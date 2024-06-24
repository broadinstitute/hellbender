package org.broadinstitute.hellbender.utils.reference;

import com.google.common.io.ByteStreams;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.PeekableIterator;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.IOException;
import java.io.Serializable;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.util.*;

/**
 * A class for reading a genomic reference stored in .2bit format. This format is defined here:
 *
 *     http://genome.ucsc.edu/FAQ/FAQformat.html#format7
 *
 * Note that this class stores the entire .2bit reference in memory (in packed form) to facilitate
 * Spark broadcasts.
 *
 * Supported public operations are:
 * {@link #getSequenceDictionary}
 * {@link #getReferenceBases(SimpleInterval)}
 * {@link #close}
 *
 * The constructor gives the option to either preserve lowercase (masked) bases, or uppercase them.
 * Ambiguity codes are not supported by the 2bit format.
 */
public class TwoBitReference implements Serializable, AutoCloseable {
    private static final long serialVersionUID = 1L;

    // Magic value at the start of a 2bit file
    public static final int TWO_BIT_SIGNATURE = 0x1A412743;

    // Version of the 2bit format we support
    public static final int TWO_BIT_SUPPORTED_VERSION = 0;

    // File extension for 2bit files
    public static final String TWO_BIT_EXTENSION = ".2bit";

    // Size of the 2bit header in bytes
    public static final int HEADER_LENGTH_IN_BYTES = 16;

    private final GATKPath referencePath;

    // If true, preserve lowercase (masked) bases, otherwise uppercase all bases
    private final boolean preserveCase;

    // The complete in-memory contents of the 2bit file. Initialized as read-only in the constructor
    private final ByteBuffer rawBytes;

    // Byte order to use, as defined in the 2bit header
    private ByteOrder byteOrder;

    // The total number of sequences in the file
    private int sequenceCount;

    // LinkedHashMap because the order of the sequence records does matter
    private LinkedHashMap<String, TwoBitSequenceRecord> sequenceRecords;

    private SAMSequenceDictionary sequenceDictionary;

    /**
     * Creates a TwoBitReference that will uppercase all bases
     *
     * @param referencePath path to the 2bit reference
     */
    public TwoBitReference( final GATKPath referencePath ) {
        this(referencePath, false);
    }

    /**
     * Creates a TwoBitReference, and specifies whether bases should be uppercased
     *
     * @param referencePath path to the 2bit reference
     * @param preserveCase if true, preserve lowercase (masked) bases, otherwise uppercase all bases
     */
    public TwoBitReference( final GATKPath referencePath, final boolean preserveCase ) {
        Utils.nonNull(referencePath);
        Utils.validateArg(referencePath.getURI().getPath().endsWith(TWO_BIT_EXTENSION), "Twobit reference must end with a " + TWO_BIT_EXTENSION + " extension");

        // Since the Java NIO APIs we're relying on use signed integers for file positions, we're
        // limited to supporting files up to MAX_INT (2^31 - 1) bytes in size, even though the 2bit
        // format should theoretically allow for files up to 2^32 bytes (ie., MAX_UNSIGNED_INT).
        try {
            final long referenceSizeInBytes = Files.size(referencePath.toPath());
            if ( referenceSizeInBytes > Integer.MAX_VALUE ) {
                throw new UserException.CouldNotReadInputFile(referencePath,
                        "We currently only support 2bit references up to " + Integer.MAX_VALUE +
                                " bytes in size, but this file is " + referenceSizeInBytes + " bytes. " +
                                "This limitation may be removed in a future release.");
            }
        } catch ( IOException e ) {
            throw new UserException.CouldNotReadInputFile(referencePath, "Error while checking the size of the 2bit reference file", e);
        }

        this.referencePath = referencePath;
        this.preserveCase = preserveCase;

        try {
            // Load the raw contents of the 2bit file into a read-only in-memory ByteBuffer
            rawBytes = ByteBuffer.wrap(ByteStreams.toByteArray(referencePath.getInputStream())).asReadOnlyBuffer();
        } catch ( IOException e ) {
            throw new UserException.CouldNotReadInputFile(referencePath, "Unable to load bytes from 2bit input file", e);
        }

        sequenceRecords = new LinkedHashMap<>();

        readHeader();
        readSequenceRecordIndexAndMetadata();

        sequenceDictionary = new SAMSequenceDictionary();
        for ( final TwoBitSequenceRecord twoBitSequenceRecord : sequenceRecords.values() ) {
            sequenceDictionary.addSequence(new SAMSequenceRecord(twoBitSequenceRecord.getSequenceName(), twoBitSequenceRecord.getDNASize()));
        }
    }

    /**
     * Reads the 2bit file header, and sets byte order for the rest of the file
     *
     * The file begins with a 16-byte header containing the following fields:
     *
     * signature - the number 0x1A412743 in the architecture of the machine that created the file
     * version - zero for now. Readers should abort if they see a version number higher than 0
     * sequenceCount - the number of sequences in the file
     * reserved - always zero for now
     */
    private void readHeader() {
        setByteOrder(ByteOrder.LITTLE_ENDIAN);
        int signature = rawBytes.getInt();
        if ( signature != TWO_BIT_SIGNATURE ) {
            setByteOrder(ByteOrder.BIG_ENDIAN);
            rawBytes.position(0);
            signature = rawBytes.getInt();

            if ( signature != TWO_BIT_SIGNATURE ) {
                throw new UserException.MalformedFile(referencePath,
                        "File does not start with the required 2bit signature value " + TWO_BIT_SIGNATURE);
            }
        }

        final int version = rawBytes.getInt();
        if ( version != TWO_BIT_SUPPORTED_VERSION ) {
            throw new UserException.CouldNotReadInputFile(referencePath,
                    "File is version " + version + ", but we only support version " + TWO_BIT_SUPPORTED_VERSION);
        }

        this.sequenceCount = rawBytes.getInt();
        if ( sequenceCount < 0 ) {
            throw new UserException.MalformedFile(referencePath,
                    "Negative sequence count in 2bit header: " + sequenceCount);
        }

        final int reservedValue = rawBytes.getInt();
        if ( reservedValue != 0 ) {
            throw new UserException.MalformedFile(referencePath,
                "Reserved header value must be 0 according to the 2bit spec, but found " + reservedValue);
        }
    }

    /**
     * Reads and parses the sequence record index and metadata
     *
     * The header is followed by a file index, which contains one entry for each sequence. Each index entry contains three fields:
     *
     * nameSize - a byte containing the length of the name field
     * name - the sequence name itself (in ASCII-compatible byte string), of variable length depending on nameSize
     * offset - the 32-bit offset of the sequence data relative to the start of the file, not aligned to any 4-byte padding boundary
     *
     * The index is followed by the sequence records, which contain nine fields:
     *
     * dnaSize - number of bases of DNA in the sequence
     * nBlockCount - the number of blocks of Ns in the file (representing unknown sequence)
     * nBlockStarts - an array of length nBlockCount of 32 bit integers indicating the (0-based) starting position of a block of Ns
     * nBlockSizes - an array of length nBlockCount of 32 bit integers indicating the length of a block of Ns
     * maskBlockCount - the number of masked (lower-case) blocks
     * maskBlockStarts - an array of length maskBlockCount of 32 bit integers indicating the (0-based) starting position of a masked block
     * maskBlockSizes - an array of length maskBlockCount of 32 bit integers indicating the length of a masked block
     * reserved - always zero for now
     * packedDna - the DNA packed to two bits per base, represented as so: T - 00, C - 01, A - 10, G - 11. The first base is in the most significant 2-bit byte; the last base is in the least significant 2 bits. For example, the sequence TCAG is represented as 00011011.
     */
    private void readSequenceRecordIndexAndMetadata() {
        final ByteBuffer bytes = independentBufferView();
        bytes.position(HEADER_LENGTH_IN_BYTES);  // seek past the header

        // First read the index, with the names of each sequence and start offsets for the sequence records
        for ( int i = 0; i < sequenceCount; i++ ) {
            final byte nameSize = rawBytes.get();
            final byte[] nameBytes = new byte[nameSize];
            rawBytes.get(nameBytes);
            final String sequenceName = new String(nameBytes);
            final int sequenceRecordStartOffset = rawBytes.getInt();

            sequenceRecords.put(sequenceName, new TwoBitSequenceRecord(sequenceName, i, sequenceRecordStartOffset));
        }

        // Next read the metadata for each sequence (DNA size, N blocks, masked blocks, and start offset for the
        // actual bases) from the sequence records
        for ( final TwoBitSequenceRecord sequenceRecord : sequenceRecords.values() ) {
            bytes.position(sequenceRecord.getSequenceRecordStartOffset());

            sequenceRecord.setDNASize(bytes.getInt());

            final int nBlockCount = bytes.getInt();
            final int[] nBlockStarts = new int[nBlockCount];
            final int[] nBlockSizes = new int[nBlockCount];
            for ( int i = 0; i < nBlockCount; i++ ) {
                nBlockStarts[i] = bytes.getInt();
            }
            for ( int i = 0; i < nBlockCount; i++ ) {
                nBlockSizes[i] = bytes.getInt();
            }
            sequenceRecord.setNBlocks(convertTwoBitBlocksToSortedMergedIntervals(sequenceRecord.getSequenceName(), nBlockStarts, nBlockSizes));

            final int maskBlockCount = bytes.getInt();
            final int[] maskBlockStarts = new int[maskBlockCount];
            final int[] maskBlockSizes = new int[maskBlockCount];
            for ( int i = 0; i < maskBlockCount; i++ ) {
                maskBlockStarts[i] = bytes.getInt();
            }
            for ( int i = 0; i < maskBlockCount; i++ ) {
                maskBlockSizes[i] = bytes.getInt();
            }
            sequenceRecord.setMaskedBlocks(convertTwoBitBlocksToSortedMergedIntervals(sequenceRecord.getSequenceName(), maskBlockStarts, maskBlockSizes));

            final int reservedValue = bytes.getInt();
            if ( reservedValue != 0 ) {
                throw new UserException.MalformedFile(referencePath,
                        "Reserved value for sequence " + sequenceRecord.getSequenceName() + " must be zero, but is non-zero");
            }

            sequenceRecord.setSequenceBasesStartOffset(bytes.position());
        }
    }

    private List<SimpleInterval> convertTwoBitBlocksToSortedMergedIntervals( final String contig, final int[] zeroBasedBlockStarts, final int[] blockSizes ) {
        final List<SimpleInterval> intervals = new ArrayList<>();
        for ( int i = 0; i < zeroBasedBlockStarts.length; i++ ) {
            final int oneBasedInclusiveStart = zeroBasedBlockStarts[i] + 1;
            final int oneBasedInclusiveEnd = oneBasedInclusiveStart + blockSizes[i] - 1;
            intervals.add(new SimpleInterval(contig, oneBasedInclusiveStart, oneBasedInclusiveEnd));
        }

        return IntervalUtils.sortAndMergeIntervalsFromSameContig(intervals, IntervalMergingRule.ALL);
    }

    private void setByteOrder( final ByteOrder byteOrder ) {
        this.byteOrder = byteOrder;
        rawBytes.order(byteOrder);
    }

    /**
     * @return The sequence dictionary for this TwoBit reference
     */
    public SAMSequenceDictionary getSequenceDictionary() {
        return sequenceDictionary;
    }

    /**
     * @param interval Genomic interval to query (must be completely on the reference)
     * @return the reference bases covering the specified interval
     */
    public ReferenceSequence getReferenceBases( final SimpleInterval interval ) {
        final TwoBitSequenceRecord sequenceRecord = sequenceRecords.get(interval.getContig());
        Utils.nonNull(sequenceRecord, "Sequence " + interval.getContig() + " not present in 2bit file " + referencePath);
        if ( interval.getEnd() > sequenceRecord.getDNASize() ) {
            throw new IllegalArgumentException("Query interval " + interval + " ends after the end of sequence " +
                    sequenceRecord.getSequenceName() + " (length " + sequenceRecord.getDNASize() + ")");
        }

        final byte[] sequenceBases = new byte[interval.getLengthOnReference()];
        final TwoBitBaseDecoder decoder = new TwoBitBaseDecoder(independentBufferView(), sequenceRecord, interval, preserveCase);
        int baseOffset = 0;

        while ( decoder.hasNext() ) {
            if ( baseOffset >= sequenceBases.length ) {
                throw new IllegalStateException("2bit base decoder returned an unexpected number of bases");
            }

            sequenceBases[baseOffset] = decoder.next();
            baseOffset++;
        }

        return new ReferenceSequence(sequenceRecord.getSequenceName(), sequenceRecord.getSequenceIndex(), sequenceBases);
    }

    /**
     * @return An independent view of the ByteBuffer containing the 2bit file contents, with its own
     * position, limit, and mark values. The actual bytes in the buffer are not copied.
     */
    private ByteBuffer independentBufferView() {
        final ByteBuffer newView = rawBytes.asReadOnlyBuffer(); // does not duplicate the content, just creates a new view of it
        newView.order(byteOrder);   // byte order is not preserved by asReadOnlyBuffer()
        newView.clear();  // resets the position/limit/mark, does not erase the buffer contents
        return newView;
    }

    @Override
    public void close() {
        // No-op for now, just here to allow seamless use in Autocloseable contexts.
        // In a future version that does not store the entire reference in memory, we will need this.
    }

    /**
     * An Iterator that decodes the raw TwoBit reference for a specified interval.
     * Each interval requires a new Iterator instance.
     */
    private static class TwoBitBaseDecoder implements Iterator<Byte>, Iterable<Byte>, Serializable {
        private static final long serialVersionUID = 1L;

        // A two-bit base encoding used as an index into this array will produce the corresponding base
        private static final Byte[] twoBitEncodingToBase = { 'T', 'C', 'A', 'G' };

        private final TwoBitSequenceRecord sequenceRecord;

        private final ByteBuffer rawBytes;

        private final SimpleInterval intervalToDecode;

        private final boolean preserveCase;

        // 1-based inclusive genomic position within our sequence's contig
        private int currentOneBasedGenomicPosition;

        // Raw value of the current packed byte we're positioned on
        private byte currentByte;

        // 0-based bit offset from the most significant bit in the current byte, marking the start of the next 2bit encoding
        // Can be 0, 2, 4, or 6
        private int currentBitOffset;

        // The base that will be returned by the next call to next(), or null if the iteration is complete
        private Byte nextBase;

        private PeekableIterator<SimpleInterval> nBlockIterator;
        private PeekableIterator<SimpleInterval> maskedBlockIterator;

        /**
         * @param rawBytes Raw byte buffer containing the complete TwoBit data. Should be an independent view of the
         *                 buffer not shared with any other concurrent code.
         * @param sequenceRecord TwoBitSequenceRecord for the contig we're accessing
         * @param intervalToDecode Query interval to decode. Must match the TwoBitSequenceRecord's contig.
         * @param preserveCase If true, preserve lowercase (masked) bases, otherwise convert all bases to uppercase.
         */
        public TwoBitBaseDecoder( final ByteBuffer rawBytes, final TwoBitSequenceRecord sequenceRecord, final SimpleInterval intervalToDecode, final boolean preserveCase ) {
            Utils.validate(intervalToDecode.getEnd() <= sequenceRecord.getDNASize(), "Interval " + intervalToDecode + " ends past the end of contig " + sequenceRecord.getSequenceName());
            Utils.validate(intervalToDecode.getContig().equals(sequenceRecord.getSequenceName()), "Interval " + intervalToDecode + " is on the wrong contig for sequence " + sequenceRecord.getSequenceName());

            this.sequenceRecord = sequenceRecord;
            this.intervalToDecode = intervalToDecode;
            this.preserveCase = preserveCase;

            this.currentOneBasedGenomicPosition = intervalToDecode.getStart();
            this.rawBytes = rawBytes;
            rawBytes.position(sequenceRecord.getSequenceBasesStartOffset() + ((currentOneBasedGenomicPosition - 1) / 4));
            this.currentByte = rawBytes.get();
            this.currentBitOffset = ((currentOneBasedGenomicPosition - 1) % 4) * 2;

            nBlockIterator = new PeekableIterator<>(sequenceRecord.getNBlocks().iterator());
            maskedBlockIterator = new PeekableIterator<>(sequenceRecord.getMaskedBlocks().iterator());
            advanceBlockIterators();

            loadNextBase();
        }

        @Override
        public Iterator<Byte> iterator() {
            return this;
        }

        @Override
        public boolean hasNext() {
            return nextBase != null;
        }

        @Override
        public Byte next() {
            if ( ! hasNext() ) {
                throw new NoSuchElementException("next() called when iterator exhausted");
            }

            Byte toReturn = nextBase;
            loadNextBase();
            return toReturn;
        }

        private void loadNextBase() {
            if ( iterationComplete() ) {
                nextBase = null;
                return;
            }

            nextBase = decodeBaseAtCurrentPosition();
            if ( isPositionedAtNBase() ) {
                nextBase = 'N';
            } else if ( isPositionedAtMaskedBase() && preserveCase ) {   // Masked 'N' is not supported
                nextBase = (byte) Character.toLowerCase(nextBase);
            }

            advancePosition();
        }

        private Byte decodeBaseAtCurrentPosition() {
            // The DNA is packed to two bits per base, represented as so: T - 00, C - 01, A - 10, G - 11.
            // The first base is in the most significant 2-bit byte; the last base is in the least significant 2 bits.
            return twoBitEncodingToBase[(currentByte >> (6 - currentBitOffset)) & 3];
        }

        private boolean iterationComplete() {
            return currentOneBasedGenomicPosition > intervalToDecode.getEnd();
        }

        private void advancePosition() {
            currentOneBasedGenomicPosition++;
            if ( iterationComplete() ) {
                return;
            }

            // Bit offsets can be 0, 2, 4, or 6, corresponding to the start bit for each of the four bases encoded in a byte
            if ( currentBitOffset <= 4 ) {
                currentBitOffset += 2;
            } else {    // Advance to the next full byte
                currentByte = rawBytes.get();
                currentBitOffset = 0;
            }

            advanceBlockIterators();
        }

        private void advanceBlockIterators() {
            advanceBlockIterator(nBlockIterator);
            advanceBlockIterator(maskedBlockIterator);
        }

        private void advanceBlockIterator( final PeekableIterator<SimpleInterval> iterator ) {
            while ( iterator.hasNext() && iterator.peek().getEnd() < currentOneBasedGenomicPosition ) {
                iterator.next(); // we've passed this interval on the genome, so discard it
            }
        }

        private boolean isPositionedAtNBase() {
            final SimpleInterval currentNBlock = nBlockIterator.peek();

            return currentNBlock != null && currentNBlock.getStart() <= currentOneBasedGenomicPosition &&
                    currentNBlock.getEnd() >= currentOneBasedGenomicPosition;
        }

        private boolean isPositionedAtMaskedBase() {
            final SimpleInterval currentMaskedBlock = maskedBlockIterator.peek();

            return currentMaskedBlock != null && currentMaskedBlock.getStart() <= currentOneBasedGenomicPosition &&
                    currentMaskedBlock.getEnd() >= currentOneBasedGenomicPosition;
        }
    }

    /**
     * Stores the raw metadata for a specific sequence from the TwoBit file
     */
    private static class TwoBitSequenceRecord implements Serializable {
        private static final long serialVersionUID = 1L;
        private final String sequenceName;

        private final int sequenceIndex;

        private final int sequenceRecordStartOffset;

        private int sequenceBasesStartOffset;

        private int dnaSize;

        private List<SimpleInterval> nBlocks;

        private List<SimpleInterval> maskedBlocks;

        public TwoBitSequenceRecord( final String sequenceName, final int sequenceIndex, final int sequenceRecordStartOffset ) {
            this.sequenceName = sequenceName;
            this.sequenceIndex = sequenceIndex;
            this.sequenceRecordStartOffset = sequenceRecordStartOffset;
        }

        public String getSequenceName() {
            return sequenceName;
        }

        public int getSequenceIndex() {
            return sequenceIndex;
        }

        public int getSequenceRecordStartOffset() {
            return sequenceRecordStartOffset;
        }

        public void setDNASize( final int dnaSize ) {
            this.dnaSize = dnaSize;
        }

        public int getDNASize() {
            return dnaSize;
        }

        public void setSequenceBasesStartOffset( final int sequenceBasesStartOffset ) {
            this.sequenceBasesStartOffset = sequenceBasesStartOffset;
        }

        public int getSequenceBasesStartOffset() {
            return sequenceBasesStartOffset;
        }

        public void setNBlocks( final List<SimpleInterval> nBlocks ) {
            this.nBlocks = nBlocks;
        }

        public List<SimpleInterval> getNBlocks() {
            return nBlocks;
        }

        public void setMaskedBlocks( final List<SimpleInterval> maskedBlocks ) {
            this.maskedBlocks = maskedBlocks;
        }

        public List<SimpleInterval> getMaskedBlocks() {
            return maskedBlocks;
        }
    }
}
