package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.*;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.TemplateFragmentOrdinal;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembler;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Memory-economical utilities for producing a FASTQ file.
 */
public class SVFastqUtils {

    public static final char HEADER_PREFIX_CHR = '@';
    public static final char FRAGMENT_NUMBER_SEPARATOR_CHR = '/';
    public static final char HEADER_FIELD_SEPARATOR_CHR = '\t';
    public static final String HEADER_FIELD_SEPARATOR_REGEXP = "\\t"; // make sure its consistent with {@link #HEADER_FIELD_SEPARATOR_CHR}
    public static final char HEADER_FIELD_EQUAL_CHR = '=';
    public static final String HEADER_FIELD_LIST_SEPARATOR_STR = AlignmentInterval.SA_TAG_INTERVAL_SEPARATOR_STR;
    private static final char LINE_SEPARATOR_CHR = '+';
    private static final String MAPPING_FIELD_NAME = "mapping";
    private static final String UNMAPPED_STR = "*";
    private static final String HEADER_FIELD_SEPARATOR_STR = "" + HEADER_FIELD_SEPARATOR_CHR;
    private static final String MAPPING_FIELD_EQUAL_TO = MAPPING_FIELD_NAME + HEADER_FIELD_EQUAL_CHR;
    private static final String UNMAPPED_DESCRIPTION_STR = MAPPING_FIELD_EQUAL_TO + UNMAPPED_STR;
    private static final Pattern FASTQ_READ_HEADER_PATTERN = Pattern.compile("^" + HEADER_PREFIX_CHR +
            "([^" + HEADER_FIELD_SEPARATOR_REGEXP +  "]+)(" + HEADER_FIELD_SEPARATOR_REGEXP + "(.*))?$");

    private static final Pattern SA_TAG_ALN_INTERVAL_SEPARATOR_PATTERN = Pattern.compile(HEADER_FIELD_LIST_SEPARATOR_STR);

    public static final class Mapping implements Locatable {

        final List<AlignmentInterval> intervals;

        public Mapping(final String str) {
            Utils.nonNull(str);
            if (str.equals(UNMAPPED_STR)) {
                intervals = Collections.emptyList();
            } else {
                intervals = Collections.unmodifiableList(SA_TAG_ALN_INTERVAL_SEPARATOR_PATTERN.splitAsStream(str)
                    .filter(s -> !s.isEmpty())
                    .map(AlignmentInterval::new)
                    .collect(Collectors.toList()));
            }
        }

        public Mapping(final GATKRead read) {
            Utils.nonNull(read);
            if (read.isUnmapped()) {
                intervals = Collections.emptyList();
            } else if (!read.hasAttribute(SAMTag.SA.name())) {
                if (read.isSupplementaryAlignment()) {
                    throw new GATKException("a supplementary alignment read record must supply a SA tag. Perhaps this constraint isn't " +
                            "satisfied by the aligner used to produce the record");
                }
                intervals = Collections.singletonList(new AlignmentInterval(read));
            } else if (!read.isSupplementaryAlignment()) {
                if (!read.hasAttribute(SAMTag.SA.name())) {
                    intervals = Collections.singletonList(new AlignmentInterval(read));
                } else {
                    intervals = Collections.unmodifiableList(Stream.concat(
                            Stream.of(new AlignmentInterval(read)),
                            SA_TAG_ALN_INTERVAL_SEPARATOR_PATTERN.splitAsStream(read.getAttributeAsString(SAMTag.SA.name()))
                                    .filter(s -> !s.isEmpty())
                                    .map(AlignmentInterval::new))
                            .collect(Collectors.toList()));
                }
            } else {
                final List<AlignmentInterval> otherIntervals = SA_TAG_ALN_INTERVAL_SEPARATOR_PATTERN.splitAsStream(read.getAttributeAsString(SAMTag.SA.name()))
                        .filter(s -> !s.isEmpty())
                        .map(AlignmentInterval::new)
                        .collect(Collectors.toList());
                if (otherIntervals.isEmpty()) {
                    throw new GATKException("the SA tag on a supplementary read record does not have alignments, perhaps this constraint is not satisfied by the aligner used to generate these records");
                }
                final List<AlignmentInterval> newIntervals = new ArrayList<>(otherIntervals.size() + 1);
                final AlignmentInterval readInterval = new AlignmentInterval(read);
                final AlignmentInterval primaryInterval = otherIntervals.get(0);
                newIntervals.add(primaryInterval);
                newIntervals.add(readInterval);
                newIntervals.addAll(otherIntervals.subList(1, otherIntervals.size()));
                intervals = Collections.unmodifiableList(newIntervals);
            }
            if (!intervals.isEmpty() && intervals.get(0).cigarAlong5to3DirectionOfContig.containsOperator(CigarOperator.H)) {
                throw new GATKException("the interval/cigar of a non-supplementary record must not contain hard-clips, perhaps this constraint is not satisfied by the aligner used to generate these records");
            }
        }

        public boolean isMapped() {
            return !intervals.isEmpty();
        }

        /**
         * Returns the contig as per the primary interval, or {@code null} if unmapped.
         * @return may return a {@code null}.
         */
        @Override
        public String getContig() {
            return isMapped() ? intervals.get(0).referenceSpan.getContig() : null;
        }

        /**
         * Returns the primary interval start position, or {@link SAMRecord#NO_ALIGNMENT_START} if unmapped.
         * @return {@link SAMRecord#NO_ALIGNMENT_START} or a strictly positive value.
         */
        @Override
        public int getStart() {
            return isMapped() ? intervals.get(0).referenceSpan.getStart() : SAMRecord.NO_ALIGNMENT_START;
        }

        /**
         * Returns the primary interval end position, or {@link SAMRecord#NO_ALIGNMENT_START} if unmapped.
         * @return {@link SAMRecord#NO_ALIGNMENT_START} or a strictly positive value.
         */
        @Override
        public int getEnd() {
            return isMapped() ? intervals.get(0).referenceSpan.getEnd() : SAMRecord.NO_ALIGNMENT_START;
        }

        /**
         * Returns the cigar of the primary interval, or the empty cigar if {@code unmapped}.
         * @return
         */
        public Cigar getCigar() {
            return isMapped() ? intervals.get(0).cigarAlongReference() : new Cigar();
        }

        public boolean isForwardStrand() {
            if (isMapped()) {
                return intervals.get(0).forwardStrand;
            } else {
                throw new UnsupportedOperationException("no forward or backward strand if unmapped");
            }
        }

        public AlignmentInterval getPrimaryInterval() {
            return intervals.isEmpty() ? null : intervals.get(0);
        }

        public List<AlignmentInterval> getSupplementaryIntervals() {
            return intervals.isEmpty() ? intervals : intervals.subList(1, intervals.size());
        }

        public List<AlignmentInterval> getAllIntervals() {
            return intervals;
        }

        /**
         * Generates the mapping string for a read. This contains the read own coordinates
         * plus the content of the SA tag if present.
         * @param read the source read.
         * @return never null. {@link #UNMAPPED_DESCRIPTION_STR} if the read is unmapped.
         */
        public static String toString(final GATKRead read) {
            Utils.nonNull(read);
            if ( read.isUnmapped() || read.getCigar().getPaddedReferenceLength() == 0 ) {
                return UNMAPPED_STR;
            } else {
                final StringBuilder builder = new StringBuilder(100);
                new AlignmentInterval(read).appendSATagString(builder);
                if (read.hasAttribute(SAMTag.SA.name())) {
                    builder.append(HEADER_FIELD_LIST_SEPARATOR_STR)
                           .append(read.getAttributeAsString(SAMTag.SA.name()));
                }
                if (builder.lastIndexOf(HEADER_FIELD_LIST_SEPARATOR_STR) == builder.length() - 1) {
                    builder.setLength(builder.length() - 1);
                }
                return builder.toString();
            }
        }

        /**
         * Generates the mapping string as it would be included in the Fastq read header.
         * @return never null. {@link #UNMAPPED_DESCRIPTION_STR} if the read is unmapped.
         */
        public String toString() {
            if (isMapped()) {
                final StringBuilder builder = new StringBuilder(100);
                for (final AlignmentInterval interval : intervals) {
                    interval.appendSATagString(builder).append(HEADER_FIELD_LIST_SEPARATOR_STR);
                }
                builder.setLength(builder.length() - 1);
                return builder.toString();
            } else {
                return UNMAPPED_STR;
            }
        }
    }

    @DefaultSerializer(FastqRead.Serializer.class)
    public static final class FastqRead implements FermiLiteAssembler.BasesAndQuals {
        private final String header; // looks just like a line in a FASTQ file: "@readName"
        private final byte[] bases;  // contains bytes with values like 'A', 'C', 'G', and 'T'
        private final byte[] quals;  // contains phred-scaled bytes without any bias --
                                     //   a bias of 33 is added when writing, and removed when reading

        public FastqRead(final GATKRead read) {
            this(read, true);
        }

        public FastqRead( final GATKRead read, final boolean includeMappingLocation ) {
            this.header = composeHeaderLine(read, includeMappingLocation);
            this.bases = read.getBases();
            this.quals = read.getBaseQualities();
            if (!read.isUnmapped() && read.isReverseStrand()) {
                SequenceUtil.reverseComplement(this.bases);
                SequenceUtil.reverseQualities(this.quals);
            }
        }

        private static String composeHeaderLine(final GATKRead read, final boolean includeMappingLocation) {
            Utils.nonNull(read);
            return HEADER_PREFIX_CHR + read.getName() + TemplateFragmentOrdinal.forRead(read)
                    + (includeMappingLocation ?  (HEADER_FIELD_SEPARATOR_STR + MAPPING_FIELD_EQUAL_TO + Mapping.toString(read)) : "");
        }

        @VisibleForTesting
        FastqRead(final String header, final byte[] bases, final byte[] quals) {
            this.header = header;
            this.bases = bases;
            this.quals = quals;
        }

        private FastqRead( final Kryo kryo, final Input input ) {
            header = input.readString();
            final int nBases = input.readInt();
            bases = new byte[nBases];
            input.readBytes(bases);
            quals = new byte[nBases];
            input.readBytes(quals);
        }

        /**
         * Returns the header line of this Fastq read starting with '@' followed by the read id and description separated by tab characters.
         * @return never {@code null}.
         */
        public String getHeader() {
            return header;
        }

        public String getId() {
            final String[] headerParts = header.split(HEADER_FIELD_SEPARATOR_STR);
            return headerParts[0].substring(1); // skip the '@'.
        }

        public String getName() { final String id = getId();
            final int fragmentNumberSeparatorIndex = id.lastIndexOf(FRAGMENT_NUMBER_SEPARATOR_CHR);
            return fragmentNumberSeparatorIndex >= 0 ? id.substring(0, fragmentNumberSeparatorIndex) : id;
        }

        public String getDescription() {
            final int tabIndex = header.indexOf(HEADER_FIELD_SEPARATOR_CHR);
            return tabIndex >= 0 ? header.substring(tabIndex + 1) : "";
        }

        public Mapping getMapping() {
            final int mappingEqualToIndex = header.indexOf(MAPPING_FIELD_EQUAL_TO);
            if (mappingEqualToIndex < 0) {
                return null;
            } else {
                final int nextFieldSeperator = header.indexOf(HEADER_FIELD_SEPARATOR_STR, mappingEqualToIndex);
                final String mappingString = nextFieldSeperator < 0
                        ? header.substring(mappingEqualToIndex + MAPPING_FIELD_EQUAL_TO.length())
                        : header.substring(mappingEqualToIndex + MAPPING_FIELD_EQUAL_TO.length(), nextFieldSeperator);
                return new Mapping(mappingString);
            }
        }

        /** returns bases with values like 'A', 'C', 'G', and 'T' */
        @Override public byte[] getBases() { return bases; }
        /** returns phred-scaled quals with no bias */
        @Override public byte[] getQuals() { return quals; }

        private void serialize( final Kryo kryo, final Output output ) {
            output.writeAscii(header);
            output.writeInt(bases.length);
            output.writeBytes(bases);
            output.writeBytes(quals);
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<FastqRead> {
            @Override
            public void write( final Kryo kryo, final Output output, final FastqRead read ) {
                read.serialize(kryo, output);
            }

            @Override
            public FastqRead read( final Kryo kryo, final Input input, final Class<FastqRead> type ) {
                return new FastqRead(kryo, input);
            }
        }
    }

    public static List<FastqRead> readFastqFile( final String fileName ) {
        final List<FastqRead> reads;
        try ( final BufferedReader reader = new BufferedReader(new InputStreamReader(BucketUtils.openFile(fileName))) ) {
            reads = readFastqStream(reader, fileName);
        }
        catch ( final IOException ioe ) {
            throw new GATKException("Can't read "+fileName, ioe);
        }
        return reads;
    }

    public static List<FastqRead> readFastqStream( final BufferedReader reader, final String fileName ) throws IOException {
        final int INITIAL_CAPACITY = 10000; // a wild guess, just something not too crazy small
        final List<FastqRead> reads = new ArrayList<>(INITIAL_CAPACITY);
        String header;
        int lineNo = 0;
        while ( (header = reader.readLine()) != null ) {
            lineNo += 1;
            if (!FASTQ_READ_HEADER_PATTERN.matcher(header).find()) {
                throw new GATKException("In FASTQ file "+fileName+" sequence identifier line does not start with @ on line "+lineNo);
            }
            final String callLine = reader.readLine();
            lineNo += 1;
            if ( callLine == null ) {
                throw new GATKException("In FASTQ file "+fileName+" file truncated: missing calls.");
            }
            final String sepLine = reader.readLine();
            lineNo += 1;
            if ( sepLine == null ) {
                throw new GATKException("In FASTQ file "+fileName+" file truncated: missing + line.");
            }
            if ( sepLine.length() < 1 || sepLine.charAt(0) != LINE_SEPARATOR_CHR ) {
                throw new GATKException("In FASTQ file " + fileName + " separator line does not start with + on line " + lineNo);
            }
            final String qualLine = reader.readLine();
            lineNo += 1;
            if ( qualLine == null ) {
                throw new UserException.BadInput("In FASTQ file "+fileName+" file truncated: missing quals.");
            }
            if ( callLine.length() != qualLine.length() ) {
                throw new UserException.BadInput("In FASTQ file "+fileName+" there are "+qualLine.length()+
                        " quality scores on line "+lineNo+" but there are "+callLine.length()+" base calls.");
            }
            final byte[] quals = qualLine.getBytes();
            SAMUtils.fastqToPhred(quals);
            final byte[] calls = callLine.getBytes();
            reads.add(new FastqRead(header, calls, quals));
        }
        return reads;
    }

    /** Convert a read's name into a FASTQ record sequence ID */
    public static String readToFastqSeqId( final GATKRead read, final boolean includeMappingLocation ) {
        final String nameSuffix = TemplateFragmentOrdinal.forRead(read).nameSuffix();
        final String description;
        if ( includeMappingLocation ) {
            final Mapping mapping = new Mapping(read);
            description = HEADER_FIELD_SEPARATOR_CHR + Mapping.toString(read);
        } else {
            description = "";
        }
        return read.getName() + nameSuffix + HEADER_FIELD_SEPARATOR_CHR + description;
    }

    /** Write a list of FASTQ records into a file. */
    public static void writeFastqFile(final String fileName, final Iterator<FastqRead> fastqReadItr ) {
        try ( final OutputStream writer =
                      new BufferedOutputStream(BucketUtils.createFile(fileName)) ) {
            writeFastqStream(writer, fastqReadItr);
        } catch ( final IOException ioe ) {
            throw new GATKException("Can't write "+fileName, ioe);
        }
    }

    public static void writeFastqStream( final OutputStream writer, final Iterator<FastqRead> fastqReadItr )
        throws IOException {
        int index = 0;
        while ( fastqReadItr.hasNext() ) {
            final FastqRead read = fastqReadItr.next();
            final String header = read.getHeader();
            if (header.contains(" ")) {
                throw new IllegalStateException("Blank found: " + header);
            }
            if ( header == null ) writer.write(Integer.toString(++index).getBytes());
            else writer.write(header.getBytes());
            writer.write('\n');
            writer.write(read.getBases());
            writer.write('\n');
            writer.write(SVFastqUtils.LINE_SEPARATOR_CHR);
            writer.write('\n');
            final byte[] quals = read.getQuals();
            final int nQuals = quals.length;
            final byte[] fastqQuals = new byte[nQuals];
            for ( int idx = 0; idx != nQuals; ++idx ) {
                fastqQuals[idx] = (byte)SAMUtils.phredToFastq(quals[idx]);
            }
            writer.write(fastqQuals);
            writer.write('\n');
        }
    }
}
