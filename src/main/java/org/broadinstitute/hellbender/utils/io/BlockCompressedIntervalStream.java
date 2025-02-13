package org.broadinstitute.hellbender.utils.io;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.seekablestream.SeekableBufferedStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.BlockCompressedStreamConstants;
import htsjdk.tribble.*;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SVInterval;
import org.broadinstitute.hellbender.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.codecs.FeatureSink;
import org.broadinstitute.hellbender.tools.sv.SVFeaturesHeader;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.channels.SeekableByteChannel;
import java.nio.file.Path;
import java.util.*;
import java.util.function.Function;

public class BlockCompressedIntervalStream {

    // our final empty block adds to the extra information
    // it has a virtual file pointer to the start of the index
    public static final byte[] EMPTY_GZIP_BLOCK_WITH_INDEX_POINTER = {
            BlockCompressedStreamConstants.GZIP_ID1,
            (byte)BlockCompressedStreamConstants.GZIP_ID2,
            BlockCompressedStreamConstants.GZIP_CM_DEFLATE,
            BlockCompressedStreamConstants.GZIP_FLG,
            0, 0, 0, 0, // modification time
            BlockCompressedStreamConstants.GZIP_XFL,
            (byte)BlockCompressedStreamConstants.GZIP_OS_UNKNOWN,
            BlockCompressedStreamConstants.GZIP_XLEN + 12, 0,
            BlockCompressedStreamConstants.BGZF_ID1,
            BlockCompressedStreamConstants.BGZF_ID2,
            BlockCompressedStreamConstants.BGZF_LEN, 0,
            39, 0, // Total block size - 1 as little-endian short
            (byte)'I', (byte)'P', 8, 0, // index pointer extra data

            // 8-byte little-endian long representing file-pointer to beginning of index
            // (this data gets overwritten each time a stream is closed)
            1, 2, 3, 4, 5, 6, 7, 8,

            3, 0, // empty payload
            0, 0, 0, 0, // crc
            0, 0, 0, 0, // uncompressedSize
    };
    public static final int FILE_POINTER_OFFSET = 22;

    public static final String BCI_FILE_EXTENSION = ".bci";

    // each compressed block of data will have (at least) one of these as a part of the index
    // for each contig that appears in a compressed block, the SVInterval tracks the smallest
    //   starting coordinate and largest end coordinate of any object in the block
    // the filePosition member is the virtual file offset of the first object in the block (or, the
    //   first object for a new contig, if there are multiple contigs represented within the block)
    public final static class IndexEntry {
        final SVInterval interval;
        long filePosition;

        public IndexEntry( final SVInterval interval, final long filePosition ) {
            this.interval = interval;
            this.filePosition = filePosition;
        }

        public IndexEntry( final DataInputStream dis ) throws IOException {
            this.interval = new SVInterval(dis.readInt(), dis.readInt(), dis.readInt());
            this.filePosition = dis.readLong();
        }

        public SVInterval getInterval() { return interval; }
        public long getFilePosition() { return filePosition; }

        public void write( final DataOutputStream dos ) throws IOException {
            dos.writeInt(interval.getContig());
            dos.writeInt(interval.getStart());
            dos.writeInt(interval.getEnd());
            dos.writeLong(filePosition);
        }
    }

    @FunctionalInterface
    public interface WriteFunc<F extends Feature> {
        void write( F feature, Writer<F> writer ) throws IOException;
    }

    // a class for writing arbitrary objects to a block compressed stream with a self-contained index
    // the only restriction is that you must supply a lambda that writes enough state to a DataOutputStream
    //   to allow you to reconstitute the object when you read it back in later
    public static class Writer <F extends Feature> implements FeatureSink<F> {
        final String path;
        final SVFeaturesHeader header;
        final WriteFunc<F> writeFunc;
        final OutputStream os;
        final BlockCompressedOutputStream bcos;
        final DataOutputStream dos;
        Feature lastInterval;
        final List<IndexEntry> indexEntries;
        long blockFilePosition;
        int blockContig;
        int blockStart;
        int blockEnd;
        boolean firstBlockMember;

        public final static int DEFAULT_COMPRESSION_LEVEL = 6;

        public Writer( final GATKPath path,
                       final SVFeaturesHeader header,
                       final WriteFunc<F> writeFunc ) {
            this(path, header, writeFunc, DEFAULT_COMPRESSION_LEVEL);
        }

        public Writer( final GATKPath path,
                       final SVFeaturesHeader header,
                       final WriteFunc<F> writeFunc,
                       final int compressionLevel ) {
            this.path = path.toString();
            this.header = header;
            this.writeFunc = writeFunc;
            this.os = path.getOutputStream();
            this.bcos = new BlockCompressedOutputStream(os, (Path)null, compressionLevel);
            this.dos = new DataOutputStream(bcos);
            this.lastInterval = null;
            this.indexEntries = new ArrayList<>();
            this.firstBlockMember = true;
            writeHeader();
        }

        @VisibleForTesting
        public Writer( final String streamSource,
                final OutputStream os,
                final SVFeaturesHeader header,
                final WriteFunc<F> writeFunc ) {
            this.path = streamSource;
            this.header = header;
            this.writeFunc = writeFunc;
            this.os = os;
            this.bcos = new BlockCompressedOutputStream(os, (Path)null, DEFAULT_COMPRESSION_LEVEL);
            this.dos = new DataOutputStream(bcos);
            this.lastInterval = null;
            this.indexEntries = new ArrayList<>();
            this.firstBlockMember = true;
            writeHeader();
        }

        private void writeHeader() {
            try {
                writeClassAndVersion(header.getClassName(), header.getVersion());
                writeSamples(header.getSampleNames());
                writeDictionary(header.getDictionary());
                dos.flush();
            } catch ( final IOException ioe ) {
                throw new UserException("can't write header to " + path, ioe);
            }
        }

        private void writeClassAndVersion( final String className, final String version )
                throws IOException {
            dos.writeUTF(className);
            dos.writeUTF(version);
        }

        private void writeSamples( final Collection<String> sampleNames ) throws IOException {
            final int nSampleNames = sampleNames.size();
            if ( nSampleNames == 0 ) {
                throw new UserException("Attempting to create " + path +
                        " with an empty list of sample names.  It won't go well.");
            }
            dos.writeInt(nSampleNames);
            for ( final String sampleName : sampleNames ) {
                dos.writeUTF(sampleName);
            }
        }

        private void writeDictionary( final SAMSequenceDictionary dictionary ) throws IOException {
            dos.writeInt(dictionary.size());
            for ( final SAMSequenceRecord rec : dictionary.getSequences() ) {
                dos.writeInt(rec.getSequenceLength());
                dos.writeUTF(rec.getSequenceName());
            }
        }

        public DataOutputStream getStream() { return dos; }

        public int getSampleIndex( final String sampleName ) {
            final Integer sampleIndex = header.getSampleIndex(sampleName);
            if ( sampleIndex == null ) {
                throw new UserException("Can't find index for sampleName " + sampleName +
                        ".  Supply a complete list of sample names.");
            }
            return sampleIndex;
        }

        public int getContigIndex( final String contigName ) {
            final SAMSequenceRecord contig = header.getDictionary().getSequence(contigName);
            if ( contig == null ) {
                throw new UserException("can't find contig name " + contigName +
                        ".  Please check your dictionaries for completeness and consistency.");
            }
            return contig.getSequenceIndex();
        }

        @Override
        public void write( final F feature ) {
            final long prevFilePosition = bcos.getPosition();
            // write the object
            try {
                writeFunc.write(feature, this);
            } catch ( final IOException ioe ) {
                throw new UserException("can't write to " + path, ioe);
            }

            // if this is the first interval we've seen in a block, just capture the block-start data
            if ( firstBlockMember || lastInterval == null ) {
                startBlock(prevFilePosition, feature);
                return;
            }

            // if the contig changes emit a new index entry (for the previous contig) and
            //   restart tracking of the block
            if ( !feature.contigsMatch(lastInterval) ) {
                addIndexEntry();
                startBlock(prevFilePosition, feature);
                return;
            }

            // extend the tracked interval, as necessary
            blockEnd = Math.max(blockEnd, feature.getEnd());
            lastInterval = feature;

            // if writing this element caused a new block to be compressed and added to the file
            if ( isNewBlock(prevFilePosition, bcos.getPosition()) ) {
                addIndexEntry();
                firstBlockMember = true;
            }
        }

        @Override
        public void close() {
            // take care of any pending index entry, if necessary
            if ( !firstBlockMember ) {
                addIndexEntry();
            }

            try {
                dos.flush(); // complete the data block

                long indexPosition = bcos.getPosition(); // current position is the start of the index

                // write the index entries
                dos.writeInt(indexEntries.size());
                for ( final IndexEntry indexEntry : indexEntries ) {
                    indexEntry.write(dos);
                }
                dos.flush(); // and complete the block

                // write a 0-length terminator block at the end that captures the index position
                final byte[] emptyBlockWithIndexPointer =
                        Arrays.copyOf(EMPTY_GZIP_BLOCK_WITH_INDEX_POINTER,
                                EMPTY_GZIP_BLOCK_WITH_INDEX_POINTER.length);
                for ( int idx = FILE_POINTER_OFFSET; idx != FILE_POINTER_OFFSET + 8; ++idx ) {
                    emptyBlockWithIndexPointer[idx] = (byte)indexPosition;
                    indexPosition >>>= 8;
                }
                os.write(emptyBlockWithIndexPointer);

                bcos.close(false); // we've already handled the terminator block
            } catch ( final IOException ioe ) {
                throw new UserException("unable to add index and close " + path, ioe);
            }
        }

        private void startBlock( final long filePosition, final Feature interval ) {
            blockFilePosition = filePosition;
            lastInterval = interval;
            blockContig = header.getDictionary().getSequenceIndex(interval.getContig());
            blockStart = interval.getStart();
            blockEnd = interval.getEnd();
            firstBlockMember = false;
        }

        private void addIndexEntry() {
            final SVInterval blockInterval = new SVInterval(blockContig, blockStart, blockEnd);
            indexEntries.add(new IndexEntry(blockInterval, blockFilePosition));
        }
    }

    // a class for reading arbitrary objects from a block compressed stream with a self-contained index
    // the only restriction is that you must supply a lambda that reads from a DataInputStream
    //   to reconstitute the object.
    public static final class Reader <T extends Feature> implements FeatureReader<T> {
        final String path;
        final FeatureCodec<T, Reader<T>> codec;
        final long indexFilePointer;
        final BlockCompressedInputStream bcis;
        final DataInputStream dis;
        final SVFeaturesHeader header;
        final long dataFilePointer;
        SVIntervalTree<Long> index;
        boolean usedByIterator;

        public Reader( final FeatureInput<T> inputDescriptor,
                       final FeatureCodec<T, Reader<T>> codec,
                       final int cloudPrefetchBufferSize,
                       final int cloudIndexPrefetchBufferSize ) {
            this.path = inputDescriptor.toPath().toUri().toString();
            this.codec = codec;
            final SeekableStream ss;
            final boolean prefetchable = BucketUtils.isEligibleForPrefetching(inputDescriptor);
            final Function<SeekableByteChannel, SeekableByteChannel> wrapper =
                    BucketUtils.getPrefetchingWrapper(prefetchable ? cloudPrefetchBufferSize : 0);
            try {
                ss = SeekableStreamFactory.getInstance().getStreamFor(path, wrapper);
            } catch ( final IOException ioe ) {
                throw new UserException("unable to open " + path, ioe);
            }
            this.indexFilePointer = findIndexFilePointer(ss);
            this.bcis = new BlockCompressedInputStream(new SeekableBufferedStream(ss));
            this.dis = new DataInputStream(bcis);
            this.header = readHeader();
            this.dataFilePointer = bcis.getPosition(); // having read header, we're pointing at the data
            this.index = null;
            this.usedByIterator = false;
            final String expectedClassName = codec.getFeatureType().getSimpleName();
            if ( !header.getClassName().equals(expectedClassName) ) {
                throw new UserException("can't use " + path + " to read " + expectedClassName +
                        " features -- it contains " + header.getClassName() + " features");
            }
        }

        public Reader( final Reader<T> reader ) {
            this.path = reader.path;
            this.codec = reader.codec;
            this.indexFilePointer = reader.indexFilePointer;
            try {
                this.bcis = new BlockCompressedInputStream(
                        new SeekableBufferedStream(
                                SeekableStreamFactory.getInstance().getStreamFor(path)));
            } catch ( final IOException ioe ) {
                throw new UserException("unable to clone stream for " + path, ioe);
            }
            this.dis = new DataInputStream(bcis);
            this.header = reader.header;
            this.dataFilePointer = reader.dataFilePointer;
            this.index = reader.index;
            this.usedByIterator = true;
        }

        @VisibleForTesting
        public Reader( final String inputStreamName,
                       final SeekableStream ss,
                       final FeatureCodec<T, Reader<T>> codec ) {
            this.path = inputStreamName;
            this.codec = codec;
            this.indexFilePointer = findIndexFilePointer(ss);
            this.bcis = new BlockCompressedInputStream(new SeekableBufferedStream(ss));
            this.dis = new DataInputStream(bcis);
            this.header = readHeader();
            this.dataFilePointer = bcis.getPosition();
            this.index = null;
            this.usedByIterator = false;
        }

        public FeatureCodecHeader getFeatureCodecHeader() {
            return new FeatureCodecHeader(header, dataFilePointer);
        }
        public DataInputStream getStream() { return dis; }
        public List<String> getSampleNames() { return header.getSampleNames(); }
        public SAMSequenceDictionary getDictionary() { return header.getDictionary(); }
        public String getVersion() { return header.getVersion(); }

        public boolean hasNext() {
            final long position = bcis.getPosition();
            // A BlockCompressedInputStream returns a 0 position when closed
            return position > 0 && position < indexFilePointer;
        }

        @Override
        public CloseableTribbleIterator<T> query( final String chr, final int start, final int end )
                throws IOException {
            if ( index == null ) {
                loadIndex(bcis);
            }
            final SVInterval interval =
                    new SVInterval(getDictionary().getSequenceIndex(chr), start, end);
            return new OverlapIterator<>(interval, this);
        }

        @Override public CloseableTribbleIterator<T> iterator() {
            return new CompleteIterator<>(this);
        }

        @Override public void close() {
            try {
                dis.close();
            } catch ( final IOException ioe ) {
                throw new UserException("unable to close " + path, ioe);
            }
        }

        @Override public List<String> getSequenceNames() {
            final SAMSequenceDictionary dict = getDictionary();
            final List<String> names = new ArrayList<>(dict.size());
            for ( final SAMSequenceRecord rec : dict.getSequences() ) {
                names.add(rec.getSequenceName());
            }
            return names;
        }

        @Override public Object getHeader() { return header; }

        @Override public boolean isQueryable() { return true; }

        public long getPosition() { return bcis.getPosition(); }

        public void seekStream( final long filePointer ) {
            try {
                bcis.seek(filePointer);
            } catch ( final IOException ioe ) {
                throw new UserException("unable to position stream for " + path, ioe);
            }
        }

        public T readStream() {
            try {
                return codec.decode(this);
            } catch ( final IOException ioe ) {
                throw new GATKException("can't read " + path, ioe);
            }
        }

        private long findIndexFilePointer( final SeekableStream ss ) {
            final int finalBlockLen = EMPTY_GZIP_BLOCK_WITH_INDEX_POINTER.length;
            final byte[] finalBlock = new byte[finalBlockLen];
            try {
                ss.seek(ss.length() - finalBlockLen);
                ss.readFully(finalBlock);
                ss.seek(0);
            } catch ( final IOException ioe ) {
                throw new UserException("unable to read final bgzip block from " + path, ioe);
            }
            for ( int idx = 0; idx != FILE_POINTER_OFFSET; ++idx ) {
                if ( EMPTY_GZIP_BLOCK_WITH_INDEX_POINTER[idx] != finalBlock[idx] ) {
                    throw new UserException(
                            "unable to recover index pointer from final block of " + path);
                }
            }
            for ( int idx = FILE_POINTER_OFFSET + 8; idx != finalBlockLen; ++idx ) {
                if ( EMPTY_GZIP_BLOCK_WITH_INDEX_POINTER[idx] != finalBlock[idx] ) {
                    throw new UserException(
                            "unable to recover index pointer from final block of " + path);
                }
            }
            long indexFilePointer = 0;
            int idx = FILE_POINTER_OFFSET + 8;
            while ( --idx >= FILE_POINTER_OFFSET ) {
                indexFilePointer <<= 8;
                indexFilePointer |= finalBlock[idx] & 0xFFL;
            }
            return indexFilePointer;
        }

        private SVFeaturesHeader readHeader() {
            try {
                final String className = dis.readUTF();
                final String version = dis.readUTF();
                final List<String> sampleNames = readSampleNames(dis);
                final SAMSequenceDictionary dictionary = readDictionary(dis);
                return new SVFeaturesHeader(className, version, dictionary, sampleNames);
            } catch ( final IOException ioe ) {
                throw new UserException("can't read header from " + path, ioe);
            }
        }

        private List<String> readSampleNames( final DataInputStream dis ) throws IOException {
            final int nSamples = dis.readInt();
            final List<String> sampleNames = new ArrayList<>(nSamples);
            for ( int sampleId = 0; sampleId < nSamples; ++sampleId ) {
                sampleNames.add(dis.readUTF());
            }
            return sampleNames;
        }

        private SAMSequenceDictionary readDictionary( final DataInputStream dis ) throws IOException {
            final int nRecs = dis.readInt();
            final List<SAMSequenceRecord> seqRecs = new ArrayList<>(nRecs);
            for ( int idx = 0; idx != nRecs; ++idx ) {
                final int contigSize = dis.readInt();
                final String contigName = dis.readUTF();
                seqRecs.add(new SAMSequenceRecord(contigName, contigSize));
            }
            return new SAMSequenceDictionary(seqRecs);
        }

        private void loadIndex( final BlockCompressedInputStream bcis ) {
            final SVIntervalTree<Long> intervalTree = new SVIntervalTree<>();
            try {
                bcis.seek(indexFilePointer);
                final DataInputStream dis = new DataInputStream(bcis);
                int nEntries = dis.readInt();
                while ( nEntries-- > 0 ) {
                    final IndexEntry entry = new IndexEntry(dis);
                    intervalTree.put(entry.getInterval(), entry.getFilePosition());
                }
                bcis.seek(dataFilePointer);
            } catch ( final IOException ioe ) {
                throw new UserException("unable to read index from " + path, ioe);
            }
            index = intervalTree;
        }

        private Reader<T> getReaderForIterator() {
            if ( !usedByIterator ) {
                usedByIterator = true;
                return this;
            }
            return new Reader<>(this);
        }

        private static class CompleteIterator <T extends Feature>
                implements CloseableTribbleIterator<T> {
            final Reader<T> reader;
            public CompleteIterator( final Reader<T> reader ) {
                this.reader = reader.getReaderForIterator();
                reader.seekStream(reader.dataFilePointer);
            }

            @Override public Iterator<T> iterator() {
                return new CompleteIterator<>(reader);
            }

            @Override public boolean hasNext() {
                return reader.hasNext();
            }

            @Override public T next() {
                if ( !hasNext() ) {
                    throw new NoSuchElementException("feature iterator has no next element");
                }
                return reader.readStream();
            }

            @Override public void close() { reader.close(); }
        }
        // find all the objects in the stream inflating just those blocks that might have relevant objects
        private static class OverlapIterator <T extends Feature>
                implements CloseableTribbleIterator<T> {
            final SVInterval interval;
            final Reader<T> reader;
            final Iterator<SVIntervalTree.Entry<Long>> indexEntryIterator;
            long blockStartPosition;
            T nextT;

            public OverlapIterator( final SVInterval interval, final Reader<T> reader ) {
                this.interval = interval;
                this.reader = reader.getReaderForIterator();
                this.indexEntryIterator = reader.index.overlappers(interval);
                this.blockStartPosition = -1;
                advance();
            }

            @Override public boolean hasNext() { return nextT != null; }

            @Override public T next() {
                final T result = nextT;
                if ( result == null ) {
                    throw new NoSuchElementException("overlapper iterator has no next element");
                }
                advance();
                return result; }

            @Override public void close() { reader.close(); nextT = null; }

            @Override public CloseableTribbleIterator<T> iterator() {
                return new OverlapIterator<>(interval, reader);
            }

            private void advance() {
                do {
                    if ( isNewBlock(blockStartPosition, reader.getPosition()) ) {
                        if ( !indexEntryIterator.hasNext() ) {
                            nextT = null;
                            return;
                        }
                        blockStartPosition = indexEntryIterator.next().getValue();
                        reader.seekStream(blockStartPosition);
                    }
                    nextT = reader.readStream();
                    final int nextContigId =
                            reader.getDictionary().getSequenceIndex(nextT.getContig());
                    if ( interval.getContig() != nextContigId ||
                            interval.getEnd() < nextT.getStart() ) {
                        nextT = null;
                        return;
                    }
                } while ( nextT.getEnd() < interval.getStart() );
            }
        }
    }

    public static boolean isNewBlock( final long filePosition1, final long filePosition2 ) {
        // upper 48 bits contain the block offset
        // check to see if there are any bit differences in those upper 48 bits
        return ((filePosition1 ^ filePosition2) & ~0xffffL) != 0;
    }
}
