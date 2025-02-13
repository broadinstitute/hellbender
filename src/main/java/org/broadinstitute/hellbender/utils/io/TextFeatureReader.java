package org.broadinstitute.hellbender.utils.io;

import htsjdk.samtools.seekablestream.ISeekableStreamFactory;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.*;
import htsjdk.tribble.index.Block;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.util.ParsingUtils;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.*;
import java.nio.channels.SeekableByteChannel;
import java.nio.file.Path;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.function.Function;

public final class TextFeatureReader<F extends Feature> implements FeatureReader<F> {
    private final Path path;
    private final String indexPath;
    private final FeatureCodec<F, Iterator<String>> codec;
    private final Function<SeekableByteChannel, SeekableByteChannel> featurePrefetchWrapper;
    private final Function<SeekableByteChannel, SeekableByteChannel> indexPrefetchWrapper;
    private FeatureCodecHeader header;
    private Boolean queryable;
    private Index index;

    private static final long BOGUS_FILE_OFFSET = -1L;

    public TextFeatureReader( final GATKPath path,
                              final FeatureCodec<F, Iterator<String>> codec,
                              final int cloudPrefetchBufferSize,
                              final int cloudIndexPrefetchBufferSize ) {
        this.path = path.toPath();
        final String indexPath = path.getTagAttributes().get("index");
        this.indexPath = indexPath != null ?
                indexPath :
                ParsingUtils.appendToPath(this.path.toUri().toString(), FileExtensions.TABIX_INDEX);
        this.codec = codec;
        final boolean prefetchable = BucketUtils.isEligibleForPrefetching(path);
        this.featurePrefetchWrapper =
                BucketUtils.getPrefetchingWrapper(prefetchable ? cloudPrefetchBufferSize : 0);
        this.indexPrefetchWrapper =
                BucketUtils.getPrefetchingWrapper(prefetchable ? cloudIndexPrefetchBufferSize : 0);
    }

    @Override
    public FeatureIterator<F> iterator() {
        return new FeatureIterator<>(path, codec, featurePrefetchWrapper);
    }

    @Override
    public boolean isQueryable() {
        if ( queryable == null ) {
            try {
                queryable = AbstractFeatureReader.isTabix(path.toUri().toString(), indexPath);
            } catch ( final IOException ioe ) {
                throw new UserException("Error while checking for existence of index for " + path, ioe);
            }
        }
        return queryable;
    }

    @Override
    public FeatureIterator<F> query( final String chr, final int start, final int end ) {
        if ( index == null ) {
            if ( !isQueryable() ) {
                throw new UserException("Querying non-queryable feature reader for " + path);
            }
            index = IndexFactory.loadIndex(indexPath, indexPrefetchWrapper);
            if ( index == null ) {
                throw new UserException("Despite its appearing to exist, couldn't load index for " + path);
            }
        }
        final List<Block> blocks = index.getBlocks(chr, start, end);
        return new FeatureIterator<>(path,
                codec,
                featurePrefetchWrapper,
                new SimpleInterval(chr, start, end),
                blocks == null || blocks.isEmpty() ? BOGUS_FILE_OFFSET : blocks.get(0).getStartPosition());
    }

    @Override
    public void close() {}

    @Override
    public List<String> getSequenceNames() {
        return Collections.emptyList();
    }

    @Override
    public Object getHeader() {
        if ( header == null ) {
            final LineIterator lineIterator = new LineIterator(path, featurePrefetchWrapper);
            try {
                header = codec.readHeader(lineIterator);
            } catch ( final IOException ioe ) {
                throw new UserException("Can't read header from " + path, ioe);
            }
            lineIterator.close();
        }
        return header.getHeaderValue();
    }

    private final static class LineIterator implements Iterator<String> {
        private final Path path;
        private final BufferedReader reader;
        private String nextLine;

        public LineIterator( final Path path,
                             final Function<SeekableByteChannel, SeekableByteChannel> wrapper ) {
            this.path = path;
            try {
                final ISeekableStreamFactory ssf = SeekableStreamFactory.getInstance();
                final SeekableStream ss =
                        ssf.getBufferedStream(ssf.getStreamFor(path.toUri().toString(), wrapper));
                if ( IOUtil.hasGzipFileExtension(path) ) {
                    reader = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(ss)));
                } else {
                    reader = new BufferedReader(new InputStreamReader(ss));
                }
            } catch ( final IOException ioe ) {
                throw new UserException("Can't open " + path, ioe);
            }
            advance();
        }

        public LineIterator( final Path path,
                             final Function<SeekableByteChannel, SeekableByteChannel> wrapper,
                             final long virtualFileOffset ) {
            this.path = path;
            try {
                final ISeekableStreamFactory ssf = SeekableStreamFactory.getInstance();
                final BlockCompressedInputStream bcis =
                        new BlockCompressedInputStream(
                                ssf.getBufferedStream(ssf.getStreamFor(path.toUri().toString(), wrapper)));
                bcis.seek(virtualFileOffset);
                reader = new BufferedReader(new InputStreamReader(bcis));
            } catch ( final IOException ioe ) {
                throw new UserException("Can't open seekable " + path, ioe);
            }
            advance();
        }

        @Override
        public boolean hasNext() {
            return nextLine != null;
        }

        @Override
        public String next() {
            if ( nextLine == null ) {
                throw new NoSuchElementException("LineIterator exhausted for " + path);
            }
            final String result = nextLine;
            advance();
            return result;
        }

        public void close() {
            nextLine = null;
            try {
                reader.close();
            } catch ( final IOException ioe ) {
                throw new UserException("Can't close " + path, ioe);
            }
        }

        private void advance() {
            try {
                nextLine = reader.readLine();
                if ( nextLine == null ) {
                    close();
                }
            } catch ( final IOException ioe ) {
                throw new UserException("Can't read " + path, ioe);
            }
        }
    }

    public final static class FeatureIterator<F extends Feature> implements CloseableTribbleIterator<F> {
        private final Path path;
        private final FeatureCodec<F, Iterator<String>> codec;
        private final Function<SeekableByteChannel, SeekableByteChannel> wrapper;
        private final SimpleInterval interval;
        private final LineIterator lineIterator;
        private F nextFeature;

        public FeatureIterator( final Path path,
                                final FeatureCodec<F, Iterator<String>> codec,
                                final Function<SeekableByteChannel, SeekableByteChannel> wrapper ) {
            this.path = path;
            this.codec = codec;
            this.wrapper = wrapper;
            this.interval = null;
            lineIterator = new LineIterator(path, wrapper);
            advance();
        }

        public FeatureIterator( final Path path,
                                final FeatureCodec<F, Iterator<String>> codec,
                                final Function<SeekableByteChannel, SeekableByteChannel> wrapper,
                                final SimpleInterval interval,
                                final long virtualFileOffset ) {
            this.path = path;
            this.codec = codec;
            this.wrapper = wrapper;
            this.interval = interval;
            if ( virtualFileOffset == BOGUS_FILE_OFFSET ) {
                lineIterator = null;
                nextFeature = null;
            } else {
                lineIterator = new LineIterator(path, wrapper, virtualFileOffset);
                advance();
            }
        }

        @Override
        public FeatureIterator<F> iterator() {
            return new FeatureIterator<>(path, codec, wrapper);
        }

        @Override
        public boolean hasNext() {
            return nextFeature != null;
        }

        @Override
        public F next() {
            if ( nextFeature == null ) {
                throw new NoSuchElementException("No more features in " + path);
            }
            final F result = nextFeature;
            advance();
            return result;
        }

        @Override
        public void close() {
            nextFeature = null;
            if ( lineIterator != null ) {
                lineIterator.close();
            }
        }

        private void advance() {
            try {
                do {
                    nextFeature = codec.decode(lineIterator);
                    if ( nextFeature == null ) {
                        close();
                    } else if ( interval == null ) {
                        break;
                    } else if ( !interval.getContig().equals(nextFeature.getContig()) ||
                                interval.getEnd() < nextFeature.getStart() ) {
                        close();
                    } else if ( nextFeature.getEnd() >= interval.getStart() ) {
                        break;
                    }
                } while ( nextFeature != null );
            } catch ( final IOException ioe ) {
                throw new UserException("Can't read feature from " + path, ioe);
            }
        }
    }
}
