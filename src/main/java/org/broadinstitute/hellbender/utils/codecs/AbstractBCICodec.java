package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.samtools.util.LocationAware;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.FeatureCodecHeader;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Reader;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Writer;

import java.io.IOException;
import java.io.InputStream;

public abstract class AbstractBCICodec<F extends Feature>
        implements FeatureOutputCodec<F, Writer<F>>, FeatureCodec<F, Reader<F>>, FeatureReaderFactory<F> {

    @Override
    public Reader<F> getReader( final FeatureInput<F> input,
                                final int cloudPrefetchBufferSize,
                                final int cloudIndexPrefetchBufferSize ) {
        return new Reader<F>(input, this, cloudPrefetchBufferSize, cloudIndexPrefetchBufferSize);
    }

    @Override
    public Feature decodeLoc( final Reader<F> reader ) throws IOException {
        return decode(reader);
    }

    @Override
    public FeatureCodecHeader readHeader( final Reader<F> reader ) throws IOException {
        return reader.getFeatureCodecHeader();
    }

    @Override
    public Reader<F> makeSourceFromStream( final InputStream is ) {
        throw new GATKException("unimplemented method");
    }

    @Override
    public LocationAware makeIndexableSourceFromStream( final InputStream is ) {
        throw new GATKException("unimplemented method");
    }

    @Override
    public boolean isDone( final Reader<F> reader ) { return !reader.hasNext(); }

    @Override
    public void close( final Reader<F> reader ) { reader.close(); }
}
