package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.samtools.util.LocationAware;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.io.TextFeatureReader;

import java.io.IOException;
import java.io.InputStream;
import java.util.Iterator;

/** Classes extending this one must implement:
 * F decode( Iterator<String> reader );
 * FeatureCodecHeader readHeader( Iterator<String> reader );
 * Class<F> getFeatureType();
 * boolean canDecode(final String path);
 */
public abstract class AbstractTextCodec<F extends Feature> implements FeatureCodec<F, Iterator<String>>,
                                                                        FeatureReaderFactory<F> {
    @Override
    public TextFeatureReader<F> getReader( final FeatureInput<F> input,
                                           final int cloudPrefetchBufferSize,
                                           final int cloudIndexPrefetchBufferSize ) {
        return new TextFeatureReader<>(input, this, cloudPrefetchBufferSize, cloudIndexPrefetchBufferSize);
    }

    @Override
    public Feature decodeLoc( final Iterator<String> itr ) throws IOException {
        return decode(itr);
    }

    @Override
    public Iterator<String> makeSourceFromStream( final InputStream is ) {
        throw new GATKException("unimplemented method");
    }

    @Override
    public LocationAware makeIndexableSourceFromStream( final InputStream is ) {
        throw new GATKException("unimplemented method");
    }

    @Override
    public boolean isDone( final Iterator<String> itr ) {
        return !itr.hasNext();
    }

    @Override
    public void close( final Iterator<String> itr ) {}
}
