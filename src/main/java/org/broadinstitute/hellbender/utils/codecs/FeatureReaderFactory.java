package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureReader;
import org.broadinstitute.hellbender.engine.FeatureInput;

/** Something that can produce a FeatureReader */
public interface FeatureReaderFactory<F extends Feature> {
    FeatureReader<F> getReader( FeatureInput<F> path,
                                int cloudPrefetchBufferSize,
                                int cloudIndexPrefetchBufferSize );
}
