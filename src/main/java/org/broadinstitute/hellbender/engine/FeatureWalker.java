package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBOptions;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * A FeatureWalker is a tool that processes a {@link Feature} at a time from a source of Features, with
 * optional contextual information from a reference, sets of reads, and/or supplementary sources
 * of Features.

 * Subclasses must implement the {@link #apply(Feature, ReadsContext, ReferenceContext, FeatureContext)} method to process each Feature,
 * as well as {@link #isAcceptableFeatureType(Class)} and {@link #getDrivingFeaturePath()}, and may optionally implement
 * {@link #onTraversalStart()}, {@link #onTraversalSuccess()}, and/or {@link #closeTool()}.
 *
 * @param <F> the driving feature type.
 */
public abstract class FeatureWalker<F extends Feature> extends WalkerBase {

    private FeatureDataSource<F> drivingFeatures;
    private Object header;

    @Override
    public boolean requiresFeatures(){
        return true;
    }
    
    @Override
    public String getProgressMeterRecordLabel() { return "features"; }

    @Override
    void initializeFeatures() {
        features = new FeatureManager(this, FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
                                      getGenomicsDBOptions());
        initializeDrivingFeatures();
    }

    /**
     * Set the intervals for traversal in the driving features.
     *
     * Marked final so that subclasses don't override it. Subclasses should override {@link #onTraversalStart} instead.
     */
    @Override
    protected final void onStartup() {
        super.onStartup();
        // set the intervals for the feature here, because they are not initialized when initialize features is set
        if ( hasUserSuppliedIntervals() ) {
            drivingFeatures.setIntervalsForTraversal(userIntervals);
        }
    }

    @SuppressWarnings("unchecked")
    private void initializeDrivingFeatures() {
        final GATKPath drivingPath = getDrivingFeaturePath();
        final FeatureCodec<? extends Feature, ?> codec = FeatureManager.getCodecForFile(drivingPath.toPath());
        final Class<? extends Feature> featureType = codec.getFeatureType();
        if (isAcceptableFeatureType(featureType)) {
            final GenomicsDBOptions options = new GenomicsDBOptions(referenceArguments.getReferencePath());
            final FeatureInput<F> drivingFeatureInput = new FeatureInput<>(drivingPath);
            drivingFeatureInput.setFeatureCodecClass((Class<FeatureCodec<F, ?>>)codec.getClass());
            drivingFeatures =
                    new FeatureDataSource<>(drivingFeatureInput, FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES,
                            featureType, cloudPrefetchBuffer, cloudIndexPrefetchBuffer, options, false);
            header = drivingFeatures.getHeader();

            final FeatureInput<F> featureInput = new FeatureInput<>(drivingPath, "drivingFeatureFile");
            featureInput.setFeatureCodecClass((Class<FeatureCodec<F, ?>>)codec.getClass());
            features.addToFeatureSources(featureInput,
                    new FeatureDataSource<>(featureInput, 0, featureType,
                        cloudPrefetchBuffer, cloudIndexPrefetchBuffer, options, false));
        } else {
            throw new UserException("File " + drivingPath.getRawInputString() + " contains features of the wrong type.");
        }
    }

    /**
     * Returns whether the given class of features is acceptable for this walker.
     */
    protected abstract boolean isAcceptableFeatureType(Class<? extends Feature> featureType);

    /**
     * {@inheritDoc}
     *
     * Implementation of Feature-based traversal.
     *
     * NOTE: You should only override {@link #traverse()} if you are writing a new walker base class in the
     * engine package that extends this class. It is not meant to be overridden by tools outside of the engine
     * package.
     */
    @Override
    public void traverse() {
        CountingReadFilter readFilter = makeReadFilter();
        // Process each feature in the input stream.
        Utils.stream(drivingFeatures).forEach(feature -> {
                    final SimpleInterval featureInterval = makeFeatureInterval(feature);
                    apply(feature,
                            new ReadsContext(reads, featureInterval, readFilter),
                            new ReferenceContext(reference, featureInterval),
                            new FeatureContext(features, featureInterval));
                    progressMeter.update(feature);
                });
    }

    /**
     * This method can be overridden if you need to customize the interval for a given feature.
     *
     * @param feature {@link Feature} to derive the interval.
     * @param <T> Class that extends feature.
     * @return Interval for the given feature.  Typically, this is just the extents of the feature itself.
     * Never {@code null}
     */
    protected <T extends Feature> SimpleInterval makeFeatureInterval(final T feature) {
        return new SimpleInterval(feature);
    }

    /**
     * Process an individual feature.
     * In general, subclasses should simply stream their output from apply(), and maintain as little internal state
     * as possible.
     *
     * @param feature Current Feature being processed.
     * @param readsContext Reads overlapping the current feature. Will be an empty, but non-null, context object
     *                     if there is no backing source of reads data (in which case all queries on it will return
     *                     an empty array/iterator)
     * @param referenceContext Reference bases spanning the current feature. Will be an empty, but non-null, context object
     *                         if there is no backing source of reference data (in which case all queries on it will return
     *                         an empty array/iterator). Can request extra bases of context around the current feature's interval
     *                         by invoking {@link ReferenceContext#setWindow}
     *                         on this object before calling {@link ReferenceContext#getBases}
     * @param featureContext Features spanning the current feature. Will be an empty, but non-null, context object
     *                       if there is no backing source of Feature data (in which case all queries on it will return an
     *                       empty List).
     */
    public abstract void apply(final F feature, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext );

    /**
     * Close the reads and reference data sources.
     *
     * Marked final so that subclasses don't override it. Subclasses should override {@link #onTraversalSuccess()} instead.
     */
    @Override
    protected final void onShutdown() {
        super.onShutdown();

        if ( drivingFeatures != null ) {
            drivingFeatures.close();
        }
    }

    /**
     * Returns the file that contains the driving features.
     *
     * @return never {@code null}.
     */
    public abstract GATKPath getDrivingFeaturePath();


    /**
     * Returns the header of the driving features file.
     */
    public Object getDrivingFeaturesHeader() {
        return header;
    }
}
