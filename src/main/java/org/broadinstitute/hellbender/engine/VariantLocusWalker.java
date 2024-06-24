package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.CountingVariantFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.transformers.VariantTransformer;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.iterators.IntervalLocusIterator;
import org.broadinstitute.hellbender.utils.iterators.ShardedIntervalIterator;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * VariantLocusWalker processes variants from a single source, grouped by locus overlap, or optionally one
 * at a time in order, with optional contextual information from a reference, sets of reads, and/or supplementary sources
 * of Features. By-variant traversal is opt in, via {@link #changeTraversalModeToByVariant()}, otherwise only loci with
 * overlapping variants are traversed.
 *
 * VariantLocusWalker authors must implement the {@link #apply} method to process each variant, and may optionally implement
 * {@link #onTraversalStart}, {@link #onTraversalSuccess} and/or {@link #closeTool}.
 */
public abstract class VariantLocusWalker extends VariantWalkerBase {

    // NOTE: using String rather than FeatureInput<VariantContext> here so that we can keep this driving source
    //       of variants separate from any other potential sources of Features
    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "A VCF file containing variants", common = false, optional = false)
    public String drivingVariantFile;

    // NOTE: keeping the driving source of variants separate from other, supplementary FeatureInputs in our FeatureManager in GATKTool
    //we do add the driving source to the Feature manager but we do need to treat it differently and thus this field.
    private FeatureDataSource<VariantContext> drivingVariants;
    private FeatureInput<VariantContext> drivingVariantsFeatureInput;

    private boolean traverseByVariant = false;

    /**
     * Traverse by variant instead of by grouped by locus. When true, the {@link #apply} method will be called for each
     * individual variant in the driving variants input, instead of by locus, with all variants overlapping that locus.
     * This mehtod must be called in OnTraversalStart();
     */
    protected void changeTraversalModeToByVariant() { traverseByVariant = true; }

    @Override
    protected SAMSequenceDictionary getSequenceDictionaryForDrivingVariants() { return drivingVariants.getSequenceDictionary(); }

    @Override
    protected Spliterator<VariantContext> getSpliteratorForDrivingVariants() { return drivingVariants.spliterator(); }

    /**
     * Marked final so that tool authors don't override it. Tool authors should override {@link #onTraversalStart} instead.
     */
    @Override
    protected final void onStartup() {
        super.onStartup();
        if ( hasUserSuppliedIntervals() ) {
            drivingVariants.setIntervalsForTraversal(userIntervals);
        }
    }

    @Override
    protected void initializeDrivingVariants() {
        drivingVariantsFeatureInput = new FeatureInput<>(drivingVariantFile, "drivingVariantFile");

        // This is the data source for the driving source of variants,
        // which uses a cache lookahead of getDrivingVariantCacheLookAheadBases()
        drivingVariants = new FeatureDataSource<>(drivingVariantsFeatureInput, getDrivingVariantCacheLookAheadBases(), VariantContext.class, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
                getGenomicsDBOptions());

        // Also add the driving datasource to the feature manager so that it can be queried. Setting cache lookahead
        // to 0 to avoid caching. Note: we are disabling lookahead here because of windowed queries that need to "look behind" as well.
        features.addToFeatureSources(0, drivingVariantsFeatureInput, VariantContext.class, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
                getGenomicsDBOptions());

        // Note: the intervals for the driving variants are set in onStartup()
    }

    /**
     * Returns the feature input for the driving variants file.
     */
    protected final FeatureInput<VariantContext> getDrivingVariantsFeatureInput() {
        return drivingVariantsFeatureInput;
    }

    /**
     * Gets the header associated with our driving source of variants as a VCFHeader.
     *
     * @return VCFHeader for our driving source of variants
     */
    public final VCFHeader getHeaderForVariants() {
        final Object header = drivingVariants.getHeader();

        if ( ! (header instanceof VCFHeader) ) {
            throw new GATKException("Header for " + drivingVariantFile + " is not in VCF header format");
        }

        return (VCFHeader)header;
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    /**
     * {@inheritDoc}
     *
     * Implementation of variant-locus traversal.
     *
     * NOTE: You should only override {@link #traverse()} if you are writing a new walker base class in the
     * engine package that extends this class. It is not meant to be overridden by tools outside of the engine
     * package.
     */
    @Override
    public void traverse() {
        final CountingReadFilter readFilter = makeReadFilter();
        final CountingVariantFilter variantFilter = makeVariantFilter();
        final VariantTransformer preTransformer  = makePreVariantFilterTransformer();
        final VariantTransformer postTransformer = makePostVariantFilterTransformer();

        if (traverseByVariant) {
            // Process each variant in the input stream, one at a time.
            getTransformedVariantStream( getSpliteratorForDrivingVariants(), preTransformer, variantFilter, postTransformer )
                    .forEachOrdered(variant -> {

                        final SimpleInterval variantInterval = new SimpleInterval(variant);

                        try{
                            apply(variant,
                                    Collections.singletonList(variant),
                                    new ReadsContext(reads, variantInterval, readFilter),
                                    new ReferenceContext(reference, variantInterval),
                                    new FeatureContext(features, variantInterval));
                        } catch (final IllegalStateException e) {
                            throw new GATKException("Exception thrown at " + variant.getContig() + ":" + variant.getStart()
                                    + " " + variant.toString(), e);
                        }

                        progressMeter.update(variantInterval);
                    });
        } else {
            // Traverse loci in shards. For any shard with overlapping variants, drop down to per-locus iteration,
            // calling apply for a single locus, only if there are overlapping variants, passing all such variants
            // as a group.
            Utils.stream(new ShardedIntervalIterator(getTraversalIntervals().iterator(), getDrivingVariantCacheLookAheadBases()))
                    .forEachOrdered (shard -> {
                        if (drivingVariants.query(shard).hasNext()) {
                            getLocusStream(new SimpleInterval(shard))
                                    .forEachOrdered(locus -> {
                                        final Iterator<VariantContext> overlappingVariants = drivingVariants.query(locus);
                                        if (overlappingVariants.hasNext()) {
                                            final List<VariantContext> filteredVariants = getTransformedVariantStream(
                                                    Spliterators.spliteratorUnknownSize(overlappingVariants, 0),
                                                    preTransformer,
                                                    variantFilter,
                                                    postTransformer)
                                                    .collect(Collectors.toList());
                                            if (!filteredVariants.isEmpty()) {
                                                try {
                                                    apply(locus,
                                                            filteredVariants,
                                                            new ReadsContext(reads, locus, readFilter),
                                                            new ReferenceContext(reference, locus),
                                                            new FeatureContext(features, locus));
                                                } catch (final IllegalStateException e) {
                                                    throw new GATKException("Exception thrown at first variant start " + filteredVariants.get(0).getContig() + ":" + filteredVariants.get(0).getStart()
                                                            + " " + filteredVariants.get(0).toString(), e);
                                                }

                                                progressMeter.update(locus);
                                            }
                                        }
                                    });
                        }
                    });
        }
    }

    // Return a Stream of SimpleInterval covering the entire territory sketched out by requestedInterval
    private Stream<SimpleInterval> getLocusStream(final SimpleInterval requestedInterval) {
        return Utils.stream(new IntervalLocusIterator(Collections.singletonList(requestedInterval).iterator()));
    }

    /**
     * Process by locus, with all variants overlapping the current locus, or by individual variant when in by-variants
     * traversal mode. Must be implemented by tool authors. In general, tool authors should simply stream their output
     * from apply(), and maintain as little internal state as possible.
     *
     * @param loc the current locus being traversed (in group by locus), or the span of the current variant (in by-variants
     *            traversal)
     * @param variants The current variant(s) being processed. NOTE: When {@link #changeTraversalModeToByVariant()}
     *                 has been called, the {@link #apply} method will be called for each individual variant in the
     *                 driving variants input. Otherwise, this method is called once for every locus that has
     *                 overlapping variants, with all variants overlapping that locus.
     *
     * @param readsContext Reads overlapping the current variant. Will be an empty, but non-null, context object
     *                     if there is no backing source of reads data (in which case all queries on it will return
     *                     an empty array/iterator)
     * @param referenceContext Reference bases spanning the current variant. Will be an empty, but non-null, context object
     *                         if there is no backing source of reference data (in which case all queries on it will return
     *                         an empty array/iterator). Can request extra bases of context around the current variant's interval
     *                         by invoking {@link ReferenceContext#setWindow}
     *                         on this object before calling {@link ReferenceContext#getBases}
     * @param featureContext Features spanning the current variant. Will be an empty, but non-null, context object
     *                       if there is no backing source of Feature data (in which case all queries on it will return an
     *                       empty List).
     */
    public abstract void apply(Locatable loc, List<VariantContext> variants, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext );

    /**
     * Close all data sources.
     *
     * Marked final so that tool authors don't override it. Tool authors should override {@link #onTraversalSuccess} and/or
     * {@link #closeTool} instead.
     */
    @Override
    protected final void onShutdown() {
        super.onShutdown();

        if ( drivingVariants != null )
            drivingVariants.close();
    }
}
