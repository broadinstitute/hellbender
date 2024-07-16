package org.broadinstitute.hellbender.engine;

import htsjdk.beta.io.bundle.Bundle;
import htsjdk.beta.io.bundle.BundleJSON;
import htsjdk.beta.io.bundle.BundleResource;
import htsjdk.beta.io.bundle.BundleResourceType;
import htsjdk.beta.plugin.IOUtils;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MultiVariantInputArgumentCollection;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.ArrayList;
import java.util.List;
import java.util.SortedSet;
import java.util.Spliterator;

/**
 * A MultiVariantWalker is a tool that processes one variant at a time, in position order, from multiple sources of
 * variants, with optional contextual information from a reference, sets of reads, and/or supplementary sources of
 * Features.
 *
 * VariantWalker authors must implement the {@link #apply} method to process each variant, and may optionally implement
 * {@link #onTraversalStart}, {@link #onTraversalSuccess} and/or {@link #closeTool}.
 */
public abstract class MultiVariantWalker extends VariantWalkerBase {

    @ArgumentCollection
    protected MultiVariantInputArgumentCollection multiVariantInputArgumentCollection = getMultiVariantInputArgumentCollection();

    // NOTE: keeping the driving source of variants separate from other, supplementary FeatureInputs in our FeatureManager
    // in GATKTool we do add the driving source to the Feature manager but we do need to treat it differently and thus this
    // field.
    private MultiVariantDataSource drivingVariants;
    private List<FeatureInput<VariantContext>> drivingVariantsFeatureInputs = new ArrayList<>(2);

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

    /**
     * Return an argument collection that provides the driving variants.  This allows subclasses to override and use a different
     * argument pattern besides the default -V
     */
    protected MultiVariantInputArgumentCollection getMultiVariantInputArgumentCollection() {
        return new MultiVariantInputArgumentCollection.DefaultMultiVariantInputArgumentCollection();
    }

    /**
     * Should this tool cross-validate the sequence dictionaries for input feature files? If a large number of inputs are expected,
     * the tools should override to return false and require a reference.
     *
     * @return true if this tool should check all input feature files' contigs against each other, otherwise false
     */
    public boolean doDictionaryCrossValidation() {
        return true;
    }

    @Override
    protected void initializeDrivingVariants() {
        multiVariantInputArgumentCollection.getDrivingVariantPaths().stream().forEach(
                gatkPath -> {
                    if (gatkPath.hasExtension(BundleJSON.BUNDLE_EXTENSION)) {
                        // expand any bundle(s) into one or more FeatureInputs (depending on whether it is a single
                        // bundle or list)
                        final List<Bundle> bundles = BundleJSON.toBundleList(IOUtils.getStringFromPath(gatkPath), GATKPath::new);
                        for (final Bundle bundle : bundles) {
                            if (bundle.getPrimaryContentType().equals(BundleResourceType.CT_VARIANT_CONTEXTS)) {
                                // use the bundle primary resource as the FeatureInput URI, and tear off and attach the
                                // individual bundle the bundle to the FI as the parent bundle so downstream code can
                                // extract other resources from it on demand
                                // note that if the original value from the user has a tag, we can't use it unless there
                                // is only one input, since FIs have to be unique
                                final FeatureInput<VariantContext> bundleFI = new FeatureInput<>(
                                        new GATKPath(bundle.getPrimaryResource().getIOPath().get().getURIString()),
                                        bundle,
                                        bundles.size() > 1 ? gatkPath.getTag() : "drivingVariants"
                                );
                                if (drivingVariantsFeatureInputs.contains(bundleFI)) {
                                    throw new UserException.BadInput("Feature inputs must be unique: " + gatkPath);
                                }
                                drivingVariantsFeatureInputs.add(bundleFI);
                                // Add each driving variants FeatureInput to the feature manager so that it can be queried, using a lookahead value
                                // of 0 to avoid caching because of windowed queries that need to "look behind" as well.
                                features.addToFeatureSources(0, bundleFI, VariantContext.class, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
                                        referenceArguments.getReferencePath());
                            } else {
                                final BundleResource br = bundle.getPrimaryResource();
                                throw new UserException.BadInput(
                                        String.format(
                                                "Feature input bundles must have a primary resource of type %s. Found content type %s with path %s",
                                                BundleResourceType.CT_VARIANT_CONTEXTS,
                                                bundle.getPrimaryContentType(),
                                                br.getIOPath().get()));
                            }
                        }
                    } else {
                        final FeatureInput<VariantContext> featureInput = new FeatureInput<>(gatkPath);
                        drivingVariantsFeatureInputs.add(featureInput);
                        // Add each driving variants FeatureInput to the feature manager so that it can be queried, using a lookahead value
                        // of 0 to avoid caching because of windowed queries that need to "look behind" as well.
                        features.addToFeatureSources(0, featureInput, VariantContext.class, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
                                referenceArguments.getReferencePath());
                    }
                }
        );

        final boolean skipDictionaryValidation = !seqValidationArguments.performSequenceDictionaryValidation() || !doDictionaryCrossValidation();
        if (skipDictionaryValidation && !hasReference() && getMasterSequenceDictionary() == null) {
            logger.warn("Input feature file dictionaries cannot be validated.  For safety, please specify a reference or " +
                    "a master sequence dictionary with `--" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME + "`");
        }

        // Create a (MultiVariantDataSource) FeatureDataSource for the driving variants inputs using the
        // cache lookahead value from getDrivingVariantCacheLookAheadBases()
        drivingVariants = new MultiVariantDataSource(drivingVariantsFeatureInputs, getDrivingVariantCacheLookAheadBases(), cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
                                                     referenceArguments.getReferencePath(), skipDictionaryValidation);

        // Note: the intervals for the driving variants are set in onStartup()
    }

    /**
     * Returns a list of feature inputs used for the driving variants for this source.
     */
    protected final List<FeatureInput<VariantContext>> getDrivingVariantsFeatureInputs() {
        return drivingVariantsFeatureInputs;
    }

    /**
     * Gets the header associated with our driving source of variants as a VCFHeader.
     *
     * @return VCFHeader for our driving source of variants
     */
    @Override
    public final VCFHeader getHeaderForVariants() {
        return drivingVariants.getHeader();
    }

    /**
     * {@inheritDoc}
     *
     * Implementation of multi-variant traversal.
     *
     * NOTE: You should only override {@link #traverse()} if you are writing a new walker base class in the
     * engine package that extends this class. It is not meant to be overridden by tools outside of the
     * engine package.
     */
    @Override
    public void traverse() {
        final CountingReadFilter readFilter = makeReadFilter();
        // Process each variant in the input stream.
        getTransformedVariantStream( makeVariantFilter() )
                .forEach(variant -> {
                    final SimpleInterval variantInterval = new SimpleInterval(variant);
                    try {
                        apply(variant,
                                new ReadsContext(reads, variantInterval, readFilter),
                                new ReferenceContext(reference, variantInterval),
                                new FeatureContext(features, variantInterval));
                    } catch (final Exception e) {
                        throw new GATKException("Exception thrown at " + variant.getContig() + ":" + variant.getStart()
                                + " " + variant.toString(), e);
                    }

                    progressMeter.update(variantInterval);
                });
    }

    /**
     * Process an individual variant. Must be implemented by tool authors.
     * In general, tool authors should simply stream their output from apply(), and maintain as little internal state
     * as possible.
     *
     * @param variant Current variant being processed.
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
    public abstract void apply( VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext );

    /**
     * Close all data sources.
     *
     * Marked final so that tool authors don't override it. Tool authors should override {@link #onTraversalSuccess} and/or
     * {@link #closeTool} instead.
     */
    @Override
    protected final void onShutdown() {
        super.onShutdown();
        if (drivingVariants != null) {
            drivingVariants.close();
        }
    }

    /**
     * Get in name-sorted order a list of samples merged from the driving variants files.
     *
     * @return SampleSet merged requiring unique names from the drivingVariants
     */
    public final SortedSet<String> getSamplesForVariants() {
        return drivingVariants.getSamples();
    }
}
