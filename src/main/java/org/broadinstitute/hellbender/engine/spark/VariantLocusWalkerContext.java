package org.broadinstitute.hellbender.engine.spark;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;

import java.util.List;

/**
 * Encapsulates a collection of {@link VariantContext} objects that overlap a given locus, with
 * its {@link ReferenceContext} and {@link FeatureContext}.
 */
public class VariantLocusWalkerContext {
    private Locatable locus;
    private final List<VariantContext> variants;
    private final ReadsContext readsContext;
    private final ReferenceContext referenceContext;
    private final FeatureContext featureContext;

    public VariantLocusWalkerContext(Locatable locus, List<VariantContext> variants, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        this.locus = locus;
        this.variants = variants;
        this.readsContext = readsContext;
        this.referenceContext = referenceContext;
        this.featureContext = featureContext;
    }

    public Locatable getLocus() {
        return locus;
    }

    public List<VariantContext> getVariants() {
        return variants;
    }

    public ReadsContext getReadsContext() {
        return readsContext;
    }

    public ReferenceContext getReferenceContext() {
        return referenceContext;
    }

    public FeatureContext getFeatureContext() {
        return featureContext;
    }
}
