package org.broadinstitute.hellbender.engine;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.validation.ConcordanceState;

/**
 * Base class for concordance walkers which define concordance as agreement between records on the ref and first alt allele.
 *
 * Subclasses must implement the {@link #apply} method to process each variant
 * and may optionally implement {@link #onTraversalStart} and/or {@link #onTraversalSuccess}.
 */
public abstract class AlleleConcordanceWalker extends AbstractConcordanceWalker {

    @Override
    protected void apply(final TruthVersusEval truthVersusEval, final ReadsContext readsContext, final ReferenceContext refContext) {
        apply(convertToAlleleTruthVersusEval(truthVersusEval), readsContext, refContext);
    }

    // the primary work of the walker.  Must be overridden in implementing classes.
    protected abstract void apply(final AlleleTruthVersusEval truthVersusEval, final ReadsContext readsContext, final ReferenceContext refContext);

    /**
     * store a truth vc in case of a false negative, an eval vc in case of a false positive, or a concordance pair of
     * truth and eval in case of a true positive.
     */
    protected static class AlleleTruthVersusEval extends TruthVersusEval {
        private final ConcordanceState concordanceState;

        public AlleleTruthVersusEval(final VariantContext truth, final VariantContext eval,
                                     final ConcordanceState concordanceState) {
            super(truth, eval);
            this.concordanceState = concordanceState;
        }

        public ConcordanceState getConcordance() {
            return concordanceState;
        }

    }

    protected AlleleTruthVersusEval convertToAlleleTruthVersusEval(final TruthVersusEval truthVersusEval) {
        if (!truthVersusEval.hasTruth()) {
            if (truthVersusEval.getEval().isFiltered()) {
                return new AlleleTruthVersusEval(truthVersusEval.getTruth(), truthVersusEval.getEval(), ConcordanceState.FILTERED_TRUE_NEGATIVE);
            } else {
                return new AlleleTruthVersusEval(truthVersusEval.getTruth(), truthVersusEval.getEval(), ConcordanceState.FALSE_POSITIVE);
            }
        }

        if (!truthVersusEval.hasEval()) {
            return new AlleleTruthVersusEval(truthVersusEval.getTruth(), truthVersusEval.getEval(), ConcordanceState.FALSE_NEGATIVE);
        }

        if (truthVersusEval.getEval().isFiltered()) {
            return new AlleleTruthVersusEval(truthVersusEval.getTruth(), truthVersusEval.getEval(), ConcordanceState.FILTERED_FALSE_NEGATIVE);
        } else {
            return new AlleleTruthVersusEval(truthVersusEval.getTruth(), truthVersusEval.getEval(), ConcordanceState.TRUE_POSITIVE);
        }
    }

    protected boolean doRefAndAltAllelesMatch(final VariantContext truth, final VariantContext eval) {
        final boolean sameRefAllele = truth.getReference().equals(eval.getReference());

        // We make sure that the truth has at least one alt allele.
        // If it does, we pick the first for comparison:
        final boolean containsAltAllele =
                (truth.getAlternateAlleles().size() == eval.getAlternateAlleles().size()) &&
                        ((truth.getAlternateAlleles().size() > 0) &&
                                eval.getAlternateAlleles().contains(truth.getAlternateAllele(0)));

        return sameRefAllele && containsAltAllele;
    }

    @Override
    protected boolean shouldVariantsBeMatched(final VariantContext truth, final VariantContext eval) {
        return doRefAndAltAllelesMatch(truth, eval);
    }
}
