package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.collections4.Predicate;
import org.apache.commons.collections4.iterators.FilterIterator;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.validation.ConcordanceState;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.File;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Set;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.stream.StreamSupport;

/**
 * Base class for concordance walkers, which process one variant at a time from one or more sources of variants,
 * with optional contextual information from a reference, sets of reads, and/or supplementary sources of
 * Features.
 *
 * Subclasses must implement the {@link #apply} and {@link #shouldVariantsBeMatched} methods to process each variant
 * and may optionally implement {@link #onTraversalStart} and/or {@link #onTraversalSuccess}.
 *
 * Created by Takuto Sato 1/30/17, abstractified by David Benjamin on 2/22/17.
 * {@link #onTraversalStart}, {@link #onTraversalSuccess} and/or {@link #closeTool}.
 */
public abstract class AbstractConcordanceWalker extends WalkerBase {

    public static final String TRUTH_VARIANTS_LONG_NAME = "truth";
    public static final String EVAL_VARIANTS_SHORT_NAME = "eval";
    public static final String EVAL_VARIANTS_LONG_NAME = "evaluation";

    public static final String CONFIDENCE_REGION_LONG_NAME = "confidence";
    public static final String CONFIDENCE_REGION_SHORT_NAME = "C";
    
    // The distance in bases to look ahead and cache when querying feature sources.
    public static final int CACHE_LOOKAHEAD = 100_000;

    @Argument(shortName = TRUTH_VARIANTS_LONG_NAME, fullName = TRUTH_VARIANTS_LONG_NAME,
            doc = "A VCF containing truth variants", optional = false)
    public String truthVariantsFile;

    @Argument(shortName = EVAL_VARIANTS_SHORT_NAME, fullName = EVAL_VARIANTS_LONG_NAME,
            doc = "A VCF containing variants to be compared to the truth", optional = false)
    public String evalVariantsFile;

    @Argument(doc = "TO BE IMPLEMENTED",
            fullName= CONFIDENCE_REGION_LONG_NAME,
            shortName = CONFIDENCE_REGION_SHORT_NAME,
            optional = true)
    protected File highConfidenceRegion;

    private FeatureDataSource<VariantContext> truthVariants;
    private FeatureDataSource<VariantContext> evalVariants;
    protected SAMSequenceDictionary dict;
    private VariantContextComparator variantContextComparator;

    // for indices, false -> 0, true -> 1
    private ConcordanceState[][][][] concordanceStatesTable =
            {// isFiltered = false
                {// genotypedAgree = false
                    {//truth-
                        //eval- , eval+
                        {null, ConcordanceState.FALSE_POSITIVE},
                     //truth+
                        //eval - , eval +
                        {ConcordanceState.FALSE_NEGATIVE, ConcordanceState.FALSE_POSITIVE}
                    },
                 // genotypesAgree = true
                    {//truth-
                        //eval- , eval+
                        {ConcordanceState.TRUE_NEGATIVE, null},
                     //truth+
                        //eval-, eval+
                        {null, ConcordanceState.TRUE_POSITIVE}
                    }
                },
             // isFiltered = true
                {// genotypedAgree = false
                    {//truth-
                        //eval- , eval+
                        {ConcordanceState.TRUE_NEGATIVE, ConcordanceState.FILTERED_TRUE_NEGATIVE},
                    //truth+
                        //eval - , eval +
                        {ConcordanceState.FALSE_NEGATIVE, ConcordanceState.FILTERED_FALSE_NEGATIVE}
                    },
                // genotypesAgree = true
                    {//truth-
                        //eval- , eval+
                        {ConcordanceState.TRUE_NEGATIVE, null},
                     //truth+
                        //eval-, eval+
                        {null, ConcordanceState.FILTERED_FALSE_NEGATIVE}
                    }
                }
            };

    // Overriding the superclass method to favor the truth variants' dictionary
    @Override
    public final SAMSequenceDictionary getBestAvailableSequenceDictionary() {
        initializeTruthVariantsIfNecessary();
        return truthVariants.getSequenceDictionary();
    }

    private void initializeTruthVariantsIfNecessary() {
        if (! BucketUtils.fileExists(truthVariantsFile)) {
            throw new IllegalArgumentException("Truth variants file " + truthVariantsFile + " does not exist or is not readable.");
        }

        if (truthVariants == null) {
            truthVariants = new FeatureDataSource<>(new FeatureInput<>(truthVariantsFile, "truth"), CACHE_LOOKAHEAD, VariantContext.class, cloudPrefetchBuffer, cloudIndexPrefetchBuffer);
        }
    }


    protected Predicate<VariantContext> makeTruthVariantFilter() {
        return vc -> !vc.isFiltered();
    }

    // this may often be overridden
    protected Predicate<VariantContext> makeEvalVariantFilter() { return vc -> true; }

    protected abstract boolean shouldVariantsBeMatched(final VariantContext truth, final VariantContext eval);

    private Spliterator<TruthVersusEval> getSpliteratorForDrivingVariants() {
        final Iterator<VariantContext> truthIterator = new FilterIterator<>(truthVariants.iterator(), makeTruthVariantFilter());
        final Iterator<VariantContext> evalIterator = new FilterIterator<>(evalVariants.iterator(), makeEvalVariantFilter());
        return new ConcordanceIterator(truthIterator, evalIterator).spliterator();
    }

    // ********** The basic traversal structure of GATKTool
    @Override
    protected final void onStartup() {
        super.onStartup();

        initializeTruthVariantsIfNecessary();
        evalVariants = new FeatureDataSource<>(new FeatureInput<>(evalVariantsFile, "eval"), CACHE_LOOKAHEAD, VariantContext.class, cloudPrefetchBuffer, cloudIndexPrefetchBuffer);

        if ( hasUserSuppliedIntervals() ) {
            truthVariants.setIntervalsForTraversal(userIntervals);
            evalVariants.setIntervalsForTraversal(userIntervals);
        }
        dict = getBestAvailableSequenceDictionary();
        variantContextComparator = new VariantContextComparator(dict);
    }

    // one can override onTraversalStart() for additional initialization

    // the primary work of the walker.  Must be overridden in implementing classes.
    protected abstract void apply(final TruthVersusEval truthVersusEval, final ReadsContext readsContext, final ReferenceContext refContext);

    /**
     * {@inheritDoc}
     *
     * Implementation of concordance traversal.
     *
     * NOTE: You should only override {@link #traverse()} if you are writing a new walker base class in the
     * engine package that extends this class. It is not meant to be overridden by tools outside of the
     * engine package.
     */
    @Override
    public void traverse() {
        // Process each variant in the input stream.
        StreamSupport.stream(getSpliteratorForDrivingVariants(), false)
                .forEach(truthVersusEval -> {
                    final SimpleInterval variantInterval = new SimpleInterval(truthVersusEval);
                    apply(truthVersusEval, new ReadsContext(reads, variantInterval), new ReferenceContext(reference, variantInterval));
                    progressMeter.update(variantInterval);
                });
    }

    // one can override onTraversalSuccess() for additional post-traversal tasks

    @Override
    protected final void onShutdown() {
        super.onShutdown();
        if( truthVariants != null ) {
            truthVariants.close();
        }
        if( evalVariants != null) {
            evalVariants.close();
        }
    }
    // ********** End of basic traversal methods

    public final VCFHeader getTruthHeader() { return getHeader(truthVariants); }

    public final VCFHeader getEvalHeader() { return getHeader(evalVariants); }

    private static VCFHeader getHeader(final FeatureDataSource<VariantContext> source) {
        final Object header = source.getHeader();
        if ( ! (header instanceof VCFHeader) ) {
            throw new GATKException("Header for " + source.getName() + " is not in VCF header format");
        }
        return (VCFHeader)header;
    }

    protected ConcordanceState getConcordanceState(final Genotype truth, final Genotype eval, final boolean evalWasFiltered) {
        final boolean isTruthPositive = truth != null && isPositive(truth);
        final boolean isEvalPositive = eval != null && isPositive(eval) && !evalWasFiltered;
        final boolean genotypesAfterFilteringAgree;
        if (eval == null || evalWasFiltered) {
            genotypesAfterFilteringAgree = !isTruthPositive;
        } else if (truth == null) {
            genotypesAfterFilteringAgree = !isEvalPositive;
        } else {
            genotypesAfterFilteringAgree = genotypesAgree(truth, eval);
        }

        return evaluateConcordanceState(isEvalPositive, isTruthPositive, genotypesAfterFilteringAgree, evalWasFiltered);
    }

    protected ConcordanceState evaluateConcordanceState(final boolean isPositiveEval, final boolean isPositiveTruth, final boolean genotypesAfterFilteringAgree, final boolean evalIsFiltered) {
        Utils.validate(!(evalIsFiltered && isPositiveEval), "isPositiveEval and evalIsFiltered cannot both be true");
        Utils.validate(!(genotypesAfterFilteringAgree && (isPositiveEval != isPositiveTruth)), "if genotypes agree truth and eval must both be positive or negative");
        if (genotypesAfterFilteringAgree) {
            if (isPositiveEval) {
                return ConcordanceState.TRUE_POSITIVE;
            } else {
                return evalIsFiltered ? ConcordanceState.FILTERED_TRUE_NEGATIVE : ConcordanceState.TRUE_NEGATIVE;
            }
        } else {
            if (isPositiveEval) {
                return ConcordanceState.FALSE_POSITIVE;
            } else {
                return evalIsFiltered ? ConcordanceState.FILTERED_FALSE_NEGATIVE : ConcordanceState.FALSE_NEGATIVE;
            }
        }
    }

    private boolean hasNonRefAlleles(final Collection<Allele> alleles) {
        return alleles.stream().anyMatch(Allele::isNonReference);
    }

    protected boolean isPositive(final Genotype genotype) {
        return hasNonRefAlleles(genotype.getAlleles());
    }

    protected boolean genotypesAgree(final Genotype geno1, final Genotype geno2) {
        // we ignore phasing
        final List<Allele> sortedAlleles1 = new ArrayList<>(geno1.getAlleles());
        Collections.sort(sortedAlleles1);

        final List<Allele> sortedAlleles2 = new ArrayList<>(geno2.getAlleles());
        Collections.sort(sortedAlleles2);

        return sortedAlleles1.equals(sortedAlleles2);
    }

    /**
     * encapsulate the iteration over truth and eval.  hasNext() returns true if either the truth or eval iterator have
     * any variants remaining.  next() produces variants in genomic order, along with their concordance state of
     * false positive for an eval variant without concordant truth variant, false negative for a truth variant without
     * a concordant eval variant, and true positive for concordant eval and truth variants.
     *
     * This implementation correctly consumes a single next() of the wrapped iterators for false positives and
     * false negatives and consumes a next() for both iterators in case of a true positive.
     */
    private class ConcordanceIterator implements Iterator<TruthVersusEval> {
        private final PeekableIterator<VariantContext> truthIterator;
        private final PeekableIterator<VariantContext> evalIterator;
        private final Queue<TruthVersusEval> truthVersusEvalQueue = new ArrayDeque<>();

        protected ConcordanceIterator(final  Iterator<VariantContext> truthIterator, final Iterator<VariantContext> evalIterator) {
            this.truthIterator = new PeekableIterator<>(truthIterator);
            this.evalIterator = new PeekableIterator<>(evalIterator);
        }

        public boolean hasNext() { return truthIterator.hasNext() || evalIterator.hasNext() || !truthVersusEvalQueue.isEmpty(); }

        public TruthVersusEval next() {

            if (!truthVersusEvalQueue.isEmpty()) {
                return truthVersusEvalQueue.poll();
            }

            if (!truthIterator.hasNext()) {
                return TruthVersusEval.evalOnly(evalIterator.next());
            } else if (!evalIterator.hasNext()) {
                return TruthVersusEval.truthOnly(truthIterator.next());
            }

            final int positionCompare = variantContextComparator.compare(truthIterator.peek(), evalIterator.peek());
            if (positionCompare > 0) {
                return TruthVersusEval.evalOnly(evalIterator.next());
            } else if (positionCompare < 0) {
                return TruthVersusEval.truthOnly(truthIterator.next());
            } else {
                //get all eval records at this location
                final List<VariantContext> evalVariants = getAllVariantsAtNextLocus(evalIterator);

                //get all truth records at this location
                final List<VariantContext> truthVariants = getAllVariantsAtNextLocus(truthIterator);

                //if there is only one truth and eval at this location, can just return those, otherwise need to match using matching function
                if (evalVariants.size() == 1 && truthVariants.size() == 1) {
                    final VariantContext truthVariant = truthVariants.get(0);
                    final VariantContext evalVariant = evalVariants.get(0);
                    if (shouldVariantsBeMatched(truthVariant, evalVariant)) {
                        return new TruthVersusEval(truthVariant, evalVariant);
                    } else {
                        truthVersusEvalQueue.add(TruthVersusEval.evalOnly(evalVariant));
                        return TruthVersusEval.truthOnly(truthVariant);
                    }
                } else {
                    truthVersusEvalQueue.addAll(matchVariants(truthVariants, evalVariants));
                    return truthVersusEvalQueue.poll();
                }
            }
        }

        private List<VariantContext> getAllVariantsAtNextLocus(final PeekableIterator<VariantContext> iterator) {
            Utils.validateArg(iterator.hasNext(), "AbstractConcordanceWalker asked for next variant in variant group when variant group was empty");
            final List<VariantContext> variants = new ArrayList<>();
            VariantContext currentVariant = iterator.next();
            variants.add(currentVariant);
            while(iterator.hasNext() && variantContextComparator.compare(currentVariant, iterator.peek()) == 0) {
                currentVariant = iterator.next();
                variants.add(currentVariant);
            }

            return variants;
        }

        protected Queue<TruthVersusEval> matchVariants(final List<VariantContext> truthVariants, final List<VariantContext> evalVariants) {
                //match based on shouldVariantsBeMatched
                //we expect only a small number of variants in each list, so a loop is fine
            final Queue<TruthVersusEval> queue = new ArrayDeque<>();
            Set<VariantContext> unmatchedTruth  = new HashSet<>(truthVariants);
            for (final VariantContext eval : evalVariants) {
                boolean matched = false;
                for (final VariantContext truth : truthVariants) {
                    if (shouldVariantsBeMatched(truth, eval)) {
                        if (!matched) {
                            queue.add(new TruthVersusEval(truth, eval));
                            matched = true;
                        }
                        unmatchedTruth.remove(truth);
                    }
                }
                if (!matched) {
                    queue.add(TruthVersusEval.evalOnly(eval));
                }
            }

            //deal with unmatched truth
            for (final VariantContext truth : unmatchedTruth) {
                queue.add(TruthVersusEval.truthOnly(truth));
            }

            return queue;
        }

        private Spliterator<TruthVersusEval> spliterator() {
            return Spliterators.spliteratorUnknownSize(this, 0);
        }
    }

    /**
     * store a truth vc in case of a false negative, an eval vc in case of a false positive, or a concordance pair of
     * truth and eval in case of a true positive.
     */
    protected static class TruthVersusEval implements Locatable {
        protected final VariantContext truth;
        protected final VariantContext eval;

        public TruthVersusEval(final VariantContext truth, final VariantContext eval) {
            this.truth = truth;
            this.eval = eval;
        }

        public static TruthVersusEval evalOnly(final VariantContext eval) {
            return new TruthVersusEval(null, eval);
        }

        public static TruthVersusEval truthOnly(final VariantContext truth) {
            return new TruthVersusEval(truth, null);
        }

        public VariantContext getTruth() {
            return truth;
        }

        public VariantContext getEval() {
            return eval;
        }

        public boolean hasTruth() {return truth != null; }
        public boolean hasEval() { return eval != null; }

        public VariantContext getTruthIfPresentElseEval() { return truth != null ? truth : eval; }

        public String getContig() { return getTruthIfPresentElseEval().getContig(); }
        public int getStart() { return getTruthIfPresentElseEval().getStart(); }
        public int getEnd() { return getTruthIfPresentElseEval().getEnd(); }
    }

}
