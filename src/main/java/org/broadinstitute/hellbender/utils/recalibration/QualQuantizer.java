package org.broadinstitute.hellbender.utils.recalibration;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/**
 * A general algorithm for quantizing quality score distributions to use a specific number of levels
 *
 * Takes a histogram of quality scores and a desired number of levels and produces a
 * map from original quality scores -> quantized quality scores.
 *
 * Note that this data structure is fairly heavy-weight, holding lots of debugging and
 * calculation information.  If you want to use it efficiently at scale with lots of
 * read groups the right way to do this:
 *
 * Map<ReadGroup, List<Byte>> map
 * for each read group rg:
 *   hist = getQualHist(rg)
 *   QualQuantizer qq = new QualQuantizer(hist, nLevels, minInterestingQual)
 *   map.set(rg, qq.getOriginalToQuantizedMap())
 *
 * This map would then be used to look up the appropriate original -> quantized
 * quals for each read as it comes in.
 */
public final class QualQuantizer {
    private static final Set<QualInterval> MY_EMPTY_SET = Collections.emptySet();

    private static final Logger logger = LogManager.getLogger(QualQuantizer.class);

    /**
     * Inputs to the QualQuantizer
     */
    final int nLevels, minInterestingQual;
    final List<Long> nObservationsPerQual;

    /**
     * Map from original qual (e.g., Q30) to new quantized qual (e.g., Q28).
     *
     * Has the same range as nObservationsPerQual
     */
    final List<Byte> originalToQuantizedMap;

    /** Sorted set of qual intervals.
     *
     * After quantize() this data structure contains only the top-level qual intervals
     */
    final TreeSet<QualInterval> quantizedIntervals;

    /**
     * Protected creator for testng use only
     */
    protected QualQuantizer(final int minInterestingQual) {
        this.nObservationsPerQual = Collections.emptyList();
        this.nLevels = 0;
        this.minInterestingQual = minInterestingQual;
        this.quantizedIntervals = null;
        this.originalToQuantizedMap = null;
    }

    /**
     * Creates a QualQuantizer for the histogram that has nLevels
     *
     * Note this is the only interface to the system.  After creating this object
     * the map can be obtained via getOriginalToQuantizedMap()
     *
     * @param nObservationsPerQual A histogram of counts of bases with quality scores.  Note that
     *  this histogram must start at 0 (i.e., get(0) => count of Q0 bases) and must include counts all the
     *  way up to the largest quality score possible in the reads.  OK if the histogram includes many 0
     *  count bins, as these are quantized for free.
     * @param nLevels the desired number of distinct quality scores to represent the full original range.  Must
     *  be at least 1.
     * @param minInterestingQual All quality scores <= this value are considered uninteresting and are freely
     *  merged together.  For example, if this value is 10, then Q0-Q10 are all considered free to merge, and
     *  quantized into a single value. For ILMN data with lots of Q2 bases this results in a Q2 bin containing
     *  all data with Q0-Q10.
     */
    public QualQuantizer(final List<Long> nObservationsPerQual, final int nLevels, final int minInterestingQual) {
        this.nObservationsPerQual = nObservationsPerQual;
        this.nLevels = nLevels;
        this.minInterestingQual = minInterestingQual;

        // some sanity checking
        if ( Collections.min(nObservationsPerQual) < 0 ) throw new GATKException("Quality score histogram has negative values at: " + Utils.join(", ", nObservationsPerQual));
        if ( nLevels < 0 ) throw new GATKException("nLevels must be >= 0");
        if ( minInterestingQual < 0 ) throw new GATKException("minInterestingQual must be >= 0");

        // actually run the quantizer
        this.quantizedIntervals = quantize();

        // store the map
        this.originalToQuantizedMap = intervalsToMap(quantizedIntervals);
    }

    /**
     * Represents a contiguous interval of quality scores.
     *
     * qStart and qEnd are inclusive, so qStart = qEnd = 2 is the quality score bin of 2
     */
    protected final class QualInterval implements Comparable<QualInterval> {
        final int qStart, qEnd, fixedQual, level;
        final long nObservations, nErrors;
        final Set<QualInterval> subIntervals;

        /** for debugging / visualization.  When was this interval created? */
        int mergeOrder;

        protected QualInterval(final int qStart, final int qEnd, final long nObservations, final long nErrors, final int level) {
            this(qStart, qEnd, nObservations, nErrors, level, -1, MY_EMPTY_SET);
        }

        protected QualInterval(final int qStart, final int qEnd, final long nObservations, final long nErrors, final int level, final Set<QualInterval> subIntervals) {
            this(qStart, qEnd, nObservations, nErrors, level, -1, subIntervals);
        }

        protected QualInterval(final int qStart, final int qEnd, final long nObservations, final long nErrors, final int level, final int fixedQual) {
            this(qStart, qEnd, nObservations, nErrors, level, fixedQual, MY_EMPTY_SET);
        }

        public QualInterval(final int qStart, final int qEnd, final long nObservations, final long nErrors, final int level, final int fixedQual, final Set<QualInterval> subIntervals) {
            this.qStart = qStart;
            this.qEnd = qEnd;
            this.nObservations = nObservations;
            this.nErrors = nErrors;
            this.fixedQual = fixedQual;
            this.level = level;
            this.mergeOrder = 0;
            this.subIntervals = Collections.unmodifiableSet(subIntervals);
        }

        /**
         * @return Human readable name of this interval: e.g., 10-12
         */
        public String getName() {
            return qStart + "-" + qEnd;
        }

        @Override
        public String toString() {
            return "QQ:" + getName();
        }

        /**
         * @return the error rate (in real space) of this interval, or 0 if there are no observations
         */
        public double getErrorRate() {
            if ( hasFixedQual() )
                return QualityUtils.qualToErrorProb((byte) fixedQual);
            else if ( nObservations == 0 )
                return 0.0;
            else
                return (nErrors+1) / (1.0 * (nObservations+1));
        }

        /**
         * @return the QUAL of the error rate of this interval, or the fixed qual if this interval was created with a fixed qual.
         */
        public byte getQual() {
            if ( ! hasFixedQual() )
                return QualityUtils.errorProbToQual(getErrorRate());
            else
                return (byte)fixedQual;
        }

        /**
         * @return true if this bin is using a fixed qual
         */
        public boolean hasFixedQual() {
            return fixedQual != -1;
        }

        @Override
        public int compareTo(final QualInterval qualInterval) {
            return Integer.valueOf(this.qStart).compareTo(qualInterval.qStart);
        }

        /**
         * Create a interval representing the merge of this interval and toMerge
         *
         * Errors and observations are combined
         * Subintervals updated in order of left to right (determined by qStart)
         * Level is 1 + highest level of this and toMerge
         * Order must be updated elsewhere
         *
         * @param toMerge
         * @return newly created merged QualInterval
         */
        public QualInterval merge(final QualInterval toMerge) {
            final QualInterval left = this.compareTo(toMerge) < 0 ? this : toMerge;
            final QualInterval right = this.compareTo(toMerge) < 0 ? toMerge : this;

            if ( left.qEnd + 1 != right.qStart )
                throw new GATKException("Attempting to merge non-contiguous intervals: left = " + left + " right = " + right);

            final long nCombinedObs = left.nObservations + right.nObservations;
            final long nCombinedErr = left.nErrors + right.nErrors;

            final int level = Math.max(left.level, right.level) + 1;
            final Set<QualInterval> subIntervals = new LinkedHashSet<>(Arrays.asList(left, right));
            QualInterval merged = new QualInterval(left.qStart, right.qEnd, nCombinedObs, nCombinedErr, level, subIntervals);

            return merged;
        }

        public double getPenalty() {
            return calcPenalty(getErrorRate());
        }


        /**
         * Calculate the penalty of this interval, given the overall error rate for the interval
         *
         * If the globalErrorRate is e, this value is:
         *
         * sum_i |log10(e_i) - log10(e)| * nObservations_i
         *
         * each the index i applies to all leaves of the tree accessible from this interval
         * (found recursively from subIntervals as necessary)
         *
         * @param globalErrorRate overall error rate in real space against which we calculate the penalty
         * @return the cost of approximating the bins in this interval with the globalErrorRate
         */
        private double calcPenalty(final double globalErrorRate) {
            if ( globalErrorRate == 0.0 ) // there were no observations, so there's no penalty
                return 0.0;

            if ( subIntervals.isEmpty() ) {
                // this is leave node
                if ( this.qEnd <= minInterestingQual )
                    // It's free to merge up quality scores below the smallest interesting one
                    return 0;
                else {
                    return (Math.abs(Math.log10(getErrorRate()) - Math.log10(globalErrorRate))) * nObservations;
                }
            } else {
                double sum = 0;
                for ( final QualInterval interval : subIntervals )
                    sum += interval.calcPenalty(globalErrorRate);
                return sum;
            }
        }
    }

    /**
     * Main method for computing the quantization intervals.
     *
     * Invoked in the constructor after all input variables are initialized.  Walks
     * over the inputs and builds the min. penalty forest of intervals with exactly nLevel
     * root nodes.  Finds this min. penalty forest via greedy search, so is not guarenteed
     * to find the optimal combination.
     *
     * TODO: develop a smarter algorithm
     *
     * @return the forest of intervals with size == nLevels
     */
    private TreeSet<QualInterval> quantize() {
        // create intervals for each qual individually
        final TreeSet<QualInterval> intervals = new TreeSet<>();
        for ( int qStart = 0; qStart < getNQualsInHistogram(); qStart++ ) {
            final long nObs = nObservationsPerQual.get(qStart);
            final double errorRate = QualityUtils.qualToErrorProb((byte)qStart);
            final double nErrors = nObs * errorRate;
            final QualInterval qi = new QualInterval(qStart, qStart, nObs, (int) Math.floor(nErrors), 0, (byte)qStart);
            intervals.add(qi);
        }

        // greedy algorithm:
        // while ( n intervals >= nLevels ):
        //   find intervals to merge with least penalty
        //   merge it
        while ( intervals.size() > nLevels ) {
            mergeLowestPenaltyIntervals(intervals);
        }

        return intervals;
    }

    /**
     * Helper function that finds and merges together the lowest penalty pair of intervals
     * @param intervals
     */
    private void mergeLowestPenaltyIntervals(final TreeSet<QualInterval> intervals) {
        // setup the iterators
        final Iterator<QualInterval> it1 = intervals.iterator();
        final Iterator<QualInterval> it1p = intervals.iterator();
        it1p.next(); // skip one

        // walk over the pairs of left and right, keeping track of the pair with the lowest merge penalty
        QualInterval minMerge = null;
        if ( logger.isDebugEnabled() ) logger.debug("mergeLowestPenaltyIntervals: " + intervals.size());
        int lastMergeOrder = 0;
        while ( it1p.hasNext() ) {
            final QualInterval left = it1.next();
            final QualInterval right = it1p.next();
            final QualInterval merged = left.merge(right);
            lastMergeOrder = Math.max(Math.max(lastMergeOrder, left.mergeOrder), right.mergeOrder);
            if ( minMerge == null || (merged.getPenalty() < minMerge.getPenalty() ) ) {
                if ( logger.isDebugEnabled() ) logger.debug("  Updating merge " + minMerge);
                minMerge = merged; // merge two bins that when merged incur the lowest "penalty"
            }
        }

        // now actually go ahead and merge the minMerge pair
        if ( logger.isDebugEnabled() ) logger.debug("  => final min merge " + minMerge);
        intervals.removeAll(minMerge.subIntervals);
        intervals.add(minMerge);
        minMerge.mergeOrder = lastMergeOrder + 1;
        if ( logger.isDebugEnabled() ) logger.debug("updated intervals: " + intervals);
    }

    /**
     * Given a final forest of intervals constructs a list mapping
     * list.get(i) => quantized qual to use for original quality score i
     *
     * This function should be called only once to initialize the corresponding
     * cached value in this object, as the calculation is a bit costly.
     *
     * @param intervals
     * @return
     */
    private List<Byte> intervalsToMap(final TreeSet<QualInterval> intervals) {
        final List<Byte> map = new ArrayList<>(getNQualsInHistogram());
        map.addAll(Collections.nCopies(getNQualsInHistogram(), Byte.MIN_VALUE));
        for ( final QualInterval interval : intervals ) {
            for ( int q = interval.qStart; q <= interval.qEnd; q++ ) {
                map.set(q, interval.getQual());
            }
        }

        if ( Collections.min(map) == Byte.MIN_VALUE )
            throw new GATKException("quantized quality score map contains an un-initialized value");

        return map;
    }

    private int getNQualsInHistogram() {
        return nObservationsPerQual.size();
    }

    public List<Byte> getOriginalToQuantizedMap() {
        return originalToQuantizedMap;
    }
}
