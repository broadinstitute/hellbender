package org.broadinstitute.hellbender.tools.dragstr;

import org.apache.commons.lang3.Range;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParams;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParamsBuilder;

import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * Implements the DRAGstr parameter estimation routine.
 */
final class DragstrParametersEstimator {

    private final DragstrHyperParameters hyperParameters;

    /**
     * Holds the possible values for the GP parameter in Phred scale (the scale used by the end user)
     */
    private final double[] phredGpValues;

    /**
     * Holds the possible values for the API parameter in Phred scale.
     */
    private final double[] phredApiValues;

    /**
     * Holds the same values as {@link #phredGpValues} but in log10 scale.
     */
    private final double[] log10GpValues;

    /**
     * Holds the same values as {@link #phredApiValues} but in log10 scale.
     */
    private final double[] log10ApiValues;

    /**
     * Fix ratio between the het. and hom-var genotype priors in linear scale [0 .. 1].
     */
    private final double hetOverHomVar;

    /**
     * Same as {@link #hetOverHomVar} but in log10 scale.
     */
    private final double log10HetOverHomVar;

    /**
     * For a given period, holds the index on {@link #phredGpValues} (thus {@link #log10GpValues} of the minimum possible best GP
     * for that period. The array in fact indexed as period - 1.
     */
    private final int[] minGpIndexByPeriod;

    /**
     * Holds the log10 scaled probability of an error for a given GP (by its index in {@link #phredGpValues}), period
     * and repeat length. This only needs to be calculated once.
     */
    private final double[][][] log10PErrorByLength;

    /**
     * Holds the log10 scaled probability of no error for a given GP (by its index in {@link #phredGpValues}), period
     * and repeat length. This only needs to be calculated once.
     */
    private final double[][][] log10PCorrectByLength;

    DragstrParametersEstimator(final DragstrHyperParameters argumentCollection) {
        hyperParameters = argumentCollection;
        phredGpValues = argumentCollection.phredGpValues.toDoubleArray();
        phredApiValues = argumentCollection.phredApiValues.toDoubleArray();
        log10GpValues = MathUtils.applyToArray(phredGpValues, d -> -.1 * d);
        log10ApiValues = MathUtils.applyToArray(phredApiValues, d -> -0.1 * d);
        hetOverHomVar = argumentCollection.hetToHomRatio;
        log10HetOverHomVar = Math.log10(hetOverHomVar);
        log10PErrorByLength = new double[phredGpValues.length][hyperParameters.maxPeriod][hyperParameters.maxRepeatLength];
        log10PCorrectByLength = new double[phredGpValues.length][hyperParameters.maxPeriod][hyperParameters.maxRepeatLength];
        minGpIndexByPeriod = new int[hyperParameters.maxPeriod];
        for (int i = 0; i < phredGpValues.length; i++) {
            for (int k = 0; k < hyperParameters.maxPeriod; k++) {
                final int period = k + 1;
                final double log10PCorrectPerPosition = MathUtils.log10OneMinusPow10(- MathUtils.LOG10_ONE_HALF + log10GpValues[i]);
                for (int j = 0; j < hyperParameters.maxRepeatLength; j++) {
                    final int repeats = j + 1;
                    final int lengthInBases = repeats * period;
                    log10PCorrectByLength[i][k][j] = lengthInBases * log10PCorrectPerPosition;
                    log10PErrorByLength[i][k][j] = MathUtils.log10OneMinusPow10(log10PCorrectByLength[i][k][j]);
                }
            }
        }
        for (int i = 0; i < minGpIndexByPeriod.length; i++) {
            final int period = i + 1;
            // Some magic numbers here... this is a trascription of Illumina's matlab script code.
            // Instead of hyperParams.maxRepeatLength there was a 20.0 and since is multiplied to period
            // it makes sense that it reference to this max that is 20 by default.
            final double gpMin = Math.ceil(-10 * Math.log10((1 - Math.pow(0.5, (1.0 / (hyperParameters.maxRepeatLength * period)) / 2.0))));
            final int index = Arrays.binarySearch(phredGpValues, gpMin);
            // since we are looking for a double, we have to be a bit tolerant in terms of differences,
            // so if no exact was found we look the position before the insertion position in case that value
            // is close enough (less than 0.001 way).
            minGpIndexByPeriod[i] = index >= 0 ? index :
                    (index < -1 && Math.abs(gpMin - phredApiValues[-index - 2]) < 0.001) ? -index - 2 : -index - 1;
        }
    }

    /**
     * Estimates the DRAGstr model parameters given the final set of cases to be used for the task.
     * @param cases the input targets.
     * @return never {@code null}.
     */
    public DragstrParams estimate(final StratifiedDragstrLocusCases cases) {
        final DragstrParamsBuilder builder = new DragstrParamsBuilder(hyperParameters.maxPeriod, hyperParameters.maxRepeatLength);
        IntStream.range(1, hyperParameters.maxPeriod + 1).parallel().forEach(period -> {
            estimatePeriod(period, builder, cases);
        });
        return builder.make(hyperParameters.phredGopValues);
    }

    /**
     * Perform the estimation for a single period.
     * <p>
     *     Parameters are estimated per period independently of other periods.
     * </p>
     * @param period the target period.
     * @param destination the builder where to store the estimated values for that period.
     * @param cases the input data to make such estimates.
     */
    private void estimatePeriod(final int period, final DragstrParamsBuilder destination, final StratifiedDragstrLocusCases cases) {
        // We calculate the "flanks" or repeat-length counts that have too little data
        // that seat in either side of the repeat-length range.
        // These will be grouped before provide estimates for those repeat-lengths:
        int accum = 0;
        int leftFlank;
        for (leftFlank = 0; leftFlank < hyperParameters.maxRepeatLength;) {
            if ((accum += cases.get(period, ++leftFlank).size()) >= hyperParameters.minLociCount) {
                break;
            }
        }
        int rightFlank; accum = 0;
        for (rightFlank = hyperParameters.maxRepeatLength; rightFlank > 1;) {
            if ((accum += cases.get(period, --rightFlank).size()) >= hyperParameters.minLociCount) {
                break;
            }
        }

        final ArrayDeque<Range<Integer>> pending = new ArrayDeque<>(hyperParameters.maxRepeatLength);

        if (rightFlank >= leftFlank) {
            // We fill 'pending' with the repeat-length groups that will be analyze:
            // [1 .. leftFlank], leftFlank + 1, leftFlank + 2, ... , [rightFlank .. maxRepeats+]
            pending.add(Range.of(1, leftFlank));
            for (leftFlank++; leftFlank <= rightFlank; leftFlank++) {
                pending.add(Range.is(leftFlank));
            }
            pending.add(Range.of(++rightFlank, hyperParameters.maxRepeatLength));
        } else { // no enough data we simply put all repeat lengths is one single group:
            pending.add(Range.of(1, hyperParameters.maxRepeatLength));
        }

        Range<Integer> last = null;

        // Done will contain the ranges already processed.
        final ArrayDeque<Range<Integer>> done = new ArrayDeque<>(hyperParameters.maxRepeatLength);

        do {
            final Range<Integer> next = pending.pop();
            estimatePeriodRepeatInterval(period, next, destination, cases);
            final double gp1 = destination.gp(period, next.getMinimum());
            final double api1 = destination.api(period, next.getMinimum());
            // if GP and API are "decreasing" with respect those from smaller repeat length
            // then we accepted them:
            if (last == null || (destination.gp(period, last.getMaximum()) >= gp1 &&
                                   destination.api(period, last.getMaximum()) + hyperParameters.apiMonothresh >= api1)) {
                done.addLast(last = next);
            // if not, the we group back this repeat-length group/range with last one and re-estimate (next-loop).
            } else {
                pending.push(Range.of(last.getMinimum(), next.getMaximum()));
                done.removeLast();
                last = !done.isEmpty() ? done.getLast() : null;
            }
        } while (!pending.isEmpty());
    }

    // Given a observed het/hom ratio the total number of reads, the variable ones and the length of the STR in bases.
    // calculate the optimal gp and api that maximizes the likelihood.
    private void estimatePeriodRepeatInterval(final int period, final Range<Integer> repeatRange,
                                              final DragstrParamsBuilder builder,
                                              final StratifiedDragstrLocusCases cases) {
        int maxApiIdx = -1;
        int maxGpIdx = -1;
        double maxLog10Prob = Double.NEGATIVE_INFINITY;
        final int minRepeat = repeatRange.getMinimum();
        final int maxRepeat = repeatRange.getMaximum();
        final int periodIdx = period - 1;
        final double maxLog10PHet = log10HetOverHomVar - Math.log10(1 + hyperParameters.hetToHomRatio);
        for (int i = 0; i < log10ApiValues.length; i++) {
            final double log10PHet = Math.min(log10ApiValues[i], maxLog10PHet);
            final double log10PHomVar = log10PHet - log10HetOverHomVar;
            // Pr homref = 1 - (Pr. het + Pr hom) in log10 scale:
            final double log10PHomRef = MathUtils.log10OneMinusPow10(MathUtils.log10SumLog10(log10PHet, log10PHomVar));
     gp_for: for (int j = minGpIndexByPeriod[periodIdx]; j < phredGpValues.length; j++) {
                double log10ProbAccumulator = 0;
                for (int r = minRepeat, repeatIdx = r - 1; r <= maxRepeat; r++, repeatIdx++) {
                    final DragstrLocusCases repeatCases = cases.get(period, r);
                    final double log10PError = log10PErrorByLength[j][periodIdx][repeatIdx];
                    final double log10PCorrect = log10PCorrectByLength[j][periodIdx][repeatIdx];
                    for (final DragstrLocusCase caze : repeatCases) {
                        if ((log10ProbAccumulator += log10ProbFunc(caze.getDepth(), caze.getIndels(), log10PError, log10PCorrect,
                                  log10PHomRef, log10PHet, log10PHomVar)) < maxLog10Prob) {
                            continue gp_for; // abort if is already smaller than the maximum so far.
                        };
                    }
                }
                if (log10ProbAccumulator > maxLog10Prob) {
                    maxApiIdx = i;
                    maxGpIdx = j;
                    maxLog10Prob = log10ProbAccumulator;
                }
            }
        }
        builder.set(period, repeatRange, phredGpValues[maxGpIdx], 10.0 / period, phredApiValues[maxApiIdx]);
    }

    private double log10ProbFunc(final int n, final int k, final double log10PError, final double log10PCorrect,
                                 final double log10PHomRef, final double log10PHet, final double log10PHomVar) {
        return MathUtils.log10SumLog10(
                log10PHomRef + k * log10PError + (n - k) * log10PCorrect,
                log10PHet + n * MathUtils.LOG10_ONE_HALF,
                   n == k ? log10PHomVar : Double.NEGATIVE_INFINITY); // this last hom_ref lk part is not accurate as is possible that an error that reversts to the reference does occur.
    }
}
