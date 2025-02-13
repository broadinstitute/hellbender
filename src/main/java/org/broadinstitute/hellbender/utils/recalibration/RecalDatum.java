package org.broadinstitute.hellbender.utils.recalibration;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMUtils;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;

import java.io.Serializable;

/**
 * The container for the 4-tuple
 *
 *   ( reported quality, empirical quality, num observations, num mismatches/errors )
 *
 * for a given set of covariates.
 *
 */
public final class RecalDatum implements Serializable {
    public static final byte MAX_RECALIBRATED_Q_SCORE = SAMUtils.MAX_PHRED_SCORE;
    private static final int UNINITIALIZED_EMPIRICAL_QUALITY = -1;
    private static final long serialVersionUID = 1L;
    private static final double MULTIPLIER = 100000.0;  //See discussion in numMismatches about what the multiplier is.

    /**
     * Estimated reported quality score based on combined data's individual q-reporteds and number of observations.
     * The estimating occurs when collapsing counts across different reported qualities.
     * TODO: investigate whether this estimated reported quality score can be eliminated.
     *
     */
    private double reportedQuality;

    /**
     * The empirical quality for datums that have been collapsed together (by read group and reported quality, for example).
     *
     * This variable was historically a double, but {@link #bayesianEstimateOfEmpiricalQuality} has always returned an integer qual score.
     * Thus the type has been changed to integer in February 2025 to highlight this implementation detail. It does not change the output.
     */
    private int empiricalQuality;

    /**
     * Number of bases seen in total
     */
    private long numObservations;

    /**
     * Number of bases seen that didn't match the reference
     * (actually sum of the error weights - so not necessarily a whole number)
     * Stored with an internal multiplier to keep it closer to the floating-point sweet spot and avoid numerical error
     * (see https://github.com/broadinstitute/gatk/wiki/Numerical-errors ).
     * However, the value of the multiplier influences the results.
     * For example, you get different results for 1000.0 and 10000.0
     * See MathUtilsUnitTest.testAddDoubles for a demonstration.
     * The value of the MULTIPLIER that we found to give consistent results insensitive to sorting is 10000.0;
     */
    private double numMismatches;

    /**
     * used when calculating empirical qualities to avoid division by zero
     */
    private static final int SMOOTHING_CONSTANT = 1;

    //---------------------------------------------------------------------------------------------------------------
    //
    // constructors
    //
    //---------------------------------------------------------------------------------------------------------------

    /**
     * Create a new RecalDatum with given observation and mismatch counts, and an reported quality
     *
     * @param numObservations
     * @param numMismatches
     * @param reportedQuality
     */
    public RecalDatum(final long numObservations, final double numMismatches, final byte reportedQuality) {
        if ( numObservations < 0 ) throw new IllegalArgumentException("numObservations < 0");
        if ( numMismatches < 0.0 ) throw new IllegalArgumentException("numMismatches < 0");
        if ( reportedQuality < 0 ) throw new IllegalArgumentException("reportedQuality < 0");

        this.numObservations = numObservations;
        this.numMismatches = (numMismatches*MULTIPLIER);
        this.reportedQuality = reportedQuality;
        this.empiricalQuality = UNINITIALIZED_EMPIRICAL_QUALITY;
    }

    /**
     * Copy copy into this recal datum, overwriting all of this objects data
     * @param copy  RecalDatum to copy
     */
    public RecalDatum(final RecalDatum copy) {
        this.numObservations = copy.numObservations;
        this.numMismatches = copy.numMismatches;
        this.reportedQuality = copy.reportedQuality;
        this.empiricalQuality = copy.empiricalQuality;
    }

    /**
     * Add in all of the data from other into this object, updating the reported quality from the expected
     * error rate implied by the two reported qualities.
     *
     * For example (the only example?), this method is called when collapsing the counts across reported quality scores within
     * the same read group.
     *
     * @param other  RecalDatum to combine
     */
    public void combine(final RecalDatum other) {
        // this is the *expected* (or theoretical) number of errors given the reported qualities and the number of observations.
        final double expectedNumErrors = this.calcExpectedErrors() + other.calcExpectedErrors();

        // increment the counts
        increment(other.getNumObservations(), other.getNumMismatches());
        
        // we use the theoretical count above to compute the "estimated" reported quality
        // after combining two datums with different reported qualities.
        reportedQuality = -10 * Math.log10(expectedNumErrors / getNumObservations());
        empiricalQuality = UNINITIALIZED_EMPIRICAL_QUALITY;
    }

    public void setReportedQuality(final double reportedQuality) {
        if ( reportedQuality < 0 ) throw new IllegalArgumentException("estimatedQReported < 0");
        if ( Double.isInfinite(reportedQuality) ) throw new IllegalArgumentException("estimatedQReported is infinite");
        if ( Double.isNaN(reportedQuality) ) throw new IllegalArgumentException("estimatedQReported is NaN");

        this.reportedQuality = reportedQuality;
        empiricalQuality = UNINITIALIZED_EMPIRICAL_QUALITY;
    }

    public double getReportedQuality() {
        return reportedQuality;
    }
    public byte getReportedQualityAsByte() {
        return (byte)(int)(Math.round(getReportedQuality()));
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // Empirical quality score -- derived from the num mismatches and observations
    //
    //---------------------------------------------------------------------------------------------------------------

    /**
     * Returns the error rate (in real space) of this interval, or 0 if there are no observations
     * @return the empirical error rate ~= N errors / N obs
     */
    public double getEmpiricalErrorRate() {
        if ( numObservations == 0 )
            return 0.0;
        else {
            // cache the value so we don't call log over and over again
            final double doubleMismatches = (numMismatches/MULTIPLIER) + SMOOTHING_CONSTANT;
            // smoothing is one error and one non-error observation, for example
            final double doubleObservations = numObservations + SMOOTHING_CONSTANT + SMOOTHING_CONSTANT;
            return doubleMismatches / doubleObservations;
        }
    }

    public void setEmpiricalQuality(final int empiricalQuality) {
        if ( empiricalQuality < 0 ) throw new IllegalArgumentException("empiricalQuality < 0");
        if ( Double.isInfinite(empiricalQuality) ) throw new IllegalArgumentException("empiricalQuality is infinite");
        if ( Double.isNaN(empiricalQuality) ) throw new IllegalArgumentException("empiricalQuality is NaN");

        this.empiricalQuality = empiricalQuality;
    }

    /**
     * Computes the empirical quality of the datum, using the reported quality as the prior.
     * @see #getEmpiricalQuality(double) below.
     */
    public double getEmpiricalQuality() {
        return getEmpiricalQuality(getReportedQuality());
    }

    /**
     * Computes the empirical base quality (roughly (num errors)/(num observations)) from the counts stored in this datum.
     */
    public double getEmpiricalQuality(final double priorQualityScore) {
        if (empiricalQuality == UNINITIALIZED_EMPIRICAL_QUALITY) {
            calcEmpiricalQuality(priorQualityScore);
        }
        return empiricalQuality;
    }

    public byte getEmpiricalQualityAsByte() {
        return (byte)(Math.round(getEmpiricalQuality()));
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // toString methods
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public String toString() {
        return String.format("%d,%.2f,%.2f", getNumObservations(), getNumMismatches(), getEmpiricalQuality());
    }

    public String stringForCSV() {
        return String.format("%s,%.2f,%.2f", toString(), getReportedQuality(), getEmpiricalQuality() - getReportedQuality());
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // increment methods
    //
    //---------------------------------------------------------------------------------------------------------------

    public long getNumObservations() {
        return numObservations;
    }

    public void setNumObservations(final long numObservations) {
        if ( numObservations < 0 ) throw new IllegalArgumentException("numObservations < 0");
        this.numObservations = numObservations;
        empiricalQuality = UNINITIALIZED_EMPIRICAL_QUALITY;
    }

    public double getNumMismatches() {
        return numMismatches/MULTIPLIER;
    }

    public void setNumMismatches(final double numMismatches) {
        if ( numMismatches < 0 ) throw new IllegalArgumentException("numMismatches < 0");
        this.numMismatches = (numMismatches*MULTIPLIER);
        empiricalQuality = UNINITIALIZED_EMPIRICAL_QUALITY;
    }

    public void incrementNumObservations(final long by) {
        numObservations += by;
        empiricalQuality = UNINITIALIZED_EMPIRICAL_QUALITY;
    }

    public void incrementNumMismatches(final double by) {
        numMismatches += (by*MULTIPLIER);
        empiricalQuality = UNINITIALIZED_EMPIRICAL_QUALITY;
    }

    public void increment(final long incObservations, final double incMismatches) {
        numObservations += incObservations;
        numMismatches += (incMismatches*MULTIPLIER); // the multiplier used to avoid underflow, or something like that.
        empiricalQuality = UNINITIALIZED_EMPIRICAL_QUALITY;
    }

    public void increment(final boolean isError) {
        increment(1, isError ? 1.0 : 0.0);
    }

    // -------------------------------------------------------------------------------------
    //
    // Private implementation helper functions
    //
    // -------------------------------------------------------------------------------------

    /**
     * calculate the expected number of errors given the estimated Q reported and the number of observations
     * in this datum.
     *
     * @return a positive (potentially fractional) estimate of the number of errors
     */
    private double calcExpectedErrors() {
        return getNumObservations() * QualityUtils.qualToErrorProb(reportedQuality);
    }

    /**
     * Calculate and cache the empirical quality score from mismatches and observations (expensive operation)
     */
    private void calcEmpiricalQuality(final double priorQualityScore) {
        // smoothing is one error and one non-error observation
        final long mismatches = (long)(getNumMismatches() + 0.5) + SMOOTHING_CONSTANT; // TODO: why add 0.5?
        final long observations = getNumObservations() + SMOOTHING_CONSTANT + SMOOTHING_CONSTANT;

        final int empiricalQual = bayesianEstimateOfEmpiricalQuality(observations, mismatches, priorQualityScore);

        empiricalQuality = Math.min(empiricalQual, MAX_RECALIBRATED_Q_SCORE);
    }

    /**
     * Compute the maximum a posteriori (MAP) estimate of the probability of sequencing error under the following model.
     *
     * Let
     *   X = number of sequencing errors,
     *   n = number of observations,
     *   theta = probability of sequencing error as a quality score,
     *   theta_rep = probability of sequencing error reported by the sequencing machine as a quality score.
     *
     * The prior and the likelihood are:
     *
     * P(theta|theta_rep) ~ Gaussian(theta - theta_rep| 0, 0.5) (Note this is done in log space)
     * P(X|n, theta) ~ Binom(X|n,theta)
     *
     * Note the prior is equivalent to
     *
     * P(theta|theta_rep) ~ Gaussian(theta | theta_rep, 0.5)
     *
     * TODO: use beta prior to do away with the search.
     *
     * @param nObservations n in the model above.
     * @param nErrors the observed number of sequencing errors.
     * @param priorMeanQualityScore the prior quality score, often the reported quality score.
     *
     * @return phredScale quality score that maximizes the posterior probability.
     */
    public static int bayesianEstimateOfEmpiricalQuality(final long nObservations, final long nErrors, final double priorMeanQualityScore) {
        final int numQualityScoreBins = (QualityUtils.MAX_REASONABLE_Q_SCORE + 1);

        final double[] logPosteriors = new IndexRange(0, numQualityScoreBins).mapToDouble(q ->
            getLogPrior(q, priorMeanQualityScore) + getLogBinomialLikelihood(q, nObservations, nErrors));
        return MathUtils.maxElementIndex(logPosteriors);
    }

    /**
     * Quals above this value should be capped down to this value (because they are too high)
     * in the base quality score recalibrator
     */
    public static final byte MAX_GATK_USABLE_Q_SCORE = 40;
    private static final double[] logPriorCache = new double[MAX_GATK_USABLE_Q_SCORE + 1];
    static {
        // normal distribution describing P(empiricalQuality - reportedQuality).  Its mean is zero because a priori we expect
        // no systematic bias in the reported quality score
        final double mean = 0.0;
        final double sigma = 0.5;   // with these parameters, deltas can shift at most ~20 Q points
        final NormalDistribution gaussian = new NormalDistribution(null, mean, sigma);

        for ( int i = 0; i <= MAX_GATK_USABLE_Q_SCORE; i++ ) {
            logPriorCache[i] = gaussian.logDensity(i);
        }
    }

    @VisibleForTesting
    protected static double getLogPrior(final double qualityScore, final double priorQualityScore) {
        final int difference = Math.min(Math.abs((int) (qualityScore - priorQualityScore)), MAX_GATK_USABLE_Q_SCORE);
        return logPriorCache[difference];
    }

    private static final long MAX_NUMBER_OF_OBSERVATIONS = Integer.MAX_VALUE - 1;

    /**
     * Given:
     * - n, the number of observations,
     * - k, the number of sequencing errors,
     * - p, the probability of error, encoded as the quality score.
     *
     * Return the binomial probability Bin(k|n,p).
     *
     * The method handles the case when the counts of type long are higher than the maximum allowed integer value,
     * Integer.MAX_VALUE = (2^31)-1 ~= 2*10^9, since the library we use for binomial probability expects integer input.
     *
     */
    @VisibleForTesting
    protected static double getLogBinomialLikelihood(final double qualityScore, long nObservations, long nErrors) {
        if ( nObservations == 0 )
            return 0.0;

        // the binomial code requires ints as input (because it does caching).  This should theoretically be fine because
        // there is plenty of precision in 2^31 observations, but we need to make sure that we don't have overflow
        // before casting down to an int.
        if ( nObservations > MAX_NUMBER_OF_OBSERVATIONS ) {
            // we need to decrease nErrors by the same fraction that we are decreasing nObservations
            final double fraction = (double)MAX_NUMBER_OF_OBSERVATIONS / (double)nObservations;
            nErrors = Math.round((double) nErrors * fraction);
            nObservations = MAX_NUMBER_OF_OBSERVATIONS;
        }

        // this is just a straight binomial PDF
        final double logLikelihood = MathUtils.logBinomialProbability((int) nObservations, (int) nErrors, QualityUtils.qualToErrorProb(qualityScore));
        return ( Double.isInfinite(logLikelihood) || Double.isNaN(logLikelihood) ) ? -Double.MAX_VALUE : logLikelihood;
    }
}
