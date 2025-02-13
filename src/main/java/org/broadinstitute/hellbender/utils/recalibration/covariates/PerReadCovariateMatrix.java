package org.broadinstitute.hellbender.utils.recalibration.covariates;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.recalibration.EventType;

/**
 * The object that holds BQSR covarites for a read as a matrix (table) of shape ( read length ) x ( num covariates ).
 * The covariates are encoded as their integer keys.
 *
 * Even though it's not a matrix in the linear algebra sense, we call it a matrix rather than table to
 * differentiate it from the recal table in BQSR.
 */
public final class PerReadCovariateMatrix {
    private static final Logger logger = LogManager.getLogger(PerReadCovariateMatrix.class);

    /**
     * This is where we store the pre-read covariates, also indexed by (event type) and (read position).
     * Thus the array has shape { event type } x { read position (aka cycle) } x { covariate }.
     * For instance, { covariate } is by default 4-dimensional (read group, base quality, context, cycle).
     */
    private final int[][][] covariates;

    /**
     * The index of the current covariate, used by addCovariate
     */
    private int currentCovariateIndex = 0;

    /**
     * Use an LRU cache to keep cache of keys (int[][][]) arrays for each read length we've seen.
     * The cache allows us to avoid the expense of recreating these arrays for every read.  The LRU
     * keeps the total number of cached arrays to less than LRU_CACHE_SIZE.
     */
    public PerReadCovariateMatrix(final int readLength, final int numberOfCovariates, final CovariateKeyCache keysCache) {
        Utils.nonNull(keysCache);
        final int[][][] cachedKeys = keysCache.get(readLength);
        if ( cachedKeys == null ) {
            if ( logger.isDebugEnabled() ) logger.debug("Keys cache miss for length " + readLength + " cache size " + keysCache.size());
            covariates = new int[EventType.values().length][readLength][numberOfCovariates];
            keysCache.put(readLength, covariates);
        } else {
            covariates = cachedKeys;
        }
    }

    public void setCovariateIndex(final int index) {
        currentCovariateIndex = index;
    }

    /**
     * Update the keys for mismatch, insertion, and deletion for the current covariate at read offset
     *
     * NOTE: no checks are performed on the number of covariates, for performance reasons.  If the count increases
     * after the keysCache has been accessed, this method will throw an ArrayIndexOutOfBoundsException.  This currently
     * only occurs in the testing harness, and we don't anticipate that it will become a part of normal runs.
     *
     * @param mismatch the mismatch key value
     * @param insertion the insertion key value
     * @param deletion the deletion key value
     * @param readOffset the read offset, must be >= 0 and <= the read length used to create this ReadCovariates
     */
    public void addCovariate(final int mismatch, final int insertion, final int deletion, final int readOffset) {
        covariates[EventType.BASE_SUBSTITUTION.ordinal()][readOffset][currentCovariateIndex] = mismatch;
        covariates[EventType.BASE_INSERTION.ordinal()][readOffset][currentCovariateIndex] = insertion;
        covariates[EventType.BASE_DELETION.ordinal()][readOffset][currentCovariateIndex] = deletion;
    }

    /**
     * Get the keys for all covariates at read position for error model
     *
     * @param readPosition
     * @param errorModel
     * @return
     */
    public int[] getCovariatesAtOffset(final int readPosition, final EventType errorModel) {
        return covariates[errorModel.ordinal()][readPosition];
    }
    // tsato: this keySet language needs to go. "getTableForErrorType". Must replace all "getKey"
    public int[][] getMatrixForErrorModel(final EventType errorModel) {
        return covariates[errorModel.ordinal()];
    }

    // ----------------------------------------------------------------------
    //
    // routines for testing
    //
    // ----------------------------------------------------------------------

    protected int[][] getMismatchMatrix() { return getMatrixForErrorModel(EventType.BASE_SUBSTITUTION); }
    protected int[][] getInsertionMatrix() { return getMatrixForErrorModel(EventType.BASE_INSERTION); }
    protected int[][] getDeletionMatrix() { return getMatrixForErrorModel(EventType.BASE_DELETION); }

    protected int[] getMismatchCovariatesAtOffset(final int readPosition) {
        return getCovariatesAtOffset(readPosition, EventType.BASE_SUBSTITUTION);
    }

    protected int[] getInsertionCovariatesAtOffset(final int readPosition) {
        return getCovariatesAtOffset(readPosition, EventType.BASE_INSERTION);
    }

    protected int[] getDeletionCovariatesAtOffset(final int readPosition) {
        return getCovariatesAtOffset(readPosition, EventType.BASE_DELETION);
    }
}
