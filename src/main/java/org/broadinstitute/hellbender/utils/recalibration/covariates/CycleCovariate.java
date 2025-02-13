package org.broadinstitute.hellbender.utils.recalibration.covariates;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;

/**
 * The Cycle covariate.
 * For ILLUMINA the cycle is simply the position in the read (counting backwards if it is a negative strand read)
 */
public final class CycleCovariate implements Covariate {
    private static final long serialVersionUID = 1L;

    private final int MAXIMUM_CYCLE_VALUE;
    public static final int CUSHION_FOR_INDELS = 4;

    public CycleCovariate(final RecalibrationArgumentCollection RAC){
        this.MAXIMUM_CYCLE_VALUE = RAC.MAXIMUM_CYCLE_VALUE;
    }

    // Used to pick out the covariate's value from attributes of the read
    @Override
    public void recordValues(final GATKRead read, final SAMFileHeader header, final PerReadCovariateMatrix values, final boolean recordIndelValues) {
        final int readLength = read.getLength();
        //Note: duplicate the loop to void checking recordIndelValues on every iteration
        if (recordIndelValues) {
            for (int i = 0; i < readLength; i++) {
                final int substitutionKey = cycleKey(i, read, false, MAXIMUM_CYCLE_VALUE);
                final int indelKey = cycleKey(i, read, true, MAXIMUM_CYCLE_VALUE);
                values.addCovariate(substitutionKey, indelKey, indelKey, i);
            }
        } else {
            for (int i = 0; i < readLength; i++) {
                final int substitutionKey = cycleKey(i, read, false, MAXIMUM_CYCLE_VALUE);
                values.addCovariate(substitutionKey, 0, 0, i);
            }
        }
    }

    @Override
    public String formatKey(final int key){
            return String.format("%d", cycleFromKey(key));
    }

    @Override
    public int keyFromValue(final Object value) {
        return (value instanceof String) ? keyFromCycle(Integer.parseInt((String) value), MAXIMUM_CYCLE_VALUE) : keyFromCycle((Integer) value, MAXIMUM_CYCLE_VALUE);
    }

    @Override
    public int maximumKeyValue() {
        return (MAXIMUM_CYCLE_VALUE << 1) + 1;
    }

    /**
     * Computes the encoded value of CycleCovariate's key for the given position at the read.
     * Uses keyFromCycle to do the encoding.
     * @param baseNumber index of the base to compute the key for
     * @param read the read
     * @param indel is this an indel key or a substitution key?
     * @param maxCycle max value of the base to compute the key for
     *                 (this method throws UserException if the computed absolute value of the cycle number is higher than this value).
     */
    public static int cycleKey(final int baseNumber, final GATKRead read, final boolean indel, final int maxCycle) {
        final boolean isNegStrand = read.isReverseStrand();
        final boolean isSecondInPair = read.isPaired() && read.isSecondOfPair();
        final int readLength = read.getLength();

        final int readOrderFactor = isSecondInPair ? -1 : 1;
        final int increment;
        int cycle;
        if (isNegStrand) {
            cycle = readLength * readOrderFactor;
            increment = -1 * readOrderFactor;
        } else {
            cycle = readOrderFactor;
            increment = readOrderFactor;
        }

        cycle += baseNumber * increment;

        if (!indel) {
            return CycleCovariate.keyFromCycle(cycle, maxCycle);
        }
        final int maxCycleForIndels = readLength - CycleCovariate.CUSHION_FOR_INDELS - 1;
        if (baseNumber < CycleCovariate.CUSHION_FOR_INDELS || baseNumber > maxCycleForIndels) {
            return -1;
        } else {
            return CycleCovariate.keyFromCycle(cycle, maxCycle);
        }
    }

    /**
     * Decodes the cycle number from the key.
     */
    public static int cycleFromKey(final int key) {
        int cycle = key >> 1; // shift so we can remove the "sign" bit
        if ( (key & 1) != 0 ) { // is the last bit set?
            cycle *= -1; // then the cycle is negative
        }
        return cycle;
    }

    /**
     * Encodes the cycle number as a key.
     */
    public static int keyFromCycle(final int cycle, final int maxCycle) {
        // no negative values because values must fit into the first few bits of the long
        int result = Math.abs(cycle);
        if ( result > maxCycle ) {
            throw new UserException("The maximum allowed value for the cycle is " + maxCycle + ", but a larger cycle (" + result + ") was detected.  Please use the --maximum-cycle-value argument (when creating the recalibration table in BaseRecalibrator) to increase this value (at the expense of requiring more memory to run)");
        }

        result <<= 1; // shift so we can add the "sign" bit
        if ( cycle < 0 ) {
            result++; // negative cycles get the lower-most bit set
        }
        return result;
    }
}