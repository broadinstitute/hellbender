package org.broadinstitute.hellbender.utils.recalibration.covariates;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.clipping.ClippingRepresentation;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;

/**
 * The read bases preceding and including the base in question (as opposed to reference bases).
 */
public final class ContextCovariate implements Covariate {
    private static final long serialVersionUID = 1L;
    private static final Logger logger = LogManager.getLogger(ContextCovariate.class);

    private final int mismatchesContextSize;
    private final int indelsContextSize;

    private final int mismatchesKeyMask;
    private final int indelsKeyMask;

    private static final int LENGTH_BITS = 4;
    private static final int LENGTH_MASK = 15;


    // the maximum context size (number of bases) permitted; we need to keep the leftmost base free so that values are
    // not negative and we reserve 4 more bits to represent the length of the context; it takes 2 bits to encode one base.
    private static final int MAX_DNA_CONTEXT = 13;
    private final byte lowQualTail;

    public static final int UNKNOWN_OR_ERROR_CONTEXT_CODE = -1;

    public ContextCovariate(final RecalibrationArgumentCollection RAC){
        mismatchesContextSize = RAC.MISMATCHES_CONTEXT_SIZE;
        indelsContextSize = RAC.INDELS_CONTEXT_SIZE;
        logger.debug("\t\tContext sizes: base substitution model " + mismatchesContextSize + ", indel substitution model " + indelsContextSize);

        if (mismatchesContextSize > MAX_DNA_CONTEXT) {
            throw new CommandLineException.BadArgumentValue("mismatches_context_size", String.format("context size cannot be bigger than %d, but was %d", MAX_DNA_CONTEXT, mismatchesContextSize));
        }
        if (indelsContextSize > MAX_DNA_CONTEXT) {
            throw new CommandLineException.BadArgumentValue("indels_context_size", String.format("context size cannot be bigger than %d, but was %d", MAX_DNA_CONTEXT, indelsContextSize));
        }

        lowQualTail = RAC.LOW_QUAL_TAIL;

        if (mismatchesContextSize <= 0 || indelsContextSize <= 0) {
            throw new CommandLineException(String.format("Context size must be positive. Mismatches: %d Indels: %d", mismatchesContextSize, indelsContextSize));
        }

        mismatchesKeyMask = createMask(mismatchesContextSize);
        indelsKeyMask = createMask(indelsContextSize);
    }

    @Override
    public void recordValues(final GATKRead read, final SAMFileHeader header, final PerReadCovariateMatrix values, final boolean recordIndelValues) {

        final int originalReadLength = read.getLength();

        // store the original bases and then write Ns over low quality ones
        final byte[] strandedClippedBases = getStrandedClippedBytes(read, lowQualTail);  //Note: this makes a copy of the read

        // Note: we're using a non-standard library here because boxing came up on profiling as taking 20% of time in applyBQSR.
        // IntList avoids boxing
        // This list has the same length as the read and contains the n-base context at each position.
        final IntList nBasePairContextAtEachCycle = getReadContextAtEachPosition(strandedClippedBases, mismatchesContextSize, mismatchesKeyMask);

        final int readLengthAfterClipping = strandedClippedBases.length;

        // this is necessary to ensure that we don't keep historical data in the ReadCovariates values
        // since the context covariate may not span the entire set of values in read covariates
        // due to the clipping of the low quality bases
        if ( readLengthAfterClipping != originalReadLength) {
            // don't bother zeroing out if we are going to overwrite the whole array
            for ( int i = 0; i < originalReadLength; i++ ){
                // this base has been clipped off, so zero out the covariate values here
                values.addCovariate(0, 0, 0, i);
            }
        }

        final boolean negativeStrand = read.isReverseStrand();

        //Note: duplicated the loop to avoid checking recordIndelValues on each iteration
        if (recordIndelValues) {
            final IntList indelKeys = getReadContextAtEachPosition(strandedClippedBases, indelsContextSize, indelsKeyMask);
            for (int i = 0; i < readLengthAfterClipping; i++) {
                final int readOffset = getStrandedOffset(negativeStrand, i, readLengthAfterClipping);
                final int indelKey = indelKeys.getInt(i);
                values.addCovariate(nBasePairContextAtEachCycle.getInt(i), indelKey, indelKey, readOffset);
            }
        } else {
            for (int i = 0; i < readLengthAfterClipping; i++) {
                final int readOffset = getStrandedOffset(negativeStrand, i, readLengthAfterClipping);
                values.addCovariate(nBasePairContextAtEachCycle.getInt(i), 0, 0, readOffset);
            }
        }
    }

    /**
     * Helper method: computes the correct offset to use in computations of covariate values.
     * @param isNegativeStrand is the read on the negative strand
     * @param offset 0-based index of the base in the read
     * @param readLength length of the read
     * @return
     */
    public static int getStrandedOffset(final boolean isNegativeStrand, final int offset, final int readLength) {
        return isNegativeStrand ? (readLength - offset - 1) : offset;
    }

    /**
     * Given a read, clips low quality ends (by overwriting with N) and returns the underlying bases, after
     * reverse-complementing for negative-strand reads.
     * @param read the read
     * @param lowQTail every base quality lower than or equal to this in the tail of the read will be replaced with N.
     * @return bases of the read (Could be an empty array if all bases are below lowQTail).
     */
    @VisibleForTesting
    static byte[] getStrandedClippedBytes(final GATKRead read, final byte lowQTail) {

        // Write N's over the low quality tail of the reads to avoid adding them into the context
        final GATKRead clippedRead = ReadClipper.clipLowQualEnds(read, lowQTail, ClippingRepresentation.WRITE_NS);

        final byte[] bases = clippedRead.getBases();
        if (read.isReverseStrand()) {
            return BaseUtils.simpleReverseComplement(bases);
        } else {
            return bases;
        }
    }

    @Override
    public String formatKey(final int key) {
        if (key == -1) // this can only happen in test routines because we do not propagate null keys to the csv file
            return null;

        return contextFromKey(key);
    }

    @Override
    public int keyFromValue(final Object value) {
        return keyFromContext((String) value);
    }

    private static int createMask(final int contextSize) {
        int mask = 0;
        // create 2*contextSize worth of bits
        for (int i = 0; i < contextSize; i++) {
            mask = (mask << 2) | 3;
        }
        // shift 4 bits to mask out the bits used to encode the length
        return mask << LENGTH_BITS;
    }

    /**
     * For each position of the read, calculate the n-base-pair *read* base context (as opposed to the reference context).
     *
     * For example, for the read [AGCTG], return the list
     *   [-1, "AG", "GC", "CT", "TG" ]
     * with each string context encoded as an integer.
     *
     * @param bases       the bases in the read to build the context from
     * @param contextSize context size to use building the context
     * @param mask        mask for pulling out just the context bits
     *
     * @return a list that has the same length as the read and contains the (preceding) n-base context at each position.
     *
     */
    private static IntList getReadContextAtEachPosition(final byte[] bases, final int contextSize, final int mask) {

        final int readLength = bases.length;

        //Note: we use a specialized collection to avoid the cost of boxing and unboxing that otherwise comes up on the profiler.
        final IntList keys = new IntArrayList(readLength);

        // the first contextSize-1 bases will not have enough previous context
        for (int i = 1; i < contextSize && i <= readLength; i++) {
            keys.add(UNKNOWN_OR_ERROR_CONTEXT_CODE);
        }

        if (readLength < contextSize) {
            return keys;
        }

        final int newBaseOffset = 2 * (contextSize - 1) + LENGTH_BITS;

        // get (and add) the key for the context starting at the first base
        int currentKey = keyFromContext(bases, 0, contextSize);
        keys.add(currentKey);

        // if the first key was -1 then there was an non-ACGT in the context; figure out how many more consecutive contexts it affects
        int currentNPenalty = 0;
        if (currentKey == -1) {
            currentKey = 0;
            currentNPenalty = contextSize - 1;
            int offset = newBaseOffset;
            int baseIndex;
            while ((baseIndex = BaseUtils.simpleBaseToBaseIndex(bases[currentNPenalty])) != -1) {
                currentKey |= (baseIndex << offset);
                offset -= 2;
                currentNPenalty--;
            }
        }

        for (int currentIndex = contextSize; currentIndex < readLength; currentIndex++) {
            final int baseIndex = BaseUtils.simpleBaseToBaseIndex(bases[currentIndex]);
            if (baseIndex == -1) { // ignore non-ACGT bases
                currentNPenalty = contextSize;
                currentKey = 0; // reset the key
            } else {
                // push this base's contribution onto the key: shift everything 2 bits, mask out the non-context bits, and add the new base and the length in
                currentKey = (currentKey >> 2) & mask;
                currentKey |= (baseIndex << newBaseOffset);
                currentKey |= contextSize;
            }

            if (currentNPenalty == 0) {
                keys.add(currentKey);
            } else {
                currentNPenalty--;
                keys.add(-1);
            }
        }

        return keys;
    }

    public static int keyFromContext(final String dna) {
        return keyFromContext(dna.getBytes(), 0, dna.length());
    }

    /**
     * Creates a int representation of a given dna string.
     *
     * @param dna    the dna sequence
     * @param start  the start position in the byte array (inclusive)
     * @param end    the end position in the array (exclusive)
     * @return the key representing the dna sequence
     */
    private static int keyFromContext(final byte[] dna, final int start, final int end) {

        int key = end - start;
        int bitOffset = LENGTH_BITS;
        for (int i = start; i < end; i++) {
            final int baseIndex = BaseUtils.simpleBaseToBaseIndex(dna[i]);
            if (baseIndex == -1) { // ignore non-ACGT bases
                return -1;
            }
            key |= (baseIndex << bitOffset);
            bitOffset += 2;
        }
        return key;
    }

    /**
     * Converts a key into the dna string representation.
     *
     * @param key    the key representing the dna sequence
     * @return the dna sequence represented by the key
     */
    public static String contextFromKey(final int key) {
        if (key < 0) {
            throw new GATKException("dna conversion cannot handle negative numbers. Possible overflow?");
        }

        final int length = key & LENGTH_MASK; // the first bits represent the length (in bp) of the context
        int mask = 48; // use the mask to pull out bases
        int offset = LENGTH_BITS;

        final StringBuilder dna = new StringBuilder(length);
        for (int i = 0; i < length; i++) {
            final int baseIndex = (key & mask) >> offset;
            dna.append((char)BaseUtils.baseIndexToSimpleBase(baseIndex));
            mask <<= 2; // move the mask over to the next 2 bits
            offset += 2;
        }

        return dna.toString();
    }

    @Override
    public int maximumKeyValue() {
        // the maximum value is T (11 in binary) for each base in the context
        final int length = Math.max(mismatchesContextSize, indelsContextSize);  // the length of the context
        int key = length;
        int bitOffset = LENGTH_BITS;
        for (int i = 0; i <length ; i++) {
            key |= (3 << bitOffset);
            bitOffset += 2;
        }
        return key;
    }
}
