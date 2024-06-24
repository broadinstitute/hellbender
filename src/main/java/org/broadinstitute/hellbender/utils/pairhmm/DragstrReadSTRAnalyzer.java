package org.broadinstitute.hellbender.utils.pairhmm;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Arrays;

/**
 * Utility to find short-tandem-repeats on read sequences.
 * <p>
 *     For each position of a read base sequence indicates the overlapping repeat unit with the largest repeat length.
 *     In case of a tie the shorter unit takes preference.
 * </p>
 * <p>
 *     To be a bit more specific, for a position <i>i</i> and repeat unit period length <i>p</i> the algorithm considers
 *     all the <i>p-kmer</i> overlapping units + the one starting in the next position <i>i + 1</i>.
 * </p>
 * <p>
 *     For example:
 *     <pre>
 *         seq:    AAATATATATAAA
 *         best-p: 1222222222211
 *         rep:    3444444444433
 *     </pre>
 *     For most part the values for the best period are the intuitive choice except for the second A with a best
 *     period is 2 and repeat length 4. Although it is part of the more eye catching triplet "AAA" that second A is
 *     right before the "AT" unit that has 4 repeats downstream.
 * </p>
 * <p>
 *     Notice that there is a peculiarity of this algorithm in that for periods larger than 1 the very first base
 *     has a repeat length that is actually 1 less than the expected value by the description above.
 * </p>
 * <p>
 *     For example:
 *     <pre>
 *         seq:    ATATATATATATAT
 *         best-p  22222222222222
 *         rep:    67777777777777
 *     </pre>
 * </p>
 * <p>
 *     This behavior is exhibit by DRAGEN, so here we copied it.
 * </p>
 */
public final class DragstrReadSTRAnalyzer {

    private final int[][] repeatsByPeriodAndPosition;
    private final int[] periodWithMostRepeats;
    private final int maxPeriod;
    private final int seqLength;

    private DragstrReadSTRAnalyzer(final int seqLength, final int[][] repeatsByPeriodAndPosition,
                                   final int[] periodWithMostRepeats, final int maxPeriod) {
        this.repeatsByPeriodAndPosition = repeatsByPeriodAndPosition;
        this.periodWithMostRepeats = periodWithMostRepeats;
        this.maxPeriod = maxPeriod;
        this.seqLength = seqLength;
    }

    /**
     * Creates an analyzer for a read.
     *
     * @param read the target read to analyze.
     * @param maxPeriod the maximum period or repeat unit length to consider.
     */
    public static DragstrReadSTRAnalyzer of(final GATKRead read, final int maxPeriod) {
        Utils.nonNull(read, "the input read cannot be null");
        final byte[] bases = Utils.nonNull(read.getBasesNoCopy(), "the input read must have bases");
        return of(bases, maxPeriod);
    }

    @VisibleForTesting
    static DragstrReadSTRAnalyzer of(final byte[] bases, final int maxPeriod) {
        ParamUtils.isPositive(maxPeriod, "the input max period must be 1 or greater");
        final int seqLength = bases.length;
        final int[][] repeatsByPeriodAndPosition = new int[maxPeriod][seqLength];
        final int[] periodWithMostRepeats = new int[seqLength];
        analyzeBases(bases, seqLength, repeatsByPeriodAndPosition, periodWithMostRepeats, maxPeriod);
        return new DragstrReadSTRAnalyzer(seqLength, repeatsByPeriodAndPosition, periodWithMostRepeats, maxPeriod);
    }

    /**
     * Returns the maximum number of repeats for any repeat unit of a particular period size at
     * that position.
     * @param position the query position. 0-based.
     * @param period the repeat unit length to consider.
     * @return 1 or greater.
     */
    public int numberOfRepeats(final int position, final int period) {
        if (period <= 0 || period > maxPeriod) {
            return 0;
        } else if (position < 0 || position >= seqLength) {
            throw new IllegalArgumentException("cannot query outside requested boundaries");
        } else {
            return repeatsByPeriodAndPosition[period - 1][position];
        }
    }

    /**
     * Returns the period length with the larger repeat length in a particular positions.
     * <p>
     *     In case of a tie, the shorter period takes preference.
     * </p>
     * @param position the query position. 0-based.
     * @return any value between 1  and {@code maxPeriod}.
     */
    public int mostRepeatedPeriod(final int position) {
        if (position >= 0 && position < seqLength) {
            return periodWithMostRepeats[position];
        } else {
            throw new IllegalArgumentException(String.format("cannot query outside requested boundaries [%d, %d): %d", 0, seqLength, position));
        }
    }

    /**
     * Returns the number of repeat units for the period with largest count.
     * <p>
     *     In case of a tie the shorter period takes preference.
     * </p>
     * <p>
     *     In short, this is equivalent to:
     *     <pre>
     *         numberOfRepeats(mostRepeatedPeriod(pos), pos);
     *     </pre>
     *     In case of a tie, the shorter period takes preference.
     * </p>
     * @param position the query position. 0-based.
     * @return any value between 1  and {@code maxPeriod}.
     */
    public int numberOfMostRepeats(final int position) {
        if (position >= 0 && position < seqLength) {
            return repeatsByPeriodAndPosition[periodWithMostRepeats[position] - 1][position];
        } else {
            throw new IllegalArgumentException("cannot query outside requested boundaries");
        }
    }

    /**
     * Does all the calculations on the base sequence.
     *
     * @param bases the sequence of bases of the read under analysis.
     * @param seqLength the length of the base sequence, its always equal to bases.length, here we are just reusing the precalculated value.
     * @param repeatsByPeriodAndPosition the output matrix period x position where the analysis will populate the repeat number
     *                                   for every period and read position.
     * @param periodWithMostRepeats the output array where we indicate the period with the maximum number of repeats on that position
     *                              ; ties are broken in favour of the shorter period.
     * @param maxPeriod, maximum period to be considered.
     */
    private static void analyzeBases(final byte[] bases, final int seqLength, final int[][] repeatsByPeriodAndPosition,
                             final int[] periodWithMostRepeats, final int maxPeriod) {

        calculateRepeatsForPeriodOne(bases, seqLength, repeatsByPeriodAndPosition[0]); // simpler and faster code for period-length == 0.
        // runLengthBuffer and heap are tmp memory reused between iterations for period 2 and above.
        final int[] runLengthBuffer = new int[seqLength + 1];
        final int[] heap = new int[maxPeriod + 1];
        for (int period = 2; period <= maxPeriod; period++) {
            calculateRepeatsForPeriodTwoAndAbove(bases, seqLength, period, runLengthBuffer, heap, repeatsByPeriodAndPosition[period - 1]);
        }

        // finally we update the periodWithMostRepeats array
        // we the best period at each position.
        Arrays.fill(periodWithMostRepeats, 0, seqLength, 1);
        final int[] mostRepeats = repeatsByPeriodAndPosition[0].clone();
        for (int periodIndex = 1; periodIndex < maxPeriod; periodIndex++) {
            final int periodLength = periodIndex + 1;
            final int[] periodValues = repeatsByPeriodAndPosition[periodIndex];
            for (int position = 0; position < seqLength; position++) {
                final int repeats = periodValues[position];
                if (repeats > mostRepeats[position]) {
                    mostRepeats[position] = repeats;
                    periodWithMostRepeats[position] = periodLength;
                }
            }
        }
    }


    private static void calculateRepeatsForPeriodTwoAndAbove(final byte[] bases, final int seqLength, final int period,
                                                             final int[] runLengthBuffer, final int[] heap, final int[] output) {

       // very small sequence (less than period) then repeats are all zero:2
        if (seqLength < period) {
            Arrays.fill(output, 0, seqLength, 0);
            return;
        }

        int position, matchedCycles, cycleIndex, positionPlusPeriod;

        // suffixes of zero run-lengths at the end of the sequence (period - 1).
        for (position = seqLength - 1, cycleIndex = period; cycleIndex > 1; position--, cycleIndex--) {
            runLengthBuffer[position] = 0;
        }

        // backward phase.
        // at the end of this phase runLength[x] will have the total number of equal repeats downstream with length starting at position X
        // inclusive.
        int carryBack = runLengthBuffer[position--] = 1; //prevValue holds the num of repeats reported in the previous (+1) position.
        for (positionPlusPeriod = position + period, matchedCycles = 0; position >= 0; position--, positionPlusPeriod--) {
            if (bases[position] == bases[positionPlusPeriod]) { // we keep matching repeat unit bases.
                if (++matchedCycles == period) { // we got a new full repeat matched so the run length increases:
                    runLengthBuffer[position] = ++carryBack;
                    matchedCycles = 0; // we reset the match-run-length to 0 as we start over (a new repeat).
                } else { // if we haven't completed a repeat we simply copy the run length from the +1 position base.
                    runLengthBuffer[position] = carryBack;
                }
            } else { // we bump into a mismatch that ends the run.
                carryBack = runLengthBuffer[position] = 1;
                matchedCycles = 0; // we reset the match-run-length to 0.
            }
        }
        // Now we propagate forward the number of equivalent repeats up stream:
        // So at the end runLength[X] == the repeats (up and down-stream) of the
        // unit of length period that starts at X.
        // The outside for-loop iterates over different repeat unit start offsets (here named cycles):
        for (cycleIndex = 0; cycleIndex < period; cycleIndex++) {
            // The left most repeated unit runLength[i] contains the actual run length for all the matching
            // units to the right. We copy that value forward to the other run-length units.
            // We do this by iterating over consecutive repeat runs.
            for (position = cycleIndex; position < seqLength; position += period) {
                final int totalRunLength = runLengthBuffer[position];
                for (int repeatInRun = 1; repeatInRun < totalRunLength; repeatInRun++) {
                    runLengthBuffer[position += period] = totalRunLength;
                }
            }
        }

        // NOTE: this line below is added solely to replicate an oddity of DRAGEN's algorithm that discounts one repeat
        // only at the beginning of the sequence for period 2 or above (e.g. ^CACATG would yield periods 122211 repeats
        // 12221 when the natural/intuitive solution is period 222211 reps 222211). This is true for the original Matlab
        // and latest DRAGEN so we must add this here:
        if (runLengthBuffer[0] > 1) runLengthBuffer[0]--;
        // Finally we need to propagate the best run-lengths to neighbor positions
        // so that a position has the maximum run-length of all the possible
        // units that contain the position + the unit that follow up-stream.
        final int heapSize = period + 1;
        Arrays.fill(heap, 0, heapSize, 0);
        int currentMax = heap[0] = runLengthBuffer[0];
        // the first period length prefix is trivial:
        // a position's max run-length is the largest seen so far and the following position's
        // We use this opportunity to populate the
        // run-length max-heap
        final int stop0 = Math.min(period, seqLength - 1);
        for (position = 0; position < stop0; ) {
            output[position] = currentMax = Math.max(currentMax, heap[++position] = runLengthBuffer[position]);
            fixHeap(heap, position, heapSize);
        }

        runLengthBuffer[seqLength] = 1; // we add this 1 as the runLength after the last sequence position which in fact does not exist.
                                       // this allows us to save a conditional to avoid a index-out-of-range.
        for (int outPosition = 0; position < seqLength;) {
            final int valueOut = runLengthBuffer[outPosition++]; // value leaving the heap.
            final int valueIn = runLengthBuffer[++position]; // value entering the heap.
            if (valueIn != valueOut) { // these two are the same (often the case) we don't need to do a heap updating at all.
                final int fixHeapIdx = ArrayUtils.indexOf(heap, valueOut); // O(n) search although usually n is very small.
                heap[fixHeapIdx] = valueIn;
                fixHeap(heap, fixHeapIdx, heapSize);
                currentMax = heap[0];
            }
            output[position - 1] = currentMax; // we use the heap's max as the final run-length for the prev position.
        }
    }

    private static void fixHeap(final int[] heap, final int idx, final int heapSize) {
        if (idx == 0 || heap[(idx - 1) >> 1] > heap[idx]) {
            fixHeapDown(heap, idx, heapSize);
        } else {
            fixHeapUp(heap, idx);
        }
    }

    private static void fixHeapUp(final int[] heap, int idx) {
        final int value = heap[idx];
        do {
            final int upIdx = (idx - 1) >> 1;
            final int upValue = heap[upIdx];
            if (upValue >= value) {
                break;
            } else {
                heap[idx] = upValue;
                idx = upIdx;
            }
        } while (idx > 0);
        heap[idx] = value;
    }

    /**
     * Check that the content of the heap is correct given its size.
     *
     * Used for debugging, it might be useful in the future; better to keep it
     * around.
     */
    @SuppressWarnings("unused")
    private static boolean checkHeap(final int[] heap, final int heapSize) {
        for (int i = 1; i < heapSize; i++) {
            if (heap[(i - 1) >> 1] < heap[i]) {
                return false;
            }
        }
        return true;
    }

    private static void fixHeapDown(final int[] heap, int idx, final int heapSize) {
        final int value = heap[idx];
        while (true) {
            final int idxRight = (idx + 1) << 1;
            final int rightValue = idxRight < heapSize ? heap[idxRight] : -1;
            final int idxLeft = idxRight - 1;
            final int leftValue = idxLeft < heapSize ? heap[idxLeft] : -1;
            if (rightValue > value) {
                if (rightValue > leftValue) {
                    heap[idx] = rightValue;
                    idx = idxRight;
                } else {
                    heap[idx] = leftValue;
                    idx = idxLeft;
                }
            } else if (leftValue > value) {
                heap[idx] = leftValue;
                idx = idxLeft;
            } else {
                break;
            }
        }
        heap[idx] = value;
    }

    /**
     * Faster and simpler implementation of {@link #calculateRepeatsForPeriodTwoAndAbove}
     * for period == 1.
     * <p>
     *     In this case we can use the output array as the run-length-buffer.
     * </p>
     */
    private static void calculateRepeatsForPeriodOne(final byte[] bases, final int seqLength, final int[] output) {
        final int rightMargin = seqLength - 1;
        byte last = bases[rightMargin];

        // backward phase:
        int carryBack = output[rightMargin] = 1;
        for (int position = rightMargin - 1; position >= 0; position--) {
            final byte next = bases[position];
            output[position] = next == last ? ++carryBack : (carryBack = 1);
            last = next;
        }
        // forward phase:
        // last = sequence[0]; // already true.
        int prevRunLength = output[0];
        for (int position = 1; position <= rightMargin; position++) {
            final byte next = bases[position];
            if (next == last) {
                output[position] = prevRunLength;
            } else {
                final int thisRunLength = output[position];
                if (prevRunLength < thisRunLength) { // so we propagate this position run-length to the prev position if longer.
                    output[position - 1] = thisRunLength;
                }
                last = next;
                prevRunLength = thisRunLength;
            }
        }
    }
}
