package org.broadinstitute.hellbender.utils.dragstr;

import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;

/**
 * Tool to figure out the period and repeat-length (in units) of STRs in a reference sequence.
 * <p>
 *    The period and repeat-length of a reference sequence position is determined as follow:
 *    The STR unit is determined only by the sequence from that base onwards.
 * </p>
 * <p>
 *     If the backward base sequence contains additional copies of that unit these are added to the repeat-length.
 * </p>
 * <p>
 *     However a larger repeat-length for a different STR unit upstream would effectively being ignored. Usually this
 *     rule does not come into play but there is an exception every now and then.
 * </p>
 * <p>
 *     For example:
 *     <pre>
 *         ...ATACACACACACACACAC^ACTACTACTACTGAAGA....
 *     </pre>
 *     Where (^) denotes the boundary between down and up-stream from the position of interest (in this case the A that follows)
 * </p><p>
 *     In this case the would take 'ACT' with 4 repeats although 'AC' has plenty more repeats upstream
 * </p>
 * <p>
 *     All sites' period and forward repeat-length are determined in a single pass thru the sequence (O(L * MaxPeriod)).
 *     The backward additional unit count is calculated on demand.
 * </p>
 */
public final class DragstrReferenceAnalyzer {

    private final int start;
    private final int end;
    private final byte[] bases;
    private final int[] period;
    private final int[] forwardRepeats;

    private DragstrReferenceAnalyzer(final byte[] bases, final int start, final int end, final int[] period, final int[] repeats) {
        this.start = start;
        this.end = end;
        this.bases = bases;
        this.period = period;
        this.forwardRepeats = repeats;
    }

    /**
     * Returns the number of consecutive repeated units for a position
     * @param position the target position in the loaded base sequene.
     * @return 1 or greater.
     */
    public int repeatLength(final int position) {
        // we get the forward repeat count:
        int result = lookup(position, forwardRepeats);
        // and then we need to add any further up-stream repeats for that unit.
        final int period = this.period[position - start];
        for (int i = position - 1, j = position + period, k = period; i >= 0; i--) {
            if (bases[i] != bases[--j]) {
                break;
            } else if (--k == 0) { // here we count complete unit matches (k reaches 0).
                k = period;
                result++;
            }
        }
        return result;
    }

    /**
     * Return a new array with a copy of the repeat unit bases starting at a particular position
     * @param position the query position.
     * @return never {@code null}.
     */
    public byte[] repeatUnit(final int position) {
        final int length = lookup(position, period);
        return Arrays.copyOfRange(bases, position, position + length);
    }

    public String repeatUnitAsString(final int position) {
        return new String(repeatUnit(position));
    }

    /**
     * Returns the STR period at a given position.ß
     * @param position the target position
     * @return 1 or greater.
     */
    public int period(final int position) {
        return lookup(position, period);
    }

    private int lookup(final int position, final int[] array) {
        final int offset = position - start;
        if (offset >= 0 && position <= end) {
            return array[offset];
        } else {
            throw new IllegalArgumentException("postion " + position + " is outside bounds");
        }
    }

    public static DragstrReferenceAnalyzer of(final byte[] bases, final int start, final int end, final int maxPeriod) {
        Utils.nonNull(bases, "the input bases cannot be null");
        if (end < start || start < 0 || end > bases.length) {
            throw new IndexOutOfBoundsException("bad indexes " + start  + " " + end + " " + bases.length);
        } else if (start >= end) {
            return new DragstrReferenceAnalyzer(bases, start, end, ArrayUtils.EMPTY_INT_ARRAY, ArrayUtils.EMPTY_INT_ARRAY);
        }

        // Period one's processing code is trivial, so do it separately and we initialize the best period and repeat
        // at the same time.
        final int[] repeats = processPeriodOne(bases, start, end);
        final int[] periods = new int[end - start];
        Arrays.fill(periods, 1);
        for (int period = 2; period <= maxPeriod; period++) {
            if (bases.length  < period << 1) {
                break;
            }
            int position, remainingToMatchUnit, carryBack, rightMargin, positionPlusPeriod, resultArrayOffset;

            // We first find the right-margin enclosing the last consecutive repeat in the input base sequence that would match
            // a prospective repeat unit within the target interval [start,end):
            rightMargin = Math.min(end + period, bases.length) - 1;
            final int rightMostStart = position = rightMargin - period;
            // reaminingToMatchUnit is a countdown reaching 0 when we have match yet another full repeat.
            remainingToMatchUnit = period;
            // carryBack would count the number of matching repeats outside the target interval
            carryBack = 1;
            for (; rightMargin < bases.length; rightMargin++) {
                        if (bases[position++] != bases[rightMargin]) {
                    break;
                } else if (--remainingToMatchUnit == 0) {
                    carryBack++;
                    remainingToMatchUnit = period;
                }
            }

            // No we work our way backwards carrying on number of consecutive matching units
            // and updating.
            // carryBack and remainingToMatchUnit has been updated in the forward pass so they correspond to the value
            // at rightMostStart.
            if ((resultArrayOffset = rightMostStart - start) >= 0 && carryBack > repeats[resultArrayOffset]) {
                repeats[resultArrayOffset] = carryBack;
                periods[resultArrayOffset] = period;
            }

            boolean inTragetRange = false;
            for (position = rightMostStart - 1, positionPlusPeriod = position + period,
                 resultArrayOffset = position - start; position >= start; position--, positionPlusPeriod--, resultArrayOffset--) {
                if (bases[position] == bases[positionPlusPeriod]) {
                    if (--remainingToMatchUnit == 0) { // have we matched yet another unit of length period?
                        ++carryBack;
                        remainingToMatchUnit = period;
                    }
                    if ((inTragetRange |= position < end) && carryBack > repeats[resultArrayOffset]) {
                        repeats[resultArrayOffset] = carryBack;
                        periods[resultArrayOffset] = period;
                    }
                } else {
                    carryBack = 1;
                    remainingToMatchUnit = period;
                }
            }
        }
        return new DragstrReferenceAnalyzer(bases, start, end, periods, repeats);
    }

    /**
     * Special faster code for period == 1.
     */
    private static int[] processPeriodOne(final byte[] bases, final int start, final int end) {
        final int[] repeats = new int[end - start];
        byte last = bases[end - 1];

        // backward phase:
        int rightMargin, position;
        for (rightMargin = end; rightMargin < bases.length && bases[rightMargin] == last; rightMargin++) {};
        int carryBack = rightMargin - end;
        int offset;
        for (position = end - 1, offset = position - start; position >= start; position--, offset--) {
            final byte next = bases[position];
            repeats[offset] = next == last ? ++carryBack : (carryBack = 1);
            last = next;
        }
        return repeats;
    }
}
