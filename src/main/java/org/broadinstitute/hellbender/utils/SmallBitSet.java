package org.broadinstitute.hellbender.utils;

import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAlleleCounts;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeIndexCalculator;

import java.util.Collection;
import java.util.Iterator;

/**
 * This class is very much like the standard library BitSet class in java.util, but instead of using
 * a long[] array it uses single (32-bit) ints, which are much faster.  This limits the capacity but
 * for small sets it makes sense.
 *
 * This class uses the binary representation of sets, in which set {i1, i2, . . .in} maps to the
 * integer 2^(i_1) + 2^(i_2) . . . + 2^(i_n).
 *
 * This is efficient because 2^n = 1 << n is simply the integer 1 with n bit-shifts, and because
 * addition of different powers of 2 is equivalent to the bitwise OR.
 *
 * Equivalently, the singleton set n in binary is represented as the integer where the nth bit is 1
 * and all others are zero, and the set {i1, i2, . . .in} is represented as the integer where
 * bits i1, i2. . .in are 1 and all others are zero.
 *
 */
public class SmallBitSet {
    public static final int MAX_ELEMENTS = 30;
    private int bits;

    public SmallBitSet() { bits = 0;}

    public SmallBitSet copy() {
        final SmallBitSet result = new SmallBitSet();
        result.bits = this.bits;
        return result;
    }

    // construct a singleton set
    public SmallBitSet(final int element) {
        bits = elementIndex(element);
    }

    // construct a two-element set
    public SmallBitSet(final int element1, final int element2) {
        bits = elementIndex(element1) | elementIndex(element2);
    }

    // construct a three-element set
    public SmallBitSet(final int element1, final int element2, final int element3) {
        bits = elementIndex(element1) | elementIndex(element2) | elementIndex(element3);
    }

    public SmallBitSet(final Collection<Integer> elements) {
        bits = 0;
        for (final int element : elements) {
            bits |= elementIndex(element);
        }
    }

    // create a full bit set of all 1s in binary up to a certain number of elements i.e. 00000000000111111....
    public static SmallBitSet fullSet(final int numElements) {
        final SmallBitSet result = new SmallBitSet();
        result.bits = (1 << numElements) - 1;
        return result;
    }

    public SmallBitSet increment() {
        bits++;
        return this;
    }

    public SmallBitSet decrement() {
        bits--;
        return this;
    }

    // the bits as an integer define a unique index within the set of bitsets
    // that is, bitsets can be enumerated as {}, {0}, {1}, {0, 1}, {2}, {0, 2} . . .
    public int index() { return bits; }

    // intersection is equivalent to bitwise AND
    public SmallBitSet intersection(final SmallBitSet other) {
        final SmallBitSet result = new SmallBitSet();
        result.bits = this.bits & other.bits;
        return result;
    }

    // union is equivalent to bitwise OR
    public SmallBitSet union(final SmallBitSet other) {
        final SmallBitSet result = new SmallBitSet();
        result.bits = this.bits | other.bits;
        return result;
    }

    public boolean contains(final SmallBitSet other) {
        return (this.bits & other.bits) == other.bits;
    }

    public void add(final int element) {
        bits |= elementIndex(element);
    }

    public void remove(final int element) {
        bits &= ~(elementIndex(element));
    }

    public void flip(final int element) {
        bits ^= elementIndex(element);
    }

    public boolean get(final int element) {
        return (bits & elementIndex(element)) != 0;
    }

    public boolean isEmpty() { return bits == 0; }

    public boolean hasElementGreaterThan(final int element) { return bits >= 1 << element; }

    private static int elementIndex(final int element) {
        return 1 << element;
    }

    @Override
    public boolean equals(Object o) {
        if (o == this)
            return true;
        if (!(o instanceof SmallBitSet))
            return false;
        SmallBitSet other = (SmallBitSet) o;
        return other.bits == this.bits;
    }

    @Override
    public int hashCode() {
        return bits;
    }

}