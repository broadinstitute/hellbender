 package org.broadinstitute.hellbender.utils;


import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.Serializable;

 /**
 * Minimal immutable class representing a 1-based closed ended genomic interval
 * SimpleInterval does not allow null contig names.  It cannot represent an unmapped Locatable.
 *
 *@warning 0 length intervals are NOT currently allowed, but support may be added in the future
 */
public final class SimpleInterval implements Locatable, Serializable {
    private static final long serialVersionUID = 1L;
    public static final char CONTIG_SEPARATOR = ':';
    public static final char START_END_SEPARATOR = '-';
    public static final String END_OF_CONTIG = "+"; //note: needs to be a String because it's used in an endsWith call.

    private final int start;
    private final int end;
    private final String contig;

    /**
     * Create a new immutable 1-based interval of the form [start, end]
     * @param contig the name of the contig, must not be null
     * @param start  1-based inclusive start position
     * @param end  1-based inclusive end position
     */
    public SimpleInterval(final String contig, final int start, final int end){
        validatePositions(contig, start, end);
        this.contig = contig;
        this.start = start;
        this.end = end;
    }

    /**
     * Create a new SimpleInterval from a {@link Locatable}
     * @param locatable any Locatable
     * @throws IllegalArgumentException if locatable violates any of the SimpleInterval constraints or is null
     */
    public SimpleInterval(final Locatable locatable){
        this(Utils.nonNull(locatable).getContig(),
                locatable.getStart(), locatable.getEnd());
    }

    /**
     * Test that these are valid values for constructing a SimpleInterval:
     *    contig cannot be null
     *    start must be >= 1
     *    end must be >= start
     * @throws IllegalArgumentException if it is invalid
     */
    static void validatePositions(final String contig, final int start, final int end) {
        Utils.validateArg(isValid(contig, start, end), () -> "Invalid interval. Contig:" + contig + " start:"+start + " end:" + end);
    }

     /**
      * Test that these are valid values for constructing a SimpleInterval:
      *    contig cannot be null
      *    start must be >= 1
      *    end must be >= start
      */
     public static boolean isValid(final String contig, final int start, final int end) {
         return contig != null && start > 0 && end >= start;
     }

     /**
     * Makes an interval by parsing the string.
     *
     * @warning this method does not fill in the true contig end values
     * for intervals that reach to the end of their contig,
     * uses {@link Integer#MAX_VALUE} instead.
     *
     * Semantics of start and end are defined in {@link Locatable}.
     * The format is one of:
     *
     * contig           (Represents the whole contig, from position 1 to the {@link Integer#MAX_VALUE})
     *
     * contig:start     (Represents the 1-element range start-start on the given contig)
     *
     * contig:start-end (Represents the range start-end on the given contig)
     *
     * contig:start+    (Represents the prefix of the contig starting at the given start position and ending at {@link Integer#MAX_VALUE})
     *
     * examples (note that _all_ commas in numbers are simply ignored, for human convenience):
     *
     * 'chr2', 'chr2:1000000' or 'chr2:1,000,000-2,000,000' or 'chr2:1000000+'
      *
      * @param str non-empty string to be parsed
     */
    public SimpleInterval(final String str){
        /* Note: we want to keep the class immutable. So all fields need to be final.
         * But only constructors can assign to final fields.
         * So we can either keep this parsing code in the constructor or make a static factory method
         * and make multiple objects. We chose the former.
         */
        Utils.nonNull(str);
        Utils.validateArg(!str.isEmpty(), "str should not be empty");

        final String contig;
        final int start;
        final int end;

        final int colonIndex = str.lastIndexOf(CONTIG_SEPARATOR);
        if (colonIndex == -1) {
            contig = str;  // chr1
            start = 1;
            end = Integer.MAX_VALUE;
        } else {
            contig = str.substring(0, colonIndex);
            final int dashIndex = str.indexOf(START_END_SEPARATOR, colonIndex);
            if(dashIndex == -1) {
                if(str.endsWith(END_OF_CONTIG)) {
                    start = parsePosition(str.substring(colonIndex + 1, str.length() - 1));  // chr:1+
                    end = Integer.MAX_VALUE;
                } else {
                    start = parsePosition(str.substring(colonIndex + 1));   // chr1:1
                    end = start;
                }
            } else {
                start = parsePosition(str.substring(colonIndex + 1, dashIndex));  // chr1:1-1
                end = parsePosition(str.substring(dashIndex + 1));
            }
        }

        validatePositions(contig, start, end);
        this.contig = contig;
        this.start = start;
        this.end = end;
    }

     /**
      * Parses a number like 100000 or 1,000,000 into an int. Throws NumberFormatException on parse failure.
      */
     static int parsePositionThrowOnFailure(final String pos) throws NumberFormatException {
         return Integer.parseInt(pos.replaceAll(",", "")); //strip commas
     }

     /**
     * Parses a number like 100000 or 1,000,000 into an int.
     */
    static int parsePosition(final String pos) {
        try {
            return parsePositionThrowOnFailure(pos);
        } catch (NumberFormatException e){
            throw new UserException("Problem parsing start/end value in interval string. Value was: " + pos, e);
        }
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final SimpleInterval that = (SimpleInterval) o;

        if (end != that.end) return false;
        if (start != that.start) return false;
        return contig.equals(that.contig);
    }

    @Override
    public int hashCode() {
        int result = start;
        result = 31 * result + end;
        result = 31 * result + contig.hashCode();
        return result;
    }

    @Override
    public String toString() {
        return IntervalUtils.locatableToString(this);
    }

     /**
     * @return name of the contig this is mapped to
     */
    @Override
    public String getContig(){
        return contig;
    }

    /** Gets the 1-based start position of the interval on the contig. */
    @Override
    public int getStart(){
        return start;
    }

    /**
    * @return the 0-based start position (from the GA4GH spec).
    */
    public long getGA4GHStart() {return start - 1; }

    /**
     * @return the 1-based closed-ended end position of the interval on the contig.
     */
    @Override
    public int getEnd(){
        return end;
    }

    /**
    * @return the typical end spans are [zero-start,end) (from the GA4GH spec).
    */
    public long getGA4GHEnd() { return end; }

    /**
     * @return number of bases covered by this interval (will always be > 0)
     */
    public int size() {
        return end - start + 1;
    }

    /**
     * Determines whether this interval overlaps the provided locatable.
     *
     * @param other interval to check
     * @return true if this interval overlaps other, otherwise false
     */
    public boolean overlaps( final Locatable other ) {
        return overlapsWithMargin(other, 0);
    }

     /**
      * Determines whether this interval comes within "margin" of overlapping the provided locatable.
      * This is the same as plain overlaps if margin=0.
      *
      * @param other interval to check
      * @param margin how many bases may be between the two interval for us to still consider them overlapping; must be non-negative
      * @return true if this interval overlaps other, otherwise false
      * @throws IllegalArgumentException if margin is negative
      */
     public boolean overlapsWithMargin(final Locatable other, final int margin) {
         if ( margin < 0 ) {
             throw new IllegalArgumentException("given margin is negative: " + margin +
                     "\tfor this: " + toString() + "\tand that: " + (other == null ? "other is null" : other.toString()));
         }
         if ( other == null || other.getContig() == null ) {
             return false;
         }

         return this.contig.equals(other.getContig()) && this.start <= other.getEnd() + margin && other.getStart() - margin <= this.end;
     }


     /**
     * Determines whether this interval contains the entire region represented by other
     * (in other words, whether it covers it).
     *
     * @param other interval to check
     * @return true if this interval contains all of the bases spanned by other, otherwise false
     */
    public boolean contains( final Locatable other ) {
        if ( other == null || other.getContig() == null ) {
            return false;
        }

        return this.contig.equals(other.getContig()) && this.start <= other.getStart() && this.end >= other.getEnd();
    }

     /**
      * Returns the intersection of the two intervals. The intervals must overlap or IllegalArgumentException will be thrown.
      */
     public SimpleInterval intersect( final Locatable that ) {
         Utils.validateArg(this.overlaps(that), () ->
                 "SimpleInterval::intersect(): The two intervals need to overlap " + this + " " + that);

         return new SimpleInterval(getContig(),
                 Math.max(getStart(), that.getStart()),
                 Math.min( getEnd(), that.getEnd()) );
     }

     /**
      * Get section of starting interval (this) that is not overlapped by the other interval (that)
      * @param that - interval to subtract from starting interval. Must overlap (but not fully contain) starting interval
      * @return - SimpleInterval representing the portion of starting interval (this) not overlapped by other interval (that)
      */
     @VisibleForTesting
     public SimpleInterval subtract(final Locatable that) {
         Utils.validateArg(this.overlaps(that), () ->
                 "SimpleIntervaL::subtract(): The two intervals need to overlap: " + this + ", " + that);
         Utils.validateArg(!that.contains(this), () ->
                 "SimpleIntervaL::subtract(): Interval to subtract " + that + " cannot contain starting interval " + this);
         if (this.getStart() < that.getStart()) {
             return new SimpleInterval(this.getContig(), this.getStart(), that.getStart());
         }
         else {
             return new SimpleInterval(this.getContig(), that.getEnd(), this.getEnd());
         }
     }

     /**
      * Returns a new SimpleInterval that represents the entire span of this and that.  Requires that
      * this and that SimpleInterval are contiguous.
      */
     public SimpleInterval mergeWithContiguous( final Locatable that ) {
         Utils.nonNull(that);
         if (!this.contiguous(that)) {
             throw new GATKException("The two intervals need to be contiguous: " + this + " " + that);
         }

         return new SimpleInterval(getContig(),
                 Math.min( getStart(), that.getStart() ),
                 Math.max( getEnd(), that.getEnd()) );
     }

     /**
      * Returns a new SimpleInterval that represents the region between the endpoints of this and other.
      *
      * Unlike {@link #mergeWithContiguous}, the two intervals do not need to be contiguous
      *
      * @param other the other interval with which to calculate the span
      * @return a new SimpleInterval that represents the region between the endpoints of this and other.
      */
     public SimpleInterval spanWith( final Locatable other ) {
         Utils.nonNull(other);
         Utils.validateArg(this.getContig().equals(other.getContig()), "Cannot get span for intervals on different contigs");
         return new SimpleInterval(contig, Math.min(start, other.getStart()), Math.max(end, other.getEnd()));
     }

     /**
      * Determines whether this interval is contiguous with (ie., overlaps or is directly adjacent to) the provided interval
      *
      * @param that interval to compare with
      * @return true if this is contiguous with that, otherwise false
      */
     public boolean contiguous(final Locatable that) {
         Utils.nonNull(that);
         return this.getContig().equals(that.getContig()) && this.getStart() <= that.getEnd() + 1 && that.getStart() <= this.getEnd() + 1;
     }

     /**
      * Returns a new SimpleInterval that represents this interval as expanded by the specified amount in both
      * directions, bounded by the contig start/stop if necessary.
      *
      * @param padding amount to expand this interval
      * @param contigLength length of this interval's contig
      * @return a new SimpleInterval that represents this interval as expanded by the specified amount in both
      *         directions, bounded by the contig start/stop if necessary.
      */
     public SimpleInterval expandWithinContig( final int padding, final int contigLength ) {
         Utils.validateArg(padding >= 0, "padding must be >= 0");
         return IntervalUtils.trimIntervalToContig(contig, start - padding, end + padding, contigLength);
     }

     /**
      * Returns a new SimpleInterval that represents this interval as expanded by the specified amount in both
      * directions, bounded by the contig start/stop if necessary.
      *
      * @param padding amount to expand this interval
      * @param sequenceDictionary dictionary to use to determine the length of this interval's contig
      * @return a new SimpleInterval that represents this interval as expanded by the specified amount in both
      *         directions, bounded by the contig start/stop if necessary.
      */
     public SimpleInterval expandWithinContig( final int padding, final SAMSequenceDictionary sequenceDictionary ) {
         Utils.nonNull(sequenceDictionary);
         final SAMSequenceRecord contigRecord = sequenceDictionary.getSequence(contig);
         Utils.nonNull( contigRecord, () -> "Contig " + contig + " not found in provided dictionary");

         return expandWithinContig(padding, contigRecord.getSequenceLength());
     }
 }
