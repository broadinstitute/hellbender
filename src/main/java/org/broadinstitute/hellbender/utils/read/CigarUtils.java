package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.smithwaterman.SWPairwiseAlignment.Parameters;
import org.broadinstitute.hellbender.utils.smithwaterman.SWPairwiseAlignment;

import java.util.*;

public final class CigarUtils {

    // used in the bubble state machine to apply Smith-Waterman to the bubble sequence
    // these values were chosen via optimization against the NA12878 knowledge base
    public static final Parameters NEW_SW_PARAMETERS = new Parameters(200, -150, -260, -11);

    private static final String SW_PAD = "NNNNNNNNNN";

    private CigarUtils(){}

    /**
     * Combines equal adjacent elements of a Cigar object
     *
     * @param rawCigar the cigar object
     * @return a combined cigar object
     */
    public static Cigar combineAdjacentCigarElements(final Cigar rawCigar) {
        Utils.nonNull(rawCigar);
        final Cigar combinedCigar = new Cigar();
        CigarElement lastElement = null;
        int lastElementLength = 0;
        for (final CigarElement cigarElement : rawCigar.getCigarElements()) {
            if (lastElement != null &&
                    ((lastElement.getOperator() == cigarElement.getOperator()) ||
                            (lastElement.getOperator() == CigarOperator.I && cigarElement.getOperator() == CigarOperator.D) ||
                            (lastElement.getOperator() == CigarOperator.D && cigarElement.getOperator() == CigarOperator.I)))
                lastElementLength += cigarElement.getLength();
            else
            {
                if (lastElement != null)
                    combinedCigar.add(new CigarElement(lastElementLength, lastElement.getOperator()));

                lastElement = cigarElement;
                lastElementLength = cigarElement.getLength();
            }
        }
        if (lastElement != null) {
            combinedCigar.add(new CigarElement(lastElementLength, lastElement.getOperator()));
        }

        return combinedCigar;
    }

    /**
     * Checks whether the cigar has any element that is not H or S
     * @return true the cigar has elements other than S or H, false otherwise.
     */
    public static boolean hasNonClippedBases(final Cigar cigar) {
        return Utils.nonNull(cigar).getCigarElements().stream()
                .anyMatch(el -> el.getOperator() != CigarOperator.SOFT_CLIP && el.getOperator() != CigarOperator.HARD_CLIP);
    }

    /**
     * Inverts the order of the operators in the cigar.
     * Eg 10M1D20M -> 20M1D10M
     */
    public static Cigar invertCigar (final Cigar cigar) {
        Utils.nonNull(cigar);
        final List<CigarElement>  els = new ArrayList<>(cigar.getCigarElements());
        Collections.reverse(els);
        return new Cigar(els);
    }

    /**
     * Compute the number of reference bases between the start (inclusive) and end (exclusive) cigar elements.
     * Note: The method does NOT use CigarOperator.consumesReferenceBases, since it checks something different.
     * The idea is you remove some elements from the beginning of the cigar string,
     * and want to recalculate what if the new starting reference position,
     * you want to count all the elements that indicate existing bases in the reference
     * (everything but I and P).
     * For example original position = 10. cigar: 2M3I2D1M. If you remove the 2M the new starting position is 12.
     * If you remove the 2M3I it is still 12. If you remove 2M3I2D (not reasonable cigar), you will get position 14.
     */
    @SuppressWarnings("fallthru")
    public static int countRefBasesBasedOnCigar(final GATKRead read, final int cigarStartIndex, final int cigarEndIndex){
        if (read == null){
            throw new IllegalArgumentException("null read");
        }
        final List<CigarElement> elems = read.getCigarElements();
        if (cigarStartIndex < 0 || cigarEndIndex > elems.size() || cigarStartIndex > cigarEndIndex){
            throw new IllegalArgumentException("invalid index:" + 0 + " -" + elems.size());
        }
        int result = 0;
        for(int i = cigarStartIndex; i < cigarEndIndex; i++){
            final CigarElement cigarElement = elems.get(i);
            switch (cigarElement.getOperator()) {
                case M:
                case D:
                case N:
                case EQ:
                case X:
                case S:
                case H:
                    result += cigarElement.getLength();
                    break;
                case I:
                case P:        //for these two, nothing happens.
                    break;
                default:
                    throw new GATKException("Unsupported cigar operator: " + cigarElement.getOperator());
            }
        }
        return result;
    }

    /**
     * Removes all clipping operators from the cigar.
     */
    public static Cigar trimReadToUnclippedBases(final Cigar cigar) {
        Utils.nonNull(cigar, "cigar is null");
        final List<CigarElement> elements = new ArrayList<>(cigar.numCigarElements());
        for ( final CigarElement ce : cigar.getCigarElements() ) {
            if ( !isClipOperator(ce.getOperator()) ) {
                elements.add(ce);
            }
        }
        return new Cigar(elements);
    }

    private static boolean isClipOperator(final CigarOperator op) {
        return op == CigarOperator.S || op == CigarOperator.H || op == CigarOperator.P;
    }

    private static boolean isClipOperator(final CigarElement el) {
        return isClipOperator(el.getOperator());
    }

    /**
     * Given a cigar1 and a read with cigar2,
     * this method creates cigar3 such that it has flanking clip operators from cigar2
     * and it has all operators from cigar1 in the middle.
     *
     * In other words if:
     * cigar2 = leftClip2 + noclips2 + rightClip2
     *
     * then
     * cigar3 = leftClip2 + cigar1 + rightClip2
     */
    public static Cigar reclipCigar(final Cigar cigar, final GATKRead read) {
        Utils.nonNull(cigar, "cigar");
        Utils.nonNull(read, "read");

        final List<CigarElement> elements = new ArrayList<>();
        int i = 0;
        final Cigar readCigar = read.getCigar();
        final int n = readCigar.numCigarElements();
        final List<CigarElement> readEls = readCigar.getCigarElements();

        //copy head clips
        while ( i < n && isClipOperator(readEls.get(i)) ) {
            elements.add(readEls.get(i));
            i++;
        }

        elements.addAll(cigar.getCigarElements());

        //skip over non-clips
        i++;
        while ( i < n && !isClipOperator(readEls.get(i)) ) {
            i++;
        }

        //copy tail clips
        while ( i < n && isClipOperator(readEls.get(i)) ) {
            elements.add(readEls.get(i));
            i++;
        }

        return new Cigar(elements);
    }

    /**
     * Returns whether the cigar has any N operators.
     */
    public static boolean containsNOperator(final Cigar cigar) {
        return containsNOperator(Utils.nonNull(cigar).getCigarElements());
    }

    /**
     * Returns whether the list has any N operators.
     */
    public static boolean containsNOperator(final List<CigarElement> cigarElements) {
        return Utils.nonNull(cigarElements).stream().anyMatch(el -> el.getOperator() == CigarOperator.N);
    }

    /**
     * A good Cigar object obeys the following rules:
     *  - is valid as per SAM spec {@link Cigar#isValid(String, long)}.
     *  - has no consecutive I/D elements
     *  - does not start or end with deletions (with or without preceding clips).
     */
    public static boolean isGood(final Cigar c) {
        Utils.nonNull(c, "cigar is null");

        //Note: the string comes from the SAMRecord so it must be a wellformed CIGAR (that is, in "\*|([0-9]+[MIDNSHPX=])+" as per SAM spec).
        //We don't have to check that
        if (c.isValid(null, -1) != null){  //if it's invalid, then it's not good
            return false;
        }
        final List<CigarElement> elems = c.getCigarElements();
        if (hasConsecutiveIndels(elems)){
            return false;
        }
        if (startsWithDeletionIgnoringClips(elems)){
            return false;
        }
        //revert the list and check deletions at the end
        final List<CigarElement> elemsRev = new ArrayList<>(elems);
        Collections.reverse(elemsRev);
        return !startsWithDeletionIgnoringClips(elemsRev);
    }

    /**
     * Checks if cigar has consecutive I/D elements.
     */
    private static boolean hasConsecutiveIndels(final List<CigarElement> elems) {
        boolean prevIndel = false;
        for (final CigarElement elem : elems) {
            final CigarOperator op = elem.getOperator();
            final boolean isIndel = (op == CigarOperator.INSERTION || op == CigarOperator.DELETION);
            if (prevIndel && isIndel) {
                return true;
            }
            prevIndel = isIndel;
        }
        return false;
    }

    /**
     * Checks if cigar starts with a deletion (ignoring any clips at the beginning).
     */
    private static boolean startsWithDeletionIgnoringClips(final List<CigarElement> elems) {
        final Iterator<CigarElement> iter = elems.iterator();
        boolean isClip = true;
        CigarOperator op = null;
        while(iter.hasNext() && isClip) { //consume clips at the beginning
            final CigarElement elem = iter.next();
            op = elem.getOperator();
            isClip = (op == CigarOperator.HARD_CLIP || op == CigarOperator.SOFT_CLIP);
        }
        //once all clips are consumed, is it a deletion or not?
        return op == CigarOperator.DELETION;
    }

    /**
     * Calculate the cigar elements for this path against the reference sequence
     *
     * @param refSeq the reference sequence that all of the bases in this path should align to
     * @return a Cigar mapping this path to refSeq, or null if no reasonable alignment could be found
     */
    public static Cigar calculateCigar(final byte[] refSeq, final byte[] altSeq) {
        Utils.nonNull(refSeq, "refSeq");
        Utils.nonNull(altSeq, "altSeq");
        if ( altSeq.length == 0 ) {
            // horrible edge case from the unit tests, where this path has no bases
            return new Cigar(Arrays.asList(new CigarElement(refSeq.length, CigarOperator.D)));
        }

        final Cigar nonStandard;

        final String paddedRef = SW_PAD + new String(refSeq) + SW_PAD;
        final String paddedPath = SW_PAD + new String(altSeq) + SW_PAD;
        final SWPairwiseAlignment alignment = new SWPairwiseAlignment( paddedRef.getBytes(), paddedPath.getBytes(), NEW_SW_PARAMETERS);

        if ( isSWFailure(alignment) ) {
            return null;
        }


        // cut off the padding bases
        final int baseStart = SW_PAD.length();
        final int baseEnd = paddedPath.length() - SW_PAD.length() - 1; // -1 because it's inclusive
        nonStandard = AlignmentUtils.trimCigarByBases(alignment.getCigar(), baseStart, baseEnd);

        if ( nonStandard.getReferenceLength() != refSeq.length ) {
            nonStandard.add(new CigarElement(refSeq.length - nonStandard.getReferenceLength(), CigarOperator.D));
        }

        // finally, return the cigar with all indels left aligned
        return leftAlignCigarSequentially(nonStandard, refSeq, altSeq, 0, 0);
    }

    /**
     * Make sure that the SW didn't fail in some terrible way, and throw exception if it did
     */
    private static boolean isSWFailure(final SWPairwiseAlignment alignment) {
        // check that the alignment starts at the first base, which it should given the padding
        if ( alignment.getAlignmentStart2wrt1() > 0 ) {
            return true;
        }

        // check that we aren't getting any S operators (which would be very bad downstream)
        for ( final CigarElement ce : alignment.getCigar().getCigarElements() ) {
            if ( ce.getOperator() == CigarOperator.S )
                return true;
            // soft clips at the end of the alignment are really insertions
        }

        return false;
    }

    /**
     * Left align the given cigar sequentially. This is needed because AlignmentUtils doesn't accept cigars with more than one indel in them.
     * This is a target of future work to incorporate and generalize into AlignmentUtils for use by others.
     * @param cigar     the cigar to left align
     * @param refSeq    the reference byte array
     * @param readSeq   the read byte array
     * @param refIndex  0-based alignment start position on ref
     * @param readIndex 0-based alignment start position on read
     * @return          the left-aligned cigar
     */
    public static Cigar leftAlignCigarSequentially(final Cigar cigar, final byte[] refSeq, final byte[] readSeq, int refIndex, int readIndex) {
        Utils.nonNull(cigar, "cigar null");
        Utils.nonNull(refSeq, "refSeq null");
        Utils.nonNull(readSeq, "readSeq null");

        final Cigar cigarToReturn = new Cigar();
        Cigar cigarToAlign = new Cigar();
        for (int i = 0; i < cigar.numCigarElements(); i++) {
            final CigarElement ce = cigar.getCigarElement(i);
            if (ce.getOperator() == CigarOperator.D || ce.getOperator() == CigarOperator.I) {
                cigarToAlign.add(ce);
                final Cigar leftAligned = AlignmentUtils.leftAlignSingleIndel(cigarToAlign, refSeq, readSeq, refIndex, readIndex, false);
                for ( final CigarElement toAdd : leftAligned.getCigarElements() ) { cigarToReturn.add(toAdd); }
                refIndex += cigarToAlign.getReferenceLength();
                readIndex += cigarToAlign.getReadLength();
                cigarToAlign = new Cigar();
            } else {
                cigarToAlign.add(ce);
            }
        }
        if( !cigarToAlign.isEmpty() ) {
            for( final CigarElement toAdd : cigarToAlign.getCigarElements() ) {
                cigarToReturn.add(toAdd);
            }
        }

        final Cigar result = AlignmentUtils.consolidateCigar(cigarToReturn);
        if( result.getReferenceLength() != cigar.getReferenceLength() )
            throw new IllegalStateException("leftAlignCigarSequentially failed to produce a valid CIGAR.  Reference lengths differ.  Initial cigar " + cigar + " left aligned into " + result);
        return result;
    }
}
