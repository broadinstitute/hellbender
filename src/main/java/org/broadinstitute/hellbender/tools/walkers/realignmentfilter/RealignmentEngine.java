package org.broadinstitute.hellbender.tools.walkers.realignmentfilter;


import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.spark.sv.utils.Strand;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.*;
import java.util.stream.Collectors;

public class RealignmentEngine {
    private final BwaMemAligner aligner;
    private final int numberOfRegularContigs;

    public RealignmentEngine(final RealignmentArgumentCollection rfac) {
        final BwaMemIndex index = new BwaMemIndex(rfac.bwaMemIndexImage);
        numberOfRegularContigs = rfac.numRegularContigs;
        aligner = new BwaMemAligner(index);
        aligner.setMinSeedLengthOption(rfac.minSeedLength);
        aligner.setDropRatioOption((float) rfac.dropRatio);
        aligner.setSplitFactorOption((float) rfac.splitFactor);
        aligner.setFlagOption(BwaMemAligner.MEM_F_ALL);
    }

    public static boolean supportsVariant(final GATKRead read, final VariantContext vc, int indelStartTolerance) {
        final byte[] readBases = read.getBasesNoCopy();

        final Pair<Integer, CigarOperator> offsetAndOperatorInRead = ReadUtils.getReadIndexForReferenceCoordinate(read, vc.getStart());

        if ( offsetAndOperatorInRead.getLeft() == ReadUtils.READ_INDEX_NOT_FOUND) {
            return false;
        } else if (vc.isSNP() && offsetAndOperatorInRead.getRight() == CigarOperator.DELETION) {
            return false;
        }

        final int variantPositionInRead = offsetAndOperatorInRead.getLeft();

        for (final Allele allele : vc.getAlternateAlleles()) {
            final int referenceLength = vc.getReference().length();
            if (allele.length() == referenceLength) {   // SNP or MNP -- check whether read bases match variant allele bases
                if (allele.basesMatch(ArrayUtils.subarray(readBases, variantPositionInRead, Math.min(variantPositionInRead + allele.length(), readBases.length)))) {
                    return true;
                }
            } else {    // indel -- check if the read has the right CIGAR operator near position -- we don't require exact an exact match because indel representation is non-unique
                final boolean isDeletion = allele.length() < referenceLength;
                int readPosition = 0;    // offset by 1 because the variant position is one base before the first indel base
                for (final CigarElement cigarElement : read.getCigarElements()) {
                    if (Math.abs(readPosition - variantPositionInRead) <= indelStartTolerance) {
                        if (isDeletion && mightSupportDeletion(cigarElement) || !isDeletion && mightSupportInsertion(cigarElement)) {
                            return true;
                        }
                    } else {
                        readPosition += cigarElement.getLength();
                    }
                }
            }
        }
        return false;
    }

    private static boolean mightSupportDeletion(final CigarElement cigarElement) {
        final CigarOperator cigarOperator = cigarElement.getOperator();
        return cigarOperator == CigarOperator.D || cigarOperator == CigarOperator.S;
    }

    private final static boolean mightSupportInsertion(final CigarElement cigarElement) {
        final CigarOperator cigarOperator = cigarElement.getOperator();
        return cigarOperator == CigarOperator.I || cigarOperator == CigarOperator.S;
    }

    public List<BwaMemAlignment> realign(final byte[] sequence) {
        final List<BwaMemAlignment> alignments = aligner.alignSeqs(Arrays.asList(sequence)).get(0).stream()
                .filter(a -> a.getRefId() >= 0)
                .collect(Collectors.toList());

        return alignments.size() == 1 ? alignments :
                alignments.stream().filter(a -> a.getRefId() < numberOfRegularContigs).collect(Collectors.toList());
    }

    /**
     * Find all common alignments of a set of unitigs.  That is, find all alignments in which each unitig aligns within a small distance of the others
     * @param allUnitigAlignments for each unitig, a list of possible alignments
     * @param maxReasonableFragmentLength
     * @return
     */
    public static List<List<BwaMemAlignment>> findJointAlignments(List<List<BwaMemAlignment>> allUnitigAlignments, int maxReasonableFragmentLength) {
        if (allUnitigAlignments.isEmpty()) {
            return Collections.emptyList();
        } else if (allUnitigAlignments.size() == 1) {
            return allUnitigAlignments.get(0).stream().map(Collections::singletonList).collect(Collectors.toList());
        }

        // it's more efficient to start from the unitig with the fewest alignments
        Collections.sort(allUnitigAlignments, Comparator.comparingInt(List::size));

        // make an OverlapDetector for each unitig with padding given by {@code maxReasonableFragmentLength}
        final List<OverlapDetector<BwaMemAlignment>> overlapDetectors = new ArrayList<>();
        for (final List<BwaMemAlignment> unitigAlignments : allUnitigAlignments) {
            final List<SimpleInterval> intervals = unitigAlignments.stream().map(RealignmentEngine::convertToInterval).collect(Collectors.toList());
            final OverlapDetector<BwaMemAlignment> overlapDetector = new OverlapDetector<BwaMemAlignment>(-maxReasonableFragmentLength/2, -maxReasonableFragmentLength/2);
            overlapDetector.addAll(unitigAlignments, intervals);
            overlapDetectors.add(overlapDetector);
        }

        final Set<List<BwaMemAlignment>> commonAlignments = new HashSet<>();
        for (final BwaMemAlignment alignment : allUnitigAlignments.get(0)) {
            final Strand strand = getStrand(alignment);
            final SimpleInterval interval = convertToInterval(alignment);

            // only alignments that have a same-strand overlap every unitig alignment
            if (!overlapDetectors.stream().allMatch(detector -> detector.getOverlaps(interval).stream().anyMatch(a -> getStrand(a) == strand))) {
                continue;
            }

            final List<BwaMemAlignment> alignments = overlapDetectors.stream()
                    .map(detector -> detector.getOverlaps(interval).stream()
                            .filter(a -> getStrand(a) == strand)
                            .max(Comparator.comparingInt(BwaMemAlignment::getAlignerScore)).get())
                    .collect(Collectors.toList());

            commonAlignments.add(alignments);
        }

        return new ArrayList<>(commonAlignments);
    }

    // note the conversion from exclusive end (BWAMemAlignment) to inclusive end (SimpleInterval) (subtract 1 from the end)
    // combined with the conversion from 0-based to 1-based contigs (add 1 to start and end)
    // which has the net result of adding 1 to the start and leaving the end alone
    private static SimpleInterval convertToInterval(final BwaMemAlignment alignment) {
        return new SimpleInterval(Integer.toString(alignment.getRefId()), alignment.getRefStart() + 1, alignment.getRefEnd());
    }

    private static Strand getStrand(final BwaMemAlignment alignment) {
        return SAMFlag.READ_REVERSE_STRAND.isSet(alignment.getSamFlag()) ? Strand.NEGATIVE : Strand.POSITIVE;
    }
}