package org.broadinstitute.hellbender.tools.walkers.contamination;

import htsjdk.samtools.util.OverlapDetector;
import org.apache.commons.lang3.Range;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.tools.copynumber.utils.segmentation.KernelSegmenter;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.*;
import java.util.function.BiFunction;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class ContaminationSegmenter {
    public static final Range<Double> ALT_FRACTIONS_FOR_SEGMENTATION = Range.of(0.1, 0.9);
    public static final double KERNEL_SEGMENTER_LINEAR_COST = 1.0;
    public static final double KERNEL_SEGMENTER_LOG_LINEAR_COST = 1.0;
    public static final int KERNEL_SEGMENTER_DIMENSION = 100;
    public static final int POINTS_PER_SEGMENTATION_WINDOW = 50;
    public static final int MAX_CHANGEPOINTS_PER_CHROMOSOME = 10;
    private static final double SEGMENTATION_KERNEL_VARIANCE = 0.025;

    static final BiFunction<PileupSummary, PileupSummary, Double> SEGMENTATION_KERNEL = (ps1, ps2) -> {
        final double maf1 = FastMath.min(ps1.getAltFraction(), 1 - ps1.getAltFraction());
        final double maf2 = FastMath.min(ps2.getAltFraction(), 1 - ps2.getAltFraction());
        return FastMath.exp(-MathUtils.square(maf1 - maf2)/(2 * SEGMENTATION_KERNEL_VARIANCE));
    };

    private ContaminationSegmenter() {}

    /**
     * Partition the genome into segments of allelic copy number state using kernel segmentation of likely hets.
     * @param sites a list of pileup summaries
     * @return a list of segment intervals.
     */
    public static List<List<PileupSummary>> findSegments(final List<PileupSummary> sites) {
        final Map<String, List<PileupSummary>> sitesByContig = sites.stream().collect(Collectors.groupingBy(PileupSummary::getContig));

        final OverlapDetector<PileupSummary> od = OverlapDetector.create(sites);

        return sitesByContig.values().stream()
                .flatMap(contig -> findContigSegments(contig).stream())
                .map(segment -> od.getOverlaps(segment).stream().sorted(Comparator.comparingInt(PileupSummary::getStart)).collect(Collectors.toList()))
                .collect(Collectors.toList());
    }

    private static List<SimpleInterval> findContigSegments(List<PileupSummary> sites) {
        // segment based on obvious hets
        final List<PileupSummary> hetSites = sites.stream()
                .filter(ps -> ALT_FRACTIONS_FOR_SEGMENTATION.contains(ps.getAltFraction()))
                .collect(Collectors.toList());

        if (hetSites.isEmpty()) {
            return Collections.emptyList();
        }

        final List<Integer> changepoints = new ArrayList<>();
        // when the kernel segmenter finds a changepoint at index n, that means index n belongs to the left segment, which goes
        // against the usual end-exclusive intervals of IndexRange etc.  This explains adding in the first changepoint of -1
        // instead of 0 and all the "changepoint + 1" constructions below
        changepoints.add(-1);
        changepoints.addAll(new KernelSegmenter<>(hetSites).findChangepoints(MAX_CHANGEPOINTS_PER_CHROMOSOME, SEGMENTATION_KERNEL, KERNEL_SEGMENTER_DIMENSION,
                Arrays.asList(POINTS_PER_SEGMENTATION_WINDOW), KERNEL_SEGMENTER_LINEAR_COST, KERNEL_SEGMENTER_LOG_LINEAR_COST, KernelSegmenter.ChangepointSortOrder.INDEX));
        changepoints.add(hetSites.size()-1);

        return IntStream.range(0, changepoints.size() - 1)
                .mapToObj(n -> {
                    final PileupSummary firstSiteInSegment = hetSites.get(changepoints.get(n) + 1);
                    final PileupSummary lastSiteInSegment = hetSites.get(changepoints.get(n+1));
                    return new SimpleInterval(firstSiteInSegment.getContig(), firstSiteInSegment.getStart(), lastSiteInSegment.getEnd());
                }).collect(Collectors.toList());
    }
}
