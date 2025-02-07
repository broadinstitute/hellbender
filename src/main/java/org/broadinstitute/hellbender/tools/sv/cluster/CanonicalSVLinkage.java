package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Iterator;
import java.util.List;

/**
 * <p>Main class for SV clustering. Two items are clustered together if they:</p>
 * <ul>
 *   <li>Match event types</li>
 *   <li>Match strands</li>
 *   <li>Match contigs</li>
 *   <li>Meet minimum interval reciprocal overlap (unless inter-chromosomal or a BND).</li>
 *   <li>Meet minimum sample reciprocal overlap using carrier status. Genotypes (GT fields) are used to determine carrier
 *   status first. If not found and the event is of type DEL or DUP, then copy number (CN fields) are used instead.
 *   In the latter case, it is expected that ploidy can be determined from the number of entries in GT.</li>
 *   <li>Start and end coordinates are within a defined distance.</li>
 * </ul>
 *
 * <p>Interval overlap, sample overlap, and coordinate proximity parameters are defined separately for depth-only/depth-only,
 * depth-only/PESR, and PESR/PESR item pairs using the {@link ClusteringParameters} class. Note that all criteria must be met for
 * two candidate items to cluster, unless the comparison is depth-only/depth-only, in which case only one of
 * interval overlap or break-end proximity must be met. For insertions with undefined length (SVLEN less than 1), interval
 * overlap is tested assuming the given start position and a length of
 * {@value CanonicalSVLinkage#INSERTION_ASSUMED_LENGTH_FOR_OVERLAP}.</p>
 */
public class CanonicalSVLinkage<T extends SVCallRecord> extends SVClusterLinkage<T> {

    protected final SAMSequenceDictionary dictionary;
    protected final boolean clusterDelWithDup;  // DEL/DUP multi-allelic site clustering
    protected ClusteringParameters depthOnlyParams;
    protected ClusteringParameters mixedParams;
    protected ClusteringParameters evidenceParams;

    protected boolean computeLinkageMetrics = false;

    public static final int INSERTION_ASSUMED_LENGTH_FOR_OVERLAP = 50;
    public static final int INSERTION_ASSUMED_LENGTH_FOR_SIZE_SIMILARITY = 1;

    public static final double DEFAULT_RECIPROCAL_OVERLAP_DEPTH_ONLY = 0.8;
    public static final double DEFAULT_SIZE_SIMILARITY_DEPTH_ONLY = 0;
    public static final int DEFAULT_WINDOW_DEPTH_ONLY = 10000000;
    public static final double DEFAULT_SAMPLE_OVERLAP_DEPTH_ONLY = 0;

    public static final double DEFAULT_RECIPROCAL_OVERLAP_MIXED = 0.8;
    public static final double DEFAULT_SIZE_SIMILARITY_MIXED = 0;
    public static final int DEFAULT_WINDOW_MIXED = 1000;
    public static final double DEFAULT_SAMPLE_OVERLAP_MIXED = 0;

    public static final double DEFAULT_RECIPROCAL_OVERLAP_PESR = 0.5;
    public static final double DEFAULT_SIZE_SIMILARITY_PESR = 0;
    public static final int DEFAULT_WINDOW_PESR = 500;
    public static final double DEFAULT_SAMPLE_OVERLAP_PESR = 0;

    protected static final Logger logger = LogManager.getLogger(CanonicalSVLinkage.class);

    public static final ClusteringParameters DEFAULT_DEPTH_ONLY_PARAMS =
            ClusteringParameters.createDepthParameters(
                    DEFAULT_RECIPROCAL_OVERLAP_DEPTH_ONLY,
                    DEFAULT_SIZE_SIMILARITY_DEPTH_ONLY,
                    DEFAULT_WINDOW_DEPTH_ONLY,
                    DEFAULT_SAMPLE_OVERLAP_DEPTH_ONLY);
    public static final ClusteringParameters DEFAULT_MIXED_PARAMS =
            ClusteringParameters.createMixedParameters(
                    DEFAULT_RECIPROCAL_OVERLAP_MIXED,
                    DEFAULT_SIZE_SIMILARITY_MIXED,
                    DEFAULT_WINDOW_MIXED,
                    DEFAULT_SAMPLE_OVERLAP_MIXED);
    public static final ClusteringParameters DEFAULT_PESR_PARAMS =
            ClusteringParameters.createPesrParameters(
                    DEFAULT_RECIPROCAL_OVERLAP_PESR,
                    DEFAULT_SIZE_SIMILARITY_PESR,
                    DEFAULT_WINDOW_PESR,
                    DEFAULT_SAMPLE_OVERLAP_PESR);

    /**
     * Create a new engine
     * @param dictionary sequence dictionary pertaining to clustered records
     * @param clusterDelWithDup enables clustering of DEL/DUP variants into CNVs
     */
    public CanonicalSVLinkage(final SAMSequenceDictionary dictionary,
                              final boolean clusterDelWithDup) {
        Utils.nonNull(dictionary);
        this.dictionary = dictionary;
        this.depthOnlyParams = DEFAULT_DEPTH_ONLY_PARAMS;
        this.mixedParams = DEFAULT_MIXED_PARAMS;
        this.evidenceParams = DEFAULT_PESR_PARAMS;
        this.clusterDelWithDup = clusterDelWithDup;
    }

    public void enableLinkageMetrics() {
        computeLinkageMetrics = true;
    }

    @Override
    public CanonicalLinkageResult areClusterable(final SVCallRecord a, final SVCallRecord b) {
        if (!typesMatch(a, b)) {
            return new CanonicalLinkageResult(false);
        }
        // Only require matching strands if BND or INV type
        if ((a.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.BND
                || a.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.INV)
                && !strandsMatch(a, b)) {
            return new CanonicalLinkageResult(false);
        }
        // Checks appropriate parameter set
        if (evidenceParams.isValidPair(a, b)) {
            return clusterTogetherWithParams(a, b, evidenceParams);
        } else if (mixedParams.isValidPair(a, b)) {
            return clusterTogetherWithParams(a, b, mixedParams);
        } else if (depthOnlyParams.isValidPair(a, b)) {
            return clusterTogetherWithParams(a, b, depthOnlyParams);
        } else {
            return new CanonicalLinkageResult(false);
        }
    }

    /**
     * Tests if SVTYPEs match, allowing for DEL/DUP/CNVs to match in some cases.
     */
    protected boolean typesMatch(final SVCallRecord a, final SVCallRecord b) {
        final GATKSVVCFConstants.StructuralVariantAnnotationType aType = a.getType();
        final GATKSVVCFConstants.StructuralVariantAnnotationType bType = b.getType();
        if (aType == bType && a.getComplexSubtype() == b.getComplexSubtype()) {
            return true;
        }
        // Allow CNVs to cluster with both DELs and DUPs, but only allow DEL/DUP clustering if enabled
        if (a.isSimpleCNV() && b.isSimpleCNV()) {
            if (clusterDelWithDup || (aType == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV || bType == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Test if breakend strands match. True if one has null strands.
     */
    protected boolean strandsMatch(final SVCallRecord a, final SVCallRecord b) {
        if (a.nullStrands() || b.nullStrands()) {
            return true;
        }
        return a.getStrandA() == b.getStrandA() && a.getStrandB() == b.getStrandB();
    }

    /**
     * Helper function to test for clustering between two items given a particular set of parameters.
     * Note that the records have already been subjected to the checks in {@link #areClusterable(SVCallRecord, SVCallRecord)}.
     */
    private static CanonicalLinkageResult clusterTogetherWithParams(final SVCallRecord a, final SVCallRecord b,
                                                                    final ClusteringParameters params) {
        // Contigs match
        if (!(a.getContigA().equals(b.getContigA()) && a.getContigB().equals(b.getContigB()))) {
            return new CanonicalLinkageResult(false);
        }

        // If complex, test complex intervals
        if (a.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.CPX
                && b.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.CPX) {
            return testComplexIntervals(a, b, params.getReciprocalOverlap(), params.getSizeSimilarity(), params.getWindow(), params.getSampleOverlap());
        }

        final Integer breakpointDistance1 = getFirstBreakpointProximity(a, b);
        final Integer breakpointDistance2 = getSecondBreakpointProximity(a, b);
        final Double reciprocalOverlap = computeReciprocalOverlap(a, b);
        final Double sizeSimliarity = computeSizeSimilarity(a, b);

        final boolean hasBreakendProximity = testBreakendProximity(breakpointDistance1, breakpointDistance2, params.getWindow());
        final boolean hasReciprocalOverlap = testReciprocalOverlap(reciprocalOverlap, params.getReciprocalOverlap());
        final boolean hasSizeSimilarity = testSizeSimilarity(sizeSimliarity, params.getSizeSimilarity());
        final boolean passesOverlapAndProximity = params.requiresOverlapAndProximity() ? (hasBreakendProximity && hasReciprocalOverlap) : (hasBreakendProximity || hasReciprocalOverlap);

        // Don't do expensive overlap calculation if it fails other checks or the threshold is 0
        final boolean result;
        if (passesOverlapAndProximity && hasSizeSimilarity && params.getSampleOverlap() > 0) {
            final Double sampleOverlap =  computeSampleOverlap(a, b);
            result = testSampleOverlap(sampleOverlap, params.getSampleOverlap());
        } else {
            result = passesOverlapAndProximity && hasSizeSimilarity;
        }
        return new CanonicalLinkageResult(result, reciprocalOverlap, sizeSimliarity, breakpointDistance1, breakpointDistance2);
    }

    /**
     * Performs overlap testing on each pair of complex intervals in two records, requiring each pair to be
     * sufficiently similar by reciprocal overlap, size similarity, and breakend proximity.
     */
    private static CanonicalLinkageResult testComplexIntervals(final SVCallRecord a, final SVCallRecord b, final double overlapThreshold,
                                                               final double sizeSimilarityThreshold, final int window,
                                                               final double sampleOverlapThreshold) {
        final List<SVCallRecord.ComplexEventInterval> intervalsA = a.getComplexEventIntervals();
        final List<SVCallRecord.ComplexEventInterval> intervalsB = b.getComplexEventIntervals();
        if (intervalsA.size() != intervalsB.size()) {
            return new CanonicalLinkageResult(false);
        }
        final Iterator<SVCallRecord.ComplexEventInterval> iterA = intervalsA.iterator();
        final Iterator<SVCallRecord.ComplexEventInterval> iterB = intervalsB.iterator();
        for (int i = 0; i < intervalsA.size(); i++) {
            final SVCallRecord.ComplexEventInterval cpxIntervalA = iterA.next();
            final SVCallRecord.ComplexEventInterval cpxIintervalB = iterB.next();
            if (cpxIntervalA.getIntervalSVType() != cpxIintervalB.getIntervalSVType()) {
                return new CanonicalLinkageResult(false);
            }
            final Integer breakpointDistance1 = getFirstBreakpointProximity(a, b);
            final Integer breakpointDistance2 = getSecondBreakpointProximity(a, b);
            final SimpleInterval intervalA = cpxIntervalA.getInterval();
            final SimpleInterval intervalB = cpxIintervalB.getInterval();
            final Double reciprocalOverlap = computeReciprocalOverlap(intervalA, intervalB);
            final Double sizeSimilarity = computeSizeSimilarity(intervalA.size(), intervalB.size());
            if (!(testReciprocalOverlap(reciprocalOverlap, overlapThreshold)
                    && testSizeSimilarity(sizeSimilarity, sizeSimilarityThreshold)
                    && testBreakendProximity(breakpointDistance1, breakpointDistance2, window))) {
                return new CanonicalLinkageResult(false);
            }
        }
        // Don't do expensive overlap calculation if threshold is 0
        final Double sampleOverlap = sampleOverlapThreshold > 0 ? computeSampleOverlap(a, b) : Double.valueOf(1.);
        return new CanonicalLinkageResult(testSampleOverlap(sampleOverlap, sampleOverlapThreshold));
    }

    private static Double computeReciprocalOverlap(final SVCallRecord a, final SVCallRecord b) {
        if (a.isIntrachromosomal() && b.isIntrachromosomal() && a.getContigA().equals(b.getContigA())) {
            final int lengthA = getLength(a, INSERTION_ASSUMED_LENGTH_FOR_OVERLAP);
            final int lengthB = getLength(b, INSERTION_ASSUMED_LENGTH_FOR_OVERLAP);
            final SimpleInterval intervalA = new SimpleInterval(a.getContigA(), a.getPositionA(), a.getPositionA() + lengthA - 1);
            final SimpleInterval intervalB = new SimpleInterval(b.getContigA(), b.getPositionA(), b.getPositionA() + lengthB - 1);
            return computeReciprocalOverlap(intervalA, intervalB);
        } else {
            return null;
        }
    }

    private static Double computeReciprocalOverlap(final SimpleInterval a, final SimpleInterval b) {
        if (!a.overlaps(b)) {
            return 0.;
        }
        return a.intersect(b).size() / (double) Math.max(a.size(), b.size());
    }

    private static boolean testReciprocalOverlap(final Double reciprocalOverlap, final double threshold) {
        return reciprocalOverlap == null || reciprocalOverlap >= threshold;
    }

    private static Double computeSizeSimilarity(final SVCallRecord a, final SVCallRecord b) {
        if (a.isIntrachromosomal() && b.isIntrachromosomal()) {
            final int lengthA = getLength(a, INSERTION_ASSUMED_LENGTH_FOR_SIZE_SIMILARITY);
            final int lengthB = getLength(b, INSERTION_ASSUMED_LENGTH_FOR_SIZE_SIMILARITY);
            return computeSizeSimilarity(lengthA, lengthB);
        } else {
            return null;
        }
    }

    private static double computeSizeSimilarity(final int lengthA, final int lengthB) {
        return Math.min(lengthA, lengthB) / (double) Math.max(lengthA, lengthB);
    }

    private static boolean testSizeSimilarity(final Double sizeSimilarity, final double threshold) {
        return sizeSimilarity == null || sizeSimilarity >= threshold;
    }

    private static boolean testBreakendProximity(final Integer distance1, final Integer distance2, final int window) {
        return distance1 != null && distance2 != null && distance1 <= window && distance2 <= window;
    }

    private static Integer getFirstBreakpointProximity(final SVCallRecord a, final SVCallRecord b) {
        if (a.getContigA().equals(b.getContigA())) {
            return Math.abs(a.getPositionA() - b.getPositionA());
        } else {
            return null;
        }
    }

    private static Integer getSecondBreakpointProximity(final SVCallRecord a, final SVCallRecord b) {
        if (a.getContigB().equals(b.getContigB())) {
            return Math.abs(a.getPositionB() - b.getPositionB());
        } else {
            return null;
        }
    }

    /**
     * Gets event length
     */
    private static int getLength(final SVCallRecord record, final int lengthIfMissing) {
        if (record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.BND) {
            if (record.isIntrachromosomal()) {
                // Correct for 0-length case, which is valid for BNDs
                return Math.max(record.getPositionB() - record.getPositionA(), 1);
            } else {
                return 0;
            }
        } else {
            // TODO lengths less than 1 shouldn't be valid
            return Math.max(record.getLength() == null ? lengthIfMissing : record.getLength(), 1);
        }
    }

    @Override
    public int getMaxClusterableStartingPosition(final SVCallRecord record) {
        return Math.max(
                getMaxClusterableStartingPositionWithParams(record, record.isDepthOnly() ? depthOnlyParams : evidenceParams, dictionary),
                getMaxClusterableStartingPositionWithParams(record, mixedParams, dictionary)
        );
    }

    /**
     * Returns max feasible start position of an item clusterable with the given record, given a set of clustering parameters.
     */
    private static int getMaxClusterableStartingPositionWithParams(final SVCallRecord call,
                                                                   final ClusteringParameters params,
                                                                   final SAMSequenceDictionary dictionary) {
        final String contig = call.getContigA();
        final int contigLength = dictionary.getSequence(contig).getSequenceLength();

        // Breakend proximity window
        final int maxPositionByWindow = Math.min(call.getPositionA() + params.getWindow(), contigLength);

        // Don't use overlap for inter-chromosomal events
        if (!call.isIntrachromosomal()) {
            return maxPositionByWindow;
        }

        // Reciprocal overlap window
        final int maxPositionByOverlap;
        final int maxPosition = (int) (call.getPositionA() + (1.0 - params.getReciprocalOverlap()) * getLength(call, INSERTION_ASSUMED_LENGTH_FOR_OVERLAP));
        maxPositionByOverlap = Math.min(maxPosition, contigLength);

        if (params.requiresOverlapAndProximity()) {
            return Math.min(maxPositionByOverlap, maxPositionByWindow);
        } else {
            return Math.max(maxPositionByOverlap, maxPositionByWindow);
        }
    }

    public final ClusteringParameters getDepthOnlyParams() {
        return depthOnlyParams;
    }
    public final ClusteringParameters getMixedParams() {
        return mixedParams;
    }
    public final ClusteringParameters getEvidenceParams() {
        return evidenceParams;
    }

    public final void setDepthOnlyParams(ClusteringParameters depthOnlyParams) { this.depthOnlyParams = depthOnlyParams; }

    public final void setMixedParams(ClusteringParameters mixedParams) {
        this.mixedParams = mixedParams;
    }

    public final void setEvidenceParams(ClusteringParameters evidenceParams) {
        this.evidenceParams = evidenceParams;
    }

    public static class CanonicalLinkageResult extends SVClusterLinkage.LinkageResult {
        private final Double reciprocalOverlap;
        private final Double sizeSimilarity;
        private final Integer breakpointDistance1;
        private final Integer breakpointDistance2;

        public CanonicalLinkageResult(final boolean result) {
            super(result);
            this.reciprocalOverlap = null;
            this.sizeSimilarity = null;
            this.breakpointDistance1 = null;
            this.breakpointDistance2 = null;
        }

        public CanonicalLinkageResult(final boolean result, final Double reciprocalOverlap, final Double sizeSimilarity,
                                      final Integer breakpointDistance1,
                                      final Integer breakpointDistance2) {
            super(result);
            this.reciprocalOverlap = reciprocalOverlap;
            this.sizeSimilarity = sizeSimilarity;
            this.breakpointDistance1 = breakpointDistance1;
            this.breakpointDistance2 = breakpointDistance2;
        }

        public Double getReciprocalOverlap() {
            return reciprocalOverlap;
        }

        public Double getSizeSimilarity() {
            return sizeSimilarity;
        }

        public Integer getBreakpointDistance1() {
            return breakpointDistance1;
        }

        public Integer getBreakpointDistance2() {
            return breakpointDistance2;
        }
    }

}
