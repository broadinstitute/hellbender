package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.util.List;
import java.util.Set;

import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.COPY_NUMBER_FORMAT;

/**
 * Clustering engine class for defragmenting depth-based DEL/DUP calls, such as those produced by
 * {@link org.broadinstitute.hellbender.tools.copynumber.GermlineCNVCaller GermlineCNVCaller}. Each variant is first padded by a fraction
 * of its length and then merged with any other overlapping variants that meet the minimum sample overlap. Additionally,
 * singleton variants (with only one carrier sample) will only be linked with variants of the same copy number.
 *
 * This class does not return any metadata via {@link org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVLinkage.CanonicalLinkageResult}.
 */
public class CNVLinkage extends SVClusterLinkage<SVCallRecord> {

    public static final double DEFAULT_SAMPLE_OVERLAP = 0.8;
    public static final double DEFAULT_PADDING_FRACTION = 0.25;

    protected final double minSampleOverlap;
    protected final double paddingFraction;
    protected final SAMSequenceDictionary dictionary;

    public CNVLinkage(final SAMSequenceDictionary dictionary, final double paddingFraction,
                      final double minSampleOverlap) {
        super();
        this.dictionary = dictionary;
        this.minSampleOverlap = minSampleOverlap;
        this.paddingFraction = paddingFraction;
    }

    @Override
    public LinkageResult areClusterable(final SVCallRecord a, final SVCallRecord b) {
        // Only do clustering on depth-only CNVs
        if (!a.isDepthOnly() || !b.isDepthOnly()) return new LinkageResult(false);
        if (!a.isSimpleCNV() || !b.isSimpleCNV()) return new LinkageResult(false);
        Utils.validate(a.getContigA().equals(a.getContigB()), "Variant A is a CNV but interchromosomal");
        Utils.validate(b.getContigA().equals(b.getContigB()), "Variant B is a CNV but interchromosomal");

        // Types match
        if (a.getType() != b.getType()) return new LinkageResult(false);

        // Interval overlap
        // Positions should be validated already by the SVCallRecord class - these checks are for thoroughness
        final SimpleInterval intervalA = getPaddedRecordInterval(a.getContigA(), a.getPositionA(), a.getPositionB());
        Utils.nonNull(intervalA, "Invalid interval " + new SimpleInterval(a.getContigA(), a.getPositionA(),
                a.getPositionB()) + " for record " + a.getId());
        final SimpleInterval intervalB = getPaddedRecordInterval(b.getContigA(), b.getPositionA(), b.getPositionB());
        Utils.nonNull(intervalB, "Invalid interval " + new SimpleInterval(b.getContigA(), b.getPositionA(),
                b.getPositionB()) + " for record " + b.getId());
        if (!intervalA.overlaps(intervalB)) return new LinkageResult(false);

        // Sample overlap
        if (minSampleOverlap > 0) {
            final double sampleOverlap = computeSampleOverlap(a, b);
            if (!testSampleOverlap(sampleOverlap, minSampleOverlap)) {
                return new LinkageResult(false);
            }
        }

        // In the single-sample case, match copy number strictly if we're looking at the same sample
        // TODO repeated check for CN attributes in hasSampleOverlap and getCarrierSamples
        final Set<String> carriersA = a.getCarrierSampleSet();
        final Set<String> carriersB = b.getCarrierSampleSet();
        if (carriersA.size() == 1 && carriersA.equals(carriersB)) {
            final Genotype genotypeA = a.getGenotypes().get(carriersA.iterator().next());
            final Genotype genotypeB = b.getGenotypes().get(carriersB.iterator().next());
            if (genotypeA.hasExtendedAttribute(COPY_NUMBER_FORMAT) && genotypeB.hasExtendedAttribute(COPY_NUMBER_FORMAT)) {
                final int copyNumberA = VariantContextGetters.getAttributeAsInt(genotypeA, COPY_NUMBER_FORMAT, 0);
                final int copyNumberB = VariantContextGetters.getAttributeAsInt(genotypeB, COPY_NUMBER_FORMAT, 0);
                final int copyNumberDeltaA = genotypeA.getPloidy() - copyNumberA;
                final int copyNumberDeltaB = genotypeB.getPloidy() - copyNumberB;
                if (copyNumberDeltaA != copyNumberDeltaB) {
                    return new LinkageResult(false);
                }
            } else {
                final List<Allele> sortedAllelesA = SVCallRecordUtils.sortAlleles(genotypeA.getAlleles());
                final List<Allele> sortedAllelesB = SVCallRecordUtils.sortAlleles(genotypeB.getAlleles());
                if (!sortedAllelesA.equals(sortedAllelesB)) {
                    return new LinkageResult(false);
                }
            }
        }
        return new LinkageResult(true);
    }

    /**
     * Max clusterable position depends on the other variant's length, which is unknown, so we assume the largest
     * possible event that would extend to the end of the contig.
     */
    @Override
    public int getMaxClusterableStartingPosition(final SVCallRecord call) {
        final int contigLength = dictionary.getSequence(call.getContigA()).getSequenceLength();
        if (!call.isSimpleCNV()) {
            return 0;
        }
        final int maxTheoreticalStart = (int) Math.floor((call.getPositionB() + paddingFraction * (call.getLength() + contigLength)) / (1.0 + paddingFraction));
        return Math.min(maxTheoreticalStart, dictionary.getSequence(call.getContigA()).getSequenceLength());
    }

    /**
     * Determine an overlap interval for clustering using specified padding.
     * @return padded interval
     */
    protected SimpleInterval getPaddedRecordInterval(final String contig, final int start, final int end) {
        return new SimpleInterval(contig, start, end)
                .expandWithinContig((int) (paddingFraction * (end - start + 1)), dictionary);
    }
}
