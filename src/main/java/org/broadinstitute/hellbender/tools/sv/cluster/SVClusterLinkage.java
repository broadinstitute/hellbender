package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVLocatable;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public abstract class SVClusterLinkage<T extends SVLocatable> {

    /**
     * Returns whether two given items cluster.
     * @param a first item
     * @param b second item
     */
    public abstract LinkageResult areClusterable(final T a, final T b);

    /**
     * Returns the maximum feasible starting position of any other item with the given item. That is, given item A and
     * X = getMaxClusterableStartingPosition(A), then for any item B on the current contig,
     * Y = start(B) > X => clusterTogether(A, B) == false. Note that this is an upper-bound, but tighter estimates
     * can greatly improve performance.
     * @param item item in question
     * @return max feasible clusterable start coordinate on the current contig
     */
    public abstract int getMaxClusterableStartingPosition(final T item);

    /**
     * Compute max feasible starting position of any other item for all items in the given collection. Note the items
     * must all have the same starting contig.
     */
    public int getMaxClusterableStartingPosition(final Collection<T> items) {
        final List<String> contigA = items.stream().map(T::getContigA).distinct().collect(Collectors.toList());
        if (contigA.size() > 1) {
            throw new IllegalArgumentException("Items start on multiple contigs");
        }
        return items.stream().mapToInt(item -> getMaxClusterableStartingPosition(item)).max().getAsInt();
    }

    /**
     * Returns number of overlapping items
     */
    protected static double getSampleSetOverlap(final Collection<String> a, final Set<String> b) {
        final double denom = Math.max(a.size(), b.size());
        if (denom == 0) {
            return 1;
        }
        return a.stream().filter(b::contains).count() / denom;
    }

    /**
     * Returns true if the overlap is null or exceeds threshold.
     */
    protected static boolean testSampleOverlap(final Double sampleOverlap, final double threshold) {
        return sampleOverlap == null || sampleOverlap >= threshold;
    }


    /**
     * Returns fractional carrier sample overlap in the two records. For CNVs, returns fraction of copy number states that match.
     */
    protected static Double computeSampleOverlap(final SVCallRecord a, final SVCallRecord b) {
        if (a.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV || b.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV) {
            // CNV sample overlap
            final GenotypesContext genotypesA = a.getGenotypes();
            final GenotypesContext genotypesB = b.getGenotypes();
            final Set<String> samples = new HashSet<>(SVUtils.hashMapCapacity(genotypesA.size() + genotypesB.size()));
            samples.addAll(genotypesA.getSampleNames());
            samples.addAll(genotypesB.getSampleNames());
            if (samples.isEmpty()) {
                return null;
            }
            int numMatches = 0;
            for (final String sample : samples) {
                final Genotype genotypeA = genotypesA.get(sample);
                final Genotype genotypeB = genotypesB.get(sample);
                // If one sample doesn't exist in the other set, assume reference copy state
                final int cnA = getCopyState(genotypeA, genotypeB);
                final int cnB = getCopyState(genotypeB, genotypeA);
                if (cnA == cnB) {
                    numMatches++;
                }
            }
            final int numSamples = samples.size();
            return (numMatches / (double) numSamples);
        } else {
            // Non-CNV
            final Set<String> samplesA = a.getCarrierSampleSet();
            final Set<String> samplesB = b.getCarrierSampleSet();
            return getSampleSetOverlap(samplesA, samplesB);
        }
    }

    /**
     * Tries to get the best copy state from the genotype. If the genotype is null, uses ploidy from a "backup"
     * genotype as the default. If we have no clue, just return -1 as a null default.
     */
    private static int getCopyState(final Genotype genotype, final Genotype matchedSampleGenotype) {
        if (genotype == null) {
            if (matchedSampleGenotype != null) {
                return VariantContextGetters.getAttributeAsInt(matchedSampleGenotype, GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, -1);
            } else {
                throw new IllegalArgumentException("Both genotypes are null");
            }
        } else {
            return VariantContextGetters.getAttributeAsInt(genotype, GATKSVVCFConstants.COPY_NUMBER_FORMAT,
                    VariantContextGetters.getAttributeAsInt(genotype, GATKSVVCFConstants.DEPTH_GENOTYPE_COPY_NUMBER_FORMAT, -1));
        }
    }

    /**
     * Used for storing the result of a clustering check between two record and any additional metadata
     */
    public static class LinkageResult {
        private final boolean result;
        public LinkageResult(final boolean result) { this.result = result; }
        public boolean getResult() { return result; }
    }
}
