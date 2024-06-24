package org.broadinstitute.hellbender.tools.walkers.mutect;


import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.hellbender.tools.walkers.annotator.BaseQuality;
import org.broadinstitute.hellbender.tools.walkers.annotator.ReadPosition;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.locusiterator.AlignmentStateMachine;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.Fragment;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignment;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignmentConstants;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * For each sample and for each allele a list feature vectors of supporting reads
 * In order to reduce the number of delimiter characters, we flatten featurized reads.  For example, suppose allele 1 has featurized reads
 * [1,2] and [3,4] and allele 2 has featurized reads [5,6] and [7,8], the annotation is
 * 1,2,3,4|5,6,7,8
 */
public class FeaturizedReadSets {
    private static final Logger logger = LogManager.getLogger(FeaturizedReadSets.class);

    public static final int DEFAULT_BASE_QUALITY = 25;

    private static final SmithWatermanAligner aligner = SmithWatermanAligner.getAligner(SmithWatermanAligner.Implementation.JAVA);
    private static final int FEATURES_PER_RANGE = 5;
    private static final List<Integer> RANGES = List.of(5, 10, 25, 50);
    public static final int NUM_RANGED_FEATURES = FEATURES_PER_RANGE * RANGES.size();
    private static final int VERY_BAD_QUAL_THRESHOLD = 10;
    private static final int BAD_QUAL_THRESHOLD = 20;

    private FeaturizedReadSets() { }

    public static List<List<List<Integer>>> getReadVectors(final VariantContext vc,
                                                           final Collection<String> samples,
                                                           final AlleleLikelihoods<GATKRead, Allele> likelihoods,
                                                           final AlleleLikelihoods<Fragment, Haplotype> haplotypeLikelihoods,
                                                           final int refDownsample,
                                                           final int altDownsample,
                                                           final M2ArgumentCollection.Mutect3DatasetMode mutect3DatasetMode,
                                                           final Map<String, Integer> readGroupIndices) {
        return getReadVectors(vc, samples, likelihoods, haplotypeLikelihoods, refDownsample, altDownsample, Collections.emptyMap(), mutect3DatasetMode, readGroupIndices);
    }

    // returns Lists (in allele order) of lists of read vectors supporting each allele
    public static List<List<List<Integer>>> getReadVectors(final VariantContext vc,
                                                           final Collection<String> samples,
                                                           final AlleleLikelihoods<GATKRead, Allele> likelihoods,
                                                           final AlleleLikelihoods<Fragment, Haplotype> haplotypeLikelihoods,
                                                           final int refDownsample,
                                                           final int altDownsample,
                                                           final Map<Allele, Integer> altDownsampleMap,
                                                           final M2ArgumentCollection.Mutect3DatasetMode mutect3DatasetMode,
                                                           final Map<String, Integer> readGroupIndices) {
        final Map<Allele, List<GATKRead>> readsByAllele = likelihoods.alleles().stream()
                .collect(Collectors.toMap(a -> a, a -> new ArrayList<>()));

        samples.stream().flatMap(s -> likelihoods.bestAllelesBreakingTies(s).stream())
                .filter(ba -> ba.isInformative())
                .forEach(ba -> readsByAllele.get(ba.allele).add(ba.evidence));

        // downsample if necessary
        final Allele refAllele = likelihoods.alleles().stream().filter(Allele::isReference).findFirst().get();
        for (final Allele allele : likelihoods.alleles()) {
            final int downsample = allele.isReference() ? refDownsample : altDownsampleMap.getOrDefault(allele, altDownsample);
            if (readsByAllele.get(allele).size() > downsample) {
                Collections.shuffle(readsByAllele.get(allele));
                readsByAllele.put(allele, readsByAllele.get(allele).subList(0, downsample));
            }
        }

        final Map<GATKRead, Haplotype> bestHaplotypes = new HashMap<>();
        samples.stream().flatMap(s -> haplotypeLikelihoods.bestAllelesBreakingTies(s).stream())
                .forEach(ba -> ba.evidence.getReads().forEach(read -> bestHaplotypes.put(read, ba.allele)));

        return vc.getAlleles().stream()
                .map(allele -> readsByAllele.get(allele).stream().map(read -> featurize(read, vc, bestHaplotypes, mutect3DatasetMode, readGroupIndices)).collect(Collectors.toList()))
                .collect(Collectors.toList());
    }


    private static List<Integer> featurize(final GATKRead read, final VariantContext vc,
                                           final Map<GATKRead, Haplotype> bestHaplotypes,
                                           final M2ArgumentCollection.Mutect3DatasetMode mutect3DatasetMode,
                                           final Map<String, Integer> readGroupIndices) {
        final List<Integer> result = new ArrayList<>();
        result.add(readGroupIndices.get(read.getReadGroup()));  // this is read group metadata rather than part of the tensor
        result.add(read.getMappingQuality());
        result.add(BaseQuality.getBaseQuality(read, vc).orElse(DEFAULT_BASE_QUALITY));
        result.add(read.isFirstOfPair() ? 1 : 0);
        result.add(read.isReverseStrand() ? 1 : 0);

        // distances from ends of read
        final int readPositionOfVariantStart = ReadPosition.getPosition(read, vc).orElse(0);
        result.add(readPositionOfVariantStart);
        result.add(read.getLength() - readPositionOfVariantStart);


        result.add(Math.abs(read.getFragmentLength()));

        if (mutect3DatasetMode == M2ArgumentCollection.Mutect3DatasetMode.ILLUMINA) {
            // distances from ends of fragment
            final int fragmentStart = Math.min(read.getMateStart(), read.getUnclippedStart());
            final int fragmentEnd = fragmentStart + Math.abs(read.getFragmentLength());
            result.add(vc.getStart() - fragmentStart);
            result.add(fragmentEnd - vc.getEnd());
        }

        // Ultima-specific read tags
        if (mutect3DatasetMode == M2ArgumentCollection.Mutect3DatasetMode.ULTIMA) {
            result.add(read.getAttributeAsInteger("si"));   // si is an integer on the order of 100s or 1000s
            result.add((int) (1000*read.getAttributeAsFloat("rq")));    // rq is a float from 0 and 1, so we multiply by 1000 and round
        }

        // mismatches versus best haplotype
        final Haplotype haplotype = bestHaplotypes.get(read);

        // TODO: fix this
        // I have no idea why this edge case occurs in Ultima data and maybe sometimes in Illumina
        if (!bestHaplotypes.containsKey(read)) {
            logger.warn(String.format("Best haplotypes don't contain read %s at position %s:%d.", read.getName(),
                    vc.getContig(), vc.getStart()));
            result.add(3);
            result.add(2);

            for (int n = 0; n < NUM_RANGED_FEATURES; n++) {
                result.add(0);
            }
        } else {
            byte[] haplotypeBases = haplotype.getBases();
            final SmithWatermanAlignment readToHaplotypeAlignment = aligner.align(haplotypeBases, read.getBases(), SmithWatermanAlignmentConstants.ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS, SWOverhangStrategy.SOFTCLIP);
            final GATKRead copy = read.copy();
            copy.setCigar(readToHaplotypeAlignment.getCigar());
            final int mismatchCount = AlignmentUtils.getMismatchCount(copy, haplotypeBases, readToHaplotypeAlignment.getAlignmentOffset()).numMismatches;
            result.add(mismatchCount);

            final long indelsVsBestHaplotype = readToHaplotypeAlignment.getCigar().getCigarElements().stream().filter(el -> el.getOperator().isIndel()).count();
            result.add((int) indelsVsBestHaplotype);

            final int readStartInHaplotype = readToHaplotypeAlignment.getAlignmentOffset();
            final AlignmentStateMachine asm = new AlignmentStateMachine(copy);
            asm.stepForwardOnGenome();
            final List<int[]> rangedFeatures = RANGES.stream().map(range -> new int[FEATURES_PER_RANGE]).toList();

            while (!asm.isRightEdge()) {
                final PileupElement pe = asm.makePileupElement();
                final int distanceFromVariant = Math.abs(asm.getReadOffset() - readPositionOfVariantStart);

                // pick which array's features we are accounting.  If the ranges are 5, 10, 25, 50 and the distance is, say 8, then the '<= 10' range is relevant
                final OptionalInt relevantRange = IntStream.range(0, RANGES.size()).filter(n -> distanceFromVariant <= RANGES.get(n)).findFirst();
                if (relevantRange.isPresent()) {
                    final int[] featuresToAddTo = rangedFeatures.get(relevantRange.getAsInt());
                    if (pe.isAfterInsertion()) {
                        featuresToAddTo[0]++;
                    }

                    if (pe.isDeletion()) {
                        featuresToAddTo[1]++;
                    } else {
                        final byte base = pe.getBase();
                        final byte qual = pe.getQual();
                        final byte haplotypeBase = haplotypeBases[asm.getGenomeOffset() + readStartInHaplotype];

                        if (base != haplotypeBase) {
                            featuresToAddTo[2]++;
                        }

                        if (qual < VERY_BAD_QUAL_THRESHOLD) {
                            featuresToAddTo[3]++;
                        } else if (qual < BAD_QUAL_THRESHOLD) {
                            featuresToAddTo[4]++;
                        }
                    }
                }
                asm.stepForwardOnGenome();
            }

            for (final int[] featuresToAdd : rangedFeatures) {
                for (final int val : featuresToAdd) {
                    result.add(val);
                }
            }
        }
        // the +1 is for the read group index that comes before the features
        Utils.validate(result.size() == mutect3DatasetMode.getNumReadFeatures() + 1, "Wrong number of features");

        return result;
    }

}
