package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.InverseAllele;
import org.broadinstitute.hellbender.tools.walkers.mutect.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.haplotype.Event;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.OutputStreamWriter;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Filtering haplotypes that contribute weak alleles to the genotyping. This is a version that determines if allele is weak using
 * Mutect2 model.
 *
 * @author Ilya Soifer &lt;ilya.soifer@ultimagen.com&gt;
 * @author Yossi Farjoun &lt;farjoun@broadinstitute.org&gt;
 *
 */

public class AlleleFilteringMutect extends AlleleFiltering {
    private SomaticGenotypingEngine genotypingEngine;
    private Set<String> normalSamples;

    public AlleleFilteringMutect(M2ArgumentCollection _m2args,
                                 OutputStreamWriter assemblyDebugStream,
                                 SomaticGenotypingEngine _genotypingEngine,
                                 Set<String> _normalSamples,
                                 final SAMFileHeader header){
        super(_m2args, assemblyDebugStream, header);
        genotypingEngine = _genotypingEngine;
        normalSamples = _normalSamples;
    }

    protected double getStringentQuality() { return -50; }

    /**
     * Calculate calculate genotype likelihood of requirement of an allele. Specifically, calculates the likelihood
     * of the data given that allele versus the likelihood of the data when all haplotypes containing the allele are removed
     * This is very similar to what is done in the callMutations function in MutectEngine, but here the haplotypes that do
     * not support the allele are all haplotypes that do not contain the allele rather than only the haplotypes that support reference
     * etc.
     * The log odds is calculated separately for the tumor and the normal samples and the smalles (the strongest) is returned
     *
     * @param alleleLikelihoods
     * @param allele
     * @return likelihood, expressed as phred-scaled confidence
     */

    @Override
    int getAlleleLikelihoodVsInverse(final AlleleLikelihoods<GATKRead, Allele> alleleLikelihoods, Allele allele) {

        final List<LikelihoodMatrix<GATKRead, Allele>> allMatrices = IntStream.range(0, alleleLikelihoods.numberOfSamples())
                .filter(i -> !normalSamples.contains(alleleLikelihoods.getSample(i)))
                .mapToObj(alleleLikelihoods::sampleMatrix)
                .collect(Collectors.toList());
        final AlleleList<Allele> alleleList = allMatrices.get(0);
        final LikelihoodMatrix<GATKRead, Allele> logAllMatrix = SomaticGenotypingEngine.combinedLikelihoodMatrix(allMatrices, alleleList);
        final double alleleLogOddsTumor = somaticAltLogOdds(logAllMatrix);
        logger.debug(() -> String.format("GAL:: Tumor %s: %f", allele.toString(), alleleLogOddsTumor));
        double alleleLogOdds = alleleLogOddsTumor;
        if (!normalSamples.isEmpty()) {
            final List<LikelihoodMatrix<GATKRead, Allele>> allMatricesNormal = IntStream.range(0, alleleLikelihoods.numberOfSamples())
                    .filter(i -> normalSamples.contains(alleleLikelihoods.getSample(i)))
                    .mapToObj(alleleLikelihoods::sampleMatrix)
                    .collect(Collectors.toList());
            final LikelihoodMatrix<GATKRead, Allele> logAllMatrixNormal = SomaticGenotypingEngine.combinedLikelihoodMatrix(allMatricesNormal, alleleList);
            double alleleLogOddsNormal = somaticAltLogOdds(logAllMatrixNormal);
            logger.debug(() -> String.format("GAL:: Normal %s: %f", allele.toString(), alleleLogOddsNormal));
            if (alleleLogOddsNormal < alleleLogOddsTumor) {
                alleleLogOdds = alleleLogOddsNormal;
            }
        }
        return (int)(10*alleleLogOdds);

    }


    /**
     * Calculate the log likelihoods of the ref/alt het genotype for each alt allele, then subtracts
     * these from the hom ref log likelihood to get the log-odds.
     *
     * This method will throw an exception if Allele is an InverseAllele
     *
     * @param matrix a matrix of log likelihoods
     */
    private double somaticAltLogOdds(final LikelihoodMatrix<GATKRead, Allele> matrix) {

        final LikelihoodMatrix<GATKRead, Allele> initialMatrix = matrix;

        if (matrix.getAllele(1-matrix.indexOfReference()) instanceof InverseAllele){
            throw new GATKException.ShouldNeverReachHereException("Inverse allele removed in filtering");
        }
        final LikelihoodMatrix<GATKRead, Allele> logMatrixWithoutThisAllele = SubsettedLikelihoodMatrix.excludingAllele(matrix, matrix.getAllele(1-matrix.indexOfReference()));

        final double logEvidenceWithoutThisAllele = logMatrixWithoutThisAllele.evidenceCount() == 0 ? 0 :
                SomaticLikelihoodsEngine.logEvidence(SomaticGenotypingEngine.getAsRealMatrix(logMatrixWithoutThisAllele),
                        genotypingEngine.makePriorPseudocounts(logMatrixWithoutThisAllele.numberOfAlleles()));
        final double logEvidenceWithAllAlleles= initialMatrix.evidenceCount() == 0 ? 0 :
                SomaticLikelihoodsEngine.logEvidence(SomaticGenotypingEngine.getAsRealMatrix(initialMatrix),
                        genotypingEngine.makePriorPseudocounts(initialMatrix.numberOfAlleles()));
        double tumorLogOdds = (-logEvidenceWithAllAlleles + logEvidenceWithoutThisAllele);
        return tumorLogOdds;
    }
}
