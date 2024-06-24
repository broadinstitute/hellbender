package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.mutable.MutableInt;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.Fragment;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;
import java.util.stream.Collectors;

/**
 *  Count of read pairs in the F1R2 and F2R1 configurations supporting the reference and alternate alleles
 *
 *  <p>This is an annotation that gathers information about the read pair configuration for the reads supporting each
 *  allele. It can be used along with downstream filtering steps to identify and filter out erroneous variants that occur
 *  with higher frequency in one read pair orientation.</p>
 *
 *  <h3>References</h3>
 *  <p>For more details about the mechanism of oxoG artifact generation, see <a href='http://www.ncbi.nlm.nih.gov/pubmed/23303777' target='_blank'>
 *      <i></i>Discovery and characterization of artefactual mutations in deep coverage targeted capture sequencing data due to oxidative DNA damage during sample preparation.</i>
 *  by Costello et al, doi: 10.1093/nar/gks1443</a></p>
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Count of read pairs in the F1R2 and F2R1 configurations supporting REF and ALT alleles (F1R2, F2R1)")
public final class OrientationBiasReadCounts implements JumboGenotypeAnnotation, StandardMutectAnnotation {

    private static final int MINIMUM_BASE_QUALITY = 20;

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(GATKVCFConstants.F1R2_KEY, GATKVCFConstants.F2R1_KEY);
    }

    @Override
    public void annotate(final ReferenceContext ref,
                                  final FeatureContext features,
                                  final VariantContext vc,
                                  final Genotype g,
                                  final GenotypeBuilder gb,
                                  final AlleleLikelihoods<GATKRead, Allele> readLikelihoods,
                                  final AlleleLikelihoods<Fragment, Allele> fragmentLikelihoods,
                                  final AlleleLikelihoods<Fragment, Haplotype> haplotypeLikelihoods){
        Utils.nonNull(gb, "gb is null");
        Utils.nonNull(vc, "vc is null");

        if (g == null || fragmentLikelihoods == null) {
            return;
        }

        final Map<Allele, MutableInt> f1r2Counts = fragmentLikelihoods.alleles().stream()
                .collect(Collectors.toMap(a -> a, a -> new MutableInt(0)));

        final Map<Allele, MutableInt> f2r1Counts = fragmentLikelihoods.alleles().stream()
                .collect(Collectors.toMap(a -> a, a -> new MutableInt(0)));

        // read orientation is a property of fragments; that is, if one read is F1 the other read is R2, hence
        // both reads of the fragment are F1R2.  Thus we can arbitrarily take the first read of each fragment
        Utils.stream(fragmentLikelihoods.bestAllelesBreakingTies(g.getSampleName()))
                .filter(ba -> ba.isInformative() && isUsableRead(ba.evidence.getReads().get(0)) && BaseQualityRankSumTest.getReadBaseQuality(ba.evidence.getReads().get(0), vc).orElse(0) >= MINIMUM_BASE_QUALITY)
                .forEach(ba -> (ReadUtils.isF2R1(ba.evidence.getReads().get(0)) ? f2r1Counts : f1r2Counts).get(ba.allele).increment());

        final int[] f1r2 = vc.getAlleles().stream().mapToInt(a -> f1r2Counts.get(a).intValue()).toArray();

        final int[] f2r1 = vc.getAlleles().stream().mapToInt(a -> f2r1Counts.get(a).intValue()).toArray();

        gb.attribute(GATKVCFConstants.F1R2_KEY, f1r2);
        gb.attribute(GATKVCFConstants.F2R1_KEY, f2r1);
    }

    protected static boolean isUsableRead(final GATKRead read) {
        return read.getMappingQuality() != 0 && read.getMappingQuality() != QualityUtils.MAPPING_QUALITY_UNAVAILABLE;
    }

}