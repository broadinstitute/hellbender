package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import javax.ws.rs.core.Variant;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;

public class SomaticGenotypingEngineUnitTest {
    private static final double DEFAULT_AF = 1e-7;

    private static final Allele ATref = Allele.create("AT",true);
    private static final Allele GT = Allele.create("GT",false);

    @DataProvider
    Object[][] germlineAltAlleleFrequenciesData() {
        // {called alleles, germline alleles, germline alt AFs, expected}
        return new Object[][] {
                //biallelic, same alt allele
                { Arrays.asList(Allele.REF_A, Allele.ALT_C), Arrays.asList(Allele.REF_A, Allele.ALT_C), new double[] {0.1}, new double[] {0.1}},

                //biallelic, different alt allele
                { Arrays.asList(Allele.REF_A, Allele.ALT_C), Arrays.asList(Allele.REF_A, Allele.ALT_G), new double[] {0.1}, new double[] {DEFAULT_AF}},

                //triallelic, same alt alleles in same order
                { Arrays.asList(Allele.REF_A, Allele.ALT_C, Allele.ALT_G), Arrays.asList(Allele.REF_A, Allele.ALT_C, Allele.ALT_G), new double[] {0.1, 0.2}, new double[] {0.1, 0.2}},

                //triallelic, same alt alleles in different order
                { Arrays.asList(Allele.REF_A, Allele.ALT_C, Allele.ALT_G), Arrays.asList(Allele.REF_A, Allele.ALT_G, Allele.ALT_C), new double[] {0.1, 0.2}, new double[] {0.2, 0.1}},

                //triallelic, only one alt in common
                { Arrays.asList(Allele.REF_A, Allele.ALT_C, Allele.ALT_G), Arrays.asList(Allele.REF_A, Allele.ALT_G, Allele.ALT_T), new double[] {0.1, 0.2}, new double[] {DEFAULT_AF, 0.1}},

                //biallelic, same alt allele in different representations
                { Arrays.asList(Allele.REF_A, Allele.ALT_G), Arrays.asList(ATref, GT), new double[] {0.1}, new double[] {0.1}},
                { Arrays.asList(ATref, GT), Arrays.asList(Allele.REF_A, Allele.ALT_G), new double[] {0.1}, new double[] {0.1}}
        };
    }

    @Test(dataProvider = "germlineAltAlleleFrequenciesData")
    public void testGetGermlineAltAlleleFrequencies(final List<Allele> calledAlleles, final List<Allele> germlineAlleles, final double[] germlineAltAFs, final double[] expected) {
        final int start = 1;
        final int length = germlineAlleles.stream().mapToInt(Allele::length).max().getAsInt();
        final int stop = start + length - 1;
        final VariantContext vc1 = new VariantContextBuilder("SOURCE", "1", start, stop, germlineAlleles)
                .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, germlineAltAFs).make();
        final double[] result = SomaticGenotypingEngine.getGermlineAltAlleleFrequencies(calledAlleles, List.of(vc1), DEFAULT_AF);
        Assert.assertEquals(result, expected, 1.0e-10);

        // multiallelic -- test splitting into multiple VCs
        if (germlineAlleles.size() > 2) {
            final Allele ref = germlineAlleles.get(0);
            final List<VariantContext> germlineVCs = new ArrayList<>();
            for (int n = 1; n < germlineAlleles.size(); n++) {
                final Allele alt = germlineAlleles.get(n);
                final VariantContext splitVC = new VariantContextBuilder("SOURCE", "1", start, stop, List.of(ref, alt))
                        .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, new double[]{germlineAltAFs[n-1]}).make();
                germlineVCs.add(splitVC);
            }

            final double[] splitResult = SomaticGenotypingEngine.getGermlineAltAlleleFrequencies(calledAlleles, germlineVCs, DEFAULT_AF);
            Assert.assertEquals(splitResult, expected, 1.0e-10);
        }
    }

    @DataProvider(name = "missingAFData")
    Object[][] missingAFData() {
        return new Object[][]{
                {new VariantContextBuilder("SOURCE", "1", 1, 1, Arrays.asList(Allele.REF_A, Allele.ALT_C))
                        .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, VCFConstants.MISSING_VALUE_v4).make()},
                {new VariantContextBuilder("SOURCE", "1", 1, 1, Arrays.asList(Allele.REF_A, Allele.ALT_C))
                        .make()}
        };
    }

    @Test(dataProvider = "missingAFData")
    public void testGetGermlineAltAlleleFrequenciesWithMissingAF(final VariantContext vc) {
        final double[] result = SomaticGenotypingEngine.getGermlineAltAlleleFrequencies(vc.getAlleles(), List.of(vc), DEFAULT_AF);
        Assert.assertEquals(result, new double[] {DEFAULT_AF}, 1.0e-10);
    }

    // test getting alt allele frequencies when each alt has its own VCF line and its own VariantContext
    @Test
    public void testGetGermlineAltAlleleFrequenciesFromSplitAllelesFormat() {
        final double af1 = 0.1;
        final double af2 = 0.2;
        final VariantContext vc1 = new VariantContextBuilder("SOURCE", "1", 1, 1, Arrays.asList(Allele.REF_A, Allele.ALT_C))
                .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, af1).make();
        final VariantContext vc2 = new VariantContextBuilder("SOURCE", "1", 1, 1, Arrays.asList(Allele.REF_A, Allele.ALT_T))
                .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, af2).make();
        final List<Allele> alleles = List.of(Allele.REF_A, Allele.ALT_C, Allele.ALT_T);
        final double[] result = SomaticGenotypingEngine.getGermlineAltAlleleFrequencies(alleles, List.of(vc1, vc2), DEFAULT_AF);
        Assert.assertEquals(result, new double[] {af1, af2}, 1.0e-10);
    }
}