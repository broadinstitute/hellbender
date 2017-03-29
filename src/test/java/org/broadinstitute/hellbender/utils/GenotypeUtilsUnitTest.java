package org.broadinstitute.hellbender.utils;

import htsjdk.variant.variantcontext.*;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GenotypeUtilsUnitTest {
    private static final Allele Aref = Allele.create("A", true);
    private static final Allele T = Allele.create("T");



    @DataProvider
    public Object[][] getGTsWithAndWithoutLikelihoods(){
        final GenotypeBuilder builder = new GenotypeBuilder("sample", Arrays.asList(Aref, T));
        final List<Allele> DIPLOID_ALLELES = Arrays.asList(Aref, T);
        final List<Allele> HAPLOID_ALLELES = Arrays.asList(Aref);
        final int[] SOME_PLS = {0, 4, 1};
        final int[] HAPLOID_PL = {0, 4};
        return new Object[][]{
                {builder.noPL().alleles(DIPLOID_ALLELES).make(), false},
                {builder.noPL().alleles(HAPLOID_ALLELES).make(), false},
                {builder.PL(HAPLOID_PL).alleles(HAPLOID_ALLELES).make(), false},
                {builder.PL(SOME_PLS).alleles(DIPLOID_ALLELES).make(), true},
                {GenotypeBuilder.createMissing("sample", 2), false}
        };
    }

    @Test(dataProvider = "getGTsWithAndWithoutLikelihoods")
    public void testIsDiploidWithLikelihoods(Genotype g, boolean expected) throws Exception {
        Assert.assertEquals(GenotypeUtils.isDiploidWithLikelihoods(g), expected);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIsDiploidWithLikelihoodsWithNull() {
        GenotypeUtils.isDiploidWithLikelihoods(null);
    }



    @DataProvider
    public Object[][] getGenotypeCountsParameters(){
        final VariantContext vc = new VariantContextBuilder("in memory", "1", 100, 100,
                                                                        Arrays.asList(Aref, T)).make();
        final ArrayList<Genotype> genotypesArray = new ArrayList<>();
        for(int i = 0; i<100; i++){
            final Genotype g = new GenotypeBuilder("sample" + i, Arrays.asList(Aref, T)).PL(new int[]{100,10,0}).make();
            genotypesArray.add(g);
        }
        final GenotypesContext genotypes = GenotypesContext.create(genotypesArray);

        return new Object[][]{
                {vc, genotypes, false, new GenotypeCounts(0, 9, 91)},
                {vc, genotypes, true, new GenotypeCounts(0, 0, 100)},
        };
    }

    @Test(dataProvider = "getGenotypeCountsParameters")
    public void testRounding(VariantContext vc, GenotypesContext gt, boolean round, GenotypeCounts expected) {
        final GenotypeCounts actual = GenotypeUtils.computeDiploidGenotypeCounts(vc, gt, round);
        Assert.assertEquals(actual.getRefs(), expected.getRefs());
        Assert.assertEquals(actual.getHets(), expected.getHets());
        Assert.assertEquals(actual.getHoms(), expected.getHoms());
    }


}