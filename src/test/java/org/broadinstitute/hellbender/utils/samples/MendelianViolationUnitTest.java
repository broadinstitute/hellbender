package org.broadinstitute.hellbender.utils.samples;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.*;

import java.util.*;

public final class MendelianViolationUnitTest extends GATKBaseTest {

    private static final int GQ30 = 30;
    private static final int GQ25 = 25;

    private static final Allele refAllele = Allele.create("A", true);
    private static final Allele altAllele = Allele.create("T");
    private static final Allele noCallAllele = Allele.NO_CALL;

    private static final List<Allele> homRef = Arrays.asList(refAllele, refAllele);
    private static final List<Allele> homVar = Arrays.asList(altAllele, altAllele);
    private static final List<Allele> het = Arrays.asList(refAllele, altAllele);
    private static final List<Allele> notCall = Arrays.asList(noCallAllele, noCallAllele);

    // haploid alleles
    private static final List<Allele> hapRef = Collections.singletonList(refAllele);
    private static final List<Allele> hapVar = Collections.singletonList(altAllele);
    private static final List<Allele> hapNotCall = Collections.singletonList(noCallAllele);

    private static final Genotype g00 = new GenotypeBuilder("2", homRef).DP(10).AD(new int[]{10,0}).GQ(GQ30).make();
    private static final Genotype g01 = new GenotypeBuilder("1", het).DP(10).AD(new int[]{5, 5}).GQ(GQ25).make();
    private static final Genotype g11 = new GenotypeBuilder("3", homVar).DP(10).AD(new int[]{0, 10}).GQ(GQ30).make();
    private static final Genotype gNo = new GenotypeBuilder("4", notCall).make();

    //haploid genotypes
    private static final Genotype ghap0  = new GenotypeBuilder("5", hapRef).DP(10).AD(new int[]{10}).GQ(GQ30).make();
    private static final Genotype ghap1  = new GenotypeBuilder("6", hapVar).DP(10).AD(new int[]{5}).GQ(GQ30).make();
    private static final Genotype ghapNo = new GenotypeBuilder("7", hapNotCall).make();

    private static final Sample sMom = new Sample("mom", "fam", null, null, Sex.FEMALE);
    private static final Sample sDad = new Sample("dad",     sMom.getFamilyID(), null, null, Sex.MALE);
    private static final Sample sChild = new Sample("child", sMom.getFamilyID(), sDad.getID(), sMom.getID(), Sex.MALE);

    @DataProvider(name = "isViolation")
    public Object[][] isViolation() {
        final List<Object[]> tests = new ArrayList<>();

        tests.add(new Object[]{g00, g00, g00, false});
        tests.add(new Object[]{g00, g00, g01, true});
        tests.add(new Object[]{g00, g00, g11, true});
        tests.add(new Object[]{g00, g00, gNo, false});

        tests.add(new Object[]{g00, g01, g00, false});
        tests.add(new Object[]{g00, g01, g01, false});
        tests.add(new Object[]{g00, g01, g11, true});
        tests.add(new Object[]{g00, g01, gNo, false});

        tests.add(new Object[]{g00, g11, g00, true});
        tests.add(new Object[]{g00, g11, g01, false});
        tests.add(new Object[]{g00, g11, g11, true});
        tests.add(new Object[]{g00, g11, gNo, false});

        tests.add(new Object[]{g00, gNo, g00, false});
        tests.add(new Object[]{g00, gNo, g01, false});
        tests.add(new Object[]{g00, gNo, g11, true});
        tests.add(new Object[]{g00, gNo, gNo, false});


        tests.add(new Object[]{g01, g00, g00, false});
        tests.add(new Object[]{g01, g00, g01, false});
        tests.add(new Object[]{g01, g00, g11, true});
        tests.add(new Object[]{g01, g00, gNo, false});

        tests.add(new Object[]{g01, g01, g00, false});
        tests.add(new Object[]{g01, g01, g01, false});
        tests.add(new Object[]{g01, g01, g11, false});
        tests.add(new Object[]{g01, g01, gNo, false});

        tests.add(new Object[]{g01, g11, g00, true});
        tests.add(new Object[]{g01, g11, g01, false});
        tests.add(new Object[]{g01, g11, g11, false});
        tests.add(new Object[]{g01, g11, gNo, false});

        tests.add(new Object[]{g01, gNo, g00, false});
        tests.add(new Object[]{g01, gNo, g01, false});
        tests.add(new Object[]{g01, gNo, g11, false});
        tests.add(new Object[]{g01, gNo, gNo, false});


        tests.add(new Object[]{g11, g00, g00, true});
        tests.add(new Object[]{g11, g00, g01, false});
        tests.add(new Object[]{g11, g00, g11, true});
        tests.add(new Object[]{g11, g00, gNo, false});

        tests.add(new Object[]{g11, g01, g00, true});
        tests.add(new Object[]{g11, g01, g01, false});
        tests.add(new Object[]{g11, g01, g11, false});
        tests.add(new Object[]{g11, g01, gNo, false});

        tests.add(new Object[]{g11, g11, g00, true});
        tests.add(new Object[]{g11, g11, g01, true});
        tests.add(new Object[]{g11, g11, g11, false});
        tests.add(new Object[]{g11, g11, gNo, false});

        tests.add(new Object[]{g11, gNo, g00, true});
        tests.add(new Object[]{g11, gNo, g01, false});
        tests.add(new Object[]{g11, gNo, g11, false});
        tests.add(new Object[]{g11, gNo, gNo, false});


        tests.add(new Object[]{gNo, g00, g00, false});
        tests.add(new Object[]{gNo, g00, g01, false});
        tests.add(new Object[]{gNo, g00, g11, true});
        tests.add(new Object[]{gNo, g00, gNo, false});

        tests.add(new Object[]{gNo, g01, g00, false});
        tests.add(new Object[]{gNo, g01, g01, false});
        tests.add(new Object[]{gNo, g01, g11, false});
        tests.add(new Object[]{gNo, g01, gNo, false});

        tests.add(new Object[]{gNo, g11, g00, true});
        tests.add(new Object[]{gNo, g11, g01, false});
        tests.add(new Object[]{gNo, g11, g11, false});
        tests.add(new Object[]{gNo, g11, gNo, false});

        tests.add(new Object[]{gNo, gNo, g00, false});
        tests.add(new Object[]{gNo, gNo, g01, false});
        tests.add(new Object[]{gNo, gNo, g11, false});
        tests.add(new Object[]{gNo, gNo, gNo, false});
      
        tests.add(new Object[]{ghap0, ghap0, ghap0, false});
        tests.add(new Object[]{ghap0, ghap0, ghap1, true});
        tests.add(new Object[]{g00, g00, ghap1, true});
        tests.add(new Object[]{g01, g00, ghap1, false});
        tests.add(new Object[]{g00, g00, ghap0, false});
        tests.add(new Object[]{ghap0, ghap0, ghapNo, false});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "isViolation")
    public void testIsViolation(final Genotype gMom, final Genotype gDad, final Genotype gChild, final boolean expected) throws Exception {
        final boolean actual = MendelianViolation.isViolation(gMom, gDad, gChild);
        Assert.assertEquals(expected, actual);
    }

    @DataProvider(name = "isViolationWithSample")
    private Iterator<Object[]> isViolationWithSample() {
        List<Object[]> ret = new ArrayList<>();
        for (Object[] o : isViolation()) {
            ret.add(getIsViolationWithSampleTestCase((Genotype)o[0], (Genotype)o[1], (Genotype)o[2], (boolean)o[3], 0));
            ret.add(getIsViolationWithSampleTestCase((Genotype)o[0], (Genotype)o[1], (Genotype)o[2], (boolean)o[3], GQ25+1));
            ret.add(getIsViolationWithSampleTestCase((Genotype)o[0], (Genotype)o[1], (Genotype)o[2], false, GQ30+1));
        }

        return ret.iterator();
    }

    private Object[] getIsViolationWithSampleTestCase(final Genotype gMom, final Genotype gDad, final Genotype gChild, final boolean originalExpected, final int minGenotypeQuality) {
        final boolean isLowQual = minGenotypeQuality > 0 && (gMom.getGQ() < minGenotypeQuality || gDad.getGQ() < minGenotypeQuality || gChild.getGQ() < minGenotypeQuality);
        final boolean expectedWithQual = !isLowQual && originalExpected;
        return new Object[]{gMom, gDad, gChild, minGenotypeQuality, expectedWithQual, isLowQual};
    }

    private boolean isNoCall(Genotype gMom, Genotype gDad, Genotype gChild) {
        return (!gMom.isCalled() && !gDad.isCalled() || !gChild.isCalled());
    }

    private VariantContext getVC(final Genotype gMom, final Genotype gDad, final Genotype gChild) {
        final ArrayList<Genotype> genotypes = new ArrayList<>(Arrays.<Genotype>asList(gMom, gDad, gChild));
        final Map<String, Integer> sampleNameToOffset= new LinkedHashMap<>(3);
        sampleNameToOffset.put("mom", 0);
        sampleNameToOffset.put("dad", 1);
        sampleNameToOffset.put("child", 2);
        final List<String> sampleNamesInOrder=Arrays.asList("mom", "dad", "child");
        final GenotypesContext context = GenotypesContext.create(genotypes, sampleNameToOffset, sampleNamesInOrder);

        final Collection<Allele> alleles= new LinkedHashSet<>(2);
        alleles.addAll(gMom.getAlleles());
        alleles.addAll(gDad.getAlleles());
        alleles.addAll(gChild.getAlleles());
        alleles.add(refAllele);
        alleles.remove(noCallAllele);

        return new VariantContextBuilder("test", "20", 10, 10, alleles).genotypes(context).make();
    }

    @Test(dataProvider = "isViolationWithSample")
    public void testParentsRefRefChildHet(final Genotype gMom, final Genotype gDad, final Genotype gChild, final int minGenotypeQuality, final boolean expected, final boolean isLowQual) throws Exception {
        final MendelianViolation mv = new MendelianViolation(minGenotypeQuality, false, false);
        final boolean actual = mv.isViolation(sMom, sDad, sChild, getVC(gMom, gDad, gChild));
        Assert.assertEquals(actual, expected);
        Assert.assertEquals(mv.getViolationsCount() > 0, expected);

        Assert.assertEquals(mv.getParentsRefRefChildHet() > 0, !isLowQual && (gMom.isHomRef() && gDad.isHomRef() && gChild.isHet()));
    }

    @Test(dataProvider = "isViolationWithSample")
    public void testFamilyNoCallCount(final Genotype gMom, final Genotype gDad, final Genotype gChild, final int minGenotypeQuality, final boolean expected, final boolean isLowQual) throws Exception {
        final MendelianViolation mv = new MendelianViolation(minGenotypeQuality, false, false);
        final boolean actual = mv.isViolation(sMom, sDad, sChild, getVC(gMom, gDad, gChild));
        Assert.assertEquals(actual, expected);

        Assert.assertEquals(mv.getFamilyCalledCount() > 0, !isLowQual && !isNoCall(gMom, gDad, gChild));
    }

    @Test(dataProvider = "isViolationWithSample")
    public void testFamilyLowQualsCount(final Genotype gMom, final Genotype gDad, final Genotype gChild, final int minGenotypeQuality, final boolean expected, final boolean isLowQual) throws Exception {
        final MendelianViolation mv = new MendelianViolation(minGenotypeQuality, false, false);
        final boolean actual = mv.isViolation(sMom, sDad, sChild, getVC(gMom, gDad, gChild));
        Assert.assertEquals(actual, expected);

        Assert.assertEquals(mv.getFamilyLowQualsCount() > 0, !isNoCall(gMom, gDad, gChild) && isLowQual);
    }

    @Test(dataProvider = "isViolationWithSample")
    public void testVarFamilyCalledCount(final Genotype gMom, final Genotype gDad, final Genotype gChild, final int minGenotypeQuality, final boolean expected, final boolean isLowQual) throws Exception {
        final MendelianViolation mv = new MendelianViolation(minGenotypeQuality, false, false);
        final boolean actual = mv.isViolation(sMom, sDad, sChild, getVC(gMom, gDad, gChild));
        Assert.assertEquals(actual, expected);

        Assert.assertEquals(mv.getVarFamilyCalledCount() == 0, (isLowQual || isNoCall(gMom, gDad, gChild) || (gMom.isHomRef() && gDad.isHomRef() && gChild.isHomRef())));
    }
}
