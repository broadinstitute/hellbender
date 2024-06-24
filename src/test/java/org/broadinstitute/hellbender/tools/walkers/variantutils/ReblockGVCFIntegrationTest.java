package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.CommandLineProgramTester;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class ReblockGVCFIntegrationTest extends CommandLineProgramTest {

    private static final String hg38_reference_20_21 = largeFileTestDir + "Homo_sapiens_assembly38.20.21.fasta";
    private static final String b37_reference_20_21 = largeFileTestDir + "human_g1k_v37.20.21.fasta";
    public static final String WARP_PROD_REBLOCKING_ARGS = " -do-qual-approx --floor-blocks -GQB 20 -GQB 30 -GQB 40 ";

    private static final String pf_reference = largeFileTestDir + "PlasmoDB-61_Pfalciparum3D7_Genome.fasta";

    @DataProvider(name = "getCommandLineArgsForExactTest")
    public Object[][] getCommandLineArgsForExactTest() {
        return new Object[][]{
                //covers inputs with old format "MQ" annotation
                {getTestFile("gvcfForReblocking.g.vcf"), getTestFile("testJustOneSample.expected.g.vcf"), " -L chr20:69771 -rgq-threshold 19", hg38_reference_20_21},
                //Broad production arguments on WGS data
                {getTestFile("prodWgsInput.g.vcf "), getTestFile("prodWgsOutput.g.vcf"), WARP_PROD_REBLOCKING_ARGS, hg38Reference},
                //Exome data with AS annotations and zero DP regression test
                {getTestFile("prodWesInput.g.vcf "), getTestFile("prodWesOutput.g.vcf"), WARP_PROD_REBLOCKING_ARGS, hg38Reference}
        };
    }

    @Test(dataProvider = "getCommandLineArgsForExactTest")
    public void testWithExactComparison(final File input, final File expected, final String extraArgs, final String reference) throws Exception {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " -O %s -R " + reference +
                        " -V " + input.getAbsolutePath() +
                        extraArgs +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false",
                Arrays.asList(expected.getAbsolutePath()));
        spec.executeTest("testWithExactComparison", this);
    }

    @Test
    public void testGVCFReblockingIsContiguous() throws Exception {
        final File output = createTempFile("reblockedgvcf", ".vcf");
        final File expected = new File(largeFileTestDir + "testProductionGVCF.expected.g.vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(b37_reference_20_21))
                .add("V", largeFileTestDir + "NA12878.prod.chr20snippet.g.vcf.gz")
                .add("rgq-threshold", "20")
                .add("L", "20:60001-1000000")
                .add("A", "Coverage")
                .add("A", "RMSMappingQuality")
                .add("A", "ReadPosRankSumTest")
                .add("A", "MappingQualityRankSumTest")
                .add("disable-tool-default-annotations", true)
                .addOutput(output);
        runCommandLine(args);

        final CommandLineProgramTester validator = ValidateVariants.class::getSimpleName;
        final ArgumentsBuilder args2 = new ArgumentsBuilder();
        args2.add("R", b37_reference_20_21);
        args2.add("V", output.getAbsolutePath());
        args2.add("L", "20:60001-1000000");
        args2.addRaw("-gvcf");
        validator.runCommandLine(args2);  //will throw a UserException if GVCF isn't contiguous

        try (final FeatureDataSource<VariantContext> actualVcs = new FeatureDataSource<>(output);
             final FeatureDataSource<VariantContext> expectedVcs = new FeatureDataSource<>(expected)) {
            GATKBaseTest.assertCondition(actualVcs, expectedVcs,
                    (a, e) -> VariantContextTestUtils.assertVariantContextsAreEqual(a, e,
                            Collections.emptyList(), Collections.emptyList()));
        }
    }

    @Test
    public void testContiguityWithSpanningDels() {
        final File input = new File(getToolTestDataDir() + "regressionTests.vcf");
        final File intervals = new File(getToolTestDataDir() + "regressionContinuityIntervals.list");
        final File output = createTempFile("reblockedgvcf", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(hg38Reference))
                .add("V", input.getAbsolutePath())
                .addOutput(output);
        runCommandLine(args);

        final CommandLineProgramTester validator = ValidateVariants.class::getSimpleName;
        final ArgumentsBuilder args2 = new ArgumentsBuilder();
        args2.add("R", hg38Reference);
        args2.add("V", output.getAbsolutePath());
        args2.add("L", intervals.getAbsolutePath());
        args2.addRaw("-gvcf");  //TODO: add --no-overlaps after ValidateVariants update is merged
        validator.runCommandLine(args2);  //will throw a UserException if GVCF isn't contiguous
    }

    @Test  //absolute minimal output
    public void testOneSampleAsForGnomAD() throws Exception {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                "-drop-low-quals -do-qual-approx -L chr20:69485-69791 -O %s -R " + hg38_reference_20_21 +
                        " -V " + getToolTestDataDir() + "gvcfForReblocking.g.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false" +
                        " -A Coverage -A RMSMappingQuality -A ReadPosRankSumTest -A MappingQualityRankSumTest --disable-tool-default-annotations true",
                Arrays.asList(getToolTestDataDir() + "testOneSampleAsForGnomAD.expected.g.vcf"));
        spec.executeTest("testOneSampleDropLows", this);
    }

    @Test  //covers non-ref AD and non-ref GT corrections
    public void testNonRefADCorrection() throws Exception {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                "-O %s -R " + hg38_reference_20_21 +
                        " -V " + getToolTestDataDir() + "nonRefAD.g.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false",
                Arrays.asList(getToolTestDataDir() + "testNonRefADCorrection.expected.g.vcf"));
        spec.executeTest("testNonRefADCorrection", this);
    }

    @Test //covers inputs with "RAW_MQ" annotation
    public void testRawMQInput() throws Exception {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                "-O %s -R " + hg38_reference_20_21 +
                        " -V " + getToolTestDataDir() + "prod.chr20snippet.withRawMQ.g.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false",
                Arrays.asList(getToolTestDataDir() + "prod.chr20snippet.withRawMQ.expected.g.vcf"));
        spec.executeTest("testRawMQInput", this);
    }

    @Test
    public void testASAnnotationsAndSubsetting() throws Exception {
        //some subsetting, but never dropping the first alt
        //also has multi-allelic that gets trimmed with ref block added
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                "-O %s -R " + b37_reference_20_21 +
                        " -drop-low-quals -do-qual-approx -V " + "src/test/resources/org/broadinstitute/hellbender/tools/walkers/CombineGVCFs/NA12878.AS.chr20snippet.g.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false",
                Arrays.asList(getToolTestDataDir() + "expected.NA12878.AS.chr20snippet.reblocked.g.vcf"));
        spec.executeTest("testASAnnotationsAndSubsetting", this);

        //one case where first alt is dropped
        final IntegrationTestSpec spec2 = new IntegrationTestSpec(
                "-O %s -R " + b37_reference_20_21 +
                        " -drop-low-quals -do-qual-approx -V " + "src/test/resources/org/broadinstitute/hellbender/tools/walkers/CombineGVCFs/NA12892.AS.chr20snippet.g.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false",
                Arrays.asList(getToolTestDataDir() + "expected.NA12892.AS.chr20snippet.reblocked.g.vcf"));
        spec2.executeTest("testASAnnotationsAndSubsetting2", this);

        //big test for as we ran for gnomADv3
        final File output = createTempFile("reblockedgvcf", ".vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("V", "src/test/resources/large/ReblockGVCF/spanDel.exome.chr20.vcf")
                .add("do-qual-approx", true)
                .add("drop-low-quals", true)
                .add("rgq-threshold", "10")
                .add("L", "chr20")
                .addReference(hg38Reference)
                .addOutput(output);
        runCommandLine(args);

        List<VariantContext> actual = VariantContextTestUtils.getVariantContexts(output);
        actual.stream().forEach(a -> {
            if (!a.getGenotype(0).isHomRef()) {
            VariantContextTestUtils.assertAlleleSpecificAnnotationLengthsCorrect(a, GATKVCFConstants.AS_RAW_QUAL_APPROX_KEY,
                    VCFHeaderLineCount.R);
            VariantContextTestUtils.assertAlleleSpecificAnnotationLengthsCorrect(a, GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY,
                    VCFHeaderLineCount.R);
            VariantContextTestUtils.assertAlleleSpecificAnnotationLengthsCorrect(a, GATKVCFConstants.AS_RAW_MAP_QUAL_RANK_SUM_KEY,
                    VCFHeaderLineCount.R);
            VariantContextTestUtils.assertAlleleSpecificAnnotationLengthsCorrect(a, GATKVCFConstants.AS_RAW_READ_POS_RANK_SUM_KEY,
                    VCFHeaderLineCount.R);
            VariantContextTestUtils.assertAlleleSpecificAnnotationLengthsCorrect(a, GATKVCFConstants.AS_SB_TABLE_KEY,
                    VCFHeaderLineCount.R);
            VariantContextTestUtils.assertAlleleSpecificAnnotationLengthsCorrect(a, GATKVCFConstants.AS_VARIANT_DEPTH_KEY,
                    VCFHeaderLineCount.R);
        } });

    }

    @Test
    public void testNewCompressionScheme() throws Exception {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                "-O %s -R " + b37_reference_20_21 +
                        " -drop-low-quals -do-qual-approx -V " + "src/test/resources/org/broadinstitute/hellbender/tools/walkers/CombineGVCFs/NA12878.AS.chr20snippet.g.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false" +
                        " --floor-blocks -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60",
                Arrays.asList(getToolTestDataDir() + "expected.NA12878.AS.chr20snippet.reblocked.hiRes.g.vcf"));
        spec.executeTest("testNewCompressionScheme", this);
    }

    @Test
    public void testAggressiveQualFiltering() throws Exception {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                "-O %s -R " + hg38_reference_20_21 +
                        " -drop-low-quals -do-qual-approx -V " + getToolTestDataDir() + "gvcfForReblocking.g.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false" +
                        " --floor-blocks" +
                        " --" + GenotypeCalculationArgumentCollection.CALL_CONFIDENCE_LONG_NAME + " 65.0",
                Arrays.asList(getToolTestDataDir() + "expected.aggressiveQualFiltering.g.vcf"));
        spec.executeTest("testVariantQualFiltering", this);
    }

    @Test
    public void testMQHeadersAreUpdated() throws Exception {
        final File output = createTempFile("reblockedgvcf", ".vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("V", getToolTestDataDir() + "justHeader.g.vcf")
                .addReference(hg38Reference)
                .addOutput(output);
        runCommandLine(args);

        Pair<VCFHeader, List<VariantContext>> actual = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        VCFHeader header = actual.getLeft();
        List<VCFInfoHeaderLine> infoLines = new ArrayList<>(header.getInfoHeaderLines());
        //check all the headers in case there's one old and one updated
        for (final VCFInfoHeaderLine line : infoLines) {
            if (line.getID().equals(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED)) {
                Assert.assertTrue(line.getType().equals(VCFHeaderLineType.Float));
                Assert.assertTrue(line.getDescription().contains("deprecated"));
            } else if (line.getID().equals(GATKVCFConstants.MAPPING_QUALITY_DEPTH_DEPRECATED)) {
                Assert.assertTrue(line.getDescription().contains("deprecated"));
            }
        }
    }

    @Test
    public void testReReblocking() {
        final File input = new File(getToolTestDataDir() + "alreadyReblocked.chr22snippet.vcf");
        final File output = createTempFile("rereblockedgvcf", ".vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("V", input)
                .addReference(hg38Reference)
                .addOutput(output);
        runCommandLine(args);

        final List<VariantContext> inputVCs = VariantContextTestUtils.readEntireVCFIntoMemory(input.getAbsolutePath()).getRight();
        final List<VariantContext> outputVCs = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath()).getRight();

        Assert.assertTrue(inputVCs.size() > outputVCs.size());
        Assert.assertTrue(outputVCs.size() == 19);
        //hom ref blocks change, but variants stay the same
        Assert.assertEquals(inputVCs.stream().filter(vc -> !vc.getGenotype(0).isHomRef()).count(), outputVCs.stream().filter(vc -> !vc.getGenotype(0).isHomRef()).count());
        List<String> inGenotypes= inputVCs.stream().filter(vc -> !vc.getGenotype(0).isHomRef()).map(vc -> vc.getGenotype(0)).map(Genotype::toString).collect(Collectors.toList());
        List<String> outGenotypes = outputVCs.stream().filter(vc -> !vc.getGenotype(0).isHomRef()).map(vc -> vc.getGenotype(0)).map(Genotype::toString).collect(Collectors.toList());
        Assert.assertTrue(inGenotypes.containsAll(outGenotypes)); //will check ref and alt alleles as part of genotype string representation
        Assert.assertTrue(outputVCs.get(18).isVariant());

        //all ref blocks have MIN_DP
        Assert.assertEquals(outputVCs.stream().filter(vc -> vc.getGenotype(0).hasExtendedAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY)).count(), outputVCs.size() - outGenotypes.size());
        //all variants have GQ
        Assert.assertEquals(outputVCs.stream().filter(vc -> vc.getGenotype(0).hasGQ()).count(), outputVCs.size());
        //we didn't ask to drop GQ0s, but they might get merged together
        Assert.assertEquals(inputVCs.stream().anyMatch(vc -> vc.getGenotype(0).getGQ() == 0), outputVCs.stream().anyMatch(vc -> vc.getGenotype(0).getGQ() == 0));
    }

    @Test
    public void testOverlappingDeletions() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                "-O %s -R " + hg38_reference_20_21 +
                        " -V " + getToolTestDataDir() + "overlappingDeletions.hc.g.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false",
                Arrays.asList(getToolTestDataDir() + "expected.overlappingDeletions.g.vcf"));
        spec.executeTest("testOverlappingDeletions", this);

        final File output = createTempFile("reblockedgvcf", ".vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("V", getToolTestDataDir() + "overlappingDels2.vcf" )
                .addReference(hg38Reference)
                .addOutput(output);
        runCommandLine(args);

        //regression test for dropped span del allele
        Pair<VCFHeader, List<VariantContext>> actual = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        final List<VariantContext> variants = actual.getRight();
        Assert.assertEquals(variants.size(), 3);
        Assert.assertEquals(variants.get(1).getGenotype(0).getAllele(0), Allele.SPAN_DEL);
    }

    @Test
    public void testHomRefCalls() throws IOException {
        final File input = new File(getToolTestDataDir() + "dropGQ0Dels.g.vcf");
        final File output = createTempFile("dropGQ0Dels.reblocked", ".g.vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("V", input)
                .addReference(hg38Reference)
                .addOutput(output);
        runCommandLine(args);

        final List<VariantContext> inputVCs = VariantContextTestUtils.readEntireVCFIntoMemory(input.getAbsolutePath()).getRight();
        final List<VariantContext> outputVCs = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath()).getRight();

        Assert.assertEquals(outputVCs.size(), 3);
        Assert.assertEquals(outputVCs.get(0).getStart(), inputVCs.get(0).getStart());
        Assert.assertEquals(outputVCs.get(outputVCs.size()-1).getEnd(), inputVCs.get(inputVCs.size()-1).getEnd());
        Assert.assertEquals(outputVCs.get(1).getGenotype(0).getGQ(), 0);  //there should be a GQ0 block in the middle from a crap variant in the input
    }

    @Test
    public void testMultipleInputs() {
        //run with multiple inputs split from chr20:19995000-19998999 of prod.chr20snippet.withRawMQ.g.vcf
        //note that an event is duplicated in shard1 and shard2 because it spans the boundary
        final File output = createTempFile("multi-input", ".vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("V", getToolTestDataDir() + "chr20.shard3.g.vcf")
                .add("V", getToolTestDataDir() + "chr20.shard2.g.vcf")
                .add("V", getToolTestDataDir() + "chr20.shard1.g.vcf")
                .add("V", getToolTestDataDir() + "chr20.shard0.g.vcf")
                .addReference(hg38Reference)
                .addOutput(output);
        runCommandLine(args);

        final File output2 = createTempFile("single-input",".vcf");
        final ArgumentsBuilder args2 = new ArgumentsBuilder();
        args2.add("V", getToolTestDataDir() + "prod.chr20snippet.withRawMQ.g.vcf")
                .add("L", "chr20:19995000-19998999")
                .addReference(hg38Reference)
                .addOutput(output2);
        runCommandLine(args2);

        try (final FeatureDataSource<VariantContext> actualVcs = new FeatureDataSource<>(output);
             final FeatureDataSource<VariantContext> expectedVcs = new FeatureDataSource<>(output2)) {
            GATKBaseTest.assertCondition(actualVcs, expectedVcs,
                    (a, e) -> VariantContextTestUtils.assertVariantContextsAreEqual(a, e,
                            Collections.emptyList(), Collections.emptyList()));
        }
    }

    @Test(expectedExceptions = UserException.class)
    //ReblockGVCF can take multiple inputs, but only if they're non-overlapping shards from the same sample
    public void testMixedSamples() {
        final File output = createTempFile("reblockedgvcf", ".vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("V", getToolTestDataDir() + "justHeader.g.vcf") //sample "Sample"
            .add("V", getToolTestDataDir() + "nonRefAD.g.vcf") //sample "HK017-0046"
            .addReference(hg38Reference)
            .addOutput(output);
        runCommandLine(args);
    }

    @Test
    //we had some external GVCFs that each went through CombineGVCFs for some reason, so GTs all went to ./.
    //chr4 is regression for overlapping ref blocks induced by a hom ref with a real alt
    //also test variants and ref blocks with no DP
    public void testNoCallGenotypes() {
        final File output = createTempFile("reblockedgvcf", ".vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("V", getToolTestDataDir() + "noCallGTs.g.vcf")
                .addReference(hg38Reference)
                .addOutput(output);
        runCommandLine(args);

        Pair<VCFHeader, List<VariantContext>> actual = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        final List<VariantContext> variants = actual.getRight();
        final List<String> variantKeys = variants.stream().map(VariantContextTestUtils::keyForVariant).collect(Collectors.toList());
        final Map<String, VariantContext> resultMap = new LinkedHashMap<>();
        for (int i = 0; i < variants.size(); i++) {
            resultMap.put(variantKeys.get(i), variants.get(i));
        }

        final List<String> expectedHomVarKeys = Arrays.asList(
                "chr22:10514994-10514994 G*, [<NON_REF>, A]",
                "chr22:10515170-10515170 C*, [<NON_REF>, T]",
                "chr22:10515223-10515223 G*, [<NON_REF>, C]",
                "chrY:6067982-6067982 T*, [<NON_REF>, G]");

        final List<String> expectedHetKeys = Arrays.asList(
                "chr4:24339212-24339218 GTATATA*, [<NON_REF>, G]",
                "chr22:10515120-10515120 A*, [<NON_REF>, AAAGC]",
                "chr22:10515223-10515223 G*, [<NON_REF>, C]");

        final List<String> expectedHomRefKeys = Arrays.asList(
                "chr4:24339199-24339209 T*, [<NON_REF>]",
                "chr4:24339210-24339211 G*, [<NON_REF>]",
                "chr4:24339219-24339219 T*, [<NON_REF>]",
                "chr4:24339221-24339221 T*, [<NON_REF>]",
                "chr4:24339222-24339226 A*, [<NON_REF>]",
                "chr22:10515118-10515118 G*, [<NON_REF>]",
                "chrY:6067983-6068157 A*, [<NON_REF>]");

        Assert.assertTrue(variantKeys.containsAll(expectedHomVarKeys));
        Assert.assertTrue(variantKeys.containsAll(expectedHetKeys));
        Assert.assertTrue(variantKeys.containsAll(expectedHomRefKeys));
        Assert.assertTrue(variants.size() == 31);
    }

    @Test
    //here genotype call disagrees with PLs and subsetting can be wrong and lead to no confidence hom-ref if we rely on GT
    public void testPosteriorDisagreementNeedingSubset() {
        final File funkyDragenVariant = new File(getToolTestDataDir() + "HG002.snippet.g.vcf");
        final File output = createTempFile("reblockedgvcf", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("V", funkyDragenVariant)
                .addReference(hg38Reference)
                .addOutput(output);
        runCommandLine(args);

        final Pair<VCFHeader, List<VariantContext>> outVCs = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        Assert.assertEquals(outVCs.getRight().size(), 2);  //one variant and one ref block from trimmed deletion
        final VariantContext testResult = outVCs.getRight().get(0);
        final Genotype g = testResult.getGenotype(0);
        Assert.assertFalse(testResult.hasAttribute(VCFConstants.END_KEY));
        Assert.assertTrue(g.isHetNonRef());
        Assert.assertEquals(testResult.getAlternateAlleles().size(), 3);
        Assert.assertTrue(testResult.getAlternateAlleles().contains(Allele.NON_REF_ALLELE));
        Assert.assertEquals(testResult.getReference().getBaseString().length(), 1);  //alleles are properly trimmed

        final VariantContext newRefBlock = outVCs.getRight().get(1);
        Assert.assertTrue(newRefBlock.hasAttribute(VCFConstants.END_KEY));
        Assert.assertEquals(newRefBlock.getAttributeAsInt(VCFConstants.END_KEY, 0), 103392934);
        final Genotype refG = newRefBlock.getGenotype(0);
        Assert.assertTrue(refG.isHomRef());
        Assert.assertEquals(newRefBlock.getAlternateAlleles().size(), 1);
        Assert.assertTrue(newRefBlock.getAlternateAlleles().contains(Allele.NON_REF_ALLELE));
    }

    @Test
    public void testLeftContigBoundary() {
        final File funkyDragenVariant = new File(getToolTestDataDir() + "testLeftContigBoundary.g.vcf");
        final File output = createTempFile("reblockedgvcf", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("V", funkyDragenVariant)
                .addReference(hg38Reference)
                .addOutput(output);
        runCommandLine(args);

        final List<VariantContext> outVCs = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath()).getRight();
        Assert.assertEquals(outVCs.size(), 2);
        Assert.assertEquals(outVCs.get(0).getStart(), 1);
        Assert.assertEquals(outVCs.get(1).getStart(), 1);
    }

    @Test
    public void testTreeScoreThreshold() {
        final File inputGvcf = new File(getToolTestDataDir() + "treeScoreGvcf.g.vcf");
        final File output = createTempFile("reblockedgvcf", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("V", inputGvcf)
                .addOutput(output)
                .addReference(hg38_reference_20_21)
                .add(ReblockGVCF.TREE_SCORE_THRESHOLD_LONG_NAME, 0.3)
                .add(ReblockGVCF.ANNOTATIONS_TO_KEEP_LONG_NAME, "TREE_SCORE")
                .add("do-qual-approx", true);

        runCommandLine(args);

        final List<VariantContext> outVCs = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath()).getRight();
        for (VariantContext vc : outVCs) {
            if(vc.getStart() == 69511) {
                Assert.assertEquals(vc.getGenotype(0).getGQ(), 99, "Site chr20:69511 should not have been changed from GQ 99.");
            }
            if(vc.getStart() == 69512) {
                Assert.assertEquals(vc.getGenotype(0).getGQ(), 99, "Ref block chr20:69512 should not have been changed from GQ 99.");
            }
            if(vc.getStart() == 69767) {
                Assert.assertEquals(vc.getGenotype(0).getGQ(), 0, "Ref block chr20:69767 should have been grouped with GQ 0 block.");
                Assert.assertEquals(vc.getEnd(), 69783, "Ref block chr20:69767 should have been expanded to include site 69771.");
            }
            Assert.assertNotEquals(vc.getStart(), 69771, "Site chr20:69771 should have been made into a ref block due to low TREE_SCORE.");
            if(vc.getStart() == 69785){
                Assert.assertEquals(vc.getGenotype(0).getGQ(), 0, "Site chr20:69785 should have been made GQ 0 due to missing TREE_SCORE");
            }
            if(!vc.isReferenceBlock()) {
                Assert.assertTrue(vc.hasAttribute(GATKVCFConstants.TREE_SCORE));
            }
        }
    }

    /**
     * Regression test for https://github.com/broadinstitute/gatk/issues/7884 using 2-record snippet of gVCF discussed there.
     */
    @Test
    public void testFirstPositionOnContigNotDropped() {
        final File input = getTestFile("testFirstPositionOnContigNotDropped.g.vcf");
        final File output = createTempFile("reblockedgvcf", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(hg38Reference))
                .add("V", input)
                .add("L", "chr12:18282464") // in the original gVCF, 18282464 is the first variant position on chr12 that is greater than the dropped position 18173860 on chr13
                .add("L", "chr13:18173860")
                .addOutput(output);
        runCommandLine(args);

        // we only check that records at both positions are retained
        final List<VariantContext> outVCs = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath()).getRight();
        Assert.assertEquals(outVCs.get(0).getContig(), "chr12");
        Assert.assertEquals(outVCs.get(0).getStart(), 18282464);
        Assert.assertEquals(outVCs.get(1).getContig(), "chr13");
        Assert.assertEquals(outVCs.get(1).getStart(), 18173860);
    }

    @Test
    public void testTreeScoreWithNoAnnotation() {
        final File output = createTempFile("reblockedgvcf", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("V", largeFileTestDir + "NA12878.prod.chr20snippet.g.vcf.gz")
                .addOutput(output)
                .addReference(b37_reference_20_21)
                .add(ReblockGVCF.TREE_SCORE_THRESHOLD_LONG_NAME, 0.3)
                .add("do-qual-approx", true);

        Assert.assertThrows(UserException.class, () -> runCommandLine(args));
    }

    @Test
    public void testDragenGvcfs() {
        final File input = getTestFile("dragen.g.vcf");
        final File output = createTempFile("reblockedgvcf", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(hg38Reference))
                .add("V", input)
                .add("L", "chr1:939436")
                .addOutput(output);
        runCommandLine(args);

        final VCFHeader outHeader = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath()).getLeft();
        final VariantContext outVC = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath()).getRight().get(0);
        // The extra allele was dropped so this site is now one alt allele and the <NON_REF> allele
        Assert.assertEquals(outVC.getAlternateAlleles().size(), 2);
        final Genotype g = outVC.getGenotype(0);
        Assert.assertEquals(g.getExtendedAttribute(GATKVCFConstants.ALLELE_FRACTION_KEY),"1.00,0.00");
        Assert.assertEquals(g.getExtendedAttribute(GATKVCFConstants.F1R2_KEY), "0,3,0");
        Assert.assertEquals(g.getExtendedAttribute(GATKVCFConstants.F2R1_KEY), "0,3,0");
        for (String attribute : g.getExtendedAttributes().keySet()) {
            final VCFHeaderLineCount countType = outHeader.getFormatHeaderLine(attribute).getCountType();
            if (countType.equals(VCFHeaderLineCount.A)) {
                Assert.assertEquals(((String) g.getExtendedAttribute(attribute)).split(",").length, 2);
            }
            if (countType.equals(VCFHeaderLineCount.R)) {
                Assert.assertEquals(((String) g.getExtendedAttribute(attribute)).split(",").length, 3);
            }
            if (countType.equals(VCFHeaderLineCount.G)) {
                Assert.assertEquals(((String) g.getExtendedAttribute(attribute)).split(",").length, 6);
            }
        }
    }

    @Test
    public void testFilters() throws IOException {
        final File input = getTestFile("dragen.g.vcf");
        final File output = createTempFile("reblockedgvcf", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(hg38Reference))
                .add("V", input)
                .add("L", "chr1")
                .add(ReblockGVCF.KEEP_SITE_FILTERS_LONG_NAME, true)
                .addOutput(output);
        runCommandLine(args);

        final VariantContext filteredVC = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath()).getRight().get(3); // last site in the file
        Assert.assertTrue(filteredVC.isFiltered());

        final VariantContext unfilteredVC = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath()).getRight().get(1);
        Assert.assertFalse(unfilteredVC.isFiltered());

        final VariantContext filteredRefBlockVC = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath()).getRight().get(0);
        Assert.assertFalse(filteredRefBlockVC.isFiltered()); // Ref block is unfiltered even though the input RefBlock and low qual variant were both filtered
        Assert.assertEquals(filteredRefBlockVC.getGenotype(0).getDP(), 12); // Ref block is combination of filtered variant with depth 22 and filtered ref block with depth 1
    }

    @Test
    public void testMovingFilters() throws IOException {
        final File input = getTestFile("dragen.g.vcf");
        final File output = createTempFile("reblockedgvcf", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(hg38Reference))
                .add("V", input)
                .add("L", "chr1")
                .add(ReblockGVCF.ADD_FILTERS_TO_GENOTYPE, true)
                .addOutput(output);
        runCommandLine(args);

        final VariantContext filteredVC = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath()).getRight().get(3); // last site in the file
        Assert.assertFalse(filteredVC.isFiltered());
        Assert.assertTrue(filteredVC.getGenotype(0).isFiltered());

        final VariantContext unfilteredVC = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath()).getRight().get(1);
        Assert.assertFalse(unfilteredVC.isFiltered());
        Assert.assertFalse(unfilteredVC.getGenotype(0).isFiltered());

        final VariantContext filteredRefBlockVC = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath()).getRight().get(0);
        Assert.assertFalse(filteredRefBlockVC.isFiltered()); // Ref block is unfiltered even though the input RefBlock and low qual variant were both filtered
        Assert.assertFalse(filteredRefBlockVC.getGenotype(0).isFiltered()); // Ref block genotype is also unfiltered
    }

    @Test
    public void testRemovingFormatAnnotations() {
        final File input = getTestFile("dragen.g.vcf");
        final File output = createTempFile("reblockedgvcf", ".vcf");
        final String priKey = "PRI";

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(hg38Reference))
                .add("V", input)
                .add(ReblockGVCF.ANNOTATIONS_TO_REMOVE_LONG_NAME, priKey)
                .addOutput(output);
        runCommandLine(args);

        final List<VariantContext> outVCs = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath()).getRight();
        for(VariantContext vc : outVCs){
            Assert.assertNull(vc.getGenotype(0).getExtendedAttribute(priKey));
        }
    }

    @Test
    public void testNonGVCFInput() {
        final File output = createTempFile("reblockedgvcf", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(b37_reference_20_21))
                .add("V", "src/test/resources/large/NA12878.HiSeq.WGS.b37_decoy.indel.recalibrated.chr20.vcf")
                .addOutput(output);

        Assert.assertThrows(GATKException.class, () -> runCommandLine(args));
    }

    @Test
    public void testAddedHcAsAnnotations() throws IOException {

        final File input = new File(largeFileTestDir + "reblockGVCFs_AS_test.p_falciparum.g.vcf");
        final File output = createTempFile("reblockedgvcf", ".vcf");
        final File expected = new File(largeFileTestDir + "expected.reblockGVCFs_AS_test.p_falciparum.rb.g.vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();

        args.addReference(new File(pf_reference))
                .add("V", input)
                .addFlag("do-qual-approx")
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false")
                .add("A", "AssemblyComplexity")
                .addFlag("floor-blocks")
                .add("GQB", 20)
                .add("GQB", 30)
                .add("GQB", 40)
                .addOutput(output);

        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(
                output.toPath(),
                expected.toPath(),
                "#",
                true
        );
    }
    @Test
    public void testSpanDelRegression() throws IOException {

        final File input = getTestFile("reblock_cleanup_bug_variant.g.vcf");
        final File output = createTempFile("reblockedgvcf", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();

        args.addReference(new File(pf_reference))
                .add("V", input)
                .addFlag("do-qual-approx")
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false")
                .add("A", "AssemblyComplexity")
                .addFlag("floor-blocks")
                .add("GQB", 20)
                .add("GQB", 30)
                .add("GQB", 40)
                .addOutput(output);

        runCommandLine(args);

        Pair<VCFHeader, List<VariantContext>> actual = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        List<VariantContext> variants = actual.getRight();
        Assert.assertEquals(variants.size(), 4);
        VariantContext v0 = variants.get(0);
        //crappy deletions are all collapsed into GQ0 ref block
        Assert.assertEquals(v0.getStart(), 646914);
        Assert.assertEquals(v0.getAttributeAsInt(VCFConstants.END_KEY, 0), 646953);
        Assert.assertTrue(v0.getGenotype(0).isHomRef());
        Assert.assertEquals(v0.getGenotype(0).getGQ(), 0);
        //no more star alleles because deletions are all gone
        Assert.assertFalse(variants.stream().anyMatch(v -> v.getAlternateAlleles().contains(Allele.SPAN_DEL)));
    }

    @Test
    public void testRemovingAnnotationFromEmptyAttributes() {
        final File input = getTestFile("gvcfWithNoPRI.vcf");
        final File output = createTempFile("reblockedgvcf", ".vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder();

        args.addReference(hg38Reference)
                .addVCF(input)
                .addOutput(output)
                .add("format-annotations-to-remove","PRI");

        //make sure it doesn't error
        runCommandLine(args);
    }
}
