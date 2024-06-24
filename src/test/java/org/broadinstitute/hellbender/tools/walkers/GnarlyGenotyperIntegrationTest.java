package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.GenomicsDBTestUtils;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static org.broadinstitute.hellbender.testutils.VariantContextTestUtils.getVariantContexts;

/**
 * Created by gauthier on 3/9/18.
 */
public class GnarlyGenotyperIntegrationTest extends CommandLineProgramTest {

    // If true, update the expected outputs in tests that assert an exact match vs. prior output,
    // instead of actually running the tests. Can be used with "./gradlew test -Dtest.single=GnarlyGenotyperIntegrationTest"
    // to update all of the exact-match tests at once. After you do this, you should look at the
    // diffs in the new expected outputs in git to confirm that they are consistent with expectations.
    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;

    private static final List<String> NO_EXTRA_ARGS = Collections.emptyList();

    /*
     * Make sure that someone didn't leave the UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS toggle turned on
     */
    @Test
    public void assertThatExpectedOutputUpdateToggleIsDisabled() {
        Assert.assertFalse(UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS, "The toggle to update expected outputs should not be left enabled");
    }

    @DataProvider(name="VCFdata")
    public Object[][] getVCFdata() {
        return new Object[][]{
               /*
                //chrX haploid sample plus diploid sample -- expected results validated with vcf-validator (samtools?)
                {new File[]{getTestFile("NA12891.chrX.haploid.rb.g.vcf"), getTestFile("NA12892.chrX.diploid.rb.g.vcf")},
                        getTestFile("haploidPlusDiploid.expected.vcf"), null, Arrays.asList(new SimpleInterval("chrX", 1000000, 5000000)), Arrays.asList("--merge-input-intervals", "--only-output-calls-starting-in-intervals"), b38_reference_20_21},

                //8 ALT alleles -- no PLs
                {new File[]{getTestFile("sample6.vcf"), getTestFile("sample7.vcf"), getTestFile("sample8.vcf"), getTestFile("sample9.vcf")},
                        getTestFile("lotsOfAltsNoPLs.vcf"), null, Arrays.asList(new SimpleInterval("chr20", 257008, 257008)), Arrays.asList("--merge-input-intervals", "--only-output-calls-starting-in-intervals"), b38_reference_20_21},
                //6 ALT alleles -- yes PLs
                {new File[]{getTestFile("sample6.vcf"), getTestFile("sample7.vcf"), getTestFile("sample8.vcf")},
                        getTestFile("lotsOfAltsYesPLs.vcf"), null, Arrays.asList(new SimpleInterval("chr20", 257008, 257008)), Arrays.asList("--merge-input-intervals", "--only-output-calls-starting-in-intervals"), b38_reference_20_21},
        */
        // Simple Test, spanning deletions; standard calling confidence
                //No variants outside requested intervals; no SNPs with QUAL < 60, no INDELs with QUAL < 69?; has star alleles after deletion at chr20:263497; has AC, AF, AN, DP, ExcessHet, FS, MQ, (MQRankSum), (ReadPosRankSum), SOR, QD; has called genotypes
                {new File[]{getTestFile("sample1.vcf"), getTestFile("sample2.vcf"), getTestFile("sample3.vcf"), getTestFile("sample4.vcf"), getTestFile("sample5.vcf")},
                         getTestFile("fiveSampleTest.vcf"), null, Arrays.asList(new SimpleInterval("chr20", 251370, 252000), new SimpleInterval("chr20", 263000, 265600)), Arrays.asList("--merge-input-intervals", "--only-output-calls-starting-in-intervals"), b38_reference_20_21},
                // Same as above, but with GenomicsDB using BCF2Codec for interchange
                {new File[]{getTestFile("sample1.vcf"), getTestFile("sample2.vcf"), getTestFile("sample3.vcf"), getTestFile("sample4.vcf"), getTestFile("sample5.vcf")},
                         getTestFile("fiveSampleTest.vcf"), null, Arrays.asList(new SimpleInterval("chr20", 251370, 252000), new SimpleInterval("chr20", 263000, 265600)), Arrays.asList("--merge-input-intervals", "--only-output-calls-starting-in-intervals", "--genomicsdb-use-bcf-codec"), b38_reference_20_21},
                //lower calling confidence
                //same as above except (different intervals and) with SNPs with 40 < QUAL < 60 and INDELs with 49 < QUAL < 69
                {new File[]{getTestFile("sample1.vcf"), getTestFile("sample2.vcf"), getTestFile("sample3.vcf"), getTestFile("sample4.vcf"), getTestFile("sample5.vcf")},
                         getTestFile("fiveSampleTest.lowerCallThreshold.vcf"), null, Arrays.asList(new SimpleInterval("chr20", 250865, 348163)), Arrays.asList("-stand-call-conf 10"), b38_reference_20_21},
                // Same as above, but with GenomicsDB using BCF2Codec for interchange
                {new File[]{getTestFile("sample1.vcf"), getTestFile("sample2.vcf"), getTestFile("sample3.vcf"), getTestFile("sample4.vcf"), getTestFile("sample5.vcf")},
                         getTestFile("fiveSampleTest.lowerCallThreshold.vcf"), null, Arrays.asList(new SimpleInterval("chr20", 250865, 348163)), Arrays.asList("-stand-call-conf 10", "--genomicsdb-use-bcf-codec"), b38_reference_20_21},
                //using latest reblocking output with allele-specific annotations
                //has all of the above plus AS_AltDP, AS_FS, AS_MQ, AS_MQRankSum, AS_QD, AS_ReadPosRankSum
                {new File[]{new File(getToolTestDataDir() + "/../variantutils/ReblockGVCF/expected.NA12878.AS.chr20snippet.reblocked.g.vcf"),
                         new File(getToolTestDataDir() + "/../variantutils/ReblockGVCF/expected.NA12892.AS.chr20snippet.reblocked.g.vcf")},
                         getTestFile("twoSampleAS.vcf"), getTestFile("twoSampleASDB.vcf"), Arrays.asList(new SimpleInterval("20")), NO_EXTRA_ARGS, b37_reference_20_21},
                // Same as above, but with GenomicsDB using BCF2Codec for interchange
                {new File[]{new File(getToolTestDataDir() + "/../variantutils/ReblockGVCF/expected.NA12878.AS.chr20snippet.reblocked.g.vcf"),
                        new File(getToolTestDataDir() + "/../variantutils/ReblockGVCF/expected.NA12892.AS.chr20snippet.reblocked.g.vcf")},
                        getTestFile("twoSampleAS.vcf"), getTestFile("twoSampleASDB.vcf"), Arrays.asList(new SimpleInterval("20")), Arrays.asList("--genomicsdb-use-bcf-codec"), b37_reference_20_21},
                //using legacy reblocking data with no raw GT count values
                //has ExcessHet, calculated from genotypes
                {new File[]{new File(getToolTestDataDir() + "noGTCount.sample1.chr20snippet.vcf"),
                        new File(getToolTestDataDir() + "noGTCount.sample2.chr20snippet.vcf")},
                        getTestFile("noGTCount.expected.chr20snippet.vcf"), null, Arrays.asList(new SimpleInterval("chr20")), NO_EXTRA_ARGS, hg38Reference
                }


        };
    }

    /**
     * Have a customized checker as the VariantContextUtils.assertEqualVariants version does a
     * simplistic ToString() comparison that does not account for differences in Scientific versus Decimal
     * notations.
     * @param actualVC List of actual variant contexts
     * @param expectedVC List of expected variant contexts
     */
    private static void assertEqualVariants(
            final List<VariantContext> actualVC, final List<VariantContext> expectedVC) {
        Utils.nonNull(actualVC, "actualVC");
        Utils.nonNull(expectedVC, "expectedVC");
        Assert.assertEquals(actualVC.size(), expectedVC.size());
        for (int i = 0; i < actualVC.size(); i++) {
            VariantContextTestUtils.assertVariantContextsAreEqual(actualVC.get(i), expectedVC.get(i), Arrays.asList("QD"), Collections.emptyList());
        }
    }

    @Test (dataProvider = "VCFdata")
    public void testUsingGenomicsDB(File[] inputs, File expected, File expectedDb, List<SimpleInterval> intervals, List<String> additionalArguments, String reference) throws IOException {
        Utils.resetRandomGenerator();
        //we need this for reproducibility because of random number generation in fixTooHighQD

        final File tempGenomicsDB = GenomicsDBTestUtils.createTempGenomicsDB(Arrays.asList(inputs), IntervalUtils.getSpanningInterval(intervals));
        final String genomicsDBUri = GenomicsDBTestUtils.makeGenomicsDBUri(tempGenomicsDB);

        final File output = createTempFile("GnarlyGenotyper", ".vcf");
        final File outputDatabase = createTempFile("GnarlyGenotyper.annotationDatabase", ".vcf");

        final File outputResolved = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected : output;
        final File outputDatabaseResolved = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedDb : outputDatabase;

        runTool(genomicsDBUri, intervals, reference, outputResolved, outputDatabaseResolved, additionalArguments);

        // Test for an exact match against past results
        if (!UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS) {
            final List<VariantContext> expectedVC = getVariantContexts(expected);
            final List<VariantContext> actualVC = getVariantContexts(output);
            assertEqualVariants(actualVC, expectedVC);
            if (expectedDb != null) {
                final List<VariantContext> expectedDB = getVariantContexts(expectedDb);
                final List<VariantContext> actualDB = getVariantContexts(outputDatabase);
                assertEqualVariants(actualDB, expectedDB);
            }
        }
    }

    protected void runTool(String input, List<SimpleInterval> intervals, String reference, File output, File outputDatabase, List<String> additionalArguments) {

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(reference))
        .add("V", input)
        .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");
        args.addOutput(output);
        intervals.forEach(args::addInterval);

        if (outputDatabase != null) {
            // to avoid null arguments when updating test resources using UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = true
            args.add("output-db", outputDatabase.getAbsolutePath());
        }

        additionalArguments.forEach(args::addRaw);

        runCommandLine(args);
    }

    @Test  //has a site with lots of alternate alleles, but Hail output is just sites-only
    public void testOnHailOutput() {
        final String input = getToolTestDataDir() + "hailOutput.chr20snippet.sites_only.vcf";
        final File output = createTempFile("GnarlyGenotyper", ".vcf");
        final File expected = new File(getToolTestDataDir() + "expected.hailOutput.chr20snippet.sites_only.vcf");

        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(hg38Reference))
                .add("V", input)
                .add("L", "chr20:10000000-10030000")
                .add("only-output-calls-starting-in-intervals", true)
                .add("keep-all-sites", true)
                .addOutput(outputPath)
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");
        runCommandLine(args);

        // Test for an exact match against past results
        if (!UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS) {
            //should have same variants as input with low QUAL variants marked with a monomorphic filter
            //QUAL column filled in, retains SB_TABLE, MQ_DP (deprecated), VarDP
            //does not expect ExcessHet because legacy input has neither genotypes nor GT counts
            try (final FeatureDataSource<VariantContext> actualVcs = new FeatureDataSource<>(output);
                 final FeatureDataSource<VariantContext> expectedVcs = new FeatureDataSource<>(expected)) {
                GATKBaseTest.assertCondition(actualVcs, expectedVcs,
                        (a, e) -> VariantContextTestUtils.assertVariantContextsAreEqual(a, e,
                                Collections.emptyList(), Collections.emptyList()));
            }
        }
    }

    //as for issue #7483 where two GVCFs merged with GenomicsDB produce some empty annotations fields
    @Test
    public void testOnEmptyAnnotations() {
        final String input = getToolTestDataDir() + "emptyASAnnotations.g.vcf";
        final File output = createTempFile("GnarlyGenotyper", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(hg38Reference))
                .add("V", input)
                .addOutput(output)
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");
        runCommandLine(args);

        final Pair<VCFHeader, List<VariantContext>> outputData = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        Assert.assertEquals(outputData.getRight().size(), 1);
        final VariantContext vc = outputData.getRight().get(0);
        Assert.assertTrue(vc.hasAttribute(GATKVCFConstants.AS_RMS_MAPPING_QUALITY_KEY));
        final List<String> mqs = vc.getAttributeAsStringList(GATKVCFConstants.AS_RMS_MAPPING_QUALITY_KEY,"");
        Assert.assertEquals(mqs.get(1), VCFConstants.MISSING_VALUE_v4);
        Assert.assertTrue(vc.hasAttribute(GATKVCFConstants.AS_FISHER_STRAND_KEY));
        final List<String> fss = vc.getAttributeAsStringList(GATKVCFConstants.AS_RMS_MAPPING_QUALITY_KEY,"");
        Assert.assertEquals(fss.get(1), VCFConstants.MISSING_VALUE_v4);
        Assert.assertTrue(vc.hasAttribute(GATKVCFConstants.AS_STRAND_ODDS_RATIO_KEY));
        final List<String> sors = vc.getAttributeAsStringList(GATKVCFConstants.AS_STRAND_ODDS_RATIO_KEY,"");
        Assert.assertEquals(sors.get(1), VCFConstants.MISSING_VALUE_v4);
    }

    @Test
    public void testHaploidInput() {
        final File haploidGVCF = new File(getToolTestDataDir(), "chrY_haploid_dragen.g.vcf");
        final File output = createTempFile("GnarlyGenotyper", ".vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(hg38Reference))
                .add("V", haploidGVCF)
                .addOutput(output)
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");
        runCommandLine(args);

        final Pair<VCFHeader, List<VariantContext>> outputData = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        Assert.assertEquals(outputData.getRight().size(), 1);
        final VariantContext vc = outputData.getRight().get(0);
        Assert.assertEquals(vc.getGenotypes().size(), 1);
        final Genotype g = vc.getGenotype(0);
        Assert.assertEquals(g.getPloidy(), 1);
        Assert.assertEquals(g.getAlleles().get(0), Allele.ALT_A);
        Assert.assertTrue(g.hasPL());
        Assert.assertEquals(g.getPL().length, 2);
        Assert.assertEquals(g.getPL(), new int[]{83,0});
    }

    @Test
    public void testNoReadsOutputAsNoCall() {
        //note these are b37 data
        File no_reads = new File(toolsTestDir, "/walkers/GnarlyGenotyper/testNoReads.rb.g.vcf");
        //this is an artisanal, hand-crafted VCF with a QUAL approx that's been artificially enhanced
        File fake_variant = new File(toolsTestDir, "/walkers/GnarlyGenotyper/fake_sample2.rb.g.vcf");
        final SimpleInterval interval =  new SimpleInterval("20", 10000000, 10000000);
        File tempGdb = GenomicsDBTestUtils.createTempGenomicsDB(Arrays.asList(no_reads, fake_variant), interval);
        final String genomicsDBUri = GenomicsDBTestUtils.makeGenomicsDBUri(tempGdb);

        final File output = createTempFile("checkNoCall", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(b37_reference_20_21)
                .addVCF(genomicsDBUri)
                .addInterval(interval)
                .addOutput(output);
        runCommandLine(args);

        final Pair<VCFHeader, List<VariantContext>> outputData = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        Assert.assertEquals(outputData.getRight().size(), 1);
        final VariantContext vc = outputData.getRight().get(0);
        Assert.assertEquals(vc.getGenotypes().size(), 2);
        final Genotype g = vc.getGenotype("GTEX-RVPV-0003");
        Assert.assertTrue(g.isNoCall()); //most importantly!
        Assert.assertFalse(g.hasGQ());
        Assert.assertFalse(g.hasDP());
        Assert.assertFalse(g.hasAD());
        Assert.assertFalse(g.hasPL());
    }

    @Test
    public void gnarlyHomRefNullPLs() {
        final File output = createTempFile("nullPLs", ".vcf");
        final File inputNull = new File(toolsTestDir, "walkers/GnarlyGenotyper/NA20890.rb.g.vcf");
        final File inputVariant = new File(toolsTestDir, "walkers/GnarlyGenotyper/NA20846.rb.g.vcf");
        final SimpleInterval interval =  new SimpleInterval("chr1", 13273, 13273);

        File tempGdb = GenomicsDBTestUtils.createTempGenomicsDB(Arrays.asList(inputNull, inputVariant), interval);
        final String genomicsDBUri = GenomicsDBTestUtils.makeGenomicsDBUri(tempGdb);

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(hg38Reference)
                .addVCF(genomicsDBUri)
                .addInterval(interval)
                .addOutput(output);
        runCommandLine(args);

        final Pair<VCFHeader, List<VariantContext>> outputData = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        Assert.assertEquals(outputData.getRight().size(), 1);
        final Genotype nullPlGenotype = outputData.getRight().get(0).getGenotype("NA20890");
        Assert.assertEquals(nullPlGenotype.getPL(), null);
        Assert.assertEquals(nullPlGenotype.getAlleles(), Arrays.asList(Allele.create("G", true), Allele.create("G", true)));
        Assert.assertEquals(nullPlGenotype.getGQ(), 0);
        Assert.assertEquals(nullPlGenotype.getDP(), 4);
    }
}