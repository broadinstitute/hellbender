package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.AbstractConcordanceWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.broadinstitute.hellbender.tools.sv.cluster.ClusteringParameters;
import org.broadinstitute.hellbender.tools.sv.cluster.SVClusterEngineArgumentsCollection;
import org.broadinstitute.hellbender.tools.sv.concordance.ClosestSVFinder;
import org.broadinstitute.hellbender.tools.sv.concordance.SVConcordanceAnnotator;
import org.broadinstitute.hellbender.tools.sv.concordance.SVConcordanceLinkage;
import org.broadinstitute.hellbender.tools.walkers.validation.Concordance;
import org.broadinstitute.hellbender.tools.walkers.validation.ConcordanceState;
import org.broadinstitute.hellbender.utils.MathUtilsUnitTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

public class SVConcordanceIntegrationTest extends CommandLineProgramTest {

    private static final double TOLERANCE = 0.001;

    @Test
    public void testRefPanel() {
        final File output = createTempFile("concord", ".vcf.gz");
        final String evalVcfPath = getToolTestDataDir() + "ref_panel_1kg.cleaned.gatk.chr22_chrY.vcf.gz";
        /**
         * Test file produced from raw standardized VCFs from manta, melt, wham, and cnv callers with SVCluster parameters:
         * --algorithm SINGLE_LINKAGE
         * --pesr-interval-overlap 0.9
         * --pesr-breakend-window 50
         */
        final String truthVcfPath = getToolTestDataDir() + "ref_panel_1kg.raw_calls.chr22_chrY.vcf.gz";

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(output)
                .add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, GATKBaseTest.FULL_HG38_DICT)
                .add(SVClusterEngineArgumentsCollection.DEPTH_SAMPLE_OVERLAP_FRACTION_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.DEPTH_INTERVAL_OVERLAP_FRACTION_NAME, 0.5)
                .add(SVClusterEngineArgumentsCollection.DEPTH_BREAKEND_WINDOW_NAME, 2000)
                .add(SVClusterEngineArgumentsCollection.MIXED_SAMPLE_OVERLAP_FRACTION_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.MIXED_INTERVAL_OVERLAP_FRACTION_NAME, 0.1)
                .add(SVClusterEngineArgumentsCollection.MIXED_BREAKEND_WINDOW_NAME, 2000)
                .add(SVClusterEngineArgumentsCollection.PESR_SAMPLE_OVERLAP_FRACTION_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.PESR_INTERVAL_OVERLAP_FRACTION_NAME, 0.1)
                .add(SVClusterEngineArgumentsCollection.PESR_BREAKEND_WINDOW_NAME, 500)
                .add(AbstractConcordanceWalker.TRUTH_VARIANTS_LONG_NAME, truthVcfPath)
                .add(AbstractConcordanceWalker.EVAL_VARIANTS_SHORT_NAME, evalVcfPath);


        runCommandLine(args, SVConcordance.class.getSimpleName());

        final Pair<VCFHeader, List<VariantContext>> outputVcf = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());

        final SAMSequenceDictionary dictionary = SVTestUtils.hg38Dict;
        final ClusteringParameters depthParameters = ClusteringParameters.createDepthParameters(0.5, 0, 2000, 0);
        final ClusteringParameters mixedParameters = ClusteringParameters.createMixedParameters(0.1, 0, 2000, 0);
        final ClusteringParameters pesrParameters = ClusteringParameters.createPesrParameters(0.1, 0, 500, 0);
        final SVConcordanceLinkage linkage = new SVConcordanceLinkage(dictionary);
        linkage.setDepthOnlyParams(depthParameters);
        linkage.setMixedParams(mixedParameters);
        linkage.setEvidenceParams(pesrParameters);
        final SVConcordanceAnnotator annotator = new SVConcordanceAnnotator();

        final List<SVCallRecord> inputEvalVariants = VariantContextTestUtils.readEntireVCFIntoMemory(evalVcfPath).getValue()
                .stream().map(SVCallRecordUtils::create).collect(Collectors.toList());
        final List<SVCallRecord> inputTruthVariants = VariantContextTestUtils.readEntireVCFIntoMemory(truthVcfPath).getValue()
                .stream().map(SVCallRecordUtils::create).collect(Collectors.toList());

        final ClosestSVFinder finder = new ClosestSVFinder(linkage, annotator::annotate, dictionary);
        final List<SVCallRecord> expectedRecords = new ArrayList<>(inputEvalVariants.size());
        final Set<Map.Entry<Long, SVCallRecord>> truthRecords = new HashSet<>();
        Long itemId = 0L;
        for (final SVCallRecord record : inputTruthVariants) {
            truthRecords.add(new AbstractMap.SimpleImmutableEntry<>(itemId++, record));
        }
        for (final SVCallRecord evalRecord : inputEvalVariants) {
            final Map.Entry<Long, SVCallRecord> closestRecord = finder.getClosestItem(evalRecord, truthRecords);
            final ClosestSVFinder.ClosestPair cluster =
                    new ClosestSVFinder.ClosestPair(evalRecord, closestRecord == null ? null : closestRecord.getValue());
            expectedRecords.add(annotator.annotate(cluster));
        }

        final Comparator<VariantContext> idComparator = Comparator.comparing(VariantContext::getID);
        final List<VariantContext> expectedVariants = expectedRecords.stream()
                .map(SVCallRecordUtils::getVariantBuilder)
                .map(VariantContextBuilder::make)
                .sorted(idComparator)
                .collect(Collectors.toList());
        final List<VariantContext> outputVariants = outputVcf.getValue().stream().sorted(idComparator).collect(Collectors.toList());
        final List<VariantContext> idSortedInputEvalVariants = inputEvalVariants.stream()
                .map(SVCallRecordUtils::getVariantBuilder)
                .map(VariantContextBuilder::make)
                .sorted(idComparator)
                .collect(Collectors.toList());

        Assert.assertEquals(outputVariants.size(), idSortedInputEvalVariants.size());
        Assert.assertEquals(outputVariants.size(), expectedVariants.size());
        final Set<String> checkedVariantsSet = new HashSet<>();
        for (int i = 0; i < outputVariants.size(); i++) {
            final VariantContext outputVariant = outputVariants.get(i);
            final VariantContext expectedVariant = expectedVariants.get(i);
            Assert.assertEquals(outputVariant.getID(), idSortedInputEvalVariants.get(i).getID());
            Assert.assertEquals(outputVariant.getID(), expectedVariant.getID());
            Assert.assertEquals(outputVariant.getContig(), expectedVariant.getContig());
            Assert.assertEquals(outputVariant.getStart(), expectedVariant.getStart());
            Assert.assertEquals(outputVariant.getEnd(), expectedVariant.getEnd());
            Assert.assertEquals(outputVariant.getAlleles(), expectedVariant.getAlleles());
            Assert.assertEquals(outputVariant.getFilters(), expectedVariant.getFilters());
            Assert.assertEquals(outputVariant.getPhredScaledQual(), expectedVariant.getPhredScaledQual());
            Assert.assertEquals(outputVariant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "test_default"), expectedVariant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "expected_default"));
            checkTruthVariantId(outputVariant, expectedVariant);
            Assert.assertEquals(outputVariant.getAttributeAsString(Concordance.TRUTH_STATUS_VCF_ATTRIBUTE, "test_default"), expectedVariant.getAttributeAsString(Concordance.TRUTH_STATUS_VCF_ATTRIBUTE, "expected_default"));
            if (outputVariant.getID().equals("ref_panel_1kg.chr22.final_cleanup_DEL_chr22_1")) {
                checkVariant(outputVariant, expectedVariant, "ref_panel_1kg_raw_00000062",
                        268./312, 66./312, checkedVariantsSet);
            } else if (outputVariant.getID().equals("ref_panel_1kg.chr22.final_cleanup_DEL_chr22_140")) {
                checkVariant(outputVariant, expectedVariant, "ref_panel_1kg_raw_000017fc",
                        310./312, 3./312, checkedVariantsSet);
            } else if (outputVariant.getID().equals("ref_panel_1kg.chr22.final_cleanup_INS_chr22_100")) {
                checkVariant(outputVariant, expectedVariant, "ref_panel_1kg_raw_00001b78",
                        1., 1./312, checkedVariantsSet);
            }
        }
        Assert.assertEquals(checkedVariantsSet.size(), 3);
        Assert.assertTrue(checkedVariantsSet.contains("ref_panel_1kg.chr22.final_cleanup_DEL_chr22_1"));
        Assert.assertTrue(checkedVariantsSet.contains("ref_panel_1kg.chr22.final_cleanup_DEL_chr22_140"));
        Assert.assertTrue(checkedVariantsSet.contains("ref_panel_1kg.chr22.final_cleanup_INS_chr22_100"));
    }

    private static void checkTruthVariantId(final VariantContext actual, final VariantContext expected) {
        final String expectedTruthId = expected.getAttributeAsString(GATKSVVCFConstants.TRUTH_VARIANT_ID_INFO, VCFConstants.EMPTY_INFO_FIELD);
        final String actualTruthId = actual.getAttributeAsString(GATKSVVCFConstants.TRUTH_VARIANT_ID_INFO, VCFConstants.EMPTY_INFO_FIELD);
        Assert.assertEquals(actualTruthId, expectedTruthId);
    }

    private static void checkVariant(final VariantContext actual, final VariantContext expected,
                                     final String expectedTruthId, final double expectedGenotypeConcordance,
                                     final double expectedAF,
                                     final Set<String> checkedVariantsSet) {
        Assert.assertFalse(checkedVariantsSet.contains(actual.getID()));
        checkedVariantsSet.add(actual.getID());
        Assert.assertEquals(actual.getAttributeAsString(GATKSVVCFConstants.TRUTH_VARIANT_ID_INFO, ""), expectedTruthId);
        Assert.assertEquals(expected.getAttributeAsString(GATKSVVCFConstants.TRUTH_VARIANT_ID_INFO, ""), expectedTruthId);
        MathUtilsUnitTest.assertEqualsDoubleSmart(actual.getAttributeAsDouble(GATKSVVCFConstants.GENOTYPE_CONCORDANCE_INFO, -1.), expectedGenotypeConcordance, TOLERANCE);
        MathUtilsUnitTest.assertEqualsDoubleSmart(expected.getAttributeAsDouble(GATKSVVCFConstants.GENOTYPE_CONCORDANCE_INFO, -1.), expectedGenotypeConcordance, TOLERANCE);
        MathUtilsUnitTest.assertEqualsDoubleSmart(actual.getAttributeAsDouble(GATKSVVCFConstants.GENOTYPE_CONCORDANCE_INFO, -1.), expected.getAttributeAsDouble(GATKSVVCFConstants.GENOTYPE_CONCORDANCE_INFO, -1.), TOLERANCE);
        MathUtilsUnitTest.assertEqualsDoubleSmart(expected.getAttributeAsDouble(GATKSVVCFConstants.GENOTYPE_CONCORDANCE_INFO, -1.), expected.getAttributeAsDouble(GATKSVVCFConstants.GENOTYPE_CONCORDANCE_INFO, -1.), TOLERANCE);
        MathUtilsUnitTest.assertEqualsDoubleSmart(actual.getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, -1.), expectedAF, TOLERANCE);
        MathUtilsUnitTest.assertEqualsDoubleSmart(expected.getAttributeAsDoubleList(VCFConstants.ALLELE_FREQUENCY_KEY, -1.).get(0).doubleValue(), expectedAF, TOLERANCE);
    }

    @Test
    public void testSelf() {

        // Run a vcf against itself and check that every variant matched itself

        final File output = createTempFile("concord", ".vcf.gz");
        final String vcfPath = getToolTestDataDir() + "ref_panel_1kg.cleaned.gatk.chr22_chrY.vcf.gz";

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(output)
                .add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, GATKBaseTest.FULL_HG38_DICT)
                .add(SVClusterEngineArgumentsCollection.DEPTH_SAMPLE_OVERLAP_FRACTION_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.DEPTH_INTERVAL_OVERLAP_FRACTION_NAME, 0.5)
                .add(SVClusterEngineArgumentsCollection.DEPTH_SIZE_SIMILARITY_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.DEPTH_BREAKEND_WINDOW_NAME, 2000)
                .add(SVClusterEngineArgumentsCollection.MIXED_SAMPLE_OVERLAP_FRACTION_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.MIXED_INTERVAL_OVERLAP_FRACTION_NAME, 0.1)
                .add(SVClusterEngineArgumentsCollection.MIXED_SIZE_SIMILARITY_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.MIXED_BREAKEND_WINDOW_NAME, 2000)
                .add(SVClusterEngineArgumentsCollection.PESR_SAMPLE_OVERLAP_FRACTION_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.PESR_INTERVAL_OVERLAP_FRACTION_NAME, 0.1)
                .add(SVClusterEngineArgumentsCollection.PESR_SIZE_SIMILARITY_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.PESR_BREAKEND_WINDOW_NAME, 500)
                .add(AbstractConcordanceWalker.TRUTH_VARIANTS_LONG_NAME, vcfPath)
                .add(AbstractConcordanceWalker.EVAL_VARIANTS_SHORT_NAME, vcfPath);

        runCommandLine(args, SVConcordance.class.getSimpleName());

        final Pair<VCFHeader, List<VariantContext>> outputVcf = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());

        final SAMSequenceDictionary dictionary = SVTestUtils.hg38Dict;
        final ClusteringParameters depthParameters = ClusteringParameters.createDepthParameters(0.5, 0, 2000, 0);
        final ClusteringParameters mixedParameters = ClusteringParameters.createMixedParameters(0.1, 0, 2000, 0);
        final ClusteringParameters pesrParameters = ClusteringParameters.createPesrParameters(0.1, 0, 500, 0);
        final SVConcordanceLinkage linkage = new SVConcordanceLinkage(dictionary);
        linkage.setDepthOnlyParams(depthParameters);
        linkage.setMixedParams(mixedParameters);
        linkage.setEvidenceParams(pesrParameters);

        final Comparator<VariantContext> idComparator = Comparator.comparing(VariantContext::getID);
        final List<VariantContext> inputVariants = VariantContextTestUtils.readEntireVCFIntoMemory(vcfPath).getValue()
                .stream()
                .sorted(idComparator)
                .collect(Collectors.toList());
        final List<VariantContext> outputVariants = outputVcf.getValue().stream().sorted(idComparator).collect(Collectors.toList());

        Assert.assertEquals(outputVariants.size(), inputVariants.size());
        for (int i = 0; i < outputVariants.size(); i++) {
            final VariantContext outputVariant = outputVariants.get(i);
            final VariantContext expectedVariant = inputVariants.get(i);
            Assert.assertEquals(outputVariant.getID(), expectedVariant.getID());
            Assert.assertEquals(outputVariant.getContig(), expectedVariant.getContig());
            Assert.assertEquals(outputVariant.getStart(), expectedVariant.getStart());
            Assert.assertEquals(outputVariant.getEnd(), expectedVariant.getEnd());
            Assert.assertEquals(outputVariant.getAlleles(), expectedVariant.getAlleles());
            Assert.assertEquals(outputVariant.getFilters(), expectedVariant.getFilters());
            Assert.assertEquals(outputVariant.getPhredScaledQual(), expectedVariant.getPhredScaledQual());
            final String svtype = outputVariant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "test_default");
            Assert.assertEquals(svtype, expectedVariant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "expected_default"));
            // check the variant matched itself with perfect concordance
            Assert.assertEquals(outputVariant.getAttributeAsString(Concordance.TRUTH_STATUS_VCF_ATTRIBUTE, "test_default"), ConcordanceState.TRUE_POSITIVE.getAbbreviation());
            Assert.assertEquals(outputVariant.getAttributeAsString(GATKSVVCFConstants.TRUTH_VARIANT_ID_INFO, ""), outputVariant.getID());
            if (svtype.equals(GATKSVVCFConstants.StructuralVariantAnnotationType.CNV.toString())) {
                MathUtilsUnitTest.assertEqualsDoubleSmart(outputVariant.getAttributeAsDouble(GATKSVVCFConstants.COPY_NUMBER_CONCORDANCE_INFO, -1.), 1.0, TOLERANCE);
            } else {
                MathUtilsUnitTest.assertEqualsDoubleSmart(outputVariant.getAttributeAsDouble(GATKSVVCFConstants.GENOTYPE_CONCORDANCE_INFO, -1.), 1.0, TOLERANCE);
            }
        }
    }

    @Test
    public void testSelfTruthSubset() {

        // Run a vcf against itself but subsetted to fewer samples and variants

        final File output = createTempFile("concord", ".vcf");
        final String evalVcfPath = getToolTestDataDir() + "ref_panel_1kg.cleaned.gatk.chr22_chrY.vcf.gz";
        final String truthVcfPath = getToolTestDataDir() + "ref_panel_1kg.cleaned.gatk.chr22_chrY.sample_subset.vcf.gz";

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(output)
                .add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, GATKBaseTest.FULL_HG38_DICT)
                .add(AbstractConcordanceWalker.EVAL_VARIANTS_LONG_NAME, evalVcfPath)
                .add(AbstractConcordanceWalker.TRUTH_VARIANTS_LONG_NAME, truthVcfPath);

        runCommandLine(args, SVConcordance.class.getSimpleName());
        assertPerfectConcordance(output, evalVcfPath);
    }

    @Test
    public void testSelfEvalSubset() {

        // Run a vcf against itself but with extra samples and variants

        final File output = createTempFile("concord", ".vcf.gz");
        final String evalVcfPath = getToolTestDataDir() + "ref_panel_1kg.cleaned.gatk.chr22_chrY.sample_subset.vcf.gz";
        final String truthVcfPath = getToolTestDataDir() + "ref_panel_1kg.cleaned.gatk.chr22_chrY.vcf.gz";

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(output)
                .add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, GATKBaseTest.FULL_HG38_DICT)
                .add(AbstractConcordanceWalker.EVAL_VARIANTS_LONG_NAME, evalVcfPath)
                .add(AbstractConcordanceWalker.TRUTH_VARIANTS_LONG_NAME, truthVcfPath);

        runCommandLine(args, SVConcordance.class.getSimpleName());
        assertPerfectConcordance(output, evalVcfPath);
    }

    private void assertPerfectConcordance(final File output, final String evalVcfPath) {

        final Pair<VCFHeader, List<VariantContext>> outputVcf = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());

        final SAMSequenceDictionary dictionary = SVTestUtils.hg38Dict;
        final ClusteringParameters depthParameters = ClusteringParameters.createDepthParameters(0.5, 0, 2000, 0);
        final ClusteringParameters mixedParameters = ClusteringParameters.createMixedParameters(0.1, 0, 2000, 0);
        final ClusteringParameters pesrParameters = ClusteringParameters.createPesrParameters(0.1, 0, 500, 0);
        final SVConcordanceLinkage linkage = new SVConcordanceLinkage(dictionary);
        linkage.setDepthOnlyParams(depthParameters);
        linkage.setMixedParams(mixedParameters);
        linkage.setEvidenceParams(pesrParameters);

        final Comparator<VariantContext> idComparator = Comparator.comparing(VariantContext::getID);

        final List<VariantContext> inputVariants = VariantContextTestUtils.readEntireVCFIntoMemory(evalVcfPath).getValue()
                .stream()
                .sorted(idComparator)
                .collect(Collectors.toList());
        final List<VariantContext> outputVariants = outputVcf.getValue().stream().sorted(idComparator).collect(Collectors.toList());

        Assert.assertEquals(outputVariants.size(), inputVariants.size());
        for (int i = 0; i < outputVariants.size(); i++) {
            final VariantContext outputVariant = outputVariants.get(i);
            final VariantContext expectedVariant = inputVariants.get(i);
            Assert.assertEquals(outputVariant.getID(), expectedVariant.getID());
            Assert.assertEquals(outputVariant.getContig(), expectedVariant.getContig());
            Assert.assertEquals(outputVariant.getStart(), expectedVariant.getStart());
            Assert.assertEquals(outputVariant.getEnd(), expectedVariant.getEnd());
            Assert.assertEquals(outputVariant.getAlleles(), expectedVariant.getAlleles());
            Assert.assertEquals(outputVariant.getFilters(), expectedVariant.getFilters());
            Assert.assertEquals(outputVariant.getPhredScaledQual(), expectedVariant.getPhredScaledQual());
            final String svtype = outputVariant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "test_default");
            Assert.assertEquals(svtype, expectedVariant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "expected_default"));
            // check the variant matched itself with perfect concordance
            if (outputVariant.getAttributeAsString(Concordance.TRUTH_STATUS_VCF_ATTRIBUTE, "test_default") == ConcordanceState.TRUE_POSITIVE.getAbbreviation()) {
                Assert.assertEquals(outputVariant.getAttributeAsString(GATKSVVCFConstants.TRUTH_VARIANT_ID_INFO, ""), outputVariant.getID());
                if (svtype.equals(GATKSVVCFConstants.StructuralVariantAnnotationType.CNV.toString())) {
                    MathUtilsUnitTest.assertEqualsDoubleSmart(outputVariant.getAttributeAsDouble(GATKSVVCFConstants.COPY_NUMBER_CONCORDANCE_INFO, -1.), 1.0, TOLERANCE);
                } else {
                    MathUtilsUnitTest.assertEqualsDoubleSmart(outputVariant.getAttributeAsDouble(GATKSVVCFConstants.GENOTYPE_CONCORDANCE_INFO, -1.), 1.0, TOLERANCE);
                }
            }
        }
    }

    @Test
    public void testSitesOnly() {
        final File output = createTempFile("concord_sites_only", ".vcf.gz");
        final String evalVcfPath = getToolTestDataDir() + "ref_panel_1kg.cleaned.gatk.chr22_chrY.sites_only.vcf.gz";
        final String truthVcfPath = getToolTestDataDir() + "ref_panel_1kg.raw_calls.chr22_chrY.sites_only.vcf.gz";

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(output)
                .add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, GATKBaseTest.FULL_HG38_DICT)
                .add(SVClusterEngineArgumentsCollection.DEPTH_SAMPLE_OVERLAP_FRACTION_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.DEPTH_INTERVAL_OVERLAP_FRACTION_NAME, 0.5)
                .add(SVClusterEngineArgumentsCollection.DEPTH_BREAKEND_WINDOW_NAME, 2000)
                .add(SVClusterEngineArgumentsCollection.MIXED_SAMPLE_OVERLAP_FRACTION_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.MIXED_INTERVAL_OVERLAP_FRACTION_NAME, 0.1)
                .add(SVClusterEngineArgumentsCollection.MIXED_BREAKEND_WINDOW_NAME, 2000)
                .add(SVClusterEngineArgumentsCollection.PESR_SAMPLE_OVERLAP_FRACTION_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.PESR_INTERVAL_OVERLAP_FRACTION_NAME, 0.1)
                .add(SVClusterEngineArgumentsCollection.PESR_BREAKEND_WINDOW_NAME, 500)
                .add(AbstractConcordanceWalker.TRUTH_VARIANTS_LONG_NAME, truthVcfPath)
                .add(AbstractConcordanceWalker.EVAL_VARIANTS_SHORT_NAME, evalVcfPath);


        runCommandLine(args, SVConcordance.class.getSimpleName());

        final Pair<VCFHeader, List<VariantContext>> outputVcf = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        final List<SVCallRecord> inputEvalVariants = VariantContextTestUtils.readEntireVCFIntoMemory(evalVcfPath).getValue()
                .stream().map(SVCallRecordUtils::create).collect(Collectors.toList());
        Assert.assertEquals(outputVcf.getValue().size(), inputEvalVariants.size());
        final long numTruePositives = outputVcf.getValue().stream().filter(v -> v.getAttributeAsString(Concordance.TRUTH_STATUS_VCF_ATTRIBUTE, "").equals(ConcordanceState.TRUE_POSITIVE.getAbbreviation())).count();
        Assert.assertEquals((int) numTruePositives, 1104);
        for (final VariantContext variant: outputVcf.getValue()) {
            // Concordance can't be calculated for sites-only
            Assert.assertTrue(variant.getAttributeAsString(GATKSVVCFConstants.GENOTYPE_CONCORDANCE_INFO, ".").equals("."));
        }
    }

    @Test(expectedExceptions = UserException.IncompatibleSequenceDictionaries.class)
    public void testHeaderContigsOutOfOrder() {
        final File output = createTempFile("concord_sites_only", ".vcf.gz");
        final String evalVcfPath = getToolTestDataDir() + "ref_panel_1kg.cleaned.gatk.chr22_chrY.vcf.gz";
        final String truthVcfPath = getToolTestDataDir() + "ref_panel_1kg.cleaned.gatk.chr22_chrY.contigs_out_of_order.vcf.gz";

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(output)
                .add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, GATKBaseTest.FULL_HG38_DICT)
                .add(AbstractConcordanceWalker.TRUTH_VARIANTS_LONG_NAME, truthVcfPath)
                .add(AbstractConcordanceWalker.EVAL_VARIANTS_SHORT_NAME, evalVcfPath);

        runCommandLine(args, SVConcordance.class.getSimpleName());
    }

    @Test
    public void testSelfUnsorted() {

        // Run a vcf against itself but don't sort the output

        final File output = createTempFile("concord", ".vcf");
        final String evalVcfPath = getToolTestDataDir() + "ref_panel_1kg.cleaned.gatk.chr22_chrY.vcf.gz";
        final String truthVcfPath = getToolTestDataDir() + "ref_panel_1kg.cleaned.gatk.chr22_chrY.vcf.gz";

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(output)
                .add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, GATKBaseTest.FULL_HG38_DICT)
                .add(AbstractConcordanceWalker.EVAL_VARIANTS_LONG_NAME, evalVcfPath)
                .add(AbstractConcordanceWalker.TRUTH_VARIANTS_LONG_NAME, truthVcfPath)
                .add(SVConcordance.UNSORTED_OUTPUT_LONG_NAME, true);

        runCommandLine(args, SVConcordance.class.getSimpleName());
        assertPerfectConcordance(output, evalVcfPath);
    }
}
