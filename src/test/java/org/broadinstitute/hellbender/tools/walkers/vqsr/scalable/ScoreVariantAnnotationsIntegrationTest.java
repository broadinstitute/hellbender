package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import com.google.common.collect.Lists;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.LabeledVariantAnnotationsData;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.VariantType;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.VariantAnnotationsModelBackend;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutorException;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collectors;

/**
 * See documentation for {@link ExtractVariantAnnotationsIntegrationTest} for information about how inputs and
 * expected outputs used there are related to those used here and in {@link TrainVariantAnnotationsModelIntegrationTest}.
 */
public final class ScoreVariantAnnotationsIntegrationTest extends CommandLineProgramTest {

    // If true, update the expected outputs in tests that assert an exact match vs. prior output,
    // instead of actually running the tests. Can be used with "./gradlew test -Dtest.single=ScoreVariantAnnotationsIntegrationTest"
    // to update all of the exact-match tests at once. After you do this, you should look at the
    // diffs in the new expected outputs in git to confirm that they are consistent with expectations.
    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;

    /**
     * Make sure that someone didn't leave the UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS toggle turned on.
     */
    @Test
    public void assertThatExpectedOutputUpdateToggleIsDisabled() {
        Assert.assertFalse(UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS, "The toggle to update expected outputs should not be left enabled.");
    }

    private static final double CALIBRATION_SENSITIVITY_THRESHOLD = 0.9;

    private static final File PACKAGE_TEST_FILES_DIR = new File(largeFileTestDir,
            "org/broadinstitute/hellbender/tools/walkers/vqsr/scalable/");
    private static final File TEST_FILES_DIR = new File(largeFileTestDir,
            "org/broadinstitute/hellbender/tools/walkers/vqsr/scalable/score");
    private static final File INPUT_FROM_TRAIN_EXPECTED_TEST_FILES_DIR = new File(largeFileTestDir,
            "org/broadinstitute/hellbender/tools/walkers/vqsr/scalable/train/expected");
    private static final File EXPECTED_TEST_FILES_DIR = new File(TEST_FILES_DIR, "expected");

    private static final File ISOLATION_FOREST_PYTHON_SCRIPT = IOUtils.writeTempResource(
            new Resource("isolation-forest.py", TrainVariantAnnotationsModel.class));

    private static final File INPUT_VCF = new File(PACKAGE_TEST_FILES_DIR, "input/small_callset_low_threshold.sites-only.chr1.1-10M.vcf");

    // Supplier and functions for creating and adding various arguments to an ArgumentsBuilder.
    private static final Supplier<ArgumentsBuilder> BASE_ARGUMENTS_BUILDER_SUPPLIER = () -> {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addVCF(INPUT_VCF);
        argsBuilder.addFlag(LabeledVariantAnnotationsWalker.DO_NOT_GZIP_VCF_OUTPUT_LONG_NAME);
        argsBuilder.add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, false);
        return argsBuilder;
    };
    private static final BiFunction<ArgumentsBuilder, String, ArgumentsBuilder> ADD_MODEL_PREFIX = (argsBuilder, modelPrefix) -> {
        argsBuilder.add(ScoreVariantAnnotations.MODEL_PREFIX_LONG_NAME, modelPrefix);
        return argsBuilder;
    };
    private static final BiFunction<ArgumentsBuilder, Double, ArgumentsBuilder> ADD_CALIBRATION_SENSITIVITY_THRESHOLD = (argsBuilder, calibrationSensitivityThreshold) -> {
        argsBuilder.add(ScoreVariantAnnotations.SNP_CALIBRATION_SENSITIVITY_THRESHOLD_LONG_NAME, calibrationSensitivityThreshold);
        argsBuilder.add(ScoreVariantAnnotations.INDEL_CALIBRATION_SENSITIVITY_THRESHOLD_LONG_NAME, calibrationSensitivityThreshold);
        return argsBuilder;
    };
    private static final BiFunction<ArgumentsBuilder, VariantAnnotationsModelBackend, ArgumentsBuilder> ADD_MODEL_BACKEND = (argsBuilder, modelBackendMode) -> {
        argsBuilder.add(TrainVariantAnnotationsModel.MODEL_BACKEND_LONG_NAME, modelBackendMode);
        return argsBuilder;
    };
    private static final Function<ArgumentsBuilder, ArgumentsBuilder> ADD_ISOLATION_FOREST_PYTHON_SCRIPT = argsBuilder -> {
        argsBuilder.add(ScoreVariantAnnotations.PYTHON_SCRIPT_LONG_NAME, ISOLATION_FOREST_PYTHON_SCRIPT);
        return argsBuilder;
    };

    /**
     * Exact-match tests for (non-exhaustive) configurations given by the Cartesian product of the following options:
     * 1) non-allele-specific ("nonAS") vs. allele-specific ("AS")
     * 2) model backend
     *      2a) Java Bayesian Gaussian Mixture Model (BGMM) backend TODO the BGMM has been reduced to a stub for this initial PR; subsequent PRs will cover the backend code and reconnect the stub
     *      2b) default PYTHON_IFOREST ("IF.score")
     *      2c) specified PYTHON_SCRIPT ("IF.score"); we will simply use the same script as the default PYTHON_IFOREST backend, so this is just a test of the command-line interface
     *      We should expect 2b-c to give functionally identical results.
     * 3) SNP-only ("snp") vs. SNP+INDEL ("snpIndel") (for both of these options, we use trained models that contain both SNP and INDEL scorers as input)
     *  TODO the BGMM has been reduced to a stub for this initial PR; subsequent PRs will cover the backend code and reconnect the stub
     */
    @DataProvider(name = "dataValidInputs")
    public Object[][] dataValidInputs() {
        final List<List<Pair<String, Function<ArgumentsBuilder, ArgumentsBuilder>>>> testConfigurations = Lists.cartesianProduct(
                Arrays.asList(
                        Pair.of("extract.nonAS.snpIndel.posUn.train.snpIndel.posOnly", Function.identity()),
                        Pair.of("extract.AS.snpIndel.posUn.train.snpIndel.posOnly", Function.identity())),
                Arrays.asList(
                        Pair.of("IF.score", ab -> ADD_MODEL_BACKEND.apply(ab, VariantAnnotationsModelBackend.PYTHON_IFOREST)), // this and the following case give the same results, so they are given the same IF.score tag
                        Pair.of("IF.score", ADD_ISOLATION_FOREST_PYTHON_SCRIPT
                                .andThen(ab -> ADD_MODEL_BACKEND.apply(ab, VariantAnnotationsModelBackend.PYTHON_SCRIPT)))),
                Arrays.asList(
                        Pair.of("snp", ExtractVariantAnnotationsIntegrationTest.ADD_SNP_MODE_AND_RESOURCES),
                        Pair.of("snpIndel", ExtractVariantAnnotationsIntegrationTest.ADD_SNP_MODE_AND_RESOURCES
                                .andThen(ExtractVariantAnnotationsIntegrationTest.ADD_INDEL_MODE_AND_RESOURCES))));

        return testConfigurations.stream()
                .map(tagAndAddFunctionPairs -> new Object[]{
                        tagAndAddFunctionPairs.stream().map(Pair::getLeft).collect(Collectors.joining(".")), // e.g., extract.nonAS.snpIndel.posUn.train.snpIndel.posOnly.IF.score.snp
                        tagAndAddFunctionPairs.stream().map(Pair::getRight)                                              // creates the corresponding ArgumentsBuilder
                                .reduce(Function.identity(), Function::andThen)                                          //  by stringing together functions that add the
                                .apply(BASE_ARGUMENTS_BUILDER_SUPPLIER.get())})                                          //  appropriate arguments
                .toArray(Object[][]::new);
    }

    /**
     * Checks expected outputs given a tag (e.g., "extract.nonAS.snpIndel.posUn.train.snpIndel.posOnly.IF.score.snp") and arguments corresponding to the
     * Cartesian products generated in {@link #dataValidInputs}.
     *
     * We perform exact-match tests of any HDF5 files produced using h5diff, which is insensitive to timestamps within the file.
     * We also perform exact-match tests of VCF files using diff. VCF indices may not be diff equivalent, so
     * we just check for their existence.
     */
    @Test(dataProvider = "dataValidInputs", groups = {"python"}) // python environment is required to run tool and to use h5diff for exact-match comparisons
    public void testValidInputs(final String tag,
                                final ArgumentsBuilder argsBuilder) {
        final File outputDir = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? EXPECTED_TEST_FILES_DIR : createTempDir("score");
        final String outputPrefix = String.format("%s/%s", outputDir, tag);
        argsBuilder.addOutput(outputPrefix);

        // add arguments for model prefix based on the
        // train tag (the portion of the tag preceding ".score", e.g., extract.nonAS.snpIndel.posUn.train.snp.posOnly.IF),
        // which gives the basename for the model files
        final String trainTag = tag.split(".score")[0];
        if (tag.contains("nonAS")) {
            ExtractVariantAnnotationsIntegrationTest.ADD_NON_ALLELE_SPECIFIC_ANNOTATIONS.apply(argsBuilder);
        } else {
            ExtractVariantAnnotationsIntegrationTest.ADD_ALLELE_SPECIFIC_ANNOTATIONS.apply(argsBuilder);
        }
        final String modelPrefix = new File(INPUT_FROM_TRAIN_EXPECTED_TEST_FILES_DIR, trainTag).toString();
        final Function<ArgumentsBuilder, ArgumentsBuilder> addModelPrefix = ab ->
                ADD_MODEL_PREFIX.apply(ab, modelPrefix);
        final Function<ArgumentsBuilder, ArgumentsBuilder> addCalibrationSensitivityThreshold = ab ->
                ADD_CALIBRATION_SENSITIVITY_THRESHOLD.apply(ab, CALIBRATION_SENSITIVITY_THRESHOLD);
        addModelPrefix.andThen(addCalibrationSensitivityThreshold).apply(argsBuilder);

        // TODO test use of sites-only VCF (output by extract tool) to label extracted sites

        runCommandLine(argsBuilder);

        if (!UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS) {
            assertExpectedOutputs(tag, outputPrefix);
        }
    }

    private static void assertExpectedOutputs(final String tag,
                                              final String outputPrefix) {
        // vcf.idx files are not reproducible
        SystemCommandUtilsTest.runDiff(
                String.format("%s/%s.vcf", EXPECTED_TEST_FILES_DIR, tag),
                String.format("%s.vcf", outputPrefix));

        SystemCommandUtilsTest.runH5Diff(
                String.format("%s/%s.annot.hdf5", EXPECTED_TEST_FILES_DIR, tag),
                String.format("%s.annot.hdf5", outputPrefix));

        SystemCommandUtilsTest.runH5Diff(
                String.format("%s/%s.scores.hdf5", EXPECTED_TEST_FILES_DIR, tag),
                String.format("%s.scores.hdf5", outputPrefix));
    }

    /**
     * In contrast to {@link ExtractVariantAnnotationsIntegrationTest#testNoResources}, the non-presence of
     * resources here does not really affect the output.
     */
    @Test(groups = {"python"}) // python environment is required to run tool
    public void testNoResources() {
        final File outputDir = createTempDir("score");
        final String outputPrefix = String.format("%s/test", outputDir);
        final ArgumentsBuilder argsBuilder = BASE_ARGUMENTS_BUILDER_SUPPLIER.get();
        argsBuilder.add(LabeledVariantAnnotationsWalker.MODE_LONG_NAME, VariantType.SNP)
                .addOutput(outputPrefix);
        final String modelPrefix = new File(INPUT_FROM_TRAIN_EXPECTED_TEST_FILES_DIR,
                "extract.nonAS.snpIndel.posUn.train.snp.posOnly.IF").toString();
        final Function<ArgumentsBuilder, ArgumentsBuilder> addModelPrefix = ab ->
                ADD_MODEL_PREFIX.apply(ab, modelPrefix);
        addModelPrefix
                .andThen(ExtractVariantAnnotationsIntegrationTest.ADD_NON_ALLELE_SPECIFIC_ANNOTATIONS)
                .apply(argsBuilder);
        runCommandLine(argsBuilder);
        Assert.assertTrue(new File(outputPrefix + ScoreVariantAnnotations.ANNOTATIONS_HDF5_SUFFIX).exists());
        Assert.assertTrue(new File(outputPrefix + ScoreVariantAnnotations.SCORES_HDF5_SUFFIX).exists());
        Assert.assertTrue(new File(outputPrefix + ".vcf").exists());
        Assert.assertTrue(new File(outputPrefix + ".vcf.idx").exists());
    }

    /**
     * If no variants are present in the input in the specified region, we do not create the scores or annotations HDF5 files.
     * This is because we cannot create HDF5 files with empty arrays/matrices.
     */
    @Test(groups = {"python"}) // python environment is required to run tool
    public void testNoVariantsInInput() {
        final File outputDir = createTempDir("score");
        final String outputPrefix = String.format("%s/test", outputDir);
        final ArgumentsBuilder argsBuilder = BASE_ARGUMENTS_BUILDER_SUPPLIER.get();
        argsBuilder.add(LabeledVariantAnnotationsWalker.MODE_LONG_NAME, VariantType.SNP)
                .add(StandardArgumentDefinitions.INTERVALS_LONG_NAME, "chr2") // the test input VCF does not have variants here
                .addOutput(outputPrefix);
        final String modelPrefix = new File(INPUT_FROM_TRAIN_EXPECTED_TEST_FILES_DIR,
                "extract.nonAS.snpIndel.posUn.train.snp.posOnly.IF").toString();
        final Function<ArgumentsBuilder, ArgumentsBuilder> addModelPrefix = ab ->
                ADD_MODEL_PREFIX.apply(ab, modelPrefix);
        addModelPrefix
                .andThen(ExtractVariantAnnotationsIntegrationTest.ADD_NON_ALLELE_SPECIFIC_ANNOTATIONS)
                .andThen(ExtractVariantAnnotationsIntegrationTest.ADD_SNP_MODE_AND_RESOURCES)
                .apply(argsBuilder);
        runCommandLine(argsBuilder);
        Assert.assertFalse(new File(outputPrefix + ScoreVariantAnnotations.ANNOTATIONS_HDF5_SUFFIX).exists());
        Assert.assertFalse(new File(outputPrefix + ScoreVariantAnnotations.SCORES_HDF5_SUFFIX).exists());
        Assert.assertTrue(new File(outputPrefix + ".vcf").exists());
        Assert.assertTrue(new File(outputPrefix + ".vcf.idx").exists());
    }

    @Test(expectedExceptions = PythonScriptExecutorException.class, groups = {"python"}) // python environment is required to run tool
    public void testAnnotationsDoNotMatchThoseUsedToTrainModel() {
        final File outputDir = createTempDir("score");
        final String outputPrefix = String.format("%s/test", outputDir);
        final ArgumentsBuilder argsBuilder = BASE_ARGUMENTS_BUILDER_SUPPLIER.get();
        argsBuilder.add(LabeledVariantAnnotationsWalker.MODE_LONG_NAME, VariantType.SNP)
                .addOutput(outputPrefix);
        final String modelPrefix = new File(INPUT_FROM_TRAIN_EXPECTED_TEST_FILES_DIR,
                "extract.nonAS.snpIndel.posUn.train.snp.posOnly.IF").toString();
        final Function<ArgumentsBuilder, ArgumentsBuilder> addModelPrefix = ab ->
                ADD_MODEL_PREFIX.apply(ab, modelPrefix);
        addModelPrefix
                .andThen(ExtractVariantAnnotationsIntegrationTest.ADD_NON_ALLELE_SPECIFIC_ANNOTATIONS)  // model was trained with non-AS annotations
                .andThen(ExtractVariantAnnotationsIntegrationTest.ADD_ALLELE_SPECIFIC_ANNOTATIONS)      // but we additionally specify AS annotations
                .andThen(ExtractVariantAnnotationsIntegrationTest.ADD_SNP_MODE_AND_RESOURCES)
                .apply(argsBuilder);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = UserException.class)
    public void testReservedSNPResourceLabel() {
        final File outputDir = createTempDir("extract");
        final String outputPrefix = String.format("%s/test", outputDir);
        final ArgumentsBuilder argsBuilder = ExtractVariantAnnotationsIntegrationTest.ADD_NON_ALLELE_SPECIFIC_ANNOTATIONS.apply(BASE_ARGUMENTS_BUILDER_SUPPLIER.get());
        argsBuilder.add(LabeledVariantAnnotationsWalker.MODE_LONG_NAME, VariantType.SNP)
                .add(StandardArgumentDefinitions.RESOURCE_LONG_NAME + String.format(":snp,%s=true", LabeledVariantAnnotationsData.SNP_LABEL), INPUT_VCF) // we just use the input VCF as a dummy resource
                .addOutput(outputPrefix);
        final String modelPrefix = new File(INPUT_FROM_TRAIN_EXPECTED_TEST_FILES_DIR,
                "extract.nonAS.snpIndel.posUn.train.snp.posOnly.IF").toString();
        final Function<ArgumentsBuilder, ArgumentsBuilder> addModelPrefix = ab ->
                ADD_MODEL_PREFIX.apply(ab, modelPrefix);
        addModelPrefix.apply(argsBuilder);
        runCommandLine(argsBuilder);
    }
}