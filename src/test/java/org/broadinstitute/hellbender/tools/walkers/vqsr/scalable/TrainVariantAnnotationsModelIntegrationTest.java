package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import com.google.common.collect.Lists;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.LabeledVariantAnnotationsData;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.VariantType;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.BGMMVariantAnnotationsScorer;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.PythonVariantAnnotationsScorer;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.VariantAnnotationsModelBackend;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collectors;

/**
 * See documentation for {@link ExtractVariantAnnotationsIntegrationTest} for information about how inputs and
 * expected outputs used there are related to those used here and in {@link ScoreVariantAnnotationsIntegrationTest}.
 */
public final class TrainVariantAnnotationsModelIntegrationTest extends CommandLineProgramTest {

    // If true, update the expected outputs in tests that assert an exact match vs. prior output,
    // instead of actually running the tests. Can be used with "./gradlew test -Dtest.single=TrainVariantAnnotationsIntegrationTest"
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

    private static final File TEST_FILES_DIR = new File(largeFileTestDir,
            "org/broadinstitute/hellbender/tools/walkers/vqsr/scalable/train");
    private static final File INPUT_FROM_EXTRACT_EXPECTED_TEST_FILES_DIR = new File(largeFileTestDir,
            "org/broadinstitute/hellbender/tools/walkers/vqsr/scalable/extract/expected");
    private static final File EXPECTED_TEST_FILES_DIR = new File(TEST_FILES_DIR, "expected");

    private static final File ISOLATION_FOREST_PYTHON_SCRIPT = IOUtils.writeTempResource(
            new Resource("isolation-forest.py", TrainVariantAnnotationsModel.class));
    private static final File ISOLATION_FOREST_HYPERPARAMETERS_JSON = new File(TEST_FILES_DIR,
            "isolation-forest-hyperparameters-different-seed.json");

    // Supplier and functions for creating and adding various arguments to an ArgumentsBuilder.
    private static final Supplier<ArgumentsBuilder> BASE_ARGUMENTS_BUILDER_SUPPLIER = ArgumentsBuilder::new;
    private static final BiFunction<ArgumentsBuilder, File, ArgumentsBuilder> ADD_ANNOTATIONS_HDF5 = (argsBuilder, annotationsHDF5) -> {
        argsBuilder.add(TrainVariantAnnotationsModel.ANNOTATIONS_HDF5_LONG_NAME, annotationsHDF5);
        return argsBuilder;
    };
    private static final BiFunction<ArgumentsBuilder, File, ArgumentsBuilder> ADD_UNLABELED_ANNOTATIONS_HDF5 = (argsBuilder, unlabeledAnnotationsHDF5) -> {
        argsBuilder.add(TrainVariantAnnotationsModel.UNLABELED_ANNOTATIONS_HDF5_LONG_NAME, unlabeledAnnotationsHDF5);
        return argsBuilder;
    };
    private static final Function<ArgumentsBuilder, ArgumentsBuilder> ADD_SNP_MODE = argsBuilder -> {
        argsBuilder.add(LabeledVariantAnnotationsWalker.MODE_LONG_NAME, VariantType.SNP);
        return argsBuilder;
    };
    private static final Function<ArgumentsBuilder, ArgumentsBuilder> ADD_INDEL_MODE = argsBuilder -> {
        argsBuilder.add(LabeledVariantAnnotationsWalker.MODE_LONG_NAME, VariantType.INDEL);
        return argsBuilder;
    };
    private static final BiFunction<ArgumentsBuilder, VariantAnnotationsModelBackend, ArgumentsBuilder> ADD_MODEL_BACKEND = (argsBuilder, modelBackendMode) -> {
        argsBuilder.add(TrainVariantAnnotationsModel.MODEL_BACKEND_LONG_NAME, modelBackendMode);
        return argsBuilder;
    };
    private static final Function<ArgumentsBuilder, ArgumentsBuilder> ADD_ISOLATION_FOREST_PYTHON_SCRIPT = argsBuilder -> {
        argsBuilder.add(TrainVariantAnnotationsModel.PYTHON_SCRIPT_LONG_NAME, ISOLATION_FOREST_PYTHON_SCRIPT);
        return argsBuilder;
    };
    private static final Function<ArgumentsBuilder, ArgumentsBuilder> ADD_ISOLATION_FOREST_HYPERPARAMETERS_JSON = argsBuilder -> {
        argsBuilder.add(TrainVariantAnnotationsModel.HYPERPARAMETERS_JSON_LONG_NAME, ISOLATION_FOREST_HYPERPARAMETERS_JSON);
        return argsBuilder;
    };

    /**
     * Exact-match tests for (non-exhaustive) configurations given by the Cartesian product of the following options:
     *  1) non-allele-specific ("nonAS") vs. allele-specific ("AS")
     *  2) SNP-only ("snp") vs. SNP+INDEL ("snpIndel") (for both of these options, we use extracted annotations that contain both SNP and INDEL variants as input)
     *  3) model backend
     *      3a) Java Bayesian Gaussian Mixture Model (BGMM) backend TODO the BGMM has been reduced to a stub for this initial PR; subsequent PRs will cover the backend code and reconnect the stub
     *      3b) default PYTHON_IFOREST with default hyperparameters ("IF")
     *      3c) default PYTHON_IFOREST with non-default seed hyperparameter ("IFDifferentSeed")
     *      3d) specified PYTHON_SCRIPT with non-default seed hyperparameter ("IFDifferentSeed"); we will simply use the same script as the default PYTHON_IFOREST backend, so this is just a test of the command-line interface
     *      We should expect 3c-d to give functionally identical results.
     */
    @DataProvider(name = "dataValidInputs")
    public Object[][] dataValidInputs() {
        final List<List<Pair<String, Function<ArgumentsBuilder, ArgumentsBuilder>>>> testConfigurations = Lists.cartesianProduct(
                Arrays.asList(
                        Pair.of("extract.nonAS.snpIndel.posUn.train", Function.identity()),
                        Pair.of("extract.AS.snpIndel.posUn.train", Function.identity())),
                Arrays.asList(
                        Pair.of("snp", ADD_SNP_MODE),
                        Pair.of("snpIndel", ADD_SNP_MODE.andThen(ADD_INDEL_MODE))),
                Collections.singletonList(
                        Pair.of("posOnly", Function.identity())),
                Arrays.asList(
                        Pair.of("IF", ab -> ADD_MODEL_BACKEND.apply(ab, VariantAnnotationsModelBackend.PYTHON_IFOREST)),
                        Pair.of("IFDifferentSeed", ADD_ISOLATION_FOREST_HYPERPARAMETERS_JSON
                                .andThen(ab -> ADD_MODEL_BACKEND.apply(ab, VariantAnnotationsModelBackend.PYTHON_IFOREST))), // this and the following case give the same results, so they are given the same IFDifferentSeed tag
                        Pair.of("IFDifferentSeed", ADD_ISOLATION_FOREST_PYTHON_SCRIPT
                                .andThen(ADD_ISOLATION_FOREST_HYPERPARAMETERS_JSON)
                                .andThen(ab -> ADD_MODEL_BACKEND.apply(ab, VariantAnnotationsModelBackend.PYTHON_SCRIPT)))));

        return testConfigurations.stream()
                .map(tagAndAddFunctionPairs -> new Object[]{
                        tagAndAddFunctionPairs.stream().map(Pair::getLeft).collect(Collectors.joining(".")), // e.g., extract.nonAS.snpIndel.posUn.train.snp.posOnly.IF
                        tagAndAddFunctionPairs.stream().map(Pair::getRight)                                              // creates the corresponding ArgumentsBuilder
                                .reduce(Function.identity(), Function::andThen)                                          //  by stringing together functions that add the
                                .apply(BASE_ARGUMENTS_BUILDER_SUPPLIER.get())})                                          //  appropriate arguments
                .toArray(Object[][]::new);
    }

    /**
     * Checks expected outputs given a tag (e.g., "extract.nonAS.snpIndel.posUn.train.snp.posOnly.IF") and arguments corresponding to the
     * Cartesian products generated in {@link #dataValidInputs}.
     *
     * We perform exact-match tests of any HDF5 files produced using h5diff, which is insensitive to timestamps within the file.
     * Binary serialized scorers may not be diff equivalent, so we just check for their existence.
     */
    @Test(dataProvider = "dataValidInputs", groups = {"python"}) // python environment is required to run tool and to use h5diff for exact-match comparisons
    public void testValidInputs(final String tag,
                                final ArgumentsBuilder argsBuilder) {
        final File outputDir = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? EXPECTED_TEST_FILES_DIR : createTempDir("train");
        final String outputPrefix = String.format("%s/%s", outputDir, tag);
        argsBuilder.addOutput(outputPrefix);

        // add arguments for positive/unlabeled annotations based on the
        // extract tag (the portion of the tag preceding ".train", e.g., extract.nonAS.snpIndel.posUn),
        // which gives the basename for the annotation files
        final String extractTag = tag.split(".train")[0];
        final File positiveAnnotationsHDF5 = new File(INPUT_FROM_EXTRACT_EXPECTED_TEST_FILES_DIR,
                extractTag + LabeledVariantAnnotationsWalker.ANNOTATIONS_HDF5_SUFFIX);
        final Function<ArgumentsBuilder, ArgumentsBuilder> addPositiveAnnotations = ab ->
                ADD_ANNOTATIONS_HDF5.apply(ab, positiveAnnotationsHDF5);
        addPositiveAnnotations.apply(argsBuilder);

        runCommandLine(argsBuilder);

        if (!UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS) {
            assertExpectedOutputs(tag, outputPrefix);
        }
    }

    private static void assertExpectedOutputs(final String tag,
                                              final String outputPrefix) {
        if (tag.contains("train.snp.")) {
            assertExpectedOutputsForVariantType(tag, outputPrefix, "snp");
            assertOutputsForVariantTypeDoNotExist(outputPrefix, "indel");
        } else if (tag.contains("train.snpIndel.")) {
            assertExpectedOutputsForVariantType(tag, outputPrefix, "snp");
            assertExpectedOutputsForVariantType(tag, outputPrefix, "indel");
        } else {
            Assert.fail("Unknown variant-type tag.");
        }
    }

    private static void assertExpectedOutputsForVariantType(final String tag,
                                                            final String outputPrefix,
                                                            final String variantType) {
        final String tagAndVariantType = String.format("%s.%s", tag, variantType);
        final String outputPrefixAndVariantType = String.format("%s.%s", outputPrefix, variantType);

        SystemCommandUtilsTest.runH5Diff(
                String.format("%s/%s", EXPECTED_TEST_FILES_DIR, tagAndVariantType + TrainVariantAnnotationsModel.TRAINING_SCORES_HDF5_SUFFIX),
                String.format("%s", outputPrefixAndVariantType + TrainVariantAnnotationsModel.TRAINING_SCORES_HDF5_SUFFIX));

        SystemCommandUtilsTest.runH5Diff(
                String.format("%s/%s", EXPECTED_TEST_FILES_DIR, tagAndVariantType + TrainVariantAnnotationsModel.CALIBRATION_SCORES_HDF5_SUFFIX),
                String.format("%s", outputPrefixAndVariantType + TrainVariantAnnotationsModel.CALIBRATION_SCORES_HDF5_SUFFIX));

        assertScorerExpectedOutputs(tagAndVariantType, outputPrefixAndVariantType);

        if (tag.contains("posOnly")) {
            Assert.assertFalse(new File(outputPrefixAndVariantType + TrainVariantAnnotationsModel.UNLABELED_SCORES_HDF5_SUFFIX).exists());
        }
    }

    private static void assertOutputsForVariantTypeDoNotExist(final String outputPrefix,
                                                              final String variantType) {
        final String outputPrefixAndVariantType = String.format("%s.%s", outputPrefix, variantType);
        Assert.assertFalse(new File(outputPrefixAndVariantType + TrainVariantAnnotationsModel.TRAINING_SCORES_HDF5_SUFFIX).exists());
        Assert.assertFalse(new File(outputPrefixAndVariantType + TrainVariantAnnotationsModel.CALIBRATION_SCORES_HDF5_SUFFIX).exists());
        Assert.assertFalse(new File(outputPrefixAndVariantType + TrainVariantAnnotationsModel.UNLABELED_SCORES_HDF5_SUFFIX).exists());
        Assert.assertFalse(new File(outputPrefixAndVariantType + BGMMVariantAnnotationsScorer.BGMM_SCORER_SER_SUFFIX).exists());
        Assert.assertFalse(new File(outputPrefixAndVariantType + PythonVariantAnnotationsScorer.PYTHON_SCORER_PKL_SUFFIX).exists());
    }

    /**
     * Binary serialized scorers may not be diff equivalent, so we just check for their existence.
     * We assume that checking elsewhere for equivalence of the scores that the scorers generate provides sufficient
     * coverage.
     */
    private static void assertScorerExpectedOutputs(final String tagAndVariantType,
                                                    final String outputPrefixAndVariantType) {
        final String scorerTag = outputPrefixAndVariantType;
        if (tagAndVariantType.contains("BGMM")) {
            Assert.assertTrue(new File(scorerTag + BGMMVariantAnnotationsScorer.BGMM_SCORER_SER_SUFFIX).exists());
            Assert.assertFalse(new File(scorerTag + PythonVariantAnnotationsScorer.PYTHON_SCORER_PKL_SUFFIX).exists());
        } else if (tagAndVariantType.contains("IF")) {
            Assert.assertTrue(new File(scorerTag + PythonVariantAnnotationsScorer.PYTHON_SCORER_PKL_SUFFIX).exists());
            Assert.assertFalse(new File(scorerTag + BGMMVariantAnnotationsScorer.BGMM_SCORER_SER_SUFFIX).exists());
        } else {
            Assert.fail("Unknown model-backend tag.");
        }
    }

    @Test(groups = {"python"}) // python environment is required to run tool and to use h5diff for exact-match comparisons
    public void testSNPOnlyModelsFromSNPOnlyAndSNPPlusIndelAnnotationsAreIdentical() {
        final File outputDir = createTempDir("train");

        final String outputPrefixSNPOnly = String.format("%s/test-snp", outputDir);
        final ArgumentsBuilder argsBuilderSNPOnly = BASE_ARGUMENTS_BUILDER_SUPPLIER.get();
        argsBuilderSNPOnly.addOutput(outputPrefixSNPOnly);
        final File positiveAnnotationsHDF5SNPOnly = new File(INPUT_FROM_EXTRACT_EXPECTED_TEST_FILES_DIR,
                "extract.nonAS.snp.pos" + LabeledVariantAnnotationsWalker.ANNOTATIONS_HDF5_SUFFIX);
        final Function<ArgumentsBuilder, ArgumentsBuilder> addPositiveAnnotationsSNPOnly = ab ->
                ADD_ANNOTATIONS_HDF5.apply(ab, positiveAnnotationsHDF5SNPOnly);
        addPositiveAnnotationsSNPOnly
                .andThen(ADD_SNP_MODE)
                .apply(argsBuilderSNPOnly);
        runCommandLine(argsBuilderSNPOnly);

        final String outputPrefixSNPPlusIndel = String.format("%s/test-snpIndel", outputDir);
        final ArgumentsBuilder argsBuilderSNPPlusIndel = BASE_ARGUMENTS_BUILDER_SUPPLIER.get();
        argsBuilderSNPPlusIndel.addOutput(outputPrefixSNPPlusIndel);
        final File positiveAnnotationsHDF5SNPPlusIndel = new File(INPUT_FROM_EXTRACT_EXPECTED_TEST_FILES_DIR,
                "extract.nonAS.snpIndel.pos" + LabeledVariantAnnotationsWalker.ANNOTATIONS_HDF5_SUFFIX);
        final Function<ArgumentsBuilder, ArgumentsBuilder> addPositiveAnnotationsSNPPlusIndel = ab ->
                ADD_ANNOTATIONS_HDF5.apply(ab, positiveAnnotationsHDF5SNPPlusIndel);
        addPositiveAnnotationsSNPPlusIndel
                .andThen(ADD_SNP_MODE)
                .apply(argsBuilderSNPPlusIndel);
        runCommandLine(argsBuilderSNPPlusIndel);

        SystemCommandUtilsTest.runH5Diff(
                outputPrefixSNPOnly + ".snp" + TrainVariantAnnotationsModel.TRAINING_SCORES_HDF5_SUFFIX,
                outputPrefixSNPPlusIndel + ".snp" + TrainVariantAnnotationsModel.TRAINING_SCORES_HDF5_SUFFIX);

        SystemCommandUtilsTest.runH5Diff(
                outputPrefixSNPOnly + ".snp" + TrainVariantAnnotationsModel.CALIBRATION_SCORES_HDF5_SUFFIX,
                outputPrefixSNPPlusIndel + ".snp" + TrainVariantAnnotationsModel.CALIBRATION_SCORES_HDF5_SUFFIX);
    }

    @Test(expectedExceptions = IllegalArgumentException.class) // python environment is required to run tool
    public void testPositiveAndUnlabeledAnnotationNamesAreNotIdentical() {
        final File outputDir = createTempDir("train");
        final String outputPrefix = String.format("%s/test", outputDir);
        final ArgumentsBuilder argsBuilder = BASE_ARGUMENTS_BUILDER_SUPPLIER.get();
        argsBuilder.addOutput(outputPrefix);
        final File positiveAnnotationsHDF5 = new File(INPUT_FROM_EXTRACT_EXPECTED_TEST_FILES_DIR,
                "extract.nonAS.snpIndel.posUn" + LabeledVariantAnnotationsWalker.ANNOTATIONS_HDF5_SUFFIX);                                          // non-allele-specific
        final Function<ArgumentsBuilder, ArgumentsBuilder> addPositiveAnnotations = ab ->
                ADD_ANNOTATIONS_HDF5.apply(ab, positiveAnnotationsHDF5);
        final File unlabeledAnnotationsHDF5 = new File(INPUT_FROM_EXTRACT_EXPECTED_TEST_FILES_DIR,
                "extract.AS.snpIndel.posUn" + ExtractVariantAnnotations.UNLABELED_TAG + LabeledVariantAnnotationsWalker.ANNOTATIONS_HDF5_SUFFIX);  // allele-specific
        final Function<ArgumentsBuilder, ArgumentsBuilder> addUnlabeledAnnotations = ab ->
                ADD_UNLABELED_ANNOTATIONS_HDF5.apply(ab, unlabeledAnnotationsHDF5);
        addPositiveAnnotations
                .andThen(addUnlabeledAnnotations)
                .apply(argsBuilder);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = UserException.BadInput.class, groups = {"python"}) // python environment is required to run tool
    public void testPositiveAnnotationsOfSpecifiedVariantTypesNotPresent() {
        final File outputDir = createTempDir("train");
        final String outputPrefix = String.format("%s/test", outputDir);
        final ArgumentsBuilder argsBuilder = BASE_ARGUMENTS_BUILDER_SUPPLIER.get();
        argsBuilder.addOutput(outputPrefix);
        final File positiveAnnotationsHDF5 = new File(INPUT_FROM_EXTRACT_EXPECTED_TEST_FILES_DIR,
                "extract.nonAS.snp.posUn" + LabeledVariantAnnotationsWalker.ANNOTATIONS_HDF5_SUFFIX);     // contains only SNPs, but SNP+INDEL is specified
        final Function<ArgumentsBuilder, ArgumentsBuilder> addPositiveAnnotations = ab ->
                ADD_ANNOTATIONS_HDF5.apply(ab, positiveAnnotationsHDF5);
        ADD_SNP_MODE
                .andThen(ADD_INDEL_MODE)
                .andThen(addPositiveAnnotations)
                .apply(argsBuilder);
        runCommandLine(argsBuilder);
    }

    // we will enable this once a positive-unlabeled backend is implemented
    @Test(expectedExceptions = UserException.BadInput.class, groups = {"python"}, enabled = false) // python environment is required to run tool
    public void testUnlabeledAnnotationsOfSpecifiedVariantTypesNotPresent() {
        final File outputDir = createTempDir("train");
        final String outputPrefix = String.format("%s/test", outputDir);
        final ArgumentsBuilder argsBuilder = BASE_ARGUMENTS_BUILDER_SUPPLIER.get();
        argsBuilder.addOutput(outputPrefix);
        final File positiveAnnotationsHDF5 = new File(INPUT_FROM_EXTRACT_EXPECTED_TEST_FILES_DIR,
                "extract.nonAS.snpIndel.posUn" + LabeledVariantAnnotationsWalker.ANNOTATIONS_HDF5_SUFFIX);
        final Function<ArgumentsBuilder, ArgumentsBuilder> addPositiveAnnotations = ab ->
                ADD_ANNOTATIONS_HDF5.apply(ab, positiveAnnotationsHDF5);
        final File unlabeledAnnotationsHDF5 = new File(INPUT_FROM_EXTRACT_EXPECTED_TEST_FILES_DIR,
                "extract.nonAS.snp.posUn" + ExtractVariantAnnotations.UNLABELED_TAG + LabeledVariantAnnotationsWalker.ANNOTATIONS_HDF5_SUFFIX);    // contains only SNPs, but SNP+INDEL is specified
        final Function<ArgumentsBuilder, ArgumentsBuilder> addUnlabeledAnnotations = ab ->
                ADD_UNLABELED_ANNOTATIONS_HDF5.apply(ab, unlabeledAnnotationsHDF5);
        ADD_SNP_MODE.andThen(ADD_INDEL_MODE)
                .andThen(addPositiveAnnotations)
                .andThen(addUnlabeledAnnotations)
                .apply(argsBuilder);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = UserException.BadInput.class, groups = {"python"}) // python environment is required to run tool
    public void testPositiveAnnotationForOneVariantTypeIsCompletelyMissing() { // TODO add analogous test that warning is emitted when annotation has zero variance?
        final File outputDir = createTempDir("train");
        final String outputPrefix = String.format("%s/test", outputDir);
        final ArgumentsBuilder argsBuilder = BASE_ARGUMENTS_BUILDER_SUPPLIER.get();
        argsBuilder.addOutput(outputPrefix);

        // we will dummy up an annotations file that contains 2 annotations (ANNOT_1 and ANNOT_2)
        // for 4 variants (2 SNPs and 2 INDELs); the INDELs will all have missing (i.e., NaN) ANNOT_1 values
        final List<String> annotationNames = Arrays.asList("ANNOT_1", "ANNOT_2");
        final double[][] annotations = new double[][]{
                new double[]{1, 2},             // SNP
                new double[]{3, 4},             // SNP
                new double[]{Double.NaN, 2},    // INDEL
                new double[]{Double.NaN, 4}};   // INDEL
        final List<Boolean> isSubset = Collections.nCopies(4, true);

        final File positiveAnnotationsHDF5 = LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile(
                annotationNames, annotations, isSubset);

        try (final HDF5File positiveAnnotationsHDF5File = new HDF5File(positiveAnnotationsHDF5, HDF5File.OpenMode.READ_WRITE)) {
            positiveAnnotationsHDF5File.makeDoubleArray("/labels/snp", new double[]{1, 1, 0, 0});
            positiveAnnotationsHDF5File.makeDoubleArray("/labels/training", new double[]{1, 1, 1, 1});
            positiveAnnotationsHDF5File.makeDoubleArray("/labels/calibration", new double[]{1, 1, 1, 1});
        }
        final Function<ArgumentsBuilder, ArgumentsBuilder> addPositiveAnnotations = ab ->
                ADD_ANNOTATIONS_HDF5.apply(ab, positiveAnnotationsHDF5);

        ADD_SNP_MODE.andThen(ADD_INDEL_MODE)
                .andThen(addPositiveAnnotations)
                .apply(argsBuilder);
        runCommandLine(argsBuilder);
    }
}