package org.broadinstitute.hellbender.engine;

import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCF3Codec;
import htsjdk.variant.vcf.VCFCodec;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentDefinition;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.codecs.table.TableCodec;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

public final class FeatureManagerUnitTest extends GATKBaseTest {
    private static final String FEATURE_MANAGER_TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";

    @DataProvider(name = "DetectCorrectFileFormatTestData")
    public Object[][] getDetectCorrectFileFormatTestData() {
        return new Object[][] {
                { new File(FEATURE_MANAGER_TEST_DIRECTORY + "minimal_vcf4_file.vcf"), VCFCodec.class },
                { new File(FEATURE_MANAGER_TEST_DIRECTORY + "minimal_vcf3_file.vcf"), VCF3Codec.class },
                { new File(FEATURE_MANAGER_TEST_DIRECTORY + "minimal_bcf_file.bcf"), BCF2Codec.class },
                { new File(FEATURE_MANAGER_TEST_DIRECTORY + "minimal_bed_file.bed"), BEDCodec.class},
                { new File(FEATURE_MANAGER_TEST_DIRECTORY + "minimal_table_file.table"), TableCodec.class}
        };
    }

    @Test(dataProvider = "DetectCorrectFileFormatTestData")
    public void testDetectCorrectFileFormat( final File file, final Class<? extends FeatureCodec<? extends Feature, ?>> expectedCodecClass ) throws Exception {
        Assert.assertEquals(FeatureManager.getCodecForFile(file.toPath()).getClass(), expectedCodecClass,
                            "Wrong codec selected for file " + file.getAbsolutePath());

        // We should also get the correct codec if we pass in the explicit expected Feature type to getCodecForFile()
        @SuppressWarnings("unchecked")
        final Class<? extends Feature> expectedCodecFeatureType = expectedCodecClass.getDeclaredConstructor().newInstance().getFeatureType();
        Assert.assertEquals(FeatureManager.getCodecForFile(file.toPath(), expectedCodecFeatureType).getClass(), expectedCodecClass,
                "Wrong codec selected for file " + file.getAbsolutePath() + " after subsetting to the expected Feature type");
    }

    @Test(expectedExceptions = UserException.NoSuitableCodecs.class)
    public void testDetectUnsupportedFileFormat() {
        final File unsupportedFile = new File(FEATURE_MANAGER_TEST_DIRECTORY + "unsupported_format_file");

        // Eliminate file existence/readability as a factor in the test
        Assert.assertTrue(unsupportedFile.canRead(), "Cannot test detection of unsupported file formats on an unreadable file");

        // Should throw, since the file exists and is readable, but is in an unsupported format
        FeatureManager.getCodecForFile(unsupportedFile.toPath());
    }

    @Test(expectedExceptions = UserException.NoSuitableCodecs.class)
    public void testRestrictCodecSelectionToWrongFeatureType() {
        final File vcf = new File(FEATURE_MANAGER_TEST_DIRECTORY + "minimal_vcf4_file.vcf");

        // If we require BED Features from this vcf file, we should get a no suitable codec exception
        FeatureManager.getCodecForFile(vcf.toPath(), BEDFeature.class);
    }

    @DataProvider(name = "IsFeatureFileTestData")
    public Object[][] getIsFeatureFileTestData() {
        return new Object[][] {
                { new File(FEATURE_MANAGER_TEST_DIRECTORY + "minimal_vcf4_file.vcf"), true },
                { new File(FEATURE_MANAGER_TEST_DIRECTORY + "minimal_vcf3_file.vcf"), true },
                { new File(FEATURE_MANAGER_TEST_DIRECTORY + "minimal_bcf_file.bcf"), true },
                { new File(FEATURE_MANAGER_TEST_DIRECTORY + "minimal_bed_file.bed"), true },
                { new File(FEATURE_MANAGER_TEST_DIRECTORY + "minimal_table_file.table"), true },
                { new File(FEATURE_MANAGER_TEST_DIRECTORY + "unsupported_format_file"), false },
                { new File(FEATURE_MANAGER_TEST_DIRECTORY + "reads_data_source_test1.bam"), false },
                { new File(FEATURE_MANAGER_TEST_DIRECTORY + "non_existent_file.bed"), false }
        };
    }

    @Test(dataProvider = "IsFeatureFileTestData")
    public void testIsFeatureFile( final File file, final boolean expectedIsFeatureFile ) {
        Assert.assertEquals(FeatureManager.isFeatureFile(file.toPath()), expectedIsFeatureFile, "isFeatureFile() returned incorrect result for file " + file.getAbsolutePath());
    }

    @CommandLineProgramProperties(summary = "", oneLineSummary = "", programGroup = TestProgramGroup.class)
    private static class ValidFeatureArgumentSource extends CommandLineProgram {
        // We should be able to detect the type parameter of a FeatureInput regardless of whether or
        // not it's wrapped within a Collection

        @Argument(fullName = "variantContextFeatureInput")
        public FeatureInput<VariantContext> variantContextFeatureInput;

        @Argument(fullName = "variantContextListFeatureInput")
        public List<FeatureInput<VariantContext>> variantContextListFeatureInput;

        @Argument(fullName = "bedFeatureInput")
        public FeatureInput<BEDFeature> bedFeatureInput;

        @Argument(fullName = "bedListFeatureInput")
        public List<FeatureInput<BEDFeature>> bedListFeatureInput;

        public ValidFeatureArgumentSource() {
            // The real argument parsing system ensures that Collection argument fields are non-null
            variantContextListFeatureInput = new ArrayList<>();
            bedListFeatureInput = new ArrayList<>();
        }

        @Override
        protected Object doWork() {
            return null;
        }
    }

    @CommandLineProgramProperties(summary = "", oneLineSummary = "",programGroup = TestProgramGroup.class)
    @SuppressWarnings("rawtypes")
    private static class InvalidFeatureArgumentSource extends CommandLineProgram {
        // FeatureInputs without type parameters (ie., raw types) should be detected as errors

        @Argument(fullName = "parameterlessFeatureInput")
        public FeatureInput parameterlessFeatureInput;

        @Argument(fullName = "parameterlessListFeatureInput")
        public List<FeatureInput> parameterlessListFeatureInput;

        @Override
        protected Object doWork() {
            return null;
        }
    }

    @DataProvider(name = "DetectFeatureTypeFromFeatureInputFieldTestData")
    public Object[][] getDetectFeatureTypeFromFeatureInputFieldTestData() {
        return new Object[][] {
                { "variantContextFeatureInput", VariantContext.class },
                { "variantContextListFeatureInput", VariantContext.class },
                { "bedFeatureInput", BEDFeature.class },
                { "bedListFeatureInput", BEDFeature.class }
        };
    }

    @Test(dataProvider = "DetectFeatureTypeFromFeatureInputFieldTestData")
    public void testDetectFeatureTypeFromFeatureInputField( final String fieldName, final Class<?> expectedFeatureType ) throws NoSuchFieldException {
        final ValidFeatureArgumentSource validFeatureSource = new ValidFeatureArgumentSource();
        final CommandLineArgumentParser clp = new CommandLineArgumentParser(validFeatureSource);
        final ArgumentDefinition fieldArgumentDefinition =
                clp.getNamedArgumentDefinitions().stream().filter(argDef -> argDef.getFullName().equals(fieldName)).findFirst().get();
        Assert.assertNotNull(fieldArgumentDefinition);
        Assert.assertEquals(FeatureManager.getFeatureTypeForFeatureInputArgument(fieldArgumentDefinition),
                            expectedFeatureType,
                            "Wrong Feature type detected for field " + fieldName);
    }

    @DataProvider(name = "DetectParameterlessFeatureInputsTestData")
    public Object[][] getDetectParameterlessFeatureInputsTestData() {
        return new Object[][] {
                { "parameterlessFeatureInput" },
                { "parameterlessListFeatureInput" }
        };
    }

    @Test(dataProvider = "DetectParameterlessFeatureInputsTestData", expectedExceptions = GATKException.class)
    public void testDetectParameterlessFeatureInputs( final String fieldName ) throws NoSuchFieldException {
        // FeatureInput fields that lack a type parameter should cause a GATKException to be thrown
        final InvalidFeatureArgumentSource invalidFeatureSource = new InvalidFeatureArgumentSource();
        final CommandLineArgumentParser clp = new CommandLineArgumentParser(invalidFeatureSource);
        final ArgumentDefinition fieldArgumentDefinition =
                clp.getNamedArgumentDefinitions().stream().filter(argDef -> argDef.getFullName().equals(fieldName)).findFirst().get();
        Assert.assertNotNull(fieldArgumentDefinition);
        FeatureManager.getFeatureTypeForFeatureInputArgument(fieldArgumentDefinition);
    }

    @Test
    public void testHandleRequestForValidFeatureInputs() {
        ValidFeatureArgumentSource toolInstance = new ValidFeatureArgumentSource();

        // Initialize two of the FeatureInput fields as they would be initialized by the argument-parsing
        // system to simulate a run of the tool with two FeatureInputs.
        toolInstance.variantContextFeatureInput = new FeatureInput<>(FEATURE_MANAGER_TEST_DIRECTORY + "feature_data_source_test.vcf");
        toolInstance.bedListFeatureInput.add(new FeatureInput<>(FEATURE_MANAGER_TEST_DIRECTORY + "minimal_bed_file.bed"));

        FeatureManager manager = new FeatureManager(toolInstance);
        List<VariantContext> vcFeatures = manager.getFeatures(toolInstance.variantContextFeatureInput, new SimpleInterval("1", 1, 2000));
        List<BEDFeature> bedFeatures = manager.getFeatures(toolInstance.bedListFeatureInput.get(0), new SimpleInterval("1", 1, 1));

        Assert.assertEquals(vcFeatures.size(), 14, "Wrong number of Features returned from VariantContext test Feature file");
        Assert.assertEquals(bedFeatures.size(), 1, "Wrong number of Features returned from BED test Feature file");
    }


    @Test
    public void testGetAllSequenceDictionaries() {
        ValidFeatureArgumentSource toolInstance = new ValidFeatureArgumentSource();

        // Initialize two of the FeatureInput fields as they would be initialized by the argument-parsing
        // system to simulate a run of the tool with two FeatureInputs.
        toolInstance.variantContextFeatureInput = new FeatureInput<>(FEATURE_MANAGER_TEST_DIRECTORY + "feature_data_source_test.vcf");
        toolInstance.bedListFeatureInput.add(new FeatureInput<>(FEATURE_MANAGER_TEST_DIRECTORY + "minimal_bed_file.bed"));

        FeatureManager manager = new FeatureManager(toolInstance);
        final List<SAMSequenceDictionary> dictionaries = manager.getAllSequenceDictionaries();
        Assert.assertEquals(dictionaries.size(), 2);
        Assert.assertEquals(dictionaries.stream().map(dict -> dict.size()).collect(Collectors.toSet()), Sets.newHashSet(1, 4));
    }

    @Test
    public void testSequenceDictionary() {
        final ValidFeatureArgumentSource toolInstance = new ValidFeatureArgumentSource();

        // Initialize the FeatureInput field as it would be initialized by the argument-parsing
        // system.
        toolInstance.variantContextFeatureInput = new FeatureInput<>(FEATURE_MANAGER_TEST_DIRECTORY + "feature_data_source_test_withSequenceDict.vcf");

        final FeatureManager manager = new FeatureManager(toolInstance);
        final List<SAMSequenceDictionary> dicts = manager.getVariantSequenceDictionaries();
        Assert.assertEquals(dicts.size(), 1);
        Assert.assertEquals(dicts.get(0).getSequences().size(), 84);
    }

    @Test
    public void testTwoSequenceDictionaries() {
        final ValidFeatureArgumentSource toolInstance = new ValidFeatureArgumentSource();

        // Initialize two of the FeatureInput fields as they would be initialized by the argument-parsing
        // system to simulate a run of the tool with two FeatureInputs.
        toolInstance.variantContextFeatureInput = new FeatureInput<>(FEATURE_MANAGER_TEST_DIRECTORY + "feature_data_source_test_withSequenceDict.vcf");
        toolInstance.variantContextListFeatureInput.add(new FeatureInput<>(FEATURE_MANAGER_TEST_DIRECTORY + "feature_data_source_test_with_bigHeader.vcf"));

        final FeatureManager manager = new FeatureManager(toolInstance);
        final List<SAMSequenceDictionary> dicts = manager.getVariantSequenceDictionaries();
        Assert.assertEquals(dicts.size(), 2);
    }

    @Test
    public void testEmptySequenceDictionary() {
        final ValidFeatureArgumentSource toolInstance = new ValidFeatureArgumentSource();

        // Initialize the FeatureInput field as it would be initialized by the argument-parsing
        // system.
        toolInstance.variantContextFeatureInput = new FeatureInput<>(FEATURE_MANAGER_TEST_DIRECTORY + "feature_data_source_test.vcf");

        final FeatureManager manager = new FeatureManager(toolInstance);
        final List<SAMSequenceDictionary> dicts = manager.getVariantSequenceDictionaries();
        Assert.assertTrue(dicts.isEmpty(), Objects.toString(dicts));
    }

    @Test
    public void testHandleRequestForValidFeatureInputIterator() {
        final ValidFeatureArgumentSource toolInstance = new ValidFeatureArgumentSource();

        // Initialize two of the FeatureInput fields as they would be initialized by the argument-parsing
        // system to simulate a run of the tool with two FeatureInputs.
        toolInstance.variantContextFeatureInput = new FeatureInput<>(FEATURE_MANAGER_TEST_DIRECTORY + "feature_data_source_test.vcf");
        toolInstance.bedListFeatureInput.add(new FeatureInput<>(FEATURE_MANAGER_TEST_DIRECTORY + "minimal_bed_file.bed"));

        final FeatureManager manager = new FeatureManager(toolInstance);
        final Iterator<VariantContext> vcIterator = manager.getFeatureIterator(toolInstance.variantContextFeatureInput);
        final Iterator<BEDFeature> bedIterator = manager.getFeatureIterator(toolInstance.bedListFeatureInput.get(0));

        final List<VariantContext> variants = Utils.stream(vcIterator).collect(Collectors.toList());
        final List<BEDFeature> beds = Utils.stream(bedIterator).collect(Collectors.toList());

        Assert.assertEquals(variants.size(), 26);
        Assert.assertEquals(variants.get(4).getStart(),280);
        Assert.assertEquals(beds.size(), 1);
    }

    @Test
    public void testHandleRequestForValidFeatureInputIteratorWithoutIndex() {
        final ValidFeatureArgumentSource toolInstance = new ValidFeatureArgumentSource();

        // Initialize a FeatureInput field as it would be initialized by the argument-parsing
        // system to simulate a run of the tool with that FeatureInput.
        toolInstance.variantContextFeatureInput = new FeatureInput<>(FEATURE_MANAGER_TEST_DIRECTORY + "feature_data_source_test.wo-idx.vcf");

        final FeatureManager manager = new FeatureManager(toolInstance);
        final Iterator<VariantContext> vcIterator = manager.getFeatureIterator(toolInstance.variantContextFeatureInput);

        final List<VariantContext> variants = Utils.stream(vcIterator).collect(Collectors.toList());

        Assert.assertEquals(variants.size(), 26);
        Assert.assertEquals(variants.get(4).getStart(),280);
    }

    @Test(expectedExceptions = GATKException.class)
    public void testHandleRequestForNonExistentFeatureInput() {
        FeatureManager manager = new FeatureManager(new ValidFeatureArgumentSource());

        // Requests for FeatureInputs not declared in the tool's class hierarchy (or associated ArgumentCollections)
        // should throw an exception
        manager.getFeatures(new FeatureInput<>("featureInputNotDeclaredInTool"), new SimpleInterval("1", 1, 1));
    }
}
