package org.broadinstitute.hellbender.engine;

import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.SparkTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.util.Collections;

public final class FeatureInputUnitTest extends GATKBaseTest {
    private static final String FEATURE_INPUT_TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";

    static class ArgumentContainer {
        @Argument(shortName="argName")
        public FeatureInput<Feature> fi;
    }

    private FeatureInput<Feature> runCommandLineWithTaggedFeatureInput(final String taggedFeatureArgument, final String argumentValue) {
        ArgumentContainer ac = new ArgumentContainer();
        CommandLineArgumentParser clp = new CommandLineArgumentParser(ac);
        final String[] args = {"--" + taggedFeatureArgument, argumentValue};
        clp.parseArguments(System.out, args);
        return ac.fi;
    }

    @DataProvider(name = "InvalidFeatureTagsDataProvider")
    public Object[][] getInvalidFeatureTags() {
        return new Object[][] {
                //{ "name:file:file2" },            //this is legal (argument has tag name "file:file2")
                { "name:" },
                { ":file" },
                { ",:file" },
                { "name,key=value=fred:file" },
                { "name,:file"},
                { ",key=value:file"},
                {"name,key:file"},
                {"name,key=:file"},
                {"name,=value:file"},
                {"name,=:file"},
                { ":" },
                { ",:" },
                { "::" },
                { "" },
                { "name,key=value1,key=value2:file" },   //duplicate key
                { "name,key=value,key=value:file" },      //duplicate key
                //{ "name:name:gendb://mydb" }            //this is legal (argument has tag name "name:gendb://mydb")
        };
    }

    @Test(dataProvider = "InvalidFeatureTagsDataProvider", expectedExceptions = CommandLineException.class)
    public void testInvalidFeatureTags( final String invalidFeatureArgument ) {
        runCommandLineWithTaggedFeatureInput(invalidFeatureArgument, "value");
    }

    @DataProvider(name = "ValidFileOnlyFeatureArgumentValuesDataProvider")
    public Object[][] getValidFileOnlyFeatureArgumentValues() {
        return new Object[][] {
                {"myFile"},
                {"myName,key1=value,myFile"},     //allowed - all of this is treated as a file name
                {"=myFile"},                      //allowed - all of this is treated as a file name
                {",myFile"},                    //allowed - all of this is treated as a file name
                {"=,myFile"},                    //allowed - all of this is treated as a file name
                {"key1=value,myFile"}             //allowed - all of this is treated as a file name
        };
    }

    @Test(dataProvider = "ValidFileOnlyFeatureArgumentValuesDataProvider")
    public void testNoFeatureNameSpecified(final String validFileOnlyFeatureArgumentValue) {
        FeatureInput<Feature> featureInput = runCommandLineWithTaggedFeatureInput("argName", validFileOnlyFeatureArgumentValue);

        Assert.assertEquals(featureInput.getFeaturePath(), validFileOnlyFeatureArgumentValue, "Wrong File in FeatureInput");
        // Name should default to the absolute path of the File when no name is specified
        Assert.assertEquals(featureInput.getName(), new File(validFileOnlyFeatureArgumentValue).getAbsolutePath(), "Wrong default name in FeatureInput");
    }

    @DataProvider(name = "GenDbPathAndNameData")
    public Object[][] genDbPathAndNameData() {
        return new Object[][] {
                // input arg name, input value, expected Feature path, expected logical name
                {"argName", "gendb://myJsons", "gendb://myJsons", "gendb://" + new File("myJsons").getAbsolutePath()},
                {"argName:myname", "gendb://myJsons", "gendb://myJsons", "myname"},
                {"argName:myname,key1=value1", "gendb://myJsons", "gendb://myJsons", "myname"},
                {"argName:myname//", "gendb://myJsons", "gendb://myJsons", "myname//"},

                {"argName", "gendb.gs://myBucket/myJsons", "gendb.gs://myBucket/myJsons", "gendb.gs://myBucket/myJsons"},
                {"argName:myname", "gendb.gs://myJsons", "gendb.gs://myJsons", "myname"},
                {"argName:myname,key1=value1", "gendb.gs://myJsons", "gendb.gs://myJsons", "myname"},
                {"argName:myname//", "gendb.gs://myJsons", "gendb.gs://myJsons", "myname//"},

                {"argName", "gendb.hdfs://localhost/myJsons", "gendb.hdfs://localhost/myJsons", "gendb.hdfs://localhost/myJsons"},
                {"argName:myname", "gendb.hdfs://myJsons", "gendb.hdfs://myJsons", "myname"},
                {"argName:myname,key1=value1", "gendb.hdfs://myJsons", "gendb.hdfs://myJsons", "myname"},
                {"argName:myname//", "gendb.hdfs://myJsons", "gendb.hdfs://myJsons", "myname//"},
        };
    }
    
    @DataProvider(name = "GcsPathAndNameData")
    public Object[][] gcsPathAndNameData() {
        return new Object[][] {
                // input arg name, input value, expected Feature path, expected logical name
                {"argName", "gs://bucket/user/my.vcf", "gs://bucket/user/my.vcf", "gs://bucket/user/my.vcf"},
                {"argName:myname", "gs://bucket/user/my.vcf", "gs://bucket/user/my.vcf", "myname"},
                {"argName:myname,key1=value1", "gs://bucket/user/my.vcf", "gs://bucket/user/my.vcf", "myname"},
                {"argName:myname//", "gs://bucket/user/my.vcf", "gs://bucket/user/my.vcf", "myname//"}
        };
    }
    
    @DataProvider(name = "HdfsPathAndNameData")
    public Object[][] hdfsPathAndNameData() {
        return new Object[][] {
                // input arg name, input value, expected Feature path, expected logical name
                {"argName", "hdfs://localhost/user/my.vcf", "hdfs://localhost/user/my.vcf", "hdfs://localhost:8020/user/my.vcf"},
                {"argName:myname", "hdfs://localhost/user/my.vcf", "hdfs://localhost/user/my.vcf", "myname"},
                {"argName:myname,key1=value1", "hdfs://localhost/user/my.vcf", "hdfs://localhost/user/my.vcf", "myname"},
                {"argName:myname//", "hdfs://localhost/user/my.vcf", "hdfs://localhost/user/my.vcf", "myname//"}
        };
    }

    @Test(dataProvider = "GenDbPathAndNameData")
    public void testGenDbPathAndName( final String argWithTags, final String inputValue, final String expectedFeaturePath, final String expectedLogicalName ) {
        FeatureInput<Feature> gendbInput = runCommandLineWithTaggedFeatureInput(argWithTags, inputValue);

        Assert.assertEquals(gendbInput.getFeaturePath(), expectedFeaturePath, "wrong featurePath");
        Assert.assertEquals(gendbInput.getName(), expectedLogicalName, "wrong logical name");
    }

    @Test(dataProvider = "GcsPathAndNameData", groups={"bucket"})
    public void testGcsPathAndName( final String argWithTags, final String inputValue, final String expectedFeaturePath, final String expectedLogicalName ) {
        final FeatureInput<Feature> gcsInput = runCommandLineWithTaggedFeatureInput(argWithTags, inputValue);

        Assert.assertEquals(gcsInput.getFeaturePath(), expectedFeaturePath, "wrong featurePath");
        Assert.assertEquals(gcsInput.getName(), expectedLogicalName, "wrong logical name");
    }

    @Test(dataProvider = "HdfsPathAndNameData")
    public void testHdfsPathAndName( final String argWithTags, final String inputValue, final String expectedFeaturePath, final String expectedLogicalName ) {
        final FeatureInput<Feature> hdfsInput = runCommandLineWithTaggedFeatureInput(argWithTags, inputValue);

        Assert.assertEquals(hdfsInput.getFeaturePath(), expectedFeaturePath, "wrong featurePath");
        Assert.assertEquals(hdfsInput.getName(), expectedLogicalName, "wrong logical name");
    }

    @Test
    public void testFeatureNameSpecified() {
        final FeatureInput<Feature> featureInput = runCommandLineWithTaggedFeatureInput("argName:myName", "myFile");

        Assert.assertEquals(featureInput.getFeaturePath(), "myFile", "Wrong File in FeatureInput");
        Assert.assertEquals(featureInput.getName(), "myName", "Wrong name in FeatureInput");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testRejectNullFile() {
        new FeatureInput<>((GATKPath) null, "sourceName1");
    }

    @Test
    public void testNullOKAsFeatureName() {
        final FeatureInput<Feature> featureInput = runCommandLineWithTaggedFeatureInput("argName:null", "myFile");

        Assert.assertEquals(featureInput.getFeaturePath(), "myFile", "Wrong File in FeatureInput");
        Assert.assertEquals(featureInput.getName(), "null", "Wrong name in FeatureInput");
    }

    @Test
    public void testNullOKAsFileName() {
        final FeatureInput<Feature> featureInput = runCommandLineWithTaggedFeatureInput("argName:myName", "null");

        Assert.assertEquals(featureInput.getFeaturePath(), "null", "Wrong File in FeatureInput");
        Assert.assertEquals(featureInput.getName(), "myName", "Wrong name in FeatureInput");
    }


    @Test
    public void testFeatureKeyValuePairsSpecified() {
        final FeatureInput<Feature> featureInput = runCommandLineWithTaggedFeatureInput("argName:myName,key1=value1,key2=value2,null=null", "myFile");

        Assert.assertEquals(featureInput.getAttribute("key1"), "value1", "wrong attribute value for key1");
        Assert.assertEquals(featureInput.getAttribute("key2"), "value2", "wrong attribute value for key2");
        Assert.assertEquals(featureInput.getAttribute("null"), "null", "wrong attribute value for key \"null\"");
        Assert.assertEquals(featureInput.getAttribute("key3"), null, "wrong attribute value for key3 (not present)");

        Assert.assertEquals(featureInput.getName(), "myName");
        Assert.assertEquals(featureInput.getFeaturePath(), "myFile");
    }

    @Test
    public void testFeatureKeyValuePairSpecified() {
        final FeatureInput<Feature> featureInput = runCommandLineWithTaggedFeatureInput("argName:myName,key1=value1", "myFile");

        Assert.assertEquals(featureInput.getAttribute("key1"), "value1", "wrong attribute value for key1");
        Assert.assertEquals(featureInput.getAttribute("key2"), null, "wrong attribute value for key2 (not present)");

        Assert.assertEquals(featureInput.getName(), "myName");
        Assert.assertEquals(featureInput.getFeaturePath(), "myFile");
    }

    @Test
    public void testFeatureKeyValuePairsSpecifiedSameValue() {
        final FeatureInput<Feature> featureInput = runCommandLineWithTaggedFeatureInput("argName:myName,key1=value,key2=value", "myFile");

        Assert.assertEquals(featureInput.getAttribute("key1"), "value", "wrong attribute value for key1");
        Assert.assertEquals(featureInput.getAttribute("key2"), "value", "wrong attribute value for key2");

        Assert.assertEquals(featureInput.getName(), "myName");
        Assert.assertEquals(featureInput.getFeaturePath(), "myFile");
    }

    @DataProvider(name = "KeyValuesDataProviderForTestingNull")
    public Object[][] getKeyValuesDataProviderForTestingNull() {
        return new Object[][] {
                { "argName:myName,key1=value1,key2=value2", "myFile" },
                { "argName:myName,null=value", "myFile"},
                { "argName:myName", "myFile" },
                { "argName", "myFile" },
                //{ "argName", "null" }, // "null" has special meaning to the CLP
                { "argName:null", "myFile" },
                { "argName:null", "null" },
        };
    }

    @Test(dataProvider = "KeyValuesDataProviderForTestingNull", expectedExceptions = IllegalArgumentException.class)
    public void testFeatureValuesForNullKey(final String argWithTags, final String inputValue ) {
        final FeatureInput<Feature> featureInput = runCommandLineWithTaggedFeatureInput(argWithTags, inputValue);
        featureInput.getAttribute(null);
    }

    @Test
    public void testFeatureCodecCache() {
        Assert.assertEquals(getVariantFeatureInputWithCachedCodec().getFeatureCodecClass(), VCFCodec.class);
    }

    @SuppressWarnings("unchecked")
    @Test
    public void testFeatureCodecCacheSerialization() throws IOException, ClassNotFoundException {
        final FeatureInput<VariantContext>featureInput = getVariantFeatureInputWithCachedCodec();

        final FeatureInput<VariantContext> roundTrippedFeatureInput = SparkTestUtils.roundTripThroughJavaSerialization(featureInput);
        Assert.assertNotNull(roundTrippedFeatureInput);

        // we expect to lose the cached feature codec class on serialization, but retain the feature path
        Assert.assertNull(roundTrippedFeatureInput.getFeatureCodecClass());
        Assert.assertEquals(featureInput.getFeaturePath(), roundTrippedFeatureInput.getFeaturePath());
    }

    @SuppressWarnings("unchecked")
    private FeatureInput<VariantContext> getVariantFeatureInputWithCachedCodec() {
        final File inputVCFFile = new File(FEATURE_INPUT_TEST_DIRECTORY, "minimal_vcf4_file.vcf");
        final FeatureInput<VariantContext> featureInput = new FeatureInput<>(inputVCFFile.getAbsolutePath());
        Assert.assertNull(featureInput.getFeatureCodecClass());

        final FeatureCodec<? extends Feature, ?> codec = FeatureManager.getCodecForFile(featureInput.toPath());
        featureInput.setFeatureCodecClass((Class<FeatureCodec<VariantContext, ?>>)codec.getClass());

        return featureInput;
    }

    @Test
    public void testToString() {
        final FeatureInput<Feature> namelessFeatureInput = runCommandLineWithTaggedFeatureInput("argName","file1");
        final FeatureInput<Feature> namedFeatureInput = runCommandLineWithTaggedFeatureInput("argName:name", "file1");
        final FeatureInput<Feature> namelessGenomicsDB = runCommandLineWithTaggedFeatureInput("argName", "gendb://file1");
        final FeatureInput<Feature> namedGenomicsDB = runCommandLineWithTaggedFeatureInput("argName:name", "gendb://file1");

        Assert.assertEquals(namelessFeatureInput.toString(), new File("file1").getAbsolutePath(), "String representation of nameless FeatureInput incorrect");
        Assert.assertEquals(namedFeatureInput.toString(), "name:" + new File("file1").getAbsolutePath(), "String representation of named FeatureInput incorrect");
        Assert.assertEquals(namelessGenomicsDB.toString(), "gendb://" + new File("file1").getAbsolutePath(), "String representation of nameless FeatureInput with genomicsDB path incorrect");
        Assert.assertEquals(namedGenomicsDB.toString(), "name:gendb://" + new File("file1").getAbsolutePath(), "String representation of named FeatureInput with genomicsDB path incorrect");
    }

    @DataProvider(name = "HasUserSuppliedNameData")
    public Object[][] hasUserSuppliedNameData() {
        return new Object[][] {
                {"argName", "hdfs://localhost/user/my.vcf", false},
                {"argName:myname", "hdfs://localhost/user/my.vcf", true},
                {"argName:myname,key1=value1", "hdfs://localhost/user/my.vcf", true},
                {"argName:myname//", "hdfs://localhost/user/my.vcf", true},
                {"argName:myname//", "/user/my.vcf", true},
                {"argName", "/user/my.vcf", false},
        };
    }

    @Test(dataProvider = "HasUserSuppliedNameData")
    public void testHasUserSuppliedName(final String argWithTags, final String inputValue, final boolean isUserSupplied) {
        final FeatureInput<Feature> input = runCommandLineWithTaggedFeatureInput(argWithTags, inputValue);
        Assert.assertEquals(input.hasUserSuppliedName(), isUserSupplied);
    }

}
