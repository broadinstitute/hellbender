package org.broadinstitute.hellbender.tools.examples;

import htsjdk.samtools.util.Log;
import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

public class ExampleMergingMultiFeatureWalkerIntegrationTest extends CommandLineProgramTest {

    @Test
    public void testDictionarySubsetIsOK() {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VERBOSITY_NAME, Log.LogLevel.ERROR.name());
        // full HG38 dictionary with 3366 entries
        argsBuilder.add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, FULL_HG38_DICT);
        // subset of HG38 dictionary with 13 entries
        argsBuilder.add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b38_reference_20_21);
        // bam has full HG38 dictionary
        argsBuilder.add(StandardArgumentDefinitions.INPUT_LONG_NAME, largeFileTestDir + "NA12878.alignedHg38.duplicateMarked.baseRealigned.bam");
        // no dictionary, no sample names, with a single feature
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_LONG_NAME, packageRootTestDir + "engine/tiny_hg38.baf.txt");
        final ExampleMergingMultiFeatureWalker example = new ExampleMergingMultiFeatureWalker();
        example.instanceMain(argsBuilder.getArgsArray());

        // full HG38 dictionary has 3366 entries
        Assert.assertEquals(example.getDictionary().size(), 3366);
        // feature file has 1 feature, no samples, no dictionary
        Assert.assertEquals(example.features.size(), 1);
        Assert.assertEquals(example.getSampleNames().size(), 0);
    }

    @Test(expectedExceptions = { UserException.class })
    public void testMixedDictionariesAreNotOK() {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VERBOSITY_NAME, Log.LogLevel.ERROR.name());
        argsBuilder.add(StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME, "true");
        // full HG38 dictionary with 3366 entries
        argsBuilder.add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, FULL_HG38_DICT);
        // subset of HG37 dictionary has conflicting names
        argsBuilder.add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21);
        // no dictionary, no sample names, with a single feature
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_LONG_NAME, packageRootTestDir + "engine/tiny_hg38.baf.txt");
        final ExampleMergingMultiFeatureWalker example = new ExampleMergingMultiFeatureWalker();
        example.instanceMain(argsBuilder.getArgsArray());
    }

    @Test(expectedExceptions = { UserException.class })
    public void testMisorderedDictionariesAreNotOK() {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VERBOSITY_NAME, Log.LogLevel.ERROR.name());
        // subset of HG38 dictionary with misordered contigs
        argsBuilder.add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, packageRootTestDir + "engine/wrongOrder.dict");
        // subset of HG38 dictionary with 13 entries
        argsBuilder.add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b38_reference_20_21);
        // no dictionary, no sample names, with a single feature
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_LONG_NAME, packageRootTestDir + "engine/tiny_hg38.baf.txt");
        final ExampleMergingMultiFeatureWalker example = new ExampleMergingMultiFeatureWalker();
        example.instanceMain(argsBuilder.getArgsArray());
    }

    @Test
    public void testGetDictionaryAndSamplesFromBCIFile() {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VERBOSITY_NAME, Log.LogLevel.ERROR.name());
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_LONG_NAME, packageRootTestDir + "engine/tiny_hg38.baf.bci");
        final ExampleMergingMultiFeatureWalker example = new ExampleMergingMultiFeatureWalker();
        example.instanceMain(argsBuilder.getArgsArray());

        // bci file has chr20, chr21, and alts for those contigs: 13 entries
        Assert.assertEquals(example.getDictionary().size(), 13);
        // feature file has 1 feature
        Assert.assertEquals(example.features.size(), 1);
        // feature file has 1 sample name
        Assert.assertEquals(example.getSampleNames().size(), 1);
    }

    @Test
    public void testContigsOrder() {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VERBOSITY_NAME, Log.LogLevel.ERROR.name());
        // full HG38 dictionary with 3366 entries
        argsBuilder.add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, FULL_HG38_DICT);
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_LONG_NAME, packageRootTestDir + "engine/tiny_hg38.baf.txt");
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_LONG_NAME, packageRootTestDir + "engine/tiny_chr9.baf.txt");
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_LONG_NAME, packageRootTestDir + "engine/tiny_chr10.baf.txt");
        final ExampleMergingMultiFeatureWalker example = new ExampleMergingMultiFeatureWalker();
        example.instanceMain(argsBuilder.getArgsArray());

        Assert.assertEquals(example.features.size(), 5);
        Assert.assertEquals(example.features.get(0).getContig(), "chr9");
        Assert.assertEquals(example.features.get(1).getContig(), "chr10");
        Assert.assertTrue(example.features.get(1).getStart() <= example.features.get(2).getStart());
        Assert.assertEquals(example.features.get(2).getContig(), "chr10");
        Assert.assertTrue(example.features.get(2).getStart() <= example.features.get(3).getStart());
        Assert.assertEquals(example.features.get(3).getContig(), "chr10");
        Assert.assertEquals(example.features.get(4).getContig(), "chr21");
    }

    @Test
    public void testCoordinateOrder() {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VERBOSITY_NAME, Log.LogLevel.ERROR.name());
        // full HG38 dictionary with 3366 entries
        argsBuilder.add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, FULL_HG38_DICT);
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_LONG_NAME, packageRootTestDir + "engine/tiny_chr10.baf.txt");
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_LONG_NAME, packageRootTestDir + "engine/tiny_chr10_2.baf.txt");
        final ExampleMergingMultiFeatureWalker example = new ExampleMergingMultiFeatureWalker();
        example.instanceMain(argsBuilder.getArgsArray());

        Assert.assertEquals(example.features.size(), 5);
        int lastStart = -1;
        for ( final Feature feature : example.features ) {
            Assert.assertTrue(feature.getStart() >= lastStart );
            lastStart = feature.getStart();
        }
    }

    @Test
    public void testRepeatedCoordinates() {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VERBOSITY_NAME, Log.LogLevel.ERROR.name());
        // full HG38 dictionary with 3366 entries
        argsBuilder.add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, FULL_HG38_DICT);
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_LONG_NAME, packageRootTestDir + "engine/tiny_chr10_3.pe.txt");
        final ExampleMergingMultiFeatureWalker example = new ExampleMergingMultiFeatureWalker();
        example.instanceMain(argsBuilder.getArgsArray());

        Assert.assertEquals(example.features.size(), 11);
        int lastStart = -1;
        for ( final Feature feature : example.features ) {
            Assert.assertTrue(feature.getStart() >= lastStart );
            lastStart = feature.getStart();
        }
    }

    @Test
    public void testIntervalSubset() {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VERBOSITY_NAME, Log.LogLevel.ERROR.name());
        argsBuilder.add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, FULL_HG38_DICT);
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_LONG_NAME, packageRootTestDir + "engine/tiny_chr9.baf.txt.gz");
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_LONG_NAME, packageRootTestDir + "engine/tiny_chr10.baf.txt.gz");
        argsBuilder.add(StandardArgumentDefinitions.INTERVALS_LONG_NAME, "chr10:2000-2002");
        final ExampleMergingMultiFeatureWalker example = new ExampleMergingMultiFeatureWalker();
        example.instanceMain(argsBuilder.getArgsArray());

        Assert.assertEquals(example.features.size(), 1);
        final Feature feature = example.features.get(0);
        Assert.assertEquals(feature.getContig(), "chr10");
        Assert.assertEquals(feature.getStart(), 2001);
    }
}
