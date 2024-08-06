package org.broadinstitute.hellbender.engine;

import htsjdk.beta.io.IOPathUtils;
import htsjdk.beta.io.bundle.*;
import htsjdk.io.IOPath;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public final class MultiVariantWalkerIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return MultiVariantWalker.class.getSimpleName();
    }

    public static File getTestDataDir() {
        return new File(publicTestDir + "org/broadinstitute/hellbender/engine/MultiVariantDataSource/");
    }

    @CommandLineProgramProperties(
            summary = "TestGATKToolWithFeatures",
            oneLineSummary = "TestGATKToolWithFeatures",
            programGroup = TestProgramGroup.class
    )
    private static final class TestMultiVariantWalker extends MultiVariantWalker {

        int count = 0;
        SimpleInterval locus;

        @Override
        public void apply(
                final VariantContext variant,
                final ReadsContext readsContext,
                final ReferenceContext referenceContext,
                final FeatureContext featureContext )
        {
            count++;

            //make sure we only move forward; if getSequenceDictionary returns null then there is only
            //a single input, and it has no sequence dictionary, so skip the test since we're just
            //iterating over a single input in order
            if (locus != null && getSequenceDictionaryForDrivingVariants() != null) {
                int locDiff = IntervalUtils.compareLocatables(
                        locus,
                        new SimpleInterval(variant),
                        getSequenceDictionaryForDrivingVariants());
                Assert.assertTrue(locDiff == 0 || locDiff == -1);
            }
            locus = new SimpleInterval(variant);
        }
    }

    @Test(expectedExceptions = UserException.class)
    public void testDuplicateSources() throws Exception {
        final TestMultiVariantWalker tool = new TestMultiVariantWalker();

        final List<String> args = new ArrayList<>();
        File testFile = new File(getTestDataDir(), "baseVariants.vcf");
        args.add("--variant");
        args.add(testFile.getAbsolutePath());
        args.add("--variant");
        args.add(testFile.getAbsolutePath()); // add it again
        tool.instanceMain(args.toArray(new String[args.size()]));
    }

    @DataProvider(name="variantFiles")
    public Object[][] getVariantFiles() {
        return new Object[][]
                {
                        {Arrays.asList(new File(getTestDataDir(), "baseVariants.vcf")), null, 26},
                        {Arrays.asList(new File(getTestDataDir(), "interleavedVariants_1.vcf")), null, 13},
                        {Arrays.asList(
                                new File(getTestDataDir(), "interleavedVariants_1.vcf"),
                                new File(getTestDataDir(), "interleavedVariants_2.vcf")), null, 26},
                        {Arrays.asList(
                                new File(getTestDataDir(), "splitVariants_1.vcf"),
                                new File(getTestDataDir(), "splitVariants_2.vcf")), null, 26},

                        // with intervals
                        {Arrays.asList(new File(getTestDataDir(), "baseVariants.vcf")), "1", 14},
                        {Arrays.asList(new File(getTestDataDir(), "interleavedVariants_1.vcf")), "1", 7},
                        {Arrays.asList(
                                new File(getTestDataDir(), "interleavedVariants_1.vcf"),
                                new File(getTestDataDir(), "interleavedVariants_2.vcf")), "1", 14},
                        {Arrays.asList(
                                new File(getTestDataDir(), "interleavedVariants_1.vcf"),
                                new File(getTestDataDir(), "interleavedVariants_2.vcf")), "2:200-600", 3},
            };
    }

    @Test(dataProvider = "variantFiles")
    public void testVariantOrder(final List<File> inputFiles, final String interval, final int expectedCount) throws Exception {
        final TestMultiVariantWalker tool = new TestMultiVariantWalker();

        final List<String> args = new ArrayList<>(inputFiles.size() * 2);
        inputFiles.forEach(f -> {args.add("--variant");
        args.add(f.getAbsolutePath());});
        if (interval != null) {
            args.add("-L");
            args.add(interval);
        }

        tool.instanceMain(args.toArray(new String[args.size()]));
        Assert.assertEquals(tool.count, expectedCount);
    }

    // tests using bundles
    @DataProvider(name="variantBundles")
    public Object[][] getVariantBundles() throws IOException {
        return new Object[][]
                {
                        {createRemoteBundleForFile(
                                new File(getTestDataDir(), "interleavedVariants_1.vcf"),
                                new File(getTestDataDir(), "interleavedVariants_1.vcf.idx"),
                                new File(getTestDataDir(), "interleavedVariants_2.vcf"),
                                new File(getTestDataDir(), "interleavedVariants_2.vcf.idx")
                            ), "1", 14},
                        {createRemoteBundleForFile(
                                new File(getTestDataDir(), "interleavedVariants_1.vcf"),
                                new File(getTestDataDir(), "interleavedVariants_1.vcf.idx"),
                                new File(getTestDataDir(), "interleavedVariants_2.vcf"),
                                new File(getTestDataDir(), "interleavedVariants_2.vcf.idx")
                        ), "2:200-600", 3},
                };
    };

    @Test(dataProvider = "variantBundles")
    public void testVariantsFromBundleOrder(final IOPath inputBundleFile, final String interval, final int expectedCount) {
        final TestMultiVariantWalker tool = new TestMultiVariantWalker();

        final List<String> args = new ArrayList<>();
        args.add("--variant");
        args.add(inputBundleFile.getRawInputString());
        if (interval != null) {
            args.add("-L");
            args.add(interval);
        }

        tool.instanceMain(args.toArray(new String[args.size()]));
        Assert.assertEquals(tool.count, expectedCount);
    }

    // copy  two .vcfs to a temporary remote location; create a bundle list with two bundles, each containing
    // a reference to the remote vcf and it's local companion index file, and then write the whole bundle out to
    // a temporary remote bundle file
    private static IOPath createRemoteBundleForFile(
            final File vcf1,
            final File index1,
            final File vcf2,
            final File index2) throws IOException {
        final String remotePath = BucketUtils.randomRemotePath(getGCPTestStaging() + "remoteBundles", "remote_bundle_test", "dir");
        final Path remoteDirPath = IOUtils.getPath(remotePath + "/");

        Files.createDirectory(remoteDirPath);
        Assert.assertTrue(Files.exists(remoteDirPath));
        Assert.assertTrue(Files.isDirectory(remoteDirPath));
        final Path remoteVCF1 = IOUtils.getPath(remotePath + "/" + vcf1.getName());
        final Path remoteVCF2 = IOUtils.getPath(remotePath + "/" + vcf2.getName());
        final Path c_remoteVCF1 = Files.copy(vcf1.toPath(), remoteVCF1);
        if (!c_remoteVCF1.equals(remoteVCF1)) {
            throw new IOException("Not equal " + vcf1 + " to " + remoteVCF1);
        }
        final Path c_remoteVCF2 = Files.copy(vcf2.toPath(), remoteVCF2);
        if (!c_remoteVCF2.equals(remoteVCF2)) {
            throw new IOException("Not equal " + vcf2 + " to " + remoteVCF2);
        }

        final Bundle variantsBundle1 = new BundleBuilder()
                .addPrimary(new IOPathResource(new GATKPath(remoteVCF1.toUri().toString()), BundleResourceType.CT_VARIANT_CONTEXTS))
                .addSecondary(new IOPathResource(new GATKPath(index1.toURI().toString()), BundleResourceType.CT_VARIANTS_INDEX))
                .build();
        final Bundle variantsBundle2 = new BundleBuilder()
                .addPrimary(new IOPathResource(new GATKPath(remoteVCF2.toUri().toString()), BundleResourceType.CT_VARIANT_CONTEXTS))
                .addSecondary(new IOPathResource(new GATKPath(index2.toURI().toString()), BundleResourceType.CT_VARIANTS_INDEX))
                .build();

        final List<Bundle> bundles = Arrays.asList(variantsBundle1, variantsBundle2);
        final GATKPath remoteBundlePath = new GATKPath(remoteDirPath.resolve("remote_bundle.json").toUri().toString());
        IOPathUtils.writeStringToPath(remoteBundlePath, BundleJSON.toJSON(bundles));
        return remoteBundlePath;
    }

    @CommandLineProgramProperties(
        summary = "TestGATKToolWithFeatures",
        oneLineSummary = "TestGATKToolWithFeatures",
        programGroup = TestProgramGroup.class
    )
    private static final class TestMultiVariantWalkerIterator extends MultiVariantWalker {
        String expectedIDOrder[];
        int count = 0;

        @Override
        public void apply(
                final VariantContext variant,
                final ReadsContext readsContext,
                final ReferenceContext referenceContext,
                final FeatureContext featureContext )
        {
            Assert.assertEquals(variant.getID(), expectedIDOrder[count]);
            count++;
        }
    }

    @Test
    public void testIteratorOverlapping() {
        //Test interleaved files that include some variants that start at the same position in both files
        final TestMultiVariantWalkerIterator tool = new TestMultiVariantWalkerIterator();
        tool.expectedIDOrder = new String[] {
                "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n",
                "o", "o_overlap",
                "p", "q", "r", "s", "t", "u", "v", "w",
                "x", "x_overlap",
                "y", "z"
        };

        final List<String> args = new ArrayList<>();
        File testFile1 = new File(getTestDataDir(), "interleavedVariants_1_WithOverlap.vcf");
        args.add("--variant");
        args.add(testFile1.getAbsolutePath());
        File testFile2 = new File(getTestDataDir(), "interleavedVariants_2_WithOverlap.vcf");
        args.add("--variant");
        args.add(testFile2.getAbsolutePath());

        tool.instanceMain(args.toArray(new String[args.size()]));
        Assert.assertEquals(tool.count, 28);
    }

    @Test
    public void testGetCompatibleHeader() throws Exception {
        final TestMultiVariantWalker tool = new TestMultiVariantWalker();

        final List<String> args = new ArrayList<>();
        File testFile1 = new File(getTestDataDir(), "interleavedVariants_1.vcf");
        args.add("--variant");
        args.add(testFile1.getAbsolutePath());
        File testFile2 = new File(getTestDataDir(), "interleavedVariants_2.vcf");
        args.add("--variant");
        args.add(testFile2.getAbsolutePath());

        tool.instanceMain(args.toArray(new String[args.size()]));
        VCFHeader header = tool.getHeaderForVariants();
        Assert.assertEquals(header.getSequenceDictionary().getSequences().size(), 4 );
    }

    @Test(expectedExceptions = UserException.IncompatibleSequenceDictionaries.class)
    public void testGetConflictingHeader() throws Exception {
        final TestMultiVariantWalker tool = new TestMultiVariantWalker();

        final List<String> args = new ArrayList<>();
        File testFile1 = new File(getTestDataDir(), "baseVariants.vcf");
        args.add("--variant");
        args.add(testFile1.getAbsolutePath());
        File testFile2 = new File(getTestDataDir(), "baseVariantsConflictingDictionary.vcf");
        args.add("--variant");
        args.add(testFile2.getAbsolutePath());

        tool.instanceMain(args.toArray(new String[args.size()]));
    }

    @Test
    public void testNoDictionaryForOnlyInput() throws Exception {
        final TestMultiVariantWalker tool = new TestMultiVariantWalker();

        final List<String> args = new ArrayList<>();
        File testFile1 = new File(getTestDataDir(), "interleavedVariants_1_NoDict.vcf");
        args.add("--variant");
        args.add(testFile1.getAbsolutePath());

        tool.instanceMain(args.toArray(new String[args.size()]));

        Assert.assertNull(tool.getSequenceDictionaryForDrivingVariants());
    }

    @Test
    public void testNoDictionaryForOneInput() throws Exception {
        final TestMultiVariantWalker tool = new TestMultiVariantWalker();

        final List<String> args = new ArrayList<>();
        File testFile1 = new File(getTestDataDir(), "interleavedVariants_1.vcf");
        args.add("--variant");
        args.add(testFile1.getAbsolutePath());
        File testFile2 = new File(getTestDataDir(), "interleavedVariants_2_NoDict.vcf");
        args.add("--variant");
        args.add(testFile2.getAbsolutePath());

        tool.instanceMain(args.toArray(new String[args.size()]));
        Assert.assertNotNull(tool.getSequenceDictionaryForDrivingVariants());
    }

    @Test(expectedExceptions = UserException.class)
    public void testNoDictionaryForAllInputs() throws Exception {
        final TestMultiVariantWalker tool = new TestMultiVariantWalker();

        final List<String> args = new ArrayList<>();
        File testFile1 = new File(getTestDataDir(), "interleavedVariants_1_NoDict.vcf");
        args.add("--variant");
        args.add(testFile1.getAbsolutePath());
        File testFile2 = new File(getTestDataDir(), "interleavedVariants_2_NoDict.vcf");
        args.add("--variant");
        args.add(testFile2.getAbsolutePath());

        tool.instanceMain(args.toArray(new String[args.size()]));
    }

}
