package org.broadinstitute.hellbender.tools;

import htsjdk.beta.io.IOPathUtils;
import htsjdk.beta.io.bundle.*;
import htsjdk.io.IOPath;
import htsjdk.samtools.*;
import htsjdk.samtools.cram.ref.ReferenceSource;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.ReadsPathDataSource;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;

public final class PrintReadsIntegrationTest extends AbstractPrintReadsIntegrationTest {

    @Test(dataProvider="testingData")
    public void testFileToFileWithMD5(String fileIn, String extOut, String reference) throws Exception {
        doFileToFile(fileIn, extOut, reference, true);
    }

    @Test
    public void testNoConflictPG() throws IOException {
        final File inFile = new File(TEST_DATA_DIR, "print_reads_withPG.sam");
        final File outFile = GATKBaseTest.createTempFile("testNoConflictRG", ".sam");
        final String[] args = new String[] {
                "--input" , inFile.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_SAM_PROGRAM_RECORD,
                "--output", outFile.getAbsolutePath()
        };
        runCommandLine(args);

        //Make sure contents are the same
        SamAssertionUtils.assertSamsEqual(outFile, inFile);

        //input has GATK PrintReads not NOT GATK PrintReads.1 in headers
        Assert.assertNotNull(SamReaderFactory.makeDefault().open(inFile).getFileHeader().getProgramRecord("GATK PrintReads"));
        Assert.assertNull(SamReaderFactory.makeDefault().open(inFile).getFileHeader().getProgramRecord("GATK PrintReads.1"));

        //output has both GATK PrintReads and GATK PrintReads.1 in headers
        Assert.assertNotNull(SamReaderFactory.makeDefault().open(outFile).getFileHeader().getProgramRecord("GATK PrintReads"));
        Assert.assertNotNull(SamReaderFactory.makeDefault().open(outFile).getFileHeader().getProgramRecord("GATK PrintReads.1"));
    }

    @DataProvider
    public Object[][] getHttpPaths(){
        final String bam = "gs://hellbender/test/resources/benchmark/CEUTrio.HiSeq.WEx.b37.NA12892.bam";
        final String bai = "gs://hellbender/test/resources/benchmark/CEUTrio.HiSeq.WEx.b37.NA12892.bam.bai";
        final String cram = "gs://hellbender/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.v3.0.samtools.cram";
        final String crai = "gs://hellbender/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.v3.0.samtools.cram.crai";
        final List<SimpleInterval> bamIntervals = Arrays.asList(new SimpleInterval("3",1_000_000, 1_000_001),
               new SimpleInterval("3", 1_000_003, 1_000_100),
                new SimpleInterval("20", 1_099_000, 1_100_000));

        final List<SimpleInterval> cramIntervals = Arrays.asList(new SimpleInterval("20", 9_999_902, 10_000_000));
        return new Object[][]{
                {BucketUtils.bucketPathToPublicHttpUrl(bam), BucketUtils.bucketPathToPublicHttpUrl(bai), bam, bai, bamIntervals, 528L},
                {BucketUtils.createSignedUrlToGcsObject(bam, 1L), BucketUtils.createSignedUrlToGcsObject(bai, 1L), bam, bai, bamIntervals, 528L},
                {BucketUtils.bucketPathToPublicHttpUrl(cram), BucketUtils.bucketPathToPublicHttpUrl(crai), cram, crai, cramIntervals, 112L},
                {BucketUtils.createSignedUrlToGcsObject(cram, 1L), BucketUtils.createSignedUrlToGcsObject(crai, 1L), cram, crai, cramIntervals, 112L}

        };
    }

    @Test(groups = {"cloud", "bucket"}, dataProvider = "getHttpPaths")
    public void testHttpPaths(String reads, String index, String nonHttpReads, String nonHttpIndex, List<SimpleInterval> intervals, long expectedNumberOfReads) throws IOException {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        final File out = createTempFile("out", ".bam");
        // this test reads tiny amounts of data from multiple places, if you don't set the prefetcher to a lower number
        // it loads large amounts of data that slows the test down significantly for no good reason
        args.addInput(reads)
                .add(StandardArgumentDefinitions.CLOUD_PREFETCH_BUFFER_LONG_NAME, 1)
                .add(StandardArgumentDefinitions.CLOUD_INDEX_PREFETCH_BUFFER_LONG_NAME, 1)
                .add("read-index", index)
                .addReference(GATKBaseTest.b37Reference)
                .addOutput(out);
        intervals.forEach(args::addInterval);
        runCommandLine(args);

        final ArgumentsBuilder args2 = new ArgumentsBuilder();
        final File out2 = createTempFile("out", ".bam");
        args2.addInput(nonHttpReads)
                .add("read-index", nonHttpIndex)
                .add(StandardArgumentDefinitions.CLOUD_PREFETCH_BUFFER_LONG_NAME, 1)
                .add(StandardArgumentDefinitions.CLOUD_INDEX_PREFETCH_BUFFER_LONG_NAME, 1)
                .addReference(GATKBaseTest.b37Reference)
                .addOutput(out2);
        intervals.forEach(args2::addInterval);
        runCommandLine(args2);

        try(final ReadsDataSource reader = new ReadsPathDataSource(out.toPath())){
            final long count = Utils.stream(reader).count();
            Assert.assertEquals( count, expectedNumberOfReads);
        }

        SamAssertionUtils.assertEqualBamFiles(out, out2, false, ValidationStringency.DEFAULT_STRINGENCY);
    }

    @Test(groups = {"cloud"})
    public void testPrintReadsWithReferenceBundle() throws IOException {
        // test that both reading and writing a cram work when the reference is specified via a bundle where the
        // reference, reference index, and reference dictionary are all in different buckets
        final IOPath testFastaFile = new GATKPath(getTestDataDir() + "/print_reads.fasta");
        final IOPath testIndexFile = new GATKPath(getTestDataDir() + "/print_reads.fasta.fai");
        final IOPath testDictFile = new GATKPath(getTestDataDir() + "/print_reads.dict");

        final String targetBucketName = BucketUtils.randomRemotePath(getGCPTestStaging(), "testPrintReadsWithReferenceBundle", "") + "/";
        final IOPath targetBucket = new GATKPath(targetBucketName);
        IOUtils.deleteOnExit(targetBucket.toPath());

        final Path remoteFasta = Files.copy(testFastaFile.toPath(), new GATKPath(targetBucketName + "print_reads.fasta").toPath());
        final IOPath targetIndex = new GATKPath(targetBucketName + "refindex/print_reads.fasta.fai");
        final Path remoteFastaIndex = Files.copy(testIndexFile.toPath(), targetIndex.toPath());
        final IOPath targetDict = new GATKPath(targetBucketName + "refdict/print_reads.dict");
        final Path remoteFastaDict = Files.copy(testDictFile.toPath(), targetDict.toPath());

        // create a bundle with the remote reference, index, and dict files
        final Bundle refBundle = new BundleBuilder()
                .addPrimary(new IOPathResource(new GATKPath(remoteFasta.toUri().toString()), BundleResourceType.CT_HAPLOID_REFERENCE))
                .addSecondary(new IOPathResource(new GATKPath(remoteFastaIndex.toUri().toString()), BundleResourceType.CT_REFERENCE_INDEX))
                .addSecondary(new IOPathResource(new GATKPath(remoteFastaDict.toUri().toString()), BundleResourceType.CT_REFERENCE_DICTIONARY))
                .build();
        final IOPath bundleFilePath = new GATKPath(targetBucketName + "refBundle.json");
        IOPathUtils.writeStringToPath(bundleFilePath, BundleJSON.toJSON(refBundle));

        final IOPath targetOutCRAM = new GATKPath(IOUtils.createTempFile("testReferenceSequenceForNioBundle", ".cram").getAbsolutePath());
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addInput(getTestDataDir() + "/print_reads.cram")
                .addReference(bundleFilePath.toString())
                .addOutput(targetOutCRAM.toString());
                runCommandLine(args);

        int count = 0;
        try (final SamReader in = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .referenceSource(new ReferenceSource(bundleFilePath.toPath()))
                .open(targetOutCRAM.toPath())) {
            for (@SuppressWarnings("unused") final SAMRecord rec : in) {
                count++;
            }
        }
        Assert.assertEquals(count, 8);
    }

}