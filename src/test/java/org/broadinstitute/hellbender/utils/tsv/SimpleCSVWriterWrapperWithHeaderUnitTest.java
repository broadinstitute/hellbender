package org.broadinstitute.hellbender.utils.tsv;

import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.LineReader;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;

public class SimpleCSVWriterWrapperWithHeaderUnitTest extends GATKBaseTest {

    @Test (groups = "bucket")
    public void testWriteToBucketPathEquivalentToLocalPath() throws IOException {
        Path bucketPath = IOUtils.getPath(BucketUtils.getTempFilePath(
                getGCPTestStaging() + "testWriteToBucketPathEquivalentToLocalPath", ".tsv"));
        Path localPath = IOUtils.getPath(createTempFile("testWriteToBucketPathEquivalentToLocalPath", ".tsv").getPath());

        try (SimpleXSVWriter bucketWriter = new SimpleXSVWriter(bucketPath, '\t');
             SimpleXSVWriter localWriter = new SimpleXSVWriter(localPath, '\t')) {

            String[] header = new String[]{"a", "b", "c"};
            bucketWriter.setHeaderLine(Arrays.asList(header), true);
            localWriter.setHeaderLine(Arrays.asList(header), true);

            for (int i = 0; i < 100; i++) {
                SimpleXSVWriter.LineBuilder bucketLine = bucketWriter.getNewLineBuilder();
                SimpleXSVWriter.LineBuilder localLine = localWriter.getNewLineBuilder();
                int finalI = i;
                Arrays.stream(header).forEach(column -> {
                    bucketLine.setColumn(column, Integer.toString(finalI));
                    localLine.setColumn(column, Integer.toString(finalI));
                });
            }
        }
        IntegrationTestSpec.assertEqualTextFiles(bucketPath, localPath, null);
    }

    @Test
    public void testFillingInBlankLines() throws IOException {
        Path outputPath = IOUtils.getPath(createTempFile("testWriteToBucketPathEquivalentToLocalPath", ".csv").getPath());

        try (SimpleXSVWriter localWriter = new SimpleXSVWriter(outputPath, ',')) {
            String[] header = new String[]{"a", "b", "c","d"};
            localWriter.setHeaderLine(Arrays.asList(header), true);

            for (int i = 0; i < 100; i++) {
                SimpleXSVWriter.LineBuilder localLine = localWriter.getNewLineBuilder();
                Arrays.stream(header).forEach(column -> {
                    localLine.setColumn("b", "10");
                    localLine.fill("0");
                    localLine.setColumn("d", "1");
                });
            }
        }

        try (final FileInputStream fis= new FileInputStream(outputPath.toFile());
             final BufferedLineReader br = new BufferedLineReader(fis)) {
            Assert.assertEquals(br.readLine(), "a,b,c,d");
            int lineCount = 0;
            while (lineCount++ < 100) {
                Assert.assertEquals(br.readLine(), "0,10,0,1");
            }
        }
    }

    @Test (expectedExceptions = IllegalStateException.class)
    public void testWrongNumberOfLines() throws IOException {
        Path outputPath = IOUtils.getPath(createTempFile("testWriteToBucketPathEquivalentToLocalPath", ".csv").getPath());

        try (SimpleXSVWriter localWriter = new SimpleXSVWriter(outputPath, ',')) {
            String[] header = new String[]{"a", "b", "c","d"};
            localWriter.setHeaderLine(Arrays.asList(header), true );

            localWriter.getNewLineBuilder().setRow(Arrays.asList("1","2","3","4","5"));
        }
    }

    @Test (expectedExceptions = GATKException.class)
    public void testMissingHeader() throws IOException {
        Path outputPath = IOUtils.getPath(createTempFile("testWriteToBucketPathEquivalentToLocalPath", ".csv").getPath());

        try (SimpleXSVWriter localWriter = new SimpleXSVWriter(outputPath, ',')) {
            localWriter.getNewLineBuilder();
        }
    }
}