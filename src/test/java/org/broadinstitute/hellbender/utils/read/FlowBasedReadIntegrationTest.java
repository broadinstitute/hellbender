package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.broadinstitute.hellbender.tools.walkers.variantrecalling.FlowTestConstants;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;

public class FlowBasedReadIntegrationTest extends GATKBaseTest {

    private final Logger logger = LogManager.getLogger(this.getClass());

    @DataProvider(name = "flowBasedReadFiles")
    public Object[][] getFlowBasedReadFiles() {
        final Object[][]        testData = {

                { publicTestDir + FlowTestConstants.FEATURE_MAPPING_DATA_DIR + "/snv_feature_mapper_input.bam",
                        null, "TGCA", 100, false }
        };

        return testData;
    }

    @Test(dataProvider = "flowBasedReadFiles")
    public void testReads(final String inputFile, final String outputPrefix, final String flowOrder,
                          final int limitCount, final boolean simulate) throws IOException {

        // create argument block
        FlowBasedArgumentCollection fbargs = new FlowBasedArgumentCollection();
        if ( simulate ) {
            fbargs.fillingValue = 0.01;
            //fbargs.remove_longer_than_one_indels = true;
            fbargs.lumpProbs = true;
            fbargs.onlyInsOrDel = true;
            fbargs.removeOneToZeroProbs = true;
        }

        // initialize reader
        final SamReader reader = SamReaderFactory.makeDefault().open(new File(inputFile));
        int count = 0 ;
        Iterator<SAMRecord> i;
        for (i = reader.iterator(), count=0; (i.hasNext()) && (count < limitCount); count++  ) {
            final FlowBasedRead   fbr = FlowBasedReadUtils.convertToFlowBasedRead(
                                            new SAMRecordToGATKReadAdapter(i.next()),
                                            reader.getFileHeader());
            fbr.applyAlignment();
            Assert.assertEquals(fbr.totalKeyBases(), fbr.seqLength());

            if ( limitCount < 1000 && outputPrefix != null ) {
                try ( final FileWriter fos = new FileWriter(outputPrefix + "." + Integer.toString(count) + ".key.txt") ) {
                    fbr.writeKey(fos);
                }
                try ( final FileWriter fos = new FileWriter(outputPrefix + "." + Integer.toString(count) + ".matrix.txt") ) {
                    fbr.writeMatrix(fos);
                }
            } else if ( (count % 1000) == 0 ) {
                logger.debug("record count: " + count);
            }
        }
    }
}
