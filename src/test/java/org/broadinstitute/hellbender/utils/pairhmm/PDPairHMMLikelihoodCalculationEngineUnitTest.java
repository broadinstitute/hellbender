package org.broadinstitute.hellbender.utils.pairhmm;

import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.BufferedLineReader;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeArguments;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

public class PDPairHMMLikelihoodCalculationEngineUnitTest extends GATKBaseTest  {

    static final String DRAGEN_GATK_TEST_ASSERT_FILE = largeFileTestDir + "expected.PDHMM.hmmresults.txt";
    static final double DOUBLE_ASSERTION_DELTA = 0.0001;

    @DataProvider(name="PDPairHMMResultsModes")
    public Object[][] PairHMMResultsModes() {
        return new Object[][]{
                {PDPairHMM.Implementation.LOGLESS_CACHING},
                {PDPairHMM.Implementation.AVX_LOGLESS_CACHING},
        };
    }

    // This test is intended to make it simple to take the outputs of the PDHMM from one impelementation and test its outputs over another.
    @Test(dataProvider = "PDPairHMMResultsModes")
    public void testInputReconstitutionDRAGEN_GATK_TEST_FILE(PDPairHMM.Implementation implementation) {
        final PDPairHMM pdPariHMM = implementation.makeNewHMM(new PairHMMNativeArguments());

        try (final FileInputStream fis= new FileInputStream(DRAGEN_GATK_TEST_ASSERT_FILE);
             final BufferedLineReader br = new BufferedLineReader(fis)) {

            br.lines().skip(1).forEach(line -> {

                final String[] split = line.split("\t");
                //debugOutputStream.write("\n" + new String(alleleBases) + "\t" + Arrays.toString(allelePDBases) + "\t" + new String(readBases) + "\t" + SAMUtils.phredToFastq(readQuals) + "\t" + SAMUtils.phredToFastq(readInsQuals) + "\t" + SAMUtils.phredToFastq(readDelQuals) + "\t" + SAMUtils.phredToFastq(overallGCP) + "\t" + String.format("%e",lk));
                byte[] alleleBases = split[0].getBytes(StandardCharsets.UTF_8);
                byte[] allelePDBases = ArrayUtils.toPrimitive(
                        Arrays.stream(split[1].substring(1, split[1].length() - 1).split(","))
                                .map(num -> Byte.parseByte(num.trim())).toArray(Byte[]::new));
                byte[] readBases = split[2].getBytes(StandardCharsets.UTF_8);
                byte[] readQuals = SAMUtils.fastqToPhred(split[3]);
                byte[] readInsQuals = SAMUtils.fastqToPhred(split[4]);
                byte[] readDelQuals = SAMUtils.fastqToPhred(split[5]);
                byte[] overallGCP = SAMUtils.fastqToPhred(split[6]);
                double expected = Double.parseDouble(split[7]);

                // There will be Negative infinities in here in the case that we skipped likelihood computation for a variant
                if (expected != Double.NEGATIVE_INFINITY) {
                    pdPariHMM.initialize(200, 500);// These are arbitrary initializations that should have no bearing on this operation
                    double actual = pdPariHMM.computeReadLikelihoodGivenHaplotypeLog10(alleleBases, allelePDBases, readBases, readQuals, readInsQuals, readDelQuals, overallGCP, false, null, null);
                    Assert.assertEquals(actual, expected, DOUBLE_ASSERTION_DELTA, String.format("Mismatching score actual: %e expected: %e computed on line %s", actual, expected, line));
                }
            });
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}