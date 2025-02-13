package org.broadinstitute.hellbender.utils.recalibration.covariates;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import static java.lang.Math.abs;

public final class CycleCovariateUnitTest extends GATKBaseTest {

    CycleCovariate covariate;
    RecalibrationArgumentCollection RAC;
    SAMReadGroupRecord illuminaReadGroup;

    @BeforeClass
    public void init() {
        RAC = new RecalibrationArgumentCollection();
        covariate = new CycleCovariate(RAC);
        illuminaReadGroup = new SAMReadGroupRecord("MY.ID");
        illuminaReadGroup.setPlatform("illumina");
    }

    @Test
    public void testSimpleCycles() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeaderWithReadGroup(illuminaReadGroup);

        short readLength = 10;
        GATKRead read = ArtificialReadUtils.createRandomRead(header, readLength);
        read.setIsPaired(true);
        read.setReadGroup(illuminaReadGroup.getReadGroupId());

        PerReadCovariateMatrix perReadCovariateMatrix = new PerReadCovariateMatrix(read.getLength(), 1, new CovariateKeyCache());
        covariate.recordValues(read, header, perReadCovariateMatrix, true);
        verifyCovariateArray(perReadCovariateMatrix.getMismatchMatrix(), 1, (short) 1);

        read.setIsReverseStrand(true);
        covariate.recordValues(read, header, perReadCovariateMatrix, true);
        verifyCovariateArray(perReadCovariateMatrix.getMismatchMatrix(), readLength, -1);

        read.setIsSecondOfPair();
        covariate.recordValues(read, header, perReadCovariateMatrix, true);
        verifyCovariateArray(perReadCovariateMatrix.getMismatchMatrix(), -readLength, 1);

        read.setIsReverseStrand(false);
        covariate.recordValues(read, header, perReadCovariateMatrix, true);
        verifyCovariateArray(perReadCovariateMatrix.getMismatchMatrix(), -1, -1);
    }

    private void verifyCovariateArray(int[][] values, int init, int increment) {
        for (int i = 0; i < values.length; i++) {
            short actual = Short.decode(covariate.formatKey(values[i][0]));
            int expected = init + (increment * i);
            Assert.assertEquals(actual, expected);
        }
    }

    @Test(expectedExceptions={UserException.class})
    public void testMoreThanMaxCycleFails() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeaderWithReadGroup(illuminaReadGroup);

        int readLength = RAC.MAXIMUM_CYCLE_VALUE + 1;
        GATKRead read = ArtificialReadUtils.createRandomRead(readLength);
        read.setIsPaired(true);
        read.setReadGroup(illuminaReadGroup.getReadGroupId());

        PerReadCovariateMatrix perReadCovariateMatrix = new PerReadCovariateMatrix(read.getLength(), 1, new CovariateKeyCache());
        covariate.recordValues(read, header, perReadCovariateMatrix, true);
    }

    @Test
    public void testMaxCyclePasses() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeaderWithReadGroup(illuminaReadGroup);

        int readLength = RAC.MAXIMUM_CYCLE_VALUE;
        GATKRead read = ArtificialReadUtils.createRandomRead(readLength);
        read.setIsPaired(true);
        read.setReadGroup(illuminaReadGroup.getReadGroupId());

        PerReadCovariateMatrix perReadCovariateMatrix = new PerReadCovariateMatrix(read.getLength(), 1, new CovariateKeyCache());
        covariate.recordValues(read, header, perReadCovariateMatrix, true);
    }

    public static int expectedCycle(GATKRead read, final int baseNumber, final boolean indel, final int maxCycle) {
        return CycleCovariate.cycleFromKey(CycleCovariate.cycleKey(baseNumber, read, indel, maxCycle));
    }

    @DataProvider(name="cycleKey")
    public Object[][] cycleKey() {
        return new Object[][]{
                //positive strand
                {"10M", false, false, 9, 20, false, (1+9)<<1},
                {"10M", false, false, 1, 20, false, (1+1)<<1},
                {"10M", false, false, 5, 20, false, (1+5)<<1},
                {"10M", false, false, 9, 20, true,  -1}, //indels at the ends are -1
                {"10M", false, false, 1, 20, true,  -1}, //indels at the ends are -1
                {"10M", false, false, 5, 20, true,  (1+5)<<1},

                //positive strand, second in pair
                {"10M", false, true, 9, 20, false, (abs(-1-9)<<1)+1},
                {"10M", false, true, 1, 20, false, (abs(-1-1)<<1)+1},
                {"10M", false, true, 5, 20, false, (abs(-1-5)<<1)+1},
                {"10M", false, true, 9, 20, true,  -1}, //indels at the ends are -1
                {"10M", false, true, 1, 20, true,  -1}, //indels at the ends are -1
                {"10M", false, true, 5, 20, true,  (abs(-1-5)<<1)+1},

                //negative strand. Read is size 10
                {"10M", true, false, 9, 20, false, (10-9)<<1},
                {"10M", true, false, 1, 20, false, (10-1)<<1},
                {"10M", true, false, 5, 20, false, (10-5)<<1},
                {"10M", true, false, 9, 20, true,  -1}, //indels at the ends are -1
                {"10M", true, false, 1, 20, true,  -1}, //indels at the ends are -1
                {"10M", true, false, 5, 20, true,  (10-5)<<1},

                //negative strand, second in pair. Read is size 10
                {"10M", true, true, 9, 20, false, abs((-10 + 9) << 1) +1},
                {"10M", true, true, 1, 20, false, abs((-10+1)<<1)+1},
                {"10M", true, true, 5, 20, false, abs((-10 + 5) << 1)+1},
                {"10M", true, true, 9, 20, true,  -1}, //indels at the ends are -1
                {"10M", true, true, 1, 20, true,  -1}, //indels at the ends are -1
                {"10M", true, true, 5, 20, true,  abs((-10+5)<<1)+1},
        };
    }

    @Test(dataProvider = "cycleKey")
    public void testCycleKey(final String cigar, final boolean isNegStrand, final boolean isSecondInPair, final int baseNumber, final int maxCycle, final boolean indel, final int expectedCycleKey){
        final GATKRead read = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode(cigar));
        read.setIsReverseStrand(isNegStrand);
        if ( isSecondInPair ) {
            read.setIsSecondOfPair();
        }
        else {
            read.setIsFirstOfPair();
        }

        final int cycleKey = CycleCovariate.cycleKey(baseNumber, read, indel, maxCycle);
        Assert.assertEquals(expectedCycleKey, cycleKey);
    }

    @DataProvider(name="cycleFromKey")
    public Object[][] cycleFromKey() {
        return new Object[][]{
                {0b100, 0b10},
                {0b101, -0b10},
                {0b110, 0b11},
                {0b10,  0b1},
                {0b11, -0b1},
        };
    }

    @Test(dataProvider = "cycleFromKey")
    public void testCycleFromKey(final int key, final int expectedCycle){
        final int cycle = CycleCovariate.cycleFromKey(key);
        Assert.assertEquals(expectedCycle, cycle);
    }

    @Test
    public void testCycleToKeyAndBack(){
        final int max = 1000000;
        for (int cycle = -abs(max)+1; cycle <= abs(max)-1; cycle++) {
            final int key = CycleCovariate.keyFromCycle(cycle, max);
            final int cycleBack = CycleCovariate.cycleFromKey(key);
            Assert.assertEquals(cycle, cycleBack);
        }
    }

    @Test(expectedExceptions = UserException.class)
    public void testCycleTooLow(){
        final int max = 10;
        final int cycle = -abs(max)-1;
        final int key = CycleCovariate.keyFromCycle(cycle, max);
    }

    @Test(expectedExceptions = UserException.class)
    public void testCycleTooHigh(){
        final int max = 10;
        final int cycle = abs(max)+1;
        final int key = CycleCovariate.keyFromCycle(cycle, max);
    }

}
