package org.broadinstitute.hellbender.utils.recalibration.covariates;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.recalibration.RecalUtils;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Random;

public final class PerReadCovariateMatrixUnitTest extends GATKBaseTest {

    @Test
    public void testCovariateGeneration() {
        final RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();

        final String[] readGroups = {"RG1", "RG2", "RGbla"};
        ReadGroupCovariate rgCov = new ReadGroupCovariate(Arrays.asList(readGroups));
        QualityScoreCovariate qsCov = new QualityScoreCovariate(RAC);
        ContextCovariate coCov = new ContextCovariate(RAC);
        CycleCovariate cyCov = new CycleCovariate(RAC);

        StandardCovariateList covariates = new StandardCovariateList(RAC, Arrays.asList(readGroups));

        final int NUM_READS = 100;
        final Random rnd = Utils.getRandomGenerator();
        final CovariateKeyCache keyCache= new CovariateKeyCache();

        for (int idx = 0; idx < NUM_READS; idx++) {
            for (final String readGroupID : readGroups) {
                final SAMReadGroupRecord readGroupRecord = new SAMReadGroupRecord(readGroupID);
                readGroupRecord.setPlatform("illumina");
                final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeaderWithReadGroup(readGroupRecord);

                final int length = 10 + rnd.nextInt(100); // random read length, at least 10 bp long
                final GATKRead read = ArtificialReadUtils.createRandomRead(header, length, false);
                read.setIsReverseStrand(rnd.nextBoolean());
                read.setReadGroup(readGroupID);

                final byte[] mQuals = read.getBaseQualities();
                final byte[] iQuals = ReadUtils.getBaseInsertionQualities(read);
                final byte[] dQuals = ReadUtils.getBaseDeletionQualities(read);
                PerReadCovariateMatrix rc = RecalUtils.computeCovariates(read, header, covariates, true, keyCache);

                // check that the length is correct
                Assert.assertEquals(rc.getMismatchMatrix().length, length);
                Assert.assertEquals(rc.getInsertionMatrix().length, length);
                Assert.assertEquals(rc.getDeletionMatrix().length, length);

                for (int i = 0; i < length; i++) {
                    // check that read group is always the same
                    Assert.assertEquals(rgCov.formatKey(rc.getMismatchCovariatesAtOffset(i)[0]), readGroupID);
                    Assert.assertEquals(rgCov.formatKey(rc.getInsertionCovariatesAtOffset(i)[0]), readGroupID);
                    Assert.assertEquals(rgCov.formatKey(rc.getDeletionCovariatesAtOffset(i)[0]), readGroupID);

                    // check quality score
                    Assert.assertEquals(qsCov.formatKey(rc.getMismatchCovariatesAtOffset(i)[1]), String.valueOf(mQuals[i]));
                    Assert.assertEquals(qsCov.formatKey(rc.getInsertionCovariatesAtOffset(i)[1]), String.valueOf(iQuals[i]));
                    Assert.assertEquals(qsCov.formatKey(rc.getDeletionCovariatesAtOffset(i)[1]),  String.valueOf(dQuals[i]));

                    // check context
                    Assert.assertEquals(coCov.formatKey(rc.getMismatchCovariatesAtOffset(i)[2]), ContextCovariateUnitTest.expectedContext(read, i, RAC.MISMATCHES_CONTEXT_SIZE, RAC.LOW_QUAL_TAIL), "read: " +idx  + " readGroup:" + readGroupID + " context mismatch key at position:" + i);
                    Assert.assertEquals(coCov.formatKey(rc.getInsertionCovariatesAtOffset(i)[2]), ContextCovariateUnitTest.expectedContext(read, i, RAC.INDELS_CONTEXT_SIZE, RAC.LOW_QUAL_TAIL), "read: " +idx  + " readGroup:" + readGroupID + " context insertion key at position:" + i);
                    Assert.assertEquals(coCov.formatKey(rc.getDeletionCovariatesAtOffset(i)[2]),  ContextCovariateUnitTest.expectedContext(read, i, RAC.INDELS_CONTEXT_SIZE, RAC.LOW_QUAL_TAIL), "read: " +idx  + " readGroup:" + readGroupID + " context deletion key at position:" + i);

                    // check cycle
                    final int expectedCycleMismatch = CycleCovariateUnitTest.expectedCycle(read, i, false, RAC.MAXIMUM_CYCLE_VALUE);
                    final int expectedCycleIndel = CycleCovariateUnitTest.expectedCycle(read, i, true, RAC.MAXIMUM_CYCLE_VALUE);
                    Assert.assertEquals(cyCov.formatKey(rc.getMismatchCovariatesAtOffset(i)[3]), String.valueOf(expectedCycleMismatch), "read: " + idx + " cycle mismatch key at position:" + i);
                    Assert.assertEquals(cyCov.formatKey(rc.getInsertionCovariatesAtOffset(i)[3]), String.valueOf(expectedCycleIndel),  "read: " +idx + " cycle insertion key at position:" + i);
                    Assert.assertEquals(cyCov.formatKey(rc.getDeletionCovariatesAtOffset(i)[3]), String.valueOf(expectedCycleIndel),  "read: " +idx + " cycle deletion key at position:" + i);
                }

            }

        }

    }

}
