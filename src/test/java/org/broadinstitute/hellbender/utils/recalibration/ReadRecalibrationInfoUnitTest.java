package org.broadinstitute.hellbender.utils.recalibration;

import htsjdk.samtools.SAMUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.recalibration.covariates.CovariateKeyCache;
import org.broadinstitute.hellbender.utils.recalibration.covariates.PerReadCovariateMatrix;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumMap;
import java.util.List;

public final class ReadRecalibrationInfoUnitTest extends GATKBaseTest {

    @DataProvider(name = "InfoProvider")
    public Object[][] createCombineTablesProvider() {
        List<Object[]> tests = new ArrayList<>();

        for ( final int readLength: Arrays.asList(10, 100, 1000) ) {
            for ( final boolean includeIndelErrors : Arrays.asList(true, false) ) {
                tests.add(new Object[]{readLength, includeIndelErrors});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "InfoProvider")
    public void testReadInfo(final int readLength, final boolean includeIndelErrors) {
        final PerReadCovariateMatrix covariates = new PerReadCovariateMatrix(readLength, 2, new CovariateKeyCache());

        final byte[] bases = new byte[readLength];
        final byte[] baseQuals = new byte[readLength];
        final byte[] insertionQuals = new byte[readLength];
        final byte[] deletionQuals = new byte[readLength];
        final boolean[] skips = new boolean[readLength];
        final double[] snpErrors = new double[readLength];
        final double[] insertionErrors = new double[readLength];
        final double[] deletionsErrors = new double[readLength];
        for ( int i = 0; i < readLength; i++ ) {
            bases[i] = 'A';
            baseQuals[i] = (byte)(i % SAMUtils.MAX_PHRED_SCORE);
            insertionQuals[i] = (byte)((i+1) % SAMUtils.MAX_PHRED_SCORE);
            deletionQuals[i] = (byte)((i+2) % SAMUtils.MAX_PHRED_SCORE);
            skips[i] = i % 2 == 0;
            snpErrors[i] = 1.0 / (i+1);
            insertionErrors[i] = 0.5 / (i+1);
            deletionsErrors[i] = 0.3 / (i+1);
        }

        final EnumMap<EventType, double[]> errors = new EnumMap<>(EventType.class);
        errors.put(EventType.BASE_SUBSTITUTION, snpErrors);
        errors.put(EventType.BASE_INSERTION, insertionErrors);
        errors.put(EventType.BASE_DELETION, deletionsErrors);

        final EnumMap<EventType, byte[]> quals = new EnumMap<>(EventType.class);
        quals.put(EventType.BASE_SUBSTITUTION, baseQuals);
        quals.put(EventType.BASE_INSERTION, insertionQuals);
        quals.put(EventType.BASE_DELETION, deletionQuals);

        final GATKRead read = ArtificialReadUtils.createArtificialRead(bases, baseQuals, readLength + "M");
        if ( includeIndelErrors ) {
            ReadUtils.setInsertionBaseQualities(read, insertionQuals);
            ReadUtils.setDeletionBaseQualities(read, deletionQuals);
        }

        final ReadRecalibrationInfo info = new ReadRecalibrationInfo(read, covariates, skips, snpErrors, insertionErrors, deletionsErrors);

        Assert.assertEquals(info.getCovariatesValues(), covariates);
        Assert.assertEquals(info.getRead(), read);

        for ( int i = 0; i < readLength; i++ ) {
            Assert.assertEquals(info.skip(i), skips[i]);
            for ( final EventType et : EventType.values() ) {
                Assert.assertEquals(info.getErrorFraction(et, i), errors.get(et)[i]);
                final byte expectedQual = et == EventType.BASE_SUBSTITUTION || includeIndelErrors ? quals.get(et)[i]: ReadUtils.DEFAULT_INSERTION_DELETION_QUAL;
                Assert.assertEquals(info.getQual(et, i), expectedQual);
            }
        }
    }
}
