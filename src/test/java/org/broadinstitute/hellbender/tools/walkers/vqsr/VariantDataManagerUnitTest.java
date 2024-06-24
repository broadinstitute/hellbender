package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.GATKBaseTest;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class VariantDataManagerUnitTest extends GATKBaseTest {

    @Test
    public final void testCalculateSortOrder() {
        final double passingQual = 400.0;
        final VariantRecalibratorArgumentCollection VRAC = new VariantRecalibratorArgumentCollection();

        VariantDataManager vdm = new VariantDataManager(new ArrayList<String>(), VRAC);

        final List<VariantDatum> theData = new ArrayList<>();
        final VariantDatum datum1 = new VariantDatum();
        datum1.atTrainingSite = true;
        datum1.failingSTDThreshold = false;
        datum1.originalQual = passingQual;
        datum1.annotations = new double[]{0.0,-10.0,10.0};
        datum1.isNull = new boolean[]{false, false, false};
        theData.add(datum1);

        final VariantDatum datum2 = new VariantDatum();
        datum2.atTrainingSite = true;
        datum2.failingSTDThreshold = false;
        datum2.originalQual = passingQual;
        datum2.annotations = new double[]{0.0,-9.0,15.0};
        datum2.isNull = new boolean[]{false, false, false};
        theData.add(datum2);

        final VariantDatum datum3 = new VariantDatum();
        datum3.atTrainingSite = false;
        datum3.failingSTDThreshold = false;
        datum3.originalQual = passingQual;
        datum3.annotations = new double[]{0.0,1.0,999.0};
        datum3.isNull = new boolean[]{false, false, false};
        theData.add(datum3);

        final VariantDatum datum4 = new VariantDatum();
        datum4.atTrainingSite = false;
        datum4.failingSTDThreshold = false;
        datum4.originalQual = passingQual;
        datum4.annotations = new double[]{0.015,2.0,1001.11};
        datum4.isNull = new boolean[]{false, false, false};
        theData.add(datum4);

        vdm.setData(theData);

        final double[] meanVector = new double[3];
        for( int iii = 0; iii < meanVector.length; iii++ ) {
            meanVector[iii] = vdm.mean(iii, true);
        }
        final List<Integer> order = vdm.calculateSortOrder(meanVector);
        Assert.assertTrue(Arrays.equals(new int[]{2,1,0}, ArrayUtils.toPrimitive(order.toArray(new Integer[order.size()]))));
    }

    @Test
    public final void testDownSamplingTrainingData() {
        final int MAX_NUM_TRAINING_DATA = 5000;
        final double passingQual = 400.0;
        final VariantRecalibratorArgumentCollection VRAC = new VariantRecalibratorArgumentCollection();
        VRAC.MAX_NUM_TRAINING_DATA = MAX_NUM_TRAINING_DATA;

        VariantDataManager vdm = new VariantDataManager(new ArrayList<String>(), VRAC);
        final List<VariantDatum> theData = new ArrayList<>();
        for( int iii = 0; iii < MAX_NUM_TRAINING_DATA * 10; iii++) {
            final VariantDatum datum = new VariantDatum();
            datum.atTrainingSite = true;
            datum.failingSTDThreshold = false;
            datum.originalQual = passingQual;
            theData.add(datum);
        }

        for( int iii = 0; iii < MAX_NUM_TRAINING_DATA * 2; iii++) {
            final VariantDatum datum = new VariantDatum();
            datum.atTrainingSite = false;
            datum.failingSTDThreshold = false;
            datum.originalQual = passingQual;
            theData.add(datum);
        }

        vdm.setData(theData);
        final List<VariantDatum> trainingData = vdm.getTrainingData();

        Assert.assertTrue( trainingData.size() == MAX_NUM_TRAINING_DATA );
    }

    @Test
    public final void testDropAggregateData() {
        final int MAX_NUM_TRAINING_DATA = 5000;
        final double passingQual = 400.0;
        final VariantRecalibratorArgumentCollection VRAC = new VariantRecalibratorArgumentCollection();
        VRAC.MAX_NUM_TRAINING_DATA = MAX_NUM_TRAINING_DATA;

        VariantDataManager vdm = new VariantDataManager(new ArrayList<String>(), VRAC);
        final List<VariantDatum> theData = new ArrayList<>();
        for( int iii = 0; iii < MAX_NUM_TRAINING_DATA * 10; iii++) {
            final VariantDatum datum = new VariantDatum();
            datum.atTrainingSite = true;
            datum.isAggregate = false;
            datum.failingSTDThreshold = false;
            datum.originalQual = passingQual;
            theData.add(datum);
        }

        for( int iii = 0; iii < MAX_NUM_TRAINING_DATA * 2; iii++) {
            final VariantDatum datum = new VariantDatum();
            datum.atTrainingSite = false;
            datum.isAggregate = true;
            datum.failingSTDThreshold = false;
            datum.originalQual = passingQual;
            theData.add(datum);
        }

        vdm.setData(theData);
        vdm.dropAggregateData();

        for( final VariantDatum datum : vdm.getData() ) {
            Assert.assertFalse( datum.isAggregate );
        }
    }
}
