package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.function.LongSupplier;

public class ProgressMeterUnitTest extends GATKBaseTest {

    private static class ListBasedTimeFunction implements LongSupplier {
        private final List<Long> values;
        private int currentIndex = 0;

        public ListBasedTimeFunction( final List<Long> values ) {
            this.values = values;
        }

        @Override
        public long getAsLong() {
            if ( currentIndex >= values.size() ) {
                throw new NoSuchElementException();
            }

            return values.get(currentIndex++);
        }

        public int size() {
            return values.size();
        }
    }

    private static class NullyLocatable implements Locatable {
        @Override
        public String getContig() {
            return null;
        }

        @Override
        public int getStart() {
            return 0;
        }

        @Override
        public int getEnd() {
            return 0;
        }
    }

    @DataProvider(name = "UpdateIntervalTestData")
    public Object[][] getUpdateIntervalTestData() {
        return new Object[][] {
                // Seconds between logger updates, time function, total number of records to process, and expected number of progress updates
                { 1.0, new ListBasedTimeFunction(Arrays.asList(1000l, 2000l)), ProgressMeter.DEFAULT_RECORDS_BETWEEN_TIME_CHECKS, 1 },
                { 1.0, new ListBasedTimeFunction(Arrays.asList(1000l, 2000l, 3000l)), ProgressMeter.DEFAULT_RECORDS_BETWEEN_TIME_CHECKS * 2, 2 },
                { 1.5, new ListBasedTimeFunction(Arrays.asList(1000l, 2000l)), ProgressMeter.DEFAULT_RECORDS_BETWEEN_TIME_CHECKS, 0 },
                { 2.0, new ListBasedTimeFunction(Arrays.asList(1000l, 2000l, 3000l, 4000l, 5000l)), ProgressMeter.DEFAULT_RECORDS_BETWEEN_TIME_CHECKS * 4, 2 },
                { 2.0, new ListBasedTimeFunction(Arrays.asList(1000l, 2000l, 3000l, 4000l)), ProgressMeter.DEFAULT_RECORDS_BETWEEN_TIME_CHECKS * 3, 1 },
                { 1.0, new ListBasedTimeFunction(Arrays.asList(1000l, 1500l, 2500l, 3500l)), ProgressMeter.DEFAULT_RECORDS_BETWEEN_TIME_CHECKS * 3, 2 },
                { 1.0, new ListBasedTimeFunction(Arrays.asList(1000l, 1500l, 2500l, 3400l)), ProgressMeter.DEFAULT_RECORDS_BETWEEN_TIME_CHECKS * 3, 1 },
        };
    }

    @Test(dataProvider = "UpdateIntervalTestData")
    public void testUpdateInterval( final double secondsBetweenUpdates, final ListBasedTimeFunction timeFunction, final long numRecords, final int expectedUpdates ) {
        final ProgressMeter meter = new ProgressMeter(secondsBetweenUpdates, timeFunction);
        meter.start();
        for ( int i = 1; i <= numRecords; ++i ) {
            meter.update(new SimpleInterval("1", 1, 1));
        }

        Assert.assertEquals(meter.numLoggerUpdates(), expectedUpdates, "Wrong number of logger updates given secondsBetweenUpdates = " + secondsBetweenUpdates);
    }

    @Test()
    public void testGetHeader() {
        // test with no percent remaining:
        final ProgressMeter m1 = new ProgressMeter(1.0, new ListBasedTimeFunction(Arrays.asList(1000L, 2000L)));
        Assert.assertEquals(m1.getHeader(), "       Current Locus  Elapsed Minutes     Records Processed   Records/Minute");

        final SAMSequenceDictionary sequenceDictionary =  new SAMSequenceDictionary(Arrays.asList(
                new SAMSequenceRecord("contig1", 1000),
                new SAMSequenceRecord("contig2", 1117),
                new SAMSequenceRecord("contig3", 1000)
        ));

        // test with percent remaining:
        final ProgressMeter m2 = new ProgressMeter(
                1.0,
                new ListBasedTimeFunction(Arrays.asList(1000L, 2000L)),
                sequenceDictionary
        );
        Assert.assertEquals(m2.getHeader(), "       Current Locus  Elapsed Minutes     Records Processed   Records/Minute  % Complete  Est. Time Remaining (s)");
    }

    @Test()
    public void testSetSequenceDictionary() {
        final ProgressMeter meter = new ProgressMeter(1.0, new ListBasedTimeFunction(Arrays.asList(1000L, 2000L)));
        Assert.assertNull(meter.sequenceDictionary);

        final SAMSequenceDictionary sequenceDictionary =  new SAMSequenceDictionary(Arrays.asList(
                new SAMSequenceRecord("contig1", 1000),
                new SAMSequenceRecord("contig2", 1117),
                new SAMSequenceRecord("contig3", 1000)
        ));

        // Validate that our sequence dictionary actually gets set:
        meter.setSequenceDictionary(sequenceDictionary);
        Assert.assertEquals(meter.sequenceDictionary, sequenceDictionary);

        // Validate that the cumulative bases in the dictionary are correct:
        Assert.assertEquals(meter.cumulativeBasesInRefDict, new LinkedHashMap<String, Long>() {{
            put("contig1", 0L);
            put("contig2", 1000L);
            put("contig3", 2117L);
        }});

        // Validate the above but with the constructor:
        final ProgressMeter meter2 = new ProgressMeter(1.0,
                new ListBasedTimeFunction(Arrays.asList(1000L, 2000L)),
                sequenceDictionary
        );

        Assert.assertEquals(meter2.cumulativeBasesInRefDict, new LinkedHashMap<String, Long>() {{
            put("contig1", 0L);
            put("contig2", 1000L);
            put("contig3", 2117L);
        }});
    }
    @DataProvider(name = "provideForCalculateEstimatedPercentComplete")
    public Object[][] provideForCalculateEstimatedPercentComplete() {

        final SAMSequenceDictionary sequenceDictionary =  new SAMSequenceDictionary(Arrays.asList(
                new SAMSequenceRecord("contig1", 2000),
                new SAMSequenceRecord("contig2", 3000),
                new SAMSequenceRecord("contig3", 5000)
        ));

        return new Object[][] {
                {sequenceDictionary, new SimpleInterval("contig1", 1000, 1000), 10.0},
                {sequenceDictionary, new SimpleInterval("contig1", 2000, 2000), 20.0},
                {sequenceDictionary, new SimpleInterval("contig2", 1000, 1000), 30.0},
                {sequenceDictionary, new SimpleInterval("contig2", 2000, 2000), 40.0},
                {sequenceDictionary, new NullyLocatable(), 0},
                {sequenceDictionary, null, 0},
        };
    }

    @Test(dataProvider = "provideForCalculateEstimatedPercentComplete")
    public void testCalculateEstimatedPercentComplete(final SAMSequenceDictionary sequenceDictionary,
                                                      final Locatable locus,
                                                      final double expectedPercent) {

        final double epsilon = 0.0001;

        final ProgressMeter meter = new ProgressMeter(1.0,
                new ListBasedTimeFunction(Arrays.asList(1000L, 2000L)),
                sequenceDictionary
        );

        meter.start();

        // We have to do at least DEFAULT_RECORDS_BETWEEN_TIME_CHECKS updates for the locus to be updated:
        for (int i = 0 ; i < ProgressMeter.DEFAULT_RECORDS_BETWEEN_TIME_CHECKS; ++i) {
            meter.update(locus);
        }

        Assert.assertTrue((meter.calculateEstimatedPercentComplete() - expectedPercent) < epsilon );
    }

    @DataProvider(name = "provideForTestCalculateEstimatedTimeRemaining_s")
    public Object[][] provideForTestCalculateEstimatedTimeRemaining_s() {

        final SAMSequenceDictionary sequenceDictionary =  new SAMSequenceDictionary(Arrays.asList(
                new SAMSequenceRecord("contig1", 2000),
                new SAMSequenceRecord("contig2", 3000),
                new SAMSequenceRecord("contig3", 5000)
        ));

        return new Object[][] {
                {sequenceDictionary, new SimpleInterval("contig1", 1000, 1000), 1000, 9},
                {sequenceDictionary, new SimpleInterval("contig1", 1000, 1000), 5000, 45},
                {sequenceDictionary, new SimpleInterval("contig1", 1000, 1000), 360000, 3240},
                {sequenceDictionary, new SimpleInterval("contig1", 2000, 2000), 360000, 1440},
                {sequenceDictionary, new SimpleInterval("contig2", 1000, 1000), 360000, 840},
                {sequenceDictionary, new SimpleInterval("contig2", 2000, 2000), 360000, 540},
        };
    }

    @Test(dataProvider = "provideForTestCalculateEstimatedTimeRemaining_s")
    public void testCalculateEstimatedTimeRemaining_s(final SAMSequenceDictionary sequenceDictionary,
                                                      final Locatable locus,
                                                      final long elapsedTime_ms,
                                                      final long expectedRemainingTime_s) {
        final ProgressMeter meter = new ProgressMeter(1.0,
                new ListBasedTimeFunction(Arrays.asList(1000L, 2000L)),
                sequenceDictionary
        );

        meter.start();

        // We have to do at least DEFAULT_RECORDS_BETWEEN_TIME_CHECKS updates for the locus to be updated:
        for (int i = 0 ; i < ProgressMeter.DEFAULT_RECORDS_BETWEEN_TIME_CHECKS; ++i) {
            meter.update(locus);
        }

        final double percentComplete = meter.calculateEstimatedPercentComplete();
        Assert.assertEquals((meter.calculateEstimatedTimeRemaining_s(percentComplete, elapsedTime_ms)), expectedRemainingTime_s );
    }


    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInvalidUpdateInterval() {
        final ProgressMeter meter = new ProgressMeter(0.0);
    }

    @DataProvider(name = "ElapsedTimeInMinutesTestData")
    public Object[][] getElapsedTimeInMinutesTestData() {
        return new Object[][] {
                // 600,000ms = 10 minutes
                { new ListBasedTimeFunction(Arrays.asList(1000l, 60000l * 10l + 1000l)), 10.0 },
                // 60,000ms = 1 minute
                { new ListBasedTimeFunction(Arrays.asList(1000l, 60000l + 1000l)), 1.0 },
                // 30,000ms = 0.5 minutes
                { new ListBasedTimeFunction(Arrays.asList(1000l, 30000l + 1000l)), 0.5 },
                // 0ms = 0 minutes
                { new ListBasedTimeFunction(Arrays.asList(1000l, 1000l)), 0.0 },
                // 10,000ms = 1.0/6.0 minutes
                { new ListBasedTimeFunction(Arrays.asList(1000l, 3000l, 6000l, 10000l, 11000l)), 1.0 / 6.0 }
        };
    }

    @Test(dataProvider = "ElapsedTimeInMinutesTestData")
    public void testElapsedTimeInMinutes( final ListBasedTimeFunction timeFunction, final double expectedElapsedMinutes ) {
        final ProgressMeter meter = new ProgressMeter(1l, timeFunction);

        meter.start();
        // The start() call consumes one value from the time function, so we need to process
        // (timeFunction.size() - 1) * ProgressMeter.RECORDS_BETWEEN_TIME_CHECKS in order
        // to consume the rest of the values in the time function
        for ( int i = 1; i <= (timeFunction.size() - 1) * ProgressMeter.DEFAULT_RECORDS_BETWEEN_TIME_CHECKS; ++i ) {
            meter.update(new SimpleInterval("1", 1, 1));
        }

        Assert.assertEquals(meter.elapsedTimeInMinutes(), expectedElapsedMinutes, "elapsed minutes differs from expected value");
    }

    @DataProvider(name = "ProcessingRateTestData")
    public Object[][] getProcessingRateTestData() {
        // We need to base the number of records we process on the time check interval (to ensure that
        // the progress meter pulls all values from the time function and gets updated properly).
        final long timeCheckRecords = ProgressMeter.DEFAULT_RECORDS_BETWEEN_TIME_CHECKS;
        return new Object[][] {
                // If we process timeCheckRecords records in 1 minute, rate should be timeCheckRecords
                { new ListBasedTimeFunction(Arrays.asList(1000l, 60000l + 1000l)), timeCheckRecords, timeCheckRecords },
                // If we process timeCheckRecords records in 2 minutes, rate should be timeCheckRecords / 2
                { new ListBasedTimeFunction(Arrays.asList(1000l, 60000l * 2l + 1000l)), timeCheckRecords, timeCheckRecords / 2},
                // timeCheckRecords records in 15 seconds should give timeCheckRecords * 4 records/minute
                { new ListBasedTimeFunction(Arrays.asList(1000l, 15000l + 1000l)), timeCheckRecords, timeCheckRecords * 4 }
        };
    }

    @Test(dataProvider = "ProcessingRateTestData")
    public void testProcessingRate( final ListBasedTimeFunction timeFunction, final long recordsToProcess, final double expectedProcessingRate ) {
        final ProgressMeter meter = new ProgressMeter(1.0, timeFunction);

        meter.start();
        for ( int i = 1; i <= recordsToProcess; ++i ) {
            meter.update(new SimpleInterval("1", 1, 1));
        }

        Assert.assertEquals(meter.processingRate(), expectedProcessingRate, "actual processing rate differs from expected value");
    }

    @Test
    public void testSecondsSinceLastPrint() {
        final ListBasedTimeFunction timeFunction = new ListBasedTimeFunction(Arrays.asList(1000l, 1500l, 2000l, 2500l, 3000l));
        final ProgressMeter meter = new ProgressMeter(1.0, timeFunction);
        final List<Double> expectedSecondsSinceLastPrint = Arrays.asList(0.5, 0.0, 0.5, 0.0);

        meter.start();
        // start() consumes one timestamp from our function, so we need to process size() - 1 more timestamps
        for ( int i = 1; i <= timeFunction.size() - 1; ++i ) {
            for ( int j = 1; j <= ProgressMeter.DEFAULT_RECORDS_BETWEEN_TIME_CHECKS; ++j ) {
                meter.update(new SimpleInterval("1", 1, 1));
            }

            Assert.assertEquals(meter.secondsSinceLastPrint(), (double) expectedSecondsSinceLastPrint.get(i - 1), "Wrong number of seconds reported since last print");
        }

    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testMustStartBeforeUpdate() throws Exception {
        ProgressMeter pm = new ProgressMeter();
        pm.update(new SimpleInterval("1",1,1));
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testCantStartTwice() throws Exception {
        ProgressMeter pm = new ProgressMeter();
        pm.start();
        pm.update(new SimpleInterval("1",1,1));
        pm.start();
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testCantStopBeforeStart() throws Exception {
        ProgressMeter pm = new ProgressMeter();
        pm.stop();
    }

    @Test
    public void testStartedAndStopped() throws Exception {
        ProgressMeter pm = new ProgressMeter();
        Assert.assertFalse(pm.started());
        Assert.assertFalse(pm.stopped());
        pm.start();
        Assert.assertTrue(pm.started());
        Assert.assertFalse(pm.stopped());
        pm.stop();
        Assert.assertTrue(pm.started());
        Assert.assertTrue(pm.stopped());
    }

    @Test
    public void testDisabledProgressMeter() {
        final ProgressMeter disabledPM = new ProgressMeter(ProgressMeter.DEFAULT_SECONDS_BETWEEN_UPDATES, true);

        disabledPM.start();
        Assert.assertFalse(disabledPM.started());
        disabledPM.update(new SimpleInterval("1", 1, 1));
        Assert.assertEquals(disabledPM.numRecordsProcessed(), 0);
        disabledPM.stop();
        Assert.assertFalse(disabledPM.stopped());
    }
}
