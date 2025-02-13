package org.broadinstitute.hellbender.engine;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SamConstants;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.function.LongSupplier;

/**
 * A basic progress meter to print out the number of records processed (and other metrics) during a traversal
 * at a configurable time interval.
 *
 * Clients set the update interval at construction, which controls how many seconds must elapse
 * before printing an update. Then call {@link #start} at traversal start, {@link #update(Locatable)}
 * after processing each record from the primary input, and {@link #stop} at traversal end to print
 * summary statistics.
 *
 * Note that {@link #start} must only be called once, before any {@link #update(Locatable)}.
 * Note no {@link #update(Locatable)} must be called after {@link #stop}.
 *
 * All output is made at INFO level via log4j.
 */
public final class ProgressMeter {
    protected static final Logger logger = LogManager.getLogger(ProgressMeter.class);

    /**
     * By default, we output a line to the logger after this many seconds have elapsed
     */
    public static final double DEFAULT_SECONDS_BETWEEN_UPDATES = 10.0;

    /**
     * We check the current time every time we process this many records,
     * by default (to cut down on the number of system calls)
     */
    public static final long DEFAULT_RECORDS_BETWEEN_TIME_CHECKS = 1000L;

    /**
     * By default, we use this function to get the current time
     */
    public static final LongSupplier DEFAULT_TIME_FUNCTION = System::currentTimeMillis;

    /**
     * Number of milliseconds in a second
     */
    public static final long MILLISECONDS_PER_SECOND = 1000L;

    /**
     * Number of milliseconds in a minute
     */
    public static final long MILLISECONDS_PER_MINUTE = MILLISECONDS_PER_SECOND * 60L;

    /**
     * Default label for records in logger messages. For display purposes only.
     */
    public static final String DEFAULT_RECORD_LABEL = "records";

    /**
     * We output a line to the logger after this many seconds have elapsed
     */
    private double secondsBetweenUpdates;

    /**
     * We check the current time every time we process this many records
     * (to cut down on the number of system calls)
     */
    private long recordsBetweenTimeChecks = DEFAULT_RECORDS_BETWEEN_TIME_CHECKS;

    /**
     * Total records processed
     */
    private long numRecordsProcessed = 0L;

    /**
     * Our start timestamp in milliseconds as returned by our {@link #timeFunction}
     */
    private long startTimeMs = 0L;

    /**
     * Current timestamp in milliseconds as returned by our {@link #timeFunction}.
     *
     * Updated only every {@link #recordsBetweenTimeChecks} records to cut down
     * on system calls.
     */
    private long currentTimeMs = 0L;

    /**
     * Timestamp in milliseconds as returned by our {@link #timeFunction} of the last time
     * we outputted a progress line to the logger
     */
    private long lastPrintTimeMs = 0L;

    /**
     * The genomic location of the most recently processed record, or null if the most recent record had no location.
     * Updated only when we actually output a line to the logger.
     */
    private Locatable currentLocus = null;

    /**
     * The number of times we've outputted a status line to the logger via {@link #printProgress}.
     * We keep track of this only for unit-testing purposes.
     */
    private long numLoggerUpdates = 0L;

    /**
     * Function that returns the current time in milliseconds (defaults to {@link #DEFAULT_TIME_FUNCTION}).
     * Should only be customized for unit testing purposes.
     */
    private LongSupplier timeFunction;

    /**
     * Cached copy of the sequence dictionary from the currently running tool.
     * Used to estimate percentage completion and time remaining.
     */
    @VisibleForTesting
     SAMSequenceDictionary sequenceDictionary;

    /**
     * Map containing the cumulative number of bases in each contig in the reference dictionary.
     * Used to estimate percentage completion and time remaining.
     */
    final HashMap<String, Long> cumulativeBasesInRefDict = new LinkedHashMap<>();

    /**
     * Keeps track of whether the progress meter has ever been started.
     */
    private boolean started;

    /**
     * Keeps track of whether the progress meter has been stopped.
     */
    private boolean stopped;

    /**
     * Keeps track of whether the progress meter has been disabled. When disabled,
     * operations on the ProgressMeter have no effect.
     */
    private final boolean disabled;

    /**
     * Label for records in logger messages. For display purposes only.
     */
    private String recordLabel = DEFAULT_RECORD_LABEL;

    /**
     * Create a progress meter with the default update interval of {@link #DEFAULT_SECONDS_BETWEEN_UPDATES} seconds
     * and the default time function {@link #DEFAULT_TIME_FUNCTION}.
     */
    public ProgressMeter() {
        this(DEFAULT_SECONDS_BETWEEN_UPDATES);
    }

    /**
     * Create a progress meter with a custom update interval and the default time function {@link #DEFAULT_TIME_FUNCTION}
     *
     * @param secondsBetweenUpdates number of seconds that should elapse before outputting a line to the logger
     */
    public ProgressMeter( final double secondsBetweenUpdates ) {
        this(secondsBetweenUpdates, DEFAULT_TIME_FUNCTION, false, null);
    }

    /**
     * Create a progress meter with a custom update interval and the default time function {@link #DEFAULT_TIME_FUNCTION},
     * and allows the ProgressMeter to be set as disabled
     *
     * @param secondsBetweenUpdates number of seconds that should elapse before outputting a line to the logger
     * @param disabled if true, set this ProgressMeter to be disabled, so that all operations on it are no-ops
     */
    public ProgressMeter( final double secondsBetweenUpdates, final boolean disabled ) {
        this(secondsBetweenUpdates, DEFAULT_TIME_FUNCTION, disabled, null);
    }

    /**
     * Create a progress meter with a custom update interval and the default time function {@link #DEFAULT_TIME_FUNCTION},
     * and allows the ProgressMeter to be set as disabled
     *
     * @param secondsBetweenUpdates number of seconds that should elapse before outputting a line to the logger
     * @param disabled if true, set this ProgressMeter to be disabled, so that all operations on it are no-ops
     */
    public ProgressMeter(final double secondsBetweenUpdates, final boolean disabled, final SAMSequenceDictionary sequenceDictionary) {
        this(secondsBetweenUpdates, DEFAULT_TIME_FUNCTION, disabled, sequenceDictionary);
    }


    /**
     * Create a progress meter with a custom update interval and a custom function for getting the current
     * time in milliseconds.
     *
     * Providing your own time function is only useful in unit tests -- in normal usage
     * clients should call one of the other constructors.
     *
     * @param secondsBetweenUpdates number of seconds that should elapse before outputting a line to the logger
     * @param timeFunction function that returns the current time in milliseconds.
     */
    @VisibleForTesting
    ProgressMeter( final double secondsBetweenUpdates, final LongSupplier timeFunction ) {
        this(secondsBetweenUpdates, timeFunction, false, null);
    }

    /**
     * Create a progress meter with a custom update interval and a custom function for getting the current
     * time in milliseconds.
     *
     * Providing your own time function is only useful in unit tests -- in normal usage
     * clients should call one of the other constructors.
     *
     * @param secondsBetweenUpdates number of seconds that should elapse before outputting a line to the logger
     * @param timeFunction function that returns the current time in milliseconds.
     * @param sequenceDictionary the {@link SAMSequenceDictionary} for the current tool.  Used to estimate % complete and remaining time.
     */
    @VisibleForTesting
    ProgressMeter(final double secondsBetweenUpdates, final LongSupplier timeFunction, final SAMSequenceDictionary sequenceDictionary) {
        this(secondsBetweenUpdates, timeFunction, false, sequenceDictionary);
    }

    /**
     * Create a progress meter with a custom update interval and a custom function for getting the current
     * time in milliseconds.
     *
     * Providing your own time function is only useful in unit tests -- in normal usage
     * clients should call one of the other constructors.
     *
     * @param secondsBetweenUpdates number of seconds that should elapse before outputting a line to the logger
     * @param timeFunction function that returns the current time in milliseconds.
     * @param disabled if true, set this ProgressMeter to be disabled, so that all operations on it are no-ops
     * @param sequenceDictionary the {@link SAMSequenceDictionary} for the current tool.  Used to estimate % complete and remaining time.
     */
    @VisibleForTesting
    ProgressMeter( final double secondsBetweenUpdates, final LongSupplier timeFunction, final boolean disabled, final SAMSequenceDictionary sequenceDictionary ) {
        Utils.nonNull(timeFunction);
        Utils.validateArg(secondsBetweenUpdates > 0, "secondsBetweenUpdates must be > 0.0");
        this.started = false;
        this.stopped = false;
        this.disabled = disabled;
        this.secondsBetweenUpdates = secondsBetweenUpdates;
        this.timeFunction = timeFunction;

        // If we have a sequence dictionary, process it for estimated time remaining and percentage of progress:
        if (sequenceDictionary != null) {
            setSequenceDictionary(sequenceDictionary);
        }
    }

    /**
     * Set the sequence dictionary for this progress meter and allow it to estimate time remaining.
     * @param sequenceDictionary the {@link SAMSequenceDictionary} for the current tool.  Used to estimate % complete and remaining time.
     */
    public void setSequenceDictionary(final SAMSequenceDictionary sequenceDictionary) {

        // Process the sequence dictionary so that we can get estimated progress:
        this.sequenceDictionary = sequenceDictionary;

        // Generate the cumulative number of bases in each contig from the reference dict:
        long cSum = 0;
        for ( final SAMSequenceRecord r : sequenceDictionary.getSequences() ) {
            cumulativeBasesInRefDict.put(r.getSequenceName(), cSum);
            cSum += r.getSequenceLength();
        }
    }

    /**
     * Set the number of records we need to process before we check the current time
     *
     * @param recordsBetweenTimeChecks number of records we need to process before we check the current time
     */
    public void setRecordsBetweenTimeChecks( final long recordsBetweenTimeChecks ) {
        this.recordsBetweenTimeChecks = recordsBetweenTimeChecks;
    }

    /**
     * Change the label used for records in logger messages. Default label is {@link #DEFAULT_RECORD_LABEL}
     *
     * @param label Label to use for records in logger messages. Not null.
     */
    public void setRecordLabel( final String label ) {
        Utils.nonNull(label);
        this.recordLabel = label;
    }

    /**
     * Start the progress meter and produce preliminary output such as column headings.
     * @throws IllegalStateException if the meter has been started before or has been stopped already
     */
    public void start() {
        if ( disabled ) {
            return;
        }

        Utils.validate( !started, "the progress meter has been started already");
        Utils.validate( !stopped, "the progress meter has been stopped already");
        started = true;
        logger.info("Starting traversal");
        printHeader();

        startTimeMs = timeFunction.getAsLong();
        currentTimeMs = startTimeMs;
        lastPrintTimeMs = startTimeMs;
        numRecordsProcessed = 0L;
        numLoggerUpdates = 0L;
        currentLocus = null;
    }

    /**
     * Signal to the progress meter that an additional record has been processed. Will output
     * statistics to the logger roughly every {@link #secondsBetweenUpdates} seconds.
     *
     * @param currentLocus the genomic location of the record just processed or null if the most recent record had no location.
     * @throws IllegalStateException if the meter has not been started yet or has been stopped already
     */
    public void update( final Locatable currentLocus ) {
        if ( disabled ) {
            return;
        }

        Utils.validate(started, "the progress meter has not been started yet");
        Utils.validate( !stopped, "the progress meter has been stopped already");
        ++numRecordsProcessed;
        if ( numRecordsProcessed % recordsBetweenTimeChecks == 0 ) {
            currentTimeMs = timeFunction.getAsLong();
            this.currentLocus = currentLocus;

            if ( secondsSinceLastPrint() >= secondsBetweenUpdates ) {
                printProgress();
                lastPrintTimeMs = currentTimeMs;
            }
        }
    }

    /**
     * Signal to the progress meter that several additional records hav been processed. Will output
     * statistics to the logger roughly every {@link #secondsBetweenUpdates} seconds.
     * @param currentLocus the genomic location of the last record just processed or null if the most recent
     *                     record had no location.
     * @param recordCountIncrease number of new records processed.
     * @throws IllegalStateException if the meter has not been started yet or has been stopped already
     */
    public void update( final Locatable currentLocus, final long recordCountIncrease ) {
        if ( disabled ) {
            return;
        }

        Utils.validate(started, "the progress meter has not been started yet");
        Utils.validate( !stopped, "the progress meter has been stopped already");
        final long previousNumRecordsProcessed = numRecordsProcessed;
        numRecordsProcessed += recordCountIncrease;
        if ( previousNumRecordsProcessed / recordsBetweenTimeChecks != (numRecordsProcessed / recordsBetweenTimeChecks) ) {
            currentTimeMs = timeFunction.getAsLong();
            this.currentLocus = currentLocus;
            if ( secondsSinceLastPrint() >= secondsBetweenUpdates ) {
                printProgress();
                lastPrintTimeMs = currentTimeMs;
            }
        }
    }

    /**
     * Stop the progress meter and output summary statistics to the logger
     * @throws IllegalStateException if the meter has not been started yet or has been stopped already
     */
    public void stop() {
        if ( disabled ) {
            return;
        }

        Utils.validate(started, "the progress meter has not been started yet");
        Utils.validate( !stopped, "the progress meter has been stopped already");
        this.stopped = true;
        currentTimeMs = timeFunction.getAsLong();
        // Output progress a final time at the end
        printProgress();
        logger.info(String.format("Traversal complete. Processed %d total %s in %.1f minutes.", numRecordsProcessed, recordLabel, elapsedTimeInMinutes()));
    }

    /**
     * Print column headings labelling the output from {@link #printProgress}
     */
    private void printHeader() {
        logger.info(getHeader());
    }

    @VisibleForTesting
    String getHeader() {
        if (sequenceDictionary != null) {
            return String.format("%20s  %15s  %20s  %15s  %10s  %23s",
                    "Current Locus", "Elapsed Minutes",
                    StringUtils.capitalize(recordLabel) + " Processed",
                    StringUtils.capitalize(recordLabel) + "/Minute",
                    "% Complete",
                    "Est. Time Remaining (s)"
            );
        }
        else{
            return String.format("%20s  %15s  %20s  %15s",
                    "Current Locus", "Elapsed Minutes",
                    StringUtils.capitalize(recordLabel) + " Processed",
                    StringUtils.capitalize(recordLabel) + "/Minute"
            );
        }
    }

    /**
     * Output traversal statistics to the logger.
     */
    private void printProgress() {
        ++numLoggerUpdates;

        // If we can display percentage of time remaining, we do so here:
        if (sequenceDictionary != null) {

            // We must check if the contig for the current locus is set:
            if ((currentLocus == null) || (currentLocus.getContig() == null)) {
                logger.info(String.format("%20s  %15.1f  %20d  %15.1f  %12s  %23s",
                        currentLocusString(), elapsedTimeInMinutes(), numRecordsProcessed, processingRate(),
                        "UNKNOWN", "UNKNOWN"));
            }
            else {
                double percentageComplete = calculateEstimatedPercentComplete();
                long   timeRemaining_s    = calculateEstimatedTimeRemaining_s(percentageComplete, currentTimeMs - startTimeMs);

                logger.info(String.format("%20s  %15.1f  %20d  %15.1f  %10.2f  %23d",
                        currentLocusString(), elapsedTimeInMinutes(), numRecordsProcessed, processingRate(),
                        percentageComplete, timeRemaining_s));
            }
        }
        else {
            logger.info(String.format("%20s  %15.1f  %20d  %15.1f",
                    currentLocusString(), elapsedTimeInMinutes(), numRecordsProcessed, processingRate()));
        }
    }

    @VisibleForTesting
    double calculateEstimatedPercentComplete() {
        if ((currentLocus == null) || (currentLocus.getContig() == null)) {
            return 0;
        }

        final long currentOverallPos = cumulativeBasesInRefDict.get(currentLocus.getContig()) + currentLocus.getStart();
        return (100.0 * currentOverallPos) / sequenceDictionary.getReferenceLength();
    }

    @VisibleForTesting
    long calculateEstimatedTimeRemaining_s(final double percentComplete, final long elapsedTime_ms) {

        // 100 in numerator to cancel out the percentage units and get a proportion:
        double totalTime_ms = (100 * elapsedTime_ms) / percentComplete;

        // Subtract percentage complete from 100 to get percent remaining.
        // Divide by 100 to convert percentage to proportion.
        return (long)((((100.0 - percentComplete) * totalTime_ms)/100.0) / MILLISECONDS_PER_SECOND );
    }

    /**
     * @return the total minutes elapsed since we called {@link #start}
     *
     * This is only accurate at set polling intervals and should not be
     * called directly except in tests.
     */
    @VisibleForTesting
    double elapsedTimeInMinutes() {
        return (currentTimeMs - startTimeMs) / (double)MILLISECONDS_PER_MINUTE;
    }

    /**
     * @return the number of seconds that have elapsed since our last progress output to the logger
     *
     * This is only accurate at set polling intervals and should not be
     * called directly except in tests.
     */
    @VisibleForTesting
    double secondsSinceLastPrint() {
        return (currentTimeMs - lastPrintTimeMs) / (double)MILLISECONDS_PER_SECOND;
    }

    /**
     * @return number of records we're processing per minute, on average
     *
     * This is only accurate at set polling intervals and should not be
     * called directly except in tests.
     */
    @VisibleForTesting
    double processingRate() {
        return numRecordsProcessed / elapsedTimeInMinutes();
    }

    /**
     * @return number of times we've outputted a progress line to the logger (for unit testing purposes)
     */
    @VisibleForTesting
    long numLoggerUpdates() {
        return numLoggerUpdates;
    }

    /**
     * @return total number of records processed so far
     */
    @VisibleForTesting
    long numRecordsProcessed() {
        return numRecordsProcessed;
    }

    /**
     * @return genomic location of the most recent record formatted for output to the logger
     */
    private String currentLocusString() {
        return currentLocus != null ? currentLocus.getContig() + ":" + currentLocus.getStart() :
                                      "unmapped";
    }

    /**
     * Returns whether the meter has been started. It returns false before the call to {@link #start} and true forever after.
     */
    public boolean started() {
        return started;
    }

    /**
     * Returns whether the meter has been stopped. It returns false before the call to {@link #stop} and true forever after.
     */
    public boolean stopped() {
        return stopped;
    }
}
