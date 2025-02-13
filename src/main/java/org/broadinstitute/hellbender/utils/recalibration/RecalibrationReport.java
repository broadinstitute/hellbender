package org.broadinstitute.hellbender.utils.recalibration;


import org.apache.commons.collections.CollectionUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.collections.NestedIntegerArray;
import org.broadinstitute.hellbender.utils.recalibration.covariates.Covariate;
import org.broadinstitute.hellbender.utils.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.utils.report.GATKReport;
import org.broadinstitute.hellbender.utils.report.GATKReportTable;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * This class has all the static functionality for reading a recalibration report file into memory. 
 */
public final class RecalibrationReport {
    private static final Logger logger = LogManager.getLogger(RecalibrationReport.class);
    private QuantizationInfo quantizationInfo; // histogram containing the counts for qual quantization (calculated after recalibration is done)
    private final RecalibrationTables recalibrationTables; // quick access reference to the tables
    private final StandardCovariateList covariates; // list of all covariates to be used in this calculation

    private final GATKReportTable argumentTable; // keep the argument table untouched just for output purposes
    private final RecalibrationArgumentCollection RAC; // necessary for quantizing qualities with the same parameter

    public RecalibrationReport(final File recalFile) {
        this(new GATKReport(recalFile));
    }

    public RecalibrationReport(final InputStream recalibrationTableStream){
        this(new GATKReport(recalibrationTableStream));
    }

    public RecalibrationReport(final GATKReport report){
        this(report, report.getReadGroups());
    }

    public RecalibrationReport(final GATKReport report, final SortedSet<String> allReadGroups) {
        argumentTable = report.getTable(RecalUtils.ARGUMENT_REPORT_TABLE_TITLE);
        RAC = initializeArgumentCollectionTable(argumentTable);

        final GATKReportTable quantizedTable = report.getTable(RecalUtils.QUANTIZED_REPORT_TABLE_TITLE);
        quantizationInfo = initializeQuantizationTable(quantizedTable);

        covariates = new StandardCovariateList(RAC, new ArrayList<>(allReadGroups));

        recalibrationTables = new RecalibrationTables(covariates, allReadGroups.size());

        initializeReadGroupCovariates(allReadGroups);

        parseReadGroupTable(report.getTable(RecalUtils.READGROUP_REPORT_TABLE_TITLE), recalibrationTables.getReadGroupTable());

        parseQualityScoreTable(report.getTable(RecalUtils.QUALITY_SCORE_REPORT_TABLE_TITLE), recalibrationTables.getQualityScoreTable());

        parseAllCovariatesTable(report.getTable(RecalUtils.ALL_COVARIATES_REPORT_TABLE_TITLE), recalibrationTables);

    }

    /**
     * Gather multiple {@link RecalibrationReport}s into a single file
     * @param inputs a list of {@link RecalibrationReport} files to gather
     * @param output a file to write the recalibration reports to
     */
    public static void gatherReportsIntoOneFile(final List<File> inputs, final File output) {
        Utils.nonNull(inputs, "inputs");
        Utils.nonNull(output, "output");
        try (final PrintStream outputFile = new PrintStream(output)){
            final GATKReport report = gatherReports(inputs);
            report.print(outputFile);
        } catch(final FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(output, e);
        }
    }

    /**
     * Gathers a set of files containing {@link RecalibrationReport}s into a single {@link GATKReport}.
     *
     * @param inputs a list of files containing {@link RecalibrationReport}s
     * @return gathered recalibration GATK report
     */
    public static GATKReport gatherReports(final List<File> inputs) {
        Utils.nonNull(inputs);
        Utils.nonEmpty(inputs, "Cannot gather an empty list of inputs");

        final SortedSet<String> allReadGroups = new TreeSet<>();
        final Map<File, Set<String>> inputReadGroups = new LinkedHashMap<>();

        // Get the read groups from each input report
        for (final File input : inputs) {
            final GATKReport report = new GATKReport(input);
            final Set<String> readGroups = report.getReadGroups();
            inputReadGroups.put(input, readGroups);
            allReadGroups.addAll(readGroups);
        }

        logTablesWithMissingReadGroups(allReadGroups, inputReadGroups);

        final RecalibrationReport result = inputs.stream()
                .map(i -> new RecalibrationReport(new GATKReport(i), allReadGroups))
                .reduce(RecalibrationReport::combine)
                .filter(r -> !r.isEmpty())
                .orElseThrow(() -> new GATKException("there is no usable data in any input file") );

        result.quantizationInfo = new QuantizationInfo(result.recalibrationTables, result.RAC.QUANTIZING_LEVELS);
        return result.createGATKReport();
    }

    /**
     * helper function to output log messages if there are tables that are missing read groups that were seen in the other tables
     * @param allReadGroups a set of all of the read groups across inputs
     * @param inputReadGroups a map from file to the read groups that file contains
     */
    private static void logTablesWithMissingReadGroups(SortedSet<String> allReadGroups, Map<File, Set<String>> inputReadGroups) {
        // Log the read groups that are missing from specific inputs
        for (final Map.Entry<File, Set<String>> entry: inputReadGroups.entrySet()) {
            final File input = entry.getKey();
            final Set<String> readGroups = entry.getValue();
            if (allReadGroups.size() != readGroups.size()) {
                // Since this is not completely unexpected, more than debug, but less than a proper warning.
                logger.info("Missing read group(s)" + ": " + input.getAbsolutePath());
                for (final Object readGroup: CollectionUtils.subtract(allReadGroups, readGroups)) {
                    logger.info("  " + readGroup);
                }
            }
        }
    }

    /**
    * Combines two recalibration reports by adding all observations and errors
    *
    * Note: This method DOES NOT recalculate the empirical qualities and quantized qualities. You have to recalculate
    * them after combining. The reason for not calculating it is because this function is intended for combining a
    * series of recalibration reports, and it only makes sense to calculate the empirical qualities and quantized
    * qualities after all the recalibration reports have been combined. Having the user recalculate when appropriate,
    * makes this method faster
    *
    * Note2: The empirical quality reported, however, is recalculated given its simplicity.
    *
    * @param other the recalibration report to combine with this one
    */
    public RecalibrationReport combine(final RecalibrationReport other) {
        if (other.isEmpty()){
            return this;
        }
        for ( int tableIndex = 0; tableIndex < recalibrationTables.numTables(); tableIndex++ ) {
            final NestedIntegerArray<RecalDatum> myTable = recalibrationTables.getTable(tableIndex);
            final NestedIntegerArray<RecalDatum> otherTable = other.recalibrationTables.getTable(tableIndex);
            RecalUtils.combineTables(myTable, otherTable);
        }
        return this;
    }

    public QuantizationInfo getQuantizationInfo() {
        return quantizationInfo;
    }

    public RecalibrationTables getRecalibrationTables() {
        return recalibrationTables;
    }

    public StandardCovariateList getCovariates() {
        return covariates;
    }

    /**
     * Initialize read group keys using the shared list of all the read groups.
     *
     * By using the same sorted set of read groups across all recalibration reports, even if
     * one report is missing a read group, all the reports use the same read group keys.
     *
     * @param allReadGroups The list of all possible read groups
     */
    private void initializeReadGroupCovariates(final SortedSet<String> allReadGroups) {
        for (final String readGroup: allReadGroups) {
            covariates.getReadGroupCovariate().keyFromValue(readGroup);
        }
    }

    /**
     * Compiles the list of keys for the Covariates table and uses the shared parsing utility to produce the actual table
     *
     * @param reportTable            the GATKReport table containing data for this table
     * @param recalibrationTables    the recalibration tables
     */
    private void parseAllCovariatesTable(final GATKReportTable reportTable, final RecalibrationTables recalibrationTables) {
        final int[] tempCOVarray = new int[4];

        for ( int i = 0; i < reportTable.getNumRows(); i++ ) {
            final Object rg = reportTable.get(i, RecalUtils.READGROUP_COLUMN_NAME);
            tempCOVarray[0] = covariates.getReadGroupCovariate().keyFromValue(rg);
            final Object qual = reportTable.get(i, RecalUtils.QUALITY_SCORE_COLUMN_NAME);
            tempCOVarray[1] = covariates.getQualityScoreCovariate().keyFromValue(qual);

            final String covName = (String)reportTable.get(i, RecalUtils.COVARIATE_NAME_COLUMN_NAME);
            final Object covValue = reportTable.get(i, RecalUtils.COVARIATE_VALUE_COLUMN_NAME);
            final Covariate cov = covariates.getCovariateByParsedName(covName);
            tempCOVarray[2] = cov.keyFromValue(covValue);

            final EventType event = EventType.eventFrom((String)reportTable.get(i, RecalUtils.EVENT_TYPE_COLUMN_NAME));
            tempCOVarray[3] = event.ordinal();

            recalibrationTables.getTableForCovariate(cov).put(getRecalDatum(reportTable, i, false), tempCOVarray);
        }
    }

    /**
     *
     * Compiles the list of keys for the QualityScore table and uses the shared parsing utility to produce the actual table
     * @param reportTable            the GATKReport table containing data for this table
     * @param qualTable               the map representing this table
     */
    private void parseQualityScoreTable(final GATKReportTable reportTable, final NestedIntegerArray<RecalDatum> qualTable) {
        final int[] tempQUALarray = new int[3];
        for ( int i = 0; i < reportTable.getNumRows(); i++ ) {
            final Object rg = reportTable.get(i, RecalUtils.READGROUP_COLUMN_NAME);
            tempQUALarray[0] = covariates.getReadGroupCovariate().keyFromValue(rg);
            final Object qual = reportTable.get(i, RecalUtils.QUALITY_SCORE_COLUMN_NAME);
            tempQUALarray[1] = covariates.getQualityScoreCovariate().keyFromValue(qual);
            final EventType event = EventType.eventFrom((String)reportTable.get(i, RecalUtils.EVENT_TYPE_COLUMN_NAME));
            tempQUALarray[2] = event.ordinal();

            qualTable.put(getRecalDatum(reportTable, i, false), tempQUALarray);
        }
    }

    /**
     * Compiles the list of keys for the ReadGroup table and uses the shared parsing utility to produce the actual table
     *
     * @param reportTable            the GATKReport table containing data for this table
     * @param rgTable                the map representing this table
     */
    private void parseReadGroupTable(final GATKReportTable reportTable, final NestedIntegerArray<RecalDatum> rgTable) {
        final int[] tempRGarray = new int[2];
        for ( int i = 0; i < reportTable.getNumRows(); i++ ) {
            final Object rg = reportTable.get(i, RecalUtils.READGROUP_COLUMN_NAME);
            tempRGarray[0] = covariates.getReadGroupCovariate().keyFromValue(rg);
            final EventType event = EventType.eventFrom((String)reportTable.get(i, RecalUtils.EVENT_TYPE_COLUMN_NAME));
            tempRGarray[1] = event.ordinal();

            rgTable.put(getRecalDatum(reportTable, i, true), tempRGarray);
        }
    }

    private static double asDouble(final Object o) {
        if ( o instanceof Double)
            return (Double)o;
        else if ( o instanceof Integer)
            return (Integer)o;
        else if ( o instanceof Long)
            return (Long)o;
        else
            throw new GATKException("Object " + o + " is expected to be either a double, long or integer but it's not either: " + o.getClass());
    }

    private static long asLong(final Object o) {
        if ( o instanceof Long)
            return (Long)o;
        else if ( o instanceof Integer)
            return ((Integer)o).longValue();
        else if ( o instanceof Double)
            return ((Double)o).longValue();
        else
            throw new GATKException("Object " + o + " is expected to be a long but it's not: " + o.getClass());
    }

    private static RecalDatum getRecalDatum(final GATKReportTable reportTable, final int row, final boolean hasEstimatedQReportedColumn) {
        final long nObservations = asLong(reportTable.get(row, RecalUtils.NUMBER_OBSERVATIONS_COLUMN_NAME));
        final double nErrors = asDouble(reportTable.get(row, RecalUtils.NUMBER_ERRORS_COLUMN_NAME));

        // the estimatedQreported column only exists in the ReadGroup table
        final double estimatedQReported = hasEstimatedQReportedColumn ?
                (Double) reportTable.get(row, RecalUtils.ESTIMATED_Q_REPORTED_COLUMN_NAME) : // we get it if we are in the read group table
                decodeByte(reportTable.get(row, RecalUtils.QUALITY_SCORE_COLUMN_NAME)); // or we use the reported quality if we are in any other table

        final RecalDatum datum = new RecalDatum(nObservations, nErrors, (byte)1);
        datum.setReportedQuality(estimatedQReported);
        //datum.setEmpiricalQuality(empiricalQuality); // don't set the value here because we will want to recompute with a different conditional Q score prior value
        return datum;
    }

    /**
     * Parses the quantization table from the GATK Report and turns it into a map of original => quantized quality scores
     *
     * @param table the GATKReportTable containing the quantization mappings
     * @return an ArrayList with the quantization mappings from 0 to MAX_SAM_QUAL_SCORE
     */
    private static QuantizationInfo initializeQuantizationTable(GATKReportTable table) {
        final Byte[] quals  = new Byte[QualityUtils.MAX_SAM_QUAL_SCORE + 1];
        final Long[] counts = new Long[QualityUtils.MAX_SAM_QUAL_SCORE + 1];
        for ( int i = 0; i < table.getNumRows(); i++ ) {
            final byte originalQual = (byte)i;
            final Object quantizedObject = table.get(i, RecalUtils.QUANTIZED_VALUE_COLUMN_NAME);
            final Object countObject = table.get(i, RecalUtils.QUANTIZED_COUNT_COLUMN_NAME);
            final byte quantizedQual = Byte.parseByte(quantizedObject.toString());
            final long quantizedCount = Long.parseLong(countObject.toString());
            quals[originalQual] = quantizedQual;
            counts[originalQual] = quantizedCount;
        }
        return new QuantizationInfo(new ArrayList<>(Arrays.asList(quals)), new ArrayList<>(Arrays.asList(counts)));
    }

    /**
     * Parses the arguments table from the GATK Report and creates a RAC object with the proper initialization values
     *
     * @param table the GATKReportTable containing the arguments and its corresponding values
     * @return a RAC object properly initialized with all the objects in the table
     */
    private static RecalibrationArgumentCollection initializeArgumentCollectionTable(GATKReportTable table) {
        final RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();

        final List<String> standardCovariateClassNames = new StandardCovariateList(RAC, Collections.emptyList()).getStandardCovariateClassNames();

        for ( int i = 0; i < table.getNumRows(); i++ ) {
            final String argument = table.get(i, "Argument").toString();
            Object value = table.get(i, RecalUtils.ARGUMENT_VALUE_COLUMN_NAME);
            if (value.equals("null")) {
                value = null; // generic translation of null values that were printed out as strings | todo -- add this capability to the GATKReport
            }

            if (argument.equals("covariate") && value != null) {
                final List<String> covs = new ArrayList<>(Arrays.asList(value.toString().split(",")));
                if (!covs.equals(standardCovariateClassNames)) {
                    throw new UserException("Non-standard covariates are not supported. Only the following are supported " + standardCovariateClassNames + " but was " + covs);
                }
            } else if (argument.equals("no_standard_covs")) {
                final boolean no_standard_covs = decodeBoolean(value);
                if (no_standard_covs){
                    throw new UserException("Non-standard covariates are not supported. Only the following are supported " + standardCovariateClassNames + " but no_standard_covs was true");
                }
            } else if (argument.equals("solid_recal_mode")) {
                final String solid_recal_mode = (String) value;
                if (!RecalibrationArgumentCollection.SOLID_RECAL_MODE.equals(solid_recal_mode)){
                    throw new UserException("Solid is not supported. Only " + RecalibrationArgumentCollection.SOLID_RECAL_MODE + " is allowed as value for solid_recal_mode");
                }
            }
            else if (argument.equals("solid_nocall_strategy")) {
                final String solid_nocall_strategy = (String) value;
                if (!RecalibrationArgumentCollection.SOLID_NOCALL_STRATEGY.equals(solid_nocall_strategy)){
                    throw new UserException("Solid is not supported. Only " + RecalibrationArgumentCollection.SOLID_NOCALL_STRATEGY + " is allowed as value for solid_nocall_strategy");
                }
            }
            else if (argument.equals("mismatches_context_size"))
                RAC.MISMATCHES_CONTEXT_SIZE = decodeInteger(value);

            else if (argument.equals("indels_context_size"))
                RAC.INDELS_CONTEXT_SIZE = decodeInteger(value);

            else if (argument.equals("mismatches_default_quality"))
                RAC.MISMATCHES_DEFAULT_QUALITY = decodeByte(value);

            else if (argument.equals("insertions_default_quality"))
                RAC.INSERTIONS_DEFAULT_QUALITY = decodeByte(value);

            else if (argument.equals("deletions_default_quality"))
                RAC.DELETIONS_DEFAULT_QUALITY = decodeByte(value);

            else if (argument.equals("maximum_cycle_value"))
                RAC.MAXIMUM_CYCLE_VALUE = decodeInteger(value);

            else if (argument.equals("low_quality_tail"))
                RAC.LOW_QUAL_TAIL = decodeByte(value);

            else if (argument.equals("default_platform"))
                RAC.DEFAULT_PLATFORM = (String) value;

            else if (argument.equals("force_platform"))
                RAC.FORCE_PLATFORM = (String) value;

            else if (argument.equals("quantizing_levels"))
                RAC.QUANTIZING_LEVELS = decodeInteger(value);

            else if (argument.equals("recalibration_report"))
                RAC.existingRecalibrationReport = (value == null) ? null : new File((String) value);

            else if (argument.equals("binary_tag_name"))
                RAC.BINARY_TAG_NAME = (value == null) ? null : (String) value;
        }

        return RAC;
    }

    private static byte decodeByte(final Object value) {
        if (value instanceof Byte){
            return (Byte)value;
        } else if ( value instanceof String ) {
            return Byte.parseByte((String)value);
        } else if ( value instanceof Long ) {
            return ((Long) value).byteValue();
        } else {
            throw new IllegalArgumentException("expected a Byte, String, or Long, but got " + value);
        }
    }

    private static int decodeInteger(final Object value) {
        Utils.validateArg(value instanceof Integer || value instanceof String, () -> "expected an Integer or a String but got " + value);
        return value instanceof Integer ? (Integer) value: Integer.parseInt((String)value);
    }

    private static boolean decodeBoolean(final Object value) {
        Utils.validateArg(value instanceof Boolean || value instanceof String, () -> "expected a Boolean or a String but got " + value);
        return value instanceof Boolean ? (Boolean) value: Boolean.parseBoolean((String) value);
    }

    /**
     * Creates the recalibration report.  Report can then be written to a stream via GATKReport.print(PrintStream).
     *
     * @return newly created recalibration report
     */
    public GATKReport createGATKReport() {
        return RecalUtils.createRecalibrationGATKReport(argumentTable, quantizationInfo, recalibrationTables, covariates);
    }

    public RecalibrationArgumentCollection getRAC() {
        return RAC;
    }

    /**
     * @return true if the report has no data
     */
    public boolean isEmpty() {
        return recalibrationTables.isEmpty();
    }

}
