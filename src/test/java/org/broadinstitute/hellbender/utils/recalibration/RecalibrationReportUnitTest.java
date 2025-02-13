package org.broadinstitute.hellbender.utils.recalibration;

import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.collections.NestedIntegerArray;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.covariates.*;
import org.broadinstitute.hellbender.utils.report.GATKReport;
import org.broadinstitute.hellbender.utils.report.GATKReportTable;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.PrintStream;
import java.util.*;

public final class RecalibrationReportUnitTest extends GATKBaseTest {

    private static final String testDir = GATKBaseTest.publicTestDir + "/org/broadinstitute/hellbender/utils/recalibration/";
    private static final File recal1 = new File(testDir + "HiSeq.1mb.1RG.sg1.table");
    private static final File recal2 = new File(testDir + "HiSeq.1mb.1RG.sg2.table");
    private static final File recal3 = new File(testDir + "HiSeq.1mb.1RG.sg3.table");
    private static final File recal4 = new File(testDir + "HiSeq.1mb.1RG.sg4.table");
    private static final File recal5 = new File(testDir + "HiSeq.1mb.1RG.sg5.table");

    private static final ImmutableList<File> recalFiles = ImmutableList.of(recal1, recal2, recal3, recal4, recal5);

    private static final File recalEmpty = new File(testDir + "HiSeq.1mb.1RG.empty.table");

    //Note: the expected output was created by running GATK3, fixing table headers and sorting in unix "sort -k 2,2 -k 3,3 -k 4,4"
    private static final File recal_original = new File(testDir + "HiSeq.1mb.1RG.noSG.table");

    @DataProvider(name = "tables")
    public Object[][] getTables(){
        final List<File> filesWithEmpty = new ArrayList<>(recalFiles);
        filesWithEmpty.add(recalEmpty);

        return new Object[][] {
                {Collections.nCopies(5, new File(testDir + "bqsr.manyObservations.piece.table")), new File(testDir + "bqsr.manyObservations.full.table")}, //gather files with multiple observations
                {recalFiles, recal_original}, //gather several files
                {filesWithEmpty, recal_original} //gather files with an empty table
        };
    }

    @Test(dataProvider = "tables")
    public void testGatherBQSR(List<File> inputTables, File expectedOutputTable) {
        testGatherReports(inputTables, expectedOutputTable);
    }

    public static void testGatherReports(List<File> inputFiles, File expectedResult) {
        final File output = GATKBaseTest.createTempFile("BQSRgathererTest", ".table");
        RecalibrationReport.gatherReportsIntoOneFile(inputFiles, output);

        final GATKReport originalReport = new GATKReport(expectedResult);
        final GATKReport calculatedReport = new GATKReport(output);

        assertReportsAreEquivalent(originalReport, calculatedReport);
    }

    private static void assertReportsAreEquivalent(final GATKReport originalReport, final GATKReport calculatedReport) {

        final ByteArrayOutputStream original = new ByteArrayOutputStream();
        originalReport.print(new PrintStream(original));

        final ByteArrayOutputStream calculated = new ByteArrayOutputStream();
        calculatedReport.print(new PrintStream(calculated));

        Assert.assertEquals(calculated.toString(),original.toString());

        // test the Arguments table
        final List<String> argumentTableColumns = Arrays.asList(RecalUtils.ARGUMENT_COLUMN_NAME, RecalUtils.ARGUMENT_VALUE_COLUMN_NAME);
        testReportsForTable(originalReport, calculatedReport, argumentTableColumns, RecalUtils.ARGUMENT_REPORT_TABLE_TITLE);

        // test the Quantized table
        final List<String> quantizedTableColumns = Arrays.asList(RecalUtils.QUALITY_SCORE_COLUMN_NAME, RecalUtils.QUANTIZED_COUNT_COLUMN_NAME, RecalUtils.QUANTIZED_VALUE_COLUMN_NAME);
        testReportsForTable(originalReport, calculatedReport, quantizedTableColumns, RecalUtils.QUANTIZED_REPORT_TABLE_TITLE);

        // test the RecalTable0 table
        final List<String> readGroupTableColumns = Arrays.asList(RecalUtils.READGROUP_COLUMN_NAME, RecalUtils.EVENT_TYPE_COLUMN_NAME, RecalUtils.ESTIMATED_Q_REPORTED_COLUMN_NAME, RecalUtils.NUMBER_OBSERVATIONS_COLUMN_NAME, RecalUtils.NUMBER_ERRORS_COLUMN_NAME);
        testReportsForTable(originalReport, calculatedReport, readGroupTableColumns, RecalUtils.READGROUP_REPORT_TABLE_TITLE);

        // test the RecalTable1 table
        final List<String> qualityScoreTableColumns = Arrays.asList(RecalUtils.READGROUP_COLUMN_NAME, RecalUtils.QUALITY_SCORE_COLUMN_NAME, RecalUtils.EVENT_TYPE_COLUMN_NAME, RecalUtils.NUMBER_OBSERVATIONS_COLUMN_NAME, RecalUtils.NUMBER_ERRORS_COLUMN_NAME);
        testReportsForTable(originalReport, calculatedReport, qualityScoreTableColumns, RecalUtils.QUALITY_SCORE_REPORT_TABLE_TITLE);

        // test the RecalTable2 table
        final List<String> allCovariatesTableColumns = Arrays.asList(RecalUtils.READGROUP_COLUMN_NAME, RecalUtils.QUALITY_SCORE_COLUMN_NAME, RecalUtils.COVARIATE_VALUE_COLUMN_NAME, RecalUtils.COVARIATE_NAME_COLUMN_NAME, RecalUtils.EVENT_TYPE_COLUMN_NAME, RecalUtils.NUMBER_OBSERVATIONS_COLUMN_NAME, RecalUtils.NUMBER_ERRORS_COLUMN_NAME, RecalUtils.EMPIRICAL_QUALITY_COLUMN_NAME);
        testReportsForTable(originalReport, calculatedReport, allCovariatesTableColumns, RecalUtils.ALL_COVARIATES_REPORT_TABLE_TITLE);

    }

    private static void testReportsForTable(GATKReport originalReport, GATKReport calculatedReport, List<String> columnsToTest, String argumentReportTableTitle) {
        final GATKReportTable originalTable = originalReport.getTable(argumentReportTableTitle);
        final GATKReportTable calculatedTable = calculatedReport.getTable(argumentReportTableTitle);
        testTablesWithColumns(originalTable, calculatedTable, columnsToTest);
    }

    /**
     * Common testing functionality given the columns to test
     *
     * @param original the original table
     * @param calculated the calculated table
     * @param columnsToTest list of columns to test. All columns will be tested with the same criteria
     */
    private static void testTablesWithColumns(GATKReportTable original, GATKReportTable calculated, List<String> columnsToTest) {
        for (int row = 0; row < original.getNumRows(); row++ ) {
            for (String column : columnsToTest) {
                final Object actual = calculated.get(row, column);
                final Object expected = original.get(row, column);
                //if ( !actual.equals(expected) )
                //    System.out.println("Row=" + row + " Table=" + original.getTableName() + " Column=" + column + " Expected=" + expected + " Actual=" + actual);
                Assert.assertEquals(actual, expected, "Row: " + row + " Original Table: " + original.getTableName() + " Column=" + column);
            }
        }

    }

    @Test
    public void testGatherMissingReadGroup() {
        // Hand modified subset of problematic gather inputs submitted by George Grant
        final File input1 = new File(testDir + "NA12878.rg_subset.chr1.recal_data.table");
        final File input2 = new File(testDir + "NA12878.rg_subset.chrY_Plus.recal_data.table");

        final GATKReport report12 = RecalibrationReport.gatherReports(Arrays.asList(input1, input2));
        final GATKReport report21 = RecalibrationReport.gatherReports(Arrays.asList(input2, input1));

        Assert.assertTrue(report12.equals(report21), "GATK reports are different when gathered in a different order.");
    }

    private static RecalDatum createRandomRecalDatum(int maxObservations, int maxErrors) {
        final Random random = new Random();
        final int nObservations = random.nextInt(maxObservations);
        final int nErrors = Math.min(random.nextInt(maxErrors), nObservations);
        final int qual = random.nextInt(QualityUtils.MAX_SAM_QUAL_SCORE);
        return new RecalDatum((long)nObservations, (double)nErrors, (byte)qual);
    }

    @Test(expectedExceptions = UserException.class)
    public void testUnsupportedCovariates(){
        File file = new File(toolsTestDir + "unsupported-covariates.table.gz");
        new RecalibrationReport(file);
    }

    @Test
    public void testOutput() {
        final int length = 100;

        List<Byte> quals = new ArrayList<>(QualityUtils.MAX_SAM_QUAL_SCORE + 1);
        List<Long> counts = new ArrayList<>(QualityUtils.MAX_SAM_QUAL_SCORE + 1);

        for (int i = 0;  i<= QualityUtils.MAX_SAM_QUAL_SCORE; i++) {
            quals.add((byte) i);
            counts.add(1L);
        }

        final QuantizationInfo quantizationInfo = new QuantizationInfo(quals, counts);
        final RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();

        quantizationInfo.noQuantization();
        final String readGroupID = "id";
        final StandardCovariateList covariateList = new StandardCovariateList(RAC, Collections.singletonList(readGroupID));

        final SAMReadGroupRecord rg = new SAMReadGroupRecord(readGroupID);
        rg.setPlatform("illumina");
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeaderWithReadGroup(rg);

        final GATKRead read = ArtificialReadUtils.createRandomRead(header, length, false);
        read.setReadGroup(rg.getReadGroupId());

        final byte [] readQuals = new byte[length];
        for (int i = 0; i < length; i++)
            readQuals[i] = 20;
        read.setBaseQualities(readQuals);

        final int expectedKeys = expectedNumberOfKeys(length, RAC.INDELS_CONTEXT_SIZE, RAC.MISMATCHES_CONTEXT_SIZE);
        int nKeys = 0;  // keep track of how many keys were produced
        final PerReadCovariateMatrix rc = RecalUtils.computeCovariates(read, header, covariateList, true, new CovariateKeyCache());

        final RecalibrationTables recalibrationTables = new RecalibrationTables(covariateList);
        final NestedIntegerArray<RecalDatum> rgTable = recalibrationTables.getReadGroupTable();
        final NestedIntegerArray<RecalDatum> qualTable = recalibrationTables.getQualityScoreTable();

        for (int offset = 0; offset < length; offset++) {

            for (EventType errorMode : EventType.values()) {

                final int[] covariates = rc.getCovariatesAtOffset(offset, errorMode);
                final int randomMax = errorMode == EventType.BASE_SUBSTITUTION ? 10000 : 100000;

                rgTable.put(createRandomRecalDatum(randomMax, 10), covariates[0], errorMode.ordinal());
                qualTable.put(createRandomRecalDatum(randomMax, 10), covariates[0], covariates[1], errorMode.ordinal());
                nKeys += 2;
                for (NestedIntegerArray<RecalDatum> covTable : recalibrationTables.getAdditionalTables()){
                    Covariate cov = recalibrationTables.getCovariateForTable(covTable);
                    final int covValue = covariates[covariateList.indexByClass(cov.getClass())];
                    if ( covValue >= 0 ) {
                        covTable.put(createRandomRecalDatum(randomMax, 10), covariates[0], covariates[1], covValue, errorMode.ordinal());
                        nKeys++;
                    }
                }
            }
        }
        Assert.assertEquals(nKeys, expectedKeys);
    }

    private static int expectedNumberOfKeys (int readLength, int indelContextSize, int mismatchesContextSize) {
        final int numCovariates = 4;
        final int numTables = 3;
        final int mismatchContextPadding = mismatchesContextSize - 1;
        final int indelContextPadding = 2 * (indelContextSize - 1);
        final int indelCyclePadding = 2 * (2 * CycleCovariate.CUSHION_FOR_INDELS);

        return (numCovariates * numTables * readLength) - mismatchContextPadding - indelContextPadding - indelCyclePadding;
    }

}
