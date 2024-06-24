package org.broadinstitute.hellbender.tools.funcotator.dataSources;

import org.apache.commons.lang3.RandomStringUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.time.LocalDate;
import java.time.Month;
import java.time.YearMonth;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Matcher;

/**
 * Test class for {@link DataSourceUtils}
 * Created by jonn on 3/22/18.
 */
public class DataSourceUtilsUnitTest extends GATKBaseTest {

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Private Members:

    //==================================================================================================================
    // Helper Methods:

    /**
     * Add a test case to the given {@code testArgs} for version-based tests.
     *
     * Calculates the expected value based on the validity of the input arguments and the valid data ranges in
     * {@link DataSourceUtils}.
     *
     * @param testArgs {@code List<Object[]>} test arguments list to which to add the test case.
     * @param maj {@code int} representing the major version to test.
     * @param min {@code int} representing the minor version to test.
     * @param date {@link LocalDate} represeting the date of the version to test.
     */
    private void addVersionTestCase(final List<Object[]> testArgs,
                                    final int maj,
                                    final int min,
                                    final LocalDate date) {
        testArgs.add(
                new Object[]{
                        maj, min, date,
                        (
                            (maj >= DataSourceUtils.MIN_MAJOR_VERSION_NUMBER) &&
                            (maj <= DataSourceUtils.MAX_MAJOR_VERSION_NUMBER) &&
                            (maj != DataSourceUtils.MIN_MAJOR_VERSION_NUMBER || min >= DataSourceUtils.MIN_MINOR_VERSION_NUMBER) &&
                            (maj != DataSourceUtils.MAX_MAJOR_VERSION_NUMBER || min <= DataSourceUtils.MAX_MINOR_VERSION_NUMBER)
                        )
                        &&
                        (
                            (date.isAfter(DataSourceUtils.MIN_DATE) || date.isEqual(DataSourceUtils.MIN_DATE)) &&
                            (date.isBefore(DataSourceUtils.MAX_DATE) || date.isEqual(DataSourceUtils.MAX_DATE))
                        )
                }
        );
    }

    /**
     * Add a test case to the given {@code testArgs} for version-based tests.
     * @param testArgs {@code List<Object[]>} test arguments list to which to add the test case.
     * @param maj {@code int} representing the major version to test.
     * @param min {@code int} representing the minor version to test.
     * @param date {@link LocalDate} represeting the date of the version to test.
     * @param expected The expected {@code boolean} value for the test case given the input version data.
     */
    private void addVersionTestCase(final List<Object[]> testArgs,
                                    final int maj,
                                    final int min,
                                    final LocalDate date,
                                    final boolean expected) {
        testArgs.add(new Object[]{maj, min, date, expected});
    }

    private List<Object[]>  createBaseTestVersionData() {

        final ArrayList<Object[]> testArgs = new ArrayList<>();

        // =============================================================
        // First do tests that check the major / minor version numbers:
        // -------------------------------------------------------------

        // major < MIN_MAJOR (minor shouldn't matter)
        addVersionTestCase(testArgs, DataSourceUtils.MIN_MAJOR_VERSION_NUMBER - 1, 0, DataSourceUtils.MIN_DATE, false);
        addVersionTestCase(testArgs, DataSourceUtils.MIN_MAJOR_VERSION_NUMBER - 1, DataSourceUtils.MIN_MINOR_VERSION_NUMBER - 1, DataSourceUtils.MIN_DATE, false);
        addVersionTestCase(testArgs, DataSourceUtils.MIN_MAJOR_VERSION_NUMBER - 1, DataSourceUtils.MIN_MINOR_VERSION_NUMBER, DataSourceUtils.MIN_DATE, false);
        addVersionTestCase(testArgs, DataSourceUtils.MIN_MAJOR_VERSION_NUMBER - 1, DataSourceUtils.MIN_MINOR_VERSION_NUMBER + 1, DataSourceUtils.MIN_DATE, false);
        addVersionTestCase(testArgs, DataSourceUtils.MIN_MAJOR_VERSION_NUMBER - 1, DataSourceUtils.MAX_MINOR_VERSION_NUMBER - 1, DataSourceUtils.MIN_DATE, false);
        addVersionTestCase(testArgs, DataSourceUtils.MIN_MAJOR_VERSION_NUMBER - 1, DataSourceUtils.MAX_MINOR_VERSION_NUMBER, DataSourceUtils.MIN_DATE, false);
        addVersionTestCase(testArgs, DataSourceUtils.MIN_MAJOR_VERSION_NUMBER - 1, DataSourceUtils.MAX_MINOR_VERSION_NUMBER + 1, DataSourceUtils.MIN_DATE, false);

        // major == MIN_MAJOR (minor does matter)
        // Truth computed in `addVersionTestCase` method:
        addVersionTestCase(testArgs, DataSourceUtils.MIN_MAJOR_VERSION_NUMBER, 0, DataSourceUtils.MIN_DATE);
        addVersionTestCase(testArgs, DataSourceUtils.MIN_MAJOR_VERSION_NUMBER, DataSourceUtils.MIN_MINOR_VERSION_NUMBER - 1, DataSourceUtils.MIN_DATE);
        addVersionTestCase(testArgs, DataSourceUtils.MIN_MAJOR_VERSION_NUMBER, DataSourceUtils.MIN_MINOR_VERSION_NUMBER, DataSourceUtils.MIN_DATE);
        addVersionTestCase(testArgs, DataSourceUtils.MIN_MAJOR_VERSION_NUMBER, DataSourceUtils.MIN_MINOR_VERSION_NUMBER + 1, DataSourceUtils.MIN_DATE);
        addVersionTestCase(testArgs, DataSourceUtils.MIN_MAJOR_VERSION_NUMBER, DataSourceUtils.MAX_MINOR_VERSION_NUMBER - 1, DataSourceUtils.MIN_DATE);
        addVersionTestCase(testArgs, DataSourceUtils.MIN_MAJOR_VERSION_NUMBER, DataSourceUtils.MAX_MINOR_VERSION_NUMBER, DataSourceUtils.MIN_DATE);
        addVersionTestCase(testArgs, DataSourceUtils.MIN_MAJOR_VERSION_NUMBER, DataSourceUtils.MAX_MINOR_VERSION_NUMBER + 1, DataSourceUtils.MIN_DATE);

        // major > MIN_MAJOR and major < MAX_MAJOR (minor shouldn't matter)
        if ( DataSourceUtils.MAX_MAJOR_VERSION_NUMBER - DataSourceUtils.MIN_MAJOR_VERSION_NUMBER > 1 ) {
            addVersionTestCase(testArgs, DataSourceUtils.MIN_MAJOR_VERSION_NUMBER + 1, 0, DataSourceUtils.MIN_DATE, true);
            addVersionTestCase(testArgs, DataSourceUtils.MIN_MAJOR_VERSION_NUMBER + 1, DataSourceUtils.MIN_MINOR_VERSION_NUMBER - 1, DataSourceUtils.MIN_DATE, true);
            addVersionTestCase(testArgs, DataSourceUtils.MIN_MAJOR_VERSION_NUMBER + 1, DataSourceUtils.MIN_MINOR_VERSION_NUMBER, DataSourceUtils.MIN_DATE, true);
            addVersionTestCase(testArgs, DataSourceUtils.MIN_MAJOR_VERSION_NUMBER + 1, DataSourceUtils.MIN_MINOR_VERSION_NUMBER + 1, DataSourceUtils.MIN_DATE, true);
            addVersionTestCase(testArgs, DataSourceUtils.MIN_MAJOR_VERSION_NUMBER + 1, DataSourceUtils.MAX_MINOR_VERSION_NUMBER - 1, DataSourceUtils.MIN_DATE, true);
            addVersionTestCase(testArgs, DataSourceUtils.MIN_MAJOR_VERSION_NUMBER + 1, DataSourceUtils.MAX_MINOR_VERSION_NUMBER, DataSourceUtils.MIN_DATE, true);
            addVersionTestCase(testArgs, DataSourceUtils.MIN_MAJOR_VERSION_NUMBER + 1, DataSourceUtils.MAX_MINOR_VERSION_NUMBER + 1, DataSourceUtils.MIN_DATE, true);
        }

        // major == MAX_MAJOR (minor does matter)
        // Truth computed in `addVersionTestCase` method:
        addVersionTestCase(testArgs, DataSourceUtils.MAX_MAJOR_VERSION_NUMBER, 0, DataSourceUtils.MIN_DATE);
        addVersionTestCase(testArgs, DataSourceUtils.MAX_MAJOR_VERSION_NUMBER, DataSourceUtils.MIN_MINOR_VERSION_NUMBER - 1, DataSourceUtils.MIN_DATE);
        addVersionTestCase(testArgs, DataSourceUtils.MAX_MAJOR_VERSION_NUMBER, DataSourceUtils.MIN_MINOR_VERSION_NUMBER, DataSourceUtils.MIN_DATE);
        addVersionTestCase(testArgs, DataSourceUtils.MAX_MAJOR_VERSION_NUMBER, DataSourceUtils.MIN_MINOR_VERSION_NUMBER + 1, DataSourceUtils.MIN_DATE);
        addVersionTestCase(testArgs, DataSourceUtils.MAX_MAJOR_VERSION_NUMBER, DataSourceUtils.MAX_MINOR_VERSION_NUMBER - 1, DataSourceUtils.MIN_DATE);
        addVersionTestCase(testArgs, DataSourceUtils.MAX_MAJOR_VERSION_NUMBER, DataSourceUtils.MAX_MINOR_VERSION_NUMBER, DataSourceUtils.MIN_DATE);
        addVersionTestCase(testArgs, DataSourceUtils.MAX_MAJOR_VERSION_NUMBER, DataSourceUtils.MAX_MINOR_VERSION_NUMBER + 1, DataSourceUtils.MIN_DATE);

        // major > MAX_MAJOR (minor shouldn't matter)
        addVersionTestCase(testArgs, DataSourceUtils.MAX_MAJOR_VERSION_NUMBER + 1, 0, DataSourceUtils.MAX_DATE, false);
        addVersionTestCase(testArgs, DataSourceUtils.MAX_MAJOR_VERSION_NUMBER + 1, DataSourceUtils.MIN_MINOR_VERSION_NUMBER - 1, DataSourceUtils.MAX_DATE, false);
        addVersionTestCase(testArgs, DataSourceUtils.MAX_MAJOR_VERSION_NUMBER + 1, DataSourceUtils.MIN_MINOR_VERSION_NUMBER, DataSourceUtils.MAX_DATE, false);
        addVersionTestCase(testArgs, DataSourceUtils.MAX_MAJOR_VERSION_NUMBER + 1, DataSourceUtils.MIN_MINOR_VERSION_NUMBER + 1, DataSourceUtils.MAX_DATE, false);
        addVersionTestCase(testArgs, DataSourceUtils.MAX_MAJOR_VERSION_NUMBER + 1, DataSourceUtils.MAX_MINOR_VERSION_NUMBER - 1, DataSourceUtils.MAX_DATE, false);
        addVersionTestCase(testArgs, DataSourceUtils.MAX_MAJOR_VERSION_NUMBER + 1, DataSourceUtils.MAX_MINOR_VERSION_NUMBER, DataSourceUtils.MAX_DATE, false);
        addVersionTestCase(testArgs, DataSourceUtils.MAX_MAJOR_VERSION_NUMBER + 1, DataSourceUtils.MAX_MINOR_VERSION_NUMBER + 1, DataSourceUtils.MAX_DATE, false);

        // =============================================================
        // Next do tests that check the dates themselves:
        // -------------------------------------------------------------

        // ----------------------------------------
        // Dates in range:
        // ----------------------------------------

        // Min version with minimally / maximally acceptable dates:
        testArgs.add(
                new Object[]{
                        DataSourceUtils.MIN_MAJOR_VERSION_NUMBER,
                        DataSourceUtils.MIN_MINOR_VERSION_NUMBER,
                        DataSourceUtils.MIN_DATE,
                        true
                }
        );
        testArgs.add(
                new Object[]{
                        DataSourceUtils.MIN_MAJOR_VERSION_NUMBER,
                        DataSourceUtils.MIN_MINOR_VERSION_NUMBER,
                        DataSourceUtils.MAX_DATE,
                        true
                }
        );

        // Max version with minimally / maximally acceptable dates:
        testArgs.add(
                new Object[]{
                        DataSourceUtils.MAX_MAJOR_VERSION_NUMBER,
                        DataSourceUtils.MAX_MINOR_VERSION_NUMBER,
                        DataSourceUtils.MIN_DATE,
                        true
                }
        );
        testArgs.add(
                new Object[]{
                        DataSourceUtils.MAX_MAJOR_VERSION_NUMBER,
                        DataSourceUtils.MAX_MINOR_VERSION_NUMBER,
                        DataSourceUtils.MAX_DATE,
                        true
                }
        );

        // Year / Month inside acceptable dates with any valid day:
        final int validYear = DataSourceUtils.MIN_DATE.getMonthValue() == 12 ? DataSourceUtils.MIN_DATE.getYear() + 1 : DataSourceUtils.MIN_DATE.getYear();
        final int validMonth = DataSourceUtils.MIN_DATE.getMonthValue() == 12 ? 1 : DataSourceUtils.MIN_DATE.getMonthValue() + 1;
        for (int day = 1; day < YearMonth.of(validYear, validMonth).lengthOfMonth() + 1; ++day) {
            testArgs.add(
                    new Object[]{
                            DataSourceUtils.MIN_MAJOR_VERSION_NUMBER,
                            DataSourceUtils.MIN_MINOR_VERSION_NUMBER,
                            LocalDate.of(validYear, validMonth, day),
                            true
                    }
            );
        }

        // 1 day inside the acceptable release window:
        testArgs.add(
                new Object[]{
                        DataSourceUtils.MIN_MAJOR_VERSION_NUMBER,
                        DataSourceUtils.MIN_MINOR_VERSION_NUMBER,
                        LocalDate.of(
                                DataSourceUtils.MIN_DATE.getYear(),
                                DataSourceUtils.MIN_DATE.getMonth(),
                                DataSourceUtils.MIN_DATE.getDayOfMonth() + 1
                        ),
                        true
                }
        );
        testArgs.add(
                new Object[]{
                        DataSourceUtils.MAX_MAJOR_VERSION_NUMBER,
                        DataSourceUtils.MAX_MINOR_VERSION_NUMBER,
                        LocalDate.of(
                                DataSourceUtils.MAX_DATE.getYear(),
                                DataSourceUtils.MAX_DATE.getMonth(),
                                DataSourceUtils.MAX_DATE.getDayOfMonth() - 1
                        ),
                        true
                }
        );

        // ----------------------------------------
        // Dates out of range:
        // ----------------------------------------

        // 1 day outside the acceptable release window:
        testArgs.add(
                new Object[]{
                        DataSourceUtils.MIN_MAJOR_VERSION_NUMBER,
                        DataSourceUtils.MIN_MINOR_VERSION_NUMBER,
                        LocalDate.of(
                                DataSourceUtils.MIN_DATE.getYear(),
                                DataSourceUtils.MIN_DATE.getMonth(),
                                DataSourceUtils.MIN_DATE.getDayOfMonth() - 1
                        ),
                        false
                }
        );
        testArgs.add(
                new Object[]{
                        DataSourceUtils.MAX_MAJOR_VERSION_NUMBER,
                        DataSourceUtils.MAX_MINOR_VERSION_NUMBER,
                        LocalDate.of(
                                DataSourceUtils.MAX_DATE.getYear(),
                                DataSourceUtils.MAX_DATE.getMonth(),
                                DataSourceUtils.MAX_DATE.getDayOfMonth() + 1
                        ),
                        false
                }
        );

        // 1 month outside the acceptable release window:
        final int yearForPrevMonth = DataSourceUtils.MIN_DATE.getMonthValue() == 1 ? DataSourceUtils.MIN_DATE.getYear() - 1 : DataSourceUtils.MIN_DATE.getYear();
        final int invalidMonthTooEarly = DataSourceUtils.MIN_DATE.getMonthValue() == 12 ? 1 : DataSourceUtils.MIN_DATE.getMonthValue() + 1;

        final int yearForNextMonth = DataSourceUtils.MAX_DATE.getMonthValue() == 12 ? DataSourceUtils.MAX_DATE.getYear() + 1 : DataSourceUtils.MAX_DATE.getYear();
        final int invalidMonthTooLate = DataSourceUtils.MAX_DATE.getMonthValue() == 12 ? 1 : DataSourceUtils.MAX_DATE.getMonthValue() + 1;

        for (int day = 1; day < YearMonth.of(yearForPrevMonth, invalidMonthTooEarly).lengthOfMonth() + 1; ++day) {
            testArgs.add(
                    new Object[]{
                            DataSourceUtils.MIN_MAJOR_VERSION_NUMBER,
                            DataSourceUtils.MIN_MINOR_VERSION_NUMBER,
                            LocalDate.of(yearForNextMonth, invalidMonthTooLate, day),
                            false
                    }
            );
        }

        for (int day = 1; day < YearMonth.of(yearForNextMonth, invalidMonthTooLate).lengthOfMonth() + 1; ++day) {
            testArgs.add(
                    new Object[]{
                            DataSourceUtils.MAX_MAJOR_VERSION_NUMBER,
                            DataSourceUtils.MAX_MINOR_VERSION_NUMBER,
                            LocalDate.of(yearForNextMonth, invalidMonthTooLate, day),
                            false
                    }
            );
        }

        // 1 year outside the acceptable release window:
        testArgs.add(
                new Object[]{
                        DataSourceUtils.MIN_MAJOR_VERSION_NUMBER,
                        DataSourceUtils.MIN_MINOR_VERSION_NUMBER,
                        LocalDate.of(
                                DataSourceUtils.MIN_DATE.getYear() - 1,
                                DataSourceUtils.MIN_DATE.getMonth(),
                                DataSourceUtils.MIN_DATE.getDayOfMonth()
                        ),
                        false
                }
        );
        testArgs.add(
                new Object[]{
                        DataSourceUtils.MAX_MAJOR_VERSION_NUMBER,
                        DataSourceUtils.MAX_MINOR_VERSION_NUMBER,
                        LocalDate.of(
                                DataSourceUtils.MAX_DATE.getYear() + 1,
                                DataSourceUtils.MAX_DATE.getMonth(),
                                DataSourceUtils.MAX_DATE.getDayOfMonth()
                        ),
                        false
                }
        );

        return testArgs;
    }

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    private Iterator<Object[]> provideForValidateVersionInformation() {
        return createBaseTestVersionData().iterator();
    }

    @DataProvider
    private Iterator<Object[]> provideForTestOldVersionRegex() {

        final ArrayList<Object[]> testArgs = new ArrayList<>();

        final List<Object[]> baseArgs = createBaseTestVersionData();

        for ( final Object[] args : baseArgs ) {
            for ( int whitespace = 0; whitespace < 2; ++whitespace ) {
                for ( int decoratorCount = 0; decoratorCount < 10; ++decoratorCount ) {

                    // Some sanity checks here for proper version numbers:
                    if (((Integer)args[0]) < 0 || ((Integer)args[1]) < 0){
                        continue;
                    }

                    final String whitespaceString = whitespace != 0 ? "\t \t \t " : " ";
                    final String decoratorString  = decoratorCount  !=0 ? RandomStringUtils.randomAlphanumeric(decoratorCount) : "";

                    testArgs.add(
                            new Object[] {
                                    args[0],args[1],args[2],
                                    decoratorString,
                                    whitespaceString
                            }
                    );
                }
            }
        }

        return testArgs.iterator();

    }

    @DataProvider
    private Iterator<Object[]> provideForTestNewVersionRegex() {

        final ArrayList<Object[]> testArgs = new ArrayList<>();

        final List<Object[]> baseArgs = createBaseTestVersionData();
        final String[] refVersions = new String[] { "hg19", "hg38"};

        for ( final Object[] args : baseArgs ) {
            for ( int whitespace = 0; whitespace < 2; ++whitespace ) {
                for ( int decoratorCount = 0; decoratorCount < 10; ++decoratorCount ) {
                    for ( int refidx = 0; refidx < refVersions.length; ++refidx ) {

                        // Some sanity checks here for proper version numbers:
                        if (((Integer) args[0]) < 0 || ((Integer) args[1]) < 0) {
                            continue;
                        }

                        final String whitespaceString = whitespace != 0 ? "\t \t \t " : " ";
                        final String decoratorString = decoratorCount != 0 ? RandomStringUtils.randomAlphanumeric(decoratorCount) : "";

                        testArgs.add(
                                new Object[]{
                                        args[0], args[1], args[2],
                                        refVersions[refidx],
                                        decoratorString,
                                        whitespaceString
                                }
                        );
                    }
                }
            }
        }

        return testArgs.iterator();

    }

    @DataProvider
    private Object[][] provideForGetDataSourceVersionString() {
        return new Object[][] {
                {
                        1,
                        6,
                        LocalDate.of(2019, Month.JANUARY, 24),
                        "v1.6.20190124"
                },
                {
                        1,
                        7,
                        LocalDate.of(2020, Month.MAY, 21),
                        "v1.7.20200521"
                },
                {
                        52,
                        13,
                        LocalDate.of(71999, Month.DECEMBER, 7),
                        "v52.13.719991207"
                },
        };
    }

    //==================================================================================================================
    // Tests:

    @Test(dataProvider = "provideForGetDataSourceVersionString")
    public void testGetDataSourceVersionString(final int major,
                                               final int minor,
                                               final LocalDate date,
                                               final String expected) {
        Assert.assertEquals(
                DataSourceUtils.getDataSourceVersionString(major, minor, date),
                expected
        );
    }

    @Test(dataProvider = "provideForValidateVersionInformation")
    public void testValidateVersionInformation(final Integer major,
                                               final Integer minor,
                                               final LocalDate releaseDate,
                                               final Boolean expected) {
        Assert.assertEquals(
                DataSourceUtils.validateVersionInformation(
                    major,
                    minor,
                    releaseDate.getYear(),
                    releaseDate.getMonthValue(),
                    releaseDate.getDayOfMonth()
                ),
            expected.booleanValue()
        );
    }

    @Test
    public void testGetDataSourceMaxVersionString() {
        Assert.assertEquals(
                DataSourceUtils.getDataSourceMaxVersionString(38),
                DataSourceUtils.getNewDataSourceVersionString(
                        DataSourceUtils.MAX_MAJOR_VERSION_NUMBER,
                        DataSourceUtils.MAX_MINOR_VERSION_NUMBER,
                        38,
                        DataSourceUtils.MAX_DATE
                )
        );
    }

    @Test
    public void testGetDataSourceMinVersionString() {
        Assert.assertEquals(
                DataSourceUtils.getDataSourceMinVersionString(),
                DataSourceUtils.getDataSourceVersionString(
                        DataSourceUtils.MIN_MAJOR_VERSION_NUMBER,
                        DataSourceUtils.MIN_MINOR_VERSION_NUMBER,
                        DataSourceUtils.MIN_DATE
                )
        );
    }

    @Test(dataProvider = "provideForTestNewVersionRegex")
    public void testNewVersionRegex(final Integer major,
                                    final Integer minor,
                                    final LocalDate releaseDate,
                                    final String reference,
                                    final String decorator,
                                    final String leadingWhitespace ) {

        // Construct the string:
        final String versionString = String.format(
                "%s%s%d.%d.%s.%4d%02d%02d%s",
                DataSourceUtils.MANIFEST_VERSION_LINE_START,
                leadingWhitespace,
                major,
                minor,
                reference,
                releaseDate.getYear(),
                releaseDate.getMonthValue(),
                releaseDate.getDayOfMonth(),
                decorator
        );

        final Matcher matcher = DataSourceUtils.NEW_VERSION_PATTERN.matcher(versionString);

        Assert.assertTrue(matcher.matches());

        final Integer versionMajor     = Integer.valueOf(matcher.group(1));
        final Integer versionMinor     = Integer.valueOf(matcher.group(2));
        final Integer versionRef       = Integer.valueOf(matcher.group(3));
        final Integer versionYear      = Integer.valueOf(matcher.group(4));
        final Integer versionMonth     = Integer.valueOf(matcher.group(5));
        final Integer versionDay       = Integer.valueOf(matcher.group(6));
        final String  versionDecorator = matcher.group(7);

        Assert.assertEquals( versionMajor, major );
        Assert.assertEquals( versionMinor, minor );
        Assert.assertEquals( "hg"+versionRef, reference );
        Assert.assertEquals( versionYear.intValue(), releaseDate.getYear() );
        Assert.assertEquals( versionMonth.intValue(), releaseDate.getMonthValue() );
        Assert.assertEquals( versionDay.intValue(), releaseDate.getDayOfMonth() );
        Assert.assertEquals( versionDecorator, decorator );
    }

    @Test(dataProvider = "provideForTestOldVersionRegex")
    public void testOldVersionRegex(final Integer major,
                                 final Integer minor,
                                 final LocalDate releaseDate,
                                 final String decorator,
                                 final String leadingWhitespace ) {

        // Construct the string:
        final String versionString = String.format(
                "%s%s%d.%d.%4d%02d%02d%s",
                DataSourceUtils.MANIFEST_VERSION_LINE_START,
                leadingWhitespace,
                major,
                minor,
                releaseDate.getYear(),
                releaseDate.getMonthValue(),
                releaseDate.getDayOfMonth(),
                decorator
        );

        final Matcher matcher = DataSourceUtils.OLD_VERSION_PATTERN.matcher(versionString);

        Assert.assertTrue(matcher.matches());

        final Integer versionMajor     = Integer.valueOf(matcher.group(1));
        final Integer versionMinor     = Integer.valueOf(matcher.group(2));
        final Integer versionYear      = Integer.valueOf(matcher.group(3));
        final Integer versionMonth     = Integer.valueOf(matcher.group(4));
        final Integer versionDay       = Integer.valueOf(matcher.group(5));
        final String  versionDecorator = matcher.group(6);

        Assert.assertEquals( versionMajor, major );
        Assert.assertEquals( versionMinor, minor );
        Assert.assertEquals( versionYear.intValue(), releaseDate.getYear() );
        Assert.assertEquals( versionMonth.intValue(), releaseDate.getMonthValue() );
        Assert.assertEquals( versionDay.intValue(), releaseDate.getDayOfMonth() );
        Assert.assertEquals( versionDecorator, decorator );
    }

    @Test(groups = {"bucket", "cloud"})
    public void testCurrentDataSourcesAvailable() {

        // Make sure that we can actually access the current min and max versions of the data sources.
        final List<String> dataSourcesBasePaths = Arrays.asList(
                DataSourceUtils.DATA_SOURCES_BUCKET_PATH +
                        DataSourceUtils.DATA_SOURCES_NAME_PREFIX +
                        "." + DataSourceUtils.getDataSourceMinVersionString(),
                DataSourceUtils.DATA_SOURCES_BUCKET_PATH +
                        DataSourceUtils.DATA_SOURCES_NAME_PREFIX +
                        "." + DataSourceUtils.getDataSourceMaxVersionString(38),
                DataSourceUtils.DATA_SOURCES_BUCKET_PATH +
                        DataSourceUtils.DATA_SOURCES_NAME_PREFIX +
                        "." + DataSourceUtils.getDataSourceMaxVersionString(19)
        );
        for (final String basePath : dataSourcesBasePaths) {
            for (final String useCaseModifier : Arrays.asList(DataSourceUtils.DS_SOMATIC_NAME_MODIFIER, DataSourceUtils.DS_GERMLINE_NAME_MODIFIER)) {
                // Data Sources File:
                IOUtils.assertFileIsReadable(
                        IOUtils.getPath(
                            basePath + useCaseModifier + DataSourceUtils.DS_EXTENSION
                        )
                );

                // Checksum for file:
                IOUtils.assertFileIsReadable(
                        IOUtils.getPath(
                                basePath + useCaseModifier + DataSourceUtils.DS_CHECKSUM_EXTENSION
                        )
                );
            }
        }
    }

}
