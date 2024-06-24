package org.broadinstitute.hellbender.tools.funcotator;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.codec.digest.MessageDigestAlgorithms;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.nio.NioFileCopierWithProgressMeter;
import org.broadinstitute.hellbender.utils.nio.NioFileCopierWithProgressMeterResults;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * {@link FuncotatorDataSourceDownloader} is a tool to download the latest data sources for <b><i>{@link Funcotator}</i></b>.
 *
 * <h3>General Information</h3>
 * <p>
 * This tool can download pre-packaged data sources for both the <strong>somatic</strong> and <strong>germline</strong> use cases.
 * The data sources downloaded by this tool correspond to the latest / current maximum of the data sources supported as defined in <b><i>{@link org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils#CURRENT_MAXIMUM_DATA_SOURCE_VERSION}</i></b>.
 * </p>
 *
 * <p>
 * To download and extract the data sources, you can invoke {@link FuncotatorDataSourceDownloader} in the following ways:
 *     <ul>
 *         <li>For <strong>somatic</strong> data sources:<br /><pre>{@code ./gatk FuncotatorDataSourceDownloader --somatic --validate-integrity --hg38 --extract-after-download}</pre></li>
 *         <li>For <strong>germline</strong> data sources:<br /><pre>{@code ./gatk FuncotatorDataSourceDownloader --germline --validate-integrity --hg19 --extract-after-download}</pre></li>
 *     </ul>
 * </p>
 *
 * <h3>Notes</h3>
 * <ul>
 *     <li>By default {@link FuncotatorDataSourceDownloader} will not overwrite data sources if they already exist locally.</li>
 *     <li>It is recommended to run the validation step after downloading the data sources to ensure download integrity.</li>
 *     <li>Using this tool will result in longer download times for data sources than if {@code gsutil} were invoked directly to copy down the data sources.
 *     <br />This is a known issue and is due to the tool printing the progress of the download.</li>
 * </ul>
 */
@CommandLineProgramProperties(
        summary = "Download the Broad Institute pre-packaged data sources for the somatic or germline use case for Funcotator.",
        oneLineSummary = "Data source downloader for Funcotator.",
        programGroup = VariantEvaluationProgramGroup.class
)
@DocumentedFeature
public class FuncotatorDataSourceDownloader extends CommandLineProgram {

    private static final Logger logger = LogManager.getLogger(FuncotatorDataSourceDownloader.class);

    //==================================================================================================================
    // Public Static Members:

    public static final String VALIDATE_INTEGRITY_ARG_LONG_NAME = "validate-integrity";
    public static final String SOMATIC_ARG_LONG_NAME            = "somatic";
    public static final String GERMLINE_ARG_LONG_NAME           = "germline";
    public static final String OVERWRITE_ARG_LONG_NAME          = "overwrite-output-file";
    public static final String EXTRACT_AFTER_DOWNLOAD           = "extract-after-download";
    public static final String HG38_ARG_LONG_NAME               = "hg38";
    public static final String HG19_ARG_LONG_NAME               = "hg19";

    //==================================================================================================================
    // Private Static Members:
    static final String TESTING_OVERRIDE_PATH_FOR_DATA_SOURCES_ARG         = "testing-override-path-for-datasources";
    static final String TESTING_OVERRIDE_PATH_FOR_DATA_SOURCES_SHA256_ARG  = "testing-override-path-for-datasources-sha256";

    //==================================================================================================================
    // Private Static Members:

    // Set to always get the latest version of the data sources:
    private static final String HG38_BASE_URL = DataSourceUtils.DATA_SOURCES_BUCKET_PATH +
            DataSourceUtils.DATA_SOURCES_NAME_PREFIX + "." + DataSourceUtils.getDataSourceMaxVersionString(38);
    private static final String HG19_BASE_URL = DataSourceUtils.DATA_SOURCES_BUCKET_PATH +
            DataSourceUtils.DATA_SOURCES_NAME_PREFIX + "." + DataSourceUtils.getDataSourceMaxVersionString(19);

    private static final String HG38_GERMLINE_GCLOUD_DATASOURCES_BASEURL     = HG38_BASE_URL + DataSourceUtils.DS_GERMLINE_NAME_MODIFIER;
    private static final String HG19_GERMLINE_GCLOUD_DATASOURCES_BASEURL     = HG19_BASE_URL + DataSourceUtils.DS_GERMLINE_NAME_MODIFIER;

    @VisibleForTesting
    static final Path   HG38_GERMLINE_GCLOUD_DATASOURCES_PATH        = IOUtils.getPath(HG38_GERMLINE_GCLOUD_DATASOURCES_BASEURL + DataSourceUtils.DS_EXTENSION);
    static final Path   HG19_GERMLINE_GCLOUD_DATASOURCES_PATH        = IOUtils.getPath(HG19_GERMLINE_GCLOUD_DATASOURCES_BASEURL + DataSourceUtils.DS_EXTENSION);
    private static final Path   HG38_GERMLINE_GCLOUD_DATASOURCES_SHA256_PATH = IOUtils.getPath(HG38_GERMLINE_GCLOUD_DATASOURCES_BASEURL + DataSourceUtils.DS_CHECKSUM_EXTENSION);
    private static final Path   HG19_GERMLINE_GCLOUD_DATASOURCES_SHA256_PATH = IOUtils.getPath(HG19_GERMLINE_GCLOUD_DATASOURCES_BASEURL + DataSourceUtils.DS_CHECKSUM_EXTENSION);

    public static final String HG38_SOMATIC_GCLOUD_DATASOURCES_BASEURL     = HG38_BASE_URL + DataSourceUtils.DS_SOMATIC_NAME_MODIFIER;;
    public static final String HG19_SOMATIC_GCLOUD_DATASOURCES_BASEURL     = HG19_BASE_URL + DataSourceUtils.DS_SOMATIC_NAME_MODIFIER;;

    public static final Path   HG38_SOMATIC_GCLOUD_DATASOURCES_PATH        = IOUtils.getPath(HG38_SOMATIC_GCLOUD_DATASOURCES_BASEURL + DataSourceUtils.DS_EXTENSION);
    public static final Path   HG19_SOMATIC_GCLOUD_DATASOURCES_PATH        = IOUtils.getPath(HG19_SOMATIC_GCLOUD_DATASOURCES_BASEURL + DataSourceUtils.DS_EXTENSION);
    private static final Path   HG38_SOMATIC_GCLOUD_DATASOURCES_SHA256_PATH = IOUtils.getPath(HG38_SOMATIC_GCLOUD_DATASOURCES_BASEURL + DataSourceUtils.DS_CHECKSUM_EXTENSION);
    private static final Path   HG19_SOMATIC_GCLOUD_DATASOURCES_SHA256_PATH = IOUtils.getPath(HG19_SOMATIC_GCLOUD_DATASOURCES_BASEURL + DataSourceUtils.DS_CHECKSUM_EXTENSION);

    //==================================================================================================================
    // Private Members:

    @Argument(fullName = VALIDATE_INTEGRITY_ARG_LONG_NAME,
            shortName = VALIDATE_INTEGRITY_ARG_LONG_NAME,
            doc = "Validate the integrity of the data sources after downloading them using sha256.",
            optional = true)
    private boolean doValidateIntegrity = false;

    @Argument(fullName = SOMATIC_ARG_LONG_NAME,
            shortName = SOMATIC_ARG_LONG_NAME,
            mutex = {GERMLINE_ARG_LONG_NAME, TESTING_OVERRIDE_PATH_FOR_DATA_SOURCES_SHA256_ARG},
            doc = "Download the latest pre-packaged datasources for somatic functional annotation.",
            optional = true)
    private boolean getSomaticDataSources = false;

    @Argument(fullName = GERMLINE_ARG_LONG_NAME,
            shortName = GERMLINE_ARG_LONG_NAME,
            mutex = {SOMATIC_ARG_LONG_NAME, TESTING_OVERRIDE_PATH_FOR_DATA_SOURCES_SHA256_ARG},
            doc = "Download the latest pre-packaged datasources for germline functional annotation.",
            optional = true)
    private boolean getGermlineDataSources = false;

    @Argument(fullName = OVERWRITE_ARG_LONG_NAME,
            shortName = OVERWRITE_ARG_LONG_NAME,
            doc = "Overwrite output file if it exists already.",
            optional = true)
    private boolean overwriteOutputFile = false;

    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName  = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "Output location for the data sources.",
            optional = true)
    protected File outputFile;

    @Argument(
            shortName = EXTRACT_AFTER_DOWNLOAD,
            fullName  = EXTRACT_AFTER_DOWNLOAD,
            doc = "Extract the data sources to a sibling folder after they have been downloaded.",
            optional = true)
    protected boolean extractDataSourcesAfterDownload = false;

    @Argument(
            shortName = HG38_ARG_LONG_NAME,
            fullName  = HG38_ARG_LONG_NAME,
            mutex = {HG19_ARG_LONG_NAME, TESTING_OVERRIDE_PATH_FOR_DATA_SOURCES_SHA256_ARG},
            doc = "If set, will extract data from the HG38 data sources bucket.",
            optional = true)
    protected boolean getHg38Datasources = false;

    @Argument(
            //TODO should these be MUTEX or should one be allowed to download either?
            shortName = HG19_ARG_LONG_NAME,
            fullName  = HG19_ARG_LONG_NAME,
            mutex = {HG38_ARG_LONG_NAME, TESTING_OVERRIDE_PATH_FOR_DATA_SOURCES_SHA256_ARG},
            doc = "If set, will extract data from the HG19 data sources bucket.",
            optional = true)
    protected boolean getHg19Datasources = false;

    // Testing arguments:
    @Hidden
    @Advanced
    @Argument(
            shortName = TESTING_OVERRIDE_PATH_FOR_DATA_SOURCES_ARG,
            fullName  = TESTING_OVERRIDE_PATH_FOR_DATA_SOURCES_ARG,
            mutex = {SOMATIC_ARG_LONG_NAME, GERMLINE_ARG_LONG_NAME},
            doc = "FOR TESTING ONLY: Override path to data sources file with another path.",
            optional = true)
    private String testingOverrideDataSourcesPath;

    @Hidden
    @Advanced
    @Argument(
            shortName = TESTING_OVERRIDE_PATH_FOR_DATA_SOURCES_SHA256_ARG,
            fullName  = TESTING_OVERRIDE_PATH_FOR_DATA_SOURCES_SHA256_ARG,
            mutex = {SOMATIC_ARG_LONG_NAME, GERMLINE_ARG_LONG_NAME},
            doc = "FOR TESTING ONLY: Override path to data sources sha256sum file with another path.",
            optional = true)
    private String testingOverrideDataSourcesSha256Path;

    //==================================================================================================================
    // Constructors:

    //==================================================================================================================
    // Override Methods:

    @Override
    protected void onStartup() {

        // Make sure the user specified at least one data source to download:
        if ((!getSomaticDataSources) && (!getGermlineDataSources) && (testingOverrideDataSourcesPath == null)) {
            throw new UserException("Must select either somatic or germline datasources.");
        }

        // Make sure the user specified at least one reference source to download:
        if ((!getHg38Datasources) && (!getHg19Datasources) && (testingOverrideDataSourcesPath == null)) {
            throw new UserException("Must select either HG19 or HG38 datasources.");
        }

        // Make sure the testing inputs are correct:
        if ( ((testingOverrideDataSourcesPath == null) && (testingOverrideDataSourcesSha256Path != null)) ||
             ((testingOverrideDataSourcesSha256Path == null) && (testingOverrideDataSourcesPath != null)) ) {
            throw new UserException("Must specify both a test data sources path and a test data sources sha256sum path.");
        }

        if ( overwriteOutputFile ) {
            logger.info("Overwrite ENABLED.  Will overwrite existing data sources download.");
        }
    }

    @Override
    protected Object doWork() {

        final String dataSourceDescription;
        final Path dataSourcesPath;
        final Path dataSourcesSha256Path;

        // Get the correct data source:
        if ( getSomaticDataSources ) {
            if (getHg38Datasources) {
                dataSourceDescription = "HG38_Somatic";
                dataSourcesPath = HG38_SOMATIC_GCLOUD_DATASOURCES_PATH;
                dataSourcesSha256Path = HG38_SOMATIC_GCLOUD_DATASOURCES_SHA256_PATH;
            } else { // Okay because HG38 and HG19 datasources are currently MUTEX and at least one is required
                dataSourceDescription = "HG19_Somatic";
                dataSourcesPath = HG19_SOMATIC_GCLOUD_DATASOURCES_PATH;
                dataSourcesSha256Path = HG19_SOMATIC_GCLOUD_DATASOURCES_SHA256_PATH;
            }
        }
        else if ( getGermlineDataSources ) {
            if (getHg38Datasources) {
                dataSourceDescription = "HG38_Germline";
                dataSourcesPath = HG38_GERMLINE_GCLOUD_DATASOURCES_PATH;
                dataSourcesSha256Path = HG38_GERMLINE_GCLOUD_DATASOURCES_SHA256_PATH;
            } else {
                dataSourceDescription = "HG19_Germline";
                dataSourcesPath = HG19_GERMLINE_GCLOUD_DATASOURCES_PATH;
                dataSourcesSha256Path = HG19_GERMLINE_GCLOUD_DATASOURCES_SHA256_PATH;
            }
        }
        else {
            // Test case:
            dataSourceDescription = "TESTING";
            dataSourcesPath = IOUtils.getPath(testingOverrideDataSourcesPath);
            dataSourcesSha256Path = IOUtils.getPath(testingOverrideDataSourcesSha256Path);
        }

        downloadAndValidateDataSources(dataSourceDescription, dataSourcesPath, dataSourcesSha256Path);

        // Token return value:
        return true;
    }

    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    private NioFileCopierWithProgressMeter createNioDownloader(final Path dsPath) {
        // Get the data sources file:
        final Path outputDestination = getOutputLocation(dsPath);

        // Set up and initiate our download:
        return NioFileCopierWithProgressMeter.create(dsPath, outputDestination, overwriteOutputFile);
    }

    private void downloadAndValidateDataSources(final String dataSourceType, final Path dsPath, final Path checksumPath) {
        logger.info(dataSourceType + " data sources selected.");

        // Create our downloader:
        final NioFileCopierWithProgressMeter xerox = createNioDownloader(dsPath);

        // Get the datasources file:
        final NioFileCopierWithProgressMeterResults results = downloadDataSources(xerox, checksumPath);

        // Confirm file integrity if requested:
        if ( doValidateIntegrity ) {
            validateIntegrity(results);
        }

        // Extract data sources if requested:
        if ( extractDataSourcesAfterDownload ) {
            IOUtils.extractTarGz(results.getDestination(), results.getDestination().getParent(), overwriteOutputFile);
        }
        else {
            logger.info("IMPORTANT: You must unzip the downloaded data sources prior to using them with Funcotator.");
        }
    }

    private NioFileCopierWithProgressMeterResults downloadDataSources(final NioFileCopierWithProgressMeter xerox,
                                                                      final Path checksumPath) {

        // Set up the validity check while the download occurs if requested:
        if ( doValidateIntegrity ) {
            // Read the sha256sum into memory:
            final String expectedSha256Sum = readSha256SumFromPath(checksumPath);

            // Setup the copier to calculate the checksum:
            xerox.setChecksumAlgorithmAndExpectedChecksum(MessageDigestAlgorithms.SHA_256, expectedSha256Sum);
        }

        // Initiate the copy and return our results:
        return xerox.initiateCopy();
    }

    private void validateIntegrity(final NioFileCopierWithProgressMeterResults results) {

        // verify the hashes are the same:
        if ( !results.isDestFileValid() ) {
            throw new UserException("ERROR: downloaded data sources are corrupt!  Unexpected checksum: " + results.getChecksum() + " != " + results.getExpectedChecksum());
        }
        else {
            logger.info("Integrity check on downloaded data sources succeeded.");
        }
    }

    private String readSha256SumFromPath(final Path sha256Path) {
        final String expectedSha256Sum;
        try {
            logger.info("Collecting expected checksum...");
            expectedSha256Sum = Files.lines(sha256Path).findFirst().orElse(null);

            if ( expectedSha256Sum == null ) {
                throw new UserException("Unable to retrieve expected checksum from: " + sha256Path.toUri());
            }

            logger.info("Collection complete!");
        }
        catch ( final IOException ex ) {
            throw new UserException("Could not read in sha256sum from file: " + sha256Path.toUri().toString(), ex);
        }

        // Clean up and return the checksum:
        return cleanExpectedSha256SumString(expectedSha256Sum);
    }

    private Path getOutputLocation(final Path dataSourcesPath) {
        if ( outputFile == null ) {
            return IOUtils.getPath(dataSourcesPath.getFileName().toString());
        }
        else {
            return outputFile.toPath();
        }
    }

    private String cleanExpectedSha256SumString(final String expectedSha256SumString) {

        String cleanString = expectedSha256SumString.trim().toLowerCase();

        // The format of the file can contain the checksum, followed by the file name.
        // If this is the case, we need to truncate the string:
        if ( cleanString.contains(" ") ) {
            cleanString = cleanString.substring(0, cleanString.indexOf(" "));
        }
        if ( cleanString.contains("\t")) {
            cleanString = cleanString.substring(0, cleanString.indexOf("\t"));
        }

        return cleanString;
    }

    //==================================================================================================================
    // Helper Data Types:

}
