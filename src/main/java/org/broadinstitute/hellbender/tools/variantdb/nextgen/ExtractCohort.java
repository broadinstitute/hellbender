package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.tools.variantdb.CommonCode;
import org.broadinstitute.hellbender.tools.variantdb.SampleList;

import java.util.*;


@CommandLineProgramProperties(
        summary = "(\"ExtractCohort\") - Filter and extract variants out of big query.",
        oneLineSummary = "Tool to extract variants out of big query for a subset of samples",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
public class ExtractCohort extends ExtractTool {
    private static final Logger logger = LogManager.getLogger(ExtractCohort.class);
    private ExtractCohortEngine engine;

   @Argument(
            fullName = "variant-filter-table",
            doc = "Fully qualified name of the filtering table to use for cohort extraction",
            optional = true
    )
    private String filteringFQTableName = null;

    @Argument(
            fullName = "cohort-extract-table",
            doc = "Fully qualified name of the table where the cohort data exists (already subsetted)",
            optional = false
    )
    private String cohortTable = null;

    @Argument(
            fullName = "min-location",
            doc = "When extracting data, only include locations >= this value",
            optional = true
    )
    private Long minLocation = null;

    @Argument(
        fullName = "max-location",
        doc = "When extracting data, only include locations <= this value",
        optional = true
    )
    private Long maxLocation = null;

    @Override
    protected void onStartup() {
        super.onStartup();

        SampleList sampleList = new SampleList(sampleTableName, null, printDebugInformation);
        Set<String> sampleNames = new HashSet<>(sampleList.getSampleNames());

        VCFHeader header = CommonCode.generateVcfHeader(sampleNames, reference.getSequenceDictionary());

        engine = new ExtractCohortEngine(
                projectID,
                vcfWriter,
                header,
                annotationEngine,
                reference,
                sampleNames,
                mode,
                cohortTable,
                minLocation,
                maxLocation,
                filteringFQTableName,
                localSortMaxRecordsInRam,
                printDebugInformation,
                vqsLodSNPThreshold,
                vqsLodINDELThreshold,
                progressMeter,
                queryMode);
        vcfWriter.writeHeader(header);
    }

    @Override
    // maybe think about creating a BigQuery Row walker?
    public void traverse() {
        progressMeter.setRecordsBetweenTimeChecks(100L);

        if ( filteringFQTableName == null || filteringFQTableName.equals("") ) {
            logger.warn("--variant-filter-table is not specified, no filtering of cohort! ");
        }

        engine.traverse();
     }

    @Override
    protected void onShutdown() {
        super.onShutdown();

        if ( engine != null ) {
            logger.info(String.format("***Processed %d total sites", engine.getTotalNumberOfSites()));
            logger.info(String.format("***Processed %d total variants", engine.getTotalNumberOfVariants()));
        }

        // Close up our writer if we have to:
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }

}
