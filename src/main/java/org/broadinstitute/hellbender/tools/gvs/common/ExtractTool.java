package org.broadinstitute.hellbender.tools.gvs.common;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.*;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public abstract class ExtractTool extends GATKTool {
    public static final int DEFAULT_LOCAL_SORT_MAX_RECORDS_IN_RAM = 1000000;
    protected VariantAnnotatorEngine annotationEngine;

    @Argument(
            fullName = "project-id",
            doc = "ID of the Google Cloud project to use when executing queries",
            optional = true
    )
    protected String projectID = null;

    @Argument(
            fullName = "dataset-id",
            doc = "ID of the Google Cloud dataset to use when executing queries",
            optional = true // I guess, but won't it break otherwise or require that a dataset be created with the name temp_tables?
    )
    protected String datasetID = null;

    @Argument(
            fullName = "sample-table",
            doc = "Fully qualified name of a bigquery table containing a single column `sample` that describes the full list of samples to extract",
            optional = true,
            mutex={"sample-file"}
    )
    protected String sampleTableName = null;

    @Argument(
            fullName = "sample-file",
            doc = "Alternative to `sample-table`. Pass in a (sample_id,sample_name) CSV that describes the full list of samples to extract. No header",
            optional = true,
            mutex={"sample-table"}

    )
    protected File sampleFileName = null;

    @Argument(
            fullName = "print-debug-information",
            doc = "If true, print extra debugging output",
            optional = true)
    protected boolean printDebugInformation = false;

    @Argument(
            fullName = "local-sort-max-records-in-ram",
            doc = "When doing local sort, store at most this many records in memory at once",
            optional = true
    )
    protected int localSortMaxRecordsInRam = DEFAULT_LOCAL_SORT_MAX_RECORDS_IN_RAM;

    @Argument(
            fullName = "ref-version",
            doc = "Remove this option!!!! only for ease of testing. Valid options are 37 or 38",
            optional = true
    )
    protected String refVersion = "37";

    @Argument(
        fullName = "min-location",
        doc = "When extracting data, only include locations >= this value",
        optional = true
    )
    protected Long minLocation = null;

    @Argument(
        fullName = "max-location",
        doc = "When extracting data, only include locations <= this value",
        optional = true
    )
    protected Long maxLocation = null;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public boolean useVariantAnnotations() {
        return true;
    }

    @Override
    public List<Annotation> getDefaultVariantAnnotations() {
        return Arrays.asList(
                // All the `AS_StandardAnnotation` implementers minus `AS_InbreedingCoeff`.
                new AS_FisherStrand(),
                new AS_StrandOddsRatio(),
                new AS_BaseQualityRankSumTest(),
                new AS_MappingQualityRankSumTest(),
                new AS_ReadPosRankSumTest(),
                new AS_RMSMappingQuality(),
                new AS_QualByDepth(),

                // All the `StandardAnnotation` implementers minus `InbreedingCoeff` and `ExcessHet`.
                new BaseQualityRankSumTest(),
                new ChromosomeCounts(),
                new Coverage(),
                new DepthPerAlleleBySample(),
                new FisherStrand(),
                new MappingQualityRankSumTest(),
                new QualByDepth(),
                new RMSMappingQuality(),
                new ReadPosRankSumTest(),
                new StrandOddsRatio()
        );
    }

    @Override
    protected void onStartup() {
        super.onStartup();

        //TODO verify what we really need here
        annotationEngine = new VariantAnnotatorEngine(makeVariantAnnotations(), null, Collections.emptyList(), false, false);

        ChromosomeEnum.setRefVersion(refVersion);
    }
}
