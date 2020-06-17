package org.broadinstitute.hellbender.tools.variantdb;

import com.google.cloud.bigquery.FieldValueList;
import com.google.cloud.bigquery.TableResult;
import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.avro.generic.GenericRecord;
import org.apache.commons.lang.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.ProgressMeter;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.variantdb.RawArrayData.ArrayGenotype;
import org.broadinstitute.hellbender.tools.walkers.ReferenceConfidenceVariantContextMerger;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bigquery.*;
import org.broadinstitute.hellbender.utils.localsort.AvroSortingCollectionCodec;
import org.broadinstitute.hellbender.utils.localsort.SortingCollection;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.Collectors;
import static org.broadinstitute.hellbender.tools.variantdb.ExtractCohortBQ.*;


public class ArrayExtractCohortEngine {
    private final DecimalFormat df = new DecimalFormat();
    private final String DOT = ".";
    
    static {
    }
    
    private static final Logger logger = LogManager.getLogger(ExtractCohortEngine.class);

    private final VariantContextWriter vcfWriter;

    private final boolean useCompressedData;
    private final boolean printDebugInformation;
    private final int localSortMaxRecordsInRam;
    private final TableReference cohortTableRef;
    private final ReferenceDataSource refSource;

    private final ProgressMeter progressMeter;
    private final String projectID;

    /** List of sample names seen in the variant data from BigQuery. */
    private final Map<Integer, String> sampleIdMap;
    private final Set<String> sampleNames;

    private final Map<Long, ProbeInfo> probeIdMap;
    private final ReferenceConfidenceVariantContextMerger variantContextMerger;

    private int totalNumberOfVariants = 0;
    private int totalNumberOfSites = 0;

    /**
     * The conf threshold above which variants are not included in the position tables.
     * This value is used to construct the genotype information of those missing samples
     * when they are merged together into a {@link VariantContext} object
     */
    public static int MISSING_CONF_THRESHOLD = 60;


    public ArrayExtractCohortEngine(final String projectID,
                                    final VariantContextWriter vcfWriter,
                                    final VCFHeader vcfHeader,
                                    final VariantAnnotatorEngine annotationEngine,
                                    final ReferenceDataSource refSource,
                                    final Map<Integer, String> sampleIdMap,
                                    final Map<Long, ProbeInfo> probeIdMap,
                                    final String cohortTableName,
                                    final int localSortMaxRecordsInRam,
                                    final boolean useCompressedData,
                                    final boolean printDebugInformation,
                                    final ProgressMeter progressMeter) {

        this.df.setMaximumFractionDigits(3);
        this.df.setGroupingSize(0);
                                
        this.localSortMaxRecordsInRam = localSortMaxRecordsInRam;

        this.projectID = projectID;
        this.vcfWriter = vcfWriter;
        this.refSource = refSource;
        this.sampleIdMap = sampleIdMap;
        this.sampleNames = new HashSet<>(sampleIdMap.values());

        this.probeIdMap = probeIdMap;

        this.cohortTableRef = new TableReference(cohortTableName, useCompressedData?SchemaUtils.RAW_ARRAY_COHORT_FIELDS_COMPRESSED:SchemaUtils.RAW_ARRAY_COHORT_FIELDS_UNCOMPRESSED);

        this.useCompressedData = useCompressedData;
        this.printDebugInformation = printDebugInformation;
        this.progressMeter = progressMeter;

        // KCIBUL: what is the right variant context merger for arrays?
        this.variantContextMerger = new ReferenceConfidenceVariantContextMerger(annotationEngine, vcfHeader);

    }

    int getTotalNumberOfVariants() { return totalNumberOfVariants; }
    int getTotalNumberOfSites() { return totalNumberOfSites; }

    public void traverse() {
        if (printDebugInformation) {
            logger.debug("using storage api with local sort");
        }
        final StorageAPIAvroReader storageAPIAvroReader = new StorageAPIAvroReader(cohortTableRef);
        createVariantsFromUngroupedTableResult(storageAPIAvroReader);
    }


    private void createVariantsFromUngroupedTableResult(final GATKAvroReader avroReader) {

        // stream out the data and sort locally
        final org.apache.avro.Schema schema = avroReader.getSchema();
        final Set<String> columnNames = new HashSet<>();
        schema.getFields().forEach(field -> columnNames.add(field.name()));

        Comparator<GenericRecord> comparator = this.useCompressedData ? COMPRESSED_PROBE_ID_COMPARATOR : UNCOMPRESSED_PROBE_ID_COMPARATOR;
        SortingCollection<GenericRecord> sortingCollection =  getAvroProbeIdSortingCollection(schema, localSortMaxRecordsInRam, comparator);
        for ( final GenericRecord queryRow : avroReader ) {
            sortingCollection.add(queryRow);
        }

        sortingCollection.printTempFileStats();

        // iterate through records and process them
        final List<GenericRecord> currentPositionRecords = new ArrayList<>(sampleIdMap.size() * 2);
        long currentProbeId = -1;

        for ( final GenericRecord sortedRow : sortingCollection ) {
            long probeId;
            if (useCompressedData) {
                final long rawData = (Long) sortedRow.get(SchemaUtils.RAW_ARRAY_DATA_FIELD_NAME);
                RawArrayData data = RawArrayData.decode(rawData);
                probeId = data.probeId;
            } else {
                probeId = (Long) sortedRow.get("probe_id");
            }

            if ( probeId != currentProbeId && currentProbeId != -1 ) {
                ++totalNumberOfSites;
                processSampleRecordsForLocation(currentProbeId, currentPositionRecords, columnNames);
                currentPositionRecords.clear();
            }

            currentPositionRecords.add(sortedRow);
            currentProbeId = probeId;
        }

        if ( ! currentPositionRecords.isEmpty() ) {
            ++totalNumberOfSites;
            processSampleRecordsForLocation(currentProbeId, currentPositionRecords, columnNames);
        }
    }

    private void processSampleRecordsForLocation(final long probeId, final Iterable<GenericRecord> sampleRecordsAtPosition, final Set<String> columnNames) {
        final List<VariantContext> unmergedCalls = new ArrayList<>();
        final Set<String> currentPositionSamplesSeen = new HashSet<>();
        boolean currentPositionHasVariant = false;

        final ProbeInfo probeInfo = probeIdMap.get(probeId);
        if (probeInfo == null) {
            throw new RuntimeException("Unable to find probeInfo for " + probeId);
        }

        final String contig = probeInfo.contig;
        final long position = probeInfo.position;
        final Allele refAllele = Allele.create(refSource.queryAndPrefetch(contig, position, position).getBaseString(), true);

        int numRecordsAtPosition = 0;

        for ( final GenericRecord sampleRecord : sampleRecordsAtPosition ) {
            final long sampleId = (Long) sampleRecord.get(SchemaUtils.SAMPLE_ID_FIELD_NAME);

            // TODO: handle missing values
            String sampleName = sampleIdMap.get((int) sampleId);            
            currentPositionSamplesSeen.add(sampleName);

            ++numRecordsAtPosition;

            if ( printDebugInformation ) {
                logger.info("\t" + contig + ":" + position + ": found record for sample " + sampleName + ": " + sampleRecord);
            }

            ++totalNumberOfVariants;
            unmergedCalls.add(createVariantContextFromSampleRecord(probeInfo, sampleRecord, columnNames, contig, position, sampleName));

        }

        if ( printDebugInformation ) {
            logger.info(contig + ":" + position + ": processed " + numRecordsAtPosition + " total sample records");
        }

        finalizeCurrentVariant(unmergedCalls, currentPositionSamplesSeen, contig, position, refAllele);
    }

    private void finalizeCurrentVariant(final List<VariantContext> unmergedCalls, final Set<String> currentVariantSamplesSeen, final String contig, final long start, final Allele refAllele) {

        // TODO: this is where we infer missing data points... once we know what we want to drop
        // final Set<String> samplesNotEncountered = Sets.difference(sampleNames, currentVariantSamplesSeen);
        // for ( final String missingSample : samplesNotEncountered ) {
        //         unmergedCalls.add(createRefSiteVariantContext(missingSample, contig, start, refAllele));
        // }

        final VariantContext mergedVC = variantContextMerger.merge(
                unmergedCalls,
                new SimpleInterval(contig, (int) start, (int) start),
                refAllele.getBases()[0],
                true,
                false,
                true);


        final VariantContext finalVC = mergedVC;

        // TODO: this was commented out... probably need to re-enable
//        final VariantContext annotatedVC = enableVariantAnnotator ?
//                variantAnnotator.annotateContext(finalizedVC, new FeatureContext(), null, null, a -> true): finalVC;

//        if ( annotatedVC != null ) {
//            vcfWriter.add(annotatedVC);
//            progressMeter.update(annotatedVC);
//        }

        if ( finalVC != null ) {
            vcfWriter.add(finalVC);
            progressMeter.update(finalVC);
        } else {
            // TODO should i print a warning here?
            vcfWriter.add(mergedVC);
            progressMeter.update(mergedVC);
        }
    }

    private String formatFloatForVcf(final Float value) {
        if (value == null || Double.isNaN(value)) {
            return DOT;
        }
        return df.format(value);
    }

    private Float getNullableFloatFromDouble(Object d) {
        return d == null ? null : (float)  ((Double) d).doubleValue();
    }

    private VariantContext createVariantContextFromSampleRecord(final ProbeInfo probeInfo, final GenericRecord sampleRecord, final Set<String> columnNames, final String contig, final long startPosition, final String sample) {
        final VariantContextBuilder builder = new VariantContextBuilder();
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder();

        builder.chr(contig);
        builder.start(startPosition);

        
        final List<Allele> alleles = new ArrayList<>();
        Allele ref = Allele.create(probeInfo.ref, true);        
        alleles.add(ref);

        Allele alleleA = Allele.create(probeInfo.alleleA, false);
        Allele alleleB = Allele.create(probeInfo.alleleB, false);

        boolean alleleAisRef = probeInfo.ref.equals(probeInfo.alleleA);
        boolean alleleBisRef = probeInfo.ref.equals(probeInfo.alleleB);

        if (alleleAisRef) {
            alleleA = ref;
        } else {
            alleles.add(alleleA);
        }

        if (alleleBisRef) {
            alleleB = ref;
        } else {
            alleles.add(alleleB);
        }

        builder.alleles(alleles);
        builder.stop(startPosition + alleles.get(0).length() - 1);

        Float normx;
        Float normy;
        Float baf;
        Float lrr;
        List<Allele> genotypeAlleles = new ArrayList<Allele>();

        if (this.useCompressedData) {
            final RawArrayData data = RawArrayData.decode((Long) sampleRecord.get(SchemaUtils.RAW_ARRAY_DATA_FIELD_NAME));
            normx = data.normx;
            normy = data.normy;
            lrr = data.lrr;
            baf = data.baf;

            // Genotype -- what about no-call?
            if (data.genotype == ArrayGenotype.AA) {
                genotypeAlleles.add(alleleA);
                genotypeAlleles.add(alleleA);
            } else if (data.genotype == ArrayGenotype.AB) {
                genotypeAlleles.add(alleleA);
                genotypeAlleles.add(alleleB);
            } else if  (data.genotype == ArrayGenotype.BB) {
                genotypeAlleles.add(alleleB);
                genotypeAlleles.add(alleleB);
            }
        } else {
            // TODO: constantize
            try {
                normx = getNullableFloatFromDouble(sampleRecord.get("call_NORMX"));
                normy = getNullableFloatFromDouble(sampleRecord.get("call_NORMY"));            
                baf = getNullableFloatFromDouble(sampleRecord.get("call_BAF"));
                lrr = getNullableFloatFromDouble(sampleRecord.get("call_LRR"));
            } catch (NullPointerException npe) {
                System.out.println("NPE on " + sampleRecord);
                System.out.println("NPE on BAF " + sampleRecord.get("call_BAF"));
                System.out.println("NPE on LRR " +sampleRecord.get("call_LRR"));
                throw npe;
            }

            Object gt = sampleRecord.get("call_GT_encoded)");
            // TODO: Genotype -- what about no-call?
            if (gt == null || gt.toString().length() == 0) {
                genotypeAlleles.add(alleleA);
                genotypeAlleles.add(alleleA);
            } else if ("X".equals(gt.toString())) {
                genotypeAlleles.add(alleleA);
                genotypeAlleles.add(alleleB);
            } else if ("B".equals(gt.toString())) {
                genotypeAlleles.add(alleleB);
                genotypeAlleles.add(alleleB);
            } else {
                System.out.println("Processing getnotype " + gt.toString());
            }
        }
        genotypeBuilder.alleles(genotypeAlleles);

        genotypeBuilder.attribute(CommonCode.NORMX, formatFloatForVcf(normx));
        genotypeBuilder.attribute(CommonCode.NORMY, formatFloatForVcf(normy));
        genotypeBuilder.attribute(CommonCode.BAF, formatFloatForVcf(baf));
        genotypeBuilder.attribute(CommonCode.LRR, formatFloatForVcf(lrr));      

        genotypeBuilder.name(sample);

        builder.genotypes(genotypeBuilder.make());

        try {
            VariantContext vc = builder.make();
            return vc;
        } catch (Exception e) {
            System.out.println("Error: "+ e.getMessage() + " processing " + sampleRecord + " and ref: " +ref + " PI: " + probeInfo.alleleA + "/" +probeInfo.alleleB + " with ga " + genotypeAlleles + " and alleles " + alleles);
            throw e;
        }
    }
}
