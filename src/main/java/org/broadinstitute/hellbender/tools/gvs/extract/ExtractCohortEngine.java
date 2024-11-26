package org.broadinstitute.hellbender.tools.gvs.extract;

import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.avro.generic.GenericRecord;
import org.apache.commons.lang.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.gvs.common.*;
import org.broadinstitute.hellbender.tools.walkers.ReferenceConfidenceVariantContextMerger;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gvs.bigquery.AvroFileReader;
import org.broadinstitute.hellbender.utils.gvs.bigquery.StorageAPIAvroReader;
import org.broadinstitute.hellbender.utils.gvs.bigquery.TableReference;
import org.broadinstitute.hellbender.utils.gvs.localsort.AvroSortingCollectionCodec;
import org.broadinstitute.hellbender.utils.gvs.localsort.SortingCollection;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;
import java.util.function.Consumer;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class ExtractCohortEngine {
    private final Logger logger;

    private final boolean printDebugInformation;
    private final int localSortMaxRecordsInRam;
    private final List<SimpleInterval> traversalIntervals;
    private final Long minLocation;
    private final Long maxLocation;
    private final TableReference filterSetInfoTableRef;
    private final TableReference filterSetSiteTableRef;
    private final TableReference samplePloidyTableRef;
    private final ReferenceDataSource refSource;
    private final Double vqScoreSNPThreshold;
    private final Double vqScoreINDELThreshold;
    private final ExtractCohort.VQScoreFilteringType vqScoreFilteringType;
    private final boolean convertFilteredGenotypesToNoCalls;
    private final OptionalLong maxAlternateAlleleCount;

    private final String projectID;
    private final boolean emitPLs;
    private final boolean emitADs;

    // Set of sample ids seen in the variant data from BigQuery.
    private final SortedSet<Long> sampleIdsToExtract;
    private final BitSet sampleIdsToExtractBitSet;

    private final Map<Long, String> sampleIdToName;

    private final ReferenceConfidenceVariantContextMerger variantContextMerger;
    private final VariantAnnotatorEngine annotationEngine;

    private int totalNumberOfVariants = 0;
    private int totalNumberOfSites = 0;
    private int totalRangeRecords = 0;
    private int totalIrrelevantRangeRecords = 0;
    private long totalEstimatedBytesScanned = 0;

    private final String vetRangesFQDataSet;

    private final String fqRangesExtractVetTable;
    private final List<String> extractVetFields;
    private final String fqRangesExtractRefTable;

    private final GATKPath vetAvroFileName;
    private final GATKPath refRangesAvroFileName;
    private final String filterSetName;

    private final GQStateEnum inferredReferenceState;
    private final InferredReferenceRecord inferredReferenceRecord;

    private final boolean presortedAvroFiles;

    private final Consumer<VariantContext> variantContextConsumer;

    private static class VariantIterables {
        public Iterable<GenericRecord> vets;
        public Iterable<GenericRecord> refRanges;
        public VariantIterables(Iterable<GenericRecord> vets, Iterable<GenericRecord> refRanges) {
            this.vets = vets;
            this.refRanges = refRanges;
        }
    }

    List<String> getFilterSetInfoTableFields() {
        return SchemaUtils.YNG_FIELDS;
    }

    String getScoreFieldName() {
        return null;
    }

    String getScoreKey() {
        return null;
    }

    String getVQScoreFieldName() {
        return SchemaUtils.VQSLOD;
    }

    String getAlleleSpecificVQSScoreKey() {
        return GATKVCFConstants.AS_VQS_LOD_KEY;
    }

    String getVqScoreSNPFailureFilterName() {
        return GATKVCFConstants.VQSR_FAILURE_SNP;
    }

    String getVqScoreINDELFailureFilterName() {
        return GATKVCFConstants.VQSR_FAILURE_INDEL;
    }

    public ExtractCohortEngine(final String projectID,
                               final VCFHeader vcfHeader,
                               final VariantAnnotatorEngine annotationEngine,
                               final ReferenceDataSource refSource,
                               final Map<Long, String> sampleIdToName,
                               final String vetRangesFQDataSet,
                               final String fqRangesExtractVetTable,
                               final VetRangesExtractVersionEnum vetRangesExtractTableVersion,
                               final String fqRangesExtractRefTable,
                               final GATKPath vetAvroFileName,
                               final GATKPath refRangesAvroFileName,
                               final List<SimpleInterval> traversalIntervals,
                               final Long minLocation,
                               final Long maxLocation,
                               final String filterSetInfoTableName,
                               final String filterSetSiteTableName,
                               final String samplePloidyTableName,
                               final int localSortMaxRecordsInRam,
                               final boolean printDebugInformation,
                               final Double vqScoreSNPThreshold,
                               final Double vqScoreINDELThreshold,
                               final String filterSetName,
                               final boolean emitPLs,
                               final boolean emitADs,
                               final ExtractCohort.VQScoreFilteringType vqScoreFilteringType,
                               final boolean convertFilteredGenotypesToNoCalls,
                               final OptionalLong maxAlternateAlleleCount,
                               final GQStateEnum inferredReferenceState,
                               final boolean presortedAvroFiles,
                               final Consumer<VariantContext> variantContextConsumer
    ) {
        this.localSortMaxRecordsInRam = localSortMaxRecordsInRam;

        this.projectID = projectID;
        this.refSource = refSource;
        this.sampleIdToName = sampleIdToName;
        this.sampleIdsToExtract = new TreeSet<>(this.sampleIdToName.keySet());

        long maxSampleId = sampleIdsToExtract.last();

        // Constraint due to easily available BitSet implementations, could
        // be extended in future work
        if (maxSampleId > Integer.MAX_VALUE) {
            throw new GATKException("Sample Ids > " + Integer.MAX_VALUE + " are not supported");
        }

        this.sampleIdsToExtractBitSet = new BitSet((int) maxSampleId + 1);
        for (Long id : sampleIdsToExtract) {
            sampleIdsToExtractBitSet.set(id.intValue());
        }

        this.emitPLs = emitPLs;
        this.emitADs = emitADs;

        this.vetRangesFQDataSet = vetRangesFQDataSet;
        this.fqRangesExtractVetTable = fqRangesExtractVetTable;

        if (vetRangesExtractTableVersion == VetRangesExtractVersionEnum.V2) {
            this.extractVetFields = SchemaUtils.EXTRACT_VET_V2_FIELDS;
        }
        else if (vetRangesExtractTableVersion == VetRangesExtractVersionEnum.V1) {
            this.extractVetFields = SchemaUtils.EXTRACT_VET_V1_FIELDS;
        }
        else {
            // This can't happen.
            throw new GATKException("Unknown vet_version " + vetRangesExtractTableVersion + " supplied");
        }
        this.fqRangesExtractRefTable = fqRangesExtractRefTable;

        this.vetAvroFileName = vetAvroFileName;
        this.refRangesAvroFileName = refRangesAvroFileName;

        this.traversalIntervals = traversalIntervals;
        this.minLocation = minLocation;
        this.maxLocation = maxLocation;

        this.printDebugInformation = printDebugInformation;
        this.vqScoreSNPThreshold = vqScoreSNPThreshold;
        this.vqScoreINDELThreshold = vqScoreINDELThreshold;
        this.vqScoreFilteringType = vqScoreFilteringType;
        this.convertFilteredGenotypesToNoCalls = convertFilteredGenotypesToNoCalls;
        this.maxAlternateAlleleCount = maxAlternateAlleleCount;

        this.filterSetSiteTableRef = vqScoreFilteringType.equals(ExtractCohort.VQScoreFilteringType.NONE) ? null : new TableReference(filterSetSiteTableName, SchemaUtils.FILTER_SET_SITE_FIELDS);
        this.filterSetInfoTableRef = vqScoreFilteringType.equals(ExtractCohort.VQScoreFilteringType.NONE) ? null : new TableReference(filterSetInfoTableName, getFilterSetInfoTableFields());
        this.samplePloidyTableRef = (samplePloidyTableName == null? null : new TableReference(samplePloidyTableName, SchemaUtils.SAMPLE_PLOIDY_FIELDS));

        this.filterSetName = filterSetName;

        this.annotationEngine = annotationEngine;
        this.variantContextMerger = new ReferenceConfidenceVariantContextMerger(
                annotationEngine, vcfHeader, false, false, true);

        this.inferredReferenceState = inferredReferenceState;

        this.presortedAvroFiles = presortedAvroFiles;

        this.inferredReferenceRecord = new InferredReferenceRecord(inferredReferenceState);
        this.variantContextConsumer = variantContextConsumer;
        logger = LogManager.getLogger(ExtractCohortEngine.class);
    }

    int getTotalNumberOfVariants() {
        return totalNumberOfVariants;
    }

    int getTotalNumberOfSites() {
        return totalNumberOfSites;
    }

    long getTotalEstimatedBytesScanned() {
        return totalEstimatedBytesScanned;
    }

    private void processBytesScanned(StorageAPIAvroReader reader) {
        long bytes = reader.getEstimatedTotalBytesScanned();
        totalEstimatedBytesScanned += bytes;
    }

    public void traverse() {

        logger.info("minLocation: " + minLocation + ", maxLocation: " + maxLocation);

        if (!SchemaUtils.decodeContig(minLocation).equals(SchemaUtils.decodeContig(maxLocation))) {
            throw new GATKException("Can not process cross-contig boundaries for Ranges implementation");
        }

        String ploidyTableRestriction = null;
        final Map<String, Integer> samplePloidyMap = new HashMap<>();

        if (minLocation != null && maxLocation != null) {
            // This block of code already requires minLocation and maxLocation to be on the same chromosome, so we can rely on that here
            ploidyTableRestriction = "chromosome = " + ((minLocation/SchemaUtils.chromAdjustment) * SchemaUtils.chromAdjustment);
        }

        if (samplePloidyTableRef != null) {
            try (StorageAPIAvroReader reader = new StorageAPIAvroReader(samplePloidyTableRef, ploidyTableRestriction, projectID)) {
                logger.info("Found ploidy lookup table.  Reading it into memory");
                for (final GenericRecord queryRow : reader) {
                    // the key will be a basic joining of chromosome (represented as a location) with sample name
                    String chromosome = queryRow.get(SchemaUtils.CHROMOSOME).toString();
                    String sampleId = queryRow.get(SchemaUtils.SAMPLE_ID).toString();
                    Integer ploidy = Integer.parseInt(queryRow.get(SchemaUtils.PLOIDY).toString());
                    samplePloidyMap.put(makePloidyLookupKey(chromosome, sampleId), ploidy);
                }
                processBytesScanned(reader);
                logger.info("Finished reading ploidy table into memory. " + samplePloidyMap.size() + " entries read.");
            }
        }

        SortedSet<Long> sampleIdsToExtract = new TreeSet<>(this.sampleIdToName.keySet());
        VariantBitSet vbs = new VariantBitSet(minLocation, maxLocation);
        VariantIterables variantIterables;
        if (fqRangesExtractVetTable != null) {
            variantIterables = createVariantIterablesFromUnsortedExtractTableBigQueryRanges(fqRangesExtractVetTable, fqRangesExtractRefTable, vbs);
        } else if (vetRangesFQDataSet != null) {
            variantIterables = createVariantIterablesFromUnsortedBigQueryRanges(vetRangesFQDataSet, sampleIdsToExtract, vbs);
        } else {
            variantIterables = createVariantIterablesFromUnsortedAvroRanges(vetAvroFileName, refRangesAvroFileName, vbs, presortedAvroFiles);
        }

        // First allele here is the ref, followed by the alts associated with that ref. We need this because at this
        // point the alleles haven't been joined and remapped to one reference allele.
        final Map<Long, Map<Allele, Map<Allele, Double>>> fullScoreMap = new HashMap<>();
        final Map<Long, Map<Allele, Map<Allele, Double>>> fullVQScoreMap = new HashMap<>();
        final Map<Long, Map<Allele, Map<Allele, String>>> fullYngMap = new HashMap<>();
        final Map<Long, List<String>> siteFilterMap = new HashMap<>();

        String rowRestriction = null;
        if (minLocation != null && maxLocation != null) {
            rowRestriction = "location >= " + minLocation + " AND location <= " + maxLocation;
        }
        final String rowRestrictionWithFilterSetName = rowRestriction + " AND " + SchemaUtils.FILTER_SET_NAME + " = '" + filterSetName + "'";

        boolean noVQScoreFilteringRequested = vqScoreFilteringType.equals(ExtractCohort.VQScoreFilteringType.NONE);
        if (!noVQScoreFilteringRequested) {
            // ensure vqScore filters (vqsLod or Sensitivity) are defined. this really shouldn't ever happen, said the engineer.
            if (vqScoreSNPThreshold == null || vqScoreINDELThreshold == null) {
                throw new UserException(getVQScoreFieldName() + " filtering thresholds for SNPs and INDELs must be defined.");
            }

            // get filter info (vqslod/sensitivity & yng values)
            try (StorageAPIAvroReader reader = new StorageAPIAvroReader(filterSetInfoTableRef, rowRestrictionWithFilterSetName, projectID)) {
                long filterRecordsProcessed = 0;
                long filterRecordsDropped = 0;
                for (final GenericRecord queryRow : reader) {
                    if (++filterRecordsProcessed % 100000 == 0) {
                        logger.info("Processed " + filterRecordsProcessed + " filter set info records, dropped " + filterRecordsDropped + ".");
                    }
                    final ExtractCohortFilterRecord filterRow = new ExtractCohortFilterRecord(queryRow, getVQScoreFieldName(), getScoreFieldName());

                    final long location = filterRow.getLocation();
                    final Allele ref = Allele.create(filterRow.getRefAllele(), true);
                    final Allele alt = Allele.create(filterRow.getAltAllele(), false);

                    if (!vbs.containsVariant(location, location + Math.max(ref.length(), alt.length()))) {
                        ++filterRecordsDropped;
                        continue;
                    }
                    final Double score = filterRow.getScore();
                    final Double vqsScore = filterRow.getVqScore();
                    final String yng = filterRow.getYng();

                    fullScoreMap.putIfAbsent(location, new HashMap<>());
                    fullScoreMap.get(location).putIfAbsent(ref, new HashMap<>());
                    fullScoreMap.get(location).get(ref).put(alt, score);
                    fullVQScoreMap.putIfAbsent(location, new HashMap<>());
                    fullVQScoreMap.get(location).putIfAbsent(ref, new HashMap<>());
                    fullVQScoreMap.get(location).get(ref).put(alt, vqsScore);
                    fullYngMap.putIfAbsent(location, new HashMap<>());
                    fullYngMap.get(location).putIfAbsent(ref, new HashMap<>());
                    fullYngMap.get(location).get(ref).put(alt, yng);
                }
                logger.info("Processed " + filterRecordsProcessed + " filter set info records, dropped " + filterRecordsDropped + ".");
                processBytesScanned(reader);
            }
        }

        // load site-level filter data into data structure
        if (filterSetSiteTableRef != null) {
            try (StorageAPIAvroReader reader = new StorageAPIAvroReader(filterSetSiteTableRef, rowRestrictionWithFilterSetName, projectID)) {
                long filterRecordsProcessed = 0;
                long filterRecordsDropped = 0;
                for (final GenericRecord queryRow : reader) {
                    if (++filterRecordsProcessed % 10000 == 0) {
                        logger.info("Processed " + filterRecordsProcessed + " filter set sites records, dropped " + filterRecordsDropped + ".");
                    }
                    long location = Long.parseLong(queryRow.get(SchemaUtils.LOCATION_FIELD_NAME).toString());
                    if (!vbs.containsVariant(location, location + 1)) {
                        ++filterRecordsDropped;
                        continue;
                    }
                    List<String> filters = Arrays.asList(queryRow.get(SchemaUtils.FILTERS).toString().split(","));
                    siteFilterMap.put(location, filters);
                }
                processBytesScanned(reader);
            }
        }

        createVariantsFromSortedRanges(sampleIdsToExtract, variantIterables, fullScoreMap, fullVQScoreMap, fullYngMap, samplePloidyMap, siteFilterMap, noVQScoreFilteringRequested);

        logger.info("Processed " + totalRangeRecords + " range records and rejected " + totalIrrelevantRangeRecords + " irrelevant ones ");
    }


    public static SortingCollection<GenericRecord> getAvroSortingCollection(org.apache.avro.Schema schema, int localSortMaxRecordsInRam) {
        final SortingCollection.Codec<GenericRecord> sortingCollectionCodec = new AvroSortingCollectionCodec(schema);
        final Comparator<GenericRecord> sortingCollectionComparator = (o1, o2) -> {
            final long firstPosition = (Long) o1.get(SchemaUtils.LOCATION_FIELD_NAME);
            final long secondPosition = (Long) o2.get(SchemaUtils.LOCATION_FIELD_NAME);

            final int result = Long.compare(firstPosition, secondPosition);
            if (result != 0) {
                return result;
            } else {
                final long firstSample = (Long) o1.get(SchemaUtils.SAMPLE_ID_FIELD_NAME);
                final long secondSample = (Long) o2.get(SchemaUtils.SAMPLE_ID_FIELD_NAME);
                return Long.compare(firstSample, secondSample);
            }
        };
        return SortingCollection.newInstance(GenericRecord.class, sortingCollectionCodec, sortingCollectionComparator, localSortMaxRecordsInRam, true);
    }

    private void addToVetSortingCollection(final SortingCollection<GenericRecord> sortingCollection,
                                           final Iterable<GenericRecord> avroReader,
                                           final VariantBitSet vbs) {
        int recordsProcessed = 0;
        long startTime = System.currentTimeMillis();

        // NOTE: if OverlapDetector takes too long, try using RegionChecker from tws_sv_local_assembler
        // need to manually add the upstream padding to the set of intervals
        List<SimpleInterval> overlapIntervals = new ArrayList<>(traversalIntervals);
        // Related to the problem in createSortedVetCollectionFromExtractTableBigQuery with having a minLocation within
        // IngestConstants.MAX_DELETION_SIZE of the start of a chromosome
        // These are 1 indexed, not 0 indexed.  1 is the first valid location.
        int adjustedStartingLocation = Math.max(1, (SchemaUtils.decodePosition(minLocation) - IngestConstants.MAX_DELETION_SIZE + 1));

        overlapIntervals.add(new SimpleInterval(SchemaUtils.decodeContig(minLocation),
                adjustedStartingLocation,
                SchemaUtils.decodePosition(minLocation)));
        final OverlapDetector<SimpleInterval> intervalsOverlapDetector = OverlapDetector.create(overlapIntervals);

        for (final GenericRecord queryRow : avroReader) {
            long location = (Long) queryRow.get(SchemaUtils.LOCATION_FIELD_NAME);
            int position = SchemaUtils.decodePosition(location);
            SimpleInterval simpleInterval = new SimpleInterval(SchemaUtils.decodeContig(location), position, position);

            if (intervalsOverlapDetector.overlapsAny(simpleInterval)) {
                vbs.setVariant(location);
                sortingCollection.add(queryRow);
                if (++recordsProcessed % 1000000 == 0) {
                    long endTime = System.currentTimeMillis();
                    logger.info("Processed " + recordsProcessed + " VET records in " + (endTime - startTime) + " ms");
                    startTime = endTime;
                }
            }
        }

        long endTime = System.currentTimeMillis();
        logger.info("Processed " + recordsProcessed + " VET records in " + (endTime - startTime) + " ms");

        sortingCollection.printTempFileStats();
    }

    private void addToRefSortingCollection(final SortingCollection<GenericRecord> sortingCollection, final Iterable<GenericRecord> avroReader, final VariantBitSet vbs) {
        int recordsProcessed = 0;
        long startTime = System.currentTimeMillis();

        for (final GenericRecord queryRow : avroReader) {
            long location = (Long) queryRow.get(SchemaUtils.LOCATION_FIELD_NAME);
            int length = ((Long) queryRow.get(SchemaUtils.LENGTH_FIELD_NAME)).intValue();

            if (vbs.containsVariant(location, location + length)) {
                sortingCollection.add(queryRow);
            }

            if (++recordsProcessed % 1000000 == 0) {
                long endTime = System.currentTimeMillis();
                logger.info("Processed " + recordsProcessed + " Reference Ranges records in " + (endTime - startTime) + " ms");
                startTime = endTime;
            }
        }

        long endTime = System.currentTimeMillis();
        logger.info("Processed " + recordsProcessed + " Reference Ranges records in " + (endTime - startTime) + " ms");

        sortingCollection.printTempFileStats();
    }

    private ExtractCohortRecord mergeSampleRecord(ExtractCohortRecord r1, ExtractCohortRecord r2) {
        if (printDebugInformation) {
            logger.info("In SampleRecord Merge Logic for " + r1 + " and " + r2);
        }

        final String r1State = r1.getState();
        final String r2State = r2.getState();

        // Valid states are m, 1-6 (ref), v, *
        // For now, just handle cases where '*' should be dropped in favor anything besides missing...
        if (r1State.equals("*") && !r2State.equals("m")) {
            return r2;
        } else if (r2State.equals("*") && !r1State.equals("m")) {
            return r1;
        } else {
            return r2;
        }
    }

    protected VariantContext processSampleRecordsForLocation(final long location,
                                                             final Iterable<ExtractCohortRecord> sampleRecordsAtPosition,
                                                             final Map<Long, Map<Allele, Map<Allele, Double>>> fullScoreMap,
                                                             final Map<Long, Map<Allele, Map<Allele, Double>>> fullVQScoreMap,
                                                             final Map<Long, Map<Allele, Map<Allele, String>>> fullYngMap,
                                                             final Map<String, Integer> samplePloidyMap,
                                                             final boolean noVQScoreFilteringRequested,
                                                             final Map<Long, List<String>> siteFilterMap,
                                                             final ExtractCohort.VQScoreFilteringType vqScoreFilteringType) {

        final List<VariantContext> unmergedCalls = new ArrayList<>();
        final List<ReferenceGenotypeInfo> refCalls = new ArrayList<>();

        int maxSampleIdToExtract = sampleIdsToExtract.last().intValue();
        final BitSet samplesSeen = new BitSet(maxSampleIdToExtract + 1);

        boolean currentPositionHasVariant = false;
        final int currentPosition = SchemaUtils.decodePosition(location);
        final String contig = SchemaUtils.decodeContig(location);
        final Allele refAllele = getReferenceAllele(refSource, location);
        int numRecordsAtPosition = 0;

        final Map<Allele, Map<Allele, Double>> scoreMap;
        final Map<Allele, Map<Allele, Double>> vqScoreMap;
        final Map<Allele, Map<Allele, String>> yngMap;

        // TODO: optimize in the case where noVQScoreFilteringRequested == true, no need to populate this

        // If there's no yng/score(vqslod/sensitivity) for this site, then we'll treat these as NAYs because VETS dropped them (they have no alt reads).
        if (fullVQScoreMap.get(SchemaUtils.encodeLocation(contig, currentPosition)) == null) {
            scoreMap = new HashMap<>();
            scoreMap.put(refAllele, new HashMap<>());
            vqScoreMap = new HashMap<>();
            vqScoreMap.put(refAllele, new HashMap<>());
            yngMap = new HashMap<>();
            yngMap.put(refAllele, new HashMap<>());
        } else {
            scoreMap = fullScoreMap.get(SchemaUtils.encodeLocation(contig, currentPosition));
            vqScoreMap = fullVQScoreMap.get(SchemaUtils.encodeLocation(contig, currentPosition));
            yngMap = fullYngMap.get(SchemaUtils.encodeLocation(contig, currentPosition));
        }

        for (final ExtractCohortRecord sampleRecord : sampleRecordsAtPosition) {
            final String sampleName = sampleIdToName.get(sampleRecord.getSampleId());

            if (sampleName == null) {
                throw new GATKException("Unable to translate sample id " + sampleRecord.getSampleId() + " to sample name");
            }

            // PERF: BOTTLENECK (~13%)
            // Note: we've already confirmed that max sample id is an int
            samplesSeen.set(sampleRecord.getSampleId().intValue());
            ++numRecordsAtPosition;

            if (printDebugInformation) {
                logger.info("\t" + contig + ":" + currentPosition + ": found record for sample " + sampleName + ": " + sampleRecord);
            }

            switch (sampleRecord.getState()) {
                case "v":   // Variant
                    ++totalNumberOfVariants;
                    VariantContext vc = createVariantContextFromSampleRecord(sampleRecord, vqScoreMap, yngMap, vqScoreFilteringType);
                    unmergedCalls.add(vc);

                    currentPositionHasVariant = true;
                    break;
                case "0":   // Non-Variant Block with GQ < 10
                    // Reference calls with GQ 0 should be rendered as no-call (#271)
                    // Nothing to do here -- just needed to mark the sample as seen so it doesn't get put in the high confidence ref band
                    break;
                case "1":  // Non-Variant Block with 10 <=  GQ < 20
                    refCalls.add(new ReferenceGenotypeInfo(sampleName, 10, sampleRecord.getSampleId().intValue()));
                    break;
                case "2":  // Non-Variant Block with 20 <= GQ < 30
                    refCalls.add(new ReferenceGenotypeInfo(sampleName, 20, sampleRecord.getSampleId().intValue()));
                    break;
                case "3":  // Non-Variant Block with 30 <= GQ < 40
                    refCalls.add(new ReferenceGenotypeInfo(sampleName, 30, sampleRecord.getSampleId().intValue()));
                    break;
                case "4":  // Non-Variant Block with 40 <= GQ < 50
                    refCalls.add(new ReferenceGenotypeInfo(sampleName, 40, sampleRecord.getSampleId().intValue()));
                    break;
                case "5":  // Non-Variant Block with 50 <= GQ < 60
                    refCalls.add(new ReferenceGenotypeInfo(sampleName, 50, sampleRecord.getSampleId().intValue()));
                    break;
                case "6":  // Non-Variant Block with 60 <= GQ (usually omitted from tables)
                    refCalls.add(new ReferenceGenotypeInfo(sampleName, 60, sampleRecord.getSampleId().intValue()));
                    break;
                case "*":   // Spanning Deletion - do nothing. just mark the sample as seen
                    break;
                case "m":   // Missing
                    // Nothing to do here -- just needed to mark the sample as seen, so that it doesn't get put in the high confidence ref band
                    break;
                default:
                    throw new GATKException("Unrecognized state: " + sampleRecord.getState());
            }

        }

        if (printDebugInformation) {
            logger.info(contig + ":" + currentPosition + ": processed " + numRecordsAtPosition + " total sample records");
        }

        return finalizeCurrentVariant(unmergedCalls, refCalls, samplesSeen, currentPositionHasVariant, location, contig, currentPosition, refAllele, scoreMap, vqScoreMap, yngMap, samplePloidyMap, noVQScoreFilteringRequested, siteFilterMap);
    }

    VariantContext finalizeCurrentVariant(final List<VariantContext> unmergedVariantCalls,
                                          final List<ReferenceGenotypeInfo> referenceCalls,
                                          final BitSet samplesSeen,
                                          final boolean currentPositionHasVariant,
                                          final long location,
                                          final String contig,
                                          final long start,
                                          final Allele refAllele,
                                          final Map<Allele, Map<Allele, Double>> scoreMap,
                                          final Map<Allele, Map<Allele, Double>> vqScoreMap,
                                          final Map<Allele, Map<Allele, String>> yngMap,
                                          final Map<String, Integer> samplePloidyMap,
                                          final boolean noVQScoreFilteringRequested,
                                          final Map<Long, List<String>> siteFilterMap) {
        // If there were no variants at this site, we don't emit a record and there's nothing to do here
        if (!currentPositionHasVariant) {
            return null;
        }

        VariantContext mergedVC = variantContextMerger.mergeWithRemapping(
                unmergedVariantCalls,
                new SimpleInterval(contig, (int) start, (int) start),
                refAllele.getBases()[0],
                true,
                false);


        // Reference Sites -- first create a single VC Builder
        final VariantContextBuilder vcWithRef = new VariantContextBuilder(mergedVC);

        // create alleles (same for all genotypes)
        List<Allele> gtAlleles = Arrays.asList(mergedVC.getReference(), mergedVC.getReference());
        List<Allele> gtHaploidAlleles = Arrays.asList(mergedVC.getReference());

        final GenotypesContext genotypes = GenotypesContext.copy(vcWithRef.getGenotypes());
        GenotypeBuilder genotypeBuilder = new GenotypeBuilder();

        // add known ref genotypes
        for (final ReferenceGenotypeInfo info : referenceCalls) {
            genotypeBuilder.reset(false);
            genotypeBuilder.name(info.getSampleName());

            long chromAsPosition = location - SchemaUtils.decodePosition(location);
            String key = makePloidyLookupKey(Long.toString(chromAsPosition), Integer.toString(info.getSampleId()));

            // Logic for determining the correct ploidy for reference data
            // If we have no info in the table, the ploidy is explicitly 2, OR we are in a PAR, use diploid reference.
            // If we have looked up the ploidy in our table and it says 1, use a haploid reference
            // Otherwise, if we have a ploidy that is neither 1 nor 2, throw a user exception because we haven't coded for this case
            int ploidy = (samplePloidyMap.containsKey(key) ? samplePloidyMap.get(key) : 2);
            if (ploidy == 2 || PloidyUtils.isLocationInPAR(location)) {
                genotypeBuilder.alleles(gtAlleles);
            } else if (ploidy == 1) {
                genotypeBuilder.alleles(gtHaploidAlleles);
            } else {
                throw new UserException("GVS cannot currently handle extracting with a ploidy of "+ploidy+" as seen at "+SchemaUtils.decodeContig(location)+": "+SchemaUtils.decodePosition(location)+".");
            }

            genotypeBuilder.GQ(info.getGQ());
            genotypes.add(genotypeBuilder.make());
        }

        // Find samples for inferred state and synthesize

        // mutates in place, so rename to be clear.
        // we don't use this again as "samplesSeen" so we don't need a full copy
        BitSet samplesNotEncountered = samplesSeen;
        samplesNotEncountered.xor(sampleIdsToExtractBitSet);

        // Iterate through the samples not encountered
        for (int sampleId = samplesNotEncountered.nextSetBit(0); sampleId >= 0; sampleId = samplesNotEncountered.nextSetBit(sampleId + 1)) {
            genotypeBuilder.reset(false);
            genotypeBuilder.name(sampleIdToName.get((long) sampleId));

            long chromAsPosition = location - SchemaUtils.decodePosition(location);
            String key = makePloidyLookupKey(Long.toString(chromAsPosition), Integer.toString(sampleId));

            // Logic for determining the correct ploidy for reference data
            // If we have no info in the table, the ploidy is explicitly 2, OR we are in a PAR, use diploid reference.
            // If we have looked up the ploidy in our table and it says 1, use a haploid reference
            // Otherwise, if we have a ploidy that is neither 1 nor 2, throw a user exception because we haven't coded for this case
            int ploidy = (samplePloidyMap.containsKey(key) ? samplePloidyMap.get(key) : 2);
            if (ploidy == 2 || PloidyUtils.isLocationInPAR(location)) {
                genotypeBuilder.alleles(gtAlleles);
            } else if (ploidy == 1) {
                genotypeBuilder.alleles(gtHaploidAlleles);
            } else {
                throw new UserException("GVS cannot currently handle extracting with a ploidy of "+ploidy+" as seen at "+SchemaUtils.decodeContig(location)+": "+SchemaUtils.decodePosition(location)+".");
            }
            genotypeBuilder.GQ(inferredReferenceState.getReferenceGQ());
            genotypes.add(genotypeBuilder.make());
        }

        vcWithRef.genotypes(genotypes);
        mergedVC = vcWithRef.make();

        if (maxAlternateAlleleCount.isPresent()) {
            if (mergedVC.getAlternateAlleles().size() > maxAlternateAlleleCount.getAsLong()) {
                    logger.warn("Dropping site on contig {} pos [{}, {}] with {} alt alleles versus maximum {}",
                            mergedVC.getContig(), mergedVC.getStart(), mergedVC.getEnd(),
                            mergedVC.getAlternateAlleles().size(), maxAlternateAlleleCount.getAsLong());
                return null;
            }
        }

        ReferenceContext referenceContext = new ReferenceContext(refSource, new SimpleInterval(mergedVC));

        VariantContext genotypedVC = annotationEngine.annotateContext(mergedVC, new FeatureContext(), referenceContext, null, a -> true);

        // apply VQScore (vqslod/calibration sensitivity)-based filters
        VariantContext filteredVC =
                noVQScoreFilteringRequested ? genotypedVC : filterSiteByAlleleSpecificVQScore(genotypedVC, scoreMap, vqScoreMap, yngMap, vqScoreFilteringType);

        // apply SiteQC-based filters, if they exist
        if (siteFilterMap.containsKey(location)) {
            final VariantContextBuilder sfBuilder = new VariantContextBuilder(filteredVC);

            Set<String> newFilters = new HashSet<>();

            // include existing filters, if any
            if (sfBuilder.getFilters() != null) {
                newFilters.addAll(sfBuilder.getFilters());
            }

            newFilters.addAll(siteFilterMap.get(location));
            filteredVC = sfBuilder.filters(newFilters).make();
        }

        // clean up extra annotations
        return removeAnnotations(filteredVC);
    }


    private String makePloidyLookupKey(final String chromosome, final String sampleId) {
        // This is almost stupidly simple, but we do this in more than one place and need to be consistent about it
        return chromosome+"-"+sampleId;
    }

    private VariantContext filterSiteByAlleleSpecificVQScore(VariantContext mergedVC,
                                                             Map<Allele, Map<Allele, Double>> scoreMap,
                                                             Map<Allele, Map<Allele, Double>> vqScoreMap,
                                                             Map<Allele, Map<Allele, String>> yngMap,
                                                             ExtractCohort.VQScoreFilteringType vqScoreFilteringType) {
        final Map<Allele, Double> remappedScoreMap = remapAllelesInMap(mergedVC, scoreMap, Double.NaN);
        final Map<Allele, Double> remappedVQScoreMap = remapAllelesInMap(mergedVC, vqScoreMap, Double.NaN);
        final Map<Allele, String> remappedYngMap = remapAllelesInMap(mergedVC, yngMap, VCFConstants.EMPTY_INFO_FIELD);

        final Map<Allele, Double> relevantScoreMap = new LinkedHashMap<>();
        mergedVC.getAlternateAlleles().forEach(key -> Optional.ofNullable(remappedScoreMap.get(key)).ifPresent(value -> relevantScoreMap.put(key, value)));
        final Map<Allele, Double> relevantVQScoreMap = new LinkedHashMap<>();
        mergedVC.getAlternateAlleles().forEach(key -> Optional.ofNullable(remappedVQScoreMap.get(key)).ifPresent(value -> relevantVQScoreMap.put(key, value)));

        final VariantContextBuilder builder = new VariantContextBuilder(mergedVC);

        if (getScoreFieldName() != null) {
            builder.attribute(getScoreKey(), relevantScoreMap.values().stream().map(val -> val.equals(Double.NaN) ? VCFConstants.EMPTY_INFO_FIELD : val.toString()).collect(Collectors.toList()));
        }
        builder.attribute(getAlleleSpecificVQSScoreKey(), relevantVQScoreMap.values().stream().map(val -> val.equals(Double.NaN) ? VCFConstants.EMPTY_INFO_FIELD : val.toString()).collect(Collectors.toList()));

        if (vqScoreFilteringType.equals(ExtractCohort.VQScoreFilteringType.SITES)) { // Note that these filters are not used with Genotype VQSLOD/Sensitivity Filtering
            int refLength = mergedVC.getReference().length();

            // if there are any Yays, the site is PASS
            if (remappedYngMap.containsValue("Y")) {
                builder.passFilters();
            } else if (remappedYngMap.containsValue("N")) {
                builder.filter(GATKVCFConstants.NAY_FROM_YNG);
            } else {
                // if it doesn't trigger any of the filters below, we assume it passes.
                builder.passFilters();
                //noinspection StatementWithEmptyBody
                if (remappedYngMap.containsValue("G")) {
                    if (isFailingSite(relevantVQScoreMap.entrySet().stream()
                            .filter(entry -> entry.getKey().length() == refLength)
                            .map(Map.Entry::getValue), vqScoreSNPThreshold)) {
                        builder.filter(getVqScoreSNPFailureFilterName());
                    }
                    if (isFailingSite(relevantVQScoreMap.entrySet().stream()
                            .filter(entry -> entry.getKey().length() != refLength)
                            .map(Map.Entry::getValue), vqScoreINDELThreshold)) {
                        builder.filter(getVqScoreINDELFailureFilterName());
                    }
                } else {
                    // per-conversation with Laura, if there is no information we let the site pass (ie no data does not imply failure)
                }
            }
        }

        // TODO: add in other annotations we need in output (like AF, etc?)
        return builder.make();
    }

    boolean isFailingSite(final Stream<Double> vqScores, final Double vqScoreThreshold) {
        Optional<Double> maxVal = vqScores
                .filter(d -> !(d.isNaN() || d.isInfinite()))
                .max(Double::compareTo);
        // It's a failing site if the maximum vqlod (if found) is less than the threshold
        return maxVal.isPresent() && maxVal.get() < vqScoreThreshold;
    }


    protected static VariantContext removeAnnotations(VariantContext vc) {

        final VariantContextBuilder builder = new VariantContextBuilder(vc);
        List<String> rmAnnotationList = new ArrayList<>(Arrays.asList(
                GATKVCFConstants.STRAND_ODDS_RATIO_KEY,
                GATKVCFConstants.AS_QUAL_BY_DEPTH_KEY,
                GATKVCFConstants.FISHER_STRAND_KEY));

        builder.rmAttributes(rmAnnotationList);
        return builder.make();
    }

    private <T> Map<Allele, T> remapAllelesInMap(VariantContext vc, Map<Allele, Map<Allele, T>> datamap, T emptyVal) {
        return remapAllelesInMap(vc.getReference(), vc.getAlternateAlleles(), vc.getContig(), vc.getStart(), datamap, emptyVal);
    }

    /*
     * Alleles from the filtering table need to be remapped to use the same ref allele that the will exist in the joined variant context.
     * This method changes the alleles in the datamap to match the representation that's in the vc.
     */
    private <T> Map<Allele, T> remapAllelesInMap(Allele ref, List<Allele> alternateAlleles, String contig, int start, Map<Allele, Map<Allele, T>> datamap, T emptyVal) {
        // create ordered results map
        Map<Allele, T> results = new LinkedHashMap<>();
        alternateAlleles.stream().forEachOrdered(allele -> results.put(allele, emptyVal));

        datamap.entrySet().stream().forEachOrdered(entry -> {
            if (entry.getKey().equals(ref)) {
                // reorder
                entry.getValue().entrySet().stream().forEach(altMapEntry -> results.put(altMapEntry.getKey(), altMapEntry.getValue()));
            } else {
                // remap
                int refLength = entry.getKey().length();
                List<Allele> allAlleles = new ArrayList<>();
                allAlleles.add(entry.getKey());
                allAlleles.addAll(entry.getValue().keySet());
                VariantContextBuilder vcb = new VariantContextBuilder("unused", contig, start, start + refLength - 1, allAlleles);
                VariantContext newvc = vcb.make();

                // If the length of the reference from the filtering table is longer than the reference in the variantContext, then that allele is not present in the extracted samples and we don't need the data
                if (refLength < ref.length()) {
                    Map<Allele, Allele> alleleMapping = GATKVariantContextUtils.createAlleleMapping(ref, newvc);
                    alleleMapping.entrySet().stream().forEach(mapped -> results.put(mapped.getValue(), entry.getValue().get(mapped.getKey())));
                }
            }
        });

        return results;
    }


    // vqScoreMap (vqsLod / calibration sensitivity), and yngMap are in/out parameters for this method. i.e. they are modified by this method
    private VariantContext createVariantContextFromSampleRecord(final ExtractCohortRecord sampleRecord,
                                                                Map<Allele, Map<Allele, Double>> vqScoreMap,
                                                                Map<Allele, Map<Allele, String>> yngMap,
                                                                ExtractCohort.VQScoreFilteringType vqScoreFilteringType) {
        final VariantContextBuilder builder = new VariantContextBuilder();
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder();

        final String contig = sampleRecord.getContig();
        final long startPosition = sampleRecord.getStart();
        final String sampleName = sampleIdToName.get(sampleRecord.getSampleId());

        if (sampleName == null) {
            throw new GATKException("Unable to translate sample id " + sampleRecord.getSampleId() + " to sample name");
        }

        builder.chr(contig);
        builder.start(startPosition);

        final List<Allele> alleles = new ArrayList<>();
        Allele ref = Allele.create(sampleRecord.getRefAllele(), true);
        alleles.add(ref);
        List<Allele> altAlleles = Arrays.stream(sampleRecord.getAltAllele().split(SchemaUtils.MULTIVALUE_FIELD_DELIMITER))
                .map(altAllele -> Allele.create(altAllele, false)).collect(Collectors.toList());

        // NOTE: gnarly needs this it seems? Should it be dropped now that we wont be using gnarly?
        altAlleles.add(Allele.NON_REF_ALLELE);

        alleles.addAll(altAlleles);
        builder.alleles(alleles);

        builder.stop(startPosition + alleles.get(0).length() - 1);

        genotypeBuilder.name(sampleName);
        vqScoreMap.putIfAbsent(ref, new HashMap<>());
        yngMap.putIfAbsent(ref, new HashMap<>());

        // need to re-prepend the leading "|" to AS_QUALapprox for use in gnarly
        if (sampleRecord.getAsQUALApprox() != null) {
            builder.attribute(SchemaUtils.AS_QUALapprox, "|" + sampleRecord.getAsQUALApprox());
        }

        if (sampleRecord.getQUALApprox() != null) {
            builder.attribute(SchemaUtils.QUALapprox, sampleRecord.getQUALApprox());
        }

        final String callGT = sampleRecord.getCallGT();
        genotypeBuilder.phased(callGT.contains("|"));

        String[] splitGT = callGT.split("[/|]");
        // This should match against anything like ".", "./.", ".|.", "././.", ".|.|.", etc
        if (callGT.matches("\\.([/|]\\.)*")) {
            genotypeBuilder.alleles(
                Arrays.stream(splitGT)
                    .map(alleleIndex -> Allele.NO_CALL)
                    .collect(Collectors.toList()));
        } else {
            List<Allele> genotypeAlleles =
                    Arrays.stream(splitGT)
                            .map(Integer::parseInt)
                            .map(alleles::get)
                            .collect(Collectors.toList());
            boolean noCallTheGenotype = false;

            if (vqScoreFilteringType.equals(ExtractCohort.VQScoreFilteringType.GENOTYPE)) {
                final List<Allele> nonRefAlleles =
                        genotypeAlleles.stream()
                                .filter(Allele::isNonReference)
                                .distinct()
                                .collect(Collectors.toList());

                final Map<Allele, Double> remappedVQScoreMap = remapAllelesInMap(ref, nonRefAlleles, contig, (int) startPosition, vqScoreMap, Double.NaN);
                final Map<Allele, String> remappedYngMap = remapAllelesInMap(ref, nonRefAlleles, contig, (int) startPosition, yngMap, VCFConstants.EMPTY_INFO_FIELD);

                // see https://github.com/broadinstitute/dsp-spec-ops/issues/291 for rationale
                // take "worst" outcome for yng/vqsScore (vqslod/sensitivity), evaluate each allele separately
                // if any allele is "N"ay, the genotype is filtered
                // if any allele is "Y"ay and the rest are "G"rey, the genotype is passed
                // if all alleles are "G"ray, the vqScore (VQSLod / calibration_sensitivity) is evaluated
                boolean anyNays = nonRefAlleles.stream().map(remappedYngMap::get).anyMatch("N"::equals);
                boolean anyYays = nonRefAlleles.stream().map(remappedYngMap::get).anyMatch("Y"::equals);

                // if there are any "N"s, the genotype is filtered
                if (anyNays) {
                    genotypeBuilder.filter(GATKVCFConstants.NAY_FROM_YNG);
                    if (convertFilteredGenotypesToNoCalls) {
                        noCallTheGenotype = true;
                    }
                } else //noinspection StatementWithEmptyBody
                    if (anyYays) {
                        // the genotype is passed, nothing to do here as non-filtered is the default
                    } else {
                        if (isFailingGenotype(nonRefAlleles.stream().filter(a -> a.length() == ref.length()),
                                remappedVQScoreMap, vqScoreSNPThreshold)) {
                            genotypeBuilder.filter(getVqScoreSNPFailureFilterName());
                            if (convertFilteredGenotypesToNoCalls) {
                                noCallTheGenotype = true;
                            }
                        }
                        if (isFailingGenotype(nonRefAlleles.stream().filter(a -> a.length() != ref.length()),
                                remappedVQScoreMap, vqScoreINDELThreshold)) {
                            genotypeBuilder.filter(getVqScoreINDELFailureFilterName());
                            if (convertFilteredGenotypesToNoCalls) {
                                noCallTheGenotype = true;
                            }
                        }
                    }
            }
            if (noCallTheGenotype) {
                // The genotype is to be no-called - that is, it would be FILTERED, but the option to no-call FT'd genotypes is enabled
                genotypeAlleles = Arrays.stream(splitGT)
                        .map(alleleIndex -> Allele.NO_CALL)
                        .collect(Collectors.toList());

            }
            genotypeBuilder.alleles(genotypeAlleles);
        }

        final String callGQ = sampleRecord.getCallGQ();
        if (callGQ != null) {
            genotypeBuilder.GQ(Integer.parseInt(callGQ));
        }

        final String callPGT = sampleRecord.getCallPGT();
        if (callPGT != null) {
            genotypeBuilder.attribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY, callPGT);
        }

        final String callPID = sampleRecord.getCallPID();
        if (callPID != null) {
            genotypeBuilder.attribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY, callPID);
        }

        final String callPS = sampleRecord.getCallPS();
        if (callPS != null) {
            genotypeBuilder.attribute(GATKVCFConstants.PHASE_SET_KEY, Integer.parseInt(callPS));
        }

        final String callPL = sampleRecord.getCallPL();
        if (this.emitPLs && callPL != null) {
            genotypeBuilder.PL(Arrays.stream(callPL.split(SchemaUtils.MULTIVALUE_FIELD_DELIMITER)).mapToInt(Integer::parseInt).toArray());
        }

        final String callAD = sampleRecord.getCallAD();
        if (this.emitADs && callAD != null) {
            genotypeBuilder.AD(Arrays.stream(callAD.split(SchemaUtils.MULTIVALUE_FIELD_DELIMITER)).mapToInt(Integer::parseInt).toArray());
        }

        final String callRGQ = sampleRecord.getCallRGQ();
        if (callRGQ != null) {
            genotypeBuilder.attribute(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY, callRGQ);
        }
        builder.genotypes(genotypeBuilder.make());

        return builder.make();
    }

    boolean isFailingGenotype(final Stream<Allele> nonRefAlleles,
                              final Map<Allele, Double> remappedVQScoreMap,
                              final Double vqScoreThreshold) {
        // get the max (best) vqslod of the alleles in this genotype
        Optional<Double> maxVal =
                nonRefAlleles
                        .map(remappedVQScoreMap::get)
                        .filter(Objects::nonNull)
                        .max(Double::compareTo);
        // It's a failing site if the maximum vqlod (if found) is less than the threshold
        return maxVal.isPresent() && maxVal.get() < vqScoreThreshold;
    }

    private SortingCollection<GenericRecord> createSortedVetCollectionFromBigQuery(final String projectID,
                                                                                   final String fqDatasetName,
                                                                                   final Set<Long> sampleIdsToExtract,
                                                                                   final Long minLocation,
                                                                                   final Long maxLocation,
                                                                                   final int localSortMaxRecordsInRam,
                                                                                   final VariantBitSet vbs
    ) {
        SortingCollection<GenericRecord> sortedVet = null;

        Map<Integer, LinkedList<Set<Long>>> tableMap = SampleList.mapSampleIdsToTableIndexes(sampleIdsToExtract);
        for (int tableIndex : tableMap.keySet()) {
            TableReference vetTableRef =
                    new TableReference(fqDatasetName + ".vet_" + String.format("%03d", tableIndex), this.extractVetFields);

            // TODO: Comment as to why (specific sample list clause)
            for (Set<Long> chunkSampleIds : tableMap.get(tableIndex)) {
                String sampleRestriction = " AND sample_id IN (" + StringUtils.join(chunkSampleIds, ",") + ")";

                // We need to look upstream MAX_DELETION_SIZE bases in case there is a deletion that begins before
                // the requested range, but spans into our processing range.  We don't use a "length" or end position
                // because it would break the clustering indexing

                // We want to look upstream... but don't want to go past the beginning of a chromosome.  Check for underflow by
                // calculating the start of the current chromosome.  Math.max ensures that it'll work as a floor, stopping
                // the subtraction from going too far back.
                long adjustedStartingLocation = (minLocation - IngestConstants.MAX_DELETION_SIZE + 1);
                long startOfChromosome = (minLocation / SchemaUtils.chromAdjustment) * SchemaUtils.chromAdjustment;
                adjustedStartingLocation = Math.max(startOfChromosome, adjustedStartingLocation);

                final String vetRowRestriction =
                        "location >= " + adjustedStartingLocation + " AND location <= " + maxLocation + sampleRestriction;
                try (StorageAPIAvroReader vetReader = new StorageAPIAvroReader(vetTableRef, vetRowRestriction, projectID)) {
                    if (sortedVet == null) {
                        sortedVet = getAvroSortingCollection(vetReader.getSchema(), localSortMaxRecordsInRam);
                    }

                    addToVetSortingCollection(sortedVet, vetReader, vbs);
                    processBytesScanned(vetReader);
                }
            }
        }
        return sortedVet;
    }

    private SortingCollection<GenericRecord> createSortedReferenceRangeCollectionFromBigQuery(final String projectID,
                                                                                              final String fqDatasetName,
                                                                                              final Set<Long> sampleIdsToExtract,
                                                                                              final Long minLocation,
                                                                                              final Long maxLocation,
                                                                                              final int localSortMaxRecordsInRam,
                                                                                              final VariantBitSet vbs
    ) {
        Map<Integer, LinkedList<Set<Long>>> tableMap = SampleList.mapSampleIdsToTableIndexes(sampleIdsToExtract);

        SortingCollection<GenericRecord> sortedReferenceRange = null;
        for (int tableIndex : tableMap.keySet()) {
            TableReference refTableRef =
                    new TableReference(fqDatasetName + ".ref_ranges_" + String.format("%03d", tableIndex), SchemaUtils.EXTRACT_REF_FIELDS);

            for (Set<Long> chunkSampleIds : tableMap.get(tableIndex)) {
                String sampleRestriction = " AND sample_id IN (" + StringUtils.join(chunkSampleIds, ",") + ")";

                // NOTE: MUST be written as location >= minLocation - MAX_REFERENCE_BLOCK_BASES +1 to not break cluster pruning in BigQuery
                final String refRowRestriction =
                        "location >= " + (minLocation - IngestConstants.MAX_REFERENCE_BLOCK_BASES + 1) + " AND location <= " + maxLocation + sampleRestriction;

                try (StorageAPIAvroReader refReader = new StorageAPIAvroReader(refTableRef, refRowRestriction, projectID)) {
                    if (sortedReferenceRange == null) {
                        sortedReferenceRange = getAvroSortingCollection(refReader.getSchema(), localSortMaxRecordsInRam);
                    }
                    addToRefSortingCollection(sortedReferenceRange, refReader, vbs);
                    processBytesScanned(refReader);
                }
            }
        }

        return sortedReferenceRange;
    }


    private VariantIterables createVariantIterablesFromUnsortedBigQueryRanges(
            final String fqDatasetName,
            final SortedSet<Long> sampleIdsToExtract,
            VariantBitSet vbs) {

        // We could handle this by making a map of BitSets or something, but it seems unnecessary to support this
        if (!SchemaUtils.decodeContig(minLocation).equals(SchemaUtils.decodeContig(maxLocation))) {
            throw new GATKException("Can not process cross-contig boundaries");
        }

        SortingCollection<GenericRecord> sortedVet = createSortedVetCollectionFromBigQuery(projectID,
                fqDatasetName,
                sampleIdsToExtract,
                minLocation,
                maxLocation,
                localSortMaxRecordsInRam,
                vbs);


        SortingCollection<GenericRecord> sortedReferenceRange = createSortedReferenceRangeCollectionFromBigQuery(projectID,
                fqDatasetName,
                sampleIdsToExtract,
                minLocation,
                maxLocation,
                localSortMaxRecordsInRam,
                vbs);

        return new VariantIterables(sortedVet, sortedReferenceRange);
    }

    //
    // BEGIN REF RANGES COHORT EXTRACT
    //
    private VariantIterables createVariantIterablesFromUnsortedExtractTableBigQueryRanges(
            final String fqVetTable,
            final String fqRefTable,
            VariantBitSet vbs) {

        // We could handle this by making a map of BitSets or something, but it seems unnecessary to support this
        if (!SchemaUtils.decodeContig(minLocation).equals(SchemaUtils.decodeContig(maxLocation))) {
            throw new GATKException("Can not process cross-contig boundaries");
        }

        SortingCollection<GenericRecord> sortedVet = createSortedVetCollectionFromExtractTableBigQuery(projectID,
                fqVetTable,
                minLocation,
                maxLocation,
                localSortMaxRecordsInRam,
                vbs);


        SortingCollection<GenericRecord> sortedReferenceRange = createSortedReferenceRangeCollectionFromExtractTableBigQuery(projectID,
                fqRefTable,
                minLocation,
                maxLocation,
                localSortMaxRecordsInRam,
                vbs);

        return new VariantIterables(sortedVet, sortedReferenceRange);
    }

    private SortingCollection<GenericRecord> createSortedVetCollectionFromExtractTableBigQuery(final String projectID,
                                                                                               final String fqVetTable,
                                                                                               final Long minLocation,
                                                                                               final Long maxLocation,
                                                                                               final int localSortMaxRecordsInRam,
                                                                                               final VariantBitSet vbs
    ) {

        TableReference tableRef = new TableReference(fqVetTable, this.extractVetFields);

        // We need to look upstream MAX_DELETION_SIZE bases in case there is a deletion that begins before
        // the requested range, but spans into our processing range.  We don't use a "length" or end position
        // because it would break the clustering indexing

        // We want to look upstream... but don't want to go past the beginning of a chromosome.  Check for underflow by
        // calculating the start of the current chromosome.  Math.max ensures that it'll work as a floor, stopping
        // the subtraction from going too far back.
        long adjustedStartingLocation = (minLocation - IngestConstants.MAX_DELETION_SIZE + 1);
        long startOfChromosome = (minLocation / SchemaUtils.chromAdjustment) * SchemaUtils.chromAdjustment;
        adjustedStartingLocation = Math.max(startOfChromosome, adjustedStartingLocation);

        final String vetRowRestriction =
                "location >= " + adjustedStartingLocation + " AND location <= " + maxLocation;
        try (StorageAPIAvroReader vetReader = new StorageAPIAvroReader(tableRef, vetRowRestriction, projectID)) {
            SortingCollection<GenericRecord> sortedVet = getAvroSortingCollection(vetReader.getSchema(), localSortMaxRecordsInRam);
            addToVetSortingCollection(sortedVet, vetReader, vbs);
            processBytesScanned(vetReader);
            return sortedVet;
        }
    }

    private SortingCollection<GenericRecord> createSortedReferenceRangeCollectionFromExtractTableBigQuery(final String projectID,
                                                                                                          final String fqRefTable,
                                                                                                          final Long minLocation,
                                                                                                          final Long maxLocation,
                                                                                                          final int localSortMaxRecordsInRam,
                                                                                                          final VariantBitSet vbs
    ) {

        TableReference tableRef = new TableReference(fqRefTable, SchemaUtils.EXTRACT_REF_FIELDS);

        // NOTE: MUST be written as location >= minLocation - MAX_REFERENCE_BLOCK_BASES +1 to not break cluster pruning in BigQuery
        final String refRowRestriction =
                "location >= " + (minLocation - IngestConstants.MAX_REFERENCE_BLOCK_BASES + 1) + " AND location <= " + maxLocation;

        try (StorageAPIAvroReader refReader = new StorageAPIAvroReader(tableRef, refRowRestriction, projectID)) {
            SortingCollection<GenericRecord> sortedReferenceRange = getAvroSortingCollection(refReader.getSchema(), localSortMaxRecordsInRam);
            addToRefSortingCollection(sortedReferenceRange, refReader, vbs);
            processBytesScanned(refReader);
            return sortedReferenceRange;
        }
    }

    //
    // END REF RANGES COHORT EXTRACT
    //
    private VariantIterables createVariantIterablesFromUnsortedAvroRanges(
            final GATKPath vetAvroFileName,
            final GATKPath refRangesAvroFileName,
            VariantBitSet vbs,
            final boolean presortedAvroFiles) {

        final AvroFileReader vetReader = new AvroFileReader(vetAvroFileName);
        final AvroFileReader refRangesReader = new AvroFileReader(refRangesAvroFileName);

        Iterable<GenericRecord> sortedVet;
        Iterable<GenericRecord> sortedReferenceRange;

        if (presortedAvroFiles) {
            sortedVet = vetReader;
            sortedReferenceRange = refRangesReader;
        } else {
            SortingCollection<GenericRecord> localSortedVet = getAvroSortingCollection(vetReader.getSchema(), localSortMaxRecordsInRam);
            addToVetSortingCollection(localSortedVet, vetReader, vbs);

            SortingCollection<GenericRecord> localSortedReferenceRange = getAvroSortingCollection(refRangesReader.getSchema(), localSortMaxRecordsInRam);
            addToRefSortingCollection(localSortedReferenceRange, refRangesReader, vbs);

            sortedVet = localSortedVet;
            sortedReferenceRange = localSortedReferenceRange;
        }

        return new VariantIterables(sortedVet, sortedReferenceRange);
    }

    void createVariantsFromSortedRanges(final SortedSet<Long> sampleIdsToExtract,
                                        final VariantIterables variantIterables,
                                        final Map<Long, Map<Allele, Map<Allele, Double>>> fullScoreMap,
                                        final Map<Long, Map<Allele, Map<Allele, Double>>> fullVQScoreMap,
                                        final Map<Long, Map<Allele, Map<Allele, String>>> fullYngMap,
                                        final Map<String, Integer> samplePloidyMap,
                                        final Map<Long, List<String>> siteFilterMap,
                                        final boolean noVQScoreFilteringRequested
    ) {

        long maxSampleId = sampleIdsToExtract.stream().max(Long::compare).orElseThrow(
                () -> new GATKException("Unable to calculate max sample id, sample list may be empty")
        );

        long minSampleId = sampleIdsToExtract.stream().min(Long::compare).orElseThrow(
                () -> new GATKException("Unable to calculate min sample id, sample list may be empty")
        );

        Long lastSample = null;
        Long lastPosition = null;

        final Map<Long, ExtractCohortRecord> currentPositionRecords = new HashMap<>(sampleIdToName.size() * 2);
        final Map<Long, Set<ReferenceRecord>> referenceCache = new HashMap<>(sampleIdToName.size());

        // Initialize cache
        for (Long sampleId : sampleIdsToExtract) {
            referenceCache.put(sampleId, new TreeSet<>());
        }

        Iterator<GenericRecord> sortedReferenceRangeIterator = variantIterables.refRanges.iterator();

        for (final GenericRecord sortedRow : variantIterables.vets) {
            final ExtractCohortRecord vetRow = new ExtractCohortRecord(sortedRow);
            long variantLocation = vetRow.getLocation();
            long variantSample = vetRow.getSampleId();

            // it's possible this variant is actually an upstream deletion before our region of interest,
            // so make sure this is not the case before we the actual variant record and discard
            if (variantLocation < minLocation) {
                handlePotentialSpanningDeletion(vetRow, referenceCache);
                continue;
            }

            // new position, fill in remainder of last position data
            // before continuing in to new position
            if (lastPosition != null && variantLocation != lastPosition) {
                // if the last VET was the last sample we're done
                if (lastSample != maxSampleId) {
                    processReferenceData(currentPositionRecords, sortedReferenceRangeIterator, referenceCache, lastPosition, lastSample + 1, maxSampleId, sampleIdsToExtract);
                }

                ++totalNumberOfSites;
                VariantContext vc = processSampleRecordsForLocation(lastPosition, currentPositionRecords.values(), fullScoreMap, fullVQScoreMap, fullYngMap, samplePloidyMap, noVQScoreFilteringRequested, siteFilterMap, vqScoreFilteringType);
                variantContextConsumer.accept(vc);
                currentPositionRecords.clear();

                lastSample = null;
            }

            long startingSample = (lastSample == null) ? minSampleId : lastSample + 1;
            processReferenceData(currentPositionRecords, sortedReferenceRangeIterator, referenceCache, variantLocation, startingSample, variantSample - 1, sampleIdsToExtract);

            // handle the actual variant record
            currentPositionRecords.merge(variantSample, vetRow, this::mergeSampleRecord);

            // if the variant record was a deletion, fabricate a spanning deletion row for the cache
            // TODO: is it possible that we will have two records now in the reference cache, and will we need logic to get "all" the records?
            // TODO: should we really build a VariantContext here, and let that sort out the length of the deletion for now, get the shortest of the alternates (biggest deletion)
            // TODO: use the genotypes of this specific sample (e.g. 0/1 vs 1/2) to decide how big the spanning deletion is.  The current logic matches what we do on ingest though
            // TODO: really, really, really think this through!!!
            handlePotentialSpanningDeletion(vetRow, referenceCache);

            lastPosition = variantLocation;
            lastSample = variantSample;
        }

        // finish writing out reference rows
        if ((lastSample != null) && (lastSample != maxSampleId)) {
            processReferenceData(currentPositionRecords, sortedReferenceRangeIterator, referenceCache, lastPosition, lastSample + 1, maxSampleId, sampleIdsToExtract);
        }

        if (!currentPositionRecords.isEmpty()) {
            ++totalNumberOfSites;
            VariantContext vc = processSampleRecordsForLocation(lastPosition, currentPositionRecords.values(), fullScoreMap, fullVQScoreMap, fullYngMap, samplePloidyMap, noVQScoreFilteringRequested, siteFilterMap, vqScoreFilteringType);
            variantContextConsumer.accept(vc);
        }
    }

    private void handlePotentialSpanningDeletion(ExtractCohortRecord vetRow, Map<Long, Set<ReferenceRecord>> referenceCache) {
        long position = vetRow.getLocation();
        long sample = vetRow.getSampleId();

        int smallestAltLength = Arrays.stream(vetRow.getAltAllele().split(",")).map(String::length).min(Integer::compare).orElse(0);
        int refLength = vetRow.getRefAllele().length();

        if (refLength > smallestAltLength) {
            referenceCache.get(sample).add(
                    new ReferenceRecord(position + 1, sample, refLength - smallestAltLength, "*")
            );
        }
    }

    private void processReferenceData(Map<Long, ExtractCohortRecord> currentPositionRecords, Iterator<GenericRecord> sortedReferenceRangeIterator, Map<Long, Set<ReferenceRecord>> referenceCache, long location, long fromSampleId, long toSampleId, SortedSet<Long> sampleIdsToExtract) {
        // in the case where there are two adjacent samples with variants, this method is called where from is greater than to
        // this is ok, there is just no reference data to process but subSet will throw an exception so we handle it with this if block
        if (toSampleId >= fromSampleId) {
            SortedSet<Long> samples = sampleIdsToExtract.subSet(fromSampleId, toSampleId + 1); // subset is start-inclusive, end-exclusive

            for (Long s : samples) {
                ExtractCohortRecord e = processReferenceData(sortedReferenceRangeIterator, referenceCache, location, s);
                currentPositionRecords.merge(s, e, this::mergeSampleRecord);
            }
        }
    }

    private ExtractCohortRecord processReferenceData(Iterator<GenericRecord> sortedReferenceRangeIterator, Map<Long, Set<ReferenceRecord>> referenceCache, long location, long sampleId) {
        ReferenceRecord r = processReferenceDataFromCache(referenceCache, location, sampleId);

        if (r == null) {
            r = processReferenceDataFromStream(sortedReferenceRangeIterator, referenceCache, location, sampleId);
        }

        return new ExtractCohortRecord(location, sampleId, r.getState());
    }

    // Refactoring opportunity:  Although this class is used as a singleton
    // the superclass ReferenceRecord requires a location, sample and length
    // as part of its construction.  This could possibly be more cleanly refactored
    // into a set of interfaces instead of an inheritance hierarchy
    public static class InferredReferenceRecord extends ReferenceRecord {
        public InferredReferenceRecord(GQStateEnum inferredState) {
            super(SchemaUtils.chromAdjustment + 1, -1, -1, inferredState.getValue());
        }
    }

    private ReferenceRecord processReferenceDataFromCache(Map<Long, Set<ReferenceRecord>> referenceCache, long location, long sampleId) {

        // PERF: bottleneck (iterator and remove) 6%
        Iterator<ReferenceRecord> iter = referenceCache.get(sampleId).iterator();
        while (iter.hasNext()) {
            ReferenceRecord row = iter.next();

            // if it ends before the current location, it's no longer relevant
            if (row.getEndLocation() < location) {
                iter.remove();
                continue;
            }

            // if it overlaps our location, use it!  Don't worry about removing it, if it
            // becomes irrelevant the next time we process it we will remove it
            if (row.getLocation() <= location && row.getEndLocation() >= location) {
                return row;
            }

            // completely after position, inferred state
            if (row.getLocation() > location) {
                return inferredReferenceRecord;
            }
        }

        // still here means QUEUE should be empty and we return a null state to indicate
        // no results were found.
        if (!referenceCache.get(sampleId).isEmpty()) {
            throw new GATKException("should have an empty queue: " + referenceCache.get(sampleId));
        } else {
            return null;
        }
    }

    private ReferenceRecord processReferenceDataFromStream(Iterator<GenericRecord> sortedReferenceRangeIterator, Map<Long, Set<ReferenceRecord>> referenceCache, long location, long sampleId) {
        while (sortedReferenceRangeIterator.hasNext()) {
            final ReferenceRecord refRow = new ReferenceRecord(sortedReferenceRangeIterator.next());
            totalRangeRecords++;

            // skip irrelevant entries
            if (refRow.getEndLocation() < location) {
                totalIrrelevantRangeRecords++;
                continue;
            }

            // is overlapping
            if (refRow.getLocation() <= location) {

                // always add to the appropriate cache for use downstream
                referenceCache.get(refRow.getSampleId()).add(refRow);

                // if this is for the requested sample, return the value
                if (refRow.getSampleId() == sampleId) {
                    return refRow;
                }
            }

            // we are now past the position, put this one entry in the cache and return the inferred state
            if (refRow.getLocation() > location) {
                referenceCache.get(refRow.getSampleId()).add(refRow);
                return inferredReferenceRecord;
            }
        }

        // if we are still here... use the inferred state
        return inferredReferenceRecord;
    }

    private static Allele getReferenceAllele(final ReferenceDataSource refSource, final long location) {
        String contig = SchemaUtils.decodeContig(location);
        long position = SchemaUtils.decodePosition(location);
        return Allele.create(refSource.queryAndPrefetch(contig, position, position).getBaseString(), true);
    }
}
