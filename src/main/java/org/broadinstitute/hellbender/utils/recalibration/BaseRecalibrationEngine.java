package org.broadinstitute.hellbender.utils.recalibration;

import htsjdk.samtools.CigarOperator;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.baq.BAQ;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.collections.NestedIntegerArray;
import org.broadinstitute.hellbender.utils.read.CigarBuilder;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.recalibration.covariates.Covariate;
import org.broadinstitute.hellbender.utils.recalibration.covariates.CovariateKeyCache;
import org.broadinstitute.hellbender.utils.recalibration.covariates.PerReadCovariateMatrix;
import org.broadinstitute.hellbender.utils.recalibration.covariates.StandardCovariateList;

import java.io.Serializable;
import java.util.Arrays;

public final class BaseRecalibrationEngine implements Serializable {
    private static final long serialVersionUID = 1L;

    protected static final Logger logger = LogManager.getLogger(BaseRecalibrationEngine.class);
    private final CovariateKeyCache keyCache;

    /*
     * Every call to EventType.values() (or any enum type) creates a new array instance but they are all equal (ie contain identical elements).
     * This is very expensive and wasteful when this array is created billions of times as in the case of BQSR.
     * The solution is to cache this array here.
     */
    private final EventType[] cachedEventTypes;

    /**
     * Reference window function for BQSR. For each read, returns an interval representing the span of
     * reference bases required by the BQSR algorithm for that read.
     *
     * Implemented as a static class rather than an anonymous class or lambda due to serialization issues in spark.
     */
    public static final class BQSRReferenceWindowFunction implements SerializableFunction<GATKRead, SimpleInterval> {
        private static final long serialVersionUID = 1L;

        @Override
        public SimpleInterval apply( GATKRead read ) {
            return BAQ.getReferenceWindowForRead(read, BAQ.DEFAULT_BANDWIDTH);
        }
    }
    public static final SerializableFunction<GATKRead, SimpleInterval> BQSR_REFERENCE_WINDOW_FUNCTION = new BQSRReferenceWindowFunction();

    private RecalibrationArgumentCollection recalArgs;

    private RecalibrationTables recalTables;

    private SAMFileHeader readsHeader;

    /**
     * list to hold the all the covariate objects that were requested (required + standard + experimental)
     */
    private StandardCovariateList covariates;

    private BAQ baq; // BAQ the reads on the fly to generate the alignment uncertainty vector
    private static final byte NO_BAQ_UNCERTAINTY = (byte)'@';

    private long numReadsProcessed = 0L;

    /**
     * Has finalizeData() been called?
     */
    private boolean finalized = false;

    public BaseRecalibrationEngine( final RecalibrationArgumentCollection recalArgs, final SAMFileHeader readsHeader ) {
        this.recalArgs = recalArgs;
        this.readsHeader = readsHeader;

        if (recalArgs.enableBAQ) {
            baq = new BAQ(recalArgs.BAQGOP); // setup the BAQ object with the provided gap open penalty
        } else {
            baq = null;
        }

        covariates = new StandardCovariateList(recalArgs, readsHeader);

        final int numReadGroups = readsHeader.getReadGroups().size();
        if ( numReadGroups < 1 ) {
            throw new UserException("Number of read groups must be >= 1, but is " + numReadGroups);
        }
        recalTables = new RecalibrationTables(covariates, numReadGroups);
        keyCache = new CovariateKeyCache();
        cachedEventTypes = recalArgs.computeIndelBQSRTables ? EventType.values() : new EventType[]{EventType.BASE_SUBSTITUTION};
    }

    public void logCovariatesUsed() {
        logger.info("The covariates being used here: ");
        for (final Covariate cov : covariates) { // list all the covariates being used
            logger.info('\t' + cov.getClass().getSimpleName());
        }
    }

    /**
     * For each read at this locus get the various covariate values and increment that location in the map based on
     * whether or not the base matches the reference at this particular location
     */
    public void processRead( final GATKRead originalRead, final ReferenceDataSource refDS, final Iterable<? extends Locatable> knownSites ) {
        final ReadTransformer transform = makeReadTransform();
        final GATKRead read = transform.apply(originalRead);

        if( read.isEmpty() ) {
            return; // the whole read was inside the adaptor so skip it
        }

        RecalUtils.updatePlatformForRead(read, readsHeader, recalArgs);

        int[] isSNP = new int[read.getLength()];
        int[] isInsertion = new int[isSNP.length];
        int[] isDeletion = new int[isSNP.length];

        //Note: this function modifies the isSNP, isInsertion and isDeletion arguments so it can't be skipped, BAQ or no BAQ
        final int nErrors = calculateIsSNPOrIndel(read, refDS, isSNP, isInsertion, isDeletion);

        // note for efficiency reasons we don't compute the BAQ array unless we actually have
        // some error to marginalize over.  For ILMN data ~85% of reads have no error
        final byte[] baqArray = (nErrors == 0 || !recalArgs.enableBAQ) ? flatBAQArray(read) : calculateBAQArray(read, refDS);
        // by default, baqArray is the array of constant value 64 i.e. [64, 64, .... , 64]
        if( baqArray != null ) { // some reads just can't be BAQ'ed
            final PerReadCovariateMatrix covariates = RecalUtils.computeCovariates(read, readsHeader, this.covariates, true, keyCache);
            final boolean[] skip = calculateSkipArray(read, knownSites); // skip known sites of variation as well as low quality and non-regular bases
            final double[] snpErrors = calculateFractionalErrorArray(isSNP, baqArray);
            final double[] insertionErrors = calculateFractionalErrorArray(isInsertion, baqArray);
            final double[] deletionErrors = calculateFractionalErrorArray(isDeletion, baqArray);

            // aggregate all of the info into our info object, and update the data
            final ReadRecalibrationInfo info = new ReadRecalibrationInfo(read, covariates, skip, snpErrors, insertionErrors, deletionErrors);
            updateRecalTablesForRead(info);
        }

        numReadsProcessed++;
    }

    /**
     * Finalize, if appropriate, all derived data in recalibrationTables.
     *
     * Called once after all calls to processRead have been issued.
     *
     * Assumes that all of the principal tables (by quality score) have been completely updated,
     * and walks over this data to create summary data tables like by read group table.
     */
    public void finalizeData() {
        Utils.validate(!finalized, "FinalizeData() has already been called");
        collapseQualityScoreTableToReadGroupTable(recalTables.getQualityScoreTable(), recalTables.getReadGroupTable());
        roundTableValues(recalTables);
        finalized = true;
    }


    /**
     * Populate the read group table, whose elements have been null up to this point,
     * by collapsing (marginalizing) the datum table index by both the read group and the reported quality score
     * over the reported quality score.
     *
     * Called once after all calls to updateDataForRead have been issued.
     *
     * @param byQualTable the RecalDatum table indexed by (read group, reported quality score)
     * @param byReadGroupTable the empty RecalDatum table to be populated by this method.
     *
     * TODO: this method should take a qual table and return a read group table.
     */
    public static void collapseQualityScoreTableToReadGroupTable(final NestedIntegerArray<RecalDatum> byQualTable,
                                                                 final NestedIntegerArray<RecalDatum> byReadGroupTable) {
        // the read group table has shape: (num read groups) x (num error modes)
        // the qual table has shape:       (num read groups) x (num reported qualities [default = MAX_PHRED_SCORE + 1 = 94]) x (num error modes)

        // iterate over all values in the qual table
        final int readGroupIndex = 0;
        final int errorModeIndex = 2;

        for ( final NestedIntegerArray.Leaf<RecalDatum> leaf : byQualTable.getAllLeaves() ) {
            final int rgKey = leaf.keys[readGroupIndex];
            final int eventIndex = leaf.keys[errorModeIndex];
            final RecalDatum rgDatum = byReadGroupTable.get(rgKey, eventIndex);
            final RecalDatum qualDatum = leaf.value;

            if ( rgDatum == null ) {
                // create a copy of qualDatum, and initialize byReadGroup table with it
                byReadGroupTable.put(new RecalDatum(qualDatum), rgKey, eventIndex);
            } else {
                // combine the qual datum with the existing datum in the byReadGroup table
                rgDatum.combine(qualDatum);
            }
        }
    }

    /**
     * To replicate the results of BQSR whether or not we save tables to disk (which we need in Spark),
     * we need to trim the numbers to a few decimal placed (that's what writing and reading does).
     */
    public static void roundTableValues(final RecalibrationTables rt) {
        for (int i = 0; i < rt.numTables(); i++) {
            for (final NestedIntegerArray.Leaf<RecalDatum> leaf : rt.getTable(i).getAllLeaves()) {
                // Empirical quality is implemented an integer qual score and does not need rounding.
                leaf.value.setNumMismatches(MathUtils.roundToNDecimalPlaces(leaf.value.getNumMismatches(), RecalUtils.NUMBER_ERRORS_DECIMAL_PLACES));
                leaf.value.setReportedQuality(MathUtils.roundToNDecimalPlaces(leaf.value.getReportedQuality(), RecalUtils.REPORTED_QUALITY_DECIMAL_PLACES));
            }
        }
    }

    /**
     * Get a possibly not-final recalibration table, to deal with distributed execution.
     */
    public RecalibrationTables getRecalibrationTables() {
        return recalTables;
    }

    /**
     * Get the final recalibration tables, after finalizeData() has been called
     *
     * This returns the finalized recalibration table collected by this engine.
     *
     * It is an error to call this function before finalizeData has been called
     *
     * @return the finalized recalibration table collected by this engine
     */
    public RecalibrationTables getFinalRecalibrationTables() {
        Utils.validate(finalized, "Cannot get final recalibration tables until finalizeData() has been called");
        return recalTables;
    }

    public StandardCovariateList getCovariates() {
        return covariates;
    }

    public long getNumReadsProcessed() {
        return numReadsProcessed;
    }

    /**
     * Update the recalibration statistics using the information in recalInfo.
     *
     * Implementation detail: we only populate the quality score table and the optional tables.
     * The read group table will be populated later by collapsing the quality score table.
     *
     * @param recalInfo data structure holding information about the recalibration values for a single read
     */
    private void updateRecalTablesForRead( final ReadRecalibrationInfo recalInfo ) {
        Utils.validate(!finalized, "FinalizeData() has already been called");

        final GATKRead read = recalInfo.getRead();
        final PerReadCovariateMatrix perReadCovariateMatrix = recalInfo.getCovariatesValues();
        final NestedIntegerArray<RecalDatum> qualityScoreTable = recalTables.getQualityScoreTable();

        final int nCovariates = covariates.size();
        final int readLength = read.getLength();
        for( int offset = 0; offset < readLength; offset++ ) {
            if( ! recalInfo.skip(offset) ) {
                for (int idx = 0; idx < cachedEventTypes.length; idx++) { //Note: we loop explicitly over cached values for speed
                    final EventType eventType = cachedEventTypes[idx];
                    final int[] covariatesAtOffset = perReadCovariateMatrix.getCovariatesAtOffset(offset, eventType);
                    final int eventIndex = eventType.ordinal();
                    final byte qual = recalInfo.getQual(eventType, offset);
                    final double isError = recalInfo.getErrorFraction(eventType, offset);

                    final int readGroup = covariatesAtOffset[StandardCovariateList.READ_GROUP_COVARIATE_DEFAULT_INDEX];
                    final int baseQuality = covariatesAtOffset[StandardCovariateList.BASE_QUALITY_COVARIATE_DEFAULT_INDEX];

                    RecalUtils.incrementDatum3keys(qualityScoreTable, qual, isError, readGroup, baseQuality, eventIndex);

                    for (int i = RecalUtils.NUM_REQUIRED_COVARIATES; i < nCovariates; i++) {
                        final int specialCovariate = covariatesAtOffset[i];
                        if (specialCovariate >= 0) {
                            RecalUtils.incrementDatum4keys(recalTables.getTable(i), qual, isError,
                                    readGroup, baseQuality, specialCovariate, eventIndex);
                        }
                    }
                }
            }
        }
    }

    private ReadTransformer makeReadTransform() {
        ReadTransformer f0 = BaseRecalibrationEngine::consolidateCigar;

        ReadTransformer f = f0.andThen(this::setDefaultBaseQualities)
                .andThen(this::resetOriginalBaseQualities)
                .andThen(ReadClipper::hardClipAdaptorSequence)
                .andThen(ReadClipper::hardClipSoftClippedBases);

        return f;
    }

    private static GATKRead consolidateCigar( final GATKRead read ) {
        // Always consolidate the cigar string into canonical form, collapsing zero-length / repeated cigar elements.
        // Downstream code cannot necessarily handle non-consolidated cigar strings.
        read.setCigar(new CigarBuilder().addAll(read.getCigar()).make());
        return read;
    }

    private GATKRead resetOriginalBaseQualities( final GATKRead read ) {
        if (! recalArgs.useOriginalBaseQualities) {
            return read;
        }
        return ReadUtils.resetOriginalBaseQualities(read);
    }

    private GATKRead setDefaultBaseQualities( final GATKRead read ) {
        // if we are using default quals, check if we need them, and add if necessary.
        // 1. we need if reads are lacking or have incomplete quality scores
        // 2. we add if defaultBaseQualities has a positive value
        if (recalArgs.defaultBaseQualities < 0) {
            return read;
        }
        byte[] reads = read.getBases();
        byte[] quals = read.getBaseQualities();
        if (quals == null || quals.length < reads.length) {
            byte[] new_quals = new byte[reads.length];
            Arrays.fill(new_quals, recalArgs.defaultBaseQualities);
            read.setBaseQualities(new_quals);
        }
        return read;
    }

    /**
     * Outputs a boolean array that has the same length as the read.
     * The array contains true at each index if the position meets one of the following criteria:
     *   1) not a regular base
     *   2) base quality is less than 6
     *   3) is a known site.
     */
    private boolean[] calculateSkipArray( final GATKRead read, final Iterable<? extends Locatable> knownSites ) {
        final int readLength = read.getLength();
        final boolean[] skip = new boolean[readLength];
        final boolean[] knownSitesArray = calculateKnownSites(read, knownSites);
        for(int i = 0; i < readLength; i++ ) {
            skip[i] = !BaseUtils.isRegularBase(read.getBase(i)) || read.getBaseQuality(i) < recalArgs.PRESERVE_QSCORES_LESS_THAN || knownSitesArray[i];
        }
        return skip;
    }

    /**
     * Outputs a boolean array that has the same length as the read and contains true at positions where known events
     * occur, as determined by the knownSites variable.
     */
    private static boolean[] calculateKnownSites( final GATKRead read, final Iterable<? extends Locatable> knownSites ) {
        final int readLength = read.getLength();
        final boolean[] knownSitesArray = new boolean[readLength];//initializes to all false
        final int softStart = read.getSoftStart();
        final int softEnd = read.getSoftEnd();
        for ( final Locatable knownSite : knownSites ) {
            if (knownSite.getEnd() < softStart || knownSite.getStart() > softEnd) {
                // knownSite is outside clipping window for the read, ignore
                continue;
            }
            final Pair<Integer, CigarOperator> featureStartAndOperatorOnRead = ReadUtils.getReadIndexForReferenceCoordinate(read, knownSite.getStart());
            int featureStartOnRead = featureStartAndOperatorOnRead.getLeft() == ReadUtils.READ_INDEX_NOT_FOUND ? 0 : featureStartAndOperatorOnRead.getLeft();
            if (featureStartAndOperatorOnRead.getRight() == CigarOperator.DELETION) {
                featureStartOnRead--;
            }

            final Pair<Integer, CigarOperator> featureEndAndOperatorOnRead = ReadUtils.getReadIndexForReferenceCoordinate(read, knownSite.getEnd());
            int featureEndOnRead = featureEndAndOperatorOnRead.getLeft() == ReadUtils.READ_INDEX_NOT_FOUND ? readLength : featureEndAndOperatorOnRead.getLeft();

            if( featureStartOnRead > readLength ) {
                featureStartOnRead = featureEndOnRead = readLength;
            }

            Arrays.fill(knownSitesArray, Math.max(0, featureStartOnRead), Math.min(readLength, featureEndOnRead + 1), true);
        }
        return knownSitesArray;
    }

    /**
     * Locates all SNP and indel events, storing them in the provided snp, isIns, and isDel arrays, and returns
     * the total number of SNP/indel events.
     *
     * @param read read to inspect
     * @param ref source of reference bses
     * @param snp storage for snp events (must be of length read.getBases().length and initialized to all 0's)
     * @param isIns storage for insertion events (must be of length read.getBases().length and initialized to all 0's)
     * @param isDel storage for deletion events (must be of length read.getBases().length and initialized to all 0's)
     * @return the total number of SNP and indel events
     */
    protected static int calculateIsSNPOrIndel(final GATKRead read, final ReferenceDataSource ref, int[] snp, int[] isIns, int[] isDel) {
        final byte[] refBases = ref.queryAndPrefetch(read.getContig(), read.getStart(), read.getEnd()).getBases();
        int readPos = 0;
        int refPos = 0;
        int nEvents = 0;

        for (final CigarElement ce : read.getCigarElements()) {
            final int elementLength = ce.getLength();
            switch (ce.getOperator()) {
                case M:
                case EQ:
                case X:
                    for (int i = 0; i < elementLength; i++) {
                        int snpInt = (BaseUtils.basesAreEqual(read.getBase(readPos), refBases[refPos]) ? 0 : 1);
                        snp[readPos] = snpInt;
                        nEvents += snpInt;
                        readPos++;
                        refPos++;
                    }
                    break;
                case D: {
                    final int index = (read.isReverseStrand() ? readPos : readPos - 1);
                    updateIndel(isDel, index);
                    refPos += elementLength;
                    break;
                }
                case N:
                    refPos += elementLength;
                    break;
                case I: {
                    final boolean forwardStrandRead = !read.isReverseStrand();
                    if (forwardStrandRead) {
                        updateIndel(isIns, readPos - 1);
                    }
                    readPos += elementLength;
                    if (!forwardStrandRead) {
                        updateIndel(isIns, readPos);
                    }
                    break;
                }
                case S: // ReferenceContext doesn't have the soft clipped bases!
                    readPos += elementLength;
                    break;
                case H:
                case P:
                    break;
                default:
                    throw new GATKException("Unsupported cigar operator: " + ce.getOperator());
            }
        }
        // we don't sum those as we go because they might set the same place to 1 twice
        nEvents += MathUtils.sum(isDel) + MathUtils.sum(isIns);
        return nEvents;
    }

    private static void updateIndel(final int[] indel, final int index) {
        if (index >= 0 && index < indel.length) {
            // protect ourselves from events at the start or end of the read (1D3M or 3M1D)
            indel[index] = 1;
        }
    }

    public static double[] calculateFractionalErrorArray( final int[] errorArray, final byte[] baqArray ) {
        if ( errorArray.length != baqArray.length ) {
            throw new GATKException("Array length mismatch detected. Malformed read?");
        }

        final int BLOCK_START_UNSET = -1;

        final double[] fractionalErrors = new double[baqArray.length];
        boolean inBlock = false;
        int blockStartIndex = BLOCK_START_UNSET;
        int i;
        for( i = 0; i < fractionalErrors.length; i++ ) {
            if( baqArray[i] == NO_BAQ_UNCERTAINTY ) {
                if( !inBlock ) {
                    fractionalErrors[i] = errorArray[i];
                } else {
                    calculateAndStoreErrorsInBlock(i, blockStartIndex, errorArray, fractionalErrors);
                    inBlock = false; // reset state variables
                    blockStartIndex = BLOCK_START_UNSET; // reset state variables
                }
            } else {
                inBlock = true;
                if( blockStartIndex == BLOCK_START_UNSET ) { blockStartIndex = i; }
            }
        }
        if( inBlock ) {
            calculateAndStoreErrorsInBlock(i-1, blockStartIndex, errorArray, fractionalErrors);
        }
        if( fractionalErrors.length != errorArray.length ) {
            throw new GATKException("Output array length mismatch detected. Malformed read?");
        }
        return fractionalErrors;
    }

    private static void calculateAndStoreErrorsInBlock( final int i,
                                                        final int blockStartIndex,
                                                        final int[] errorArray,
                                                        final double[] fractionalErrors ) {
        int totalErrors = 0;
        for( int j = Math.max(0, blockStartIndex - 1); j <= i; j++ ) {
            totalErrors += errorArray[j];
        }
        for( int j = Math.max(0, blockStartIndex - 1); j <= i; j++ ) {
            fractionalErrors[j] = ((double) totalErrors) / ((double)(i - Math.max(0, blockStartIndex - 1) + 1));
        }
    }

    /**
     * Create a BAQ style array that indicates no alignment uncertainty
     * @param read the read for which we want a BAQ array
     * @return a BAQ-style non-null byte[] counting NO_BAQ_UNCERTAINTY values
     * // TODO -- could be optimized avoiding this function entirely by using this inline if the calculation code above
     */
    protected static byte[] flatBAQArray(final GATKRead read) {
        final byte[] baq = new byte[read.getLength()];
        Arrays.fill(baq, NO_BAQ_UNCERTAINTY);
        return baq;
    }

    /**
     * Compute an actual BAQ array for read, based on its quals and the reference sequence
     * @param read the read to BAQ
     * @return a non-null BAQ tag array for read
     */
    private byte[] calculateBAQArray( final GATKRead read, final ReferenceDataSource refDS ) {
        baq.baqRead(read, refDS, BAQ.CalculationMode.RECALCULATE, BAQ.QualityMode.ADD_TAG);
        return BAQ.getBAQTag(read);
    }
}