package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeArguments;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParams;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.haplotype.PartiallyDeterminedHaplotype;
import org.broadinstitute.hellbender.utils.pairhmm.*;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.OutputStreamWriter;
import java.util.*;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;

/**
 * Classic likelihood computation: full pair-hmm all haplotypes vs all reads.
 *
 * Note this is the "Partially Determined" Pair hmm engine, which means the haplotypes it sees are expected to be
 * PartiallyDeterminedHaplotype objects. These are special in that they all have a subset of "determined" or present
 * bases corresponding to variants and everything else is "undetermined." The likelihood scores are computed by taking
 * the maximum possible score for the read vs the un-determined bases and behaving normally over the determined ones
 * that have a real subset of alleles from the assembly region.
 *
 * Since the operations and methods differ from the normal PairHMM significantly we have opted to create an entirely seperate
 * codepath of PD___HMM classes to handle the differences.
 */
public final class PDPairHMMLikelihoodCalculationEngine implements ReadLikelihoodCalculationEngine {


    private static final int MIN_ADJUSTED_QSCORE = 10;

    @VisibleForTesting
    static final double INITIAL_QSCORE = 40.0;
    public static final String HMM_BASE_QUALITIES_TAG = "HMMQuals";

    private final byte constantGCP;

    private final double log10globalReadMismappingRate;

    private final PDPairHMM pdPairHMM;
    public static final String UNCLIPPED_ORIGINAL_SPAN_ATTR = "originalAlignment";

    // DRAGEN-GATK related parameters
    private final DragstrParams dragstrParams;
    private final boolean dynamicDisqualification;
    private final double readDisqualificationScale;
    private final double expectedErrorRatePerBase;
    private final boolean disableCapReadQualitiesToMapQ;
    private final boolean symmetricallyNormalizeAllelesToReference;
    private final boolean modifySoftclippedBases;
    private final int rangeForReadOverlapToDeterminedBases;

    private final PairHMMLikelihoodCalculationEngine.PCRErrorModel pcrErrorModel;

    private final byte baseQualityScoreThreshold;

    /**
     * Create a new PairHMMLikelihoodCalculationEngine using provided parameters and hmm to do its calculations
     *
     * @param constantGCP the gap continuation penalty to use with the PairHMM
     * @param hmmType the type of the HMM to use
     * @param resultsFile output file to dump per-read, per-haplotype inputs and outputs for debugging purposes (null if not enabled).
     * @param log10globalReadMismappingRate the global mismapping probability, in log10(prob) units.  A value of
     *                                      -3 means that the chance that a read doesn't actually belong at this
     *                                      location in the genome is 1 in 1000.  The effect of this parameter is
     *                                      to cap the maximum likelihood difference between the reference haplotype
     *                                      and the best alternative haplotype by -3 log units.  So if the best
     *                                      haplotype is at -10 and this parameter has a value of -3 then even if the
     *                                      reference haplotype gets a score of -100 from the pairhmm it will be
     *                                      assigned a likelihood of -13.
     * @param pcrErrorModel model to correct for PCR indel artifacts
     * @param baseQualityScoreThreshold Base qualities below this threshold will be reduced to the minimum usable base
     *                                  quality.
     */
    public PDPairHMMLikelihoodCalculationEngine(final byte constantGCP,
                                                final DragstrParams dragstrParams,
                                                final PairHMMNativeArguments arguments,
                                                final PairHMM.Implementation hmmType,
                                                final GATKPath resultsFile,
                                                final double log10globalReadMismappingRate,
                                                final PairHMMLikelihoodCalculationEngine.PCRErrorModel pcrErrorModel,
                                                final byte baseQualityScoreThreshold,
                                                final boolean dynamicReadDisqualificaiton,
                                                final double readDisqualificationScale,
                                                final double expectedErrorRatePerBase,
                                                final boolean symmetricallyNormalizeAllelesToReference,
                                                final boolean disableCapReadQualitiesToMapQ,
                                                final boolean modifySoftclippedBases,
                                                final int rangeForReadOverlapToDeterminedBases) {
        Utils.nonNull(hmmType, "hmmType is null");
        Utils.nonNull(pcrErrorModel, "pcrErrorModel is null");
        if (constantGCP < 0){
            throw new IllegalArgumentException("gap continuation penalty must be non-negative");
        }
        if (log10globalReadMismappingRate > 0){
            throw new IllegalArgumentException("log10globalReadMismappingRate must be negative");
        }
        this.dragstrParams = dragstrParams;
        this.constantGCP = constantGCP;
        this.log10globalReadMismappingRate = log10globalReadMismappingRate;
        this.pcrErrorModel = this.dragstrParams == null ? pcrErrorModel : PairHMMLikelihoodCalculationEngine.PCRErrorModel.NONE;
        // TODO later we probably need a LOG and LOGLESS version for parsimony with DRAGEN
        this.pdPairHMM = new VectorLoglessPairPDHMM(VectorLoglessPairPDHMM.Implementation.OMP , new PairHMMNativeArguments());
        if (resultsFile != null) {
            pdPairHMM.setAndInitializeDebugOutputStream(new OutputStreamWriter(resultsFile.getOutputStream()));
        }
        this.dynamicDisqualification = dynamicReadDisqualificaiton;
        this.readDisqualificationScale = readDisqualificationScale;
        this.symmetricallyNormalizeAllelesToReference = symmetricallyNormalizeAllelesToReference;
        this.expectedErrorRatePerBase = expectedErrorRatePerBase;
        this.disableCapReadQualitiesToMapQ = disableCapReadQualitiesToMapQ;
        this.modifySoftclippedBases = modifySoftclippedBases;
        this.rangeForReadOverlapToDeterminedBases = rangeForReadOverlapToDeterminedBases;

        initializePCRErrorModel();

        if (baseQualityScoreThreshold < QualityUtils.MIN_USABLE_Q_SCORE) {
            throw new IllegalArgumentException("baseQualityScoreThreshold must be greater than or equal to " + QualityUtils.MIN_USABLE_Q_SCORE + " (QualityUtils.MIN_USABLE_Q_SCORE)");
        }
        this.baseQualityScoreThreshold = baseQualityScoreThreshold;
    }


    @Override
    public void close() {
        pdPairHMM.close();
    }

    @Override
    public void filterPoorlyModeledEvidence(AlleleLikelihoods<GATKRead, ? extends Haplotype> result, boolean dynamicDisqualification, double expectedErrorRatePerBase, double readDisqualificationScale) {
        ReadLikelihoodCalculationEngine.super.filterPoorlyModeledEvidence(result, dynamicDisqualification, expectedErrorRatePerBase, readDisqualificationScale);
    }

    @Override
    public ToDoubleFunction<GATKRead> log10MinTrueLikelihood(double maximumErrorPerBase, boolean capLikelihoods) {
        return ReadLikelihoodCalculationEngine.super.log10MinTrueLikelihood(maximumErrorPerBase, capLikelihoods);
    }

    @Override
    public AlleleLikelihoods<GATKRead, Haplotype> computeReadLikelihoods(List<Haplotype> haplotypeList, SAMFileHeader hdr, SampleList samples, Map<String, List<GATKRead>> perSampleReadList, boolean filterPoorly) {
        throw new GATKException.ShouldNeverReachHereException("We should never get here, this HMM engine exists only for PD haplotype computation");
    }


    @Override
    @SuppressWarnings("unchecked")
    // NOTE that the PairHMM doesn't need to interrogate the header so we skip checking it for this version
    public AlleleLikelihoods<GATKRead, Haplotype> computeReadLikelihoods(AssemblyResultSet assemblyResultSet, SampleList samples,
                                                                         Map<String, List<GATKRead>> perSampleReadList,
                                                                         boolean filterPoorly) {
        Utils.nonNull(assemblyResultSet, "assemblyResultSet is null");
        if (assemblyResultSet.getHaplotypeList().stream().anyMatch(haplotype -> !haplotype.isPartiallyDetermined())) {
            throw new GATKException("PDHMM engine requires PartiallyDeterminedHaplotype objects as input");
        }
        final List<PartiallyDeterminedHaplotype> haplotypeList = assemblyResultSet.getHaplotypeList().stream().map(h -> (PartiallyDeterminedHaplotype)h).collect(Collectors.toList());

        AlleleLikelihoods<GATKRead, Haplotype> resut = (AlleleLikelihoods<GATKRead, Haplotype>) computeReadLikelihoodsPartiallyDetermined(haplotypeList, null, samples, perSampleReadList, filterPoorly);
        return resut;
    }

    public AlleleLikelihoods<GATKRead, ? extends Haplotype> computeReadLikelihoodsPartiallyDetermined(final List<PartiallyDeterminedHaplotype> haplotypeList,
                                                                                                      final SAMFileHeader hdr,
                                                                                                      final SampleList samples,
                                                                                                      final Map<String, List<GATKRead>> perSampleReadList, final boolean filterPoorly) {
        Utils.nonNull(samples, "samples is null");
        Utils.nonNull(perSampleReadList, "perSampleReadList is null");
        Utils.nonNull(haplotypeList, "haplotypeList is null");

        final AlleleList<PartiallyDeterminedHaplotype> haplotypes = new IndexedAlleleList<>(haplotypeList);

        initializePairHMM(haplotypeList, perSampleReadList);

        // Add likelihoods for each sample's reads to our result
        final AlleleLikelihoods<GATKRead, PartiallyDeterminedHaplotype> result = new AlleleLikelihoods<>(samples, haplotypes, perSampleReadList);
        final int sampleCount = result.numberOfSamples();
        for (int i = 0; i < sampleCount; i++) {
            computeReadLikelihoods(result.sampleMatrix(i));
        }

        result.normalizeLikelihoods(log10globalReadMismappingRate, symmetricallyNormalizeAllelesToReference);
        filterPoorlyModeledEvidence(result, dynamicDisqualification, expectedErrorRatePerBase, readDisqualificationScale);
        return result;
    }


    /**
     * Creates a new GATKRead with the source read's header, read group and mate
     * information, but with the following fields set to user-supplied values:
     *  - Read Bases
     *  - Base Qualities
     *  - Base Insertion Qualities
     *  - Base Deletion Qualities
     *
     *  Cigar string is empty (not-null)
     *
     * Use this method if you want to create a new GATKRead based on
     * another GATKRead, but with modified bases and qualities
     *
     * @param read a read to copy the header from
     * @param readBases an array containing the new bases you wish use in place of the originals
     * @param baseQualities an array containing the new base qualities you wish use in place of the originals
     * @param baseInsertionQualities an array containing the new base insertion qaulities
     * @param baseDeletionQualities an array containing the new base deletion qualities
     * @return a read with modified bases and qualities, safe for the GATK
     */
    private static GATKRead createQualityModifiedRead(final GATKRead read,
                                                      final byte[] readBases,
                                                      final byte[] baseQualities,
                                                      final byte[] baseInsertionQualities,
                                                      final byte[] baseDeletionQualities) {
        Utils.validateArg( baseQualities.length == readBases.length && baseInsertionQualities.length == readBases.length && baseDeletionQualities.length == readBases.length,
                () -> String.format("Read bases and read quality arrays aren't the same size: Bases: %d vs Base Q's: %d vs Insert Q's: %d vs Delete Q's: %d.",
                        readBases.length, baseQualities.length, baseInsertionQualities.length, baseDeletionQualities.length));

        final GATKRead processedRead = ReadUtils.emptyRead(read);
        processedRead.setBases(readBases);
        processedRead.setBaseQualities(baseQualities);
        ReadUtils.setInsertionBaseQualities(processedRead, baseInsertionQualities);
        ReadUtils.setDeletionBaseQualities(processedRead, baseDeletionQualities);
        return processedRead;
    }

    /**
     * Initialize our pairHMM with parameters appropriate to the haplotypes and reads we're going to evaluate
     *
     * After calling this routine the PairHMM will be configured to best evaluate all reads in the samples
     * against the set of haplotypes
     *
     * @param haplotypes a non-null list of haplotypes
     * @param perSampleReadList a mapping from sample -> reads
     */
    private void initializePairHMM(final List<PartiallyDeterminedHaplotype> haplotypes, final Map<String, List<GATKRead>> perSampleReadList) {
        final int readMaxLength = perSampleReadList.entrySet().stream().flatMap(e -> e.getValue().stream()).mapToInt(GATKRead::getLength).max().orElse(0);
        final int haplotypeMaxLength = haplotypes.stream().mapToInt(h -> h.getBases().length).max().orElse(0);

        // initialize arrays to hold the probabilities of being in the match, insertion and deletion cases
        pdPairHMM.initialize(haplotypes, perSampleReadList, readMaxLength, haplotypeMaxLength);
    }

    private void computeReadLikelihoods(final LikelihoodMatrix<GATKRead, PartiallyDeterminedHaplotype> likelihoods) {
        // Modify the read qualities by applying the PCR error model and capping the minimum base,insertion,deletion qualities
        final List<GATKRead> processedReads = modifyReadQualities(likelihoods.evidence());

        if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
            for(int counter = 0; counter < processedReads.size(); counter++) {
                GATKRead read = processedReads.get(counter);
                HaplotypeCallerGenotypingDebugger.println("Range for Overlaps to Variants for consideration: "+rangeForReadOverlapToDeterminedBases);
                HaplotypeCallerGenotypingDebugger.println("read "+counter +": "+read.getName()+" cigar: "+read.getCigar()+" mapQ: "+read.getMappingQuality()+" loc: ["+read.getStart() +"-"+ read.getEnd()+"] unclippedloc: ["+read.getUnclippedStart()+"-"+read.getUnclippedEnd()+"]");
                HaplotypeCallerGenotypingDebugger.println(Arrays.toString(read.getBaseQualitiesNoCopy()));
            }
        }
        // Run the PairHMM to calculate the log10 likelihood of each (processed) reads' arising from each haplotype
        pdPairHMM.computeLog10Likelihoods(likelihoods, processedReads, inputScoreImputator, rangeForReadOverlapToDeterminedBases);
    }

    /**
     * Pre-processing of the reads to be evaluated at the current location from the current sample.
     * We apply the PCR Error Model, and cap the minimum base, insertion, and deletion qualities of each read.
     * Modified copies of reads are packed into a new list, while original reads are retained for downstream use
     *
     * @param reads The original list of unmodified reads
     * @return processedReads. A new list of reads, in the same order, whose qualities have been altered by PCR error model and minimal quality thresholding
     */
    private List<GATKRead> modifyReadQualities(final List<GATKRead> reads) {
        final List<GATKRead> result = new ArrayList<>(reads.size());

        for (final GATKRead read : reads) {
            final GATKRead maybeUnclipped = modifySoftclippedBases ? read : ReadClipper.hardClipSoftClippedBases(read);    //Clip the bases here to remove hap
            final byte[] readBases = maybeUnclipped.getBases();
            //For FRD we want to have scores for reads that don't overlap
            final SimpleInterval readSpan = new SimpleInterval(read.getContig(), read.getUnclippedStart(), read.getUnclippedEnd());

            // NOTE -- must clone anything that gets modified here so we don't screw up future uses of the read
            //Using close here is justified - it's an array of primitives.
            final byte[] readQuals = maybeUnclipped.getBaseQualities().clone();
            final byte[] readInsQuals = ReadUtils.getBaseInsertionQualities(maybeUnclipped).clone();
            final byte[] readDelQuals = ReadUtils.getBaseDeletionQualities(maybeUnclipped).clone();

            applyPCRErrorModel(readBases, readInsQuals, readDelQuals);
            capMinimumReadQualities(maybeUnclipped, readQuals, readInsQuals, readDelQuals, baseQualityScoreThreshold, disableCapReadQualitiesToMapQ);

            // Store the actual qualities
            read.setTransientAttribute(HMM_BASE_QUALITIES_TAG, readQuals);
            // Create a new copy of the read and sets its base qualities to the modified versions.
            GATKRead qualityModifiedRead = createQualityModifiedRead(maybeUnclipped, readBases, readQuals, readInsQuals, readDelQuals);
            qualityModifiedRead.setTransientAttribute(UNCLIPPED_ORIGINAL_SPAN_ATTR, readSpan);
            result.add(qualityModifiedRead);
        }
        return result;
    }

    private static void capMinimumReadQualities(final GATKRead read, final byte[] readQuals, final byte[] readInsQuals, final byte[] readDelQuals, final byte baseQualityScoreThreshold, final boolean disableCapReadQualitiesToMapQ) {
        for( int i = 0; i < readQuals.length; i++ ) {
            if (!disableCapReadQualitiesToMapQ) {
                readQuals[i] = (byte) Math.min(0xff & readQuals[i], read.getMappingQuality()); // cap base quality by mapping quality, as in UG
            }
            readQuals[i] =    setToFixedValueIfTooLow( readQuals[i],    baseQualityScoreThreshold,             QualityUtils.MIN_USABLE_Q_SCORE );
            readInsQuals[i] = setToFixedValueIfTooLow( readInsQuals[i], QualityUtils.MIN_USABLE_Q_SCORE,       QualityUtils.MIN_USABLE_Q_SCORE );
            readDelQuals[i] = setToFixedValueIfTooLow( readDelQuals[i], QualityUtils.MIN_USABLE_Q_SCORE,       QualityUtils.MIN_USABLE_Q_SCORE );
        }
    }

    private static byte setToFixedValueIfTooLow(final byte currentVal, final byte minQual, final byte fixedQual){
        return currentVal < minQual ? fixedQual : currentVal;
    }

    private static Map<GATKRead, byte[]> buildGapContinuationPenalties(final List<GATKRead> reads, final byte gapPenalty) {
        final Map<GATKRead, byte[]> result = new HashMap<>(reads.size());
        reads.stream().forEach(read -> result.put(read, Utils.dupBytes(gapPenalty, read.getLength())));
        return result;
    }

    /* --------------------------------------------------------------------------------
    *
    * Experimental attempts at PCR error rate modeling
    *
    -------------------------------------------------------------------------------- */

    private byte[] pcrIndelErrorModelCache;
    private PairHMMInputScoreImputator inputScoreImputator;

    private void initializePCRErrorModel() {

        inputScoreImputator = dragstrParams == null
                ? StandardPairHMMInputScoreImputator.newInstance(constantGCP)
                : DragstrPairHMMInputScoreImputator.of(dragstrParams) ;

        if ( !pcrErrorModel.hasRateFactor() ) {
            return;
        }

        pcrIndelErrorModelCache = new byte[ReadLikelihoodCalculationEngine.MAX_REPEAT_LENGTH + 1];

        final double rateFactor = pcrErrorModel.getRateFactor();

        for(int i = 0; i <= ReadLikelihoodCalculationEngine.MAX_REPEAT_LENGTH; i++ ) {
            pcrIndelErrorModelCache[i] = getErrorModelAdjustedQual(i, rateFactor);
        }
    }

    static byte getErrorModelAdjustedQual(final int repeatLength, final double rateFactor) {
        return (byte) Math.max(MIN_ADJUSTED_QSCORE, MathUtils.fastRound(INITIAL_QSCORE - Math.exp(repeatLength / (rateFactor * Math.PI)) + 1.0));
    }

    @VisibleForTesting
    void applyPCRErrorModel( final byte[] readBases, final byte[] readInsQuals, final byte[] readDelQuals ) {
        if ( pcrErrorModel == PairHMMLikelihoodCalculationEngine.PCRErrorModel.NONE ) {
            return;
        }

        for ( int i = 1; i < readBases.length; i++ ) {
            final int repeatLength = ReadLikelihoodCalculationEngine.findTandemRepeatUnits(readBases, i-1).getRight();
            readInsQuals[i-1] = (byte) Math.min(0xff & readInsQuals[i - 1], 0xff & pcrIndelErrorModelCache[repeatLength]);
            readDelQuals[i-1] = (byte) Math.min(0xff & readDelQuals[i - 1], 0xff & pcrIndelErrorModelCache[repeatLength]);
        }
    }

}
