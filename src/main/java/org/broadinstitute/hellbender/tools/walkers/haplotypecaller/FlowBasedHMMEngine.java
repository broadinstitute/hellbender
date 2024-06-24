package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParams;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.pairhmm.FlowBasedPairHMM;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMMInputScoreImputation;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMMInputScoreImputator;
import org.broadinstitute.hellbender.utils.read.FlowBasedReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.haplotype.FlowBasedHaplotype;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.function.ToDoubleFunction;

/**
 *  Flow Based HMM, intended to incorporate the scoring model of the {@link FlowBasedAlignmentLikelihoodEngine} while allowing for frame-shift insertions
 *  and deletions for better genotyping.
 */
public class FlowBasedHMMEngine implements ReadLikelihoodCalculationEngine {
    private final double readDisqualificationScale;
    private final boolean dynamicReadDisqualification;
    private double log10globalReadMismappingRate;
    private final double expectedErrorRatePerBase;
    private PairHMMLikelihoodCalculationEngine.PCRErrorModel pcrErrorModel;
    final FlowBasedArgumentCollection fbargs;

    @VisibleForTesting
    public static final double INITIAL_QSCORE = 40.0;

    private final FlowBasedPairHMM flowPairHMM;

    public static final byte MIN_USABLE_Q_SCORE_DEFAULT = 6;
    private static final int MIN_ADJUSTED_QSCORE = 10;

    private byte minUsableIndelScoreToUse;

    private PairHMMInputScoreImputator inputScoreImputator;
    private final DragstrParams dragstrParams;
    private byte constantGCP;
    private final byte flatInsertionPenalty;
    private final byte flatDeletionPenalty;


    /**
     * Default constructor
     * @param fbargs - arguments
     * @param log10globalReadMismappingRate - probability for wrong mapping (maximal contribution of the read to data likelihood)
     * @param expectedErrorRatePerBase - the expected rate of random sequencing errors for a read originating from its true haplotype.
     * @param pcrErrorModel
     */
    public FlowBasedHMMEngine(final FlowBasedArgumentCollection fbargs, final byte constantGCP, final double log10globalReadMismappingRate, final double expectedErrorRatePerBase,
                              final PairHMMLikelihoodCalculationEngine.PCRErrorModel pcrErrorModel, final DragstrParams dragstrParams, final boolean dynamicReadDisqualification, final double readDisqualificationScale,
                              final int minUsableIndelScoreToUse, final byte flatDeletionPenalty, final byte flatInsertionPenalty) {
        this.fbargs = fbargs;
        this.log10globalReadMismappingRate = log10globalReadMismappingRate;
        this.expectedErrorRatePerBase = expectedErrorRatePerBase;
        this.readDisqualificationScale = readDisqualificationScale;
        this.dynamicReadDisqualification = dynamicReadDisqualification;
        this.pcrErrorModel = pcrErrorModel;
        this.flowPairHMM = new FlowBasedPairHMM();
        this.dragstrParams = dragstrParams;
        this.constantGCP = constantGCP;
        this.flatDeletionPenalty = flatDeletionPenalty;
        this.flatInsertionPenalty = flatInsertionPenalty;
        this.minUsableIndelScoreToUse = (byte)minUsableIndelScoreToUse;

        // validity checks
        if (fbargs.flowOrderCycleLength != 4) {
            throw new GATKException("FlowBasedHMMEngine requires flow order of 4 elements but cycle length is specified as " + fbargs.flowOrderCycleLength);
        }

        initializePCRErrorModel();
    }

    public AlleleLikelihoods<GATKRead, Haplotype> computeReadLikelihoods(final List<Haplotype> haplotypeList,
                                                                         final SAMFileHeader hdr,
                                                                         final SampleList samples,
                                                                         final Map<String, List<GATKRead>> perSampleReadList,
                                                                         final boolean filterPoorly){
        return computeReadLikelihoods(haplotypeList, hdr, samples, perSampleReadList, filterPoorly, true);
    }


    public AlleleLikelihoods<GATKRead, Haplotype> computeReadLikelihoods(final List<Haplotype> haplotypeList,
                                                                         final SAMFileHeader hdr,
                                                                         final SampleList samples,
                                                                         final Map<String, List<GATKRead>> perSampleReadList,
                                                                         final boolean filterPoorly,
                                                                         final boolean normalizeLikelihoods) {
        Utils.nonNull(samples, "samples is null");
        Utils.nonNull(perSampleReadList, "perSampleReadList is null");
        Utils.nonNull(haplotypeList, "haplotypeList is null");

        final AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(haplotypeList);

        // Add likelihoods for each sample's reads to our result
        final AlleleLikelihoods<GATKRead, Haplotype> result = new AlleleLikelihoods<>(samples, haplotypes, perSampleReadList);
        final int sampleCount = result.numberOfSamples();
        for (int i = 0; i < sampleCount; i++) {
            computeReadLikelihoods(result.sampleMatrix(i), hdr);
        }
        if (normalizeLikelihoods) {
            result.normalizeLikelihoods(log10globalReadMismappingRate, true);
        }
        if ( filterPoorly ) {
            filterPoorlyModeledEvidence(result, dynamicReadDisqualification, expectedErrorRatePerBase, readDisqualificationScale);
        }

        return result;
    }


    /** Calculate minimal likelihood that reasonably matching read can get. We divide errors into "expected", e.g.
     * hmer indels without 0->1/1->0 errors. Those on average occur with expectedErrorRate and "catastrophic' e.g.
     * 1->0 and 0->1 errors (or large hmer indels). Those occur with probability fbargs.fillingValue.
     * If the read has more than 3 expected and more than 2 "catastrophic" errors it will probably be deemed unfit to the
     * haplotype
     *
     * @param expectedErrorRate  error rate for expected errors.
     * @return minimal likelihood for the read to be considered not poorly modeled
     */
    @Override
    public ToDoubleFunction<GATKRead> log10MinTrueLikelihood(final double expectedErrorRate, final boolean capLikelihoods) {
        final double log10ErrorRate = Math.log10(expectedErrorRate);
        final double largeEventErrorRate = 0.001; // error rate for non-hmer/snv errors that are not seq. errors.
        final double log10catastrophicErrorRate = Math.log10(largeEventErrorRate);

        return read -> {
            final double maxErrorsForRead = Math.max(3.0, Math.ceil(read.getLength() * expectedErrorRate));
            final double maxCatastrophicErrorsForRead = Math.max(2.0, Math.ceil(read.getLength() * largeEventErrorRate));
            return maxErrorsForRead * log10ErrorRate + maxCatastrophicErrorsForRead*log10catastrophicErrorRate;
        };
    }

    private byte[] pcrIndelErrorModelCache;

    private void initializePCRErrorModel() {

        inputScoreImputator = dragstrParams == null
                ? NonSymmetricalPairHMMInputScoreImputator.newInstance(constantGCP, flatInsertionPenalty, flatDeletionPenalty)
                : DragstrPairHMMInputScoreImputator.of(dragstrParams) ;

        if ( !pcrErrorModel.hasRateFactor() ) {
            return;
        }

        pcrIndelErrorModelCache = new byte[ReadLikelihoodCalculationEngine.MAX_REPEAT_LENGTH + 1];

        final double rateFactor = pcrErrorModel.getRateFactor();

        for( int i = 0; i <= ReadLikelihoodCalculationEngine.MAX_REPEAT_LENGTH; i++ ) {
            pcrIndelErrorModelCache[i] = getErrorModelAdjustedQual(i, rateFactor, minUsableIndelScoreToUse);
        }

    }


    // TODO these methods are ripped whole cloth from the PairHMMLikelihoodsCalculationEngine and should be uinfied/fixed if possible
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

    private static void capMinimumReadIndelQualities(final byte[] readInsQuals, final byte[] readDelQuals, final byte minUsableQualScore) {
        for( int i = 0; i < readInsQuals.length; i++ ) {
            readInsQuals[i] = (byte)Math.max(readInsQuals[i], minUsableQualScore);
            readDelQuals[i] = (byte)Math.max(readDelQuals[i], minUsableQualScore);
        }
    }

    static byte getErrorModelAdjustedQual(final int repeatLength, final double rateFactor, final byte minUsableIndelScoreToUse) {
        return (byte) Math.max(minUsableIndelScoreToUse, MathUtils.fastRound(INITIAL_QSCORE - Math.exp(repeatLength / (rateFactor * Math.PI)) + 1.0));
    }

    /**
     * Initialize our flowPairHMM with parameters appropriate to the haplotypes and reads we're going to evaluate
     *
     * After calling this routine the PairHMM will be configured to best evaluate all reads in the samples
     * against the set of haplotypes
     *  @param haplotypes a non-null list of haplotypes
     * @param perSampleReadList a mapping from sample -> reads
     */
    private void initializeFlowPairHMM(final List<FlowBasedHaplotype> haplotypes, final List<FlowBasedRead> perSampleReadList) {
        initializePCRErrorModel();

        final int readMaxLength = perSampleReadList.stream().mapToInt(FlowBasedRead::getKeyLength).max().orElse(0);
        final int haplotypeMaxLength = haplotypes.stream().mapToInt(h -> h.getKeyLength()).max().orElse(0);

        // initialize arrays to hold the probabilities of being in the match, insertion and deletion cases
        flowPairHMM.initialize(readMaxLength, haplotypeMaxLength);;
    }

    /**
     * Compute read likelihoods for a single sample
     * @param likelihoods Single sample likelihood matrix
     * @param hdr SAM header that corresponds to the sample
     */
    private void computeReadLikelihoods(final LikelihoodMatrix<GATKRead, Haplotype> likelihoods,
                                        final SAMFileHeader hdr) {

        final List<FlowBasedRead> processedReads = new ArrayList<>(likelihoods.evidenceCount());
        final List<FlowBasedHaplotype> processedHaplotypes = new ArrayList<>(likelihoods.numberOfAlleles());

        // establish flow order based on the first evidence. Note that all reads belong to the same sample (group)
        final FlowBasedReadUtils.ReadGroupInfo rgInfo = (likelihoods.evidenceCount() != 0)
                ? FlowBasedReadUtils.getReadGroupInfo(hdr, likelihoods.evidence().get(0))
                : null;
        final String flowOrder = (rgInfo != null)
                ? rgInfo.flowOrder.substring(0, fbargs.flowOrderCycleLength)
                : FlowBasedReadUtils.findFirstUsableFlowOrder(hdr, fbargs);

        //convert all reads to FlowBasedReads (i.e. parse the matrix of P(call | haplotype) for each read from the BAM)
        for (int i = 0 ; i < likelihoods.evidenceCount(); i++) {
            final GATKRead rd = likelihoods.evidence().get(i);

            // create a flow based read
            final FlowBasedRead fbRead = new FlowBasedRead(rd, flowOrder, rgInfo.maxClass, fbargs);
            fbRead.applyAlignment();

            //TODO This currently supports any ScoreImputator in GATK but will probably need a custom one to handle flow based data in the future:
            final PairHMMInputScoreImputation inputScoreImputation = inputScoreImputator.impute(fbRead);
            final byte[] readInsQuals = inputScoreImputation.insOpenPenalties();
            final byte[] readDelQuals = inputScoreImputation.delOpenPenalties();
            final byte[] overallGCP = inputScoreImputation.gapContinuationPenalties();

            applyPCRErrorModel(fbRead.getBases(), readInsQuals, readDelQuals);
            capMinimumReadIndelQualities(readInsQuals, readDelQuals, minUsableIndelScoreToUse);

            fbRead.setReadInsQuals(readInsQuals);
            fbRead.setReadDelQuals(readDelQuals);
            fbRead.setOverallGCP(overallGCP);

            processedReads.add(fbRead);
        }

        //same for the haplotypes - each haplotype is converted to FlowBasedHaplotype

        for (int i = 0; i < likelihoods.numberOfAlleles(); i++){
            final FlowBasedHaplotype fbh = new FlowBasedHaplotype(likelihoods.alleles().get(i), flowOrder);
            processedHaplotypes.add(fbh);
        }

        initializeFlowPairHMM(processedHaplotypes, processedReads);


        flowPairHMM.computeLog10LikelihoodsFlowBased(likelihoods, processedReads, processedHaplotypes);
    }

    @Override
    public void close() {}

}
