package org.broadinstitute.hellbender.utils.pairhmm;

import com.intel.gkl.pdhmm.IntelPDHMM;
import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.gatk.nativebindings.pairhmm.HaplotypeDataHolder;
import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeArguments;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.PDPairHMMLikelihoodCalculationEngine;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.haplotype.PartiallyDeterminedHaplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Class for performing the pair HMM for global alignment using AVX instructions contained in a native shared library.
 */
public final class VectorLoglessPairPDHMM extends LoglessPDPairHMM {

    //TODO what should this size be? make this configurable
    private final int GKL_COMP_BATCH_SIZE = 1000;

    /**
     * Type for implementation of VectorLoglessPairHMM
     */
    public enum Implementation {
//        /**
//         * AVX-accelerated version of PairHMM
//         */
//        AVX,
        /**
         * OpenMP multi-threaded AVX-accelerated version of PairHMM
         */
        OMP,
    }

    private static final Logger logger = LogManager.getLogger(VectorLoglessPairPDHMM.class);
    private static long pairHMMTime = 0;
    private static long pairHMMSetupTime = 0;

    private final IntelPDHMM pairPDHmm;

    /**
     * Create a VectorLoglessPairHMM
     *
     * @param implementation    which implementation to use (AVX or OMP)
     * @param args              arguments to the native GKL implementation
     */
    public VectorLoglessPairPDHMM(Implementation implementation, PairHMMNativeArguments args) throws UserException.HardwareFeatureException {
        final boolean isSupported;

        switch (implementation) {
//            case AVX:
//                pairHmm = new IntelPairHmm();
//                isSupported = pairHmm.load(null);
//                if (!isSupported) {
//                    throw new UserException.HardwareFeatureException("Machine does not support AVX PairHMM.");
//                }
//                break;

            case OMP:
                pairPDHmm = new IntelPDHMM();
                isSupported = pairPDHmm.load(null);
                if (!isSupported) {
                    throw new UserException.HardwareFeatureException("Machine does not support OpenMP AVX PairHMM.");
                }
                break;

            default:
                throw new UserException.HardwareFeatureException("Unknown PairHMM implementation.");
        }
        //TODO no initialize argument?
        //pairPDHmm.initialize(args);
    }

//    /**
//     * {@inheritDoc}
//     */
////    @Override
//    public void initialize(final int readMaxLength, final int haplotypeMaxLength ) {
//        super.initialize(readMaxLength, haplotypeMaxLength);
//
//        branchMatchMatrix = new double[paddedMaxReadLength][paddedMaxHaplotypeLength];
//        branchInsertionMatrix = new double[paddedMaxReadLength][paddedMaxHaplotypeLength];
//        branchDeletionMatrix = new double[paddedMaxReadLength][paddedMaxHaplotypeLength];
//        // do not need to call super.initialize()
//        int numHaplotypes = haplotypes.size();
//        mHaplotypeDataArray = new HaplotypeDataHolder[numHaplotypes];
//        int idx = 0;
//        haplotypeToHaplotypeListIdxMap.clear();
//        for (final Haplotype currHaplotype : haplotypes) {
//            mHaplotypeDataArray[idx] = new HaplotypeDataHolder();
//            mHaplotypeDataArray[idx].haplotypeBases = currHaplotype.getBases();
//            haplotypeToHaplotypeListIdxMap.put(currHaplotype, idx);
//            ++idx;
//        }
//    }

//    /**
//     * {@inheritDoc}
//     */
//    @Override
//    public void initialize(final List<Haplotype> haplotypes, final Map<String, List<GATKRead>> perSampleReadList,
//                           final int readMaxLength, final int haplotypeMaxLength) {
//
//    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void computeLog10Likelihoods(final LikelihoodMatrix<GATKRead, PartiallyDeterminedHaplotype> logLikelihoods,
                                        final List<GATKRead> processedReads,
                                        final PairHMMInputScoreImputator inputScoreImputator,
                                        final int rangeForReadOverlapToDeterminedBases) {
        if (processedReads.isEmpty()) {
            return;
        }
        if(doProfiling) {
            startTime = System.nanoTime();
        }
        // (re)initialize the pairHMM only if necessary
        final int readMaxLength = findMaxReadLength(processedReads);
        final int haplotypeMaxLength = findMaxAlleleLength(logLikelihoods.alleles());
        final int readCount = processedReads.size();
        final List<PartiallyDeterminedHaplotype> alleles = logLikelihoods.alleles();
        final int alleleCount = alleles.size();
        final int totalComps = readCount * alleleCount;
        if (!initialized || readMaxLength > maxReadLength || haplotypeMaxLength > maxHaplotypeLength) {
            initialize(readMaxLength, haplotypeMaxLength);
        }

        final GKLSubmissionBatchManager submissionManager = new GKLSubmissionBatchManager(readMaxLength, haplotypeMaxLength);
        final List<PDHMMComparisonContainer> submittedComps = new ArrayList<>();
        final List<PDHMMComparisonContainer> nullComps = new ArrayList<>();

        // Iterate over all reads and alleles, and submit the comparisons to the submission manager (or record which ones were missed)
        int currentTestcase = 0;
        for(int readIndex = 0; readIndex < readCount; readIndex++) {
            GATKRead read = processedReads.get(readIndex);
            final PairHMMInputScoreImputation inputScoreImputation = inputScoreImputator.impute(read);
            // peek at the next haplotype in the list (necessary to get nextHaplotypeBases, which is required for caching in the array implementation)
            for (int a = 0; a < alleleCount; a++) {
                final PartiallyDeterminedHaplotype allele = alleles.get(a);

                PDHMMComparisonContainer comp = new PDHMMComparisonContainer(read, readIndex, inputScoreImputation, allele, a);
                // if we aren't apply the optimization (range == -1) or if the read overlaps the determined region for calling:
                if (rangeForReadOverlapToDeterminedBases < 0 || allele.getMaximumExtentOfSiteDeterminedAlleles()
                        .overlapsWithMargin((Locatable) read.getTransientAttribute(PDPairHMMLikelihoodCalculationEngine.UNCLIPPED_ORIGINAL_SPAN_ATTR), rangeForReadOverlapToDeterminedBases + 1)) { // the +1 here is us erring on the side of caution
                    // append individual comparision case to the full arrays
                    submissionManager.addComp(comp);
                    submittedComps.add(comp);
                } else {
                    nullComps.add(comp);
                }
                currentTestcase++;
            }
        }

        mLogLikelihoodArray = new double[readCount * alleleCount];
        double[] mLogLikelihoodArray = submissionManager.getCompMatrixOutput();

        assert(mLogLikelihoodArray.length == submittedComps.size()); //If this fails, something went wrong with array filling indexing and we should not proceed
        assert(submittedComps.size() + nullComps.size() == totalComps); //If this fails, something went wrong with array filling indexing and we should not proceed
        assert(currentTestcase == totalComps); //If this fails, something went wrong with array filling indexing and we should not proceed

        // Overwrite the likelihoods for reads that have invalid likelihoods due to not spanning the determined region
        for(int i = 0; i < submittedComps.size(); i++) {
            PDHMMComparisonContainer comp = submittedComps.get(i);
            double lk = mLogLikelihoodArray[i];
            logLikelihoods.set(comp.alleleIndex, comp.readIndex, lk);
        }

        // Blank out the null comparisions that we skipped over anyway
        for (PDHMMComparisonContainer comp : nullComps) {
            logLikelihoods.set(comp.alleleIndex, comp.readIndex, Double.NEGATIVE_INFINITY);
        }
    }

    @Override
    public void close() {
        //pairPDHmm.done(); // TODO evidently this doesn't have a close command?
        if (doProfiling) {
            logger.info("Time spent in setup for JNI call : " + (pairHMMSetupTime * 1e-9));
            logger.info("Time spent in JNI call : " + (pairHMMTime * 1e-9));
        }
        super.close();
    }



    private class GKLSubmissionBatchManager {

        private static long perepTimeStart = doProfiling ? System.nanoTime() : 0;
        private static long nativeTimeStart = 0;

        // (re)initialize the pairHMM only if necessary
        final int readMaxLength;// = findMaxReadLength(processedReads);
        final int haplotypeMaxLength;// = findMaxAlleleLength(logLikelihoods.alleles());
        private final List<PDHMMComparisonContainer> submissionQueue = new ArrayList<>();

        // Storing the internal submission chunk outputs to hide from consumer
        private List<double[]> compMatrixOutputChunks = new ArrayList<>(GKL_COMP_BATCH_SIZE);

        GKLSubmissionBatchManager(final int readMaxLength, final int haplotypeMaxLength) {
            this.readMaxLength = readMaxLength;
            this.haplotypeMaxLength = haplotypeMaxLength;
        }

        /**
         * Method that adds a comparison to the batch to be submitted to the GKL native PairHMM code.
         */
        void addComp(final PDHMMComparisonContainer container) {
            if (submissionQueue.size() >= GKL_COMP_BATCH_SIZE) {
                submitBatch();
            }
            submissionQueue.add(container);
        }

        /**
         * Method that returns the output of the GKL native PairHMM code as a single array.
         * @return the output of the GKL native PairHMM likelihood calculations as a merged array in the order they were submitted
         */
        double[] getCompMatrixOutput() {
            // empty the queue if there are any reads left
            if (!submissionQueue.isEmpty()) {
                submitBatch();
            }

            int totalLength = compMatrixOutputChunks.stream().mapToInt(arr -> arr.length).sum();
            // Create the merged array
            double[] mergedArray = new double[totalLength];
            // Copy each array from the list into the merged array
            int startIndex = 0;
            for (double[] array : compMatrixOutputChunks) {
                System.arraycopy(array, 0, mergedArray, startIndex, array.length);
                startIndex += array.length;
            }
            return mergedArray;
        }

        /**
         * Method that clears the Queue and actually handles constructing and submitting the arrays to the GKL native PairHMM code.
         */
        private void submitBatch() {
            final int currentBatchSize = submissionQueue.size();

            final int hapArraySize = currentBatchSize * haplotypeMaxLength;
            final int readArraySize = currentBatchSize * readMaxLength;
            byte[] alleleBasesFull = new byte[hapArraySize];
            byte[] allelePDBasesFull = new byte[hapArraySize];
            byte[] readBasesFull = new byte[readArraySize];
            byte[] readQualsFull = new byte[readArraySize];
            byte[] readInsQualsFull = new byte[readArraySize];
            byte[] readDelQualsFull = new byte[readArraySize];
            byte[] overallGCPFull = new byte[readArraySize];
            long[] hapLengthsFull = new long[currentBatchSize];
            long[] readLengthsFull = new long[currentBatchSize];

            // Populate the array information from the queued read/haplotype information
            for (int i = 0; i < submissionQueue.size(); i++) {
                final PDHMMComparisonContainer container = submissionQueue.get(i);
                // pull information form the queued reads and alleles (NOTE: the intention is that none of these operations are copy operations and are fairly light-weight)
                final GATKRead read = container.read;
                final PartiallyDeterminedHaplotype allele = container.allele;
                final PairHMMInputScoreImputation inputScoreImputation = container.inputScoreImputation;
                final byte[] readBases = read.getBases();
                final byte[] readQuals = read.getBaseQualities();
                final byte[] readInsQuals = inputScoreImputation.insOpenPenalties();
                final byte[] readDelQuals = inputScoreImputation.delOpenPenalties();
                final byte[] overallGCP = inputScoreImputation.gapContinuationPenalties();
                final byte[] alleleBases = allele.getBases();

                // Copy the information into the full arrays that will be submitted to the GKL
                System.arraycopy(alleleBases, 0, alleleBasesFull, i * haplotypeMaxLength,
                        alleleBases.length);
                System.arraycopy(allele.getAlternateBases(), 0, allelePDBasesFull, i * haplotypeMaxLength,
                        allele.getAlternateBases().length);
                System.arraycopy(readBases, 0, readBasesFull, i * readMaxLength, readBases.length);
                System.arraycopy(readQuals, 0, readQualsFull, i * readMaxLength, readQuals.length);
                System.arraycopy(readInsQuals, 0, readInsQualsFull, i * readMaxLength,
                        readInsQuals.length);
                System.arraycopy(readDelQuals, 0, readDelQualsFull, i * readMaxLength,
                        readDelQuals.length);
                System.arraycopy(overallGCP, 0, overallGCPFull, i * readMaxLength, overallGCP.length);
                hapLengthsFull[i] = alleleBases.length;
                readLengthsFull[i] = readBases.length;
            }

            if (doProfiling) {
                pairHMMSetupTime += (System.nanoTime() - perepTimeStart);
                nativeTimeStart = System.nanoTime();
            }

            // Call the JNI PairPDHMM
            compMatrixOutputChunks.add(pairPDHmm.computePDHMM(alleleBasesFull, allelePDBasesFull,
                    readBasesFull, readQualsFull, readInsQualsFull, readDelQualsFull,
                    overallGCPFull, hapLengthsFull, readLengthsFull, currentBatchSize, haplotypeMaxLength, readMaxLength));

            if (doProfiling) {
                pairHMMTime += (System.nanoTime() - nativeTimeStart);
                startTime = System.nanoTime(); // count all of the time between batched calls as setup time
            }

            // Empty the Queues for the next batch
            submissionQueue.clear();
        }
    }

    // Container for holding submission information
    record PDHMMComparisonContainer(GATKRead read, int readIndex,
                                    PairHMMInputScoreImputation inputScoreImputation,
                                    PartiallyDeterminedHaplotype allele,
                                    int alleleIndex) {
    }
}
