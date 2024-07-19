package org.broadinstitute.hellbender.utils.pairhmm;

import com.intel.gkl.pairhmm.IntelPairHmm;
import com.intel.gkl.pairhmm.IntelPairHmmOMP;
import com.intel.gkl.pdhmm.IntelPDHMM;
import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.gatk.nativebindings.pairhmm.HaplotypeDataHolder;
import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeArguments;
import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeBinding;
import org.broadinstitute.gatk.nativebindings.pairhmm.ReadDataHolder;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.PDPairHMMLikelihoodCalculationEngine;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.haplotype.PartiallyDeterminedHaplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Class for performing the pair HMM for global alignment using AVX instructions contained in a native shared library.
 */
public final class VectorLoglessPairPDHMM extends LoglessPDPairHMM {

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
    private long threadLocalSetupTimeDiff = 0;
    private long pairHMMSetupTime = 0;

    private final IntelPDHMM pairPDHmm;

    //Hold the mapping between haplotype and index in the list of Haplotypes passed to initialize
    //Use this mapping in computeLikelihoods to find the likelihood value corresponding to a given Haplotype
    private final Map<Haplotype, Integer> haplotypeToHaplotypeListIdxMap = new LinkedHashMap<>();
    private HaplotypeDataHolder[] mHaplotypeDataArray;

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

        int hapArraySize = totalComps * haplotypeMaxLength;
        int readArraySize = totalComps * readMaxLength;

        // The following arrays are used to pass data to the native code, they expect the data to be in a contiguous arrays
        byte[] alleleBasesFull = new byte[hapArraySize];
        byte[] allelePDBasesFull = new byte[hapArraySize];
        byte[] readBasesFull = new byte[readArraySize];
        byte[] readQualsFull = new byte[readArraySize];
        byte[] readInsQualsFull = new byte[readArraySize];
        byte[] readDelQualsFull = new byte[readArraySize];
        byte[] overallGCPFull = new byte[readArraySize];
        double[] expectedFull = new double[totalComps];
        long[] hapLength = new long[totalComps];
        long[] readLength = new long[totalComps];

        int idx = 0;
        int readIndex = 0;
        int currentTestcase = 0;
        for(final GATKRead read : processedReads){
            final PairHMMInputScoreImputation inputScoreImputation = inputScoreImputator.impute(read);
            final byte[] readBases = read.getBases();

            final byte[] readQuals = read.getBaseQualities();
            final byte[] readInsQuals = inputScoreImputation.insOpenPenalties();
            final byte[] readDelQuals = inputScoreImputation.delOpenPenalties();
            final byte[] overallGCP = inputScoreImputation.gapContinuationPenalties();

            // peek at the next haplotype in the list (necessary to get nextHaplotypeBases, which is required for caching in the array implementation)
            for (int a = 0; a < alleleCount; a++) {
                final PartiallyDeterminedHaplotype allele = alleles.get(a);
                final byte[] alleleBases = allele.getBases();
                // if we aren't apply the optimization (range == -1) or if the read overlaps the determined region for calling:
//                if (rangeForReadOverlapToDeterminedBases < 0 || allele.getMaximumExtentOfSiteDeterminedAlleles()
//                        .overlapsWithMargin((Locatable) read.getTransientAttribute(PDPairHMMLikelihoodCalculationEngine.UNCLIPPED_ORIGINAL_SPAN_ATTR), rangeForReadOverlapToDeterminedBases + 1)) { // the +1 here is us erring on the side of caution
                    // append individual comparision case to the full arrays
                    System.arraycopy(alleleBases, 0, alleleBasesFull, currentTestcase * haplotypeMaxLength,
                            alleleBases.length);
                    System.arraycopy(allele.getAlternateBases(), 0, allelePDBasesFull, currentTestcase * haplotypeMaxLength,
                            allele.getAlternateBases().length);
                    System.arraycopy(readBases, 0, readBasesFull, currentTestcase * readMaxLength, readBases.length);
                    System.arraycopy(readQuals, 0, readQualsFull, currentTestcase * readMaxLength, readQuals.length);
                    System.arraycopy(readInsQuals, 0, readInsQualsFull, currentTestcase * readMaxLength,
                            readInsQuals.length);
                    System.arraycopy(readDelQuals, 0, readDelQualsFull, currentTestcase * readMaxLength,
                            readDelQuals.length);
                    System.arraycopy(overallGCP, 0, overallGCPFull, currentTestcase * readMaxLength, overallGCP.length);

                    // TODO we should in future check that we can apply the no-op optimzation for non-overlapping reads
                    hapLength[currentTestcase] = alleleBases.length;
                    readLength[currentTestcase] = readBases.length;
                    currentTestcase++;

                    // Otherwise we record that the likelihood array is bogus here for later validation and set it to -infinity
//                } else {
//                    lk = Double.NEGATIVE_INFINITY;
//                }

            }
            readIndex++;
        }

        assert(currentTestcase == totalComps); //If this fails, something went wrong with array filling indexing and we should not proceed

        if (doProfiling) {
            threadLocalSetupTimeDiff = (System.nanoTime() - startTime);
        }
        mLogLikelihoodArray = new double[readCount * alleleCount];
        double[] mLogLikelihoodArray = pairPDHmm.computePDHMM(alleleBasesFull, allelePDBasesFull,
                readBasesFull, readQualsFull, readInsQualsFull, readDelQualsFull,
                overallGCPFull, hapLength, readLength, totalComps, maxHaplotypeLength, maxReadLength);


        // Overwrite the likelihoods for reads that have invalid likelihoods due to not spanning the determined region
        //TODO this is a hack, move it
        idx = 0;
        readIndex = 0;
        int overwriteInx = 0;
        for(final GATKRead read : processedReads){
            for (int a = 0; a < alleleCount; a++) {
                final PartiallyDeterminedHaplotype allele = alleles.get(a);
                double lk;
                if (rangeForReadOverlapToDeterminedBases < 0 || allele.getMaximumExtentOfSiteDeterminedAlleles()
                        .overlapsWithMargin((Locatable) read.getTransientAttribute(PDPairHMMLikelihoodCalculationEngine.UNCLIPPED_ORIGINAL_SPAN_ATTR), rangeForReadOverlapToDeterminedBases + 1)) { // the +1 here is us erring on the side of caution
                    lk = mLogLikelihoodArray[overwriteInx];
                } else {
                    lk = Double.NEGATIVE_INFINITY;
                }
                logLikelihoods.set(a, readIndex, lk);
                mLogLikelihoodArray[idx++] = lk;
                // TODO we don't have these elements available at this point so we can't write them out
                //writeToResultsFileIfApplicable(readBases, readQuals, readInsQuals, readDelQuals, overallGCP, alleleBases, allele.getAlternateBases(), lk);
            }

        }

        if (doProfiling) {
            threadLocalPairHMMComputeTimeDiff = (System.nanoTime() - startTime);
            pairHMMComputeTime += threadLocalPairHMMComputeTimeDiff;
            pairHMMSetupTime += threadLocalSetupTimeDiff;
        }
    }


    @Override
    public void close() {
        //pairPDHmm.done(); // TODO evidently this doesn't have a close command?
        if (doProfiling)
            logger.info("Time spent in setup for JNI call : " + (pairHMMSetupTime * 1e-9));
        super.close();
    }
}
