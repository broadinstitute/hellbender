package org.broadinstitute.hellbender.tools.walkers.groundtruth;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.FlowBasedProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import org.broadinstitute.hellbender.utils.read.FlowBasedReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

/**
 * Convert quality string that reports probability of indel errors (in flow chemistry) to the quality string that approximates probability of mismatch error.
 *
 * <p>
 * </p>
 *
 * <p>
 * The reference is strictly required when handling CRAM files.
 * </p>
 *
 * <h3> Input </h3>
 * <ul>
 *     <li> Coordinate-sorted and indexed SAM/BAM/CRAM </li>
 * </ul>
 *
 * <h3> Output </h3>
 * <ul>
 *     <li> Coordinate-sorted and indexed SAM/BAM/CRAM </li>
 * </ul>
 *
 * {@GATK.walkertype ReadWalker}
 */
@CommandLineProgramProperties(
        summary = "Prints reads from the input SAM/BAM/CRAM file to the SAM/BAM/CRAM file while adding a base quality attribute.",
        oneLineSummary = "Add base quality attribute to reads in in the SAM/BAM/CRAM file",
        programGroup = FlowBasedProgramGroup.class
)
@ExperimentalFeature
@WorkflowProperties
public final class AddFlowBaseQuality extends ReadWalker {

    public static final String MINIMAL_ERROR_RATE_LONG_NAME = "minimal-error-rate";
    public static final String MAXIMAL_QUALITY_SCORE_LONG_NAME = "maximal-quality-score";
    public static final String REPLACE_QUALITY_MODE_LONG_NAME = "replace-quality-mode";
    public static final String BASE_QUALITY_ATTRIBUTE_NAME = "XQ";
    public static final String OLD_QUALITY_ATTRIBUTE_NAME = "OQ";
    public static final char PHRED_ASCII_BASE = '!';

    public static final int ERROR_PROB_BAND_1LESS = 0;
    public static final int ERROR_PROB_BAND_KEY = 1;
    public static final int ERROR_PROB_BAND_1MORE = 2;
    public static final int ERROR_PROB_BANDS = 3;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="Write output to this file")
    @WorkflowOutput(optionalCompanions={StandardArgumentDefinitions.OUTPUT_INDEX_COMPANION})
    public GATKPath output;
    private SAMFileGATKReadWriter outputWriter;

    @Argument(fullName = MINIMAL_ERROR_RATE_LONG_NAME, doc = "lower floor for flow error rate values")
    public double minErrorRate = 1e-3;

    @Argument(fullName = MAXIMAL_QUALITY_SCORE_LONG_NAME, doc = "clip quality score to the given value)")
    public int maxQualityScore = 93;

    @Argument(fullName = REPLACE_QUALITY_MODE_LONG_NAME, doc = "replace existing base qualities while saving previous qualities to OQ (when true) or simply write to BQ (when false) ")
    public boolean replaceQualityMode = false;

    @ArgumentCollection
    public FlowBasedArgumentCollection fbargs = new FlowBasedArgumentCollection();

    @Override
    public void onTraversalStart() {
        outputWriter = createSAMWriter(output, true);

        // do not trim the boundary homopolymers as the default
        fbargs.keepBoundaryFlows = true;
    }

    @Override
    public void apply(final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        outputWriter.addRead(addBaseQuality(read));
    }

    @Override
    public void closeTool() {
        if ( outputWriter != null ) {
            outputWriter.close();
        }
    }

    /**
     * this is the actual 'worker' of the tool
     */
    private GATKRead addBaseQuality(final GATKRead read) {

        // convert to a flow base read
        final FlowBasedReadUtils.ReadGroupInfo rgInfo = FlowBasedReadUtils.getReadGroupInfo(getHeaderForReads(), read);
        final FlowBasedRead   fbRead = new FlowBasedRead(read, rgInfo.flowOrder, rgInfo.maxClass, fbargs);
        final int flowOrderLength = calcFlowOrderLength(rgInfo.flowOrder);

        // generate base quality
        final double[]      baseErrorProb = generateBaseErrorProbability(fbRead, flowOrderLength);
        final byte[]        phred = convertErrorProbToPhred(baseErrorProb);

        // install in read
        if ( !replaceQualityMode ) {
            read.setAttribute(BASE_QUALITY_ATTRIBUTE_NAME, convertPhredToString(phred));
        } else {
            read.setAttribute(OLD_QUALITY_ATTRIBUTE_NAME, convertPhredToString(read.getBaseQualitiesNoCopy()));
            read.setBaseQualities(phred);
        }

        // return reused read
        return read;
    }

    private byte[] convertErrorProbToPhred(double[] errorProb) {

        final byte[] phred = new byte[errorProb.length];
        for ( int i = 0 ; i < errorProb.length ; i++ ) {
            if ( errorProb[i] == 0 ) {
                phred[i] = (byte)(maxQualityScore);
            } else {
                phred[i] = (byte)Math.min((maxQualityScore), (int) (-10 * Math.log10(errorProb[i])));
            }
        }
        return phred;
    }

    private String convertPhredToString(final byte[] phred) {

        final byte[] s = new byte[phred.length];
        for ( int i = 0 ; i < phred.length ; i++ ) {
            s[i] = (byte)(phred[i] + PHRED_ASCII_BASE);
        }
        return new String(s);
    }

    private int calcFlowOrderLength(String flowOrder) {

        final int i = flowOrder.indexOf(flowOrder.charAt(0), 1);

        return (i < 0) ? flowOrder.length() : i;
    }

    private double[] generateBaseErrorProbability(final FlowBasedRead fbRead, final int flowOrderLength) {

        /**
         * access key and error probabilities
         * for a description of the flow probabilities see {@link FlowBasedRead#flowMatrix}
         */
        final int[]       key = fbRead.getKey();
        final double[][]  errorProbBands = extractErrorProbBands(fbRead, minErrorRate);
        final double[]    result = new double[fbRead.getBasesNoCopy().length];

        // loop over hmers via flow key
        int               base = 0;
        for ( int flow = 0 ; flow < key.length ; flow++ ) {
            if ( key[flow] != 0 ) {

                // establish hmer quality
                final int         hmerLength = key[flow];
                final double[]    hmerBaseErrorProbs = generateHmerBaseErrorProbabilities(key, errorProbBands, flow, flowOrderLength);

                // fill in

                // first base of the read gets the original error probability of the hmer, otherwise use computed error probability
                if ( base == 0 ) {
                    result[base++] = errorProbBands[ERROR_PROB_BAND_KEY][flow];
                } else {
                    result[base++] = hmerBaseErrorProbs[0];  // first base, or only base in case of a single base hmer
                }

                // for hmers longer than 1
                if ( hmerLength > 1 ) {

                    // skip all but last (leave with zero error probability)
                    base += (hmerLength - 2);

                    // fill last base from computed error probability
                    result[base++] = hmerBaseErrorProbs[1]; // last base, if hmer is longer than 1
                }

                // override result for the last base with the orignal hmer error probability
                if ( base == result.length ) {
                    result[base - 1] = errorProbBands[ERROR_PROB_BAND_KEY][flow];
                }
            }
        }

        return result;
    }

    // extract error probability bands. middle (1) band is the key prob.
    // lower (0) and high (2) are corresponding to -1 and +1 in hmer lengths
    private static double[][] extractErrorProbBands(final FlowBasedRead flowRead, final double minValue) {

        // access key
        final int[] key = flowRead.getKey();

        // allocate result
        double[][] result = new double[ERROR_PROB_BANDS][];
        for ( int i = 0 ; i < result.length ; i++ ) {
            result[i] = new double[key.length];
        }

        for ( int i = 0 ; i < key.length ; i++ ) {

            // extract key probability
            result[ERROR_PROB_BAND_KEY][i] = Math.max(flowRead.getProb(i, key[i]), minValue);

            // -1
            if ( key[i] > 0 ) {
                result[ERROR_PROB_BAND_1LESS][i] = Math.max(flowRead.getProb(i, key[i] - 1), minValue);
            } else {
                result[ERROR_PROB_BAND_1LESS][i] = minValue;
            }

            // +1
            if ( key[i] < flowRead.getMaxHmer() ) {
                result[ERROR_PROB_BAND_1MORE][i] = Math.max(flowRead.getProb(i, key[i] + 1), minValue);
            } else {
                result[ERROR_PROB_BAND_1MORE][i] = minValue;
            }
        }

        return result;
    }

    /**
     * The following functions estimate the error probability for an hmer. specifically two error
     * probability values are generated: one for the first base of the hmer and another for the
     * rest of its bases.
     *
     * The computation itself is performed in a subsequent function: generateSidedHmerBaseErrorProbability
     * It iterates over the possible valid combinations of errors and sums them up.
     *
     * @param key - key (hmer length) in flow space
     * @param errorProbBands - for each flow (position in the key) three error probabilities are provided:
     *                       [0] - for the hmer being one base shorter
     *                       [1] - for the hmer to be at its length
     *                       [2] - for the hmer to be one base longer
     * @param flow - the flow (index) for which to generate the probabilities (0 <= flow < key.length)
     * @param flowOrderLength - the cycle length of of the flow order (usually 4)
     * @return an array of two probabilities:
     * [0] - probability for the first base of the hmer
     * [1] - probability for the rest of the bases of the hmer
     */
    @VisibleForTesting
    protected static double[] generateHmerBaseErrorProbabilities(final int[] key, final double[][] errorProbBands, final int flow, final int flowOrderLength) {

        // result is left/right error probabilities
        final double[]          errorProbs = new double[2];
        final int               hmerLength = key[flow];

        errorProbs[0] = generateSidedHmerBaseErrorProbability(key, errorProbBands, flow, -1, flowOrderLength);
        if ( hmerLength != 1 ) {
            errorProbs[1] = generateSidedHmerBaseErrorProbability(key, errorProbBands, flow, 1, flowOrderLength);
        }

        return errorProbs;
    }

    private static double generateSidedHmerBaseErrorProbability(final int[] key, final double[][] errorProbBands, final int flow, final int sideIncr, final int flowOrderLength) {

        // create a key slice of the area around the flow/hmer.
        final int minIndex = Math.max(flow - flowOrderLength + 1, 0);
        final int maxIndex = Math.min(flow + flowOrderLength - 1, key.length - 1);
        final int[] slice = Arrays.copyOfRange(key, minIndex, maxIndex + 1);
        final int hmerLength = key[flow];

        // walk the flows towards the side until (and including) the first non-zero key
        // on hmers of length 1 we walk both sides
        final List<int[]>     slices = new LinkedList<>();
        final int[]           incrs = (hmerLength != 1)
                                            ? new int[] { sideIncr }
                                            : new int[] { sideIncr, -sideIncr};
        for (int incr : incrs) {
            for (int sideFlow = flow + incr; sideFlow >= 0 && sideFlow < key.length; sideFlow += incr) {

                // create a alternative key slice by incrementing sideFlow and decrementing flow
                final int[] altSlice = Arrays.copyOf(slice, slice.length);
                altSlice[sideFlow - minIndex] += 1;
                altSlice[flow - minIndex] -= 1;
                if (sliceIsValid(altSlice, flowOrderLength)) {
                    slices.add(altSlice);
                }

                // is the sideFlow non-zero? if so, break
                if (key[sideFlow] != 0) {
                    break;
                }
            }
        }

        // at this point, we have a list of valid slices. figure out the error probability for each of them
        // and compute the base quality
        final double keyP = sliceProb(slice, minIndex, key, errorProbBands);
        double sumP = keyP;
        for ( int[] s : slices ) {
            sumP += sliceProb(s, minIndex, key, errorProbBands);
        }
        final double ep = 1 - (keyP / sumP);

        return ep;
    }

    // compute probability for a slice
    private static double sliceProb(int[] slice, int minIndex, int[] key, double[][] errorProbBands) {

        double p = 1.0;
        for ( int i = 0 ; i < slice.length ; i++ ) {
            final int band;
            if ( slice[i] < key[i + minIndex] ) {
                band = ERROR_PROB_BAND_1LESS;
            } else if ( slice[i] > key[i + minIndex] ) {
                band = ERROR_PROB_BAND_1MORE;
            } else {
                band = ERROR_PROB_BAND_KEY;
            }
            p *= errorProbBands[band][i + minIndex];
        }
        return p;
    }

    private static boolean sliceIsValid(final int[] slice, final int flowOrderLength) {

        // look for strings of consecutive zeros in length of flowOrderLength - 1
        int     consecutiveZeros = 0;
        for ( int key : slice ) {
            if ( key != 0 ) {
                consecutiveZeros = 0;
            } else {
                consecutiveZeros++;
                if ( consecutiveZeros >= (flowOrderLength - 1) ) {
                    return false;
                }
            }
        }

        // if here, not found -> valid
        return true;
    }
}
