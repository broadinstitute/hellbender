package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMUtils;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.FlowBasedProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.groundtruth.SeriesStats;
import org.broadinstitute.hellbender.utils.read.*;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

@CommandLineProgramProperties(
        summary = "This program converts the flow qualities that Ultima Genomics CRAM reports to more conventional base qualities. " +
                "Specifically, the generated quality will report the probability that a specific base is a sequencing error mismatch, " +
                "while auxilary tags qa, qt, qg, qc report specific probability that a specific base X is a A->X error. " +
                "Since mismatch error in flow-based chemistries can only occur as a result of several indel errors, " +
                "we implemented several strategies to estimate the probability of a mismatch which can be specified" +
                "using the svnq-mode parameter: " +
                "Legacy - the quality value from flow matrix is used. " +
                "Optimistic - assuming that the probability of the indel errors are p1 and p2, then snvq=p1*p2 - assuming they always coincide. " +
                "Pessimistic - snvq=(1-p1)*(1-p2) - assuming they never coincide. " +
                "Geometric - snvq=sqrt(Optimistic*Pessimistic) - i.e. the geometric mean of the optimistic and Pessimistic modes. " +
                "The Geometric is set as the default mode",
        oneLineSummary = "Add SNV Quality to the flow-based CRAM",
        programGroup = FlowBasedProgramGroup.class
)

@DocumentedFeature
@ExperimentalFeature
public final class AddFlowSNVQuality extends ReadWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "File to which reads should be written")
    @WorkflowOutput(optionalCompanions={StandardArgumentDefinitions.OUTPUT_INDEX_COMPANION})
    public GATKPath output = null;
    private SAMFileGATKReadWriter outputWriter;

    @ArgumentCollection
    public FlowBasedArgumentCollection fbargs = new FlowBasedArgumentCollection();

    @ArgumentCollection
    public AddFlowSNVQualityArgumentCollection aqArgs = new AddFlowSNVQualityArgumentCollection();

    public static final char PHRED_ASCII_BASE = '!';

    public static final int ERROR_PROB_BAND_1LESS = 0;
    public static final int ERROR_PROB_BAND_KEY = 1;
    public static final int ERROR_PROB_BAND_1MORE = 2;
    public static final int ERROR_PROB_BANDS = 3;

    public double minLikelihoodProbRate = 1e-6;
    public int maxQualityScore = 60;

    // locals
    private SeriesStats                         inputQualStats = new SeriesStats();
    private SeriesStats                         outputBQStats = new SeriesStats();
    private SeriesStats                         outputQAltStats = new SeriesStats();
    private SeriesStats                         outputQCalledStats = new SeriesStats();
    private SeriesStats                         outputSumPStats = new SeriesStats();

    // private class to hold the base probabilities and SNVQ probabilties for a read
    class ReadProbs {
        double[] baseProbs;
        double[][] snvqProbs; // length of first dimension is flow order length
    }

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        outputWriter = createSAMWriter(output, true);
    }

    @Override
    public void closeTool() {
        super.closeTool();
        if ( outputWriter != null ) {
            outputWriter.close();
        }

        try {
            if ( aqArgs.debugCollectStatsInto != null )
                printStats(aqArgs.debugCollectStatsInto);
        } catch (IOException e) {
            throw new GATKException("", e);
        }
    }

    @Override
    public void apply(final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        // include supplementary alignments?
        if ( read.isSupplementaryAlignment() && !aqArgs.keepSupplementaryAlignments ) {
            return;
        }

        // include qc-failed reads?
        if (  read.failsVendorQualityCheck() && !aqArgs.includeQcFailedReads ) {
            return;
        }

        // collect input stats
        if ( aqArgs.debugCollectStatsInto != null ) {
            collectInputStats(read);
        }

        // add SNVQ attributes
        addBaseQuality(read, getHeaderForReads(), aqArgs.maxPhredScore, fbargs);

        // collect output stats
        if ( aqArgs.debugCollectStatsInto != null ) {
            collectOutputStats(read);
            if ( aqArgs.debugReadName.size() != 0 && aqArgs.debugReadName.contains(read.getName()) ) {
                dumpOutputRead(read);
            }
        }

        // write read to output
        outputWriter.addRead(read);
    }

    private void collectInputStats(GATKRead read) {
        for ( byte q : read.getBaseQualitiesNoCopy() ) {
            inputQualStats.add(q);
        }
    }

    private void collectOutputStats(GATKRead read) {
        if ( aqArgs.outputQualityAttribute != null ) {
            if (read.hasAttribute(aqArgs.outputQualityAttribute)) {
                for (byte q : read.getAttributeAsString(aqArgs.outputQualityAttribute).getBytes()) {
                    outputBQStats.add(SAMUtils.fastqToPhred((char) q));
                }
            }
        } else {
            for (byte q : read.getBaseQualitiesNoCopy()) {
                outputBQStats.add(q);
            }
        }
        final FlowBasedReadUtils.ReadGroupInfo rgInfo = FlowBasedReadUtils.getReadGroupInfo(getHeaderForReads(), read);
        final byte[] bases = read.getBasesNoCopy();
        final double[] sumP = new double[bases.length];
        for ( int i = 0 ; i < 4 ; i++ ) {
            byte altBase = rgInfo.flowOrder.getBytes()[i];
            String attrValue = read.getAttributeAsString(attrNameForNonCalledBase(altBase));
            int ofs = 0;
            for ( byte q : attrValue.getBytes() ) {
                if ( bases[ofs] != altBase ) {
                    outputQAltStats.add(SAMUtils.fastqToPhred((char)q));
                } else {
                    outputQCalledStats.add(SAMUtils.fastqToPhred((char)q));
                }
                sumP[ofs] += Math.pow(10.0, SAMUtils.fastqToPhred((char)q) / -10.0);
                ofs++;

            }
        }
        for ( double p : sumP ) {
            outputSumPStats.add(p);
        }
    }

    // dump read as a csv for easier analysis
    private void dumpOutputRead(GATKRead read) {

        try {
            // open file
            final String fname = aqArgs.debugCollectStatsInto + "." + read.getName() + ".csv";
            logger.info("dumping read into: " + fname);
            final PrintWriter pw = new PrintWriter(fname);

            // write header
            final StringBuilder hdr = new StringBuilder();
            hdr.append("pos,base,qual,tp,t0,bq");
            final FlowBasedReadUtils.ReadGroupInfo rgInfo = FlowBasedReadUtils.getReadGroupInfo(getHeaderForReads(), read);
            for (int i = 0; i < 4; i++) {
                hdr.append(",");
                hdr.append(attrNameForNonCalledBase(rgInfo.flowOrder.charAt(i)));
            }
            hdr.append(",qCalled");
            pw.println(hdr);

            // access data
            final byte[] bases = read.getBasesNoCopy();
            final byte[] quals = read.getBaseQualitiesNoCopy();
            final byte[] tp = read.getAttributeAsByteArray(FlowBasedRead.FLOW_MATRIX_TAG_NAME);
            final byte[] t0 = read.getAttributeAsByteArray(FlowBasedRead.FLOW_MATRIX_T0_TAG_NAME);
            final byte[] bq = (aqArgs.outputQualityAttribute != null)
                    ? read.getAttributeAsString(aqArgs.outputQualityAttribute).getBytes()
                    : null;
            final byte[][] qX = new byte[4][];
            for (int i = 0; i < 4; i++) {
                qX[i] = read.getAttributeAsString(attrNameForNonCalledBase(rgInfo.flowOrder.charAt(i))).getBytes();
            }

            // write lines
            List<String> line = new LinkedList<>();
            for (int pos = 0; pos < bases.length; pos++) {
                line.clear();

                // position
                line.add(Integer.toString(pos));

                // base, qual
                line.add(Character.toString(bases[pos]));
                line.add(Integer.toString(quals[pos]));

                // tp,t0,bq
                line.add(Integer.toString(tp[pos]));
                line.add(Integer.toString(SAMUtils.fastqToPhred((char)t0[pos])));
                if ( bq != null ) {
                    line.add(Integer.toString(SAMUtils.fastqToPhred((char) bq[pos])));
                }

                // qX
                int calledIndex = -1;
                for (int i = 0; i < 4; i++) {
                    line.add(Integer.toString(SAMUtils.fastqToPhred((char)qX[i][pos])));
                    if ( bases[pos] == rgInfo.flowOrder.charAt(i) ) {
                        calledIndex = i;
                    }
                }

                // qCalled
                if ( calledIndex >= 0 ) {
                    line.add(Integer.toString(SAMUtils.fastqToPhred((char)qX[calledIndex][pos])));
                } else {
                    line.add("-1");
                }

                // write the line
                pw.println(StringUtils.join(line, ','));
            }

            // close file
            pw.close();
        } catch (IOException e) {
            throw new GATKException("", e);
        }
    }

    private void printStats(final String fname) throws IOException {

        inputQualStats.csvWrite(fname + ".inputQual.csv");
        outputBQStats.csvWrite(fname + ".outputBQ.csv");
        outputQAltStats.csvWrite(fname + ".outputQAlt.csv");
        outputQCalledStats.csvWrite(fname + ".outputQCalled.csv");
        outputSumPStats.csvWrite(fname + ".outputSumP.csv");
    }

    static public String attrNameForNonCalledBase(byte nonCalledBase) {
        return attrNameForNonCalledBase((char)nonCalledBase);
    }

    static public String attrNameForNonCalledBase(char nonCalledBase) {
        return "q" + Character.toLowerCase(nonCalledBase);
    }

    public void addBaseQuality(final GATKRead read, final SAMFileHeader hdr, double maxPhredScore, FlowBasedArgumentCollection fbargs) {

        // take in phred score limit
        if ( !Double.isNaN(maxPhredScore) ) {
            maxQualityScore = (int)maxPhredScore;
            minLikelihoodProbRate = Math.pow(10, -maxPhredScore / 10.0);
        }

        // convert to a flow base read
        final FlowBasedReadUtils.ReadGroupInfo rgInfo = FlowBasedReadUtils.getReadGroupInfo(hdr, read);
        final FlowBasedRead fbRead = new FlowBasedRead(read, rgInfo.flowOrder, rgInfo.maxClass, fbargs);
        final int flowOrderLength = FlowBasedReadUtils.calcFlowOrderLength(rgInfo.flowOrder);

        // generate base and snvq probabilities for the read
        final ReadProbs readProbs = generateFlowReadBaseAndSNVQErrorProbabilities(fbRead, flowOrderLength, rgInfo.flowOrder.getBytes());

        // install in read
        if ( aqArgs.outputQualityAttribute != null ) {
            read.setAttribute(aqArgs.outputQualityAttribute, new String(convertErrorProbToFastq(readProbs.baseProbs)));
        } else {
            read.setBaseQualities(convertErrorProbToPhred(readProbs.baseProbs));
        }
        for ( int i = 0 ; i < flowOrderLength ; i++ ) {
            final String name = AddFlowSNVQuality.attrNameForNonCalledBase(rgInfo.flowOrder.charAt(i));
            read.setAttribute(name, new String(convertErrorProbToFastq(readProbs.snvqProbs[i])));
        }
    }

    // Not using SamUtils function since normally an error probability can not be zero.
    // still, this method is called to convert base quality as well as snvq, which is computed.
    // the following check is a safety, in case snvq produces a zero.
    private char[] convertErrorProbToFastq(double[] errorProb) {

        byte[] phred = convertErrorProbToPhred(errorProb);
        return SAMUtils.phredToFastq(phred).toCharArray();
    }

    // Not using SamUtils function since normally an error probability can not be zero.
    // still, this method is called to convert base quality as well as snvq, which is computed.
    // the following check is a safety, in case snvq produces a zero.
    private byte[] convertErrorProbToPhred(double[] errorProb) {

        final byte[] phred = new byte[errorProb.length];
        for ( int i = 0 ; i < errorProb.length ; i++ ) {

            if ( errorProb[i] == 0 ) {
                phred[i] = (byte)maxQualityScore;
            } else {
                final double p = errorProb[i];
                phred[i] =  (byte)Math.round(-10 * Math.log10(p));
            }
        }
        return phred;
    }

    /**
     * generate base and snvq probabilties for a read.
     *
     * @param fbRead a flow based read
     * @param flowOrderLength number of bases in flow order (essentially number of valid base values)
     * @param flowOrder the flow order itself (which can be the size of flowOrderLength or a repeat of it
     *
     * @return an instance of a private class containing the base probabilities as well as the snvq probabilities
     */
    private ReadProbs generateFlowReadBaseAndSNVQErrorProbabilities(final FlowBasedRead fbRead, final int flowOrderLength, byte[] flowOrder) {

        /**
         * access key and error probabilities
         * for a description of the flow probabilities see {@link FlowBasedRead#flowMatrix}
         */
        final int[]       key = fbRead.getKey();
        final double[][]  errorProbBands = extractErrorProbBands(fbRead, minLikelihoodProbRate);

        // allocate returned prob arrays
        final double[]    baseProbs = new double[fbRead.getBasesNoCopy().length];
        final double[][]  snvqProbs = new double[flowOrderLength][];
        for ( int i = 0 ; i < snvqProbs.length ; i++ ) {
            snvqProbs[i] = new double[baseProbs.length];
        }

        // loop over hmers via flow key
        int               base = 0;
        Map<Byte, Double> allBaseProb0 = new LinkedHashMap<>();
        Map<Byte, Double> allBaseProb1 = new LinkedHashMap<>();

        for ( int flow = 0 ; flow < key.length ; flow++ ) {
            if ( key[flow] != 0 ) {

                // establish initial stat
                allBaseProb0.clear();
                allBaseProb1.clear();
                int flow_i = (flow % flowOrderLength);

                // establish hmer quality
                final int         hmerLength = key[flow];
                final double[]    hmerBaseErrorProbs = generateHmerBaseErrorProbabilities(key, errorProbBands, flow, flowOrderLength, flowOrder, allBaseProb0, allBaseProb1);

                // install value in first byte of the hmer
                baseProbs[base++] = hmerBaseErrorProbs[0];  // first base, or only base in case of a single base hmer
                for ( int i = 0 ; i < flowOrderLength ; i++ ) {
                    if ( allBaseProb0.containsKey(flowOrder[i]) ) {
                        snvqProbs[i][base - 1] = allBaseProb0.get(flowOrder[i]);
                    } else if ( i != flow_i ) {
                        snvqProbs[i][base - 1] = minLikelihoodProbRate;
                    }
                }

                // for hmers longer than 1
                if ( hmerLength > 1 ) {

                    // skip all but last (leave with zero error probability)
                    base += (hmerLength - 2);

                    // fill last base from computed error probability
                    baseProbs[base++] = hmerBaseErrorProbs[1]; // last base, if hmer is longer than 1

                    for ( int i = 0 ; i < flowOrderLength ; i++ ) {
                        if ( allBaseProb1.containsKey(flowOrder[i]) ) {
                            final double p = allBaseProb1.get(flowOrder[i]);
                            for ( int j = 0 ; j < hmerLength - 1 ; j++ ) {
                                snvqProbs[i][base - 1 - j] = (j == 0) ? p : minLikelihoodProbRate; // all but last get the min prob
                            }
                        } else if ( i != flow_i ) {
                            for ( int j = 0 ; j < hmerLength - 1 ; j++ ) {
                                snvqProbs[i][base - 1 - j] = minLikelihoodProbRate;
                            }
                        }
                    }
                }

                // override result for the last base with the original hmer error probability
                if ( base == baseProbs.length ) {
                    baseProbs[base - 1] = errorProbBands[ERROR_PROB_BAND_KEY][flow];
                }
            }
        }

        // adjust probability of called bases so that sum will be 1, also enforce min prob
        final byte[] bases = fbRead.getBasesNoCopy();
        for ( int ofs = 0 ; ofs < bases.length ; ofs++ ) {

            // go through alt bases and accumulate p, find out index of called bases (in flow order)
            final byte calledBase = bases[ofs];
            double altP = 0;
            int calledIndex = -1;
            for (int i = 0; i < flowOrderLength; i++) {
                if ( calledBase != flowOrder[i] ) {
                    snvqProbs[i][ofs] = Math.max(minLikelihoodProbRate, snvqProbs[i][ofs]);
                    altP += snvqProbs[i][ofs];
                } else {
                    calledIndex = i;
                }
            }
            if ( calledBase < 0 ) {
                throw new GATKException(String.format("failed to locate called base %c in flow order %s", (char)calledBase, flowOrder));
            }

            // install probability in called base slot
            snvqProbs[calledIndex][ofs] = Math.max(0, 1 - altP);

            // at this point, bq becomes trivial (?)
            baseProbs[ofs] = 1 - snvqProbs[calledIndex][ofs];
        }

        // build return value
        ReadProbs readProbs = new ReadProbs();
        readProbs.baseProbs = baseProbs;
        readProbs.snvqProbs = snvqProbs;
        return readProbs;
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

    @VisibleForTesting
    protected double[] generateHmerBaseErrorProbabilities(final int[] key, final double[][] errorProbBands, final int flow,
                                                          final int flowOrderLength, byte[] flowOrder,
                                                          Map<Byte, Double> allBaseProb0, Map<Byte, Double> allBaseProb1) {

        // result is left/right error probabilities
        final double[]          errorProbs = new double[2];
        final int               hmerLength = key[flow];

        errorProbs[0] = generateSidedHmerBaseErrorProbability(key, errorProbBands, flow, -1, flowOrderLength, flowOrder, allBaseProb0);
        if ( hmerLength != 1 ) {
            errorProbs[1] = generateSidedHmerBaseErrorProbability(key, errorProbBands, flow, 1, flowOrderLength, flowOrder, allBaseProb1);
        }

        return errorProbs;
    }

    private double generateSidedHmerBaseErrorProbability(final int[] key, final double[][] errorProbBands, final int flow, final int sideIncr,
                                                         final int flowOrderLength, final byte[] flowOrder, final Map<Byte, Double> allBaseProb) {

        // create a key slice of the area around the flow/hmer.
        final int minIndex = Math.max(flow - (flowOrderLength - 1), 0);
        final int maxIndex = Math.min(flow + (flowOrderLength - 1), key.length - 1);
        final int[] slice = Arrays.copyOfRange(key, minIndex, maxIndex + 1);
        final int hmerLength = key[flow];

        // walk the flows towards the side until (and including) the first non-zero key
        // on hmers of length 1 we walk both sides
        final class SliceInfo {
            int[] slice;
            byte altByte;
            int sideFlow;
        }
        final List<SliceInfo> slices = new LinkedList<>();
        final int[]           incrs = (hmerLength != 1)
                ? new int[] { sideIncr }
                : new int[] { sideIncr, -sideIncr};
        for (int incr : incrs) {
            for (int sideFlow = flow + incr; sideFlow >= 0 && sideFlow < key.length; sideFlow += incr) {

                // side flow can no overflow the slice
                if ( sideFlow < minIndex || sideFlow > maxIndex ) {
                    break;
                }

                // create a alternative key slice by incrementing sideFlow and decrementing flow
                final int[] altSlice = Arrays.copyOf(slice, slice.length);
                altSlice[sideFlow - minIndex] += 1;
                altSlice[flow - minIndex] -= 1;
                if ( sliceIsValidForConsideration(altSlice, flowOrderLength) ) {
                    SliceInfo si = new SliceInfo();
                    si.slice = altSlice;
                    si.altByte = flowOrder[sideFlow % flowOrderLength];
                    si.sideFlow = sideFlow;
                    slices.add(si);
                }

                // is the sideFlow (the first encountered) non-zero? if so, break
                if (key[sideFlow] != 0) {
                    break;
                }
            }
        }

        // at this point, we have a list of valid slices. figure out the error probability for each of them
        // and compute the base quality
        final double keyP = sliceProbs(slice, minIndex, key, errorProbBands, flow, flow)[0];
        double sumP = keyP;
        for ( final SliceInfo si : slices ) {
            final double[] sliceP = sliceProbs(si.slice, minIndex, key, errorProbBands, flow, si.sideFlow);
            if ( allBaseProb != null ) {
                allBaseProb.put(si.altByte, getSnvq(sliceP[0], sliceP[1], sliceP[2], aqArgs.snvMode));
            }
            sumP += sliceP[0];
        }
        final double ep = 1 - (keyP / sumP);

        return ep;
    }

    static double getSnvq(final double sliceP, final double p1, final double p2, AddFlowSNVQualityArgumentCollection.SnvqModeEnum snvMode) {
        if ( snvMode == AddFlowSNVQualityArgumentCollection.SnvqModeEnum.Legacy ) {
            return sliceP;
        } else if ( snvMode == AddFlowSNVQualityArgumentCollection.SnvqModeEnum.Optimistic ) {
            return (p1 * p2);
        } else if ( snvMode == AddFlowSNVQualityArgumentCollection.SnvqModeEnum.Pessimistic ) {
            return (1 - (1 - p1) * (1 - p2));
        } else if ( snvMode == AddFlowSNVQualityArgumentCollection.SnvqModeEnum.Geometric ) {
            return Math.sqrt((p1 * p2) * (1 - (1 - p1) * (1 - p2)));
        } else {
            throw new GATKException("unknown snvqMode: " +  snvMode);
        }
    }

    // compute probability for a slice
    private static double[] sliceProbs(final int[] slice, final int minIndex, final int[] key, final double[][] errorProbBands,
                                     final int flow, final int sideFlow) {

        double accumulatedP = 1.0;
        double p1 = 0.0;
        double p2 = 0.0;
        int key_i = minIndex;
        for ( int i = 0 ; i < slice.length ; i++, key_i++ ) {
            final int hmer = key[key_i];
            final int band;
            if ( slice[i] == (hmer - 1) ) {
                band = ERROR_PROB_BAND_1LESS;
            } else if ( slice[i] == (hmer + 1) ) {
                band = ERROR_PROB_BAND_1MORE;
            } else if ( slice[i] == hmer ){
                band = ERROR_PROB_BAND_KEY;
            } else {
                throw new GATKException("slice[i] and hmer are too far apart: " + slice[i] + " " + hmer);
            }
            final double p = errorProbBands[band][key_i];
            accumulatedP *= p;

            // collect p1/p2 (flow and sideFlow probs)
            if ( key_i == flow ) {
                p1 = p;
            }
            if ( key_i == sideFlow ) {
                p2 = p;
            }
        }

        return new double[] {accumulatedP, p1, p2};
    }

    static boolean sliceIsValidForConsideration(final int[] slice, final int flowOrderLength) {

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

