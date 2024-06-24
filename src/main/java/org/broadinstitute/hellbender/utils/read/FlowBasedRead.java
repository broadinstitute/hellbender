package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.*;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.FlowBasedHMMEngine;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Tail;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;

import java.io.*;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;
/**
 * Adds flow information to the usual GATKRead. In addition to the usual read data this class keeps flowMatrix,
 * that contains probabilities for alternative hmer calls.
 *
 * Main function deals with parsing flow-specific QUAL representation readFlowMatrix.
 * Note that there is a lot of code that deals with other various formats of the representation (e.g. when the matrix
 * is coded in the tags of the BAM and is given in flow space). This code is not used in production, but was used in
 * development and testing
 *
 * A common usage pattern is to covert a GATKRead into a FlowBasedRead. Additionally
 * a SAMRecord can also be converted into a FlowBasedRead. Follows a common usage pattern:
 *
 * For a self contained example of a usage pattern, see {@link FlowBasedReadUtils#convertToFlowBasedRead(GATKRead, SAMFileHeader)}
 *
 **/

public class FlowBasedRead extends SAMRecordToGATKReadAdapter implements GATKRead, Serializable {

    public static final int  MAX_CLASS = 12; //note - this is a historical value to support files with max class is not defined in the header, it is expected that mc tag exists in the CRAM
    public static final String DEFAULT_FLOW_ORDER = "TGCA";
    private static final long serialVersionUID = 42L;
    private static final Logger logger = LogManager.getLogger(FlowBasedRead.class);
    private static final OneShotLogger vestigialOneShotLogger = new OneShotLogger(FlowBasedRead.class);

    // constants
    static private final int MINIMAL_READ_LENGTH = 10; // check if this is the right number
    private final double MINIMAL_CALL_PROB = 0.1;

    // constants for clippingTagContains.
    // The tag is present when the end of the read was clipped at base calling.
    // The value of the tag is a string consisting of any one or more of the following:
    // A - adaptor clipped
    // Q - quality clipped
    // Z - "three zeros" clipped
    final  public static String  CLIPPING_TAG_NAME = "tm";

    final public static String FLOW_MATRIX_TAG_NAME = "tp";
    final public static String FLOW_MATRIX_T0_TAG_NAME = "t0";
    final public static String FLOW_MATRiX_OLD_TAG_KR = "kr";
    final public static String FLOW_MATRiX_OLD_TAG_TI = "ti";
    public static final String FLOW_MATRiX_OLD_TAG_KH = "kh";
    public static final String FLOW_MATRiX_OLD_TAG_KF = "kf";
    public static final String FLOW_MATRiX_OLD_TAG_KD = "kd";

    final public static String MAX_CLASS_READ_GROUP_TAG = "mc";

    final public static String[] HARD_CLIPPED_TAGS = new String[] {
            FLOW_MATRIX_TAG_NAME, FLOW_MATRiX_OLD_TAG_KR, FLOW_MATRiX_OLD_TAG_TI,
            FLOW_MATRiX_OLD_TAG_KH, FLOW_MATRiX_OLD_TAG_KF, FLOW_MATRiX_OLD_TAG_KD,
            FLOW_MATRIX_T0_TAG_NAME
    };

    /**
     * The sam record from which this flow based read originated
     */
    private SAMRecord samRecord;

    /**
     * The read's sequence, always in forward direction
     */
    private byte[] forwardSequence;

    /**
     * The flow key for the read - i.e. lengths of hmers in an flow order.
     *
     * For example, assuming a flow order of TGCA, and a forward sequence of GGAAT, the key will be 0,2,0,2,1
     */
    private int[] key;

    /**
     * the maping of key elements to their origin locations in the sequence. Entry n contains the offset in the sequence
     * where the hmer described by this key element starts.
     */
    private int [] flow2base;

    /**
     * The maximal length of an hmer that can be encoded (normally in the 10-12 range)
     */
    private int maxHmer;

    /**
     * The value to fill the flow matrix with. Normally 0.001
     */
    private double totalMinErrorProb;
    private double perHmerMinErrorProb;

    /**
     * The order in which flow key in encoded (See decription for key field). Flow order may be wrapped if a longer one
     * needed.
     */
    private byte[] flowOrder;

    /**
     * The probability matrix for this read. [n][m] position represents that probablity that an hmer of n length will be
     * present at the m key position. Therefore, the first dimention is in the maxHmer order, where the second dimension
     * is length(key).
     */
    private double[][] flowMatrix;

    /**
     * The validity status of the key. Certain operations may produce undefined/errornous results. This is signaled by
     * the read being marked with a validKey == false
     */
    private boolean validKey = true;

    /**
     * The direction of this read. After being red, the direction will also swing to be REFERENCE
     */
    private Direction direction = Direction.SYNTHESIS;

    /**
     * Was base clipping applied to this read? (normally to trim a read to the span of a haplotype)
     */
    private boolean baseClipped = false;

    /**
     * If applyBaseClipping was called, the left trimming that was actually applied to the read
     */
    private int trimLeftBase = 0 ;

    /**
     * If applyBaseClipping was called, the right trimming that was actually applied to the read
     */
    private int trimRightBase = 0 ;

    /**
     * The flow based argument collection under which this read was created
     */
    private final FlowBasedArgumentCollection fbargs;

    /**
     * This allows tools to reduce/enlarge the lower limit of read size for clipping operations.
     */
    static private int minimalReadLength = MINIMAL_READ_LENGTH;


    /**
     * Various methods for storing arrays of pre-computed qualities to be used in {@link FlowBasedHMMEngine}.
     */
    public byte[] getReadInsQuals() {
        return readInsQuals;
    }
    public void setReadInsQuals(byte[] readInsQuals) {
        this.readInsQuals = readInsQuals;
    }
    public byte[] getReadDelQuals() {
        return readDelQuals;
    }
    public void setReadDelQuals(byte[] readDelQuals) {
        this.readDelQuals = readDelQuals;
    }
    public byte[] getOverallGCP() {
        return overallGCP;
    }
    public void setOverallGCP(byte[] overallGCP) {
        this.overallGCP = overallGCP;
    }
    private byte[] readInsQuals;
    private byte[] readDelQuals;
    private byte[] overallGCP;

    public enum Direction {
        REFERENCE, SYNTHESIS
    }

    /**

     * Constructor from GATKRead. flow order, hmer and arguments
     * @param read GATK read
     * @param flowOrder flow order string (one cycle)
     * @param maxHmer maximal hmer to keep in the flow matrix
     * @param fbargs arguments that control resolution etc. of the flow matrix
     */
    public FlowBasedRead(final GATKRead read, final String flowOrder, final int maxHmer, final FlowBasedArgumentCollection fbargs) {
        this(read.convertToSAMRecord(null), flowOrder, maxHmer, fbargs);
    }

    /**
     * Same as above but constructs from SAMRecord
     * @param samRecord record from SAM file
     * @param flowOrder flow order (single cycle)
     * @param maxHmer maximal hmer to keep in the flow matrix
     * @param fbargs arguments that control resoltion of the flow matrix
     */
    public FlowBasedRead(final SAMRecord samRecord, final String flowOrder, final int maxHmer, final FlowBasedArgumentCollection fbargs) {
        super(samRecord);
        Utils.nonNull(fbargs);
        Utils.validate(FlowBasedReadUtils.hasFlowTags(samRecord), "FlowBasedRead can only be used on flow reads. failing read: " + samRecord);
        this.fbargs = fbargs;
        this.maxHmer = maxHmer;
        this.samRecord = samRecord;
        forwardSequence = getForwardSequence();

        // read flow matrix in. note that below code contains accomodates for old formats
        if ( samRecord.hasAttribute(FLOW_MATRIX_TAG_NAME) ) {
            perHmerMinErrorProb = fbargs.fillingValue;
            totalMinErrorProb = perHmerMinErrorProb;
            // this path is the production path. A flow read should contain a FLOW_MATRIX_TAG_NAME tag
            readFlowMatrix(flowOrder);

        } else {

            // NOTE: this path is vestigial and deals with old formats of the matrix
            if ( samRecord.hasAttribute(FLOW_MATRiX_OLD_TAG_KR) ) {
                readVestigialFlowMatrixFromKR(flowOrder);
            } else if ( samRecord.hasAttribute(FLOW_MATRiX_OLD_TAG_TI) ) {
                readVestigialFlowMatrixFromTI(flowOrder);
            } else {
                throw new GATKException("read missing flow matrix attribute: " + FLOW_MATRIX_TAG_NAME);
            }
        }

        implementMatrixMods(FlowBasedReadUtils.getFlowMatrixModsInstructions(fbargs.flowMatrixMods, maxHmer));

        //Spread boundary flow probabilities when the read is unclipped
        //in this case the value of the hmer is uncertain
        if ( !fbargs.keepBoundaryFlows ) {
            if (samRecord.getReadUnmappedFlag() || CigarUtils.countClippedBases(samRecord.getCigar(), Tail.LEFT, CigarOperator.HARD_CLIP) == 0) {
                spreadFlowLengthProbsAcrossCountsAtFlow(findFirstNonZero(key));
            }
            if (samRecord.getReadUnmappedFlag() || CigarUtils.countClippedBases(samRecord.getCigar(), Tail.RIGHT, CigarOperator.HARD_CLIP) == 0) {
                spreadFlowLengthProbsAcrossCountsAtFlow(findLastNonZero(key));
            }
        }

        if ( logger.isDebugEnabled() ) {
            logger.debug("cons: name: " + samRecord.getReadName()
                    + " len: " + samRecord.getReadLength()
                    + " loc: " + samRecord.getStart() + "-" + samRecord.getEnd()
                    + " rev: " + isReverseStrand()
                    + " cigar:" + samRecord.getCigarString());
            logger.debug("     bases: " + new String(samRecord.getReadBases()));
            logger.debug("       key: " + FlowBasedKeyCodec.keyAsString(key));
        }

    }

    //since the last unclipped flow is uncertain (we give high probabilities to
    //also hmers higher than the called hmer)
    private void spreadFlowLengthProbsAcrossCountsAtFlow(final int flowToSpread) {
        if (flowToSpread<0) //boundary case when all the key is zero
            return;

        final int call = key[flowToSpread];
        if (call==0){
            throw new IllegalStateException("Boundary key value should not be zero for the spreading");
        }

        final int numberToFill = maxHmer - call+1;
        double total = 0;
        for (int i = call; i < maxHmer+1; i++)
            total += flowMatrix[i][flowToSpread];
        final double fillProb = Math.max(total / numberToFill, perHmerMinErrorProb);
        for (int i = call; i < maxHmer+1; i++){
            flowMatrix[i][flowToSpread] = fillProb;
        }
    }



    // This is the code for parsing the current/production BAM format (with TP tag)
    private void readFlowMatrix(final String _flowOrder) {

        // generate key (base to flow space)
        setDirection(Direction.REFERENCE);  // base is always in reference/alignment direction

        key = FlowBasedKeyCodec.baseArrayToKey(samRecord.getReadBases(), _flowOrder);
        flow2base = FlowBasedKeyCodec.getKeyToBase(key);
        flowOrder = FlowBasedKeyCodec.getFlowToBase(_flowOrder, key.length);

        if (perHmerMinErrorProb == 0){
            // estimate the lowest possible error probability and fill the matrix with it
            // it is not recommended to do this in the M2 or HaplotypeCaller, since obviously the error are not
            // equally distributed over the homopolymers
            totalMinErrorProb = estimateMinErrorProb();
            perHmerMinErrorProb = totalMinErrorProb/maxHmer;
        }

        // initialize matrix. fill first line, copy subsequent lines from first
        flowMatrix = new double[maxHmer+1][key.length];
        Arrays.fill(flowMatrix[0], perHmerMinErrorProb);
        for (int i = 1 ; i < maxHmer+1; i++) {
            System.arraycopy(flowMatrix[0], 0, flowMatrix[i], 0, key.length);
        }

        // access qual, convert to flow representation
        final byte[]      quals = samRecord.getBaseQualities();
        final byte[]      tp = samRecord.getSignedByteArrayAttribute(FLOW_MATRIX_TAG_NAME);

        boolean specialTreatmentForZeroCalls = false;
        final byte[]      t0 = SAMUtils.fastqToPhred(samRecord.getStringAttribute(FLOW_MATRIX_T0_TAG_NAME));
        final double[]     t0probs = new double[quals.length];
        if ((t0!=null) && fbargs.useT0Tag){
            specialTreatmentForZeroCalls = true;

            if (t0.length!=tp.length){
                throw new GATKException("Illegal read len(t0)!=len(qual): " + samRecord.getReadName());
            }
        }

        final double[]    probs = new double[quals.length];
        for ( int i = 0 ; i < quals.length ; i++ ) {
            final double q = quals[i];
            final double p = Math.pow(10, -q/10);
            probs[i] = p;
            if (specialTreatmentForZeroCalls){
                final double qq = t0[i];
                t0probs[i] = Math.pow(10,-qq/10);
            }
        }

        // apply key and qual/tp to matrix
        int     qualOfs = 0; //converts between base -> flow
        for ( int i = 0 ; i < key.length ; i++ ) {
            final int        run = key[i];
            if (run > 0) {
                parseSingleHmer(probs, tp, i, run, qualOfs);
            }

            if ((run == 0)&&(specialTreatmentForZeroCalls)){
                parseZeroQuals(t0probs, i, qualOfs);
            }

            double totalErrorProb = 0;

            for (int k=0; k < maxHmer; k++ ){
                totalErrorProb += flowMatrix[k][i];
            }
            final double callProb = Math.max(MINIMAL_CALL_PROB, 1-totalErrorProb);
            // the probability in the recalibration is not divided by two for hmers of length 1
            flowMatrix[Math.min(run, maxHmer)][i] = callProb;
            qualOfs+=run;
        }
        applyFilteringFlowMatrix();
    }


    //convert qualities from the single hmer to a column in a flow matrix
    private void parseSingleHmer(final double[] probs, final byte[] tp, final int flowIdx,
                                 final int flowCall, final int qualOfs){
        for (int i = qualOfs ; i < qualOfs+flowCall; i++) {
            if (tp[i]!=0) {
                final int loc = Math.max(Math.min(flowCall+tp[i], maxHmer),0);
                if (flowMatrix[loc][flowIdx] == perHmerMinErrorProb) {
                    flowMatrix[loc][flowIdx] = probs[i];
                } else {
                    flowMatrix[loc][flowIdx]    += probs[i];
                }
            }
        }
    }

    // convert qualities from the t0 tag to the probabilities of 1->0 error.
    // This function deals with t0 tag that encodes the probability of 1->0 error
    // in this case there is no nucleotide to place the error probability on, so we
    // place it on the neighboring bases and choose the **lower** error probability between the
    // neighbors (that's how T0 encoding works). The error is placed only on the 1-mer error assuming
    // that 2->0 errors are negligibly rare.
    private void parseZeroQuals(final double[] probs, final int flowIdx, final int qualOfs){
        if ((qualOfs == 0) | (qualOfs==probs.length)){ // do not report zero error probability on the edge of the read
            return;
        }

        double prob0 = Math.min(probs[qualOfs-1],probs[qualOfs]);
        // cases where t0 actually does not report anything. This awkward situation comes because
        // if fbargs.fillingValue is zero, the empty cell probability is maximalQuality / getMaxHmer()
        // while if not, the empty cell probability is fbargs.fillingValue.
        if (prob0 <= totalMinErrorProb * 3) {
            prob0 = 0;
        }
        flowMatrix[1][flowIdx] = Math.max(flowMatrix[1][flowIdx], prob0);
    }

    public String getFlowOrder() {
        return new String(Arrays.copyOfRange(flowOrder, 0, Math.min(fbargs.flowOrderCycleLength,flowOrder.length)));
    }


    public int getMaxHmer() {
        return maxHmer;
    }

    public int getNFlows() {
        return key.length;
    }

    public Direction getDirection(){
        return direction;
    }


    private byte[] getForwardSequence(){
        if (!isReverseStrand()) {
            return samRecord.getReadBases();
        } else {
            final byte[] result = new byte[samRecord.getReadBases().length];
            System.arraycopy(samRecord.getReadBases(), 0, result, 0, result.length);
            SequenceUtil.reverseComplement(result);
            return result;
        }
    }


    private int[] getAttributeAsIntArray(final String attributeName, final boolean isSigned) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);
        final Object attributeValue = this.samRecord.getAttribute(attributeName);

        if (attributeValue == null) {
            return null;
        } else if (attributeValue instanceof byte[]) {
            final byte[] byteArrayAttributeValue = (byte[]) attributeValue;
            final int[] ret = new int[byteArrayAttributeValue.length];
            for (int i = 0; i < ret.length; i++)
                if (!isSigned) ret[i] = byteArrayAttributeValue[i]&0xff;
                else ret[i]=byteArrayAttributeValue[i]; //converting signed byte to unsigned
            return Arrays.copyOf(ret, ret.length);
        } else if ((attributeValue instanceof int[])) {
            final int[] ret = (int[]) attributeValue;
            return Arrays.copyOf(ret, ret.length);
        } else if  (attributeValue instanceof short[]) {
            final short [] shortArrayAttributeValue = (short[]) attributeValue;
            final int[] ret = new int[shortArrayAttributeValue.length];
            for (int i = 0 ; i < shortArrayAttributeValue.length; i++ )
                ret[i] = shortArrayAttributeValue[i];
            return Arrays.copyOf(ret, ret.length);
        }else {
            throw new GATKException.ReadAttributeTypeMismatch(attributeName, "integer array");
        }
    }

    public boolean isValid() {
        return validKey;
    }

    /**
     * get a specific cell from the flow matrix. Each cell contains the probability
     * for an hmer of the given length to appear the given position in the flow key
     *
     * @param flow - position in the flow key (index into key[])
     * @param hmer - length of the hmer
     * @return
     */
    public double getProb(final int flow, final int hmer) {
        double prob = flowMatrix[hmer < maxHmer ? hmer : maxHmer][flow];
        return (prob <= 1) ? prob : 1;
    }

    /*
    * Legacy function from the time when the error probability were in flow space and when the read was clipped we had to
    * translate the clipping in the base space to the clipping in the flow space. Now does nothing (isBaseFormat() is true)
    * but is kept for R&D cases.
     */
    public void applyAlignment(){
        if ((getDirection() == Direction.SYNTHESIS) && ( isReverseStrand() )) {
            flipMatrix();
            ArrayUtils.reverse(key);
            flow2base = FlowBasedKeyCodec.getKeyToBase(key);
            SequenceUtil.reverseComplement(flowOrder);

        }

        final boolean isBase = isBaseFormat();
        final int[] basePair = {0, 0};
        final int[] clipLeftPair = !isBase ? findLeftClippingFromCigar() : basePair;
        final int[] clipRightPair = !isBase ? findRightClippingFromCigar() : basePair;
        final int clipLeft = clipLeftPair[0];
        final int leftHmerClip = clipLeftPair[1];
        final int clipRight = clipRightPair[0];
        final int rightHmerClip = clipRightPair[1];

        applyClipping(clipLeft, leftHmerClip, clipRight, rightHmerClip, false);

        setDirection(Direction.REFERENCE);

    }

    private boolean isBaseFormat() {
       return samRecord.hasAttribute(FLOW_MATRiX_OLD_TAG_TI) || samRecord.hasAttribute(FLOW_MATRIX_TAG_NAME);
    }

    private void fillFlowMatrix(final int [] kh, final int [] kf,
                                final double [] kdProbs ) {
        for ( int i = 0 ; i < kh.length; i++ ) {
            // normal matrix filling
            final int     pos = kf[i];
            final int     hmer = kh[i] & 0xff;
            if (hmer > maxHmer){
                flowMatrix[maxHmer][pos] = Math.max(flowMatrix[maxHmer][pos], kdProbs[i]);
            } else {
                flowMatrix[hmer][pos] = Math.max(flowMatrix[hmer][pos], kdProbs[i]);
            }
        }

    }

    // execute the matrix modifications
    private void implementMatrixMods(final int[] flowMatrixModsInstructions) {

        if ( flowMatrixModsInstructions != null ) {
            for (int hmer = 0; hmer < flowMatrixModsInstructions.length; hmer++) {
                final int hmer2 = flowMatrixModsInstructions[hmer];
                if (hmer2 != 0) {
                    for (int pos = 0; pos < flowMatrix[0].length; pos++) {

                        if (flowMatrix[hmer][pos] > flowMatrix[hmer2][pos]) {
                            flowMatrix[hmer2][pos] = flowMatrix[hmer][pos];
                        }

                        // if we are copying bacwards, zero out source
                        if (hmer > hmer2)
                            flowMatrix[hmer][pos] = 0;
                    }
                }
            }
        }
    }


    private void flipMatrix() {
        for ( int i = 0 ; i < flowMatrix.length; i++) {
            ArrayUtils.reverse(flowMatrix[i]);
        }
    }

    private static int findFirstNonZero(final int[] array){
        int result = -1;
        for (int i = 0 ; i < array.length; i++){
            if (array[i]!=0) {
                result = i;
                break;
            }
        }
        return result;
    }

    private static int findLastNonZero(final int[] array){
        int result = -1;
        for (int i = array.length-1 ; i >= 0; i--){
            if (array[i]!=0) {
                result = i;
                break;
            }
        }
        return result;
    }

    private static void shiftColumnUp(final double[][] matrix, final int colnum, final int shift) {
        for (int i = 0; i < matrix.length - shift; i ++ ) {
            matrix[i][colnum] = matrix[i+shift][colnum];
        }
        for (int i = matrix.length - shift; i < matrix.length; i ++ ) {
            matrix[i][colnum] = 0;
        }

    }

    public void setDirection(final Direction dir ) {
        direction = dir;
    }

    //trims base-spaced reads. Usually not needed, but kept for completeness
    public void applyBaseClipping(final int clipLeftBase, final int clipRightBase, boolean spread){
        final int[] clipLeftPair = findLeftClipping(clipLeftBase);
        final int[] clipRightPair = findRightClipping(clipRightBase);
        final int clipLeft = clipLeftPair[0];
        final int leftHmerClip = clipLeftPair[1];
        final int clipRight = clipRightPair[0];
        final int rightHmerClip = clipRightPair[1];
        if (getLength() - clipLeftBase - clipRightBase < minimalReadLength) {
            baseClipped = true;
            validKey = false;
            trimLeftBase = clipLeftBase;
            trimRightBase = clipRightBase;
        } else {
            applyClipping(clipLeft, leftHmerClip, clipRight, rightHmerClip, spread);
            baseClipped = true;
            validKey = true;
            trimLeftBase = clipLeftBase;
            trimRightBase = clipRightBase;
        }
    }

    private void applyClipping(int clipLeft, final int leftHmerClip, int clipRight, final int rightHmerClip, final boolean spread){
        if ((clipLeft < 0) || (clipRight < 0)  || (clipLeft >= getKeyLength() ) || ( clipRight >= getKeyLength())) {
            throw new IllegalStateException(String.format("Weird read clip calculated: left/right/keyLength %d/%d/%d", clipLeft, clipRight, getKeyLength()));
        }

        if ((leftHmerClip < 0) || (rightHmerClip < 0)  || (leftHmerClip >= (getMaxHmer() + 2) ) || ( rightHmerClip >= (getMaxHmer() + 2))) {
            throw new IllegalStateException(String.format("Weird read clip calculated: left/right/maxHmer+2 %d/%d/%d", leftHmerClip, rightHmerClip, getMaxHmer() + 2));
        }

        final int originalLength = key.length;

        key[clipLeft]-=leftHmerClip;
        boolean shiftLeft = true;
        while (key[clipLeft] == 0) {
            clipLeft += 1 ;
            shiftLeft = false;
        }

        key[key.length - clipRight-1] -= rightHmerClip;
        boolean shiftRight = true;
        while (key[originalLength - 1- clipRight] == 0) {
            clipRight += 1 ;
            shiftRight = false;
        }

        key = Arrays.copyOfRange(key, clipLeft, originalLength - clipRight);
        flow2base = FlowBasedKeyCodec.getKeyToBase(key);
        flowOrder = Arrays.copyOfRange(flowOrder, clipLeft, originalLength - clipRight);

        final double [][] newFlowMatrix = new double[flowMatrix.length][originalLength - clipLeft - clipRight] ;
        for ( int i = 0 ; i < newFlowMatrix.length; i++) {
            newFlowMatrix[i] = Arrays.copyOfRange(flowMatrix[i], clipLeft, originalLength - clipRight);
        }

        flowMatrix = newFlowMatrix;
        if (shiftLeft) {
            shiftColumnUp(flowMatrix, 0, leftHmerClip);
        }

        if (shiftRight) {
            shiftColumnUp(flowMatrix, flowMatrix[0].length-1, rightHmerClip);
        }

        //Spread boundary flow probabilities for the boundary hmers of the read
        //in this case the value of the genome hmer is uncertain
        if ( spread ) {
            spreadFlowLengthProbsAcrossCountsAtFlow(findFirstNonZero(key));
            spreadFlowLengthProbsAcrossCountsAtFlow(findLastNonZero(key));
        }
    }

    private int[] findLeftClippingFromCigar() {
        final List<CigarElement> cigar = getCigarElements();
        final int[] result = new int[2];
        if (cigar.size() == 0 ) {
            return result;
        }

        final CigarElement start = cigar.get(0);
        if (start.getOperator() != CigarOperator.H) {
            return result;
        }

        final int basesClipped = start.getLength();
        return findLeftClipping(basesClipped);
    }


    private int[] findLeftClipping(final int basesClipped){
        return FlowBasedReadUtils.findLeftClipping(basesClipped, flow2base, key);
    }

    private int[] findRightClippingFromCigar() {
        final List<CigarElement> cigar = getCigarElements();
        final int[] result = new int[2];
        if (cigar.size() == 0 ) {
            result[0] = 0;
            result[1] = 0;
            return result;
        }

        final CigarElement end = cigar.get(cigar.size()-1);
        if (end.getOperator() != CigarOperator.H) {
            result[0] = 0;
            result[1] = 0;
            return result;
        }

        final int basesClipped = end.getLength();

        return findRightClipping(basesClipped);
    }

    private int[] findRightClipping(final int basesClipped) {
        final int[] rkey = new int[key.length];
        for (int i = 0 ; i < key.length; i++ ){
            rkey[i] = key[key.length-1-i];
        }
        final int[] rflow2base = FlowBasedKeyCodec.getKeyToBase(rkey);

        return FlowBasedReadUtils.findRightClipping(basesClipped, rflow2base, rkey);
    }


    public void writeKey(final FileWriter oos)
            throws IOException {
        for (int i = 0; i < key.length; i++)
            oos.write(key[i]+"\n");

    }

    /**
     * Flow matrix logger
     *
     * This is used exclusively for testing
     *
     * @param oos
     * @throws IOException
     */
    protected void writeMatrix(final OutputStreamWriter oos)
            throws IOException {
        final DecimalFormat formatter = new DecimalFormat("0.000000", DecimalFormatSymbols.getInstance(Locale.ENGLISH));

        final byte[]      bases = samRecord.getReadBases();
        int         basesOfs = 0;
        final byte[]      quals = samRecord.getBaseQualities();
        final byte[]      ti = samRecord.hasAttribute(FLOW_MATRiX_OLD_TAG_TI) ? samRecord.getByteArrayAttribute(FLOW_MATRiX_OLD_TAG_TI) : (new byte[key.length*3]);
        if ( isReverseStrand() )
        {
            ArrayUtils.reverse(quals);
            ArrayUtils.reverse(ti);
        }

        for ( int col = 0 ; col < key.length ; col++ ) {
            oos.write("C,R,F,B,Bi,Q,ti\n");
            final byte base = (key[col] != 0) ? (basesOfs < bases.length ? bases[basesOfs] : (byte)'?') : (byte)'.';
            final String bi = (key[col] != 0) ? Integer.toString(basesOfs) : ".";
            final String q = (key[col] != 0) ? Integer.toString(quals[basesOfs]) : ".";
            final String Ti = (key[col] != 0) ? Integer.toString(ti[basesOfs]) : ".";
            for (int row = 0; row < flowMatrix.length; row++) {
                final String s = formatter.format(flowMatrix[row][col]);
                oos.write(""
                        + col + ","
                        + row + ","
                        + key[col] + ","
                        + (char)base + ","
                        + bi + ","
                        + q + ","
                        + Ti + ","
                        + (isReverseStrand() ? "r" : ".") + " "
                        + s + "\n");
            }
            if ( key[col] != 0 )
                basesOfs +=  key[col];
            oos.write("\n");
        }

    }


    public byte [] getFlowOrderArray() {
        return flowOrder;
    }

    public int getKeyLength() {
        return key.length;
    }

    public int[] getKey() {
        return key;
    }

    /**
     * Number of total bases that the flow based key generates
     * @return number of bases
     */
    public int totalKeyBases()  {
        int sum = 0 ;
        for (int i = 0 ; i < key.length; i++){
            sum += key[i];
        }
        return sum;
    }

    public int seqLength(){
        return forwardSequence.length;
    }
    public boolean isBaseClipped() {
        return baseClipped;
    }

    public int getTrimmedStart() {
        return trimLeftBase + getStart();
    }
    public int getTrimmedEnd() {
        return getEnd() - trimRightBase;
    }



    //functions that take care of simulating base format
    //they perform modifications on the flow matrix that are defined in applyFilteringFlowMatrix
    //this function was only applied when we tested what is the necessary information to be reported in the flow matrix
    private void applyFilteringFlowMatrix(){

        if (fbargs.disallowLargerProbs) {
            removeLargeProbs();
        }

        if (fbargs.removeLongerThanOneIndels) {
            removeLongIndels( key );
        }

        if (fbargs.removeOneToZeroProbs) {
            removeOneToZeroProbs(key);
        }

        if ((fbargs.lumpProbs)) {
            lumpProbs();
        }
        clipProbs();

        if (fbargs.symmetricIndels) {
            smoothIndels(key);
        }
        if (fbargs.onlyInsOrDel) {
            reportInsOrDel(key);
        }

        if ((fbargs.retainMaxNProbs)){
            reportMaxNProbsHmer(key);
        }

    }

    /**
     * clip probability values in a way that probability between perHmerMinErrorProb and 3*perHmerMinErrorProb
     * is automatically clipped to perHmerMinErrorProb. This is done to avoid issues with rounding of
     * small error probabilities in the basecalling
     */
    private void clipProbs() {
        double probabilityThreshold = perHmerMinErrorProb*3;
        for ( int i = 0 ; i < getMaxHmer(); i++ ) {
            for ( int j =0; j < getNFlows(); j++) {
                if ((flowMatrix[i][j] <= probabilityThreshold) &&
                        (key[j]!=i)) {
                    flowMatrix[i][j] = perHmerMinErrorProb;
                }
            }
        }
    }

    /**
     * remove probabilities larger than 1
     */
    private void removeLargeProbs(){
        for (int i = 0; i < getNFlows(); i++){
            for (int j = 0 ; j < getMaxHmer()+1; j++) {
                if (flowMatrix[j][i] > 1) {
                    flowMatrix[j][i] = 1;
                }
            }
        }
    }

    /**
     * This is vestigial and applies only to old formats
     * @param key_kh
     */
    private void removeLongIndels(final int [] key_kh ){
        for ( int i = 0 ; i < getNFlows(); i++ ) {
            for (int j = 0; j < getMaxHmer()+1; j++){
                if (Math.abs(j-key_kh[i])>1){
                    flowMatrix[j][i] = perHmerMinErrorProb;
                }
            }
        }
    }

    /**
     * This is vestigial and applies only to old formats
     * @param key_kh
     */
    private void removeOneToZeroProbs(final int [] key_kh) {
        for (int i = 0 ; i < getNFlows(); i++){
            if (key_kh[i] == 0){
                for (int j = 1; j < getMaxHmer()+1; j++){
                    flowMatrix[j][i] = perHmerMinErrorProb;
                }
            }
        }
    }


    /**
     * Quantize probability values according to fbargs.probabilityQuantization and fbargs.probabilityScalingFactor
     * @param kd_probs
     */

    private void quantizeProbs(final int [] kd_probs ) {
        final int nQuants = fbargs.probabilityQuantization;
        final double bin_size = 6*fbargs.probabilityScalingFactor/(float)nQuants;
        for ( int i = 0 ; i < kd_probs.length; i++) {
            if (kd_probs[i] <=0)
                continue;
            else {
                kd_probs[i] = (byte)(bin_size * (byte)(kd_probs[i]/bin_size)+1);
            }
        }
    }

    /**
     * Smooth out probabilities by averaging with neighbours
     * @param kr
     */
    private void smoothIndels(final int [] kr ) {
        for ( int i = 0 ; i < kr.length; i++ ){
            final int idx = kr[i];
            if (( idx > 1 ) && ( idx < maxHmer) ) {
                final double prob = (flowMatrix[idx - 1][i] + flowMatrix[idx + 1][i]) / 2;
                flowMatrix[idx - 1][i] = prob;
                flowMatrix[idx + 1][i] = prob;
            }
        }
    }

    /**
     * Only report probability of insertions or of deletions, never both
     * @param kr
     */
    private void reportInsOrDel(final int [] kr ) {
        for ( int i = 0 ; i < kr.length; i++ ){
            final int idx = kr[i];
            if (( idx > 1 ) && ( idx < maxHmer) ) {
                if ((flowMatrix[idx-1][i] > perHmerMinErrorProb) && (flowMatrix[idx+1][i] > perHmerMinErrorProb)) {
                    final int fixCell = flowMatrix[idx-1][i] > flowMatrix[idx+1][i] ? idx+1 : idx-1;
                    flowMatrix[fixCell][i] = perHmerMinErrorProb;
                }
            }
        }
    }

    /**
     * Combine all probabilities of insertions together and report them as probabilities of 1-mer insertion
     * Combine all probabilities of deletions together and report them as probabilities of 1-mer deletion
     */
    private void lumpProbs() {

        for (int i = 0; i < getMaxHmer(); i++) {
            for (int j = 0 ; j < getNFlows(); j ++ ) {
                final int fkey = key[j];
                if (flowMatrix[i][j]<=perHmerMinErrorProb) {
                    continue;
                } else {
                    if ( (i - fkey) < -1 ){
                        flowMatrix[fkey-1][j]+=flowMatrix[i][j];
                        flowMatrix[i][j] = perHmerMinErrorProb;
                    } else if ((i-fkey) > 1) {
                        flowMatrix[fkey+1][j]+=flowMatrix[i][j];
                        flowMatrix[i][j] = perHmerMinErrorProb;
                    }

                }

            }
        }

    }

    /*
    Given full vector of error probabilities retain only the probabilities that can be reported in the base format
    (N+1/2 highest error probabilities)
     */
    private void reportMaxNProbsHmer(final int [] key) {
        final double [] tmpContainer = new double[maxHmer];
        for (int i = 0 ; i < key.length;i++){

            for (int j = 0 ; j < tmpContainer.length; j++) {
                tmpContainer[j] = flowMatrix[j][i];
            }
            final int k = (key[i]+1)/2;
            final double kth_highest = findKthLargest(tmpContainer, k+1);
            for (int j = 0 ; j < maxHmer; j++)
                if (flowMatrix[j][i] < kth_highest)
                    flowMatrix[j][i] = perHmerMinErrorProb;
        }

    }


    private static double findKthLargest(final double[] nums, final int k) {
        final PriorityQueue<Double> q = new PriorityQueue<Double>(k);
        for(final double i: nums){
            q.offer(i);

            if(q.size()>k){
                q.poll();
            }
        }

        return q.peek();
    }

    private double[] phredToProb(final int [] kq) {
        final double [] result = new double[kq.length];
        for (int i = 0 ; i < kq.length; i++ ) {
            //disallow probabilities below fillingValue
            result[i] = Math.max(Math.pow(10, ((double)-kq[i])/fbargs.probabilityScalingFactor), perHmerMinErrorProb);
        }
        return result;
    }


    public static void setMinimalReadLength(int minimalReadLength) {
        FlowBasedRead.minimalReadLength = minimalReadLength;
    }

    // convert qualities and ti tag to flow matrix
    private void readVestigialFlowMatrixFromTI(final String _flowOrder) {

        vestigialOneShotLogger.warn("Vestigial read format detected: " + samRecord);

        // generate key (base to flow space)
        setDirection(Direction.REFERENCE);  // base is always in reference/alignment direction
        key = FlowBasedKeyCodec.baseArrayToKey(samRecord.getReadBases(), _flowOrder);
        flow2base = FlowBasedKeyCodec.getKeyToBase(key);
        flowOrder = FlowBasedKeyCodec.getFlowToBase(_flowOrder, key.length);

        // initialize matrix
        flowMatrix = new double[maxHmer+1][key.length];
        for (int i = 0 ; i < maxHmer+1; i++) {
            for (int j = 0 ; j < key.length; j++ ){
                flowMatrix[i][j] = perHmerMinErrorProb;
            }
        }

        // access qual, convert to flow representation
        final byte[]      quals = samRecord.getBaseQualities();
        final byte[]      ti = samRecord.getByteArrayAttribute(FLOW_MATRiX_OLD_TAG_TI);
        final double[]    probs = new double[quals.length];
        for ( int i = 0 ; i < quals.length ; i++ ) {
            final double q = quals[i];
            final double p = QualityUtils.qualToErrorProb(q);
            probs[i] = p*2;
        }

        // apply key and qual/ti to matrix
        int     qualOfs = 0;
        for ( int i = 0 ; i < key.length ; i++ ) {
            final int        run = key[i];

            // the probability is not divided by two for hmers of length 1
            if ( run == 1 ) {
                probs[qualOfs] = probs[qualOfs]/2;
            }

            //filling the probability for the called hmer (not reported by the quals
            if ( run <= maxHmer ) {
                flowMatrix[run][i] = (run > 0) ? (1 - probs[qualOfs]) : 1;
                //require a prob. at least 0.1
                flowMatrix[run][i] = Math.max(MINIMAL_CALL_PROB, flowMatrix[run][i]);

            }

            if ( run != 0 ) {
                if ( quals[qualOfs] != 40 ) {
                    final int     run1 = (ti[qualOfs] == 0) ? (run - 1) : (run + 1);
                    if (( run1 <= maxHmer ) && (run <= maxHmer)){
                        flowMatrix[run1][i] = probs[qualOfs] / flowMatrix[run][i];
                    }
                    if (run <= maxHmer) {
                        flowMatrix[run][i] /= flowMatrix[run][i]; // for comparison to the flow space - probabilities are normalized by the key's probability
                    }
                }
                qualOfs += run;
            }

        }

        //this is just for tests of all kinds of
        applyFilteringFlowMatrix();
    }

    // code for reading BAM format where the flow matrix is stored in sparse representation in kr,kf,kh and kd tags
    // used for development of the new basecalling, but not in production code
    private void readVestigialFlowMatrixFromKR(final String _flowOrder) {

        vestigialOneShotLogger.warn("Vestigial read format detected: " + samRecord);

        key = getAttributeAsIntArray(FLOW_MATRiX_OLD_TAG_KR, true);

        // creates a translation from flow # to base #
        flow2base = FlowBasedKeyCodec.getKeyToBase(key);

        // create a translation from
        flowOrder = FlowBasedKeyCodec.getFlowToBase(_flowOrder, key.length);

        flowMatrix = new double[maxHmer+1][key.length];
        for (int i = 0 ; i < maxHmer+1; i++) {
            for (int j = 0 ; j < key.length; j++ ){
                flowMatrix[i][j] = perHmerMinErrorProb;
            }
        }

        int [] kh = getAttributeAsIntArray(FLOW_MATRiX_OLD_TAG_KH, true);
        int [] kf = getAttributeAsIntArray(FLOW_MATRiX_OLD_TAG_KF, false);
        int [] kd = getAttributeAsIntArray(FLOW_MATRiX_OLD_TAG_KD, true);

        final int [] key_kh = key;
        final int [] key_kf = new int[key.length];
        for ( int i = 0 ; i < key_kf.length ; i++)
            key_kf[i] = i;
        final int [] key_kd = new int[key.length];

        kh = ArrayUtils.addAll(kh, key_kh);
        kf = ArrayUtils.addAll(kf, key_kf);
        kd = ArrayUtils.addAll(kd, key_kd);

        quantizeProbs(kd);

        final double [] kdProbs = phredToProb(kd);
        fillFlowMatrix( kh, kf, kdProbs);
        applyFilteringFlowMatrix();
    }
    //Finds the quality that is being set when the probability of error is very low
    private double estimateMinErrorProb(){
        final byte[] quals = samRecord.getBaseQualities();
        double maxQual = 0;

        for (int i = 0; i < quals.length; i++) {
            if (quals[i] > maxQual){
                maxQual = quals[i];
            }
        }

        // in the very rare case when there is no qualities
        if (maxQual==0){
            maxQual = 40;
        }
        return Math.pow(10, -maxQual / 10);
    }
}


