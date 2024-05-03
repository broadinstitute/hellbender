package org.broadinstitute.hellbender.tools.walkers.annotator.flow;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import org.broadinstitute.hellbender.utils.read.FlowBasedKeyCodec;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Base class for flow based annotations
 *
 * Some flow based annotations depend on the results from other annotations, regardless
 * if they were called for by user arguments. To overcome this, this class contains all shared
 * code to compute flow based annotations.
 *
 * Each (specific) annotation is implemented as a subclass of this class. It then invokes the
 * annotation calculation methods contained here (the shared code) to compute its prerequisite
 * and itself.
 *
 * State between such calls is kept in a LocalContext, a local class. Its is there were annotations
 * are accumulated as well.
 */
public abstract class FlowAnnotatorBase implements InfoFieldAnnotation {
    private final static Logger logger = LogManager.getLogger(FlowAnnotatorBase.class);
    protected final OneShotLogger flowMissingOneShotLogger = new OneShotLogger(FlowAnnotatorBase.class);


    // additional constants
    protected static final String   C_INSERT = "ins";
    protected static final String   C_DELETE = "del";
    protected static final String   C_NA = "NA";
    protected static final String   C_CSS_CS = "cycle-skip";
    protected static final String   C_CSS_PCS = "possible-cycle-skip";
    protected static final String   C_CSS_NS = "non-skip";

    protected static final String   C_SNP = "snp";
    protected static final String   C_NON_H_MER = "non-h-indel";
    protected static final String   C_H_MER = "h-indel";


    protected static final int      MOTIF_SIZE = 5;
    protected static final int      GC_CONTENT_SIZE = 10;
    protected static final int      BASE_TYPE_COUNT = 4;

    private List<String>            flowOrder;


    protected class LocalContext {
        ReferenceContext ref;
        AlleleLikelihoods<GATKRead, Allele> likelihoods;
        String      flowOrder;

        List<String> indel;
        List<Integer> indelLength;
        List<Integer> hmerIndelLength;
        List<String>  leftMotif;
        List<String>  rightMotif;

        Map<String, Object> attributes = new LinkedHashMap<>();

        boolean     generateAnnotation;

        protected LocalContext(final ReferenceContext ref,
                               final VariantContext vc,
                               final AlleleLikelihoods<GATKRead, Allele> likelihoods,
                               final boolean needsRef) {
            Utils.nonNull(vc);
            if ( needsRef ) {
                Utils.validate(ref == null || ref.hasBackingDataSource(), "-R (reference) argument must be provided");
            }

            // some annotators share results
            this.ref = ref;
            this.likelihoods = likelihoods;

            // annotation will be generated by default
            this.generateAnnotation = !needsRef || (ref != null);
        }

        protected Map<String, Object> asAttributes() {

            if ( !generateAnnotation ) {
                return Collections.emptyMap();
            } else {
                return attributes.entrySet().stream()
                        .filter(x -> getKeyNames().contains(x.getKey()))
                        .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
            }
        }
    }

    /*
    This function establishes the flow order to be used for manipulating reads in flow space.

    The most natural source for the flow order is are the reads themselves. Alas reads will
    not always be sourced from a bam file with a flow order. In these cases, we can either get it from a
    --flow-order parameter (VariantAnnotator tool) or the default input source (bam)
     */
    private String establishFlowOrder(final LocalContext localContext, final AlleleLikelihoods<GATKRead, Allele> likelihoods) {

        // extract from a read
        if ( (likelihoods != null) && (likelihoods.numberOfSamples() > 0) ) {
            final List<GATKRead>  reads = likelihoods.sampleEvidence(0);
            if ( reads.size() > 0 ) {
                GATKRead        read = reads.get(0);
                if ( read instanceof FlowBasedRead ) {
                    return ((FlowBasedRead)read).getFlowOrder();
                } else if ( flowOrder != null )  {
                    establishReadGroupFlowOrder(localContext, read.getReadGroup());
                }
            }
        }

        // use global
        return establishReadGroupFlowOrder(localContext, null);
    }

    /*
        the flow order might be different for each read group.
        provided flow order can be a list of [group:]flowOrder separated by a comma
        no group: means all/rest
     */
    private String establishReadGroupFlowOrder(final LocalContext localContext, final String readGroup) {

        // find flow order for the readGroup
        if ( flowOrder != null ) {
            for (String elem : flowOrder) {
                final String toks[] = elem.split(":");
                if (toks.length == 1) {
                    return toks[0];
                } else if (toks[0].equals(readGroup)) {
                    return toks[1];
                }
            }
        }

        // if here, no flow order was found. may we use a default?
        if ( isActualFlowOrderRequired() ) {
            localContext.generateAnnotation = false;
            flowMissingOneShotLogger.warn(this.getClass().getSimpleName() + " annotation will not be calculated, no '" + StandardArgumentDefinitions.FLOW_ORDER_FOR_ANNOTATIONS + "' argument provided");
        }

        return FlowBasedRead.DEFAULT_FLOW_ORDER;
    }

    protected boolean isActualFlowOrderRequired() {
        return false;
    }

    // "indel_classify" and "indel_length"
    protected void indelClassify(final VariantContext vc, final LocalContext localContext) {

        final List<String>      indelClassify = new LinkedList<>();
        final List<Integer>     indelLength = new LinkedList<>();
        final int               refLength = vc.getReference().length();
        for ( Allele a : vc.getAlleles() ) {
            if ( !a.isReference() ) {
                indelClassify.add(refLength == a.length() ? C_NA : (refLength < a.length() ? C_INSERT : C_DELETE));
                if ( !isSpecial(a) && (a.length() != refLength) ) {
                    indelLength.add(Math.abs(refLength - a.length()));
                } else {
                    indelLength.add(null);
                }
            }
        }
        localContext.attributes.put(GATKVCFConstants.FLOW_INDEL_CLASSIFY, localContext.indel = indelClassify);
        localContext.attributes.put(GATKVCFConstants.FLOW_INDEL_LENGTH, localContext.indelLength = indelLength);
    }

    // "indel_classify" and "indel_length"
    protected void variantType(final VariantContext vc, final LocalContext localContext) {
        List<Allele> alleles = vc.getAlternateAlleles();
        boolean isSnp = true;
        for (int i = 0;  i < alleles.size(); i++){
            if (isSpecial(alleles.get(i))){
                continue;
            }
            if (!localContext.indel.get(i).equals(C_NA)){
                isSnp=false;
            }
        }
        if (isSnp){
            localContext.attributes.put(GATKVCFConstants.FLOW_VARIANT_TYPE, C_SNP);
            return;
        }

        boolean isHmer = true;
        for (int i = 0;  i < alleles.size(); i++){
            if (isSpecial(alleles.get(i))){
                continue;
            }
            if ((localContext.hmerIndelLength==null) || (localContext.hmerIndelLength.get(i)==null) || (localContext.hmerIndelLength.get(i)==0)){
                isHmer=false;
            }
        }

        if (isHmer){
            localContext.attributes.put(GATKVCFConstants.FLOW_VARIANT_TYPE, C_H_MER);
            return;
        }

        localContext.attributes.put(GATKVCFConstants.FLOW_VARIANT_TYPE, C_NON_H_MER);
    }


    /*
    This function determines if the vc is an hmer indel. If so, it marks it as such
     */
    protected void isHmerIndel(final VariantContext vc, final LocalContext localContext) {

        // loop over all allels
        final List<Integer>     hmerIndelLength = new LinkedList<>();
        final List<String>      hmerIndelNuc = new LinkedList<>();
        final List<String>      rightMotif = new LinkedList<>();
        for ( Allele a : vc.getAlleles() ) {

            // skip reference
            if ( a.isReference() ) {
                continue;
            }

            // establish flow order
            if ( localContext.flowOrder == null ) {
                localContext.flowOrder = establishFlowOrder(localContext, localContext.likelihoods);}

            // assume no meaningful result
            hmerIndelLength.add(0);
            hmerIndelNuc.add(null);
            rightMotif.add(null);

            // access alleles
            final Allele      ref = vc.getReference();
            final Allele      alt = a;
            if ( isSpecial(a) )
                continue;;

            // get byte before and after
            final byte        before = getReferenceNucleotide(localContext, vc.getStart() - 1);
            final byte[]      after = getReferenceHmerPlus(localContext, vc.getEnd() + 1, MOTIF_SIZE);
            if (after.length==0){
                logger.warn("Failed to get hmer from reference, isHmerIndel and RightMotif annotations will not be calculated. " +
                        "Start index: " + vc.getEnd() + 1 + " Reference length: " + localContext.ref.getBases().length);
                return;
            }
            // build two haplotypes. add byte before and after
            final byte[]      refHap = buildHaplotype(before, ref.getBases(), after);
            final byte[]      altHap = buildHaplotype(before, alt.getBases(), after);

            // convert to flow space
            final int[]       refKey = FlowBasedKeyCodec.baseArrayToKey(refHap, localContext.flowOrder);
            final int[]       altKey = FlowBasedKeyCodec.baseArrayToKey(altHap, localContext.flowOrder);
            if ( refKey == null || altKey == null ) {
                throw new GATKException("failed to generate key from reference or alternate sequence");
            }

            // key must be the same length to begin with
            if ( refKey.length != altKey.length ) {
                continue;
            }

            // key must have only one difference, which should not be between a zero and something
            int     diffIndex = -1;
            int     refBasesCountUpInclHmer = 0;
            for ( int n = 0 ; n < refKey.length ; n++ ) {
                // count ref bases up to and including difference key
                if ( diffIndex < 0 ) {
                    refBasesCountUpInclHmer += refKey[n];
                }

                // is this the (one) difference key?
                if ( refKey[n] != altKey[n] ) {
                    if ( diffIndex >= 0 ) {
                        // break away
                        diffIndex = -1;
                        break;
                    } else {
                        diffIndex = n;
                    }
                }
            }

            // check if we've actually encountered a significant different key
            if ( diffIndex < 0 ) {
                continue;
            }
            if ( Math.max(refKey[diffIndex], altKey[diffIndex]) == 0 ) {
                continue;
            }

            // if here, we found the difference. replace last element of list
            final int             length = Math.max(refKey[diffIndex], altKey[diffIndex]);
            final byte            nuc = localContext.flowOrder.getBytes()[diffIndex % localContext.flowOrder.length()];
            hmerIndelLength.set(hmerIndelLength.size() - 1, length);
            hmerIndelNuc.set(hmerIndelNuc.size() - 1, Character.toString((char)nuc));

            // at this point, we can generate the right motif (for the hmer indel) as we already have the location
            // of the hmer-indel and the bases following it
            if ( a.length() != ref.length() ) {
                final String    motif = new String(Arrays.copyOfRange(refHap, refBasesCountUpInclHmer, Math.min(refHap.length, refBasesCountUpInclHmer + MOTIF_SIZE)));
                rightMotif.set(rightMotif.size() - 1, motif);
            }
        }

        // reflect back to attributs and context
        localContext.attributes.put(GATKVCFConstants.FLOW_HMER_INDEL_LENGTH, localContext.hmerIndelLength = hmerIndelLength);
        localContext.attributes.put(GATKVCFConstants.FLOW_HMER_INDEL_NUC, hmerIndelNuc);
        localContext.rightMotif = rightMotif;
    }

    private byte[] buildHaplotype(final byte before, final byte[] bases, final byte[] after) {

        final byte[]  hap = new byte[1 + bases.length + after.length];

        hap[0] = before;
        System.arraycopy(bases, 0, hap, 1, bases.length);
        System.arraycopy(after, 0, hap, 1 + bases.length, after.length);

        return hap;
    }

    protected void getLeftMotif(final VariantContext vc, final LocalContext localContext) {

        final int         refLength = vc.getReference().length();
        final List<String> leftMotif = new LinkedList<>();

        for ( Allele a : vc.getAlleles() ) {
            if ( a.isReference() ) {
                continue;
            }

            String  motif = getRefMotif(localContext, vc.getStart() - MOTIF_SIZE, MOTIF_SIZE);
            if (motif.length() != MOTIF_SIZE){
                logger.warn("Failed to get motif from reference, getLeftMotif annotation will not be calculated. " +
                        "Start index: " + vc.getStart() + " Reference length: " + localContext.ref.getBases().length);
                return;
            }
            if ( a.length() != refLength ) {
                motif = motif.substring(1) + vc.getReference().getBaseString().substring(0, 1);
            }
            leftMotif.add(motif);
        }

        localContext.attributes.put(GATKVCFConstants.FLOW_LEFT_MOTIF, localContext.leftMotif = leftMotif);
    }

    protected void getRightMotif(final VariantContext vc, final LocalContext localContext) {

        final int         refLength = vc.getReference().length();
        String      motif = getRefMotif(localContext, vc.getStart() + refLength, MOTIF_SIZE);
        if (motif.length() != MOTIF_SIZE){
            logger.warn("Failed to get motif from reference, getRightMotif annotation will not be calculated. " +
                    "Start index: " + vc.getStart() + refLength + " Reference length: " + localContext.ref.getBases().length);
            return;
        }
        // fill empty entries (non indel alelles)
        for ( int i = 0 ; i < localContext.rightMotif.size() ; i++ ) {
            if ( localContext.rightMotif.get(i) == null ) {
                localContext.rightMotif.set(i, motif);
            }
        }

        localContext.attributes.put(GATKVCFConstants.FLOW_RIGHT_MOTIF, localContext.rightMotif);
    }

    protected void gcContent(final VariantContext vc, final LocalContext localContext) {

        final int         begin = vc.getStart() - (GC_CONTENT_SIZE / 2);
        final String      seq = getRefMotif(localContext, begin + 1, GC_CONTENT_SIZE);
        if ( seq.length() != GC_CONTENT_SIZE ) {
            logger.warn("gcContent will not be calculated at position " + vc.getContig() + ":" + vc.getStart() +
                    " too close to the edge of the reference");
            return;
        }
        int         gcCount = 0;
        for ( byte b : seq.getBytes() ) {
            if ( b == 'G' || b == 'C' ) {
                gcCount++;
            }
        }
        localContext.attributes.put(GATKVCFConstants.FLOW_GC_CONTENT, (float)gcCount / seq.length());
    }

    protected void cycleSkip(final VariantContext vc, final LocalContext localContext) {

        // establish flow order
        if ( localContext.flowOrder == null ) {
            localContext.flowOrder = establishFlowOrder(localContext, localContext.likelihoods);
        }

        // loop over alleles
        final List<String>  css = new LinkedList<>();
        final int           refLength = vc.getReference().length();
        for ( Allele a : vc.getAlleles() ) {
            if ( a.isReference() ) {
                continue;
            }

            // meaningful only for non indels
            if ( isSpecial(a) || (a.length() != refLength) ) {
                css.add(C_NA);
            } else {

                // access alleles
                final Allele      ref = vc.getReference();
                final Allele      alt = a;

                // convert to flow space
                final int         i = css.size();   // always working on the last
                final int[]       refKey = FlowBasedKeyCodec.baseArrayToKey((localContext.leftMotif.get(i) + ref.getBaseString() + localContext.rightMotif.get(i)).getBytes(), localContext.flowOrder);
                final int[]       altKey = FlowBasedKeyCodec.baseArrayToKey((localContext.leftMotif.get(i) + (!isSpecial(alt) ? alt.getBaseString() : "") + localContext.rightMotif.get(i)).getBytes(), localContext.flowOrder);

                // assign initial css
                String            cssValue = (refKey.length != altKey.length) ? C_CSS_CS : C_CSS_NS;

                // if same length (NS) then see if it is possible-cycle-skip
                if ( cssValue == C_CSS_NS ) {
                    for ( int n = 0 ; n < refKey.length ; n++ ) {
                        if ( (refKey[n] == 0) ^ (altKey[n] == 0) ) {
                            cssValue = C_CSS_PCS;
                            break;
                        }
                    }
                }

                css.add(cssValue);
            }
        }

        localContext.attributes.put(GATKVCFConstants.FLOW_CYCLESKIP_STATUS, css);
    }

    // get a single nucleotid from reference
    private byte getReferenceNucleotide(final LocalContext localContext, final int start) {
        final int         index = start - localContext.ref.getWindow().getStart();
        final byte[]      bases = localContext.ref.getBases();
        Utils.validIndex(index, bases.length); // do not catch, if here the location is outside of the reference, there is a problem!
        return bases[index];
    }

    // get an hmer from reference plus a number of additional bases
    private byte[] getReferenceHmerPlus(final LocalContext localContext, final int start, final int additional) {
        int               index = start - localContext.ref.getWindow().getStart();
        final byte[]      bases = localContext.ref.getBases();
        try {
            Utils.validIndex(index, bases.length);
        } catch (IllegalArgumentException e) {
            logger.warn("Failed to get hmer from reference. Start index: " + index + " Reference length: " + bases.length);
            return new byte[0];
        }

        // get hmer
        final StringBuilder sb = new StringBuilder();
        final byte          base0 = bases[index++];
        sb.append((char)base0);
        for ( ; index < bases.length && bases[index] == base0 ; index++ ) {
            sb.append((char) bases[index]);
        }

        // get additional
        for ( int n = 0 ; n < additional && index < bases.length ; n++, index++ ) {
            sb.append((char) bases[index]);
        }

        return sb.toString().getBytes();
    }
    // get motif from reference
    private String getRefMotif(final LocalContext localContext, final int start, final int length) {
        final byte[]      bases = localContext.ref.getBases();
        final int         startIndex = start - localContext.ref.getWindow().getStart();
        final int         endIndex = startIndex + length;
        try {
            Utils.validIndex(startIndex, bases.length);
            Utils.validIndex(endIndex-1, bases.length);
        } catch (IllegalArgumentException e) {
            return "";
        }
        return new String(Arrays.copyOfRange(bases, startIndex, endIndex));
    }

    public void setFlowOrder(final List<String> flowOrder) {
        this.flowOrder = flowOrder;
    }

    private boolean isSpecial(Allele a) {
        return a.equals(Allele.SPAN_DEL) || a.equals(Allele.NON_REF_ALLELE);
    }
}
