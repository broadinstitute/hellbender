package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotationData;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.*;
import java.util.stream.IntStream;


/**
 * Root Mean Square of the mapping quality of reads across all samples.
 *
 * <p>This annotation provides an estimation of the overall mapping quality of reads supporting a variant call, averaged over all samples in a cohort.</p>
 *
 * <h3>Statistical notes</h3>
 * <p>The root mean square is equivalent to the mean of the mapping qualities plus the standard deviation of the mapping qualities.</p>
 *
 * <h3>Related annotations</h3>
 * <ul>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_MappingQualityRankSumTest.php">MappingQualityRankSumTest</a></b> compares the mapping quality of reads supporting the REF and ALT alleles.</li>
 * </ul>
 *
 */
public final class RMSMappingQuality extends InfoFieldAnnotation implements StandardAnnotation, ReducibleAnnotation {

    @Override
    public String getRawKeyName() { return GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_KEY;}

    /**
     * Generate the raw data necessary to calculate the annotation. Raw data is the final endpoint for gVCFs.
     */
    @Override
    public Map<String, Object> annotateRawData(final ReferenceContext ref,
                                               final VariantContext vc,
                                               final ReadLikelihoods<Allele> likelihoods){
        Utils.nonNull(vc);
        if (likelihoods == null || likelihoods.readCount() == 0) {
            return Collections.emptyMap();
        }

        final Map<String, Object> annotations = new HashMap<>();
        final ReducibleAnnotationData<Number> myData = new ReducibleAnnotationData<>(null);
        calculateRawData(vc, likelihoods, myData);
        final String annotationString = formatedValue((double) myData.getAttributeMap().get(Allele.NO_CALL));
        annotations.put(getRawKeyName(), annotationString);
        return annotations;
    }

    @Override
    @SuppressWarnings({"unchecked", "rawtypes"})//FIXME
    public void calculateRawData(final VariantContext vc,
                                 final ReadLikelihoods<Allele> likelihoods,
                                 final ReducibleAnnotationData rawAnnotations){
        //put this as a double, like GATK3.5
        final double squareSum = IntStream.range(0, likelihoods.numberOfSamples()).boxed()
                .flatMap(s -> likelihoods.sampleReads(s).stream())
                .map(GATKRead::getMappingQuality)
                .filter(mq -> mq != QualityUtils.MAPPING_QUALITY_UNAVAILABLE)
                .mapToDouble(mq -> mq * mq).sum();

        rawAnnotations.putAttribute(Allele.NO_CALL, squareSum);
    }

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final ReadLikelihoods<Allele> likelihoods) {
        return annotateRawData(ref, vc, likelihoods);
    }

    @VisibleForTesting
    static String formatedValue(double rms) {
        return String.format("%.2f", rms);
    }

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(VCFConstants.RMS_MAPPING_QUALITY_KEY, getRawKeyName());
    }


    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(VCFStandardHeaderLines.getInfoLine(getKeyNames().get(0)), GATKVCFHeaderLines.getInfoLine(getRawKeyName()));
    }

    /**
     * converts {@link GATKVCFConstants#RAW_RMS_MAPPING_QUALITY_KEY} into  {@link VCFConstants#RMS_MAPPING_QUALITY_KEY}  annotation if present
     * @param vc which potentially contains rawMQ
     * @return if vc contained {@link GATKVCFConstants#RAW_RMS_MAPPING_QUALITY_KEY} it will be replaced with {@link VCFConstants#RMS_MAPPING_QUALITY_KEY}
     * otherwise return the original vc
     */
    public VariantContext finalizeRawMQ(final VariantContext vc) {
        final String rawMQdata = vc.getAttributeAsString(getRawKeyName(), null);
        if (rawMQdata == null) {
            return vc;
        } else {
            final double squareSum = parseRawDataString(rawMQdata);
            final String finalizedRMSMAppingQuality = makeFinalizedAnnotationString(vc, squareSum);
            final VariantContext result = new VariantContextBuilder(vc)
                    .rmAttribute(getRawKeyName())
                    .attribute(getKeyNames().get(0), finalizedRMSMAppingQuality)
                    .make();
            return result;
        }
    }

    private static double parseRawDataString(String rawDataString) {
        try {
            final double squareSum = Double.parseDouble(rawDataString.split(",")[0]);
            return squareSum;
        } catch (final NumberFormatException e){
            throw new UserException.BadInput("malformed " + GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_KEY +" annotation: " + rawDataString);
        }

    }


    private static String makeFinalizedAnnotationString(final VariantContext vc, final double squareSum) {
        final int numOfReads = getNumOfReads(vc);
        return formatedValue(Math.sqrt(squareSum /numOfReads));
    }

    /**
     *
     * @return the number of reads at the vc position (-1 if all read data is null)
     */
    private static int getNumOfReads(final VariantContext vc) {
        //don't use the full depth because we don't calculate MQ for reference blocks
        int numOfReads = vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, -1);
        if(vc.hasGenotypes()) {
            for(final Genotype gt : vc.getGenotypes()) {
                if(gt.isHomRef()) {
                    //site-level DP contribution will come from MIN_DP for gVCF-called reference variants or DP for BP resolution
                    if (gt.hasExtendedAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY)) {
                        numOfReads -= Integer.parseInt(gt.getExtendedAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY).toString());
                    }
                    else if (gt.hasDP()) {
                        numOfReads -= gt.getDP();
                    }
                }
            }
        }
        if (numOfReads <= 0){
            throw new UserException.BadInput("Cannot calculate Root Mean Square Mapping Quality if there are 0 or less reads." +
                                                     "\nNumber of reads recorded as :" +numOfReads +
                                                     "\nIn VariantContext: "+ vc.toStringDecodeGenotypes());
        }
        return numOfReads;
    }

}
