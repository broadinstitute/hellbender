package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.StringUtils;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.*;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.*;

public final class AnnotationUtils {

    public static final String ALLELE_SPECIFIC_ANNOTATION_KEY_PREFIX = "AS_";
    public static final String ALLELE_SPECIFIC_RAW_DELIM = "|";
    public static final String ALLELE_SPECIFIC_REDUCED_DELIM = ",";
    public static final String ALLELE_SPECIFIC_SPLIT_REGEX = "\\|"; //String.split takes a regex, so we need to escape the pipe
    public static final String BRACKET_REGEX = "\\[|\\]";
    public static final String LIST_DELIMITER = ",";
    public static final String MISSING_VALUE = ".";

    private AnnotationUtils(){}

    /**
     * Helper function to parse the list into the annotation string
     * @param valueList the ArrayList returned from StrandBiasBySample.annotate()
     * @return the array used by the per-sample Strand Bias annotation
     */
    public static String encodeValueList(final List<Double> valueList, final String precisionFormat ) {
        List<String> outputList = new ArrayList<>();
        for (Double d : valueList) {
            outputList.add(String.format(precisionFormat, d));
        }
        return StringUtils.join(outputList, LIST_DELIMITER);
    }

    /**
     * Helper function to convert a List of Strings to a comma-separated String
     * @param stringList the ArrayList with String data
     * @return a comma-separated String
     */
    public static String encodeStringList( final List<String> stringList) {
        return StringUtils.join(stringList, LIST_DELIMITER);
    }

    /**
     * Helper function to convert a List of Strings to a @{value ALLELE_SPECIFIC_RAW_DELIM)-separated String, as for raw annotations
     * @param somethingList the ArrayList with String data
     * @return a delimited String
     */
    public static String encodeAnyASListWithRawDelim(final List<?> somethingList) {
        return StringUtils.join(somethingList, ALLELE_SPECIFIC_RAW_DELIM).replaceAll(BRACKET_REGEX, "");  //Who actually wants brackets at the ends of their string?  Who???
    }

    /**
     * Helper method to split a "raw" annotation string delimited with {@value ALLELE_SPECIFIC_RAW_DELIM}
     * @param somethingList a String, possibly read from a VCF
     * @return a List of Strings
     */
    public static List<String> decodeAnyASListWithRawDelim(final String somethingList) {
        return Arrays.asList(StringUtils.splitByWholeSeparatorPreserveAllTokens(somethingList.replaceAll(BRACKET_REGEX, ""), ALLELE_SPECIFIC_RAW_DELIM));
    }

    /**
     * Helper function to convert a comma-separated String (such as a vc.getAttrbute().toString() output) to a List of Strings
     * @param somethingList the allele-specific annotations string; may have brackets
     * @return a list of allele-specific annotation entries
     */
    public static List<String> decodeAnyASList( final String somethingList) {
        return Arrays.asList(StringUtils.splitByWholeSeparatorPreserveAllTokens(somethingList.replaceAll(BRACKET_REGEX, ""), ALLELE_SPECIFIC_REDUCED_DELIM));
    }

    /**
     * Helper function to determine if an annotation is allele-specific
     * @param annotation the annotation to be tested
     * @return true if the annotation is expected to have values per-allele
     */
    public static boolean isAlleleSpecific(final VariantAnnotation annotation) {
        return annotation instanceof AlleleSpecificAnnotation;
    }

    public static boolean isAlleleSpecificGatkKey(final String annotationKey) {
        final VCFInfoHeaderLine header = GATKVCFHeaderLines.getInfoLine(annotationKey);
        return header.getCountType().equals(VCFHeaderLineCount.A) ||
                header.getCountType().equals(VCFHeaderLineCount.R) ||
                annotationKey.startsWith(ALLELE_SPECIFIC_ANNOTATION_KEY_PREFIX);
    }

    /**
     * Handles all the Java and htsjdk parsing shenanigans from getAttributeAsString
     * @param rawDataString may have surrounding brackets, with raw delimiter
     * @return
     */
    public static List<String> getAlleleLengthListOfStringFromRawData(String rawDataString) {
        if (rawDataString == null) {
            return Collections.emptyList();
        }
        if (rawDataString.startsWith("[")) {
            rawDataString = rawDataString.substring(1, rawDataString.length() - 1).replaceAll("\\s", "");
        }
        return Arrays.asList(rawDataString.split(ALLELE_SPECIFIC_SPLIT_REGEX, -1)); //-1 to keep empty data
    }

    /**
     * Handles all the Java and htsjdk parsing shenanigans from getAttributeAsString
     * @param dataString may have surrounding brackets, with reduced delimieter
     * @return
     */
    public static List<String> getAlleleLengthListOfString(String dataString) {
        if (dataString == null) {
            return Collections.emptyList();
        }
        if (dataString.startsWith("[")) {
            dataString = dataString.substring(1, dataString.length() - 1).replaceAll("\\s", "");
        }
        return Arrays.asList(dataString.split(ALLELE_SPECIFIC_REDUCED_DELIM, -1)); //-1 to keep empty data
    }

    static public String generateMissingDataWarning(final VariantContext vc, final Genotype g, final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        final StringBuilder outString = new StringBuilder("Annotation will not be calculated at position " + vc.getContig() + ":" + vc.getStart() +
                " and possibly subsequent");
        if (!g.isCalled()) {
            outString.append("; genotype for sample " + g.getSampleName() + " is not called");
        }
        if (likelihoods == null) {
            outString.append("; alleleLikelihoodMap is null");
        }
        return outString.toString();
    }
}
