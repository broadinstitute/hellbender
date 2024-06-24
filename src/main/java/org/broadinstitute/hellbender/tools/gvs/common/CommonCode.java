package org.broadinstitute.hellbender.tools.gvs.common;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

//TODO rename this or get rid of it. a place holder for now
public class CommonCode {
    public static String CALL_RATE = "CALL_RATE";
    public static String INVARIANT = "INVARIANT";
    public static String HWE = "HWE";

    public static String HWE_FILTER = "HWE";
    public static String CALL_RATE_FILTER = "CALL_RATE";
    public static String INVARIANT_FILTER = "INVARIANT";

    public static final int EXCESS_ALLELES_THRESHOLD = 6;

    /**
     * If the alleles are "missing", -1 will be returned as the index
     * @param variant
     * @return
     */
    public static List<Integer> getGTAlleleIndexes(final VariantContext variant) {
        IndexedAlleleList<Allele> alleleList = new IndexedAlleleList<>(variant.getAlleles());
        ArrayList<Integer> allele_indices = new ArrayList<Integer>();

        for (Allele allele : variant.getGenotype(0).getAlleles()) {
            allele_indices.add(alleleList.indexOfAllele(allele));
        }
        return allele_indices;
    }

    public static String getGTString(final VariantContext variant) {
        List<Integer> allele_indices = getGTAlleleIndexes(variant);

        // As of VS-910, we aren't going to consider any number of alleles to necessarily be a UserException
        // This code used to throw an error if allele_indices.size() != 2

        List<String> gsStrings = allele_indices.stream().map(index -> index == -1 ? "." : index.toString()).collect(Collectors.toList());
        String separator = variant.getGenotype(0).isPhased() ? VCFConstants.PHASED : VCFConstants.UNPHASED;
        return StringUtils.join(gsStrings, separator);
    }

    public enum OutputType {
        TSV,
        ORC,
        AVRO,
        PARQUET,
        TSV2,
        BQ,
        NONE
    }
}
