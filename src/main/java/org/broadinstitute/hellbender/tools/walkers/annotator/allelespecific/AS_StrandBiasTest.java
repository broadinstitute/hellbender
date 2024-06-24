package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.AnnotationUtils;
import org.broadinstitute.hellbender.tools.walkers.annotator.StrandBiasTest;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;

/**
 * Allele-specific implementation of strand bias annotations
 */
public abstract class AS_StrandBiasTest extends StrandBiasTest implements ReducibleAnnotation, AlleleSpecificAnnotation {
    private final static Logger logger = LogManager.getLogger(AS_StrandBiasTest.class);
    public static final String SPLIT_DELIM = "\\|"; //String.split takes a regex, so we need to escape the pipe
    public static final String PRINT_DELIM = "|";
    public static final String REDUCED_DELIM = ",";
    public static final int MIN_COUNT = 2;
    public static final double MIN_PVALUE = 1.0E-320;
    public static final int FORWARD = 0;
    public static final int REVERSE = 1;

    @Override
    public String getPrimaryRawKey() { return GATKVCFConstants.AS_SB_TABLE_KEY; }

    /**
     * @return true if annotation has secondary raw keys
     */
    @Override
    public boolean hasSecondaryRawKeys() {
        return false;
    }

    /**
     * Get additional raw key strings that are not the primary key
     *
     * @return may be null
     */
    @Override
    public List<String> getSecondaryRawKeys() {
        return null;
    }

    @Override
    public String getEmptyRawValue() {
        return "0,0";
    }

    /**
     * Method which determines how the Strand Bias read direction allele data must be combined into a final annotation
     * Must be overridden by client methods.
     *
     * @param combinedData
     * @return
     */
    protected abstract Map<Allele,Double> calculateReducedData(final AlleleSpecificAnnotationData<List<Integer>> combinedData );


    /**
     * Uses the likelihoods map to generate a 2x2 strand contingency table by counting the total read support for each
     * allele in either the forward or reverse direction.
     *
     * @param ref the reference context for this annotation
     * @param vc the variant context to annotate
     * @param likelihoods likelihoods indexed by sample, allele, and read within sample
     * @return
     */
    @Override
    public Map<String, Object> annotateRawData(final ReferenceContext ref,
                                               final VariantContext vc,
                                               final AlleleLikelihoods<GATKRead, Allele> likelihoods ) {

        //for allele-specific annotations we only call from HC and we only use likelihoods
        if ( likelihoods == null) {
            return Collections.emptyMap();
        }
        return StrandBiasUtils.computeSBAnnotation(vc, likelihoods, getPrimaryRawKey());
    }

    protected String makeReducedAnnotationString(VariantContext vc, Map<Allele,Double> perAltsStrandCounts) {
        String annotationString = "";
        for (Allele a : vc.getAlternateAlleles()) {
            if (!annotationString.isEmpty()) {
                annotationString += REDUCED_DELIM;
            }
            if (!perAltsStrandCounts.containsKey(a)) {
                logger.warn("ERROR: VC allele not found in annotation alleles -- maybe there was trimming?");
            } else {
                annotationString += perAltsStrandCounts.get(a) == null ? VCFConstants.MISSING_VALUE_v4 : String.format("%.3f", perAltsStrandCounts.get(a));
            }
        }
        return annotationString;
    }

    /**
     * Method which combines the per allele contingency tables from the underlying variant contexts by totaling
     * supported values for both forward and reverse data and outputting it as a new contingency table.
     *
     * @param vcAlleles
     * @param annotationList
     * @return
     */
    @Override
    @SuppressWarnings({"unchecked", "rawtypes"})//FIXME generics here blow up
    public Map<String, Object> combineRawData(final List<Allele> vcAlleles, final List<ReducibleAnnotationData<?>>  annotationList) {
        //VC already contains merged alleles from ReferenceConfidenceVariantContextMerger
        ReducibleAnnotationData combinedData = new AlleleSpecificAnnotationData(vcAlleles, null);

        for (final ReducibleAnnotationData currentValue : annotationList) {
            parseRawDataString(currentValue);
            StrandBiasUtils.combineAttributeMap(currentValue, combinedData);
        }
        final String annotationString = StrandBiasUtils.makeRawAnnotationString(vcAlleles, combinedData.getAttributeMap());
        return Collections.singletonMap(getPrimaryRawKey(), annotationString);
    }

    /**
     * Parses the raw data stings of combined contingency matrix data and calls client methods calculateReducedData(myData)
     * implementation to generate double digest of provided allele information which is stored in '|' delineated lists.
     *
     * @param vc -- contains the final set of alleles, possibly subset by GenotypeGVCFs
     * @param originalVC -- used to get all the alleles for all gVCFs
     * @return
     */
    @Override
    public  Map<String, Object> finalizeRawData(final VariantContext vc, final VariantContext originalVC) {
        if (!vc.hasAttribute(getPrimaryRawKey())) {
            return new HashMap<>();
        }
        String rawContingencyTableData = vc.getAttributeAsString(getPrimaryRawKey(),null);
        if (rawContingencyTableData == null) {
            return new HashMap<>();
        }
        AlleleSpecificAnnotationData<List<Integer>> myData = new AlleleSpecificAnnotationData<>(originalVC.getAlleles(), rawContingencyTableData);
        parseRawDataString(myData);

        Map<Allele, Double> perAltRankSumResults = calculateReducedData(myData);

        String annotationString = makeReducedAnnotationString(vc, perAltRankSumResults);
        String rawAnnotationsString = StrandBiasUtils.makeRawAnnotationString(vc.getAlleles(), myData.getAttributeMap());
        Map<String, Object> returnMap = new HashMap<>();
        returnMap.put(getKeyNames().get(0), annotationString);
        returnMap.put(getPrimaryRawKey(), rawAnnotationsString);  //this is in case raw annotations are requested
        return returnMap;
    }

    protected void parseRawDataString(ReducibleAnnotationData<List<Integer>> myData) {
        List<String> values = AnnotationUtils.getAlleleLengthListOfStringFromRawData(myData.getRawData());
        if (values.size() != myData.getAlleles().size()) {
            throw new IllegalStateException("Number of alleles and number of allele-specific entries do not match.  " +
                    "Allele-specific annotations should have an entry for each allele including the reference.");
        }

        Map<Allele, List<Integer>> perAlleleValues = new HashMap<>();
        for (int i = 0; i < values.size(); i++) {
            List<Integer> perAlleleList = new ArrayList<>();
            String[] rawListEntriesAsStringVector = values.get(i).split(",");
            //Read counts will only ever be integers
            for (String s : rawListEntriesAsStringVector) {
                if (!s.isEmpty()) {
                    perAlleleList.add(Integer.parseInt(s.trim()));
                }
            }
            perAlleleValues.put(myData.getAlleles().get(i), perAlleleList);
        }

        myData.setAttributeMap(perAlleleValues);
    }

    @Override
    //Allele-specific annotations cannot be called from walkers other than HaplotypeCaller
    protected Map<String, Object> calculateAnnotationFromGTfield(final GenotypesContext genotypes){
        return Collections.emptyMap();
    }

    public static String rawValueAsString(int[][] table) {
        return table[0][0]+","+table[0][1]+ PRINT_DELIM +table[1][0]+","+table[1][1];
    }
}
