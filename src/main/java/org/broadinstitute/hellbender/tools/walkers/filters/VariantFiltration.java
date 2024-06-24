package org.broadinstitute.hellbender.tools.walkers.filters;

import com.google.common.collect.Sets;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.VariantContextUtils.JexlVCMatchExp;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.AlleleFilterUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;
import java.util.stream.Collectors;

import static java.util.Collections.singleton;
import static org.broadinstitute.hellbender.utils.variant.GATKVCFConstants.AS_FILTER_STATUS_KEY;


/**
 * Filter variant calls based on INFO and/or FORMAT annotations
 *
 * <p>
 * This tool is designed for hard-filtering variant calls based on certain criteria. Records are hard-filtered by
 * changing the value in the FILTER field to something other than PASS. Filtered records will be preserved in the output
 * unless their removal is requested in the command line. </p>
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li>A VCF of variant calls to filter.</li>
 *     <li>One or more filtering expressions and corresponding filter names.</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <p>
 * A filtered VCF in which passing variants are annotated as PASS and failing variants are annotated with the name(s) of
 * the filter(s) they failed.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 *   gatk VariantFiltration \
 *   -R reference.fasta \
 *   -V input.vcf.gz \
 *   -O output.vcf.gz \
 *   --filter-name "my_filter1" \
 *   --filter-expression "AB < 0.2" \
 *   --filter-name "my_filter2" \
 *   --filter-expression "MQ0 > 50"
 * </pre>
 *
 * <h3>Note</h3>
 * <p>
 * Composing filtering expressions can range from very simple to extremely complicated depending on what you're
 * trying to do.
 * <p>
 * Compound expressions (ones that specify multiple conditions connected by &&, AND, ||, or OR, and reference
 * multiple attributes) require special consideration. By default, variants that are missing one or more of the
 * attributes referenced in a compound expression are treated as PASS for the entire expression, even if the variant
 * would satisfy the filter criteria for another part of the expression. This can lead to unexpected results if any
 * of the attributes referenced in a compound expression are present for some variants, but missing for others.
 * <p>
 * It is strongly recommended that such expressions be provided as individual arguments, each referencing a
 * single attribute and specifying a single criteria. This ensures that all of the individual expression are
 * applied to each variant, even if a given variant is missing values for some of the expression conditions.
 * <p>
 * As an example, multiple individual expressions provided like this:
 * <pre>
 *   gatk VariantFiltration \
 *   -R reference.fasta \
 *   -V input.vcf.gz \
 *   -O output.vcf.gz \
 *   --filter-name "my_filter1" \
 *   --filter-expression "AB < 0.2" \
 *   --filter-name "my_filter2" \
 *   --filter-expression "MQ0 > 50"
 * </pre>
 *
 * are preferable to a single compound expression such as this:
 *
 *  <pre>
 *    gatk VariantFiltration \
 *    -R reference.fasta \
 *    -V input.vcf.gz \
 *    -O output.vcf.gz \
 *    --filter-name "my_filter" \
 *    --filter-expression "AB < 0.2 || MQ0 > 50"
 *  </pre>
 * See this <a href="https://gatk.broadinstitute.org/hc/en-us/articles/360035891011-JEXL-filtering-expressions">article about using JEXL expressions</a>
 * for more information.
 */
@CommandLineProgramProperties(
        summary = "Filter variant calls based on INFO and/or FORMAT annotations.",
        oneLineSummary = "Filter variant calls based on INFO and/or FORMAT annotations",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
public final class VariantFiltration extends VariantWalker {

    public static final String FILTER_EXPRESSION_LONG_NAME = "filter-expression";
    public static final String FILTER_NAME_LONG_NAME = "filter-name";
    public static final String GENOTYPE_FILTER_EXPRESSION_LONG_NAME = "genotype-filter-expression";
    public static final String GENOTYPE_FILTER_NAME_LONG_NAME = "genotype-filter-name";
    public static final String CLUSTER_SIZE_LONG_NAME = "cluster-size";
    public static final String CLUSTER_WINDOW_SIZE_LONG_NAME = "cluster-window-size";
    public static final String MASK_EXTENSION_LONG_NAME = "mask-extension";
    public static final String MASK_NAME_LONG_NAME = "mask-name";
    public static final String MASK_DESCRIPTION_LONG_NAME = "mask-description";
    public static final String FILTER_NOT_IN_MASK_LONG_NAME = "filter-not-in-mask";
    public static final String MISSING_VAL_LONG_NAME = "missing-values-evaluate-as-failing";
    public static final String INVERT_LONG_NAME = "invert-filter-expression";
    public static final String INVERT_GT_LONG_NAME = "invert-genotype-filter-expression";
    public static final String NO_CALL_GTS_LONG_NAME = "set-filtered-genotype-to-no-call";
    public static final String ALLELE_SPECIFIC_LONG_NAME = "apply-allele-specific-filters";

    private static final String FILTER_DELIMITER = ";";

    /**
     * Any variant which overlaps entries from the provided mask file will be filtered. If the user wants logic to be reversed,
     * i.e. filter variants that do not overlap with provided mask, then argument --filter-not-in-mask can be used.
     * Note that it is up to the user to adapt the name of the mask to make it clear that the reverse logic was used
     * (e.g. if masking against Hapmap, use --mask-name=hapmap for the normal masking and --mask-name=not_hapmap for the reverse masking).
     */
    @Argument(fullName="mask", shortName="mask", doc="Input mask", optional=true)
    public FeatureInput<Feature> mask;

    @Argument(doc="File to which variants should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public GATKPath out = null;

    /**
     * VariantFiltration accepts any number of JEXL expressions (so you can have two named filters by using
     * --filter-name One --filter-expression "X < 1" --filter-name Two --filter-expression "X > 2").
     *
     * It is preferable to use multiple expressions, each specifying an individual filter criteria, to a single
     * compound expression that specifies multiple filter criteria.
     */
    @Argument(fullName=FILTER_EXPRESSION_LONG_NAME, shortName="filter", doc="One or more expressions used with INFO fields to filter", optional=true)
    public List<String> filterExpressions = new ArrayList<>();

    /**
     * This name is put in the FILTER field for variants that get filtered.  Note that there must be a 1-to-1 mapping between filter expressions and filter names.
     */
    @Argument(fullName=FILTER_NAME_LONG_NAME, doc="Names to use for the list of filters", optional=true)
    public List<String> filterNames = new ArrayList<>();

    /**
     * Similar to the INFO field based expressions, but used on the FORMAT (genotype) fields instead.
     * VariantFiltration will add the sample-level FT tag to the FORMAT field of filtered samples (this does not affect the record's FILTER tag).
     * One can filter normally based on most fields (e.g. "GQ < 5.0"), but the GT (genotype) field is an exception. We have put in convenience
     * methods so that one can now filter out hets ("isHet == 1"), refs ("isHomRef == 1"), or homs ("isHomVar == 1"). Also available are
     * expressions isCalled, isNoCall, isMixed, and isAvailable, in accordance with the methods of the Genotype object.
     *
     * It is preferable to use multiple expressions, each specifying an individual filter criteria, to a single compound expression
     * that specifies multiple filter criteria.
     */
    @Argument(fullName=GENOTYPE_FILTER_EXPRESSION_LONG_NAME, shortName="G-filter", doc="One or more expressions used with FORMAT (sample/genotype-level) fields to filter (see documentation guide for more info)", optional=true)
    public List<String> genotypeFilterExpressions = new ArrayList<>();

    /**
     * Similar to the INFO field based expressions, but used on the FORMAT (genotype) fields instead.
     */
    @Argument(fullName=GENOTYPE_FILTER_NAME_LONG_NAME, shortName="G-filter-name", doc="Names to use for the list of sample/genotype filters (must be a 1-to-1 mapping); this name is put in the FILTER field for variants that get filtered", optional=true)
    public List<String> genotypeFilterNames = new ArrayList<>();

    /**
     * Works together with the --cluster-window-size argument.
     */
    @Argument(fullName=CLUSTER_SIZE_LONG_NAME, shortName="cluster", doc="The number of SNPs which make up a cluster. Must be at least 2", optional=true)
    public Integer clusterSize = 3;

    /**
     * Works together with the --cluster-size argument.  To disable the clustered SNP filter, set this value to less than 1.
     */
    @Argument(fullName=CLUSTER_WINDOW_SIZE_LONG_NAME, shortName="window", doc="The window size (in bases) in which to evaluate clustered SNPs", optional=true)
    public Integer clusterWindow = 0;

    @Argument(fullName=MASK_EXTENSION_LONG_NAME, doc="How many bases beyond records from a provided 'mask' should variants be filtered", optional=true)
    public Integer maskExtension = 0;

    /**
     * When using the --mask argument, the mask-name will be annotated in the variant record.
     * Note that when using the --filter-not-in-mask argument to reverse the masking logic,
     * it is up to the user to adapt the name of the mask to make it clear that the reverse logic was used
     * (e.g. if masking against Hapmap, use --mask-name=hapmap for the normal masking and --mask-name=not_hapmap for the reverse masking).
     */
    @Argument(fullName=MASK_NAME_LONG_NAME, doc="The text to put in the FILTER field if a 'mask' is provided and overlaps with a variant call", optional=true)
    public String maskName = "Mask";

    /**
     * By default, if the --mask argument is used, any variant falling in a mask will be filtered.
     * If this argument is used, logic is reversed, and variants falling outside a given mask will be filtered.
     * Use case is, for example, if we have an interval list or BED file with "good" sites.
     * Note that it is up to the user to adapt the name of the mask to make it clear that the reverse logic was used
     * (e.g. if masking against Hapmap, use --mask-name=hapmap for the normal masking and --mask-name=not_hapmap for the reverse masking).
     */
    @Argument(fullName=FILTER_NOT_IN_MASK_LONG_NAME, doc="Filter records NOT in given input mask.", optional=true)
    public boolean filterRecordsNotInMask = false;

    /**
     * By default, if JEXL cannot evaluate your expression for a particular record because one of the annotations is not present, the whole expression evaluates as PASSing.
     * Use this argument to have it evaluate as failing filters instead for these cases.
     */
    @Argument(fullName=MISSING_VAL_LONG_NAME, doc="When evaluating the JEXL expressions, missing values should be considered failing the expression", optional=true)
    public Boolean failMissingValues = false;

    /**
     * Invalidate previous filters applied to the VariantContext, applying only the filters here
     */
    @Argument(fullName=StandardArgumentDefinitions.INVALIDATE_PREVIOUS_FILTERS_LONG_NAME, doc="Remove previous filters applied to the VCF", optional=true)
    boolean invalidatePreviousFilters = false;

    /**
     * Invert the selection criteria for --filter-expression
     */
    @Argument(fullName=INVERT_LONG_NAME, shortName="invfilter", doc="Invert the selection criteria for --filter-expression", optional=true)
    public boolean invertFilterExpression = false;

    /**
     * Invert the selection criteria for --genotype-filter-expression
     */
    @Argument(fullName=INVERT_GT_LONG_NAME, shortName="invG-filter", doc="Invert the selection criteria for --genotype-filter-expression", optional=true)
    public boolean invertGenotypeFilterExpression = false;

    /**
     * If this argument is provided, set filtered genotypes to no-call (./.).
     */
    @Argument(fullName=NO_CALL_GTS_LONG_NAME, optional=true, doc="Set filtered genotypes to no-call")
    public boolean setFilteredGenotypesToNocall = false;

    @Argument(fullName=ALLELE_SPECIFIC_LONG_NAME, optional=true, doc="Set mask at the allele level. This option is not compatible with clustering.")
    public boolean applyForAllele = false;

    /**
     * If a mask interval list is provided, then set the description of the filter in the VCF header to this String.
     * Note that if spaces are needed, then the entire description should be enclosed in quotes. Also note that if
     * --filter-not-in-mask is used, the description should be adapted to reflect the reverse logic.
     */
    @Argument(fullName=MASK_DESCRIPTION_LONG_NAME, optional=true, doc="Description to add to the FILTER field in VCF header for the mask filter.")
    public String maskDescription;

    // JEXL expressions for the filters
    private List<JexlVCMatchExp> filterExps;
    private List<JexlVCMatchExp> genotypeFilterExps;

    private JexlMissingValueTreatment howToTreatMissingValues;

    public static final String CLUSTERED_SNP_FILTER_NAME = "SnpCluster";

    private VariantContextWriter writer;

    private final List<Allele> diploidNoCallAlleles = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);

    /**
     * Invert logic if specified
     *
     * @param logic boolean logical operation value
     * @param invert whether to invert logic
     * @return invert logic if invert flag is true, otherwise leave the logic
     */
    private static boolean invertLogic(final boolean logic, final boolean invert){
        return invert ? !logic : logic;
    }

    /**
     * Prepend inverse phrase to description if --invert-filter-expression
     *
     * @param description the description
     * @return the description with inverse prepended if --invert_filter_expression
     */
    private String possiblyInvertFilterExpression( final String description ){
        return invertFilterExpression ? "Inverse of: " + description : description;
    }

    private void initializeVcfWriter() {
        writer = createVCFWriter(out);

        // setup the header fields
        final Set<VCFHeaderLine> hInfo = new LinkedHashSet<>();
        hInfo.addAll(getHeaderForVariants().getMetaDataInInputOrder());
        if (applyForAllele) {
            hInfo.add(new VCFInfoHeaderLine(AS_FILTER_STATUS_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.String, "Filter status for each allele, as assessed by ApplyVQSR. Note that the VCF filter field will reflect the most lenient/sensitive status across all alleles."));
        }

        // need AC, AN and AF since output if set filtered genotypes to no-call
        // If setting filtered genotypes to no-call, then allele counts (AC, AN and AF ) will be recomputed and these annotations
        // need to be included in the header
        if ( setFilteredGenotypesToNocall ) {
            GATKVariantContextUtils.addChromosomeCountsToHeader(hInfo);
        }

        if ( clusterWindow > 0 ) {
            hInfo.add(new VCFFilterHeaderLine(CLUSTERED_SNP_FILTER_NAME, "SNPs found in clusters"));
        }

        if ( !genotypeFilterExps.isEmpty() ) {
            hInfo.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_FILTER_KEY));
        }

        try {
            for ( final JexlVCMatchExp exp : filterExps ) {
                hInfo.add(new VCFFilterHeaderLine(exp.name, possiblyInvertFilterExpression(exp.exp.toString())));
            }
            for ( final JexlVCMatchExp exp : genotypeFilterExps ) {
                hInfo.add(new VCFFilterHeaderLine(exp.name, possiblyInvertFilterExpression(exp.exp.toString())));
            }

            if ( mask != null ) {
                if (maskDescription != null) {
                    hInfo.add(new VCFFilterHeaderLine(maskName, maskDescription));
                } else if (filterRecordsNotInMask) {
                    hInfo.add(new VCFFilterHeaderLine(maskName, "Doesn't overlap a user-input mask"));
                } else {
                    hInfo.add(new VCFFilterHeaderLine(maskName, "Overlaps a user-input mask"));
                }
            }
        } catch (final IllegalArgumentException e) {
            throw new UserException.BadInput(e.getMessage());
        }

        hInfo.addAll(getDefaultToolVCFHeaderLines());
        writer.writeHeader(new VCFHeader(hInfo, getHeaderForVariants().getGenotypeSamples()));
    }

    @Override
    public void onTraversalStart() {
        if (clusterSize <= 1){
            throw new CommandLineException.BadArgumentValue(CLUSTER_SIZE_LONG_NAME, "values lower than 2 are not allowed");
        }
        if ( maskExtension < 0 ) {
            throw new CommandLineException.BadArgumentValue(MASK_EXTENSION_LONG_NAME, "negative values are not allowed");
        }

        if (filterRecordsNotInMask && mask == null) {
            throw new CommandLineException.BadArgumentValue(FILTER_NOT_IN_MASK_LONG_NAME, "argument not allowed if mask argument is not provided");
        }
        if (maskDescription != null && mask == null) {
            throw new CommandLineException.BadArgumentValue(MASK_DESCRIPTION_LONG_NAME, "argument not allowed if mask argument is not provided");
        }
        filterExps = VariantContextUtils.initializeMatchExps(filterNames, filterExpressions);
        genotypeFilterExps = VariantContextUtils.initializeMatchExps(genotypeFilterNames, genotypeFilterExpressions);
        howToTreatMissingValues = failMissingValues ? JexlMissingValueTreatment.TREAT_AS_MATCH : JexlMissingValueTreatment.TREAT_AS_MISMATCH;

        VariantContextUtils.engine.get().setSilent(true);

        initializeVcfWriter();
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext ref, final FeatureContext featureContext) {
        if (applyForAllele) {
            final List<VariantContext> filtered = splitMultiAllelics(variant).stream().map(vc -> filter(vc, new FeatureContext(featureContext, new SimpleInterval(vc.getContig(), vc.getStart(), vc.getEnd())))).collect(Collectors.toList());
            // get filters for each allele
            final List<Set<String>> alleleFilters = filtered.stream().map(filteredvc -> filteredvc.getFilters()).collect(Collectors.toList());
            // add in the AS_FilterStatus and set the variant filters
            final VariantContext filteredVC = AlleleFilterUtils.addAlleleAndSiteFilters(variant, alleleFilters, invalidatePreviousFilters);
            writer.add(filteredVC);
        } else {
            writer.add(filter(variant, featureContext));
        }
    }

    private List<VariantContext> splitMultiAllelics(VariantContext vc) {
        final List<VariantContext> results = new ArrayList<>();
        final VariantContextBuilder vcb = new VariantContextBuilder("SimpleSplit", vc.getContig(), vc.getStart(), vc.getEnd(),
                    Arrays.asList(vc.getReference(), Allele.NO_CALL));
        vc.getAlternateAlleles().forEach(allele -> results.add(GATKVariantContextUtils.trimAlleles(
                    vcb.alleles(Arrays.asList(vc.getReference(), allele)).make(true), true, true)));
        return results;
    }

    /**
     * Add mask to variant context filters if it covers its location
     * @return VariantContext with the mask added if the VariantContext is within the extended mask area
     */
    private VariantContext addMaskIfCoversVariant(final VariantContext vc, final FeatureContext featureContext) {
        final List<Feature> maskVariants = featureContext.getValues(mask, maskExtension, maskExtension);

        final boolean variantsMasked = maskVariants.isEmpty() == filterRecordsNotInMask;
        if (variantsMasked) {
            final Set<String> oldFiltersPlusNewOne = Sets.union(vc.getFilters(), singleton(maskName));
            return new VariantContextBuilder(vc).filters(oldFiltersPlusNewOne).make();
        } else {
            return vc;
        }
    }

    private boolean isMaskFilterPresent(final VariantContext vc) {
        return vc.getFilters() != null && vc.getFilters().contains(maskName);
    }

    private VariantContext filter(final VariantContext variant, final FeatureContext featureContext) {
        final VariantContext vcModFilters = invalidatePreviousFilters ? (new VariantContextBuilder(variant)).unfiltered().make() : variant;
        final VariantContext vc = isMaskFilterPresent(vcModFilters) ? vcModFilters: addMaskIfCoversVariant(vcModFilters, featureContext);
        final VariantContextBuilder builder = new VariantContextBuilder(vc);

        // make new Genotypes based on filters
        if ( !genotypeFilterExps.isEmpty() || setFilteredGenotypesToNocall ) {
            GATKVariantContextUtils.setFilteredGenotypeToNocall(builder, vc, setFilteredGenotypesToNocall, this::getGenotypeFilters);
        }

        // make a new variant context based on filters
        //Note that making this empty set effectively converts the VC to PASS, whereas an unfiltered VC has null filters
        final Set<String> filters = new LinkedHashSet<>(vc.getFilters());

        // test for clustered SNPs if requested
        if (areClusteredSNPs(featureContext, vc)){
            filters.add(CLUSTERED_SNP_FILTER_NAME);
        }

        for ( final JexlVCMatchExp exp : filterExps ) {
            // Note that g is set to null since filterExps contains INFO-level filters
            if ( matchesFilter(vc, null, exp, invertFilterExpression) ) {
                filters.add(exp.name);
            }
        }

        //even if the original filters we null, we created a Set<String> filters above
        if ( filters.isEmpty() ) {
            if (!invalidatePreviousFilters) {
                builder.passFilters();
            } else {
                builder.unfiltered();
            }
        } else {
            builder.filters(filters);
        }

        return builder.make();
    }

    /**
     * Get the genotype filters
     *
     * @param vc the variant context
     * @param g the genotype
     * @return list of genotype filter names
     */
    private List<String> getGenotypeFilters(final VariantContext vc, final Genotype g) {
        final List<String> filters = new ArrayList<>();
        if (g.isFiltered()) {
            filters.addAll(Arrays.asList(g.getFilters().split(FILTER_DELIMITER)));
        }

        // Add if expression filters the variant context
        for (final JexlVCMatchExp exp : genotypeFilterExps) {
            // Give the genotype g to mathchesFilter enables JEXL matching on genotype fields
            if (matchesFilter(vc, g, exp, invertGenotypeFilterExpression)) {
                filters.add(exp.name);
            }
        }

        return filters;
    }

    /**
     * Return true if matches the filter expression
     */
    private boolean matchesFilter(final VariantContext vc, final Genotype g, final VariantContextUtils.JexlVCMatchExp exp, final boolean invertVCfilterExpression) {
        return invertLogic(VariantContextUtils.match(vc, g, exp, howToTreatMissingValues), invertVCfilterExpression);
    }

    /**
     * Return true if there is a window of size {@link #clusterWindow} that contains as least {@link #clusterSize} SNPs.
     */
    private boolean areClusteredSNPs(final FeatureContext featureContext, final VariantContext current) {
        if (clusterWindow < 1){ //as per argument doc, snpsInVicinity < 1 imply no clustering
            return false;
        }
        if (!current.isSNP()){
            return false;
        }
        //Need to fetch SNPs from left and right of the current SNPs.
        //Note: The "vicinity" here is the region around the current variant (up and down n bases, where n = clusterWindow).
        final List<VariantContext> snpsInVicinity = featureContext.getValues(getDrivingVariantsFeatureInput(), clusterWindow, clusterWindow)
                                                    .stream().filter(v -> v.isSNP()).collect(Collectors.toList());

        if (snpsInVicinity.size() < clusterSize){  //not enough variants - will never be a cluster no matter what.
            return false;
        }
        snpsInVicinity.sort(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR); //Note: by construction, they will be from the same contig.

        //Scan the list of variants and check if any stretch of clusterSize variants
        //that contains the 'current variant' has more variants than clusterWindow.
        final int currentStartPos = current.getStart();
        int firstIndex = 0;

        //start with the lowest pos
        int firstPos = snpsInVicinity.get(firstIndex).getStart();//there is at least 1 element so it's ok
        final int n = clusterSize - 1;
        while(firstPos <= currentStartPos && firstIndex + n < snpsInVicinity.size()) {
            final int firstPlusNPos = snpsInVicinity.get(firstIndex + n).getStart();

            if (firstPlusNPos - firstPos < clusterWindow){
                return true;
            }

            firstIndex++;
            firstPos = snpsInVicinity.get(firstIndex).getStart();//there is at least 1 element so it's ok
        }
        return false;
    }

    @Override
    public void closeTool() {
        if (writer != null){
            writer.close();
        }
    }

}
