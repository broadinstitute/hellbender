package org.broadinstitute.hellbender.utils.help;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.HashMap;
import java.util.Map;

public final class HelpConstants {

    private HelpConstants() {};

    public static final String LEGACY_GATK_FORUM_URL = "http://gatkforums.broadinstitute.org/";
    public static final String GATK_FORUM_URL = "https://gatk.broadinstitute.org/hc/en-us/community/topics";
    public static final String GATK_MAIN_SITE = "https://software.broadinstitute.org/gatk/";
    public static final String GATK_DOCS_SITE = "https://gatk.broadinstitute.org/hc/";

    //TODO anything using this is a dead link, we should move any remaining references to this into some technical docs for the future.
    public static String forumPost(String post) {
        return LEGACY_GATK_FORUM_URL + post;
    }

    public static String forumArticle(String article) {
        return GATK_DOCS_SITE + "articles/" + article;
    }

    /**
     * Definition of the group names / descriptions for documentation/help purposes.
     *
     */

    // Start GATK Program groups
    public final static String DOC_CAT_ANNOTATORS = "Variant Annotations";
    public final static String DOC_CAT_ANNOTATORS_SUMMARY = "Available to HaplotypeCaller, Mutect2, VariantAnnotator and GenotypeGVCFs. " +
            "See https://software.broadinstitute.org/gatk/documentation/article?id=10836";

    public final static String DOC_CAT_FLOW_ANNOTATORS = "Flow Annotations";
    public final static String DOC_CAT_FLOW_ANNOTATORS_SUMMARY = "Flow annotataions available to HaplotypeCaller, Mutect2, VariantAnnotator and GenotypeGVCFs";


    public final static String DOC_CAT_COVERAGE_ANALYSIS = "Coverage Analysis";
    public final static String DOC_CAT_COVERAGE_ANALYSIS_SUMMARY = "Tools that count coverage, e.g. depth per allele";

    public final static String DOC_CAT_CNV = "Copy Number Variant Discovery";
    public final static String DOC_CAT_CNV_SUMMARY = "Tools that analyze read coverage to detect copy number variants.";

    public final static String DOC_CAT_LONGREAD = "Long Read Tools";
    public final static String DOC_CAT_LONGREAD_SUMMARY = "Tools that are designed for use with long reads.";


    public final static String DOC_CAT_EXAMPLE = "Example Tools";
    public final static String DOC_CAT_EXAMPLE_SUMMARY = "Example tools that show developers how to implement new tools";

    public final static String DOC_CAT_METAGENOMICS = "Metagenomics";
    public final static String DOC_CAT_METAGENOMICS_SUMMARY = "Tools that perform metagenomic analysis, e.g. microbial community composition and pathogen detection";

    public final static String DOC_CAT_READFILTERS = "Read Filters";
    public final static String DOC_CAT_READFILTERS_SUMMARY = "Applied by engine to select reads for analysis";

    public final static String DOC_CAT_SHORT_VARIANT_DISCOVERY = "Short Variant Discovery";
    public final static String DOC_CAT_SHORT_VARIANT_DISCOVERY_SUMMARY = "Tools that perform variant calling and genotyping for short variants (SNPs, SNVs and Indels)";

    public final static String DOC_CAT_SV_DISCOVERY = "Structural Variant Discovery";
    public final static String DOC_CAT_SV_DISCOVERY_SUMMARY = "Tools that detect structural variants";

    public final static String DOC_CAT_TEST = "Test Tools";
    public final static String DOC_CAT_TEST_SUMMARY = "Tools for internal test purposes";

    public final static String DOC_CAT_RNA = "RNA-Specific Tools";
    public final static String DOC_CAT_RNA_SUMMARY = "Tools intended to be used for processing RNA data.";

    public final static String DOC_CAT_METHYLATION_DISCOVERY = "Methylation-Specific Tools";
    public final static String DOC_CAT_METHYLATION_DISCOVERY_SUMMARY = "Tools that perform methylation calling, processing bisulfite sequenced, methylation-aware aligned BAM";

    public final static String DOC_CAT_SHORT_FLOW_BASED = "Flow Based Tools";
    public final static String DOC_CAT_SHORT_FLOW_BASED_SUMMARY = "Tools designed specifically to operate on flow based data";

    // End GATK Program groups

    /**
     * List of "supercategory" values used for doc purposes. Every doc group name can/should be put into
     * one of the following super-categories.
     */
    private final static String DOC_SUPERCAT_TOOLS = "tools";
    private final static String DOC_SUPERCAT_UTILITIES = "utilities";
    private final static String DOC_SUPERCAT_EXCLUDE = "exclude";

    private static Map<String, String> groupToSuperCategory;

    private static Map<String, String> getSuperCategoryMap() {
        if (groupToSuperCategory == null) {

            // do this only on demand since we only need it during docgen
            groupToSuperCategory = new HashMap<>();

            // supercat Tools

            // start Picard program groups
            groupToSuperCategory.put(new picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup().getName(), DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(new picard.cmdline.programgroups.IntervalsManipulationProgramGroup().getName(), DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(new picard.cmdline.programgroups.OtherProgramGroup().getName(), DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(new picard.cmdline.programgroups.ReadDataManipulationProgramGroup().getName(), DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(new picard.cmdline.programgroups.ReferenceProgramGroup().getName(), DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(new picard.cmdline.programgroups.VariantEvaluationProgramGroup().getName(), DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(new picard.cmdline.programgroups.VariantFilteringProgramGroup().getName(), DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(new picard.cmdline.programgroups.VariantManipulationProgramGroup().getName(), DOC_SUPERCAT_TOOLS);
            // end Picard program groups

            groupToSuperCategory.put(DOC_CAT_CNV, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_COVERAGE_ANALYSIS, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_METAGENOMICS, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_SHORT_VARIANT_DISCOVERY, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_SV_DISCOVERY, DOC_SUPERCAT_TOOLS);

            // supercat Utilities
            groupToSuperCategory.put(DOC_CAT_READFILTERS, DOC_SUPERCAT_UTILITIES);
            groupToSuperCategory.put(DOC_CAT_ANNOTATORS, DOC_SUPERCAT_UTILITIES);
            groupToSuperCategory.put(DOC_CAT_FLOW_ANNOTATORS, DOC_SUPERCAT_UTILITIES);

            // supercat Exclude
            groupToSuperCategory.put(DOC_CAT_EXAMPLE, DOC_SUPERCAT_EXCLUDE);
            groupToSuperCategory.put(DOC_CAT_TEST, DOC_SUPERCAT_EXCLUDE);
        }
        return groupToSuperCategory;
    }

    /**
     * Given a group name, return a supercategory string for use by the online doc system to determine which
     * supercateogry the group is in. The strings returned by this method should match those used in the
     * corresponding help template.
     *
     * @param groupName
     * @return supercategory string corresponding to {@code groupName} for use in for use in determining
     * which supercateogry the group is in for online doc purposes.
     */
    public static String getSuperCategoryProperty(final String groupName) {
        Utils.nonNull(groupName);
        return getSuperCategoryMap().getOrDefault(groupName, "other");
    }

}
