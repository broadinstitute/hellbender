package org.broadinstitute.hellbender.tools.walkers.genotyper;

/**
 * Created by davidben on 6/10/16.
 */
public enum GenotypeAssignmentMethod {

    /**
     * set all of the genotype GT values to NO_CALL
     */
    SET_TO_NO_CALL,

    /**
     * Use the subsetted PLs to greedily assign genotypes
     */
    USE_PLS_TO_ASSIGN,

    /**
     * Use the existing, subsetted posteriors array to assign genotypes
     */
    USE_POSTERIORS_ANNOTATION,

    /**
     * set all of the genotype GT values to NO_CALL and remove annotations
     */
    SET_TO_NO_CALL_NO_ANNOTATIONS,

    /**
     * Try to match the original GT calls, if at all possible
     *
     * Suppose I have 3 alleles: A/B/C and the following samples:
     *
     *       original_GT best_match to A/B best_match to A/C
     * S1 => A/A A/A A/A
     * S2 => A/B A/B A/A
     * S3 => B/B B/B A/A
     * S4 => B/C A/B A/C
     * S5 => C/C A/A C/C
     *
     * Basically, all alleles not in the subset map to ref.  It means that het-alt genotypes
     * when split into 2 bi-allelic variants will be het in each, which is good in some cases,
     * rather than the undetermined behavior when using the PLs to assign, which could result
     * in hom-var or hom-ref for each, depending on the exact PL values.
     */
    BEST_MATCH_TO_ORIGINAL,

    /**
     * do not even bother changing the GTs
     */
    DO_NOT_ASSIGN_GENOTYPES,

    /**
     * Calculate posterior probabilities and use those to assign genotypes
     */
    USE_POSTERIOR_PROBABILITIES,

    /**
     * Use PLs unless they are unavailable, in which case use best match to original
     * GQ0 hom-refs will be converted to no-calls
     */
    PREFER_PLS
}
