package org.broadinstitute.hellbender.tools.walkers.validation;

/**
 * Created by davidben on 3/2/17.
 */
public enum ConcordanceState {
    TRUE_POSITIVE("TP"), FALSE_POSITIVE("FP"), FALSE_NEGATIVE("FN"), FILTERED_TRUE_NEGATIVE("FTN"), FILTERED_FALSE_NEGATIVE("FFN"), TRUE_NEGATIVE("TN");

    private final String abbreviation;

    ConcordanceState(final String abbreviation) {
        this.abbreviation = abbreviation;
    }

    public String getAbbreviation() { return abbreviation; }
}
