package org.broadinstitute.hellbender.cmdline;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public final class ReadFilterArgumentDefinitions {
    private ReadFilterArgumentDefinitions(){}

    public static final String AMBIGUOUS_FILTER_FRACTION_NAME = "ambig-filter-frac";
    public static final String AMBIGUOUS_FILTER_BASES_NAME = "ambig-filter-bases";

    public static final String MAX_FRAGMENT_LENGTH_NAME = "max-fragment-length";

    public static final String LIBRARY_NAME = "library";

    public static final String MINIMUM_MAPPING_QUALITY_NAME = "minimum-mapping-quality";
    public static final String MAXIMUM_MAPPING_QUALITY_NAME = "maximum-mapping-quality";

    public static final String FILTER_TOO_SORT_NAME = "filter-too-short";
    public static final String DONT_REQUIRE_SOFT_CLIPS_BOTH_ENDS = "dont-require-soft-clips-both-ends";

    public static final String PL_FILTER_NAME_LONG_NAME = "platform-filter-name";

    public static final String BLACK_LISTED_LANES_LONG_NAME = "black-listed-lanes";

    public static final String BLACK_LIST_LONG_NAME = "read-group-black-list";

    public static final String KEEP_READ_GROUP_LONG_NAME = "keep-read-group";

    public static final String MAX_LENGTH_ARG_NAME = "max-read-length";
    public static final String MIN_LENGTH_ARG = "min-read-length";

    public static final String READ_NAME_LONG_NAME = "read-name";

    public static final String KEEP_REVERSE_STRAND_ONLY_NAME = "keep-reverse-strand-only";

    public static final String SAMPLE_NAME = "sample";
}
