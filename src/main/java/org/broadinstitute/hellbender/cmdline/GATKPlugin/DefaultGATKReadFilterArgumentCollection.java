package org.broadinstitute.hellbender.cmdline.GATKPlugin;

import com.google.common.collect.Lists;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;

import java.util.*;

/**
 * Default {@link GATKReadFilterArgumentCollection} applied in GATK for optional read filters in the command line.
 * It contains arguments that allow the user to:
 *
 * - Provide a list of read filters to apply.
 * - Disable some and/or all read filters.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class DefaultGATKReadFilterArgumentCollection extends GATKReadFilterArgumentCollection {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME,
            shortName = ReadFilterArgumentDefinitions.READ_FILTER_SHORT_NAME,
            doc="Read filters to be applied before analysis", optional=true, common = true)
    public final List<String> userEnabledReadFilterNames = new ArrayList<>(); // preserve order

    @Argument(fullName = ReadFilterArgumentDefinitions.INVERTED_READ_FILTER_LONG_NAME,
            shortName = ReadFilterArgumentDefinitions.INVERTED_READ_FILTER_SHORT_NAME,
            doc="Inverted (with flipped acceptance/failure conditions) read filters applied before analysis (after regular read filters).", optional=true, common = true)
    public final List<String> userEnabledInvertedReadFilterNames = new ArrayList<>(); // preserve order

    @Argument(fullName = ReadFilterArgumentDefinitions.DISABLE_READ_FILTER_LONG_NAME,
            shortName = ReadFilterArgumentDefinitions.DISABLE_READ_FILTER_SHORT_NAME,
            doc="Read filters to be disabled before analysis", optional=true, common = true)
    public final List<String> userDisabledReadFilterNames = new ArrayList<>();

    @Advanced
    @Argument(fullName = ReadFilterArgumentDefinitions.DISABLE_TOOL_DEFAULT_READ_FILTERS,
            shortName = ReadFilterArgumentDefinitions.DISABLE_TOOL_DEFAULT_READ_FILTERS,
            doc = "Disable all tool default read filters (WARNING: many tools will not function correctly without their default read filters on)", common = true, optional = true)
    public boolean disableToolDefaultReadFilters = false;

    private List<String> fullUserEnabledReadFilterNames = new ArrayList<>();
    /** Returns the list with the union of inverted and non-inverted read filters provided by the user. Includes user default filters followed by the user inverted filters.  */
    public List<String> getAllUserEnabledReadFilterNames() {
        if (fullUserEnabledReadFilterNames.isEmpty()) {
            fullUserEnabledReadFilterNames.addAll(userEnabledReadFilterNames);
            fullUserEnabledReadFilterNames.addAll(userEnabledInvertedReadFilterNames);
        }
        return fullUserEnabledReadFilterNames;
    }

    /** Returns the list with the read filters provided by the user, preserving the order. */
    @Override
    public List<String> getUserEnabledReadFilterNames() {
        return userEnabledReadFilterNames;
    }

    /** Returns the list with the inverted read filters provided by the user, preserving the order. */
    @Override
    public List<String> getUserEnabledInvertedReadFilterNames() {
        return userEnabledInvertedReadFilterNames;
    }

    /** Returns the set of filters disabled by the user. */
    @Override
    public List<String> getUserDisabledReadFilterNames() {
        return userDisabledReadFilterNames;
    }

    /** {@inheritDoc}. */
    @Override
    public boolean getDisableToolDefaultReadFilters() {
        return disableToolDefaultReadFilters;
    }

}
