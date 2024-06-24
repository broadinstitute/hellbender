package org.broadinstitute.hellbender.cmdline.GATKPlugin;

import java.io.Serializable;
import java.util.List;

/**
 * An abstract ArgumentCollection for defining the set of read filter descriptor plugin arguments that are exposed to the user on the command line.
 *
 * Subclasses should provide {@link org.broadinstitute.barclay.argparser.Argument} annotations for the arguments that should be exposed.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 * @see org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKReadFilterPluginDescriptor
 */
public abstract class GATKReadFilterArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    /**
     * Returns the enabled filter names. Order should be honored.
     */
    public abstract List<String> getUserEnabledReadFilterNames();

    /**
     * Returns the inverted enabled read filters. Applied after the regular read filters. Order should be honored.
     */
    public abstract List<String> getUserEnabledInvertedReadFilterNames();

    /**
     * Returns the union of all of the user enabled read filters and the user enabled inverted read filters.
     */
    public abstract List<String> getAllUserEnabledReadFilterNames();

    /**
     * Returns the disabled read filter names. Order should be honored.
     */
    public abstract List<String> getUserDisabledReadFilterNames();

    /**
     * Returns {@code true} if all default filters are disabled; {@code false} otherwise.
     */
    public abstract boolean getDisableToolDefaultReadFilters();

}
