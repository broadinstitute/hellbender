package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

/**
 * Tools that are designed for use with long reads.
 * @author Jonn Smith
 */
public final class LongReadProgramGroup implements CommandLineProgramGroup {
    @Override
    public String getName() { return HelpConstants.DOC_CAT_LONGREAD; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_LONGREAD; }
}
