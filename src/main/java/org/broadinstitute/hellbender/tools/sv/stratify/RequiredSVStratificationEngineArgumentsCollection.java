package org.broadinstitute.hellbender.tools.sv.stratify;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

public class RequiredSVStratificationEngineArgumentsCollection extends SVStratificationEngineArgumentsCollection {
    private static final long serialVersionUID = 1L;

    /**
     * Expected format is tab-delimited and contains columns NAME, SVTYPE, MIN_SIZE, MAX_SIZE, track. First line must
     * be a header with column names. Comment lines starting with {@link TableUtils#COMMENT_PREFIX} are ignored.
     */
    @Argument(
            doc = "Stratification configuration file (.tsv)",
            fullName = STRATIFY_CONFIG_FILE_LONG_NAME
    )
    public GATKPath configFile;

}
