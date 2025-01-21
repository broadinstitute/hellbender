package org.broadinstitute.hellbender.tools.sv.stratify;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.GATKPath;

import java.io.Serializable;
import java.util.List;

/**
 * Arguments for use with {@link SVStratificationEngine}.
 */
public class SVStratificationEngineArgumentsCollection implements Serializable {
    // Command-line arguments
    public static final String STRATIFY_CONFIG_FILE_LONG_NAME = "stratify-config";
    public static final String TRACK_NAME_FILE_LONG_NAME = "track-name";
    public static final String TRACK_INTERVAL_FILE_LONG_NAME = "track-intervals";
    public static final String OVERLAP_FRACTION_LONG_NAME = "stratify-overlap-fraction";
    public static final String NUM_BREAKPOINT_OVERLAPS_LONG_NAME = "stratify-num-breakpoint-overlaps";
    public static final String NUM_BREAKPOINT_INTERCHROM_OVERLAPS_LONG_NAME = "stratify-num-breakpoint-overlaps-interchromosomal";
    private static final long serialVersionUID = 1L;

    @Argument(
            doc = "Track intervals file. Can be specified multiple times.",
            fullName = TRACK_INTERVAL_FILE_LONG_NAME,
            optional = true
    )
    public List<GATKPath> trackFileList;

    @Argument(
            doc = "Track names. Must be once for each --" + TRACK_INTERVAL_FILE_LONG_NAME,
            fullName = TRACK_NAME_FILE_LONG_NAME,
            optional = true
    )
    public List<String> trackNameList;

    @Argument(
            doc = "Minimum overlap fraction for tracks",
            minValue = 0,
            maxValue = 1,
            fullName = OVERLAP_FRACTION_LONG_NAME
    )
    public double overlapFraction = 0;

    @Argument(
            doc = "Minimum number of variant endpoint overlaps for tracks",
            minValue = 0,
            maxValue = 2,
            fullName = NUM_BREAKPOINT_OVERLAPS_LONG_NAME
    )
    public int numBreakpointOverlaps = 1;

    @Argument(
            doc = "Minimum number of breakpoint overlaps for tracks for interchromosomal variants (e.g. BNDs)",
            minValue = 1,
            maxValue = 2,
            fullName = NUM_BREAKPOINT_INTERCHROM_OVERLAPS_LONG_NAME
    )
    public int numBreakpointOverlapsInterchrom = 1;
}
