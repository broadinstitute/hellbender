package org.broadinstitute.hellbender.tools.sv;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.PairWalker;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.GATKReadWriter;

import java.io.IOException;

@CommandLineProgramProperties(
        summary = "Shards BAM for expansion hunter.",
        oneLineSummary = "Shards BAM for expansion hunter.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@ExperimentalFeature
public class ShardBAMForExpansionHunter extends PairWalker {
    @Argument(
            doc = "Output path for relevant reads.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME )
    private GATKPath outputPath;

    private GATKReadWriter writer;

    @Override public void onTraversalStart() {
        super.onTraversalStart();
        writer = createSAMWriter(outputPath, false);
    }

    @Override public void apply( final GATKRead read, final GATKRead mate ) {
        writer.addRead(read);
        writer.addRead(mate);
    }

    @Override public void applyUnpaired( final GATKRead read ) {
        if ( !read.isSecondaryAlignment() && !read.isSupplementaryAlignment() ) {
            writer.addRead(read);
        }
    }

    @Override public void closeTool() {
        if ( writer != null ) {
            try {
                writer.close();
            } catch ( final IOException ioe ) {
                throw new UserException("failed to close output file", ioe);
            }
        }
        super.closeTool();
    }
}
