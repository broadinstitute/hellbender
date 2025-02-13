package org.broadinstitute.hellbender.tools.walkers.longreads;

import htsjdk.samtools.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

/**
 * Recover uncorrected reads
 *
 * <h3>Input</h3>
 * <ul>
 *     <li>The uncorrected BAM file</li>
 *     <li>The CCS-corrected BAM file</li>
 * </ul>
 *
 * <h3>Outputs</h3>
 * <ul>
 *     <li>A BAM file of remaining subreads that were not corrected by the CCS step</li>
 * </ul>
 *
 * <h3>Usage Example</h3>
 * <h4>Quickly count errors</h4>
 * <pre>
 *   gatk RecoverUncorrectedReads \
 *     -I input.bam \
 *     -C ccs.bam \
 *     -O output.bam
 * </pre>
 */
@DocumentedFeature
@ExperimentalFeature
@CommandLineProgramProperties(
        summary = "Recover uncorrected reads",
        oneLineSummary = "Recover uncorrected reads",
        programGroup = ReadDataManipulationProgramGroup.class
)
public final class RecoverUncorrectedReads extends ReadWalker {
    @Argument(fullName = "corrected",
            shortName = "C",
            doc="CCS-corrected reads")
    public String corrected;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="Write output to this file")
    public String output;
    private SAMFileGATKReadWriter outputWriter;

    private final Set<String> correctedReadNames = new HashSet<>();
    private int uncorrectedReads = 0;

    @Override
    public void onTraversalStart() {
        final SamReaderFactory srf = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);

        try (final SamReader srs = srf.open(IOUtils.getPath(corrected))) {
            for ( final SAMRecord sr : srs ) {
                String[] name        = sr.getReadName().split("/");
                String   readZmwName = name[ 0 ] + "/" + name[ 1 ];
                correctedReadNames.add(readZmwName);
            }

        } catch (final IOException e) {
            throw new UserException("Unable to open path: " + corrected, e);
        }

        outputWriter = createSAMWriter(new GATKPath(IOUtils.getPath(output).toUri().toString()), true);
    }

    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        String[] name = read.getName().split("/");
        if (name.length >= 2) {
            String readZmwName = name[0] + "/" + name[1];
            if (!correctedReadNames.contains(readZmwName)) {
                outputWriter.addRead(read);
                uncorrectedReads++;

                logger.debug("uncorrected: {}", readZmwName);
            } else {
                logger.debug("corrected: {}", readZmwName);
            }
        } else {
            outputWriter.addRead(read);
            uncorrectedReads++;

            logger.debug("uncorrected: {}", read.getName());
        }
    }

    @Override
    public void closeTool() {
        logger.info("corrected: {} ; uncorrected: {}", correctedReadNames.size(), uncorrectedReads);

        if ( outputWriter != null ) {
            outputWriter.close();
        }
    }
}
