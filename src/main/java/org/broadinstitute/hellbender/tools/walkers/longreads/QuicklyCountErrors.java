package org.broadinstitute.hellbender.tools.walkers.longreads;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.LongReadProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.FileNotFoundException;
import java.io.PrintStream;

/**
 * Quickly count errors
 *
 * <h3>Input</h3>
 * <ul>
 *     <li>A single BAM file</li>
 * </ul>
 *
 * <h3>Outputs</h3>
 * <ul>
 *     <li>A collection of BAM files where care has been keep reads from the same ZMW in the same file</li>
 * </ul>
 *
 * <h3>Usage Example</h3>
 * <h4>Quickly count errors</h4>
 * <pre>
 *   gatk QuicklyCountErrors \
 *     -I input.bam \
 *     -O output.txt
 * </pre>
 */
@DocumentedFeature
@ExperimentalFeature
@CommandLineProgramProperties(
        summary = "Quickly count errors",
        oneLineSummary = "Quickly count errors",
        programGroup = LongReadProgramGroup.class
)
public final class QuicklyCountErrors extends ReadWalker {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="Path to which variants should be written")
    public String outPathName;
    //final Path outPath = IOUtils.getPath(outPathName);

    private long numMismatches = 0;
    private long numInsertions = 0;
    private long numDeletions = 0;
    private long numBases = 0;

    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        boolean xseen = false;
        for (CigarElement ce : read.getCigarElements()) {
            if (ce.getOperator().equals(CigarOperator.X)) {
                xseen = true;
                numMismatches++;
            } else if (ce.getOperator().equals(CigarOperator.D)) {
                numDeletions += ce.getLength();
            } else if (ce.getOperator().equals(CigarOperator.I)) {
                numInsertions += ce.getLength();
            }
        }

        if (!xseen && read.hasAttribute("NM")) {
            numMismatches += read.getAttributeAsInteger("NM");
        }

        numBases += read.getLength();
    }

    @Override
    public void closeTool() {
        try ( final PrintStream out = new PrintStream(outPathName) ) {
            out.println(String.join(" ", "numMismatches", "numDeletions", "numInsertions", "numBases"));
            out.println(String.join(" ", Long.valueOf(numMismatches).toString(), Long.valueOf(numDeletions).toString(), Long.valueOf(numInsertions).toString(), Long.valueOf(numBases).toString()));
        } catch (final FileNotFoundException e) {
            throw new UserException("Cannot open given outfile: " + outPathName, e);
        }

        logger.info("numMismatches: " + numMismatches + " numDeletions: " + numDeletions + " numInsertions: " + numInsertions + " numBases:" + numBases);
    }
}
