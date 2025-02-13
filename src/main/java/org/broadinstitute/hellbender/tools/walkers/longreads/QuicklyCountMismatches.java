package org.broadinstitute.hellbender.tools.walkers.longreads;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

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
 *   gatk QuicklyCountMismatches \
 *     -I input.bam \
 *     -O output.txt
 * </pre>
 */
@DocumentedFeature
@ExperimentalFeature
@CommandLineProgramProperties(
        summary = "Quickly count errors",
        oneLineSummary = "Quickly count errors",
        programGroup = ReadDataManipulationProgramGroup.class
)
public final class QuicklyCountMismatches extends LocusWalker {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="Path to which variants should be written")
    public String outPathName;

    @Argument(fullName = "insHist",
            shortName = "insHist",
            doc="Path to which variants should be written")
    public String insPathName;

    @Argument(fullName = "delHist",
            shortName = "delHist",
            doc="Path to which variants should be written")
    public String delPathName;

    private final Map<Integer, Long> numMismatchesMap = new TreeMap<>();
    private final Map<Integer, Long> numInsertionsMap = new TreeMap<>();
    private final Map<Integer, Long> numDeletionsMap = new TreeMap<>();
    private final Map<Integer, Long> numBasesMap   = new TreeMap<>();

    private final Map<Integer, Integer> insHist = new TreeMap<>();
    private final Map<Integer, Integer> delHist = new TreeMap<>();

    @Override
    public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        for (PileupElement pe : alignmentContext.getBasePileup()) {
            int npr = (pe.getRead().hasAttribute("np")) ? pe.getRead().getAttributeAsInteger("np") : 0;

            for (int np : Arrays.asList(-1, npr)) {
                if (!numMismatchesMap.containsKey(np)) { numMismatchesMap.put(np, 0L); }
                if (!numInsertionsMap.containsKey(np)) { numInsertionsMap.put(np, 0L); }
                if (!numDeletionsMap.containsKey(np)) { numDeletionsMap.put(np, 0L); }
                if (!numBasesMap.containsKey(np)) { numBasesMap.put(np, 0L); }

                if (pe.isDeletion() && pe.atStartOfCurrentCigar()) {
                    int len = pe.getCurrentCigarElement().getLength();
                    if (!delHist.containsKey(len)) { delHist.put(len, 0); }
                    delHist.put(len, delHist.get(len) + 1);

                    if (len > 1) {
                        numDeletionsMap.put(np, numDeletionsMap.get(np) + 1L);
                    }

                    //logger.info("{} {} {} {} {}", pe.getRead(), pe.getCurrentCigarElement(), pe.getCurrentCigarOffset(), pe.getOffsetInCurrentCigar(), pe);
                } else if (pe.isBeforeInsertion() && pe.getBasesOfImmediatelyFollowingInsertion() != null) {
                    int len = pe.getBasesOfImmediatelyFollowingInsertion().length();

                    if (!insHist.containsKey(len)) { insHist.put(len, 0); }
                    insHist.put(len, insHist.get(len) + 1);

                    if (len > 1) {
                        numInsertionsMap.put(np, numInsertionsMap.get(np) + 1L);
                    }
                } else {
                    if (pe.getBase() != referenceContext.getBase()) {
                        numMismatchesMap.put(np, numMismatchesMap.get(np) + 1L);
                    }
                }

                numBasesMap.put(np, numBasesMap.get(np) + 1L);
            }
        }
    }

    @Override
    public void closeTool() {
            try (final PrintStream out = new PrintStream(outPathName)) {

                for ( final int np : numMismatchesMap.keySet() ) {
                    out.println(String.join(" ", numMismatchesMap.get(np).toString(), numInsertionsMap.get(np).toString(), numDeletionsMap.get(np).toString(), numBasesMap.get(np).toString()));
                }
            } catch (final FileNotFoundException e) {
                throw new UserException("Cannot open given outfile: " + outPathName, e);
            }

            try (PrintStream outIns = new PrintStream(insPathName)) {
                for ( final int len : insHist.keySet() ) {
                    outIns.println(String.join(" ", Integer.valueOf(len).toString(), insHist.get(len).toString()));
                }
            } catch (final FileNotFoundException e) {
                throw new UserException("Cannot open given outfile: " + insPathName, e);
            }

            try (final PrintStream outDel = new PrintStream(delPathName)) {
                for ( final int len : delHist.keySet() ) {
                    outDel.println(String.join(" ", Integer.valueOf(len).toString(), delHist.get(len).toString()));
                }
            } catch (final FileNotFoundException e) {
                throw new UserException("Cannot open given outfile: " + delPathName, e);
            }
    }
}
