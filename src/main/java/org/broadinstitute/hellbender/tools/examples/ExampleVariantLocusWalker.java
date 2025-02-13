package org.broadinstitute.hellbender.tools.examples;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.List;

/**
 * Example/toy program that shows how to implement the VariantLocusWalker interface. Prints locus-based coverage from
 * supplied variants, and reference bases if provided
 */
@CommandLineProgramProperties(
        summary = "Example tool that prints locus-based coverage from supplied variants to the specified output file (stdout if none provided), along with overlapping reference bases/features (if provided)",
        oneLineSummary = "Example tool that prints locus-based variants coverage with optional contextual data",
        programGroup = ExampleProgramGroup.class,
        omitFromCommandLine = true
)
public class ExampleVariantLocusWalker extends VariantLocusWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private File outputFile = null;

    private PrintStream outputStream = null;

    @Override
    public void onTraversalStart() {
        try {
            outputStream = outputFile != null ? new PrintStream(outputFile) : System.out;
        }
        catch ( final FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(outputFile, e);
        }
    }

    @Override
    public void apply(Locatable loc, List<VariantContext> variants, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        // print the locus
        outputStream.printf("Current locus %s:%d\n", loc.getContig(), loc.getStart());
        // print the reference context if available
        if ( referenceContext.hasBackingDataSource() ) {
            outputStream.println("\tReference base(s): " + new String(referenceContext.getBases()));
        }
        // print the variants
        for (final VariantContext variant : variants) {
            outputStream.printf("\tVariant at %s:%d-%d. Ref: %s Alt(s): %s\n",
                    variant.getContig(), variant.getStart(), variant.getEnd(), variant.getReference(), variant.getAlternateAlleles());
        }

        outputStream.println();
    }

    @Override
    public void closeTool() {
        if ( outputStream != null ) {
            outputStream.close();
        }
    }
}
