package org.broadinstitute.hellbender.tools.examples;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.spark.VariantLocusWalkerContext;
import org.broadinstitute.hellbender.engine.spark.VariantLocusWalkerSpark;

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
public final class ExampleVariantLocusWalkerSpark extends VariantLocusWalkerSpark {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private String outputFile = null;

    @Override
    protected void processVariants(JavaRDD<VariantLocusWalkerContext> rdd, JavaSparkContext ctx) {
        rdd.map(variantFunction()).saveAsTextFile(outputFile);
    }

    private static Function<VariantLocusWalkerContext, String> variantFunction() {
        return (Function<VariantLocusWalkerContext, String>) context -> {
            Locatable loc = context.getLocus();
            List<VariantContext> variants = context.getVariants();
            ReadsContext readsContext = context.getReadsContext();
            ReferenceContext referenceContext = context.getReferenceContext();
            FeatureContext featureContext = context.getFeatureContext();

            StringBuilder sb = new StringBuilder();
            sb.append(String.format("Current locus %s:%d\n", loc.getContig(), loc.getStart()));

            if ( referenceContext.hasBackingDataSource() ) {
                sb.append(String.format("\tReference base(s): %s\n", new String(referenceContext.getBases())));
            }

            for (final VariantContext variant : variants) {
                sb.append(String.format("\tVariant at %s:%d-%d. Ref: %s Alt(s): %s\n",
                        variant.getContig(), variant.getStart(), variant.getEnd(), variant.getReference(), variant.getAlternateAlleles()));
            }
            return sb.toString();
        };
    }
}
