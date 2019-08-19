package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.MultiplePassVariantWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.contamination.CalculateContamination;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;
import org.broadinstitute.hellbender.tools.walkers.readorientation.LearnReadOrientationModel;

import java.io.File;
import java.nio.file.Path;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * <p>Filter variants in a Mutect2 VCF callset.</p>
 *
 * <p>
 *     FilterMutectCalls applies filters to the raw output of {@link Mutect2}.
 *     Parameters are contained in {@link M2FiltersArgumentCollection} and described in
 *     <a href='https://github.com/broadinstitute/gatk/tree/master/docs/mutect/mutect.pdf' target='_blank'>https://github.com/broadinstitute/gatk/tree/master/docs/mutect/mutect.pdf</a>.
 *     To filter based on sequence context artifacts, specify the --orientation-bias-artifact-priors [artifact priors tar.gz file] argument
 *     one or more times.  This input is generated by {@link LearnReadOrientationModel}.
 * </p>
 * <p>
 *     If given a --contamination-table file, e.g. results from {@link CalculateContamination}, the tool will additionally
 *     filter variants due to contamination. This argument may be specified with a table for one or more tumor samples.
 *     Alternatively, provide an estimate of the contamination with the --contamination argument.
 *
 *     FilterMutectCalls can also be given one or more --tumor-segmentation files, which are also output by {@link CalculateContamination}.
 * </p>
 * <p>
 *     This tool is featured in the Somatic Short Mutation calling Best Practice Workflow.
 *     See <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11136">Tutorial#11136</a> for a
 *     step-by-step description of the workflow and <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11127">Article#11127</a>
 *     for an overview of what traditional somatic calling entails. For the latest pipeline scripts, see the
 *     <a href="https://github.com/broadinstitute/gatk/tree/master/scripts/mutect2_wdl">Mutect2 WDL scripts directory</a>.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk FilterMutectCalls \
 *   -R reference.fasta \
 *   -V somatic.vcf.gz \
 *   --contamination-table contamination.table \
 *   --tumor-segmentation segments.tsv \
 *   -O filtered.vcf.gz
 * </pre>
 *
 *
 * When running on unfiltered output of Mutect2 in --mitochondria mode, setting the advanced option --autosomal-coverage
 * argument (default 0) activates a recommended filter against likely erroneously mapped  <a href="https://en.wikipedia.org/wiki/NUMT">NuMTs (nuclear mitochondrial DNA segments)</a>.
 * For the value, provide the median coverage expected in autosomal regions with coverage.
 *
 */
@CommandLineProgramProperties(
        summary = "Filter somatic SNVs and indels called by Mutect2",
        oneLineSummary = "Filter somatic SNVs and indels called by Mutect2",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
public final class FilterMutectCalls extends MultiplePassVariantWalker {

    public static final String FILTERING_STATS_LONG_NAME = "filtering-stats";

    public static final String FILTERING_STATUS_VCF_KEY = "filtering_status";

    public static final String FILTERING_STATS_EXTENSION = ".filteringStats.tsv";

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName =StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output filtered VCF file", optional=false)
    private final String outputVcf = null;

    @Argument(fullName = Mutect2.MUTECT_STATS_SHORT_NAME, doc="The Mutect stats file output by Mutect2", optional=true)
    private final String statsTable = null;

    @Argument(fullName = FILTERING_STATS_LONG_NAME, doc="The output filtering stats file", optional=true)
    private final String filteringStatsOutput = null;

    @ArgumentCollection
    protected M2FiltersArgumentCollection MTFAC = new M2FiltersArgumentCollection();

    private VariantContextWriter vcfWriter;

    private Mutect2FilteringEngine filteringEngine;

    private static final int NUMBER_OF_LEARNING_PASSES = 2;

    @Override
    protected int numberOfPasses() { return NUMBER_OF_LEARNING_PASSES + 2; }    // {@coode NUMBER_OF_LEARNING_PASSES} passes for learning, one for the threshold, and one for calling

    @Override
    public boolean requiresReference() { return true;}

    @Override
    public void onTraversalStart() {
        final VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> headerLines = inputHeader.getMetaDataInSortedOrder().stream()
                .filter(line -> !line.getKey().equals(FILTERING_STATUS_VCF_KEY)) //remove header line from Mutect2 stating that calls are unfiltered.
                .collect(Collectors.toSet());
        headerLines.add(new VCFHeaderLine(FILTERING_STATUS_VCF_KEY, "These calls have been filtered by " + FilterMutectCalls.class.getSimpleName() + " to label false positives with a list of failed filters and true positives with PASS."));

        GATKVCFConstants.MUTECT_FILTER_NAMES.stream().map(GATKVCFHeaderLines::getFilterLine).forEach(headerLines::add);

        headerLines.addAll(getDefaultToolVCFHeaderLines());

        final VCFHeader vcfHeader = new VCFHeader(headerLines, inputHeader.getGenotypeSamples());
        vcfWriter = createVCFWriter(new File(outputVcf));
        vcfWriter.writeHeader(vcfHeader);


        final File mutect2StatsTable = new File(statsTable == null ? drivingVariantFile + Mutect2.DEFAULT_STATS_EXTENSION : statsTable);
        filteringEngine = new Mutect2FilteringEngine(MTFAC, vcfHeader, mutect2StatsTable);
        if (!mutect2StatsTable.exists()) {
            throw new UserException.CouldNotReadInputFile("Mutect stats table " + mutect2StatsTable + " not found.  When Mutect2 outputs a file calls.vcf it also creates" +
                    " a calls.vcf" + Mutect2.DEFAULT_STATS_EXTENSION + " file.  Perhaps this file was not moved along with the vcf, or perhaps it was not delocalized from a" +
                    " virtual machine while running in the cloud." );
        }
    }

    @Override
    protected void nthPassApply(final VariantContext variant,
                                final ReadsContext readsContext,
                                final ReferenceContext referenceContext,
                                final FeatureContext featureContext,
                                final int n) {
        ParamUtils.isPositiveOrZero(n, "Passes must start at the 0th pass.");
        if (n <= NUMBER_OF_LEARNING_PASSES) {
            filteringEngine.accumulateData(variant, referenceContext);
        } else if (n == NUMBER_OF_LEARNING_PASSES + 1) {
            vcfWriter.add(filteringEngine.applyFiltersAndAccumulateOutputStats(variant, referenceContext));
        } else {
            throw new GATKException.ShouldNeverReachHereException("This walker should never reach (zero-indexed) pass " + n);
        }
    }

    @Override
    protected void afterNthPass(final int n) {
        if (n < NUMBER_OF_LEARNING_PASSES) {
            filteringEngine.learnParameters();
        } else if (n == NUMBER_OF_LEARNING_PASSES) {
            // it's important for filter parameters to stay the same and only learn the threshold in the final pass so that the
            // final threshold used corresponds exactly to the filters
            filteringEngine.learnThreshold();
        }else if (n == NUMBER_OF_LEARNING_PASSES + 1) {
            final Path filteringStats = IOUtils.getPath(
                filteringStatsOutput != null ? filteringStatsOutput
                    : outputVcf + FILTERING_STATS_EXTENSION);
            filteringEngine.writeFilteringStats(filteringStats);
        } else {
            throw new GATKException.ShouldNeverReachHereException("This walker should never reach (zero-indexed) pass " + n);
        }
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }

}
