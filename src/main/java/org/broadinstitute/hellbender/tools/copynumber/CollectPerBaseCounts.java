package org.broadinstitute.hellbender.tools.copynumber;

import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.LocusWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.MappingQualityReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.datacollection.PerBaseCountCollector;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.PerBaseCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.Metadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.MetadataUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Collects per-base counts at specified sites where bases are one of those specified in {@link PerBaseCountCollector#BASES}.
 * If two reads share the same query name and, therefore, originate from the same template of DNA, this
 * only counts the higher base quality nucleotide of the two reads. If the base qualities of the two reads are equal,
 * the nucleotide from read one will be counted.
 * Only reads that pass the specified read filters, have mapping qualities greater than {@code minimum-mapping-quality},
 * and bases that exceed the specified {@code minimum-base-quality} will be counted.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         SAM format read data
 *     </li>
 *     <li>
 *         Reference FASTA file
 *     </li>
 *     <li>
 *         Sites at which per-base counts will be collected
 *     </li>
 * </ul>
 *
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         PerBase-counts file.
 *         This is a tab-separated values (TSV) file with a SAM-style header containing a read group sample name, a sequence dictionary,
 *         a row specifying the column headers contained in {@link PerBaseCountCollection.PerBaseCountTableColumn},
 *         and the corresponding entry rows.
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk CollectPerBaseCounts \
 *          -I sample.bam \
 *          -R reference.fa \
 *          -L sites.interval_list \
 *          -O sample.perBaseCounts.tsv
 * </pre>
 *
 * @author Lee Lichtenstein &lt;lichtens@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 * @author Robert Klein &lt;rklein@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Collects per-base counts at specified sites",
        oneLineSummary = "Collects per-base counts at specified sites",
        programGroup = CoverageAnalysisProgramGroup.class
)
@DocumentedFeature
public final class CollectPerBaseCounts extends LocusWalker {
    private static final int DEFAULT_MINIMUM_MAPPING_QUALITY = 30;
    static final int DEFAULT_MINIMUM_BASE_QUALITY = 20;

    static final List<ReadFilter> DEFAULT_ADDITIONAL_READ_FILTERS = ImmutableList.of(
            ReadFilterLibrary.MAPPED,
            ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT,
            ReadFilterLibrary.NOT_DUPLICATE,
            ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK,
            ReadFilterLibrary.NOT_SECONDARY_ALIGNMENT,
            ReadFilterLibrary.NOT_SUPPLEMENTARY_ALIGNMENT,
            new MappingQualityReadFilter(DEFAULT_MINIMUM_MAPPING_QUALITY));

    public static final String MINIMUM_BASE_QUALITY_LONG_NAME = "minimum-base-quality";

    @Argument(
            doc = "Output file for per-base counts.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputPerBaseCountsFile;

    @Argument(
            doc = "Minimum base quality.  Base calls with lower quality will be filtered out of pileups.",
            fullName = MINIMUM_BASE_QUALITY_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private int minimumBaseQuality = DEFAULT_MINIMUM_BASE_QUALITY;

    private PerBaseCountCollector perBaseCountCollector;

    @Override
    public boolean emitEmptyLoci() {
        return true;
    }

    @Override
    public boolean requiresReference() {
        return false;
    }

    @Override
    public boolean requiresIntervals() {
        return true;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> readFilters = new ArrayList<>(super.getDefaultReadFilters());
        readFilters.addAll(DEFAULT_ADDITIONAL_READ_FILTERS);
        return readFilters;
    }

    @Override
    public void onTraversalStart() {
        validateArguments();

        final SampleLocatableMetadata metadata = MetadataUtils.fromHeader(getHeaderForReads(), Metadata.Type.SAMPLE_LOCATABLE);
        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        //this check is currently redundant, since the master dictionary is taken from the reads;
        //however, if any other dictionary is added in the future, such a check should be performed
        if (!CopyNumberArgumentValidationUtils.isSameDictionary(metadata.getSequenceDictionary(), sequenceDictionary)) {
            logger.warn("Sequence dictionary in BAM does not match the master sequence dictionary.");
        }
        perBaseCountCollector = new PerBaseCountCollector(metadata);
        logger.info("Collecting per-base counts...");
    }

    private void validateArguments() {
        CopyNumberArgumentValidationUtils.validateOutputFiles(outputPerBaseCountsFile);
    }

    @Override
    public Object onTraversalSuccess() {
        logger.info(String.format("Writing allelic counts to %s...", outputPerBaseCountsFile.getAbsolutePath()));
        perBaseCountCollector.getPerBaseCounts().write(outputPerBaseCountsFile);

        logger.info(String.format("%s complete.", getClass().getSimpleName()));

        return null;
    }

    @Override
    public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        perBaseCountCollector.collectAtLocus(alignmentContext.getBasePileup(), alignmentContext.getLocation(), minimumBaseQuality);
    }
}

