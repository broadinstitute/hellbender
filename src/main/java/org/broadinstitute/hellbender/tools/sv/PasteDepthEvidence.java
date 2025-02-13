package org.broadinstitute.hellbender.tools.sv;

import htsjdk.tribble.Feature;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.codecs.*;

import java.util.ArrayList;
import java.util.List;

/**
 * <p>Merges depth evidence files into a single output file.</p>
 * <p>The locus-sorted input files must have the same number of rows and identical bins (i.e., intervals).</p>
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         One or more depth evidence files (extensions .rd.txt, .rd.txt.gz, or .rd.bci).
 *         These must be locus-sorted, and must all contain the same bins.
 *         <br />Or a .list file containing a list of such depth evidence files, one per line.
 *     </li>
 *     <li>
 *         A sequence dictionary (necessary only when all inputs are .rd.txt or .rd.txt.gz format).
 *     </li>
 *     <li>
 *         An optional bin-size.  Only intervals of this size will be output.
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         An output file containing merged depth evidence from the inputs.
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk PasteDepthEvidence \
 *       -F file1.rd.txt.gz [-F file2.rd.txt.gz ...] \
 *       --sequence-dictionary ref.dict \
 *       -O merged.rd.bci \
 *       --bin-size 100
 * </pre>
 *
 * @author Ted Sharpe &lt;tsharpe@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Merges multiple sources of depth evidence into a single output file.",
        oneLineSummary = "Merges DepthEvidence records.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@ExperimentalFeature
public class PasteDepthEvidence extends TabularMultiFeatureWalker<DepthEvidence> {
    public static final String EVIDENCE_FILE_LONG_ARGUMENT_NAME = "depth-evidence-file";
    public static final String COMPRESSION_LEVEL_ARGUMENT_NAME = "compression-level";
    public static final String BIN_SIZE_ARGUMENT_NAME = "bin-size";

    @Argument(
            doc = "Input DepthEvidence file URI(s) with extension .rd.txt, .rd.txt.gz, or rd.bci",
            fullName = EVIDENCE_FILE_LONG_ARGUMENT_NAME,
            shortName = StandardArgumentDefinitions.FEATURE_SHORT_NAME
    )
    private List<FeatureInput<DepthEvidence>> inputPaths;

    @Argument(
            doc = "Output file for features of a type matching the input. Will be indexed if it " +
                    "has a block-compressed extension (e.g. '.gz' or '.bci').",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private GATKPath outputFilePath;

    @Argument(
            doc = "Optional bin size.  If not specified, the \"correct\" bin size " +
            "will be taken from the first row.",
            fullName = BIN_SIZE_ARGUMENT_NAME,
            optional = true
    )
    private int binSize;

    @Argument(
            doc = "Output compression level",
            fullName = COMPRESSION_LEVEL_ARGUMENT_NAME,
            minValue = 0, maxValue = 9, optional = true
    )
    private int compressionLevel = 4;

    private FeatureSink<DepthEvidence> outputSink;
    private int[] aggregatedCounts;

    @Override
    @SuppressWarnings("unchecked")
    public void onTraversalStart() {
        final FeatureOutputCodec<? extends Feature, ? extends FeatureSink<? extends Feature>> codec =
                FeatureOutputCodecFinder.find(outputFilePath, DepthEvidence.class);
        final Class<? extends Feature> outputClass = codec.getFeatureType();
        if ( !DepthEvidence.class.isAssignableFrom(outputClass) ) {
            throw new UserException("Output file " + outputFilePath + " implies Feature subtype " +
                    outputClass.getSimpleName() + " but this tool writes DepthEvidence.");
        }
        outputSink = (FeatureSink<DepthEvidence>)codec.makeSortMerger(outputFilePath,
                            getDictionary(), new ArrayList<>(getSampleNames()), compressionLevel);
        aggregatedCounts = new int[getSampleNames().size()];
    }

    @Override
    public void apply( final List<DepthEvidence> features,
                       final List<Object> headers,
                       final ReadsContext readsContext,
                       final ReferenceContext referenceContext ) {
        final DepthEvidence firstFeature = features.get(0);
        if ( binSize == 0 ) {
            binSize = firstFeature.getLengthOnReference();
        }
        if ( firstFeature.getLengthOnReference() != binSize ) {
            return;
        }
        int idx = 0;
        for ( final DepthEvidence depthEvidence : features ) {
            if ( !firstFeature.getContig().equals(depthEvidence.getContig()) ||
                    firstFeature.getStart() != depthEvidence.getStart() ||
                    firstFeature.getEnd() != depthEvidence.getEnd() ) {
                throw new UserException("inputs do not have identical intervals at " + firstFeature);
            }
            final int[] counts = depthEvidence.getCounts();
            System.arraycopy(counts, 0, aggregatedCounts, idx, counts.length);
            idx += counts.length;
        }
        if ( idx != aggregatedCounts.length ) {
            throw new GATKException("didn't fill output counts array");
        }
        outputSink.write(new DepthEvidence(firstFeature, aggregatedCounts));
    }

    @Override
    public Object onTraversalSuccess() {
        outputSink.close();
        return null;
    }
}
