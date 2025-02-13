package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.Feature;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.codecs.FeatureOutputCodec;
import org.broadinstitute.hellbender.utils.codecs.FeatureOutputCodecFinder;
import org.broadinstitute.hellbender.utils.codecs.FeatureSink;

/**
 * <p>Combines adjacent intervals in DepthEvidence files.</p>
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         A locus-sorted DepthEvidence file (name ends with ".rd.txt", ".rd.txt.gz", or "rd.bci").
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         An output file containing merged evidence from the inputs.
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk CondenseDepthEvidence \
 *       -F input.rd.txt.gz \
 *       -O merged.rd.txt.gz
 * </pre>
 *
 * @author Ted Sharpe &lt;tsharpe@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Merges adjacent DepthEvidence records.",
        oneLineSummary = "Merges adjacent DepthEvidence records.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@ExperimentalFeature
public class CondenseDepthEvidence extends FeatureWalker<DepthEvidence> {
    public static final String INPUT_ARGNAME = "depth-evidence";
    public static final String COMPRESSION_LEVEL_ARGNAME = "compression-level";
    public static final String MAX_INTERVAL_SIZE_ARGNAME = "max-interval-size";
    public static final String MIN_INTERVAL_SIZE_ARGNAME = "min-interval-size";

    @Argument(
            doc = "DepthEvidence input file",
            fullName = INPUT_ARGNAME,
            shortName = StandardArgumentDefinitions.FEATURE_SHORT_NAME
    )
    private GATKPath inputPath;

    @Argument(
            doc = "Merged DepthEvidence output file",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private GATKPath outputPath;

    @Argument(
            doc = "Maximum interval size",
            fullName = MAX_INTERVAL_SIZE_ARGNAME,
            optional = true
    )
    private int maxIntervalLength = 1000;

    @Argument(
            doc = "Minimum interval size",
            fullName = MIN_INTERVAL_SIZE_ARGNAME,
            optional = true
    )
    private int minIntervalLength = 0;

    @Argument(
            doc = "Output compression level",
            fullName = COMPRESSION_LEVEL_ARGNAME,
            minValue = 0, maxValue = 9, optional = true
    )
    private int compressionLevel = 4;

    private FeatureSink<DepthEvidence> outputSink;
    private DepthEvidence accumulator;


    @Override
    public GATKPath getDrivingFeaturePath() {
        return inputPath;
    }

    @Override
    protected boolean isAcceptableFeatureType( final Class<? extends Feature> featureType ) {
        return featureType.equals(DepthEvidence.class);
    }

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();

        if ( minIntervalLength > maxIntervalLength ) {
            throw new UserException("Minimum interval length exceeds maximum interval length.");
        }

        final FeatureOutputCodec<DepthEvidence, ? extends FeatureSink<DepthEvidence>> codec =
                FeatureOutputCodecFinder.find(outputPath, DepthEvidence.class);
        final SVFeaturesHeader header = (SVFeaturesHeader)getDrivingFeaturesHeader();
        final SAMSequenceDictionary dict =
                header.getDictionary() != null ? header.getDictionary() : getBestAvailableSequenceDictionary();
        outputSink = codec.makeSink(outputPath, dict, header.getSampleNames(), compressionLevel);
    }

    @Override
    public void apply( final DepthEvidence feature,
                       final ReadsContext readsContext,
                       final ReferenceContext referenceContext,
                       final FeatureContext featureContext ) {
        if ( accumulator == null ) {
            accumulator = feature;
            return;
        }
        final int intervalLength = accumulator.getLengthOnReference();
        if ( !isAdjacent(accumulator, feature) || intervalLength >= maxIntervalLength ) {
            if ( intervalLength >= minIntervalLength ) {
                outputSink.write(accumulator);
            }
            accumulator = feature;
            return;
        }
        final int[] accumCounts = accumulator.getCounts();
        MathUtils.addToArrayInPlace(accumCounts, feature.getCounts());
        accumulator = new DepthEvidence(feature.getContig(), accumulator.getStart(),
                                            feature.getEnd(), accumCounts);
    }

    @Override
    public Object onTraversalSuccess() {
        super.onTraversalSuccess();
        if ( accumulator != null && accumulator.getLengthOnReference() >= minIntervalLength ) {
            outputSink.write(accumulator);
        }
        outputSink.close();
        return null;
    }

    private static boolean isAdjacent( final DepthEvidence ev1, final DepthEvidence ev2 ) {
        return ev1.getContig().equals(ev2.getContig()) && ev1.getEnd() + 1 == ev2.getStart();
    }
}
