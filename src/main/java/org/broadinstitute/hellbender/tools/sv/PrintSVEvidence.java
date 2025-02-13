package org.broadinstitute.hellbender.tools.sv;

import com.google.common.annotations.VisibleForTesting;
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

import java.util.*;
/**
 * <p>Merges locus-sorted files of evidence for structural variation into a single output file.</p>
 * <p>The tool can also subset the inputs to specified genomic intervals, or to a specified list of samples.</p>
 * <p>The evidence types and their files extensions are:</p>
 * <dl>
 *     <dt>BafEvidence</dt>
 *     <dd>The biallelic frequency of a SNP in some sample at some locus.
 *          File extensions are *.baf.txt, *.baf.txt.gz, or *.baf.bci.</dd>
 *     <dt>DepthEvidence</dt>
 *     <dd>The read counts of any number of samples on some interval.
 *          File extensions are *.rd.txt, *.rd.txt.gz, or *.rd.bci.</dd>
 *     <dt>DiscordantPairEvidence</dt>
 *     <dd>Evidence of a read pair that spans a genomic distance that's too large or too small.
 *          File extensions are *.pe.txt, *.pe.txt.gz, or *.pe.bci.</dd>
 *     <dt>SiteDepth</dt>
 *     <dd>The read counts of each base call for some sample at some locus.
 *          File extensions are *.sd.txt, *.sd.txt.gz, or *.sd.bci.</dd>
 *     <dt>SplitReadEvidence</dt>
 *     <dd>The number of chimeric reads in some sample at some locus.
 *          File extensions are *.sr.txt, *.sr.txt.gz, or *.sr.bci.</dd>
 * </dl>
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         One or more evidence files.
 *         These must be locus-sorted, and must all contain the same type of evidence.
 *         <br />Or a file containing a list of evidence files, one per line.
 *     </li>
 *     <li>Optional:  A list of sample names to extract.</li>
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
 *     gatk PrintSVEvidence \
 *       -F file1.baf.txt.gz [-F file2.baf.txt.gz ...] \
 *       -O merged.baf.bci \
 *       --sample-names sample1 [--sample-names sample2 ...]
 * </pre>
 *
 * @author Ted Sharpe &lt;tsharpe@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Merges multiple sources of SV evidence records of some particular feature type" +
        " into a single output file.  Inputs must be locus-sorted." +
        "  Can also subset by regions or samples.",
        oneLineSummary = "Merges SV evidence records.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@ExperimentalFeature
public class PrintSVEvidence extends MergingMultiFeatureWalker<SVFeature> {
    public static final String EVIDENCE_FILE_NAME = "evidence-file";
    public static final String SAMPLE_NAMES_NAME = "sample-names";
    public static final String COMPRESSION_LEVEL_NAME = "compression-level";

    @Argument(
            doc = "Input feature file URI(s) with extension '"
                    + SplitReadEvidenceCodec.FORMAT_SUFFIX + "', '"
                    + DiscordantPairEvidenceCodec.FORMAT_SUFFIX + "', '"
                    + SiteDepthCodec.FORMAT_SUFFIX + "', '"
                    + BafEvidenceCodec.FORMAT_SUFFIX + "', or '"
                    + DepthEvidenceCodec.FORMAT_SUFFIX + "' (may be gzipped). "
                    + "Can also handle bci rather than txt files.",
            fullName = EVIDENCE_FILE_NAME,
            shortName = StandardArgumentDefinitions.FEATURE_SHORT_NAME
    )
    private List<FeatureInput<SVFeature>> inputPaths;

    @Argument(doc = "List of sample names to extract from the sources (either as a .list file or " +
            " as repeated arguments).  If not specified, all samples will be merged.",
              fullName = SAMPLE_NAMES_NAME, optional = true)
    @VisibleForTesting
    Set<String> sampleNames = new LinkedHashSet<>();

    @Argument(
            doc = "Output file for features of a type matching the input. Will be indexed if it " +
                    "has a block-compressed extension (e.g. '.gz' or '.bci').",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private GATKPath outputFilePath;

    @Argument(
            doc = "Output compression level",
            fullName = COMPRESSION_LEVEL_NAME,
            minValue = 0, maxValue = 9, optional = true
    )
    private int compressionLevel = 4;

    private boolean noSampleFiltering = false;
    private FeatureSink<SVFeature> outputSink;

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();

        final FeatureOutputCodec<SVFeature,FeatureSink<SVFeature>> codec =
                FeatureOutputCodecFinder.find(outputFilePath, SVFeature.class);

        final Class<? extends SVFeature> outputClass = codec.getFeatureType();
        for ( final FeatureInput<SVFeature> input : inputPaths ) {
            try {
                final Class<? extends Feature> inputClass =
                        input.getFeatureCodecClass().getDeclaredConstructor().newInstance().getFeatureType();
                if ( !outputClass.isAssignableFrom(inputClass) ) {
                    throw new UserException("Incompatible feature input " + input.getFeaturePath() +
                            " produces features of type " + inputClass.getSimpleName() +
                            " rather than features of type " + outputClass.getSimpleName() +
                            " as dictated by the output path " + outputFilePath);
                }
            } catch ( final ReflectiveOperationException roe ) {
                throw new GATKException("Failed to instantiate codec " +
                                            input.getFeatureCodecClass().getSimpleName());
            }
        }
        if ( sampleNames.isEmpty() ) {
            // use the complete set of sample names we found in the headers of the feature files
            sampleNames.addAll(getSampleNames());
            if ( sampleNames.isEmpty() ) {
                noSampleFiltering = true;
            }
        }

        outputSink = codec.makeSortMerger(outputFilePath, getDictionary(),
                                            new ArrayList<>(sampleNames), compressionLevel);
    }

    @Override
    public void apply( final SVFeature featureArg,
                       final Object header,
                       final ReadsContext readsContext,
                       final ReferenceContext referenceContext ) {
        final SVFeature feature;
        if ( noSampleFiltering ) {
            feature = featureArg;
        } else {
            feature = featureArg.extractSamples(sampleNames, header);
            if ( feature == null ) {
                return;
            }
        }
        outputSink.write(feature);
    }

    @Override
    public Object onTraversalSuccess() {
        super.onTraversalSuccess();
        outputSink.close();
        return null;
    }
}
