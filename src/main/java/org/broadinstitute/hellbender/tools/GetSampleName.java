package org.broadinstitute.hellbender.tools;

import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Emit a single sample name from the bam header into an output file. The sample name is that in the read group (RG) sample (SM) field
 *
 * <p>
 *     Note: If the bam has zero or more than one sample names in the header, this tool will error, by design.
 *     This tool has not been tested extensively.  Most options supported by the GATK are irrelevant for this tool.
 * </p>
 *
 * <h3>Input</h3>
 * <ul>
 *     <li>A BAM file with a single sample name in the header</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>A file with a single sample name in it</li>
 * </ul>
 *
 * <h3>Example Usage</h3>
 * <pre>
 *   gatk GetSampleName \
 *     -I input.bam \
 *     -O sample_name.txt
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Emit a single sample name from the bam header into an output file. " +
                "The sample name is that in the read group (RG) sample (SM) field",
        oneLineSummary = "Emit a single sample name",
        programGroup = DiagnosticsAndQCProgramGroup.class
)
final public class GetSampleName extends GATKTool {

    public static final String STANDARD_ENCODING = "UTF-8";

    @Argument(
            doc = "Output file with only the sample name in it.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    protected GATKPath outputSampleNameFile;

    public static final String URL_ENCODING_LONG_NAME = "use-url-encoding";
    public static final String URL_ENCODING_SHORT_NAME = "encode";

    @Argument(
            doc = "Apply URL encoding to convert spaces and other special characters in sample name.",
            fullName = URL_ENCODING_LONG_NAME,
            shortName = URL_ENCODING_SHORT_NAME
    )
    protected boolean urlEncode;

    @Override
    public void traverse() {
        // Do nothing!
    }

    @Override
    public boolean requiresReads() {return true;}

    @Override
    public void onTraversalStart() {

        // Grab the header info
        if ((getHeaderForReads() == null) || (getHeaderForReads().getReadGroups() == null)) {
            throw new UserException.BadInput("The given input bam has no header or no read groups.  Cannot determine a sample name.");
        }

        final List<String> sampleNames = getHeaderForReads().getReadGroups().stream().map(s -> s.getSample()).distinct().collect(Collectors.toList());

        if (sampleNames.size() == 0) {
            throw new UserException.BadInput("The given bam input has no sample names.");
        }

        try (final OutputStreamWriter fileWriter = new OutputStreamWriter(outputSampleNameFile.getOutputStream())) {
            final String outputSamplesOnSeparateLines = sampleNames.stream()
                    .map(rawSample -> urlEncode ? IOUtils.urlEncode(rawSample) : rawSample)
                    .collect(Collectors.joining("\n"));
            fileWriter.write(outputSamplesOnSeparateLines);
        } catch (final IOException ioe) {
            throw new UserException(String.format("Could not write to output file %s.", outputSampleNameFile), ioe);
        }
    }

}
