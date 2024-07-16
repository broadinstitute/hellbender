package org.broadinstitute.hellbender.tools;

import htsjdk.beta.io.bundle.*;
import htsjdk.beta.plugin.variants.VariantsBundle;
import htsjdk.samtools.util.FileExtensions;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.StandardOpenOption;
import java.util.*;

/**
 * Create a bundle (JSON) file for use with a GATK tool.
 *
 * other inputs are NEVER inferred, and must always be provided with a content type tag.
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Create a bundle (JSON) file for use with a GATK tool",
        oneLineSummary = "Create a bundle (JSON) file for use with a GATK tool",
        programGroup = OtherProgramGroup.class
)
public class CreateBundle extends CommandLineProgram {
    protected static final Logger logger = LogManager.getLogger(CreateBundle.class);

    public static final String SUPPRESS_INDEX_RESOLUTION_FULL_NAME = "suppress-index-resolution";
    public static final String OTHER_INPUT_FULL_NAME = "other-input";

    @Argument(fullName = StandardArgumentDefinitions.PRIMARY_INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.PRIMARY_INPUT_SHORT_NAME,
            doc="Path to the primary bundle input (content type will be inferred if no content type tag is specified)")
    GATKPath primaryInput;

    @Argument(fullName = StandardArgumentDefinitions.SECONDARY_INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.SECONDARY_INPUT_SHORT_NAME,
            doc = "Path to a secondary bundle input for" + StandardArgumentDefinitions.PRIMARY_INPUT_LONG_NAME
                    + "(usually an index). The type will be inferred if no content type tag is specified. If no "+
                    " secondary input is specified, an index for the primary bundle file will be automatically inferred unless " +
                    SUPPRESS_INDEX_RESOLUTION_FULL_NAME + " is specified.",
            optional = true)
    GATKPath secondaryInput;

    @Argument(fullName = OTHER_INPUT_FULL_NAME,
            shortName = OTHER_INPUT_FULL_NAME,
            doc = "Path to other bundle inputs for " + StandardArgumentDefinitions.PRIMARY_INPUT_LONG_NAME
                    + " A content type tag MUST be provided for each other input.",
            optional = true)
    List<GATKPath> otherInputs;

    @Argument(fullName = SUPPRESS_INDEX_RESOLUTION_FULL_NAME,
            doc ="Don't attempt to resolve the primary input's index file (defaults to false) if no secondary input is provided",
            optional = true)
    boolean suppressIndexResolution = false;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Path the output bundle file (must end with the suffix .json.")
    GATKPath outputBundlePath;

    private enum BundleType {
        VCF,
        REFERENCE,
        OTHER
    }
    private BundleType outputBundleType;

    @Override
    protected String[] customCommandLineValidation() {
        if (!outputBundlePath.toString().endsWith(".json")){
            return new String[]{"Output bundle path must end with the suffix .json"};
        }
        outputBundleType = determinePrimaryContentType();
        return super.customCommandLineValidation();
    }

    @Override
    protected Object doWork() {
        try (final BufferedWriter writer = Files.newBufferedWriter(outputBundlePath.toPath(), StandardOpenOption.CREATE)) {
            final Bundle bundle = switch (outputBundleType) {
                case VCF -> createVCFBundle();
                case REFERENCE ->
                    throw new IllegalArgumentException ("Reference bundles are not yet supported");
                case OTHER -> createOtherBundle();
            };
            writer.write(BundleJSON.toJSON(bundle));
        } catch (final IOException e) {
            throw new RuntimeException(String.format("Failed writing bundle to output %s", outputBundlePath), e);
        }
        return null;
    }

    private BundleType determinePrimaryContentType() {
        BundleType bundleType;

        // determine the type of bundle to create; consult the tag attributes if any, otherwise try to infer from the
        // primary input file extension
        final String primaryContentTag = primaryInput.getTag();
        if (primaryContentTag != null && !primaryContentTag.isEmpty()) {
            if (primaryContentTag.equals(BundleResourceType.CT_VARIANT_CONTEXTS)) {
                bundleType = BundleType.VCF;
            } else if (primaryContentTag.equals(BundleResourceType.CT_HAPLOID_REFERENCE)) {
                bundleType = BundleType.REFERENCE;
            } else {
                logger.info(String.format("Primary input content type %s for %s not recognized. A bundle will be created using content typse from the provided argument tags.",
                        primaryContentTag,
                        primaryInput));
                bundleType = BundleType.OTHER;
            }
        } else {
            logger.info(String.format("A content type for the primary input was not provided. Attempting to infer the content type from the %s extension.", primaryInput));
            bundleType = inferPrimaryContentType(primaryInput);
        }
        return bundleType;
    }

    private BundleType inferPrimaryContentType(final GATKPath primaryInput) {
        logger.info("Attempting to infer bundle content type from file extension.");
        if (FileExtensions.VCF_LIST.stream().anyMatch(ext -> primaryInput.hasExtension(ext))) {
            return BundleType.VCF;
        } else if (FileExtensions.FASTA.stream().anyMatch(ext -> primaryInput.hasExtension(ext))) {
            return BundleType.REFERENCE;
        } else {
            throw new IllegalArgumentException(String.format("Unable to infer bundle content type from file extension %s. A content type must be provided as part of the argument.", primaryInput));
        }
    }

    private Bundle createVCFBundle() {
        final Collection<BundleResource> bundleResources = new ArrayList<>();

        bundleResources.add(new IOPathResource(primaryInput, BundleResourceType.CT_VARIANT_CONTEXTS));
        if (secondaryInput != null) {
            final String secondaryContentType = secondaryInput.getTag();
            if (secondaryContentType == null) {
                logger.info(String.format("A content type for the secondary input was not provided. Assuming %s is an index.", secondaryInput));
                bundleResources.add(new IOPathResource(secondaryInput, BundleResourceType.CT_VARIANTS_INDEX));
            } else {
                bundleResources.add(new IOPathResource(secondaryInput, secondaryContentType));
            }
        } else if (!suppressIndexResolution) {
            // secondary input is null, and index resolution suppression is off
            final Optional<GATKPath> indexPath = VariantsBundle.resolveIndex(primaryInput, GATKPath::new);
            if (indexPath.isEmpty()) {
                throw new IllegalArgumentException(
                        String.format(
                                "Could not infer an index for %s, you must either specify the index path as a secondary input on the command line or specify the %s argument.",
                                primaryInput.getRawInputString(),
                                SUPPRESS_INDEX_RESOLUTION_FULL_NAME));
            }
            bundleResources.add(new IOPathResource(indexPath.get(), BundleResourceType.CT_VARIANTS_INDEX));
        }
        if (otherInputs != null) {
            for (final GATKPath otherInput : otherInputs) {
                final String otherContentType = otherInput.getTag();
                if (otherContentType == null) {
                    throw new IllegalArgumentException(
                            String.format(
                                    "A content must be provided for \"other\" input %s.",
                                    otherInput.getRawInputString()));
                } else {
                    bundleResources.add(new IOPathResource(otherInput, otherContentType));
                }
            }
        }
        return new VariantsBundle(bundleResources);
    }

    private Bundle createOtherBundle() {
        final Collection<BundleResource> bundleResources = new ArrayList<>();
        bundleResources.add(new IOPathResource(primaryInput, primaryInput.getTag()));
        if (secondaryInput != null) {
            final String secondaryContentType = secondaryInput.getTag();
            if (secondaryContentType == null) {
                throw new IllegalArgumentException(String.format("A content type for the secondary input must be provided."));
            } else {
                bundleResources.add(new IOPathResource(secondaryInput, secondaryContentType));
            }
        }
        if (otherInputs != null) {
            for (final GATKPath otherInput : otherInputs) {
                final String otherContentType = otherInput.getTag();
                if (otherContentType == null) {
                    throw new IllegalArgumentException(
                            String.format(
                                    "A content type must be provided for \"other\" input %s.",
                                    otherInput.getRawInputString()));
                } else {
                    bundleResources.add(new IOPathResource(otherInput, otherContentType));
                }
            }
        }
        return new Bundle(primaryInput.getTag(), bundleResources);
    }
}
