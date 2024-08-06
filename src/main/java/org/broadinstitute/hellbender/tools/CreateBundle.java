package org.broadinstitute.hellbender.tools;

import htsjdk.beta.io.bundle.*;
import htsjdk.beta.plugin.registry.HaploidReferenceResolver;
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
 * <p>
 * Create a single bundle (JSON) file for one or more related inputs, for use as input to a GATK tool.
 * <h2>Bundle Structure</h2>
 * A bundle is a JSON file that contains references to one or more related resources (i.e., a VCF file
 * and it's associated index file, or a .fasta reference file and it's associated index and dictionary files).
 * Bundle files can be supplied as inputs for many GATK tools. The primary advantage to using a bundle input
 * over an individual file input is that a bundle allows the individual resources to be located in different
 * directories (ie., when providing a VCF input, it's corresponding index file generally need to be a sibling
 * file in the same directory as the VCF, whereas using a bundle file, because you reference each resource
 * explicitly, the index can be located in a different directory from the VCF).
 * <p>
 * Each resource in a bundle has an associated content type, which is a string that identifies the type of data
 * in that resource. One resource in the bundle is always designated as the "primary" resource, which determines
 * the type of the bundle.
 *
 * <p>An example bundle JSON file is show below. The bundle has 3 resources, with content types
 * "HAPLOID_REFERENCE", "REFERENCE_DICTIONARY", and "REFERENCE_INDEX" respectively, with the "HAPLOID_REFERENCE"
 * resource designated as the primary resource. This bundle is a reference bundle, and can be used as an
 * input where ever a reference argument is required.
 * <pre>
 *  {
 *    "schema": {
 *     "schemaVersion": "0.1.0",
 *     "schemaName": "htsbundle"
 *    },
 *    "HAPLOID_REFERENCE": {"path": "file:///projects/print_reads.fasta"},
 *    "REFERENCE_DICTIONARY": {"path": "file:///projects/print_reads.dict"},
 *    "REFERENCE_INDEX": {"path": "file:///projects/print_reads.fasta.fai"},
 *    "primary": "HAPLOID_REFERENCE"
 *  }
 * </pre>
 * <h2>Using CreateBundle</h2>
 * <p>
 * CreateBundle requires at least one primary resource input. One or more optional secondary resource inputs
 * may also be supplied.
 * <p>
 * The simplest way to use CreateBundle is to specify only the primary resource. In this case, with no secondary
 * resources are explicitly provided, CreateBundle will attempt to locate and infer "standard" secondary resources
 * (see "Standard Secondary Resources" below) for the primary resource, as long as:
 * <ul>
 *     <li>the content type of the primary resource can be determined from the resources's file extension</li>
 *     <li>the primary resource content type is a well known type (i.e., a variants or a fasta file)</li>
 *     <li>the secondary resources are siblings (in the same parent directory) as the primary resource</li>
 * </ul>
 * The "--suppress-resource-resolution" argument can be used to suppress this secondary
 * resource inference behavior. If the type of the primary resource is not expliclty provided, or cannot be
 * determined from the file extensions, or if the standard secondary resources for a standard primary resource
 * cannot be found in the same director as the primary, an exception will be thrown.
 * <p>
 * Alternatively, you can explicitly specify all of the resources and their content types. This is useful when
 * the resources are not in the same directory, or when the content types are not standard.
 * <p>
 * For the primary resource, if the content type is not specified on the command line (content types are
 * supplied on the command line using argument tags - see "Standard Content Types" and the examples below),
 * CreateBundle will attempt to  determine the content type of the primary resource. Content types for explicitly
 * supplied secondary resources are never inferred, and must always be supplied explicitly.
 * <p>
 * Bundle output file names must end with the suffix ".json".
 * <p>
 * <h2>Standard Content Types</h2>
 * In general, bundle content types can be any string, but many tools expect bundles to use standard, well known
 * content types that are pre-defined, such as content types for a VCF, a VCF index, a .fasta file, or a reference
 * dictionary file. The common well known content types are:
 * <h4>Standard VCF Content Types:</h4>
 * <ul>
 *  <li>"CT_VARIANT_CONTEXTS": a VCF file</li>
 *  <li>"CT_VARIANTS_INDEX": a VCF index file</li>
 * </ul>
 * <h4>Standard Reference Content Types</h4>
 * <ul>
 *  <li>"CT_HAPLOID_REFERENCE": a fasta reference file</li>
 *  <li>"CT_HAPLOID_REFERENCE_INDEX": a fasta index file</li>
 *  <li>"CT_HAPLOID_REFERENCE_DICTIONARY": a fasta dictionary file</li>
 * </ul><p>
 * <h2>Standard Secondary Resources</h2>
 * <h3>For VCFS, the standard (inferred) secondary resources are:</h3>
 * <ul>
 * <li>an index file</li>
 * </ul>
 *<h3>For references, the standard secondary resources are:</h3>
 * <ul>
 * <li>an index file</li>
 *<li>a dictionary file</li>
 * </ul>
 * <h3>Common bundle creation examples:</h3>
 * </p>
 * <h3>VCF Bundle Examples</h3>
 * <p>
 * 1) Create a resource bundle for a VCF from just the VCF, letting the tool resolve the secondary (index) resource by
 * automatically finding the sibling index file, and letting the tool determine the content types. If the sibling index
 * file cannot be found, an exception will be thrown. The resulting bundle contains the VCF and associated index.
 * <pre>
 *    CreateBundle \
 *      --primary path/to/my.vcf \
 *      --output mybundle.json
 * </pre><p>
 * The exact same bundle could be created manually by specifying both the resources and the content types explicitly:
 * <pre>
 *     CreateBundle \
 *      --primary:CT_VARIANT_CONTEXTS path/to/my.vcf \
 *      --secondary:CT_VARIANTS_INDEX path/to/my.vcf.idx \
 *      --output mybundle.json
 * </pre><p>
 * <p>
 * 2) Create a resource bundle for a VCF from just the VCF, but suppress automatic resolution of the secondary
 * resources. Let the tool determine the content type. The resulting bundle will contain only the VCF resource:
 * <pre>
 *    CreateBundle \
 *      --primary path/to/my.vcf \
 *      --suppress-resource-resolution \
 *      --output mybundle.json
 * </pre><p>
 *  3) Create a resource bundle for a VCF, but specify the VCF AND the secondary index resource explicitly (which
 *  suppresses automatic secondary resolution). This is useful when the VCF and index are not in the same directory.
 *  Let the tool determine the primary content type. The resulting bundle will contain the VCF and index resources:
 * <pre>
 *     CreateBundle \
 *      --primary path/to/my.vcf \
 *      --secondary:CT_VARIANTS_INDEX some/other/path/to/vcd.idx \
 *      --output mybundle.json
 * </pre><p>
 *  4) Create a resource bundle for a VCF, but specify the VCF AND the secondary index resource explicitly (this
 *  is useful when the VCF and index are not in the same directory), and specify the content types explicitly via
 *  command line argument tags. The resulting bundle will contain the VCF and index resources.
 * <pre>
 *     CreateBundle \
 *      --primary:CT_VARIANT_CONTEXTS path/to/my.vcf \
 *      --secondary:CT_VARIANTS_INDEX some/other/path/to/vcd.idx \
 *      --output mybundle.json
 * </pre><p>
 *<h3>Reference Bundle Examples</h3>
 *<p>
 * 1) Create a resource bundle for a reference from just the .fasta, letting the tool resolve the secondary
 * (index and dictionary) resource by automatically finding the sibling files, and determining the content types.
 * If the sibling index file cannot be found, an exception will be thrown. The resulting bundle will contain the
 * reference, index, and dictionary.
 * <pre>
 *    CreateBundle \
 *      --primary path/to/my.fasta \
 *      --output mybundle.json
 * </pre><p>
 * 2) Create a resource bundle for a reference from just the .fasta, but suppress resolution of the secondary index and
 * dictionary resources. Let the tool determine the content type. The resulting bundle will contain only the .fasta
 * resource:
 * <pre>
 *    CreateBundle \
 *      --primary path/to/my.fasta \
 *      --suppress-resource-resolution \
 *      --output mybundle.json
 * </pre><p>
 *  3) Create a resource bundle for a fasta, but specify the fasta AND the secondary index and dictionary resources
 *  explicitly (which suppresses automatic secondary resolution). Let the tool determine the content types. The
 *  resulting bundle will contain the fasta, index  and dictionary resources:
 * <pre>
 *     CreateBundle \
 *      --primary path/to/my.fasta \
 *      --secondary:CT_HAPLOID_REFERENCE_INDEX some/other/path/to/my.fai \
 *      --secondary:CT_HAPLOID_REFERENCE_DICTIONARY some/other/path/to/my.dict \
 *      --output mybundle.json
 * </pre><p>
 *  4) Create a resource bundle for a fasta, but specify the fasta, index and dictionary resources and the content
 *  types explicitly. The resulting bundle will contain the fasta, index  and dictionary resources:
 * <pre>
 *     CreateBundle \
 *      --primary:CT_HAPLOID_REFERENCE path/to/my.fasta \
 *      --secondary:CT_HAPLOID_REFERENCE_INDEX some/other/path/to/my.fai \
 *      --secondary:CT_HAPLOID_REFERENCE_DICTIONARY some/other/path/to/my.dict \
 *      --output mybundle.json
 * </pre><p>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Create a bundle (JSON) file for use with a GATK tool",
        oneLineSummary = "Create a bundle (JSON) file for use with a GATK tool",
        programGroup = OtherProgramGroup.class
)
public class CreateBundle extends CommandLineProgram {
    protected static final Logger logger = LogManager.getLogger(CreateBundle.class);

    @Argument(fullName = StandardArgumentDefinitions.PRIMARY_RESOURCE_LONG_NAME,
            shortName = StandardArgumentDefinitions.PRIMARY_RESOURCE_SHORT_NAME,
            doc="Path to the primary resource (content type will be inferred if no content type tag is specified)")
    GATKPath primaryResource;

    @Argument(fullName = StandardArgumentDefinitions.SECONDARY_RESOURCE_LONG_NAME,
            shortName = StandardArgumentDefinitions.SECONDARY_RESOURCE_SHORT_NAME,
            doc = "Path to a secondary resource for" + StandardArgumentDefinitions.PRIMARY_RESOURCE_LONG_NAME +
                    " (usually an index). If no secondary resource is specified, standard secondary resources" +
                    " for the primary bundle resource will be automatically inferred unless " +
                    SUPPRESS_SECONDARY_RESOURCE_RESOLUTION_FULL_NAME +
                    " is specified.",
            optional = true)
    GATKPath secondaryResource;

    public static final String OTHER_RESOURCE_FULL_NAME = "other-resource";
    @Argument(fullName = OTHER_RESOURCE_FULL_NAME,
            shortName = OTHER_RESOURCE_FULL_NAME,
            doc = "Path to other bundle resources for " + StandardArgumentDefinitions.PRIMARY_RESOURCE_LONG_NAME
                    + ". The content is not inferred for \"other\" (non-primary and non-secondary) resources, so a" +
                    " content type tag MUST be provided for each \"other\" resource.",
            optional = true)
    List<GATKPath> otherResources;

    public static final String SUPPRESS_SECONDARY_RESOURCE_RESOLUTION_FULL_NAME = "suppress-resource-resolution";
    @Argument(fullName = SUPPRESS_SECONDARY_RESOURCE_RESOLUTION_FULL_NAME,
            doc ="Don't attempt to resolve the primary resource's secondary resources if no secondary resource is explicitly" +
                    " provided (defaults to false) ",
            optional = true)
    boolean suppressResourceResolution = false;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Path the output bundle file (must end with the suffix " + BundleJSON.BUNDLE_EXTENSION + ")")
    GATKPath outputBundlePath;

    private enum BundleType {
        VCF,
        REFERENCE,
        CUSTOM
    }
    private BundleType outputBundleType;

    @Override
    protected String[] customCommandLineValidation() {
        if (!outputBundlePath.toString().endsWith(BundleJSON.BUNDLE_EXTENSION)){
            return new String[]{ String.format("Output bundle path must end with the suffix %s", BundleJSON.BUNDLE_EXTENSION) };
        }
        outputBundleType = determinePrimaryContentType();
        return super.customCommandLineValidation();
    }

    @Override
    protected Object doWork() {
        try (final BufferedWriter writer = Files.newBufferedWriter(outputBundlePath.toPath(), StandardOpenOption.CREATE)) {
            final Bundle bundle = switch (outputBundleType) {
                case VCF -> createVCFBundle();
                case REFERENCE -> createHaploidReferenceBundle();
                case CUSTOM -> createOtherBundle();
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
        // primary resource's file extension
        final String primaryContentTag = primaryResource.getTag();
        if (primaryContentTag != null) {
            if (primaryContentTag.equals(BundleResourceType.CT_VARIANT_CONTEXTS)) {
                bundleType = BundleType.VCF;
            } else if (primaryContentTag.equals(BundleResourceType.CT_HAPLOID_REFERENCE)) {
                bundleType = BundleType.REFERENCE;
            } else {
                logger.info(String.format("Primary input content type %s for %s not recognized. A bundle will be created using content types from the provided argument tags.",
                        primaryContentTag,
                        primaryResource));
                bundleType = BundleType.CUSTOM;
            }
        } else {
            logger.info(String.format("A content type for the primary input was not provided. Attempting to infer the primary content type from the %s extension.", primaryResource));
            bundleType = inferPrimaryContentType(primaryResource);
        }
        return bundleType;
    }

    private BundleType inferPrimaryContentType(final GATKPath primaryResource) {
        logger.info("Attempting to infer bundle content type from file extension.");
        if (FileExtensions.VCF_LIST.stream().anyMatch(ext -> primaryResource.hasExtension(ext))) {
            return BundleType.VCF;
        } else if (FileExtensions.FASTA.stream().anyMatch(ext -> primaryResource.hasExtension(ext))) {
            return BundleType.REFERENCE;
        } else {
            throw new IllegalArgumentException(String.format("Unable to infer bundle content type from file extension %s. A content type must be provided as part of the argument.", primaryResource));
        }
    }

    private Bundle createVCFBundle() {
        final Collection<BundleResource> bundleResources = new ArrayList<>();

        bundleResources.add(new IOPathResource(primaryResource, BundleResourceType.CT_VARIANT_CONTEXTS));
        if (secondaryResource != null) {
            final String secondaryContentType = secondaryResource.getTag();
            if (secondaryContentType == null) {
                throw new IllegalArgumentException(String.format("A content type for the secondary input %s must be provided.",
                        secondaryResource.getRawInputString()));
            } else {
                bundleResources.add(new IOPathResource(secondaryResource, secondaryContentType));
            }
        } else if (!suppressResourceResolution) {
            // secondary resources is null, and resource resolution suppression is off
            final Optional<GATKPath> indexPath = VariantsBundle.resolveIndex(primaryResource, GATKPath::new);
            if (indexPath.isEmpty()) {
                throw new IllegalArgumentException(
                        String.format(
                                "Could not infer an index for %s, you must either specify the index path as a secondary input on the command line or specify the %s argument.",
                                primaryResource.getRawInputString(),
                                SUPPRESS_SECONDARY_RESOURCE_RESOLUTION_FULL_NAME));
            }
            bundleResources.add(new IOPathResource(indexPath.get(), BundleResourceType.CT_VARIANTS_INDEX));
        }
        if (otherResources != null && !otherResources.isEmpty()) {
            for (final GATKPath otherInput : otherResources) {
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

    private Bundle createHaploidReferenceBundle() {
        final Collection<BundleResource> bundleResources = new ArrayList<>();
        if (secondaryResource == null && (otherResources == null || otherResources.isEmpty())) {
            // infer dictionary and index
            return HaploidReferenceResolver.referenceBundleFromFastaPath(primaryResource, GATKPath::new);
        }

        bundleResources.add(new IOPathResource(primaryResource, primaryResource.getTag()));
        if (secondaryResource != null) {
            final String secondaryContentType = secondaryResource.getTag();
            if (secondaryContentType == null) {
                throw new IllegalArgumentException(String.format("A content type for the secondary input %s must be provided.",
                        secondaryResource.getRawInputString()));
            } else {
                bundleResources.add(new IOPathResource(secondaryResource, secondaryContentType));
            }
        }
        if (otherResources != null && !otherResources.isEmpty()) {
            for (final GATKPath otherInput : otherResources) {
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
        return new Bundle(primaryResource.getTag(), bundleResources);
    }

    private Bundle createOtherBundle() {
        final Collection<BundleResource> bundleResources = new ArrayList<>();
        bundleResources.add(new IOPathResource(primaryResource, primaryResource.getTag()));
        if (secondaryResource != null) {
            final String secondaryContentType = secondaryResource.getTag();
            if (secondaryContentType == null) {
                throw new IllegalArgumentException(String.format("A content type for the secondary input %s must be provided.",
                        secondaryResource.getRawInputString()));
            } else {
                bundleResources.add(new IOPathResource(secondaryResource, secondaryContentType));
            }
        }
        if (otherResources != null && !otherResources.isEmpty()) {
            for (final GATKPath otherInput : otherResources) {
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
        return new Bundle(primaryResource.getTag(), bundleResources);
    }
}
