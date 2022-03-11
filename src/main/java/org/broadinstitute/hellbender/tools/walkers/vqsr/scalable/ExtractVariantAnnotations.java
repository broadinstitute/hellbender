package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.LabeledVariantAnnotationsData;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.VariantType;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.io.File;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

/**
 * TODO
 *
 * DEVELOPER NOTE: See documentation in {@link LabeledVariantAnnotationsWalker}.
 */
@CommandLineProgramProperties(
        // TODO
        summary = "",
        oneLineSummary = "",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
public final class ExtractVariantAnnotations extends LabeledVariantAnnotationsWalker {

    /**
     * TODO
     */
    @Argument(
            fullName = "maximum-number-of-unlabeled-variants",
            doc = "", // TODO
            minValue = 0)
    private int maximumNumberOfUnlabeledVariants = 10_000_000;

    @Argument(
            fullName = "reservoir-sampling-random-seed",
            doc = "") // TODO
    private int reservoirSamplingRandomSeed = 0;

    private RandomGenerator rng;
    private LabeledVariantAnnotationsData unlabeledDataReservoir; // will not be sorted in genomic order
    private int unlabeledIndex = 0;

    @Override
    public void afterOnTraversalStart() {
        rng = RandomGeneratorFactory.createRandomGenerator(new Random(reservoirSamplingRandomSeed));
        unlabeledDataReservoir = maximumNumberOfUnlabeledVariants == 0
                ? null
                : new LabeledVariantAnnotationsData(annotationNames, resourceLabels, useASAnnotations, maximumNumberOfUnlabeledVariants);
    }

    @Override
    protected void nthPassApply(final VariantContext variant,
                                final ReadsContext readsContext,
                                final ReferenceContext referenceContext,
                                final FeatureContext featureContext,
                                final int n) {
        if (n == 0) {
            final List<Triple<List<Allele>, VariantType, TreeSet<String>>> metadata = extractVariantMetadata(
                    variant, featureContext, unlabeledDataReservoir != null);
            final boolean isVariantExtracted = !metadata.isEmpty();
            if (isVariantExtracted) {
                final boolean isUnlabeled = metadata.stream().map(Triple::getRight).allMatch(Set::isEmpty);
                if (!isUnlabeled) {
                    addExtractedVariantToData(data, variant, metadata);
                    writeExtractedVariantToVCF(variant, metadata);
                } else {
                    // Algorithm R for reservoir sampling: https://en.wikipedia.org/wiki/Reservoir_sampling#Simple_algorithm
                    if (unlabeledIndex < maximumNumberOfUnlabeledVariants) {
                        addExtractedVariantToData(unlabeledDataReservoir, variant, metadata);
                    } else {
                        final int j = rng.nextInt(unlabeledIndex);
                        if (j < maximumNumberOfUnlabeledVariants) {
                            setExtractedVariantInData(unlabeledDataReservoir, variant, metadata, j);
                        }
                    }
                    unlabeledIndex++;
                }
            }
        }
    }

    @Override
    protected void afterNthPass(final int n) {
        if (n == 0) {
            writeAnnotationsToHDF5AndClearData();
            if (unlabeledDataReservoir != null) {
                writeUnlabeledAnnotationsToHDF5AndClearData();
                // TODO write extracted unlabeled variants to VCF, which can be used to mark extraction in scoring step
            }
            if (vcfWriter != null) {
                vcfWriter.close();
            }
        }
    }

    @Override
    public Object onTraversalSuccess() {

        // TODO FAIL if annotations that are all NaN
        // TODO WARN if annotations that have zero variance

        logger.info(String.format("%s complete.", getClass().getSimpleName()));

        return null;
    }

    private static void setExtractedVariantInData(final LabeledVariantAnnotationsData data,
                                                  final VariantContext variant,
                                                  final List<Triple<List<Allele>, VariantType, TreeSet<String>>> metadata,
                                                  final int index) {
        data.set(index, variant,
                metadata.stream().map(Triple::getLeft).collect(Collectors.toList()),
                metadata.stream().map(Triple::getMiddle).collect(Collectors.toList()),
                metadata.stream().map(Triple::getRight).collect(Collectors.toList()));
    }

    // TODO clean up
    private void writeUnlabeledAnnotationsToHDF5AndClearData() {
        final File outputUnlabeledAnnotationsFile = new File(outputPrefix + ".unlabeled" + ANNOTATIONS_HDF5_SUFFIX);
        if (unlabeledDataReservoir.size() == 0) {
            throw new GATKException("No unlabeled variants were present in the input VCF.");
        }
        for (final VariantType variantType : variantTypesToExtract) {
            logger.info(String.format("Extracted unlabeled annotations for %d variants of type %s.",
                    unlabeledDataReservoir.getVariantTypeFlat().stream().mapToInt(t -> t == variantType ? 1 : 0).sum(), variantType));
        }
        logger.info(String.format("Extracted unlabeled annotations for %s total variants.", unlabeledDataReservoir.size()));

        logger.info("Writing unlabeled annotations...");
        // TODO coordinate sort
        unlabeledDataReservoir.writeHDF5(outputUnlabeledAnnotationsFile, omitAllelesInHDF5);
        logger.info(String.format("Unlabeled annotations and metadata written to %s.", outputUnlabeledAnnotationsFile.getAbsolutePath()));

        unlabeledDataReservoir.clear();
    }
}