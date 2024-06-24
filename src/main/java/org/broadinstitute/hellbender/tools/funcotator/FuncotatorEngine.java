package org.broadinstitute.hellbender.tools.funcotator;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.util.ParsingUtils;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.broadinstitute.hellbender.engine.filters.VariantFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.compositeoutput.CompositeOutputRenderer;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.genelistoutput.GeneListOutputRenderer;
import org.broadinstitute.hellbender.tools.funcotator.mafOutput.MafOutputRenderer;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadata;
import org.broadinstitute.hellbender.tools.funcotator.simpletsvoutput.SimpleTsvOutputRenderer;
import org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer;
import org.broadinstitute.hellbender.transformers.VariantTransformer;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Class that performs functional annotation of variants.
 *
 * Requires a set of data sources ({@link DataSourceFuncotationFactory}) from which to create {@link Funcotation}s.
 */
public final class FuncotatorEngine implements AutoCloseable {

    /** Obligatory logger. */
    private static final Logger logger = LogManager.getLogger(FuncotatorEngine.class);
    private static final String SIMPLE_TSV_SEG_FILE_CONFIG = "org/broadinstitute/hellbender/tools/funcotator/simple_funcotator_seg_file.config";

    @VisibleForTesting
    static final String GENE_LIST_FILE_SUFFIX = ".gene_list.txt";

    /**
     * Information about what kinds of {@link Funcotation}s are going to be created by this {@link FuncotatorEngine}.
     */
    private final FuncotationMetadata inputMetadata;

    /**
     * The {@link DataSourceFuncotationFactory} that will create {@link Funcotation}s for this {@link FuncotatorEngine}.
     */
    private final List<DataSourceFuncotationFactory> dataSourceFactories;

    /**
     * The arguments given to the instance of the {@link GATKTool} running this {@link FuncotatorEngine}.
     */
    private final BaseFuncotatorArgumentCollection funcotatorArgs;

    /**
     * The {@link SAMSequenceDictionary} for the driving variants (i.e. the input variant file).
     */
    private final SAMSequenceDictionary sequenceDictionaryForDrivingVariants;

    /**
     * Whether the input variant contigs must be converted to hg19.
     * This is only the case when the input reference is b37 AND when
     * the reference version is hg19 (i.e. {@link FuncotatorVariantArgumentCollection#referenceVersion} == {@link FuncotatorVariantArgumentCollection.FuncotatorReferenceVersion#hg19}).
     */
    private final boolean mustConvertInputContigsToHg19;

    /**
     * Whether the output variant contigs must be converted back to B37 from hg19 before being returned.
     * (NOTE: this means that the output contigs will continue to use B37 contig names even if internally we converted them to hg19)
     */
    private boolean mustRevertVariantContigsFromHg19ToB37 = false;

    /**
     * Whether this {@link FuncotatorEngine} has only produced annotations on variants that have been labeled by the
     * {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationFactory} as {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification#IGR}.
     */
    private boolean onlyProducedIGRs = true;

    /**
     * Create a {@link FuncotatorEngine} using the given {@code metadata} and {@code funcotationFactories} representing
     * the kinds of {@link Funcotation}s to be created and the data sources from which they should be created,
     * respectively.
     * @param metadata {@link FuncotationMetadata} containing information on the kinds of {@link Funcotation}s this {@link FuncotatorEngine} will create to represent the input file.
     *          For example, this could be based on the existing annotations for an input VCF.
     * @param funcotationFactories A {@link List<DataSourceFuncotationFactory>} which can create the desired {@link Funcotation}s.
     */
    public FuncotatorEngine(final BaseFuncotatorArgumentCollection funcotatorArgs,
                            final SAMSequenceDictionary sequenceDictionaryForDrivingVariants,
                            final FuncotationMetadata metadata,
                            final List<DataSourceFuncotationFactory> funcotationFactories) {

        this.sequenceDictionaryForDrivingVariants = sequenceDictionaryForDrivingVariants;
        this.funcotatorArgs = funcotatorArgs;
        inputMetadata = metadata;

        dataSourceFactories = funcotationFactories;
        // Note: The dataSourceFactories must be sorted to ensure that as we iterate through them
        // to create funcotations, the inherent dependencies between different funcotation types are preserved.
        // For example, most FuncotationFactories require that a GencodeFuncotation is present before they can
        // create their annotations.   This sorting enables such dependencies.
        dataSourceFactories.sort(DataSourceUtils::datasourceComparator);

        // Determine whether we have to convert given variants from B37 to HG19:
        mustConvertInputContigsToHg19 = determineReferenceAndDatasourceCompatibility();

        // Read in the custom variant classification order file here so that it can be shared across all engines:
        if (funcotatorArgs.customVariantClassificationOrderFile != null) {
            FuncotatorUtils.setVariantClassificationCustomSeverity(funcotatorArgs.customVariantClassificationOrderFile);
        }
    }

    /**
     * @return An unmodifiable {@link List<DataSourceFuncotationFactory>} being used by this {@link FuncotatorEngine} to create {@link Funcotation}s.
     */
    public List<DataSourceFuncotationFactory> getFuncotationFactories() {
        return Collections.unmodifiableList(dataSourceFactories);
    }

    /**
     * Creates a {@link FuncotationMap} for the given {@code variantContext}.
     *
     * @param variantContext   {@link VariantContext} to annotate.  Never {@code null}.
     * @param referenceContext {@link ReferenceContext} corresponding to the given {@code variantContext}.  Never {@code null}.
     * @param featureContext {@link FeatureContext} corresponding to the given {@code variantContext}.  Never {@code null}.
     * @return an instance of FuncotationMap that maps transcript IDs to lists of funcotations for the given variantContext context.
     */
    public FuncotationMap createFuncotationMapForVariant(final VariantContext variantContext,
                                                         final ReferenceContext referenceContext,
                                                         final FeatureContext featureContext) {

        Utils.nonNull(variantContext);
        Utils.nonNull(referenceContext);
        Utils.nonNull(featureContext);

        //==============================================================================================================
        // First create only the transcript (Gencode) funcotations:

        if (retrieveGencodeFuncotationFactoryStream().count() > 1) {
            logger.warn("Attempting to annotate with more than one GENCODE datasource.  If these have overlapping transcript IDs, errors may occur.");
        }

        final List<GencodeFuncotation> transcriptFuncotations = retrieveGencodeFuncotationFactoryStream()
                .map(gf -> gf.createFuncotations(variantContext, referenceContext, featureContext))
                .flatMap(List::stream)
                .map(f -> {
                        final GencodeFuncotation gf = (GencodeFuncotation) f;
                        if (onlyProducedIGRs && (gf.getVariantClassification() != GencodeFuncotation.VariantClassification.IGR)) {
                            onlyProducedIGRs = false;
                        }
                        return gf;
                    }
                )
                .collect(Collectors.toList());

        //==============================================================================================================
        // Create the funcotations for non-Gencode data sources:

        // Create a place to keep our funcotations:
        final FuncotationMap funcotationMap = FuncotationMap.createFromGencodeFuncotations(transcriptFuncotations);

        // Perform the rest of the annotation.  Note that this code manually excludes the Gencode Funcotations.
        for (final DataSourceFuncotationFactory funcotationFactory : dataSourceFactories ) {

            // Note that this guarantees that we do not add GencodeFuncotations a second time.
            if (!funcotationFactory.getType().equals(FuncotatorArgumentDefinitions.DataSourceType.GENCODE)) {
                final List<String> txIds = funcotationMap.getTranscriptList();

                for (final String txId: txIds) {
                    funcotationMap.add(txId, funcotationFactory.createFuncotations(variantContext, referenceContext,
                            featureContext, funcotationMap.getGencodeFuncotations(txId)));
                }
            }
        }

        //==============================================================================================================
        // Create the funcotations for the input and add to all txID mappings.

        final List<String> txIds = funcotationMap.getTranscriptList();

        for (final String txId: txIds) {
            funcotationMap.add(txId, FuncotatorUtils.createFuncotations(variantContext, inputMetadata, FuncotatorConstants.DATASOURCE_NAME_FOR_INPUT_VCFS));
        }

        return funcotationMap;
    }

    /**
     * Creates a {@link FuncotationMap} for the given {@code variantContext} using the datasources initialized with this
     *  engine.
     *
     * @param segmentAsVariantContext   {@link VariantContext} to annotate.  Never {@code null}.
     * @param referenceContext {@link ReferenceContext} corresponding to the given {@code variantContext}.  Never {@code null}.
     * @param featureContext {@link FeatureContext} corresponding to the given {@code variantContext}.  Never {@code null}.
     * @return an instance of FuncotationMap that maps transcript IDs to lists of funcotations for the given variantContext context.
     */
     FuncotationMap createFuncotationMapForSegment(final VariantContext segmentAsVariantContext,
                                                         final ReferenceContext referenceContext,
                                                         final FeatureContext featureContext) {

        Utils.nonNull(segmentAsVariantContext);
        Utils.nonNull(referenceContext);
        Utils.nonNull(featureContext);

        //==============================================================================================================
        // We do not need to treat the Gencode funcotations as a special case

        if (retrieveGencodeFuncotationFactoryStream().count() > 1) {
            logger.error("Attempting to annotate with more than one GENCODE datasource.  If these have overlapping transcript IDs, errors may occur, so it has been disallowed.");
            throw new UserException.BadInput("Attempting to funcotate segments with more than one GENCODE datasource.  This is currently not supported.  Please post to the forum if you would like to see support for this.");
        }

        final List<Funcotation> funcotations = dataSourceFactories.stream()
                .filter(DataSourceFuncotationFactory::isSupportingSegmentFuncotation)
                .map(ff -> ff.createFuncotations(segmentAsVariantContext, referenceContext,
                        featureContext))
                .flatMap(List::stream)
                .collect(Collectors.toList());

        // Create a place to keep our funcotations:
        return FuncotationMap.createNoTranscriptInfo(funcotations);
    }

    /**
     * Create an output renderer for the data created by this instance of {@link FuncotatorEngine}.
     * @param annotationDefaultsMap {@link LinkedHashMap<String,String>} of annotation names and their default values.
     * @param annotationOverridesMap {@link LinkedHashMap<String,String>} of annotation names and the values for these fields overridden by the user.
     * @param headerForVariants {@link VCFHeader} for the input VCF file containing the variants to annotate.
     * @param defaultToolVcfHeaderLines {@link Set<VCFHeaderLine>} containing the default {@link VCFHeaderLine}s for the given {@code gatkToolInstance}.
     * @param gatkToolInstance {@link GATKTool} instance from which we will be using this {@link FuncotatorEngine}.
     * @return The requested {@link OutputRenderer} based on the given {@code funcotatorArgs}.
     */
    OutputRenderer createOutputRenderer(final LinkedHashMap<String, String> annotationDefaultsMap,
                                        final LinkedHashMap<String, String> annotationOverridesMap,
                                        final VCFHeader headerForVariants,
                                        final Set<VCFHeaderLine> defaultToolVcfHeaderLines,
                                        final GATKTool gatkToolInstance) {

        final OutputRenderer outputRenderer;

        // Determine which annotations are accounted for (by the funcotation factories) and which are not.
        final LinkedHashMap<String, String> unaccountedForDefaultAnnotations = getUnaccountedForAnnotations(  getFuncotationFactories(), annotationDefaultsMap );
        final LinkedHashMap<String, String> unaccountedForOverrideAnnotations = getUnaccountedForAnnotations( getFuncotationFactories(), annotationOverridesMap );

        // Set up our output renderer:
        switch (funcotatorArgs.outputFormatType) {
            case MAF:
                outputRenderer = new MafOutputRenderer(
                        funcotatorArgs.outputFile.toPath(),
                        getFuncotationFactories(),
                        headerForVariants,
                        unaccountedForDefaultAnnotations,
                        unaccountedForOverrideAnnotations,
                        defaultToolVcfHeaderLines.stream().map(Object::toString).collect(Collectors.toCollection(LinkedHashSet::new)),
                        funcotatorArgs.referenceVersion,
                        funcotatorArgs.excludedFields,
                        gatkToolInstance.getVersion()
                );
                break;

            case VCF:
                outputRenderer = new VcfOutputRenderer(
                        gatkToolInstance.createVCFWriter(funcotatorArgs.outputFile),
                        getFuncotationFactories(),
                        headerForVariants,
                        unaccountedForDefaultAnnotations,
                        unaccountedForOverrideAnnotations,
                        defaultToolVcfHeaderLines,
                        funcotatorArgs.excludedFields,
                        gatkToolInstance.getVersion()
                );
                break;

            case SEG:
                // Create an output renderer that will actually write multiple files.
                outputRenderer =  new CompositeOutputRenderer(
                            Arrays.asList(
                                    SimpleTsvOutputRenderer.createFromResource(funcotatorArgs.outputFile.toPath(),
                                        unaccountedForDefaultAnnotations,
                                        unaccountedForOverrideAnnotations, funcotatorArgs.excludedFields,
                                        Paths.get(SIMPLE_TSV_SEG_FILE_CONFIG),
                                        gatkToolInstance.getVersion(), true),

                                    new GeneListOutputRenderer(new File(funcotatorArgs.outputFile.getAbsolutePath() + GENE_LIST_FILE_SUFFIX).toPath(),
                                        unaccountedForDefaultAnnotations, unaccountedForOverrideAnnotations,
                                        funcotatorArgs.excludedFields, gatkToolInstance.getVersion(), funcotatorArgs.minNumBasesForValidSegment)
                        ), gatkToolInstance.getVersion());
                break;
            default:
                throw new GATKException("Unsupported output format type specified: " + funcotatorArgs.outputFormatType.toString());
        }

        return outputRenderer;
    }

    /**
     * @return A {@link VariantFilter} that will ignore any variants that have been filtered (if the user requested that the filter is turned on).  Otherwise returns a no-op filter.
     */
    public static VariantFilter makeVariantFilter(final boolean isRemovingFilteredVariants) {
        // Ignore variants that have been filtered if the requested:
        return isRemovingFilteredVariants ?
                VariantFilterLibrary.PASSES_FILTERS :
                VariantFilterLibrary.ALLOW_ALL_VARIANTS;
    }

    /**
     * Create a new {@link VariantContext} which will match the given Reference if there is a mismatch for input between the B37 reference and the HG19 reference.
     * @param variant A {@link VariantContext} object containing the variant to convert.
     * @return A {@link VariantContext} whose contig has been transformed to HG19 if requested by the user.  Otherwise, an identical variant.
     */
    private VariantContext getCorrectVariantContextForReference(final VariantContext variant) {
        if ( mustConvertInputContigsToHg19 ) {
            final VariantContextBuilder vcb = new VariantContextBuilder(variant);
            vcb.chr(FuncotatorUtils.convertB37ContigToHg19Contig(variant.getContig()));
            return vcb.make();
        }
        else {
            return variant;
        }
    }

    /**
     * Create a new {@link VariantContext} which will match the given Reference if there is a mismatch for input between the B37 reference and the HG19 reference.
     * @param variant A {@link VariantContext} object containing the variant to convert.
     * @return A {@link VariantContext} whose contig has been transformed to HG19 if requested by the user.  Otherwise, an identical variant.
     */
    VariantContext getCorrectVariantContextForOutput(final VariantContext variant) {
        if ( mustRevertVariantContigsFromHg19ToB37 ) {
            final VariantContextBuilder vcb = new VariantContextBuilder(variant);
            vcb.chr(FuncotatorUtils.convertHG19ContigToB37Contig(variant.getContig()));
            return vcb.make();
        }
        else {
            return variant;
        }
    }

    /**
     * @return The default {@link VariantTransformer} which will automatically convert from the B37 reference standard to the HG19 reference standard for contig names.
     */
    public VariantTransformer getDefaultVariantTransformer() {
        return variantContext -> getCorrectVariantContextForReference(variantContext);
    }

    /**
     * Shutdown the engine.  Closes all datasource factories.
     */
    public void close() {
        for ( final DataSourceFuncotationFactory factory : dataSourceFactories ) {
            if ( factory != null ) {
                factory.close();
            }
        }
    }

    /**
     * Processes the given {@link Set} into a list of transcript IDs.
     * This is necessary because the command-line input argument is overloaded to be either a file containing transcript
     * IDs (1 per line) OR as a list of transcript IDs.
     * @param rawTranscriptSet {@link Set} of {@link String}s from which to create a list of Transcript IDs.  If of size 1, will try to open as a file.
     * @return A {@link Set} of {@link String} contianing Transcript IDs in which the user is interested.
     */
    public static Set<String> processTranscriptList(final Set<String> rawTranscriptSet) {
        if ( rawTranscriptSet.size() == 1 ) {
            final String filePathString = rawTranscriptSet.iterator().next();
            try ( final BufferedReader bufferedReader = Files.newBufferedReader(IOUtils.getPath(filePathString)) ) {
                logger.info("Opened transcript file: " + filePathString);

                // Create a place to put our output:
                final Set<String> transcriptIdSet = new HashSet<>();

                String line = bufferedReader.readLine();
                while ( line != null ) {
                    logger.info("    Adding transcript ID to transcript set: " + line);
                    transcriptIdSet.add(line);
                    line = bufferedReader.readLine();
                }
                logger.info("Transcript parsing complete.");

                return transcriptIdSet;
            }
            catch ( final IOException ex ) {
                logger.warn("Could not open transcript selection list as a file.  Using it as a singleton list of transcript IDs: [" + filePathString + "]");
                return rawTranscriptSet;
            }
        }
        else {
            return rawTranscriptSet;
        }
    }


    /**
     * Gets the correct {@link ReferenceContext} for the {@code variant} being processed based on if the B37->HG19 conversion is required.
     * @param variant {@link VariantContext} to check for B37/HG19 compliance.
     * @param referenceContext {@link ReferenceContext} on which the given {@code variant} was originally based before the variant transformation.
     * @return A {@link ReferenceContext} that is guaranteed to match the given {@code variant} for HG19/B37 compliance.
     */
    public ReferenceContext getCorrectReferenceContext(final VariantContext variant, final ReferenceContext referenceContext) {

        final ReferenceContext correctReferenceContext;

        // Check to see if we need to revert the ReferenceContext's interval to the original variant interval
        // (This would only happen in the case where we were given b37 variants with hg19 data sources):
        if ( mustConvertInputContigsToHg19 ) {

            // Convert our contig back to B37 here so it matches the variant:
            final SimpleInterval interval = new SimpleInterval(
                    FuncotatorUtils.convertHG19ContigToB37Contig(variant.getContig()), variant.getStart(), variant.getEnd()
            );

            correctReferenceContext = new ReferenceContext(referenceContext, interval);
        }
        else {
            correctReferenceContext = referenceContext;
        }

        return correctReferenceContext;
    }

    /**
     * Returns whether this {@link FuncotatorEngine} has only produced annotations on variants that have been labeled by the
     * {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationFactory} as {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification#IGR}.
     * @return {@code true} IFF this {@link FuncotatorEngine} has only produced IGR annotations.
     */
    public boolean onlyProducedIGRs() {
        return onlyProducedIGRs;
    }

    // =================================================================================================================

    /**
     * Creates a {@link LinkedHashMap} of annotations in the given {@code annotationMap} that do not occur in the given {@code dataSourceFactories}.
     * @param dataSourceFactories {@link List} of {@link DataSourceFuncotationFactory} to check for whether each annotation in the {@code annotationMap} is handled.
     * @param annotationMap {@link Map} (of ANNOTATION_NAME : ANNOTATION_VALUE) to check
     * @return A {@link LinkedHashMap} of annotations in the given {@code annotationMap} that do not occur in the given {@code dataSourceFactories}.
     */
    public LinkedHashMap<String, String> getUnaccountedForAnnotations( final List<DataSourceFuncotationFactory> dataSourceFactories,
                                                                        final Map<String, String> annotationMap ) {
        final LinkedHashMap<String, String> outAnnotations = new LinkedHashMap<>();

        // Check each field in each factory:
        for ( final String field : annotationMap.keySet() ) {
            boolean accountedFor = false;
            for ( final DataSourceFuncotationFactory funcotationFactory : dataSourceFactories ) {

                if ( funcotationFactory.getSupportedFuncotationFields().contains(field) ) {
                    accountedFor = true;
                    break;
                }
            }
            if ( !accountedFor ) {
                outAnnotations.put(field, annotationMap.get(field));
            }
        }

        return outAnnotations;
    }

    /**
     * Split each element of the given {@link List} into a key and value.
     * Assumes each element of the given {@link List} is formatted as follows:
     *     KEY:VALUE
     * @param annotationArgs {@link List} of strings formatted KEY:VALUE to turn into a {@link Map}.
     * @return A {@link LinkedHashMap} of KEY:VALUE pairs corresponding to entries in the given list.
     */
    public static LinkedHashMap<String, String> splitAnnotationArgsIntoMap( final List<String> annotationArgs ) {

        final LinkedHashMap<String, String> annotationMap = new LinkedHashMap<>();

        for ( final String s : annotationArgs ) {
            final List<String> keyVal = ParsingUtils.split(s, FuncotatorArgumentDefinitions.MAP_NAME_VALUE_DELIMITER);
            if ( keyVal.size() != 2) {
                throw new UserException.BadInput( "Argument annotation incorrectly formatted: " + s );
            }

            annotationMap.put( keyVal.get(0), keyVal.get(1) );
        }

        return annotationMap;
    }

    private boolean determineReferenceAndDatasourceCompatibility() {

        boolean mustConvertInputContigsToHg19 = false;

        // Do individual checks here so we can have a helpful log message for each case:
        if ( funcotatorArgs.forceB37ToHg19ContigNameConversion ) {
            logger.info("Forcing B37 -> HG19 Variant conversion.");
            mustConvertInputContigsToHg19 = true;
        }
        else if ( funcotatorArgs.referenceVersion.equals(BaseFuncotatorArgumentCollection.FuncotatorReferenceVersionHg19) &&
                        FuncotatorUtils.isSequenceDictionaryUsingB37Reference(sequenceDictionaryForDrivingVariants) ) {
            logger.info("VCF sequence dictionary detected as B37 in HG19 annotation mode.  Performing conversion. (NOTE: the output VCF will still be B37)");
            mustConvertInputContigsToHg19 = true;
        }
        else {
            logger.info("Using given VCF and Reference.  No conversion required.");
        }

        if (mustConvertInputContigsToHg19) {
            // NOTE AND WARNING:
            // hg19 is from ucsc. b37 is from the genome reference consortium.
            // ucsc decided the grc version had some bad data in it, so they blocked out some of the bases, aka "masked" them
            // so the lengths of the contigs are the same, the bases are just _slightly_ different.
            // ALSO, the contig naming convention is different between hg19 and hg38:
            //      hg19 uses contigs of the form "chr1"
            //      b37 uses contigs of the form  "1"
            // This naming convention difference causes a LOT of issues and was a bad idea.

            logger.warn("WARNING: You are using B37 as a reference.  " +
                    "Funcotator will convert your variants to GRCh37, and this will be fine in the vast majority of cases.  " +
                    "There MAY be some errors (e.g. in the Y chromosome, but possibly in other places as well) due to changes between the two references.");
        }

        // Record whether we need to revert the contigs back to B37 after annotation:
        if (mustConvertInputContigsToHg19 && FuncotatorUtils.isSequenceDictionaryUsingB37Reference(sequenceDictionaryForDrivingVariants)) {
            this.mustRevertVariantContigsFromHg19ToB37 = true;
        }

        return mustConvertInputContigsToHg19;
    }

    private Stream<DataSourceFuncotationFactory> retrieveGencodeFuncotationFactoryStream() {
        return dataSourceFactories.stream()
                .filter(f -> f.getType().equals(FuncotatorArgumentDefinitions.DataSourceType.GENCODE));
    }
}
