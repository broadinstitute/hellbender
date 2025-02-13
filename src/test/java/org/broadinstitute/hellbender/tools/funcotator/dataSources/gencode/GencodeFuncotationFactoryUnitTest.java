package org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode;

import com.google.common.collect.Sets;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.testutils.FuncotatorReferenceTestUtils;
import org.broadinstitute.hellbender.tools.funcotator.*;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.gtf.*;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.FuncotatorTestUtils;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Unit test class for the {@link GencodeFuncotationFactory} class.
 * Created by jonn on 9/1/17.
 */
public class GencodeFuncotationFactoryUnitTest extends GATKBaseTest {
    // Gencode v19
    private static final String CNTN4_GENCODE_ANNOTATIONS_FILE_NAME = toolsTestDir + "funcotator/gencode.v19.CNTN4.annotation.gtf";
    private static final String CNTN4_GENCODE_TRANSCRIPT_FASTA_FILE = toolsTestDir + "funcotator/gencode.v19.CNTN4.pc_transcripts.fasta";

    private static final String GENCODE_HG19_BIG_GTF_FILE = FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER + "gencode/hg19/gencode.v19.testVariantSubset.gtf";
    private static final String GENCODE_HG19_BIG_FASTA_FILE = FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER + "gencode/hg19/gencode.v19.testVariantSubset.pc_transcripts.fa";

    //==================================================================================================================
    // Multi-Test Static Variables:

    private static final double doubleEqualsEpsilon = 0.000001;

    private static final FeatureReader<GencodeGtfFeature> gencodeHg19FeatureReader;
    private static final FeatureReader<GencodeGtfFeature> gencodeHg19TP53FeatureReader;
    private static final FeatureReader<GencodeGtfFeature> gencodeHg38FeatureReader;
    private static final FeatureReader<GencodeGtfFeature> muc16NonBasicFeatureReader;

    private static final ReferenceDataSource refDataSourceHg19Ch19;
    private static final ReferenceDataSource refDataSourceHg19Ch17;
    private static final ReferenceDataSource refDataSourceHg19Ch3;
    private static final ReferenceDataSource refDataSourceHg38;

    private static GencodeFuncotationFactory testMuc16SnpCreateFuncotationsFuncotationFactory;

    private static final List<AutoCloseable> autoCloseableList = new ArrayList<>();

    // Initialization of static variables:
    static {
        muc16NonBasicFeatureReader = AbstractFeatureReader.getFeatureReader(FuncotatorTestConstants.MUC16_GENCODE_NON_BASIC_ANNOTATIONS_FILE_NAME, new GencodeGtfCodec());
        gencodeHg19FeatureReader = AbstractFeatureReader.getFeatureReader(FuncotatorTestConstants.GENCODE_DATA_SOURCE_GTF_PATH_HG19, new GencodeGtfCodec());
        gencodeHg19TP53FeatureReader = AbstractFeatureReader.getFeatureReader(FuncotatorTestConstants.TP53_REVERSE_ORDER_GENCODE_ANNOTATIONS_FILE_NAME, new GencodeGtfCodec());
        gencodeHg38FeatureReader = AbstractFeatureReader.getFeatureReader(FuncotatorTestConstants.GENCODE_DATA_SOURCE_GTF_PATH_HG38, new GencodeGtfCodec());
        refDataSourceHg19Ch19 = ReferenceDataSource.of(IOUtils.getPath(FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref()));
        refDataSourceHg19Ch17 = ReferenceDataSource.of(IOUtils.getPath(FuncotatorReferenceTestUtils.retrieveHg19Chr17Ref()));
        refDataSourceHg19Ch3 = ReferenceDataSource.of( IOUtils.getPath(FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref()) );
        refDataSourceHg38 = ReferenceDataSource.of( IOUtils.getPath(hg38Reference) );

        // Gets cleaned up in `cleanupAfterTests()`
        // NOTE: This is initialized here to save time in testing.
        testMuc16SnpCreateFuncotationsFuncotationFactory = new GencodeFuncotationFactory(
                IOUtils.getPath(FuncotatorTestConstants.GENCODE_TRANSCRIPT_FASTA_FILE_NAME),
                "VERSION",
                GencodeFuncotationFactory.DEFAULT_NAME,
                FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_DEFAULT_VALUE,
                new HashSet<>(),
                new LinkedHashMap<>(),
                createFeatureInputForMuc16Ds(GencodeFuncotationFactory.DEFAULT_NAME),
                "TEST");

        // Add all to the closeable list:
        autoCloseableList.add( muc16NonBasicFeatureReader );
        autoCloseableList.add( gencodeHg19FeatureReader );
        autoCloseableList.add( gencodeHg19TP53FeatureReader );
        autoCloseableList.add( gencodeHg38FeatureReader );
        autoCloseableList.add( refDataSourceHg19Ch19 );
        autoCloseableList.add( refDataSourceHg19Ch17 );
        autoCloseableList.add( refDataSourceHg19Ch3 );
        autoCloseableList.add( refDataSourceHg38 );
        autoCloseableList.add( testMuc16SnpCreateFuncotationsFuncotationFactory );
    }

    //==================================================================================================================
    // Setup and Breakdown Methods with Hooks:

    @BeforeClass
    public static void setupBeforeTests() {
        System.out.println("Setting up before tests...");
    }

    @AfterClass
    public static void cleanupAfterTests() {
        System.out.println("Cleaning up after tests...");

        for ( final AutoCloseable cl : autoCloseableList ) {

            try {
                cl.close();
            }
            catch ( final Exception ex ) {
                throw new GATKException("Could not close " + cl.toString(), ex);
            }
        }
    }

    //==================================================================================================================
    // Helper Methods:

    private static Object[] helpProvideForTestGetTranscriptEndPaddingBases( final String refAllele,
                                                                            final String altAllele,
                                                                            final int startPos,
                                                                            final List<? extends Locatable> exonPositionList,
                                                                            final String refSequence,
                                                                            final String expected) {
        final String contig = "TEST";

        final VariantContext vc = new VariantContextBuilder(
                "TEST_SOURCE",
                contig,
                startPos,
                startPos + refAllele.length() - 1,
                Arrays.asList(Allele.create(refAllele, true), Allele.create(altAllele, false))
        ).make();

        return new Object[] {
            vc,
            vc.getAlleles().get(1),
            exonPositionList,
            FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(refSequence, vc.getContig(),  vc.getEnd(),  vc.getEnd()),
            expected
        };
    }

    private static GencodeGtfExonFeature getExonForVariant( final CloseableTribbleIterator<GencodeGtfFeature> gtfFeatureIterator,
                                                            final SimpleInterval variantLocation ) {

        final Optional<GencodeGtfExonFeature> exonOption = gtfFeatureIterator.stream()
                                                    .map( f -> ((GencodeGtfGeneFeature) f))
                                                    .flatMap(g -> g.getTranscripts().stream())
                                                    .flatMap(t -> t.getExons().stream())
                                                    .filter( x -> x.overlaps(variantLocation) )
                                                    .findFirst();

        return exonOption.orElseThrow( () -> new GATKException("Could not get an exon associated with this variant's interval: " + variantLocation.toString()) );
    }

    private static GencodeGtfExonFeature getExonForVariant( final GencodeGtfGeneFeature gene,
                                                            final SimpleInterval variantLocation ) {

        final Optional<GencodeGtfExonFeature> exonOption = gene.getTranscripts().stream()
                                                                .flatMap(t -> t.getExons().stream())
                                                                .filter( x -> x.overlaps(variantLocation) )
                                                                .findFirst();

        return exonOption.orElseThrow( () -> new GATKException("Could not get an exon associated with this variant's interval: " + variantLocation.toString()) );
    }

    private static GencodeGtfTranscriptFeature getMuc16Transcript(final GencodeGtfGeneFeature gene) {
        final Optional<GencodeGtfTranscriptFeature> transcriptOption = gene.getTranscripts().stream()
                .filter( x -> x.getTranscriptId().equals(FuncotatorTestConstants.MUC16_TRANSCRIPT) )
                .findFirst();

        return transcriptOption.orElseThrow( () -> new GATKException("Could not get the MUC16 transcript from the MUC16 gene!  The test has mutated!  No good!") );
    }

    private static List<Object[]> addReferenceDataToUnitTestData(final List<Object[]> unitTestData,
                                                                 final String referenceFileName,
                                                                 final FeatureReader<GencodeGtfFeature> featureReader,
                                                                 final ReferenceDataSource referenceDataSource,
                                                                 final String transcriptFastaFile,
                                                                 final String transcriptGtfFile) {

        // unitTestData should include:
//        final String expectedGeneName,
//        final int chromosomeNumber,
//        final int start,
//        final int end,
//        final GencodeFuncotation.VariantClassification expectedVariantClassification,
//        final GencodeFuncotation.VariantType expectedVariantType,
//        final String ref,
//        final String alt,
//        final String expectedGenomeChange,
//        final String expectedStrand,
//        final String expectedCDnaChange,
//        final String expectedCodonChange,
//        final String expectedProteinChange,

        final List<Object[]> outList = new ArrayList<>(unitTestData.size());

        for ( final Object[] rawData : unitTestData ) {
            final Object[] dataWithReference = new Object[rawData.length + 5];
            for ( int i = 0; i < rawData.length; ++i ) {
                dataWithReference[i] = rawData[i];
            }
            dataWithReference[dataWithReference.length-5] = referenceFileName;
            dataWithReference[dataWithReference.length-4] = featureReader;
            dataWithReference[dataWithReference.length-3] = referenceDataSource;
            dataWithReference[dataWithReference.length-2] = transcriptFastaFile;
            dataWithReference[dataWithReference.length-1] = transcriptGtfFile;
            outList.add(dataWithReference);
        }

        return outList;
    }

    private Set<String> getValidTranscriptsForGene(final String expectedGeneName) {

        final Set<String> requestedTranscriptIds = new HashSet<>();
        if ( expectedGeneName != null ){
            if ( expectedGeneName.equals("PIK3CA") ) {
                requestedTranscriptIds.add( FuncotatorTestConstants.PIK3CA_TRANSCRIPT );
            }
            else if ( expectedGeneName.equals("MUC16") ) {
                requestedTranscriptIds.add( FuncotatorTestConstants.MUC16_TRANSCRIPT );
            }
        }
        return requestedTranscriptIds;
    }

    private static GencodeFuncotation createFuncotationForTestGencodeFuncotationComparatorUnitTest(
            final String transcriptId,
            final GencodeGTFFieldConstants.FeatureTag apprisLevel,
            final Integer locusLevel,
            final GencodeFuncotation.VariantClassification variantClassification,
            final Integer transcriptLength
    ) {

        final GencodeFuncotationBuilder builder = new GencodeFuncotationBuilder();

        builder.setAnnotationTranscript( transcriptId );
        builder.setApprisRank( apprisLevel );
        builder.setLocusLevel( locusLevel );
        builder.setVariantClassification( variantClassification );
        builder.setTranscriptLength( transcriptLength );

        return builder.build();
    }

    private GencodeGtfFeatureBaseData createGtfBaseDataForTestIs5Prime(final SimpleInterval interval) {
        return new GencodeGtfFeatureBaseData(GencodeGtfCodec.GTF_FILE_TYPE_STRING, 
                1,
                interval.getContig(),
                GencodeGtfFeature.ANNOTATION_SOURCE_ENSEMBL,
                GencodeGtfFeature.FeatureType.GENE,
                interval.getStart(),
                interval.getEnd(),
                Strand.POSITIVE,
                GencodeGtfFeature.GenomicPhase.DOT,
                "TEST-GENE-ID",
                "TEST-TX-ID",
                GencodeGTFFieldConstants.KnownGeneBiotype.PROTEIN_CODING.toString(),
                GencodeGTFFieldConstants.GeneTranscriptStatus.PUTATIVE.toString(),
                "TEST-GENE",
                GencodeGTFFieldConstants.KnownGeneBiotype.PROTEIN_CODING.toString(),
                GencodeGTFFieldConstants.GeneTranscriptStatus.PUTATIVE.toString(),
                "TEST-TX",
                1,
                "",
                GencodeGTFFieldConstants.LocusLevel.AUTOMATICALLY_ANNOTATED.toString(),
                null
        );
    }

    private GencodeGtfExonFeature helpCreateExonFeature(final SimpleInterval interval,
                                                        final GencodeGtfStartCodonFeature startCodonFeature) {

        final GencodeGtfFeatureBaseData baseData = createGtfBaseDataForTestIs5Prime(interval);
        baseData.featureType = GencodeGtfFeature.FeatureType.EXON;
        final GencodeGtfExonFeature exon = (GencodeGtfExonFeature)GencodeGtfExonFeature.create(baseData);

        if ( startCodonFeature != null ) {
            exon.setStartCodon(startCodonFeature);
        }

        return exon;
    }

    private GencodeGtfTranscriptFeature provideForTestIs5PrimeUtrTranscriptHelper(
                                            final SimpleInterval interval,
                                            final Strand strand,
                                            final GencodeGtfStartCodonFeature startCodonFeature) {
        //    All that is needed for the Transcript is:
        //      Strand
        //      Exons
        //          start codon within an exon

        final GencodeGtfFeatureBaseData baseData = createGtfBaseDataForTestIs5Prime(interval);
        baseData.featureType = GencodeGtfFeature.FeatureType.TRANSCRIPT;
        baseData.genomicStrand = strand;

        final GencodeGtfTranscriptFeature transcriptFeature = (GencodeGtfTranscriptFeature) GencodeGtfTranscriptFeature.create(baseData);

        if ( startCodonFeature == null ) {
            // Always create a list of 3 exons:
            final int NUM_EXONS = 3;

            final int step = (interval.getEnd()-interval.getStart())/NUM_EXONS;
            for (int i = interval.getStart() ; i < interval.getEnd() ; i += step ) {
                transcriptFeature.addExon(
                        helpCreateExonFeature(new SimpleInterval(interval.getContig(), i, i+step-1), null)
                );
            }
        }
        else {
            // Create 2 exons - 1 to contain the start codon, the other to be a placeholder:
            final int NUM_EXONS = 2;

            boolean firstExon = true;
            final int step = (interval.getEnd()-interval.getStart())/NUM_EXONS;
            for (int i = interval.getStart() ; i < interval.getEnd() ; i += step ) {
                GencodeGtfStartCodonFeature startCodonLoopFeature = null;
                if (firstExon) {
                    startCodonLoopFeature = startCodonFeature;
                    firstExon = false;
                }

                transcriptFeature.addExon(
                        helpCreateExonFeature(new SimpleInterval(interval.getContig(), i, i+step-1), startCodonLoopFeature)
                );
            }
        }

        return transcriptFeature;
    }

    private GencodeGtfUTRFeature provideForTestIs5PrimeUtrUTRHelper(final SimpleInterval interval) {
        //    All that is needed for the UTR is:
        //      start
        //      end
        final GencodeGtfFeatureBaseData baseData = createGtfBaseDataForTestIs5Prime(interval);
        baseData.featureType = GencodeGtfFeature.FeatureType.UTR;
        return (GencodeGtfUTRFeature) GencodeGtfUTRFeature.create(baseData);
    }

    private GencodeGtfStartCodonFeature provideForTestIs5PrimeUtrStartCodonHelper(final SimpleInterval interval) {
        //    All that is needed for the StartCodon is:
        //      start
        //      end
        final GencodeGtfFeatureBaseData baseData = createGtfBaseDataForTestIs5Prime(interval);
        baseData.featureType = GencodeGtfFeature.FeatureType.START_CODON;
        return (GencodeGtfStartCodonFeature)GencodeGtfStartCodonFeature.create(baseData);
    }

    private GencodeGtfTranscriptFeature provideArtificialTranscriptForTestMANEExtractMode(final SimpleInterval interval,
                                                                                          final boolean isBasic,
                                                                                          final boolean isMANESelect,
                                                                                          final boolean isMANEPlusClinical ) {
        List<GencodeGtfFeature.OptionalField<?>> optionalFields = new ArrayList<>();
        if ( isBasic ) {
            optionalFields.add(new GencodeGtfFeature.OptionalField<String>("tag", GencodeGTFFieldConstants.FeatureTag.BASIC.toString()));
        }
        if ( isMANEPlusClinical ) {
            optionalFields.add(new GencodeGtfFeature.OptionalField<String>("tag", GencodeGTFFieldConstants.FeatureTag.MANE_PLUS_CLINICAL.toString()));
        }
        if ( isMANESelect ) {
            optionalFields.add(new GencodeGtfFeature.OptionalField<String>("tag", GencodeGTFFieldConstants.FeatureTag.MANE_SELECT.toString()));
        }

        final GencodeGtfTranscriptFeature baseData = (GencodeGtfTranscriptFeature) GencodeGtfTranscriptFeature.create( new GencodeGtfFeatureBaseData(GencodeGtfCodec.GTF_FILE_TYPE_STRING,
                1,
                interval.getContig(),
                GencodeGtfFeature.ANNOTATION_SOURCE_HAVANA,
                GencodeGtfFeature.FeatureType.GENE,
                interval.getStart(),
                interval.getEnd(),
                Strand.POSITIVE,
                GencodeGtfFeature.GenomicPhase.DOT,
                "TEST-GENE-ID",
                "TEST-TX-ID-"+Math.random(), // Randomize the transcript ID to ensure uniqueness for later comparisons
                GencodeGTFFieldConstants.KnownGeneBiotype.PROTEIN_CODING.toString(),
                GencodeGTFFieldConstants.GeneTranscriptStatus.PUTATIVE.toString(),
                "TEST-GENE",
                GencodeGTFFieldConstants.KnownGeneBiotype.PROTEIN_CODING.toString(),
                GencodeGTFFieldConstants.GeneTranscriptStatus.PUTATIVE.toString(),
                "TEST-TX",
                1,
                "",
                GencodeGTFFieldConstants.LocusLevel.AUTOMATICALLY_ANNOTATED.toString(),
                optionalFields.isEmpty() ? null : optionalFields
        ));
        return baseData;
    }


    //==================================================================================================================
    // Data Providers:

    @DataProvider
    Object[][] provideTestForMANESelectMode() {
        final GencodeGtfTranscriptFeature NONBASIC_A = provideArtificialTranscriptForTestMANEExtractMode(new SimpleInterval("chr1", 99, 101), false, false, false);
        final GencodeGtfTranscriptFeature BASIC_A = provideArtificialTranscriptForTestMANEExtractMode(new SimpleInterval("chr1", 99, 101), true, false, false);
        final GencodeGtfTranscriptFeature BASIC_B = provideArtificialTranscriptForTestMANEExtractMode(new SimpleInterval("chr1", 99, 102), true, false, false);
        final GencodeGtfTranscriptFeature BASIC_MANESELECT_A = provideArtificialTranscriptForTestMANEExtractMode(new SimpleInterval("chr1", 99, 101), true, true, false);
        final GencodeGtfTranscriptFeature BASIC_MANESELECT_B = provideArtificialTranscriptForTestMANEExtractMode(new SimpleInterval("chr1", 99, 102), true, true, false);
        final GencodeGtfTranscriptFeature BASIC_MANEPLUSCLINICAL_A = provideArtificialTranscriptForTestMANEExtractMode(new SimpleInterval("chr1", 99, 101), true, false, true);
        final GencodeGtfTranscriptFeature BASIC_MANEPLUSCLINICAL_B = provideArtificialTranscriptForTestMANEExtractMode(new SimpleInterval("chr1", 99, 102), true, false, true);
        final GencodeGtfTranscriptFeature BASIC_MANESELECT_MANEPLUSCLINICAL = provideArtificialTranscriptForTestMANEExtractMode(new SimpleInterval("chr1", 99, 101), true, true, true);
        final GencodeGtfTranscriptFeature NONBASIC_MANESELECT = provideArtificialTranscriptForTestMANEExtractMode(new SimpleInterval("chr1", 99, 101), false, true, false);
        final GencodeGtfTranscriptFeature NONBASIC_MANEPLUSCLINICAL = provideArtificialTranscriptForTestMANEExtractMode(new SimpleInterval("chr1", 99, 101), false, false, true);

        return new Object[][]{
                // trivial cases where only one transcript is present
                {List.of(NONBASIC_A), List.of()},
                {List.of(BASIC_A), List.of(BASIC_A)},
                {List.of(BASIC_MANESELECT_A), List.of(BASIC_MANESELECT_A)},
                {List.of(BASIC_MANEPLUSCLINICAL_A), List.of(BASIC_MANEPLUSCLINICAL_A)},
                {List.of(BASIC_MANESELECT_MANEPLUSCLINICAL), List.of(BASIC_MANESELECT_MANEPLUSCLINICAL)},
                // NOTE that these are (NOT BASIC TRANSCRIPTS), we still return them in MANE select mode
                {List.of(NONBASIC_MANESELECT), List.of(NONBASIC_MANESELECT)},
                {List.of(NONBASIC_MANESELECT), List.of(NONBASIC_MANESELECT)},

                // More complicated cases with multiple transcripts
                {List.of(NONBASIC_A, BASIC_B), List.of(BASIC_B)}, // filter out non-basic transcripts
                {List.of(BASIC_A, BASIC_B), List.of(BASIC_A, BASIC_B)}, // return all basic transcripts
                {List.of(BASIC_A, BASIC_MANESELECT_A), List.of(BASIC_MANESELECT_A)}, // return only MANE select transcript
                {List.of(BASIC_A, BASIC_MANESELECT_A, BASIC_MANESELECT_B), List.of(BASIC_MANESELECT_A, BASIC_MANESELECT_B)}, // return ALL MANE_SELECT transcripts present
                {List.of(BASIC_A, BASIC_MANESELECT_A, NONBASIC_MANESELECT), List.of(BASIC_MANESELECT_A, NONBASIC_MANESELECT)}, // return ALL MANE_SELECT transcripts present, even if one is non-basic
                {List.of(BASIC_A, BASIC_MANEPLUSCLINICAL_B), List.of(BASIC_MANEPLUSCLINICAL_B)}, // return only MANE+Clinical transcript
                {List.of(BASIC_A, BASIC_MANEPLUSCLINICAL_A, BASIC_MANEPLUSCLINICAL_B), List.of(BASIC_MANEPLUSCLINICAL_A, BASIC_MANEPLUSCLINICAL_B)}, // return multiple MANE+Clinical transcripts
                {List.of(BASIC_A, BASIC_MANESELECT_A, BASIC_MANEPLUSCLINICAL_B), List.of(BASIC_MANEPLUSCLINICAL_B)}, // MANE+Clinical Overrides MANE select

                //edge cases
                {List.of(BASIC_A, BASIC_MANESELECT_A, NONBASIC_MANEPLUSCLINICAL), List.of(NONBASIC_MANEPLUSCLINICAL)}, // MANE+Clinical Overrides MANE select EVEN if it is non-basic and MANE select is basic (probably doesn't happen in practice)
                {List.of(BASIC_A, BASIC_MANESELECT_A, BASIC_MANESELECT_B, BASIC_MANESELECT_MANEPLUSCLINICAL), List.of(BASIC_MANESELECT_MANEPLUSCLINICAL)}, // MANE+Clinical trumps MANE_Select if both are present (should not happen in practice)
        };
    }

    @DataProvider
    Object[][] provideTranscriptForGetSortedCdsAndStartStopPositions() {
        return new Object[][] {
                {
                        DataProviderForExampleGencodeGtfGene.createGencodeGtfGeneFeature().getTranscripts().get(0),
                        Arrays.asList(
                                new SimpleInterval("chr1", 99, 101),
                                new SimpleInterval("chr1", 201, 400),
                                new SimpleInterval("chr1", 401, 600),
                                new SimpleInterval("chr1", 601, 800),
                                new SimpleInterval("chr1", 900, 902)
                        )
                },
                {
                        DataProviderForExampleGencodeGtfGene.createGencodeGtfGeneFeature().getTranscripts().get(1),
                        Arrays.asList(
                                new SimpleInterval("chr1", 1099, 1200),
                                new SimpleInterval("chr1", 1201, 1400),
                                new SimpleInterval("chr1", 1401, 1600),
                                new SimpleInterval("chr1", 1601, 1800),
                                new SimpleInterval("chr1", 1801, 1899),
                                new SimpleInterval("chr1", 1900, 1902)
                        )
                },
                {
                        DataProviderForExampleGencodeGtfGene.createGencodeGtfGeneFeature().getTranscripts().get(2),
                        Collections.emptyList()
                },
        };
    }

    @DataProvider
    Object[][] provideMuc16SnpDataForGetVariantClassification() {
        final List<Object[]> l = new ArrayList<>();

        l.addAll( DataProviderForSnpsOnMuc16.provideSnpDataForGetVariantClassification_0() );
        l.addAll( DataProviderForSnpsOnMuc16.provideSnpDataForGetVariantClassification_1() );
        l.addAll( DataProviderForSnpsOnMuc16.provideSnpDataForGetVariantClassification_2() );
        l.addAll( DataProviderForSnpsOnMuc16.provideSnpDataForGetVariantClassification_3() );
        
        return l.toArray(new Object[][]{{}});
    }

    @DataProvider
    Object[][] provideForTestSpliceSiteWindowSettings() {
        return new Object[][] {
                { 19, 8987107, 8987107, GencodeFuncotation.VariantType.SNP, "A", "G", 0, GencodeFuncotation.VariantClassification.SILENT },
                { 19, 8987107, 8987107, GencodeFuncotation.VariantType.SNP, "A", "G", 1, GencodeFuncotation.VariantClassification.SILENT },
                { 19, 8987107, 8987107, GencodeFuncotation.VariantType.SNP, "A", "G", 2, GencodeFuncotation.VariantClassification.SILENT },
                { 19, 8987107, 8987107, GencodeFuncotation.VariantType.SNP, "A", "G", 3, GencodeFuncotation.VariantClassification.SILENT },
                { 19, 8987107, 8987107, GencodeFuncotation.VariantType.SNP, "A", "G", 4, GencodeFuncotation.VariantClassification.SILENT },
                { 19, 8987107, 8987107, GencodeFuncotation.VariantType.SNP, "A", "G", 5, GencodeFuncotation.VariantClassification.SILENT },
                { 19, 8987107, 8987107, GencodeFuncotation.VariantType.SNP, "A", "G", 6, GencodeFuncotation.VariantClassification.SPLICE_SITE },
                { 19, 8987107, 8987107, GencodeFuncotation.VariantType.SNP, "A", "G", 7, GencodeFuncotation.VariantClassification.SPLICE_SITE },
                { 19, 8987107, 8987107, GencodeFuncotation.VariantType.SNP, "A", "G", 8, GencodeFuncotation.VariantClassification.SPLICE_SITE },
                { 19, 8987049, 8987049, GencodeFuncotation.VariantType.SNP, "G", "A", 0, GencodeFuncotation.VariantClassification.SILENT },
                { 19, 8987049, 8987049, GencodeFuncotation.VariantType.SNP, "G", "A", 1, GencodeFuncotation.VariantClassification.SILENT },
                { 19, 8987049, 8987049, GencodeFuncotation.VariantType.SNP, "G", "A", 2, GencodeFuncotation.VariantClassification.SILENT },
                { 19, 8987049, 8987049, GencodeFuncotation.VariantType.SNP, "G", "A", 3, GencodeFuncotation.VariantClassification.SILENT },
                { 19, 8987049, 8987049, GencodeFuncotation.VariantType.SNP, "G", "A", 4, GencodeFuncotation.VariantClassification.SILENT },
                { 19, 8987049, 8987049, GencodeFuncotation.VariantType.SNP, "G", "A", 5, GencodeFuncotation.VariantClassification.SPLICE_SITE },
                { 19, 8987049, 8987049, GencodeFuncotation.VariantType.SNP, "G", "A", 6, GencodeFuncotation.VariantClassification.SPLICE_SITE },
                { 19, 8987049, 8987049, GencodeFuncotation.VariantType.SNP, "G", "A", 7, GencodeFuncotation.VariantClassification.SPLICE_SITE },
        };
    }

    @DataProvider
    Object[][] provideMuc16SnpDataForGetVariantClassificationWithOutOfCdsData() {
        final List<Object[]> l = new ArrayList<>();

        l.addAll( DataProviderForSnpsOnMuc16.provideSnpDataForGetVariantClassification_0() );
        l.addAll( DataProviderForSnpsOnMuc16.provideSnpDataForGetVariantClassification_1() );
        l.addAll( DataProviderForSnpsOnMuc16.provideSnpDataForGetVariantClassification_2() );
        l.addAll( DataProviderForSnpsOnMuc16.provideSnpDataForGetVariantClassification_3() );
        l.addAll( DataProviderForSnpsOnMuc16.provideSnpDataForGetVariantClassification_NotInCDS() );

        return l.toArray(new Object[][]{{}});
    }

    @DataProvider
    Object[][] provideForCreateNonBasicFuncotations() {
        return new Object[][] {
                { 9002502, 9002502 },
                { 9002154, 9002154 },
                { 9001832, 9001832 },
                { 9000443, 9000443 },
                { 9000148, 9000148 },
                { 8999393, 8999393 },
                { 8999026, 8999026 },
                { 8998699, 8998699 },
                { 8997413, 8997413 },
                { 8997119, 8997119 },
                { 8996322, 8996322 },
                { 8995955, 8995955 },
                { 8995636, 8995636 },
                { 8994418, 8994418 },
                { 8994143, 8994143 },
                { 8993374, 8993374 },
                { 8993008, 8993008 },
                { 8987211, 8987211 },
                { 8987046, 8987046 },
                { 8982158, 8982158 },
                { 8973586, 8973586 },
                { 9006643, 9006643 },
                { 9006345, 9006345 },
                { 9005560, 9005560 },
                { 9005195, 9005195 },
                { 9004870, 9004870 },
                { 9003567, 9003567 },
                { 9003294, 9003294 },
                { 9002502, 9002502 },
                { 9002154, 9002154 },
                { 9001832, 9001832 },
                { 9000443, 9000443 },
                { 9000148, 9000148 },
                { 8999393, 8999393 },
                { 8999026, 8999026 },
                { 8998699, 8998699 },
                { 8997413, 8997413 },
                { 8997119, 8997119 },
                { 8996322, 8996322 },
                { 8995955, 8995955 },
                { 8994418, 8994418 },
                { 8994143, 8994143 },
                { 8993374, 8993374 },
                { 8993008, 8993008 },
                { 8987211, 8987211 },
                { 8987046, 8987046 },
                { 8982158, 8982158 },
                { 8979218, 8979218 },
                { 8977641, 8977641 },
                { 8976740, 8976740 },
                { 8976582, 8976582 },
                { 8976261, 8976261 },
                { 8973973, 8973973 },
                { 8973550, 8973550 },
                { 8971677, 8971677 },
                { 8969263, 8969263 },
                { 8968881, 8968881 },
                { 8966750, 8966750 },
                { 9006643, 9006643 },
                { 9006345, 9006345 },
                { 9005560, 9005560 },
                { 9005195, 9005195 },
                { 9004870, 9004870 },
                { 9003567, 9003567 },
                { 9003294, 9003294 },
                { 9002502, 9002502 },
                { 9002154, 9002154 },
                { 9001832, 9001832 },
                { 9000443, 9000443 },
                { 9000148, 9000148 },
                { 8999393, 8999393 },
                { 8999026, 8999026 },
                { 8998699, 8998699 },
                { 8997413, 8997413 },
                { 8997119, 8997119 },
                { 8996322, 8996322 },
                { 8995955, 8995955 },
                { 8995636, 8995636 },
                { 8994418, 8994418 },
                { 8994143, 8994143 },
                { 8993374, 8993374 },
                { 8993008, 8993008 },
                { 8987211, 8987211 },
                { 8987046, 8987046 },
                { 8982158, 8982158 },
                { 8979218, 8979218 },
                { 8976740, 8976740 },
                { 8976582, 8976582 },
                { 8976261, 8976261 },
                { 8973973, 8973973 },
                { 8973550, 8973550 },
                { 8971677, 8971677 },
                { 8969277, 8969277 },
                { 8968881, 8968881 },
                { 8966651, 8966651 },
                { 8962355, 8962355 },
                { 8961953, 8961953 },
                { 8959612, 8959612 },

                // Include some test cases outside the bounds of the MUC16 gene:
                { 9092019, 9092019 },
                { 9092050, 9092050 },
                { 8959519, 8959519 },
                { 8959510, 8959510 },
        };
    }

    @DataProvider
    Object[][] provideDataForCreateFuncotations() {

        final List<Object[]> outList = new ArrayList<>();

        // MUC16 SNPs / DNPs:
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForMuc16MnpFullData.provideMnpDataForMuc16_1(),       FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(), gencodeHg19FeatureReader, refDataSourceHg19Ch19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_GTF_PATH_HG19 ) );
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForMuc16MnpFullData.provideMnpDataForMuc16_2(),       FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(), gencodeHg19FeatureReader, refDataSourceHg19Ch19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_GTF_PATH_HG19  ) );
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForMuc16MnpFullData.provideMnpDataForMuc16_3(),       FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(), gencodeHg19FeatureReader, refDataSourceHg19Ch19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_GTF_PATH_HG19  ) );
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForMuc16MnpFullData.provideMnpDataForMuc16_4(),       FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(), gencodeHg19FeatureReader, refDataSourceHg19Ch19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_GTF_PATH_HG19  ) );
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForMuc16MnpFullData.provideMnpDataForMuc16_5(),       FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(), gencodeHg19FeatureReader, refDataSourceHg19Ch19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_GTF_PATH_HG19  ) );
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForMuc16MnpFullData.provideEdgeCasesForMUC16Data_1(), FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(), gencodeHg19FeatureReader, refDataSourceHg19Ch19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_GTF_PATH_HG19  ) );

        // MUC16 INDELs:
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForMuc16IndelData.provideIndelDataForMuc16(), FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(), gencodeHg19FeatureReader, refDataSourceHg19Ch19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_GTF_PATH_HG19 ) );

        // PIK3CA SNPs / DNPs:
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForPik3caTestData.providePik3caMnpData(), FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(), gencodeHg19FeatureReader, refDataSourceHg19Ch3, FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_GTF_PATH_HG19 ) );

        // PIK3CA INDELs:
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForPik3caTestData.providePik3caInDelData(), FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(), gencodeHg19FeatureReader, refDataSourceHg19Ch3, FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_GTF_PATH_HG19 ) );

        // PIK3CA Other Indels:
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForPik3caTestData.providePik3caInDelData2(), FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(), gencodeHg19FeatureReader, refDataSourceHg19Ch3, FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19, FuncotatorTestConstants.GENCODE_DATA_SOURCE_GTF_PATH_HG19 ) );

        // Trouble variants:
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForTroubleVariants.provideSymbolicAllelesAndMaskedBasesForHg38(), hg38Reference, gencodeHg38FeatureReader, refDataSourceHg38, FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG38, FuncotatorTestConstants.GENCODE_DATA_SOURCE_GTF_PATH_HG38 ) );
        outList.addAll( addReferenceDataToUnitTestData(DataProviderForTroubleVariants.provideTP53AllelesForOrderTestsHg19(), FuncotatorReferenceTestUtils.retrieveHg19Chr17Ref(), gencodeHg19TP53FeatureReader, refDataSourceHg19Ch17, FuncotatorTestConstants.TP53_REVERSE_ORDER_GENCODE_TRANSCRIPT_FASTA_FILE, FuncotatorTestConstants.TP53_REVERSE_ORDER_GENCODE_TRANSCRIPT_FASTA_FILE ) );

        return outList.toArray(new Object[][]{{}});
    }

    @DataProvider
    Object[][] provideForTestGetTranscriptEndPaddingBases() {

                          // Base Position (read from bottom to top):
                          //
                          //1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
                          //0000000001111111111222222222233333333334444444444555555555566666666667777777777888888888899999999990
                          //0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001
        final String seq = "AAAAATTTTTGGGGGCCCCCAAAAATTTTTGGGGGCCCCCAAAAATTTTTGGGGGCCCCCAAAAATTTTTGGGGGCCCCCAAAAATTTTTGGGGGCCCCC";
        final String contig = "test";

        final List<Locatable> exonPositionList = new ArrayList<>();
        exonPositionList.add( new SimpleInterval(contig, 5,15) );
        exonPositionList.add( new SimpleInterval(contig, 17,25) );
        exonPositionList.add( new SimpleInterval(contig, 35,55) );
        exonPositionList.add( new SimpleInterval(contig, 56,90) );

//        final VariantContext variant,
//        final Allele altAllele,
//        final List<? extends htsjdk.samtools.util.Locatable> exonPositionList,
//        final ReferenceContext reference,
//        final String expected

        return new Object[][] {
                // Not Indel:
                helpProvideForTestGetTranscriptEndPaddingBases("A", "G", 5, exonPositionList, seq, ""),
                helpProvideForTestGetTranscriptEndPaddingBases(seq.substring(seq.length()-1, seq.length()), "A", seq.length(), exonPositionList, seq, ""),
                // Insertion 1 base, out of range:
                helpProvideForTestGetTranscriptEndPaddingBases("A", "AG", 5, exonPositionList, seq, ""),
                // Insertion 1 base, border of range:
                helpProvideForTestGetTranscriptEndPaddingBases("C", "CG", 80, exonPositionList, seq, ""),
                // Insertion 1 base, in range:
                helpProvideForTestGetTranscriptEndPaddingBases("C", "CG", 82, exonPositionList, seq, "GGGGGC"),
                // Insertion 2 bases, out of range:
                helpProvideForTestGetTranscriptEndPaddingBases("A", "AGC", 5, exonPositionList, seq, ""),
                // Insertion 2 bases, border of range:
                helpProvideForTestGetTranscriptEndPaddingBases("C", "CGC", 80, exonPositionList, seq, ""),
                // Insertion 2 bases, in range:
                helpProvideForTestGetTranscriptEndPaddingBases("C", "CGA", 82, exonPositionList, seq, "GGGGGC"),
                // Insertion 3 bases, out of range:
                helpProvideForTestGetTranscriptEndPaddingBases("A", "AGCG", 5, exonPositionList, seq, ""),
                // Insertion 3 bases, border of range:
                helpProvideForTestGetTranscriptEndPaddingBases("C", "CGCT", 80, exonPositionList, seq, ""),
                // Insertion 3 bases, in range:
                helpProvideForTestGetTranscriptEndPaddingBases("C", "CGAT", 82, exonPositionList, seq, "GGGGGC"),
                // Insertion 4 bases, out of range:
                helpProvideForTestGetTranscriptEndPaddingBases("A", "AGCGT", 5, exonPositionList, seq, ""),
                // Insertion 4 bases, border of range:
                helpProvideForTestGetTranscriptEndPaddingBases("C", "CGCTG", 80, exonPositionList, seq, ""),
                // Insertion 4 bases, in range:
                helpProvideForTestGetTranscriptEndPaddingBases("C", "CGATT", 82, exonPositionList, seq, "GGGGGCCCC"),
                // Deletion 1 base, out of range:
                helpProvideForTestGetTranscriptEndPaddingBases("AT", "A", 5, exonPositionList, seq, ""),
                // Deletion 1 base, border of range:
                helpProvideForTestGetTranscriptEndPaddingBases("CA", "C", 80, exonPositionList, seq, ""),
                // Deletion 1 base, in range:
                helpProvideForTestGetTranscriptEndPaddingBases("AA", "A", 82, exonPositionList, seq, "GGGGGC"),
                // Deletion 2 bases, out of range:
                helpProvideForTestGetTranscriptEndPaddingBases("ATT", "A", 5, exonPositionList, seq, ""),
                // Deletion 2 bases, border of range:
                helpProvideForTestGetTranscriptEndPaddingBases("CAA", "C", 80, exonPositionList, seq, ""),
                // Deletion 2 bases, in range:
                helpProvideForTestGetTranscriptEndPaddingBases("AAA", "A", 82, exonPositionList, seq, "GGGGGC"),
                // Deletion 3 bases, out of range:
                helpProvideForTestGetTranscriptEndPaddingBases("ATTT", "A", 5, exonPositionList, seq, ""),
                // Deletion 3 bases, border of range:
                helpProvideForTestGetTranscriptEndPaddingBases("CAAA", "C", 80, exonPositionList, seq, ""),
                // Deletion 3 bases, in range:
                helpProvideForTestGetTranscriptEndPaddingBases("AAAA", "A", 82, exonPositionList, seq, "GGGGGC"),
                // Deletion 3 bases, out of range:
                helpProvideForTestGetTranscriptEndPaddingBases("ATTTT", "A", 5, exonPositionList, seq, ""),
                // Deletion 3 bases, border of range:
                helpProvideForTestGetTranscriptEndPaddingBases("CAAAA", "C", 80, exonPositionList, seq, ""),
                // Deletion 3 bases, in range:
                helpProvideForTestGetTranscriptEndPaddingBases("AAAAT", "A", 82, exonPositionList, seq, "GGGGGCCCC"),
        };
    }

    @DataProvider
    Object[][] provideForTestCreateGencodeFuncotationBuilderWithTrivialFieldsPopulated() {

        //Line 1229 or so

        // Just simple field checks are OK.
        // Doesn't need to be fancy.

//        final VariantContext variant,
//        final Allele altAllele,
//        final GencodeGtfGeneFeature gtfFeature,
//        final GencodeGtfTranscriptFeature transcript

        final SimpleInterval variantInterval =  new SimpleInterval("chr3", 178921515, 178921517);
        final Allele refAllele = Allele.create("GCA", true);
        final Allele altAllele = Allele.create("TTG");
        final VariantContext variant = new VariantContextBuilder(
                    FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
                    variantInterval.getContig(),
                    variantInterval.getStart(),
                    variantInterval.getEnd(),
                    Arrays.asList(refAllele, altAllele)
                ).make();


        final GencodeGtfGeneFeature gene = DataProviderForExampleGencodeGtfGene.createGencodeGtfGeneFeature();

        return new Object[][] {
                { variant, altAllele, gene, gene.getTranscripts().get(0) }
        };
    }

    @DataProvider
    Object[][] provideForCreateDefaultFuncotationsOnProblemVariant() {

//        variant, altAllele, gtfFeature, reference, transcript, version

        final SimpleInterval variantInterval =  new SimpleInterval("chr3", 178921515, 178921517);
        final Allele refAllele = Allele.create("GCA", true);
        final Allele altAllele = Allele.create("TTG");
        final VariantContext variant = new VariantContextBuilder(
                FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
                variantInterval.getContig(),
                variantInterval.getStart(),
                variantInterval.getEnd(),
                Arrays.asList(refAllele, altAllele)
        ).make();

        final String versionString = "VERSION";
        final String dataSourceName = "TEST_GENCODE_NAME";

        // ======================
        // Create the GencodeGtfFeature:
        GencodeGtfFeatureBaseData data;

        data = new GencodeGtfFeatureBaseData(GencodeGtfCodec.GTF_FILE_TYPE_STRING, 1, variantInterval.getContig(), GencodeGtfFeature.ANNOTATION_SOURCE_ENSEMBL, GencodeGtfFeature.FeatureType.GENE,
                variantInterval.getStart()-2000, variantInterval.getEnd()+2000, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", null, GencodeGTFFieldConstants.KnownGeneBiotype.PROTEIN_CODING.toString(),
                null, "TEST_GENE", null, null, null, -1, null, GencodeGTFFieldConstants.LocusLevel.AUTOMATICALLY_ANNOTATED.toString(), null);
        final GencodeGtfGeneFeature gene = (GencodeGtfGeneFeature)GencodeGtfFeature.create(data);

        // ======================

        data = new GencodeGtfFeatureBaseData(GencodeGtfCodec.GTF_FILE_TYPE_STRING, 2, variantInterval.getContig(), GencodeGtfFeature.ANNOTATION_SOURCE_ENSEMBL, GencodeGtfFeature.FeatureType.TRANSCRIPT,
                variantInterval.getStart()-1000, variantInterval.getEnd()+1000, Strand.POSITIVE, GencodeGtfFeature.GenomicPhase.DOT, "TEST_GENE1", "TEST_TRANSCRIPT1", GencodeGTFFieldConstants.KnownGeneBiotype.PROTEIN_CODING.toString(),
                null, "TEST_GENE", GencodeGTFFieldConstants.KnownGeneBiotype.PROTEIN_CODING.toString(), null, "TEST_TRANSCRIPT1", -1, null, GencodeGTFFieldConstants.LocusLevel.AUTOMATICALLY_ANNOTATED.toString(),
                Collections.emptyList()
        );
        final GencodeGtfTranscriptFeature transcript1 = (GencodeGtfTranscriptFeature) GencodeGtfFeature.create(data);
        gene.addTranscript(transcript1);

        // ======================

        return new Object[][] {
                {
                    variant,
                    variant.getAlleles().get(1),
                    gene,
                    new ReferenceContext( refDataSourceHg19Ch3, variantInterval ),
                    gene.getTranscripts().get(0),
                    versionString,
                    dataSourceName,
                    new GencodeFuncotationBuilder()
                            .setDataSourceName(dataSourceName)
                            .setHugoSymbol(gene.getGeneName())
                            .setChromosome(variant.getContig())
                            .setStart(variant.getStart())
                            .setEnd(variant.getEnd())
                            .setVariantClassification(GencodeFuncotation.VariantClassification.COULD_NOT_DETERMINE)
                            .setVariantType(GencodeFuncotation.VariantType.TNP)
                            .setRefAllele(variant.getReference())
                            .setTumorSeqAllele2(variant.getAlternateAllele(0).getBaseString())
                            .setGenomeChange("g.chr3:178921515_178921517GCA>TTG")
                            .setAnnotationTranscript(transcript1.getTranscriptName())
                            .setStrand(gene.getGenomicStrand())
                            .setReferenceContext("TATAAATAGTGCACTCAGAATAA")
                            .setGcContent(0.3399503722084367)
                            .setLocusLevel(Integer.valueOf(gene.getLocusLevel().toString()))
                            // This is OK because there are no exons:
                            .setTranscriptLength(0)
                            .setVersion(versionString)
                            .setGeneTranscriptType(transcript1.getTranscriptType())
                            .setNcbiBuild("TEST")
                            .build()
                },
        };
    }

    @DataProvider
    Object[][] provideForTestIs5PrimeUtr() {

        // Need:
        //    Transcript:
        //      Strand
        //      Exons
        //          start codon within an exon
        //
        //    StartCodon:
        //      start
        //      end
        //
        //    UTR:
        //      start
        //      end

        final SimpleInterval transcriptInterval = new SimpleInterval("T", 1, 2000);

        final GencodeGtfStartCodonFeature startCodonFeatureFront = provideForTestIs5PrimeUtrStartCodonHelper(new SimpleInterval("T", 500, 502));
        final GencodeGtfStartCodonFeature startCodonFeatureBack = provideForTestIs5PrimeUtrStartCodonHelper(new SimpleInterval("T", 1498, 1500));

        final GencodeGtfUTRFeature utr = provideForTestIs5PrimeUtrUTRHelper(new SimpleInterval("T", 750, 1250));

        return new Object[][] {
                // Strand +
                    // Start Codon is Null
                { utr, provideForTestIs5PrimeUtrTranscriptHelper(transcriptInterval, Strand.POSITIVE, null), false },
                    // Start Codon non-null
                { utr, provideForTestIs5PrimeUtrTranscriptHelper(transcriptInterval, Strand.POSITIVE, startCodonFeatureFront), false },
                { utr, provideForTestIs5PrimeUtrTranscriptHelper(transcriptInterval, Strand.POSITIVE, startCodonFeatureBack), true },
                // Strand -
                    // Start Codon is Null
                { utr, provideForTestIs5PrimeUtrTranscriptHelper(transcriptInterval, Strand.NEGATIVE, null), false },
                    // Start Codon non-null
                { utr, provideForTestIs5PrimeUtrTranscriptHelper(transcriptInterval, Strand.NEGATIVE, startCodonFeatureFront), true },
                { utr, provideForTestIs5PrimeUtrTranscriptHelper(transcriptInterval, Strand.NEGATIVE, startCodonFeatureBack), false },
        };
    }

    @DataProvider
    Object[][] provideDataForTestIsVariantInCodingRegion() {
        return new Object[][] {

                // Trivially false cases:
                {GencodeFuncotation.VariantClassification.INTRON, null, false},
                {GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantClassification.INTRON, false},
                {GencodeFuncotation.VariantClassification.INTRON, GencodeFuncotation.VariantClassification.MISSENSE, false},
                {GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR, null, false},
                {GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR, GencodeFuncotation.VariantClassification.INTRON, false},
                {GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR, GencodeFuncotation.VariantClassification.MISSENSE, false},
                {GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, null, false},
                {GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, GencodeFuncotation.VariantClassification.INTRON, false},
                {GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, GencodeFuncotation.VariantClassification.MISSENSE, false},
                {GencodeFuncotation.VariantClassification.IGR, null, false},
                {GencodeFuncotation.VariantClassification.IGR, GencodeFuncotation.VariantClassification.INTRON, false},
                {GencodeFuncotation.VariantClassification.IGR, GencodeFuncotation.VariantClassification.MISSENSE, false},
                {GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK, null, false},
                {GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK, GencodeFuncotation.VariantClassification.INTRON, false},
                {GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK, GencodeFuncotation.VariantClassification.MISSENSE, false},
                {GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME, null, false},
                {GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME, GencodeFuncotation.VariantClassification.INTRON, false},
                {GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME, GencodeFuncotation.VariantClassification.MISSENSE, false},
                {GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME, null, false},
                {GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME, GencodeFuncotation.VariantClassification.INTRON, false},
                {GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME, GencodeFuncotation.VariantClassification.MISSENSE, false},
                {GencodeFuncotation.VariantClassification.RNA, null, false},
                {GencodeFuncotation.VariantClassification.RNA, GencodeFuncotation.VariantClassification.INTRON, false},
                {GencodeFuncotation.VariantClassification.RNA, GencodeFuncotation.VariantClassification.MISSENSE, false},
                {GencodeFuncotation.VariantClassification.LINCRNA, null, false},
                {GencodeFuncotation.VariantClassification.LINCRNA, GencodeFuncotation.VariantClassification.INTRON, false},
                {GencodeFuncotation.VariantClassification.LINCRNA, GencodeFuncotation.VariantClassification.MISSENSE, false},

                // Trivially true cases:
                {GencodeFuncotation.VariantClassification.MISSENSE, null, true},
                {GencodeFuncotation.VariantClassification.MISSENSE, GencodeFuncotation.VariantClassification.INTRON, true},
                {GencodeFuncotation.VariantClassification.MISSENSE, GencodeFuncotation.VariantClassification.MISSENSE, true},
                {GencodeFuncotation.VariantClassification.NONSENSE, null, true},
                {GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.INTRON, true},
                {GencodeFuncotation.VariantClassification.NONSENSE, GencodeFuncotation.VariantClassification.MISSENSE, true},
                {GencodeFuncotation.VariantClassification.NONSTOP, null, true},
                {GencodeFuncotation.VariantClassification.NONSTOP, GencodeFuncotation.VariantClassification.INTRON, true},
                {GencodeFuncotation.VariantClassification.NONSTOP, GencodeFuncotation.VariantClassification.MISSENSE, true},
                {GencodeFuncotation.VariantClassification.SILENT, null, true},
                {GencodeFuncotation.VariantClassification.SILENT, GencodeFuncotation.VariantClassification.INTRON, true},
                {GencodeFuncotation.VariantClassification.SILENT, GencodeFuncotation.VariantClassification.MISSENSE, true},
                {GencodeFuncotation.VariantClassification.IN_FRAME_DEL, null, true},
                {GencodeFuncotation.VariantClassification.IN_FRAME_DEL, GencodeFuncotation.VariantClassification.INTRON, true},
                {GencodeFuncotation.VariantClassification.IN_FRAME_DEL, GencodeFuncotation.VariantClassification.MISSENSE, true},
                {GencodeFuncotation.VariantClassification.IN_FRAME_INS, null, true},
                {GencodeFuncotation.VariantClassification.IN_FRAME_INS, GencodeFuncotation.VariantClassification.INTRON, true},
                {GencodeFuncotation.VariantClassification.IN_FRAME_INS, GencodeFuncotation.VariantClassification.MISSENSE, true},
                {GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS, null, true},
                {GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS, GencodeFuncotation.VariantClassification.INTRON, true},
                {GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS, GencodeFuncotation.VariantClassification.MISSENSE, true},
                {GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL, null, true},
                {GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL, GencodeFuncotation.VariantClassification.INTRON, true},
                {GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL, GencodeFuncotation.VariantClassification.MISSENSE, true},
                {GencodeFuncotation.VariantClassification.START_CODON_SNP, null, true},
                {GencodeFuncotation.VariantClassification.START_CODON_SNP, GencodeFuncotation.VariantClassification.INTRON, true},
                {GencodeFuncotation.VariantClassification.START_CODON_SNP, GencodeFuncotation.VariantClassification.MISSENSE, true},
                {GencodeFuncotation.VariantClassification.START_CODON_INS, null, true},
                {GencodeFuncotation.VariantClassification.START_CODON_INS, GencodeFuncotation.VariantClassification.INTRON, true},
                {GencodeFuncotation.VariantClassification.START_CODON_INS, GencodeFuncotation.VariantClassification.MISSENSE, true},
                {GencodeFuncotation.VariantClassification.START_CODON_DEL, null, true},
                {GencodeFuncotation.VariantClassification.START_CODON_DEL, GencodeFuncotation.VariantClassification.INTRON, true},
                {GencodeFuncotation.VariantClassification.START_CODON_DEL, GencodeFuncotation.VariantClassification.MISSENSE, true},

                // Splice site special cases:
                {GencodeFuncotation.VariantClassification.SPLICE_SITE, GencodeFuncotation.VariantClassification.INTRON, false},
                {GencodeFuncotation.VariantClassification.SPLICE_SITE, GencodeFuncotation.VariantClassification.MISSENSE, true},
        };
    }

    @DataProvider
    Object[][] provideDataForTestGencodeFuncotationComparatorUnitTest() {

        final String transcriptId1 = "ENST0000123456";
        final String transcriptId2 = "ENST0000987654";
        final Set<String> transcriptSet = new HashSet<>( Arrays.asList(transcriptId1, transcriptId2) );

        return new Object[][] {
                // ==================================================================================================
                // CanonicalGencodeFuncotationComparator
                // Transcript list:
                {
                        // Goes down to Locus Level:
                        TranscriptSelectionMode.CANONICAL.getComparator(new HashSet<>()),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        -1
                },
                {
                        // Goes down to Locus Level:
                        TranscriptSelectionMode.CANONICAL.getComparator(new HashSet<>(Collections.singletonList("TEST"))),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest("", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest("", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        1
                },
                // Locus Level:
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL,  null, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, null, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        1
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        1
                },
                // IGR / non-IGR:
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.IGR, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.IGR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        1
                },
                // Variant Classification:
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.SILENT, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.SILENT, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        1
                },
                // Appris Annotation:
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", null, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", null, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        1
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_CANDIDATE_HIGHEST_SCORE, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL.ordinal() - GencodeGTFFieldConstants.FeatureTag.APPRIS_CANDIDATE_HIGHEST_SCORE.ordinal()
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_CANDIDATE_HIGHEST_SCORE, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        GencodeGTFFieldConstants.FeatureTag.APPRIS_CANDIDATE_HIGHEST_SCORE.ordinal() - GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL.ordinal()
                },
                // Transcript Length:
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, null),
                        -1
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, null),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        1
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 500),
                        -1
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 500),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        1
                },
                // ABC Order of Transcript ID:
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(null, GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(null, GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        1
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        transcriptId1.compareTo(transcriptId2)
                },
                {
                        TranscriptSelectionMode.CANONICAL.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        transcriptId2.compareTo(transcriptId1)
                },


                // ==================================================================================================
                // BestEffectGencodeFuncotationComparator
                // Transcript list:
                {
                        // Goes down to Locus Level:
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(new HashSet<>()),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        -1
                },
                {
                        // Goes down to Locus Level:
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(new HashSet<>(Collections.singletonList("TEST"))),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest("NOT_USER_REQUESTED", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest("NOT_USER_REQUESTED", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        1
                },
                // IGR / non-IGR:
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.IGR, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.IGR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        1
                },
                // Variant Classification:
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.SILENT, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.SILENT, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        1
                },
                // Locus Level:
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL,  null, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, null, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        1
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                        1
                },
                // Appris Annotation:
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", null, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", null, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        1
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_CANDIDATE_HIGHEST_SCORE, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL.ordinal() - GencodeGTFFieldConstants.FeatureTag.APPRIS_CANDIDATE_HIGHEST_SCORE.ordinal()
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_CANDIDATE_HIGHEST_SCORE, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        GencodeGTFFieldConstants.FeatureTag.APPRIS_CANDIDATE_HIGHEST_SCORE.ordinal() - GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL.ordinal()
                },
                // Transcript Length:
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, null),
                        -1
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, null),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        1
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 500),
                        -1
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 500),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        1
                },
                // ABC Order of Transcript ID:
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(null, GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        -1
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(null, GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        1
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        transcriptId1.compareTo(transcriptId2)
                },
                {
                        TranscriptSelectionMode.BEST_EFFECT.getComparator(transcriptSet),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1 + ".5", GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.NONSTOP, 5000),
                        transcriptId2.compareTo(transcriptId1)
                }
        };
    }

    @DataProvider
    Object[][] provideDataForTestCalculateGcContent() {

                          // Base Position:
                          //
                          //1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
                          //0000000001111111111222222222233333333334444444444555555555566666666667777777777888888888899999999990
                          //0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001
        final String seq = "AAAAATTTTTGGGGGCCCCCAAAAATTTTTGGGGGCCCCCAAAAATTTTTGGGGGCCCCCAAAAATTTTTGGGGGCCCCCAAAAATTTTTGGGGGCCCCC";
        final String contig = "test";

        return new Object[][] {
                // MNPs:
                { Allele.create("A", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,  1,  1),  1, 0.0 },
                { Allele.create("A", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,  1,  1),  9, 0.0 },
                { Allele.create("A", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,  1,  1), 10, 1.0/11.0 },
                { Allele.create("C", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,100,100),  1, 1 },
                { Allele.create("C", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,100,100),  9, 1 },
                { Allele.create("C", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,100,100), 10, 10.0/11.0 },
                { Allele.create("T", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 50, 50),  1, 1.0/3.0 },
                { Allele.create("T", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 50, 50),  9, 9.0/19.0 },
                { Allele.create("T", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 50, 50), 10, 11.0/21.0 },
                { Allele.create("T", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 50, 50), 50, 50.0/100.0 },
                { Allele.create("AAAAA",  true), Allele.create("GGGGG"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 1,  5),   1, 0.0 },
                { Allele.create("AAAAA",  true), Allele.create("GGGGG"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 1,  5),   9, 4.0 / 14.0 },
                { Allele.create("AAAAA",  true), Allele.create("GGGGG"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 1,  5),  10, 5.0/15.0 },
                { Allele.create("GCCCCC", true), Allele.create("AGGGGG"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,95,100),   1, 1 },
                { Allele.create("GCCCCC", true), Allele.create("AGGGGG"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,95,100),   9, 10.0/15.0 },
                { Allele.create("GCCCCC", true), Allele.create("AGGGGG"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,95,100),  10, 10.0/16.0 },
                { Allele.create("TGGGGG", true), Allele.create("AAAAAA"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,50, 55),   1, 6.0/8.0 },
                { Allele.create("TGGGGG", true), Allele.create("AAAAAA"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,50, 55),   9, 10.0/24.0 },
                { Allele.create("TGGGGG", true), Allele.create("AAAAAA"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,50, 55),  10, 11.0/26.0 },
                { Allele.create("TGGGGG", true), Allele.create("AAAAAA"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,50, 55),  50, 50.0/100.0 },

                // Insertions:
                { Allele.create("A", true), Allele.create("AG"),    FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,  1,  1),  1, 0.0 },
                { Allele.create("A", true), Allele.create("AGG"),   FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,  1,  1),  9, 0.0 },
                { Allele.create("A", true), Allele.create("AGGG"),  FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,  1,  1), 10, 1.0/11.0 },
                { Allele.create("C", true), Allele.create("CG"),    FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,100,100),  1, 1 },
                { Allele.create("C", true), Allele.create("CGG"),   FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,100,100),  9, 1 },
                { Allele.create("C", true), Allele.create("CGGG"),  FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,100,100), 10, 1 },
                { Allele.create("T", true), Allele.create("TG"),    FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 50, 50),  1, 1.0/2.0 },
                { Allele.create("T", true), Allele.create("TGG"),   FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 50, 50),  9, 9.0/18.0 },
                { Allele.create("T", true), Allele.create("TGGG"),  FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 50, 50), 10, 10.0/20.0 },
                { Allele.create("T", true), Allele.create("TGGGG"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 50, 50), 49, 49.0/98.0 },

                // Deletions:
                { Allele.create("AAAAA",  true), Allele.create("A"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig, 1,  5),  10,  5.0/15.0 },
                { Allele.create("GCCCCC", true), Allele.create("G"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,95,100),  10, 10.0/15.0  },
                { Allele.create("TGGGGG", true), Allele.create("T"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,50, 55),  10, 10.0/25.0 },
                { Allele.create("TG", true),     Allele.create("T"), FuncotatorTestUtils.createReferenceContextFromBasesAndLocation(seq, contig,50, 51),  10, 10.0/21.0 },
        };
    }

    //==================================================================================================================
    // Tests:

    @Test ( dataProvider = "provideTestForMANESelectMode" )
    void testMANESelectTranscriptSelectionCriteria(final List<GencodeGtfTranscriptFeature> inputTranscripts, final List<GencodeGtfTranscriptFeature> expectedTranscripts) {
        final List<GencodeGtfTranscriptFeature> selectedTranscripts = GencodeFuncotationFactory.retreiveMANESelectModeTranscriptsCriteria(inputTranscripts);
        Assert.assertEquals(selectedTranscripts.size(), expectedTranscripts.size());
        for (int i = 0; i < selectedTranscripts.size(); i++) {
            Assert.assertEquals(selectedTranscripts.get(i), expectedTranscripts.get(i));
        }
    }

    @Test ( dataProvider = "provideTranscriptForGetSortedCdsAndStartStopPositions")
    void testGetSortedExonAndStartStopPositions(final GencodeGtfTranscriptFeature transcript, final List<? extends Locatable> expected) {

        final List<? extends Locatable> exons = GencodeFuncotationFactory.getSortedCdsAndStartStopPositions(transcript);

        Assert.assertEquals(exons.size(), expected.size());

        for( int i = 0; i < exons.size() ; ++i ) {
            final SimpleInterval interval = new SimpleInterval( exons.get(i).getContig(), exons.get(i).getStart(), exons.get(i).getEnd());
            Assert.assertEquals( interval, expected.get(i) );
        }
    }

    @Test (dataProvider = "provideMuc16SnpDataForGetVariantClassification")
    void testGetVariantClassificationForCodingRegions(final int chromosomeNumber,
                                      final int start,
                                      final int end,
                                      final GencodeFuncotation.VariantType variantType,
                                      final String ref,
                                      final String alt,
                                      final GencodeFuncotation.VariantClassification expectedVariantClassification) {

        // This test can only deal with variants in coding regions.
        // So we ignore any tests that are expected outside of coding regions.
        // i.e. expectedVariantClassification is one of:
        //     { INTRON, FIVE_PRIME_UTR, THREE_PRIME_UTR, IGR, FIVE_PRIME_FLANK, DE_NOVO_START_IN_FRAME, DE_NOVO_START_OUT_FRAME, RNA, LINCRNA }
        // We test these cases in another unit test.
        if ((expectedVariantClassification == GencodeFuncotation.VariantClassification.INTRON) ||
            (expectedVariantClassification == GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR) ||
            (expectedVariantClassification == GencodeFuncotation.VariantClassification.THREE_PRIME_UTR) ||
            (expectedVariantClassification == GencodeFuncotation.VariantClassification.IGR) ||
            (expectedVariantClassification == GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK) ||
            (expectedVariantClassification == GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME) ||
            (expectedVariantClassification == GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME) ||
            (expectedVariantClassification == GencodeFuncotation.VariantClassification.RNA) ||
            (expectedVariantClassification == GencodeFuncotation.VariantClassification.LINCRNA) )
        {
            return;
        }

        final String contig = "chr" + Integer.toString(chromosomeNumber);
        final SimpleInterval variantInterval = new SimpleInterval( contig, start, end );

        final Allele refAllele = Allele.create(ref, true);
        final Allele altAllele = Allele.create(alt);

        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder(
                FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(),
                contig,
                start,
                end,
                Arrays.asList(refAllele, altAllele)
        );
        final VariantContext variantContext = variantContextBuilder.make();

        // Get our gene feature iterator:
        final CloseableTribbleIterator<GencodeGtfFeature> gtfFeatureIterator;
        try {
            gtfFeatureIterator = gencodeHg19FeatureReader.query(contig, start, end);
        }
        catch (final IOException ex) {
            throw new GATKException("Could not finish the test!", ex);
        }

        // Get the gene.
        // We know the first gene is the right one - the gene in question is the MUC16 gene:
        final GencodeGtfGeneFeature             gene = (GencodeGtfGeneFeature) gtfFeatureIterator.next();
        final GencodeGtfTranscriptFeature transcript = getMuc16Transcript(gene);
        final GencodeGtfExonFeature             exon = getExonForVariant( gene, variantInterval );

        final ReferenceContext referenceContext = new ReferenceContext( refDataSourceHg19Ch19, variantInterval );

        final List<? extends Locatable> exonPositionList = GencodeFuncotationFactory.getSortedCdsAndStartStopPositions(transcript);

        final ReferenceDataSource muc16TranscriptDataSource = ReferenceDataSource.of(new File(FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19).toPath());
        final Map<String, GencodeFuncotationFactory.MappedTranscriptIdInfo> muc16TranscriptIdMap = GencodeFuncotationFactory. createTranscriptIdMap(muc16TranscriptDataSource);

        final SequenceComparison seqComp =
                GencodeFuncotationFactory.createSequenceComparison(
                        variantContext,
                        altAllele,
                        referenceContext,
                        transcript,
                        exonPositionList,
                        muc16TranscriptIdMap,
                        false);

        final GencodeFuncotation.VariantClassification varClass = GencodeFuncotationFactory.createVariantClassification(
                variantContext,
                altAllele,
                variantType,
                exon,
                transcript.getExons().size(),
                seqComp,
                FuncotatorUtils.DEFAULT_SPLICE_SITE_WINDOW_SIZE
        );

        Assert.assertEquals( varClass, expectedVariantClassification );
    }

    @Test (dataProvider = "provideForTestSpliceSiteWindowSettings")
    void testSpliceSiteWindowSettings(final int chromosomeNumber,
                                      final int start,
                                      final int end,
                                      final GencodeFuncotation.VariantType variantType,
                                      final String ref,
                                      final String alt,
                                      final int spliceSiteWindow,
                                      final GencodeFuncotation.VariantClassification expectedVariantClassification) {

        // This test can only deal with variants in coding regions.
        // So we ignore any tests that are expected outside of coding regions.
        // i.e. expectedVariantClassification is one of:
        //     { INTRON, FIVE_PRIME_UTR, THREE_PRIME_UTR, IGR, FIVE_PRIME_FLANK, DE_NOVO_START_IN_FRAME, DE_NOVO_START_OUT_FRAME, RNA, LINCRNA }
        // We test these cases in another unit test.
        if ((expectedVariantClassification == GencodeFuncotation.VariantClassification.INTRON) ||
                (expectedVariantClassification == GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR) ||
                (expectedVariantClassification == GencodeFuncotation.VariantClassification.THREE_PRIME_UTR) ||
                (expectedVariantClassification == GencodeFuncotation.VariantClassification.IGR) ||
                (expectedVariantClassification == GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK) ||
                (expectedVariantClassification == GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME) ||
                (expectedVariantClassification == GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME) ||
                (expectedVariantClassification == GencodeFuncotation.VariantClassification.RNA) ||
                (expectedVariantClassification == GencodeFuncotation.VariantClassification.LINCRNA) )
        {
            return;
        }

        final String contig = "chr" + Integer.toString(chromosomeNumber);
        final SimpleInterval variantInterval = new SimpleInterval( contig, start, end );

        final Allele refAllele = Allele.create(ref, true);
        final Allele altAllele = Allele.create(alt);

        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder(
                FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(),
                contig,
                start,
                end,
                Arrays.asList(refAllele, altAllele)
        );
        final VariantContext variantContext = variantContextBuilder.make();

        // Get our gene feature iterator:
        final CloseableTribbleIterator<GencodeGtfFeature> gtfFeatureIterator;
        try {
            gtfFeatureIterator = gencodeHg19FeatureReader.query(contig, start, end);
        }
        catch (final IOException ex) {
            throw new GATKException("Could not finish the test!", ex);
        }

        // Get the gene.
        // We know the first gene is the right one - the gene in question is the MUC16 gene:
        final GencodeGtfGeneFeature             gene = (GencodeGtfGeneFeature) gtfFeatureIterator.next();
        final GencodeGtfTranscriptFeature transcript = getMuc16Transcript(gene);
        final GencodeGtfExonFeature             exon = getExonForVariant( gene, variantInterval );

        final ReferenceContext referenceContext = new ReferenceContext( refDataSourceHg19Ch19, variantInterval );

        final List<? extends Locatable> exonPositionList = GencodeFuncotationFactory.getSortedCdsAndStartStopPositions(transcript);

        final ReferenceDataSource muc16TranscriptDataSource = ReferenceDataSource.of(new File(FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19).toPath());
        final Map<String, GencodeFuncotationFactory.MappedTranscriptIdInfo> muc16TranscriptIdMap = GencodeFuncotationFactory. createTranscriptIdMap(muc16TranscriptDataSource);

        final SequenceComparison seqComp =
                GencodeFuncotationFactory.createSequenceComparison(
                        variantContext,
                        altAllele,
                        referenceContext,
                        transcript,
                        exonPositionList,
                        muc16TranscriptIdMap,
                        false);

        final GencodeFuncotation.VariantClassification varClass = GencodeFuncotationFactory.createVariantClassification(
                variantContext,
                altAllele,
                variantType,
                exon,
                transcript.getExons().size(),
                seqComp,
                spliceSiteWindow
        );

        Assert.assertEquals( varClass, expectedVariantClassification );
    }

    @Test (dataProvider = "provideMuc16SnpDataForGetVariantClassificationWithOutOfCdsData")
    void testMuc16SnpCreateFuncotations(final int chromosomeNumber,
                                        final int start,
                                        final int end,
                                        final GencodeFuncotation.VariantType expectedVariantType,
                                        final String ref,
                                        final String alt,
                                        final GencodeFuncotation.VariantClassification expectedVariantClassification) {

        final String contig = "chr" + Integer.toString(chromosomeNumber);
        final SimpleInterval variantInterval = new SimpleInterval( contig, start, end );

        final Allele refAllele = Allele.create(ref, true);
        final Allele altAllele = Allele.create(alt);

        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder(
                FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(),
                contig,
                start,
                end,
                Arrays.asList(refAllele, altAllele)
        );
        final VariantContext variantContext = variantContextBuilder.make();

        // Get our gene feature iterator:
        final CloseableTribbleIterator<GencodeGtfFeature> gtfFeatureIterator;
        try {
            gtfFeatureIterator = gencodeHg19FeatureReader.query(contig, start, end);
        }
        catch (final IOException ex) {
            throw new GATKException("Could not finish the test!", ex);
        }

        // Get the gene.
        // We know the first gene is the right one - the gene in question is the MUC16 gene:
        final GencodeGtfGeneFeature gene = (GencodeGtfGeneFeature) gtfFeatureIterator.next();
        final ReferenceContext referenceContext = new ReferenceContext(refDataSourceHg19Ch19, variantInterval );
        final Set<String> requestedTranscriptIds = getValidTranscriptsForGene("MUC16");

        // Create a factory for our funcotations:
        try (final GencodeFuncotationFactory funcotationFactory = new GencodeFuncotationFactory(
                IOUtils.getPath(FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19),
                "VERSION",
                GencodeFuncotationFactory.DEFAULT_NAME,
                FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_DEFAULT_VALUE,
                requestedTranscriptIds,
                new LinkedHashMap<>(), createFeatureInputForMuc16Ds(GencodeFuncotationFactory.DEFAULT_NAME), "HG19")) {

            // Generate our funcotations:
            final List<Feature> featureList = new ArrayList<>();
            featureList.add( gene );
            final List<Funcotation> funcotations = funcotationFactory.createFuncotationsOnVariant(variantContext, referenceContext, featureList);

            // Make sure we get what we expected:
            Assert.assertEquals(funcotations.size(), 1);

            final GencodeFuncotation funcotation = (GencodeFuncotation)funcotations.get(0);

            Assert.assertEquals(funcotation.getVariantClassification(), expectedVariantClassification);
            Assert.assertEquals(funcotation.getVariantType(), expectedVariantType);
        }
    }

    /*
     * Test the case where we have no basic transcripts at all. In these cases, we should get a single IGR annotation
     */
    @Test(dataProvider = "provideForCreateNonBasicFuncotations")
    void createNonBasicFuncotations(final int start, final int end) {

        final String contig = "chr19";

        final SimpleInterval variantInterval = new SimpleInterval( contig, start, end );

        final Allele refAllele = Allele.create("C", true);
        final Allele altAllele = Allele.create("G", false);

        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder(
                FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(),
                contig,
                start,
                end,
                Arrays.asList(refAllele, altAllele)
        );
        final VariantContext variantContext = variantContextBuilder.make();

        // Get our gene feature iterator:
        final CloseableTribbleIterator<GencodeGtfFeature> gtfFeatureIterator;
        try {
            gtfFeatureIterator = muc16NonBasicFeatureReader.query(contig, start, end);
        }
        catch (final IOException ex) {
            throw new GATKException("Could not finish the test!", ex);
        }

        // Get the gene.
        // We know the first gene is the right one - the gene in question is the MUC16 gene:
        final GencodeGtfGeneFeature gene = (GencodeGtfGeneFeature) gtfFeatureIterator.next();
        final ReferenceContext referenceContext = new ReferenceContext(refDataSourceHg19Ch19, variantInterval );

        // Create a factory for our funcotations:
        try (final GencodeFuncotationFactory funcotationFactory = new GencodeFuncotationFactory(
                                                                    IOUtils.getPath(FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19),
                                                                    "VERSION",
                                                                    GencodeFuncotationFactory.DEFAULT_NAME,
                                                                    FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_DEFAULT_VALUE,
                                                                    new HashSet<>(),
                                                                    new LinkedHashMap<>(),
                                                                    createFeatureInputForMuc16Ds(GencodeFuncotationFactory.DEFAULT_NAME),
                                                                    "HG19")) {

            // Generate our funcotations:
            final List<Feature> featureList = new ArrayList<>();
            featureList.add( gene );
            final List<Funcotation> funcotations = funcotationFactory.createFuncotationsOnVariant(variantContext, referenceContext, featureList);

            // Make sure we get what we expected (a single IGR funcotation):
            Assert.assertEquals(funcotations.size(), 1);
            Assert.assertTrue(funcotations.get(0) instanceof GencodeFuncotation);

            final GencodeFuncotation gencodeFuncotation = (GencodeFuncotation)funcotations.get(0);
            Assert.assertEquals(gencodeFuncotation.getVariantClassification(), GencodeFuncotation.VariantClassification.IGR);
            Assert.assertNull(gencodeFuncotation.getHugoSymbol());
        }
    }

    @Test (dataProvider = "provideDataForCreateFuncotations")
    public void testCreateFuncotations(final String expectedGeneName,
                                final int chromosomeNumber,
                                final int start,
                                final int end,
                                final GencodeFuncotation.VariantClassification expectedVariantClassification,
                                final GencodeFuncotation.VariantType expectedVariantType,
                                final String ref,
                                final String alt,
                                final String expectedGenomeChange,
                                final String expectedStrand,
                                final String expectedCDnaChange,
                                final String expectedCodonChange,
                                final String expectedProteinChange,
                                final String referenceFileName,
                                final FeatureReader<GencodeGtfFeature> featureReader,
                                final ReferenceDataSource referenceDataSource,
                                final String transcriptFastaFile, final String transcriptGtfFile) {

        final String contig = "chr" + Integer.toString(chromosomeNumber);
        final SimpleInterval variantInterval = new SimpleInterval( contig, start, end );

        final Allele refAllele = Allele.create(ref, true);
        final Allele altAllele = Allele.create(alt);

        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder(
                referenceFileName,
                contig,
                start,
                end,
                Arrays.asList(refAllele, altAllele)
        );
        final VariantContext variantContext = variantContextBuilder.make();

        // Get our gene feature iterator:
        final CloseableTribbleIterator<GencodeGtfFeature> gtfFeatureIterator;
        try {
            gtfFeatureIterator = featureReader.query(contig, start, end);
        }
        catch (final IOException ex) {
            throw new GATKException("Could not finish the test!", ex);
        }

        // Get the gene.
        // We know the first gene is the right one - the gene in question is the MUC16 gene:
        final GencodeGtfGeneFeature gene = (GencodeGtfGeneFeature) gtfFeatureIterator.next();
        final ReferenceContext referenceContext = new ReferenceContext(referenceDataSource, variantInterval );
        final Set<String> requestedTranscriptIds = getValidTranscriptsForGene(expectedGeneName);

        // Run this test with flanking turned on for both ends, to make sure that we don't get
        // any false positive flank calls. The flanks feature is tested independently below in
        // testCreateFuncotationsWithFlanks().
        final FlankSettings flankSettingsForTesting = new FlankSettings(5000, 5000);

        // Create a factory for our funcotations:
        try (final GencodeFuncotationFactory funcotationFactory = new GencodeFuncotationFactory(
                IOUtils.getPath(transcriptFastaFile),
                "VERSION",
                GencodeFuncotationFactory.DEFAULT_NAME,
                FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_DEFAULT_VALUE,
                requestedTranscriptIds,
                new LinkedHashMap<>(),
                new FeatureInput<>(transcriptGtfFile, GencodeFuncotationFactory.DEFAULT_NAME, Collections.emptyMap()),
                flankSettingsForTesting,
                "TEST")) {

            final List<Feature> featureList = new ArrayList<>();
            featureList.add( gene );
            final List<Funcotation> funcotations = funcotationFactory.createFuncotationsOnVariant(variantContext, referenceContext, featureList);

            // Make sure we get what we expected:
            Assert.assertEquals(funcotations.size(), 1);

            final GencodeFuncotation funcotation = (GencodeFuncotation)funcotations.get(0);

            final boolean geneNameCorrect              = Objects.equals( funcotation.getHugoSymbol(), expectedGeneName );
            final boolean variantClassificationCorrect = Objects.equals( funcotation.getVariantClassification(), expectedVariantClassification );
            final boolean variantTypeCorrect           = Objects.equals( funcotation.getVariantType(), expectedVariantType );
            final boolean genomeChangeCorrect          = Objects.equals( funcotation.getGenomeChange(), expectedGenomeChange );
            final boolean strandCorrect                = Objects.equals( funcotation.getTranscriptStrand(), expectedStrand );
            final boolean cDnaChangeCorrect            = Objects.equals( funcotation.getcDnaChange(), expectedCDnaChange );
            final boolean codonChangeCorrect           = Objects.equals( funcotation.getCodonChange(), expectedCodonChange );
            final boolean proteinChangeCorrect         = Objects.equals( funcotation.getProteinChange(), expectedProteinChange );

            final StringBuilder errorMessageStringBuilder = new StringBuilder();

            if (!geneNameCorrect) {
                errorMessageStringBuilder.append("\n\tGene Name is not correct!\n\t\tExpected:  [");
                errorMessageStringBuilder.append(expectedGeneName);
                errorMessageStringBuilder.append("]\n\t\tBut found: [");
                errorMessageStringBuilder.append(funcotation.getHugoSymbol());
                errorMessageStringBuilder.append("]");
            }
            if (!variantClassificationCorrect) {
                errorMessageStringBuilder.append("\n\tVariant Classification is not correct!\n\t\tExpected:  [");
                errorMessageStringBuilder.append(expectedVariantClassification);
                errorMessageStringBuilder.append("]\n\t\tBut found: [");
                errorMessageStringBuilder.append(funcotation.getVariantClassification());
                errorMessageStringBuilder.append("]");
            }
            if (!variantTypeCorrect) {
                errorMessageStringBuilder.append("\n\tVariant Type is not correct!\n\t\tExpected:  [");
                errorMessageStringBuilder.append(expectedVariantType);
                errorMessageStringBuilder.append("]\n\t\tBut found: [");
                errorMessageStringBuilder.append(funcotation.getVariantType());
                errorMessageStringBuilder.append("]");
            }
            if (!genomeChangeCorrect) {
                errorMessageStringBuilder.append("\n\tGenome Change is not correct!\n\t\tExpected:  [");
                errorMessageStringBuilder.append(expectedGenomeChange);
                errorMessageStringBuilder.append("]\n\t\tBut found: [");
                errorMessageStringBuilder.append(funcotation.getGenomeChange());
                errorMessageStringBuilder.append("]");
            }
            if (!strandCorrect) {
                errorMessageStringBuilder.append("\n\tStrand is not correct!\n\t\tExpected:  [");
                errorMessageStringBuilder.append(expectedStrand);
                errorMessageStringBuilder.append("]\n\t\tBut found: [");
                errorMessageStringBuilder.append(funcotation.getTranscriptStrand());
                errorMessageStringBuilder.append("]");
            }
            if (!cDnaChangeCorrect) {
                errorMessageStringBuilder.append("\n\tCDna Change is not correct!\n\t\tExpected:  [");
                errorMessageStringBuilder.append(expectedCDnaChange);
                errorMessageStringBuilder.append("]\n\t\tBut found: [");
                errorMessageStringBuilder.append(funcotation.getcDnaChange());
                errorMessageStringBuilder.append("]");
            }
            if (!codonChangeCorrect) {
                errorMessageStringBuilder.append("\n\tCodon Change is not correct!\n\t\tExpected:  [");
                errorMessageStringBuilder.append(expectedCodonChange);
                errorMessageStringBuilder.append("]\n\t\tBut found: [");
                errorMessageStringBuilder.append(funcotation.getCodonChange());
                errorMessageStringBuilder.append("]");
            }
            if (!proteinChangeCorrect) {
                errorMessageStringBuilder.append("\n\tProtein Change is not correct!\n\t\tExpected:  [");
                errorMessageStringBuilder.append(expectedProteinChange);
                errorMessageStringBuilder.append("]\n\t\tBut found: [");
                errorMessageStringBuilder.append(funcotation.getProteinChange());
                errorMessageStringBuilder.append("]");
            }
            errorMessageStringBuilder.append("\n");
            errorMessageStringBuilder.append(funcotation.getAnnotationTranscript());

            Assert.assertTrue(
                    geneNameCorrect && variantClassificationCorrect && variantTypeCorrect &&
                    genomeChangeCorrect && strandCorrect && cDnaChangeCorrect && codonChangeCorrect && proteinChangeCorrect,
                    errorMessageStringBuilder.toString() + "\n"
            );

            Assert.assertNotNull(funcotation.getMetadata(), "Metadata was null.");
            Assert.assertTrue(funcotation.getMetadata().retrieveAllHeaderInfo().size() > 0, "Metadata was empty.");

            final Set<String> metaDataFields =  funcotation.getMetadata().retrieveAllHeaderInfo().stream()
                    .map(f -> f.getID()).collect(Collectors.toSet());
            final Set<String> symmetricDifference = Sets.symmetricDifference(metaDataFields, funcotation.getFieldNames());

            Assert.assertEquals(symmetricDifference.size(), 0, "Metadata fields did not match exactly the funcotation field names: " +
                symmetricDifference.stream().collect(Collectors.joining(", ")));
        }
    }

    @Test(dataProvider = "provideForTestGetTranscriptEndPaddingBases")
    void testGetTranscriptEndPaddingBases(final VariantContext variant,
                                          final Allele altAllele,
                                          final List<? extends htsjdk.samtools.util.Locatable> exonPositionList,
                                          final ReferenceContext reference,
                                          final String expected) {
        Assert.assertEquals(GencodeFuncotationFactory.getTranscriptEndPaddingBases(variant, altAllele, exonPositionList, reference), expected);
    }

    @Test(dataProvider = "provideForTestCreateGencodeFuncotationBuilderWithTrivialFieldsPopulated")
    void testCreateGencodeFuncotationBuilderWithTrivialFieldsPopulated(final VariantContext variant,
                                                                       final Allele altAllele,
                                                                       final GencodeGtfGeneFeature gtfFeature,
                                                                       final GencodeGtfTranscriptFeature transcript ) {

        final GencodeFuncotationBuilder builder =
                GencodeFuncotationFactory.createGencodeFuncotationBuilderWithTrivialFieldsPopulated(variant, altAllele,
                        transcript, "TEST");
        final GencodeFuncotation gf = builder.gencodeFuncotation;

        // Ultra-trivial checks:
        Assert.assertEquals( gf.getRefAllele(), variant.getReference().getBaseString() );
        Assert.assertEquals( gf.getTranscriptStrand(), transcript.getGenomicStrand().encode() );
        Assert.assertEquals( gf.getHugoSymbol(), gtfFeature.getGeneName() );
        Assert.assertEquals( gf.getNcbiBuild(), "TEST" );
        Assert.assertEquals( gf.getChromosome(), gtfFeature.getChromosomeName() );
        Assert.assertEquals( gf.getStart(), variant.getStart() );
        Assert.assertEquals( gf.getGeneTranscriptType(), transcript.getTranscriptType() );
        Assert.assertEquals( gf.getEnd(), variant.getStart() + altAllele.length() - 1 );
        Assert.assertEquals( gf.getTumorSeqAllele2(), altAllele.getBaseString() );
        Assert.assertEquals( gf.getAnnotationTranscript(), transcript.getTranscriptId() );
        Assert.assertEquals( gf.getLocusLevel(), Integer.valueOf(transcript.getLocusLevel().toString()) );
        Assert.assertEquals( gf.getEnd(), variant.getEnd() );

        // Still simple, but computational, checks:
        Assert.assertEquals( gf.getVariantType(), GencodeFuncotationFactory.getVariantType(variant.getReference(), altAllele));
        Assert.assertEquals( gf.getGenomeChange(), GencodeFuncotationFactory.getGenomeChangeString(variant, altAllele) );
        Assert.assertEquals( gf.getApprisRank(), GencodeFuncotationFactory.getApprisRank( transcript ) );
        Assert.assertTrue( gf.getTranscriptLength() == transcript.getExons().stream().mapToInt(Locatable::getLengthOnReference).sum() );
    }

    @Test(dataProvider = "provideForTestIs5PrimeUtr")
    void testIs5PrimeUtr(final GencodeGtfUTRFeature utr, final GencodeGtfTranscriptFeature transcript, final boolean expected) {
        Assert.assertEquals(GencodeFuncotationFactory.is5PrimeUtr( utr, transcript), expected);
    }

    @Test (dataProvider = "provideDataForTestIsVariantInCodingRegion")
    void testIsVariantInCodingRegion(final GencodeFuncotation.VariantClassification varClass, final GencodeFuncotation.VariantClassification secondaryVarClass, final boolean expected) {
        Assert.assertEquals( GencodeFuncotationFactory.isVariantInCodingRegion(varClass, secondaryVarClass), expected );
    }

    @Test ( dataProvider = "provideForCreateDefaultFuncotationsOnProblemVariant")
    void testCreateDefaultFuncotationsOnProblemVariant(final VariantContext variant,
                                                       final Allele altAllele,
                                                       final GencodeGtfGeneFeature gtfFeature,
                                                       final ReferenceContext reference,
                                                       final GencodeGtfTranscriptFeature transcript,
                                                       final String version,
                                                       final String dataSourceName,
                                                       final GencodeFuncotation expected) {

        final GencodeFuncotation funcotation = GencodeFuncotationFactory.createDefaultFuncotationsOnProblemVariant(variant, altAllele, reference, transcript, version, dataSourceName, "TEST");
        Assert.assertEquals( funcotation, expected );
    }

    @Test ( dataProvider = "provideDataForTestGencodeFuncotationComparatorUnitTest")
    void testGencodeFuncotationComparatorUnitTest( final Comparator<GencodeFuncotation> comparator,
                                                   final GencodeFuncotation a,
                                                   final GencodeFuncotation b,
                                                   final int expected) {
        Assert.assertEquals( comparator.compare(a,b), expected );
    }

    @Test (dataProvider = "provideDataForTestCalculateGcContent")
    void testCalculateGcContent(final Allele refAllele,
                                final Allele altAllele,
                                final ReferenceContext referenceContext,
                                final int windowSize,
                                final double expected) {
        Assert.assertEquals( GencodeFuncotationFactory.calculateGcContent( refAllele, altAllele, referenceContext, windowSize ), expected, doubleEqualsEpsilon);
    }

    @DataProvider
    public Object[][] provideSortingOfUserRequestedTranscripts() {
        final String transcriptId1 = "ENST0000001.4";
        final String transcriptId2 = "ENST0000002.3";
        final String transcriptId3 = "ENST0000003.8";
        return new Object[][] {
                {       Arrays.asList(
                            createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1, GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                            createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2, GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                            createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId3, GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.MISSENSE, 5000)
                        ),
                        Collections.singleton(transcriptId2),
                        transcriptId2
                },
                {       Arrays.asList(
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1, GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2, GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId3, GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.MISSENSE, 5000)
                        ),
                        Collections.singleton(transcriptId1),
                        transcriptId1
                },
                {       Arrays.asList(
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1, GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.MISSENSE, 5000),
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2, GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId3, GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.MISSENSE, 5000)
                        ),
                        new HashSet<>(Arrays.asList(transcriptId2, transcriptId1)),
                        transcriptId1
                },
                {       Arrays.asList(
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1, GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.MISSENSE, 5000),
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2, GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000),
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId3, GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.MISSENSE, 5000)
                        ),
                        new HashSet<>(Arrays.asList(transcriptId2, transcriptId1)),
                        transcriptId1
                },
                {       Arrays.asList(
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1, GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.MISSENSE, 5000),
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2, GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, 5000)
                        ),
                        new HashSet<>(Arrays.asList(transcriptId2, transcriptId1)),
                        transcriptId1
                },
                {       Arrays.asList(
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1, GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.MISSENSE, 5000),
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2, GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.MISSENSE, 5000)
                        ),
                        new HashSet<>(Arrays.asList(transcriptId2, transcriptId1)),
                        transcriptId1
                },
                {       Arrays.asList(
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1, GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 2, GencodeFuncotation.VariantClassification.MISSENSE, 5000),
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2, GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.MISSENSE, 5000)
                        ),
                        new HashSet<>(Arrays.asList(transcriptId2, transcriptId1)),
                        transcriptId2
                },
                {       Arrays.asList(
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1, GencodeGTFFieldConstants.FeatureTag.APPRIS_CANDIDATE, 1, GencodeFuncotation.VariantClassification.MISSENSE, 5000),
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2, GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.MISSENSE, 5000)
                        ),
                        new HashSet<>(Arrays.asList(transcriptId2, transcriptId1)),
                        transcriptId2
                },
                {       Arrays.asList(
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1, GencodeGTFFieldConstants.FeatureTag.APPRIS_CANDIDATE_HIGHEST_SCORE, 1, GencodeFuncotation.VariantClassification.MISSENSE, 5000),
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2, GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.MISSENSE, 5000)
                        ),
                        new HashSet<>(Arrays.asList(transcriptId2, transcriptId1)),
                        transcriptId2
                },
                {       Arrays.asList(
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1, GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.MISSENSE, 2000),
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2, GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.MISSENSE, 5000)
                        ),
                        new HashSet<>(Arrays.asList(transcriptId2, transcriptId1)),
                        transcriptId2
                },
                {       Arrays.asList(
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId1, GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.MISSENSE, 2000),
                                createFuncotationForTestGencodeFuncotationComparatorUnitTest(transcriptId2, GencodeGTFFieldConstants.FeatureTag.APPRIS_PRINCIPAL, 1, GencodeFuncotation.VariantClassification.MISSENSE, 5000)
                        ),Collections.singleton(transcriptId3),
                        transcriptId2
                }
        };
    }

    /**
     * Also, tests some basic sorting independent of the user defined list.
     */
    @Test(dataProvider = "provideSortingOfUserRequestedTranscripts")
    public void testSortingOfUserRequestedTranscripts(final List<GencodeFuncotation> gencodeFuncotations, final Set<String> userRequestedTranscripts, final String gtFirstTranscript) {

        final List<TranscriptSelectionMode> transcriptSelectionModes = Arrays.asList(TranscriptSelectionMode.ALL, TranscriptSelectionMode.BEST_EFFECT, TranscriptSelectionMode.CANONICAL);
        for (final TranscriptSelectionMode transcriptSelectionMode : transcriptSelectionModes) {
            final Comparator<GencodeFuncotation> comparator = transcriptSelectionMode.getComparator(userRequestedTranscripts);
            gencodeFuncotations.sort(comparator);
            Assert.assertEquals(gencodeFuncotations.get(0).getAnnotationTranscript(), gtFirstTranscript, " Failed on " + transcriptSelectionMode.toString());
        }
    }


    /**
     * This test (of {@link GencodeFuncotationFactory#createFuncotationsOnVariant}) makes sure that if multiple gene features are detected, there is still only one transcript returned
     *  when in the BEST_EFFECT or CANONICAL selection mode.  And tests that ALL has multiple transcripts returned.
     *  The selected test data will ensure that ALL returns more than one transcript.
     * @throws IOException
     */
    @Test
    public void testMultipleGeneFeaturesOnlyProduceOneTranscript() throws IOException {

        // Note the "chr" here to make this work.
        final SimpleInterval variantInterval = new SimpleInterval("chr3", 2944600, 2944600);
        final ReferenceContext referenceContext = new ReferenceContext(refDataSourceHg19Ch3, variantInterval);
        final VariantContext vc = new VariantContextBuilder()
                .alleles(Arrays.asList(Allele.create("T", true), Allele.create("AT", false)))
                .chr(variantInterval.getContig()).start(variantInterval.getStart()).stop(variantInterval.getEnd())
                .make();

        final List<Feature> gencodeFeatures = getCntn4Features();

        try (final GencodeFuncotationFactory funcotationFactory = new GencodeFuncotationFactory(
                IOUtils.getPath(CNTN4_GENCODE_TRANSCRIPT_FASTA_FILE),
                "VERSION",
                GencodeFuncotationFactory.DEFAULT_NAME,
                TranscriptSelectionMode.CANONICAL,
                Collections.emptySet(),
                new LinkedHashMap<>(), createFeatureInputForCntn4Ds(GencodeFuncotationFactory.DEFAULT_NAME),
                "TEST")) {
            final List<Funcotation> gencodeFuncotations = funcotationFactory.createFuncotationsOnVariant(vc, referenceContext, gencodeFeatures);
            Assert.assertEquals(gencodeFuncotations.size(), 1);
        }

        try (final GencodeFuncotationFactory funcotationFactory = new GencodeFuncotationFactory(
                IOUtils.getPath(CNTN4_GENCODE_TRANSCRIPT_FASTA_FILE),
                "VERSION",
                GencodeFuncotationFactory.DEFAULT_NAME,
                TranscriptSelectionMode.BEST_EFFECT,
                Collections.emptySet(),
                new LinkedHashMap<>(), createFeatureInputForCntn4Ds(GencodeFuncotationFactory.DEFAULT_NAME),
                "TEST")) {
            final List<Funcotation> gencodeFuncotations = funcotationFactory.createFuncotationsOnVariant(vc, referenceContext, gencodeFeatures);
            Assert.assertEquals(gencodeFuncotations.size(), 1);
        }

        try (final GencodeFuncotationFactory funcotationFactory = new GencodeFuncotationFactory(
                IOUtils.getPath(CNTN4_GENCODE_TRANSCRIPT_FASTA_FILE),
                "VERSION",
                GencodeFuncotationFactory.DEFAULT_NAME,
                TranscriptSelectionMode.ALL,
                Collections.emptySet(),
                new LinkedHashMap<>(), createFeatureInputForCntn4Ds(GencodeFuncotationFactory.DEFAULT_NAME),
                "TEST")) {
            final List<Funcotation> gencodeFuncotations = funcotationFactory.createFuncotationsOnVariant(vc, referenceContext, gencodeFeatures);
            Assert.assertTrue(gencodeFuncotations.size() > 1);
        }
    }

    private List<Feature> getCntn4Features() throws IOException {
        final GencodeGtfCodec gencodeGtfCodec = new GencodeGtfCodec();
        Assert.assertTrue(gencodeGtfCodec.canDecode(CNTN4_GENCODE_ANNOTATIONS_FILE_NAME));

        final List<Feature> gencodeFeatures = new ArrayList<>();
        try (BufferedInputStream bufferedInputStream =
                     new BufferedInputStream(
                             new FileInputStream(CNTN4_GENCODE_ANNOTATIONS_FILE_NAME)
                     )
        ) {
            // Get the line iterator:
            final LineIterator lineIterator = gencodeGtfCodec.makeSourceFromStream(bufferedInputStream);

            // Get the header (required for the read to work correctly):
            gencodeGtfCodec.readHeader(lineIterator);

            while (lineIterator.hasNext()) {
                gencodeFeatures.add(gencodeGtfCodec.decode(lineIterator));
            }
            Assert.assertTrue(gencodeFeatures.size() > 1);
        }
        return gencodeFeatures;
    }

    private static FeatureInput<? extends Feature> createFeatureInputForMuc16Ds(final String dsName) {
        return new FeatureInput<>(FuncotatorTestConstants.GENCODE_DATA_SOURCE_FASTA_PATH_HG19, dsName, Collections.emptyMap());
    }

    private static FeatureInput<? extends Feature> createFeatureInputForCntn4Ds(final String dsName) {
        return new FeatureInput<>(CNTN4_GENCODE_ANNOTATIONS_FILE_NAME, dsName, Collections.emptyMap());
    }

    @DataProvider
    public Object[][] provideSimpleGenesField() {
        return new Object[][] {

                // Two genes found
                {new VariantContextBuilder()
                        .chr("chr3").start(2613200).stop(3200000)
                        .attribute(VCFConstants.END_KEY, 3200000)
                        .alleles(Arrays.asList(
                                Allele.create("T", true),
                                Allele.create(SimpleSVType.createBracketedSymbAlleleString(StructuralVariantType.DEL.name()), false))
                        )
                        .make(), "CNTN4,CNTN4-AS1"},

                // Three genes found.  Note sorting is alphabetical, not order seen in the genome.  This is how users
                //  expect it.
                {new VariantContextBuilder()
                        .chr("chr3").start(2100000).stop(3200000)
                        .attribute(VCFConstants.END_KEY, 3200000)
                        .alleles(Arrays.asList(
                                Allele.create("T", true),
                                Allele.create(SimpleSVType.createBracketedSymbAlleleString(StructuralVariantType.DEL.name()), false))
                        )
                        .make(), "CNTN4,CNTN4-AS1,CNTN4-AS2"},

                // We use no call for copy neutral
                {new VariantContextBuilder()
                        .chr("chr3").start(2100000).stop(3200000)
                        .attribute(VCFConstants.END_KEY, 3200000)
                        .alleles(Arrays.asList(
                                Allele.create("T", true),
                                AnnotatedIntervalToSegmentVariantContextConverter.COPY_NEUTRAL_ALLELE
                        ))
                        .make(), "CNTN4,CNTN4-AS1,CNTN4-AS2"}
        };
    }

    @Test(dataProvider = "provideSimpleGenesField")
    public void testSimpleGenesField(final VariantContext vc, final String genesTruth) throws IOException {

        try (final GencodeFuncotationFactory funcotationFactory = new GencodeFuncotationFactory(
                IOUtils.getPath(CNTN4_GENCODE_TRANSCRIPT_FASTA_FILE),
                "VERSION",
                GencodeFuncotationFactory.DEFAULT_NAME,
                TranscriptSelectionMode.CANONICAL,
                Collections.emptySet(),
                new LinkedHashMap<>(), createFeatureInputForCntn4Ds(GencodeFuncotationFactory.DEFAULT_NAME), "TEST")) {
            final ReferenceContext referenceContext = new ReferenceContext(refDataSourceHg19Ch3, new SimpleInterval(vc) );

            final FeatureContext featureContext = FuncotatorTestUtils.createFeatureContext(
                    Collections.singletonList(funcotationFactory), "TEST", new SimpleInterval(vc),
                    0,0,0,null);

            final List<Funcotation> funcotations = funcotationFactory.createFuncotations(
                    vc, referenceContext, featureContext);
            Assert.assertEquals(funcotations.size(), 1);
            Assert.assertEquals(funcotations.get(0).getField(funcotationFactory.getName() + "_" + funcotationFactory.getVersion() + "_" + "genes"), genesTruth);
        }
    }

    @DataProvider
    public Object[][] provideForTestIsFivePrimeFlankAndIsThreePrimeFlank() {
        final GencodeGtfTranscriptFeature positiveStrandTranscript =
                FuncotatorTestUtils.createArtificialGencodeGtfTranscriptFeatureForTesting("chr1", 1000, 2000, Strand.POSITIVE);

        final GencodeGtfTranscriptFeature negativeStrandTranscript =
                FuncotatorTestUtils.createArtificialGencodeGtfTranscriptFeatureForTesting("chr1", 1000, 2000, Strand.NEGATIVE);

        return new Object[][] {
                // Variant, transcript, 5' flank size, 3' flank size, expected isFivePrimeFlank, expected isThreePrimeFlank

                // SNP, + strand, in/around the 5' flank:
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 988, 988, "A", "G"),
                        positiveStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 989, 989, "A", "G"),
                        positiveStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 990, 990, "A", "G"),
                        positiveStrandTranscript, 10, 10, true, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 995, 995, "A", "G"),
                        positiveStrandTranscript, 10, 10, true, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 999, 999, "A", "G"),
                        positiveStrandTranscript, 10, 10, true, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 1000, 1000, "A", "G"),
                        positiveStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 1001, 1001, "A", "G"),
                        positiveStrandTranscript, 10, 10, false, false },

                // SNP, + strand, in/around the 3' flank:
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 1999, 1999, "A", "G"),
                        positiveStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2000, 2000, "A", "G"),
                        positiveStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2001, 2001, "A", "G"),
                        positiveStrandTranscript, 10, 10, false, true },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2005, 2005, "A", "G"),
                        positiveStrandTranscript, 10, 10, false, true },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2010, 2010, "A", "G"),
                        positiveStrandTranscript, 10, 10, false, true },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2011, 2011, "A", "G"),
                        positiveStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2012, 2012, "A", "G"),
                        positiveStrandTranscript, 10, 10, false, false },

                // SNP, - strand, in/around the 3' flank:
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 988, 988, "A", "G"),
                        negativeStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 989, 989, "A", "G"),
                        negativeStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 990, 990, "A", "G"),
                        negativeStrandTranscript, 10, 10, false, true },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 995, 995, "A", "G"),
                        negativeStrandTranscript, 10, 10, false, true },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 999, 999, "A", "G"),
                        negativeStrandTranscript, 10, 10, false, true },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 1000, 1000, "A", "G"),
                        negativeStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 1001, 1001, "A", "G"),
                        negativeStrandTranscript, 10, 10, false, false },

                // SNP, - strand, in/around the 5' flank:
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 1999, 1999, "A", "G"),
                        negativeStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2000, 2000, "A", "G"),
                        negativeStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2001, 2001, "A", "G"),
                        negativeStrandTranscript, 10, 10, true, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2005, 2005, "A", "G"),
                        negativeStrandTranscript, 10, 10, true, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2010, 2010, "A", "G"),
                        negativeStrandTranscript, 10, 10, true, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2011, 2011, "A", "G"),
                        negativeStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2012, 2012, "A", "G"),
                        negativeStrandTranscript, 10, 10, false, false },

                // Insertion, + strand, in/around the 5' flank:
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 988, 988, "A", "GT"),
                        positiveStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 989, 989, "A", "GT"),
                        positiveStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 990, 990, "A", "GT"),
                        positiveStrandTranscript, 10, 10, true, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 995, 995, "A", "GT"),
                        positiveStrandTranscript, 10, 10, true, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 999, 999, "A", "GT"),
                        positiveStrandTranscript, 10, 10, true, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 1000, 1000, "A", "GT"),
                        positiveStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 1001, 1001, "A", "GT"),
                        positiveStrandTranscript, 10, 10, false, false },

                // Insertion, + strand, in/around the 3' flank:
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 1999, 1999, "A", "GT"),
                        positiveStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2000, 2000, "A", "GT"),
                        positiveStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2001, 2001, "A", "GT"),
                        positiveStrandTranscript, 10, 10, false, true },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2005, 2005, "A", "GT"),
                        positiveStrandTranscript, 10, 10, false, true },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2010, 2010, "A", "GT"),
                        positiveStrandTranscript, 10, 10, false, true },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2011, 2011, "A", "GT"),
                        positiveStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2012, 2012, "A", "GT"),
                        positiveStrandTranscript, 10, 10, false, false },

                // Insertion, - strand, in/around the 3' flank:
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 988, 988, "A", "GT"),
                        negativeStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 989, 989, "A", "GT"),
                        negativeStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 990, 990, "A", "GT"),
                        negativeStrandTranscript, 10, 10, false, true },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 995, 995, "A", "GT"),
                        negativeStrandTranscript, 10, 10, false, true },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 999, 999, "A", "GT"),
                        negativeStrandTranscript, 10, 10, false, true },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 1000, 1000, "A", "GT"),
                        negativeStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 1001, 1001, "A", "GT"),
                        negativeStrandTranscript, 10, 10, false, false },

                // Insertion, - strand, in/around the 5' flank:
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 1999, 1999, "A", "GT"),
                        negativeStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2000, 2000, "A", "GT"),
                        negativeStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2001, 2001, "A", "GT"),
                        negativeStrandTranscript, 10, 10, true, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2005, 2005, "A", "GT"),
                        negativeStrandTranscript, 10, 10, true, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2010, 2010, "A", "GT"),
                        negativeStrandTranscript, 10, 10, true, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2011, 2011, "A", "GT"),
                        negativeStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2012, 2012, "A", "GT"),
                        negativeStrandTranscript, 10, 10, false, false },

                // Deletion, + strand, in/around the 5' flank:
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 986, 988, "ACT", "G"),
                        positiveStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 987, 989, "ACT", "G"),
                        positiveStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 988, 990, "ACT", "G"),
                        positiveStrandTranscript, 10, 10, true, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 989, 991, "ACT", "G"),
                        positiveStrandTranscript, 10, 10, true, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 990, 992, "ACT", "G"),
                        positiveStrandTranscript, 10, 10, true, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 994, 996, "ACT", "G"),
                        positiveStrandTranscript, 10, 10, true, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 997, 999, "ACT", "G"),
                        positiveStrandTranscript, 10, 10, true, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 998, 1000, "ACT", "G"),
                        positiveStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 999, 1001, "ACT", "G"),
                        positiveStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 1000, 1002, "ACT", "G"),
                        positiveStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 1001, 1003, "ACT", "G"),
                        positiveStrandTranscript, 10, 10, false, false },

                // Deletion, + strand, in/around the 3' flank:
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 1997, 1999, "ACT", "G"),
                        positiveStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 1998, 2000, "ACT", "G"),
                        positiveStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 1999, 2001, "ACT", "G"),
                        positiveStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2000, 2002, "ACT", "G"),
                        positiveStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2001, 2003, "ACT", "G"),
                        positiveStrandTranscript, 10, 10, false, true },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2004, 2006, "ACT", "G"),
                        positiveStrandTranscript, 10, 10, false, true },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2008, 2010, "ACT", "G"),
                        positiveStrandTranscript, 10, 10, false, true },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2009, 2011, "ACT", "G"),
                        positiveStrandTranscript, 10, 10, false, true },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2010, 2012, "ACT", "G"),
                        positiveStrandTranscript, 10, 10, false, true },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2011, 2013, "ACT", "G"),
                        positiveStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2012, 2014, "ACT", "G"),
                        positiveStrandTranscript, 10, 10, false, false },

                // Deletion, - strand, in/around the 3' flank:
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 986, 988, "ACT", "G"),
                        negativeStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 987, 989, "ACT", "G"),
                        negativeStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 988, 990, "ACT", "G"),
                        negativeStrandTranscript, 10, 10, false, true },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 989, 991, "ACT", "G"),
                        negativeStrandTranscript, 10, 10, false, true },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 990, 992, "ACT", "G"),
                        negativeStrandTranscript, 10, 10, false, true },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 994, 996, "ACT", "G"),
                        negativeStrandTranscript, 10, 10, false, true },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 997, 999, "ACT", "G"),
                        negativeStrandTranscript, 10, 10, false, true },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 998, 1000, "ACT", "G"),
                        negativeStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 999, 1001, "ACT", "G"),
                        negativeStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 1000, 1002, "ACT", "G"),
                        negativeStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 1001, 1003, "ACT", "G"),
                        negativeStrandTranscript, 10, 10, false, false },

                // Deletion, - strand, in/around the 5' flank:
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 1997, 1999, "ACT", "G"),
                        negativeStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 1998, 2000, "ACT", "G"),
                        negativeStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 1999, 2001, "ACT", "G"),
                        negativeStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2000, 2002, "ACT", "G"),
                        negativeStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2001, 2003, "ACT", "G"),
                        negativeStrandTranscript, 10, 10, true, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2004, 2006, "ACT", "G"),
                        negativeStrandTranscript, 10, 10, true, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2008, 2010, "ACT", "G"),
                        negativeStrandTranscript, 10, 10, true, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2009, 2011, "ACT", "G"),
                        negativeStrandTranscript, 10, 10, true, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2010, 2012, "ACT", "G"),
                        negativeStrandTranscript, 10, 10, true, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2011, 2013, "ACT", "G"),
                        negativeStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2012, 2014, "ACT", "G"),
                        negativeStrandTranscript, 10, 10, false, false },

                // With flank size == 0, we should never see a true return value:
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 999, 999, "A", "G"),
                        positiveStrandTranscript, 0, 0, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 999, 999, "A", "G"),
                        negativeStrandTranscript, 0, 0, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 1000, 1000, "A", "G"),
                        positiveStrandTranscript, 0, 0, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 1000, 1000, "A", "G"),
                        negativeStrandTranscript, 0, 0, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2001, 2001, "A", "G"),
                        positiveStrandTranscript, 0, 0, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2001, 2001, "A", "G"),
                        negativeStrandTranscript, 0, 0, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2000, 2000, "A", "G"),
                        positiveStrandTranscript, 0, 0, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2000, 2000, "A", "G"),
                        negativeStrandTranscript, 0, 0, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 997, 999, "ACT", "G"),
                        positiveStrandTranscript, 0, 0, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 997, 999, "ACT", "G"),
                        negativeStrandTranscript, 0, 0, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2001, 2003, "ACT", "G"),
                        positiveStrandTranscript, 0, 0, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2001, 2003, "ACT", "G"),
                        negativeStrandTranscript, 0, 0, false, false },

                // Test cases to ensure that the 5' flank size affects only 5' flank detection, and
                // the 3' flank size affects only 3' flank detection:
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 999, 999, "A", "G"),
                        positiveStrandTranscript, 1, 0, true, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 999, 999, "A", "G"),
                        negativeStrandTranscript, 1, 0, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 999, 999, "A", "G"),
                        positiveStrandTranscript, 0, 1, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 999, 999, "A", "G"),
                        negativeStrandTranscript, 0, 1, false, true },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2001, 2001, "A", "G"),
                        positiveStrandTranscript, 1, 0, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2001, 2001, "A", "G"),
                        negativeStrandTranscript, 1, 0, true, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2001, 2001, "A", "G"),
                        positiveStrandTranscript, 0, 1, false, true },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr1", 2001, 2001, "A", "G"),
                        negativeStrandTranscript, 0, 1, false, false },

                // Variants on a different contig from the transcript should never produce a true return value:
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr2", 999, 999, "A", "G"),
                        positiveStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr2", 2001, 2001, "A", "G"),
                        positiveStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr2", 999, 999, "A", "G"),
                        negativeStrandTranscript, 10, 10, false, false },
                { FuncotatorTestUtils.createSimpleVariantContext(null, "chr2", 2001, 2001, "A", "G"),
                        negativeStrandTranscript, 10, 10, false, false },
        };
    }

    @Test(dataProvider = "provideForTestIsFivePrimeFlankAndIsThreePrimeFlank")
    public void testIsFivePrimeFlankAndIsThreePrimeFlank(final VariantContext variant,
                                                         final GencodeGtfTranscriptFeature transcript,
                                                         final int fivePrimeFlankSize,
                                                         final int threePrimeFlankSize,
                                                         final boolean expectedIsFivePrimeFlank,
                                                         final boolean expectedIsThreePrimeFlank) {

        Assert.assertEquals(GencodeFuncotationFactory.isFivePrimeFlank(variant, transcript, fivePrimeFlankSize), expectedIsFivePrimeFlank, "isFivePrimeFlank() returned the wrong value");
        Assert.assertEquals(GencodeFuncotationFactory.isThreePrimeFlank(variant, transcript, threePrimeFlankSize), expectedIsThreePrimeFlank, "isThreePrimeFlank() returned the wrong value");
    }

    @DataProvider
    public Object[][] provideDataForCreateFuncotationsWithFlanks() {
        // MUC16 basic transcripts span: chr19:8959520-9092018
        // MUC16 basic transcripts: chr19:8959520-9092018, chr19:8959522-9003586, chr19:8973992-8977690

        // PIK3CA basic transcripts span: chr3:178866311-178957881
        // PIK3CA basic transcripts: chr3: 178866311-178957881

        final String MUC16_CONTIG = FuncotatorTestConstants.MUC16_HG19_BASIC_TRANSCRIPTS_SPAN.getContig();
        final int MUC16_START = FuncotatorTestConstants.MUC16_HG19_BASIC_TRANSCRIPTS_SPAN.getStart();
        final int MUC16_END = FuncotatorTestConstants.MUC16_HG19_BASIC_TRANSCRIPTS_SPAN.getEnd();

        final String PIK3CA_CONTIG = FuncotatorTestConstants.PIK3CA_HG19_BASIC_TRANSCRIPTS_SPAN.getContig();
        final int PIK3CA_START = FuncotatorTestConstants.PIK3CA_HG19_BASIC_TRANSCRIPTS_SPAN.getStart();
        final int PIK3CA_END = FuncotatorTestConstants.PIK3CA_HG19_BASIC_TRANSCRIPTS_SPAN.getEnd();

        return new Object[][] {
                /* ***********************************
                 * MUC16 Test Cases (reverse strand):
                 *************************************/

                /*
                 * SNP 5' reverse strand POSITIVE test cases:
                 */

                // SNP 1 base past the 5' end, 5' flank size = 3000
                { "MUC16", GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK, MUC16_CONTIG, MUC16_END + 1, MUC16_END + 1, "C", "A", 3000, 0, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // SNP 10 bases past the 5' end, 5' flank size = 10
                { "MUC16", GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK, MUC16_CONTIG, MUC16_END + 10, MUC16_END + 10, "A", "T", 10, 0, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // SNP 10 bases past the 5' end, 5' flank size = 10 and 3' flank size = 9
                { "MUC16", GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK, MUC16_CONTIG, MUC16_END + 10, MUC16_END + 10, "A", "T", 10, 9, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // SNP 5 bases past the 5' end, 5' flank size = 10
                { "MUC16", GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK, MUC16_CONTIG, MUC16_END + 5, MUC16_END + 5, "T", "A", 10, 0, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },

                /*
                 * SNP 5' reverse strand NEGATIVE test cases:
                 */

                // SNP 1 base past the 5' end, 5' flank size = 0
                { null, GencodeFuncotation.VariantClassification.IGR, MUC16_CONTIG, MUC16_END + 1, MUC16_END + 1, "C", "A", 0, 0, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // SNP 1 base past the 5' end, 5' flank size = 0 and 3' flank size = 10
                { null, GencodeFuncotation.VariantClassification.IGR, MUC16_CONTIG, MUC16_END + 1, MUC16_END + 1, "C", "A", 0, 10, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // SNP 10 bases past the 5' end, 5' flank size = 9
                { null, GencodeFuncotation.VariantClassification.IGR, MUC16_CONTIG, MUC16_END + 10, MUC16_END + 10, "A", "T", 9, 0, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // SNP on first base of the 5' end, 5' flank size = 10
                { "MUC16", GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR, MUC16_CONTIG, MUC16_END, MUC16_END, "T", "A", 10, 10, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },

                /*
                 * SNP 3' reverse strand POSITIVE test cases:
                 */

                // SNP 1 base past the 3' end, 3' flank size = 100
                { "MUC16", GencodeFuncotation.VariantClassification.THREE_PRIME_FLANK, MUC16_CONTIG, MUC16_START - 1, MUC16_START - 1, "G", "A", 0, 100, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // SNP 10 bases past the 3' end, 3' flank size = 10
                { "MUC16", GencodeFuncotation.VariantClassification.THREE_PRIME_FLANK, MUC16_CONTIG, MUC16_START - 10, MUC16_START - 10, "A", "T", 0, 10, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // SNP 10 bases past the 3' end, 3' flank size = 10 and 5' flank size = 9
                { "MUC16", GencodeFuncotation.VariantClassification.THREE_PRIME_FLANK, MUC16_CONTIG, MUC16_START - 10, MUC16_START - 10, "A", "T", 9, 10, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // SNP 5 bases past the 3' end, 3' flank size = 10
                { "MUC16", GencodeFuncotation.VariantClassification.THREE_PRIME_FLANK, MUC16_CONTIG, MUC16_START - 5, MUC16_START - 5, "C", "T", 0, 10, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },

                /*
                 * SNP 3' reverse strand NEGATIVE test cases:
                 */

                // SNP 1 base past the 3' end, 3' flank size = 0
                { null, GencodeFuncotation.VariantClassification.IGR, MUC16_CONTIG, MUC16_START - 1, MUC16_START - 1, "G", "A", 0, 0, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // SNP 1 base past the 3' end, 3' flank size = 0 and 5' flank size = 10
                { null, GencodeFuncotation.VariantClassification.IGR, MUC16_CONTIG, MUC16_START - 1, MUC16_START - 1, "G", "A", 10, 0, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // SNP 10 bases past the 3' end, 3' flank size = 9
                { null, GencodeFuncotation.VariantClassification.IGR, MUC16_CONTIG, MUC16_START - 10, MUC16_START - 10, "A", "T", 0, 9, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // SNP on last base of the 3' end, 3' flank size = 10
                { "MUC16", GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, MUC16_CONTIG, MUC16_START, MUC16_START, "T", "A", 10, 10, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },

                /*
                 * Multi-base Deletion 5' reverse strand POSITIVE test cases:
                 */

                // 3-base deletion in the first 3 bases of the 5' flank region, 5' flank size = 10
                { "MUC16", GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK, MUC16_CONTIG, MUC16_END + 1, MUC16_END + 3, "CTA", "A", 10, 0, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // 3-base deletion in the middle of the 5' flank region, 5' flank size = 10
                { "MUC16", GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK, MUC16_CONTIG, MUC16_END + 4, MUC16_END + 6, "ATT", "A", 10, 0, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // 3-base deletion in the last 3 bases of the 5' flank region, 5' flank size = 10
                { "MUC16", GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK, MUC16_CONTIG, MUC16_END + 8, MUC16_END + 10, "AGA", "A", 10, 0, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // 3-base deletion, 2 bases in the 5' flank and 1 base outside of it, 5' flank size = 10
                { "MUC16", GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK, MUC16_CONTIG, MUC16_END + 9, MUC16_END + 11, "GAA", "G", 10, 0, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // 3-base deletion, 1 base in the 5' flank and 2 bases outside of it, 5' flank size = 10
                { "MUC16", GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK, MUC16_CONTIG, MUC16_END + 10, MUC16_END + 12, "AAG", "A", 10, 0, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },

                /*
                 * Multi-base Deletion 5' reverse strand NEGATIVE test cases:
                 */

                // 3-base deletion entirely in the 5' UTR, 5' flank size = 10
                { "MUC16", GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR, MUC16_CONTIG, MUC16_END - 3, MUC16_END - 1, "GCT", "T", 10, 0, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // 3-base deletion at the very end of the 5' UTR, flank size = 10
                { "MUC16", GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR, MUC16_CONTIG, MUC16_END - 2, MUC16_END, "CTT", "C", 10, 0, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // 3-base deletion with 2 bases in the 5' UTR and 1 base in the 5' flank, flank size = 10
                { "MUC16", GencodeFuncotation.VariantClassification.COULD_NOT_DETERMINE, MUC16_CONTIG, MUC16_END - 1, MUC16_END + 1, "TTC", "T", 10, 0, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // 3-base deletion with 1 base in the 5' UTR and 2 bases in the 5' flank, flank size = 10
                { "MUC16", GencodeFuncotation.VariantClassification.COULD_NOT_DETERMINE, MUC16_CONTIG, MUC16_END, MUC16_END + 2, "TCT", "T", 10, 0, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // 3-base deletion 1 base past the end of the 5' flank, flank size = 10
                { null, GencodeFuncotation.VariantClassification.IGR, MUC16_CONTIG, MUC16_END + 11, MUC16_END + 13, "AGA", "A", 10, 0, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // 3-base deletion 2 bases past the end of the 5' flank, flank size = 10
                { null, GencodeFuncotation.VariantClassification.IGR, MUC16_CONTIG, MUC16_END + 12, MUC16_END + 14, "GAC", "G", 10, 0, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },

                /*
                 * Multi-base Deletion 3' reverse strand POSITIVE test cases:
                 */

                // 3-base deletion in the first 3 bases of the 3' flank region, 3' flank size = 10
                { "MUC16", GencodeFuncotation.VariantClassification.THREE_PRIME_FLANK, MUC16_CONTIG, MUC16_START - 3, MUC16_START - 1, "GTG", "G", 0, 10, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // 3-base deletion in the middle of the 3' flank region, 3' flank size = 10
                { "MUC16", GencodeFuncotation.VariantClassification.THREE_PRIME_FLANK, MUC16_CONTIG, MUC16_START - 6, MUC16_START - 4, "ACT", "A", 0, 10, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // 3-base deletion in the last 3 bases of the 3' flank region, 3' flank size = 10
                { "MUC16", GencodeFuncotation.VariantClassification.THREE_PRIME_FLANK, MUC16_CONTIG, MUC16_START - 10, MUC16_START - 8, "AGA", "A", 0, 10, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // 3-base deletion, 2 bases in the 3' flank and 1 base outside of it, 3' flank size = 10
                { "MUC16", GencodeFuncotation.VariantClassification.THREE_PRIME_FLANK, MUC16_CONTIG, MUC16_START - 11, MUC16_START - 9, "CAG", "C", 0, 10, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // 3-base deletion, 1 base in the 3' flank and 2 bases outside of it, 3' flank size = 10
                { "MUC16", GencodeFuncotation.VariantClassification.THREE_PRIME_FLANK, MUC16_CONTIG, MUC16_START - 12, MUC16_START - 10, "TCA", "T", 0, 10, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },

                /*
                 * Multi-base Deletion 3' reverse strand NEGATIVE test cases:
                 */

                // 3-base deletion entirely in the 3' UTR, 3' flank size = 10
                { "MUC16", GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, MUC16_CONTIG, MUC16_START + 1, MUC16_START + 3, "CCG", "C", 0, 10, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // 3-base deletion at the very end of the 3' UTR, flank size = 10
                { "MUC16", GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, MUC16_CONTIG, MUC16_START, MUC16_START + 2, "TCC", "T", 0, 10, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // 3-base deletion with 2 bases in the 3' UTR and 1 base in the 3' flank, flank size = 10
                { "MUC16", GencodeFuncotation.VariantClassification.COULD_NOT_DETERMINE, MUC16_CONTIG, MUC16_START - 1, MUC16_START + 1, "GTC", "G", 0, 10, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // 3-base deletion with 1 base in the 3' UTR and 2 bases in the 5' flank, flank size = 10
                { "MUC16", GencodeFuncotation.VariantClassification.COULD_NOT_DETERMINE, MUC16_CONTIG, MUC16_START - 2, MUC16_START, "TGT", "T", 0, 10, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // 3-base deletion 1 base past the end of the 3' flank, flank size = 10
                { null, GencodeFuncotation.VariantClassification.IGR, MUC16_CONTIG, MUC16_START - 13, MUC16_START - 11, "CTC", "C", 0, 10, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },
                // 3-base deletion 2 bases past the end of the 3' flank, flank size = 10
                { null, GencodeFuncotation.VariantClassification.IGR, MUC16_CONTIG, MUC16_START - 14, MUC16_START - 12, "GCT", "G", 0, 10, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.MUC16_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch19 },

                /* ***********************************
                 * PIK3CA Test Cases (forward strand):
                 *************************************/

                /*
                 * SNP 5' forward strand POSITIVE test cases:
                 */

                // SNP 1 base past the 5' end, 5' flank size = 3000
                { "PIK3CA", GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK, PIK3CA_CONTIG, PIK3CA_START - 1, PIK3CA_START - 1, "C", "A", 3000, 0, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // SNP 10 bases past the 5' end, 5' flank size = 10
                { "PIK3CA", GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK, PIK3CA_CONTIG, PIK3CA_START - 10, PIK3CA_START - 10, "C", "A", 10, 0, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // SNP 10 bases past the 5' end, 5' flank size = 10 and 3' flank size = 9
                { "PIK3CA", GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK, PIK3CA_CONTIG, PIK3CA_START - 10, PIK3CA_START - 10, "C", "A", 10, 9, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // SNP 5 bases past the 5' end, 5' flank size = 10
                { "PIK3CA", GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK, PIK3CA_CONTIG, PIK3CA_START - 5, PIK3CA_START - 5, "T", "A", 10, 0, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },

                /*
                 * SNP 5' forward strand NEGATIVE test cases:
                 */

                // SNP 1 base past the 5' end, 5' flank size = 0
                { null, GencodeFuncotation.VariantClassification.IGR, PIK3CA_CONTIG, PIK3CA_START - 1, PIK3CA_START - 1, "C", "A", 0, 0, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // SNP 1 base past the 5' end, 5' flank size = 0 and 3' flank size = 10
                { null, GencodeFuncotation.VariantClassification.IGR, PIK3CA_CONTIG, PIK3CA_START - 1, PIK3CA_START - 1, "C", "A", 0, 10, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // SNP 10 bases past the 5' end, 5' flank size = 9
                { null, GencodeFuncotation.VariantClassification.IGR, PIK3CA_CONTIG, PIK3CA_START - 10, PIK3CA_START - 10, "C", "A", 9, 0, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // SNP on first base of the 5' end, 5' flank size = 10
                { "PIK3CA", GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR, PIK3CA_CONTIG, PIK3CA_START, PIK3CA_START, "T", "A", 10, 10, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },

                /*
                 * SNP 3' forward strand POSITIVE test cases:
                 */

                // SNP 1 base past the 3' end, 3' flank size = 100
                { "PIK3CA", GencodeFuncotation.VariantClassification.THREE_PRIME_FLANK, PIK3CA_CONTIG, PIK3CA_END + 1, PIK3CA_END + 1, "C", "A", 0, 100, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // SNP 10 bases past the 3' end, 3' flank size = 10
                { "PIK3CA", GencodeFuncotation.VariantClassification.THREE_PRIME_FLANK, PIK3CA_CONTIG, PIK3CA_END + 10, PIK3CA_END + 10, "T", "A", 0, 10, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // SNP 10 bases past the 3' end, 3' flank size = 10 and 5' flank size = 9
                { "PIK3CA", GencodeFuncotation.VariantClassification.THREE_PRIME_FLANK, PIK3CA_CONTIG, PIK3CA_END + 10, PIK3CA_END + 10, "T", "A", 9, 10, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // SNP 5 bases past the 3' end, 3' flank size = 10
                { "PIK3CA", GencodeFuncotation.VariantClassification.THREE_PRIME_FLANK, PIK3CA_CONTIG, PIK3CA_END + 5, PIK3CA_END + 5, "T", "A", 0, 10, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },

                /*
                 * SNP 3' forward strand NEGATIVE test cases:
                 */

                // SNP 1 base past the 3' end, 3' flank size = 0
                { null, GencodeFuncotation.VariantClassification.IGR, PIK3CA_CONTIG, PIK3CA_END + 1, PIK3CA_END + 1, "C", "A", 0, 0, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // SNP 1 base past the 3' end, 3' flank size = 0 and 5' flank size = 10
                { null, GencodeFuncotation.VariantClassification.IGR, PIK3CA_CONTIG, PIK3CA_END + 1, PIK3CA_END + 1, "C", "A", 10, 0, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // SNP 10 bases past the 3' end, 3' flank size = 9
                { null, GencodeFuncotation.VariantClassification.IGR, PIK3CA_CONTIG, PIK3CA_END + 10, PIK3CA_END + 10, "T", "A", 0, 9, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // SNP on last base of the 3' end, 3' flank size = 10
                { "PIK3CA", GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, PIK3CA_CONTIG, PIK3CA_END, PIK3CA_END, "A", "C", 10, 10, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },

                /*
                 * Multi-base Deletion 5' forward strand POSITIVE test cases:
                 */

                // 3-base deletion in the last 3 bases of the 5' flank region, 5' flank size = 10
                { "PIK3CA", GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK, PIK3CA_CONTIG, PIK3CA_START - 3, PIK3CA_START - 1, "TCC", "T", 10, 0, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // 3-base deletion in the middle of the 5' flank region, 5' flank size = 10
                { "PIK3CA", GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK, PIK3CA_CONTIG, PIK3CA_START - 6, PIK3CA_START - 4, "CTC", "C", 10, 0, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // 3-base deletion in the first 3 bases of the 5' flank region, 5' flank size = 10
                { "PIK3CA", GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK, PIK3CA_CONTIG, PIK3CA_START - 10, PIK3CA_START - 8, "CCC", "C", 10, 0, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // 3-base deletion, 2 bases in the 5' flank and 1 base outside of it, 5' flank size = 10
                { "PIK3CA", GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK, PIK3CA_CONTIG, PIK3CA_START - 11, PIK3CA_START - 9, "GCC", "G", 10, 0, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // 3-base deletion, 1 base in the 5' flank and 2 bases outside of it, 5' flank size = 10
                { "PIK3CA", GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK, PIK3CA_CONTIG, PIK3CA_START - 12, PIK3CA_START - 10, "CGC", "C", 10, 0, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },

                /*
                 * Multi-base Deletion 5' forward strand NEGATIVE test cases:
                 */

                // 3-base deletion entirely in the 5' UTR, 5' flank size = 10
                { "PIK3CA", GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR, PIK3CA_CONTIG, PIK3CA_START + 1, PIK3CA_START + 3, "CTC", "C", 10, 0, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // 3-base deletion at the very start of the 5' UTR, flank size = 10
                { "PIK3CA", GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR, PIK3CA_CONTIG, PIK3CA_START, PIK3CA_START + 2, "TCT", "T", 10, 0, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // 3-base deletion with 2 bases in the 5' UTR and 1 base in the 5' flank, flank size = 10
                { "PIK3CA", GencodeFuncotation.VariantClassification.COULD_NOT_DETERMINE, PIK3CA_CONTIG, PIK3CA_START - 1, PIK3CA_START + 1, "CTC", "C", 10, 0, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // 3-base deletion with 1 base in the 5' UTR and 2 bases in the 5' flank, flank size = 10
                { "PIK3CA", GencodeFuncotation.VariantClassification.COULD_NOT_DETERMINE, PIK3CA_CONTIG, PIK3CA_START - 2, PIK3CA_START, "CCT", "C", 10, 0, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // 3-base deletion 1 base past the end of the 5' flank, flank size = 10
                { null, GencodeFuncotation.VariantClassification.IGR, PIK3CA_CONTIG, PIK3CA_START - 13, PIK3CA_START - 11, "CCG", "C", 10, 0, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // 3-base deletion 2 bases past the end of the 5' flank, flank size = 10
                { null, GencodeFuncotation.VariantClassification.IGR, PIK3CA_CONTIG, PIK3CA_START - 14, PIK3CA_START - 12, "GCC", "G", 10, 0, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },

                /*
                 * Multi-base Deletion 3' forward strand POSITIVE test cases:
                 */

                // 3-base deletion in the first 3 bases of the 3' flank region, 3' flank size = 10
                { "PIK3CA", GencodeFuncotation.VariantClassification.THREE_PRIME_FLANK, PIK3CA_CONTIG, PIK3CA_END + 1, PIK3CA_END + 3, "CTT", "C", 0, 10, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // 3-base deletion in the middle of the 3' flank region, 3' flank size = 10
                { "PIK3CA", GencodeFuncotation.VariantClassification.THREE_PRIME_FLANK, PIK3CA_CONTIG, PIK3CA_END + 4, PIK3CA_END + 6, "TTT", "T", 0, 10, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // 3-base deletion in the last 3 bases of the 3' flank region, 3' flank size = 10
                { "PIK3CA", GencodeFuncotation.VariantClassification.THREE_PRIME_FLANK, PIK3CA_CONTIG, PIK3CA_END + 8, PIK3CA_END + 10, "TGT", "T", 0, 10, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // 3-base deletion, 2 bases in the 3' flank and 1 base outside of it, 3' flank size = 10
                { "PIK3CA", GencodeFuncotation.VariantClassification.THREE_PRIME_FLANK, PIK3CA_CONTIG, PIK3CA_END + 9, PIK3CA_END + 11, "GTC", "G", 0, 10, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // 3-base deletion, 1 base in the 3' flank and 2 bases outside of it, 3' flank size = 10
                { "PIK3CA", GencodeFuncotation.VariantClassification.THREE_PRIME_FLANK, PIK3CA_CONTIG, PIK3CA_END + 10, PIK3CA_END + 12, "TCC", "T", 0, 10, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },

                /*
                 * Multi-base Deletion 3' forward strand NEGATIVE test cases:
                 */

                // 3-base deletion entirely in the 3' UTR, 3' flank size = 10
                { "PIK3CA", GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, PIK3CA_CONTIG, PIK3CA_END - 3, PIK3CA_END - 1, "ATT", "A", 0, 10, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // 3-base deletion at the very end of the 3' UTR, flank size = 10
                { "PIK3CA", GencodeFuncotation.VariantClassification.THREE_PRIME_UTR, PIK3CA_CONTIG, PIK3CA_END - 2, PIK3CA_END, "TTA", "T", 0, 10, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // 3-base deletion with 2 bases in the 3' UTR and 1 base in the 3' flank, flank size = 10
                { "PIK3CA", GencodeFuncotation.VariantClassification.COULD_NOT_DETERMINE, PIK3CA_CONTIG, PIK3CA_END - 1, PIK3CA_END + 1, "TAC", "T", 0, 10, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // 3-base deletion with 1 base in the 3' UTR and 2 bases in the 3' flank, flank size = 10
                { "PIK3CA", GencodeFuncotation.VariantClassification.COULD_NOT_DETERMINE, PIK3CA_CONTIG, PIK3CA_END, PIK3CA_END + 2, "ACT", "A", 0, 10, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // 3-base deletion 1 base past the end of the 3' flank, flank size = 10
                { null, GencodeFuncotation.VariantClassification.IGR, PIK3CA_CONTIG, PIK3CA_END + 11, PIK3CA_END + 13, "CCA", "C", 0, 10, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },
                // 3-base deletion 2 bases past the end of the 3' flank, flank size = 10
                { null, GencodeFuncotation.VariantClassification.IGR, PIK3CA_CONTIG, PIK3CA_END + 12, PIK3CA_END + 14, "CAA", "C", 0, 10, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_ANNOTATIONS_FILE_NAME, FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE, refDataSourceHg19Ch3 },

        };
    }

    @Test(dataProvider = "provideDataForCreateFuncotationsWithFlanks")
    public void testCreateFuncotationsWithFlanks(final String expectedGeneName,
                                          final GencodeFuncotation.VariantClassification expectedVariantClassification,
                                          final String contig,
                                          final int start,
                                          final int end,
                                          final String ref,
                                          final String alt,
                                          final int fivePrimeFlankSize,
                                          final int threePrimeFlankSize,
                                          final String transcriptGtfFile,
                                          final String transcriptFastaFile,
                                          final ReferenceDataSource referenceDataSource) {

        final Allele refAllele = Allele.create(ref, true);
        final Allele altAllele = Allele.create(alt);

        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder(
                "",
                contig,
                start,
                end,
                Arrays.asList(refAllele, altAllele)
        );
        final VariantContext variantContext = variantContextBuilder.make();
        final SimpleInterval variantInterval = new SimpleInterval(contig, start, end);

        final ReferenceContext referenceContext = new ReferenceContext(referenceDataSource, variantInterval);

        final FeatureInput<GencodeGtfFeature> gencodeFeatureInput = new FeatureInput<>(transcriptGtfFile, GencodeFuncotationFactory.DEFAULT_NAME, Collections.emptyMap());
        final Map<FeatureInput<? extends Feature>, Class<? extends Feature>> featureInputMap = new HashMap<>();
        featureInputMap.put(gencodeFeatureInput, GencodeGtfFeature.class);
        final FeatureContext featureContext = FeatureContext.createFeatureContextForTesting(featureInputMap, "dummyName", variantInterval, VariantWalker.DEFAULT_DRIVING_VARIANTS_LOOKAHEAD_BASES, 0, 0, null);

        // Create a factory for our funcotations:
        try (final GencodeFuncotationFactory funcotationFactory = new GencodeFuncotationFactory(
                IOUtils.getPath(transcriptFastaFile),
                "VERSION",
                GencodeFuncotationFactory.DEFAULT_NAME,
                FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_DEFAULT_VALUE,
                Collections.emptySet(),
                new LinkedHashMap<>(),
                gencodeFeatureInput,
                new FlankSettings(fivePrimeFlankSize,threePrimeFlankSize),
                "TEST")) {

            // We test against createFuncotations() rather than createFuncotationsOnVariant() here because
            // the flanking feature relies on the query done in createFuncotations(), and we need to test
            // that the query is padded appropriately.
            final List<Funcotation> funcotations = funcotationFactory.createFuncotations(variantContext, referenceContext, featureContext);

            // Make sure we get what we expected:
            Assert.assertEquals(funcotations.size(), 1);

            final GencodeFuncotation funcotation = (GencodeFuncotation)funcotations.get(0);

            Assert.assertEquals(funcotation.getVariantClassification(), expectedVariantClassification, "Variant classification not correct");
            Assert.assertEquals(funcotation.getHugoSymbol(), expectedGeneName, "Gene name not correct");
        }
    }

    @Test
    public void testRequiresFeatures() {
        Assert.assertTrue(testMuc16SnpCreateFuncotationsFuncotationFactory.requiresFeatures());
    }

    @Test
    public void testFivePrimeFlankSorting() {
        // This variant context would get a sorting that was different than the ground truth in FuncotatorIntegrationTest#nonTrivialLargeDataValidationTest
        final VariantContext variantContext = new VariantContextBuilder("", "1", 871276, 871276,
                Arrays.asList(Allele.create("G", true), Allele.create("A"))).make();
        final String transcriptGtfFile = GENCODE_HG19_BIG_GTF_FILE;
        final String transcriptFastaFile = GENCODE_HG19_BIG_FASTA_FILE;
        final SimpleInterval variantInterval = new SimpleInterval(variantContext);
        final ReferenceContext referenceContext = new ReferenceContext(ReferenceDataSource.of(Paths.get(b37Reference)), variantInterval);
        final VariantContext vcHg19 = new VariantContextBuilder(variantContext).chr(FuncotatorUtils.convertB37ContigToHg19Contig(variantContext.getContig())).make();
        final SimpleInterval vcHg19Interval = new SimpleInterval(vcHg19.getContig(), vcHg19.getStart(), vcHg19.getEnd());

        final FeatureInput<GencodeGtfFeature> gencodeFeatureInput = new FeatureInput<>(transcriptGtfFile, GencodeFuncotationFactory.DEFAULT_NAME, Collections.emptyMap());
        final Map<FeatureInput<? extends Feature>, Class<? extends Feature>> featureInputMap = new HashMap<>();
        featureInputMap.put(gencodeFeatureInput, GencodeGtfFeature.class);
        final FeatureContext featureContext = FeatureContext.createFeatureContextForTesting(featureInputMap, "dummyName", vcHg19Interval, VariantWalker.DEFAULT_DRIVING_VARIANTS_LOOKAHEAD_BASES, 0, 0, null);

        // Sorts canonically
        try (final GencodeFuncotationFactory funcotationFactory = new GencodeFuncotationFactory(
                IOUtils.getPath(transcriptFastaFile),
                "VERSION",
                GencodeFuncotationFactory.DEFAULT_NAME,
                FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_DEFAULT_VALUE,
                Collections.emptySet(),
                new LinkedHashMap<>(),
                gencodeFeatureInput,
                new FlankSettings(5000, 0),
                "TEST")) {

            // We have to set the contig conversion because of the datasources / B37 mismatch:
            funcotationFactory.setDoExonContigConversionToB37ForTranscripts(true);

            final List<Funcotation> gencodeFuncotationList = funcotationFactory.createFuncotations(vcHg19, referenceContext, featureContext);

            System.out.println(gencodeFuncotationList.size());
        }
    }

//    @Test
//    public void testGetCodingSequenceFromTranscriptFasta() {
//
//
//        getCodingSequenceFromReferenceContext();
//
//        final String transcriptFastaFile    = FuncotatorTestConstants.PIK3CA_ALL_TRANSCRIPTS_GENCODE_TRANSCRIPT_FASTA_FILE;
//        final String transcriptId           = "ENST00000263967.1";
//        final String expectedCodingSequence = "ATGAGTGGGAGCTGCT";
//    }
}
