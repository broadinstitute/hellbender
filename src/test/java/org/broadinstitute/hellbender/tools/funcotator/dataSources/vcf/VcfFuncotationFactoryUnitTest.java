package org.broadinstitute.hellbender.tools.funcotator.dataSources.vcf;

import com.google.common.collect.ImmutableMap;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.lang3.RandomStringUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.lang3.tuple.Triple;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.testutils.FuncotatorReferenceTestUtils;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorArgumentDefinitions;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorTestConstants;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.FuncotatorTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Class to test the {@link VcfFuncotationFactory}.
 * Created by jonn on 3/26/18.
 */
public class VcfFuncotationFactoryUnitTest extends GATKBaseTest {

    //==================================================================================================================
    // Private Static Members:

    private static final String FACTORY_NAME = "TestFactory";
    private static final String FACTORY_VERSION = "TEST_VERSION";
    private static final String EXAC_SNIPPET = toolsTestDir + "funcotator/test_exac.vcf";

    private static final String DEFAULT_FILTER_STRING = "TODAY;A;Variant;Was;FILTERED";

    //==================================================================================================================
    // Private Members:
    private static final LinkedHashMap<String, Object> FIELD_DEFAULT_MAP = new LinkedHashMap<>();

    static {
        FIELD_DEFAULT_MAP.put("ASP", "false");
        FIELD_DEFAULT_MAP.put("ASS", "false");
        FIELD_DEFAULT_MAP.put("CAF", "");
        FIELD_DEFAULT_MAP.put("CDA", "false");
        FIELD_DEFAULT_MAP.put("CFL", "false");
        FIELD_DEFAULT_MAP.put("COMMON", "");
        FIELD_DEFAULT_MAP.put("DSS", "false");
        FIELD_DEFAULT_MAP.put("G5", "false");
        FIELD_DEFAULT_MAP.put("G5A", "false");
        FIELD_DEFAULT_MAP.put("GENEINFO", "");
        FIELD_DEFAULT_MAP.put("GNO", "false");
        FIELD_DEFAULT_MAP.put("HD", "false");
        FIELD_DEFAULT_MAP.put("INT", "false");
        FIELD_DEFAULT_MAP.put("KGPhase1", "false");
        FIELD_DEFAULT_MAP.put("KGPhase3", "false");
        FIELD_DEFAULT_MAP.put("LSD", "false");
        FIELD_DEFAULT_MAP.put("MTP", "false");
        FIELD_DEFAULT_MAP.put("MUT", "false");
        FIELD_DEFAULT_MAP.put("NOC", "false");
        FIELD_DEFAULT_MAP.put("NOV", "false");
        FIELD_DEFAULT_MAP.put("NSF", "false");
        FIELD_DEFAULT_MAP.put("NSM", "false");
        FIELD_DEFAULT_MAP.put("NSN", "false");
        FIELD_DEFAULT_MAP.put("OM", "false");
        FIELD_DEFAULT_MAP.put("OTH", "false");
        FIELD_DEFAULT_MAP.put("PM", "false");
        FIELD_DEFAULT_MAP.put("PMC", "false");
        FIELD_DEFAULT_MAP.put("R3", "false");
        FIELD_DEFAULT_MAP.put("R5", "false");
        FIELD_DEFAULT_MAP.put("REF", "false");
        FIELD_DEFAULT_MAP.put("RS", "");
        FIELD_DEFAULT_MAP.put("RSPOS", "");
        FIELD_DEFAULT_MAP.put("RV", "false");
        FIELD_DEFAULT_MAP.put("S3D", "false");
        FIELD_DEFAULT_MAP.put("SAO", "");
        FIELD_DEFAULT_MAP.put("SLO", "false");
        FIELD_DEFAULT_MAP.put("SSR", "");
        FIELD_DEFAULT_MAP.put("SYN", "false");
        FIELD_DEFAULT_MAP.put("TOPMED", "");
        FIELD_DEFAULT_MAP.put("TPA", "false");
        FIELD_DEFAULT_MAP.put("U3", "false");
        FIELD_DEFAULT_MAP.put("U5", "false");
        FIELD_DEFAULT_MAP.put("VC", "");
        FIELD_DEFAULT_MAP.put("VLD", "false");
        FIELD_DEFAULT_MAP.put("VP", "");
        FIELD_DEFAULT_MAP.put("WGT", "");
        FIELD_DEFAULT_MAP.put("WTD", "false");
        FIELD_DEFAULT_MAP.put("dbSNPBuildID", "");
        FIELD_DEFAULT_MAP.put("ID", "");
        FIELD_DEFAULT_MAP.put("FILTER", "");
    }

    //==================================================================================================================
    // Helper Methods:



    private Object[] helpProvideForTestCreateFuncotations(final String variantFeatureFileName,
                                                          final String contig,
                                                          final int start,
                                                          final int end,
                                                          final String refAlleleString,
                                                          final String altAlleleString,
                                                          final List<Funcotation> expected) {
        return new Object[]{
                variantFeatureFileName,
                FuncotatorTestUtils.createSimpleVariantContext(FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(), contig, start, end, refAlleleString, altAlleleString),
                new ReferenceContext(ReferenceDataSource.of(IOUtils.getPath(b37Reference)), new SimpleInterval(contig, start, end)),
                expected
        };
    }


    //==================================================================================================================
    // Data Providers:

    @DataProvider
    private Iterator<Object[]> provideForTestGetName() {
        final ArrayList<Object[]> data = new ArrayList<>();

        for (int i = 0; i < 20; ++i) {
            data.add(new Object[]{RandomStringUtils.randomAlphanumeric(i)});
        }

        return data.iterator();
    }

    @DataProvider
    private Object[][] provideForTestCreateFuncotationsOnVariant() {

        return new Object[][]{
                // Trivial Case: No overlapping features:
                helpProvideForTestCreateFuncotations(FuncotatorTestConstants.DBSNP_HG19_SNIPPET_FILE_PATH,
                        "3", 61650, 61650, "T", "C",
                        Collections.singletonList(
                                TableFuncotation.create(FIELD_DEFAULT_MAP.keySet().stream().map(s -> FACTORY_NAME + "_" + s).collect(Collectors.toList()),
                                        FIELD_DEFAULT_MAP.values().stream().map(Object::toString).collect(Collectors.toList()),
                                        Allele.create("C"), FACTORY_NAME, null)
                        )
                ),
                // One overlapping VCF feature:
                helpProvideForTestCreateFuncotations(FuncotatorTestConstants.DBSNP_HG19_SNIPPET_FILE_PATH,
                        "3", 61662, 61662, "T", "C",
                        Collections.singletonList(
                                TableFuncotation.create(FIELD_DEFAULT_MAP.keySet().stream().map(s -> FACTORY_NAME + "_" + s).collect(Collectors.toList()),
                                        Arrays.asList("true", "false", "0.9744,0.02556", "false", "false", "1", "false", "true", "false", "", "false", "true", "false", "true", "true", "false", "false", "false", "false", "false", "false", "false", "false", "false", "false", "false", "false", "false", "false", "false", "73009205", "61662", "false", "false", "0", "true", "0", "false", "0.954392,0.0456075", "false", "false", "false", "SNV", "true", "0x05010000000515043e000100", "1", "false", "130", "rs73009205", ""),
                                        Allele.create("C"), FACTORY_NAME, null)
                        )
                ),
                // No matching VCF features (three overlap by position only), since there are no indels in dbSNP (the test datasource), so the ground truth should be a default entry, which was constructed here manually:
                helpProvideForTestCreateFuncotations(FuncotatorTestConstants.DBSNP_HG19_SNIPPET_FILE_PATH,
                        "3", 64157, 64166, "AGAAAGGTCA", "TCTTTCCAGT",
                        Collections.singletonList(TableFuncotation.create(FIELD_DEFAULT_MAP.keySet().stream().map(s -> FACTORY_NAME + "_" + s).collect(Collectors.toList()),
                                Arrays.asList("false", "false", "", "false", "false", "", "false", "false", "false", "", "false", "false", "false", "false", "false", "false", "false", "false", "false", "false", "false", "false", "false", "false", "false", "false", "false", "false", "false", "false", "", "", "false", "false", "", "false", "", "false", "", "false", "false", "false", "", "false", "", "", "false", "", "", ""),
                                Allele.create("TCTTTCCAGT"), FACTORY_NAME, null))
                ),
                // One overlapping VCF feature, non-empty FILTER field:
                helpProvideForTestCreateFuncotations(FuncotatorTestConstants.DBSNP_HG19_SNIPPET_WITH_FILTERS_FILE_PATH,
                        "3", 61662, 61662, "T", "C",
                        Collections.singletonList(
                                TableFuncotation.create(FIELD_DEFAULT_MAP.keySet().stream().map(s -> FACTORY_NAME + "_" + s).collect(Collectors.toList()),
                                        Arrays.asList("true", "false", "0.9744,0.02556", "false", "false", "1", "false", "true", "false", "", "false", "true", "false", "true", "true", "false", "false", "false", "false", "false", "false", "false", "false", "false", "false", "false", "false", "false", "false", "false", "73009205", "61662", "false", "false", "0", "true", "0", "false", "0.954392,0.0456075", "false", "false", "false", "SNV", "true", "0x05010000000515043e000100", "1", "false", "130", "rs73009205", "FILTER_73"),
                                        Allele.create("C"), FACTORY_NAME, null)
                        )
                ),
        };
    }

    //==================================================================================================================
    // Tests:

    @Test
    public void testGetAnnotationFeatureClass() {
        final VcfFuncotationFactory vcfFuncotationFactory = createVcfFuncotationFactory(FACTORY_NAME, FACTORY_VERSION, IOUtils.getPath(FuncotatorTestConstants.VARIANT_FILE_HG19_CHR3));
        Assert.assertEquals(vcfFuncotationFactory.getAnnotationFeatureClass(), VariantContext.class);
    }

    @Test
    public void testRequiresFeatures() {
        Assert.assertTrue(createVcfFuncotationFactory(FACTORY_NAME, FACTORY_VERSION, IOUtils.getPath(FuncotatorTestConstants.VARIANT_FILE_HG19_CHR3)).requiresFeatures());
    }

    @Test
    public void testGetType() {
        final VcfFuncotationFactory vcfFuncotationFactory = createVcfFuncotationFactory(FACTORY_NAME, FACTORY_VERSION, IOUtils.getPath(FuncotatorTestConstants.VARIANT_FILE_HG19_CHR3));
        Assert.assertEquals(vcfFuncotationFactory.getType(), FuncotatorArgumentDefinitions.DataSourceType.VCF);
    }

    @Test(dataProvider = "provideForTestGetName")
    public void testGetName(final String name) {
        final VcfFuncotationFactory vcfFuncotationFactory = createVcfFuncotationFactory(name, FACTORY_VERSION, IOUtils.getPath(FuncotatorTestConstants.VARIANT_FILE_HG19_CHR3));
        Assert.assertEquals(vcfFuncotationFactory.getName(), name);
    }

    @Test
    public void testGetSupportedFuncotationFields() {
        final VcfFuncotationFactory vcfFuncotationFactory =
                createVcfFuncotationFactory(FACTORY_NAME, FACTORY_VERSION, IOUtils.getPath(FuncotatorTestConstants.DBSNP_HG19_SNIPPET_FILE_PATH));

        final LinkedHashSet<String> expectedFieldNames = new LinkedHashSet<>();

        for (final String field : FIELD_DEFAULT_MAP.keySet()) {
            expectedFieldNames.add(FACTORY_NAME + "_" + field);
        }

        Assert.assertEquals(vcfFuncotationFactory.getSupportedFuncotationFields(), expectedFieldNames);
    }

    @Test(dataProvider = "provideForTestCreateFuncotationsOnVariant")
    public void testCreateFuncotationsOnVariant(final String variantFeatureDataFileName,
                                                final VariantContext variant,
                                                final ReferenceContext referenceContext,
                                                final List<Funcotation> expected) {
        // Make our factory:
        final VcfFuncotationFactory vcfFuncotationFactory =
                createVcfFuncotationFactory(FACTORY_NAME, FACTORY_VERSION, IOUtils.getPath(variantFeatureDataFileName));

        // Create features from the file:
        final List<Feature> vcfFeatures;
        try (final VCFFileReader vcfReader = new VCFFileReader(IOUtils.getPath(variantFeatureDataFileName))) {
            vcfFeatures = vcfReader.query(variant.getContig(), variant.getStart(), variant.getEnd()).stream().collect(Collectors.toList());
        }

        // Check to see if we're in business:
        Assert.assertEquals(
                vcfFuncotationFactory.createFuncotationsOnVariant(
                        variant,
                        referenceContext,
                        vcfFeatures
                ),
                expected
        );

        Assert.assertEquals(
                vcfFuncotationFactory.createFuncotationsOnVariant(
                        variant,
                        referenceContext,
                        vcfFeatures,
                        Collections.emptyList()
                ),
                expected
        );
    }

    @Test(dataProvider = "provideForTestCreateFuncotationsOnVariant")
    public void testCreateFuncotationMetadata(final String variantFeatureDataFileName,
                                              final VariantContext variant,
                                              final ReferenceContext referenceContext,
                                              final List<Funcotation> expected) {
        // Don't need the expected gt for this test, but useful to reuse the data provider.
        // Make our factory:
        final VcfFuncotationFactory vcfFuncotationFactory =
                createVcfFuncotationFactory(FACTORY_NAME, FACTORY_VERSION, IOUtils.getPath(variantFeatureDataFileName));

        // Create features from the file:
        final List<Feature> vcfFeatures;
        try (final VCFFileReader vcfReader = new VCFFileReader(IOUtils.getPath(variantFeatureDataFileName))) {
            vcfFeatures = vcfReader.query(variant.getContig(), variant.getStart(), variant.getEnd()).stream().collect(Collectors.toList());
        }

        // test the metadata
        final List<Funcotation> funcotations = vcfFuncotationFactory.createFuncotationsOnVariant(
                variant,
                referenceContext,
                vcfFeatures,
                Collections.emptyList()
        );

        Assert.assertEquals(funcotations.stream().map(f -> f.getMetadata().retrieveAllHeaderInfo()).collect(Collectors.toSet()).size(), 1);
        final Pair<VCFHeader, List<VariantContext>> vcfInfo = VariantContextTestUtils.readEntireVCFIntoMemory(variantFeatureDataFileName);
        final List<VCFInfoHeaderLine> gtOutputVcfInfoHeaderLines = vcfFuncotationFactory.createFuncotationVcfInfoHeaderLines(vcfInfo.getLeft());

        // Get the info headers that are in the VCF and make sure that these are also present in the metadata
        final Set<String> headerInfoLines = funcotations.get(0).getFieldNames();
        final Set<String> metadataFields = funcotations.get(0).getMetadata().retrieveAllHeaderInfo().stream()
                .map(f -> f.getID())
                .collect(Collectors.toSet());
        Assert.assertEquals(metadataFields, headerInfoLines);
        Assert.assertEquals(metadataFields, vcfFuncotationFactory.getSupportedFuncotationFields());
        Assert.assertEquals(funcotations.get(0).getMetadata().retrieveAllHeaderInfo(), gtOutputVcfInfoHeaderLines);
    }

    @DataProvider
    public Object[][] provideMultiallelicTest() {
        // These were chosen to correspond to test cases in the test exac datasource VCF.
        return new Object[][] {
                // 3	69521	.	T	A,C AC_AMR=2,0; AC_Het=0,3,0 AC=2,3 -- DP_HIST=4891|699|176|41|7229|10522|4675|4512|4936|3378|1833|885|500|250|131|64|34|24|15|139,0|0|0|0|0|0|0|0|0|1|0|0|0|0|0|0|0|0|0|0,0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|3;
                //  Note that AC Het is of type=., so we should test that we return the entire string.
                {new SimpleInterval("3", 69521, 69521), Arrays.asList("T", "C"),
                        Collections.singletonList(ImmutableMap.of("_AC_AMR", "0", "_AC_Het", "0,3,0", "_AC", "3", "_DP_HIST", "4891|699|176|41|7229|10522|4675|4512|4936|3378|1833|885|500|250|131|64|34|24|15|139,0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|3"))},
                {new SimpleInterval("3", 69521, 69521), Arrays.asList("T", "A"),
                        Collections.singletonList(ImmutableMap.of("_AC_AMR", "2", "_AC_Het", "0,3,0", "_AC", "2","_DP_HIST", "4891|699|176|41|7229|10522|4675|4512|4936|3378|1833|885|500|250|131|64|34|24|15|139,0|0|0|0|0|0|0|0|0|1|0|0|0|0|0|0|0|0|0|0"))},
                // 3	69552	rs55874132	G	T,A,C  AC_AMR=3,0,0 AC_Het=1,1,0,0,0,0 AC=3,3,5  4764|1048|70|7|7472|10605|4702|4511|4937|3377|1835|886|500|250|128|63|35|22|13|117,0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1|0|0|1,0|0|0|0|0|0|0|0|0|0|0|1|0|0|0|0|0|0|0|1,3|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0
                {new SimpleInterval("3", 69552, 69552), Arrays.asList("G", "A"),
                        Collections.singletonList(ImmutableMap.of("_AC_AMR", "0", "_AC_Het", "1,1,0,0,0,0", "_AC", "3","_DP_HIST", "4764|1048|70|7|7472|10605|4702|4511|4937|3377|1835|886|500|250|128|63|35|22|13|117,0|0|0|0|0|0|0|0|0|0|0|1|0|0|0|0|0|0|0|1"))},
                {new SimpleInterval("3", 69552, 69552), Arrays.asList("G", "T"),
                        Collections.singletonList(ImmutableMap.of("_AC_AMR", "3", "_AC_Het", "1,1,0,0,0,0", "_AC", "3","_DP_HIST", "4764|1048|70|7|7472|10605|4702|4511|4937|3377|1835|886|500|250|128|63|35|22|13|117,0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1|0|0|1"))},
                {new SimpleInterval("3", 69552, 69552), Arrays.asList("G", "C"),
                        Collections.singletonList(ImmutableMap.of("_AC_AMR", "0", "_AC_Het", "1,1,0,0,0,0", "_AC", "5","_DP_HIST", "4764|1048|70|7|7472|10605|4702|4511|4937|3377|1835|886|500|250|128|63|35|22|13|117,3|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0"))},
                // 3	324682	.	ACCAGGCCCAGCTCATGCTTCTTTGCAGCCTCT	TCCAGGCCCAGCTCATGCTTCTTTGCAGCCTCT,A  AC=7,2; AC_AMR=0,0 ;AC_Het=1,0,0  DP_HIST=428|427|186|183|1953|705|127|19|1|2|1|0|0|0|0|0|0|0|0|0,0|0|1|0|1|1|0|1|0|0|0|0|0|0|0|0|0|0|0|0,0|1|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;
                {new SimpleInterval("3", 324682, 324714), Arrays.asList("ACCAGGCCCAGCTCATGCTTCTTTGCAGCCTCT", "A"),
                        Collections.singletonList(ImmutableMap.of("_AC_AMR", "0", "_AC_Het", "1,0,0", "_AC", "2","_DP_HIST", "428|427|186|183|1953|705|127|19|1|2|1|0|0|0|0|0|0|0|0|0,0|1|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0"))},
                {new SimpleInterval("3", 324682, 324714), Arrays.asList("ACCAGGCCCAGCTCATGCTTCTTTGCAGCCTCT", "TCCAGGCCCAGCTCATGCTTCTTTGCAGCCTCT", "A"),
                        Arrays.asList(ImmutableMap.of("_AC_AMR", "0", "_AC_Het", "1,0,0", "_AC", "7","_DP_HIST", "428|427|186|183|1953|705|127|19|1|2|1|0|0|0|0|0|0|0|0|0,0|0|1|0|1|1|0|1|0|0|0|0|0|0|0|0|0|0|0|0"),
                                ImmutableMap.of("_AC_AMR", "0", "_AC_Het", "1,0,0", "_AC", "2","_DP_HIST", "428|427|186|183|1953|705|127|19|1|2|1|0|0|0|0|0|0|0|0|0,0|1|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0"))},
                //HARD!!  Same as the previous test
                {new SimpleInterval("3", 324682, 324682), Arrays.asList("A", "T"),
                        Collections.singletonList(ImmutableMap.of("_AC_AMR", "0", "_AC_Het", "1,0,0", "_AC", "7","_DP_HIST", "428|427|186|183|1953|705|127|19|1|2|1|0|0|0|0|0|0|0|0|0,0|0|1|0|1|1|0|1|0|0|0|0|0|0|0|0|0|0|0|0"))},

                // Control case (no multiallelics)
                //3	13372	.	G	C AC=3; AC_AMR=0 AC_Het=0 DP_HIST=14728|2455|2120|518|121|499|534|314|111|21|10|2|2|0|0|0|0|0|0|0,1|0|0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0
                {new SimpleInterval("3", 13372, 13372), Arrays.asList("G", "C"),
                        Collections.singletonList(ImmutableMap.of("_AC_AMR", "0", "_AC_Het", "0", "_AC", "3","_DP_HIST", "14728|2455|2120|518|121|499|534|314|111|21|10|2|2|0|0|0|0|0|0|0,1|0|0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0"))},

                // Control case (no multiallelics in datasource, but multiallelic query) -- Should produce a funcotation for both alleles, but second one should be blank
                //3	13372	.	G	C AC=3; AC_AMR=0 AC_Het=0 DP_HIST=14728|2455|2120|518|121|499|534|314|111|21|10|2|2|0|0|0|0|0|0|0,1|0|0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0
                {new SimpleInterval("3", 13372, 13372), Arrays.asList("G", "C", "T"),
                        Arrays.asList(ImmutableMap.of("_AC_AMR", "0", "_AC_Het", "0", "_AC", "3","_DP_HIST", "14728|2455|2120|518|121|499|534|314|111|21|10|2|2|0|0|0|0|0|0|0,1|0|0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0"),
                                ImmutableMap.of("_AC_AMR", "", "_AC_Het", "", "_AC", "","_DP_HIST", ""))},
                // Control case (no multiallelics in datasource, but multiallelic query) -- Should produce a funcotation for both alleles, but first one should be blank
                //3	13372	.	G	C AC=3; AC_AMR=0 AC_Het=0 DP_HIST=14728|2455|2120|518|121|499|534|314|111|21|10|2|2|0|0|0|0|0|0|0,1|0|0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0
                {new SimpleInterval("3", 13372, 13372), Arrays.asList("G", "T", "C"),
                        Arrays.asList(ImmutableMap.of("_AC_AMR", "", "_AC_Het", "", "_AC", "","_DP_HIST", ""),
                                ImmutableMap.of("_AC_AMR", "0", "_AC_Het", "0", "_AC", "3","_DP_HIST", "14728|2455|2120|518|121|499|534|314|111|21|10|2|2|0|0|0|0|0|0|0,1|0|0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0"))},
        };
    }

    /**
     * This also tests the VCFFuncotationFactory Caching.
     */
    @Test(dataProvider = "provideMultiallelicTest")
    public void testQueryIntoMultiallelic(final SimpleInterval variantInterval, final List<String> alleles,
                                          final List<Map<String, String>> gtAttributes) {

        // Make the factory
        final VcfFuncotationFactory vcfFuncotationFactory =
                createVcfFuncotationFactory(FACTORY_NAME, FACTORY_VERSION, IOUtils.getPath(EXAC_SNIPPET));

        final ReferenceContext referenceContext = new ReferenceContext(ReferenceDataSource.of(Paths.get(FuncotatorReferenceTestUtils.retrieveB37Chr3Ref())), variantInterval);

        final List<Feature> vcfFeatures;
        try (final VCFFileReader vcfReader = new VCFFileReader(IOUtils.getPath(EXAC_SNIPPET))) {
            vcfFeatures = vcfReader.query(variantInterval.getContig(), variantInterval.getStart(), variantInterval.getEnd()).stream().collect(Collectors.toList());
        }

        final VariantContext variant = new VariantContextBuilder()
                .chr(variantInterval.getContig()).start(variantInterval.getStart()).stop(variantInterval.getEnd())
                .alleles(alleles)
                .make();

        final List<Funcotation> funcotations = vcfFuncotationFactory.createFuncotationsOnVariant(
                variant,
                referenceContext,
                vcfFeatures,
                Collections.emptyList()
        );
        Assert.assertEquals(funcotations.size(), alleles.size() - 1);
        Assert.assertEquals(funcotations.size(), gtAttributes.size());

        IntStream.range(0, funcotations.size()).forEach( j ->
                gtAttributes.get(j).keySet().forEach(k ->
                        Assert.assertEquals(funcotations.get(j).getField(vcfFuncotationFactory.getName() + k), gtAttributes.get(j).get(k), "Mismatch with " + k)
            )
        );

        Assert.assertEquals(vcfFuncotationFactory.cacheHits, 0);
        Assert.assertEquals(vcfFuncotationFactory.cacheMisses, 1);  // Should match the number of times createFuncotationOnVariant was called.

        // Funcotate again, so that we should get a cache hit.
        final List<Funcotation> funcotations2 = vcfFuncotationFactory.createFuncotationsOnVariant(
                variant,
                referenceContext,
                vcfFeatures,
                Collections.emptyList()
        );

        Assert.assertEquals(vcfFuncotationFactory.cacheHits, 1);
        Assert.assertEquals(vcfFuncotationFactory.cacheMisses, 1);
        Assert.assertEquals(funcotations, funcotations2, "Even though there was a cache hit, the funcotations are not equal.");

        // Sanity check that we get the same funcotations whether we use the three parameter or the four parameter
        //  version of createFuncotationsOnVariant.
        Assert.assertEquals(vcfFuncotationFactory.createFuncotationsOnVariant(
                variant,
                referenceContext,
                vcfFeatures), funcotations);
    }

    @Test
    public void testCacheOnObjectReference(){
        // This code is a bit complex, since the cache is based exclusively on object references.  That works great in Funcotator,
        //  but not as great in the general case (incl. autotests)
        // We do not care so much about the content of each variant context.  We change the position to control whether
        //  there is a cache hit or not.
        // Please note that this test does not actually test the content of the funcotations.  Just whether the cache
        //  was set to the appropriate size and that the hit/miss counters are being maintained properly.

        // Create dummy data.  Remember that since the cache is based on reference, we always have to index into this list.
        final List<String> alleles = Arrays.asList("G", "C", "T");
        final List<Triple<VariantContext, ReferenceContext, List<Feature>>> dummyTriples = IntStream.range(0, VcfFuncotationFactory.LRUCache.MAX_ENTRIES + 1)
                .boxed().map(i -> createDummyCacheTriples(alleles, i)).collect(Collectors.toList());

        // Create our funcotation factory to test
        final VcfFuncotationFactory vcfFuncotationFactory =
                createVcfFuncotationFactory(FACTORY_NAME, FACTORY_VERSION, IOUtils.getPath(EXAC_SNIPPET));

        for (int i = 0; i < VcfFuncotationFactory.LRUCache.MAX_ENTRIES; i++) {
            funcotateForCacheTest(vcfFuncotationFactory, dummyTriples.get(i));
            Assert.assertEquals(vcfFuncotationFactory.cacheHits, 0);
            Assert.assertEquals(vcfFuncotationFactory.cacheMisses, i+1);  // Should match the number of times createFuncotationOnVariant was called.
        }
        // We will get one more miss in this loop, since [0] will have been purged from the cache.  We will test this below.
        for (int i = 0; i < (VcfFuncotationFactory.LRUCache.MAX_ENTRIES + 1); i++) {
            funcotateForCacheTest(vcfFuncotationFactory, dummyTriples.get(i));
        }
        Assert.assertEquals(vcfFuncotationFactory.cacheHits, VcfFuncotationFactory.LRUCache.MAX_ENTRIES);

        // This should be another miss, since the variant at index = 0 should no longer be in the cache.
        funcotateForCacheTest(vcfFuncotationFactory, dummyTriples.get(0));
        Assert.assertEquals(vcfFuncotationFactory.cacheMisses, (VcfFuncotationFactory.LRUCache.MAX_ENTRIES + 2));
    }

    // Performs a dummy funcotation with an offset for controlling the cache.
    private void funcotateForCacheTest(final VcfFuncotationFactory vcfFuncotationFactory, final Triple<VariantContext, ReferenceContext, List<Feature>> cacheTriple) {
        vcfFuncotationFactory.createFuncotationsOnVariant(
                cacheTriple.getLeft(),
                cacheTriple.getMiddle(),
                cacheTriple.getRight(),
                Collections.emptyList()
        );
    }

    // Creates dummy triples for testing the cache.  If offset is zero, and ref="G" and alt(s) are in {"C", "T"}, there
    //  will be a hit on the exac snippet.  This is to make sure that one of the dummy triplets was a hit in a VCF
    //  funcotation factory.
    private Triple<VariantContext, ReferenceContext, List<Feature>> createDummyCacheTriples(final List<String> alleles, final int offset) {

        // Create an interval for a variant that overlaps an entry in the exac snippet test file when offset = 0.
        final SimpleInterval variantInterval = new SimpleInterval("3", 13372+offset, 13372+offset);
        final VariantContext vc =  new VariantContextBuilder()
                .chr(variantInterval.getContig()).start(variantInterval.getStart()).stop(variantInterval.getEnd())
                .alleles(alleles)
                .make();
        final ReferenceContext referenceContext = new ReferenceContext(ReferenceDataSource.of(Paths.get(FuncotatorReferenceTestUtils.retrieveB37Chr3Ref())), variantInterval);
        final List<Feature> vcfFeatures;
        try (final VCFFileReader vcfReader = new VCFFileReader(IOUtils.getPath(EXAC_SNIPPET))) {
            vcfFeatures = vcfReader.query(variantInterval.getContig(), variantInterval.getStart(), variantInterval.getEnd()).stream().collect(Collectors.toList());
        }

        return Triple.of(vc, referenceContext, vcfFeatures);
    }

    private VcfFuncotationFactory createVcfFuncotationFactory(final String name,
                                                              final String version,
                                                              final Path sourceFilePath) {
        return new VcfFuncotationFactory(name, version, sourceFilePath, new LinkedHashMap<>(), new FeatureInput<VariantContext>(sourceFilePath.toString(), name, new HashMap<>()));
    }

    @Test
    public void testNoSupportOfSegments() {
        final VcfFuncotationFactory factory = createVcfFuncotationFactory("FAKE", "MOCK1", IOUtils.getPath(EXAC_SNIPPET));

        Assert.assertFalse(factory.isSupportingSegmentFuncotation());
        Assert.assertEquals(factory.getSupportedFuncotationFieldsForSegments(), Collections.emptyList());
    }
}
