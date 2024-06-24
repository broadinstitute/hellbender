package org.broadinstitute.hellbender.tools.funcotator.mafOutput;

import com.google.common.collect.ImmutableMap;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.tools.funcotator.DataSourceFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorConstants;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorTestConstants;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.vcf.VcfFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadataUtils;
import org.broadinstitute.hellbender.tools.funcotator.metadata.TumorNormalPair;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.FuncotatorTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.funcotator.mafOutput.CustomMafFuncotationCreator.MAF_DBSNP_VAL_STATUS_DELIMITER;
import static org.broadinstitute.hellbender.tools.funcotator.mafOutput.CustomMafFuncotationCreator.MAF_DBSNP_VAL_STATUS_FIELD;
import static org.broadinstitute.hellbender.tools.funcotator.mafOutput.MafOutputRendererConstants.DBSNP_DS_NAME;

public class CustomMafFuncotationCreatorUnitTest extends GATKBaseTest {

    private static final String HG19_CHR1_1M_FASTA = publicTestDir + "Homo_sapiens_assembly19_chr1_1M.fasta";

    /**
     * The ground truth field is ordered by {@link CustomMafFuncotationCreator#COUNT_FIELD_NAMES}
     * @return array of objects to be used in tests.  Never {@code null}
     */
    @DataProvider
    public Object[][] provideVariantsAndPairs() {
        return new Object[][] {
                {new VariantContextBuilder("foo", "3", 1000000, 1000000, Arrays.asList(Allele.create("A", true), Allele.create("T", false)))
                        .genotypes(
                                Arrays.asList(
                                        new GenotypeBuilder("T1",
                                                Arrays.asList(Allele.create("A", true), Allele.create("T", false)))
                                                .AD(new int[]{10,2})
                                                .attribute("AF", "0.166").make(),
                                        new GenotypeBuilder("N1",
                                                Arrays.asList(Allele.create("A", true), Allele.create("T", false)))
                                                .AD(new int[]{20,0}).make()
                                        )
                        )
                        .make(),
                        Collections.singletonList(new TumorNormalPair("T1", "N1")),
                        ImmutableMap.of(Allele.create("T", false), Arrays.asList("A", "T", "2", "10", "0", "20", "0.166"))
                },{new VariantContextBuilder("foo", "3", 1000000, 1000000, Arrays.asList(Allele.create("A", true), Allele.create("T", false)))
                    .genotypes(
                            Arrays.asList(
                                    new GenotypeBuilder("T1",
                                            Arrays.asList(Allele.create("A", true), Allele.create("T", false)))
                                            .AD(new int[]{10,2})
                                            .attribute("AF", "0.166").make()
                                    )
                            )
                    .make(),
                    Collections.singletonList(new TumorNormalPair("T1", "")),
                    ImmutableMap.of(Allele.create("T", false), Arrays.asList(FuncotatorConstants.UNKNOWN_VALUE_STRING, FuncotatorConstants.UNKNOWN_VALUE_STRING, "2", "10", "", "", "0.166"))
                },
                {new VariantContextBuilder("foo", "3", 1000000, 1000000, Arrays.asList(Allele.create("A", true), Allele.create("T", false), Allele.create("AT", false)))
                        .genotypes(
                                Arrays.asList(
                                        new GenotypeBuilder("T1",
                                                Arrays.asList(Allele.create("A", true), Allele.create("T", false), Allele.create("AT", false)))
                                                .AD(new int[]{10,2,4})
                                                .attribute("AF", "0.166,0.32").make(),
                                        new GenotypeBuilder("N1",
                                                Arrays.asList(Allele.create("A", true), Allele.create("T", false), Allele.create("AT", false)))
                                                .AD(new int[]{20,0,1}).make()
                                )
                        )
                        .make(),
                        Collections.singletonList(new TumorNormalPair("T1", "N1")),
                        ImmutableMap.of(Allele.create("T", false), Arrays.asList("A", "T", "2", "10", "0", "20", "0.166"),
                                Allele.create("AT", false), Arrays.asList("A", "T", "4", "10", "1", "20", "0.32"))
                },

                // Test when we have no recognizable pairs that the fields are generated but empty.
                {new VariantContextBuilder("foo", "3", 1000000, 1000000, Arrays.asList(Allele.create("A", true), Allele.create("T", false), Allele.create("AT", false)))
                        .genotypes(
                                Arrays.asList(
                                        new GenotypeBuilder("GERMLINE1",
                                                Arrays.asList(Allele.create("A", true), Allele.create("T", false), Allele.create("AT", false)))
                                                .AD(new int[]{10,2,4})
                                                .attribute("AF", "0.166,0.32").make(),
                                        new GenotypeBuilder("GERMLINE2",
                                                Arrays.asList(Allele.create("A", true), Allele.create("T", false), Allele.create("AT", false)))
                                                .AD(new int[]{20,0,1}).make()
                                )
                        )
                        .make(),
                        Collections.emptyList(),
                        ImmutableMap.of(Allele.create("T", false), Arrays.asList(FuncotatorConstants.UNKNOWN_VALUE_STRING, FuncotatorConstants.UNKNOWN_VALUE_STRING, "", "", "", "",""),
                                Allele.create("AT", false),  Arrays.asList(FuncotatorConstants.UNKNOWN_VALUE_STRING, FuncotatorConstants.UNKNOWN_VALUE_STRING, "", "", "", "",""))
                }

        };
    }

    @Test(dataProvider = "provideVariantsAndPairs")
    public void testCreationOfCustomCountFields(final VariantContext variant, final List<TumorNormalPair> tnPair, final Map<Allele, List<String>> gtOrderedCustomFieldsByAllele) {
        final List<Funcotation> guess = CustomMafFuncotationCreator.createCustomMafCountFields(variant, tnPair);
        Assert.assertEquals(guess.size(), variant.getAlternateAlleles().size());
        for (final Allele allele: variant.getAlternateAlleles()) {
            final Funcotation funcotationOfInterest = guess.stream().filter(f -> f.getAltAllele().equals(allele)).findFirst().orElse(null);
            Assert.assertNotNull(funcotationOfInterest);
            final List<String> gtValues = gtOrderedCustomFieldsByAllele.getOrDefault(allele, Collections.emptyList());

            for (int i = 0; i < CustomMafFuncotationCreator.COUNT_FIELD_NAMES.size(); i++) {
                final String countField = CustomMafFuncotationCreator.COUNT_FIELD_NAMES.get(i);
                Assert.assertEquals(funcotationOfInterest.getField(countField), gtValues.get(i));
            }
        }
    }

    @DataProvider
    public Object[][] provideDbSnpVariants(){
        return new Object[][] {
                {new VariantContextBuilder("testing", "1", 10177, 10177,
                        Arrays.asList(Allele.create("A", true), Allele.create("AC"))).make(), 1, "byFrequency" + MAF_DBSNP_VAL_STATUS_DELIMITER + "by1000genomes"},
                {new VariantContextBuilder("testing", "1", 10352, 10352,
                        Arrays.asList(Allele.create("T", true), Allele.create("TA"))).make(), 2, "byFrequency" + MAF_DBSNP_VAL_STATUS_DELIMITER + "by1000genomes"},
        };
    }

    /**
     * As a side effect, this also tests (https://github.com/broadinstitute/gatk/issues/4972)
     */
    @Test(dataProvider = "provideDbSnpVariants")
    public void testCreateDbSnpCustomFields(final VariantContext variant, final int gtNumHits, final String gtDbSnpValStatusField) {

        final Path sourceFilePath = IOUtils.getPath(FuncotatorTestConstants.DBSNP_HG19_SNIPPET_FILE_PATH);
        final DataSourceFuncotationFactory vcfFuncotationFactory =
                new VcfFuncotationFactory(DBSNP_DS_NAME,
                        "snippetTest",
                        sourceFilePath,
                        new LinkedHashMap<>(),
                        new FeatureInput<VariantContext>(sourceFilePath.toString(), DBSNP_DS_NAME, new HashMap<>())
                );

        /* dbSNP records of relevance.
        1	10177	rs367896724	A	AC	.	.	RS=367896724;RSPOS=10177;dbSNPBuildID=138;SSR=0;SAO=0;VP=0x050000020005170026000200;GENEINFO=DDX11L1:100287102;WGT=1;VC=DIV;R5;ASP;VLD;G5A;G5;KGPhase3;CAF=0.5747,0.4253;COMMON=1
        1	10352	rs555500075	T	TA	.	.	RS=555500075;RSPOS=10352;dbSNPBuildID=142;SSR=0;SAO=0;VP=0x050000020005170026000200;GENEINFO=DDX11L1:100287102;WGT=1;VC=DIV;R5;ASP;VLD;G5A;G5;KGPhase3;CAF=0.5625,0.4375;COMMON=1
        1	10352	rs145072688	T	TA	.	.	RS=145072688;RSPOS=10353;dbSNPBuildID=134;SSR=0;SAO=0;VP=0x050000020005000002000200;GENEINFO=DDX11L1:100287102;WGT=1;VC=DIV;R5;ASP;CAF=0.5625,0.4375;COMMON=1
         */
        // Create features from the file:
        final List<Feature> vcfFeatures;
        try (final VCFFileReader vcfReader = new VCFFileReader(IOUtils.getPath(FuncotatorTestConstants.DBSNP_HG19_SNIPPET_FILE_PATH))) {
            vcfFeatures = vcfReader.query(variant.getContig(), variant.getStart(), variant.getEnd()).stream().collect(Collectors.toList());
        }

        Assert.assertEquals(vcfFeatures.size(), gtNumHits);

        final Map<String, List<Feature>> vcfFuncotationSourceMap = ImmutableMap.of(DBSNP_DS_NAME, vcfFeatures);
        final ReferenceContext referenceContext = new ReferenceContext(ReferenceDataSource.of(Paths.get(HG19_CHR1_1M_FASTA)),
                new SimpleInterval(variant.getContig(),
                variant.getStart(), variant.getEnd()));

        final FeatureContext featureContext = FuncotatorTestUtils.createFeatureContext(Collections.singletonList(vcfFuncotationFactory),
                "TEST_CREATE_DB_SNP_CUSTOM_FIELDS",
                new SimpleInterval(variant.getContig(), variant.getStart(), variant.getEnd()),
                0, 0, 0, null);

        final List<Funcotation> funcotations = vcfFuncotationFactory.createFuncotations(variant, referenceContext, featureContext);
        Assert.assertTrue(funcotations.size() > 0);
        for (final Funcotation f : funcotations) {
            Assert.assertEquals(Utils.split(f.getField(DBSNP_DS_NAME + "_VLD"), "|").size(), vcfFuncotationSourceMap.get(DBSNP_DS_NAME).size());
        }

        final List<Funcotation> customDbSnpFuncotations = CustomMafFuncotationCreator.createCustomMafDbSnpFields(funcotations);
        Assert.assertEquals(customDbSnpFuncotations.stream().map(f -> f.getField(MAF_DBSNP_VAL_STATUS_FIELD)).collect(Collectors.toList()),
                Collections.singletonList(gtDbSnpValStatusField));

        // Now add some dummy (non-DbSNP) funcotations and make sure that we are not getting any additional custom dbsnp maf fields.
        final List<String> dummyFieldNames = Arrays.asList("foo_field", "bar_field");
        final Funcotation dummyFuncotation = TableFuncotation.create(dummyFieldNames, Arrays.asList("1", "2"), Allele.create("AA"), "DUMMY", FuncotationMetadataUtils.createWithUnknownAttributes(dummyFieldNames));
        funcotations.add(dummyFuncotation);

        final List<Funcotation> customDbSnpFuncotationsWithoutDummies = CustomMafFuncotationCreator.createCustomMafDbSnpFields(funcotations);
        Assert.assertEquals(customDbSnpFuncotationsWithoutDummies, customDbSnpFuncotations);
    }

    @Test
    public void testDummyOnlyCustomDbSNPField() {
        final List<String> dummyFieldNames = Arrays.asList("foo_field", "bar_field");
        final Funcotation dummyFuncotation = TableFuncotation.create(dummyFieldNames, Arrays.asList("1", "2"), Allele.create("AA"), "DUMMY", FuncotationMetadataUtils.createWithUnknownAttributes(dummyFieldNames));
        final Funcotation dummyFuncotation2 = TableFuncotation.create(dummyFieldNames, Arrays.asList("4", "5"), Allele.create("AA"), "DUMMY", FuncotationMetadataUtils.createWithUnknownAttributes(dummyFieldNames));
        final List<Funcotation> funcotations = Arrays.asList(dummyFuncotation, dummyFuncotation2);

        // Should be empty
        Assert.assertEquals(CustomMafFuncotationCreator.createCustomMafDbSnpFields(funcotations).size(), 0);
    }
}
