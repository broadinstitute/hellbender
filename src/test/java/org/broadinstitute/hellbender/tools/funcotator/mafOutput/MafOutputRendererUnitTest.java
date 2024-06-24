package org.broadinstitute.hellbender.tools.funcotator.mafOutput;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.collections.MapUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.DummyPlaceholderGatkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.testutils.FuncotatorReferenceTestUtils;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.tools.funcotator.*;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.FuncotatorTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.*;

/**
 * Unit test class for the {@link MafOutputRenderer}.
 * Created by jonn on 1/22/18.
 */
public class MafOutputRendererUnitTest extends GATKBaseTest {

    //==================================================================================================================
    // Private Static Members:

    private final String OTHER_TRANSCRIPTS_RAW_NAMES =
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString()            + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString()            + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.IN_FRAME_INS.toString()            + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS.toString()         + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL.toString()         + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.MISSENSE.toString()                + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.NONSENSE.toString()                + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.SILENT.toString()                  + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.SPLICE_SITE.toString()             + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME.toString()  + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME.toString() + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.START_CODON_SNP.toString()         + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.START_CODON_INS.toString()         + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.START_CODON_DEL.toString()         + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.NONSTOP.toString()                 + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR.toString()          + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.THREE_PRIME_UTR.toString()         + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK.toString()        + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.INTRON.toString()                  + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.LINCRNA.toString();

    private final String OTHER_TRANSCRIPTS_MAF_READY_NAMES =
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString())     + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString())     + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.IN_FRAME_INS.toString())     + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS.toString())  + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL.toString())  + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.MISSENSE.toString())         + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.NONSENSE.toString())         + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.SILENT.toString())           + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.SPLICE_SITE.toString())      + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME.toString()                                                    + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME.toString()                                                   + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.START_CODON_SNP.toString()                                                           + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + GencodeFuncotation.VariantClassification.START_CODON_INS.toString()                                                           + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.START_CODON_DEL.toString())  + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.NONSTOP.toString())          + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR.toString())   + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.THREE_PRIME_UTR.toString())  + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK.toString()) + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.INTRON.toString())           + VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER +
            "MITD1_ENST00000466880.1_" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.LINCRNA.toString());

    //==================================================================================================================
    // Private Members:

    //==================================================================================================================
    // Helper Methods:

    private List<String> createFieldValuesFromNameList(final String prefix,
                                                       final List<String> baseFieldList,
                                                       final int fieldSize) {
        final List<String> outList = new ArrayList<>(baseFieldList.size());

        for ( int i = 0; i < baseFieldList.size() ; ++i ) {
            final String formatString = "%s%0" +
                    ((fieldSize - prefix.length()) > 0 ? fieldSize - prefix.length() : "") +
                    "d";
            outList.add(String.format(formatString, prefix, i+1));
        }

        return outList;
    }

    private MafOutputRenderer createMafOutputRenderer(final File outputFile, final String referenceVersion) {
        return createMafOutputRenderer(outputFile, referenceVersion, new HashSet<>());
    }

    private MafOutputRenderer createMafOutputRenderer(final File outputFile, final String referenceVersion, final Set<String> excludedFields) {

        final Map<Path, Properties> configData =
                DataSourceUtils.getAndValidateDataSourcesFromPaths(
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        Collections.singletonList(FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER)
                );
        
        final List<DataSourceFuncotationFactory> dataSourceFuncotationFactories = DataSourceUtils.createDataSourceFuncotationFactoriesForDataSources(
                configData,
                new LinkedHashMap<>(),
                TranscriptSelectionMode.BEST_EFFECT,
                new HashSet<>(),
                new DummyPlaceholderGatkTool(),
                FuncotatorArgumentDefinitions.LOOKAHEAD_CACHE_IN_BP_DEFAULT_VALUE,
                new FlankSettings(0, 0),
                false,
                FuncotatorUtils.DEFAULT_MIN_NUM_BASES_FOR_VALID_SEGMENT,
                FuncotatorUtils.DEFAULT_SPLICE_SITE_WINDOW_SIZE
        );

        // Sort the datasources to ensure the same order every time:
        dataSourceFuncotationFactories.sort( Comparator.comparing(DataSourceFuncotationFactory::getName) );

        return new MafOutputRenderer(
                IOUtils.getPath(outputFile.toURI().toString()),
                dataSourceFuncotationFactories,
                new VCFHeader(),
                new LinkedHashMap<>(),
                new LinkedHashMap<>(),
                new HashSet<>(),
                referenceVersion,
                excludedFields,
                "Unknown"
        );
    }

    private MafOutputRenderer createMafOutputRenderer(final String referenceVersion) {
        return createMafOutputRenderer(getSafeNonExistentFile("TestMafOutputFile"), referenceVersion);
    }

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    private Object[][] provideForTestReplaceFuncotationValuesWithMafCompliantValues() {
        return new Object[][] {
                // Empty maps:
                {
                    new HashMap<>(),
                    new LinkedHashMap<>()
                },
                // Singleton map that doesn't contain a replaceable element:
                {
                        MapUtils.putAll(new HashMap<String, Object>(),
                                new Object[][] {
                                        {"Unreplaceable_Key", "Unreplaceable_Value"},
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        {"Unreplaceable_Key", "Unreplaceable_Value"},
                                }
                        )
                },
                // Map with multiple elements, none of which are replaced:
                {
                        MapUtils.putAll(new HashMap<String, Object>(),
                                new Object[][] {
                                        {"Unreplaceable_Key1", "Unreplaceable_Value1"},
                                        {"Unreplaceable_Key2", "Unreplaceable_Value2"},
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        {"Unreplaceable_Key1", "Unreplaceable_Value1"},
                                        {"Unreplaceable_Key2", "Unreplaceable_Value2"},
                                }
                        )
                },
                // Singleton map containing a replaced value:
                {
                        MapUtils.putAll(new HashMap<String, Object>(),
                                new Object[][] {
                                        { MafOutputRendererConstants.FieldName_Other_Transcripts, OTHER_TRANSCRIPTS_RAW_NAMES },
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        { MafOutputRendererConstants.FieldName_Other_Transcripts, OTHER_TRANSCRIPTS_MAF_READY_NAMES.replaceAll(VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER, MafOutputRendererConstants.OTHER_TRANSCRIPT_DELIMITER)},
                                }
                        )
                },
                // Map with multiple elements, some of which are replaced:
                {
                        MapUtils.putAll(new HashMap<String, Object>(),
                                new Object[][] {
                                        {"Unreplaceable_Key1", "Unreplaceable_Value1"},
                                        { MafOutputRendererConstants.FieldName_Other_Transcripts, OTHER_TRANSCRIPTS_RAW_NAMES },
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        {"Unreplaceable_Key1", "Unreplaceable_Value1"},
                                        { MafOutputRendererConstants.FieldName_Other_Transcripts, OTHER_TRANSCRIPTS_MAF_READY_NAMES.replaceAll(VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER, MafOutputRendererConstants.OTHER_TRANSCRIPT_DELIMITER)},
                                }
                        )
                },
                // Map with multiple elements, all of which are replaced:
                {
                        MapUtils.putAll(new HashMap<String, Object>(),
                                new Object[][] {
                                        { MafOutputRendererConstants.FieldName_Other_Transcripts, OTHER_TRANSCRIPTS_RAW_NAMES },
                                        { MafOutputRendererConstants.FieldName_Variant_Type,      GencodeFuncotation.VariantType.INS.toString() },
                                        { MafOutputRendererConstants.FieldName_Reference_Allele,  "A" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1, "A" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, "AT" },
                                        { MafOutputRendererConstants.FieldName_Start_Position,    "100" },
                                        { MafOutputRendererConstants.FieldName_End_Position,      "100" },
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        { MafOutputRendererConstants.FieldName_Other_Transcripts, OTHER_TRANSCRIPTS_MAF_READY_NAMES.replaceAll(VcfOutputRenderer.OTHER_TRANSCRIPT_DELIMITER, MafOutputRendererConstants.OTHER_TRANSCRIPT_DELIMITER)},
                                        { MafOutputRendererConstants.FieldName_Variant_Type,      GencodeFuncotation.VariantType.INS.toString() },
                                        { MafOutputRendererConstants.FieldName_Reference_Allele,  "-" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1, "-" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, "T" },
                                        { MafOutputRendererConstants.FieldName_Start_Position,    "100" },
                                        { MafOutputRendererConstants.FieldName_End_Position,      "101" },
                                }
                        )
                }
        };
    }

    @DataProvider
    private Object[][] provideForAdjustIndelAlleleInformationForMafOutput() {
        return new Object[][] {
                // Empty maps:
                {
                        new LinkedHashMap<>(),
                        new LinkedHashMap<>()
                },
                // Singleton map that doesn't contain a replaceable element:
                {
                        MapUtils.putAll(new LinkedHashMap<String, Object>(),
                                new Object[][] {
                                        {"Unreplaceable_Key", "Unreplaceable_Value"},
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        {"Unreplaceable_Key", "Unreplaceable_Value"},
                                }
                        )
                },
                // Singleton map that contains a variant type that is not replaced:
                {
                        MapUtils.putAll(new LinkedHashMap<String, Object>(),
                                new Object[][] {
                                        {MafOutputRendererConstants.FieldName_Variant_Type, GencodeFuncotation.VariantType.SNP.toString()},
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        {MafOutputRendererConstants.FieldName_Variant_Type, GencodeFuncotation.VariantType.SNP.toString()},
                                }
                        )
                },
                // Map with multiple elements, none of which are replaced:
                {
                        MapUtils.putAll(new LinkedHashMap<String, Object>(),
                                new Object[][] {
                                        {"Unreplaceable_Key1", "Unreplaceable_Value1"},
                                        {MafOutputRendererConstants.FieldName_Variant_Type, GencodeFuncotation.VariantType.SNP.toString()},
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        {"Unreplaceable_Key1", "Unreplaceable_Value1"},
                                        {MafOutputRendererConstants.FieldName_Variant_Type, GencodeFuncotation.VariantType.SNP.toString()},
                                }
                        )
                },
                // Map contains an INSERTION of size 1 with associated fields:
                {
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        { MafOutputRendererConstants.FieldName_Variant_Type,      GencodeFuncotation.VariantType.INS.toString() },
                                        { MafOutputRendererConstants.FieldName_Reference_Allele,  "A" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1, "A" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, "AT" },
                                        { MafOutputRendererConstants.FieldName_Start_Position,    "100" },
                                        { MafOutputRendererConstants.FieldName_End_Position,      "100" },
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        { MafOutputRendererConstants.FieldName_Variant_Type,      GencodeFuncotation.VariantType.INS.toString() },
                                        { MafOutputRendererConstants.FieldName_Reference_Allele,  "-" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1, "-" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, "T" },
                                        { MafOutputRendererConstants.FieldName_Start_Position,    "100" },
                                        { MafOutputRendererConstants.FieldName_End_Position,      "101" },
                                }
                        )
                },
                // Map contains a DELETION of size 1 with associated fields:
                {
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        { MafOutputRendererConstants.FieldName_Variant_Type,      GencodeFuncotation.VariantType.DEL.toString() },
                                        { MafOutputRendererConstants.FieldName_Reference_Allele,  "AA" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1, "AA" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, "A" },
                                        { MafOutputRendererConstants.FieldName_Start_Position,    "100" },
                                        { MafOutputRendererConstants.FieldName_End_Position,      "100" },
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        { MafOutputRendererConstants.FieldName_Variant_Type,      GencodeFuncotation.VariantType.DEL.toString() },
                                        { MafOutputRendererConstants.FieldName_Reference_Allele,  "A" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1, "A" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, "-" },
                                        { MafOutputRendererConstants.FieldName_Start_Position,    "101" },
                                        { MafOutputRendererConstants.FieldName_End_Position,      "101" },
                                }
                        )
                },
                // Map contains an INSERTION of size 5 with associated fields:
                {
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        { MafOutputRendererConstants.FieldName_Variant_Type,      GencodeFuncotation.VariantType.INS.toString() },
                                        { MafOutputRendererConstants.FieldName_Reference_Allele,  "AAAAA" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1, "AAAAA" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, "AAAAATTTTT" },
                                        { MafOutputRendererConstants.FieldName_Start_Position,    "100" },
                                        { MafOutputRendererConstants.FieldName_End_Position,      "104" },
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        { MafOutputRendererConstants.FieldName_Variant_Type,      GencodeFuncotation.VariantType.INS.toString() },
                                        { MafOutputRendererConstants.FieldName_Reference_Allele,  "-" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1, "-" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, "TTTTT" },
                                        { MafOutputRendererConstants.FieldName_Start_Position,    "100" },
                                        { MafOutputRendererConstants.FieldName_End_Position,      "101" },
                                }
                        )
                },
                // Map contains a DELETION of size 5 with associated fields:
                {
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        { MafOutputRendererConstants.FieldName_Variant_Type,      GencodeFuncotation.VariantType.DEL.toString() },
                                        { MafOutputRendererConstants.FieldName_Reference_Allele,  "ATTTTT" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1, "ATTTTT" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, "A" },
                                        { MafOutputRendererConstants.FieldName_Start_Position,    "100" },
                                        { MafOutputRendererConstants.FieldName_End_Position,      "100" },
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        { MafOutputRendererConstants.FieldName_Variant_Type,      GencodeFuncotation.VariantType.DEL.toString() },
                                        { MafOutputRendererConstants.FieldName_Reference_Allele,  "TTTTT" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1, "TTTTT" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, "-" },
                                        { MafOutputRendererConstants.FieldName_Start_Position,    "101" },
                                        { MafOutputRendererConstants.FieldName_End_Position,      "105" },
                                }
                        )
                },
                // Map contains a DELETION of size 3 with associated fields:
                {
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        { MafOutputRendererConstants.FieldName_Variant_Type,      GencodeFuncotation.VariantType.DEL.toString() },
                                        { MafOutputRendererConstants.FieldName_Reference_Allele,  "CAGG" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1, "CAGG" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, "C" },
                                        { MafOutputRendererConstants.FieldName_Start_Position,    "2222103" },
                                        { MafOutputRendererConstants.FieldName_End_Position,      "2222103" },
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        { MafOutputRendererConstants.FieldName_Variant_Type,      GencodeFuncotation.VariantType.DEL.toString() },
                                        { MafOutputRendererConstants.FieldName_Reference_Allele,  "AGG" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele1, "AGG" },
                                        { MafOutputRendererConstants.FieldName_Tumor_Seq_Allele2, "-" },
                                        { MafOutputRendererConstants.FieldName_Start_Position,    "2222104" },
                                        { MafOutputRendererConstants.FieldName_End_Position,      "2222106" },
                                }
                        )
                },
        };
    }

    @DataProvider
    private Object[][] provideForMafTransform() {
        return new Object[][] {
                // Empty Strings:
                {
                    "", "", FuncotatorTestConstants.REFERENCE_VERSION_HG19, ""
                },

                // Field has no replacements
                {
                        "SPECIAL_NO_REPLACEY_FIELD", "GARBAGE", FuncotatorTestConstants.REFERENCE_VERSION_HG19, "GARBAGE"
                },

                // Field containing replaceable string from another field:
                {
                    MafOutputRendererConstants.FieldName_Variant_Classification,
                        MafOutputRendererConstants.FieldValue_Gencode_Chromosome_Mito,
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        MafOutputRendererConstants.FieldValue_Gencode_Chromosome_Mito
                },

                // Variant classification tests:
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        "In_Frame_Del"
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.IN_FRAME_INS.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        "In_Frame_Ins"
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        "Frame_Shift_Ins"
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        "Frame_Shift_Del"
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.MISSENSE.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        "Missense_Mutation"
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.NONSENSE.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        "Nonsense_Mutation"
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.SILENT.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        "Silent"
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.SPLICE_SITE.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        "Splice_Site"
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.START_CODON_DEL.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        "Translation_Start_Site"
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.NONSTOP.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        "Nonstop_Mutation"
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        "5'UTR"
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.THREE_PRIME_UTR.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        "3'UTR"
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        "5'Flank"
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.INTRON.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        "Intron"
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        GencodeFuncotation.VariantClassification.LINCRNA.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        "RNA"
                },

                // Chromosome tests:
                {
                    MafOutputRendererConstants.FieldName_Chromosome,
                    MafOutputRendererConstants.FieldValue_Gencode_Chromosome_Mito,
                    FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                    MafOutputRendererConstants.FieldValue_Chromosome_Mito
                },
                {
                        MafOutputRendererConstants.FieldName_Chromosome,
                        "chr7",
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        "7"
                },
                {
                        MafOutputRendererConstants.FieldName_Chromosome,
                        "ChR7",
                        "hG19",
                        "7"
                },
                {
                        MafOutputRendererConstants.FieldName_Chromosome,
                        "chracteristic",
                        "Hg19",
                        "chracteristic"
                },
                {
                        MafOutputRendererConstants.FieldName_Chromosome,
                        "chrZ",
                        "HG19",
                        "chrZ"
                },
                {
                        MafOutputRendererConstants.FieldName_Chromosome,
                        "chr7",
                        "b37",
                        "7"
                },
                {
                        MafOutputRendererConstants.FieldName_Chromosome,
                        "ChR7",
                        "B37",
                        "7"
                },
                {
                        MafOutputRendererConstants.FieldName_Chromosome,
                        "chr7",
                        FuncotatorTestConstants.REFERENCE_VERSION_HG38,
                        "chr7"
                },
                {
                        MafOutputRendererConstants.FieldName_Chromosome,
                        "ChR7",
                        "hG38",
                        "ChR7"
                },
                {
                        MafOutputRendererConstants.FieldName_Chromosome,
                        "chracteristic",
                        "Hg38",
                        "chracteristic"
                },
                {
                        MafOutputRendererConstants.FieldName_Chromosome,
                        "chrZ",
                        "HG38",
                        "chrZ"
                },
                {
                        MafOutputRendererConstants.FieldName_Chromosome,
                        "chrQ",
                        "Xenomorph A1",
                        "chrQ"
                },

                // Other transcripts:
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.IN_FRAME_INS.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.IN_FRAME_INS.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.MISSENSE.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.MISSENSE.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.NONSENSE.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.NONSENSE.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.SILENT.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.SILENT.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.SPLICE_SITE.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.SPLICE_SITE.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.START_CODON_DEL.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.START_CODON_DEL.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.NONSTOP.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.NONSTOP.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.THREE_PRIME_UTR.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.THREE_PRIME_UTR.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.INTRON.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.INTRON.toString())
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.LINCRNA.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.LINCRNA.toString())
                },

                // Fields that we do not replace:
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME.toString()
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME.toString()
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.START_CODON_SNP.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.START_CODON_SNP.toString()
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.START_CODON_INS.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.START_CODON_INS.toString()
                },

                // Multiple other transcripts:
                {
                    MafOutputRendererConstants.FieldName_Other_Transcripts,
                        OTHER_TRANSCRIPTS_RAW_NAMES,
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        OTHER_TRANSCRIPTS_MAF_READY_NAMES
                }
        };
    }

    @DataProvider
    private Object[][] provideForMafTransformInvert() {
        return new Object[][] {
                // Empty Strings:
                {
                        "", "", FuncotatorTestConstants.REFERENCE_VERSION_HG19, ""
                },

                // Field has no replacements
                {
                        "SPECIAL_NO_REPLACEY_FIELD", "GARBAGE", FuncotatorTestConstants.REFERENCE_VERSION_HG19, "GARBAGE"
                },

                // Field containing replaceable string from another field:
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        MafOutputRendererConstants.FieldValue_Gencode_Chromosome_Mito,
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        MafOutputRendererConstants.FieldValue_Gencode_Chromosome_Mito
                },

                // Variant classification tests:
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        "In_Frame_Del",
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        "In_Frame_Ins",
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.IN_FRAME_INS.toString()
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        "Frame_Shift_Ins",
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        "Frame_Shift_Del",
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        "Missense_Mutation",
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.MISSENSE.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        "Nonsense_Mutation",
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.NONSENSE.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        "Silent",
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.SILENT.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        "Splice_Site",
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.SPLICE_SITE.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        "Translation_Start_Site",
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.START_CODON_DEL.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        "Nonstop_Mutation",
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.NONSTOP.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        "5'UTR",
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        "3'UTR",
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.THREE_PRIME_UTR.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        "5'Flank",
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        "Intron",
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.INTRON.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Variant_Classification,
                        "RNA",
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.LINCRNA.toString(),
                },

                // Chromosome tests:
                {
                        MafOutputRendererConstants.FieldName_Chromosome,
                        MafOutputRendererConstants.FieldValue_Chromosome_Mito,
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        MafOutputRendererConstants.FieldValue_Gencode_Chromosome_Mito,
                },
                {
                        MafOutputRendererConstants.FieldName_Chromosome,
                        "7",
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        "chr7",
                },
                {
                        MafOutputRendererConstants.FieldName_Chromosome,
                        "acteristic",
                        "Hg19",
                        "acteristic",
                },
                {
                        MafOutputRendererConstants.FieldName_Chromosome,
                        "7",
                        "b37",
                        "chr7",
                },
                {
                        MafOutputRendererConstants.FieldName_Chromosome,
                        "chr7",
                        FuncotatorTestConstants.REFERENCE_VERSION_HG38,
                        "chr7"
                },
                {
                        MafOutputRendererConstants.FieldName_Chromosome,
                        "ChR7",
                        "hG38",
                        "ChR7"
                },
                {
                        MafOutputRendererConstants.FieldName_Chromosome,
                        "chracteristic",
                        "Hg38",
                        "chracteristic"
                },
                {
                        MafOutputRendererConstants.FieldName_Chromosome,
                        "chrZ",
                        "HG38",
                        "chrZ"
                },
                {
                        MafOutputRendererConstants.FieldName_Chromosome,
                        "chrQ",
                        "Xenomorph A1",
                        "chrQ"
                },

                // Other transcripts:
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString()),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString()),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.IN_FRAME_DEL.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.IN_FRAME_INS.toString()),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.IN_FRAME_INS.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS.toString()),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL.toString()),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.MISSENSE.toString()),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.MISSENSE.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.NONSENSE.toString()),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.NONSENSE.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.SILENT.toString()),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.SILENT.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.SPLICE_SITE.toString()),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.SPLICE_SITE.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.START_CODON_DEL.toString()),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.START_CODON_DEL.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.NONSTOP.toString()),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.NONSTOP.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR.toString()),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.FIVE_PRIME_UTR.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.THREE_PRIME_UTR.toString()),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.THREE_PRIME_UTR.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK.toString()),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.FIVE_PRIME_FLANK.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.INTRON.toString()),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.INTRON.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.LINCRNA.toString()),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.LINCRNA.toString(),
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        "LINCRNA",
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        "LINCRNA",
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        "LINCRNA|" + MafOutputRendererConstants.VariantClassificationMap.get(GencodeFuncotation.VariantClassification.LINCRNA.toString()),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        "LINCRNA|" + GencodeFuncotation.VariantClassification.LINCRNA.toString(),
                },


                // Fields that we do not replace:
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.DE_NOVO_START_IN_FRAME.toString()
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.DE_NOVO_START_OUT_FRAME.toString()
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.START_CODON_SNP.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.START_CODON_SNP.toString()
                },
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        GencodeFuncotation.VariantClassification.START_CODON_INS.toString(),
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        GencodeFuncotation.VariantClassification.START_CODON_INS.toString()
                },

                // Multiple other transcripts:
                {
                        MafOutputRendererConstants.FieldName_Other_Transcripts,
                        OTHER_TRANSCRIPTS_MAF_READY_NAMES,
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                        OTHER_TRANSCRIPTS_RAW_NAMES,
                }
        };
    }

    @DataProvider
    private Object[][] provideForWrite() {

        final VariantContext variant1 = new VariantContextBuilder().chr("chr1").start(25).stop(25).alleles("A", "T").make();
        final VariantContext variant2 = new VariantContextBuilder().chr("chrX").start(1759).stop(1759).alleles("G", "C").make();
        final VariantContext variant3 = new VariantContextBuilder().chr("chr3").start(89247).stop(89248).alleles("AA", "GG").make();
        final VariantContext variant4 = new VariantContextBuilder().chr("chr6").start(888996).stop(888999).alleles("TCCC", "T").make();
        final VariantContext variant5 = new VariantContextBuilder().chr("chr9").start(1067897).stop(1067899).alleles("AGG", "A", "AG", "AT").make();

        final List<VariantContext> variantList = Arrays.asList(variant1, variant2, variant3, variant4, variant5);

        final List<String> baseFieldNameList = new ArrayList<>(
                createMafOutputRenderer( getSafeNonExistentFile("GARBAGE"), FuncotatorTestConstants.REFERENCE_VERSION_HG19 ).getDefaultMap().keySet()
        );

        final int fieldSize = 10;
        
        // NOTE: The data field names must match data sources that are checked in for this to work in an expected way:
        final List<List<Funcotation>> funcotationList = Arrays.asList(
                Collections.singletonList(
                        TableFuncotation.create(
                                baseFieldNameList,
                                createFieldValuesFromNameList("A", baseFieldNameList, fieldSize),
                                Allele.create("T"),
                                GencodeFuncotationFactory.DEFAULT_NAME, null
                        )
                ),
                Collections.singletonList(
                        TableFuncotation.create(
                                baseFieldNameList,
                                createFieldValuesFromNameList("B", baseFieldNameList, fieldSize),
                                Allele.create("C"),
                                GencodeFuncotationFactory.DEFAULT_NAME, null
                        )
                ),
                Collections.singletonList(
                        TableFuncotation.create(
                                baseFieldNameList,
                                createFieldValuesFromNameList("C", baseFieldNameList, fieldSize),
                                Allele.create("GG"),
                                GencodeFuncotationFactory.DEFAULT_NAME, null
                        )
                ),
                Collections.singletonList(
                        TableFuncotation.create(
                                baseFieldNameList,
                                createFieldValuesFromNameList("D", baseFieldNameList, fieldSize),
                                Allele.create("T"),
                                "TestDataSource4", null
                        )
                ),
                Arrays.asList(
                        TableFuncotation.create(
                                baseFieldNameList,
                                createFieldValuesFromNameList("E", baseFieldNameList, fieldSize),
                                Allele.create("A"),
                                "TestDataSource5", null
                        ),
                        TableFuncotation.create(
                                baseFieldNameList,
                                createFieldValuesFromNameList("F", baseFieldNameList, fieldSize),
                                Allele.create("AG"),
                                "TestDataSource5", null
                        ),
                        TableFuncotation.create(
                                baseFieldNameList,
                                createFieldValuesFromNameList("G", baseFieldNameList, fieldSize),
                                Allele.create("AT"),
                                "TestDataSource5", null
                        )
                )
        );

        return new Object[][] {
                {
                    variantList,
                    funcotationList,
                    new File(FuncotatorTestConstants.FUNCOTATOR_TEST_DIR + File.separator + "ExampleMafFileForTests.maf")
                },
        };
    }

    //==================================================================================================================
    // Tests:

    @Test(dataProvider = "provideForTestReplaceFuncotationValuesWithMafCompliantValues")
    public void testReplaceFuncotationValuesWithMafCompliantValues(final Map<String, Object> outputMap, final LinkedHashMap<String, String> expected ) {
        final LinkedHashMap<String, String> compliantMap = createMafOutputRenderer(FuncotatorTestConstants.REFERENCE_VERSION_HG19).replaceFuncotationValuesWithMafCompliantValues(outputMap);
        Assert.assertEquals(compliantMap, expected);
    }

    @Test(dataProvider = "provideForMafTransform")
    public void testMafTransform(final String key, final String value, final String referenceVersion, final String expectedValue) {
        final String transformedValue = MafOutputRenderer.mafTransform(key, value, referenceVersion);
        Assert.assertEquals(transformedValue, expectedValue);
    }

    @Test(dataProvider = "provideForMafTransformInvert")
    public void testMafTransformInvert(final String key, final String value, final String referenceVersion, final String expectedValue) {
        final String transformedValue = MafOutputRenderer.mafTransformInvert(key, value, referenceVersion);
        Assert.assertEquals(transformedValue, expectedValue);
    }

    @Test(dataProvider = "provideForAdjustIndelAlleleInformationForMafOutput")
    public void testAdjustIndelAlleleInformationForMafOutput(final LinkedHashMap<String, String> outputMap, final LinkedHashMap<String, String> expected) {
        createMafOutputRenderer(FuncotatorTestConstants.REFERENCE_VERSION_HG19).adjustIndelAlleleInformationForMafOutput(outputMap);
        Assert.assertEquals(outputMap, expected);
    }

    @Test(dataProvider = "provideForWrite")
    public void testWrite(final List<VariantContext> variants, final List<List<Funcotation>> funcotations, final File expectedFile) {

        final File outFile = getSafeNonExistentFile("TestMafOutputFile");
        try ( final MafOutputRenderer mafOutputRenderer = createMafOutputRenderer( outFile, FuncotatorTestConstants.REFERENCE_VERSION_HG19 ) ) {
            for ( int i = 0 ; i < variants.size(); ++i ) {
                final FuncotationMap funcotationMap = FuncotationMap.createNoTranscriptInfo(funcotations.get(i));
                mafOutputRenderer.write(variants.get(i), funcotationMap);
            }
        }

        // Make sure our files are as we expect them to be:
        try {
            // For this test it is important to check that there is no trailing whitespace, hence
            // we are running it with doTrimWhitespace set to false:
            IntegrationTestSpec.assertEqualTextFiles(outFile, expectedFile, MafOutputRendererConstants.COMMENT_STRING, false);
        }
        catch (final IOException ex) {
            throw new GATKException("ERROR comparing text files: " + outFile.toURI().toString() + " and " + expectedFile.toURI().toString(), ex);
        }
    }

    /**
     * Creatse a dummy gencode funcotation and a dummy funcotation from a fake datasource named FAKEDATA.
     * Then exclude one of the FAKEDATA fields and make sure that it does not appear in the output.
     */
    @Test
    public void testCreateMafCompliantOutputMapExclusion() {
        final File outFile = getSafeNonExistentFile("TestMafOutputExclusionFile.maf");
        final String dummyTranscriptName = "FAKE00001.1";
        final VariantContext dummyVariantContext = FuncotatorTestUtils.createSimpleVariantContext(FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),"3", 1000000, 1000000, "C", "T");
        final GencodeFuncotation dummyGencodeFuncotation = (GencodeFuncotation) FuncotatorTestUtils.createDummyGencodeFuncotation(dummyTranscriptName, dummyVariantContext);
        final Set<String> excludedFields = Collections.singleton("FAKEDATA_FOO");
        try ( final MafOutputRenderer mafOutputRenderer = createMafOutputRenderer( outFile, FuncotatorTestConstants.REFERENCE_VERSION_HG19, excludedFields) ) {
            final FuncotationMap funcotationMap = FuncotationMap.createFromGencodeFuncotations(Collections.singletonList(dummyGencodeFuncotation));
            funcotationMap.add(dummyTranscriptName, FuncotatorTestUtils.createDummyTableFuncotation());
            mafOutputRenderer.write(dummyVariantContext, funcotationMap);
        }

        final AnnotatedIntervalCollection maf = AnnotatedIntervalCollection.create(outFile.toPath(), null);
        Assert.assertTrue(maf.getRecords().size() > 0);
        maf.getRecords().forEach(r -> Assert.assertFalse(r.hasAnnotation("FAKEDATA_FOO")));
        maf.getRecords().forEach(r -> Assert.assertTrue(r.hasAnnotation("FAKEDATA_BAR")));
    }

    @Test
    public void testCreateMafCompliantOutputMapSanitized() {
        final File outFile = getSafeNonExistentFile("TestMafOutputSanitized.maf");
        final String dummyTranscriptName = "FAKE00001.1";
        final VariantContext dummyVariantContext = FuncotatorTestUtils.createSimpleVariantContext(FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),"3", 1000000, 1000000, "C", "T");
        final GencodeFuncotation dummyGencodeFuncotation = (GencodeFuncotation) FuncotatorTestUtils.createDummyGencodeFuncotation(dummyTranscriptName, dummyVariantContext);
        final Set<String> excludedFields = Collections.emptySet();
        try ( final MafOutputRenderer mafOutputRenderer = createMafOutputRenderer( outFile, FuncotatorTestConstants.REFERENCE_VERSION_HG19, excludedFields) ) {
            final FuncotationMap funcotationMap = FuncotationMap.createFromGencodeFuncotations(Collections.singletonList(dummyGencodeFuncotation));
            funcotationMap.add(dummyTranscriptName, FuncotatorTestUtils.createDummyTableFuncotation());
            mafOutputRenderer.write(dummyVariantContext, funcotationMap);
        }

        final AnnotatedIntervalCollection maf = AnnotatedIntervalCollection.create(outFile.toPath(), null);
        Assert.assertTrue(maf.getRecords().size() > 0);
        maf.getRecords().forEach(r -> Assert.assertTrue(r.hasAnnotation("FAKEDATA_FOO")));
        maf.getRecords().forEach(r -> Assert.assertTrue(r.hasAnnotation("FAKEDATA_BAR")));
        maf.getRecords().forEach(r -> Assert.assertTrue(r.hasAnnotation("FAKEDATA_BAZ")));
        maf.getRecords().forEach(r -> Assert.assertEquals(r.getAnnotationValue("FAKEDATA_BAZ"), "_%09_YES_%0A_"));
    }

    @Test
    // Asserting that we can write a MAF file with a version 43 GencodeFuncotation without error
    public void testVersion43GencodeFuncotation() {
        final File outFile = getSafeNonExistentFile("TestMafOutputSanitized.maf");
        final String dummyTranscriptName = "FAKE00001.1";
        final VariantContext dummyVariantContext = FuncotatorTestUtils.createSimpleVariantContext(FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),"3", 1000000, 1000000, "C", "T");
        final GencodeFuncotation dummyGencodeFuncotation = (GencodeFuncotation) FuncotatorTestUtils.createDummyGencodeFuncotation(dummyTranscriptName, dummyVariantContext, "43");
        final Set<String> excludedFields = Collections.emptySet();
        try ( final MafOutputRenderer mafOutputRenderer = createMafOutputRenderer( outFile, FuncotatorTestConstants.REFERENCE_VERSION_HG19, excludedFields) ) {
            final FuncotationMap funcotationMap = FuncotationMap.createFromGencodeFuncotations(Collections.singletonList(dummyGencodeFuncotation));
            funcotationMap.add(dummyTranscriptName, FuncotatorTestUtils.createDummyTableFuncotation());
            mafOutputRenderer.write(dummyVariantContext, funcotationMap);
        }

        final AnnotatedIntervalCollection maf = AnnotatedIntervalCollection.create(outFile.toPath(), null);
        Assert.assertTrue(maf.getRecords().size() > 0);
        maf.getRecords().forEach(r -> Assert.assertTrue(r.hasAnnotation("FAKEDATA_FOO")));
        maf.getRecords().forEach(r -> Assert.assertTrue(r.hasAnnotation("FAKEDATA_BAR")));
        maf.getRecords().forEach(r -> Assert.assertTrue(r.hasAnnotation("FAKEDATA_BAZ")));
        maf.getRecords().forEach(r -> Assert.assertEquals(r.getAnnotationValue("FAKEDATA_BAZ"), "_%09_YES_%0A_"));
    }
}
