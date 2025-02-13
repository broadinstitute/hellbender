package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.testutils.FuncotatorReferenceTestUtils;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationFactory;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class TranscriptUtilsUnitTest extends GATKBaseTest {

    //==================================================================================================================
    // Private Static Members:

    private static final ReferenceDataSource refDataSourceHg19Ch3;
    private static final ReferenceDataSource refDataSourceB37;

    private static final List<AutoCloseable> autoCloseableList = new ArrayList<>();
    static {
        refDataSourceHg19Ch3 = ReferenceDataSource.of( IOUtils.getPath(FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref()) );
        refDataSourceB37 = ReferenceDataSource.of( IOUtils.getPath(b37Reference) );
        autoCloseableList.add( refDataSourceHg19Ch3 );
        autoCloseableList.add( refDataSourceB37 );
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
    // Private Members:

    //==================================================================================================================
    // Helper Methods:

    //==================================================================================================================
    // Data Providers:

    @DataProvider(name = "provideForTestExtractTrascriptFromReference")
    Object[][] provideForTestExtractTrascriptFromReference() {

        final SimpleInterval cntn4Interval = new SimpleInterval("chr3", 2140497, 3099645);
        final ReferenceContext refContextCntn4 = new ReferenceContext( refDataSourceHg19Ch3, cntn4Interval);

        return new Object[][] {
                // Trivial tests:
                { refContextCntn4, List.of(new SimpleInterval("chr3", 2140497, 2140497)), "A" },
                { refContextCntn4, List.of(new SimpleInterval("chr3", 2140497, 2140498)), "AG" },
                { refContextCntn4, List.of(new SimpleInterval("chr3", 2140497, 2140499)), "AGC" },
                { refContextCntn4, List.of(new SimpleInterval("chr3", 2140497, 2140500)), "AGCC" },
                { refContextCntn4, List.of(new SimpleInterval("chr3", 2140497, 2140501)), "AGCCG" },
                { refContextCntn4, List.of(new SimpleInterval("chr3", 2140498, 2140501)), "GCCG" },
                { refContextCntn4, List.of(new SimpleInterval("chr3", 2140499, 2140501)), "CCG" },
                { refContextCntn4, List.of(new SimpleInterval("chr3", 2140500, 2140501)), "CG" },
                { refContextCntn4, List.of(new SimpleInterval("chr3", 2140501, 2140501)), "G" },

                // Tests with real transcripts:
                {
                    refContextCntn4,
                    List.of(
                            new SimpleInterval("chr3", 2613188, 2613242),
                            new SimpleInterval("chr3", 2777899, 2778025),
                            new SimpleInterval("chr3", 2787206, 2787371),
                            new SimpleInterval("chr3", 2140497, 2140662),
                            new SimpleInterval("chr3", 2380862, 2380917),
                            new SimpleInterval("chr3", 2613100, 2613187)
                    ),
                    "AGCCGCCGCGAGCGCCGGGAGGGAGGGCAGCGCCGCGCCAGACGGAGCCCGGGGTTCGACGCCAGGATTGGCTGCAAGTAGGGAGCTTTCGCCGCCGCCCCGGGCCCCTCGGACTGTGCCGGCGCCGCACCCGAGGCTCTCGCCAGCCCGGCGCCCCGGTGCTGAGGAATCATTGACATAGAGTAACTCCACAGCATGTGTCTTCAAGAGCTTCCCTAAAAGATTAAAGGTTATACAAAACTTAAAAGAAGCAGCAATTCTATTCGCTTGTTATTGGACTTGAAACTCCCTTTGACCTCGGAAACTGAAGATGAGGTTGCCATGGGAACTGCTGGTACTGCAATCATTCATTTTGTGCCTTGCAGATGATTCCACACTGCATGGCCCGATTTTTATTCAAGAACCAAGTCCTGTAATGTTCCCTTTGGATTCTGAGGAGAAAAAAGTGAAGCTCAATTGTGAAGTTAAAGGAAATCCAAAACCTCATATCAGGTGGAAGTTAAATGGAACAGATGTTGACACTGGTATGGATTTCCGCTACAGTGTTGTTGAAGGGAGCTTGTTGATCAATAACCCCAATAAAACCCAAGATGCTGGAACGTACCAGTGCACAGCGACAAACTCGTTTGGAACAATTGTTAGCAGAGAAGCAAAGCTT"
                },
                {
                        refContextCntn4,
                        List.of(
                                new SimpleInterval("chr3", 2613188, 2613242),
                                new SimpleInterval("chr3", 2777899, 2778025),
                                new SimpleInterval("chr3", 2787206, 2787381),
                                new SimpleInterval("chr3", 2861170, 2861265),
                                new SimpleInterval("chr3", 2908436, 2908633),
                                new SimpleInterval("chr3", 2928724, 2928767),
                                new SimpleInterval("chr3", 2140550, 2140662),
                                new SimpleInterval("chr3", 2140958, 2141333),
                                new SimpleInterval("chr3", 2380862, 2380917),
                                new SimpleInterval("chr3", 2613100, 2613187),
                                new SimpleInterval("chr3", 2928768, 2928908),
                                new SimpleInterval("chr3", 2942369, 2942505),
                                new SimpleInterval("chr3", 2944560, 2944689),
                                new SimpleInterval("chr3", 2967313, 2967463),
                                new SimpleInterval("chr3", 3030029, 3030156),
                                new SimpleInterval("chr3", 3067786, 3067961),
                                new SimpleInterval("chr3", 3072542, 3072659),
                                new SimpleInterval("chr3", 3076316, 3076474),
                                new SimpleInterval("chr3", 3078863, 3079012),
                                new SimpleInterval("chr3", 3080617, 3080687),
                                new SimpleInterval("chr3", 3081721, 3081955),
                                new SimpleInterval("chr3", 3083994, 3084106),
                                new SimpleInterval("chr3", 3084661, 3084662)
                        ),
                        "GTTCGACGCCAGGATTGGCTGCAAGTAGGGAGCTTTCGCCGCCGCCCCGGGCCCCTCGGACTGTGCCGGCGCCGCACCCGAGGCTCTCGCCAGCCCGGCGCCCCGGTGCTGAGCCGGAAAATAAGTTTGTTGCGCTGCGAGGCAGCCACAAAACAAGGAACCGAGAGCCCGGAATGCTGCGGGAAGCCTTCAAGTCAGCTCCTCCGACTGGTTCGGGCTACTGCCCCCTCTCCGTGCGCCCTGGCCTCTGGCGCCGGGTTCCCGGCGGGGCTTTTCTTCTGACAGCCCAGTCACAGCCCGCAGCAGAGGGACGCGAACCTGGGGAGTGGAGGGACCTGGGACTAAAGGAACAGGAGCCCGTAGCCGTGGTGGAAGGAGCCGCGTGGAGACGGAGGCTGATGTCTGTGGCGCCCGCTGGGTGCCGGGCTGGCTGCTGAGCGCTGAGGCTGCGGCGGCGAGCGACAGGCCAGGTGCCTGCTCTTAGGGAAGGAATCATTGACATAGAGTAACTCCACAGCATGTGTCTTCAAGAGCTTCCCTAAAAGATTAAAGGTTATACAAAACTTAAAAGAAGCAGCAATTCTATTCGCTTGTTATTGGACTTGAAACTCCCTTTGACCTCGGAAACTGAAGATGAGGTTGCCATGGGAACTGCTGGTACTGCAATCATTCATTTTGTGCCTTGCAGATGATTCCACACTGCATGGCCCGATTTTTATTCAAGAACCAAGTCCTGTAATGTTCCCTTTGGATTCTGAGGAGAAAAAAGTGAAGCTCAATTGTGAAGTTAAAGGAAATCCAAAACCTCATATCAGGTGGAAGTTAAATGGAACAGATGTTGACACTGGTATGGATTTCCGCTACAGTGTTGTTGAAGGGAGCTTGTTGATCAATAACCCCAATAAAACCCAAGATGCTGGAACGTACCAGTGCACAGCGACAAACTCGTTTGGAACAATTGTTAGCAGAGAAGCAAAGCTTCAGTTTGCTTATCTTGACAACTTTAAAACAAGAACAAGAAGCACTGTGTCTGTCCGTCGAGGTCAAGGAATGGTGCTACTGTGTGGCCCGCCACCCCATTCTGGAGAGCTGAGTTATGCCTGGATCTTCAATGAATACCCTTCCTATCAGGATAATCGCCGCTTTGTTTCTCAAGAGACTGGGAATCTGTATATTGCCAAAGTAGAAAAATCAGATGTTGGGAATTATACCTGTGTGGTTACCAATACCGTGACAAACCACAAGGTCCTGGGGCCACCTACACCACTAATATTGAGAAATGATGTCCAGTACCAACTATTATCTGGCGAAGAGCTGATGGAAAGCCAATAGCAAGGAAAGCCAGAAGACACAAGTCAAATGGAATTCTTGAGATCCCTAATTTTCAGCAGGAGGATGCTGGTTTATATGAATGTGTAGCTGAAAATTCCAGAGGGAAAAATGTAGCAAGGGGACAGCTAACTTTCTATGCTCAACCTAATTGGATTCAAAAAATAAATGATATTCACGTGGCCATGGAAGAAAATGTCTTTTGGGAATGTAAAGCAAATGGAAGGCCTAAGCCTACATACAAGTGGCTAAAAAATGGCGAACCTCTGCTAACTCGGGATAGAATTCAAATTGAGCAAGGAACACTCAACATAACAATAGTGAACCTCTCAGATGCTGGCATGTATCAGTGTTTGGCAGAGAATAAACATGGAGTTATCTTTTCCAACGCAGAGCTTAGTGTTATAGCTGTAGGTCCAGATTTTTCAAGAACACTCTTGAAAAGAGTAACTCTTGTCAAAGTGGGAGGTGAAGTTGTCATTGAGTGTAAGCCAAAAGCGTCTCCAAAACCTGTTTACACCTGGAAGAAAGGAAGGGATATATTAAAAGAAAATGAAAGAATTACCATTTCTGAAGATGGAAACCTCAGAATCATCAACGTTACTAAATCAGACGCTGGGAGTTATACCTGTATAGCCACTAACCATTTTGGAACTGCTAGCAGTACTGGAAACTTGGTAGTGAAAGATCCAACAAGGGTAATGGTACCCCCTTCCAGTATGGATGTCACTGTTGGAGAGAGTATTGTTTTACCGTGCCAGGTAACGCATGATCACTCGCTAGACATCGTGTTTACTTGGTCATTTAATGGACACCTGATAGACTTTGACAGAGATGGGGACCACTTTGAAAGAGTTGGAGGGGATTCAGCTGGTGATTTGATGATCCGAAACATCCAACTGAAGCATGCTGGGAAATATGTCTGCATGGTCCAAACAAGTGTGGACAGGCTATCTGCTGCTGCAGACCTGATTGTAAGAGGTCCTCCAGGTCCCCCAGAGGCTGTGACAATAGACGAAATCACAGATACCACTGCTCAGCTCTCCTGGAGACCCGGGCCTGACAACCACAGCCCCATCACCATGTATGTCATTCAAGCCAGGACTCCATTCTCCGTGGGCTGGCAAGCAGTCAGTACAGTCCCAGAACTCATTGATGGGAAGACATTCACAGCGACCGTGGTGGGTTTGAACCCTTGGGTTGAATATGAATTCCGCACAGTTGCAGCCAACGTGATTGGGATTGGGGAGCCCAGCCGCCCCTCAGAGAAACGGAGAACAGAAGAAGCTCTCCCCGAAGTCACACCAGCGAATGTCAGTGGTGGCGGAGGCAGCAAATCTGAACTGGTTATAACCTGGGAGACGGTCCCTGAGGAATTACAGAATGGTCGAGGCTTTGGTTATGTGGTGGCCTTCCGGCCCTACGGTAAAATGATCTGGATGCTGACAGTGCTGGCCTCAGCTGATGCCTCTAGATACGTGTTCAGGAATGAGAGCGTGCACCCCTTCTCTCCCTTTGAGGTTAAAGTAGGTGTCTTCAACAACAAAGGAGAAGGCCCTTTCAGTCCCACCACGGTGGTGTATTCTGCAGAAGAAGAACCCACCAAACCACCAGCCAGTATCTTTGCCAGAAGTCTTTCTGCCACAGATATTGAAGTTTTCTGGGCCTCCCCACTGGAGAAGAATAGAGGACGAATACAAGGTTATGAGGT"
                },
                {
                        refContextCntn4,
                        List.of(
                                new SimpleInterval("chr3", 2613188, 2613242),
                                new SimpleInterval("chr3", 2777899, 2778007),
                                new SimpleInterval("chr3", 2140572, 2140662),
                                new SimpleInterval("chr3", 2142242, 2142323),
                                new SimpleInterval("chr3", 2304030, 2304114),
                                new SimpleInterval("chr3", 2380862, 2380917),
                                new SimpleInterval("chr3", 2613100, 2613187)
                        ),
                        "AAGTAGGGAGCTTTCGCCGCCGCCCCGGGCCCCTCGGACTGTGCCGGCGCCGCACCCGAGGCTCTCGCCAGCCCGGCGCCCCGGTGCTGAGTGGGTGAAAAAGAACAGTGTGTCATGAAGACAGGCACCTGGAGGTGATTTGGGTGGCATTCATGAGAAAATTCACGTTACCCACTTCGCTAACTGATGAAACACTGCCAGAGCCCACTTGGACTGTTTCGCAACTTTACTTAAGGGAATAAAAGCTCTCATGCTGAGGAATCATTGACATAGAGTAACTCCACAGCATGTGTCTTCAAGAGCTTCCCTAAAAGATTAAAGGTTATACAAAACTTAAAAGAAGCAGCAATTCTATTCGCTTGTTATTGGACTTGAAACTCCCTTTGACCTCGGAAACTGAAGATGAGGTTGCCATGGGAACTGCTGGTACTGCAATCATTCATTTTGTGCCTTGCAGATGATTCCACACTGCATGGCCCGATTTTTATTCAAGAACCAAGTCCTGTAATGTTCCCTTTGGATTCTGAGGAGAAAAAAGTGAAGCTCAATTGTGAAGTTAAAGGAAA"
                },
        };
    }

    @DataProvider(name = "provideForTestExtractTrascriptFromReferenceWithContigNameConversion")
    Object[][] provideForTestExtractTrascriptFromReferenceWithContigNameConversion() {

        final SimpleInterval cntn4Interval = new SimpleInterval("3", 2140497, 3099645);
        final ReferenceContext refContextCntn4 = new ReferenceContext( refDataSourceB37, cntn4Interval);

        return new Object[][] {
                // Trivial tests:
                { refContextCntn4, List.of(new SimpleInterval("chr3", 2140497, 2140497)), "A" },
                { refContextCntn4, List.of(new SimpleInterval("chr3", 2140497, 2140498)), "AG" },
                { refContextCntn4, List.of(new SimpleInterval("chr3", 2140497, 2140499)), "AGC" },
                { refContextCntn4, List.of(new SimpleInterval("chr3", 2140497, 2140500)), "AGCC" },
                { refContextCntn4, List.of(new SimpleInterval("chr3", 2140497, 2140501)), "AGCCG" },
                { refContextCntn4, List.of(new SimpleInterval("chr3", 2140498, 2140501)), "GCCG" },
                { refContextCntn4, List.of(new SimpleInterval("chr3", 2140499, 2140501)), "CCG" },
                { refContextCntn4, List.of(new SimpleInterval("chr3", 2140500, 2140501)), "CG" },
                { refContextCntn4, List.of(new SimpleInterval("chr3", 2140501, 2140501)), "G" },

                // Tests with real transcripts:
                {
                        refContextCntn4,
                        List.of(
                                new SimpleInterval("chr3", 2613188, 2613242),
                                new SimpleInterval("chr3", 2777899, 2778025),
                                new SimpleInterval("chr3", 2787206, 2787371),
                                new SimpleInterval("chr3", 2140497, 2140662),
                                new SimpleInterval("chr3", 2380862, 2380917),
                                new SimpleInterval("chr3", 2613100, 2613187)
                        ),
                        "AGCCGCCGCGAGCGCCGGGAGGGAGGGCAGCGCCGCGCCAGACGGAGCCCGGGGTTCGACGCCAGGATTGGCTGCAAGTAGGGAGCTTTCGCCGCCGCCCCGGGCCCCTCGGACTGTGCCGGCGCCGCACCCGAGGCTCTCGCCAGCCCGGCGCCCCGGTGCTGAGGAATCATTGACATAGAGTAACTCCACAGCATGTGTCTTCAAGAGCTTCCCTAAAAGATTAAAGGTTATACAAAACTTAAAAGAAGCAGCAATTCTATTCGCTTGTTATTGGACTTGAAACTCCCTTTGACCTCGGAAACTGAAGATGAGGTTGCCATGGGAACTGCTGGTACTGCAATCATTCATTTTGTGCCTTGCAGATGATTCCACACTGCATGGCCCGATTTTTATTCAAGAACCAAGTCCTGTAATGTTCCCTTTGGATTCTGAGGAGAAAAAAGTGAAGCTCAATTGTGAAGTTAAAGGAAATCCAAAACCTCATATCAGGTGGAAGTTAAATGGAACAGATGTTGACACTGGTATGGATTTCCGCTACAGTGTTGTTGAAGGGAGCTTGTTGATCAATAACCCCAATAAAACCCAAGATGCTGGAACGTACCAGTGCACAGCGACAAACTCGTTTGGAACAATTGTTAGCAGAGAAGCAAAGCTT"
                },
                {
                        refContextCntn4,
                        List.of(
                                new SimpleInterval("chr3", 2613188, 2613242),
                                new SimpleInterval("chr3", 2777899, 2778025),
                                new SimpleInterval("chr3", 2787206, 2787381),
                                new SimpleInterval("chr3", 2861170, 2861265),
                                new SimpleInterval("chr3", 2908436, 2908633),
                                new SimpleInterval("chr3", 2928724, 2928767),
                                new SimpleInterval("chr3", 2140550, 2140662),
                                new SimpleInterval("chr3", 2140958, 2141333),
                                new SimpleInterval("chr3", 2380862, 2380917),
                                new SimpleInterval("chr3", 2613100, 2613187),
                                new SimpleInterval("chr3", 2928768, 2928908),
                                new SimpleInterval("chr3", 2942369, 2942505),
                                new SimpleInterval("chr3", 2944560, 2944689),
                                new SimpleInterval("chr3", 2967313, 2967463),
                                new SimpleInterval("chr3", 3030029, 3030156),
                                new SimpleInterval("chr3", 3067786, 3067961),
                                new SimpleInterval("chr3", 3072542, 3072659),
                                new SimpleInterval("chr3", 3076316, 3076474),
                                new SimpleInterval("chr3", 3078863, 3079012),
                                new SimpleInterval("chr3", 3080617, 3080687),
                                new SimpleInterval("chr3", 3081721, 3081955),
                                new SimpleInterval("chr3", 3083994, 3084106),
                                new SimpleInterval("chr3", 3084661, 3084662)
                        ),
                        "GTTCGACGCCAGGATTGGCTGCAAGTAGGGAGCTTTCGCCGCCGCCCCGGGCCCCTCGGACTGTGCCGGCGCCGCACCCGAGGCTCTCGCCAGCCCGGCGCCCCGGTGCTGAGCCGGAAAATAAGTTTGTTGCGCTGCGAGGCAGCCACAAAACAAGGAACCGAGAGCCCGGAATGCTGCGGGAAGCCTTCAAGTCAGCTCCTCCGACTGGTTCGGGCTACTGCCCCCTCTCCGTGCGCCCTGGCCTCTGGCGCCGGGTTCCCGGCGGGGCTTTTCTTCTGACAGCCCAGTCACAGCCCGCAGCAGAGGGACGCGAACCTGGGGAGTGGAGGGACCTGGGACTAAAGGAACAGGAGCCCGTAGCCGTGGTGGAAGGAGCCGCGTGGAGACGGAGGCTGATGTCTGTGGCGCCCGCTGGGTGCCGGGCTGGCTGCTGAGCGCTGAGGCTGCGGCGGCGAGCGACAGGCCAGGTGCCTGCTCTTAGGGAAGGAATCATTGACATAGAGTAACTCCACAGCATGTGTCTTCAAGAGCTTCCCTAAAAGATTAAAGGTTATACAAAACTTAAAAGAAGCAGCAATTCTATTCGCTTGTTATTGGACTTGAAACTCCCTTTGACCTCGGAAACTGAAGATGAGGTTGCCATGGGAACTGCTGGTACTGCAATCATTCATTTTGTGCCTTGCAGATGATTCCACACTGCATGGCCCGATTTTTATTCAAGAACCAAGTCCTGTAATGTTCCCTTTGGATTCTGAGGAGAAAAAAGTGAAGCTCAATTGTGAAGTTAAAGGAAATCCAAAACCTCATATCAGGTGGAAGTTAAATGGAACAGATGTTGACACTGGTATGGATTTCCGCTACAGTGTTGTTGAAGGGAGCTTGTTGATCAATAACCCCAATAAAACCCAAGATGCTGGAACGTACCAGTGCACAGCGACAAACTCGTTTGGAACAATTGTTAGCAGAGAAGCAAAGCTTCAGTTTGCTTATCTTGACAACTTTAAAACAAGAACAAGAAGCACTGTGTCTGTCCGTCGAGGTCAAGGAATGGTGCTACTGTGTGGCCCGCCACCCCATTCTGGAGAGCTGAGTTATGCCTGGATCTTCAATGAATACCCTTCCTATCAGGATAATCGCCGCTTTGTTTCTCAAGAGACTGGGAATCTGTATATTGCCAAAGTAGAAAAATCAGATGTTGGGAATTATACCTGTGTGGTTACCAATACCGTGACAAACCACAAGGTCCTGGGGCCACCTACACCACTAATATTGAGAAATGATGTCCAGTACCAACTATTATCTGGCGAAGAGCTGATGGAAAGCCAATAGCAAGGAAAGCCAGAAGACACAAGTCAAATGGAATTCTTGAGATCCCTAATTTTCAGCAGGAGGATGCTGGTTTATATGAATGTGTAGCTGAAAATTCCAGAGGGAAAAATGTAGCAAGGGGACAGCTAACTTTCTATGCTCAACCTAATTGGATTCAAAAAATAAATGATATTCACGTGGCCATGGAAGAAAATGTCTTTTGGGAATGTAAAGCAAATGGAAGGCCTAAGCCTACATACAAGTGGCTAAAAAATGGCGAACCTCTGCTAACTCGGGATAGAATTCAAATTGAGCAAGGAACACTCAACATAACAATAGTGAACCTCTCAGATGCTGGCATGTATCAGTGTTTGGCAGAGAATAAACATGGAGTTATCTTTTCCAACGCAGAGCTTAGTGTTATAGCTGTAGGTCCAGATTTTTCAAGAACACTCTTGAAAAGAGTAACTCTTGTCAAAGTGGGAGGTGAAGTTGTCATTGAGTGTAAGCCAAAAGCGTCTCCAAAACCTGTTTACACCTGGAAGAAAGGAAGGGATATATTAAAAGAAAATGAAAGAATTACCATTTCTGAAGATGGAAACCTCAGAATCATCAACGTTACTAAATCAGACGCTGGGAGTTATACCTGTATAGCCACTAACCATTTTGGAACTGCTAGCAGTACTGGAAACTTGGTAGTGAAAGATCCAACAAGGGTAATGGTACCCCCTTCCAGTATGGATGTCACTGTTGGAGAGAGTATTGTTTTACCGTGCCAGGTAACGCATGATCACTCGCTAGACATCGTGTTTACTTGGTCATTTAATGGACACCTGATAGACTTTGACAGAGATGGGGACCACTTTGAAAGAGTTGGAGGGGATTCAGCTGGTGATTTGATGATCCGAAACATCCAACTGAAGCATGCTGGGAAATATGTCTGCATGGTCCAAACAAGTGTGGACAGGCTATCTGCTGCTGCAGACCTGATTGTAAGAGGTCCTCCAGGTCCCCCAGAGGCTGTGACAATAGACGAAATCACAGATACCACTGCTCAGCTCTCCTGGAGACCCGGGCCTGACAACCACAGCCCCATCACCATGTATGTCATTCAAGCCAGGACTCCATTCTCCGTGGGCTGGCAAGCAGTCAGTACAGTCCCAGAACTCATTGATGGGAAGACATTCACAGCGACCGTGGTGGGTTTGAACCCTTGGGTTGAATATGAATTCCGCACAGTTGCAGCCAACGTGATTGGGATTGGGGAGCCCAGCCGCCCCTCAGAGAAACGGAGAACAGAAGAAGCTCTCCCCGAAGTCACACCAGCGAATGTCAGTGGTGGCGGAGGCAGCAAATCTGAACTGGTTATAACCTGGGAGACGGTCCCTGAGGAATTACAGAATGGTCGAGGCTTTGGTTATGTGGTGGCCTTCCGGCCCTACGGTAAAATGATCTGGATGCTGACAGTGCTGGCCTCAGCTGATGCCTCTAGATACGTGTTCAGGAATGAGAGCGTGCACCCCTTCTCTCCCTTTGAGGTTAAAGTAGGTGTCTTCAACAACAAAGGAGAAGGCCCTTTCAGTCCCACCACGGTGGTGTATTCTGCAGAAGAAGAACCCACCAAACCACCAGCCAGTATCTTTGCCAGAAGTCTTTCTGCCACAGATATTGAAGTTTTCTGGGCCTCCCCACTGGAGAAGAATAGAGGACGAATACAAGGTTATGAGGT"
                },
                {
                        refContextCntn4,
                        List.of(
                                new SimpleInterval("chr3", 2613188, 2613242),
                                new SimpleInterval("chr3", 2777899, 2778007),
                                new SimpleInterval("chr3", 2140572, 2140662),
                                new SimpleInterval("chr3", 2142242, 2142323),
                                new SimpleInterval("chr3", 2304030, 2304114),
                                new SimpleInterval("chr3", 2380862, 2380917),
                                new SimpleInterval("chr3", 2613100, 2613187)
                        ),
                        "AAGTAGGGAGCTTTCGCCGCCGCCCCGGGCCCCTCGGACTGTGCCGGCGCCGCACCCGAGGCTCTCGCCAGCCCGGCGCCCCGGTGCTGAGTGGGTGAAAAAGAACAGTGTGTCATGAAGACAGGCACCTGGAGGTGATTTGGGTGGCATTCATGAGAAAATTCACGTTACCCACTTCGCTAACTGATGAAACACTGCCAGAGCCCACTTGGACTGTTTCGCAACTTTACTTAAGGGAATAAAAGCTCTCATGCTGAGGAATCATTGACATAGAGTAACTCCACAGCATGTGTCTTCAAGAGCTTCCCTAAAAGATTAAAGGTTATACAAAACTTAAAAGAAGCAGCAATTCTATTCGCTTGTTATTGGACTTGAAACTCCCTTTGACCTCGGAAACTGAAGATGAGGTTGCCATGGGAACTGCTGGTACTGCAATCATTCATTTTGTGCCTTGCAGATGATTCCACACTGCATGGCCCGATTTTTATTCAAGAACCAAGTCCTGTAATGTTCCCTTTGGATTCTGAGGAGAAAAAAGTGAAGCTCAATTGTGAAGTTAAAGGAAA"
                },
        };
    }

    //==================================================================================================================
    // Tests:

    @Test(dataProvider = "provideForTestExtractTrascriptFromReference")
    public void TestExtractTrascriptFromReference(final ReferenceContext refContext, final List<Locatable> exons, final String expectedTranscript) {
        final String transcript = TranscriptUtils.extractTrascriptFromReference(refContext, exons);

        if ( !transcript.equals(expectedTranscript) ) {
            System.out.println(transcript);
            System.out.println(expectedTranscript);
        }
        Assert.assertEquals(transcript, expectedTranscript);
    }

    @Test(dataProvider = "provideForTestExtractTrascriptFromReferenceWithContigNameConversion")
    public void TestExtractTrascriptFromReferenceWithContigNameConversion(final ReferenceContext refContext, final List<Locatable> exons, final String expectedTranscript) {
        final String transcript = TranscriptUtils.extractTrascriptFromReference(refContext, exons, true);

        if ( !transcript.equals(expectedTranscript) ) {
            System.out.println(transcript);
            System.out.println(expectedTranscript);
        }
        Assert.assertEquals(transcript, expectedTranscript);
    }

}
