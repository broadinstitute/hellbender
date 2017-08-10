package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype.InternalFilterLongReadAlignmentsSAMSpark.getCanonicalChromosomes;
import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertTrue;

public class InternalFilterLongReadAlignmentsSAMSparkUnitTest extends BaseTest {

    private static final Set<String> canonicalChromosomes = getCanonicalChromosomes(null, null);

    @DataProvider(name = "contigAlignmentsHeuristicFilter")
    private Object[][] createTestData() throws IOException {

        final List<Object[]> data = new ArrayList<>(20);

        AlignmentInterval intervalOne = new AlignmentInterval(new SimpleInterval("chr1", 1948156, 1948936),
                1, 787, TextCigarCodec.decode("257M4I182M2I342M361S"), true, 60, 8, 733, false, false);
        AlignmentInterval intervalTwo = new AlignmentInterval(new SimpleInterval("chr1", 1948935, 1949190),
                893, 1148, TextCigarCodec.decode("892H256M"), true, 60, 3, 241, false, false);
        AlignedContig contig = new AlignedContig("asm000063:tig00003", "CCACTGTGCCCGGCCAAGGGTCCCCGGTTCTGAAAGTGGAAGGGGTGCGGCTGCCTCAGGAGTCACCACGGCAACAAGAACCTGGACCTGAGCGCAGGTGGTCAGATTCTGGGGCCAGCAGCTTTTTGGTTTTTAGAGACGAGGTCTCACTCTGTTGCCCAGGCTGGAGTGCAGTGGTGCGATCACTGCACCCTGCAGCCTCGGCCTCCTGGTTTCAAGTGACCACAGATGCATGCAGCCATGCTTGGCATATATAAATATATATATATATATATATTTATGTGTATATTGGTAGAGACATGGTCTTGTTATATTGCCCAGGCTGATCGCAAACATCTGCTTAAGCGATCCTCCTGCGTTGGCTCTCCAAAGTATTGGGATTATAGGCATGAGCTACCATGGCCTGGCCTCCTTATTCTAGTCTTTTCTTTCCTTTCTTCTTGTTTTTTTTTTTTTTTTGGCAGGGTCTCACTCTGTCACCCAGGCTGCAGTGCAGTGGTGTGATCACAGCTCACTGCAGCCTCAACTTCCCAGGCTCAAGCGATCCTCCCGGCTCAGCATCCTGAGTAGCTGGGACTACAGATGCATGTCACCACGCCTGGCTAAATTTTCTTCTTTGTAGATATGGGGTCTCACCATGTAGTACTTTTCAATGTATTAAGCATCCTTATTTGATATTTGATGCCTGATAATACCCATGTCTGAACCATGCAAGATTGCTGCAATTCCTTCCTTCCTTCCCTCCCTCCTTCCCTTCCTTCCTTCCCTTTCCTTCCTTCCTCTTTCCCTCCCTTCTTTCCTTCCCTTTCCCTCCCTCCCTTCCTTCCTCTTTCCTTCCTTCCTTTCCCTCCCTTACTCCTTCCTTCCCTTCCCCTTCCTTCTTCCTTCTCTCCCTCCCTCCCTTCCCCTCCCTTACTCCCTTCCTTCCTCCTTCCCTCCCTCCTTTCCTTCATTCCCTTCCTTCCCCTTCCCCTTCCTTCCTTCTCTCCCTCCCTCCTTCCTTCCCTCCTTTCCTTCCTTCCTTCCTTTCCTTTCCCTCCTTCCTCCCTCCCTCCTTTCCTTCCTTCCTTTCCTTTCCTCCCTTCCCTCCCTCCCTCCCTCCCTTCCTTCCCCTCCCTCCCTCCTTTCCTTCTTTCGACAGAGTCTTG".getBytes(),
                Arrays.asList(intervalOne, intervalTwo), false);
        data.add(new Object[]{contig, Arrays.asList(intervalOne), Arrays.asList(intervalOne, intervalTwo), 1, 2});

        intervalOne = new AlignmentInterval(new SimpleInterval("chr2", 1422222, 1422435),
                1, 270,  TextCigarCodec.decode("75M56I139M"), false, 60, 56, 142, false, false);
        intervalTwo = new AlignmentInterval(new SimpleInterval("chr2_KI270774v1_alt", 105288, 105557),
                1, 270,  TextCigarCodec.decode("114M1I27M1I127M"), false, 56, 13, 179, false, false);
        contig = new AlignedContig("asm002608:tig00001", "ATGCTGGGGAATTTGTGTGCTCCTTGGGTGGGGACGAGCATGGAAGGCGCGTGGGACTGAAGCCTTGAAGACCCCGCAGGCGCCTCTCCTGGACAGACCTCGTGCAGGCGCCTCTCCTGGACCGACCTCGTGCAGGCGCCTCTCCTGGACAGACCTCGTGCAGGCGCCTCTCCTGGACCGACCTCGTGCAGGCGCCGCGCTGGACCGACCTCGTGCAGGCGCCGCGCTGGGCCATGGGGAGAGCGAGAGCCTGGTGTGCCCCTCAGGGAC".getBytes(),
                Arrays.asList(intervalOne, intervalTwo), true);
        data.add(new Object[]{contig, Arrays.asList(intervalTwo), Arrays.asList(intervalOne), 3, 1});

        intervalOne = new AlignmentInterval(new SimpleInterval("chr1", 30374719, 30375721),
                1, 1002,  TextCigarCodec.decode("966M1D36M2362H"), true, 60, 6, 960, false, false);
        intervalTwo = new AlignmentInterval(new SimpleInterval("chr1", 30375922, 30378473),
                826, 3364,  TextCigarCodec.decode("825S33M1D1047M7D553M5D906M"), true, 60, 24, 2423, false, false);
        AlignmentInterval intervalThree = new AlignmentInterval(new SimpleInterval("chr1_KI270760v1_alt", 22529, 23531),
                1, 1002,  TextCigarCodec.decode("966M1D36M2362H"), true, 14, 3, 975, false, false);
        AlignmentInterval intervalFour = new AlignmentInterval(new SimpleInterval("chr1_KI270760v1_alt", 23681, 26220),
                826, 3364,  TextCigarCodec.decode("825H33M1D2506M"), true, 60, 2, 2517, false, false);
        contig = new AlignedContig("asm027070:tig00000", "GAGCCCATCTCCTTGACTGTGGCTCTGATGCTGCCTCCACACTGGGATCTCTCTGCTCTCTTCACCTCATACCTCCTTCCCCCCACCTCACCCCATCGCCCCCGTTCTTGATCCTGCAATTGTAGAAACAGAAAGTTGGCTGATTTCTTGGGCCCGCAAATTGCCCAACAGGGAGACTGGGTGGGCGGCCCCCGCTTCCACTCCATCGCCCACCCTGATGCATCGTCTGACACTTTCAATTTATTTTTCAATTCCTCTACCATCAGAAATGACGATTAGATTTCCAGCATAAATACCGCCTTACCAAACTGAATTAATCACGGCAAGGAGGGGCACACACAGGCTCCAGCAGCCTGGGCAGAACATCCCCAGCATTAACCCTTCCGTCCTCACCCAGGCCCCCACCAGCAGGACGGAGGCTCCAGGCCTCACAGAAGACGCCACTCAAAATATCACTGGGGTCACCTAATCCCATCCCCCTTACCCTTTGCAGCCTCCCTCCTGTGGGAGTTCCTAGGAAGTGTCTTGCCCAAAGCCATCCACTCCATCAGGGCAGAGTCAGAGACACTGGCCCCTCATCTCCAGCCCCATCAGGGAAGGAGGCTCCATCCACATCCAGGACAAGATGTGGGAGTATCCGGGGTTTGGCGTTGTCCAGGACACATACGGGACGGGACTCCTGCAGACCCGAGGGTGGGGGCACCCAGTGATCACAGGGCCTGAACTGAAAGGGGTCTTGGAGAGACCTGGAGGCAGGTTCCAACCCTTGCCCCACAAACAAGACCATCACCCCTCTTTGCTGAGACTGTTCATTGCTCAGTCCAACAACCACAGCTCAGGTTGACCTCCAGCCTCCCCACTTCTCCACCTCCCTGACTCCAACCACAGCTCAGGGTGACCTCCAGCCTCCCCACTTCTCCACCTCCCTGACTCCAACCACAGCTCAGGGTGACATCCAGCCTCCCCACTTCTCCACCTCCCTGACTCCAGCCACAGCTCAGGCTCCTTCCTATGAGACCCCCATGGCCTCTCACAGCCTCTCCACTTCTATGCCTGTTCTCACCCAATCCCCATCCCTCAGCAGTCATCACCTCAAAATGCAAACACTGTCCTATGGTTTCCTGGCTCAGAACCCATCGGGCCCTCCTCTGCTCTCAAATCAGGCCCCCACCCTTCAAGGCCATGAGGACTGGGCTGGCCTGGCCCCTACCGGTCAGTGCACTCCCCCATCCTGGCTGGGTTGTCTCCTCTTTCTCCTTCAAGTTTTTCTATTTAAAATTCCCCTCCTCAGAGAACCTTCTCTGGCCACCATCCCCCAATCTAAATTAGGTTCTCCCTCCTAAGGTTCTTTCTCAAATCCATTTCCTTTCCTTCTGAGCACTTAAGCGAGCGATAATTACACACTAACTTGTGTAATTTGTTTAATAGGATCTTTGGGACAGAGACTTTATCTGACTCGCTTGATGCTGCAGCTGCTAGAACCCAGACCGTAATGTAGTGGGAGCTCAGTGCAGACTTTTGAAGGAGTAAGTGAGTAAAAGAACAACAAGCCCCTCTTGGTGCCCACCAAGTGCCAAGCTGAGACTGGGCCCTGGAGCTGGAGTCAAGATGTGGACCTGGCCTTGGTGTGCTGGGCCCTAACAGATGAGTAGGAGTTTGCCGAGCACTGAAGGTGGGGTTGACATGACCAACTTCTGAGAGGCACTCTTTGCCTCTGGATGGCCCCTTCCCAGTCACCCCAAAAGGAAGCCCTTGCCCTTTCAAAAGTGGTGAATGTGGTGGTTCAGATCGGTAGGTGTTCCTATGAATAGGTGAGGGGCCAGGCTTCAGGTCAGTTGAACCTGGGTTTGAATCCTGATTTTGCTCTTGGTACTAGGGCAGGTCACTGAGACGCTCTGAGCCTCTCTGCTCCAGGATGAGGATCCCTTCATCCATGCTCACTCAAAGTCCTGCCCACCAGGATGGAGGCAGACAGGCTGCAATGCCCTCCCCTCTCAGTGGGGGAAAAATACCAGGTCAGGCAGCCAGCAGCCGAGAATGCCAGGCAGAGCAAAGGTGTCCTAAGGGATGGACAGAATAAGGGCTTGAGAGCCTAGCCAAGGGTGAGGCTAGGAGAGGCTTCCCGGAGGACGAGGCAAGTCAGAGCTCTTTGCCTCTTACTCCCATGACTGTGGGTGCCTTTCTCCTCCTCCTCTCATTCTCTCTCCTTTCCAGCTCCTGCTCTGCTCATTTCTTCACCTCAGTCTCTCTGCCCCGACAGGAGCCCTGAGGGACACAACCCCGTCCCGAGGAATGTATCTGCCCACTTCCAGCAGGTTCCTGGAGGCCCTCTAAATTCCCCTTCCCCCCAAAGTCATCTCCCAACACTGCTGCTCCCAGGGTGGGACGCCTGCTGCTGCACCTCCACACACGTGCACACACCCAGCCAGGTGCAGACAGCGTGGGCAGTGCAGAGGGGAGGGCTGGGGATTAAGGAGTTCGTGTTCTTGAGCAGCCTGGAAAGCAGCAGGGCTTCCACAGGAGCCGCCCCTGCCCTCACCCCTGCCCAGTAGGGTTAAGGGGCTGGCTTAGATGTCACCCCAAGCCAAGGCTGTCCTTCTCAGAGGCTCCTTCCCAGCTCCCCTGAGTGGGTCAGTCCCTTCCCCTCTCTGAGCCCCTCTTTCCTCTTCTGTAAAGCAGACTCAGTGATGTTGCTCAGAGGATTGAAGGACAAAGAAAAGCAACACAATGGACAGCAGGGATTTGCAAACAGCCGGGTGCTGTACCCAAGACAGGGTATTGCTGGTGATGTCTGATGGATGGGGAGTTGAAAGACTCAGCTGTCACTGGGCAGCTGGGTCTGGTTCCCCTGAGTCATTCGTAATTCACCAACCCAGTCTATAGAAGCTTATTAAGCACTTATTGTGTGCCATGCTCCATGCAAGGGCCAAAGACACCATGAGCAGAGCCAGACCCCACCCTCAGGTTCCCCCATGGGATGGGGTTAGCCAGATGACCTGAAGGCCTCTCCAGCCAGCTCAACCCCCTTAATCCAGAATTACTCCCTGTGCCAGGCTGACGGTGTGGCCAGAGAGGCCAGGGCCTGGGAGGGGGCCTGGCAGTGGGTGGTGGGAAGAGATGGAGTGGCTGTGTCAGGGGAAGGAGAGAGCAGGTTGTTCCTGTACAGGTTTCGCTCCTCGGATAGGGGGCTGCAATGACAGCTTCCAGGAAAGACCAGGCAAGTGCCTCACCCCATCCATTCTTGCTCACCCCTGCGGCCTCTTGGCCAATGGCTGCTGTGACCCTGTCCTCCTCTGGGAATCTGGTCTCGGGGAGGAGCCCTGGACCCTGACATTGACTAGAAACCTGACCCCATGTCTGAGCA".getBytes(),
                Arrays.asList(intervalOne, intervalTwo, intervalThree, intervalFour), false);
        data.add(new Object[]{contig, Arrays.asList(intervalOne, intervalFour), Arrays.asList(intervalThree, intervalFour), 2, 2});

        return data.toArray(new Object[data.size()][]);
    }

    @Test(dataProvider = "contigAlignmentsHeuristicFilter", groups = "sv")
    public void testScoreComputation(final AlignedContig contig,
                                     final List<AlignmentInterval> configuration,
                                     final List<AlignmentInterval> configurationEquallyGoodOrBetter,
                                     final int expectedConfigurationCount,
                                     final int expectedAICount) {

        final double scoreOne = InternalFilterLongReadAlignmentsSAMSpark.computeScoreOfConfiguration(configuration, canonicalChromosomes, 60);
        final double equallyGoodOrBetterScore = InternalFilterLongReadAlignmentsSAMSpark.computeScoreOfConfiguration(configurationEquallyGoodOrBetter, canonicalChromosomes, 60);
        assertTrue( scoreOne <= equallyGoodOrBetterScore);
    }

    @Test(dataProvider = "contigAlignmentsHeuristicFilter", groups = "sv")
    public void testNewFiltering(final AlignedContig contig,
                                 final List<AlignmentInterval> configurationBetter,
                                 final List<AlignmentInterval> configurationWorse,
                                 final int expectedConfigurationCount,
                                 final int expectedAICount) {

        assertEquals(InternalFilterLongReadAlignmentsSAMSpark.pickBestConfigurations(contig, canonicalChromosomes).size(), expectedConfigurationCount);
    }

    @Test(dataProvider = "contigAlignmentsHeuristicFilter", groups = "sv")
    public void testIntegrative(final AlignedContig contig,
                                final List<AlignmentInterval> configurationBetter,
                                final List<AlignmentInterval> configurationWorse,
                                final int expectedConfigurationCount,
                                final int expectedAICount) {
        final AlignedContig tig =
                InternalFilterLongReadAlignmentsSAMSpark.filterAndSplitGappedAI(
                        SparkContextFactory.getTestSparkContext().parallelize(Collections.singletonList(contig))
                , null, null).collect().get(0);
        assertEquals(tig.alignmentIntervals.size(), expectedAICount,
                tig.alignmentIntervals.stream().map(AlignmentInterval::toPackedString).collect(Collectors.toList()).toString());
    }
}
