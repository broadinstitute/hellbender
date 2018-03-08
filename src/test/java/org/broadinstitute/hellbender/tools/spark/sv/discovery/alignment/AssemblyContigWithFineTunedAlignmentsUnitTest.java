package org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SVTestUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoverFromLocalAssemblyContigAlignmentsSpark.SAMFormattedContigAlignmentParser.parseReadsAndOptionallySplitGappedAlignments;

public class AssemblyContigWithFineTunedAlignmentsUnitTest extends GATKBaseTest {

    @Test(groups = "sv", dataProvider = "createData")
    void test(final AlignedContig alignedContig, final boolean expectedHasIncompletePicture) {

        Assert.assertEquals(new AssemblyContigWithFineTunedAlignments(alignedContig).hasIncompletePicture(), expectedHasIncompletePicture);
    }

    @DataProvider(name = "createData")
    private Object[][] createTestData() {
        final List<Object[]> data = new ArrayList<>(20);

        final List<String> bases1 = Arrays.asList("ACCAAGGGCAAAATTGTTCCACAGCATGAAAGAAATCCATAAGTTTTTCTGTATCAACTTTTACCCTACCATGCTTCAAGAGCTGCTGTAGCAAGCTCAAATACATGATGTACTTACTTTCAGTTTGTCCCATTTGTGTCCCTAGCTTTCTCCGAGTGCCCCGCTTACCTGCAGAGCTTGAAACTTTTTCATCCTTGGGAGTCCTTTGTCTGTTGGTCCTCTGTTTCACACACTTGAGTGTTCCTTCACCGGATTCTTTCAGGCCCCACGTTGGGCGCCAGAATGTTGGGGACCAGCCTCAACACCACCCGTAGGGTACCCAAAGTCCAATGGTGACAAAGGAATGAGAAAAGACAGGTTAAGAGTTCATAAAGGTGGGAGCCAGGGGACCAGTTGCAAAATGGAGGCTGCAAAAGGCTCAGAGCTCTGGTCTCCACACTATTTATTGGGTACAATCACTTAGATGTAAAAAGCAGATGTTCAGGGTGAAACAGTGAAAGGGTGGCAGTGCATCATAGGTGTAATTTATAGCAATAGTAGTTTAAATGAATCTCCTTGTGCTCAGTGTATCTTTAACTTATTGGAGAGTAGCTAGTGGGAGTGGGCTTAACTAGGAGCCTGCATGTCTGTCCGCATTCCCATGCTTCAAAGGAGTGTCTTTCTCCTGGAACACAGTGTTTACAAATAAGAGAGCGGGTCTCGCTCTGAGCATGGGAACATGATGGCAATTAGGAGGCTTTCCTCCTCAGAGGCCTTTTGTGGCTTTTCACAACTTATTTTCCTATATTTTTATGGCCAGTTTATACAGGCACCCCACAAGCCCTTTTCCCAACAATTTAAGTTTCCTGAGGCCTCCTGAGAAGCAGAAGCCACTATACTTCTTGTACAGCCTATAGAACCTCTTTTCTTTATAAATTACTCAGTTTCAGGTATTCCTTTATAACAACGCGAGAATGGACTAATAAAATCACCACATCATCATCCTCATCACCATTATCATCACCATCAACATCACTGTCATCATCACCATCATCATTATCACCATCACCATCATTATAGTTACCATCACCATAATTATCACCATCAACATCACCATCGGTGCCATCAACATCATCATCATCACCATCACCATCATCATCATCATTACCATTACCATCACCATCATCATCTTCACCATCACCATTGCACCATCACCATTATCATCATCATAATCACCATCAACATCACCATAATAATCATCACCATTATCATCATTATCACCATCACCATCATCATCATCATTACCATTACCATCACCATCATCATCTTCACCATCACCATCACCATTGCACCATCACCATTATCATCATCTTAATCACTATCAACATCACCATAATAATCACCATCACCATCACCATTATCATCATTATCACCATCACCATTATCACCATCATAATCATCATGAAATTACCATAATCATACTCACCATTACTATCACTGTCATCACCATCACAATCATTACCACTAACACCACCACCATCATCATCATCATCACCACCACCATCATCATCATTATCCTAACAATAAAGATGGCAGAACAAATGAATCTCATTTGTTAATGCCCCACATGGTGGCACTGTGCTGAGGGCCCTTATCTGCTACATCTCATCACTCAGCCTTATGTTGCTTCCCTCAGACCTTTAAGGATTTCTAAGCAAAATGGGATAACCTTATTCCTGGAAGAATTAATTGCTTCTTGAGTCAATAAACTATATTGAGGACCTAGTATGTGGTGGGCATCCATCAGAGAGCAAAACCAGGTGTGGTTCCTGCCCTCGTGGAGCTTACAGTCCAGTGGGGAGACAGATATTACTCATTACATACCGAATGGACACTTACAGATAGAGGTAAGTAGCTTGATAGAAAGTTCCATGGGGTTGGCCAGGTGTGGTGGCTCATGCCTGTAATCCCAGCATTTTGGGAGGCTGAGGTGGGTGGATCACAAGGTCAGGAGTTTGAGACCAGCCTGGCCAATATGGTGAAACCCTGTCTCTACTAAAAATACAAAAATTACCTGGGTGTGGTGGTGCGGGCCTGTAGTCCCAGTTACTTGGGAGGCTGAGGCAGAAGAATCGCTTAAACCCGGGAGGCGGAGGTTGCAGTTAGCCAAGATTGCACCACTGCACTCCAGCCTGGGTGACAGAGTGAGACTCCATCTCAA",
                "ATCACCATCACCATTATCATCATTATCACCATCACCATTATCACCATCATAATCATCATGAAATTACCATAATCATACTCACCATTACTATCACTGTCATCACCATCACAATCATTACCACTAACACCACCACCATCATCATCATCATCACCACCACCATCATCATCATTATCCTAACAATAAAGATGGCAGAACAAATGAATCTCATTTGTTAATGCCCCACATGGTGGCACTGTGCTGAGGGCCCTTATCTGCTACATCTCATCACTCAGCCTTATGTTGCTTCCCTCAGACCTTTAAGGATTTCTAAGCAAAATGGGATAACCTTATTCCTGGAAGAATTAATTGCTTCTTGAGTCAATAAACTATATTGAGGACCTAGTATGTGGTGGGCATCCATCAGAGAGCAAAACCAGGTGTGGTTCCTGCCCTCGTGGAGCTTACAGTCCAGTGGGGAGACAGATATTACTCATTACATACCGAATGGACACTTACAGATAGAGGTAAGTAGCTTGATAGAAAGTTCCATGGGGTTGGCCAGGTGTGGTGGCTCATGCCTGTAATCCCAGCATTTTGGGAGGCTGAGGTGGGTGGATCACAAGGTCAGGAGTTTGAGACCAGCCTGGCCAATATGGTGAAACCCTGTCTCTACTAAAAATACAAAAATTACCTGGGTGTGGTGGTGCGGGCCTGTAGTCCCAGTTACTTGGGAGGCTGAGGCAGAAGAATCGCTTAAACCCGGGAGGCGGAGGTTGCAGTTAGCCAAGATTGCACCACTGCACTCCAGCCTGGGTGACAGAGTGAGACTCCATCTCAA",
                "CCATCAACATCATCATCATCACCATCACCATCATCATCATCATTACCATTACCATCACCATCATCATC");
        final List<String> cigars1 = Arrays.asList("1017M3I52M1128S", "1383H817M", "1101H68M1031H");
        final List<String> chromosomes1 = Arrays.asList("chr1", "chr1", "chr7");
        final List<Integer> positions1 = Arrays.asList(19931877, 19932932, 355947);
        final List<Boolean> reverseStrands1 = Arrays.asList(false, false, false);
        final List<Boolean> suppStatus1 = Arrays.asList(false, true, true);

        final List<SAMRecord> reads1 = createSAMRecordsWithEssentialInfo(bases1, cigars1, chromosomes1, positions1, reverseStrands1, suppStatus1);

        final AlignedContig alignedContig1 = parseReadsAndOptionallySplitGappedAlignments(reads1, GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY, true);

        data.add(new Object[]{alignedContig1, false});



        final List<String> bases2 = Arrays.asList("CTGGAGTGCAATGGCATGATATTGGCTCACTGCAACCTCCACCTCCTGGGTTCAAGCAATTCTCCTGCCTCAGCTTCTAGAGTAGCTGGGATTACAGGTGCACACTACCACGCCCAGCTAATTTTTGTATTTTTATTAGAGATGGGGTTTCATCATGTTGGTCAGGCTGGTCTCGAACTCCTGACCTCAGGTGATTGTCCTGCCTCAGCATCCCAAAGTGCTGGGATTACAGGCATGAGGCACCGCGCCCAGCCAGCATGGAGGTATTTGAGAGCAACAGTGATCAGAACCATTTGGTTCAAGCAGCGGTTTTAAAACGGAAGTGGAGAAGGAATTAGCAGATCCCTGACATCCTCTTCAATCAGAGTTCCTCCATTGTGAACTGGTTTACATGTCAGCATTATGGATTTTGGTGCAACACCTGCCCCCAACAGGAAGAAAAGAAGAAAAAGAAAGAAGAGGAAGGAAGAAGAGAAAGACAAAGAAGAAGAAGGAGGAGGAGGCGGCAGGAGGAGTAGGAGGGAGGAAGAAGGAGGAGGAAGAAAAAGAGGAAGAAGAAAGGAGGAAGGAAGAAGAAAGAAGAATTGGGAGGCTGAGACAGGCAGATCACAAGGTCAGGAGTTCAAGACCAGCCTGGCCAACATGGTGAAACCCCATCACTTCTAAAAATACAAAAATTAGCTGGGCGTGTTGGCACATGCCTATAATCCCAGCTACTTGGGAGGCTGAGGCAGGAGAATTGCTTGAACCTGGGAGGCGGAGGTTGCAGTGAGCAGAGAGCTCGCCACTGCACTCCAGTCTGGGCAATGAGCGAGACTGTCTTGAAAAAAAAAAAAGGAAGGAAAGAAGGAGGAGGAGGAGAAGAAAGAAGAAGTTTTATTATTGTTATTTTTTGAGATGGAGTCTTGCTCTGTAGCCCAGGCTGGAATACAGTGGCACACTCTTGACTCCCTGCCACCTCTGACTGCTGGGTTCAAAGGAGAAGGCGGAAGAAGAAGGAAGAAGAAGAAAGAAGGAGAAAGGCTGGGTGCAGTGGCTCACACCTGTAATCTCAGCACTTTGGGAGGCCGAGGCAGGTGAATCACAAGGTAAGGAGTTCGAGACCAGCCTGGCCAACTGGTGAAACCCTGCCTCTACTAAAAGTACAAAAATTAGCCGGGCGTGGTCTGAGGCAGGAGAATCGCTTTAACCTGGGAGGAGGAGGTTGCAGTCAGTCAAGATGGCGCCACTGCACTCCAGCCTGGGTGACAGAGTGAGACTTTGTCTCAAAAAAAAAAAAGGAAAGAAGGAGAAGAAAGAAGGAGGAGGAGGAGAAGAAGAAGAGGAAAAGGAGGGGGAGGAGAAGGGGAGGGGGAGGAAAGAGGAGGAGAAGAAAGAAGAAGGAGAAGGAAGAAAGGAGGGAGGAGGAGGAGGAGAAGGGGGAGGGGGAGGAAGGAGGAGGAGAAGAATGAGGAGAAGGAGAAGAAGAAGGAAGAAGAAAGAAGGAAGAAGAAGGAGGAGGAGGAGGGGGAGGAGGGGGAGGGGGAGGAGGAGGAGGAGAGAAGGAGGAGAAAAGTAGTTGAGGCCCAAACACCAAGAGGGAGCAAAGATTGAAAAGATGAGATGAGCCATGAAAGCAAGTACAGGAGTTACTGATGGTACTGGGGAGCCCGTGTAGGTTTAGGATCTGAGCTTTTGGAAGATTGATTGGGTAGCCTTTGAGCCACCTGATAAGTGGAAAGAACAAGAGAGGCTGGATCTGTGTTTTCAGGAAGCATATGTTGGCCCAGCAATTGCTGGCTTGTAGTGAGGAGGCACACAACTGGCCTAGGACAGTGGTCATGAAAATGCAGAGGAGGTAAAGTCCCTGCACTCCTAGGGAGACTAGTCCTGATGTCAGTCTGGAGTCAGTCAGAATGGTGTCCTCTCCCTCCCTGCACTACCCAGCCCAGTCAGTGGGAGGACTTCCTCAATTCCAGTAGCCATTCAAGTCCCTGGAATTGGTGGCTGTCACTTGCAAACTATAGCCACTTGAGCAGAAATGGGCCAGGATTACTTATCTTTAATCTGCATATCATTGGGAGGCACTTACCTGCTAGCTCTGGCTAAAAACTAGAGCAACCCTGGCCTGCCGTAGCTCCTGCTGCCCAGACAACTCCTCCAATATGAAAGGGATGAGGGGAACTCAAAGTTACAATGTCCTACTTGGAGCAGTAAGTTCAGTAGACATATCACTTGCCTCATTAACATCAAGCATCCCAAAACCCAGTCTGGGTCAGTTTTGCCCAGAGTGGGGTTTGTAGAACACGGGTTCTCCTGGGATCCTATACCTAGCCCAGAATCAGTTGCAAAAGCCAGGCCATAGCGAATTGTCCTGCCAGCCAGATAGCAGAGAATCTGACGGCAGCAGGCAGAAGGAGCCGCTCCATTGCAGTAAGCCAAGATCGCGCCACTTGCCTCATTACATCAAGCATCCCAAAACCCAGTCTGGGTCAGTTTTGCCCA",
                "GGAGGGGGAGGAGGAGGAGGAGAGAAGGAGGAGAAAAGTAGTTGAGGCCCAAACACCAAGAGGGAGCAAAGATTGAAAAGATGAGATGAGCCATGAAAGCAAGTACAGGAGTTACTGATGGTACTGGGGAGCCCGTGTAGGTTTAGGATCTGAGCTTTTGGAAGATTGATTGGGTAGCCTTTGAGCCACCTGATAAGTGGAAAGAACAAGAGAGGCTGGATCTGTGTTTTCAGGAAGCATATGTTGGCCCAGCAATTGCTGGCTTGTAGTGAGGAGGCACACAACTGGCCTAGGACAGTGGTCATGAAAATGCAGAGGAGGTAAAGTCCCTGCACTCCTAGGGAGACTAGTCCTGATGTCAGTCTGGAGTCAGTCAGAATGGTGTCCTCTCCCTCCCTGCACTACCCAGCCCAGTCAGTGGGAGGACTTCCTCAATTCCAGTAGCCATTCAAGTCCCTGGAATTGGTGGCTGTCACTTGCAAACTATAGCCACTTGAGCAGAAATGGGCCAGGATTACTTATCTTTAATCTGCATATCATTGGGAGGCACTTACCTGCTAGCTCTGGCTAAAAACTAGAGCAACCCTGGCCTGCCGTAGCTCCTGCTGCCCAGACAACTCCTCCAATATGAAAGGGATGAGGGGAACTCAAAGTTACAATGTCCTACTTGGAGCAGTAAGTTCAGTAGACATATCACTTGCCTCATTAACATCAAGCATCCCAAAACCCAGTCTGGGTCAGTTTTGCCCAGAGTGGGGTTTGTAGAACACGGGTTCTCCTGGGATCCTATACCTAGCCCAGAATCAGTTGCAAAAGCCAGGCCATAGCGAATTGTCCTGCCAGCCAGATAGCAGAGAATCTGACGGCAGCAGGCAGAAGGAGCCGCTCCATTGCAGTAAGCCAAGATCGCGCCACTTGCCTCATTACATCAAGCATCCCAAAACCCAGTCTGGGTCAGTTTTGCCCA",
                "CCTCCTCCTCCTCCCCCTCCCCCTCCTCCCCCTCCTCCTCCTCCT");
        final List<String> cigars2 = Arrays.asList("1313M1171S", "1517H967M", "947H45M1492H");
        final List<String> chromosomes2 = Arrays.asList("chr1", "chr1", "chr14");
        final List<Integer> positions2 = Arrays.asList(39043258, 39044558, 70910125);
        final List<Boolean> reverseStrands2 = Arrays.asList(false, false, true);
        final List<Boolean> suppStatus2 = Arrays.asList(false, true, true);

        final List<SAMRecord> reads2 = createSAMRecordsWithEssentialInfo(bases2, cigars2, chromosomes2, positions2, reverseStrands2, suppStatus2);

        final AlignedContig alignedContig2 = parseReadsAndOptionallySplitGappedAlignments(reads2, GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY, true);
        data.add(new Object[]{alignedContig2, false});


        final List<String> bases3 = Arrays.asList("CCCGCCTCAGCCTCCCTAAGTGCTGAGATTACAGGCCTGAGCCACTGCGCCAGGCCTGGTTTTTTGGTTTCAAACCACAATAGACATTGCTGGAGAATCAAGCTCATAGTTTCTTTTTACTCTGCATGATATCCCTCCAAAAGCTTGTCTATTCTCATGACTTCATGACAGTTCTTTGCCAATGATTTGCAAACATTATCTCCAGTCTTGATTCTTTCCTAGGCTTTATTTCTAAATACTCACTAGTCATTTCCACTTAGAAGCTTTGTCTTCTTTTCAAACTCAGCATATCCAAAACTGACCTCATCTTGTTGCCCTAACTAGAACATGCCACAACCCATTTCTGCCTATAATGTCATTATTCTCTTAAGGCCCCCATATTTGCAAATCAGGAGTCATCTTTGCCTCCCTTCTCTGAATGCTGCTGTGCTCACAGTCAGTACCATACAGTTAGTACTTGTTCTCAGTTTAATGCACACTTGTTTCATATGACCTCCTAGGGCAGAACCAAGATTCGTGGGTGGAAGTTGCAGGGAGGCAGATTTGCCTCCATATAAGGAAGAACTTTTTAATACTTTGAGCTGTCTGAATGGAATGGACTGCCTCCAGAAGTCTCGGGTTCTCCATCACTGAAGGTGTTTGAGCAGAGGCTGCCTGACTAATTGCTAATTGAAAGGACTGTTCTGGCATAAGATGGTGTTTTTACAAGATGATTTCTAGAATCCCTCTAATCCTGAGAGCCAGTGAGTCGATAGAAGGTAGCTTTGTCTCTCCTGCTAGACTCCCTTAGGACAGGGAGACTATTTTACCTTTCTTTTATATTCTGTACAGCACTTAATTCAGGTGCTGGTCTCTTAATTGCCTAAAGATGATTATTTACAGGTTAATTGATTCTTTTCATTTTGTTCCAATATTTGGTTAAACACCAAATATTGTGGATTTTTTTCCTTTGAAATATCTTCTGTAGTCTGGGCACGGTGGCTCATGCCTGTAATCCCAGCACTTGGCTCACCGCAATTACAGGCATGAGCCACCGTGCCCAGACTACAGAAGATATTTCAAAGGAAATATCTTGGGGAGACCAAGGCAGGTGGATCACCTGAGATCAGGAGTTCGAGACCAGCCAGGCCAACATGGCGAAACCTTGTCTCTACTAAAAATACAAAAATTCGCCGGGCGTGGTGGTGCATGCCTGTAGTTCCAGCTACTCGGGAGGCTGAGGCAGGAGAATCGCTTGAACCTGGAAGGTGGAGGTTGCAGTGAGCTGAGATCATGCCATCACACTCCAGCCTGGGCAACAGAGTGAGACTCCATCTCAAAAAATTACAGTAATAATAAATAAATACAAATATCTTTTGTAGTAGCATTATTTTAGATGAATAGCAGCTTTTAACCTAACTTTCTAGACCCTAATATCATCACTTGTCTCCACATCCCATTCTGTTTATCAGACAGTACTTGTAAGCTCTTTGCTTTCATCATATTATTCTTTCACTCAAAACCCTTGAATGGCTTCACATTGATCACAACGTCAAATTTAAATTTTCCCTGGCTTTCAGTGCCCTCCATCATCTGGCCCAGCTTGCTCATGCATCCTTGTTCCCTATAATACTCCATGTTCCATGCAGGCTAACCCACTCACATTTTCTGTACATAGCCTGCTTGGTCTCCTTTTTTGCTTTTGACTTGCTCTGGAATGGCTTCCCTTTTTTCTCTTGCCTCTTCAAGATGCCTCACTTCCCTCCTGAAACCAAGCTACCGCCAGTCCTCATTGATTTCCTCTGAACTCTCAGAGCATGTAGTAATTTATGTAATTCAGAACCGTAGCAACAACATTCTGTGAAGAAAAATCTGCAAGAATAGGCTGATAATTTAACTTTCCCTAATCCAACTGGATATTCCCATAATAAAACTTTTAAAAATATAGGCTGGCTGTCATGGTTCACATCTGTAATCCCAGCATTTTGGGAGGCTGAGGCAAGAGGACTGTTTGAGCCCAGGAGTTTGAGACTAGCCTGAGCAACATAGTGAGACTCTGTCTCTATCACACACATA",
                "AAGTGCTGGGATTACAGGCATGAGCCACCGTGCCCAGACTAC",
                "CCCGCCTCAGCCTCCCTAAGTGCTGAGATTACAGGCCTGAGCCACTGCGCCAGGCCTGGTTTTTTGGTTTCAAACCACAATAGACATTGCTGGAGAATCAAGCTCATAGTTTCTTTTTACTCTGCATGATATCCCTCCAAAAGCTTGTCTATTCTCATGACTTCATGACAGTTCTTTGCCAATGATTTGCAAACATTATCTCCAGTCTTGATTCTTTCCTAGGCTTTATTTCTAAATACTCACTAGTCATTTCCACTTAGAAGCTTTGTCTTCTTTTCAAACTCAGCATATCCAAAACTGACCTCATCTTGTTGCCCTAACTAGAACATGCCACAACCCATTTCTGCCTATAATGTCATTATTCTCTTAAGGCCCCCATATTTGCAAATCAGGAGTCATCTTTGCCTCCCTTCTCTGAATGCTGCTGTGCTCACAGTCAGTACCATACAGTTAGTACTTGTTCTCAGTTTAATGCACACTTGTTTCATATGACCTCCTAGGGCAGAACCAAGATTCGTGGGTGGAAGTTGCAGGGAGGCAGATTTGCCTCCATATAAGGAAGAACTTTTTAATACTTTGAGCTGTCTGAATGGAATGGACTGCCTCCAGAAGTCTCGGGTTCTCCATCACTGAAGGTGTTTGAGCAGAGGCTGCCTGACTAATTGCTAATTGAAAGGACTGTTCTGGCATAAGATGGTGTTTTTACAAGATGATTTCTAGAATCCCTCTAATCCTGAGAGCCAGTGAGTCGATAGAAGGTAGCTTTGTCTCTCCTGCTAGACTCCCTTAGGACAGGGAGACTATTTTACCTTTCTTTTATATTCTGTACAGCACTTAATTCAGGTGCTGGTCTCTTAATTGCCTAAAGATGATTATTTACAGGTTAATTGATTCTTTTCATTTTGTTCCAATATTTGGTTAAACACCAAATATTGTGGATTTTTTTCCTTTGAAATATCTT");
        final List<String> cigars3 = Arrays.asList("1064S991M", "1050H42M963H", "961M1094H");
        final List<String> chromosomes3 = Arrays.asList("chr1", "chr16", "chr1");
        final List<Integer> positions3 = Arrays.asList(40050406, 48151, 40049455);
        final List<Boolean> reverseStrands3 = Arrays.asList(true, false, true);
        final List<Boolean> suppStatus3 = Arrays.asList(false, true, true);

        final List<SAMRecord> reads3 = createSAMRecordsWithEssentialInfo(bases3, cigars3, chromosomes3, positions3, reverseStrands3, suppStatus3);

        final AlignedContig alignedContig3 = parseReadsAndOptionallySplitGappedAlignments(reads3, GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY, true);
        data.add(new Object[]{alignedContig3, false});

        return data.toArray(new Object[data.size()][]);
    }

    private static List<SAMRecord> createSAMRecordsWithEssentialInfo(final List<String> readBases,
                                                                     final List<String> cigars,
                                                                     final List<String> chromosomes,
                                                                     final List<Integer> positions,
                                                                     final List<Boolean> isReverseStrand,
                                                                     final List<Boolean> isSupplementary) {
        Utils.validateArg(readBases.size() == cigars.size() &&
                        cigars.size() == chromosomes.size() &&
                        chromosomes.size() == positions.size() &&
                        positions.size() == isSupplementary.size() &&
                        isReverseStrand.size() == isSupplementary.size(),
                        "input of different sizes");

        return IntStream.range(0, readBases.size())
                .mapToObj(i -> {
                    final byte[] dummyQuals = new byte[readBases.get(i).length()];
                    Arrays.fill(dummyQuals, (byte)'A');
                    final GATKRead artificialRead = ArtificialReadUtils.createArtificialRead(readBases.get(i).getBytes(), dummyQuals, cigars.get(i));
                    artificialRead.setPosition(chromosomes.get(i), positions.get(i));
                    artificialRead.setIsReverseStrand(isReverseStrand.get(i));
                    artificialRead.setIsSupplementaryAlignment(isSupplementary.get(i));
                    return artificialRead;
                })
                .map(read -> read.convertToSAMRecord(null))
                .collect(Collectors.toList());
    }

    @Test(groups = "sv")
    public void testSerialization() {
        final AlignmentInterval one = SVTestUtils.fromSAMRecordString("asm002362:tig00002\t16\tchr2\t1422222\t60\t75M56I139M\t*\t0\t0\tATGCTGGGGAATTTGTGTGCTCCTTGGGTGGGGACGAGCATGGAAGGCGCGTGGGACTGAAGCCTTGAAGACCCCGCAGGCGCCTCTCCTGGACAGACCTCGTGCAGGCGCCTCTCCTGGACCGACCTCGTGCAGGCGCCTCTCCTGGACAGACCTCGTGCAGGCGCCTCTCCTGGACCGACCTCGTGCAGGCGCCGCGCTGGACCGACCTCGTGCAGGCGCCGCGCTGGGCCATGGGGAGAGCGAGAGCCTGGTGTGCCCCTCAGGGAC\t*\tSA:Z:chr2_KI270774v1_alt,105288,-,114M1I27M1I127M,56,13;\tMD:Z:214\tRG:Z:GATKSVContigAlignments\tNM:i:56\tAS:i:142\tXS:i:0\n",
                true);
        final AlignmentInterval two = SVTestUtils.fromSAMRecordString("asm002362:tig00002\t2064\tchr2_KI270774v1_alt\t105288\t56\t114M1I27M1I127M\t*\t0\t0\tATGCTGGGGAATTTGTGTGCTCCTTGGGTGGGGACGAGCATGGAAGGCGCGTGGGACTGAAGCCTTGAAGACCCCGCAGGCGCCTCTCCTGGACAGACCTCGTGCAGGCGCCTCTCCTGGACCGACCTCGTGCAGGCGCCTCTCCTGGACAGACCTCGTGCAGGCGCCTCTCCTGGACCGACCTCGTGCAGGCGCCGCGCTGGACCGACCTCGTGCAGGCGCCGCGCTGGGCCATGGGGAGAGCGAGAGCCTGGTGTGCCCCTCAGGGAC\t*\tSA:Z:chr2,1422222,-,75M56I139M,60,56;\tMD:Z:94C17G1G6T13T3G1G34A3T9T68T8\tRG:Z:GATKSVContigAlignments\tNM:i:13\tAS:i:179\tXS:i:142",
                true);
        final AlignedContig sourceTig = new AlignedContig("asm002362:tig00002", "GTCCCTGAGGGGCACACCAGGCTCTCGCTCTCCCCATGGCCCAGCGCGGCGCCTGCACGAGGTCGGTCCAGCGCGGCGCCTGCACGAGGTCGGTCCAGGAGAGGCGCCTGCACGAGGTCTGTCCAGGAGAGGCGCCTGCACGAGGTCGGTCCAGGAGAGGCGCCTGCACGAGGTCTGTCCAGGAGAGGCGCCTGCGGGGTCTTCAAGGCTTCAGTCCCACGCGCCTTCCATGCTCGTCCCCACCCAAGGAGCACACAAATTCCCCAGCAT".getBytes(),
                Arrays.asList(one, two), false);
        final List<AssemblyContigAlignmentsConfigPicker.GoodAndBadMappings> config = AssemblyContigAlignmentsConfigPicker.pickBestConfigurations(sourceTig, new HashSet<>(Collections.singletonList("chr2")), 0.);
        final AssemblyContigWithFineTunedAlignments tig = AssemblyContigAlignmentsConfigPicker.reConstructContigFromPickedConfiguration(new Tuple2<>(new Tuple2<>(sourceTig.contigName, sourceTig.contigSequence),
                config)).next();

        final ByteArrayOutputStream bos = new ByteArrayOutputStream();
        final Output out = new Output(bos);
        final Kryo kryo = new Kryo();
        kryo.writeClassAndObject(out, tig);
        out.flush();

        final ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
        final Input in = new Input(bis);
        @SuppressWarnings("unchecked")
        final AssemblyContigWithFineTunedAlignments roundTrip = (AssemblyContigWithFineTunedAlignments) kryo.readClassAndObject(in);
        Assert.assertEquals(tig, roundTrip);
    }
}
