package org.broadinstitute.hellbender.utils.read;


import htsjdk.samtools.*;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class GATKReadAdaptersUnitTest extends GATKBaseTest {

    public static final String BASIC_READ_NAME = "basic_read";
    public static final String BASIC_READ_CONTIG = "1";
    public static final int BASIC_READ_START = 5;
    public static final int BASIC_READ_MAPPING_QUALITY = 20;
    public static final byte[] BASIC_READ_BASES = {'A', 'C', 'G', 'T'};
    public static final byte[] BASIC_READ_BASE_QUALITIES = {30, 40, 30, 50};
    public static final String BASIC_READ_CIGAR = "1M1I2M";
    public static final String BASIC_READ_GROUP = "Foo";
    public static final String BASIC_PROGRAM = "x";
    public static final int BASIC_READ_END = BASIC_READ_START + TextCigarCodec.decode(BASIC_READ_CIGAR).getReferenceLength() - 1;
    public static final String BASIC_READ_MATE_CONTIG = "1";
    public static final int BASIC_READ_MATE_START = 125;
    private static final String BASIC_SAMPLE_NAME = "Sample1";

    private static GATKRead basicReadBackedBySam() {
        return new SAMRecordToGATKReadAdapter(basicSAMRecord());
    }

    private static SAMFileHeader getSAMHeader() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(2, 1, 1000000);
        final SAMReadGroupRecord readGroupRecord = new SAMReadGroupRecord(BASIC_READ_GROUP);
        readGroupRecord.setSample(BASIC_SAMPLE_NAME);
        header.addReadGroup(readGroupRecord);
        return header;
    }

    private static SAMRecord basicSAMRecord() {
        final SAMFileHeader header = getSAMHeader();

        final SAMRecord read = ArtificialReadUtils.createArtificialSAMRecord(
                header,
                BASIC_READ_NAME,
                header.getSequenceIndex(BASIC_READ_CONTIG),
                BASIC_READ_START,
                BASIC_READ_BASES,
                BASIC_READ_BASE_QUALITIES,
                BASIC_READ_CIGAR
        );
        read.setAttribute(SAMTag.RG.name(), BASIC_READ_GROUP);
        read.setMappingQuality(BASIC_READ_MAPPING_QUALITY);
        read.setMateReferenceName(BASIC_READ_MATE_CONTIG);
        read.setMateAlignmentStart(BASIC_READ_MATE_START);
        read.setReadPairedFlag(true);
        read.setFirstOfPairFlag(true);
        read.setReferenceName(BASIC_READ_CONTIG);
        // ArtificialReadUtils adds this tag, but explicitly add it for symmetry with basicGoogleGenomicsRead
        read.setAttribute(SAMTag.PG.name(), BASIC_PROGRAM);

        return read;
    }

    private static SAMRecord basicUnmappedSAMRecord() {
        SAMRecord sam = basicSAMRecord();
        sam.setReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
        return sam;
    }

    private static List<GATKRead> getUnmappedReads() {
        List<GATKRead> unmappedReads = new ArrayList<>();

        final SAMRecord unmappedFlagSam = basicSAMRecord();
        unmappedFlagSam.setReadUnmappedFlag(true);
        unmappedReads.add(new SAMRecordToGATKReadAdapter(unmappedFlagSam));

        final SAMRecord unmappedContigSam = basicSAMRecord();
        unmappedContigSam.setReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
        unmappedReads.add(new SAMRecordToGATKReadAdapter(unmappedContigSam));

        final SAMRecord noAlignmentStartSam = basicSAMRecord();
        noAlignmentStartSam.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
        unmappedReads.add(new SAMRecordToGATKReadAdapter(noAlignmentStartSam));

        return unmappedReads;
    }

    @DataProvider(name = "GetAndSetPositionData")
    public Object[][] getAndSetPositionData() {
        List<Object[]> testCases = new ArrayList<>();

        testCases.add(new Object[]{basicReadBackedBySam(), BASIC_READ_CONTIG, BASIC_READ_START, BASIC_READ_END});

        for ( GATKRead unmappedRead : getUnmappedReads() ) {
            testCases.add(new Object[]{unmappedRead, null, ReadConstants.UNSET_POSITION, ReadConstants.UNSET_POSITION});
        }

        return testCases.toArray(new Object[][]{});
    }

    @Test(dataProvider = "GetAndSetPositionData")
    public void testGetAndSetPosition( final GATKRead read, final String expectedContig, final int expectedStart, final int expectedEnd ) {
        Assert.assertEquals(read.getContig(), expectedContig, "Wrong contig for read");
        Assert.assertEquals(read.getStart(), expectedStart, "Wrong start position for read");
        Assert.assertEquals(read.getEnd(), expectedEnd, "Wrong end position for read");

        read.setPosition("2", 17);
        Assert.assertFalse(read.isUnmapped(), "Read should be mapped after setPosition() call");
        Assert.assertEquals(read.getContig(), "2", "Wrong contig after setPosition()");
        Assert.assertEquals(read.getStart(), 17, "Wrong start position after setPosition()");
    }

    @DataProvider(name = "InvalidSetPositionData")
    public Object[][] invalidSetPositionData() {
        return new Object[][]{
                { null, 1 },
                { SAMRecord.NO_ALIGNMENT_REFERENCE_NAME, 1 },
                { "1", SAMRecord.NO_ALIGNMENT_START }
        };
    }

    @Test(dataProvider = "InvalidSetPositionData", expectedExceptions = IllegalArgumentException.class)
    public void testSetInvalidPositionOnSamBackedRead( final String contig, final int start ) {
        final GATKRead read = new SAMRecordToGATKReadAdapter(basicSAMRecord());
        read.setPosition(contig, start);
    }

    @DataProvider(name = "GetAssignedPositionData")
    public Object[][] getAssignedPositionData() {
        List<Object[]> testCases = new ArrayList<>();

        testCases.add(new Object[]{basicReadBackedBySam(), BASIC_READ_CONTIG, BASIC_READ_START});

        // SAMRecord, unmapped flag set, but has an assigned position
        final SAMRecord unmappedFlagSam = basicSAMRecord();
        unmappedFlagSam.setReadUnmappedFlag(true);
        testCases.add(new Object[]{new SAMRecordToGATKReadAdapter(unmappedFlagSam), BASIC_READ_CONTIG, BASIC_READ_START});

        // SAMRecord, unmapped flag not set, but no contig
        final SAMRecord unmappedContigSam = basicSAMRecord();
        unmappedContigSam.setReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
        testCases.add(new Object[]{new SAMRecordToGATKReadAdapter(unmappedContigSam), SAMRecord.NO_ALIGNMENT_REFERENCE_NAME, BASIC_READ_START});

        // SAMRecord, unmapped flag not set, but no start
        final SAMRecord noAlignmentStartSam = basicSAMRecord();
        noAlignmentStartSam.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
        testCases.add(new Object[]{new SAMRecordToGATKReadAdapter(noAlignmentStartSam), BASIC_READ_CONTIG, ReadConstants.UNSET_POSITION});

        return testCases.toArray(new Object[][]{});
    }

    @Test(dataProvider = "GetAssignedPositionData")
    public void testGetAssignedPosition( final GATKRead read, final String expectedAssignedContig, final int expectedAssignedStart ) {
        Assert.assertEquals(read.getAssignedContig(), expectedAssignedContig, "Wrong assigned contig for read");
        Assert.assertEquals(read.getAssignedStart(), expectedAssignedStart, "Wrong assigned start position for read");
    }

    @DataProvider(name = "GetAndSetNameData")
    public Object[][] getAndSetNameData() {
        final SAMRecord namelessSam = basicSAMRecord();
        namelessSam.setReadName(null);

        return new Object[][]{
                { basicReadBackedBySam(), BASIC_READ_NAME },
                { new SAMRecordToGATKReadAdapter(namelessSam), null },
        };
    }

    @Test(dataProvider = "GetAndSetNameData")
    public void testGetAndSetName( final GATKRead read, final String expectedName ) {
        Assert.assertEquals(read.getName(), expectedName, "Wrong name for read");

        read.setName("NEWNAME");
        Assert.assertEquals(read.getName(), "NEWNAME", "Wrong name for read after setName()");

        read.setName(null);
        Assert.assertEquals(read.getName(), null, "Wrong name for read after setName(null)");
    }

    @DataProvider(name = "GetLengthData")
    public Object[][] getLengthData() {
        final SAMRecord baselessSam = basicSAMRecord();
        baselessSam.setReadBases(SAMRecord.NULL_SEQUENCE);


        return new Object[][]{
                { basicReadBackedBySam(), BASIC_READ_BASES.length },
                { new SAMRecordToGATKReadAdapter(baselessSam), 0 },
        };
    }

    @Test(dataProvider = "GetLengthData")
    public void testGetLength( final GATKRead read, final int expectedLength ) {
        Assert.assertEquals(read.getLength(), expectedLength, "Wrong length for read");
    }

    @DataProvider(name = "GetUnclippedStartAndEndData")
    public Object[][] getUnclippedStartAndEndData() {
        final SAMRecord softClippedSam = basicSAMRecord();
        softClippedSam.setCigarString("1S2M1S");

        final SAMRecord hardClippedSam = basicSAMRecord();
        hardClippedSam.setCigarString("3H2M2H");



        return new Object[][]{
                { new SAMRecordToGATKReadAdapter(softClippedSam), BASIC_READ_START - 1, BASIC_READ_START + 2 },
                { new SAMRecordToGATKReadAdapter(hardClippedSam), BASIC_READ_START - 3, BASIC_READ_START + 3 },
        };
    }

    @Test(dataProvider = "GetUnclippedStartAndEndData")
    public void testGetAndSetUnclippedStartAndEnd( final GATKRead read, final int expectedUnclippedStart, final int expectedUnclippedEnd ) {
        Assert.assertEquals(read.getUnclippedStart(), expectedUnclippedStart, "Wrong unclipped start for read");
        Assert.assertEquals(read.getUnclippedEnd(), expectedUnclippedEnd, "Wrong unclipped end for read");
    }

    @Test
    public void testGetSoftStartAndEnd() {
        final GATKRead read = basicReadBackedBySam();

        Assert.assertEquals(read.getSoftStart(), 5);
        Assert.assertEquals(read.getSoftEnd(), 7);

        // Should invalidate the cached soft start/end
        read.setCigar("1S2M1S1H");

        Assert.assertEquals(read.getSoftStart(), 4);
        Assert.assertEquals(read.getSoftEnd(), 7);

        // Should invalidate the cached soft start/end
        read.setCigar("4M");

        Assert.assertEquals(read.getSoftStart(), 5);
        Assert.assertEquals(read.getSoftEnd(), 8);
    }

    @Test
    public void testGetAdaptorBoundary() {
        final GATKRead read = basicReadBackedBySam();
        read.setFragmentLength(2);
        read.setIsReverseStrand(false);
        read.setMateIsReverseStrand(true);

        Assert.assertEquals(read.getAdaptorBoundary(), ReadUtils.getAdaptorBoundary(read));

        // Should invalidate the cached adaptor boundary
        read.setIsReverseStrand(true);

        Assert.assertEquals(read.getAdaptorBoundary(), ReadUtils.getAdaptorBoundary(read));
    }

    @DataProvider(name = "GetAndSetMatePositionData")
    public Object[][] getAndSetMatePositionData() {
        final SAMRecord samWithUnmappedMate = basicSAMRecord();
        samWithUnmappedMate.setMateUnmappedFlag(true);

        return new Object[][]{
                { basicReadBackedBySam(), BASIC_READ_MATE_CONTIG, BASIC_READ_MATE_START },
                { new SAMRecordToGATKReadAdapter(samWithUnmappedMate), null, ReadConstants.UNSET_POSITION },
        };
    }

    @Test(dataProvider = "GetAndSetMatePositionData")
    public void testGetAndSetMatePosition( final GATKRead read, final String expectedMateContig, final int expectedMateStart ) {
        Assert.assertEquals(read.getMateContig(), expectedMateContig, "Wrong contig for mate");
        Assert.assertEquals(read.getMateStart(), expectedMateStart, "Wrong start position for mate");

        read.setMatePosition("2", 52);
        Assert.assertFalse(read.mateIsUnmapped(), "Mate should be mapped after setMatePosition()");
        Assert.assertEquals(read.getMateContig(), "2", "Wrong contig for mate after setMatePosition()");
        Assert.assertEquals(read.getMateStart(), 52, "Wrong start position for mate after setMatePosition()");

        read.setMatePosition(new SimpleInterval("1", 40, 40));
        Assert.assertEquals(read.getMateContig(), "1", "Wrong contig for mate after setMatePosition()");
        Assert.assertEquals(read.getMateStart(), 40, "Wrong start position for mate after setMatePosition()");

        // Setting mate position should have the additional effect of marking the read as paired
        read.setIsPaired(false);
        read.setMatePosition("1", 1);
        Assert.assertTrue(read.isPaired(), "Read should be marked as paired after setMatePosition()");

        read.setIsPaired(false);
        read.setMatePosition(new SimpleInterval("2", 1, 1));
        Assert.assertTrue(read.isPaired(), "Read should be marked as paired after setMatePosition()");
    }

    @Test(dataProvider = "InvalidSetPositionData", expectedExceptions = IllegalArgumentException.class)
    public void testSetInvalidMatePositionOnSamBackedRead( final String contig, final int start ) {
        final GATKRead read = new SAMRecordToGATKReadAdapter(basicSAMRecord());
        read.setMatePosition(contig, start);
    }

    @DataProvider(name = "GetAndSetFragmentLengthData")
    public Object[][] getAndSetFragmentLengthData() {
        final SAMRecord samWithISize = basicSAMRecord();
        samWithISize.setInferredInsertSize(120);

        return new Object[][]{
                { basicReadBackedBySam(), 0 },
                { new SAMRecordToGATKReadAdapter(samWithISize), 120 },
        };
    }

    @Test(dataProvider = "GetAndSetFragmentLengthData")
    public void testGetAndSetFragmentLength( final GATKRead read, final int expectedFragmentLength ) {
        Assert.assertEquals(read.getFragmentLength(), expectedFragmentLength, "Wrong fragment length for read");

        read.setFragmentLength(50);
        Assert.assertEquals(read.getFragmentLength(), 50, "Wrong fragment length for read after setFragmentLength()");

        // Negative fragment lengths explicitly allowed
        read.setFragmentLength(-50);
        Assert.assertEquals(read.getFragmentLength(), -50, "Wrong fragment length for read after setFragmentLength()");
    }

    @DataProvider(name = "GetAndSetMappingQualityData")
    public Object[][] getAndSetMappingQualityData() {
        final SAMRecord samWithMappingQualityZero = basicSAMRecord();
        samWithMappingQualityZero.setMappingQuality(0);


        return new Object[][]{
                { basicReadBackedBySam(), BASIC_READ_MAPPING_QUALITY },
                { new SAMRecordToGATKReadAdapter(samWithMappingQualityZero), 0 },
        };
    }

    @Test(dataProvider = "GetAndSetMappingQualityData")
    public void testGetAndSetMappingQuality( final GATKRead read, final int expectedMappingQuality ) {
        Assert.assertEquals(read.getMappingQuality(), expectedMappingQuality, "Wrong mapping quality for read");

        read.setMappingQuality(50);
        Assert.assertEquals(read.getMappingQuality(), 50, "Wrong mapping quality for read after setMappingQuality()");
    }

    @DataProvider(name = "InvalidMappingQualityData")
    public Object[][] invalidMappingQualityData() {
        return new Object[][]{
                { -1 },
                { -10 },
                { 256 },
                { 1000 }
        };
    }

    @Test(dataProvider = "InvalidMappingQualityData", expectedExceptions = IllegalArgumentException.class)
    public void testSetInvalidMappingQualityOnSamBackedRead( final int invalidMappingQuality ) {
        final GATKRead read = basicReadBackedBySam();
        read.setMappingQuality(invalidMappingQuality);
    }

    @DataProvider(name = "GetAndSetBasesData")
    public Object[][] getAndSetBasesData() {
        final SAMRecord baselessSam = basicSAMRecord();
        baselessSam.setReadBases(null);

        final SAMRecord noAlignedSequenceSam = basicSAMRecord();
        noAlignedSequenceSam.setReadBases(SAMRecord.NULL_SEQUENCE);


        return new Object[][]{
                { basicReadBackedBySam(), BASIC_READ_BASES, "ACGT" },
                { new SAMRecordToGATKReadAdapter(baselessSam), new byte[0], "*" },
                { new SAMRecordToGATKReadAdapter(noAlignedSequenceSam), new byte[0], "*" },
        };
    }

    @Test(dataProvider = "GetAndSetBasesData")
    public void testGetAndSetBases( final GATKRead read, final byte[] expectedBases, final String expectedBasesString ) {
        Assert.assertEquals(read.getBases(), expectedBases, "Wrong bases for read");
        Assert.assertEquals(read.getBasesString(), expectedBasesString, "Wrong base string for read");

        final byte[] newBases = {'G', 'C', 'G', 'G'};
        read.setBases(newBases);
        Assert.assertEquals(read.getBases(), newBases, "Wrong bases for read after setBases()");
        Assert.assertEquals(read.getBasesString(), "GCGG", "Wrong base string for read after setBases()");
        for (int i = 0; i < newBases.length; i++) {
            Assert.assertEquals(read.getBase(i), newBases[i], "Wrong base string for read after setBases()");
        }
    }

    @DataProvider(name = "GetAndSetBaseQualitiesData")
    public Object[][] getAndSetBaseQualitiesData() {
        final SAMRecord noQualsSam = basicSAMRecord();
        noQualsSam.setBaseQualities(null);

        final SAMRecord emptyQualsSam = basicSAMRecord();
        emptyQualsSam.setBaseQualities(new byte[0]);

        return new Object[][]{
                { basicReadBackedBySam(), BASIC_READ_BASE_QUALITIES },
                { new SAMRecordToGATKReadAdapter(noQualsSam), new byte[0] },
                { new SAMRecordToGATKReadAdapter(emptyQualsSam), new byte[0] },
        };
    }

    @Test(dataProvider = "GetAndSetBaseQualitiesData")
    public void testGetAndSetBaseQualities( final GATKRead read, final byte[] expectedQuals ) {
        Assert.assertEquals(read.getBaseQualities(), expectedQuals, "Wrong base qualities for read");
        Assert.assertEquals(read.getBaseQualityCount(), expectedQuals.length, "Wrong number of base qualities for read");

        final byte[] newQuals = {1, 2, 3, 4};
        read.setBaseQualities(newQuals);
        Assert.assertEquals(read.getBaseQualities(), newQuals, "Wrong base qualities for read after setBaseQualities()");
        Assert.assertEquals(read.getBaseQualityCount(), newQuals.length, "Wrong number of base qualities for read after setBaseQualities()");
        for (int i = 0; i < newQuals.length; i++) {
            Assert.assertEquals(read.getBaseQuality(i), newQuals[i], "Wrong base quality for read after setBaseQualities()");
        }
    }

    @DataProvider(name = "GetAndSetCigarData")
    public Object[][] getAndSetCigarData() {

        SAMRecord noCigarSam = basicSAMRecord();
        noCigarSam.setCigar(null);

        SAMRecord emptyCigarSam = basicSAMRecord();
        emptyCigarSam.setCigar(new Cigar());

        return new Object[][]{
            { basicReadBackedBySam(), TextCigarCodec.decode(BASIC_READ_CIGAR) },
            { new SAMRecordToGATKReadAdapter(noCigarSam), new Cigar() },
            { new SAMRecordToGATKReadAdapter(emptyCigarSam), new Cigar() },
        };
    }

    @Test(dataProvider = "GetAndSetCigarData")
    public void testGetAndSetCigar( final GATKRead read, final Cigar expectedCigar ) {
        Assert.assertEquals(read.getCigar(), expectedCigar, "Wrong cigar for read");
        Assert.assertEquals(read.numCigarElements(), expectedCigar.numCigarElements(), "Wrong numCigarElements for read");

        final Cigar newCigar = TextCigarCodec.decode("4M");
        read.setCigar(newCigar);
        Assert.assertEquals(read.getCigar(), newCigar, "Wrong cigar for read after setCigar()");
        Assert.assertEquals(read.numCigarElements(), newCigar.numCigarElements(), "Wrong numCigarElements for read");
        Assert.assertEquals(read.getCigarElement(0), newCigar.getCigarElement(0), "Wrong numCigarElement for read");

        read.setCigar("2M2I");
        Assert.assertEquals(read.getCigar(), TextCigarCodec.decode("2M2I"), "Wrong cigar for read after setCigar()");

        read.setCigar(new Cigar());
        Assert.assertEquals(read.getCigar(), new Cigar(), "Wrong cigar for read after setCigar()");

        read.setCigar((Cigar)null);
        Assert.assertEquals(read.getCigar(), new Cigar(), "Wrong cigar for read after setCigar()");
    }

    @Test
    public void testCigarLengthCaching() {
        final GATKRead read = basicReadBackedBySam();

        Assert.assertEquals(read.numCigarElements(), 3);

        // Should invalidate the cached cigar length
        read.setCigar("4M");

        Assert.assertEquals(read.numCigarElements(), 1);
    }

    @DataProvider(name = "GetAndSetReadGroupData")
    public Object[][] getAndSetReadGroupData() {
        SAMRecord noRGSam = basicSAMRecord();
        noRGSam.clearAttributes();

        return new Object[][] {
            { basicReadBackedBySam(), BASIC_READ_GROUP },
            { new SAMRecordToGATKReadAdapter(noRGSam), null },
        };
    }

    @Test(dataProvider = "GetAndSetReadGroupData")
    public void testGetAndSetReadGroup( final GATKRead read, final String expectedReadGroup ) {
        Assert.assertEquals(read.getReadGroup(), expectedReadGroup, "Wrong read group for read");

        read.setReadGroup("NewReadGroup");
        Assert.assertEquals(read.getReadGroup(), "NewReadGroup", "Wrong read group for read after setReadGroup()");

        read.setReadGroup(null);
        Assert.assertEquals(read.getReadGroup(), null, "Read group should be null after setReadGroup(null)");
    }

    @DataProvider(name = "IsPairedData")
    public Object[][] isPairedData() {
        SAMRecord unpairedSAM = basicSAMRecord();
        unpairedSAM.setReadPairedFlag(false);

        SAMRecord properlyPairedSAM = basicSAMRecord();
        properlyPairedSAM.setProperPairFlag(true);

        SAMRecord unpairedProperlyPairedSAM = basicSAMRecord();
        unpairedProperlyPairedSAM.setReadPairedFlag(false);
        unpairedProperlyPairedSAM.setProperPairFlag(true);


        return new Object[][] {
                { basicReadBackedBySam(), true, false },
                { new SAMRecordToGATKReadAdapter(unpairedSAM), false, false },
                { new SAMRecordToGATKReadAdapter(properlyPairedSAM), true, true },

                // We only consider reads to be properly paired if they are also marked as paired
                { new SAMRecordToGATKReadAdapter(unpairedProperlyPairedSAM), false, false },
        };
    }

    @Test(dataProvider = "IsPairedData")
    public void testIsPaired( final GATKRead read, final boolean expectedIsPaired, final boolean expectedIsProperlyPaired ) {
        Assert.assertEquals(read.isPaired(), expectedIsPaired, "Read paired status incorrect");
        Assert.assertEquals(read.isProperlyPaired(), expectedIsProperlyPaired, "Read properly paired status incorrect");

        read.setIsPaired(false);
        Assert.assertFalse(read.isPaired(), "Read should be unpaired after setIsPaired(false)");
        Assert.assertFalse(read.isProperlyPaired(), "Read should not be marked as properly paired after setIsPaired(false)");

        read.setIsPaired(true);
        Assert.assertTrue(read.isPaired(), "Read should be paired after setIsPaired(true)");
        Assert.assertFalse(read.isProperlyPaired(), "Read incorrectly marked as properly paired");

        read.setIsProperlyPaired(true);
        Assert.assertTrue(read.isProperlyPaired(), "Read should be marked as properly paired after setIsProperlyPaired(true)");
        Assert.assertTrue(read.isPaired(), "Read should still be paired after setIsProperlyPaired(true)");

        read.setIsProperlyPaired(false);
        Assert.assertFalse(read.isProperlyPaired(), "Read should not be marked as properly paired after setIsProperlyPaired(false)");
        Assert.assertTrue(read.isPaired(), "Read should still be paired after setIsProperlyPaired(false)");

        read.setIsPaired(false);
        read.setIsProperlyPaired(true);
        Assert.assertTrue(read.isPaired(), "setIsProperlyPaired(true) should have also marked the read as paired");
    }

    @DataProvider(name = "IsUnmappedData")
    public Object[][] isUnmappedData() {
        List<Object[]> testCases = new ArrayList<>();

        for ( final GATKRead unmappedRead : getUnmappedReads() ) {
            testCases.add(new Object[]{ unmappedRead, true });
        }

        testCases.add(new Object[]{ basicReadBackedBySam(), false });

        return testCases.toArray(new Object[][]{});
    }

    @Test(dataProvider = "IsUnmappedData")
    public void testIsUnmapped( final GATKRead read, final boolean expectedIsUnmapped ) {
        Assert.assertEquals(read.isUnmapped(), expectedIsUnmapped, "Read mapped/unmapped status reported by isUnmapped() is incorrect");

        read.setIsUnmapped();
        Assert.assertTrue(read.isUnmapped(), "Read should be unmapped after setIsUnmapped()");

        read.setPosition("1", 1);
        Assert.assertFalse(read.isUnmapped(), "Read should be mapped after setPosition()");

        read.setIsUnmapped();
        Assert.assertTrue(read.isUnmapped(), "Read should be unmapped after setIsUnmapped()");
    }

    @DataProvider(name = "MateIsUnmappedData")
    public Object[][] mateIsUnmappedData() {
        SAMRecord samWithUnmappedMate = basicSAMRecord();
        samWithUnmappedMate.setMateUnmappedFlag(true);

        SAMRecord samWithUnmappedMate2 = basicSAMRecord();
        samWithUnmappedMate2.setMateReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);

        SAMRecord samWithUnmappedMate3 = basicSAMRecord();
        samWithUnmappedMate3.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);

        return new Object[][] {
                { basicReadBackedBySam(), false },
                { new SAMRecordToGATKReadAdapter(samWithUnmappedMate), true },
                { new SAMRecordToGATKReadAdapter(samWithUnmappedMate2), true },
                { new SAMRecordToGATKReadAdapter(samWithUnmappedMate3), true },
         };
    }

    @Test(dataProvider = "MateIsUnmappedData")
    public void testMateIsUnmapped( final GATKRead read, final boolean expectedMateIsUnmapped ) {
        Assert.assertEquals(read.mateIsUnmapped(), expectedMateIsUnmapped, "Mate unmapped status incorrect");

        read.setMateIsUnmapped();
        Assert.assertTrue(read.mateIsUnmapped(), "Mate should be unmapped after setMateIsUnmapped()");

        read.setMatePosition("1", 1);
        Assert.assertFalse(read.mateIsUnmapped(), "Mate should be mapped after setMatePosition()");

        read.setMateIsUnmapped();
        Assert.assertTrue(read.mateIsUnmapped(), "Mate should be unmapped after setMateIsUnmapped()");

        // Calling setMateIsUnmapped() should have the side effect of marking the read as paired
        read.setIsPaired(false);
        read.setMateIsUnmapped();
        Assert.assertTrue(read.isPaired(), "Read should be marked as paired after call to setMateIsUnmapped()");
    }

    @DataProvider(name = "InvalidMateIsUnmappedData")
    public Object[][] invalidMateIsUnmappedData() {
        SAMRecord unpairedSAM = basicSAMRecord();
        unpairedSAM.setReadPairedFlag(false);

        return new Object[][] {
                { new SAMRecordToGATKReadAdapter(unpairedSAM) },
        };
    }

    @Test(dataProvider = "InvalidMateIsUnmappedData", expectedExceptions = IllegalStateException.class)
    public void testInvalidMateIsUnmapped( final GATKRead read ) {
        // Calling mateIsUnmapped() on an unpaired read should throw an error
        read.mateIsUnmapped();
    }

    @DataProvider(name = "IsReverseStrandData")
    public Object[][] isReverseStrandData() {
        SAMRecord reverseStrandSam = basicSAMRecord();
        reverseStrandSam.setReadNegativeStrandFlag(true);

        return new Object[][] {
                { basicReadBackedBySam(), false },
                { new SAMRecordToGATKReadAdapter(reverseStrandSam), true },
        };
    }

    @Test(dataProvider = "IsReverseStrandData")
    public void testIsReverseStrand( final GATKRead read, final boolean expectedIsReverseStrand ) {
        Assert.assertEquals(read.isReverseStrand(), expectedIsReverseStrand, "Read strand incorrect");

        read.setIsReverseStrand(true);
        Assert.assertTrue(read.isReverseStrand(), "Read should be marked as reverse strand after setIsReverseStrand(true)");

        read.setIsReverseStrand(false);
        Assert.assertFalse(read.isReverseStrand(), "Read should not be marked as reverse strand after setIsReverseStrand(false)");
    }

    @DataProvider(name = "MateIsReverseStrandData")
    public Object[][] mateIsReverseStrandData() {
        SAMRecord samWithReverseStrandMate = basicSAMRecord();
        samWithReverseStrandMate.setMateNegativeStrandFlag(true);

        return new Object[][] {
                { basicReadBackedBySam(), false },
                { new SAMRecordToGATKReadAdapter(samWithReverseStrandMate), true },
        };
    }

    @Test(dataProvider = "MateIsReverseStrandData")
    public void testMateIsReverseStrand( final GATKRead read, final boolean expectedMateIsReverseStrand ) {
        Assert.assertEquals(read.mateIsReverseStrand(), expectedMateIsReverseStrand, "Mate strandedness incorrect");

        read.setMateIsReverseStrand(true);
        Assert.assertTrue(read.mateIsReverseStrand(), "Mate should be reverse strand after setMateIsReverseStrand(true)");

        read.setMateIsReverseStrand(false);
        Assert.assertFalse(read.mateIsReverseStrand(), "Mate should not be reverse strand after setMateIsReverseStrand(false)");

        // Calling setMateIsReverseStrand() should have the side effect of marking the read as paired.
        read.setIsPaired(false);
        read.setMateIsReverseStrand(true);
        Assert.assertTrue(read.isPaired(), "Read should be marked as paired after call to setMateIsReverseStrand()");
    }

    @DataProvider(name = "InvalidMateIsReverseStrandData")
    public Object[][] invalidMateIsReverseStrandData() {
        SAMRecord unpairedSAM = basicSAMRecord();
        unpairedSAM.setReadPairedFlag(false);

        return new Object[][] {
                { new SAMRecordToGATKReadAdapter(unpairedSAM) },
        };
    }

    @Test(dataProvider = "InvalidMateIsReverseStrandData", expectedExceptions = IllegalStateException.class)
    public void testInvalidMateIsReverseStrand( final GATKRead read ) {
        // Should throw when called on an unpaired read
        read.mateIsReverseStrand();
    }

    @DataProvider(name = "ReadNumberTestData")
    public Object[][] readNumberTestData() {
        SAMRecord secondOfPairSam = basicSAMRecord();
        secondOfPairSam.setSecondOfPairFlag(true);
        secondOfPairSam.setFirstOfPairFlag(false);

        return new Object[][] {
                { basicReadBackedBySam(), true, false },
                { new SAMRecordToGATKReadAdapter(secondOfPairSam), false, true },
        };
    }

    @Test(dataProvider = "ReadNumberTestData")
    public void testReadNumber( final GATKRead read, final boolean expectedIsFirstOfPair, final boolean expectedIsSecondOfPair ) {
        Assert.assertEquals(read.isFirstOfPair(), expectedIsFirstOfPair, "Read first of pair setting incorrect");
        Assert.assertEquals(read.isSecondOfPair(), expectedIsSecondOfPair, "Read second of pair setting incorrect");

        read.setIsFirstOfPair();
        Assert.assertTrue(read.isFirstOfPair(), "Read should be marked first of pair after setIsFirstOfPair()");
        Assert.assertFalse(read.isSecondOfPair(), "Read should not be marked second of pair after setIsFirstOfPair()");

        read.setIsSecondOfPair();
        Assert.assertFalse(read.isFirstOfPair(), "Read should not be marked first of pair after setIsSecondOfPair()");
        Assert.assertTrue(read.isSecondOfPair(), "Read should be marked second of pair after setIsSecondOfPair()");
    }

    @DataProvider(name = "GetReadNumberOfUnpairedReadData")
    public Object[][] getReadNumberOfUnpairedReadData() {
        SAMRecord unpairedSAM = basicSAMRecord();
        unpairedSAM.setFirstOfPairFlag(true);
        unpairedSAM.setSecondOfPairFlag(false);
        unpairedSAM.setReadPairedFlag(false);

        SAMRecord unpairedSAM2 = basicSAMRecord();
        unpairedSAM2.setSecondOfPairFlag(true);
        unpairedSAM2.setFirstOfPairFlag(false);
        unpairedSAM2.setReadPairedFlag(false);

        return new Object[][] {
                { new SAMRecordToGATKReadAdapter(unpairedSAM) },
                { new SAMRecordToGATKReadAdapter(unpairedSAM2) },
        };
    }

    @Test(dataProvider = "GetReadNumberOfUnpairedReadData")
    public void testGetReadNumberOfUnpairedRead( final GATKRead unpairedRead ) {
        Assert.assertFalse(unpairedRead.isFirstOfPair(), "Unpaired read should not be marked as first of pair");
        Assert.assertFalse(unpairedRead.isSecondOfPair(), "Unpaired read should not be marked as second of pair");
    }

    @DataProvider(name = "GetAndSetSimpleFlagsData")
    public Object[][] getAndSetSimpleFlagsData() {
        return new Object[][] {
                { basicReadBackedBySam() },
        };
    }

    @Test(dataProvider = "GetAndSetSimpleFlagsData")
    public void testGetAndSetSimpleFlags( final GATKRead read ) {
        read.setIsSecondaryAlignment(true);
        Assert.assertTrue(read.isSecondaryAlignment(), "Read should be marked as secondary alignment after setIsSecondaryAlignment(true)");

        read.setIsSecondaryAlignment(false);
        Assert.assertFalse(read.isSecondaryAlignment(), "Read should not be marked as secondary alignment after setIsSecondaryAlignment(false)");

        read.setIsSupplementaryAlignment(true);
        Assert.assertTrue(read.isSupplementaryAlignment(), "Read should be marked as supplementary alignment after setIsSupplementaryAlignment(true)");

        read.setIsSupplementaryAlignment(false);
        Assert.assertFalse(read.isSupplementaryAlignment(), "Read should not be marked as supplementary alignment after setIsSupplementaryAlignment(false)");

        read.setFailsVendorQualityCheck(true);
        Assert.assertTrue(read.failsVendorQualityCheck(), "Read should be marked as failing vendor quality checks after setFailsVendorQualityCheck(true)");

        read.setFailsVendorQualityCheck(false);
        Assert.assertFalse(read.failsVendorQualityCheck(), "Read should not be marked as failing vendor quality checks after setFailsVendorQualityCheck(false)");

        read.setIsDuplicate(true);
        Assert.assertTrue(read.isDuplicate(), "Read should be marked as a duplicate after setIsDuplicate(true)");

        read.setIsDuplicate(false);
        Assert.assertFalse(read.isDuplicate(), "Read should not be marked as a duplicate after setIsDuplicate(false)");
    }

    @DataProvider(name = "GetAndSetAttributesData")
    public Object[][] getAndSetAttributesData() {
        return new Object[][] {
                { basicReadBackedBySam() },
        };
    }

    @Test(dataProvider = "GetAndSetAttributesData")
    public void testGetAndSetIntegerAttribute( final GATKRead read ) {
        Assert.assertNull(read.getAttributeAsInteger("DR"), "Attribute DR should be null");

        read.setAttribute("DR", 5);
        Assert.assertEquals(read.getAttributeAsInteger("DR").intValue(), 5, "Wrong value for attribute DR");

        // test type coercion
        read.setAttribute("DR", "6");
        Assert.assertEquals(read.getAttributeAsInteger("DR").intValue(), 6, "Wrong value for attribute DR");

        read.setAttribute("DR", 10);
        Assert.assertEquals(read.getAttributeAsInteger("DR").intValue(), 10, "Wrong value for attribute DR");

        read.clearAttribute("DR");
        Assert.assertNull(read.getAttributeAsInteger("DR"), "Attribute DR should be null");
    }

    @Test(dataProvider = "GetAndSetAttributesData", expectedExceptions = GATKException.ReadAttributeTypeMismatch.class)
    public void testGetNonIntegerAttributeAsInteger( final GATKRead read ) {
        read.setAttribute("DR", "notaninteger");
        read.getAttributeAsInteger("DR");
    }

    @Test(dataProvider = "GetAndSetAttributesData")
    public void testGetAndSetStringAttribute( final GATKRead read ) {
        Assert.assertNull(read.getAttributeAsString("DR"), "Attribute DR should be null");

        read.setAttribute("DR", "foo");
        Assert.assertEquals(read.getAttributeAsString("DR"), "foo", "Wrong value for attribute DR");

        // Test type coercion
        read.setAttribute("DR", 5);
        Assert.assertEquals(read.getAttributeAsString("DR"), "5", "Wrong value for attribute DR");

        read.setAttribute("DR", "bar");
        Assert.assertEquals(read.getAttributeAsString("DR"), "bar", "Wrong value for attribute DR");

        read.setAttribute("DR", new byte[] {'A', 'C', 'T', 'G'});
        Assert.assertEquals(read.getAttributeAsString("DR"), "ACTG", "Wrong value for attribute DR");

        read.clearAttribute("DR");
        Assert.assertNull(read.getAttributeAsString("DR"), "Attribute DR should be null");
    }

    @Test(dataProvider = "GetAndSetAttributesData")
    public void testGetAndSetByteArrayAttribute( final GATKRead read ) {
        Assert.assertNull(read.getAttributeAsByteArray("DR"), "Attribute DR should be null");

        read.setAttribute("DR", new byte[]{1, 2, 3});
        Assert.assertEquals(read.getAttributeAsByteArray("DR"), new byte[]{1, 2, 3}, "Wrong value for attribute DR");

        // Test type coercion
        read.setAttribute("DR", "abc");
        Assert.assertEquals(read.getAttributeAsByteArray("DR"), new byte[]{'a', 'b', 'c'}, "Wrong value for attribute DR");

        read.setAttribute("DR", new byte[]{4, 5, 6, 7});
        Assert.assertEquals(read.getAttributeAsByteArray("DR"), new byte[]{4, 5, 6, 7}, "Wrong value for attribute DR");

        read.clearAttribute("DR");
        Assert.assertNull(read.getAttributeAsByteArray("DR"), "Attribute DR should be null");
    }

    @Test(expectedExceptions = GATKException.ReadAttributeTypeMismatch.class)
    public void testGetNonByteArrayAttributeAsByteArray() {
        // Only SAMRecord-backed reads can encounter a type mismatch for byte array attributes
        SAMRecord sam = basicSAMRecord();
        sam.setAttribute("DR", 5);
        final GATKRead samBackedRead = new SAMRecordToGATKReadAdapter(sam);

        samBackedRead.getAttributeAsByteArray("DR");
    }

    @Test(dataProvider = "GetAndSetAttributesData")
    public void testClearAllAttributes( final GATKRead read ) {
        read.setAttribute("DR", "foo");
        read.setAttribute("LT", 5);

        Assert.assertNotNull(read.getAttributeAsString("DR"), "Attribute DR should not be null");
        Assert.assertNotNull(read.getAttributeAsInteger("LT"), "Attribute LT should not be null");

        read.clearAttributes();
        Assert.assertNull(read.getAttributeAsString("DR"), "Attribute DR should be null");
        Assert.assertNull(read.getAttributeAsInteger("LT"), "Attribute LT should be null");
    }

    @Test
    public void testBasicGetSAMString() {
        // 1. GATKRead backed by a SAM record
        final SAMRecord samRecord = basicSAMRecord();
        final String samRecordString = samRecord.getSAMString();
        final GATKRead samBackedRead = new SAMRecordToGATKReadAdapter(samRecord);
        Assert.assertEquals(samRecordString, samBackedRead.getSAMString(), "SAM-backed GATKRead string should match wrapped SAM record string");
    }

    @Test(dataProvider = "GetAndSetPositionData")
    public void testSAMStringPosition(final GATKRead read, final String expectedContig, final int expectedStart, final int expectedEnd) {
        // just calling setPosition leaves the read in a state where isReverseStrand results in a missing field
        // exception for "strand" when calling getSAMString(), so we need to set the reverse strand as well
        read.setPosition("2", 17);
        read.setIsReverseStrand(false);
        final SAMRecord samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: GetAndSetPositionData");
    }

    @Test(dataProvider = "GetAndSetNameData")
    public void testSAMStringName( final GATKRead read, final String expectedName ) {
        read.setName("NEWNAME");
        SAMRecord samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: GetAndSetNameData");

        read.setName(null);
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: GetAndSetNameData");
    }

    @Test(dataProvider = "GetLengthData")
    public void testSAMStringLength(final GATKRead read, final int expectedLength) {
        final SAMRecord samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: GetLengthData");
    }

    @Test(dataProvider = "GetUnclippedStartAndEndData")
    public void testSAMStringUnclippedStartAndEnd(final GATKRead read, final int expectedUnclippedStart, final int expectedUnclippedEnd) {
        final SAMRecord samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: GetAndSetMatePositionData");
    }

    @Test(dataProvider = "GetAndSetMatePositionData")
    public void testSAMStringMatePosition(final GATKRead read, final String expectedMateContig, final int expectedMateStart) {
        SAMRecord samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: GetAndSetMatePositionData");

        read.setMatePosition("2", 52);
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: GetAndSetMatePositionData");

        read.setMatePosition(new SimpleInterval("1", 40, 40));
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: GetAndSetMatePositionData");

        // Setting mate position should have the additional effect of marking the read as paired
        read.setIsPaired(false);
        read.setMatePosition("1", 1);
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: GetAndSetMatePositionData");

        read.setIsPaired(false);
        samRec = read.convertToSAMRecord(getSAMHeader());
        final int flagInt = parseSAMStringFlags(read.getSAMString());
        Assert.assertTrue(((flagInt & ReadUtils.SAM_READ_PAIRED_FLAG) != 0) == false);
    }

    @Test(dataProvider = "GetAndSetFragmentLengthData")
    public void testSAMStringFragmentLengthString(final GATKRead read, final int expectedFragmentLength) {
        SAMRecord samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: GetAndSetMappingQualityData");

        read.setFragmentLength(50);
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: GetAndSetMappingQualityData");

        // Negative fragment lengths explicitly allowed
        read.setFragmentLength(-50);
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: GetAndSetMappingQualityData");
    }

    @Test(dataProvider = "GetAndSetMappingQualityData")
    public void testSAMStringMappingQuality(final GATKRead read, final int expectedMappingQuality) {
        SAMRecord samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: GetAndSetMappingQualityData");

        read.setMappingQuality(50);
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: GetAndSetMappingQualityData");
    }

    @Test(dataProvider = "GetAndSetBasesData")
    public void testSAMStringBases(final GATKRead read, final byte[] expectedBases, final String expectedBasesString) {
        final byte[] newBases = {'G', 'C', 'G', 'G'};
        read.setBases(newBases);

        final SAMRecord samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: GetAndSetBasesData");
    }

    @Test(dataProvider = "GetAndSetBaseQualitiesData")
    public void testSAMStringBaseQualities(final GATKRead read, final byte[] expectedQuals) {
        SAMRecord samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: GetAndSetBaseQualitiesData");

        final byte[] newQuals = {1, 2, 3, 4};
        read.setBaseQualities(newQuals);
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: GetAndSetBaseQualitiesData");
    }

    @Test(dataProvider = "GetAndSetCigarData")
    public void testSAMStringCigar(final GATKRead read, final Cigar expectedCigar) {
        SAMRecord samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: GetAndSetCigarData");

        final Cigar newCigar = TextCigarCodec.decode("4M");
        read.setCigar(newCigar);
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: GetAndSetCigarData");

        read.setCigar("2M2I");
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: GetAndSetCigarData");

        read.setCigar(new Cigar());
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: GetAndSetCigarData");

        read.setCigar((Cigar) null);
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: GetAndSetCigarData");
    }

    @Test(dataProvider = "GetAndSetReadGroupData")
    public void testSAMStringReadGroup(final GATKRead read, final String expectedReadGroup) {
        SAMRecord samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: IsPairedData");

        read.setReadGroup("NewReadGroup");
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: IsPairedData");

        read.setReadGroup(null);
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: IsPairedData");
    }

    @Test
    public void testSAMStringUnmappedMateContig() {
        // test the special case of an unpaired read to make sure we get the proper mate contig name
        SAMFileHeader samHeader = getSAMHeader();
        SAMRecord samRec = ArtificialReadUtils.createArtificialSAMRecord(
                samHeader,
                BASIC_READ_NAME,
                samHeader.getSequenceIndex(BASIC_READ_CONTIG),
                BASIC_READ_START,
                BASIC_READ_BASES,
                BASIC_READ_BASE_QUALITIES,
                BASIC_READ_CIGAR
        );

        final GATKRead read = new SAMRecordToGATKReadAdapter(samRec);
        Assert.assertTrue(samRec.getSAMString().equals(read.getSAMString()));
    }

    @Test(dataProvider = "IsPairedData")
    public void testSAMStringIsPaired(final GATKRead read, final boolean expectedIsPaired, final boolean expectedIsProperlyPaired) {
        int flagInt = parseSAMStringFlags(read.getSAMString());
        Assert.assertTrue(((flagInt & ReadUtils.SAM_READ_PAIRED_FLAG) != 0) == expectedIsPaired);

        read.setIsPaired(false);
        flagInt = parseSAMStringFlags(read.getSAMString());
        Assert.assertTrue(((flagInt & ReadUtils.SAM_READ_PAIRED_FLAG) != 0) == false);

        read.setIsProperlyPaired(true);
        flagInt = parseSAMStringFlags(read.getSAMString());
        Assert.assertTrue(((flagInt & ReadUtils.SAM_PROPER_PAIR_FLAG) != 0) == true);

        read.setIsProperlyPaired(false);
        flagInt = parseSAMStringFlags(read.getSAMString());
        Assert.assertTrue(((flagInt & ReadUtils.SAM_PROPER_PAIR_FLAG) != 0) == false);
    }

    @Test(dataProvider = "IsUnmappedData")
    public void testSAMStringIsUnmapped(final GATKRead read, final boolean expectedIsUnmapped) {
        read.setIsUnmapped();
        SAMRecord samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: MateIsUnmappedData");

        read.setPosition("1", 1);
        read.setIsReverseStrand(false);
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: MateIsUnmappedData");

        read.setIsUnmapped();
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: MateIsUnmappedData");
    }

    @Test(dataProvider = "MateIsUnmappedData")
    public void testSAMStringMateIsUnmapped(final GATKRead read, final boolean expectedMateIsUnmapped) {
        read.setMatePosition("1", 1);
        final SAMRecord samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: MateIsUnmappedData");

        // Calling setMateIsUnmapped() should have the side effect of marking the read as paired
        read.setMateIsUnmapped();
        final int flagInt = parseSAMStringFlags(read.getSAMString());
        Assert.assertTrue(((flagInt & ReadUtils.SAM_MATE_UNMAPPED_FLAG) != 0) == true);
        Assert.assertTrue(((flagInt & ReadUtils.SAM_READ_PAIRED_FLAG )!= 0) == true);
    }

    @Test(dataProvider = "IsReverseStrandData")
    public void testSAMStringIsReverseStrand(final GATKRead read, final boolean expectedIsReverseStrand) {
        SAMRecord samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: IsReverseStrandData");

        read.setIsReverseStrand(true);
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: IsReverseStrandData");

        read.setIsReverseStrand(false);
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: IsReverseStrandData");
    }

    @Test(dataProvider = "MateIsReverseStrandData")
    public void testSAMStringMateIsReverseStrand(final GATKRead read, final boolean expectedMateIsReverseStrand) {
        SAMRecord samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: ReadNumberTestData");

        read.setMateIsReverseStrand(true);
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: ReadNumberTestData");

        read.setMateIsReverseStrand(false);
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: ReadNumberTestData");

        // Calling setMateIsReverseStrand() should have the side effect of marking the read as paired.
        read.setIsPaired(false);
        read.setMateIsReverseStrand(true);
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: ReadNumberTestData");
    }

    @Test(dataProvider = "ReadNumberTestData")
    public void testSAMStringReadNumber(final GATKRead read, final boolean expectedIsFirstOfPair, final boolean expectedIsSecondOfPair) {
        SAMRecord samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: ReadNumberTestData");

        read.setIsFirstOfPair();
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: ReadNumberTestData");

        read.setIsSecondOfPair();
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "SAM string comparison failure: ReadNumberTestData");
    }

    @Test(dataProvider = "GetAndSetAttributesData")
    public void testSAMStringAttributes( final GATKRead read ) {
        Assert.assertNull(read.getAttributeAsInteger("DR"), "Attribute DR should be null");

        read.setAttribute("DR", 5);
        String samString = read.getSAMString();

        // We lose type information for custom attributes for Google Genomics reads, but the data provider
        // for this test includes both SAM-backed and Google read-backed items so we need to test for either
        Assert.assertTrue(samString.contains("DR:i:5") || samString.contains("DR:Z:5"), "SAM string comparison failure: GetAndSetAttributesData");

        // test type coercion
        read.setAttribute("DR", "6");
        samString = read.getSAMString();
        Assert.assertTrue(samString.contains("DR:Z:6"), "SAM string comparison failure: GetAndSetAttributesData");

        read.clearAttribute("DR");
        samString = read.getSAMString();
        Assert.assertTrue(!samString.contains("DR:"), "SAM string comparison failure: GetAndSetAttributesData");

        // We lose type information for custom attributes for Google Genomics reads, but the data provider
        // for this test includes both SAM-backed and Google read-backed items so we need to test for either
        read.setAttribute("DR", new byte[]{1, 2, 3});
        samString = read.getSAMString();
        Assert.assertTrue(samString.contains("DR:"), "SAM string comparison failure: GetAndSetAttributesData");

        read.clearAttribute("DR");
        samString = read.getSAMString();
        Assert.assertTrue(!samString.contains("DR:"), "SAM string comparison failure: GetAndSetAttributesData");
    }

    @Test(dataProvider = "GetAndSetSimpleFlagsData")
    public void testSAMStringSimpleFlags(final GATKRead read) {
        SAMRecord samRec;

        read.setIsSecondaryAlignment(true);
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "Invalid SAM string after setIsSecondaryAlignment(true)");

        read.setIsSecondaryAlignment(false);
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "Invalid SAM string after setIsSecondaryAlignment(false)");

        read.setIsSupplementaryAlignment(true);
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "Invalid SAM string after setIsSupplementaryAlignment(true)");

        read.setIsSupplementaryAlignment(false);
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "Invalid SAM string after setIsSupplementaryAlignment(false)");

        read.setFailsVendorQualityCheck(true);
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "Invalid SAM string after setFailsVendorQualityCheck(true)");

        read.setFailsVendorQualityCheck(false);
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "Invalid SAM string after setFailsVendorQualityCheck(false)");

        read.setIsDuplicate(true);
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "Invalid SAM string after setIsDuplicate(true)");

        read.setIsDuplicate(false);
        samRec = read.convertToSAMRecord(getSAMHeader());
        Assert.assertEquals(read.getSAMString(), samRec.getSAMString(), "Invalid SAM string after setIsDuplicate(false)");
    }

    //  pull the flags field out of a string produced by getSAMString()
    private int parseSAMStringFlags(String samString) {
        final Pattern p = Pattern.compile("(\\t)(\\d+)"); // find the int field in the second column, which are the flags
        final  Matcher m = p.matcher(samString);
        Assert.assertTrue(m.find());
        return Integer.parseInt(m.group(2)); // use the second capture group in the reg ex above since we don't want the tab
    }

    @DataProvider(name="copyData")
    public Object[][] getCopyData() {
        List<Object[]> testCases = new ArrayList<>();

        testCases.add(new Object[]{basicReadBackedBySam()});

        return testCases.toArray(new Object[][]{});
    }

    @Test(dataProvider="copyData")
    public void testShallowCopy(final GATKRead read) {
        Assert.assertEquals(read, read.copy());

        read.setIsReverseStrand(false);
        final GATKRead shallowCopy = read.copy();
        Assert.assertEquals(read.isReverseStrand(), shallowCopy.isReverseStrand());
        read.setIsReverseStrand(true);
        Assert.assertNotEquals(read.isReverseStrand(), shallowCopy.isReverseStrand());
        Assert.assertNotEquals(read, shallowCopy);

        Assert.assertEquals(read, read.deepCopy());
    }

    @Test(dataProvider="copyData")
    public void testDeepCopyReverseStrand(final GATKRead read) {
        Assert.assertEquals(read, read.deepCopy());

        // make sure we've copied deeply; i.e. for GoogleRead backed reads,
        // reverseStrand is nested as LinearAlignment->Position->reverseStrand
        read.setIsReverseStrand(false);
        final GATKRead deepCopy = read.deepCopy();
        Assert.assertEquals(read.isReverseStrand(), deepCopy.isReverseStrand());
        read.setIsReverseStrand(true);
        Assert.assertNotEquals(read.isReverseStrand(), deepCopy.isReverseStrand());
        Assert.assertNotEquals(read, deepCopy);

        Assert.assertEquals(read, read.deepCopy());
    }

    @Test(dataProvider="copyData")
    public void testDeepCopyAttributes(final GATKRead read) {
        Assert.assertEquals(read, read.deepCopy());

        GATKRead deepCopy = read.deepCopy();
        byte attr[] = new byte[]{'B', 'I'};
        read.setAttribute("BI", attr);
        Assert.assertEquals(read.getAttributeAsByteArray("BI"), attr);
        Assert.assertNull(deepCopy.getAttributeAsByteArray("BI"));
        Assert.assertNotEquals(read, deepCopy);

        deepCopy = read.deepCopy();
        Assert.assertEquals(read.getAttributeAsByteArray("BI"), deepCopy.getAttributeAsByteArray("BI"));
    }

    @Test(dataProvider="copyData")
    public void testDeepCopyPosition(final GATKRead read) {
        Assert.assertEquals(read, read.deepCopy());

        // position
        final GATKRead deepCopy = read.deepCopy();
        read.setPosition("2", 1);
        Assert.assertNotEquals(read.getContig(), deepCopy.getContig());
        Assert.assertNotEquals(read, deepCopy);

        Assert.assertEquals(read, read.deepCopy());
    }

    @DataProvider
    public Object[][] unmappedReadDataProvider() {
        List<Object[]> testCases = new ArrayList<>();

        SAMRecord unmappedUnplacedRead = new SAMRecord(hg19Header);
        SAMRecord unmappedPlacedRead = new SAMRecord(hg19Header);
        SAMRecord mappedUnplacedRead = new SAMRecord(hg19Header);
        SAMRecord mappedPlacedRead = new SAMRecord(hg19Header);

        unmappedUnplacedRead.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
        unmappedUnplacedRead.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
        unmappedUnplacedRead.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
        unmappedUnplacedRead.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
        unmappedUnplacedRead.setReadUnmappedFlag(true);
        unmappedUnplacedRead.setMateUnmappedFlag(true);
        unmappedUnplacedRead.setReadPairedFlag(true);
        testCases.add(new Object[]{new SAMRecordToGATKReadAdapter(unmappedUnplacedRead), true, true});


        unmappedPlacedRead.setReferenceIndex(0);
        unmappedPlacedRead.setMateReferenceIndex(0);
        unmappedPlacedRead.setAlignmentStart(100);
        unmappedPlacedRead.setMateAlignmentStart(100);
        unmappedPlacedRead.setReadUnmappedFlag(true);
        unmappedPlacedRead.setMateUnmappedFlag(true);
        unmappedPlacedRead.setReadPairedFlag(true);
        testCases.add(new Object[]{new SAMRecordToGATKReadAdapter(unmappedPlacedRead), true, false});


        mappedUnplacedRead.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
        mappedUnplacedRead.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
        mappedUnplacedRead.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
        mappedUnplacedRead.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
        mappedUnplacedRead.setReadUnmappedFlag(false);
        mappedUnplacedRead.setMateUnmappedFlag(false);
        mappedUnplacedRead.setReadPairedFlag(true);
        testCases.add(new Object[]{new SAMRecordToGATKReadAdapter(mappedUnplacedRead), true, true}); // Unmapped is true here because it tests alignment info second


        mappedPlacedRead.setReferenceIndex(0);
        mappedPlacedRead.setMateReferenceIndex(0);
        mappedPlacedRead.setAlignmentStart(100);
        mappedPlacedRead.setMateAlignmentStart(100);
        mappedPlacedRead.setReadUnmappedFlag(false);
        mappedPlacedRead.setMateUnmappedFlag(false);
        mappedPlacedRead.setReadPairedFlag(true);
        testCases.add(new Object[]{new SAMRecordToGATKReadAdapter(mappedPlacedRead), false, false});

        return testCases.toArray(new Object[][]{});
    }

    @Test(dataProvider="unmappedReadDataProvider")
    public void testUnmappedUnpairedDistinction(final GATKRead read, final boolean unmapped, final boolean unplaced) {
        Assert.assertEquals(read.isUnmapped(), unmapped);
        Assert.assertEquals(read.isUnplaced(), unplaced);

        read.setIsUnmapped();

        Assert.assertEquals(read.isUnmapped(), true);
        Assert.assertEquals(read.isUnplaced(), unplaced);

        read.setIsUnplaced();
        Assert.assertEquals(read.isUnmapped(), true);
        Assert.assertEquals(read.isUnplaced(), true);

        read.convertToSAMRecord(hg19Header).setReadUnmappedFlag(false); // have to unset unmapped to allow SAMRecord to show its underlying data
        Assert.assertEquals(read.convertToSAMRecord(hg19Header).getStart(), SAMRecord.NO_ALIGNMENT_START);
        Assert.assertEquals(read.convertToSAMRecord(hg19Header).getMappingQuality(), SAMRecord.NO_MAPPING_QUALITY);
    }

    @Test(dataProvider="unmappedReadDataProvider")
    public void testUnmappedUnpairedDistinctionmates(final GATKRead read, final boolean unmapped, final boolean unplaced) {
        Assert.assertEquals(read.mateIsUnmapped(), unmapped);
        Assert.assertEquals(read.mateIsUnplaced(), unplaced);

        read.setMateIsUnmapped();

        Assert.assertEquals(read.mateIsUnmapped(), true);
        Assert.assertEquals(read.mateIsUnplaced(), unplaced);

        read.setMateIsUnplaced();
        Assert.assertEquals(read.mateIsUnmapped(), true);
        Assert.assertEquals(read.mateIsUnplaced(), true);

        read.convertToSAMRecord(hg19Header).setReadUnmappedFlag(false); // have to unset unmapped to allow SAMRecord to show its underlying data
        Assert.assertEquals(read.convertToSAMRecord(hg19Header).getMateAlignmentStart(), SAMRecord.NO_ALIGNMENT_START);
    }

    @DataProvider(name = "pairOrientationEmptyResults")
    public Object[][] pairOrientationEmpty() {
        final GATKRead unmappedRead = ArtificialReadUtils.createArtificialRead("100M");
        unmappedRead.setIsPaired(true);
        unmappedRead.setIsUnmapped();

        final GATKRead mateUnmappedRead = ArtificialReadUtils.createArtificialRead("100M");
        mateUnmappedRead.setIsPaired(true);
        mateUnmappedRead.setIsUnmapped();

        final GATKRead mateOnDifferentContigRead = ArtificialReadUtils.createArtificialRead("100M");
        mateOnDifferentContigRead.setIsPaired(true);
        mateOnDifferentContigRead.setIsReverseStrand(false);
        mateOnDifferentContigRead.setMateIsReverseStrand(true);
        mateOnDifferentContigRead.setPosition(new SimpleInterval("1", 100, 200));
        mateOnDifferentContigRead.setMatePosition(new SimpleInterval("2", 200, 300));

        return new Object[][]{
                {unmappedRead},
                {mateUnmappedRead},
                {mateOnDifferentContigRead}
        };
    }

    @Test(dataProvider = "pairOrientationEmptyResults")
    public void testGetPairOrientationEmpty(final GATKRead read) {
        Assert.assertFalse(read.getPairOrientation().isPresent());
    }

    @DataProvider(name = "pairOrientationFailures")
    public Object[][] pairOrientationFailures() {
        final GATKRead unpairedRead = ArtificialReadUtils.createArtificialRead("100M");
        unpairedRead.setIsPaired(false);

        final GATKRead unmappedRead = ArtificialReadUtils.createArtificialRead("100M");
        unmappedRead.setIsPaired(true);
        unmappedRead.setIsUnmapped();

        final GATKRead mateUnmappedRead = ArtificialReadUtils.createArtificialRead("100M");
        mateUnmappedRead.setIsPaired(true);
        mateUnmappedRead.setIsUnmapped();

        return new Object[][]{
                {unpairedRead}
        };
    }

    @Test(dataProvider = "pairOrientationFailures", expectedExceptions = IllegalArgumentException.class)
    public void testGetPairOrientationFailures(final GATKRead read) {
        final Optional<SamPairUtil.PairOrientation> pairOrientation = read.getPairOrientation();
    }

    @DataProvider(name = "pairOrientation")
    public Object[][] pairOrientation() {
        final GATKRead forwardStrandInnie = ArtificialReadUtils.createArtificialRead("100M");
        forwardStrandInnie.setIsPaired(true);
        forwardStrandInnie.setIsReverseStrand(false);
        forwardStrandInnie.setMateIsReverseStrand(true);
        forwardStrandInnie.setPosition(new SimpleInterval("1", 100, 200));
        forwardStrandInnie.setMatePosition(new SimpleInterval("1", 200, 300));
        forwardStrandInnie.setFragmentLength(201);

        final GATKRead reverseStrandInnie = ArtificialReadUtils.createArtificialRead("100M");
        reverseStrandInnie.setIsPaired(true);
        reverseStrandInnie.setIsReverseStrand(true);
        reverseStrandInnie.setMateIsReverseStrand(false);
        reverseStrandInnie.setPosition(new SimpleInterval("1", 200, 300));
        reverseStrandInnie.setMatePosition(new SimpleInterval("1", 100, 200));
        reverseStrandInnie.setFragmentLength(201);

        final GATKRead forwardStrandOutie = ArtificialReadUtils.createArtificialRead("100M");
        forwardStrandOutie.setIsPaired(true);
        forwardStrandOutie.setIsReverseStrand(false);
        forwardStrandOutie.setMateIsReverseStrand(true);
        forwardStrandOutie.setPosition(new SimpleInterval("1", 200, 300));
        forwardStrandOutie.setMatePosition(new SimpleInterval("1", 100, 200));
        forwardStrandOutie.setFragmentLength(-201);

        final GATKRead reverseStrandOutie = ArtificialReadUtils.createArtificialRead("100M");
        reverseStrandOutie.setIsPaired(true);
        reverseStrandOutie.setIsReverseStrand(true);
        reverseStrandOutie.setMateIsReverseStrand(false);
        reverseStrandOutie.setPosition(new SimpleInterval("1", 100, 200));
        reverseStrandOutie.setMatePosition(new SimpleInterval("1", 200, 300));
        reverseStrandOutie.setFragmentLength(-201);

        final GATKRead forwardStrandTandem = ArtificialReadUtils.createArtificialRead("100M");
        forwardStrandTandem.setIsPaired(true);
        forwardStrandTandem.setIsReverseStrand(false);
        forwardStrandTandem.setMateIsReverseStrand(false);
        forwardStrandTandem.setPosition(new SimpleInterval("1", 100, 200));
        forwardStrandTandem.setMatePosition(new SimpleInterval("1", 200, 300));

        final GATKRead reverseStrandTandem = ArtificialReadUtils.createArtificialRead("100M");
        reverseStrandTandem.setIsPaired(true);
        reverseStrandTandem.setIsReverseStrand(true);
        reverseStrandTandem.setMateIsReverseStrand(true);
        reverseStrandTandem.setPosition(new SimpleInterval("1", 100, 200));
        reverseStrandTandem.setMatePosition(new SimpleInterval("1", 200, 300));

        return new Object[][]{
                {forwardStrandInnie, SamPairUtil.PairOrientation.FR},
                {reverseStrandInnie, SamPairUtil.PairOrientation.FR},
                {forwardStrandOutie, SamPairUtil.PairOrientation.RF},
                {reverseStrandOutie, SamPairUtil.PairOrientation.RF},
                {forwardStrandTandem, SamPairUtil.PairOrientation.TANDEM},
                {reverseStrandTandem, SamPairUtil.PairOrientation.TANDEM},
        };
    }

    @Test(dataProvider = "pairOrientation")
    public void testGetPairOrientation(final GATKRead read, final SamPairUtil.PairOrientation expected) {
        final SamPairUtil.PairOrientation pairOrientation = read.getPairOrientation().get();
        Assert.assertEquals(pairOrientation, expected);
    }

    private GATKRead makeFlowBasedRead(final byte[] bases, final boolean isReverse) {

        byte[] quals = new byte[bases.length];

        final String cigar = String.format("%dM", bases.length);
        GATKRead read = ArtificialReadUtils.createArtificialRead(bases, quals, cigar);

        final SAMFileHeader header = getSAMHeader();
        SAMReadGroupRecord record = header.getReadGroup(BASIC_READ_GROUP);
        record.setAttribute(SAMReadGroupRecord.FLOW_ORDER_TAG, FlowBasedRead.DEFAULT_FLOW_ORDER);
        ((SAMRecordToGATKReadAdapter)read).setHeader(header);
        read.setAttribute(SAMTag.RG.name(), BASIC_READ_GROUP);

        read.setPosition("chr1",100);
        read.setIsReverseStrand(isReverse);
        // convert to a flow based read
        ArtificialReadUtils.makeIntoFlowBased(read);
        return read;
    }

    @DataProvider(name="flowBasedReads")
    public Object[][] makeMultipleReads(){
        List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[]{makeFlowBasedRead(new byte[]{'T','A','G','C','G','A'}, false), 1,3,6});
        tests.add(new Object[]{makeFlowBasedRead(new byte[]{'A','C','C','G','A','T'}, false), 0,4,6});
        tests.add(new Object[]{makeFlowBasedRead(new byte[]{'T','A','G','C','G','A'}, true), 0,6,6});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "flowBasedReads")
    public void testHardClipAttributes(final GATKRead read, final int newStart, final int newLength, final int originalLength){

        //Add an attribute that is not suppsed to be clipped
        byte[]      ua = new byte[read.getBasesNoCopy().length];
        Arrays.fill(ua, (byte)1);
        String untrimmableAttribute = "ua";
        read.setAttribute(untrimmableAttribute, ua);

        read.hardClipAttributes(newStart, newLength, originalLength);
        for (final String tag : FlowBasedRead.HARD_CLIPPED_TAGS){
            if (read.hasAttribute(tag)){
                Object attribute = ((SAMRecordToGATKReadAdapter)read).getEncapsulatedSamRecord().getAttribute(tag);
                if (attribute instanceof String) {
                    Assert.assertEquals(((String)attribute).length(), newLength);
                } else {
                    Assert.assertEquals(read.getAttributeAsByteArray(tag).length, newLength);
                }
            }
        }

        Assert.assertEquals(read.getAttributeAsByteArray(untrimmableAttribute).length, originalLength);

    }

}
