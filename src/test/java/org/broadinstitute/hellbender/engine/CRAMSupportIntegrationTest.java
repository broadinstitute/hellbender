package org.broadinstitute.hellbender.engine;

import com.google.common.jimfs.Configuration;
import com.google.common.jimfs.Jimfs;
import htsjdk.samtools.*;
import htsjdk.samtools.cram.ref.ReferenceSource;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.PrintReads;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.FileSystem;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.function.Function;

public final class CRAMSupportIntegrationTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File("src/test/resources/org/broadinstitute/hellbender/engine");

    @Override
    public String getTestedClassName() {
        return PrintReads.class.getSimpleName();
    }

    @DataProvider(name="roundTripCRAMTests")
    public Object[][] getRoundTripCRAMTests() {
        return new Object[][] {
                // we need to use lenient equality for this test because this bam has ~14k reads that fail full
                // read equality; at least some of which are because they are unmapped/unplaced, but have cigar
                // strings that both samtools and htsjdk drop when roundtripping
                {NA12878_20_21_WGS_bam, b37_reference_20_21, true, false},
                // this cram is the result of converting the above bam to cram using samtools; once the file is
                // converted, we can use full read equality when roundtripping through cram, so we don't need to
                // be lenient
                {NA12878_20_21_WGS_cram, b37_reference_20_21, false, false},
                // roundtrip a v3.0 file
                { NA12878_20_21_WGS_cram, b37_reference_20_21, false, false },
                // roundtrip a v3.1 file
                { NA12878_20_21_WGS_cram_31, b37_reference_20_21, false, false },
        };
    }

    @Test(dataProvider="roundTripCRAMTests")
    public void testRoundTripToCRAM(
            final String sourceBamOrCRAM,
            final String reference,
            final boolean lenientEquality,
            final boolean emitDetail) throws IOException {
        final File outputCram = createTempFile("testRoundTripCRAM", ".cram");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(sourceBamOrCRAM);
        args.addOutput(outputCram);
        args.addReference(reference);
        args.add(ReadFilterArgumentDefinitions.DISABLE_TOOL_DEFAULT_READ_FILTERS, true);

        runCommandLine(args, "PrintReads");
        assertRoundTripFidelity(new File(sourceBamOrCRAM), outputCram, new File(reference), lenientEquality, emitDetail);
    }

    public void assertRoundTripFidelity(
            final File sourceFile,
            final File targetCRAMFile,
            final File referenceFile,
            final boolean lenientEquality,
            final boolean emitDetail) throws IOException {
        try (final SamReader sourceReader = SamReaderFactory.makeDefault()
                .referenceSequence(referenceFile)
                .validationStringency((ValidationStringency.SILENT))
                .open(sourceFile);
             final CRAMFileReader copyReader = new CRAMFileReader(targetCRAMFile, new ReferenceSource(referenceFile))) {
            final SAMRecordIterator sourceIterator = sourceReader.iterator();
            final SAMRecordIterator targetIterator = copyReader.getIterator();

            while (sourceIterator.hasNext() && targetIterator.hasNext()) {
                if (lenientEquality) {
                    final SAMRecord sourceRec = sourceIterator.next();
                    final SAMRecord targetRec = targetIterator.next();
                    Assert.assertEquals(targetRec.getReadName(), sourceRec.getReadName());
                    Assert.assertEquals(targetRec.getReadBases(), sourceRec.getReadBases());
                    Assert.assertEquals(targetRec.getBaseQualities(), sourceRec.getBaseQualities());
                    Assert.assertEquals(targetRec.getAlignmentStart(), sourceRec.getAlignmentStart());
                    Assert.assertEquals(targetRec.getAlignmentEnd(), sourceRec.getAlignmentEnd());
                } else if (emitDetail) {
                    final SAMRecord sourceRec = sourceIterator.next();
                    final SAMRecord targetRec = targetIterator.next();
                    if (!sourceRec.equals(targetRec)) {
                        System.out.println("Difference found:");
                        System.out.println(sourceRec.getSAMString());
                        System.out.println(targetRec.getSAMString());
                    }
                    Assert.assertEquals(targetRec, sourceRec);
                } else {
                    Assert.assertEquals(targetIterator.next(), sourceIterator.next());
                }
            }
            Assert.assertEquals(sourceIterator.hasNext(), targetIterator.hasNext());
        }
    }


    @DataProvider(name = "ReadEntireCramTestData")
    public Object[][] readEntireCramTestData() {
        final File ref = new File(hg19MiniReference);
        final List<Object[]> testCases = new ArrayList<>();
        for ( final String outputExtension : Arrays.asList(".sam", ".bam", ".cram") ) {
            // cram, reference, output file extension, expected read names
            testCases.add(new Object[]{new File(TEST_DATA_DIR, "cram_with_bai_index.cram"), ref, outputExtension, Arrays.asList("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k")});
            testCases.add(new Object[]{new File(TEST_DATA_DIR, "cram_with_crai_index.cram"), ref, outputExtension, Arrays.asList("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k")});
        }
        return testCases.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ReadEntireCramTestData")
    public void testReadEntireCram(final File cramFile, final File reference, final String outputExtension, final List<String> expectedReadNames ) throws IOException {
        final File outputFile = createTempFile("testReadEntireCram", outputExtension);
        final String[] args = new String[] {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, cramFile.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, reference.getAbsolutePath()
        };
        runCommandLine(args);

        SamAssertionUtils.assertCRAMContentsIfCRAM(outputFile);
        checkReadNames(outputFile, reference, expectedReadNames);
    }

    @Test
    public void testSamtoolsGeneratedCRAMSliceMD5Calculation() throws IOException {
        // Note: The input CRAM used for this test was generated using samtools, with the "ambiguityCodes.fasta" file
        // as the reference. Since the reference contains ambiguity codes (essential, since part of what we're trying
        // to validate is that htsjdk calculates the same MD5 value for a slice who's reference spans ambiguity codes
        // as samtools does), and since GATK wants an accompanying sequence dictionary, the .dict was generated
        // with GATK because samtools has a bug in dictionary generation when the reference has ambiguity codes.
        // See https://github.com/samtools/samtools/issues/704, and gatk tracking issue:
        // https://github.com/broadinstitute/gatk/issues/3306
        final File samtoolsGeneratedCRAM = new File(TEST_DATA_DIR, "samtoolsSliceMD5WithAmbiguityCodesTest.cram");
        final File referenceWithAmbiguityCodes = new File(TEST_DATA_DIR, "ambiguityCodes.fasta");
        final File outputFile = createTempFile("testReadSamtoolsGeneratedCRAM", ".cram");
        final String[] args = new String[] {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, samtoolsGeneratedCRAM.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, referenceWithAmbiguityCodes.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME
        };
        runCommandLine(args);  // no assert, just make sure we don't throw
    }


    @DataProvider(name = "ReadCramWithIntervalsIndexTestData")
    public Object[][] readCramWithIntervalsBAIIndexTestData() {
        final File ref = new File(hg19MiniReference);
        final List<Object[]> testCases = new ArrayList<>();
        for ( final String outputExtension : Arrays.asList(".sam", ".bam", ".cram") ) {
            // cram, reference, output file extension, intervals, expected read names
            testCases.add(new Object[]{ new File(TEST_DATA_DIR, "cram_with_bai_index.cram"), ref, outputExtension, Arrays.asList("1:1000-2000"), Arrays.asList("d", "e") });
            testCases.add(new Object[]{ new File(TEST_DATA_DIR, "cram_with_bai_index.cram"), ref, outputExtension, Arrays.asList("1:1000-1099"), Arrays.asList("d") });
            testCases.add(new Object[]{ new File(TEST_DATA_DIR, "cram_with_bai_index.cram"), ref, outputExtension, Arrays.asList("1:1-200", "1:1000-2000"), Arrays.asList("a", "d", "e") });
            testCases.add(new Object[]{ new File(TEST_DATA_DIR, "cram_with_bai_index.cram"), ref, outputExtension, Arrays.asList("1:1000-2000", "3:300-400"), Arrays.asList("d", "e", "i", "j") });
            testCases.add(new Object[]{ new File(TEST_DATA_DIR, "cram_with_bai_index.cram"), ref, outputExtension, Arrays.asList("1:1000-2000", "3:300-399"), Arrays.asList("d", "e", "i") });

            testCases.add(new Object[]{ new File(TEST_DATA_DIR, "cram_with_crai_index.cram"), ref, outputExtension, Arrays.asList("1:1000-2000"), Arrays.asList("d", "e") });
            testCases.add(new Object[]{ new File(TEST_DATA_DIR, "cram_with_crai_index.cram"), ref, outputExtension, Arrays.asList("1:1000-1099"), Arrays.asList("d") });
            testCases.add(new Object[]{ new File(TEST_DATA_DIR, "cram_with_crai_index.cram"), ref, outputExtension, Arrays.asList("1:1-200", "1:1000-2000"), Arrays.asList("a", "d", "e") });
            testCases.add(new Object[]{ new File(TEST_DATA_DIR, "cram_with_crai_index.cram"), ref, outputExtension, Arrays.asList("1:1000-2000", "3:300-400"), Arrays.asList("d", "e", "i", "j") });
            testCases.add(new Object[]{ new File(TEST_DATA_DIR, "cram_with_crai_index.cram"), ref, outputExtension, Arrays.asList("1:1000-2000", "3:300-399"), Arrays.asList("d", "e", "i") });
        }
        return testCases.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ReadCramWithIntervalsIndexTestData")
    public void testReadCramWithIntervalsWithBAIIndex( final File cramFile, final File reference, final String outputExtension,
                                                       final List<String> intervalArgs, final List<String> expectedReadNames ) throws IOException {
        final File outputFile = createTempFile("testReadCramWithIntervalsWithBAIIndex", outputExtension);
        final List<String> args = new ArrayList<>();
        args.addAll(Arrays.asList(
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, cramFile.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, reference.getAbsolutePath()
        ));
        intervalArgs.stream().forEach(intervalArg -> { args.add("-L"); args.add(intervalArg); });

        runCommandLine(args);

        SamAssertionUtils.assertCRAMContentsIfCRAM(outputFile);
        checkReadNames(outputFile, reference, expectedReadNames);
    }

    private void checkReadNames( final File outputFile, final File reference, final List<String> expectedReadNames ) throws IOException {
        List<String> actualReadNames = new ArrayList<>();
        try ( final SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).referenceSequence(reference).open(outputFile) ) {
            for ( SAMRecord read : reader ) {
                actualReadNames.add(read.getReadName());
            }
        }
        Assert.assertEquals(actualReadNames, expectedReadNames, "Read names in output do not match expected read names");
    }

    @Test(dataProvider="testingDataNoRef", expectedExceptions = UserException.MissingReference.class)
    public void testNoRef(String fileIn, String extOut) throws Exception {
        final File outFile = GATKBaseTest.createTempFile(fileIn + ".", extOut);
        File readInput = new File(TEST_DATA_DIR, fileIn);
        final String[] args = new String[]{
                "--input" , readInput.getAbsolutePath(),
                "--output", outFile.getAbsolutePath()
        };
        runCommandLine(args);
    }

    @DataProvider(name="testingDataNoRef")
    public Object[][] testingDataNoRef() {
        return new String[][]{
                {"cramtest.cram", ".sam"}
        };
    }

    @Test(dataProvider="testingDataNoRefMultipleInputs", expectedExceptions = UserException.MissingReference.class)
    public void testNoRefMulti(String fileIn1, String fileIn2, String extOut) throws Exception {
        final File outFile = GATKBaseTest.createTempFile(fileIn1 + ".", extOut);
        File readInput1 = new File(TEST_DATA_DIR, fileIn1);
        File readInput2 = new File(TEST_DATA_DIR, fileIn2);
        final String[] args = new String[]{
                "--input" , readInput1.getAbsolutePath(),
                "--input" , readInput2.getAbsolutePath(),
                "--output", outFile.getAbsolutePath()
        };
        runCommandLine(args);
    }

    @DataProvider(name="testingDataNoRefMultipleInputs")
    public Object[][] testingDataNoRefMulti() {
        return new String[][]{
                {"cramtest.sam", "cramtest.cram", ".bam"}
        };
    }

    // This test case shows that when a CRAM input is provided with a reference that does not have all of the contigs
    // from the CRAM in its sequence dictionary, we throw a UserException.
    @Test(dataProvider="testingDataWrongRef", expectedExceptions = UserException.IncompatibleSequenceDictionaries.class)
    public void testWrongRef(String fileIn, String extOut, String referenceFile) throws Exception {
        final File outFile = GATKBaseTest.createTempFile(fileIn + ".", extOut);
        File readInput = new File(TEST_DATA_DIR, fileIn);
        File reference = new File(TEST_DATA_DIR, referenceFile);
        final String[] args = new String[]{
                "--input" , readInput.getAbsolutePath(),
                "--output", outFile.getAbsolutePath(),
                "-R", reference.getAbsolutePath()
        };
        runCommandLine(args);
    }

    @DataProvider(name="testingDataWrongRef")
    public Object[][] testingDataWrongRef() {
        return new String[][]{
                {"cramtest.cram", ".sam", "cramtestWrongRef.fasta"},
        };
    }

    @DataProvider(name = "serialQueriesOnRemoteFileTest")
    public Object[][] serialQueriesForRemoteFileTest() {
        final File localCRAMWithCrai = new File(publicTestDir +
                "org/broadinstitute/hellbender/engine/CEUTrio.HiSeq.WGS.b37.NA12878.snippet_with_unmapped.cram");

        return new Object[][]{
                { localCRAMWithCrai, new File(b37_reference_20_21), Arrays.asList("20:10000009-10000011", "unmapped"),
                        Arrays.asList("a", "b", "c", "d", "e", "g", "h", "h", "i", "i"),
                        Arrays.asList("g", "h", "h", "i", "i") },
                { localCRAMWithCrai, new File(b37_reference_20_21), Arrays.asList("20:10000009-10000013", "unmapped"),
                        Arrays.asList("a", "b", "c", "d", "e", "f", "f", "g", "h", "h", "i", "i"),
                        Arrays.asList("g", "h", "h", "i", "i") }
        };
    }

    // regression test for https://github.com/broadinstitute/gatk/issues/6475
    @Test(dataProvider = "serialQueriesOnRemoteFileTest")
    public void testSerialQueriesOnRemoteFile(
            final File cramFile,
            final File referenceFile,
            final List<String> intervalStrings,
            final List<String> expectedMappedNames,
            final List<String> expectedUnmappedNames) throws IOException {

        final File cramIndex = SamFiles.findIndex(cramFile);
        try (final FileSystem jimfs = Jimfs.newFileSystem(Configuration.unix())) {
            final Path jimfsCRAM = jimfs.getPath(cramFile.getName());
            final Path jimfsCRAI = jimfs.getPath(cramIndex.getName());
            Files.copy(cramFile.toPath(), jimfsCRAM);
            Files.copy(cramIndex.toPath(), jimfsCRAI);
            final File tempBAMOutFile = createTempFile("testSerialQueriesOnRemoteFile", "bam");

            final ArgumentsBuilder args = new ArgumentsBuilder();
            args.addRaw("-I"); args.addRaw(jimfsCRAM.toUri().toString());
            args.addRaw("-O"); args.addRaw(tempBAMOutFile.getAbsolutePath());
            args.addRaw("-R"); args.addRaw(referenceFile);
            for ( final String intervalString : intervalStrings ) {
                args.addRaw("-L"); args.addRaw(intervalString);
            }
            runCommandLine(args);

            Assert.assertEquals(getReadNames(tempBAMOutFile, s -> s.iterator()), expectedMappedNames);
            Assert.assertEquals(getReadNames(tempBAMOutFile, s -> s.queryUnmapped()), expectedUnmappedNames);
        }
    }

    private final List<String> getReadNames(final File bamFile, final Function<SamReader, SAMRecordIterator> getIt) throws IOException {
        final List<String> allReadNames = new ArrayList<>();
        try (final SamReader samReader = SamReaderFactory.makeDefault().open(bamFile);
             final SAMRecordIterator it = getIt.apply(samReader)) {
            while (it.hasNext()) {
                allReadNames.add(it.next().getReadName());
            }
        }
        return allReadNames;
    }

}