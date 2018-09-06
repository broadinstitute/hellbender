package org.broadinstitute.hellbender.tools.walkers.variantutils;

import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.nio.file.*;
import java.util.*;


public class LeftAlignAndTrimVariantsIntegrationTest extends CommandLineProgramTest {
    final Path testDataDir = Paths.get(getToolTestDataDir());

    private Object[] getTestSet(String expectedVcf, String Options) {
        return new Object[]{testDataDir.resolve("test_left_align_hg38.vcf"), Paths.get(b38_reference_20_21), testDataDir.resolve(expectedVcf), Options};
    }

    @DataProvider(name = "LeftAlignDataProvider")
    public Object[][] LeftAlignTestData() {
        return new Object[][]{getTestSet("expected_left_align_hg38.vcf", ""),
                getTestSet("expected_left_align_hg38_split_multiallelics.vcf", " --split-multi-allelics"),
                getTestSet("expected_left_align_hg38_notrim.vcf", " --dont-trim-alleles"),
                getTestSet("expected_left_align_hg38_notrim_split_multiallelics.vcf", " --dont-trim-alleles --split-multi-allelics"),
                getTestSet("expected_left_align_hg38_split_multiallelics_keepOrigAC.vcf", " --split-multi-allelics --keep-original-ac")
        };
    }

    @Test(dataProvider = "LeftAlignDataProvider")
    public void testLeftAlignment(Path inputFile, Path ref, Path expectedOutputFile, String options) throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + ref.toString()
                        + " -V " + inputFile
                        + " -O %s"
                        + " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false"
                        + " --suppress-reference-path "
                        + options,
                Collections.singletonList(expectedOutputFile.toString())
        );
        spec.executeTest("testLeftAlignment--" + expectedOutputFile.toString(), this);
    }

    @DataProvider(name = "LeftAlignRequireReferenceDataProvider")
    public Object[][] LeftAlignRequireReferenceData() {
        return new Object[][]{{testDataDir.resolve("test_left_align_hg38.vcf")}};
    }

    @Test(dataProvider = "LeftAlignRequireReferenceDataProvider")
    public void testLefAlignRequireReference(Path inputFile) throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " -V " + inputFile
                        + " -O %s"
                        + " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false"
                        + " --suppress-reference-path ",
                1, CommandLineException.MissingArgument.class
        );
        spec.executeTest("testLeftAlignment--requireReference", this);
    }
}
