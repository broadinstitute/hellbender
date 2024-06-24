package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.tribble.TribbleException;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;

import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

public class FilterVariantTranchesIntegrationTest  extends CommandLineProgramTest {
    private static final String trancheVCF = largeFileTestDir + "VQSR/g94982_20_1m_10m_python_2dcnn.vcf.gz";
    private static final String snpInput = largeFileTestDir + "VQSR/g94982_20_1m_10m_python_2dcnn.snps.vcf.gz";
    private static final String indelInput = largeFileTestDir + "VQSR/g94982_20_1m_10m_python_2dcnn.indels.vcf.gz";
    private final static String indelTruthVCF = largeFileTestDir + "VQSR/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.sites.20.1M-10M.vcf";
    private final static String snpTruthVCF = largeFileTestDir + "VQSR/Omni25_sites_1525_samples.b37.20.1M-10M.vcf";
    private final static String emptyVCF = largeFileTestDir + "VQSR/emptyB37.vcf";
    private final static String badRef = toolsTestDir + "VQSR/chr20.badRefAlleles.snps.vcf";
    private static final List<String> NO_EXTRA_ARGS = Collections.emptyList();

    @DataProvider(name="getFilteringArgs")
    public Object[][] getFilteringArgs() {
        return new Object[][]{
                //default tranches
                {
                        Arrays.asList(snpTruthVCF, indelTruthVCF), NO_EXTRA_ARGS,
                        largeFileTestDir + "VQSR/expected/g94982_20_1m_10m_tranche_defaults.vcf"
                },
                //with tranches
                {
                        Arrays.asList(snpTruthVCF, indelTruthVCF), Arrays.asList("-indel-tranche", "99.0", "-snp-tranche", "99.5", "--invalidate-previous-filters", "true"),
                        largeFileTestDir + "VQSR/expected/g94982_20_1m_10m_tranched_99_99.5.vcf"
                },
                //duplicate tranches specified
                {
                        Arrays.asList(snpTruthVCF, indelTruthVCF), Arrays.asList("-snp-tranche", "99.5", "--indel-tranche", "99.0", "-snp-tranche", "99.5", "-indel-tranche", "99.0", "--invalidate-previous-filters", "true"),
                        largeFileTestDir + "VQSR/expected/g94982_20_1m_10m_tranched_99_99.5.vcf"
                },
                //multiple tranches
                {
                        Arrays.asList(snpTruthVCF, indelTruthVCF),
                        Arrays.asList("-snp-tranche", "95.0",
                        "-snp-tranche", "99.0",
                        "-snp-tranche", "99.9",
                        "-indel-tranche", "95.0",
                        "-indel-tranche", "99.0",
                        "-indel-tranche", "99.9",
                        "--invalidate-previous-filters", "true"),
                        largeFileTestDir + "VQSR/expected/g94982_20_1m_10m_tranched_95_99_99.9.vcf"
                },
                //different tranche ordering -- expected is same as above
                {
                        Arrays.asList(snpTruthVCF, indelTruthVCF),
                        Arrays.asList("-snp-tranche", "99.9",
                        "-snp-tranche", "95.0",
                        "-snp-tranche", "99.0",
                        "-indel-tranche", "99.9",
                        "-indel-tranche", "95.0",
                        "-indel-tranche", "99.0",
                        "--invalidate-previous-filters", "true"),
                        largeFileTestDir + "VQSR/expected/g94982_20_1m_10m_tranched_95_99_99.9.vcf"
                },
                //the old tranche values are lower than the defaults
                {
                        Arrays.asList(snpTruthVCF, indelTruthVCF),
                        Arrays.asList("-snp-tranche", "99.0","-indel-tranche", "99.0"),
                        largeFileTestDir + "VQSR/expected/g94982_20_1m_10m_tranched_99_with_old_filters.vcf"
                }
        };
    }

    @DataProvider(name="getBadResourceCombos")
    public Object[][] getBadResourceCombos() {
        return new Object[][]{
                //default tranches
                {snpInput, indelTruthVCF}, {indelInput, snpTruthVCF}, {emptyVCF, snpTruthVCF}, {trancheVCF, emptyVCF},
                {badRef, snpTruthVCF}
        };
    }

    @DataProvider(name="getGoodResourceCombos")
    public Object[][] getGoodResourceCombos() {
        return new Object[][]{
                //default tranches
                {snpInput, snpTruthVCF}, {indelInput, indelTruthVCF}
        };
    }

    /**
     * Run the tool on a small test VCF.
     */
    @Test(dataProvider = "getFilteringArgs")
    public void runTrancheFiltering(final List<String> resources, final List<String> extraArgs,
                                    final String expectedOutput) throws IOException {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VARIANT_LONG_NAME, trancheVCF)
                .add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "%s")
                .add("info-key", "MIX_SMALL_2D_W_DROPOUT")
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false")
                .addRaw(StringUtils.join(extraArgs.toArray(), " "));
       resources.stream().forEach(r -> argsBuilder.add(StandardArgumentDefinitions.RESOURCE_LONG_NAME, r));

        final IntegrationTestSpec spec = new IntegrationTestSpec(argsBuilder.toString(),
                Arrays.asList(expectedOutput));
        spec.executeTest("testTrancheFiltering", this);
    }

    @Test(dataProvider = "getBadResourceCombos", expectedExceptions = UserException.BadInput.class)
    public void runTrancheFilteringOnBadInputs(final String inputVcf, final String resourceVcf) {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVcf)
                .add(StandardArgumentDefinitions.RESOURCE_LONG_NAME, resourceVcf)
                .add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "%s")
                .add("info-key", "MIX_SMALL_2D_W_DROPOUT")
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        runCommandLine(argsBuilder);
    }

    @Test(dataProvider = "getGoodResourceCombos")
    public void runTrancheFilteringOnGoodInputs(final String inputVcf, final String resourceVcf) {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVcf)
                .add(StandardArgumentDefinitions.RESOURCE_LONG_NAME, resourceVcf)
                .add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "%s")
                .add("info-key", "MIX_SMALL_2D_W_DROPOUT")
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        runCommandLine(argsBuilder);
    }
}
