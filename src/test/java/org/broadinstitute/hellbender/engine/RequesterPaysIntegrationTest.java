package org.broadinstitute.hellbender.engine;

import com.google.cloud.storage.StorageException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.examples.ExampleReadWalkerWithVariants;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class RequesterPaysIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedToolName() {
        return ExampleReadWalkerWithVariants.class.getSimpleName();
    }
    
    @DataProvider
    public Object[][] getRequesterPaysPaths() {

        //Files here are paid for by the bucket owners
        final String OWNER = getGCPTestInputPath() + "org/broadinstitute/hellbender/engine/RequesterPaysIntegrationTest/";

        //Files here are requester pays
        final String REQUESTER = getGCPRequesterPaysBucket() + "test/resources/nio/";

        return new Object[][]{
                {OWNER, OWNER, OWNER, OWNER, false},
                {OWNER, OWNER, OWNER, REQUESTER, true},
                {OWNER, OWNER, REQUESTER, OWNER, true},
                {OWNER, OWNER, REQUESTER, REQUESTER, true},
                {OWNER, REQUESTER, OWNER, OWNER, true},
                {OWNER, REQUESTER, OWNER, REQUESTER, true},
                {OWNER, REQUESTER, REQUESTER, OWNER, true},
                {OWNER, REQUESTER, REQUESTER, REQUESTER, true},
                {REQUESTER, OWNER, OWNER, OWNER, true},
                {REQUESTER, OWNER, OWNER, REQUESTER, true},
                {REQUESTER, OWNER, REQUESTER, OWNER, true},
                {REQUESTER, OWNER, REQUESTER, REQUESTER, true},
                {REQUESTER, REQUESTER, OWNER, OWNER, true},
                {REQUESTER, REQUESTER, OWNER, REQUESTER, true},
                {REQUESTER, REQUESTER, REQUESTER, OWNER, true},
                {REQUESTER, REQUESTER, REQUESTER, REQUESTER, true},
        };
    }

    @Test(dataProvider = "getRequesterPaysPaths", groups="cloud")
    public void testMixedNormalAndRequesterPays(String referenceBase, String bamBase, String vcfBase,
                                                  String intervalBase, boolean requiresRequesterPays) throws IOException {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        final File output = IOUtils.createTempFile("out", ".txt");
        args.addReference(referenceBase + "hg19mini.fasta")
                        .addInput(bamBase + "reads_data_source_test1.bam"  )
                        .addVCF(vcfBase + "example_variants_withSequenceDict.vcf")
                        .addInterval(intervalBase + "hg19mini.all.interval_list")
                        .addOutput(output)
                        .add(StandardArgumentDefinitions.NIO_PROJECT_FOR_REQUESTER_PAYS_LONG_NAME, getGCPTestProject());
                runCommandLine(args);
        IntegrationTestSpec.assertEqualTextFiles(output, new File(packageRootTestDir+"engine/RequesterPaysIntegrationTest/expected_ExampleReadWalkerWithVariantsIntegrationTest_output.txt"));
    }

    @Test(dataProvider = "getRequesterPaysPaths", groups="cloud")
    public void testWithoutRequesterPaysArgument(String referenceBase, String bamBase, String vcfBase,
                                                  String intervalBase, boolean requiresRequesterPays) {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(referenceBase + "hg19mini.fasta")
                .addInput(bamBase + "reads_data_source_test1.bam"  )
                .addVCF(vcfBase + "example_variants_withSequenceDict.vcf")
                .addInterval(intervalBase + "hg19mini.all.interval_list")
                .addOutput(IOUtils.createTempFile("out", ".txt"));
        try{
            runCommandLine(args);
            Assert.assertFalse(requiresRequesterPays, "This shouldn't have reached here because it should if thrown.");
        } catch (final UserException.CouldNotReadInputFile | StorageException ex){
            if( !requiresRequesterPays){
                Assert.fail("This should have thrown an exception.");
            }
        }
    }
}
