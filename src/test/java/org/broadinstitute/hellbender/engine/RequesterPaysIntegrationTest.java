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

    //Files here are paid for by the bucket owners
    private static final String NOT_REQUESTER = getGCPTestInputPath() + "org/broadinstitute/hellbender/engine/RequesterPaysIntegrationTest/";

    //Files here are requester pays
    private static final String REQUESTER = getGCPRequesterPaysBucket() + "test/resources/nio/";

    @DataProvider
    public Object[][] getRequesterPaysPaths(){
        return new Object[][]{
                {NOT_REQUESTER, NOT_REQUESTER, NOT_REQUESTER, NOT_REQUESTER, false},
                {NOT_REQUESTER, NOT_REQUESTER, NOT_REQUESTER, REQUESTER, true},
                {NOT_REQUESTER, NOT_REQUESTER, REQUESTER, NOT_REQUESTER, true},
                {NOT_REQUESTER, NOT_REQUESTER, REQUESTER, NOT_REQUESTER, true},
                {NOT_REQUESTER, REQUESTER, NOT_REQUESTER, NOT_REQUESTER, true},
                {NOT_REQUESTER, REQUESTER, NOT_REQUESTER, REQUESTER, true},
                {NOT_REQUESTER, REQUESTER, REQUESTER, NOT_REQUESTER, true},
                {NOT_REQUESTER, REQUESTER, REQUESTER, REQUESTER, true},
                {REQUESTER, NOT_REQUESTER, NOT_REQUESTER, NOT_REQUESTER, true},
                {REQUESTER, NOT_REQUESTER, NOT_REQUESTER, REQUESTER, true},
                {REQUESTER, NOT_REQUESTER, REQUESTER, NOT_REQUESTER, true},
                {REQUESTER, NOT_REQUESTER, REQUESTER, NOT_REQUESTER, true},
                {REQUESTER, REQUESTER, NOT_REQUESTER, NOT_REQUESTER, true},
                {REQUESTER, REQUESTER, NOT_REQUESTER, REQUESTER, true},
                {REQUESTER, REQUESTER, REQUESTER, NOT_REQUESTER, true},
                {REQUESTER, REQUESTER, REQUESTER, REQUESTER, true},
        }   ;
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
