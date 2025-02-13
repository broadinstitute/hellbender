package org.broadinstitute.hellbender.tools.permutect;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class PermutectArgumentConstantsUnitTest extends GATKBaseTest {

    private class DummyPermutectArgCollection {
        @Argument(fullName = PermutectArgumentConstants.NUM_EPOCHS_NAME, doc = "argument from an argument collection", optional = true)
        private String Arg3 = null;
    }

    private class dummyPermutectWrapper extends CommandLineProgram {

        @Argument(fullName = "dummy-argument",doc = "not in python argument list", optional = true)
        private String Arg1 = null;

        // TMP_DIR_NAME = "tmp_dir" // this is a representative inhereited argument that is present in the python argument list

        @Argument(fullName = PermutectArgumentConstants.OUTPUT_NAME, doc = "a standard permutect argument", optional = false)
        private String Arg2 = null;

        @Argument(fullName = PermutectArgumentConstants.INFO_LAYERS_NAME, doc = "in python argument list", optional = true)
        private String Arg3 = null;

        @Argument(fullName = PermutectArgumentConstants.BASE_MODEL_NAME, doc = "in python argument list, has GATK defined default value is overwritten", optional = true)
        private String Arg4 = "THIS_SHOULD_NOT_BE_HERE";

        @Argument(fullName = PermutectArgumentConstants.BATCH_SIZE_NAME, doc = "in python argument list, has GATK defined default value, but is not specified on the cli", optional = true)
        private String Arg4b = "THIS_SHOULD_NOT_BE_HERE";

        @Argument(fullName = PermutectArgumentConstants.DROPOUT_P_NAME, doc = "flag argument", optional = true)
        private boolean Arg5 = false;

        @Argument(fullName = PermutectArgumentConstants.NUM_READ_FEATURES_NAME, doc = "integer arguments", optional = true)
        private int Arg6 = 3;

        @Argument(fullName = PermutectArgumentConstants.AGGREGATION_LAYERS_NAME, doc = "list argument, optional", optional = true)
        private List<String> Arg7 = new ArrayList<>();

        @ArgumentCollection
        DummyPermutectArgCollection args = new DummyPermutectArgCollection();

        @Override
        protected Object doWork() { return null; }
    }

    @Test
    public void testGetPtyhonClassArgumentsFromToolParser() {
        ArgumentsBuilder builder = new ArgumentsBuilder();
        builder.add(PermutectArgumentConstants.OUTPUT_NAME, "output");
        builder.add("dummy-argument", "THIS_SHOULD_NOT_BE_HERE");
        builder.add(PermutectArgumentConstants.INFO_LAYERS_NAME, "info_layers");
        builder.add(PermutectArgumentConstants.BASE_MODEL_NAME, "base_model");
        builder.addFlag(PermutectArgumentConstants.DROPOUT_P_NAME);
        builder.add(PermutectArgumentConstants.AGGREGATION_LAYERS_NAME, "agg1");
        builder.add(PermutectArgumentConstants.AGGREGATION_LAYERS_NAME, "agg2");
        builder.add(PermutectArgumentConstants.AGGREGATION_LAYERS_NAME, "agg3");
        builder.add(PermutectArgumentConstants.AGGREGATION_LAYERS_NAME, "agg4");
        builder.add(PermutectArgumentConstants.NUM_READ_FEATURES_NAME, "2");
        builder.add(PermutectArgumentConstants.NUM_EPOCHS_NAME, "num_epochs");
        CommandLineParser parser = new dummyPermutectWrapper().getCommandLineParser();
        final boolean conversionMap = parser.parseArguments(new PrintStream(System.err), builder.getArgsArray());

        List<String> pyArgs = PermutectArgumentConstants.getPtyhonClassArgumentsFromToolParser(parser);
        Assert.assertTrue(pyArgs.contains("--output"));
        Assert.assertTrue(pyArgs.contains("output"));
        Assert.assertEquals(pyArgs.indexOf("output") - 1, pyArgs.indexOf("--output"));

        Assert.assertTrue(pyArgs.contains("--info_layers"));
        Assert.assertTrue(pyArgs.contains("info_layers"));
        Assert.assertEquals(pyArgs.indexOf("info_layers") - 1, pyArgs.indexOf("--info_layers"));

        Assert.assertTrue(pyArgs.contains("--dropout_p"));

        Assert.assertTrue(pyArgs.contains("--aggregation_layers"));
        Assert.assertTrue(pyArgs.contains("agg1"));
        Assert.assertTrue(pyArgs.contains("agg2"));
        Assert.assertTrue(pyArgs.contains("agg3"));
        Assert.assertTrue(pyArgs.contains("agg4"));
        Assert.assertEquals(pyArgs.indexOf("agg1") - 1, pyArgs.indexOf("--aggregation_layers"));
        Assert.assertEquals(pyArgs.indexOf("agg2") - 1, pyArgs.indexOf("agg1"));
        Assert.assertEquals(pyArgs.indexOf("agg3") - 1, pyArgs.indexOf("agg2"));
        Assert.assertEquals(pyArgs.indexOf("agg4") - 1, pyArgs.indexOf("agg3"));

        Assert.assertTrue(pyArgs.contains("--num_read_features"));
        Assert.assertTrue(pyArgs.contains("2"));

        Assert.assertTrue(pyArgs.contains("--num_epochs"));
        Assert.assertTrue(pyArgs.contains("num_epochs"));
        Assert.assertEquals(pyArgs.indexOf("num_epochs") - 1, pyArgs.indexOf("--num_epochs"));

        Assert.assertFalse(pyArgs.contains("--dummy-argument"));
        Assert.assertFalse(pyArgs.contains("THIS_SHOULD_NOT_BE_HERE"));

        Assert.assertTrue(pyArgs.contains("--base_model"));
        Assert.assertTrue(pyArgs.contains("base_model"));
        Assert.assertEquals(pyArgs.indexOf("base_model") - 1, pyArgs.indexOf("--base_model"));

        Assert.assertFalse(pyArgs.contains("--tmp_dir"));
        Assert.assertFalse(pyArgs.contains("tmp_dir"));

        Assert.assertFalse(pyArgs.contains("--batch_size"));
    }

    @Test
    public void testGenerateArgumentMap() {
        final Map<String, String> conversionMap = PermutectArgumentConstants.PERMUTECT_PYTHON_ARGUMENT_MAP;

        Assert.assertNotNull(conversionMap);
        Assert.assertTrue(conversionMap.entrySet().size() > 30); // assert that the map is not empty and that it reflectively picked up a lot of arguments, the exact number will be subject to change

        for (Map.Entry<String, String> entry : conversionMap.entrySet()) {
            Assert.assertNotNull(entry.getKey());
            Assert.assertFalse(entry.getKey().contains("_"));

            // ptyhon arguments should not contain hyphens
            Assert.assertNotNull(entry.getValue());
            Assert.assertFalse(entry.getValue().contains("-"));
        }

        // various illegal fields that could have snuck into the reflection by acciedent that we want to make sure didn't
        Assert.assertFalse(conversionMap.containsKey("PERMUTECT_PYTHON_ARGUMENT_MAP"));
        Assert.assertFalse(conversionMap.containsKey("dragen-mode"));
        Assert.assertFalse(conversionMap.containsKey("getPythonClassArgumentsFromToolParser"));
        Assert.assertFalse(conversionMap.containsKey("serialVersionUID"));
    }
}