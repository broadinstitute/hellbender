package org.broadinstitute.hellbender.tools.permutect;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.util.List;

@CommandLineProgramProperties(
        summary = "train the Mutect3 artifact model.", //TODO this needs to be properly labeled
        oneLineSummary = "train the Mutect3 artifact modela",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public class PermutectEditDataset extends CommandLineProgram {
    public static final String PERMUTECT_EDIT_DATASET = "edit_dataset.py";

    @Argument(
            doc = "Size in bytes of output binary data files.",
            fullName = PermutectArgumentConstants.CHUNK_SIZE_NAME,
            optional = true
    )
    public String chunkSize = String.valueOf((int) 2e9);

    @Argument(
            doc = "How to modify the dataset.",
            fullName = PermutectArgumentConstants.DATASET_EDIT_TYPE_NAME,
            optional = false
    )
    public String datasetEditType;

    @Argument(
            doc = "New source integer to apply.",
            fullName = PermutectArgumentConstants.SOURCE_NAME,
            optional = true
    )
    public String source;

    @Argument(
            doc = "Tarfile(s) of training/validation datasets produced by preprocess_dataset.py.",
            fullName = PermutectArgumentConstants.TRAIN_TAR_NAME,
            optional = false
    )
    public List<String> trainTar;

    @Argument(
            doc = "Path to pruned dataset file.",
            fullName = PermutectArgumentConstants.OUTPUT_NAME,
            optional = false
    )
    public String output;

    @Override
    protected Object doWork() {
        PythonScriptExecutor executor = new PythonScriptExecutor(true);
        List<String> pythonifiedArguments = PermutectArgumentConstants.getPtyhonClassArgumentsFromToolParser(getCommandLineParser());

        return executor.executeScript(
                new Resource(PERMUTECT_EDIT_DATASET, PermutectTrainBaseModel.class),
                null,
                pythonifiedArguments);
    }
}
