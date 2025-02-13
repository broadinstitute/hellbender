package org.broadinstitute.hellbender.tools.permutect;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.util.ArrayList;
import java.util.List;

@CommandLineProgramProperties(
        summary = "train the Permutect read set representation model.",
        oneLineSummary = "train the Permutect read set representation model",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public class PermutectTrainBaseModel extends CommandLineProgram {

    public static final String TRAIN_BASE_MODEL_PY = "train_base_model.py";

    @Argument(
            doc = "Options [SUPERVISED, SEMISUPERVISED, SUPERVISED_CLUSTERING, AFFINE, MASK_PREDICTION, AUTOENCODER, DEEPSAD, MARS].",
            fullName = PermutectArgumentConstants.LEARNING_METHOD_NAME,
            optional = true
    )
    public String trainingDatasetName = null;

    @Argument(
            doc = "Tarfile of training/validation datasets produced by preprocess_dataset.",
            fullName = PermutectArgumentConstants.TRAIN_TAR_NAME,
            optional = false
    )
    public String chunkSizeName = null;

    @Argument(
            doc = "Output location for the saved model file.",
            fullName = PermutectArgumentConstants.OUTPUT_NAME,
            optional = false
    )
    public String sources = null;

    @Argument(
            doc = "output tensorboard directory.",
            fullName = PermutectArgumentConstants.TENSORBOARD_DIR_NAME,
            optional = true
    )
    public String outputTarGz = null;

    // Shared argument collections to include in arguments
    @ArgumentCollection
    PermutectBaseModelArgumentCollection baseArgumentCollection = new PermutectBaseModelArgumentCollection();
    @ArgumentCollection
    PermutectTrainingParamsArgumentCollection trainingParamsArgumentCollection = new PermutectTrainingParamsArgumentCollection();

    @Override
    protected Object doWork() {
        PythonScriptExecutor executor = new PythonScriptExecutor(true);
        List<String> pythonifiedArguments = PermutectArgumentConstants.getPtyhonClassArgumentsFromToolParser(getCommandLineParser());

        return executor.executeScript(
                new Resource(TRAIN_BASE_MODEL_PY, PermutectTrainBaseModel.class),
                null,
                pythonifiedArguments);
    }


}