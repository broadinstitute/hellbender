package org.broadinstitute.hellbender.tools.permutect;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.util.List;

@CommandLineProgramProperties(
        summary = "train the Permutect read set representation model.",
        oneLineSummary = "train the Permutect read set representation model",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public class PermutectTrainArtifactModel extends CommandLineProgram {

    public static final String TRAIN_BASE_MODEL_PY = "train_model.py";

    @Argument(
            doc = "Flag to include artifact priors and allele fraction spectra in saved output. This is worth doing if labeled training data is available but might work poorly when Mutect3 generates weak labels based on allele fractions.",
            fullName = PermutectArgumentConstants.LEARN_ARTIFACT_SPECTRA_NAME,
            optional = true
    )
    public boolean learnArtifactSpectra = false;

    @Argument(
            doc = "Total number of sites considered by Mutect2 in all training data, including those lacking variation or artifacts, hence absent from input datasets. Necessary for learning priors since otherwise rates of artifacts and variants would be overinflated. Only required if learning artifact log priors.",
            fullName = PermutectArgumentConstants.GENOMIC_SPAN_NAME,
            optional = true
    )
    public String genomicSpan;

    @Argument(
            doc = "Tarfile of training/validation datasets produced by preprocess_dataset.py.",
            fullName = PermutectArgumentConstants.TRAIN_TAR_NAME,
            optional = false
    )
    public String trainTarName;

    @Argument(
            doc = "Base model from train_base_model.py.",
            fullName = PermutectArgumentConstants.BASE_MODEL_NAME,
            optional = true
    )
    public String baseModelName;

    @Argument(
            doc = "Path to output saved model file.",
            fullName = PermutectArgumentConstants.OUTPUT_NAME,
            optional = false
    )
    public String outputName;

    @Argument(
            doc = "Path to output tensorboard directory.",
            fullName = PermutectArgumentConstants.TENSORBOARD_DIR_NAME,
            optional = true
    )
    public String tensorboardDirName = "tensorboard";

    // Shared argument collections to include in arguments
    @ArgumentCollection
    PermutectArtifactModelArgumentCollection artifactModelArgs = new PermutectArtifactModelArgumentCollection();
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
