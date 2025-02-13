package org.broadinstitute.hellbender.tools.permutect;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.tools.copynumber.GermlineCNVCaller;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.util.ArrayList;
import java.util.List;

@CommandLineProgramProperties(
        summary = "Preprocess plain text training dataset into tarfile of normalized binary data for permutect.",
        oneLineSummary = "Preprocess plain text training dataset into tarfile of normalized binary data for permutect",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public class PermutectPreprocessDataset extends CommandLineProgram {

    public static final String PERMUTECT_PREPREOCESS_DATASET_SCRIPT = "preprocess_dataset.py";

    //TODO handle lists for this? Make it a gatk list?
    @Argument(
            doc = "List of plain text data files.",
            fullName = PermutectArgumentConstants.TRAINING_DATASETS_NAME
    )
    public String trainingDatasetName = null;

    @Argument(
            doc = "Size in bytes of output binary data files. Default is 2e9.",
            fullName = PermutectArgumentConstants.CHUNK_SIZE_NAME,
            optional = true
    )
    public String chunkSizeName = null;

    @Argument(
            doc = "Integer sources corresponding to plain text data files for distinguishing different sequencing conditions.",
            fullName = PermutectArgumentConstants.SOURCES_NAME,
            optional = true
    )
    public String sources = null;

    @Argument(
            doc = "Path to output tarfile of training data.",
            fullName = PermutectArgumentConstants.OUTPUT_NAME
    )
    public String outputTarGz = null;


    @Override
    protected Object doWork() {
        PythonScriptExecutor executor = new PythonScriptExecutor(true);
        List<String> pythonifiedArguments = PermutectArgumentConstants.getPtyhonClassArgumentsFromToolParser(getCommandLineParser());

        return executor.executeScript(
                new Resource(PERMUTECT_PREPREOCESS_DATASET_SCRIPT, PermutectTrainBaseModel.class),
                null,
                pythonifiedArguments);
    }

}