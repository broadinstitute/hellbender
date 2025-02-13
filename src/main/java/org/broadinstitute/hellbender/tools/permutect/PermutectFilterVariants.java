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
public class PermutectFilterVariants extends CommandLineProgram {

    public static final String FILTER_VARIANTS = "permutect_filter_variants.py";

    @Argument(
            doc = "Unfiltered input Mutect2 VCF.",
            fullName = PermutectArgumentConstants.INPUT_NAME,
            optional = false
    )
    public String inputName;

    @Argument(
            doc = "Plain text dataset file corresponding to variants in input VCF.",
            fullName = PermutectArgumentConstants.TEST_DATASET_NAME,
            optional = false
    )
    public String testDatasetName;

    @Argument(
            doc = "Trained Permutect model from train_model.py.",
            fullName = PermutectArgumentConstants.M3_MODEL_NAME,
            optional = false
    )
    public String m3ModelName;

    @Argument(
            doc = "Table of contig names vs integer indices.",
            fullName = PermutectArgumentConstants.CONTIGS_TABLE_NAME,
            optional = false
    )
    public String contigsTableName;

    @Argument(
            doc = "Path to output filtered VCF.",
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

    @Argument(
            doc = "Batch size.",
            fullName = PermutectArgumentConstants.BATCH_SIZE_NAME,
            optional = true
    )
    public String batchSize = "64";

    @Argument(
            doc = "Number of subprocesses devoted to data loading, including reading from memory map, collating batches, and transferring to GPU.",
            fullName = PermutectArgumentConstants.NUM_WORKERS_NAME,
            optional = true
    )
    public String numWorkers = "0";

    @Argument(
            doc = "Size in bytes of intermediate binary datasets.",
            fullName = PermutectArgumentConstants.CHUNK_SIZE_NAME,
            optional = true
    )
    public String chunkSize = "100000";

    @Argument(
            doc = "Number of epochs for fitting allele fraction spectra.",
            fullName = PermutectArgumentConstants.NUM_SPECTRUM_ITERATIONS_NAME,
            optional = true
    )
    public String numSpectrumIterations = "10";

    @Argument(
            doc = "Learning rate for fitting allele fraction spectra.",
            fullName = PermutectArgumentConstants.SPECTRUM_LEARNING_RATE_NAME,
            optional = true
    )
    public String spectrumLearningRate = "0.001";

    @Argument(
            doc = "Initial value for natural log prior of somatic variants.",
            fullName = PermutectArgumentConstants.INITIAL_LOG_VARIANT_PRIOR_NAME,
            optional = true
    )
    public String initialLogVariantPrior = "-10.0";

    @Argument(
            doc = "Initial value for natural log prior of artifacts.",
            fullName = PermutectArgumentConstants.INITIAL_LOG_ARTIFACT_PRIOR_NAME,
            optional = true
    )
    public String initialLogArtifactPrior = "-10.0";

    @Argument(
            doc = "Number of sites considered by Mutect2, including those lacking variation or artifacts, hence absent from input dataset. Necessary for learning priors since otherwise rates of artifacts and variants would be overinflated.",
            fullName = PermutectArgumentConstants.GENOMIC_SPAN_NAME,
            optional = false
    )
    public String genomicSpan;

    @Argument(
            doc = "Copy-number segmentation file from GATK containing minor allele fractions. Useful for modeling germline variation as the minor allele fraction determines the distribution of germline allele counts.",
            fullName = PermutectArgumentConstants.MAF_SEGMENTS_NAME,
            optional = true
    )
    public String mafSegmentsName;

    @Argument(
            doc = "Copy-number segmentation file from GATK containing minor allele fractions in the normal/control sample.",
            fullName = PermutectArgumentConstants.NORMAL_MAF_SEGMENTS_NAME,
            optional = true
    )
    public String normalMafSegmentsName;

    @Argument(
            doc = "Flag for genotyping both somatic and somatic variants distinctly but considering both as non-errors (true positives), which affects the posterior threshold set by optimal F1 score.",
            fullName = PermutectArgumentConstants.GERMLINE_MODE_NAME,
            optional = true
    )
    public boolean germlineMode = false;

    @Argument(
            doc = "Beta shape parameter for germline spectrum beta binomial if we want to override binomial.",
            fullName = PermutectArgumentConstants.HET_BETA_NAME,
            optional = true
    )
    public String hetBeta;

    @Argument(
            doc = "Flag for not genotyping germline events so that the only possibilities considered are somatic, artifact, and sequencing error. This is useful for certain validation where pseudo-somatic events are created by mixing germline events at varying fractions.",
            fullName = PermutectArgumentConstants.NO_GERMLINE_MODE_NAME,
            optional = true
    )
    public boolean noGermlineMode = false;

    @Override
    protected Object doWork() {
        PythonScriptExecutor executor = new PythonScriptExecutor(true);
        List<String> pythonifiedArguments = PermutectArgumentConstants.getPtyhonClassArgumentsFromToolParser(getCommandLineParser());

        return executor.executeScript(
                new Resource(FILTER_VARIANTS, PermutectTrainBaseModel.class),
                null,
                pythonifiedArguments);
    }
}
