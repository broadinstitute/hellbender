package org.broadinstitute.hellbender.tools.permutect;

import org.broadinstitute.barclay.argparser.Argument;

public class PermutectTrainingParamsArgumentCollection {
    @Argument(
            doc = "Learning rate for the model.",
            fullName = PermutectArgumentConstants.LEARNING_RATE_NAME,
            optional = true
    )
    public String learningRate = "0.001";

    @Argument(
            doc = "Weight decay for the optimizer.",
            fullName = PermutectArgumentConstants.WEIGHT_DECAY_NAME,
            optional = true
    )
    public String weightDecay = "0.0";

    @Argument(
            doc = "Batch size for training.",
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
            doc = "Number of epochs for primary training loop.",
            fullName = PermutectArgumentConstants.NUM_EPOCHS_NAME,
            optional = false
    )
    public String numEpochs;

    @Argument(
            doc = "Number of calibration-only epochs.",
            fullName = PermutectArgumentConstants.NUM_CALIBRATION_EPOCHS_NAME,
            optional = true
    )
    public String numCalibrationEpochs = "0";

    @Argument(
            doc = "Batch size when performing model inference (not training).",
            fullName = PermutectArgumentConstants.INFERENCE_BATCH_SIZE_NAME,
            optional = true
    )
    public String inferenceBatchSize = "8192";

}
