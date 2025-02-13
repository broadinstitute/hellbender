package org.broadinstitute.hellbender.tools.permutect;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.barclay.argparser.NamedArgumentDefinition;

import java.lang.reflect.Field;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static java.util.Map.entry;

public class PermutectArgumentConstants {

    // Java-style (kebab case) without _K suffix
    public static final String STATE_DICT_NAME = "model-state-dict";
    public static final String ARTIFACT_LOG_PRIORS_NAME = "artifact-log-priors";
    public static final String ARTIFACT_SPECTRA_STATE_DICT_NAME = "artifact-spectra-state-dict";
    public static final String HYPERPARAMS_NAME = "hyperparams";
    public static final String NUM_READ_FEATURES_NAME = "num-read-features";
    public static final String NUM_INFO_FEATURES_NAME = "num-info-features";
    public static final String REF_SEQUENCE_LENGTH_NAME = "ref-sequence-length";
    public static final String HIDDEN_LAYERS_NAME = "hidden-layers";
    public static final String NUM_BASE_FEATURES_NAME = "num-base-features";
    public static final String NUM_REF_ALT_FEATURES_NAME = "num-ref-alt-features";

    public static final String SOURCES_NAME = "sources";
    public static final String SOURCE_NAME = "source";

    public static final String INPUT_NAME = "input";
    public static final String OUTPUT_NAME = "output";
    public static final String OUTPUT_DIR_NAME = "output-dir";

    public static final String READ_LAYERS_NAME = "read-layers";
    public static final String SELF_ATTENTION_HIDDEN_DIMENSION_NAME = "self-attention-hidden-dimension";
    public static final String NUM_SELF_ATTENTION_LAYERS_NAME = "num-self-attention-layers";

    public static final String LEARNING_METHOD_NAME = "learning-method";

    public static final String INFO_LAYERS_NAME = "info-layers";
    public static final String AGGREGATION_LAYERS_NAME = "aggregation-layers";
    public static final String CALIBRATION_LAYERS_NAME = "calibration-layers";
    public static final String REF_SEQ_LAYER_STRINGS_NAME = "ref-seq-layer-strings";
    public static final String DROPOUT_P_NAME = "dropout-p";
    public static final String LEARNING_RATE_NAME = "learning-rate";
    public static final String WEIGHT_DECAY_NAME = "weight-decay";
    public static final String BATCH_NORMALIZE_NAME = "batch-normalize";
    public static final String LEARN_ARTIFACT_SPECTRA_NAME = "learn-artifact-spectra";

    public static final String TRAINING_DATASETS_NAME = "training-datasets";
    public static final String TRAIN_TAR_NAME = "train-tar";
    public static final String EVALUATION_TAR_NAME = "evaluation-tar";
    public static final String TEST_DATASET_NAME = "test-dataset";
    public static final String NORMAL_ARTIFACT_DATASETS_NAME = "normal-artifact-datasets";
    public static final String REWEIGHTING_RANGE_NAME = "reweighting-range";
    public static final String BATCH_SIZE_NAME = "batch-size";
    public static final String CHUNK_SIZE_NAME = "chunk-size";
    public static final String NUM_EPOCHS_NAME = "num-epochs";
    public static final String NUM_CALIBRATION_EPOCHS_NAME = "num-calibration-epochs";
    public static final String INFERENCE_BATCH_SIZE_NAME = "inference-batch-size";
    public static final String NUM_WORKERS_NAME = "num-workers";
    public static final String NUM_SPECTRUM_ITERATIONS_NAME = "num-spectrum-iterations";
    public static final String SPECTRUM_LEARNING_RATE_NAME = "spectrum-learning-rate";

    public static final String DATASET_EDIT_TYPE_NAME = "dataset-edit";

    public static final String TENSORBOARD_DIR_NAME = "tensorboard-dir";

    public static final String INITIAL_LOG_VARIANT_PRIOR_NAME = "initial-log-variant-prior";
    public static final String INITIAL_LOG_ARTIFACT_PRIOR_NAME = "initial-log-artifact-prior";
    public static final String CONTIGS_TABLE_NAME = "contigs-table";
    public static final String GENOMIC_SPAN_NAME = "genomic-span";
    public static final String MAF_SEGMENTS_NAME = "maf-segments";
    public static final String NORMAL_MAF_SEGMENTS_NAME = "normal-maf-segments";
    public static final String GERMLINE_MODE_NAME = "germline-mode";
    public static final String NO_GERMLINE_MODE_NAME = "no-germline-mode";
    public static final String HET_BETA_NAME = "het-beta";

    public static final String BASE_MODEL_NAME = "base-model";
    public static final String M3_MODEL_NAME = "permutect-model";
    public static final String PRETRAINED_MODEL_NAME = "pretrained-model";

    @VisibleForTesting
    static final Map<String, String> PERMUTECT_PYTHON_ARGUMENT_MAP = Collections.unmodifiableMap(generateArgumentMap());


    /**
     * Takes in the command line parser for a permutect tool and converts and returns a string list of all of the appropriate arguments
     * for the wrapped python script that are A) actually present for the tool and B) have been set by the user.
     *
     * @param parser the command line parser for the tool in question from which to generate python arguments
     */
    //TODO this might be easier done by directly taking the input arguments directly
    public static List<String> getPtyhonClassArgumentsFromToolParser(CommandLineParser parser) {
        if (parser instanceof CommandLineArgumentParser argParser) {
            List<String> pythonArgs = new ArrayList<>();
            for (Map.Entry<String, String> entry : PERMUTECT_PYTHON_ARGUMENT_MAP.entrySet()) {
                NamedArgumentDefinition arg = argParser.getNamedArgumentDefinitionByAlias(entry.getKey());
                if (arg != null && arg.getHasBeenSet()) { // arg can be null if it is not actually a valid argument for the tool in question
                    pythonArgs.add("--" + entry.getValue());

                    //TODO double check the toString() method for the argument value
                    if (arg.isFlag()) {
                        continue; // flags don't have values
                    } else if (arg.isCollection()) {
                        // The python argument code for permutect expects a sequenctial list of strings following the list argument
                        ((Collection<?>) arg.getArgumentValue()).forEach(value -> pythonArgs.add(value.toString()));
                    } else {
                        pythonArgs.add(arg.getArgumentValue().toString());
                    }
                }
            }
            return pythonArgs;

        } else {
            throw new IllegalArgumentException("command line parser is not CommandLineArgumentParser");
        }
    }

    /**
     * A number of utilities to make converting from the java wrappers to the python methods as easy as possible.
     */
    private static String convertToPythonStyle(String javaStyle) {
        return javaStyle.replace('-', '_');
    }

    /**
     * Generate the static map using reflection.
     */
    public static Map<String, String> generateArgumentMap() {
        return Stream.of(PermutectArgumentConstants.class.getDeclaredFields())
                .filter(field -> java.lang.reflect.Modifier.isStatic(field.getModifiers())
                        && java.lang.reflect.Modifier.isFinal(field.getModifiers())
                        && field.getType().equals(String.class))
                .collect(Collectors.toMap(
                        PermutectArgumentConstants::getFieldValue, // Java-style name
                        field -> convertToPythonStyle(getFieldValue(field)) // Python-style name
                ));
    }

    /**
     * Safely get the value of a static final field.
     */
    private static String getFieldValue(Field field) {
        try {
            return (String) field.get(null);
        } catch (IllegalAccessException e) {
            throw new RuntimeException("Unable to access field: " + field.getName(), e);
        }
    }

}
