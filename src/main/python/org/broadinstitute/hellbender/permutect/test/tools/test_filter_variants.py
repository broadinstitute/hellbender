import tempfile
from argparse import Namespace

from permutect import constants
from permutect.tools import filter_variants


def test_filtering_on_dream1_chr20():
    # Inputs
    artifact_model = '/Users/davidben/mutect3/permutect/integration-tests/singular-10-Mb/artifact-model.pt'

    mutect2_vcf = '/Users/davidben/mutect3/permutect/integration-tests/dream1-chr20/mutect2_chr20.vcf'
    maf_segments = '/Users/davidben/mutect3/permutect/integration-tests/dream1-chr20/segments.table'
    contigs_table = '/Users/davidben/mutect3/permutect/integration-tests/dream1-chr20/contigs.table'
    filtering_dataset = '/Users/davidben/mutect3/permutect/integration-tests/dream1-chr20/test_chr20.dataset'

    # Outputs
    permutect_vcf = tempfile.NamedTemporaryFile()
    tensorboard_dir = tempfile.TemporaryDirectory()

    filtering_args = Namespace()
    setattr(filtering_args, constants.INPUT_NAME, mutect2_vcf)
    setattr(filtering_args, constants.TEST_DATASET_NAME, filtering_dataset)
    setattr(filtering_args, constants.M3_MODEL_NAME, artifact_model)
    setattr(filtering_args, constants.OUTPUT_NAME, permutect_vcf.name)
    setattr(filtering_args, constants.TENSORBOARD_DIR_NAME, tensorboard_dir.name)
    setattr(filtering_args, constants.BATCH_SIZE_NAME, 64)
    setattr(filtering_args, constants.NUM_WORKERS_NAME, 0)
    setattr(filtering_args, constants.CHUNK_SIZE_NAME, 100000)
    setattr(filtering_args, constants.NUM_SPECTRUM_ITERATIONS_NAME, 2)
    setattr(filtering_args, constants.HET_BETA_NAME, 10)
    setattr(filtering_args, constants.SPECTRUM_LEARNING_RATE_NAME, 0.001)
    setattr(filtering_args, constants.INITIAL_LOG_VARIANT_PRIOR_NAME, -10.0)
    setattr(filtering_args, constants.INITIAL_LOG_ARTIFACT_PRIOR_NAME, -10.0)
    setattr(filtering_args, constants.GENOMIC_SPAN_NAME, 60000000)
    setattr(filtering_args, constants.MAF_SEGMENTS_NAME, None)
    setattr(filtering_args, constants.CONTIGS_TABLE_NAME, contigs_table)
    setattr(filtering_args, constants.NORMAL_MAF_SEGMENTS_NAME, None)
    setattr(filtering_args, constants.GERMLINE_MODE_NAME, False)
    setattr(filtering_args, constants.NO_GERMLINE_MODE_NAME, False)

    filter_variants.main_without_parsing(filtering_args)
    h = 9