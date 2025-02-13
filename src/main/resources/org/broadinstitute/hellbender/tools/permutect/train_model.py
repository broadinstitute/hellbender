import argparse

import psutil
import torch
from torch.utils.tensorboard import SummaryWriter

from permutect import constants, utils
from permutect.architecture.artifact_model import ArtifactModel
from permutect.architecture.artifact_spectra import ArtifactSpectra
from permutect.architecture.posterior_model import plot_artifact_spectra
from permutect.architecture.base_model import load_base_model
from permutect.data.base_dataset import BaseDataset
from permutect.data.artifact_dataset import ArtifactDataset
from permutect.data.base_datum import ArtifactDatum
from permutect.parameters import TrainingParameters, add_training_params_to_parser, parse_training_params, \
    ArtifactModelParameters, parse_artifact_model_params, add_artifact_model_params_to_parser
from permutect.utils import Variation, Label


def train_artifact_model(hyperparams: ArtifactModelParameters, training_params: TrainingParameters, summary_writer: SummaryWriter, dataset: ArtifactDataset):
    model = ArtifactModel(params=hyperparams, num_base_features=dataset.num_base_features, num_ref_alt_features=dataset.num_ref_alt_features, device=utils.gpu_if_available())
    # TODO: magic constant
    model.learn(dataset, training_params, summary_writer=summary_writer, epochs_per_evaluation=10)

    for n, var_type in enumerate(Variation):
        cal_fig, cal_axes = model.calibration[n].plot_calibration()
        summary_writer.add_figure("calibration by count for " + var_type.name, cal_fig)

    return model


def learn_artifact_priors_and_spectra(artifact_dataset: ArtifactDataset, genomic_span_of_data: int):
    artifact_counts = torch.zeros(len(utils.Variation))
    types_list, depths_list, alt_counts_list = [], [], []

    artifact_datum: ArtifactDatum
    for artifact_datum in artifact_dataset:
        if artifact_datum.get_label() != Label.ARTIFACT:
            continue
        variant_type = artifact_datum.get_variant_type()
        artifact_counts[variant_type] += 1
        types_list.append(variant_type)
        counts_and_seq_lks = artifact_datum.one_dimensional_data.get_counts_and_seq_lks()
        depths_list.append(counts_and_seq_lks.depth)
        alt_counts_list.append(counts_and_seq_lks.alt_count)

    # turn the lists into tensors
    types_tensor = torch.LongTensor(types_list)
    depths_tensor = torch.Tensor(depths_list).float()
    alt_counts_tensor = torch.Tensor(alt_counts_list).float()

    log_artifact_priors = torch.log(artifact_counts / genomic_span_of_data)
    artifact_spectra = ArtifactSpectra(num_components=2)

    # TODO: hard-coded num epochs!!!
    artifact_spectra.fit(num_epochs=10, types_b=types_tensor, depths_1d_tensor=depths_tensor,
                         alt_counts_1d_tensor=alt_counts_tensor, batch_size=64)

    return log_artifact_priors, artifact_spectra


def parse_arguments():
    parser = argparse.ArgumentParser(description='train the Permutect artifact model')

    add_artifact_model_params_to_parser(parser)
    add_training_params_to_parser(parser)

    parser.add_argument('--' + constants.LEARN_ARTIFACT_SPECTRA_NAME, action='store_true',
                        help='flag to include artifact priors and allele fraction spectra in saved output.  '
                             'This is worth doing if labeled training data is available but might work poorly '
                             'when Mutect3 generates weak labels based on allele fractions.')
    parser.add_argument('--' + constants.GENOMIC_SPAN_NAME, type=float, required=False,
                        help='Total number of sites considered by Mutect2 in all training data, including those lacking variation or artifacts, hence absent from input datasets.  '
                             'Necessary for learning priors since otherwise rates of artifacts and variants would be overinflated. '
                             'Only required if learning artifact log priors')

    # inputs and outputs
    parser.add_argument('--' + constants.TRAIN_TAR_NAME, type=str, required=True,
                        help='tarfile of training/validation datasets produced by preprocess_dataset.py')
    parser.add_argument('--' + constants.BASE_MODEL_NAME, type=str, help='Base model from train_base_model.py')
    parser.add_argument('--' + constants.OUTPUT_NAME, type=str, required=True, help='path to output saved model file')
    parser.add_argument('--' + constants.TENSORBOARD_DIR_NAME, type=str, default='tensorboard', required=False,
                        help='path to output tensorboard directory')

    return parser.parse_args()


def main_without_parsing(args):
    params = parse_artifact_model_params(args)
    training_params = parse_training_params(args)
    learn_artifact_spectra = getattr(args, constants.LEARN_ARTIFACT_SPECTRA_NAME)
    genomic_span = getattr(args, constants.GENOMIC_SPAN_NAME)

    tensorboard_dir = getattr(args, constants.TENSORBOARD_DIR_NAME)
    summary_writer = SummaryWriter(tensorboard_dir)

    base_model = load_base_model(getattr(args, constants.BASE_MODEL_NAME))
    print(f"Memory usage percent before creating BaseDataset: {psutil.virtual_memory().percent:.1f}")
    base_dataset = BaseDataset(data_tarfile=getattr(args, constants.TRAIN_TAR_NAME), num_folds=10)
    print(f"Memory usage percent before creating ArtifactDataset: {psutil.virtual_memory().percent:.1f}")
    artifact_dataset = ArtifactDataset(base_dataset,
                                       base_model,
                                       base_loader_num_workers=training_params.num_workers,
                                       base_loader_batch_size=training_params.inference_batch_size)
    print(f"Memory usage percent after creating ArtifactDataset: {psutil.virtual_memory().percent:.1f}")

    model = train_artifact_model(hyperparams=params, training_params=training_params, summary_writer=summary_writer, dataset=artifact_dataset)
    print(f"Memory usage percent after training artifact model: {psutil.virtual_memory().percent:.1f}")

    artifact_log_priors, artifact_spectra = learn_artifact_priors_and_spectra(artifact_dataset, genomic_span) if learn_artifact_spectra else (None, None)
    if artifact_spectra is not None:
        art_spectra_fig, art_spectra_axs = plot_artifact_spectra(artifact_spectra, depth=50)
        summary_writer.add_figure("Artifact AF Spectra", art_spectra_fig)

    summary_writer.close()
    model.save_with_base_model(base_model, getattr(args, constants.OUTPUT_NAME), artifact_log_priors, artifact_spectra)


def main():
    args = parse_arguments()
    main_without_parsing(args)


if __name__ == '__main__':
    main()
