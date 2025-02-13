from matplotlib import pyplot as plt
from matplotlib.figure import Figure
import math
from typing import List


# one or more simple plots of y data vs x data on shared axes
def simple_plot(x_y_lab_tuples, x_label, y_label, title):
    fig = plt.figure()
    curve = fig.gca()
    labels_present = False
    for (x, y, lab) in x_y_lab_tuples:
        if lab is not None:
            curve.plot(x, y, label=lab)
            labels_present = True
        else:
            curve.plot(x, y)
    curve.set_title(title)
    curve.set_xlabel(x_label)
    curve.set_ylabel(y_label)
    if labels_present:
        curve.legend()
    return fig, curve


def simple_plot_on_axis(ax, x_y_lab_tuples, x_label, y_label):
    labels_present = False
    for (x, y, lab) in x_y_lab_tuples:
        if lab is not None:
            ax.plot(x, y, label=lab)
            labels_present = True
        else:
            ax.plot(x, y)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    if labels_present:
        ax.legend()


def simple_histograms_on_axis(ax, list_of_histogram_data, list_of_labels, num_bins):
    ax.hist(list_of_histogram_data, bins=num_bins, alpha=0.5, label=list_of_labels)


# apply grouped bar plot to an axis (subplot) object
# heights by category is a dict of category to bar heights, where the nth bar height
# corresponds to the nth x label
def grouped_bar_plot_on_axis(ax, heights_by_category, x_labels, y_label):
    spacing = 5
    bar_width = 0.7 * spacing / len(heights_by_category)

    for n, (category, heights) in enumerate(heights_by_category.items()):
        offset = n * bar_width
        x_positions = [offset + spacing*i for i in range(len(heights))]
        ax.bar(x_positions, heights, width=bar_width, edgecolor='white', label=category)

    # Add xticks on the middle of the group bars
    # plt.xlabel('group', fontweight='bold')
    ticks_offset = bar_width * len(heights_by_category)/2
    ax.set_xticks([ticks_offset + spacing*i for i in range(len(x_labels))], labels=x_labels)
    plt.setp(ax.get_xticklabels(), rotation=90)
    ax.set_ylabel(y_label)
    ax.legend()


# heights by category is a dict of category to bar heights, where the nth bar height
# corresponds to the nth x label
def grouped_bar_plot(heights_by_category, x_labels, y_label):
    fig, ax = plt.subplots()
    grouped_bar_plot_on_axis(ax, heights_by_category, x_labels, y_label)
    return fig, ax


# labels are 0 for non-artifact, 1 for artifact
# predictions_and_labels has form [[(pred, label), (pred, label). . . for roc 1], [likewise for roc 2] etc]
# predictions are logits, not probabilities!!
def plot_accuracy_vs_accuracy_roc_on_axis(lists_of_predictions_and_labels, curve_labels, axis, given_threshold: float = None, sens_prec: bool = False):
    x_y_lab_tuples = []
    small_dots = []
    big_dots = []
    for predictions_and_labels, curve_label in zip(lists_of_predictions_and_labels, curve_labels):
        thresh_and_accs, _ = get_roc_data(predictions_and_labels, given_threshold, sens_prec, use_harmonic_mean=True)
        x_y_lab_tuples.append(([x[1] for x in thresh_and_accs], [x[2] for x in thresh_and_accs], curve_label))

        for threshold, art_acc, non_art_acc in thresh_and_accs:
            if threshold == 0:
                big_dots.append((art_acc, non_art_acc, 'rs'))   # red square
            elif given_threshold is not None and abs(threshold - given_threshold) < 0.001:
                big_dots.append((art_acc, non_art_acc, 'kd'))  # black diamond
            else:
                small_dots.append((art_acc, non_art_acc, 'go'))     # green circle

    simple_plot_on_axis(axis, x_y_lab_tuples, "precision" if sens_prec else "artifact accuracy", "sensitivity" if sens_prec else "non-artifact accuracy")
    for x, y, spec in small_dots:
        axis.plot(x, y, spec, markersize=2,label="")  # point
    for x, y, spec in big_dots:
        axis.plot(x, y, spec, markersize=6,label="")  # point


# similar to the above, but labels are not known and we just have the predicted error probabilities
# we generate a theoretical ROC curve i.e. what the ROC curve would be if these predicted probabilities were
# perfectly calibrated
# labels are 0 for non-artifact, 1 for artifact
# predicted_error_probs is a list of list of floats between 0 and 1
# curve labels is a list of strings
def plot_theoretical_roc_on_axis(predicted_error_probs, curve_labels, axis):
    x_y_lab_tuples = []
    dots = []
    best_thresholds = []
    for error_probs, curve_label in zip(predicted_error_probs, curve_labels):
        thresh_and_accs, best_threshold = get_theoretical_roc_data(error_probs) # best threshold is (threshold, art accuracay, non-art accuracy)
        x_y_lab_tuples.append(([x[1] for x in thresh_and_accs], [x[2] for x in thresh_and_accs], curve_label))

        for threshold, art_acc, non_art_acc in thresh_and_accs:
            dots.append((art_acc, non_art_acc, 'go'))
        dots.append((best_threshold[1], best_threshold[2], 'ro'))
        best_thresholds.append(best_threshold)

    simple_plot_on_axis(axis, x_y_lab_tuples, "precision", "sensitivity")
    for x, y, spec in dots:
        axis.plot(x, y, spec, markersize=2,label="")  # point
    return best_thresholds


# input is list of (predicted artifact logit, binary artifact/non-artifact label) tuples
# 1st output is (threshold, accuracy on artifacts, accuracy on non-artifacts) tuples
# 2nd output is the threshold that maximizes mean (by default harmonic mean, otherwise artithmetic) of accuracy on artifacts and accuracy on non-artifacts
def get_roc_data(predictions_and_labels, given_threshold: float = None, sens_prec: bool = False, use_harmonic_mean: bool = True):
    predictions_and_labels.sort(key=lambda p_and_l: p_and_l[0])  # sort from least to greatest artifact logit
    num_calls = len(predictions_and_labels)
    total_artifact = sum([label for _, label in predictions_and_labels]) + 0.0001
    total_non_artifact = num_calls - total_artifact + 0.0002
    # start at threshold = -infinity; that is, everything is called an artifact, and pick up one variant at a time
    thresh_and_accs = []  # tuples of threshold, accuracy on artifacts, accuracy on non-artifacts
    art_found, non_art_found = total_artifact, 0
    best_art_acc, best_non_art_acc = 1, 0
    next_threshold = -10
    best_threshold, best_mean = -99999, 0
    given_threshold_reached = (given_threshold is None)
    for pred_logit, label in predictions_and_labels:
        art_found -= label  # if labeled as artifact, one artifact has slipped below threshold
        non_art_found += (1 - label)  # if labeled as non-artifact, one non-artifact has been gained
        art_acc, non_art_acc = art_found / total_artifact, non_art_found / total_non_artifact

        # stuff for sensitivity-precision mode
        tp = non_art_found  # non-artifacts that pass threshold are true positives
        fp = total_artifact - art_found  # artifacts that do not fail threshold are false positives

        # in sensitivity-precision mode we care about the precision, not the absolute accuracy of artifact calls
        art_metric = (tp/(tp + fp)) if sens_prec else art_acc

        harmonic_mean = 0 if (art_metric == 0 or non_art_acc == 0) else 1 / ((1 / art_metric) + (1 / non_art_acc))
        arithmetic_mean = (art_acc + non_art_acc) / 2
        mean = harmonic_mean if use_harmonic_mean else arithmetic_mean

        if mean > best_mean:
            best_mean = mean
            best_threshold = pred_logit
            best_art_acc, best_non_art_acc = art_acc, non_art_acc

        if pred_logit > next_threshold:
            thresh_and_accs.append((next_threshold, art_metric, non_art_acc))
            next_threshold = math.ceil(pred_logit)

        # the first time we reach a logit greater than the given threshold, we are basically *at* that threshold
        if not given_threshold_reached and pred_logit > given_threshold:
            thresh_and_accs.append((given_threshold, art_metric, non_art_acc))
            given_threshold_reached = True

    return thresh_and_accs, best_threshold


# input is list of artifact probabilities
# NOTE: this actually includes all errors, such as germline and seq error, but for later possible
# fixing code duplication with the above method we'll call it "artifact"
# 1st output is (threshold, accuracy on non-errors, accuracy on errors) tuples
# 2nd output is the threshold that maximizes harmonic mean of these two accuracies
def get_theoretical_roc_data(artifact_probs):
    artifact_probs.sort(key=lambda p: p)  # sort from least to greatest error probability
    num_calls = len(artifact_probs)
    total_artifact = sum([prob for prob in artifact_probs]) + 0.0001
    total_non_artifact = num_calls - total_artifact + 0.0002
    # start at threshold = 0; that is, everything is called an artifact, and pick up one variant at a time
    # by increasing the probability threshold
    thresh_and_accs = []  # tuples of threshold, accuracy on artifacts, accuracy on non-artifacts
    art_found, non_art_found = total_artifact, 0
    next_threshold = -1
    best_threshold, best_harmonic_mean = (0, 1, 0), 0 # best threshold is threshold, precision, sensitivity
    for prob in artifact_probs:
        art_found -= prob  # lose a fractional artifact
        non_art_found += (1 - prob)  # gain a fractional non-artifact

        tp = non_art_found  # non-artifacts that pass threshold are true positives
        fp = total_artifact - art_found     # artifacts that do not fail threshold are false positives

        # in sensitivity-precision mode we care about the precision, not the absolute accuracy of artifact calls
        sensitivity = tp / total_non_artifact
        precision = tp / (tp + fp)

        harmonic_mean = 0 if (precision == 0 or sensitivity == 0) else 1 / ((1 / sensitivity) + (1 / precision))

        if harmonic_mean > best_harmonic_mean:
            best_harmonic_mean = harmonic_mean
            best_threshold = (prob, precision, sensitivity)

        if prob > next_threshold:
            thresh_and_accs.append((next_threshold, precision, sensitivity))
            next_threshold = math.ceil(prob*20)/20  # we are basically having thresholds of 0.05, 0.1, 0.15. . .
    return thresh_and_accs, best_threshold


def tidy_subplots(figure: Figure, axes, x_label: str = None, y_label: str = None,
                  column_labels: List[str] = None, row_labels: List[str] = None, keep_axes_tick_labels=False):
    """
    Combines various tidying operations on figures with subplots
    1.  Removes the individual axis legends and replaces with a single figure legend.  This assumes
        that all axes have the same lines.
    2.  Show x (y) labels and tick labels only in bottom row (leftmost column)
    3.  Apply column headings and row labels
    4.  Apply overall x and y labels to the figure as a whole

    figure matplotlib.figure.Figure
    axes:   2D array of matplotlib.axes.Axes

    We assume these have been generated together via figure, axes = plt.subplots(. . .)

    """
    handles, labels = figure.get_axes()[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc='upper center')

    for ax in figure.get_axes():
        if not keep_axes_tick_labels:
            ax.label_outer()  # y tick labels only shown in leftmost column, x tick labels only shown on bottom row
        ax.legend().set_visible(False)  # hide the redundant identical subplot legends

        # remove the subplot labels and title -- these will be given manually to the whole figure and to the outer rows
        ax.set_xlabel(None)
        ax.set_ylabel(None)
        ax.set_title(None)

    if x_label is not None:
        figure.supxlabel(x_label)

    if y_label is not None:
        figure.supylabel(y_label)

    if row_labels is not None:
        assert len(row_labels) == len(axes)
        for row_idx, label in enumerate(row_labels):
            axes[row_idx][0].set_ylabel(label)   # note that we use row 0 and set_title to make this a column heading

    if column_labels is not None:
        assert len(column_labels) == len(axes[0])
        for col_idx, label in enumerate(column_labels):
            axes[0][col_idx].set_title(label)

    figure.tight_layout()

