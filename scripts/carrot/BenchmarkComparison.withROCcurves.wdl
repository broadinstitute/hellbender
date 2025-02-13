version 1.0

import "https://raw.githubusercontent.com/broadinstitute/palantir-workflows/2e039c7795e0b0dfdad23ab56bf9be18247632a5/BenchmarkVCFs/BenchmarkVCFs.wdl" as Benchmark

workflow BenchmarkComparison {
    input{
        String? analysisRegion
        File evalVcf
        String evalLabel
        File evalVcfIndex
        Array[File] evalRuntimeSummaries
        File evalMonitoringExample

        File controlVcf
        String controlLabel
        File controlVcfIndex
        Array[File] controlRuntimeSummaries
        File controlMonitoringExample


        File truthVcf
        File confidenceInterval
        String truthLabel
        File truthVcfIndex
        File reference
        File refIndex
        File refDict
        File hapMap
        Array[File] stratIntervals = []
        Array[String] stratLabels = []
        Array[String]? jexlVariantSelectors
        Array[String]? variantSelectorLabels
        String referenceVersion
        Int? threadsVcfEval=2
        Boolean doIndelLengthStratification=true
        Int? preemptible
        String gatkTag="4.0.11.0"
        Boolean requireMatchingGenotypes=true
        File? gatkJarForAnnotation
        Array[String]? annotationNames
        Boolean enableRefOverlap = false
        Boolean passingOnly=true
        String? vcfScoreField
        String? dummyInputForTerraCallCaching

        File monitoring_log_process_script = "gs://emeryj-testing/plot_monitoring_data.R"
    }

    meta {
        description: "A workflow to calculate sensitivity and precision of a germline variant calling pipeline by comparing a 'call' vcf produced by the pipeline to a gold standard 'truth' vcf.  Allows for stratification based on interval lists, bed files, or variant types defined according to GATK SelectVariants."
    }

    parameter_meta {
        evalVcf: {description: "vcfs to be evaluated"}
        evalLabel: {description: "label to identify vcf to be evaluated"}
        evalVcfIndex: {description: "vcf index for evalVcf"}
        truthVcf: {description: "truth vcf against which to evaluate"}
        truthLabel: {description: "label by which to indentify truth set"}
        confidenceInterval: {description: "confidence interval for truth set (can be bed or picard interval_list)"}
        reference: {description: "reference fasta"}
        refIndex: {description: "reference index"}
        refDict: {description: "reference dict"}
        hapMap: {description: "reference haplotype map for CrosscheckFingerprints"}
        referenceVersion: {description: "reference version used, for igv xml session (must be either 'hg19' or 'hg38')"}
        stratIntervals: {description: "intervals for stratifiction (can be picard interval_list or bed format)"}
        stratLabels: {description: "labels by which to identify stratification intervals (must be same length as stratIntervals)"}
        jexlVariantSelectors: {description: "variant types to select over (defined by jexl fed to GATK SelectVariants)"}
        variantSelectorLabels: {description: "labels by which to identify variant selectors (must be same length as jexlVariantSelectors)"}
        doIndelLengthStratification: {description: "whether or not to perform stratification by indel length"}
        requireMatchingGenotypes: {description: "whether to require genotypes to match in order to be a true positive"}
        gatkTag: {description: "version of gatk docker to use.  Defaults to 4.0.11.0"}
        analysisRegion: {description: "if provided (gatk format, single interval e.g., 'chr20', or 'chr20:1-10') all the analysis will be performed only within the region."}
        passingOnly: {description:"Have vcfEval only consider the passing variants"}
        vcfScoreField: {description:"Have vcfEval use this field for making the roc-plot. If this is an info field (like VSQLOD) it should be provided as INFO.VQSLOD, otherewise it is assumed to be a format field."}
        gatkJarForAnnotation: {description:"GATK jar that can calculate necessary annotations for jexl Selections when using VCFEval."}
        annotationNames: {description:"Annotation arguments to GATK (-A argument, multiple OK)"}
        dummyInputForTerraCallCaching: {description:"When running on Terra, use workspace.name as this input to ensure that all tasks will only cache hit to runs in your own workspace. This will prevent call caching from failing with 'Cache Miss (10 failed copy attempts)'. Outside of Terra this can be left empty. This dummy input is only needed for tasks that have no inputs specific to the sample being run (such as CreateIntervalList which does not take in any sample data)."}
    }

    call HandleRuntimesTask as CONTROLRuntimeTask {
        input:
            time_summaries = controlRuntimeSummaries,
            representative_monitoring_log = controlMonitoringExample,
            monitoring_log_process_script = monitoring_log_process_script
    }

    call HandleRuntimesTask as EVALRuntimeTask {
        input:
            time_summaries = evalRuntimeSummaries,
            representative_monitoring_log = evalMonitoringExample,
            monitoring_log_process_script = monitoring_log_process_script
    }

    call Benchmark.Benchmark as BenchmarkVCFTestSample{
        input:
            analysisRegion = analysisRegion,
            evalVcf = evalVcf,
            evalLabel = evalLabel,
            evalVcfIndex = evalVcfIndex,
            truthVcf = truthVcf,
            confidenceInterval = confidenceInterval,
            truthLabel = truthLabel,
            truthVcfIndex = truthVcfIndex,
            reference = reference,
            refIndex = refIndex,
            refDict = refDict,
            hapMap = hapMap,
            stratIntervals = stratIntervals,
            stratLabels = stratLabels,
            jexlVariantSelectors = jexlVariantSelectors,
            variantSelectorLabels = variantSelectorLabels,
            referenceVersion = referenceVersion,
            doIndelLengthStratification = doIndelLengthStratification,
            gatkTag = gatkTag,
            requireMatchingGenotypes = requireMatchingGenotypes,
            passingOnly = passingOnly,
            vcfScoreField = vcfScoreField,
            gatkJarForAnnotation = gatkJarForAnnotation,
            annotationNames = annotationNames
    }

    call Benchmark.Benchmark as BenchmarkVCFControlSample{
        input:
            analysisRegion = analysisRegion,
            evalVcf = controlVcf,
            evalLabel = controlLabel,
            evalVcfIndex = controlVcfIndex,
            truthVcf = truthVcf,
            confidenceInterval = confidenceInterval,
            truthLabel = truthLabel,
            truthVcfIndex = truthVcfIndex,
            reference = reference,
            refIndex = refIndex,
            refDict = refDict,
            hapMap = hapMap,
            stratIntervals = stratIntervals,
            stratLabels = stratLabels,
            jexlVariantSelectors = jexlVariantSelectors,
            variantSelectorLabels = variantSelectorLabels,
            referenceVersion = referenceVersion,
            doIndelLengthStratification = doIndelLengthStratification,
            gatkTag = gatkTag,
            requireMatchingGenotypes = requireMatchingGenotypes,
            passingOnly = passingOnly,
            vcfScoreField = vcfScoreField,
            gatkJarForAnnotation = gatkJarForAnnotation,
            annotationNames = annotationNames
    }

    Array[File] roc_tables_1 = flatten([
        select_all(BenchmarkVCFTestSample.snpRocs),
        select_all(BenchmarkVCFTestSample.nonSnpRocs),
        select_all(BenchmarkVCFControlSample.snpRocs),
        select_all(BenchmarkVCFControlSample.nonSnpRocs)
    ])

    call PlotROCTask as ROCPlot {
        input:
            sample_id = truthLabel,
            roc_tables = roc_tables_1,
            tool1_label = evalLabel,
            tool2_label = controlLabel,
            stratifiers = select_first([stratLabels, []]),
            preemptible = preemptible
    }
    Array[File] roc_plots = ROCPlot.plots
    call MergePNGs as MergeROC {
        input:
            pngs = roc_plots
    }

    output {
        File? evalsummary=BenchmarkVCFTestSample.summary
        Float evalsnpPrecision=BenchmarkVCFTestSample.snpPrecision
        Float evalindelPrecision=BenchmarkVCFTestSample.indelPrecision
        Float evalsnpRecall=BenchmarkVCFTestSample.snpRecall
        Float evalindelRecall=BenchmarkVCFTestSample.indelRecall
        Float evalsnpF1Score=BenchmarkVCFTestSample.snpF1Score
        Float evalindelF1Score=BenchmarkVCFTestSample.indelF1Score

        File? controlsummary=BenchmarkVCFControlSample.summary
        Float controlsnpPrecision=BenchmarkVCFControlSample.snpPrecision
        Float controlindelPrecision=BenchmarkVCFControlSample.indelPrecision
        Float controlsnpRecall=BenchmarkVCFControlSample.snpRecall
        Float controlindelRecall=BenchmarkVCFControlSample.indelRecall
        Float controlsnpF1Score=BenchmarkVCFControlSample.snpF1Score
        Float controlindelF1Score=BenchmarkVCFControlSample.indelF1Score

        Array[File] rocPlots = ROCPlot.plots
        File rocPlotsmerged = MergeROC.plots

        Float controlHCwallclockhours=CONTROLRuntimeTask.HC_wallclock_hours
        Float controlHCwallclockmax=CONTROLRuntimeTask.HC_wallclock_max
        Float controlHCprocesshours=CONTROLRuntimeTask.HC_process_hours
        Float controlHCsystemhours=CONTROLRuntimeTask.HC_system_hours
        File controlMonitoringLogs=CONTROLRuntimeTask.monitoring_logs

        Float evalHCwallclockhours=EVALRuntimeTask.HC_wallclock_hours
        Float evalHCwallclockmax=EVALRuntimeTask.HC_wallclock_max
        Float evalHCprocesshours=EVALRuntimeTask.HC_process_hours
        Float evalHCsystemhours=EVALRuntimeTask.HC_system_hours
        File evalMonitoringLogs=EVALRuntimeTask.monitoring_logs
    }
}

task HandleRuntimesTask {
    input {
        Array[File] time_summaries
        File representative_monitoring_log
        File monitoring_log_process_script
    }

    command <<<
        set -xeuo pipefail

        cat <<'EOF' > sumRuntime.py
import os
import sys

files = ["~{sep='",\n"' time_summaries}"]

systemRuntimes = []
usrRuntimes = []
osRuntimes = []
print(files)

for file in files:
    lines = open(file, 'r').readlines()
    systemRuntimes.append(float(lines[0]))
    usrRuntimes.append(float(lines[1]))
    osRuntimes.append(float(lines[2]))

print(systemRuntimes)
print("wallclock hours for gatk: " +str(sum(systemRuntimes)/3600))
print("max shard wallclock for gatk: " +str(max(systemRuntimes)/3600))
open("wallclocksum.txt", 'w').write(str((sum(systemRuntimes)/3600)))
open("wallclockmax.txt", 'w').write(str((max(systemRuntimes)/3600)))
print(usrRuntimes)
print("process cpu hours for gatk: " +str(sum(usrRuntimes)/3600))
open("processcpusum.txt", 'w').write(str((sum(usrRuntimes)/3600)))
print(osRuntimes)
print("system cpu hours for gatk: " +str(sum(osRuntimes)/3600))
open("systemcpusum.txt", 'w').write(str((sum(osRuntimes)/3600)))
EOF
        python sumRuntime.py

        DEBIAN_FRONTEND=noninteractive
        TZ=America/Boston
        wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
        echo "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" >> /etc/apt/sources.list
        apt update && apt install r-base -y
        chmod a+x ~{monitoring_log_process_script}
        ~{monitoring_log_process_script} ~{representative_monitoring_log} monitoring.pdf

    >>>

    output {
        Float HC_wallclock_hours = read_float("wallclocksum.txt")
        Float HC_wallclock_max = read_float("wallclockmax.txt")
        Float HC_process_hours = read_float("processcpusum.txt")
        Float HC_system_hours = read_float("systemcpusum.txt")
        File monitoring_logs = "monitoring.pdf"
    }

    runtime {
        docker: "docker.io/broadinstitute/gatk:gatkbase-3.2.0"
        disks: "local-disk 200 HDD"
        maxRetries: 3
    }
 }

task PlotROCTask {
    input {
        String sample_id
        Array[File] roc_tables
        String tool1_label
        String tool2_label
        Array[String] stratifiers
        String? additional_label
        Int? mem_gb
        Int? preemptible
    }

    Int machine_mem_gb = select_first([mem_gb, 8])

    String additional_label_arg = if defined(additional_label) then "--additional-label \"" + additional_label + "\"" else ""

    command <<<
        set -xeuo pipefail

        source activate fe_evaluation

        cat <<'EOF' > script.py
import matplotlib
import matplotlib.pyplot as plt
import csv
import gzip
import argparse
import os

matplotlib.rcParams['text.usetex'] = False
matplotlib.rcParams['mathtext.default'] = 'regular'
matplotlib.rcParams['font.family'] = 'serif'


def f1(tp, fp, fn):
    return tp / (tp + 0.5 * (fp + fn))

def read_roc_data(filename):
    data = dict()
    max_f1 = 0
    best_qual = None
    with gzip.open(filename, 'rt') as roc_file:
        for line in roc_file:
            if line.startswith('#score field'):
                break
        for line in csv.DictReader(roc_file, delimiter='\t'):
            qual = float(line['#score'])
            f1_score = f1(float(line['true_positives_baseline']), float(line['false_positives']), float(line['false_negatives']))
            if f1_score > max_f1:
                best_qual = (qual, f1_score)
            false_positives = float(line['false_positives'])
            sensitivity = float(line['sensitivity'])
            data[qual] = (false_positives, sensitivity)
    return data, best_qual


def read_file(data, best_qual, filename, stratifiers):
    basename = os.path.basename(filename)
    info = basename.split('_')
    print(info)
    tool = "tool1" if "TEST" in info[0] else "tool2"
    print(tool)

    if info[4] == 'vcfeval':
        region = 'all'
    else:
        region = info[4]
    if region not in stratifiers:
        raise RuntimeError('Invalid stratifier {} in file {}'.format(region, filename))

    if info[5 if region == 'all' else 6] == 'non':
        var_type = 'indel'
    else:
        var_type = 'snp'

    print((region, var_type, tool))
    data[(region, var_type, tool)], best_qual[(region, var_type, tool)] = read_roc_data(filename)

def plot_roc(ax, data, best_qual, region, var_type):
    print(data)
    X1 = [x[0] for x in data[(region, var_type, 'tool1')].values()]
    Y1 = [x[1] for x in data[(region, var_type, 'tool1')].values()]
    X2 = [x[0] for x in data[(region, var_type, 'tool2')].values()]
    Y2 = [x[1] for x in data[(region, var_type, 'tool2')].values()]

    ax.plot(X1, Y1, c='C0', marker='.', markersize=2, linestyle=' ')
    ax.plot(X2, Y2, c='C1', marker='.', markersize=2, linestyle=' ')

    best_qual_tool1 = best_qual[(region, var_type, 'tool1')]
    best_qual_tool2 = best_qual[(region, var_type, 'tool2')]

    best_coordinates_tool1 = data[(region, var_type, 'tool1')][best_qual_tool1[0]]
    best_coordinates_tool2 = data[(region, var_type, 'tool2')][best_qual_tool2[0]]

    ax.plot([best_coordinates_tool1[0]], [best_coordinates_tool1[1]], c='C0', marker='x', markersize=15, linestyle=' ')
    ax.plot([best_coordinates_tool2[0]], [best_coordinates_tool2[1]], c='C1', marker='x', markersize=15, linestyle=' ')

    ax.annotate(r'$F_{1, max}$ = '+'{:.3f} @ Q{:.0f}'.format(best_qual_tool1[1], best_qual_tool1[0]), (0.95, 0.05), xycoords='axes fraction', xytext=(0, 12), textcoords='offset points', c='C0', ha='right', va='bottom')
    ax.annotate(r'$F_{1, max}$ = '+'{:.3f} @ Q{:.0f}'.format(best_qual_tool2[1], best_qual_tool2[0]), (0.95, 0.05), xycoords='axes fraction', xytext=(0, 0), textcoords='offset points', c='C1', ha='right', va='bottom')

    ax.set_xlabel(r'FP')
    ax.set_ylabel(r'Sensitivity')

    ax.set_title(r'{} {}'.format(var_type, region), zorder=0)

def plot_data(data, best_qual, stratifiers, sample_id, tool1_label, tool2_label, additional_label):
    num_columns = max(len(stratifiers), 3)
    fig, axes = plt.subplots(3, num_columns, figsize=(6*num_columns,16))
    for row, var_type in enumerate(['snp', 'indel']):
        for col, region in enumerate(stratifiers):
            column_to_plot = col if len(stratifiers) > 1 else 1
            ax = axes[row, column_to_plot]
            plot_roc(ax, data, best_qual, region, var_type)

    # Clear axes for legend
    for i in range(num_columns):
        axes[2, i].axis('off')

    # Clear non-used axes if plotting less than 3 columns
    if len(stratifiers) == 1:
        axes[0, 0].axis('off')
        axes[0, 2].axis('off')
        axes[1, 0].axis('off')
        axes[1, 2].axis('off')
    if len(stratifiers) == 2:
        axes[0, 2].axis('off')
        axes[1, 2].axis('off')


    fig.suptitle('Sample: {}'.format(sample_id) + ('' if not additional_label else ', {}'.format(additional_label)) + '\n')

    legend_tool1_line = matplotlib.lines.Line2D([], [], color='C0', label=tool1_label)
    legend_tool2_line = matplotlib.lines.Line2D([], [], color='C1', label=tool2_label)
    fig.legend(bbox_to_anchor=(0.5, 0.2), loc='center', handles=[legend_tool1_line, legend_tool2_line])
    plt.tight_layout()
    fig.savefig('roc_plot_{}.png'.format(sample_id), dpi=100)

def main(roc_tables, sample_id, tool1_label, tool2_label, stratifiers, additional_label):
    if stratifiers is None:
        stratifiers = ['all']
    else:
        stratifiers.insert(0, 'all')

    data = dict()
    best_qual = dict()
    for filename in roc_tables:
        read_file(data, best_qual, filename, stratifiers)

    plot_data(data, best_qual, stratifiers, sample_id, tool1_label, tool2_label, additional_label)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create F1 functional equivalence plots.')
    parser.add_argument('--additional-label', type=str)
    parser.add_argument('--stratifiers', type=str, nargs='*')
    required_named = parser.add_argument_group('Required named arguments')
    required_named.add_argument('--sample-id', type=str)
    required_named.add_argument('--tool1', required=True, type=str)
    required_named.add_argument('--tool2', required=True, type=str)
    required_named.add_argument('--roc-tables', type=str, nargs='+')
    args = parser.parse_args()
    main(args.roc_tables, args.sample_id, args.tool1, args.tool2, args.stratifiers, args.additional_label)
EOF

        python script.py --sample-id "~{sample_id}" --tool1 "~{tool1_label}" --tool2 "~{tool2_label}" ~{additional_label_arg} --roc-tables ~{sep=' ' roc_tables} --stratifiers ~{sep=' ' stratifiers}
    >>>

    output {
        Array[File] plots = glob("*.png")
    }

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/functionalequivalence/fe_evaluation:1.0.0"
        preemptible: select_first([preemptible, 0])
        disks: "local-disk 200 HDD"
        maxRetries: 3
    }
}

task MergePNGs {
    input {
        Array[File] pngs
    }

    command {
        convert -append ~{sep=" " pngs} plots.png
    }

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/functionalequivalence/merge_pngs:1.0.0"
        memory: "8 GB"
        disks: "local-disk 20 HDD"
        maxRetries: 3
    }

    output{
        File plots = "plots.png"
    }
}
