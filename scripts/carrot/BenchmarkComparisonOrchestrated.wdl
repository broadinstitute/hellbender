version 1.0

import "https://raw.githubusercontent.com/broadinstitute/gatk/je_carrotHCTestsScripts/scripts/carrot/BenchmarkComparison.withROCcurves.wdl" as BenchmarkAndRocCurves

workflow BenchmarkVCFsHeadToHeadOrchestrated {
       input{
            File CHM_evalVcf
            String CHM_evalLabel
            File CHM_evalVcfIndex
            Array[File] CHM_evalRuntimeSummaries
            File CHM_evalMonitoringExample
            File CHM_controlVcf
            String CHM_controlLabel
            File CHM_controlVcfIndex
            Array[File] CHM_controlRuntimeSummaries
            File CHM_controlMonitoringExample

            File CHM_truthVcf
            File CHM_confidenceInterval
            String CHM_truthLabel
            File CHM_truthVcfIndex


            File NIST_evalVcf
            String NIST_evalLabel
            File NIST_evalVcfIndex
            Array[File] NIST_evalRuntimeSummaries
            File NIST_evalMonitoringExample
            File NIST_controlVcf
            String NIST_controlLabel
            File NIST_controlVcfIndex
            Array[File] NIST_controlRuntimeSummaries
            File NIST_controlMonitoringExample

            File NIST_truthVcf
            File NIST_confidenceInterval
            String NIST_truthLabel
            File NIST_truthVcfIndex


            File EXOME1_evalVcf
            String EXOME1_evalLabel
            File EXOME1_evalVcfIndex
            Array[File] EXOME1_evalRuntimeSummaries
            File EXOME1_evalMonitoringExample
            File EXOME1_controlVcf
            String EXOME1_controlLabel
            File EXOME1_controlVcfIndex
            Array[File] EXOME1_controlRuntimeSummaries
            File EXOME1_controlMonitoringExample

            File EXOME1_truthVcf
            File EXOME1_confidenceInterval
            String EXOME1_truthLabel
            File EXOME1_truthVcfIndex

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

           File monitoring_log_process_script = "gs://dsp-methods-carrot-data/test_data/haplotypecaller_tests/plot_monitoring_data.R"
       }

        call BenchmarkAndRocCurves.BenchmarkComparison as CHMSampleHeadToHead {
            input:
                evalVcf = CHM_evalVcf,
                evalLabel = CHM_evalLabel,
                evalVcfIndex = CHM_evalVcfIndex,
                evalRuntimeSummaries = CHM_evalRuntimeSummaries,
                evalMonitoringExample = CHM_evalMonitoringExample,

                controlVcf = CHM_controlVcf,
                controlLabel = CHM_controlLabel,
                controlVcfIndex = CHM_controlVcfIndex,
                controlRuntimeSummaries = CHM_controlRuntimeSummaries,
                controlMonitoringExample = CHM_controlMonitoringExample,

                truthVcf = CHM_truthVcf,
                confidenceInterval = CHM_confidenceInterval,
                truthLabel = CHM_truthLabel,
                truthVcfIndex = CHM_truthVcfIndex,

                reference = reference,
                refIndex = refIndex,
                refDict = refDict,
                hapMap = hapMap,
                referenceVersion = referenceVersion,
                stratIntervals = stratIntervals,
                stratLabels = stratLabels,

                monitoring_log_process_script = monitoring_log_process_script
        }

        call BenchmarkAndRocCurves.BenchmarkComparison as NISTSampleHeadToHead {
            input:
                evalVcf = NIST_evalVcf,
                evalLabel = NIST_evalLabel,
                evalVcfIndex = NIST_evalVcfIndex,
                evalRuntimeSummaries = NIST_evalRuntimeSummaries,
                evalMonitoringExample = NIST_evalMonitoringExample,

                controlVcf = NIST_controlVcf,
                controlLabel = NIST_controlLabel,
                controlVcfIndex = NIST_controlVcfIndex,
                controlRuntimeSummaries = NIST_controlRuntimeSummaries,
                controlMonitoringExample = NIST_controlMonitoringExample,

                truthVcf = NIST_truthVcf,
                confidenceInterval = NIST_confidenceInterval,
                truthLabel = NIST_truthLabel,
                truthVcfIndex = NIST_truthVcfIndex,

                reference = reference,
                refIndex = refIndex,
                refDict = refDict,
                hapMap = hapMap,
                referenceVersion = referenceVersion,
                stratIntervals = stratIntervals,
                stratLabels = stratLabels,

                monitoring_log_process_script = monitoring_log_process_script
        }


        call BenchmarkAndRocCurves.BenchmarkComparison as EXOME1SampleHeadToHead {
            input:
                evalVcf = EXOME1_evalVcf,
                evalLabel = EXOME1_evalLabel,
                evalVcfIndex = EXOME1_evalVcfIndex,
                evalRuntimeSummaries = EXOME1_evalRuntimeSummaries,
                evalMonitoringExample = EXOME1_evalMonitoringExample,

                controlVcf = EXOME1_controlVcf,
                controlLabel = EXOME1_controlLabel,
                controlVcfIndex = EXOME1_controlVcfIndex,
                controlRuntimeSummaries = EXOME1_controlRuntimeSummaries,
                controlMonitoringExample = EXOME1_controlMonitoringExample,

                truthVcf = EXOME1_truthVcf,
                confidenceInterval = EXOME1_confidenceInterval,
                truthLabel = EXOME1_truthLabel,
                truthVcfIndex = EXOME1_truthVcfIndex,

                reference = reference,
                refIndex = refIndex,
                refDict = refDict,
                hapMap = hapMap,
                referenceVersion = referenceVersion,

                monitoring_log_process_script = monitoring_log_process_script
        }

        call CreateHTMLReport {
            input:
                merged_chm_plots = CHMSampleHeadToHead.rocPlotsmerged,
                merged_nist_plots = NISTSampleHeadToHead.rocPlotsmerged,
                merged_exome1_plots = EXOME1SampleHeadToHead.rocPlotsmerged,
        }

        output {
            File ROC_Plots_Reported = CreateHTMLReport.report

            File? CHMevalsummary=CHMSampleHeadToHead.evalsummary
            Float CHMevalsnpPrecision=CHMSampleHeadToHead.evalsnpPrecision
            Float CHMevalindelPrecision=CHMSampleHeadToHead.evalindelPrecision
            Float CHMevalsnpRecall=CHMSampleHeadToHead.evalsnpRecall
            Float CHMevalindelRecall=CHMSampleHeadToHead.evalindelRecall
            Float CHMevalsnpF1Score=CHMSampleHeadToHead.evalsnpF1Score
            Float CHMevalindelF1Score=CHMSampleHeadToHead.evalindelF1Score

            File? CHMcontrolsummary=CHMSampleHeadToHead.controlsummary
            Float CHMcontrolsnpPrecision=CHMSampleHeadToHead.controlsnpPrecision
            Float CHMcontrolindelPrecision=CHMSampleHeadToHead.controlindelPrecision
            Float CHMcontrolsnpRecall=CHMSampleHeadToHead.controlsnpRecall
            Float CHMcontrolindelRecall=CHMSampleHeadToHead.controlindelRecall
            Float CHMcontrolsnpF1Score=CHMSampleHeadToHead.controlsnpF1Score
            Float CHMcontrolindelF1Score=CHMSampleHeadToHead.controlindelF1Score

            Float CHMcontrolHCwallclockhours=CHMSampleHeadToHead.controlHCwallclockhours
            Float CHMcontrolHCwallclockmax=CHMSampleHeadToHead.controlHCwallclockmax
            Float CHMcontrolHCprocesshours=CHMSampleHeadToHead.controlHCprocesshours
            Float CHMcontrolHCsystemhours=CHMSampleHeadToHead.controlHCsystemhours
            File CHMcontrolMonitoringLogs=CHMSampleHeadToHead.controlMonitoringLogs

            Float CHMevalHCwallclockhours=CHMSampleHeadToHead.evalHCwallclockhours
            Float CHMevalHCwallclockmax=CHMSampleHeadToHead.evalHCwallclockmax
            Float CHMevalHCprocesshours=CHMSampleHeadToHead.evalHCprocesshours
            Float CHMevalHCsystemhours=CHMSampleHeadToHead.evalHCsystemhours
            File CHMevalMonitoringLogs=CHMSampleHeadToHead.evalMonitoringLogs


            File? NISTevalsummary=NISTSampleHeadToHead.evalsummary
            Float NISTevalsnpPrecision=NISTSampleHeadToHead.evalsnpPrecision
            Float NISTevalindelPrecision=NISTSampleHeadToHead.evalindelPrecision
            Float NISTevalsnpRecall=NISTSampleHeadToHead.evalsnpRecall
            Float NISTevalindelRecall=NISTSampleHeadToHead.evalindelRecall
            Float NISTevalsnpF1Score=NISTSampleHeadToHead.evalsnpF1Score
            Float NISTevalindelF1Score=NISTSampleHeadToHead.evalindelF1Score

            File? NISTcontrolsummary=NISTSampleHeadToHead.controlsummary
            Float NISTcontrolsnpPrecision=NISTSampleHeadToHead.controlsnpPrecision
            Float NISTcontrolindelPrecision=NISTSampleHeadToHead.controlindelPrecision
            Float NISTcontrolsnpRecall=NISTSampleHeadToHead.controlsnpRecall
            Float NISTcontrolindelRecall=NISTSampleHeadToHead.controlindelRecall
            Float NISTcontrolsnpF1Score=NISTSampleHeadToHead.controlsnpF1Score
            Float NISTcontrolindelF1Score=NISTSampleHeadToHead.controlindelF1Score

            Float NISTcontrolHCwallclockhours=NISTSampleHeadToHead.controlHCwallclockhours
            Float NISTcontrolHCwallclockmax=NISTSampleHeadToHead.controlHCwallclockmax
            Float NISTcontrolHCprocesshours=NISTSampleHeadToHead.controlHCprocesshours
            Float NISTcontrolHCsystemhours=NISTSampleHeadToHead.controlHCsystemhours
            File NISTcontrolMonitoringLogs=NISTSampleHeadToHead.controlMonitoringLogs

            Float NISTevalHCwallclockhours=NISTSampleHeadToHead.evalHCwallclockhours
            Float NISTevalHCwallclockmax=NISTSampleHeadToHead.evalHCwallclockmax
            Float NISTevalHCprocesshours=NISTSampleHeadToHead.evalHCprocesshours
            Float NISTevalHCsystemhours=NISTSampleHeadToHead.evalHCsystemhours
            File NISTevalMonitoringLogs=NISTSampleHeadToHead.evalMonitoringLogs



            File? EXOME1evalsummary=EXOME1SampleHeadToHead.evalsummary
            Float EXOME1evalsnpPrecision=EXOME1SampleHeadToHead.evalsnpPrecision
            Float EXOME1evalindelPrecision=EXOME1SampleHeadToHead.evalindelPrecision
            Float EXOME1evalsnpRecall=EXOME1SampleHeadToHead.evalsnpRecall
            Float EXOME1evalindelRecall=EXOME1SampleHeadToHead.evalindelRecall
            Float EXOME1evalsnpF1Score=EXOME1SampleHeadToHead.evalsnpF1Score
            Float EXOME1evalindelF1Score=EXOME1SampleHeadToHead.evalindelF1Score

            File? EXOME1controlsummary=EXOME1SampleHeadToHead.controlsummary
            Float EXOME1controlsnpPrecision=EXOME1SampleHeadToHead.controlsnpPrecision
            Float EXOME1controlindelPrecision=EXOME1SampleHeadToHead.controlindelPrecision
            Float EXOME1controlsnpRecall=EXOME1SampleHeadToHead.controlsnpRecall
            Float EXOME1controlindelRecall=EXOME1SampleHeadToHead.controlindelRecall
            Float EXOME1controlsnpF1Score=EXOME1SampleHeadToHead.controlsnpF1Score
            Float EXOME1controlindelF1Score=EXOME1SampleHeadToHead.controlindelF1Score

          }
}


task CreateHTMLReport {
    input {
        File merged_chm_plots
        File merged_nist_plots
        File merged_exome1_plots
        String additional_label = ""
    }

    command <<<
        set -xeuo pipefail

        source activate fe_evaluation

        chm_plots_base64=$(base64 -w 0 ~{merged_chm_plots})
        nist_plots_base64=$(base64 -w 0 ~{merged_nist_plots})
        exome1_plots_base64=$(base64 -w 0 ~{merged_exome1_plots})

        cat <<EOF > report.html
<!DOCTYPE html>
<html>
    <head>
        <meta charset="UTF-8">
        <title>HaplotypeCaller CARROT Test Report ~{additional_label}</title>
        <style>
            body {
                font-family: sans-serif;
                font-size: 14px;
                padding: 0 26px;
                line-height: 1.6;
            }
            img {
                max-width: 100%;
                max-height: 100%;
            }
            table {
                border-collapse: collapse;
            }
            th, td {
                padding: 5px 10px;
            }
            table, td {
                border: 1px solid black;
            }
        </style>
    </head>
    <body>
        <h2>CHM plots</h2>
        <table>
            <tr>
                <th style="text-align: center;">~{additional_label}</th>
            </tr>
            <tr>
                <td><img src="data:image/png;base64,$chm_plots_base64" /></td>
            </tr>
        </table>
        <h2>HGOO2 plots</h2>
        <table>
            <tr>
                <th style="text-align: center;">~{additional_label}</th>
            </tr>
            <tr>
                <td><img src="data:image/png;base64,$nist_plots_base64" /></td>
            </tr>
        </table>
        <h2>NA12878 Twist Exome plots</h2>
        <table>
            <tr>
                <th style="text-align: center;">~{additional_label}</th>
            </tr>
            <tr>
                <td><img src="data:image/png;base64,$exome1_plots_base64" /></td>
            </tr>
        </table>
    </body>
</html>
EOF
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/functionalequivalence/fe_evaluation:1.0.0"
        memory: "2 GB"
        disks: "local-disk 20 HDD"
        maxRetries: 3
    }

    output {
        File report = "report.html"
    }
}