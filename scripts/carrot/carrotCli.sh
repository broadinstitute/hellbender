#!/bin/bash

#carrot_cli software create --name example-test --description "An example CARROT test software" --repository_url https://github.com/broadinstitute/carrot-example-test

# report ROC plot file
carrot_cli result create --name "ROC_Plots_Reported" --description "Report file for ROC plots" --result_type file

#CHM outputs
carrot_cli result create --name "CHM evalsummary" --description "VCF Eval output summary file" --result_type file
carrot_cli result create --name "CHM evalsnpPrecision" --description "VCF eval snp precision at max f1 point for CHM sample" --result_type numeric
carrot_cli result create --name "CHM evalindelPrecision" --description "VCF eval indel precision at max f1 point for CHM sample" --result_type numeric
carrot_cli result create --name "CHM evalsnpRecall" --description "VCF eval snp recall at max f1 point for CHM sample" --result_type numeric
carrot_cli result create --name "CHM evalsnpF1Score" --description "VCF eval snp F1Score at max f1 point for CHM sample" --result_type numeric
carrot_cli result create --name "CHM evalindelF1Score" --description "VCF eval indel F1Score at max f1 point for CHM sample" --result_type numeric
carrot_cli result create --name "CHM controlsummary" --description "VCF control output summary file" --result_type file
carrot_cli result create --name "CHM controlsnpPrecision" --description "VCF control snp precision at max f1 point for CHM sample" --result_type numeric
carrot_cli result create --name "CHM controlindelPrecision" --description "VCF control indel precision at max f1 point for CHM sample" --result_type numeric
carrot_cli result create --name "CHM controlsnpRecall" --description "VCF control snp recall at max f1 point for CHM sample" --result_type numeric
carrot_cli result create --name "CHM controlsnpF1Score" --description "VCF control snp F1Score at max f1 point for CHM sample" --result_type numeric
carrot_cli result create --name "CHM controlindelF1Score" --description "VCF control indel F1Score at max f1 point for CHM sample" --result_type numeric
# wallclocks
carrot_cli result create --name "CHM evalHCwallclockhours" --description "Total Haplotype caller wallclock hours for eval CHM" --result_type numeric
carrot_cli result create --name "CHM evalHCwallclockmax" --description "Maximum shard Haplotype caller wallclock hours for eval CHM" --result_type numeric
carrot_cli result create --name "CHM evalHCprocesshours" --description "Total Haplotype caller processor cpu hours for eval CHM" --result_type numeric
carrot_cli result create --name "CHM evalHCsystemhours" --description "Total Haplotype caller system cpu hours for eval CHM" --result_type numeric
carrot_cli result create --name "CHM evalMonitoringLogs" --description "PDF file with monitoring logs outputs from representative eval shard for CHM" --result_type file
carrot_cli result create --name "CHM controlHCwallclockhours" --description "Total Haplotype caller wallclock hours for control CHM" --result_type numeric
carrot_cli result create --name "CHM controlHCwallclockmax" --description "Maximum shard Haplotype caller wallclock hours for control CHM" --result_type numeric
carrot_cli result create --name "CHM controlHCprocesshours" --description "Total Haplotype caller processor cpu hours for control CHM" --result_type numeric
carrot_cli result create --name "CHM controlHCsystemhours" --description "Total Haplotype caller system cpu hours for control CHM" --result_type numeric
carrot_cli result create --name "CHM controlMonitoringLogs" --description "PDF file with monitoring logs outputs from representative control shard for CHM" --result_type file

#HG002 outputs
carrot_cli result create --name "NIST evalsummary" --description "VCF Eval output summary file" --result_type file
carrot_cli result create --name "NIST evalsnpPrecision" --description "VCF eval snp precision at max f1 point for NIST sample" --result_type numeric
carrot_cli result create --name "NIST evalindelPrecision" --description "VCF eval indel precision at max f1 point for NIST sample" --result_type numeric
carrot_cli result create --name "NIST evalsnpRecall" --description "VCF eval snp recall at max f1 point for NIST sample" --result_type numeric
carrot_cli result create --name "NIST evalsnpF1Score" --description "VCF eval snp F1Score at max f1 point for NIST sample" --result_type numeric
carrot_cli result create --name "NIST evalindelF1Score" --description "VCF eval indel F1Score at max f1 point for NIST sample" --result_type numeric
carrot_cli result create --name "NIST controlsummary" --description "VCF control output summary file" --result_type file
carrot_cli result create --name "NIST controlsnpPrecision" --description "VCF control snp precision at max f1 point for NIST sample" --result_type numeric
carrot_cli result create --name "NIST controlindelPrecision" --description "VCF control indel precision at max f1 point for NIST sample" --result_type numeric
carrot_cli result create --name "NIST controlsnpRecall" --description "VCF control snp recall at max f1 point for NIST sample" --result_type numeric
carrot_cli result create --name "NIST controlsnpF1Score" --description "VCF control snp F1Score at max f1 point for NIST sample" --result_type numeric
carrot_cli result create --name "NIST controlindelF1Score" --description "VCF control indel F1Score at max f1 point for NIST sample" --result_type numeric
# wallclocks
carrot_cli result create --name "NIST evalHCwallclockhours" --description "Total Haplotype caller wallclock hours for eval NIST" --result_type file
carrot_cli result create --name "NIST evalHCwallclockmax" --description "Maximum shard Haplotype caller wallclock hours for eval NIST" --result_type numeric
carrot_cli result create --name "NIST evalHCprocesshours" --description "Total Haplotype caller processor cpu hours for eval NIST" --result_type numeric
carrot_cli result create --name "NIST evalHCsystemhours" --description "Total Haplotype caller system cpu hours for eval NIST" --result_type numeric
carrot_cli result create --name "NIST evalMonitoringLogs" --description "PDF file with monitoring logs outputs from representative eval shard for NIST" --result_type file
carrot_cli result create --name "NIST controlHCwallclockhours" --description "Total Haplotype caller wallclock hours for control NIST" --result_type numeric
carrot_cli result create --name "NIST controlHCwallclockmax" --description "Maximum shard Haplotype caller wallclock hours for control NIST" --result_type numeric
carrot_cli result create --name "NIST controlHCprocesshours" --description "Total Haplotype caller processor cpu hours for control NIST" --result_type numeric
carrot_cli result create --name "NIST controlHCsystemhours" --description "Total Haplotype caller system cpu hours for control NIST" --result_type numeric
carrot_cli result create --name "NIST controlMonitoringLogs" --description "PDF file with monitoring logs outputs from representative control shard for NIST" --result_type file

#Exome outputs
carrot_cli result create --name "EXOME1 evalsummary" --description "VCF Eval output summary file" --result_type file
carrot_cli result create --name "EXOME1 evalsnpPrecision" --description "VCF eval snp precision at max f1 point for EXOME1 sample" --result_type numeric
carrot_cli result create --name "EXOME1 evalindelPrecision" --description "VCF eval indel precision at max f1 point for EXOME1 sample" --result_type numeric
carrot_cli result create --name "EXOME1 evalsnpRecall" --description "VCF eval snp recall at max f1 point for EXOME1 sample" --result_type numeric
carrot_cli result create --name "EXOME1 evalsnpF1Score" --description "VCF eval snp F1Score at max f1 point for EXOME1 sample" --result_type numeric
carrot_cli result create --name "EXOME1 evalindelF1Score" --description "VCF eval indel F1Score at max f1 point for EXOME1 sample" --result_type numeric
carrot_cli result create --name "EXOME1 controlsummary" --description "VCF control output summary file" --result_type file
carrot_cli result create --name "EXOME1 controlsnpPrecision" --description "VCF control snp precision at max f1 point for EXOME1 sample" --result_type numeric
carrot_cli result create --name "EXOME1 controlindelPrecision" --description "VCF control indel precision at max f1 point for EXOME1 sample" --result_type numeric
carrot_cli result create --name "EXOME1 controlsnpRecall" --description "VCF control snp recall at max f1 point for EXOME1 sample" --result_type numeric
carrot_cli result create --name "EXOME1 controlsnpF1Score" --description "VCF control snp F1Score at max f1 point for EXOME1 sample" --result_type numeric
carrot_cli result create --name "EXOME1 controlindelF1Score" --description "VCF control indel F1Score at max f1 point for EXOME1 sample" --result_type numeric

#HaplotypeCaller
carrot_cli pipeline create --name "HaplotypeCaller_eval_on_test_samples" --description "Pipeline for testing for results regressions in HaplotypeCaller on standard samples"

carrot_cli template create --pipeline "HaplotypeCaller_eval_on_test_samples" --name "HaplotypeCaller orchestrated template" --test_wdl "gs://dsp-methods-carrot-data/test_data/haplotypecaller_tests/wdls/HaplotypeCallerComparisonForCarrot.SamplesOrchestrated.wdl" --eval_wdl "gs://dsp-methods-carrot-data/test_data/haplotypecaller_tests/wdls/BenchmarkComparisonOrchestrated.wdl"

## BIND THE RESULTS TO OUR TEMPLATE
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "CHM evalsummary" BenchmarkVCFsHeadToHeadOrchestrated.CHMevalsummary
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "CHM evalsnpPrecision" BenchmarkVCFsHeadToHeadOrchestrated.CHMevalsnpPrecision
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "CHM evalindelPrecision" BenchmarkVCFsHeadToHeadOrchestrated.CHMevalindelPrecision
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "CHM evalsnpRecall" BenchmarkVCFsHeadToHeadOrchestrated.CHMevalsnpRecall
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "CHM evalsnpF1Score" BenchmarkVCFsHeadToHeadOrchestrated.CHMevalsnpF1Score
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "CHM evalindelF1Score" BenchmarkVCFsHeadToHeadOrchestrated.CHMevalindelF1Score
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "CHM controlsummary" BenchmarkVCFsHeadToHeadOrchestrated.CHMcontrolsummary
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "CHM controlsnpPrecision" BenchmarkVCFsHeadToHeadOrchestrated.CHMcontrolsnpPrecision
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "CHM controlindelPrecision" BenchmarkVCFsHeadToHeadOrchestrated.CHMcontrolindelPrecision
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "CHM controlsnpRecall" BenchmarkVCFsHeadToHeadOrchestrated.CHMcontrolsnpRecall
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "CHM controlsnpF1Score" BenchmarkVCFsHeadToHeadOrchestrated.CHMcontrolsnpF1Score
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "CHM controlindelF1Score" BenchmarkVCFsHeadToHeadOrchestrated.CHMcontrolindelF1Score
# wallclocks
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "CHM evalHCwallclockhours" BenchmarkVCFsHeadToHeadOrchestrated.CHMevalHCwallclockhours
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "CHM evalHCwallclockmax" BenchmarkVCFsHeadToHeadOrchestrated.CHMevalHCwallclockmax
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "CHM evalHCprocesshours" BenchmarkVCFsHeadToHeadOrchestrated.CHMevalHCprocesshours
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "CHM evalHCsystemhours" BenchmarkVCFsHeadToHeadOrchestrated.CHMevalHCsystemhours
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "CHM evalMonitoringLogs" BenchmarkVCFsHeadToHeadOrchestrated.CHMevalMonitoringLogs
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "CHM controlHCwallclockhours" BenchmarkVCFsHeadToHeadOrchestrated.CHMcontrolHCwallclockhours
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "CHM controlHCwallclockmax" BenchmarkVCFsHeadToHeadOrchestrated.CHMcontrolHCwallclockmax
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "CHM controlHCprocesshours" BenchmarkVCFsHeadToHeadOrchestrated.CHMcontrolHCprocesshours
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "CHM controlHCsystemhours" BenchmarkVCFsHeadToHeadOrchestrated.CHMcontrolHCsystemhours
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "CHM controlMonitoringLogs" BenchmarkVCFsHeadToHeadOrchestrated.CHMcontrolMonitoringLogs

#HG002 outputs
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "NIST evalsummary" BenchmarkVCFsHeadToHeadOrchestrated.NISTevalsummary
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "NIST evalsnpPrecision" BenchmarkVCFsHeadToHeadOrchestrated.NISTevalsnpPrecision
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "NIST evalindelPrecision" BenchmarkVCFsHeadToHeadOrchestrated.NISTevalindelPrecision
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "NIST evalsnpRecall" BenchmarkVCFsHeadToHeadOrchestrated.NISTevalsnpRecall
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "NIST evalsnpF1Score" BenchmarkVCFsHeadToHeadOrchestrated.NISTevalsnpF1Score
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "NIST evalindelF1Score" BenchmarkVCFsHeadToHeadOrchestrated.NISTevalindelF1Score
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "NIST controlsummary" BenchmarkVCFsHeadToHeadOrchestrated.NISTcontrolsummary
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "NIST controlsnpPrecision" BenchmarkVCFsHeadToHeadOrchestrated.NISTcontrolsnpPrecision
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "NIST controlindelPrecision" BenchmarkVCFsHeadToHeadOrchestrated.NISTcontrolindelPrecision
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "NIST controlsnpRecall" BenchmarkVCFsHeadToHeadOrchestrated.NISTcontrolsnpRecall
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "NIST controlsnpF1Score" BenchmarkVCFsHeadToHeadOrchestrated.NISTcontrolsnpF1Score
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "NIST controlindelF1Score" BenchmarkVCFsHeadToHeadOrchestrated.NISTcontrolindelF1Score
# wallclocks
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "NIST evalHCwallclockhours" BenchmarkVCFsHeadToHeadOrchestrated.NISTevalHCwallclockhours
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "NIST evalHCwallclockmax" BenchmarkVCFsHeadToHeadOrchestrated.NISTevalHCwallclockmax
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "NIST evalHCprocesshours" BenchmarkVCFsHeadToHeadOrchestrated.NISTevalHCprocesshours
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "NIST evalHCsystemhours" BenchmarkVCFsHeadToHeadOrchestrated.NISTevalHCsystemhours
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "NIST evalMonitoringLogs" BenchmarkVCFsHeadToHeadOrchestrated.NISTevalMonitoringLogs
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "NIST controlHCwallclockhours" BenchmarkVCFsHeadToHeadOrchestrated.NISTcontrolHCwallclockhours
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "NIST controlHCwallclockmax" BenchmarkVCFsHeadToHeadOrchestrated.NISTcontrolHCwallclockmax
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "NIST controlHCprocesshours" BenchmarkVCFsHeadToHeadOrchestrated.NISTcontrolHCprocesshours
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "NIST controlHCsystemhours" BenchmarkVCFsHeadToHeadOrchestrated.NISTcontrolHCsystemhours
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "NIST controlMonitoringLogs" BenchmarkVCFsHeadToHeadOrchestrated.NISTcontrolMonitoringLogs

#Exome outputs
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "EXOME1 evalsummary" BenchmarkVCFsHeadToHeadOrchestrated.EXOME1evalsummary
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "EXOME1 evalsnpPrecision" BenchmarkVCFsHeadToHeadOrchestrated.EXOME1evalsnpPrecision
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "EXOME1 evalindelPrecision" BenchmarkVCFsHeadToHeadOrchestrated.EXOME1evalindelPrecision
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "EXOME1 evalsnpRecall" BenchmarkVCFsHeadToHeadOrchestrated.EXOME1evalsnpRecall
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "EXOME1 evalsnpF1Score" BenchmarkVCFsHeadToHeadOrchestrated.EXOME1evalsnpF1Score
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "EXOME1 evalindelF1Score" BenchmarkVCFsHeadToHeadOrchestrated.EXOME1evalindelF1Score
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "EXOME1 controlsummary" BenchmarkVCFsHeadToHeadOrchestrated.EXOME1controlsummary
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "EXOME1 controlsnpPrecision" BenchmarkVCFsHeadToHeadOrchestrated.EXOME1controlsnpPrecision
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "EXOME1 controlindelPrecision" BenchmarkVCFsHeadToHeadOrchestrated.EXOME1controlindelPrecision
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "EXOME1 controlsnpRecall" BenchmarkVCFsHeadToHeadOrchestrated.EXOME1controlsnpRecall
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "EXOME1 controlsnpF1Score" BenchmarkVCFsHeadToHeadOrchestrated.EXOME1controlsnpF1Score
carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "EXOME1 controlindelF1Score" BenchmarkVCFsHeadToHeadOrchestrated.EXOME1controlindelF1Score

carrot_cli template map_to_result "HaplotypeCaller orchestrated template" "ROC_Plots_Reported" BenchmarkVCFsHeadToHeadOrchestrated.ROC_Plots_Reported


carrot_cli test create --name "HaplotypeCaller CARROT Regression Tests" --template "HaplotypeCaller orchestrated template" --description "Standard set of scientific tests to check for haplotype caller releases" --test_input_defaults "HC_carrot_test_default_inputs.json" --eval_input_defaults "HC_carrot_eval_default_inputs.json"