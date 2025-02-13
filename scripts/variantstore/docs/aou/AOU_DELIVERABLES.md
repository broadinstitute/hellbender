# Running the Genome Variant Store (GVS) Pipelines for AoU

## Setup
- Create a Terra workspace
  - Request a new AoU workspace as described in Section 3.D. of [AoU DRC Protocols](https://docs.google.com/document/d/1ooK0wbHLgSueiepjTTyLLI6Zz7vi1-GhKTjmCd8ZHwU/edit?usp=sharing).
  - As described in the "Getting Started" of [Operational concerns for running Hail in Terra Cromwell/WDL](https://docs.google.com/document/d/1_OY2rKwZ-qKCDldSZrte4jRIZf4eAw2d7Jd-Asi50KE/edit?usp=sharing), this workspace will need permission in Terra to run Hail dataproc clusters within WDL. Contact Emily to request this access as part of setting up the new workspace.
  - There is a quota that needs to be upgraded for the process of Bulk Ingest.
    When we ingest data, we use the Write API, which is part of BQ’s Storage API. Since we are hitting this API with so much data all at once, we want to increase our CreateWriteStream quota. Follow the [Quota Request Template](workspace/CreateWriteStreamRequestIncreasedQuota.md).
    Once that quota has been increased, the `load_data_scatter_width` value needs to be updated based on that new quota (for information on what we did for Echo, see the "Calculate Quota To be Requested" section in the [Quota Request Template](workspace/CreateWriteStreamRequestIncreasedQuota.md) doc).
  - Create and push a feature branch (e.g. `EchoCallset`) based off the `ah_var_store` branch to the GATK GitHub repo.
    - Update the .dockstore.yml file on that feature branch to add the feature branch for all the WDLs that will be loaded into the workspace in the next step.
- Once the requested workspace has been created and permissioned, populate with the following WDLs:
  - [GvsBulkIngestGenomes](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsBulkIngestGenomes) workflow
  - [GvsWithdrawSamples](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsWithdrawSamples) workflow
  - [GvsPopulateAltAllele](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsPopulateAltAllele) workflow
  - [GvsCreateFilterSet](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsCreateFilterSet) workflow
  - [GvsExtractAvroFilesForHail](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsExtractAvroFilesForHail) workflow
  - [GvsPrepareRangesCallset](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsPrepareRangesCallset) workflow
  - [GvsExtractCallset](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsExtractCallset) workflow
  - [GvsExtractCallsetPgenMerged](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsExtractCallsetPgenMerged) workflow
  - [GvsCallsetStatistics](https://dockstore.org/workflows/github.com/broadinstitute/gatk/GvsCallsetStatistics) workflow
  - [GvsCalculatePrecisionAndSensitivity](https://dockstore.org/workflows/github.com/broadinstitute/gatk/GvsCalculatePrecisionAndSensitivity) workflow
  - [GvsCallsetCost](https://dockstore.org/workflows/github.com/broadinstitute/gatk/GvsCallsetCost) workflow
  - [GvsCreateVDS](https://dockstore.org/workflows/github.com/broadinstitute/gatk/GvsCreateVDS) workflow
  - [GvsCreateVATfromVDS](https://dockstore.org/workflows/github.com/broadinstitute/gatk/GvsCreateVATfromVDS) workflow
  - [GvsValidateVat](https://dockstore.org/my-workflows/github.com/broadinstitute/gatk/GvsValidateVat) workflow
- Install the [Fetch WGS metadata for samples from list](./workspace/Fetch%20WGS%20metadata%20for%20samples%20from%20list.ipynb) python notebook in the workspace that has been created.
  - Place the file with the list of the new samples to ingest in a GCS location the notebook (running with your @pmi-ops account) will have access to.  This will grab the samples from the workspace where they were reblocked and bring them into this callset workspace.
  - Run the steps in the notebook:
      - In section 1.2, Set the `sample_list_file_path` variable in that notebook to the path of the file
      - Run the cells up to and through, section 1.4 ("Perform the copy, in batches")
      - If you want to automatically break up the new samples into smaller sample sets, then run the "now that the data have been copied, you can make sample sets if you wish" step. Set the `SUBSET_SIZE` and `set_name` variables to customize.
- **NOTE** If you want to create a large sample set after you have run the notebook, Terra provides (and recommends you use) this python [script](https://github.com/broadinstitute/firecloud-tools/tree/master/scripts/import_large_tsv) which allows you to upload a sample set to the workspace.
- Create a dataset in the Google project. Make sure that when you are creating the dataset that you set the `location type` to be `Multi-Region`.
- Make a note of the Google project ID ("aou-genomics-curation-prod"), dataset name (e.g. "aou_wgs" — if it does not exist be sure to create one before running any workflows) and callset identifier (e.g. "echo") as these will be inputs (`project_id`, `dataset_name` and `call_set_identifier`) to all or most of the GVS workflows. The [naming conventions for other aspects of GVS datasets are outlined here](https://docs.google.com/document/d/1pNtuv7uDoiOFPbwe4zx5sAGH7MyxwKqXkyrpNmBxeow).
- Once the **non-control** samples have been fully ingested into BQ using the `GvsBulkIngestGenomes` workflow, the **control** samples can be manually added to the workspace and loaded in separately

## The Main Pipeline
1. `GvsBulkIngestGenomes` workflow
   - For use with **non-control** samples only! To ingest control samples (required for running `GvsCalculatePrecisionAndSensitivity`), use the`GvsAssignIds` and `GvsImportGenomes` workflows described below.
   - Set `sample_id_column_name` to "research_id" to use the shorter unique ID from AoU for the `sample_name` column.
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
   - This workflow will be run twice: first to load only VCF headers for validation purposes, then a second time to load variant and reference data.
   1. `GvsBulkIngestGenomes` header ingest and validation
      - Set `load_vcf_headers` to `true` and `load_vet_and_ref_ranges` to `false` to load VCF header data only.
      - Note that when run with these parameters the workflow is expected to call `TerminateWorkflow` with a message that header data has been successfully loaded. In the Terra UI this will look like a failure but in this particular case is the desired outcome.
      - Once these headers have been loaded, run a sanity checking query on the DRAGEN version ```
  SELECT
  REGEXP_EXTRACT(vcf_header_lines, r'SW: [0-9\.]+') AS version,
  COUNT(*)
FROM
  `aou-genomics-curation-prod.echo.vcf_header_lines` AS header_lines
WHERE
  header_lines.is_expected_unique = TRUE
  AND CONTAINS_SUBSTR(vcf_header_lines, 'DRAGENCommandLine=<ID=dragen,')
GROUP BY
  version```. The version string here appears to be a mix of hardware and software versions. What matters for us is that the last triplet is `3.7.8`. In the Echo callset this query currently returns two rows with `version` values of `SW: 05.021.604.3.7.8` and `SW: 07.021.604.3.7.8`. Assuming this query returns only rows with `3.7.8` as the final triplet, proceed with the second invocation of `GvsBulkIngestGenomes` documented below.
   1. `GvsBulkIngestGenomes` variant and reference data ingest
      - If and only if the header ingest described above completed successfully, proceed with the loading of variant and reference data.
      - Set `load_vcf_headers` to `false` and `load_vet_and_ref_ranges` to `true` to load variant and reference data.
      - **NOTE** Be sure to set the input `drop_state` to ZERO (this will have the effect of dropping GQ0 reference blocks) and set `use_compressed_references` to true (this will further compress the reference data).
   - Note: In case of mistakenly ingesting a large number of bad samples, instructions for removing them can be found in [this Jira ticket](https://broadworkbench.atlassian.net/browse/VS-1206)
1. `GvsWithdrawSamples` workflow
   - Run if there are any samples to withdraw from the last callset.
   - When you run the `GvsWithdrawSamples` workflow, you should inspect the output of the workflow.
     - The output `num_samples_withdrawn` indicates the number of samples that have been withdrawn. This number should agree with that which you expect.
     - If the workflow fails, it may have failed if the list of samples that was supplied to it includes samples that have not yet been ingested. To determine if this is the case, inspect the output (STDOUT) of the workflow and if it includes a list of samples that need to be ingested, then do so (or investigate the discrepancy). Note that there is a boolean variable `allow_uningested_samples` for this workflow that will allow it to pass if this condition occurs.
1. `GvsPopulateAltAllele` workflow
   - This step loads data into the `alt_allele` table from the `vet_*` tables (which were populated during ingest) in preparation for running the filtering step.
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
1. `GvsCreateFilterSet` workflow
   - This step calculates features from the `alt_allele` table, and trains the VETS filtering model along with site-level QC filters and loads them into BigQuery into a series of `filter_set_*` tables.
   - See [naming conventions doc](https://docs.google.com/document/d/1pNtuv7uDoiOFPbwe4zx5sAGH7MyxwKqXkyrpNmBxeow) for guidance on what to use for `filter_set_name`, which you will need to keep track of for the `GvsExtractAvroFilesForHail` WDL. If, for some reason, this step needs to be run multiple times, be sure to use a different `filter_set_name` (the doc has guidance for this, as well).
   - Run this workflow with the workflow option "Retry with more memory" and choose a "Memory retry factor" of 1.5
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
1. `GvsExtractAvroFilesForHail` workflow
   - This workflow extracts the data in BigQuery and transforms it into Avro files in a Google bucket, incorporating the VETS filter set data.
   - The extracted Avro files will then be used as an input for `GvsCreateVDS` workflow described below.
   - This workflow needs to be run with the `filter_set_name` from `GvsCreateFilterSet` step.
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
1. `GvsCreateVDS` workflow
   - This step creates a VDS based on the Avro files generated from the `GvsExtractAvroFilesForHail` workflow above.
   - You can find what the `avro_path` input should be by going to the `GvsExtractAvroFilesForHail` run in Job Manager; the output `avro_path` is the location of the files created by that workflow.
   - The `vds_path` path input to this workflow represents the output path for the VDS. VDSes should be written under the AoU delivery bucket `gs://prod-drc-broad/`. Ask Lee for the exact path to use for the VDS in the `#dsp-variants` slack channel.
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
   - Once a VDS has been created the Variants team will also generate callset statistics using `GvsCallsetStatistics` as described below. The Variants team then forwards both the path to the VDS and the output callset statistics TSV to Lee to quality check the VDS.
   - If you are debugging a Hail-related issue, you may want to set `leave_hail_cluster_running_at_end` to `true` and refer to [the suggestions for debugging issues with Hail](HAIL_DEBUGGING.md).
1. `GvsCallsetStatistics` workflow
    - You will need to run `GvsPrepareRangesCallset` workflow for callset statistics first, if it has not been run already
       - This workflow transforms the data in the vet tables into a schema optimized for callset stats creation and for calculating sensitivity and precision.
       - The `only_output_vet_tables` input should be set to `true` (the default value is `false`).
       - The `enable_extract_table_ttl` input should be set to `true` (the default value is `false`), which will add a TTL of two weeks to the tables it creates.
       - See [naming conventions doc](https://docs.google.com/document/d/1pNtuv7uDoiOFPbwe4zx5sAGH7MyxwKqXkyrpNmBxeow) for guidance on what to use for `extract_table_prefix` or cohort prefix, which you will need to keep track of for the callset stats.
       - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
    - You will need to have the "BigQuery Data Viewer" role for your @pmi-ops proxy group on the `spec-ops-aou:gvs_public_reference_data.gnomad_v3_sites` table
    - This workflow needs to be run with the `extract_table_prefix` input from `GvsPrepareRangesCallset` step.
    - This workflow needs to be run with the `filter_set_name` input from `GvsCreateFilterSet` step.
    - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
1. `GvsCalculatePrecisionAndSensitivity` workflow
    - Please see the detailed instructions for running the Precision and Sensitivity workflow [here](../../tieout/AoU_PRECISION_SENSITIVITY.md).
1. `GvsCallsetCost` workflow
    - This workflow calculates the total BigQuery cost of generating this callset (which is not represented in the Terra UI total workflow cost) using the above GVS workflows; it's used to calculate the cost as a whole and by sample.


## Internal sign-off protocol

The Variants team currently has the following VDS internal sign-off protocol:

1. Generate a VDS for the candidate callset into the "delivery" bucket.
1. Open up the VDS in a [beefy](vds/cluster/AoU%20VDS%20Cluster%20Configuration.md) notebook and confirm the "shape" looks right.
1. Run `GvsPrepareRangesCallset.wdl` to generate a prepare table of VET data.
1. Run `GvsCallsetStatistics.wdl` to generate callset statistics for the candidate callset using the prepare VET table created in the preceding step.
1. Copy the output of `GvsCallsetStatistics.wdl` into the "delivery" bucket.
1. Email the paths to the VDS and callset statistics to Lee/Wail for QA / approval.


## Main Deliverables (via email to stakeholders once the above steps are complete)

The Callset Stats and S&P files can be simply `gsutil cp`'ed to the AoU delivery bucket since they are so much smaller.
1. GCS location of the VDS in the AoU delivery bucket
1. Fully qualified name of the BigQuery dataset (composed of the `project_id` and `dataset_name` inputs from the workflows)
1. GCS location of the CSV output from `GvsCallsetStatistics` workflow in the AoU delivery bucket
1. GCS location of the TSV output from `GvsCalculatePrecisionAndSensitivity` in the AoU delivery bucket

## Running the VAT pipeline
To create a BigQuery table of variant annotations, you may follow the instructions here:
[process to create variant annotations table](../../variant-annotations-table/README.md)
The pipeline takes in the VDS and outputs a variant annotations table in BigQuery.

Once the VAT table is created and a tsv is exported, the AoU research workbench team should be notified of its creation and permission should be granted so that several members of the team have view permission.

- Grant `BigQuery Data Viewer` permission to specific people's PMI-OPS accounts. This will include members of the AoU research workbench team.
- Copy the tarred and bgzipped export of the VAT into the pre-delivery bucket.
- Send an email out notifying the AoU research workbench team of the readiness of the VAT. Additionally, a RW Jira ticket will be made by project management to request copying the VAT to pre-prod.
- A document describing how this information was shared (for previous callsets) is located [here](https://docs.google.com/document/d/1caqgCS1b_dDJXQT4L-tRxjOkLGDgRNkO9eac1xd9ib0/edit)

## Additional Deliverables

### Smaller Interval Lists

1. You will need to run the `GvsPrepareRangesCallset` workflow for each "[Region](https://support.researchallofus.org/hc/en-us/articles/14929793660948-Smaller-Callsets-for-Analyzing-Short-Read-WGS-SNP-Indel-Data-with-Hail-MT-VCF-and-PLINK)" (interval list) for which a PGEN or VCF deliverable is required for the callset.
   - This workflow transforms the data in the vet, ref_ranges, and samples tables into a schema optimized for extract.
   - The `enable_extract_table_ttl` input should be set to `true` (the default value is `false`), which will add a TTL of two weeks to the tables it creates.
   - `extract_table_prefix` should be set to a name that is unique to the given Region / interval list. See the [naming conventions doc](https://docs.google.com/document/d/1pNtuv7uDoiOFPbwe4zx5sAGH7MyxwKqXkyrpNmBxeow) for guidance on what to use.
   - Specify the `interval_list` appropriate for the PGEN / VCF extraction run you are performing.
   - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
1. `GvsExtractCallset` / `GvsExtractCallsetPgenMerged` workflows ("small callset" Exome, Clinvar, and ACAF threshold extracts in VCF and PGEN formats respectively)
    - Specify the same `call_set_identifier`, `dataset_name`, `project_id`, `extract_table_prefix`, and `interval_list` that were used in the `GvsPrepareRangesCallset` run documented above.
    - Specify the `interval_weights_bed` appropriate for the PGEN / VCF extraction run you are performing. `gs://gvs_quickstart_storage/weights/gvs_full_vet_weights_1kb_padded_orig.bed` is the interval weights BED used for Quickstart.
    - For both `GvsExtractCallset` and `GvsExtractCallsetPgenMerged`, select the workflow option "Retry with more memory" and choose a "Memory retry factor" of 1.5
    - `GvsExtractCallset` currently defaults to 1000 alt alleles maximum, which means that any sites having more than that number of alt alleles will be dropped. For AoU callsets make sure to specify the AoU-appropriate `maximum_alternate_alleles` value, currently 100.
    - `GvsExtractCallsetPgen` currently defaults to 100 alt alleles maximum, which means that any sites having more than that number of alt alleles will be dropped.
    - For both `GvsExtractCallset` and `GvsExtractCallsetPgenMerged`, be sure to set the `output_gcs_dir` to the proper path in the AoU delivery bucket so you don't need to copy them there yourself once the workflows have finished.
    - For `GvsExtractCallset`, you will probably (check the requirements to confirm) want to set the input `bgzip_output_vcfs` to `true`.
    - For `GvsExtractCallset`, make sure to set the `ploidy_table_name` to 'sample_chromosome_ploidy' (the default ploidy table name created during ingest)
    - If the overall workflow appears to fail due to angry cloud issues in `PgenExtractTask` (PGEN) or `ExtractTask` (VCF), you can re-run the workflow with the exact same inputs with call caching turned on; the successful shards will cache and only the failed ones will re-run.
    - The ACAF and Clinvar interval lists appear to be particularly challenging for the GVS extract workflows and have some special configuration requirements as described below.
      - For both PGEN and VCF extracts of ACAF and Clinvar the following considerations apply:
        - Specify a `split_intervals_disk_size_override` of 1000 (GiB).
        - For `GvsExtractCallsetPgenMerged` only, specify a `scatter_count` of 25000 (overrides the default 30000 which gets inflated by this particular interval list to result in too many shards for Cromwell).
      - For both PGEN and VCF extracts of ACAF only:
        - Specify an `extract_overhead_memory_override_gib` of 10 (GiB, up from the default of 3 GiB).
        - Specify a `y_bed_weight_scaling` of 8 (up from the default of 4).
        - When re-running the extract workflow with call caching enabled, it will be necessary to increase memory in the `ExtractTask` / `PgenExtractTask` tasks. Due to the way call caching works in Cromwell (i.e. the `memory` runtime attribute is not part of the call caching hashes), it is possible to edit the value of the `memory` runtime attribute of a task _in the WDL_ without breaking call caching. However, do *not* alter the value of the `memory_gib` input parameter as changing that absolutely will break call caching and will cause tens of thousands of shards to re-run needlessly! Both VCF and PGEN extracts can have their memory set to `"50 GiB"` for the call-caching re-run. Most extract shards should finish on the first re-run attempt, but a few stragglers will likely OOM and automatically re-run with more memory.
    - These workflows do not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.

### Smaller Sample Lists

#### VCF Outputs
You can take advantage of our existing sub-cohort WDL, `GvsExtractCohortFromSampleNames.wdl`, to create VCFs for a subset of callset samples.
- Specify the same `call_set_identifier`, `gvs_dataset` (same as `dataset_name` in other runs), `gvs_project` (same as `project_id` in other runs), and `filter_set_name` that were used in the creation of the main callset. 
- Specify a unique `cohort_table_prefix` to this subset of samples so as to not overwrite the prepare tables for the full callset.
- You will need to either fill out `cohort_sample_names` with a GCS path to a newline-delimited list of the sample names or `cohort_sample_names_array` with an Array of sample name Strings.  The `cohort_sample_names_array` input will take precedence over `cohort_sample_names` if both are set. 
- Be sure to set the `output_gcs_dir` to the proper path in the AoU delivery bucket so you don't need to copy the output files there yourself once the workflow has finished.
- This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.

#### PGEN Outputs
1. You will need to run the `GvsPrepareRangesCallset` workflow for each subset of samples for which a PGEN or VCF deliverable is required for the callset.
    - This workflow transforms the data in the vet, ref_ranges, and samples tables into a schema optimized for extract.
    - The `enable_extract_table_ttl` input should be set to `true` (the default value is `false`), which will add a TTL of two weeks to the tables it creates.
    - `extract_table_prefix` should be set to a name that is unique to the given sample list. See the [naming conventions doc](https://docs.google.com/document/d/1pNtuv7uDoiOFPbwe4zx5sAGH7MyxwKqXkyrpNmBxeow) for guidance on what to use.
    - Specify the `interval_list` appropriate for the PGEN / VCF extraction run you are performing.
    - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.
1. `GvsExtractCallsetPgenMerged` workflow
    - Specify the same `call_set_identifier`, `dataset_name`, `project_id`, `extract_table_prefix`, and `interval_list` that were used in the `GvsPrepareRangesCallset` run documented above.
    - Specify the `interval_weights_bed` appropriate for the PGEN extraction run you are performing. `gs://gvs_quickstart_storage/weights/gvs_full_vet_weights_1kb_padded_orig.bed` is the interval weights BED used for Quickstart.
    - Select the workflow option "Retry with more memory" and choose a "Memory retry factor" of 1.5
    - Set the `extract_maxretries_override` input to 5, `split_intervals_disk_size_override` to 1000, `scatter_count` to 25000, and `y_bed_weight_scaling` to 8 to start; you will likely have to adjust one or more of these values in subsequent attempts.
    - `GvsExtractCallsetPgen` currently defaults to 100 alt alleles maximum, which means that any sites having more than that number of alt alleles will be dropped.
    - Be sure to set the `output_gcs_dir` to the proper path in the AoU delivery bucket so you don't need to copy the output files there yourself once the workflow has finished.
    - For `GvsExtractCallsetPgen` (which is called by `GvsExtractCallsetPgenMerged`), if one (or several) of the `PgenExtractTask` shards fail because of angry cloud, you can re-run the workflow with the exact same inputs with call caching turned on; the successful shards will cache and only the failed ones will re-run.
    - If you want to collect the monitoring logs from a large number of `Extract` shards, the `summarize_task_monitor_logs.py` script will not work if the task is scattered too wide.  Use the `summarize_task_monitor_logs_from_file.py` script, instead, which takes a FOFN of GCS paths instead of a space-separated series of localized files.
    - This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.


### VID to Participant ID Mapping Table.
Once the VAT has been created, you will need to create a database table mapping the VIDs (Variant IDs) from that table to all the participants in the dataset that share that VID. This table is used by the AoU Researcher Workbench, and will need to be copied over to a location specified by them. 

1. Create the database table. Using (for instance) the BigQuery cloud user interface, run the query below. Note that you should redirect the output of this query to a new database table in the same dataset, for instance by using the 'query settings' feature in the BigQuery cloud user interface. Also note that you will need to specify the `project`, `dataset`, and `vat_table_name` fields before running the query. Further note that this query might take an hour or two to run to completion:
    ```
   CREATE TEMP FUNCTION vidToLocation(vid string)
    RETURNS int64
    AS (
        (CASE SPLIT(vid, '-')[OFFSET(0)]
                            WHEN 'X' THEN 23
                            WHEN 'Y' THEN 24
                            ELSE CAST(SPLIT(vid, '-')[OFFSET(0)] AS int64) END) * 1000000000000 +
                    CAST(SPLIT(vid, '-')[OFFSET(1)] AS int64)
    );
    
    SELECT vat.vid as vid, ARRAY_AGG(SAFE_CAST(si.sample_name as INT64) IGNORE NULLS) AS person_ids
        FROM `<project>.<dataset>.alt_allele` AS aa
                JOIN `<project>.<dataset>.sample_info` AS si
                    ON aa.sample_id = si.sample_id
                JOIN
            (SELECT vid,
                vidToLocation(vid) AS location,
                SPLIT(vid, '-')[OFFSET(2)] AS ref_allele,
                SPLIT(vid, '-')[OFFSET(3)] AS alt_allele
            FROM `<project>.<dataset>.<vat_table_name>`
            GROUP BY vid, location) AS vat
        ON
            vat.ref_allele = aa.ref AND
            vat.alt_allele = aa.allele AND
            vat.location = aa.location
    GROUP BY vat.vid
    ORDER BY
        vidToLocation(vat.vid),
        SPLIT(vat.vid, '-')[OFFSET(2)],
        SPLIT(vat.vid, '-')[OFFSET(3)]
   ```
1. Once the query has successfully finished, you should cluster it on the field `vid`. This can be accomplished using the command below. Note that you will need to specify the `project`, `dataset`, and `mapping_table_name` fields before running the command:
    ```
    bq update --project_id=<project> --clustering_fields=vid <dataset>.<mapping_table_name>
    ```

1. Copy the created mapping table to the dataset specified by the All of Us DRC. I specifically reached out to Justin Cook and Brian Freeman for the dataset to copy to.

1. Note well that there will be a small difference in the number of vids in the VAT and that in the new mapping table that you have just created. For the Echo callset there are 3595 vids in the VAT that don't exist in the new mapping table (use the query below to determine them). The difference occurs because when we generate annotations using Nirvana, Nirvana left align and truncates the indels. For some of these indels, they now are remapped to a slightly different location than we had defined them as in the `alt_allele` table. This is a known issue that will be resolved in a future release.

    ```
   select distinct vid from `<dataset>.<vat_table_name>` where vid not in (select vid from `<dataset>.<mapping_table_name>`) ;
    ```
