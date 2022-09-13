version 1.0

workflow GvsCallsetStatistics {
    input {
        String project_id
        String dataset_name
        String filter_set_name
        String extract_prefix
        String metrics_table = "~{extract_prefix}_sample_metrics"
        String aggregate_metrics_table = "~{extract_prefix}_sample_metrics_agg"
        String statistics_table = "~{extract_prefix}_stats"
    }

    call CreateTables {
        input:
            project_id = project_id,
            dataset_name = dataset_name,
            metrics_table = metrics_table,
            aggregate_metrics_table = aggregate_metrics_table,
            statistics_table = statistics_table
    }

    # Only collect statistics for the autosomal chromosomes, the first 22 in our location scheme.
    scatter(chrom in range(22)) {
        call CollectStatisticsForChromosome {
            input:
                go = CreateTables.done,
                project_id = project_id,
                dataset_name = dataset_name,
                filter_set_name = filter_set_name,
                extract_prefix = extract_prefix,
                metrics_table = metrics_table,
                chromosome = chrom + 1 # 0-based ==> 1-based
        }
    }

    call AggregateStatisticsAcrossChromosomes {
        input:
            go = CollectStatisticsForChromosome.done[0],
            project_id = project_id,
            dataset_name = dataset_name,
            filter_set_name = filter_set_name,
            extract_prefix = extract_prefix,
            metrics_table = metrics_table,
            aggregate_metrics_table = aggregate_metrics_table
    }

    call CollectStatistics {
        input:
            go = AggregateStatisticsAcrossChromosomes.done,
            project_id = project_id,
            dataset_name = dataset_name,
            filter_set_name = filter_set_name,
            extract_prefix = extract_prefix,
            metrics_table = metrics_table,
            aggregate_metrics_table = aggregate_metrics_table,
            statistics_table = statistics_table
    }
}

task CreateTables {
    input {
        String project_id
        String dataset_name
        String metrics_table
        String aggregate_metrics_table
        String statistics_table
    }
    command <<<
        set -o errexit -o nounset -o xtrace -o pipefail

        set +o errexit
        bq --project_id=~{project_id} show ~{dataset_name}.~{metrics_table}
        BQ_SHOW_RC=$?
        set -o errexit

        if [ $BQ_SHOW_RC -eq 0 ]; then
            echo "Output metrics table already exists, exiting."
            exit 1
        fi

        set +o errexit
        bq --project_id=~{project_id} show ~{dataset_name}.~{aggregate_metrics_table}
        BQ_SHOW_RC=$?
        set -o errexit

        if [ $BQ_SHOW_RC -eq 0 ]; then
            echo "Output aggregate metrics table already exists, exiting."
            exit 1
        fi

        set +o errexit
        bq --project_id=~{project_id} show ~{dataset_name}.~{statistics_table}
        BQ_SHOW_RC=$?
        set -o errexit

        if [ $BQ_SHOW_RC -eq 0 ]; then
            echo "Output statistics table already exists, exiting."
            exit 1
        fi

        # Schemas extracted programatically: https://stackoverflow.com/a/66987934
        #
        # After cleaning up header and quotes:
        #
        # cat raw.json | jq -M '[.[]|{name,type,mode}' > clean.json

        cat > metrics_schema.json <<FIN
        [
          {
            "name": "filter_set_name",
            "type": "STRING",
            "mode": "NULLABLE"
          },
          {
            "name": "sample_id",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "variant_entries",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "del_count",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "ins_count",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "snp_count",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "ti_count",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "tv_count",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "snp_het_count",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "snp_homvar_count",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "indel_het_count",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "indel_homvar_count",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "singleton",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "pass_qc",
            "type": "INT64",
            "mode": "NULLABLE"
          }
        ]
        FIN

        cat > statistics_schema.json <<FIN
        [
          {
            "name": "sample_id",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "sample_name",
            "type": "STRING",
            "mode": "NULLABLE"
          },
          {
            "name": "del_count",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "m_del_count",
            "type": "FLOAT64",
            "mode": "NULLABLE"
          },
          {
            "name": "mad_del_count",
            "type": "FLOAT64",
            "mode": "NULLABLE"
          },
          {
            "name": "pass_del_count",
            "type": "BOOL",
            "mode": "NULLABLE"
          },
          {
            "name": "ins_count",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "m_ins_count",
            "type": "FLOAT64",
            "mode": "NULLABLE"
          },
          {
            "name": "mad_ins_count",
            "type": "FLOAT64",
            "mode": "NULLABLE"
          },
          {
            "name": "pass_ins_count",
            "type": "BOOL",
            "mode": "NULLABLE"
          },
          {
            "name": "snp_count",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "m_snp_count",
            "type": "FLOAT64",
            "mode": "NULLABLE"
          },
          {
            "name": "mad_snp_count",
            "type": "FLOAT64",
            "mode": "NULLABLE"
          },
          {
            "name": "pass_snp_count",
            "type": "BOOL",
            "mode": "NULLABLE"
          },
          {
            "name": "singleton",
            "type": "INT64",
            "mode": "NULLABLE"
          },
          {
            "name": "m_singleton",
            "type": "FLOAT64",
            "mode": "NULLABLE"
          },
          {
            "name": "mad_singleton",
            "type": "FLOAT64",
            "mode": "NULLABLE"
          },
          {
            "name": "pass_singleton",
            "type": "BOOL",
            "mode": "NULLABLE"
          },
          {
            "name": "ins_del_ratio",
            "type": "FLOAT64",
            "mode": "NULLABLE"
          },
          {
            "name": "m_ins_del_ratio",
            "type": "FLOAT64",
            "mode": "NULLABLE"
          },
          {
            "name": "mad_ins_del_ratio",
            "type": "FLOAT64",
            "mode": "NULLABLE"
          },
          {
            "name": "pass_ins_del_ratio",
            "type": "BOOL",
            "mode": "NULLABLE"
          },
          {
            "name": "ti_tv_ratio",
            "type": "FLOAT64",
            "mode": "NULLABLE"
          },
          {
            "name": "m_ti_tv_ratio",
            "type": "FLOAT64",
            "mode": "NULLABLE"
          },
          {
            "name": "mad_ti_tv_ratio",
            "type": "FLOAT64",
            "mode": "NULLABLE"
          },
          {
            "name": "pass_ti_tv_ratio",
            "type": "BOOL",
            "mode": "NULLABLE"
          },
          {
            "name": "snp_het_homvar_ratio",
            "type": "FLOAT64",
            "mode": "NULLABLE"
          },
          {
            "name": "m_snp_het_homvar_ratio",
            "type": "FLOAT64",
            "mode": "NULLABLE"
          },
          {
            "name": "mad_snp_het_homvar_ratio",
            "type": "FLOAT64",
            "mode": "NULLABLE"
          },
          {
            "name": "pass_snp_het_homvar_ratio",
            "type": "BOOL",
            "mode": "NULLABLE"
          },
          {
            "name": "indel_het_homvar_ratio",
            "type": "FLOAT64",
            "mode": "NULLABLE"
          },
          {
            "name": "m_indel_het_homvar_ratio",
            "type": "FLOAT64",
            "mode": "NULLABLE"
          },
          {
            "name": "mad_indel_het_homvar_ratio",
            "type": "FLOAT64",
            "mode": "NULLABLE"
          },
          {
            "name": "pass_indel_het_homvar_ratio",
            "type": "BOOL",
            "mode": "NULLABLE"
          }
        ]
        FIN

        bq mk --table ~{project_id}:~{dataset_name}.~{metrics_table} metrics_schema.json
        bq mk --table ~{project_id}:~{dataset_name}.~{aggregate_metrics_table} metrics_schema.json
        bq mk --table ~{project_id}:~{dataset_name}.~{statistics_table} statistics_schema.json
    >>>
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_2022_08_22"
    }
    output {
        Boolean done = true
    }
}

task CollectStatisticsForChromosome {
    input {
        Boolean go = true
        String project_id
        String dataset_name
        String filter_set_name
        String extract_prefix
        String metrics_table
        Int chromosome
    }
    command <<<
        set -o errexit -o nounset -o xtrace -o pipefail

        echo "project_id = ~{project_id}" > ~/.bigqueryrc

        bq query --location=US --project_id=~{project_id} --format=csv --use_legacy_sql=false '
        CREATE TEMPORARY FUNCTION titv(ref STRING, allele STRING)
        RETURNS STRING
            LANGUAGE js AS """
                if ( ref.length > 1 || allele.length > 1) {
                    return "other";
                } else if ( (ref == "A" && allele == "G") ||
                            (ref == "G" && allele == "A") ||
                            (ref == "C" && allele == "T") ||
                            (ref == "T" && allele == "C") ) {
                    return "ti";
                } else {
                    return "tv";
                }
        """;

        CREATE TEMPORARY FUNCTION type(ref STRING, allele STRING, gt_str STRING)
        RETURNS STRING
            LANGUAGE js AS """

        alts = allele.split(",")

        // get the the non-reference allele indexes
        ai = gt_str.replace("|","/").split("/").filter(i => i != "0");

        // the the distinct set of lengths of the alternates
        alt_lengths = new Set(ai.map(i => alts[parseInt(i)-1].length))

        if (alt_lengths.size > 1) {
            return "complex"
        } else {
            // get first (only) element
            al = alt_lengths.keys().next().value

            if ( ref.length == al && al == 1) {
                return "snp"
            } else if (ref.length > al) {
                return "del"
            } else if (ref.length < al) {
                return "ins"
            } else {
                return "other"
            }
        }
        """;

        INSERT `~{project_id}.~{dataset_name}.~{metrics_table}` (
            filter_set_name,
            sample_id,
            variant_entries,
            del_count,
            ins_count,
            snp_count,
            ti_count,
            tv_count,
            snp_het_count,
            snp_homvar_count,
            indel_het_count,
            indel_homvar_count,
            singleton,
            pass_qc
        )
        SELECT "~{filter_set_name}" filter_set_name,
               sample_id,
               count(1) variant_entries,
               SUM(CASE WHEN type = "del" THEN 1 ELSE 0 END) del_count,
               SUM(CASE WHEN type = "ins" THEN 1 ELSE 0 END) ins_count,
               SUM(CASE WHEN type = "snp" THEN 1 ELSE 0 END) snp_count,
               SUM(CASE WHEN type = "snp" AND titv = "ti" THEN 1 ELSE 0 END) ti_count, # TODO: minimize alleles
               SUM(CASE WHEN type = "snp" AND titv = "tv" THEN 1 ELSE 0 END) tv_count, # TODO: minimize alleles
               SUM(CASE WHEN type = "snp" AND gt_type = "het" THEN 1 ELSE 0 END) snp_het_count,
               SUM(CASE WHEN type = "snp" AND gt_type = "homvar" THEN 1 ELSE 0 END) snp_homvar_count,
               SUM(CASE WHEN type IN ("ins","del") AND gt_type = "het" THEN 1 ELSE 0 END) indel_het_count,
               SUM(CASE WHEN type IN ("ins","del") AND gt_type = "homvar" THEN 1 ELSE 0 END) indel_homvar_count,
               COUNTIF(not in_gnomad) singleton,
               null AS pass_qc
        FROM (
            SELECT sample_id,
                   type(ref, alt, call_GT) as type,
                   CASE WHEN INSTR(call_GT, "0") > 0 THEN "het" ELSE "homvar" END as gt_type,
                   titv(ref, alt) as titv,
                   CASE WHEN gnomad.location IS NULL THEN false ELSE true END in_gnomad
            FROM `~{project_id}.~{dataset_name}.~{extract_prefix}__VET_DATA` v
            LEFT JOIN `spec-ops-aou.gvs_public_reference_data.gnomad_v3_sites` gnomad ON (v.location = gnomad.location)
            WHERE call_GT != "./."
            AND v.location >= ~{chromosome}000000000000
            AND v.location < ~{chromosome + 1}000000000000) GROUP BY 1,2

        '

    >>>
    output {
        Boolean done = true
    }
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_2022_08_22"
    }
}

task AggregateStatisticsAcrossChromosomes {
    input {
        Boolean go
        String project_id
        String dataset_name
        String filter_set_name
        String extract_prefix
        String metrics_table
        String aggregate_metrics_table
    }
    command <<<
        bq query --location=US --project_id=~{project_id} --format=csv --use_legacy_sql=false '
        INSERT `~{project_id}.~{dataset_name}.~{aggregate_metrics_table}` (
            filter_set_name,
            sample_id,
            variant_entries,
            del_count,
            ins_count,
            snp_count,
            ti_count,
            tv_count,
            snp_het_count,
            snp_homvar_count,
            indel_het_count,
            indel_homvar_count,
            singleton,
            pass_qc
        )
        SELECT "~{filter_set_name}" filter_set_name,
            sample_id,
            SUM(variant_entries) variant_entries,
            SUM(del_count) del_count,
            SUM(ins_count) ins_count,
            SUM(snp_count) snp_count,
            SUM(ti_count) ti_count,
            SUM(tv_count) tv_count,
            SUM(snp_het_count) snp_het_count,
            SUM(snp_homvar_count) snp_homvar_count,
            SUM(indel_het_count) indel_het_count,
            SUM(indel_homvar_count) indel_homvar_count,
            SUM(singleton) singleton,
            null AS pass_qc
        FROM `~{project_id}.~{dataset_name}.~{metrics_table}` GROUP BY 1,2

        '
    >>>
    output {
        Boolean done = true
    }
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_2022_08_22"
    }
}

task CollectStatistics {
    input {
        Boolean go
        String project_id
        String dataset_name
        String filter_set_name
        String extract_prefix
        String metrics_table
        String aggregate_metrics_table
        String statistics_table
    }
    command <<<

        bq query --location=US --project_id=~{project_id} --format=csv --use_legacy_sql=false '

        WITH fss AS (
          SELECT *,
                 (ins_count / del_count) as ins_del_ratio,
                 (ti_count / tv_count) as ti_tv_ratio,
                 (snp_het_count / snp_homvar_count) snp_het_homvar_ratio,
                 (indel_het_count / indel_homvar_count) as indel_het_homvar_ratio
          FROM `~{project_id}.~{dataset_name}.~{aggregate_metrics_table}`
          WHERE filter_set_name = "~{filter_set_name}"),
        medians AS (
            SELECT
                `bqutil`.fn.median(ARRAY_AGG(del_count IGNORE NULLS)) as m_del_count,
                `bqutil`.fn.median(ARRAY_AGG(ins_count IGNORE NULLS)) as m_ins_count,
                `bqutil`.fn.median(ARRAY_AGG(snp_count IGNORE NULLS)) as m_snp_count,
                `bqutil`.fn.median(ARRAY_AGG(singleton IGNORE NULLS)) as m_singleton,
                `bqutil`.fn.median(ARRAY_AGG(ins_del_ratio IGNORE NULLS)) as m_ins_del_ratio,
                `bqutil`.fn.median(ARRAY_AGG(ti_tv_ratio IGNORE NULLS)) as m_ti_tv_ratio,
                `bqutil`.fn.median(ARRAY_AGG(snp_het_homvar_ratio IGNORE NULLS)) as m_snp_het_homvar_ratio,
                `bqutil`.fn.median(ARRAY_AGG(indel_het_homvar_ratio IGNORE NULLS)) as m_indel_het_homvar_ratio
            FROM fss),
        mads AS (
            SELECT
                `bqutil`.fn.median(ARRAY_AGG(ABS(del_count - m_del_count) IGNORE NULLS)) as mad_del_count,
                `bqutil`.fn.median(ARRAY_AGG(ABS(ins_count - m_ins_count) IGNORE NULLS)) as mad_ins_count,
                `bqutil`.fn.median(ARRAY_AGG(ABS(snp_count - m_snp_count) IGNORE NULLS)) as mad_snp_count,
                `bqutil`.fn.median(ARRAY_AGG(ABS(singleton - m_singleton) IGNORE NULLS)) as mad_singleton,
                `bqutil`.fn.median(ARRAY_AGG(ABS(ins_del_ratio - m_ins_del_ratio) IGNORE NULLS)) as mad_ins_del_ratio,
                `bqutil`.fn.median(ARRAY_AGG(ABS(ti_tv_ratio - m_ti_tv_ratio) IGNORE NULLS)) as mad_ti_tv_ratio,
                `bqutil`.fn.median(ARRAY_AGG(ABS(snp_het_homvar_ratio - m_snp_het_homvar_ratio) IGNORE NULLS)) as mad_snp_het_homvar_ratio,
                `bqutil`.fn.median(ARRAY_AGG(ABS(indel_het_homvar_ratio - m_indel_het_homvar_ratio) IGNORE NULLS)) as mad_indel_het_homvar_ratio
            FROM fss
            CROSS JOIN medians
            WHERE filter_set_name = "~{filter_set_name}")

        SELECT
            fss.sample_id,
            si.sample_name,
            del_count, m_del_count, mad_del_count,
            CASE WHEN del_count BETWEEN m_del_count - 4*mad_del_count AND m_del_count + 4*mad_del_count THEN true ELSE false END pass_del_count,

            ins_count, m_ins_count, mad_ins_count,
            CASE WHEN ins_count BETWEEN m_ins_count - 4*mad_ins_count AND m_ins_count + 4*mad_ins_count THEN true ELSE false END pass_ins_count,

            snp_count, m_snp_count, mad_snp_count,
            CASE WHEN snp_count BETWEEN m_snp_count - 4*mad_snp_count AND m_snp_count + 4*mad_snp_count THEN true ELSE false END pass_snp_count,

            singleton, m_singleton, mad_singleton,
            CASE WHEN singleton BETWEEN m_singleton - 8*mad_singleton AND m_singleton + 8*mad_singleton THEN true ELSE false END pass_singleton,

            ins_del_ratio, m_ins_del_ratio, mad_ins_del_ratio,
            CASE WHEN ins_del_ratio BETWEEN m_ins_del_ratio - 4*mad_ins_del_ratio AND m_ins_del_ratio + 4*mad_ins_del_ratio THEN true ELSE false END pass_ins_del_ratio,

            ti_tv_ratio, m_ti_tv_ratio, mad_ti_tv_ratio,
            CASE WHEN ti_tv_ratio BETWEEN m_ti_tv_ratio - 4*mad_ti_tv_ratio AND m_ti_tv_ratio + 4*mad_ti_tv_ratio THEN true ELSE false END pass_ti_tv_ratio,

            snp_het_homvar_ratio, m_snp_het_homvar_ratio, mad_snp_het_homvar_ratio,
            CASE WHEN snp_het_homvar_ratio BETWEEN m_snp_het_homvar_ratio - 4*mad_snp_het_homvar_ratio AND m_snp_het_homvar_ratio + 4*mad_snp_het_homvar_ratio THEN true ELSE false END pass_snp_het_homvar_ratio,

            indel_het_homvar_ratio, m_indel_het_homvar_ratio, mad_indel_het_homvar_ratio,
            CASE WHEN indel_het_homvar_ratio BETWEEN m_indel_het_homvar_ratio - 4*mad_indel_het_homvar_ratio AND m_indel_het_homvar_ratio + 4*mad_indel_het_homvar_ratio THEN true ELSE false END pass_indel_het_homvar_ratio,
        FROM fss
        JOIN `~{project_id}.~{dataset_name}.~{extract_prefix}__SAMPLES` si ON (fss.sample_id = si.sample_id)
        CROSS JOIN medians
        CROSS JOIN mads
        order by 1

        '
    >>>
    output {
        Boolean done = true
    }
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_2022_08_22"
    }
}
