#!/usr/bin/env bash

# Downloads the HGNC data source from the HGNC website.

curl 'https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&amp;col=gd_app_sym&amp;col=gd_app_name&amp;col=gd_status&amp;col=gd_locus_type&amp;col=gd_locus_group&amp;col=gd_prev_sym&amp;col=gd_prev_name&amp;col=gd_aliases&amp;col=gd_name_aliases&amp;col=gd_pub_chrom_map&amp;col=gd_date_mod&amp;col=gd_date_sym_change&amp;col=gd_date_name_change&amp;col=gd_pub_acc_ids&amp;col=gd_enz_ids&amp;col=gd_pub_eg_id&amp;col=gd_pub_ensembl_id&amp;col=gd_pubmed_ids&amp;col=gd_pub_refseq_ids&amp;col=family.id&amp;col=family.name&amp;col=gd_ccds_ids&amp;col=gd_vega_ids&amp;col=md_eg_id&amp;col=md_mim_id&amp;col=md_refseq_id&amp;col=md_prot_id&amp;col=md_ensembl_id&amp;col=md_ucsc_id&amp;status=Approved&amp;status_opt=2&amp;where=&amp;order_by=gd_app_sym_sort&amp;format=text&amp;limit=&amp;hgnc_dbtag=on&amp;submit=submi' > hgnc_download_$(date +%b%d%Y).tsv

