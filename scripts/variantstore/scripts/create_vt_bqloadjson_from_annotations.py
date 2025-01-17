from datetime import datetime

import json
import ijson
import gzip
import argparse
import logging
import sys
import Bio.Seq

vat_nirvana_positions_dictionary = {
    "position": "position", # required
    "ref_allele": "refAllele", # required
    "alt_allele": "altAlleles", # required
}

vat_nirvana_variants_dictionary = {
    "vid": "vid", # required
    "contig": "chromosome", # required
    "variant_type": "variantType", # required
    "genomic_location": "hgvsg", # required
    "dbsnp_rsid": "dbsnp",  # nullable
}

vat_nirvana_transcripts_dictionary = {
    "transcript": "transcript", # nullable
    "gene_symbol": "hgnc", # nullable
    "transcript_source": "source", # nullable
    "aa_change": "hgvsp", # nullable
    "consequence": "consequence", # nullable
    "dna_change_in_transcript": "hgvsc", # nullable
    "exon_number": "exons", # nullable
    "intron_number": "introns", # nullable
    "gene_id": "geneId", # nullable
    "is_canonical_transcript": "isCanonical" # nullable
}

vat_nirvana_gvs_alleles_dictionary = {
    "gvs_all_an": "AN", # required
    "gvs_all_ac": "AC", # required
    "gvs_all_af": "AF" # required
}

vat_nirvana_revel_dictionary = {
    "revel": "score" # nullable
}

vat_nirvana_splice_ai_dictionary = {
    "splice_ai_acceptor_gain_score": "acceptorGainScore", # nullable
    "splice_ai_acceptor_gain_distance": "acceptorGainDistance", # nullable
    "splice_ai_acceptor_loss_score": "acceptorLossScore", # nullable
    "splice_ai_acceptor_loss_distance": "acceptorLossDistance", # nullable
    "splice_ai_donor_gain_score": "donorGainScore", # nullable
    "splice_ai_donor_gain_distance": "donorGainDistance", # nullable
    "splice_ai_donor_loss_score": "donorLossScore", # nullable
    "splice_ai_donor_loss_distance": "donorLossDistance" # nullable
}

vat_nirvana_clinvar_dictionary = {
    "clinvar_classification": "significance", # nullable
    "clinvar_last_updated": "lastUpdatedDate", # nullable
    "clinvar_phenotype": "phenotypes" # nullable
}

vat_nirvana_gnomad_dictionary = {
    "gnomad_all_af": "allAf", # nullable
    "gnomad_all_ac": "allAc", # nullable
    "gnomad_all_an": "allAn", # nullable
    "gnomad_failed_filter": "failedFilter" #nullable
}

significance_ordering = [
    "benign",
    "likely benign",
    "uncertain significance",
    "likely pathogenic",
    "pathogenic",
    "drug response",
    "association",
    "risk factor",
    "protective",
    "affects",
    "conflicting data from submitters",
    "other",
    "not provided",
    "'-'"
]

gnomad_ordering = [
    "afr",
    "amr",
    "asj",
    "eas",
    "fin",
    "nfe",
    "oth",
    "sas"
]

gvs_subpopulations = [
    "afr",
    "amr",
    "eas",
    "eur",
    "mid",
    "oth",
    "sas"
]

def check_filtering(variant):
    # skip any row (with a warning) if no gvsAnnotations exist
    if variant.get("gvsAnnotations") == None: # <-- enum since we need this to be in tandem with the custom annotations header / template
        logging.warning("WARNING: There has been an error in creating custom annotations for AC/AF/AN", variant.get("vid"))
        return False
    # skip any row (with a warning) if the AC value is 0
    elif variant["gvsAnnotations"].get("AC") == 0:
        logging.warning("WARNING: Its AC is 0 so we are dropping this variant", variant.get("vid"))
        return False
    # skip any row (with a warning) if AC, AN or AF is missing
    elif variant["gvsAnnotations"].get("AC") == None:
        logging.warning("WARNING: There has been an error -- there is no AC value -- should AN be 0 for this variant?", variant.get("vid"))
        return False
    elif variant["gvsAnnotations"].get("AN") == None:
        logging.warning("WARNING: There has been an error -- there is an AC value -- but no AN value", variant.get("vid"))
        return False
    elif variant["gvsAnnotations"].get("AF") == None:
        logging.warning("WARNING: There has been an error -- there is an AC value -- but no AF value", variant.get("vid"))
        return False
    else:
        return True

def get_gnomad_subpop(gnomad_obj):
    row = {}
    max_ac = None
    max_an = None
    max_af = None
    max_subpop = ""
    ## TODO defend against possible unexpected values appearing in gnomad_subpop (values not in gnomad_ordering)
    for gnomad_subpop in gnomad_ordering: # since we cycle through this in order, if there is a tie, we just ignore it because we wouldn't chose it anyway based on order
        subpop_af_key = "".join([gnomad_subpop, "Af"])
        subpop_ac_key = "".join([gnomad_subpop, "Ac"])
        subpop_an_key = "".join([gnomad_subpop, "An"])
        subpop_af_val = gnomad_obj.get(subpop_af_key) # note that these can be null if there is no value in the annotations. They will be null in the VAT
        subpop_ac_val = gnomad_obj.get(subpop_ac_key)
        subpop_an_val = gnomad_obj.get(subpop_an_key)
        # here we set the subpopulation ac/an/af values
        row["_".join(["gnomad", gnomad_subpop, "an"])] = subpop_an_val
        row["_".join(["gnomad", gnomad_subpop, "ac"])] = subpop_ac_val
        row["_".join(["gnomad", gnomad_subpop, "af"])] = subpop_af_val
        if subpop_af_val != None and (max_af == None or subpop_af_val > max_af): # this will also set the first max_af value
            max_subpop = gnomad_subpop
            max_ac = subpop_ac_val
            max_an = subpop_an_val
            max_af = subpop_af_val
    # here we set the MAX subpopulation ac/an/af values
    row["gnomad_max_subpop"] = max_subpop
    row["gnomad_max_ac"] = max_ac
    row["gnomad_max_an"] = max_an
    row["gnomad_max_af"] = max_af
    return row

def get_subpopulation_calculations(subpop_annotations):
    row = {}
    max_ac = None
    max_an = None
    max_af = None
    max_sc = None
    max_subpop = ""
    for gvs_subpop in gvs_subpopulations:
        subpop_ac_val = subpop_annotations.get("_".join(["AC", gvs_subpop]), 0)
        subpop_an_val = subpop_annotations.get("_".join(["AN", gvs_subpop]), 0)
        subpop_af_val = subpop_annotations.get("_".join(["AF", gvs_subpop]), None)
        # From the VDS, we get the homozygote_count and subtract it from the AC value
        subpop_sc_val = subpop_annotations.get("_".join(["AC", gvs_subpop]), 0) - subpop_annotations.get("_".join(["Hom", gvs_subpop]), 0)

        row["_".join(["gvs", gvs_subpop, "ac"])] = subpop_ac_val
        row["_".join(["gvs", gvs_subpop, "an"])] = subpop_an_val
        row["_".join(["gvs", gvs_subpop, "af"])] = subpop_af_val
        row["_".join(["gvs", gvs_subpop, "sc"])] = subpop_sc_val
        if subpop_af_val != None and (max_af == None or subpop_af_val > max_af): # this will also set the first max_af value
            max_subpop = gvs_subpop
            max_ac = subpop_ac_val
            max_an = subpop_an_val
            max_af = subpop_af_val
            max_sc = subpop_sc_val
    # here we set the MAX subpopulation ac/an/af values
    row["gvs_max_subpop"] = max_subpop
    row["gvs_max_ac"] = max_ac
    row["gvs_max_an"] = max_an
    row["gvs_max_af"] = max_af
    row["gvs_max_sc"] = max_sc
    return row

def make_annotated_json_row(row_position, row_ref, row_alt, variant_line, transcript_line):
    row = {}
    row["position"] = row_position # this is a required field
    row["ref_allele"] = row_ref # this is a required field
    row["alt_allele"] = row_alt # this is a required field

    for vat_variants_fieldname in vat_nirvana_variants_dictionary.keys():  # like "contig"
        nirvana_variants_fieldname = vat_nirvana_variants_dictionary.get(vat_variants_fieldname)
        variant_fieldvalue = variant_line.get(nirvana_variants_fieldname)
        row[vat_variants_fieldname] = variant_fieldvalue

    if transcript_line != None:
        for vat_transcripts_fieldname in vat_nirvana_transcripts_dictionary.keys():  # like "transcript"
            nirvana_transcripts_fieldname = vat_nirvana_transcripts_dictionary.get(vat_transcripts_fieldname)
            transcript_fieldvalue = transcript_line.get(nirvana_transcripts_fieldname)
            row[vat_transcripts_fieldname] = transcript_fieldvalue

        if (variant_line.get("spliceAI") != None) and (transcript_line.get("hgnc") != None):
            splice_ai_list = variant_line["spliceAI"]
            for splice_ai_obj in splice_ai_list:
                # get the splice AI value that matches to the transcript_line transcripts.hgnc to "spliceAI.hgnc"
                if splice_ai_obj.get("hgnc") == transcript_line["hgnc"]:
                    for vat_splice_ai_fieldname in vat_nirvana_splice_ai_dictionary.keys():  # like "splice_ai_acceptor_gain_score"
                        nirvana_splice_ai_fieldname = vat_nirvana_splice_ai_dictionary.get(vat_splice_ai_fieldname)
                        splice_ai_fieldvalue = splice_ai_obj.get(nirvana_splice_ai_fieldname)
                        row[vat_splice_ai_fieldname] = splice_ai_fieldvalue

    if variant_line.get("gnomad") != None:
        for vat_gnomad_fieldname in vat_nirvana_gnomad_dictionary.keys():  # like "gnomad_all_af"
            nirvana_gnomad_fieldname = vat_nirvana_gnomad_dictionary.get(vat_gnomad_fieldname)
            gnomad_fieldvalue = variant_line.get("gnomad").get(nirvana_gnomad_fieldname)
            row[vat_gnomad_fieldname] = gnomad_fieldvalue
        gnomad_row = get_gnomad_subpop(variant_line["gnomad"])
        row.update(gnomad_row)

    if variant_line.get("clinvar") != None:
        clinvar_objs = variant_line["clinvar"]
        var_ref = variant_line["refAllele"]
        var_alt = variant_line["altAllele"]
        significance_values = [] # ordered by Benign, Likely Benign, Uncertain significance, Likely pathogenic, Pathogenic # https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/
        updated_dates = [] # grab the most recent
        phenotypes = [] # ordered alphabetically
        clinvar_ids = [] # For easy validation downstream
        # Note that inside the clinvar array, are multiple objects that may or may not be the one we are looking for.
        # We check by making sure the ref and alt are the same including any reverse complements
        ## TODO add clinvar star rating!!!
        for clinvar_obj in clinvar_objs:
            # get only the clinvar objs with right variant and the id that starts with RCV
            if (((clinvar_obj.get("refAllele") == var_ref) and (clinvar_obj.get("altAllele") == var_alt)) or
            ((clinvar_obj.get("refAllele") == Bio.Seq.reverse_complement(variant_line["refAllele"])) and (clinvar_obj.get("altAllele") == Bio.Seq.reverse_complement(variant_line["altAllele"])))) and (clinvar_obj.get("id")[:3] == "RCV"):
                clinvar_ids.append(clinvar_obj.get("id"))
                significance_values.extend([x.lower() for x in clinvar_obj.get("significance")])
                updated_dates.append(clinvar_obj.get("lastUpdatedDate"))
                phenotypes.extend(clinvar_obj.get("phenotypes"))
                ## TODO add the ("variationId") and the ("reviewStatus")--note that the reviewStatus will need to maintain the ordering of the significance arrays
                ## we need to do this with a tuple so that the reviewStatus lines up with the significance (since significance seems to be an array, while star is a single value)
        if len(clinvar_ids) > 0:
            ordered_significance_values = []
            # We want to collect all the significance values and order them by the significance_ordering list
            # So I will loop through the significance_ordering values and check for matching values in the significance_values list and put them in a new list
            # And then anything that did not match will get added at the end of the list (sorted alpha)
            for value in significance_ordering:
                if value in significance_values:
                    ordered_significance_values.append(value) # this adds the id to the end of the list
            values_not_accounted_for = list(set(significance_values).difference(ordered_significance_values))
            values_not_accounted_for.sort() # alphabetize this so it is deterministic
            ordered_significance_values.extend(values_not_accounted_for) # add any values that aren't in significance_ordering to the end
            row["clinvar_id"] = clinvar_ids # array
            row["clinvar_classification"] = ordered_significance_values # special sorted array
            updated_dates.sort(key=lambda date: datetime.strptime(date, "%Y-%m-%d")) # note: method is in-place, and returns None
            row["clinvar_last_updated"] = updated_dates[-1] # most recent date
            row["clinvar_phenotype"] = sorted(phenotypes) # union of all phenotypes

    if variant_line.get("revel") != None:
        row["revel"] = variant_line.get("revel").get("score")

    gvs_annotations = variant_line["gvsAnnotations"]
    for vat_gvs_alleles_fieldname in vat_nirvana_gvs_alleles_dictionary.keys():  # like "gvs_all_ac"
        nirvana_gvs_alleles_fieldname = vat_nirvana_gvs_alleles_dictionary.get(vat_gvs_alleles_fieldname)
        gvs_alleles_fieldvalue = gvs_annotations.get(nirvana_gvs_alleles_fieldname)
        row[vat_gvs_alleles_fieldname] = gvs_alleles_fieldvalue
        row["gvs_all_sc"] = gvs_annotations.get("AC") - gvs_annotations.get('Hom') ## TODO: should this include the hemi value?


    subpopulation_info = get_subpopulation_calculations(gvs_annotations)
    row.update(subpopulation_info)

    return row


def make_positions_json(annotated_json, output_json):
    output_file=gzip.open(output_json, 'w')

    if annotated_json.endswith(".gz"):
        json_data = gzip.open(annotated_json, 'rb')
    else:
        json_data = open(annotated_json, 'rb')

    positions = ijson.items(json_data, 'item', use_float=True)

    position_count = 0
    for p in positions:
        position_count += 1
        position=p['position']
        refAllele=p['refAllele'] # ref_allele
        altAlleles=p['altAlleles'] # this is an array that we need to split into each alt_allele
        variants=p['variants']
        if '*' in altAlleles:
            altAlleles.remove('*')
        ## there should be the same number of altAlleles as variants, right?
        assert len(altAlleles)==len(variants)
        # row for each variant - transcript
        # so let's start with each variant
        for index, variant in enumerate(variants):
            # check if it's a row we care about
            if check_filtering(variant) is False:
                continue
            # remember that we want one for each variant-transcript and variant-null for variants without transcripts
            if variant.get("transcripts") == None:
                # then we make a special row
                row = make_annotated_json_row(position, refAllele, altAlleles[index], variant, None)
                json_str = json.dumps(row) + "\n"
                json_bytes = json_str.encode('utf-8')
                output_file.write(json_bytes)
            else:
                transcript_lines = variant.get("transcripts")
                # Collect all the transcript sources and check for if they contain Ensembl
                sources = [transcript.get('source') for transcript in transcript_lines]
                if "Ensembl" in sources:
                    for transcript in transcript_lines:
                        if transcript.get('source') == "Ensembl":
                            row = make_annotated_json_row(position, refAllele, altAlleles[index], variant, transcript)
                            json_str = json.dumps(row) + "\n"
                            json_bytes = json_str.encode('utf-8')
                            output_file.write(json_bytes)
                else:
                    # if there are transcripts, but they are not Ensembl, we now only want one row in the VAT, not one row per transcript
                    row = make_annotated_json_row(position, refAllele, altAlleles[index], variant, None)
                    json_str = json.dumps(row) + "\n"
                    json_bytes = json_str.encode('utf-8')
                    output_file.write(json_bytes)
    output_file.close()
    json_data.close()

    if position_count == 0:
        logging.warning(f"WARNING: Found no items in annotated json file: {annotated_json}")
        sys.exit(0)

def make_annotation_json(annotated_json, output_json, logging):
    logging.info("Making the positions json")
    make_positions_json(annotated_json, output_json)
    logging.info("Done")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Create BQ load friendly json for VAT VT table creation')
    parser.add_argument('--annotated_json', type=str, help='nirvana created annotation json', required=True)
    parser.add_argument('--output_vt_json', type=str, help='name of the vt json', required=True)

    args = parser.parse_args()

    make_annotation_json(args.annotated_json,
                         args.output_vt_json,
                         logging)
