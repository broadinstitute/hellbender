import argparse
import hail as hl


def vds_mt(vds_path):
    vds = hl.vds.read_vds(vds_path)
    # This is currently required as otherwise GT takes priority but AC will be wrong as alleles will be in different
    # orders in the VCF versus VDS representations.
    vds.variant_data = vds.variant_data.drop('GT')
    mt = hl.vds.to_dense_mt(vds)
    dense_vds_path = 'vds_dense.mt'
    mt.write(dense_vds_path, overwrite=True)
    mt = hl.read_matrix_table(dense_vds_path).key_rows_by('locus')
    return mt


def vcf_mt(vcf_paths):
    # Import a VCF that we will use a truth data for this test
    # Setting array_elements_required to false is done as a workaround because Hail has a hard time with unconventional fields with empty values e.g. AS_YNG=.,.,. Avoiding explicitly acknowledging the use of missing elements in arrays requires Hail to make a decision in several ambiguous cases
    mt = hl.import_vcf(vcf_paths, force_bgz=True, reference_genome='GRCh38', array_elements_required=False).key_rows_by('locus')
    return mt


def joined_mt(mt_path):
    joined = hl.experimental.full_outer_join_mt(vcf_mt, vds_mt)
    joined = joined.rename({
        'left_col' : 'vcf_col',
        'right_col': 'vds_col',
        'left_row' : 'vcf_row',
        'right_row': 'vds_row',
        'left_entry' : 'vcf_entry',
        'right_entry': 'vds_entry',
    })

    joined.write(mt_path, overwrite=True)
    mt = hl.read_matrix_table(mt_path)
    return mt


if __name__ == '__main__':oh
    parser = argparse.ArgumentParser()
    parser.add_argument('--vds-path', required=True,
                        help='Input VDS for tieout')
    parser.add_argument('--joined-matrix-table-path', required=True, help='Output joined MatrixTable')
    parser.add_argument('vcf_paths', nargs='+')

    args = parser.parse_args()

    vds_mt = vds_mt(args.vds_path)
    vcf_mt = vcf_mt(args.vcf_paths)
    joined_mt(args.joined_matrix_table_path)


