package org.broadinstitute.hellbender.tools.gvs.ingest;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.gvs.common.CommonCode;
import org.broadinstitute.hellbender.tools.walkers.annotator.AnnotationUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import static htsjdk.variant.vcf.VCFConstants.DEPTH_KEY;
import static org.broadinstitute.hellbender.utils.variant.GATKVCFConstants.AS_RAW_QUAL_APPROX_KEY;
import static org.broadinstitute.hellbender.utils.variant.GATKVCFConstants.AS_VARIANT_DEPTH_KEY;
import static org.broadinstitute.hellbender.utils.variant.GATKVCFConstants.RAW_QUAL_APPROX_KEY;
import static org.broadinstitute.hellbender.utils.variant.GATKVCFConstants.VARIANT_DEPTH_KEY;

/**
 * Expected headers for the Variant Table (VET)
 *     sample_id, // req
 *     location, // req
 *     reference_bases, // req
 *     alternate_bases_alt, // req
 *     alternate_bases_AS_RAW_MQ,
 *     alternate_bases_AS_RAW_MQRankSum,
 *     alternate_bases_AS_QUALapprox,
 *     alternate_bases_AS_RAW_ReadPosRankSum,
 *     alternate_bases_AS_SB_TABLE,
 *     alternate_bases_AS_VarDP,
 *     call_genotype, // req
 *     call_AD,
 *     call_DP, // Laura says consider removing for now-- so similar to AS_VarDP
 *     call_GQ, // req
 *     call_PGT,
 *     call_PID,
 *     call_PS,
 *     call_PL
 *
 */
public enum VetFieldEnum {
    // This where the validation step (required vs not) lives  -- fail if there is missing data for a required field
    // and just leave it empty if not required

    sample_id, // Required-- sample Id for sample
    location, // Required-- encoded chromosome and position

    ref { // Required
        public String getColumnValue(final VariantContext variant, final boolean forceLoadingFromNonAlleleSpecific) {
            final String referenceBase = variant.getReference().getBaseString();
            if (referenceBase == null) {
                throw new IllegalArgumentException("Cannot be missing required value for reference_bases"); // TODO, should this be UserException too?
            }
            return referenceBase;
        }
    },

    alt { // remove "<NON_REF>"
        //TODO what if this field is null and if <NON_REF> is not there--throw an error
        public String getColumnValue(final VariantContext variant, final boolean forceLoadingFromNonAlleleSpecific) {
            List<String> outList = new ArrayList<>();
            for(Allele a : variant.getAlternateAlleles()) {
                if (!a.isNonRefAllele()) { // TODO unit test this
                    outList.add(a.getDisplayString());
                }
            }
            return String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, outList);
        }
    },

    AS_RAW_MQ {
        // Required for VQSR Data
        // can strip off the first one?
        // TODO sci notation?
        public String getColumnValue(final VariantContext variant, final boolean forceLoadingFromNonAlleleSpecific) {
            String out = getAttribute(variant, GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY, null);
            String outNotAlleleSpecific = getAttribute(variant, GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY, null);
            String outNotAlleleSpecificAndOld = getAttribute(variant, GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED, null);
            if (out == null && outNotAlleleSpecific == null && outNotAlleleSpecificAndOld == null) {
                throw new UserException("Cannot be missing required value for alternate_bases.AS_RAW_MQ, RAW_MQandDP or RAW_MQ.");
            }
            if (!forceLoadingFromNonAlleleSpecific && out != null) {
// KC: we are seeing a TON of these!
//                if (!out.endsWith("|0.00")) {
//                    logger.warn("Expected AS_RAW_MQ value to end in |0.00. value is: " + out + " for variant " + variant.toString());
//                }
                out = out.substring(0, out.lastIndexOf("|"));
                String[] outValues = out.split("\\|");
                out = Arrays
                        .stream(outValues)
                        .map(val -> val.endsWith(".00") ? val.substring(0, val.length() - 3) : val)
                        .collect(Collectors.joining(AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM));
                return out;
            // If we have gvcfs that are not allele specific from GATK4 we'll get RAW_MQandDP.
            // We can drop DP here and use AS_VarDP when finalizing RMS Mapping Quality
            } else {
                String outValue;
                if (outNotAlleleSpecific != null) {
                    String[] outValues = outNotAlleleSpecific.split(",");
                    if (outValues.length !=2) {
                        throw new UserException("Expected RAW_MQandDP to be two values separated by a comma.");
                    }
                    // First value is MQ the second is DP. Use the only MQ value we have for all alleles since we're faking allele specific annotations.
                    outValue = outValues[0];
                } else {
                    outValue = outNotAlleleSpecificAndOld;
                }

                double mq = Double.parseDouble(outValue);
                if (variant.getAlleles().size() == 3) {
                    outNotAlleleSpecific = 0 + AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM + (int) mq;
                } else if (variant.getAlleles().size() == 4) {
                    outNotAlleleSpecific = 0 + AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM + (int) mq + AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM + (int) mq;
                } else {
                    throw new UserException("Expected diploid sample to either have 3 alleles (ref, alt, non-ref) or 4 alleles (ref, alt 1, alt 2, non-ref)");
                }
                return outNotAlleleSpecific;
            }
        }
        public boolean isVqsrSpecificField() {
            return true;
        }
    },

    AS_RAW_MQRankSum { // TODO -- maybe rely on 1/1 for call_GT, also get rid of the | at the beginning
        // Required for VQSR Data
        public String getColumnValue(final VariantContext variant, final boolean forceLoadingFromNonAlleleSpecific) {
            // in the case where neither allele is reference, don't return a value
            if (isGenotypeAllNonRef(variant.getGenotype(0))) {
                return "";
            }

            // e.g. AS_RAW_MQRankSum=|1.4,1|NaN;
            String out =  getAttribute(variant, GATKVCFConstants.AS_RAW_MAP_QUAL_RANK_SUM_KEY, null);

            if (forceLoadingFromNonAlleleSpecific || out == null) {
                // Try to use non-AS version
                // TODO: it looks like the AS_RAW version also trims to a single decimal point??
                // e.g. MQRankSum=1.465 and turn it into |1.465,1|
                String outNotAlleleSpecific = getAttribute(variant, GATKVCFConstants.MAP_QUAL_RANK_SUM_KEY, null);


                if ( outNotAlleleSpecific == null || "".equals(outNotAlleleSpecific) || outNotAlleleSpecific.contentEquals("||") || outNotAlleleSpecific.contentEquals("|||") ) {
                    return "";
                }

                int numVariants = variant.getAlleles().size() - 2;

                // As of VS-910, we aren't going to consider any number of alleles to necessarily be a UserException
                // Moving this TODO note up here as an open question for how we want to handle this quantity
                // with more than one allele
                // TODO: just replicate rather than distribute, is this right?
                out = "|";
                for (int i = 0; i < numVariants; ++i) {
                    out += outNotAlleleSpecific + ",1|";
                }
            }

            if ( out == null || out.contentEquals("||") || out.contentEquals("|||") ) {
                out = "";
                return out;
            }
            if (out.startsWith("|")) {
                out = out.substring(1);
            } else {
                throw new UserException("Expected AS_RAW_MQRankSum value to begin with a |");
            }
            if (out.endsWith("|NaN")) {
                out = out.substring(0, out.length() - 4);
            } else {
                // for now remove the last value even if not NaN
                out = out.substring(0, out.lastIndexOf("|"));
                //throw new UserException("Expected AS_RAW_MQRankSum value to be ||, ||| or to end in |NaN");
            }
            return out;
        }
        public boolean isVqsrSpecificField() {
            return true;
        }
    },

    QUALapprox {
        // Required for VQSR Data
        public String getColumnValue(final VariantContext variant, final boolean forceLoadingFromNonAlleleSpecific) {
            String out = getAttribute(variant, RAW_QUAL_APPROX_KEY, null);
            if (out == null) {
                throw new UserException("Cannot be missing required value for QUALapprox at site: " + variant);
            }
            return out;
        }
        public boolean isVqsrSpecificField() {
            return true;
        }
    },

    AS_QUALapprox {
        // Required for VQSR Data
        public String getColumnValue(final VariantContext variant, final boolean forceLoadingFromNonAlleleSpecific) {
            String out = getAttribute(variant, AS_RAW_QUAL_APPROX_KEY, null);
            if (out == null) {
                String outNotAlleleSpecific = getAttribute(variant, RAW_QUAL_APPROX_KEY, null);
                if (outNotAlleleSpecific == null) {
                    throw new UserException("Cannot be missing required value for AS_QUALapprox or QUALapprox at site: " + variant);
                }
                return outNotAlleleSpecific;
            }
            if (out.startsWith("|")) {
                out = out.substring(1);
            } else {
                throw new UserException("Expected AS_RAW_MQRankSum value to begin with a |");
            }
            // check to see if there are two or three values, make sure last is smallest, throw out last
            List<String> outList = Arrays.asList(out.split("\\|"));
            if (outList.size() == 2 | outList.size() == 3) { // check length of array -- needs to be 2 or 3
                if (outList.lastIndexOf(Collections.min(outList)) == outList.size() - 1) { // this should be the smallest value
                    out = StringUtils.join(outList.subList(0, outList.size() - 1), AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM);
                } else {
                    throw new UserException(String.format("Expected the final value of AS_QUALapprox to be the smallest at %d", variant.getStart()));
                }
            } else {
                throw new UserException("Expected AS_QUALapprox to have two or three values");
            }
            return out;
        }
        public boolean isVqsrSpecificField() {
            return true;
        }
    },

    AS_RAW_ReadPosRankSum {  // TODO -- maybe rely on 1/1 for call_GT
        public String getColumnValue(final VariantContext variant, final boolean forceLoadingFromNonAlleleSpecific) {
            // in the case where neither allele is reference, don't return a value
            if (isGenotypeAllNonRef(variant.getGenotype(0))) {
                return "";
            }

            // e.g. AS_RAW_ReadPosRankSum=|-0.3,1|0.6,1
            String out =  getAttribute(variant, GATKVCFConstants.AS_RAW_READ_POS_RANK_SUM_KEY, null);

            if (forceLoadingFromNonAlleleSpecific || out == null) {
                // Try to use non-AS version
                // TODO: it looks like the AS_RAW version also trims to a single decimal point??
                // e.g. ReadPosRankSum=-0.511 and turn it into |-0.511,1|
                String outNotAlleleSpecific = getAttribute(variant, GATKVCFConstants.READ_POS_RANK_SUM_KEY, null);

                if ( outNotAlleleSpecific == null || "".equals(outNotAlleleSpecific) || outNotAlleleSpecific.contentEquals("||") || outNotAlleleSpecific.contentEquals("|||") ) {
                    return "";
                }

                int numVariants = variant.getAlleles().size() - 2;
                // Moving this TODO note up here as an open question for how we want to handle this quantity
                // with more than one allele
                // TODO: just replicate rather than distribute, is this right?
                out = "|";
                for (int i = 0; i < numVariants; ++i) {
                    out += outNotAlleleSpecific + ",1|";
                }
                // As of VS-910, we aren't going to consider any number of alleles to necessarily be a UserException
            }

            if (out == null || out.contentEquals("||") || out.contentEquals("|||") ) {
                out = "";
                return out;
            }
            if (out.startsWith("|")) {
                out = out.substring(1);
            } else {
                throw new UserException("Expected AS_RAW_ReadPosRankSum value to begin with a |");
            }
            if (out.endsWith("|NaN")) {
                out = out.substring(0, out.length() - 4);
            } else {
                // for now remove the last value even if not NaN
                out = out.substring(0, out.lastIndexOf("|"));
                //throw new UserException("Expected AS_RAW_ReadPosRankSum value to be ||, ||| or to end in |NaN");
            }
            return out;
        }
        public boolean isVqsrSpecificField() {
            return true;
        }
    },

    AS_SB_TABLE {
        // Required for VQSR Data
        public String getColumnValue(final VariantContext variant, final boolean forceLoadingFromNonAlleleSpecific) {
            String out = getAttribute(variant, GATKVCFConstants.AS_SB_TABLE_KEY, null);
            if (forceLoadingFromNonAlleleSpecific || out == null) {
                String outNotAlleleSpecific = variant.getGenotype(0).getExtendedAttribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY, null).toString();
                String[] outValues = outNotAlleleSpecific.split(",");

                int numVariants = variant.getAlleles().size() - 2;
                out = outValues[0] + "," + outValues[1];
                for (int i = 0; i < numVariants; ++i) {
                    out += "|" + Integer.parseInt(outValues[2])/numVariants + "," + Integer.parseInt(outValues[3])/numVariants;
                }
                // As of VS-910, we aren't going to consider any number of alleles to necessarily be a UserException

                return out;
            }
            if (out.endsWith("|0,0")) {
                out = out.substring(0, out.length() - 4);
            } else {
                // for now remove the last value even if not NaN
                out = out.substring(0, out.lastIndexOf("|"));
                //throw new UserException("Expected AS_SB_TABLE value to end in |0,0");
            }
            return out;
        }
        public boolean isVqsrSpecificField() {
            return true;
        }
    },

    AS_VarDP {
        // Required for VQSR Data
        public String getColumnValue(final VariantContext variant, final boolean forceLoadingFromNonAlleleSpecific) {
            String out = getAttribute(variant, AS_VARIANT_DEPTH_KEY, null);
            if (out == null) {
                String varDP = getAttribute(variant, VARIANT_DEPTH_KEY, null);
                String dP = getAttribute(variant, DEPTH_KEY, null);
                if (varDP == null || dP == null) {
                    throw new UserException("Cannot be missing required value for AS_VarDP, or VarDP and DP, at site:" + variant);
                }
                int refDP = Integer.parseInt(dP) - Integer.parseInt(varDP);
                out = refDP + "|" + varDP + "|";
            }
            if (out.endsWith("|0")) {
                out = out.substring(0, out.length() - 2);
            } else {
                // for now remove the last value even if not NaN
                out = out.substring(0, out.lastIndexOf("|"));
                //throw new UserException("Expected AS_VarDP value to end in |0");
            }
            return out;
        }
        public boolean isVqsrSpecificField() {
            return true;
        }
    },

    call_GT {
        public String getColumnValue(final VariantContext variant, final boolean forceLoadingFromNonAlleleSpecific) {
            if (!variant.hasGenotypes()) {
                throw new UserException("Cannot be missing required value for call.GT");
            }
            // TODO how is missing handled?
            return CommonCode.getGTString(variant);
        }
    },

    call_AD {
        public String getColumnValue(final VariantContext variant, final boolean forceLoadingFromNonAlleleSpecific) {
            String out = variant.getGenotype(0).hasAD() ? Arrays.stream(variant.getGenotype(0).getAD())
                    .mapToObj(String::valueOf)
                    .collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)) : "";
            if (out.equals("")) {
                return out;
            }
            if (out.endsWith(",0")) {
                out = out.substring(0, out.length() - 2);
            } else {
                // for now remove the last value even if not NaN
                out = out.substring(0, out.lastIndexOf("|"));
                //throw new UserException("Expected call_AD to have a final value of 0");
            }
            return out;
        }
    },

    // call_DP { // TODO we can drop whole column since it is similar to AS_VarDP
    // TODO come up with a check-- looks like we can drop this whole one!!!!!
    // public String getColumnValue(final VariantContext variant) {
    // return variant.getGenotype(0).hasDP() ? String.valueOf(variant.getGenotype(0).getDP()): "";
    // }
    // },

    call_GQ { // Required
        public String getColumnValue(final VariantContext variant, final boolean forceLoadingFromNonAlleleSpecific) {
            if (!variant.getGenotype(0).hasGQ()) {
                throw new UserException("Cannot be missing required value for call.GQ");
            }
            return  String.valueOf(variant.getGenotype(0).getGQ());
        }
    },

    call_PGT {
        public String getColumnValue(final VariantContext variant, final boolean forceLoadingFromNonAlleleSpecific) {
            return variant.getGenotype(0).hasAnyAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY) ? String.valueOf(variant.getGenotype(0).getAnyAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY)) : "";
        }
    },

    call_PID {
        public String getColumnValue(final VariantContext variant, final boolean forceLoadingFromNonAlleleSpecific) {
            return variant.getGenotype(0).hasAnyAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY) ? String.valueOf(variant.getGenotype(0).getAnyAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY)) : "";
        }
    },

    call_PS {
        public String getColumnValue(final VariantContext variant, final boolean forceLoadingFromNonAlleleSpecific) {
            return variant.getGenotype(0).hasAnyAttribute(VCFConstants.PHASE_SET_KEY) ? String.valueOf(variant.getGenotype(0).getAnyAttribute(VCFConstants.PHASE_SET_KEY)) : "";
        }
    },

    call_PL {
        public String getColumnValue(final VariantContext variant, final boolean forceLoadingFromNonAlleleSpecific) {
            return variant.getGenotype(0).hasPL() ? Arrays.stream(variant.getGenotype(0).getPL())
                    .mapToObj(String::valueOf)
                    .collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)) : "";
        }
    };


    public String getColumnValue(final VariantContext variant, final boolean forceLoadingFromNonAlleleSpecific) {
        throw new IllegalArgumentException("Not implemented");
    }

    public boolean isVqsrSpecificField() {
        return false;
    }

    public static boolean isGenotypeAllNonRef(Genotype g) {
        // True iff all alleles are non-reference
        return g.getAlleles().stream().allMatch(Allele::isNonReference);
    }

    static final Logger logger = LogManager.getLogger(CreateVariantIngestFiles.class);

    private static String getAttribute(VariantContext vc, String key, String defaultValue){
        Object attr = vc.getAttribute(key);
        if ( attr == null ) return defaultValue;
        if ( attr instanceof String ) return (String)attr;
        if ( attr instanceof List) return StringUtils.join((List)attr, VCFConstants.INFO_FIELD_ARRAY_SEPARATOR);
        return String.valueOf(attr); // throws an exception if this isn't a string
    }
}

