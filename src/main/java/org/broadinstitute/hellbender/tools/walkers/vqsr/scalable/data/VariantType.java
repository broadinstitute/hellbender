package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

/**
 * Logic for determining variant types was retained from VQSR.
 */
public enum VariantType {
    SNP,
    INDEL;

    public static boolean checkVariantType(final VariantContext vc,
                                           final VariantContext resourceVC) {
        switch (resourceVC.getType()) {
            case SNP:
            case MNP:
                return getVariantType(vc) == SNP;
            case INDEL:
            case MIXED:
            case SYMBOLIC:
                return getVariantType(vc) == INDEL;
            default:
                return false;
        }
    }

    public static VariantType getVariantType(final VariantContext vc) {
        if (vc.isSNP() || vc.isMNP()) {
            return SNP;
        } else if (vc.isStructuralIndel() || vc.isIndel() || vc.isMixed() || vc.isSymbolic()) {
            return INDEL;
        } else {
            throw new IllegalStateException("Encountered unknown variant type: " + vc.getType());
        }
    }

    public static VariantType getVariantType(final VariantContext vc,
                                             final Allele allele) {
        if (vc.getReference().length() == allele.length()) {
            //note that spanning deletions are considered SNPs by this logic
            return SNP;
        } else if ((vc.getReference().length() != allele.length()) || allele.isSymbolic()) {
            return INDEL;
        } else {
            throw new IllegalStateException("Encountered unknown variant type: " + vc.getType());
        }
    }
}
