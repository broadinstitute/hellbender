package org.broadinstitute.hellbender.tools.sv.concordance;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVLinkage;

public class SVConcordanceLinkage extends CanonicalSVLinkage<SVCallRecord> {

    public SVConcordanceLinkage(final SAMSequenceDictionary dictionary) {
        super(dictionary, false);
    }

    @Override
    public CanonicalLinkageResult areClusterable(final SVCallRecord a, final SVCallRecord b) {
        final GATKSVVCFConstants.StructuralVariantAnnotationType aType = a.getType();
        final GATKSVVCFConstants.StructuralVariantAnnotationType bType = b.getType();
        // Don't allow CNV/DEL or CNV/DUP matching, which is problematic for concordance calculations
        if ((aType == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV || bType != GATKSVVCFConstants.StructuralVariantAnnotationType.CNV) && aType != bType) {
            return new CanonicalLinkageResult(false);
        }
        return super.areClusterable(a, b);
    }
}
