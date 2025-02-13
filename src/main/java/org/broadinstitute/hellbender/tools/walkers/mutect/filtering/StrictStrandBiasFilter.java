package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.StrandBiasUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

public class StrictStrandBiasFilter extends HardAlleleFilter {
    private final int minReadsOnEachStrand;

    public StrictStrandBiasFilter(final int minReadsOnEachStrand) {
        this.minReadsOnEachStrand = minReadsOnEachStrand;
    }

    @Override
    public ErrorType errorType() { return ErrorType.ARTIFACT; }

    @Override
    public List<Boolean> areAllelesArtifacts(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext) {
        List<List<Integer>> sbs = StrandBiasUtils.getSBsForAlleles(vc);
        if (minReadsOnEachStrand == 0 || sbs == null || sbs.isEmpty() || sbs.size() <= 1) {
            return Collections.emptyList();
        }
        // remove symbolic alleles
        if (vc.hasSymbolicAlleles()) {
            sbs = GATKVariantContextUtils.removeDataForSymbolicAlleles(vc, sbs);
        }
        // skip the reference
        return sbs.subList(1, sbs.size()).stream().map(altList -> altList.stream().anyMatch(x -> x <= minReadsOnEachStrand)).collect(Collectors.toList());
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.STRICT_STRAND_BIAS_FILTER_NAME;
    }

    @Override
    protected List<String> requiredInfoAnnotations() { return Collections.singletonList(GATKVCFConstants.AS_SB_TABLE_KEY); }
}
