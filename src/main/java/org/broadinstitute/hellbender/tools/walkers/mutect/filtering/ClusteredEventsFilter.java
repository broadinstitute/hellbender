package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.List;

public class ClusteredEventsFilter extends HardFilter {
    private final int maxEventsInRegion;
    private final int maxEventsInHaplotype;

    public ClusteredEventsFilter(final int maxEventsInRegion, final int maxEventsInHaplotype) {
        this.maxEventsInRegion = maxEventsInRegion;
        this.maxEventsInHaplotype = maxEventsInHaplotype;
    }

    @Override
    public ErrorType errorType() { return ErrorType.ARTIFACT; }

    @Override
    public boolean isArtifact(final VariantContext vc, final Mutect2FilteringEngine filteringEngine) {
        final List<Integer> haplotypeEventCounts = vc.getAttributeAsIntList(GATKVCFConstants.EVENT_COUNT_IN_HAPLOTYPE_KEY, 0);
        final int regionEventCounts = vc.getAttributeAsInt(GATKVCFConstants.EVENT_COUNT_IN_REGION_KEY, 0);
        return haplotypeEventCounts.stream().mapToInt(n -> n).max().getAsInt() > maxEventsInHaplotype || regionEventCounts > maxEventsInRegion;
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.CLUSTERED_EVENTS_FILTER_NAME;
    }

    @Override
    protected List<String> requiredInfoAnnotations() { return List.of(GATKVCFConstants.EVENT_COUNT_IN_REGION_KEY, GATKVCFConstants.EVENT_COUNT_IN_HAPLOTYPE_KEY); }
}
