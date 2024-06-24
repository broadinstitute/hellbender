package org.broadinstitute.hellbender.utils.downsampling;

import org.apache.commons.lang3.mutable.MutableInt;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.*;


public final class MutectDownsampler extends ReadsDownsampler {
    private static final int SUSPICIOUS_MAPPING_QUALITY = 50;
    private final MutableInt suspiciousReadCount;
    private final int maxSuspiciousReadsPerStride;
    private boolean rejectAllReadsInStride;

    private final int stride;
    private final int maxCoverage;
    private final List<GATKRead> pendingReads;
    private List<GATKRead> finalizedReads;

    private GATKRead firstReadInStride;


    /**
     * @param maxReadsPerAlignmentStart Maximum number of reads per alignment start position. Must be > 0
     * @param stride Length in bases constituting a single pool of reads to downsample
     */
    public MutectDownsampler(final int maxReadsPerAlignmentStart,
                             final int maxSuspiciousReadsPerAlignmentStart,
                             final int stride) {
        // convert coverage per base to coverage per stride
        maxCoverage = maxReadsPerAlignmentStart <= 0 ? Integer.MAX_VALUE : (maxReadsPerAlignmentStart * stride);
        this.stride = ParamUtils.isPositive(stride, "stride must be > 0");
        maxSuspiciousReadsPerStride = maxSuspiciousReadsPerAlignmentStart <= 0 ? Integer.MAX_VALUE : stride * ParamUtils.isPositive(maxSuspiciousReadsPerAlignmentStart, "maxSuspiciousReadsPerAlignmentStart must be > 0");

        pendingReads = new ArrayList<>();
        finalizedReads = new ArrayList<>();
        rejectAllReadsInStride = false;
        suspiciousReadCount = new MutableInt(0);

        clearItems();
        resetStats();
    }

    @Override
    public void submit( final GATKRead newRead ) {
        Utils.nonNull(newRead);
        if (ReadUtils.readHasNoAssignedPosition(newRead)) {
            finalizedReads.add(newRead);
            return;
        } else {
            handlePositionalChange(newRead);

            if (rejectAllReadsInStride) {
                return;
            }

            if (newRead.getMappingQuality() <= SUSPICIOUS_MAPPING_QUALITY) {
                suspiciousReadCount.increment();
            }

            rejectAllReadsInStride |= suspiciousReadCount.intValue() >= maxSuspiciousReadsPerStride;

            pendingReads.add(newRead);
        }
    }

    private void handlePositionalChange( final GATKRead newRead ) {
        if (ReadUtils.readHasNoAssignedPosition(newRead)) {
            return;
        } else if (firstReadInStride == null) {
            firstReadInStride = newRead;
        } else {
            final boolean newContig = !newRead.getAssignedContig().equals(firstReadInStride.getAssignedContig());
            final boolean newStride = newContig || newRead.getAssignedStart() >= firstReadInStride.getAssignedStart() + stride;
            if (newStride) {
                finalizePendingReads();
                firstReadInStride = newRead;
                rejectAllReadsInStride = false;
                suspiciousReadCount.setValue(0);
            }
        }
    }

    private void finalizePendingReads() {
        if (!rejectAllReadsInStride) {
            // most common case: no downsampling necessary
            if (pendingReads.size() <= maxCoverage) {
                finalizedReads.addAll(pendingReads);
            } else {
                // if we exceed the max coverage, just use well-mapped reads.  Maybe the number of such reads won't reach
                // the desired coverage, but if the region is decently mappable the shortfall will be minor.
                final ReservoirDownsampler wellMappedDownsampler = new ReservoirDownsampler(maxCoverage, false);
                pendingReads.stream().filter(read -> read.getMappingQuality() > SUSPICIOUS_MAPPING_QUALITY).forEach(wellMappedDownsampler::submit);
                wellMappedDownsampler.signalEndOfInput();
                final List<GATKRead> readsToFinalize = wellMappedDownsampler.consumeFinalizedItems();
                if (stride > 1) {
                    Collections.sort(readsToFinalize, Comparator.comparingInt(GATKRead::getAssignedStart));
                }
                finalizedReads.addAll(readsToFinalize);
            }
        }
        pendingReads.clear();

    }

    @Override
    public boolean hasFinalizedItems() {
        return ! finalizedReads.isEmpty();
    }

    @Override
    public boolean hasPendingItems() { return !pendingReads.isEmpty(); }

    @Override
    public GATKRead peekFinalized() { return finalizedReads.isEmpty() ? null : finalizedReads.get(0); }

    @Override
    public GATKRead peekPending() { return pendingReads.isEmpty() ? null : pendingReads.get(0); }

    @Override
    public List<GATKRead> consumeFinalizedItems() {
        Collections.sort(finalizedReads, Comparator.comparingInt(GATKRead::getAssignedStart));
        final List<GATKRead> toReturn = finalizedReads;
        finalizedReads = new ArrayList<>();
        return toReturn;
    }

    @Override
    public int size() { return finalizedReads.size() + pendingReads.size(); }

    @Override
    public void signalEndOfInput() {
        finalizePendingReads();
    }

    @Override
    public void clearItems() {
        pendingReads.clear();
        finalizedReads.clear();
        firstReadInStride = null;
        rejectAllReadsInStride = false;
        suspiciousReadCount.setValue(0);
    }

    @Override
    public boolean requiresCoordinateSortOrder() {
        return true;
    }

    @Override
    public void signalNoMoreReadsBefore( final GATKRead read ) {
        handlePositionalChange(read);
    }
}
