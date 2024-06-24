package org.broadinstitute.hellbender.tools.sv.concordance;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Efficiently clusters a set of evaluation ("eval") SVs with their closest truth SVs.
 *
 * "Closest" is defined in the {@link #getClosestItem} method and selects first on total breakpoint distance
 * (sum of both ends). As a tiebreaker, it then considers the distance of the closest breakend, number of
 * matching genotypes, and finally the variant ID.
 *
 * The clustering algorithm is based on that of the {@link org.broadinstitute.hellbender.tools.sv.cluster.SVClusterEngine}
 * which makes efficient use of memory by dumping variants which are out of range of the current genomic position
 * and therefore irrelevant, assuming records are fed in sorted order.
 *
 * Furthermore, completed variants are held in a second buffer than ensures outputs are produced in sorted order as well.
 *
 * The output should be flushed frequently to avoid excessive memory usage, and a "forced" flush must be performed
 * between contigs. Developers should see {@link org.broadinstitute.hellbender.tools.walkers.sv.SVConcordance} as an
 * example use case.
 */
public class ClosestSVFinder {

    protected final boolean sortOutput;
    protected final Map<Long, SVCallRecord> truthIdToItemMap;
    protected final Map<Long, ActiveClosestPair> idToClusterMap;
    private final SVConcordanceLinkage linkage;
    private final Function<ClosestPair, SVCallRecord> collapser;

    private final PriorityQueue<SVCallRecord> outputBuffer;

    private Long nextItemId;
    private Integer lastItemStart;
    private String lastItemContig;

    /**
     * @param linkage linkage to use for matching non-truth and truth variants
     * @param collapser output collapsing function
     * @param dictionary reference dict
     */
    public ClosestSVFinder(final SVConcordanceLinkage linkage,
                           final Function<ClosestPair, SVCallRecord> collapser,
                           final boolean sortOutput,
                           final SAMSequenceDictionary dictionary) {
        this.sortOutput = sortOutput;
        this.linkage = Utils.nonNull(linkage);
        this.collapser = Utils.nonNull(collapser);
        outputBuffer = new PriorityQueue<>(SVCallRecordUtils.getCallComparator(dictionary));
        truthIdToItemMap = new HashMap<>();
        idToClusterMap = new HashMap<>();
        nextItemId = 0L;
        lastItemStart = null;
        lastItemContig = null;
    }

    /**
     * Sorts output by default
     */
    public ClosestSVFinder(final SVConcordanceLinkage linkage,
                           final Function<ClosestPair, SVCallRecord> collapser,
                           final SAMSequenceDictionary dictionary) {
        this(linkage, collapser, true, dictionary);
    }

    /**
     * Should be run frequently to reduce memory usage. Forced flushes must be run when a new contig is encountered.
     * @param force flushes all variants in the output buffer regardless of position
     * @return finalized variants
     */
    public List<SVCallRecord> flush(final boolean force) {
        final List<SVCallRecord> collapsedRecords = flushClusters(force).stream()
                .map(c -> new ClosestPair(c.getItem(), c.getClosest()))
                .map(collapser)
                .collect(Collectors.toList());
        if (sortOutput) {
            outputBuffer.addAll(collapsedRecords);
            return flushBuffer(force);
        } else {
            return collapsedRecords;
        }
    }

    /**
     * Flushes output buffer
     */
    private List<SVCallRecord> flushBuffer(final boolean force) {
        if (force) {
            final List<SVCallRecord> output = new ArrayList<>(outputBuffer.size());
            while (!outputBuffer.isEmpty()) {
                output.add(outputBuffer.poll());
            }
            return output;
        } else {
            final ArrayList<SVCallRecord> output = new ArrayList<>();
            final Integer minActiveStartPosition = minActiveStartPosition();
            while (!outputBuffer.isEmpty() &&
                    (minActiveStartPosition == null || outputBuffer.peek().getPositionA() <= minActiveStartPosition)) {
                output.add(outputBuffer.poll());
            }
            output.trimToSize();
            return output;
        }
    }

    private Integer minActiveStartPosition() {
        return idToClusterMap.isEmpty() ? null : idToClusterMap.values().stream().mapToInt(c -> c.getItem().getPositionA()).min().getAsInt();
    }

    /**
     * Flushes active clusters
     */
    private List<ActiveClosestPair> flushClusters(final boolean force) {
        if (force) {
            final List<ActiveClosestPair> output = new ArrayList<>(idToClusterMap.values());
            truthIdToItemMap.clear();
            idToClusterMap.clear();
            lastItemStart = null;
            lastItemContig = null;
            return output;
        } else {
            // Remove finalized ref items
            truthIdToItemMap.values().removeIf(v -> linkage.getMaxClusterableStartingPosition(v) < lastItemStart);
            // Find and remove finalized clusters
            final List<Map.Entry<Long, ActiveClosestPair>> finalizedClusters = idToClusterMap.entrySet().stream()
                    .filter(e -> e.getValue().getMaxClusterableStartingPosition() < lastItemStart)
                    .collect(Collectors.toList());
            finalizedClusters.forEach(e -> idToClusterMap.remove(e.getKey()));
            return finalizedClusters.stream().map(Map.Entry::getValue).collect(Collectors.toList());
        }
    }

    /**
     * Adds and clusters a new variant. Variants must be added in dictionary-sorted order.
     */
    public void add(final SVCallRecord item, final boolean isTruthVariant) {
        Utils.validateArg(lastItemContig == null || lastItemContig.equals(item.getContigA()), "Attempted to add item on a new contig; please run a force flush beforehand");
        Utils.validateArg(lastItemStart == null || lastItemStart <= item.getPositionA(), "Items must be added in dictionary-sorted order");
        lastItemContig = item.getContigA();
        lastItemStart = item.getPositionA();
        if (isTruthVariant) {
            truthIdToItemMap.put(nextItemId, item);
            idToClusterMap.values().stream()
                .filter(other -> linkage.areClusterable(other.getItem(), item))
                .forEach(cluster -> cluster.update(nextItemId, item));
        } else {
            final int maxStart = linkage.getMaxClusterableStartingPosition(item);
            final Map.Entry<Long, SVCallRecord> minEntry = getClosestItem(item, truthIdToItemMap.entrySet());
            final Long closestId = minEntry == null ? null : minEntry.getKey();
            final SVCallRecord closest = minEntry == null ? null : minEntry.getValue();
            final ActiveClosestPair cluster = new ActiveClosestPair(nextItemId, item, closestId, closest, maxStart);
            idToClusterMap.put(nextItemId, cluster);
        }
        nextItemId++;
    }

    /**
     * Compares a set of records to the given item and returns the "closest" one.
     */
    @VisibleForTesting
    public Map.Entry<Long, SVCallRecord> getClosestItem(final SVCallRecord evalRecord, final Set<Map.Entry<Long, SVCallRecord>> candidates) {
        final Comparator<Map.Entry<Long, SVCallRecord>> distanceComparator = Comparator.comparingInt(o -> totalDistance(evalRecord, o.getValue()));
        final Comparator<Map.Entry<Long, SVCallRecord>> minDistanceComparator = Comparator.comparingInt(o -> minDistance(evalRecord, o.getValue()));
        final Comparator<Map.Entry<Long, SVCallRecord>> genotypeDistanceComparator = Comparator.comparingInt(o -> genotypeDistance(evalRecord, o.getValue()));
        // For consistency, in case all other criteria are equal
        final Comparator<Map.Entry<Long, SVCallRecord>> idEqualComparator = Comparator.comparing(o -> !o.getValue().getId().equals(evalRecord.getId()));
        final Comparator<Map.Entry<Long, SVCallRecord>> idOrderComparator = Comparator.comparing(o -> o.getValue().getId());
        final Optional<Map.Entry<Long, SVCallRecord>> result = candidates.stream()
                .filter(other -> linkage.areClusterable(evalRecord, other.getValue()))
                .min(distanceComparator.thenComparing(minDistanceComparator).thenComparing(genotypeDistanceComparator)
                        .thenComparing(idEqualComparator).thenComparing(idOrderComparator));
        return result.orElseGet(() -> null);
    }

    /**
     * Total distance between breakends, or {@link Integer#MAX_VALUE} if one record is null. Asserts
     * that non-null records have the same start and end contigs.
     */
    public static int totalDistance(final SVCallRecord a, final SVCallRecord b) {
        if (a == null || b == null) {
            return Integer.MAX_VALUE;
        } else {
            Utils.validate(a.getContigA().equals(b.getContigA()), "Start is on different contigs");
            Utils.validate(a.getContigB().equals(b.getContigB()), "End is on different contigs");
            return Math.abs(a.getPositionA() - b.getPositionA()) + Math.abs(a.getPositionB() - b.getPositionB());
        }
    }

    /**
     * Distance between closest breakends, or {@link Integer#MAX_VALUE} if one record is null. Asserts
     * that non-null records have the same start and end contigs.
     */
    public static int minDistance(final SVCallRecord a, final SVCallRecord b) {
        if (a == null || b == null) {
            return Integer.MAX_VALUE;
        } else {
            Utils.validate(a.getContigA().equals(b.getContigA()), "Start is on different contigs");
            Utils.validate(a.getContigB().equals(b.getContigB()), "End is on different contigs");
            return Math.min(Math.abs(a.getPositionA() - b.getPositionA()), Math.abs(a.getPositionB() - b.getPositionB()));
        }
    }

    /**
     * Genotype distance, defined as the negative fraction of matching genotypes, or {@link Integer#MAX_VALUE} if one
     * record is null. Asserts that non-null records have the same start and end contigs.
     */
    public static int genotypeDistance(final SVCallRecord a, final SVCallRecord b) {
        if (a == null || b == null) {
            return Integer.MAX_VALUE;
        } else {
            final GenotypesContext genotypesA = a.getGenotypes();
            final GenotypesContext genotypesB = b.getGenotypes();
            int matches = 0;
            for (final Genotype genotypeA : genotypesA) {
                final Genotype genotypeB = genotypesB.get(genotypeA.getSampleName());
                if (genotypeB != null && genotypeA.sameGenotype(genotypeB)) {
                    matches++;
                }
            }
            return -matches;
        }
    }

    /**
     * Output container for an evaluation record and its closest truth record.
     */
    public static class ClosestPair {
        final SVCallRecord evalItem;
        final SVCallRecord closest;
        public ClosestPair(final SVCallRecord evalItem, final SVCallRecord closest) {
            this.evalItem = evalItem;
            this.closest = closest;
        }

        public SVCallRecord getEvalItem() {
            return evalItem;
        }

        public SVCallRecord getClosest() {
            return closest;
        }
    }

    /**
     * Internal representation of a eval-truth pair.
     */
    private static class ActiveClosestPair {

        final Long itemId;
        final SVCallRecord item;
        Long closestId;
        SVCallRecord closest;
        final int maxClusterableStartingPosition;

        ActiveClosestPair(final Long itemId, final SVCallRecord item,
                          final Long closestId, final SVCallRecord closest,
                          final int maxClusterableStartingPosition) {
            this.itemId = Utils.nonNull(itemId);
            this.item = Utils.nonNull(item);
            this.closestId = closestId;
            this.closest = closest;
            this.maxClusterableStartingPosition = maxClusterableStartingPosition;
        }

        /**
         * Compares the given new paired record with the eval record and updates if it's closer than the current one.
         */
        boolean update(final Long newClosestId, final SVCallRecord newClosest) {
            Utils.nonNull(newClosest);
            final Comparator<Map.Entry<Long, SVCallRecord>> distanceComparator = Comparator.comparingInt(o -> totalDistance(item, o.getValue()));
            final Comparator<Map.Entry<Long, SVCallRecord>> minDistanceComparator = Comparator.comparingInt(o -> minDistance(item, o.getValue()));
            final Comparator<Map.Entry<Long, SVCallRecord>> genotypeDistanceComparator = Comparator.comparingInt(o -> genotypeDistance(item, o.getValue()));
            final Comparator<Map.Entry<Long, SVCallRecord>> idEqualComparator = Comparator.comparing(o -> !o.getValue().getId().equals(item.getId()));
            final Comparator<Map.Entry<Long, SVCallRecord>> idOrderComparator = Comparator.comparing(o -> o.getValue().getId());
            final List<Map.Entry<Long, SVCallRecord>> candidates = new ArrayList<>(2);
            candidates.add(new AbstractMap.SimpleImmutableEntry<>(closestId, closest));
            candidates.add(new AbstractMap.SimpleImmutableEntry<>(newClosestId, newClosest));
            final Map.Entry<Long, SVCallRecord> result = candidates.stream()
                    .min(distanceComparator.thenComparing(minDistanceComparator).thenComparing(genotypeDistanceComparator)
                            .thenComparing(idEqualComparator).thenComparing(idOrderComparator)).get();
            if (result.getKey() == closestId) {
                return false;
            }
            closestId = Utils.nonNull(newClosestId);
            closest = newClosest;
            return true;
        }

        Long getItemId() {
            return itemId;
        }

        SVCallRecord getItem() {
            return item;
        }

        Long getClosestId() {
            return closestId;
        }

        SVCallRecord getClosest() {
            return closest;
        }

        int getMaxClusterableStartingPosition() {
            return maxClusterableStartingPosition;
        }

        String getContig() {
            return item.getContigA();
        }
    }
}
