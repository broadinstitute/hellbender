package org.broadinstitute.hellbender.tools.sv.concordance;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVLinkage;
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

    private final PriorityQueue<LinkageConcordanceRecord> outputBuffer;

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
        outputBuffer = new PriorityQueue<>(new Comparator<>() {
            private final Comparator<SVCallRecord> callComparator = SVCallRecordUtils.getCallComparator(dictionary);
            @Override
            public int compare(LinkageConcordanceRecord o1, LinkageConcordanceRecord o2) {
                return callComparator.compare(o1.record, o2.record);
            }
        });
        truthIdToItemMap = new HashMap<>();
        idToClusterMap = new HashMap<>();
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
    public List<LinkageConcordanceRecord> flush(final boolean force) {
        final List<LinkageConcordanceRecord> collapsedRecords = flushClusters(force).stream()
                .map(c -> new ClosestPair(c.getItemId(), c.getItem(), c.getClosestRecord(), c.getClosestLinkage()))
                .map(p -> new LinkageConcordanceRecord(p.getId(), collapser.apply(p), p.getLinkage()))
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
    private List<LinkageConcordanceRecord> flushBuffer(final boolean force) {
        if (force) {
            final List<LinkageConcordanceRecord> output = new ArrayList<>(outputBuffer.size());
            while (!outputBuffer.isEmpty()) {
                output.add(outputBuffer.poll());
            }
            return output;
        } else {
            final ArrayList<LinkageConcordanceRecord> output = new ArrayList<>();
            final Integer minActiveStartPosition = minActiveStartPosition();
            while (!outputBuffer.isEmpty() &&
                    (minActiveStartPosition == null || outputBuffer.peek().record.getPositionA() <= minActiveStartPosition)) {
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
    public void add(final SVCallRecord item, final Long id, final boolean isTruthVariant) {
        Utils.validateArg(lastItemContig == null || lastItemContig.equals(item.getContigA()), "Attempted to add item on a new contig; please run a force flush beforehand");
        Utils.validateArg(lastItemStart == null || lastItemStart <= item.getPositionA(), "Items must be added in dictionary-sorted order");
        Utils.validateArg(!idToClusterMap.containsKey(id), "ID already in use: " + id);
        lastItemContig = item.getContigA();
        lastItemStart = item.getPositionA();
        if (isTruthVariant) {
            Utils.validateArg(!truthIdToItemMap.containsKey(id), "ID already in use: " + id);
            truthIdToItemMap.put(id, item);
            for (final ActiveClosestPair activeClosestPair : idToClusterMap.values()) {
                final CanonicalSVLinkage.CanonicalLinkageResult result = linkage.areClusterable(activeClosestPair.getItem(), item);
                if (result.getResult()) {
                    activeClosestPair.update(new LinkageConcordanceRecord(id, item, result));
                }
            }
        } else {
            final int maxStart = linkage.getMaxClusterableStartingPosition(item);
            final LinkageConcordanceRecord closest = getClosestItem(item, truthIdToItemMap);
            final ActiveClosestPair cluster = new ActiveClosestPair(id, item, closest, maxStart);
            idToClusterMap.put(id, cluster);
        }
    }

    /**
     * Compares a set of records to the given item and returns the "closest" one.
     */
    @VisibleForTesting
    public LinkageConcordanceRecord getClosestItem(final SVCallRecord evalRecord, final Map<Long, SVCallRecord> candidates) {
        final Comparator<LinkageConcordanceRecord> distanceComparator = Comparator.comparingInt(o -> totalDistance(evalRecord, o.record()));
        final Comparator<LinkageConcordanceRecord> minDistanceComparator = Comparator.comparingInt(o -> minDistance(evalRecord, o.record()));
        final Comparator<LinkageConcordanceRecord> genotypeDistanceComparator = Comparator.comparingInt(o -> genotypeDistance(evalRecord, o.record()));
        // For consistency, in case all other criteria are equal
        final Comparator<LinkageConcordanceRecord> idEqualComparator = Comparator.comparing(o -> !o.record().getId().equals(evalRecord.getId()));
        final Comparator<LinkageConcordanceRecord> idOrderComparator = Comparator.comparing(o -> o.record().getId());
        final Optional<LinkageConcordanceRecord> result = candidates.entrySet().stream()
                .map(other -> new LinkageConcordanceRecord(other.getKey(), other.getValue(), linkage.areClusterable(evalRecord, other.getValue())))
                .filter(o -> o.linkage().getResult())
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

    public String getLastItemContig() {
        return lastItemContig;
    }

    public record ConcordanceRecord (Long id, SVCallRecord record) {}

    public record LinkageConcordanceRecord(Long id, SVCallRecord record, CanonicalSVLinkage.CanonicalLinkageResult linkage) {}

    /**
     * Output container for an evaluation record and its closest truth record.
     */
    public static class ClosestPair {
        final Long id;
        final SVCallRecord evalItem;
        final SVCallRecord closest;
        final CanonicalSVLinkage.CanonicalLinkageResult linkage;

        public ClosestPair(final Long id, final SVCallRecord evalItem, final SVCallRecord closest, final CanonicalSVLinkage.CanonicalLinkageResult linkage) {
            this.id = id;
            this.evalItem = evalItem;
            this.closest = closest;
            this.linkage = linkage;
        }

        /**
         * Constructor with linkage metrics set to null and a "true" result
         */
        public ClosestPair(final Long id, final SVCallRecord evalItem, final SVCallRecord closest) {
            this(id, evalItem, closest, new CanonicalSVLinkage.CanonicalLinkageResult(true));
        }

        public Long getId() {
            return id;
        }

        public SVCallRecord getEvalItem() {
            return evalItem;
        }

        public SVCallRecord getClosest() {
            return closest;
        }

        public CanonicalSVLinkage.CanonicalLinkageResult getLinkage() {
            return linkage;
        }
    }

    /**
     * Internal representation of a eval-truth pair.
     */
    private static class ActiveClosestPair {

        final Long itemId;
        final SVCallRecord item;
        LinkageConcordanceRecord closest;
        final int maxClusterableStartingPosition;

        ActiveClosestPair(final Long itemId, final SVCallRecord item,
                          final LinkageConcordanceRecord closest,
                          final int maxClusterableStartingPosition) {
            this.itemId = Utils.nonNull(itemId);
            this.item = Utils.nonNull(item);
            this.closest = closest;
            this.maxClusterableStartingPosition = maxClusterableStartingPosition;
        }

        /**
         * Compares the given new paired record with the eval record and updates if it's closer than the current one.
         */
        void update(final LinkageConcordanceRecord newClosest) {
            Utils.nonNull(newClosest);
            if (closest == null) {
                closest = newClosest;
            }
            final Comparator<LinkageConcordanceRecord> distanceComparator = Comparator.comparingInt(o -> totalDistance(item, o.record()));
            final Comparator<LinkageConcordanceRecord> minDistanceComparator = Comparator.comparingInt(o -> minDistance(item, o.record()));
            final Comparator<LinkageConcordanceRecord> genotypeDistanceComparator = Comparator.comparingInt(o -> genotypeDistance(item, o.record()));
            final Comparator<LinkageConcordanceRecord> idEqualComparator = Comparator.comparing(o -> !o.record().getId().equals(item.getId()));
            final Comparator<LinkageConcordanceRecord> idOrderComparator = Comparator.comparing(o -> o.record().getId());
            final List<LinkageConcordanceRecord> candidates = Arrays.asList(closest, newClosest);
            closest = candidates.stream().min(distanceComparator
                            .thenComparing(minDistanceComparator)
                            .thenComparing(genotypeDistanceComparator)
                            .thenComparing(idEqualComparator)
                            .thenComparing(idOrderComparator))
                            .get();
        }

        Long getItemId() {
            return itemId;
        }

        SVCallRecord getItem() {
            return item;
        }

        SVCallRecord getClosestRecord() {
            return closest == null ? null : closest.record();
        }

        CanonicalSVLinkage.CanonicalLinkageResult getClosestLinkage() {
            return closest == null ? null : closest.linkage();
        }

        int getMaxClusterableStartingPosition() {
            return maxClusterableStartingPosition;
        }

        String getContig() {
            return item.getContigA();
        }
    }
}
