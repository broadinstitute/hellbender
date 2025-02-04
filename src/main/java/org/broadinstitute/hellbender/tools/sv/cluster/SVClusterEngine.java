package org.broadinstitute.hellbender.tools.sv.cluster;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.SVLocatable;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * <p>Base class for clustering items that possess start/end genomic coordinates. Efficient algorithms are implemented for
 * single-linkage and max-clique clustering that leverage the low dimensionality of position-based clustering criteria.
 * These algorithms are suitable for clustering problems testing for overlapping events in coordinate-sorted order, with
 * additional possible criteria, when the maximum feasible starting position for an item can be easily estimated (e.g.
 * reciprocal overlap, end-point distance).</p>
 *
 * <p>Based on the method from Marschall T, Costa VG, Canzar S, et al. CLEVER: Clique-enumerating variant finder.
 * Bioinformatics. 2012;28(22):2875-2882.</p>
 *
 * <p>NOTE: precise implementation of {@link SVClusterLinkage#getMaxClusterableStartingPosition(SVLocatable)}
 * is important for efficiency because it determines when a cluster can be finalized and omitted from further clustering tests.</p>
 */
public class SVClusterEngine {

    /**
     * Available clustering algorithms
     */
    public enum CLUSTERING_TYPE {
        SINGLE_LINKAGE,
        MAX_CLIQUE
    }

    private final Function<OutputCluster, SVCallRecord> collapser; // Flattens clusters into a single representative item for output
    private final SVClusterLinkage<SVCallRecord> linkage;
    private Map<Integer, Cluster> idToClusterMap; // Active clusters
    private final Map<Integer, SVCallRecord> idToItemMap; // Active items
    protected final CLUSTERING_TYPE clusteringType;
    private final Comparator<SVCallRecord> itemComparator;

    private String currentContig;
    private int nextItemId;
    private int nextClusterId;
    private int lastStart;
    private Integer minActiveStartingPositionItemId;

    /**
     * @param clusteringType algorithm choice
     * @param collapser function that ingests a collection of clustered items and returns a single representative item
     */
    public SVClusterEngine(final CLUSTERING_TYPE clusteringType,
                           final Function<OutputCluster, SVCallRecord> collapser,
                           final SVClusterLinkage<SVCallRecord> linkage,
                           final SAMSequenceDictionary dictionary) {
        this.clusteringType = clusteringType;
        this.collapser = Utils.nonNull(collapser);
        this.linkage = Utils.nonNull(linkage);
        idToClusterMap = new HashMap<>();
        currentContig = null;
        idToItemMap = new HashMap<>();
        itemComparator = SVCallRecordUtils.getSVLocatableComparator(dictionary);
        nextItemId = 0;
        nextClusterId = 0;
        lastStart = 0;
        minActiveStartingPositionItemId = null;
    }

    @VisibleForTesting
    public Function<OutputCluster, SVCallRecord> getCollapser() {
        return collapser;
    }

    @VisibleForTesting
    public SVClusterLinkage<SVCallRecord> getLinkage() {
        return linkage;
    }

    public SVCallRecord getMinActiveStartingPositionItem() {
        Utils.validate(minActiveStartingPositionItemId == null || idToItemMap.containsKey(minActiveStartingPositionItemId),
                "Unregistered item id " + minActiveStartingPositionItemId);
        return idToItemMap.get(minActiveStartingPositionItemId);
    }

    /**
     * Returns true if there are any active or finalized clusters.
     */
    public final boolean isEmpty() {
        return idToClusterMap.isEmpty();
    }

    /**
     * Adds and clusters the given item. Note that items must be added in order of increasing start position.
     * @param item item to cluster
     */
    public final List<SVCallRecord> addAndFlush(final SVCallRecord item) {
        // Start a new cluster if on a new contig
        if (!item.getContigA().equals(currentContig)) {
            final List<SVCallRecord> result = flush();
            currentContig = item.getContigA();
            lastStart = 0;
            seedCluster(registerItem(item));
            return result;
        } else {
            final int itemId = registerItem(item);
            final List<Integer> clusterIdsToProcess = cluster(itemId);
            return processClusters(clusterIdsToProcess);
        }
    }

    private final int registerItem(final SVCallRecord item) {
        Utils.validate(item.getPositionA() >= lastStart, "Items must be added in order of increasing start coordinate");
        lastStart = item.getPositionA();
        final int itemId = nextItemId++;
        idToItemMap.put(itemId, item);
        if (minActiveStartingPositionItemId == null || item.getPositionA() < getMinActiveStartingPositionItem().getPositionA()) {
            minActiveStartingPositionItemId = itemId;
        }
        return itemId;
    }

    private final int getMaxClusterableStartingPositionByIds(final Collection<Integer> itemIds) {
        Utils.nonNull(itemIds);
        Utils.nonEmpty(itemIds);
        return linkage.getMaxClusterableStartingPosition(itemIds.stream().map(this::getItem).collect(Collectors.toList()));
    }

    /**
     * Add a new {@param <T>} to the current clusters and determine which are complete
     * @param itemId id of registered item to add
     * @return the IDs for clusters that are complete and ready for processing
     */
    private final List<Integer> cluster(final Integer itemId) {
        final SVCallRecord item = getItem(itemId);
        // Get list of item IDs from active clusters that cluster with this item
        final Set<Integer> linkedItems = idToClusterMap.values().stream().map(Cluster::getItemIds)
                .flatMap(List::stream)
                .distinct()
                .filter(other -> !other.equals(itemId) && linkage.areClusterable(item, getItem(other)).getResult())
                .collect(Collectors.toCollection(LinkedHashSet::new));

        // Find clusters to which this item belongs, and which active clusters we're definitely done with
        final List<Integer> clusterIdsToProcess = new ArrayList<>();
        // Clusters to which we simply add the item
        final List<Integer> clustersToAugment = new ArrayList<>();
        // New clusters, formed from subsets of currently active clusters, to which we will add the item
        final Set<List<Integer>> clustersToSeedWith = new HashSet<>();    // Use set to prevent creating duplicate clusters
        for (final Map.Entry<Integer, Cluster> entry : idToClusterMap.entrySet()) {
            final Integer clusterIndex = entry.getKey();
            final Cluster cluster = entry.getValue();
            final List<Integer> clusterItems = cluster.getItemIds();
            if (item.getPositionA() > cluster.getMaxClusterableStart()) {
                clusterIdsToProcess.add(clusterIndex);  //this cluster is complete -- process it when we're done
            } else {
                if (clusteringType.equals(CLUSTERING_TYPE.MAX_CLIQUE)) {
                    final List<Integer> linkedClusterItems = clusterItems.stream().filter(linkedItems::contains).collect(Collectors.toList());
                    final int numLinkedItems = linkedClusterItems.size();
                    if (numLinkedItems == clusterItems.size()) {
                        clustersToAugment.add(clusterIndex);
                    } else if (numLinkedItems > 0) {
                        clustersToSeedWith.add(linkedClusterItems);
                    }
                } else if (clusteringType.equals(CLUSTERING_TYPE.SINGLE_LINKAGE)) {
                    final boolean matchesCluster = clusterItems.stream().anyMatch(linkedItems::contains);
                    if (matchesCluster) {
                        clustersToAugment.add(clusterIndex);
                    }
                } else {
                    throw new IllegalArgumentException("Clustering algorithm for type " + clusteringType.name() + " not implemented");
                }
            }
        }

        // Create new clusters from subsets (max-clique only)
        if (!clustersToSeedWith.isEmpty()) {
            // Currently existing clusters to which we will add the current item
            final List<Set<Integer>> augmentedClusterItemLists = clustersToAugment.stream()
                    .map(idToClusterMap::get)
                    .map(Cluster::getItemIds)
                    .map(HashSet::new)
                    .collect(Collectors.toList());
            final List<Set<Integer>> triggeredClusterItemSets = new ArrayList<>(clustersToSeedWith.size() + augmentedClusterItemLists.size());
            // New clusters formed from subsets of the currently active clusters
            triggeredClusterItemSets.addAll(clustersToSeedWith.stream().map(HashSet::new).collect(Collectors.toList()));
            triggeredClusterItemSets.addAll(augmentedClusterItemLists);
            triggeredClusterItemSets.sort(Comparator.comparingInt(Set::size));
            for (int i = 0; i < triggeredClusterItemSets.size(); i++) {
                final Set<Integer> seedItems = triggeredClusterItemSets.get(i);
                // Check that this cluster is not a sub-cluster of any of the others being created
                boolean isSubset = false;
                for (int j = i + 1; j < triggeredClusterItemSets.size(); j++) {
                    if (triggeredClusterItemSets.get(j).containsAll(seedItems)) {
                        isSubset = true;
                        break;
                    }
                }
                if (!isSubset) {
                    seedWithExistingCluster(itemId, seedItems);
                }
            }
        }

        // Add to or merge existing clusters
        if (clusteringType.equals(CLUSTERING_TYPE.SINGLE_LINKAGE)) {
            if (!clustersToAugment.isEmpty()) {
                combineClusters(clustersToAugment, itemId);
            }
        } else {
            for (final Integer clusterId : clustersToAugment) {
                addToCluster(clusterId, itemId);
            }
        }
        
        // If there weren't any matches, create a new singleton cluster
        if (clustersToAugment.isEmpty() && clustersToSeedWith.isEmpty()) {
            seedCluster(itemId);
        }
        return clusterIdsToProcess;
    }

    /**
     * Creates a new cluster by agglomerating clusters with the given ids together (and deleting them from the currently
     * active set), along with an additional item.
     * @param clusterIds ids of clusters to combine
     * @param itemId id of item to add to new cluster
     */
    private final void combineClusters(final Collection<Integer> clusterIds, final Integer itemId) {
        final List<Cluster> clusters = clusterIds.stream().map(this::getCluster).collect(Collectors.toList());
        clusterIds.stream().forEach(idToClusterMap::remove);
        final List<Integer> clusterItems = clusters.stream()
                .map(Cluster::getItemIds)
                .flatMap(List::stream)
                .distinct()
                .collect(Collectors.toList());
        final List<Integer> newClusterItems = new ArrayList<>(clusterItems.size() + 1);
        newClusterItems.addAll(clusterItems);
        newClusterItems.add(itemId);
        idToClusterMap.put(nextClusterId++, new Cluster(getMaxClusterableStartingPositionByIds(newClusterItems), newClusterItems));
    }

    /**
     * Finalizes a single cluster, removing it from the currently active set and adding it to the output buffer.
     */
    private final SVCallRecord processCluster(final int clusterIndex) {
        final Cluster cluster = getCluster(clusterIndex);
        idToClusterMap.remove(clusterIndex);
        final List<Integer> clusterItemIds = cluster.getItemIds();
        final OutputCluster outputCluster = new OutputCluster(clusterItemIds.stream().map(idToItemMap::get).collect(Collectors.toList()));
        final SVCallRecord result = collapser.apply(outputCluster);
        // Clean up item id map
        if (clusterItemIds.size() == 1) {
            // Singletons won't be present in any other clusters
            idToItemMap.remove(clusterItemIds.get(0));
        } else {
            // Need to check that items aren't present in any other clusters
            final Set<Integer> activeItemIds = idToClusterMap.values().stream()
                    .map(Cluster::getItemIds)
                    .flatMap(List::stream)
                    .collect(Collectors.toSet());
            final List<Integer> itemsToRemove = idToItemMap.keySet().stream().filter(i -> !activeItemIds.contains(i))
                    .collect(Collectors.toList());
            for (final Integer i : itemsToRemove) {
                idToItemMap.remove(i);
            }
        }
        // Update min active start position
        if (clusterItemIds.contains(minActiveStartingPositionItemId)) {
            findAndSetMinActiveStart();
        }
        return result;
    }

    /**
     * Scans active items for the current min active starting position.
     */
    private final void findAndSetMinActiveStart() {
        minActiveStartingPositionItemId = null;
        SVCallRecord minActiveStartingPositionItem = null;
        for (final Integer itemId : idToItemMap.keySet()) {
            final SVCallRecord item = idToItemMap.get(itemId);
            if (minActiveStartingPositionItemId == null || itemComparator.compare(item, minActiveStartingPositionItem) < 0) {
                minActiveStartingPositionItemId = itemId;
                minActiveStartingPositionItem = idToItemMap.get(itemId);
            }
        }
    }

    /**
     * Finalizes a set of clusters.
     */
    private final List<SVCallRecord> processClusters(final List<Integer> clusterIdsToProcess) {
        final List<SVCallRecord> result = new ArrayList<>(clusterIdsToProcess.size());
        for (final Integer clusterId : clusterIdsToProcess) {
            result.add(processCluster(clusterId));
        }
        return result;
    }

    /**
     * Finalizes all active clusters and adds them to the output buffer. Also clears the currently active set of clusters
     * and items.
     */
    public final List<SVCallRecord> flush() {
        final List<Integer> clustersToFlush = new ArrayList<>(idToClusterMap.keySet());
        final List<SVCallRecord> result = new ArrayList<>(clustersToFlush.size());
        for (final Integer clusterId : clustersToFlush) {
            result.add(processCluster(clusterId));
        }
        idToItemMap.clear();
        minActiveStartingPositionItemId = null;
        nextItemId = 0;
        nextClusterId = 0;
        return result;
    }

    /**
     * Creates a new singleton cluster containing the given item.
     */
    private final void seedCluster(final Integer item) {
        final List<Integer> newClusters = new ArrayList<>(1);
        newClusters.add(item);
        idToClusterMap.put(nextClusterId++, new Cluster(linkage.getMaxClusterableStartingPosition(getItem(item)), newClusters));
    }

    /**
     * Create a new cluster containing a new item and a set of existing items, and add it to the set of active clusters.
     * Note that if there exists an identical activate cluster, it will not be added.
     * @param item new item (assumed registered)
     * @param seedItems existing items
     */
    private final void seedWithExistingCluster(final Integer item, final Collection<Integer> seedItems) {
        final List<Integer> newClusterItems = new ArrayList<>(1 + seedItems.size());
        newClusterItems.addAll(seedItems);
        newClusterItems.add(item);
        idToClusterMap.put(nextClusterId++,
                new Cluster(getMaxClusterableStartingPositionByIds(newClusterItems), newClusterItems));
    }

    private final Cluster getCluster(final int id) {
        Utils.validateArg(idToClusterMap.containsKey(id), "Cluster ID " + id + " does not exist.");
        return idToClusterMap.get(id);
    }

    private final SVCallRecord getItem(final int id) {
        Utils.validateArg(idToItemMap.containsKey(id), "Item ID " + id + " does not exist.");
        return idToItemMap.get(id);
    }

    /**
     * Add the item specified by {@param itemId} to the cluster specified by {@param clusterIndex}
     * and expand the clustering interval
     * @param clusterId
     * @param itemId
     */
    private final void addToCluster(final int clusterId, final Integer itemId) {
        final Cluster cluster = getCluster(clusterId);
        final List<Integer> clusterItems = cluster.getItemIds();
        clusterItems.add(itemId);
        final SVCallRecord item = getItem(itemId);
        final int itemClusterableStartPosition = linkage.getMaxClusterableStartingPosition(item);
        cluster.setMaxClusterableStart(Math.max(cluster.getMaxClusterableStart(), itemClusterableStartPosition));
    }

    public static final class OutputCluster {
        final List<SVCallRecord> items;
        public OutputCluster(final List<SVCallRecord> items) {
            this.items = items;
        }

        public List<SVCallRecord> getItems() {
            return items;
        }
    }

    /**
     * Container class for clustered items
     */
    private static final class Cluster {
        private int maxClusterableStart;
        private final List<Integer> itemIds;

        public Cluster(final int maxClusterableStart, final List<Integer> itemIds) {
            Utils.nonNull(itemIds);
            this.maxClusterableStart = maxClusterableStart;
            this.itemIds = itemIds;
        }

        public int getMaxClusterableStart() {
            return maxClusterableStart;
        }

        public void setMaxClusterableStart(final int position) {
            maxClusterableStart = position;
        }

        public List<Integer> getItemIds() {
            return itemIds;
        }

        /**
         * Note we do not check for equality on max clusterable start position, which could be dependent on the
         * state of the engine.
         */
        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof Cluster)) return false;
            Cluster cluster = (Cluster) o;
            return itemIds.equals(cluster.itemIds);
        }

        @Override
        public int hashCode() {
            return Objects.hash(itemIds);
        }
    }
}
