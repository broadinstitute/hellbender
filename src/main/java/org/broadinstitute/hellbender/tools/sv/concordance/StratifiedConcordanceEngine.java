package org.broadinstitute.hellbender.tools.sv.concordance;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.cluster.SVClusterEngineArgumentsCollection;
import org.broadinstitute.hellbender.tools.sv.stratify.SVStratificationEngine;
import org.broadinstitute.hellbender.tools.sv.stratify.SVStratificationEngineArgumentsCollection;
import org.broadinstitute.hellbender.tools.walkers.sv.SVStratify;
import org.broadinstitute.hellbender.tools.walkers.validation.Concordance;
import org.broadinstitute.hellbender.tools.walkers.validation.ConcordanceState;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

public class StratifiedConcordanceEngine {

    protected final Map<String, ClosestSVFinder> clusterEngineMap;
    protected final SVStratificationEngine stratificationEngine;
    protected final Map<Long, ItemTracker> variantStatusMap = new HashMap<>();
    protected List<VariantContext> outputBuffer;
    protected ClosestSVFinder defaultEngine;
    protected SAMSequenceDictionary dictionary;
    protected final SVStratificationEngineArgumentsCollection stratArgs;
    protected Long nextItemId = 0L;

    public StratifiedConcordanceEngine(final Map<String, ClosestSVFinder> clusterEngineMap,
                                       final SVStratificationEngine stratificationEngine,
                                       final SVStratificationEngineArgumentsCollection stratArgs,
                                       final SVClusterEngineArgumentsCollection defaultClusteringArgs,
                                       final SVConcordanceAnnotator defaultCollapser,
                                       final SAMSequenceDictionary dictionary) {
        Utils.validate(stratificationEngine.getStrata().size() == clusterEngineMap.size(),
                "Stratification and clustering configurations have a different number of groups.");
        for (final SVStratificationEngine.Stratum stratum : stratificationEngine.getStrata()) {
            Utils.validate(clusterEngineMap.containsKey(stratum.getName()),
                    "Could not find group " + stratum.getName() + " in clustering configuration.");
        }

        this.clusterEngineMap = clusterEngineMap;
        this.stratificationEngine = stratificationEngine;
        this.stratArgs = stratArgs;
        this.outputBuffer = new ArrayList<>();
        this.dictionary = dictionary;

        final SVConcordanceLinkage defaultLinkage = new SVConcordanceLinkage(dictionary);
        defaultLinkage.setDepthOnlyParams(defaultClusteringArgs.getDepthParameters());
        defaultLinkage.setMixedParams(defaultClusteringArgs.getMixedParameters());
        defaultLinkage.setEvidenceParams(defaultClusteringArgs.getPESRParameters());
        this.defaultEngine = new ClosestSVFinder(defaultLinkage, defaultCollapser::annotate, false, dictionary);
    }

    public Collection<VariantContext> flush(boolean force) {
        for (final String name : clusterEngineMap.keySet()) {
            flushEngineToBuffer(clusterEngineMap.get(name), force, name);
        }
        flushEngineToBuffer(defaultEngine, force, SVStratify.DEFAULT_STRATUM);
        final Collection<VariantContext> result = outputBuffer;
        outputBuffer = new ArrayList<>();
        return result;
    }

    public boolean isEmpty() {
        return variantStatusMap.isEmpty() && outputBuffer.isEmpty();
    }

    public void addTruthVariant(final SVCallRecord record) {
        final Long id = nextItemId++;
        // truth variants can be matched across all stratification groups
        for (final String name : clusterEngineMap.keySet()) {
            addToEngine(record, id, true, clusterEngineMap.get(name), name);
        }
        addToEngine(record, id, true, defaultEngine, SVStratify.DEFAULT_STRATUM);
    }

    public void addEvalVariant(final SVCallRecord record) {
        final Long id = nextItemId++;
        // eval variants get stratified
        final List<String> stratifications = stratificationEngine.getMatches(record,
                        stratArgs.overlapFraction, stratArgs.numBreakpointOverlaps, stratArgs.numBreakpointOverlapsInterchrom)
                .stream().map(SVStratificationEngine.Stratum::getName).collect(Collectors.toUnmodifiableList());
        for (final String stratum : stratifications) {
            Utils.validate(clusterEngineMap.containsKey(stratum), "Group undefined: " + stratum);
            addToEngine(record, id, false, clusterEngineMap.get(stratum), stratum);
        }
        // default stratum
        addToEngine(record, id, false, defaultEngine, SVStratify.DEFAULT_STRATUM);
        // Keep track of which groups it was added to
        Utils.validate(!variantStatusMap.containsKey(id), "Attempted to add duplicate id");
        final List<String> stratList = new ArrayList<>(stratifications.size() + 1);
        stratList.addAll(stratifications);
        stratList.add(SVStratify.DEFAULT_STRATUM);
        variantStatusMap.put(id, new ItemTracker(stratList));
    }

    /**
     * Adds a record to the given concordance engine, flushing if hit a new contig
     */
    protected void addToEngine(final SVCallRecord record, final Long id, final boolean isTruth, final ClosestSVFinder engine, final String name) {
        if (engine.getLastItemContig() != null && !record.getContigA().equals(engine.getLastItemContig())) {
            flushEngineToBuffer(engine, true, name);
        }
        engine.add(record, id, isTruth);
        flushEngineToBuffer(engine, false, name);
    }

    protected void flushEngineToBuffer(final ClosestSVFinder engine, final boolean force, final String name) {
        for (final ClosestSVFinder.ConcordanceRecord record: engine.flush(force)) {
            final ItemTracker tracker = variantStatusMap.get(record.id());
            Utils.validate(tracker != null, "Unregistered variant id: " + record.id());
            tracker.eject(name, record.record());
            if (tracker.allEjected()) {
                Utils.nonNull(tracker.getOutput(), "No tracker output");
                variantStatusMap.remove(record.id());
                final List<String> strata = tracker.getGroups().stream().sorted().collect(Collectors.toUnmodifiableList());
                final VariantContextBuilder builder = SVCallRecordUtils.getVariantBuilder(tracker.getOutput())
                        .attribute(GATKSVVCFConstants.STRATUM_INFO_KEY, strata);
                outputBuffer.add(builder.make());
            }
        }
    }

    protected static class ItemTracker {

        private final Map<String, Integer> groupPriorityMap;
        private SVCallRecord out;
        private Integer outPriority;
        private List<String> groups;

        public ItemTracker(final List<String> submittedTo) {
            this.groups = submittedTo;
            this.groupPriorityMap = new HashMap<>(SVUtils.hashMapCapacity(submittedTo.size()));
            for (int i = 0; i < submittedTo.size(); i++) {
                this.groupPriorityMap.put(submittedTo.get(i), i);
            }
        }

        public void eject(final String name, final SVCallRecord record) {
            final Integer priority = groupPriorityMap.remove(name);
            Utils.nonNull(priority, "Unregistered name: " + name);
            if (outPriority == null || (isTruePositive(record)) && (priority < outPriority || (priority > outPriority && !isTruePositive(out)))) {
                out = record;
                outPriority = priority;
            }
        }

        private boolean isTruePositive(final SVCallRecord record) {
            return ConcordanceState.TRUE_POSITIVE.getAbbreviation().equals(record.getAttributes().get(Concordance.TRUTH_STATUS_VCF_ATTRIBUTE));
        }

        public boolean allEjected() {
            return groupPriorityMap.isEmpty();
        }

        public SVCallRecord getOutput() {
            return out;
        }

        public List<String> getGroups() {
            return groups;
        }
    }
}
