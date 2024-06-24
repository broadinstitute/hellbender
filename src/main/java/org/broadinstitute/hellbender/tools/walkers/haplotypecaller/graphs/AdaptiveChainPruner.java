package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;
import java.util.stream.Collectors;

public class AdaptiveChainPruner<V extends BaseVertex, E extends BaseEdge> extends ChainPruner<V,E> {
    private final double initialErrorProbability;
    private final double logOddsThreshold;
    private final double seedingLogOddsThreshold;   // threshold for seeding subgraph of good vertices
    private final int maxUnprunedVariants;

    public AdaptiveChainPruner(final double initialErrorProbability, final double logOddsThreshold, final double seedingLogOddsThreshold, final int maxUnprunedVariants) {
        ParamUtils.isPositive(initialErrorProbability, "Must have positive error probability");
        this.initialErrorProbability = initialErrorProbability;
        this.logOddsThreshold = logOddsThreshold;
        this.seedingLogOddsThreshold = seedingLogOddsThreshold;
        this.maxUnprunedVariants = maxUnprunedVariants;
    }

    @Override
    protected Collection<Path<V,E>> chainsToRemove(final List<Path<V, E>> chains) {
        if (chains.isEmpty()) {
            return Collections.emptyList();
        }

        final BaseGraph<V,E> graph = chains.get(0).getGraph();



        Collection<Path<V,E>> probableErrorChains = likelyErrorChains(chains, graph, initialErrorProbability);
        final int errorCount = probableErrorChains.stream().mapToInt(c -> c.getLastEdge().getMultiplicity()).sum();
        final int totalBases = chains.stream().mapToInt(c -> c.getEdges().stream().mapToInt(E::getMultiplicity).sum()).sum();
        final double errorRate = (double) errorCount / totalBases;

        return likelyErrorChains(chains, graph, errorRate).stream().filter(c -> !c.getEdges().stream().anyMatch(BaseEdge::isRef)).collect(Collectors.toList());
    }

    private Collection<Path<V,E>> likelyErrorChains(final List<Path<V, E>> chains, final BaseGraph<V,E> graph, final double errorRate) {
        // pre-compute the left and right log odds of each chain
        final Map<Path<V,E>, Pair<Double, Double>> chainLogOdds = chains.stream()
                .collect(Collectors.toMap(c -> c, c-> chainLogOdds(c, graph, errorRate)));

        // compute correspondence of vertices to incident chains with log odds above the seeding and extending thresholds
        final Multimap<V, Path<V,E>> vertexToSeedableChains = ArrayListMultimap.create();
        final Multimap<V, Path<V,E>> vertexToGoodIncomingChains = ArrayListMultimap.create();
        final Multimap<V, Path<V,E>> vertexToGoodOutgoingChains = ArrayListMultimap.create();

        for (final Path<V,E> chain : chains) {
            if (chainLogOdds.get(chain).getRight() >= logOddsThreshold || chain.getEdges().get(0).isRef()) {
                vertexToGoodIncomingChains.put(chain.getLastVertex(), chain);
            }

            if (chainLogOdds.get(chain).getLeft() >= logOddsThreshold || chain.getEdges().get(0).isRef()) {
                vertexToGoodOutgoingChains.put(chain.getFirstVertex(), chain);
            }

            // seed-worthy chains must pass the more stringent seeding log odds threshold on both sides
            // in addition to that, we only seed from vertices with multiple such chains incoming or outgoing (see below)
            if (chainLogOdds.get(chain).getRight() >= seedingLogOddsThreshold && chainLogOdds.get(chain).getLeft() >= seedingLogOddsThreshold) {
                vertexToSeedableChains.put(chain.getFirstVertex(), chain);
                vertexToSeedableChains.put(chain.getLastVertex(), chain);
            }
        }


        // We have a priority queue of chains to add to the graph, with priority given by the log odds (higher first)
        // for determinism we have a tie breaker based on chains' first vertex
        // Note that chains can we added twice to the queue, once for each side
        final PriorityQueue<Pair<Path<V,E>, Double>> chainsToAdd = new PriorityQueue<>(
                Comparator.comparingDouble((Pair<Path<V,E>, Double> p) -> -p.getRight())
                .thenComparing((Pair<Path<V,E>, Double> p) -> p.getLeft().getFirstVertex().getSequence(), BaseUtils.BASES_COMPARATOR)
                .thenComparing((Pair<Path<V,E>, Double> p) -> p.getLeft().getBases(), BaseUtils.BASES_COMPARATOR)); // Handle rare edge case for non-determinism where two chains with equivalent score start on the same vertex

        // seed the subgraph of good chains by starting with some definitely-good chains.  These include the max-weight chain
        // and chains emanating from vertices with two incoming or two outgoing chains (plus one outgoing or incoming for a total of 3 or more) with good log odds
        // The idea is that a high-multiplicity error chain A that branches into a second error chain B and a continuation-of-the-original-error chain A'
        // may have a high log odds for A'.  However, only in the case of true variation will multiple branches leaving the same vertex have good log odds.
        final Path<V,E> maxWeightChain = getMaxWeightChain(chains);
        chainsToAdd.add(ImmutablePair.of(maxWeightChain, Double.POSITIVE_INFINITY));
        final Set<V> processedVertices = new LinkedHashSet<>(); // vertices whose incident chains have already been enqueued
        for (final V vertex : vertexToSeedableChains.keySet()) {
            if (vertexToSeedableChains.get(vertex).size() > 2) {
                vertexToGoodOutgoingChains.get(vertex).forEach(chain -> chainsToAdd.add(ImmutablePair.of(chain, chainLogOdds.get(chain).getLeft())));
                vertexToGoodIncomingChains.get(vertex).forEach(chain -> chainsToAdd.add(ImmutablePair.of(chain, chainLogOdds.get(chain).getRight())));
                processedVertices.add(vertex);
            }
        }

        final Set<Path<V,E>> goodChains = new LinkedHashSet<>();
        final Set<V> verticesThatAlreadyHaveOutgoingGoodChains = new LinkedHashSet<>();
        final MutableInt variantCount = new MutableInt(0);

        // starting from the high-confidence seed vertices, grow the "good" subgraph along chains with above-threshold log odds,
        // discovering good chains as we go.
        while (!chainsToAdd.isEmpty() && variantCount.intValue() <= maxUnprunedVariants) {
            final Path<V,E> chain = chainsToAdd.poll().getLeft();

            if (!goodChains.add(chain)) {
                continue;
            }

            // When we add an outgoing chain starting from a vertex with other good outgoing chains, we add a variant
            final boolean newVariant = !verticesThatAlreadyHaveOutgoingGoodChains.add(chain.getFirstVertex());
            if (newVariant) {
                variantCount.increment();
            }

            // check whether we've already added this chain from the other side or we've exceeded the variant count limit
            if (newVariant && variantCount.intValue() > maxUnprunedVariants) {
                continue;
            }

            for (final V vertex : Arrays.asList(chain.getFirstVertex(), chain.getLastVertex())) {
                if (!processedVertices.contains(vertex)) {
                    vertexToGoodOutgoingChains.get(vertex).forEach(c -> chainsToAdd.add(ImmutablePair.of(c, chainLogOdds.get(c).getLeft())));
                    vertexToGoodIncomingChains.get(vertex).forEach(c -> chainsToAdd.add(ImmutablePair.of(c, chainLogOdds.get(c).getRight())));
                    processedVertices.add(vertex);
                }
            }
        }

        return chains.stream().filter(c -> !goodChains.contains(c)).collect(Collectors.toList());
    }

    // find the chain containing the edge of greatest weight, taking care to break ties deterministically
    private Path<V, E> getMaxWeightChain(final Collection<Path<V, E>> chains) {
        return chains.stream()
                .max(Comparator.comparingInt((Path<V, E> chain) -> chain.getEdges().stream().mapToInt(BaseEdge::getMultiplicity).max().orElse(0))
                        .thenComparingInt(Path::length)
                        .thenComparing((Path<V,E> c) -> c.getFirstVertex().getSequence(), BaseUtils.BASES_COMPARATOR))
                .get();
    }

    // left and right chain log odds
    private Pair<Double, Double> chainLogOdds(final Path<V,E> chain, final BaseGraph<V,E> graph, final double errorRate) {
        final int leftTotalMultiplicity = MathUtils.sumIntFunction(graph.outgoingEdgesOf(chain.getFirstVertex()), E::getMultiplicity);
        final int rightTotalMultiplicity = MathUtils.sumIntFunction(graph.incomingEdgesOf(chain.getLastVertex()), E::getMultiplicity);

        final int leftMultiplicity = chain.getEdges().get(0).getMultiplicity();
        final int rightMultiplicity = chain.getLastEdge().getMultiplicity();

        final double leftLogOdds = graph.isSource(chain.getFirstVertex()) ? 0.0 :
                Mutect2Engine.logLikelihoodRatio(leftTotalMultiplicity - leftMultiplicity, leftMultiplicity, errorRate);
        final double rightLogOdds = graph.isSink(chain.getLastVertex()) ? 0.0 :
                Mutect2Engine.logLikelihoodRatio(rightTotalMultiplicity - rightMultiplicity, rightMultiplicity, errorRate);

        return ImmutablePair.of(leftLogOdds, rightLogOdds);
    }

}
