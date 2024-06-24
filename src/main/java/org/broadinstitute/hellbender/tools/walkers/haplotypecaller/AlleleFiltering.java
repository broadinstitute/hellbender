package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.StrandOddsRatio;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.InverseAllele;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.haplotype.Event;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.jgrapht.graph.DefaultDirectedWeightedGraph;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.io.ComponentNameProvider;
import org.jgrapht.io.DOTExporter;
import org.jgrapht.io.IntegerComponentNameProvider;

import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.google.common.annotations.VisibleForTesting;

/**
 * Filtering haplotypes that contribute weak alleles to the genotyping.
 *
 * @author Ilya Soifer &lt;ilya.soifer@ultimagen.com&gt;
 * @author Yossi Farjoun &lt;farjoun@broadinstitute.org&gt;
 *
 */

public abstract class AlleleFiltering {

    final protected static Logger logger = LogManager.getLogger(AlleleFiltering.class);
    final protected AssemblyBasedCallerArgumentCollection assemblyArgs;
    final private OutputStreamWriter assemblyDebugOutStream;
    AlleleFiltering(final AssemblyBasedCallerArgumentCollection assemblyArgs, final OutputStreamWriter assemblyDebugOutStream){
        this.assemblyArgs = assemblyArgs;
        this.assemblyDebugOutStream = assemblyDebugOutStream;
    }

    /**
     * Finds alleles that are likely not contributing much to explaining the data and remove the haplotypes
     * that contribute them.
     *
     * The alleles from the active region are divided into clusters of alleles that likely "compete" with each
     * other, where compete means that they are the same allele up to a sequencing error although they might
     * be assigned to a different genomic location. In each cluster we iteratively calculate the quality of
     * each allele relative to other alleles in the cluster and remove the allele with the lowest quality.
     * We then also select in each cluster alleles with high SOR and remove them.
     *
     * Every haplotype that contributes a filtered allele is filtered out.
     *
     * @param readLikelihoods unfiltered read x haplotype likelihood matrix
     * @param activeWindowStart location of the active windows (assemblyResult.getPaddedReferenceLoc().getStart()
     * @param suspiciousLocations set of suspicious locations for further marking in genotyping
     * @return Subsetted read x haplotype where only the haplotypes that do not contribute filtered alleles show. Also
     * locations of filtered alleles on the genome added to `suspiciousLocations` list
     */

    public AlleleLikelihoods<GATKRead, Haplotype> filterAlleles(final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
                                                                final int activeWindowStart, Set<Integer> suspiciousLocations){

        logger.debug("SHA:: filter alleles - start");
        final AlleleLikelihoods<GATKRead, Haplotype> subsettedReadLikelihoodsFinal = subsetHaplotypesByAlleles(readLikelihoods, assemblyArgs, activeWindowStart, suspiciousLocations);
        logger.debug("SHA:: filter alleles - end");

        readLikelihoods.setFilteredHaplotypeCount(readLikelihoods.numberOfAlleles() - subsettedReadLikelihoodsFinal.numberOfAlleles());

        if (assemblyDebugOutStream != null) {
            try {
                assemblyDebugOutStream.write("\nThere were " + subsettedReadLikelihoodsFinal.alleles().size() + " haplotypes found after subsetting by alleles. Here they are:\n");
                subsettedReadLikelihoodsFinal.alleles().forEach(h -> {
                    try {
                        assemblyDebugOutStream.write(h.toString());
                        assemblyDebugOutStream.append("\n");
                    } catch (IOException e) {
                        throw new UserException("Error writing to debug output stream", e);
                    }
                });
            } catch (IOException e) {
                throw new UserException("Error writing to debug output stream", e);
            }
        }

        return subsettedReadLikelihoodsFinal;
    }

    /**
     * Main function that filters haplotypes that contribute weak alleles
     * @param readLikelihoods read x haplotype matrix
     * @param assemblyArgs HaplotypeCaller/Mutect2 parameters
     * @param activeWindowStart Genomic location of the start
     * @param suspiciousLocations set of positions where the alleles are being filtered (modified)
     * @return read x haplotype matrix where the filtered haplotypes are removed
     * @throws IOException if output file can't be written
     */
    private AlleleLikelihoods<GATKRead, Haplotype> subsetHaplotypesByAlleles(final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
                                                                         final AssemblyBasedCallerArgumentCollection assemblyArgs,
                                                                         final int activeWindowStart, Set<Integer> suspiciousLocations) {
        // 1. Collect all alleles in the active region
        final Set<Haplotype> disabledHaplotypes = new HashSet<>();
        final Map<Haplotype, Collection<Event>> haplotypeAlleleMap  = new CollectionUtil.DefaultingMap<>((k) -> new ArrayList<>(), true);
        readLikelihoods.alleles().forEach(h -> h.getEventMap().getEvents().stream().forEach(jh -> haplotypeAlleleMap.get(h).add(jh)));

        // 2. Split them into sets to genotype together. The goal is to cluster true allele with all its variants from
        // biased seq. error.
        // The alleles split into clusters of alleles that potentially interact (compete with each other for reads)
        // First we generate a graph with edge for each pair of alleles that do not occur in the same haplotype
        // Then we only keep the edges where the alleles are close or up to hmer indel from each other
        // the connected components of the graph are genotyped together
        final OccurrenceMatrix<Haplotype, Event> occm = new OccurrenceMatrix<>(haplotypeAlleleMap);
        List<Pair<Event, Event>> nonCoOcurringAlleles = occm.nonCoOcurringColumns();
        final List<Pair<Event, Event>> closeNonCoOccurringAlleles = filterByDistance(nonCoOcurringAlleles, 0, 3);
        nonCoOcurringAlleles = filterSameUpToHmerPairs(filterByDistance(nonCoOcurringAlleles,0,20),
                findReferenceHaplotype(readLikelihoods.alleles()), activeWindowStart);
        nonCoOcurringAlleles.addAll(closeNonCoOccurringAlleles);
        final List<Set<Event>> independentAlleles = occm.getIndependentSets(nonCoOcurringAlleles);

        // 3. For each cluster - remove weak alleles
        for (final Set<Event> alleleSet : independentAlleles) {

            // debugging - write the interaction map of the location (we will keep this function from the unused approach
            // where we only attempted to filter alleles that strongly affect an another allele's quality. This approach
            // failed to deliver a significant improvement and thus is not used.
            // interaction map is the graph of how much quality of each allele is improved when another allele is removed
            if (assemblyArgs.writeFilteringGraphs) {
                if (alleleSet.size() > 1 ) {
                    final List<Event> alleleSetAsList = new ArrayList<>(alleleSet);
                    final Map<Event, Integer> initialRPLsMap = new HashMap<>();
                    final DefaultDirectedWeightedGraph<Event, DefaultWeightedEdge> intm =
                            interactionMatrixToGraph(getInteractionMatrix(alleleSetAsList, haplotypeAlleleMap,
                                    readLikelihoods, initialRPLsMap), initialRPLsMap);
                    printInteractionGraph(intm, initialRPLsMap, alleleSet);
                }
            }

            boolean removedAlleles = true;
            final Set<Haplotype> activeHaplotypes = new HashSet<>(readLikelihoods.alleles());

            while (removedAlleles) {
                removedAlleles = false;
                // b. Marginalize: calculate quality of each allele relative to all other alleles
                logger.debug("GAL::start of iteration");
                final List<Event> activeAlleles = activeHaplotypes.stream()
                        .flatMap(h -> h.getEventMap().getEvents().stream().filter(alleleSet::contains))
                        .distinct()
                        .collect(Collectors.toList());;

                final Map<Event, List<Haplotype>> alleleHaplotypeMap = new CollectionUtil.DefaultingMap<>((k) -> new ArrayList<>(), true);
                readLikelihoods.alleles().stream().filter(activeHaplotypes::contains)
                        .forEach(h ->
                                h.getEventMap().getEvents().stream()
                                            .filter(alleleSet::contains)
                                            .forEach(jh -> alleleHaplotypeMap.get(jh).add(h))

                        );

                logger.debug("AHM::printout start");
                for (final Event al : alleleHaplotypeMap.keySet()) {
                    logger.debug("AHM::allele block ---> ");
                    for (final Allele h : alleleHaplotypeMap.get(al)) {
                        logger.debug(() -> String.format("AHM:: (%d) %s/%s: %s", al.getStart(), al.altAllele().getBaseString(), al.refAllele().getBaseString(), h.getBaseString()));
                    }
                    logger.debug("AHM::allele block ---< ");

                }
                logger.debug("AHM::printout end");



                final List<AlleleLikelihoods<GATKRead, Allele>> alleleLikelihoods =
                        activeAlleles.stream().map(al -> getAlleleLikelihoodMatrix(readLikelihoods, al,
                                haplotypeAlleleMap, activeHaplotypes)).collect(Collectors.toList());
                //    c. Calculate SOR and RPL
                // Note that the QUAL is calculated as a PL, that is -10*log likelihood. This means that high PL is low quality allele
                final List<Integer> collectedRPLs = IntStream.range(0, activeAlleles.size()).mapToObj(i -> getAlleleLikelihoodVsInverse(alleleLikelihoods.get(i), activeAlleles.get(i).altAllele())).collect(Collectors.toList());
                final List<Double> collectedSORs = IntStream.range(0, activeAlleles.size()).mapToObj(i -> getAlleleSOR(alleleLikelihoods.get(i), activeAlleles.get(i).altAllele())).collect(Collectors.toList());

                //    d. Generate variants that are below SOR threshold and below RPL threshold
                final List<Event> filteringCandidates = identifyBadAlleles(collectedRPLs,
                                                                                collectedSORs,
                                                                                activeAlleles,
                                                                                assemblyArgs.prefilterQualThreshold,
                                                                                assemblyArgs.prefilterSorThreshold);


                //very weak candidates are filtered out in any case, even if they are alone (they will be filtered anyway even in the GVCF mode)
                // the very weak quality is hardcoded
                final List<Event> filteringCandidatesStringent = identifyBadAlleles(collectedRPLs,
                        collectedSORs,
                        activeAlleles,
                        1,
                        Integer.MAX_VALUE);


                //for now we just mark all locations with close alleles, one of which is weak.
                //We write them in suspiciousLocations and they will be then annotated as SUSP_NOISY... in the VCF
                if ((filteringCandidates.size() > 0 ) && (alleleSet.size()>0)) {
                    activeAlleles.forEach(laa -> suspiciousLocations.add(laa.getStart()));
                }

                //    e. For every variant - calculate what is the effect of its deletion and if higher than threshold - delete and continue

                // (This is a currently disabled code from the approach that would disable only the candidates that strongly
                // affect other alleles
                //Event candidateToDisable = identifyStrongInteractingAllele(filteringCandidates,
                //        hcargs.prefilterQualThreshold, activeAlleles, collectedRPLs, readLikelihoods, haplotypeAlleleMap, alleleHaplotypeMap); )

                // if weak candidate had been identified - add its haplotypes into blacklist, remove the allele from the
                // current cluster and genotype again.
                if ((filteringCandidates.size()>0 && activeAlleles.size()>1) ||
                (activeAlleles.size()==1 && filteringCandidatesStringent.size()>0) ||
                        (filteringCandidates.size()>0 && this.assemblyArgs.filterLoneAlleles)) {

                    if ((filteringCandidatesStringent.size()>0) && (filteringCandidates.size() == 0 )) {
                        throw new GATKException.ShouldNeverReachHereException("The thresholds for stringent allele " +
                                "filtering should always be higher than for the relaxed one");
                    }

                    final Event candidateToDisable = filteringCandidates.get(0);
                    logger.debug(() -> String.format("GAL:: Remove %s", candidateToDisable.toString()));
                    removedAlleles = true;
                    final List<Haplotype> haplotypesToRemove = alleleHaplotypeMap.get(candidateToDisable);
                    disabledHaplotypes.addAll(haplotypesToRemove);
                    activeHaplotypes.removeAll(haplotypesToRemove);
                }
                logger.debug("GAL::end of iteration");

            }
        }

        // finalizing: remove all disabled genotypes
        logger.debug("----- SHA list of removed haplotypes start ----");
        for (Haplotype hap : disabledHaplotypes) {
            logger.debug(() -> String.format("SHA :: Removed haplotype : %s ", hap.toString()));
        }
        logger.debug("----- SHA list of removed haplotypes end ----");

        final Set<Haplotype> eventualAlleles = new HashSet<>();
        readLikelihoods.alleles().stream().filter(al -> !disabledHaplotypes.contains(al)).forEach(eventualAlleles::add);
        logger.debug("----- SHA list of remaining haplotypes start ----");
        for (Haplotype hap : eventualAlleles) {
            logger.debug(() -> String.format("SHA :: Remaining haplotype : %s ", hap.toString()));
        }
        logger.debug("----- SHA list of remaining haplotypes end ----");



        final AlleleLikelihoods<GATKRead, Haplotype> currentReadLikelihoods = readLikelihoods.removeAllelesToSubset(eventualAlleles);
        logger.debug("----- SHA list of remaining alleles start ----");
        final Set<Event> locAllele = new HashSet<>();
        currentReadLikelihoods.alleles().forEach(h -> h.getEventMap().getEvents().stream().forEach(locAllele::add));
        for (final Event al: locAllele) {
            logger.debug(() -> String.format("---- SHA :: %s ", al.toString()));
        }
        logger.debug("----- SHA list of remaining alleles end ----");

        return currentReadLikelihoods;
    }


    /**
     * Finds a list of alleles that are candidate for removal in the order of precedence (first - the best candidate to be removed)
     *
     * @param collectedRPLs list of each allele qualities (collected by {@link AlleleFiltering#getAlleleLikelihoodVsInverse}
     * @param collectedSORs list of each allele SORs (collected by {@link AlleleFiltering#getAlleleSOR(AlleleLikelihoods, Allele)}
     * @param alleles list of alleles in the same order as in collectedRPLs/collectedSORs
     * @param qualThreshold only variants with quality below qualThreshold will be considered
     * @param sorThreshold only variants with SOR above threshold will be considered
     * @return list of alleles that can be removed
     */
    @VisibleForTesting
    List<Event> identifyBadAlleles(final List<Integer> collectedRPLs, final List<Double> collectedSORs,
                                                      final List<Event> alleles,
                                                      final double qualThreshold,
                                                      final double sorThreshold) {

        //collected RPLs are the -10*QUAL of the alleles. high RPL means low quality.
        // SORs are regular: high SOR - strongly biased
        final int[] rplsIndices = getSortedIndexList(collectedRPLs);
        final int[] sorIndices = rplsIndices; // the variants that have high sor are ordered according to their quality


        //this list will contain all alleles that should be filtered in the order of priority
        final List<Event> result = new ArrayList<>();
        final double THRESHOLD = -1 * qualThreshold; // quality threshold is like in GATK (GL) and we collected PL, so QUAL 30 will appear as -30.
        final double SOR_THRESHOLD = sorThreshold;

        //note that high value is a poor quality allele, so the worst allele is the highest collectedRPL
        //we first collect all allleles with low quality: from the lowest
        for (int i = rplsIndices.length-1 ; (i >= 0) && (collectedRPLs.get(rplsIndices[i])>THRESHOLD) ; i--) {
            result.add(alleles.get(rplsIndices[i]));
        }
        int rplCount = result.size();
        //we then add alleles with high SOR. Note that amongh all allleles with the SOR higher than the SOR_THRESHOLD
        //we will first filter the one with the lowest QUAL.
        logger.debug(() -> String.format("SHA:: Have %d candidates with low QUAL", rplCount));
        for (int i = sorIndices.length-1 ; (i >= 0) ; i--) {
             if (collectedSORs.get(sorIndices[i])>SOR_THRESHOLD){
                if (!result.contains(alleles.get(sorIndices[i]))) {
                    result.add(alleles.get(sorIndices[i]));
                }
            }
        }
        logger.debug(() -> String.format("SHA:: Have %d candidates with high SOR", result.size() - rplCount));
        return result;
    }


    /**
     * Generates from read x haplotype matrix a read x allele matrix with two alleles: Allele and ~Allele where Allele
     * is supported by all haplotypes that contain this allele and ~Allele is supported by all other haplotypes.
     * Similar to {@link AlleleLikelihoods#marginalize(Map)}.
     *
     * @param readLikelihoods read x haplotype matrix
     * @param allele allele to consider
     * @param haplotypeAlleleMap map between alleles and haplotypes
     * @param enabledHaplotypes list of haplotypes that are currently inactive
     * @return read x allele matrix
     */
    private AlleleLikelihoods<GATKRead, Allele> getAlleleLikelihoodMatrix(final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
                                                                                     final Event allele,
                                                                                     final Map<Haplotype, Collection<Event>> haplotypeAlleleMap,
                                                                                     final Set<Haplotype> enabledHaplotypes
                                                                                     ){
        final Map<Allele,List<Haplotype>> alleleHaplotypeMap = new CollectionUtil.DefaultingMap<>((k) -> new ArrayList<>(), true);

        final Allele notAllele= InverseAllele.of(allele.altAllele(), true);
        readLikelihoods.alleles().stream().filter(enabledHaplotypes::contains)
                .filter(h->haplotypeAlleleMap.get(h).contains(allele))
                .forEach(alleleHaplotypeMap.get(allele.altAllele())::add);
        readLikelihoods.alleles().stream().filter(enabledHaplotypes::contains)
                .filter(h -> !haplotypeAlleleMap.get(h).contains(allele))
                .forEach(alleleHaplotypeMap.get(notAllele)::add);

        final AlleleLikelihoods<GATKRead, Allele> alleleLikelihoods = readLikelihoods.marginalize(alleleHaplotypeMap);
        logger.debug(() -> String.format("GALM: %s %d %d", allele, alleleHaplotypeMap.get(allele.altAllele()).size(), alleleHaplotypeMap.get(notAllele).size()));
        return alleleLikelihoods;
    }

    //functions to get allele likelihoods and SOR. Differ between the mutect and the HC implementations
    abstract int getAlleleLikelihoodVsInverse(final AlleleLikelihoods<GATKRead, Allele> alleleLikelihoods, final Allele allele);

    private double getAlleleSOR(final AlleleLikelihoods<GATKRead, Allele> alleleLikelihoods, final Allele allele) {
        final Allele notAllele = InverseAllele.of(allele, true);
        final int [][] contingency_table = StrandOddsRatio.getContingencyTableWrtAll(alleleLikelihoods, notAllele, Collections.singletonList(allele), 1);
        final double sor = StrandOddsRatio.calculateSOR(contingency_table);
        logger.debug(() -> String.format("GAS:: %s: %f (%d %d %d %d)", allele.toString(), sor, contingency_table[0][0], contingency_table[0][1], contingency_table[1][0], contingency_table[1][1]));
        return sor;

    }

    //filters pairs of alleles by distance
    private List<Pair<Event, Event>> filterByDistance(
            final List<Pair<Event, Event>> allelePairs,
            final int minDist, final int maxDist) {
        logger.debug(() -> String.format("FBD: input %d pairs ", allelePairs.size()));
        final List<Pair<Event, Event>> result = new ArrayList<>(allelePairs);
        result.removeIf(v -> Math.abs(v.getLeft().getStart() - v.getRight().getStart())>maxDist);
        result.removeIf(v -> Math.abs(v.getLeft().getStart() - v.getRight().getStart())<minDist);
        logger.debug(() -> String.format("FBD: output %d pairs ", allelePairs.size()));

        return result;
    }

    //filters pairs of alleles that are not same up to hmer indel
    private List<Pair<Event, Event>> filterSameUpToHmerPairs(final List<Pair<Event,
            Event>> allelePairs, final Haplotype refHaplotype, final int activeWindowStart) {

        final List<Pair<Event, Event>> result = new ArrayList<>();
        for (final Pair<Event, Event> allelePair: allelePairs) {
            final Pair<Haplotype, Haplotype> modifiedHaplotypes = new ImmutablePair<>(
                    refHaplotype.insertAllele(
                            allelePair.getLeft().refAllele(),
                            allelePair.getLeft().altAllele(),
                            allelePair.getLeft().getStart()),
                    refHaplotype.insertAllele(
                            allelePair.getRight().refAllele(),
                            allelePair.getRight().altAllele(),
                            allelePair.getRight().getStart()));

            if ( BaseUtils.equalUpToHmerChange(modifiedHaplotypes.getLeft().getBases(), modifiedHaplotypes.getRight().getBases()) ) {
                result.add(allelePair);
            }

        }

        return result;
    }


    // find (the) reference haplotype within a list of haplotypes
    static Haplotype findReferenceHaplotype(final List<Haplotype> haplotypeList) {
        for (final Haplotype h: haplotypeList ) {
            if (h.isReference()) {
                return h;
            }
        }
        return null;
    }

    // if alleles are different in length, return their the length of their (potentially) common prefix, otherwise return 0
    private int getCommonPrefixLength(final Allele al1, final Allele al2){
        if (al1.length()!=al2.length()){
            return Math.min(al1.length(), al2.length());
        } else {
            return 0;
        }
    }

    // sort an integer list
    private int[] getSortedIndexList(final List<Integer> values) {
        return IntStream.range(0, values.size()).
                mapToObj(i -> new ImmutablePair<>(i, values.get(i)))
                .sorted(Comparator.comparingInt( v -> (int)v.getRight()))
                .mapToInt(v-> v.getLeft()).toArray();

    }

    // the functions below are currently unused but kept for potential future uses.
    // The goal of these functions is to look at how one allele affects the other and make decisions
    // only for the alleles that really affect others. The approach did not currently work that well
    @SuppressWarnings("unused")
    private Event identifyStrongInteractingAllele(final List<Event> candidateList,
                                                             final float prefilterThreshold,
                                                             final List<Event> allAlleles,
                                                             final List<Integer> rpls,
                                                             final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
                                                             final Map<Haplotype, Collection<Event>> haplotypeAlleleMap,
                                                             final Map<Event, List<Haplotype>> alleleHaplotypeMap
                                                              ){


        logger.debug("ISIA :: start");
        final Map<Event, Integer> initialRPLsMap = new HashMap<>();
        IntStream.range(0, allAlleles.size()).forEach(i -> initialRPLsMap.put(allAlleles.get(i), rpls.get(i)));

        for (final Event cand: candidateList){
            logger.debug(() -> String.format("ISIA :: test %s", cand.toString()));
            if ( initialRPLsMap.get(cand) > (-1)*prefilterThreshold){
                logger.debug( String.format("ISIA:: selected %s due to low QUAL", cand));
                return cand;
            }

            if (allAlleles.size() <=1) {
                return null;
            }

            final Map<Event, Integer> interactionVector = getInteractionVector(cand,
                    haplotypeAlleleMap, alleleHaplotypeMap, readLikelihoods, initialRPLsMap);
            for (final Event allele: interactionVector.keySet()){
                logger.debug(() -> String.format(" --- %s: %d", allele.toString(), initialRPLsMap.get(allele) - interactionVector.get(allele)));
                if (initialRPLsMap.get(allele) - interactionVector.get(allele) > prefilterThreshold ){
                    logger.debug(String.format("ISIA:: selected %s", cand));
                    return cand;
                }
            }
        }
        logger.debug("ISIA :: end");

        return null;

    }


    // function to calculate interactions matrix between the alleles
    private Map<Event, Map<Event, Integer>> getInteractionMatrix(
            final List<Event> alleles,
            final Map<Haplotype, Collection<Event>> haplotypeAlleleMap,
            final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
            final Map<Event, Integer> initialRPLsMap) {

        final Map<Event, List<Haplotype>> alleleHaplotypeMap = new CollectionUtil.DefaultingMap<>((k) -> new ArrayList<>(), true);
        final Set<Haplotype> haplotypes = new HashSet<>(readLikelihoods.alleles());
        readLikelihoods.alleles().stream().forEach(h -> h.getEventMap().getEvents().stream().filter(al -> alleles.contains(al)).forEach(
                        jh -> alleleHaplotypeMap.get(jh).add(h))
                );

        final List<Event> allAlleles = new ArrayList<>(alleleHaplotypeMap.keySet());

        final List<AlleleLikelihoods<GATKRead, Allele>> initialAlleleLikelihoods =
                allAlleles.stream().map(c -> getAlleleLikelihoodMatrix(readLikelihoods, c, haplotypeAlleleMap, haplotypes)).collect(Collectors.toList());

        final List<Integer> initialRPLs = IntStream.range(0, allAlleles.size()).mapToObj(i -> getAlleleLikelihoodVsInverse(initialAlleleLikelihoods.get(i),
                allAlleles.get(i).altAllele())).collect(Collectors.toList());

        for (int i = 0 ; i < allAlleles.size(); i++) {
            initialRPLsMap.put(allAlleles.get(i), initialRPLs.get(i));
        }

        final Map<Event, Map<Event, Integer>> result = new HashMap<>();
        for ( final Event alleleToDisable : allAlleles) {
            Map<Event, Integer> rplsWithoutAlleleMap = getInteractionVector(alleleToDisable, haplotypeAlleleMap, alleleHaplotypeMap, readLikelihoods, initialRPLsMap);
            result.put(alleleToDisable, rplsWithoutAlleleMap);
        }

        return result;
    }

    // function to create interaction of a single allele with other alleles
    private Map<Event, Integer> getInteractionVector(
            final Event alleleToDisable,
            final Map<Haplotype, Collection<Event>> haplotypeAlleleMap,
            final Map<Event, List<Haplotype>> alleleHaplotypeMap,
            final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
            final Map<Event, Integer> initialRPLsMap) {


        final Set<Event> allAlleles = initialRPLsMap.keySet();
        final List<Event> allelesWithoutDisabledAllele = allAlleles.stream().filter(al -> al!=alleleToDisable).collect(Collectors.toList());
        final Set<Haplotype> haplotypes = haplotypeAlleleMap.keySet();
        final Set<Haplotype> haplotypesWithoutDisabledAllele = haplotypes.stream().filter( h -> !alleleHaplotypeMap.get(alleleToDisable).contains(h)).collect(Collectors.toSet());

        final List<AlleleLikelihoods<GATKRead, Allele>> disabledAlleleLikelihood =
                allelesWithoutDisabledAllele.stream().map(c -> getAlleleLikelihoodMatrix(readLikelihoods, c, haplotypeAlleleMap, haplotypesWithoutDisabledAllele)).collect(Collectors.toList());

        final List<Integer> rplsWithoutAllele = IntStream.range(0, allelesWithoutDisabledAllele.size()).mapToObj(i -> getAlleleLikelihoodVsInverse(disabledAlleleLikelihood.get(i),
                allelesWithoutDisabledAllele.get(i).altAllele())).collect(Collectors.toList());

        final Map<Event, Integer> rplsWithoutAlleleMap = new HashMap<>();
        IntStream.range(0, allelesWithoutDisabledAllele.size()).forEach( i -> rplsWithoutAlleleMap.put(allelesWithoutDisabledAllele.get(i), rplsWithoutAllele.get(i)));

        return rplsWithoutAlleleMap;
    }


    private DefaultDirectedWeightedGraph<Event, DefaultWeightedEdge> interactionMatrixToGraph(final Map<Event, Map<Event, Integer>> interactionMatrix,
                                                                                                 final Map<Event, Integer> initialRPL ){
        final DefaultDirectedWeightedGraph<Event, DefaultWeightedEdge> result = new DefaultDirectedWeightedGraph<>(DefaultWeightedEdge.class);
        initialRPL.keySet().stream().forEach(x -> result.addVertex(x));


        for ( final Event loc1 : interactionMatrix.keySet() ) {
            for ( final Event loc2 : interactionMatrix.get(loc1).keySet()){
                final int diff = interactionMatrix.get(loc1).get(loc2) - initialRPL.get(loc2);
                if (diff < 0){
                    final DefaultWeightedEdge edge = result.addEdge(loc1, loc2);
                    result.setEdgeWeight(edge, diff);
                }
            }
        }
        return result;
    }

    //debug function - prints dot file with edges between the alleles that affect each other
    void printInteractionGraph(final DefaultDirectedWeightedGraph<Event, DefaultWeightedEdge>  intm,
                               final Map<Event, Integer> rpls,
                               final Set<Event> alleleSet ){
        final IntegerComponentNameProvider<Event> p1 = new IntegerComponentNameProvider<>();
        final ComponentNameProvider<Event> p2 = (v -> v.toString() + " = " + rpls.get(v)) ;
        final ComponentNameProvider<DefaultWeightedEdge> p4 = (e -> String.valueOf(intm.getEdgeWeight(e)));

        final DOTExporter<Event, DefaultWeightedEdge> dotExporter = new DOTExporter<>(p1, p2, p4,
                null, null);
        final String contig = alleleSet.iterator().next().getContig();
        final int rangeStart = alleleSet.stream().mapToInt(al -> al.getStart()).min().getAsInt();
        final int rangeEnd = alleleSet.stream().mapToInt(al -> al.getStart()).max().getAsInt();
        try {
            final Writer outfile = new FileWriter(String.format("allele.interaction.%s.%d-%d.dot", contig, rangeStart, rangeEnd));
            dotExporter.exportGraph(intm, outfile);
        }
        catch (IOException e) {
            throw new RuntimeException("Unable to write a DOT file" + String.format("allele.interaction.%s.%d-%d.dot", contig, rangeStart, rangeEnd));
        }
    }
}

