package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.tribble.Feature;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.IntervalCoverageFinder.CandidateCoverageInterval;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.BreakpointEvidence.ExternalEvidence;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.BreakpointEvidence.ReadEvidence;
import org.broadinstitute.hellbender.tools.spark.sv.utils.*;
import org.broadinstitute.hellbender.tools.spark.utils.FlatMapGluer;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchUniqueMultiMapSpark;
import org.broadinstitute.hellbender.utils.SVInterval;
import org.broadinstitute.hellbender.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndexCache;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.spark.SparkUtils;
import scala.Tuple2;

import java.io.*;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * (Internal) Produces local assemblies of genomic regions that may harbor structural variants
 *
 * <p>This tool is used in development and should not be of interest to most researchers.  It packages the
 * identification of genomic regions that might contain structural variation and the generation of local assemblies
 * of these regions as a separate tool, independent of the calling of structural variations from these assemblies.
 * Most researchers will run StructuralVariationDiscoveryPipelineSpark, which both generates local assemblies
 * of interesting genomic regions, and then calls structural variants from these assemblies.</p>
 * <p>This tool identifies genomic regions that may harbor structural variants by integrating evidence from split reads,
 * discordant read pairs, template-length anomalies, and copy-number variation.  It then prepares local assemblies of
 * these regions for structural variant calling.  In addition to the reads that align to these regions, reads sharing
 * kmers (fixed-length subsequences) with the reads aligned in these regions are extracted to produce the local
 * assemblies.</p>
 * <p>The local assemblies are done with <a href="https://github.com/lh3/fermi-lite">FermiLite</a>,
 * and the assembled contigs are aligned to reference with <a href="https://arxiv.org/abs/1303.3997">BWA-MEM</a>.</p>
 * <p>The output is a file of aligned contigs from local assemblies to be used in calling structural variants.</p>
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li>A SAM/BAM/CRAM file of paired-end, aligned and coordinate-sorted reads.</li>
 *     <li>A BWA index image for the reference.
 *         You can use BwaMemIndexImageCreator to create the index image file.</li>
 *     <li>A list of ubiquitous kmers to ignore.
 *         You can use FindBadGenomicGenomicKmersSpark to create the list of kmers to ignore.</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>A file of aligned contigs.</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 * <pre>
 *   gatk FindBreakpointEvidenceSpark \
 *     -I input_reads.bam \
 *     --aligner-index-image reference.img \
 *     --kmers-to-ignore ignored_kmers.txt \
 *     -O assemblies.sam
 * </pre>
 * <p>This tool can be run without explicitly specifying Spark options. That is to say, the given example command
 * without Spark options will run locally. See
 * <a href ="https://software.broadinstitute.org/gatk/documentation/article?id=10060">Tutorial#10060</a>
 * for an example of how to set up and run a Spark tool on a cloud Spark cluster.</p>
 *
 * <h3>Caveats</h3>
 * <p>Expected input is a paired-end, coordinate-sorted BAM with around 30x coverage.
 * Coverage much lower than that probably won't work well.</p>
 */
@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "(Internal) Produces local assemblies of genomic regions that may harbor structural variants",
        summary =
        "This tool is used in development and should not be of interest to most researchers.  It packages the" +
        " identification of genomic regions that might contain structural variation and the generation of local assemblies" +
        " of these regions as a separate tool, independent of the calling of structural variations from these assemblies" +
        " Most researchers will run StructuralVariationDiscoveryPipelineSpark, which both generates local assemblies" +
        " of interesting genomic regions, and then calls structural variants from these assemblies." +
        " This tool identifies genomic regions that may harbor structural variants by integrating evidence from split reads," +
        " discordant read pairs, template-length anomalies, and copy-number variation.  It then prepares local assemblies of" +
        " these regions for structural variant calling.  In addition to the reads that align to these regions, reads sharing" +
        " kmers (fixed-length subsequences) with the reads aligned in these regions are extracted to produce the local" +
        " assemblies." +
        " The local assemblies are done with FermiLite, and the assembled contigs are aligned to reference with BWA-MEM." +
        " The output is a file of aligned contigs from local assemblies to be used in calling structural variants.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class)
public final class FindBreakpointEvidenceSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    // window size for computing per-base depth
    public static final int DEPTH_WINDOW_SIZE = 100000;

    @ArgumentCollection
    private final FindBreakpointEvidenceSparkArgumentCollection params =
            new FindBreakpointEvidenceSparkArgumentCollection();


    @Argument(doc = "sam file for aligned contigs", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputAssemblyAlignments;

    @Override
    public boolean requiresReads()
    {
        return true;
    }

    @Override
    protected void runTool( final JavaSparkContext ctx ) {

        validateParams();

        gatherEvidenceAndWriteContigSamFile(ctx, params, getHeaderForReads(), getUnfilteredReads(),
                outputAssemblyAlignments, logger);

    }

    /**
     * Gathers evidence reads and outputs them in a directory where reads are written out as interleaved FASTQ's.
     * Also produces the SAM records of contigs locally-assembled from such reads for downstream variant discovery.
     *
     * @return the in-memory representation of assembled contigs alignments, whose length equals the number of local assemblies (regardless of success of failure status)
     */
    public static AssembledEvidenceResults gatherEvidenceAndWriteContigSamFile(
            final JavaSparkContext ctx,
            final FindBreakpointEvidenceSparkArgumentCollection params,
            final SAMFileHeader header,
            final JavaRDD<GATKRead> unfilteredReads,
            final String outputAssemblyAlignments,
            final Logger logger) {

        final SVReadFilter filter = new SVReadFilter(params);
        final ReadMetadata readMetadata = buildMetadata(params, header, unfilteredReads, filter, logger);
        log("Metadata retrieved.", logger);

        // develop evidence, intervals, and, finally, a set of template names for each interval
        final EvidenceScanResults
                evidenceScanResults = getMappedQNamesSet(params, readMetadata, ctx, header, unfilteredReads, filter, logger);
        final List<SVInterval> intervals = evidenceScanResults.intervals;
        if ( intervals.isEmpty() ) return new AssembledEvidenceResults(
                evidenceScanResults.readMetadata,
                intervals,
                new ArrayList<>(),
                evidenceScanResults.evidenceTargetLinks);

        final HopscotchUniqueMultiMapSpark<String, Integer, QNameAndInterval> qNamesMultiMap = evidenceScanResults.qNamesForAssemblyMultiMap;

        // supplement the template names with other reads that share kmers
        final List<AlignedAssemblyOrExcuse> alignedAssemblyOrExcuseList;
        if ( params.intervalOnlyAssembly ) {
            alignedAssemblyOrExcuseList = new ArrayList<>();
        } else {
            alignedAssemblyOrExcuseList = addAssemblyQNames(params, readMetadata, ctx, qNamesMultiMap, intervals.size(),
                    unfilteredReads, filter, logger);
        }

        // write a FASTQ file for each interval
        final FermiLiteAssemblyHandler fermiLiteAssemblyHandler =
                new FermiLiteAssemblyHandler(params.alignerIndexImageFile, params.maxFASTQSize,
                                                params.fastqDir, params.writeGFAs,
                                                params.popVariantBubbles, params.removeShadowedContigs,
                                                params.expandAssemblyGraph, params.zDropoff);
        alignedAssemblyOrExcuseList.addAll(
                handleAssemblies(ctx, qNamesMultiMap, unfilteredReads, filter, intervals.size(),
                        params.includeMappingLocation, fermiLiteAssemblyHandler));

        alignedAssemblyOrExcuseList.sort(Comparator.comparingInt(AlignedAssemblyOrExcuse::getAssemblyId));

        // record the intervals
        if ( params.intervalFile != null ) {
            AlignedAssemblyOrExcuse.writeIntervalFile(params.intervalFile, header, intervals, alignedAssemblyOrExcuseList);
        }

        // write alignments of the assembled contigs
        if ( outputAssemblyAlignments != null ) {
            AlignedAssemblyOrExcuse.writeAssemblySAMFile(outputAssemblyAlignments, alignedAssemblyOrExcuseList, header, params.assembliesSortOrder);
            log("Wrote SAM file of aligned contigs.", logger);
        }

        return new AssembledEvidenceResults(evidenceScanResults.readMetadata, intervals, alignedAssemblyOrExcuseList,
                                            evidenceScanResults.evidenceTargetLinks);
    }

    private void validateParams() {
        if( !(outputAssemblyAlignments.endsWith(".bam") || outputAssemblyAlignments.endsWith(".sam")) )
                throw new UserException("Output assembly alignments does not end with \".bam\" or \".sam\": " + outputAssemblyAlignments);

        params.validate();
    }

    public static final class AssembledEvidenceResults {
        final ReadMetadata readMetadata;
        final List<SVInterval> assembledIntervals;
        final List<AlignedAssemblyOrExcuse> alignedAssemblyOrExcuseList;
        final List<EvidenceTargetLink> evidenceTargetLinks;

        public AssembledEvidenceResults(final ReadMetadata readMetadata,
                                        final List<SVInterval> assembledIntervals,
                                        final List<AlignedAssemblyOrExcuse> alignedAssemblyOrExcuseList,
                                        final List<EvidenceTargetLink> evidenceTargetLinks) {
            this.readMetadata = readMetadata;
            this.assembledIntervals = assembledIntervals;
            this.alignedAssemblyOrExcuseList = alignedAssemblyOrExcuseList;
            this.evidenceTargetLinks = evidenceTargetLinks;
        }

        public ReadMetadata getReadMetadata() {
            return readMetadata;
        }

        public List<SVInterval> getAssembledIntervals() { return assembledIntervals; }

        public List<AlignedAssemblyOrExcuse> getAlignedAssemblyOrExcuseList() {
            return alignedAssemblyOrExcuseList;
        }

        public List<EvidenceTargetLink> getEvidenceTargetLinks() {
            return evidenceTargetLinks;
        }
    }

    public static ReadMetadata buildMetadata( final FindBreakpointEvidenceSparkArgumentCollection params,
                                              final SAMFileHeader header,
                                              final JavaRDD<GATKRead> unfilteredReads,
                                              final SVReadFilter filter,
                                              final Logger logger ) {
        Utils.validate(header.getSortOrder() == SAMFileHeader.SortOrder.coordinate,
                "The reads must be coordinate sorted.");

        final Set<Integer> crossContigsToIgnoreSet;
        if ( params.crossContigsToIgnoreFile == null ) crossContigsToIgnoreSet = Collections.emptySet();
        else crossContigsToIgnoreSet = readCrossContigsToIgnoreFile(params.crossContigsToIgnoreFile,
                header.getSequenceDictionary());
        final ReadMetadata readMetadata =
                new ReadMetadata(crossContigsToIgnoreSet, header, params.maxTrackedFragmentLength, unfilteredReads, filter, logger);
        if ( params.metadataFile != null ) {
            ReadMetadata.writeMetadata(readMetadata, params.metadataFile);
        }
        return readMetadata;
    }

    /**
     * Find the breakpoint evidence,
     * cluster the evidence into intervals,
     * find the template names mapped into each interval,
     * kmerize each of these templates,
     * clean up by removing some intervals that are bogus as evidenced by ubiquitous kmers,
     * and return a set of template names and the intervals to which they belong.
     */
    private static EvidenceScanResults getMappedQNamesSet(
            final FindBreakpointEvidenceSparkArgumentCollection params,
            final ReadMetadata readMetadata,
            final JavaSparkContext ctx,
            final SAMFileHeader header,
            final JavaRDD<GATKRead> unfilteredReads,
            final SVReadFilter filter,
            final Logger logger)
    {
        final Broadcast<ReadMetadata> broadcastMetadata = ctx.broadcast(readMetadata);
        final List<List<BreakpointEvidence>> externalEvidence =
                readExternalEvidence(params.externalEvidenceFile, readMetadata,
                                        params.externalEvidenceWeight, params.externalEvidenceUncertainty);
        log("External evidence retrieved.", logger);

        final SVIntervalTree<SVInterval> highCoverageSubintervalTree =
                findGenomewideHighCoverageIntervalsToIgnore(params, readMetadata, ctx, header, unfilteredReads, filter, logger, broadcastMetadata);

        final Broadcast<SVIntervalTree<SVInterval>> broadcastHighCoverageSubIntervals = ctx.broadcast(highCoverageSubintervalTree);

        final Broadcast<List<List<BreakpointEvidence>>> broadcastExternalEvidence = ctx.broadcast(externalEvidence);
        final Tuple2<List<SVInterval>, List<EvidenceTargetLink>> intervalsAndEvidenceTargetLinks =
                getIntervalsAndEvidenceTargetLinks(params, broadcastMetadata, broadcastExternalEvidence, header,
                        unfilteredReads, filter, logger, broadcastHighCoverageSubIntervals);
        List<SVInterval> intervals = intervalsAndEvidenceTargetLinks._1();

        SparkUtils.destroyBroadcast(broadcastExternalEvidence, "external evidence");

        final int nIntervals = intervals.size();
        log("Discovered " + nIntervals + " intervals.", logger);

        if ( nIntervals == 0 ) return new EvidenceScanResults(readMetadata, intervals, intervalsAndEvidenceTargetLinks._2(), null);

        if ( params.exclusionIntervalsFile != null ) {
            intervals = removeIntervalsNearGapsAndLog(intervals, params.exclusionIntervalPadding, readMetadata,
                    params.exclusionIntervalsFile, logger);
        }

        final int nIntervalsAfterGapRemoval = intervals.size();

        // remove any intervals that happen to be completely contained in a high-depth region
        final Iterator<SVInterval> intervalIterator = intervals.iterator();
        while (intervalIterator.hasNext()) {
            final SVInterval interval = intervalIterator.next();
            final Iterator<SVIntervalTree.Entry<SVInterval>> overlappers = highCoverageSubintervalTree.overlappers(interval);
            while (overlappers.hasNext()) {
                final SVIntervalTree.Entry<SVInterval> next = overlappers.next();
                if (next.getInterval().overlapLen(interval) == interval.getLength()) {
                    intervalIterator.remove();
                    break;
                }
            }
        }

        final int nIntervalsAfterDepthCleaning = intervals.size();
        log("Removed " + (nIntervalsAfterGapRemoval - nIntervalsAfterDepthCleaning) + " intervals that were entirely high-depth.", logger);

        final HopscotchUniqueMultiMapSpark<String, Integer, QNameAndInterval> qNamesMultiMap =
                getQNames(params, ctx, broadcastMetadata, intervals, unfilteredReads, filter, broadcastHighCoverageSubIntervals);

        SparkUtils.destroyBroadcast(broadcastHighCoverageSubIntervals, "high-coverage subintervals");
        SparkUtils.destroyBroadcast(broadcastMetadata, "read metadata");

        if ( params.qNamesMappedFile != null ) {
            QNameAndInterval.writeQNames(params.qNamesMappedFile, qNamesMultiMap);
        }
        log("Discovered " + qNamesMultiMap.size() + " mapped template names.", logger);

        return new EvidenceScanResults(readMetadata, intervals, intervalsAndEvidenceTargetLinks._2(), qNamesMultiMap);
    }

    static SVIntervalTree<SVInterval> findGenomewideHighCoverageIntervalsToIgnore(final FindBreakpointEvidenceSparkArgumentCollection params,
                                                                                  final ReadMetadata readMetadata,
                                                                                  final JavaSparkContext ctx,
                                                                                  final SAMFileHeader header,
                                                                                  final JavaRDD<GATKRead> unfilteredReads,
                                                                                  final SVReadFilter filter,
                                                                                  final Logger logger,
                                                                                  final Broadcast<ReadMetadata> broadcastMetadata) {
        final int capacity = header.getSequenceDictionary().getSequences().stream()
                .mapToInt(seqRec -> (seqRec.getSequenceLength() + DEPTH_WINDOW_SIZE - 1)/DEPTH_WINDOW_SIZE).sum();
        final List<SVInterval> depthIntervals = new ArrayList<>(capacity);
        for (final SAMSequenceRecord sequenceRecord : header.getSequenceDictionary().getSequences()) {
            final int contigID = readMetadata.getContigID(sequenceRecord.getSequenceName());
            final int contigLength = sequenceRecord.getSequenceLength();
            for (int i = 1; i < contigLength; i = i + DEPTH_WINDOW_SIZE) {
                depthIntervals.add(new SVInterval(contigID, i, Math.min(contigLength, i + DEPTH_WINDOW_SIZE)));
            }
        }

        final List<SVInterval> highCoverageSubintervals = findHighCoverageSubintervalsAndLog(
                params, ctx, broadcastMetadata, depthIntervals, unfilteredReads, filter, logger);
        final SVIntervalTree<SVInterval> highCoverageSubintervalTree = new SVIntervalTree<>();
        highCoverageSubintervals.forEach(i -> highCoverageSubintervalTree.put(i, i));

        return highCoverageSubintervalTree;
    }

    static final class EvidenceScanResults {
        final ReadMetadata readMetadata;
        final List<SVInterval> intervals;
        final List<EvidenceTargetLink> evidenceTargetLinks;
        final HopscotchUniqueMultiMapSpark<String, Integer, QNameAndInterval> qNamesForAssemblyMultiMap;

        public EvidenceScanResults(final ReadMetadata readMetadata,
                                   final List<SVInterval> intervals,
                                   final List<EvidenceTargetLink> evidenceTargetLinks,
                                   final HopscotchUniqueMultiMapSpark<String, Integer, QNameAndInterval> qNamesForAssemblyMultiMap) {
            this.readMetadata = readMetadata;
            this.intervals = intervals;
            this.evidenceTargetLinks = evidenceTargetLinks;
            this.qNamesForAssemblyMultiMap = qNamesForAssemblyMultiMap;
        }
    }

    /** Read a file of contig names that will be ignored when checking for inter-contig pairs. */
    private static Set<Integer> readCrossContigsToIgnoreFile( final String crossContigsToIgnoreFile,
                                                              final SAMSequenceDictionary dictionary ) {
        final Set<Integer> ignoreSet = new HashSet<>();
        try ( final BufferedReader rdr =
                      new BufferedReader(
                              new InputStreamReader(BucketUtils.openFile(crossContigsToIgnoreFile))) ) {
            String line;
            while ( (line = rdr.readLine()) != null ) {
                final int tigId = dictionary.getSequenceIndex(line);
                if ( tigId == -1 ) {
                    throw new UserException("crossContigToIgnoreFile contains an unrecognized contig name: "+line);
                }
                ignoreSet.add(tigId);
            }
        }
        catch ( final IOException ioe ) {
            throw new UserException("Can't read crossContigToIgnore file "+crossContigsToIgnoreFile, ioe);
        }
        return ignoreSet;
    }

    @VisibleForTesting
    static List<List<BreakpointEvidence>> readExternalEvidence( final String path,
                                                                final ReadMetadata readMetadata,
                                                                final int externalEvidenceWeight,
                                                                final int externalEvidenceUncertainty ) {
        final int nPartitions = readMetadata.getNPartitions();
        final List<List<BreakpointEvidence>> evidenceByPartition = new ArrayList<>(nPartitions);
        for ( int idx = 0; idx != nPartitions; ++idx ) {
            evidenceByPartition.add(new ArrayList<>());
        }

        if ( path == null ) {
            return evidenceByPartition;
        }

        final SVLocation[] partitionBoundaries = new SVLocation[nPartitions+1];
        for ( int idx = 0; idx != nPartitions; ++idx ) {
            final ReadMetadata.PartitionBounds bounds = readMetadata.getPartitionBounds(idx);
            partitionBoundaries[idx] = new SVLocation(bounds.getFirstContigID(), bounds.getFirstStart());
        }
        partitionBoundaries[nPartitions] = new SVLocation(Integer.MAX_VALUE, Integer.MAX_VALUE);
        for ( int idx = 0; idx != nPartitions; ++idx ) {
            if ( partitionBoundaries[idx].compareTo(partitionBoundaries[idx+1]) > 0 ) {
                throw new GATKException("Partition boundaries are not coordinate sorted.");
            }
        }

        final Map<String, Integer> contigNameMap = readMetadata.getContigNameMap();

        SVLocation prevLocation = new SVLocation(0, 0);
        int partitionIdx = 0;
        try ( final FeatureDataSource<Feature> dataSource = new FeatureDataSource<>(path, null, 0, Feature.class) ) {
            for ( final Feature feature : dataSource ) {
                final Integer contigID = contigNameMap.get(feature.getContig());
                if ( contigID == null ) {
                    throw new UserException(path + " contains a contig name not present in the BAM header: " + feature.getContig());
                }
                final int featureStart = feature.getStart();
                final SVLocation featureLocation = new SVLocation(contigID, featureStart);
                if ( prevLocation.compareTo(featureLocation) > 0 ) {
                    throw new UserException("Features in " + path + " are not coordinate sorted.");
                }
                prevLocation = featureLocation;
                while ( featureLocation.compareTo(partitionBoundaries[partitionIdx+1]) >= 0 ) {
                    partitionIdx += 1;
                }
                final List<BreakpointEvidence> partitionEvidence = evidenceByPartition.get(partitionIdx);
                final int featureEnd = feature.getEnd();
                if ( featureEnd - featureStart <= 2*externalEvidenceUncertainty ) {
                    partitionEvidence.add(new ExternalEvidence(contigID, featureStart, featureEnd, externalEvidenceWeight));
                } else {
                    final SVInterval interval1 =
                            BreakpointEvidence.fixedWidthInterval(contigID, featureStart, externalEvidenceUncertainty);
                    partitionEvidence.add(new ExternalEvidence(interval1, externalEvidenceWeight));
                    final SVInterval interval2 =
                            BreakpointEvidence.fixedWidthInterval(contigID, featureEnd, externalEvidenceUncertainty);
                    final SVLocation startLocation2 = interval2.getStartLocation();
                    int partitionIdx2 = partitionIdx;
                    if ( startLocation2.compareTo(partitionBoundaries[partitionIdx+1]) >= 0 ) {
                        partitionIdx2 = Arrays.binarySearch(partitionBoundaries, interval2.getStartLocation());
                        // Binary search returns the complement of the insert position when the exact key is not found.
                        // I.e., ~partitionIdx2 will point to the smallest entry larger than the sought key.
                        // We want the partition before that one, if it exists.
                        if ( partitionIdx2 < 0 ) {
                            partitionIdx2 = ~partitionIdx2;
                            if ( partitionIdx2 > 0 ) partitionIdx2 -= 1;
                        }
                    }
                    evidenceByPartition.get(partitionIdx2).add(new ExternalEvidence(interval2, externalEvidenceWeight));
                }
            }
        }
        return evidenceByPartition;
    }

    /**
     * Kmerize each read mapped into a breakpoint interval,
     * get the template names of all reads sharing these kmers (regardless of where or if they're mapped),
     * and add these template names to the set of names for each interval.
     * Intervals having too many reads are killed.
     * The return is a description (as intervalId and explanatory String) of the intervals that were killed.
     */
    private static List<AlignedAssemblyOrExcuse> addAssemblyQNames(
            final FindBreakpointEvidenceSparkArgumentCollection params,
            final ReadMetadata readMetadata,
            final JavaSparkContext ctx,
            final HopscotchUniqueMultiMapSpark<String, Integer, QNameAndInterval> qNamesMultiMap,
            final int nIntervals,
            final JavaRDD<GATKRead> unfilteredReads,
            final SVReadFilter filter,
            final Logger logger)
    {
        final Tuple2<List<AlignedAssemblyOrExcuse>, HopscotchUniqueMultiMapSpark<SVKmer, Integer, KmerAndInterval>> kmerIntervalsAndDispositions =
                getKmerAndIntervalsSet(params, readMetadata, ctx, qNamesMultiMap, nIntervals,
                                        unfilteredReads, filter, logger);

        final HopscotchUniqueMultiMapSpark<SVKmer, Integer, KmerAndInterval> kmersAndIntervals =
                removeUbiquitousKmers(params, readMetadata, ctx, kmerIntervalsAndDispositions._2(), unfilteredReads, filter, logger);

        qNamesMultiMap.addAll(getAssemblyQNames(params, ctx, kmersAndIntervals, unfilteredReads, filter));

        if ( params.qNamesAssemblyFile != null ) {
            QNameAndInterval.writeQNames(params.qNamesAssemblyFile, qNamesMultiMap);
        }

        log("Discovered "+qNamesMultiMap.size()+" unique template names for assembly.", logger);
        return kmerIntervalsAndDispositions._1();
    }

    /**
     * Kmerize reads having template names in a given set,
     * filter out low complexity kmers and kmers that appear too often in the genome to be helpful in localizing reads,
     * kill intervals that have too few surviving kmers.
     * The return is a Tuple2 in which
     * _1 describes the intervals that have been killed for having too few kmers (as a map from intervalId onto an explanatory string),
     * and _2 describes the good kmers that we want to use in local assemblies (as a multimap from kmer onto intervalId).
     */
    private static Tuple2<List<AlignedAssemblyOrExcuse>, HopscotchUniqueMultiMapSpark<SVKmer, Integer, KmerAndInterval>> getKmerAndIntervalsSet(
            final FindBreakpointEvidenceSparkArgumentCollection params,
            final ReadMetadata readMetadata,
            final JavaSparkContext ctx,
            final HopscotchUniqueMultiMapSpark<String, Integer, QNameAndInterval> qNamesMultiMap,
            final int nIntervals,
            final JavaRDD<GATKRead> unfilteredReads,
            final SVReadFilter filter,
            final Logger logger)
    {
        final Set<SVKmer> kmerKillSet =
                SVFileUtils.readKmersFile(params.kmersToIgnoreFile, params.kSize);
        if ( params.adapterSequence != null ) {
            SVKmerizer.stream(params.adapterSequence, params.kSize, 0, new SVKmerLong())
                    .forEach(kmer -> kmerKillSet.add(kmer.canonical(params.kSize)));
        }
        log("Ignoring " + kmerKillSet.size() + " genomically common kmers.", logger);

        final Tuple2<List<AlignedAssemblyOrExcuse>, List<KmerAndInterval>> kmerIntervalsAndDispositions =
                getKmerIntervals(params, readMetadata, ctx, qNamesMultiMap, nIntervals, kmerKillSet,
                                    unfilteredReads, filter, logger);
        final HopscotchUniqueMultiMapSpark<SVKmer, Integer, KmerAndInterval> kmerMultiMap =
                new HopscotchUniqueMultiMapSpark<>(kmerIntervalsAndDispositions._2());
        log("Discovered " + kmerMultiMap.size() + " kmers.", logger);

        return new Tuple2<>(kmerIntervalsAndDispositions._1(), kmerMultiMap);
    }

    /**
     * Functional interface that consumes the raw materials of an assembly to be aligned (i.e., a Tuple2 of assemblyId
     * and list of sequences) and returns an aligned assembly or an excuse for not producing one.
     */
    @VisibleForTesting interface LocalAssemblyHandler
            extends Serializable, Function<Tuple2<Integer,List<SVFastqUtils.FastqRead>>, AlignedAssemblyOrExcuse> {
    }

    /**
     * Transform all the reads for a supplied set of template names in each interval into FASTQ records
     * for each interval, and do something with the list of FASTQ records for each interval (like write it to a file).
     */
    @VisibleForTesting static List<AlignedAssemblyOrExcuse> handleAssemblies(
            final JavaSparkContext ctx,
            final HopscotchUniqueMultiMapSpark<String, Integer, QNameAndInterval> qNamesMultiMap,
            final JavaRDD<GATKRead> unfilteredReads,
            final SVReadFilter filter,
            final int nIntervals,
            final boolean includeMappingLocation,
            final LocalAssemblyHandler localAssemblyHandler ) {

        final int[] counts = new int[nIntervals];
        for ( final QNameAndInterval qNameAndInterval : qNamesMultiMap ) {
            counts[qNameAndInterval.getIntervalId()] += 1;
        }
        final ComplexityPartitioner partitioner = new ComplexityPartitioner(counts);

        final Broadcast<HopscotchUniqueMultiMapSpark<String, Integer, QNameAndInterval>> broadcastQNamesMultiMap =
                ctx.broadcast(qNamesMultiMap);
        final List<AlignedAssemblyOrExcuse> intervalDispositions =
            unfilteredReads
                .mapPartitionsToPair(readItr ->
                        new ReadsForQNamesFinder(broadcastQNamesMultiMap.value(), nIntervals,
                                includeMappingLocation, readItr, filter).iterator(), false)
                .combineByKey(x -> x,
                                SVUtils::concatenateLists,
                                SVUtils::concatenateLists,
                                partitioner, false, null)
                .map(localAssemblyHandler::apply)
                .collect();

        SparkUtils.destroyBroadcast(broadcastQNamesMultiMap, "QNames multi map");
        BwaMemIndexCache.closeAllDistributedInstances(ctx);

        return intervalDispositions;
    }

    public static final class IntPair {
        private final int int1;
        private final int int2;

        public IntPair( final int int1, final int int2 ) {
            this.int1 = int1;
            this.int2 = int2;
        }

        public int int1() { return int1; }
        public int int2() { return int2; }

        public static IntPair reduce( final IntPair pair1, final IntPair pair2 ) {
            return new IntPair(pair1.int1 + pair2.int1, pair1.int2 + pair2.int2);
        }
    }

    /**
     * For a set of interesting kmers, count occurrences of each over all reads, and remove those
     * that appear too frequently from the set.
     */
    private static HopscotchUniqueMultiMapSpark<SVKmer, Integer, KmerAndInterval> removeUbiquitousKmers(
            final FindBreakpointEvidenceSparkArgumentCollection params,
            final ReadMetadata readMetadata,
            final JavaSparkContext ctx,
            final HopscotchUniqueMultiMapSpark<SVKmer, Integer, KmerAndInterval> kmersAndIntervals,
            final JavaRDD<GATKRead> unfilteredReads,
            final SVReadFilter filter,
            final Logger logger ) {
        final Broadcast<HopscotchUniqueMultiMapSpark<SVKmer, Integer, KmerAndInterval>> broadcastKmersAndIntervals =
                ctx.broadcast(kmersAndIntervals);

        final int kmersPerPartition = kmersAndIntervals.size();
        final int kSize = params.kSize;
        final int partitionSpan = readMetadata.getMedianPartitionSpan();
        // Especially if reads are long and partitions small, the number of partitions in which we observe a kmer will
        // be larger than the actual number of copies at independent loci.
        // This calculation computes the expected inflator given the typical partition size, read length, and kmer size.
        final float partitionsPerLocus = (partitionSpan + readMetadata.getAvgReadLen() - params.kSize) / partitionSpan;
        final int maxPartitionCount = (int)(params.cleanerMaxCopyNumber * partitionsPerLocus);
        logger.info("Cleanup: maxPartitions=" + maxPartitionCount);
        final int maxCount = (int)(params.cleanerMaxCopyNumber * readMetadata.getAccurateKmerCoverage(kSize));
        // The IntPair class is used to keep two counts:
        //  int1 is the number of partitions in which we observe the kmer.
        //  int2 is the total number of observations of the kmer.
        final List<SVKmer> ubiquitousKmers =
                unfilteredReads
                        .filter(filter::notJunk)
                        .filter(filter::isPrimaryLine)
                        .mapPartitions(readItr ->
                            new KmerCounter(kSize,kmersPerPartition,broadcastKmersAndIntervals.getValue()).apply(readItr))
                        .mapToPair(kmerAndCount -> new Tuple2<>(kmerAndCount.getKey(), new IntPair(1, kmerAndCount.getValue())))
                        .reduceByKey(IntPair::reduce)
                        .filter(pair -> {
                            final IntPair intPair = pair._2();
                            return intPair.int1() > maxPartitionCount || intPair.int2() > maxCount;
                        })
                        .map(Tuple2::_1)
                        .collect();

        for ( final SVKmer kmer : ubiquitousKmers ) {
            final Iterator<KmerAndInterval> entryItr = kmersAndIntervals.findEach(kmer);
            while ( entryItr.hasNext() ) {
                entryItr.next();
                entryItr.remove();
            }
        }

        SparkUtils.destroyBroadcast(broadcastKmersAndIntervals, "kmers and intervals");

        log("Removed "+ubiquitousKmers.size()+" ubiquitous kmers.", logger);

        return kmersAndIntervals;
    }

    /**
     * Grab template names for all reads that contain kmers associated with a given breakpoint.
     */
    @VisibleForTesting static List<QNameAndInterval> getAssemblyQNames(
            final FindBreakpointEvidenceSparkArgumentCollection params,
            final JavaSparkContext ctx,
            final HopscotchUniqueMultiMapSpark<SVKmer, Integer, KmerAndInterval> kmerMultiMap,
            final JavaRDD<GATKRead> unfilteredReads,
            final SVReadFilter filter ) {
        final Broadcast<HopscotchUniqueMultiMapSpark<SVKmer, Integer, KmerAndInterval>> broadcastKmersAndIntervals =
                ctx.broadcast(kmerMultiMap);

        final int kSize = params.kSize;
        final List<QNameAndInterval> qNamesAndIntervals =
            unfilteredReads
                .filter(filter::notJunk)
                .filter(filter::isPrimaryLine)
                .mapPartitions(readItr ->
                        new FlatMapGluer<>(new QNameIntervalFinder(kSize,broadcastKmersAndIntervals.getValue()), readItr))
                .collect();

        SparkUtils.destroyBroadcast(broadcastKmersAndIntervals, "cleaned kmers and intervals");

        return qNamesAndIntervals;
    }

    /** find kmers for each interval */
    @VisibleForTesting static Tuple2<List<AlignedAssemblyOrExcuse>, List<KmerAndInterval>> getKmerIntervals(
            final FindBreakpointEvidenceSparkArgumentCollection params,
            final ReadMetadata readMetadata,
            final JavaSparkContext ctx,
            final HopscotchUniqueMultiMapSpark<String, Integer, QNameAndInterval> qNamesMultiMap,
            final int nIntervals,
            final Set<SVKmer> kmerKillSet,
            final JavaRDD<GATKRead> unfilteredReads,
            final SVReadFilter filter,
            final Logger logger ) {

        final Broadcast<Set<SVKmer>> broadcastKmerKillSet = ctx.broadcast(kmerKillSet);
        final Broadcast<HopscotchUniqueMultiMapSpark<String, Integer, QNameAndInterval>> broadcastQNameAndIntervalsMultiMap =
                ctx.broadcast(qNamesMultiMap);

        // given a set of template names with interval IDs and a kill set of ubiquitous kmers,
        // produce a set of interesting kmers for each interval ID
        final int kSize = params.kSize;
        final int kmersPerPartition = (int)(2 * readMetadata.getNRefBases() / readMetadata.getNPartitions());
        // diploid coverage from read metadata is:
        //  {number of sequenced bases}/{size of genome}*{a deflation factor that accounts for the number of erroneous
        //  kmers that will be produced given the mean base quality}
        final float effectiveDiploidCoverage = readMetadata.getAccurateKmerCoverage(kSize);
        // haploid coverage is half the diploid coverage
        final float effectiveHaploidCoverage = effectiveDiploidCoverage / 2;
        // expected distribution of the coverage is expected to be poisson with np (==coverage) large enough
        // that it'll look pretty similar to a normal distribution with mean and variance equal to the coverage.
        final float effectiveHaploidCoverageSD = (float)Math.sqrt(effectiveHaploidCoverage);
        // take mean - 3*sigma as the minimum we expect to see by chance
        final float minCoverage = effectiveHaploidCoverage - 3 * effectiveHaploidCoverageSD;
        // but don't go smaller than the specified min count
        final int minKmers = Math.max(params.cleanerMinKmerCount, Math.round(minCoverage));
        // expected number of observations given the copy number and coverage
        final int maxKmers = Math.round(params.cleanerMaxCopyNumber * effectiveDiploidCoverage);
        logger.info("Cleanup: minKmers=" + minKmers + " maxKmers=" + maxKmers);
        final int maxIntervals = params.cleanerMaxIntervals;
        final int maxDUSTScore = params.maxDUSTScore;
        final List<KmerAndInterval> kmerIntervals =
            unfilteredReads
                .mapPartitionsToPair(readItr ->
                        new FlatMapGluer<>(
                                new QNameKmerizer(
                                        broadcastQNameAndIntervalsMultiMap.value(),
                                        broadcastKmerKillSet.value(), kSize, maxDUSTScore, filter),
                                readItr), false)
                .reduceByKey(Integer::sum)
                .mapPartitions(itr ->
                        new KmerCleaner(itr, kmersPerPartition, minKmers, maxKmers, maxIntervals).iterator())
                .collect();

        SparkUtils.destroyBroadcast(broadcastQNameAndIntervalsMultiMap, "QNames and intervals");
        SparkUtils.destroyBroadcast(broadcastKmerKillSet, "kmer kill set");

        final int[] intervalKmerCounts = new int[nIntervals];
        for ( final KmerAndInterval kmerAndInterval : kmerIntervals ) {
            intervalKmerCounts[kmerAndInterval.getIntervalId()] += 1;
        }
        final Set<Integer> intervalsToKill = new HashSet<>();
        final List<AlignedAssemblyOrExcuse> intervalDispositions = new ArrayList<>();
        for ( int idx = 0; idx != intervalKmerCounts.length; ++idx ) {
            if ( intervalKmerCounts[idx] < params.minKmersPerInterval ) {
                intervalsToKill.add(idx);
                intervalDispositions.add(new AlignedAssemblyOrExcuse(idx,"FASTQ not written -- too few kmers"));
            }
        }

        qNamesMultiMap.removeIf( qNameAndInterval -> intervalsToKill.contains(qNameAndInterval.getIntervalId()) );

        final List<KmerAndInterval> filteredKmerIntervals = kmerIntervals.stream()
                .filter(kmerAndInterval -> !intervalsToKill.contains(kmerAndInterval.getIntervalId()))
                .collect(SVUtils.arrayListCollector(kmerIntervals.size()));

        // record the kmers with their interval IDs
        if ( params.kmerFile != null ) {
            try (final OutputStreamWriter writer = new OutputStreamWriter(new BufferedOutputStream(
                    BucketUtils.createFile(params.kmerFile)))) {
                for (final KmerAndInterval kmerAndInterval : filteredKmerIntervals) {
                    writer.write(kmerAndInterval.toString(kSize) + " " + kmerAndInterval.getIntervalId() + "\n");
                }
            } catch (final IOException ioe) {
                throw new GATKException("Can't write kmer intervals file " + params.kmerFile, ioe);
            }
        }

        return new Tuple2<>(intervalDispositions, filteredKmerIntervals);
    }

    private static List<SVInterval> removeIntervalsNearGapsAndLog( final List<SVInterval> intervals,
                                                                   final int minDistanceToGap,
                                                                   final ReadMetadata readMetadata,
                                                                   final String exclusionIntervalsFile,
                                                                   final Logger logger ) {
        final List<SVInterval> result = removeIntervalsNearGaps(intervals, minDistanceToGap,
                readMetadata.getContigNameMap(), exclusionIntervalsFile);
        final int nKilledIntervals = intervals.size() - result.size();
        log("Killed " + nKilledIntervals + " intervals that were near reference gaps.", logger);
        return result;
    }

    /** remove intervals that are near gaps */
    @VisibleForTesting static List<SVInterval> removeIntervalsNearGaps( final List<SVInterval> intervals,
                                                                        final int minDistanceToGap,
                                                                        final Map<String, Integer> contigNameMap,
                                                                        final String exclusionIntervalsFile ) {
        if ( exclusionIntervalsFile == null ) return intervals;
        final SortedSet<SVInterval> gaps =
                new TreeSet<>(SVFileUtils.readIntervalsFile(exclusionIntervalsFile, contigNameMap));
        return intervals.stream()
                .filter(interval -> {
                    final SortedSet<SVInterval> headSet = gaps.headSet(interval);
                    if ( !headSet.isEmpty() ) {
                        final SVInterval gapInterval = headSet.last();
                        if ( gapInterval.gapLen(interval) < minDistanceToGap ) return false;
                    }
                    final SortedSet<SVInterval> tailSet = gaps.tailSet(interval);
                    if ( !tailSet.isEmpty() ) {
                        final SVInterval gapInterval = tailSet.first();
                        if ( interval.gapLen(gapInterval) < minDistanceToGap ) return false;
                    }
                    return true;
                })
                .collect(Collectors.toCollection(() -> new ArrayList<>(intervals.size())));
    }

    private static List<SVInterval> findHighCoverageSubintervalsAndLog(
            final FindBreakpointEvidenceSparkArgumentCollection params,
            final JavaSparkContext ctx,
            final Broadcast<ReadMetadata> broadcastMetadata,
            final List<SVInterval> intervals,
            final JavaRDD<GATKRead> unfilteredReads,
            final SVReadFilter filter,
            final Logger logger) {

        final int minFlankingHighCovFactor = params.highDepthCoverageFactor;
        final int minPeakHighCovFactor = params.highDepthCoveragePeakFactor;

        final ReadMetadata shortReadMetadata = broadcastMetadata.getValue();
        int minFlankingHighCoverageValue = (int) (minFlankingHighCovFactor * shortReadMetadata.getCoverage());
        int minPeakHighCoverageValue = (int) (minPeakHighCovFactor * shortReadMetadata.getCoverage());
        final List<SVInterval> result =
                findHighCoverageSubIntervals(ctx, broadcastMetadata, intervals, unfilteredReads,
                        filter,
                        minFlankingHighCoverageValue,
                        minPeakHighCoverageValue);
        log("Found " + result.size() + " sub-intervals with coverage over " + minFlankingHighCoverageValue +
                " and a peak coverage of over " + minPeakHighCoverageValue + ".", logger);

        final String intervalFile = params.highCoverageIntervalsFile;
        if (intervalFile != null) {
            try (final OutputStreamWriter writer =
                         new OutputStreamWriter(new BufferedOutputStream(BucketUtils.createFile(intervalFile)))) {
                for (final SVInterval svInterval : result) {
                    final String bedLine = shortReadMetadata.getContigName(svInterval.getContig()) + "\t" + (svInterval.getStart() - 1) + "\t" + svInterval.getEnd() + "\n";
                    writer.write(bedLine);
                }
            } catch (final IOException ioe) {
                throw new UserException.CouldNotCreateOutputFile("Can't write high coverage intervals file " + intervalFile, ioe);
            }
        }
        return result;
    }

    @VisibleForTesting static List<SVInterval> findHighCoverageSubIntervals(
            final JavaSparkContext ctx,
            final Broadcast<ReadMetadata> broadcastMetadata,
            final List<SVInterval> intervals,
            final JavaRDD<GATKRead> unfilteredReads,
            final SVReadFilter filter,
            final int minFlankingHighCoverageValue,
            final int minPeakHighCoverageValue) {
        final Broadcast<List<SVInterval>> broadcastIntervals = ctx.broadcast(intervals);
        final List<CandidateCoverageInterval> highCoverageIntervalRegions =
                unfilteredReads
                    .mapPartitionsToPair(readItr ->
                        new IntervalCoverageFinder(
                                broadcastMetadata.value(),broadcastIntervals.value(),readItr,filter).iterator())
                    .reduceByKey((a, b) -> {
                        final int[] result = new int[a.length];
                        for (int i = 0; i < a.length; i++) {
                            result[i] = a[i] + b[i];
                        }
                        return result;
                    })
                    .flatMap(kv -> findHighCoverageSubintervals(broadcastIntervals, kv, minFlankingHighCoverageValue, minPeakHighCoverageValue))
                        .collect();

        final List<CandidateCoverageInterval> sortedRegions = new ArrayList<>(highCoverageIntervalRegions);
        sortedRegions.sort(Comparator.comparing(CandidateCoverageInterval::getInterval, SVInterval::compareTo));
        // we do two passes of merging together high coverage intervals:
        // first, merge two intervals if (gapLen == 0 and at least one has a max cov peak), and discard unmerged intervals w/o max cov peaks
        //    this takes care of boundary-crossing issues where part of a high-cov interval without a peak is in a different depth window
        // then, merge intervals where (gap len < read length and both have max cov peaks)
        final List<SVInterval> windowBoundaryMergedHighCoverageIntervalRegions = new ArrayList<>(sortedRegions.size());
        if (sortedRegions.size() > 0) {
            CandidateCoverageInterval prev = sortedRegions.get(0);
            for (int i = 1; i < sortedRegions.size(); i++) {
                CandidateCoverageInterval curr = sortedRegions.get(i);

                if (prev.getInterval().gapLen(curr.getInterval()) == 0) {
                    prev = new CandidateCoverageInterval(prev.getInterval().join(curr.getInterval()), prev.containsMaxCoveragePeak() || curr.containsMaxCoveragePeak());
                } else {
                    if (prev.containsMaxCoveragePeak()) {
                        windowBoundaryMergedHighCoverageIntervalRegions.add(prev.getInterval());
                    }
                    prev = curr;
                }
            }
            if (prev.containsMaxCoveragePeak()) {
                windowBoundaryMergedHighCoverageIntervalRegions.add(prev.getInterval());
            }
        }

        final List<SVInterval> gapMergedHighCoverageIntervalRegions = new ArrayList<>(windowBoundaryMergedHighCoverageIntervalRegions.size());
        if (windowBoundaryMergedHighCoverageIntervalRegions.size() > 0) {
            SVInterval prev = windowBoundaryMergedHighCoverageIntervalRegions.get(0);
            for (int i = 1; i < windowBoundaryMergedHighCoverageIntervalRegions.size(); i++) {
                SVInterval curr = windowBoundaryMergedHighCoverageIntervalRegions.get(i);

                if (prev.gapLen(curr) <= broadcastMetadata.getValue().getAvgReadLen()) {
                    prev = prev.join(curr);
                } else {
                    gapMergedHighCoverageIntervalRegions.add(prev);
                    prev = curr;
                }
            }
            gapMergedHighCoverageIntervalRegions.add(prev);
        }

        SparkUtils.destroyBroadcast(broadcastIntervals, "intervals");

        return gapMergedHighCoverageIntervalRegions;
    }

    static Iterator<CandidateCoverageInterval> findHighCoverageSubintervals(final Broadcast<List<SVInterval>> broadcastIntervals,
                                                                                                   final Tuple2<Integer, int[]> intervalCoverage,
                                                                                                   final int minFlankingHighCoverageValue,
                                                                                                   final int minPeakHighCoverageValue) {
        final SVInterval interval = broadcastIntervals.getValue().get(intervalCoverage._1);
        final int[] coverageArray = intervalCoverage._2;
        return IntervalCoverageFinder.getHighCoverageIntervalsInWindow(minFlankingHighCoverageValue, minPeakHighCoverageValue, interval, coverageArray);
    }

    /** find template names for reads mapping to each interval */
    @VisibleForTesting static HopscotchUniqueMultiMapSpark<String, Integer, QNameAndInterval> getQNames(
            final FindBreakpointEvidenceSparkArgumentCollection params,
            final JavaSparkContext ctx,
            final Broadcast<ReadMetadata> broadcastMetadata,
            final List<SVInterval> intervals,
            final JavaRDD<GATKRead> unfilteredReads,
            final SVReadFilter filter,
            final Broadcast<SVIntervalTree<SVInterval>> broadcastHighCoverageSubIntervals) {
        final Broadcast<List<SVInterval>> broadcastIntervals = ctx.broadcast(intervals);
        final List<QNameAndInterval> qNameAndIntervalList =
                unfilteredReads
                    .mapPartitions(readItr ->
                            new FlatMapGluer<>(
                                    new QNameFinder(broadcastMetadata.value(), broadcastIntervals.value(), filter,
                                            broadcastHighCoverageSubIntervals.value()),
                                    readItr), false)
                    .collect();

        SparkUtils.destroyBroadcast(broadcastIntervals, "intervals");

        final HopscotchUniqueMultiMapSpark<String, Integer, QNameAndInterval> qNamesMultiMap =
                new HopscotchUniqueMultiMapSpark<>(params.assemblyToMappedSizeRatioGuess*qNameAndIntervalList.size());
        qNamesMultiMap.addAll(qNameAndIntervalList);
        return qNamesMultiMap;
    }

    /**
     * Identify funky reads that support a hypothesis of a breakpoint in the vicinity, group the reads,
     * and declare a breakpoint interval where there is sufficient density of evidence.
     */
    @VisibleForTesting static Tuple2<List<SVInterval>, List<EvidenceTargetLink>> getIntervalsAndEvidenceTargetLinks(
            final FindBreakpointEvidenceSparkArgumentCollection params,
            final Broadcast<ReadMetadata> broadcastMetadata,
            final Broadcast<List<List<BreakpointEvidence>>> broadcastExternalEvidenceByPartition,
            final SAMFileHeader header,
            final JavaRDD<GATKRead> unfilteredReads,
            final SVReadFilter filter,
            final Logger logger, final Broadcast<SVIntervalTree<SVInterval>> highCoverageSubintervalTree) {
        // find all breakpoint evidence, then filter for pile-ups
        final int nContigs = header.getSequenceDictionary().getSequences().size();
        final int allowedOverhang = params.allowedShortFragmentOverhang;
        final int minEvidenceMapQ = params.minEvidenceMapQ;

        // 1) identify well-mapped reads
        // 2) that look like they support a hypothesis of a breakpoint in the vicinity
        // 3a) filter out those that lack supporting evidence from a sufficient number of other reads, except
        // 3b) pass through everything within a fragment length of partition boundaries
        final JavaRDD<BreakpointEvidence> evidenceRDD = unfilteredReads
                .mapPartitions(readItr -> {
                    final GATKRead sentinel = new SAMRecordToGATKReadAdapter(null);
                    return FlatMapGluer.applyMapFunc(
                            new ReadClassifier(broadcastMetadata.value(), sentinel, allowedOverhang, filter, highCoverageSubintervalTree.getValue()),
                            readItr, sentinel);
                }, true);
        evidenceRDD.cache();

        // record the evidence
        if ( params.unfilteredEvidenceDir != null ) {
            evidenceRDD.map(e -> e.stringRep(broadcastMetadata.getValue(), minEvidenceMapQ))
                       .saveAsTextFile(params.unfilteredEvidenceDir);
        }

        final JavaRDD<EvidenceTargetLink> evidenceTargetLinkJavaRDD = evidenceRDD.mapPartitions(
                itr -> {
                    final ReadMetadata readMetadata = broadcastMetadata.getValue();
                    final EvidenceTargetLinkClusterer clusterer =
                            new EvidenceTargetLinkClusterer(readMetadata, minEvidenceMapQ);
                    return clusterer.cluster(itr);
                }).filter(link -> link.readPairs >= 2 || link.splitReads >= 1);

        final List<EvidenceTargetLink> evidenceTargetLinks = evidenceTargetLinkJavaRDD.collect();

        log("Collected " + evidenceTargetLinks.size() + " evidence target links", logger);

        writeTargetLinks(broadcastMetadata, evidenceTargetLinks, params.targetLinkFile);

        final JavaRDD<BreakpointEvidence> filteredEvidenceRDD =
                evidenceRDD
                    .mapPartitionsWithIndex( (idx, evidenceItr1) -> {
                            final ReadMetadata readMetadata = broadcastMetadata.value();
                            final PartitionCrossingChecker xChecker =
                                    new PartitionCrossingChecker(idx, readMetadata,
                                                                 readMetadata.getMaxMedianFragmentSize());
                            final Iterator<BreakpointEvidence> evidenceItr2 =
                                    broadcastExternalEvidenceByPartition.value().get(idx).iterator();
                            final List<Iterator<BreakpointEvidence>> evidenceItrList = new ArrayList<>(2);
                            evidenceItrList.add(evidenceItr1);
                            evidenceItrList.add(evidenceItr2);
                            final Iterator<BreakpointEvidence> evidenceItr =
                                    FlatMapGluer.concatIterators(evidenceItrList.iterator());
                            return getFilter(evidenceItr,readMetadata, params, xChecker);
                        }, true);

        filteredEvidenceRDD.cache();

        // record the evidence
        if ( params.evidenceDir != null ) {
            filteredEvidenceRDD
                    .filter(BreakpointEvidence::isValidated)
                    .saveAsTextFile(params.evidenceDir);
        }

        // replace clumps of evidence with a single aggregate indicator of evidence (to save memory), except
        // don't clump anything within two fragment lengths of a partition boundary.
        // collect the whole mess in the driver.
        final int maxFragmentSize = broadcastMetadata.value().getMaxMedianFragmentSize();
        final List<BreakpointEvidence> collectedEvidence =
                filteredEvidenceRDD
                        .mapPartitionsWithIndex( (idx, readEvidenceItr) ->
                                new FlatMapGluer<>(
                                        new BreakpointEvidenceClusterer(maxFragmentSize,
                                                new PartitionCrossingChecker(idx,broadcastMetadata.value(),2*maxFragmentSize)),
                                        readEvidenceItr,
                                        new BreakpointEvidence(new SVInterval(nContigs,1,1),0,false)), true)
                        .collect();

        filteredEvidenceRDD.unpersist();
        evidenceRDD.unpersist();

        // reapply the density filter (all data collected -- no more worry about partition boundaries).
        final Iterator<BreakpointEvidence> evidenceIterator =
                getFilter(collectedEvidence.iterator(), broadcastMetadata.value(), params, new PartitionCrossingChecker() );
        final List<BreakpointEvidence> allEvidence = new ArrayList<>(collectedEvidence.size());
        while ( evidenceIterator.hasNext() ) {
            allEvidence.add(evidenceIterator.next());
        }

        // write additional validated read evidence from partition-boundary reads
        if ( params.evidenceDir != null ) {
            final String crossPartitionFile = params.evidenceDir+"/part-xxxxx";
            try ( final OutputStreamWriter writer =
                      new OutputStreamWriter(new BufferedOutputStream(BucketUtils.createFile(crossPartitionFile))) ) {
                for ( final BreakpointEvidence ev : allEvidence ) {
                    // Only tell 'em about the ReadEvidence that was validated in the driver's pass over the stream.
                    // (The validated ReadEvidence instances well away from the partition boundaries that have already
                    // been reported will have been replaced by a generic BreakpointEvidence object that represents
                    // a whole cluster.  A few bits of validated ReadEvidence near the partition boundaries will
                    // get double reported, unfortunately.)
                    if ( ev instanceof ReadEvidence ) {
                        writer.write(ev.toString());
                        writer.write('\n');
                    }
                }
            }
            catch ( final IOException ioe ) {
                throw new GATKException("Can't write cross-partition evidence to "+crossPartitionFile, ioe);
            }
        }

        // re-cluster the new evidence -- we can now glue across partition boundaries
        final Iterator<BreakpointEvidence> evidenceIterator2 =
            new FlatMapGluer<>(new BreakpointEvidenceClusterer(maxFragmentSize,new PartitionCrossingChecker()),
                               allEvidence.iterator(),
                               new BreakpointEvidence(new SVInterval(nContigs,1,1),0,false));

        final List<SVInterval> intervals = new ArrayList<>(allEvidence.size());
        while ( evidenceIterator2.hasNext() ) {
            intervals.add(evidenceIterator2.next().getLocation());
        }

        return new Tuple2<>(intervals, evidenceTargetLinks);
    }

    private static Iterator<BreakpointEvidence> getFilter(
            final Iterator<BreakpointEvidence> evidenceItr,
            final ReadMetadata readMetadata,
            final StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection params,
            final PartitionCrossingChecker partitionCrossingChecker
    ) {
        return new BreakpointDensityFilter(
            evidenceItr, readMetadata, params.minEvidenceWeightPerCoverage,
            params.minCoherentEvidenceWeightPerCoverage, partitionCrossingChecker, params.minEvidenceMapQ
        );
    }

    private static void writeTargetLinks(final Broadcast<ReadMetadata> broadcastMetadata,
                                         final List<EvidenceTargetLink> targetLinks, final String targetLinkFile) {
        if ( targetLinkFile != null ) {
            try (final OutputStreamWriter writer =
                         new OutputStreamWriter(new BufferedOutputStream(BucketUtils.createFile(targetLinkFile)))) {
                targetLinks.iterator().forEachRemaining(entry -> {
                    final String bedpeRecord = entry.toBedpeString(broadcastMetadata.getValue());
                    try {
                        writer.write(bedpeRecord + "\n");
                    } catch (final IOException ioe) {
                        throw new GATKException("Can't write target links to "+ targetLinkFile, ioe);
                    }
                });
            } catch ( final IOException ioe ) {
                throw new GATKException("Can't write target links to "+ targetLinkFile, ioe);
            }
        }
    }

    private static void log(final String message, final Logger logger) {
        logger.info(message);
    }
}
