package org.broadinstitute.hellbender.tools.walkers.readorientation;

import com.google.common.primitives.Ints;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.metrics.StringHeader;
import htsjdk.samtools.util.Histogram;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class F1R2CountsCollector {

    public static final String ALT_TABLE_EXTENSION = ".alt_table";
    public static final String ALT_HIST_EXTENSION = ".alt_histogram";
    public static final String REF_HIST_EXTENSION = ".ref_histogram";

    private static final Logger logger = LogManager.getLogger(F1R2CountsCollector.class);

    private final Set<String> samples;

    private final SAMFileHeader header;

    private final CollectF1R2CountsArgumentCollection CF1R2Args;

    // For each sample and for reference context, count ref sites in a histogram keyed by depth
    private final Map<String, Map<String, Histogram<Integer>>> refSiteHistograms;

    // For each sample store the total depths of alt sites with alt depth = 1 separately to save memory
    private Map<String, DepthOneHistograms> depthOneAltHistograms;

    // alt table writers for each sample and 3-base context
    // This is by far the largest quantity of data (for each sample histograms involve 64 3-base contexts x 3
    // substitution errors x 2 orientations x (by default) 200 histogram count bins, which comes out to maybe 1MB of RAM.
    // In a genome of 3 billion bases where every hundredth site or so has more than 1 alt base the alt records could come out
    // to a gigabyte or more.  Thus it makes sense to split the alt table writers by sample *and* by 3-base context
    // so that downstream tools (LearnReadOrientationModel) need only load the records for one context at a time.
    private Map<String, Map<String, AltSiteRecord.AltSiteRecordTableWriter>> altTableWriters;

    // We put all pre-tarred output files in a temporary directory.  These include three files for each sample:
    // an alt records file produced by that sample's altTableWriter, a depth one alt histogram file, and a ref site
    // histogram file.  Each one of these files contains information for all 3-base contexts.
    private final File tmpDir = IOUtils.createTempDir("untarred");

    // the final output file is produced by archiving every file in the temporary directory
    private final File outputTarGzFile;

    public F1R2CountsCollector(final CollectF1R2CountsArgumentCollection CF1R2Args, final SAMFileHeader header, final File outputTarGzFile, final Collection<String> samples) {
        this.CF1R2Args = CF1R2Args;
        this.samples = samples.size() == 1 ? Collections.singleton(samples.iterator().next()) : new HashSet<>(samples);
        this.header = header;
        this.outputTarGzFile = outputTarGzFile;

        refSiteHistograms = new HashMap<>(samples.size());
        depthOneAltHistograms = new HashMap<>(samples.size());
        altTableWriters = new HashMap<>(samples.size());

        for (final String sample : samples) {
            final Map<String, Histogram<Integer>> refSitesHistogramsForThisSample = new HashMap<>(F1R2FilterConstants.ALL_KMERS.size());

            altTableWriters.put(sample, new HashMap<>());
            // Initialize for each reference the histogram of the counts of reference sites by depth
            for (final String context : F1R2FilterConstants.ALL_KMERS) {
                Histogram<Integer> emptyRefHistogram = F1R2FilterUtils.createRefHistogram(context, CF1R2Args.maxDepth);
                refSitesHistogramsForThisSample.put(context, emptyRefHistogram);

                // Intentionally not use try-with-resources so that the writer stays open outside of the try block
                final File altTableFile = new File(tmpDir, IOUtils.urlEncode(sample) + "." + context + ALT_TABLE_EXTENSION);
                try {
                    altTableWriters.get(sample).put(context, new AltSiteRecord.AltSiteRecordTableWriter(IOUtils.fileToPath(altTableFile), sample));
                } catch (IOException e) {
                    throw new UserException(String.format("Encountered an IO exception creating a writer for %s", altTableFile), e);
                }
            }

            refSiteHistograms.put(sample, refSitesHistogramsForThisSample);
            depthOneAltHistograms.put(sample, new DepthOneHistograms(CF1R2Args.maxDepth));
        }
    }

    public void process(final ReadPileup pileup, final ReferenceContext referenceContext) {
        final int position = referenceContext.getInterval().getStart();
        final String refContext = referenceContext.getKmerAround(position, F1R2FilterConstants.REF_CONTEXT_PADDING);
        if (refContext == null) {
            return;
        }
        final Nucleotide refBase = F1R2FilterUtils.getMiddleBase(refContext);

        if (refContext.contains("N") || refContext.length() != F1R2FilterConstants.REFERENCE_CONTEXT_SIZE) {
            return;
        }

        if (refContext == null) {
            logger.warn(String.format("Skipped a site with null reference at interval %s, k-mer = %s",
                    referenceContext.getInterval().toString(), refContext));
            return;
        }

        // optimize the common case of a single tumor sample
        final String onlySample = samples.size() == 1 ? samples.iterator().next() : null;
        final Map<String, ReadPileup> splitPileup = samples.size() == 1 ? Collections.singletonMap(onlySample,
                        pileup.makeFilteredPileup(pe -> Objects.equals(ReadUtils.getSampleName(pe.getRead(), header), onlySample) && pe.getQual() > CF1R2Args.minBaseQuality))
                : pileup.splitBySample(header, null);

        for (final Map.Entry<String, ReadPileup> entry : splitPileup.entrySet()) {
            final String sample = entry.getKey();

            // if size == 1 we filtered by base quality already
            final ReadPileup samplePileup = samples.size() == 1 ? entry.getValue() :
                    entry.getValue().makeFilteredPileup(pe -> pe.getQual() > CF1R2Args.minBaseQuality);
            final int[] baseCounts = samplePileup.getBaseCounts();
            final int depth = (int) MathUtils.sum(baseCounts);

            if (!isPileupGood(samplePileup)) {
                return;
            }

            // Make a copy of base counts and update the counts of ref to -1. Now the maxElementIndex of the array gives us
            // the alt base.
            final int[] baseCountsCopy = Arrays.copyOf(baseCounts, baseCounts.length);
            baseCountsCopy[refBase.ordinal()] = -1;
            final int altBaseIndex = MathUtils.maxElementIndex(baseCountsCopy);
            final boolean referenceSite = baseCounts[altBaseIndex] == 0;

            // If the site is ref, we simply update the coverage histogram
            if (referenceSite) {
                refSiteHistograms.get(sample).get(refContext).increment(Math.min(depth, CF1R2Args.maxDepth));
                return;
            }

            // If we got here, we have an alt site with a single alt base
            final Nucleotide altBase = Nucleotide.decode(BaseUtils.baseIndexToSimpleBase(altBaseIndex));

            final int refCount = baseCounts[refBase.ordinal()];
            final int altCount = baseCounts[altBaseIndex];
            Utils.validate(altCount > 0, "We must have a nonzero alt read but got " + altCount);

            final int refF1R2 = samplePileup.getNumberOfElements(pe -> Nucleotide.decode(pe.getBase()) == refBase && ReadUtils.isF1R2(pe.getRead()));
            final int altF1R2 = samplePileup.getNumberOfElements(pe -> Nucleotide.decode(pe.getBase()) == altBase && ReadUtils.isF1R2(pe.getRead()));

            if (altCount == 1) {
                final ReadOrientation type = altF1R2 == 1 ? ReadOrientation.F1R2 : ReadOrientation.F2R1;
                depthOneAltHistograms.get(sample).increment(refContext, altBase, type, depth);
                return;
            }

            try {
                altTableWriters.get(sample).get(refContext).writeRecord(new AltSiteRecord(refContext, refCount, altCount, refF1R2, altF1R2, altBase));
            } catch (IOException e) {
                throw new UserException("Encountered an IO Exception writing to the alt data table", e);
            }
        }
    }

    public void writeHistograms() {
        for (final String sample : samples) {
            final MetricsFile<?, Integer> refMetricsFile = new MetricsFile<>();
            refMetricsFile.addHeader(new StringHeader(sample));
            refSiteHistograms.get(sample).values().forEach(refMetricsFile::addHistogram);
            refMetricsFile.write(new File(tmpDir,IOUtils.urlEncode(sample) + REF_HIST_EXTENSION));

            final MetricsFile<?, Integer> altMetricsFile = new MetricsFile<>();
            altMetricsFile.addHeader(new StringHeader(sample));
            depthOneAltHistograms.get(sample).getHistograms().forEach(altMetricsFile::addHistogram);
            altMetricsFile.write(new File(tmpDir, IOUtils.urlEncode(sample) + ALT_HIST_EXTENSION));

        }
    }

    public void closeAndArchiveFiles() {
        if (altTableWriters != null) {
            altTableWriters.values().stream().flatMap(map -> map.values().stream()).forEach( writer -> {
                if (writer != null) {
                    try {
                        writer.close();
                    } catch (IOException e) {
                        throw new UserException("Encountered an IO exception while closing the alt table writer", e);
                    }
                }
            });
        }

        try {
            IOUtils.writeTarGz(outputTarGzFile.getAbsolutePath(), tmpDir.listFiles());
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile("Could not create tar.gz file", ex);
        }
    }

    /**
     * Use a series of heuristics to detect a bad pileup.
     */
    private boolean isPileupGood(final ReadPileup pileup){
        final int[] baseCounts = pileup.getBaseCounts();
        final int depth = (int) MathUtils.sum(baseCounts);

        List<Integer> mappingQualities = Ints.asList(pileup.getMappingQuals());

        // If more than 1% of the reads is indel then consider this site an indel
        final int indelThreshold = depth/100;
        boolean isIndel = pileup.getNumberOfElements(pe -> pe.isDeletion() || pe.isAfterInsertion() || pe.isBeforeDeletionStart()) > indelThreshold;

        // If depth (the sum of base counts) is 0 but the pileup is non-empty, that means all the reads
        // have deleted bases at this particular locus
        isIndel = isIndel || depth == 0 && pileup.size() > 0;

        return depth > 0 && ! isIndel && MathUtils.median(mappingQualities) >= CF1R2Args.minMedianMapQual;
    }

    public static List<File> getRefHistogramsFromExtractedTar(final File extractedTarDir) {
        return Arrays.stream(extractedTarDir.listFiles()).filter(file -> file.getAbsolutePath().endsWith(REF_HIST_EXTENSION)).collect(Collectors.toList());
    }

    public static List<File> getAltHistogramsFromExtractedTar(final File extractedTarDir) {
        return Arrays.stream(extractedTarDir.listFiles()).filter(file -> file.getAbsolutePath().endsWith(ALT_HIST_EXTENSION)).collect(Collectors.toList());
    }

    public static List<File> getAltTablesFromExtractedTar(final File extractedTarDir) {
        return Arrays.stream(extractedTarDir.listFiles()).filter(file -> file.getAbsolutePath().endsWith(ALT_TABLE_EXTENSION)).collect(Collectors.toList());
    }

    public static List<File> getAltTablesFromExtractedTar(final File extractedTarDir, final String refContext) {
        return Arrays.stream(extractedTarDir.listFiles()).filter(file -> file.getAbsolutePath().endsWith(refContext + ALT_TABLE_EXTENSION)).collect(Collectors.toList());
    }
}
