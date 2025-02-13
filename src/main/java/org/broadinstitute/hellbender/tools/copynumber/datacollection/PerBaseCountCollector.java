package org.broadinstitute.hellbender.tools.copynumber.datacollection;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.PerBaseCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.PerBaseCount;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

import java.util.*;

/**
 * Collects per-base counts at specified sites where bases are one of those specified in {@link PerBaseCountCollector#BASES}.
 * If two reads share the same query name and, therefore, originate from the same template of DNA, this
 * only counts the higher base quality nucleotide of the two reads. If the base qualities of the two reads are equal,
 * the nucleotide from read one will be counted.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 * @author Robert Klein &lt;rklein@broadinstitute.org&gt;
 */
public class PerBaseCountCollector {
    public static final List<Nucleotide> BASES = Collections.unmodifiableList(Arrays.asList(Nucleotide.A, Nucleotide.C, Nucleotide.G, Nucleotide.T, Nucleotide.N));

    private final SampleLocatableMetadata metadata;
    private final List<PerBaseCount> perBaseCounts = new ArrayList<>();

    public PerBaseCountCollector(final SampleLocatableMetadata metadata) {
        this.metadata = Utils.nonNull(metadata);
    }

    /**
     * Add counts to this class for a specific locus.
     *
     * @param pileup associated pileup at the locus.  Not {@code null}
     * @param locus position in genome to collect alellic counts.  Not {@code null}
     * @param minBaseQuality minimum base quality in the read for that read to count at that position.  Must be greater than or equal to 0.
     */
    public void collectAtLocus(final ReadPileup pileup, final Locatable locus, final int minBaseQuality) {
        Utils.nonNull(pileup);
        Utils.nonNull(locus);
        ParamUtils.isPositiveOrZero(minBaseQuality, "Minimum base quality must be zero or higher.");

        final HashMap<Nucleotide, Integer> perBaseCount = new HashMap<>();
        for(Nucleotide base : BASES) {
            perBaseCount.put(base, 0);
        }

        // Track the base quality of the first encountered read of a pair until the pair is found or iteration ends
        final HashMap<String, Integer> unpairedTemplateQuality = new HashMap<>();
        // Track the nucleotide of the first encountered read of a pair until the pair is found or iteration ends
        final HashMap<String, Nucleotide> unpairedTemplateBase = new HashMap<>();

        for (final PileupElement element : pileup) {
            final String queryName = element.getRead().getName();
            final int baseQuality = element.getQual();
            final Nucleotide base = Nucleotide.decode(element.getBase());
            final boolean firstOfPair = element.getRead().isFirstOfPair();

            if (element.isDeletion() || baseQuality < minBaseQuality) {
                continue;
            }

            // If there is no query name or the read is not paired
            if(queryName == null || !element.getRead().isPaired()) {
                final int count = perBaseCount.get(base) + 1;
                perBaseCount.put(base, count);
                continue;
            }

            if (unpairedTemplateQuality.containsKey(queryName)) {
                // If we've seen the other paired read from this template
                final int pairBaseQuality = unpairedTemplateQuality.get(queryName);
                final Nucleotide pairBase = unpairedTemplateBase.get(queryName);
                if (pairBaseQuality > baseQuality) {
                    // If the other read of the pair has a higher base quality we count its base
                    final int count = perBaseCount.get(pairBase) + 1;
                    perBaseCount.put(pairBase, count);
                } else if(pairBaseQuality == baseQuality) {
                    // If the paired bases have equal quality, count the base from the first of pair read
                    if(firstOfPair) {
                        final int count = perBaseCount.get(base) + 1;
                        perBaseCount.put(base, count);
                    } else {
                        final int count = perBaseCount.get(pairBase) + 1;
                        perBaseCount.put(pairBase, count);
                    }
                } else if(baseQuality > pairBaseQuality) {
                    // Otherwise we count this reads base
                    final int count = perBaseCount.get(base) + 1;
                    perBaseCount.put(base, count);
                }
                // We've found the read pair, so remove the first encountered read from the unpaired HashMaps
                unpairedTemplateBase.remove(queryName);
                unpairedTemplateQuality.remove(queryName);
            } else {
                // This is the first read of the pair, so add it to the unpaired HashMaps
                unpairedTemplateBase.put(queryName, base);
                unpairedTemplateQuality.put(queryName, baseQuality);
            }
        }

        // add all the unpaired bases to the perBaseCount
        for(final Nucleotide base : unpairedTemplateBase.values()) {
            final int count = perBaseCount.get(base) + 1;
            perBaseCount.put(base, count);
        }

        perBaseCounts.add(new PerBaseCount(
                new SimpleInterval(locus.getContig(), locus.getStart(), locus.getEnd()),
                perBaseCount));
    }

    /**
     * Get the per-base counts gathered so far.
     *
     * @return a <em>reference</em> to the PerBaseCountCollection
     */
    public PerBaseCountCollection getPerBaseCounts() {
        return new PerBaseCountCollection(metadata, perBaseCounts);
    }

    /**
     * Reminder that any additional information used through this method will not be able to enforce the minBaseQuality.
     *
     * @param perBaseCountCollector input data to combine with this.
     */
    public void collectFromCollector(final PerBaseCountCollector perBaseCountCollector) {
        if (perBaseCountCollector != null) {
            this.perBaseCounts.addAll(perBaseCountCollector.getPerBaseCounts().getRecords());
        }
    }

    /** TODO: Docs and input parameter checking
     *
     * @param perBaseCountCollector1
     * @param perBaseCountCollector2
     * @return a new per-base count collector with the combined contents of the two inputs
     */
    public static PerBaseCountCollector combine(final PerBaseCountCollector perBaseCountCollector1, final PerBaseCountCollector perBaseCountCollector2,
                                                final SampleLocatableMetadata sampleMetadata) {
        final PerBaseCountCollector result = new PerBaseCountCollector(sampleMetadata);
        result.collectFromCollector(perBaseCountCollector1);
        result.collectFromCollector(perBaseCountCollector2);
        return result;
    }
}
