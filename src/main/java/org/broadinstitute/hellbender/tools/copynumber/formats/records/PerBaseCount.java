package org.broadinstitute.hellbender.tools.copynumber.formats.records;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.tools.copynumber.datacollection.PerBaseCountCollector;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Per-base counts at a site specified by an interval.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 * @author Robert Klein &lt;rklein@broadinstitute.org&gt;
 */
public class PerBaseCount implements Locatable {

    /* these are mandatory */
    private final SimpleInterval interval;
    private final Map<Nucleotide, Integer> perBaseCount;
    public static final List<Nucleotide> BASES = PerBaseCountCollector.BASES;

    /**
     * Construct the per-base count object.
     * @param interval the genomic interval
     * @param perBaseCount a mapping of nucleotides to counts
     */
    public PerBaseCount(final SimpleInterval interval,
                        final Map<Nucleotide, Integer> perBaseCount) {
        Utils.nonNull(interval);
        Utils.nonNull(perBaseCount);

        Utils.validateArg(interval.getStart() == interval.getEnd(), "PerBaseCount interval must have same start and end.");

        for(Integer value : perBaseCount.values()) {
            ParamUtils.isPositiveOrZero(value, "Cannot construct PerBaseCount with negative base counts.");
        }

        boolean hasAllBases = perBaseCount.keySet().containsAll(PerBaseCountCollector.BASES);
        Utils.validateArg(hasAllBases, "Constructing PerBaseCount requires perBaseCount contain all bases in {A, C, G, T, N}");
        boolean hasNoExtraBases = PerBaseCountCollector.BASES.containsAll(perBaseCount.keySet());
        Utils.validateArg(hasNoExtraBases, "Constructing PerBaseCount requires perBaseCount contain only bases in {A, C, G, T, N}");

        this.interval = interval;
        this.perBaseCount = perBaseCount;
    }

    @Override
    public String getContig() { return interval.getContig(); }

    @Override
    public int getStart() { return interval.getStart(); }

    @Override
    public int getEnd() { return interval.getEnd(); }

    public SimpleInterval getInterval() { return interval; }

    public Map<Nucleotide, Integer> getPerBaseCount() { return perBaseCount; }

    public List<Nucleotide> getBases() { return BASES; }

    public static HashMap<Nucleotide, Integer> getPerBaseCountFromCounts(
            int aCounts, int cCounts, int gCounts, int tCounts, int nCounts) {
        final HashMap<Nucleotide, Integer> perBaseCount = new HashMap<>();
        perBaseCount.put(Nucleotide.A, ParamUtils.isPositiveOrZero(aCounts, "Counts must be positive integers or zero."));
        perBaseCount.put(Nucleotide.C, ParamUtils.isPositiveOrZero(cCounts, "Counts must be positive integers or zero."));
        perBaseCount.put(Nucleotide.G, ParamUtils.isPositiveOrZero(gCounts, "Counts must be positive integers or zero."));
        perBaseCount.put(Nucleotide.T, ParamUtils.isPositiveOrZero(tCounts, "Counts must be positive integers or zero."));
        perBaseCount.put(Nucleotide.N, ParamUtils.isPositiveOrZero(nCounts, "Counts must be positive integers or zero."));
        return(perBaseCount);
    }

    public int getTotalBaseCount() {
        return perBaseCount.values().stream().mapToInt(Integer::intValue).sum();
    }

    /**
     * Check that intervals are equal and all counts are equal. We assume BASES are the same between instances.
     */
    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof PerBaseCount)) {
            return false;
        }

        final PerBaseCount count = (PerBaseCount) o;
        boolean isEqual = interval.equals(count.interval);
        for(Nucleotide base : BASES) {
            isEqual = isEqual && perBaseCount.get(base).equals(count.perBaseCount.get(base));
        }
        return isEqual;
    }

    @Override
    public int hashCode() {
        int result = interval.hashCode();
        for(Nucleotide base : BASES) {
            result = 31 * result + perBaseCount.get(base);
        }
        return result;
    }

    @Override
    public String toString() {
        String str = "PerBaseCount{interval=" + interval;
        for(Nucleotide base : BASES) {
            str = str + ", " + base.toString() + "=" + perBaseCount.get(base);
        }
        str = str + "}";
        return str;
    }
}
