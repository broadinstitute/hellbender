package org.broadinstitute.hellbender.tools.walkers.longreads;

import com.google.common.base.Joiner;
import htsjdk.samtools.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.LongReadProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.*;

/**
 * Compute simple metrics on long read data
 * <h3>Input</h3>
 * <ul>
 *     <li>An unaligned or aligned BAM file</li>
 * </ul>
 *
 * <h3>Outputs</h3>
 * <ul>
 *     <li>A table with metrics</li>
 * </ul>
 *
 * <h3>Usage Example</h3>
 * <pre>
 *   gatk ComputeLongReadMetrics \
 *     -I my.bam \
 *     -O metrics
 * </pre>
 */
@DocumentedFeature
@ExperimentalFeature
@CommandLineProgramProperties(
        summary = "Compute simple metrics on long read data",
        oneLineSummary = "Compute simple metrics on long read data",
        programGroup = LongReadProgramGroup.class
)
public final class ComputeLongReadMetrics extends ReadWalker {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc="Write output to this file prefix")
    public String OUTPUT;

    @Argument(fullName = "roundLengthsToNearest", shortName="rn", doc="Round read lengths to nearest <value>")
    public Integer ROUND_TO_NEAREST = 100;

    private long num_reads = 0;
    private long yield = 0;
    private long yieldTimesPasses = 0;
    private final Map<Integer, Long> prlHist      = new TreeMap<>();
    private final Map<Integer, Long> prlYieldHist  = new TreeMap<>();
    private final Map<Integer, Long> rlHist       = new TreeMap<>();
    private final Map<Integer, Long> rlYieldHist = new TreeMap<>();
    private final Map<Long, Long>    zmwNumHist     = new TreeMap<>();
    private final Map<Integer, Long> npHist       = new TreeMap<>();
    private final Map<Integer, Long> rangeGapHist = new TreeMap<>();

    private final List<Integer> polymeraseReadLengths = new ArrayList<>();
    private final List<Integer> readLengths           = new ArrayList<>();

    private long lastZmw = -1;
    private int prl           = 0;
    private int lastRangeStop = 0;

    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        SAMRecord sr = read.convertToSAMRecord(this.getHeaderForReads());

        if (!sr.isSecondaryOrSupplementary() || sr.getReadUnmappedFlag()) {
            int np = sr.hasAttribute("np") ? sr.getIntegerAttribute("np") : 1;
            int readLength = sr.getReadLength();
            long zmw = lastZmw + 1;

            int range_start = 0;
            int range_stop = 0;

            if (sr.getReadName().contains("/")) {
                String[] pieces = sr.getReadName().split("/");
                if (pieces.length == 3) {
                    zmw = Long.parseLong(pieces[1]);

                    if (pieces[2].contains("_")) {
                        String[] range = pieces[2].split("_");
                        range_start = Integer.parseInt(range[0]);
                        range_stop = Integer.parseInt(range[1]);
                    }
                }
            }

            num_reads++;
            yield += readLength;
            yieldTimesPasses += (long) readLength * np;

            increment(rlHist, roundToNearest(readLength, ROUND_TO_NEAREST));
            increment(zmwNumHist, zmw);
            increment(npHist, np);

            add(rlYieldHist, roundToNearest(readLength, ROUND_TO_NEAREST), (long) readLength);

            if ( lastRangeStop != 0 && lastZmw == zmw) {
                increment(rangeGapHist, range_start - lastRangeStop);
            }

            readLengths.add(readLength);

            if ( lastZmw != -1 && lastZmw != zmw) {
                increment(prlHist, roundToNearest(prl, ROUND_TO_NEAREST));
                add(prlYieldHist, roundToNearest(prl, ROUND_TO_NEAREST), (long) prl);

                polymeraseReadLengths.add(prl);

                prl = 0;
            }

            prl = (zmw == lastZmw) ? prl + readLength : readLength;

            lastZmw = zmw;
            lastRangeStop = range_stop;
        }
    }

    private int roundToNearest(int num, int nearest) {
        return ((num + (nearest - 1)) / nearest) * nearest;
    }

    private void increment(Map<Long, Long> h, Long k) {
        if (!h.containsKey(k)) {
            h.put(k, 1L);
        } else {
            h.put(k, h.get(k) + 1L);
        }
    }

    private void increment(Map<Integer, Long> h, Integer k) {
        if (!h.containsKey(k)) {
            h.put(k, 1L);
        } else {
            h.put(k, h.get(k) + 1L);
        }
    }

    private void add(Map<Integer, Long> h, Integer k, Long v) {
        if (!h.containsKey(k)) {
            h.put(k, v);
        } else {
            h.put(k, h.get(k) + v);
        }
    }

    private int NX(List<Integer> l, int X) {
        double sum = 0.0;
        for (int i : l) { sum += i; }

        double h = (((double) X)/100.0)*sum;

        double cs = 0.0;
        int last = l.get(0);
        for (int i : l) {
            cs += i;
            if (cs > h) {
                break;
            }

            last = i;
        }

        return last;
    }

    @Override
    public void closeTool() {
        readLengths.sort(Collections.reverseOrder());
        polymeraseReadLengths.sort(Collections.reverseOrder());

        final PrintStream polymeraseReadLengthCountsOut;
        final PrintStream polymeraseReadLengthHistOut;
        final PrintStream polymeraseReadLengthYieldHistOut;
        final PrintStream polymeraseReadLengthNxOut;

        final PrintStream readLengthCountsOut;
        final PrintStream readLengthHistOut;
        final PrintStream readLengthYieldHistOut;
        final PrintStream readLengthNxOut;

        final PrintStream zmwHistOut;
        final PrintStream npHistOut;
        final PrintStream rangeGapHistOut;

        // TODO: add data for number of passes vs. read length boxplot

        try {
            polymeraseReadLengthCountsOut = new PrintStream(new FileOutputStream(OUTPUT + ".prl_counts.txt"));
            polymeraseReadLengthHistOut = new PrintStream(new FileOutputStream(OUTPUT + ".prl_hist.txt"));
            polymeraseReadLengthYieldHistOut = new PrintStream(new FileOutputStream(OUTPUT + ".prl_yield_hist.txt"));
            polymeraseReadLengthNxOut = new PrintStream(new FileOutputStream(OUTPUT + ".prl_nx.txt"));

            readLengthCountsOut = new PrintStream(new FileOutputStream(OUTPUT + ".rl_counts.txt"));
            readLengthHistOut = new PrintStream(new FileOutputStream(OUTPUT + ".rl_hist.txt"));
            readLengthYieldHistOut = new PrintStream(new FileOutputStream(OUTPUT + ".rl_yield_hist.txt"));
            readLengthNxOut = new PrintStream(new FileOutputStream(OUTPUT + ".rl_nx.txt"));

            zmwHistOut = new PrintStream(new FileOutputStream(OUTPUT + ".zmw_hist.txt"));
            npHistOut = new PrintStream(new FileOutputStream(OUTPUT + ".np_hist.txt"));
            rangeGapHistOut = new PrintStream(new FileOutputStream(OUTPUT + ".range_gap_hist.txt"));
        } catch (FileNotFoundException e) {
            throw new GATKException("Couldn't open output file.");
        }

        final Map<Long, Long> zmw_hist = new TreeMap<>();
        for (long zmw : zmwNumHist.keySet()) {
            increment(zmw_hist, zmwNumHist.get(zmw));
        }

        final Double polymerase_read_length_mean = polymeraseReadLengths.stream().mapToInt(Integer::intValue).average().getAsDouble();

        polymeraseReadLengthCountsOut.println(Joiner.on("\t").join( "num_polymerase_reads", "polymerase_read_yield_bp", "polymerase_read_yield_bp_times_passes", "polymerase_read_min_length_bp", "polymerase_read_max_length_bp", "polymerase_read_mean_length_bp", "polymerase_read_sd_length_bp" ));
        writeReadStatsToFile(polymeraseReadLengthCountsOut, polymerase_read_length_mean, polymeraseReadLengths);

        polymeraseReadLengthHistOut.println("polymerase_read_length\tlength_bp");
        for (int i : prlHist.keySet()) {
            polymeraseReadLengthHistOut.println(i + "\t" + prlHist.get(i));
        }

        polymeraseReadLengthYieldHistOut.println("polymerase_read_length\tyield_bp");
        for (int i : prlYieldHist.keySet()) {
            polymeraseReadLengthYieldHistOut.println(i + "\t" + prlYieldHist.get(i));
        }

        polymeraseReadLengthNxOut.println("NX\tvalue");
        readLengthNxOut.println("NX\tvalue");
        for (int X = 1; X < 100; X++) {
            polymeraseReadLengthNxOut.println(X + "\t" + NX(polymeraseReadLengths, X));
            readLengthNxOut.println(X + "\t" + NX(readLengths, X));
        }

        final Double read_length_mean = readLengths.stream().mapToInt(Integer::intValue).average().getAsDouble();

        readLengthCountsOut.println(Joiner.on("\t").join( "num_reads", "read_yield_bp", "read_yield_bp_times_passes", "read_min_length_bp", "read_max_length_bp", "read_mean_length_bp", "read_sd_length_bp" ));
        writeReadStatsToFile(readLengthCountsOut, read_length_mean, readLengths);

        readLengthHistOut.println("read_length_bp\tcount");
        for (int i : rlHist.keySet()) {
            readLengthHistOut.println(i + "\t" + rlHist.get(i));
        }

        readLengthYieldHistOut.println("read_length\tyield_bp");
        for (int i : rlYieldHist.keySet()) {
            readLengthYieldHistOut.println(i + "\t" + rlYieldHist.get(i));
        }

        zmwHistOut.println("zmw_count\tcount");
        for (long i : zmw_hist.keySet()) {
            zmwHistOut.println(i + "\t" + zmw_hist.get(i));
        }

        npHistOut.println("np\tcount");
        for (int i : npHist.keySet()) {
            npHistOut.println(i + "\t" + npHist.get(i));
        }

        rangeGapHistOut.println("range_gap\tcount");
        for (int i : rangeGapHist.keySet()) {
            rangeGapHistOut.println(i + "\t" + rangeGapHist.get(i));
        }

        polymeraseReadLengthCountsOut.close();
        polymeraseReadLengthHistOut.close();
        polymeraseReadLengthYieldHistOut.close();
        polymeraseReadLengthNxOut.close();
        readLengthCountsOut.close();
        readLengthHistOut.close();
        readLengthYieldHistOut.close();
        readLengthNxOut.close();
        zmwHistOut.close();
        npHistOut.close();
        rangeGapHistOut.close();
    }

    /**
     * Write the read length stats to a file
     * @param outStream The output stream to which to write.
     * @param readLengthMean The mean read length.
     * @param readLengthList The list of read lengths.
     */
    private void writeReadStatsToFile(final PrintStream outStream, final Double readLengthMean, final List<Integer> readLengthList) {
        outStream.println(Joiner.on("\t").join(
                readLengthList.size(),
                readLengthList.stream().mapToInt(Integer::intValue).sum(),
                yieldTimesPasses,
                readLengthList.stream().min(Integer::compareTo).get(),
                readLengthList.stream().max(Integer::compareTo).get(),
                readLengthMean,
                Math.sqrt(readLengthList.stream()
                        .map((integer) -> Math.pow(integer - readLengthMean, 2))
                        .mapToDouble(Double::doubleValue).sum() / (readLengthList.size() - 1))
        ));
    }
}
