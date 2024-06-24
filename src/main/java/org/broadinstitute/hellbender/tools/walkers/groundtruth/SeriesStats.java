package org.broadinstitute.hellbender.tools.walkers.groundtruth;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.concurrent.atomic.AtomicInteger;

public class SeriesStats {

    private static final Logger logger = LogManager.getLogger(SeriesStats.class);

    // local state
    private double last = Double.NaN;
    private int count = 0;
    private double sum = 0;
    private double min = Double.NaN;
    private double max = Double.NaN;
    private SortedMap<Double, AtomicInteger> bins = new TreeMap<>();
    private int intCount = 0;
    private Map<Double, SeriesStats> auxBins = new LinkedHashMap<>();

    public void csvWrite(final String path) throws IOException {
        logger.info("Writing SeriesStats " + toDigest() + " into " + path);
        PrintWriter pw = new PrintWriter(path);
        pw.println("value,count");
        boolean intKeys = isIntKeys();
        for (Map.Entry<Double, AtomicInteger> entry : bins.entrySet() ) {
            if ( intKeys ) {
                pw.println(String.format("%d,%d", entry.getKey().intValue(), entry.getValue().get()));
            } else {
                pw.println(String.format("%f,%d", entry.getKey(), entry.getValue().get()));
            }
        }
        pw.close();
    }

    public void add(int v) {
        add((double)v);
        intCount++;
    }

    public void add(double v) {

        // save in simple values
        last = v;
        sum += v;
        if ( count > 0 ) {
            min = Math.min(min, v);
            max = Math.max(max, v);
        } else {
            min = max = v;
        }
        count++;

        // save in bins
        final Double key = v;
        if ( bins.containsKey(key) ) {
            bins.get(key).incrementAndGet();
        } else {
            bins.put(key, new AtomicInteger(1));
        }
    }

    public double getLast() {
        return last;
    }

    public int getCount() {
        return count;
    }

    public double getMin() {
        return (count != 0) ? min : Double.NaN;
    }

    public double getMax() {
        return (count != 0) ? max : Double.NaN;
    }

    public int getUniq() {
        return bins.size();
    }

    public double getMean() {
        return (count != 0) ? (sum / count) : Double.NaN;
    }

    public double getMedian() {
        return getPercentile(50);
    }

    public double getPercentile(double precentile) {
        if ( count == 0 ) {
            return Double.NaN;
        } else if ( count == 1 ) {
            return last;
        } else {

            int percentileIndex = (int)(count * precentile / 100);
            int index = 0;
            for (Map.Entry<Double, AtomicInteger> entry : bins.entrySet() ) {
                int binSize = entry.getValue().get();
                if ( percentileIndex >= index && (percentileIndex < (index + binSize)) ) {
                    return entry.getKey();
                }
                index += binSize;
            }

            // if here, we need the highest entry
            return bins.lastKey();
        }
    }

    public double getStd() {

        if (count == 0) {
            return Double.NaN;
        }

        // calculate mean
        double mean = getMean();

        // calculate sum of sqr(deviations)
        double vsum = 0;
        for (Map.Entry<Double, AtomicInteger> entry : bins.entrySet()) {
            int binSize = entry.getValue().get();
            vsum += (Math.pow(entry.getKey() - mean, 2) * binSize);
        }

        // calculate variance and std
        double variance = vsum / count;
        return Math.sqrt(variance);
    }

    public Map<Double, AtomicInteger> getBins() {
        return this.bins;
    }

    public Map<Double, SeriesStats> getAuxBins() {
        return this.auxBins;
    }

    public String toDigest() {
        if ( isIntKeys() ) {
            return String.format("count=%d, min=%d, max=%d, median=%d, bin.count=%d", getCount(), (int)getMin(), (int)getMax(), (int)getMedian(), getBins().size());
        } else {
            return String.format("count=%d, min=%f, max=%f, median=%f, bin.count=%d", getCount(), getMin(), getMax(), getMedian(), getBins().size());
        }
    }

    private boolean isIntKeys() {
        return (count == intCount);
    }
}
