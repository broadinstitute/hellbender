package org.broadinstitute.hellbender.utils.dragstr;

import org.apache.commons.lang3.Range;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Partial mutable collection of Dragstr Parameters used to compose the final immutable {@link DragstrParams}.
 */
public final class DragstrParamsBuilder {

    private final int maxPeriod;
    private final int maxRepeats;

    public final double[][] gp;
    public final double[][] gcp;
    public final double[][] api;

    public DragstrParamsBuilder(final int maxPeriod, final int maxRepeats) {
        this.maxPeriod = maxPeriod;
        this.maxRepeats = maxRepeats;
        this.gp =  new double[maxPeriod][maxRepeats];
        this.gcp = new double[maxPeriod][maxRepeats];
        this.api = new double[maxPeriod][maxRepeats];
    }

    public void set(final int period, final Range<Integer> repeatRange, final double gp, final double gcp, final double api) {
        for (int i = repeatRange.getMinimum(); i <= repeatRange.getMaximum(); i++) {
            set(period, i, gp, gcp, api);
        }
    }

    public void set(final int period, final int repeats, final double gp, final double gcp, final double api) {
        final int p = period - 1;
        final int r = repeats - 1;
        Utils.validIndex(p, maxPeriod);
        Utils.validIndex(r, maxRepeats);
        this.gp[p][r] = gp;
        this.gcp[p][r] = gcp;
        this.api[p][r] = api;
    }

    public DragstrParams make(final DoubleSequence gopValues) {
        return make(gopValues.min(), gopValues.max(), gopValues.step());
    }

    public DragstrParams make(final double minGop, final double maxGop, final double stepwise) {

        final double[][] gop = new double[gp.length][gp[0].length];
        for (int i = 0; i < gop.length; i++) {
            for (int j = 0; j < gop[i].length; j++) {
                // seems strange to use as min-gop 0.0 and then make sure that is no less than minGop (usually 10).
                // That is the way it is in DRAGstr (matlab scripts) and the code is not totally equivalent [0..10] g
                gop[i][j] = Math.max(minGop, gopCalculation(gp[i][j], gcp[i][j], i + 1, 0.0, maxGop, stepwise));
            }
        }
        return DragstrParams.of(maxPeriod, maxRepeats, gop, gcp, api);
    }

    private double gopCalculation(final double gp, final double gcp, final int period, final double minGop, final double maxGop, final double stepwise) {

        double gpProb = Math.pow(10, -.1 * gp);
        double bestGop = minGop;
        final double c = Math.pow(10, -0.1 * gcp);
        double bestCost = Double.POSITIVE_INFINITY;
        for (double gop = minGop; maxGop - gop > -0.001; gop += stepwise) {
            final double g = Math.pow(10, -0.1 * gop);
            final double prGap = g * Math.pow(c, period - 1) * (1.0 - c);
            final double prNoGap = Math.pow(1 - 2 * g, period + 1);
            final double cost = Math.abs(prGap / prNoGap - gpProb);
            if (cost < bestCost) {
                bestGop = gop;
                bestCost = cost;
            }
        }
        return bestGop;
    }

    public double gp(final int period, final int repeatLength) {
        return gp[period - 1][repeatLength - 1];
    }

    public double api(final int period, final int repeatLength) {
        return api[period - 1][repeatLength - 1];
    }
}
