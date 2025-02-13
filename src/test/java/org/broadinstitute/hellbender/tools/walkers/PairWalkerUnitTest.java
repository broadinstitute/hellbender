package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.PrintDistantMates;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

public class PairWalkerUnitTest {
    final static SAMFileHeader hdr =
            ArtificialReadUtils.createArtificialSamHeader(3, 1, 5000);
    final static SAMSequenceDictionary dict = hdr.getSequenceDictionary();
    final static List<SimpleInterval> intervalList = Arrays.asList(
            new SimpleInterval(dict.getSequence(0).getSequenceName(), 1001, 2000),
            new SimpleInterval(dict.getSequence(1).getSequenceName(), 1001, 2000),
            new SimpleInterval(dict.getSequence(1).getSequenceName(), 3001, 4000),
            new SimpleInterval(dict.getSequence(2).getSequenceName(), 1001, 2000)
    );

    private static GATKRead createRead( final int refIndex, final int alignStart ) {
        return ArtificialReadUtils.createArtificialRead(hdr, "r1", refIndex, alignStart, 100);
    }

    @Test
    public void testRegionCheckerBasics() {
        final PairWalker.RegionChecker regionChecker = new PairWalker.RegionChecker(intervalList, dict);
        Assert.assertFalse(regionChecker.isInInterval(createRead(0, 501)));
        Assert.assertTrue(regionChecker.isInInterval(createRead(0, 1501)));
        Assert.assertFalse(regionChecker.isInInterval(createRead(0, 2501)));
        Assert.assertFalse(regionChecker.isInInterval(createRead(1, 501)));
        Assert.assertTrue(regionChecker.isInInterval(createRead(1, 1501)));
        Assert.assertTrue(regionChecker.isInInterval(createRead(1, 1701)));
        Assert.assertTrue(regionChecker.isInInterval(createRead(1, 1901)));
        Assert.assertFalse(regionChecker.isInInterval(createRead(1, 2501)));
        Assert.assertFalse(regionChecker.isInInterval(createRead(1, 2901)));
        Assert.assertTrue(regionChecker.isInInterval(createRead(1, 3501)));
        Assert.assertFalse(regionChecker.isInInterval(createRead(1, 4501)));
        Assert.assertFalse(regionChecker.isInInterval(createRead(2, 501)));
        Assert.assertTrue(regionChecker.isInInterval(createRead(2, 1501)));
        Assert.assertFalse(regionChecker.isInInterval(createRead(2, 2501)));
    }

    @Test
    public void testRegionCheckerEdges() {
        final PairWalker.RegionChecker regionChecker = new PairWalker.RegionChecker(intervalList, dict);
        Assert.assertFalse(regionChecker.isInInterval(createRead(0, 901)));
        Assert.assertTrue(regionChecker.isInInterval(createRead(0, 902)));
        Assert.assertTrue(regionChecker.isInInterval(createRead(0, 2000)));
        Assert.assertFalse(regionChecker.isInInterval(createRead(0, 2001)));
    }

    @Test
    public void testRegionCheckerIntervalSkipping() {
        final PairWalker.RegionChecker regionChecker = new PairWalker.RegionChecker(intervalList, dict);
        Assert.assertTrue(regionChecker.isInInterval(createRead(1, 1501)));
        Assert.assertTrue(regionChecker.isInInterval(createRead(2, 1501)));
    }

    private final static class PairWalkerTester extends PairWalker {
        int nPairs = 0;
        public void apply( final GATKRead read, final GATKRead mate ) { nPairs += 1; }
        public void applyUnpaired( final GATKRead read ) {
            throw new GATKException("not expecting unpaired reads during testing");
        }
    }

    private static int countValidPairs( final List<GATKRead> reads ) {
        final PairWalkerTester tester = new PairWalkerTester();
        tester.transformTraversalIntervals(intervalList, dict);
        for ( final GATKRead read : reads ) {
            tester.apply(read, null, null);
        }
        return tester.nPairs;
    }

    private static List<GATKRead> createReadPair( final int refIndex,
                                                  final int leftStart,
                                                  final int rightStart ) {
        return ArtificialReadUtils.createPair(hdr, "r1", 100, refIndex,
                leftStart, rightStart, true, false);
    }

    @Test
    public void testPairWalkerBasics() {
        Assert.assertEquals(countValidPairs(createReadPair(0, 501, 1501)), 1);
        Assert.assertEquals(countValidPairs(createReadPair(0, 1501, 1701)), 1);
        Assert.assertEquals(countValidPairs(createReadPair(0, 1701, 2501)), 1);
        Assert.assertEquals(countValidPairs(createReadPair(0, 2501, 3501)), 0);
        Assert.assertEquals(countValidPairs(createReadPair(0, 501, 3501)), 0);
    }

    @Test
    public void testUnmapped() {
        final List<GATKRead> pair = createReadPair(0, 1501, 3501);
        pair.get(1).setIsUnplaced();
        Assert.assertEquals(countValidPairs(pair), 1);
    }

    @Test
    public void testPairWalkerInterleavedPair() {
        final List<GATKRead> reads = new ArrayList<>(4);
        reads.addAll(createReadPair(0, 501, 1701));
        for ( final GATKRead read : reads ) {
            read.setName("r2");
        }
        reads.addAll(createReadPair(0, 1501, 2501));
        reads.sort(Comparator.comparingInt(GATKRead::getStart));
        Assert.assertEquals(countValidPairs(reads), 2);
    }

    private static List<GATKRead> createDistantReadPair( final int refIndex,
                                                         final int leftStart,
                                                         final int rightStart ) {
        final List<GATKRead> distantPair = createReadPair(refIndex, leftStart, rightStart);
        final GATKRead leftRead = distantPair.get(0);
        final GATKRead rightRead = distantPair.get(1);
        final GATKRead leftDistantMate =
                PrintDistantMates.doDistantMateAlterations(leftRead);
        final GATKRead rightDistantMate =
                PrintDistantMates.doDistantMateAlterations(rightRead);
        return Arrays.asList(leftRead, rightDistantMate, rightRead, leftDistantMate);
    }

    @Test
    public void testPairWalkerDistantPairs() {
        Assert.assertEquals(countValidPairs(createDistantReadPair(1, 1501, 3501)), 1);
        Assert.assertEquals(countValidPairs(createDistantReadPair(1, 501, 3501)), 1);
        Assert.assertEquals(countValidPairs(createDistantReadPair(1, 1501, 4501)), 1);
    }
}
