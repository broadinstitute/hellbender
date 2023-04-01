package org.broadinstitute.hellbender.utils.assembly;

import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.utils.assembly.LocalAssembler.*;
import org.broadinstitute.hellbender.utils.read.UnalignedRead;
import org.broadinstitute.hellbender.utils.read.ByteSequence;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;

public class LocalAssemblerUnitTest {
    private static final byte QMIN = LocalAssembler.QMIN_DEFAULT;
    private static final int MIN_THIN_OBS = LocalAssembler.MIN_THIN_OBS_DEFAULT;
    private static final int MIN_GAPFILL_COUNT = LocalAssembler.MIN_GAPFILL_COUNT_DEFAULT;
    private static final int TOO_MANY_TRAVERSALS = LocalAssembler.TOO_MANY_TRAVERSALS_DEFAULT;
    private static final int TOO_MANY_SCAFFOLDS = LocalAssembler.TOO_MANY_SCAFFOLDS_DEFAULT;
    private static final int MIN_SV_SIZE = LocalAssembler.MIN_SV_SIZE_DEFAULT;
    private static final int KMER_SET_CAPACITY = 200;

    private static final String[] SEQS_FOR_DOGBONE_GRAPH = new String[] {
            "ACGCGCCGGCGCAGGCGCAGAGACACATGCTACCGCGTCCAGGGGTGGAGGCATGGCGCAGGCGCAGAGA",
            "CCGCGCCGGCGCAGGCGCAGAGACACATGCTACCGCGTCCAGGGGTGGAGGCATGGCGCAGGCGCAGAGT"
    };

    // produces a graph that looks like this:
    //    -------------|                                                    |---------------
    //                 .----------------------------------------------------.
    //    -------------|                                                    |---------------
    // (i.e., 5 contigs, with a branch at each end of a longer contig)
    public static List<ContigImpl> makeDogbone( final KmerSet<KmerAdjacency> kmers ) {
        for ( final String seq : SEQS_FOR_DOGBONE_GRAPH ) {
            KmerAdjacency.kmerize(seq, MIN_THIN_OBS, kmers);
        }
        final List<ContigImpl> contigs = LocalAssembler.buildContigs(kmers);
        LocalAssembler.connectContigs(contigs);
        return contigs;
    }

    private static final String SEQ_FOR_LARIAT =
            "ACGCGCCGGCGCAGGCGCAGAGACACATGCTACCGCGTCCAGGGGTGGAGGCATGGCGCAGGCGCAGAGTCGCGCCGGCGCAGGCGCAGAGACACATGCTA";

    private static byte[] makeQuals( final int length ) {
        final byte[] qualBytes = new byte[length];
        Arrays.fill(qualBytes, QMIN);
        return qualBytes;
    }

    // produces a graph that cycles back on itself
    public static List<ContigImpl> makeLariat( final KmerSet<KmerAdjacency> kmers ) {
        KmerAdjacency.kmerize(SEQ_FOR_LARIAT, MIN_THIN_OBS, kmers);
        final List<ContigImpl> contigs = LocalAssembler.buildContigs(kmers);
        LocalAssembler.connectContigs(contigs);
        return contigs;
    }


    @Test
    public void testKmers() {
        Assert.assertTrue(Kmer.KSIZE < 32);
        Assert.assertTrue((Kmer.KSIZE & 1) != 0);
        final Kmer kmer = new Kmer(0x0BADF00DDEADBEEFL);
        for ( int call = 0; call != 4; ++call ) {
            final Kmer pred = new Kmer(kmer.getPredecessorVal(call));
            Assert.assertEquals(pred.getInitialCall(), call);
            Assert.assertEquals(new Kmer(pred.getSuccessorVal(kmer.getFinalCall())), kmer);
            final Kmer succ = new Kmer(kmer.getSuccessorVal(call));
            Assert.assertEquals(succ.getFinalCall(), call);
            Assert.assertEquals(new Kmer(succ.getPredecessorVal(kmer.getInitialCall())), kmer);
        }
    }

    @Test
    public void testKmerAdjacencies() {
        final KmerSet<KmerAdjacency> kmers = new KmerSet<>(KMER_SET_CAPACITY);
        // canonicality call is here -----v
        final String seq = "ACGTACGTACGTACACACGTACGTACGTACG";
        KmerAdjacency.kmerize(seq, 0, kmers);
        Assert.assertEquals(kmers.size(), 1);
        final KmerAdjacency kmer = kmers.iterator().next();
        Assert.assertTrue(kmer.isCanonical());
        Assert.assertEquals(kmer.getInitialCall(), 0);
        Assert.assertEquals(kmer.getFinalCall(), 2);
        Assert.assertEquals(kmer.getPredecessorCount(), 0);
        Assert.assertEquals(kmer.getPredecessorMask(), 0);
        Assert.assertNull(kmer.getSolePredecessor());
        Assert.assertEquals(kmer.getSuccessorCount(), 0);
        Assert.assertEquals(kmer.getSuccessorMask(), 0);
        Assert.assertNull(kmer.getSoleSuccessor());
        Assert.assertNull(kmer.getContig());
        Assert.assertEquals(kmer.getNObservations(), 0);
        Assert.assertSame(kmer.rc().rc(), kmer);
        Assert.assertSame(kmer.canonical(), kmer);
        Assert.assertEquals(kmer.toString(), seq);

        final ByteSequence calls = new ByteSequence(("A" + seq + "T").getBytes());
        final ByteSequence quals = new ByteSequence(makeQuals(calls.length()));
        KmerAdjacency.kmerize(calls, quals, QMIN, kmers);
        Assert.assertEquals(kmer.getNObservations(), 1);
        Assert.assertSame(KmerAdjacency.find(kmer.getKVal(), kmers), kmer);
        Assert.assertSame(KmerAdjacency.find(kmer.rc().getKVal(), kmers), kmer.rc());
        Assert.assertEquals(kmer.getPredecessorCount(), 1);
        Assert.assertEquals(kmer.getPredecessorMask(), 1); // predecessor is "A"
        final KmerAdjacency prev = kmer.getSolePredecessor();
        Assert.assertNotNull(prev);
        Assert.assertSame(KmerAdjacency.find(prev.getKVal(), kmers), prev);
        Assert.assertEquals(prev.getNObservations(), 1);
        Assert.assertSame(prev.getSoleSuccessor(), kmer);
        Assert.assertSame(prev.rc().getSolePredecessor(), kmer.rc());
        Assert.assertSame(KmerAdjacency.find(prev.rc().getKVal(), kmers), prev.rc());
        Assert.assertEquals(kmer.getSuccessorCount(), 1);
        Assert.assertEquals(kmer.getSuccessorMask(), 8); // successor is "T"
        final KmerAdjacency next = kmer.getSoleSuccessor();
        Assert.assertNotNull(next);
        Assert.assertSame(KmerAdjacency.find(next.getKVal(), kmers), next);
        Assert.assertEquals(next.getNObservations(), 1);
        Assert.assertSame(next.getSolePredecessor(), kmer);
        Assert.assertSame(next.rc().getSoleSuccessor(), kmer.rc());
        Assert.assertSame(KmerAdjacency.find(next.rc().getKVal(), kmers), next.rc());

        final ByteSequence calls2 = new ByteSequence(("C" + seq + "C").getBytes());
        KmerAdjacency.kmerize(calls2, quals, QMIN, kmers);
        Assert.assertEquals(kmer.getNObservations(), 2);
        Assert.assertEquals(kmer.getPredecessorCount(), 2);
        Assert.assertNull(kmer.getSolePredecessor());
        Assert.assertEquals(kmer.getSuccessorCount(), 2);
        Assert.assertNull(kmer.getSoleSuccessor());
    }

    @Test
    public void testSingleContigBuild() {
        final KmerSet<KmerAdjacency> kmers = new KmerSet<>(KMER_SET_CAPACITY);
        final String seq = "AACGTACGTACGTACACACGTACGTACGTACGT";
        final ByteSequence calls = new ByteSequence(seq.getBytes());
        final ByteSequence quals = new ByteSequence(makeQuals(calls.length()));
        KmerAdjacency.kmerize(calls, quals, QMIN, kmers);
        final List<ContigImpl> contigs = LocalAssembler.buildContigs(kmers);
        Assert.assertEquals(contigs.size(), 1);

        final Contig contigAsBuilt = contigs.get(0);
        Assert.assertTrue(contigAsBuilt.isCanonical()); // by definition, the one that got built is canonical
        Assert.assertFalse(contigAsBuilt.rc().isCanonical());
        Assert.assertSame(contigAsBuilt.canonical(), contigAsBuilt);
        Assert.assertSame(contigAsBuilt.rc().canonical(), contigAsBuilt);

        // which strand the contig is built for depends on minute details of the implementation
        // just grab the right one
        final String contigSeq = contigAsBuilt.getSequence().toString();
        final Contig contig = seq.equals(contigSeq) ? contigAsBuilt : contigAsBuilt.rc();
        Assert.assertEquals(contig.getSequence().toString(), seq);
        Assert.assertEquals(contig.getMaxObservations(), 1);
        Assert.assertEquals(contig.getFirstKmer().toString(), seq.substring(0, seq.length()-2));
        Assert.assertEquals(contig.getLastKmer().toString(), seq.substring(2));
        Assert.assertEquals(contig.size(), seq.length());

        Assert.assertSame(contig.rc().rc(), contig);

        Assert.assertEquals(SequenceUtil.reverseComplement(contig.getSequence().toString()),
                            contig.rc().getSequence().toString());
    }

    @Test
    public void testContigConnections() {
        final List<ContigImpl> contigs = makeDogbone(new KmerSet<>(KMER_SET_CAPACITY));
        Assert.assertEquals(contigs.size(), 5);
        for ( final Contig contig : contigs ) {
            final List<Contig> predecessors = contig.getPredecessors();
            final int nPreds = predecessors.size();
            for ( int idx = 0; idx != nPreds; ++idx ) {
                final Contig pred = predecessors.get(idx);
                Assert.assertTrue(pred.getSuccessors().contains(contig));
                Assert.assertTrue(pred.rc().getPredecessors().contains(contig.rc()));
                Assert.assertSame(contig.rc().getSuccessors().get(nPreds - idx - 1).rc(), pred);
            }
            final List<Contig> successors = contig.getSuccessors();
            final int nSuccs = successors.size();
            for ( int idx = 0; idx != nSuccs; ++idx ) {
                final Contig succ = successors.get(idx);
                Assert.assertTrue(succ.getPredecessors().contains(contig));
                Assert.assertTrue(succ.rc().getSuccessors().contains(contig.rc()));
                Assert.assertSame(contig.rc().getPredecessors().get(nSuccs - idx - 1).rc(), succ);
            }
        }
    }

    @Test
    public void testPaths() {
        final KmerSet<KmerAdjacency> kmers = new KmerSet<>(KMER_SET_CAPACITY);
        makeDogbone(kmers);
        final PathBuilder pathBuilder = new PathBuilder(kmers);
        final byte[] path1CallBytes = SEQS_FOR_DOGBONE_GRAPH[0].getBytes();
        final ByteSequence path1Calls = new ByteSequence(path1CallBytes);
        final Path path1 = new Path(path1Calls, pathBuilder);
        Assert.assertEquals(path1.getParts().size(), 3);
        Assert.assertTrue(pathsEquivalent(path1, path1.rc().rc()));

        // exercise our rudimentary error correction by changing one call in the middle of a contig
        final int pos = path1Calls.length() / 2;
        path1CallBytes[pos] = (byte)(path1CallBytes[pos] == 'A' ? 'C' : 'A');
        Assert.assertTrue(pathsEquivalent(path1, new Path(path1Calls, pathBuilder)));

        // path of RC sequence ought to be equivalent to RC of path
        SequenceUtil.reverseComplement(path1CallBytes);
        Assert.assertTrue(pathsEquivalent(path1.rc(), new Path(path1Calls, pathBuilder)));

        final Path path2 = new Path(new ByteSequence(SEQS_FOR_DOGBONE_GRAPH[1].getBytes()), pathBuilder);
        Assert.assertEquals(path2.getParts().size(), 3);
        Assert.assertFalse(pathsEquivalent(path1, path2));
    }

    // this looks like it should just be an equals method on Path and their parts,
    // but Paths are never tested for equality by the LocalAssembler, just by the unit tests.
    // so it's pulled out and implemented here.
    private static boolean pathsEquivalent( final Path path1, final Path path2 ) {
        final List<PathPart> parts1 = path1.getParts();
        final List<PathPart> parts2 = path2.getParts();
        if ( parts1.size() != parts2.size() ) {
            return false;
        }
        final int nParts = parts1.size();
        for ( int idx = 0; idx != nParts; ++idx ) {
            if ( !pathPartsEquivalent(parts1.get(idx), parts2.get(idx)) ) {
                return false;
            }
        }
        return true;
    }

    private static boolean pathPartsEquivalent( final PathPart part1, final PathPart part2 ) {
        final String part1Seq = part1.getSequence() == null ? null : part1.getSequence().toString();
        final String part2Seq = part2.getSequence() == null ? null : part2.getSequence().toString();
        return part1.getContig() == part2.getContig() &&
                Objects.equals(part1Seq, part2Seq) &&
                part1.getStart() == part2.getStart() &&
                part1.getStop() == part2.getStop();
    }

    // test removal of poorly attested contigs
    @Test
    public void testGraphTrimming() {
        final KmerSet<KmerAdjacency> kmers = new KmerSet<>(KMER_SET_CAPACITY);
        final List<ContigImpl> contigs = makeDogbone(kmers);
        Assert.assertEquals(contigs.size(), 5);

        // tear the ears off the dog bone -- the gap-fill kmerization left them at 0 max-observations
        LocalAssembler.removeThinContigs(contigs, MIN_THIN_OBS, kmers);
        Assert.assertEquals(contigs.size(), 1);
    }

    // test suppression of removing cut-point contigs
    @Test
    public void testCutPoints() {
        final KmerSet<KmerAdjacency> kmers = new KmerSet<>(KMER_SET_CAPACITY);
        makeDogbone(kmers);

        // adjust the number of observations so that the central contig has no observations
        // and the ears have MIN_THIN_OBS observations
        for ( final KmerAdjacency kmer : kmers ) {
            if ( kmer.getNObservations() > 0 ) {
                kmer.observe(null, null, -kmer.getNObservations());
            } else {
                kmer.observe(null, null, MIN_THIN_OBS);
            }
            kmer.clearContig();
        }

        // re-create the graph, and make sure central contig isn't removed despite having no
        // observations (because it's a cut point)
        final List<ContigImpl> contigs2 = LocalAssembler.buildContigs(kmers);
        LocalAssembler.connectContigs(contigs2);
        LocalAssembler.removeThinContigs(contigs2, MIN_THIN_OBS, kmers);
        Assert.assertEquals(contigs2.size(), 5);
    }

    // test joining adjacent contigs that have no branching
    @Test
    public void testWeldPipes() {
        final KmerSet<KmerAdjacency> kmers = new KmerSet<>(KMER_SET_CAPACITY);
        final List<ContigImpl> contigs = makeDogbone(kmers);
        for ( final Contig contig : contigs ) {
            // find the central contig
            if ( contig.getPredecessors().size() > 1 && contig.getSuccessors().size() > 1 ) {
                final KmerAdjacency prevKmer = contig.getPredecessors().get(0).getLastKmer();
                prevKmer.observe(null, null, MIN_THIN_OBS);
                final KmerAdjacency nextKmer = contig.getSuccessors().get(0).getFirstKmer();
                nextKmer.observe(null, null, MIN_THIN_OBS);
                break;
            }
        }
        for ( final KmerAdjacency kmer : kmers ) {
            kmer.clearContig();
        }
        final List<ContigImpl> contigs2 = LocalAssembler.buildContigs(kmers);
        LocalAssembler.connectContigs(contigs2);
        LocalAssembler.removeThinContigs(contigs2, MIN_THIN_OBS, kmers);
        Assert.assertEquals(contigs2.size(), 3);
        LocalAssembler.weldPipes(contigs2);
        Assert.assertEquals(contigs2.size(), 1);
        Assert.assertEquals(contigs2.get(0).getSequence().length(), SEQS_FOR_DOGBONE_GRAPH[0].length());
    }

    @Test
    public void testCycleDetection() {
        final KmerSet<KmerAdjacency> kmers = new KmerSet<>(KMER_SET_CAPACITY);
        final List<ContigImpl> contigs = makeLariat(kmers);
        Assert.assertEquals(contigs.size(), 2);
        final Contig c1 = contigs.get(0);
        final Contig c2 = contigs.get(1);
        final Contig longer = c1.getSequence().length() > c2.getSequence().length() ? c1 : c2;

        // longer contig is a self-cycle in the lariat sequence
        Assert.assertTrue(contigInList(longer, longer.getSuccessors()));
        Assert.assertTrue(contigInList(longer, longer.getPredecessors()));

        LocalAssembler.markCycles(contigs);
        Assert.assertTrue(longer.isCycleMember());
        final Contig shorter = longer == c1 ? c2 : c1;
        Assert.assertFalse(shorter.isCycleMember());
    }

    private boolean contigInList( final Contig target, final List<Contig> contigs ) {
        for ( final Contig contig : contigs ) {
            if ( contig == target ) return true;
        }
        return false;
    }

    @Test
    public void testGapFilling() {
        final ByteSequence calls = new ByteSequence(SEQ_FOR_LARIAT.getBytes());
        final byte[] qualBytes = makeQuals(calls.length());
        qualBytes[qualBytes.length / 2] = 0; // one bad quality score in the middle
        final ByteSequence quals = new ByteSequence(qualBytes);
        final UnalignedRead read = new UnalignedRead("noName", calls, quals);
        final KmerSet<KmerAdjacency> kmers = new KmerSet<>(KMER_SET_CAPACITY);
        final List<UnalignedRead> reads = new ArrayList<>(MIN_GAPFILL_COUNT);
        for ( int iii = 0; iii != MIN_GAPFILL_COUNT; ++iii ) {
            reads.add(read);
            KmerAdjacency.kmerize(calls, quals, QMIN, kmers);
        }
        final List<ContigImpl> contigs = LocalAssembler.buildContigs(kmers);
        LocalAssembler.connectContigs(contigs);
        LocalAssembler.markCycles(contigs);

        // broken lariat
        Assert.assertEquals(contigs.size(), 3);
        Assert.assertFalse(contigs.get(0).isCycleMember() || contigs.get(1).isCycleMember());

        // can we find and fill a gap?
        Assert.assertTrue(LocalAssembler.fillGaps(kmers, MIN_GAPFILL_COUNT, reads));

        final List<ContigImpl> contigs2 = LocalAssembler.buildContigs(kmers);
        LocalAssembler.connectContigs(contigs2);
        LocalAssembler.markCycles(contigs2);

        // lariat healed by gap fill
        Assert.assertEquals(contigs2.size(), 2);
        Assert.assertTrue(contigs2.get(0).isCycleMember() || contigs2.get(1).isCycleMember());
    }

    @Test
    public void testTraversalPhasing() {
        final ByteSequence calls1 = new ByteSequence(SEQS_FOR_DOGBONE_GRAPH[0].getBytes());
        final ByteSequence calls2 = new ByteSequence(SEQS_FOR_DOGBONE_GRAPH[1].getBytes());
        final ByteSequence quals = new ByteSequence(makeQuals(Math.max(calls1.length(),calls2.length())));
        final UnalignedRead read1 = new UnalignedRead("read1", calls1, quals);
        final UnalignedRead read2 = new UnalignedRead("read2", calls2, quals);
        final List<UnalignedRead> reads = new ArrayList<>(2);
        reads.add(read1);
        reads.add(read2);
        final KmerSet<KmerAdjacency> kmers = new KmerSet<>(KMER_SET_CAPACITY);
        LocalAssembler.kmerizeReads(reads, QMIN, kmers);
        final List<ContigImpl> contigs = LocalAssembler.buildContigs(kmers);
        LocalAssembler.connectContigs(contigs);
        final List<Path> readPaths = LocalAssembler.pathReads(kmers, reads);
        Assert.assertEquals(readPaths.size(), 2);
        Assert.assertEquals(readPaths.get(0).getParts().size(), 3);
        Assert.assertEquals(readPaths.get(1).getParts().size(), 3);
        Assert.assertTrue(pathPartsEquivalent(readPaths.get(0).getParts().get(1),
                                                readPaths.get(1).getParts().get(1)));
        final Map<Contig,List<TransitPairCount>> contigTransitsMap =
                LocalAssembler.collectTransitPairCounts(contigs, readPaths);
        Assert.assertEquals(contigTransitsMap.size(), 2);

        // we put transits for both strands into the map, so the two entries ought to be RC of each other
        final Iterator<Map.Entry<Contig, List<TransitPairCount>>> mapItr =
                contigTransitsMap.entrySet().iterator();
        final Map.Entry<Contig, List<TransitPairCount>> entry1 = mapItr.next();
        final Map.Entry<Contig, List<TransitPairCount>> entry2 = mapItr.next();
        Assert.assertSame(entry1.getKey().rc(), entry2.getKey());
        Assert.assertEquals(entry1.getValue().size(), 2);
        Assert.assertEquals(entry2.getValue().size(), 2);
        final List<TransitPairCount> list1 = entry1.getValue();
        final List<TransitPairCount> list2 = entry2.getValue();
        Assert.assertEquals(list1.get(0).getRC(), list2.get(0));
        Assert.assertEquals(list1.get(1).getRC(), list2.get(1));

        final Set<Traversal> allTraversals =
                    LocalAssembler.traverseAllPaths(contigs, readPaths, TOO_MANY_TRAVERSALS, contigTransitsMap);
        Assert.assertEquals(allTraversals.size(), 2);
        final Iterator<Traversal> travItr = allTraversals.iterator();
        final String trav1Seq = travItr.next().getSequence();
        final String trav2Seq = travItr.next().getSequence();
        if ( trav1Seq.equals(SEQS_FOR_DOGBONE_GRAPH[0]) ||
                SequenceUtil.reverseComplement(trav1Seq).equals(SEQS_FOR_DOGBONE_GRAPH[0]) ) {
            Assert.assertTrue(trav2Seq.equals(SEQS_FOR_DOGBONE_GRAPH[1]) ||
                    SequenceUtil.reverseComplement(trav2Seq).equals(SEQS_FOR_DOGBONE_GRAPH[1]));
        } else {
            Assert.assertTrue(trav1Seq.equals(SEQS_FOR_DOGBONE_GRAPH[1]) ||
                    SequenceUtil.reverseComplement(trav1Seq).equals(SEQS_FOR_DOGBONE_GRAPH[1]));
            Assert.assertTrue(trav2Seq.equals(SEQS_FOR_DOGBONE_GRAPH[0]) ||
                    SequenceUtil.reverseComplement(trav2Seq).equals(SEQS_FOR_DOGBONE_GRAPH[0]));
        }
    }

    @Test
    public void testScaffolds() {
        final ByteSequence calls = new ByteSequence(SEQS_FOR_DOGBONE_GRAPH[0].getBytes());
        final int length = calls.length() - 1;
        final ByteSequence quals = new ByteSequence(makeQuals(length));

        // this time no read transits the central contig
        final UnalignedRead read1 = new UnalignedRead("read1", calls.subSequence(1, length), quals);
        final UnalignedRead read2 = new UnalignedRead("read2", calls.subSequence(0, length), quals);

        final ByteSequence calls2 = new ByteSequence(SEQS_FOR_DOGBONE_GRAPH[1].getBytes());
        final UnalignedRead read3 = new UnalignedRead("read3", calls2.subSequence(1, length), quals);
        final UnalignedRead read4 = new UnalignedRead("read4", calls2.subSequence(0, length), quals);

        final List<UnalignedRead> reads = new ArrayList<>(4);
        reads.add(read1);
        reads.add(read2);
        reads.add(read3);
        reads.add(read4);

        final KmerSet<KmerAdjacency> kmers = new KmerSet<>(KMER_SET_CAPACITY);
        LocalAssembler.kmerizeReads(reads, QMIN, kmers);
        final List<ContigImpl> contigs = LocalAssembler.buildContigs(kmers);
        Assert.assertEquals(contigs.size(), 5); // same dogbone as before
        LocalAssembler.connectContigs(contigs);
        final List<Path> readPaths = LocalAssembler.pathReads(kmers, reads);
        final Map<Contig,List<TransitPairCount>> contigTransitsMap =
                LocalAssembler.collectTransitPairCounts(contigs, readPaths);
        Assert.assertEquals(contigTransitsMap.size(), 0);
        final List<Traversal> allTraversals = new ArrayList<>(
                LocalAssembler.traverseAllPaths(contigs, readPaths, TOO_MANY_TRAVERSALS, contigTransitsMap));

        // this should reconstruct the original dogbone sequences, but also another pair
        // with the phasing of the initial and final SNP swapped (since we have no transits
        // to establish phasing)
        final Collection<Traversal> scaffolds =
                LocalAssembler.createScaffolds(allTraversals, TOO_MANY_SCAFFOLDS, MIN_SV_SIZE);
        Assert.assertEquals(scaffolds.size(), 4);
        final List<String> scaffoldSeqs =
            scaffolds.stream().map(Traversal::getSequence).collect(Collectors.toList());
        Assert.assertTrue(containsSeqOrRC(SEQS_FOR_DOGBONE_GRAPH[0], scaffoldSeqs));
        Assert.assertTrue(containsSeqOrRC(SEQS_FOR_DOGBONE_GRAPH[1], scaffoldSeqs));
        final String phaseRev1 =
                SEQS_FOR_DOGBONE_GRAPH[1].charAt(0) + SEQS_FOR_DOGBONE_GRAPH[0].substring(1);
        Assert.assertTrue(containsSeqOrRC(phaseRev1, scaffoldSeqs));
        final String phaseRev2 =
                SEQS_FOR_DOGBONE_GRAPH[0].charAt(0) + SEQS_FOR_DOGBONE_GRAPH[1].substring(1);
        Assert.assertTrue(containsSeqOrRC(phaseRev2, scaffoldSeqs));
    }

    public boolean containsSeqOrRC( final String seq, final List<String> scaffoldSeqs ) {
        return scaffoldSeqs.contains(seq) ||
                scaffoldSeqs.contains(SequenceUtil.reverseComplement(seq));
    }

    @Test
    public void testCycles() {
        final String[] seqs = {
                "ACGCGCCGGCGCAGGCGCAGAGACACATGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACTACCGCGTCCAGGGGTGGAGGCATGGCGCAGGCGCAGAGA",
                "TCGCGCCGGCGCAGGCGCAGAGACACATGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        };
        final ByteSequence calls1 = new ByteSequence(seqs[0].getBytes());
        final ByteSequence quals1 = new ByteSequence(makeQuals(calls1.length()));
        final UnalignedRead read1 = new UnalignedRead("read1", calls1, quals1);

        final ByteSequence calls2 = new ByteSequence(seqs[1].getBytes());
        final ByteSequence quals2 = new ByteSequence(makeQuals(calls2.length()));
        final UnalignedRead read2 = new UnalignedRead("read2", calls2, quals2);

        final List<UnalignedRead> reads = new ArrayList<>(2);
        reads.add(read1);
        reads.add(read2);

        final KmerSet<KmerAdjacency> kmers = new KmerSet<>(KMER_SET_CAPACITY);
        LocalAssembler.kmerizeReads(reads, QMIN, kmers);
        final List<ContigImpl> contigs = LocalAssembler.buildContigs(kmers);
        Assert.assertEquals(contigs.size(), 5); // same dogbone as before
        LocalAssembler.connectContigs(contigs);
        LocalAssembler.markCycles(contigs);
        final List<Path> readPaths = LocalAssembler.pathReads(kmers, reads);
        final Map<Contig,List<TransitPairCount>> contigTransitsMap =
                LocalAssembler.collectTransitPairCounts(contigs, readPaths);
        // the polyA contig, the one just upstream, and their RCs make 4
        Assert.assertEquals(contigTransitsMap.size(), 4);
        final List<Traversal> traversals = new ArrayList<>(
                LocalAssembler.traverseAllPaths(contigs, readPaths, TOO_MANY_TRAVERSALS, contigTransitsMap));
        Assert.assertEquals(traversals.size(), 2);

        // this should reconstruct the original dogbone sequences, but also another pair
        // with the phasing of the initial and final SNP swapped (since we have no transits
        // to establish phasing)
        final Collection<Traversal> scaffolds =
                LocalAssembler.createScaffolds(traversals, TOO_MANY_SCAFFOLDS, MIN_SV_SIZE);
        Assert.assertEquals(scaffolds.size(), 2);
    }
}
