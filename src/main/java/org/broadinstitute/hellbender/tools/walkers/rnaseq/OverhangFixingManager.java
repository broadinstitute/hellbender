package org.broadinstitute.hellbender.tools.walkers.rnaseq;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.Tuple;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.read.*;

import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;

/**
 * The class manages reads and splices and tries to apply overhang clipping when appropriate. Overhangs correspond to
 * places where a read is aligned partially over Intronic splits generated by SplitNCigarReads. The class is designed
 * perform two passes identical passes, where in the first pass every place that the manager changes the relevant
 * mate strand information is recorded to be repaired with setPredictedMateInformation(). Running the activateWriting()
 * switches the tool to emit output to the writer.
 * Important note: although for efficiency the manager does try to send reads to the underlying writer in coordinate
 * sorted order, it does NOT guarantee that it will do so in every case!  So unless there's a good reason not to,
 * methods that instantiate this manager should pass in a writer that does not assume the reads are pre-sorted.
 */
public class OverhangFixingManager {

    private final Map<String, Tuple<Integer, String>> mateChangedReads;

    protected static final Logger logger = LogManager.getLogger(OverhangFixingManager.class);
    private static final boolean DEBUG = false;

    // how many reads should we store in memory before flushing the queue?
    private final int maxRecordsInMemory;

    // how many mismatches do we tolerate in the overhangs?
    private final int maxMismatchesInOverhang;

    // how many bases do we tolerate in the overhang before deciding not to clip?
    private final int maxBasesInOverhang;

    // should we not bother fixing overhangs?
    private final boolean doNotFixOverhangs;

    // should we process secondary reads at all
    private final boolean processSecondaryReads;

    // should reads be written to file output or recorded to repair mate information
    private boolean outputToFile;

    // header for the reads
    private final SAMFileHeader header;

    // where we ultimately write out our records
    private final GATKReadWriter writer;

    // fasta reference reader to check overhanging edges in the exome reference sequence
    private final ReferenceSequenceFile referenceReader;

    // the genome unclippedLoc parser
    private final GenomeLocParser genomeLocParser;

    // the read cache
    private static final int INITIAL_CAPACITY = 5000;
    private final PriorityQueue<List<SplitRead>> waitingReadGroups;
    private int waitingReads;

    // the set of current splices to use
    private final Set<Splice> splices = new TreeSet<>(new SpliceComparator());

    protected static final int MAX_SPLICES_TO_KEEP = 1000;


    /**
     *
     * @param header                   header for the reads
     * @param writer                   actual writer
     * @param genomeLocParser          the GenomeLocParser object
     * @param referenceReader          the reference reader
     * @param maxRecordsInMemory       max records to keep in memory
     * @param maxMismatchesInOverhangs max number of mismatches permitted in the overhangs before requiring clipping
     * @param maxBasesInOverhangs      max number of bases permitted in the overhangs before deciding not to clip
     * @param doNotFixOverhangs        if true, don't clip overhangs at all
     * @param processSecondaryReads    if true, allow secondary reads to be overhang clipped. If false, secondary reads
     *                                 will still be written out to the writer but will not be edited by the clipper.
     */
    public OverhangFixingManager(final SAMFileHeader header,
                                 final GATKReadWriter writer,
                                 final GenomeLocParser genomeLocParser,
                                 final ReferenceSequenceFile referenceReader,
                                 final int maxRecordsInMemory,
                                 final int maxMismatchesInOverhangs,
                                 final int maxBasesInOverhangs,
                                 final boolean doNotFixOverhangs,
                                 final boolean processSecondaryReads) {
        this.header = header;
        this.writer = writer;
        this.genomeLocParser = genomeLocParser;
        this.referenceReader = referenceReader;
        this.maxRecordsInMemory = maxRecordsInMemory;
        this.maxMismatchesInOverhang = maxMismatchesInOverhangs;
        this.maxBasesInOverhang = maxBasesInOverhangs;
        this.doNotFixOverhangs = doNotFixOverhangs;
        this.waitingReadGroups = new PriorityQueue<List<SplitRead>>(INITIAL_CAPACITY, new SplitReadComparator());
        this.outputToFile = false;
        this.mateChangedReads = new HashMap<>();
        this.processSecondaryReads = processSecondaryReads;
    }

    final int getNReadsInQueue() { return waitingReads; }

    /**
     * For testing purposes only
     *
     * @return the list of reads currently in the queue
     */
    List<List<SplitRead>> getReadsInQueueForTesting() {
        return new ArrayList<List<SplitRead>>(waitingReadGroups);
    }

    /**
     * For testing purposes only
     *
     * @return the list of splices currently in the queue
     */
    public List<Splice> getSplicesForTesting() {
        return new ArrayList<>(splices);
    }

    /**
     * Add a new observed split to the list to use
     *
     * @param contig  the contig
     * @param start   the start of the split, inclusive
     * @param end     the end of the split, inclusive
     * @return the splice created, for testing purposes
     */
    @VisibleForTesting
    public Splice addSplicePosition(final String contig, final int start, final int end) {
        if ( doNotFixOverhangs ) {
            return null;
        }

        // is this a new splice?  if not, we are done
        final Splice splice = new Splice(contig, start, end);
        if ( splices.contains(splice) ) {
            return null;
        }

        // initialize it with the reference context
        // we don't want to do this until we know for sure that it's a new splice position
        splice.initialize(referenceReader);

        // clear the set of old split positions seen if we hit a new contig
        final boolean sameContig = splices.isEmpty() || splices.iterator().next().loc.getContig().equals(contig);
        if ( !sameContig ) {
            splices.clear();
        }

        // run this position against the existing reads
        waitingReadGroups.parallelStream().forEach( readGroup -> {
            final int size = readGroup.size();
            for (int i = 0; i < size; i++ ) {
                fixSplit(readGroup.get(i), splice);
            }
        } );

        splices.add(splice);

        if ( splices.size() > MAX_SPLICES_TO_KEEP ) {
            cleanSplices();
        }
        return splice;
    }

    /**
     * Add a family of split reads to the manager
     *
     * @param readGroup the family of reads to add to the manager (families are assumed be supplementary alignments to each other)
     */
    public void addReadGroup(final List<GATKRead> readGroup) {
        Utils.nonEmpty(readGroup, "readGroup added to manager is empty, which is not allowed");

        // if the new read is on a different contig or we have too many reads, then we need to flush the queue and clear the map
        final boolean tooManyReads = getNReadsInQueue() >= maxRecordsInMemory;
        GATKRead topRead = ((getNReadsInQueue()>0)? waitingReadGroups.peek().get(0).read : null);
        GATKRead firstNewGroup = readGroup.get(0);
        final boolean encounteredNewContig = getNReadsInQueue() > 0
                && !topRead.isUnmapped()
                && !firstNewGroup.isUnmapped()
                && !topRead.getContig().equals(firstNewGroup.getContig());

        if ( tooManyReads || encounteredNewContig ) {
            if ( DEBUG ) {
                logger.warn("Flushing queue on " + (tooManyReads ? "too many reads" : ("move to new contig: " + firstNewGroup.getContig() + " from " + topRead.getContig())) + " at " + firstNewGroup.getStart());
            }

            final int targetQueueSize = encounteredNewContig ? 0 : maxRecordsInMemory / 2;
            writeReads(targetQueueSize);
        }

        List<SplitRead> newReadGroup = readGroup.stream().map(this::getSplitRead).collect(Collectors.toList());

        // Check every stored read for an overhang with the new splice
        for ( final Splice splice : splices) {
            for (int i = 0; i < newReadGroup.size(); i++) {
                fixSplit(newReadGroup.get(i), splice);
            }
        }
        // add the new reads to the queue
        waitingReadGroups.add(newReadGroup);
        waitingReads = waitingReads + newReadGroup.size();

    }

    /**
     * Clean up the list of splices by removing the lowest half of sequential splices
     */
    private void cleanSplices() {
        final int targetQueueSize = splices.size() / 2;
        final Iterator<Splice> iter = splices.iterator();
        for ( int i = 0; i < targetQueueSize; i++ ) {
            iter.next();
            iter.remove();
        }
    }

    /**
     * Try to fix the given splitRead using the given split
     *
     * @param splitRead        the splitRead to fix
     * @param splice      the split (bad region to clip out)
     */
    @VisibleForTesting
    void fixSplit(final SplitRead splitRead, final Splice splice) {
        // if the splitRead doesn't even overlap the split position then we can just exit
        if ( splitRead.unclippedLoc == null || !splice.loc.overlapsP(splitRead.unclippedLoc) ) {
            return;
        }
        // if processSecondaryReads == false, filter out clipping of secondary alignments
        if (!processSecondaryReads && splitRead.read.isSecondaryAlignment()) {
            return;
        }

        GenomeLoc readLoc = splitRead.unclippedLoc;
        GATKRead read = splitRead.read;
        // Compute the number of non-clipped bases consumed according to the cigar
        final int readBasesLength = Utils.stream(read.getCigar().getCigarElements())
                                         .filter(c -> c.getOperator().consumesReadBases() && !c.getOperator().isClipping())
                                         .mapToInt(CigarElement::getLength)
                                         .sum();
        if ( isLeftOverhang(readLoc, splice.loc) ) {
            final int overhang = splice.loc.getStop() - read.getStart() + 1;
            if ( overhangingBasesMismatch(read.getBases(), read.getStart() - readLoc.getStart(), readBasesLength, splice.reference, splice.reference.length - overhang, overhang) ) {
                final GATKRead clippedRead = ReadClipper.softClipByReadCoordinates(read, 0, splice.loc.getStop() - readLoc.getStart());
                splitRead.setRead(clippedRead);
            }
        }
        else if ( isRightOverhang(readLoc, splice.loc) ) {
            final int overhang = readLoc.getStop() - splice.loc.getStart() + 1;
            if ( overhangingBasesMismatch(read.getBases(), read.getLength() - overhang, readBasesLength, splice.reference, 0, read.getEnd() - splice.loc.getStart() + 1) ) {
                final GATKRead clippedRead = ReadClipper.softClipByReadCoordinates(read, read.getLength() - overhang, read.getLength() - 1);
                splitRead.setRead(clippedRead);
            }
        }
    }

    /**
     * Is this a proper overhang on the left side of the read?
     *
     * @param readLoc    the read's unclippedLoc
     * @param spliceLoc   the split's unclippedLoc
     * @return true if it's a left side overhang
     */
    protected static boolean isLeftOverhang(final GenomeLoc readLoc, final GenomeLoc spliceLoc) {
        return readLoc.getStart() <= spliceLoc.getStop() && readLoc.getStart() > spliceLoc.getStart() && readLoc.getStop() > spliceLoc.getStop();
    }

    /**
     * Is this a proper overhang on the right side of the read?
     *
     * @param readLoc    the read's unclippedLoc
     * @param spliceLoc   the split's unclippedLoc
     * @return true if it's a right side overhang
     */
    protected static boolean isRightOverhang(final GenomeLoc readLoc, final GenomeLoc spliceLoc) {
        return readLoc.getStop() >= spliceLoc.getStart() && readLoc.getStop() < spliceLoc.getStop() && readLoc.getStart() < spliceLoc.getStart();
    }

    /**
     * Are there too many mismatches to the reference among the overhanging bases?
     *
     * @param read                  the read bases
     * @param readStartIndex        where to start on the read
     * @param readLength            the length of the read according to the reference, used to prevent overclipping
     *                              soft-clipped output from SplitNCigarReads)
     * @param reference             the reference bases
     * @param referenceStartIndex   where to start on the reference
     * @param spanToTest            how many bases to test
     * @return true if too many overhanging bases mismatch, false otherwise
     */
    protected boolean overhangingBasesMismatch(final byte[] read,
                                               final int readStartIndex,
                                               final int readLength,
                                               final byte[] reference,
                                               final int referenceStartIndex,
                                               final int spanToTest) {
        // don't process too small a span, too large a span, or a span that is most of a read
        if ( spanToTest < 1 || spanToTest > maxBasesInOverhang || spanToTest > readLength / 2 ) {
            return false;
        }

        int numMismatchesSeen = 0;
        for ( int i = 0; i < spanToTest; i++ ) {
            if ( read[readStartIndex + i] != reference[referenceStartIndex + i] ) {
                if ( ++numMismatchesSeen > maxMismatchesInOverhang) {
                    return true;
                }
            }
        }

        // we can still mismatch overall if at least half of the bases mismatch

        return numMismatchesSeen >= ((spanToTest+1)/2);
    }

    /**
     * Close out the manager stream by clearing the read cache
     */
    public void flush() {
        writeReads(0);
    }

    /**
     * Writes read groups off the top of waitingReads until the total number of reads in the set is less than the
     * target queue size. If outputToFile == false, then the it will instead mark the first item to setMateChanged in
     * mateChagnedReads
     */
    private void writeReads(int targetQueueSize) {
        // write out all of the remaining reads
        while ( getNReadsInQueue() > targetQueueSize ) {
            List<SplitRead> waitingGroup = waitingReadGroups.poll();
            waitingReads = waitingReads - waitingGroup.size();

            // Repair the supplementary groups together and add them into the writer
            if (outputToFile) {
                SplitNCigarReads.repairSupplementaryTags(waitingGroup.stream()
                        .map( r -> r.read )
                        .collect(Collectors.toList()), header);
                for (SplitRead splitRead : waitingGroup) {
                    writer.addRead(splitRead.read);
                }

                // On the first traversal we want to store reads that would be the mate
            } else {
                // Don't mark the readgroup if it is secondary (mate information should always point to the primary alignment)
                if (!waitingGroup.get(0).read.isSecondaryAlignment() && waitingGroup.get(0).hasBeenOverhangClipped()) {
                    waitingGroup.get(0).setMateChanged();
                }
            }
        }
    }

    /**
     * Activates output writing for the Overhang Fixing Manager. This command is used to allow the manager to write
     * clipped and unclipped reads to the underlying file writer.
      */
    public void activateWriting() {
        if (outputToFile) {
            throw new GATKException("Cannot activate writing for OverhangClippingManager multiple times");
        }
        flush();
        splices.clear();
        logger.info("Overhang Fixing Manager saved "+mateChangedReads.size()+" reads in the first pass");
        outputToFile = true;
    }

    // class to represent the reads with their soft-clip-included GenomeLocs
    public final class SplitRead {

        // Relevant information to determine if the read has been clipped by the manager
        private final Cigar oldCigar;
        private final int oldStart;

        public GATKRead read;
        public GenomeLoc unclippedLoc;

        public SplitRead(final GATKRead read) {
            oldCigar = read.getCigar();
            oldStart = read.getStart();
            setRead(read);
        }

        public void setRead(final GATKRead read) {
            this.read = read;
            if ( ! read.isUnmapped() ) {
                unclippedLoc = genomeLocParser.createGenomeLoc(read.getContig(), read.getSoftStart(), read.getSoftEnd());
            }
        }

        // Returns true if either of the required mate information fields have been changed by the clipper
        public boolean hasBeenOverhangClipped() {
            return (!oldCigar.equals(read.getCigar())) || (oldStart != read.getStart());
        }

        // Adds the relevant information for repairing the mate to setMateChanged keyed on a string composed of the start position
        public void setMateChanged() {
            if (!read.isUnmapped()) {
                mateChangedReads.put( makeKey(read.getName(), !read.isFirstOfPair(), oldStart ),
                        new Tuple<>(read.getStart(), TextCigarCodec.encode(read.getCigar())));
            }
        }
    }

    // Generates the string key to be used for
    private static String makeKey(String name, boolean firstOfPair, int mateStart) {
        return name + (firstOfPair ? 1 : 0) + mateStart;
    }
    /**
     * Will edit the mate MC tag and mate start position field for given read if that read has been recorded being edited
     * by the OverhangFixingManager before. Returns true if the read was edited by the tool, false otherwise.
     *
     * @param read the read to be edited
     */
    public boolean setPredictedMateInformation(GATKRead read) {
        if (!outputToFile) {
            return false;
        }
        if (!read.isEmpty() && read.isPaired()) {
            String keystring = makeKey(read.getName(), read.isFirstOfPair(), read.getMateStart());
            if (mateChangedReads.containsKey(keystring)) {
                Tuple<Integer, String> value = mateChangedReads.get(keystring);

                // update the start position so it is accurate
                read.setMatePosition(read.getMateContig(), value.a);

                // if the MC tag is present, update it too
                if (read.hasAttribute("MC")) {read.setAttribute("MC", value.b);}
                return true;
            }
        }
        return false;
    }

    // class to represent the comparator for the split reads
    private final class SplitReadComparator implements Comparator<List<SplitRead>>, Serializable {

        private static final long serialVersionUID = 7956407034441782842L;
        private final ReadCoordinateComparator readComparator;

        public SplitReadComparator() {
            readComparator = new ReadCoordinateComparator(header);
        }

        @Override
        public int compare(final List<SplitRead> readgroup1, final List<SplitRead> readgroup2) {
            return readComparator.compare(readgroup1.get(0).read, readgroup2.get(0).read);
        }
    }

    @VisibleForTesting
    SplitRead getSplitRead(final GATKRead read){
        return new SplitRead(read);
    }

    // class to represent the split positions
    protected final class Splice {

        public final GenomeLoc loc;
        public byte[] reference;

        public Splice(final String contig, final int start, final int end) {
            loc = genomeLocParser.createGenomeLoc(contig, start, end);
        }

        public void initialize(final ReferenceSequenceFile referenceReader) {
            reference = referenceReader.getSubsequenceAt(loc.getContig(), loc.getStart(), loc.getStop()).getBases();
        }

        @Override
        public boolean equals(final Object other) {
            return other != null && (other instanceof Splice) && this.loc.equals(((Splice)other).loc);
        }

        @Override
        public int hashCode() {
            return loc.hashCode();
        }
    }

    // class to represent the comparator for the split reads
    private final class SpliceComparator implements Comparator<Splice>, Serializable {
        private static final long serialVersionUID = -7783679773557594065L;

        @Override
        public int compare(final Splice position1, final Splice position2) {
            return position1.loc.compareTo(position2.loc);
        }
    }
}
