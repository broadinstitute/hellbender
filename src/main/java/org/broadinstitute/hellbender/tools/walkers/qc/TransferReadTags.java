package org.broadinstitute.hellbender.tools.walkers.qc;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.PeekableIterator;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.QCProgramGroup;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.ReadsPathDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadQueryNameComparator;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * This tool is designed for a pair of SAM files sharing the same read names (e.g.
 * SAM files before and after alignment), where one SAM lacks some read tags that are
 * present in the other.
 *
 * This situation may happen if for instance an unaligned bam is converted to fastq files,
 * and some read tags get lost during this conversion. The converting to fastq is sometimes unavoidable,
 * due to the fact that some tools, such as the adapter clipping tools, are written for fastq only,
 * while others are written specifically for SAM files.
 *
 * This tools behaves similarly to Picard MergeBamAlignment (MBA). The difference is that whereas
 * MBA merges the alignment information to the unaligned bam, TransferReadTags uses the aligned bam as the
 * base and adds the read tags from the unaligned bam to the aligned bam.
 *
 * Currently, the tool is implemented for the specific case of transfering UMI read tags (RX)
 * from an unaligned bam.
 */
@CommandLineProgramProperties(
        summary = "Incorporate read tags in a SAM file to that of a matching SAM file",
        oneLineSummary = "Incorporate read tags in a SAM file to that of a matching SAM file",
        programGroup = QCProgramGroup.class // sato: change
)
public class TransferReadTags extends GATKTool {
    final static String UMI_TAG_NAME = "RX";

    @Argument(fullName = "unmapped-sam", doc = "query-name sorted unmapped sam file containing the read tag of interest")
    public GATKPath unmappedSamFile;
    ReadsDataSource unmappedSam;

    @Argument(shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, doc = "")
    public File outSamFile;
    SAMFileGATKReadWriter writer;

    @Argument(fullName = "read-tags", doc = "read tag names to transfer")
    public List<String> readTags = new ArrayList<>();


    PeekableIterator<GATKRead> alignedSamIterator;
    PeekableIterator<GATKRead> unmappedSamIterator;

    GATKRead currentRead1;
    GATKRead currentRead2;
    ReadQueryNameComparator queryNameComparator;

    public void onTraversalStart(){
        Utils.nonEmpty(readTags, "read tags may not be empty");
        queryNameComparator = new ReadQueryNameComparator();
        alignedSamIterator = new PeekableIterator<>(directlyAccessEngineReadsDataSource().iterator());
        SAMFileHeader.SortOrder sortOrder1 = directlyAccessEngineReadsDataSource().getHeader().getSortOrder();
        Utils.validate(sortOrder1 == SAMFileHeader.SortOrder.queryname, "aligned sam must be sorted by queryname");

        unmappedSam = new ReadsPathDataSource(unmappedSamFile.toPath());

        // We would check that the unmapped bam is query-name sorted, but
        // the SortOrder field is often not populated. Thus we simply assume that
        // the unmapped sam is query-name sorted and rely on the traversal to throw an
        // error when it is not.
        unmappedSamIterator = new PeekableIterator<>(unmappedSam.iterator());

        // Initialize Read2. Read1 will be initialized in traverse()
        if (unmappedSamIterator.hasNext()){
            currentRead2 = unmappedSamIterator.next();
        } else {
            throw new UserException("unaligned sam iterator is empty.");
        }

        writer = createSAMWriter(new GATKPath(outSamFile.getAbsolutePath()), false);
    }

    @Override
    public void traverse() {
        // will need to use this later
        while (alignedSamIterator.hasNext()){
            currentRead1 = alignedSamIterator.next();
            int diff = queryNameComparator.compareReadNames(currentRead1, currentRead2);

            if (diff == 0) {
                // The query names match.
                GATKRead udpatedRead = updateReadTags(currentRead1, currentRead2);
                writer.addRead(udpatedRead);
                progressMeter.update(currentRead1);
                continue;
            } else if (diff > 0){
                // Read1 is ahead. Play unmapped reads foward until it catches up.
                while (unmappedSamIterator.hasNext()){
                    currentRead2 = unmappedSamIterator.next();
                    diff = queryNameComparator.compareReadNames(currentRead1, currentRead2);
                    if (diff > 0){
                        continue;
                    } else if (diff == 0){
                        // caught up: star moving the aligned reads forward
                        GATKRead udpatedRead = updateReadTags(currentRead1, currentRead2);
                        writer.addRead(udpatedRead);
                        break;
                    } else {
                        throw new IllegalStateException("Aligned read is lexicographically smaller than the unmapped read: " +
                                "aligned read = " + currentRead1.getName() + ", unmapped read = " + currentRead2.getName());
                    }
                }
            } else {
                throw new IllegalStateException("Aligned read is lexicographically smaller than the unmapped read: " +
                        "aligned read = " + currentRead1.getName() + ", unmapped read = " + currentRead2.getName());
            }
        }
    }

    /**
     *
     * @param targetRead We use this read as the template
     * @param originRead The read from which we extract the requested fields
     * @return A new instance of GATKRead derived from targetRead, updated with the requested fields from originRead
     */
    private GATKRead updateReadTags(final GATKRead targetRead, final GATKRead originRead){
        final GATKRead updatedRead = targetRead.copy();
        for (String tagName : readTags){
            final String tagValue = originRead.getAttributeAsString(tagName);
            Utils.nonNull(tagValue, "The attribute is empty: read " + currentRead2.getName());
            updatedRead.setAttribute(tagName, tagValue);
        }
        return updatedRead;
    }

    @Override
    public Object onTraversalSuccess(){
        Utils.validate(!alignedSamIterator.hasNext(), "aligned sam iterator has to have iterated through");
        if (writer != null) {
            writer.close();
        }
        return "SUCCESS";
    }
}

