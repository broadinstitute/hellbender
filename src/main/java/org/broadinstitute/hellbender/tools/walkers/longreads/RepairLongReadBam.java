package org.broadinstitute.hellbender.tools.walkers.longreads;

import htsjdk.samtools.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.util.List;

/**
 * Restore annotations from unaligned BAM files that are discarded during conversion by samtools fastq
 *
 * <h3>Input</h3>
 * <ul>
 *     <li>A unaligned BAM file</li>
 *     <li>A aligned BAM file</li>
 * </ul>
 *
 * <h3>Outputs</h3>
 * <ul>
 *     <li>An aligned BAM file with annotations restored</li>
 * </ul>
 *
 * <h3>Usage Example</h3>
 * <h4>Restore annotations from unaligned BAM files that may be discarded during conversion by samtools fastq</h4>
 * <pre>
 *   gatk RepairLongReadBam \
 *     -I unaligned.bam \
 *     -A aligned.bam \
 *     -O restored.bam
 * </pre>
 */
@DocumentedFeature
@ExperimentalFeature
@CommandLineProgramProperties(
        summary = "Restore annotations from unaligned BAM files that are discarded during conversion by samtools fastq",
        oneLineSummary = "Restore annotations from unaligned BAM files that are discarded during conversion by samtools fastq",
        programGroup = ReadDataManipulationProgramGroup.class
)
public final class RepairLongReadBam extends ReadWalker {
    @Argument(fullName = "aligned", shortName = "A", doc="aligned reads")
    public String aligned;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc="Write output to this file")
    public String output;

    @Argument(fullName = "sort", shortName = "s", doc="Sort output", optional = true)
    public boolean sort = false;

    private SAMFileWriter writer;
    private SAMRecordIterator it;
    private SAMRecord currentAlignedRead;

    @Override
    public void onTraversalStart() {
        if ( readArguments.getReadPathSpecifiers().size() != 1 ) {
            throw new UserException("Specify a single unaligned BAM file to -I and a single aligned BAM file to -A");
        }

        SamReaderFactory srf = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
        SamReader srs = srf.open(IOUtils.getPath(aligned));

        if (getHeaderForReads() != null && getHeaderForReads().getReadGroups() != null && getHeaderForReads().getReadGroups().size() != 1) {
            throw new UserException("One (and only one) read group per aligned/unaligned BAM file required");
        }

        it = srs.iterator();
        currentAlignedRead = it.hasNext() ? it.next() : null;

        writer = createWriter(output, srs.getFileHeader(), srs.getFileHeader().getSequenceDictionary(), sort);
    }

    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        SAMRecord unalignedRead = read.convertToSAMRecord(this.getHeaderForReads());

        do {
            if (currentAlignedRead != null) {
                if (!unalignedRead.getReadName().equals(currentAlignedRead.getReadName())) {
                    throw new UserException(
                            "Expected aligned read to have name '" +
                            unalignedRead.getReadName() +
                            "', but saw '" +
                            currentAlignedRead.getReadName() +
                            "'.  Are these unaligned and aligned versions of exactly the same data?"
                    );
                } else {
                    List<SAMRecord.SAMTagAndValue> attrs = unalignedRead.getAttributes();
                    attrs.addAll(currentAlignedRead.getAttributes());

                    for (SAMRecord.SAMTagAndValue tv : attrs) {
                        currentAlignedRead.setAttribute(tv.tag, tv.value);
                    }

                    writer.addAlignment(currentAlignedRead);
                }
            }

            currentAlignedRead = it.hasNext() ? it.next() : null;
        } while (currentAlignedRead != null && unalignedRead.getReadName().equals(currentAlignedRead.getReadName()));
    }

    @Override
    public void closeTool() {
        if ( writer != null ) {
            writer.close();
        }
    }

    /**
     * Creates SAMFileWriter
     * @return A SAMFileWriter.
     */
    private SAMFileWriter createWriter(String out, SAMFileHeader sfh, SAMSequenceDictionary ssd, boolean sortOutput) {
        sfh.setSortOrder(sortOutput ? SAMFileHeader.SortOrder.coordinate : SAMFileHeader.SortOrder.unsorted);
        sfh.setSequenceDictionary(ssd);

        return new SAMFileWriterFactory().makeSAMOrBAMWriter(sfh, !sortOutput, new File(out));
    }
}
