package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.aeonbits.owner.util.Collections;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.argparser.WorkflowProperties;
import org.broadinstitute.barclay.argparser.WorkflowOutput;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.AlignmentAgreesWithHeaderReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.util.ArrayList;
import java.util.List;

/**
 * Replace bases in reads with reference bases.
 *
 * Used to anonymize reads with information from the reference.
 * This tool is useful in the case where you want to use data for analysis,
 * but cannot publish the data without anonymizing the sequence information.
 *
 * Reads are processed, then emitted in a new file.  For a read to be processed
 * it must have valid start/end positions, sequence information, and consistent
 * lengths between the sequence, base qualities, and CIGAR string.  In addition
 * reads must be consistent with the read file header's sequence dictionary.
 *
 * For each aligned read, any base that does not match the reference is transformed
 * to the reference base.  For any transformed bases, the quality is set to a constant
 * value ({@link #refQual}).  For bases not transformed, the quality is preserved.
 * Reads transformed in this way also have their CIGARs rewritten to match the new
 * sequence information.
 */
@CommandLineProgramProperties(
        summary = "Replace bases in reads with reference bases.",
        oneLineSummary = "Replace bases in reads with reference bases.",
        programGroup = OtherProgramGroup.class
)
@DocumentedFeature
@ExperimentalFeature
@WorkflowProperties
public final class ReadAnonymizer extends ReadWalker {

    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="Output bam file.")
    @WorkflowOutput(optionalCompanions={StandardArgumentDefinitions.OUTPUT_INDEX_COMPANION})
    public GATKPath output;

    @Argument(
            fullName = "ref-base-quality",
            shortName = "ref-base-quality",
            doc = "Quality for bases that are set to the reference base.",
            minValue = 0,
            maxValue = 60,
            optional = true
    )
    public int refQual = 60;

    @Argument(
            fullName = "use-simple-cigar",
            shortName = "use-simple-cigar",
            doc = "If true, will produce a simplified cigar string (without `=` and `X`).",
            optional = true
    )
    public boolean useSimpleCigar = false;

    private SAMFileGATKReadWriter outputWriter;

    @Override
    public boolean requiresReference() { return true; }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        // Because of the operations we perform on the reads,
        // we must guarantee that certain information is contained in each read:
        return Collections.list(
            ReadFilterLibrary.VALID_ALIGNMENT_START,
            ReadFilterLibrary.VALID_ALIGNMENT_END,
            ReadFilterLibrary.READLENGTH_EQUALS_CIGARLENGTH,
            ReadFilterLibrary.SEQ_IS_STORED,
            ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS,
            ReadFilterLibrary.MAPPED,
            new AlignmentAgreesWithHeaderReadFilter()
        );
    }

    @Override
    public void onTraversalStart() {
        outputWriter = createSAMWriter(output, false);
    }

    @Override
    public void apply(final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext ) {

        final GATKRead sanitizedRead = anonymizeRead(read, referenceContext);

        // Write the read to the output file:
        outputWriter.addRead(sanitizedRead);
    }

    private GATKRead anonymizeRead(final GATKRead read, final ReferenceContext referenceContext) {
        final SimpleInterval readInterval  = new SimpleInterval(read.getContig(), read.getStart(), read.getEnd());
        final byte[]         readBases     = read.getBasesNoCopy();
        final byte[]         readbaseQuals = read.getBaseQualitiesNoCopy();
        final byte[]         refBases      = referenceContext.getBases(readInterval);

        final List<Byte>         newReadBases     = new ArrayList<>();
        final List<Byte>         newBaseQualities = new ArrayList<>();
        final List<CigarElement> newCigar         = new ArrayList<>();

        int readIndex = 0;
        int refIndex = 0;

        // Track our current cigar operator so we can accumulate cigars:
        CigarOperator currentNewCigarOp = null;
        int currentNewCigarOpCount = 0;

        CigarOperator iterCigarOp = null;
        int iterCigarOpCount = 0;

        for ( final CigarElement cigarElement :  read.getCigar().getCigarElements() ) {

            switch (cigarElement.getOperator()) {
                // For these elements we don't have to do anything special:
                case H:
                case N:
                case P:
                    iterCigarOp = cigarElement.getOperator();
                    iterCigarOpCount = cigarElement.getLength();
                    break;
                case S:
                case EQ:
                    for (int i = 0; i < cigarElement.getLength(); ++i) {
                        newReadBases.add(readBases[readIndex + i]);
                        newBaseQualities.add(readbaseQuals[readIndex + i]);
                    }
                    iterCigarOp = cigarElement.getOperator();
                    iterCigarOpCount = cigarElement.getLength();
                    break;
                // For the rest of the elements, we have to do something more:
                case M:
                    for (int i = 0; i < cigarElement.getLength(); ++i) {
                        newReadBases.add(refBases[refIndex + i]);

                        if (readBases[readIndex + i] == refBases[refIndex + i]) {
                            // Since the base is the same as the reference anyway, we keep the qual:
                            newBaseQualities.add(readbaseQuals[readIndex + i]);
                        }
                        else {
                            // Since we replaced the reference base we set the quality to our default:
                            newBaseQualities.add((byte)refQual);
                        }
                    }
                    iterCigarOp = useSimpleCigar ? CigarOperator.M : CigarOperator.EQ;
                    iterCigarOpCount = cigarElement.getLength();
                    break;
                case X:
                case D:
                    // For these operators we need to add in the reference bases as matches:
                    for (int i = 0; i < cigarElement.getLength(); ++i) {
                        newReadBases.add(refBases[refIndex + i]);
                        // Since we know it's the reference base, we set it to the max quality (60):
                        newBaseQualities.add((byte)refQual);
                    }
                    iterCigarOp = useSimpleCigar ? CigarOperator.M : CigarOperator.EQ;
                    iterCigarOpCount = cigarElement.getLength();
                    break;
                case I:
                    // Inserted bases are simply removed and ignored.
                    iterCigarOp = currentNewCigarOp;
                    iterCigarOpCount = 0;
                    break;
                default:
                    throw new UserException.MalformedFile("Unexpected cigar operation: " + cigarElement.toString());
            }

            // Update Cigar:
            if ( iterCigarOp == currentNewCigarOp ) {
                currentNewCigarOpCount += iterCigarOpCount;
            }
            else {
                if ( currentNewCigarOp != null ) {
                    newCigar.add(new CigarElement(currentNewCigarOpCount, currentNewCigarOp));
                }
                currentNewCigarOp = iterCigarOp;
                currentNewCigarOpCount = iterCigarOpCount;
            }

            // Update indices:
            if (cigarElement.getOperator().consumesReferenceBases()) {
                refIndex += cigarElement.getLength();
            }
            if (cigarElement.getOperator().consumesReadBases()) {
                readIndex += cigarElement.getLength();
            }
        }

        // Add in the last cigar element now that we're done iterating:
        newCigar.add(new CigarElement(currentNewCigarOpCount, currentNewCigarOp));

        // Replace the old read data with the new data:
        read.setCigar(new Cigar(newCigar));

        final Byte[] newBases = newReadBases.toArray(new Byte[newReadBases.size()]);
        read.setBases(ArrayUtils.toPrimitive(newBases));

        final Byte[] newQualitites = newBaseQualities.toArray(new Byte[newBaseQualities.size()]);
        read.setBaseQualities(ArrayUtils.toPrimitive(newQualitites));

        // Clean the attributes except the read group:
        final String readGroup = read.getReadGroup();
        read.clearAttributes();
        read.setReadGroup(readGroup);

        return read;
    }

    @Override
    public void closeTool() {
        if ( outputWriter != null ) {
            outputWriter.close();
        }
    }
}
