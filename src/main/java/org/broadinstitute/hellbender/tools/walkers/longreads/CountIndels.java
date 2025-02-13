package org.broadinstitute.hellbender.tools.walkers.longreads;

import com.google.common.base.Joiner;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.LongReadProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.*;
import java.util.*;

/**
 * Quickly count errors
 *
 * <h3>Input</h3>
 * <ul>
 *     <li>A single BAM file</li>
 * </ul>
 *
 * <h3>Outputs</h3>
 * <ul>
 *     <li>A collection of BAM files where care has been keep reads from the same ZMW in the same file</li>
 * </ul>
 *
 * <h3>Usage Example</h3>
 * <h4>Quickly count errors</h4>
 * <pre>
 *   gatk CountIndels \
 *     -I input.bam \
 *     -L region.bed \
 *     -O output.txt
 * </pre>
 */
@DocumentedFeature
@ExperimentalFeature
@CommandLineProgramProperties(
        summary = "Count indels",
        oneLineSummary = "Count indels",
        programGroup = LongReadProgramGroup.class
)
public final class CountIndels extends IntervalWalker {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="Path to which per-read information should be written")
    public String perReadPathName;

    @Argument(fullName = "outputIntervalInfo",
            shortName = "OI",
            doc="Path to which per-interval information should be written")
    public String perIntervalPathName;

    private static class DataEntry {
        private final String  readName;
        private final boolean isNegativeStrand;
        private final CigarOperator type;
        private final int           refPos;
        private final int           alleleLength;
        private final String allele;

        public DataEntry(String readName, boolean isNegativeStrand, CigarOperator type, int refPos, int alleleLength, String allele) {
            this.readName = readName;
            this.isNegativeStrand = isNegativeStrand;
            this.type = type;
            this.refPos = refPos;
            this.alleleLength = alleleLength;
            this.allele = allele;
        }

        public String getReadName() { return readName; }
        public boolean isNegativeStrand() { return isNegativeStrand; }
        public CigarOperator getType() { return type; }
        public int getRefPos() { return refPos; }
        public int getAlleleLength() { return alleleLength; }
        public String getAllele() { return allele; }

        @Override
        public String toString() {
            return "DataEntry{" +
                    "readName='" + readName + '\'' +
                    ", isNegativeStrand=" + isNegativeStrand +
                    ", type=" + type +
                    ", refPos=" + refPos +
                    ", alleleLength=" + alleleLength +
                    ", allele='" + allele + '\'' +
                    '}';
        }
    }

    private SimpleInterval curInterval = null;
    private Set<DataEntry> d = null;

    BufferedWriter readStatsWriter, intervalStatsWriter;

    @Override
    public void onTraversalStart() {
        try {
            readStatsWriter = new BufferedWriter(new FileWriter(perReadPathName));
            readStatsWriter.write(Joiner.on("\t").join(
                    "interval",
                    "readName",
                    "isNegativeStrand",
                    "lengthOnReference",
                    "numDeletedBases",
                    "numInsertedBases",
                    "delLengths",
                    "insLengths",
                    "insAlleles"
            ) + "\n");

            intervalStatsWriter = new BufferedWriter(new FileWriter(perIntervalPathName));
            intervalStatsWriter.write(Joiner.on("\t").join(
                    "interval",
                    "numReads",
                    "lengthOnReference",
                    "allNumDeletedBasesFwd",
                    "allNumDeletedBasesRev",
                    "allNumInsertedBasesFwd",
                    "allNumInsertedBasesRev",
                    "allDelLengthsFwd",
                    "allDelLengthsRev",
                    "allInsAllelesFwd",
                    "allInsAllelesRev",
                    "allInsAllelesFwd",
                    "allInsAllelesRev"
            ) + "\n");
        } catch (IOException e) {
            throw new GATKException(e.getMessage());
        }
    }

    @Override
    public void apply(SimpleInterval interval, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        if (curInterval == null || interval.getStart() != curInterval.getStart()) {
            if (d != null) {
                summarizeInterval(curInterval, d);
            }

            curInterval = interval;
            d = new HashSet<>();

            logger.info("{}", interval);
        }

        for (GATKRead read : readsContext) {
            int curOffset = 1;
            for (CigarElement ce : read.getCigar()) {
                int opOffset = curOffset;
                int refOffset = 0;
                if (ce.getOperator().isIndel()) {
                    opOffset = curOffset - 1;
                    if (ce.getOperator().equals(CigarOperator.D)) {
                        refOffset = 1;
                    }
                }
                int refPos = read.convertToSAMRecord(getHeaderForReads()).getReferencePositionAtReadPosition(opOffset) + refOffset;

                int alleleLength = ce.getLength();
                String allele = ce.getOperator().equals(CigarOperator.D) ? "." : read.getBasesString().substring(curOffset - 1, curOffset - 1 + alleleLength);

                if (ce.getOperator().isIndel() && alleleLength > 3 && interval.overlapsWithMargin(new SimpleInterval(interval.getContig(), refPos, refPos + alleleLength), 10)) {
                    d.add(new DataEntry(read.getName(), read.isReverseStrand(), ce.getOperator(), refPos, alleleLength, allele));
                }

                curOffset += ce.getOperator().consumesReadBases() ? ce.getLength() : 0;
            }
        }
    }

    private void summarizeInterval(SimpleInterval interval, Set<DataEntry> d) {
        Map<String, DataEntry> readNames = new HashMap<>();
        if ( !d.isEmpty() ) {
            for (DataEntry de : d) {
                readNames.put(de.getReadName(), de);
            }

            int allNumDeletedBasesFwd = 0;
            int allNumDeletedBasesRev = 0;
            int allNumInsertedBasesFwd = 0;
            int allNumInsertedBasesRev = 0;
            List<Integer> allDelLengthsFwd = new ArrayList<>();
            List<Integer> allDelLengthsRev = new ArrayList<>();
            List<Integer> allInsLengthsFwd = new ArrayList<>();
            List<Integer> allInsLengthsRev = new ArrayList<>();
            List<String> allInsAllelesFwd = new ArrayList<>();
            List<String> allInsAllelesRev = new ArrayList<>();

            for (String readName : readNames.keySet()) {
                int numDeletedBases = 0;
                int numInsertedBases = 0;
                boolean isNegativeStrand = false;
                List<Integer> delLengths = new ArrayList<>();
                List<Integer> insLengths = new ArrayList<>();
                List<String> insAlleles = new ArrayList<>();

                DataEntry de = readNames.get(readName);

                numDeletedBases += de.getType().equals(CigarOperator.D) ? de.alleleLength : 0;
                numInsertedBases += de.getType().equals(CigarOperator.I) ? de.alleleLength : 0;
                isNegativeStrand = de.isNegativeStrand();

                if (de.getType().equals(CigarOperator.D)) {
                    delLengths.add(de.getAlleleLength());
                } else {
                    insLengths.add(de.getAlleleLength());
                    insAlleles.add(de.allele);
                }

                if (isNegativeStrand) {
                    allNumDeletedBasesRev += numDeletedBases;
                    allNumInsertedBasesRev += numInsertedBases;
                    allDelLengthsRev.addAll(delLengths);
                    allInsLengthsRev.addAll(insLengths);
                    allInsAllelesRev.addAll(insAlleles);
                } else {
                    allNumDeletedBasesFwd += numDeletedBases;
                    allNumInsertedBasesFwd += numInsertedBases;
                    allDelLengthsFwd.addAll(delLengths);
                    allInsLengthsFwd.addAll(insLengths);
                    allInsAllelesFwd.addAll(insAlleles);
                }

                try {
                    readStatsWriter.write(Joiner.on("\t").join(
                            interval,
                            readName,
                            isNegativeStrand ? "-" : "+",
                            interval.getLengthOnReference(),
                            numDeletedBases,
                            numInsertedBases,
                            Joiner.on(",").join(delLengths),
                            Joiner.on(",").join(insLengths),
                            !insAlleles.isEmpty() ? Joiner.on(",").join(insAlleles) : "."
                    ) + "\n");
                } catch (IOException e) {
                    throw new GATKException(e.getMessage());
                }
            }

            try {
                intervalStatsWriter.write(Joiner.on("\t").join(
                        interval,
                        readNames.size(),
                        interval.getLengthOnReference(),
                        allNumDeletedBasesFwd,
                        allNumDeletedBasesRev,
                        allNumInsertedBasesFwd,
                        allNumInsertedBasesRev,
                        Joiner.on(",").join(allDelLengthsFwd),
                        Joiner.on(",").join(allDelLengthsRev),
                        Joiner.on(",").join(allInsAllelesFwd),
                        Joiner.on(",").join(allInsAllelesRev),
                        !allInsAllelesFwd.isEmpty() ? Joiner.on(",").join(allInsAllelesFwd) : ".",
                        !allInsAllelesRev.isEmpty() ? Joiner.on(",").join(allInsAllelesRev) : "."
                ) + "\n");
            } catch (IOException e) {
                throw new GATKException(e.getMessage());
            }
        }
    }

    @Override
    public void closeTool() {
        summarizeInterval(curInterval, d);

        try {
            readStatsWriter.close();
            intervalStatsWriter.close();
        } catch (IOException e) {
            throw new GATKException(e.getMessage());
        }
    }
}
