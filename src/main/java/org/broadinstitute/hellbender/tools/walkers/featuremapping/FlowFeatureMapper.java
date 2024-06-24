package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import htsjdk.samtools.*;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.util.Precision;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.FlowBasedProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.FlowBasedReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.writers.GVCFWriter;
import org.broadinstitute.hellbender.utils.haplotype.FlowBasedHaplotype;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;

import java.util.*;


/**
 * Finds specific features in reads, scores the confidence of each feature relative to the
 * reference in each read and writes them into a VCF file.
 *
 * The sense of what a 'feature' is left somewhat open. In the most general sense, it is a haplotype
 * located in a specific location on the read. It is not necessarily defined as a deviation from the reference.
 *
 * A feature is indeed scored against the reference (in terms of its deviation).
 *
 * The current version implements a single type of feature: a SNP (aka SNV).
 *
 * <p>
 * At this point, this tool finds SNVs
 * </p>
 *
 * <h3> Input </h3>
 * <ul>
 *     <li> Coordinate-sorted and indexed SAM/BAM/CRAM </li>
 * </ul>
 *
 * <h3> Output </h3>
 * <ul>
 *     <li> Coordinate-sorted and indexed VCF </li>
 * </ul>
 *
 * <h3>Usage examples</h3>
 * Find SNVs in chromosome 20.
 * <pre>
 * gatk FlowFeatureMapper \
 *   -I input.bam \
 *   -L 20 \
 *   -O chr20_snv.vcf
 * </pre>
 *
 * {@GATK.walkertype ReadWalker}
 */
@CommandLineProgramProperties(
        summary = "Mapping features (flow space processing)",
        oneLineSummary = "Map/find features in BAM file, output VCF. Initially mapping SNVs",
        programGroup = FlowBasedProgramGroup.class
)


@DocumentedFeature
@ExperimentalFeature
public final class FlowFeatureMapper extends ReadWalker {

    static class CopyAttrInfo {
        public final String name;
        public final VCFHeaderLineType type;
        public final String desc;

        public CopyAttrInfo(final String spec) {
            final String[] toks = spec.split(",");
            name = toks[0];
            type = toks.length > 1 ? VCFHeaderLineType.valueOf(toks[1]) : VCFHeaderLineType.String;
            desc = toks.length > 2 ? StringUtils.join(Arrays.copyOfRange(toks, 2, toks.length), ",") : ("copy-attr: " + name);
        }
    }

    private static final Logger     logger = LogManager.getLogger(FlowFeatureMapper.class);

    private static final String     VCB_SOURCE = "fm";

    private static final String     VCF_READ_NAME = "X_RN";
    private static final String     VCF_SCORE = "X_SCORE";
    private static final String     VCF_FLAGS = "X_FLAGS";
    private static final String     VCF_MAPQ = "X_MAPQ";
    private static final String     VCF_CIGAR = "X_CIGAR";
    private static final String     VCF_READ_COUNT = "X_READ_COUNT";
    private static final String     VCF_FILTERED_COUNT = "X_FILTERED_COUNT";
    private static final String     VCF_FC1 = "X_FC1";
    private static final String     VCF_FC2 = "X_FC2";
    private static final String     VCF_LENGTH = "X_LENGTH";
    private static final String     VCF_EDIST = "X_EDIST";
    private static final String     VCF_INDEX = "X_INDEX";
    private static final String     VCF_SMQ_LEFT = "X_SMQ_LEFT";
    private static final String     VCF_SMQ_RIGHT = "X_SMQ_RIGHT";
    private static final String     VCF_SMQ_LEFT_MEAN = "X_SMQ_LEFT_MEAN";
    private static final String     VCF_SMQ_RIGHT_MEAN = "X_SMQ_RIGHT_MEAN";
    private static final String     VCF_ADJACENT_REF_DIFF = "X_ADJACENT_REF_DIFF";


    private static final Double     LOWEST_PROB = 0.0001;
    private static final int        VENDOR_QUALITY_CHECK_FLAG = 0x200;

    private static final String     INCLUDE_QC_FAILED_READS_FULL_NAME = "include-qc-failed-reads";

    final private List<CopyAttrInfo> copyAttrInfo = new LinkedList<>();

    // order here is according to SequenceUtil.VALID_BASES_UPPER
    final private String scoreForBaseNames[] = new String[SequenceUtil.VALID_BASES_UPPER.length];

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "File to which variants should be written")
    public GATKPath outputVCF = null;

    @ArgumentCollection
    private FlowFeatureMapperArgumentCollection fmArgs = new FlowFeatureMapperArgumentCollection();

    @Advanced
    @Argument(fullName= AssemblyBasedCallerArgumentCollection.EMIT_REF_CONFIDENCE_LONG_NAME, shortName= AssemblyBasedCallerArgumentCollection.EMIT_REF_CONFIDENCE_SHORT_NAME, doc="Mode for emitting reference confidence scores (For Mutect2, this is a BETA feature)", optional = true)
    public ReferenceConfidenceMode emitReferenceConfidence = ReferenceConfidenceMode.NONE;

    @Advanced
    @Argument(fullName = HaplotypeCallerArgumentCollection.GQ_BAND_LONG_NAME, shortName = HaplotypeCallerArgumentCollection.GQ_BAND_SHORT_NAME, doc= "Exclusive upper bounds for reference confidence GQ bands " +
            "(must be in [1, 100] and specified in increasing order)", optional = true)
    public List<Integer> GVCFGQBands = new ArrayList<>(70);
    {
        for (int i=1; i<=60; ++i) {
            GVCFGQBands.add(i);
        }
        GVCFGQBands.add(70); GVCFGQBands.add(80); GVCFGQBands.add(90); GVCFGQBands.add(99);
    };

    @Advanced
    @Argument(fullName=HaplotypeCallerArgumentCollection.OUTPUT_BLOCK_LOWER_BOUNDS, doc = "Output the band lower bound for each GQ block regardless of the data it represents", optional = true)
    public boolean floorBlocks = false;

    @Advanced
    @Argument(fullName=INCLUDE_QC_FAILED_READS_FULL_NAME, doc = "include reads with QC failed flag", optional = true)
    public boolean includeQcFailedReads = false;

    @ArgumentCollection
    public FlowBasedArgumentCollection fbargs = new FlowBasedArgumentCollection();

    protected static class ReadContext implements Comparable<ReadContext> {
        final GATKRead         read;
        final ReferenceContext referenceContext;

        ReadContext(final GATKRead read, final ReferenceContext referenceContext) {
            this.read = read;
            this.referenceContext = referenceContext;
        }

        @Override
        public int compareTo(ReadContext o) {
            int     delta = read.getContig().compareTo(o.read.getContig());

            delta = (delta != 0) ? delta : Integer.compare(read.getStart(), o.read.getStart());
            delta = (delta != 0) ? delta : Integer.compare(read.getEnd(), o.read.getEnd());

            return delta;
        }
    }

    protected static class MappedFeature implements Comparable<MappedFeature> {

        GATKRead    read;
        FlowFeatureMapperArgumentCollection.MappingFeatureEnum  type;
        byte[]      readBases;
        byte[]      refBases;
        int         readBasesOffset; // offset of read bases array
        int         start;      // location (on reference)
        int         offsetDelta;
        double      score;
        int         readCount;
        int         filteredCount;
        int         nonIdentMBasesOnRead;
        int         featuresOnRead;
        int         refEditDistance;
        int         index;
        int         smqLeft;
        int         smqRight;
        int         smqLeftMean;
        int         smqRightMean;

        double      scoreForBase[];
        boolean     adjacentRefDiff;

        public MappedFeature(GATKRead read, FlowFeatureMapperArgumentCollection.MappingFeatureEnum  type, byte[] readBases,
                             byte[] refBases, int readBasesOffset, int start, int offsetDelta) {
            this.read = read;
            this.type = type;
            this.readBases = readBases;
            this.refBases = refBases;
            this.readBasesOffset = readBasesOffset;
            this.start = start;
            this.offsetDelta = offsetDelta;
        }

        static MappedFeature makeSNV(GATKRead read, int offset, byte refBase, int start, int offsetDelta) {
            byte[]      readBases = {read.getBasesNoCopy()[offset]};
            byte[]      refBases = {refBase};
            return new MappedFeature(
                    read,
                    FlowFeatureMapperArgumentCollection.MappingFeatureEnum.SNV,
                    readBases,
                    refBases,
                    offset,
                    start,
                    offsetDelta);
        }

        @Override
        public String toString() {
            return "Feature{" +
                    "read=" + read +
                    ", type=" + type +
                    ", readBases=" + Arrays.toString(readBases) +
                    ", refBases=" + Arrays.toString(refBases) +
                    ", readBasesOffset=" + readBasesOffset +
                    ", start=" + start +
                    '}';
        }

        @Override
        public int compareTo(MappedFeature o) {

            int     delta = this.read.getContig().compareTo(o.read.getContig());
            if ( delta != 0 ) {
                return delta;
            } else {
                return this.start - o.start;
            }
        }
    }

    // locals
    private VariantContextWriter                vcfWriter;
    final private PriorityQueue<MappedFeature>        featureQueue = new PriorityQueue<>();
    final private PriorityQueue<ReadContext>    readQueue = new PriorityQueue<>();
    private FeatureMapper                       mapper;

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        mapper = buildMapper();

        // enforce requirement for sorted input
        if ( getHeaderForReads().getSortOrder() != SAMFileHeader.SortOrder.coordinate ) {
            throw new IllegalArgumentException("input file must be coordinated sorted");
        }

        // open output vcf
        // The HC engine will make the right kind (VCF or GVCF) of writer for us
        final SAMSequenceDictionary sequenceDictionary = getHeaderForReads().getSequenceDictionary();
        vcfWriter = makeVCFWriter(outputVCF, sequenceDictionary, createOutputVariantIndex, createOutputVariantMD5, outputSitesOnlyVCFs);
        vcfWriter.writeHeader(makeVCFHeader(sequenceDictionary, getDefaultToolVCFHeaderLines()));
    }

    @Override
    public void closeTool() {
        flushQueue(null, null);
        super.closeTool();
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }

    public VariantContextWriter makeVCFWriter( final GATKPath outputVCF, final SAMSequenceDictionary readsDictionary,
                                               final boolean createOutputVariantIndex, final boolean  createOutputVariantMD5,
                                               final boolean sitesOnlyMode ) {
        Utils.nonNull(outputVCF);
        Utils.nonNull(readsDictionary);

        final List<Options> options = new ArrayList<>(2);
        if (createOutputVariantIndex) {options.add(Options.INDEX_ON_THE_FLY);}
        if (sitesOnlyMode) {options.add(Options.DO_NOT_WRITE_GENOTYPES);}

        VariantContextWriter writer = GATKVariantContextUtils.createVCFWriter(
                outputVCF.toPath(),
                readsDictionary,
                createOutputVariantMD5,
                options.toArray(new Options[options.size()])
        );

        if ( emitReferenceConfidence == ReferenceConfidenceMode.GVCF ) {
            try {
                writer = new GVCFWriter(writer, new ArrayList<Number>(GVCFGQBands), floorBlocks);
            } catch ( IllegalArgumentException e ) {
                throw new CommandLineException.BadArgumentValue("GQBands", "are malformed: " + e.getMessage());
            }
        }

        return writer;
    }

    public VCFHeader makeVCFHeader(final SAMSequenceDictionary sequenceDictionary, final Set<VCFHeaderLine>  defaultToolHeaderLines ) {
        final Set<VCFHeaderLine> headerInfo = new HashSet<>();
        headerInfo.addAll(defaultToolHeaderLines);

        // all callers need to add these standard annotation header lines
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));

        // all callers need to add these standard FORMAT field header lines
        VCFStandardHeaderLines.addStandardFormatLines(headerInfo, true,
                VCFConstants.GENOTYPE_KEY,
                VCFConstants.GENOTYPE_QUALITY_KEY,
                VCFConstants.DEPTH_KEY,
                VCFConstants.GENOTYPE_PL_KEY);

        // add our own headers
        headerInfo.add(new VCFInfoHeaderLine(VCF_READ_NAME, 1, VCFHeaderLineType.String, "Read name"));
        headerInfo.add(new VCFInfoHeaderLine(VCF_SCORE, 1, VCFHeaderLineType.Float, "Mapping score"));
        headerInfo.add(new VCFInfoHeaderLine(VCF_FLAGS, 1, VCFHeaderLineType.Integer, "Read flags"));
        headerInfo.add(new VCFInfoHeaderLine(VCF_MAPQ, 1, VCFHeaderLineType.Integer, "Read mapqe"));
        headerInfo.add(new VCFInfoHeaderLine(VCF_CIGAR, 1, VCFHeaderLineType.String, "Read CIGAR"));
        headerInfo.add(new VCFInfoHeaderLine(VCF_READ_COUNT, 1, VCFHeaderLineType.Integer, "Number of reads containing this location"));
        headerInfo.add(new VCFInfoHeaderLine(VCF_FILTERED_COUNT, 1, VCFHeaderLineType.Integer, "Number of reads containing this location that agree with reference according to fitler"));
        headerInfo.add(new VCFInfoHeaderLine(VCF_FC1, 1, VCFHeaderLineType.Integer, "Number of M bases different on read from references"));
        headerInfo.add(new VCFInfoHeaderLine(VCF_FC2, 1, VCFHeaderLineType.Integer, "Number of features before score threshold filter"));
        headerInfo.add(new VCFInfoHeaderLine(VCF_LENGTH, 1, VCFHeaderLineType.Integer, "Read length"));
        headerInfo.add(new VCFInfoHeaderLine(VCF_EDIST, 1, VCFHeaderLineType.Integer, "Read Levenshtein edit distance from reference"));
        headerInfo.add(new VCFInfoHeaderLine(VCF_INDEX, 1, VCFHeaderLineType.Integer, "Ordinal index, from start of the read, where the feature was found"));
        headerInfo.add(new VCFInfoHeaderLine(VCF_SMQ_LEFT, 1, VCFHeaderLineType.Integer, "Ordinal Median quality of N bases to the left of the feature"));
        headerInfo.add(new VCFInfoHeaderLine(VCF_SMQ_RIGHT, 1, VCFHeaderLineType.Integer, "Ordinal Median quality of N bases to the right of the feature"));
        headerInfo.add(new VCFInfoHeaderLine(VCF_SMQ_LEFT_MEAN, 1, VCFHeaderLineType.Integer, "Ordinal Mean quality of N bases to the left of the feature"));
        headerInfo.add(new VCFInfoHeaderLine(VCF_SMQ_RIGHT_MEAN, 1, VCFHeaderLineType.Integer, "Ordinal Mean quality of N bases to the right of the feature"));
        for ( String spec : fmArgs.copyAttr ) {
            final CopyAttrInfo info = new CopyAttrInfo(spec);
            headerInfo.add(new VCFInfoHeaderLine(fmArgs.copyAttrPrefix + info.name, 1, info.type, info.desc));
            copyAttrInfo.add(info);
        }

        // validation mode?
        if ( fmArgs.reportAllAlts ) {
            for ( int baseIndex = 0 ; baseIndex < SequenceUtil.VALID_BASES_UPPER.length ; baseIndex++ ) {
                headerInfo.add(new VCFInfoHeaderLine(scoreNameForBase(baseIndex), 1, VCFHeaderLineType.Float, "Base specific mapping score"));
            }
        }
        if ( fmArgs.reportAllAlts || fmArgs.tagBasesWithAdjacentRefDiff ) {
            headerInfo.add(new VCFInfoHeaderLine(VCF_ADJACENT_REF_DIFF, 1, VCFHeaderLineType.Flag, "Adjacent base filter indication: indel in the adjacent 5 bases to the considered base on the read"));
        }

        final VCFHeader vcfHeader = new VCFHeader(headerInfo);
        vcfHeader.setSequenceDictionary(sequenceDictionary);
        return vcfHeader;
    }

    private String scoreNameForBase(int baseIndex) {
        if ( scoreForBaseNames[baseIndex] == null ) {
            scoreForBaseNames[baseIndex] = VCF_SCORE + "_" + new String(new byte[]{SequenceUtil.VALID_BASES_UPPER[baseIndex]});
        }
        return scoreForBaseNames[baseIndex];
    }

    @Override
    public void apply(final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        // include dups?
        if ( read.isDuplicate() && !fmArgs.includeDupReads ) {
            return;
        }

        // include supplementary alignments?
        if ( read.isSupplementaryAlignment() && !fmArgs.keepSupplementaryAlignments ) {
            return;
        }

        // include qc-failed reads?
        if ( ((read.getFlags() & VENDOR_QUALITY_CHECK_FLAG) != 0) && !includeQcFailedReads ) {
            return;
        }

        // flush qeues up to this read
        flushQueue(read, referenceContext);

        // find features in read
        mapper.forEachOnRead(read, referenceContext, fr -> {
            if ( logger.isDebugEnabled() ) {
                logger.debug("fr: " + fr);
            }

            // score the feature
            fr.score = scoreFeature(fr);

            // score for validation mode?
            if ( fmArgs.reportAllAlts) {
                fr.scoreForBase = new double[SequenceUtil.VALID_BASES_UPPER.length];
                for ( int baseIndex = 0 ; baseIndex < fr.scoreForBase.length ; baseIndex++ ) {
                    final byte base = SequenceUtil.VALID_BASES_UPPER[baseIndex];
                    boolean incl = base != fr.readBases[0];
                    if ( incl ) {
                        fr.scoreForBase[baseIndex] = scoreFeature(fr, base);
                    } else {
                        fr.scoreForBase[baseIndex] = Double.NaN;
                    }
                }
            }

            // emit feature if filters in
            if ( filterFeature(fr) ) {
                featureQueue.add(fr);
            }
        });
    }

    private void flushQueue(final GATKRead read, final ReferenceContext referenceContext) {

        // emit all?
        if ( read == null ) {
            while ( featureQueue.size() != 0 ) {
                final MappedFeature fr = featureQueue.poll();
                enrichFeature(fr);
                emitFeature(fr);
            }
        } else {
            // enter read into the queue
            readQueue.add(new ReadContext(read, referenceContext));

            // emit all features that start before this read
            while ( featureQueue.size() != 0 ) {
                MappedFeature fr = featureQueue.peek();
                if ( !fr.read.getContig().equals(read.getContig())
                            || (fr.start < read.getStart()) ) {
                    fr = featureQueue.poll();
                    enrichFeature(fr);
                    emitFeature(fr);
                }
                else {
                    break;
                }
            }

            // remove all reads that start before this read
            while ( readQueue.size() != 0 ) {
                ReadContext rc = readQueue.peek();

                if ( !rc.read.getContig().equals(read.getContig())
                        || (rc.read.getEnd() < read.getStart()) ) {
                    rc = readQueue.poll();
                }
                else {
                    break;
                }
            }
        }
    }

    private void enrichFeature(final MappedFeature fr) {

        // loop on queued reads, count and check if should be counted as filtered
        final Locatable   loc = new SimpleInterval(fr.read.getContig(), fr.start, fr.start);
        for ( ReadContext rc : readQueue ) {
            if ( rc.read.contains(loc) ) {
                fr.readCount++;
                if ( mapper.noFeatureButFilterAt(rc.read, rc.referenceContext, fr.start) == FeatureMapper.FilterStatus.Filtered ) {
                    fr.filteredCount++;
                }
            }
        }
    }

    private double scoreFeature(final MappedFeature fr) {
        return scoreFeature(fr, (byte)0);
    }
    private double scoreFeature(final MappedFeature fr, byte altBase) {

        // build haplotypes
        final FlowBasedReadUtils.ReadGroupInfo rgInfo = FlowBasedReadUtils.getReadGroupInfo(getHeaderForReads(), fr.read);
        final FlowBasedHaplotype[]    haplotypes = buildHaplotypes(fr, rgInfo.flowOrder, altBase);

        // create flow read
        final FlowBasedRead flowRead = new FlowBasedRead(fr.read, rgInfo.flowOrder, rgInfo.maxClass, fbargs);

        final int diffLeft = haplotypes[0].getStart() - flowRead.getStart() + fr.offsetDelta;
        final int diffRight = flowRead.getEnd() - haplotypes[0].getEnd();
        flowRead.applyBaseClipping(Math.max(0, diffLeft), Math.max(diffRight, 0), false);

        if ( !flowRead.isValid() ) {
            return -1;
        }

        // compute alternative score
        final int         hapKeyLength = Math.min(haplotypes[0].getKeyLength(), haplotypes[1].getKeyLength());
        final double      readScore = computeLikelihoodLocal(flowRead, haplotypes[0], hapKeyLength, false);
        final double      refScore = computeLikelihoodLocal(flowRead, haplotypes[1], hapKeyLength, false);
        double            score = readScore - refScore;
        if ( !Double.isNaN(fmArgs.limitScore) ) {
            score = Math.min(score, fmArgs.limitScore);
        }

        if ( ((Double.isNaN(score) || (score < 0)) && fmArgs.debugNegatives)
             || (fmArgs.debugReadName != null && fmArgs.debugReadName.contains(fr.read.getName())) ) {
            logger.info("**** debug read: " + fr.read);
            logger.info("readBases: " + fr.read.getBasesString());
            logger.info("flowRead: " + flowRead);
            logger.info("flowBases: " + flowRead.getBasesString());
            logger.info("flowOrder: " + flowRead.getFlowOrder());
            logger.info("flowKey: " + flowRead.getKeyLength() + " " + Arrays.toString(flowRead.getKey()));
            logger.info("readHaplotype: " + haplotypes[0]);
            logger.info("readHapKey: " + haplotypes[0].getKeyLength() + " " + Arrays.toString(haplotypes[0].getKey()));
            computeLikelihoodLocal(flowRead, haplotypes[0], hapKeyLength, true);
            logger.info("refrHaplotype: " + haplotypes[1]);
            logger.info("refrHapKey: " + haplotypes[1].getKeyLength() + " " + Arrays.toString(haplotypes[1].getKey()));
            computeLikelihoodLocal(flowRead, haplotypes[1], hapKeyLength, true);
            logger.info("score: " + score);

            // analyze read
            final FlowBasedRead flowRead2 = new FlowBasedRead(fr.read, rgInfo.flowOrder, rgInfo.maxClass, fbargs);
            final int[]        key2 = flowRead2.getKey();
            for ( int i = 0 ; i < key2.length ; i++ ) {
                final double      p1 = flowRead2.getProb(i, key2[i]);
                for ( int j = 0 ; j < rgInfo.maxClass ; j++ ) {
                    final double      p2 = flowRead2.getProb(i, j);
                    if ( p2 > p1 )
                        logger.info(String.format("prob at %s key[%d]=%d, %f is lower than at %d which is %f",
                                                                flowRead2.getName(), i, key2[i], p1, j, p2));
                }
            }
        }

        if ( score < 0 && !fmArgs.keepNegatives && score != -1.0 ) {
            score = 0;
        }

        return score;
    }

    public static double computeLikelihoodLocal(final FlowBasedRead read, final FlowBasedHaplotype haplotype, final int hapKeyLength, final boolean debug) {

        final byte[] flowOrder = haplotype.getFlowOrderArray();
        final byte   readFlowOrder0 = read.getFlowOrderArray()[0];
        int startingPoint = 0;
        for (int i = 0; i < flowOrder.length; i++) {
            if (flowOrder[i] == readFlowOrder0) {
                startingPoint = i;
                break;
            }
        }
        final int[]         key = haplotype.getKey();

        // debug support
        StringBuffer        debugMessage = null;
        if ( debug )
            debugMessage = new StringBuffer(Integer.toString(startingPoint) + " hmer prob |");
        double              result = 0 ;
        for (int i = 0; i < read.getKeyLength(); i++) {
            int     index = i + startingPoint;
            double  prob = 0;
            int     locationToFetch = 0;
            if ( index < hapKeyLength ) {
                locationToFetch = Math.min(key[index] & 0xff, read.getMaxHmer() + 1);
                prob = read.getProb(i, locationToFetch);
            } else {
                if ( debug ) {
                    debugMessage.append(" clip");
                }
                break;
            }
            if ( Precision.equals(prob, 0.0) ) {
                prob = LOWEST_PROB;
            }
            result += Math.log10(prob);

            if ( debug ) {
                debugMessage.append(String.format(" %d %.4f", locationToFetch, prob));
            }
        }

        if ( debug ) {
            debugMessage.append(" | " + result);
            logger.info("debugMessage: " + debugMessage);
        }

        return result;
    }

    private FlowBasedHaplotype[] buildHaplotypes(final MappedFeature fr, final String flowOrder, byte altBase) {

        // build bases for flow haplotypes
        // NOTE!!!: this code assumes length of feature on read and reference is the same
        // this is true for SNP but not for INDELs - it will have to be re-written!
        // TODO: write for INDEL
        byte[] bases = fr.read.getBasesNoCopy();
        int         offset = fr.readBasesOffset;
        int         refStart = fr.start;
        int         refModOfs = 0;

        // install alt base?
        byte orgBase = 0;
        if ( altBase != 0 ) {
            orgBase = fr.refBases[0];
            fr.refBases[0] = altBase;
        }

        if ( offset > 0 ) {
            // reach into hmer before
            offset--;
            refModOfs++;
            refStart--;

            // extend until start of hmer
            final byte        hmerBase = bases[offset];
            while ( offset > 0 && bases[offset-1] == hmerBase ) {
                offset--;
                refModOfs++;
                refStart--;
            }
        }
        final byte[]      sAltBases = Arrays.copyOfRange(bases, offset, bases.length);
        final byte[]      sRefBases = Arrays.copyOf(sAltBases, sAltBases.length);
        System.arraycopy(fr.refBases, 0, sRefBases, refModOfs, fr.refBases.length);

        // construct haplotypes
        final SimpleInterval genomeLoc = new SimpleInterval(fr.read.getContig(), refStart, refStart + sAltBases.length - 1);
        final Cigar          cigar = new Cigar();
        cigar.add(new CigarElement(sAltBases.length, CigarOperator.M));
        final Haplotype      altHaplotype = new Haplotype(sAltBases, false);
        final Haplotype      refHaplotype = new Haplotype(sRefBases, true);
        altHaplotype.setGenomeLocation(genomeLoc);
        refHaplotype.setGenomeLocation(genomeLoc);
        altHaplotype.setCigar(cigar);
        refHaplotype.setCigar(cigar);

        // prepare flow based haplotypes
        final FlowBasedHaplotype[] result = {
                                new FlowBasedHaplotype(altHaplotype, flowOrder),
                                new FlowBasedHaplotype(refHaplotype, flowOrder)
                            };

        // restore changes
        if ( altBase != 0 ) {
            fr.refBases[0] = orgBase;
        }

        // return
        return result;
    }

    private boolean filterFeature(final MappedFeature fr) {

        if ( fmArgs.excludeNaNScores && Double.isNaN(fr.score) ) {
            return false;
        } else if ( fr.score > fmArgs.maxScore ) {
            return false;
        } else if ( fr.score < fmArgs.minScore ) {
            return false;
        }

        return true;
    }

    private void emitFeature(final MappedFeature fr) {

        // create alleles
        final Collection<Allele>          alleles = new LinkedList<>();
        if ( fmArgs.reportAllAlts && Arrays.equals(fr.readBases, fr.refBases) ) {
            alleles.add(Allele.create("*".getBytes(), false));
        } else {
            alleles.add(Allele.create(fr.readBases, false));
        }
        alleles.add(Allele.create(fr.refBases, true));

        // create variant context builder
        final VariantContextBuilder       vcb = new VariantContextBuilder(
                                                    VCB_SOURCE,
                                                    fr.read.getContig(),
                                                    fr.start,
                                                    fr.start + fr.refBases.length - 1,
                                                    alleles);

        // copy attributes
        vcb.attribute(VCF_READ_NAME, fr.read.getName());
        vcb.attribute(VCF_SCORE, String.format("%.5f", fr.score));
        vcb.attribute(VCF_FLAGS, fr.read.getFlags());
        vcb.attribute(VCF_MAPQ, fr.read.getMappingQuality());
        vcb.attribute(VCF_CIGAR, fr.read.getCigar().toString());
        vcb.attribute(VCF_READ_COUNT, fr.readCount);
        vcb.attribute(VCF_FILTERED_COUNT, fr.filteredCount);
        vcb.attribute(VCF_FC1, fr.nonIdentMBasesOnRead);
        vcb.attribute(VCF_FC2, fr.featuresOnRead);
        vcb.attribute(VCF_LENGTH, fr.read.getLength());
        vcb.attribute(VCF_EDIST, fr.refEditDistance);
        vcb.attribute(VCF_INDEX, fr.index);

        // median/mean quality on?
        if ( fmArgs.surroundingMediaQualitySize != null ) {
            vcb.attribute(VCF_SMQ_LEFT, fr.smqLeft);
            vcb.attribute(VCF_SMQ_RIGHT, fr.smqRight);
        }
        if ( fmArgs.surroundingMeanQualitySize != null ) {
            vcb.attribute(VCF_SMQ_LEFT_MEAN, fr.smqLeftMean);
            vcb.attribute(VCF_SMQ_RIGHT_MEAN, fr.smqRightMean);
        }

        for ( CopyAttrInfo info : copyAttrInfo ) {
            if ( fr.read.hasAttribute(info.name) ) {
                final String attrName = fmArgs.copyAttrPrefix + info.name;
                if ( info.type == VCFHeaderLineType.Integer ) {
                    vcb.attribute(attrName, fr.read.getAttributeAsInteger(info.name));
                } else if ( info.type == VCFHeaderLineType.Float ) {
                    vcb.attribute(attrName, fr.read.getAttributeAsFloat(info.name));
                } else {
                    vcb.attribute(attrName, fr.read.getAttributeAsString(info.name));
                }
            }
        }

        // validation mode?
        if ( fmArgs.reportAllAlts ) {
            if ( fr.scoreForBase != null ) {
                for (int baseIndex = 0; baseIndex < SequenceUtil.VALID_BASES_UPPER.length; baseIndex++) {
                    if (!Double.isNaN(fr.scoreForBase[baseIndex])) {
                        vcb.attribute(scoreNameForBase(baseIndex), String.format("%.5f", fr.scoreForBase[baseIndex]));
                    }
                }
            }
        }
        if ( fr.adjacentRefDiff ) {
            vcb.attribute(VCF_ADJACENT_REF_DIFF, true);
        }

        // build it!
        final VariantContext      vc = vcb.make();

        // write to file
        vcfWriter.add(vc);
    }

    private FeatureMapper buildMapper() {

        // build appropriate mapper
        if ( fmArgs.mappingFeature == FlowFeatureMapperArgumentCollection.MappingFeatureEnum.SNV ) {
            return new SNVMapper(fmArgs);
        } else {
            throw new GATKException("unsupported mappingFeature: " + fmArgs.mappingFeature);
        }
    }
}

