package org.broadinstitute.hellbender.tools.walkers.sv;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.tribble.Feature;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.bed.FullBEDFeature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.compress.utils.Sets;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.codecs.gtf.*;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.*;

import static java.util.Objects.isNull;
import static org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata.buildContigNameToIDMap;

/**
 * Adds gene overlap and variant consequence annotations to SV VCF from GATK-SV pipeline
 * Input files are an SV VCF and a GTF file containing primary or canonical transcripts
 * Output file is an annotated SV VCF
 */
@CommandLineProgramProperties(
        summary = "Adds gene overlap and variant consequence annotations to SV VCF from GATK-SV pipeline." +
                "Input files are an SV VCF and a GTF file containing primary or canonical transcripts." +
                "Output file is an annotated SV VCF.",
        oneLineSummary = "Adds gene overlap and variant consequence annotations to SV VCF from GATK-SV pipeline",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
public final class SVAnnotate extends VariantWalker {

    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output file (if not provided, defaults to STDOUT)",
            common = false,
            optional = true
    )
    private GATKPath outputFile = null;

    @Argument(
            fullName="proteinCodingGTF",
            doc="protein-coding GTF file (canonical only)",
            optional=true
    )
    private File proteinCodingGTFFile;

    @Argument(
            fullName="promoterBed",
            doc="BED file (with header) containing promoter regions. Columns: chrom, start, end, name (gene to which the promoter corresponds), score (.), strand",
            optional=true
    )
    private File promoterBedFile;

    @Argument(
            fullName="nonCodingBed",
            doc="BED file (with header) containing non-coding features. Columns: chrom, start, end, name, score (.), strand",
            optional=true
    )
    private File nonCodingBedFile;

    private VariantContextWriter vcfWriter = null;
    private OverlapDetector<GencodeGtfGeneFeature> gtfOverlapDetector;
    private OverlapDetector<FullBEDFeature> promoterOverlapDetector;
    private OverlapDetector<FullBEDFeature> nonCodingOverlapDetector;
    private SVIntervalTree<String> transcriptionStartSiteTree;
    private final Set<String> MSVExonOverlapClassifications = Sets.newHashSet(GATKSVVCFConstants.LOF, GATKSVVCFConstants.INT_EXON_DUP, GATKSVVCFConstants.DUP_PARTIAL, GATKSVVCFConstants.PARTIAL_EXON_DUP, GATKSVVCFConstants.COPY_GAIN);
    private Map<String, Integer> contigNameToID;
    private SAMSequenceDictionary sequenceDictionary;
    private int maxContigLength;
    private enum StructuralVariantAnnotationType {
        DEL,
        DUP,
        INS,
        INV,
        CPX,
        BND,
        CTX,
        CNV
    }

    // mini class for SV intervals (type and segment) within CPX events
    private static final class SVSegment {
        private final StructuralVariantAnnotationType intervalSVType;
        private final SimpleInterval interval;
        private SVSegment(final StructuralVariantAnnotationType svType, final SimpleInterval interval) {
            this.intervalSVType = svType;
            this.interval = interval;
        }
    }

    @Override
    public void onTraversalStart() {
        final VCFHeader header = getHeaderForVariants();
        sequenceDictionary = header.getSequenceDictionary();
        maxContigLength = sequenceDictionary.getSequences().stream().mapToInt(SAMSequenceRecord::getSequenceLength).max().getAsInt() + 1;
        contigNameToID = buildContigNameToIDMap(sequenceDictionary);

        final FeatureDataSource<GencodeGtfGeneFeature> proteinCodingGTFSource = new FeatureDataSource<>(proteinCodingGTFFile);
        gtfOverlapDetector = OverlapDetector.create(Lists.newArrayList(proteinCodingGTFSource));

        buildTranscriptionStartSiteTree(proteinCodingGTFSource);

        final FeatureDataSource<FullBEDFeature> promoterSource = new FeatureDataSource<>(promoterBedFile);
        promoterOverlapDetector = OverlapDetector.create(Lists.newArrayList(promoterSource));

        final FeatureDataSource<FullBEDFeature> nonCodingSource = new FeatureDataSource<>(nonCodingBedFile);
        nonCodingOverlapDetector = OverlapDetector.create(Lists.newArrayList(nonCodingSource));


        vcfWriter = createVCFWriter(outputFile);
        updateAndWriteHeader(header);
    }

    private void buildTranscriptionStartSiteTree(FeatureDataSource<GencodeGtfGeneFeature> proteinCodingGTFSource) {
        transcriptionStartSiteTree = new SVIntervalTree<>();
        for (final GencodeGtfGeneFeature gene : proteinCodingGTFSource) {
            final List<GencodeGtfTranscriptFeature> transcriptsForGene = gene.getTranscripts();
            for (GencodeGtfTranscriptFeature transcript : transcriptsForGene) {
                final int start = transcript.getGenomicStrand().equals(Strand.decode("-")) ? transcript.getEnd() : transcript.getStart();
                final int end = start + 1;
                transcriptionStartSiteTree.put(new SVInterval(contigNameToID.get(transcript.getContig()), start, end), transcript.getGeneName());
            }
        }
    }

    private void addAnnotationInfoKeysToHeader(final VCFHeader header) {
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.LOF, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) on which the SV is predicted to have a loss-of-function effect."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.INT_EXON_DUP, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) on which the SV is predicted to result in intragenic exonic duplication."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.COPY_GAIN, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) on which the SV is predicted to have a copy-gain effect."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.DUP_PARTIAL, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) which are partially overlapped by an SV's duplication."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.INTRONIC, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) where the SV was found to lie entirely within an intron."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.PARTIAL_EXON_DUP, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) where the SV was found to partially overlap a single exon."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.INV_SPAN, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) which are entirely spanned by an SV's inversion."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.UTR, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) for which the SV is predicted to disrupt a UTR."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.MSV_EXON_OVERLAP, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) on which the multiallelic SV would be predicted to have a LOF, INTRAGENIC_EXON_DUP, COPY_GAIN, DUP_PARTIAL, or PARTIAL_EXON_DUP annotation if the SV were biallelic."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.PROMOTER, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) for which the SV is predicted to overlap the promoter region."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.BREAKEND_EXON, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) for which the SV breakend is predicted to fall in an exon."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.INTERGENIC, 0, VCFHeaderLineType.Flag, "SV does not overlap coding sequence."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.NONCODING_SPAN, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Class(es) of noncoding elements spanned by SV."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.NONCODING_BREAKPOINT, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Class(es) of noncoding elements disrupted by SV breakpoint."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.NEAREST_TSS, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Nearest transcription start site to intragenic variants."));

    }

    private void updateAndWriteHeader(VCFHeader header) {
        addAnnotationInfoKeysToHeader(header);
        vcfWriter.writeHeader(header);
    }

    protected static boolean variantSpansFeature(final SimpleInterval variantInterval, final SimpleInterval featureInterval) {
        return variantInterval.contains(featureInterval);
    }

    protected static int countBreakendsInsideFeature(final SimpleInterval variantInterval, final SimpleInterval featureInterval) {
        int count = 0;
        if (variantInterval.getContig().equals(featureInterval.getContig())) {
            if (variantInterval.getStart() >= featureInterval.getStart() && variantInterval.getStart() <= featureInterval.getEnd()) {
                count++;
            }
            if (variantInterval.getEnd() >= featureInterval.getStart() && variantInterval.getEnd() <= featureInterval.getEnd()) {
                count++;
            }
        }
        return count;
    }

    protected static boolean variantOverlapsFeature(final SimpleInterval variantInterval, final SimpleInterval featureInterval) {
        return IntervalUtils.overlaps(variantInterval, featureInterval);
    }

    private static void updateVariantConsequenceDict(final Map<String, Set<String>> variantConsequenceDict, final String key, final String value) {
        variantConsequenceDict.putIfAbsent(key, new HashSet<>());
        variantConsequenceDict.get(key).add(value);
    }

    private String annotateDeletionOrInsertion(final SimpleInterval variantInterval, final GencodeGtfTranscriptFeature gtfTranscript) {
        final List<GencodeGtfFeature> gtfFeaturesForTranscript = gtfTranscript.getAllFeatures();
        String consequence = GATKSVVCFConstants.INTRONIC;
        for (GencodeGtfFeature gtfFeature : gtfFeaturesForTranscript) {
            final SimpleInterval featureInterval = new SimpleInterval(gtfFeature);
            if (!variantOverlapsFeature(variantInterval, featureInterval)) {
                continue;
            }
            if (gtfFeature.getFeatureType() == GencodeGtfFeature.FeatureType.CDS) {
                consequence = GATKSVVCFConstants.LOF;
                break;
            } else if (gtfFeature.getFeatureType() == GencodeGtfFeature.FeatureType.UTR) {
                consequence = GATKSVVCFConstants.UTR;
            }
        }
        return consequence;
    }

    private String annotateDuplication(final SimpleInterval variantInterval, final GencodeGtfTranscriptFeature gtfTranscript) {
        String consequence = GATKSVVCFConstants.INTRONIC;
        final SimpleInterval transcriptInterval = new SimpleInterval(gtfTranscript);
        if (variantSpansFeature(variantInterval, transcriptInterval)) {
            consequence = GATKSVVCFConstants.COPY_GAIN;
        } else if (countBreakendsInsideFeature(variantInterval, transcriptInterval) == 1) {
            consequence = GATKSVVCFConstants.DUP_PARTIAL;
        } else {
            // both breakpoints inside transcript
            final List<GencodeGtfFeature> gtfFeaturesForTranscript = gtfTranscript.getAllFeatures();
            int numBreakpointsInExon = 0;  // TODO: CDS or exon?
            int numBreakpointsInUTR = 0;
            int numExonsSpanned = 0;
            for (GencodeGtfFeature gtfFeature : gtfFeaturesForTranscript) {
                final SimpleInterval featureInterval = new SimpleInterval(gtfFeature);
                if (!variantOverlapsFeature(variantInterval, featureInterval)) {
                    continue;
                }
                if (gtfFeature.getFeatureType() == GencodeGtfFeature.FeatureType.EXON) {
                    if (variantSpansFeature(variantInterval, featureInterval)) {
                        numExonsSpanned++;  // TODO: CDS or exon? may differ from breakpoints
                    } else {
                        numBreakpointsInExon = numBreakpointsInExon + countBreakendsInsideFeature(variantInterval, featureInterval);
                    }
                } else if (gtfFeature.getFeatureType() == GencodeGtfFeature.FeatureType.UTR) {
                    numBreakpointsInUTR = numBreakpointsInUTR + countBreakendsInsideFeature(variantInterval, featureInterval);
                }
            }
            if (numBreakpointsInExon == 2) {
                consequence = GATKSVVCFConstants.LOF;
            } else if (numExonsSpanned > 0) {
                consequence = GATKSVVCFConstants.INT_EXON_DUP;  // formerly DUP_LOF - consider INTERNAL
            } else if (numBreakpointsInExon == 1) {
                consequence = GATKSVVCFConstants.PARTIAL_EXON_DUP;  // new category - could collapse with DUP_PARTIAL
            } else if (numBreakpointsInUTR > 0) {
                consequence = GATKSVVCFConstants.UTR;
            }
        }
        return consequence;
    }

    private String annotateCopyNumberVariant(final SimpleInterval variantInterval, final GencodeGtfTranscriptFeature gtfTranscript) {
        String consequence = annotateDuplication(variantInterval, gtfTranscript);
        if (MSVExonOverlapClassifications.contains(consequence)) {
            return GATKSVVCFConstants.MSV_EXON_OVERLAP;
        } else {
            return consequence;  // TODO: MCNV classifications ???
        }
    }

    private String annotateInversion(final SimpleInterval variantInterval, final GencodeGtfTranscriptFeature gtfTranscript) {
        String consequence = GATKSVVCFConstants.INTRONIC;
        final SimpleInterval transcriptInterval = new SimpleInterval(gtfTranscript);
        if (variantSpansFeature(variantInterval, transcriptInterval)) {
            consequence = GATKSVVCFConstants.INV_SPAN;
        } else if (countBreakendsInsideFeature(variantInterval, transcriptInterval) == 1) {
            consequence = GATKSVVCFConstants.LOF;
        } else {
            // both breakpoints inside transcript
            final List<GencodeGtfFeature> gtfFeaturesForTranscript = gtfTranscript.getAllFeatures();
            for (GencodeGtfFeature gtfFeature : gtfFeaturesForTranscript) {
                final SimpleInterval featureInterval = new SimpleInterval(gtfFeature);
                if (!variantOverlapsFeature(variantInterval, featureInterval)) {
                    continue;
                }
                // TODO: if overlaps exon, it's LOF unless both breakpoints are in the same UTR?
                if (gtfFeature.getFeatureType() == GencodeGtfFeature.FeatureType.EXON) {
                    consequence = GATKSVVCFConstants.LOF;  // TODO: CDS or exon here?
                } else if (gtfFeature.getFeatureType() == GencodeGtfFeature.FeatureType.UTR) {
                    if (countBreakendsInsideFeature(variantInterval, featureInterval) == 2) {
                        consequence = GATKSVVCFConstants.UTR;
                    }
                }
            }
        }
        return consequence;
    }

    private String annotateTranslocation(final SimpleInterval variantInterval, final GencodeGtfTranscriptFeature gtfTranscript) {
        // already checked for transcript overlap, and if a translocation breakpoint falls inside a gene it's automatically LOF
        // TODO: eliminate unnecessary function or keep for aesthetics/future flexibility?
        return GATKSVVCFConstants.LOF;
    }

    private String annotateBreakend(final SimpleInterval variantInterval, final GencodeGtfTranscriptFeature gtfTranscript) {
        String consequence = annotateDeletionOrInsertion(variantInterval, gtfTranscript);
        if (consequence.equals(GATKSVVCFConstants.LOF)) {
            consequence = GATKSVVCFConstants.BREAKEND_EXON;
        }
        return consequence;
    }

    private void annotateTranscript(final SimpleInterval variantInterval, final StructuralVariantAnnotationType svType, final GencodeGtfTranscriptFeature transcript, final Map<String, Set<String>> variantConsequenceDict) {
        if (!variantOverlapsFeature(variantInterval, new SimpleInterval(transcript))) {
            return;
        }
        String consequence = null;
        if (svType.equals(StructuralVariantAnnotationType.DEL) || svType.equals(StructuralVariantAnnotationType.INS)) {
            consequence = annotateDeletionOrInsertion(variantInterval, transcript);
        } else if (svType.equals(StructuralVariantAnnotationType.DUP)) {
            consequence = annotateDuplication(variantInterval, transcript);
        } else if (svType.equals(StructuralVariantAnnotationType.CNV)) {
            consequence = annotateCopyNumberVariant(variantInterval,transcript);
        } else if (svType.equals(StructuralVariantAnnotationType.INV)) {
            consequence = annotateInversion(variantInterval, transcript);
        } else if (svType.equals(StructuralVariantAnnotationType.CTX)) {
            consequence = annotateTranslocation(variantInterval, transcript);
        } else if (svType.equals(StructuralVariantAnnotationType.BND)) {
            consequence = annotateBreakend(variantInterval, transcript);
        }

        if (consequence != null) {
            updateVariantConsequenceDict(variantConsequenceDict, consequence, transcript.getGeneName());
        }
    }

    /**
     Annotate promoter overlaps and return a boolean: true if there are any promoter overlaps, false if not
     */
    private boolean annotatePromoterOverlaps(final SimpleInterval variantInterval, final Map<String, Set<String>> variantConsequenceDict) {
        boolean anyPromoterOverlaps = false;
        final Set<FullBEDFeature> promotersForVariant = promoterOverlapDetector.getOverlaps(variantInterval);
        final Set<String> codingAnnotationGenes = new HashSet<>();
        variantConsequenceDict.values().forEach(codingAnnotationGenes::addAll);
        for (FullBEDFeature promoter : promotersForVariant) {
            if (!codingAnnotationGenes.contains(promoter.getName())) {
                updateVariantConsequenceDict(variantConsequenceDict, GATKSVVCFConstants.PROMOTER, promoter.getName());
                anyPromoterOverlaps = true;
            }
        }
        return anyPromoterOverlaps;
    }

    private void annotateNonCodingOverlaps(final SimpleInterval variantInterval, final Map<String, Set<String>> variantConsequenceDict) {
        final Set<FullBEDFeature> nonCodingFeaturesForVariant = nonCodingOverlapDetector.getOverlaps(variantInterval);
        for (FullBEDFeature feature : nonCodingFeaturesForVariant) {
            String consequence = GATKSVVCFConstants.NONCODING_BREAKPOINT;
            if (variantSpansFeature(variantInterval, new SimpleInterval(feature))) {
                consequence = GATKSVVCFConstants.NONCODING_SPAN;
            }
            updateVariantConsequenceDict(variantConsequenceDict, consequence, feature.getName());
        }
    }


    protected static void annotateNearestTranscriptionStartSite(final SimpleInterval variantInterval, final Map<String, Set<String>> variantConsequenceDict, SVIntervalTree<String> transcriptionStartSiteTree, int maxContigLength, int variantContigID) {
        // TODO: keep all nearest TSS for dispersed CPX / CTX or choose closest?
        // TODO: will start < end ever? Shouldn't at this point in the pipeline
        SVIntervalTree.Entry<String> nearestBefore = transcriptionStartSiteTree.max(new SVInterval(variantContigID, variantInterval.getStart(), variantInterval.getEnd()));
        SVIntervalTree.Entry<String> nearestAfter = transcriptionStartSiteTree.min(new SVInterval(variantContigID, variantInterval.getStart(), variantInterval.getEnd()));
        // nearest TSS only "valid" for annotation if non-null and on the same contig as the variant
        boolean beforeInvalid = (isNull(nearestBefore) || nearestBefore.getInterval().getContig() != variantContigID );
        boolean afterInvalid = (isNull(nearestAfter) || nearestAfter.getInterval().getContig() != variantContigID );
        // if at least one result is valid, keep one with shorter distance
        if (!(beforeInvalid && afterInvalid)) {
            // if result is invalid, set distance to longest contig length so that other TSS will be kept
            int distanceBefore = beforeInvalid ? maxContigLength : variantInterval.getStart() - nearestBefore.getInterval().getEnd();
            int distanceAfter = afterInvalid ? maxContigLength : nearestAfter.getInterval().getStart() - variantInterval.getEnd();
            String nearestTSSGeneName = (distanceBefore < distanceAfter) ? nearestBefore.getValue() : nearestAfter.getValue();
            updateVariantConsequenceDict(variantConsequenceDict, GATKSVVCFConstants.NEAREST_TSS, nearestTSSGeneName);
        }
        // TODO: return consequence instead?
    }

    private StructuralVariantAnnotationType getSVType(final VariantContext variant) {
        // TODO: haha majorly clean this up
        // return variant.getStructuralVariantType().name();
        final Allele alt = variant.getAlternateAllele(0); // TODO: any chance of multiallelic alt field for SV?
        if (alt.isBreakpoint()) {
            return StructuralVariantAnnotationType.BND;
        } else if (alt.isSingleBreakend()) {
            throw new IllegalArgumentException("what even is single breakend??: " + alt);
        } else if (alt.isSymbolic()) {
            if (alt.toString().contains("INS")) {
                return StructuralVariantAnnotationType.INS;
            } else {
                return StructuralVariantAnnotationType.valueOf(alt.toString().substring(1, alt.toString().length()-1));  // assume <SVTYPE>
            }
        } else {
            throw new IllegalArgumentException("Unexpected ALT allele: " + alt);
        }
    }

    private void annotateSVSegment(final SimpleInterval variantInterval, final StructuralVariantAnnotationType svType, final Map<String, Set<String>> variantConsequenceDict) {
        final Set<GencodeGtfGeneFeature> gtfGenesForVariant = gtfOverlapDetector.getOverlaps(variantInterval);
        for (GencodeGtfGeneFeature geneOverlapped : gtfGenesForVariant) {
            final List<GencodeGtfTranscriptFeature> transcriptsForGene = geneOverlapped.getTranscripts();
            for (GencodeGtfTranscriptFeature transcript : transcriptsForGene) {
                annotateTranscript(variantInterval, svType, transcript, variantConsequenceDict);
            }
        }
    }

    private List<SVSegment> getSVSegments(VariantContext variant, StructuralVariantAnnotationType overallSVType) {
        final List<SVSegment> intervals = new ArrayList<>();
        if (overallSVType.equals(StructuralVariantAnnotationType.CPX)) {
            final List<String> cpxIntervalsString = variant.getAttributeAsStringList(GATKSVVCFConstants.CPX_INTERVALS, "NONE");
            for (String cpxInterval : cpxIntervalsString) {
                final String[] parsed = cpxInterval.split("_");
                final StructuralVariantAnnotationType svTypeForInterval = StructuralVariantAnnotationType.valueOf(parsed[0]);
                final SimpleInterval interval = new SimpleInterval(parsed[1]);
                intervals.add(new SVSegment(svTypeForInterval, interval));
            }
        } else if (overallSVType.equals(StructuralVariantAnnotationType.CTX)) {
            intervals.add(new SVSegment(overallSVType, new SimpleInterval(variant)));
            // annotate both breakpoints of translocation - CHR2:END2-END2
            intervals.add(new SVSegment(overallSVType, new SimpleInterval(variant.getAttributeAsString(GATKSVVCFConstants.END_CONTIG_ATTRIBUTE, "NONE"), variant.getAttributeAsInt(GATKSVVCFConstants.END_CONTIG_POSITION, 0), variant.getAttributeAsInt(GATKSVVCFConstants.END_CONTIG_POSITION, 0))));
        } else {
            intervals.add(new SVSegment(overallSVType, new SimpleInterval(variant)));
        }

        return intervals;
    }

    private Map<String,String> formatVariantConsequenceDict(Map<String,Set<String>> variantConsequenceDict) {
        Map<String,String> formatted = new HashMap<>();
        for (String consequence : variantConsequenceDict.keySet()) {
            List<String> sortedGenes = new ArrayList<>(variantConsequenceDict.get(consequence));
            Collections.sort(sortedGenes);
            formatted.put(consequence, String.join(",", sortedGenes));
        }
        return formatted;
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        final Map<String, Set<String>> variantConsequenceDict = new HashMap<>();
        final StructuralVariantAnnotationType overallSVType = getSVType(variant);
        final List<SVSegment> svSegments = getSVSegments(variant, overallSVType);
        for (SVSegment svSegment : svSegments) {
            annotateSVSegment(svSegment.interval, svSegment.intervalSVType, variantConsequenceDict);
        }

        // if variant consequence dictionary is empty (no protein-coding annotations), apply INTERGENIC flag
        boolean noCodingAnnotations = variantConsequenceDict.isEmpty();

        // then annotate promoter overlaps and non-coding feature overlaps
        boolean anyPromoterOverlaps = false;
        for (SVSegment svSegment : svSegments) {
            anyPromoterOverlaps = anyPromoterOverlaps || annotatePromoterOverlaps(svSegment.interval, variantConsequenceDict);
            annotateNonCodingOverlaps(svSegment.interval, variantConsequenceDict);
        }

        // annotate nearest TSS for intergenic variants with no promoter overlaps
        if (!anyPromoterOverlaps && noCodingAnnotations) {
            for (SVSegment svSegment : svSegments) {
                annotateNearestTranscriptionStartSite(svSegment.interval, variantConsequenceDict, transcriptionStartSiteTree, maxContigLength, contigNameToID.get(svSegment.interval.getContig()));
            }
        }

        VariantContextBuilder vcb = new VariantContextBuilder(variant);
        vcb.putAttributes(formatVariantConsequenceDict(variantConsequenceDict));
        vcb.attribute(GATKSVVCFConstants.INTERGENIC, noCodingAnnotations);
        vcfWriter.add(vcb.make());
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }
}
