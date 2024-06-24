package org.broadinstitute.hellbender.tools.walkers.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.bed.FullBEDFeature;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SVInterval;
import org.broadinstitute.hellbender.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.gtf.*;

import java.io.File;
import java.util.*;

/**
 * Adds gene overlap, predicted functional consequence, and noncoding element overlap annotations to
 * a structural variant (SV) VCF from the GATK-SV pipeline.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         VCF containing structural variant (SV) records from the GATK-SV pipeline
 *     </li>
 *     <li>
 *         GTF file containing primary or canonical transcripts (optional; required for protein-coding gene,
 *         transcription start site, and promoter overlap annotations)
 *     </li>
 *     <li>
 *         BED file of noncoding elements in which the fourth column specifies the type of element
 *         (optional; required for noncoding element overlap annotations)
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Annotated VCF
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk SVAnnotate \
 *       -V structural.vcf.gz \
 *       --protein-coding-gtf canonical.gtf \
 *       --non-coding-bed noncoding.bed \
 *       -O annotated.vcf.gz
 * </pre>
 *
 * <h3>Annotation categories</h3>
 * <p>
 *     If a variant overlaps a gene, promoter, or noncoding element, the predicted functional impact will be annotated
 *     and the gene or noncoding element name listed. The list below describes the functional consequence categories
 *     in more detail. For complex variants (CPX), the annotations represent the union of the independent annotations
 *     of each component of the complex event according to its coordinates and SV type.
 * </p>
 *
 * <ul>
 *     <li><p><i>PREDICTED_LOF</i><br />
 *     Gene(s) on which the SV is predicted to have a loss-of-function (LOF) effect.
 *     The conditions for a LOF consequence depend on the SV type:
 *     <ul>
 *         <li>
 *             Deletion (DEL): the deletion overlaps any coding sequence (CDS) or the TSS
 *         </li>
 *         <li>
 *             Duplication (DUP): the duplication has both breakpoints in CDS, or one breakpoint in CDS and
 *             another in the 3' or 5' untranslated region (UTR)
 *         </li>
 *         <li>
 *             Insertion (INS): the insertion falls within CDS
 *         </li>
 *         <li>
 *             Inversion (INV): the inversion overlaps any coding sequence (CDS) or the TSS, except if it spans
 *             the entire gene (<i>PREDICTED_INV_SPAN</i>)
 *         </li>
 *         <li>
 *             Translocation (CTX): any translocation breakpoint falls within the transcript
 *         </li>
 *         <li>
 *             Multiallelic copy number variant (CNV): not annotated as <i>PREDICTED_LOF</i>.
 *             See <i>PREDICTED_MSV_EXON_OVERLAP</i>
 *         </li>
 *         <li>
 *             Breakend (BND): not annotated as <i>PREDICTED_LOF</i>. See <i>PREDICTED_BREAKEND_EXONIC</i>
 *         </li>
 *     </ul>
 *     </p></li>
 *     <li><p><i>PREDICTED_COPY_GAIN</i><br />
 *     Gene(s) on which the SV is predicted to have a copy-gain effect. This occurs when a duplication spans the entire
 *     transcript, from the first base of the 5' UTR to the last base of the 3' UTR. </p></li>
 *     <li><p><i>PREDICTED_INTRAGENIC_EXON_DUP</i><br />
 *     Gene(s) on which the SV is predicted to result in intragenic exonic duplication without breaking any coding
 *     sequences. This occurs when a duplication spans at least one coding exon and neither breakpoint is in CDS; both
 *     breakpoints are in UTR or intron. The result is that intact exons are duplicated within the boundaries of the
 *     gene body. </p></li>
 *     <li><p><i>PREDICTED_PARTIAL_EXON_DUP</i><br />
 *     Gene(s) where the duplication SV has one breakpoint in the coding sequence. This occurs when a duplication
 *     has exactly one breakpoint in CDS and the other breakpoint is in intron or UTR. When the duplication is in
 *     tandem, the result is that the endogenous copy of the breakpoint-harboring exon remains intact and a partial
 *     duplicate copy of that exon is also found elsewhere in the same gene.</p></li>
 *     <li><p><i>PREDICTED_TSS_DUP</i><br />
 *     Gene(s) for which the SV is predicted to duplicate the transcription start site (TSS). This occurs when a
 *     duplication has one breakpoint before the start of a transcript and the other breakpoint within the transcript.
 *     When the duplication is in tandem, the result is that there is one intact copy of the full endogenous gene, but
 *     an additional transcription start site is duplicated upstream (5') of the endogenous TSS. </p></li>
 *     <li><p><i>PREDICTED_DUP_PARTIAL</i><br />
 *     Gene(s) which are partially overlapped by an SV's duplication, but the transcription start site is not
 *     duplicated. The partial duplication occurs when a duplication has one breakpoint within the transcript and one
 *     breakpoint after the end of the transcript. When the duplication is in tandem, the result is that there is one
 *     intact copy of the full endogenous gene.</p></li>
 *     <li><p><i>PREDICTED_PARTIAL_DISPERSED_DUP</i><br />
 *     Gene(s) which are partially overlapped by the duplicated segment involved in an SV's dispersed duplication.
 *     This annotation is applied to a dispersed (non-tandem) duplication segment that is part of a complex SV if the
 *     duplicated segment overlaps part of a transcript but not the entire transcript (which would be a
 *     PREDICTED_COPY_GAIN event).</p></li>
 *     <li><p><i>PREDICTED_INV_SPAN</i><br />
 *     Gene(s) which are entirely spanned by an SV's inversion. A whole-gene inversion occurs when an inversion spans
 *     the entire transcript, from the first base of the 5' UTR to the last base of the 3' UTR. </p></li>
 *     <li><p><i>PREDICTED_MSV_EXON_OVERLAP</i><br />
 *     Gene(s) on which the multiallelic CNV would be predicted to have a LOF, INTRAGENIC_EXON_DUP, COPY_GAIN,
 *     DUP_PARTIAL, TSS_DUP, or PARTIAL_EXON_DUP annotation if the SV were biallelic. The functional impact of the
 *     multiallelic CNV on an individual sample depends on the copy number of the individual. </p></li>
 *     <li><p><i>PREDICTED_UTR</i><br />
 *     Gene(s) for which the SV is predicted to disrupt a UTR. This occurs when the SV has at least one breakpoint in
 *     a gene's 5' or 3' UTR but does not meet any of the criteria for a different gene-disrupting categorization
 *     above.</p></li>
 *     <li><p><i>PREDICTED_INTRONIC</i><br />
 *     Gene(s) where the SV was found to lie entirely within an intron. </p></li>
 *     <li><p><i>PREDICTED_BREAKEND_EXONIC</i><br />
 *     Gene(s) for which the SV breakend is predicted to fall in an exon. This category is reserved for breakend (BND)
 *     SVs with a breakpoint in CDS. </p></li>
 *     <li><p><i>PREDICTED_INTERGENIC</i><br />
 *     SV does not overlap any protein-coding gene transcripts in the GTF. </p></li>
 *     <li><p><i>PREDICTED_PROMOTER</i><br />
 *     Gene(s) for which the SV is predicted to overlap the promoter region. This occurs when the variant overlaps the
 *     predicted promoter region but does not overlap the transcript. The promoter region is inferred from the GTF as a
 *     window upstream of the TSS. The size of the window can be altered with the --promoter-window-length
 *     argument. </p></li>
 *     <li><p><i>PREDICTED_NEAREST_TSS</i><br />
 *     Nearest transcription start site to an intergenic variant. The gene with the nearest TSS to either side of the SV
 *     is annotated for intergenic variants that do not overlap any promoter regions. </p></li>
 *     <li><p><i>PREDICTED_NONCODING_SPAN</i><br />
 *     Class(es) of noncoding elements spanned by SV. </p></li>
 *     <li><p><i>PREDICTED_NONCODING_BREAKPOINT</i><br />
 *     Class(es) of noncoding elements disrupted by SV breakpoint. </p></li>
 *
 * </ul>
 */
@CommandLineProgramProperties(
        summary = "Adds predicted functional consequence, gene overlap, and noncoding element overlap annotations " +
                "to SV VCF from GATK-SV pipeline. " +
                "Input files are an SV VCF, a GTF file containing primary or canonical transcripts, " +
                "and a BED file containing noncoding elements. " +
                "Output file is an annotated SV VCF.",
        oneLineSummary = "Adds gene overlap and variant consequence annotations to SV VCF from GATK-SV pipeline",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
public final class SVAnnotate extends VariantWalker {
    public static final String PROTEIN_CODING_GTF_NAME = "protein-coding-gtf";
    public static final String PROMOTER_WINDOW_NAME = "promoter-window-length";
    public static final String NON_CODING_BED_NAME = "non-coding-bed";
    public static final String MAX_BND_LEN_NAME = "max-breakend-as-cnv-length";

    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output file (if not provided, defaults to STDOUT)",
            common = false,
            optional = true
    )
    private GATKPath outputFile = null;

    @Argument(
            fullName=PROTEIN_CODING_GTF_NAME,
            doc="Protein-coding GTF file containing primary or canonical transcripts (1-2 transcripts per gene only)",
            optional=true
    )
    private File proteinCodingGTFFile;

    @Argument(
            fullName=PROMOTER_WINDOW_NAME,
            doc="Promoter window (bp) upstream of TSS. Promoters will be inferred as the {window} bases upstream of the TSS. Default: 1000",
            minValue=0, optional=true
    )
    private int promoterWindow = 1000;

    @Argument(
            fullName=NON_CODING_BED_NAME,
            doc="BED file (with header) containing non-coding features. Columns: chrom, start, end, name, score (.), strand",
            optional=true
    )
    private File nonCodingBedFile;

    @Argument(
            fullName=MAX_BND_LEN_NAME,
            doc="Length in bp. Provide to annotate BNDs smaller than this size as deletions or duplications if applicable. Recommended value: < 2000000",
            minValue=0, optional=true
    )
    private int maxBreakendLen = -1;

    private VariantContextWriter vcfWriter = null;
    private SVIntervalTree<String> nonCodingIntervalTree;
    private SVAnnotateEngine.GTFIntervalTreesContainer gtfIntervalTrees;
    private SAMSequenceDictionary sequenceDictionary;
    private SVAnnotateEngine svAnnotateEngine;

    @Override
    public void onTraversalStart() {
        final VCFHeader header = getHeaderForVariants();
        // get contigs from VCF
        sequenceDictionary = header.getSequenceDictionary();
        // TODO: more elegant way to make reference inputs optional than checking 10x?
        // Load protein-coding GTF data into memory as interval tree of transcripts if GTF provided
        if (proteinCodingGTFFile != null) {
            final FeatureDataSource<GencodeGtfGeneFeature> proteinCodingGTFSource = new FeatureDataSource<>(proteinCodingGTFFile);
            gtfIntervalTrees = buildIntervalTreesFromGTF(proteinCodingGTFSource, sequenceDictionary, promoterWindow);
        }

        // Load noncoding BED file into memory as interval tree of noncoding elements if BED provided
        if (nonCodingBedFile != null) {
            final FeatureDataSource<FullBEDFeature> nonCodingSource = new FeatureDataSource<>(nonCodingBedFile);
            nonCodingIntervalTree = buildIntervalTreeFromBED(nonCodingSource, sequenceDictionary);
        }

        vcfWriter = createVCFWriter(outputFile);
        updateAndWriteHeader(header);

        svAnnotateEngine = new SVAnnotateEngine(gtfIntervalTrees, nonCodingIntervalTree, sequenceDictionary,
                maxBreakendLen);
    }

    /**
     * Check if strand of transcript from GTF is negative
     * @param transcript - protein-coding GTF transcript to check
     * @return - boolean True if negative strand, False if positive
     */
    @VisibleForTesting
    protected static boolean isNegativeStrand(final GencodeGtfTranscriptFeature transcript) {
        return transcript.getGenomicStrand().equals(Strand.decode("-"));
    }

    /**
     * Get transcription start site from GTF transcript: first (5') base of transcript
     * @param transcript - protein-coding GTF transcript
     * @return - int position of TSS on transcript (start of transcript if + strand, end if - strand)
     */
    @VisibleForTesting
    protected static int getTranscriptionStartSite(final GencodeGtfTranscriptFeature transcript) {
        final boolean negativeStrand = isNegativeStrand(transcript);
        return negativeStrand ? transcript.getEnd() : transcript.getStart();  // GTF codec reassigns start to be earlier in contig;
    }

    /**
     * Get promoter interval for GTF transcript: user-defined window upstream (5') of TSS
     * @param transcript - protein-coding GTF transcript
     * @param promoterWindow - size of promoter window in bp
     * @return - SimpleInterval representing promoter region of size promoterWindow upstream of TSS
     */
    @VisibleForTesting
    protected static SimpleInterval getPromoterInterval(final GencodeGtfTranscriptFeature transcript,
                                                        final int promoterWindow) {
        final int tss = getTranscriptionStartSite(transcript);
        final boolean negativeStrand = isNegativeStrand(transcript);
        final int promoterLeft = negativeStrand ? tss + 1 : Math.max(tss - promoterWindow, 1);
        final int promoterRight = negativeStrand ? tss + promoterWindow : tss - 1;
        return new SimpleInterval(transcript.getContig(), promoterLeft, promoterRight);
    }



    /**
     * Builds transcript, TSS, and promoter interval trees from protein-coding GTF
     * @param proteinCodingGTFSource - GTF as FeatureDataSource
     * @param sequenceDictionary - SAMSequenceDictionary for VCF
     * @param promoterWindow - size of promoter window in bp
     * @return - container class packaging transcript, promoter, and TSS interval trees for annotation
     */
    @VisibleForTesting
    protected static SVAnnotateEngine.GTFIntervalTreesContainer buildIntervalTreesFromGTF(
            final FeatureDataSource<GencodeGtfGeneFeature> proteinCodingGTFSource,
            final SAMSequenceDictionary sequenceDictionary, final int promoterWindow
    ) {
        final SVIntervalTree<GencodeGtfTranscriptFeature> transcriptIntervalTree = new SVIntervalTree<>();
        final SVIntervalTree<String> promoterIntervalTree = new SVIntervalTree<>();
        final SVIntervalTree<String> transcriptionStartSiteTree = new SVIntervalTree<>();
        for (final GencodeGtfGeneFeature gene : proteinCodingGTFSource) {
            final List<GencodeGtfTranscriptFeature> transcriptsForGene = gene.getTranscripts();
            for (GencodeGtfTranscriptFeature transcript : transcriptsForGene) {
                final int contigID = sequenceDictionary.getSequenceIndex(transcript.getContig());
                if (contigID < 0) {
                    continue; // if GTF input contains chromosome not in VCF sequence dictionary, just ignore it
                }
                transcriptIntervalTree.put(SVUtils.locatableToSVInterval(transcript, sequenceDictionary), transcript);
                final String geneName = transcript.getGeneName();
                final int tss = getTranscriptionStartSite(transcript);
                transcriptionStartSiteTree.put(new SVInterval(contigID, tss, tss + 1), geneName);
                SimpleInterval promoterInterval = getPromoterInterval(transcript, promoterWindow);
                promoterIntervalTree.put(SVUtils.locatableToSVInterval(promoterInterval, sequenceDictionary), geneName);
            }
        }
        return new SVAnnotateEngine.GTFIntervalTreesContainer(transcriptIntervalTree, promoterIntervalTree, transcriptionStartSiteTree);
    }

    /**
     * Builds interval tree of noncoding elements to annotate from BED file input
     * @param BEDSource - noncoding element BED file as FeatureDataSource
     * @param sequenceDictionary - SAMSequenceDictionary for VCF
     * @return - SVIntervalTree of nonocoding elements for annotation
     */
    @VisibleForTesting
    protected static SVIntervalTree<String> buildIntervalTreeFromBED(final FeatureDataSource<FullBEDFeature> BEDSource,
                                                            final SAMSequenceDictionary sequenceDictionary) {
        final SVIntervalTree<String> BEDIntervalTree = new SVIntervalTree<>();
        for (final FullBEDFeature feature : BEDSource) {
            // BED feature class already does start+1 conversion to 1-based closed interval
            try {
                BEDIntervalTree.put(SVUtils.locatableToSVInterval(feature, sequenceDictionary), feature.getName());
            } catch (IllegalArgumentException e) {
                continue;  // if BED input contains chromosome not in VCF sequence dictionary, just ignore it
            }
        }
        return BEDIntervalTree;
    }

    /**
     * Adds SV functional annotation INFO keys to VCF header
     * @param header - starting VCF header to which to add INFO keys
     */
    private void addAnnotationInfoKeysToHeader(final VCFHeader header) {
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.LOF, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) on which the SV is predicted to have a loss-of-function effect."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.INT_EXON_DUP, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) on which the SV is predicted to result in intragenic exonic duplication without breaking any coding sequences."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.COPY_GAIN, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) on which the SV is predicted to have a copy-gain effect."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TSS_DUP, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) for which the SV is predicted to duplicate the transcription start site."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.DUP_PARTIAL, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) which are partially overlapped by an SV's duplication, but the transcription start site is not duplicated."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.INTRONIC, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) where the SV was found to lie entirely within an intron."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.PARTIAL_EXON_DUP, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) where the duplication SV has one breakpoint in the coding sequence."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.INV_SPAN, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) which are entirely spanned by an SV's inversion."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.UTR, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) for which the SV is predicted to disrupt a UTR."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.MSV_EXON_OVERLAP, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) on which the multiallelic SV would be predicted to have a LOF, INTRAGENIC_EXON_DUP, COPY_GAIN, DUP_PARTIAL, TSS_DUP, or PARTIAL_EXON_DUP annotation if the SV were biallelic."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.PROMOTER, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) for which the SV is predicted to overlap the promoter region."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.BREAKEND_EXON, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) for which the SV breakend is predicted to fall in an exon."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.INTERGENIC, 0, VCFHeaderLineType.Flag, "SV does not overlap any protein-coding genes."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.NONCODING_SPAN, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Class(es) of noncoding elements spanned by SV."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.NONCODING_BREAKPOINT, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Class(es) of noncoding elements disrupted by SV breakpoint."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.NEAREST_TSS, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Nearest transcription start site to an intergenic variant."));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.PARTIAL_DISPERSED_DUP, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) overlapped partially by the duplicated interval involved in a dispersed duplication event in a complex SV."));

    }

    /**
     * Updates VCF header with SV annotation INFO keys and writes header
     * @param header - starting VCF header
     */
    private void updateAndWriteHeader(final VCFHeader header) {
        addAnnotationInfoKeysToHeader(header);
        vcfWriter.writeHeader(header);
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext,
                      final FeatureContext featureContext) {
        final VariantContext annotated = svAnnotateEngine.createAnnotatedStructuralVariantContext(variant);
        vcfWriter.add(annotated);
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }
}
