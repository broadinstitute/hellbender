package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang3.mutable.MutableInt;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.HashSet;
import java.util.Set;

/**
 * Annotate every variant in a VCF with the depth at that locus in a bam.  Note that this bam is *not* the bam
 * from which the vcf was derived, otherwise we would simply use the DP INFO field.
 *
 * <p>
 * In the CRSP sensitivity validation, we have a bam derived from a pool of 5, 10, or 20 samples
 * and a vcf of all known variants in those samples.  The pooled bam is a simulated tumor and
 * the vcf of individual variants is our truth data.  We annotate the truth data with the depth
 * in the pooled bam in order to bin the results of our sensitivity analysis by depth.
 * </p>
 *
 * <h3>Example</h3>
 *
 * <pre>
 * gatk --java-options "-Xmx4g" AnnotateVcfWithBamDepth \
 *   -V input.vcf \
 *   -I reads.bam \
 *   -O output.vcf
 * </pre>
 *
 * Created by David Benjamin on 1/30/17.
 */
@CommandLineProgramProperties(
        summary = "Annotate a vcf with a bam's read depth at each variant locus",
        oneLineSummary = "(Internal) Annotate a vcf with a bam's read depth at each variant locus",
        programGroup = VariantEvaluationProgramGroup.class
)
@DocumentedFeature
public class AnnotateVcfWithBamDepth extends VariantWalker {

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output filtered VCF file",
            optional=false)
    private final GATKPath outputVcf = null;

    private VariantContextWriter vcfWriter;

    public static final String POOLED_BAM_DEPTH_ANNOTATION_NAME = "BAM_DEPTH";

    @Override
    public void onTraversalStart() {
        final VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> headerLines = new HashSet<>(inputHeader.getMetaDataInSortedOrder());
        headerLines.add(new VCFInfoHeaderLine(POOLED_BAM_DEPTH_ANNOTATION_NAME, 1, VCFHeaderLineType.Integer, "pooled bam depth"));
        headerLines.addAll(getDefaultToolVCFHeaderLines());
        final VCFHeader vcfHeader = new VCFHeader(headerLines, inputHeader.getGenotypeSamples());
        vcfWriter = createVCFWriter(outputVcf);
        vcfWriter.writeHeader(vcfHeader);
    }

    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext refContext, final FeatureContext fc) {
        final MutableInt depth = new MutableInt(0);
        for (final GATKRead read : readsContext) {
            if (!read.failsVendorQualityCheck() && !read.isDuplicate() && !read.isUnmapped() && read.getEnd() > read.getStart() && new SimpleInterval(read).contains(vc) ) {
                depth.increment();
            }
        }
        vcfWriter.add(new VariantContextBuilder(vc).attribute(POOLED_BAM_DEPTH_ANNOTATION_NAME, depth.intValue()).make());
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }
}
