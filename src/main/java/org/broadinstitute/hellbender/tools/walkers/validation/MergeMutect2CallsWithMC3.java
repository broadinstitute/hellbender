package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.collections4.Predicate;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.AlleleConcordanceWalker;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.util.*;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "UNSUPPORTED.  FOR EVALUATION ONLY. Merge M2 (eval) calls with MC3 (truth)",
        oneLineSummary = "UNSUPPORTED.  FOR EVALUATION ONLY. Merge M2 calls with MC",
        programGroup = VariantEvaluationProgramGroup.class
)
@ExperimentalFeature
public class MergeMutect2CallsWithMC3 extends AlleleConcordanceWalker {

    @Argument(doc = "Merged vcf.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    protected GATKPath outputVcf;

    public static final String CENTERS_KEY = "CENTERS";
    public static final String M2_CENTER_NAME = "M2";
    public static final String MC3_REF_COUNT_KEY = "NREF";
    public static final String MC3_ALT_COUNT_KEY = "NALT";

    public static final String M2_FILTERS_KEY = "M2_FILTERS";
    public static final VCFInfoHeaderLine M2_FILTERS_HEADER_LINE =
            new VCFInfoHeaderLine(M2_FILTERS_KEY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "M2 filters applied to variant.");

    private VariantContextWriter vcfWriter;
    private String tumorSample;

    @Override
    protected Predicate<VariantContext> makeTruthVariantFilter() {
        return vc -> true;
    }

    @Override
    protected Predicate<VariantContext> makeEvalVariantFilter() { return vc -> true; }

    @Override
    public void onTraversalStart() {
        final Set<VCFHeaderLine> headerLines = new HashSet<>(getTruthHeader().getMetaDataInSortedOrder());

        VCFStandardHeaderLines.addStandardFormatLines(headerLines, true,
                VCFConstants.GENOTYPE_KEY,
                VCFConstants.GENOTYPE_ALLELE_DEPTHS);

        headerLines.addAll(getDefaultToolVCFHeaderLines());
        headerLines.add(M2_FILTERS_HEADER_LINE);
        tumorSample = getEvalHeader().getMetaDataLine(Mutect2Engine.TUMOR_SAMPLE_KEY_IN_VCF_HEADER).getValue();
        final VCFHeader mergedHeader = new VCFHeader(headerLines, Collections.singletonList(tumorSample));
        vcfWriter = createVCFWriter(outputVcf);
        vcfWriter.writeHeader(mergedHeader);
    }

    @Override
    protected void apply(final AlleleTruthVersusEval truthVersusEval, final ReadsContext readsContext, final ReferenceContext refContext) {
        final ConcordanceState concordanceState = truthVersusEval.getConcordance();

        // get AD from M2 if available, otherwise from NALT and NREF MC3 info fields
        final int[] tumorAlleleCounts = truthVersusEval.hasEval() ? truthVersusEval.getEval().getGenotype(tumorSample).getAD() :
                new int[] {truthVersusEval.getTruth().getAttributeAsInt(MC3_REF_COUNT_KEY, 0),
                        truthVersusEval.getTruth().getAttributeAsInt(MC3_ALT_COUNT_KEY, 0)};

        final Genotype tumorGenotype = new GenotypeBuilder(tumorSample, truthVersusEval.getTruthIfPresentElseEval().getAlleles())
                .AD(tumorAlleleCounts)
                .make();
        final List<Genotype> genotypes = new ArrayList<>(Arrays.asList(tumorGenotype));

        switch (concordanceState) {
            case TRUE_POSITIVE: // variant present in MC3 and unfiltered in M2.  Add M2 to CENTERS field.
                vcfWriter.add(makeVariantContextBuilderWithM2Center(truthVersusEval.getTruth()).genotypes(genotypes).make());
                break;
            case FALSE_POSITIVE: // variant absent in MC3 and unfiltered in M2.  Take M2 site and alleles and add M2 center INFO field
                final VariantContext m2 = truthVersusEval.getEval();
                vcfWriter.add(new VariantContextBuilder(m2.getSource(), m2.getContig(), m2.getStart(), m2.getEnd(), m2.getAlleles())
                        .attribute(CENTERS_KEY, M2_CENTER_NAME).genotypes(genotypes).make());
                break;
            case FALSE_NEGATIVE:    // variant present in MC3 and absent in M2.  Emit unchanged except for genotype
                vcfWriter.add(new VariantContextBuilder(truthVersusEval.getTruth()).genotypes(genotypes).make());
                break;
            case FILTERED_TRUE_NEGATIVE:    // absent in MC3, filtered in M2.  Skip.
                break;
            case FILTERED_FALSE_NEGATIVE:   // present in MC3, filtered in M2.  Add M2 center and M2 filters INFO field to MC3.
                final VariantContextBuilder vcb = makeVariantContextBuilderWithM2Center(truthVersusEval.getTruth())
                        .attribute(M2_FILTERS_KEY, truthVersusEval.getEval().getFilters().stream().collect(Collectors.toList()))
                        .genotypes(genotypes);
                vcfWriter.add(vcb.make());
                break;
            default:
                throw new IllegalStateException("Unexpected ConcordanceState: " + concordanceState.toString());
        }

    }

    private VariantContextBuilder makeVariantContextBuilderWithM2Center(final VariantContext mc3) {
        final List<String> centers = mc3.hasAttribute(CENTERS_KEY) ?
                mc3.getAttributeAsStringList(CENTERS_KEY, "") : new ArrayList<>();
        centers.add(M2_CENTER_NAME);
        return new VariantContextBuilder(mc3).attribute(CENTERS_KEY, centers);
    }

    @Override
    public void closeTool(){
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }
}
