package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * Local reference context at a variant position.
 *
 * </p>The annotation gives ten reference bases each to the left and right of the variant start and the start base for a total of 21 reference bases.
 * Start position is defined as one base before indels.  For example, the reference context AAAAAAAAAACTTTTTTTTTT would apply to a SNV variant
 * context with ref allele C and alt allele G as well as to a deletion variant context with ref allele CT and alt allele C.</p>
 *
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Annotate with local reference bases (REF_BASES)")
public class ReferenceBases implements InfoFieldAnnotation {

    private static final int NUM_BASES_ON_EITHER_SIDE = 10;
    private static final int REFERENCE_CONTEXT_LENGTH = 2*NUM_BASES_ON_EITHER_SIDE + 1;

    protected final OneShotLogger warning = new OneShotLogger(this.getClass());

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(GATKVCFConstants.REFERENCE_BASES_KEY); }

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        if (ref==null)  {
            warning.warn("REF_BASES requires the reference to annotate, none was provided");
            return Collections.emptyMap();
        }

        final String bases = annotate(ref, vc);
        return Collections.singletonMap(GATKVCFConstants.REFERENCE_BASES_KEY, bases );

    }

    public static String annotate(final ReferenceContext ref, final VariantContext vc) {
        Utils.nonNull(ref);
        final int basesToDiscardInFront = Math.max(vc.getStart() - ref.getWindow().getStart() - NUM_BASES_ON_EITHER_SIDE, 0);
        final String allBases = new String(ref.getBases());
        final int endIndex = Math.min(basesToDiscardInFront + 2 * NUM_BASES_ON_EITHER_SIDE + 1, allBases.length());
        String localBases = allBases.substring(basesToDiscardInFront, endIndex);
        if (localBases.length() < REFERENCE_CONTEXT_LENGTH) {
            localBases = String.join("", localBases, StringUtils.repeat("N", REFERENCE_CONTEXT_LENGTH - localBases.length()));
        }

        return localBases;
    }

    public static String getNMiddleBases(final String bases, final int n){
        Utils.validateArg(bases.length() >= n, "bases must have n or more bases. bases = " + bases);
        Utils.validateArg( bases.length() % 2 == 1, "the length of bases must be an odd number");
        Utils.validateArg( n % 2 == 1, "n must be odd");

        final int numBasesOnEachSide = n/2;
        final int middleIndex = bases.length()/2;
        return bases.substring(middleIndex - numBasesOnEachSide, middleIndex + numBasesOnEachSide + 1);

    }
}
