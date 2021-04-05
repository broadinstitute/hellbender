package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.primitives.Ints;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.List;
import java.util.OptionalInt;
import java.util.stream.Collectors;

import static java.lang.Math.min;

/**
 * Median distance of variant starts from ends of reads supporting each alt allele.
 *
 * </p>The output is an array containing, for each alt allele, the median distance of the variant start from the closest read end over all reads that best match that allele.</p>
 * </p>For example, for variant context with ref allele A and alt allele C the read position for alt-supporting read GGGGCTT is 2 because the A to C
 * substitution is 2 bases from the right end of the read, which is less than its distance from the left end.
 * For variant context with ref allele AG and alt allele A (deletion) the read position of alt-supporting read ATTTTT is 0.
 * For variant context with ref allele A and alt allele AG (insertion) the read position of alt-supporting read TTTTAG is 1.</p>
 * <p>The annotation considers only the read's bases themselves and not the position they map to with respect to the reference.  For example,
 * suppose a substitution is preceded by 80 matching bases and followed by 10 matching bases, a 10-base deletion, and 10 more matching bases.  Its distance from the end of the read
 * is 20 bases, not 30 bases, because the deleted bases belong to the reference, not the read.  Similarly soft-clipped bases are counted in the distance.</p>
 * <p>This annotation is useful for filtering alignment artifacts.</p>
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Median distance of variant starts from ends of reads supporting each allele (MPOS)")
public class MaxStartPositionDifference extends PerAlleleAnnotation implements StandardMutectAnnotation {

    // we don't want a GGA mode allele with no reads to prejudice us against a site so we assign a non-suspicious value
    private static final int VALUE_FOR_NO_READS = 50;

    @Override
    protected int aggregate(final List<Integer> values) {
        List<Integer> positiveValues = values.stream()
                .filter(v -> v != 0)
                .collect(Collectors.toList());

        if (positiveValues.isEmpty()) {
            return 0;
        }

        System.out.println("Max start positions " + positiveValues);

        int positiveMax = MathUtils.arrayMax(Ints.toArray(positiveValues)) - MathUtils.arrayMin(Ints.toArray(positiveValues));

        return values.isEmpty() ? VALUE_FOR_NO_READS : positiveMax;
    }

    @Override
    protected String getVcfKey() { return GATKVCFConstants.READ_START_POSITION_MAX_DIFF_KEY; }

    @Override
    protected String getDescription() { return "Maximum difference in start position"; }

    @Override
    protected OptionalInt getValueForRead(final GATKRead read, final VariantContext vc) {
        if (vc.getStart() < read.getStart() || read.getEnd() < vc.getStart()) {
            return OptionalInt.empty();
        }

        OptionalInt valueAsInt = OptionalInt.empty();

        // Assume innies.  If they aren't, then we won't return a value
        if(read.isReverseStrand() && !read.mateIsReverseStrand()) {
            valueAsInt = OptionalInt.of(read.getEnd());
        }

        return valueAsInt.isPresent() ? valueAsInt : OptionalInt.empty();
    }
}