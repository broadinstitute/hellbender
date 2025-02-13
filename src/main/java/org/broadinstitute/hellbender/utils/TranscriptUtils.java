package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.util.Locatable;
import javassist.Loader;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorUtils;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Class with utility methods for working with transcripts.
 */
public class TranscriptUtils {

    //==================================================================================================================
    // Public Static Members:

    /**
     * Extracts the transcript sequence from the reference context given the exons.
     * The given exons are assumed to be all on the same contig.
     * The exons will be sorted in the order they appear in the genome before extracting the sequences.
     * @param refContext the reference context.  Must not be {@code null}.
     * @param exons the exons of the transcript.  Must not be {@code null}.
     * @return the transcript sequence as coded in the given reference context.
     */
    public static final String extractTrascriptFromReference(final ReferenceContext refContext,
                                                             final List<Locatable> exons) {
        return extractTrascriptFromReference(refContext, exons, false);
    }

    /**
     * Extracts the transcript sequence from the reference context given the exons.
     * The given exons are assumed to be all on the same contig.
     * The exons will be sorted in the order they appear in the genome before extracting the sequences.
     * @param refContext the reference context.  Must not be {@code null}.
     * @param exons the exons of the transcript.  Must not be {@code null}.
     * @param convertHg19ContigToB37Contig whether to convert the contig names from hg19 to b37.
     * @return the transcript sequence as coded in the given reference context.
     */
    public static final String extractTrascriptFromReference(final ReferenceContext refContext,
                                                             final List<Locatable> exons,
                                                             final boolean convertHg19ContigToB37Contig ) {

        // TODO: THIS CONVERSION IS BAD.  WE MUST REFACTOR THIS EVENTUALLY TO REMOVE THIS AND MAKE SOME B37 DATA SOURCES.
        Utils.nonNull(refContext);
        Utils.nonNull(exons);

        final StringBuilder transcript = new StringBuilder();

        // We should iterate through the list of exons in sorted order so we can simply append them together.
        for (final Locatable exon : exons.stream().sorted(Comparator.comparingInt(Locatable::getStart).thenComparing(Locatable::getEnd)).toList() ) {
            final SimpleInterval exonInterval;

            if ( convertHg19ContigToB37Contig ) {
                exonInterval = new SimpleInterval(FuncotatorUtils.convertHG19ContigToB37Contig(exon.getContig()), exon.getStart(), exon.getEnd());
            }
            else {
                exonInterval = new SimpleInterval(exon.getContig(), exon.getStart(), exon.getEnd());
            }
            transcript.append(new String(refContext.getBases(exonInterval)));
        }
        return transcript.toString();
    }

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Private Members:

    //==================================================================================================================
    // Constructors:

    //==================================================================================================================
    // Override Methods:

    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    //==================================================================================================================
    // Helper Data Types:

}
