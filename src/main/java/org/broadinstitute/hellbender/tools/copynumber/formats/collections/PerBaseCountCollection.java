package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.PerBaseCount;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;

/**
 * Simple data structure to pass and read/write a List of {@link PerBaseCount} objects.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 * @author Robert Klein &lt;rklein@broadinstitute.org&gt;
 */
public final class PerBaseCountCollection extends AbstractSampleLocatableCollection<PerBaseCount> {
    //note to developers: repeat the column headers in Javadoc so that they are viewable when linked
    /**
     * CONTIG, POSITION, A, C, G, T, N
     */
    enum PerBaseCountTableColumn {
        CONTIG,
        POSITION,
        A,
        C,
        G,
        T,
        N;

        static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

    private static final Function<DataLine, PerBaseCount> PER_BASE_COUNT_RECORD_FROM_DATA_LINE_DECODER = dataLine -> {
        final String contig = dataLine.get(PerBaseCountTableColumn.CONTIG);
        final int position = dataLine.getInt(PerBaseCountTableColumn.POSITION);
        final int a = dataLine.getInt(PerBaseCountTableColumn.A);
        final int c = dataLine.getInt(PerBaseCountTableColumn.C);
        final int g = dataLine.getInt(PerBaseCountTableColumn.G);
        final int t = dataLine.getInt(PerBaseCountTableColumn.T);
        final int n = dataLine.getInt(PerBaseCountTableColumn.N);
        final HashMap<Nucleotide, Integer> perBaseCount = new HashMap<>();
        perBaseCount.put(Nucleotide.A, a);
        perBaseCount.put(Nucleotide.C, c);
        perBaseCount.put(Nucleotide.G, g);
        perBaseCount.put(Nucleotide.T, t);
        perBaseCount.put(Nucleotide.N, n);
        final SimpleInterval interval = new SimpleInterval(contig, position, position);
        return new PerBaseCount(interval, perBaseCount);
    };

    private static final BiConsumer<PerBaseCount, DataLine> PER_BASE_COUNT_RECORD_TO_DATA_LINE_ENCODER = (allelicCount, dataLine) ->
            dataLine.append(allelicCount.getInterval().getContig())
                    .append(allelicCount.getInterval().getEnd())
                    .append(allelicCount.getPerBaseCount().get(Nucleotide.A))
                    .append(allelicCount.getPerBaseCount().get(Nucleotide.C))
                    .append(allelicCount.getPerBaseCount().get(Nucleotide.G))
                    .append(allelicCount.getPerBaseCount().get(Nucleotide.T))
                    .append(allelicCount.getPerBaseCount().get(Nucleotide.N));

    public PerBaseCountCollection(final File inputFile) {
        super(inputFile, PerBaseCountCollection.PerBaseCountTableColumn.COLUMNS, PER_BASE_COUNT_RECORD_FROM_DATA_LINE_DECODER, PER_BASE_COUNT_RECORD_TO_DATA_LINE_ENCODER);
    }

    public PerBaseCountCollection(final SampleLocatableMetadata metadata,
                                  final List<PerBaseCount> PerBaseCounts) {
        super(metadata, PerBaseCounts, PerBaseCountCollection.PerBaseCountTableColumn.COLUMNS, PER_BASE_COUNT_RECORD_FROM_DATA_LINE_DECODER, PER_BASE_COUNT_RECORD_TO_DATA_LINE_ENCODER);
    }
}

