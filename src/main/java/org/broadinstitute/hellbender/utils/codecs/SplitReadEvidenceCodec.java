package org.broadinstitute.hellbender.utils.codecs;

import com.google.common.base.Splitter;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.index.tabix.TabixFormat;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.tools.sv.SplitReadEvidence;
import org.broadinstitute.hellbender.tools.sv.SplitReadEvidenceSortMerger;
import org.broadinstitute.hellbender.utils.io.FeatureOutputStream;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

/** Codec to handle SplitReadEvidence in tab-delimited text files */
public class SplitReadEvidenceCodec extends AbstractTextCodec<SplitReadEvidence>
        implements FeatureOutputCodec<SplitReadEvidence, FeatureOutputStream<SplitReadEvidence>> {

    public static final String FORMAT_SUFFIX = ".sr.txt";
    public static final String COL_DELIMITER = "\t";
    public static final String DIRECTION_RIGHT = "right";
    public static final String DIRECTION_LEFT = "left";
    private static final Splitter splitter = Splitter.on(COL_DELIMITER);

    @Override
    public boolean canDecode( final String path ) {
        String toDecode = path.toLowerCase();
        if ( IOUtil.hasBlockCompressedExtension(toDecode) ) {
            toDecode = toDecode.substring(0, toDecode.lastIndexOf('.'));
        }
        return toDecode.endsWith(FORMAT_SUFFIX);
    }

    @Override
    public Class<SplitReadEvidence> getFeatureType() {
        return SplitReadEvidence.class;
    }

    @Override
    public FeatureCodecHeader readHeader( final Iterator<String> itr ) {
        return FeatureCodecHeader.EMPTY_HEADER;
    }

    @Override
    public SplitReadEvidence decode( final Iterator<String> itr ) {
        return itr.hasNext() ? decode(itr.next()) : null;
    }

    @Override
    public TabixFormat getTabixFormat() {
        return new TabixFormat(TabixFormat.ZERO_BASED, 1, 2, 0, '#', 0);
    }

    public SplitReadEvidence decode( final String line ) {
        final List<String> tokens = splitter.splitToList(line);
        if (tokens.size() != 5) {
            throw new IllegalArgumentException("Invalid number of columns: " + tokens.size());
        }
        final String contig = tokens.get(0);
        final int position = Integer.parseUnsignedInt(tokens.get(1)) + 1; // Adjust for 0-based indexing
        if (!tokens.get(2).equals(DIRECTION_LEFT) && !tokens.get(2).equals(DIRECTION_RIGHT)) {
            throw new IllegalArgumentException("Invalid direction: " + tokens.get(2));
        }
        final boolean strand = tokens.get(2).equals(DIRECTION_RIGHT);
        final int count = Integer.parseUnsignedInt(tokens.get(3));
        final String sample = tokens.get(4);
        return new SplitReadEvidence(sample, contig, position, count, strand);
    }

    @Override
    public FeatureOutputStream<SplitReadEvidence> makeSink( final GATKPath path,
                                                            final SAMSequenceDictionary dict,
                                                            final List<String> sampleNames,
                                                            final int compressionLevel ) {
        return new FeatureOutputStream<>(path,
                                        getTabixFormat(),
                                        SplitReadEvidenceCodec::encode,
                                        dict,
                                        compressionLevel);
    }

    @Override
    public void encode( final SplitReadEvidence ev,
                        final FeatureOutputStream<SplitReadEvidence> os ) {
        os.write(ev);
    }

    @Override
    public FeatureSink<SplitReadEvidence> makeSortMerger( final GATKPath path,
                                                          final SAMSequenceDictionary dict,
                                                          final List<String> sampleNames,
                                                          final int compressionLevel ) {
        return new SplitReadEvidenceSortMerger(dict, makeSink(path, dict, sampleNames, compressionLevel));
    }

    public static String encode(final SplitReadEvidence ev) {
        final List<String> columns = Arrays.asList(
                ev.getContig(),
                Integer.toString(ev.getStart() - 1),
                ev.getStrand() ? SplitReadEvidenceCodec.DIRECTION_RIGHT : SplitReadEvidenceCodec.DIRECTION_LEFT,
                Integer.toString(ev.getCount()),
                ev.getSample()
        );
        return String.join(COL_DELIMITER, columns);
    }
}
