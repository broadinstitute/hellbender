package org.broadinstitute.hellbender.utils.codecs;

import com.google.common.base.Splitter;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.index.tabix.TabixFormat;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.tools.sv.DiscordantPairEvidence;
import org.broadinstitute.hellbender.utils.io.FeatureOutputStream;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

/** Codec to handle DiscordantPairEvidence in tab-delimited text files */
public class DiscordantPairEvidenceCodec extends AbstractTextCodec<DiscordantPairEvidence>
        implements FeatureOutputCodec<DiscordantPairEvidence, FeatureOutputStream<DiscordantPairEvidence>> {

    public static final String FORMAT_SUFFIX = ".pe.txt";
    public static final String COL_DELIMITER = "\t";
    private static final Splitter splitter = Splitter.on(COL_DELIMITER);
    public static final String STRAND_PLUS = "+";
    public static final String STRAND_MINUS = "-";

    @Override
    public boolean canDecode( final String path ) {
        String toDecode = path.toLowerCase();
        if ( IOUtil.hasBlockCompressedExtension(toDecode) ) {
            toDecode = toDecode.substring(0, toDecode.lastIndexOf('.'));
        }
        return toDecode.endsWith(FORMAT_SUFFIX);
    }

    @Override
    public Class<DiscordantPairEvidence> getFeatureType() {
        return DiscordantPairEvidence.class;
    }

    @Override
    public FeatureCodecHeader readHeader( final Iterator<String> itr ) {
        return FeatureCodecHeader.EMPTY_HEADER;
    }

    @Override
    public DiscordantPairEvidence decode( final Iterator<String> itr ) {
        return itr.hasNext() ? decode(itr.next()) : null;
    }

    @Override
    public TabixFormat getTabixFormat() {
        return new TabixFormat(TabixFormat.ZERO_BASED, 1, 2, 0, '#', 0);
    }

    public DiscordantPairEvidence decode( final String line ) {
        final List<String> tokens = splitter.splitToList(line);
        if (tokens.size() != 7) {
            throw new IllegalArgumentException("Invalid number of columns: " + tokens.size());
        }
        final String startContig = tokens.get(0);
        final int start = Integer.parseUnsignedInt(tokens.get(1)) + 1; // Adjust for 0-based indexing
        final boolean startStrand = tokens.get(2).equals(STRAND_PLUS);
        final String endContig = tokens.get(3);
        final int end = Integer.parseUnsignedInt(tokens.get(4)) + 1; // Adjust for 0-based indexing
        final boolean endStrand = tokens.get(5).equals(STRAND_PLUS);
        final String sample = tokens.get(6);
        return new DiscordantPairEvidence(sample, startContig, start, startStrand, endContig, end, endStrand);
    }

    @Override
    public FeatureOutputStream<DiscordantPairEvidence> makeSink( final GATKPath path,
                                                                 final SAMSequenceDictionary dict,
                                                                 final List<String> sampleNames,
                                                                 final int compressionLevel ) {
        return new FeatureOutputStream<>(path,
                                        getTabixFormat(),
                                        DiscordantPairEvidenceCodec::encode,
                                        dict,
                                        compressionLevel);
    }

    @Override
    public void encode( final DiscordantPairEvidence ev,
                        final FeatureOutputStream<DiscordantPairEvidence> os ) {
        os.write(ev);
    }

    @Override
    public FeatureSink<DiscordantPairEvidence> makeSortMerger( final GATKPath path,
                                                               final SAMSequenceDictionary dict,
                                                               final List<String> sampleNames,
                                                               final int compressionLevel ) {
        return makeSink(path, dict, sampleNames, compressionLevel);
    }

    public static String encode(final DiscordantPairEvidence ev) {
        final List<String> columns = Arrays.asList(
                ev.getContig(),
                Integer.toString(ev.getStart() - 1),
                ev.getStartStrand() ? STRAND_PLUS : STRAND_MINUS,
                ev.getEndContig(),
                Integer.toString(ev.getEndPosition() - 1),
                ev.getEndStrand() ? STRAND_PLUS : STRAND_MINUS,
                ev.getSample()
        );
        return String.join(COL_DELIMITER, columns);
    }
}
