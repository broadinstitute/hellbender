package org.broadinstitute.hellbender.utils.codecs;

import com.google.common.base.Splitter;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.index.tabix.TabixFormat;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.sv.DepthEvidence;
import org.broadinstitute.hellbender.tools.sv.DepthEvidenceSortMerger;
import org.broadinstitute.hellbender.tools.sv.SVFeaturesHeader;
import org.broadinstitute.hellbender.utils.io.FeatureOutputStream;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/** Codec to handle DepthEvidence in tab-delimited text files */
public class DepthEvidenceCodec extends AbstractTextCodec<DepthEvidence>
        implements FeatureOutputCodec<DepthEvidence, FeatureOutputStream<DepthEvidence>> {

    public static final String FORMAT_SUFFIX = ".rd.txt";
    public static final String COL_DELIMITER = "\t";
    public static final String COMMENT_LINE_START = "#";
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
    public Class<DepthEvidence> getFeatureType() {
        return DepthEvidence.class;
    }

    @Override
    public FeatureCodecHeader readHeader( final Iterator<String> itr ) {
        if ( itr.hasNext() ) {
            final String nextLine = itr.next();
            if ( nextLine.startsWith(COMMENT_LINE_START) ) {
                final List<String> headerCols = splitter.splitToList(nextLine);
                final SVFeaturesHeader svFeaturesHeader =
                        new SVFeaturesHeader(DepthEvidence.class.getSimpleName(),
                                            "unknown",
                                            null,
                                            headerCols.subList(3, headerCols.size()));
                return new FeatureCodecHeader(svFeaturesHeader, FeatureCodecHeader.NO_HEADER_END);
            }
        }
        throw new UserException.BadInput("Depth evidence file did not have a header line");
    }

    @Override
    public DepthEvidence decode( final Iterator<String> itr ) {
        while ( itr.hasNext() ) {
            final String nextLine = itr.next();
            if ( !nextLine.startsWith(COMMENT_LINE_START) ) {
                return decode(nextLine);
            }
        }
        return null;
    }

    @Override
    public TabixFormat getTabixFormat() {
        return new TabixFormat(TabixFormat.ZERO_BASED, 1, 2, 3, COMMENT_LINE_START.charAt(0), 1);
    }

    public DepthEvidence decode( final String line ) {
        final List<String> tokens = splitter.splitToList(line);
        if (tokens.size() < 3) {
            throw new IllegalArgumentException("Expected at least 3 columns but found " + tokens.size());
        }
        final String contig = tokens.get(0);
        final int start = Integer.parseUnsignedInt(tokens.get(1)) + 1; // Adjust for 0-based indexing
        final int end = Integer.parseUnsignedInt(tokens.get(2)); // Adjust for 0-based indexing and inclusive intervals
        final int numCounts = tokens.size() - 3;
        final int[] counts = new int[numCounts];
        for (int i = 3; i < tokens.size(); i++) {
            counts[i - 3] = Integer.parseUnsignedInt(tokens.get(i));
        }
        return new DepthEvidence(contig, start, end, counts);
    }

    @Override
    public FeatureOutputStream<DepthEvidence> makeSink( final GATKPath path,
                                                        final SAMSequenceDictionary dict,
                                                        final List<String> sampleNames,
                                                        final int compressionLevel ) {
        final FeatureOutputStream<DepthEvidence> foStream =
                new FeatureOutputStream<>(path,
                                        getTabixFormat(),
                                        DepthEvidenceCodec::encode,
                                        dict,
                                        compressionLevel);
        final StringBuilder sb = new StringBuilder("#Chr\tStart\tEnd");
        for ( final String sampleName : sampleNames ) {
            sb.append('\t').append(sampleName);
        }
        foStream.writeHeader(sb.toString());
        return foStream;
    }

    @Override
    public void encode( final DepthEvidence ev, final FeatureOutputStream<DepthEvidence> os ) {
        os.write(ev);
    }

    @Override
    public FeatureSink<DepthEvidence> makeSortMerger( final GATKPath path,
                                                              final SAMSequenceDictionary dict,
                                                              final List<String> sampleNames,
                                                              final int compressionLevel ) {
        return new DepthEvidenceSortMerger(dict, makeSink(path, dict, sampleNames, compressionLevel));
    }

    public static String encode(final DepthEvidence ev) {
        final int[] counts = ev.getCounts();
        final int numCounts = counts.length;
        final List<String> columns = new ArrayList<>(3 + numCounts);
        columns.add(ev.getContig());
        columns.add(Integer.toString(ev.getStart() - 1));
        columns.add(Integer.toString(ev.getEnd()));
        for ( final int count : counts ) {
            columns.add(Integer.toString(count));
        }
        return String.join(COL_DELIMITER, columns);
    }
}
