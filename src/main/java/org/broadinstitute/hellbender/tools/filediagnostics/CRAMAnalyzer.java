package org.broadinstitute.hellbender.tools.filediagnostics;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.cram.build.CramIO;
import htsjdk.samtools.cram.encoding.CRAMEncoding;
import htsjdk.samtools.cram.encoding.EncodingFactory;
import htsjdk.samtools.cram.io.ITF8;
import htsjdk.samtools.cram.structure.*;
import htsjdk.samtools.seekablestream.SeekablePathStream;
import htsjdk.samtools.util.RuntimeIOException;
import org.broadinstitute.hellbender.engine.GATKPath;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Base64;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Analyzer for CRAM files. Displays metadata for each CRAM container,
 * slice, and block.
 *
 * Note: the analyzer does not require a reference for the since it only
 * enumerates metadata and doesn't attempt to dereference the reads.
 */
public class CRAMAnalyzer extends HTSAnalyzer {

    final Map<DataSeries, Long> externalDataSeriesDataSizes = new LinkedHashMap<>();
    final Map<Integer, Long> externalTagDataSizes = new LinkedHashMap<>();
    long coreBlocksDataSize = 0L;
    long recordCount = 0;
    final long countLimit;
    final FileOutputStream fos;

    public CRAMAnalyzer(final GATKPath inputPathName, final File outputFile, final long countLimit) {
        super(inputPathName, outputFile);
        this.countLimit = countLimit;
        try {
            fos = new FileOutputStream(outputFile);
        } catch (final IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    @Override
    public void close() throws IOException {
        if (fos != null) {
            fos.close();
        }
    }

    protected void emitln(final String s) {
        try {
            fos.write(s.getBytes());
            fos.write('\n');
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Run the analyzer for the file.
     */
    protected void doAnalysis() {
        int containerCount = 0;
        try (final SeekablePathStream seekableStream = new SeekablePathStream(this.inputPath.toPath())) {
            final CramHeader cramHeader = analyzeCRAMHeader(seekableStream);
            boolean isEOF = false;
            while (!isEOF && containerCount < countLimit) {
                long containerOffset = seekableStream.position();
                final Container container = new Container(
                        cramHeader.getCRAMVersion(),
                        seekableStream,
                        containerOffset);
                 isEOF = analyzeContainer(container, ++containerCount) ;
            }
        }
        catch (IOException e) {
            throw new RuntimeIOException(e);
        }

        emitln("\nTotal Record Count: " + recordCount);
        emitDataDistribution();
    }

    /**
     * Display metadata for a CRAM file header.
     *
     */
    public CramHeader analyzeCRAMHeader(InputStream is) {
        final CramHeader cramHeader = CramIO.readCramHeader(is);
        emitln("\nCRAM File: " + inputPath);
        emitln("CRAM Version: " + cramHeader.getCRAMVersion().toString());
        emitln("CRAM ID Contents: " + String.format("%s", Base64.getEncoder().encodeToString(cramHeader.getId())));

        final SAMFileHeader samHeader = Container.readSAMFileHeaderContainer(cramHeader.getCRAMVersion(), is, inputPath.getRawInputString());
        emitln("\n" + samHeader.toString());
        final SAMSequenceDictionary dict = samHeader.getSequenceDictionary();
        emitln(dict.toString());
        dict.getSequences().forEach(e -> emitln(e.toString()));
        return cramHeader;
    }

    /**
     * Display metadata for a CRAM file container.
     * return true if container is EOF container
     */
    public boolean analyzeContainer(Container container, int containerCount) {
        final ContainerHeader containerHeader = container.getContainerHeader();
        emitln(String.format(
                "\n***Container #:%d %s byteOffset=%d",
                containerCount, containerHeader.toString(), container.getContainerByteOffset()));
        if (container.isEOF()) {
            return true;
        }
        // contentID to DataSeries (excludes tags ?)
        final Map<Integer, DataSeries> dataSeriesByContentID = analyzeCompressionHeader(container.getCompressionHeader());
        int sliceCount = 0;
        for (final Slice slice : container.getSlices()) {
            analyzeSlice(slice, ++sliceCount, dataSeriesByContentID);
        }
        return false;
    }

    public Map<Integer, DataSeries> analyzeCompressionHeader(final CompressionHeader compressionHeader) {
        //preservation map, data series encoding map, and tag encoding map
        analyzePreservationMap(compressionHeader);
        final Map<Integer, DataSeries> dataSeriesByContentID = analyzeDataSeriesEncodingMap(compressionHeader);
        analyzeTagEncodingMap(compressionHeader);
        return dataSeriesByContentID;
    }

    public void analyzePreservationMap(final CompressionHeader compressionHeader) {
        emitln(String.format(
                "Requires reference (%b); Preserved read names (%b); APDelta (%b)",
                    compressionHeader.isReferenceRequired(),
                    compressionHeader.isPreserveReadNames(),
                    compressionHeader.isAPDelta()));
    }

    public Map<Integer, DataSeries> analyzeDataSeriesEncodingMap(final CompressionHeader compressionHeader) {
        // Since each container can in theory use a different set of ids for the data series in that
        // container, we need to reconstruct the mapping for each container, and return it so that as
        // we process Slice blocks and track data, we can calculate the data series to which we should
        // accrue that block's data
        final Map<Integer, DataSeries> dataSeriesByContentID = new LinkedHashMap<>();

        emitln("\nData Series Encodings:\n");
        final CompressionHeaderEncodingMap encodingMap = compressionHeader.getEncodingMap();
        for (final DataSeries dataSeries: DataSeries.values()) {
            final EncodingDescriptor encodingDescriptor = encodingMap.getEncodingDescriptorForDataSeries(dataSeries);
            if (encodingDescriptor == null) {
                emitln(String.format("%-50s not present",
                        String.format("DataSeries (%s/%s)",
                                dataSeries.getCanonicalName(),
                                dataSeries.name())));
            } else {
                // the encoding map has an encoding descriptor for this data series; determine if its external
                // and if so, update the contentID to data series map so we can track how much data is used for
                // the blocks for this data series
                final Integer externalContentID = contentIDFromExternalDescriptor(dataSeries.getType(), encodingDescriptor);
                if (externalContentID != null) {
                    dataSeriesByContentID.put(externalContentID,  dataSeries);
                }
                emitln(String.format("%-50s %s",
                        String.format("DataSeries (%s/%s)",
                                dataSeries.getCanonicalName(),
                                dataSeries.name()),
                        encodingDescriptorAsString(
                                dataSeries.getType(),
                                encodingDescriptor)));
            }
        }
        return dataSeriesByContentID;
    }

    public void analyzeTagEncodingMap(final CompressionHeader compressionHeader) {
        emitln("\nTag Encodings:");
        for (final Map.Entry<Integer, EncodingDescriptor> entry : compressionHeader.getTagEncodingMap().entrySet()) {
            final Integer contentID = entry.getKey(); // is this content ID ?
            final EncodingDescriptor ep = entry.getValue();
            emitln(String.format("%-50s %s",
                    String.format("Content ID/Tag (%s/%s)",
                            contentID,
                            decomposeTagNameAndType(contentID)),
                     encodingDescriptorAsString(DataSeriesType.BYTE_ARRAY, ep)));
       }
    }

    /**
     * Display metadata for a CRAM container slice.
     *
     */
    public void analyzeSlice(final Slice slice, final int sliceCount, final Map<Integer, DataSeries> dataSeriesByContentID) {
        emitln(String.format("\n******Slice #: %d %s",
                sliceCount,
                slice.toString()));
        emitln(String.format("%-50s %s",
                "Header block ",
                slice.getSliceHeaderBlock()));
        emitln(String.format("%-50s %s",
                "Core block ",
                slice.getSliceBlocks().getCoreBlock()));

        if (slice.getEmbeddedReferenceContentID() != Slice.EMBEDDED_REFERENCE_ABSENT_CONTENT_ID) {
            emitln(String.format("Embedded reference block ID %d", slice.getEmbeddedReferenceContentID()));
        }

        slice.getSliceBlocks().getExternalContentIDs().forEach((id) -> {
            final String blockID = dataSeriesByContentID.get(id) == null ?
                    Integer.toString(id) : // not a fixed data series (i.e., tag block)
                    dataSeriesByContentID.get(id).getCanonicalName();
            emitln(String.format("%-50s %s",
                    String.format("External Block (%s):", blockID),
                    slice.getSliceBlocks().getExternalBlock(id).toString()));
        });

        updateDataDistribution(slice, dataSeriesByContentID);
        recordCount += slice.getNumberOfRecords();
    }

    final void updateDataDistribution(final Slice slice, final Map<Integer, DataSeries> dataSeriesByContentID) {
        final SliceBlocks sliceBlocks = slice.getSliceBlocks();

        coreBlocksDataSize += sliceBlocks.getCoreBlock().getCompressedContentSize();
        final Map<Integer, EncodingDescriptor> tagContentIDs = slice.getCompressionHeader().getTagEncodingMap();
        for (final Integer contentID : sliceBlocks.getExternalContentIDs()) {
            //if (tagContentIDs.containsKey(contentID)) {
            final DataSeries ds = dataSeriesByContentID.get(contentID);
            if (ds == null) {
                // accrue to tag data
                externalTagDataSizes.merge(
                        contentID,
                        Long.valueOf(slice.getSliceBlocks().getExternalBlock(contentID).getCompressedContentSize()),
                        (oldValue, increment) -> oldValue + increment);
            } else {
                // accrue to fixed DataSeries ID
                externalDataSeriesDataSizes.merge(
                        ds,
                        Long.valueOf(sliceBlocks.getExternalBlock(contentID).getCompressedContentSize()),
                        (oldValue, increment) -> oldValue + increment);
            }
        }
    }

    public void emitDataDistribution() {
        emitln("\nCore Block(s) Total: " + String.format("%,d\n", coreBlocksDataSize));
        emitln("External Data Series Totals (external block resolution - all core data encodings accrue to the core block):\n");
        for (final DataSeries ds : externalDataSeriesDataSizes.keySet()) {
            emitln(String.format("%s: %,d", ds, externalDataSeriesDataSizes.get(ds)));
        }

        emitln("\nTag Series Distribution:\n");
        for (final Map.Entry<Integer, Long> externalEntry : externalTagDataSizes.entrySet()) {
            final Integer contentID = externalEntry.getKey();
            final String tagName = String.format("%d (%s)", contentID, decomposeTagNameAndType(externalEntry.getKey()));
            emitln(String.format("%s: %,d", tagName, externalEntry.getValue()));
        }
    }

    private String encodingDescriptorAsString(final DataSeriesType dsType, final EncodingDescriptor descriptor) {
        final String encodingIDString = EncodingID.values()[descriptor.getEncodingID().getId()].toString();
        final CRAMEncoding<?> encoding = EncodingFactory.createCRAMEncoding(
                dsType, descriptor.getEncodingID(), descriptor.getEncodingParameters());
        return String.format("%s (%s)", encodingIDString, encoding.toString());
    }

    private Integer contentIDFromExternalDescriptor(final DataSeriesType dsType, final EncodingDescriptor descriptor) {
        final CRAMEncoding<?> encoding = EncodingFactory.createCRAMEncoding(
                dsType, descriptor.getEncodingID(), descriptor.getEncodingParameters());
        switch (descriptor.getEncodingID()) {
            case EXTERNAL:
                return ITF8.readUnsignedITF8(encoding.toSerializedEncodingParams());

            case BYTE_ARRAY_STOP:
                final ByteBuffer buf = ByteBuffer.wrap(encoding.toSerializedEncodingParams());
                buf.order(ByteOrder.LITTLE_ENDIAN);
                final byte stopByte = buf.get(); // discard it
                return ITF8.readUnsignedITF8(buf);

            // Everything else is either not external, or is hybrid (BYTE_ARRAY_LEN), so we can't
            // track the data for these at the block level since the data foes either partially or
            // fully into the core block
            case BYTE_ARRAY_LEN:
            case NULL:
            case GOLOMB:
            case HUFFMAN:
            case BETA:
            case SUBEXPONENTIAL:
            case GOLOMB_RICE:
            case GAMMA:
            default:
                return null;
        }
    }

    private String decomposeTagNameAndType(final int contentID) {
        return ReadTag.intToNameType4Bytes(contentID);
    }

}

