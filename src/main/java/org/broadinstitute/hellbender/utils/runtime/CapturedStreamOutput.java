package org.broadinstitute.hellbender.utils.runtime;

import org.apache.commons.io.IOUtils;
import org.apache.commons.io.output.NullOutputStream;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.HardThresholdingOutputStream;

import java.io.*;
import java.util.EnumMap;

/**
 * Stream output captured from a stream.
 */
public class CapturedStreamOutput extends StreamOutput {
    protected final InputStream processStream;
    private final EnumMap<StreamLocation, OutputStream> outputStreams = new EnumMap<>(StreamLocation.class);

    // Size of buffers used to transfer data read from process streams
    public final static int STREAM_BLOCK_TRANSFER_SIZE = 4096;

    /**
     * The byte stream to capture content or null if no output string content was requested.
     */
    protected final ByteArrayOutputStream bufferStream;

    /**
     * True if the buffer is truncated.
     */
    private boolean bufferTruncated = false;

    /**
     * @param settings       Settings that define what to capture.
     * @param processStream  Stream to capture output.
     * @param standardStream Stream to write debug output.
     */
    public CapturedStreamOutput(OutputStreamSettings settings, InputStream processStream, PrintStream standardStream) {
        this.processStream = processStream;
        int bufferSize = settings.getBufferSize();
        this.bufferStream = (bufferSize < 0) ? new ByteArrayOutputStream() : new ByteArrayOutputStream(bufferSize);

        for (StreamLocation location : settings.getStreamLocations()) {
            OutputStream outputStream;
            switch (location) {
                case Buffer:
                    if (bufferSize < 0) {
                        outputStream = this.bufferStream;
                    } else {
                        outputStream = new HardThresholdingOutputStream(bufferSize) {
                            @Override
                            protected OutputStream getStream() {
                                return bufferTruncated ? NullOutputStream.INSTANCE : bufferStream;
                            }

                            @Override
                            protected void thresholdReached() {
                                bufferTruncated = true;
                            }
                        };
                    }
                    break;
                case File:
                    try {
                        outputStream = new FileOutputStream(settings.getOutputFile(), settings.isAppendFile());
                    } catch (IOException e) {
                        throw new UserException.BadInput(e.getMessage());
                    }
                    break;
                case Standard:
                    outputStream = standardStream;
                    break;
                default:
                    throw new GATKException("Unexpected stream location: " + location);
            }
            this.outputStreams.put(location, outputStream);
        }
    }

    @Override
    public byte[] getBufferBytes() {
        return bufferStream.toByteArray();
    }

    @Override
    public boolean isBufferTruncated() {
        return bufferTruncated;
    }

    /**
     * Drain the input stream to keep the process from backing up until it's empty.
     * File streams will be closed automatically when this method returns.
     *
     * @throws IOException When unable to read or write.
     */
    @SuppressWarnings("deprecation")
    public void read() throws IOException {
        int readCount = 0;
        try {
            // read until eof
            byte[] buf = new byte[STREAM_BLOCK_TRANSFER_SIZE];
            while ((readCount = processStream.read(buf)) >= 0)
                for (OutputStream outputStream : this.outputStreams.values()) {
                    outputStream.write(buf, 0, readCount);
                }
        } finally {
            for (StreamLocation location : this.outputStreams.keySet()) {
                OutputStream outputStream = this.outputStreams.get(location);
                outputStream.flush();
                if (location != StreamLocation.Standard)
                    IOUtils.closeQuietly(outputStream);
            }
        }
    }
}
