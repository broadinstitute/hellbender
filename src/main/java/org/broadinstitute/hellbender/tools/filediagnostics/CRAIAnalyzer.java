package org.broadinstitute.hellbender.tools.filediagnostics;

import htsjdk.samtools.CRAMCRAIIndexer;
import htsjdk.samtools.cram.CRAIIndex;
import htsjdk.samtools.util.RuntimeIOException;
import org.broadinstitute.hellbender.engine.GATKPath;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;

/**
 * Analyzer for CRAM (.crai) index files.
 */
public class CRAIAnalyzer extends HTSAnalyzer {

    final FileOutputStream fos;

    public CRAIAnalyzer(final GATKPath inputPath, final File outputFile) {
        super(inputPath, outputFile);
        try {
            fos = new FileOutputStream(outputFile);
        } catch (final IOException e) {
            throw new RuntimeIOException(e);
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
        try (final InputStream is = inputPath.getInputStream()) {
            final CRAIIndex craiIndex = CRAMCRAIIndexer.readIndex(is);
            emitln("\nSeqId AlignmentStart AlignmentSpan ContainerOffset SliceOffset SliceSize\n");
            craiIndex.getCRAIEntries().stream().forEach(e -> emitln(e.toString()));
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public void close() throws IOException {
        if (fos != null) {
            fos.close();
        }
    }

}

