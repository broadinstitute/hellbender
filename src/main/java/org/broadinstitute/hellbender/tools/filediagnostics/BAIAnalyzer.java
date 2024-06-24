package org.broadinstitute.hellbender.tools.filediagnostics;

import htsjdk.samtools.BAMIndexer;
import org.broadinstitute.hellbender.engine.GATKPath;

import java.io.File;
import java.io.IOException;

/**
 * Analyzer for BAI files.
 */
public class BAIAnalyzer extends HTSAnalyzer {

    public BAIAnalyzer(final GATKPath inputPath, final File outputFile) {
        super(inputPath, outputFile);
    }

    /**
     * Run the analyzer for the file.
     */
    protected void doAnalysis() {
        System.out.println(String.format("\nOutput written to %s\n", outputFile));
        BAMIndexer.createAndWriteIndex(inputPath.toPath().toFile(), outputFile, true);
    }

    @Override
    public void close() throws IOException {
    }

}

