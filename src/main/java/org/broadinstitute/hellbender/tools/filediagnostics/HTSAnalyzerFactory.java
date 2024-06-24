package org.broadinstitute.hellbender.tools.filediagnostics;

import htsjdk.samtools.util.FileExtensions;
import org.broadinstitute.hellbender.engine.GATKPath;

import java.io.File;

/**
 * Class for creating an analyzer based on an alignment file type.
 */
public class HTSAnalyzerFactory {

    public static HTSAnalyzer getFileAnalyzer(final GATKPath inputPath, final File outputFile, final long countLimit) {
        System.out.println(inputPath.getRawInputString());
        if (inputPath.isCram()) {
            return new CRAMAnalyzer(inputPath, outputFile, countLimit);
        } else if (inputPath.hasExtension(FileExtensions.CRAM_INDEX)) {
            return new CRAIAnalyzer(inputPath, outputFile);
        } else if (inputPath.hasExtension(FileExtensions.BAI_INDEX)) {
            return new BAIAnalyzer(inputPath, outputFile);
        } else {
            throw new RuntimeException("Unsupported diagnostic file type: " + inputPath.getRawInputString());
        }
    }
}
