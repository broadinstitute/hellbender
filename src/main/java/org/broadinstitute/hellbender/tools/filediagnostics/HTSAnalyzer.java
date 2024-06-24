package org.broadinstitute.hellbender.tools.filediagnostics;

import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.engine.GATKPath;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;

/**
 * Base class for alignment file analyzers.
 */
public abstract class HTSAnalyzer implements Closeable {

    protected GATKPath inputPath;
    protected File outputFile;

    public HTSAnalyzer(final GATKPath filePath, final File outputFile) {
        this.inputPath = filePath;
        this.outputFile = outputFile;
    }

    /**
     * Run the analyzer for the file specified by fileName.
     */
    public void analyze() {
        // set and then reset the global log level
        Log.setGlobalLogLevel(Log.LogLevel.ERROR);
        doAnalysis();
        emitln("");
        Log.setGlobalLogLevel(Log.LogLevel.DEBUG);
    }

    /**
     * Emit a string followed by a newline to the output destination.
     *
     * @param s
     */
    protected void emitln(final String s) {
        System.out.println(s);
    }

    protected abstract void doAnalysis();

    public abstract void close() throws IOException;
}
