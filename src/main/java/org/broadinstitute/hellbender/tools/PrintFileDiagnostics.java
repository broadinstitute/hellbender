package org.broadinstitute.hellbender.tools;

import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.tools.filediagnostics.HTSAnalyzer;
import org.broadinstitute.hellbender.tools.filediagnostics.HTSAnalyzerFactory;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.io.File;

/**
 * A diagnostic tool that prints meta information about a GATK input file.
 *
 * Works on files ending in .cram, .crai, and .bai.
 *
 * Sample Usage:
 *
 * gatk PrintFileDiagnostics \
 *   -I input.cram \
 *   -count-limit 10
 */
@ExperimentalFeature
@WorkflowProperties
@CommandLineProgramProperties(
        summary = "Print diagnostic information about a genomics file to stdout",
        oneLineSummary = "Print diagnostic information about a genomics file to stdout",
        programGroup = OtherProgramGroup.class
)
public class PrintFileDiagnostics extends CommandLineProgram {

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc = "Input path for diagnostics",
            optional = false,
            common = true)
    @WorkflowInput
    public GATKPath inputPath;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Outut file for diagnostics (must be a local file)",
            optional = false,
            common = true)
    @WorkflowInput
    public File outputFile;

    @Argument(shortName="count-limit",
            fullName="count-limit",
            doc="Limit on how much output to emit (.cram only)")
    private long countLimit = 1000;

    private HTSAnalyzer htsAnalyzer;

    @Override
    protected void onStartup() {
        super.onStartup();
        htsAnalyzer = HTSAnalyzerFactory.getFileAnalyzer(inputPath, outputFile, countLimit);
    }

    @Override
    protected Object doWork() {
        htsAnalyzer.analyze();
        return 0;
    }

    @Override
    protected void onShutdown() {
        if ( htsAnalyzer != null ) {
            try {
                htsAnalyzer.close();
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
        }
    }

}
