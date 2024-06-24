package org.broadinstitute.hellbender.testutils;

import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.logging.BunnyLog;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * Utility interface for CommandLine Program testing. API users that have their own Main implementation
 * should override {@link #makeCommandLineArgs(List)} to make their tests work with {@link IntegrationTestSpec}.
 */
public interface CommandLineProgramTester {

    /**
     * Returns the name for the tested tool.
     */
    public String getTestedToolName();

    /**
     * For testing support. Given a name of a Main CommandLineProgram and it's arguments, builds the arguments appropriate for calling the
     * program through Main.
     *
     * Default behaviour uses {@link #makeCommandLineArgs(List, String)} with the tool name provided by {@link #getTestedToolName()}.
     *
     * @param args List of command line arguments
     * @return String[] of command line arguments
     */
    default String[] makeCommandLineArgs(final List<String> args) {
        return makeCommandLineArgs(args, getTestedToolName());
    }

    /**
     * For testing support. Given a name of a Main CommandLineProgram and it's arguments, builds the arguments appropriate for calling the
     * program through Main.
     *
     * Default behaviour generates a command line in the form "toolname args", with the verbosity parameter returned by {@link #injectDefaultVerbosity(List)}.
     *
     * @param args     List of command line arguments
     * @param toolname name of the tool to test
     * @return String[] of command line arguments
     */
    default String[] makeCommandLineArgs(final List<String> args, final String toolname) {
        List<String> curatedArgs = injectDefaultVerbosity(args);
        final String[] commandLineArgs = new String[curatedArgs.size() + 1];
        commandLineArgs[0] = toolname;
        int i = 1;
        for (final String arg : curatedArgs) {
            commandLineArgs[i++] = arg;
        }
        return commandLineArgs;
    }

    /**
     * Inject the verbosity parameter into the list.
     *
     * Default behaviour look for --verbosity argument; if not found, supply a default value that minimizes the amount of logging output.
     */
    default List<String> injectDefaultVerbosity(final List<String> args) {

        // global toggle for BunnyLog output.
        BunnyLog.setEnabled(false);

        for (String arg : args) {
            if (arg.equalsIgnoreCase("--" + StandardArgumentDefinitions.VERBOSITY_NAME) || arg.equalsIgnoreCase("-" + StandardArgumentDefinitions.VERBOSITY_NAME)) {
                return args;
            }
        }
        List<String> argsWithVerbosity = new ArrayList<>(args);
        argsWithVerbosity.add("--" + StandardArgumentDefinitions.VERBOSITY_NAME);
        argsWithVerbosity.add(Log.LogLevel.ERROR.name());
        return argsWithVerbosity;
    }

    /**
     * Runs the command line implemented by this test.
     *
     * Default behavior uses {@link Main} with the command line arguments created by {@link #makeCommandLineArgs(List)}.
     */
    default Object runCommandLine(final List<String> args) {
        return new Main().instanceMain(makeCommandLineArgs(args));
    }

    /**
     * Lets you explicitly specify a tool to run with the provided arguments
     *
     * Default behavior uses {@link Main} with the command line arguments created by {@link #makeCommandLineArgs(List, String)}.
     */
    default Object runCommandLine(final List<String> args, final String toolName) {
        return new Main().instanceMain(makeCommandLineArgs(args, toolName));
    }

    /**
     * Lets you explicitly specify a tool to run with the provided arguments
     *
     * Default behavior uses {@link Main} with the command line arguments created by {@link #makeCommandLineArgs(List, String)}.
     */
    default Object runCommandLine(final ArgumentsBuilder args, final String toolName) {
        return new Main().instanceMain(makeCommandLineArgs(args.getArgsList(), toolName));
    }

    default Object runCommandLine(final String[] args) {
        return runCommandLine(Arrays.asList(args));
    }

    default Object runCommandLine(final ArgumentsBuilder args) {
        return runCommandLine(args.getArgsList());
    }

}
