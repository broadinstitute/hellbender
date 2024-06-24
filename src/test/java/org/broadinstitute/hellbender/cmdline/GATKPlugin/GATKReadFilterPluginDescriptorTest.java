package org.broadinstitute.hellbender.cmdline.GATKPlugin;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.TextCigarCodec;
import org.apache.commons.io.output.NullOutputStream;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.testpluggables.TestReadFilter;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.filters.*;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.mockito.internal.util.io.IOUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import javax.validation.constraints.AssertTrue;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.Consumer;
import java.util.stream.Collectors;


public class GATKReadFilterPluginDescriptorTest extends GATKBaseTest {

    // null print stream for the tests
    @SuppressWarnings("deprecation")
    private static final PrintStream nullMessageStream = new PrintStream(new NullOutputStream());

    @DataProvider(name = "defaultFiltersForAllowedValues")
    public Object[][] defaultFiltersForAllowedValues() {
        return new Object[][] {
                {Collections.emptyList()},
                {Collections.singletonList(ReadFilterLibrary.GOOD_CIGAR)},
                {Arrays.asList(ReadFilterLibrary.MAPPED, ReadFilterLibrary.GOOD_CIGAR)},
                {Arrays.asList(ReadFilterLibrary.GOOD_CIGAR, ReadFilterLibrary.MAPPED)}
        };
    }

    @Test(dataProvider = "defaultFiltersForAllowedValues")
    public void testGetAllowedValuesForDescriptorArgument(final List<ReadFilter> defaultFilters) {
        final CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(defaultFilters)),
                Collections.emptySet());
        clp.parseArguments(nullMessageStream, new String[]{});

        final GATKReadFilterPluginDescriptor pluginDescriptor = clp.getPluginDescriptor(GATKReadFilterPluginDescriptor.class);

        // the help for disable-read-filter should point out to the default filters in the order provided
        Assert.assertEquals(pluginDescriptor.getAllowedValuesForDescriptorHelp(ReadFilterArgumentDefinitions.DISABLE_READ_FILTER_LONG_NAME),
                defaultFilters.stream().map(rf -> rf.getClass().getSimpleName()).collect(Collectors.toSet()));

        // test if the help for readFilter is not empty after parsing: if custom validation throws, the help should print the readFilter available
        // the complete set could not checked because that requires to discover all the implemented filters
        Assert.assertFalse(pluginDescriptor.getAllowedValuesForDescriptorHelp(ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME).isEmpty());
    }

    //Filters with arguments to verify filter test method actually filters
    @DataProvider(name="filtersWithFilteringArguments")
    public Object[][] filtersWithGoodArguments(){
        return new Object[][]{
                { LibraryReadFilter.class.getSimpleName(),
                        (Consumer<SetupTest>) this::setupLibraryTest, "--library", "Foo" },
                { MappingQualityReadFilter.class.getSimpleName(),
                        (Consumer<SetupTest>) this::setupMappingQualityTest, "--minimum-mapping-quality", "255" },
                { MappingQualityReadFilter.class.getSimpleName(),
                        (Consumer<SetupTest>) this::setupMappingQualityTest, "--minimum-mapping-quality", "60" },
                { MappingQualityReadFilter.class.getSimpleName(),
                        (Consumer<SetupTest>) this::setupMappingQualityTest, "--maximum-mapping-quality", "40"},
                { ReadGroupReadFilter.class.getSimpleName(),
                        (Consumer<SetupTest>) this::setupReadGroupTest, "--keep-read-group", "fred" },
                { ReadNameReadFilter.class.getSimpleName(),
                        (Consumer<SetupTest>) this::setupReadNameTest, "--read-name", "fred" },
                { SampleReadFilter.class.getSimpleName(),
                        (Consumer<SetupTest>) this::setupSampleTest, "--sample", "fred" }
        };
    }

    @Test(dataProvider = "filtersWithFilteringArguments")
    public void testFilterArguments(
            final String filter,
            final Consumer<SetupTest> setup,
            final String argName,
            final String argValue)
    {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);

        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(null)),
                Collections.emptySet());

        String[] args = {
                "--" + ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME, filter,
                argName, argValue
        };
        clp.parseArguments(nullMessageStream, args);
        ReadFilter rf = instantiateFilter(clp, header);

        // to ensure that the filter is actually working, verify that the test record
        // we're using fails the filter test *before* we set it up to pass the filter
        Assert.assertFalse(rf.test(read));

        // setup the header and read for the test
        setup.accept(new SetupTest(header, read, argValue));

        Assert.assertTrue(rf.test(read));
    }

    @DataProvider(name="filtersWithRequiredArguments")
    public Object[][] filtersWithRequiredArguments(){
        return new Object[][]{
                { LibraryReadFilter.class.getSimpleName(), "--library", "fakeLibrary" },
                { PlatformReadFilter.class.getSimpleName(), "--" + ReadFilterArgumentDefinitions.PL_FILTER_NAME_LONG_NAME, "fakePlatform" },
                { PlatformUnitReadFilter.class.getSimpleName(), "--" + ReadFilterArgumentDefinitions.BLACK_LISTED_LANES_LONG_NAME, "fakeUnit" },
                { ReadGroupReadFilter.class.getSimpleName(), "--" + ReadFilterArgumentDefinitions.KEEP_READ_GROUP_LONG_NAME, "fakeGroup" },
                { ReadGroupBlackListReadFilter.class.getSimpleName(), "--" + ReadFilterArgumentDefinitions.READ_GROUP_BLACK_LIST_LONG_NAME, "tg:sub"},
                { ReadNameReadFilter.class.getSimpleName(), "--" + ReadFilterArgumentDefinitions.READ_NAME_LONG_NAME, "fakeRead" },
                { ReadLengthReadFilter.class.getSimpleName(), "--" + ReadFilterArgumentDefinitions.MAX_READ_LENGTH_ARG_NAME, "10" },
                { SampleReadFilter.class.getSimpleName(), "--sample", "fakeSample" }
        };
    }

    // fail if a filter with required arguments is specified without corresponding arguments
    @Test(dataProvider = "filtersWithRequiredArguments", expectedExceptions = CommandLineException.MissingArgument.class)
    public void testDependentFilterArguments(
            final String filter,
            final String argName,   //unused
            final String argValue)  //unused
    {
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(null)),
                Collections.emptySet());

        String[] args = {
                "--" + ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME, filter  // no args, just enable filters
        };

        clp.parseArguments(nullMessageStream, args);
    }

    // fail if a filter's arguments are passed but the filter itself is not enabled
    @Test(dataProvider = "filtersWithRequiredArguments", expectedExceptions = CommandLineException.class)
    public void testDanglingFilterArguments(
            final String filter,
            final String argName,
            final String argValue)
    {
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(null)),
                Collections.emptySet());

        String[] args = { argName, argValue }; // no read filter set

        // no need to instantiate the filters - dependency errors are caught by the command line parser
        clp.parseArguments(nullMessageStream, args);
    }

    @Test
    public void testNoFiltersSpecified() {
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(null)),
                Collections.emptySet());

        clp.parseArguments(nullMessageStream, new String[]{});

        // get the command line read filters
        final GATKReadFilterPluginDescriptor readFilterPlugin = clp.getPluginDescriptor(GATKReadFilterPluginDescriptor.class);
        final List<ReadFilter> readFilters = readFilterPlugin.getResolvedInstances();
        Assert.assertEquals(readFilters.size(), 0);
    }

    @Test
    public void testEnableMultipleFilters() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);

        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(null)),
                Collections.emptySet());

        String[] args = {
                "--" + ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME, ReadLengthReadFilter.class.getSimpleName(),
                "--" + ReadFilterArgumentDefinitions.MIN_READ_LENGTH_ARG_NAME, "10",
                "--" + ReadFilterArgumentDefinitions.MAX_READ_LENGTH_ARG_NAME, "102",
                "--" + ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME, ReadNameReadFilter.class.getSimpleName(),
                "--" + ReadFilterArgumentDefinitions.READ_NAME_LONG_NAME, "fred"
        };
        clp.parseArguments(nullMessageStream, args);
        ReadFilter rf = instantiateFilter(clp, header);

        // default read name is rejected by the readName filter
        Assert.assertFalse(rf.test(read));

        String readName = read.getName();
        read.setName("fred");
        Assert.assertTrue(rf.test(read)); // accepted

        // trigger the read length filter to reject
        read.setBases(new byte[150]);
        Assert.assertFalse(rf.test(read));
    }

    @Test(expectedExceptions = CommandLineException.class)
    public void testEnableNonExistentFilter() {
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(null)),
                Collections.emptySet());
        clp.parseArguments(nullMessageStream, new String[] {"--RF", "fakeFilter"});
    }

    @DataProvider(name = "duplicateFilters")
    public Object[][] duplicateFilters() {
        return new Object[][] {
                {new String[] {
                        "--RF", ReadLengthReadFilter.class.getSimpleName(),
                        "--RF", ReadLengthReadFilter.class.getSimpleName()}},
                {new String[] {
                        "--RF", ReadFilterLibrary.MAPPED.getClass().getSimpleName(),
                        "--RF", ReadLengthReadFilter.class.getSimpleName(),
                        "--RF", ReadLengthReadFilter.class.getSimpleName()}
                }
        };
    }

    @Test(dataProvider = "duplicateFilters", expectedExceptions = CommandLineException.BadArgumentValue.class)
    public void testEnableDuplicateFilter(final String[] arguments) {
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(null)),
                Collections.emptySet());
        clp.parseArguments(nullMessageStream, arguments);
    }

    @DataProvider(name = "duplicateDisabledFilters")
    public Object[][] duplicateDisabledFilters() {
        return new Object[][] {
                {Collections.singletonList(ReadFilterLibrary.MAPPED),
                        new String[]{
                                "--DF", ReadFilterLibrary.MAPPED.getClass().getSimpleName(),
                                "--DF", ReadFilterLibrary.MAPPED.getClass().getSimpleName()}
                }
        };
    }

    @Test(dataProvider = "duplicateDisabledFilters", expectedExceptions = CommandLineException.BadArgumentValue.class)
    public void testDisableDuplicateFilter(final List<ReadFilter> defaults, final String[] arguments) {
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(defaults)),
                Collections.emptySet());
        clp.parseArguments(nullMessageStream, arguments);
    }

    @Test
    public void testDisableOneFilter() {
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(
                        Collections.singletonList(ReadFilterLibrary.GOOD_CIGAR)
                )),
                Collections.emptySet());
        clp.parseArguments(nullMessageStream, new String[] {
                "--RF", ReadFilterLibrary.MAPPED.getClass().getSimpleName(),
                "--RF", ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS.getClass().getSimpleName(),
                "--" + ReadFilterArgumentDefinitions.DISABLE_READ_FILTER_LONG_NAME, ReadFilterLibrary.GOOD_CIGAR.getClass().getSimpleName()});

        // get the command line read filters
        GATKReadFilterPluginDescriptor readFilterPlugin = clp.getPluginDescriptor(GATKReadFilterPluginDescriptor.class);
        List<ReadFilter> readFilters = readFilterPlugin.getResolvedInstances();

        Assert.assertEquals(readFilters.size(), 2);
        Assert.assertEquals(readFilters.get(0).getClass().getSimpleName(),
                ReadFilterLibrary.MAPPED.getClass().getSimpleName());
        Assert.assertEquals(readFilters.get(1).getClass().getSimpleName(),
                ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS.getClass().getSimpleName());

        Assert.assertEquals(readFilterPlugin.userArgs.getUserDisabledReadFilterNames().size(), 1);
        Assert.assertTrue(readFilterPlugin.userArgs.getUserDisabledReadFilterNames().contains(
                ReadFilterLibrary.GOOD_CIGAR.getClass().getSimpleName()));
        Assert.assertTrue(readFilterPlugin.isDisabledFilter(
                ReadFilterLibrary.GOOD_CIGAR.getClass().getSimpleName()));
    }

    @Test
    public void testDisableNotEnabledFilter() {
        GATKReadFilterPluginDescriptor rfDesc = new GATKReadFilterPluginDescriptor(null);
        final CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(rfDesc),
                Collections.emptySet());
        final String filterName = ReadFilterLibrary.MAPPED.getClass().getSimpleName();
        clp.parseArguments(nullMessageStream, new String[]{
                "--" + ReadFilterArgumentDefinitions.DISABLE_READ_FILTER_LONG_NAME, filterName
        });
        // Make sure mapped filter got disabled with no exception
        Assert.assertTrue(rfDesc.userArgs.getUserDisabledReadFilterNames().contains(filterName));
        Assert.assertTrue(rfDesc.isDisabledFilter(filterName));
    }

    @Test(expectedExceptions = CommandLineException.class)
    public void tesDisableNonExistentReadFilter() {
        final CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(null)),
                Collections.emptySet());
        clp.parseArguments(nullMessageStream, new String[] {
                "--" + ReadFilterArgumentDefinitions.DISABLE_TOOL_DEFAULT_READ_FILTERS, "asdfasdf"
        });
    }

    @Test
    public void testDisableMultipleFilters() {
        List<ReadFilter> defaultFilters = new ArrayList<>();
        defaultFilters.add(ReadFilterLibrary.GOOD_CIGAR);
        defaultFilters.add(ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS);
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(defaultFilters)),
                Collections.emptySet());
        clp.parseArguments(nullMessageStream, new String[] {
                "--RF", ReadFilterLibrary.MAPPED.getClass().getSimpleName(),
                "--" + ReadFilterArgumentDefinitions.DISABLE_READ_FILTER_LONG_NAME, ReadFilterLibrary.GOOD_CIGAR.getClass().getSimpleName(),
                "--" + ReadFilterArgumentDefinitions.DISABLE_READ_FILTER_LONG_NAME, ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS.getClass().getSimpleName()});

        // get the command line read filters
        GATKReadFilterPluginDescriptor readFilterPlugin = clp.getPluginDescriptor(GATKReadFilterPluginDescriptor.class);
        List<ReadFilter> readFilters = readFilterPlugin.getResolvedInstances();

        Assert.assertEquals(readFilters.size(), 1);
        Assert.assertEquals(readFilters.get(0).getClass().getSimpleName(),
                ReadFilterLibrary.MAPPED.getClass().getSimpleName());

        Assert.assertEquals(readFilterPlugin.userArgs.getUserDisabledReadFilterNames().size(), 2);
        Assert.assertTrue(readFilterPlugin.userArgs.getUserDisabledReadFilterNames().contains(
                ReadFilterLibrary.GOOD_CIGAR.getClass().getSimpleName()));
        Assert.assertTrue(readFilterPlugin.userArgs.getUserDisabledReadFilterNames().contains(
                ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS.getClass().getSimpleName()));

        ReadFilter rf = instantiateFilter(clp, createHeaderWithReadGroups());
        Assert.assertEquals(
                rf.getClass().getSimpleName(),
                ReadFilterLibrary.MAPPED.getClass().getSimpleName());
    }

    @DataProvider(name = "disableToolDefaulFiltersArguments")
    public Object[][] disableToolDefaulFiltersArguments() {
        List<ReadFilter> defaultFilters = new ArrayList<>();
        defaultFilters.add(ReadFilterLibrary.GOOD_CIGAR);
        defaultFilters.add(ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS);
        return new Object[][] {
                {new String[]{
                        "--" + ReadFilterArgumentDefinitions.DISABLE_TOOL_DEFAULT_READ_FILTERS},
                        defaultFilters, null},
                {new String[]{
                        "--" + ReadFilterArgumentDefinitions.DISABLE_TOOL_DEFAULT_READ_FILTERS,
                        "--RF", ReadFilterLibrary.MAPPED.getClass().getSimpleName()},
                        defaultFilters, ReadFilterLibrary.MAPPED},
                {new String[]{
                        "--" + ReadFilterArgumentDefinitions.DISABLE_TOOL_DEFAULT_READ_FILTERS,
                        "--RF", ReadFilterLibrary.GOOD_CIGAR.getClass().getSimpleName()},
                        defaultFilters, ReadFilterLibrary.GOOD_CIGAR}
        };
    }

    @Test(dataProvider = "disableToolDefaulFiltersArguments")
    public void testDisableToolDefaultFilters(final String[] args, final List<ReadFilter> defaultFilters, final ReadFilter expectedFilter) {
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(defaultFilters)),
                Collections.emptySet());

        clp.parseArguments(nullMessageStream, args);

        GATKReadFilterPluginDescriptor readFilterPlugin = clp.getPluginDescriptor(GATKReadFilterPluginDescriptor.class);
        Assert.assertTrue(readFilterPlugin.userArgs.getDisableToolDefaultReadFilters());

        // no instances because no readFilter was provided
        List<ReadFilter> readFilters = readFilterPlugin.getResolvedInstances();
        Assert.assertEquals(readFilters.size(), (expectedFilter == null) ? 0 : 1);

        // all the default filters returns true for isDisabledFilter if it is not the expected
        defaultFilters.stream().map(df -> df.getClass().getSimpleName()).filter(df -> expectedFilter == null || expectedFilter.getClass().getSimpleName().equals(df))
                .forEach(df -> Assert.assertTrue(readFilterPlugin.isDisabledFilter(df.getClass().getSimpleName())));

        ReadFilter rf = instantiateFilter(clp, createHeaderWithReadGroups());
        Assert.assertEquals(
                rf.getClass().getSimpleName(),
                (expectedFilter == null) ? ReadFilterLibrary.ALLOW_ALL_READS.getClass().getSimpleName()
                                         : expectedFilter.getClass().getSimpleName());
    }

    // Disabled due to https://github.com/broadinstitute/barclay/issues/23
    // For now this will just generate a warning, but it should throw
    @Test(expectedExceptions = CommandLineException.class, enabled=false)
    public void testDisabledDefaultWithArgsProvided() {
        // test for arguments provided for a default filter that is also disabled
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(
                        Collections.singletonList(new SampleReadFilter())
                )),
                Collections.emptySet());
        clp.parseArguments(nullMessageStream, new String[] {
                "--sample", "fred",
                "--disable-fead-filter", SampleReadFilter.class.getSimpleName()});

        // get the command line read filters
        clp.getPluginDescriptor(GATKReadFilterPluginDescriptor.class);
    }

    @Test
    public void testDisabledDefaultWithArgsNotProvided() {
        // test for disabling a default filter that has arguments, but they are not provided
        GATKReadFilterPluginDescriptor rfDesc = new GATKReadFilterPluginDescriptor(
                Collections.singletonList(new ReadLengthReadFilter(1, 10)));
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(rfDesc),
                Collections.emptySet());
        final String filterName = ReadLengthReadFilter.class.getSimpleName();
        clp.parseArguments(nullMessageStream, new String[] {
                "--" + ReadFilterArgumentDefinitions.DISABLE_READ_FILTER_LONG_NAME, filterName});

        // Make sure ReadLengthReadFilter got disabled without an exception
        Assert.assertTrue(rfDesc.userArgs.getUserDisabledReadFilterNames().contains(filterName));
        Assert.assertTrue(rfDesc.isDisabledFilter(filterName));
    }

    @Test(expectedExceptions = CommandLineException.class)
    public void testDisabledNonDefaultFilterWithArgsProvided() {
        // test for arguments provided for a non-default filter that is also disabled
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(Collections.emptyList())),
                Collections.emptySet());
        clp.parseArguments(nullMessageStream, new String[] {
                "--" + ReadFilterArgumentDefinitions.DISABLE_TOOL_DEFAULT_READ_FILTERS, SampleReadFilter.class.getSimpleName(),
                "--sample", "fred"
        });
        // get the command line read filters
        clp.getPluginDescriptor(GATKReadFilterPluginDescriptor.class);
    }

    @Test
    public void testDisableDefaultsAndReenableWithDifferentArgs() {
        GATKReadFilterPluginDescriptor rfDesc = new GATKReadFilterPluginDescriptor(
                Collections.singletonList(new ReadLengthReadFilter(1, 10)));
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(rfDesc),
                Collections.emptySet());
        clp.parseArguments(nullMessageStream, new String[] {
                "--" + ReadFilterArgumentDefinitions.DISABLE_TOOL_DEFAULT_READ_FILTERS,
                "--" + ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME, ReadLengthReadFilter.class.getSimpleName(),
                "--" + ReadFilterArgumentDefinitions.MAX_READ_LENGTH_ARG_NAME, "13"}
        );
        List<ReadFilter> allFilters = rfDesc.getResolvedInstances();
        ReadLengthReadFilter rf = (ReadLengthReadFilter) allFilters.get(0);
        Assert.assertEquals(rf.maxReadLength.intValue(), 13);
    }

    // TODO this test is enforcing the expected behavior that if a tool declares a default annotation with an optional argument
    // TODO and if the user then disables and reenables that argument, that the tools default optional argument value will be
    // TODO used over the global default, which seems wrong. See https://github.com/broadinstitute/gatk/issues/3848
    @Test (enabled = false)
    public void testDisableDefaultsAndReenableWithDifferentArgsDefaultValue() {
        GATKReadFilterPluginDescriptor rfDesc = new GATKReadFilterPluginDescriptor(
                Collections.singletonList(new ReadLengthReadFilter(1, 10)));
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(rfDesc),
                Collections.emptySet());
        clp.parseArguments(nullMessageStream, new String[] {
                "--disableToolDefaultReadFilters",
                "--readFilter", ReadLengthReadFilter.class.getSimpleName(),
                "--maxReadLength", "13"}
        );
        List<ReadFilter> allFilters = rfDesc.getResolvedInstances();
        ReadLengthReadFilter rf = (ReadLengthReadFilter) allFilters.get(0);
        Assert.assertEquals(rf.maxReadLength.intValue(), 13);
        Assert.assertEquals(rf.minReadLength, 1);
    }

    @Test
    public void testDisableDefaultsAndReorderWithUserEnabledFirts() {
        final ReadLengthReadFilter rlrf = new ReadLengthReadFilter(2, 10);
        GATKReadFilterPluginDescriptor rfDesc = new GATKReadFilterPluginDescriptor(
                Arrays.asList(rlrf, ReadFilterLibrary.MAPPED));
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(rfDesc),
                Collections.emptySet());
        clp.parseArguments(nullMessageStream, new String[] {
                "--" + ReadFilterArgumentDefinitions.DISABLE_TOOL_DEFAULT_READ_FILTERS,
                "--" + ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME, ReadFilterLibrary.GOOD_CIGAR.getClass().getSimpleName(),
                "--" + ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME, ReadFilterLibrary.MAPPED.getClass().getSimpleName(),
                "--" + ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME, ReadLengthReadFilter.class.getSimpleName()}
        );
        List<ReadFilter> allFilters = rfDesc.getResolvedInstances();
        // User filter is first according to the command line
        Assert.assertSame(allFilters.get(0).getClass(), ReadFilterLibrary.GOOD_CIGAR.getClass());
        // default filters are exactly the same object
        Assert.assertSame(allFilters.get(1), ReadFilterLibrary.MAPPED);
        ReadLengthReadFilter rf = (ReadLengthReadFilter) allFilters.get(2);
        Assert.assertSame(rf, rlrf);
        // and the state is the same as in the provided default
        Assert.assertEquals(rf.minReadLength, 2);
        Assert.assertEquals(rf.maxReadLength.intValue(), 10);
    }

    @Test(expectedExceptions = CommandLineException.class)
    public void testEnableDisableConflict() {
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(null)),
                Collections.emptySet());
        clp.parseArguments(nullMessageStream, new String[] {
                "--RF", "GoodCigarReadFilter",
                "--" + ReadFilterArgumentDefinitions.DISABLE_READ_FILTER_LONG_NAME, "GoodCigarReadFilter"});
    }

    @Test
    public void testDisableToolDefaultFiltersWithInvertedFilter() {
        GATKReadFilterPluginDescriptor rfDesc = new GATKReadFilterPluginDescriptor(
                Collections.singletonList(new ReadLengthReadFilter(1, 10)));
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(rfDesc),
                Collections.emptySet());
        clp.parseArguments(nullMessageStream, new String[] {
                "--" + ReadFilterArgumentDefinitions.DISABLE_TOOL_DEFAULT_READ_FILTERS,
                "--" + ReadFilterArgumentDefinitions.INVERTED_READ_FILTER_LONG_NAME, ReadLengthReadFilter.class.getSimpleName()}
        );
        List<ReadFilter> allFilters = rfDesc.getResolvedInstances();
        Assert.assertEquals(allFilters.size(), 1);
        Assert.assertEquals(allFilters.get(0).getClass(), ReadFilter.ReadFilterNegate.class);
    }

    @Test(expectedExceptions = CommandLineException.class)
    public void testEnableInvertConflict() {
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(null)),
                Collections.emptySet());
        clp.parseArguments(nullMessageStream, new String[] {
                "--RF", "GoodCigarReadFilter",
                "--" + ReadFilterArgumentDefinitions.INVERTED_READ_FILTER_LONG_NAME, "GoodCigarReadFilter"});
    }

    @Test(expectedExceptions = CommandLineException.class)
    public void testInvertToolDefaultConflict() {
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(List.of(new ReadFilterLibrary.GoodCigarReadFilter()))),
                Collections.emptySet());
        clp.parseArguments(nullMessageStream, new String[] {
                "--" + ReadFilterArgumentDefinitions.INVERTED_READ_FILTER_LONG_NAME, "GoodCigarReadFilter"});
    }

    @Test
    public void testPreserveCommandLineOrder() {
        List<ReadFilter> orderedDefaults = new ArrayList<>();
        orderedDefaults.add(new WellformedReadFilter());
        orderedDefaults.add(ReadFilterLibrary.HAS_READ_GROUP);
        orderedDefaults.add(ReadFilterLibrary.MAPPING_QUALITY_NOT_ZERO);

        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(orderedDefaults)),
                Collections.emptySet());
        clp.parseArguments(nullMessageStream, new String[] {
                "--" + ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME, ReadFilterLibrary.MAPPED.getClass().getSimpleName(),
                "--" + ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME, ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS.getClass().getSimpleName(),
                "--" + ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME, ReadFilterLibrary.GOOD_CIGAR.getClass().getSimpleName()});

        GATKReadFilterPluginDescriptor readFilterPlugin = clp.getPluginDescriptor(GATKReadFilterPluginDescriptor.class);

        // get and verify the final resolved list of filters enabled on the command line, resolved with the defaults
        List<ReadFilter> orderedFilters = readFilterPlugin.getResolvedInstances();
        Assert.assertEquals(orderedFilters.size(), 6);
        Assert.assertEquals(orderedFilters.get(0).getClass().getSimpleName(),
                WellformedReadFilter.class.getSimpleName());
        Assert.assertEquals(orderedFilters.get(1).getClass().getSimpleName(),
                ReadFilterLibrary.HAS_READ_GROUP.getClass().getSimpleName());
        Assert.assertEquals(orderedFilters.get(2).getClass().getSimpleName(),
                ReadFilterLibrary.MAPPING_QUALITY_NOT_ZERO.getClass().getSimpleName());
        Assert.assertEquals(orderedFilters.get(3).getClass().getSimpleName(),
                ReadFilterLibrary.MAPPED.getClass().getSimpleName());
        Assert.assertEquals(orderedFilters.get(4).getClass().getSimpleName(),
                ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS.getClass().getSimpleName());
        Assert.assertEquals(orderedFilters.get(5).getClass().getSimpleName(),
                ReadFilterLibrary.GOOD_CIGAR.getClass().getSimpleName());

        // Now get the final merged read filter and verify the execution order. We need to ensure that
        // getMergedReadFilter creates a composite filter that honors the filter test execution
        // order rules (tool defaults first, in order, followed by command line-specified, in order
        // listed), so validate by reaching inside the filter and navigating down the tree.
        ReadFilter rf = instantiateFilter(clp, createHeaderWithReadGroups());
        String expectedOrder[] = {
                WellformedReadFilter.class.getSimpleName(),
                ReadFilterLibrary.HAS_READ_GROUP.getClass().getSimpleName(),
                ReadFilterLibrary.MAPPING_QUALITY_NOT_ZERO.getClass().getSimpleName(),
                ReadFilterLibrary.MAPPED.getClass().getSimpleName(),
                ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS.getClass().getSimpleName(),
                ReadFilterLibrary.GOOD_CIGAR.getClass().getSimpleName()
        };

        int count = ReadFilterUnitTest.verifyAndFilterOrder(rf, expectedOrder);
        Assert.assertEquals(count, 6);
    }

    @Test
    public void testPreserveToolDefaultFilterOrder() {
        List<ReadFilter> orderedDefaults = new ArrayList<>();
        orderedDefaults.add(new WellformedReadFilter());
        orderedDefaults.add(ReadFilterLibrary.MAPPED);
        orderedDefaults.add(ReadFilterLibrary.HAS_READ_GROUP);
        orderedDefaults.add(ReadFilterLibrary.MAPPING_QUALITY_NOT_ZERO);
        orderedDefaults.add(ReadFilterLibrary.PAIRED);
        orderedDefaults.add(ReadFilterLibrary.NONZERO_FRAGMENT_LENGTH_READ_FILTER);
        orderedDefaults.add(ReadFilterLibrary.FIRST_OF_PAIR);
        orderedDefaults.add(ReadFilterLibrary.PROPERLY_PAIRED);
        orderedDefaults.add(ReadFilterLibrary.NOT_DUPLICATE);
        orderedDefaults.add(ReadFilterLibrary.NOT_SECONDARY_ALIGNMENT);
        orderedDefaults.add(ReadFilterLibrary.NOT_SUPPLEMENTARY_ALIGNMENT);

        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(orderedDefaults)),
                Collections.emptySet());
        clp.parseArguments(nullMessageStream, new String[] {
                //disable one just to mix things up
                "-" + ReadFilterArgumentDefinitions.DISABLE_READ_FILTER_SHORT_NAME, ReadFilterLibrary.MAPPED.getClass().getSimpleName(),
                "-" + ReadFilterArgumentDefinitions.READ_FILTER_SHORT_NAME, ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS.getClass().getSimpleName(),
                "-" + ReadFilterArgumentDefinitions.READ_FILTER_SHORT_NAME, ReadFilterLibrary.GOOD_CIGAR.getClass().getSimpleName()});

        // Now get the final merged read filter and verify the execution order. We need to ensure that
        // getMergedReadFilter creates a composite filter that honors the filter test execution
        // order rules (tool defaults first, in order, followed by command line-specified, in order
        // listed). So reach inside the filter and navigate down the tree.
        ReadFilter rf = instantiateFilter(clp, createHeaderWithReadGroups());
        String expectedOrder[] = {
                WellformedReadFilter.class.getSimpleName(),
                ReadFilterLibrary.HAS_READ_GROUP.getClass().getSimpleName(),
                ReadFilterLibrary.MAPPING_QUALITY_NOT_ZERO.getClass().getSimpleName(),
                ReadFilterLibrary.PAIRED.getClass().getSimpleName(),
                ReadFilterLibrary.NONZERO_FRAGMENT_LENGTH_READ_FILTER.getClass().getSimpleName(),
                ReadFilterLibrary.FIRST_OF_PAIR.getClass().getSimpleName(),
                ReadFilterLibrary.PROPERLY_PAIRED.getClass().getSimpleName(),
                ReadFilterLibrary.NOT_DUPLICATE.getClass().getSimpleName(),
                ReadFilterLibrary.NOT_SECONDARY_ALIGNMENT.getClass().getSimpleName(),
                ReadFilterLibrary.NOT_SUPPLEMENTARY_ALIGNMENT.getClass().getSimpleName(),
                ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS.getClass().getSimpleName(),
                ReadFilterLibrary.GOOD_CIGAR.getClass().getSimpleName()
        };

        int count = ReadFilterUnitTest.verifyAndFilterOrder(rf, expectedOrder);
        Assert.assertEquals(count, 12);
    }

    @Test
    public void testReadLengthFilter() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);

        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(null)),
                Collections.emptySet());
        String[] args = {
                "--" + ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME, ReadLengthReadFilter.class.getSimpleName(),
                "--" + ReadFilterArgumentDefinitions.MIN_READ_LENGTH_ARG_NAME, "10",
                "--" + ReadFilterArgumentDefinitions.MAX_READ_LENGTH_ARG_NAME, "20"
        };
        clp.parseArguments(nullMessageStream, args);
        ReadFilter rf = instantiateFilter(clp, header);

        read.setBases(new byte[5]);
        Assert.assertFalse(rf.test(read));
        read.setBases(new byte[25]);
        Assert.assertFalse(rf.test(read));

        read.setBases(new byte[15]);
        Assert.assertTrue(rf.test(read));
    }

    @Test
    public void testInvertReadLengthFilter() {
        final SAMFileHeader header = createHeaderWithReadGroups();
        final GATKRead read = simpleGoodRead(header);

        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKReadFilterPluginDescriptor(null)),
                Collections.emptySet());
        String[] args = {
                "--" + ReadFilterArgumentDefinitions.INVERTED_READ_FILTER_LONG_NAME, ReadLengthReadFilter.class.getSimpleName(),
                "--" + ReadFilterArgumentDefinitions.MIN_READ_LENGTH_ARG_NAME, "10",
                "--" + ReadFilterArgumentDefinitions.MAX_READ_LENGTH_ARG_NAME, "20"
        };
        clp.parseArguments(nullMessageStream, args);
        ReadFilter rf = instantiateFilter(clp, header);

        read.setBases(new byte[5]);
        Assert.assertTrue(rf.test(read));
        read.setBases(new byte[25]);
        Assert.assertTrue(rf.test(read));

        read.setBases(new byte[15]);
        Assert.assertFalse(rf.test(read));
    }

    final private static String readgroupName = "ReadGroup1";

    @DataProvider(name="testDefaultFilters")
    public Object[][] defaultFilters() {
        return new Object[][]{
                // ReadGroupReadFilter has a required "keepReadGroup" arg; provide it
                {new String[]{"--" + ReadFilterArgumentDefinitions.KEEP_READ_GROUP_LONG_NAME, readgroupName}},

                // ReadGroupReadFilter has a required "keepReadGroup" arg; provide it
                // *and* specify it on the command line
                {new String[]{
                        "--" + ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME, "ReadGroupReadFilter",
                        "--" + ReadFilterArgumentDefinitions.KEEP_READ_GROUP_LONG_NAME, readgroupName}
                }
        };
    }

    @Test(dataProvider = "testDefaultFilters")
    public void testToolHasDefaultRequiredArgsPositive(final String[] args) {
        CommandLineParser clp = new CommandLineArgumentParser(
            new Object(),
            Collections.singletonList(
                    new GATKReadFilterPluginDescriptor(Collections.singletonList(new ReadGroupReadFilter()))
            ),
            Collections.emptySet()
        );
        clp.parseArguments(nullMessageStream, args);

        SAMFileHeader samHeader = createHeaderWithReadGroups();
        ReadFilter rf = instantiateFilter(clp, samHeader);
        final GATKRead read = simpleGoodRead(samHeader);

        // make sure the test read actually has this read group
        Assert.assertEquals(read.getReadGroup(), readgroupName);
        Assert.assertTrue(rf.test(read));
    }

    @Test(expectedExceptions = CommandLineException.MissingArgument.class)
    public void testToolHasDefaultRequiredArgsNegative() {
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(
                        new GATKReadFilterPluginDescriptor(Collections.singletonList(new ReadGroupReadFilter()))
                ),
                Collections.emptySet());

        // ReadGroupReadFilter has a required "keepReadGroup" arg; don't provide it and fail
        String[] args = {};
        clp.parseArguments(nullMessageStream, args);
    }

    private ReadFilter instantiateFilter(
            final CommandLineParser clp,
            final SAMFileHeader header)
    {
        GATKReadFilterPluginDescriptor readFilterPlugin = clp.getPluginDescriptor(GATKReadFilterPluginDescriptor.class);
        return readFilterPlugin.getMergedReadFilter(header);
    }

    private static final int CHR_COUNT = 2;
    private static final int CHR_START = 1;
    private static final int CHR_SIZE = 1000;
    private static final int GROUP_COUNT = 5;

    private SAMFileHeader createHeaderWithReadGroups() {
        return ArtificialReadUtils.createArtificialSamHeaderWithGroups(CHR_COUNT, CHR_START, CHR_SIZE, GROUP_COUNT);
    }

    /**
     * Creates a read record.
     *
     * @param header header for the new record
     * @param cigar the new record CIGAR.
     * @param group the new record group index that must be in the range \
     *              [0,{@link #GROUP_COUNT})
     * @param reference the reference sequence index (0-based)
     * @param start the start position of the read alignment in the reference
     *              (1-based)
     * @return never <code>null</code>
     */
    private GATKRead createRead(final SAMFileHeader header, final Cigar cigar, final int group,
                                final int reference, final int start ) {
        final GATKRead record = ArtificialReadUtils.createArtificialRead(header, cigar);
        record.setPosition(header.getSequence(reference).getSequenceName(), start);
        record.setReadGroup(header.getReadGroups().get(group).getReadGroupId());
        return record;
    }

    private GATKRead simpleGoodRead( final SAMFileHeader header ) {
        final String cigarString = "101M";
        final Cigar cigar = TextCigarCodec.decode(cigarString);
        GATKRead read = createRead(header, cigar, 1, 0, 10);
        read.setMappingQuality(50);
        return read;
    }

    // object for holding test setup params to use as a type for test lambdas
    public class SetupTest {
        public SAMFileHeader hdr;
        public GATKRead read;
        public String argValue;

        public SetupTest(SAMFileHeader hdr, GATKRead read, String argValue) {
            this.hdr = hdr;
            this.read = read;
            this.argValue = argValue;
        }
    }

    // setup methods to initialize reads for individual filter tests
    private void setupLibraryTest(SetupTest setup) {
        setup.hdr.getReadGroup(setup.read.getReadGroup()).setLibrary(setup.argValue);
    }

    private void setupMappingQualityTest(SetupTest setup) {
        setup.read.setMappingQuality(Integer.parseInt(setup.argValue));
    }

    private void setupReadGroupTest(SetupTest setup) {
        setup.read.setReadGroup(setup.argValue);
    }

    private void setupReadNameTest(SetupTest setup) {
        setup.read.setName(setup.argValue);
    }

    private void setupSampleTest(SetupTest setup) {
        setup.hdr.getReadGroup(setup.read.getReadGroup()).setSample(setup.argValue);
    }
    
    @CommandLineProgramProperties(summary="test tool to check ReadFilter loading",
            oneLineSummary = "test tool to check ReadFilter loading",
            programGroup = TestProgramGroup.class,
            omitFromCommandLine = true)
    public static class TestReadFiltersTool extends GATKTool {
        @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
        String output;

        @Override
        public void traverse() {
            String readFilters = getCommandLineParser()
                    .getPluginDescriptor(GATKReadFilterPluginDescriptor.class)
                    .getResolvedInstances()
                    .stream()
                    .map(a -> a.getClass().getSimpleName())
                    .collect(Collectors.joining("\n"));
            IOUtil.writeText(readFilters, new File(output));
        }
    }

    @Test
    public void testConfigFileControlsReadFiltersPackages() throws IOException {
        final File output = createTempFile("annotations", "txt");
        final File configFile = new File(packageRootTestDir + "cmdline/GATKPlugin/changePluginPackages.properties");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add(ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME, TestReadFilter.class.getSimpleName())
                .add(StandardArgumentDefinitions.GATK_CONFIG_FILE_OPTION, configFile)
                .addOutput(output);
        runToolInNewJVM("TestReadFiltersTool", args);
        Assert.assertEquals(Files.readAllLines(output.toPath()), Arrays.asList(WellformedReadFilter.class.getSimpleName(), TestReadFilter.class.getSimpleName()));
    }

}
