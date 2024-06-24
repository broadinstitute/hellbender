package org.broadinstitute.hellbender.cmdline.GATKPlugin;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLinePluginDescriptor;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.config.ConfigFactory;
import org.broadinstitute.hellbender.utils.config.GATKConfig;

import java.io.IOException;
import java.util.*;
import java.util.function.BiFunction;
import java.util.stream.Collectors;

/**
 * A CommandLinePluginDescriptor for ReadFilter plugins
 */
public class GATKReadFilterPluginDescriptor extends CommandLinePluginDescriptor<ReadFilter> {

    protected transient Logger logger = LogManager.getLogger(this.getClass());

    private void readObject(java.io.ObjectInputStream in)
            throws IOException, ClassNotFoundException {
        in.defaultReadObject();
        logger = LogManager.getLogger(this.getClass()); // Logger is not serializable (even by Kryo)
    }

    /**
     * At startup, set the plugin package name to the one(s) in the configuration file.
     */
    private static final List<String> PLUGIN_PACKAGE_NAMES;
    static {
        // Get our configuration:
        final GATKConfig config = ConfigFactory.getInstance().getGATKConfig();
        // Exclude abstract classes and interfaces from the list of discovered codec classes
        PLUGIN_PACKAGE_NAMES = Collections.unmodifiableList(config.read_filter_packages());
    }

    private static final Class<ReadFilter> pluginBaseClass = org.broadinstitute.hellbender.engine.filters.ReadFilter.class;
    private static final String READ_PLUGIN_DISPLAY_NAME = "readFilter";

    // the purpose of this argument collection is to allow the caller to control the exposure of the command line arguments
    @VisibleForTesting
    @ArgumentCollection
    final GATKReadFilterArgumentCollection userArgs;

    // Map of read filter (simple) class names to the corresponding discovered plugin instance
    private final Map<String, ReadFilter> allDiscoveredReadFilters = new HashMap<>();

    // Map of read filter (simple) class names to the corresponding default plugin instance
    // it is a LinkedHashMap because we want to remember the order in which these were provided, and also keep the
    // actual instances in case they have any additional state provided by the tool
    // when they were created
    private final Map<String, ReadFilter> toolDefaultReadFilters = new LinkedHashMap<>();

    // Set of predecessor readFilters for which we've seen arguments that must exist either as a tool default or be supplied by the user
    // (eg. ReadLengthReadFilter if we see "--maxReadLength" on the command line)
    private final Set<String> requiredPredecessors = new HashSet<>();

    /**
     * @param userArgs           Argument collection to control the exposure of the command line arguments.
     * @param toolDefaultFilters Default filters that may be supplied with arguments
     *                           on the command line. May be null.
     */
    public GATKReadFilterPluginDescriptor(final GATKReadFilterArgumentCollection userArgs, final List<ReadFilter> toolDefaultFilters) {
        this.userArgs = userArgs;
        if (null != toolDefaultFilters) {
            toolDefaultFilters.forEach(f -> {
                final Class<? extends ReadFilter> rfClass = f.getClass();
                // anonymous classes have a 0-length simple name, and thus cannot be accessed or
                // controlled by the user via the command line, but they should still be valid
                // as default filters, so use the full name to ensure that their map entries
                // don't clobber each other
                String className = rfClass.getSimpleName();
                if (className.length() == 0) {
                    className = rfClass.getName();
                }
                toolDefaultReadFilters.put(className, f);
            });
        }
    }

    /**
     * @param toolDefaultFilters Default filters that may be supplied with arguments
     *                           on the command line. May be null.
     */
    public GATKReadFilterPluginDescriptor(final List<ReadFilter> toolDefaultFilters) {
        this(new DefaultGATKReadFilterArgumentCollection(), toolDefaultFilters);
    }

    /////////////////////////////////////////////////////////
    // GATKCommandLinePluginDescriptor implementation methods

    /**
     * Return a display name to identify this plugin to the user
     * @return A short user-friendly name for this plugin.
     */
    @Override
    public String getDisplayName() {
        // The value returned by this method is placed into the freemarker property map by the docgen system,
        // and becomes a variable in the freemarker template language. We can't use the actual argument name
        // here because the kebabified version ("read-filter") isn't a valid freemarker variable (and is
        // parsed as an expression with an embedded '-'). So use the pre-kebabified version, which matches
        // what is used in the template).
        return READ_PLUGIN_DISPLAY_NAME;
    }

    /**
     * @return the class object for the base class of all plugins managed by this descriptor
     */
    @Override
    public Class<ReadFilter> getPluginBaseClass() {return pluginBaseClass;}

    /**
     * A list of package names which will be searched for plugins managed by the descriptor.
     * @return
     */
    @Override
    public List<String> getPackageNames() {return PLUGIN_PACKAGE_NAMES;}

    @Override
    public boolean includePluginClass(Class<?> c) {
        // don't use the ReadFilter base class, it's inner classes, the CountingReadFilter,
        // or the unit tests
        return !c.getName().equals(this.getPluginBaseClass().getName()) &&
                !c.getName().startsWith(CountingReadFilter.class.getName()) &&
                !c.getName().startsWith(this.getPluginBaseClass().getName() + "$") &&
                !c.getName().contains("UnitTest$");
    }

    // Instantiate a new ReadFilter derived object and save it in the list
    @Override
    @SuppressWarnings("deprecation")
    public ReadFilter createInstanceForPlugin(final Class<?> pluggableClass) throws IllegalAccessException, InstantiationException {
        ReadFilter readFilter = null;
        final String simpleName = pluggableClass.getSimpleName();

        if (allDiscoveredReadFilters.containsKey(simpleName)) {
            // we found a plugin class with a name that collides with an existing class;
            // plugin names must be unique even across packages
            throw new IllegalArgumentException(
                    String.format("A plugin class name collision was detected (%s/%s). " +
                            "Simple names of plugin classes must be unique across packages.",
                            pluggableClass.getName(),
                            allDiscoveredReadFilters.get(simpleName).getClass().getName())
            );
        } else if (toolDefaultReadFilters.containsKey(simpleName)) {
            // an instance of this class was provided by the tool as one of it's default filters;
            // use the default instance as the target for command line argument values
            // rather than creating a new one, in case it has state provided by the tool
            readFilter = toolDefaultReadFilters.get(simpleName);
        } else {
            readFilter = (ReadFilter) pluggableClass.newInstance();
        }

        // Add all filters to the allDiscoveredReadFilters list, even if the instance came from the
        // tool defaults list (we want the actual instances to be shared to preserve state)
        allDiscoveredReadFilters.put(simpleName, readFilter);
        return readFilter;
    }

    @Override
    public boolean isDependentArgumentAllowed(final Class<?> predecessorClass) {
        // Make sure the predecessor for a dependent argument was either specified on the command line or
        // is a tool default, otherwise reject it.
        // NOTE: This method is called by the CLP during parsing at the time the dependet argument is seen
        // on the command line. Even if this check passes at the time this method is called, its possible
        // for the user to subsequently disable the required predecessor. That case is caught during final
        // validation done by the validateArguments method.
        String predecessorName = predecessorClass.getSimpleName();
        boolean isAllowed = userArgs.getAllUserEnabledReadFilterNames().contains(predecessorName)
                || (toolDefaultReadFilters.get(predecessorName) != null);
        if (isAllowed) {
            // Keep track of the ones we allow so we can validate later that they weren't subsequently disabled
            requiredPredecessors.add(predecessorName);
        }
        return isAllowed;
    }

    /**
     * @return an ordered list of ReadFilter instances that result from resolving all command line
     * arguments with any default plugins that have been provided to this descriptor. This list
     * represents the ReadFilters that will actually be used by the consumer.
     *
     * NOTE: calling this method before argument parsing (and thus before {@link #validateAndResolvePlugins}
     * has been called) may return a different list than calling it after parsing, since the resolution
     * policy is implementation-dependent and may depend on the actual arguments specified by the user
     * on the command line.
     */
    @Override
    public List<ReadFilter> getResolvedInstances() {
        // start with the tool's default filters in the order they were specified, and remove any that
        // were disabled on the command line (if --disableToolDefaultReadFilters is specified, just initialize
        // an empty list with initial capacity of user filters)
        final List<ReadFilter> finalFilters =
                userArgs.getDisableToolDefaultReadFilters() ?
                        new ArrayList<>(userArgs.getUserEnabledReadFilterNames().size()) :
                        toolDefaultReadFilters.entrySet()
                                .stream()
                                .filter(e -> !isDisabledFilter(e.getKey()))
                                .map(e -> e.getValue())
                                .collect(Collectors.toList());

        // now add in any additional filters enabled on the command line (preserving order)
        final List<ReadFilter> userEnabledFilters = getUserEnabledInstances();
        if (userEnabledFilters != null) {
            userEnabledFilters.stream()
                    .filter(f -> !finalFilters.contains(f)) // remove redundant filters
                    .forEach(f -> finalFilters.add(f));
        }
        return finalFilters;
    }

    // Return a list with a read filter instance for each read filter enabled by the user.
    private List<ReadFilter> getUserEnabledInstances() {
        final ArrayList<ReadFilter> filters = new ArrayList<>(userArgs.getUserEnabledReadFilterNames().size());
        userArgs.getUserEnabledReadFilterNames().forEach(s -> {
            ReadFilter rf = allDiscoveredReadFilters.get(s);
            filters.add(rf);
        });
        userArgs.getUserEnabledInvertedReadFilterNames().forEach(s -> {
            ReadFilter rf = allDiscoveredReadFilters.get(s).negate();
            filters.add(rf);
        });
        return filters;
    }

    /**
     * Get the list of default plugins used for this instance of this descriptor. Used for help/doc generation.
     *
     * NOTE: this method does not account for disabled default filters and just return ALL default instances.
     * The refactored interface in Barclay changes it's contract to allows returning a list with only 'enabled' default
     * instances. We'll change the implementation when we integrate the updated interface.
     */
    @Override
    public List<ReadFilter> getDefaultInstances() { return new ArrayList<>(toolDefaultReadFilters.values()); }

    /**
     * Return the allowed values for read-filter/disable-read-filter names for use by the help system.
     * @param longArgName long name of the argument for which help is requested
     * @return
     */
    @Override
    public Set<String> getAllowedValuesForDescriptorHelp(final String longArgName) {
        if (longArgName.equals(ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME)) {
            return allDiscoveredReadFilters.keySet();
        }
        if (longArgName.equals(ReadFilterArgumentDefinitions.DISABLE_READ_FILTER_LONG_NAME)) {
            return toolDefaultReadFilters.keySet();
        }
        return null;
    }

    /**
     * Return the class object for the plugin with simple class name {@code pluginName}
     * Used for help/usage and documentation generation.
     *
     * @param pluginName Name of the plugin requested
     * @return Class object for the plugin instance requested, upper bounded by type {@code T}
     */
    @Override
    public Class<?> getClassForPluginHelp(final String pluginName) {
        if (allDiscoveredReadFilters.get(pluginName) != null) {
            return allDiscoveredReadFilters.get(pluginName).getClass();
        } else {
            throw new IllegalArgumentException(String.format("Can't resolve ReadFilter plugin for name %s", pluginName));
        }
    }

    /**
     * Validate the list of arguments and reduce the list of read filters to those
     * actually seen on the command line. This is called by the command line parser
     * after all arguments have been parsed.
     */
    @Override
    public void validateAndResolvePlugins() {
        // throw if a filter is *enabled* more than once by the user
        final Set<String> duplicateUserEnabledFilterNames = Utils.getDuplicatedItems(userArgs.getUserEnabledReadFilterNames());
        if (!duplicateUserEnabledFilterNames.isEmpty()) {
            throw new CommandLineException.BadArgumentValue(
                    String.format("The read filter(s) are enabled more than once: %s",
                            Utils.join(", ", duplicateUserEnabledFilterNames)));
        }

        // throw if a filter is *disabled* more than once by the user
        final Set<String> duplicateDisabledUserFilterNames = Utils.getDuplicatedItems(userArgs.getUserDisabledReadFilterNames());
        if (!duplicateDisabledUserFilterNames.isEmpty()) {
            throw new CommandLineException.BadArgumentValue(
                    String.format("The read filter(s) are disabled more than once: %s",
                            Utils.join(", ", duplicateDisabledUserFilterNames)));
        }

        // throw if a filter is both enabled *and* disabled by the user
        final Set<String> enabledAndDisabled = new HashSet<>(userArgs.getAllUserEnabledReadFilterNames());
        enabledAndDisabled.retainAll(userArgs.getUserDisabledReadFilterNames());
        if (!enabledAndDisabled.isEmpty()) {
            final String badFiltersList = Utils.join(", ", enabledAndDisabled);
            throw new CommandLineException(
                    String.format("The read filter(s): %s are both enabled and disabled", badFiltersList));
        }

        // throw if a filter is enabled and inverted by the user
        final Set<String> enabledAndInverted = new HashSet<>(userArgs.getUserEnabledReadFilterNames());
        enabledAndInverted.retainAll(userArgs.getUserEnabledInvertedReadFilterNames());
        if (!enabledAndInverted.isEmpty()) {
            final String badFiltersList = Utils.join(", ", enabledAndDisabled);
            throw new CommandLineException(
                    String.format("The read filter(s): %s are both enabled and inverted", badFiltersList));
        }

        // throw if a disabled filter doesn't exist; warn if it wasn't enabled by the tool in the first place
        userArgs.getUserDisabledReadFilterNames().forEach(s -> {
            if (!allDiscoveredReadFilters.containsKey(s)) {
                throw new CommandLineException.BadArgumentValue(String.format("Disabled filter (%s) does not exist", s));
            } else if (!toolDefaultReadFilters.containsKey(s)) {
                logger.warn(String.format("Disabled filter (%s) is not enabled by this tool", s));
            }
        });

        // warn if a filter is both default and enabled by the user
        final Set<String> redundant = new HashSet<>(toolDefaultReadFilters.keySet());
        redundant.retainAll(userArgs.getUserEnabledReadFilterNames());
        redundant.forEach(
            s -> {
                logger.warn(String.format("Redundant enabled filter (%s) is enabled for this tool by default", s));
            });

        // throw if a filter is both default and inverted by the user as we do not want to inadvertently filter all reads from the input
        final Set<String> redundantInv = new HashSet<>(toolDefaultReadFilters.keySet());
        redundantInv.retainAll(userArgs.getUserEnabledInvertedReadFilterNames());
        if (!userArgs.getDisableToolDefaultReadFilters() && redundantInv.size() > 0) {
            throw new CommandLineException(
                    String.format("The read filter(s): %s exist as tool default filters and were inverted, disable tool default read filters in order to use this filter", redundantInv));
        }

        // Throw if args were specified for a filter that was also disabled, or that was not enabled by the
        // tool by default.
        //
        // Note that this is also checked during command line argument parsing, but needs to be checked again
        // here. Whenever the command line parser sees a dependent argument on the command line, it delegates
        // back to the descriptor's isDependentArgumentAllowed method to allow it to validate that the predecessor
        // for that dependent argument has been supplied, either by a default read filter, or by an explicitly
        // enabled read filter. However, its possible for the user to subsequently try to disable that
        // predecessor, which is what we want to catch here.
        //
        userArgs.getUserDisabledReadFilterNames().forEach(s -> {
            if (requiredPredecessors.contains(s)) {
                String message = String.format("Values were supplied for (%s) that is also disabled", s);
                if (toolDefaultReadFilters.containsKey(s)) {
                    // NOTE: https://github.com/broadinstitute/barclay/issues/23
                    // This is a special case to work around the issue where we can't really tell if the
                    // predecessor was added as a result of a user-provided value, or a default value. The
                    // CLP doesn't distinguish, so we only warn here for now.
                    logger.warn(message);
                } else {
                    throw new CommandLineException(message);
                }
            }
        });

        // throw if a filter name was specified that has no corresponding instance
        final Map<String, ReadFilter> requestedReadFilters = new HashMap<>();
        userArgs.getAllUserEnabledReadFilterNames().forEach(s -> {
            ReadFilter trf = allDiscoveredReadFilters.get(s);
            if (null == trf) {
                if (!toolDefaultReadFilters.containsKey(s)) {
                    throw new CommandLineException("Unrecognized read filter name: " + s);
                }
            } else {
                requestedReadFilters.put(s, trf);
            }
        });
    }

    /////////////////////////////////////////////////////////
    // ReadFilter plugin-specific helper methods

    /**
     * Determine if a particular ReadFilter was disabled on the command line, either directly of by disabling all
     * tool defaults.
     * @param filterName name of the filter to query.
     * @return {@code true} if the name appears in the list of disabled filters, or is a tool default not provided by
     * the user and all tool defaults are disabled; {@code false} otherwise.
     */
    public boolean isDisabledFilter(final String filterName) {
        return userArgs.getUserDisabledReadFilterNames().contains(filterName)
                || (userArgs.getDisableToolDefaultReadFilters() && !userArgs.getAllUserEnabledReadFilterNames().contains(filterName));
    }

    /**
     * Merge the default filters with the users's command line read filter requests, then initialize
     * the resulting filters.
     *
     * @param samHeader - a SAMFileHeader to use to initialize read filter instances
     * @return Single merged read filter.
     */
    public final ReadFilter getMergedReadFilter(final SAMFileHeader samHeader) {
        Utils.nonNull(samHeader);
        return getMergedReadFilter(
                samHeader,
                ReadFilter::fromList
        );
    }

    /**
     * Merge the default filters with the users's command line read filter requests, then initialize
     * the resulting filters.
     *
     * @param samHeader - a SAMFileHeader to use to initialize read filter instances
     * @return Single merged counting read filter.
     */
    public final CountingReadFilter getMergedCountingReadFilter(final SAMFileHeader samHeader) {
        Utils.nonNull(samHeader);
        return getMergedReadFilter(
                samHeader,
                CountingReadFilter::fromList
        );
    }

    /**
     * Merge the default filters with the users's command line read filter requests, then initialize
     * the resulting filters.
     *
     * @param samHeader a SAMFileHeader to initialize read filter instances. May not be null.
     * @param aggregateFunction function to use to merge ReadFilters, usually ReadFilter::fromList. The function
     *                          must return the ALLOW_ALL_READS filter wrapped in the appropriate type when passed
     *                          a null or empty list.
     * @param <T> extends ReadFilter, type returned by the wrapperFunction
     * @return Single merged read filter.
     */
    public <T extends ReadFilter> T getMergedReadFilter(
            final SAMFileHeader samHeader,
            final BiFunction<List<ReadFilter>, SAMFileHeader, T> aggregateFunction) {

        Utils.nonNull(samHeader);
        Utils.nonNull(aggregateFunction);

        final List<ReadFilter> finalFilters = getResolvedInstances();
        return aggregateFunction.apply(finalFilters, samHeader);
    }

}
