package org.broadinstitute.hellbender.testutils;

import htsjdk.samtools.util.Locatable;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Builder for command line argument lists with lots of convenience methods for adding standard GATK arguments such as
 * input, output, variants, reference, and intervals.
 * 
 * Use this only in test code.
 */
public final class ArgumentsBuilder {
    private final List<String> args= new ArrayList<>();

    public ArgumentsBuilder(){}

    /**
     * static factory to allow fluent style creation
     * create().add()...
     * @return new ArgumentsBuilder
     */
    public static ArgumentsBuilder create(){
        return new ArgumentsBuilder();
    }

    public ArgumentsBuilder(Object[] args){
        for (Object arg: args){
            if (arg instanceof String){
                addRaw((String) arg);
            } else {
                addRaw(arg);
            }
        }
    }

    /**
     * Concatenate the arguments from other onto the end of this builder
     * @return this builder combined with other
     */
    public ArgumentsBuilder concat(ArgumentsBuilder other){
        this.args.addAll(other.args);
        return this;
    }

    /**
     * @return a copy of this builder
     */
    public ArgumentsBuilder copy(){
        return ArgumentsBuilder.create().concat(this);
    }

    /**
     * Add a string to the arguments list
     * Strings are processed specially, they are reformatted to match the new unix style arguments
     *
     * NOTE: In general this method should be avoided in favor of other methods that handle adding dashes.
     * @param arg A string representing one or more arguments
     * @return the ArgumentsBuilder
     */
    public ArgumentsBuilder addRaw(String arg){
        List<String> chunks = Arrays.asList(StringUtils.split(arg.trim()));
        for (String chunk : chunks){
            args.add(chunk);
        }
        return this;
    }

    /**
     * Add any object's string representation to the arguments list
     */
    public ArgumentsBuilder addRaw(Object arg) {
        args.add(arg.toString());
        return this;
    }

    // ARGUMENT/VALUE METHODS

    /**
     * add an argument with a given value to this builder.
     *
     * This is the fundamental add method that others invoke.  It adds dashes to the argument name.
     */
    public ArgumentsBuilder add(final String argumentName, final String argumentValue) {
        Utils.nonNull(argumentValue);
        Utils.nonNull(argumentName);
        addRaw("--" + argumentName);
        addRaw(argumentValue);
        return this;
    }

    public ArgumentsBuilder add(final String argumentName, final File file){
        Utils.nonNull(file);
        return add(argumentName, file.getAbsolutePath());
    }

    public ArgumentsBuilder add(final String argumentName, final Path path){
        Utils.nonNull(path);
        return add(argumentName, path.toString());
    }

    public ArgumentsBuilder add(final String argumentName, final boolean yes){
        return add(argumentName, String.valueOf(yes));
    }

    public ArgumentsBuilder add(final String argumentName, final Number value){
        Utils.nonNull(value);
        return add(argumentName, value.toString());
    }

    public ArgumentsBuilder add(final String argumentName, final Enum<?> enummerationValue){
        Utils.nonNull(enummerationValue);
        return add(argumentName, enummerationValue.name());
    }

    // CONVENIENCE METHODS WITH BUILT-IN STANDARD ARGUMENTS

    // INPUT

    public ArgumentsBuilder addInput(final File input) {
        return add(StandardArgumentDefinitions.INPUT_LONG_NAME, input);
    }

    public ArgumentsBuilder addInput(final String input) {
        return add(StandardArgumentDefinitions.INPUT_LONG_NAME, input);
    }

    public ArgumentsBuilder addInput(final Path input) {
        return add(StandardArgumentDefinitions.INPUT_LONG_NAME, input);
    }

    // OUTPUT

    public ArgumentsBuilder addOutput(final File output) {
        return add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, output);
    }

    public ArgumentsBuilder addOutput(final String output) {
        return add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, output);
    }

    public ArgumentsBuilder addOutput(final Path output) {
        return add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, output.toString());
    }

    // REFERENCE

    public ArgumentsBuilder addReference(final File reference){
        return add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, reference);
    }

    public ArgumentsBuilder addReference(final String reference){
        return add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, reference);
    }

    public ArgumentsBuilder addReference(final Path reference){
        return add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, reference.toString());
    }

    // VCF

    public ArgumentsBuilder addVCF(final File vcf) {
        return add(StandardArgumentDefinitions.VARIANT_LONG_NAME, vcf);
    }

    public ArgumentsBuilder addVCF(final String vcf) {
        return add(StandardArgumentDefinitions.VARIANT_LONG_NAME, vcf);
    }

    public ArgumentsBuilder addVCF(final Path vcf) {
        return add(StandardArgumentDefinitions.VARIANT_LONG_NAME, vcf);
    }


    //FLAG

    public ArgumentsBuilder addFlag(final String argumentName) {
        Utils.nonNull(argumentName);
        return addRaw("--" + argumentName);
    }

    // INTERVALS

    public ArgumentsBuilder addInterval(final String interval){
        Utils.nonNull(interval);
        return add(StandardArgumentDefinitions.INTERVALS_LONG_NAME, interval);
    }

    public ArgumentsBuilder addInterval(final Locatable interval){
        Utils.nonNull(interval);
        return add(StandardArgumentDefinitions.INTERVALS_LONG_NAME, IntervalUtils.locatableToString(interval));
    }

    public ArgumentsBuilder addIntervals(final File interval){
        Utils.nonNull(interval);
        return add(StandardArgumentDefinitions.INTERVALS_LONG_NAME, interval);
    }

    public ArgumentsBuilder addMask(final File mask){
        return add(IntervalArgumentCollection.EXCLUDE_INTERVALS_LONG_NAME, mask);
    }

    /**
     * @return the arguments as List
     */
    public List<String> getArgsList(){
        return args;
    }

    /**
     * @return the arguments as String[]
     */
    public String[] getArgsArray(){
        return args.toArray(new String[this.args.size()]);
    }


    /**
     * @return the arguments as a single String
     */
    public String getString() {
        return String.join(" ", args);
    }

    @Override
    public String toString(){
        return getString();
    }
}
