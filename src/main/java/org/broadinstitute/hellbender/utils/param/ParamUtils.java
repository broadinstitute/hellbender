package org.broadinstitute.hellbender.utils.param;

import org.apache.commons.lang3.Range;
import org.apache.commons.math3.exception.NotFiniteNumberException;
import org.apache.commons.math3.util.MathUtils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.*;

/**
 * This class should eventually be merged into Utils, which is in hellbender, and then this class should be deleted.
 *
 * Double.NaN will generally cause a check to fail.  Note that any comparison with a NaN yields false.
 *
 * Created by lichtens on 8/25/15.
 */
public class ParamUtils {

    public static final double INV_LOG_2 = 1/Math.log(2.0);

    private ParamUtils () {}

    /**
     * Checks that the  input is within range and returns the same value or throws an {@link IllegalArgumentException}
     *
     * <p>Note that min can be greater than max, but that will guarantee that an exception is thrown.</p>
     *
     * @param val value to check
     * @param min minimum value for val
     * @param max maximum value for val
     * @param message the text message that would be pass to the exception thrown when {@code o == null}.
     * @return the same value
     * @throws IllegalArgumentException
     */
    public static double inRange(final double val, final double min, final double max, final String message) {
        Utils.validateArg(val >= min && val <= max, message);
        return val;
    }

    /**
     * Checks that the  input is within range and returns the same value or throws an {@link IllegalArgumentException}
     * @param val value to check
     * @param min minimum value for val
     * @param max maximum value for val
     * @param message the text message that would be pass to the exception thrown when val gt min or val lt max.
     * @return the same value
     * @throws IllegalArgumentException
     */
    public static long inRange(final long val, final double min, final double max, final String message) {
        Utils.validateArg(val >= min && val <= max, message);
        return val;
    }

    /**
     * Checks that the  input is within range and returns the same value or throws an {@link IllegalArgumentException}
     * @param val value to check
     * @param min minimum value for val (inclusive)
     * @param max maximum value for val (inclusive)
     * @param message the text message that would be pass to the exception thrown when val gt min or val lt max.
     * @return the same value
     * @throws IllegalArgumentException
     */
    public static long inRange(final long val, final long min, final long max, final String message) {
        Utils.validateArg(val >= min && val <= max, message);
        return val;
    }

    /**
     * Checks that the input int value is within range and returns the same value or throws an {@link IllegalArgumentException}
     * @param val value to check
     * @param min minimum value for val (inclusive)
     * @param max maximum value for val (inclusive)
     * @param message the text message that would be pass to the exception thrown when val gt min or val lt max.
     * @return the same input int value
     * @throws IllegalArgumentException
     */
    public static int inRange(final int val, final int min, final int max, final String message) {
        Utils.validateArg(val >= min && val <= max, message);
        return val;
    }

    /**
     * Checks that the  input is positive or zero and returns the same value or throws an {@link IllegalArgumentException}
     * @param val value to check
     * @param message the text message that would be pass to the exception thrown
     * @return the same value
     * @throws IllegalArgumentException
     */
    public static double isPositiveOrZero(final double val, final String message) {
        Utils.validateArg(val >= 0, message);
        return val;
    }

    /**
     * Checks that the  input is positive or zero and returns the same value or throws an {@link IllegalArgumentException}
     * @param val value to check
     * @param message the text message that would be pass to the exception thrown
     * @return the same value
     * @throws IllegalArgumentException
     */
    public static double isNegativeOrZero(final double val, final String message) {
        Utils.validateArg(val <= 0, message);
        return val;
    }

    /**
     * Checks that the  input is positive or zero and returns the same value or throws an {@link IllegalArgumentException}
     * @param val value to check
     * @param message the text message that would be pass to the exception thrown
     * @return the same value
     * @throws IllegalArgumentException
     */
    public static long isPositiveOrZero(final long val, final String message) {
        Utils.validateArg(val >= 0, message);
        return val;
    }

    /**
     * Checks that the  input is positive or zero and returns the same value or throws an {@link IllegalArgumentException}
     * @param val value to check
     * @param message the text message that would be pass to the exception thrown
     * @return the same value
     * @throws IllegalArgumentException
     */
    public static int isPositiveOrZero(final int val, final String message) {
        Utils.validateArg(val >= 0, message);
        return val;
    }

    /**
     * Checks that the  input is greater than zero and returns the same value or throws an {@link IllegalArgumentException}
     * @param val value to check
     * @param message the text message that would be pass to the exception thrown
     * @return the same value
     * @throws IllegalArgumentException
     */
    public static long isPositive(final long val, final String message) {
        Utils.validateArg(val > 0, message);
        return val;
    }

    /**
     * Checks that the input {@code int} value is greater than zero and returns the same value or
     * throws an {@link IllegalArgumentException}
     * @param val value to check
     * @param message the text message that would be pass to the exception thrown
     * @return the same value.
     * @throws IllegalArgumentException if the input value is 0 or less.
     */
    public static int isPositive(final int val, final String message) {
        Utils.validateArg(val > 0, message);
        return val;
    }

    /**
     * Checks that the  input is greater than zero and returns the same value or throws an {@link IllegalArgumentException}
     * @param val value to check
     * @param message the text message that would be pass to the exception thrown
     * @return the same value
     * @throws IllegalArgumentException
     */
    public static double isPositive(final double val, final String message) {
        Utils.validateArg(val > 0, message);
        return val;
    }

    /**
     * Checks that the  input is not infinity nor NaN or throws an {@link IllegalArgumentException}
     * @param val value to check
     * @param message the text message that would be pass to the exception thrown
     * @return the same value
     * @throws IllegalArgumentException
     */
    public static double isFinite(final double val, final String message) {
        try {
            MathUtils.checkFinite(val);
        } catch (final NotFiniteNumberException ne) {
            throw new IllegalArgumentException(message);
        }
        return val;
    }

    /**
     * log the value in base b
     *
     * @param val -- The value to be logged.
     * @param base -- The base of the logging
     * @return the logged value.  If val is zero, returns Double.NEGATIVE_INFINITY.  If base is zero, returns -0.0.  If any input
     *  is negative, returns Double.NaN.
     */
    public static double logb(final double val, final double base) {
        return Math.log(val)/Math.log(base);
    }

    /**
     * log the value in base 2
     *
     * @param val -- The value to be logged.
     * @return the logged value.  If val is zero, returns Double.NEGATIVE_INFINITY.  If val
     *  is negative, returns Double.NaN.
     */
    public static double log2(final double val) {
        return Math.log(val) * INV_LOG_2;
    }

    /**
     *  Writes a list of strings to the given filename.  The file will be overwritten if it exists.
     *
     * @param stringIterator Strings that will be written (one to a line)
     * @param outputFile File to write
     */
    public static void writeStringListToFile(final Iterable<String> stringIterator, final File outputFile) {
        Utils.nonNull(stringIterator, "String list cannot be null");
        try (final FileWriter writer = new FileWriter(outputFile.getAbsolutePath())) {
            for (String str : stringIterator) {
                writer.write(str + "\n");
            }
        } catch (final IOException ioe) {
            throw new GATKException("Cannot write to file.", ioe);
        }
    }

    /**
     *  Writes a double array to the given filename.  The file will be overwritten if it exists.
     *
     * @param data doubles that will be written (one to a line)
     * @param outputFile File to write
     */
    public static void writeValuesToFile(final double[] data, final File outputFile) {
        Utils.nonNull(data, "input data cannot be null");
        try (final FileWriter writer = new FileWriter(outputFile.getAbsolutePath())) {
            for (final double d : data) {
                writer.write(String.valueOf(d) + System.lineSeparator());
            }
        } catch (final IOException ioe) {
            throw new GATKException("Cannot write to file.", ioe);
        }
    }

    /**
     *  Reads a double array from the given filename.
     *
     * @param inputFile doubles that will be read (one to a line)
     */
    public static double[] readValuesFromFile(final File inputFile) {
        Utils.nonNull(inputFile, "input file cannot be null");
        try (final BufferedReader reader = new BufferedReader(new FileReader(inputFile.getAbsolutePath()))) {
            return reader.lines().mapToDouble(Double::parseDouble).toArray();
        } catch (final IOException ioe) {
            throw new GATKException("Cannot write to file.", ioe);
        }
    }

    /**
     * Validates the value of a parameter.
     * <p>
     * An invalid value will result in an {@link IllegalArgumentException}.
     * </p>
     * <p>
     *   Notice that {@link Double#NaN nan} are always considered invalid whereas negative or infinity values
     *   will depends on the declared input range extreme values.
     * </p>
     *
     * @param validRange the valid range for the parameter value.
     * @param value the parameter value itself.
     * @param definition a human friendly description of the parameter to be used in an explanatory exception.
     * @return the input value.
     * @throws IllegalArgumentException if the value provided is in-valid.
     */
    public static double inRange(final Range<Double> validRange, final double value, final String definition) {
        Utils.nonNull(validRange);
        Utils.validateArg(validRange.contains(value), String.format("invalid value for %s: %g is not in [%g, %g]",
                definition, value, validRange.getMinimum(), validRange.getMaximum()));
        return value;
    }

    /**
     * Validates the value of a parameter.
     * <p>
     * An invalid value will result in an {@link IllegalArgumentException}.
     * </p>
     *
     * @param validRange the valid range for the parameter value.
     * @param value the parameter value itself.
     * @param definition a human friendly description of the parameter to be used in an explanatory exception.
     * @return the input value.
     * @throws IllegalArgumentException if the value provided is in-valid.
     */
    public static int inRange(final Range<Integer> validRange, final int value, final String definition) {
        Utils.nonNull(validRange);
        if (!validRange.contains(value)) {
            final String prefix = definition == null ? "invalid value" : "invalid value for " + definition;
            throw new IllegalArgumentException(String.format("%s: %d is not in [%d, %d]",
                    prefix, value, validRange.getMinimum(), validRange.getMaximum()));
        }
        return value;
    }

    /**
     * Checks that a range in an array specified using the first element index/offset and the length of the range is valid.
     * @param from  the first position of th range.
     * @param length the length of the range.
     * @param arrayLength the target array length, thus indexes from 0 to arrayLength - 1.
     * @param definition the description of the parameter.
     */
    public static void isValidArrayOffsetAndRangeLength(final int from, final int length,
                                                        final int arrayLength, final String definition) {
        if (length < 0) {
            throw new IllegalArgumentException("range length cannot be negative for " + definition);
        } else if (from < 0) {
            throw new IllegalArgumentException("range offset cannot be negative for " + definition);
        } else if (from + length > arrayLength) {
            throw new IllegalArgumentException("range goes beyond the end of the array for " + definition);
        }
    }

}
