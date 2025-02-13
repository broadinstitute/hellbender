package org.broadinstitute.hellbender.utils.tsv;

import com.opencsv.CSVWriter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.Closeable;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

/**
 * A simple TSV/CSV/XSV writer with support for writing in the cloud with configurable delimiter.
 *
 * The expected use case for this class is that first {@link #setHeaderLine} is called with a list of the column names
 * which will be used to determine the number of columns per line as well as how the header is indexed. Then in order to
 * construct a new line call {@link #getNewLineBuilder} to get a line builder for each line, which then has convienent
 * methods for individually assigning column values based on the header line etc. Once a line is finished being mutated
 * one simply needs to call write() on the line to validate and finalize the line.
 *
 * Header lines are encoded in the same format as each row, a single row of delimeted column titles as the first row in the table.
 *
 * Note: this class is intended for creating XSV files with loosely defined input types. If there exists a well defined object
 * that summarizes your table data points then consider using {@link TableWriter}.
 */
public class SimpleXSVWriter implements Closeable {
    private int expectedNumColumns;
    private Map<String, Integer> headerMap = null;
    private CSVWriter outputWriter;

    // The current incomplete line in the writer.
    private LineBuilder currentLineBuilder = null;

    /**
     * Creates a new table writer given the file and column names.
     *
     * @param path         the destination path. This could be a cloud uri (ex. gs://...)
     * @param separator    separator to use for the XSV file
     * @throws IOException              if one was raised when opening the the destination file for writing.
     */
    public SimpleXSVWriter(final Path path, final char separator) throws IOException {
        this( new OutputStreamWriter(
                Files.newOutputStream(Utils.nonNull(path, "The path cannot be null."))),
        separator);
    }

    /**
     * Creates a new table writer given an initialized writer and column names.
     *
     * @param writer       the destination writer.
     * @param separator    separator to use for the TSV file
     * @throws IOException              if one was raised when opening the the destination file for writing.
     */
    public SimpleXSVWriter(final Writer writer, final char separator) {
        Utils.validate(separator!='\n', "Column separator cannot be a newline character");
        outputWriter = new CSVWriter(writer, separator);
    }

    /**
     * Provides a header line to the XSV output file. Note that this will throw an exception if all header lines
     * are not unique as it attempts to create an index for the provided header lines for convenience when building
     * rows of the XSV.
     *
     * NOTE: This can only be set once, XSV output files are expected to only have a single row as header.
     *
     * @param columns Ordered list of header lines to be built into the XSV
     */
    public void setHeaderLine(final List<String> columns, final boolean displayHeader) {
        if (headerMap != null) {
            throw new GATKException("Cannot modify header line once set");
        }
        if (displayHeader) {
            outputWriter.writeNext(columns.toArray(new String[0]), false);
        }
        expectedNumColumns = columns.size();

        // Create the mapping between header and column
        headerMap = new HashMap<>();
        for (int i = 0; i < columns.size(); i++) {
            Utils.nonNull(columns.get(i), "Provided header had null column at position: " + i);
            if (headerMap.putIfAbsent(columns.get(i), i) != null) {
                throw new GATKException("Column names must be unique, but found a duplicate name: " + columns.get(i));
            }
        }
    }

    private void writeLine(String[] line) {
        outputWriter.writeNext(line, false);
        currentLineBuilder = null;
    }

    /**
     * Builds a new LineBuilder and writes out the previous line if it exists.
     *
     * @return a blank LineBuilder to allow for defining the next line
     */
    public LineBuilder getNewLineBuilder() {
        if (headerMap == null) {
            throw new GATKException("Cannot construct line without first setting the header line");
        }
        if (currentLineBuilder != null) {
            currentLineBuilder.write();
        }
        currentLineBuilder = new LineBuilder(expectedNumColumns);
        return currentLineBuilder;
    }

    /**
     * @param column header line to get index for
     * @return zero based index corresponding to that header string, throws an exception if the headerline doesn't exist
     */
    public Integer getIndexForColumn(String column) {
        Utils.nonNull(headerMap, "Cannot request column index if the header has not been specified");
        Integer index = headerMap.get(column);
        Utils.nonNull(index, "Requested column " + column + " does not exist in the provided header");
        return index;
    }

    @Override
    public void close() throws IOException {
        if (currentLineBuilder != null) {
            currentLineBuilder.write();
        }
        outputWriter.close();
    }

    /**
     * Helper to allow for incremental construction of a body line using either indexes or column headings
     * <p>
     * Calling build() will cause the line to be written out into the underlying CSV writer in its current state. Doing
     * so will result in a validation call where an exception will be thrown if any columns of the current line have
     * not been defined. fill() can be used to provide a default value for undefined columns.
     */
    public class LineBuilder {
        String[] lineToBuild;
        boolean hasBuilt = false;

        LineBuilder(int lineLength) {
            lineToBuild = new String[lineLength];
        }

        /**
         * @param row complete line corresponding to this row of the tsv
         */
        public LineBuilder setRow(final String[] row) {
            checkAlterationAfterWrite();
            Utils.validate(row.length == lineToBuild.length, "Provided line must have the correct number of columns");
            for (int i = 0; i < row.length; i++) {
                lineToBuild[i] = row[i];
            }
            return this;
        }

        /**
         * @param row complete line corresponding to this row of the tsv
         */
        public LineBuilder setRow(final List<String> row) {
            checkAlterationAfterWrite();
            Utils.validate(row.size() == lineToBuild.length, "Provided line must have the correct number of columns");
            for (int i = 0; i < row.size(); i++) {
                lineToBuild[i] = row.get(i);
            }
            return this;
        }

        /**
         * @param index Column index to be set
         * @param value Value to be placed into the line
         */
        public LineBuilder setColumn(final int index, final String value) {
            checkAlterationAfterWrite();
            lineToBuild[index] = value;
            return this;
        }

        /**
         * @param heading Column heading to be set
         * @param value   Value to be placed into the line
         */
        public LineBuilder setColumn(final String heading, final String value) {
            int index = getIndexForColumn(heading);
            return setColumn(index, value);
        }

        /**
         * Fills in every empty column of the pending line with the provided value
         */
        public LineBuilder fill(final String filling) {
            checkAlterationAfterWrite();
            for (int i = 0; i < lineToBuild.length; i++) {
                if (lineToBuild[i] == null) {
                    lineToBuild[i] = filling;
                }
            }
            return this;
        }

        /**
         * Constructs the line and writes it out to the output
         */
        public void write() {
            Utils.validate(Arrays.stream(lineToBuild).noneMatch(Objects::isNull), () -> "Attempted to construct an incomplete line, make sure all columns are filled: " + Arrays.toString(lineToBuild));
            writeLine(lineToBuild);
            hasBuilt = true;
        }

        // Throw an exception if we try to alter an already written out line
        private void checkAlterationAfterWrite() {
            Utils.validate(!hasBuilt, "Cannot make alterations to an already written out CSV line");
        }
    }
}