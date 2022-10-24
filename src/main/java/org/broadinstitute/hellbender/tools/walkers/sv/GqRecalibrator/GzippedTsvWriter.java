package org.broadinstitute.hellbender.tools.walkers.sv.GqRecalibrator;

import net.minidev.json.JSONObject;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;
import java.util.zip.GZIPOutputStream;

public class GzippedTsvWriter implements AutoCloseable {
    private boolean open;
    private final String name;
    private final OutputStream outputStream;
    private String typeName = null;
    private int numRows;
    private int numSamples;
    static final int defaultBufferSize = 1000000;
    private final Map<Set<String>, Integer> stringSetEncoderMap = new HashMap<>();
    private final Map<String, Integer> stringEncoderMap = new HashMap<>();
    final static String stringType = "String";
    final static String stringSetType = "StringSet";
    final static String booleanType = "boolean";
    final static String byteType = "byte";
    final static String shortType = "short";
    final static String intType = "int";
    final static String floatType = "float";
    final static String doubleType = "double";

    GzippedTsvWriter(final String variableName, final Path outputFolder) {
        this(variableName, outputFolder, defaultBufferSize);
    }

    GzippedTsvWriter(final String variableName, final Path outputFolder, final int bufferSize) {
        this.name = variableName;
        final Path tsvPath = outputFolder.resolve(name + ".tsv.gz");
        outputStream = getGZipOutputStream(tsvPath, bufferSize);
        open = true;
    }

    static OutputStream getGZipOutputStream(final Path outputPath, final int bufferSize) {
        try {
            return new BufferedOutputStream(
                new GZIPOutputStream(
                    Files.newOutputStream(outputPath), bufferSize, true
                ),
                bufferSize
            );
        } catch (IOException ioException) {
            throw new RuntimeException("Unable to open a gzipped output stream at " + outputPath, ioException);
        }
    }

    @Override
    public void close() {
        if(open) {
            open = false;
            try {
                outputStream.close();
            } catch (IOException ioException) {
                throw new RuntimeException("Error closing " + name, ioException);
            }
        }
    }

    int getNumRows() { return numRows; }

    int getNumSamples() { return numSamples; }

    String getName() { return name; }

    String getTypeName() { return typeName; }

    boolean hasLabels() {
        return (typeName.equals(stringType) && !stringEncoderMap.isEmpty()) ||
            (typeName.equals(stringSetType) && !stringSetEncoderMap.isEmpty());
    }

    static void savePropertiesSummaryJson(Collection<GzippedTsvWriter> tsvWriters, final Path jsonPath) {
        try(final OutputStream outputStream = Files.newOutputStream(jsonPath)) {
            outputStream.write(getPropertiesSummaryJson(tsvWriters).toString().getBytes());
        } catch(IOException ioException) {
            throw new RuntimeException("Error writing GzippedTsvWriter properties summary", ioException);
        }
    }

    static JSONObject getPropertiesSummaryJson(Collection<GzippedTsvWriter> tsvWriters) {
        final JSONObject summariesObject = new JSONObject();
        for(final GzippedTsvWriter tsvWriter : tsvWriters) {
            summariesObject.put(tsvWriter.getName(), tsvWriter.getPropertySummaryJson());
        }
        return summariesObject;
    }

    JSONObject getPropertySummaryJson() {
        final JSONObject summaryObject = new JSONObject();
        summaryObject.put("type", typeName);
        summaryObject.put("num_rows", numRows);
        summaryObject.put("num_samples", numSamples);
        if(hasLabels()) {
            summaryObject.put("codes", getEncodings());
        }
        return summaryObject;
    }

    List<String> getEncodings() {
        switch (typeName) {
            case stringType:
                return stringEncoderMap
                    .entrySet()
                    .stream()
                    .sorted(Map.Entry.comparingByValue())
                    .map(Map.Entry::getKey)
                    .collect(Collectors.toList());
            case stringSetType:
                return stringSetEncoderMap
                    .entrySet()
                    .stream()
                    .sorted(Map.Entry.comparingByValue())
                    .map(Map.Entry::getKey)
                    .map(stringSet -> stringSet.stream().sorted().collect(Collectors.joining(",")))
                    .collect(Collectors.toList());
            default:
                return null;
        }
    }

    List<String> getAllLabels() {
        switch (typeName) {
            case stringType:
                return stringEncoderMap.keySet().stream().sorted().distinct().collect(Collectors.toList());
            case stringSetType:
                return stringSetEncoderMap.keySet().stream().flatMap(Set::stream).sorted().distinct()
                    .collect(Collectors.toList());
            default:
                return null;
        }
    }

    private void setOrCheckType(final String typeName, final int numSamples)  {
        if(this.typeName == null) {
            // this is the first line of data. Set typename and write header line
            this.typeName = typeName;
            this.numSamples = numSamples;
            numRows = 0;
        } else if(!Objects.equals(typeName, this.typeName)) {
            throw new IllegalArgumentException("Value type changed from " + this.typeName + " to " + typeName);
        } else if(numSamples != this.numSamples) {
            throw new IllegalArgumentException("Number of samples changed from " + this.numSamples + " to " + numSamples);
        }
    }

    private void write(final String string) {
        try {
            outputStream.write(string.getBytes(StandardCharsets.UTF_8));
        } catch (IOException ioException) {
            throw new RuntimeException("Error writing to gzipped TSV for " + name);
        }
    }

    private void write(final int value) {
        write(String.format("%d", value));
    }

    private void newline() {
        try {
            outputStream.write("\n".getBytes(StandardCharsets.UTF_8));
        } catch (IOException ioException) {
            throw new RuntimeException("Error writing to gzipped TSV for " + name);
        }
        ++numRows;
    }

    void append(final String value) {
        append(value, true);
    }

    private int encode(final String value) {
        Integer encoded = stringEncoderMap.getOrDefault(value, null);
        if(encoded == null) {
            encoded = stringEncoderMap.size();
            stringEncoderMap.put(value, encoded);
        }
        return encoded;
    }

    void append(final String value, final boolean ordinalEncode) {
        setOrCheckType(stringType, 0);
        if(ordinalEncode) {
            write(encode(value));
        } else {
            write(value);
        }
        newline();
    }

    void append(final String[] values) {
        setOrCheckType(stringType, values.length);
        if(values.length > 0) {
            write(encode(values[0]));
            for (int i = 1; i < values.length; ++i) {
                write("\t" + encode(values[i]));
            }
        }
        newline();
    }

    private int encode(final Set<String> stringSet) {
        Integer encoded = stringSetEncoderMap.getOrDefault(stringSet, null);
        if(encoded == null) {
            encoded = stringSetEncoderMap.size();
            stringSetEncoderMap.put(stringSet, encoded);
        }
        return encoded;
    }

    void append(final Set<String> value) {
        setOrCheckType(stringSetType, 0);
        write(encode(value));
        newline();
    }

    void append(final Set<String>[] values) {
        setOrCheckType(stringSetType, values.length);
        if(values.length > 0) {
            write(encode(values[0]));
            for(int i = 1; i < values.length; ++i) {
                write("\t" + encode(values[i]));
            }
        }
        newline();
    }

    void append(final boolean value) {
        setOrCheckType(booleanType, 0);
        write(value ? "true" : "false");
        newline();
    }

    void append(final boolean[] values) {
        setOrCheckType(booleanType, values.length);
        if(values.length > 0) {
            write(values[0] ? "true" : "false");
            for(int i = 1; i < values.length; ++i) {
                write(values[0] ? "\ttrue" : "\tfalse");
            }
        }
        newline();
    }

    void append(final byte value) {
        setOrCheckType(byteType, 0);
        write((int)value);
        newline();
    }

    void append(final byte[] values) {
        setOrCheckType(byteType, values.length);
        if(values.length > 0) {
            write((int)values[0]);
            for(int i = 1; i < values.length; ++i) {
                write("\t" + (int)values[i]);
            }
        }
        newline();
    }

    void append(final short value) {
        setOrCheckType(shortType, 0);
        write(value);
        newline();
    }

    void append(final short[] values) {
        setOrCheckType(shortType, values.length);
        if(values.length > 0) {
            write(values[0]);
            for(int i = 1; i < values.length; ++i) {
                write("\t" + values[i]);
            }
        }
        newline();
    }

    void append(final int value) {
        setOrCheckType(intType, 0);
        write(value);
        newline();
    }

    void append(final int[] values) {
        setOrCheckType(intType, values.length);
        if(values.length > 0) {
            write(values[0]);
            for(int i = 1; i < values.length; ++i) {
                write("\t" + values[i]);
            }
        }
        newline();
    }

    void append(final float value) {
        setOrCheckType(floatType, 0);
        write(String.format("%.7f", value));
        newline();
    }

    void append(final float[] values) {
        setOrCheckType(floatType, values.length);
        if(values.length > 0) {
            write(String.format("%.7f", values[0]));
            for(int i = 1; i < values.length; ++i) {
                write(String.format("\t%.7f", values[i]));
            }
        }
        newline();
    }

    void append(final double value) {
        setOrCheckType(doubleType, 0);
        write(String.format("%.16f", value));
        newline();
    }

    void append(final double[] values) {
        setOrCheckType(doubleType, values.length);
        if(values.length > 0) {
            write(String.format("%.16f", values[0]));
            for(int i = 1; i < values.length; ++i) {
                write(String.format("\t%.16f", values[i]));
            }
        }
        newline();
    }
}
