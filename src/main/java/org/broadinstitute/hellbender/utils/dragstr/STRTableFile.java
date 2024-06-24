package org.broadinstitute.hellbender.utils.dragstr;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceDictionaryCodec;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.Lazy;
import htsjdk.samtools.util.LineReader;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.output.NullWriter;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.dragstr.DragstrLocus;
import org.broadinstitute.hellbender.tools.dragstr.DragstrLocusUtils;
import org.broadinstitute.hellbender.utils.BinaryTableReader;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.ZipUtils;
import org.broadinstitute.hellbender.tools.dragstr.STRDecimationTable;

import java.io.*;
import java.nio.file.Files;

/**
 * Class to create and access STR table file contents.
 */
public final class STRTableFile implements AutoCloseable {

    public static final String SITES_FILE_NAME = "sites.bin";
    public static final String SITES_INDEX_FILE_NAME = "sites.idx";
    public static final String SITES_TEXT_FILE_NAME = "sites.txt";
    public static final String SUMMARY_FILE_NAME = "summary.txt";
    public static final String REF_DICTIONARY_FILE_NAME = "reference.dict";
    public static final String DECIMATION_TABLE_FILE_NAME = "decimation.txt";

    private final File dir;
    private boolean closed;
    private SAMSequenceDictionary dictionary;
    private Lazy<DragstrLocusUtils.BinaryTableIndex> lociIndex;

    /**
     * Array with the name of the minimum content in the str-table file zip.
     */
    private static final String[] ESSENTIAL_FILES = { SITES_FILE_NAME, SITES_INDEX_FILE_NAME,
                                                     REF_DICTIONARY_FILE_NAME, DECIMATION_TABLE_FILE_NAME};

    private STRTableFile(final File dir) {
        this.dir = dir;
        final String source = new File(dir, REF_DICTIONARY_FILE_NAME).toString();
        try (final LineReader lineReader = new BufferedLineReader(new FileInputStream(source))) {
            this.dictionary = new SAMSequenceDictionaryCodec(NullWriter.INSTANCE).decode(lineReader, source);
        } catch (FileNotFoundException e) {
            throw new GATKException("cannot read dictionary for str-table-file in " + dir);
        }
        lociIndex = new Lazy<>(() -> DragstrLocusUtils.BinaryTableIndex.load(new File(this.dir, SITES_INDEX_FILE_NAME).toString()));
    }

    /**
     * Return a loci reader for a given interval in the genome.
     * @param interval the target interval.
     *
     * @return never {@code null}, but perhaps an "empty" reader that won't return any locus.
     */
    public BinaryTableReader<DragstrLocus> locusReader(final SimpleInterval interval) {
        checkIsNotClose();
        try {
           return DragstrLocusUtils.binaryReader(new File(dir, SITES_FILE_NAME).toString(), lociIndex.get(),
                    dictionary.getSequenceIndex(interval.getContig()), interval.getStart(), interval.getEnd());
        } catch (final IOException ex) {
            throw new GATKException("problems accessing to " + interval + " in str-table-file " + dir);
        }
    }

    /**
     * Return a loci reader across the whole genome.
     *
     * @return never {@code null}, but perhaps an "empty" reader that won't return any locus.
     */
    public BinaryTableReader<DragstrLocus> locusReader() {
        checkIsNotClose();
        try {
            return DragstrLocusUtils.binaryReader(new File(dir, SITES_FILE_NAME));
        } catch (final IOException ex) {
            throw new GATKException("problems accessing to in str-table-file " + dir);
        }
    }

    /**
     * Instantiates an STR table given the location of the file it is stored in.
     * @param path path to the table zip file.
     * @return never {@code null}.
     */
    public static STRTableFile open(final GATKPath path) {
        File dir;
        try {
            dir = Files.createTempDirectory("STRTableFile").toFile();
        }
        catch ( IOException e ) {
            throw new GATKException("Unable to create temp directory for STRTableFile", e);
        }

        try {
            ZipUtils.unzip(path, dir, ESSENTIAL_FILES);
        } catch (final Exception ex) {
            throw new GATKException("issues loading str-table-file at " + path, ex);
        }
        for (final String essentialFileName : ESSENTIAL_FILES) {
            if (!new File(dir, essentialFileName).exists()) {
                throw new GATKException("missing files in the str-table zip: e.g. " + essentialFileName);
            }
        }
        return new STRTableFile(dir);
    }

    /**
     * @throws GATKException if any issue arises during closure.
     */
    @Override
    public void close() {
        if (!closed) {
            closed = true;
            try {
                FileUtils.deleteDirectory(dir);
            } catch (final IOException e) {
                throw new GATKException("issues closing the str-table-file at " + dir, e);
            }
        }
    }

    private void checkIsNotClose() {
        if (closed) {
            throw new IllegalStateException("already closed");
        }
    }

    /**
     * Returns a copy of the decimation table in this STR table file.
     * @return never {@code null}.
     */
    public STRDecimationTable decimationTable() {
        checkIsNotClose();
        return new STRDecimationTable(new File(dir, DECIMATION_TABLE_FILE_NAME).toString());
    }

    /**
     * Returns the dictionary of the reference this str-table is based on.
     * @return never {@code null}.
     */
    public SAMSequenceDictionary dictionary() {
        checkIsNotClose(); // in theory we could return the dictionary anyway since is preloaded but lest fail if
                           // closed to be consistent with other accessors.
        return dictionary;
    }

}
