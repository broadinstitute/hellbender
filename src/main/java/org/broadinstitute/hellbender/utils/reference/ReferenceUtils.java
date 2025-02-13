package org.broadinstitute.hellbender.utils.reference;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.Iterator;

/**
 * A collection of static methods for dealing with references.
 */
public final class ReferenceUtils {

    // Private so that no one will instantiate this class.
    private ReferenceUtils() {}

    /**
     * Given a fasta filename, return the name of the corresponding index file.
     * (This also works if the file is in gs://)
     */
    public static String getFastaIndexFileName(String fastaFilename) {
        return fastaFilename + ".fai";
    }

    /**
     * Given a fasta filename, return the name of the corresponding dictionary file.
     * (This also works if the file is in gs://)
     */
    public static String getFastaDictionaryFileName(String fastaFilename) {
        int lastDot = fastaFilename.lastIndexOf('.');
        return fastaFilename.substring(0, lastDot) + ".dict";
    }

    /**
     * Given a fasta dictionary file, returns its sequence dictionary
     *
     * @param fastaDictionaryFile fasta dictionary file
     * @return the SAMSequenceDictionary from fastaDictionaryFile
     */
    public static SAMSequenceDictionary loadFastaDictionary( final GATKPath fastaDictionaryFile ) {
        try ( final InputStream fastaDictionaryStream = fastaDictionaryFile.getInputStream() ) {
            return loadFastaDictionary(fastaDictionaryStream);
        }
        catch ( IOException e ) {
            throw new UserException.CouldNotReadInputFile("Error loading fasta dictionary file " + fastaDictionaryFile, e);
        }
        catch ( UserException.MalformedFile e ) {
            throw new UserException.MalformedFile(
                    "Could not read sequence dictionary from given fasta file " +
                            fastaDictionaryFile
            );
        }
    }

    /**
     * Given a fasta dictionary file, returns its sequence dictionary
     *
     * @param fastaDictionaryFile fasta dictionary file
     * @return the SAMSequenceDictionary from fastaDictionaryFile
     */
    public static SAMSequenceDictionary loadFastaDictionary( final File fastaDictionaryFile ) {
        return loadFastaDictionary(new GATKPath(fastaDictionaryFile.getAbsolutePath()));
    }

    /**
     * Given an InputStream connected to a fasta dictionary, returns its sequence dictionary
     *
     * Note: does not close the InputStream it's passed
     *
     * @param fastaDictionaryStream InputStream connected to a fasta dictionary
     * @return the SAMSequenceDictionary from the fastaDictionaryStream
     */
    public static SAMSequenceDictionary loadFastaDictionary( final InputStream fastaDictionaryStream ) {
        // Don't close the reader when we're done, since we don't want to close the client's InputStream for them
        final BufferedLineReader reader = new BufferedLineReader(fastaDictionaryStream);

        final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
        final SAMFileHeader header = codec.decode(reader, fastaDictionaryStream.toString());

        // Make sure we have a valid sequence dictionary before continuing:
        if (header.getSequenceDictionary() == null || header.getSequenceDictionary().isEmpty()) {
            throw new UserException.MalformedFile (
                    "Could not read sequence dictionary from given fasta stream " +
                            fastaDictionaryStream
            );
        }

        return header.getSequenceDictionary();
    }

    public static CachingIndexedFastaSequenceFile createReferenceReader(final GATKPath referenceInput) {
        // fasta reference reader to supplement the edges of the reference sequence
        return new CachingIndexedFastaSequenceFile(referenceInput.toPath());
    }

    public static byte[] getRefBaseAtPosition(final ReferenceSequenceFile reference, final String contig, final int start) {
        return reference.getSubsequenceAt(contig, start, start).getBases();
    }

    public static byte[] getRefBasesAtPosition(final ReferenceSequenceFile reference, final String contig, final int start, final int length) {
        return reference.getSubsequenceAt(contig, start, start+length-1).getBases();
    }

    /**
     * Given a reference path and a sequence interval, calculates the MD5 for the given sequence.
     *
     * Note: MD5 calculation mimics the MD5 calculation in {@link picard.sam.CreateSequenceDictionary}: allows IUPAC bases
     *       but uppercases all bases.
     *
     * @param referencePath The path to the reference.
     * @param interval The interval of the sequence.
     * @return the sequence's MD5 as a String.
     */
    public static String calculateMD5(GATKPath referencePath, SimpleInterval interval){
        MessageDigest md5;
        try {
            md5 = MessageDigest.getInstance("MD5");
        }
        catch(NoSuchAlgorithmException exception){
            throw new GATKException("Incorrect MessageDigest algorithm specified in calculateMD5()", exception);
        }

        // pass in true as the second argument ReferenceDataSource.of() to prevent modification of IUPAC bases or capitalization
        try(final ReferenceDataSource source = ReferenceDataSource.of(referencePath.toPath(), true)) {
            Iterator<Byte> baseIterator = source.query(interval);
            while (baseIterator.hasNext()) {
                Byte b = baseIterator.next();
                md5.update(StringUtil.toUpperCase(b));
            }

            String hash = new BigInteger(1, md5.digest()).toString(16);
            if (hash.length() != 32) {
                final String zeros = "00000000000000000000000000000000";
                hash = zeros.substring(0, 32 - hash.length()) + hash;
            }
            return hash;
        }
    }
}
