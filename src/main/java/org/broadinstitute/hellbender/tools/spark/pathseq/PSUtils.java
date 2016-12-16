package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

/**
 * Common functions for PathSeq
 */
public final class PSUtils {

    public static JavaRDD<GATKRead> primaryReads(final JavaRDD<GATKRead> reads) {
        return reads.filter(read -> !(read.isSecondaryAlignment() || read.isSupplementaryAlignment()));
    }

    /**
     * Groups pairs using read names
     */
    public static JavaRDD<Iterable<GATKRead>> groupReadPairs(final JavaRDD<GATKRead> reads) {
        return reads.groupBy(GATKRead::getName).values();
    }

    /**
     * Gets flattened RDD of reads that have a mate
     */
    public static JavaRDD<GATKRead> pairedReads(final JavaRDD<GATKRead> reads) {

        return groupReadPairs(reads).filter(p -> {
            final Iterator<GATKRead> itr = p.iterator();
            itr.next();
            return itr.hasNext();
        }).flatMap(Iterable::iterator);
    }

    /**
     * Gets flattened RDD of reads that do not have a mate
     */
    protected static JavaRDD<GATKRead> unpairedReads(final JavaRDD<GATKRead> reads) {
        return PSUtils.groupReadPairs(reads).filter(p -> {
            final Iterator<GATKRead> itr = p.iterator();
            itr.next();
            return !itr.hasNext();
        }).flatMap(Iterable::iterator);
    }

    public static String[] parseCommaDelimitedArgList(final String arg) {
        if (arg == null || arg.isEmpty()) {
            return new String[0];
        }
        return arg.split(",");
    }

    /**
     * Parses command-line option for specifying kmer spacing masks
     */
    public static byte[] parseMask(final String maskArg, final int kSize) {

        final String[] kmerMaskString = parseCommaDelimitedArgList(maskArg);
        final byte[] kmerMask = new byte[kmerMaskString.length];
        Utils.validateArg((kmerMaskString.length & 1) == 0, "Needed kmer mask index list of even length, but found " + kmerMaskString.length);
        for (int i = 0; i < kmerMaskString.length; i++) {
            kmerMask[i] = (byte) Integer.parseInt(kmerMaskString[i]);
            Utils.validateArg(kmerMask[i] >= 0 && kmerMask[i] < kSize, "Invalid kmer mask index: " + kmerMaskString[i]);
        }
        return kmerMask;
    }

    protected static void addReferenceSequencesToHeader( final SAMFileHeader header,
                                                         final String referencePath,
                                                         final SerializableFunction<GATKRead, SimpleInterval> windowFunction,
                                                         final PipelineOptions options) {
        final List<SAMSequenceRecord> refSeqs = PSUtils.getReferenceSequences(referencePath,
                windowFunction, options);
        for (final SAMSequenceRecord rec : refSeqs) {
            if (header.getSequence(rec.getSequenceName()) == null) {
                header.addSequence(rec);
            }
        }
    }

    private static List<SAMSequenceRecord> getReferenceSequences(final String referencePath,
                                                                   final SerializableFunction<GATKRead, SimpleInterval> windowFunction,
                                                                   final PipelineOptions options) {
        final ReferenceMultiSource referenceSource = new ReferenceMultiSource(options, referencePath, windowFunction);
        final SAMSequenceDictionary referenceDictionary = referenceSource.getReferenceSequenceDictionary(null);
        if (referenceDictionary == null) {
            throw new UserException.MissingReferenceDictFile(referencePath);
        }
        return referenceDictionary.getSequences();
    }

    /**
     * Prints warning message followed by a list of relevant items
     */
    public static void logItemizedWarning(final Logger logger, final Collection<String> items, final String warning) {
        if (!items.isEmpty()) {
            String str = "";
            for (final String acc : items) str += acc + ", ";
            str = str.substring(0,str.length() - 2);
            logger.warn(warning + " : " + str);
        }
    }

    /**
     * Writes two objects using Kryo to specified local file path.
     * NOTE: using setReferences(false), which must also be set when reading the file. Does not work with nested
     * objects that reference its parent.
     */
    public static void writeKryoTwo(final String filePath, final Object obj1, final Object obj2) {
        try {
            final Kryo kryo = new Kryo();
            kryo.setReferences(false);
            Output output = new Output(new FileOutputStream(filePath));
            kryo.writeClassAndObject(output, obj1);
            kryo.writeClassAndObject(output, obj2);
            output.close();
        } catch (final FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile("Could not serialize objects to file", e);
        }
    }

    /**
     * Reads taxonomy database that has been written using writeKryoTwo()
     */
    @SuppressWarnings("unchecked")
    public static PSTaxonomyDatabase readTaxonomyDatabase(final String filePath, final PipelineOptions options) {
        final Kryo kryo = new Kryo();

        //There is an odd bug where a ClassNotFound exception gets thrown in Spark mode when Kryo-deserializing the library
        //This has something to do with the default system ClassLoader not loading GATK classes when in Spark mode
        //Working solution: explicitly replace the Kryo system class loader with GATKTool's
        final ClassLoader loader = GATKTool.class.getClassLoader();

        kryo.setClassLoader(loader);
        kryo.setReferences(false);
        final Input input = new Input(BucketUtils.openFile(filePath, options));
        final PSTree tree = (PSTree) kryo.readClassAndObject(input);
        final HashMap<String, String> map = (HashMap<String, String>) kryo.readClassAndObject(input);
        input.close();
        return new PSTaxonomyDatabase(tree, map);
    }

    /**
     * Same as GATKSparkTool's getRecommendedNumReducers(), but can specify input BAM path (for when --input is not used)
     */
    public static int pathseqGetRecommendedNumReducers(final String inputPath, final int numReducers,
                                                 final PipelineOptions options, final int targetPartitionSize) {
        if (numReducers != 0) {
            return numReducers;
        }
        return 1 + (int) (BucketUtils.dirSize(inputPath, options) / targetPartitionSize);
    }
}
