package org.broadinstitute.hellbender.engine.spark.datasources;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SBIIndexWriter;
import htsjdk.samtools.util.FileExtensions;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.HeaderlessSAMRecordCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import org.broadinstitute.hellbender.utils.spark.SparkUtils;
import org.disq_bio.disq.*;

import java.io.IOException;
import java.util.Comparator;

/**
 * ReadsSparkSink writes GATKReads to a file. This code lifts from the HadoopGenomics/Hadoop-BAM
 * read writing code as well as from bigdatagenomics/adam.
 */
public final class ReadsSparkSink {

    private final static Logger logger = LogManager.getLogger(ReadsSparkSink.class);

    /**
     * writeReads writes rddReads to outputFile with header as the file header.
     * @param ctx the JavaSparkContext to write.
     * @param outputFile path to the output bam.
     * @param referencePathSpecifier GATKPath to the reference. required for cram output, otherwise may be null.
     * @param reads reads to write.
     * @param header the header to put at the top of the files
     * @param format should the output be a single file, sharded, ADAM, etc.
     */
    public static void writeReads(
            final JavaSparkContext ctx, final String outputFile, final GATKPath referencePathSpecifier, final JavaRDD<GATKRead> reads,
            final SAMFileHeader header, ReadsWriteFormat format) throws IOException {
        writeReads(ctx, outputFile, referencePathSpecifier, reads, header, format, 0, null, true, SBIIndexWriter.DEFAULT_GRANULARITY);
    }

    /**
     * writeReads writes rddReads to outputFile with header as the file header.
     * @param ctx the JavaSparkContext to write.
     * @param outputFile path to the output bam.
     * @param referencePathSpecifier GATKPath to the reference. required for cram output, otherwise may be null.
     * @param reads reads to write.
     * @param header the header to put at the top of the files
     * @param format should the output be a single file, sharded, ADAM, etc.
     * @param numReducers the number of reducers to use when writing a single file. A value of zero indicates that the default
     *                    should be used.
     * @param outputPartsDir directory for temporary files for SINGLE output format, should be null for default value of filename + .output
     * @param sortReadsToHeader if true, the writer will perform a sort of reads according to the sort order of the header before writing
     * @param splittingIndexGranularity the granularity of the splitting index
     */
    public static void writeReads(
            final JavaSparkContext ctx, final String outputFile, final GATKPath referencePathSpecifier, final JavaRDD<GATKRead> reads,
            final SAMFileHeader header, ReadsWriteFormat format, final int numReducers, final String outputPartsDir, final boolean sortReadsToHeader, final long splittingIndexGranularity) throws IOException {
        writeReads(ctx, outputFile, referencePathSpecifier, reads, header, format, numReducers, outputPartsDir, true, true, sortReadsToHeader, splittingIndexGranularity);
    }

    /**
     * writeReads writes rddReads to outputFile with header as the file header.
     * @param ctx the JavaSparkContext to write.
     * @param outputFile path to the output bam.
     * @param referencePathSpecifier GATKPath to the reference. required for cram output, otherwise may be null.
     * @param reads reads to write.
     * @param header the header to put at the top of the files
     * @param format should the output be a single file, sharded, ADAM, etc.
     * @param numReducers the number of reducers to use when writing a single file. A value of zero indicates that the default
     *                    should be used.
     * @param outputPartsDir directory for temporary files for SINGLE output format, should be null for default value of filename + .output
     * @param writeBai whether to write a BAI file (when writing BAM format)
     * @param writeSbi whether to write an SBI file (when writing BAM format)
     * @param sortReadsToHeader whether to sort the reads in the underlying RDD to match the header sort order option before writing
     * @param splittingIndexGranularity  the granularity of the splitting index if one is created
     */
    public static void writeReads(
            final JavaSparkContext ctx, final String outputFile, final GATKPath referencePathSpecifier, final JavaRDD<GATKRead> reads,
            final SAMFileHeader header, ReadsWriteFormat format, final int numReducers, final String outputPartsDir,
            final boolean writeBai, final boolean writeSbi, final boolean sortReadsToHeader, final long splittingIndexGranularity) throws IOException {

        String absoluteOutputFile = BucketUtils.makeFilePathAbsolute(outputFile);
        ReadsSparkSource.checkCramReference(ctx, new GATKPath(absoluteOutputFile), referencePathSpecifier);

        // The underlying reads are required to be in SAMRecord format in order to be
        // written out, so we convert them to SAMRecord explicitly here. If they're already
        // SAMRecords, this will effectively be a no-op. The SAMRecords will be headerless
        // for efficient serialization.
        final JavaRDD<SAMRecord> samReads = reads.map(read -> read.convertToSAMRecord(null));
        final JavaRDD<SAMRecord> readsToOutput = sortReadsToHeader ? sortSamRecordsToMatchHeader(samReads, header, numReducers) : samReads;

        if (format == ReadsWriteFormat.SINGLE) {
            FileCardinalityWriteOption fileCardinalityWriteOption = FileCardinalityWriteOption.SINGLE;
            final String outputPartsDirectory = (outputPartsDir == null)? getDefaultPartsDirectory(outputFile)  : outputPartsDir;
            TempPartsDirectoryWriteOption tempPartsDirectoryWriteOption = new TempPartsDirectoryWriteOption(outputPartsDirectory);
            BaiWriteOption baiWriteOption = BaiWriteOption.fromBoolean(writeBai);
            SbiWriteOption sbiWriteOption = SbiWriteOption.fromBoolean(writeSbi);
            if (absoluteOutputFile.endsWith(FileExtensions.BAM) ||
                    absoluteOutputFile.endsWith(FileExtensions.CRAM) ||
                    absoluteOutputFile.endsWith(FileExtensions.SAM)) {
                // don't specify a write option for format since it is inferred from the extension in the path
                writeReads(ctx, absoluteOutputFile, referencePathSpecifier, readsToOutput, header,
                        splittingIndexGranularity, fileCardinalityWriteOption, tempPartsDirectoryWriteOption, baiWriteOption, sbiWriteOption);
            } else {
                // default to BAM
                ReadsFormatWriteOption formatWriteOption = ReadsFormatWriteOption.BAM;
                writeReads(ctx, absoluteOutputFile, referencePathSpecifier, readsToOutput, header, splittingIndexGranularity, formatWriteOption,
                        fileCardinalityWriteOption, tempPartsDirectoryWriteOption, baiWriteOption, sbiWriteOption);
            }
        } else if (format == ReadsWriteFormat.SHARDED) {
            if (outputPartsDir!=null) {
                throw new  GATKException(String.format("You specified the bam output parts directory %s, but requested a sharded output format which does not use this option",outputPartsDir));
            }
            ReadsFormatWriteOption formatWriteOption = ReadsFormatWriteOption.BAM; // use BAM if output file is a directory
            FileCardinalityWriteOption fileCardinalityWriteOption = FileCardinalityWriteOption.MULTIPLE;
            writeReads(ctx, absoluteOutputFile, referencePathSpecifier, readsToOutput, header, splittingIndexGranularity, formatWriteOption, fileCardinalityWriteOption);
        }
    }

    private static void writeReads(
            final JavaSparkContext ctx, final String outputFile, final GATKPath referencePathSpecifier, final JavaRDD<SAMRecord> reads,
            final SAMFileHeader header, final long sbiIndexGranularity, final WriteOption... writeOptions) throws IOException {

        Broadcast<SAMFileHeader> headerBroadcast = ctx.broadcast(header);
        final JavaRDD<SAMRecord> sortedReadsWithHeader = reads.map(read -> {
            read.setHeaderStrict(headerBroadcast.getValue());
            return read;
        });
        HtsjdkReadsRdd htsjdkReadsRdd = new HtsjdkReadsRdd(header, sortedReadsWithHeader);
        HtsjdkReadsRddStorage.makeDefault(ctx)
                .referenceSourcePath(referencePathSpecifier == null ? null : referencePathSpecifier.getRawInputString())
                .sbiIndexGranularity(sbiIndexGranularity)
                .write(htsjdkReadsRdd, outputFile, writeOptions);
    }

    private static void deleteHadoopFile(String fileToObliterate, Configuration conf) throws IOException {
        final Path pathToDelete = new Path(fileToObliterate);
        pathToDelete.getFileSystem(conf).delete(pathToDelete, true);
    }

    /**
     * Gets the default parts directory for a given file by appending .parts/ to the end of it
     */
    public static String getDefaultPartsDirectory(String file) {
        return file + ".parts/";
    }

    /**
     * Sorts the given reads according to the sort order in the header.
     * @param reads the reads to sort
     * @param header the header specifying the sort order,
     *               if the header specifies {@link SAMFileHeader.SortOrder#unsorted} or {@link SAMFileHeader.SortOrder#unknown}
     *               then no sort will be performed
     * @param numReducers the number of reducers to use; a value of 0 means use the default number of reducers
     * @return a sorted RDD of reads
     */
    private static JavaRDD<SAMRecord> sortSamRecordsToMatchHeader(final JavaRDD<SAMRecord> reads, final SAMFileHeader header, final int numReducers) {
        final Comparator<SAMRecord> comparator = getSAMRecordComparator(header);
        if ( comparator == null ) {
            return reads;
        } else {
            return SparkUtils.sortUsingElementsAsKeys(reads, comparator, numReducers);
        }
    }

    //Returns the comparator to use or null if no sorting is required.
    private static Comparator<SAMRecord> getSAMRecordComparator(final SAMFileHeader header) {
        switch (header.getSortOrder()){
            case coordinate: return new HeaderlessSAMRecordCoordinateComparator(header);
            //duplicate isn't supported because it doesn't work right on headerless SAMRecords
            case duplicate: throw new UserException.UnimplementedFeature("The sort order \"duplicate\" is not supported in Spark.");
            case queryname:
            case unsorted:   return header.getSortOrder().getComparatorInstance();
            default:         return null; //NOTE: javac warns if you have this (useless) default BUT it errors out if you remove this default.
        }
    }

}
