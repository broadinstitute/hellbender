package org.broadinstitute.hellbender.engine.spark.datasources;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.ReadsPathDataSource;
import org.broadinstitute.hellbender.engine.TraversalParameters;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.*;
import org.broadinstitute.hellbender.utils.spark.SparkUtils;
import org.disq_bio.disq.HtsjdkReadsRdd;
import org.disq_bio.disq.HtsjdkReadsRddStorage;
import org.disq_bio.disq.HtsjdkReadsTraversalParameters;

import java.io.IOException;
import java.io.Serializable;
import java.nio.file.Files;
import java.util.Collections;
import java.util.List;
import java.util.Objects;

/** Loads the reads from disk either serially (using samReaderFactory) or in parallel using Hadoop-BAM.
 * The parallel code is a modified version of the example writing code from Hadoop-BAM.
 */
public final class ReadsSparkSource implements Serializable {
    private static final long serialVersionUID = 1L;

    private transient final JavaSparkContext ctx;
    private ValidationStringency validationStringency = ReadConstants.DEFAULT_READ_VALIDATION_STRINGENCY;

    public ReadsSparkSource(final JavaSparkContext ctx) { this.ctx = ctx; }

    public ReadsSparkSource(final JavaSparkContext ctx, final ValidationStringency validationStringency)
    {
        this.ctx = ctx;
        this.validationStringency = validationStringency;
    }

    /**
     * Loads Reads using Hadoop-BAM. For local files, readFileName must have the fully-qualified path,
     * i.e., file:///path/to/bam.bam.
     * @param readsPathSpecifier file to load
     * @param referencePathSpecifier GATKPath for reference or null if not available. Reference is required for CRAM files.
     * @param traversalParameters parameters controlling which reads to include. If <code>null</code> then all the reads (both mapped and unmapped) will be returned.
     * @return RDD of (SAMRecord-backed) GATKReads from the file.
     */
    public JavaRDD<GATKRead> getParallelReads(final GATKPath readsPathSpecifier, final GATKPath referencePathSpecifier, final TraversalParameters traversalParameters) {
        return getParallelReads(readsPathSpecifier, referencePathSpecifier, traversalParameters, 0);
    }

    /**
     * Loads Reads using Hadoop-BAM. For local files, bam must have the fully-qualified path,
     * i.e., file:///path/to/bam.bam.
     * @param readsPathSpecifier file to load
     * @param referencePathSpecifier GATKPath for reference or null if not available. Reference is required for CRAM files.
     * @param traversalParameters parameters controlling which reads to include. If <code>null</code> then all the reads (both mapped and unmapped) will be returned.
     * @param splitSize maximum bytes of bam file to read into a single partition, increasing this will result in fewer partitions. A value of zero means
     *                  use the default split size (determined by the Hadoop input format, typically the size of one HDFS block).
     * @return RDD of (SAMRecord-backed) GATKReads from the file.
     */
    public JavaRDD<GATKRead> getParallelReads(final GATKPath readsPathSpecifier, final GATKPath referencePathSpecifier, final TraversalParameters traversalParameters, final long splitSize) {
        return getParallelReads(readsPathSpecifier, referencePathSpecifier, traversalParameters, splitSize, false);
    }

    /**
     * Loads Reads using Hadoop-BAM. For local files, bam must have the fully-qualified path,
     * i.e., file:///path/to/bam.bam.
     * @param readPathSpecifier file to load
     * @param referencePathSpecifier GATKPath for reference or null if not available. Reference is required for CRAM files.
     * @param traversalParameters parameters controlling which reads to include. If <code>null</code> then all the reads (both mapped and unmapped) will be returned.
     * @param splitSize maximum bytes of bam file to read into a single partition, increasing this will result in fewer partitions. A value of zero means
     *                  use the default split size (determined by the Hadoop input format, typically the size of one HDFS block).
     * @param useNio whether to use NIO or the Hadoop filesystem for reading files
     * @return RDD of (SAMRecord-backed) GATKReads from the file.
     */
    public JavaRDD<GATKRead> getParallelReads(final GATKPath readPathSpecifier, final GATKPath referencePathSpecifier, final TraversalParameters traversalParameters, final long splitSize, final boolean useNio) {
        try {
            final GATKPath cramReferencePathSpec = checkCramReference(ctx, readPathSpecifier, referencePathSpecifier);
            HtsjdkReadsTraversalParameters<SimpleInterval> tp = traversalParameters == null ? null :
                    new HtsjdkReadsTraversalParameters<>(traversalParameters.getIntervalsForTraversal(), traversalParameters.traverseUnmappedReads());
            HtsjdkReadsRdd htsjdkReadsRdd = HtsjdkReadsRddStorage.makeDefault(ctx)
                    .useNio(useNio)
                    .splitSize((int) splitSize)
                    .validationStringency(validationStringency)
                    .referenceSourcePath(cramReferencePathSpec == null ? null : cramReferencePathSpec.getRawInputString())
                    .read(readPathSpecifier.getRawInputString(), tp);
            JavaRDD<GATKRead> reads = htsjdkReadsRdd.getReads()
                    .map(read -> (GATKRead) SAMRecordToGATKReadAdapter.headerlessReadAdapter(read))
                    .filter(Objects::nonNull);
            return fixPartitionsIfQueryGrouped(ctx, htsjdkReadsRdd.getHeader(), reads);
        } catch (IOException | IllegalArgumentException e) {
            throw new UserException("Failed to load reads from " + readPathSpecifier.getRawInputString() + "\n Caused by:" + e.getMessage(), e);
        }
    }

    private static JavaRDD<GATKRead> fixPartitionsIfQueryGrouped(JavaSparkContext ctx, SAMFileHeader header, JavaRDD<GATKRead> reads) {
        if( ReadUtils.isReadNameGroupedBam(header)) {
            return SparkUtils.putReadsWithTheSameNameInTheSamePartition(header, reads, ctx);
        } else {
            return reads;
        }
    }

    /**
     * Loads Reads using Hadoop-BAM. For local files, readFileName must have the fully-qualified path,
     * i.e., file:///path/to/bam.bam.
     * @param readsPathSpecifier file to load
     * @param referencePathSpecifier Reference path or null if not available. Reference is required for CRAM files.
     * @return RDD of (SAMRecord-backed) GATKReads from the file.
     */
    public JavaRDD<GATKRead> getParallelReads(final GATKPath readsPathSpecifier, final GATKPath referencePathSpecifier) {
        return getParallelReads(readsPathSpecifier, referencePathSpecifier, 0);
    }

    /**
     * Loads Reads using Hadoop-BAM. For local files, readFileName must have the fully-qualified path,
     * i.e., file:///path/to/bam.bam.
     * @param readsPathSpecifier file to load
     * @param referencePathSpecifier Reference path or null if not available. Reference is required for CRAM files.
     * @param splitSize maximum bytes of bam file to read into a single partition, increasing this will result in fewer partitions. A value of zero means
     *                  use the default split size (determined by the Hadoop input format, typically the size of one HDFS block).
     * @return RDD of (SAMRecord-backed) GATKReads from the file.
     */
    public JavaRDD<GATKRead> getParallelReads(final GATKPath readsPathSpecifier, final GATKPath referencePathSpecifier, int splitSize) {
        return getParallelReads(readsPathSpecifier, referencePathSpecifier, null /* all reads */, splitSize);
    }

    /**
     * Loads the header using Hadoop-BAM.
     * @param filePathSpecifier path to the bam.
     * @param referencePathSpecifier Reference path or null if not available. Reference is required for CRAM files.
     * @return the header for the bam.
     */
    public SAMFileHeader getHeader(final GATKPath filePathSpecifier, final GATKPath referencePathSpecifier) {
        final GATKPath cramReferencePathSpec = checkCramReference(ctx, filePathSpecifier, referencePathSpecifier);

        // GCS case
        if (BucketUtils.isGcsUrl(filePathSpecifier)) {
            final SamReaderFactory factory = SamReaderFactory.makeDefault()
                    .validationStringency(validationStringency)
                    .referenceSequence(cramReferencePathSpec == null ? null : referencePathSpecifier.toPath());
            try (final ReadsDataSource readsDataSource =
                         new ReadsPathDataSource(Collections.singletonList(filePathSpecifier.toPath()), factory)) {
                 return readsDataSource.getHeader();
            }
        }

        // local file or HDFs case
        try {
            return HtsjdkReadsRddStorage.makeDefault(ctx)
                    .validationStringency(validationStringency)
                    .referenceSourcePath(cramReferencePathSpec == null ? null : cramReferencePathSpec.getRawInputString())
                    .read(filePathSpecifier.getRawInputString())
                    .getHeader();
        } catch (IOException | IllegalArgumentException e) {
            throw new UserException("Failed to read bam header from " + filePathSpecifier.getRawInputString() + "\n Caused by:" + e.getMessage(), e);
        }
    }

    /**
     * Check that for CRAM the reference is set to a file that exists and is not 2bit.
     * @return a <code>GATKPath</code> or <code>null</code> if not CRAM
     */
    static GATKPath checkCramReference(final JavaSparkContext ctx, final GATKPath filePathSpecifier, final GATKPath referencePathSpecifier) {
        if (filePathSpecifier.isCram()) {
            if (referencePathSpecifier == null) {
                throw new UserException.MissingReference("A reference is required for CRAM input");
            } else if (ReferenceTwoBitSparkSource.isTwoBit(referencePathSpecifier)) { // htsjdk can't handle 2bit reference files
                throw new UserException("A 2bit file cannot be used as a CRAM file reference");
            } else if (referencePathSpecifier.isHadoopURL()) {
                // For Hadoop file system, use a org.apache.hadoop.fs.Path
                if (!SparkUtils.hadoopPathExists(ctx, referencePathSpecifier.getURI())) {
                    throw new UserException.MissingReference("The specified fasta file (" + referencePathSpecifier + ") does not exist.");
                }
            } else {
                // For local or GCS, use nio Path
                final java.nio.file.Path nioReferencePath = referencePathSpecifier.toPath();
                if (!Files.exists(nioReferencePath)) {
                    throw new UserException.MissingReference("The specified fasta file (" + referencePathSpecifier + ") does not exist.");
                }
            }
            return referencePathSpecifier;
        }
        return null;
    }

    /**
     * Tests if a given SAMRecord overlaps any interval in a collection. This is only used as a fallback option for
     * formats that don't support query-by-interval natively at the Hadoop-BAM layer.
     */
    //TODO: use OverlapDetector, see https://github.com/broadinstitute/gatk/issues/1531
    private static boolean samRecordOverlaps(final SAMRecord record, final TraversalParameters traversalParameters ) {
        if (traversalParameters == null) {
            return true;
        }
        if (traversalParameters.traverseUnmappedReads() && record.getReadUnmappedFlag() && record.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START) {
            return true; // include record if unmapped records should be traversed and record is unmapped
        }
        List<SimpleInterval> intervals = traversalParameters.getIntervalsForTraversal();
        if (intervals == null || intervals.isEmpty()) {
            return false; // no intervals means 'no mapped reads'
        }
        for (SimpleInterval interval : intervals) {
            if (record.getReadUnmappedFlag() && record.getAlignmentStart() != SAMRecord.NO_ALIGNMENT_START) {
                // This follows the behavior of htsjdk's SamReader which states that "an unmapped read will be returned
                // by this call if it has a coordinate for the purpose of sorting that is in the query region".
                int start = record.getAlignmentStart();
                return interval.getStart() <= start && interval.getEnd() >= start;
            } else  if (interval.overlaps(record)) {
                return true;
            }
        }
        return false;
    }
}
