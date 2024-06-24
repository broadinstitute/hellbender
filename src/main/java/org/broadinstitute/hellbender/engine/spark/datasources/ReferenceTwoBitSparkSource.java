package org.broadinstitute.hellbender.engine.spark.datasources;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.reference.TwoBitReference;

import java.io.IOException;
import java.io.Serializable;


/**
 * A ReferenceSource impl that is backed by a .2bit representation of a reference genome.  This loads an entire .2bit
 * file into a byte array that is encapsulated by this object.  This is particularly useful for fast reference queries
 * if the entire reference can fit into memory.
 */
public class ReferenceTwoBitSparkSource implements ReferenceSparkSource, Serializable {
    private static final long serialVersionUID = 1L;

    public static final String TWO_BIT_EXTENSION = ".2bit";

    private final String referenceURL;
    private final TwoBitReference twoBitFile;

    public ReferenceTwoBitSparkSource( GATKPath referencePathSpecifier) throws IOException {
        // It would simplify this class if we could cache the GATKPath, but ReferenceFileSparkSource
        // objects are used as Spark broadcast variables, and caching GATKPath here triggers a known
        // issue during broadcast with the Java 11 GATK build. See https://issues.apache.org/jira/browse/SPARK-26963.
        this.referenceURL = referencePathSpecifier.getRawInputString();
        Utils.validateArg(isTwoBit(referencePathSpecifier), "ReferenceTwoBitSource can only take .2bit files");
        this.twoBitFile = new TwoBitReference(new GATKPath(referenceURL));
    }

    /**
     * Gets the reference bases spanning the requested interval. If the interval ends beyond the end of its
     * contig according to our reference source's dictionary, it will be truncated at the contig end.
     *
     * @param interval query interval
     * @return A ReferenceBases containing the reference bases spanning the requested interval, cropped at the
     *         contig end if necessary
     */
    @Override
    public ReferenceBases getReferenceBases(SimpleInterval interval) throws IOException {
        final SimpleInterval queryInterval = cropIntervalAtContigEnd(interval);
        final ReferenceSequence result = twoBitFile.getReferenceBases(queryInterval);
        return new ReferenceBases(result.getBases(), queryInterval);
    }

    @Override
    public SAMSequenceDictionary getReferenceSequenceDictionary(SAMSequenceDictionary optReadSequenceDictionaryToMatch) throws IOException {
        return twoBitFile.getSequenceDictionary();
    }

    public static boolean isTwoBit(final GATKPath referenceSpecifier) {
        return referenceSpecifier.getURI().getPath().endsWith(TWO_BIT_EXTENSION);
    }

    private SimpleInterval cropIntervalAtContigEnd( final SimpleInterval interval ) {
        // The 2bit query API does not support queries beyond the ends of contigs, so we need
        // to truncate our interval at the contig end if necessary.
        final SAMSequenceRecord contigRecord = twoBitFile.getSequenceDictionary().getSequence(interval.getContig());
        Utils.nonNull(contigRecord, () -> "Contig " + interval.getContig() + " not found in reference dictionary");
        return new SimpleInterval(interval.getContig(), interval.getStart(), Math.min(interval.getEnd(), contigRecord.getSequenceLength()));
    }
}
