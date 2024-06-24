package org.broadinstitute.hellbender.utils.read;

/**
 * Possible output formats when writing reads.
 *
 * (Currently used only in the context of Spark.)
 */
public enum ReadsWriteFormat {
    /**
     * Write reads to a single BAM file
     */
    SINGLE,

    /**
     * Write reads to a sharded set of BAM files
     */
    SHARDED
}
