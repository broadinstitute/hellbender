package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.SVInterval;
import org.broadinstitute.hellbender.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.util.ArrayList;

public class CNVInputReader {

    /**
     * Loads an external cnv call list and returns the results in an SVIntervalTree. NB: the contig indices in the SVIntervalTree
     * are based on the sequence indices in the SAM header, _NOT_ the ReadMetadata (which we might not have access to at this
     * time).
     */
    public static SVIntervalTree<VariantContext> loadCNVCalls( final String cnvCallsFile,
                                                               final SAMFileHeader headerForReads) {
        Utils.validate(cnvCallsFile != null, "Can't load null CNV calls file");
        try ( final FeatureDataSource<VariantContext> dataSource = new FeatureDataSource<>(cnvCallsFile, null, 0, VariantContext.class) ) {

            final String sampleId = SVUtils.getSampleId(headerForReads);
            validateCNVcallDataSource(headerForReads, sampleId, dataSource);

            final SVIntervalTree<VariantContext> cnvCallTree = new SVIntervalTree<>();
            Utils.stream(dataSource.iterator())
                    .map(vc -> new VariantContextBuilder(vc).genotypes(vc.getGenotype(sampleId)).make()) // forces a decode of the genotype for serialization purposes
                    .map(vc -> new Tuple2<>(new SVInterval(headerForReads.getSequenceIndex(vc.getContig()), vc.getStart(), vc.getEnd()),vc))
                    .forEach(pv -> cnvCallTree.put(pv._1(), pv._2()));
            return cnvCallTree;
        }
    }

    private static void validateCNVcallDataSource(final SAMFileHeader headerForReads,
                                                  final String sampleId,
                                                  final FeatureDataSource<VariantContext> dataSource) {
        final VCFHeader cnvCallHeader = (VCFHeader) dataSource.getHeader();
        final ArrayList<String> sampleNamesInOrder = cnvCallHeader.getSampleNamesInOrder();
        Utils.validate(sampleNamesInOrder.size() == 1, "CNV call VCF should be single sample");
        Utils.validate(sampleNamesInOrder.contains(sampleId), ("CNV call VCF does not contain calls for sample " + sampleId));
        Utils.validate(cnvCallHeader.getSequenceDictionary() != null,
                "CNV calls file does not have a valid sequence dictionary");
        Utils.validate(cnvCallHeader.getSequenceDictionary().isSameDictionary(headerForReads.getSequenceDictionary()),
                "CNV calls file does not have the same sequence dictionary as the read evidence");
    }
}
