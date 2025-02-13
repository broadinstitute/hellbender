package org.broadinstitute.hellbender.tools.spark.transforms;

import com.google.common.collect.Iterators;
import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.SparkFiles;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.*;
import org.broadinstitute.hellbender.utils.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;

public final class BaseRecalibratorSparkFn {

    /**
     * Run the {@link BaseRecalibrationEngine} on reads and overlapping variants.
     * @param readsWithVariants the RDD of reads with overlapping variants
     * @param header the reads header
     * @param referenceFileName the name of the reference file added via {@code SparkContext#addFile()}
     * @param recalArgs arguments to use during recalibration
     * @return the recalibration report object
     */
    public static RecalibrationReport apply(final JavaPairRDD<GATKRead, Iterable<GATKVariant>> readsWithVariants, final SAMFileHeader header, final String referenceFileName, final RecalibrationArgumentCollection recalArgs) {
        JavaRDD<RecalibrationTables> unmergedTables = readsWithVariants.mapPartitions(readsWithVariantsIterator -> {
            String pathOnExecutor = SparkFiles.get(referenceFileName);
            ReferenceDataSource referenceDataSource = new ReferenceFileSource(IOUtils.getPath(pathOnExecutor));
            final BaseRecalibrationEngine bqsr = new BaseRecalibrationEngine(recalArgs, header);
            bqsr.logCovariatesUsed();
            Utils.stream(readsWithVariantsIterator).forEach(t -> bqsr.processRead(t._1, referenceDataSource, t._2));
            return Iterators.singletonIterator(bqsr.getRecalibrationTables());
        });

        final RecalibrationTables emptyRecalibrationTable = new RecalibrationTables(new StandardCovariateList(recalArgs, header));
        final RecalibrationTables combinedTables = unmergedTables.treeAggregate(emptyRecalibrationTable,
                RecalibrationTables::inPlaceCombine,
                RecalibrationTables::inPlaceCombine,
                Math.max(1, (int)(Math.log(unmergedTables.partitions().size()) / Math.log(2))));

        BaseRecalibrationEngine.collapseQualityScoreTableToReadGroupTable(combinedTables.getQualityScoreTable(), combinedTables.getReadGroupTable());
        BaseRecalibrationEngine.roundTableValues(combinedTables);

        final QuantizationInfo quantizationInfo = new QuantizationInfo(combinedTables, recalArgs.QUANTIZING_LEVELS);

        final StandardCovariateList covariates = new StandardCovariateList(recalArgs, header);
        return RecalUtils.createRecalibrationReport(recalArgs.generateReportTable(covariates.covariateNames()), quantizationInfo.generateReportTable(), RecalUtils.generateReportTables(combinedTables, covariates));
    }
}
