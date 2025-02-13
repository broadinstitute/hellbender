package org.broadinstitute.hellbender.utils.spark;

import com.google.common.collect.Iterators;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.spark.SparkFiles;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.function.PairFlatMapFunction;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.collections.AutoCloseableCollection;
import org.broadinstitute.hellbender.utils.config.ConfigFactory;
import org.broadinstitute.hellbender.utils.iterators.CloseAtEndIterator;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import org.broadinstitute.hellbender.utils.variant.VariantContextVariantAdapter;
import scala.Tuple2;

import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Joins an RDD of GATKReads to variant data by copying the variants files to every node, using Spark's file
 * copying mechanism.
 */
public final class JoinReadsWithVariants {
    private static final int DEFAULT_QUERY_LOOKAHEAD_BASES = 100000;

    private JoinReadsWithVariants() {
    }

    /**
     * Joins each read of an RDD<GATKRead> with overlapping variants from a list of variants files.
     *
     * @param reads the RDD of reads, in coordinate-sorted order
     * @param variantsFileNames the names of the variants files added via {@code SparkContext#addFile()}
     * @return an RDD that contains each read along with the overlapping variants
     */
    public static JavaPairRDD<GATKRead, Iterable<GATKVariant>> join(final JavaRDD<GATKRead> reads, final List<String> variantsFileNames) {
        return reads.mapPartitionsToPair((PairFlatMapFunction<Iterator<GATKRead>, GATKRead, Iterable<GATKVariant>>) gatkReadIterator -> {
            List<FeatureDataSource<VariantContext>> variantSources = variantsFileNames.stream().map(fileName -> openFeatureSource(SparkFiles.get(fileName))).collect(Collectors.toList());
            Iterator<Tuple2<GATKRead, Iterable<GATKVariant>>> iterator = Iterators.transform(gatkReadIterator, read -> getVariantsOverlappingRead(read, variantSources));
            return new CloseAtEndIterator<>(iterator, new AutoCloseableCollection<>(variantSources)); // close FeatureDataSource at end of iteration
        });
    }

    private static Tuple2<GATKRead, Iterable<GATKVariant>> getVariantsOverlappingRead(final GATKRead read, final List<FeatureDataSource<VariantContext>> variantSources) {
        if (SimpleInterval.isValid(read.getContig(), read.getStart(), read.getEnd())) {
            return new Tuple2<>(read, getVariantsOverlappingInterval(variantSources, new SimpleInterval(read)));
        } else {
            //Sometimes we have reads that do not form valid intervals (reads that do not consume any ref bases, eg CIGAR 61S90I
            //In those cases, we'll just say that nothing overlaps the read
            return new Tuple2<>(read, Collections.emptyList());
        }
    }

    private static FeatureDataSource<VariantContext> openFeatureSource(String path) {
        int cloudPrefetchBuffer = ConfigFactory.getInstance().getGATKConfig().cloudPrefetchBuffer();
        int cloudIndexPrefetchBuffer = ConfigFactory.getInstance().getGATKConfig().cloudIndexPrefetchBuffer();
        return new FeatureDataSource<>(path, null, DEFAULT_QUERY_LOOKAHEAD_BASES, VariantContext.class, cloudPrefetchBuffer, cloudIndexPrefetchBuffer);
    }

    private static List<GATKVariant> getVariantsOverlappingInterval(FeatureDataSource<VariantContext> variantSource, SimpleInterval interval) {
        return Utils.stream(variantSource.query(interval)).map(VariantContextVariantAdapter::sparkVariantAdapter).collect(Collectors.toList());
    }

    private static List<GATKVariant> getVariantsOverlappingInterval(List<FeatureDataSource<VariantContext>> variantSources, SimpleInterval interval) {
        return variantSources.stream().map(variantSource -> getVariantsOverlappingInterval(variantSource, interval)).flatMap(List::stream).collect(Collectors.toList());
    }
}