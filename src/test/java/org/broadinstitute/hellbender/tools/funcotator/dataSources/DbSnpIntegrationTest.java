package org.broadinstitute.hellbender.tools.funcotator.dataSources;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.funcotator.BaseFuncotatorArgumentCollection;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorDataSourceDownloader;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.nio.file.Path;
import java.util.List;

/**
 * Class to hold integration tests for the dbSNP data source.
 */
public class DbSnpIntegrationTest extends CommandLineProgramTest {

    private final Path DB_SNP_HG19_FILE_PATH       = IOUtils.getPath(
            FuncotatorDataSourceDownloader.HG19_SOMATIC_GCLOUD_DATASOURCES_BASEURL + "/"
                    + "dbsnp/hg19/" + "hg19_All_20180423.vcf.gz"
    );
    private final Path DB_SNP_HG19_INDEX_FILE_PATH = IOUtils.getPath(
            FuncotatorDataSourceDownloader.HG19_SOMATIC_GCLOUD_DATASOURCES_BASEURL + "/"
                    + "dbsnp/hg19/" + "hg19_All_20180423.vcf.gz.tbi"
    );

    private final Path DB_SNP_HG38_FILE_PATH       = IOUtils.getPath(
            FuncotatorDataSourceDownloader.HG38_SOMATIC_GCLOUD_DATASOURCES_BASEURL + "/"
                    + "dbsnp/hg38/" + "hg38_All_20180418.vcf.gz"
    );
    private final Path DB_SNP_HG38_INDEX_FILE_PATH = IOUtils.getPath(
            FuncotatorDataSourceDownloader.HG38_SOMATIC_GCLOUD_DATASOURCES_BASEURL + "/"
                    + "dbsnp/hg38/" + "hg38_All_20180418.vcf.gz.tbi"
    );

    @DataProvider
    private Object[][] provideFortestDbSnpDataSourceParsing() {
        return new Object[][] {
                // HG19 tests:
                {
                    BaseFuncotatorArgumentCollection.FuncotatorReferenceVersionHg19,
                    new SimpleInterval("chr1", 10318704, 10318704),
                    "rs746945770",
                    Allele.create("G", true),
                    Allele.create("A")
                },
                {
                    BaseFuncotatorArgumentCollection.FuncotatorReferenceVersionHg19,
                    new SimpleInterval("chrX", 31213723, 31213723),
                    "rs5972332",
                    Allele.create("C", true),
                    Allele.create("T")
                },
                {
                    BaseFuncotatorArgumentCollection.FuncotatorReferenceVersionHg19,
                    new SimpleInterval("chrY", 8551842, 8551842),
                    "rs562075277",
                    Allele.create("G", true),
                    Allele.create("A")
                },
                {
                    BaseFuncotatorArgumentCollection.FuncotatorReferenceVersionHg19,
                    new SimpleInterval("chrM", 5005, 5005),
                    "rs879008075",
                    Allele.create("T", true),
                    Allele.create("C")
                },
                // HG38 tests:
                {
                    BaseFuncotatorArgumentCollection.FuncotatorReferenceVersionHg38,
                    new SimpleInterval("chr1", 84349785, 84349785),
                    "rs17131617",
                    Allele.create("T", true),
                    Allele.create("C")
                },
                {
                    BaseFuncotatorArgumentCollection.FuncotatorReferenceVersionHg38,
                    new SimpleInterval("chrX", 80688070, 80688070),
                    "rs3122407",
                    Allele.create("T", true),
                    Allele.create("C")
                },
                {
                    BaseFuncotatorArgumentCollection.FuncotatorReferenceVersionHg38,
                    new SimpleInterval("chrY", 13355944, 13355944),
                    "rs2032654",
                    Allele.create("A", true),
                    Allele.create("G")
                },
                {
                    BaseFuncotatorArgumentCollection.FuncotatorReferenceVersionHg38,
                    new SimpleInterval("chrM", 5131, 5133),
                    "rs199476116",
                    Allele.create("TAA", true),
                    Allele.create("T")
                }
        };
    }

    @Test(
            groups={"cloud"},
            dataProvider = "provideFortestDbSnpDataSourceParsing"
    )
    public void testDbSnpDataSourceParsing( final String refVersion,
                                            final Locatable interval,
                                            final String expectedID,
                                            final Allele expectedRefAllele,
                                            final Allele expectedAltAllele) {
        // 1 - Get the correct version of dbSNP from the funcotator data sources bucket:
        final Path dbSnpFile = (refVersion.equals(BaseFuncotatorArgumentCollection.FuncotatorReferenceVersionHg38))
                ? DB_SNP_HG38_FILE_PATH
                : DB_SNP_HG19_FILE_PATH;

        // 2 - Create a FeatureDataSource from the dbSNP VCF:
        try ( final FeatureDataSource<VariantContext> dbSnpDataSource = new FeatureDataSource<>(dbSnpFile.toUri().toString()) ) {

            // Do a dummy check here:
            Assert.assertNotNull(dbSnpDataSource);

            // 3 - Attempt to read sites and features from the FeatureDataSource:
            final List<VariantContext> features = dbSnpDataSource.queryAndPrefetch(interval);
            Assert.assertEquals(features.size(), 1);

            final VariantContext dbSnpVariant = features.get(0);
            Assert.assertEquals(dbSnpVariant.getContig(), interval.getContig());
            Assert.assertEquals(dbSnpVariant.getStart(), interval.getStart());
            Assert.assertEquals(dbSnpVariant.getEnd(), interval.getEnd());
            Assert.assertEquals(dbSnpVariant.getID(), expectedID);
            Assert.assertEquals(dbSnpVariant.getAlleles().size(), 2);
            Assert.assertEquals(dbSnpVariant.getAlleles().get(0), expectedRefAllele, "Variant has incorrect ref allele: " + dbSnpVariant.getAlleles().get(0) + " != " + expectedRefAllele + " [" + interval + " in " + dbSnpFile + "]");
            Assert.assertEquals(dbSnpVariant.getAlleles().get(1), expectedAltAllele);
        }
    }
}
