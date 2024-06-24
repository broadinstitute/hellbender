package org.broadinstitute.hellbender.tools.walkers.mutect;

import com.google.common.collect.ImmutableSet;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFUtils;
import org.apache.commons.collections4.SetUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.walkers.annotator.AnnotationUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.PileupDetectionArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReadThreadingAssemblerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceMode;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.FilterMutectCalls;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.M2FiltersArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.NuMTFilterTool;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.ReadOrientationFilter;
import org.broadinstitute.hellbender.tools.walkers.readorientation.LearnReadOrientationModel;
import org.broadinstitute.hellbender.tools.walkers.validation.Concordance;
import org.broadinstitute.hellbender.tools.walkers.validation.ConcordanceSummaryRecord;
import org.broadinstitute.hellbender.tools.walkers.variantutils.ValidateVariants;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Created by davidben on 9/1/16.
 */
@Test(groups = {"variantcalling"})
public class Mutect2IntegrationTest extends CommandLineProgramTest {
    // positions 10,000,000 - 11,000,000 of chr 20 and with most annotations removed
    private static final File GNOMAD = new File(largeFileTestDir, "very-small-gnomad.vcf");
    private static final String DREAM_BAMS_DIR = largeFileTestDir + "mutect/dream_synthetic_bams/";
    private static final File DREAM_4_NORMAL = new File(DREAM_BAMS_DIR, "normal_4.bam");
    private static final File DREAM_3_NORMAL = new File(DREAM_BAMS_DIR, "normal_3.bam");
    private static final File DREAM_4_TUMOR = new File(DREAM_BAMS_DIR, "tumor_4.bam");
    private static final File DREAM_3_TUMOR = new File(DREAM_BAMS_DIR, "tumor_3.bam");
    private static final File DREAM_2_NORMAL = new File(DREAM_BAMS_DIR, "normal_2.bam");
    private static final File DREAM_1_NORMAL = new File(DREAM_BAMS_DIR, "normal_1.bam");
    private static final File DREAM_2_TUMOR = new File(DREAM_BAMS_DIR, "tumor_2.bam");
    private static final File DREAM_1_TUMOR = new File(DREAM_BAMS_DIR, "tumor_1.bam");

    private static final String DREAM_VCFS_DIR = toolsTestDir + "mutect/dream/vcfs/";
    private static final File DREAM_4_TRUTH = new File(DREAM_VCFS_DIR, "sample_4.vcf");
    private static final File DREAM_3_TRUTH = new File(DREAM_VCFS_DIR, "sample_3.vcf");
    private static final File DREAM_2_TRUTH = new File(DREAM_VCFS_DIR, "sample_2.vcf");
    private static final File DREAM_1_TRUTH = new File(DREAM_VCFS_DIR, "sample_1.vcf");

    private static final String DREAM_MASKS_DIR = toolsTestDir + "mutect/dream/masks/";
    private static final File DREAM_4_MASK = new File(DREAM_MASKS_DIR, "mask4.list");
    private static final File DREAM_3_MASK = new File(DREAM_MASKS_DIR, "mask3.list");
    private static final File DREAM_2_MASK = new File(DREAM_MASKS_DIR, "mask2.list");
    private static final File DREAM_1_MASK = new File(DREAM_MASKS_DIR, "mask1.list");

    private static final File NO_CONTAMINATION_TABLE = new File(toolsTestDir, "mutect/no-contamination.table");
    private static final File FIVE_PCT_CONTAMINATION_TABLE = new File(toolsTestDir, "mutect/five-pct-contamination.table");
    private static final File TEN_PCT_CONTAMINATION_TABLE = new File(toolsTestDir, "mutect/ten-pct-contamination.table");

    private static final File NA12878_MITO_BAM = new File(toolsTestDir, "mutect/mito/NA12878.bam");
    private static final File NA12878_MITO_VCF = new File(toolsTestDir, "mutect/mito/unfiltered-with-assb.vcf");
    private static final File NA12878_MITO_GVCF = new File(toolsTestDir, "mitochondria/NA12878.MT.g.vcf");
    private static final File NA12878_MITO_INITIAL_FILTERED_VCF = new File(toolsTestDir, "mutect/mito/initialFiltered.vcf");
    private static final File MITO_REF = new File(toolsTestDir, "mutect/mito/Homo_sapiens_assembly38.mt_only.fasta");
    private static final File DEEP_MITO_BAM = new File(largeFileTestDir, "mutect/highDPMTsnippet.bam");
    private static final String DEEP_MITO_SAMPLE_NAME = "mixture";

    private static final File FILTERING_DIR = new File(toolsTestDir, "mutect/filtering");

    private static final File FLOW_BASED_INPUT_BAM = new File(publicTestDir + "/large", "input_jukebox_for_test.bam");

    private static final File GNOMAD_WITHOUT_AF_SNIPPET = new File(toolsTestDir, "mutect/gnomad-without-af.vcf");
    private static final File UNPARSIMONIOUS_GNOMAD_SNIPPET = new File(toolsTestDir, "mutect/unparsimonious_germline.vcf");

    private static final double TLOD_MATCH_EPSILON = 0.05;
    private static final double VARIANT_TLOD_MATCH_PCT = 0.01;

    private static final String CHROMOSOME_20 = "20";

    // tumor bams, normal bams, truth vcf, mask, required sensitivity, whether to error-correct
    @DataProvider(name = "dreamSyntheticData")
    public Object[][] dreamSyntheticData() {
        return new Object[][]{
                {DREAM_1_TUMOR, Optional.of(DREAM_1_NORMAL), DREAM_1_TRUTH, DREAM_1_MASK, 0.98, false},
                {DREAM_2_TUMOR, Optional.of(DREAM_2_NORMAL), DREAM_2_TRUTH, DREAM_2_MASK, 0.98, false},
                {DREAM_2_TUMOR, Optional.empty(), DREAM_2_TRUTH, DREAM_2_MASK, 0.98, false},
                {DREAM_2_TUMOR, Optional.empty(), DREAM_2_TRUTH, DREAM_2_MASK, 0.98, true},
                {DREAM_3_TUMOR, Optional.of(DREAM_3_NORMAL), DREAM_3_TRUTH, DREAM_3_MASK, 0.95, false},
                {DREAM_4_TUMOR, Optional.of(DREAM_4_NORMAL), DREAM_4_TRUTH, DREAM_4_MASK, 0.8, false},
                {DREAM_4_TUMOR, Optional.of(DREAM_4_NORMAL), DREAM_4_TRUTH, DREAM_4_MASK, 0.7, true},
                {DREAM_4_TUMOR, Optional.empty(), DREAM_4_TRUTH, DREAM_4_MASK, 0.7, false},
        };
    }

    /**
     * Several DREAM challenge bams with synthetic truth data.  In order to keep file sizes manageable, bams are restricted
     * to chromosome 20, leaving ~100-200 variants, and then further restricted to 400-bp intervals centered around
     * these variants.
     * <p>
     * Because this truth data is synthetic, ideal performance is perfect concordance with the truth vcfs.
     * <p>
     * The characteristics of the data are as follows (note that we removed structural variants from the original
     * DREAM challenge vcfs):
     * <p>
     * Sample 1: pure monoclonal sample, SNVs only
     * Sample 2: 80% pure monoclonal sample, SNVs only
     * Sample 3: pure triclonal sample, subclone minor allele frequencies are 1/2, 1/3, and 1/5, SNVs and indels
     * Sample 4: 80% biclonal sample, subclone minor allele fractions are 50% and 35%, SNVs and indels
     *
     */
    @Test(dataProvider = "dreamSyntheticData")
    public void testDreamTumorNormal(final File tumor, final Optional<File> normal, final File truth, final File mask,
                                     final double requiredSensitivity, final boolean errorCorrectReads) throws Exception {
        Utils.resetRandomGenerator();
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");
        final File filteredVcf = createTempFile("filtered", ".vcf");
        final File f1r2Counts = createTempFile("f1r2", ".tar.gz");
        final File orientationModel = createTempFile("orientation", ".tar.gz");
        final File dataset = createTempFile("dataset", ".txt");

        final List<File> normals = normal.isPresent() ? Collections.singletonList(normal.get()) : Collections.emptyList();
        runMutect2(Collections.singletonList(tumor), normals, unfilteredVcf, CHROMOSOME_20, b37Reference, Optional.of(GNOMAD),
                args -> args.addMask(mask).add(M2ArgumentCollection.F1R2_TAR_GZ_NAME, f1r2Counts),
                args -> args.add(M2ArgumentCollection.MUTECT3_DATASET_LONG_NAME, dataset),
                args -> args.addFlag(M2ArgumentCollection.MUTECT3_TRAINING_MODE_LONG_NAME),
                args -> errorCorrectReads ? args.add(ReadThreadingAssemblerArgumentCollection.PILEUP_ERROR_CORRECTION_LOG_ODDS_LONG_NAME, 3.0) : args
        );

        // verify that alleles contained in likelihoods matrix but dropped from somatic calls do not show up in annotations
        // also check that alleles have been properly clipped after dropping any non-called alleles, i.e. if we had AAA AA A
        // and A got dropped, we need AAA AA -> AA A.  The condition we don't want is that all alleles share a common first base
        // and no allele has length 1.
        VariantContextTestUtils.streamVcf(unfilteredVcf)
                .forEach(vc -> {
                    for (final Genotype genotype : vc.getGenotypes()) {
                        final int[] f1r2 = ReadOrientationFilter.getF1R2(genotype);
                        Assert.assertEquals(f1r2.length, vc.getNAlleles());
                        if (vc.getAlleles().stream().filter(a -> !a.isSymbolic()).map(a -> a.getBases()[0]).distinct().count() == 1) {
                            Assert.assertTrue(vc.getAlleles().stream().anyMatch(a -> a.getBases().length == 1));
                        }
                    }
                });

        final ArgumentsBuilder orientationBiasArgs = new ArgumentsBuilder().addInput(f1r2Counts).addOutput(orientationModel);
        runCommandLine(orientationBiasArgs, LearnReadOrientationModel.class.getSimpleName());

        for (final boolean runOrientationFilter : new boolean[] { true, false}) {

            runFilterMutectCalls(unfilteredVcf, filteredVcf, b37Reference,
                    args -> runOrientationFilter ? args.add(M2FiltersArgumentCollection.ARTIFACT_PRIOR_TABLE_NAME, orientationModel) : args);

            final File concordanceSummary = createTempFile("concordance", ".txt");
            final File truePositivesFalseNegatives = createTempFile("tpfn", ".vcf");
            runConcordance(truth, filteredVcf,concordanceSummary, CHROMOSOME_20, mask, Optional.of (truePositivesFalseNegatives));

            final List<ConcordanceSummaryRecord> summaryRecords = new ConcordanceSummaryRecord.Reader(concordanceSummary.toPath()).toList();
            summaryRecords.forEach(rec -> {
                if (rec.getTruePositives() + rec.getFalseNegatives() > 0) {
                    Assert.assertTrue(rec.getSensitivity() > requiredSensitivity);
                    // tumor-only will have germline variants sneak in
                    if (!normals.isEmpty()) {
                        Assert.assertTrue(rec.getPrecision() > 0.5);
                    }
                }
            });
        }
    }

    @Test
    public void testNA12878NormalNormalFiltering() {
        Utils.resetRandomGenerator();
        final File unfilteredVcf = new File(FILTERING_DIR, "NA12878.vcf");
        final File contamination = new File(FILTERING_DIR, "contamination.table");
        final File segments = new File(FILTERING_DIR, "segments.table");

        final File filteredVcf = createTempFile("filtered", ".vcf");

        runFilterMutectCalls(unfilteredVcf, filteredVcf, b37Reference,
                args -> args.add(M2FiltersArgumentCollection.TUMOR_SEGMENTATION_LONG_NAME, segments),
                args -> args.add(M2FiltersArgumentCollection.CONTAMINATION_TABLE_LONG_NAME, contamination));

        final long numPassVariants = VariantContextTestUtils.streamVcf(filteredVcf)
                .filter(vc -> vc.getFilters().isEmpty()).count();

        Assert.assertTrue(numPassVariants < 13);
    }

    // tumorBams, normalBam, truthVcf, mask, requiredSensitivity
    @DataProvider(name = "twoTumorData")
    public Object[][] twoTumorData() {
        return new Object[][]{
                {Arrays.asList(DREAM_1_TUMOR, DREAM_2_TUMOR), Collections.singletonList(DREAM_1_NORMAL), DREAM_1_TRUTH, DREAM_1_MASK, 0.97},
                {Arrays.asList(DREAM_3_TUMOR, DREAM_4_TUMOR), Collections.singletonList(DREAM_3_NORMAL), DREAM_3_TRUTH, DREAM_3_MASK, 0.90}
        };
    }

    @Test(dataProvider = "twoTumorData")
    public void testTwoDreamTumorSamples(final List<File> tumors, final List<File> normals,
                                         final File truth, final File mask, final double requiredSensitivity) throws Exception {
        Utils.resetRandomGenerator();
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");
        final File filteredVcf = createTempFile("filtered", ".vcf");

        runMutect2(tumors, normals, unfilteredVcf, CHROMOSOME_20, b37Reference, Optional.of(GNOMAD), args -> args.addMask(mask));
        runFilterMutectCalls(unfilteredVcf, filteredVcf, b37Reference);

        final File concordanceSummary = createTempFile("concordance", ".txt");
        runConcordance(truth, filteredVcf, concordanceSummary, CHROMOSOME_20, mask, Optional.empty());

        final List<ConcordanceSummaryRecord> summaryRecords = new ConcordanceSummaryRecord.Reader(concordanceSummary.toPath()).toList();
        summaryRecords.forEach(rec -> {
            if (rec.getTruePositives() + rec.getFalseNegatives() > 0) {
                Assert.assertTrue(rec.getSensitivity() > requiredSensitivity);
                // tumor-only will have germline variants sneak in
                if (!normals.isEmpty()) {
                    //Assert.assertTrue(rec.getPrecision() > 0.5);
                }
            }
        });
    }

    // make a pon with a tumor and then use this pon to call somatic variants on the same tumor
    // if the pon is doing its job all calls should be filtered by this pon
    @Test
    public void testPon() {
        Utils.resetRandomGenerator();
        final File tumor = DREAM_1_TUMOR;
        final File normal = DREAM_1_NORMAL;
        final File ponVcf = createTempFile("pon", ".vcf");
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");
        final File filteredVcf = createTempFile("filtered", ".vcf");

        runMutect2(tumor, normal, ponVcf, CHROMOSOME_20, b37Reference, Optional.empty());
        runMutect2(tumor, normal, unfilteredVcf, CHROMOSOME_20, b37Reference, Optional.empty(),
                args -> args.add(M2ArgumentCollection.PANEL_OF_NORMALS_LONG_NAME, ponVcf));
        runFilterMutectCalls(unfilteredVcf, filteredVcf, b37Reference);

        final long numVariants = VariantContextTestUtils.streamVcf(filteredVcf)
                .filter(vc -> vc.getFilters().isEmpty()).count();

        Assert.assertEquals(numVariants, 0);
    }

    // run tumor-normal mode using the original DREAM synthetic sample 1 tumor and normal restricted to
    // 1/3 of our dbSNP interval, in which there is only one true positive.
    // we want to see that the number of false positives is small
    @Test
    public void testTumorNormal()  {
        Utils.resetRandomGenerator();
        final File unfilteredVcf = createTempFile("output", ".vcf");
        final File filteredVcf = createTempFile("filtered", ".vcf");
        final List<File> tumor = Collections.singletonList(new File(DREAM_BAMS_DIR, "tumor.bam"));
        final List<File> normals = Arrays.asList(new File(DREAM_BAMS_DIR, "normal.bam"), DREAM_2_NORMAL);

        runMutect2(tumor, normals, unfilteredVcf, "20:10000000-10100000", b37Reference, Optional.empty());
        runFilterMutectCalls(unfilteredVcf, filteredVcf, b37Reference);

        VariantContextTestUtils.streamVcf(unfilteredVcf).flatMap(vc -> vc.getGenotypes().stream()).forEach(g -> Assert.assertTrue(g.hasAD()));
        final long numVariants = VariantContextTestUtils.streamVcf(filteredVcf)
                .filter(VariantContext::isNotFiltered)
                .count();
        Assert.assertTrue(numVariants < 4);
    }

    // run tumor-only using our mini gnomAD on NA12878, which is not a tumor
    @Test
    public void testTumorOnly() {
        Utils.resetRandomGenerator();
        final File tumor = new File(NA12878_20_21_WGS_bam);
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");
        final File filteredVcf = createTempFile("filtered", ".vcf");

        runMutect2(tumor, unfilteredVcf, "20:10000000-10010000", b37Reference, Optional.of(GNOMAD));
        runFilterMutectCalls(unfilteredVcf, filteredVcf, b37Reference);

        final long numVariantsBeforeFiltering = VariantContextTestUtils.streamVcf(unfilteredVcf).count();

        final long numVariantsPassingFilters = VariantContextTestUtils.streamVcf(filteredVcf)
                .filter(vc -> vc.getFilters().isEmpty()).count();

        // just a sanity check that this bam has some germline variants on this interval so that our test doesn't pass trivially!
        Assert.assertTrue(numVariantsBeforeFiltering > 15);

        // every variant on this interval in this sample is in gnomAD
        Assert.assertTrue(numVariantsPassingFilters < 2);
    }

    @Test
    public void testMutect3Dataset() {
        Utils.resetRandomGenerator();
        final File tumor = new File(NA12878_20_21_WGS_bam);
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");
        final File mutect3Dataset = createTempFile("mutect3", ".data");

        runMutect2(tumor, unfilteredVcf, "20:10000000-10010000", b37Reference, Optional.of(GNOMAD),
                   args -> args.addFlag(ReadThreadingAssemblerArgumentCollection.LINKED_DE_BRUIJN_GRAPH_LONG_NAME),
                   args -> args.add(M2ArgumentCollection.MUTECT3_DATASET_LONG_NAME, mutect3Dataset));
    }

    // make sure we can call tumor alts when the normal has a different alt at the same site
    // regression test for https://github.com/broadinstitute/gatk/issues/6901
    @Test
    public void testDifferentAltsInTumorAndNormal() {
        Utils.resetRandomGenerator();
        final File tumor = new File(toolsTestDir, "mutect/different-alts-tumor.bam");
        final File normal = new File(toolsTestDir, "mutect/different-alts-normal.bam");
        final File output = createTempFile("output", ".vcf");


        runMutect2(tumor, normal, output, "20:10020000-10021200", b37Reference, Optional.empty(),
                args -> args.add(M2ArgumentCollection.NORMAL_LOG_10_ODDS_LONG_NAME, 0.0));

        final Map<Integer, Allele> altAllelesByPosition = VariantContextTestUtils.streamVcf(output)
                .collect(Collectors.toMap(VariantContext::getStart, VariantContext::getAltAlleleWithHighestAlleleCount));

        Assert.assertTrue(altAllelesByPosition.get(10020042).basesMatch(Allele.ALT_C)); //tumor G->C, normal G->A
        Assert.assertTrue(altAllelesByPosition.get(10020124).basesMatch(Allele.ALT_G)); //tumor A->G, normal A->T
    }

    // test on an artificial bam with several contrived MNPs
    @Test
    public void testMnps() {
        Utils.resetRandomGenerator();
        final File bam = new File(toolsTestDir, "mnp.bam");

        for (final int maxMnpDistance : new int[]{0, 1, 2, 3, 5}) {
            final File outputVcf = createTempFile("unfiltered", ".vcf");

            runMutect2(bam, outputVcf, "20:10019000-10022000", b37Reference, Optional.empty(),
                    args -> args.add(M2ArgumentCollection.EMISSION_LOG_SHORT_NAME, 15),
                    args -> args.add(AssemblyBasedCallerArgumentCollection.MAX_MNP_DISTANCE_SHORT_NAME, maxMnpDistance));

            // note that for testing HaplotypeCaller GVCF mode we will always have the symbolic <NON REF> allele
            final Map<Integer, List<String>> alleles = VariantContextTestUtils.streamVcf(outputVcf)
                    .collect(Collectors.toMap(VariantContext::getStart, vc -> vc.getAlternateAlleles().stream().filter(a -> !a.isSymbolic()).map(Allele::getBaseString).collect(Collectors.toList())));

            // phased, two bases apart
            if (maxMnpDistance < 2) {
                Assert.assertEquals(alleles.get(10019968), Collections.singletonList("G"));
                Assert.assertEquals(alleles.get(10019970), Collections.singletonList("G"));
            } else {
                Assert.assertEquals(alleles.get(10019968), Collections.singletonList("GAG"));
                Assert.assertTrue(!alleles.containsKey(10019970));
            }

            // adjacent and out of phase
            Assert.assertEquals(alleles.get(10020229), Collections.singletonList("A"));
            Assert.assertEquals(alleles.get(10020230), Collections.singletonList("G"));

            // 4-substitution MNP w/ spacings 2, 3, 4
            if (maxMnpDistance < 2) {
                Assert.assertEquals(alleles.get(10020430), Collections.singletonList("G"));
                Assert.assertEquals(alleles.get(10020432), Collections.singletonList("G"));
                Assert.assertEquals(alleles.get(10020435), Collections.singletonList("G"));
                Assert.assertEquals(alleles.get(10020439), Collections.singletonList("G"));
            } else if (maxMnpDistance < 3) {
                Assert.assertEquals(alleles.get(10020430), Collections.singletonList("GAG"));
                Assert.assertEquals(alleles.get(10020435), Collections.singletonList("G"));
                Assert.assertEquals(alleles.get(10020439), Collections.singletonList("G"));
            } else if (maxMnpDistance < 4) {
                Assert.assertEquals(alleles.get(10020430), Collections.singletonList("GAGTTG"));
                Assert.assertEquals(alleles.get(10020439), Collections.singletonList("G"));
            } else {
                Assert.assertEquals(alleles.get(10020430), Collections.singletonList("GAGTTGTCTG"));
            }

            // two out of phase DNPs that overlap and have a base in common
            if (maxMnpDistance > 0) {
                Assert.assertEquals(alleles.get(10020680), Collections.singletonList("TA"));
                Assert.assertEquals(alleles.get(10020681), Collections.singletonList("AT"));
            }
        }
    }

    @Test
    public void testForceCalling() {
        Utils.resetRandomGenerator();
        final File tumor = new File(NA12878_20_21_WGS_bam);

        // The kmerSize = 1 case is a ridiculous setting that forces assembly to fail.
        for (final int kmerSize : new int[] {1, 20}) {
            final File unfilteredVcf = createTempFile("unfiltered", ".vcf");
            final File forceCalls = new File(toolsTestDir, "mutect/gga_mode.vcf");

            runMutect2(tumor, unfilteredVcf, "20:9998500-10010000", b37Reference, Optional.empty(),
                    args -> args.add(AssemblyBasedCallerArgumentCollection.FORCE_CALL_ALLELES_LONG_NAME, forceCalls),
                    args -> args.add(ReadThreadingAssemblerArgumentCollection.KMER_SIZE_LONG_NAME, kmerSize),
                    args -> args.add(ReadThreadingAssemblerArgumentCollection.DONT_INCREASE_KMER_SIZE_LONG_NAME, true));

            final Map<Integer, List<Allele>> altAllelesByPosition = VariantContextTestUtils.streamVcf(unfilteredVcf)
                    .collect(Collectors.toMap(VariantContext::getStart, VariantContext::getAlternateAlleles));
            for (final VariantContext vc : new FeatureDataSource<VariantContext>(forceCalls)) {
                final List<Allele> altAllelesAtThisLocus = altAllelesByPosition.get(vc.getStart());
                vc.getAlternateAlleles().stream().filter(a-> a.length() > 0 && BaseUtils.isNucleotide(a.getBases()[0])).forEach(a -> Assert.assertTrue(altAllelesAtThisLocus.contains(a)));
            }
        }
    }

    // regression test for PR: https://github.com/broadinstitute/gatk/pull/8717, which fixed an issue wherein germline events were
    // incorrectly contributing to the bad haplotype and clustered events filters.  In order to make the test more stringent we turn
    // on the genotype-germline-sites flag
    //  The test is on two variants that caused particular trouble previously in the DREAM 2 sample
    @Test
    public void testFilteredHaplotypeAndClusteredEventsFilters() throws Exception {
        Utils.resetRandomGenerator();
        final File tumor = DREAM_2_TUMOR;
        final File normal = DREAM_2_NORMAL;
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");
        final File filteredVcf = createTempFile("filtered", ".vcf");

        for (final int locus : List.of(4567628, 20870771)) {
            final String interval = "20:" + (locus - 500) + "-" + (locus + 500);
            runMutect2(List.of(tumor), List.of(normal), unfilteredVcf, interval, b37Reference, Optional.of(GNOMAD),
                    args -> args.addFlag(M2ArgumentCollection.GENOTYPE_GERMLINE_SITES_LONG_NAME));
            runFilterMutectCalls(unfilteredVcf, filteredVcf, b37Reference);
            final Optional<VariantContext> vcShouldPass = VariantContextTestUtils.streamVcf(filteredVcf).filter(vc -> vc.getStart() == locus).findFirst();
            Assert.assertTrue(vcShouldPass.isPresent());
            Assert.assertFalse(vcShouldPass.get().isFiltered());
        }
    }

    // test that the dont-use-soft-clips option actually does something
    @Test
    public void testDontUseSoftClips() {
        Utils.resetRandomGenerator();
        final File tumor = new File(NA12878_20_21_WGS_bam);
        final int start = 10050000;

        final SimpleInterval interval = new SimpleInterval("20", start, start + 25000);

        final File calls1 = createTempFile("unfiltered", ".vcf");
        runMutect2(tumor, calls1, interval.toString(), b37Reference, Optional.of(GNOMAD));

        final File calls2 = createTempFile("unfiltered", ".vcf");
        runMutect2(tumor, calls2, interval.toString(), b37Reference, Optional.of(GNOMAD),
                args -> args.addFlag(AssemblyBasedCallerArgumentCollection.DONT_USE_SOFT_CLIPPED_BASES_LONG_NAME));

        final List<VariantContext> indelsWithSoftClips = VariantContextTestUtils.streamVcf(calls1).filter(VariantContext::isIndel).collect(Collectors.toList());
        final List<VariantContext> indelsWithoutSoftClips = VariantContextTestUtils.streamVcf(calls2).filter(VariantContext::isIndel).collect(Collectors.toList());

        Assert.assertTrue(indelsWithoutSoftClips.size() < indelsWithSoftClips.size());

        final int startOfDroppedVariant = 10068160;
        final int endOfDroppedVariant = 10068174;
        Assert.assertTrue(indelsWithSoftClips.stream().anyMatch(vc -> vc.getStart() == startOfDroppedVariant && vc.getEnd() == endOfDroppedVariant));
        Assert.assertFalse(indelsWithoutSoftClips.stream().anyMatch(vc -> vc.getStart() == startOfDroppedVariant && vc.getEnd() == endOfDroppedVariant));

    }


    // make sure that force calling with given alleles that normally wouldn't be called due to complete lack of coverage
    // doesn't run into any edge case bug involving empty likelihoods matrices
    @Test
    public void testForceCallingZeroCoverage() {
        Utils.resetRandomGenerator();
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");
        final File forceCalls = new File(toolsTestDir, "mutect/gga_mode_2.vcf");

        runMutect2(DREAM_3_TUMOR, unfilteredVcf, "20:1119000-1120000", b37Reference, Optional.empty(),
                args -> args.add(AssemblyBasedCallerArgumentCollection.FORCE_CALL_ALLELES_LONG_NAME, forceCalls));
    }

    // make sure we have fixed a bug where germline resources with AF=. throw errors
    @Test
    public void testMissingAF() {
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");
        runMutect2(DREAM_4_TUMOR, unfilteredVcf, "20:1119000-1120000", b37Reference, Optional.of(GNOMAD_WITHOUT_AF_SNIPPET),
                args -> args.addInterval(new SimpleInterval("20:10837425-10837426")));
    }

    // make sure we have fixed a bug where germline resources with less parsimonious representations (eg AC -> TC when M2 calls A->T) still
    // yield a correct annotation
    @Test
    public void testUnparsimoniousGermlineResource() {
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");
        final int position1 = 1366255;
        final int position2 = 1401389;
        final int position3 = 2289449;
        final double defaultAF = 9.78e-7;
        final String interval1 = "20:" + (position1 - 1000) + "-" + (position1 + 1000);
        final String interval2 = "20:" + (position2 - 1000) + "-" + (position2 + 1000);
        final String interval3 = "20:" + (position3 - 1000) + "-" + (position3 + 1000);
        runMutect2(DREAM_4_TUMOR, unfilteredVcf, interval1, b37Reference, Optional.of(UNPARSIMONIOUS_GNOMAD_SNIPPET),
                args -> args.addInterval(interval2).addInterval(interval3).add(M2ArgumentCollection.DEFAULT_AF_LONG_NAME, defaultAF)
                        .addFlag(M2ArgumentCollection.GENOTYPE_GERMLINE_SITES_LONG_NAME));

        // when we call on these intervals we should get variants 20:1366255 G->A, 20:1401389 G->C, and 20:2289449 TCTGGGGACAAA->T
        // In the germline resource file the alleles are represented as GT -> AT (equivalent to G->A), GT -> CA (not equivalent to G->C),
        // and TCTGGGGACAAA->TAA (equivalent to TCTGGGGACA->T) with an extraneous TCTGGGGACAAA->T, respectively
        // In the first case, we should get AF = 0.1 from GG -> AG = G -> A; in the second case we should get the default AF; in the third case
        // we should get AF = 0.1 from TCTGGGGACAAA->TAA and NOT AF = 0.2 from TCTGGGGACAAA->T

        final Map<Integer, double[]> startToAFArray = VariantContextTestUtils.streamVcf(unfilteredVcf).collect(Collectors.toMap(VariantContext::getStart,
                vc -> MathUtils.applyToArray(VariantContextGetters.getAttributeAsDoubleArray(vc, GATKVCFConstants.POPULATION_AF_KEY), x -> Math.pow(10, -x))));

        Assert.assertTrue(startToAFArray.containsKey(position1));
        Assert.assertTrue(startToAFArray.containsKey(position2));
        Assert.assertTrue(startToAFArray.containsKey(position3));

        Assert.assertEquals(startToAFArray.get(position1), new double[] {0.1}, 1e-8);
        Assert.assertEquals(startToAFArray.get(position2), new double[] {defaultAF}, 1e-8);
        Assert.assertEquals(startToAFArray.get(position3), new double[] {0.1}, 1e-8);
    }

    @Test
    public void testContaminationFilter() {
        Utils.resetRandomGenerator();
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");
        runMutect2(new File(NA12878_20_21_WGS_bam), unfilteredVcf, "20:10000000-20010000", b37Reference, Optional.of(GNOMAD));

        final Map<Integer, Set<VariantContext>> filteredVariants = Arrays.stream(new int[] {0, 5, 10}).boxed().collect(Collectors.toMap(pct -> pct, pct -> {
            final File filteredVcf = createTempFile("filtered-" + pct, ".vcf");
            final File contaminationTable = pct == 0 ? NO_CONTAMINATION_TABLE :
                    (pct == 5 ? FIVE_PCT_CONTAMINATION_TABLE : TEN_PCT_CONTAMINATION_TABLE);
            runFilterMutectCalls(unfilteredVcf, filteredVcf, b37Reference,
                    args -> args.add(M2FiltersArgumentCollection.CONTAMINATION_TABLE_LONG_NAME, contaminationTable));

            return VariantContextTestUtils.streamVcf(filteredVcf).collect(Collectors.toSet());
        }));


        final int variantsFilteredAtZeroPercent = (int) filteredVariants.get(0).stream()
                .filter(vc -> vc.getFilters().contains(GATKVCFConstants.CONTAMINATION_FILTER_NAME))
                .count();

        final List<VariantContext> variantsFilteredAtFivePercent = filteredVariants.get(5).stream()
                .filter(vc -> vc.getFilters().contains(GATKVCFConstants.CONTAMINATION_FILTER_NAME)).collect(Collectors.toList());
        Assert.assertEquals(variantsFilteredAtZeroPercent, 0);
        Assert.assertTrue(variantsFilteredAtFivePercent.size() <
                filteredVariants.get(10).stream().filter(vc -> vc.getFilters().contains(GATKVCFConstants.CONTAMINATION_FILTER_NAME)).count());

        final List<VariantContext> missedObviousVariantsAtTenPercent = filteredVariants.get(10).stream()
                .filter(vc -> !vc.isFiltered())
                .filter(VariantContext::isBiallelic)
                .filter(vc -> {
                    final int[] AD = vc.getGenotype(0).getAD();
                    return AD[1] < 0.15 * AD[0];
                }).collect(Collectors.toList());

        Assert.assertTrue(missedObviousVariantsAtTenPercent.isEmpty());

        // If the filter is smart, it won't filter variants with allele fraction much higher than the contamination
        final List<VariantContext> highAlleleFractionFilteredVariantsAtFivePercent = variantsFilteredAtFivePercent.stream()
                .filter(VariantContext::isBiallelic)
                .filter(vc -> vc.getAttributeAsDouble(GATKVCFConstants.CONTAMINATION_QUAL_KEY, 100) < 30)
                .filter(vc -> {
                    final int[] AD = vc.getGenotype(0).getAD();
                    return MathUtils.sum(AD) > 30 && AD[1] > AD[0];
                }).collect(Collectors.toList());

        Assert.assertTrue(highAlleleFractionFilteredVariantsAtFivePercent.isEmpty());
    }

    // test that ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT removes reads that consume zero reference bases
    // e.g. read name HAVCYADXX150109:1:2102:20528:2129 with cigar 23S53I
    @Test
    public void testReadsThatConsumeZeroReferenceReads()  {
        final File bam = new File(publicTestDir + "org/broadinstitute/hellbender/tools/mutect/na12878-chr20-consumes-zero-reference-bases.bam");
        final File outputVcf = createTempFile("output", ".vcf");

        runMutect2(bam, outputVcf, CHROMOSOME_20, b37Reference, Optional.empty());
    }

    // make sure that unpaired reads that pass filtering do not cause errors
    // in particular, the read HAVCYADXX150109:1:1109:11610:46575 with SAM flag 16 fails without the patch
    @Test
    public void testUnpairedReads()  {
        final File bam = new File(toolsTestDir + "unpaired.bam");
        final File outputVcf = createTempFile("output", ".vcf");

        runMutect2(bam, outputVcf, CHROMOSOME_20, b37Reference, Optional.empty());
    }

    // some bams from external pipelines use faulty adapter trimming programs that introduce identical repeated reads
    // into bams.  Although these bams fail the Picard tool ValidateSamFile, one can run HaplotypeCaller and Mutect on them
    // and get fine results.  This test ensures that this remains the case.  The test bam is a small chunk of reads surrounding
    // a germline SNP in NA12878, where we have duplicated 40 of the reads. (In practice bams of this nature might have one bad read
    // per megabase).
    @Test
    public void testBamWithRepeatedReads() {
        final File bam = new File(toolsTestDir + "mutect/repeated_reads.bam");
        final File outputVcf = createTempFile("output", ".vcf");

        runMutect2(bam, outputVcf, "20:10018000-10020000", b37Reference, Optional.empty());
    }

    // In the rare case that a particular fragment is only supported by supplementary reads, that should not
    // result in an exception.  This test ensures that Mutect2 does not fail to finish in that case.
    @Test
    public void testBamWithOnlySupplementaryReads() {
        final File bam = new File(toolsTestDir + "mutect/only_supplementary_reads.bam");
        final File outputVcf = createTempFile("output", ".vcf");

        runMutect2(bam, outputVcf, "20:10018000-10020000", b37Reference, Optional.empty());
    }

    // basic test on a small chunk of NA12878 mitochondria.  This is not a validation, but rather a sanity check
    // that M2 makes obvious calls, doesn't trip up on the beginning of the circular chromosome, and can handle high depth
    @Test
    public void testMitochondria()  {
        Utils.resetRandomGenerator();
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");

        runMutect2(NA12878_MITO_BAM, unfilteredVcf, "chrM:1-1000", MITO_REF.getAbsolutePath(), Optional.empty(),
                args -> args.add(M2ArgumentCollection.MITOCHONDRIA_MODE_LONG_NAME, true));

        final List<VariantContext> variants = VariantContextTestUtils.streamVcf(unfilteredVcf).collect(Collectors.toList());
        final List<String> variantKeys = variants.stream().map(VariantContextTestUtils::keyForVariant).collect(Collectors.toList());

        final List<String> expectedKeys = Arrays.asList(
                "chrM:152-152 T*, [C]",
                "chrM:263-263 A*, [G]",
                "chrM:302-302 A*, [AC, ACC, ACCCCCCCCCCCCC, C]",
                "chrM:310-310 T*, [C, TC]",
                "chrM:750-750 A*, [G]");
        Assert.assertTrue(variantKeys.containsAll(expectedKeys));

        Assert.assertEquals(variants.get(0).getAttributeAsInt(GATKVCFConstants.ORIGINAL_CONTIG_MISMATCH_KEY, 0), 1741);
    }

    // basic test on a small chunk of NA12878 mitochondria. This is intended to be a sanity check to ensure that the pileup-caller
    // is A) Hooked up and functioning in Mutect and B) its arguments can be used to make reasonable calls when assembly doesn't detect variants
    @Test
    public void testPileupDetection()  {
        Utils.resetRandomGenerator();
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");

        runMutect2(NA12878_MITO_BAM, unfilteredVcf, "chrM:1-1000", MITO_REF.getAbsolutePath(), Optional.empty(),
                args -> args.add(M2ArgumentCollection.MITOCHONDRIA_MODE_LONG_NAME, true)
                        // Pileup-Detection arguments
                        .add(PileupDetectionArgumentCollection.PILEUP_DETECTION_LONG_NAME, true)
                        .add(PileupDetectionArgumentCollection.PILEUP_DETECTION_ENABLE_INDELS, true)
                        //.add(PileupDetectionArgumentCollection.PILEUP_DETECTION_ACTIVE_REGION_LOD_THRESHOLD_LONG_NAME, 0) // This is currently incompatible with Mutect 2 (isActive() doesn't store the original activity scores in M2Engine)

                        // NOTE: These arguments are intended to force the assembly to fail at all sites in order to test what gets recovered by the pileup code
                        .add(ReadThreadingAssemblerArgumentCollection.KMER_SIZE_LONG_NAME, 1)
                        .add(ReadThreadingAssemblerArgumentCollection.DONT_INCREASE_KMER_SIZE_LONG_NAME, true)
        );

        final List<VariantContext> variants = VariantContextTestUtils.streamVcf(unfilteredVcf).collect(Collectors.toList());
        final List<String> variantKeys = variants.stream().map(VariantContextTestUtils::keyForVariant).collect(Collectors.toList());
        System.out.print(variantKeys);

        // NOTE: This list is un-principled, its simply asserting that we found SOMETHING with the pileupcaller we should because assembly is broken
        //       this is not set in stone, future changes to the pileupcaller/Mutect2 might cause this list to change
        final List<String> expectedKeys = Arrays.asList(
                "chrM:152-152 T*, [C]",
                "chrM:263-263 A*, [G]",
                "chrM:302-302 A*, [AC]", //Pileup-calling doesn't detect overlapping indels as well as assembly
                "chrM:310-310 T*, [C, TC]",
                "chrM:750-750 A*, [G]");
        Assert.assertTrue(variantKeys.containsAll(expectedKeys));
    }

    /**
     * Several difficult force calling sites including regression test for a thorny force calling bug involving
     * T -> A at chrM:8316.  Previously this call was lost while trimming the assemblyResultSet because every haplotype
     * containing the allele contained a deletion that terminated exactly where trimming occurred, thereby causing the
     * GATK to discard the trimmed haplotype.
     */
    @Test
    public void testDifficultForceCalling() {
        final File dir = new File(largeFileTestDir + "mutect/mito/");
        final File ref = new File(dir, "mito_shifted_8000.fasta");
        final File bam = new File(dir, "mito.bam");
        final File force = new File(dir, "alleles.vcf");
        final String interval = "chrM:8023-9140";

        final File output = createTempFile("output", ".bam");

        runMutect2(bam, output, interval, ref.getAbsolutePath(), Optional.empty(),
                args -> args.add(AssemblyBasedCallerArgumentCollection.FORCE_CALL_ALLELES_LONG_NAME, force),
                args -> args.addFlag(M2ArgumentCollection.MITOCHONDRIA_MODE_LONG_NAME),
                args -> args.add(AssemblyRegionArgumentCollection.MAX_STARTS_LONG_NAME, 75));

        final SimpleInterval callingInterval = new SimpleInterval(interval);

        final Map<Integer, VariantContext> callsByLocus = VariantContextTestUtils.getVariantContexts(output).stream()
                .collect(Collectors.toMap(vc -> vc.getStart(), vc -> vc));

        final List<VariantContext> forceCallingVariants = VariantContextTestUtils.getVariantContexts(force).stream()
                .filter(vc -> vc.overlaps(callingInterval))
                .collect(Collectors.toList());

        for (final VariantContext vc : forceCallingVariants) {
            Assert.assertTrue(callsByLocus.containsKey(vc.getStart()));
            Assert.assertTrue(callsByLocus.get(vc.getStart()).hasAllele(vc.getAlternateAllele(0)));
        }
    }

    @DataProvider(name = "vcfsForFiltering")
    public Object[][] vcfsForFiltering() {
        return new Object[][]{
                {NA12878_MITO_VCF, 0.5, Collections.emptyList(), Arrays.asList(
                        ImmutableSet.of(GATKVCFConstants.STRAND_ARTIFACT_FILTER_NAME, GATKVCFConstants.STRICT_STRAND_BIAS_FILTER_NAME),
                        Collections.emptySet(),
                        ImmutableSet.of( GATKVCFConstants.TUMOR_EVIDENCE_FILTER_NAME,
                                GATKVCFConstants.ALLELE_FRACTION_FILTER_NAME),
                        Collections.emptySet(),
                        Collections.emptySet(),
                        ImmutableSet.of(GATKVCFConstants.DUPLICATED_EVIDENCE_FILTER_NAME),
                        ImmutableSet.of(GATKVCFConstants.FAIL)),
                        Arrays.asList(
                                Arrays.asList(GATKVCFConstants.STRAND_ARTIFACT_FILTER_NAME + ", " + GATKVCFConstants.STRICT_STRAND_BIAS_FILTER_NAME), // strand_bias, strict_stand
                                Arrays.asList(GATKVCFConstants.SITE_LEVEL_FILTERS), // SITE
                                Arrays.asList(GATKVCFConstants.TUMOR_EVIDENCE_FILTER_NAME + ", " + GATKVCFConstants.ALLELE_FRACTION_FILTER_NAME), // weak_evidence, low_allele_frac
                                Arrays.asList(GATKVCFConstants.SITE_LEVEL_FILTERS, GATKVCFConstants.TUMOR_EVIDENCE_FILTER_NAME + ", " + GATKVCFConstants.STRAND_ARTIFACT_FILTER_NAME + ", " + GATKVCFConstants.ALLELE_FRACTION_FILTER_NAME, GATKVCFConstants.STRAND_ARTIFACT_FILTER_NAME + ", " + GATKVCFConstants.STRICT_STRAND_BIAS_FILTER_NAME + ", " + GATKVCFConstants.ALLELE_FRACTION_FILTER_NAME), // SITE|weak_evidence, strand_bias, low_allele_frac|strand_bias, strict_strand, low_allele_frac
                                Arrays.asList(GATKVCFConstants.SITE_LEVEL_FILTERS),  // SITE
                                Arrays.asList(GATKVCFConstants.DUPLICATED_EVIDENCE_FILTER_NAME), // duplicate
                                Arrays.asList(GATKVCFConstants.TUMOR_EVIDENCE_FILTER_NAME + ", " + GATKVCFConstants.STRAND_ARTIFACT_FILTER_NAME + ", " + GATKVCFConstants.STRICT_STRAND_BIAS_FILTER_NAME, GATKVCFConstants.ALLELE_FRACTION_FILTER_NAME) // weak_evidence, strand_bias, strict_stand|low_allele_frac

                        )},
                {NA12878_MITO_GVCF, .0009, Arrays.asList("MT:1", "MT:37", "MT:40", "MT:152", "MT:157"), Arrays.asList(
                        Collections.emptySet(),
                        ImmutableSet.of(GATKVCFConstants.MEDIAN_BASE_QUALITY_FILTER_NAME, GATKVCFConstants.TUMOR_EVIDENCE_FILTER_NAME, GATKVCFConstants.STRAND_ARTIFACT_FILTER_NAME),
                        ImmutableSet.of(GATKVCFConstants.TUMOR_EVIDENCE_FILTER_NAME,GATKVCFConstants.STRAND_ARTIFACT_FILTER_NAME, GATKVCFConstants.STRICT_STRAND_BIAS_FILTER_NAME),
                        Collections.emptySet(),
                        ImmutableSet.of(GATKVCFConstants.MEDIAN_BASE_QUALITY_FILTER_NAME, GATKVCFConstants.CONTAMINATION_FILTER_NAME, GATKVCFConstants.ALLELE_FRACTION_FILTER_NAME,
                                GATKVCFConstants.TUMOR_EVIDENCE_FILTER_NAME, GATKVCFConstants.READ_POSITION_FILTER_NAME, GATKVCFConstants.MEDIAN_MAPPING_QUALITY_FILTER_NAME, GATKVCFConstants.STRICT_STRAND_BIAS_FILTER_NAME)),
                        Arrays.asList(
                                Arrays.asList(GATKVCFConstants.SITE_LEVEL_FILTERS), // SITE,
                                Arrays.asList(GATKVCFConstants.TUMOR_EVIDENCE_FILTER_NAME + ", " + GATKVCFConstants.MEDIAN_BASE_QUALITY_FILTER_NAME + ", " + GATKVCFConstants.STRAND_ARTIFACT_FILTER_NAME, GATKVCFConstants.SITE_LEVEL_FILTERS), //"weak_evidence, base_qual, strand_bias|SITE",
                                Arrays.asList(GATKVCFConstants.TUMOR_EVIDENCE_FILTER_NAME + ", " + GATKVCFConstants.STRAND_ARTIFACT_FILTER_NAME + ", " + GATKVCFConstants.STRICT_STRAND_BIAS_FILTER_NAME, GATKVCFConstants.SITE_LEVEL_FILTERS), // "weak_evidence, strict_strand, strand_bias|SITE",
                                Arrays.asList(GATKVCFConstants.SITE_LEVEL_FILTERS, GATKVCFConstants.TUMOR_EVIDENCE_FILTER_NAME + ", " + GATKVCFConstants.MEDIAN_BASE_QUALITY_FILTER_NAME + ", " + GATKVCFConstants.STRAND_ARTIFACT_FILTER_NAME + ", " + GATKVCFConstants.ALLELE_FRACTION_FILTER_NAME, GATKVCFConstants.SITE_LEVEL_FILTERS), //".|weak_evidence, base_qual, strand_bias, low_allele_frac|SITE",
                                Arrays.asList(GATKVCFConstants.TUMOR_EVIDENCE_FILTER_NAME + ", " + GATKVCFConstants.MEDIAN_BASE_QUALITY_FILTER_NAME + ", " + GATKVCFConstants.MEDIAN_MAPPING_QUALITY_FILTER_NAME + ", " + GATKVCFConstants.CONTAMINATION_FILTER_NAME + ", " + GATKVCFConstants.STRICT_STRAND_BIAS_FILTER_NAME + ", " + GATKVCFConstants.READ_POSITION_FILTER_NAME + ", " + GATKVCFConstants.ALLELE_FRACTION_FILTER_NAME, GATKVCFConstants.SITE_LEVEL_FILTERS) // "weak_evidence, base_qual, map_qual, contamination, strand_artifact, position, low_allele_frac|SITE"
                        )},
                {NA12878_MITO_GVCF, .0009, Arrays.asList("MT:15326"), Arrays.asList(Collections.emptySet()), Arrays.asList(Arrays.asList(GATKVCFConstants.SITE_LEVEL_FILTERS, GATKVCFConstants.SITE_LEVEL_FILTERS))}
        };
    }

    @Test(dataProvider = "vcfsForFiltering")
    public void testFilterMitochondria(File unfiltered, final double minAlleleFraction, final List<String> intervals, List<Set<String>> expectedFilters, List<List<String>> expectedASFilters)  {
        final File filteredVcf = createTempFile("filtered", ".vcf");

        // vcf sequence dicts don't match ref
        runFilterMutectCalls(unfiltered, filteredVcf, MITO_REF.getAbsolutePath(),
                args -> args.add(M2ArgumentCollection.MITOCHONDRIA_MODE_LONG_NAME, true),
                args -> args.add(StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME, true),
                args -> args.add(M2FiltersArgumentCollection.MIN_AF_LONG_NAME, minAlleleFraction),
                args -> args.add(M2FiltersArgumentCollection.MIN_READS_ON_EACH_STRAND_LONG_NAME, 1),
                args -> args.add(M2FiltersArgumentCollection.UNIQUE_ALT_READ_COUNT_LONG_NAME, 2),
                args -> {
                    intervals.stream().map(SimpleInterval::new).forEach(args::addInterval);
                    return args;
                });

        final List<Set<String>> actualFilters = VariantContextTestUtils.streamVcf(filteredVcf)
                .map(VariantContext::getFilters).collect(Collectors.toList());

        final List<List<String>> actualASFilters = VariantContextTestUtils.streamVcf(filteredVcf)
                .map(vc -> AnnotationUtils.decodeAnyASListWithRawDelim(vc.getCommonInfo().getAttributeAsString(GATKVCFConstants.AS_FILTER_STATUS_KEY, ""))).collect(Collectors.toList());
        Assert.assertEquals(actualASFilters, expectedASFilters);

        Assert.assertEquals(actualFilters.size(), expectedFilters.size());
        for (int n = 0; n < actualFilters.size(); n++) {
            Assert.assertTrue(actualFilters.get(n).containsAll(expectedFilters.get(n)), "Actual filters missing some expected filters: " + SetUtils.difference(expectedFilters.get(n), actualFilters.get(n)));
            Assert.assertTrue(expectedFilters.get(n).containsAll(actualFilters.get(n)), "Expected filters missing some actual filters: " + SetUtils.difference(actualFilters.get(n), expectedFilters.get(n)));
        }

        Assert.assertEquals(actualFilters, expectedFilters);
    }

    @Test(dataProvider = "vcfsForFiltering")
    public void testFilterMicrobial(File unfiltered, final double minAlleleFraction, final List<String> intervals, List<Set<String>> expectedFilters, List<List<String>> expectedASFilters)  {
        final File filteredVcf = createTempFile("filtered", ".vcf");

        // vcf sequence dicts don't match ref
        runFilterMutectCalls(unfiltered, filteredVcf, MITO_REF.getAbsolutePath(),
                args -> args.add(M2ArgumentCollection.MICROBIAL_MODE_LONG_NAME, true),
                args -> args.add(StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME, true),
                args -> args.add(M2FiltersArgumentCollection.MIN_AF_LONG_NAME, minAlleleFraction),
                args -> args.add(M2FiltersArgumentCollection.MIN_READS_ON_EACH_STRAND_LONG_NAME, 1),
                args -> args.add(M2FiltersArgumentCollection.UNIQUE_ALT_READ_COUNT_LONG_NAME, 2),
                args -> {
                    intervals.stream().map(SimpleInterval::new).forEach(args::addInterval);
                    return args;
                });

        final List<Set<String>> actualFilters = VariantContextTestUtils.streamVcf(filteredVcf)
                .map(VariantContext::getFilters).collect(Collectors.toList());

        final List<List<String>> actualASFilters = VariantContextTestUtils.streamVcf(filteredVcf)
                .map(vc -> AnnotationUtils.decodeAnyASListWithRawDelim(vc.getCommonInfo().getAttributeAsString(GATKVCFConstants.AS_FILTER_STATUS_KEY, ""))).collect(Collectors.toList());
        Assert.assertEquals(actualASFilters, expectedASFilters);

        Assert.assertEquals(actualFilters.size(), expectedFilters.size());
        for (int n = 0; n < actualFilters.size(); n++) {
            Assert.assertTrue(actualFilters.get(n).containsAll(expectedFilters.get(n)), "Actual filters missing some expected filters: " + SetUtils.difference(expectedFilters.get(n), actualFilters.get(n)));
            Assert.assertTrue(expectedFilters.get(n).containsAll(actualFilters.get(n)), "Expected filters missing some actual filters: " + SetUtils.difference(actualFilters.get(n), expectedFilters.get(n)));
        }

        Assert.assertEquals(actualFilters, expectedFilters);
    }

    @DataProvider(name = "vcfsForNuMTFiltering")
    public Object[][] vcfsForNuMTFiltering() {
        return new Object[][]{
                {NA12878_MITO_INITIAL_FILTERED_VCF, 30, Collections.emptyList(), Arrays.asList(
                        ImmutableSet.of(GATKVCFConstants.STRAND_ARTIFACT_FILTER_NAME, GATKVCFConstants.STRICT_STRAND_BIAS_FILTER_NAME),
                        Collections.emptySet(),
                        ImmutableSet.of( GATKVCFConstants.TUMOR_EVIDENCE_FILTER_NAME,
                                GATKVCFConstants.POSSIBLE_NUMT_FILTER_NAME,
                                GATKVCFConstants.ALLELE_FRACTION_FILTER_NAME),
                        Collections.emptySet(),
                        Collections.emptySet(),
                        ImmutableSet.of(GATKVCFConstants.DUPLICATED_EVIDENCE_FILTER_NAME),
                        ImmutableSet.of(GATKVCFConstants.FAIL)),
                        Arrays.asList(
                                Arrays.asList(GATKVCFConstants.STRAND_ARTIFACT_FILTER_NAME + ", " + GATKVCFConstants.STRICT_STRAND_BIAS_FILTER_NAME), // strand_bias, strict_stand
                                Arrays.asList(GATKVCFConstants.SITE_LEVEL_FILTERS), // SITE
                                Arrays.asList(GATKVCFConstants.TUMOR_EVIDENCE_FILTER_NAME + ", " + GATKVCFConstants.ALLELE_FRACTION_FILTER_NAME + ", "  + GATKVCFConstants.POSSIBLE_NUMT_FILTER_NAME), // weak_evidence, low_allele_frac, possible_numt
                                Arrays.asList(GATKVCFConstants.SITE_LEVEL_FILTERS, GATKVCFConstants.TUMOR_EVIDENCE_FILTER_NAME + ", " + GATKVCFConstants.STRAND_ARTIFACT_FILTER_NAME + ", " + GATKVCFConstants.ALLELE_FRACTION_FILTER_NAME + ", " + GATKVCFConstants.POSSIBLE_NUMT_FILTER_NAME, GATKVCFConstants.STRAND_ARTIFACT_FILTER_NAME + ", " + GATKVCFConstants.STRICT_STRAND_BIAS_FILTER_NAME + ", " + GATKVCFConstants.ALLELE_FRACTION_FILTER_NAME + ", " + GATKVCFConstants.POSSIBLE_NUMT_FILTER_NAME), // SITE|weak_evidence, strand_bias, low_allele_frac, possible_numt|strand_bias, strict_strand, low_allele_frac, possible_numt
                                Arrays.asList(GATKVCFConstants.SITE_LEVEL_FILTERS),  // SITE
                                Arrays.asList(GATKVCFConstants.DUPLICATED_EVIDENCE_FILTER_NAME), // duplicate
                                Arrays.asList(GATKVCFConstants.TUMOR_EVIDENCE_FILTER_NAME + ", " + GATKVCFConstants.STRAND_ARTIFACT_FILTER_NAME + ", " + GATKVCFConstants.STRICT_STRAND_BIAS_FILTER_NAME, GATKVCFConstants.ALLELE_FRACTION_FILTER_NAME + ", " + GATKVCFConstants.POSSIBLE_NUMT_FILTER_NAME) // weak_evidence, strand_bias, strict_stand|low_allele_frac, possible_numt

                        )}
        };
    }

    @Test(dataProvider = "vcfsForNuMTFiltering")
    public void testNuMTFilterMitochondria(File initialFilters, final double autosomalCoverage, final List<String> intervals, List<Set<String>> expectedFilters, List<List<String>> expectedASFilters)  {
        final File filteredVcf = createTempFile("filtered", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addVCF(initialFilters)
                .addOutput(filteredVcf)
                .addReference(MITO_REF.getAbsolutePath())
                .add(NuMTFilterTool.MEDIAN_AUTOSOMAL_COVERAGE_LONG_NAME, autosomalCoverage)
                .add(NuMTFilterTool.MAX_NUMT_COPIES_IN_AUTOSOME_LONG_NAME, 4.0);

        intervals.stream().map(SimpleInterval::new).forEach(args::addInterval);

        // vcf sequence dicts don't match ref
        runCommandLine(args, NuMTFilterTool.class.getSimpleName());

        final List<Set<String>> actualFilters = VariantContextTestUtils.streamVcf(filteredVcf)
                .map(VariantContext::getFilters).collect(Collectors.toList());

        final List<List<String>> actualASFilters = VariantContextTestUtils.streamVcf(filteredVcf)
                .map(vc -> AnnotationUtils.decodeAnyASListWithRawDelim(vc.getCommonInfo().getAttributeAsString(GATKVCFConstants.AS_FILTER_STATUS_KEY, ""))).collect(Collectors.toList());
        Assert.assertEquals(actualASFilters, expectedASFilters);

        Assert.assertEquals(actualFilters.size(), expectedFilters.size());
        for (int n = 0; n < actualFilters.size(); n++) {
            Assert.assertTrue(actualFilters.get(n).containsAll(expectedFilters.get(n)), "Actual filters missing some expected filters: " + SetUtils.difference(expectedFilters.get(n), actualFilters.get(n)));
            Assert.assertTrue(expectedFilters.get(n).containsAll(actualFilters.get(n)), "Expected filters missing some actual filters: " + SetUtils.difference(actualFilters.get(n), expectedFilters.get(n)));
        }

        Assert.assertEquals(actualFilters, expectedFilters);
    }

    @Test
    public void testMitochondrialRefConf()  {
        Utils.resetRandomGenerator();
        final File standardVcf = createTempFile("standard", ".vcf");
        final File unthresholded = createTempFile("unthresholded", ".vcf");
        final double minAF = 0.01;

        runMutect2(NA12878_MITO_BAM, standardVcf, "chrM:1-1000", MITO_REF.getAbsolutePath(), Optional.empty(),
                args -> args.add(AssemblyBasedCallerArgumentCollection.EMIT_REF_CONFIDENCE_LONG_NAME, ReferenceConfidenceMode.GVCF.toString()),
                args -> args.add(M2ArgumentCollection.MINIMUM_ALLELE_FRACTION_LONG_NAME, minAF),
                args -> args.add(M2ArgumentCollection.LOD_BAND_LONG_NAME, -2.0),
                args -> args.add(M2ArgumentCollection.LOD_BAND_LONG_NAME, 0.0));

        //check ref conf-specific headers are output
        final Pair<VCFHeader, List<VariantContext>> result = VariantContextTestUtils.readEntireVCFIntoMemory(standardVcf.getAbsolutePath());
        Assert.assertTrue(result.getLeft().hasFormatLine(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY));
        Assert.assertTrue(result.getLeft().getMetaDataLine(GATKVCFConstants.SYMBOLIC_ALLELE_DEFINITION_HEADER_TAG) != null);

        final List<VariantContext> variants = result.getRight();
        final Map<String, VariantContext> variantMap = variants.stream().collect(Collectors.toMap(VariantContextTestUtils::keyForVariant, Function.identity()));
        final List<String> variantKeys = new ArrayList<>(variantMap.keySet());

        final List<String> expectedKeys = Arrays.asList(
                "chrM:152-152 T*, [<NON_REF>, C]",
                "chrM:263-263 A*, [<NON_REF>, G]",
                //"chrM:297-297 A*, [<NON_REF>, AC, C]",
                //"chrM:301-301 A*, [<NON_REF>, AC, ACC, ACCC]",
                "chrM:302-302 A*, [<NON_REF>, AC, ACC, ACCC, C]",  //one of these commented out variants has an allele that only appears in debug mode
                "chrM:310-310 T*, [<NON_REF>, TC]",
                "chrM:750-750 A*, [<NON_REF>, G]");
        Assert.assertTrue(variantKeys.containsAll(expectedKeys));
        //First entry should be a homRef block
        Assert.assertTrue(VariantContextTestUtils.keyForVariant(variants.get(0)).contains("*, [<NON_REF>]"));

        final ArgumentsBuilder validateVariantsArgs = new ArgumentsBuilder()
                .add("R", MITO_REF.getAbsolutePath())
                .add("V", standardVcf.getAbsolutePath())
                .add("L", IntervalUtils.locatableToString(new SimpleInterval("chrM:1-1000")))
                .addRaw("-gvcf");
        runCommandLine(validateVariantsArgs, ValidateVariants.class.getSimpleName());

        runMutect2(NA12878_MITO_BAM, unthresholded, "chrM:1-1000", MITO_REF.getAbsolutePath(), Optional.empty(),
                args -> args.add(AssemblyBasedCallerArgumentCollection.EMIT_REF_CONFIDENCE_LONG_NAME, ReferenceConfidenceMode.GVCF.toString()),
                args -> args.add(M2ArgumentCollection.MINIMUM_ALLELE_FRACTION_LONG_NAME, 0.00),
                args -> args.add(M2ArgumentCollection.LOD_BAND_LONG_NAME, -2.0),
                args -> args.add(M2ArgumentCollection.LOD_BAND_LONG_NAME, 0.0));

        final Pair<VCFHeader, List<VariantContext>> result_noThreshold = VariantContextTestUtils.readEntireVCFIntoMemory(unthresholded.getAbsolutePath());

        final Map<String, VariantContext> variantMap2 = result_noThreshold.getRight().stream().collect(Collectors.toMap(VariantContextTestUtils::keyForVariant, Function.identity()));

        //TLODs for variants should not change too much for variant allele, should change significantly for non-ref
        // however, there are edge cases where this need not be true (this might indicate the need to fix our
        // LOD calculation for the NON-REF allele), so we allow one anomalous site
        final long changedRegularAlleleLodCount = expectedKeys.stream()
                .filter(key -> !onlyNonRefTlodsChange(variantMap.get(key), variantMap2.get(key), minAF))
                .count();

        Assert.assertTrue(changedRegularAlleleLodCount <= 1);

        final List<String> expectedRefKeys = Arrays.asList(
                //ref blocks will be dependent on TLOD band values
                "chrM:218-218 A*, [<NON_REF>]",
                "chrM:264-266 C*, [<NON_REF>]",
                "chrM:475-483 A*, [<NON_REF>]",
                "chrM:488-492 T*, [<NON_REF>]");

        //ref block boundaries aren't particularly stable, so try a few and make sure we check at least one
        final List<String> refBlockKeys = expectedRefKeys.stream()
                .filter(key -> variantMap.containsKey(key) && variantMap2.containsKey(key))
                .collect(Collectors.toList());
        Assert.assertFalse(refBlockKeys.isEmpty());

        refBlockKeys.forEach(key -> Assert.assertTrue(onlyNonRefTlodsChange(variantMap.get(key), variantMap2.get(key), minAF)));
    }

    private boolean onlyNonRefTlodsChange(final VariantContext v1, final VariantContext v2, final double minAF) {
        if (v1 == null || v2 == null || !v1.getReference().equals(v2.getReference()) ||
                !(v1.getAlternateAlleles().size() == v2.getAlternateAlleles().size())) {
            return false;
        }

        //ref blocks have TLOD in format field
        final boolean isHomRef = v1.getGenotype(0).isHomRef();
        final double[] tlods1 = !isHomRef ? VariantContextGetters.getAttributeAsDoubleArray(v1, GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY)
                : new double[]{VariantContextGetters.getAttributeAsDouble(v1.getGenotype(0), GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, 0)};
        final double[] tlods2 = !isHomRef ? VariantContextGetters.getAttributeAsDoubleArray(v2, GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY)
                : new double[]{VariantContextGetters.getAttributeAsDouble(v2.getGenotype(0), GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, 0)};

        final double[] af1 = VariantContextGetters.getAttributeAsDoubleArray(v1, GATKVCFConstants.ALLELE_FRACTION_KEY, () -> null, 0);
        final double[] af2 = VariantContextGetters.getAttributeAsDoubleArray(v2, GATKVCFConstants.ALLELE_FRACTION_KEY, () -> null, 0);

        for (int i = 0; i < v1.getAlternateAlleles().size(); i++) {
            if (!v1.getAlternateAllele(i).equals(v2.getAlternateAllele(i))) {
                return false;
            }
            //we expect the AF threshold to have a significant effect on the NON_REF TLOD, but only for that allele
            if (!v1.getAlternateAllele(i).equals(Allele.NON_REF_ALLELE)) {
                if (tlods1[i] > 0) {
                    if (Math.abs(tlods1[i] - tlods2[i]) / tlods1[i] > VARIANT_TLOD_MATCH_PCT) {
                        return false;
                    }
                } else if (af2 != null && af2[i] > minAF && Math.abs(tlods1[i] - tlods2[i]) > TLOD_MATCH_EPSILON) {
                    return false;
                }
            } else {
                if (Math.abs(tlods1[i] - tlods2[i]) < TLOD_MATCH_EPSILON) {
                    return false;
                }
            }
        }
        return true;
    }

    @Test
    @SuppressWarnings("deprecation")
    public void testAFAtHighDP()  {
        Utils.resetRandomGenerator();
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");

        runMutect2(DEEP_MITO_BAM, unfilteredVcf, "chrM:1-1018", MITO_REF.getAbsolutePath(), Optional.empty(),
                args -> args.add(IntervalArgumentCollection.INTERVAL_PADDING_LONG_NAME, 300));

        final List<VariantContext> variants = VariantContextTestUtils.streamVcf(unfilteredVcf).collect(Collectors.toList());

        for (final VariantContext vc : variants) {
            Assert.assertTrue(vc.isBiallelic()); //I do some lazy parsing below that won't hold for multiple alternate alleles
            final Genotype g = vc.getGenotype(DEEP_MITO_SAMPLE_NAME);
            Assert.assertTrue(g.hasAD());
            final int[] ADs = g.getAD();
            Assert.assertTrue(g.hasExtendedAttribute(GATKVCFConstants.ALLELE_FRACTION_KEY));
            //Assert.assertEquals(Double.parseDouble(String.valueOf(vc.getGenotype(DEEP_MITO_SAMPLE_NAME).getExtendedAttribute(GATKVCFConstants.ALLELE_FRACTION_KEY,"0"))), (double)ADs[1]/(ADs[0]+ADs[1]), 1e-6);
            Assert.assertEquals(Double.parseDouble(String.valueOf(vc.getGenotype(DEEP_MITO_SAMPLE_NAME).getAttributeAsString(GATKVCFConstants.ALLELE_FRACTION_KEY, "0"))), (double) ADs[1] / (ADs[0] + ADs[1]), 2e-3);
        }
    }

    // check that the somatic clustering model works with high-depth, low-AF cfDNA clustering
    @Test
    public void testBloodBiopsyFiltering() {
        final File unfiltered = new File(toolsTestDir, "mutect/cfdna/cfdna-unfiltered.vcf");
        final File filtered = createTempFile("filtered", ".vcf");

        runFilterMutectCalls(unfiltered, filtered, b37Reference);

        final Map<Integer, Set<String>> filtersBySite = VariantContextTestUtils.streamVcf(filtered).collect(Collectors.toMap(VariantContext::getStart, VariantContext::getFilters));

        // these are sites that caused trouble for a previous version of the somatic clustering model
        final List<Integer> lowAFSitesWeShouldNotCall = Arrays.asList(25963056, 47162531, 142178205, 151841902, 31325209, 41521982);
        for (final int site : lowAFSitesWeShouldNotCall) {
            Assert.assertTrue(filtersBySite.get(site).contains(GATKVCFConstants.TUMOR_EVIDENCE_FILTER_NAME));
        }
    }

    @Test
    public void testBamout() {
        final File outputVcf = createTempFile("output", ".vcf");
        final File bamout = createTempFile("bamout", ".bam");

        runMutect2(DREAM_1_TUMOR, outputVcf, "20:10000000-13000000", b37Reference, Optional.empty(),
                args -> args.add(AssemblyBasedCallerArgumentCollection.BAM_OUTPUT_LONG_NAME, bamout));
        Assert.assertTrue(bamout.exists());
    }

    // this test is here for basic sanity of flow mode
    @Test
    public void testBamoutFlowMode() {
        final File outputVcf = createTempFile("output", ".vcf");
        final File bamout = createTempFile("bamout", ".bam");

        runMutect2(FLOW_BASED_INPUT_BAM, outputVcf, "chr9:81149486-81177047", hg38Reference, Optional.empty(),
                args -> {
                    args.add(M2ArgumentCollection.FLOW_M2_MODE_LONG_NAME, "ADVANCED");
                    return args.add(AssemblyBasedCallerArgumentCollection.BAM_OUTPUT_LONG_NAME, bamout);

                });
        Assert.assertTrue(bamout.exists());
    }



    @Test
    public void testFilteringHeaders() {
        Utils.resetRandomGenerator();
        final File tumor = new File(NA12878_20_21_WGS_bam);
        final File unfilteredVcf = createTempFile("unfiltered", ".vcf");
        final File filteredVcf = createTempFile("filtered", ".vcf");

        runMutect2(tumor, unfilteredVcf, "20:10000000-10010000", b37Reference, Optional.of(GNOMAD));
        runFilterMutectCalls(unfilteredVcf, filteredVcf, b37Reference);

        final VCFHeader filteredHeader = VariantContextTestUtils.getVCFHeader(filteredVcf.getAbsolutePath());
        //explicit check for PASS header
        Assert.assertTrue(filteredHeader.hasFilterLine(VCFConstants.PASSES_FILTERS_v4));
        for (final String header : GATKVCFConstants.MUTECT_FILTER_NAMES){
            Assert.assertTrue(filteredHeader.hasFilterLine(header));
        }
    }


    @DataProvider(name="ExposureOfSmithWatermanParametersTestInputs")
    public Object[][] getExposureOfSmithWatermanParametersTestInputs() {
        return new Object[][] {
                {NA12878_20_21_WGS_bam, b37_reference_20_21}
        };
    }

    /**
     * See {@link org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerIntegrationTest#testExposureOfSmithWatermanParametersIsConsistentWithPastResults}.
     */
    @Test(dataProvider="ExposureOfSmithWatermanParametersTestInputs")
    public void testExposureOfSmithWatermanParametersIsConsistentWithPastResults(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testExposureOfSmithWatermanParametersIsConsistentWithPastResults", ".vcf");
        final File expected = new File(toolsTestDir + "mutect", "expected.testExposureOfSmithWatermanParameters.M2.gatk4.vcf");

        final String outputPath = output.getAbsolutePath();

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", outputPath,
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--" + AssemblyBasedCallerArgumentCollection.SMITH_WATERMAN_DANGLING_END_MATCH_VALUE_LONG_NAME, "1",
                "--" + AssemblyBasedCallerArgumentCollection.SMITH_WATERMAN_DANGLING_END_MISMATCH_PENALTY_LONG_NAME, "-2",
                "--" + AssemblyBasedCallerArgumentCollection.SMITH_WATERMAN_DANGLING_END_GAP_OPEN_PENALTY_LONG_NAME, "-3",
                "--" + AssemblyBasedCallerArgumentCollection.SMITH_WATERMAN_DANGLING_END_GAP_EXTEND_PENALTY_LONG_NAME, "-4",
                "--" + AssemblyBasedCallerArgumentCollection.SMITH_WATERMAN_HAPLOTYPE_TO_REFERENCE_MATCH_VALUE_LONG_NAME, "5",
                "--" + AssemblyBasedCallerArgumentCollection.SMITH_WATERMAN_HAPLOTYPE_TO_REFERENCE_MISMATCH_PENALTY_LONG_NAME, "-6",
                "--" + AssemblyBasedCallerArgumentCollection.SMITH_WATERMAN_HAPLOTYPE_TO_REFERENCE_GAP_OPEN_PENALTY_LONG_NAME, "-7",
                "--" + AssemblyBasedCallerArgumentCollection.SMITH_WATERMAN_HAPLOTYPE_TO_REFERENCE_GAP_EXTEND_PENALTY_LONG_NAME, "-8",
                "--" + AssemblyBasedCallerArgumentCollection.SMITH_WATERMAN_READ_TO_HAPLOTYPE_MATCH_VALUE_LONG_NAME, "9",
                "--" + AssemblyBasedCallerArgumentCollection.SMITH_WATERMAN_READ_TO_HAPLOTYPE_MISMATCH_PENALTY_LONG_NAME, "-10",
                "--" + AssemblyBasedCallerArgumentCollection.SMITH_WATERMAN_READ_TO_HAPLOTYPE_GAP_OPEN_PENALTY_LONG_NAME, "-11",
                "--" + AssemblyBasedCallerArgumentCollection.SMITH_WATERMAN_READ_TO_HAPLOTYPE_GAP_EXTEND_PENALTY_LONG_NAME, "-12",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
        };

        runCommandLine(args);

        // This used to be an exact text match, which cause hours of aggravation when the test passed locally and
        // failed on the cloud.
        final double[] outputTlods = VariantContextTestUtils.streamVcf(output)
                        .flatMap(vc -> vc.getAttributeAsDoubleList(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, 0).stream())
                .mapToDouble(x -> x).toArray();

        final double[] expectedTlods = VariantContextTestUtils.streamVcf(expected)
                .flatMap(vc -> vc.getAttributeAsDoubleList(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, 0).stream())
                .mapToDouble(x -> x).toArray();

        Assert.assertEquals(outputTlods, expectedTlods);
    }

    @SafeVarargs
    final private void runMutect2(final List<File> tumors, final List<File> normals, final File output,
                                  final String interval, final String reference,
                                  final Optional<File> gnomad, final Function<ArgumentsBuilder, ArgumentsBuilder>... appendExtraArguments) {
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(output)
                .addReference(reference);

        tumors.forEach(args::addInput);

        normals.forEach(normal -> {
            args.addInput(normal);
            args.add(M2ArgumentCollection.NORMAL_SAMPLE_LONG_NAME, getSampleName(normal));
        });

        gnomad.ifPresent(g -> args.add(M2ArgumentCollection.GERMLINE_RESOURCE_LONG_NAME, g));

        args.addInterval(new SimpleInterval(interval));

        ArgumentsBuilder argsWithAdditions = args;

        for (final Function<ArgumentsBuilder, ArgumentsBuilder> extraArgument : appendExtraArguments) {
            argsWithAdditions = extraArgument.apply(args);
        }

        runCommandLine(argsWithAdditions);
    }

    @SafeVarargs
    final private void runMutect2(final File tumor, final File normal, final File output, final String interval, final String reference,
                                  final Optional<File> gnomad, final Function<ArgumentsBuilder, ArgumentsBuilder>... appendExtraArguments) {
        runMutect2(Collections.singletonList(tumor), Collections.singletonList(normal), output, interval, reference, gnomad, appendExtraArguments);
    }

    @SafeVarargs
    final private void runMutect2(final File tumor, final File output, final String interval, final String reference,
                                  final Optional<File> gnomad, final Function<ArgumentsBuilder, ArgumentsBuilder>... appendExtraArguments) {
        runMutect2(Collections.singletonList(tumor), Collections.emptyList(), output, interval, reference, gnomad, appendExtraArguments);
    }

    @SafeVarargs
    final private void runFilterMutectCalls(final File unfilteredVcf, final File filteredVcf, final String reference,
                                            final Function<ArgumentsBuilder, ArgumentsBuilder>... appendExtraArguments) {
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addVCF(unfilteredVcf)
                .addOutput(filteredVcf)
                .addReference(reference);

        ArgumentsBuilder argsWithAdditions = args;

        for (final Function<ArgumentsBuilder, ArgumentsBuilder> extraArgument : appendExtraArguments) {
            argsWithAdditions = extraArgument.apply(args);
        }

        runCommandLine(argsWithAdditions, FilterMutectCalls.class.getSimpleName());
    }

    private void runConcordance(final File truth, final File eval, final File summary, final String interval, final File mask, final Optional<File> truePositivesFalseNegatives) {
        final ArgumentsBuilder concordanceArgs = new ArgumentsBuilder()
                .add(Concordance.TRUTH_VARIANTS_LONG_NAME, truth)
                .add(Concordance.EVAL_VARIANTS_LONG_NAME, eval)
                .addInterval(new SimpleInterval(interval))
                .add(IntervalArgumentCollection.EXCLUDE_INTERVALS_LONG_NAME, mask)
                .add(Concordance.SUMMARY_LONG_NAME, summary);

        truePositivesFalseNegatives.ifPresent(file -> concordanceArgs.add(Concordance.TRUE_POSITIVES_AND_FALSE_NEGATIVES_SHORT_NAME, file));
        runCommandLine(concordanceArgs, Concordance.class.getSimpleName());
    }

    private String getSampleName(final File bam)  {
        try {
            final File nameFile = createTempFile("sample_name", ".txt");
            new Main().instanceMain(makeCommandLineArgs(Arrays.asList("-I", bam.getAbsolutePath(), "-O", nameFile.getAbsolutePath(), "-encode"), "GetSampleName"));
            return Files.readAllLines(nameFile.toPath()).get(0);
        } catch (final IOException ex) {
            throw new IllegalArgumentException(ex);
        }
    }
}
