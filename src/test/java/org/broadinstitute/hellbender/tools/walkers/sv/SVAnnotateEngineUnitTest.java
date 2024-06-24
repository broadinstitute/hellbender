package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.bed.FullBEDFeature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.commons.compress.utils.Sets;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SVInterval;
import org.broadinstitute.hellbender.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.gtf.GencodeGtfGeneFeature;
import org.broadinstitute.hellbender.utils.codecs.gtf.GencodeGtfTranscriptFeature;
import org.testng.Assert;
import org.testng.TestException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

public class SVAnnotateEngineUnitTest extends GATKBaseTest {
    private final File TOY_GTF_FILE = new File(getToolTestDataDir().replaceFirst("Engine", "") + "unittest.gtf");
    private final File TINY_NONCODING_BED_FILE = new File(getToolTestDataDir().replaceFirst("Engine", "") + "noncoding.unittest.bed");


    // Pairs of intervals with different relationships to check if first (variant) interval spans second (feature)
    // interval (boolean) and count how many of the variant breakpoints (start and end of first interval) are within
    // the feature (second) interval (integer)
    @DataProvider(name="variantFeatureComparisons")
    public Object[][] getVariantFeatureComparisonsData() {
        return new Object[][] {
                { new SimpleInterval("chr1", 1, 2000), new SimpleInterval("chr2", 100, 200), false, 0 },
                { new SimpleInterval("chrX", 6000, 8000), new SimpleInterval("chrX", 30, 400), false, 0 },
                { new SimpleInterval("chr3", 10000, 30000), new SimpleInterval("chr3", 9000, 20000), false, 1 },
                { new SimpleInterval("chrY", 500, 700), new SimpleInterval("chrY", 600, 800), false, 1 },
                { new SimpleInterval("chr22", 40, 90), new SimpleInterval("chr22", 25, 700), false, 2 },
                { new SimpleInterval("chr19", 33, 33), new SimpleInterval("chr19", 22, 44), false, 2 },
                { new SimpleInterval("chr10", 900, 4000), new SimpleInterval("chr10", 2000, 2200), true, 0 },
        };
    }

    // Check variant / feature interval comparison functions for span and number of breakpoints
    @Test(dataProvider = "variantFeatureComparisons")
    public void testVariantFeatureComparisons(
            final SimpleInterval variantInterval,
            final SimpleInterval featureInterval,
            final boolean expectedVariantSpansFeature,
            final int expectedNumBreakpointsInsideFeature)
    {
        final boolean actualVariantSpansFeature = SVAnnotateEngine.variantSpansFeature(variantInterval, featureInterval);
        Assert.assertEquals(actualVariantSpansFeature, expectedVariantSpansFeature);

        final int actualNumBreakpointsInsideFeature =
                SVAnnotateEngine.countBreakendsInsideFeature(variantInterval, featureInterval);
        Assert.assertEquals(actualNumBreakpointsInsideFeature, expectedNumBreakpointsInsideFeature);
    }

    // Toy transcription start sites and initTree() method for testing annotateNearestTranscriptionStartSite()
    private static final SVInterval[] transcriptionStartSites = {
            new SVInterval(0, 100, 101),
            new SVInterval(1, 150, 151),
            new SVInterval(1, 200, 201),
            new SVInterval(1, 250, 251),
            new SVInterval(2, 1, 2)
    };

    // initialize toy TSS SVIntervalTree to test nearest TSS annotation
    private static SVIntervalTree<String> initTree() {
        final SVIntervalTree<String> tree = new SVIntervalTree<>();
        final String[] genes = {"A", "B", "C", "D", "E"};
        if (transcriptionStartSites.length != genes.length) {
            throw new TestException("Transcription start sites list and genes array are not the same length");
        }
        for ( int idx = 0; idx < genes.length; ++idx ) {
            tree.put(transcriptionStartSites[idx], genes[idx]);
        }
        return tree;
    }

    // Toy variants and expected nearest TSS for testing annotateNearestTranscriptionStartSite()
    @DataProvider(name="nearestTSS")
    public Object[][] getNearestTSSData() {
        return new Object[][] {
                { new SimpleInterval("chr1", 1, 50), "A" },
                { new SimpleInterval("chr2", 105, 110), "B" },  // Before first gene on contig
                { new SimpleInterval("chr2", 155, 160), "B" },  // Closer left
                { new SimpleInterval("chr2", 160, 195), "C" },  // Closer right
                { new SimpleInterval("chr2", 3000, 4000), "D" },  // After last gene on contig
                { new SimpleInterval("chr3", 33, 33), "E" },
                { new SimpleInterval("chr4", 900, 4000), null },  // On different contig from TSS data
        };
    }

    // Test annotateNearestTranscriptionStartSite() with toy data
    @Test(dataProvider = "nearestTSS")
    public void testAnnotateNearestTranscriptionStartSite(
            final SimpleInterval variantInterval,
            final String expectedNearestTSSGene)
    {
        final SAMSequenceDictionary sequenceDictionary =
                SVAnnotateUnitTest.createSequenceDictionary(Arrays.asList("chr1", "chr2", "chr3", "chr4"));
        final SVIntervalTree<String> transcriptionStartSiteTree = initTree();
        SVAnnotateEngine svAnnotateEngine = new SVAnnotateEngine(
                new SVAnnotateEngine.GTFIntervalTreesContainer(null, null, transcriptionStartSiteTree),
                null, sequenceDictionary, -1);
        final Map<String, Set<String>> variantConsequenceDict = new HashMap<>();
        svAnnotateEngine.annotateNearestTranscriptionStartSite(variantInterval, variantConsequenceDict);
        if (expectedNearestTSSGene == null) {
            Assert.assertNull(variantConsequenceDict.get(GATKSVVCFConstants.NEAREST_TSS));
        } else {
            Assert.assertEquals(variantConsequenceDict.get(GATKSVVCFConstants.NEAREST_TSS),
                    new HashSet<>(Arrays.asList(expectedNearestTSSGene)));
        }
    }


    /**
     * Get the first transcript from the first gene in the toy GTF
     * @return - first transcript from the first gene in the toy GTF
     */
    private GencodeGtfTranscriptFeature loadToyGtfTranscript() {
        final FeatureDataSource<GencodeGtfGeneFeature> toyGTFSource = SVAnnotateUnitTest.loadToyGTFSource(TOY_GTF_FILE);
        // get only first gene, EMMA1, which has only one transcript
        final GencodeGtfGeneFeature toyGene = toyGTFSource.iterator().next();
        return toyGene.getTranscripts().get(0);
    }

    // Toy variants and protein-coding annotations if the interval is a DUP, DEL, CNV, or INV
    @DataProvider(name = "toyIntervalVariants")
    public Object[][] getToyIntervalVariantTestData() {
        return new Object[][] {
                { new SimpleInterval("chr1", 350, 550), GATKSVVCFConstants.LOF, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.LOF },
                { new SimpleInterval("chr1", 320, 360), GATKSVVCFConstants.LOF, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.LOF },
                { new SimpleInterval("chr1", 550, 750), GATKSVVCFConstants.LOF, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.LOF },
                { new SimpleInterval("chr1", 450, 650), GATKSVVCFConstants.INT_EXON_DUP, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.LOF },
                { new SimpleInterval("chr1", 250, 850), GATKSVVCFConstants.INT_EXON_DUP, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.LOF },
                { new SimpleInterval("chr1", 150, 450), GATKSVVCFConstants.INT_EXON_DUP, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.LOF },
                { new SimpleInterval("chr1", 250, 350), GATKSVVCFConstants.PARTIAL_EXON_DUP, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.LOF },
                { new SimpleInterval("chr1", 350, 650), GATKSVVCFConstants.PARTIAL_EXON_DUP, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.LOF },
                { new SimpleInterval("chr1", 150, 250), GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR },
                { new SimpleInterval("chr1", 650, 850), GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR },
                { new SimpleInterval("chr1", 150, 750), GATKSVVCFConstants.INT_EXON_DUP, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.LOF },
                { new SimpleInterval("chr1", 220, 280), GATKSVVCFConstants.INTRONIC, GATKSVVCFConstants.INTRONIC, GATKSVVCFConstants.INTRONIC, GATKSVVCFConstants.INTRONIC },
                { new SimpleInterval("chr1", 750, 950), GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR },
                { new SimpleInterval("chr1", 50, 150), GATKSVVCFConstants.TSS_DUP, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.LOF },
                { new SimpleInterval("chr1", 50, 250), GATKSVVCFConstants.TSS_DUP, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.LOF },
                { new SimpleInterval("chr1", 50, 550), GATKSVVCFConstants.TSS_DUP, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.LOF },
                { new SimpleInterval("chr1", 950, 1500), GATKSVVCFConstants.DUP_PARTIAL, GATKSVVCFConstants.UTR, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.UTR },
                { new SimpleInterval("chr1", 650, 950), GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR },
                { new SimpleInterval("chr1", 50, 2000), GATKSVVCFConstants.COPY_GAIN, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.INV_SPAN },
                { new SimpleInterval("chr1", 750, 790), GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR },
                { new SimpleInterval("chr1", 350, 1200), GATKSVVCFConstants.DUP_PARTIAL, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.LOF},
                { new SimpleInterval("chr1", 450, 1300), GATKSVVCFConstants.DUP_PARTIAL, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.LOF}
        };
    }

    // Test annotation of interval variants (DUP, DEL, CNV, INV)
    @Test(dataProvider = "toyIntervalVariants")
    public void testAnnotateIntervalSVTypes(
            final SimpleInterval toyVariant,
            final String expectedDuplicationConsequence,
            final String expectedDeletionConsequence,
            final String expectedCopyNumberVariantConsequence,
            final String expectedInversionConsequence
    ) {
        final GencodeGtfTranscriptFeature toyTranscript = loadToyGtfTranscript();

        final String actualDuplicationConsequence = SVAnnotateEngine.annotateDuplication(toyVariant, toyTranscript, false);
        Assert.assertEquals(actualDuplicationConsequence, expectedDuplicationConsequence);

        final String actualDeletionConsequence = SVAnnotateEngine.annotateDeletion(toyVariant, toyTranscript);
        Assert.assertEquals(actualDeletionConsequence, expectedDeletionConsequence);

        final String actualCopyNumberVariantConsequence =
                SVAnnotateEngine.annotateCopyNumberVariant(toyVariant, toyTranscript,
                        SVAnnotateEngine.MSV_EXON_OVERLAP_CLASSIFICATIONS);
        Assert.assertEquals(actualCopyNumberVariantConsequence, expectedCopyNumberVariantConsequence);

        final String actualInversionConsequence = SVAnnotateEngine.annotateInversion(toyVariant, toyTranscript);
        Assert.assertEquals(actualInversionConsequence, expectedInversionConsequence);
    }

    // Toy "point" variants with 2-bp interval (INS, sometimes CTX) and 1-bp interval (BND, sometimes CTX)
    // and expected consequences when annotated as INS, BND, or CTX
    @DataProvider(name = "toyPointVariants")
    public Object[][] getToyPointVariantTestData() {
        return new Object[][] {
                // in UTR
                { new SimpleInterval("chr1", 150, 151),
                        new SimpleInterval("chr1", 150, 150),
                        GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR, GATKSVVCFConstants.LOF },
                // in intron
                { new SimpleInterval("chr1", 250, 251),
                        new SimpleInterval("chr1", 250, 250),
                        GATKSVVCFConstants.INTRONIC, GATKSVVCFConstants.INTRONIC, GATKSVVCFConstants.LOF },
                // in CDS
                { new SimpleInterval("chr1", 350, 351),
                        new SimpleInterval("chr1", 350, 350),
                        GATKSVVCFConstants.LOF, GATKSVVCFConstants.BREAKEND_EXON, GATKSVVCFConstants.LOF },
                // overlap last base of a feature (CDS) only
                { new SimpleInterval("chr1", 600, 601),
                        new SimpleInterval("chr1", 600, 600),
                        GATKSVVCFConstants.LOF, GATKSVVCFConstants.BREAKEND_EXON, GATKSVVCFConstants.LOF },
                // overlap first base of a feature (UTR) only
                { new SimpleInterval("chr1", 100, 101),
                        new SimpleInterval("chr1", 100, 100),
                        GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR, GATKSVVCFConstants.LOF }
        };
    }

    // Test annotation of "point" variants (INS, BND, CTX)
    @Test(dataProvider = "toyPointVariants")
    public void testAnnotatePointSVTypes(
            final SimpleInterval toyTwoBaseVariant,
            final SimpleInterval toyPointVariant,
            final String expectedInsertionConsequence,
            final String expectedBreakendConsequence,
            final String expectedTranslocationVariantConsequence
    ) {
        final GencodeGtfTranscriptFeature toyTranscript = loadToyGtfTranscript();

        final String actualInsertionConsequence = SVAnnotateEngine.annotateInsertion(toyTwoBaseVariant, toyTranscript);
        Assert.assertEquals(actualInsertionConsequence, expectedInsertionConsequence);

        // BND and CTX are annotated one breakpoint at a time
        final String actualBreakendConsequence = SVAnnotateEngine.annotateBreakend(toyPointVariant, toyTranscript);
        Assert.assertEquals(actualBreakendConsequence, expectedBreakendConsequence);

        final String actualPointTranslocationConsequence =
                SVAnnotateEngine.annotateTranslocation(toyPointVariant, toyTranscript);
        Assert.assertEquals(actualPointTranslocationConsequence, expectedTranslocationVariantConsequence);

        final String actualTwoBaseTranslocationConsequence =
                SVAnnotateEngine.annotateTranslocation(toyTwoBaseVariant, toyTranscript);
        Assert.assertEquals(actualTwoBaseTranslocationConsequence, expectedTranslocationVariantConsequence);
    }

    /**
     * Create list of SV segments with SAME SVTYPE - convenience function for testing getSVSegments
     * @param svType - SV type for all segments
     * @param intervals - list of intervals
     * @return - list of SV segments with provided SV type, one for each interval
     */
    private List<SVAnnotateEngine.SVSegment> createListOfSVSegments(final GATKSVVCFConstants.StructuralVariantAnnotationType svType,
                                                                    final SimpleInterval[] intervals) {
        final List<SVAnnotateEngine.SVSegment> segments = new ArrayList<>(intervals.length);
        for (final SimpleInterval interval : intervals) {
            segments.add(new SVAnnotateEngine.SVSegment(svType, interval));
        }
        return segments;
    }

    /**
     * Create list of SV segments with different SVTYPEs - convenience function
     * @param svTypes - list of SV types
     * @param intervals - list of intervals
     * @return - list of SV segments
     */
    private List<SVAnnotateEngine.SVSegment> createListOfSVSegmentsDifferentTypes(final GATKSVVCFConstants.StructuralVariantAnnotationType[] svTypes,
                                                                                  final SimpleInterval[] intervals) {
        Assert.assertEquals(svTypes.length, intervals.length);
        final List<SVAnnotateEngine.SVSegment> segments = new ArrayList<>(intervals.length);
        for (int i = 0; i < svTypes.length; i++) {
            segments.add(new SVAnnotateEngine.SVSegment(svTypes[i], intervals[i]));
        }
        return segments;
    }

    // CPX_TYPE & CPX_INTERVALS INFO fields specifying complex variant intervals, and expected SV segments
    @DataProvider(name = "complexVariantIntervals")
    public Object[][] getComplexVariantIntervalsTestData() {
        return new Object[][] {
                { GATKSVVCFConstants.ComplexVariantSubtype.dDUP, "DUP_chr1:280-420",
                        createListOfSVSegmentsDifferentTypes(new GATKSVVCFConstants.StructuralVariantAnnotationType[]{
                                GATKSVVCFConstants.StructuralVariantAnnotationType.DUP },
                                new SimpleInterval[]{new SimpleInterval("chr1", 280, 420)}) },
                { GATKSVVCFConstants.ComplexVariantSubtype.dupINV, "DUP_chr1:44355904-44356327,INV_chr1:44355904-44357498",
                        createListOfSVSegmentsDifferentTypes(new GATKSVVCFConstants.StructuralVariantAnnotationType[]{
                                GATKSVVCFConstants.StructuralVariantAnnotationType.DUP,
                                        GATKSVVCFConstants.StructuralVariantAnnotationType.INV},
                        new SimpleInterval[]{new SimpleInterval("chr1", 44355904, 44356327),
                                new SimpleInterval("chr1", 44356327, 44357498)}) },
                { GATKSVVCFConstants.ComplexVariantSubtype.dupINVdup, "DUP_chr1:247660974-248129213,INV_chr1:247660974-248587216,DUP_chr1:248520217-248587216",
                        createListOfSVSegmentsDifferentTypes(new GATKSVVCFConstants.StructuralVariantAnnotationType[]{
                                GATKSVVCFConstants.StructuralVariantAnnotationType.DUP,
                                        GATKSVVCFConstants.StructuralVariantAnnotationType.DUP,
                                        GATKSVVCFConstants.StructuralVariantAnnotationType.INV},
                        new SimpleInterval[]{new SimpleInterval("chr1", 247660974, 248129213),
                                new SimpleInterval("chr1", 248520217, 248587216),
                                new SimpleInterval("chr1", 248129213, 248520217)}) },
                { GATKSVVCFConstants.ComplexVariantSubtype.dDUP_iDEL, "DUP_chr2:131488885-131489335,INV_chr2:131488885-131489335,DEL_chr2:130185450-130185720",
                        createListOfSVSegmentsDifferentTypes(new GATKSVVCFConstants.StructuralVariantAnnotationType[]{
                                GATKSVVCFConstants.StructuralVariantAnnotationType.DUP,
                                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL},
                        new SimpleInterval[]{new SimpleInterval("chr2", 131488885,131489335),
                                new SimpleInterval("chr2", 130185450,130185720)}) },
                { GATKSVVCFConstants.ComplexVariantSubtype.dDUP_iDEL, "DUP_chr3:95751919-95752156,DEL_chr3:95746923-95749272",
                        createListOfSVSegmentsDifferentTypes(new GATKSVVCFConstants.StructuralVariantAnnotationType[]{
                                GATKSVVCFConstants.StructuralVariantAnnotationType.DUP,
                                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL},
                        new SimpleInterval[]{new SimpleInterval("chr3", 95751919,95752156),
                                new SimpleInterval("chr3", 95746923,95749272)}) },
                { GATKSVVCFConstants.ComplexVariantSubtype.INS_iDEL, "DEL_chr3:60521333-60521483",
                        createListOfSVSegmentsDifferentTypes(new GATKSVVCFConstants.StructuralVariantAnnotationType[]{
                                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL},
                        new SimpleInterval[]{new SimpleInterval("chr3", 60521333,60521483)}) },
                { GATKSVVCFConstants.ComplexVariantSubtype.delINV, "DEL_chr2:120379742-120383130,INV_chr2:120383130-120384123",
                        createListOfSVSegmentsDifferentTypes(new GATKSVVCFConstants.StructuralVariantAnnotationType[]{
                                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                                GATKSVVCFConstants.StructuralVariantAnnotationType.INV},
                        new SimpleInterval[]{new SimpleInterval("chr2", 120379742,120383130),
                                new SimpleInterval("chr2", 120383130,120384123)}) },
                { GATKSVVCFConstants.ComplexVariantSubtype.INVdel, "INV_chr2:122719025-122719803,DEL_chr2:122719803-122724929",
                        createListOfSVSegmentsDifferentTypes(new GATKSVVCFConstants.StructuralVariantAnnotationType[]{
                                GATKSVVCFConstants.StructuralVariantAnnotationType.INV,
                                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL},
                        new SimpleInterval[]{new SimpleInterval("chr2", 122719025,122719803),
                                new SimpleInterval("chr2", 122719803,122724929)}) },
                { GATKSVVCFConstants.ComplexVariantSubtype.delINVdup, "DEL_chr2:54002577-54003019,INV_chr2:54003019-54006204,DUP_chr2:54006057-54006204",
                        createListOfSVSegmentsDifferentTypes(new GATKSVVCFConstants.StructuralVariantAnnotationType[]{
                                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                                GATKSVVCFConstants.StructuralVariantAnnotationType.DUP,
                                        GATKSVVCFConstants.StructuralVariantAnnotationType.INV},
                        new SimpleInterval[]{new SimpleInterval("chr2", 54002577,54003019),
                                new SimpleInterval("chr2", 54006057,54006204),
                                new SimpleInterval("chr2", 54003019,54006057)}) },
                { GATKSVVCFConstants.ComplexVariantSubtype.dupINVdel, "DUP_chr2:4157678-4157846,INV_chr2:4157678-4165085,DEL_chr2:4165085-4175888",
                        createListOfSVSegmentsDifferentTypes(new GATKSVVCFConstants.StructuralVariantAnnotationType[]{
                                GATKSVVCFConstants.StructuralVariantAnnotationType.DUP,
                                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                                GATKSVVCFConstants.StructuralVariantAnnotationType.INV},
                        new SimpleInterval[]{new SimpleInterval("chr2", 4157678,4157846),
                                new SimpleInterval("chr2", 4165085,4175888),
                                new SimpleInterval("chr2", 4157846,4165085)}) },
                { GATKSVVCFConstants.ComplexVariantSubtype.delINVdel, "DEL_chr2:62384663-62387814,INV_chr2:62387814-62388322,DEL_chr2:62388322-62390272",
                        createListOfSVSegmentsDifferentTypes(new GATKSVVCFConstants.StructuralVariantAnnotationType[]{
                                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                                GATKSVVCFConstants.StructuralVariantAnnotationType.INV,
                                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL},
                        new SimpleInterval[]{new SimpleInterval("chr2", 62384663,62387814),
                                new SimpleInterval("chr2", 62387814,62388322),
                                new SimpleInterval("chr2", 62388322,62390272)}) }
        };
    }

    // Test getComplexAnnotationIntervals()
    @Test(dataProvider = "complexVariantIntervals")
    public void testGetSegmentsFromComplexIntervals(
            final GATKSVVCFConstants.ComplexVariantSubtype complexType,
            final String cpxIntervalsString,
            final List<SVAnnotateEngine.SVSegment> expectedSVSegments
    ) {
        final List<SVAnnotateEngine.SVSegment> cpxIntervals = SVAnnotateEngine.parseComplexIntervals(Arrays.asList(cpxIntervalsString.split(",")));
        final List<SVAnnotateEngine.SVSegment> actualSegments = SVAnnotateEngine.getComplexAnnotationIntervals(cpxIntervals,
                complexType);
        assertSegmentListEqual(actualSegments, expectedSVSegments);
    }

    // CPX_INTERVALS INFO field string specifying complex variant intervals, and expected annotation(s)
    @DataProvider(name = "toyComplexVariants")
    public Object[][] getToyComplexVariantTestData() {
        return new Object[][] {
                { GATKSVVCFConstants.ComplexVariantSubtype.dDUP, "DUP_chr1:280-420", Sets.newHashSet(GATKSVVCFConstants.PARTIAL_DISPERSED_DUP) },
                { GATKSVVCFConstants.ComplexVariantSubtype.INVdup, "INV_chr1:90-1010,DUP_chr1:890-1010",
                        Sets.newHashSet(GATKSVVCFConstants.LOF, GATKSVVCFConstants.PARTIAL_DISPERSED_DUP) },
                { GATKSVVCFConstants.ComplexVariantSubtype.delINVdup, "DEL_chr1:250-450,INV_chr1:450-650,DUP_chr1:610-650",
                        Sets.newHashSet(GATKSVVCFConstants.LOF, GATKSVVCFConstants.PARTIAL_DISPERSED_DUP) }
        };
    }

    // Test annotation of CPX events for gene overlaps from CPX_INTERVALS string
    @Test(dataProvider = "toyComplexVariants")
    public void testAnnotateComplexEvents(
            final GATKSVVCFConstants.ComplexVariantSubtype complexType,
            final String cpxIntervalsString,
            final Set<String> expectedConsequences
    ) {
        final FeatureDataSource<GencodeGtfGeneFeature> toyGTFSource = SVAnnotateUnitTest.loadToyGTFSource(TOY_GTF_FILE);
        final Map<String, Set<String>> variantConsequenceDict = new HashMap<>();
        final SAMSequenceDictionary sequenceDictionary = SVAnnotateUnitTest.createSequenceDictionary(Arrays.asList("chr1"));
        final int promoterWindow = 1000;

        final SVAnnotateEngine.GTFIntervalTreesContainer gtfTrees =
                SVAnnotate.buildIntervalTreesFromGTF(toyGTFSource, sequenceDictionary, promoterWindow);
        SVAnnotateEngine svAnnotateEngine = new SVAnnotateEngine(gtfTrees, null, sequenceDictionary,
                -1);

        final List<SVAnnotateEngine.SVSegment> cpxIntervals = SVAnnotateEngine.parseComplexIntervals(Arrays.asList(cpxIntervalsString.split(",")));
        final List<SVAnnotateEngine.SVSegment> cpxSegments = SVAnnotateEngine.getComplexAnnotationIntervals(cpxIntervals, complexType);
        for (final SVAnnotateEngine.SVSegment cpxSegment: cpxSegments) {
            svAnnotateEngine.annotateGeneOverlaps(cpxSegment.getInterval(), cpxSegment.getIntervalSVType(), true,
                    variantConsequenceDict);
        }

        Assert.assertEquals(variantConsequenceDict.keySet(), expectedConsequences);
    }

    private Map<String, Set<String>> createVariantConsequenceDict(final String[] consequences,
                                                                  final String[] features) {
        Assert.assertEquals(consequences.length, features.length);
        final Map<String, Set<String>> map = new HashMap<>();
        for (int i = 0; i < consequences.length; i++) {
            SVAnnotateEngine.updateVariantConsequenceDict(map, consequences[i], features[i]);
        }
        return map;
    }


    // Test sortVariantConsequenceDict() sorts lists of genes in variant consequence map
    @Test
    public void testSortVariantConsequenceDict() {
        final Map<String, Set<String>> before = createVariantConsequenceDict(
                new String[]{GATKSVVCFConstants.LOF, GATKSVVCFConstants.LOF, GATKSVVCFConstants.LOF,
                        GATKSVVCFConstants.LOF, GATKSVVCFConstants.DUP_PARTIAL, GATKSVVCFConstants.LOF,
                        GATKSVVCFConstants.TSS_DUP},
                new String[]{"NOC2L", "KLHL17", "PLEKHN1", "PERM1", "SAMD11", "HES4", "ISG15"}
        );

        final Map<String, Object> expectedAfter = new HashMap<>();
        expectedAfter.put(GATKSVVCFConstants.DUP_PARTIAL, Arrays.asList("SAMD11"));
        expectedAfter.put(GATKSVVCFConstants.TSS_DUP, Arrays.asList("ISG15"));
        expectedAfter.put(GATKSVVCFConstants.LOF, Arrays.asList("HES4", "KLHL17", "NOC2L", "PERM1", "PLEKHN1"));

        Assert.assertEquals(SVAnnotateEngine.sortVariantConsequenceDict(before), expectedAfter);
    }

    @DataProvider(name = "proteinCodingConsequences")
    public Object[][] getProteinCodingConsequences() {
        return new Object[][] {
                // protein-coding consequence -> not intergenic
                { createVariantConsequenceDict(
                        new String[]{GATKSVVCFConstants.TSS_DUP, GATKSVVCFConstants.PARTIAL_DISPERSED_DUP},
                        new String[]{"HES4", "SAMD11"}),
                        false
                },
                // partial dispersed dup does not count towards protein-coding consequences for this
                { createVariantConsequenceDict(
                        new String[]{GATKSVVCFConstants.PARTIAL_DISPERSED_DUP},
                        new String[]{"SAMD11"}),
                        true
                },
                // ignore noncoding consequences
                { createVariantConsequenceDict(
                        new String[]{GATKSVVCFConstants.NONCODING_BREAKPOINT, GATKSVVCFConstants.PROMOTER},
                        new String[]{"Enhancer", "RP1"}),
                        true
                },
                { new HashMap<>(),
                        true
                }
        };
    }

    @Test(dataProvider = "proteinCodingConsequences")
    public void testIsIntergenic(
            final Map<String, Set<String>> variantConsequenceDict,
            final boolean expectedIsIntergenic
    ){
        final boolean actualIsIntergenic = SVAnnotateEngine.isIntergenic(variantConsequenceDict,
                SVAnnotateEngine.PROTEIN_CODING_CONSEQUENCES);
        Assert.assertEquals(actualIsIntergenic, expectedIsIntergenic);
    }

    @DataProvider(name = "complexSubtypes")
    public Object[][] getComplexSubtypes() {
        return new Object[][] {
                { null, false },
                { GATKSVVCFConstants.ComplexVariantSubtype.delINV, false },
                { GATKSVVCFConstants.ComplexVariantSubtype.INVdel, false },
                { GATKSVVCFConstants.ComplexVariantSubtype.dupINV, true },
                { GATKSVVCFConstants.ComplexVariantSubtype.INVdup, true },
                { GATKSVVCFConstants.ComplexVariantSubtype.delINVdel, false },
                { GATKSVVCFConstants.ComplexVariantSubtype.dupINVdup, true },
                { GATKSVVCFConstants.ComplexVariantSubtype.delINVdup, true },
                { GATKSVVCFConstants.ComplexVariantSubtype.dupINVdel, true },
                { GATKSVVCFConstants.ComplexVariantSubtype.dDUP, true },
                { GATKSVVCFConstants.ComplexVariantSubtype.dDUP_iDEL, true },
                { GATKSVVCFConstants.ComplexVariantSubtype.INS_iDEL, false },
                { GATKSVVCFConstants.ComplexVariantSubtype.CTX_PP_QQ, false },
                { GATKSVVCFConstants.ComplexVariantSubtype.CTX_PQ_QP, false },
                { GATKSVVCFConstants.ComplexVariantSubtype.CTX_INV, false }
        };
    }

    @Test(dataProvider = "complexSubtypes")
    public void testIncludesDispersedDuplication(
            final GATKSVVCFConstants.ComplexVariantSubtype complexType,
            final boolean expectedIncludesDispersedDuplication
    ){
        final boolean actualIncludesDispersedDuplication = SVAnnotateEngine.includesDispersedDuplication(complexType,
                SVAnnotateEngine.COMPLEX_SUBTYPES_WITH_DISPERSED_DUP);
        Assert.assertEquals(actualIncludesDispersedDuplication, expectedIncludesDispersedDuplication);
    }


    /**
     * Assert two lists of SVAnnotate.SVSegment objects are equal in contents
     * Lists must be same size and contain equal SVSegments in the same order
     * @param segmentsA - first list of SVSegments to compare
     * @param segmentsB - second list of SVSegments to compare
     */
    private void assertSegmentListEqual(final List<SVAnnotateEngine.SVSegment> segmentsA,
                                        final List<SVAnnotateEngine.SVSegment> segmentsB) {
        final int lengthA = segmentsA.size();
        if (lengthA != segmentsB.size()) {
            Assert.fail("Segment lists differ in length");
        }
        for (int i = 0; i < lengthA; i++) {
            SVAnnotateEngine.SVSegment segmentA = segmentsA.get(i);
            SVAnnotateEngine.SVSegment segmentB = segmentsB.get(i);
            if (!segmentA.equals(segmentB)) {
                Assert.fail("Segment items differ");
            }
        }
    }

    /**
     * Utility function to create a VariantContext object with relevant attributes for testing
     */
    private VariantContext createVariantContext(final String chrom, final int pos, final int end, final String chr2,
                                                final Integer end2, final String ref, final String alt,
                                                final Integer svLen, final String strands, final String cpxType,
                                                final List<String> cpxIntervals) {
        final Map<String, Object> attributes = new HashMap<>();
        if (chr2 != null) {
            attributes.put(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, chr2);
        }
        if (end2 != null) {
            attributes.put(GATKSVVCFConstants.END2_ATTRIBUTE, end2);
        }
        if (svLen != null) {
            attributes.put(GATKSVVCFConstants.SVLEN, svLen);
        }
        if (strands != null) {
            attributes.put(GATKSVVCFConstants.STRANDS_ATTRIBUTE, strands);
        }
        if (cpxType != null) {
            attributes.put(GATKSVVCFConstants.CPX_TYPE, cpxType);
        }
        if (cpxIntervals != null) {
            attributes.put(GATKSVVCFConstants.CPX_INTERVALS, cpxIntervals);
        }
        return new VariantContextBuilder()
                .source("source")
                .id("id")
                .chr(chrom)
                .start(pos)
                .stop(end)
                .alleles(Arrays.asList(ref != null ? Allele.create(ref, true) : Allele.REF_N,
                        alt != null ? Allele.create(alt, false) : Allele.ALT_N))
                .attributes(attributes)
                .make();
    }

    // VariantContext objects with all types and representations to test getSVSegments and getSVType
    @DataProvider(name = "typesAndSegments")
    public Object[][] getSVTypesAndSegmentsTestData() {
        return new Object[][] {
                { createVariantContext("chr2", 86263976, 86263977, "chr19", 424309, "N",
                        "<CTX>", null, null, "CTX_PP/QQ", null),
                        GATKSVVCFConstants.ComplexVariantSubtype.CTX_PP_QQ,
                        GATKSVVCFConstants.StructuralVariantAnnotationType.CTX,
                        createListOfSVSegments(GATKSVVCFConstants.StructuralVariantAnnotationType.CTX,
                                new SimpleInterval[]{ new SimpleInterval("chr2", 86263976, 86263976),
                                        new SimpleInterval("chr2", 86263977, 86263977),
                                        new SimpleInterval("chr19", 424309, 424309),
                                        new SimpleInterval("chr19", 424310, 424310) }),
                        null },
                { createVariantContext("chr2", 86263976, 86263976, null, 424309, "G",
                        "G]chr19:424309]", null, null,"CTX_PP/QQ", null),
                        GATKSVVCFConstants.ComplexVariantSubtype.CTX_PP_QQ,
                        GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                        Arrays.asList(new SVAnnotateEngine.SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.CTX,
                                new SimpleInterval("chr2", 86263976, 86263976))),
                        null},
                { createVariantContext("chr2", 86263977, 86263977, null, 424309, "A",
                        "[chr19:424310[A", -1, null, "CTX_PP/QQ", null),
                        GATKSVVCFConstants.ComplexVariantSubtype.CTX_PP_QQ,
                        GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                        Arrays.asList(new SVAnnotateEngine.SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.CTX,
                                new SimpleInterval("chr2", 86263977, 86263977))),
                        null },
                { createVariantContext("chr4", 21923233, 22506191, "chr8", 107912008, "N",
                        "<CTX>", null, null, "CTX_INV",
                        Collections.singletonList("INV_chr4:21923233-22506191")),
                        GATKSVVCFConstants.ComplexVariantSubtype.CTX_INV,
                        GATKSVVCFConstants.StructuralVariantAnnotationType.CTX,
                        Arrays.asList(new SVAnnotateEngine.SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.INV,
                                        new SimpleInterval("chr4", 21923233, 22506191)),
                                new SVAnnotateEngine.SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.CTX,
                                        new SimpleInterval("chr4", 21923233, 21923233)),
                                new SVAnnotateEngine.SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.CTX,
                                        new SimpleInterval("chr4", 22506191, 22506191)),
                                new SVAnnotateEngine.SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.CTX,
                                        new SimpleInterval("chr8", 107912008, 107912008)),
                                new SVAnnotateEngine.SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.CTX,
                                        new SimpleInterval("chr8", 107912009, 107912009)) ),
                        null },
                { createVariantContext("chr2", 205522308, 205522384, "chr2", null, "N",
                        "<INV>", 76, null, null, null),
                        null,
                        GATKSVVCFConstants.StructuralVariantAnnotationType.INV,
                        Arrays.asList(new SVAnnotateEngine.SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.INV,
                                new SimpleInterval("chr2", 205522308, 205522384))),
                        null },
                { createVariantContext("chr19", 424309, 424309, null, 424309, "T",
                        "T]chr2:86263976]", null, null, "CTX_PP/QQ", null),
                        GATKSVVCFConstants.ComplexVariantSubtype.CTX_PP_QQ,
                        GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                        Arrays.asList(new SVAnnotateEngine.SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.CTX,
                                new SimpleInterval("chr19", 424309, 424309))),
                        null },
                { createVariantContext("chr19", 424310, 424310, null, 424309, "C",
                        "[chr2:86263977[C", null, null, "CTX_PP/QQ", null),
                        GATKSVVCFConstants.ComplexVariantSubtype.CTX_PP_QQ,
                        GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                        Arrays.asList(new SVAnnotateEngine.SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.CTX,
                                new SimpleInterval("chr19", 424310, 424310))),
                        null },
                { createVariantContext("chr22", 10510000, 10694100, "chr22", null, "N",
                        "<DEL>", 184100, null, null, null),
                        null,
                        GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                        createListOfSVSegments(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 10510000, 10694100)}),
                        null},
                { createVariantContext("chr22", 10510000, 10694100, "chr22", null, "N",
                        "<DEL:ME:LINE1>", 184100, null, null, null),
                        null,
                        GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                        createListOfSVSegments(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 10510000, 10694100)}),
                        null},
                { createVariantContext("chr22", 10524000, 10710000, "chr22", null, "N",
                        "<DUP>", 186000, null, null, null),
                        null,
                        GATKSVVCFConstants.StructuralVariantAnnotationType.DUP,
                        createListOfSVSegments(GATKSVVCFConstants.StructuralVariantAnnotationType.DUP,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 10524000, 10710000)}),
                        null },
                { createVariantContext("chr22", 10532563, 10532611, "chr22", null, "N",
                        "<INS:ME:ALU>", 245, null, null, null),
                        null,
                        GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                        createListOfSVSegments(GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 10532563, 10532564)}),
                        null },
                { createVariantContext("chr22", 10572758, 10572788, "chr22", null, "N",
                        "<INS>", 57, null, null, null),
                        null,
                        GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                        createListOfSVSegments(GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 10572758, 10572759)}),
                        null },
                { createVariantContext("chr22", 10717890, 10717890, "chr22", null, "N",
                        "<BND>", 5170, "-+", null, null),
                        null,
                        GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                        createListOfSVSegments(GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 10717890, 10717890),
                                        new SimpleInterval("chr22", 10723060, 10723060) }),
                        createListOfSVSegments(GATKSVVCFConstants.StructuralVariantAnnotationType.DUP,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 10717890, 10723060)}) },
                { createVariantContext("chr22", 10774600, 10784500, "chr22", null, "N",
                        "<CNV>", 9900, null, null, null),
                        null,
                        GATKSVVCFConstants.StructuralVariantAnnotationType.CNV,
                        createListOfSVSegments(GATKSVVCFConstants.StructuralVariantAnnotationType.CNV,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 10774600, 10784500)}),
                        null },
                { createVariantContext("chr22", 10930458, 10930458, "chr22", 11564561, "N",
                        "<BND>", 634103, "--", null, null),
                        null,
                        GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                        createListOfSVSegments(GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 10930458, 10930458),
                                        new SimpleInterval("chr22", 11564561, 11564561) }),
                        null },
                { createVariantContext("chr22", 17636024, 17636024, "chr22", null, "N",
                        "<BND>", 10709, "+-", null, null),
                        null,
                        GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                        createListOfSVSegments(GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 17636024, 17636024),
                                        new SimpleInterval("chr22", 17646733, 17646733) }),
                        createListOfSVSegments(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 17636024, 17646733)}) },
                { createVariantContext("chr22", 18971159, 18971435, "chr22", null, "N",
                        "<CPX>", 386, null, "dDUP",
                        Arrays.asList("INV_chr22:20267228-20267614","DUP_chr22:20267228-20267614")),
                        GATKSVVCFConstants.ComplexVariantSubtype.dDUP,
                        GATKSVVCFConstants.StructuralVariantAnnotationType.CPX,
                        Arrays.asList(
                                new SVAnnotateEngine.SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.DUP,
                                        new SimpleInterval("chr22", 20267228, 20267614)),
                                new SVAnnotateEngine.SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                                        new SimpleInterval("chr22", 18971159, 18971160))),
                        null },
                { createVariantContext("chr22", 22120897, 22120897, "chrX", 126356858, "N",
                        "<BND>", -1, "++", null, null),
                        null,
                        GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                        createListOfSVSegments(GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 22120897, 22120897),
                                        new SimpleInterval("chrX", 126356858, 126356858) }),
                        null },
                { createVariantContext("chr22", 22196261, 22196261, "chr22", null, "N",
                        "<BND>", 708725, "+-", null, null),
                        null,
                        GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                        createListOfSVSegments(GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 22196261, 22196261),
                                        new SimpleInterval("chr22", 22904986, 22904986) }),
                        null },
                { createVariantContext("chr22", 22196261, 22196261, null, null, "A",
                        "A[chr22:22904986[", null, "+-", null, null),
                        null,
                        GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                        createListOfSVSegments(GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 22196261, 22196261) }),
                        null },
                { createVariantContext("chr22", 22904986, 22904986, null, null, "T",
                        "]chr22:22196261]T", null, "+-", null, null),
                        null,
                        GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                        createListOfSVSegments(GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 22904986, 22904986) }),
                        null },
                { createVariantContext("chr22", 36533058, 36538234, "chr22", null, "N",
                        "<CPX>", 5176, null, "dupINV",
                        Arrays.asList("DUP_chr22:36533058-36533299","INV_chr22:36533058-36538234")),
                        GATKSVVCFConstants.ComplexVariantSubtype.dupINV,
                        GATKSVVCFConstants.StructuralVariantAnnotationType.CPX,
                        Arrays.asList(
                                new SVAnnotateEngine.SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.DUP,
                                        new SimpleInterval("chr22", 36533058, 36533299)),
                                new SVAnnotateEngine.SVSegment(GATKSVVCFConstants.StructuralVariantAnnotationType.INV,
                                        new SimpleInterval("chr22", 36533299, 36538234))),
                        null }
        };
    }

    // Test getSVType() and getSVSegments() on representative list of VariantContext objects
    @Test(dataProvider = "typesAndSegments")
    public void testGetSVTypeAndSegments(
            final VariantContext variant,
            final GATKSVVCFConstants.ComplexVariantSubtype complexType,
            final GATKSVVCFConstants.StructuralVariantAnnotationType expectedSVType,
            final List<SVAnnotateEngine.SVSegment> expectedSVSegments,
            final List<SVAnnotateEngine.SVSegment> expectedSVSegmentsWithBNDOverlap
    ) {
        GATKSVVCFConstants.StructuralVariantAnnotationType actualSVType = SVAnnotateEngine.getSVType(variant);
        Assert.assertEquals(actualSVType, expectedSVType);

        final List<SVAnnotateEngine.SVSegment> actualSegments = SVAnnotateEngine.getSVSegments(variant,
                actualSVType, -1, complexType);
        assertSegmentListEqual(actualSegments, expectedSVSegments);

        final List<SVAnnotateEngine.SVSegment> actualSegmentsWithBNDOverlap = SVAnnotateEngine.getSVSegments(variant,
                actualSVType, 15000, complexType);
        assertSegmentListEqual(actualSegmentsWithBNDOverlap,
                expectedSVSegmentsWithBNDOverlap != null ? expectedSVSegmentsWithBNDOverlap : expectedSVSegments);
    }

    /**
     * Utility function to create key -> value Map for each key,value in lists
     * For testing INFO field annotations
     * Visible for use in integration tests as well
     * @param keys - list of keys
     * @param values - list of values
     * @return - key -> value map for each key,value in input lists, paired based on order of lists
     */
    protected static Map<String, Object> createAttributesMap(final List<String> keys, final List<Object> values) {
        final Map<String, Object> attributes = new HashMap<>();
        final int len = keys.size();
        if (len != values.size()) {
            throw new TestException("Length of keys list != length of values list");
        }
        for (int i = 0; i < len; i++) {
            final Object value = values.get(i);
            attributes.put(keys.get(i), value.getClass() == String.class ? Arrays.asList(value) : value);
        }
        return attributes;
    }

    // VariantContext objects with expected annotations to test full process
    // with focus on noncoding, promoter, and TSS annotation not otherwise covered
    @DataProvider(name = "toyIntegratedVariants")
    public Object[][] getAnnotateStructuralVariantTestData() {
        return new Object[][]{
                { createVariantContext("chr1", 10, 50, null, null, null,
                        "<CNV>", 40, null, null, null),
                        createAttributesMap(
                                Arrays.asList(GATKSVVCFConstants.PROMOTER, GATKSVVCFConstants.INTERGENIC),
                                Arrays.asList("EMMA1", true)) },
                // CPX with dDUP CG and promoter annotations
                { createVariantContext("chr1", 10, 11, null, null, null,
                        "<CPX>", 1500, null, "dDUP",
                        Collections.singletonList("DUP_chr1:2000-3500")),
                        createAttributesMap(
                                Arrays.asList(GATKSVVCFConstants.PROMOTER, GATKSVVCFConstants.COPY_GAIN,
                                        GATKSVVCFConstants.INTERGENIC),
                                Arrays.asList("EMMA1", "EMMA2", false)) },
                // dupINVdup with CG
                { createVariantContext("chr1", 20, 3500, null, null, null,
                        "<CPX>", 3480, null, "dupINVdup",
                        Arrays.asList("DUP_chr1:20-70", "INV_chr1:20-3500", "DUP_chr1:2000-3500")),
                        createAttributesMap(
                                Arrays.asList(GATKSVVCFConstants.INV_SPAN, GATKSVVCFConstants.COPY_GAIN,
                                        GATKSVVCFConstants.NONCODING_SPAN, GATKSVVCFConstants.INTERGENIC),
                                Arrays.asList("EMMA1", "EMMA2", Arrays.asList("DNase", "Enhancer"), false)) },
                // ignore INV for dDUP; ignore dDUP for promoter; CPX noncoding span; CPX intergenic
                { createVariantContext("chr1", 1101, 1102, null, null, null,
                        "<CPX>", 700, null, "dDUP_iDEL",
                        Arrays.asList("DUP_chr1:10-60", "INV_chr1:10-60","DEL_chr1:1100-1700")),
                        createAttributesMap(
                                Arrays.asList(GATKSVVCFConstants.NONCODING_SPAN, GATKSVVCFConstants.NEAREST_TSS,
                                        GATKSVVCFConstants.INTERGENIC),
                                Arrays.asList("Enhancer", "EMMA1", true)) },
                // Ignore INV and source INS in CPX_INTERVALS for INS_iDEL
                { createVariantContext("chr1", 1101, 1102, null, null, null,
                        "<CPX>", 700, null, "INS_iDEL",
                        Arrays.asList("INV_chr1:10-60","DEL_chr1:1100-1700","INS_chr1:10-60")),
                        createAttributesMap(
                                Arrays.asList(GATKSVVCFConstants.NONCODING_SPAN, GATKSVVCFConstants.NEAREST_TSS,
                                        GATKSVVCFConstants.INTERGENIC),
                                Arrays.asList("Enhancer", "EMMA1", true)) },
                // no noncoding breakpoint for DUP segment of CPX; modify INV interval; multiple genes for a consequence
                { createVariantContext("chr1", 450, 3100, null, null, null,
                        "<CPX>", 2650, null, "dupINV",
                        Arrays.asList("DUP_chr1:450-2200", "INV_chr1:450-3100")),
                        createAttributesMap(
                                Arrays.asList(GATKSVVCFConstants.PARTIAL_DISPERSED_DUP, GATKSVVCFConstants.LOF,
                                        GATKSVVCFConstants.INTERGENIC),
                                Arrays.asList(Arrays.asList("EMMA1", "EMMA2"), "EMMA2", false)) },
                // intergenic with partial dispersed dup
                { createVariantContext("chr1", 1100, 1100, null, null, null,
                        "<CPX>", 550, null, "dDUP",
                        Arrays.asList("DUP_chr1:450-1000")),
                        createAttributesMap(
                                Arrays.asList(GATKSVVCFConstants.PARTIAL_DISPERSED_DUP, GATKSVVCFConstants.NEAREST_TSS,
                                        GATKSVVCFConstants.INTERGENIC),
                                Arrays.asList("EMMA1", "EMMA1", true)) },
                // merge INV + DEL for nearest TSS; noncoding breakpoint for CPX
                { createVariantContext("chr1", 1400, 1900, null, null, null,
                        "<CPX>", 500, null, "delINV",
                        Arrays.asList("DEL_chr1:1400-1500", "INV_chr1:1500-1900")),
                        createAttributesMap(
                                Arrays.asList(GATKSVVCFConstants.NONCODING_BREAKPOINT, GATKSVVCFConstants.NEAREST_TSS,
                                        GATKSVVCFConstants.INTERGENIC),
                                Arrays.asList("Enhancer", "EMMA2", true)) },
                { createVariantContext("chr1", 50, 450, null, null, null,
                        "<DEL>", 400, null, null, null),
                        createAttributesMap(
                                Arrays.asList(GATKSVVCFConstants.LOF, GATKSVVCFConstants.NONCODING_BREAKPOINT,
                                        GATKSVVCFConstants.INTERGENIC),
                                Arrays.asList("EMMA1", "DNase", false)) },
                { createVariantContext("chr1", 1100, 1700, null, null, null,
                        "<INV>", 600, null, null, null),
                        createAttributesMap(
                                Arrays.asList(GATKSVVCFConstants.NONCODING_SPAN, GATKSVVCFConstants.NEAREST_TSS,
                                        GATKSVVCFConstants.INTERGENIC),
                                Arrays.asList("Enhancer", "EMMA1", true)) },
                { createVariantContext("chr1", 1300, 1800, null, null, null,
                        "<INS:LINE1>", 42, null, null, null),
                        createAttributesMap(
                                Arrays.asList(GATKSVVCFConstants.NONCODING_BREAKPOINT, GATKSVVCFConstants.NEAREST_TSS,
                                        GATKSVVCFConstants.INTERGENIC),
                                Arrays.asList("Enhancer", "EMMA1", true)) },
                { createVariantContext("chr1", 2000, 2200, null, null, null,
                        "<DUP>", 42, null, null, null),
                        createAttributesMap(
                                Arrays.asList(GATKSVVCFConstants.DUP_PARTIAL,GATKSVVCFConstants.INTERGENIC),
                                Arrays.asList("EMMA2", false)) },
                { createVariantContext("chr1", 3100, 3300, null, null, null,
                        "<DUP>", 200, null, null, null),
                        createAttributesMap(
                                Arrays.asList(GATKSVVCFConstants.PROMOTER,GATKSVVCFConstants.INTERGENIC),
                                Arrays.asList("EMMA2", true)) },
                // CTX_INV - this test case is not interchromosomal because the test GTF only has chr1
                { createVariantContext("chr1", 10, 1010, "chr1", 2300, null,
                        "<CTX>", -1, null, "CTX_INV", Arrays.asList("INV_chr1:10-1010")),
                        createAttributesMap(
                                Arrays.asList(GATKSVVCFConstants.INV_SPAN, GATKSVVCFConstants.NONCODING_SPAN,
                                        GATKSVVCFConstants.LOF, GATKSVVCFConstants.INTERGENIC),
                                Arrays.asList("EMMA1", "DNase", "EMMA2", false))},
                // check annotate promoter for all segments in multi-segment SV
                { createVariantContext("chr1", 30, 30, "chr1", 3030, null,
                        "<BND>", null, "-+", null, null),
                        createAttributesMap(
                                Arrays.asList(GATKSVVCFConstants.PROMOTER, GATKSVVCFConstants.INTERGENIC),
                                Arrays.asList(Arrays.asList("EMMA1", "EMMA2"), true)) }
        };
    }

    /**
     * Tests full integrated annotateStructuralVariant process, with a focus on noncoding & promoter annotation
     * and logic controlling when to annotate TSS and intergenic, since protein-coding and TSS annotation
     * is tested elsewhere
     */
    @Test(dataProvider = "toyIntegratedVariants")
    public void testAnnotateStructuralVariant(
            final VariantContext variant,
            final Map<String, Object> expectedAttributes
    ) {
        final SAMSequenceDictionary sequenceDictionary =
                SVAnnotateUnitTest.createSequenceDictionary(Arrays.asList("chr1"));
        final int promoterWindow = 200;
        final int maxBreakendLen = -1;

        final FeatureDataSource<GencodeGtfGeneFeature> toyGTFSource =
                SVAnnotateUnitTest.loadToyGTFSource(TOY_GTF_FILE);
        final SVAnnotateEngine.GTFIntervalTreesContainer gtfTrees =
                SVAnnotate.buildIntervalTreesFromGTF(toyGTFSource, sequenceDictionary, promoterWindow);

        final FeatureDataSource<FullBEDFeature> tinyNoncodingBedSource =
                SVAnnotateUnitTest.loadTinyNoncodingBEDSource(TINY_NONCODING_BED_FILE);
        final SVIntervalTree<String> nonCodingIntervalTree =
                SVAnnotate.buildIntervalTreeFromBED(tinyNoncodingBedSource, sequenceDictionary);

        SVAnnotateEngine svAnnotateEngine = new SVAnnotateEngine(gtfTrees, nonCodingIntervalTree,
                sequenceDictionary, maxBreakendLen);
        final Map<String, Object> actualAttributes =
                svAnnotateEngine.annotateStructuralVariant(variant);

        Assert.assertEquals(actualAttributes, expectedAttributes);
    }
}
