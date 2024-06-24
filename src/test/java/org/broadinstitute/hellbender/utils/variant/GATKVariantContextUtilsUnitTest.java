package org.broadinstitute.hellbender.utils.variant;

import com.google.common.jimfs.Configuration;
import com.google.common.jimfs.Jimfs;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.FeatureReader;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceMemorySource;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAssignmentMethod;
import org.broadinstitute.hellbender.tools.walkers.variantutils.LeftAlignAndTrimVariants;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.nio.file.FileSystem;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public final class GATKVariantContextUtilsUnitTest extends GATKBaseTest {

    private static final int LENGTH_OF_ARTIFICIAL_READ = 50;

    Allele Aref, Cref, Gref, Tref, A, T, C, G, ATC, ATCATC;
    Allele ATCATCT;
    Allele ATref;
    Allele ATCref;
    Allele Anoref;
    Allele GT;
    SimpleInterval baseLoc = new SimpleInterval("20", 1000, 1000);

    final String refBases1 = "GCAGAGCTGACCCTCCCTCCCCTCTCCCAGTGCAACAGCACGGGCGGCGACTGCTTTTACCGAGGCTACACGTCAGGCGTGGCGGCTGTCCAGGACTGGTACCACTTCCACTATGTGGATCTCTGCTGAGGACCAGGAAAGCCAGCACCCGCAGAGACTCTTCCCCAGTGCTCCATACGATCACCATTCTCTGCAGAAGG";
    final String longStr = "AGTCGCTCGAGCTCGAGCTCGAGTGTGCGCTCTACAGCTCAGCTCGCTCGCACACAT";
    final List<String> longPieces = Arrays.asList("AAAAAAAAAAAAAAAAAAAAAAAAAAAA", "TCTCTCTCTCTCTC", "CAGTCCAGTCCAGTCCAGTCCAGTCCAGTCCAGTCCAGTCCAGTCCAGTCCAGTC", Utils.dupString("A", 120), Utils.dupString("TCACGTC", 60), Utils.dupString(longStr, 2)); // where we'll perform tests
    final List<Integer> strLengths = Arrays.asList(1, 2, 5, 1, 7, 57); //
    final List<String> refBasesStrings = longPieces.stream().map(l -> refBases1 + l + refBases1).collect(Collectors.toList());

    final List<Integer> contigStops = refBasesStrings.stream().map(String::length).collect(Collectors.toList());
    final List<SAMFileHeader> headers = contigStops.stream().map(c ->
            ArtificialReadUtils.createArtificialSamHeader(1, 1, c)).collect(Collectors.toList());
    final String artificialContig = "1";
    final int repeatStart = refBases1.length();
    final List<SimpleInterval> refIntervals = contigStops.stream().map(c ->
            new SimpleInterval(artificialContig, 1, c)).collect(Collectors.toList());
    final List<ReferenceBases> refBases = IntStream.range(0, longPieces.size()).mapToObj(i -> new ReferenceBases(refBasesStrings.get(i).getBytes(), refIntervals.get(i)))
            .collect(Collectors.toList());
    final List<ReferenceMemorySource> refSources = IntStream.range(0, longPieces.size()).mapToObj(i -> new ReferenceMemorySource(refBases.get(i), headers.get(i).getSequenceDictionary()))
            .collect(Collectors.toList());

    @BeforeClass
    public void setup() throws IOException {
        // alleles
        Aref = Allele.create("A", true);
        Cref = Allele.create("C", true);
        Gref = Allele.create("G", true);
        Tref = Allele.create("T", true);
        A = Allele.create("A");
        T = Allele.create("T");
        C = Allele.create("C");
        G = Allele.create("G");
        ATC = Allele.create("ATC");
        ATCref = Allele.create("ATC", true);
        ATCATC = Allele.create("ATCATC");
        ATCATCT = Allele.create("ATCATCT");
        ATref = Allele.create("AT",true);
        Anoref = Allele.create("A",false);
        GT = Allele.create("GT",false);
    }

    private Genotype makeG(String sample, Allele a1, Allele a2, double log10pError, int... pls) {
        return new GenotypeBuilder(sample, Arrays.asList(a1, a2)).log10PError(log10pError).PL(pls).make();
    }


    private Genotype makeG(String sample, Allele a1, Allele a2, double log10pError) {
        return new GenotypeBuilder(sample, Arrays.asList(a1, a2)).log10PError(log10pError).make();
    }

    private VariantContext makeVC(String source, List<Allele> alleles) {
        return makeVC(source, alleles, null, null);
    }

    private VariantContext makeVC(String source, List<Allele> alleles, Genotype... g1) {
        return makeVC(source, alleles, Arrays.asList(g1));
    }

    private VariantContext makeVC(String source, List<Allele> alleles, String filter) {
        return makeVC(source, alleles, filter.equals(".") ? null : new LinkedHashSet<String>(Collections.singletonList(filter)));
    }

    private VariantContext makeVC(String source, List<Allele> alleles, Set<String> filters) {
        return makeVC(source, alleles, null, filters);
    }

    private VariantContext makeVC(String source, List<Allele> alleles, Collection<Genotype> genotypes) {
        return makeVC(source, alleles, genotypes, null);
    }

    private VariantContext makeVC(String source, List<Allele> alleles, Collection<Genotype> genotypes, Set<String> filters) {
        int start = 10;
        int stop = start + alleles.get(0).length() - 1; // alleles.contains(ATC) ? start + 3 : start;
        return new VariantContextBuilder(source, "1", start, stop, alleles).genotypes(genotypes).filters(filters).make();
    }

    private VariantContext makeVC(final List<Allele> alleles, final int start) {
        int stop = start + alleles.get(0).length() - 1;
        return new VariantContextBuilder("source", "1", start, stop, alleles).make();
    }

    private VariantContext makeVC(final Allele allele, final int start) {
        return makeVC(Collections.singletonList(allele), start);
    }

    private VariantContext makeVC(final Allele ref, final Allele alt, final int start) {
        return makeVC(Arrays.asList(ref, alt), start);
    }

    @Test
    public void testHomozygousAlleleList() throws Exception {
        final List<Allele> alleles = GATKVariantContextUtils.homozygousAlleleList(T, 2);
        Assert.assertEquals(alleles, Arrays.asList(T, T));
    }

    final static private Locatable START_AT_1 = new SimpleInterval("1",1,1);
    final static private Locatable START_AT_2 = new SimpleInterval("1",2,2);

    @DataProvider
    Object[][] provideDataForDetermineReferenceAllele() {
        // {list of vcs, locus, expected ref - contig is "1" throughout
        return new Object[][] {
                { Collections.emptyList(), null, null },
                { Collections.emptyList(), START_AT_1, null },
                { Collections.singletonList(makeVC(Aref,1)), START_AT_1, Aref },
                { Collections.singletonList(makeVC(Aref,1)), START_AT_2, null },

                {Arrays.asList(makeVC(Aref,1), makeVC(ATref,1)), START_AT_1, ATref},
                {Arrays.asList(makeVC(ATref,1), makeVC(Aref,1)), START_AT_1, ATref},
                {Arrays.asList(makeVC(Aref,1), makeVC(ATref,1)), null, ATref},
                {Arrays.asList(makeVC(ATref,1), makeVC(Aref,1)), START_AT_2, null},
                {Arrays.asList(makeVC(Aref,1), makeVC(ATref,2)), START_AT_1, Aref},

                {Arrays.asList(makeVC(Aref, C,1), makeVC(ATref,ATCATC,1)), START_AT_1, ATref},
                {Arrays.asList(makeVC(Aref,1), makeVC(ATCref,1), makeVC(ATref,1)), START_AT_1, ATCref},
        };
    }

    @Test(dataProvider = "provideDataForDetermineReferenceAllele")
    public void testDetermineReferenceAllele(final List<VariantContext> vcs, final Locatable loc, final Allele expectedRef) {

        final Allele ref = GATKVariantContextUtils.determineReferenceAllele(vcs, loc);

        if (expectedRef == null) {
            Assert.assertNull(ref);
        } else {
            Assert.assertTrue(ref.isReference());
            Assert.assertEquals(ref.getBaseString(), expectedRef.getBaseString());
        }
    }

    @Test
    public void testMatchingAlleles() {
        Assert.assertTrue(GATKVariantContextUtils.isAlleleInList(Allele.REF_A, Allele.ALT_T, ATCref, Arrays.asList(Allele.ALT_A, Allele.create("TTC",false))));
        Assert.assertTrue(GATKVariantContextUtils.isAlleleInList(Allele.REF_A, Allele.ALT_T, ATCref, Arrays.asList(Allele.ALT_T, Allele.create("TTC",false))));
    }

    @DataProvider
    Object[][] alleleIndicesData() {
        // {alleles1, alleles2, expected indices}
        return new Object[][] {
                // ref only
                { Arrays.asList(Allele.REF_A), Arrays.asList(Allele.REF_A), Arrays.asList(OptionalInt.of(0))},

                // extra alt
                { Arrays.asList(Allele.REF_A), Arrays.asList(Allele.REF_A, Allele.ALT_C), Arrays.asList(OptionalInt.of(0))},

                // same biallelic
                { Arrays.asList(Allele.REF_A, Allele.ALT_C), Arrays.asList(Allele.REF_A, Allele.ALT_C), Arrays.asList(OptionalInt.of(0), OptionalInt.of(1))},

                // missing alt
                { Arrays.asList(Allele.REF_A, Allele.ALT_C), Arrays.asList(Allele.REF_A), Arrays.asList(OptionalInt.of(0), OptionalInt.empty())},

                // different alt
                { Arrays.asList(Allele.REF_A, Allele.ALT_C), Arrays.asList(Allele.REF_A, Allele.ALT_G), Arrays.asList(OptionalInt.of(0), OptionalInt.empty())},

                // triallelic, order switched
                { Arrays.asList(Allele.REF_A, Allele.ALT_C, Allele.ALT_G), Arrays.asList(Allele.REF_A, Allele.ALT_G, Allele.ALT_C), Arrays.asList(OptionalInt.of(0), OptionalInt.of(2), OptionalInt.of(1))},

                // longer representation
                { Arrays.asList(Allele.REF_A, Allele.ALT_G), Arrays.asList(ATref, GT), Arrays.asList(OptionalInt.of(0), OptionalInt.of(1))},

                // shorter representation
                { Arrays.asList(ATref, GT), Arrays.asList(Allele.REF_A, Allele.ALT_G), Arrays.asList(OptionalInt.of(0), OptionalInt.of(1))}
        };
    }

    @Test(dataProvider = "alleleIndicesData")
    public void testAlleleIndices(final List<Allele> alleles1, final List<Allele> alleles2, final List<OptionalInt> expected ) {
        final List<OptionalInt> result = GATKVariantContextUtils.alleleIndices(alleles1, alleles2);
        Assert.assertEquals(result, expected);
    }

    // --------------------------------------------------------------------------------
    //
    // Test allele merging
    //
    // --------------------------------------------------------------------------------

    private class MergeAllelesTest extends TestDataProvider {
        List<List<Allele>> inputs;
        List<Allele> expected;

        @SafeVarargs
        @SuppressWarnings("varargs")
        private MergeAllelesTest(List<Allele>... arg) {
            super(MergeAllelesTest.class);
            LinkedList<List<Allele>> all = new LinkedList<>(Arrays.asList(arg));
            expected = all.pollLast();
            inputs = all;
        }

        public String toString() {
            return String.format("MergeAllelesTest input=%s expected=%s", inputs, expected);
        }
    }
    @DataProvider(name = "mergeAlleles")
    public Object[][] mergeAllelesData() {
        // first, do no harm
        new MergeAllelesTest(Collections.singletonList(Aref), Collections.singletonList(Aref));

        new MergeAllelesTest(Collections.singletonList(Aref), Collections.singletonList(Aref), Collections.singletonList(Aref));

        new MergeAllelesTest(Collections.singletonList(Aref), Arrays.asList(Aref, T), Arrays.asList(Aref, T));

        new MergeAllelesTest(Arrays.asList(Aref, C), Arrays.asList(Aref, T), Arrays.asList(Aref, C, T));

        new MergeAllelesTest(Arrays.asList(Aref, T),
                Arrays.asList(Aref, C),
                Arrays.asList(Aref, T, C)); // in order of appearence

        new MergeAllelesTest(Arrays.asList(Aref, C, T),
                Arrays.asList(Aref, C),
                Arrays.asList(Aref, C, T));

        new MergeAllelesTest(Arrays.asList(Aref, C, T), Arrays.asList(Aref, C, T));

        new MergeAllelesTest(Arrays.asList(Aref, T, C), Arrays.asList(Aref, T, C));

        new MergeAllelesTest(Arrays.asList(Aref, T, C),
                Arrays.asList(Aref, C),
                Arrays.asList(Aref, T, C)); // in order of appearence

        new MergeAllelesTest(Arrays.asList(Aref),
                Arrays.asList(Aref, ATC),
                Arrays.asList(Aref, ATC));

        new MergeAllelesTest(Arrays.asList(Aref),
                Arrays.asList(Aref, ATC, ATCATC),
                Arrays.asList(Aref, ATC, ATCATC));

        // alleles in the order we see them
        new MergeAllelesTest(Arrays.asList(Aref, ATCATC),
                Arrays.asList(Aref, ATC, ATCATC),
                Arrays.asList(Aref, ATCATC, ATC));

        // same
        new MergeAllelesTest(Arrays.asList(Aref, ATC),
                Arrays.asList(Aref, ATCATC),
                Arrays.asList(Aref, ATC, ATCATC));

        new MergeAllelesTest(Arrays.asList(ATref, ATC, Anoref, G),
                Arrays.asList(Aref, ATCATC, G),
                Arrays.asList(ATref, ATC, Anoref, G, ATCATCT, GT));

        return MergeAllelesTest.getTests(MergeAllelesTest.class);
    }

    @Test( dataProvider = "mergeAlleles")
    public void testMergeAlleles(MergeAllelesTest cfg) {
        final List<VariantContext> inputs = new ArrayList<>();

        int i = 0;
        for ( final List<Allele> alleles : cfg.inputs ) {
            final String name = "vcf" + ++i;
            inputs.add(makeVC(name, alleles));
        }

        final List<String> priority = vcs2priority(inputs);

        final VariantContext merged = GATKVariantContextUtils.simpleMerge(
                inputs, priority,
                GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
                GATKVariantContextUtils.GenotypeMergeType.PRIORITIZE, false);

        Assert.assertEquals(merged.getAlleles().size(),cfg.expected.size());
        Assert.assertEquals(new LinkedHashSet<>(merged.getAlleles()), new LinkedHashSet<>(cfg.expected));   //HACK this is a hack to get around a bug in the htsjdk.  The method returns a list with an unspecified order.
    }

    // --------------------------------------------------------------------------------
    //
    // Test rsID merging
    //
    // --------------------------------------------------------------------------------

    private class SimpleMergeRSIDTest extends TestDataProvider {
        List<String> inputs;
        String expected;

        private SimpleMergeRSIDTest(String... arg) {
            super(SimpleMergeRSIDTest.class);
            LinkedList<String> allStrings = new LinkedList<>(Arrays.asList(arg));
            expected = allStrings.pollLast();
            inputs = allStrings;
        }

        public String toString() {
            return String.format("SimpleMergeRSIDTest vc=%s expected=%s", inputs, expected);
        }
    }

    @DataProvider(name = "simplemergersiddata")
    public Object[][] createSimpleMergeRSIDData() {
        new SimpleMergeRSIDTest(".", ".");
        new SimpleMergeRSIDTest(".", ".", ".");
        new SimpleMergeRSIDTest("rs1", "rs1");
        new SimpleMergeRSIDTest("rs1", "rs1", "rs1");
        new SimpleMergeRSIDTest(".", "rs1", "rs1");
        new SimpleMergeRSIDTest("rs1", ".", "rs1");
        new SimpleMergeRSIDTest("rs1", "rs2", "rs1,rs2");
        new SimpleMergeRSIDTest("rs1", "rs2", "rs1", "rs1,rs2"); // duplicates
        new SimpleMergeRSIDTest("rs2", "rs1", "rs2,rs1");
        new SimpleMergeRSIDTest("rs2", "rs1", ".", "rs2,rs1");
        new SimpleMergeRSIDTest("rs2", ".", "rs1", "rs2,rs1");
        new SimpleMergeRSIDTest("rs1", ".", ".", "rs1");
        new SimpleMergeRSIDTest("rs1", "rs2", "rs3", "rs1,rs2,rs3");

        return SimpleMergeRSIDTest.getTests(SimpleMergeRSIDTest.class);
    }

    @Test(dataProvider = "simplemergersiddata")
    public void testRSIDMerge(SimpleMergeRSIDTest cfg) {
        VariantContext snpVC1 = makeVC("snpvc1", Arrays.asList(Aref, T));
        final List<VariantContext> inputs = cfg.inputs.stream()
                .map(id -> new VariantContextBuilder(snpVC1).id(id).make())
                .collect(Collectors.toList());

        final VariantContext merged = GATKVariantContextUtils.simpleMerge(
                inputs,null,
                GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
                GATKVariantContextUtils.GenotypeMergeType.UNSORTED, false);
        Assert.assertEquals(merged.getID(), cfg.expected);
    }

    // --------------------------------------------------------------------------------
    //
    // Test filtered merging
    //
    // --------------------------------------------------------------------------------

    private class MergeFilteredTest extends TestDataProvider {
        List<VariantContext> inputs;
        VariantContext expected;
        String setExpected;
        GATKVariantContextUtils.FilteredRecordMergeType type;


        private MergeFilteredTest(String name, VariantContext input1, VariantContext input2, VariantContext expected, String setExpected) {
            this(name, input1, input2, expected, GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED, setExpected);
        }

        private MergeFilteredTest(String name, VariantContext input1, VariantContext input2, VariantContext expected, GATKVariantContextUtils.FilteredRecordMergeType type, String setExpected) {
            super(MergeFilteredTest.class, name);
            LinkedList<VariantContext> all = new LinkedList<>(Arrays.asList(input1, input2));
            this.expected = expected;
            this.type = type;
            inputs = all;
            this.setExpected = setExpected;
        }

        public String toString() {
            return String.format("%s input=%s expected=%s", super.toString(), inputs, expected);
        }
    }

    @DataProvider(name = "mergeFiltered")
    public Object[][] mergeFilteredData() {
        new MergeFilteredTest("AllPass",
                makeVC("1", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS),
                makeVC("2", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS),
                makeVC("3", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS),
                GATKVariantContextUtils.MERGE_INTERSECTION);

        new MergeFilteredTest("noFilters",
                makeVC("1", Arrays.asList(Aref, T), "."),
                makeVC("2", Arrays.asList(Aref, T), "."),
                makeVC("3", Arrays.asList(Aref, T), "."),
                GATKVariantContextUtils.MERGE_INTERSECTION);

        new MergeFilteredTest("oneFiltered",
                makeVC("1", Arrays.asList(Aref, T), "."),
                makeVC("2", Arrays.asList(Aref, T), "FAIL"),
                makeVC("3", Arrays.asList(Aref, T), "."),
                String.format("1-%s2", GATKVariantContextUtils.MERGE_FILTER_PREFIX));

        new MergeFilteredTest("onePassOneFail",
                makeVC("1", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS),
                makeVC("2", Arrays.asList(Aref, T), "FAIL"),
                makeVC("3", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS),
                String.format("1-%s2", GATKVariantContextUtils.MERGE_FILTER_PREFIX));

        new MergeFilteredTest("AllFiltered",
                makeVC("1", Arrays.asList(Aref, T), "FAIL"),
                makeVC("2", Arrays.asList(Aref, T), "FAIL"),
                makeVC("3", Arrays.asList(Aref, T), "FAIL"),
                GATKVariantContextUtils.MERGE_FILTER_IN_ALL);

        // test ALL vs. ANY
        new MergeFilteredTest("FailOneUnfiltered",
                makeVC("1", Arrays.asList(Aref, T), "FAIL"),
                makeVC("2", Arrays.asList(Aref, T), "."),
                makeVC("3", Arrays.asList(Aref, T), "."),
                GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
                String.format("%s1-2", GATKVariantContextUtils.MERGE_FILTER_PREFIX));

        new MergeFilteredTest("OneFailAllUnfilteredArg",
                makeVC("1", Arrays.asList(Aref, T), "FAIL"),
                makeVC("2", Arrays.asList(Aref, T), "."),
                makeVC("3", Arrays.asList(Aref, T), "FAIL"),
                GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ALL_UNFILTERED,
                String.format("%s1-2", GATKVariantContextUtils.MERGE_FILTER_PREFIX));

        // test excluding allele in filtered record
        new MergeFilteredTest("DontIncludeAlleleOfFilteredRecords",
                makeVC("1", Arrays.asList(Aref, T), "."),
                makeVC("2", Arrays.asList(Aref, T), "FAIL"),
                makeVC("3", Arrays.asList(Aref, T), "."),
                String.format("1-%s2", GATKVariantContextUtils.MERGE_FILTER_PREFIX));

        // promotion of site from unfiltered to PASSES
        new MergeFilteredTest("UnfilteredPlusPassIsPass",
                makeVC("1", Arrays.asList(Aref, T), "."),
                makeVC("2", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS),
                makeVC("3", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS),
                GATKVariantContextUtils.MERGE_INTERSECTION);

        new MergeFilteredTest("RefInAll",
                makeVC("1", Arrays.asList(Aref), VariantContext.PASSES_FILTERS),
                makeVC("2", Arrays.asList(Aref), VariantContext.PASSES_FILTERS),
                makeVC("3", Arrays.asList(Aref), VariantContext.PASSES_FILTERS),
                GATKVariantContextUtils.MERGE_REF_IN_ALL);

        new MergeFilteredTest("RefInOne",
                makeVC("1", Arrays.asList(Aref), VariantContext.PASSES_FILTERS),
                makeVC("2", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS),
                makeVC("3", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS),
                "2");

        return MergeFilteredTest.getTests(MergeFilteredTest.class);
    }

    @Test(dataProvider = "mergeFiltered")
    public void testMergeFiltered(MergeFilteredTest cfg) {
        final List<String> priority = vcs2priority(cfg.inputs);
        final VariantContext merged = GATKVariantContextUtils.simpleMerge(
                cfg.inputs, priority, cfg.type, GATKVariantContextUtils.GenotypeMergeType.PRIORITIZE, false);

        // test alleles are equal
        Assert.assertEquals(merged.getAlleles(), cfg.expected.getAlleles());

        // test filter field
        Assert.assertEquals(merged.getFilters(), cfg.expected.getFilters());
    }

    // --------------------------------------------------------------------------------
    //
    // Test genotype merging
    //
    // --------------------------------------------------------------------------------

    private class MergeGenotypesTest extends TestDataProvider {
        List<VariantContext> inputs;
        VariantContext expected;
        List<String> priority;

        private MergeGenotypesTest(String name, String priority, VariantContext... arg) {
            super(MergeGenotypesTest.class, name);
            LinkedList<VariantContext> all = new LinkedList<>(Arrays.asList(arg));
            this.expected = all.pollLast();
            inputs = all;
            this.priority = Arrays.asList(priority.split(","));
        }

        public String toString() {
            return String.format("%s input=%s expected=%s", super.toString(), inputs, expected);
        }
    }

    @DataProvider(name = "mergeGenotypes")
    public Object[][] mergeGenotypesData() {
        new MergeGenotypesTest("TakeGenotypeByPriority-1,2", "1,2",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, -1)),
                makeVC("2", Arrays.asList(Aref, T), makeG("s1", Aref, T, -2)),
                makeVC("3", Arrays.asList(Aref, T), makeG("s1", Aref, T, -1)));

        new MergeGenotypesTest("TakeGenotypeByPriority-1,2-nocall", "1,2",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Allele.NO_CALL, Allele.NO_CALL, -1)),
                makeVC("2", Arrays.asList(Aref, T), makeG("s1", Aref, T, -2)),
                makeVC("3", Arrays.asList(Aref, T), makeG("s1", Allele.NO_CALL, Allele.NO_CALL, -1)));

        new MergeGenotypesTest("TakeGenotypeByPriority-2,1", "2,1",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, -1)),
                makeVC("2", Arrays.asList(Aref, T), makeG("s1", Aref, T, -2)),
                makeVC("3", Arrays.asList(Aref, T), makeG("s1", Aref, T, -2)));

        new MergeGenotypesTest("NonOverlappingGenotypes", "1,2",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, -1)),
                makeVC("2", Arrays.asList(Aref, T), makeG("s2", Aref, T, -2)),
                makeVC("3", Arrays.asList(Aref, T), makeG("s1", Aref, T, -1), makeG("s2", Aref, T, -2)));

        new MergeGenotypesTest("PreserveNoCall", "1,2",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Allele.NO_CALL, Allele.NO_CALL, -1)),
                makeVC("2", Arrays.asList(Aref, T), makeG("s2", Aref, T, -2)),
                makeVC("3", Arrays.asList(Aref, T), makeG("s1", Allele.NO_CALL, Allele.NO_CALL, -1), makeG("s2", Aref, T, -2)));

        new MergeGenotypesTest("PerserveAlleles", "1,2",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, -1)),
                makeVC("2", Arrays.asList(Aref, C), makeG("s2", Aref, C, -2)),
                makeVC("3", Arrays.asList(Aref, T, C), makeG("s1", Aref, T, -1), makeG("s2", Aref, C, -2)));

        new MergeGenotypesTest("TakeGenotypePartialOverlap-1,2", "1,2",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, -1)),
                makeVC("2", Arrays.asList(Aref, T), makeG("s1", Aref, T, -2), makeG("s3", Aref, T, -3)),
                makeVC("3", Arrays.asList(Aref, T), makeG("s1", Aref, T, -1), makeG("s3", Aref, T, -3)));

        new MergeGenotypesTest("TakeGenotypePartialOverlap-2,1", "2,1",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, -1)),
                makeVC("2", Arrays.asList(Aref, T), makeG("s1", Aref, T, -2), makeG("s3", Aref, T, -3)),
                makeVC("3", Arrays.asList(Aref, T), makeG("s1", Aref, T, -2), makeG("s3", Aref, T, -3)));

        //
        // merging genothpes with PLs
        //

        // first, do no harm
        new MergeGenotypesTest("OrderedPLs", "1",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, -1, 1, 2, 3)),
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, -1, 1, 2, 3)));

        // first, do no harm
        new MergeGenotypesTest("OrderedPLs-3Alleles", "1",
                makeVC("1", Arrays.asList(Aref, C, T), makeG("s1", Aref, T, -1, 1, 2, 3, 4, 5, 6)),
                makeVC("1", Arrays.asList(Aref, C, T), makeG("s1", Aref, T, -1, 1, 2, 3, 4, 5, 6)));

        // first, do no harm
        new MergeGenotypesTest("OrderedPLs-3Alleles-2", "1",
                makeVC("1", Arrays.asList(Aref, T, C), makeG("s1", Aref, T, -1, 1, 2, 3, 4, 5, 6)),
                makeVC("1", Arrays.asList(Aref, T, C), makeG("s1", Aref, T, -1, 1, 2, 3, 4, 5, 6)));

        // first, do no harm
        new MergeGenotypesTest("OrderedPLs-3Alleles-2", "1",
                makeVC("1", Arrays.asList(Aref, T, C), makeG("s1", Aref, T, -1, 1, 2, 3, 4, 5, 6)),
                makeVC("1", Arrays.asList(Aref, T, C), makeG("s2", Aref, C, -1, 1, 2, 3, 4, 5, 6)),
                makeVC("1", Arrays.asList(Aref, T, C), makeG("s1", Aref, T, -1, 1, 2, 3, 4, 5, 6), makeG("s2", Aref, C, -1, 1, 2, 3, 4, 5, 6)));

        new MergeGenotypesTest("TakeGenotypePartialOverlapWithPLs-2,1", "2,1",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, -1,5,0,3)),
                makeVC("2", Arrays.asList(Aref, T), makeG("s1", Aref, T, -2,4,0,2), makeG("s3", Aref, T, -3,3,0,2)),
                makeVC("3", Arrays.asList(Aref, T), makeG("s1", Aref, T, -2,4,0,2), makeG("s3", Aref, T, -3,3,0,2)));

        new MergeGenotypesTest("TakeGenotypePartialOverlapWithPLs-1,2", "1,2",
                makeVC("1", Arrays.asList(Aref,ATC), makeG("s1", Aref, ATC, -1,5,0,3)),
                makeVC("2", Arrays.asList(Aref, T), makeG("s1", Aref, T, -2,4,0,2), makeG("s3", Aref, T, -3,3,0,2)),
                // no likelihoods on result since type changes to mixed multiallelic
                makeVC("3", Arrays.asList(Aref, ATC, T), makeG("s1", Aref, ATC, -1), makeG("s3", Aref, T, -3)));

        new MergeGenotypesTest("MultipleSamplePLsDifferentOrder", "1,2",
                makeVC("1", Arrays.asList(Aref, C, T), makeG("s1", Aref, C, -1, 1, 2, 3, 4, 5, 6)),
                makeVC("2", Arrays.asList(Aref, T, C), makeG("s2", Aref, T, -2, 6, 5, 4, 3, 2, 1)),
                // no likelihoods on result since type changes to mixed multiallelic
                makeVC("3", Arrays.asList(Aref, C, T), makeG("s1", Aref, C, -1), makeG("s2", Aref, T, -2)));

        return MergeGenotypesTest.getTests(MergeGenotypesTest.class);
    }

    @Test(dataProvider = "mergeGenotypes")
    public void testMergeGenotypes(MergeGenotypesTest cfg) {
        final VariantContext merged = GATKVariantContextUtils.simpleMerge(
                cfg.inputs, cfg.priority, GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
                GATKVariantContextUtils.GenotypeMergeType.PRIORITIZE, false);

        // test alleles are equal
        Assert.assertEquals(merged.getAlleles(), cfg.expected.getAlleles());

        // test genotypes
        assertGenotypesAreMostlyEqual(merged.getGenotypes(), cfg.expected.getGenotypes());
    }

    // necessary to not overload equals for genotypes
    private void assertGenotypesAreMostlyEqual(GenotypesContext actual, GenotypesContext expected) {
        if (actual == expected) {
            return;
        }

        if (actual == null || expected == null) {
            Assert.fail("Maps not equal: expected: " + expected + " and actual: " + actual);
        }

        if (actual.size() != expected.size()) {
            Assert.fail("Maps do not have the same size:" + actual.size() + " != " + expected.size());
        }

        for (Genotype value : actual) {
            Genotype expectedValue = expected.get(value.getSampleName());

            Assert.assertEquals(value.getAlleles(), expectedValue.getAlleles(), "Alleles in Genotype aren't equal");
            Assert.assertEquals(value.getGQ(), expectedValue.getGQ(), "GQ values aren't equal");
            Assert.assertEquals(value.hasLikelihoods(), expectedValue.hasLikelihoods(), "Either both have likelihoods or both not");
            if ( value.hasLikelihoods() )
                Assert.assertEquals(value.getLikelihoods().getAsVector(), expectedValue.getLikelihoods().getAsVector(), "Genotype likelihoods aren't equal");
        }
    }

    @Test
    public void testMergeGenotypesUniquify() {
        final VariantContext vc1 = makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, -1));
        final VariantContext vc2 = makeVC("2", Arrays.asList(Aref, T), makeG("s1", Aref, T, -2));

        final VariantContext merged = GATKVariantContextUtils.simpleMerge(
                Arrays.asList(vc1, vc2), null, GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
                GATKVariantContextUtils.GenotypeMergeType.UNIQUIFY, false);

        // test genotypes
        Assert.assertEquals(merged.getSampleNames(), new LinkedHashSet<>(Arrays.asList("s1.1", "s1.2")));
    }

// TODO: remove after testing
//    @Test(expectedExceptions = IllegalStateException.class)
//    public void testMergeGenotypesRequireUnique() {
//        final VariantContext vc1 = makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, -1));
//        final VariantContext vc2 = makeVC("2", Arrays.asList(Aref, T), makeG("s1", Aref, T, -2));
//
//        final VariantContext merged = VariantContextUtils.simpleMerge(
//                Arrays.asList(vc1, vc2), null, GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
//                GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE, false, false, "set", false, false);
//    }

    // --------------------------------------------------------------------------------
    //
    // Misc. tests
    //
    // --------------------------------------------------------------------------------
    
    private static List<String> vcs2priority(final Collection<VariantContext> vcs) {
        return vcs.stream()
                .map(VariantContext::getSource)
                .collect(Collectors.toList());
    }

    // --------------------------------------------------------------------------------
    //
    // basic allele clipping test
    //
    // --------------------------------------------------------------------------------

    private class ReverseClippingPositionTestProvider extends TestDataProvider {
        final String ref;
        final List<Allele> alleles = new ArrayList<>();
        final int expectedClip;

        private ReverseClippingPositionTestProvider(final int expectedClip, final String ref, final String... alleles) {
            super(ReverseClippingPositionTestProvider.class);
            this.ref = ref;
            for ( final String allele : alleles )
                this.alleles.add(Allele.create(allele));
            this.expectedClip = expectedClip;
        }

        @Override
        public String toString() {
            return String.format("ref=%s allele=%s reverse clip %d", ref, alleles, expectedClip);
        }
    }

    // --------------------------------------------------------------------------------
    //
    // test splitting into bi-allelics
    //
    // --------------------------------------------------------------------------------

    @DataProvider(name = "SplitBiallelics")
    public Object[][] makeSplitBiallelics() throws CloneNotSupportedException {
        List<Object[]> tests = new ArrayList<>();

        final VariantContextBuilder root = new VariantContextBuilder("x", "20", 10, 10, Arrays.asList(Aref, C));

        // biallelic -> biallelic
        tests.add(new Object[]{root.make(), Collections.singletonList(root.make())});

        // monos -> empty
        root.alleles(Collections.singletonList(Aref));
        tests.add(new Object[]{root.make(), Collections.emptyList()});

        root.alleles(Arrays.asList(Aref, C, T));
        tests.add(new Object[]{root.make(),
                Arrays.asList(
                        root.alleles(Arrays.asList(Aref, C)).make(),
                        root.alleles(Arrays.asList(Aref, T)).make())});

        root.alleles(Arrays.asList(Aref, C, T, G));
        tests.add(new Object[]{root.make(),
                Arrays.asList(
                        root.alleles(Arrays.asList(Aref, C)).make(),
                        root.alleles(Arrays.asList(Aref, T)).make(),
                        root.alleles(Arrays.asList(Aref, G)).make())});

        final Allele C      = Allele.create("C");
        final Allele CA      = Allele.create("CA");
        final Allele CAA     = Allele.create("CAA");
        final Allele CAAAA   = Allele.create("CAAAA");
        final Allele CAAAAA  = Allele.create("CAAAAA");
        final Allele Cref      = Allele.create("C", true);
        final Allele CAref     = Allele.create("CA", true);
        final Allele CAAref    = Allele.create("CAA", true);
        final Allele CAAAref   = Allele.create("CAAA", true);

        root.alleles(Arrays.asList(Cref, CA, CAA));
        tests.add(new Object[]{root.make(),
                Arrays.asList(
                        root.alleles(Arrays.asList(Cref, CA)).make(),
                        root.alleles(Arrays.asList(Cref, CAA)).make())});

        root.alleles(Arrays.asList(CAAref, C, CA)).stop(12);
        tests.add(new Object[]{root.make(),
                Arrays.asList(
                        root.alleles(Arrays.asList(CAAref, C)).make(),
                        root.alleles(Arrays.asList(CAref, C)).stop(11).make())});

        root.alleles(Arrays.asList(CAAAref, C, CA, CAA)).stop(13);
        tests.add(new Object[]{root.make(),
                Arrays.asList(
                        root.alleles(Arrays.asList(CAAAref, C)).make(),
                        root.alleles(Arrays.asList(CAAref, C)).stop(12).make(),
                        root.alleles(Arrays.asList(CAref, C)).stop(11).make())});

        root.alleles(Arrays.asList(CAAAref, CAAAAA, CAAAA, CAA, C)).stop(13);
        tests.add(new Object[]{root.make(),
                Arrays.asList(
                        root.alleles(Arrays.asList(Cref, CAA)).stop(10).make(),
                        root.alleles(Arrays.asList(Cref, CA)).stop(10).make(),
                        root.alleles(Arrays.asList(CAref, C)).stop(11).make(),
                        root.alleles(Arrays.asList(CAAAref, C)).stop(13).make())});

        final Allele threeCopies = Allele.create("GTTTTATTTTATTTTA", true);
        final Allele twoCopies = Allele.create("GTTTTATTTTA", true);
        final Allele zeroCopies = Allele.create("G", false);
        final Allele oneCopies = Allele.create("GTTTTA", false);
        tests.add(new Object[]{root.alleles(Arrays.asList(threeCopies, zeroCopies, oneCopies)).stop(25).make(),
                Arrays.asList(
                        root.alleles(Arrays.asList(threeCopies, zeroCopies)).stop(25).make(),
                        root.alleles(Arrays.asList(twoCopies, zeroCopies)).stop(20).make())});

        return tests.toArray(new Object[][]{});
    }


    @DataProvider(name = "GQLog10PosteriorsTest")
    public Object[][] makeGQLog10PosteriorsTest() {
        List<Object[]> tests = new ArrayList<>();

        // testing the 3 allele case
        tests.add(new Object[]{0, new double[]{-1.0, -2.0, -2.2}, -1.787});
        tests.add(new Object[]{1, new double[]{-1.0, -2.0, -2.2}, -0.973});
        tests.add(new Object[]{2, new double[]{-1.0, -2.0, -2.2}, -0.958});

        // testing in the 3 allele case where the choice between two genotypes is ambiguous
        tests.add(new Object[]{0, new double[]{0.0, 0.0, -0.2}, 0});
        tests.add(new Object[]{1, new double[]{0.0, 0.0, -0.2}, 0});
        tests.add(new Object[]{2, new double[]{0.0, 0.0, -0.2}, 0});

        // testing in the 4+ allele case where the choice between two genotypes is ambiguous (if not careful this might have resulted in a negative GQ)
        tests.add(new Object[]{0, new double[]{0.0, 0.0, -0.2, -0.2, -0.2, 0.0}, 0});
        tests.add(new Object[]{1, new double[]{0.0, 0.0, -0.2, -0.2, -0.2, 0.0}, 0});
        tests.add(new Object[]{2, new double[]{0.0, 0.0, -0.2, -0.2, -0.2, 0.0}, 0});
        tests.add(new Object[]{3, new double[]{0.0, 0.0, -0.2, -0.2, -0.2, 0.0}, 0});
        tests.add(new Object[]{4, new double[]{0.0, 0.0, -0.2, -0.2, -0.2, 0.0}, 0});
        tests.add(new Object[]{5, new double[]{0.0, 0.0, -0.2, -0.2, -0.2, 0.0}, 0});


        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "GQLog10PosteriorsTest")
    public void testGetGQLog10FromPosteriors(final int bestGenotypeIndex, final double[] genotypeArray, final double expectedResult)  {
        final double actualResult = GATKVariantContextUtils.getGQLog10FromPosteriors(bestGenotypeIndex, genotypeArray);
        assertEqualsDoubleSmart(actualResult, expectedResult);
    }

    // --------------------------------------------------------------------------------
    //
    // Test repeats
    //
    // --------------------------------------------------------------------------------

    private class RepeatDetectorTest extends TestDataProvider {
        String ref;
        boolean isTrueRepeat;
        VariantContext vc;

        private RepeatDetectorTest(boolean isTrueRepeat, String ref, String refAlleleString, String ... altAlleleStrings) {
            super(RepeatDetectorTest.class);
            this.isTrueRepeat = isTrueRepeat;
            this.ref = ref;

            List<Allele> alleles = new LinkedList<>();
            final Allele refAllele = Allele.create(refAlleleString, true);
            alleles.add(refAllele);
            for ( final String altString: altAlleleStrings) {
                final Allele alt = Allele.create(altString, false);
                alleles.add(alt);
            }

            VariantContextBuilder builder = new VariantContextBuilder("test", "chr1", 1, refAllele.length(), alleles);
            this.vc = builder.make();
        }

        public String toString() {
            return String.format("%s refBases=%s trueRepeat=%b vc=%s", super.toString(), ref, isTrueRepeat, vc);
        }
    }

    @DataProvider(name = "RepeatDetectorTest")
    public Object[][] makeRepeatDetectorTest() {
        new RepeatDetectorTest(true,  "NAAC", "N", "NA");
        new RepeatDetectorTest(true,  "NAAC", "NA", "N");
        new RepeatDetectorTest(false, "NAAC", "NAA", "N");
        new RepeatDetectorTest(false, "NAAC", "N", "NC");
        new RepeatDetectorTest(false, "AAC", "A", "C");

        // running out of ref bases => false
        new RepeatDetectorTest(false, "NAAC", "N", "NCAGTA");

        // complex repeats
        new RepeatDetectorTest(true,  "NATATATC", "N", "NAT");
        new RepeatDetectorTest(true,  "NATATATC", "N", "NATA");
        new RepeatDetectorTest(true,  "NATATATC", "N", "NATAT");
        new RepeatDetectorTest(true,  "NATATATC", "NAT", "N");
        new RepeatDetectorTest(false, "NATATATC", "NATA", "N");
        new RepeatDetectorTest(false, "NATATATC", "NATAT", "N");

        // multi-allelic
        new RepeatDetectorTest(true,  "NATATATC", "N", "NAT", "NATAT");
        new RepeatDetectorTest(true,  "NATATATC", "N", "NAT", "NATA");
        new RepeatDetectorTest(true,  "NATATATC", "NAT", "N", "NATAT");
        new RepeatDetectorTest(true,  "NATATATC", "NAT", "N", "NATA"); // two As
        new RepeatDetectorTest(false, "NATATATC", "NAT", "N", "NATC"); // false
        new RepeatDetectorTest(false, "NATATATC", "NAT", "N", "NCC"); // false
        new RepeatDetectorTest(false, "NATATATC", "NAT", "NATAT", "NCC"); // false

        return RepeatDetectorTest.getTests(RepeatDetectorTest.class);
    }

    @Test( dataProvider = "RepeatDetectorTest")
    public void testRepeatDetectorTest(RepeatDetectorTest cfg) {

        // test alleles are equal
        Assert.assertEquals(GATKVariantContextUtils.isTandemRepeat(cfg.vc, cfg.ref.getBytes()), cfg.isTrueRepeat);
    }

    @Test
    public void testFindNumberOfRepetitions() throws Exception {
        /*
         *    GATAT has 0 leading repeats of AT but 2 trailing repeats of AT
         *    ATATG has 1 leading repeat of A but 2 leading repeats of AT
         *    CCCCCCCC has 2 leading and 2 trailing repeats of CCC
         */
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("AT".getBytes(), "GATAT".getBytes(), false), 2);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("AT".getBytes(), "GATAT".getBytes(), true), 0);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("A".getBytes(), "ATATG".getBytes(), true), 1);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("AT".getBytes(), "ATATG".getBytes(), true), 2);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("CCC".getBytes(), "CCCCCCCC".getBytes(), true),2);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("CCC".getBytes(), "CCCCCCCC".getBytes(), false),2);

        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("ATG".getBytes(), "ATGATGATGATG".getBytes(), true),4);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("G".getBytes(), "ATGATGATGATG".getBytes(), true),0);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("T".getBytes(), "T".getBytes(), true),1);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("AT".getBytes(), "ATGATGATCATG".getBytes(), true),1);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("CCCCCCCC".getBytes(), "CCC".getBytes(), true),0);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("AT".getBytes(), "AT".getBytes(), true), 1);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("AT".getBytes(), "".getBytes(), true), 0); //empty test string

        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("ATG".getBytes(), "ATGATGATGATG".getBytes(), false),4);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("G".getBytes(), "ATGATGATGATG".getBytes(), false),1);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("T".getBytes(), "T".getBytes(), false),1);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("AT".getBytes(), "ATGATGATCATG".getBytes(), false),0);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("CCCCCCCC".getBytes(), "CCC".getBytes(), false),0);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("AT".getBytes(), "AT".getBytes(), false), 1);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("AT".getBytes(), "".getBytes(), false), 0); //empty test string
    }

    @Test
    public void testFindNumberOfRepetitionsFullArray() throws Exception {
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("XXXATG".getBytes(),3 ,3, "ATGATGATGATGYYY".getBytes(), 0, 12, true),4);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("GGGG".getBytes(), 0 ,1 , "GGGGATGATGATGATG".getBytes(), 4, 12, true),0);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("T".getBytes(), 0, 1, "TTTTT".getBytes(), 0, 1, true),1);

        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("AT".getBytes(), 0, 2, "AT".getBytes(), 0, 0, true), 0); //empty test string
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("AT".getBytes(), 0, 2, "AT".getBytes(), 1, 0, true), 0); //empty test string
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("AT".getBytes(), 0, 2, "".getBytes(), 0, 0, true), 0); //empty test string

        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("XXXAT".getBytes(), 3, 2, "XXXGATAT".getBytes(), 4, 4, false), 2);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("AT".getBytes(), 0, 2, "GATAT".getBytes(), 0, 5, false), 2);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("ATG".getBytes(), 0, 3, "ATGATGATGATG".getBytes(), 0, 12, false),4);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("ATG".getBytes(), 0, 3, "ATGATGATGATGATG".getBytes(), 3, 12, false),4);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("G".getBytes(), 0, 1, "ATGATGATGATG".getBytes(), 0, 12, false),1);
        Assert.assertEquals(GATKVariantContextUtils.findNumberOfRepetitions("G".getBytes(), 0, 1, "ATGATGATGATGATG".getBytes(), 0, 12, false),1);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testFindNumberOfRepetitionsNullArg1() throws Exception {
        GATKVariantContextUtils.findNumberOfRepetitions(null, "AT".getBytes(), false);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testFindNumberOfRepetitionsNullArg2() throws Exception {
        GATKVariantContextUtils.findNumberOfRepetitions("AT".getBytes(), null, false);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testFindNumberOfRepetitionsEmptyArg1() throws Exception {
        GATKVariantContextUtils.findNumberOfRepetitions("".getBytes(), "AT".getBytes(), false);
    }

    @Test
    public void testRepeatAllele() {
        Allele nullR = Aref;
        Allele nullA = Allele.create("A", false);
        Allele atc   = Allele.create("AATC", false);
        Allele atcatc   = Allele.create("AATCATC", false);
        Allele ccccR = Allele.create("ACCCC", true);
        Allele cc   = Allele.create("ACC", false);
        Allele cccccc   = Allele.create("ACCCCCC", false);
        Allele gagaR   = Allele.create("AGAGA", true);
        Allele gagagaga   = Allele.create("AGAGAGAGA", false);

        // - / ATC [ref] from 20-22
        String delLoc = "chr1";

        // - [ref] / ATC from 20-20
        String insLoc = "chr1";
        int insLocStart = 20;
        int insLocStop = 20;

        Pair<List<Integer>,byte[]> result;
        byte[] refBytes = "ATCATCATCGGA".getBytes();    // excludes leading match base common to VC's ref and alt alleles

        Assert.assertEquals(GATKVariantContextUtils.findRepeatedSubstring("ATG".getBytes()),3);
        Assert.assertEquals(GATKVariantContextUtils.findRepeatedSubstring("AAA".getBytes()),1);
        Assert.assertEquals(GATKVariantContextUtils.findRepeatedSubstring("CACACAC".getBytes()),7);
        Assert.assertEquals(GATKVariantContextUtils.findRepeatedSubstring("CACACA".getBytes()),2);
        Assert.assertEquals(GATKVariantContextUtils.findRepeatedSubstring("CATGCATG".getBytes()),4);
        Assert.assertEquals(GATKVariantContextUtils.findRepeatedSubstring("AATAATA".getBytes()),7);


        // A*,ATC, context = ATC ATC ATC : (ATC)3 -> (ATC)4
        VariantContext vc = new VariantContextBuilder("foo", insLoc, insLocStart, insLocStop, Arrays.asList(nullR,atc)).make();
        result = GATKVariantContextUtils.getNumTandemRepeatUnits(vc, refBytes);
        Assert.assertEquals(result.getLeft().toArray()[0],3);
        Assert.assertEquals(result.getLeft().toArray()[1],4);
        Assert.assertEquals(result.getRight().length,3);

        // ATC*,A,ATCATC
        vc = new VariantContextBuilder("foo", insLoc, insLocStart, insLocStart+3, Arrays.asList(Allele.create("AATC", true),nullA,atcatc)).make();
        result = GATKVariantContextUtils.getNumTandemRepeatUnits(vc, refBytes);
        Assert.assertEquals(result.getLeft().toArray()[0],3);
        Assert.assertEquals(result.getLeft().toArray()[1],2);
        Assert.assertEquals(result.getLeft().toArray()[2],4);
        Assert.assertEquals(result.getRight().length,3);

        // simple non-tandem deletion: CCCC*, -
        refBytes = "CCCCCCCCATG".getBytes();    // excludes leading match base common to VC's ref and alt alleles
        vc = new VariantContextBuilder("foo", delLoc, 10, 14, Arrays.asList(ccccR,nullA)).make();
        result = GATKVariantContextUtils.getNumTandemRepeatUnits(vc, refBytes);
        Assert.assertEquals(result.getLeft().toArray()[0],8);
        Assert.assertEquals(result.getLeft().toArray()[1],4);
        Assert.assertEquals(result.getRight().length,1);

        // CCCC*,CC,-,CCCCCC, context = CCC: (C)7 -> (C)5,(C)3,(C)9
        refBytes = "CCCCCCCAGAGAGAG".getBytes();    // excludes leading match base common to VC's ref and alt alleles
        vc = new VariantContextBuilder("foo", insLoc, insLocStart, insLocStart+4, Arrays.asList(ccccR,cc, nullA,cccccc)).make();
        result = GATKVariantContextUtils.getNumTandemRepeatUnits(vc, refBytes);
        Assert.assertEquals(result.getLeft().toArray()[0],7);
        Assert.assertEquals(result.getLeft().toArray()[1],5);
        Assert.assertEquals(result.getLeft().toArray()[2],3);
        Assert.assertEquals(result.getLeft().toArray()[3],9);
        Assert.assertEquals(result.getRight().length,1);

        // GAGA*,-,GAGAGAGA
        refBytes = "GAGAGAGAGATTT".getBytes();  // excludes leading match base common to VC's ref and alt alleles
        vc = new VariantContextBuilder("foo", insLoc, insLocStart, insLocStart+4, Arrays.asList(gagaR, nullA,gagagaga)).make();
        result = GATKVariantContextUtils.getNumTandemRepeatUnits(vc, refBytes);
        Assert.assertEquals(result.getLeft().toArray()[0],5);
        Assert.assertEquals(result.getLeft().toArray()[1],3);
        Assert.assertEquals(result.getLeft().toArray()[2],7);
        Assert.assertEquals(result.getRight().length,2);

    }

    // --------------------------------------------------------------------------------
    //
    // test forward clipping
    //
    // --------------------------------------------------------------------------------

    @DataProvider(name = "ClipAlleleTest")
    public Object[][] makeClipAlleleTest() {
        List<Object[]> tests = new ArrayList<>();

        // this functionality can be adapted to provide input data for whatever you might want in your data
        tests.add(new Object[]{Arrays.asList("ACC", "AC", "<NON_REF>"), Arrays.asList("AC", "A", "<NON_REF>"), 0});
        tests.add(new Object[]{Arrays.asList("ACC", "AC", "*"), Arrays.asList("AC", "A", "*"), 0});

        tests.add(new Object[]{Arrays.asList("ACC", "AC"), Arrays.asList("AC", "A"), 0});
        tests.add(new Object[]{Arrays.asList("ACGC", "ACG"), Arrays.asList("GC", "G"), 2});
        tests.add(new Object[]{Arrays.asList("ACGC", "ACGA"), Arrays.asList("C", "A"), 3});
        tests.add(new Object[]{Arrays.asList("ACGC", "AGC"), Arrays.asList("AC", "A"), 0});
        tests.add(new Object[]{Arrays.asList("AT", "AC", "AG"), Arrays.asList("T", "C", "G"), 1});
        tests.add(new Object[]{Arrays.asList("AT", "AC", "ACG"), Arrays.asList("T", "C", "CG"), 1});
        tests.add(new Object[]{Arrays.asList("AC", "ACT", "ACG"), Arrays.asList("C", "CT", "CG"), 1});
        tests.add(new Object[]{Arrays.asList("ACG", "ACGT", "ACGTA"), Arrays.asList("G", "GT", "GTA"), 2});
        tests.add(new Object[]{Arrays.asList("ACG", "ACGT", "ACGCA"), Arrays.asList("G", "GT", "GCA"), 2});

        // trims from left and right
        tests.add(new Object[]{Arrays.asList("ACGTT", "ACCTT"), Arrays.asList("G", "C"), 2});
        tests.add(new Object[]{Arrays.asList("ACGTT", "ACCCTT"), Arrays.asList("G", "CC"), 2});
        tests.add(new Object[]{Arrays.asList("ACGTT", "ACGCTT"), Arrays.asList("G", "GC"), 2});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ClipAlleleTest")
    public void testClipAlleles(final List<String> alleleStrings, final List<String> expected, final int numLeftClipped) {
        final int start = 10;
        final VariantContext unclipped = GATKVariantContextUtils.makeFromAlleles("test", "20", start, alleleStrings);
        final VariantContext clipped = GATKVariantContextUtils.trimAlleles(unclipped, true, true);

        Assert.assertEquals(clipped.getStart(), unclipped.getStart() + numLeftClipped);
        for ( int i = 0; i < unclipped.getAlleles().size(); i++ ) {
            final Allele trimmed = clipped.getAlleles().get(i);
            Assert.assertEquals(trimmed.getDisplayString(), expected.get(i));  //note that getBaseString doesn't work for <NON_REF>
        }
    }

    // --------------------------------------------------------------------------------
    //
    // test allele remapping
    //
    // --------------------------------------------------------------------------------

    @DataProvider(name = "AlleleRemappingData")
    public Object[][] makeAlleleRemappingData() {
        List<Object[]> tests = new ArrayList<>();

        final Allele originalBase1 = Allele.create((byte)'A');
        final Allele originalBase2 = Allele.create((byte)'T');

        for ( final byte base1 : BaseUtils.BASES ) {
            for ( final byte base2 : BaseUtils.BASES ) {
                for ( final int numGenotypes : Arrays.asList(0, 1, 2, 5) ) {
                    Map<Allele, Allele> map = new LinkedHashMap<>(2);
                    map.put(originalBase1, Allele.create(base1));
                    map.put(originalBase2, Allele.create(base2));

                    tests.add(new Object[]{map, numGenotypes});
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "AlleleRemappingData")
    public void testAlleleRemapping(final Map<Allele, Allele> alleleMap, final int numGenotypes) {

        final GATKVariantContextUtils.AlleleMapper alleleMapper = new GATKVariantContextUtils.AlleleMapper(alleleMap);

        @SuppressWarnings("unchecked")
        final GenotypesContext originalGC = createGenotypesContext(numGenotypes, new ArrayList<>(alleleMap.keySet()));

        final GenotypesContext remappedGC = GATKVariantContextUtils.updateGenotypesWithMappedAlleles(originalGC, alleleMapper);

        for ( int i = 0; i < numGenotypes; i++ ) {
            final Genotype originalG = originalGC.get(String.format("%d", i));
            final Genotype remappedG = remappedGC.get(String.format("%d", i));

            Assert.assertEquals(originalG.getAlleles().size(), remappedG.getAlleles().size());
            for ( int j = 0; j < originalG.getAlleles().size(); j++ )
                Assert.assertEquals(remappedG.getAllele(j), alleleMap.get(originalG.getAllele(j)));
        }
    }

    private static GenotypesContext createGenotypesContext(final int numGenotypes, final List<Allele> alleles) {
        Utils.resetRandomGenerator();
        final Random random = Utils.getRandomGenerator();

        final GenotypesContext gc = GenotypesContext.create();
        for ( int i = 0; i < numGenotypes; i++ ) {
            // choose alleles at random
            final List<Allele> myAlleles = new ArrayList<>();
            myAlleles.add(alleles.get(random.nextInt(2)));
            myAlleles.add(alleles.get(random.nextInt(2)));

            final Genotype g = new GenotypeBuilder(String.format("%d", i)).alleles(myAlleles).make();
            gc.add(g);
        }

        return gc;
    }

    // --------------------------------------------------------------------------------
    //
    // Test subsetDiploidAlleles
    //
    // --------------------------------------------------------------------------------

    @DataProvider(name = "SubsetAllelesData")
    public Object[][] makesubsetAllelesData() {
        List<Object[]> tests = new ArrayList<>();

        final List<Allele> AA = Arrays.asList(Aref,Aref);
        final List<Allele> AC = Arrays.asList(Aref,C);
        final List<Allele> CC = Arrays.asList(C,C);
        final List<Allele> AG = Arrays.asList(Aref,G);
        final List<Allele> GG = Arrays.asList(G,G);
        final List<Allele> ACG = Arrays.asList(Aref,C,G);

        final VariantContext vcBase = new VariantContextBuilder("test", "20", 10, 10, AC).make();

        // haploid, one alt allele
        final double[] haploidRefPL = MathUtils.normalizeSumToOne(new double[]{0.9, 0.1});
        final double[] haploidAltPL = MathUtils.normalizeSumToOne(new double[]{0.1, 0.9});
        final double[] haploidUninformative = new double[]{0, 0};

        // diploid, one alt allele
        final double[] homRefPL = MathUtils.normalizeSumToOne(new double[]{0.9, 0.09, 0.01});
        final double[] hetPL = MathUtils.normalizeSumToOne(new double[]{0.09, 0.9, 0.01});
        final double[] homVarPL = MathUtils.normalizeSumToOne(new double[]{0.01, 0.09, 0.9});
        final double[] uninformative = new double[]{0, 0, 0};

        final Genotype base = new GenotypeBuilder("NA12878").DP(10).GQ(50).make();

        // the simple case where no selection occurs
        final Genotype aHaploidGT = new GenotypeBuilder(base).alleles(Collections.singletonList(Aref)).AD(new int[]{10,2}).PL(haploidRefPL).GQ(8).make();
        final Genotype cHaploidGT = new GenotypeBuilder(base).alleles(Collections.singletonList(C)).AD(new int[]{10,2}).PL(haploidAltPL).GQ(8).make();
        final Genotype aaGT = new GenotypeBuilder(base).alleles(AA).AD(new int[]{10,2}).PL(homRefPL).GQ(8).make();
        final Genotype acGT = new GenotypeBuilder(base).alleles(AC).AD(new int[]{10,2}).PL(hetPL).GQ(8).make();
        final Genotype ccGT = new GenotypeBuilder(base).alleles(CC).AD(new int[]{10,2}).PL(homVarPL).GQ(8).make();

        // haploid
        tests.add(new Object[]{new VariantContextBuilder(vcBase).genotypes(aHaploidGT).make(), AC, Collections.singletonList(new GenotypeBuilder(aHaploidGT).make())});
        tests.add(new Object[]{new VariantContextBuilder(vcBase).genotypes(cHaploidGT).make(), AC, Collections.singletonList(new GenotypeBuilder(cHaploidGT).make())});
        // diploid
        tests.add(new Object[]{new VariantContextBuilder(vcBase).genotypes(aaGT).make(), AC, Collections.singletonList(new GenotypeBuilder(aaGT).make())});
        tests.add(new Object[]{new VariantContextBuilder(vcBase).genotypes(acGT).make(), AC, Collections.singletonList(new GenotypeBuilder(acGT).make())});
        tests.add(new Object[]{new VariantContextBuilder(vcBase).genotypes(ccGT).make(), AC, Collections.singletonList(new GenotypeBuilder(ccGT).make())});

        // uninformative test cases
        // diploid
        final Genotype uninformativeGT = new GenotypeBuilder(base).alleles(CC).PL(uninformative).GQ(0).make();
        final Genotype emptyGT = new GenotypeBuilder(base).alleles(GATKVariantContextUtils.noCallAlleles(2)).noPL().noGQ().make();
        tests.add(new Object[]{new VariantContextBuilder(vcBase).genotypes(uninformativeGT).make(), AC, Collections.singletonList(emptyGT)});
        // haploid
        final Genotype haploidUninformativeGT = new GenotypeBuilder(base).alleles(Collections.singletonList(C)).PL(haploidUninformative).GQ(0).make();
        final Genotype haplpoidEmptyGT = new GenotypeBuilder(base).alleles(GATKVariantContextUtils.noCallAlleles(1)).noPL().noGQ().make();
        tests.add(new Object[]{new VariantContextBuilder(vcBase).genotypes(haploidUninformativeGT).make(), AC, Collections.singletonList(haplpoidEmptyGT)});

        // subsetting from 3 to 2 alleles
        // diploid PL order is: AA, AC, CC, AG, CG, GG (00, 01, 11, 02, 12, 22)
        final double[] homRef3AllelesPL = new double[]{0, -10, -20, -30, -40, -50};
        final double[] hetRefC3AllelesPL = new double[]{-10, 0, -20, -30, -40, -50};
        final double[] homC3AllelesPL = new double[]{-20, -10, 0, -30, -40, -50};
        final double[] hetRefG3AllelesPL = new double[]{-20, -10, -30, 0, -40, -50};
        final double[] hetCG3AllelesPL = new double[]{-20, -10, -30, -40, 0, -50};
        final double[] homG3AllelesPL = new double[]{-20, -10, -30, -40, -50, 0};
        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(AA).PL(homRef3AllelesPL).make()).make(),
                AC,
                Collections.singletonList(new GenotypeBuilder(base).alleles(AA).PL(new double[]{0, -10, -20}).GQ(100).make())});
        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(AA).PL(hetRefC3AllelesPL).make()).make(),
                AC,
                Collections.singletonList(new GenotypeBuilder(base).alleles(AC).PL(new double[]{-10, 0, -20}).GQ(100).make())});
        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(AA).PL(homC3AllelesPL).make()).make(),
                AC,
                Collections.singletonList(new GenotypeBuilder(base).alleles(CC).PL(new double[]{-20, -10, 0}).GQ(100).make())});
        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(AA).PL(hetRefG3AllelesPL).make()).make(),
                AG,
                Collections.singletonList(new GenotypeBuilder(base).alleles(AG).PL(new double[]{-20, 0, -50}).GQ(200).make())});
        // wow, scary -- bad output but discussed with Eric and we think this is the only thing that can be done
        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(AA).PL(hetCG3AllelesPL).make()).make(),
                AG,
                Collections.singletonList(new GenotypeBuilder(base).alleles(AA).PL(new double[]{0, -20, -30}).GQ(200).make())});

        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(AA).PL(homG3AllelesPL).make()).make(),
                AG,
                Collections.singletonList(new GenotypeBuilder(base).alleles(GG).PL(new double[]{-20, -40, 0}).GQ(200).make())});

        // haploid PL order is: A, C, G (0, 1, 2)
        final double[] haploidRef3AllelesPL = new double[]{0, -10, -20};
        final double[] haploidAltC3AllelesPL = new double[]{-10, 0, -20};
        final double[] haploidAltG3AllelesPL = new double[]{-20, -10, 0};
        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(Collections.singletonList(Aref)).PL(haploidRef3AllelesPL).make()).make(),
                AC,
                Collections.singletonList(new GenotypeBuilder(base).alleles(Collections.singletonList(Aref)).PL(new double[]{0, -10}).GQ(100).make())});
        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(Collections.singletonList(C)).PL(haploidAltC3AllelesPL).make()).make(),
                AC,
                Collections.singletonList(new GenotypeBuilder(base).alleles(Collections.singletonList(C)).PL(new double[]{-10, 0}).GQ(100).make())});
        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(Collections.singletonList(G)).PL(haploidAltG3AllelesPL).make()).make(),
                AG,
                Collections.singletonList(new GenotypeBuilder(base).alleles(Collections.singletonList(G)).PL(new double[]{-20, 0}).GQ(200).make())});

        return tests.toArray(new Object[][]{});
    }

    /*@Test(dataProvider = "SubsetAllelesData")
    public void testSubsetAlleles(final VariantContext inputVC,
                                  final List<Allele> allelesToUse,
                                  final List<Genotype> expectedGenotypes) {
        // initialize cache of allele anyploid indices
        for (final Genotype genotype : inputVC.getGenotypes()) {
            GenotypeLikelihoods.initializeAnyploidPLIndexToAlleleIndices(inputVC.getNAlleles() - 1, genotype.getPloidy());
        }

        final GenotypesContext actual = AlleleSubsettingUtils.subsetAlleles(inputVC, allelesToUse, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN);

        Assert.assertEquals(actual.size(), expectedGenotypes.size());
        for ( final Genotype expected : expectedGenotypes ) {
            final Genotype actualGT = actual.get(expected.getSampleName());
            Assert.assertNotNull(actualGT);
            VariantContextTestUtils.assertGenotypesAreEqual(actualGT, expected);
        }
    }*/

    @DataProvider(name = "UpdateGenotypeAfterSubsettingData")
    public Object[][] makeUpdateGenotypeAfterSubsettingData() {
        final List<Object[]> tests = new ArrayList<>();

        final List<Allele> AA = Arrays.asList(Aref, Aref);
        final List<Allele> AC = Arrays.asList(Aref, C);
        final List<Allele> CC = Arrays.asList(C, C);
        final List<Allele> AG = Arrays.asList(Aref, G);
        final List<Allele> CG = Arrays.asList(C, G);
        final List<Allele> GG = Arrays.asList(G, G);
        final List<Allele> AAA = Arrays.asList(Aref, Aref, Aref);
        final List<Allele> AAC = Arrays.asList(Aref, Aref, C);
        final List<Allele> ACC = Arrays.asList(Aref, C, C);
        final List<Allele> CCC = Arrays.asList(C, C, C);
        final List<Allele> AAG = Arrays.asList(Aref, Aref, G);
        final List<Allele> ACG = Arrays.asList(Aref, C, G);
        final List<Allele> CCG = Arrays.asList(C, C, G);
        final List<Allele> AGG = Arrays.asList(Aref, G, G);
        final List<Allele> CGG = Arrays.asList(C, G, G);
        final List<Allele> GGG = Arrays.asList(G, G, G);
        final List<List<Allele>> allDiploidSubsetAlleles = Arrays.asList(AC, AG, ACG);
        final List<List<Allele>> allTriploidSubsetAlleles = Arrays.asList(AAA, AAC, ACC, CCC, AAG, ACG, CCG, AGG, CGG, GGG);

        // for P=1, the index of the genotype a is a
        final double[] aRefPL = new double[]{0.9, 0.09, 0.01};
        final double[] cPL = new double[]{0.09, 0.9, 0.01};
        final double[] gPL = new double[]{0.01, 0.09, 0.9};
        final List<double[]> allHaploidPLs = Arrays.asList(aRefPL, cPL, gPL);
        final List<List<Allele>> allHaploidSubsetAlleles = Arrays.asList(Arrays.asList(Aref), Arrays.asList(G));

        // for P=2 and N=1, the ordering is 00,01,11
        final double[] homRefPL = new double[]{0.9, 0.09, 0.01};
        final double[] hetPL = new double[]{0.09, 0.9, 0.01};
        final double[] homVarPL = new double[]{0.01, 0.09, 0.9};
        final double[] uninformative = new double[]{0.33, 0.33, 0.33};
        final List<double[]> allDiploidPLs = Arrays.asList(homRefPL, hetPL, homVarPL, uninformative);

        // for P=3 and N=2, the ordering is 000, 001, 011, 111, 002, 012, 112, 022, 122, 222
        final double[] aaaPL = new double[]{0.9, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09};
        final double[] aacPL = new double[]{0.01, 0.9, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09};
        final double[] accPL = new double[]{0.01, 0.02, 0.9, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09};
        final double[] cccPL = new double[]{0.01, 0.02, 0.03, 0.9, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09};
        final double[] aagPL = new double[]{0.01, 0.02, 0.03, 0.04, 0.9, 0.05, 0.06, 0.07, 0.08, 0.09};
        final double[] acgPL = new double[]{0.01, 0.02, 0.03, 0.04, 0.05, 0.9, 0.06, 0.07, 0.08, 0.09};
        final double[] ccgPL = new double[]{0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.9, 0.07, 0.08, 0.09};
        final double[] aggPL = new double[]{0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.9, 0.08, 0.09};
        final double[] cggPL = new double[]{0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.9, 0.09};
        final double[] gggPL = new double[]{0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.9};
        final double[] uninformativeTriploid = new double[]{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
        final List<double[]> allTriploidPLs = Arrays.asList(homRefPL, hetPL, homVarPL, uninformativeTriploid);


        for (final List<Allele> alleles : allHaploidSubsetAlleles) {
            tests.add(new Object[]{1, GenotypeAssignmentMethod.SET_TO_NO_CALL, allHaploidPLs.get(0), Arrays.asList(Aref), alleles, GATKVariantContextUtils.noCallAlleles(1)});
        }

        for (final List<Allele> alleles : allDiploidSubsetAlleles) {
            tests.add(new Object[]{2, GenotypeAssignmentMethod.SET_TO_NO_CALL, allDiploidPLs.get(0), AA, alleles, GATKVariantContextUtils.noCallAlleles(2)});
        }

        for (final List<Allele> alleles : allTriploidSubsetAlleles) {
            tests.add(new Object[]{3, GenotypeAssignmentMethod.SET_TO_NO_CALL, allTriploidPLs.get(0), AAA, alleles, GATKVariantContextUtils.noCallAlleles(3)});
        }

        final List<Allele> originalHaploidGT = Arrays.asList(Aref, C, G);
        tests.add(new Object[]{1, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, aRefPL, originalHaploidGT, originalHaploidGT, Arrays.asList(Aref)});
        tests.add(new Object[]{1, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, cPL, originalHaploidGT, originalHaploidGT, Arrays.asList(C)});
        tests.add(new Object[]{1, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, gPL, originalHaploidGT, originalHaploidGT, Arrays.asList(G)});

        for (final List<Allele> originalGT : Arrays.asList(AA, AC, CC, AG, CG, GG)) {
            tests.add(new Object[]{2, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, homRefPL, originalGT, AC, AA});
            tests.add(new Object[]{2, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, hetPL, originalGT, AC, AC});
            tests.add(new Object[]{2, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, homVarPL, originalGT, AC, CC});
        }

        for (final List<Allele> originalGT : allTriploidSubsetAlleles) {
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, aaaPL, originalGT, ACG, AAA});
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, aacPL, originalGT, ACG, AAC});
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, accPL, originalGT, ACG, ACC});
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, cccPL, originalGT, ACG, CCC});
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, aagPL, originalGT, ACG, AAG});
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, acgPL, originalGT, ACG, ACG});
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, ccgPL, originalGT, ACG, CCG});
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, aggPL, originalGT, ACG, AGG});
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, cggPL, originalGT, ACG, CCG});
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, gggPL, originalGT, ACG, GGG});
        }

        tests.add(new Object[]{1, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, allHaploidPLs.get(0), Arrays.asList(Aref, C, G), Arrays.asList(Aref), Arrays.asList(Aref)});
        tests.add(new Object[]{1, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, allHaploidPLs.get(0), Arrays.asList(Aref, C, G), Arrays.asList(C), Arrays.asList(C)});
        tests.add(new Object[]{1, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, allHaploidPLs.get(0), Arrays.asList(Aref, C, G), Arrays.asList(G), Arrays.asList(G)});

        tests.add(new Object[]{2, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, allDiploidPLs.get(0), AA, AC, AA});
        tests.add(new Object[]{2, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, allDiploidPLs.get(0), AC, AC, AC});
        tests.add(new Object[]{2, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, allDiploidPLs.get(0), CC, AC, CC});
        tests.add(new Object[]{2, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, allDiploidPLs.get(0), CG, AC, AC});

        tests.add(new Object[]{2, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, allDiploidPLs.get(0), AA, AG, AA});
        tests.add(new Object[]{2, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, allDiploidPLs.get(0), AC, AG, AA});
        tests.add(new Object[]{2, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, allDiploidPLs.get(0), CC, AG, AA});
        tests.add(new Object[]{2, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, allDiploidPLs.get(0), CG, AG, AG});

        tests.add(new Object[]{2, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, allDiploidPLs.get(0), AA, ACG, AA});
        tests.add(new Object[]{2, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, allDiploidPLs.get(0), AC, ACG, AC});
        tests.add(new Object[]{2, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, allDiploidPLs.get(0), CC, ACG, CC});
        tests.add(new Object[]{2, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, allDiploidPLs.get(0), AG, ACG, AG});
        tests.add(new Object[]{2, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, allDiploidPLs.get(0), CG, ACG, CG});
        tests.add(new Object[]{2, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, allDiploidPLs.get(0), GG, ACG, GG});

        tests.add(new Object[]{3, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, allTriploidPLs.get(0), AAA, AAC, AAA});
        tests.add(new Object[]{3, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, allTriploidPLs.get(0), ACC, AAC, ACC});
        tests.add(new Object[]{3, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, allTriploidPLs.get(0), AAC, AAC, AAC});
        tests.add(new Object[]{3, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, allTriploidPLs.get(0), AAC, ACG, AAC});
        tests.add(new Object[]{3, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, allTriploidPLs.get(0), GGG, AAA, AAA});

        return tests.toArray(new Object[][]{});
    }

    @DataProvider(name = "MakeGenotypeCallData")
    public Object[][] makeGenotypeCallData() {
        final List<Object[]> tests = new ArrayList<>();

        final List<Allele> AA = Arrays.asList(Aref,Aref);
        final List<Allele> AC = Arrays.asList(Aref,C);
        final List<Allele> CC = Arrays.asList(C,C);
        final List<Allele> AG = Arrays.asList(Aref,G);
        final List<Allele> CG = Arrays.asList(C,G);
        final List<Allele> GG = Arrays.asList(G,G);
        final List<Allele> AAA = Arrays.asList(Aref,Aref,Aref);
        final List<Allele> AAC = Arrays.asList(Aref,Aref,C);
        final List<Allele> ACC = Arrays.asList(Aref,C,C);
        final List<Allele> CCC = Arrays.asList(C,C,C);
        final List<Allele> AAG = Arrays.asList(Aref,Aref,G);
        final List<Allele> ACG = Arrays.asList(Aref,C,G);
        final List<Allele> CCG = Arrays.asList(C,C,G);
        final List<Allele> AGG = Arrays.asList(Aref,G,G);
        final List<Allele> CGG = Arrays.asList(C,G,G);
        final List<Allele> GGG = Arrays.asList(G,G,G);
        final List<List<Allele>> allDiploidSubsetAlleles = Arrays.asList(AC,AG,ACG);
        final List<List<Allele>> allTriploidSubsetAlleles = Arrays.asList(AAA,AAC,ACC,CCC,AAG,ACG,CCG,AGG,CGG,GGG);

        // for P=1, the index of the genotype a is a
        final double[] aRefPL = new double[]{0.9, 0.09, 0.01};
        final double[] cPL = new double[]{0.09, 0.9, 0.01};
        final double[] gPL = new double[]{0.01, 0.09, 0.9};
        final double[] nonRefPL = gPL;
        final List<double[]> allHaploidPLs = Arrays.asList(aRefPL, cPL, gPL);
        final List<List<Allele>> allHaploidSubsetAlleles = Arrays.asList(Collections.singletonList(Aref), Collections.singletonList(G));

        // for P=2 and N=1, the ordering is 00,01,11
        final double[] homRefPL = new double[]{0.9, 0.09, 0.01};
        final double[] hetPL = new double[]{0.09, 0.9, 0.01};
        final double[] homVarPL = new double[]{0.01, 0.09, 0.9};
        final double[] uninformative = new double[]{0.33, 0.33, 0.33};
        final List<double[]> allDiploidPLs = Arrays.asList(homRefPL, hetPL, homVarPL, uninformative);

        // for P=3 and N=2, the ordering is 000, 001, 011, 111, 002, 012, 112, 022, 122, 222
        final double[] aaaPL = new double[]{0.9, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09};
        final double[] aacPL = new double[]{0.01, 0.9, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09};
        final double[] accPL = new double[]{0.01, 0.02, 0.9, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09};
        final double[] cccPL = new double[]{0.01, 0.02, 0.03, 0.9, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09};
        final double[] aagPL = new double[]{0.01, 0.02, 0.03, 0.04, 0.9, 0.05, 0.06, 0.07, 0.08, 0.09};
        final double[] acgPL = new double[]{0.01, 0.02, 0.03, 0.04, 0.05, 0.9, 0.06, 0.07, 0.08, 0.09};
        final double[] ccgPL = new double[]{0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.9, 0.07, 0.08, 0.09};
        final double[] aggPL = new double[]{0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.9, 0.08, 0.09};
        final double[] cggPL = new double[]{0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.9, 0.09};
        final double[] gggPL = new double[]{0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.9};
        final double[] uninformativeTriploid = new double[]{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
        final List<double[]> allTriploidPLs = Arrays.asList(homRefPL, hetPL, homVarPL, uninformativeTriploid);

        for ( final List<Allele> alleles : allHaploidSubsetAlleles ) {
            tests.add(new Object[]{1, GenotypeAssignmentMethod.SET_TO_NO_CALL, allHaploidPLs.get(0), Collections.singletonList(Aref), alleles, GATKVariantContextUtils.noCallAlleles(1)});
        }

        for ( final List<Allele> alleles : allDiploidSubsetAlleles ) {
            tests.add(new Object[]{2, GenotypeAssignmentMethod.SET_TO_NO_CALL, allDiploidPLs.get(0), AA, alleles, GATKVariantContextUtils.noCallAlleles(2)});
        }

        for ( final List<Allele> alleles : allTriploidSubsetAlleles ) {
            tests.add(new Object[]{3, GenotypeAssignmentMethod.SET_TO_NO_CALL, allTriploidPLs.get(0), AAA, alleles, GATKVariantContextUtils.noCallAlleles(3)});
        }

        final List<Allele> originalHaploidGT = Collections.singletonList(Aref);
        final List<Allele> haploidAllelesToUse = Arrays.asList(Aref, C, G );
        final List<Allele> haploidAllelesWithNonRef = Arrays.asList(Aref, C, Allele.NON_REF_ALLELE );
        tests.add(new Object[]{1, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, aRefPL, originalHaploidGT, haploidAllelesToUse, Collections.singletonList(Aref)});
        tests.add(new Object[]{1, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, cPL, originalHaploidGT, haploidAllelesToUse, Collections.singletonList(C)});
        tests.add(new Object[]{1, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, gPL, originalHaploidGT, haploidAllelesToUse, Collections.singletonList(G)});
        tests.add(new Object[]{1, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, nonRefPL, originalHaploidGT, haploidAllelesWithNonRef, Collections.singletonList(Aref)});  //this used to be no-call, but new convention is no call means NO DATA

        for ( final List<Allele> originalGT : Arrays.asList(AA, AC, CC, AG, CG, GG) ) {
            tests.add(new Object[]{2, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, homRefPL, originalGT, AC, AA});
            tests.add(new Object[]{2, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, hetPL, originalGT, AC, AC});
            tests.add(new Object[]{2, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, homVarPL, originalGT, AC, CC});
        }

        for ( final List<Allele> originalGT : allTriploidSubsetAlleles) {
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, aaaPL, originalGT, ACG, AAA});
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, aacPL, originalGT, ACG, AAC});
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, accPL, originalGT, ACG, ACC});
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, cccPL, originalGT, ACG, CCC});
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, aagPL, originalGT, ACG, AAG});
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, acgPL, originalGT, ACG, ACG});
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, ccgPL, originalGT, ACG, CCG});
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, aggPL, originalGT, ACG, AGG});
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, cggPL, originalGT, ACG, CCG});
            tests.add(new Object[]{3, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, gggPL, originalGT, ACG, GGG});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "MakeGenotypeCallData")
    public void testMakeGenotypeCall(final int ploidy,
                                     final GenotypeAssignmentMethod mode,
                                     final double[] likelihoods,
                                     final List<Allele> originalGT,
                                     final List<Allele> allelesToUse,
                                     final List<Allele> expectedAlleles) {
        Utils.validateArg(originalGT.size() == ploidy, "original call must be consistent with ploidy");

        final GenotypeBuilder gb = new GenotypeBuilder("test");
        final double[] logLikelhoods = MathUtils.normalizeLog10(likelihoods);

        GATKVariantContextUtils.makeGenotypeCall(originalGT.size(), gb, mode, logLikelhoods, allelesToUse, gb.make(), null);

        final Genotype g = gb.make();
        Assert.assertEquals(g.getAlleles().size(), expectedAlleles.size());
        Assert.assertEquals(new LinkedHashSet<>(g.getAlleles()), new LinkedHashSet<>(expectedAlleles));
    }

    @Test()
    public void testSubsetToRef() {
        final Map<Genotype, Genotype> tests = new LinkedHashMap<>();

        for ( final List<Allele> alleles : Arrays.asList(Collections.singletonList(Aref), Collections.singletonList(C), Arrays.asList(Aref, C), Arrays.asList(Aref, C, C) ) ) {
            for ( final String name : Arrays.asList("test1", "test2") ) {
                final GenotypeBuilder builder = new GenotypeBuilder(name, alleles);
                builder.DP(10);
                builder.GQ(30);
                builder.AD(alleles.size() == 1 ? new int[]{1} : (alleles.size() == 2 ? new int[]{1, 2} : new int[]{1, 2, 3}));
                builder.PL(alleles.size() == 1 ? new int[]{1} : (alleles.size() == 2 ? new int[]{1, 2} : new int[]{1, 2, 3}));
                builder.attribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY,
                        alleles.size() == 1 ? new int[]{1, 2}  : (alleles.size() == 2 ? new int[]{1, 2, 3, 4} : new int[]{1, 2, 3, 4, 5, 6}));
                final List<Allele> refs = Collections.nCopies(alleles.size(), Aref);
                tests.put(builder.make(), builder.alleles(refs).noAD().noPL().attribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY, null).make());
            }
        }

        for ( final int n : Arrays.asList(1, 2, 3) ) {
            for ( final List<Genotype> genotypes : Utils.makePermutations(new ArrayList<>(tests.keySet()), n, false) ) {
                final VariantContext vc = new VariantContextBuilder("test", "20", 1, 1, Arrays.asList(Aref, C)).genotypes(genotypes).make();
                final GenotypesContext gc = GATKVariantContextUtils.subsetToRefOnly(vc, 2);

                Assert.assertEquals(gc.size(), genotypes.size());
                for ( int i = 0; i < genotypes.size(); i++ ) {
                    VariantContextTestUtils.assertGenotypesAreEqual(gc.get(i), tests.get(genotypes.get(i)));
                }
            }
        }
    }

    // --------------------------------------------------------------------------------
    //
    // Test methods for merging reference confidence VCs
    //
    // --------------------------------------------------------------------------------

    @Test(dataProvider = "totalPloidyData")
    public void testTotalPloidy(final int[] ploidies, final int defaultPloidy, final int expected) {
        final Genotype[] genotypes = new Genotype[ploidies.length];
        final List<Allele> vcAlleles = Arrays.asList(Aref,C);
        for (int i = 0; i < genotypes.length; i++)
            genotypes[i] = new GenotypeBuilder().alleles(GATKVariantContextUtils.noCallAlleles(ploidies[i])).make();
        final VariantContext vc = new VariantContextBuilder().chr("seq1").genotypes(genotypes).alleles(vcAlleles).make();
        Assert.assertEquals(GATKVariantContextUtils.totalPloidy(vc,defaultPloidy),expected," " + defaultPloidy + " " + Arrays.toString(ploidies));
    }

    @DataProvider(name="totalPloidyData")
    public Object[][] totalPloidyData() {
        final Random rdn = Utils.getRandomGenerator();
        final List<Object[]> resultList = new ArrayList<>();
        for (int i = 0; i < 100; i++) {
            final int sampleCount = rdn.nextInt(10);

            int expected = 0;
            final int defaultPloidy = rdn.nextInt(10) + 1;
            final int[] plodies = new int[sampleCount];
            for (int j = 0; j < sampleCount; j++) {
                plodies[j] = rdn.nextInt(10);
                expected += plodies[j] == 0 ? defaultPloidy : plodies[j];
            }
            resultList.add(new Object[] { plodies, defaultPloidy, expected });
        }
        return resultList.toArray(new Object[100][]);
    }

    private byte[] randomBases(final int length, final boolean reference) {
        final byte[] bases = new byte[length];
        bases[0] = (byte) (reference  ? 'A' : 'C');
        BaseUtils.fillWithRandomBases(bases, 1, bases.length);
        return bases;
    }

    @Test
    public void testCreateVariantContextWriterNoOptions() {
        final File tmpDir = createTempDir("createVCFTest");
        final File outFile = new File(tmpDir.getAbsolutePath(), "createVCFTest" + ".vcf");

        final VariantContextWriter vcw = GATKVariantContextUtils.createVCFWriter(outFile.toPath(), null, false);
        vcw.close();

        final File outFileIndex = new File(outFile.getAbsolutePath() + ".idx");
        final File outFileMD5 = new File(outFile.getAbsolutePath() + ".md5");

        Assert.assertTrue(outFile.exists(), "No output file was not created");
        Assert.assertFalse(outFileIndex.exists(), "The createIndex argument was not honored");
        Assert.assertFalse(outFileMD5.exists(), "The createMD5 argument was not honored");
    }

    @DataProvider(name="createVCFWriterData")
    public Object[][] createVCFWriterData() {
        return new Object[][]{
                {".vcf", ".idx", true, true},
                {".vcf", ".idx", false, true},
                {".vcf", ".idx", true, false},

                {".bcf", ".idx", true, true},
                {".bcf", ".idx", false, true},
                {".bcf", ".idx", true, false},

                {".vcf.bgz", ".tbi", true, true},
                // AbstractFeatureReader fails to recognize this as block compressed unless it has an accompanying index
                // Is that correct behavior ?
                //".vcf.bgz", ".tbi", false, true},
                {".vcf.gz", ".tbi", false, true},
                {".vcf.bgz", ".tbi", true, false},

                // defaults to .vcf
                {".tmp", ".idx", false, true},
        };
    }
    @Test(dataProvider = "createVCFWriterData")
    public void testCreateVCFWriterWithOptions(
            final String outputExtension,
            final String indexExtension,
            final boolean createIndex,
            final boolean createMD5) throws IOException {

        final File tmpDir = createTempDir("createVCFTest");
        final File outputFile = new File(tmpDir.getAbsolutePath(), "createVCFTest" + outputExtension);

        Options options[] = createIndex ?
                new Options[] {Options.INDEX_ON_THE_FLY} :
                new Options[] {};
        try (final VariantContextWriter writer = GATKVariantContextUtils.createVCFWriter(
                outputFile.toPath(),
                makeSimpleSequenceDictionary(),
                createMD5,
                options)) {
            writeHeader(writer);
        }

        final File outFileIndex = new File(outputFile.getAbsolutePath() + indexExtension);
        final File outFileMD5 = new File(outputFile.getAbsolutePath() + ".md5");

        Assert.assertTrue(outputFile.exists(), "No output file was not created");
        Assert.assertEquals(outFileIndex.exists(), createIndex, "The createIndex argument was not honored");
        Assert.assertEquals(outFileMD5.exists(), createMD5, "The createMD5 argument was not honored");

        verifyFileType(outputFile, outputExtension);
    }

    // just make sure we can read the file with the corresponding codec
    private void verifyFileType(
            final File resultVCFFile,
            final String outputExtension) {
        final FeatureCodec<? extends Feature, ?> featureCodec = FeatureManager.getCodecForFile(resultVCFFile.toPath());

        if (outputExtension.equals(".vcf") ||
            outputExtension.equals(".vcf.bgz") ||
            outputExtension.equals(".vcf.gz") ||
            outputExtension.equals(".tmp"))
        {
            Assert.assertEquals(featureCodec.getClass(), VCFCodec.class,
                    "Wrong codec selected for file " + resultVCFFile.getAbsolutePath());
        }
        else if (outputExtension.equals(".bcf")) {
            Assert.assertEquals(featureCodec.getClass(), BCF2Codec.class,
                    "Wrong codec selected for file " + resultVCFFile.getAbsolutePath());
        }
        else {
            throw new IllegalArgumentException("Unknown file extension in createVCFWriter test validation");
        }
    }

    @DataProvider(name="createVCFWriterLenientData")
    public Object[][] createVCFWriterLenientData() {
        return new Object[][]{
                {".vcf", ".idx", true, true},
                {".vcf", ".idx", false, true},
                {".vcf", ".idx", true, false}
        };
    }

    @Test(dataProvider = "createVCFWriterLenientData")
    public void testCreateVCFWriterLenientTrue(
            final String outputExtension,
            final String indexExtension,
            final boolean createIndex,
            final boolean createMD5) throws IOException {

        final File tmpDir = createTempDir("createVCFTest");
        final File outputFile = new File(tmpDir.getAbsolutePath(), "createVCFTest" + outputExtension);

        Options options[] = createIndex ?
                new Options[] {Options.ALLOW_MISSING_FIELDS_IN_HEADER, Options.INDEX_ON_THE_FLY} :
                new Options[] {Options.ALLOW_MISSING_FIELDS_IN_HEADER};
        try (final VariantContextWriter vcw = GATKVariantContextUtils.createVCFWriter(
                outputFile.toPath(),
                makeSimpleSequenceDictionary(),
                createMD5,
                options)) {
            writeHeader(vcw);
            writeBadVariant(vcw);  // verify leniency by writing a bogus attribute
        }

        final File outFileIndex = new File(outputFile.getAbsolutePath() + indexExtension);
        final File outFileMD5 = new File(outputFile.getAbsolutePath() + ".md5");

        Assert.assertTrue(outputFile.exists(), "No output file was not created");
        Assert.assertEquals(outFileIndex.exists(), createIndex, "The createIndex argument was not honored");
        Assert.assertEquals(outFileMD5.exists(), createMD5, "The createMD5 argument was not honored");
    }

    @Test(dataProvider = "createVCFWriterLenientData", expectedExceptions = IllegalStateException.class)
    public void testCreateVCFWriterLenientFalse(
            final String outputExtension,
            final String indexExtension, // unused
            final boolean createIndex,
            final boolean createMD5) throws IOException {

        final File tmpDir = createTempDir("createVCFTest");
        final File outputFile = new File(tmpDir.getAbsolutePath(), "createVCFTest" + outputExtension);

        Options options[] = createIndex ?
                new Options[] {Options.INDEX_ON_THE_FLY} :
                new Options[] {};
        try (final VariantContextWriter vcw = GATKVariantContextUtils.createVCFWriter(
                outputFile.toPath(),
                makeSimpleSequenceDictionary(),
                createMD5,
                options)) {
            writeHeader(vcw);
            writeBadVariant(vcw); // write a bad attribute and throw...
        }
    }

    @Test(expectedExceptions=IllegalArgumentException.class)
    public void testCreateVariantContextWriterNoReference() {
        // should throw due to lack of reference
        final File outFile = createTempFile("testVCFWriter", ".vcf");
        final VariantContextWriter vcw =
                     GATKVariantContextUtils.createVCFWriter(
                             outFile.toPath(),
                             null,
                             true,
                             Options.INDEX_ON_THE_FLY);
        vcw.close();
    }

    // This test verifies that we can write valid .gz/.tbi pair using a .gz input file with a large header.
    // Specifically, we want to make sure that index queries on the result return the first variants emitted into
    // the file, and that we don't encounter https://github.com/broadinstitute/gatk/issues/2821 and/or
    // https://github.com/broadinstitute/gatk/issues/2801.
    @Test
    public void testOnTheFlyTabixCreation() throws IOException {
        final File inputGZIPFile = new File(publicTestDir + "org/broadinstitute/hellbender/engine/8_mutect2_sorted.vcf.gz");
        final File outputGZIPFile = createTempFile("testOnTheFlyTabixCreation", ".vcf.gz");

        long recordCount = 0;
        try (final VariantContextWriter vcfWriter = GATKVariantContextUtils.createVCFWriter(
                 outputGZIPFile.toPath(),
                 null,
                 false,
                 Options.INDEX_ON_THE_FLY);
            final FeatureReader<VariantContext> inputFileReader =
                     AbstractFeatureReader.getFeatureReader(inputGZIPFile.getAbsolutePath(), new VCFCodec(), false))
        {
            vcfWriter.writeHeader((VCFHeader)inputFileReader.getHeader());
            final Iterator<VariantContext> it = inputFileReader.iterator();
            while (it.hasNext()) {
                vcfWriter.add(it.next());
                recordCount++;
            }
        }

        // make sure we got a tabix index
        final File tabixIndexFile = new File(outputGZIPFile.getAbsolutePath() + FileExtensions.TABIX_INDEX);
        Assert.assertTrue(tabixIndexFile.exists());
        Assert.assertTrue(tabixIndexFile.length() > 0);

        // verify the index via query
        final SimpleInterval queryInterval = new SimpleInterval("chr6:33414233-118314029");
        try (final FeatureReader<VariantContext> outputFileReader = AbstractFeatureReader.getFeatureReader(
                outputGZIPFile.getAbsolutePath(),
                new VCFCodec(),
                true)) // require index
        {
            final long actualCount = outputFileReader.query(
                    queryInterval.getContig(), queryInterval.getStart(), queryInterval.getEnd()).stream().count();
            Assert.assertEquals(actualCount, recordCount);
        }
    }

    @Test
    public void testCreateVcfWriterOnNio() throws IOException {
        final File inputGZIPFile = new File(
            publicTestDir + "org/broadinstitute/hellbender/engine/8_mutect2_sorted.vcf.gz");
        try (FileSystem jimfs = Jimfs.newFileSystem(Configuration.unix())) {
            final Path outputGZIP = jimfs.getPath("testCreateVcfWriterOnNio.vcf.gz");
            final Path tabixIndex = outputGZIP.resolveSibling(outputGZIP.getFileName().toString() + FileExtensions.TABIX_INDEX);
            long recordCount = 0;

            try (final VariantContextWriter vcfWriter = GATKVariantContextUtils.createVCFWriter(
                outputGZIP,
                null,
                false,
                Options.INDEX_ON_THE_FLY);
                final FeatureReader<VariantContext> inputFileReader =
                    AbstractFeatureReader
                        .getFeatureReader(inputGZIPFile.getAbsolutePath(), new VCFCodec(), false)) {
                vcfWriter.writeHeader((VCFHeader) inputFileReader.getHeader());
                final Iterator<VariantContext> it = inputFileReader.iterator();
                while (it.hasNext()) {
                    vcfWriter.add(it.next());
                    recordCount++;
                }
            }

            // make sure we got an output and tabix index
            Assert.assertTrue(Files.exists(outputGZIP));
            Assert.assertTrue(Files.exists(tabixIndex));
            Assert.assertTrue(Files.size(tabixIndex) > 0);

            // make sure we got an output and queryable index
            long roundTripRecordCount = 0;
            try (final VCFFileReader outputFileReader = new VCFFileReader(outputGZIP, true)) {
                final Iterator<VariantContext> it = outputFileReader.query(new SimpleInterval("chr6", 1, 999999999));
                while (it.hasNext()) {
                    it.next();
                    roundTripRecordCount++;
                }
            }
            Assert.assertEquals(roundTripRecordCount, recordCount);
        }
    }

    private void writeHeader(final VariantContextWriter writer) {
        final Set<VCFHeaderLine> metaData = new HashSet<>();
        metaData.add(new VCFHeaderLine(
                VCFHeaderVersion.VCF4_2.getFormatString(),
                VCFHeaderVersion.VCF4_2.getVersionString()));
        final VCFHeader vcfHeader = new VCFHeader(metaData, Collections.emptyList());
        vcfHeader.setSequenceDictionary(makeSimpleSequenceDictionary());
        writer.writeHeader(vcfHeader);
    }

    private void writeBadVariant(final VariantContextWriter writer) {
        //write a variant with a (bad) attribute that doesn't appear in the header to the output
        final VariantContextBuilder vcBuilder = new VariantContextBuilder("","chr1", 1, 1, Arrays.asList(Aref));
        vcBuilder.attribute("fake", new Object());
        final VariantContext vc = vcBuilder.make();
        writer.add(vc);
    }

    private SAMSequenceDictionary makeSimpleSequenceDictionary() {
        final SAMSequenceDictionary seqDictionary = new SAMSequenceDictionary();
        seqDictionary.addSequence(new SAMSequenceRecord("chr1", 10));
        return seqDictionary;
    }

    @DataProvider(name = "setFilteredGenotypeToNocallTest")
    public Object[][] makeSetFilteredGenotypeToNocallTest() {
        List<Object[]> tests = new ArrayList<>();

        tests.add(new Object[]{new ArrayList<>(Arrays.asList(Aref, C)), true, Arrays.asList(0), 0, Arrays.asList(0.5)});
        tests.add(new Object[]{new ArrayList<>(Arrays.asList(Aref, C)), false, Arrays.asList(1), 2, Arrays.asList(0.5)});
        tests.add(new Object[]{new ArrayList<>(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL)), true, Arrays.asList(1), 2, Arrays.asList(0.5)});

        return tests.toArray(new Object[][]{});
    }

    private List<String> getGenotypeFilters(final VariantContext vc, final Genotype g) {
        final List<String> filters = new ArrayList<>();
        if (g.isFiltered()) {
            filters.add(g.getFilters());
        }

        return filters;
    }

    @Test(dataProvider="setFilteredGenotypeToNocallTest")
    public void testSetFilteredGenotypeToNocall(final List<Allele> genotypeAlleles, final boolean setFilteredGenotypesToNocall,
                                                final List<Integer> calledAltAlleles, final int calledAlleles, final List<Double> alleleFrequency) {
        final List<Allele> alleles = new ArrayList<>(Arrays.asList(Aref, C));
        final Genotype genotype = new GenotypeBuilder("NA12878").alleles(genotypeAlleles).filter("FILTER").make();
        final VariantContextBuilder builder = new VariantContextBuilder("test", "chr1", 1, Aref.length(), alleles).genotypes(genotype).
                attribute(VCFConstants.ALLELE_COUNT_KEY, new int[]{1} ).
                attribute(VCFConstants.ALLELE_NUMBER_KEY, 2).
                attribute(VCFConstants.ALLELE_FREQUENCY_KEY, new double[]{0.5});
        GATKVariantContextUtils.setFilteredGenotypeToNocall(builder, builder.make(), setFilteredGenotypesToNocall, this::getGenotypeFilters);
        VariantContext vc = builder.make();
        if ( vc.hasAttribute((VCFConstants.ALLELE_COUNT_KEY)) ) {
            final int[] array = (int[]) vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY);
            Assert.assertEquals(vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY), calledAltAlleles.toArray());
        }
        if ( vc.hasAttribute((VCFConstants.ALLELE_NUMBER_KEY)) ) {
            Assert.assertEquals(vc.getAttribute(VCFConstants.ALLELE_NUMBER_KEY), calledAlleles);
        }
        if ( vc.hasAttribute((VCFConstants.ALLELE_FREQUENCY_KEY)) ) {
            Assert.assertEquals(vc.getAttribute(VCFConstants.ALLELE_FREQUENCY_KEY), alleleFrequency.toArray());
        }
    }

    @DataProvider(name="gqFromPLsData")
    public Object[][] gqFromPLsData() {
        return new Object[][]{
                {new int[]{0, 15}, 15},
                {new int[]{15, 0}, 15},
                {new int[]{0, 10, 20}, 10},
                {new int[]{20, 10, 0}, 10},
                {new int[]{0, 10, 20, 30, 40}, 10},
                {new int[]{30, 40, 20, 10, 0}, 10},
                {new int[]{-10, 20, 35}, 30},
                {new int[]{35, 40, -10, 15, 20}, 25},
                {new int[]{0, 10, 20, 30, 40, 50, 5}, 5},
                {new int[]{15, 15, 0, 5}, 5},
                {new int[]{15, 15, 0, 25}, 15},
                {new int[]{0, 15, 0, 25}, 0}
        };
    }

    @Test(dataProvider = "gqFromPLsData")
    public void testCalculateGQFromPLs(final int[] plValues, final int expectedGQ) {
        Assert.assertEquals(GATKVariantContextUtils.calculateGQFromPLs(plValues), expectedGQ);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testCalculateGQFromShortPLArray() {
        final int[] plValues = new int[]{0};
        GATKVariantContextUtils.calculateGQFromPLs(plValues);
    }

    @Test
    public void testCreateFilterListWithAppend() {
        final List<Allele> alleles = new ArrayList<>();
        alleles.add(Cref);
        alleles.add(Allele.create("A", false));
        final VariantContextBuilder vcBuilder = new VariantContextBuilder("","chr1", 1, 1, alleles);
        vcBuilder.filters("F1", "F2", "Y");
        final VariantContext vc = vcBuilder.make();

        final String testFilterString = "TEST_FILTER";
        final List<String> filterResult = GATKVariantContextUtils.createFilterListWithAppend(vc, testFilterString);
        Assert.assertEquals(filterResult.size(), 4);
        Assert.assertTrue(vc.getFilters().stream().allMatch(f -> filterResult.contains(f)));
        Assert.assertTrue(filterResult.contains(testFilterString));

        // Check that it is in the proper position
        Assert.assertEquals(filterResult.get(2), testFilterString);
    }

    @Test
    public void testCreateFilterListWithAppendToEmpty() {
        final List<Allele> alleles = new ArrayList<>();
        alleles.add(Cref);
        alleles.add(Allele.create("A", false));
        final VariantContextBuilder vcBuilder = new VariantContextBuilder("","chr1", 1, 1, alleles);
        final VariantContext vc = vcBuilder.make();

        final String testFilterString = "TEST_FILTER";
        final List<String> filterResult = GATKVariantContextUtils.createFilterListWithAppend(vc, testFilterString);
        Assert.assertEquals(filterResult.size(), 1);
        Assert.assertTrue(filterResult.contains(testFilterString));
    }

    @DataProvider
    Object[][] provideAllelesAndFrameshiftResults() {
        return new Object[][] {
                { Allele.create((byte)'A'), Allele.create((byte)'A'), false },
                { Allele.create((byte)'A'), Allele.create((byte)'T'), false },
                {
                        Allele.create(new byte[] {(byte)'A',(byte)'A'}),
                        Allele.create(new byte[] {(byte)'A',(byte)'A'}),
                        false
                },
                {
                        Allele.create(new byte[] {(byte)'A',(byte)'A'}),
                        Allele.create(new byte[] {(byte)'A',(byte)'T'}),
                        false
                },
                {
                        Allele.create(new byte[] {(byte)'A',(byte)'A',(byte)'A'}),
                        Allele.create(new byte[] {(byte)'A',(byte)'T',(byte)'T'}),
                        false
                },
                {
                        Allele.create(new byte[] {(byte)'A',(byte)'A',(byte)'A',(byte)'A'}),
                        Allele.create(new byte[] {(byte)'A'}),
                        false
                },
                {
                        Allele.create(new byte[] {(byte)'A'}),
                        Allele.create(new byte[] {(byte)'A',(byte)'A',(byte)'A',(byte)'A'}),
                        false
                },

                // ======================
                {
                        Allele.create(new byte[] {(byte)'A',(byte)'A'}),
                        Allele.create(new byte[] {(byte)'A'}),
                        true
                },
                {
                        Allele.create(new byte[] {(byte)'A',(byte)'A',(byte)'A'}),
                        Allele.create(new byte[] {(byte)'A'}),
                        true
                },
                {
                        Allele.create(new byte[] {(byte)'A',(byte)'A',(byte)'A'}),
                        Allele.create(new byte[] {(byte)'A',(byte)'A'}),
                        true
                },
                {
                        Allele.create(new byte[] {(byte)'A',(byte)'A',(byte)'A',(byte)'A'}),
                        Allele.create(new byte[] {(byte)'A',(byte)'A'}),
                        true
                },
                {
                        Allele.create(new byte[] {(byte)'A',(byte)'A',(byte)'A',(byte)'A'}),
                        Allele.create(new byte[] {(byte)'A',(byte)'A',(byte)'A'}),
                        true
                },

                {
                        Allele.create(new byte[] {(byte)'A'}),
                        Allele.create(new byte[] {(byte)'A',(byte)'A'}),
                        true
                },
                {
                        Allele.create(new byte[] {(byte)'A'}),
                        Allele.create(new byte[] {(byte)'A',(byte)'A',(byte)'A'}),
                        true
                },
                {
                        Allele.create(new byte[] {(byte)'A',(byte)'A'}),
                        Allele.create(new byte[] {(byte)'A',(byte)'A',(byte)'A'}),
                        true
                },
                {
                        Allele.create(new byte[] {(byte)'A',(byte)'A'}),
                        Allele.create(new byte[] {(byte)'A',(byte)'A',(byte)'A',(byte)'A'}),
                        true
                },
                {
                        Allele.create(new byte[] {(byte)'A',(byte)'A',(byte)'A'}),
                        Allele.create(new byte[] {(byte)'A',(byte)'A',(byte)'A',(byte)'A'}),
                        true
                },
        };
    }

    @DataProvider
    Object[][] providePositionsAndFrameshiftResults() {
        return new Object[][] {
                { 1,1,1, false },
                { 1,3,1, true },
                { 1,3,2, true },
                { 1,3,3, false },
                { 1,3,233, true },
                { 1,3,234, false },
                { 1,3,235, true },
                { 8,9,8, true },
                { 8,9,9, false },
                { 8,9,10, true },
                { 8,9,11, true },
                { 8,9,12, false },
        };
    }

    @DataProvider
    Object[][] provideDataForTestIsInsertion() {
        return new Object[][] {
                { Allele.create("A", true),     Allele.create("T"),     false },
                { Allele.create("A", true),     Allele.create("TT"),    true },
                { Allele.create("AA", true),    Allele.create("TT"),    false },
                { Allele.create("AA", true),    Allele.create("T"),     false },
                { Allele.create("A", true),     Allele.create("TTTTT"), true },
                { Allele.create("AAAAA", true), Allele.create("T"),     false },
                { Allele.create("AAAAA", true), Allele.create("TTTTT"), false },
        };
    }

    @DataProvider
    Object[][] provideDataForTestIsDeletion() {
        return new Object[][] {
                { Allele.create("A", true),     Allele.create("T"),     false },
                { Allele.create("A", true),     Allele.create("TT"),    false },
                { Allele.create("AA", true),    Allele.create("TT"),    false },
                { Allele.create("AA", true),    Allele.create("T"),     true },
                { Allele.create("A", true),     Allele.create("TTTTT"), false },
                { Allele.create("AAAAA", true), Allele.create("T"),     true },
                { Allele.create("AAAAA", true), Allele.create("TTTTT"), false },
        };
    }

    @DataProvider
    Object[][] provideDataForTestIsOnp() {
        return new Object[][] {
                { Allele.create("A", true),     Allele.create("T"),     true },
                { Allele.create("A", true),     Allele.create("TT"),    false },
                { Allele.create("AA", true),    Allele.create("TT"),    true },
                { Allele.create("AA", true),    Allele.create("T"),     false },
                { Allele.create("A", true),     Allele.create("TTTTT"), false },
                { Allele.create("AAAAA", true), Allele.create("T"),     false },
                { Allele.create("AAAAA", true), Allele.create("TTTTT"), true },
        };
    }

    @DataProvider
    Object[][] provideDataForTestIsIndel() {
        return new Object[][] {
                { Allele.create("A", true),     Allele.create("T"),     false },
                { Allele.create("A", true),     Allele.create("TT"),    true },
                { Allele.create("AA", true),    Allele.create("TT"),    false },
                { Allele.create("AA", true),    Allele.create("T"),     true },
                { Allele.create("A", true),     Allele.create("TTTTT"), true },
                { Allele.create("AAAAA", true), Allele.create("T"),     true },
                { Allele.create("AAAAA", true), Allele.create("TTTTT"), false },
        };
    }

    @DataProvider
    Object[][] provideDataForIsTransition() {
        return new Object[][] {
                { makeVC("source", Arrays.asList(Aref, G)), true},
                { makeVC("source", Arrays.asList(Cref, T)), true},
                { makeVC("source", Arrays.asList(Gref, A)), true},
                { makeVC("source", Arrays.asList(Tref, C)), true},

                { makeVC("source", Arrays.asList(Aref, C)), false},
                { makeVC("source", Arrays.asList(Aref, T)), false},

                { makeVC("source", Arrays.asList(Cref, A)), false},
                { makeVC("source", Arrays.asList(Cref, G)), false},

                { makeVC("source", Arrays.asList(Gref, C)), false},
                { makeVC("source", Arrays.asList(Gref, T)), false},

                { makeVC("source", Arrays.asList(Tref, A)), false},
                { makeVC("source", Arrays.asList(Tref, G)), false}

        };
    }


    @Test(dataProvider = "provideAllelesAndFrameshiftResults")
    void testIsFrameshift(final Allele ref, final Allele alt, final boolean expected) {
        Assert.assertEquals( GATKVariantContextUtils.isFrameshift(ref, alt), expected );
        Assert.assertEquals( GATKVariantContextUtils.isFrameshift(ref.getBaseString(), alt.getBaseString()), expected );
    }

    @Test(dataProvider = "providePositionsAndFrameshiftResults")
    void testIsFrameshiftByPositions(final int refStart, final int refEnd, final int altEnd, final boolean expected) {
        Assert.assertEquals( GATKVariantContextUtils.isFrameshift(refStart, refEnd, altEnd), expected );
    }

    @Test(dataProvider = "provideDataForTestIsInsertion")
    void testIsInsertion(final Allele ref, final Allele alt, final boolean expected) {
        Assert.assertEquals( GATKVariantContextUtils.isInsertion(ref, alt), expected );
        Assert.assertEquals( GATKVariantContextUtils.isInsertion(ref.getBaseString(), alt.getBaseString()), expected );
    }

    @Test(dataProvider = "provideDataForTestIsDeletion")
    void testIsDeletion(final Allele ref, final Allele alt, final boolean expected) {
        Assert.assertEquals( GATKVariantContextUtils.isDeletion(ref, alt), expected );
        Assert.assertEquals( GATKVariantContextUtils.isDeletion(ref.getBaseString(), alt.getBaseString()), expected );
    }

    @Test(dataProvider = "provideDataForTestIsOnp")
    void testIsOnp(final Allele ref, final Allele alt, final boolean expected) {
        Assert.assertEquals( GATKVariantContextUtils.isXnp(ref, alt), expected );
        Assert.assertEquals( GATKVariantContextUtils.isXnp(ref.getBaseString(), alt.getBaseString()), expected );
    }

    @Test(dataProvider = "provideDataForTestIsIndel")
    void testIsIndel(final Allele ref, final Allele alt, final boolean expected ) {
        Assert.assertEquals( GATKVariantContextUtils.isIndel(ref, alt), expected );
        Assert.assertEquals( GATKVariantContextUtils.isIndel(ref.getBaseString(), alt.getBaseString()), expected );
    }

    @Test(dataProvider = "provideDataForIsTransition")
    public void testIsTransition(final VariantContext vc, final boolean isTransition) {
        Assert.assertEquals(GATKVariantContextUtils.isTransition(vc), isTransition);
    }


    @DataProvider
    public Object[][] provideMatchAlleles() {
        // These were chosen to correspond to test cases in the test exac datasource VCF.
        return new Object[][] {
                {new SimpleInterval("3", 69521, 69521), Arrays.asList("T", "C"),
                        new SimpleInterval("3", 69521, 69521), Arrays.asList("T", "A", "C"),
                        new int[]{1}},
                {new SimpleInterval("3", 69521, 69521), Arrays.asList("T", "A"),
                        new SimpleInterval("3", 69521, 69521), Arrays.asList("T", "A", "C"),
                        new int[]{0}},
                {new SimpleInterval("3", 69521, 69521), Arrays.asList("T", "C", "A"),
                        new SimpleInterval("3", 69521, 69521), Arrays.asList("T", "A", "C"),
                        new int[]{1, 0}},
                {new SimpleInterval("3", 69521, 69521), Arrays.asList("T", "C", "A", "G"),
                        new SimpleInterval("3", 69521, 69521), Arrays.asList("T", "A", "C"),
                        new int[]{1, 0, -1}},
                {new SimpleInterval("3", 69552, 69552), Arrays.asList("G", "A"),
                        new SimpleInterval("3", 69552, 69552), Arrays.asList("G", "T", "A", "C"),
                        new int[]{1}},
                {new SimpleInterval("3", 69521, 69521), Arrays.asList("G", "T"),
                        new SimpleInterval("3", 69521, 69521), Arrays.asList("G", "T", "A", "C"),
                        new int[]{0}},
                {new SimpleInterval("3", 69550, 69550), Arrays.asList("G", "C"),
                        new SimpleInterval("3", 69550, 69550), Arrays.asList("G", "T", "A", "C"),
                        new int[]{2}},
                {new SimpleInterval("3", 324682, 324714), Arrays.asList("ACCAGGCCCAGCTCATGCTTCTTTGCAGCCTCT", "A"),
                        new SimpleInterval("3", 324682, 324714), Arrays.asList("ACCAGGCCCAGCTCATGCTTCTTTGCAGCCTCT", "TCCAGGCCCAGCTCATGCTTCTTTGCAGCCTCT", "A"),
                        new int[]{1}},
                {new SimpleInterval("3", 324682, 324714), Arrays.asList("ACCAGGCCCAGCTCATGCTTCTTTGCAGCCTCT", "TCCAGGCCCAGCTCATGCTTCTTTGCAGCCTCT", "A"),
                        new SimpleInterval("3", 324682, 324714), Arrays.asList("ACCAGGCCCAGCTCATGCTTCTTTGCAGCCTCT", "TCCAGGCCCAGCTCATGCTTCTTTGCAGCCTCT", "A"),
                        new int[]{0, 1}},
                //HARD!!  Same as the previous test
                {new SimpleInterval("3", 324682, 324682), Arrays.asList("A", "T"),
                        new SimpleInterval("3", 324682, 324714), Arrays.asList("ACCAGGCCCAGCTCATGCTTCTTTGCAGCCTCT", "TCCAGGCCCAGCTCATGCTTCTTTGCAGCCTCT", "A"),
                        new int[]{0}},
                {new SimpleInterval("3", 324683, 324683), Arrays.asList("C", "T"), //(See second base in ref and alt)
                        new SimpleInterval("3", 324682, 324714), Arrays.asList("TCCAGGCCCAGCTCATGCTTCTTTGCAGCCTCT", "TTCAGGCCCAGCTCATGCTTCTTTGCAGCCTCT", "A"),
                        new int[]{0}},
                {new SimpleInterval("3", 13372, 13372), Arrays.asList("G", "C"),
                        new SimpleInterval("3", 13372, 13372), Arrays.asList("G", "C"),
                        new int[]{0}},
                {new SimpleInterval("3", 13372, 13372), Arrays.asList("G", "C", "T"),
                        new SimpleInterval("3", 13372, 13372), Arrays.asList("G", "C"),
                        new int[]{0, -1}},
                {new SimpleInterval("3", 13372, 13372), Arrays.asList("G", "T", "C"),
                        new SimpleInterval("3", 13372, 13372), Arrays.asList("G", "C"),
                        new int[]{-1, 0}},
                {new SimpleInterval("3", 13372, 13372), Arrays.asList("G", "GTT", "GT"),
                        new SimpleInterval("3", 13371, 13372), Arrays.asList("AG", "AC", "AGTT"),
                        new int[]{1, -1}},

        };
    }
    @Test(dataProvider = "provideMatchAlleles")
    public void testMatchAlleles(final SimpleInterval variant1Interval, final List<String> variant1Alleles,
                                 final SimpleInterval variant2Interval, final List<String> variant2Alleles,
                                 final int[] gtMatch) {
        final VariantContext variant1 = new VariantContextBuilder()
                .chr(variant1Interval.getContig()).start(variant1Interval.getStart()).stop(variant1Interval.getEnd())
                .alleles(variant1Alleles)
                .make();

        final VariantContext variant2 = new VariantContextBuilder()
                .chr(variant1Interval.getContig()).start(variant2Interval.getStart()).stop(variant2Interval.getEnd())
                .alleles(variant2Alleles)
                .make();

        final int[] matches = GATKVariantContextUtils.matchAllelesAndStart(variant1, variant2);
        Assert.assertTrue(Arrays.equals(matches, gtMatch), "Failed");
    }

    @DataProvider(name = "multiAllelic")
    public Object[][] multiAllelicDataProvider() {
        Genotype g1 = makeG("sample1", Aref, T, -2.5,25,0,25,25,25,25);
        Genotype g1Split=makeG("sample1",Aref,T,-2.5,25,0,25);
        Genotype g1Ref = makeG("sample1", Aref, Aref, 0,0,0,0);

        Genotype g2 = makeG("sample2", Aref, ATC, -2.5,25,25,25,0,25,25);
        Genotype g2Split = makeG("sample2", Aref, ATC, -2.5,25,0,25);
        Genotype g2Ref = makeG("sample2", Aref, Aref, 0,0,0,0);

        Genotype g3 = makeG("sample3", Aref, T, -2.5,25,0,25,25,25,25);
        Genotype g3Split = makeG("sample3", Aref, T, -2.5,25,0,25);
        Genotype g3Ref = makeG("sample3", Aref, Aref, 0,0,0,0);

        Genotype g4 = makeG("sample3", Aref, Aref, -2.5,0,25,25,25,25,25);
        Genotype g4Split = makeG("sample3", Aref, Aref, -2.5,0,25,25);

        GenotypesContext gc1 = GenotypesContext.create(g1, g2, g3, g4);

        GenotypesContext gc2 = GenotypesContext.create(g1Split, g2Ref, g3Split, g4Split);
        GenotypesContext gc3 = GenotypesContext.create(g1Ref, g2Split, g3Ref, g4Split);

        VariantContext vcIn = new VariantContextBuilder("source", "1", 10, 10, Arrays.asList(Aref, T, ATC)).genotypes(gc1)
                .attribute("AC", Arrays.asList(2, 1)).attribute("AF", Arrays.asList(0.25, 0.125)).attribute("AN", 8).attribute("DP",10).make();

        VariantContext expectedVc1 = new VariantContextBuilder("source", "1", 10, 10, Arrays.asList(Aref, T)).genotypes(gc2)
                .attribute("AC", 2).attribute("AF", 0.25).attribute("AN", 8).make();
        VariantContext expectedVc2 = new VariantContextBuilder("source", "1", 10, 10, Arrays.asList(Aref, ATC)).genotypes(gc3)
                .attribute("AC", 1).attribute("AF", 0.125).attribute("AN", 8).make();

        VariantContext expectedVc1Keep = new VariantContextBuilder("source", "1", 10, 10, Arrays.asList(Aref, T)).genotypes(gc2)
                .attribute("AC", 2).attribute("AF", 0.25).attribute("AN", 8).attribute("AC_Orig", 2).attribute("AF_Orig", 0.25).attribute("AN_Orig", 8).make();

        VariantContext expectedVc2Keep = new VariantContextBuilder("source", "1", 10, 10, Arrays.asList(Aref, ATC)).genotypes(gc3)
                .attribute("AC", 1).attribute("AF", 0.125).attribute("AN", 8).attribute("AC_Orig", 1).attribute("AF_Orig", 0.125).attribute("AN_Orig", 8).make();


        return new Object[][]{
                {vcIn, new LinkedList<VariantContext>(Arrays.asList(expectedVc1, expectedVc2)), false},
                {vcIn, new LinkedList<VariantContext>(Arrays.asList(expectedVc1Keep, expectedVc2Keep)), true}
        };
    }

    @Test(dataProvider = "multiAllelic")
    public void testSplitMultiAllelic(final VariantContext vcToSplit, final List<VariantContext> expectedVariants, Boolean keepOriginalChrCounts) {
        final List<VariantContext> outVcs = GATKVariantContextUtils.splitVariantContextToBiallelics(vcToSplit, true, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, keepOriginalChrCounts);
        Assert.assertEquals(outVcs.size(), expectedVariants.size());
        for (int i = 0; i < outVcs.size(); i++) {
            VariantContextTestUtils.assertVariantContextsAreEqual(outVcs.get(i), expectedVariants.get(i), new ArrayList<String>(), Collections.emptyList());
        }
    }

    @Test(dataProvider = "SplitBiallelics")
    public void testSplitBiallelicsNoGenotypes(final VariantContext vc, final List<VariantContext> expectedBiallelics) {
        final List<VariantContext> biallelics = GATKVariantContextUtils.splitVariantContextToBiallelics(vc, false, GenotypeAssignmentMethod.SET_TO_NO_CALL, false);
        Assert.assertEquals(biallelics.size(), expectedBiallelics.size());
        for (int i = 0; i < biallelics.size(); i++) {
            final VariantContext actual = biallelics.get(i);
            final VariantContext expected = expectedBiallelics.get(i);
            VariantContextTestUtils.assertVariantContextsAreEqual(actual, expected, new ArrayList<String>(), Collections.emptyList());
        }
    }

    @Test(dataProvider = "SplitBiallelics", dependsOnMethods = "testSplitBiallelicsNoGenotypes")
    public void testSplitBiallelicsGenotypes(final VariantContext vc, final List<VariantContext> expectedBiallelics) {
        final List<Genotype> genotypes = new ArrayList<Genotype>();

        int sampleI = 0;
        for (final List<Allele> alleles : Utils.makePermutations(vc.getAlleles(), 2, true)) {
            genotypes.add(GenotypeBuilder.create("sample" + sampleI++, alleles));
        }
        genotypes.add(GenotypeBuilder.createMissing("missing", 2));

        final VariantContext vcWithGenotypes = new VariantContextBuilder(vc).genotypes(genotypes).make();

        final List<VariantContext> biallelics = GATKVariantContextUtils.splitVariantContextToBiallelics(vcWithGenotypes, false, GenotypeAssignmentMethod.SET_TO_NO_CALL, false);
        for (int i = 0; i < biallelics.size(); i++) {
            final VariantContext actual = biallelics.get(i);
            Assert.assertEquals(actual.getNSamples(), vcWithGenotypes.getNSamples()); // not dropping any samples

            for (final Genotype inputGenotype : genotypes) {
                final Genotype actualGenotype = actual.getGenotype(inputGenotype.getSampleName());
                Assert.assertNotNull(actualGenotype);
                if (!vc.isVariant() || vc.isBiallelic())
                    Assert.assertEquals(actualGenotype, vcWithGenotypes.getGenotype(inputGenotype.getSampleName()));
                else
                    Assert.assertTrue(actualGenotype.isNoCall());
            }
        }
    }

    /** Makes a list of alleles, with the first one being marked as reference. */
    private static List<Allele> makeAlleles(String... ss) {
        final List<Allele> alleles = new ArrayList<>();
        for (int i=0; i<ss.length; ++i) {
            alleles.add(Allele.create(ss[i], i== 0));
        }
        return alleles;
    }

    @Test()
    public void testIsMnpWithoutNonRef() {
        final String NON_REF = Allele.NON_REF_STRING;
        Assert.assertTrue(GATKVariantContextUtils.isUnmixedMnpIgnoringNonRef(makeVC(makeAlleles("ACT", "GCG"), 1)));
        Assert.assertTrue(GATKVariantContextUtils.isUnmixedMnpIgnoringNonRef(makeVC(makeAlleles("AA", "GG", "CT"), 1)));
        Assert.assertTrue(GATKVariantContextUtils.isUnmixedMnpIgnoringNonRef(makeVC(makeAlleles("ACT", "GCG", NON_REF), 1)));
        Assert.assertTrue(GATKVariantContextUtils.isUnmixedMnpIgnoringNonRef(makeVC(makeAlleles("AA", "GG", "CT", NON_REF), 1)));

        Assert.assertFalse(GATKVariantContextUtils.isUnmixedMnpIgnoringNonRef(makeVC(makeAlleles("A", "C"), 1)));
        Assert.assertFalse(GATKVariantContextUtils.isUnmixedMnpIgnoringNonRef(makeVC(makeAlleles("A", "C", NON_REF), 1)));
        Assert.assertFalse(GATKVariantContextUtils.isUnmixedMnpIgnoringNonRef(makeVC(makeAlleles("A", "AC", NON_REF), 1)));
        Assert.assertFalse(GATKVariantContextUtils.isUnmixedMnpIgnoringNonRef(makeVC(makeAlleles("ACT", "A", NON_REF), 1)));
        Assert.assertFalse(GATKVariantContextUtils.isUnmixedMnpIgnoringNonRef(makeVC(makeAlleles("ACT", "A", "AGG", NON_REF), 1)));
    }

    @Test(dataProvider = "variantTypes")
    public void testVariantTypesAndIsComplex(final String ref, final String alt, final VariantContext.Type gtType, boolean isComplexIndel) {
        Assert.assertEquals(GATKVariantContextUtils.typeOfVariant(Allele.create(ref), Allele.create(alt)), gtType);
        Assert.assertEquals(GATKVariantContextUtils.isComplexIndel(Allele.create(ref), Allele.create(alt)), isComplexIndel);
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testSymbolicRef() {
        GATKVariantContextUtils.typeOfVariant(Allele.NON_REF_ALLELE, Allele.create("C"));
    }

    @DataProvider(name = "variantTypes")
    public Object[][] variantTypes() {
        return new Object[][]{
                // ref, alt, type, isComplex?
                {"CCTTGGCTTATTCCA", "C", VariantContext.Type.INDEL, false},
                {"C", "CCTTGGCTTATTCCA", VariantContext.Type.INDEL, false},
                {"ACTAG", "A", VariantContext.Type.INDEL, false},
                {"ATT", "AT", VariantContext.Type.INDEL, false},
                {"AT", "ATT", VariantContext.Type.INDEL, false},
                {"CT", "CAGG", VariantContext.Type.INDEL, true},
                {"CTTT", "CAGG", VariantContext.Type.MNP, false},
                {"CTTT", "CAGGG", VariantContext.Type.INDEL, true},
                {"T", "T", VariantContext.Type.NO_VARIATION, false},
                {"CTAG", "CTAG", VariantContext.Type.NO_VARIATION, false},
                {"A", "AAGAAGCATGC", VariantContext.Type.INDEL, false},
                {"A", "C", VariantContext.Type.SNP, false},
                {"AG", "CA", VariantContext.Type.MNP, false},
                {"AGAAGG", "CATTCC", VariantContext.Type.MNP, false},
                {"GC", "GA", VariantContext.Type.SNP, false},
                {"GA", "<NON_REF>", VariantContext.Type.SYMBOLIC, false},
                {"GA", "*", VariantContext.Type.NO_VARIATION, false},

                // There are two MNPs here
                {"AGAAGG", "CATACC", VariantContext.Type.MNP, false},

                // Note that this is technically a simple AT insertion, but the isComplex cannot handle this properly.
                {"CT", "CATT", VariantContext.Type.INDEL, true},
        };
    }

    @Test(dataProvider = "bestAlleleSNP")
    public void testChooseBestAlleleSNP(final Allele referenceAllele, final List<Allele> altAlleles, int offsetIntoRead, int minBaseQualityCutoff, int gtChosenInAlt, boolean isGtRef) {

        // No indels
        if (isGtRef) {
            final PileupElement pileupElement = ArtificialReadUtils.createNonIndelPileupElement(offsetIntoRead, referenceAllele, GATKVariantContextUtilsUnitTest.LENGTH_OF_ARTIFICIAL_READ);
            final Allele chosen = GATKVariantContextUtils.chooseAlleleForRead(pileupElement, referenceAllele, altAlleles, minBaseQualityCutoff);
            Assert.assertEquals(chosen, referenceAllele);
        } else {
            final PileupElement pileupElement = ArtificialReadUtils.createNonIndelPileupElement(offsetIntoRead, altAlleles.get(gtChosenInAlt), GATKVariantContextUtilsUnitTest.LENGTH_OF_ARTIFICIAL_READ);
            final Allele chosen = GATKVariantContextUtils.chooseAlleleForRead(pileupElement, referenceAllele, altAlleles, minBaseQualityCutoff);
            Assert.assertEquals(chosen, altAlleles.get(gtChosenInAlt));
        }
    }

    @DataProvider(name = "bestAlleleSNP")
    public Object[][] createSNPAlleles() {
        return new Object[][] {
                // final Allele referenceAllele, final List<Allele> altAlleles, int offsetIntoRead, int minBaseQualityCutoff, int gtChosenInAlt, boolean isGtRef
                {
                        Allele.create("A", true), Arrays.asList(Allele.create("T", false), Allele.create("CA", false)), 0, 0,
                        0, false
                }, {
                Allele.create("A", true), Arrays.asList(Allele.create("T", false), Allele.create("CA", false)), 0,  0,
                0, true
        }, {
                Allele.create("A", true), Arrays.asList(Allele.create("C", false), Allele.create("T", false)), 1, 0,
                0, false
        }, {
                Allele.create("AA", true), Arrays.asList(Allele.create("CC", false), Allele.create("T", false)), 1, 0,
                0, false
        }, {
                Allele.create("AA", true), Arrays.asList(Allele.create("T", false), Allele.create("CC", false)), 1, 0,
                1, false
        }, {
                Allele.create("GA", true), Arrays.asList(Allele.create("T", false), Allele.create("CC", false)), 1, 0,
                1, true
        }
        };
    }

    @Test(dataProvider = "testInsertions")
    public void testChooseBestAlleleInsertion(final int offsetIntoRead, final String refBases, final String altBases,
                                              final boolean isSpliceInAlt) {

        // We assume that we are on the base right before the insertion.  We are pretending each read is 50 bases long.
        final Allele referenceAllele = Allele.create(refBases, true);
        final Allele insertionAllele = Allele.create(altBases, false);

        if (isSpliceInAlt) {
            final PileupElement pileupElement = ArtificialReadUtils.createSplicedInsertionPileupElement(offsetIntoRead, insertionAllele, GATKVariantContextUtilsUnitTest.LENGTH_OF_ARTIFICIAL_READ);
            Assert.assertEquals(
                    GATKVariantContextUtils.chooseAlleleForRead(pileupElement, referenceAllele,
                            Collections.singletonList(insertionAllele), 0),
                    insertionAllele);
        } else {
            final PileupElement pileupElement = ArtificialReadUtils.createNonIndelPileupElement(offsetIntoRead, referenceAllele, GATKVariantContextUtilsUnitTest.LENGTH_OF_ARTIFICIAL_READ);
            Assert.assertEquals(
                    GATKVariantContextUtils.chooseAlleleForRead(pileupElement, referenceAllele,
                            Collections.singletonList(insertionAllele), 0),
                    referenceAllele);
        }
    }

    @Test(dataProvider = "testDeletions")
    public void testChooseBestAlleleDeletion(final int offsetIntoRead, final String refBases, final String altBases,
                                             final boolean isSpliceInAlt) {

        // We assume that we are on the base right before the insertion.  We are pretending each read is 50 bases long.
        final Allele referenceAllele = Allele.create(refBases, true);
        final Allele deletionAllele = Allele.create(altBases, false);

        if (isSpliceInAlt) {
            final PileupElement pileupElement = ArtificialReadUtils.createSplicedDeletionPileupElement(offsetIntoRead, referenceAllele, GATKVariantContextUtilsUnitTest.LENGTH_OF_ARTIFICIAL_READ);
            Assert.assertEquals(
                    GATKVariantContextUtils.chooseAlleleForRead(pileupElement, referenceAllele,
                            Collections.singletonList(deletionAllele), 0),
                    deletionAllele);
        } else {
            final PileupElement pileupElement = ArtificialReadUtils.createNonIndelPileupElement(offsetIntoRead, referenceAllele, GATKVariantContextUtilsUnitTest.LENGTH_OF_ARTIFICIAL_READ);
            Assert.assertEquals(
                    GATKVariantContextUtils.chooseAlleleForRead(pileupElement, referenceAllele,
                            Collections.singletonList(deletionAllele), 0),
                    referenceAllele);
        }
    }

    @Test
    public void testChooseBestAlleleNull() {

        final int offsetIntoRead = 10;

        // We assume that we are on the base right before the insertion.  We are pretending each read is 50 bases long.
        final Allele referenceAllele = Allele.create("T", true);
        final Allele deletionAllele = Allele.create("A", false);

        final PileupElement pileupElement = ArtificialReadUtils.createNonIndelPileupElement(offsetIntoRead, referenceAllele, GATKVariantContextUtilsUnitTest.LENGTH_OF_ARTIFICIAL_READ);
        final byte[] newBases = pileupElement.getRead().getBases();
        newBases[offsetIntoRead] = (byte) 'C';
        pileupElement.getRead().setBases(newBases);

        Assert.assertEquals(
                GATKVariantContextUtils.chooseAlleleForRead(pileupElement, referenceAllele,
                        Collections.singletonList(deletionAllele), 0),
                null);
    }

    @DataProvider(name="testInsertions")
    public Object[][] createTestInsertions() {
        return new Object[][] {
                {0, "A", "ATT", true},
                {0, "A", "ATT", false},
                {10, "A", "ATT", true},
                {10, "A", "ATT", false}
        };
    }
    @DataProvider(name="testDeletions")
    public Object[][] createTestDeletions() {
        return new Object[][] {
                {0, "ATT", "A", true},
                {0, "ATT", "A", false},
                {10, "ATT", "A", true},
                {10, "ATT", "A", false}
        };
    }

    @Test(dataProvider = "doesReadContainAllele")
    public void testDoesReadContainAllele(final int offsetIntoRead, final String actualBases, final String searchAlleleBases,
                                          final Trilean gt) {

        final PileupElement pileupElement = ArtificialReadUtils.createNonIndelPileupElement(offsetIntoRead, Allele.create(actualBases, true), GATKVariantContextUtilsUnitTest.LENGTH_OF_ARTIFICIAL_READ);
        Assert.assertEquals(GATKVariantContextUtils.doesReadContainAllele(pileupElement, Allele.create(searchAlleleBases, false)),
                gt);
        Assert.assertEquals(GATKVariantContextUtils.doesReadContainAllele(pileupElement, Allele.create(actualBases, false)),
                Trilean.TRUE);
    }

    @DataProvider(name="doesReadContainAllele")
    public Object[][] createDoesReadContainAlelle() {
        return new Object[][] {
                {10, "ATT", "C", Trilean.FALSE},
                {10, "AT", "AT", Trilean.TRUE},
                {49, "A", "ATT", Trilean.UNKNOWN},
                {48, "AT", "ATT", Trilean.UNKNOWN},
        };
    }

    @DataProvider(name = "LeftAlignDataProvider")
    public Object[][] LeftAlignTestData() throws UnsupportedEncodingException {
        List<Object[]> tests = new ArrayList<Object[]>();
        for (int iLongPiece = 0; iLongPiece < longPieces.size(); iLongPiece++) {
            //longPiece is repeat of str
            final String longPiece = longPieces.get(iLongPiece);
            final ReferenceMemorySource refSource = refSources.get(iLongPiece);
            final int strLength = strLengths.get(iLongPiece);
            //final int nStrs=longPiece.length()/strLength; //number of strs in the repeat
            final byte[] theseRefBases = refBases.get(iLongPiece).getBases();
            Integer indelIndexStart, indelIndexStop;
            if (longPiece.length() < 90) {
                //for shorter repeats, we will place the indel at each location from just before the repeat to 20 bases into the repeat
                indelIndexStart = repeatStart - 1;
                indelIndexStop = repeatStart + Math.min(longPiece.length(), 20);
            } else {
                //for longer repeats, we will place the indel at the 61st base of the repeat, to test that automatic expansion works
                indelIndexStart = repeatStart + 60;
                indelIndexStop = indelIndexStart + 1;
            }

            for (int indelIndex = indelIndexStart; indelIndex < indelIndexStop; indelIndex++) {
                //we will delete up to 50 units of the str (if possible), and insert up to 50 units of the str
                for (int indelRepeats = Math.max(-(longPiece.length() - (indelIndex - repeatStart + 1)), -50) / strLength; indelRepeats < 50; indelRepeats++) {
                    if (indelRepeats == 0) {
                        continue;
                    }
                    final List<Allele> alleles = new ArrayList<Allele>();
                    if (indelRepeats < 0) { // deletion
                        byte[] basesRef = new byte[Math.abs(indelRepeats) * strLength + 1];
                        byte[] basesAlt = new byte[1];
                        System.arraycopy(theseRefBases, indelIndex, basesRef, 0, Math.abs(indelRepeats) * strLength + 1);
                        System.arraycopy(theseRefBases, indelIndex, basesAlt, 0, 1);
                        alleles.add(Allele.create(new String(basesRef, "UTF-8"), true));
                        alleles.add(Allele.create(new String(basesAlt, "UTF-8"), false));
                    } else { //insertion
                        byte[] basesRef = new byte[1];
                        byte[] basesAlt = new byte[Math.abs(indelRepeats) * strLength + 1];
                        byte[] basesRepeat = new byte[strLength];
                        System.arraycopy(theseRefBases, indelIndex, basesRef, 0, 1);
                        System.arraycopy(theseRefBases, indelIndex, basesAlt, 0, 1);
                        if (indelIndex < repeatStart + longPiece.length() - strLength) {
                            //look at next strLength bases to find str
                            System.arraycopy(theseRefBases, indelIndex + 1, basesRepeat, 0, strLength);
                        } else {
                            //look at previous strLength bases (including this one) to find str
                            System.arraycopy(theseRefBases, indelIndex - strLength + 1, basesRepeat, 0, strLength);
                        }
                        System.arraycopy(Utils.dupString(new String(basesRepeat, "UTF-8"), indelRepeats).getBytes(), 0, basesAlt, 1, Math.abs(indelRepeats) * strLength);
                        alleles.add(Allele.create(new String(basesRef, "UTF-8"), true));
                        alleles.add(Allele.create(new String(basesAlt, "UTF-8"), false));
                    }
                    final SimpleInterval interval = new SimpleInterval(artificialContig, indelIndex + 1, indelIndex + alleles.get(0).length());
                    ReferenceContext ref = new ReferenceContext(refSource, interval);

                    final VariantContext vc = new VariantContextBuilder("test", artificialContig, indelIndex + 1, indelIndex + alleles.get(0).length(), alleles).make();
                    final boolean expectRealigned = (indelIndex != repeatStart - 1);
                    final int expectedStart = expectRealigned ? repeatStart : indelIndex + 1;
                    tests.add(new Object[]{vc, ref, expectRealigned, expectedStart});
                }
            }
            //also add an insertion of 1 str one base past end of repeat.  This insertion should not be moved
            byte[] basesRef = new byte[1];
            byte[] basesAlt = new byte[strLength + 1];
            byte[] basesRepeat = new byte[strLength];
            final int indelIndex = repeatStart + longPiece.length();
            System.arraycopy(theseRefBases, indelIndex, basesRef, 0, 1);
            System.arraycopy(theseRefBases, indelIndex, basesAlt, 0, 1);
            // insert first strLength bases of repeat
            System.arraycopy(theseRefBases, repeatStart, basesRepeat, 0, strLength);
            System.arraycopy(Utils.dupString(new String(basesRepeat, "UTF-8"), 1).getBytes(), 0, basesAlt, 1, strLength);
            final List<Allele> alleles = new ArrayList<Allele>();
            alleles.add(Allele.create(new String(basesRef, "UTF-8"), true));
            alleles.add(Allele.create(new String(basesAlt, "UTF-8"), false));
            final SimpleInterval interval = new SimpleInterval(artificialContig, indelIndex + 1, indelIndex + alleles.get(0).length());
            ReferenceContext ref = new ReferenceContext(refSource, interval);

            final VariantContext vc = new VariantContextBuilder("test", artificialContig, indelIndex + 1, indelIndex + alleles.get(0).length(), alleles).make();
            tests.add(new Object[]{vc, ref, false, indelIndex + 1});

        }
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "LeftAlignDataProvider")
    public void testLeftAlign(final VariantContext vc, final ReferenceContext ref, final boolean expectRealigned, final int expectedStart) {
        final VariantContext realignedV = GATKVariantContextUtils.leftAlignAndTrim(vc, ref, LeftAlignAndTrimVariants.DEFAULT_MAX_LEADING_BASES, true);
        Assert.assertEquals(realignedV != vc, expectRealigned);
        Assert.assertEquals(realignedV.getStart(), expectedStart);
    }

    @Test
    public void testLeftAlignWithVaryingMaxDistances() {

        final byte[] refSequence = new byte[200];
        Arrays.fill(refSequence, 0, 100, (byte) 'C');
        Arrays.fill(refSequence, 100, 200, (byte) 'A');

        final String contig = "1";
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, refSequence.length);
        final SimpleInterval interval = new SimpleInterval(contig, 1, refSequence.length);


        final ReferenceBases refBases = new ReferenceBases(refSequence, interval);
        final ReferenceDataSource referenceDataSource = new ReferenceMemorySource(refBases, header.getSequenceDictionary());
        final ReferenceContext referenceContext = new ReferenceContext(referenceDataSource, interval);

        final List<Allele> alleles = Arrays.asList(Allele.create("AA", true), Allele.create("A", false));
        for (int maxShift : new int[] {0, 1, 5, 10, 30}) {
            for (int start : new int[]{101, 105, 109, 110, 111, 115, 120, 130, 141}) {
                final VariantContext vc = new VariantContextBuilder("source", contig, start, start + 1, alleles).make();
                final VariantContext realignedV = GATKVariantContextUtils.leftAlignAndTrim(vc, referenceContext, maxShift, true);

                Assert.assertEquals(realignedV.getStart(), Math.max(start - maxShift, 100));
            }
        }
    }
}
