package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.CoordMath;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVCollapser;
import org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVLinkage;
import org.broadinstitute.hellbender.tools.sv.cluster.ClusteringParameters;
import org.broadinstitute.hellbender.tools.sv.cluster.SVClusterEngine;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.testng.Assert;
import org.testng.TestException;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.CPX_SV_SYB_ALT_ALLELE_STR;

public class SVTestUtils {

    public static final ReferenceSequenceFile hg38Reference = ReferenceUtils.createReferenceReader(new GATKPath(GATKBaseTest.hg38Reference));
    public final static SAMSequenceDictionary hg38Dict = hg38Reference.getSequenceDictionary();
    public final static int chr1Length = hg38Dict.getSequence("chr1").getSequenceLength();

    private final static GenomeLocParser glParser = new GenomeLocParser(SVTestUtils.hg38Dict);
    public static final CanonicalSVCollapser defaultCollapser =
            new CanonicalSVCollapser(
                    hg38Reference,
                    CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE,
                    CanonicalSVCollapser.BreakpointSummaryStrategy.MEDIAN_START_MEDIAN_END,
                    CanonicalSVCollapser.FlagFieldLogic.OR);

    public static final String PESR_ALGORITHM = "pesr";
    public static final List<String> DEPTH_ONLY_ALGORITHM_LIST = Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM);
    public static final List<String> PESR_ONLY_ALGORITHM_LIST = Collections.singletonList(PESR_ALGORITHM);
    public static final Allele CPX_ALLELE = Allele.create(SimpleSVType.createBracketedSymbAlleleString(CPX_SV_SYB_ALT_ALLELE_STR));
    public static final Allele BND_ALLELE = Allele.create(SimpleSVType.createBracketedSymbAlleleString(GATKSVVCFConstants.BREAKEND_STR));

    public static CanonicalSVLinkage<SVCallRecord> getNewDefaultLinkage() {
        final CanonicalSVLinkage<SVCallRecord> linkage = new CanonicalSVLinkage<>(SVTestUtils.hg38Dict, false);
        linkage.setDepthOnlyParams(defaultDepthOnlyParameters);
        linkage.setMixedParams(defaultMixedParameters);
        linkage.setEvidenceParams(defaultEvidenceParameters);
        return linkage;
    }

    public static SVClusterEngine getNewDefaultSingleLinkageEngine() {
        return new SVClusterEngine(SVClusterEngine.CLUSTERING_TYPE.SINGLE_LINKAGE, defaultCollapser::collapse, getNewDefaultLinkage(), hg38Dict);
    }

    public static SVClusterEngine getNewDefaultMaxCliqueEngine() {
        return new SVClusterEngine(SVClusterEngine.CLUSTERING_TYPE.MAX_CLIQUE, defaultCollapser::collapse, getNewDefaultLinkage(), hg38Dict);
    }

    public static final ClusteringParameters defaultDepthOnlyParameters = ClusteringParameters.createDepthParameters(0.8, 0, 10000000, 0);
    public static final ClusteringParameters defaultMixedParameters = ClusteringParameters.createMixedParameters(0.8, 0, 1000, 0);
    public static final ClusteringParameters defaultEvidenceParameters = ClusteringParameters.createPesrParameters(0.5, 0, 500, 0);

    public static final SVClusterEngine defaultSingleLinkageEngine = getNewDefaultSingleLinkageEngine();
    public static final SVClusterEngine defaultMaxCliqueEngine = getNewDefaultMaxCliqueEngine();

    public final static int start = 10001;

    public final static int length = 10000;

    public final static int length1b = (int)Math.round(defaultDepthOnlyParameters.getReciprocalOverlap() * length);

    public final static int start1b = start + length - length1b;

    //separated from end of call1 by defragmenter padding (in bin space, according to bins defined by targetIntervals below)
    public final static int start2 = start + length*14/9;

    public final static int start3 = start + (int)Math.round(Math.floor((1- defaultDepthOnlyParameters.getReciprocalOverlap())*length));

    //make intervals like xxxx----xxxx----xxxx----xxxx----xxxx
    public final static List<GenomeLoc> targetIntervals = new ArrayList<>(
            Arrays.asList(glParser.createGenomeLoc("chr1", 1, length/2),  //left edge call
                    glParser.createGenomeLoc("chr1", start, start + length/9),  //start of call1
                    glParser.createGenomeLoc("chr1", start + length*2/9, start + length*3/9),
                    glParser.createGenomeLoc("chr1", start + length*4/9, start + length*5/9),
                    glParser.createGenomeLoc("chr1", start + length*6/9, start + length*7/9),
                    glParser.createGenomeLoc("chr1", start + length*8/9, start + length),  //end of call1
                    glParser.createGenomeLoc("chr1", start + length*10/9, start + length*11/9),  //padding for call1
                    glParser.createGenomeLoc("chr1", start + length*12/9, start + length*13/9),  //padding for call1
                    glParser.createGenomeLoc("chr1", start2, start2 + length/9),
                    glParser.createGenomeLoc("chr1", start2 + length*2/9, start2 + length*3/9),
                    glParser.createGenomeLoc("chr1", start2 + length*4/9, start2 + length*5/9),
                    glParser.createGenomeLoc("chr1", start2 + length*6/9, start2 + length*7/9),
                    glParser.createGenomeLoc("chr1", start2 + length*8/9, start2 + length),
                    glParser.createGenomeLoc("chr1", start2 + length, start2 + length + length/9),
                    glParser.createGenomeLoc("chr1", start2 + length+ length*2/9, start2 + length+ length*3/9),
                    glParser.createGenomeLoc("chr1", start2 + length+ length*4/9, start2 + length+ length*5/9),
                    glParser.createGenomeLoc("chr1", start2 + length+ length*6/9, start2 + length+ length*7/9),
                    glParser.createGenomeLoc("chr1", start2 + length+ length*8/9, start2 + length+ length),
                    glParser.createGenomeLoc("chr1", chr1Length - 99, chr1Length),
                    glParser.createGenomeLoc("chr10", start, start + length - 1)))
            ;

    public final static GenotypeBuilder sample1 = new GenotypeBuilder("sample1",
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL))
            .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 1);

    public final static GenotypeBuilder sample1Ploidy1 = new GenotypeBuilder("sample1",
            Lists.newArrayList(Allele.SV_SIMPLE_DEL))
            .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 1);

    public final static GenotypeBuilder sample1_CN0 = new GenotypeBuilder("sample1",
            Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL))
            .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 0);

    public final static GenotypeBuilder sample2 = new GenotypeBuilder("sample2",
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP))
            .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 3);

    public final static GenotypeBuilder sample3 = new GenotypeBuilder("sample3",
            Lists.newArrayList(Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DEL))
            .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 0);

    public final static SVCallRecord makeRecord(final String id,
                                                final String contigA,
                                                final int positionA,
                                                final Boolean strandA,
                                                final String contigB,
                                                final int positionB,
                                                final Boolean strandB,
                                                final GATKSVVCFConstants.StructuralVariantAnnotationType type,
                                                final Integer length,
                                                final List<String> algorithms,
                                                final List<Allele> alleles,
                                                final List<GenotypeBuilder> genotypeBuilders) {
        return makeRecordWithEvidenceAndQuality(id, contigA, positionA, strandA, contigB, positionB, strandB, type, length, Collections.emptyList(), algorithms, alleles, genotypeBuilders, null);
    }

    public final static SVCallRecord makeRecordWithEvidenceAndQuality(final String id,
                                                 final String contigA,
                                                 final int positionA,
                                                 final Boolean strandA,
                                                 final String contigB,
                                                 final int positionB,
                                                 final Boolean strandB,
                                                 final GATKSVVCFConstants.StructuralVariantAnnotationType type,
                                                 final Integer length,
                                                 final List<GATKSVVCFConstants.EvidenceTypes> evidence,
                                                 final List<String> algorithms,
                                                 final List<Allele> alleles,
                                                 final List<GenotypeBuilder> genotypeBuilders,
                                                 final Double log10PError) {
        final Allele refAllele = Allele.create(ReferenceUtils.getRefBaseAtPosition(hg38Reference, contigA, positionA), true);
        final List<Allele> newAlleles = replaceRefAlleles(alleles, refAllele);
        final List<Genotype> genotypes = new ArrayList<>(genotypeBuilders.size());
        for (final GenotypeBuilder builder : genotypeBuilders) {
            genotypes.add(makeGenotypeWithRefAllele(builder, refAllele));
        }
        return new SVCallRecord(id, contigA, positionA, strandA, contigB, positionB, strandB, type, null, Collections.emptyList(), length, evidence, algorithms,
                newAlleles, genotypes, Collections.emptyMap(), Collections.emptySet(), log10PError, hg38Dict);
    }

    public static final Genotype makeGenotypeWithRefAllele(final GenotypeBuilder builder, final Allele refAllele) {
        final List<Allele> alleles = replaceRefAlleles(builder.make().getAlleles(), refAllele);
        builder.alleles(alleles);
        builder.attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, alleles.size());
        return builder.make();
    }

    public static final SVCallRecord makeRecordWithCarriers(final List<String> allSamples, final Set<String> carrierSamples) {
        final List<GenotypeBuilder> genotypes = makeDeletionGenotypesWithCarriers(allSamples, carrierSamples, Allele.SV_SIMPLE_DEL);
        return makeRecord("", "chr1", 100, Boolean.TRUE, "chr1", 200, Boolean.FALSE,
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                List.of(Allele.REF_N, Allele.SV_SIMPLE_DEL), genotypes);
    }

    public static final SVCallRecord makeComplexRecordWithCarriers(final List<String> allSamples, final Set<String> carrierSamples) {
        final List<Genotype> genotypes = makeDeletionGenotypesWithCarriers(allSamples, carrierSamples, SVTestUtils.CPX_ALLELE)
                .stream().map(GenotypeBuilder::make).collect(Collectors.toUnmodifiableList());
        return new SVCallRecord("cpx1", "chr1", 1000, null,
                "chr1", 2000, null,
                GATKSVVCFConstants.StructuralVariantAnnotationType.CPX,
                GATKSVVCFConstants.ComplexVariantSubtype.delINV,
                Arrays.asList(SVCallRecord.ComplexEventInterval.decode("DEL_chr1:1100-1500", SVTestUtils.hg38Dict), SVCallRecord.ComplexEventInterval.decode("INV_chr1:1600-1900", SVTestUtils.hg38Dict)),
                1000, Collections.emptyList(), Collections.singletonList(SVTestUtils.PESR_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, SVTestUtils.CPX_ALLELE),
                genotypes, Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
    }

    private static List<GenotypeBuilder> makeDeletionGenotypesWithCarriers(final List<String> allSamples, final Set<String> carrierSamples, final Allele altAllele) {
        final List<GenotypeBuilder> genotypes = new ArrayList<>(allSamples.size());
        for (final String sample : allSamples) {
            final GenotypeBuilder builder = new GenotypeBuilder(sample);
            if (carrierSamples.contains(sample)) {
                builder.alleles(List.of(Allele.REF_N, altAllele));
            } else {
                builder.alleles(List.of(Allele.REF_N, Allele.REF_N));
            }
            builder.attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2);
            genotypes.add(builder);
        }
        return genotypes;
    }

    public final static SVCallRecord rightEdgeCall = makeRecord("rightEdgeCall", "chr1", chr1Length - 99, null,
            "chr1", chr1Length, null,
            GATKSVVCFConstants.StructuralVariantAnnotationType.CNV, 100,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
            Arrays.asList(sample1, sample2));

    public final static SVCallRecord leftEdgeCall = makeRecord("leftEdgeCall", "chr1", 1, null,
            "chr1", length/2, null,
            GATKSVVCFConstants.StructuralVariantAnnotationType.CNV, length/2,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
            Arrays.asList(sample1, sample2));

    public final static SVCallRecord nonDepthOnly = makeRecord("nonDepthOnly", "chr1", start, null,
            "chr1", start + length -1, null,
            GATKSVVCFConstants.StructuralVariantAnnotationType.CNV, length,
            Arrays.asList(GATKSVVCFConstants.DEPTH_ALGORITHM, "PE"),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
            Arrays.asList(sample1, sample2));

    public final static SVCallRecord call1 = makeRecord("call1", "chr1", start, null,
            "chr1", start + length -1, null,
            GATKSVVCFConstants.StructuralVariantAnnotationType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
            Arrays.asList(sample1, sample2));

    public final static SVCallRecord call1b = makeRecord("call1b", "chr1", start1b, null,
            "chr1", start1b + length1b - 1, null,
            GATKSVVCFConstants.StructuralVariantAnnotationType.CNV, length1b,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
            Arrays.asList(sample1, sample2));

    public final static SVCallRecord call2 = makeRecord("call2", "chr1", start2, null,
            "chr1", start2 + length - 1, null,
            GATKSVVCFConstants.StructuralVariantAnnotationType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
            Arrays.asList(sample1, sample2));

    public final static SVCallRecord call3 = makeRecord("call3", "chr1", start2+length, null,
            "chr1", start2 + length + length - 1, null,
            GATKSVVCFConstants.StructuralVariantAnnotationType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
            Arrays.asList(sample1, sample2));

    public final static SVCallRecord call4_chr10 = makeRecord("call4_chr10", "chr10", start, null,
            "chr10", start + length - 1, null,
            GATKSVVCFConstants.StructuralVariantAnnotationType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
            Arrays.asList(sample1, sample2));

    public final static SVCallRecord call1_CN1 = makeRecord("call1_CN1", "chr1", start, null,
            "chr1", start + length - 1, null,
            GATKSVVCFConstants.StructuralVariantAnnotationType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
            Arrays.asList(sample1));

    public final static SVCallRecord call2_CN0 = makeRecord("call2_CN0", "chr1", start2, null,
            "chr1", start2 + length - 1, null,
            GATKSVVCFConstants.StructuralVariantAnnotationType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
            Arrays.asList(sample1_CN0));

    public final static SVCallRecord sameBoundsSampleMismatch = makeRecord("sameBoundsSampleMismatch", "chr1", start, null,
            "chr1", start + length - 1, null,
            GATKSVVCFConstants.StructuralVariantAnnotationType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
            Arrays.asList(sample1, sample3));

    public final static List<Genotype> threeGenotypes = Arrays.asList(sample1.make(), sample2.make(), sample3.make());

    public final static SVCallRecord inversion = makeRecord("inversion", "chr1", 10000, true, "chr1", 20000, true,
            GATKSVVCFConstants.StructuralVariantAnnotationType.INV, 10001, Arrays.asList("SR", "PE"), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP), Collections.singletonList(sample2));

    public static final SVCallRecord depthOnly = makeRecord("depthOnly", "chr1", 10000, null, "chr1", 20000, null,
            GATKSVVCFConstants.StructuralVariantAnnotationType.CNV, 10001, Collections.singletonList("depth"), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL), Collections.singletonList(sample1));

    public static final SVCallRecord depthAndStuff = makeRecord("depthAndStuff", "chr1", 10000, null, "chr1", 20000, null,
                    GATKSVVCFConstants.StructuralVariantAnnotationType.CNV, 10001, Arrays.asList("depth", "PE"), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP), Collections.singletonList(sample2));

    public static final SVCallRecord depthAndStuff2 = makeRecord("depthAndStuff2", "chr1", 10000 - defaultEvidenceParameters.getWindow(), null, "chr1", 20000, null,
            GATKSVVCFConstants.StructuralVariantAnnotationType.CNV, 10001 + defaultEvidenceParameters.getWindow(), Arrays.asList("depth", "PE"), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP), Collections.singletonList(sample2));

    public static final SVCallRecord depthAndStuff3 = makeRecord("depthAndStuff3", "chr1", 10000 + defaultEvidenceParameters.getWindow(), null, "chr1", 20000, null,
            GATKSVVCFConstants.StructuralVariantAnnotationType.CNV, 10001 - defaultEvidenceParameters.getWindow(), Arrays.asList("depth", "PE"), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP), Collections.singletonList(sample2));

    public final static SVCallRecord overlapsCall1 = makeRecord("overlapsCall1", "chr1", start3, null,
            "chr1", start3 + length -1, null,
            GATKSVVCFConstants.StructuralVariantAnnotationType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
            Arrays.asList(sample1, sample2));

    public static void assertEqualsExceptMembership(final SVCallRecord one, final SVCallRecord two) {
        assertEqualsExceptExcludedAttributes(one, two, Collections.singletonList(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY));
    }

    public static void assertEqualsExceptMembershipAndGT(final SVCallRecord one, final SVCallRecord two) {
        assertEqualsExceptExcludedAttributes(one, two, Collections.singletonList(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY), true);
    }

    private static Genotype nullGT(final Genotype gt) {
        return new GenotypeBuilder(gt).alleles(null).make();
    }

    public static void assertEqualsExceptExcludedAttributes(final SVCallRecord one, final SVCallRecord two,
                                                            final Collection<String> excludedAttributeKeys) {
        assertEqualsExceptExcludedAttributes(one, two, excludedAttributeKeys, false);
    }

    public static void assertEqualsExceptExcludedAttributes(final SVCallRecord one, final SVCallRecord two,
                                                            final Collection<String> excludedAttributeKeys,
                                                            final boolean ignoreGT) {
        if (one == two) return;
        Assert.assertEquals(one.getAlgorithms(), two.getAlgorithms());
        Assert.assertEquals(one.getAlleles(), two.getAlleles());

        // Exclude MEMBERS field
        final Map<String, Object> attributesOne = new HashMap<>(one.getAttributes());
        final Map<String, Object> attributesTwo = new HashMap<>(two.getAttributes());
        for (final String key : excludedAttributeKeys) {
            attributesOne.remove(key);
            attributesTwo.remove(key);
        }
        Assert.assertEquals(attributesOne, attributesTwo);

        Assert.assertEquals(one.getPositionAInterval(), two.getPositionAInterval());
        Assert.assertEquals(one.getPositionBInterval(), two.getPositionBInterval());
        Assert.assertEquals(one.getStrandA(), two.getStrandA());
        Assert.assertEquals(one.getStrandB(), two.getStrandB());
        Assert.assertEquals(one.getId(), two.getId());
        Assert.assertEquals(one.getLength(), two.getLength());
        Assert.assertEquals(one.getType(), two.getType());
        Assert.assertEquals(one.getComplexSubtype(), two.getComplexSubtype());
        Assert.assertEquals(one.getComplexEventIntervals(), two.getComplexEventIntervals());
        Assert.assertEquals(one.getGenotypes().size(), two.getGenotypes().size());
        for (int i = 0; i < one.getGenotypes().size(); i++) {
            final Genotype g1 = ignoreGT ? nullGT(one.getGenotypes().get(i)) : one.getGenotypes().get(i);
            final Genotype g2 = ignoreGT ? nullGT(two.getGenotypes().get(i)) : two.getGenotypes().get(i);
            VariantContextTestUtils.assertGenotypesAreEqual(g1, g2);
        }
    }

    public static void assertContainsAllIgnoreRefAlleleBase(final GenotypesContext one, final GenotypesContext two,
                                                            final boolean ignoreGT) {
        Assert.assertTrue(one.size() >= two.size());
        if ( !one.isEmpty()) {
            Assert.assertTrue(one.getSampleNames().containsAll(two.getSampleNames()), "sample names set");
            final Set<String> samples = two.getSampleNames();
            for ( final String sample : samples ) {
                assertSameGenotypeExceptRefBase(one.get(sample), two.get(sample), ignoreGT);
            }
        }
    }

    public static void assertSameGenotypeExceptRefBase(final Genotype g1, final Genotype g2, final boolean ignoreGT) {
        final GenotypeBuilder builderOne = new GenotypeBuilder(g1);
        final GenotypeBuilder builderTwo = new GenotypeBuilder(g2);
        if (ignoreGT) {
            builderOne.alleles(null);
            builderTwo.alleles(null);
        } else {
            builderOne.alleles(replaceRefAlleles(g1.getAlleles(), Allele.REF_N));
            builderTwo.alleles(replaceRefAlleles(g2.getAlleles(), Allele.REF_N));
        }
        VariantContextTestUtils.assertGenotypesAreEqual(builderOne.make(), builderTwo.make());
    }

    public static List<Allele> replaceRefAlleles(final List<Allele> alleles, final Allele replace) {
        return alleles.stream().map(a -> a.isReference() ? replace : a).collect(Collectors.toList());
    }

    public static Genotype buildHomGenotypeWithPloidy(final Allele allele, final int ploidy) {
        return new GenotypeBuilder()
                .alleles(buildHomAlleleListWithPloidy(allele, ploidy))
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, ploidy)
                .make();
    }

    public static List<Allele> buildHomAlleleListWithPloidy(final Allele allele, final int ploidy) {
        return Collections.nCopies(ploidy, allele);
    }

    public static List<Allele> buildBiallelicListWithPloidy(final Allele altAllele, final Allele refAllele, final int copyNumber, final int ploidy) {
        final int numAlt;
        if (altAllele.equals(Allele.SV_SIMPLE_DEL)) {
            numAlt = Math.max(ploidy - copyNumber, 0);
        } else if (altAllele.equals(Allele.SV_SIMPLE_DUP)) {
            numAlt = Math.max(copyNumber - ploidy, 0);
        } else {
            throw new TestException("Unsupported alt allele: " + altAllele.getDisplayString());
        }
        final List<Allele> alleles = new ArrayList<>(ploidy);
        final int numRef = ploidy - numAlt;
        for (int i = 0; i < numRef; i++) {
            alleles.add(refAllele);
        }
        for (int i = 0; i < numAlt; i++) {
            alleles.add(altAllele);
        }
        return alleles;
    }

    public static SVCallRecord newCallRecordWithAlleles(final List<Allele> genotypeAlleles, final List<Allele> variantAlleles,
                                                        final GATKSVVCFConstants.StructuralVariantAnnotationType svtype, final Integer expectedCopyNumber,
                                                        final Integer copyNumber) {
        return newCallRecordWithAllelesAndSampleName("sample", genotypeAlleles, variantAlleles, svtype, expectedCopyNumber, copyNumber);
    }

    public static SVCallRecord newCallRecordWithAllelesAndSampleName(final String sample,
                                                                     final List<Allele> genotypeAlleles, final List<Allele> variantAlleles,
                                                                     final GATKSVVCFConstants.StructuralVariantAnnotationType svtype, final Integer expectedCopyNumber,
                                                                     final Integer copyNumber) {
        GenotypeBuilder builder = new GenotypeBuilder(sample).alleles(genotypeAlleles);
        if (expectedCopyNumber != null) {
            builder = builder.attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, expectedCopyNumber);
        }
        if (copyNumber != null) {
            builder = builder.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, copyNumber);
        }
        return new SVCallRecord("", "chr1", 100, getValidTestStrandA(svtype), "chr1", 199, getValidTestStrandB(svtype),
                svtype, null, Collections.emptyList(), 100, Collections.emptyList(), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                variantAlleles,
                Collections.singletonList(builder.make()),
                Collections.emptyMap(), Collections.emptySet(), null);
    }

    public static SVCallRecord newDeletionRecordWithCoordinates(final String id, final String contig, final int positionA, final int positionB) {
        return new SVCallRecord(id, contig, positionA, true, contig, positionB, false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                getLength(positionA, positionB, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL),
                Collections.emptyList(), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Collections.emptyList(),
                Collections.emptyList(),
                Collections.emptyMap(), Collections.emptySet(), null);
    }

    public static SVCallRecord newDeletionRecordWithAttributes(final Map<String, Object> attributes) {
        return new SVCallRecord("", "chr1", 100, true, "chr1", 199, false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                100, Collections.emptyList(), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Collections.emptyList(),
                Collections.emptyList(),
                attributes, Collections.emptySet(), null);
    }

    public static SVCallRecord newNamedDeletionRecordWithAttributesAndGenotypes(final String id,
                                                                                final List<Genotype> genotypes,
                                                                                final Map<String, Object> attributes) {
        return new SVCallRecord(id, "chr1", 100, true, "chr1", 199, false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                100, Collections.emptyList(), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                genotypes,
                attributes, Collections.emptySet(), null);
    }

    public static final Map<String, Object> keyValueArraysToMap(final String[] keys, final Object[] values) {
        final Map<String, Object> map = new HashMap<>();
        for (int i = 0; i < keys.length; i++) {
            map.put(keys[i], values[i]);
        }
        return map;
    }

    // Note that strands may not be valid
    public static SVCallRecord newCallRecordWithLengthAndType(final Integer length, final GATKSVVCFConstants.StructuralVariantAnnotationType svtype) {
        final int positionB = length == null ? 1 : CoordMath.getEnd(1, length);
        return new SVCallRecord("", "chr1", 1, getValidTestStrandA(svtype), "chr1", positionB, getValidTestStrandB(svtype),
                svtype, null, Collections.emptyList(), length, Collections.emptyList(), PESR_ONLY_ALGORITHM_LIST, Collections.emptyList(), Collections.emptyList(),
                Collections.emptyMap(), Collections.emptySet(), null);
    }

    public static SVCallRecord newDeletionCallRecordWithIdAndAlgorithms(final String id, final List<String> algorithms) {
        return new SVCallRecord(id, "chr1", 1, true, "chr1", 100, false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(), 100, Collections.emptyList(), algorithms, Collections.emptyList(),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null);
    }

    // Note strands and length may not be set properly
    public static SVCallRecord newPESRCallRecordWithIntervalAndType(final int start, final int end, final GATKSVVCFConstants.StructuralVariantAnnotationType svtype) {
        return new SVCallRecord("", "chr1", start, getValidTestStrandA(svtype), "chr1", end, getValidTestStrandB(svtype),
                svtype, null, Collections.emptyList(), getLength(start, end, svtype), Collections.emptyList(), PESR_ONLY_ALGORITHM_LIST, Collections.emptyList(),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null);
    }

    // Note strands and length may not be set properly
    public static SVCallRecord newInsertionWithPositionAndLength(final int start, final int length) {
        return new SVCallRecord("", "chr1", start, true, "chr1", start + 1, false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.INS, null, Collections.emptyList(), length, Collections.emptyList(), PESR_ONLY_ALGORITHM_LIST, Collections.emptyList(),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null);
    }

    public static SVCallRecord newDepthCallRecordWithIntervalAndType(final int start, final int end, final GATKSVVCFConstants.StructuralVariantAnnotationType svtype) {
        return new SVCallRecord("", "chr1", start, getValidTestStrandA(svtype), "chr1", end, getValidTestStrandB(svtype),
                svtype, null, Collections.emptyList(), getLength(start, end, svtype), Collections.emptyList(), DEPTH_ONLY_ALGORITHM_LIST, Collections.emptyList(),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null);
    }

    // Note strands and length may not be set properly
    public static SVCallRecord newCallRecordWithContigsIntervalAndType(final String startContig, final int start, final String endContig, final int end, final GATKSVVCFConstants.StructuralVariantAnnotationType svtype) {
        return new SVCallRecord("", startContig, start, getValidTestStrandA(svtype), endContig, end, getValidTestStrandB(svtype),
                svtype, null, Collections.emptyList(), getLength(start, end, svtype), Collections.emptyList(), PESR_ONLY_ALGORITHM_LIST, Collections.emptyList(),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null);
    }

    public static Integer getLength(final int start, final int end, final GATKSVVCFConstants.StructuralVariantAnnotationType type) {
        if (type == GATKSVVCFConstants.StructuralVariantAnnotationType.BND || type == GATKSVVCFConstants.StructuralVariantAnnotationType.CTX) {
            return null;
        } else if (type == GATKSVVCFConstants.StructuralVariantAnnotationType.INS || type == GATKSVVCFConstants.StructuralVariantAnnotationType.CPX) {
            return 100;
        }
        return end - start + 1;
    }

    public static SVCallRecord newBndCallRecordWithStrands(final boolean strandA, final boolean strandB) {
        return new SVCallRecord("", "chr1", 1000, strandA, "chr1", 1000, strandB, GATKSVVCFConstants.StructuralVariantAnnotationType.BND, null, Collections.emptyList(), null, Collections.emptyList(),
                Collections.singletonList(PESR_ALGORITHM),
                Collections.emptyList(),
                Collections.emptyList(),
                Collections.emptyMap(), Collections.emptySet(), null);
    }

    public static SVCallRecord newCtxCallRecord() {
        return new SVCallRecord("", "chr1", 1000, null, "chr1", 1000, null, GATKSVVCFConstants.StructuralVariantAnnotationType.CTX, null, Collections.emptyList(), null, Collections.emptyList(),
                Collections.singletonList(PESR_ALGORITHM),
                Collections.emptyList(),
                Collections.emptyList(),
                Collections.emptyMap(), Collections.emptySet(), null);
    }

    public static SVCallRecord newCpxCallRecordWithLength(final int length) {
        return new SVCallRecord("", "chr1", 1000, null, "chr1", 1000, null, GATKSVVCFConstants.StructuralVariantAnnotationType.CPX, null, Collections.emptyList(), length, Collections.emptyList(),
                Collections.singletonList(PESR_ALGORITHM),
                Collections.emptyList(),
                Collections.emptyList(),
                Collections.emptyMap(), Collections.emptySet(), null);
    }

    public static SVCallRecord newCnvCallRecordWithStrands(final Boolean strandA, final Boolean strandB) {
        return new SVCallRecord("", "chr1", 1000, strandA, "chr1", 1999, strandB, GATKSVVCFConstants.StructuralVariantAnnotationType.CNV, null, Collections.emptyList(), 1000, Collections.emptyList(),
                Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Collections.emptyList(),
                Collections.emptyList(),
                Collections.emptyMap(), Collections.emptySet(), null);
    }

    public static SVCallRecord newCallRecordWithCoordinates(final String id, final String chrA, final int posA, final String chrB, final int posB) {
        return new SVCallRecord(id, chrA, posA, true, chrB, posB, false, GATKSVVCFConstants.StructuralVariantAnnotationType.BND, null, Collections.emptyList(), null, Collections.emptyList(),
                Collections.singletonList("peser"),
                Collections.emptyList(),
                Collections.emptyList(),
                Collections.emptyMap(), Collections.emptySet(), null);
    }

    public static SVCallRecord newCallRecordWithCoordinatesAndType(final String id, final String chrA, final int posA, final String chrB, final int posB, final GATKSVVCFConstants.StructuralVariantAnnotationType type) {
        return new SVCallRecord(id, chrA, posA, true, chrB, posB, false, type, null, Collections.emptyList(), getLength(posA, posB, type), Collections.emptyList(),
                Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Collections.emptyList(),
                Collections.emptyList(),
                Collections.emptyMap(), Collections.emptySet(), null);
    }

    public static SVCallRecord newCallRecordWithAlgorithms(final List<String> algorithms) {
        return new SVCallRecord("", "chr1", 1000, true, "chr1", 1000, false, GATKSVVCFConstants.StructuralVariantAnnotationType.INS, null, Collections.emptyList(), length, Collections.emptyList(),
                algorithms,
                Collections.emptyList(),
                Collections.emptyList(),
                Collections.emptyMap(), Collections.emptySet(), null);
    }

    public static SVCallRecord newCallRecordInsertionWithLength(final Integer length) {
        return new SVCallRecord("", "chr1", 1000, true, "chr1", 1000, false, GATKSVVCFConstants.StructuralVariantAnnotationType.INS, null, Collections.emptyList(), length, Collections.emptyList(),
                PESR_ONLY_ALGORITHM_LIST,
                Collections.emptyList(),
                Collections.emptyList(),
                Collections.emptyMap(), Collections.emptySet(), null);
    }

    public static SVCallRecord newCallRecordInsertionWithLengthAndCoordinates(final String chrA, final int posA, final Integer length) {
        return new SVCallRecord("", chrA, posA, true, chrA, posA, false, GATKSVVCFConstants.StructuralVariantAnnotationType.INS, null, Collections.emptyList(), length, Collections.emptyList(),
                PESR_ONLY_ALGORITHM_LIST,
                Collections.emptyList(),
                Collections.emptyList(),
                Collections.emptyMap(), Collections.emptySet(), null);
    }

    public static VariantContext newVariantContext(final String id, final String chrA, final int posA,
                                                   final int end, final List<Allele> alleles,
                                                   final List<Genotype> genotypes, final Integer svlen,
                                                   final String strands, final GATKSVVCFConstants.StructuralVariantAnnotationType svtype,
                                                   final List<String> algorithms, final String chr2, final Integer end2,
                                                   final Map<String, Object> extendedAttributes, final Double log10PError) {
        final Map<String, Object> attributes = new HashMap<>();
        attributes.put(GATKSVVCFConstants.SVLEN, svlen);
        attributes.put(GATKSVVCFConstants.SVTYPE, svtype);
        attributes.put(GATKSVVCFConstants.STRANDS_ATTRIBUTE, strands);
        attributes.put(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, algorithms);
        if (chr2 != null) {
            attributes.put(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, chr2);
        }
        if (end2 != null) {
            attributes.put(GATKSVVCFConstants.END2_ATTRIBUTE, end2);
        }
        attributes.putAll(extendedAttributes);
        VariantContextBuilder builder = new VariantContextBuilder()
                .id(id)
                .start(posA)
                .chr(chrA)
                .stop(end)
                .alleles(alleles)
                .genotypes(genotypes)
                .attributes(attributes);
        if (log10PError != null) {
            builder = builder.log10PError(log10PError);
        }
        return builder.make();
    }

    public static Boolean getValidTestStrandA(final GATKSVVCFConstants.StructuralVariantAnnotationType type) {
        if (type == GATKSVVCFConstants.StructuralVariantAnnotationType.DUP) {
            return Boolean.FALSE;
        } else if (type == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV) {
            return null;
        } else if (type == GATKSVVCFConstants.StructuralVariantAnnotationType.INV) {
            return Boolean.TRUE;
        }
        return Boolean.TRUE;
    }

    public static Boolean getValidTestStrandB(final GATKSVVCFConstants.StructuralVariantAnnotationType type) {
        if (type == GATKSVVCFConstants.StructuralVariantAnnotationType.DUP) {
            return Boolean.TRUE;
        } else if (type == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV) {
            return null;
        } else if (type == GATKSVVCFConstants.StructuralVariantAnnotationType.INV) {
            return Boolean.TRUE;
        }
        return Boolean.FALSE;
    }

    public static List<Allele> getCNVAlleles(final GATKSVVCFConstants.StructuralVariantAnnotationType type) {
        if (type == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV) {
            return Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP);
        } else if (type == GATKSVVCFConstants.StructuralVariantAnnotationType.DEL) {
            return Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL);
        } else if (type == GATKSVVCFConstants.StructuralVariantAnnotationType.DUP) {
            return Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP);
        } else {
            throw new TestException("Not a CNV type:" + type);
        }
    }

    public static GenotypeBuilder getDiploidCNVGenotypeBuilder(final String sample, final int copyNumber) {
        return new GenotypeBuilder(sample, Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL))
                .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 2)
                .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, copyNumber);
    }

    public static Map<String, Object> buildMapFromArrays(final String[] keys, final Object[] values) {
        if (keys.length != values.length) {
            throw new TestException("Keys and values have different lengths: " + keys.length + " and " + values.length);
        }
        final Map<String, Object> map = new HashMap<>();
        for (int i = 0; i < keys.length; i++) {
            if (keys[i] == null) {
                throw new TestException("Encountered null key");
            }
            map.put(keys[i], values[i]);
        }
        return map;
    }

    public static SVCallRecord getCNVRecordWithCN(final int ploidy, List<Allele> alleles, final GATKSVVCFConstants.StructuralVariantAnnotationType svtype,
                                           final int[] copyNumbers, final String cnField) {
        // Create genotypes with copy number attribute (and no GT)
        final List<Genotype> genotypesWithCopyNumber = IntStream.range(0, copyNumbers.length)
                .mapToObj(i -> new GenotypeBuilder(String.valueOf(i))
                        .attribute(cnField, copyNumbers[i])
                        .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, ploidy)
                        .alleles(SVTestUtils.buildHomAlleleListWithPloidy(Allele.NO_CALL, ploidy))
                        .make())
                .collect(Collectors.toList());
        return new SVCallRecord("", "chr1", 1000, SVTestUtils.getValidTestStrandA(svtype),
                "chr1", 1999, SVTestUtils.getValidTestStrandB(svtype), svtype, null, Collections.emptyList(),
                1000, Collections.emptyList(), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                alleles, GenotypesContext.copy(genotypesWithCopyNumber), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
    }
}
