package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CoordMath;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.COPY_NUMBER_FORMAT;

public class SVCallRecord implements SVLocatable {

    public static final String STRAND_PLUS = "+";
    public static final String STRAND_MINUS = "-";
    public static final int UNDEFINED_LENGTH = -1;
    public static final List<String> INVALID_ATTRIBUTES = Lists.newArrayList(
            VCFConstants.END_KEY,
            GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE,
            GATKSVVCFConstants.SVLEN,
            GATKSVVCFConstants.CONTIG2_ATTRIBUTE,
            GATKSVVCFConstants.END2_ATTRIBUTE,
            GATKSVVCFConstants.STRANDS_ATTRIBUTE,
            GATKSVVCFConstants.SVTYPE,
            GATKSVVCFConstants.CPX_TYPE
    );

    private final String id;
    private final String contigA;
    private final int positionA;
    private final Boolean strandA;
    private final String contigB;
    private final int positionB;
    private final Boolean strandB;
    private final GATKSVVCFConstants.StructuralVariantAnnotationType type;
    private final Integer length;
    private final List<String> algorithms;
    private final List<Allele> alleles;
    private final Allele refAllele;
    private final List<Allele> altAlleles;
    private final GenotypesContext genotypes;
    private final Map<String,Object> attributes;
    private final Set<String> filters;
    private final Double log10PError;

    // CPX related fields
    private final GATKSVVCFConstants.ComplexVariantSubtype cpxSubtype;

    public SVCallRecord(final String id,
                        final String contigA,
                        final int positionA,
                        final Boolean strandA,
                        final String contigB,
                        final int positionB,
                        final Boolean strandB,
                        final GATKSVVCFConstants.StructuralVariantAnnotationType type,
                        final GATKSVVCFConstants.ComplexVariantSubtype cpxSubtype,
                        final Integer length,
                        final List<String> algorithms,
                        final List<Allele> alleles,
                        final List<Genotype> genotypes,
                        final Map<String,Object> attributes,
                        final Set<String> filters,
                        final Double log10PError,
                        final SAMSequenceDictionary dictionary) {
        this(id, contigA, positionA, strandA, contigB, positionB, strandB, type, cpxSubtype, length, algorithms, alleles, genotypes, attributes, filters, log10PError);
        validateCoordinates(dictionary);
    }

    protected SVCallRecord(final String id,
                           final String contigA,
                           final int positionA,
                           final Boolean strandA,
                           final String contigB,
                           final int positionB,
                           final Boolean strandB,
                           final GATKSVVCFConstants.StructuralVariantAnnotationType type,
                           final GATKSVVCFConstants.ComplexVariantSubtype cpxSubtype,
                           final Integer length,
                           final List<String> algorithms,
                           final List<Allele> alleles,
                           final List<Genotype> genotypes,
                           final Map<String, Object> attributes,
                           final Set<String> filters,
                           final Double log10PError) {
        Utils.nonNull(algorithms);
        Utils.nonNull(alleles);
        Utils.nonNull(genotypes);
        Utils.nonNull(attributes);
        Utils.nonNull(filters);
        this.id = Utils.nonNull(id);
        this.contigA = contigA;
        this.positionA = positionA;
        this.contigB = contigB;
        this.positionB = positionB;
        this.type = Utils.nonNull(type);
        this.cpxSubtype = cpxSubtype;
        this.algorithms = Collections.unmodifiableList(algorithms);
        this.alleles = Collections.unmodifiableList(alleles);
        this.altAlleles = alleles.stream().filter(allele -> !allele.isNoCall() && !allele.isReference()).collect(Collectors.toList());
        final List<Allele> refAllelesList = alleles.stream().filter(allele -> !allele.isNoCall() && allele.isReference()).collect(Collectors.toList());
        Utils.validate(refAllelesList.size() <= 1, "Encountered multiple reference alleles");
        this.refAllele = refAllelesList.isEmpty() ? null : refAllelesList.get(0);
        this.genotypes = GenotypesContext.copy(genotypes).immutable();
        this.attributes = validateAttributes(attributes);
        this.length = inferLength(type, positionA, positionB, length);
        final Pair<Boolean, Boolean> strands = inferStrands(type, strandA, strandB);
        this.strandA = strands.getLeft();
        this.strandB = strands.getRight();
        this.filters = filters;
        this.log10PError = log10PError;
    }

    /**
     * Ensures start/end loci are valid and ordered
     */
    private void validateCoordinates(final SAMSequenceDictionary dictionary) {
        Utils.nonNull(contigA);
        Utils.nonNull(contigB);
        Utils.nonNull(dictionary);
        validatePosition(contigA, positionA, dictionary);
        validatePosition(contigB, positionB, dictionary);
        // CPX types may have position B precede A, such as dispersed duplications where A is the insertion point and
        // B references the source sequence.
        if (type != GATKSVVCFConstants.StructuralVariantAnnotationType.CPX) {
            Utils.validateArg(IntervalUtils.compareLocatables(getPositionAInterval(), getPositionBInterval(), dictionary) <= 0,
                    "End coordinate cannot precede start");
        }
    }

    private static void validatePosition(final String contig, final int position, final SAMSequenceDictionary dictionary) {
        final SAMSequenceRecord seq = dictionary.getSequence(contig);
        Utils.validateArg(seq != null, "Contig " + contig + " not found in dictionary");
        Utils.validateArg(position > 0 && position <= seq.getSequenceLength(), "Invalid position " + contig + ":" + position);
    }

    private static Map<String, Object> validateAttributes(final Map<String, Object> attributes) {
        for (final String key : INVALID_ATTRIBUTES) {
            Utils.validateArg(!attributes.containsKey(key), "Attempted to create record with invalid key: " + key);
        }
        return attributes;
    }

    /**
     * Determines length annotation for the variant. For CNV and INV, length is calculated based on positions and
     * cross-checked with the length if provided. For INS, the provided length is used. In all other cases the
     * inferred length is null.
     * @param type SV type
     * @param positionA first position
     * @param positionB second position
     * @param inputLength assumed length, may be null
     * @return variant length, may be null
     */
    private static Integer inferLength(final GATKSVVCFConstants.StructuralVariantAnnotationType type,
                                       final int positionA,
                                       final int positionB,
                                       final Integer inputLength) {
        if (type == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV || type == GATKSVVCFConstants.StructuralVariantAnnotationType.DEL
                || type == GATKSVVCFConstants.StructuralVariantAnnotationType.DUP || type == GATKSVVCFConstants.StructuralVariantAnnotationType.INV) {
            // Intrachromosomal classes
            final int length = CoordMath.getLength(positionA, positionB);
            if (inputLength != null) {
                Utils.validateArg(inputLength.intValue() == length, "Input length does not match calculated length");
            }
            return length;
        } else {
            if ((type == GATKSVVCFConstants.StructuralVariantAnnotationType.BND
                    || type == GATKSVVCFConstants.StructuralVariantAnnotationType.CTX) && inputLength != null) {
                throw new IllegalArgumentException("Input length should be null for type " + type.name() + " but found " + inputLength);
            }
            return inputLength;
        }
    }

    /**
     * Determines strands for the variant. For SV classes with implicit strand orientations (e.g. DEL, DUP) or that have
     * multiple strand types on each end (i.e. CNVs which are either DEL/DUP which are +/- and -/+), inferred strands
     * are null.
     * @param type SV type
     * @param inputStrandA assumed first strand, may be null
     * @param inputStrandB assumed second strand, may be null
     * @return pair containing (first strand, second strand), which may contain null values
     */
    private static Pair<Boolean, Boolean> inferStrands(final GATKSVVCFConstants.StructuralVariantAnnotationType type,
                                                       final Boolean inputStrandA,
                                                       final Boolean inputStrandB) {
        if (type == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV) {
            Utils.validateArg(inputStrandA == null && inputStrandB == null, "Attempted to create CNV with non-null strands");
            return Pair.of(null, null);
        } else if (type == GATKSVVCFConstants.StructuralVariantAnnotationType.DEL) {
            if (inputStrandA != null) {
                Utils.validateArg(inputStrandA.booleanValue() == true, "Attempted to create DEL with negative first strand");
            }
            if (inputStrandB != null) {
                Utils.validateArg(inputStrandB.booleanValue() == false, "Attempted to create DEL with positive second strand");
            }
            return Pair.of(Boolean.TRUE, Boolean.FALSE);
        } else if (type == GATKSVVCFConstants.StructuralVariantAnnotationType.DUP) {
            if (inputStrandA != null) {
                Utils.validateArg(inputStrandA.booleanValue() == false, "Attempted to create DUP with positive first strand");
            }
            if (inputStrandB != null) {
                Utils.validateArg(inputStrandB.booleanValue() == true, "Attempted to create DUP with negative second strand");
            }
            return Pair.of(Boolean.FALSE, Boolean.TRUE);
        } else {
            return Pair.of(inputStrandA, inputStrandB);
        }
    }


    /**
     * Determines whether the given genotype represents a carrier of the alt allele. Specifically, returns false iff:
     *   (a) there are no alt alleles for this variant,
     *   (b) the sample is ploidy 0,
     *   (c) the GT field contains calls but no alt alleles,
     *   (d) the GT field is not called but the variant is a biallelic CNV with a non-ref CN field (e.g. less than the
     *       ploidy for deletions).
     * Throws an exception if the carrier status cannot be determined from available information.
     */
    private boolean isCarrier(final Genotype genotype) {

        // No alts exist, so can't be a carrier
        if (altAlleles.isEmpty()) {
            return false;
        }

        // Retrieve and check non-zero ploidy
        final int expectedCopyNumber = getExpectedCopyNumber(genotype);
        if (expectedCopyNumber == 0) {
            return false;
        }

        // If a genotype call exists, use it
        if (genotype.isCalled()) {
            return genotype.getAlleles().stream().filter(SVCallRecordUtils::isAltAllele).count() > 0;
        }

        // Otherwise, try to infer status if it's a biallelic CNV with a copy number call
        final int copyNumber = VariantContextGetters.getAttributeAsInt(genotype, COPY_NUMBER_FORMAT, expectedCopyNumber);
        if (type == GATKSVVCFConstants.StructuralVariantAnnotationType.DEL) {
            return copyNumber < expectedCopyNumber;
        } else if (type == GATKSVVCFConstants.StructuralVariantAnnotationType.DUP) {
            return copyNumber > expectedCopyNumber;
        } else if (type == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV) {
            return copyNumber != expectedCopyNumber;
        }

        // No-calls (that aren't bi-allelic CNVs)
        return false;
    }

    public static int getExpectedCopyNumber(final Genotype genotype) {
        Utils.validateArg(genotype.hasExtendedAttribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT),
                "Genotype missing required field " + GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT);
        return VariantContextGetters.getAttributeAsInt(genotype, GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 0);
    }

    public Set<String> getCarrierSampleSet() {
        return getCarrierSampleStream().collect(Collectors.toSet());
    }

    public List<Genotype> getCarrierGenotypeList() {
        return getCarrierGenotypeStream().collect(Collectors.toList());
    }

    private Stream<String> getCarrierSampleStream() {
        return getCarrierGenotypeStream().map(Genotype::getSampleName);
    }

    private Stream<Genotype> getCarrierGenotypeStream() {
        return genotypes.stream().filter(this::isCarrier);
    }

    public boolean isDepthOnly() {
        return algorithms.size() == 1 && algorithms.get(0).equals(GATKSVVCFConstants.DEPTH_ALGORITHM);
    }

    public boolean isSimpleCNV() {
        return type == GATKSVVCFConstants.StructuralVariantAnnotationType.DEL
                || type == GATKSVVCFConstants.StructuralVariantAnnotationType.DUP
                || type == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV;
    }

    public boolean nullStrands() {
        return strandA == null && strandB == null;
    }

    public Map<String, Object> getAttributes() {
        return attributes;
    }

    public String getId() {
        return id;
    }

    @Override
    public String getContigA() {
        return contigA;
    }

    @Override
    public int getPositionA() {
        return positionA;
    }

    @Override
    public String getContigB() {
        return contigB;
    }

    @Override
    public int getPositionB() {
        return positionB;
    }

    public Boolean getStrandA() {
        return strandA;
    }

    public Boolean getStrandB() {
        return strandB;
    }

    @Override
    public GATKSVVCFConstants.StructuralVariantAnnotationType getType() {
        return type;
    }

    public GATKSVVCFConstants.ComplexVariantSubtype getComplexSubtype() {
        return cpxSubtype;
    }

    public Integer getLength() {
        return length;
    }

    public List<String> getAlgorithms() {
        return algorithms;
    }

    public List<Allele> getAlleles() {
        return alleles;
    }

    public List<Allele> getAltAlleles() { return altAlleles; }

    public Allele getRefAllele() { return refAllele; }

    public Set<String> getAllSamples() {
        return genotypes.stream().map(Genotype::getSampleName)
                .collect(Collectors.toCollection(LinkedHashSet::new));
    }

    public GenotypesContext getGenotypes() {
        return genotypes;
    }

    public SimpleInterval getPositionAInterval() {
        return new SimpleInterval(contigA, positionA, positionA);
    }

    public SimpleInterval getPositionBInterval() {
        return new SimpleInterval(contigB, positionB, positionB);
    }

    public Set<String> getFilters() {
        return filters;
    }

    public Double getLog10PError() {
        return log10PError;
    }

    public GATKSVVCFConstants.ComplexVariantSubtype getCpxSubtype() {
        return cpxSubtype;
    }
}
