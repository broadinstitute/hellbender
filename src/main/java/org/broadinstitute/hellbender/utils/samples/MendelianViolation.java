package org.broadinstitute.hellbender.utils.samples;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.*;

/**
 * Class for the identification and tracking of mendelian violation. It can be used in 2 distinct ways:
 * - Either using an instance of the MendelianViolation class to track mendelian violations for each of the families while
 * walking over the variants
 * - Or using the static methods to directly get information about mendelian violation in a family at a given locus
 *
 */
public class MendelianViolation {
    //List of families with violations
    private List<String> violationFamilies = new ArrayList<>();

    //Call information
    private int noCall = 0;
    private int familyCalled = 0;
    private int varFamilyCalled = 0;
    private int lowQual = 0;

    //Stores occurrences of inheritance
    private EnumMap<GenotypeType, EnumMap<GenotypeType,EnumMap<GenotypeType,Integer>>> inheritance;

    private int violations_total = 0;

    private final boolean allCalledOnly;
    private final double minGenotypeQuality;
    private final boolean abortOnSampleNotFound;

    /**
     * @param minGenotypeQualityP - the minimum phred scaled genotype quality score necessary to asses mendelian violation
     */
    public MendelianViolation(final double minGenotypeQualityP) {
        this(minGenotypeQualityP, true, false);
    }

    /**
     * @param minGenotypeQualityP - the minimum phred scaled genotype quality score necessary to asses mendelian violation
     * @param abortOnSampleNotFound - Whether to stop execution if a family is passed but no relevant genotypes are found. If false, then the family is ignored.
     * @param completeTriosOnly true if only complete trios are considered, false to include parent/child pairs are
     */
    public MendelianViolation(final double minGenotypeQualityP, final boolean abortOnSampleNotFound, boolean completeTriosOnly) {
        minGenotypeQuality = minGenotypeQualityP;
        this.abortOnSampleNotFound = abortOnSampleNotFound;
        createInheritanceMap();
        allCalledOnly = completeTriosOnly;
    }

    //Count of violations of the type HOM_REF/HOM_REF -> HET
    public int getParentsRefRefChildHet() {
        return inheritance.get(GenotypeType.HOM_REF).get(GenotypeType.HOM_REF).get(GenotypeType.HET);
    }

    private void resetCounts() {
        noCall = 0;
        lowQual = 0;
        familyCalled = 0;
        varFamilyCalled = 0;
        violations_total = 0;
        violationFamilies.clear();

        clearInheritanceMap();
    }

    /**
     * Tests whether there is a mendelian violation between the supplied samples.  Note: this will reset any accumulated stats.
     * @param mother
     * @param father
     * @param child
     * @param vc
     * @return Whether a violation is present
     */
    public boolean isViolation(final Sample mother, final Sample father, final Sample child, final VariantContext vc){
        resetCounts();
        updateViolations(mother.getFamilyID(), mother.getID(), father.getID(), child.getID(),vc);
        return violations_total > 0;
    }

    protected void updateViolations(final String familyId, final String motherId, final String fatherId, final String childId, final VariantContext vc){
        final Genotype gMom = vc.getGenotype(motherId);
        final Genotype gDad = vc.getGenotype(fatherId);
        final Genotype gChild = vc.getGenotype(childId);

        if (gMom == null || gDad == null || gChild == null){
            if(abortOnSampleNotFound) {
                throw new IllegalArgumentException(String.format("Variant %s:%d: Missing genotypes for family %s: mom=%s dad=%s family=%s", vc.getContig(), vc.getStart(), familyId, motherId, fatherId, childId));
            }
            else {
                return;
            }
        }

        //Count No calls
        if(allCalledOnly && (!gMom.isCalled() || !gDad.isCalled() || !gChild.isCalled())){
            noCall++;
        }
        else if (!gMom.isCalled() && !gDad.isCalled() || !gChild.isCalled()){
            noCall++;
        }
        //Count lowQual. Note that if min quality is set to 0, even values with no quality associated are returned
        else if (minGenotypeQuality > 0 && (
            gMom.getGQ()   < minGenotypeQuality ||
            gDad.getGQ()   < minGenotypeQuality ||
            gChild.getGQ() < minGenotypeQuality )) {
            //no call
            lowQual++;
        }
        else {
            //Count all families per loci called
            familyCalled++;

            if (!(gMom.isHomRef() && gDad.isHomRef() && gChild.isHomRef())) {
                varFamilyCalled++;
            }

            if(isViolation(gMom, gDad, gChild)){
                violationFamilies.add(familyId);
                violations_total++;
            }
            final int count = inheritance.get(gMom.getType()).get(gDad.getType()).get(gChild.getType());
            inheritance.get(gMom.getType()).get(gDad.getType()).put(gChild.getType(),count+1);
        }
    }

    /**
     * Counts the number of mendelian violations in the provided sample DB.  Note: this method resets any stats accumulated prior to calling it.
     * @param sampleDB contains the database of samples containing samples for families to be checked for Mendelian violations
     * @param sampleIDs ids of the subset of the samples contained in <code>sampleDB</code> who's family's violations will be checked
     * @param vc the variant context to extract the genotypes and alleles for mom, dad and child.
     * @return whether or not there is a mendelian violation at the site.
     */
    public int countFamilyViolations(SampleDB sampleDB, Set<String> sampleIDs, VariantContext vc) {
        resetCounts();
        Map<String, Set<Sample>> families = sampleDB.getFamilies(sampleIDs);

        for (final Set<Sample> family : families.values()) {
            final Iterator<Sample> sampleIterator = family.iterator();
            Sample sample;
            while (sampleIterator.hasNext()) {
                sample = sampleIterator.next();
                if (!sampleDB.getParents(sample).isEmpty()) {
                    updateViolations(sample.getFamilyID(), sample.getMaternalID(), sample.getPaternalID(), sample.getID(), vc);
                }
            }
        }
        return violations_total;
    }

    /**
     * Evaluate the genotypes of mom, dad, and child to detect Mendelian violations
     *
     * @param gMom
     * @param gDad
     * @param gChild
     * @return true if the three genotypes represent a Mendelian violation; false otherwise
     */
    public static boolean isViolation(final Genotype gMom, final Genotype gDad, final Genotype gChild) {
        switch (gChild.getPloidy()) {
            case 1 : {
                final Allele child_allele = gChild.getAlleles().get(0);
                if (child_allele.isNoCall()){ //cannot possibly be a violation if allele child is no call
                    return false;
                    }
                //one parents is no "call
                if(!gMom.isCalled() || !gDad.isCalled()){
                    return false;
                    }
                return  !gMom.getAlleles().contains(child_allele) && !gDad.getAlleles().contains(child_allele);
                }
            case 2 : {
                if (gChild.isNoCall()){ //cannot possibly be a violation if child is no call
                    return false;
                }
                if(gMom.isHomRef() && gDad.isHomRef() && gChild.isHomRef()) {
                    return false;
                }

                //1 parent is no "call
                if(!gMom.isCalled()){
                    return (gDad.isHomRef() && gChild.isHomVar()) || (gDad.isHomVar() && gChild.isHomRef());
                }
                else if(!gDad.isCalled()){
                    return (gMom.isHomRef() && gChild.isHomVar()) || (gMom.isHomVar() && gChild.isHomRef());
                }
                //Both parents have genotype information
                final Allele childRef = gChild.getAlleles().get(0);
                return !(gMom.getAlleles().contains(childRef) && gDad.getAlleles().contains(gChild.getAlleles().get(1)) ||
                    gMom.getAlleles().contains(gChild.getAlleles().get(1)) && gDad.getAlleles().contains(childRef));
                }
            default: {
                return false;
                }
        }
    }

    private void createInheritanceMap(){
        inheritance = new EnumMap<>(GenotypeType.class);
        for(GenotypeType mType : GenotypeType.values()){
            inheritance.put(mType, new EnumMap<>(GenotypeType.class));
            for(GenotypeType dType : GenotypeType.values()){
                inheritance.get(mType).put(dType, new EnumMap<>(GenotypeType.class));
                for(GenotypeType cType : GenotypeType.values()){
                    inheritance.get(mType).get(dType).put(cType, 0);
                }
            }
        }
    }

    private void clearInheritanceMap(){
        for(GenotypeType mType : GenotypeType.values()){
            for(GenotypeType dType : GenotypeType.values()){
                for(GenotypeType cType : GenotypeType.values()){
                    inheritance.get(mType).get(dType).put(cType, 0);
                }
            }
        }
    }

    /**
     * @return Number of families with genotype information for all members
     */
    public int getFamilyCalledCount(){
        return familyCalled;
    }

    /**
     * @return Number of families where at least one member is not homozygous ref
     */
    public int getVarFamilyCalledCount(){
        return varFamilyCalled;
    }

    /**
     * @return Number of families missing genotypes for one or more of their members
     */
    public int getFamilyNoCallCount(){
        return noCall;
    }

    /**
     * @return Number of families with genotypes below the set quality threshold
     */
    public int getFamilyLowQualsCount(){
        return lowQual;
    }

    /**
     * @return Number of violations identified
     */
    public int getViolationsCount(){
        return violations_total;
    }

    protected EnumMap<GenotypeType, EnumMap<GenotypeType, EnumMap<GenotypeType, Integer>>> getInheritance() {
        return inheritance;
    }
}
