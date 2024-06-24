package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.utils.GeneralUtils;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAssignmentMethod;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.samples.Sample;
import org.broadinstitute.hellbender.utils.samples.SampleDB;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;

/**
 * Utility to compute genotype posteriors given family priors.
 */
public final class FamilyLikelihoods {

    private static final Logger logger = LogManager.getLogger(FamilyLikelihoods.class);

    //Matrix of priors for all genotype combinations
    private final EnumMap<GenotypeType,EnumMap<GenotypeType,EnumMap<GenotypeType,Integer>>> mvCountMatrix =
            new EnumMap<>(GenotypeType.class);

    static final int NUM_CALLED_GENOTYPETYPES = 3; //HOM_REF, HET, and HOM_VAR

    double[] configurationLikelihoodsMatrix = new double[NUM_CALLED_GENOTYPETYPES*NUM_CALLED_GENOTYPETYPES*NUM_CALLED_GENOTYPETYPES];

    private List<Sample> trios = new ArrayList<>();

    public final double NO_JOINT_VALUE = -1.0;

    private double deNovoPrior = 1e-8;

    private static final double ONE_THIRD = 0.333333333333333333;
    private static final double LOG10_OF_ONE_THIRD = -0.4771213;

    private enum FamilyMember {
        MOTHER,
        FATHER,
        CHILD
    }

    public FamilyLikelihoods(final SampleDB sampleDB, final double DNprior, final Set<String> vcfSamples, final Map<String,Set<Sample>> families){
        this.deNovoPrior = DNprior;
        Arrays.fill(configurationLikelihoodsMatrix,0);
        buildMatrices();
        trios = setTrios(sampleDB, vcfSamples, families);
    }

    /**
     * Applies the trio genotype combination to the given trio.
     * @param motherGenotype: Original genotype of the mother
     * @param fatherGenotype: Original genotype of the father
     * @param childGenotype: Original genotype of the child
     * @param updatedGenotypes: An ArrayList<Genotype> to which the newly updated genotypes are added in the following order: Mother, Father, Child
     */
    public void getUpdatedGenotypes(final VariantContext vc, final Genotype motherGenotype, final Genotype fatherGenotype, final Genotype childGenotype, final ArrayList<Genotype> updatedGenotypes){
        //genotypes here can be no call
        final boolean fatherIsCalled = fatherGenotype != null && hasCalledGT(fatherGenotype.getType()) && fatherGenotype.hasLikelihoods();
        final boolean motherIsCalled = motherGenotype != null && hasCalledGT(motherGenotype.getType()) && motherGenotype.hasLikelihoods();
        final boolean childIsCalled = childGenotype != null && hasCalledGT(childGenotype.getType()) && childGenotype.hasLikelihoods();

        //default to posteriors equal to likelihoods (flat priors) in case input genotypes are not called
        final double[] uninformativeLikelihoods = {ONE_THIRD, ONE_THIRD, ONE_THIRD};

        final double[] motherLikelihoods = motherIsCalled? GeneralUtils.normalizeFromLog10(motherGenotype.getLikelihoods().getAsVector()) : uninformativeLikelihoods;
        final double[] fatherLikelihoods = fatherIsCalled? GeneralUtils.normalizeFromLog10(fatherGenotype.getLikelihoods().getAsVector()) : uninformativeLikelihoods;
        final double[] childLikelihoods = childIsCalled? GeneralUtils.normalizeFromLog10(childGenotype.getLikelihoods().getAsVector()) : uninformativeLikelihoods;

        //these are also in log10 space
        final double[] motherLog10Posteriors = getPosteriors(FamilyMember.MOTHER);
        final double[] fatherLog10Posteriors = getPosteriors(FamilyMember.FATHER);
        final double[] childLog10Posteriors = getPosteriors(FamilyMember.CHILD);

        final double[] motherPosteriors = GeneralUtils.normalizeFromLog10(motherLog10Posteriors);
        final double[] fatherPosteriors = GeneralUtils.normalizeFromLog10(fatherLog10Posteriors);
        final double[] childPosteriors = GeneralUtils.normalizeFromLog10(childLog10Posteriors);


        double jointPosteriorProbability =  -1;
        //jointTrioLikelihood is combined likelihoods (before prior) of best configuration after applying prior
        double jointTrioLikelihood = -1;
        if(childIsCalled && motherIsCalled && fatherIsCalled) {
            jointTrioLikelihood = motherLikelihoods[MathUtils.maxElementIndex(motherPosteriors)]*fatherLikelihoods[MathUtils.maxElementIndex(fatherPosteriors)]*childLikelihoods[MathUtils.maxElementIndex(childPosteriors)];
            jointPosteriorProbability = MathUtils.arrayMax(motherPosteriors)*MathUtils.arrayMax(fatherPosteriors)*MathUtils.arrayMax(childPosteriors);
        }

        updatedGenotypes.add(getUpdatedGenotype(vc, motherGenotype, jointTrioLikelihood, jointPosteriorProbability, motherLog10Posteriors));
        updatedGenotypes.add(getUpdatedGenotype(vc, fatherGenotype, jointTrioLikelihood, jointPosteriorProbability, fatherLog10Posteriors));
        updatedGenotypes.add(getUpdatedGenotype(vc, childGenotype, jointTrioLikelihood, jointPosteriorProbability, childLog10Posteriors));
    }

    private Genotype getUpdatedGenotype(final VariantContext vc, final Genotype genotype, final double jointLikelihood, final double jointPosteriorProb, final double[] log10Posteriors){
        //Don't update null, missing or unavailable genotypes
        if(genotype == null || !hasCalledGT(genotype.getType())) {
            return genotype;
        }

        int phredScaledJL = -1;
        int phredScaledJP = -1;
        if(jointLikelihood != NO_JOINT_VALUE){
            final double dphredScaledJL = QualityUtils.phredScaleLog10ErrorRate(Math.log10(1-jointLikelihood));
            phredScaledJL = dphredScaledJL < Byte.MAX_VALUE ? (byte)dphredScaledJL : Byte.MAX_VALUE;
        }
        if(jointPosteriorProb != NO_JOINT_VALUE){
            final double dphredScaledJP = QualityUtils.phredScaleLog10ErrorRate(Math.log10(1-jointPosteriorProb));
            phredScaledJP = dphredScaledJP < Byte.MAX_VALUE ? (byte)dphredScaledJP : Byte.MAX_VALUE;
        }

        //Add the joint trio calculations
        final Map<String, Object> genotypeAttributes = new LinkedHashMap<>();
        genotypeAttributes.putAll(genotype.getExtendedAttributes());
        genotypeAttributes.put(GATKVCFConstants.JOINT_LIKELIHOOD_TAG_NAME, phredScaledJL);
        genotypeAttributes.put(GATKVCFConstants.JOINT_POSTERIOR_TAG_NAME, phredScaledJP);

        final GenotypeBuilder builder = new GenotypeBuilder(genotype);

        //update genotype types based on posteriors
        GATKVariantContextUtils.makeGenotypeCall(vc.getMaxPloidy(2), builder, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, log10Posteriors, vc.getAlleles(), genotype, null);

        builder.attribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY,
                Utils.listFromPrimitives(GenotypeLikelihoods.fromLog10Likelihoods(log10Posteriors).getAsPLs()));
        builder.attributes(genotypeAttributes);
        return builder.make();
    }

    //marginalize over the configurationLikelihoodsMatrix and normalize to get the posteriors
    private double[] getPosteriors(final FamilyMember recalcInd) {
        final double[] marginalOverChangedHR = new double[NUM_CALLED_GENOTYPETYPES*NUM_CALLED_GENOTYPETYPES];
        final double[] marginalOverChangedHET = new double[NUM_CALLED_GENOTYPETYPES*NUM_CALLED_GENOTYPETYPES];
        final double[] marginalOverChangedHV = new double[NUM_CALLED_GENOTYPETYPES*NUM_CALLED_GENOTYPETYPES];
        final double[] recalcPosteriors = new double[NUM_CALLED_GENOTYPETYPES];

        final GenotypeType[] calledTypes = {GenotypeType.HOM_REF, GenotypeType.HET, GenotypeType.HOM_VAR};
        int counter = 0;

        switch (recalcInd) {
            case MOTHER:
                for(final GenotypeType father : calledTypes) {
                    for(final GenotypeType child : calledTypes) {
                        GenotypeType mother;
                        mother = GenotypeType.HOM_REF;
                        marginalOverChangedHR[counter] = configurationLikelihoodsMatrix[getLikelihoodMatrixIndex(mother, father, child)];
                        mother = GenotypeType.HET;
                        marginalOverChangedHET[counter] = configurationLikelihoodsMatrix[getLikelihoodMatrixIndex(mother, father, child)];
                        mother = GenotypeType.HOM_VAR;
                        marginalOverChangedHV[counter] = configurationLikelihoodsMatrix[getLikelihoodMatrixIndex(mother, father, child)];
                        counter++;
                    }
                }
                break;
            case FATHER:
                for(final GenotypeType mother : calledTypes){
                    for (final GenotypeType child : calledTypes){
                        GenotypeType father;
                        father = GenotypeType.HOM_REF;
                        marginalOverChangedHR[counter] = configurationLikelihoodsMatrix[getLikelihoodMatrixIndex(mother, father, child)];
                        father = GenotypeType.HET;
                        marginalOverChangedHET[counter] = configurationLikelihoodsMatrix[getLikelihoodMatrixIndex(mother, father, child)];
                        father = GenotypeType.HOM_VAR;
                        marginalOverChangedHV[counter] = configurationLikelihoodsMatrix[getLikelihoodMatrixIndex(mother, father, child)];
                        counter++;
                    }
                }
                break;
            case CHILD:
                for(final GenotypeType mother : calledTypes){
                    for (final GenotypeType father: calledTypes){
                        GenotypeType child;
                        child = GenotypeType.HOM_REF;
                        marginalOverChangedHR[counter] = configurationLikelihoodsMatrix[getLikelihoodMatrixIndex(mother, father, child)];
                        child = GenotypeType.HET;
                        marginalOverChangedHET[counter] = configurationLikelihoodsMatrix[getLikelihoodMatrixIndex(mother, father, child)];
                        child = GenotypeType.HOM_VAR;
                        marginalOverChangedHV[counter] = configurationLikelihoodsMatrix[getLikelihoodMatrixIndex(mother, father, child)];
                        counter++;
                    }
                }
                break;
            default:
                throw new UserException(String.format("%d does not indicate a valid trio FamilyMember -- use 0 for mother, 1 for father, 2 for child",recalcInd.ordinal()));
        }

        recalcPosteriors[0] = MathUtils.log10sumLog10(marginalOverChangedHR,0);
        recalcPosteriors[1] = MathUtils.log10sumLog10(marginalOverChangedHET,0);
        recalcPosteriors[2] = MathUtils.log10sumLog10(marginalOverChangedHV,0);

        return MathUtils.scaleLogSpaceArrayForNumericalStability(recalcPosteriors);
    }

    /**
     * Computes phred-scaled genotype posteriors given the data in the given variant context and family priors given by this object.
     */
    public GenotypesContext calculatePosteriorGLs(final VariantContext vc){
        Utils.nonNull(vc);
        final GenotypesContext genotypesContext = GenotypesContext.copy(vc.getGenotypes());

        for (final Sample sample : trios) {
            final Genotype mother = vc.getGenotype(sample.getMaternalID());
            final Genotype father = vc.getGenotype(sample.getPaternalID());
            final Genotype child = vc.getGenotype(sample.getID());

            //Keep only trios and parent/child pairs
            if(mother == null && father == null || child == null) {
                logger.warn("Null genotypes in variant: "+vc.toStringDecodeGenotypes());
                continue;
            }

            final ArrayList<Genotype> trioGenotypes = new ArrayList<>(3);
            updateFamilyGenotypes(vc, mother, father, child, trioGenotypes);

            //replace uses sample names to match genotypes, so order doesn't matter
            if (!trioGenotypes.isEmpty()) {
                genotypesContext.replace(trioGenotypes.get(0));
                genotypesContext.replace(trioGenotypes.get(1));
                genotypesContext.replace(trioGenotypes.get(2));
            }
        }

        return genotypesContext;
    }

    /**
     * Select trios and parent/child pairs only
     */
    private List<Sample> setTrios(final SampleDB sampleDB, final Set<String> vcfSamples, final Map<String, Set<Sample>> families){
        final List<Sample> trios = new ArrayList<>();
        for(final Map.Entry<String,Set<Sample>> familyEntry : families.entrySet()){
            Set<Sample> family = familyEntry.getValue();

            // Since getFamilies(vcfSamples) above still returns parents of samples in the VCF even if those parents are not in the VCF, need to subset down here:
            final Set<Sample> familyMembersInVCF = new TreeSet<>();
            for(final Sample familyMember : family){
                if (vcfSamples.contains(familyMember.getID())) {
                    familyMembersInVCF.add(familyMember);
                }
            }
            family = familyMembersInVCF;

            if(family.size() == 3){
                for(final Sample familyMember : family){
                    final List<Sample> parents = sampleDB.getParents(familyMember);
                    if(parents.size()==2){
                        if(family.containsAll(parents)) {
                            trios.add(familyMember);
                        }
                    }
                }
            }

        }
        return trios;
    }

    //Create a lookup matrix to find the number of MVs for each family genotype combination
    private void buildMatrices(){
        for(final GenotypeType mother : GenotypeType.values()){
            mvCountMatrix.put(mother, new EnumMap<>(GenotypeType.class));
            for(final GenotypeType father : GenotypeType.values()){
                mvCountMatrix.get(mother).put(father, new EnumMap<>(GenotypeType.class));
                for(final GenotypeType child : GenotypeType.values()){
                    mvCountMatrix.get(mother).get(father).put(child, getCombinationMVCount(mother, father, child));
                }
            }
        }
    }

    //Returns the number of Mendelian Violations for a given genotype combination.
    //If one of the parents' genotypes is missing, it will consider it as a parent/child pair
    //If the child genotype or both parents genotypes are missing, 0 is returned.
    private int getCombinationMVCount(final GenotypeType mother, final GenotypeType father, final GenotypeType child){

        //Child is no call => No MV
        if(child == GenotypeType.NO_CALL || child == GenotypeType.UNAVAILABLE) {
            return 0;
        }
        //Add parents with genotypes for the evaluation
        final ArrayList<GenotypeType> parents = new ArrayList<>();
        if (!(mother == GenotypeType.NO_CALL || mother == GenotypeType.UNAVAILABLE)) {
            parents.add(mother);
        }
        if (!(father == GenotypeType.NO_CALL || father == GenotypeType.UNAVAILABLE)) {
            parents.add(father);
        }

        //Both parents no calls => No MV
        if (parents.isEmpty()) {
            return 0;
        }

        //If at least one parent had a genotype, then count the number of ref and alt alleles that can be passed
        int parentsNumRefAlleles = 0;
        int parentsNumAltAlleles = 0;

        for(final GenotypeType parent : parents){
            if(parent == GenotypeType.HOM_REF){
                parentsNumRefAlleles++;
            }
            else if(parent == GenotypeType.HET){
                parentsNumRefAlleles++;
                parentsNumAltAlleles++;
            }
            else if(parent == GenotypeType.HOM_VAR){
                parentsNumAltAlleles++;
            }
        }

        //Case Child is HomRef
        if(child == GenotypeType.HOM_REF){
            if(parentsNumRefAlleles == parents.size()) {
                return 0;
            } else {
                return (parents.size() - parentsNumRefAlleles);
            }
        }

        //Case child is HomVar
        if(child == GenotypeType.HOM_VAR){
            if(parentsNumAltAlleles == parents.size()) {
                return 0;
            } else {
                return parents.size() - parentsNumAltAlleles;
            }
        }

        //Case child is Het
        if(child == GenotypeType.HET && ((parentsNumRefAlleles > 0 && parentsNumAltAlleles > 0) || parents.size()<2)) {
            return 0;
        }

        //MV
        return 1;
    }

    /**
     * Updates the genotypes of the given trio. If one of the parents is null, it is considered a parent/child pair.
     * @param vc: Input variant context
     * @param mother: Mother's genotype from vc input
     * @param father: Father's genotype from vc input
     * @param child: Child's genotype from vc input
     * @param finalGenotypes: An ArrayList<Genotype> containing the updated genotypes
     */
    private void updateFamilyGenotypes(final VariantContext vc, final Genotype mother, final Genotype father, final Genotype child, final ArrayList<Genotype> finalGenotypes) {

        //If one of the parents is not called, fill in with uninformative likelihoods
        final Map<GenotypeType,Double> motherLikelihoods = getLikelihoodsAsMapSafeNull(mother);
        final Map<GenotypeType,Double> fatherLikelihoods = getLikelihoodsAsMapSafeNull(father);
        final Map<GenotypeType,Double> childLikelihoods = getLikelihoodsAsMapSafeNull(child);

        //if the child isn't called or neither parent is called, there's no extra inheritance information in that trio so return
        if (!hasCalledGT(child.getType()) || (!hasCalledGT(mother.getType()) && !hasCalledGT(father.getType()))) {
            return;
        }

        //Fill the configurationLikelihoodsMatrix for each genotype combination
        int matInd;
        int mvCount;
        double jointLikelihood;
        double mvCoeff;
        double configurationLikelihood;
        for(final Map.Entry<GenotypeType,Double> childGenotype :
                childLikelihoods.entrySet()){
            for(final Map.Entry<GenotypeType,Double> motherGenotype :
                    motherLikelihoods.entrySet()){
                for(final Map.Entry<GenotypeType,Double> fatherGenotype :
                        fatherLikelihoods.entrySet()){
                    mvCount = mvCountMatrix.get(motherGenotype.getKey()).get(fatherGenotype.getKey()).get(childGenotype.getKey());
                    jointLikelihood = motherGenotype.getValue()+fatherGenotype.getValue()+childGenotype.getValue();
                    mvCoeff = mvCount>0 ? Math.pow(deNovoPrior,mvCount) : (1.0-10*deNovoPrior-deNovoPrior*deNovoPrior);
                    configurationLikelihood =  Math.log10(mvCoeff) + jointLikelihood;
                    matInd = getLikelihoodMatrixIndex(motherGenotype.getKey(), fatherGenotype.getKey(), childGenotype.getKey());
                    configurationLikelihoodsMatrix[matInd] = configurationLikelihood;
                }
            }
        }

        getUpdatedGenotypes(vc, mother, father, child, finalGenotypes);
    }

    //Get a Map of genotype (log10)likelihoods
    private EnumMap<GenotypeType,Double> getLikelihoodsAsMapSafeNull(final Genotype genotype){
        final EnumMap<GenotypeType,Double> likelihoodsMap = new EnumMap<>(GenotypeType.class);
        final double[] likelihoods;

        if (genotype != null && hasCalledGT(genotype.getType()) && genotype.hasExtendedAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY)) {
            final Object GPfromVCF = genotype.getExtendedAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY);
            //parse the GPs into a vector of probabilities
            final String[] likelihoodsAsStringVector = ((String)GPfromVCF).split(",");
            final double[] likelihoodsAsVector = new double[likelihoodsAsStringVector.length];
            for ( int i = 0; i < likelihoodsAsStringVector.length; i++ ) {
                likelihoodsAsVector[i] = Double.parseDouble(likelihoodsAsStringVector[i]) / -10.0;
            }
            //keep in log10 space for large GQs
            likelihoods = GeneralUtils.normalizeFromLog10(likelihoodsAsVector, true, true);
        }

        //In case of null, unavailable or no call, all likelihoods are log10(1/3)
        else if(genotype == null || !hasCalledGT(genotype.getType()) || genotype.getLikelihoods() == null){
            likelihoods = new double[NUM_CALLED_GENOTYPETYPES];
            likelihoods[0] = LOG10_OF_ONE_THIRD;
            likelihoods[1] = LOG10_OF_ONE_THIRD;
            likelihoods[2] = LOG10_OF_ONE_THIRD;
        }

        //No posteriors in VC, use PLs
        else {
            likelihoods = GeneralUtils.normalizeFromLog10(genotype.getLikelihoods().getAsVector(), true, true);
        }

        if (likelihoods.length != NUM_CALLED_GENOTYPETYPES) {
            final String key = genotype.hasExtendedAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY) ?
                    GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY : VCFConstants.GENOTYPE_PL_KEY;
            throw new UserException(genotype + " has " + likelihoods.length + " " + key + " values, should be " + NUM_CALLED_GENOTYPETYPES +
                    " since only the diploid case is supported when applying family priors.");
        }

        likelihoodsMap.put(GenotypeType.HOM_REF,likelihoods[genotypeTypeToValue(GenotypeType.HOM_REF)]);
        likelihoodsMap.put(GenotypeType.HET,likelihoods[genotypeTypeToValue(GenotypeType.HET)]);
        likelihoodsMap.put(GenotypeType.HOM_VAR, likelihoods[genotypeTypeToValue(GenotypeType.HOM_VAR)]);
        return likelihoodsMap;
    }

    private int getLikelihoodMatrixIndex(final GenotypeType mother, final GenotypeType father, final GenotypeType child){
        final int childInd = genotypeTypeToValue(child);
        final int motherInd;
        final int fatherInd;
        final int INVALID = -1;
        motherInd = genotypeTypeToValue(mother);
        fatherInd = genotypeTypeToValue(father);

        if (childInd == INVALID || motherInd == INVALID || fatherInd == INVALID) //any of the genotypes are NO_CALL, UNAVAILABLE or MIXED
        {
            return INVALID;
        }

        //index into array playing the part of a 3x3x3 matrix (where 3=NUM_CALLED_GENOTYPETYPES)
        return motherInd*NUM_CALLED_GENOTYPETYPES*NUM_CALLED_GENOTYPETYPES + fatherInd*NUM_CALLED_GENOTYPETYPES + childInd;
    }

    private int genotypeTypeToValue(final GenotypeType input){
        if (input == GenotypeType.HOM_REF) {
            return 0;
        }
        if (input == GenotypeType.HET) {
            return 1;
        }
        if (input == GenotypeType.HOM_VAR) {
            return 2;
        }
        return -1;
    }

    //this excludes mixed genotypes, whereas the htsjdk Genotype.isCalled() will return true if the GenotypeType is mixed
    private boolean hasCalledGT(final GenotypeType genotype){
        return genotype == GenotypeType.HOM_REF || genotype == GenotypeType.HET || genotype == GenotypeType.HOM_VAR;
    }

}