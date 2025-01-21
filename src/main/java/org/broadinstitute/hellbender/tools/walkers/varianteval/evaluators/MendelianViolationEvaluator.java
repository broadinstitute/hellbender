package org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators;

import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.Analysis;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.DataPoint;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;
import org.broadinstitute.hellbender.utils.samples.MendelianViolation;
import org.broadinstitute.hellbender.utils.samples.SampleDB;

/**
 * Mendelian violation detection and counting
 * <p/>
 * a violation looks like:
 * Suppose dad = A/B and mom = C/D
 * The child can be [A or B] / [C or D].
 * If the child doesn't match this, the site is a violation
 * <p/>
 * Some examples:
 * <p/>
 * mom = A/A, dad = C/C
 * child can be A/C only
 * <p/>
 * mom = A/C, dad = C/C
 * child can be A/C or C/C
 * <p/>
 * mom = A/C, dad = A/C
 * child can be A/A, A/C, C/C
 * <p/>
 * The easiest way to do this calculation is to:
 * <p/>
 * Get alleles for mom => A/B
 * Get alleles for dad => C/D
 * Make allowed genotypes for child: A/C, A/D, B/C, B/D
 * Check that the child is one of these.
 */
@Analysis(name = "Mendelian Violation Evaluator", description = "Mendelian Violation Evaluator")
public class MendelianViolationEvaluator extends VariantEvaluator {
    public MendelianViolationEvaluator(VariantEvalEngine engine) {
        super(engine);

        mv = new ExtendedMendelianViolation(getEngine().getVariantEvalArgs().getMendelianViolationQualThreshold());
    }

    @DataPoint(description = "Number of variants found with at least one family having genotypes", format = "%d")
    public long nVariants;
    @DataPoint(description = "Number of variants found with no family having genotypes -- these sites do not count in the nNoCall", format = "%d")
    public long nSkipped;
    @DataPoint(description="Number of variants x families called (no missing genotype or lowqual)", format = "%d")
    public long nFamCalled;
    @DataPoint(description="Number of variants x families called (no missing genotype or lowqual) that contain at least one var allele.", format = "%d")
    public long nVarFamCalled;
    @DataPoint(description="Number of variants x families discarded as low quality", format = "%d")
    public long nLowQual;
    @DataPoint(description="Number of variants x families discarded as no call", format = "%d")
    public long nNoCall;
    @DataPoint(description="Number of loci with mendelian violations", format = "%d")
    public long nLociViolations;
    @DataPoint(description = "Number of mendelian violations found", format = "%d")
    public long nViolations;

    @DataPoint(description="Number of mendelian violations of the type HOM_REF/HOM_REF -> HOM_VAR", format = "%d")
    public long mvRefRef_Var;
    @DataPoint(description="Number of mendelian violations of the type HOM_REF/HOM_REF -> HET", format = "%d")
    public long mvRefRef_Het;
    @DataPoint(description="Number of mendelian violations of the type HOM_REF/HET -> HOM_VAR", format = "%d")
    public long mvRefHet_Var;
    @DataPoint(description="Number of mendelian violations of the type HOM_REF/HOM_VAR -> HOM_VAR", format = "%d")
    public long mvRefVar_Var;
    @DataPoint(description="Number of mendelian violations of the type HOM_REF/HOM_VAR -> HOM_REF", format = "%d")
    public long mvRefVar_Ref;
    @DataPoint(description="Number of mendelian violations of the type HOM_VAR/HET -> HOM_REF", format = "%d")
    public long mvVarHet_Ref;
    @DataPoint(description="Number of mendelian violations of the type HOM_VAR/HOM_VAR -> HOM_REF", format = "%d")
    public long mvVarVar_Ref;
    @DataPoint(description="Number of mendelian violations of the type HOM_VAR/HOM_VAR -> HET", format = "%d")
    public long mvVarVar_Het;

    @DataPoint(description="Number of HomRef/HomRef/HomRef trios", format = "%d")
    public long HomRefHomRef_HomRef;
    @DataPoint(description="Number of Het/Het/Het trios", format = "%d")
    public long HetHet_Het;
    @DataPoint(description="Number of Het/Het/HomRef trios", format = "%d")
    public long HetHet_HomRef;
    @DataPoint(description="Number of Het/Het/HomVar trios", format = "%d")
    public long HetHet_HomVar;
    @DataPoint(description="Number of HomVar/HomVar/HomVar trios", format = "%d")
    public long HomVarHomVar_HomVar;
    @DataPoint(description="Number of HomRef/HomVar/Het trios", format = "%d")
    public long HomRefHomVAR_Het;
    @DataPoint(description="Number of ref alleles inherited from het/het parents", format = "%d")
    public long HetHet_inheritedRef;
    @DataPoint(description="Number of var alleles inherited from het/het parents", format = "%d")
    public long HetHet_inheritedVar;
    @DataPoint(description="Number of ref alleles inherited from homRef/het parents", format = "%d")
    public long HomRefHet_inheritedRef;
    @DataPoint(description="Number of var alleles inherited from homRef/het parents", format = "%d")
    public long HomRefHet_inheritedVar;
    @DataPoint(description="Number of ref alleles inherited from homVar/het parents", format = "%d")
    public long HomVarHet_inheritedRef;
    @DataPoint(description="Number of var alleles inherited from homVar/het parents", format = "%d")
    public long HomVarHet_inheritedVar;

    ExtendedMendelianViolation mv;

    public String getName() {
        return "mendelian_violations";
    }

    public int getComparisonOrder() {
        return 1;   // we only need to see each eval track
    }

    @Override
    public void update1(final VariantContext vc, final VariantEvalContext context) {
        if (vc.isBiallelic() && vc.hasGenotypes()) { // todo -- currently limited to biallelic loci

            if (mv.countFamilyViolations(getEngine().getSampleDB(), null, vc) > 0){
                nLociViolations++;
                nViolations += mv.getViolationsCount();
                mvRefRef_Var += mv.getParentsRefRefChildVar();
                mvRefRef_Het += mv.getParentsRefRefChildHet();
                mvRefHet_Var += mv.getParentsRefHetChildVar();
                mvRefVar_Var += mv.getParentsRefVarChildVar();
                mvRefVar_Ref += mv.getParentsRefVarChildRef();
                mvVarHet_Ref += mv.getParentsVarHetChildRef();
                mvVarVar_Ref += mv.getParentsVarVarChildRef();
                mvVarVar_Het += mv.getParentsVarVarChildHet();

            }
            HomRefHomRef_HomRef += mv.getRefRefRef();
            HetHet_Het += mv.getHetHetHet();
            HetHet_HomRef += mv.getHetHetHomRef();
            HetHet_HomVar += mv.getHetHetHomVar();
            HomVarHomVar_HomVar += mv.getVarVarVar();
            HomRefHomVAR_Het += mv.getRefVarHet();
            HetHet_inheritedRef += mv.getParentsHetHetInheritedRef();
            HetHet_inheritedVar += mv.getParentsHetHetInheritedVar();
            HomRefHet_inheritedRef += mv.getParentsRefHetInheritedRef();
            HomRefHet_inheritedVar += mv.getParentsRefHetInheritedVar();
            HomVarHet_inheritedRef += mv.getParentsVarHetInheritedRef();
            HomVarHet_inheritedVar += mv.getParentsVarHetInheritedVar();

            if(mv.getFamilyCalledCount()>0 || mv.getFamilyLowQualsCount()>0 || mv.getFamilyNoCallCount()>0){
                nVariants++;
                nFamCalled += mv.getFamilyCalledCount();
                nLowQual += mv.getFamilyLowQualsCount();
                nNoCall += mv.getFamilyNoCallCount();
                nVarFamCalled += mv.getVarFamilyCalledCount();
            }
            else{
                nSkipped++;
            }
        }
    }

    private static class ExtendedMendelianViolation extends MendelianViolation
    {
        public ExtendedMendelianViolation(double threshold)
        {
            super(threshold, false, false);
        }

        /**
         * Count of violations of the type HOM_REF/HOM_REF -> HOM_VAR
         */
        public int getParentsRefRefChildVar(){
            return getInheritance().get(GenotypeType.HOM_REF).get(GenotypeType.HOM_REF).get(GenotypeType.HOM_VAR);
        }

        /**
         * Count of violations of the type HOM_REF/HOM_REF -> HET
         */
        public int getParentsRefRefChildHet(){
            return getInheritance().get(GenotypeType.HOM_REF).get(GenotypeType.HOM_REF).get(GenotypeType.HET);
        }

        /**
         * Count of violations of the type HOM_REF/HET -> HOM_VAR
         */
        public int getParentsRefHetChildVar(){
            return getInheritance().get(GenotypeType.HOM_REF).get(GenotypeType.HET).get(GenotypeType.HOM_VAR) + getInheritance().get(GenotypeType.HET).get(GenotypeType.HOM_REF).get(GenotypeType.HOM_VAR);
        }

        /**
         * Count of violations of the type HOM_REF/HOM_VAR -> HOM_VAR
         */
        public int getParentsRefVarChildVar(){
            return getInheritance().get(GenotypeType.HOM_REF).get(GenotypeType.HOM_VAR).get(GenotypeType.HOM_VAR) + getInheritance().get(GenotypeType.HOM_VAR).get(GenotypeType.HOM_REF).get(GenotypeType.HOM_VAR);
        }

        /**
         * Count of violations of the type HOM_REF/HOM_VAR -> HOM_REF
         */
        public int getParentsRefVarChildRef(){
            return getInheritance().get(GenotypeType.HOM_REF).get(GenotypeType.HOM_VAR).get(GenotypeType.HOM_REF)+ getInheritance().get(GenotypeType.HOM_VAR).get(GenotypeType.HOM_REF).get(GenotypeType.HOM_REF);
        }

        /**
         * Count of violations of the type HOM_VAR/HET -> HOM_REF
         */
        public int getParentsVarHetChildRef(){
            return getInheritance().get(GenotypeType.HET).get(GenotypeType.HOM_VAR).get(GenotypeType.HOM_REF) + getInheritance().get(GenotypeType.HOM_VAR).get(GenotypeType.HET).get(GenotypeType.HOM_REF);
        }

        /**
         * Count of violations of the type HOM_VAR/HOM_VAR -> HOM_REF
         */
        public int getParentsVarVarChildRef(){
            return getInheritance().get(GenotypeType.HOM_VAR).get(GenotypeType.HOM_VAR).get(GenotypeType.HOM_REF);
        }

        /**
         * Count of violations of the type HOM_VAR/HOM_VAR -> HET
         */
        public int getParentsVarVarChildHet(){
            return getInheritance().get(GenotypeType.HOM_VAR).get(GenotypeType.HOM_VAR).get(GenotypeType.HET);
        }

        /**
         * Count of HomRef/HomRef/HomRef trios
         */
        public int getRefRefRef(){
            return getInheritance().get(GenotypeType.HOM_REF).get(GenotypeType.HOM_REF).get(GenotypeType.HOM_REF);
        }

        /**
         * Count of HomVar/HomVar/HomVar trios
         */
        public int getVarVarVar(){
            return getInheritance().get(GenotypeType.HOM_VAR).get(GenotypeType.HOM_VAR).get(GenotypeType.HOM_VAR);
        }

        /**
         * Count of HomRef/HomVar/Het trios
         */
        public int getRefVarHet(){
            return getInheritance().get(GenotypeType.HOM_REF).get(GenotypeType.HOM_VAR).get(GenotypeType.HET) + getInheritance().get(GenotypeType.HOM_VAR).get(GenotypeType.HOM_REF).get(GenotypeType.HET);
        }

        /**
         * Count of Het/Het/Het trios
         */
        public int getHetHetHet(){
            return getInheritance().get(GenotypeType.HET).get(GenotypeType.HET).get(GenotypeType.HET);
        }

        /**
         * Count of Het/Het/HomRef trios
         */
        public int getHetHetHomRef(){
            return getInheritance().get(GenotypeType.HET).get(GenotypeType.HET).get(GenotypeType.HOM_REF);
        }

        /**
         * Count of Het/Het/HomVar trios
         */
        public int getHetHetHomVar(){
            return getInheritance().get(GenotypeType.HET).get(GenotypeType.HET).get(GenotypeType.HOM_VAR);
        }

        /**
         * Count of ref alleles inherited from Het/Het parents (no violation)
         */
        public int getParentsHetHetInheritedRef(){
            return getInheritance().get(GenotypeType.HET).get(GenotypeType.HET).get(GenotypeType.HET) + 2*getInheritance().get(GenotypeType.HET).get(GenotypeType.HET).get(GenotypeType.HOM_REF);
        }

        /**
         * Count of var alleles inherited from Het/Het parents (no violation)
         */
        public int getParentsHetHetInheritedVar(){
            return getInheritance().get(GenotypeType.HET).get(GenotypeType.HET).get(GenotypeType.HET) + 2*getInheritance().get(GenotypeType.HET).get(GenotypeType.HET).get(GenotypeType.HOM_VAR);
        }

        /**
         * Count of ref alleles inherited from HomRef/Het parents (no violation)
         */
        public int getParentsRefHetInheritedRef(){
            return getInheritance().get(GenotypeType.HOM_REF).get(GenotypeType.HET).get(GenotypeType.HOM_REF) + getInheritance().get(GenotypeType.HET).get(GenotypeType.HOM_REF).get(GenotypeType.HOM_REF);
        }

        /**
         * Count of var alleles inherited from HomRef/Het parents (no violation)
         */
        public int getParentsRefHetInheritedVar(){
            return getInheritance().get(GenotypeType.HOM_REF).get(GenotypeType.HET).get(GenotypeType.HET) + getInheritance().get(GenotypeType.HET).get(GenotypeType.HOM_REF).get(GenotypeType.HET);
        }

        /**
         * Count of ref alleles inherited from HomVar/Het parents (no violation)
         */
        public int getParentsVarHetInheritedRef(){
            return getInheritance().get(GenotypeType.HOM_VAR).get(GenotypeType.HET).get(GenotypeType.HET) + getInheritance().get(GenotypeType.HET).get(GenotypeType.HOM_VAR).get(GenotypeType.HET);
        }

        /**
         * Count of var alleles inherited from HomVar/Het parents (no violation)
         */
        public int getParentsVarHetInheritedVar(){
            return getInheritance().get(GenotypeType.HOM_VAR).get(GenotypeType.HET).get(GenotypeType.HOM_VAR) + getInheritance().get(GenotypeType.HET).get(GenotypeType.HOM_VAR).get(GenotypeType.HOM_VAR);
        }
    }
}
