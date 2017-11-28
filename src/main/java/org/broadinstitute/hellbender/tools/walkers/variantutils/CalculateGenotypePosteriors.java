package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.samples.*;
import org.broadinstitute.hellbender.utils.variant.*;

import java.io.File;
import java.util.*;

/**
 * Calculate genotype posterior probabilities given family and/or known population genotypes
 *
 * <p>
 * This tool calculates the posterior genotype probability for each sample genotype in a VCF of input variant calls,
 * based on the genotype likelihoods from the samples themselves and, optionally, from input VCFs describing allele
 * frequencies in related populations. The input variants must possess genotype likelihoods generated by
 * HaplotypeCaller, UnifiedGenotyper or another source that provides <b>unbiased</b> genotype likelihoods.</p>
 *
 * <h4>Statistical notes</h4>
 * <p>The AF field is not used in the calculation as it does not provide a way to estimate the confidence
 * interval or uncertainty around the allele frequency, unlike AN which does provide this necessary information. This
 * uncertainty is modeled by a Dirichlet distribution: that is, the frequency is known up to a Dirichlet distribution
 * with parameters AC1+q,AC2+q,...,(AN-AC1-AC2-...)+q, where "q" is the global frequency prior (typically q << 1). The
 * genotype priors applied then follow a Dirichlet-Multinomial distribution, where 2 alleles per sample are drawn
 * independently. This assumption of independent draws follows from the assumption of Hardy-Weinberg equilibrium (HWE).
 * Thus, HWE is imposed on the likelihoods as a result of CalculateGenotypePosteriors.</p>
 *
 * <h3>Inputs</h3>
 * <p>
 *     <ul>
 *         <li>A VCF with genotype likelihoods, and optionally genotypes, AC/AN fields, or MLEAC/AN fields.</li>
 *         <li>(Optional) A PED pedigree file containing the description of the relationships between individuals.</li>
 *     </ul>
 * </p>
 *
 * <p>
 * Optionally, a collection of VCFs can be provided for the purpose of informing allele frequency priors. Each of
 * these resource VCFs must satisfy at least one of the following requirement sets:
 * </p>
 * <ul>
 *     <li>AC field and AN field</li>
 *     <li>MLEAC field and AN field</li>
 *     <li>Genotypes</li>
 * </ul>
 * </p>
 *
 * <h3>Output</h3>
 * <p>A new VCF with the following information:</p>
 * <ul>
 *     <li>Genotype posteriors added to the FORMAT fields ("PP")</li>
 *     <li>Genotypes and GQ assigned according to these posteriors (note that the original genotype and GQ may change)</li>
 *     <li>Per-site genotype priors added to the INFO field ("PG")</li>
 *     <li>(Optional) Per-site, per-trio joint likelihoods (JL) and joint posteriors (JL) given as Phred-scaled probability
 *  of all genotypes in the trio being correct based on the PLs for JL and the PPs for JP. These annotations are added to
 *  the FORMAT fields.</li>
 * </ul>
 *
 * <h3>Notes</h3>
 * <p>
 * By default, priors will be applied to each variant separately, provided each variant features data from at least
 * 10 called samples (no-calls do not count). SNP sites in the input callset that have a SNP at the matching site in
 * the supporting VCF will have priors applied based on the AC from the supporting samples and the input callset
 * unless the --ignoreInputSamples flag is used. If a site is not called in the supporting VCF, priors will be
 * applied using the discovered AC from the input samples unless the --discoveredACpriorsOff flag is used.
 * For any non-SNP sites in the input callset, flat priors are applied.
 * </p>
 *
 * <h3>Usage examples</h3>
 *
 * <h4>Refine genotypes based on the discovered allele frequency in an input VCF containing many samples</h4>
 * <pre>
 * gatk --java-options "-Xmx4g" CalculateGenotypePosteriors \
 *   -V multisample_input.vcf.gz \
 *   -O output.vcf.gz
 * </pre>
 *
 * <h4>Inform the genotype assignment of a single sample using the 1000G phase 3 samples</h4>
 * <pre>
 * gatk --java-options "-Xmx4g" CalculateGenotypePosteriors \
 *   -V sample_input.vcf.gz \
 *   -O sample_output.1000G_PPs.vcf.gz \
 *   -supporting 1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf.gz
 * </pre>
 *
 * <h4>Apply only family priors to a callset</h4>
 * <pre>
 * gatk --java-options "-Xmx4g" CalculateGenotypePosteriors \
 *   -V input.vcf.gz \
 *   -O output.vcf.gz \
 *   -ped family.ped \
 *   --skipPopulationPriors
 * </pre>
 *
 * <h4>Apply frequency and HWE-based priors to the genotypes of a family without including the family allele counts
 * in the allele frequency estimates</h4>
 * <pre>
 * gatk --java-options "-Xmx4g" CalculateGenotypePosteriors \
 *   -V input.vcf.gz \
 *   -O output.vcf.gz \
 *   --ignoreInputSamples
 * </pre>
 *
 * <h4>Calculate the posterior genotypes of a callset, and impose that a variant *not seen* in the external panel
 * is tantamount to being AC=0, AN=5008 within that panel</h4>
 * <pre>
 * gatk --java-options "-Xmx4g" CalculateGenotypePosteriors \
 *   -V input.vcf.gz \
 *   -O output.vcf.gz \
 *   -supporting 1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf.gz \
 *   --num-reference-samples-if-no-call 2504
 * </pre>
 *
 * <h3>Caveat</h3>
 * <p>If applying family priors, only diploid family genotypes are supported</p>
 */
@CommandLineProgramProperties(
        summary = "This tool calculates the posterior genotype probability for each sample genotype in a VCF of input variant calls,\n" +
                " based on the genotype likelihoods from the samples themselves and, optionally, from input VCFs describing allele\n" +
                " frequencies in related populations. The input variants must possess genotype likelihoods generated by\n" +
                " HaplotypeCaller, UnifiedGenotyper or another source that provides *unbiased* genotype likelihoods.",
        oneLineSummary = "Calculate genotype posterior probabilities given family and/or known population genotypes",
        programGroup = VariantProgramGroup.class
)
@DocumentedFeature
public final class CalculateGenotypePosteriors extends VariantWalker {

    private static final Logger logger = LogManager.getLogger(CalculateGenotypePosteriors.class);

    /**
     * Supporting external panels. Allele counts from these panels (taken from AC,AN or MLEAC,AN or raw genotypes) will
     * be used to inform the frequency distribution underlying the genotype priors. These files must be VCF 4.2 spec or later.
     */
    @Argument(fullName="supporting-callsets", shortName = "supporting", doc="Other callsets to use in generating genotype posteriors", optional=true)
    public List<FeatureInput<VariantContext>> supportVariants = new ArrayList<>();

    @Argument(doc="File to which variants should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public String out = null;

    /**
     * The global prior of a variant site -- i.e. the expected allele frequency distribution knowing only that N alleles
     * exist, and having observed none of them. This is the "typical" 1/x trend, modeled here as not varying
     * across alleles. The calculation for this parameter is (Effective population size) * (steady state mutation rate)
     *
     */
     @Argument(fullName="global-prior", doc="Global Dirichlet prior parameters for the allele frequency",optional=true)
     public double globalPrior = HomoSapiensConstants.SNP_HETEROZYGOSITY;

    /**
     * The mutation prior -- i.e. the probability that a new mutation occurs. Sensitivity analysis on known de novo 
     * mutations suggests a default value of 10^-6.
     *
     */
    @Argument(fullName="de-novo-prior", doc="Prior for de novo mutations",optional=true)
    public double deNovoPrior = 1e-6;

    /**
     * When a variant is not seen in a panel, this argument controls whether to infer (and with what effective strength)
     * that only reference alleles were observed at that site. E.g. "If not seen in 1000Genomes, treat it as AC=0,
     * AN=2000". This is applied across all external panels, so if numRefIsMissing = 10, and the variant is absent in
     * two panels, this confers evidence of AC=0,AN=20.
     */
    @Argument(fullName="num-reference-samples-if-no-call",doc="Number of hom-ref genotypes to infer at sites not present in a panel",optional=true)
    public int numRefIfMissing = 0;

    /**
     * By default the tool looks for MLEAC first, and then falls back to AC if MLEAC is not found. When this
     * flag is set, the behavior is flipped and the tool looks first for the AC field and then fall back to MLEAC or
     * raw genotypes.
     */
    @Argument(fullName="default-to-allele-count",doc="Use AC rather than MLEAC",optional=true)
    public boolean defaultToAC = false;

    /**
     * When this flag is set, only the AC and AN calculated from external sources will be used, and the calculation
     * will not use the discovered allele frequency in the callset whose posteriors are being calculated. Useful for
     * callsets containing related individuals.
     */
    @Argument(fullName="ignore-input-samples",doc="Use external information only",optional=true)
    public boolean ignoreInputSamples = false;

    /**
     * Calculate priors for missing external variants from sample data -- default behavior is to apply flat priors
     */
    @Argument(fullName="discovered-allele-count-priors-off",doc="Do not use discovered allele count in the input callset " +
            "for variants that do not appear in the external callset. ", optional=true)
    public boolean useACoff = false;

    /**
     * Skip application of population-based priors
     */
    @Argument(fullName="skip-population-priors",doc="Skip application of population-based priors", optional=true)
    public boolean skipPopulationPriors = false;

    /**
     * Skip application of family-based priors. Note: if pedigree file is absent, family-based priors will always be skipped.
     */
    @Argument(fullName="skip-family-priors",doc="Skip application of family-based priors", optional=true)
    public boolean skipFamilyPriors = false;

    /**
     * See https://software.broadinstitute.org/gatk/documentation/article.php?id=7696 for more details on the PED
     * format. Note that each -ped argument can be tagged with NO_FAMILY_ID, NO_PARENTS, NO_SEX, NO_PHENOTYPE to
     * tell the GATK PED parser that the corresponding fields are missing from the ped file.
     *
     */
    @Argument(fullName=StandardArgumentDefinitions.PEDIGREE_FILE_LONG_NAME, shortName=StandardArgumentDefinitions.PEDIGREE_FILE_SHORT_NAME, doc="Pedigree file for samples", optional=true)
    private File pedigreeFile = null;

    private FamilyLikelihoods famUtils;
    private SampleDB sampleDB = null;

    private VariantContextWriter vcfWriter;

    @Override
    public void onTraversalStart() {
        final VariantContextWriterBuilder builder = new VariantContextWriterBuilder().setOutputFile(out).setOutputFileType(VariantContextWriterBuilder.OutputType.VCF);
        if (hasReference()){
            vcfWriter = builder.setReferenceDictionary(getBestAvailableSequenceDictionary()).setOption(Options.INDEX_ON_THE_FLY).build();
        } else {
            vcfWriter = builder.unsetOption(Options.INDEX_ON_THE_FLY).build();
            logger.info("Can't make an index for output file " + out + " because a reference dictionary is required for creating Tribble indices on the fly");
        }

        sampleDB = initializeSampleDB();

        // Get list of samples to include in the output
        final Map<String, VCFHeader> vcfHeaders = Collections.singletonMap(getDrivingVariantsFeatureInput().getName(), getHeaderForVariants());
        final Set<String> vcfSamples = VcfUtils.getSortedSampleSet(vcfHeaders, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);

        //Get the trios from the families passed as ped
        if (!skipFamilyPriors){
            final Set<Trio> trios = sampleDB.getTrios();
            if(trios.isEmpty()) {
                logger.info("No PED file passed or no *non-skipped* trios found in PED file. Skipping family priors.");
                skipFamilyPriors = true;
            }
        }

        final VCFHeader header = vcfHeaders.values().iterator().next();
        if ( ! header.hasGenotypingData() ) {
            throw new UserException("VCF has no genotypes");
        }

        if ( header.hasInfoLine(GATKVCFConstants.MLE_ALLELE_COUNT_KEY) ) {
            final VCFInfoHeaderLine mleLine = header.getInfoHeaderLine(GATKVCFConstants.MLE_ALLELE_COUNT_KEY);
            if ( mleLine.getCountType() != VCFHeaderLineCount.A ) {
                throw new UserException("VCF does not have a properly formatted MLEAC field: the count type should be \"A\"");
            }

            if ( mleLine.getType() != VCFHeaderLineType.Integer ) {
                throw new UserException("VCF does not have a properly formatted MLEAC field: the field type should be \"Integer\"");
            }
        }

        // Initialize VCF header
        final Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfHeaders.values(), true);
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.GENOTYPE_PRIOR_KEY));
        if (!skipFamilyPriors) {
            headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.JOINT_LIKELIHOOD_TAG_NAME));
            headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.JOINT_POSTERIOR_TAG_NAME));
        }
        headerLines.addAll(getDefaultToolVCFHeaderLines());

        vcfWriter.writeHeader(new VCFHeader(headerLines, vcfSamples));

        final Map<String,Set<Sample>> families = sampleDB.getFamilies(vcfSamples);
        famUtils = new FamilyLikelihoods(sampleDB, deNovoPrior, vcfSamples, families);
    }

    /**
     * Entry-point function to initialize the samples database from input data
     */
    private SampleDB initializeSampleDB() {
        final SampleDBBuilder sampleDBBuilder = new SampleDBBuilder(PedigreeValidationType.STRICT);
        if (pedigreeFile != null) {
            sampleDBBuilder.addSamplesFromPedigreeFiles(Collections.singletonList(pedigreeFile));
        }
        return sampleDBBuilder.getFinalSampleDB();
    }

    @Override
    public void apply(final VariantContext variant,
                      final ReadsContext readsContext,
                      final ReferenceContext referenceContext,
                      final FeatureContext featureContext) {
        final Collection<VariantContext> vcs = featureContext.getValues(getDrivingVariantsFeatureInput());

        final Collection<VariantContext> otherVCs = featureContext.getValues(supportVariants);

        final int missing = supportVariants.size() - otherVCs.size();

        for ( final VariantContext vc : vcs ) {
            final VariantContext vc_familyPriors;
            final VariantContext vc_bothPriors;

            //do family priors first (if applicable)
            final VariantContextBuilder builder = new VariantContextBuilder(vc);
            //only compute family priors for biallelelic sites
            if (!skipFamilyPriors && vc.isBiallelic()){
                final GenotypesContext gc = famUtils.calculatePosteriorGLs(vc);
                builder.genotypes(gc);
            }
            VariantContextUtils.calculateChromosomeCounts(builder, false);
            vc_familyPriors = builder.make();

            if (!skipPopulationPriors) {
                vc_bothPriors = PosteriorProbabilitiesUtils.calculatePosteriorProbs(vc_familyPriors, otherVCs, missing * numRefIfMissing, globalPrior, !ignoreInputSamples, defaultToAC, useACoff);
            } else {
                final VariantContextBuilder builder2 = new VariantContextBuilder(vc_familyPriors);
                VariantContextUtils.calculateChromosomeCounts(builder, false);
                vc_bothPriors = builder2.make();
            }
            vcfWriter.add(vc_bothPriors);
        }
    }

    @Override
    public void closeTool(){
        vcfWriter.close();
    }
}

