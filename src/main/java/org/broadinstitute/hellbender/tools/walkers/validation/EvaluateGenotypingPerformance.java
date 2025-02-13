package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.util.ParsingUtils;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFUtils;
import org.apache.commons.collections4.Predicate;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.AbstractConcordanceWalker;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.broadinstitute.hellbender.engine.filters.VariantFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.VariantIDsVariantFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Evaluate genotyping accuracy against truth.  Both truth and eval vcfs may be single or multi-sample.  Produces two output files:
 * <ol>
 *     <li>Pearson correlations as a function of alt allele frequency.  Sites are binned by alt allele fraction.  Pearson's r is calculated between the genotypes (or other specified annotation) in the truth and eval data, for each bin.</li>
 *     <li>Accuracy metrics.  True positives, false positive, true negatives, and false negatives are counted for each sample, and recall, precision, and accuracy are calculated.  These are not stratified by alt allele frequency.</li>
 * </ol>
 * All results are stratified by sample and SNP or INDEL.
 *
 * <h3>Usage examples</h3>
 *
 * <h4>Sub-population allele frequencies for sample1 are given by --af-annotation.</h4>
 * <pre>
 *     gatk EvaluateGenotypingPerformance
 *          -eval eval.vcf
 *          -truth truth.vcf
 *          --af-annotation sample1:af1
 *          -O output.correlations.tsv
 *          -OA output.accuracy.tsv
 * </pre>
 *
 * <h4>Individuals from different subpopulations should have different annotations containing the sub-population allele frequencies for their ethnicity.</h4>
 * <pre>
 *     gatk EvaluateGenotypingPerformance
 *          -eval eval.vcf
 *          -truth truth.vcf
 *          --af-annotation sample1:af1
 *          --af-annotation sample2:af2
 *          -O output.correlations.tsv
 *          -OA output.accuracy.tsv
 * </pre>
 *
 * <h4>If sample names in eval and truth differ, they can be mapped to one another using --sample-map.</h4>
 * <pre>
 *     gatk EvaluateGenotypingPerformance
 *          -eval eval.vcf
 *          -truth truth.vcf
 *          --af-annotation sample1_eval:af1
 *          --af-annotation sample2_eval:af2
 *          --sample-map sample1_eval:sample1_truth
 *          --sample-map sample2_eval:sample2_truth
 *          -O output.correlations.tsv
 *          -OA output.accuracy.tsv
 * </pre>
 *
 * <h4>If sample names differ, and both truth and eval contain only a single sample, comparison can be forced between the individual samples regardless of name using --force-compare-single-sample</h4>
 *  <pre>
 *     gatk EvaluateGenotypingPerformance
 *          -eval eval.vcf
 *          -truth truth.vcf
 *          --af-annotation sample1:af1
 *          --force-compare-single-sample
 *          -O output.correlations.tsv
 *          -OA output.accuracy.tsv
 * </pre>
 *
 * <h4>Allele fractions can be given in a resource vcf instead of the eval vcf using --resource.</h4>
 *  <pre>
 *     gatk EvaluateGenotypingPerformance
 *          -eval eval.vcf
 *          -truth truth.vcf
 *          --af-annotation sample1:af1
 *          --resource af_resource.vcf
 *          -O output.correlations.tsv
 *          -OA output.accuracy.tsv
 * </pre>
 */
@CommandLineProgramProperties(
        summary = EvaluateGenotypingPerformance.USAGE_SUMMARY,
        oneLineSummary = EvaluateGenotypingPerformance.USAGE_ONE_LINE_SUMMARY,
        programGroup = VariantEvaluationProgramGroup.class
)
@DocumentedFeature
public class EvaluateGenotypingPerformance extends AbstractConcordanceWalker {
    final static String USAGE_SUMMARY = "Evaluate genotyping in a (potentially multisample) vcf against truth genotypes.  This is intended generally for imputed data, where we are evaluating genotype accuracy, not full variant calling accuracy." +
            " Outputs are a table of Pearson correlations vs allele frequency, and well as a table of sensitivity, precision, and accuracy.  Note that this tool evaluates only genotyping accuracy.  That is, if a variant is called in the truth" +
            " but there is no corresponding record in the evaluation vcf, there is no contribution to the result.  If, however, there is a corresponding record in the evaluation vcf, which is genotyped as a hom-ref, while the truth has " +
            "a het, this site would be counted as a false negative.  This is generally useful for evaluating array or imputed data, where the objective is to genotype known sites, but is generally not desirable for sequencing variant calling.";
    final static String USAGE_ONE_LINE_SUMMARY = "Evaluate genotyping in a (potentially multisample) vcf against truth genotypes.";

    /**
     * A resource file from which to extract alt allele frequencies.  If not included, then allele frequencies should be annotated in the vcf being evaluated.
     */
    @Argument(fullName= StandardArgumentDefinitions.RESOURCE_LONG_NAME, doc="External resource VCF file containing alt allele frequencies.", optional=true)
    private FeatureInput<VariantContext> af_resource;

    /**
     * A resource file from which to extract alt allele counts in the reference panel.
     */
    @Argument(fullName = "ac-resource", doc="External resource VCF file containing alt allele counts from reference panel.", optional=true)
    private FeatureInput<VariantContext> ac_resource;

    /**
     * max ac for accuracy stratification
     */
    @Argument(fullName = "max-ac-bin", optional = true)
    private int maxACBin = 20;
    /**
     * List of IDs (or a .list file containing ids) to select. The tool will only select variants whose ID
     * field is present in this list of IDs. The matching is done by exact string matching. If a file, the file
     * name must end in ".list", and the expected file format is simply plain text with one ID per line.
     */
    @Argument(fullName="keep-ids", shortName="ids", doc="List of variant rsIDs to select", optional=true)
    private Set<String> rsIDsToKeep = new HashSet<>();

    /**
     * Genotype annotation to use for correlation calculation.  Annotation to use for calculating correlation between truth and eval data, instead of genotype.
     * If a genotype does not include this annotation, the the genotype call will be converted into a dosage (this is likely to happen if the truth does not include a dosage field).
     * Usually this will be "DS", if used.
     */
    @Argument(fullName="dosage-field", doc="Genotype annotation to use for correlation calculation. Generally this will be 'DS', if used.", optional=true)
    String dosageField = null;

    @Argument(
            doc = "Output file for correlation.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private GATKPath outputFile = null;

    @Argument(
            doc = "Output file for accuracy.",
            fullName = "output-accuracy",
            shortName = "OA"
    )
    private GATKPath outputAccuracyFile = null;

    @Argument(
            doc = "Output file for accuracy by AF",
            fullName = "output-accuracy-af"
    )
    private GATKPath outputAccuracyAFFile = null;

    @Argument(
            doc = "Output file for accuracy by AC",
            fullName = "output-accuracy-ac",
            optional = true
    )
    private GATKPath outputAccuracyACFile = null;

    @Argument(
            doc = "Output file for GP calibration",
            fullName = "output-gp-calibration",
            optional = true
    )
    private GATKPath outputCalibrationFile = null;

    @Argument(
            fullName = "n-calibration-bins",
            optional = true
    )
    private int nCalibrationBins = 10;

    /**
     * List of strings specifying the annotations from which alt allele frequencies for each sample should be read.  Each string should be the sample name, and corresponding annotation, separated by a colon.
     */
    @Argument(
        fullName = "af-annotations",
        doc="List of mappings from sample name to corresponding annotation storing alt allele frequencies.  Mappings use colon as separator."
    )
    private List<String> afAnnotations;

    /**
     * List of mappings between sample names in the evaluation vcf and in the truth vcf.  Corresponding sample names should be separated by colons.  If no mappings is passed for a particular sample, it is
     * assumed that the same sample name is also in the truth vcf.
     */
    @Argument(
            fullName = "sample-map",
            doc = "List of mappings between corresponding eval and truth sample names.  Mappings use colon as separator.",
            optional = true
    )
    private List<String> sampleMappings;


    /**
     * Number of bins to bin sites by alt allele frequency for correlation calculation.  (nBins - 1) bins will be spaced equally in log space between firstBinRightEdge and 1.  The first bin will span 0 to firstBinRightEdge.
     */
    @Argument(
            fullName = "nbins",
            optional = true,
            doc = "Number of bins to bin sites by alt allele frequency for correlation calculation.",
            minValue = 1
    )
    int nBins = 14;

    @Argument(
            fullName = "first-bin-right-edge",
            optional = true,
            doc = "The right edge of the first bin for binning sites by alt allele frequency for correlation calculation.",
            minValue = 0,
            maxValue = 1
    )
    double firstBinRightEdge = 0.0005;

    @Argument(
            fullName = "min-af-for-accuracy-metrics",
            optional = true,
            doc = "Minimum alt allele frequency to include a site in accuracy calculations. Sites with lower alt allele frequency are not included.",
            minValue = 0,
            maxValue = 1
    )
    double minAfForAccuracyMetrics = 0;

    @Argument(
            fullName = "allow-differing-ploidies",
            optional = true,
            doc = "Allow ploidy to differ at a site between truth and eval.  This should generally only occur if chromosome X is being treated differently in males.  If this option is true " +
                    "and different ploidies are found at a site, the genotypes are considered to agree if the set of observed alleles is the same in both the truth and eval."
    )
    boolean allowDifferingPloidy = false;

    @Argument(
            fullName = "force-compare-single-sample",
            optional = true,
            doc = "Compare single samples in truth and eval vcfs, regardless of sample names.  Both vcfs must contain only one sample, or else error will be thrown.  " +
                    "Cannot be used with --sample-map argument."
    )
    boolean forceCompareSingleSample = false;

    final List<List<AFCorrelationAggregator>> aggregators = new ArrayList<>();
    final List<AccuracyMetrics> snpMetrics = new ArrayList<>();
    final List<AccuracyMetrics> indelMetrics = new ArrayList<>();
    final List<List<AFAccuracyMetrics>> afSnpMetrics = new ArrayList<>();
    final List<List<AFAccuracyMetrics>> afIndelMetrics = new ArrayList<>();
    final List<List<ACAccuracyMetrics>> acSnpMetrics = new ArrayList<>();
    final List<List<ACAccuracyMetrics>> acIndelMetrics = new ArrayList<>();
    final Map<String, String> afAnnotationsMap = new HashMap<>();
    final Set<String> afAnnotationSet = new HashSet<>();
    final Map<String, String> sampleMap = new HashMap<>();
    final List<List<List<GPCalibrationMetrics>>> snpCalibrations = new ArrayList<>();
    final List<List<List<GPCalibrationMetrics>>> indelCalibrations = new ArrayList<>();

    final List<GenotypeType> genotypeTypes = Arrays.asList(GenotypeType.HOM_REF, GenotypeType.HET, GenotypeType.HOM_VAR);
    @Override
    public void onTraversalStart() {
        if (forceCompareSingleSample && sampleMappings.size() > 0) {
            throw new GATKException("Cannot use --force-compare-single-sample and --sample-mappings arguments together");
        }

        final VCFHeader header = getEvalHeader();
        final VCFHeader truthHeader = getTruthHeader();
        final List<String> samples = header.getSampleNamesInOrder();
        final Set<String> truthSamples = new HashSet<>(truthHeader.getSampleNamesInOrder());
        loadMapping(afAnnotations, afAnnotationsMap);
        //create a set of the unique annotations that will need to be extracted for each site
        afAnnotationSet.addAll(afAnnotationsMap.values());

        if (sampleMappings.size() > 0) {
            loadMapping(sampleMappings, sampleMap);
        }

        if (forceCompareSingleSample) {
            if (samples.size() != 1 || truthSamples.size() != 1) {
                throw new GATKException("When using --force-compare-single-sample buth eval and truth must contain only one sample.  Eval contains " + samples.size() + " samples" +
                        " and truth contains " + truthSamples.size() + " samples");
            }
            sampleMap.put(samples.get(0), truthSamples.iterator().next());
        }

        //initialize the pearson correlation aggregators and accuracy metrics for each sample

        for (final String sample : samples) {
            final String mappedSample = getMappedSample(sample);
            if ((!afAnnotationsMap.containsKey(sample) && !afAnnotationsMap.containsKey(mappedSample))) {
                //samples which are not included in the annotations map are ignored.
                logger.warn("Sample " + sample + " which maps to " + mappedSample + " does not have corresponding afAnnotation given.  " + sample + " will not be included in results.");
                continue;
            }

            if (!truthSamples.contains(mappedSample)) {
                //samples which are not matched to a sample in the truth vcf are ignored
                logger.warn("Sample " + sample + " maps to " + mappedSample + ", but " + mappedSample + " is not found in truth vcf.  " + sample + " will not be included in results.");
                continue;
            }
            //each sample has a list of correlation, snp & indel accuracy aggregators, one for each af bin
            final List<AFCorrelationAggregator> theseAggregators = new ArrayList<>();
            final List<AFAccuracyMetrics> theseSnpMetrics = new ArrayList<>();
            final List<AFAccuracyMetrics> theseIndelMetrics = new ArrayList<>();
            theseAggregators.add(new AFCorrelationAggregator(firstBinRightEdge / 2, sample));
            theseSnpMetrics.add(new AFAccuracyMetrics(firstBinRightEdge / 2, VariantContext.Type.SNP, sample));
            theseIndelMetrics.add(new AFAccuracyMetrics(firstBinRightEdge / 2, VariantContext.Type.INDEL, sample));
            for (int i = 1; i < nBins; i++) {
                //bins are spaced evenly in log space
                final double log10_bin_width = (-Math.log10(firstBinRightEdge)) / (double) (nBins - 1);
                final double logBinCenter = Math.log10(firstBinRightEdge) + (i - 0.5) * log10_bin_width;
                final double binCenter = Math.pow(10, logBinCenter);
                theseAggregators.add(new AFCorrelationAggregator(binCenter, sample));
                theseSnpMetrics.add(new AFAccuracyMetrics(binCenter, VariantContext.Type.SNP, sample));
                theseIndelMetrics.add(new AFAccuracyMetrics(binCenter, VariantContext.Type.INDEL, sample));
            }
            final List<ACAccuracyMetrics> theseACSnpMetrics = new ArrayList<>();
            final List<ACAccuracyMetrics> theseACIndelMetrics = new ArrayList<>();
            for (int i = 0; i < maxACBin; i++) {
                theseACSnpMetrics.add(new ACAccuracyMetrics(i+1, VariantContext.Type.SNP, sample));
                theseACIndelMetrics.add(new ACAccuracyMetrics(i+1, VariantContext.Type.INDEL, sample));
            }

            final List<List<GPCalibrationMetrics>> theseSnpGPCalibrationMetrics = new ArrayList<>();
            final List<List<GPCalibrationMetrics>> theseIndelGPCalibrationMetrics = new ArrayList<>();
            for (int i =0; i < nCalibrationBins; i++) {
                final List<GPCalibrationMetrics> thisBinSnpCalibrationMetrics = new ArrayList<>();
                final List<GPCalibrationMetrics> thisBinIndelCalibrationMetrics = new ArrayList<>();
                for (final GenotypeType genotypeType : genotypeTypes) {
                    thisBinSnpCalibrationMetrics.add(new GPCalibrationMetrics(sample, (i + 0.5) / (float) nCalibrationBins, VariantContext.Type.SNP, genotypeType));
                    thisBinIndelCalibrationMetrics.add(new GPCalibrationMetrics(sample, (i + 0.5) / (float) nCalibrationBins, VariantContext.Type.INDEL, genotypeType));
                }
                theseSnpGPCalibrationMetrics.add(thisBinSnpCalibrationMetrics);
                theseIndelGPCalibrationMetrics.add(thisBinIndelCalibrationMetrics);
            }
            aggregators.add(theseAggregators);
            snpMetrics.add(new AccuracyMetrics(VariantContext.Type.SNP, sample));
            indelMetrics.add(new AccuracyMetrics(VariantContext.Type.INDEL, sample));
            afSnpMetrics.add(theseSnpMetrics);
            afIndelMetrics.add(theseIndelMetrics);
            acSnpMetrics.add(theseACSnpMetrics);
            acIndelMetrics.add(theseACIndelMetrics);
            snpCalibrations.add(theseSnpGPCalibrationMetrics);
            indelCalibrations.add(theseIndelGPCalibrationMetrics);
        }

        if (aggregators.size() == 0) {
            throw new GATKException("There are no valid comparisons to preform, there is likely a mistake in sample or af annotation mappings.");
        }
    }

    /**
     * Load mappings in mappingStrings in map
     * @param mappingStrings list of string which map keys to values
     * @param map map into which the mappings are placed
     */
    private void loadMapping(final List<String> mappingStrings, final Map<String, String> map) {
        for (final String mappingString : mappingStrings) {
            String[] tokens = new String[2];
            StringUtil.split(mappingString, tokens,':');
            map.put(tokens[0], tokens[1]);
        }
    }

    /**
     * get corresponding truth sample name for a particular eval sample.
     */
    private String getMappedSample(final String sample) {
        return sampleMap.getOrDefault(sample, sample);
    }

    /**
     * get the annotation from which to extract alt allele frequency information for a particular sample
     */
    private String getAfAnnotation(final String sample) {
        String afAnnotation = afAnnotationsMap.get(sample);
        if (afAnnotation == null) {
            afAnnotation = afAnnotationsMap.get(getMappedSample(sample));
        }

        return afAnnotation;
    }

    @Override
    protected Predicate<VariantContext> makeTruthVariantFilter() {
        return VariantFilterLibrary.ALLOW_ALL_VARIANTS::test;
    }

    @Override
    protected Predicate<VariantContext> makeEvalVariantFilter() {
        final VariantFilter filter;
        if (rsIDsToKeep != null && !rsIDsToKeep.isEmpty()) {
            filter = new VariantFilterLibrary.PassesFiltersVariantFilter().and(new VariantIDsVariantFilter(rsIDsToKeep));
        } else {
            filter = new VariantFilterLibrary.PassesFiltersVariantFilter();
        }
        return filter::test;
    }

    @Override
    protected boolean shouldVariantsBeMatched(final VariantContext truth, final VariantContext eval) {
        return truth.getAlleles().containsAll(eval.getAlleles());
    }

    @Override
    protected void apply(final TruthVersusEval truthVersusEval, final ReadsContext readsContext, final ReferenceContext refContext) {

        if (!truthVersusEval.hasEval()) {
            //only include sites which have both an eval record
            return;
        }

        final VariantContext evalVC = truthVersusEval.getEval();

        final VariantContext truthVC = truthVersusEval.getTruth();

        //since multiple samples are likely to share the same afAnnotation, we build a map to store the values for each afAnnotation, instead of extracting the value for each sample
        final Map<String, Double> afMap = buildAFMap(evalVC);
        final Integer ac = getRefPanelAC(evalVC);

        for (int i=0; i<aggregators.size(); i++) {
            //loop through samples
            final String sample = snpMetrics.get(i).SAMPLE;
            final String mappedSample = getMappedSample(sample);
            final String afAnnotation = getAfAnnotation(sample);
            final Double af = afMap.get(afAnnotation);
            if (outputCalibrationFile == null && af == null && ac == null) {
                // if we don't have alt allele frac information for this site, skip it
                continue;
            }
            final Integer bin = af!=null? getBin(af) : null;
            final Genotype evalGenotype = evalVC.getGenotype(sample);
            final Genotype truthGenotype = truthVC != null? truthVC.getGenotype(mappedSample) : null;
            if ((truthGenotype != null && truthGenotype.isNoCall()) || (evalGenotype != null && evalGenotype.isNoCall())) {
                //not confident in truth, or eval is a noCall, so skip
                continue;
            }
            final double evalDosageFrac = getDosageFrac(evalGenotype, evalVC.getReference(), dosageField);
            final double truthDosageFrac = getDosageFrac(truthGenotype, truthVC != null ? truthVC.getReference() : null, dosageField);
            final ConcordanceState concordanceState = getConcordanceState(truthGenotype, evalGenotype, evalVC.isFiltered());
            if (evalVC.isSNP()) {
                if (bin != null) {
                    aggregators.get(i).get(bin).snp_pearsonCorrelationAggregator.addEntry(evalDosageFrac, truthDosageFrac);
                    afSnpMetrics.get(i).get(bin).incrementMetrics(concordanceState);
                    if (af >= minAfForAccuracyMetrics) {
                        snpMetrics.get(i).incrementMetrics(concordanceState);
                    }
                }

                if (ac != null) {
                    final Integer acBin = ac < maxACBin ? ac - 1 : maxACBin -1;
                    acSnpMetrics.get(i).get(acBin).incrementMetrics(concordanceState);
                }

                // calibration
                if (evalGenotype.hasExtendedAttribute("GP")) {
                    final List<String> gpStrings = ParsingUtils.split((String)evalGenotype.getExtendedAttribute("GP"), ',');
                    if(gpStrings.size()== 3) {
                        for (int j=0; j<3; j++) {
                            final double gp = VCFUtils.parseVcfDouble(gpStrings.get(j));
                            final int gp_bin = gp == 1? nCalibrationBins - 1 : (int) (gp * nCalibrationBins);
                            final GenotypeType truthGenotypeType = truthGenotype == null? GenotypeType.HOM_REF : truthGenotype.getType();
                            final GenotypeType thisGenotypeType = genotypeTypes.get(j);
                            snpCalibrations.get(i).get(gp_bin).get(j).increment(thisGenotypeType == truthGenotypeType);
                        }
                    }
                }

            } else if (evalVC.isIndel()) {
                if (bin != null) {
                    aggregators.get(i).get(bin).indel_pearsonCorrelationAggregator.addEntry(evalDosageFrac, truthDosageFrac);
                    afIndelMetrics.get(i).get(bin).incrementMetrics(concordanceState);
                    if (af >= minAfForAccuracyMetrics) {
                        indelMetrics.get(i).incrementMetrics(concordanceState);
                    }
                }

                if (ac != null) {
                    final Integer acBin = ac < maxACBin ? ac - 1 : maxACBin -1;
                    acIndelMetrics.get(i).get(acBin).incrementMetrics(concordanceState);
                }

                // calibration
                if (evalGenotype.hasExtendedAttribute("GP")) {
                    final List<String> gpStrings = ParsingUtils.split((String)evalGenotype.getExtendedAttribute("GP"), ',');
                    if(gpStrings.size()== 3) {
                        for (int j=0; j<3; j++) {
                            final double gp = VCFUtils.parseVcfDouble(gpStrings.get(j));
                            final int gp_bin = gp == 1? nCalibrationBins - 1 : (int) (gp * nCalibrationBins);
                            final GenotypeType truthGenotypeType = truthGenotype == null? GenotypeType.HOM_REF : truthGenotype.getType();
                            final GenotypeType thisGenotypeType = genotypeTypes.get(j);
                            indelCalibrations.get(i).get(gp_bin).get(j).increment(thisGenotypeType == truthGenotypeType);
                        }
                    }
                }
            }
        }
    }

    /**
     * builds a map storing the alt allele frequency values for all required af annotations for the given eval variant
     * @param vc eval variant context
     * @return map from annotation to extracted alt allele context
     */
    private Map<String, Double> buildAFMap(final VariantContext vc) {
        if (af_resource != null) {
            final FeatureContext featureContext = new FeatureContext(features, new SimpleInterval(vc));
            final List<VariantContext> resourceFeatures = featureContext.getValues(af_resource, vc.getStart());
            for (final VariantContext feature : resourceFeatures) {
                //find the overlapping resource feature which starts at the same position and has the same reference allele as the eval vc
                if (feature.getReference().equals(vc.getReference())) {
                    //find the index in the resource feature which corresponds to the first alt allele of the eval vc
                    final int altAlleleIndex = feature.getAlleleIndex(vc.getAlternateAllele(0));
                    if (altAlleleIndex > 0) {
                        return buildAFMapForIndex(feature, afAnnotationSet,altAlleleIndex - 1);
                    }
                }
            }
            //haven't found a feature to extract afs from, so return empty map
            return Collections.emptyMap();
        } else {
            //if no allele frequency resource has been passed, annotations should be stored in eval vc
            return buildAFMapForIndex(vc, afAnnotationSet, 0);
        }
    }

    private Integer getRefPanelAC(final VariantContext vc) {
        if (ac_resource != null) {
            final FeatureContext featureContext = new FeatureContext(features, new SimpleInterval(vc));
            final List<VariantContext> resourceFeatures = featureContext.getValues(ac_resource, vc.getStart());
            for (final VariantContext feature : resourceFeatures) {
                if (feature.getReference().equals(vc.getReference())) {
                    final int altAlleleIndex = feature.getAlleleIndex(vc.getAlternateAllele(0));
                    if (altAlleleIndex > 0) {
                        final List<Integer> ac_list = feature.getAttributeAsIntList("AC", 0);
                        return ac_list.get(altAlleleIndex - 1);
                    }
                }
            }
        }
        return null;
    }

    /**
     * build a map storing the alt allele frequency values for all required af annotations from the given variant
     * @param vc variant context containing af annotations
     * @param index this index of the alt allele frequency to extract
     * @return map from annotation to extracted alt allele context
     */
    static Map<String, Double> buildAFMapForIndex(final VariantContext vc, final Collection<String> annotations, final int index) {
        final Map<String, Double> afMap = new HashMap<>();
        for (final String afAnnotation : annotations) {
            final List<Double> afList = vc.getAttributeAsDoubleList(afAnnotation, -99);
            if (afList.size() < index + 1) {
                continue;
            }
            final Double af = afList.get(index);
            if (af < 0 || af > 1) {
                throw new GATKException("Invalid AF value " + af + " at " + vc.getContig() + ":" + vc.getStart() + " for allele " + vc.getAlternateAllele(index-1));
            }
            afMap.put(afAnnotation, af);
        }

        return afMap;
    }


    /**
     * Get the dosage frequency from the given genotype.
     */
     private double getDosageFrac(final Genotype geno, final Allele refAllele, final String dosageField) {
         if (geno == null) {
             return 0;
         }
         if (dosageField != null) {
             final String dosageString = (String)geno.getExtendedAttribute(dosageField);
             return dosageString != null ? VCFUtils.parseVcfDouble(dosageString)/(double)geno.getPloidy() : getDosageFracFromGenotypeCall(geno, refAllele);
         }

         return getDosageFracFromGenotypeCall(geno, refAllele);
    }

    private double getDosageFracFromGenotypeCall(final Genotype geno, final Allele refAllele) {
        return 1d - (double)geno.countAllele(refAllele)/(double)geno.getPloidy();
    }

    @Override
    protected boolean genotypesAgree(final Genotype geno1, final Genotype geno2) {
        if (geno1 != null && geno1.getPloidy() != geno2.getPloidy()) {
            //situation could arise from haploid vs diploid on X
            if (!allowDifferingPloidy) {
                throw new GATKException("sample " + geno1.getSampleName() + " is ploidy " + geno1.getPloidy() + " while truth sample " + geno2.getSampleName() + " is ploidy " + geno2.getPloidy() + "." +
                        "  This may be due to haploid vs diploid representation on X.  If you would like to allow for this TYPE of data, use the allowDifferingPloidy argument.");
            }

            //build a set of all alleles in truth and eval genotypes.  If sets are the same, then genotypes "agree", if different, then genotypes "agree".  so 0/0 and 0 agree, 1/1 and 1 agree, but 0/1 and 1 disagree.
            final Set<Allele> truthAlleles = new HashSet<>(geno1.getAlleles());
            final Set<Allele> evalAlleles = new HashSet<>(geno2.getAlleles());

            return truthAlleles.equals(evalAlleles);
        } else {
            return super.genotypesAgree(geno1, geno2);
        }
    }

    @Override
    protected ConcordanceState getConcordanceState(final Genotype truth, final Genotype eval, final boolean evalWasFiltered) {
        if (truth != null && truth.getPloidy() != eval.getPloidy()) {
            //situation could arise from haploid vs diploid on X
            if (!allowDifferingPloidy) {
                throw new GATKException("sample " + eval.getSampleName() + " is ploidy " + eval.getPloidy() + " while truth sample " + truth.getSampleName() + " is ploidy " + truth.getPloidy() + "." +
                        "  This may be due to haploid vs diploid representation on X.  If you would like to allow for this TYPE of data, use the allowDifferingPloidy argument.");
            }

            //build a set of all alleles in truth and eval genotypes.  If sets are the same, then genotypes "agree", if different, then genotypes "agree".  so 0/0 and 0 agree, 1/1 and 1 agree, but 0/1 and 1 disagree.
            final Set<Allele> truthAlleles = new HashSet<>(truth.getAlleles());
            final Set<Allele> evalAlleles = new HashSet<>(eval.getAlleles());

            final boolean isPositiveEval = isPositive(eval);
            final boolean isPositiveTruth = isPositive(truth);
            final boolean genotypesAgree = truthAlleles.equals(evalAlleles);

            return evaluateConcordanceState(isPositiveEval, isPositiveTruth, genotypesAgree, evalWasFiltered);
        }

        //now normal situation, agree on ploidy.
        return super.getConcordanceState(truth, eval, evalWasFiltered);
    }



    @Override
    public Object onTraversalSuccess() {
        final MetricsFile<AFCorrelationMetric, Integer> correlationWriter = getMetricsFile();
        for (final List<AFCorrelationAggregator> theseAggregators : aggregators ) {
            for (final AFCorrelationAggregator aggregator : theseAggregators) {
                correlationWriter.addMetric(new AFCorrelationMetric(aggregator));
            }
        }
        correlationWriter.write(outputFile.toPath().toFile());

        //write out accuracy results
        final MetricsFile<AccuracyMetrics, Integer> accuracyWriter = getMetricsFile();
        snpMetrics.forEach(AccuracyMetrics::calculateDerivedMetrics);
        indelMetrics.forEach(AccuracyMetrics::calculateDerivedMetrics);
        accuracyWriter.addAllMetrics(snpMetrics);
        accuracyWriter.addAllMetrics(indelMetrics);
        accuracyWriter.write(outputAccuracyFile.toPath().toFile());

        final MetricsFile<AFAccuracyMetrics, Integer> afAccuracryWriter = getMetricsFile();
        for (final List<AFAccuracyMetrics> metrics : afSnpMetrics) {
            metrics.forEach(AccuracyMetrics::calculateDerivedMetrics);
            afAccuracryWriter.addAllMetrics(metrics);
        }
        for (final List<AFAccuracyMetrics> metrics : afIndelMetrics) {
            metrics.forEach(AccuracyMetrics::calculateDerivedMetrics);
            afAccuracryWriter.addAllMetrics(metrics);
        }
        afAccuracryWriter.write(outputAccuracyAFFile.toPath().toFile());

        if (outputAccuracyACFile != null) {
            final MetricsFile<ACAccuracyMetrics, Integer> acAccuracyMetricsWriter = getMetricsFile();
            for (final List<ACAccuracyMetrics> metrics : acSnpMetrics) {
                metrics.forEach(AccuracyMetrics::calculateDerivedMetrics);
                acAccuracyMetricsWriter.addAllMetrics(metrics);
            }
            for (final List<ACAccuracyMetrics> metrics : acIndelMetrics) {
                metrics.forEach(AccuracyMetrics::calculateDerivedMetrics);
                acAccuracyMetricsWriter.addAllMetrics(metrics);
            }
            acAccuracyMetricsWriter.write(outputAccuracyACFile.toPath().toFile());
        }

        if (outputCalibrationFile != null) {
            final MetricsFile<GPCalibrationMetrics, Integer> calibrationMetricsWriter = getMetricsFile();
            for (final List<List<GPCalibrationMetrics>> metricsList : snpCalibrations) {
                for (final List<GPCalibrationMetrics> metrics : metricsList) {
                    calibrationMetricsWriter.addAllMetrics(metrics);
                }
            }
            for (final List<List<GPCalibrationMetrics>> metricsList : indelCalibrations) {
                for (final List<GPCalibrationMetrics> metrics : metricsList) {
                    calibrationMetricsWriter.addAllMetrics(metrics);
                }
            }
            calibrationMetricsWriter.write(outputCalibrationFile.toPath().toFile());
        }
        return null;
    }

    /**
     * get the af bin to assign a particular site to based on its allele frequency
     */
    int getBin(double af) {
        if (af < firstBinRightEdge) {
            return 0;
        }

        return (int)Math.ceil((1-Math.log10(af)/Math.log10(firstBinRightEdge)) * (nBins - 1));
    }

    /**
     * A class for holding Pearson correlation aggregators for a particular sample and af bin.
     */
    private static final class AFCorrelationAggregator {
        final double binCenter;
        final String sampleName;
        final PearsonCorrelationAggregator snp_pearsonCorrelationAggregator = new PearsonCorrelationAggregator();
        final PearsonCorrelationAggregator indel_pearsonCorrelationAggregator = new PearsonCorrelationAggregator();

        AFCorrelationAggregator(final double binCenter, final String sampleName) {
            this.binCenter = binCenter;
            this.sampleName = sampleName;
        }
    }

     public static final class AFCorrelationMetric extends MetricBase {
        public double BIN_CENTER;
        public String SAMPLE;
        public double SNP_CORRELATION;
        public double SNP_SITES;
        public double INDEL_CORRELATION;
        public double INDEL_SITES;

        AFCorrelationMetric(final AFCorrelationAggregator afCorrelationAggregator) {
            BIN_CENTER = afCorrelationAggregator.binCenter;
            SAMPLE = afCorrelationAggregator.sampleName;
            SNP_CORRELATION = afCorrelationAggregator.snp_pearsonCorrelationAggregator.getCorrelation();
            SNP_SITES = afCorrelationAggregator.snp_pearsonCorrelationAggregator.n;
            INDEL_CORRELATION = afCorrelationAggregator.indel_pearsonCorrelationAggregator.getCorrelation();
            INDEL_SITES = afCorrelationAggregator.indel_pearsonCorrelationAggregator.n;
        }
    }

    /**
     * A class used to calculate Pearson correlation in a single pass with minimal memory usage
     */
    static final class PearsonCorrelationAggregator {
        /*
        We will calculate r using the formula r = (<xy> - <x><y>)/(sqrt(<x^2> - <x>^2) * sqrt(<y^2> - <y>^2)).  So we need to store
        1) the sums of:
            xy
            x
            y
            x^2
            y^2
        2) the total number of entries

        From these we can calculate all the required expectation values.
         */
        private double sum_xy;
        private double sum_x;
        private double sum_y;
        private double sum_x2;
        private double sum_y2;
        private int n;

        void addEntry(final double x, final double y) {
            sum_xy += x*y;
            sum_x += x;
            sum_y += y;
            sum_x2 += x*x;
            sum_y2 += y*y;
            n++;
        }

        double getCorrelation() {
            // r = (<xy> - <x><y>)/(sqrt(<x^2> - <x>^2) * sqrt(<y^2> - <y>^2))
            final double n_d = n;
            final double e_xy = sum_xy/n_d;
            final double e_x = sum_x/n_d;
            final double e_y = sum_y/n_d;
            final double e_x2 = sum_x2/n_d;
            final double e_y2 = sum_y2/n_d;

            return (e_xy - e_x * e_y)/(Math.sqrt(e_x2 - e_x*e_x)* Math.sqrt(e_y2 - e_y*e_y));
        }
    }

    public static class GPCalibrationMetrics extends MetricBase {
        public String SAMPLE;
        public double GP_BIN_CENTER;
        public VariantContext.Type TYPE;
        public GenotypeType GTYPE;
        public int CORRECT;
        public int TOTAL;

        GPCalibrationMetrics(final String sample, final double bin_center, final VariantContext.Type type,
                             final GenotypeType gtype) {
            this.SAMPLE = sample;
            this.GP_BIN_CENTER = bin_center;
            this.TYPE = type;
            this.GTYPE = gtype;
        }

        void increment(final boolean isCorrect) {
            TOTAL++;
            if (isCorrect) {
                CORRECT++;
            }
        }
    }

    /**
     * A class to store accuracy metrics
     */
    public static class AccuracyMetrics extends MetricBase {
        public String SAMPLE;
        public VariantContext.Type TYPE;
        public int TRUE_POSITIVES;
        public int TRUE_NEGATIVES;
        public int FALSE_POSITIVES;
        public int FALSE_NEGATIVES;
        public double RECALL;
        public double PRECISION;
        public double ACCURACY;


        AccuracyMetrics(final VariantContext.Type type, final String sampleName) {
            this.TYPE = type;
            this.SAMPLE = sampleName;
        }

        void incrementMetrics(final ConcordanceState concordanceState) {
            switch (concordanceState) {
                case TRUE_POSITIVE:
                    TRUE_POSITIVES++;
                    break;
                case TRUE_NEGATIVE:
                case FILTERED_TRUE_NEGATIVE:
                    TRUE_NEGATIVES++;
                    break;
                case FALSE_POSITIVE:
                    FALSE_POSITIVES++;
                    break;
                case FALSE_NEGATIVE:
                case FILTERED_FALSE_NEGATIVE:
                    FALSE_NEGATIVES++;
            }
        }

        void calculateDerivedMetrics() {
            // recall = TP/(TP + FN)
            RECALL = TRUE_POSITIVES /((double) TRUE_POSITIVES + (double) FALSE_NEGATIVES);

            // precision = TP/(TP + FP)
            PRECISION = TRUE_POSITIVES /((double) TRUE_POSITIVES + (double) FALSE_POSITIVES);

            // accuracy = (TP + FN) / (TP + FP + TN + FN)
            // this is the fraction of sites correctly genotyped
            ACCURACY = (TRUE_POSITIVES + TRUE_NEGATIVES)/(double)(TRUE_POSITIVES + FALSE_POSITIVES + TRUE_NEGATIVES + FALSE_NEGATIVES);
        }
    }

    /**
     * A class to store accuracy metrics for a particular frequency bin
     */
    public static final class AFAccuracyMetrics extends AccuracyMetrics {
        public double BIN_CENTER;

        AFAccuracyMetrics(final double binCenter, final VariantContext.Type type, final String sampleName) {
            super(type, sampleName);
            BIN_CENTER = binCenter;
        }
    }

    /**
     * A class to store accuracy metrics for a particular AC bin
     */
    public static final class ACAccuracyMetrics extends AccuracyMetrics {
        public double AC;

        ACAccuracyMetrics(final double ac, final VariantContext.Type type, final String sampleName) {
            super(type, sampleName);
            AC = ac;
        }
    }

}
