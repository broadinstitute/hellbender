package org.broadinstitute.hellbender.utils.variant;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.utils.ValidationUtils;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.*;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.tuple.MutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.annotator.AnnotationUtils;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.StrandBiasUtils;
import org.broadinstitute.hellbender.tools.walkers.genotyper.AlleleSubsettingUtils;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAlleleCounts;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAssignmentMethod;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypesCache;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.GenotypePriorCalculator;
import org.broadinstitute.hellbender.utils.haplotype.Event;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;

import java.io.Serializable;
import java.nio.file.Path;
import java.util.*;
import java.util.function.BiFunction;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public final class GATKVariantContextUtils {

    private static final Logger logger = LogManager.getLogger(GATKVariantContextUtils.class);

    public static final String MERGE_FILTER_PREFIX = "filterIn";
    public static final String MERGE_REF_IN_ALL = "ReferenceInAlln";
    public static final String MERGE_FILTER_IN_ALL = "FilteredInAll";
    public static final String MERGE_INTERSECTION = "Intersection";

    public static final int DEFAULT_PLOIDY = HomoSapiensConstants.DEFAULT_PLOIDY;

    public static final double SUM_GL_THRESH_NOCALL = -0.1; // if sum(gl) is bigger than this threshold, we treat GL's as non-informative and will force a no-call.

    public static boolean isInformative(final double[] gls) {
        return MathUtils.sum(gls) < GATKVariantContextUtils.SUM_GL_THRESH_NOCALL;
    }

    /**
     * A method that constructs a mutable set of VCF header lines containing the tool name, version, date and command line,
     * and return to the requesting tool.
     * @param toolkitShortName  short name for the tool kit, e.g. "gatk"
     * @param toolName          name of the tool
     * @param versionString     version of the tool
     * @param dateTime          date and time at invocation of the tool
     * @param cmdLine           command line (e.g. options, flags) when the tool is invoked.
     * @return A mutable set of VCF header lines containing the tool name, version, date and command line.
     */
    public static Set<VCFHeaderLine> getDefaultVCFHeaderLines(final String toolkitShortName, final String toolName,
                                                              final String versionString, final String dateTime,
                                                              final String cmdLine) {
        final Set<VCFHeaderLine> defaultVCFHeaderLines = new HashSet<>();
        final Map<String, String> simpleHeaderLineMap = new HashMap<>(4);
        simpleHeaderLineMap.put("ID", toolName);
        simpleHeaderLineMap.put("Version", versionString);
        simpleHeaderLineMap.put("Date", dateTime);
        simpleHeaderLineMap.put("CommandLine", cmdLine);
        defaultVCFHeaderLines.add(new VCFHeaderLine("source", toolName));
        defaultVCFHeaderLines.add(new VCFSimpleHeaderLine(String.format("%sCommandLine", toolkitShortName), simpleHeaderLineMap));
        return defaultVCFHeaderLines;
    }

    /**
     * Creates a VariantContextWriter whose outputFile type is based on the extension of the output file name.
     * The default options set by VariantContextWriter are cleared before applying ALLOW_MISSING_FIELDS_IN_HEADER (if
     * <code>lenientProcessing</code> is set), followed by the set of options specified by any <code>options</code> args.
     *
     * @param outPath output Path for this writer. May not be null.
     * @param referenceDictionary required if on the fly indexing is set, otherwise can be null
     * @param createMD5 true if an md5 file should be created
     * @param options variable length list of additional Options to be set for this writer
     * @returns VariantContextWriter must be closed by the caller
     */
    public static VariantContextWriter createVCFWriter(
            final Path outPath,
            final SAMSequenceDictionary referenceDictionary,
            final boolean createMD5,
            final Options... options)
    {
        Utils.nonNull(outPath);

        VariantContextWriterBuilder vcWriterBuilder =
                new VariantContextWriterBuilder().clearOptions().setOutputPath(outPath);

        if (VariantContextWriterBuilder.OutputType.UNSPECIFIED == VariantContextWriterBuilder.determineOutputTypeFromFile(outPath)) {
            // the only way the user has to specify an output type is by file extension, and htsjdk
            // throws if it can't map the file extension to a known vcf type, so fallback to a default
            // of VCF
            logger.warn(String.format(
                    "Can't determine output variant file format from output file extension \"%s\". Defaulting to VCF.",
                    FilenameUtils.getExtension(outPath.getFileName().toString())));
            vcWriterBuilder = vcWriterBuilder.setOutputFileType(VariantContextWriterBuilder.OutputType.VCF);
        }

        if (createMD5) {
            vcWriterBuilder.setCreateMD5();
        }

        if (null != referenceDictionary) {
            vcWriterBuilder = vcWriterBuilder.setReferenceDictionary(referenceDictionary);
        }

        for (Options opt : options) {
            vcWriterBuilder = vcWriterBuilder.setOption(opt);
        }

        return vcWriterBuilder.build();
    }

    /**
     * Diploid NO_CALL allele list...
     *
     * @deprecated you should use {@link #noCallAlleles(int)} instead. It indicates the presence of a hardcoded diploid assumption which is bad.
     */
    @Deprecated
    public final static List<Allele> NO_CALL_ALLELES = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);

    private static boolean hasPLIncompatibleAlleles(final Collection<Allele> alleleSet1, final Collection<Allele> alleleSet2) {
        final Iterator<Allele> it1 = alleleSet1.iterator();
        final Iterator<Allele> it2 = alleleSet2.iterator();

        while ( it1.hasNext() && it2.hasNext() ) {
            final Allele a1 = it1.next();
            final Allele a2 = it2.next();
            if ( ! a1.equals(a2) )
                return true;
        }

        // by this point, at least one of the iterators is empty.  All of the elements
        // we've compared are equal up until this point.  But it's possible that the
        // sets aren't the same size, which is indicated by the test below.  If they
        // are of the same size, though, the sets are compatible
        return it1.hasNext() || it2.hasNext();
    }

    /**
     *  Does an allele match any in a list of alleles?
     *  We don't assume that ref1/alt1 are in their minimal representation
     * @param ref1
     * @param alt1
     * @param ref2
     * @param altList
     * @return
     */
    public static boolean isAlleleInList(final Allele ref1, final Allele alt1, final Allele ref2, final List<Allele> altList) {
        final Allele commonRef;
        if (ref1.equals(ref2)) {
            return altList.contains(alt1);
        } else {
                commonRef = determineReferenceAllele(ref1, ref2);
        }
        final Map<Allele, Allele> alleleMap;
        if (ref1.equals(commonRef)) {
            alleleMap = GATKVariantContextUtils.createAlleleMapping(commonRef, ref2, altList);
            return alleleMap.values().contains(alt1);
        } else if (ref2.equals(commonRef)) {
            alleleMap = GATKVariantContextUtils.createAlleleMapping(commonRef, ref1, Arrays.asList(alt1));
            return altList.contains(alleleMap.get(alt1));
        } else {
            throw new IllegalStateException("Reference alleles " + ref1 + " and " + ref2 + " have common reference allele " + commonRef + " which is equal to neither.");
        }
    }

    /**
     *  Find the indices in one allele list (or OptionalInt.empty() if none exists) in another list of alleles
     *  We assume that reference alleles are the first element of each list and don't assume that either list is in its minimal representation
     * @param alleles1  The alleles whose indices we want to find within {@code alleles2}
     * @param alleles2  The alleles where we find the occurrence of alleles1
     * @return  Example: alleles1 = {A, G, T}; alleles2 = {AA, GA, CC} output = {0,1,EMPTY}
     */
    public static List<OptionalInt> alleleIndices(final List<Allele> alleles1, final List<Allele> alleles2) {
        Utils.validateArg(!alleles1.isEmpty() && !alleles2.isEmpty(), "alleles lists must at least contain a reference allele.");
        final Allele ref1 = alleles1.get(0);
        final Allele ref2 = alleles2.get(0);
        Utils.validateArg(ref1.isReference() && ref2.isReference(), "First allele in each list must be reference.");

        final Allele commonRef = determineReferenceAllele(ref1, ref2);
        final Map<Allele, Allele> alleles1ToCommon = createAlleleMapping(commonRef, ref1, alleles1.subList(1, alleles1.size()));
        final Map<Allele, Allele> alleles2ToCommon = createAlleleMapping(commonRef, ref2, alleles2.subList(1, alleles2.size()));

        final Map<Allele, Integer> commonToIndex2 = IntStream.range(0, alleles2.size()).boxed()
                .collect(Collectors.toMap(n -> alleles2ToCommon.get(alleles2.get(n)), n -> n));

        return alleles1.stream()
                .map(a -> commonToIndex2.getOrDefault(alleles1ToCommon.get(a), -1))
                .map(n -> n < 0 ? OptionalInt.empty() : OptionalInt.of(n))
                .collect(Collectors.toList());
    }

    /**
     * Determines the common reference allele
     *
     * @param VCs    the list of VariantContexts
     * @param loc    if not null, ignore records that do not begin at this start location
     * @return possibly null Allele
     */
    public static Allele determineReferenceAllele(final List<VariantContext> VCs, final Locatable loc) {
        Allele ref = null;

        for ( final VariantContext vc : VCs ) {
            if ( contextMatchesLoc(vc, loc) ) {
                final Allele myRef = vc.getReference();
                try {
                    ref = determineReferenceAllele(ref, myRef);
                } catch (IllegalStateException e) {
                    throw new IllegalStateException(String.format("The provided variant file(s) have inconsistent references " +
                            "for the same position(s) at %s:%d, %s vs. %s", vc.getContig(), vc.getStart(), ref, myRef));
                }
            }
        }
        return ref;
    }

    public static Allele determineReferenceAllele(final Allele ref1, final Allele ref2) {
        if ( ref1 == null || ref1.length() < ref2.length() ) {
            return ref2;
        } else if ( ref2 == null || ref2.length() < ref1.length()) {
            return ref1;
        }
        else if ( ref1.length() == ref2.length() && ! ref1.equals(ref2) ) {
            throw new IllegalStateException(String.format("The provided reference alleles do not appear to represent the same position, %s vs. %s", ref1, ref2));
        } else {  //the lengths are the same and they're equal, so we could return ref1 or ref2
            return ref1;
        }
    }

    /**
     * Calculates the total ploidy of a variant context as the sum of all plodies across genotypes.
     * @param vc the target variant context.
     * @param defaultPloidy the default ploidy to be assume when there is no ploidy information for a genotype.
     * @return never {@code null}.
     */
    public static int totalPloidy(final VariantContext vc, final int defaultPloidy) {
        Utils.nonNull(vc, "the vc provided cannot be null");
        Utils.validateArg(defaultPloidy >= 0, "the default ploidy must 0 or greater");
        return vc.getGenotypes().stream().mapToInt(Genotype::getPloidy)
                .map(p -> p <= 0 ? defaultPloidy : p).sum();
    }

    /**
     * Returns the type of the substitution (transition vs transversion) represented by this variant. Only applicable to bi allelic SNPs.
     */
    public static BaseUtils.BaseSubstitutionType getSNPSubstitutionType(final VariantContext context) {
        Utils.nonNull(context);
        if (!context.isSNP() || !context.isBiallelic()) {
            throw new IllegalArgumentException("Requested SNP substitution type for bialleic non-SNP " + context);
        }
        return BaseUtils.SNPSubstitutionType(context.getReference().getBases()[0], context.getAlternateAllele(0).getBases()[0]);
    }

    /**
     * If this is a BiAllelic SNP, is it a transition?
     */
    public static boolean isTransition(final VariantContext context) {
        Utils.nonNull(context);
        return getSNPSubstitutionType(context) == BaseUtils.BaseSubstitutionType.TRANSITION;
    }


    /**
     * Returns a homozygous call allele list given the only allele and the ploidy.
     *
     * @param allele the only allele in the allele list.
     * @param ploidy the ploidy of the resulting allele list.
     *
     * @throws IllegalArgumentException if {@code allele} is {@code null} or ploidy is negative.
     *
     * @return never {@code null}.
     */
    public static List<Allele> homozygousAlleleList(final Allele allele, final int ploidy) {
        Utils.nonNull(allele);
        Utils.validateArg(ploidy >= 0, "ploidy cannot be negative");

        // Use a tailored inner class to implement the list:
        return Collections.nCopies(ploidy,allele);
    }

    /**
     * Add the genotype call (GT) field to GenotypeBuilder using the requested {@link GenotypeAssignmentMethod}
     *
     * @param ploidy output ploidy for no-call GTs and likelihood array length
     * @param gb the builder where we should put our newly called alleles, cannot be null
     * @param assignmentMethod the method to use to do the assignment, cannot be null
     * @param genotypeLikelihoods a vector of likelihoods to use if the method requires PLs, should be log10 likelihoods, cannot be null
     * @param allelesToUse the alleles with respect to which the likelihoods are defined
     * @param originalGT Genotype that includes GQ when available
     * @param gpc utility class to help with likelihood calculations
     */
    public static void makeGenotypeCall(final int ploidy,
                                        final GenotypeBuilder gb,
                                        final GenotypeAssignmentMethod assignmentMethod,
                                        final double[] genotypeLikelihoods,
                                        final List<Allele> allelesToUse,
                                        final Genotype originalGT,
                                        final GenotypePriorCalculator gpc) {
        if(originalGT == null && assignmentMethod == GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL) {
            throw new IllegalArgumentException("original GT cannot be null if assignmentMethod is BEST_MATCH_TO_ORIGINAL");
        }
        //only do this if we're not no-calling so CombineGVCFs and GenomicsDB keep attributes
        if (assignmentMethod != GenotypeAssignmentMethod.SET_TO_NO_CALL && (originalGT.isHomRef() || originalGT.isNoCall()) && originalGT.hasGQ() && originalGT.getGQ() == 0) {
            gb.alleles(noCallAlleles(ploidy));
            //in the "no data" case set to plain no-call (./.)
            if (originalGT.hasDP() && originalGT.getDP() == 0) {
                gb.noPL().noDP().noAD().noGQ().noAttributes();
                return;
            }
        }
        if (assignmentMethod == GenotypeAssignmentMethod.SET_TO_NO_CALL) {
            gb.alleles(noCallAlleles(ploidy));
        } else if (assignmentMethod == GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN ||
                    assignmentMethod == GenotypeAssignmentMethod.PREFER_PLS) {
            if ((genotypeLikelihoods == null || !isInformative(genotypeLikelihoods)) && assignmentMethod == GenotypeAssignmentMethod.PREFER_PLS) {
                if (originalGT == null) {
                    throw new IllegalArgumentException("original GT cannot be null if assignmentMethod is PREFER_PLS");
                } else {
                    gb.alleles(bestMatchToOriginalGT(allelesToUse, originalGT.getAlleles()));
                }
            } else {
                if (genotypeLikelihoods == null) {
                    gb.alleles(bestMatchToOriginalGT(allelesToUse, originalGT.getAlleles()));
                } else {
                    final int maxLikelihoodIndex = MathUtils.maxElementIndex(genotypeLikelihoods);
                    final GenotypeAlleleCounts alleleCounts = GenotypesCache.get(ploidy, maxLikelihoodIndex);
                    final List<Allele> finalAlleles = alleleCounts.asAlleleList(allelesToUse);
                    final double gq = GenotypeLikelihoods.getGQLog10FromLikelihoods(maxLikelihoodIndex, genotypeLikelihoods);
                    if (finalAlleles.contains(Allele.NON_REF_ALLELE)) {
                        final Allele ref = allelesToUse.stream().filter(Allele::isReference).collect(Collectors.toList()).get(0);
                        gb.alleles(Collections.nCopies(ploidy, ref));
                        gb.PL(new int[genotypeLikelihoods.length]).log10PError(0);
                    } else if (maxLikelihoodIndex == 0 && gq > GATKVariantContextUtils.SUM_GL_THRESH_NOCALL) { //GLs are log10s, so they're negative, GQ is also negative, so less negative is greater than
                        gb.alleles(Collections.nCopies(ploidy, Allele.NO_CALL));
                    } else {
                        gb.alleles(finalAlleles);
                    }
                    final int numAltAlleles = allelesToUse.size() - 1;
                    if (numAltAlleles > 0) {
                        gb.log10PError(gq);
                    }
                }
            }
        } else if (assignmentMethod == GenotypeAssignmentMethod.SET_TO_NO_CALL_NO_ANNOTATIONS) {
            gb.alleles(noCallAlleles(ploidy)).noGQ().noAD().noPL().noAttributes();
        } else if (assignmentMethod == GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL) {
            if (originalGT.hasGQ() && originalGT.getGQ() == 0 &&
                    (!originalGT.hasPL() || originalGT.getPL()[0] == 0)) {
                gb.alleles(noCallAlleles(ploidy));
            } else {
                gb.alleles(bestMatchToOriginalGT(allelesToUse, originalGT.getAlleles()));
            }
        } else if (assignmentMethod == GenotypeAssignmentMethod.USE_POSTERIOR_PROBABILITIES) {
            if (gpc == null) {
                throw new GATKException("cannot uses posteriors without an genotype prior calculator present");
            } else {
                // Calculate posteriors.
                final double[] log10Priors = gpc.getLog10Priors(ploidy, allelesToUse);
                final double[] log10Posteriors = MathUtils.ebeAdd(log10Priors, genotypeLikelihoods);
                final double[] normalizedLog10Posteriors = MathUtils.scaleLogSpaceArrayForNumericalStability(log10Posteriors);
                // Update GP and PG annotations:
                gb.attribute(VCFConstants.GENOTYPE_POSTERIORS_KEY, Arrays.stream(normalizedLog10Posteriors)
                        .map(v -> v == 0.0 ? 0.0 : v * -10) // the reason for the == 0.0 is to avoid a signed 0 output "-0.0"
                        .mapToObj(GATKVariantContextUtils::formatGP).toArray());
                gb.attribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY, Arrays.stream(log10Priors)
                        .map(v -> v == 0.0 ? 0.0 : v * -10)
                        .mapToObj(GATKVariantContextUtils::formatGP).toArray());
                // Set the GQ accordingly
                final int maxPosteriorIndex = MathUtils.maxElementIndex(log10Posteriors);
                if ( allelesToUse.size() > 0 ) {
                    gb.log10PError(getGQLog10FromPosteriors(maxPosteriorIndex, normalizedLog10Posteriors));
                }
                // Finally we update the genotype alleles.
                gb.alleles(GenotypesCache.get(ploidy, maxPosteriorIndex).asAlleleList(allelesToUse));
            }
        }
    }

    private static List<Allele> bestMatchToOriginalGT(final List<Allele> allelesToUse, final List<Allele> originalGT) {
        final List<Allele> best = new LinkedList<>();
        final Allele ref = allelesToUse.get(0);
        for (final Allele originalAllele : originalGT) {
            best.add((allelesToUse.contains(originalAllele) || originalAllele.isNoCall()) ? originalAllele : ref);
        }
        return best;
    }

    @VisibleForTesting
    static double getGQLog10FromPosteriors(final int bestGenotypeIndex, final double[] /**/log10Posteriors) {
        if (bestGenotypeIndex < 0) {
            return CommonInfo.NO_LOG10_PERROR;
        } else {
            switch (log10Posteriors.length) {
                case 0:
                case 1: return CommonInfo.NO_LOG10_PERROR;
                case 2: return bestGenotypeIndex == 0 ? log10Posteriors[1] : log10Posteriors[0];
                case 3: return Math.min(0.0, MathUtils.log10SumLog10(
                                 log10Posteriors[ bestGenotypeIndex == 0 ? 2 : bestGenotypeIndex - 1],
                                 log10Posteriors[ bestGenotypeIndex == 2 ? 0 : bestGenotypeIndex + 1]));
                default:
                    if (bestGenotypeIndex == 0) {
                        return Math.min(0.0, MathUtils.log10SumLog10(log10Posteriors, 1, log10Posteriors.length));
                    } else if (bestGenotypeIndex == log10Posteriors.length - 1) {
                        return Math.min(0.0, MathUtils.log10SumLog10(log10Posteriors, 0, bestGenotypeIndex));
                    } else {
                        return Math.min(0.0, MathUtils.log10SumLog10(
                                MathUtils.log10sumLog10(log10Posteriors, 0, bestGenotypeIndex),
                                MathUtils.log10sumLog10(log10Posteriors, bestGenotypeIndex + 1, log10Posteriors.length)
                        ));
                    }
            }
        }
    }

    private static String formatGP(final double gp) {
        final String formatted = String.format("%.2f", gp);
        int last = formatted.length() - 1;
        if (formatted.charAt(last) == '0') {
            if (formatted.charAt(--last) == '0') {
                return formatted.substring(0, --last); // exclude the '.' as-well.
            } else {
                return formatted.substring(0, ++last);
            }
        } else {
            return formatted;
        }
    }

    /**
     *
     * @param ploidy
     * @param allele
     * @return (return a single no-call allele for ploidy 0, as on Y for females)
     */
    public static List<Allele> makePloidyLengthAlleleList(final int ploidy, final Allele allele) {
        if (ploidy == 0) {
            return Arrays.asList(Allele.NO_CALL);
        }
        final List<Allele> repeatedList = new ArrayList<>();
        for (int i = 0; i < ploidy; i++) {
            repeatedList.add(allele);
        }
        return repeatedList;
    }

    /**
     * Return the rightmost variant context in maybeOverlapping that overlaps curPos
     *
     * @param curPos non-null genome loc
     * @param maybeOverlapping a collection of variant contexts that might overlap curPos
     * @return the rightmost VariantContext, or null if none overlaps
     */
    public static VariantContext getOverlappingVariantContext(final Locatable curPos, final Collection<VariantContext> maybeOverlapping) {
        VariantContext overlaps = null;
        for ( final VariantContext vc : maybeOverlapping ) {
            if ( curPos.overlaps(vc) ) {
                if ( overlaps == null || vc.getStart() > overlaps.getStart() ) {
                    overlaps = vc;
                }
            }
        }
        return overlaps;
    }

    /**
     * Determines whether the provided VariantContext has real alternate alleles.
     *
     * @param vc  the VariantContext to evaluate
     * @return true if it has proper alternate alleles, false otherwise
     */
    public static boolean isProperlyPolymorphic(final VariantContext vc) {
        //obvious cases
        if (vc == null || vc.getAlternateAlleles().isEmpty()) {
            return false;
        } else if (vc.isBiallelic()) {
            return !(GATKVCFConstants.isSpanningDeletion(vc.getAlternateAllele(0)) || vc.isSymbolic());
        } else if (GATKVCFConstants.isSpanningDeletion(vc.getAlternateAllele(0)) && vc.getAlternateAllele(1).equals(Allele.NON_REF_ALLELE)){
            return false;
        } else {
            return true;
        }
    }

    /** This is lifted directly from htsjdk with some minor modifications!  However, it is a private method there.
     *
     * This method cannot return {@link VariantContext.Type} MIXED
     *
     * Please see https://github.com/samtools/htsjdk/issues/999
     *
     * <p>Here are some cases that will not work properly, though this may not be an issue in practice:  </p>
     *  <ul>
     *      <li>"CGT" --> "GGA" this will be a MNP, but really it is two SNPs.</li>
     *      <li>Spanning deletions for alternate will show as {@link VariantContext.Type} NO_VARIATION</li>
     *      <li>Spanning deletions for reference will throw exception. </li>
     *      <li>Reference that is symbolic will throw an exception.</li>
     *  </ul>
     *
     * @param ref reference allele. Never {@code null}
     * @param allele alternate allele to compare. Never {@code null}
     * @return
     */
    public static VariantContext.Type typeOfVariant(final Allele ref, final Allele allele) {
        Utils.nonNull(ref);
        Utils.nonNull(allele);

        if ( ref.isSymbolic() )
            throw new IllegalStateException("Unexpected error: encountered a record with a symbolic reference allele");

        if ( allele.isSymbolic() )
            return VariantContext.Type.SYMBOLIC;

        if (allele.equals(Allele.SPAN_DEL)) {
            return VariantContext.Type.NO_VARIATION;
        }

        if ( ref.equals(Allele.SPAN_DEL) )
            throw new IllegalStateException("Unexpected error: encountered a record with a spanning deletion reference allele");

        if ( ref.length() == allele.length() ) {
            if (ref.basesMatch(allele)) {
                return VariantContext.Type.NO_VARIATION;
            } else if ( allele.length() == 1 )
                return VariantContext.Type.SNP;

            // If the two alleles are the same length and only differ by one base, then still a SNP.
            else if (IntStream.range(0, ref.length()).filter(i -> ref.getBases()[i] != allele.getBases()[i]).count() == 1) {
                return VariantContext.Type.SNP;
            } else
                return VariantContext.Type.MNP;
        }

        // Important note: previously we were checking that one allele is the prefix of the other.  However, that's not an
        // appropriate check as can be seen from the following example:
        // REF = CTTA and ALT = C,CT,CA
        // This should be assigned the INDEL type but was being marked as a MIXED type because of the prefix check.
        // In truth, it should be absolutely impossible to return a MIXED type from this method because it simply
        // performs a pairwise comparison of a single alternate allele against the reference allele (whereas the MIXED type
        // is reserved for cases of multiple alternate alleles of different types).  Therefore, if we've reached this point
        // in the code (so we're not a SNP, MNP, or symbolic allele), we absolutely must be an INDEL.

        return VariantContext.Type.INDEL;

        // old incorrect logic:
        // if (oneIsPrefixOfOther(ref, allele))
        //     return Type.INDEL;
        // else
        //     return Type.MIXED;
    }

    /**
     *  This method should only be run on variants that are known to be indels.  See {@code typeOfVariant}
     *
     *<p>Here are some cases that will not work properly, though this may not be an issue in practice:  </p>
     *  <ul>
     *      <li>"CT" --> "CATT" this is really just a simple AT insertion, but this will show up as complex.</li>
     *  </ul>
     * @param ref reference allele. Never {@code null}
     * @param allele alternate allele to compare. Never {@code null}
     * @return true if the indel is complex (for example, also includes a SNP), false if simple indel.  If the input alleles define a variant that is not
     *  an indel, then the behavior of this method is undefined (though will probably just return false).
     *
     */
    public static boolean isComplexIndel(final Allele ref, final Allele allele) {

        Utils.nonNull(ref);
        Utils.nonNull(allele);

        // Symbolic --> false
        if (ref.isSymbolic() || (ref.length() == 0)) {
            return false;
        }
        if (allele.isSymbolic() || (allele.length() == 0)) {
            return false;
        }

        // SNP, MNP, or no variation --> false
        if ( ref.length() == allele.length() ) {
            return false;
        }

        // obvious simple del or simple indel
        if ((allele.length() == 1) || (ref.length() == 1)) {
            return false;
        }

        // If the ref starts with the alt or vice versa, this is still simple.
        if (allele.length() > ref.length()) {
            final boolean isAltStartsWithRef = IntStream.range(0, ref.length()).allMatch(i -> ref.getBases()[i] == allele.getBases()[i]);
            return !isAltStartsWithRef;
        } else {
            final boolean isRefStartsWithAlt = IntStream.range(0, allele.length()).allMatch(i -> ref.getBases()[i] == allele.getBases()[i]);
            return !isRefStartsWithAlt;
        }
    }

    /**
     * Given a set of alleles (reference and alternate), choose the allele that is the best match for the given read (and offset)
     * TODO: This method cannot recognize equivalent alleles (See https://github.com/broadinstitute/gatk/issues/5061)
     * @param pileupElement read and offset.  Never {@code null}
     * @param referenceAllele Reference allele.  Never {@code null}
     * @param altAlleles List of candidate alternate alleles.  Never {@code null}
     * @param minBaseQualityCutoff minimum base quality for the bases that match the allele in order to be counted.
     *                             Must be positive or zero.  If you do not want any filtering, specify 0.
     * @return The allele (reference or from the altAlleles) that matches.  {@code null} if none are a match or the base qualities
     *      corresponding to the allele don't all exceed the minimum.
     */
    public static Allele chooseAlleleForRead(final PileupElement pileupElement, final Allele referenceAllele, final List<Allele> altAlleles, int minBaseQualityCutoff) {
        Utils.nonNull(pileupElement);
        Utils.nonNull(referenceAllele);
        Utils.nonNull(altAlleles);
        ParamUtils.isPositiveOrZero(minBaseQualityCutoff, "Minimum base quality must be positive or zero.");

        final boolean isRef = referenceAllele.basesMatch(getBasesForAlleleInRead(pileupElement, referenceAllele))
                && !pileupElement.isBeforeDeletionStart() && !pileupElement.isBeforeInsertion();

        Allele pileupAllele = null;
        if (!isRef) {

            for (Allele altAllele : altAlleles) {
                final VariantContext.Type variantType = typeOfVariant(referenceAllele, altAllele);

                if (variantType == VariantContext.Type.INDEL) {
                    if (isIndelInThePileupElement(pileupElement, referenceAllele, altAllele)) {
                        pileupAllele = altAllele;
                    }

                } else if (variantType == VariantContext.Type.MNP || variantType == VariantContext.Type.SNP) {
                    if (doesReadContainAllele(pileupElement, altAllele) == Trilean.TRUE) {
                        pileupAllele = altAllele;
                    }
                }
            }
        } else {
            pileupAllele = referenceAllele;
        }

        if ((pileupAllele != null) && (getMinBaseQualityForAlleleInRead(pileupElement, pileupAllele) < minBaseQualityCutoff)) {
            pileupAllele = null;
        }

        return pileupAllele;
    }

    private static boolean isIndelInThePileupElement(final PileupElement pileupElement, final Allele referenceAllele, final Allele altAllele) {
        boolean isAltAlleleInThePileup = false;

        // Check insertion
        if (pileupElement.isBeforeInsertion()) {
            final String insertionBases = pileupElement.getBasesOfImmediatelyFollowingInsertion();
            // edge case: ignore a deletion immediately preceding an insertion as p.getBasesOfImmediatelyFollowingInsertion() returns null [EB]
            if ((insertionBases != null) && (Allele.extend(referenceAllele, insertionBases.getBytes()).basesMatch(altAllele))) {
                isAltAlleleInThePileup = true;
            }
        } else if (pileupElement.isBeforeDeletionStart()) {
            final int deletionLength = pileupElement.getLengthOfImmediatelyFollowingIndel();
            if ((referenceAllele.getBases().length - altAllele.getBases().length) == deletionLength) {
                isAltAlleleInThePileup = true;
            }
        }
        return isAltAlleleInThePileup;
    }

    /**
     * @param pileupElement pileup element representing the read.  Never {@code null}
     * @param allele allele to get overlapping bases in read.  Never {@code null}
     * @return array of the bytes that correspond to the allele in the pileup element.  Note that, if the read ends, this
     * list can be smaller than the length of the allele.
     */
    private static byte[] getBasesForAlleleInRead(final PileupElement pileupElement, final Allele allele) {
        Utils.nonNull(pileupElement);
        Utils.nonNull(allele);
        return ArrayUtils.subarray(pileupElement.getRead().getBases(), pileupElement.getOffset(), pileupElement.getOffset() + allele.getBases().length);
    }

    /**
     * TODO: Test.  And make sure to test with reference alleles to make sure "*" is not included.
     * @param pileupElement pileup element representing the read.  Never {@code null}
     * @param allele query allele.  Never {@code null}
     * @return Whether the read contains the allele.  Note that unknown can occur as well.
     */
    public static Trilean doesReadContainAllele(final PileupElement pileupElement, final Allele allele) {
        Utils.nonNull(pileupElement);
        Utils.nonNull(allele);

        final byte[] readBases = ArrayUtils.subarray(pileupElement.getRead().getBases(), pileupElement.getOffset(), pileupElement.getOffset() + allele.getBases().length);

        if (readBases.length < allele.getBases().length) {
            return Trilean.UNKNOWN;
        }

        if (allele.basesMatch(readBases)) {
            return Trilean.TRUE;
        } else {
            return Trilean.FALSE;
        }
    }

    /**
     * Find the minimum base quality for all bases in a read that correspond to a given allele.
     *
     * @param pileupElement pileup element representing the read.  Never {@code null}
     * @param allele query allele.  Never {@code null}
     * @return lowest base quality seen in the corresponding bases of the read
     */
    private static int getMinBaseQualityForAlleleInRead(final PileupElement pileupElement, final Allele allele) {
        Utils.nonNull(pileupElement);
        Utils.nonNull(allele);
        final byte[] alleleBases = allele.getBases();
        final byte[] pileupBaseQualities = ArrayUtils.subarray(pileupElement.getRead().getBaseQualities(), pileupElement.getOffset(), pileupElement.getOffset() + alleleBases.length);
        final OptionalInt minQuality = IntStream.range(0, pileupBaseQualities.length).map(i -> Byte.toUnsignedInt(pileupBaseQualities[i])).min();
        return minQuality.orElse(-1);
    }

    public static boolean containsInlineIndel(final VariantContext vc) {
        return containsInlineIndel(vc.getReference(), vc.getAlternateAlleles());
    }

    public static boolean containsInlineIndel(final Allele ref, final Collection<Allele> altAlleles) {
        final int refLength = ref.length();
        for (final Allele alt : altAlleles) {
            if (!alt.isSymbolic() && alt != Allele.SPAN_DEL && alt.length() != refLength) {
                return true;
            }
        }
        return false;
    }

    public enum GenotypeMergeType {
        /**
         * Make all sample genotypes unique by file. Each sample shared across RODs gets named sample.ROD.
         */
        UNIQUIFY,
        /**
         * Take genotypes in priority order (see the priority argument).
         */
        PRIORITIZE,
        /**
         * Take the genotypes in any order.
         */
        UNSORTED,
        /**
         * Require that all samples/genotypes be unique between all inputs.
         */
        REQUIRE_UNIQUE
    }

    public enum FilteredRecordMergeType {
        /**
         * Union - leaves the record if any record is unfiltered.
         */
        KEEP_IF_ANY_UNFILTERED,
        /**
         * Requires all records present at site to be unfiltered. VCF files that don't contain the record don't influence this.
         */
        KEEP_IF_ALL_UNFILTERED,
        /**
         * If any record is present at this site (regardless of possibly being filtered), then all such records are kept and the filters are reset.
         */
        KEEP_UNCONDITIONAL
    }

    /**
     * Returns true iff VC is an non-complex indel where every allele represents an expansion or
     * contraction of a series of identical bases in the reference.
     *
     * For example, suppose the ref bases are CTCTCTGA, which includes a 3x repeat of CTCTCT
     *
     * If VC = -/CT, then this function returns true because the CT insertion matches exactly the
     * upcoming reference.
     * If VC = -/CTA then this function returns false because the CTA isn't a perfect match
     *
     * Now consider deletions:
     *
     * If VC = CT/- then again the same logic applies and this returns true
     * The case of CTA/- makes no sense because it doesn't actually match the reference bases.
     *
     * The logic of this function is pretty simple.  Take all of the non-null alleles in VC.  For
     * each insertion allele of n bases, check if that allele matches the next n reference bases.
     * For each deletion allele of n bases, check if this matches the reference bases at n - 2 n,
     * as it must necessarily match the first n bases.  If this test returns true for all
     * alleles you are a tandem repeat, otherwise you are not.
     *
     * @param vc
     * @param refBasesStartingAtVCWithPad not this is assumed to include the PADDED reference
     * @return
     */
    public static boolean isTandemRepeat(final VariantContext vc, final byte[] refBasesStartingAtVCWithPad) {
        final String refBasesStartingAtVCWithoutPad = new String(refBasesStartingAtVCWithPad).substring(1);
        if ( ! vc.isIndel() ) // only indels are tandem repeats
            return false;

        final Allele ref = vc.getReference();

        for ( final Allele allele : vc.getAlternateAlleles() ) {
            if ( ! isRepeatAllele(ref, allele, refBasesStartingAtVCWithoutPad) )
                return false;
        }

        // we've passed all of the tests, so we are a repeat
        return true;
    }

    /**
     *
     * @param refBasesStartingAtVCWithoutPad    Ref bases excluding the initial base of the variant context where the alt matches the ref.
     *                                          For example, if the reference sequence is GATCCACCACCAGTCGA and we have a deletion
     *                                          of one STR unit CCA, it is represented as a variant context TCCA -> T, where the 'T' is
     *                                          the padding base.  In this case, {@code refBasesStartingAtVCWithoutPad} is CCACCACCAGTCGA.
     * @return
     */
    public static Pair<List<Integer>, byte[]> getNumTandemRepeatUnits(final Allele refAllele, final List<Allele> altAlleles, final byte[] refBasesStartingAtVCWithoutPad) {
        Utils.nonNull(refAllele);
        Utils.nonNull(altAlleles);
        Utils.nonNull(refBasesStartingAtVCWithoutPad);

        if ( altAlleles.stream().allMatch(a -> a.length() == refAllele.length()) ){ // only indels are tandem repeats
            return null;
        }

        final byte[] refAlleleBases = Arrays.copyOfRange(refAllele.getBases(), 1, refAllele.length());

        byte[] repeatUnit = null;
        final List<Integer> lengths = new ArrayList<>();

        for ( final Allele allele : altAlleles ) {
            Pair<int[],byte[]> result = getNumTandemRepeatUnits(refAlleleBases, Arrays.copyOfRange(allele.getBases(), 1, allele.length()), refBasesStartingAtVCWithoutPad);

            final int[] repetitionCount = result.getLeft();
            // repetition count = 0 means allele is not a tandem expansion of context
            if (repetitionCount[0] == 0 || repetitionCount[1] == 0)
                return null;

            if (lengths.isEmpty()) {
                lengths.add(repetitionCount[0]); // add ref allele length only once
            }
            lengths.add(repetitionCount[1]);  // add this alt allele's length

            repeatUnit = result.getRight();
        }

        return new MutablePair<>(lengths,repeatUnit);
    }

    public static Pair<List<Integer>, byte[]> getNumTandemRepeatUnits(final VariantContext vc, final byte[] refBasesStartingAtVCWithoutPad) {
        return getNumTandemRepeatUnits(vc.getReference(), vc.getAlternateAlleles(), refBasesStartingAtVCWithoutPad);
    }

    public static Pair<int[],byte[]> getNumTandemRepeatUnits(final byte[] refBases, final byte[] altBases, final byte[] remainingRefContext) {
         /* we can't exactly apply same logic as in basesAreRepeated() to compute tandem unit and number of repeated units.
           Consider case where ref =ATATAT and we have an insertion of ATAT. Natural description is (AT)3 -> (AT)2.
         */

        byte[] longB;
        // find first repeat unit based on either ref or alt, whichever is longer
        if (altBases.length > refBases.length)
            longB = altBases;
        else
            longB = refBases;

        // see if non-null allele (either ref or alt, whichever is longer) can be decomposed into several identical tandem units
        // for example, -*,CACA needs to first be decomposed into (CA)2
        final int repeatUnitLength = findRepeatedSubstring(longB);
        final byte[] repeatUnit = Arrays.copyOf(longB, repeatUnitLength);

        final int[] repetitionCount = new int[2];
        // look for repetitions forward on the ref bases (i.e. starting at beginning of ref bases)
        int repetitionsInRef = findNumberOfRepetitions(repeatUnit, refBases, true);
        repetitionCount[0] = findNumberOfRepetitions(repeatUnit, ArrayUtils.addAll(refBases, remainingRefContext), true)-repetitionsInRef;
        repetitionCount[1] = findNumberOfRepetitions(repeatUnit, ArrayUtils.addAll(altBases, remainingRefContext), true)-repetitionsInRef;

        return new MutablePair<>(repetitionCount, repeatUnit);

    }

    /**
     * Find out if a string can be represented as a tandem number of substrings.
     * For example ACTACT is a 2-tandem of ACT,
     * but ACTACA is not.
     *
     * @param bases                 String to be tested
     * @return                      Length of repeat unit, if string can be represented as tandem of substring (if it can't
     *                              be represented as one, it will be just the length of the input string)
     */
    public static int findRepeatedSubstring(byte[] bases) {

        int repLength;
        for (repLength=1; repLength <=bases.length; repLength++) {
            final byte[] candidateRepeatUnit = Arrays.copyOf(bases,repLength);
            boolean allBasesMatch = true;
            for (int start = repLength; start < bases.length; start += repLength ) {
                // check that remaining of string is exactly equal to repeat unit
                final byte[] basePiece = Arrays.copyOfRange(bases,start,start+candidateRepeatUnit.length);
                if (!Arrays.equals(candidateRepeatUnit, basePiece)) {
                    allBasesMatch = false;
                    break;
                }
            }
            if (allBasesMatch)
                return repLength;
        }

        return repLength;
    }

    /**
     * Finds number of repetitions a string consists of.
     * For example, for string ATAT and repeat unit AT, number of repetitions = 2
     * @param repeatUnit             Non-empty substring represented by byte array
     * @param testString             String to test (represented by byte array), may be empty
     * @param leadingRepeats         Look for leading (at the beginning of string) or trailing (at end of string) repetitions
     * For example:
     *    GATAT has 0 leading repeats of AT but 2 trailing repeats of AT
     *    ATATG has 1 leading repeat of A but 2 leading repeats of AT
     *    CCCCCCCC has 2 leading and 2 trailing repeats of CCC
     *
     * @return  Number of repetitions (0 if testString is not a concatenation of n repeatUnit's, including the case of empty testString)
     */
    public static int findNumberOfRepetitions(byte[] repeatUnit, byte[] testString, boolean leadingRepeats) {
        Utils.nonNull(repeatUnit, "repeatUnit");
        Utils.nonNull(testString, "testString");
        Utils.validateArg(repeatUnit.length != 0, "empty repeatUnit");
        if (testString.length == 0){
            return 0;
        }
        return findNumberOfRepetitions(repeatUnit, 0, repeatUnit.length, testString, 0, testString.length, leadingRepeats);
    }

    /**
     * Finds number of repetitions a string consists of.
     * Same as {@link #findNumberOfRepetitions} but operates on subarrays of a bigger array to save on copying.
     * For example, for string ATAT and repeat unit AT, number of repetitions = 2
     * @param repeatUnitFull             Non-empty substring represented by byte array
     * @param offsetInRepeatUnitFull     the offset in repeatUnitFull from which to read the repeat unit
     * @param repeatUnitLength           length of the repeat unit
     * @param testStringFull             string to test (represented by byte array), may be empty
     * @param offsetInTestStringFull     the offset in offsetInRepeatUnitFull from which to read the test string
     * @param testStringLength           length of the test string
     * @param leadingRepeats         Look for leading (at the beginning of string) or trailing (at end of string) repetitions
     * For example:
     *    GATAT has 0 leading repeats of AT but 2 trailing repeats of AT
     *    ATATG has 1 leading repeat of A but 2 leading repeats of AT
     *    CCCCCCCC has 2 leading and 2 trailing repeats of CCC
     * @return  Number of repetitions (0 if testString is not a concatenation of n repeatUnit's, including the case of empty testString)
     */
    public static int findNumberOfRepetitions(final byte[] repeatUnitFull, final int offsetInRepeatUnitFull, final int repeatUnitLength, final byte[] testStringFull, final int offsetInTestStringFull, final int testStringLength, final boolean leadingRepeats) {
        Utils.nonNull(repeatUnitFull, "repeatUnit");
        Utils.nonNull(testStringFull, "testString");
        Utils.validIndex(offsetInRepeatUnitFull, repeatUnitFull.length);
        Utils.validateArg(repeatUnitLength >= 0 && repeatUnitLength <= repeatUnitFull.length, "repeatUnitLength");
        if (testStringLength == 0){
            return 0;
        }
        Utils.validIndex(offsetInTestStringFull, testStringFull.length);
        Utils.validateArg(testStringLength >= 0 && testStringLength <= testStringFull.length, "testStringLength");
        final int lengthDifference = testStringLength - repeatUnitLength;

        if (leadingRepeats) {
            int numRepeats = 0;
            // look forward on the test string
            for (int start = 0; start <= lengthDifference; start += repeatUnitLength) {
                if(Utils.equalRange(testStringFull, start + offsetInTestStringFull, repeatUnitFull, offsetInRepeatUnitFull, repeatUnitLength)) {
                    numRepeats++;
                } else {
                    return numRepeats;
                }
            }
            return numRepeats;
        } else {
            // look backward. For example, if repeatUnit = AT and testString = GATAT, number of repeat units is still 2
            int numRepeats = 0;
            // look backward on the test string
            for (int start = lengthDifference; start >= 0; start -= repeatUnitLength) {
                if (Utils.equalRange(testStringFull, start + offsetInTestStringFull, repeatUnitFull, offsetInRepeatUnitFull, repeatUnitLength)) {
                    numRepeats++;
                } else {
                    return numRepeats;
                }
            }
            return numRepeats;
        }
    }

    /**
     * Helper function for isTandemRepeat that checks that allele matches somewhere on the reference
     */
    protected static boolean isRepeatAllele(final Allele ref, final Allele alt, final String refBasesStartingAtVCWithoutPad) {
        if ( ! Allele.oneIsPrefixOfOther(ref, alt) )
            return false; // we require one allele be a prefix of another

        if ( ref.length() > alt.length() ) { // we are a deletion
            return basesAreRepeated(ref.getBaseString(), alt.getBaseString(), refBasesStartingAtVCWithoutPad, 2);
        } else { // we are an insertion
            return basesAreRepeated(alt.getBaseString(), ref.getBaseString(), refBasesStartingAtVCWithoutPad, 1);
        }
    }

    protected static boolean basesAreRepeated(final String l, final String s, final String ref, final int minNumberOfMatches) {
        final String potentialRepeat = l.substring(s.length()); // skip s bases

        for ( int i = 0; i < minNumberOfMatches; i++) {
            final int start = i * potentialRepeat.length();
            final int end = (i+1) * potentialRepeat.length();
            if ( ref.length() < end )
                return false; // we ran out of bases to test
            final String refSub = ref.substring(start, end);
            if ( ! refSub.equals(potentialRepeat) )
                return false; // repeat didn't match, fail
        }

        return true; // we passed all tests, we matched
    }

    /**
     * Subset the samples in VC to reference only information with ref call alleles
     *
     * Preserves DP if present
     *
     * @param vc the variant context to subset down to
     * @param defaultPloidy defaultPloidy to use if a genotype doesn't have any alleles
     * @return a GenotypesContext
     */
    public static GenotypesContext subsetToRefOnly(final VariantContext vc, final int defaultPloidy) {
        Utils.nonNull(vc == null, "vc cannot be null");
        Utils.validateArg(defaultPloidy >= 1, () -> "defaultPloidy must be >= 1 but got " + defaultPloidy);

        final GenotypesContext oldGTs = vc.getGenotypes();
        if (oldGTs.isEmpty()) return oldGTs;
        final GenotypesContext newGTs = GenotypesContext.create(oldGTs.size());

        final Allele ref = vc.getReference();
        final List<Allele> diploidRefAlleles = Arrays.asList(ref, ref);

        // create the new genotypes
        for ( final Genotype g : vc.getGenotypes() ) {
            final int gPloidy = g.getPloidy() == 0 ? defaultPloidy : g.getPloidy();
            final List<Allele> refAlleles = gPloidy == 2 ? diploidRefAlleles : Collections.nCopies(gPloidy, ref);
            final GenotypeBuilder gb = new GenotypeBuilder(g.getSampleName(), refAlleles);
            if ( g.hasDP() ) gb.DP(g.getDP());
            if ( g.hasGQ() ) gb.GQ(g.getGQ());
            newGTs.add(gb.make());
        }

        return newGTs;
    }

    public static Genotype removePLsAndAD(final Genotype g) {
        return ( g.hasLikelihoods() || g.hasAD() ) ? new GenotypeBuilder(g).noPL().noAD().make() : g;
    }

    //TODO consider refactor variant-context merging code so that we share as much as possible between
    //TODO simpleMerge and referenceConfidenceMerge
    //TODO likely using a separate helper class or hierarchy.
    /**
     * Merges VariantContexts into a single hybrid.  Takes genotypes for common samples in priority order, if provided.
     * If uniquifySamples is true, the priority order is ignored and names are created by concatenating the VC name with
     * the sample name
     *
     * @param unsortedVCs               collection of unsorted VCs
     * @param priorityListOfVCs         priority list detailing the order in which we should grab the VCs
     * @param filteredRecordMergeType   merge type for filtered records
     * @param genotypeMergeOptions      merge option for genotypes
     * @param filteredAreUncalled       are filtered records uncalled?
     * @return new VariantContext       representing the merge of unsortedVCs
     */
    public static VariantContext simpleMerge(final Collection<VariantContext> unsortedVCs,
                                             final List<String> priorityListOfVCs,
                                             final FilteredRecordMergeType filteredRecordMergeType,
                                             final GenotypeMergeType genotypeMergeOptions,
                                             final boolean filteredAreUncalled) {
        int originalNumOfVCs = priorityListOfVCs == null ? 0 : priorityListOfVCs.size();
        return simpleMerge(unsortedVCs, priorityListOfVCs, originalNumOfVCs, filteredRecordMergeType, genotypeMergeOptions, filteredAreUncalled);
    }

    /**
     * Merges VariantContexts into a single hybrid.  Takes genotypes for common samples in priority order, if provided.
     * If uniquifySamples is true, the priority order is ignored and names are created by concatenating the VC name with
     * the sample name.
     * simpleMerge does not verify any more unique sample names EVEN if genotypeMergeOptions == GenotypeMergeType.REQUIRE_UNIQUE. One should use
     * SampleUtils.verifyUniqueSamplesNames to check that before using simpleMerge.
     *
     * For more information on this method see: http://www.thedistractionnetwork.com/programmer-problem/
     *
     * @param unsortedVCs               collection of unsorted VCs
     * @param priorityListOfVCs         priority list detailing the order in which we should grab the VCs
     * @param filteredRecordMergeType   merge type for filtered records
     * @param genotypeMergeOptions      merge option for genotypes
     * @param filteredAreUncalled       are filtered records uncalled?
     * @return new VariantContext       representing the merge of unsortedVCs
     */
    public static VariantContext simpleMerge(final Collection<VariantContext> unsortedVCs,
                                             final List<String> priorityListOfVCs,
                                             final int originalNumOfVCs,
                                             final FilteredRecordMergeType filteredRecordMergeType,
                                             final GenotypeMergeType genotypeMergeOptions,
                                             final boolean filteredAreUncalled) {
        if ( unsortedVCs == null || unsortedVCs.isEmpty() )
            return null;

        if (priorityListOfVCs != null && originalNumOfVCs != priorityListOfVCs.size())
            throw new IllegalArgumentException("the number of the original VariantContexts must be the same as the number of VariantContexts in the priority list");

        final List<VariantContext> preFilteredVCs = sortVariantContextsByPriority(unsortedVCs, priorityListOfVCs, genotypeMergeOptions);
        // Make sure all variant contexts are padded with reference base in case of indels if necessary
        List<VariantContext> VCs = preFilteredVCs.stream()
                .filter(vc -> !filteredAreUncalled || vc.isNotFiltered())
                .collect(Collectors.toList());

        if ( VCs.isEmpty() ) // everything is filtered out and we're filteredAreUncalled
            return null;

        // establish the baseline info from the first VC
        final VariantContext first = VCs.get(0);
        final String name = first.getSource();
        final Allele refAllele = determineReferenceAllele(VCs);

        final LinkedHashSet<Allele> alleles = new LinkedHashSet<>();
        final Set<String> filters = new LinkedHashSet<>();
        final Map<String, Object> attributes = new LinkedHashMap<>();
        final Set<String> inconsistentAttributes = new LinkedHashSet<>();
        final Set<String> variantSources = new LinkedHashSet<>(); // contains the set of sources we found in our set of VCs that are variant
        final Set<String> rsIDs = new LinkedHashSet<>(1); // most of the time there's one id

        VariantContext longestVC = first;
        int depth = 0;
        double log10PError = CommonInfo.NO_LOG10_PERROR;
        boolean anyVCHadFiltersApplied = false;
        GenotypesContext genotypes = GenotypesContext.create();

        // counting the number of filtered and variant VCs
        int nFiltered = 0;

        // cycle through and add info from the other VCs, making sure the loc/reference matches
        for ( final VariantContext vc : VCs ) {
            Utils.validate(longestVC.getStart() == vc.getStart(), () -> "BUG: attempting to merge VariantContexts with different start sites: first="+ first.toString() + " second=" + vc.toString());

            if ( VariantContextUtils.getSize(vc) > VariantContextUtils.getSize(longestVC) )
                longestVC = vc; // get the longest location

            nFiltered += vc.isFiltered() ? 1 : 0;
            if ( vc.isVariant() ) variantSources.add(vc.getSource());

            AlleleMapper alleleMapping = resolveIncompatibleAlleles(refAllele, vc);

            alleles.addAll(alleleMapping.values());

            mergeGenotypes(genotypes, vc, alleleMapping, genotypeMergeOptions == GenotypeMergeType.UNIQUIFY);

            // We always take the QUAL of the first VC with a non-MISSING qual for the combined value
            if ( log10PError == CommonInfo.NO_LOG10_PERROR )
                log10PError =  vc.getLog10PError();

            filters.addAll(vc.getFilters());
            anyVCHadFiltersApplied |= vc.filtersWereApplied();

            //
            // add attributes
            //
            // special case DP (add it up) and ID (just preserve it)
            //
            if (vc.hasAttribute(VCFConstants.DEPTH_KEY))
                depth += vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0);
            if ( vc.hasID() ) rsIDs.add(vc.getID());

            for (final Map.Entry<String, Object> p : vc.getAttributes().entrySet()) {
                final String key = p.getKey();
                final Object value = p.getValue();
                // only output annotations that have the same value in every input VC
                // if we don't like the key already, don't go anywhere
                if ( ! inconsistentAttributes.contains(key) ) {
                    final boolean alreadyFound = attributes.containsKey(key);
                    final Object boundValue = attributes.get(key);
                    final boolean boundIsMissingValue = alreadyFound && boundValue.equals(VCFConstants.MISSING_VALUE_v4);

                    if ( alreadyFound && ! boundValue.equals(value) && ! boundIsMissingValue ) {
                        // we found the value but we're inconsistent, put it in the exclude list
                        inconsistentAttributes.add(key);
                        attributes.remove(key);
                    } else if ( ! alreadyFound || boundIsMissingValue )  { // no value
                        attributes.put(key, value);
                    }
                }
            }
        }

        // if we have more alternate alleles in the merged VC than in one or more of the
        // original VCs, we need to strip out the GL/PLs (because they are no longer accurate), as well as allele-dependent attributes like AC,AF, and AD
        for ( final VariantContext vc : VCs ) {
            if (vc.getAlleles().size() == 1)
                continue;
            if ( hasPLIncompatibleAlleles(alleles, vc.getAlleles())) {
                if ( ! genotypes.isEmpty() ) {
                    logger.debug(String.format("Stripping PLs at %s:%d-%d due to incompatible alleles merged=%s vs. single=%s",
                            vc.getContig(), vc.getStart(), vc.getEnd(), alleles, vc.getAlleles()));
                }
                genotypes = stripPLsAndAD(genotypes);
                // this will remove stale AC,AF attributed from vc
                VariantContextUtils.calculateChromosomeCounts(vc, attributes, true);
                break;
            }
        }

        // if at least one record was unfiltered and we want a union, clear all of the filters
        if ( (filteredRecordMergeType == FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED && nFiltered != VCs.size()) || filteredRecordMergeType == FilteredRecordMergeType.KEEP_UNCONDITIONAL )
            filters.clear();

        if ( depth > 0 )
            attributes.put(VCFConstants.DEPTH_KEY, String.valueOf(depth));

        final String ID = rsIDs.isEmpty() ? VCFConstants.EMPTY_ID_FIELD : Utils.join(",", rsIDs);

        final VariantContextBuilder builder = new VariantContextBuilder().source(name).id(ID);
        builder.loc(longestVC.getContig(), longestVC.getStart(), longestVC.getEnd());
        builder.alleles(alleles);
        builder.genotypes(genotypes);
        builder.log10PError(log10PError);
        if ( anyVCHadFiltersApplied ) {
            builder.filters(filters.isEmpty() ? filters : new TreeSet<>(filters));
        }
        builder.attributes(new TreeMap<>(attributes));

        // Trim the padded bases of all alleles if necessary
        final VariantContext merged = builder.make();
        return merged;
    }


    //TODO as part of a larger refactoring effort remapAlleles can be merged with createAlleleMapping.

    public static GenotypesContext stripPLsAndAD(final GenotypesContext genotypes) {
        final GenotypesContext newGs = GenotypesContext.create(genotypes.size());
        for ( final Genotype g : genotypes ) {
            newGs.add(removePLsAndAD(g));
        }
        return newGs;
    }

    private static Allele determineReferenceAllele(final List<VariantContext> VCs) {
        return determineReferenceAllele(VCs, null);
    }

    public static boolean contextMatchesLoc(final VariantContext vc, final Locatable loc) {
        return loc == null || loc.getStart() == vc.getStart();
    }

    public static AlleleMapper resolveIncompatibleAlleles(final Allele refAllele, final VariantContext vc) {
        if ( refAllele.equals(vc.getReference()) )
            return new AlleleMapper(vc);
        else {
            final Map<Allele, Allele> map = createAlleleMapping(refAllele, vc);
            map.put(vc.getReference(), refAllele);
            return new AlleleMapper(map);
        }
    }

    public static Map<Allele, Allele> createAlleleMapping(final Allele refAllele, final VariantContext oneVC) {
        return createAlleleMapping(refAllele, oneVC.getReference(), oneVC.getAlternateAlleles());
    }

    //TODO as part of a larger refactoring effort {@link #createAlleleMapping} can be merged with {@link ReferenceConfidenceVariantContextMerger#remapAlleles}.
    /**
     * Create an allele mapping for the given context where its reference allele must (potentially) be extended to the given allele
     *
     * The refAllele is the longest reference allele seen at this start site.
     * So imagine it is:
     * refAllele: ACGTGA
     * myRef:     ACGT
     * myAlt:     A
     *
     * We need to remap all of the alleles in vc to include the extra GA so that
     * myRef => refAllele and myAlt => AGA
     *
     * @param refAllele          the new (extended) reference allele
     * @param inputRef           the reference allele that may need to be extended
     * @param inputAlts          the alternate alleles that may need to be extended
     * @return a non-null mapping of original alleles to new (extended) ones
     */
    public static Map<Allele, Allele> createAlleleMapping(final Allele refAllele, final Allele inputRef, final List<Allele> inputAlts) {
        Utils.validate( refAllele.length() >= inputRef.length(), () -> "BUG: inputRef="+inputRef+" is longer than refAllele="+refAllele);

        // frequent simple case where there is already a common reference
        if (refAllele.length() == inputRef.length()) {
            final Map<Allele, Allele> map = new LinkedHashMap<>();
            inputAlts.forEach(a -> map.put(a,a));
            return map;
        }

        final byte[] extraBases = Arrays.copyOfRange(refAllele.getBases(), inputRef.length(), refAllele.length());

        final Map<Allele, Allele> map = new LinkedHashMap<>();
        for ( final Allele a : inputAlts ) {
            if ( isNonSymbolicExtendableAllele(a) ) {
                Allele extended = Allele.extend(a, extraBases);
                map.put(a, extended);
            } else if (a.equals(Allele.SPAN_DEL)) {
                map.put(a, a);
            }
        }

        return map;
    }

    private static boolean isNonSymbolicExtendableAllele(final Allele allele) {
        return ! (allele.isReference() || allele.isSymbolic() || allele.equals(Allele.SPAN_DEL));
    }

    public static List<VariantContext> sortVariantContextsByPriority(Collection<VariantContext> unsortedVCs, List<String> priorityListOfVCs, GenotypeMergeType mergeOption ) {
        if ( mergeOption == GenotypeMergeType.PRIORITIZE && priorityListOfVCs == null )
            throw new IllegalArgumentException("Cannot merge calls by priority with a null priority list");

        if ( priorityListOfVCs == null || mergeOption == GenotypeMergeType.UNSORTED )
            return new ArrayList<>(unsortedVCs);
        else {
            ArrayList<VariantContext> sorted = new ArrayList<>(unsortedVCs);
            Collections.sort(sorted, new CompareByPriority(priorityListOfVCs));
            return sorted;
        }
    }

    private static void mergeGenotypes(GenotypesContext mergedGenotypes, VariantContext oneVC, AlleleMapper alleleMapping, boolean uniquifySamples) {
        //TODO: should we add a check for cases when the genotypeMergeOption is REQUIRE_UNIQUE
        for ( final Genotype g : oneVC.getGenotypes() ) {
            final String name = mergedSampleName(oneVC.getSource(), g.getSampleName(), uniquifySamples);
            if ( ! mergedGenotypes.containsSample(name) ) {
                // only add if the name is new
                Genotype newG = g;

                if ( uniquifySamples || alleleMapping.needsRemapping() ) {
                    final List<Allele> alleles = alleleMapping.needsRemapping() ? alleleMapping.remap(g.getAlleles()) : g.getAlleles();
                    newG = new GenotypeBuilder(g).name(name).alleles(alleles).make();
                }

                mergedGenotypes.add(newG);
            }
        }
    }

    /**
     * Cached NO_CALL immutable lists where the position ith contains the list with i elements.
     */
    @SuppressWarnings({"rawtypes", "unchecked"})
    private static List<Allele>[] NOCALL_LISTS = new List[] {
            Collections.emptyList(),
            Collections.singletonList(Allele.NO_CALL),
            Collections.nCopies(2,Allele.NO_CALL)
    };

    /**
     * Code to ensure that {@link #NOCALL_LISTS} has enough entries beyond the requested ploidy
     * @param capacity the requested ploidy.
     */
    private static void ensureNoCallListsCapacity(final int capacity) {
        final int currentCapacity = NOCALL_LISTS.length - 1;
        if (currentCapacity >= capacity)
            return;
        NOCALL_LISTS = Arrays.copyOf(NOCALL_LISTS,Math.max(capacity,currentCapacity << 1) + 1);
        for (int i = currentCapacity + 1; i < NOCALL_LISTS.length; i++)
            NOCALL_LISTS[i] = Collections.nCopies(i,Allele.NO_CALL);
    }

    /**
     * Returns a {@link Allele#NO_CALL NO_CALL} allele list provided the ploidy.
     *
     * @param ploidy the required ploidy.
     *
     * @return never {@code null}, but an empty list if {@code ploidy} is equal or less than 0. The returned list
     *   might or might not be mutable.
     */
    public static List<Allele> noCallAlleles(final int ploidy) {
        if (NOCALL_LISTS.length <= ploidy)
            ensureNoCallListsCapacity(ploidy);
        return NOCALL_LISTS[ploidy];
    }


    public static String mergedSampleName(String trackName, String sampleName, boolean uniquify ) {
        return uniquify ? sampleName + "." + trackName : sampleName;
    }

    /**
     * Trim the alleles in inputVC from the reverse direction
     *
     * @param inputVC a non-null input VC whose alleles might need a haircut
     * @return a non-null VariantContext (may be == to inputVC) with alleles trimmed up
     */
    public static VariantContext reverseTrimAlleles( final VariantContext inputVC ) {
        return trimAlleles(inputVC, false, true);
    }

    /**
     * Trim the alleles in inputVC forward and reverse, as requested
     *
     * @param inputVC a non-null input VC whose alleles might need a haircut
     * @param trimForward should we trim up the alleles from the forward direction?
     * @param trimReverse should we trim up the alleles from the reverse direction?
     * @return a non-null VariantContext (may be == to inputVC) with trimmed up alleles
     */
    public static VariantContext trimAlleles(final VariantContext inputVC, final boolean trimForward, final boolean trimReverse) {
        Utils.nonNull(inputVC);

        if ( inputVC.getNAlleles() <= 1 || inputVC.getAlleles().stream().anyMatch(a -> a.length() == 1 && !a.equals(Allele.SPAN_DEL)) ) {
            return inputVC;
        }

        final List<byte[]> sequences = inputVC.getAlleles().stream().filter(a -> !a.isSymbolic() && !a.equals(Allele.SPAN_DEL)).map(Allele::getBases).collect(Collectors.toList());
        final List<IndexRange> ranges = inputVC.getAlleles().stream().filter(a -> !a.isSymbolic() && !a.equals(Allele.SPAN_DEL)).map(a -> new IndexRange(0, a.length())).collect(Collectors.toList());

        final Pair<Integer, Integer> shifts = AlignmentUtils.normalizeAlleles(sequences, ranges, 0, true);
        final int endTrim = shifts.getRight();
        final int startTrim = -shifts.getLeft();

        final boolean emptyAllele = ranges.stream().anyMatch(r -> r.size() == 0);
        final boolean restoreOneBaseAtEnd = emptyAllele && startTrim == 0;
        final boolean restoreOneBaseAtStart = emptyAllele && startTrim > 0;

        // if the end trimming consumed all the bases, leave one base
        final int endBasesToClip = restoreOneBaseAtEnd ? endTrim - 1 : endTrim;
        final int startBasesToClip = restoreOneBaseAtStart ? startTrim - 1 : startTrim;

        return trimAlleles(inputVC, (trimForward ? startBasesToClip : 0) - 1, trimReverse ? endBasesToClip : 0);
    }

    /**
     * Trim up alleles in inputVC, cutting out all bases up to fwdTrimEnd inclusive and
     * the last revTrim bases from the end
     *
     * @param inputVC a non-null input VC
     * @param fwdTrimEnd bases up to this index (can be -1) will be removed from the start of all alleles
     * @param revTrim the last revTrim bases of each allele will be clipped off as well
     * @return a non-null VariantContext (may be == to inputVC) with trimmed up alleles
     */
    protected static VariantContext trimAlleles(final VariantContext inputVC,
                                                final int fwdTrimEnd,
                                                final int revTrim) {
        if( fwdTrimEnd == -1 && revTrim == 0 ) // nothing to do, so just return inputVC unmodified
            return inputVC;

        final List<Allele> alleles = new LinkedList<>();
        final Map<Allele, Allele> originalToTrimmedAlleleMap = new LinkedHashMap<>();

        for (final Allele a : inputVC.getAlleles()) {
            if (a.isSymbolic() || a.equals(Allele.SPAN_DEL)) {
                alleles.add(a);
                originalToTrimmedAlleleMap.put(a, a);
            } else {
                // get bases for current allele and create a new one with trimmed bases
                final byte[] newBases = Arrays.copyOfRange(a.getBases(), fwdTrimEnd+1, a.length()-revTrim);
                final Allele trimmedAllele = Allele.create(newBases, a.isReference());
                alleles.add(trimmedAllele);
                originalToTrimmedAlleleMap.put(a, trimmedAllele);
            }
        }

        // now we can recreate new genotypes with trimmed alleles
        final AlleleMapper alleleMapper = new AlleleMapper(originalToTrimmedAlleleMap);
        final GenotypesContext genotypes = updateGenotypesWithMappedAlleles(inputVC.getGenotypes(), alleleMapper);

        final int start = inputVC.getStart() + (fwdTrimEnd + 1);
        final VariantContextBuilder builder = new VariantContextBuilder(inputVC);
        builder.start(start);
        builder.stop(start + alleles.get(0).length() - 1);
        builder.alleles(alleles);
        builder.genotypes(genotypes);
        return builder.make();
    }

    protected static GenotypesContext updateGenotypesWithMappedAlleles(final GenotypesContext originalGenotypes, final AlleleMapper alleleMapper) {
        final GenotypesContext updatedGenotypes = GenotypesContext.create(originalGenotypes.size());

        for ( final Genotype genotype : originalGenotypes ) {
            final List<Allele> updatedAlleles = alleleMapper.remap(genotype.getAlleles());
            updatedGenotypes.add(new GenotypeBuilder(genotype).alleles(updatedAlleles).make());
        }

        return updatedGenotypes;
    }

    protected static class AlleleMapper {
        private VariantContext vc = null;
        private Map<Allele, Allele> map = null;
        public AlleleMapper(VariantContext vc)          { this.vc = vc; }
        public AlleleMapper(Map<Allele, Allele> map)    { this.map = map; }
        public boolean needsRemapping()                 { return this.map != null; }
        public Collection<Allele> values()              { return map != null ? map.values() : vc.getAlleles(); }
        public Allele remap(Allele a)                   { return map != null && map.containsKey(a) ? map.get(a) : a; }

        public List<Allele> remap(List<Allele> as) {
            List<Allele> newAs = as.stream()
                    .map(this::remap)
                    .collect(Collectors.toList());
            return newAs;
        }

    }

    private static class CompareByPriority implements Comparator<VariantContext>, Serializable {
        private static final long serialVersionUID = 0L;

        List<String> priorityListOfVCs;
        public CompareByPriority(List<String> priorityListOfVCs) {
            this.priorityListOfVCs = priorityListOfVCs;
        }

        private int getIndex(VariantContext vc) {
            return Utils.validIndex(priorityListOfVCs.indexOf(vc.getSource()), priorityListOfVCs.size());
        }

        @Override
        public int compare(VariantContext vc1, VariantContext vc2) {
            return Integer.valueOf(getIndex(vc1)).compareTo(getIndex(vc2));
        }
    }

    /**
     * For testing purposes only.  Create a site-only VariantContext at contig:start containing alleles
     *
     * @param name the name of the VC
     * @param contig the contig for the VC
     * @param start the start of the VC
     * @param alleleStrings a non-null, non-empty list of strings for the alleles.  The first will be the ref allele, and others the
     *                      alt.  Will compute the stop of the VC from the length of the reference allele
     * @return a non-null VariantContext
     */
    public static VariantContext makeFromAlleles(final String name, final String contig, final int start, final List<String> alleleStrings) {
        Utils.nonNull(alleleStrings, "alleleStrings is null");
        Utils.nonEmpty(alleleStrings, "alleleStrings is empty");

        final List<Allele> alleles = new LinkedList<>();
        final int length = alleleStrings.get(0).length();

        boolean first = true;
        for ( final String alleleString : alleleStrings ) {
            alleles.add(Allele.create(alleleString, first));
            first = false;
        }
        return new VariantContextBuilder(name, contig, start, start+length-1, alleles).make();
    }

    public static List<VariantContext> splitSomaticVariantContextToBiallelics(final VariantContext vc, final boolean trimLeft, final VCFHeader outputHeader) {
        Utils.nonNull(vc);

        if (!vc.isVariant() || vc.isBiallelic()) {
            // non variant or biallelics already satisfy the contract
            return Collections.singletonList(vc);
        } else {
            final List<VariantContext> biallelics = new LinkedList<>();

            List<String> attrsSpecialFormats = new ArrayList<String>(Arrays.asList(GATKVCFConstants.AS_FILTER_STATUS_KEY, GATKVCFConstants.AS_SB_TABLE_KEY));
            List<Map<String, Object>> attributesByAllele = splitAttributesIntoPerAlleleLists(vc, attrsSpecialFormats, outputHeader);
            splitASSBTable(vc, attributesByAllele);
            splitASFilters(vc, attributesByAllele);

            ListIterator<Map<String, Object>> attributesByAlleleIterator = attributesByAllele.listIterator();

            for (final Allele alt : vc.getAlternateAlleles()) {
                final VariantContextBuilder builder = new VariantContextBuilder(vc);

                // make biallelic alleles
                final List<Allele> alleles = Arrays.asList(vc.getReference(), alt);
                builder.alleles(alleles);
                Map<String, Object> attributes = attributesByAlleleIterator.next();
                builder.attributes(attributes);

                // now add the allele specific filters to the variant context
                String filters = (String) attributes.get(GATKVCFConstants.AS_FILTER_STATUS_KEY);
                // checking for . and PASS should be able to be removed. these were temporarily used to indicate no allele specific filter
                if (filters != null && !filters.isEmpty() && !filters.equals(VCFConstants.EMPTY_INFO_FIELD) && !filters.equals(GATKVCFConstants.SITE_LEVEL_FILTERS) && !filters.equals((VCFConstants.PASSES_FILTERS_v4))) {
                    AnnotationUtils.decodeAnyASList(filters).stream().forEach(filter -> builder.filter(filter));
                }

                int alleleIndex = vc.getAlleleIndex(alt);
                builder.genotypes(AlleleSubsettingUtils.subsetSomaticAlleles(outputHeader, vc.getGenotypes(), alleles, new int[]{0, alleleIndex}));
                final VariantContext trimmed = trimAlleles(builder.make(), trimLeft, true);
                biallelics.add(trimmed);
            }
            return biallelics;
        }
    }

    public static void splitASSBTable(VariantContext vc, List<Map<String, Object>> attrsByAllele) {
        List<String> sbs = StrandBiasUtils.getSBsForAlleles(vc).stream().map(ints -> StrandBiasUtils.encode(ints)).collect(Collectors.toList());
        new IndexRange(1, sbs.size()).forEach(i -> {
            String newattrs = String.join(AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM, new ArrayList<String>(Arrays.asList(sbs.get(0), sbs.get(i))));
            attrsByAllele.get(i - 1).put(GATKVCFConstants.AS_SB_TABLE_KEY, newattrs);
        });
    }

    public static void splitASFilters(VariantContext vc, List<Map<String, Object>> attrsByAllele) {
        // the reason we are getting as list and then joining on , is because the default getAttributeAsString for a list will add spaces between items which we don't
        // want to have to trim out later in the code
        String asfiltersStr = String.join(",", vc.getCommonInfo().getAttributeAsStringList(GATKVCFConstants.AS_FILTER_STATUS_KEY, GATKVCFConstants.SITE_LEVEL_FILTERS));
        List<String> filtersList = AnnotationUtils.decodeAnyASListWithRawDelim(asfiltersStr);
        new IndexRange(0, filtersList.size()).forEach(i -> attrsByAllele.get(i).put(GATKVCFConstants.AS_FILTER_STATUS_KEY, filtersList.get(i)));
    }

    public static List<Map<String, Object>> splitAttributesIntoPerAlleleLists(VariantContext vc, List<String> skipAttributes, VCFHeader outputHeader) {
        List<Map<String, Object>> results = new ArrayList<>(vc.getNAlleles()-1);
        vc.getAlternateAlleles().forEach(alt -> results.add(new HashMap<>()));

        Map<String, Object> attributes = vc.getAttributes();
        attributes.entrySet().stream().filter(entry -> !skipAttributes.contains(entry.getKey())).forEachOrdered(entry -> {
            String key = entry.getKey();
            // default to unbounded in case header is not found
            VCFHeaderLineCount countType = VCFHeaderLineCount.UNBOUNDED;
            try {
                VCFInfoHeaderLine header = outputHeader.getInfoHeaderLine(key);
                countType = header.getCountType();
            } catch (IllegalStateException ex) {
                // this happens for DP if we use GATKVCFHeaderLines.getInfoLine(key)
                // shouldn't happen now that we use the generated output header
                logger.warn("Could not find header info for key " + key);
            }
            // override count type for this attribute
            if (key.equals(GATKVCFConstants.REPEATS_PER_ALLELE_KEY)) {
                countType = VCFHeaderLineCount.R;
            }
            List<Object> attr;
            switch (countType) {
                case A:
                    attr = vc.getCommonInfo().getAttributeAsList(key);
                    ValidationUtils.validateArg(attr.size() == results.size(), "Incorrect attribute size for " + key);
                    new IndexRange(0, attr.size()).forEach(i -> results.get(i).put(key, attr.get(i)));
                    break;
                case R:
                    attr = vc.getCommonInfo().getAttributeAsList(key);
                    ValidationUtils.validateArg(attr.size() == vc.getNAlleles(), "Incorrect attribute size for " + key);
                    new IndexRange(1, attr.size()).forEach(i -> {
                        List<Object> newattrs = new ArrayList<Object>(Arrays.asList(attr.get(0), attr.get(i)));
                        results.get(i-1).put(key, newattrs);
                    });
                    break;
                default:
                    results.forEach(altMap -> altMap.put(key, entry.getValue()));

            }

        });
        return results;
    }

    /**
     * Split variant context into its biallelic components if there are more than 2 alleles
     * <p>
     * For VC has A/B/C alleles, returns A/B and A/C contexts.
     * Alleles are right trimmed to satisfy VCF conventions
     * <p>
     * If vc is biallelic or non-variant it is just returned
     * <p>
     * Chromosome counts are updated (but they are by definition 0)
     *
     * @param vc                       a potentially multi-allelic variant context
     * @param trimLeft                 if true, we will also left trim alleles, potentially moving the resulting vcs forward on the genome
     * @param genotypeAssignmentMethod assignment strategy for the (subsetted) PLs
     * @param keepOriginalChrCounts    keep the orignal chromosome counts before subsetting
     * @return a list of bi-allelic (or monomorphic) variant context
     */
    public static List<VariantContext> splitVariantContextToBiallelics(final VariantContext vc, final boolean trimLeft, final GenotypeAssignmentMethod genotypeAssignmentMethod,
                                                          final boolean keepOriginalChrCounts) {
        Utils.nonNull(vc);
        if (!vc.isVariant()) {
            return Collections.emptyList();
        } else if (vc.isBiallelic())
            // non variant or biallelics already satisfy the contract
            return Collections.singletonList(vc);
        else {
            final List<VariantContext> biallelics = new LinkedList<>();

            // if any of the genotypes are het-not-ref (i.e. 1/2), set all of them to no-call
            final GenotypeAssignmentMethod genotypeAssignmentMethodUsed = hasHetNonRef(vc.getGenotypes()) ? GenotypeAssignmentMethod.SET_TO_NO_CALL_NO_ANNOTATIONS : genotypeAssignmentMethod;

            for (final Allele alt : vc.getAlternateAlleles()) {
                final VariantContextBuilder builder = new VariantContextBuilder(vc);

                // make biallelic alleles
                final List<Allele> alleles = Arrays.asList(vc.getReference(), alt);
                builder.alleles(alleles);

                // since the VC has been subset, remove the invalid attributes
                for (final String key : vc.getAttributes().keySet()) {
                    if (!(key.equals(VCFConstants.ALLELE_COUNT_KEY) || key.equals(VCFConstants.ALLELE_FREQUENCY_KEY) || key.equals(VCFConstants.ALLELE_NUMBER_KEY)) ||
                            genotypeAssignmentMethodUsed == GenotypeAssignmentMethod.SET_TO_NO_CALL_NO_ANNOTATIONS) {
                        builder.rmAttribute(key);
                    }
                }

                // subset INFO field annotations if available if genotype is called
                if (genotypeAssignmentMethodUsed != GenotypeAssignmentMethod.SET_TO_NO_CALL_NO_ANNOTATIONS &&
                        genotypeAssignmentMethodUsed != GenotypeAssignmentMethod.SET_TO_NO_CALL)
                    AlleleSubsettingUtils.addInfoFieldAnnotations(vc, builder, keepOriginalChrCounts);

                builder.genotypes(AlleleSubsettingUtils.subsetAlleles(vc.getGenotypes(),2,vc.getAlleles(), alleles, null, genotypeAssignmentMethodUsed));
                final VariantContext trimmed = trimAlleles(builder.make(), trimLeft, true);
                biallelics.add(trimmed);
            }

            return biallelics;
        }
    }

    /**
     * Split variant context into its biallelic components if there are more than 2 alleles
     * <p>
     * For VC has A/B/C alleles, returns A/B and A/C contexts.
     * Alleles are right trimmed to satisfy VCF conventions
     * <p>
     * If vc is biallelic or non-variant it is just returned
     * <p>
     * Chromosome counts are updated (but they are by definition 0)
     *
     * @param vc                       a potentially multi-allelic variant context
     * @param trimLeft                 if true, we will also left trim alleles, potentially moving the resulting vcs forward on the genome
     * @param genotypeAssignmentMethod assignment strategy for the (subsetted) PLs
     * @param keepOriginalChrCounts    keep the orignal chromosome counts before subsetting
     * @return a list of bi-allelic (or monomorphic) variant context
     */
    public static List<Event> splitVariantContextToEvents(final VariantContext vc, final boolean trimLeft, final GenotypeAssignmentMethod genotypeAssignmentMethod,
                                                          final boolean keepOriginalChrCounts) {

            return splitVariantContextToBiallelics(vc, trimLeft, genotypeAssignmentMethod, keepOriginalChrCounts).stream()
                    .map(Event::ofWithoutAttributes).collect(Collectors.toList());
    }

    /**
     * Check if any of the genotypes is heterozygous, non-reference (i.e. 1/2)
     *
     * @param genotypesContext genotype information
     * @return true if any of the genotypes are heterozygous, non-reference, false otherwise
     */
    private static boolean hasHetNonRef(final GenotypesContext genotypesContext) {
        for (final Genotype gt : genotypesContext) {
            if (gt.isHetNonRef())
                return true;
        }
        return false;
    }

    /**
     * Add chromosome counts (AC, AN and AF) to the VCF header lines
     *
     * @param headerLines the VCF header lines
     */
    public static void addChromosomeCountsToHeader(final Set<VCFHeaderLine> headerLines) {
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY));
    }

    /**
     * Set the builder's filtered genotypes to no-call and update AC, AN and AF
     *
     * @param builder builder for variant context
     * @param vc the VariantContext record to set filtered genotypes to no-call
     * @param setFilteredGenotypesToNocall flag to set filtered genotype to NO CALL
     * @param filters the filters for each genotype
     */
    public static void setFilteredGenotypeToNocall(final VariantContextBuilder builder, final VariantContext vc,
                                                   final boolean setFilteredGenotypesToNocall, BiFunction<VariantContext, Genotype, List<String>> filters) {
        Utils.nonNull(vc);
        Utils.nonNull(builder);
        Utils.nonNull(filters);

        final GenotypesContext genotypes = GenotypesContext.create(vc.getGenotypes().size());

        // update genotypes if filtered genotypes are set to no-call
        boolean haveFilteredNoCallAlleles = false;
        for (final Genotype g : vc.getGenotypes()) {
            if (g.isCalled()) {
                List<String> filterNames = filters.apply(vc,g);
                if (!filterNames.isEmpty() && setFilteredGenotypesToNocall) {
                    haveFilteredNoCallAlleles = true;
                    genotypes.add(new GenotypeBuilder(g).filters(filterNames).alleles(GATKVariantContextUtils.noCallAlleles((g.getPloidy()))).make());
                } else {
                    genotypes.add(new GenotypeBuilder(g).filters(filterNames).make());
                }
            } else {
                genotypes.add(g);
            }
        }

        // if filtered genotypes are set to no-call, output recomputed AC, AN, AF
        if ( haveFilteredNoCallAlleles ) {
            final Map<String, Object> attributes = new LinkedHashMap<>(vc.getAttributes()); // need to make mutable
            VariantContextUtils.calculateChromosomeCounts(builder.genotypes(genotypes).make(), attributes, true, vc.getSampleNames());
            builder.attributes(attributes);
        }

        // update genotypes
        builder.genotypes(genotypes);
    }

    /**
     * Calculate the Genotype Quality (GQ) by subtracting the smallest Phred-scaled likelihoods (PL) from the second smallest.
     *
     * @param plValues  array of PL values
     * @return          the genotype quality corresponding to the PL values
     */
    public static int calculateGQFromPLs(final int[] plValues) {
        Utils.nonNull(plValues);
        Utils.validateArg(plValues.length >= 2, () -> "Array of PL values must contain at least two elements.");

        int first = plValues[0];
        int second = plValues[1];
        if (first > second) {
            second = first;
            first = plValues[1];
        }
        for (int i = 2; i < plValues.length; i++) {
            final int candidate = plValues[i];
            if (candidate >= second) {
                continue;
            }
            if (candidate <= first) {
                second = first;
                first = candidate;
            } else
                second = candidate;
        }
        return second - first;
    }

    /**
     *  Create a string with existing filters plus the one to append
     * @param vc VariantContext to base initial filter values.  Not {@code null}
     * @param filterToAppend the filter value to append to the strings Not {@code null}
     * @return String of filter values.  Sorted alphabetically.
     */
    public static List<String> createFilterListWithAppend(VariantContext vc, String filterToAppend) {
        Utils.nonNull(vc);
        Utils.nonNull(filterToAppend);

        final List<String> filtersAsList = new ArrayList<>(vc.getFilters());
        filtersAsList.add(filterToAppend);
        filtersAsList.sort(String::compareToIgnoreCase);
        return filtersAsList;
    }

    /**
     * Determines whether the given reference and alternate alleles constitute a frameshift mutation.
     * @param reference The reference {@link Allele}.  Must not be {@code null}.
     * @param alternate The alternate / variant {@link Allele}.  Must not be {@code null}.
     * @return {@code true} if replacing the reference with the alternate results in a frameshift.  {@code false} otherwise.
     */
    public static boolean isFrameshift(final Allele reference, final Allele alternate) {

        Utils.nonNull(reference);
        Utils.nonNull(alternate);

        // We know it's a frameshift if we have a replacement that is not of a
        // length evenly divisible by 3 because that's how many bases are read at once:
        return ((Math.abs( reference.length() - alternate.length() ) % 3) != 0);
    }

    /**
     * Determines whether the given reference and alternate alleles constitute a frameshift mutation.
     * @param reference The {@link String} containing the bases of the reference allele.  Must not be {@code null}.
     * @param alternate The {@link String} containing the bases of the alternate allele.  Must not be {@code null}.
     * @return {@code true} if replacing the reference with the alternate results in a frameshift.  {@code false} otherwise.
     */
    public static boolean isFrameshift(final String reference, final String alternate) {

        Utils.nonNull(reference);
        Utils.nonNull(alternate);

        final String refComparable = reference.replaceAll("\\*", "");
        final String altComperable = alternate.replaceAll("\\*", "");

        // We know it's a frameshift if we have a replacement that is not of a
        // length evenly divisible by 3 because that's how many bases are read at once:
        return ((Math.abs( refComparable.length() - altComperable.length() ) % 3) != 0);
    }

    /**
     * Determines whether the given reference and alternate alleles constitute a frameshift mutation.
     * This does not take into account the cases where either or both of the alleles overlap splice sites.
     * @param startPos Genomic start position (1-based, inclusive) of the variant.  Must be > 0.
     * @param refEnd Genomic end position (1-based, inclusive) of the reference allele.  Must be > 0.
     * @param altEnd Genomic end position (1-based, inclusive) of the alternate allele.  Must be > 0.
     * @return {@code true} if replacing the reference with the alternate results in a frameshift.  {@code false} otherwise.
     */
    public static boolean isFrameshift(final int startPos, final int refEnd, final int altEnd) {

        ParamUtils.isPositive( startPos, "Genomic positions must be > 0." );
        ParamUtils.isPositive( refEnd, "Genomic positions must be > 0." );
        ParamUtils.isPositive( altEnd, "Genomic positions must be > 0." );

        final int refLength = refEnd - startPos + 1;
        final int altLength = altEnd - startPos + 1;

        // We know it's a frameshift if we have a replacement that is not of a
        // length evenly divisible by 3 because that's how many bases are read at once:
        return ((Math.abs( refLength - altLength ) % 3) != 0);
    }

    /**
     * Determines whether the given reference and alternate alleles constitute an insertion mutation.
     * @param reference The reference {@link Allele}.   Must not be {@code null}.
     * @param alternate The alternate / variant {@link Allele}.   Must not be {@code null}.
     * @return {@code true} if replacing the reference with the alternate results in an insertion.  {@code false} otherwise.
     */
    public static boolean isInsertion(final Allele reference, final Allele alternate) {

        Utils.nonNull(reference);
        Utils.nonNull(alternate);

        // If we have more bases in the alternate, we have an insertion:
        return reference.length() < alternate.length();
    }

    /**
     * Determines whether the given reference and alternate alleles constitute an insertion mutation.
     * @param reference A {@link String} containing the bases of the reference allele.  Must not be {@code null}.
     * @param alternate A {@link String} containing the bases of the alternate / variant allele.  Must not be {@code null}.
     * @return {@code true} if replacing the reference with the alternate results in an insertion.  {@code false} otherwise.
     */
    public static boolean isInsertion(final String reference, final String alternate) {

        Utils.nonNull(reference);
        Utils.nonNull(alternate);

        final String refComparable = reference.replaceAll("\\*", "");
        final String altComperable = alternate.replaceAll("\\*", "");

        // If we have more bases in the alternate, we have an insertion:
        return refComparable.length() < altComperable.length();
    }

    /**
     * Determines whether the given reference and alternate alleles constitute a deletion mutation.
     * @param reference A {@link String} containing the bases of the reference allele.  Must not be {@code null}.
     * @param alternate A {@link String} containing the bases of the alternate / variant allele.  Must not be {@code null}.
     * @return {@code true} if replacing the reference with the alternate results in a deletion.  {@code false} otherwise.
     */
    public static boolean isDeletion(final String reference, final String alternate) {

        Utils.nonNull(reference);
        Utils.nonNull(alternate);

        final String refComparable = reference.replaceAll("\\*", "");
        final String altComperable = alternate.replaceAll("\\*", "");

        // If we have fewer bases in the alternate, we have a deletion:
        return refComparable.length() > altComperable.length();
    }

    /**
     * Determines whether the given reference and alternate alleles constitute a deletion mutation.
     * @param reference The reference {@link Allele}.  Must not be {@code null}.
     * @param alternate The alternate / variant {@link Allele}.  Must not be {@code null}.
     * @return {@code true} if replacing the reference with the alternate results in a deletion.  {@code false} otherwise.
     */
    public static boolean isDeletion(final Allele reference, final Allele alternate) {

        Utils.nonNull(reference);
        Utils.nonNull(alternate);

        // If we have fewer bases in the alternate, we have a deletion:
        return reference.length() > alternate.length();
    }

    /**
     * Determines whether the given reference and alternate alleles constitute an insertion or deletion mutation.
     * @param reference A {@link String} containing the bases of the reference allele.  Must not be {@code null}.
     * @param alternate A {@link String} containing the bases of the alternate / variant allele.  Must not be {@code null}.
     * @return {@code true} if replacing the reference with the alternate results in an insertion or deletion.  {@code false} otherwise.
     */
    public static boolean isIndel(final String reference, final String alternate) {

        Utils.nonNull(reference);
        Utils.nonNull(alternate);

        final String refComparable = reference.replaceAll("\\*", "");
        final String altComperable = alternate.replaceAll("\\*", "");

        // If we do not have the same number of bases in the reference and alternate alleles,
        // then we have an indel:
        return refComparable.length() != altComperable.length();
    }

    /**
     * Determines whether the given reference and alternate alleles constitute an insertion or deletion mutation.
     * @param reference The reference {@link Allele}.  Must not be {@code null}.
     * @param alternate The alternate / variant {@link Allele}.  Must not be {@code null}.
     * @return {@code true} if replacing the reference with the alternate results in an insertion or deletion.  {@code false} otherwise.
     */
    public static boolean isIndel(final Allele reference, final Allele alternate) {

        Utils.nonNull(reference);
        Utils.nonNull(alternate);

        // If we do not have the same number of bases in the reference and alternate alleles,
        // then we have an indel:
        return reference.length() != alternate.length();
    }

    /**
     * Determines whether the given reference and alternate alleles constitute a polymorphism in one or more nucleotides (XNP).
     * @param reference The reference {@link Allele}.  Must not be {@code null}.
     * @param alternate The alternate / variant {@link Allele}.  Must not be {@code null}.
     * @return {@code true} if replacing the reference with the alternate results in an XNP.  {@code false} otherwise.
     */
    public static boolean isXnp(final Allele reference, final Allele alternate) {

        Utils.nonNull(reference);
        Utils.nonNull(alternate);

        // If we have an equal number of bases in the reference and the alternate, we have an ONP:
        return (reference.length() == alternate.length()) && (!reference.equals(alternate));
    }

    /**
     * Determines whether the given reference and alternate alleles constitute a polymorphism in one or more nucleotides (XNP).
     * @param reference A {@link String} containing the bases of the reference allele.  Must not be {@code null}.
     * @param alternate A {@link String} containing the bases of the alternate / variant allele.  Must not be {@code null}.
     * @return {@code true} if replacing the reference with the alternate results in an XNP.  {@code false} otherwise.
     */
    public static boolean isXnp(final String reference, final String alternate) {

        Utils.nonNull(reference);
        Utils.nonNull(alternate);

        final String refComparable = reference.replaceAll("\\*", "");
        final String altComperable = alternate.replaceAll("\\*", "");

        // If we have an equal number of bases in the reference and the alternate, we have an ONP:
        return ((refComparable.length() == altComperable.length()) && (!refComparable.equals(altComperable)));
    }

    /**
     * @param vc {@link VariantContext to test}
     * @return true if the only alternate allele for this VariantContext is a spanning deletion, otherwise false.
     */
    public static boolean isSpanningDeletionOnly(final VariantContext vc){
        return vc.getAlternateAlleles().size() == 1 && GATKVCFConstants.isSpanningDeletion(vc.getAlternateAllele(0));
    }

    /**
     *
     * Attempt to match allele ref/alt pairs, even if the allele pairs in the given variant contexts are equivalent,
     *  but not exact.
     *
     * For example, if variant 1 and variant2 have the same position, but
     * Variant 1: "A", T", "C"
     * Variant 2: "ACCAGGCCCAGCTCATGCTTCTTTGCAGCCTCT", "TCCAGGCCCAGCTCATGCTTCTTTGCAGCCTCT", "A", "AC"
     *
     * Then the returned array would be:  {0, -1}
     * Since A>T matches in variant1 matches the the first alt allele in variant 2.  And A>C does not match any allele
     *  in variant 2.
     *
     * Do not use this method for doing a full split of a variant context into biallelic components.  This method
     *  ignores (and drops) genotypes and INFO attributes.  It is really meant just for alleles, but works on a
     *  VariantContext to better leverage existing code in
     *  {@link GATKVariantContextUtils#trimAlleles(VariantContext, boolean, boolean)}
     *
     * @param variant1 Never {@code null}
     * @param variant2 Never {@code null}
     * @return Returns indexes into the alternate alleles of variant2.  Note that this method assumes that (when biallelic) the variant
     *  contexts are already trimmed.  See {@link GATKVariantContextUtils#trimAlleles(VariantContext, boolean, boolean)}
     *  Never {@code null}.  The length will be equal to the number of alternate alleles in variant1.  A value of -1
     *  indicates no match in variant2.  If the reference alleles do not match, the output array will be populated
     *  exclusively with -1.
     */
    public static int[] matchAllelesAndStart(final VariantContext variant1, final VariantContext variant2) {
        Utils.nonNull(variant1);
        Utils.nonNull(variant2);

        // First split and trim all variant contexts into events
        final List<Event> events1 = simpleSplitIntoBiallelics(variant1);
        final List<Event> events2 = simpleSplitIntoBiallelics(variant2);

        // result[i] = j such that variant1's ith alt and variant2's jth alt are equivalent, or -1 if no such j exists
        return IntStream.range(0, events1.size())
                .map(i -> IntStream.range(0, events2.size()).filter(j -> events1.get(i).equals(events2.get(j))).findFirst().orElse(-1))
                .toArray();
    }

    /**
     * Do not use this method for doing a full split of a variant context into biallelic components.  This method
     *  ignores (and drops) genotypes and INFO attributes.  It is really meant just for alleles, but works on a
     *  VariantContext to better leverage existing code in
     *  {@link GATKVariantContextUtils#trimAlleles(VariantContext, boolean, boolean)}
     *
     * This method is trying to be fast, otherwise.
     *
     * NOTE: this will attempt to trim alleles to a minimal representation for all variant contexts.
     *
     * @param vc variant context to split into simple biallelics.  Never {@code null}
     * @return a list of variant contexts.  Each will be biallelic.  Length will be the number of alt alleles in the input vc.
     * Note that the variant contexts are usually stripped of attributes and genotypes.  Never {@code null}.  Empty list
     * if variant context has no alt alleles.
     */
    private static List<Event> simpleSplitIntoBiallelics(final VariantContext vc) {
        Utils.nonNull(vc);
        final List<Event> result = new ArrayList<>();

        if (vc.isBiallelic()) {
            return Collections.singletonList(Event.ofWithoutAttributes(GATKVariantContextUtils.trimAlleles(vc, true, true)));
        } else {
            // Since variant context builders are slow to keep re-creating.  Just create one and spew variant contexts from it, since
            //  the only difference will be the alternate allele.  Initialize the VCB with a dummy alternate allele,
            //  since it will be overwritten in all cases.
            final VariantContextBuilder vcb = new VariantContextBuilder("SimpleSplit", vc.getContig(), vc.getStart(), vc.getEnd(),
                    Arrays.asList(vc.getReference(), Allele.NO_CALL));
            vc.getAlternateAlleles().forEach(allele -> result.add(Event.ofWithoutAttributes(GATKVariantContextUtils.trimAlleles(
                    vcb.alleles(Arrays.asList(vc.getReference(), allele)).make(true), true, true))));
        }

        return result;
    }

    /**
     * Returns true if the context represents a MNP. If the context supplied contains the GVCF NON_REF symbolic
     * allele, then the determination is performed after discarding the NON_REF symbolic allele.
     *
     * @param vc a Variant context
     * @return true if the context represents an unmixed MNP (i.e. all alleles are non-symbolic and of the same
     *         length > 1), false otherwise
     */
    public static boolean isUnmixedMnpIgnoringNonRef(final VariantContext vc) {
        final List<Allele> alleles = vc.getAlleles();
        int length = vc.getReference().length(); // ref allele cannot be symbolic
        if (length < 2) { return false; }

        for (final Allele a : alleles) {
            if (a.isSymbolic() && !Allele.NON_REF_ALLELE.equals(a)) { return false; }
            else if (!a.isSymbolic() && a.length() != length) { return false; }
        }

        return true;
    }

    public static <T> List<T> removeDataForSymbolicAltAlleles(VariantContext vc, List<T> data) {
        return removeDataForSymbolicAlleles(vc, data, false);
    }

    public static <T> List<T> removeDataForSymbolicAlleles(VariantContext vc, List<T> data) {
        return removeDataForSymbolicAlleles(vc, data, true);
    }

    protected static <T> List<T> removeDataForSymbolicAlleles(VariantContext vc, List<T> data, boolean dataContainsReference) {
        if (vc.hasSymbolicAlleles()) {
            List<Allele> symbolicAlleles = vc.getAlternateAlleles().stream().filter(allele -> allele.isSymbolic()).collect(Collectors.toList());
            // convert allele index to index for data
            int offset = dataContainsReference ? 0 : 1;
            List<Integer> symAltIndexes = vc.getAlleleIndices(symbolicAlleles).stream().map(i -> i-offset).collect(Collectors.toList());
            return removeItemsByIndex(data, symAltIndexes);
        } else {
            return data;
        }
    }

    public static <T> List<T> removeItemsByIndex(List<T> data, List<Integer> indexesToRemove) {
        List<T> updated = new ArrayList<>();
        new IndexRange(0, data.size()).forEach(i -> {
            if (!indexesToRemove.contains(i)) {
                updated.add(data.get(i));
            }
        });
        return updated;
    }

    /**
     * Main routine workhorse for {@link org.broadinstitute.hellbender.tools.walkers.variantutils.LeftAlignAndTrimVariants}.
     * By definition, it will only take biallelic vc's. Splitting into multiple alleles has to be
     * handled by calling routine.
     *
     * @param vc  Input VC with variants to left align
     * @param ref Reference context
     * @return new VC.
     */
    public static VariantContext leftAlignAndTrim(final VariantContext vc, final ReferenceContext ref, final int maxLeadingBases, final boolean trim) {
        if (!vc.isIndel() || maxLeadingBases <= 0) {
            return vc;
        }


        for(int leadingBases = Math.min(maxLeadingBases, 10); leadingBases <= maxLeadingBases; leadingBases = Math.min(2*leadingBases, maxLeadingBases)) {
            final int refStart = Math.max(vc.getStart() - leadingBases, 1);

            // reference sequence starting before the variant (to give space for left-alignment) and ending at the variant end
            final byte[] refSeq = ref.getBases(new SimpleInterval(vc.getContig(), refStart, vc.getEnd()));

            final int variantOffsetInRef = vc.getStart() - refStart;

            final List<byte[]> sequences = vc.getAlleles().stream().map(a -> {
                final byte[] result = new byte[variantOffsetInRef + a.length()];
                System.arraycopy(refSeq, 0, result, 0, variantOffsetInRef);
                System.arraycopy(a.getBases(), 0, result, variantOffsetInRef, a.length());
                return result;
            }).collect(Collectors.toList());

            final List<IndexRange> alleleRanges = vc.getAlleles().stream()
                    .map(a -> new IndexRange(variantOffsetInRef + 1, variantOffsetInRef + a.length()))  // +1 to ignore the shared base in front
                    .collect(Collectors.toList());

            // note that this also shifts the index ranges as a side effect, so below they can be used to output allele bases
            final Pair<Integer, Integer> shifts = AlignmentUtils.normalizeAlleles(sequences, alleleRanges, variantOffsetInRef, trim);

            if (shifts.getLeft() == 0 && shifts.getRight() == 0) {
                return vc;
            } else if (shifts.getLeft() == variantOffsetInRef && leadingBases < maxLeadingBases) {
                continue;
            }

            final Map<Allele, Allele> alleleMap = IntStream.range(0, alleleRanges.size()).boxed()
                    .collect(Collectors.toMap(
                            n -> vc.getAlleles().get(n),
                            n -> Allele.create(Arrays.copyOfRange(sequences.get(n), variantOffsetInRef - shifts.getLeft(), variantOffsetInRef - shifts.getRight() + vc.getAlleles().get(n).length()), n == 0)));

            final GenotypesContext newGenotypes = GenotypesContext.create(vc.getNSamples());
            for (final Genotype genotype : vc.getGenotypes()) {
                final List<Allele> newAlleles = genotype.getAlleles().stream().map(a -> alleleMap.getOrDefault(a, Allele.NO_CALL)).collect(Collectors.toList());
                newGenotypes.add(new GenotypeBuilder(genotype).alleles(newAlleles).make());
            }

            return new VariantContextBuilder(vc).start(vc.getStart() - shifts.getLeft()).stop(vc.getEnd() - shifts.getRight())
                    .alleles(alleleMap.values()).genotypes(newGenotypes).make();
        }

        return vc;
    }


}

