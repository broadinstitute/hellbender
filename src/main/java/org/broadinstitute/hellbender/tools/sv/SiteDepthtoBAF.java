package org.broadinstitute.hellbender.tools.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.sv.CollectSVEvidence;
import org.broadinstitute.hellbender.utils.MathUtils.RunningAverage;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.codecs.FeatureOutputCodec;
import org.broadinstitute.hellbender.utils.codecs.FeatureOutputCodecFinder;
import org.broadinstitute.hellbender.utils.codecs.FeatureSink;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

/**
 * <p>Merges locus-sorted SiteDepth evidence files, and calculates the bi-allelic frequency (baf)
 * for each sample and site, and writes these values as a BafEvidence output file.</p>
 * <p>Samples at sites that have too few depth counts (controlled by --min-total-depth), and samples
 * failing a Pearson's chi square test for goodness of fit to a bi-allelic model are excluded. (I.e.,
 * only samples that are judged to be heterozygous at a site are included.)</p>
 * <p>Across all evaluable samples, sites where the samples have bafs with too much variance
 * (controlled by --max-std) are excluded.  The median baf value at each site is adjusted to be
 * exactly 0.5.</p>
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         One or more SiteDepth evidence files, or a file containing a list of evidence files, one per line.
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         An output file containing merged BafEvidence from the inputs.
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk SiteDepthtoBAF \
 *       -F file1.sd.txt.gz [-F file2.sd.txt.gz ...] \
 *       -O merged.baf.bci
 * </pre>
 *
 * @author Ted Sharpe &lt;tsharpe@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Convert SiteDepth to BafEvidence",
        oneLineSummary = "Convert SiteDepth to BafEvidence",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@ExperimentalFeature
@DocumentedFeature
public class SiteDepthtoBAF extends MergingMultiFeatureWalker<SiteDepth> {
    public static final String LOCUS_DEPTH_FILE_NAME = "site-depth";
    public static final String BAF_SITES_VCF_LONG_NAME = "baf-sites-vcf";
    public static final String SAMPLE_NAMES_NAME = "sample-names";
    public static final String BAF_EVIDENCE_FILE_NAME = "baf-evidence-output";
    public static final String MIN_TOTAL_DEPTH_NAME = "min-total-depth";
    public static final String MAX_STD_NAME = "max-std";
    public static final String MIN_HET_PROBABILITY = "min-het-probability";
    public static final String COMPRESSION_LEVEL_NAME = "compression-level";
    public static ChiSquaredDistribution chiSqDist = new ChiSquaredDistribution(null, 1.);

    @Argument(
            doc = "SiteDepth files to process",
            fullName = LOCUS_DEPTH_FILE_NAME,
            shortName = StandardArgumentDefinitions.FEATURE_SHORT_NAME )
    private List<FeatureInput<SiteDepth>> siteDepthFiles;

    @Argument(fullName = BAF_SITES_VCF_LONG_NAME,
            doc = "Input VCF of SNPs marking loci for allele counts")
    public GATKPath bafSitesFile;

    @Argument(fullName = SAMPLE_NAMES_NAME,
            doc = "Sample names. Necessary only when writing a *.baf.bci output file.", optional = true)
    private List<String> sampleNames;

    @Argument(
            doc = "BAF output file",
            fullName = BAF_EVIDENCE_FILE_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME )
    private GATKPath bafOutputFile;

    @Argument(
            doc = "Maximum standard deviation across bafs at a locus (otherwise locus will be ignored)",
            fullName = MAX_STD_NAME, optional = true )
    private double maxStdDev = .2;

    @Argument(
            doc = "Minimum total depth at a locus (otherwise locus will be ignored)",
            fullName = MIN_TOTAL_DEPTH_NAME, optional = true )
    private int minTotalDepth = 10;

    @Argument(
            doc = "Minimum probability of site being biallelic and heterozygous",
            fullName = MIN_HET_PROBABILITY, optional = true )
    private double minHetProbability = .5;

    @Argument(
            doc = "Output compression level",
            fullName = COMPRESSION_LEVEL_NAME,
            minValue = 0, maxValue = 9, optional = true )
    private int compressionLevel = 4;

    private FeatureDataSource<VariantContext> bafSitesSource;
    private Iterator<VariantContext> bafSitesIterator;
    private FeatureSink<BafEvidence> writer;
    private final List<SiteDepth> sameLocusBuffer = new ArrayList<>();

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        bafSitesSource = new FeatureDataSource<>(bafSitesFile.getRawInputString());
        bafSitesIterator = new CollectSVEvidence.BAFSiteIterator(bafSitesSource.iterator());
        final SAMSequenceDictionary dict = getDictionary();
        dict.assertSameDictionary(bafSitesSource.getSequenceDictionary());
        final FeatureOutputCodec<BafEvidence,FeatureSink<BafEvidence>> codec =
                FeatureOutputCodecFinder.find(bafOutputFile, BafEvidence.class);
        writer = codec.makeSink(bafOutputFile, dict, sampleNames, compressionLevel);
    }

    @Override
    public void apply( final SiteDepth siteDepth, final Object header,
                       final ReadsContext reads, final ReferenceContext ref ) {
        if ( !sameLocus(siteDepth) ) {
            processBuffer();
        }
        sameLocusBuffer.add(siteDepth);
    }

    @Override
    public Object onTraversalSuccess() {
        super.onTraversalSuccess();
        if ( !sameLocusBuffer.isEmpty() ) {
            processBuffer();
        }
        bafSitesSource.close();
        writer.close();
        return null;
    }

    /** Performs a Pearson's chi square test for goodness of fit to a biallelic model to determine
     * heterozygosity.  If this fails, null is returned.  If the test succeeds,
     * the depth of the alt call as a fraction of the total depth is returned as BafEvidence.
     */
    @VisibleForTesting BafEvidence calcBAF( final SiteDepth sd,
                                            final int refIndex,
                                            final int altIndex ) {
        final int totalDepth = sd.getTotalDepth();
        if ( totalDepth < minTotalDepth ) {
            return null;
        }
        final double expectRefAlt = totalDepth / 2.;
        final double altDepth = sd.getDepth(altIndex);
        final double refDiff = sd.getDepth(refIndex) - expectRefAlt;
        final double altDiff = altDepth - expectRefAlt;
        final double chiSq = (refDiff * refDiff + altDiff * altDiff) / expectRefAlt;
        final double fitProb = 1. - chiSqDist.cumulativeProbability(chiSq);
        if ( fitProb < minHetProbability ) {
            return null;
        }
        return new BafEvidence(sd.getSample(), sd.getContig(), sd.getStart(), altDepth / totalDepth);
    }

    private boolean sameLocus( final Locatable locus ) {
        if ( sameLocusBuffer.isEmpty() ) {
            return true;
        }
        final SiteDepth curLocus = sameLocusBuffer.get(0);
        return curLocus.getStart() == locus.getStart() &&
                curLocus.getContig().equals(locus.getContig());
    }

    private void processBuffer() {
        final VariantContext vc = nextSite();
        if ( !sameLocus(vc) ) {
            final Locatable loc = sameLocusBuffer.get(0);
            throw new UserException("expecting locus " + loc.getContig() + ":" + loc.getStart() +
                    ", but found locus " + vc.getContig() + ":" + vc.getStart() + " in baf sites vcf");
        }
        final List<Allele> alleles = vc.getAlleles();
        final int refIndex = Nucleotide.decode(alleles.get(0).getBases()[0]).ordinal();
        if ( refIndex > 3 ) {
            throw new UserException("ref call is not [ACGT] in vcf at " + vc.getContig() + ":" + vc.getStart());
        }
        final int altIndex = Nucleotide.decode(alleles.get(1).getBases()[0]).ordinal();
        if ( altIndex > 3 ) {
            throw new UserException("alt call is not [ACGT] in vcf at " + vc.getContig() + ":" + vc.getStart());
        }
        final List<BafEvidence> beList = new ArrayList<>(sameLocusBuffer.size());
        for ( final SiteDepth sd : sameLocusBuffer ) {
            final BafEvidence bafEvidence = calcBAF(sd, refIndex, altIndex);
            if ( bafEvidence != null ) {
                beList.add(bafEvidence);
            }
        }
        sameLocusBuffer.clear();
        final int nBafs = beList.size();
        if ( nBafs <= 1 ) {
            if ( nBafs == 1 ) {
                writer.write(new BafEvidence(beList.get(0), .5));
            }
            return;
        }
        final RunningAverage ra = new RunningAverage();
        for ( final BafEvidence bafEvidence : beList ) {
            ra.add(bafEvidence.getValue());
        }
        final double stddev = ra.stddev();
        if ( stddev <= maxStdDev ) {
            final double[] vals = new double[nBafs];
            for ( int idx = 0; idx != nBafs; ++idx ) {
                vals[idx] = beList.get(idx).getValue();
            }
            Arrays.sort(vals);
            final int midIdx = nBafs / 2;
            final double median =
                    (nBafs % 2) != 0 ? vals[midIdx] : (vals[midIdx] + vals[midIdx - 1])/2.;
            final double adj = .5 - median;
            for ( final BafEvidence bafEvidence : beList ) {
                writer.write(new BafEvidence(bafEvidence, bafEvidence.getValue()+adj));
            }
        }
    }

    private VariantContext nextSite() {
        if ( !bafSitesIterator.hasNext() ) {
            throw new UserException("baf sites vcf exhausted before site depth data");
        }
        return bafSitesIterator.next();
    }
}
