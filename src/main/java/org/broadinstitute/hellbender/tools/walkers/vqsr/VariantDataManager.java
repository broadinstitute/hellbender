package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.AnnotationUtils;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;

import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;
import java.util.stream.Collectors;


public class VariantDataManager {
    private List<VariantDatum> data = Collections.emptyList();
    private double[] meanVector;
    private double[] varianceVector; // this is really the standard deviation
    public List<String> annotationKeys;
    private final VariantRecalibratorArgumentCollection VRAC;
    protected final static Logger logger = LogManager.getLogger(VariantDataManager.class);
    protected final List<TrainingSet> trainingSets;
    private static final double SAFETY_OFFSET = 0.01;     //To use for example as 1/(X + SAFETY_OFFSET) to protect against dividing or taking log of X=0.
    private static final double PRECISION = 0.01;         //To use mainly with MathUtils.compareDoubles(a,b,PRECISION)

    public VariantDataManager( final List<String> annotationKeys, final VariantRecalibratorArgumentCollection VRAC ) {
        this.data = Collections.emptyList();
        final List<String> uniqueAnnotations = annotationKeys.stream().distinct().collect(Collectors.toList());
        if (annotationKeys.size() != uniqueAnnotations.size()) {
            logger.warn("Ignoring duplicate annotations for recalibration {}", Utils.getDuplicatedItems(annotationKeys));
        }
        this.annotationKeys = new ArrayList<>( uniqueAnnotations );
        this.VRAC = VRAC;
        meanVector = new double[this.annotationKeys.size()];
        varianceVector = new double[this.annotationKeys.size()];
        trainingSets = new ArrayList<>();
    }

    public void setData( final List<VariantDatum> data ) {
        this.data = data;
    }

    public void setNormalization(final Map<String, Double> anMeans, final Map<String, Double> anStdDevs) {
        for (int i = 0; i < this.annotationKeys.size(); i++) {
            meanVector[i] = anMeans.get(annotationKeys.get(i));
            varianceVector[i] = anStdDevs.get(annotationKeys.get(i));
        }
    }

    public List<VariantDatum> getData() {
        return data;
    }

    /**
     * Normalize annotations to mean 0 and standard deviation 1.
     * Order the variant annotations by the provided list {@code theOrder} or standard deviation.
     *
     * @param calculateMeans Boolean indicating whether or not to calculate the means
     * @param theOrder a list of integers specifying the desired annotation order. If this is null
     *                 annotations will get sorted in decreasing size of their standard deviations.
     */
    public void normalizeData(final boolean calculateMeans, List<Integer> theOrder) {
        boolean foundZeroVarianceAnnotation = false;
        for( int iii = 0; iii < meanVector.length; iii++ ) {
            final double theMean, theSTD;
            if (calculateMeans) {
                theMean = mean(iii, true);
                theSTD = standardDeviation(theMean, iii, true);
                if (Double.isNaN(theMean)) {
                    throw new UserException.BadInput("Values for " + annotationKeys.get(iii) + " annotation not detected for ANY training variant in the input callset. VariantAnnotator may be used to add these annotations.");
                }

                foundZeroVarianceAnnotation = foundZeroVarianceAnnotation || (theSTD < 1E-5);
                meanVector[iii] = theMean;
                varianceVector[iii] = theSTD;
            }
            else {
                theMean = meanVector[iii];
                theSTD = varianceVector[iii];
            }
            logger.info(annotationKeys.get(iii) + String.format(": \t mean = %.2f\t standard deviation = %.2f", theMean, theSTD));
            for( final VariantDatum datum : data ) {
                // Transform each data point via: (x - mean) / standard deviation
                datum.annotations[iii] = ( datum.isNull[iii] ? 0.1 * Utils.getRandomGenerator().nextGaussian() : ( datum.annotations[iii] - theMean ) / theSTD );
            }
        }
        if( foundZeroVarianceAnnotation ) {
            throw new UserException.BadInput( "Found annotations with zero variance. They must be excluded before proceeding." );
        }

        // trim data by standard deviation threshold and mark failing data for exclusion later
        for( final VariantDatum datum : data ) {
            boolean remove = false;
            for( final double val : datum.annotations ) {
                remove = remove || (Math.abs(val) > VRAC.STD_THRESHOLD);
            }
            datum.failingSTDThreshold = remove;
        }

        // re-order the data by increasing standard deviation so that the results don't depend on the order things were specified on the command line
        // standard deviation over the training points is used as a simple proxy for information content, perhaps there is a better thing to use here
        // or use the serialized report's annotation order via the argument theOrder
        if (theOrder == null){
            theOrder = calculateSortOrder(meanVector);
        }
        annotationKeys = reorderList(annotationKeys, theOrder);
        varianceVector = ArrayUtils.toPrimitive(reorderArray(ArrayUtils.toObject(varianceVector), theOrder));
        meanVector = ArrayUtils.toPrimitive(reorderArray(ArrayUtils.toObject(meanVector), theOrder));
        for( final VariantDatum datum : data ) {
            datum.annotations = ArrayUtils.toPrimitive(reorderArray(ArrayUtils.toObject(datum.annotations), theOrder));
            datum.isNull = ArrayUtils.toPrimitive(reorderArray(ArrayUtils.toObject(datum.isNull), theOrder));
        }
        logger.info("Annotation order is: " + annotationKeys.toString());
    }

    public double[] getMeanVector() {
        return meanVector;
    }

    public double[] getVarianceVector() {
        return varianceVector;
    }

    /**
     * Get a list of indices which give the ascending sort order of the data array
     * @param inputVector the data to consider
     * @return a non-null list of integers with length matching the length of the input array
     */
    protected List<Integer> calculateSortOrder(final double[] inputVector) {
        final List<Integer> theOrder = new ArrayList<>(inputVector.length);
        final List<MyDoubleForSorting> toBeSorted = new ArrayList<>(inputVector.length);
        int count = 0;
        for( int iii = 0; iii < inputVector.length; iii++ ) {
            toBeSorted.add(new MyDoubleForSorting(-1.0 * Math.abs(inputVector[iii] - mean(iii, false)), count++));
        }
        Collections.sort(toBeSorted);
        for( final MyDoubleForSorting d : toBeSorted ) {
            theOrder.add(d.originalIndex); // read off the sort order by looking at the index field
        }
        return theOrder;
    }

    // small private class to assist in reading off the new ordering of the annotation array
    private class MyDoubleForSorting implements Comparable<MyDoubleForSorting> {
        final Double myData;
        final int originalIndex;

        public MyDoubleForSorting(final double myData, final int originalIndex) {
            this.myData = myData;
            this.originalIndex = originalIndex;
        }

        @Override
        public int compareTo(final MyDoubleForSorting other) {
            return myData.compareTo(other.myData);
        }
    }

    /**
     * Convenience connector method to work with arrays instead of lists. See ##reorderList##
     */
    private <T> T[] reorderArray(final T[] data, final List<Integer> order) {
        return reorderList(Arrays.asList(data), order).toArray(data);
    }

    /**
     * Reorder the given data list to be in the specified order
     * @param data the data to reorder
     * @param order the new order to use
     * @return a reordered list of data
     */
    private <T> List<T> reorderList(final List<T> data, final List<Integer> order) {
        final List<T> returnList = new ArrayList<>(data.size());
        for( final int index : order ) {
            returnList.add( data.get(index) );
        }
        return returnList;
    }

    /**
     * Convert a normalized point to it's original annotation value
     *
     * norm = (orig - mu) / sigma
     * orig = norm * sigma + mu
     *
     * @param normalizedValue the normalized value of the ith annotation
     * @param annI the index of the annotation value
     * @return the denormalized value for the annotation
     */
    public double denormalizeDatum(final double normalizedValue, final int annI) {
        final double mu = meanVector[annI];
        final double sigma = varianceVector[annI];
        return normalizedValue * sigma + mu;
    }

    public void addTrainingSet( final TrainingSet trainingSet ) {
        trainingSets.add( trainingSet );
    }

    public List<String> getAnnotationKeys() {
        return annotationKeys;
    }

    public boolean checkHasTrainingSet() {
        for( final TrainingSet trainingSet : trainingSets ) {
            if( trainingSet.isTraining ) { return true; }
        }
        return false;
    }

    public boolean checkHasTruthSet() {
        for( final TrainingSet trainingSet : trainingSets ) {
            if( trainingSet.isTruth ) { return true; }
        }
        return false;
    }

    public List<VariantDatum> getTrainingData() {
        final List<VariantDatum> trainingData = new ArrayList<>();
        for( final VariantDatum datum : data ) {
            if( datum.atTrainingSite && !datum.failingSTDThreshold ) {
                trainingData.add( datum );
            } else if (datum.failingSTDThreshold && VRAC.debugStdevThresholding) {
                logger.warn("Datum at " + datum.loc + " with ref " + datum.referenceAllele + " and alt " + datum.alternateAllele + " failing std thresholding: " + Arrays.toString(datum.annotations));
            }
        }
        logger.info( "Training with " + trainingData.size() + " variants after standard deviation thresholding." );
        if( trainingData.size() < VRAC.MIN_NUM_BAD_VARIANTS ) {
            logger.warn( "WARNING: Training with very few variant sites! Please check the model reporting PDF to ensure the quality of the model is reliable." );
        } else if( trainingData.size() > VRAC.MAX_NUM_TRAINING_DATA ) {
            logger.warn( "WARNING: Very large training set detected. Downsampling to " + VRAC.MAX_NUM_TRAINING_DATA + " training variants." );
            Collections.shuffle(trainingData, Utils.getRandomGenerator());
            return trainingData.subList(0, VRAC.MAX_NUM_TRAINING_DATA);
        }
        return trainingData;
    }

    public List<VariantDatum> selectWorstVariants() {
        final List<VariantDatum> trainingData = new ArrayList<>();

        for( final VariantDatum datum : data ) {
            if( datum != null && !datum.failingSTDThreshold && !Double.isInfinite(datum.lod) && datum.lod < VRAC.BAD_LOD_CUTOFF ) {
                datum.atAntiTrainingSite = true;
                trainingData.add( datum );
            }
        }

        logger.info( "Selected worst " + trainingData.size() + " scoring variants --> variants with LOD <= " + String.format("%.4f", VRAC.BAD_LOD_CUTOFF) + "." );

        return trainingData;
    }

    public List<VariantDatum> getEvaluationData() {
        final List<VariantDatum> evaluationData = new ArrayList<>();

        for( final VariantDatum datum : data ) {
            if( datum != null && !datum.failingSTDThreshold && !datum.atTrainingSite && !datum.atAntiTrainingSite ) {
                evaluationData.add( datum );
            }
        }

        return evaluationData;
    }

    /**
     * Remove all VariantDatum's from the data list which are marked as aggregate data
     */
    public void dropAggregateData() {
        final Iterator<VariantDatum> iter = data.iterator();
        while (iter.hasNext()) {
            final VariantDatum datum = iter.next();
            if( datum.isAggregate ) {
                iter.remove();
            }
        }
    }

    public List<VariantDatum> getRandomDataForPlotting( final int numToAdd, final List<VariantDatum> trainingData, final List<VariantDatum> antiTrainingData, final List<VariantDatum> evaluationData ) {
        final List<VariantDatum> returnData = new ArrayList<>();
        Collections.shuffle(trainingData, Utils.getRandomGenerator());
        Collections.shuffle(antiTrainingData, Utils.getRandomGenerator());
        Collections.shuffle(evaluationData, Utils.getRandomGenerator());
        returnData.addAll(trainingData.subList(0, Math.min(numToAdd, trainingData.size())));
        returnData.addAll(antiTrainingData.subList(0, Math.min(numToAdd, antiTrainingData.size())));
        returnData.addAll(evaluationData.subList(0, Math.min(numToAdd, evaluationData.size())));
        Collections.shuffle(returnData, Utils.getRandomGenerator());
        return returnData;
    }

    protected double mean( final int index, final boolean trainingData ) {
        double sum = 0.0;
        int numNonNull = 0;
        for( final VariantDatum datum : data ) {
            if( (trainingData == datum.atTrainingSite) && !datum.isNull[index] ) {
                sum += datum.annotations[index];
                numNonNull++;
            }
        }
        return sum / ((double) numNonNull);
    }

    protected double standardDeviation( final double mean, final int index, final boolean trainingData ) {
        double sum = 0.0;
        int numNonNull = 0;
        for( final VariantDatum datum : data ) {
            if( (trainingData == datum.atTrainingSite) && !datum.isNull[index] ) { sum += ((datum.annotations[index] - mean)*(datum.annotations[index] - mean)); numNonNull++; }
        }
        return Math.sqrt( sum / ((double) numNonNull) );
    }

    public void decodeAnnotations( final VariantDatum datum, final VariantContext vc, final boolean jitter ) {
        final double[] annotations = new double[annotationKeys.size()];
        final boolean[] isNull = new boolean[annotationKeys.size()];
        int iii = 0;
        for( final String key : annotationKeys ) {
            isNull[iii] = false;
            annotations[iii] = decodeAnnotation( key, vc, jitter, VRAC, datum );
            if( Double.isNaN(annotations[iii]) ) { isNull[iii] = true; }
            iii++;
        }
        datum.annotations = annotations;
        datum.isNull = isNull;
    }
    /** Transforms an interval [xmin, xmax] to (-inf, +inf) **/
    private static double logitTransform( final double x, final double xmin, final double xmax) {
        return Math.log((x - xmin)/(xmax - x));
    }

    private static double decodeAnnotation( final String annotationKey, final VariantContext vc, final boolean jitter, final VariantRecalibratorArgumentCollection vrac, final VariantDatum datum ) {
        double value;

        final double LOG_OF_TWO = 0.6931472;

        try {
            //if we're in allele-specific mode and an allele-specific annotation has been requested, parse the appropriate value from the list
            if(vrac.useASannotations && annotationKey.startsWith(GATKVCFConstants.ALLELE_SPECIFIC_PREFIX)) {
                final List<Object> valueList = vc.getAttributeAsList(annotationKey);
                //FIXME: we need to look at the ref allele here too
                if (vc.hasAllele(datum.alternateAllele)) {
                    final int altIndex = vc.getAlleleIndex(datum.alternateAllele)-1; //-1 is to convert the index from all alleles (including reference) to just alternate alleles
                    value = Double.parseDouble((String)valueList.get(altIndex));
                }
                //if somehow our alleles got mixed up
                else
                    throw new IllegalStateException("VariantDatum allele " + datum.alternateAllele + " is not contained in the input VariantContext.");
            }
            else
                value = vc.getAttributeAsDouble( annotationKey, Double.NaN );
            if( Double.isInfinite(value) ) { value = Double.NaN; }
            if( jitter && annotationKey.equalsIgnoreCase(GATKVCFConstants.HAPLOTYPE_SCORE_KEY) && MathUtils.compareDoubles(value, 0.0, PRECISION) == 0 ) { value += 0.01 * Utils.getRandomGenerator().nextGaussian(); }
            if( jitter && (annotationKey.equalsIgnoreCase(GATKVCFConstants.FISHER_STRAND_KEY) || annotationKey.equalsIgnoreCase(GATKVCFConstants.AS_FILTER_STATUS_KEY)) && MathUtils.compareDoubles(value, 0.0, PRECISION) == 0 ) { value += 0.01 * Utils.getRandomGenerator().nextGaussian(); }
            if( jitter && annotationKey.equalsIgnoreCase(GATKVCFConstants.INBREEDING_COEFFICIENT_KEY) && MathUtils.compareDoubles(value, 0.0, PRECISION) == 0 ) { value += 0.01 * Utils.getRandomGenerator().nextGaussian(); }
            if( jitter && (annotationKey.equalsIgnoreCase(GATKVCFConstants.STRAND_ODDS_RATIO_KEY) || annotationKey.equalsIgnoreCase(GATKVCFConstants.AS_STRAND_ODDS_RATIO_KEY)) && MathUtils.compareDoubles(value, LOG_OF_TWO, PRECISION) == 0 ) { value += 0.01 * Utils.getRandomGenerator().nextGaussian(); }   //min SOR is 2.0, then we take ln
            if( jitter && (annotationKey.equalsIgnoreCase(VCFConstants.RMS_MAPPING_QUALITY_KEY))) {
                if( vrac.MQ_CAP > 0) {
                    value = logitTransform(value, -SAFETY_OFFSET, vrac.MQ_CAP + SAFETY_OFFSET);
                    if (MathUtils.compareDoubles(value, logitTransform(vrac.MQ_CAP, -SAFETY_OFFSET, vrac.MQ_CAP + SAFETY_OFFSET), PRECISION) == 0 ) {
                        value += vrac.MQ_JITTER * Utils.getRandomGenerator().nextGaussian();
                    }
                } else if( MathUtils.compareDoubles(value, vrac.MQ_CAP, PRECISION) == 0 ) {
                    value += vrac.MQ_JITTER * Utils.getRandomGenerator().nextGaussian();
                }
            }
            if( jitter && (annotationKey.equalsIgnoreCase(GATKVCFConstants.AS_RMS_MAPPING_QUALITY_KEY))){
                value += vrac.MQ_JITTER * Utils.getRandomGenerator().nextGaussian();
            }
        } catch( NumberFormatException e ) {
            value = Double.NaN; // VQSR works with missing data by marginalizing over the missing dimension when evaluating the Gaussian mixture model
        }

        return value;
    }

    public void parseTrainingSets(
            final FeatureContext featureContext,
            final VariantContext evalVC,
            final VariantDatum datum,
            final boolean TRUST_ALL_POLYMORPHIC ) {
        datum.isKnown = false;
        datum.atTruthSite = false;
        datum.atTrainingSite = false;
        datum.atAntiTrainingSite = false;
        datum.prior = 2.0;

        for( final TrainingSet trainingSet : trainingSets ) {
            List<VariantContext> vcs = featureContext.getValues(trainingSet.variantSource, featureContext.getInterval().getStart());
            for( final VariantContext trainVC : vcs ) {
                if (VRAC.useASannotations && !doAllelesMatch(trainVC, datum))
                    continue;
                if( isValidVariant( evalVC, trainVC, TRUST_ALL_POLYMORPHIC ) ) {
                    datum.isKnown = datum.isKnown || trainingSet.isKnown;
                    datum.atTruthSite = datum.atTruthSite || trainingSet.isTruth;
                    datum.atTrainingSite = datum.atTrainingSite || trainingSet.isTraining;
                    datum.prior = Math.max( datum.prior, trainingSet.prior );
                }
                if( trainVC != null ) {
                    datum.atAntiTrainingSite = datum.atAntiTrainingSite || trainingSet.isAntiTraining;
                }
            }
        }
    }

    private boolean isValidVariant( final VariantContext evalVC, final VariantContext trainVC, final boolean TRUST_ALL_POLYMORPHIC) {
        return trainVC != null && trainVC.isNotFiltered() && trainVC.isVariant() && checkVariationClass( evalVC, trainVC ) &&
                (TRUST_ALL_POLYMORPHIC || !trainVC.hasGenotypes() || trainVC.isPolymorphicInSamples());
    }

    private boolean doAllelesMatch(final VariantContext trainVC, final VariantDatum datum) {
        //only do this check in the allele-specific case, where each datum represents one allele
        if (datum.alternateAllele == null) {
            return true;
        }
        try {
            return GATKVariantContextUtils.isAlleleInList(datum.referenceAllele, datum.alternateAllele, trainVC.getReference(), trainVC.getAlternateAlleles());
        } catch (final IllegalStateException e) {
            throw new IllegalStateException("Reference allele mismatch at position " + trainVC.getContig() + ":" + trainVC.getStart() + " : ", e);
        }
    }

    protected static boolean checkVariationClass( final VariantContext evalVC, final VariantContext trainVC ) {
        switch( trainVC.getType() ) {
            case SNP:
            case MNP:
                return checkVariationClass( evalVC, VariantRecalibratorArgumentCollection.Mode.SNP );
            case INDEL:
            case MIXED:
            case SYMBOLIC:
                return checkVariationClass( evalVC, VariantRecalibratorArgumentCollection.Mode.INDEL );
            default:
                return false;
        }
    }

    protected static boolean checkVariationClass( final VariantContext evalVC, final VariantRecalibratorArgumentCollection.Mode mode ) {
        switch( mode ) {
            case SNP:
                return evalVC.isSNP() || evalVC.isMNP();
            case INDEL:
                return evalVC.isStructuralIndel() || evalVC.isIndel() || evalVC.isMixed() || evalVC.isSymbolic();
            case BOTH:
                return true;
            default:
                throw new IllegalStateException( "Encountered unknown recal mode: " + mode );
        }
    }

    protected static boolean checkVariationClass( final VariantContext evalVC, final Allele allele, final VariantRecalibratorArgumentCollection.Mode mode ) {
        switch( mode ) {
            case SNP:
                //note that spanning deletions are considered SNPs by this logic
                return evalVC.getReference().length() == allele.length();
            case INDEL:
                return (evalVC.getReference().length() != allele.length()) || allele.isSymbolic();
            case BOTH:
                return true;
            default:
                throw new IllegalStateException( "Encountered unknown recal mode: " + mode );
        }
    }

    public void writeOutRecalibrationTable(final VariantContextWriter recalWriter, final SAMSequenceDictionary seqDictionary) {
        // we need to sort in coordinate order in order to produce a valid VCF
        Collections.sort( data, VariantDatum.getComparator(seqDictionary) );

        // create dummy alleles to be used
        List<Allele> alleles = Arrays.asList(Allele.create("N", true), Allele.create("<VQSR>", false));

        for( final VariantDatum datum : data ) {
            if (VRAC.useASannotations)
                alleles = Arrays.asList(datum.referenceAllele, datum.alternateAllele); //use the alleles to distinguish between multiallelics in AS mode
            VariantContextBuilder builder = new VariantContextBuilder("VQSR", datum.loc.getContig(), datum.loc.getStart(), datum.loc.getEnd(), alleles);
            builder.attribute(VCFConstants.END_KEY, datum.loc.getEnd());
            builder.attribute(GATKVCFConstants.VQS_LOD_KEY, String.format("%.4f", datum.lod));
            builder.attribute(GATKVCFConstants.CULPRIT_KEY, (datum.worstAnnotation != -1 ? annotationKeys.get(datum.worstAnnotation) : "NULL"));

            if ( datum.atTrainingSite ) builder.attribute(GATKVCFConstants.POSITIVE_LABEL_KEY, true);
            if ( datum.atAntiTrainingSite ) builder.attribute(GATKVCFConstants.NEGATIVE_LABEL_KEY, true);

            recalWriter.add(builder.make());
        }
    }
}
