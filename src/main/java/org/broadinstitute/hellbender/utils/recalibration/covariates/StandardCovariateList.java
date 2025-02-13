package org.broadinstitute.hellbender.utils.recalibration.covariates;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Represents the list of standard BQSR covariates.
 *
 * Note: the first two covariates ({@link ReadGroupCovariate} and {@link QualityScoreCovariate})
 * are special in the way that they are represented in the BQSR recalibration table.
 *
 * The remaining covariates are called "additional covariates".
 */
public final class StandardCovariateList implements Iterable<Covariate>, Serializable {
    private static final long serialVersionUID = 1L;
    private final ReadGroupCovariate readGroupCovariate;
    private final QualityScoreCovariate qualityScoreCovariate;
    private final List<Covariate> additionalCovariates;
    private final List<Covariate> allCovariates;

    private final Map<Class<? extends Covariate>, Integer> indexByClass;

    public static final int READ_GROUP_COVARIATE_DEFAULT_INDEX = 0;
    public static final int BASE_QUALITY_COVARIATE_DEFAULT_INDEX = 1;
    public static final int CONTEXT_COVARIATE_DEFAULT_INDEX = 2;
    public static final int CYCLE_COVARIATE_DEFAULT_INDEX = 3;
    public static final int NUM_REQUIRED_COVARITES = 2;

    /**
     * Creates a new list of standard BQSR covariates and initializes each covariate.
     */
    public StandardCovariateList(final RecalibrationArgumentCollection rac, final List<String> allReadGroups) {
        readGroupCovariate = new ReadGroupCovariate(allReadGroups);
        qualityScoreCovariate = new QualityScoreCovariate(rac);
        final ContextCovariate contextCovariate = new ContextCovariate(rac);
        final CycleCovariate cycleCovariate = new CycleCovariate(rac);

        additionalCovariates = Collections.unmodifiableList(Arrays.asList(contextCovariate, cycleCovariate));
        allCovariates = Collections.unmodifiableList(Arrays.asList(readGroupCovariate, qualityScoreCovariate, contextCovariate, cycleCovariate));
        //precompute for faster lookup (shows up on profile)
        indexByClass = new LinkedHashMap<>();
        for(int i = 0; i < allCovariates.size(); i++){
            indexByClass.put(allCovariates.get(i).getClass(), i);
        }
    }

    /**
     * Creates a new list of standard BQSR covariates and initializes each covariate.
     */
    public StandardCovariateList(final RecalibrationArgumentCollection rac, final SAMFileHeader header){
        this(rac, ReadGroupCovariate.getReadGroupIDs(header)); // tsato: but how does ApplyBQSR get the list of covariates?
    }

    /**
     * Returns 2. ReadGroupCovariate and QualityScoreCovariate are special
     */
    public int numberOfSpecialCovariates() {
        return 2;
    }

    /**
     * Returns the list of simple class names of standard covariates. The returned list is unmodifiable.
     * For example CycleCovariate.
     */
    public List<String> getStandardCovariateClassNames() {
        return Collections.unmodifiableList(allCovariates.stream().map(cov -> cov.getClass().getSimpleName()).collect(Collectors.toList()));
    }

    /**
     * Returns the size of the list of standard covariates.
     */
    public int size(){
        return allCovariates.size();
    }

    /**
     * Returns a new iterator over all covariates in this list.
     * Note: the list is unmodifiable and the iterator does not support modifying the list.
     */
    @Override
    public Iterator<Covariate> iterator() {
        return allCovariates.iterator();
    }

    public ReadGroupCovariate getReadGroupCovariate() {
        return readGroupCovariate;
    }

    public QualityScoreCovariate getQualityScoreCovariate() {
        return qualityScoreCovariate;
    }

    /**
     * returns an unmodifiable view of the additional covariates stored in this list.
     */
    public Iterable<Covariate> getAdditionalCovariates() {
        return additionalCovariates;
    }

    /**
     * Return a human-readable string representing the used covariates
     *
     * @return a non-null comma-separated string
     */
    public String covariateNames() {
        return String.join(",", getStandardCovariateClassNames());
    }

    /**
     * Get the covariate by the index.
     * @throws IndexOutOfBoundsException if the index is out of range
     *         (<tt>index &lt; 0 || index &gt;= size()</tt>)
     */
    public Covariate get(final int covIndex) {
        return allCovariates.get(covIndex);
    }

    /**
     * Returns the index of the covariate by class name or -1 if not found.
     */
    public int indexByClass(final Class<? extends Covariate> clazz){
        return indexByClass.getOrDefault(clazz, -1);
    }

    /**
     * For each covariate compute the values for all positions in this read and
     * record the values in the provided storage object.
     */
    public void populatePerReadCovariateMatrix(final GATKRead read, final SAMFileHeader header, final PerReadCovariateMatrix perReadCovariateMatrix, final boolean recordIndelValues) {
        for (int i = 0, n = allCovariates.size(); i < n; i++) {
            final Covariate cov = allCovariates.get(i);
            perReadCovariateMatrix.setCovariateIndex(i); // TODO: avoid this pattern. "cov" should already know which index it belongs to.
            cov.recordValues(read, header, perReadCovariateMatrix, recordIndelValues);
        }
    }

    /**
     * Retrieves a covariate by the parsed name {@link Covariate#parseNameForReport()} or null
     * if no covariate with that name exists in the list.
     */
    public Covariate getCovariateByParsedName(final String covName) {
        return allCovariates.stream().filter(cov -> cov.parseNameForReport().equals(covName)).findFirst().orElse(null);
    }

}
