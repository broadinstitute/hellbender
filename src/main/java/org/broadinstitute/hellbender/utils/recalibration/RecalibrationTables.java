package org.broadinstitute.hellbender.utils.recalibration;

import org.broadinstitute.hellbender.utils.recalibration.covariates.Covariate;
import org.broadinstitute.hellbender.utils.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.collections.NestedIntegerArray;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Utility class to facilitate base quality score recalibration.
 */
public final class RecalibrationTables implements Serializable, Iterable<NestedIntegerArray<RecalDatum>> {
    private static final long serialVersionUID = 1L;
    private final int qualDimension;
    private final int eventDimension = EventType.values().length;
    private final int numReadGroups;

    final StandardCovariateList covariates;  // save the covariates this was created with

    //These two tables are special
    private final NestedIntegerArray<RecalDatum> readGroupTable;
    private final NestedIntegerArray<RecalDatum> qualityScoreTable;

    private final List<NestedIntegerArray<RecalDatum>> additionalTables;
    private final List<NestedIntegerArray<RecalDatum>> allTables;    //special + additional

    //bidirectional mappings
    private final Map<Covariate, NestedIntegerArray<RecalDatum>> covariateToTable;
    private final Map<NestedIntegerArray<RecalDatum>, Covariate> tableToCovariate;



    public RecalibrationTables(final StandardCovariateList covariates) {
        this(covariates, covariates.getReadGroupCovariate().maximumKeyValue() + 1);
    }

    public RecalibrationTables(StandardCovariateList covariates, final int numReadGroups) {
        this.covariates = covariates;
        this.additionalTables = new ArrayList<>();
        this.allTables = new ArrayList<>();
        this.covariateToTable = new LinkedHashMap<>();
        this.tableToCovariate = new LinkedHashMap<>();

        this.qualDimension = covariates.getQualityScoreCovariate().maximumKeyValue() + 1;
        this.numReadGroups = numReadGroups;

        // two special (i.e. required) tables: read group table and reported quality score table
        this.readGroupTable = new NestedIntegerArray<>(numReadGroups, eventDimension);
        allTables.add(readGroupTable);
        covariateToTable.put(covariates.getReadGroupCovariate(), readGroupTable);
        tableToCovariate.put(readGroupTable, covariates.getReadGroupCovariate());

        this.qualityScoreTable = makeQualityScoreTable();
        allTables.add(qualityScoreTable);
        covariateToTable.put(covariates.getQualityScoreCovariate(), qualityScoreTable);
        tableToCovariate.put(qualityScoreTable, covariates.getQualityScoreCovariate());

        //Non-special tables
        for (Covariate cov : covariates.getAdditionalCovariates()){
            final NestedIntegerArray<RecalDatum> table = new NestedIntegerArray<>(numReadGroups, qualDimension, cov.maximumKeyValue() + 1, eventDimension);
            additionalTables.add(table);
            allTables.add(table);
            covariateToTable.put(cov, table);
            tableToCovariate.put(table, cov);
        }
    }

    public NestedIntegerArray<RecalDatum> getTableForCovariate(Covariate cov) {
        return covariateToTable.get(cov);
    }

    public Covariate getCovariateForTable(NestedIntegerArray<RecalDatum> table) {
        return tableToCovariate.get(table);
    }

    public boolean isReadGroupTable(NestedIntegerArray<RecalDatum> table) {
        //Note: NestedIntegerArray does not implement equals so we use reference identity to check equality
        //to explicitly check identity (and future-proof against an equals method in NestedIntegerArray).
        return table == getReadGroupTable();
    }

    public boolean isQualityScoreTable(NestedIntegerArray<RecalDatum> table) {
        //Note: NestedIntegerArray does not implement equals so we use reference identity to check equality
        //to explicitly check identity (and future-proof against an equals method in NestedIntegerArray).
        return table == getQualityScoreTable();
    }

    public boolean isAdditionalCovariateTable(NestedIntegerArray<RecalDatum> table) {
        return additionalTables.contains(table);
    }

    public NestedIntegerArray<RecalDatum> getReadGroupTable() {
        return readGroupTable;
    }

    public NestedIntegerArray<RecalDatum> getQualityScoreTable() {
        return qualityScoreTable;
    }

    public int numTables() {
        return allTables.size();
    }

    @Override
    public Iterator<NestedIntegerArray<RecalDatum>> iterator() {
        return allTables.iterator();
    }

    /**
     * @return true if all the tables contain no RecalDatums
     */
    public boolean isEmpty() {
        for( final NestedIntegerArray<RecalDatum> table : allTables ) {
            if( !table.getAllValues().isEmpty() ) { return false; }
        }
        return true;
    }

    /**
     * Allocate a new quality score table, based on requested parameters
     * in this set of tables, without any data in it.
     * @return a newly allocated, empty read group x quality score table
     */
    public NestedIntegerArray<RecalDatum> makeQualityScoreTable() {
        return new NestedIntegerArray<>(numReadGroups, qualDimension, eventDimension);
    }

    /**
     * Merge all of the tables from toMerge into this set of tables
     */
    public RecalibrationTables combine(final RecalibrationTables toMerge) {
        if ( numTables() != toMerge.numTables() )
            throw new IllegalArgumentException("Attempting to merge RecalibrationTables with different sizes");

        for ( int i = 0; i < numTables(); i++ ) {
            final NestedIntegerArray<RecalDatum> myTable = this.allTables.get(i);
            final NestedIntegerArray<RecalDatum> otherTable = toMerge.allTables.get(i);
            RecalUtils.combineTables(myTable, otherTable);
        }

        return this;
    }

    /**
     * Combines the two tables into a new table (allocating a new table in the process)
     *
     * @param left first table to combine
     * @param right second table to combine
     * @return a new table with the merged contents of left and right
     */
    public static RecalibrationTables safeCombine(final RecalibrationTables left, final RecalibrationTables right){
        Utils.nonNull(left);
        Utils.nonNull(right);

        final RecalibrationTables newTable = new RecalibrationTables(left.covariates, left.numReadGroups);
        newTable.combine(left);
        newTable.combine(right);
        return newTable;
    }

    /**
     * Combines the right table into the left table, in-place (without making a copy)
     *
     * @param left first table to combine
     * @param right second table to combine
     * @return modified version of left with the contents of right incorporated into it
     */
    public static RecalibrationTables inPlaceCombine(final RecalibrationTables left, final RecalibrationTables right){
        Utils.nonNull(left);
        Utils.nonNull(right);

        return left.combine(right);
    }

    //XXX this should not be accessible by index
    public NestedIntegerArray<RecalDatum> getTable(int index) {
        return allTables.get(index);
    }

    public List<NestedIntegerArray<RecalDatum>> getAdditionalTables() {
        return additionalTables;
    }
}
