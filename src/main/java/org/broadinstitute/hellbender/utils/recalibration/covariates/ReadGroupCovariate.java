package org.broadinstitute.hellbender.utils.recalibration.covariates;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * The Read Group covariate.
 */
public final class ReadGroupCovariate implements Covariate {
    private static final long serialVersionUID = 1L;

    public static final int MISSING_READ_GROUP_KEY = -1;

    //Note: these maps are initialized and made umodifiable at construction so the whole covariate is an immutable object once it's constructed.

    /*
     * Stores the mapping from read group id to a number.
     */
    private final Map<String, Integer> readGroupLookupTable;

    /*
     * Stores the reverse mapping, from number to read group id.
     */
    private final Map<Integer, String> readGroupReverseLookupTable;

    public ReadGroupCovariate(final List<String> readGroups){
        final Map<String, Integer> rgLookupTable = new LinkedHashMap<>();
        final Map<Integer, String> rgReverseLookupTable = new LinkedHashMap<>();

        readGroups.forEach(
                rg -> {
                    if (!rgLookupTable.containsKey(rg)) {
                        final int nextId = rgLookupTable.size();
                        rgLookupTable.put(rg, nextId); // read group index starts at 0 TODO: refactor to use a for loop
                        rgReverseLookupTable.put(nextId, rg);
                    }
                }
        );
        readGroupLookupTable = Collections.unmodifiableMap(rgLookupTable);
        readGroupReverseLookupTable = Collections.unmodifiableMap(rgReverseLookupTable);
    }

    @Override
    public void recordValues(final GATKRead read, final SAMFileHeader header, final PerReadCovariateMatrix perReadCovariateMatrix, final boolean recordIndelValues) {
        final SAMReadGroupRecord rg = ReadUtils.getSAMReadGroupRecord(read, header);
        final String readGroupIdentifier = getReadGroupIdentifier(rg); // note that the identifier by default is PU, not the ID.

        final int key = keyForReadGroup(readGroupIdentifier);
        Utils.validate(key >= -1, "key must be a nonnegative integer or the error code -1, but is " + key);

        final int readLength = read.getLength();
        for (int i = 0; i < readLength; i++) {
            perReadCovariateMatrix.addCovariate(key, key, key, i); // (mismatch, insertion, deletion) for the first three args
        }
    }

    /**
     * If present, we use the Platform Unit (PU) as the identifier of a read group, rather than the read group ID.
     * PU has the format {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE}.
     *
     * See https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
     *
     */
    public static String getReadGroupIdentifier(final SAMReadGroupRecord rg) {
        final String pu = rg.getPlatformUnit();
        return pu == null ? rg.getId() : pu;
    }

    @Override
    public String formatKey(final int key) {
        Utils.validate(readGroupReverseLookupTable.containsKey(key), () -> "missing key " + key);
        return readGroupReverseLookupTable.get(key);
    }

    @Override
    public int keyFromValue(final Object value) {
        return keyForReadGroup((String) value);
    }

    /**
     * Given the identifier (PU) for the read group, return that integer code (key) that represents it
     * in the perReadCovariateMatrix.
     *
     * @param readGroupIdentifier PU by default, read group ID if PU not available.
     * @return The integer code/key if the read group exists. -1 if the read group does not exist in the recal table.
     */
    private int keyForReadGroup(final String readGroupIdentifier) {
        if (readGroupLookupTable.containsKey(readGroupIdentifier)) {
            return readGroupLookupTable.get(readGroupIdentifier);
        } else {
            // ApplyBQSR is responsible for handling this error code appropriately; if --allow-missing-reads is set to false,
            // which is the default, it will throw an error.
            // TODO: throw an error here if this part is reached from BaseRecalibrator
            return MISSING_READ_GROUP_KEY;
        }
    }

    @Override
    public int maximumKeyValue() {
        return readGroupLookupTable.size() - 1;
    }

    public static List<String> getReadGroupIDs(final SAMFileHeader header) {
        return header.getReadGroups().stream().map(rg -> getReadGroupIdentifier(rg)).collect(Collectors.toList());
    }
}
