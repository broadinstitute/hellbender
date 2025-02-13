package org.broadinstitute.hellbender.utils.recalibration.covariates;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

/**
 * The Reported Quality Score covariate.
 */
public final class QualityScoreCovariate implements Covariate {
    private static final long serialVersionUID = 1L;

    public QualityScoreCovariate(final RecalibrationArgumentCollection RAC){
        //nothing to initialize
    }

    @Override
    public void recordValues(final GATKRead read, final SAMFileHeader header, final PerReadCovariateMatrix values, final boolean recordIndelValues) {
        final int baseQualityCount = read.getBaseQualityCount();
        final byte[] baseInsertionQualities = recordIndelValues ? ReadUtils.getBaseInsertionQualities(read) : null;
        final byte[] baseDeletionQualities = recordIndelValues ? ReadUtils.getBaseDeletionQualities(read) : null;

        //note: duplicate the loop to avoid checking recordIndelValues on every iteration
        if (recordIndelValues) {
            for (int i = 0; i < baseQualityCount; i++) {
                values.addCovariate(read.getBaseQuality(i), baseInsertionQualities[i], baseDeletionQualities[i], i);
            }
        } else {
            for (int i = 0; i < baseQualityCount; i++) {
                values.addCovariate(read.getBaseQuality(i), 0, 0, i);
            }
        }
    }

    @Override
    public String formatKey(final int key) {
        return String.format("%d", key);
    }

    @Override
    public int keyFromValue(final Object value) {
        if ( value instanceof String ){
            return Byte.parseByte((String) value);
        } else if ( value instanceof Long ) {
            return ((Long) value).intValue();
        } else {
            return (int)(Byte) value;
        }
    }

    @Override
    public int maximumKeyValue() {
        return QualityUtils.MAX_SAM_QUAL_SCORE;
    }
}
