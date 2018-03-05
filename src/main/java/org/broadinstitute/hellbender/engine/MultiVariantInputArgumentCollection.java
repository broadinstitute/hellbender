package org.broadinstitute.hellbender.engine;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * Class that defines the variant arguments used for a MultiVariantWalker.  DefaultArgumentCollection uses the standard --variant argument; however,
 * subclasses of MultiVariantWalker can override MultiVariantWalker.getMultiVariantInputArgumentCollection() and provide their own argument pattern.
 */
public abstract class MultiVariantInputArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    abstract protected List<String> getDrivingVariantPaths();

    public static class DefaultArgumentCollection extends MultiVariantInputArgumentCollection {
        private static final long serialVersionUID = 1L;

        // NOTE: using List<String> rather than List<FeatureInput> here so that we can initialize the driving source of variants separately
        // from any other potential sources of Features.
        @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME,
                doc = "One or more VCF files containing variants", common = false, optional = false)
        public List<String> drivingVariantFiles = new ArrayList<>();

        @Override
        protected List<String> getDrivingVariantPaths() {
            return drivingVariantFiles;
        }
    }
}
