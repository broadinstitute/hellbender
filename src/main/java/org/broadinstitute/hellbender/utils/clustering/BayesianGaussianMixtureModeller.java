package org.broadinstitute.hellbender.utils.clustering;

import org.apache.commons.lang3.NotImplementedException;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.io.Serializable;

public final class BayesianGaussianMixtureModeller implements Serializable {
    private static final long serialVersionUID = 1L;

    public enum InitMethod {
        K_MEANS_PLUS_PLUS, RANDOM, TEST
    }

    private BayesianGaussianMixtureModeller(final int nComponents,
                                            final double tol,
                                            final double regCovar,
                                            final int maxIter,
                                            final int nInit,
                                            final InitMethod initMethod,
                                            final double weightConcentrationPrior,
                                            final double meanPrecisionPrior,
                                            final RealVector meanPrior,
                                            final Double degreesOfFreedomPrior,
                                            final RealMatrix covariancePrior,
                                            final int seed,
                                            final boolean warmStart,
                                            final int verboseInterval,
                                            final double relativeSymmetryThreshold,
                                            final double absolutePositivityThreshold,
                                            final double epsilon) {
        throw new NotImplementedException("BGMM module implemented in separate PR.");
    }
}