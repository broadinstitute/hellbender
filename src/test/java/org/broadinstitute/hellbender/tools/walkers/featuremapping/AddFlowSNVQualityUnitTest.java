package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import org.broadinstitute.hellbender.testutils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class AddFlowSNVQualityUnitTest extends BaseTest {

    @DataProvider(name = "sliceIsValidDataProvider")
    Object[][] sliceIsValidDataProvider() {
        return new Object[][] {
                { new int[] {0, 0, 0, 0, 0, 0}, 4, false },
                { new int[] {0, 0, 1, 0, 0, 1}, 4, true },
                { new int[] {1, 0, 0, 0, 1, 1}, 4, false },
        };
    }

    @Test(dataProvider = "sliceIsValidDataProvider")
    void testSliceIsValid(final int[] slice, final int flowOrderLength, final boolean isValid) {
        Assert.assertEquals(AddFlowSNVQuality.sliceIsValidForConsideration(slice, flowOrderLength), isValid);
    }

    @DataProvider(name = "getSnvqDataProvider")
    Object[][] getSnvqDataProvider() {
        return new Object[][] {
                { 0.5, 0.2, 0.3, AddFlowSNVQualityArgumentCollection.SnvqModeEnum.Legacy, 0.5 },
                { 0.5, 0.2, 0.3, AddFlowSNVQualityArgumentCollection.SnvqModeEnum.Optimistic, 0.06 },
                { 0.5, 0.2, 0.3, AddFlowSNVQualityArgumentCollection.SnvqModeEnum.Pessimistic, 0.44 },
                { 0.5, 0.2, 0.3, AddFlowSNVQualityArgumentCollection.SnvqModeEnum.Geometric, Math.sqrt(0.06 * 0.44) },
        };
    }

    @Test(dataProvider = "getSnvqDataProvider")
    void testGetSnvq(final double sliceP, final double p1, final double p2, AddFlowSNVQualityArgumentCollection.SnvqModeEnum snvMode, final double expected) {
        Assert.assertEquals(AddFlowSNVQuality.getSnvq(sliceP, p1, p2, snvMode), expected, 0.0001);
    }

}
