package org.broadinstitute.hellbender.utils.param;

import org.apache.commons.lang3.Range;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

/**
 * Created by lichtens on 8/25/15.
 */
public class ParamUtilsUnitTest {

    @Test
    public void testInRangeSuccess(){
        Assert.assertTrue(4 == ParamUtils.inRange(4L, 3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(4 == ParamUtils.inRange(4, 3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(4 == ParamUtils.inRange(Range.<Integer>of(3, 6), 4, "error"));
        Assert.assertTrue(4.1 == ParamUtils.inRange(4.1, 3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(4.1 == ParamUtils.inRange(4.1, -3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(4.1 == ParamUtils.inRange(Range.<Double>of(-3.0, 6.0), 4.1, "error"));
        Assert.assertTrue(0 == ParamUtils.inRange(0L, -3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(0 == ParamUtils.inRange(0, -3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(0 == ParamUtils.inRange(Range.<Integer>of(-3, 6), 0, "error"));
        Assert.assertTrue(0.0 == ParamUtils.inRange(0.0, -3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(0.0 == ParamUtils.inRange(Range.<Double>of(-3.0, 6.0), 0.0, "error"));
        Assert.assertTrue(0 == ParamUtils.inRange(0L, -3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(0 == ParamUtils.inRange(0, -3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(0 == ParamUtils.inRange(Range.<Integer>of(-3, 6), 0, "error"));
        Assert.assertTrue(-1 == ParamUtils.inRange(-1L, -3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(-1 == ParamUtils.inRange(-1, -3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(-1 == ParamUtils.inRange(Range.<Integer>of(-3, 6), -1, "error"));
        Assert.assertTrue(-1.5 == ParamUtils.inRange(-1.5, -3, 6, "Range calculation did not work properly"), "Did not return proper value");
        Assert.assertTrue(-1.5 == ParamUtils.inRange(Range.<Double>of(-3.0, 6.0), -1.5, "error"));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureAllPositiveLong(){
        ParamUtils.inRange(4L, 7, 10, "Range calculation did not work properly");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureAllPositiveInt(){
        ParamUtils.inRange(4, 7, 10, "Range calculation did not work properly");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureAllPositiveIntRange(){
        ParamUtils.inRange(Range.<Integer>of(7, 10), 4, "Range calculation did not work properly");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureAllPositiveDoubleRange(){
        ParamUtils.inRange(Range.<Double>of(7.0, 10.0), 4, "Range calculation did not work properly");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureNaN(){
        ParamUtils.inRange(Double.NaN, 7, 10, "Range calculation did not work properly");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureNaNHigh(){
        ParamUtils.inRange(7, 10, Double.NaN, "Range calculation did not work properly");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureNaNDoubleRange(){
        ParamUtils.inRange(Range.<Double>of(7.0, 10.0), Double.NaN, "Range calculation did not work properly");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureNaNLow(){
        ParamUtils.inRange(7, Double.NaN, 10, "Range calculation did not work properly");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureNaNTwo1(){
        ParamUtils.inRange(7, Double.NaN, Double.NaN, "Range calculation did not work properly");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureNaNTwo2(){
        ParamUtils.inRange(Double.NaN, 10, Double.NaN, "Range calculation did not work properly");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureNaNThree(){
        ParamUtils.inRange(Double.NaN, Double.NaN, Double.NaN, "Range calculation did not work properly");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureAllPositiveDouble(){
        ParamUtils.inRange(4.1, 7.2, 10, "Range calculation did not work properly");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureMinNegativeDouble(){
        ParamUtils.inRange(45.2, -7.2, 10, "Range calculation did not work properly");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureValNegativeDouble(){
        ParamUtils.inRange(-10.3, -7.2, 10, "Range calculation did not work properly");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureValNegativeLong(){
        ParamUtils.inRange(-10L, -7, 10, "Range calculation did not work properly");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureValNegativeInt(){
        ParamUtils.inRange(-10, -7, 10, "Range calculation did not work properly");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureValNegativeInt2(){
        ParamUtils.inRange(-10, -7.2, 10, "Range calculation did not work properly");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureIntValMinMaxUnreasonable(){
        ParamUtils.inRange(5, 10, -7, "Range calculation did not work properly.  Min was greater than max, so will always throw exception.");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInRangeFailureValMinMaxUnreasonable(){
        ParamUtils.inRange(5L, 10, -7, "Range calculation did not work properly.  Min was greater than max, so will always throw exception.");
    }

    @Test
    public void testInRangeLongDoubleDouble(){
        ParamUtils.inRange(5, -1.2, 7.1, "Range calculation did not work properly!");
    }

    @Test(expectedExceptions = IllegalArgumentException.class, expectedExceptionsMessageRegExp = "That is infinite!")
    public void testInfinite(){
        ParamUtils.isFinite(Double.POSITIVE_INFINITY, "That is infinite!");
    }

    @Test(expectedExceptions = IllegalArgumentException.class, expectedExceptionsMessageRegExp = "That is infinite!")
    public void testNegInfinite(){
        ParamUtils.isFinite(Double.NEGATIVE_INFINITY, "That is infinite!");
    }

    @Test(expectedExceptions = IllegalArgumentException.class, expectedExceptionsMessageRegExp = "That is negative!")
    public void testNegativeIsNotPositive(){
        ParamUtils.isPositiveOrZero(Double.NEGATIVE_INFINITY, "That is negative!");
    }

    @Test(expectedExceptions = IllegalArgumentException.class, expectedExceptionsMessageRegExp = "That is negative!")
    public void testNegativeIsGreaterThanZero(){
        ParamUtils.isPositive(Double.NEGATIVE_INFINITY, "That is negative!");
    }

    @Test(expectedExceptions = IllegalArgumentException.class, expectedExceptionsMessageRegExp = "ERROR")
    public void testNaNIsGreaterThanZero(){
        ParamUtils.isPositive(Double.NaN, "ERROR");
    }

    @Test(expectedExceptions = IllegalArgumentException.class, expectedExceptionsMessageRegExp = "That is negative!")
    public void testNegativeIsNotPositive2(){
        ParamUtils.isPositiveOrZero(-4.2, "That is negative!");
    }

    @Test(expectedExceptions = IllegalArgumentException.class, expectedExceptionsMessageRegExp = "ERROR")
    public void testIsPositiveOnNan(){
        ParamUtils.isPositiveOrZero(Double.NaN, "ERROR");
    }


    @Test(expectedExceptions = IllegalArgumentException.class, expectedExceptionsMessageRegExp = "That is negative!")
    public void testNegativeIsGreaterThanZero2(){
        ParamUtils.isPositive(-4.2, "That is negative!");
    }

    @Test
    public void testPositiveIsGreaterThanZero(){
        ParamUtils.isPositive(4.2, "That is negative!");
    }

    @Test
    public void testPositiveIsGreaterThanZero2(){
        Assert.assertEquals(ParamUtils.isPositive(4L, "That is negative!"), 4L);
    }

    @Test
    public void testPositiveIntIsGreaterThanZero2() {
        Assert.assertEquals(ParamUtils.isPositive(4, "That is negative!"), 4);
    }

    @Test
    public void testPositiveIsGreaterThanZeroWithInf(){
        ParamUtils.isPositive(Double.POSITIVE_INFINITY, "That is negative!");
    }

    @Test(expectedExceptions = IllegalArgumentException.class, expectedExceptionsMessageRegExp = "That is zero!")
    public void testZeroIsGreaterThanZero2(){
        ParamUtils.isPositive(0L, "That is zero!");
    }

    @Test(expectedExceptions = IllegalArgumentException.class, expectedExceptionsMessageRegExp = "That is zero!")
    public void testIntZeroIsGreaterThanZero2(){
        ParamUtils.isPositive(0, "That is zero!");
    }

    @Test(expectedExceptions = IllegalArgumentException.class, expectedExceptionsMessageRegExp = "That is zero!")
    public void testZeroIsGreaterThanZero3(){
        ParamUtils.isPositive(0.0, "That is zero!");
    }

    @Test
    public void testLogB(){
        Assert.assertEquals(ParamUtils.logb(4, 2), 2.0);
        Assert.assertEquals(ParamUtils.logb(4.1, 2), 2.0356, 1e-3);
        Assert.assertEquals(ParamUtils.logb(4.1, 2.1), 1.9018, 1e-3);
        Assert.assertEquals(ParamUtils.logb(4, 2.1), 1.8685, 1e-3);

    }

    @Test
    public void testLog2(){
        Assert.assertEquals(2.0, ParamUtils.log2(4));
        Assert.assertEquals(10.0, ParamUtils.log2(1024));
        Assert.assertEquals(-13.2877, ParamUtils.log2(.0001), 1e-3);
    }

    @Test
    public void testLogBOfBaseZero(){
        Assert.assertEquals(-0.0, ParamUtils.logb(337, -0));
        Assert.assertEquals(-0.0, ParamUtils.logb(337, 0));
        Assert.assertTrue(Double.isNaN(ParamUtils.logb(-337, 0)));
        Assert.assertTrue(Double.isNaN(ParamUtils.logb(-337, -0)));
        Assert.assertTrue(Double.isNaN(ParamUtils.logb(0, 0)));
        Assert.assertTrue(Double.isNaN(ParamUtils.logb(-0, -0)));
        Assert.assertTrue(Double.isNaN(ParamUtils.logb(-0, 0)));
        Assert.assertTrue(Double.isNaN(ParamUtils.logb(0, -0)));
    }

    @Test
    public void testLogBOfZero(){
        Assert.assertEquals(ParamUtils.logb(0, 337), Double.NEGATIVE_INFINITY);
        Assert.assertEquals(ParamUtils.logb(-0, 337), Double.NEGATIVE_INFINITY);
    }

    @Test
    public void testLogBOfNegativeNumber() {
        // Note you cannot do Double.NaN == val  see: http://stackoverflow.com/questions/1456566/how-do-you-test-to-see-if-a-double-is-equal-to-nan
        Assert.assertTrue(Double.isNaN(ParamUtils.logb(-4, 5)));
        Assert.assertTrue(Double.isNaN(ParamUtils.logb(4, -5.1)));
    }

    @Test
    public void testLog2OfNegativeNumber() {
        Assert.assertTrue(Double.isNaN(ParamUtils.log2(-4)));
        Assert.assertTrue(Double.isNaN(ParamUtils.log2(-5.1)));
    }

    @Test
    public void testWriteAndReadDoublesToFile() {
        final double[] testData = new double[]{0.0, 1.234552345, Math.PI, Double.NaN, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, -3.91};
        final File outputFile = IOUtils.createTempFile("param-test-write-doubles-", ".txt");
        ParamUtils.writeValuesToFile(testData, outputFile);

        // Test values
        final double[] valuesInFile = ParamUtils.readValuesFromFile(outputFile);

        Assert.assertEquals(valuesInFile.length, testData.length);

        for (int i = 0; i < testData.length; i++) {
            if (Double.isNaN(testData[i])) {
                Assert.assertTrue(Double.isNaN(valuesInFile[i]));
            } else {
                Assert.assertEquals(valuesInFile[i], testData[i], 1e-20);
            }
        }
    }
}
