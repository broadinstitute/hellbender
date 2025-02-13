package org.broadinstitute.hellbender.utils.recalibration.covariates;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

public final class StandardCovariateListUnitTest extends GATKBaseTest {

    public StandardCovariateList makeCovariateList() {
        return new StandardCovariateList(new RecalibrationArgumentCollection(), Collections.singletonList("readGroup"));
    }

    @Test
    public void testSize() {
        StandardCovariateList scl = makeCovariateList();
        Assert.assertEquals(scl.size(), 4);
    }

    @Test
    public void testCovariateNames() {
        StandardCovariateList scl = makeCovariateList();
        Assert.assertEquals(scl.covariateNames(), "ReadGroupCovariate,QualityScoreCovariate,ContextCovariate,CycleCovariate");
    }

    @Test
    public void testIterator() {
        StandardCovariateList scl = makeCovariateList();
        Assert.assertEquals(Utils.stream(scl).count(), 4);
    }

    @Test
    public void testGetCovariates() {
        StandardCovariateList scl = makeCovariateList();
        Assert.assertEquals(scl.getReadGroupCovariate().parseNameForReport(), "ReadGroup");
        Assert.assertEquals(scl.getQualityScoreCovariate().parseNameForReport(), "QualityScore");
        final List<Covariate> additionalCovars = Utils.stream(scl.getAdditionalCovariates()).collect(Collectors.toList());
        Assert.assertEquals(additionalCovars.get(0).parseNameForReport(), "Context");
        Assert.assertEquals(additionalCovars.get(1).parseNameForReport(), "Cycle");
    }

    @Test
    public void testGetCovariatesByIndex() {
        StandardCovariateList scl = makeCovariateList();
        Assert.assertEquals(scl.get(0).parseNameForReport(), "ReadGroup");
        Assert.assertEquals(scl.get(1).parseNameForReport(), "QualityScore");
        Assert.assertEquals(scl.get(2).parseNameForReport(), "Context");
        Assert.assertEquals(scl.get(3).parseNameForReport(), "Cycle");
    }

    @Test(expectedExceptions = IndexOutOfBoundsException.class)
    public void testGetCovariatesByIndexInvalid() {
        StandardCovariateList scl = makeCovariateList();
        scl.get(4);
    }

    @Test
    public void testGetCovariatesByIndexClass() {
        StandardCovariateList scl = makeCovariateList();
        Assert.assertEquals(scl.indexByClass(ReadGroupCovariate.class), 0);
        Assert.assertEquals(scl.indexByClass(QualityScoreCovariate.class), 1);
        Assert.assertEquals(scl.indexByClass(ContextCovariate.class), 2);
        Assert.assertEquals(scl.indexByClass(CycleCovariate.class), 3);

        //finally, test an anonymous subclass
        Assert.assertEquals(scl.indexByClass(new Covariate() {
            private static final long serialVersionUID = 1L;

            @Override
            public void recordValues(GATKRead read, SAMFileHeader header, PerReadCovariateMatrix values, boolean recordIndels) {

            }

            @Override
            public String formatKey(int key) {
                return null;
            }

            @Override
            public int keyFromValue(Object value) {
                return 0;
            }

            @Override
            public int maximumKeyValue() {
                return 0;
            }
        }.getClass()), -1);
    }

    @Test
    public void testGetCovariatesByParsedName() {
        StandardCovariateList scl = makeCovariateList();
        final String[] parsedNames = {"ReadGroup", "QualityScore", "Context", "Cycle"};
        for (String parsedName : parsedNames) {
            Assert.assertEquals(scl.getCovariateByParsedName(parsedName).parseNameForReport(), parsedName);
        }
        Assert.assertEquals(scl.getCovariateByParsedName("fred"), null);
    }

    @Test
    public void testCovariateInitialize() {
        StandardCovariateList scl = makeCovariateList();
        //this just tests non blowing up.
    }
}
