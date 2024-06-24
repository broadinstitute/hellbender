package org.broadinstitute.hellbender.tools.walkers.annotator.flow;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;

import java.util.*;

public class FlowAnnotatorUnitTest extends GATKBaseTest {

    private FlowAnnotatorBase[]     allAnnotators = {
            new IndelClassify(),
            new IndelLength(),
            new HmerIndelLength(),
            new HmerIndelNuc(),
            new HmerMotifs(),
            new GcContent(),
            new CycleSkipStatus(),
            new VariantType()
    };

    @DataProvider(name = "testData")
    public Object[][] getTestData() {

        final Object[][]        testData = {
                // order:
                // refbases, altAllele (include a space before and after refAllele
                // indel-class, indel-length, hmer-indel-lenfth, hmer-indel-nuc
                // left-motif, right-motif, gc-content, cycleskip-status
                // value that starts with "!" means ignore
                {
                        // a simple SNP
                        "GTATC A ACATCGGA", "C",
                        "NA", "", "0", "", "GTATC", "ACATC", "0.3", "non-skip", "snp"
                },
                {
                        // a possible-cycle-skip SNP
                        "GTATC A TCATCGGA", "C",
                        "NA", "", "0", "", "GTATC", "TCATC", "0.3", "possible-cycle-skip", "snp"
                },
                {
                        // a cycle-skip SNP
                        "GTATC A ACATCGGA", "T",
                        "NA", "", "0", "", "GTATC", "ACATC", "0.3", "cycle-skip", "snp"
                },
                {
                        // hmer indel (1->0)
                        "TATCT CA TTGACCAA", "C",
                        "del", "1", "1", "A", "ATCTC", "TTGAC", "0.3", "NA", "h-indel"
                },
                {
                        // del hmer indel
                        "TATCTC AT TGACCAA", "A",
                        "del", "1", "2", "T", "TCTCA", "GACCA", "0.4", "NA", "h-indel"
                },
                {
                        // ins hmer indel
                        "TATCT C ATTGACCAA", "CA",
                        "ins", "1", "2", "A", "ATCTC", "TTGAC", "0.3", "NA", "h-indel"
                },
                //The following tests test the edge case of variant so close to the edge of the chromosome
                //that some of the annotations can't be calculated
                {
                        // Close to 5' edge of the chromosome, expect not to crash on finding the left motif
                        "ATCT C ATTGACCAA", "CA",
                        "ins", "1", "2", "A", null, "TTGAC", "0.3", "NA", "h-indel"
                },
                {
                        // close to 3' edge of the chromosome, should not find left motif or GC content
                        "TATCT C ATT", "CA",
                        "ins", "1", "2", "A", "ATCTC", null, null, "NA", "h-indel"
                },
                {
                        // very close to 3' edge of the chromosome, should not find left motif, GC content or if it is
                        // an hmer indel
                        "TATCT C ", "CA",
                        "ins", "1", null, null, "ATCTC", null, null, "NA", "non-h-indel"
                }



        };

        return testData;
    }

    @Test(dataProvider = "testData")
    public void testBasic(final String refBases, final String altAllele,
                          final String indelClass, final String indelLength,
                          final String hmerIndelLength, final String hmerIndelNuc,
                          final String leftMotif, final String rightMotif,
                          final String gcContent, final String cycleskipStatus, final String variantType) {

        // should be in same order as test data!!!!
        final List<String>      expectedAttrs = allKeys();

        // prepare
        final int        refAlleleStart = refBases.indexOf(' ');
        final int        refAlleleEnd = refBases.indexOf(' ', refAlleleStart + 1);
        final String     refAllele = refBases.substring(refAlleleStart + 1, refAlleleEnd);
        final ReferenceContext ref = buildReferenceContext(refBases.replace(" ", ""), refAlleleStart + 1, refAlleleEnd - 1);
        final VariantContext vc = buildVariantContext(ref, refAllele, altAllele);
        String          msg = "on " + refBases + " " + altAllele;

        // invoke
        final Map<String, Object> attrs = allAnnotate(ref, vc);
        Assert.assertNotNull(attrs, msg);

        // check that all expected attributes are there
        String[]        testResults = {indelClass, indelLength, hmerIndelLength, hmerIndelNuc,
                leftMotif, rightMotif, gcContent, cycleskipStatus, variantType};
        for ( int n = 0 ; n < 9 ; n++ ) {
            String       key = expectedAttrs.get(n);
            String       elem = testResults[n];
            String       keyMsg = "on " + key + " " + msg;
            if ( elem != null && !elem.startsWith("!") ) {
                Object v = attrs.get(key);
                if (v instanceof List) {
                    v = StringUtils.join((List) v, ",");
                }
                Assert.assertEquals(v.toString(), elem, keyMsg);
            } else if ( elem == null ) {
                Assert.assertFalse(attrs.containsKey(key), keyMsg);
            }
        }
    }

    private Map<String, Object> allAnnotate(final ReferenceContext ref, final VariantContext vc) {

        final Map<String, Object>     attrs = new LinkedHashMap<>();

        for ( FlowAnnotatorBase a : allAnnotators ) {
            attrs.putAll(a.annotate(ref, vc, null));
        }

        return attrs;
    }

    private List<String> allKeys() {

        List<String>     keys = new LinkedList<>();

        for ( FlowAnnotatorBase a : allAnnotators ) {
            keys.addAll(a.getKeyNames());
            a.setFlowOrder(Collections.singletonList(FlowBasedRead.DEFAULT_FLOW_ORDER));
        }

        return keys;

    }

    private ReferenceContext buildReferenceContext(String refBases, int start, int stop) {

        // note that locations here are 1 based
        final String                insLoc = "chr1";
        final SimpleInterval        interval = new SimpleInterval(insLoc, start, stop);
        final byte[]                refBytes = refBases.getBytes();
        final SimpleInterval        interval1 = new SimpleInterval(insLoc, 1, refBytes.length);
        final ReferenceBases        ref1 = new ReferenceBases(refBytes, interval1);
        final SAMSequenceDictionary dict = new SAMSequenceDictionary(Arrays.asList(new SAMSequenceRecord(insLoc, refBytes.length)));
        final ReferenceContext      ref = new ReferenceContext(ReferenceDataSource.of(ref1, dict), interval, start - 1, 20);

        return ref;
    }

    private VariantContext buildVariantContext(ReferenceContext ref, String refBases, String altBases) {

        final Allele refAllele = Allele.create(refBases, true);
        final Allele altAllele = Allele.create(altBases, false);
        final VariantContext vc = new VariantContextBuilder("foo",
                ref.getContig(), ref.getStart(), ref.getEnd(),
                Arrays.asList(refAllele, altAllele)).make();

        return vc;
    }
}