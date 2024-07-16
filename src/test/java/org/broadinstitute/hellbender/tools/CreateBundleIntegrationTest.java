package org.broadinstitute.hellbender.tools;

import htsjdk.beta.io.bundle.*;
import htsjdk.beta.plugin.IOUtils;
import htsjdk.beta.plugin.variants.VariantsBundle;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CreateBundleIntegrationTest extends CommandLineProgramTest {

    // force our local paths to use absolute path names to make BundleResource and IOPath equality checks easier,
    // since all serialized JSON bundles always contain absolute path names for local files
    private final static String LOCAL_VCF = new GATKPath(getTestDataDir() + "/count_variants_withSequenceDict.vcf").getURIString();
    private final static String LOCAL_VCF_IDX = new GATKPath(getTestDataDir() + "/count_variants_withSequenceDict.vcf.idx").getURIString();
    private final static String LOCAL_VCF_GZIP = new GATKPath("src/test/resources/large/NA24385.vcf.gz").getURIString();
    private final static String LOCAL_VCF_TBI = new GATKPath("src/test/resources/large/NA24385.vcf.gz.tbi").getURIString();
    private final static String LOCAL_VCF_WITH_NO_INDEX = new GATKPath("src/test/resources/org/broadinstitute/hellbender/tools/count_variants_withSequenceDict_noIndex.vcf").getURIString();
    private final static String CLOUD_VCF = "gs://hellbender/test/resources/large/1000G.phase3.broad.withGenotypes.chr20.10100000.vcf";
    private final static String CLOUD_VCF_IDX = "gs://hellbender/test/resources/large/1000G.phase3.broad.withGenotypes.chr20.10100000.vcf.idx";

    private final static String PRIMARY_CT = "primary_ct";
    private final static String SECONDARY_CT = "secondary_ct";
    private final static String OTHER_CT = "other_ct";

    @DataProvider(name = "bundleCases")
    public Object[][] bundleCases() {
        return new Object[][] {
                // primary, primary tag, secondary, secondary tag, other(s), other tag(s), suppressIndexResolution, expectedBundle

                // common VCF bundle cases, with inferred content types
                {LOCAL_VCF, null, null, null, null, null, true, new VariantsBundle(new GATKPath(LOCAL_VCF))},
                {LOCAL_VCF, null, null, null, null, null, false, new VariantsBundle(new GATKPath(LOCAL_VCF), new GATKPath(LOCAL_VCF_IDX))},
                {LOCAL_VCF, null, LOCAL_VCF_IDX, null, null, null, true, new VariantsBundle(new GATKPath(LOCAL_VCF), new GATKPath(LOCAL_VCF_IDX))},
                {LOCAL_VCF, null, LOCAL_VCF_IDX, null, null, null, false, new VariantsBundle(new GATKPath(LOCAL_VCF), new GATKPath(LOCAL_VCF_IDX))},
                {LOCAL_VCF_WITH_NO_INDEX, null, null, null, null, null, true, new VariantsBundle(new GATKPath(LOCAL_VCF_WITH_NO_INDEX))},
                {LOCAL_VCF_GZIP, null, null, null, null, null, false, new VariantsBundle(new GATKPath(LOCAL_VCF_GZIP), new GATKPath(LOCAL_VCF_TBI))},
                {LOCAL_VCF_GZIP, null, LOCAL_VCF_TBI, null, null, null, true, new VariantsBundle(new GATKPath(LOCAL_VCF_GZIP), new GATKPath(LOCAL_VCF_TBI))},

                {CLOUD_VCF, null, null, null, null, null, true, new VariantsBundle(new GATKPath(CLOUD_VCF))},
                {CLOUD_VCF, null, null, null, null, null, false, new VariantsBundle(new GATKPath(CLOUD_VCF), new GATKPath(CLOUD_VCF_IDX))},

                {CLOUD_VCF, null, CLOUD_VCF_IDX, null, null, null, true, new VariantsBundle(new GATKPath(CLOUD_VCF), new GATKPath(CLOUD_VCF_IDX))},
                {CLOUD_VCF, null, CLOUD_VCF_IDX, null, null, null, false, new VariantsBundle(new GATKPath(CLOUD_VCF), new GATKPath(CLOUD_VCF_IDX))},

                // common vcf bundle cases, with explicit content types
                {LOCAL_VCF, BundleResourceType.CT_VARIANT_CONTEXTS, null, null, null, null, true, new VariantsBundle(new GATKPath(LOCAL_VCF))},
                {LOCAL_VCF, BundleResourceType.CT_VARIANT_CONTEXTS, null, null, null, null, false, new VariantsBundle(new GATKPath(LOCAL_VCF), new GATKPath(LOCAL_VCF_IDX))},
                {LOCAL_VCF, BundleResourceType.CT_VARIANT_CONTEXTS, LOCAL_VCF_IDX, BundleResourceType.CT_VARIANTS_INDEX, null, null, true, new VariantsBundle(new GATKPath(LOCAL_VCF), new GATKPath(LOCAL_VCF_IDX))},
                {LOCAL_VCF, BundleResourceType.CT_VARIANT_CONTEXTS, LOCAL_VCF_IDX, BundleResourceType.CT_VARIANTS_INDEX, null, null, false, new VariantsBundle(new GATKPath(LOCAL_VCF), new GATKPath(LOCAL_VCF_IDX))},

                // vcf bundle with a vcf, an index, and some other resource with an explicit content type
                {LOCAL_VCF, BundleResourceType.CT_VARIANT_CONTEXTS, LOCAL_VCF_IDX, BundleResourceType.CT_VARIANTS_INDEX, Arrays.asList("someVariantsCompanion.txt"), Arrays.asList("someVariantsCT"), false,
                    new BundleBuilder()
                        .addPrimary(new IOPathResource(new GATKPath(LOCAL_VCF), BundleResourceType.CT_VARIANT_CONTEXTS))
                        .addSecondary(new IOPathResource(new GATKPath(LOCAL_VCF_IDX), BundleResourceType.CT_VARIANTS_INDEX))
                        .addSecondary(new IOPathResource(new GATKPath(new GATKPath("someVariantsCompanion.txt").getURIString()), "someVariantsCT"))
                        .build()},

                // "other" bundles
                {
                        LOCAL_VCF, PRIMARY_CT, null, null, null, null, true,
                        new BundleBuilder()
                                .addPrimary(new IOPathResource(new GATKPath(LOCAL_VCF), PRIMARY_CT))
                                .build()
                },
                {
                        LOCAL_VCF, PRIMARY_CT, LOCAL_VCF_IDX, SECONDARY_CT, null, null, true,
                        new BundleBuilder()
                                .addPrimary(new IOPathResource(new GATKPath(LOCAL_VCF), PRIMARY_CT))
                                .addSecondary(new IOPathResource(new GATKPath(LOCAL_VCF_IDX), SECONDARY_CT))
                                .build()
                },
                {
                        // frankenbundle with multiple resources
                        LOCAL_VCF, PRIMARY_CT, LOCAL_VCF_IDX, SECONDARY_CT, Arrays.asList(LOCAL_VCF_TBI), Arrays.asList(OTHER_CT), true,
                        new BundleBuilder()
                                .addPrimary(new IOPathResource(new GATKPath(LOCAL_VCF), PRIMARY_CT))
                                .addSecondary(new IOPathResource(new GATKPath(LOCAL_VCF_IDX), SECONDARY_CT))
                                .addSecondary(new IOPathResource(new GATKPath(LOCAL_VCF_TBI), OTHER_CT))
                                .build()
                },
        };
    }

    @DataProvider(name = "negativeBundleCases")
    public Object[][] negativeBundleCases() {
        return new Object[][] {
                // primary, primary tag, secondary, secondary tag, other(s), other tag(s), suppressIndexResolution, expectedBundle

                // no index file can be inferred
                {LOCAL_VCF_WITH_NO_INDEX, null, null, null, null, null, false, new VariantsBundle(new GATKPath(LOCAL_VCF_WITH_NO_INDEX))},
                // primary content type not provided and cannot be inferred from the extension
                {"primaryFile.ext", null, null, null, null, null, false, null},
                // secondary content type not provided
                {"primaryFile.ext", PRIMARY_CT, "secondaryFile.ext", null, null, null, false, null},

                // other bundle with other content type not provided
                {"primaryFile.ext", PRIMARY_CT, "secondaryFile.ext", SECONDARY_CT, Arrays.asList("other.txt"), null, false, null},
                // vcf bundle with other content type not provided
                {LOCAL_VCF, BundleResourceType.CT_VARIANT_CONTEXTS, null, null, Arrays.asList("other.txt"), null, false, null},
        };
    }

    @Test(dataProvider = "bundleCases")
    public void testBundleCases(
            final String primaryInput,
            final String primaryInputTag,
            final String secondaryInput,
            final String secondaryInputTag,
            final List<String> otherInputs,
            final List<String> otherInputTags,
            final boolean resolveIndex,
            final Bundle expectedBundle) {
        doCreateBundleTest (primaryInput, primaryInputTag, secondaryInput, secondaryInputTag, otherInputs, otherInputTags, resolveIndex, expectedBundle);
    }

    @Test(dataProvider = "negativeBundleCases", expectedExceptions = IllegalArgumentException.class)
    public void testNegativeBundleCases(
            final String primaryInput,
            final String primaryInputTag,
            final String secondaryInput,
            final String secondaryInputTag,
            final List<String> otherInputs,
            final List<String> otherInputTags,
            final boolean resolveIndex,
            final Bundle expectedBundle) {
        doCreateBundleTest (primaryInput, primaryInputTag, secondaryInput, secondaryInputTag, otherInputs, otherInputTags, resolveIndex, expectedBundle);
    }

    private void doCreateBundleTest(
            final String primaryInput,
            final String primaryInputTag,
            final String secondaryInput,
            final String secondaryInputTag,
            final List<String> otherInputs,
            final List<String> otherInputTags,
            final boolean suppressIndexResolution,
            final Bundle expectedBundle) {
        final GATKPath outputPath = new GATKPath(createTempFile("test", ".bundle.json").getAbsolutePath().toString());

        final List<String> args = new ArrayList<>();

        args.add("--" + StandardArgumentDefinitions.PRIMARY_INPUT_LONG_NAME + (primaryInputTag != null ? ":" + primaryInputTag : ""));
        args.add(primaryInput);
        if (secondaryInput != null) {
            args.add("--" + StandardArgumentDefinitions.SECONDARY_INPUT_LONG_NAME + (secondaryInputTag != null ? ":" + secondaryInputTag : ""));
            args.add(secondaryInput);
        }
        if (otherInputs != null) {
            for (int i = 0; i < otherInputs.size(); i++) {
                args.add("--" + CreateBundle.OTHER_INPUT_FULL_NAME + ((otherInputTags != null && otherInputTags.get(i) != null) ? ":" + otherInputTags.get(i) : ""));
                args.add(otherInputs.get(i) != null ? otherInputs.get(i) : "");
            }
        }
        if (suppressIndexResolution == true) {
            args.add("--" + CreateBundle.SUPPRESS_INDEX_RESOLUTION_FULL_NAME);
        }
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.add(outputPath.toString());

        runCommandLine(args);

        final Bundle actualBundle = BundleJSON.toBundle(IOUtils.getStringFromPath(outputPath), GATKPath::new);

        // bundle resource order is not preserved when roundtripping through JSON, so compare ignoring order
        Assert.assertTrue(Bundle.equalsIgnoreOrder(actualBundle, expectedBundle));
    }
}
