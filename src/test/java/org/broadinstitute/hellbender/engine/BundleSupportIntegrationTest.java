package org.broadinstitute.hellbender.engine;

import htsjdk.beta.io.IOPathUtils;
import htsjdk.beta.io.bundle.Bundle;
import htsjdk.beta.io.bundle.BundleJSON;
import htsjdk.io.IOPath;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.IOException;

public class BundleSupportIntegrationTest extends GATKBaseTest {

    // this test uses a serialized bundle file to ensure that we don't unintentionally pick up any
    // code (like, from htsjdk) that introduces backward compatibility issues
    @Test
    public void testReadWriteSerializedReferenceBundle() throws IOException {
        // This test file contains absolute paths to files on a local dev machine, so it shouldn't really be used
        // for anything other than this test, since the absolute paths are unlikely to work on any other machine.
        // But here we just want to make sure we can consume and roundtrip it without error
        final IOPath testBundleFilePath = new GATKPath("src/test/resources/org/broadinstitute/hellbender/engine/print_reads_bundle_do_not_use.json");

        // get our test bundle from the file (ensure we canparse it), then write it out to a temp file, read it back
        // in, and compare
        final Bundle testBundle = BundleJSON.toBundle(IOPathUtils.getStringFromPath(testBundleFilePath));
        final IOPath roundTrippedBundleFilePath = new GATKPath(
                createTempPath("testReadWriteSerializedReferenceBundle", ".json").toString());
        IOPathUtils.writeStringToPath(roundTrippedBundleFilePath, BundleJSON.toJSON(testBundle));
        final Bundle roundTrippedBundle = BundleJSON.toBundle(IOPathUtils.getStringFromPath(testBundleFilePath));
        Assert.assertTrue(Bundle.equalsIgnoreOrder(roundTrippedBundle, testBundle));
    }

}
