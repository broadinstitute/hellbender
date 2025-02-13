package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.util.ArrayList;
import java.util.List;

/**
 * This class knows about all FeatureOutputCodec implementations, and allows you to find an
 * appropriate codec to create a given file type.
 */
public final class FeatureOutputCodecFinder {

    private static final List<FeatureOutputCodec<? extends Feature, ? extends FeatureSink<?>>> outputCodecs =
            new ArrayList<>(10);
    static {
        outputCodecs.add(new BafEvidenceCodec());
        outputCodecs.add(new DepthEvidenceCodec());
        outputCodecs.add(new DiscordantPairEvidenceCodec());
        outputCodecs.add(new SiteDepthCodec());
        outputCodecs.add(new SplitReadEvidenceCodec());
        outputCodecs.add(new BafEvidenceBCICodec());
        outputCodecs.add(new DepthEvidenceBCICodec());
        outputCodecs.add(new DiscordantPairEvidenceBCICodec());
        outputCodecs.add(new SiteDepthBCICodec());
        outputCodecs.add(new SplitReadEvidenceBCICodec());
    }

    @SuppressWarnings("unchecked")
    public static <F extends Feature>
            FeatureOutputCodec<F, FeatureSink<F>> find( final GATKPath outputFilePath,
                                                        final Class<F> featureClass ) {
        FeatureOutputCodec<F, FeatureSink<F>> result = null;
        final String outputFileName = outputFilePath.getRawInputString();
        for ( final FeatureOutputCodec<?, ?> codec : outputCodecs ) {
            if ( codec.canDecode(outputFileName) ) {
                if ( result != null ) {
                    throw new GATKException("Found multiple output codecs for " + outputFileName);
                }
                if ( featureClass.isAssignableFrom(codec.getFeatureType()) ) {
                    result = (FeatureOutputCodec<F, FeatureSink<F>>)codec;
                }
            }
        }
        if ( result == null ) {
            throw new UserException("No feature output codec that writes " +
                    featureClass.getSimpleName() + " was found for " + outputFileName);
        }
        return result;
    }
}
