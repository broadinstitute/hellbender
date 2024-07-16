package org.broadinstitute.hellbender.engine;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.beta.io.bundle.Bundle;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.Serializable;
import java.util.*;

/**
 * Class to represent a Feature-containing input file. Tools should declare @Argument-annotated fields of
 * this type (or Collections of this type), and the Feature management system will automatically discover
 * them at runtime (provided that they are declared in the tool itself, a superclass of the tool, or an
 * ArgumentCollection of the tool).
 *
 * DO NOT ATTEMPT TO INSTANTIATE THIS CLASS DIRECTLY! FeatureInputs must be instantiated by the argument-parsing
 * system only in order to be recognized by the Feature management system. This is why the constructor is
 * marked as protected.
 *
 * If you still want to instantiate this class directly, you will have to call {@link GATKTool#addFeatureInputsAfterInitialization(String, String, Class, int)}
 *  in order to register the FeatureInput with the engine.
 *
 * FeatureInputs can be assigned logical names on the command line using the syntax:
 *
 *     --argument_name:logical_name feature_file
 *
 * These logical names can then be retrieved by the tool at runtime via {@link #getName}
 *
 * Furthermore, a list of comma-separated key=value pairs may be provided as follows:
 *
 *     --argument_name:logical_name,key1=value1,key2=value2 feature_file
 *
 * the string value provided for a given key can be retrieved via {@link #getAttribute(String)}. Keys must be unique.
 *
 * @param <T> the type of Feature that this FeatureInput file contains (eg., VariantContext, BEDFeature, etc.)
 */
public final class FeatureInput<T extends Feature> extends GATKPath implements Serializable {
    final protected static Logger logger = LogManager.getLogger(FeatureInput.class);

    private static final long serialVersionUID = 1L;

    /**
     * File containing Features as specified by the user on the command line
     */

    /**
     * Cache the codec for this feature input the first time we discover it, so we only do it once
     */
    private transient Class<FeatureCodec<T, ?>> featureCodecClass;

    /**
     * retain any containing bundle in case we need to extract other resources from it
     */
    private Bundle parentBundle;

    /**
     * Delimiter between the logical name and the file name in the --argument_name logical_name:feature_file syntax
     */
    public static final String FEATURE_ARGUMENT_TAG_DELIMITER = ":";

    /**
     * Construct a FeatureInput from a raw String argument value. To specify a logical name or tags, use
     * {@link #FeatureInput(String, String)} or {@link #FeatureInput( String, String, Map<String, String>)}.
     *
     * Only meant to be called by the argument parsing system, and therefore marked as package-visible --
     * FeatureInputs constructed some other way will not be recognized by the engine.
     *
     * Note: cannot delegate to another constructor because Java only allows a call to "this" on the first line of a constructor.
     *
     * @param rawArgumentValue String of the form "logical_name:feature_file" or "feature_file"
     */
    FeatureInput(final String rawArgumentValue) {
        super(rawArgumentValue);
        Utils.nonNull(rawArgumentValue, "rawArgumentValue");
    }

    /**
     * Construct a FeatureInput from a GATKPath argument value.
     *
     * @param gatkPath an existing GATKPath object
     */
    public FeatureInput(final GATKPath gatkPath) {
        super(gatkPath);
    }

    /**
     * Construct a FeatureInput from a GATKPath argument value.
     *
     * @param gatkPath an existing GATKPath object
     * @param name the tag name for this feature input
     */
    public FeatureInput(final GATKPath gatkPath, final String name) {
        super(gatkPath);
        Utils.nonNull(name, "name");
        if (gatkPath.getTag() != null) {
            logger.warn(String.format(
                    "FeatureInput: user-provided tag name %s will be replaced with %s",
                    gatkPath.getTag(),
                    name));
        }

        setTag(name);
    }

    /**
     * Construct a FeatureInput from a path and a name
     *
     * This constructor is meant to be called only by the engine and test classes,
     * which is why it has package access.
     */
    FeatureInput(final String featurePath, final String name) {
        this(featurePath, name, Collections.emptyMap());
    }

    /**
     * Construct a FeatureInput from raw components: name, key value pairs and the file.
     *
     * This constructor is meant to be called by the engine and test classes --
     * FeatureInputs constructed some other way will not be recognized by the engine.
     */
    @VisibleForTesting
    public FeatureInput(final String rawInputSpecifier, final String name, final Map<String, String> keyValueMap) {
        super(rawInputSpecifier);

        Utils.nonNull(name, "name");
        Utils.nonNull(keyValueMap, "kevValueMap");
        Utils.nonNull(rawInputSpecifier, "feature-file");

        setTag(name);
        setTagAttributes(keyValueMap);
    }

    /**
     * Construct a FeatureInput from a Bundle.
     *
     * @param primaryResourcePath the path for the primary feature resource for this bundle
     * @param featureBundle an existing Bundle object; resources in this bundle MUST be IOPathBundleResources (that is,
     *                      they must be backed by an IOPath, not an in-memory object)
     * @param name the tag name for this feature input - may be null
     */
    public FeatureInput(
            final GATKPath primaryResourcePath,
            final Bundle featureBundle,
            final String name) {
        super(primaryResourcePath);
        // retain the containing bundle for later so we can interrogate it for other resources, like the index
        this.parentBundle = featureBundle;
        if (name != null) {
            if (primaryResourcePath.getTag() != null) {
                logger.warn(String.format(
                        "FeatureInput: user-provided tag name %s will be replaced with %s",
                        primaryResourcePath.getTag(),
                        name));
            }
            setTag(name);
        }

    }


    /**
     * Remember the FeatureCodec class for this input the first time it is discovered so we can bypass dynamic codec
     * discovery when multiple FeatureDataSources are created for the same input.
     */
    public void setFeatureCodecClass(final Class<FeatureCodec<T, ?>> featureCodecClass) {
        this.featureCodecClass = featureCodecClass;
    }

    /**
     * @return The previously established FeatureCodec class to use for this input, if any. May return {@code null}.
     */
    public Class<FeatureCodec<T, ?>> getFeatureCodecClass() {
        return this.featureCodecClass;
    }

    /**
     * @return the parent bundle for this FeatureInput, if this input was derived from a Bundle.  May
     * return {@code null}. The returned bundle can be interrogated for companion resources.
     */
    public Bundle getParentBundle() {
        return parentBundle;
    }

    /**
     * creates a name from the given filePath by finding the absolute path of the given input
     */
    private String makeIntoAbsolutePath() {
        if(IOUtils.isGenomicsDBPath(this)){
            return IOUtils.getAbsolutePathWithGenomicsDBURIScheme(this);
        } else if (getScheme() != null && !getScheme().equals("file")) { // local files always have a "file" scheme
            return toPath().toAbsolutePath().toUri().toString();
        } else {
            return getURI().getPath();
        }
    }

    /**
     * Gets the value for the given key associated with this Feature source or {@code null}
     * if no value is associated with a given key.
     * @throws IllegalArgumentException if the key is {@code null}.
     */
    public String getAttribute(final String key) {
        Utils.nonNull(key);
        return getTagAttributes().get(key);
    }

    /**
     * Gets the logical name of this Feature source. This will be a user-provided value if the
     * --argument_name logical_name:feature_file was used on the command line, otherwise it will
     * default to the absolute path of the backing file
     *
     * @return logical name of this source of Features
     */
    public String getName() {
        return getTag() != null ? getTag() : makeIntoAbsolutePath();
    }

    /**
     *
     * @return true if the value for name does not match the default, indicating it was a user supplied name (i.e. foo:file.vcf)
     */
    public boolean hasUserSuppliedName() { return getTag() != null;
    }

    /**
     * Gets the file backing this source of Features
     *
     * @return file backing this source of Features
     */
    //TODO: this should go away; most consumers of this method assume this is a path to a File;
    public String getFeaturePath() {
        return getRawInputString();
    }

    /**
     * FeatureInputs will be hashed by the engine, so make an effort to produce a reasonable hash code
     *
     * @return hash code for this FeatureInput (combination of hash code of the name and file)
     */
    @Override
    public int hashCode() {
        return super.hashCode() + 31 * getRawInputString().hashCode();
    }

    /**
     * Test FeatureInputs for equality
     *
     * @param other object to compare this FeatureInput with
     * @return true if this FeatureInput equals other, otherwise false
     */
    @Override
    public boolean equals(final Object other) {
        if (! (other instanceof FeatureInput)) {
            return false;
        }

        final FeatureInput<?> otherFeature = (FeatureInput<?>)other;
        return super.equals(otherFeature) &&
                Objects.equals(getRawInputString(), otherFeature.getRawInputString());
    }

    /**
     * Returns a String representation of this FeatureInput. Will be the absolute path to
     * the featureFile if we have no logical name, or a String of the form
     * "logical_name:absolute_path_to_featureFile" if we do have a logical name.
     *
     * @return String representation of this FeatureInput
     */
    @Override
    public String toString() {
        final String featureFilePath = makeIntoAbsolutePath();
        return hasUserSuppliedName() ?
                String.format("%s%s%s", getTag(), FEATURE_ARGUMENT_TAG_DELIMITER, featureFilePath):
                featureFilePath;
    }
}
