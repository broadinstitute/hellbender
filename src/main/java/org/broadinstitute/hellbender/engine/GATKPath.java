package org.broadinstitute.hellbender.engine;

import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem;

import htsjdk.io.HtsPath;

import org.broadinstitute.barclay.argparser.TaggedArgument;
import org.broadinstitute.barclay.argparser.TaggedArgumentParser;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.Serializable;
import java.nio.file.FileSystemNotFoundException;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;
import java.util.Objects;

/**
 * GATK tool command line arguments that are input or output resources. These can
 * have an optional name supplied on the command line, as well as one or more optional tag/value pairs.
 */
public class GATKPath extends HtsPath implements TaggedArgument, Serializable {
    private static final long serialVersionUID = 1L;

    public static final String HDFS_SCHEME = "hdfs";
    public static final String S3A_SCHEME = "s3a";

    private String tagName;
    private Map<String, String> tagAttributes;

    /**
     * The raw value for this specifier as provided by the user. Can be a local file or valid URI.
     * @param uriString
     */
    public GATKPath(final String uriString) {
        super(uriString);
    }

    /**
     * Create a GATKPath from an existing GATKPath. Propagates tag and tag attributes.
     * @param sourceGATKPath
     */
    public GATKPath(final GATKPath sourceGATKPath) {
        super(sourceGATKPath);

        if (sourceGATKPath.getTag() != null) {
            this.setTag(sourceGATKPath.getTag());
        }

        final Map<String, String> sourceTagMap = sourceGATKPath.getTagAttributes();
        if (sourceTagMap != null) {
            this.setTagAttributes(new HashMap<>(sourceTagMap));
        }
    }

    @Override
    public Path toPath() {
        // special case GCS, in case the filesystem provider wasn't installed properly but is available.
        if (CloudStorageFileSystem.URI_SCHEME.equals(getURI().getScheme())) {
            final Path tempPath = BucketUtils.getPathOnGcs(getURIString());
            setCachedPath(tempPath);
            return tempPath;
        } else {
            try {
                return super.toPath();
            } catch (FileSystemNotFoundException e) {
                try {
                    ClassLoader cl = Thread.currentThread().getContextClassLoader();
                    if (cl == null) {
                        throw e;
                    }
                    return FileSystems.newFileSystem(getURI(), new HashMap<>(), cl).provider().getPath(getURI());
                } catch (final IOException ioe) {
                    throw new GATKException("IOException loading new file system", ioe);
                }
            }
        }
    }

    @Override
    public InputStream getInputStream() {
        if (!isPath()) {
            throw new UserException(getToPathFailureReason());
        }

        try {
            InputStream inputStream;
            if (isHadoopURL()) {
                org.apache.hadoop.fs.Path file = new org.apache.hadoop.fs.Path(getURIString());
                org.apache.hadoop.fs.FileSystem fs = file.getFileSystem(new org.apache.hadoop.conf.Configuration());
                inputStream = fs.open(file);
            } else {
                inputStream = super.getInputStream();
            }

            return inputStream;
        } catch (IOException e) {
            throw new UserException.CouldNotReadInputFile(this, "Can't create input stream", e);
        }
    }

    /**
     * Open a path for writing regardless of whether it's on GCS, HDFS or local disk.
     * For writing to GCS it'll use the application/octet-stream MIME type.
     *
     * @return an OutputStream that writes to the path specified by this UR.
     */
    @Override
    public OutputStream getOutputStream() {
        if (!isPath()) {
            throw new UserException(getToPathFailureReason());
        }

        try {
            OutputStream outputStream;
            if (isHadoopURL()) {
                org.apache.hadoop.fs.Path fsPath = new org.apache.hadoop.fs.Path(getURIString());
                org.apache.hadoop.fs.FileSystem fs = fsPath.getFileSystem(new org.apache.hadoop.conf.Configuration());
                return fs.create(fsPath);
            } else {
                outputStream = super.getOutputStream();
            }

            return outputStream;
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile(getRawInputString(), e);
        }
    }

    /**
     * Returns true if this is a HDFS (Hadoop filesystem) URL.
     */
    public boolean isHadoopURL() {
        return getScheme().equals(HDFS_SCHEME) || getScheme().equals(S3A_SCHEME);
    }

    @Override
    public void setTag(String tagName) {
        this.tagName = tagName;
    }

    @Override
    public String getTag() {
        return tagName;
    }

    @Override
    public void setTagAttributes(final Map<String, String> attributes) {
            this.tagAttributes = attributes;
    }

    @Override
    public Map<String, String> getTagAttributes() {
        return tagAttributes;
    }

    @Override
    public String toString() {
        if (getTagAttributes() != null) {
            return super.toString() + TaggedArgumentParser.getDisplayString("", this);
        } else {
            return super.toString();
        }
    }

    /**
     * Return true if the tag name and attributes for two {@code TaggedArgument} objects are the equal. For use
     * by TaggedArgument implementations. Canbe replaced with the Barclay version once
     * https://github.com/broadinstitute/barclay/pull/152 is in.
     *
     * @param first first {@code TaggedArgument} to compare
     * @param second second {@code TaggedArgument}
     * @return true if the tag name and attributes are equal
     */
    static boolean tagNameAndAttributesAreEqual(final TaggedArgument first, final TaggedArgument second) {
        return Objects.equals(first.getTag(), second.getTag()) &&
                Objects.equals(first.getTagAttributes(), second.getTagAttributes());
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof GATKPath)) return false;
        if (!super.equals(o)) return false;

        return tagNameAndAttributesAreEqual(this, (GATKPath) o);
    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        result = tagName == null ?
                result :
                31 * result + tagName.hashCode();
        result = getTagAttributes() == null ?
                result :
                31 * result + getTagAttributes().hashCode();
        return result;
    }
}
