package org.broadinstitute.hellbender.utils.gcs;

import com.google.cloud.hadoop.gcsio.GoogleCloudStorageFileSystem;
import com.google.cloud.storage.BlobInfo;
import com.google.cloud.storage.HttpMethod;
import com.google.cloud.storage.Storage;
import com.google.cloud.storage.StorageOptions;
import com.google.cloud.storage.contrib.nio.CloudStorageConfiguration;
import com.google.cloud.storage.contrib.nio.CloudStorageConfiguration.Builder;
import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem;
import com.google.cloud.storage.contrib.nio.CloudStorageFileSystemProvider;
import com.google.cloud.storage.contrib.nio.SeekableByteChannelPrefetcher;
import com.google.common.base.Strings;
import com.google.common.io.ByteStreams;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.http.nio.HttpFileSystemProvider;
import org.broadinstitute.http.nio.HttpFileSystemProviderSettings;
import org.broadinstitute.http.nio.HttpsFileSystemProvider;
import com.google.api.gax.retrying.RetrySettings;
import com.google.auth.oauth2.GoogleCredentials;
import com.google.cloud.http.HttpTransportOptions;
import org.broadinstitute.http.nio.RetryHandler;
import org.threeten.bp.Duration;

import java.io.*;
import java.net.URL;
import java.nio.channels.SeekableByteChannel;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.UUID;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;

/**
 * Utilities for dealing with google buckets.
 */
public final class BucketUtils {
    public static final String GCS_PREFIX = GoogleCloudStorageFileSystem.SCHEME + "://";
    public static final String HTTP_PREFIX = HttpFileSystemProvider.SCHEME + "://";
    public static final String HTTPS_PREFIX = HttpsFileSystemProvider.SCHEME +"://";
    public static final String HDFS_SCHEME = "hdfs";
    public static final String HDFS_PREFIX = HDFS_SCHEME + "://";

    // slashes omitted since hdfs paths seem to only have 1 slash which would be weirder to include than no slashes
    public static final String FILE_PREFIX = "file:";
    public static final int MAX_TIMEOUT = 120000;

    private BucketUtils(){} //private so that no one will instantiate this class

    /**
     * @param path path to inspect
     * @return true if this path represents a gcs location
     */
    public static boolean isGcsUrl(final String path) {
        Utils.nonNull(path);
        return path.startsWith(GCS_PREFIX);
    }

    /**
     * Return true if this {@code GATKPath} represents a gcs URI.
     * @param pathSpec specifier to inspect
     * @return true if this {@code GATKPath} represents a gcs URI.
     */
    public static boolean isGcsUrl(final GATKPath pathSpec) {
        Utils.nonNull(pathSpec);
        return pathSpec.getScheme().equals(GoogleCloudStorageFileSystem.SCHEME);
    }

    /**
     * @param pathSpec specifier to inspect
     * @return true if this {@code GATKPath} represents a remote storage system which may benefit from prefetching (gcs or http(s))
     */
    public static boolean isEligibleForPrefetching(final GATKPath pathSpec) {
        Utils.nonNull(pathSpec);
        return isEligibleForPrefetching(pathSpec.getScheme());
     }

    /**
     * @param path path to inspect
     * @return true if this {@code Path} represents a remote storage system which may benefit from prefetching (gcs or http(s))
     */
    public static boolean isEligibleForPrefetching(final java.nio.file.Path path) {
        Utils.nonNull(path);
        return isEligibleForPrefetching(path.toUri().getScheme());
    }

    private static boolean isEligibleForPrefetching(final String scheme){
        return scheme != null
                && (scheme.equals(GoogleCloudStorageFileSystem.SCHEME)
                || scheme.equals(HttpFileSystemProvider.SCHEME)
                || scheme.equals(HttpsFileSystemProvider.SCHEME));
    }

    /**
     * @return true if the given path is an http or https Url.
     */
    public static boolean isHttpUrl(String path){
        return path.startsWith(HTTP_PREFIX) || path.startsWith(HTTPS_PREFIX);
    }

    /**
     * Returns true if the given path is a HDFS (Hadoop filesystem) URL.
     */
    public static boolean isHadoopUrl(String path) {
        return path.startsWith(HDFS_PREFIX);
    }

    /**
     * Returns true if the given path is a GCS, HDFS (Hadoop filesystem), or Http(s) URL.
     */
    public static boolean isRemoteStorageUrl(String path) {
        return isGcsUrl(path) || isHadoopUrl(path) || isHttpUrl(path);
    }

    /**
     * Changes relative local file paths to be absolute file paths. Paths with a scheme are left unchanged.
     * @param path the path
     * @return an absolute file path if the original path was a relative file path, otherwise the original path
     */
    //TODO: get rid of this..
    public static String makeFilePathAbsolute(String path){
        if (isGcsUrl(path) || isHadoopUrl(path) || isFileUrl(path) || isHttpUrl(path)){
            return path;
        } else {
            return new File(path).getAbsolutePath();
        }
    }

    /**
     * Open a file for reading regardless of whether it's on GCS, HDFS or local disk.
     *
     * If the file ends with .gz will attempt to wrap it in an appropriate unzipping stream
     *
     * @param path the GCS, HDFS or local path to read from. If GCS, it must start with "gs://", or "hdfs://" for HDFS.
     * @return an InputStream that reads from the specified file.
     */
    public static InputStream openFile(String path) {
        try {
            Utils.nonNull(path);
            InputStream inputStream;
            if (isGcsUrl(path)) {
                java.nio.file.Path p = getPathOnGcs(path);
                inputStream = Files.newInputStream(p);
            } else if (isHadoopUrl(path)) {
                Path file = new org.apache.hadoop.fs.Path(path);
                FileSystem fs = file.getFileSystem(new Configuration());
                inputStream = fs.open(file);
            } else {
                inputStream = new FileInputStream(path);
            }

            if(IOUtil.hasBlockCompressedExtension(path)){
                return IOUtils.makeZippedInputStream(new BufferedInputStream(inputStream));
            } else {
                return inputStream;
            }
        } catch (IOException x) {
            throw new UserException.CouldNotReadInputFile(path, x);
        }
    }

    /**
     * Open a binary file for writing regardless of whether it's on GCS, HDFS or local disk.
     * For writing to GCS it'll use the application/octet-stream MIME type.
     *
     * @param path the GCS or local path to write to. If GCS, it must start with "gs://", or "hdfs://" for HDFS.
     * @return an OutputStream that writes to the specified file.
     */
    public static OutputStream createFile(String path) {
        Utils.nonNull(path);
        try {
            if (isGcsUrl(path)) {
                java.nio.file.Path p = getPathOnGcs(path);
                return Files.newOutputStream(p);
            } else if (isHadoopUrl(path)) {
                Path file = new Path(path);
                FileSystem fs = file.getFileSystem(new Configuration());
                return fs.create(file);
            } else {
                return new FileOutputStream(path);
            }
        } catch (IOException x) {
            throw new UserException.CouldNotCreateOutputFile("Could not create file at path:" + path + " due to " + x.getMessage(), x);
        }
    }

    /**
     * Copies a file. Can be used to copy e.g. from GCS to local.
     *
     * @param sourcePath the path to read from. If GCS, it must start with "gs://", or "hdfs://" for HDFS.
     * @param destPath the path to copy to. If GCS, it must start with "gs://", or "hdfs://" for HDFS.
     * @throws IOException
     */
    public static void copyFile(String sourcePath, String destPath) throws IOException {
        try (
            InputStream in = openFile(sourcePath);
            OutputStream fout = createFile(destPath)) {
            ByteStreams.copy(in, fout);
        }
    }

    /**
     * Deletes a file: local, GCS or HDFS.
     *  @param pathToDelete the path to delete. If GCS, it must start with "gs://", or "hdfs://" for HDFS.
     *
     */
    public static void deleteFile(String pathToDelete) throws IOException {
        if (isGcsUrl(pathToDelete)) {
            java.nio.file.Path p = getPathOnGcs(pathToDelete);
            Files.delete(p);
        } else if (isHadoopUrl(pathToDelete)) {
            Path file = new Path(pathToDelete);
            FileSystem fs = file.getFileSystem(new Configuration());
            fs.delete(file, false);
        } else {
            boolean ok = new File(pathToDelete).delete();
            if (!ok) throw new IOException("Unable to delete '"+pathToDelete+"'");
        }
    }

    /**
     * Get a temporary file path based on the prefix and extension provided.
     * This file (and possible indexes associated with it) will be scheduled for deletion on shutdown
     *
     * @param prefix a prefix for the file name
     *               for remote paths this should be a valid URI to root the temporary file in (ie. gs://hellbender/staging/)
     *               there is no guarantee that this will be used as the root of the tmp file name, a local prefix may be placed in the tmp folder for example
     * @param extension and extension for the temporary file path, the resulting path will end in this
     * @return a path to use as a temporary file, on remote file systems which don't support an atomic tmp file reservation a path is chosen with a long randomized name
     *
     */
    public static String getTempFilePath(String prefix, String extension){
        if (isGcsUrl(prefix) || (isHadoopUrl(prefix))){
            final String path = randomRemotePath(prefix, "", extension);
            IOUtils.deleteOnExit(IOUtils.getPath(path));
            IOUtils.deleteOnExit(IOUtils.getPath(path + FileExtensions.TRIBBLE_INDEX));
            IOUtils.deleteOnExit(IOUtils.getPath(path + FileExtensions.TABIX_INDEX));
            IOUtils.deleteOnExit(IOUtils.getPath(path + ".bai"));
            IOUtils.deleteOnExit(IOUtils.getPath(path + ".md5"));
            IOUtils.deleteOnExit(IOUtils.getPath(path.replaceAll(extension + "$", ".bai"))); //if path ends with extension, replace it with .bai
            return path;
        } else {
            return IOUtils.createTempFile(prefix, extension).getAbsolutePath();
        }
    }

    /**
     * Picks a random name, by putting some random letters between "prefix" and "suffix".
     *
     * @param stagingLocation The folder where you want the file to be. Must start with "gs://" or "hdfs://"
     * @param prefix The beginning of the file name
     * @param suffix The end of the file name, e.g. ".tmp"
     */
    public static String randomRemotePath(String stagingLocation, String prefix, String suffix) {
        if (isGcsUrl(stagingLocation)) {
            // Go through URI because Path.toString isn't guaranteed to include the "gs://" prefix.
            return getPathOnGcs(stagingLocation).resolve(prefix + UUID.randomUUID().toString() + suffix).toUri().toString();
        } else if (isHadoopUrl(stagingLocation)) {
            return new Path(stagingLocation, prefix + UUID.randomUUID().toString() + suffix).toString();
        } else {
            throw new IllegalArgumentException("Staging location is not remote: " + stagingLocation);
        }
    }

    /**
     * Returns true if we can read the first byte of the file.
     *  @param path The folder where you want the file to be (local, GCS or HDFS).
     *
     */
    public static boolean fileExists(String path) {
        final boolean MAYBE = false;
        try (InputStream inputStream = openFile(path)) {
            int ignored = inputStream.read();
        } catch (UserException.CouldNotReadInputFile notthere) {
            // file isn't there
            return false;
        } catch (FileNotFoundException x) {
            // file isn't there
            return false;
        } catch (IOException x) {
            // unexpected problem while reading the file. The file may exist, but it's not accessible.
            return MAYBE;
        }
        return true;
    }

    /**
     * Returns the file size of a file pointed to by a GCS/HDFS/local path
     *
     * @param path The URL to the file whose size to return
     * @return the file size in bytes
     * @throws IOException
     */
    public static long fileSize(String path) throws IOException {
        if (isGcsUrl(path)) {
            java.nio.file.Path p = getPathOnGcs(path);
            return Files.size(p);
        } else if (isHadoopUrl(path)) {
            Path hadoopPath = new Path(path);
            FileSystem fs = hadoopPath.getFileSystem(new Configuration());
            return fs.getFileStatus(hadoopPath).getLen();
        } else {
            return new File(path).length();
        }
    }

    /**
     * Returns the total file size of all files in a directory, or the file size if the path specifies a file.
     * Note that sub-directories are ignored - they are not recursed into.
     * Only supports HDFS and local paths.
     *
     * @param pathSpecifier The URL to the file or directory whose size to return
     * @return the total size of all files in bytes
     */
    public static long dirSize(final GATKPath pathSpecifier) {
        try {
            // GCS case (would work with local too)
            if (isGcsUrl(pathSpecifier)) {
                java.nio.file.Path p = getPathOnGcs(pathSpecifier.getRawInputString());
                if (Files.isRegularFile(p)) {
                    return Files.size(p);
                }
                return Files.list(p).mapToLong(
                    q -> {
                        try {
                            return (Files.isRegularFile(q) ? Files.size(q) : 0);
                        } catch (IOException e) {
                            throw new RuntimeIOException(e);
                        }
                    }
                ).sum();
            }
            // local file or HDFS case
            Path hadoopPath = new Path(pathSpecifier.getURIString());
            FileSystem fs = new Path(pathSpecifier.getURIString()).getFileSystem(new Configuration());
            FileStatus status = fs.getFileStatus(hadoopPath);
            if (status == null) {
                throw new UserException.CouldNotReadInputFile(pathSpecifier.getRawInputString(), "File not found.");
            }
            long size = 0;
            if (status.isDirectory()) {
                for (FileStatus st : fs.listStatus(status.getPath())) {
                    if (st.isFile()) {
                        size += st.getLen();
                    }
                }
            } else {
                size += status.getLen();
            }
            return size;
        } catch (RuntimeIOException | IOException e) {
            throw new UserException("Failed to determine total input size of " + pathSpecifier.getRawInputString() + "\n Caused by:" + e.getMessage(), e);
        }
    }

    public static boolean isFileUrl(String path) {
        return path.startsWith(FILE_PREFIX);
    }

    /**
     * Given a path of the form "gs://bucket/folder/folder/file", returns "bucket".
     */
    public static String getBucket(String path) {
        return path.split("/")[2];
    }

    /**
     * Given a path of the form "gs://bucket/folder/folder/file", returns "folder/folder/file".
     */
    public static String getPathWithoutBucket(String path) {
        final String[] split = path.split("/");
        return String.join("/", Arrays.copyOfRange(split, 3, split.length));
    }

    /**
     * Sets max_reopens, requester_pays, and generous timeouts as the global default.
     * These will apply even to library code that creates its own paths to access with NIO.
     *
     * @param maxReopens If the channel errors out, how many times it will attempt to
     *                   re-initiate the connection.
     * @param requesterProject Project to bill when accessing "requester pays" buckets. If unset,
     *                         these buckets cannot be accessed.
     */
    public static void setGlobalNIODefaultOptions(int maxReopens, String requesterProject) {
        CloudStorageFileSystemProvider.setDefaultCloudStorageConfiguration(getCloudStorageConfiguration(maxReopens, requesterProject));
        CloudStorageFileSystemProvider.setStorageOptions(setGenerousTimeouts(StorageOptions.newBuilder()).build());
        setGlobalHttpNioDefaultOptions(maxReopens);
    }

    /**
     * Set the options for http-nio to have generous timeouts and defaults
     */
    public static void setGlobalHttpNioDefaultOptions(final int maxRetries) {
        final HttpFileSystemProviderSettings.RetrySettings retrySettings = new HttpFileSystemProviderSettings.RetrySettings(
                maxRetries,
                RetryHandler.DEFAULT_RETRYABLE_HTTP_CODES,
                RetryHandler.DEFAULT_RETRYABLE_EXCEPTIONS,
                RetryHandler.DEFALT_RETRYABLE_MESSAGES,
                e -> false);
        final HttpFileSystemProviderSettings settings = new HttpFileSystemProviderSettings(
                java.time.Duration.ofMillis(MAX_TIMEOUT),
                HttpFileSystemProviderSettings.DEFAULT_SETTINGS.redirect(),
                retrySettings);

        HttpFileSystemProvider.setSettings(settings);
    }

    /**
     * String -> Path. This *should* not be necessary (use Paths.get(URI.create(...)) instead) , but it currently is
     * on Spark because using the fat, shaded jar breaks the registration of the GCS FilesystemProvider.
     * To transform other types of string URLs into Paths, use IOUtils.getPath instead.
     */
    public static java.nio.file.Path getPathOnGcs(String gcsUrl) {
        // use a split limit of -1 to preserve empty split tokens, especially trailing slashes on directory names
        final String[] split = gcsUrl.split("/", -1);
        final String BUCKET = split[2];
        final String pathWithoutBucket = String.join("/", Arrays.copyOfRange(split, 3, split.length));
        return CloudStorageFileSystem.forBucket(BUCKET).getPath(pathWithoutBucket);
    }

    /**
     * The config we want to use.
     *
     * @param maxReopens If the GCS bucket channel errors out, how many times it will attempt to
     *                   re-initiate the connection.
     * @param requesterProject Project to bill when accessing "requester pays" buckets. If unset,
     *                         these buckets cannot be accessed.
     *
     **/
    public static CloudStorageConfiguration getCloudStorageConfiguration(int maxReopens, String requesterProject) {
        Builder builder = CloudStorageConfiguration.builder()
            // if the channel errors out, re-open up to this many times
            .maxChannelReopens(maxReopens);
        if (!Strings.isNullOrEmpty(requesterProject)) {
            // enable requester pays and indicate who pays
            builder = builder.autoDetectRequesterPays(true).userProject(requesterProject);
        }

        //this causes the gcs filesystem to treat files that end in a / as a directory
        //true is the default but this protects against future changes in behavior
        builder.usePseudoDirectories(true);
        return builder.build();
    }

    private static StorageOptions.Builder setGenerousTimeouts(StorageOptions.Builder builder) {
        return builder
            .setTransportOptions(HttpTransportOptions.newBuilder()
                .setConnectTimeout(MAX_TIMEOUT)
                .setReadTimeout(MAX_TIMEOUT)
                .build())
            .setRetrySettings(RetrySettings.newBuilder()
                .setMaxAttempts(15)
                .setMaxRetryDelay(Duration.ofMillis(256_000L))
                .setTotalTimeout(Duration.ofMillis(4000_000L))
                .setInitialRetryDelay(Duration.ofMillis(1000L))
                .setRetryDelayMultiplier(2.0)
                .setInitialRpcTimeout(Duration.ofMillis(180_000L))
                .setRpcTimeoutMultiplier(1.0)
                .setMaxRpcTimeout(Duration.ofMillis(180_000L))
                .build());
    }

    /**
     * Get an authenticated GCS-backed NIO FileSystem object representing the selected projected and bucket.
     * Credentials are found automatically when running on Compute/App engine, logged into gcloud, or
     * if the GOOGLE_APPLICATION_CREDENTIALS env. variable is set. In that case leave credentials null.
     * Otherwise, you must pass the contents of the service account credentials file.
     * See https://github.com/GoogleCloudPlatform/gcloud-java#authentication
     *
     * Note that most of the time it's enough to just open a file via
     * Files.newInputStream(Paths.get(URI.create( path ))).
     **/
    public static java.nio.file.FileSystem getAuthenticatedGcs(String projectId, String bucket, byte[] credentials) throws IOException {
        StorageOptions.Builder builder = StorageOptions.newBuilder()
                .setProjectId(projectId);
        if (null != credentials) {
            builder = builder.setCredentials(GoogleCredentials.fromStream(new ByteArrayInputStream(credentials)));
        }
        // generous timeouts, to avoid tests failing when not warranted.
        StorageOptions storageOptions = setGenerousTimeouts(builder).build();

        // 2. Create GCS filesystem object with those credentials
        return CloudStorageFileSystem.forBucket(bucket, CloudStorageConfiguration.DEFAULT, storageOptions);
    }

    /**
     * Wrap a SeekableByteChannel with a prefetcher.
     * @param bufferSizeMB buffer size in mb which the prefetcher should fetch ahead.
     * @param channel a channel that needs prefetching
     */
    public static SeekableByteChannel addPrefetcher(final int bufferSizeMB, final SeekableByteChannel channel) {
        try {
            return SeekableByteChannelPrefetcher.addPrefetcher(bufferSizeMB, channel);
        } catch (final IOException ex) {
            throw new GATKException("Unable to initialize the prefetcher: " + ex);
        }
    }

    /**
     * Creates a wrapping function which adds a prefetcher if the buffer size is > 0 if it's <= 0 then this wrapper returns the
     * original channel.
     * @param cloudPrefetchBuffer the prefetcher buffer size in MB
     */
    public static Function<SeekableByteChannel, SeekableByteChannel> getPrefetchingWrapper(final int cloudPrefetchBuffer) {
        return cloudPrefetchBuffer > 0 ? rawChannel -> addPrefetcher(cloudPrefetchBuffer, rawChannel) : Utils.identityFunction();
    }

    /**
     * Take a GCS path and return a signed url to the same resource which allows unauthenticated users to access the file.
     * @param path String representing a GCS path
     * @param hoursToLive how long in hours the url will remain valid
     * @return A signed url which provides access to the bucket location over http allowing unauthenticated users to access it
     */
    public static String createSignedUrlToGcsObject(String path, final long hoursToLive) {
        final Storage storage = StorageOptions.getDefaultInstance().getService();
        final BlobInfo info = BlobInfo.newBuilder(getBucket(path), getPathWithoutBucket(path)).build();
        final URL signedUrl = storage.signUrl(info, hoursToLive, TimeUnit.HOURS, Storage.SignUrlOption.httpMethod(HttpMethod.GET));
        return signedUrl.toString();
    }

    /**
     * Convert a GCS bucket location into the equivalent public http url.  This doesn't do any validation checking
     * to be sure that the location actually exists or is accessible.  It's just a string -> string conversion
     * @param path String representing the gs:// path to an object in a public bucket
     * @return  String representing the https:// path to the same object
     */
    public static String bucketPathToPublicHttpUrl(String path){
        return String.format("https://storage.googleapis.com/%s/%s", getBucket(path), getPathWithoutBucket(path));
    }
}
