How to enable s3a:// file accesses
===

## Why do we need s3a://?

Genome datasets are often too big to migrate.
The data migration causes redundant data copies and requirements for additional disk space.
Genome datasets can be accumulated in any place such as AWS S3.
However, current GATK provides direct accesses only for gs:// files.

We have several choices to directly access S3, but I believe [s3a://](https://hadoop.apache.org/docs/current/hadoop-aws/tools/hadoop-aws/index.html) perform well for GATK including Spark-based tools.
s3a:// enables us to directly & efficiently access files on S3-compatible object storage such as AWS S3, IBM COS, and Ceph Object Storage.

## How do we enable s3a://?

Example test scripts

```bash
$ ./gradelew localJar
$ cd scripts/s3a
$ ./setup.sh

$ ./test-local.sh
or
$ ./test-spark.sh
```

Detailed steps

1. Download [hadoop-2.8.2](https://archive.apache.org/dist/hadoop/core/hadoop-2.8.2/hadoop-2.8.2.tar.gz) and [jsr203-s3a-0.0.2.jar](https://repo1.maven.org/maven2/net/fnothaft/jsr203-s3a/0.0.2/jsr203-s3a-0.0.2.jar)
2. Add class libraries for `hadoop-core`, `hadoop-aws`, `aws-java-sdk`, `net.fnothaft:jsr203-s3a`, and gatk-local to Java CLASSPATH.
3. Add `core-site.xml` under a directory in CLASSPATH with properties such as `fs.s3a.access_key` and `fs.s3a.secret_key` ([documentation](https://hadoop.apache.org/docs/current/hadoop-aws/tools/hadoop-aws/index.html#Authenticating_with_S3)).

The above `test-local.sh` uses anonymous accesses to public example datasets in AWS S3.
If you use non-AWS S3 storage like IBM COS, you need to add property `fs.s3a.endpoint` to specify the endpoint URL.

Note that we need to carefully exclude libraries conflicting with existing libraries packaged into `gatk-*.jar`.
For example, the above `test-local.sh` excludes `apache-commons` from CLASSPATH.

The methodology here is leveraging the fact that GATK and other relevant genomics toolkits use java.nio.file.Filesystem API to read files.
In other words, we can also add supports for other kinds of file URI in the same manner here.
