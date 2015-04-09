package org.broadinstitute.hellbender.tools.dataflow.pipelines;

import com.google.api.services.genomics.model.Read;
import com.google.api.services.storage.Storage;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.io.TextIO;
import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.dataflow.sdk.transforms.*;
import com.google.cloud.dataflow.sdk.values.PBegin;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.readers.bam.BAMIO;
import com.google.cloud.genomics.dataflow.readers.bam.ReadBAMTransform;
import com.google.cloud.genomics.dataflow.readers.bam.ReadConverter;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import com.google.cloud.genomics.dataflow.utils.GenomicsOptions;
import com.google.cloud.genomics.utils.Contig;
import com.google.cloud.genomics.utils.GenomicsFactory;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTool;
import org.broadinstitute.hellbender.engine.dataflow.SAMSerializableFunction;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.engine.dataflow.PTransformSAM;
import org.broadinstitute.hellbender.transformers.ReadTransformer;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.logging.Logger;
import java.util.stream.Collectors;

public abstract class DataflowReadsPipeline extends DataflowTool{


    public static final String GCS_PREFIX = "gs://";
    @Argument
    private String outputFile;

    @Argument
    private File clientSecret = new File("client_secret.json");

    @Argument
    String bam;



    abstract protected PTransformSAM getTool();

    protected ImmutableList<ReadFilter> getReadFilters(){
        return ImmutableList.of();
    }

    protected ImmutableList<ReadTransformer> getReadTransformers(){
        return ImmutableList.of();
    }

    private GenomicsFactory.OfflineAuth getAuth(GCSOptions options){
        try {
            return GCSOptions.Methods.createGCSAuth(options);
        } catch (IOException e) {
            throw new GATKException("Couldn't create a dataflow auth object.", e);
        }
    }

    private String getHeaderString(GCSOptions options, String bamPath) {
        if(isCloudStorageUrl(bamPath)) {
            try {
                Storage.Objects storageClient = GCSOptions.Methods.createStorageClient(options, getAuth(options));
                final SamReader reader = BAMIO.openBAM(storageClient, bamPath);
                return reader.getFileHeader().getTextHeader();
            } catch (IOException e) {
                throw new GATKException("Failed to read bam header from " + bamPath + ".", e);
            }
        } else {
            return SamReaderFactory.makeDefault().getFileHeader(new File(bamPath)).getTextHeader();
        }
    }


    private SerializableFunction<Read,Boolean> wrapFilter(ReadFilter filter, String headerString){
        return new SAMSerializableFunction<>(headerString, filter);
    }


    @Override
    final protected void setupPipeline(Pipeline pipeline) {
        GCSOptions ops = PipelineOptionsFactory.fromArgs(new String[]{"--genomicsSecretsFile=" + clientSecret.getAbsolutePath()}).as(GCSOptions.class);
        GenomicsOptions.Methods.validateOptions(ops);


        GenomicsFactory.OfflineAuth auth = getAuth(ops);
        String headerString = getHeaderString(ops, bam);


        Iterable<Contig> contigs = intervals.stream()
                .map(i -> new Contig(i.getContig(), i.getStart(), i.getEnd()))
                .collect(Collectors.toList());

        PCollection<Read> preads;
        if( isCloudStorageUrl(bam)) {
            preads = ReadBAMTransform.getReadsFromBAMFilesSharded(pipeline, auth, contigs, ImmutableList.of(bam));
        } else {
            preads = pipeline.apply(Create.of(ImmutableList.of(bam)))
                    .apply(ParDo.of(new DoFn<String, Read>() {
                        @Override
                        public void processElement(ProcessContext c) throws Exception {
                            ReadsDataSource sams = new ReadsDataSource(new File(c.element()));
                            sams.setIntervalsForTraversal(intervals);
                            for (SAMRecord sam : sams) {
                                c.output(ReadConverter.makeRead(sam));
                            }
                        }
                    }));
        }

        for (ReadFilter filter : getReadFilters()) {
            preads = preads.apply(Filter.by(wrapFilter(filter, headerString)));
        }

        PTransformSAM f = getTool();
        f.setHeaderString(headerString);

        PCollection<String> pstrings = preads.apply(f);
        pstrings.apply(TextIO.Write.to(outputFile));
    }

    private boolean isCloudStorageUrl(String path) {
        return path.startsWith(GCS_PREFIX);
    }
}
