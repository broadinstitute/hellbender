package org.broadinstitute.hellbender.engine.dataflow.datasources;

import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.util.LinkedHashMap;
import java.util.Map;

public class ReferenceFileSource implements Serializable {
    private static final long serialVersionUID = 1L;

    private final String referencePath;

    public ReferenceFileSource(final String referencePath) {
        this.referencePath = referencePath;
    }

    public ReferenceBases getReferenceBases( final SimpleInterval interval ) throws IOException {
        ReferenceSequenceFile referenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(new File(referencePath));
        ReferenceSequence sequence = referenceSequenceFile.getSubsequenceAt(interval.getContig(), interval.getStart(), interval.getEnd());
        return new ReferenceBases(sequence.getBases(), interval);
    }

    public Map<String, ReferenceBases> getAllReferenceBases() throws IOException {
        ReferenceSequenceFile referenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(new File(referencePath));
        Map<String, ReferenceBases> bases = new LinkedHashMap<>();
        ReferenceSequence seq;
        while ((seq = referenceSequenceFile.nextSequence()) != null) {
            String name = seq.getName();
            bases.put(name, new ReferenceBases(seq.getBases(), new SimpleInterval(name, 1, seq.length())));
        }
        return bases;
    }


}
