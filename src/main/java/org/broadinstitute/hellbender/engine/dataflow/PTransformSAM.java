package org.broadinstitute.hellbender.engine.dataflow;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.exceptions.GATKException;


public abstract class PTransformSAM extends PTransform<PCollection<Read>,PCollection<String>> {
    private String headerString;
    public String getHeaderString(){
        if (headerString == null){
            throw new GATKException("You must call setHeaderString before calling getHeaderString");
        }
        return headerString;
    }

    public void setHeaderString(String headerString) {
        this.headerString = headerString;
    }
}
