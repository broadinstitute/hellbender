package org.broadinstitute.hellbender.tools.gvs.extract;

import htsjdk.samtools.util.Locatable;
import org.apache.avro.generic.GenericRecord;
import org.apache.commons.lang3.builder.ReflectionToStringBuilder;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;

import java.util.Objects;

public class ExtractCohortRecord implements Locatable {

    private final long location;
    private final long sampleId;
    private final String contig;
    private final int start;
    private final int end;

    private final String state;
    private final String refAllele;
    private final String altAllele;
    private final String callGT;
    private final String callAD;
    private final String callGQ;
    private final String callRGQ;
    private final String qualapprox;
    private final String asQualapprox;
    private final String callPL;
    private final String callPGT;

    private final String callPID;

    private final String callPS;

    // COHORT_FIELDS
//    public static final ImmutableSet<String> COHORT_FIELDS = ImmutableSet.of(
//            SchemaUtils.LOCATION_FIELD_NAME,
//            SchemaUtils.SAMPLE_ID_FIELD_NAME,
//            SchemaUtils.STATE_FIELD_NAME,
//            SchemaUtils.REF_ALLELE_FIELD_NAME,
//            SchemaUtils.ALT_ALLELE_FIELD_NAME,
//            SchemaUtils.CALL_GT,
//            SchemaUtils.CALL_GQ,
//            SchemaUtils.CALL_RGQ,
//            SchemaUtils.QUALapprox,
//            SchemaUtils.AS_QUALapprox,
//            SchemaUtils.CALL_PL,
//            SchemaUtils.CALL_PGT,
//            SchemaUtils.CALL_PID,
//            SchemaUtils.CALL_PS);//, AS_VarDP);

    // If the requested key is missing in an Avro record:
    //
    // Avro 1.11 throws: https://github.com/apache/avro/blob/release-1.11.0/lang/java/avro/src/main/java/org/apache/avro/generic/GenericData.java#L267-L269
    // Avro 1.8 returns null: https://github.com/apache/avro/blob/release-1.8.0/lang/java/avro/src/main/java/org/apache/avro/generic/GenericData.java#L208
    //
    // Most of the code here was written for Avro 1.8 behavior; this function provides the equivalent for Avro 1.11.
    private static Object getNoThrow(GenericRecord record, String key) {
        return record.hasField(key) ? record.get(key) : null;
    }

    public ExtractCohortRecord(GenericRecord genericRecord) {
        this.location = (Long) genericRecord.get(SchemaUtils.LOCATION_FIELD_NAME);
        this.sampleId = (Long) genericRecord.get(SchemaUtils.SAMPLE_ID_FIELD_NAME);
        this.contig = SchemaUtils.decodeContig(location);
        this.start = SchemaUtils.decodePosition(location);

        this.end = start;

        // if this record is being constructed from the VET data, we won't have a state so we default it to 'v'
        this.state = Objects.toString(getNoThrow(genericRecord, SchemaUtils.STATE_FIELD_NAME), "v");

        // the rest are nullable
        this.refAllele = Objects.toString(genericRecord.get(SchemaUtils.REF_ALLELE_FIELD_NAME), null);
        this.altAllele = Objects.toString(genericRecord.get(SchemaUtils.ALT_ALLELE_FIELD_NAME), null);
        this.callGT = Objects.toString(genericRecord.get(SchemaUtils.CALL_GT), null);
        this.callAD = Objects.toString(getNoThrow(genericRecord, SchemaUtils.CALL_AD), null);
        this.callGQ = Objects.toString(genericRecord.get(SchemaUtils.CALL_GQ), null);
        this.qualapprox = Objects.toString(genericRecord.get(SchemaUtils.QUALapprox), null);
        this.asQualapprox = Objects.toString(genericRecord.get(SchemaUtils.AS_QUALapprox), null);
        this.callPL = Objects.toString(genericRecord.get(SchemaUtils.CALL_PL), null);

        // Phasing-specific fields
        this.callPGT = Objects.toString(getNoThrow(genericRecord, SchemaUtils.CALL_PGT), null);
        this.callPID = Objects.toString(getNoThrow(genericRecord, SchemaUtils.CALL_PID), null);
        this.callPS  = Objects.toString(getNoThrow(genericRecord, SchemaUtils.CALL_PS), null);

        // to keep callRGQ final...
        String tmpRGQ = Objects.toString(getNoThrow(genericRecord, SchemaUtils.CALL_RGQ), null);

        // if we don't get RGQ from the database, we can calculate it from the PLs
        if (tmpRGQ != null) {
            this.callRGQ = tmpRGQ;
        } else if (this.callPL != null) {
            this.callRGQ = this.callPL.split(",")[0];
        } else {
            this.callRGQ = null;
        }
    }

    // constructor for single base reference
    public ExtractCohortRecord(long location, long sampleId, String state) {
        this(location, sampleId, state, null, null, null, null, null, null, null, null, null, null, null, null );
    }

    public ExtractCohortRecord(long location, long sampleId, String state, String refAllele, String altAllele, String callGT, String callAD, String callGQ, String callRGQ, String qualapprox, String asQualapprox, String callPL, String callPGT, String callPID, String callPS) {
        this.location = location;
        this.sampleId = sampleId;
        this.contig = SchemaUtils.decodeContig(location);
        this.start = SchemaUtils.decodePosition(location);
        this.end = start;
        this.state = state;

        this.refAllele = refAllele;
        this.altAllele = altAllele;
        this.callGT = callGT;
        this.callAD = callAD;
        this.callGQ = callGQ;
        this.callRGQ = callRGQ;
        this.qualapprox = qualapprox;
        this.asQualapprox = asQualapprox;
        this.callPL = callPL;
        this.callPGT = callPGT;
        this.callPID = callPID;
        this.callPS = callPS;
    }

    @Override
    public String getContig() { return this.contig; }

    @Override
    public int getStart() { return this.start; }

    @Override
    public int getEnd() { return this.end; }

    public long getLocation() { return this.location; }

    public Long getSampleId() { return this.sampleId; }

    public String getState() { return this.state; }

    public String getRefAllele() { return this.refAllele; }

    public String getAltAllele() { return this.altAllele; }

    public String getCallGT() { return this.callGT; }

    public String getCallAD() { return this.callAD; }

    public String getCallGQ() { return this.callGQ; }

    public String getCallRGQ() { return this.callRGQ; }

    public String getQUALApprox() { return this.qualapprox; }

    public String getAsQUALApprox() { return this.asQualapprox; }

    public String getCallPL() { return this.callPL; }

    public String getCallPGT() { return callPGT; }

    public String getCallPID() { return callPID; }

    public String getCallPS() { return callPS; }

    public String toString() {
        return ReflectionToStringBuilder.toString(this);
    }
}
