package org.broadinstitute.hellbender.tools.gvs.extract;

import htsjdk.samtools.util.Locatable;
import org.apache.avro.generic.GenericRecord;
import org.apache.avro.util.Utf8;
import org.apache.commons.lang3.builder.ReflectionToStringBuilder;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;

public class ReferenceRecord implements Locatable, Comparable<ReferenceRecord> {

    private final long location;
    private final long sampleId;
    private final String contig;
    private final Integer start;
    private final Integer end;
    private final Long endLocation;

    private final String state;

    public ReferenceRecord(GenericRecord genericRecord) {
        this.location = (Long) genericRecord.get(SchemaUtils.LOCATION_FIELD_NAME);
        this.sampleId = (Long) genericRecord.get(SchemaUtils.SAMPLE_ID_FIELD_NAME);
        this.contig = SchemaUtils.decodeContig(location);
        this.start = SchemaUtils.decodePosition(location);

        int length = Math.toIntExact((Long) genericRecord.get("length"));
        this.end = this.start + length - 1;
        this.endLocation = this.location + length - 1;
        this.state = ((Utf8) genericRecord.get(SchemaUtils.STATE_FIELD_NAME)).toString();
    }

    public ReferenceRecord(long location, long sampleId, int length, String state) {
        this.location = location;
        this.sampleId = sampleId;
        this.contig = SchemaUtils.decodeContig(location);
        this.start = SchemaUtils.decodePosition(location);

        this.end = this.start + length - 1;
        this.endLocation = this.location + length - 1;

        this.state = state;
    }

    @Override
    public String getContig() { return this.contig; }

    @Override
    public int getStart() { return this.start; }

    @Override
    public int getEnd() { return this.end; }

    public long getLocation() { return this.location; }
    public long getEndLocation() { return this.endLocation; }

    public long getSampleId() { return this.sampleId; }

    public String getState() { return this.state; }

    public String toString() {
        return ReflectionToStringBuilder.toString(this);
    }

    @Override
    public int compareTo(ReferenceRecord o) {
        final long firstPosition = this.location;
        final long secondPosition = o.location;

        final int result = Long.compare(firstPosition, secondPosition);
        if (result != 0) {
            return result;
        } else {
            final long firstSample = this.location;
            final long secondSample = o.location;
            return Long.compare(firstSample, secondSample);
        }
    }
}
