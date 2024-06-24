package org.broadinstitute.hellbender.tools.gvs.extract;

import htsjdk.samtools.util.Locatable;
import org.apache.avro.generic.GenericRecord;
import org.apache.commons.lang3.builder.ReflectionToStringBuilder;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;

public class ExtractCohortFilterRecord implements Locatable {

    private final Long location;
    private final String contig;
    private final Integer start;
    private final Integer end;

    private final Double score;
    private final Double vqScore;
    private final String yng;
    private final String refAllele;
    private final String altAllele;

    public ExtractCohortFilterRecord(GenericRecord genericRecord, final String vqScoreFieldName, final String scoreFieldName) {
        this.location = Long.parseLong(genericRecord.get(SchemaUtils.LOCATION_FIELD_NAME).toString());
        this.contig = SchemaUtils.decodeContig(location);
        this.start = SchemaUtils.decodePosition(location);
        this.end = start;

        this.score = (scoreFieldName != null) ? Double.parseDouble(genericRecord.get(SchemaUtils.SCORE).toString()) : null;
        this.vqScore = Double.parseDouble(genericRecord.get(vqScoreFieldName).toString());
        this.yng = genericRecord.get(SchemaUtils.YNG_STATUS).toString();

        this.refAllele = genericRecord.get(SchemaUtils.REF_ALLELE_FIELD_NAME).toString();
        this.altAllele = genericRecord.get(SchemaUtils.ALT_ALLELE_FIELD_NAME).toString();
    }

    @Override
    public String getContig() { return this.contig; }

    @Override
    public int getStart() { return this.start; }

    @Override
    public int getEnd() { return this.end; }

    public long getLocation() { return this.location; }

    public Double getScore() { return score; }

    public double getVqScore() { return vqScore; }

    public String getYng() { return this.yng; }

    public String getRefAllele() { return this.refAllele; }

    public String getAltAllele() { return this.altAllele; }

    public String toString() {
        return ReflectionToStringBuilder.toString(this);
    }
}
