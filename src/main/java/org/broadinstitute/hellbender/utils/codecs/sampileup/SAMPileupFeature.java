package org.broadinstitute.hellbender.utils.codecs.sampileup;

import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.Feature;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * A tribble feature representing a SAM pileup.
 *
 * Allows intake of simple mpileups. Simple pileup features will contain only basic information, no reconstructed reads.
 *
 * @author Danil Gomez-Sanchez (magiDGS)
 */
public class SAMPileupFeature implements Feature {

    // genomic location
    private final String contig;
    private final int position;
    // reference base
    private final byte refBase;

    // list of pileup elements
    private final List<SAMPileupElement> pileupElements;

    SAMPileupFeature(final String contig, final int position, final byte refBase, final List<SAMPileupElement> pileupElements) {
        Utils.nonNull(pileupElements);
        this.contig = contig;
        this.position = position;
        this.refBase = refBase;
        this.pileupElements = pileupElements;
    }

    @Override
    @Deprecated
    public String getChr() {
        return getContig();
    }

    @Override
    public String getContig() {
        return contig;
    }

    @Override
    public int getStart() {
        return position;
    }

    @Override
    public int getEnd() {
        return position;
    }

    /**
     * Returns pile of obseved qualities over the genomic location
     *
     * Note: this call costs O(n) and allocates fresh array each time
     */
    public String getQualsString() {
        return SAMUtils.phredToFastq(getBaseQuals());
    }

    /**
     * Returns pile of observed bases over the genomic location.
     *
     * Note: this call costs O(n) and allocates fresh array each time
     */
    public String getBasesString() {
        return StringUtil.bytesToString(getBases());
    }

    /**
     * Returns the reference basse
     */
    public byte getRef() {
        return refBase;
    }

    /**
     * Return the number of observed bases over the genomic location
     */
    public int size() {
        return pileupElements.size();
    }

    /**
     * Format in a samtools-like string.
     * Each line represents a genomic position, consisting of chromosome name, coordinate,
     * reference base, read bases and read qualities
     */
    public String getPileupString() {
        // In the pileup format,
        return String.format("%s %s %c %s %s",
                getContig(), getStart(),    // chromosome name and coordinate
                getRef(),                                                     // reference base
                getBasesString(),
                getQualsString());
    }

    /**
     * Gets the bases in byte array form
     *
     * Note: this call costs O(n) and allocates fresh array each time
     */
    public byte[] getBases() {
        final List<Byte> bases = getBasesStream().collect(Collectors.toList());
        return ArrayUtils.toPrimitive(bases.toArray(new Byte[bases.size()]));
    }

    /**
     * Gets the PHRED base qualities.
     *
     * Note: this call costs O(n) and allocates fresh array each time
     */
    public byte[] getBaseQuals() {
        final List<Byte> quals = getQualsStream().collect(Collectors.toList());
        return ArrayUtils.toPrimitive(quals.toArray(new Byte[quals.size()]));
    }

    /**
     * Get the bases as a stream
     */
    public Stream<Byte> getBasesStream() {
        return pileupElements.stream().map(SAMPileupElement::getBase);
    }

    /**
     * Get the qualities as a stream
     */
    public Stream<Byte> getQualsStream() {
        return pileupElements.stream().map(SAMPileupElement::getBaseQuality);
    }

}
