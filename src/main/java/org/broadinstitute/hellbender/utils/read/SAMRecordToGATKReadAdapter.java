package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.*;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.Serializable;
import java.lang.reflect.Array;
import java.nio.charset.Charset;
import java.util.*;

/**
 * Implementation of the {@link GATKRead} interface for the {@link SAMRecord} class.
 *
 * This adapter wraps a {@link SAMRecord} without making a copy, so construction is cheap,
 * but care must be exercised if the underlying read has been exposed somewhere before
 * wrapping.
 */
public class SAMRecordToGATKReadAdapter implements GATKRead, Serializable {
    private static final long serialVersionUID = 1L;

    private static final Charset DEFAULT_CHARSET = Charset.forName("UTF-8");

    private final SAMRecord samRecord;

    private transient Integer cachedSoftStart = null;
    private transient Integer cachedSoftEnd = null;
    private transient Integer cachedAdaptorBoundary = null;
    private transient Integer cachedCigarLength = null;

    private void clearCachedValues() {
        cachedSoftStart = null;
        cachedSoftEnd = null;
        cachedAdaptorBoundary = null;
        cachedCigarLength = null;
    }

    public int getFlags() {
        return samRecord.getFlags();
    }

    public SAMRecordToGATKReadAdapter( final SAMRecord samRecord ) {
        this.samRecord = samRecord;
    }

    /**
     * Produces a SAMRecordToGATKReadAdapter wrapping the provided SAMRecord,
     * and nulls out the header in the encapsulated read. This is useful for
     * Spark tools, in order to avoid serializing the SAMFileHeader for each
     * record.
     *
     * @param samRecord read to adapt (header will be stripped)
     * @return SAMRecordToGATKReadAdapter wrapping the headerless samRecord
     */
    public static SAMRecordToGATKReadAdapter headerlessReadAdapter( final SAMRecord samRecord ) {
        samRecord.setHeaderStrict(null);
        return new SAMRecordToGATKReadAdapter(samRecord);
    }

    @Override
    public String getName() {
        return samRecord.getReadName();
    }

    @Override
    public void setName( final String name ) {
        clearCachedValues();
        samRecord.setReadName(name);
    }

    @Override
    public String getContig() {
        if ( isUnmapped() ) {
            return null;
        }

        // Guaranteed not to be null or SAMRecord.NO_ALIGNMENT_REFERENCE_NAME due to the isUnmapped() check above
        return samRecord.getReferenceName();
    }

    @Override
    public int getStart() {
        if ( isUnmapped() ) {
            return ReadConstants.UNSET_POSITION;
        }

        // Guaranteed not to be SAMRecord.NO_ALIGNMENT_START due to the isUnmapped() check above
        return samRecord.getAlignmentStart();
    }

    @Override
    public int getEnd() {
        if ( isUnmapped() ) {
            return ReadConstants.UNSET_POSITION;
        }

        // Guaranteed not to be SAMRecord.NO_ALIGNMENT_START due to the isUnmapped() check above
        return samRecord.getAlignmentEnd();
    }

    @Override
    public void setPosition( final String contig, final int start ) {
        assertPositionIsValid(contig, start);

        clearCachedValues();
        samRecord.setReferenceName(contig);
        samRecord.setAlignmentStart(start);
        samRecord.setReadUnmappedFlag(false);
    }

    @Override
    public void setPosition( final Locatable locatable ) {
        if ( locatable == null ) {
            throw new IllegalArgumentException("Cannot set read position to null");
        }

        setPosition(locatable.getContig(), locatable.getStart());
    }

    @Override
    public String getAssignedContig() {
        return samRecord.getReferenceName();
    }

    @Override
    public int getAssignedStart() {
        return samRecord.getAlignmentStart();
    }

    @Override
    public int getUnclippedStart() {
        if ( isUnmapped() ) {
            return ReadConstants.UNSET_POSITION;
        }

        return samRecord.getUnclippedStart();
    }

    @Override
    public int getUnclippedEnd() {
        if ( isUnmapped() ) {
            return ReadConstants.UNSET_POSITION;
        }

        return samRecord.getUnclippedEnd();
    }

    @Override
    public int getSoftStart() {
        if ( cachedSoftStart == null ) {
            cachedSoftStart = ReadUtils.getSoftStart(this);
        }

        return cachedSoftStart;
    }

    @Override
    public int getSoftEnd() {
        if ( cachedSoftEnd == null ) {
            cachedSoftEnd = ReadUtils.getSoftEnd(this);
        }

        return cachedSoftEnd;
    }

    @Override
    public int getAdaptorBoundary() {
        if ( cachedAdaptorBoundary == null ) {
            cachedAdaptorBoundary = ReadUtils.getAdaptorBoundary(this);
        }
        return cachedAdaptorBoundary;
    }

    @Override
    public String getMateContig() {
        if ( mateIsUnmapped() ) {
            return null;
        }

        return samRecord.getMateReferenceName();
    }

    @Override
    public int getMateStart() {
        if ( mateIsUnmapped() ) {
            return ReadConstants.UNSET_POSITION;
        }

        return samRecord.getMateAlignmentStart();
    }

    @Override
    public void setMatePosition( final String contig, final int start ) {
        assertPositionIsValid(contig, start);

        clearCachedValues();

        // Calling this method has the additional effect of marking the read as paired
        setIsPaired(true);

        samRecord.setMateReferenceName(contig);
        samRecord.setMateAlignmentStart(start);
        samRecord.setMateUnmappedFlag(false);
    }

    private static void assertPositionIsValid(final String contig, final int start) {
        if ( contig == null || contig.equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME) || start < 1 ) {
            throw new IllegalArgumentException("contig must be non-null and not equal to " + SAMRecord.NO_ALIGNMENT_REFERENCE_NAME + ", and start must be >= 1"
            + "\nContig: " + contig + "\nStart: " + start);
        }
    }

    @Override
    public void setMatePosition( final Locatable locatable ) {
        Utils.nonNull(locatable, "Cannot set mate position to null");

        setMatePosition(locatable.getContig(), locatable.getStart());
    }

    @Override
    public int getFragmentLength() {
        return samRecord.getInferredInsertSize();
    }

    @Override
    public void setFragmentLength( final int fragmentLength ) {
        clearCachedValues();

        // May be negative if mate maps to lower position than read
        samRecord.setInferredInsertSize(fragmentLength);
    }

    @Override
    public int getMappingQuality() {
        return samRecord.getMappingQuality() != SAMRecord.NO_MAPPING_QUALITY ? samRecord.getMappingQuality() : ReadConstants.NO_MAPPING_QUALITY;
    }

    @Override
    public void setMappingQuality( final int mappingQuality ) {
        if ( mappingQuality < 0 || mappingQuality > 255 ) {
            throw new IllegalArgumentException("mapping quality must be >= 0 and <= 255 \nMappingQuality: " + mappingQuality);
        }

        clearCachedValues();
        samRecord.setMappingQuality(mappingQuality);
    }

    @Override
    public byte[] getBases() {
        final byte[] bases = samRecord.getReadBases();

        // Make a defensive copy to protect against direct modification of the returned array
        return bases != null && bases.length > 0 ? Arrays.copyOf(bases, bases.length) : ArrayUtils.EMPTY_BYTE_ARRAY;
    }

    @Override
    public byte[] getBasesNoCopy() {
        final byte[] bases = samRecord.getReadBases();
        return bases != null ? bases : ArrayUtils.EMPTY_BYTE_ARRAY;
    }

    //Overridden default method to avoid a call to getBases which makes a copy of data
    //Bounds checking is the caller's responsibility, as it's too expensive in this hotspot method
    @Override
    public byte getBase(final int i){
        final byte[] bases = samRecord.getReadBases();
        return bases[i];
    }

    @Override
    public int getLength() {
        final byte[] bases = samRecord.getReadBases();
        return bases == null ? 0 : bases.length;
    }

    @Override
    public void setBases( final byte[] bases ) {
        clearCachedValues();
        samRecord.setReadBases(bases);
    }

    @Override
    public byte[] getBaseQualities() {
        final byte[] baseQualities = samRecord.getBaseQualities();

        // Make a defensive copy to protect against direct modification of the returned array
        return baseQualities != null && baseQualities.length > 0 ? Arrays.copyOf(baseQualities, baseQualities.length) : ArrayUtils.EMPTY_BYTE_ARRAY;
    }

    @Override
    public byte[] getBaseQualitiesNoCopy() {
        final byte[] baseQualities = samRecord.getBaseQualities();
        return baseQualities != null ? baseQualities : ArrayUtils.EMPTY_BYTE_ARRAY;
    }

    @Override
    public int getBaseQualityCount(){
        final byte[] baseQualities = samRecord.getBaseQualities();
        return baseQualities == null ? 0 : baseQualities.length;
    }

    //Overridden default method to avoid a call to getBaseQualities which makes a copy of data
    //Bounds checking is the caller's responsibility, as it's too expensive in this hotspot method
    @Override
    public byte getBaseQuality(final int i)  {
        final byte[] baseQualities = samRecord.getBaseQualities();
        return baseQualities[i];
    }

    @Override
    public void setBaseQualities( final byte[] baseQualities ) {
        if ( baseQualities != null ) {
            for ( byte b : baseQualities ) {
                if ( b < 0 ) {
                    throw new IllegalArgumentException("Base quality score " + b + " is invalid");
                }
            }
        }

        clearCachedValues();
        samRecord.setBaseQualities(baseQualities);
    }

    @Override
    public Cigar getCigar() {
        // Make a defensive copy before returning to guard against modification of the return value,
        // since Cigar is a mutable type:
        return samRecord.getCigar() != null ? new Cigar(samRecord.getCigar().getCigarElements()) : new Cigar();
    }

    /**
     * This implementation does not make a new Cigar object but instead provides
     * an unmodifiable view of the underlying list of CigarElements.
     * This is done to reduce the amount of object allocation.
     */
    @Override
    public List<CigarElement> getCigarElements(){
        //Cigar.getCigarElements returns an unmodifiable list so we don't wrap it again
        return samRecord.getCigar() == null ? Collections.emptyList() : samRecord.getCigar().getCigarElements();
    }

    /**
     * This implementation avoids the creation of the unmodifiable view of the underlying list of CigarElements
     * and simply retrieves the element that is requested.
     *
     * Bounds checking is the caller's responsibility, as it's too expensive in this hotspot method.
     */
    @Override
    public CigarElement getCigarElement(final int index) {
        return samRecord.getCigar().getCigarElement(index);
    }

    /**
     * This implementation saves time by not making an unmodifiable view of the list of
     * elements but returns the length of the list directly (or 0 if there's no cigar).
     */
    @Override
    public int numCigarElements() {
        // It's surprising and bizarre, but profiling reveals that caching the cigar length
        // actually helps performance in some cases (eg., the HaplotypeCaller)
        if ( cachedCigarLength == null ) {
            cachedCigarLength = samRecord.getCigar() == null ? 0 : samRecord.getCigarLength();
        }
        return cachedCigarLength;
    }

    @Override
    public void setCigar( final Cigar cigar ) {
        clearCachedValues();
        samRecord.setCigar(cigar);
    }

    @Override
    public void setCigar( final String cigarString ) {
        clearCachedValues();
        samRecord.setCigarString(cigarString);
    }

    @Override
    public String getReadGroup() {
        // May return null
        return (String)samRecord.getAttribute(SAMTag.RG.getBinaryTag());
    }

    @Override
    public void setReadGroup( final String readGroupID ) {
        clearCachedValues();
        samRecord.setAttribute(SAMTag.RG.name(), readGroupID);
    }

    @Override
    public boolean isPaired() {
        return samRecord.getReadPairedFlag();
    }

    @Override
    public void setIsPaired( final boolean isPaired ) {
        clearCachedValues();

        samRecord.setReadPairedFlag(isPaired);
        if ( ! isPaired ) {
            samRecord.setProperPairFlag(false);
        }
    }

    @Override
    public boolean isProperlyPaired() {
        return isPaired() && samRecord.getProperPairFlag();
    }

    @Override
    public void setIsProperlyPaired( final boolean isProperlyPaired ) {
        clearCachedValues();

        if ( isProperlyPaired ) {
            setIsPaired(true);
        }

        samRecord.setProperPairFlag(isProperlyPaired);
    }

    @Override
    public boolean isUnmapped() {
        return samRecord.getReadUnmappedFlag() ||
               samRecord.getReferenceName() == null || samRecord.getReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME) ||
               samRecord.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START;
    }

    @Override
    public boolean isUnplaced() {
        return samRecord.getReferenceName() == null || samRecord.getReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME) ||
                samRecord.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START;
    }

    @Override
    public void setIsUnmapped() {
        clearCachedValues();

        samRecord.setReadUnmappedFlag(true);
    }


    @Override
    public void setIsUnplaced() {
        clearCachedValues();

        samRecord.setReadUnmappedFlag(true);
        samRecord.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
        samRecord.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
        samRecord.setMappingQuality(SAMRecord.NO_MAPPING_QUALITY);
    }

    @Override
    public boolean mateIsUnmapped() {
        Utils.validate(isPaired(), "Cannot get mate information for an unpaired read");

        return samRecord.getMateUnmappedFlag() ||
               samRecord.getMateReferenceName() == null || samRecord.getMateReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME) ||
               samRecord.getMateAlignmentStart() == SAMRecord.NO_ALIGNMENT_START;
    }

    @Override
    public boolean mateIsUnplaced() {
        Utils.validate(isPaired(), "Cannot get mate information for an unpaired read");

        return samRecord.getMateReferenceName() == null || samRecord.getMateReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME) ||
                samRecord.getMateAlignmentStart() == SAMRecord.NO_ALIGNMENT_START;
    }

    @Override
    public void setMateIsUnmapped() {
        clearCachedValues();

        // Calling this method has the side effect of marking the read as paired.
        setIsPaired(true);

        samRecord.setMateUnmappedFlag(true);
    }

    @Override
    public void setMateIsUnplaced() {
        clearCachedValues();

        // Calling this method has the side effect of marking the read as paired.
        setIsPaired(true);

        samRecord.setMateUnmappedFlag(true);
        samRecord.setMateNegativeStrandFlag(false);
        samRecord.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
        samRecord.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
    }

    @Override
    public boolean isReverseStrand() {
        return samRecord.getReadNegativeStrandFlag();
    }

    @Override
    public void setIsReverseStrand( final boolean isReverseStrand ) {
        clearCachedValues();

        samRecord.setReadNegativeStrandFlag(isReverseStrand);
    }

    @Override
    public boolean mateIsReverseStrand() {
        Utils.validate(isPaired(), "Cannot get mate information for an unpaired read");

        return samRecord.getMateNegativeStrandFlag();
    }

    @Override
    public void setMateIsReverseStrand( final boolean mateIsReverseStrand ) {
        clearCachedValues();

        // Calling this method has the side effect of marking the read as paired.
        setIsPaired(true);

        samRecord.setMateNegativeStrandFlag(mateIsReverseStrand);
    }

    @Override
    public boolean isFirstOfPair() {
        return isPaired() && samRecord.getFirstOfPairFlag();
    }

    @Override
    public void setIsFirstOfPair() {
        clearCachedValues();

        // Calling this method has the side effect of marking the read as paired.
        setIsPaired(true);

        samRecord.setFirstOfPairFlag(true);
        samRecord.setSecondOfPairFlag(false);
    }

    @Override
    public boolean isSecondOfPair() {
        return isPaired() && samRecord.getSecondOfPairFlag();
    }

    @Override
    public void setIsSecondOfPair() {
        clearCachedValues();

        // Calling this method has the side effect of marking the read as paired.
        setIsPaired(true);

        samRecord.setSecondOfPairFlag(true);
        samRecord.setFirstOfPairFlag(false);
    }

    @Override
    public boolean isSecondaryAlignment() {
        return samRecord.isSecondaryAlignment();
    }

    @Override
    public void setIsSecondaryAlignment( final boolean isSecondaryAlignment ) {
        clearCachedValues();

        samRecord.setSecondaryAlignment(isSecondaryAlignment);
    }

    @Override
    public boolean isSupplementaryAlignment() {
        return samRecord.getSupplementaryAlignmentFlag();
    }

    @Override
    public void setIsSupplementaryAlignment( final boolean isSupplementaryAlignment ) {
        clearCachedValues();

        samRecord.setSupplementaryAlignmentFlag(isSupplementaryAlignment);
    }

    @Override
    public boolean failsVendorQualityCheck() {
        return samRecord.getReadFailsVendorQualityCheckFlag();
    }

    @Override
    public void setFailsVendorQualityCheck( final boolean failsVendorQualityCheck ) {
        clearCachedValues();

        samRecord.setReadFailsVendorQualityCheckFlag(failsVendorQualityCheck);
    }

    @Override
    public boolean isDuplicate() {
        return samRecord.getDuplicateReadFlag();
    }

    @Override
    public void setIsDuplicate( final boolean isDuplicate ) {
        clearCachedValues();

        samRecord.setDuplicateReadFlag(isDuplicate);
    }

    @Override
    public boolean hasAttribute( final String attributeName ) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);
        return samRecord.getAttribute(attributeName) != null;
    }

    @Override
    public Integer getAttributeAsInteger( final String attributeName ) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);
        final Object attributeValue = samRecord.getAttribute(attributeName);

        if ( attributeValue == null ) {
            return null;
        }
        else if ( attributeValue instanceof Integer ) {
            return (Integer)attributeValue;
        }
        else {
            try {
                return Integer.parseInt(attributeValue.toString());
            }
            catch ( NumberFormatException e ) {
                throw new GATKException.ReadAttributeTypeMismatch(attributeName, "integer", e);
            }
        }
    }

    @Override
    public Float getAttributeAsFloat( final String attributeName ) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);
        final Object attributeValue = samRecord.getAttribute(attributeName);

        if ( attributeValue == null ) {
            return null;
        }
        else if ( attributeValue instanceof Float ) {
            return (Float)attributeValue;
        }
        else {
            try {
                return Float.parseFloat(attributeValue.toString());
            }
            catch ( NumberFormatException e ) {
                throw new GATKException.ReadAttributeTypeMismatch(attributeName, "integer", e);
            }
        }
    }

    @Override
    public String getAttributeAsString( final String attributeName ) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);
        final Object attributeValue = samRecord.getAttribute(attributeName);
        if ( attributeValue instanceof byte[]) {
            // in case that the attribute is a byte[] array, the toString method will format it as name@hashCode
            // for a good representation of the byte[] as String, it encodes the bytes with the default charset (UTF-8)
            final byte[] val = (byte[]) attributeValue;
            return (val.length == 0) ? "" : new String(val, DEFAULT_CHARSET);
        }
        // otherwise, just use the toString() method unless it is null
        return attributeValue != null ? attributeValue.toString() : null;
    }

    @Override
    public byte[] getAttributeAsByteArray( final String attributeName ) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);
        final Object attributeValue = samRecord.getAttribute(attributeName);

        if ( attributeValue == null ) {
            return null;
        }
        else if ( attributeValue instanceof byte[] ) {
            // In the case where the attribute value is already a byte[], make a defensive
            // copy before returning to guard against modification of the return value.
            final byte[] ret = (byte[])attributeValue;
            return Arrays.copyOf(ret, ret.length);
        }
        else if ( attributeValue instanceof String ) {
            return ((String)attributeValue).getBytes(DEFAULT_CHARSET);
        }
        else {
            throw new GATKException.ReadAttributeTypeMismatch(attributeName, "byte array");
        }
    }

    @Override
    public void setAttribute( final String attributeName, final Integer attributeValue ) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);

        clearCachedValues();
        samRecord.setAttribute(attributeName, attributeValue);
    }

    @Override
    public void setAttribute( final String attributeName, final String attributeValue ) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);

        clearCachedValues();
        samRecord.setAttribute(attributeName, attributeValue);
    }

    @Override
    public void setAttribute( final String attributeName, final byte[] attributeValue ) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);

        clearCachedValues();
        samRecord.setAttribute(attributeName, attributeValue);
    }

    @Override
    public void clearAttribute( final String attributeName ) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);

        clearCachedValues();
        samRecord.setAttribute(attributeName, null);
    }

    @Override
    public void clearAttributes() {
        clearCachedValues();

        samRecord.clearAttributes();
    }

    @Override
    public void clearTransientAttribute( final String attributeName ) {
        clearCachedValues();
        samRecord.removeTransientAttribute(attributeName);
    }

    @Override
    public GATKRead copy() {
        // Produces a shallow but "safe to use" copy.
        return new SAMRecordToGATKReadAdapter(ReadUtils.cloneSAMRecord(samRecord));
    }

    @Override
    public GATKRead deepCopy() {
        // Produces a true deep copy.
        return new SAMRecordToGATKReadAdapter(samRecord.deepCopy());
    }

    @Override
    public String getSAMString() {
        return samRecord.getSAMString();
    }

    @Override
    public void reverseComplement() {
        clearCachedValues();

        samRecord.reverseComplement(true);
    }

    /**
     * Clip attributes that change after a hard clipping operation
     *
     * The method clips string or array attributes that are part of a pre-set list
     * originating from FlowBasedRead (a white-list), that are at the same length
     * as the sequence.
     */
    @Override
    public void hardClipAttributes(final int newStart, final int newLength, final int originalLength)
    {
        // this is currently applicable only for flow based read groups
        if ( samRecord.getReadGroup() == null || samRecord.getReadGroup().getFlowOrder() == null ) {
            return;
        }

        // loop over whitelist
        for ( final String tag : FlowBasedRead.HARD_CLIPPED_TAGS ) {

            // attribute present?
            if ( !samRecord.hasAttribute(tag) ) {
                continue;
            }

            // ignore if not array
            Object  value = samRecord.getAttribute(tag);

            if ( !(value.getClass().isArray() || (value instanceof String) )) {
                continue;
            }

            // check if length is same as original sequence
            if ( value.getClass().isArray() && (Array.getLength(value) != originalLength) ) {
                continue;
            }

            if ((value instanceof String )&&(((String) value).length()!=originalLength)) {
                continue;
            }

            // clip
            if ( value instanceof byte[] ) {
                samRecord.setAttribute(tag, Arrays.copyOfRange((byte[])value, newStart, newStart + newLength));
            } else if ( value instanceof short[] ) {
                samRecord.setAttribute(tag, Arrays.copyOfRange((short[])value, newStart, newStart + newLength));
            } else if ( value instanceof int[] ) {
                samRecord.setAttribute(tag, Arrays.copyOfRange((int[])value, newStart, newStart + newLength));
            } else if ( value instanceof float[] ) {
                samRecord.setAttribute(tag, Arrays.copyOfRange((float[]) value, newStart, newStart + newLength));
            } else if (value instanceof String){
                samRecord.setAttribute(tag, ((String) value).substring(newStart, newStart + newLength));
            } else {
                throw new GATKException("unknown array type on attribute value: " + value.getClass().getName());
            }
        }
    }

    @Override
    public SAMRecord convertToSAMRecord( final SAMFileHeader header ) {
        samRecord.setHeaderStrict(header);
        return samRecord;
    }

    public SAMRecord getEncapsulatedSamRecord() {
        return samRecord;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        SAMRecordToGATKReadAdapter that = (SAMRecordToGATKReadAdapter) o;

        return Objects.equals(samRecord, that.samRecord);

    }

    @Override
    public int hashCode() {
        return Objects.hashCode(samRecord);
    }

    @Override
    public String toString() {
        return commonToString();
    }

    public boolean hasHeader() {
        return samRecord.getHeader() != null;
    }

    public void setHeader(SAMFileHeader header) {
        clearCachedValues();

        samRecord.setHeaderStrict(header);
    }

    /**
     * This is used to access the transient attribute store in the underlying SAMRecord.
     *
     * NOTE: This is an advanced use case for SAMRecord and you should probably use setAttribute() instead
     * @param key key whose value is to be retrived
     */
    public Object getTransientAttribute(Object key) {
        return samRecord.getTransientAttribute(key);
    }

    /**
     * This is used to access the transient attribute store in the underlying SAMRecord. This is used to store temporary
     * attributes that will not be serialized and that do not trigger the SAMRecord to parse the attributes if they are not needed.
     *
     * NOTE: This is an advanced use case for SAMRecord and you should probably use setAttribute() instead
     * @param key key under which the value will be stored
     * @param value value to be keyed
     */
    public void setTransientAttribute(Object key, Object value) {
        samRecord.setTransientAttribute(key, value);
    }
}
