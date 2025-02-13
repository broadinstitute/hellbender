package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.tribble.Feature;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.sv.SVFeaturesHeader;

import java.util.*;

/**
 * A MultiFeatureWalkerBase is a base class for tools that consume multiple sources of Features.
 *
 * There are two subclasses that may be convenient:
 * The MergingMultiFeatureWalker presents one feature at a time, in sorted order, by merging
 * multiple sources of features into a single stream.
 * The TabularMultiFeatureWalker presents one feature from each of the sources simultaneously.  It
 * requires that each source have the same number of features (e.g., if each source is conceptually
 * a column in a 2-D matrix of features).
 *
 * If neither of the existing subclasses work for you, you can extend this walker by implementing the
 * abstract {@link #traverse()} method in a class that declares a collection of FeatureInputs as an argument.
 * You may optionally implement {@link #onTraversalStart()}, {@link #onTraversalSuccess()},
 * and/or {@link #closeTool()}.
 */
public abstract class MultiFeatureWalkerBase extends WalkerBase {

    private SAMSequenceDictionary dictionary;
    private final Set<String> samples = new TreeSet<>();

    @Override
    public boolean requiresFeatures(){
        return true;
    }

    @Override
    public String getProgressMeterRecordLabel() { return "features"; }

    /**
     * Operations performed just prior to the start of traversal.
     */
    @Override
    public final void onStartup() {
        super.onStartup();
        setDictionaryAndSamples();
    }

    /**
     * Get the dictionary we settled on
     */
    public SAMSequenceDictionary getDictionary() { return dictionary; }

    /**
     * Get the list of sample names we accumulated
     */
    public Set<String> getSampleNames() { return Collections.unmodifiableSet(samples); }

    /**
     * Choose the most comprehensive dictionary available (see betterDictionary method below),
     * and concatenate the sample names available from each feature input.
     * Each feature input may have its own dictionary, and the user can specify an additional master
     * dictionary, reference dictionary, and reads dictionary.  This method makes certain that all
     * of these dictionaries are consistent with regard to contig name and order.  It's OK if one
     * dictionary is a subset of another:  we'll choose the most comprehensive dictionary.
     * (Can't use the getBestAvailableSequenceDictionary method -- it throws if there are multiple
     * dictionaries available.)
     */
    private void setDictionaryAndSamples() {
        DictSource dictSource = new DictSource(getMasterSequenceDictionary(),
                StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME);
        if ( hasReference() ) {
            final DictSource refDictSource = new DictSource(reference.getSequenceDictionary(),
                    StandardArgumentDefinitions.REFERENCE_LONG_NAME);
            dictSource = betterDictionary(refDictSource, dictSource);
        }
        if ( hasReads() ) {
            final DictSource readsDictSource = new DictSource(reads.getSequenceDictionary(), "read-source");
            dictSource = betterDictionary(readsDictSource, dictSource);
        }
        for ( final FeatureInput<? extends Feature> input : features.getAllInputs() ) {
            final Object header = features.getHeader(input);
            if ( header instanceof SVFeaturesHeader ) {
                final SVFeaturesHeader svFeaturesHeader = (SVFeaturesHeader)header;
                final DictSource featureDictSource = new DictSource(svFeaturesHeader.getDictionary(),
                        input.getName());
                dictSource = betterDictionary(featureDictSource, dictSource);
                final List<String> sampleNames = svFeaturesHeader.getSampleNames();
                if ( sampleNames != null ) {
                    samples.addAll(svFeaturesHeader.getSampleNames());
                }
            } else if (header instanceof VCFHeader ) {
                final VCFHeader vcfHeader = (VCFHeader)header;
                final DictSource featureDictSource = new DictSource(vcfHeader.getSequenceDictionary(),
                        input.getName());
                dictSource = betterDictionary(featureDictSource, dictSource);
                samples.addAll(vcfHeader.getSampleNamesInOrder());
            }
        }
        if ( dictSource.getDictionary() == null ) {
            throw new UserException("No dictionary found.  Provide one as --" +
                    StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME + " or --" +
                    StandardArgumentDefinitions.REFERENCE_LONG_NAME + ".");
        }
        dictionary = dictSource.getDictionary();
    }

    /**
     * Makes sure that the two dictionaries are consistent with regard to contig names and order.
     * Returns the more comprehensive (larger) dictionary if they're consistent.
     */
    private static DictSource betterDictionary( final DictSource newDict,
                                                final DictSource curDict ) {
        if ( curDict.getDictionary() == null ) return newDict;
        if ( newDict.getDictionary() == null ) return curDict;
        final DictSource smallDict;
        final DictSource largeDict;
        if ( newDict.getDictionary().size() <= curDict.getDictionary().size() ) {
            smallDict = newDict;
            largeDict = curDict;
        } else {
            smallDict = curDict;
            largeDict = newDict;
        }
        int lastIdx = -1;
        final SAMSequenceDictionary largeDictionary = largeDict.getDictionary();
        for ( final SAMSequenceRecord rec : smallDict.getDictionary().getSequences() ) {
            final int newIdx = largeDictionary.getSequenceIndex(rec.getContig());
            if ( newIdx == -1 ) {
                throw new UserException("Contig " + rec.getContig() + " in the dictionary read from " +
                        smallDict.getSource() + " does not appear in the larger dictionary read from " +
                        largeDict.getSource());
            }
            if ( newIdx <= lastIdx ) {
                final String prevContig = largeDictionary.getSequence(lastIdx).getContig();
                throw new UserException("Contigs out of order: Contig " + rec.getContig() +
                        " comes before contig " + prevContig + " in the dictionary read from " +
                        largeDict.getSource() + ", but follows it in the dictionary read from " +
                        smallDict.getSource());
            }
            lastIdx = newIdx;
        }
        return largeDict;
    }

    public static final class DictSource {
        private final SAMSequenceDictionary dictionary;
        private final String source;

        public DictSource( final SAMSequenceDictionary dictionary, final String source ) {
            this.dictionary = dictionary;
            this.source = source;
        }

        public SAMSequenceDictionary getDictionary() { return dictionary; }
        public String getSource() { return source; }
    }
}
