package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.*;

/**
 * A MergingMultiFeatureWalker is a base class for a tool that processes one {@link Feature} at a
 * time in sorted order from multiple sources of Features.
 * The input files for each feature must be sorted by locus.
 *
 * To use this walker you need only implement the abstract
 * {@link #apply(F, Object, ReadsContext, ReferenceContext)} method in a class that declares
 * a collection of FeatureInputs as an argument.
 * You may optionally implement {@link #onTraversalStart()}, {@link #onTraversalSuccess()},
 * and/or {@link #closeTool()}.
 */
public abstract class MergingMultiFeatureWalker<F extends Feature> extends MultiFeatureWalkerBase {
    /**
     * {@inheritDoc}
     *
     * Implementation of Feature-based traversal where the Features from the various sources are
     * presented in dictionary order.
     *
     * NOTE: You should only override {@link #traverse()} if you are writing a new walker base class
     * in the engine package that extends this class. It is not meant to be overridden by tools
     * outside the engine package.
     */
    @Override
    public void traverse() {
        final CountingReadFilter readFilter = makeReadFilter();
        final MergingIterator<F> iterator =
                new MergingIterator<>(getDictionary(), features, userIntervals);
        while ( iterator.hasNext() ) {
            final PQEntry<F> entry = iterator.next();
            final F feature = entry.getFeature();
            final SimpleInterval featureInterval = new SimpleInterval(feature);
            apply(feature,
                    entry.getHeader(),
                    new ReadsContext(reads, featureInterval, readFilter),
                    new ReferenceContext(reference, featureInterval));
            progressMeter.update(feature);
        }
    }

    /**
     * Process an individual feature.
     * In general, subclasses should simply stream their output from apply(), and maintain as little
     * internal state as possible.
     *
     * @param feature Current Feature being processed.
     * @param header Header object for the source from which the feature was drawn (may be null)
     * @param readsContext An object that allows querying for the reads the overlap the feature
     * @param referenceContext An object that allows querying for the reference sequence associated with the feature
     */
    public abstract void apply( final F feature,
                                final Object header,
                                final ReadsContext readsContext,
                                final ReferenceContext referenceContext );

    public static final class MergingIterator<F extends Feature> implements Iterator<PQEntry<F>> {
        final SAMSequenceDictionary dictionary;
        final PriorityQueue<PQEntry<F>> priorityQueue;

        @SuppressWarnings("unchecked")
        public MergingIterator( final SAMSequenceDictionary dictionary,
                                final FeatureManager featureManager,
                                final List<SimpleInterval> intervals ) {
            final Set<FeatureInput<? extends Feature>> inputs = featureManager.getAllInputs();
            this.dictionary = dictionary;
            this.priorityQueue = new PriorityQueue<>(inputs.size());
            for ( final FeatureInput<? extends Feature> input : inputs ) {
                final Iterator<F> iterator =
                        (Iterator<F>)featureManager.getFeatureIterator(input, intervals);
                final Object header = featureManager.getHeader(input);
                addEntry(new PQContext<>(iterator, dictionary, header));
            }
        }

        @Override
        public boolean hasNext() {
            return !priorityQueue.isEmpty();
        }

        @Override
        public PQEntry<F> next() {
            final PQEntry<F> entry = priorityQueue.poll();
            if ( entry == null ) {
                throw new NoSuchElementException("iterator is exhausted");
            }
            final PQEntry<F> newEntry = addEntry(entry.getContext());
            if ( newEntry != null ) {
                if ( newEntry.compareTo(entry) < 0 ) {
                    final Feature feature = newEntry.getFeature();
                    throw new UserException("inputs are not sorted at " +
                                            feature.getContig() + ":" + feature.getStart());
                }
            }
            return entry;
        }

        private PQEntry<F> addEntry( final PQContext<F> context ) {
            final Iterator<F> iterator = context.getIterator();
            if ( !iterator.hasNext() ) {
                return null;
            }
            final PQEntry<F> entry  = new PQEntry<>(context, iterator.next());
            priorityQueue.add(entry);
            return entry;
        }
    }

    public static final class PQContext<F extends Feature> {
        private final Iterator<F> iterator;
        private final SAMSequenceDictionary dictionary;
        private final Object header;

        public PQContext( final Iterator<F> iterator,
                          final SAMSequenceDictionary dictionary,
                          final Object header ) {
            this.iterator = iterator;
            this.dictionary = dictionary;
            this.header = header;
        }

        public Iterator<F> getIterator() { return iterator; }
        public SAMSequenceDictionary getDictionary() { return dictionary; }
        public Object getHeader() { return header; }
    }

    public static final class PQEntry<F extends Feature> implements Comparable<PQEntry<F>> {
        private final PQContext<F> context;
        private final int contigIndex;
        private final F feature;

        public PQEntry( final PQContext<F> context, final F feature ) {
            this.context = context;
            this.feature = feature;
            this.contigIndex = context.getDictionary().getSequenceIndex(feature.getContig());
        }

        public PQContext<F> getContext() { return context; }
        public F getFeature() { return feature; }
        public Object getHeader() { return context.getHeader(); }

        @Override
        public int compareTo( final PQEntry<F> entry ) {
            int result = Integer.compare(contigIndex, entry.contigIndex);
            if ( result == 0 ) {
                result = Integer.compare(feature.getStart(), entry.feature.getStart());
                if ( result == 0 ) {
                    result = Integer.compare(feature.getEnd(), entry.feature.getEnd());
                }
            }
            return result;
        }
    }
}
