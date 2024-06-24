package org.broadinstitute.hellbender.engine.filters;

import com.google.common.annotations.VisibleForTesting;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;
import java.util.List;
import java.util.function.Predicate;

/**
 * Filters which operate on {@link GATKRead} should subclass this by overriding {@link #test(GATKRead)}
 *
 * ReadFilter implements Predicate and Serializable.  It provides a default implementation of apply
 * based on the implementing class's implementation of test().
 *
 * To be accessible from the command line, subclasses must have a zero-arg constructor and contain
 * (optional) ArgumentCollection or Argument annotated fields for filter-specific arguments. Arguments
 * that are optional=true should have an initial value since the command line parser will not require
 * the user to provide one; arguments that are optional=false should not have an initial value (and
 * used boxed types with null default/initial values) since otherwise they will be seen by the command
 * line parser as having been provided, and the user will not be notified that the value must be provided
 * on the command line).
 */
public abstract class ReadFilter implements Predicate<GATKRead>, Serializable {

    protected SAMFileHeader samHeader = null;

    public void setHeader(SAMFileHeader samHeader) { this.samHeader = samHeader; }

    @VisibleForTesting
    public static class ReadFilterNegate extends ReadFilter {
        private static final long serialVersionUID = 1L;

        private final ReadFilter delegate;

        protected ReadFilterNegate(ReadFilter delegate) {
            Utils.nonNull(delegate, "Delegate filter cannot be null");
            this.delegate = delegate;
        }

        @Override
        public boolean test( GATKRead read ) {
            return !delegate.test(read);
        }

        @Override
        public String getName() {
            return "NOT " + delegate.getName();
        }
    }

    protected abstract static class ReadFilterBinOp extends ReadFilter {
        private static final long serialVersionUID = 1L;

        final protected ReadFilter lhs;
        final protected ReadFilter rhs;

        public ReadFilterBinOp(final ReadFilter lhs, final ReadFilter rhs) {
            Utils.nonNull(lhs, "ReadFilterBinOp lhs filter cannot be null");
            Utils.nonNull(rhs, "ReadFilterBinOp rhs filter cannot be null");
            this.lhs = lhs;
            this.rhs = rhs;
        }

        @Override
        public void setHeader(SAMFileHeader samHeader) {
            super.setHeader(samHeader);
            if ( lhs != null ) {
                lhs.setHeader(samHeader);
            }
            if ( rhs != null ) {
                rhs.setHeader(samHeader);
            }
        }

    }

    @VisibleForTesting
    protected static class ReadFilterAnd extends ReadFilterBinOp {
        private static final long serialVersionUID = 1L;

        public ReadFilterAnd(ReadFilter lhs, ReadFilter rhs) { super(lhs, rhs); }

        @Override
        public boolean test( GATKRead read ) { return lhs.test(read) && rhs.test(read); }

        @Override
        public String getName() {
            return lhs.getName() + " AND " + rhs.getName();
        }
    }

    private static class ReadFilterOr extends ReadFilterBinOp {
        private static final long serialVersionUID = 1L;

        public ReadFilterOr(ReadFilter lhs, ReadFilter rhs) { super(lhs, rhs); }

        @Override
        public boolean test( GATKRead read ) { return lhs.test(read) || rhs.test(read);}

        @Override
        public String getName() {
            return lhs.getName() + " OR " + rhs.getName();
        }
    }

    // It turns out, this is necessary. Please don't remove it.
    // Without this line, we see the following error:
    // java.io.InvalidClassException: org.broadinstitute.hellbender.engine.filters.ReadFilter; local class incompatible:
    // stream classdesc serialVersionUID = -5040289903122017748, local class serialVersionUID = 6814309376393671214
    private static final long serialVersionUID = 1L;

    /**
     * Return a composite (and) {@code ReadFilter} constructed from a list of {@code ReadFilter}. Each
     * filter in the list is first initialized with the {@code SAMFileHeader} param. The resulting filter
     * honors the order of the input list and tests the filter conditions in the same order as the iteration
     * order of the input list.
     * @param readFilters If null or empty, the ALLOW_ALL_READS read filter will be returned
     * @param samHeader {@code SAMFileHeader} used to initialize each filter. May not be null
     * @return Composite ReadFilter
     */
    public static ReadFilter fromList(final List<ReadFilter> readFilters, final SAMFileHeader samHeader) {
        Utils.nonNull(samHeader, "SAMFileHeader must not be null");
        if (readFilters == null || readFilters.isEmpty()) {
            return ReadFilterLibrary.ALLOW_ALL_READS;
        }
        readFilters.forEach(f -> f.setHeader(samHeader));
        ReadFilter compositeFilter = readFilters.get(0);
        for (int i = 1; i < readFilters.size(); i++) {
            compositeFilter = compositeFilter.and(readFilters.get(i));
        }
        return compositeFilter;
     }

    //HACK: These methods are a hack to get the type system to accept compositions of ReadFilters.
    /**
     * Specialization of {@link #and(Predicate)} so that ReadFilters anded with other ReadFilters produce a ReadFilter
     */
    public ReadFilter and( ReadFilter other ) {
        Utils.nonNull(other);
        return new ReadFilterAnd(this, other);
    }

    /**
     * Specialization of {@link #or(Predicate)} so that ReadFilters or'ed with other ReadFilters produce a ReadFilter
     */
    public ReadFilter or( ReadFilter other ) {
        Utils.nonNull(other);
        return new ReadFilterOr(this, other);
    }

    /**
     * Specialization of negate so that the resulting object is still a ReadFilter
     */
    @Override
    public ReadFilter negate(){ return new ReadFilterNegate(this); }

    @Override
    public abstract boolean test( GATKRead read );

    /**
     * Retruns the display name for logs for this filter
     */
    String getName() {
        return this.getClass().getSimpleName();
    }
}
