package org.broadinstitute.hellbender.engine.filters;

import com.google.common.annotations.VisibleForTesting;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;
import java.util.Stack;
import java.util.function.Predicate;
import java.util.stream.IntStream;

/**
 * Wrapper/adapter for {@link org.broadinstitute.hellbender.engine.filters.ReadFilter} that counts the number of
 * reads filtered, and provides a filter count summary.
 *
 * For filters with complex predicates that are composed of compound and/or operators, provides counts
 * by level for each component predicate. All counts reflect short-circuited evaluation of and/or operators
 * (ie., not all read filters in a compound predicate will necessarily be evaluated every time). Also note
 * that the count for a compound predicate does not always equal the sum of the counts of it's component
 * predicates, i.e. an "or" filter will report a filter count of 1 in the case where both component predicates
 * are evaluated and both fail, but the the individual component filters will report a count of 1 at the next
 * level.
 */
public class CountingReadFilter extends ReadFilter {

    private static final long serialVersionUID = 1L;

    // Underlying ReadFilter we delegate to if we're wrapping a simple ReadFilter.
    @VisibleForTesting
    protected final ReadFilter delegateFilter;

    // Number of reads filtered by this filter
    protected long filteredCount = 0;

    public CountingReadFilter(final ReadFilter readFilter) {
        Utils.nonNull(readFilter);
        delegateFilter = readFilter;
    }

    /**
     * Return a composite (and) {@code CountingReadFilter} constructed from a list of
     * {@link org.broadinstitute.hellbender.engine.filters.ReadFilter}. Each filter in the list is first
     * initialized with the {@code SAMFileHeader} param. The resulting filter honors the order of the input list
     * and tests the filter conditions in the same order as the iteration order of the input list.
     * @param readFilters If null or empty, the ALLOW_ALL_READS read filter will be returned
     * @param samHeader {@code SAMFileHeader} used to initialize each filter. May not be null
     * @return Composite CountingReadFilter
     */
    public static CountingReadFilter fromList(final List<ReadFilter> readFilters, final SAMFileHeader samHeader) {
        Utils.nonNull(samHeader, "SAMFileHeader must not be null");
        if (readFilters == null || readFilters.isEmpty()) {
            return new CountingReadFilter(ReadFilterLibrary.ALLOW_ALL_READS);
        }
        readFilters.forEach(f -> f.setHeader(samHeader));
        CountingReadFilter compositeFilter = new CountingReadFilter(readFilters.get(0));
        for (int i = 1; i < readFilters.size(); i++) {
            compositeFilter = compositeFilter.and(new CountingReadFilter(readFilters.get(i)));
        }
        return compositeFilter;
    }

    // Used only by the nested CountingBinopReadFilter subclass and its derivatives, which must
    // override the test method with an implementation that does not depend on delegateFilter.
    private CountingReadFilter() {
        delegateFilter = null;
    }

    // Return the number of reads filtered by this filter
    public long getFilteredCount() {
        return filteredCount;
    }

    public void resetFilteredCount() {
        filteredCount = 0;
    }

    public String getName() {return delegateFilter.getName();}

    // Returns a summary line with filter counts organized by level
    public String getSummaryLine() {return getSummaryLineForLevel(0);}

    protected String getSummaryLineForLevel(final int indentLevel) {
        return getIndentString(indentLevel) + Long.toString(filteredCount) + " read(s) filtered by: " + getName() + " \n";
    }

    protected String getIndentString(final int indentLevel) {
        final StringBuilder bldr = new StringBuilder();
        IntStream.range(0, indentLevel).forEach(i -> bldr.append("  "));
        return bldr.toString();
    }

    /**
     * Specialization of {@link #and(Predicate)} so that CountingReadFilter and'ed with other CountingReadFilter produce a CountingReadFilter
     */
    //@Override
    public CountingReadFilter and(final CountingReadFilter other) {
        Utils.nonNull(other);
        return new CountingAndReadFilter(this, other);
    }

    /**
     * Specialization of {@link #or(Predicate)} so that CountingReadFilter ored with other CountingReadFilter produce a CountingReadFilter
     */
    //@Override
    public CountingReadFilter or(final CountingReadFilter other) {
        Utils.nonNull(other);
        return new CountingOrReadFilter(this, other);
    }

    /**
     * Specialization of negate so that the resulting object is still a CountingReadFilter
     */
    @Override
    public CountingReadFilter negate() {
        return new CountingNegateReadFilter(this);
    }

    @Override
    public boolean test(final GATKRead read) {
        final boolean accept = delegateFilter.test(read);
        if (!accept) {
            filteredCount++;
        }
        return accept;
    }

    private static class CountingNegateReadFilter extends CountingReadFilter {
        private static final long serialVersionUID = 1L;

        CountingReadFilter delegateCountingFilter;

        public CountingNegateReadFilter(CountingReadFilter delegate) {
            this.delegateCountingFilter = delegate;
        }

        @Override
        public boolean test(GATKRead read) {
            final boolean accept = !delegateCountingFilter.test(read);
            if (!accept) {
                filteredCount++;
            }
            return accept;
        }

        @Override
        public String getName() {
            return "NOT " + delegateCountingFilter.getName();
        }
    }

    /**
     * Private class for Counting binary operator (and/or) filters; these keep track of how many reads are filtered at
     * each level of filter nesting.
     *
     * Subclasses must override the test method.
     */
    private static abstract class CountingBinopReadFilter extends CountingReadFilter {

        private static final long serialVersionUID = 1L;

        protected final CountingReadFilter lhs;
        protected final CountingReadFilter rhs;

        public CountingBinopReadFilter(final CountingReadFilter lhs, final CountingReadFilter rhs) {
            Utils.nonNull(lhs);
            Utils.nonNull(rhs);
            this.lhs = lhs;
            this.rhs = rhs;
        }

        @Override
        protected String getSummaryLineForLevel(final int indentLevel) {
            final String indent = getIndentString(indentLevel);
            if (0 == filteredCount) {
                return indent +  "0 read(s) filtered by: " + getName() + " \n";
            }
            else {
                return indent + Long.toString(filteredCount) + " read(s) filtered by: " + getName() + "\n"
                        + (lhs.getFilteredCount() > 0 ? lhs.getSummaryLineForLevel(indentLevel + 1) : "")
                        + (rhs.getFilteredCount() > 0 ? rhs.getSummaryLineForLevel(indentLevel + 1) : "");
            }
        }

        @Override
        public void resetFilteredCount() {
            super.resetFilteredCount();
            this.lhs.resetFilteredCount();
            this.rhs.resetFilteredCount();
        }

        @Override
        public abstract String getName();
    }

    /**
     * Private class for Counting AND filters
     */
    @VisibleForTesting
    protected static final class CountingAndReadFilter extends CountingBinopReadFilter {

        private static final long serialVersionUID = 1L;
        protected long totalCount = 0;

        private CountingAndReadFilter(final CountingReadFilter lhs, final CountingReadFilter rhs) {
            super(lhs, rhs);
        }

        @Override
        public boolean test(final GATKRead read) {
            final boolean accept = lhs.test(read) && rhs.test(read);
            if (!accept) {
                filteredCount++;
            }
            totalCount++;
            return accept;
        }

        @Override
        protected String getSummaryLineForLevel(final int indentLevel) {
            if(indentLevel == 0) {
                String simplifiedSummary = getSummaryLineForLevelAllAndsSimplified();
                if(!simplifiedSummary.isEmpty()) {
                    return simplifiedSummary;
                }
            }
            return super.getSummaryLineForLevel(indentLevel);
        }

        /**
         * @return simplified summary line if this filter consists only of filters joined by AND, otherwise empty String
         */
        private String getSummaryLineForLevelAllAndsSimplified() {
            StringBuilder summaryLine = new StringBuilder();

            Stack<CountingReadFilter> unread = new Stack<>();
            unread.push(this.rhs);
            unread.push(this.lhs);
            CountingReadFilter curFilter;

            while(!unread.empty()) {
                curFilter = unread.pop();
                if (curFilter instanceof CountingAndReadFilter) {
                    unread.push(((CountingAndReadFilter) curFilter).rhs);
                    unread.push(((CountingAndReadFilter) curFilter).lhs);
                } else if (curFilter instanceof CountingBinopReadFilter || curFilter instanceof CountingNegateReadFilter) {
                    return "";
                } else {
                    summaryLine.append(curFilter.getSummaryLineForLevel(0));
                }
            }

            summaryLine.append(this.getFilteredCount() + " total reads filtered out of "+ this.totalCount+ " reads processed");

            return summaryLine.toString();
        }

        @Override
        public String getName() {
            return "(" + lhs.getName() + " AND " + rhs.getName() + ")";
        }
    }

    /**
     * Private class for Counting OR filters
     */
    private static final class CountingOrReadFilter extends CountingBinopReadFilter {

        private static final long serialVersionUID = 1L;

        private CountingOrReadFilter(final CountingReadFilter lhs, final CountingReadFilter rhs) {
            super(lhs, rhs);
        }

        @Override
        public boolean test(final GATKRead read) {
            final  boolean accept = lhs.test(read) || rhs.test(read);
            if (!accept) {
                filteredCount++;
            }
            return accept;
        }

        @Override
        public String getName() {
            return "(" + lhs.getName() + " OR " + rhs.getName() + ")";
        }
    }
}
