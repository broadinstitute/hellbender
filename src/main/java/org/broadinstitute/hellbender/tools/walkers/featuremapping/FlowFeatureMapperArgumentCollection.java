package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Hidden;

import java.io.Serializable;
import java.util.LinkedList;
import java.util.List;

/**
 * Set of arguments for the {@link FlowFeatureMapper}
 */
public class FlowFeatureMapperArgumentCollection implements Serializable{
    private static final long serialVersionUID = 1L;

    enum MappingFeatureEnum {
        SNV
    };

    /**
     *  kind of feature we are mapping (looking for)
     **/
    @Argument(fullName = "mapping-feature", doc = "Kind of feaure being mapped", optional = true)
    public MappingFeatureEnum mappingFeature = MappingFeatureEnum.SNV;

    /**
     *  maximum value for delta in score
     **/
    @Argument(fullName = "limit-score", doc = "Limit value for score", optional = true)
    public double limitScore = Double.NaN;

    /**
     *  attributes to copy from bam
     **/
    @Argument(fullName = "copy-attr", doc = "attributes to copy from bam. format <name>,<type>,<desc>. types: Integer, Float, String, Character, Flag", optional = true)
    public List<String> copyAttr = new LinkedList<>();

    /**
     *  prefix to add to attributes to copy from bam
     **/
    @Argument(fullName = "copy-attr-prefix", doc = "prefix to add to attributes to copy from bam", optional = true)
    public String copyAttrPrefix = "";

    /**
     *  number of bases that need to be identical before the snv
     **/
    @Argument(fullName = "snv-identical-bases", doc = "number of bases that need to be identical before the snv", optional = true)
    public int snvIdenticalBases = 1;

    /**
     *  number of bases that need to be identical after the snv
     **/
    @Argument(fullName = "snv-identical-bases-after", doc = "number of bases that need to be identical after the snv. 0 means same as number of bases before", optional = true)
    public int snvIdenticalBasesAfter = 0;

    /**
     *  threshold of score delta to for emitting (will be emitted if lower)
     **/
    @Argument(fullName = "max-score", doc = "threshold of score delta to for emitting (will be emitted if lower)", optional = true)
    public double maxScore = Double.POSITIVE_INFINITY;

    /**
     *  minimal threshold of score delta to for emitting (will be emitted if higher)
     **/
    @Argument(fullName = "min-score", doc = "minimal threshold of score delta to for emitting (will be emitted if higher)", optional = true)
    public double minScore = Double.NEGATIVE_INFINITY;

    /**
     *  exclude NaN score records?
     **/
    @Argument(fullName = "exclude-nan-scores", doc = "exclude nan scores", optional = true)
    public boolean excludeNaNScores = false;

    /**
     *  include duplicate read?
     **/
    @Argument(fullName = "include-dup-reads", doc = "include duplicate reads", optional = true)
    public boolean includeDupReads = false;

    /**
     *  keep negatives?
     **/
    @Argument(fullName = "keep-negatives", doc = "keep nevative scores?", optional = true)
    public boolean keepNegatives = false;

    /**
     *  keep supplementary alignments?
     **/
    @Argument(fullName = "keep-supplementary-alignments", doc = "keep supplementary alignments ?", optional = true)
    public boolean keepSupplementaryAlignments = false;

    /**
     *  debug negatives?
     **/
    @Hidden
    @Argument(fullName = "debug-negatives", doc = "debug negative scores?", optional = true)
    public boolean debugNegatives = false;

    /**
     *  debug read names?
     **/
    @Hidden
    @Argument(fullName = "debug-read-name", doc = "debug specific reads?", optional = true)
    public List<String> debugReadName = null;

    /**
     *  surrounding media quality (size) - if not specified, this feature is off
     **/
    @Hidden
    @Argument(fullName = "surrounding-median-quality-size", doc = "number of bases around the feature to calculate surrounding median quality", optional = true)
    public Integer surroundingMediaQualitySize = null;

    /**
     *  surrounding mean quality (size) - if not specified, this feature is off
     **/
    @Hidden
    @Argument(fullName = "surrounding-mean-quality-size", doc = "number of bases around the feature to calculate surrounding mean quality", optional = true)
    public Integer surroundingMeanQualitySize = null;

    /**
     *  validation mode - if not specified, this feature is off
     **/
    @Hidden
    @Argument(fullName = "report-all-alts", doc = "In this mode (aka validation mode), every base of every read in the input CRAM and interval is reported, and an X_SCORE value is calculated for all 3 possible alts", optional = true)
    public boolean reportAllAlts = false;

    /**
     *  adjacent-ref-diff mode - if not specified, this feature is off
     **/
    @Hidden
    @Argument(fullName = "tag-bases-with-adjacent-ref-diff", doc = "In this mode bases that have an adjacent difference from the reference on the same read are not discarded, and tagged with X_ADJACENT_REF_DIFFm", optional = true)
    public boolean tagBasesWithAdjacentRefDiff = false;
}
