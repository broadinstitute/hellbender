package org.broadinstitute.hellbender.tools.walkers.conversion;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import org.broadinstitute.barclay.argparser.Argument;

import java.util.Comparator;
import java.util.Map;

public class CompareGtfInfo implements Comparator<Map.Entry<String, GtfInfo>>{

    @Argument()//TODO: figure this out
    SAMSequenceDictionary dictionary;

    CompareGtfInfo(SAMSequenceDictionary dictionary){
        this.dictionary = dictionary;
    }

    @Override
    public int compare(Map.Entry<String, GtfInfo> e1, Map.Entry<String, GtfInfo> e2) { //compare two entries of a map where key = geneId or transcriptId and value = gtfInfo object
        Interval e1Interval = e1.getValue().getInterval();
        Interval e2Interval = e2.getValue().getInterval();
        //compare by contig then start
        try {
            return Comparator.comparingInt((Interval interval) -> dictionary.getSequence(interval.getContig()).getSequenceIndex()).thenComparingInt(Interval::getStart).compare(e1Interval, e2Interval);
        }
        catch(NullPointerException e){
            throw new NullPointerException("GTF and Dictionary are different versions");
        }
    }
}
