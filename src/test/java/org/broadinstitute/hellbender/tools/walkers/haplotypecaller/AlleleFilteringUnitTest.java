package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.haplotype.Event;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;

public class AlleleFilteringUnitTest {

    @Test
    public void testNoNeedToFilter(){
        List<Haplotype> haplotypeList = new ArrayList<>();
        final byte[] fullReferenceWithPadding = "CATGCATG".getBytes();
        Haplotype haplotype = new Haplotype(fullReferenceWithPadding, true, 0, TextCigarCodec.decode("8M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 108));
        haplotype.setEventMap(EventMap.fromHaplotype(haplotype, fullReferenceWithPadding, 0));


        haplotypeList.add(haplotype);
        haplotype = new Haplotype("CAGGCATG".getBytes(), false, 0, TextCigarCodec.decode("8M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 108));

        haplotypeList.add(haplotype);
        haplotype.setEventMap(EventMap.fromHaplotype(haplotype, fullReferenceWithPadding, 0));

        AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(haplotypeList);
        SampleList samples = new IndexedSampleList(Arrays.asList("sm1"));

        List<GATKRead> readList = new ArrayList<>(30);
        Map<String, List<GATKRead>>ebs = new HashMap<>();
        ebs.put("sm1", readList);

        for (int i = 0 ; i < 30; i++) {
            readList.add(ArtificialReadUtils.createArtificialRead("20M"));
        }

        double[][] values = {{0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3},
                                {3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0}};

        AlleleLikelihoods<GATKRead, Haplotype> lks = new AlleleLikelihoods<>(samples, haplotypes, ebs);
        LikelihoodMatrix<GATKRead, Haplotype> lkm = lks.sampleMatrix(0);
        for (int i = 0; i < lks.numberOfAlleles(); i++){
            for (int j = 0 ; j < lkm.evidenceCount(); j++) {
                lkm.set(i,j,values[i][j]);
            }
        }

        HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
        HaplotypeCallerGenotypingEngine genotypingEngine = new HaplotypeCallerGenotypingEngine(hcArgs, samples, ! hcArgs.doNotRunPhysicalPhasing, false);


        AlleleFiltering alleleFiltering = new AlleleFilteringHC(hcArgs, null, genotypingEngine);
        AlleleLikelihoods<GATKRead, Haplotype> filtered_lks = alleleFiltering.filterAlleles(lks, 0, new HashSet<>());
        Assert.assertEquals(filtered_lks.alleles(), lks.alleles());
    }
    @Test
    public void testNoNeedToFilterTwoSamples(){
        List<Haplotype> haplotypeList = new ArrayList<>();
        final byte[] fullReferenceWithPadding = "CATGCATG".getBytes();
        Haplotype haplotype = new Haplotype(fullReferenceWithPadding, true, 0, TextCigarCodec.decode("8M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 108));
        haplotype.setEventMap(EventMap.fromHaplotype(haplotype, fullReferenceWithPadding, 0));
        haplotypeList.add(haplotype);

        haplotype = new Haplotype("CAGGCATG".getBytes(), false, 0, TextCigarCodec.decode("8M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 108));
        haplotypeList.add(haplotype);
        haplotype.setEventMap(EventMap.fromHaplotype(haplotype, fullReferenceWithPadding, 0));

        haplotype = new Haplotype("CATTCATG".getBytes(), false, 0, TextCigarCodec.decode("8M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 108));
        haplotypeList.add(haplotype);
        haplotype.setEventMap(EventMap.fromHaplotype(haplotype, fullReferenceWithPadding, 0));

        AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(haplotypeList);
        SampleList samples = new IndexedSampleList(Arrays.asList(new String[]{"sm1", "sm2"}));

        List<GATKRead> readList = new ArrayList<>(30);
        Map<String, List<GATKRead>>ebs = new HashMap<>();
        ebs.put("sm1", readList);

        for (int i = 0 ; i < 30; i++) {
            readList.add(ArtificialReadUtils.createArtificialRead("20M"));
        }

        List<GATKRead> readList2 = new ArrayList<>(30);
        ebs.put("sm2", readList2);
        for (int i = 0 ; i < 30; i++) {
            readList2.add(ArtificialReadUtils.createArtificialRead("20M"));
        }

        double[][] values = {{0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3},
                {3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0},
                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
        };
        AlleleLikelihoods<GATKRead, Haplotype> lks = new AlleleLikelihoods<>(samples, haplotypes, ebs);
        LikelihoodMatrix<GATKRead, Haplotype> lkm = lks.sampleMatrix(0);
        for (int i = 0; i < lks.numberOfAlleles(); i++){
            for (int j = 0 ; j < lkm.evidenceCount(); j++) {
                lkm.set(i,j,values[i][j]);
            }
        }


        double[][] values2 = {{0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3},
                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                {3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0}
        };
        LikelihoodMatrix<GATKRead, Haplotype> lkm2 = lks.sampleMatrix(1);
        for (int i = 0; i < lks.numberOfAlleles(); i++){
            for (int j = 0 ; j < lkm2.evidenceCount(); j++) {
                lkm2.set(i,j,values2[i][j]);
            }
        }

        HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
        HaplotypeCallerGenotypingEngine genotypingEngine = new HaplotypeCallerGenotypingEngine(hcArgs, samples, ! hcArgs.doNotRunPhysicalPhasing, false);


        AlleleFiltering alleleFiltering = new AlleleFilteringHC(hcArgs, null, genotypingEngine);
        AlleleLikelihoods<GATKRead, Haplotype> filtered_lks = alleleFiltering.filterAlleles(lks, 0, new HashSet<>());
        Assert.assertEquals(filtered_lks.alleles(), lks.alleles());
    }



    @Test
    public void testFilterCloseMis(){

        // create haplotypes
        List<Haplotype> haplotypeList = new ArrayList<>();
        final byte[] fullReferenceWithPadding = "CATGCATG".getBytes();
        Haplotype haplotype = new Haplotype(fullReferenceWithPadding, true, 0, TextCigarCodec.decode("8M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 108));
        haplotype.setEventMap(EventMap.fromHaplotype(haplotype, fullReferenceWithPadding, 0));
        haplotypeList.add(haplotype);

        haplotype = new Haplotype("CAGGCATG".getBytes(), false, 0, TextCigarCodec.decode("8M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 108));

        haplotype.setEventMap(EventMap.fromHaplotype(haplotype, fullReferenceWithPadding, 0));
        haplotypeList.add(haplotype);

        haplotype = new Haplotype("CAGTCATG".getBytes(), false, 0, TextCigarCodec.decode("8M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 108));

        haplotype.setEventMap(EventMap.fromHaplotype(haplotype, fullReferenceWithPadding, 0));
        haplotypeList.add(haplotype);

        AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(haplotypeList);
        SampleList samples = new IndexedSampleList(Arrays.asList("sm1"));



        List<GATKRead> readList = new ArrayList<>(30);
        Map<String, List<GATKRead>>ebs = new HashMap<>();
        ebs.put("sm1", readList);

        for (int i = 0 ; i < 30; i++) {
            readList.add(ArtificialReadUtils.createArtificialRead("20M"));
        }

        double[][] values = {{0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3},
                {3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0},
                {2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
        };

        AlleleLikelihoods<GATKRead, Haplotype> lks = new AlleleLikelihoods<>(samples, haplotypes, ebs);
        LikelihoodMatrix<GATKRead, Haplotype> lkm = lks.sampleMatrix(0);
        for (int i = 0; i < lks.numberOfAlleles(); i++){
            for (int j = 0 ; j < lkm.evidenceCount(); j++) {
                lkm.set(i,j,values[i][j]);
            }
        }

        HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
        HaplotypeCallerGenotypingEngine genotypingEngine = new HaplotypeCallerGenotypingEngine(hcArgs, samples, ! hcArgs.doNotRunPhysicalPhasing, false);


        AlleleFiltering alleleFiltering = new AlleleFilteringHC(hcArgs, null, genotypingEngine);
        AlleleLikelihoods<GATKRead, Haplotype> filtered_lks = alleleFiltering.filterAlleles(lks, 0, new HashSet<>());
        Assert.assertEquals(filtered_lks.alleles(), haplotypeList.subList(0,2));
    }

    @Test
    public void testFilterDistantHindel(){

        // create haplotypes
        List<Haplotype> haplotypeList = new ArrayList<>();
        final byte[] fullReferenceWithPadding = "CATGCATG".getBytes();
        Haplotype haplotype = new Haplotype(fullReferenceWithPadding, true, 0, TextCigarCodec.decode("8M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 108));
        haplotype.setEventMap(EventMap.fromHaplotype(haplotype, fullReferenceWithPadding, 0));
        haplotypeList.add(haplotype);

        haplotype = new Haplotype("CAGGCATG".getBytes(), false, 0, TextCigarCodec.decode("8M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 108));

        haplotype.setEventMap(EventMap.fromHaplotype(haplotype, fullReferenceWithPadding, 0));
        haplotypeList.add(haplotype);

        haplotype = new Haplotype("CAGGCATTG".getBytes(), false, 0, TextCigarCodec.decode("7M1I1M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 109));

        haplotype.setEventMap(EventMap.fromHaplotype(haplotype, fullReferenceWithPadding, 0));
        haplotypeList.add(haplotype);

        AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(haplotypeList);
        SampleList samples = new IndexedSampleList(Arrays.asList("sm1"));



        List<GATKRead> readList = new ArrayList<>(30);
        Map<String, List<GATKRead>>ebs = new HashMap<>();
        ebs.put("sm1", readList);

        for (int i = 0 ; i < 30; i++) {
            readList.add(ArtificialReadUtils.createArtificialRead("20M"));
        }

        double[][] values = {{0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3},
                {3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0},
                {2, 0, 1, 0, 2, 0, 2, 0, 2, 0, 1, 0, 4, 0, 2, 0, 2, 0, 2, 0, 2, 0, 4, 0, 2, 0, 2, 0, 2, 0}
        };

        AlleleLikelihoods<GATKRead, Haplotype> lks = new AlleleLikelihoods<>(samples, haplotypes, ebs);
        LikelihoodMatrix<GATKRead, Haplotype> lkm = lks.sampleMatrix(0);
        for (int i = 0; i < lks.numberOfAlleles(); i++){
            for (int j = 0 ; j < lkm.evidenceCount(); j++) {
                lkm.set(i,j,values[i][j]);
            }
        }

        HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
        HaplotypeCallerGenotypingEngine genotypingEngine = new HaplotypeCallerGenotypingEngine(hcArgs, samples, ! hcArgs.doNotRunPhysicalPhasing, false);


        AlleleFiltering alleleFiltering = new AlleleFilteringHC(hcArgs, null, genotypingEngine);
        AlleleLikelihoods<GATKRead, Haplotype> filtered_lks = alleleFiltering.filterAlleles(lks, 0, new HashSet<>());
        Assert.assertEquals(filtered_lks.alleles(), haplotypeList.subList(0,2));
    }

    @Test
    public void testNotFilterDistantM(){

        // create haplotypes
        List<Haplotype> haplotypeList = new ArrayList<>();
        final byte[] fullReferenceWithPadding = "CATGCATG".getBytes();
        Haplotype haplotype = new Haplotype(fullReferenceWithPadding, true, 0, TextCigarCodec.decode("8M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 108));
        haplotype.setEventMap(EventMap.fromHaplotype(haplotype, fullReferenceWithPadding, 0));
        haplotypeList.add(haplotype);

        haplotype = new Haplotype("CAGGCATG".getBytes(), false, 0, TextCigarCodec.decode("8M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 108));

        haplotype.setEventMap(EventMap.fromHaplotype(haplotype, fullReferenceWithPadding, 0));
        haplotypeList.add(haplotype);

        haplotype = new Haplotype("CATGCATC".getBytes(), false, 0, TextCigarCodec.decode("8M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 108));

        haplotype.setEventMap(EventMap.fromHaplotype(haplotype, fullReferenceWithPadding, 0));
        haplotypeList.add(haplotype);

        AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(haplotypeList);
        SampleList samples = new IndexedSampleList(Arrays.asList("sm1"));



        List<GATKRead> readList = new ArrayList<>(30);
        Map<String, List<GATKRead>>ebs = new HashMap<>();
        ebs.put("sm1", readList);

        for (int i = 0 ; i < 30; i++) {
            readList.add(ArtificialReadUtils.createArtificialRead("20M"));
        }

        double[][] values = {{0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3},
                {3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0},
                {3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0}
        };

        AlleleLikelihoods<GATKRead, Haplotype> lks = new AlleleLikelihoods<>(samples, haplotypes, ebs);
        LikelihoodMatrix<GATKRead, Haplotype> lkm = lks.sampleMatrix(0);
        for (int i = 0; i < lks.numberOfAlleles(); i++){
            for (int j = 0 ; j < lkm.evidenceCount(); j++) {
                lkm.set(i,j,values[i][j]);
            }
        }

        HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
        HaplotypeCallerGenotypingEngine genotypingEngine = new HaplotypeCallerGenotypingEngine(hcArgs, samples, ! hcArgs.doNotRunPhysicalPhasing, false);


        AlleleFiltering alleleFiltering = new AlleleFilteringHC(hcArgs, null, genotypingEngine);
        AlleleLikelihoods<GATKRead, Haplotype> filtered_lks = alleleFiltering.filterAlleles(lks, 100, new HashSet<>());
        Assert.assertEquals(filtered_lks.alleles(), haplotypeList.subList(0,1));
    }

    @Test
    public void testNotFilterLoneWeakAllele(){
        //if hcArgs.filterLoneAlleles is false - weak lone alleles should be kept

        // create haplotypes
        List<Haplotype> haplotypeList = new ArrayList<>();
        final byte[] fullReferenceWithPadding = "CATGCATG".getBytes();
        Haplotype haplotype = new Haplotype(fullReferenceWithPadding, true, 0, TextCigarCodec.decode("8M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 108));
        haplotype.setEventMap(EventMap.fromHaplotype(haplotype, fullReferenceWithPadding, 0));
        haplotypeList.add(haplotype);

        haplotype = new Haplotype("CAGGCATG".getBytes(), false, 0, TextCigarCodec.decode("8M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 108));

        haplotype.setEventMap(EventMap.fromHaplotype(haplotype, fullReferenceWithPadding, 0));
        haplotypeList.add(haplotype);

        AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(haplotypeList);
        SampleList samples = new IndexedSampleList(Arrays.asList("sm1"));

        List<GATKRead> readList = new ArrayList<>(30);
        Map<String, List<GATKRead>>ebs = new HashMap<>();
        ebs.put("sm1", readList);

        for (int i = 0 ; i < 30; i++) {
            readList.add(ArtificialReadUtils.createArtificialRead("20M"));
        }

        double[][] values = {{3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0},
                {3, 0, 4, 0, 3, 0, 4, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0}
        };

        AlleleLikelihoods<GATKRead, Haplotype> lks = new AlleleLikelihoods<>(samples, haplotypes, ebs);
        LikelihoodMatrix<GATKRead, Haplotype> lkm = lks.sampleMatrix(0);
        for (int i = 0; i < lks.numberOfAlleles(); i++){
            for (int j = 0 ; j < lkm.evidenceCount(); j++) {
                lkm.set(i,j,values[i][j]);
            }
        }

        HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
        HaplotypeCallerGenotypingEngine genotypingEngine = new HaplotypeCallerGenotypingEngine(hcArgs, samples, ! hcArgs.doNotRunPhysicalPhasing, false);


        AlleleFiltering alleleFiltering = new AlleleFilteringHC(hcArgs, null, genotypingEngine);
        AlleleLikelihoods<GATKRead, Haplotype> filtered_lks = alleleFiltering.filterAlleles(lks, 100, new HashSet<>());
        //hcArgs.filterLoneAlleles is false, so we keep the weak lone allele
        Assert.assertEquals(filtered_lks.alleles(), haplotypeList);

        hcArgs.filterLoneAlleles = true;
        alleleFiltering = new AlleleFilteringHC(hcArgs, null, genotypingEngine);
        filtered_lks = alleleFiltering.filterAlleles(lks, 100, new HashSet<>());

        //hcArgs.filterLoneAlleles is true, so we keep the remove the weak lone allele
        Assert.assertEquals(filtered_lks.alleles(), haplotypeList.subList(0,1));

    }

    @Test //check that we filter strong allele with high SOR 
    public void testFilterDistantHindelSor() {

        // create haplotypes
        List<Haplotype> haplotypeList = new ArrayList<>();
        final byte[] fullReferenceWithPadding = "CAGGCATG".getBytes();
        Haplotype haplotype = new Haplotype(fullReferenceWithPadding, true, 0, TextCigarCodec.decode("8M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 108));
        haplotype.setEventMap(EventMap.fromHaplotype(haplotype, fullReferenceWithPadding, 0));
        haplotypeList.add(haplotype);

        haplotype = new Haplotype("CAGGCATTG".getBytes(), false, 0, TextCigarCodec.decode("7M1I1M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 109));

        haplotype.setEventMap(EventMap.fromHaplotype(haplotype, fullReferenceWithPadding, 0));
        haplotypeList.add(haplotype);
        haplotype = new Haplotype("CAGGCATTTG".getBytes(), false, 0, TextCigarCodec.decode("7M2I1M"));
        haplotype.setGenomeLocation(new SimpleInterval("chr", 100, 110));

        haplotype.setEventMap(EventMap.fromHaplotype(haplotype, fullReferenceWithPadding, 0));
        haplotypeList.add(haplotype);

        AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(haplotypeList);
        SampleList samples = new IndexedSampleList(Arrays.asList("sm1"));

        List<GATKRead> readList = new ArrayList<>(30);
        Map<String, List<GATKRead>> ebs = new HashMap<>();
        ebs.put("sm1", readList);

        for (int i = 0; i < 40; i++) {
            readList.add(ArtificialReadUtils.createArtificialRead("20M"));
        }


        double[][] values = {
                { 0,  3,  0, 3,  0, 3,  0, 3,  0, 3,  0, 3,  0, 3,  0, 3,  0, 3,  0, 3,  0, 3,  0, 3,  0, 3,  0, 3,  0, 3, 0, 3, 0, 3,
                        0, 3, 0, 3, 0,
                        3 },
                { 3,  0,  3, 0,  3, 0,  3, 0,  3, 0,  3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0,
                        3, 0, 3, 0, 3,
                        0 },
                {  10, 0,  0, 0,  10, 0,  0, 0,  10, 0,  0,  0,  10, 0,  0, 0,  10, 0,  0, 0,  10, 0,  0, 0,  10, 0,  0, 0, 10, 0, 0,0,10,0,0,0,10,0,0,0}
        };
        for (int i = 0; i < 40; i++) {
            if (i % 4 == 0) {
                readList.get(i).setIsReverseStrand(true);
            } 
        }

        AlleleLikelihoods<GATKRead, Haplotype> lks = new AlleleLikelihoods<>(samples, haplotypes, ebs);
        LikelihoodMatrix<GATKRead, Haplotype> lkm = lks.sampleMatrix(0);
        for (int i = 0; i < lks.numberOfAlleles(); i++) {
            for (int j = 0; j < lkm.evidenceCount(); j++) {
                lkm.set(i, j, values[i][j]);
            }
        }

        HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
        HaplotypeCallerGenotypingEngine genotypingEngine = new HaplotypeCallerGenotypingEngine(hcArgs, samples,
                !hcArgs.doNotRunPhysicalPhasing, false);

        AlleleFiltering alleleFiltering = new AlleleFilteringHC(hcArgs, null, genotypingEngine);
        AlleleLikelihoods<GATKRead, Haplotype> filtered_lks = alleleFiltering.filterAlleles(lks, 0, new HashSet<>());
        Assert.assertEquals(filtered_lks.alleles(), haplotypeList.subList(0, 2));
    }

    @Test
    public void testIdentifyBadAlleles(){
        Event a = new Event("chr1", 10, Allele.create("A",true), Allele.create("T", false));
        Event b = new Event("chr1", 10, Allele.create("T",true), Allele.create("G", false));
        Event c = new Event("chr1", 10, Allele.create("C", true), Allele.create("G", false));

        List<Event> events = List.of(a,b,c);
        List<Integer> rpls = List.of(10,20,0);
        List<Double> sors = List.of(0.0,1.0,3.5);
        HaplotypeCallerGenotypingEngine ge = new HaplotypeCallerGenotypingEngine(new HaplotypeCallerArgumentCollection(),
                SampleList.singletonSampleList("test"), false, false);
        AlleleFiltering af = new AlleleFilteringHC(null, null,ge);
        List<Event> badAlleles = af.identifyBadAlleles(rpls, sors, events, 30, 3);
        Assert.assertEquals(badAlleles, List.of(b, a, c));
        rpls = List.of(-100, -200, 0);
        sors = List.of(0.0,1.0,3.5);
        badAlleles = af.identifyBadAlleles(rpls, sors, events, 30, 3);
        Assert.assertEquals(badAlleles, List.of(c));

        rpls = List.of(-100, -200, -300);
        sors = List.of(0.0,1.0,3.5);
        badAlleles = af.identifyBadAlleles(rpls, sors, events, 30, 3);
        Assert.assertEquals(badAlleles, List.of(c));


    }
}
