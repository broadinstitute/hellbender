package org.broadinstitute.hellbender.tools.funcotator;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.CountingVariantFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;
import org.broadinstitute.hellbender.tools.funcotator.metadata.VcfFuncotationMetadata;
import org.broadinstitute.hellbender.transformers.VariantTransformer;
import org.broadinstitute.hellbender.utils.SequenceDictionaryUtils;
import org.broadinstitute.hellbender.utils.Utils;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.nio.file.Path;
import java.util.*;

/**
 * Funcotator (FUNCtional annOTATOR) analyzes given variants for their function (as retrieved from a set of data sources) and produces the analysis in a specified output file.
 *
 * <p>
 *     This tool is a functional annotation tool that allows a user to add annotations to called variants based on a set of data sources, each with its own matching criteria.
 * </p>
 *
 * <h3>Detailed Information and Tutorial</h3>
 * <p>Detailed information and a tutorial can be found here:
 *     <ul>
 *         <li><a href="https://gatk.broadinstitute.org/hc/en-us/articles/360035889931-Funcotator-Information-and-Tutorial">https://gatk.broadinstitute.org/hc/en-us/articles/360035889931-Funcotator-Information-and-Tutorial</a></li>
 *     </ul>
 * </p>
 *
 * <h3>Data Sources</h3>
 * <p>
 *     Data sources are expected to be in folders that are specified as input arguments.  While multiple data source folders can be specified, <b>no two data sources can have the same name</b>.
 * </p>
 * <h4>Data Source Folders</h4>
 * <p>
 *     In each main data source folder, there should be sub-directories for each individual data source, with further sub-directories for a specific reference (i.e. <i>hg19</i> or <i>hg38</i>).
 *     In the reference-specific data source directory, there is a configuration file detailing information about the data source and how to match it to a variant.  This configuration file is required.
 * </p>
 * <p>
 *     An example of a data source directory is the following:
 *
 *     <pre>
 *         dataSourcesFolder/
 *              Data_Source_1/
 *                  hg19
 *                      data_source_1.config
 *                      data_source_1.data.file.one
 *                      data_source_1.data.file.two
 *                      data_source_1.data.file.three
 *                      ...
 *                   hg38
 *                      data_source_1.config
 *                      data_source_1.data.file.one
 *                      data_source_1.data.file.two
 *                      data_source_1.data.file.three
 *                      ...
 *              Data_Source_2/
 *                  hg19
 *                      data_source_2.config
 *                      data_source_2.data.file.one
 *                      data_source_2.data.file.two
 *                      data_source_2.data.file.three
 *                      ...
 *                   hg38
 *                      data_source_2.config
 *                      data_source_2.data.file.one
 *                      data_source_2.data.file.two
 *                      data_source_2.data.file.three
 *                      ...
 *               ...
 *     </pre>
 * </p>
 * <h4>Pre-packaged Data Sources</h4>
 * <p>
 *     The GATK includes two sets of pre-packaged data sources, allowing for {@link Funcotator} use without (much) additional configuration.
 *     These data source packages correspond to the <strong>germline</strong> and <strong>somatic</strong> use cases.
 *     Broadly speaking, if you have a <strong>germline VCF</strong>, the <strong>germline data sources</strong> are what you want to use to start with.
 *     Conversely, if you have a <strong>somatic VCF</strong>, the <strong>somatic data sources</strong> are what you want to use to start with.
 *
 *     <br /><b>Versioned gzip archives of data source files are provided here:</b>
 *     <ul>
 *         <li><a href="ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/funcotator/">ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/funcotator/</a></li>
 *         <li><a href="https://console.cloud.google.com/storage/browser/broad-public-datasets/funcotator">gs://broad-public-datasets/funcotator/</a></li>
 *     </ul>
 * </p>
 *
 * <h5>gnomAD</h5>
 * <p>
 *     The pre-packaged data sources include gnomAD, a large database of known variants.  gnomAD is split into two parts - one based on exome data, one based on whole genome data.
 *     Due to the size of gnomAD, it cannot be included in the data sources package directly.  Instead, the configuration data are present and point to a Google bucket in which
 *     the gnomAD data reside.  This will cause <i>{@link Funcotator}</i> to actively connect to that bucket when it is run.  <br/>
 *     For this reason, <strong>gnomAD is disabled by default</strong>.<br />
 *     To enable gnomAD, simply change directories to your data sources directory and untar the gnomAD tar.gz files:
 *     <pre>
 *         cd DATA_SOURCES_DIR
 *         tar -zxf gnomAD_exome.tar.gz
 *         tar -zxf gnomAD_genome.tar.gz
 *     </pre>
 * </p>
 * <p>
 *     Because <i>{@link Funcotator}</i> will query the Internet when gnomAD is enabled, performance will be impacted by the machine's Internet connection speed.
 *     If this degradation is significant, you can localize gnomAD to the machine running <i>{@link Funcotator}</i> to improve performance (however due to the size of gnomAD this may be impractical).
 * </p>
 *
 * <h4>Data Source Downloader Tool</h4>
 * <p>
 *     To improve ease-of-use of <b><i>{@link Funcotator}</i></b>, there is a tool to download the pre-packaged data sources to the user's machine.
 *     This tool is the <strong><i>{@link FuncotatorDataSourceDownloader}</i></strong> and can be run to retrieve the pre-packaged data sources from the google bucket and localize them to the machine on which it is run.
 *     <br />Briefly:
 *     <ul>
 *         <li>For <strong>somatic</strong> data sources:<br /><pre>{@code ./gatk FuncotatorDataSourceDownloader --somatic --validate-integrity --extract-after-download}</pre></li>
 *         <li>For <strong>germline</strong> data sources:<br /><pre>{@code ./gatk FuncotatorDataSourceDownloader --germline --validate-integrity --extract-after-download}</pre></li>
 *     </ul>
 * </p>
 * <h4>Disabling Data Sources</h4>
 * <p>
 *     A data source can be disabled by removing the folder containing the configuration file for that source.  This can be done on a per-reference basis.  If the entire data source should be disabled, the entire top-level data source folder can be removed.
 * </p>
 * <p>
 *     If it is possible that the data source will be re-enabled in the future, then we recommend zipping the data source folder and removing the folder itself, leaving only the zip file in its place.  When the time comes to enable data source again, simply unzip the file and the data source will be ready to go the next time <i>{@link Funcotator}</i> is run.
 * </p>
 * <h4>User-Defined Data Sources</h4>
 * <p>
 *     Users can define their own data sources by creating a new correctly-formatted data source sub-directory in the main data sources folder.  In this sub-directory, the user must create an additional folder for the reference for which the data source is valid.  If the data source is valid for multiple references, then multiple reference folders should be created.
 *     Inside each reference folder, the user should place the file(s) containing the data for the data source.  Additionally the user <b>must</b> create a configuration file containing metadata about the data source.
 * </p>
 * <p>
 *     <i>{@link Funcotator}</i> allows for data sources with source files that live on the cloud, enabling users to annotate with data sources that are not physically present on the machines running <i>{@link Funcotator}</i>.<br />
 *     To create a data source based on the cloud, create a configuration file for that data source and put the cloud URL in as the src_file property (see <i>Configuration File Format</i> for details).<br />
 *     E.g.:
 *     <pre>
 *         ...
 *         src_file = gs://broad-references/hg19/v0/1000G_phase1.snps.high_confidence.b37.vcf.gz
 *         ...
 *     </pre>
 * </p>
 * <p>
 *     There are several formats allowed for data sources, however the two most useful are arbitrarily separated value (XSV) files, such as comma-separated value (CSV), tab-separated value (TSV).  These files contain a table of data that can be matched to a variant by <i>gene name</i>, <i>transcript ID</i>, or <i>genome position</i>.
 *     In the case of <i>gene name</i> and <i>transcript ID</i>, one column must contain the <i>gene name</i> or <i>transcript ID</i> for each row's data.
 *     <ul>
 *         <li>For <i>gene name</i>, when a variant is annotated with a gene name that <i>exactly matches</i> an entry in the gene name column for a row, that row's other fields will be added as annotations to the variant.</li>
 *         <li>For <i>transcript ID</i>, when a variant is annotated with a transcript ID that <i>exactly matches</i> an entry in the transcript ID column for a row, that row's other fields will be added as annotations to the variant.</li>
 *         <li>For <i>genome position</i>, one column must contain the contig ID, another column must contain the start position (1-based, inclusive), and a column must contain the stop position (1-based, inclusive).  The start and stop columns may be the same column.  When a variant is annotated with a genome position that <i>overlaps</i> an entry in the three genome position columns for a row, that row's other fields will be added as annotations to the variant.</li>
 *     </ul>
 * </p>
 * <h4>Configuration File Format</h4>
 * <p>
 *     The configuration file is a standard Java properties-style configuration file with key-value pairs.  This file name <b>must end in .config</b>.
 * </p>
 * <p>
 *     The following is an example of a genome position XSV configuration file (for the ORegAnno data source):
 *     <pre>
 *         name = Oreganno
 *         version = 20160119
 *         src_file = oreganno.tsv
 *         origin_location = http://www.oreganno.org/dump/ORegAnno_Combined_2016.01.19.tsv
 *         preprocessing_script = getOreganno.py
 *
 *         # Supported types:
 *         # simpleXSV    -- Arbitrary separated value table (e.g. CSV), keyed off Gene Name OR Transcript ID
 *         # locatableXSV -- Arbitrary separated value table (e.g. CSV), keyed off a genome location
 *         # gencode      -- Custom datasource class for GENCODE
 *         # cosmic       -- Custom datasource class for COSMIC
 *         # vcf          -- Custom datasource class for Variant Call Format (VCF) files
 *         type = locatableXSV
 *
 *         # Required field for GENCODE files.
 *         # Path to the FASTA file from which to load the sequences for GENCODE transcripts:
 *         gencode_fasta_path =
 *
 *         # Required field for GENCODE files.
 *         # NCBI build version (either hg19 or hg38):
 *         ncbi_build_version =
 *
 *         # Required field for simpleXSV files.
 *         # Valid values:
 *         #     GENE_NAME
 *         #     TRANSCRIPT_ID
 *         xsv_key =
 *
 *         # Required field for simpleXSV files.
 *         # The 0-based index of the column containing the key on which to match
 *         xsv_key_column =
 *
 *         # Required field for simpleXSV AND locatableXSV files.
 *         # The delimiter by which to split the XSV file into columns.
 *         xsv_delimiter = \t
 *
 *         # Required field for simpleXSV files.
 *         # Whether to permissively match the number of columns in the header and data rows
 *         # Valid values:
 *         #     true
 *         #     false
 *         xsv_permissive_cols = true
 *
 *         # Required field for locatableXSV files.
 *         # The 0-based index of the column containing the contig for each row
 *         contig_column = 1
 *
 *         # Required field for locatableXSV files.
 *         # The 0-based index of the column containing the start position for each row
 *         start_column = 2
 *
 *         # Required field for locatableXSV files.
 *         # The 0-based index of the column containing the end position for each row
 *         end_column = 3
 *     </pre>
 * </p>
 *
 * <h3>Required Inputs</h3>
 * <ul>
 *     <li>A reference genome sequence.</li>
 *     <li>The version of the reference genome sequence being used (e.g. <i>hg19</i>, <i>hg38</i>, etc.).</li>
 *     <li>A VCF of variant calls to annotate.</li>
 *     <li>The path to a folder of data sources formatted for use by Funcotator.</li>
 *     <li>The desired output format for the annotated variants file (either <i>MAF</i> or <i>VCF</i>)</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <p>
 *     The basic output of <i>{@link Funcotator}</i> is:
 * <ul>
 *     <li>A VCF or MAF file containing all variants from the input file with added annotations corresponding to annotations from each data source that matched a given variant according to that data source's matching criteria.</li>
 * </ul>
 * </p>
 *
 * <h4>Annotations for Pre-Packaged Data Sources</h4>
 * <p>
 *     The pre-packaged data sources will create a set of baseline, or default annotations for an input data set.
 *     Most of these data sources copy and paste values from their source files into the output of <i>{@link Funcotator}</i> to create annotations.  In this sense they are trivial data sources.
 * </p>
 * <h5>Gencode</h5>
 * <p>
 *     <i>{@link Funcotator}</i> performs some processing on the input data to create the Gencode annotations.  Gencode is currently required, so <i>{@link Funcotator}</i> will create these annotations for all input variants.
 *     The order and a specification of the Gencode annotations that <i>{@link Funcotator}</i> creates is as follows:
 *     <ol>
 *         <li>
 *             <b><i>hugoSymbol</i></b><br />
 *             <i>Type: </i>{@link String}<br />
 *             The name of the gene in which the annotated variant allele occurs.  If the variant allele occurs outside of any known gene boundaries, then this field is set to "Unknown".
 *         </li>
 *         <li>
 *             <b><i>ncbiBuild</i></b><br />
 *             <i>Type: </i>{@link String}<br />
 *             The reference which was used to create this Gencode annotation.  Current valid values are: "<i>hg19</i>" or "<i>hg38</i>".
 *         </li>
 *         <li>
 *             <b><i>chromosome</i></b><br />
 *             <i>Type: </i>{@link String}<br />
 *             The contig in which the variant occurs.  Will always correspond to the contig in the variant position.
 *         </li>
 *         <li>
 *             <b><i>start</i></b><br />
 *             <i>Type: </i>{@link Integer}<br />
 *             The start position in genomic coordinates of the variant allele being annotated (1-based, inclusive).  Will always correspond to the start in the variant position.
 *         </li>
 *         <li>
 *             <b><i>end</i></b><br />
 *             <i>Type: </i>{@link Integer}<br />
 *             The end position in genomic coordinates of the variant allele being annotated (1-based, inclusive).  Will always correspond to the position last base in the variant allele.
 *         </li>
 *         <li>
 *             <b><i>variantClassification</i></b><br />
 *             <i>Type: </i>{@link String}<br />
 *             The classification of the variant being annotated.  Will always be one of the following:
 *             <ul>
 *
 *                  <li><p><i>COULD_NOT_DETERMINE</i><br />
 *                  Variant classification could not be determined.</p></li>
 *
 *                  <li><p><i>INTRON</i><br />
 *                  Variant lies between exons within the bounds of the chosen transcript.
 *                  Only valid for Introns.</p></li>
 *
 *                  <li><p><i>FIVE_PRIME_UTR</i><br />
 *                  Variant is on the 5'UTR for the chosen transcript.
 *                  Only valid for UTRs.</p></li>
 *
 *                  <li><p><i>THREE_PRIME_UTR</i><br />
 *                  Variant is on the 3'UTR for the chosen transcript
 *                  Only valid for UTRs.</p></li>
 *
 *                  <li><p><i>IGR</i><br />
 *                  Intergenic region. Does not overlap any transcript.
 *                  Only valid for IGRs.</p></li>
 *
 *                  <li><p><i>FIVE_PRIME_FLANK</i><br />
 *                  The variant is upstream of the chosen transcript
 *                  Only valid for IGRs.</p></li>
 *
 *                  <li><p><i>THREE_PRIME_FLANK</i><br />
 *                  The variant is downstream of the chosen transcript
 *                  Only valid for IGRs.</p></li>
 *
 *                  <li><p><i>MISSENSE</i><br />
 *                  The point mutation alters the protein structure by one amino acid.
 *                  Can occur in Coding regions or Introns.</p></li>
 *
 *                  <li><p><i>NONSENSE</i><br />
 *                  A premature stop codon is created by the variant.
 *                  Can occur in Coding regions or Introns.</p></li>
 *
 *                  <li><p><i>NONSTOP</i><br />
 *                  Variant removes stop codon.
 *                  Can occur in Coding regions or Introns.</p></li>
 *
 *                  <li><p><i>SILENT</i><br />
 *                  Variant is in coding region of the chosen transcript, but protein structure is identical.
 *                  Can occur in Coding regions or Introns.</p></li>
 *
 *                  <li><p><i>SPLICE_SITE</i><br />
 *                  The variant is within a configurable number of bases  of a splice site. See the secondary classification to determine if it lies on the exon or intron side.
 *                  Can occur in Coding regions or Introns.</p></li>
 *
 *                  <li><p><i>IN_FRAME_DEL</i><br />
 *                  Deletion that keeps the sequence in frame.
 *                  Can occur in Coding regions or Introns.</p></li>
 *
 *                  <li><p><i>IN_FRAME_INS</i><br />
 *                  Insertion that keeps the sequence in frame.
 *                  Can occur in Coding regions or Introns.</p></li>
 *
 *                  <li><p><i>FRAME_SHIFT_INS</i><br />
 *                  Insertion that moves the coding sequence out of frame.
 *                  Can occur in Coding regions or Introns.</p></li>
 *
 *                  <li><p><i>FRAME_SHIFT_DEL</i><br />
 *                  Deletion that moves the sequence out of frame.
 *                  Can occur in Coding regions or Introns.</p></li>
 *
 *                  <li><p><i>START_CODON_SNP</i><br />
 *                  Point mutation that overlaps the start codon.
 *                  Can occur in Coding regions or Introns.</p></li>
 *
 *                  <li><p><i>START_CODON_INS</i><br />
 *                  Insertion that overlaps the start codon.
 *                  Can occur in Coding regions or Introns.</p></li>
 *
 *                  <li><p><i>START_CODON_DEL</i><br />
 *                  Deletion that overlaps the start codon.
 *                  Can occur in Coding regions or Introns.</p></li>
 *
 *                  <li><p><i>DE_NOVO_START_IN_FRAME</i><br />
 *                  New start codon is created by the given variant using the chosen transcript.
 *                  However, it is in frame relative to the coded protein, meaning that if the coding sequence were extended
 *                  then the new start codon would be in frame with the
 *                  existing start and stop codons.
 *
 *                  This can only occur in a 5' UTR.</p></li>
 *
 *                  <li><p><i>DE_NOVO_START_OUT_FRAME</i><br />
 *                  New start codon is created by the given variant using the chosen transcript.
 *                  However, it is out of frame relative to the coded protein, meaning that if the coding sequence were extended
 *                  then the new start codon would NOT be in frame with
 *                  the existing start and stop codons.
 *
 *                  This can only occur in a 5' UTR.</p></li>
 *
 *                  <li><p><i>RNA</i><br />
 *                  Variant lies on one of the RNA transcripts.
 *                  (special catch-all case)</p></li>
 *
 *                  <li><p><i>LINCRNA</i><br />
 *                  Variant lies on one of the lincRNAs.
 *                  (special catch-all case)</p></li>
 *             </ul>
 *         </li>
 *         <li>
 *             <b><i>secondaryVariantClassification</i></b><br />
 *             <i>Type: </i>{@link String}<br />
 *             Additional variant classification information for variant alleles that have a {@code variantClassification} of {@code SPLICE_SITE}.
 *             For a variant allele with the {@code variantClassification} of {@code SPLICE_SITE}, this will indicate the specific classification of the variant.
 *             For all variants that do not have the {@code variantClassification} of {@code SPLICE_SITE}, this will be the empty string.
 *         </li>
 *         <li>
 *             <b><i>variantType</i></b><br />
 *             <i>Type: </i>{@link String}<br />
 *             Basic information about the variant allele being annotated.  Can be one of:
 *             <ul>
 *                 <li><i>INS</i> - The variant allele is some kind of insertion.</li>
 *                 <li><i>DEL</i> - The variant allele is some kind of deletion.</li>
 *                 <li><i>SNP</i> - The variant allele is a single nucleotide polymorphism.</li>
 *                 <li><i>DNP</i> - The variant allele is a di-nucleotide polymorphism.</li>
 *                 <li><i>TNP</i> - The variant allele is a tri-nucleotide polymorphism.</li>
 *                 <li><i>ONP</i> - The variant allele is an oligo-nucleotide polymorphism (Synonymous with MNP).</li>
 *                 <li><i>MNP</i> - The variant allele is a multi-nucleotide polymorphism (Synonymous with ONP).</li>
 *                 <li><i>NA</i> - The variant allele type cannot be determined.</li>
 *             </ul>
 *         </li>
 *         <li>
 *             <b><i>refAllele</i></b><br />
 *             <i>Type: </i>{@link String}<br />
 *             The reference allele for the position at which this this variant allele occurs.<br />
 *             For insertions, this will be set to <i>"-"</i>.
 *         </li>
 *         <li>
 *             <b><i>tumorSeqAllele1</i></b><br />
 *             <i>Type: </i>{@link String}<br />
 *             Always the same as the reference allele.  This field is a hold-over required for MAF annotations.<br />
 *             For insertions, this will be set to <i>"-"</i>.
 *         </li>
 *         <li>
 *             <b><i>tumorSeqAllele2</i></b><br />
 *             <i>Type: </i>{@link String}<br />
 *             The variant allele being annotated.  This field only includes the bases that are different from the reference.<br />
 *             For the input VCF records, this field may slightly differ from the alternate allele reported in the base data for the {@link VariantContext}.<br />
 *             For deletions, this will be set to <i>"-"</i>.
 *         </li>
 *         <li>
 *             <b><i>genomeChange</i></b><br />
 *             <i>Type: </i>{@link String}<br />
 *             A {@link String} summarizing the change resulting from this variant allele within the context of the whole genome sequence.<br />
 *             Generally the format of this field is:<br />
 *             <pre>
 *             g.[CONTIG]:[POSITION][BASES CHANGED]
 *             </pre>
 *
 *             The format of this field slightly varies based on {@code VariantType}:
 *             <ul>
 *                 <li>
 *                     <i><u>Insertion</u></i><br />
 *                     <pre>g.[CONTIG]:[POSITION OF BASE PRIOR TO INSERTION];_[POSITION OF BASE AFTER INSERTION]ins[BASES INSERTED]</pre><br />
 *                     E.g.:<br />
 *                     <pre>g.chr19:2018023_2018024insAATCG</pre><br />
 *                     This indicates that the bases AATCG were inserted between bases 2018023 and 2018024 on chromosome 19.
 *                 </li>
 *                 <li>
 *                     <i><u>Deletion</u></i><br />
 *                     <pre>g.[CONTIG]:[POSITION OF BASE DELETED]del[BASE DELETED]</pre><br />
 *                     E.g.:<br />
 *                     <pre>g.chr19:2018023delT</pre><br />
 *                     This indicates that the base T was deleted at position 2018023 on chromosome 19.
 *                     <i>OR</i>
 *                     <pre>g.[CONTIG]:[POSITION OF FIRST BASE DELETED]_[POSITION OF LAST BASE DELETED]del[BASES DELETED]</pre><br />
 *                     E.g.:<br />
 *                     <pre>g.chr19:2018023_2018025delTTG</pre><br />
 *                     This indicates that the bases TTG were deleted starting at position 2018023 and ending at position 2018025 on chromosome 19.
 *                 </li>
 *                 <li>
 *                     <i><u>SNP</u></i><br />
 *                     <pre>g.[CONTIG]:[POSITION OF BASE ALTERED][REFERENCE BASE]>[ALTERNATE BASE]</pre><br />
 *                     E.g.:<br />
 *                     <pre>g.chr19:2018023T>G</pre><br />
 *                     This indicates that the base T was changed to G at position 2018023 on chromosome 19.
 *                 </li>
 *                 <li>
 *                     <i><u>MNP (including DNPs, TNPs)</u></i><br />
 *                     <pre>g.[CONTIG]:[POSITION OF FIRST BASE ALTERED]_[POSITION OF LAST BASE ALTERED][REFERENCE BASES>&gt;[ALTERNATE BASES]</pre><br />
 *                     E.g.:<br />
 *                     <pre>g.chr19:2018023_2018025TTG>GAT</pre><br />
 *                     This indicates that the bases TTG were changed to GAT from position 2018023 to position 2018025 on chromosome 19.
 *                 </li>
 *             </ul>
 *         </li>
 *         <li>
 *             <b><i>annotationTranscript</i></b><br />
 *             <i>Type: </i>{@link String}<br />
 *             The ID of the transcript chosen for the detailed Gencode annotation reporting.<br />E.g.: <i>ENST00000435064.1</i><br />
 *             If the variant allele does not occur within the bounds of any transcript (e.g. is of type <i>IGR</i>), then this field is empty.
 *         </li>
 *         <li>
 *             <b><i>transcriptStrand</i></b><br />
 *             <i>Type: </i>{@link String}<br />
 *             The strand direction associated with the transcript on which this variant allele occurs.<br />
 *             Either "<i>+</i>" or "<i>-</i>".
 *         </li>
 *         <li>
 *             <b><i>transcriptExon</i></b><br />
 *             <i>Type: </i>{@link Integer} or Empty<br />
 *             The exon number on the transcript in which this variant allele occurs (1-based).  Corresponds directly to the Gencode exon number.<br />
 *             If the variant does not occur in the expressed transcript of the corresponding gene (e.g. is of type <i>INTRON</i> or <i>IGR</i>), then this field is empty.
 *         </li>
 *         <li>
 *             <b><i>transcriptPos</i></b><br />
 *             <i>Type: </i>{@link Integer} or Empty<br />
 *             Position in the chosen transcript of the variant allele.  All positions listed are 1-based and inclusive (meaning that the first base in the transcript starts at and ends at position 1).<br />
 *             For variant alleles that occur at a single base, the format is simply the position at which that variant occurs in the transcript (e.g. <i>1294</i>)<br />
 *             For variant alleles spanning multiple bases, the format is:<br />
 *             <pre>
 *             [START]_[END]
 *             </pre>
 *             E.g.: <pre>1236_1237</pre><br />
 *             If the variant does not occur in the expressed transcript of the corresponding gene (e.g. is of type <i>INTRON</i> or <i>IGR</i>), then this field is empty.
 *         </li>
 *         <li>
 *             <b><i>cDnaChange</i></b><br />
 *             <i>Type: </i>{@link String}<br />
 *             A {@link String} that summarizes the change resulting from this variant allele in the coding sequence for the transcript in which it occurs.<br />
 *             Positions in this field are <i>relative to the start of the transcript (1-based, inclusive)</i> unless otherwise noted.<br />
 *             Generally the format of this field is:<br />
 *             <pre>
 *             c.[POSITION][BASES CHANGED]
 *             </pre>
 *
 *             The format of this field slightly varies based on {@code VariantType}, the number of affected bases, and whether the variant allele is a <i>SPLICE_SITE</i>:
 *             <ul>
 *                 <li>
 *                     <i><u>Insertions</u></i><br />
 *                     <pre>c.[POSITION OF BASE PRIOR TO INSERTION]_[POSITION OF BASE AFTER INSERTION]ins[BASES INSERTED]</pre><br />
 *                     E.g.:<br />
 *                     <pre>c.2018_2019insAA</pre><br />
 *                     This indicates that the bases AA were inserted between bases 2018 and 2019 in the transcript associated with this variant allele.
 *                 </li>
 *                 <li>
 *                     <i><u>Deletions of One Base</u></i><br />
 *                     <pre>c.[POSITION OF BASE DELETED]del[BASE DELETED]</pre><br />
 *                     E.g.:<br />
 *                     <pre>c2018delT</pre><br />
 *                     This indicates that the base T was deleted at position 2018 in the transcript associated with this variant allele.
 *                 </li>
 *                 <li>
 *                     <i><u>Deletions of Multiple Bases</u></i><br />
 *                     <pre>c.[POSITION OF FIRST BASE DELETED]_[POSITION OF LAST BASE DELETED]del[BASES DELETED]</pre><br />
 *                     E.g.:<br />
 *                     <pre>c2018_2022delTTCAG</pre><br />
 *                     This indicates that the bases TTCAG were deleted from position 2018 to position 2022 in the transcript associated with this variant allele.
 *                 </li>
 *                 <li>
 *                     <i><u>SNPs</u></i><br />
 *                     <pre>c.[POSITION OF BASE CHANGED]&gt;[NEW BASE]</pre><br />
 *                     E.g.:<br />
 *                     <pre>c.1507T>G</pre><br />
 *                     This indicates that the base T was changed to G at position 1507 in the transcript associated with this variant allele.
 *                 </li>
 *                 <li>
 *                     <i><u>MNPs (including DNPs, TNPs)</u></i><br />
 *                     <pre>c.[POSITION OF FIRST BASE CHANGED]_[POSITION OF LAST BASE CHANGED]&gt;[NEW BASES]</pre><br />
 *                     E.g.:<br />
 *                     <pre>c.12899_12900AG>TA</pre><br />
 *                     This indicates that the bases AG were changed to TA from position 12899 to position 12900 in the transcript associated with this variant allele.
 *                 </li>
 *                 <li>
 *                     <i><u>INTRONIC Variants</u></i><br />
 *                     For variants occuring in INTRONs, the cDNA string position indicates the offset from the exon boundary for the start of the variant.  The whole string takes the form:<br />
 *                     <pre>c.e[EXON NUMBER][+|-][BASES FROM EXON][REF ALLELE]>[ALT ALLELE]</pre><br />
 *                     E.g.:<br />
 *                     <pre>c.e81-4TAA>A</pre><br />
 *                     This indicates that the bases TAA were changed to A starting four bases before exon 81 in the transcript associated with this variant allele.
 *                 </li>
 *             </ul>
 *             If the variant does not occur in the expressed transcript of the corresponding gene (e.g. is of type <i>IGR</i>), then this field is empty.
 *         </li>
 *         <li>
 *             <b><i>codonChange</i></b><br />
 *             <i>Type: </i>{@link String}<br />
 *             A {@link String} that representing the <i>codon-aligned change</i> resulting from this variant allele in the coding sequence for the transcript in which it occurs.<br />
 *             Positions in this field are <i>relative to the start of the transcript (1-based, inclusive) and aligned to the coding sequence</i>unless otherwise noted.<br />
 *             Unlike the cDnaChange, the bases reported in the codonChange string will always have a length evenly divisble by 3 (except for frameshifts) and represent what the codons would be if the variant alternate allele were expressed in the reference sequence.<br />
 *             Capitalized bases represent the bases changed by the variant alternate allele.  Lower-case bases represent reference bases.<br />
 *             Generally the format of this field is:<br />
 *             <pre>
 *             c.[POSITION][BASES CHANGED]
 *             </pre>
 *             The format of this field slightly varies based on {@code VariantType}, the number of affected bases, and whether the variant allele occurs in an Intron:
 *             <ul>
 *                 <li>
 *                     <i><u>In-Frame Insertions</u></i><br />
 *                     <pre>c.([POSITION OF FIRST BASE IN FIRST CODON IN THE REFERENCE AFFECTED BY THIS VARIANT]-[POSITION OF LAST BASE IN LAST CODON IN THE REFERENCE AFFECTED BY THIS VARIANT][REFERENCE CODONS]&gt;[EXPRESSED CODONS]</pre><br />
 *                     E.g.:<br />
 *                     <pre>c.(19-21)ctt>ctCGTt</pre><br />
 *                     This indicates that the bases <i>CGT</i> were inserted before the 6th codon (starting at base 19, ending at base 21) in the transcript associated with this variant allele, and the resulting expressed codons would be <i>ctCGTt</i>.
 *                 </li>
 *                 <li>
 *                     <i><u>In-Frame Deletions of Complete Codons (codon-aligned deletions)</u></i><br />
 *                     <pre>c.([POSITION OF FIRST BASE IN FIRST CODON DELETED]-[POSITION OF LAST BASE IN LAST CODON DELETED][REFERENCE CODONS]del</pre><br />
 *                     E.g.:<br />
 *                     <pre>c.(997-999)gcadel</pre><br />
 *                     This indicates that the 332nd codon (starting at base 997, ending at base 999) was deleted in the transcript associated with this variant allele, and the deleted codon bases are <i>gca</i>.
 *                 </li>
 *                 <li>
 *                     <i><u>In-Frame Deletions Spanning Multiple Codons</u></i><br />
 *                     <pre>c.([POSITION OF FIRST BASE IN FIRST CODON DELETED]-[POSITION OF LAST BASE IN LAST CODON DELETED][REFERENCE CODONS]&gt;[EXPRESSED CODONS]</pre><br />
 *                     E.g.:<br />
 *                     <pre>c.(997-1002)gcactc>gtc</pre><br />
 *                     This indicates that bases in the 332nd codon (starting at base 997) and 333rd codon (ending at base 1002) were deleted in the transcript associated with this variant allele, and the resulting expressed codon would be <i>gtc</i>.
 *                 </li>
 *                 <li>
 *                     <i><u>Frame Shift Insertions and Deletions</u></i><br />
 *                     <pre>c.([POSITION OF FIRST BASE IN LAST CORRECTLY EXPRESSED/REFERENCE CODON]-[POSITION OF LAST BASE IN LAST CORRECTLY EXPRESSED/REFERENCE CODON][REFERENCE CODONS]&gt;[EXPRESSED CODONS]</pre><br />
 *                     E.g.:<br />
 *                     <pre>c.(997-999)gcafs</pre><br />
 *                     This indicates that bases just <u><i>AFTER</i></u> the 332nd codon (starting at base 997, ending at base 999) were inserted or deleted in the transcript associated with this variant allele resulting in a frame shift, and that the last correctly transcribed codon would be codon 332 (starting at base 997, ending at base 999), '<i>gca</i>'.
 *                 </li>
 *                 <i><u>SNPs / MNPs</u></i><br />
 *                     <pre>c.([POSITION OF FIRST BASE IN FIRST CODON IN THE REFERENCE AFFECTED BY THIS VARIANT]-[POSITION OF LAST BASE IN LAST CODON IN THE REFERENCE AFFECTED BY THIS VARIANT][REFERENCE CODONS]&gt;[EXPRESSED CODONS]</pre><br />
 *                     E.g. 1:<br />
 *                     <pre>c.(39871-39873)cCC>cTT</pre><br />
 *                     This indicates that the bases <i>CC</i> were changed to <i>TT</i> in the 13290th codon (starting at base 39871, ending at base 39873) in the transcript associated with this variant allele, and the resulting expressed codon would be <i>cTT</i>.<br />
 *                     E.g. 2:<br />
 *                     <pre>c.(4-9)ctAAgc>ctGCgc</pre><br />
 *                     This indicates that the bases <i>AA</i> starting in the 2nd codon (starting at base 4) and ending in the 3rd codon (ending at base 9) were changed to <i>GC</i> in the transcript associated with this variant allele, and the resulting expressed codons would be <i>ctGCgc</i>.
 *                 </li>
 *             </ul>
 *         </li>
 *         <li>
 *             <b><i>proteinChange</i></b><br />
 *             <i>Type: </i>{@link String}<br />
 *             A short string representing the predicted amino acid sequence change in the product of the gene transcript in which this variant alternate allele occurs.<br />
 *             Positions in this field are <i>relative to the start of the amino acid sequence (1-based, inclusive) resulting from decoding the codons in the transcript in which this variant alternate allele occurs</i>unless otherwise noted.<br />
 *             Amino acid abbreviations are the standard letters as can be found in <a href="https://en.wikipedia.org/wiki/DNA_codon_table"></a>this table</a> with the exception of the stop codon, which is represented by '<i><b>*</b></i>'.<br />
 *             It is important to note that the positions and amino acids reported in this string may not directly align to the codons in which the variant alternate allele occurs.  This is most often due to the variant occuring in a set of tandem repeats which would cause the amino acid change to be "pushed" to the end of the tandem repeats.<br />
 *             For protein change strings in the Mitochondrial contig, <a href="https://en.wikipedia.org/wiki/Vertebrate_mitochondrial_code">the mitochondrial genetic code</a> is used (rather than the standard code).
 *             The format of this field takes two forms:<br />
 *             <ul>
 *                 <li>
 *                     <i>Protein Changes with One Amino Acid Changed</i><br />
 *                     <pre>p.[REFERENCE AMINO ACID][POSITION][PREDICTED EXPRESSED AMINO ACID]</pre><br />
 *                     E.g. 1:<br />
 *                     <pre>p.V5T</pre><br />
 *                     The amino acid at protein position 5 was V (Valine) in the reference and would become T (Threonine) with the variant alternate allele expressed.<br />
 *                     E.g. 2:<br />
 *                     <pre>p.R2R</pre><br />
 *                     The amino acid at protein position 2 was R (Argenine) in the reference and would become R (Argenine) with the variant alternate allele expressed (no change in amino acid sequence / silent variant classification).
 *                 </li>
 *                 <li>
 *                     <i>Protein Changes with Multiple Amino Acids Changed</i><br />
 *                     <pre>p.[FIRST AFFECTED AMINO ACID POSITION]_[LAST AFFECTED AMINO ACID POSITION][REFERENCE AMINO ACIDS]&gt;[PREDICTED EXPRESSED AMINO ACIDS]</pre><br />
 *                     E.g.:<br />
 *                     <pre>p.100_101Q*>FL</pre><br />
 *                     The amino acid sequence starting at protein position 100 and ending at protein position 101 was Q* (Glutamine,STOP) in the reference and would become FL (Phenylalanine,Leucine) with the variant alternate allele expressed.
 *                 </li>
 *             </ul>
 *             If the variant alternate allele does not occur in a coding region, this field will be empty.
 *         </li>
 *         <li>
 *             <b><i>gcContent</i></b><br />
 *             <i>Type: </i>{@link Double}<br />
 *             Represents the fraction of Guanine and Cytosine bases in a window of a given size around a variant.  This window size does not include any bases in the variant alternate allele itself.  By default the window size is 200 bases.
 *         </li>
 *         <li>
 *             <b><i>referenceContext</i></b><br />
 *             <i>Type: </i>{@link String}<br />
 *             The <u>strand-correct</u> reference coding sequence in a given window around the reference allele.  By default the window size is <i>10 bases</i>.<br />
 *             E.g. For the reference context around a variant with the reference allele '<i>C</i>' on the '+' strand:
 *
 *             <pre>
 *                      [REF ALLELE]
 *                           |
 *                           v
 *                 GAACCCACGTCGGTGAGGGCC
 *                 |________| |________|
 *                     v           v
 *                  10 bases      10 bases
 *                (window size)  (window size)
 *             </pre>
 *             Strand-correct specifically means that if the strand of this transcript is determined to be '<i>-</i>' then the sequence is reverse complemented.  <br />
 *             E.g. For the reference context around a variant with the reference allele '<i>C</i>' on the '-' strand:
 *             <pre>
 *                      [REF ALLELE]
 *                           |
 *                           v
 *                 CACGAAAGTCTTGCGGATCT
 *                 |________| |________|
 *                     v           v
 *                  10 bases      10 bases
 *                (window size)  (window size)
 *             </pre>
 *         </li>
 *         <li>
 *             <b><i>otherTranscripts</i></b><br />
 *             <i>Type: </i>{@link String}<br />
 *             A summary of the other transcripts in which this variant occurs, which were not chosen for detailed reporting due to the transcript selection scheme.<br />
 *             Each other transcript is represented by a condensed string that indicates how that transcript would be affected by this variant alternate allele.<br />
 *             Each other transcript takes the form:<br />
 *             <pre>[HUGO SYMBOL]_[TRANSCRIPT ID]_[VARIANT CLASSIFICATION]_[PROTEIN CHANGE STRING]</pre><br />
 *             E.g.:
 *             <pre>SDF4_ENST00000263741.7_MISSENSE_p.R243Q</pre>
 *             If another transcript were to be an IGR, the other transcript field would be populated with '<i>IGR_ANNOTATON</i>'<br />
 *             If the other transcript does not have a protein change string, then that part is not rendered.<br />
 *             In the event that there are multiple other transcripts, these transcripts are separated by '<i>/</i>'.<br />
 *             E.g.:
 *             <pre>SDF4_ENST00000263741.7_MISSENSE_p.R243Q/TNFRSF4_ENST00000379236.3_FIVE_PRIME_FLANK</pre>
 *             If this variant alternate allele occurs in only one transcript, this field will be empty.
 *         </li>
 *     </ol>
 * </p>
 * <p>
 *     Other annotations will follow the Gencode annotations and will be based on the data sources included in the data sources directory.
 * </p>
 * <h3>Usage example</h3>
 * <pre>
 *   ./gatk Funcotator \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   -O outputFile \
 *   --output-file-format MAF \
 *   --data-sources-path dataSourcesFolder/ \
 *   --ref-version hg19
 * </pre>
 *
 * <h3>Notes</h3>
 * <ul>
 *     <li>This tool is the spiritual successor to <a href="https://github.com/broadinstitute/oncotator">Oncotator</a>, with better support for germline data, numerous fixes for correctness, and many other features.</li>
 *     <li>REMEMBER: <strong>Funcotator is NOT Oncotator.</strong></li>
 * </ul>
 *
 * <h3>Known Issues</h3>
 * <p>A complete list of known open issues can be found on <a href="https://github.com/broadinstitute/gatk/issues?q=is%3Aopen+is%3Aissue+label%3AFuncotator">the GATK github entry for funcotator here.</a></p>
 */
@CommandLineProgramProperties(
        summary = "Create functional annotations on given variants cross-referenced by a given set of data sources.\n" +
                "A GATK functional annotation tool (similar functionality to Oncotator).",
        oneLineSummary = "Functional Annotator",
        programGroup = VariantEvaluationProgramGroup.class
)
@DocumentedFeature
public class Funcotator extends VariantWalker {
    private static final Logger logger = LogManager.getLogger(Funcotator.class);

    //==================================================================================================================
    // Arguments:

    @ArgumentCollection
    private final FuncotatorVariantArgumentCollection funcotatorArgs = new FuncotatorVariantArgumentCollection();

    //==================================================================================================================

    protected OutputRenderer outputRenderer;

    protected FuncotatorEngine funcotatorEngine;

    //==================================================================================================================

    /**
     * @return The {@link Funcotator}-specific arguments used to instantiate this {@link Funcotator} instance.
     */
    public FuncotatorVariantArgumentCollection getArguments() {
        return funcotatorArgs;
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {

        if (seqValidationArguments.performSequenceDictionaryValidation()) {
            logger.info("Validating sequence dictionaries...");
            // Ensure that the reference dictionary is a superset of the variant dictionary:
            checkReferenceDictionaryIsSupersetOfVariantDictionary();
        }
        else {
            logger.info("Skipping sequence dictionary validation.");
        }

        logger.info("Processing user transcripts/defaults/overrides...");
        Utils.validateArg(funcotatorArgs.outputFormatType != FuncotatorArgumentDefinitions.OutputFormatType.SEG,
                "This tool does not support segment output.  Please see FuncotateSegments.");


        // Next set up our transcript list:
        final Set<String> finalUserTranscriptIdSet = FuncotatorEngine.processTranscriptList(funcotatorArgs.userTranscriptIdSet);

        // Get our overrides for annotations:
        final LinkedHashMap<String, String> annotationDefaultsMap = FuncotatorEngine.splitAnnotationArgsIntoMap(funcotatorArgs.annotationDefaults);
        final LinkedHashMap<String, String> annotationOverridesMap = FuncotatorEngine.splitAnnotationArgsIntoMap(funcotatorArgs.annotationOverrides);

        // Get the header for our variants:
        final VCFHeader vcfHeader = getHeaderForVariants();

        if (!funcotatorArgs.reannotateVCF) {
            checkIfAlreadyAnnotated(vcfHeader, drivingVariantFile);
        }

        logger.info("Initializing data sources...");
        // Initialize all of our data sources:
        // Sort data sources to make them process in the same order each time:
        funcotatorArgs.dataSourceDirectories.sort(Comparator.naturalOrder());
        final Map<Path, Properties> configData = DataSourceUtils.getAndValidateDataSourcesFromPaths(funcotatorArgs.referenceVersion, funcotatorArgs.dataSourceDirectories);

        logger.info("Finalizing data sources (this step can be long if data sources are cloud-based)...");
        // Create the data sources from the input:
        // This will also create and register the FeatureInputs (created by the Data Sources)
        // with the GATK Engine, so we do not have to plumb them in after the fact.
        final List<DataSourceFuncotationFactory> dataSourceFuncotationFactories = DataSourceUtils.createDataSourceFuncotationFactoriesForDataSources(
                configData,
                annotationOverridesMap,
                funcotatorArgs.transcriptSelectionMode,
                finalUserTranscriptIdSet,
                this,
                funcotatorArgs.lookaheadFeatureCachingInBp,
                new FlankSettings(funcotatorArgs.fivePrimeFlankSize, funcotatorArgs.threePrimeFlankSize),
                false,
                funcotatorArgs.minNumBasesForValidSegment,
                funcotatorArgs.spliceSiteWindow
        );

        logger.info("Initializing Funcotator Engine...");
        // Create our engine to do our work and drive this Funcotation train!
        funcotatorEngine = new FuncotatorEngine(
                funcotatorArgs,
                getSequenceDictionaryForDrivingVariants(),
                VcfFuncotationMetadata.create(
                    new ArrayList<>(vcfHeader.getInfoHeaderLines())
                ),
                dataSourceFuncotationFactories
        );

        // Create our output renderer:
        logger.info("Creating a " + funcotatorArgs.outputFormatType + " file for output: " + funcotatorArgs.outputFile.toURI());
        outputRenderer = funcotatorEngine.createOutputRenderer(
                annotationDefaultsMap,
                annotationOverridesMap,
                vcfHeader,
                getDefaultToolVCFHeaderLines(),
                this
        );
    }

    /**
     *  Checks to see if the given vcf has already been annotated.
     *
     *  No need to annotate again.
     */
    @VisibleForTesting
    static void checkIfAlreadyAnnotated(final VCFHeader vcfHeader, GATKPath drivingVariantFile) {
        if (vcfHeader.getOtherHeaderLine(FuncotatorConstants.FUNCOTATOR_VERSION_VCF_HEADERLINE_KEY) != null) {
            throw new UserException.BadInput("Given VCF " +drivingVariantFile+ " has already been annotated!");
        }
    }

    /**
     * Checks to see that the given reference's sequence dictionary is a
     * superset of the given variant file's dictionary.
     *
     * This is a more strict check than the one found in {@link GATKTool#validateSequenceDictionaries()}.
     */
    private void checkReferenceDictionaryIsSupersetOfVariantDictionary() {

        final SAMSequenceDictionary referenceDictionary = getReferenceDictionary();
        final SAMSequenceDictionary variantDictionary = getSequenceDictionaryForDrivingVariants();

        if ( referenceDictionary == null ) {
            throw new UserException.BadInput("Reference fasta sequence dictionary is null!");
        }

        if ( variantDictionary == null ) {
            throw new UserException.BadInput("Funcotator by default requires that the variant input have a sequence dictionary in its header. To disable this safety check, use argument --" + StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME);
        }

        SequenceDictionaryUtils.validateDictionaries(
                "Reference", getReferenceDictionary(),
                "Driving Variants", getSequenceDictionaryForDrivingVariants(),
                true,
                false
                );
    }

    @Override
    protected CountingVariantFilter makeVariantFilter() {
        return new CountingVariantFilter(FuncotatorEngine.makeVariantFilter(funcotatorArgs.removeFilteredVariants));
    }

    @Override
    public VariantTransformer makePostVariantFilterTransformer(){
        return funcotatorEngine.getDefaultVariantTransformer();
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        // Get the correct reference for B37/HG19 compliance:
        // This is necessary because of the variant transformation that gets applied in VariantWalkerBase::apply.
        final ReferenceContext correctReferenceContext = funcotatorEngine.getCorrectReferenceContext(variant, referenceContext);
        // Place the variant on our queue to be funcotated:
        enqueueAndHandleVariant(variant, correctReferenceContext, featureContext);
    }

    @Override
    public Object onTraversalSuccess() {

        // If we only saw IGRs, we most likely have a configuration issue.
        // Make sure the user knows this by making a HUGE stink about it.
        if ( funcotatorEngine.onlyProducedIGRs() ) {
            logger.warn("================================================================================");
            logger.warn("\u001B[43m     _  _  _   __        __               _                   _  _  _           ");
            logger.warn("    | || || |  \\ \\      / /_ _ _ __ _ __ (_)_ __   __ _      | || || |        ");
            logger.warn("    | || || |   \\ \\ /\\ / / _` | '__| '_ \\| | '_ \\ / _` |     | || || |     ");
            logger.warn("    |_||_||_|    \\ \\V V / (_| | |  | | | | | | | | (_| |     |_||_||_|        ");
            logger.warn("    (_)(_)(_)     \\_/\\_/ \\__,_|_|  |_| |_|_|_| |_|\\__, |     (_)(_)(_)      ");
            logger.warn("                                                  |___/                         \u001B[0;0m");
            logger.warn("--------------------------------------------------------------------------------");
            logger.warn(" Only IGRs were produced for this dataset.  This STRONGLY indicates that this   ");
            logger.warn(" run was misconfigured.     ");
            logger.warn(" You MUST check your data sources to make sure they are correct for these data.");
            logger.warn("================================================================================");
        }
        return true;
    }

    @Override
    public void closeTool() {
        if ( funcotatorEngine != null) {
            funcotatorEngine.close();
        }

        if ( outputRenderer != null ) {
            outputRenderer.close();
        }
    }

    //==================================================================================================================

    /**
     * Creates an annotation on the given {@code variant} or enqueues it to be processed during a later call to this method.
     * @param variant {@link VariantContext} to annotate.
     * @param referenceContext {@link ReferenceContext} corresponding to the given {@code variant}.
     * @param featureContext {@link FeatureContext} corresponding to the given {@code variant}.
     */
    protected void enqueueAndHandleVariant(final VariantContext variant, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        final FuncotationMap funcotationMap = funcotatorEngine.createFuncotationMapForVariant(variant, referenceContext, featureContext);

        // This is necessary because we want to revert the variant contig name change if it was applied in the FuncotatorEngine::getCorrectVariantContextForReference method before outputting the vcf.
        // NOTE: this will only revert the variantContext if it was originally changed (only for B37 VCFs)
        final VariantContext variantContextForOutput = funcotatorEngine.getCorrectVariantContextForOutput(variant);

        // At this point there is only one transcript ID in the funcotation map if canonical or best effect are selected
        outputRenderer.write(variantContextForOutput, funcotationMap);
    }
}
