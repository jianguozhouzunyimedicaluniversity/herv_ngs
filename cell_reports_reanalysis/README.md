# Reanalysis of FTD-ALS RNA sequencing data

Replicate the paper's processing of RNA-Seq data, use TEToolkit to identify HERVs, and compare results with the paper. 

## Paper

Liu, E. Y., Russ, J., Cali, C. P., Phan, J. M., Amlie-Wolf, A., & Lee, E. B. (2019). Loss of Nuclear TDP-43 Is Associated with Decondensation of LINE Retrotransposons. Cell reports, 27(5), 1409–1421.e6. doi:10.1016/j.celrep.2019.04.003

## Data

RNA-Seq data were deposited in the Gene Expression Omnibus (GEO) database with accession number GSE126543.

The list of SRA accessions numbers (14 total) is included in this folder. 

RNA was extracted from the neuronal nuclei of the frontal cortex from 7 patients. Each patient had two samples, one TDP-43 positive and one TDP-43 negative. 

## Software

 * FastQC 0.11.8 - version used in paper was 0.11.3     
 * STAR 2.6.1c - version used in paper was 2.2.4  
 * SAMtools 1.9 - version not listed in paper
 * PORT v0.8.5
 
 * SRA-Toolkit 2.9.6 
 * TEToolkit 2.0.3  


## RNA-seq Pipeline

All RNA-Seq data processing steps were done on NIH's Biowulf cluster.

### 1. Download RNA-Seq data

__Load SRA-Toolkit__

`module load sratoolkit/2.9.6`    

__Run fastq-dump on SRA paired-end files__   

`fastq-dump -I --split-files SRR8571937`  


### 2. Perform QC with FastQC

__Load FastQC__    

`module load fastqc/0.11.8`  

__Run FastQC__  

`fastqc -o fastqc_reports SRR8571937_1.fastq`  

Per Dr. Edward Lee, trimming and filtering of reads were not performed. 

### 3. Map reads to reference genome with STAR

In the paper, sequencing reads were aligned against GRCh38 (GENCODE release 22) using STAR 2.2.4 with option --outFilterIntronMotifs RemoveNonCanonical. All other parameters were kept as default.

For this reanalysis, we aligned reads against GRCh38 (GENCODE release 28) instead since we already have been using this version for other analyses. In addition, the GTF file from this release was concatenated with Dr. Nath's lab's GTF file of HERVK annotations (*ALS_Annotations.gtf*). As a result, our genome indices for STAR was slightly different from the authors' genome indices. 

We added the parameter --outSAMunmapped Within KeepPairs to ensure that the data was compatible with PORT since we were using STAR >= v2.5.1a. 

We did consider using STAR 2.2.4 to reproduce the results as closely as possible, however this version was not available in Alex Dobin's STAR github, https://github.com/alexdobin/STAR. 

__Load STAR__  

`module load STAR/2.6.1c`

__Run STAR__    

<pre>
STAR --runThreadN $SLURM_CPUS_PER_TASK --genomeDir /data/ALS_Working_Grp/Star/indices/hg38 --readFilesIn SRR8571939_1.fastq SRR8571939_2.fastq --outFileNamePrefix /data/ALS_Working_Grp/Cell_Reports_reanalysis/sam/SRR8571939 --outFilterIntronMotifs RemoveNoncanonical --outSAMunmapped Within KeepPairs </pre>


Note: The most recent version of STAR, 2.7.0f, had the following error: *Genome version: 20201 is INCOMPATIBLE with running STAR version: 2.7.0f*. As a work around, we used an older version to avoid having to re-generate the genome indices from scratch. 

### 4. Filter reads with PORT

Liu et al. kept only uniquely mapping reads and filtered out ribosomal, mitochondrial, and non-standard chromosomal reads (i.e. X or Y chromosome reads). The authors modified the script called *norm_scripts/filter_sam.pl* from the PORT pipeline https://github.com/itmat/Normalization.git. 

__Download PORT__  

`git clone https://github.com/itmat/Normalization.git`

__Reorganize file directories for PORT__

The input files, FASTQ and SAM, must be organized into a specific directory structure for PORT to run properly. 

Example:

<pre>
STUDY
└── READS
    ├── Sample_1
    │   ├── Unaligned reads
    │   └── Aligned.sam/bam
    ├── Sample_2
    │   ├── Unaligned reads
    │   └── Aligned.sam/bam
    ├── Sample_3
    │   ├── Unaligned reads
    │   └── Aligned.sam/bam
    └── Sample_4
        ├── Unaligned reads
        └── Aligned.sam/bam
</pre>

__Run BLAST to find ribosomal reads in FASTQ files__

<pre> perl runblast.pl <dir> <loc> <blastdir> <query> [option] rRNA_mm9.fa </pre>

-fq: set this if the unaligned files are in fastq format
-pe \"<unlaligned1>,<unaligned2>\" : set this if the data are paired end and provide two unaligned files
 
 The developers for PORT provided a FASTA file of ribosomal RNAs that can be used for mammals, *norm_scripts/rRNA_mm9.fa*.
 
__Parse BLAST output to extract ribosomal IDs__

<pre>perl parseblastout.pl <id> <loc> </pre>

__Filter SAM file__


<pre>perl PORT/norm_scripts/filter_sam.pl SRR8571939.sam SRR8571939.filtered.sam ribosomalids.txt -u</pre>

The script was modified to filter out reads mapped to chromosomes X and Y.  

This script will remove all rows from the SAM file except those that satisfy all of the following:  
1. Unique mapper (-u option)  
2. Both forward and reverse reads map consistently  
3. Read ID not in the ribosomalids text file  
4. Chromosome is one of the numbered ones (e.g. chr1, chr2, OR 1, 2)
5. Is a forward mapper (script outputs forward mappers only)


### 5. Sort by coordinate 

__Load SAMtools__

`module load samtools/1.9`

`samtools view`

`samtools sort`

`samtools idxstats alignments.bam`

### 6. 
