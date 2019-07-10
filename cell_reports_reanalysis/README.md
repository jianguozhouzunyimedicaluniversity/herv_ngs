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
 * PORT v0.8.5 - version not listed in paper
 * BLAST 2.2.30+ - comes with PORT installation
 
 * SRA-Toolkit 2.9.6 
 * TEToolkit 2.0.3  


## RNA-seq Pipeline

All RNA-Seq data processing steps were done on NIH's Biowulf cluster.

An interactive session was created using the following command:  
`sinteractive --cpus-per-task=12 --mem=100g --gres=lscratch:100`  


### 1. Download RNA-Seq data

__Load SRA-Toolkit__

`module load sratoolkit/2.9.6`    

__Run fastq-dump on SRA paired-end files__   

Example:  
`fastq-dump -I --split-files SRR8571937`  


### 2. Perform QC with FastQC

__Load FastQC__    

`module load fastqc/0.11.8`  

__Run FastQC__  

Example:  
`fastqc -o /data/ALS_Working_Grp/Cell_Reports_reanalysis/seqs/fastqc_reports
/data/ALS_Working_Grp/Cell_Reports_reanalysis/seqs/SRR8571937_1.fastq` 

Per Dr. Edward Lee, trimming and filtering of reads were not performed. 

### 3. Map reads to reference genome with STAR

In the paper, sequencing reads were aligned against GRCh38 (GENCODE release 22) using STAR 2.2.4 with option --outFilterIntronMotifs RemoveNonCanonical. All other parameters were kept as default.

For this reanalysis, we aligned reads against GRCh38 (GENCODE release 28) instead since we already have been using this version for other analyses. In addition, the GTF file from this release was concatenated with Dr. Nath's lab's GTF file of HERVK annotations (*ALS_Annotations.gtf*). As a result, our genome indices for STAR was slightly different from the authors' genome indices. 

We added the parameter --outSAMunmapped Within KeepPairs to ensure that the data was compatible with PORT since we were using STAR >= v2.5.1a. 

We did consider using STAR 2.2.4 to reproduce the results as closely as possible, however this version was not available in Alex Dobin's STAR github, https://github.com/alexdobin/STAR. 

__Load STAR__  

`module load STAR/2.6.1c`

__Run STAR__    

Example:  
<pre>
STAR --runThreadN $SLURM_CPUS_PER_TASK --genomeDir /data/ALS_Working_Grp/Star/indices/hg38   
--readFilesIn /data/ALS_Working_Grp/Cell_Reports_reanalysis/seqs/SRR8571937_1.fastq  
/data/ALS_Working_Grp/Cell_Reports_reanalysis/seqs/SRR8571937_2.fastq 
--outFileNamePrefix /data/ALS_Working_Grp/Cell_Reports_reanalysis/sam/SRR8571937
--outFilterIntronMotifs RemoveNoncanonical --outSAMunmapped Within KeepPairs 
</pre>


Note: The most recent version of STAR, 2.7.0f, had the following error: *Genome version: 20201 is INCOMPATIBLE with running STAR version: 2.7.0f*. As a work around, we used an older version to avoid having to re-generate the genome indices from scratch. 

### 4. Filter reads with PORT

Liu et al. kept only uniquely mapping reads and filtered out ribosomal, mitochondrial, and non-standard chromosomal reads (i.e. X or Y chromosome reads). The authors modified the script called *norm_scripts/filter_sam.pl* from the PORT pipeline https://github.com/itmat/Normalization.git. 

__Download PORT__  

`git clone https://github.com/itmat/Normalization.git`

__Reorganize file directories for PORT__

The input files, FASTQ and SAM, must be organized into a specific directory structure for PORT to run properly. All alignment SAM files must have the same name across samples.

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

__Convert FASTQ files to FASTA__

Used PORT's *norm_scripts/fastq2fasta.pl* script to convert all FASTQ files to FASTA files. 

Example:  
<pre>
perl /data/ALS_Working_Grp/Cell_Reports_reanalysis/PORT/norm_scripts/fastq2fasta.pl   
/data/ALS_Working_Grp/Cell_Reports_reanalysis/reads/SRR8571937/SRR8571937_1.fastq   
/data/ALS_Working_Grp/Cell_Reports_reanalysis/reads/SRR8571937/SRR8571937_1.fasta  
</pre>

__Load BLAST__ 

`module load blast/2.2.30+`

__Build BLAST databases with FASTA files__

Example:  

Database for sample's forward reads:   
<pre>
makeblastdb -dbtype nucl -max_file_sz 300MB -in 
/data/ALS_Working_Grp/Cell_Reports_reanalysis/reads/SRR8571937/SRR8571937_1.fasta -out 
/data/ALS_Working_Grp/Cell_Reports_reanalysis/reads/SRR8571937/blastdb1.SRR8571937
</pre>

Database for sample's reverse reads:
<pre>
makeblastdb -dbtype nucl -max_file_sz 300MB -in
/data/ALS_Working_Grp/Cell_Reports_reanalysis/reads/SRR8571937/SRR8571937_2.fasta 
-out /data/ALS_Working_Grp/Cell_Reports_reanalysis/reads/SRR8571937/blastdb2.SRR8571937
</pre>

Note: We tried to use *norm_scripts/runall_runblast.pl* but received an error when using cluster option -other and Biowulf's sbatch parameters. 

__Identify ribosomal reads in BLAST databases__

Example:  
<pre>
blastn -task blastn -db 
/data/ALS_Working_Grp/Cell_Reports_reanalysis/reads/SRR8571937/blastdb1.SRR8571937 
-query /data/ALS_Working_Grp/Cell_Reports_reanalysis/PORT/norm_scripts/rRNA_mm9.fa 
-num_descriptions 1000000000 -num_alignments 1000000000 -num_threads $SLURM_CPUS_PER_TASK > 
/data/ALS_Working_Grp/Cell_Reports_reanalysis/reads/SRR8571937/blastdb1.SRR8571937_blastout
</pre>
 
 The developers for PORT provided a FASTA file of ribosomal RNAs that can be used for mammals, *norm_scripts/rRNA_mm9.fa*.
 
__Parse BLAST output to extract ribosomal IDs__

Example:
<pre>
perl /data/ALS_Working_Grp/Cell_Reports_reanalysis/PORT/norm_scripts/parseblastout.pl SRR8571937 
/data/ALS_Working_Grp/Cell_Reports_reanalysis/reads
</pre>

The *norm_scripts/parseblastout.pl* script was modified so that the BLAST databases and output are not deleted.

__Filter SAM file__

<pre>
perl /data/ALS_Working_Grp/Cell_Reports_reanalysis/PORT/norm_scripts/filter_sam.pl 
/data/ALS_Working_Grp/Cell_Reports_reanalysis/reads/SRR8571937/Aligned.sam 
/data/ALS_Working_Grp/Cell_Reports_reanalysis/filtered_sam/SRR8571937.sam 
/data/ALS_Working_Grp/Cell_Reports_reanalysis/reads/SRR8571937/SRR8571937.ribosomalids.txt -u</pre>

The *norm_scripts/filter_sam.pl* script was modified to filter out reads mapped to chromosomes X and Y. 
In addition, the script originally outputted forward mappers only. 
However, TEcount expects both forward and reverse reads since it knows that reads are paired. 
With just the forward mappers only, I got the following error with TEcount:

*If the BAM file is sorted by coordinates, please specify --sortByPos and re-run!*

This modified script will remove all rows from the SAM file except those that satisfy all of the following:  
1. Unique mapper (-u option)  
2. Both forward and reverse reads map consistently  
3. Read ID not in the ribosomalids text file  
4. Chromosome is one of the numbered ones (e.g. chr1, chr2, OR 1, 2)

__Gunzip SAM files (output)__
<pre> gunzip /data/ALS_Working_Grp/Cell_Reports_reanalysis/filtered_sam/SRR8571937_u.sam.gz </pre>

### 5. Convert SAM to BAM format and sort by coordinate 

__Load SAMtools__

`module load samtools/1.9`

__Run SAMtools faidx__

Index the reference sequence in the FASTA format. 

This step is required to add in the @SQ lines in the header of the SAM files. This header information was removed during the 
filtering SAM file step of the PORT pipeline.

<pre>
samtools faidx /data/ALS_Working_Grp/Reference/hg38.fa
</pre>

__Run SAMtools view to generate SAM file with header info__

Example:  
<pre>
samtools view -ht /data/ALS_Working_Grp//Reference/hg38.fa.fai 
/data/ALS_Working_Grp/Cell_Reports_reanalysis/filtered_sam/SRR8571937_u.sam > 
/data/ALS_Working_Grp/Cell_Reports_reanalysis/filtered_sam/SRR8571937_u_header.sam
</pre>

__Run SAMtools to convert SAM file to BAM file__

Example:  
<pre>
samtools view -S -b -o /data/ALS_Working_Grp/Cell_Reports_reanalysis/bam/SRR8571937.bam  
/data/ALS_Working_Grp/Cell_Reports_reanalysis/filtered_sam/SRR8571937_u_header.sam
</pre>


__Run SAMtools sort__

Example:  
<pre>
samtools sort -T /lscratch/$SLURM_JOB_ID/SRR8571937 
/data/ALS_Working_Grp/Cell_Reports_reanalysis/bam/SRR8571937.bam -O BAM -o
/data/ALS_Working_Grp/Cell_Reports_reanalysis/sorted_bam/SRR8571937_coordsort.bam
</pre>

### 6. Annotate reads to transposable elements with TEcount

__Load TEToolKit__

`module load tetoolkit/2.0.3`

__Run TEcount__

According to Tara, the PI told her that the reads were not stranded, so we selected the option "no" for the "stranded" argument.

Example:  
<pre>
TEcount --sortByPos --format BAM --mode multi --stranded no -i 100
-b /data/ALS_Working_Grp/Cell_Reports_reanalysis/sorted_bam/SRR8571937_coordsort.bam 
--GTF /data/ALS_Working_Grp/Gtf/hg38.gtf --TE /data/ALS_Working_Grp/Gtf/HERVK_Nath_2.gtf
--project /data/ALS_Working_Grp/Cell_Reports_reanalysis/tecounts/SRR8571937_tecounts
</pre>

Required and optional arguments:  
--GTF genic-GTF-file: GTF file for gene annotations    
--TE TE-GTF-file: GTF file for transposable element annotations    
--mode [TE counting mode]: uniq (unique mappers only) or multi (distribute among all alignments). DEFAULT: multi  
--format [input file format]: Input file format: BAM or SAM. DEFAULT: BAM  
--stranded [option] Is this a stranded library? (yes, no, or reverse). DEFAULT: yes  
--sortByPos: Input file is sorted by chromosome position  
--project [name]: Prefix used for output files (e.g. project name) DEFAULT: TEcount_out  
-i | --iteration: maximum number of iterations used to optimize multi-reads assignment. DEFAULT: 100  
-b: BAM file

For more details on how to use TEcount, go to https://github.com/mhammell-laboratory/tetoolkit.



`samtools flagstat alignments.bam`
