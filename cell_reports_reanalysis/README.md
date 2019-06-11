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
 * STAR 2.7.0f - version used in paper was 2.2.4  
 * SAMtools 1.9 - version not listed in paper
 
 * SRA-Toolkit 2.9.6 
 * TEToolkit 2.0.3  


## RNA-seq Pipeline

All RNA-Seq data processing steps were done on NIH's Biowulf cluster.

### 1. Download RNA-Seq data

Load SRA-Toolkit    
`module load sratoolkit/2.9.6`  

Run fastq-dump on SRA paired-end files    
`fastq-dump -I --split-files SRR8571937`


### 2. Perform QC with FastQC

Per Dr. Edward Lee, trimming and filtering of reads were not performed. 

### 3. Map reads to reference genome with STAR

In the paper, sequencing reads were aligned against GRCh38 (GENCODE release 22) using STAR 2.2.4 with option --outFilterIntronMotifs RemoveNonCanonical. All other parameters were kept as default.

For this reanalysis, we aligned reads against GRCh38 (GENCODE release 28) instead since we already have been using this version for other analyses.

### 4. Filter reads

Liu et al. kept only uniquely mapping reads and filtered out ribosomal, mitochondrial, and non-standard chromosomal reads (i.e. X or Y chromosome reads). The authors modified the script called *norm_scripts/filter_sam.pl* from the PORT pipeline https://github.com/itmat/Normalization.git. 

### 5. Convert SAM files to BAM files and sort by coordinate 

`samtools view`

`samtools sort`

### 6. 
