# Reanalysis of FTD-ALS RNA sequencing data

Replicate the paper's processing of RNA-Seq data, use TEToolkit to identify HERVs, and compare results with the paper. 

## Paper

Liu, E. Y., Russ, J., Cali, C. P., Phan, J. M., Amlie-Wolf, A., & Lee, E. B. (2019). Loss of Nuclear TDP-43 Is Associated with Decondensation of LINE Retrotransposons. Cell reports, 27(5), 1409–1421.e6. doi:10.1016/j.celrep.2019.04.003

## Data

RNA-Seq data was deposited in the Gene Expression Omnibus (GEO) database with accession number GSE126543.

## Software

FastQC 0.11.3 (0.11.8)
STAR 2.2.4 (2.7.0f)
TEToolkit 2.0.3
SAMtools 1.9


## RNA-seq Pipeline

All RNA-Seq data processing steps were done on NIH's Biowulf cluster.

### 1. Download RNA-Seq data

### 2. Perform QC with FastQC

Per Dr. Edward Lee, trimming and filtering of reads were not performed. 

### 3. Map reads to reference genome with STAR

In the paper, sequencing reads were aligned against GRCh38 (GENCODE release 22) using STAR 2.2.4 with option --outFilterIntronMotifs RemoveNonCanonical. All other parameters were kept at default.

For this reanalysis, we aligned reads against GRCh38 (GENCODE release 28) instead since we already have been using this version for other analyses.

### 4. Filter reads

Liu et al. kept only uniquely mapping reads and filtered out ribosomal, mitochondrial, and non-standard chromosomal reads (i.e. X or Y chromosome reads). The authors modified the script called *norm_scripts/filter_sam.pl* from the PORT pipeline https://github.com/itmat/Normalization.git. 

norm_scripts/filter_sam.pl

### 5. Convert SAM files to BAM files and sort by coordinate 

`samtools view`

`samtools sort`

### 6. 
