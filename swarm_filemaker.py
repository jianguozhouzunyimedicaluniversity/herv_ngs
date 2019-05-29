#!/usr/bin/env python3

import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
import os
from os.path import join
from glob import glob
import sys
import argparse

'''
Purpose: Generates swarm and other files for each step in the HERV NGS pipeline.
Reads in an excel file and obtains sample IDs from specified sheet/column.

To use:
python swarm_filemaker.py <sample_info EXCELFILE> <sheet STRING> <sample_ids STRING> [options]

Example:
python swarm_filemaker.py Updated_Sample_Info_v20180315_TDO.xlsx \
        Sample_Info \
        NYGC_DataFileID \

'''
__author__ = 'Sara Jones'
__author_email__ = 'jonessarae@gmail.com'
__doc__ = 'Generates swarm and other files for deployment on Biowulf cluster.'
__date_modified__ = '5/29/19'

################ path to files, change according to project ###################

# path to STAR indices for hg38
genome_dir = '/data/johnsonko/Tara/star_hg38/indices/hg38'
# path to directory containing bam files
#bam_dir = '/data/johnsonko/Tara/bam/'
bam_dir = '/Users/jonesse3/Desktop'
# path to directory containing fastq files
#seq_dir = '/data/johnsonko/Tara/Seqs/'
seq_dir = '/Users/jonesse3/Desktop'
# path to gtf file for hg38 gene annotations
hg38_gtf = '/data/ALS_Working_Grp/Gtf/hg38.gtf'
# path to gtf file for transposable element annotations
te_gtf = '/data/ALS_Working_Grp/Gtf/HERVK_Nath_2.gtf'

def star_swarm(seq_dir, sample_list, genome_dir, bam_dir):
    '''
    Writes a swarm file for mapping trimmed reads to a reference genome using
    STAR and uses samtools to generate a bam file.

    Arguments:
        seq_dir: path to directory containing trimmed fastq files
        sample_list: list of sample IDs
        genome_dir: path to directory containing STAR indices
        bam_dir: path to store bam files
    '''
    # open and write to file
    f = open('star.swarm', 'w')

    # loop through list of sample IDs
    for id in sample_list:

        # list to hold file names
        files = []

        # loop through directory to find files with 1.fastq and 2.fastq suffices
        for ext in ('*1.fastq', '*2.fastq'):
            files.extend(glob(join(seq_dir, id + ext)))

        # check that there are two files for the sample ID
        if len(files) == 0:
            print('\nMissing both R1 and R2 fastq files for {}.'.format(id))
        elif len(files) == 1:
            print('\nMissing a file for {}. Only the following was found:\n'.format(id))
            for file in files:
                print(file)
        elif len(files) > 2:
            print('\nExpecting two files for {}.'.format(id))
            print('\nThe following files were found:\n')
            for file in files:
                print(file)
        else:
            # check that we have R1 and R2 files
            if files[0].endswith('1.fastq'):
                r1 = files[0] # R1 file of paired-end run
            else:
                print('\nMissing R1 file for {}.\n'.format(id))
                continue
            if files[1].endswith('2.fastq'):
                r2 = files[1] # R2 file of paired-end run
            else:
                print('\nMissing R2 file for {}.\n'.format(id))
                continue
            # write STAR command to file
            f.write('STAR --runThreadN $SLURM_CPUS_PER_TASK '+
                '--genomeDir {} '.format(genome_dir) +
                '--sjdbOverhang 100 ' +
                '--readFilesIn {} {} '.format(r1,r2) +
                '--outFilterMultimapNmax 100 ' +
                '--winAnchorMultimapNmax 150 ' +
                '--genomeLoad LoadAndRemove ' +
                '--outSAMattributes Standard ' +
                '--outSAMunmapped None ' +
                '--outFilterType BySJout ' +
                '--outStd SAM ' +
                '--outSAMstrandField intronMotif ' +
                '--alignSJoverhangMin 8 ' +
                '--alignSJDBoverhangMin 1 ' +
                '--outFileNamePrefix {} '.format(id) +
                '--outTmpDir=/lscratch/$SLURM_JOB_ID/STARtmp | ' +
                'samtools view -S -b -o {}unsorted/{}.bam\n'.format(bam_dir, id))

    # close file
    f.close()

def samtools_sort(sample_list, bam_dir):
    '''
    Writes a swarm file that sorts a bam file by read name using samtools.

    Arguments:
        sample_list: list of sample IDs
        bam_dir: path to bam files
    '''
    # open and write to file
    f = open('ReadNameSort.swarm', 'w')

    # loop through list of sample IDs
    for id in sample_list:

        # loop through directory to find files
        file = glob(join(bam_dir, id + '.bam'))

        # check that there is a bam file for sample ID
        if not file:
            print('\nMissing bam file for {}.\n'.format(id))
            continue
        else:
            # write samtools command to file
            f.write('samtools sort -n -T '+
                '/lscratch/$SLURM_JOB_ID/{} '.format(id) +
                '{}.bam '.format(id) +
                '-O BAM -o {}/sorted/READ_NAME_SORTED_{}.bam\n'.format(bam_dir, id))

    # close file
    f.close()

def te_count(sample_list, bam_dir, hg38_gtf, te_gtf):
    '''
    Writes a swarm file that generates a count table of transposable elements
    using TEcounts.

    Arguments:
        sample_list: list of sample IDs
        bam_dir: path to bam files
        hg38_gtf: path to hg38 annotation file
        te_gtf: path to transposable element annotation file
    '''
    # open and write to file
    f = open('TE_count.swarm', 'w')

    # loop through list of sample IDs
    for id in sample_list:

        # loop through directory to find files
        file = glob(join(bam_dir, 'READ_NAME_SORTED_' + id + '.bam'))

        # check that there is a bam file for sample ID
        if not file:
            print('\nMissing name-sorted bam file for {}.\n'.format(id))
            continue
        else:
            # write TEcount command to file
            f.write('TEcount ' +
                '-b {}/sorted/READ_NAME_SORTED_{}.bam '.format(bam_dir, id) +
                '--GTF {} '.format(hg38_gtf) +
                '--TE {} '.format(te_gtf) +
                '--format BAM ' +
                '--stranded reverse ' +
                '--mode multi ' +
                '-i 100 ' +
                '--project {}/sorted/TE_COUNTS_{}\n'.format(bam_dir, id))

    # close file
    f.close()

def herv_filter(sample_list, bam_dir):
    '''
    Writes a bash script that filters out HERVs using the pattern "dup"
    in the count table generated from TEcount.

    Arguments:
        sample_list: list of sample IDs
        bam_dir: path to bam files
    '''
    # open and write to file
    f = open('herv_filter.sh', 'w')
    f.write('#!/bin/bash\n')

    # loop through list of sample IDs
    for id in sample_list:

        # loop through directory to find files
        file = glob(join(bam_dir, 'TE_COUNTS_' + id + '.cntTable'))

        # check that there is a bam file for sample ID
        if not file:
            print('\nMissing count table file for {}.\n'.format(id))
            continue
        else:
            # write linux command to file
            f.write('more {}/sorted/TE_COUNTS_{}.cntTable | '.format(bam_dir, id) +
                'grep dup > {}/sorted/HERV_COUNTS_{}.txt &\n'.format(bam_dir, id))

    # close file
    f.close()

def main(args):

    # path to excel file containing sample info
    sample_file = args.sample_info
    # name of sheet with sample info
    sheetname = args.sheet
    # column name with sample ID
    sample_id_col = args.sample_ids

    # create dataframe of excel file
    df = pd.read_excel(sample_file, sheet_name = sheetname)

    # obtain list of samples
    sample_list = df[sample_id_col]

    # create star_swarm file
    star_swarm(seq_dir, sample_list, genome_dir, bam_dir)

    # create ReadNameSort.swarm file
    samtools_sort(sample_list, bam_dir)

    # create TE_count.swarm file
    te_count(sample_list, bam_dir, hg38_gtf, te_gtf)

    # create herv_filter.sh file
    herv_filter(sample_list, bam_dir)


if __name__ == "__main__":

    # Create arguments
    p = argparse.ArgumentParser(description=__doc__,\
                                prog='swarm_filemaker.py', \
                                usage='%(prog)s ' +
                                '<sample_info FILE> ' +
                                '<sheet STRING> ' +
                                '<sample_ids STRING> ' +
                                '[options]', add_help=True)
    p.add_argument('sample_info', help='path to excel file with sample info')
    p.add_argument('sheet', help='name of excel sheet with sample info')
    p.add_argument('sample_ids', help='column name containing sample IDs')
    #p.add_argument('--top_k', default=5, help='number of the top classes to print out, default is 5', type=int)
    #p.add_argument('--gpu', action='store_true', default=False, help='enable gpu')

    # Check number of arguments entered by user
    if len(sys.argv) < 4:
        print('Please enter required info in <>.\n')
        p.print_help()
        sys.exit(0)
    else:
        # Set arguments
        args = p.parse_args()
        print(args)
        main(args)
