#!/usr/bin/env python3

import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
import os
from os.path import join
from glob import glob
import sys

genome_dir = '/data/johnsonko/Tara/star_hg38/indices/hg38'
bam_dir = 'data/johnsonko/Tara/bam/'

# path to excel file containing sample info
sample_file = '/Users/Sara/Downloads/Updated_Sample_Info_v20180315_TDO.xlsx'
# name of sheet with sample info
sheetname = 'Sample_Info'
# column name with sample ID
sample_id_col = 'NYGC_DataFileID'
# path to directory containing fastq files
seq_dir = os.path.realpath('/Users/Sara/Desktop')

#seq_dir = '/Users/Sara/Desktop/'
# create dataframe of excel file
df = pd.read_excel(sample_file, sheet_name = sheetname)

# list of samples
sample_list = df[sample_id_col]

files = []
missed_samples = []
for ext in ('*1.fastq', '*2.fastq'):
    files.extend(glob(join(seq_dir,sample_list[2] + ext)))

# check that there are two files for the sample ID
if len(files) == 0:
    print('Missing both R1 and R2 fastq files for {}.'.format(sample_list[2]))
elif len(files) == 1:
    print('Missing a file for {}. Only the following was found:\n.'.format(sample_list[2]))
    for file in files:
        print(file)
elif len(files) > 2:
    print('Sample {} has more than 2 files. They are:\n'.format(sample_list[1]))
    for file in files:
        print(file)
    print('\nThis sample will be skipped. Please check that you have the correct files.')
else:
    f = open('my_file.txt', 'w')
    r1 = files[0] # read 1 file of paired-end run
    r2 = files[1] # read 2 file of paired-end run
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
        '--outFileNamePrefix {} '.format(sample_list[0]) +
        '--outTmpDir=/lscratch/$SLURM_JOB_ID/STARtmp | ' +
        'samtools view -S -b -o {}{}.bam'.format(bam_dir, sample_list[0]))
    f.close()
'''
with open('my_path/my_file.txt', 'r') as f:
    file_data = f.read()

f = open('my_path/my_file.txt', 'w')
f.write("Hello there!")
f.close()

def create_cast_list(filename):
    cast_list = []
    with open(filename) as f:
        for line in f:
            name = line.split(",")[0]
            cast_list.append(name)

    return cast_list

cast_list = create_cast_list('flying_circus_cast.txt')
for actor in cast_list:
    print(actor)
'''
