#!/bin/bash 

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --nodes=1
#SBATCH --job-name=cutadapt
#SBATCH --error=cutadapt.err
#SBATCH --output=cutadapt.out
#SBATCH --time=0-12:00:00
#SBATCH --cpus-per-task=8

cd /projects/bgmp/wrosales/bioinfo/Bi623/QAA/
conda activate QAA
/usr/bin/time -v cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -o 11_2H_both_S9_L008_trim_outputs/11_2H_both_S9_L008_R1_cutadapt_output.fastq.gz /projects/bgmp/shared/2017_sequencing/demultiplexed/11_2H_both_S9_L008_R1_001.fastq.gz
/usr/bin/time -v cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o 11_2H_both_S9_L008_trim_outputs/11_2H_both_S9_L008_R2_cutadapt_output.fastq.gz /projects/bgmp/shared/2017_sequencing/demultiplexed/11_2H_both_S9_L008_R2_001.fastq.gz
/usr/bin/time -v cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -o 29_4E_fox_S21_L008_trim_outputs/29_4E_fox_S21_L008_R1_cutadapt_output.fastq.gz /projects/bgmp/shared/2017_sequencing/demultiplexed/29_4E_fox_S21_L008_R1_001.fastq.gz
/usr/bin/time -v cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o 29_4E_fox_S21_L008_trim_outputs/29_4E_fox_S21_L008_R2_cutadapt_output.fastq.gz /projects/bgmp/shared/2017_sequencing/demultiplexed/29_4E_fox_S21_L008_R2_001.fastq.gz
