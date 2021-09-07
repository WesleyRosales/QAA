#!/bin/bash 

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --nodes=1
#SBATCH --job-name=trimmomatic
#SBATCH --error=trimmomatic.err
#SBATCH --output=trimmomatic.out
#SBATCH --time=0-12:00:00
#SBATCH --cpus-per-task=8

cd /projects/bgmp/wrosales/bioinfo/Bi623/QAA/
conda activate QAA
trimmomatic PE -threads 8 -phred33 11_2H_both_S9_L008_trim_outputs/11_2H_both_S9_L008_R1_cutadapt_output.fastq.gz 11_2H_both_S9_L008_trim_outputs/11_2H_both_S9_L008_R2_cutadapt_output.fastq.gz \
11_2H_both_S9_L008_trim_outputs/11_2H_both_S9_L008_R1_paired.fastq.gz 11_2H_both_S9_L008_trim_outputs/11_2H_both_S9_L008_R1_unpaired.fastq.gz \
11_2H_both_S9_L008_trim_outputs/11_2H_both_S9_L008_R2_paired.fastq.gz 11_2H_both_S9_L008_trim_outputs/11_2H_both_S9_L008_R2_unpaired.fastq.gz \
LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:35
trimmomatic PE -threads 8 -phred33 29_4E_fox_S21_L008_trim_outputs/29_4E_fox_S21_L008_R1_cutadapt_output.fastq.gz 29_4E_fox_S21_L008_trim_outputs/29_4E_fox_S21_L008_R2_cutadapt_output.fastq.gz \
29_4E_fox_S21_L008_trim_outputs/29_4E_fox_S21_L008_R1_paired.fastq.gz 29_4E_fox_S21_L008_trim_outputs/29_4E_fox_S21_L008_R1_unpaired.fastq.gz \
29_4E_fox_S21_L008_trim_outputs/29_4E_fox_S21_L008_R2_paired.fastq.gz 29_4E_fox_S21_L008_trim_outputs/29_4E_fox_S21_L008_R2_unpaired.fastq.gz \
LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:35