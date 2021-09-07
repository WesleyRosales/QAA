#!/bin/bash 

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --nodes=1
#SBATCH --job-name=QAA_hist
#SBATCH --error=QAA_hist.err
#SBATCH --output=QAA_hist.out
#SBATCH --time=0-12:00:00
#SBATCH --cpus-per-task=8

cd /projects/bgmp/wrosales/bioinfo/Bi623/QAA/
conda activate QAA
/usr/bin/time -v ./read_length_dist.py -pr1 11_2H_both_S9_L008_trim_outputs/11_2H_both_S9_L008_R1_paired.fastq.gz -ur1 11_2H_both_S9_L008_trim_outputs/11_2H_both_S9_L008_R1_unpaired.fastq.gz -pr2 11_2H_both_S9_L008_trim_outputs/11_2H_both_S9_L008_R2_paired.fastq.gz -ur2 11_2H_both_S9_L008_trim_outputs/11_2H_both_S9_L008_R2_unpaired.fastq.gz -o 11_2H_both_S9_L008_trim_outputs/11_2H_both_S9_L008_read_length_dist.png
/usr/bin/time -v ./read_length_dist.py -pr1 29_4E_fox_S21_L008_trim_outputs/29_4E_fox_S21_L008_R1_paired.fastq.gz -ur1 29_4E_fox_S21_L008_trim_outputs/29_4E_fox_S21_L008_R1_unpaired.fastq.gz -pr2 29_4E_fox_S21_L008_trim_outputs/29_4E_fox_S21_L008_R2_paired.fastq.gz -ur2 29_4E_fox_S21_L008_trim_outputs/29_4E_fox_S21_L008_R2_unpaired.fastq.gz -o 29_4E_fox_S21_L008_trim_outputs/29_4E_fox_S21_L008_read_length_dist.png
