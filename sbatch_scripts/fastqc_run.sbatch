#!/bin/bash 

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --nodes=1
#SBATCH --job-name=fastqc
#SBATCH --error=fastqc.err
#SBATCH --output=fastqc.out
#SBATCH --time=0-12:00:00
#SBATCH --cpus-per-task=8

cd /projects/bgmp/wrosales/bioinfo/Bi623/QAA/
conda activate QAA
/usr/bin/time -v fastqc -o 11_2H_both_S9_L008_fastqc_outputs /projects/bgmp/shared/2017_sequencing/demultiplexed/11_2H_both_S9_L008_R1_001.fastq.gz /projects/bgmp/shared/2017_sequencing/demultiplexed/11_2H_both_S9_L008_R2_001.fastq.gz
/usr/bin/time -v fastqc -o 29_4E_fox_S21_L008_fastqc_outputs /projects/bgmp/shared/2017_sequencing/demultiplexed/29_4E_fox_S21_L008_R1_001.fastq.gz /projects/bgmp/shared/2017_sequencing/demultiplexed/29_4E_fox_S21_L008_R2_001.fastq.gz 
