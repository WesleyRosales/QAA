#!/bin/bash 

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --nodes=1
#SBATCH --job-name=htseqcount
#SBATCH --error=htseqcount.err
#SBATCH --output=htseqcount.out
#SBATCH --time=0-10:00:00
#SBATCH --cpus-per-task=8

cd /projects/bgmp/wrosales/bioinfo/Bi623/QAA
conda activate QAA
htseq-count --stranded=yes 11_2H_both_S9_L008_alignmentAligned.out.sam Mus_musculus.GRCm39.104.gtf > 11_2H_both_S9_L008_stranded_htseq_out.txt
htseq-count --stranded=no 11_2H_both_S9_L008_alignmentAligned.out.sam Mus_musculus.GRCm39.104.gtf > 11_2H_both_S9_L008_nonstranded_htseq_out.txt
htseq-count --stranded=yes 29_4E_fox_S21_L008_alignmentAligned.out.sam Mus_musculus.GRCm39.104.gtf > 29_4E_fox_S21_L008_stranded_htseq_out.txt
htseq-count --stranded=no 29_4E_fox_S21_L008_alignmentAligned.out.sam Mus_musculus.GRCm39.104.gtf > 29_4E_fox_S21_L008_nonstranded_htseq_out.txt
