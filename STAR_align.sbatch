#!/bin/bash 

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --nodes=1
#SBATCH --job-name=STAR_align
#SBATCH --error=STAR_align.err
#SBATCH --output=STAR_align.out
#SBATCH --time=0-10:00:00
#SBATCH --cpus-per-task=8

cd /projects/bgmp/wrosales/bioinfo/Bi623/QAA
conda activate QAA
/usr/bin/time -v STAR --runThreadN 8 --runMode alignReads \
--outFilterMultimapNmax 3 \
--outSAMunmapped Within KeepPairs \
--alignIntronMax 1000000 --alignMatesGapMax 1000000 \
--readFilesCommand zcat \
--readFilesIn /projects/bgmp/wrosales/bioinfo/Bi623/QAA/11_2H_both_S9_L008_trim_outputs/11_2H_both_S9_L008_R1_paired.fastq.gz /projects/bgmp/wrosales/bioinfo/Bi623/QAA/11_2H_both_S9_L008_trim_outputs/11_2H_both_S9_L008_R2_paired.fastq.gz \
--genomeDir /projects/bgmp/wrosales/bioinfo/Bi623/QAA/Mus_musculus.GRCm39.dna.ens104.STAR_2.7.9a \
--outFileNamePrefix 11_2H_both_S9_L008_alignment
/usr/bin/time -v STAR --runThreadN 8 --runMode alignReads \
--outFilterMultimapNmax 3 \
--outSAMunmapped Within KeepPairs \
--alignIntronMax 1000000 --alignMatesGapMax 1000000 \
--readFilesCommand zcat \
--readFilesIn /projects/bgmp/wrosales/bioinfo/Bi623/QAA/29_4E_fox_S21_L008_trim_outputs/29_4E_fox_S21_L008_R1_paired.fastq.gz /projects/bgmp/wrosales/bioinfo/Bi623/QAA/29_4E_fox_S21_L008_trim_outputs/29_4E_fox_S21_L008_R2_paired.fastq.gz \
--genomeDir /projects/bgmp/wrosales/bioinfo/Bi623/QAA/Mus_musculus.GRCm39.dna.ens104.STAR_2.7.9a \
--outFileNamePrefix 29_4E_fox_S21_L008_alignment