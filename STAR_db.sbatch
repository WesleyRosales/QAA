#!/bin/bash 

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --nodes=1
#SBATCH --job-name=STARdb
#SBATCH --error=STARdb.err
#SBATCH --output=STARdb.out
#SBATCH --time=0-10:00:00
#SBATCH --cpus-per-task=8

cd /projects/bgmp/wrosales/bioinfo/Bi623/QAA
conda activate QAA
/usr/bin/time -v STAR --runThreadN 8 --runMode genomeGenerate \
--genomeDir /projects/bgmp/wrosales/bioinfo/Bi623/QAA/Mus_musculus.GRCm39.dna.ens104.STAR_2.7.9a \
--genomeFastaFiles /projects/bgmp/wrosales/bioinfo/Bi623/QAA/Mus_musculus.GRCm39.dna.primary_assembly.fa \
--sjdbGTFfile /projects/bgmp/wrosales/bioinfo/Bi623/QAA/Mus_musculus.GRCm39.104.gtf