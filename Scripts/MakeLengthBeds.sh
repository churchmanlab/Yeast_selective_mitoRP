#!/bin/bash

#SBATCH -c 4
#SBATCH -t 0-4:00
#SBATCH -p short
#SBATCH --mem=5G
#SBATCH -e logs/LengthBed_%j.err
#SBATCH -o logs/LengthBed_%j.log
#SBATCH --mail-type=END
#SBATCH --mail-user=mail@gmail.com

### WRITTEN:		05/14/2020 by Mary Couvillion

### USE: 			This script will make 5' and 3' bed files with Mito position and 
###					read length
###					sbatch MakeLengthBeds.sh
###					LibName


module load python/3.7.4


LibName=$1

# BAM to SAM
samtools view -@ 3 -h -o ${LibName}_Aligned.Mito_mRNA.noDups_forLengthBed.sam  ${LibName}_Aligned.Mito_mRNA.noDups.bam


# 5' bed
python ./Scripts/SAM2length5pBED_softClip.py -f ${LibName}_Aligned.Mito_mRNA.noDups_forLengthBed.sam  -p ${LibName}_Mito_mRNA_lengths5p_P.bed

# 3' bed
python ./Scripts/SAM2length3pBED_softClip.py -f ${LibName}_Aligned.Mito_mRNA.noDups_forLengthBed.sam  -p ${LibName}_Mito_mRNA_lengths3p_P.bed

rm ${LibName}_Aligned.Mito_mRNA.noDups_forLengthBed.sam  



