#!/bin/bash

#SBATCH -c 4
#SBATCH -t 0-01:00
#SBATCH -p short
#SBATCH --mem=1G
#SBATCH -e logs/Asite_%j.err
#SBATCH -o logs/Asite_%j.log
#SBATCH --mail-type=END
#SBATCH --mail-user=mtcouvi@gmail.com

### WRITTEN:		05/1/2020 by Mary Couvillion
###					updated 11/2022 for yeast

### USE: 			This script will transform .bam yeast mito ribosome profiling output 
###					to A site bedGraphs for viewing on IGV. 
###					And will count frame preferences across mito mRNAs
###					LibName sizeRange

### REQUIREMENTS:	logs directory (in same directory)
###					.bam file (in same directory)
###					
###					SAM2hMitoFBED_softClip.py
###					chrM.chrom.sizes (for making A site bedGraphs)
###					FrameCountBed_chrM_ignore1st3.py
###					cutCodons_hMito_noOverlap.bed (used for frame count)

# load required modules
module load gcc/6.2.0
module load python/2.7.12
module load samtools/1.9
module load bedtools/2.27.1



# Input from command line. LibName is a short descriptive name you want as prefix to all files
LibName=$1
Exp=$2 # 
sizeRange=$3 # e.g. 37to41 

Ascript="./Scripts/SAM2yMitoFPBEDplus_MTC.py"


### A site bedGraph
echo "Making A site bed"
# bam to sam
samtools view -h -o ${LibName}_Aligned.Mito_mRNA.noDups_for${sizeRange}.sam  ${LibName}_Aligned.Mito_mRNA.noDups.bam 
# ## Custom script to make bed file
python $Ascript -i ${LibName}_Aligned.Mito_mRNA.noDups_for${sizeRange}.sam -s ${sizeRange}
# 
# # # convert to bedGraph
echo "Converting to bedGraph"
# # # plus
genomeCoverageBed -bga -trackline -i ${LibName}_Aligned.Mito_mRNA.noDups.Asite_${sizeRange}_P.bed -g ./SeqFiles/bowtie_index/sacCer3.chrom.sizes > ${LibName}_Mito_mRNA.noDups.Asite_${sizeRange}_P.bedGraph
# 
rm ${LibName}_Aligned.Mito_mRNA.noDups_for${sizeRange}.sam 
rm ${LibName}_Aligned.Mito_mRNA.noDups.Asite_${sizeRange}_P.bed


############ Fill all positions in bedGraph files ###########
echo "Filling positions in bedGraph"
python ./Scripts/FillMissingPositionsBedGraph.py -p ${LibName}_Mito_mRNA.noDups.Asite_${sizeRange}_P.bedGraph -s ${sizeRange} 

# # ######## Frame Count ###########
# 
echo "Counting frames"
python ./Scripts/FrameCountBed_chrM_yMito.py ${LibName}_Mito_mRNA.noDups.Asite_${sizeRange}_P.bedGraph ${LibName}_Asite_${sizeRange}


