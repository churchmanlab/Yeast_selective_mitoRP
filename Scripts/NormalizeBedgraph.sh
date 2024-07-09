#!/bin/bash

#SBATCH -c 1
#SBATCH -t 0-0:10
#SBATCH -p short
#SBATCH --mem=1G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mail@gmail.com



LibName=$1
suffix=$2 # Mito_mRNA.noDups.5p.plus Mito_mRNA.noDups.Asite_${sizeRange}_P
sizeRange=$3



# Plus
infile=${LibName}_${suffix}.bedGraph
outfile=${LibName}_${suffix}.rpm.bedGraph


# Get sum of all reads
sum=`awk 'BEGIN {OFS = "\t"} {sum+= $4*($3-$2); next;} END {print sum}' ${infile}`
# Normalize
awk 'BEGIN {OFS = "\t"} {print $1, $2, $3, $4/s*1000000}' s=$sum $infile > $outfile

############ Fill all positions in bedGraph files ###########
echo "Filling positions in bedGraph"
python ./Scripts/FillMissingPositionsBedGraph.py -p ${LibName}_${suffix}.rpm.bedGraph
python ./Scripts/FillMissingPositionsBedGraph.py -p ${LibName}_${suffix}.bedGraph

