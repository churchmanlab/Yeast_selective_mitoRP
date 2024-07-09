#!/bin/bash

#SBATCH -c 1
#SBATCH -t 0-0:10
#SBATCH -p short
#SBATCH --mem=1G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mail@gmail.com



SampName=$1 # Pet309
numReps=$2 # 2
sizeRange=$3 # 37to41

outfile1=${SampName}_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.bedGraph
outfile2=${SampName}_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.rpm.bedGraph

if [ "$numReps" = "1" ]
then
infile1=${SampName}_1_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.bedGraph

mv $infile1 $outfile1

elif [ "$numReps" = "2" ]
then
infile1=${SampName}_1_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.bedGraph
infile2=${SampName}_2_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.bedGraph

paste $infile1 $infile2 | awk 'BEGIN {OFS = "\t"} NR > 1 { print $1, $2, $3, $4+$8 }' > $outfile1

elif [ "$numReps" = "3" ]
then
infile1=${SampName}_1_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.bedGraph
infile2=${SampName}_2_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.bedGraph
infile3=${SampName}_3_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.bedGraph

paste $infile1 $infile2 $infile3 | awk 'BEGIN {OFS = "\t"} NR > 1 { print $1, $2, $3, $4+$8+$12}' > $outfile1

elif [ "$numReps" = "4" ]
then
infile1=${SampName}_1_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.bedGraph
infile2=${SampName}_2_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.bedGraph
infile3=${SampName}_3_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.bedGraph
infile4=${SampName}_4_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.bedGraph

paste $infile1 $infile2 $infile3 $infile4 | awk 'BEGIN {OFS = "\t"} NR > 1 { print $1, $2, $3, $4+$8+$12+$16}' > $outfile1

elif [ "$numReps" = "5" ]
then
infile1=${SampName}_1_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.bedGraph
infile2=${SampName}_2_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.bedGraph
infile3=${SampName}_3_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.bedGraph
infile4=${SampName}_4_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.bedGraph
infile5=${SampName}_5_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.bedGraph

paste $infile1 $infile2 $infile3 $infile4 $infile5 | awk 'BEGIN {OFS = "\t"} NR > 1 { print $1, $2, $3, $4+$8+$12+$16+$20}' > $outfile1


paste $infile1 $infile2 $infile3 $infile4 | awk 'BEGIN {OFS = "\t"} NR > 1 { print $1, $2, $3, $4+$8+$12+$16}' > $outfile1

elif [ "$numReps" = "6" ]
then
infile1=${SampName}_1_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.bedGraph
infile2=${SampName}_2_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.bedGraph
infile3=${SampName}_3_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.bedGraph
infile4=${SampName}_4_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.bedGraph
infile5=${SampName}_5_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.bedGraph
infile6=${SampName}_6_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.bedGraph

paste $infile1 $infile2 $infile3 $infile4 $infile5 | awk 'BEGIN {OFS = "\t"} NR > 1 { print $1, $2, $3, $4+$8+$12+$16+$20+$24}' > $outfile1

fi

# Get sum of all reads
sum=`awk 'BEGIN {OFS = "\t"} {sum+= $4*($3-$2); next;} END {print sum}' ${outfile1}`
# Normalize
awk 'BEGIN {OFS = "\t"} {print $1, $2, $3, $4/s*1000000}' s=$sum $outfile1 > $outfile2

sed  -i '1i\track type=bedGraph' $outfile1
sed  -i '1i\track type=bedGraph' $outfile2


