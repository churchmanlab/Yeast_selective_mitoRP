# Commands for submitting jobs to slurm job scheduler: copy-paste into command-line



####################################################################
1. Trim raw reads, align, get library makeups, footprint sizes, and make 5' and 3' bedGraphs
####################################################################

Libs=("Mrps17_1" "Mrps17_2" "Aep2_1" "Aep2_2")

UMI="3p10_5p4"
COX2oligo="no" # whether to remove sequence from a synthetic oligo sometimes added to library prep
types=("Mono" "Mono" "Mono" "Mono") # other options are 40S or Di80 for aligning longer reads

mkdir logs

script="ProcessFASTQ_yMitoRP.sh"
for i in {0..3}
do
sbatch -e logs/Map_${Libs[i]}.err -o logs/Map_${Libs[i]}.log ./Scripts/${script} ${Libs[i]} ${Libs[i]}.fastq.gz $UMI $COX2oligo ${types[i]}
done



# Combine footprint length files
Experiment="SelRP"
Libs=("Mrps17_1" "Mrps17_2" "Aep2_1" "Aep2_2") 

echo 'Length' > ${Experiment}_lengths.txt
for x in {11..55} # length range of choice
do
echo $x >> ${Experiment}_lengths.txt
done 

list='Aligned.Mito_mRNA.noDups'

for type in $list
do
for i in {0..3}
do
join -j 1 -a1 -e 0 -o 1.1,2.2 ${Experiment}_lengths.txt ${Libs[i]}_${type}_samstatsReadlength.noSoft.txt > ${Libs[i]}_${type}_samstatsReadlengthAll.noSoft.txt
done
done

for type in $list
do
paste <(awk '{print $1"\t"$2}' ${Libs[0]}_${type}_samstatsReadlengthAll.noSoft.txt) <(awk '{print $2}' ${Libs[1]}_${type}_samstatsReadlengthAll.noSoft.txt) <(awk '{print $2}' ${Libs[2]}_${type}_samstatsReadlengthAll.noSoft.txt) > ${Experiment}_${type}_noSoftLengthDist.txt
done

# Use the resulting file to plot read length distributions
# Paths, colors, names, etc need to be updated
ReadLengthPlot.R
# In R:
source('path_to_script/ReadLengthPlot.R')



####################################################################
2. Get periodicity for each length of read
####################################################################

Libs="Mrps17_1 Mrps17_2 Aep2_1 Aep2_2
min=35 # This is the range of sizes periodicity will be counted for
max=44
script="CountFramePerLength.sh" 
for lib in $Libs
do
sbatch -e logs/FramePerLength_${lib}.err -o logs/FramePerLength_${lib}.log ./Scripts/${script} $lib $min $max
done

# Combine
Experiment="SelRP"
Libs=("Mrps17_1" "Mrps17_2" "Aep2_1" "Aep2_2")

ends='5 3'
for end in $ends
do
suffix="_${end}p_FrameCounts.txt"
outsuffix="_${end}p_allFrameCounts.txt"
paste <(awk '{print $1"\t"$2}' ${Libs[0]}_${end}p_FrameCounts.txt) <(awk '{print $2}' ${Libs[1]}_${end}p_FrameCounts.txt) <(awk '{print $2}' ${Libs[2]}_${end}p_FrameCounts.txt) <(awk '{print $2}' ${Libs[3]}_${end}p_FrameCounts.txt) > ${Experiment}${outsuffix}
done

# Use the resulting files to plot the periodicity for each length in order to help decide what size ranges to use for A site transformation (copy-paste into excel or other)



####################################################################
3. Make 5' and 3' length beds for V plots
####################################################################


Experiment="SelRP"
Libs="Mrps17_1 Mrps17_2 Aep2_1 Aep2_2"

script="MakeLengthBeds.sh" # Update paths within
for lib in $Libs
do
sbatch -e logs/LengthBeds_${lib}.err -o logs/LengthBeds_${lib}.log ./Scripts/${script} $lib
done

# Use resulting files to make V plots
ScatterPlotsLengths.R 
# In R:
source('path_to_script/ScatterPlotsLengths.R')



####################################################################
4. A site transform and count frame
####################################################################
Input:
_Aligned.Mito_mRNA.noDups.bam 
Outputs:
_Mito_mRNA.noDups.Asite_${sizeRange}_P.bedGraph
_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.bedGraph # This lists out every position on chrM rather than the condensed bedGraph format

Exp="SelRP"
Libs="Mrps17_1 Mrps17_2 Aep2_1 Aep2_2"

sizerange='37to41' 
sizes="37 38 39 40 41"
offsets="15 16 17 17 17"

script="AsiteAndCountFrame.sh"
for lib in $Libs
do
sbatch /n/groups/churchman/mc348/yMitoRP/Scripts/$script $lib $Exp $sizerange "$sizes" "$offsets"
done


####################################################################
5. Normalize bedGraphs by rpm
####################################################################
Input: 
${LibName}_${suffix}.bedGraph
Outputs: 
${LibName}_${suffix}.rpm.bedGraph

Libs="Mrps17_1 Mrps17_2 Aep2_1 Aep2_2"
suffixes="Mito_mRNA.noDups.Asite_37to41_P Mito_mRNA.noDups.Asite_37to41_PAll" 

for lib in $Libs
do
for suffix in $suffixes
do
sbatch ./Scripts/NormalizeBedgraph.sh $lib $suffix
done
done





####################################################################
6. Combine replicates from Asite bedgraphs, with all positions filled
####################################################################
Input:
${SampName}_1_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.bedGraph
${SampName}_2_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.bedGraph
Outputs:
${SampName}_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.bedGraph
${SampName}_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.rpm.bedGraph

Samps=("Mrps17" "Aep2")
numReps=("2" "6")

sizeRange="37to41"

script="SumBedgraphs.sh"
for i in {0..1} 
do
sbatch ./Scripts/${script} ${Samps[i]} ${numReps[i]} $sizeRange
done




mkdir rpmNormalizedBedgraphs
mv *.rpm.bedGraph rpmNormalizedBedgraphs/.
cd rpmNormalizedBedgraphs


####################################################################
7. Normalize bedGraphs by fold change compared to control
####################################################################

Libs="Aep2"
normSamp="Mrps17"

window="9" 
libsuffix='_Mito_mRNA.noDups.Asite_37to41_PAll.rpm'
normsuffix='_Mito_mRNA.noDups.Asite_37to41_PAll.rpm'
rawfilesuffix='_Mito_mRNA.noDups.Asite_37to41_PAll'


script="NormBedgraphToControl_v5.R" # _v4 is fixed window, _v5 is rolling window
zerothresh=3 # max number of zeros in the window
for lib in $Libs
do
sbatch -p short -t 0-00:30 --mem=1G --wrap="Rscript /n/groups/churchman/mc348/yMitoRP/Scripts/${script} $normSamp $lib $window $libsuffix $normsuffix $rawfilesuffix $zerothresh" 
done



####################################################################
8. Plotting for figures
####################################################################

# Plotting across the genome using Gviz
GenomeViewer_enrichment.R
# In R:
source('path_to_script/GenomeViewer_enrichment.R')


# Plotting across coding sequences, option for adding protein features like transmembrane domains
TAEnrichment.R
# In R:
source('path_to_script/TAEnrichment.R')


# Quantifying relative enrichment and barplot
CountReadsInRegion_includeMaturases.R
# In R:
source('path_to_script/CountReadsInRegion_includeMaturases.R')











