#!/bin/bash

#SBATCH -c 4
#SBATCH -t 0-10:00
#SBATCH -p short
#SBATCH --mem=80G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mail@gmail.com

### WRITTEN:		11/1/22 by Mary Couvillion based on python script written ~2014

### USE: 			This script will trim reads, map to the genome, remove PCR duplicates, 
###					and produce files for visualizing on IGV
###					sbatch -e logs/LibName.err -o logs/LibName.log ProcessFASTQ_yMitoRP_2022.sh
###					LibName FASTQfile UMItype


### REQUIREMENTS:	logs directory (in same directory)
###					fastq.gz files (in same directory)
###					
###					extractMolecularBarcodeFrom3pr10AND5pr4.py
###					removePCRdupsFromBAM_MC20180913.py


# load required modules
module load gcc/6.2.0
module load python/2.7.12
module load samtools/1.9
module load bedtools/2.27.1
module load cutadapt/1.14
module load R/4.0.1
module load bowtie/1.2.1.1
module load bamtools/2.4.1
module load fastqc/0.11.3
module load fastx/0.0.13
module load igvtools/2.3.88
module load java/jdk-1.8u112
module load star/2.7.3a
module load picard/2.8.0

# Input from command line. LibName is a short descriptive name you want as prefix to all files
LibName=$1
InputFASTQ=$2
UMI=$3 # 3p6, 3p6_5p4, 3p10_5p4 Rooijers Pearce Li
COX2oligo=$4
type=$5
PWD=`pwd` # This is to get the path to the directory you're working in for the featureCounts_MitoRP.R script so it finds the right files to read and write


# ####### CLEAN READS ###########

if [ "${UMI}" = "3p10_5p4" ]
then
	# Cutadapt (minimum is 31-10-4=17 after trimming --untrimmed-output writes reads with     		  no adaptor found INSTEAD of writing to final output file. Do this here so that we keep only reads that have adapter, which means we know where the barcode is
	cutadapt -e 0.2 --untrimmed-output ${LibName}_noAdaptor.fastq -a CTGTAGGCACCATCAAT -m 31 ${InputFASTQ} -o ${LibName}_trimmed.fastq
	rm ${LibName}_noAdaptor.fastq
	###### To extract barcode ######
	python ./Scripts/extractMolecularBarcodeFrom3pr10AND5pr4.py ${LibName}_trimmed.fastq ${LibName}_Cleaned.fastq
	rm ${LibName}_trimmed.fastq
fi
# NOTE MTC 11/1/22: Got rid of the 5' 1 nt trim, since we're aligning with STAR now. Instead allow a 5' 1 nt mismatch in bowtie steps
if [ "${UMI}" = "none" ]
then
	cutadapt -e 0.2 -a CTGTAGGCACCATCAAT -m 17 ${InputFASTQ} -o ${LibName}_trimmed.fastq
	mv ${LibName}_trimmed.fastq ${LibName}_Cleaned.fastq
fi


####### MAP READS ###########


# Remove COX2 oligo sequence
if [ "${COX2oligo}" = "yes" ]
then
./Scripts/programs/bbmap/bbduk.sh in=${LibName}_Cleaned.fastq hdist=1 k=30 literal=TTAACAACATTCATTATGAATGATGTACCAACACCTTATGC out=${LibName}_COX2filtered.fastq outm=${LibName}_COX2oligo.fastq # Add an extra T here because it's added in RT then the longer version is not soft clipped and is not removed if I don't
fi
if [ "${COX2oligo}" = "no" ]
then
mv ${LibName}_Cleaned.fastq ${LibName}_COX2filtered.fastq
fi

# Bowtie
# To get counts (number of mapped reads) for each step, look in error file (.err)

# Nuclear rRNA
bowtie -v 1 -5 1 -S ./SeqFiles/bowtie_index/rna_coding_Nuclear_rRNA ${LibName}_COX2filtered.fastq ${LibName}_Nuc_rRNA_almts.sam --un ${LibName}_Nuc_rRNA_un.fastq
rm ${LibName}_Nuc_rRNA_almts.sam

#Nuclear tRNA, 3' 3nt mismatch to account for untemplated -CCA
bowtie -v 1 -3 3 -5 1 -S ./SeqFiles/bowtie_index/rna_coding_Nuclear_tRNA ${LibName}_Nuc_rRNA_un.fastq ${LibName}_Nuc_tRNA_almts.sam --un ${LibName}_Nuc_tRNA_un.fastq
rm ${LibName}_Nuc_tRNA_almts.sam

#Nuclear other ncRNA
bowtie -v 1 -5 1 -S ./SeqFiles/bowtie_index/rna_coding_Nuclear_ncRNAother ${LibName}_Nuc_tRNA_un.fastq ${LibName}_Nuc_ncRNA_almts.sam --un ${LibName}_Nuc_ncRNA_filtered.fastq
rm ${LibName}_Nuc_ncRNA_almts.sam

#Mito rRNA
bowtie -v 1 -5 1 -S ./SeqFiles/bowtie_index/rna_coding_Mito_rRNA_edited ${LibName}_Nuc_ncRNA_filtered.fastq ${LibName}_yMito_rRNA_almts.sam --al ${LibName}_yMito_rRNA_al.fastq --un ${LibName}_yMito_rRNA_un.fastq
rm ${LibName}_yMito_rRNA_almts.sam

#Mito tRNA, 3' 3nt mismatch to account for untemplated -CCA
bowtie -v 1 -3 3 -5 1 -S ./SeqFiles/bowtie_index/rna_coding_Mito_tRNA_edited ${LibName}_yMito_rRNA_un.fastq ${LibName}_yMito_tRNA_almts.sam --al ${LibName}_yMito_tRNA_al.fastq --un ${LibName}_yMito_tRNA_un.fastq
rm ${LibName}_yMito_tRNA_almts.sam

#Mito other ncRNA
bowtie -v 1 -5 1 -S ./SeqFiles/bowtie_index/rna_coding_chrM_ncRNAother ${LibName}_yMito_tRNA_un.fastq ${LibName}_Mito_ncRNA_almts.sam --al ${LibName}_Mito_ncRNA_al.fastq --un ${LibName}_all_ncRNA_filtered.fastq


## Use STAR to map, intron-aware
## This STAR index is different from the one used to map cytoRP reads: the --sjdbOverhang is bigger (45)

## To make index (do this ahead of time): 
# mkdir S288C_refseq_chrAll_45bp
# sbatch -p short -t 0-3:00 -c 4 --mem=10G --wrap="STAR --runMode genomeGenerate --runThreadN 4 --genomeDir S288C_refseq_chrAll_80bp --genomeFastaFiles ./SeqFiles/bowtie_index/S288C_refseq_chrAll.fa --sjdbGTFfile ./SeqFiles/bowtie_index/S288C_refseq_chrAll.gtf --sjdbOverhang 80 --genomeSAindexNbases 10"

# Align nuc ncRNA-filtered reads
if [ "${type}" = "Di80" ] || [ "${type}" = "40S" ]
then
genome="path_to_star_index/STAR_index/S288C_refseq_chrAll_80bp"
else
genome="path_to_star_index/STAR_index/S288C_refseq_chrAll_45bp"
fi

STAR --runThreadN 4 --genomeDir $genome --readFilesIn ${LibName}_Nuc_ncRNA_filtered.fastq --outFileNamePrefix ${LibName}_ --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM NM MD --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFilterMultimapNmax 1000 --outFilterMismatchNoverReadLmax 0.07 --outFilterMismatchNmax 3 --outReadsUnmapped Fastx --limitBAMsortRAM 3000000000 --alignSJDBoverhangMin 1 --alignIntronMax 5000

rm -r ${LibName}__STARtmp
rm ${LibName}_SJ.out.tab


# Map to oGAB
bowtie -v 1 -5 1 -S ./SeqFiles/bowtie_index/oGAB ${LibName}_Unmapped.out.mate1 ${LibName}_oGAB_almts.sam
rm ${LibName}_oGAB_almts.sam

# # ############ Get counts from mapping ###########
# # ## Make list of mapping references, in order
echo Input > ${LibName}_MapList.txt
echo Nuc rRNA >> ${LibName}_MapList.txt
echo Nuc tRNA >> ${LibName}_MapList.txt
echo Nuc other ncRNA >> ${LibName}_MapList.txt
echo Mito rRNA >> ${LibName}_MapList.txt
echo Mito tRNA >> ${LibName}_MapList.txt
echo Mito other ncRNA >> ${LibName}_MapList.txt
echo oGAB >> ${LibName}_MapList.txt
# Grab mapped numbers from bowtie error file
grep -m 1 'reads processed:' logs/Map_${LibName}.err > ${LibName}_bowtieCounts.txt
grep 'reads with at least one reported alignment:' logs/Map_${LibName}.err >> ${LibName}_bowtieCounts.txt
paste <(awk '{print $0}' ${LibName}_MapList.txt) <(awk -F ': ' '{print $2}' ${LibName}_bowtieCounts.txt) > ${LibName}_Counts.txt
# Get mapped reads from STAR log file
echo Mapped from STAR >> ${LibName}_Counts.txt
egrep 'Uniquely mapped reads number | Number of reads mapped to multiple loci' ${LibName}_Log.final.out | awk -F '|\t' '{print $2}' | awk '{sum += $1} END {print sum}' >> ${LibName}_Counts.txt

# Get unmapped reads from STAR log file
echo Unmapped from STAR \(oGAB counts should be subtracted from this\): >> ${LibName}_Counts.txt
egrep 'Number of reads mapped to too many loci | Number of reads unmapped: too short | Number of reads unmapped: other' ${LibName}_Log.final.out | awk -F '|\t' '{print $2}' | awk '{sum += $1} END {print sum}' >> ${LibName}_Counts.txt

# Get COX2 oligo counts
echo COX2 oligo >> ${LibName}_Counts.txt
wc -l ${LibName}_COX2oligo.fastq >> ${LibName}_Counts.txt

rm ${LibName}_MapList.txt
rm ${LibName}_bowtieCounts.txt





########## REMOVE PCR DUPLICATES ###########

## Make BAM index (samtools)
samtools index -@ 3 ${LibName}_Aligned.sortedByCoord.out.bam

if [ "${UMI}" != "none" ]
then

## Removal (updated 9/2018 to include CIGAR string as well)
python ./Scripts/removePCRdupsFromBAM_MC20180913.py ${LibName}_Aligned.sortedByCoord.out.bam ${LibName}_Aligned_noDups.bam

else 

cp ${LibName}_Aligned.sortedByCoord.out.bam ${LibName}_Aligned_noDups.bam

fi

############ MAKE MITO mRNA ONLY BAM FILE ###########

## These should be the footprints
## Make sure they are sorted by coordinate
samtools sort -@ 3 -o ${LibName}_Aligned_noDups.sort.bam ${LibName}_Aligned_noDups.bam
# Index 
samtools index -@ 3 ${LibName}_Aligned_noDups.sort.bam

## Mito only
samtools view -@ 3 -b ${LibName}_Aligned_noDups.sort.bam -o ${LibName}_chrMall_noDups.bam chrM 

## mito mRNAs only by specifying mito bed file
samtools view -@ 3 -b ${LibName}_Aligned_noDups.sort.bam -o ${LibName}_Aligned.Mito_mRNAplus.noDups.bam -L ./SeqFiles/yMito_mRNAsANDsurrounding.bed 

## And remove the rRNA intron-spanning reads that overlap SceI (actually the whole first exon of 21S rRNA
samtools view -@ 3 -h ${LibName}_Aligned.Mito_mRNAplus.noDups.bam | awk '{if ($1 ~ /^@/ || $4 < 58007 || $4 > 60724) print $0}' | samtools view -@ 3 -h -b -o ${LibName}_Aligned.Mito_mRNA.noDups.bam

rm ${LibName}_Aligned.Mito_mRNAplus.noDups.bam

# ############ Get counts from noDups ###########
echo '' >> ${LibName}_Counts.txt
echo STAR mapping locations with duplicates: >> ${LibName}_Counts.txt
samtools view -@ 3 -F 0x4 -c ${LibName}_Aligned.sortedByCoord.out.bam >> ${LibName}_Counts.txt
echo STAR mapping locations noDups: >> ${LibName}_Counts.txt
samtools view -@ 3 -F 0x4 -c ${LibName}_Aligned_noDups.bam >> ${LibName}_Counts.txt
echo Mito mRNA noDups \(this is with local alignment and mismatches allowed\): >> ${LibName}_Counts.txt
samtools view -@ 3 -F 0x4 -c ${LibName}_Aligned.Mito_mRNA.noDups.bam >> ${LibName}_Counts.txt


# Add header to Counts file
sed -i '1i\ \n'${LibName} ${LibName}_Counts.txt



# Get read length distributions after removing soft clipped bases

list='Aligned.Mito_mRNA.noDups'
for type in $list
do

# convert from .bam to .sam
samtools view -@ 3 -h -O SAM -o ${LibName}_${type}.sam ${LibName}_${type}.bam

# remove soft clipped bases from reads
java -jar ./Scripts/programs/jvarkit/biostar84452.jar -o ${LibName}_${type}.noSoft.sam ${LibName}_${type}.sam
rm ${LibName}_${type}.sam

# convert back to .bam
samtools view -@ 3 -b -h -o ${LibName}_${type}.noSoft.bam ${LibName}_${type}.noSoft.sam
rm ${LibName}_${type}.noSoft.sam

# sort by coord
samtools sort -@ 3 -o ${LibName}_${type}.noSoft.sort.bam ${LibName}_${type}.noSoft.bam  
rm ${LibName}_${type}.noSoft.bam

samtools stats ${LibName}_${type}.noSoft.sort.bam > ${LibName}_${type}_samstats.txt
grep ^RL ${LibName}_${type}_samstats.txt | cut -f 2- > ${LibName}_${type}_samstatsReadlength.noSoft.txt
# Add header (variable has to be outside of quotes for sed)
sed -i '1i\Length\t'${LibName} ${LibName}_${type}_samstatsReadlength.noSoft.txt
rm ${LibName}_${type}_samstats.txt
rm ${LibName}_${type}.noSoft.sort.bam
done


# ############ MAKE BEDGRAPHs FOR VIEWING ON IGV ###########

### Make bam index
samtools index ${LibName}_Aligned.Mito_mRNA.noDups.bam
samtools index ${LibName}_chrMall_noDups.bam

## 5' bedGraph (bedtools genomeCoverageBed)
# mito mRNA reads
# plus
genomeCoverageBed -bga -5 -strand + -trackline -ibam ${LibName}_Aligned.Mito_mRNA.noDups.bam > ${LibName}_Mito_mRNA.noDups.5p.plus.bedGraph
# minus
genomeCoverageBed -bga -5 -strand - -trackline -ibam ${LibName}_Aligned.Mito_mRNA.noDups.bam > ${LibName}_Mito_mRNA.noDups.5p.minus.bedGraph

# All mito reads
# plus
genomeCoverageBed -bga -5 -strand + -trackline -ibam ${LibName}_chrMall_noDups.bam > ${LibName}_chrMall_noDups.5p.plus.bedGraph
# minus
genomeCoverageBed -bga -5 -strand - -trackline -ibam ${LibName}_chrMall_noDups.bam > ${LibName}_chrMall_noDups.5p.minus.bedGraph


## 3' bedGraph (bedtools genomeCoverageBed)
# mito mRNA reads
# plus
genomeCoverageBed -bga -3 -strand + -trackline -ibam ${LibName}_Aligned.Mito_mRNA.noDups.bam > ${LibName}_Mito_mRNA.noDups.3p.plus.bedGraph
# minus
genomeCoverageBed -bga -3 -strand - -trackline -ibam ${LibName}_Aligned.Mito_mRNA.noDups.bam > ${LibName}_Mito_mRNA.noDups.3p.minus.bedGraph



