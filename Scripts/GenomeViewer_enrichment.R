# To install bioconductor,
# source("http://bioconductor.org/biocLite.R")
# biocLite("BiocUpgrade")
# biocLite(c("GenomicFeatures","biomaRt","rtracklayer","Gviz"))

# To choose CRAN mirror in terminal and not rely on the xquartz window popping up:
# chooseCRANmirror(81)



library(biomaRt)
library(GenomicFeatures)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(Gviz)   # plots multiple genomic tracks
library(scales)

####
# Load genomic annotations
####

windowtype = 'rolling'
maxzeros = 3
toPlot='COBfactors' # all
log = 'yes' # no yes
autofindlimits = 'no'
type = 'enrichment' #enrichment rpm
plottype = 'png'

if (toPlot == 'all') {
samples <- c('Pet309', 'Aep3', 'Atp22', 'Aep1','Aep2','Smt1','Cbs1','Cbs2','Cbp1','Pet111')
ht=6
} else if (toPlot == 'some') {
samples <- c('Mrps17','Pet309', 'Aep3', 'Atp22','Cbp1','Aep1','Aep2', 'Pet111')
ht=4.5
ymins=c(0, 0,0,0,0,0,0,0)
ymaxs=c(16000,400,400,400,10,300, 300, 200)
} else if (toPlot == 'COBfactors') {
samples <- c('Mrx4', 'Mrx4_cbp3', 'Cbs1','Cbs1_mrx4','Mrps17_mrx4')
ht=4.5
if (log == 'yes') {
	ymins=c(-2,-2,-2,-2,-2)
	ymaxs=c(6,6,6,6,6)
}
if (log == 'no') {
	ymins=c(0,0,0,0,0)
	ymaxs=c(50,50,80,80,10)
}
} else if (toPlot == 'COBTAs') {
samples <- c('Mrps17', 'Cbp1', 'Cbs1', 'Cbs2')
ht=4.5
if (log == 'yes') {
	ymins=c(0,-3,-3,-3)
	ymaxs=c(16000,4,4,4)
}
if (log == 'no') {
	ymins=c(0,0,0,0)
	ymaxs=c(16000,15,15,15)
}
} else {

samples <- toPlot

ht=2
}
wdth=12


# Exon from ENSEMBL
mart = useMart("ensembl")
ensembl = useDataset("scerevisiae_gene_ensembl",mart)
yeast_exon <- getBM(c("strand","chromosome_name","exon_chrom_start","exon_chrom_end","ensembl_transcript_id","ensembl_gene_id","ensembl_exon_id"),mart=ensembl,filters="chromosome_name",values=c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI","XVII","XVIII","Mito"))
    #"ensembl_exon_id","external_gene_id"

# RangedData
yeast_exon <- RangedData(IRanges(start=yeast_exon$exon_chrom_start,
                                 end=yeast_exon$exon_chrom_end),
                         space=yeast_exon$chromosome_name,
                         strand=yeast_exon$strand,
                         transcript=yeast_exon$ensembl_transcript_id,
                         gene=yeast_exon$ensembl_gene_id,
                         exon=yeast_exon$ensembl_exon_id)  # THIS IS NECESSARY for labeling on TrackPlots
#                        universe="sacCer3")
#,symbol=yeast_exon$external_gene_id



yeast_exon <- as(yeast_exon,"GRanges")
seqlevels(yeast_exon) <- paste("chr",seqlevels(yeast_exon),sep="")
seqlevels(yeast_exon) <- sub('chrMito','chrM', seqlevels(yeast_exon))


## Transmembrane domains from file I made from SGD and ENSEMBL
## Go to SGD genome browser, select S.c. proteome, then do a sequence dump of the track
#yeast_transdomain <- read.delim('~/Desktop/Data/Yeast_genome/TransmembranePositions_nt.txt', header=T)
#yeast_transdomain <- RangedData(IRanges(start=yeast_transdomain$Transmembrane_Domain_Start, end=yeast_transdomain$Transmembrane_Domain_End), space='chrM', strand='+', gene=yeast_transdomain$Ensembl_Gene_ID, symbol='')
#yeast_transdomain <- as(yeast_transdomain, 'GRanges')
#
##Signal Peptide sequences
#yeast_SP <- read.delim('~/Desktop/Data/Yeast_genome/ChrM_SP_nt.txt', header=T)
#yeast_SP <- RangedData(IRanges(start=yeast_SP$Signal_Peptide_start, end=yeast_SP$Signal_Peptide_end), space='chrM', strand='+', gene=yeast_SP$Ensembl_Gene_ID, symbol='')
#yeast_SP <- as(yeast_SP, 'GRanges')

# Sequences from BSgenome
yeast_seqs <- getSeq(Scerevisiae)
names(yeast_seqs)<- c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI","chrM")


#####
# Load experimental data
#####
path <- 'path_to_bedGraphs/bedGraphs/'


for (i in c(1:length(samples))) {

	# Fixed window
	if (windowtype == 'fixed') {
	tab <-import.bedGraph(paste0(path,'enrichment/',samples[i],"_Mito_mRNA.noDups.Asite_37to41_PAll.rpm_Mrps17_normWin10.bedGraph"),genome="sacCer3")
	}
	# Rolling window
	if (windowtype == 'rolling') {
		if (type == 'enrichment') {
			if (samples[i] == 'Mrx4_cbp3') {
				tab <-import.bedGraph(paste0(path,'enrichment/',samples[i],"_Mito_mRNA.noDups.Asite_37to41_PAll.rpm_Mrps17_Cbp3_normRollWin9_max",maxzeros,"zeros.bedGraph"),genome="sacCer3")
				} else if (samples[i] == 'Cbs1_mrx4') {
				tab <-import.bedGraph(paste0(path,'enrichment/',samples[i],"_Mito_mRNA.noDups.Asite_37to41_PAll.rpm_Mrps17_mrx4_normRollWin9_max",maxzeros,"zeros.bedGraph"),genome="sacCer3")
				} else if (samples[i] == 'Mrps17') {
				tab <-import.bedGraph(paste0(path,samples[i],"_Mito_mRNA.noDups.Asite_37to41_PAll.rpm.bedGraph"),genome="sacCer3")
				} else {
				tab <-import.bedGraph(paste0(path,'enrichment/',samples[i],"_Mito_mRNA.noDups.Asite_37to41_PAll.rpm_Mrps17_normRollWin9_max",maxzeros,"zeros.bedGraph"),genome="sacCer3")		
				}
		}
		
		if (type == 'rpm') {
			tab <- import.bedGraph(paste0(path,samples[i],"_Mito_mRNA.noDups.Asite_37to41_PAll.rpm.bedGraph"),genome="sacCer3")	
		}

			
	
		
	}


	if (type == 'enrichment' & samples[i] != 'Mrps17') {
		tab$NA. <- NULL
		tab$NA.1 <- NULL
		if (log == 'yes') {
		tab$NA.2 <- NULL
		tab$Forward <- tab$score 
		tab$Reverse <- tab$score
		}
		if (log == 'no') {
		tab$score <- NULL
		tab$NA.2 <- as.numeric(tab$NA.2)
		tab$Forward <- tab$NA.2
		tab$Reverse <- tab$NA.2
		}
		tab$Forward[tab$Forward < 0] <- 0
		tab$Reverse[tab$Reverse > 0] <- 0
		if (log == 'yes') {
		tab$score <- NULL
		}
		if (log == 'no') {
		tab$NA.2 <- NULL
		}
		} else {
		tab$Forward <- tab$score
		tab$Reverse <- 0
		tab$score <- NULL
		}
	
		assign(samples[i], tab)

}


# allscore=unlist(list(Pet309$score, Aep3$score, Atp22$score, Cbs1$score, Cbs2$score, Cbp1$score, Aep1$score, Aep2$score, Pet111$score))

#####
# Visualize using genome browser
#####

    
# gen <- genome(yeast_tx)
chr <- "chrM"

# Set up axis at top
trackgenome <- GenomeAxisTrack(lwd=1, reverseStrand=FALSE, showTitle=FALSE, cex=0.6, labelPos='below') #cex=0.3
#trackgene <- AnnotationTrack(yeast_tx,group=values(yeast_tx)$tx_name,chromosome=chr,name="Gene")
#track_SP <- GeneRegionTrack(range=yeast_SP, chromosome=chr, shape = 'box', fill='red', name="SP", stackHeight = 0.5)
#track_transdomain <- GeneRegionTrack(range=yeast_transdomain, chromosome=chr, shape = 'box', fill='gray35', name="Trans", stackHeight = 0.5)

# Annotations
trackexon <- GeneRegionTrack(yeast_exon,chromosome=chr,name="Exon", stackHeight=.7, size=0.3, showId=FALSE, col = 'gray30', fill='gray65', lwd=0.5)
#strack <- SequenceTrack(yeast_seqs, chromosome=chr)
    #fontcolor=c(A='green',T='red',C='blue',G='black', N='gray')
#trackcds <- AnnotationTrack(range=yeastcs,chromosome=chr,name="CDS")
#trackmrnaplus <- DataTrack(range=rna_plus,chromosome=chr,name="mRNA +")

# Data

for (i in c(1:length(samples))) {

sample = samples[i]

if (sample == 'Mrps17') {
color1 = 'darkgoldenrod2'
if (log == 'yes') {
color2 = 'gold1' 
} else {color2 = 'darkgoldenrod2'}
}
if (sample == 'Pet309' | sample == 'Pet111') {
color1 = 'indianred3'
if (log == 'yes') {
color2 = 'pink' 
} else {color2 = 'indianred3'}
}
if (sample == 'Aep3' | sample == 'Atp22' | sample == 'Smt1' | sample == 'Aep1' | sample == 'Aep2') {
color1 = 'forestgreen'
if (log == 'yes') {
color2 = 'darkseagreen1'
} else {color2 = 'forestgreen'}
}
if (sample == 'Cbs1' | sample == 'Cbs1'  | sample == 'Cbp1' | sample == 'Cbs1_mrx4') {
color1 = 'dodgerblue'
if (log == 'yes') {
color2 = 'lightblue'
} else {color2 = 'dodgerblue'}
}
if (sample == 'Mrx1') {
color1 = 'darkgoldenrod'
if (log == 'yes') {
color2 = 'wheat1'
} else {color2 = 'darkgoldenrod'}
}
if (sample == 'Smt1') {
color1 = 'seagreen3'
if (log == 'yes') {
color2 = 'seagreen3'
} else {color2 = 'seagreen3'}
}
if (sample == 'Oxa1') {
color1 = 'orange'
if (log == 'yes') {
color2 = 'tan'
} else {color2 = 'orange'}
}
if (sample == 'Mrx4' | sample == 'Mrx4_cbp3') {
color1 = 'purple'
if (log == 'yes') {
color2 = 'violet'
} else {color2 = 'purple'}
}
if (sample == 'Mrps17_mrx4') {
color1 = 'sandybrown'
if (log == 'yes') {
color2 = 'gold1'
} else {color2 = 'tan'}
}

tbl = get(sample)

# Auto-find y limits
if (autofindlimits == 'yes') {
ylimits = c(min(tbl$Reverse, na.rm=TRUE), max(tbl$Forward, na.rm=TRUE))}
# Specify y limits
if (autofindlimits == 'no') {
ylimits = c(ymins[i], ymaxs[i])}


assign(paste0('track',sample), DataTrack(range=tbl, chromosome=chr, name=sample, ylim=ylimits, col=c(color1, color2), groups=c("Forward","Reverse"), baseline=0, col.baseline='grey60', lwd.baseline=.2, cex.sampleNames=1.5))

print(ylimits)

}


colors=c('grey30', 'grey70')








# Plot




# Full genome
if (type == 'enrichment' & plottype == 'pdf') {
pdf(paste0('path_to_output/GenomeBrowser/SelRP_WholeGenome_',toPlot,'_max',maxzeros,'zeros_',log,'log.pdf'),
     width=wdth, # width=4.5, width=2, width=2.5
     height=ht, #4
     pointsize=12)
}
if (type == 'rpm' & plottype == 'pdf') {
pdf(paste0('path_to_output/GenomeBrowser/SelRP_WholeGenome_',toPlot,'rpm.pdf'),
     width=wdth, # width=4.5, width=2, width=2.5
     height=ht, #4
     pointsize=12)
}
if (type == 'enrichment' & plottype == 'png') {
png(paste0('path_to_output/GenomeBrowser/SelRP_WholeGenome_',toPlot,'_max',maxzeros,'zeros_',log,'log.png'),
     width=wdth, # width=4.5, width=2, width=2.5
     height=ht, #4
     units='in',
     res=350)
}
if (type == 'rpm' & plottype == 'png') {
png(paste0('path_to_output/GenomeBrowser/SelRP_WholeGenome_',toPlot,'rpm.png'),
     width=wdth, # width=4.5, width=2, width=2.5
     height=ht, #4
     units='in',
     res=350)
}


n=length(samples)
datatracks = c(trackgenome)
for (i in c(1:length(samples))) {
datatrack <- get(paste0('track', samples[i]))
datatracks = c(datatracks, datatrack)
}
datatracks = c(datatracks, trackexon)


plotTracks(datatracks,from=12000, 83000, sizes=c(.25,rep(.3, n),.3), type = 'histogram', background.title='gray50', cex.title=0.6, cex.axis=0.3, legend=FALSE) 

dev.off()



# Zoom
gene='COX1_ATP8_ATP6'
if (gene == 'COB') {
st=36000
end=44000
extrackh=.1
}
if (gene == 'COX1_ATP8_ATP6') {
st=12400
end=29500
extrackh=.3
}

if (type == 'enrichment') {
pdf(paste0('path_to_output/GenomeBrowser/SelRP_WholeGenome_',toPlot,'_max',maxzeros,'zeros_',log,'log_zoom',gene,'.pdf'),
     width=wdth, # width=4.5, width=2, width=2.5
     height=ht, #4
     pointsize=12)
}
if (type == 'rpm') {
pdf(paste0('path_to_output/GenomeBrowser/SelRP_WholeGenome_',toPlot,'rpm_zoom.pdf'),
     width=wdth, # width=4.5, width=2, width=2.5
     height=ht, #4
     pointsize=12)
}

plotTracks(datatracks,from=st, end, sizes=c(.25,rep(.3, n),.3), type = 'histogram', background.title='gray50', cex.title=0.6, cex.axis=0.3, legend=FALSE) 

dev.off()




# source('/Users/Mary/Desktop/Data/SelectiveRP/Scripts/GenomeViewer_enrichment.R')

