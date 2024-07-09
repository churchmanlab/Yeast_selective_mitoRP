# All read sizes used here, comparison between samples only, not set up to calculate gene expression across genes (RPK)
# Updated 6/6/24 from CountReadsInRegion to include the maturase genes (parts that don't overlap main genes), and to exclude maturase regions from main genes


library('scales')
library('data.table')

Folder='2022_Jake_TAs' # 2022_Jake_TAs 2015_NoXlink 2023_Disome Wagner
Experiment='SelRP' # SelRP NoXlink Disome Wagner


# plotName = 'PPRs' 
# sampNames <-  c('Pet309','Aep3','Atp22','Cbp1','Aep1','Aep2', 'Pet111', 'Mrx1') 
# normNames = rep('Mrps17', 8)  # 'MrpS17'

plotName = 'COB_KOs' # PPRs MrpS17smt1d
sampNames <-  c('Mrx4', 'Mrx4_cbp3', 'Cbs1', 'Cbs1_mrx4', 'Mrps17_mrx4')
normNames = c('Mrps17', 'Mrps17_cbp3', 'Mrps17', 'Mrps17_mrx4', 'Mrps17')  # 'MrpS17'

region = 'all' # all 5UTR

path <- paste0('path_to_folder/',Folder,'/MitoRP/ReadLengthPosition/') 


genes = c('COX1','ATP8','ATP6' ,'COB' ,'ATP9' ,'VAR1' ,'COX2' ,'COX3', 'AI1', 'AI2', 'AI3', 'AI4', 'AI5a', 'AI5b', 'BI2', 'BI3', 'BI4','SCE1', 'Q0255')

UTR5starts = list(c(13580, 16435, 18954, 20508, 21995, 23612,25318,26229), 27620, 28300, c(36440,37723, 39141, 40841, 42508, 43297), 46552, 48801, 73700, 79113, 13987, 16471, 18992, 20985,22247,24156,37737, 39218, 41091, 61022, 74514)

starts = list(c(13818, 16435, 18954, 20508, 21995, 23612,25318,26229), 27666, 28487, c(36440,37723, 39141, 40841, 42508, 43297), 46723, 48901, 73758, 79213, 13987, 16471, 18992, 20985, 22247, 24156,37737, 39218, 41091, 61022, 74514)

ends = list(c(13986,16470,18991,20984,22246,23746,25342, 26701), 27812, 29266, c(36954, 37736, 39217, 41090, 42558, 43647), 46953, 50097, 74513, 80022, 16322, 18830,19996,21935, 23167, 25255,38579, 40265, 42251, 61729, 75984)
colors = c('indianred3','forestgreen','seagreen3' ,'dodgerblue' ,'seagreen2' ,'grey50' ,'indianred1' ,'indianred2', rep('violet', 6,), rep('aquamarine', 3), 'gray30', 'pink')

if (region == 'all') {
ends = ends }
if (region == '5UTR') {
ends = starts }

meannormlibtotalcountss=list()


for (i in c(1:length(sampNames))) {
	normName = normNames[i]

	if (normName == 'Mrps17') {
	normlibNames=c('Mrps17_1', 'Mrps17_2')}
	if (normName == 'Mrps17_Cbp3') {
	normlibNames=c('Mrps17_Cbp3_1', 'Mrps17_Cbp3_2')}
	if (normName == 'Mrps17_mrx4') {
	normlibNames=c('Mrps17_mrx4_1', 'Mrps17_mrx4_2')}

	normlibtotalcountss=list()
	for (normlibName in normlibNames) {

		normseries5pr <- fread(paste0(path,'bedFiles/', normlibName, '_Mito_mRNA_lengths5p_P.bed'), skip=1)
		normseries3pr <- fread(paste0(path,'bedFiles/', normlibName, '_Mito_mRNA_lengths3p_P.bed'), skip=1)

		normlibtotalcounts=c()
		for (m in c(1:length(genes))) {
		GOI=genes[m]

		normlibtotalcount = nrow(normseries5pr[V2 %in% unlist(mapply(seq, UTR5starts[[m]], ends[[m]], SIMPLIFY = FALSE))])/nrow(normseries5pr)*1000
	
	
	
	
		normlibtotalcounts=c(normlibtotalcounts, normlibtotalcount)
		}
		normlibtotalcountss[[length(normlibtotalcountss)+1]] <- normlibtotalcounts
	}

	print(normlibtotalcountss)
# foldchangess = list()
# logfoldchangess = list()
# allfoldchanges = c()
# alllogfoldchanges = c()
# cols = c()

# for (i in c(1:length(sampNames))) {
	sampName = sampNames[i]

	if (sampName == 'Mrps17') {
	libNames=c('Mrps17_1', 'Mrps17_2')}
	if (sampName == 'Pet309') {
	libNames=c('Pet309_1', 'Pet309_2')}
	if (sampName == 'Atp22') {
	libNames=c('Atp22_1', 'Atp22_2', 'Atp22_3')}
	if (sampName == 'Cbs1') {
	libNames=c('Cbs1_1', 'Cbs1_2', 'Cbs1_3', 'Cbs1_4')}
	if (sampName == 'Cbs2') {
	libNames=c('Cbs2_1', 'Cbs2_2', 'Cbs2_3', 'Cbs2_4')}
	if (sampName == 'Cbp1') {
	libNames=c('Cbp1_1', 'Cbp1_2')}
	if (sampName == 'Aep1') {
	libNames=c('Aep1_1', 'Aep1_2', 'Aep1_3', 'Aep1_4')}
	if (sampName == 'Aep2') {
	libNames=c('Aep2_1', 'Aep2_2', 'Aep2_3', 'Aep2_4', 'Aep2_5', 'Aep2_6')}
	if (sampName == 'Aep3') {
	libNames=c('Aep3_1', 'Aep3_2', 'Aep3_3')}
	if (sampName == 'Pet111') {
	libNames=c('Pet111_1', 'Pet111_2')}
	if (sampName == 'Mrx4') {
	libNames=c('Mrx4_1', 'Mrx4_2')}
	if (sampName == 'Mrx4_cbp3') {
	libNames=c('Mrx4_cbp3_1', 'Mrx4_cbp3_2','Mrx4_cbp3_3')}
	if (sampName == 'Mrps17_mrx4') {
	libNames=c('Mrps17_mrx4_1', 'Mrps17_mrx4_2')}
	if (sampName == 'Mrps17_cbp3') {
	libNames=c('Mrps17_Cbp3_1', 'Mrps17_Cbp3_2')}
	if (sampName == 'Mrx4') {
	libNames=c("Mrx4_1", "Mrx4_2")}
	if (sampName == 'Cbs1_mrx4') {
	libNames=c('Cbs1_mrx4_1', 'Cbs1_mrx4_2')}

	libNum=length(libNames)

	libtotalcountss=list()
	for (libName in libNames) {
		series5pr <- fread(paste0(path,'bedFiles/', libName, '_Mito_mRNA_lengths5p_P.bed'), skip=1)
		series3pr <- fread(paste0(path,'bedFiles/', libName, '_Mito_mRNA_lengths3p_P.bed'), skip=1)

		libtotalcounts=c()
		for (k in c(1:length(genes))) {
		GOI=genes[k]

		libtotalcount = nrow(series5pr[V2 %in% unlist(mapply(seq, UTR5starts[[k]], ends[[k]], SIMPLIFY = FALSE))])/nrow(series5pr)*1000
	
		libtotalcounts=c(libtotalcounts, libtotalcount)
		}
		libtotalcountss[[length(libtotalcountss)+1]] = libtotalcounts
	} # Finish going through each replicate

	# Make a data table for each sample with replicates as columns 
	DT_samp <- data.table::copy(libtotalcountss) # make DT from list
	setDT(DT_samp)
	assign(paste0(sampName,'DT_samp'), DT_samp)




	# Get means for each gene across the norm sample replicates
	DT_norm <- data.table::copy(normlibtotalcountss) # make DT from list
	meannormlibtotalcounts = rowMeans(setDT(DT_norm))
	meannormlibtotalcountss[[length(meannormlibtotalcountss) + 1]] <- meannormlibtotalcounts
} # Finish going through each sample

print('meannormlibtotalcountss')
print(meannormlibtotalcountss)



# set up matrix

meanss=list()
sdss=list()
sess=list()
rangess=list()
log2sampss=list()
log2meanss=list()
log2sdss=list()
log2sess=list()
log2rangess=list()
for (i in c(1:length(genes))) {
	# Get value for each TA for this gene
	print(genes[i])
	samps = list()
	for (j in c(1:length(sampNames))) {
		sampName=sampNames[j]	
		print(sampName)
		print('meannormlibtotalcountss')
		print(meannormlibtotalcountss[[j]][i])
		GOI = get(paste0(sampName,'DT_samp'))[i,]/meannormlibtotalcountss[[j]][i]
		vec = as.vector(as.matrix(GOI[1]))
		print(vec)
		samps[[length(samps) + 1]] <- vec
	}
	# Means are enrichment over Mrps17
	means = sapply(samps, mean)
	meanss[[length(meanss) + 1]] <- means
	sds = sapply(samps, sd)
	sdss[[length(sdss) + 1]] <- sds
	ses = sds/sqrt(sapply(samps, length))
	sess[[length(sess) + 1]] <- ses
	ranges = sapply(samps, range)
	rangess[[length(rangess) + 1]] <- ranges 

	# log values (here have to calc each then do mean at end
	log2samps = lapply(samps, log2)
	log2sampss[[length(log2sampss) + 1]] <- log2samps
	log2sds = sapply(log2samps, sd)
	log2sdss[[length(log2sdss) + 1]] <- log2sds
	log2ses = log2sds/sqrt(sapply(log2samps, length))
	log2sess[[length(log2sess) + 1]] <- log2ses
	log2ranges = sapply(log2samps, range)
	log2rangess[[length(log2rangess) + 1]] <- log2ranges 
	log2means = sapply(log2samps, mean)
	log2meanss[[length(log2meanss) + 1]] <- log2means 
}
sampNum = length(sampNames)

matmeans = matrix(unlist(meanss), ncol = sampNum, byrow = TRUE)
matsds = matrix(unlist(sdss), ncol = sampNum, byrow = TRUE)
matses = matrix(unlist(sess), ncol = sampNum, byrow = TRUE)
matranges = matrix(unlist(rangess), ncol = sampNum*length(genes), byrow = FALSE)

log2matmeans = matrix(unlist(log2meanss), ncol = sampNum, byrow = TRUE)
log2matsds = matrix(unlist(log2sdss), ncol = sampNum, byrow = TRUE)
log2matses = matrix(unlist(log2sess), ncol = sampNum, byrow = TRUE)
log2matranges = matrix(unlist(log2rangess), ncol = sampNum*length(genes), byrow = FALSE)


# Plot
pdf(paste0(path, 'Counts/Enrichment_transcriptsum_barplot_',region,'_',plotName,'.pdf'),width=sampNum*length(genes)/10,height=4.5,pointsize=8) # height=17

if (plotName == 'PPRs') {
	if (plotName == 'Mrps17_KOs') {
	ylimits = c(-1,3.5) }
	if (region == 'all') {
	ylimits = c(-5,7) }
	if (region == '5UTR') {
	ylimits = c(-5, 8) }
	}
if (plotName == 'Mrps17_KOs') {
ylimits = c(-1,3.5) }
if (plotName == 'Smt1') {
ylimits = c(-1,2) }
if (plotName == 'COB_KOs') {
ylimits = c(-2,3) }
if (plotName == 'COBfactors') {
ylimits = c(-7,7) }

xx=barplot(log2matmeans, names.arg=paste0(sampNames), beside = TRUE, col=colors, border=NA, ylab='Enrichment over respective MrpS17 (log2)', legend.text = c(genes), args.legend = list(x = 'bottomleft', bty = 'n', fill = colors, cex = .9, border='white'), ylim = ylimits)

# Add separator line
abline(v=xx[1,][-1]-1)
# Add error bars
# se
arrows(xx, log2matmeans-log2matses,xx, log2matmeans+log2matses, length=0.05, angle=90, code=3, lwd=.5)
# range
# arrows(t(xx), log2matranges[1,], t(xx), log2matranges[2,], length=0.05, angle=90, code=3, lwd=.5)


# And the other way, with genes along the x axis
colors2 = rep(c('indianred3','forestgreen','yellowgreen' ,'dodgerblue' ,'seagreen2','grey50','indianred1' ,'indianred2' ), each=sampNum)

yy=barplot(t(log2matmeans), names.arg=paste0(genes), beside = TRUE, col=colors2, border=NA, ylab='Enrichment over MrpS17 (log2)', legend.text = c(sampNames), args.legend = list(x = 'bottomleft', bty = 'n', fill = NULL, cex = .9, border='white'), ylim = ylimits)

# Add error bars
# se
arrows(yy, t(log2matmeans)-t(log2matses),yy, t(log2matmeans)+t(log2matses), length=0.05, angle=90, code=3, lwd=.5)
# range
# arrows(t(xx), log2matranges[1,], t(xx), log2matranges[2,], length=0.05, angle=90, code=3, lwd=.5)

dev.off()
# 

# source('/Users/Mary/Desktop/Data/SelectiveRP/Scripts/CountReadsInRegion_includeMaturases.R')

