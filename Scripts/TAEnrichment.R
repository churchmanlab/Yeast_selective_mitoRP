

library('scales')
library('data.table')

Folder='2022_Jake_TAs' # 2022_Jake_TAs 2023_Disome
# Experiment='Disome' # SelRP Disome
sampName <- 'Mrps17_cbp3' # Pet309 Atp22 Smt1 Cbs1 Cbs2 Cbp1 Aep1 Aep2 InpAep1 InpAep2 Mrx1 Mrps17smt1d Atp22smt1d Mrps17 Pet111 Aep3 Di40
data1 <- 'enrich' # enrich rpm(if rpm, log2 values are taken to fit on y axis) 
normSamp <- 'Mrps17' # Mrps17 Mono
windowtype <- 'Roll' # 'Roll' for sliding window, '' for fixed
window <- 9
maxzeros <- 3
plottype = 'h' # 'h' 'l'
smooth = 'no' # yes no both
smoothness = 30 # bigger is smoother
datatrans <- 'log' # log nolog
zoom = 'none' #  none  5end CDS noneForCrossCorr
offset = 'no' # Offset plot to see where signal lines up with domains (e.g. for Oxa1)
offsetamount = 90 
drawpauses = 'no' # these would be ad-hoc pauses, for carefully-called pauses see 'plotextra' below
plothisttoo = 'no'
linecolor = 'sandybrown' # indianred3 dodgerblue should be default
linewidth = .5 # 1.5 should be default
if (datatrans == 'log') {
	ylimits = c(-3, 3) # c(-5,12) c(-2,2)
}
if (datatrans == 'nolog') {
	ylimits = NULL
}

structure = 'UP_transmem' # UP_transmem (uniprot transmembrane) TMHMM2.0
plotraw = 'no'
compare = 'no' # 
doublenorm = 'no' # Make new plot with the normalized sample normalized again to the comparison sample (e.g. (Cbs1;mrx4delta/Mrps17;mrx4delta)/(Cbs1/Mrps17)) * use the difference between log2 values
plotcompare = 'no'
toplot = 'none' #'hydrophobicity_KD' none charge hydrophobicity_HW chargeANDhydrophobicity_KD

if (compare == 'yes') {
Folder2='2022_Jake_TAs'
data2='enrich'
sampName2 <- 'Cbs1'
normSamp2 <- 'Mrps17' # Mrps17 Mono
window2 <- 9 # 9 10

ylimits2 = c(-6.5,5)
if (datatrans == 'nolog') {
	ylimits2 = c(-.1,.1)}
}

if (smooth == 'no') {opacity = c(1,0)} # opacity of lines/bars
if (smooth == 'yes') {opacity = c(0,1)} # opacity of lines/bars
if (smooth == 'both') {opacity = c(1,1)} # opacity of lines/bars


if (zoom == 'none' | zoom == 'noneForCrossCorr' ) {
u = 300 # upstream
d = 1600 # downstream
}
if (zoom == '5end') {
u = 250 # upstream
d = 200 # downstream
}
if (zoom == 'CDS') {
u = 50 # upstream
d = 800 # downstream
}


path <- paste0('path_to_folder/',Folder,'/MitoRP/bedGraphs/') 

# Fixed window:
# file <- paste0(sampName,'_Mito_mRNA.noDups.Asite_37to41_PAll.rpm_',normSamp,'_norm',windowtype,'Win',window,'.bedGraph')
if (data1 == 'rpm') {
# Not enrichment
file <- paste0(sampName,'_Mito_mRNA.noDups.Asite_37to41_PAll.rpm.bedGraph') # combined reps
# file <- paste0(sampName,'_Mito_mRNA.noDups.Asite_37to41_P.rpmAll.txt') # ind reps


bgDT <- data.table(read.table(paste0(path, file),skip=1))
bgDT$V4[bgDT$V4 !=0] <- log(bgDT$V4[bgDT$V4 !=0],2)
} else {
# Rolling window:
file <- paste0('enrichment/',sampName,'_Mito_mRNA.noDups.Asite_37to41_PAll.rpm_',normSamp,'_norm',windowtype,'Win',window,'_max',maxzeros,'zeros.bedGraph')
bgDT <- data.table(read.table(paste0(path, file),skip=1))
}


if (compare == 'yes') {
	path2 <- paste0('path_to_folder/',Folder2,'/MitoRP/bedGraphs/enrichment/') 
	if (data2 == 'rpm') {
	# Not enrichment
	file2 <- paste0(sampName2,'_Mito_mRNA.noDups.Asite_37to41_PAll.rpm.bedGraph') # combined reps
	# file <- paste0(sampName,'_Mito_mRNA.noDups.Asite_37to41_P.rpmAll.txt') # ind reps
	bgDT2 <- data.table(read.table(paste0(path, file2),skip=1))
	bgDT2$V4[bgDT2$V4 !=0] <- log(bgDT2$V4[bgDT2$V4 !=0],2)
	} else {
	# Rolling window
	file2 <- paste0(sampName2,'_Mito_mRNA.noDups.Asite_37to41_PAll.rpm_',normSamp2,'_norm',windowtype,'Win',window2,'_max',maxzeros,'zeros.bedGraph')
	# Fixed window
	# file2 <- paste0(sampName2,'_Mito_mRNA.noDups.Asite_37to41_PAll.rpm_',normSamp2,'_normWin',window2,'.bedGraph')
	bgDT2 <- data.table(read.table(paste0(path2, file2),skip=1))
	}
}

# Extra bedgraph
plotextra = 'no'
if (plotextra == 'yes') {
	eDT <- fread('anotherFeatureFile.bedGraph', skip=1)
	# Fill in all positions so intron removal and transformation works later
	fillDT=data.table(V1='chrM', V3=c(1:85779))
	exDT=merge(eDT, fillDT, by='V3', all=TRUE)
	extraDT=data.table(V1='chrM', V2=exDT$V3-3, V3=exDT$V3, V4=exDT$V4)
	extraDT$V4[is.na(extraDT$V4)] <- 0
	extracolor = 'dodgerblue'
}

# Gene structure data
genemodelDT <- data.table(read.table('./SeqFiles/S288C_refseq_chrAll_modMaturases.bed'))
genemodelDT <- genemodelDT[V1=='chrM']



genes = c('COX1','ATP8','ATP6' ,'COB' ,'OLI1' ,'VAR1' ,'COX2' ,'COX3')

# Transform coordinates to remove introns

# Get intron info for each gene: length and positions, new coordinates for each gene
sumintrons = 0
intronposs = c() # Make vector of all positions in introns for later
starts = c()
ends = c()
intronlocss=c()

for (gene in genes) { 
	tmpDT <- genemodelDT[genemodelDT$V4 == gene]
	start = min(tmpDT$V2) - sumintrons

	exonss=c()
	for (i in c(1:nrow(tmpDT))) {
		exons=c(tmpDT[i]$V2:tmpDT[i]$V3)
		exonss=c(exonss, exons)
	}
	genelen = max(tmpDT$V3) - min(tmpDT$V2) + 1
	exonlen = sum(tmpDT$V7)
	intronlen = genelen - exonlen
	`%notin%` <- Negate(`%in%`)
	intronpos = c(min(tmpDT$V2):max(tmpDT$V3))[c(min(tmpDT$V2):max(tmpDT$V3)) %notin% exonss]
	intronposs = c(intronposs, intronpos)
	sumintrons = sumintrons + intronlen
	end = max(tmpDT$V3) - sumintrons

	starts = c(starts, start)
	ends = c(ends, end)
	# Keep track of where introns were
	len=0 # running total length of exons
	intronlocs = c()
	if (nrow(tmpDT) > 1) {
		for (i in c(1:(nrow(tmpDT)-1))) {
			len=len+tmpDT[i]$V7
			intronlocs = c(intronlocs,(start + len))
		}
	} else {intronlocs = NULL}

	intronlocss = c(intronlocss, intronlocs)
}

# Original positions, without removing introns
# starts = c(13818, 27666, 28487, 36540, 46723, 48901, 73758, 79213)
# ends = c(26701, 27812, 29266, 43647, 46953, 50097, 74513, 80022)

# Take out introns from bg
trans_bgDT <- bgDT[bgDT$V3 %notin% intronposs]
trans_extraDT <- extraDT[extraDT$V3 %notin% intronposs]
# Renumber bg
trans_bgDT[, V3 := c(1:nrow(trans_bgDT))]
trans_bgDT[, V2 := c(0:(nrow(trans_bgDT)-1))]
trans_extraDT[, V3 := c(1:nrow(trans_extraDT))]
trans_extraDT[, V2 := c(0:(nrow(trans_extraDT)-1))]

if (compare == 'yes') {
	trans_bgDT2 <- bgDT2[bgDT2$V3 %notin% intronposs]
	trans_bgDT2[, V3 := c(1:nrow(trans_bgDT2))]
	trans_bgDT2[, V2 := c(0:(nrow(trans_bgDT2)-1))]
}

if (zoom == 'none') {
# 	pdf(paste0(path, 'Viz/',sampName,'_norm_',normSamp,'_Viewer_',structure,'_',windowtype,'window',window,'plotcompare_',plotcompare,'_', datatrans,'.pdf'),width=12,height=3,pointsize=8) # height=17
	pdf(paste0(path, 'Viz/',sampName,'_norm_',normSamp,'_Viewer_',structure,'_',windowtype,'window',window,'plotcompare_',plotcompare,'_', datatrans,'.pdf'),width=5,height=2,pointsize=8) # height=17
}
if (zoom == 'noneForCrossCorr') {
	pdf(paste0(path, 'Viz/',sampName,'_norm_',normSamp,'_Viewer_',structure,'_',windowtype,'window',window,'plotcompare_',plotcompare,'_', datatrans,'.pdf'),width=4,height=3,pointsize=8) # height=17
}
if (zoom == '5end') {
	pdf(paste0(path, 'Viz/',sampName,'_norm_',normSamp,'_Viewer_',structure,'_',windowtype,'window',window,'plotcompare_',plotcompare,'_5end_', datatrans,'.pdf'),width=6,height=3,pointsize=8) # height=17
}
if (zoom == 'CDS') {
	pdf(paste0(path, 'Viz/',sampName,'_norm_',normSamp,'_Viewer_',structure,'_',windowtype,'window',window,'plotcompare_',plotcompare,'_CDS_', datatrans,'.pdf'),width=6,height=3,pointsize=8) # height=17
}

# par(mfrow = c(length(genes)+1, 1))

# Set up lists to store values for plotting all together at end
xs=list()
ys=list()
if (compare == 'yes') {
	ys2=list()
}
# Set up vectors for plotting mean values in barplot at the end
aveenrichs = c()
ses = c()
sds =c()


for (i in c(1:length(genes))) { # to line 317

	print(i)

	x = c((starts[i]-u):(starts[i]+d))
	if (datatrans == 'log') {
	y = trans_bgDT$V4[trans_bgDT$V3 %in% x]
	} else if (datatrans == 'nolog') {
	y = trans_bgDT$V7[trans_bgDT$V3 %in% x]
	}

	xs[[length(xs) + 1]] <- c(1:1951) # Make x coords the same for all for plotting together
	ys[[length(ys) + 1]] <- y

	plotx = x
	main = paste0(genes[i], ' - ', sampName)

	if (offset == 'yes') {
		plotx <- x - offsetamount
		main = paste0(genes[i], ' - ', sampName, ' -',offsetamount,' offset')
	}

	if (compare == 'yes') {
	if (datatrans == 'log') {
	y2 = trans_bgDT2$V4[trans_bgDT2$V3 %in% x]
	} else if (datatrans == 'nolog') {
	y2 = trans_bgDT2$V7[trans_bgDT2$V3 %in% x]
	}
	ys2[[length(ys2) + 1]] <- y2
	}

	if (datatrans == 'log') {
	ylabel = paste0('Fold enrichment over ',normSamp,' (log2)') }
	if (datatrans == 'nolog') {
	ylabel = paste0('Enrichment over ',normSamp) }
	
	plot(plotx, y, type=plottype, col=alpha(linecolor, opacity[1]), xlab='', ylab=ylabel, ylim = ylimits, lwd=linewidth)
	title(main, line = -2)
	# Add start, stop, and intron lines
	abline(v=c(starts[i], ends[i]), col=c('green','red'), lwd = .1)
# 	abline(v=c(intronlocss), lwd = .1)
	segments(starts[i], 100, starts[i]+100, 100) # For scale bar
	
	if (smooth == 'yes' | smooth == 'both') {
		# Plot smoothed line
		# Get number of values and calculate span
		notna = which(!is.na(y))	
		spanval = smoothness/length(notna)
		# get loess
		lo <- loess(y~plotx, span=spanval, degree=2)
		smooth <- predict(lo)
		# Add back in NA values
		allsmooth = rep(NA, length(plotx))
		allsmooth[notna] <- smooth
		# Plot line
		lines(plotx,allsmooth, col=alpha('dodgerblue', opacity[2]), lwd=1.5)
	}
	
	# Highlight pauses
	# Get pause values, top 3% in gene
	if (drawpauses == 'yes') {
		pausevals = sort(y, decreasing=TRUE)[1:round(length(notna)*.03, 0)]
		pausemin = min(pausevals)
		tmp=rle(y>pausemin)
		cons = with(tmp, rep(lengths * values, lengths))
	# 	pausesind = which(y %in% pausevals) # simple values over x
		pausesind = which(cons >= 3) # at least 3 in a row > min pause value
		pauses = plotx[pausesind]
		# Draw lines
		abline(v=pauses, lwd = 2, col=alpha('dodgerblue', .3))
	}
	# Plot info from "extraDT"
	if (plotextra == 'yes') {
		abline(v=trans_extraDT$V3[trans_extraDT$V3 %in% plotx & trans_extraDT$V4 == 10], col = alpha(extracolor, .5), lwd=2)
	}
	if (compare == 'yes' & plotcompare == 'yes') {
	points(x, y2, col=alpha('red', opacity[1]), type=plottype, lwd=1.5)
	# Plot smoothed line
	# Get number of values and calculate span
	notna = which(!is.na(y2))	
	spanval = smoothness/length(notna) 
	# get loess
	lo <- loess(y2~plotx, span=spanval, degree=2)
	smooth <- predict(lo)
	# Add back in NA values
	allsmooth = rep(NA, length(plotx))
	allsmooth[notna] <- smooth
	# Plot line
	lines(plotx,allsmooth, col=alpha('red', opacity[2]), lwd=1.5)

	}
	
	if (plotraw == 'yes') {
	lines(x, trans_bgDT$V5[trans_bgDT$V3 %in% x], col=alpha('pink', .4))
	lines(x, trans_bgDT$V6[trans_bgDT$V3 %in% x], col=alpha('purple', .2))
	}
	# if (sampName == 'Oxa1' & offset == 'yes') {
	# abline(v=plotOxa, col=alpha('seagreen', .2))
	# }
	# abline(v=c(starts[i], ends[i]), col=c('green','red'))
	# abline(v=c(intronlocss), lwd = .1, col=alpha('grey60',.2))
	# fold=1.5
	abline(h=0, lty=1, lwd=.5)
	# if (datatrans == 'log') {
	# abline(h=log(fold,2), lty=2, lwd=.5)
	# }




	# Get mean values across ORF for barplot
	subx = x[x>starts[i] & x<ends[i]]
	suby = y[x>starts[i] & x<ends[i]]
	aveenrich = mean(suby, na.rm=TRUE)
	sd = sd(suby, na.rm=TRUE) # standard deviation
	se = sd(suby, na.rm=TRUE)/sqrt(na.omit(length(suby))) # standard error

	aveenrichs = c(aveenrichs, aveenrich)
	ses = c(ses, se)
	sds = c(sds, sd)

	# Now plot features
	# Protein structure data

	if (structure == 'TMHMM2.0') {

	featDT <- data.table(read.table(paste0('./SeqFiles/ProteinFeatures/TMHMM2.0/',genes[i],'.txt'),skip=6))
	# First convert coordinates to intronless genome coords
	featDT[, first := V4*3+starts[i]-3]
	featDT[, last := V5*3+starts[i]-3]

	# Get values for each category
	inside <- featDT[V3=='inside']
	ins = list()
	if (nrow(inside) >0) {
	for (j in c(1:nrow(inside))) {
	mat = c(inside[j]$first:inside[j]$last)
	ins[[ (length(ins) + 1) ]] <- mat 
	}}
	outside <- featDT[V3=='outside']
	outs = list()
	if (nrow(outside) >0) {
	for (j in c(1:nrow(outside))) {
	cyt = c(outside[j]$first:outside[j]$last)
	outs[[ (length(outs) + 1) ]] <- cyt 
	}}
	TMh <- featDT[V3=='TMhelix']
	TMs = list()
	if (nrow(TMh) >0) {
	for (j in c(1:nrow(TMh))) {
	TM = c(TMh[j]$first:TMh[j]$last)
	TMs[[ (length(TMs) + 1) ]] <- TM 
	}}

	}





	if (structure == 'UP_transmem') {
	featDTtrans <- data.table(read.table('./SeqFiles/ProteinFeatures/UniProt_UP000002311_559292_transmem.bed',sep='\t'))
	featDTtrans <- featDTtrans[featDTtrans$V1 == 'chrM']
	featDTtrans <- featDTtrans[grep('^P', featDTtrans$V4)]

	featDTtop <- data.table(read.table('./SeqFiles/ProteinFeatures/UniProt_UP000002311_559292_topo_dom.bed',sep='\t'))
	featDTtop <- featDTtop[featDTtop$V1 == 'chrM']
	featDTtop <- featDTtop[grep('^P', featDTtop$V4)]
	featDTmat <- featDTtop[grep('matrix', featDTtop$V14)]
	featDTims <- featDTtop[grep('intermembrane', featDTtop$V14)]
	featDTother <- featDTtop[grep('other', featDTtop$V14)]

	# Convert coordinates
	dts = c('featDTtrans', 'featDTmat', 'featDTims', 'featDTother')
	for (p in c(1:length(dts))) {
	dt=get(dts[p])
	firsts=c()
	lasts=c()
	for (o in c(1:nrow(dt))) {
	srg=c(1:dt[o]$V2) # firsts range
	lrg=c(1:dt[o]$V3) # lasts range
	first = length(srg[srg %notin% intronposs])
	firsts = c(firsts, first)
	last = length(lrg[lrg %notin% intronposs])
	lasts = c(lasts, last)
	}
	dt[, first := firsts]
	dt[, last := lasts]

	assign(dts[p], dt)
	}

	# Get values

	ins = list()
	for (j in c(1:nrow(featDTmat))) {
	mat = c(featDTmat[j]$first:featDTmat[j]$last)
	ins[[ (length(ins) + 1) ]] <- mat 
	}

	outs = list()
	for (j in c(1:nrow(featDTims))) {
	out = c(featDTims[j]$first:featDTims[j]$last)
	outs[[ (length(outs) + 1) ]] <- out 
	}

	others = list()
	for (j in c(1:nrow(featDTother))) {
	other = c(featDTother[j]$first:featDTother[j]$last)
	others[[ (length(others) + 1) ]] <- other 
	}

	TMs = list()
	for (j in c(1:nrow(featDTtrans))) {
	TM = c(featDTtrans[j]$first:featDTtrans[j]$last)
	TMs[[ (length(TMs) + 1) ]] <- TM 
	}
	}


	# plot on graph
	linewidth = .6
	trans = .1
	incol = 'orange'
	outcol= 'cyan3'
	TMcol = 'grey60'
	othercol = 'grey30'
	yaxis=min(ylimits)/4*3.7

	# First black line for CDS
	segments(starts[i], yaxis, ends[i], yaxis, col='black', lwd=2, lend=2)
	if (structure != 'none') {
		if (length(ins) > 0 ) {
		for (k in c(1:length(ins))) {
		segments(ins[[k]][1], yaxis, ins[[k]][length(ins[[k]])],yaxis, col=incol, lwd=2, lend=2)
		}}
		if (length(outs) > 0 ) {
		for (k in c(1:length(outs))) {
		segments(outs[[k]][1], yaxis, outs[[k]][length(outs[[k]])],yaxis, col=outcol, lwd=2,lend=2)
		}}

		if (length(others) > 0 ) {
		for (k in c(1:length(others))) {
		segments(others[[k]][1], yaxis, others[[k]][length(others[[k]])],yaxis, col=othercol, lwd=2,lend=2)
		}}

		if (length(TMs) > 0 ) {
		for (k in c(1:length(TMs))) {
		segments(TMs[[k]][1], yaxis, TMs[[k]][length(TMs[[k]])],yaxis, col=TMcol, lwd=5,lend=3)
		}}
		}
	legend('topleft', c(paste0(structure), 'matrix', 'IMS', 'TM domain', 'other'), col = c('white', c(incol, outcol, TMcol, othercol)) ,lwd = 5,seg.len=1, box.lty = 0, cex=1.3)
	if (plotraw == 'yes') {
	legend('topright', c(paste0('log10(',normSamp, ' raw counts)'),paste0('log10(',sampName, ' raw counts)')), text.col = c('purple','pink'),cex=1.3,box.lty = 0)
	}
	if (plotcompare == 'yes') {
		if (data1 == 'rpm') {
		legend('topright', c(paste0(sampName,' rpm (log2)'),paste0(sampName2,' norm to ', normSamp2)), text.col = alpha(c('dodgerblue', 'red'), max(opacity)),cex=1.3,box.lty = 0)
		} else if (data1 == 'rpm' & data2 == 'rpm') {
		legend('topright', c(paste0(sampName,' rpm (log2)'),paste0(sampName,' rpm (log2)')), text.col = alpha(c('dodgerblue', 'red'), max(opacity)),cex=1.3,box.lty = 0)
		} else if (data2 == 'rpm') {
		legend('topright', c(paste0(sampName,' norm to ', normSamp),paste0(sampName,' rpm (log2)')), text.col = alpha(c('dodgerblue', 'red'), max(opacity)),cex=1.3,box.lty = 0)
		} else {
		legend('topright', c(paste0(sampName,' norm to ', normSamp),paste0(sampName2,' norm to ', normSamp2)), text.col = alpha(c('dodgerblue', 'red'), max(opacity)),cex=1.3,box.lty = 0)
		}
		
	}

	if (toplot == 'charge' | toplot == 'chargeANDhydrophobicity_KD') {
	# From here: https://www.bioinformatics.nl/cgi-bin/emboss/charge
	chargeDT <- data.table(read.table(paste0('./SeqFiles/ProteinFeatures/Charge_',genes[i],'.txt'),sep='\t'))

	lines(c(starts[i]: (ends[i]-3)), rep(10*(chargeDT$Charge), each=3), col=alpha('gold1', .3))
	}

	if (toplot == 'hydrophobicity_KD' | toplot == 'chargeANDhydrophobicity_KD') {
	# From here: https://www.novoprolabs.com/tools/protein-hydropathy
	hydroDT <- data.table(read.table(paste0('./SeqFiles/ProteinFeatures/HydophobicityScores_',genes[i],'.csv'),sep=','))

	lines(c((starts[i]+12): (ends[i]-15)), rep(hydroDT$V2, each=3), col=alpha('black', .3))
	legend('bottomleft', c('K-D hydrophobicity'), col = alpha('black', .3) ,lwd = 1.8,seg.len=1, box.lty = 0, cex=1.3)

	if (genes[i] != 'ATP8') {
	removefirst=100
	scattery=rep(hydroDT$V2, each=3)[(removefirst+1):length(rep(hydroDT$V2, each=3))] # hydro values
	scatterx=y[(u+1+12+removefirst): (u+removefirst+length(scattery))] # TA values
		if (offset == 'yes') {
		scatterx=y[(u+1+12+offsetamount+removefirst): (u+12++removefirst+length(scattery)+offsetamount)]
		}
	# get correlation
	pearson_r = cor.test(~scatterx+scattery, method = c('pearson')) 
	r=pearson_r$estimate
	mylabel = bquote(italic(r) == .(format(pearson_r$estimate, digits = 2)))
	mylabel2 = bquote(R^2 == .(format(pearson_r$estimate^2, digits = 3)))

	plot(scatterx,scattery, pch=16, cex=.7, col=alpha('black', .7), xlab='Enrichment', ylab='K-D hydrophobicity')
	text(max(scatterx, na.rm=TRUE)*(4/5),max(scattery, na.rm=TRUE)*(4/5), labels = mylabel, cex=2)

	# And plot cross-correlation
	scatterx <- scatterx[!is.na(scatterx) & !is.na(scattery)]
	scattery <- scattery[!is.na(scatterx) & !is.na(scattery)]
	ccf(scatterx, scattery, ylim=c(-.25,.25), main=paste0('Cross-correlation: Oxa1 and hydrophobicity ', genes[i]), ylab='CCF')
	}
	}

	if (toplot == 'hydrophobicity_HW') {
	# From here: https://web.expasy.org/protscale/ window=11
	hydroDT <- fread(paste0('./SeqFiles/ProteinFeatures/HydophobicityScores_HW_',genes[i],'.txt'))

	lines(c((starts[i]+15): (ends[i]-18)), rep(hydroDT$V2, each=3), col=alpha('black', .3))
	legend('bottomleft', c('H-W hydrophilicity'), col = alpha('black', .3) ,lwd = 1.8,seg.len=1, box.lty = 0, cex=1.3)

	if (genes[i] != 'ATP8') {
	removefirst=100
	scattery=rep(hydroDT$V2, each=3)[(removefirst+1):length(rep(hydroDT$V2, each=3))] 
	scatterx=y[(u+1+15+removefirst): (u++removefirst+length(scattery))]
		if (offset == 'yes') {
		scatterx=y[(u+1+15+offsetamount+removefirst): (u+15+length(scattery)+removefirst+offsetamount)]
		}

	# get correlation
	pearson_r = cor.test(~scatterx+scattery, method = c('pearson')) 
	r=pearson_r$estimate
	mylabel = bquote(italic(r) == .(format(pearson_r$estimate, digits = 2)))
	mylabel2 = bquote(R^2 == .(format(pearson_r$estimate^2, digits = 3)))

	plot(scatterx,scattery, pch=16, cex=.7, col=alpha('black', .7), xlab='Enrichment', ylab='H-W hydrophilicity')
	text(max(scatterx, na.rm=TRUE)*(4/5),max(scattery, na.rm=TRUE)*(4/5), labels = mylabel, cex=2)


	# And plot cross-correlation
	scatterx <- scatterx[!is.na(scatterx) & !is.na(scattery)]
	scattery <- scattery[!is.na(scatterx) & !is.na(scattery)]
	ccf(scatterx, scattery, ylim=c(-.4,.4), main=paste0('Cross-correlation with hydrophilicity ', genes[i]))
	}

	}


	# Plot again with greater y limits (ylimits2) to include upstream peak
	if (plothisttoo == 'yes') {
		plot(plotx, y, type='h', col=alpha('red', max(opacity)), xlab='', ylab=ylabel, ylim = ylimits2, lwd=1.5)
		abline(v=c(starts[i], ends[i]), col=c('green','red'), lwd = .1)
		# And gene models
		segments(starts[i], yaxis, ends[i], yaxis, col='black', lwd=2, lend=2)

		if (length(ins) > 0 ) {
		for (k in c(1:length(ins))) {
		segments(ins[[k]][1], yaxis, ins[[k]][length(ins[[k]])],yaxis, col=incol, lwd=2, lend=2)
		}}
		if (length(outs) > 0 ) {
		for (k in c(1:length(outs))) {
		segments(outs[[k]][1], yaxis, outs[[k]][length(outs[[k]])],yaxis, col=outcol, lwd=2,lend=2)
		}}

		if (length(others) > 0 ) {
		for (k in c(1:length(others))) {
		segments(others[[k]][1], yaxis, others[[k]][length(others[[k]])],yaxis, col=othercol, lwd=2,lend=2)
		}}

		if (length(TMs) > 0 ) {
		for (k in c(1:length(TMs))) {
		segments(TMs[[k]][1], yaxis, TMs[[k]][length(TMs[[k]])],yaxis, col=TMcol, lwd=5,lend=3)
		}}
	}


}  # Finish going through each gene


colors=rainbow(8)




if (doublenorm == 'yes') {
		if (datatrans == 'log') {
	pdf(paste0(path, 'Viz/',sampName,'_norm_',normSamp,'_vs_',sampName2,'_norm_',normSamp2, '_DoubleNormViewer_',structure,'_',windowtype,'window',window,'_', datatrans,'.pdf'),width=9,height=3,pointsize=8) # height=17
		} else if (datatrans == 'nolog') {
	pdf(paste0(path, 'Viz/',sampName,'_norm_',normSamp,'_vs_',sampName2,'_norm_',normSamp2, '_DoubleNormViewer_',structure,'_',windowtype,'window',window,'_', datatrans,'_nolog.pdf'),width=9,height=3,pointsize=8) # height=17
	}




	for (i in c(1:length(genes))) {

	x = c((starts[i]-u):(starts[i]+d))
	if (datatrans == 'log') {
	y = trans_bgDT$V4[trans_bgDT$V3 %in% x]
	}
	if (datatrans == 'nolog') {
	y = trans_bgDT$V7[trans_bgDT$V3 %in% x]
	}
	xs[[length(xs) + 1]] <- c(1:1951) # Make x coords the same for all for plotting together
	ys[[length(ys) + 1]] <- y

	if (compare == 'yes') {
	if (datatrans == 'log') {
	y2 = trans_bgDT2$V4[trans_bgDT2$V3 %in% x]
	} else if (datatrans == 'nolog') {
	y2 = trans_bgDT2$V7[trans_bgDT2$V3 %in% x]
	}
	ys2[[length(ys2) + 1]] <- y2
	}

	# Create a logical index for positions where either vector has a value of 0
	# zero_index <- y == 0 | y2 == 0

	# Subtract vectors, marking NA where either vector has a value of 0
	if (datatrans == 'log') {
	y3 = y - y2} else if (datatrans == 'nolog') {
	y3 = y/y2
	}
	# y3[zero_index] <- NA


	if (datatrans == 'log') {
	ylabel2 = paste0('Fold enrichment over ',sampName2,':',normSamp2,' (log2)') }
	if (datatrans == 'nolog') {
	ylabel2 = paste0('Enrichment over ',sampName2,':',normSamp2)}

	plot(x, y3, type='h', col='orange', xlab='', ylab=ylabel2, ylim = ylimits2, lwd=1.5)
	title('Double normalized', line = -1)
	title(paste0(genes[i], ' - ', sampName,':',normSamp), line = -2)


	abline(v=c(starts[i], ends[i]), col=c('green','red'))
	abline(v=c(intronlocss), lwd = .1)
	fold=1.5
	abline(h=c(0,log(fold,2)), lty=c(1,2), lwd=.5)



	if (structure == 'UP_transmem') {
	# plot on graph
	linewidth = .6
	trans = .1
	incol = 'orange'
	outcol= 'cadetblue1'
	TMcol = 'grey50'
	othercol = 'grey30'
	yaxis=min(ylimits)/4*3.7
	if (length(ins) > 0 ) {
	for (k in c(1:length(ins))) {
	segments(ins[[k]][1], yaxis, ins[[k]][length(ins[[k]])],yaxis, col=incol, lwd=2, lend=2)
	}}
	if (length(outs) > 0 ) {
	for (k in c(1:length(outs))) {
	segments(outs[[k]][1], yaxis, outs[[k]][length(outs[[k]])],yaxis, col=outcol, lwd=2,lend=2)
	}}

	if (length(others) > 0 ) {
	for (k in c(1:length(others))) {
	segments(others[[k]][1], yaxis, others[[k]][length(others[[k]])],yaxis, col=othercol, lwd=2,lend=2)
	}}

	if (length(TMs) > 0 ) {
	for (k in c(1:length(TMs))) {
	segments(TMs[[k]][1], yaxis, TMs[[k]][length(TMs[[k]])],yaxis, col=TMcol, lwd=5,lend=3)
	}}

	legend('topleft', c(paste0(structure), 'matrix', 'IMS', 'TM domain', 'other'), col = c('white', c(incol, outcol, TMcol, othercol)) ,lwd = 5,seg.len=1, box.lty = 0, cex=1.3)

	}
	}
	dev.off()
} # Finish with double norm plotting

# And "metaplot"
# plot(xs[[1]], ys[[1]], type='l', col=alpha('dodgerblue', .4), xlab='', ylab=paste0('Fold enrichment over ',normSamp,' (log2)'), ylim = ylimits, lwd=1.5)
# title('metaplot', line = -2)
# for (i in c(2:8)) {
# lines(xs[[i]], ys[[i]], col=alpha('dodgerblue', .2))
# }
# fold=1.5
# abline(h=c(0,log(fold,2)), lty=c(1,2), lwd=.5)


# And barplot with mean enrichments

# Artificially expand vectors so I can plot in the same dimensions and it will look ok
add=11
aveenrichs=c(aveenrichs, rep(0,add))
sds=c(sds, rep(0,add))
ses=c(ses, rep(0,add))
genes_edit=c(genes, rep('',add))

xx=barplot(aveenrichs, ylim=c(min(aveenrichs, na.rm=TRUE)-.2, max(aveenrichs, na.rm=TRUE)+.2), names.arg=genes_edit, col=c('indianred', 'forestgreen', 'forestgreen','dodgerblue','forestgreen', 'gold2', 'indianred'), border=NA, ylab='Mean within ORF')
title(sampName, line = -1, adj=.2)

# standard error
arrows(xx, aveenrichs-ses/2, xx, aveenrichs+ses/2, length=0.05, angle=90, code=3, lwd=.2)
# standard deviation
# arrows(xx, aveenrichs-sds/2, xx, aveenrichs+sds/2, length=0.05, angle=90, code=3, lwd=.2)

dev.off()




# source('/Users/Mary/Desktop/Data/SelectiveRP/Scripts/TAEnrichment.R')

