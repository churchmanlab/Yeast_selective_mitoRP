

library('scales')
library('data.table')

Folder='TAs'
Experiment='SelRP' 
sampName <- 'Mrps17' 
downsample <- 'yes' # no
normSamp <- 'none' # Mrps17 none # This is only used if plotting counts at the end of script is uncommented
zoom = 'UTR' # start none stop widestart 5end ProxUTR UTR
ylimits = c(25,50)
if (sampName == 'Di80') {
ylimits = c(25,99)
}
if (sampName == '40S' | sampName == '80S') {
ylimits = c(15,60)
}
filetype = 'pdf'


if (zoom == 'start' | zoom == 'stop') {
u = 80 # upstream
d = 150 # downstream
}
if (zoom == 'widestart') {
u = 350 # upstream
d = 100 # downstream
}
if (zoom == 'none') {
u = 180 # upstream
d = 2000 # downstream
}
if (zoom == '5end') {
u = 55 # upstream
d = 200 # downstream
}
if (zoom == 'UTR') {
u = 1000 # upstream
d = 10 # downstream
maxcount = 6000
}
if (zoom == 'ProxUTR') {
u = 55 # upstream
d = 40 # downstream
maxcount = 6000
}

path <- paste0('/Users/Mary/Desktop/Data/SelectiveRP/',Folder,'/MitoRP/ReadLengthPosition/') 

genes = c('COX1','ATP8','ATP6' ,'COB' ,'ATP9' ,'VAR1' ,'COX2' ,'COX3')
starts = c(13818, 27666, 28487, 36540, 46723, 48901, 73758, 79213)
ends = c(26701, 27812, 29266, 43647, 46953, 50097, 74513, 80022)

ranges=list()

if (zoom == 'start' | zoom == 'none' | zoom == 'widestart' | zoom == '5end' | zoom == 'UTR' | zoom == 'ProxUTR') {
for (i in c(1:length(starts))) {
range = c((starts[i]-u):(starts[i]+d))
ranges[[length(ranges) + 1]] <- range
}
}
if (zoom == 'stop') {
for (i in c(1:length(starts))) {
range = c((ends[i]-u):(ends[i]+d))
ranges[[length(ranges) + 1]] <- range
}
}

if (sampName == 'Mrps17') {
libNames=c('Mrps17_1', 'Mrps17_2')}
if (sampName == 'Aep1') {
libNames=c('Aep1_1', 'Aep1_2', 'Aep1_3', 'Aep1_4')}
if (sampName == 'Aep2') {
libNames=c('Aep2_1', 'Aep2_2', 'Aep2_3', 'Aep2_4', 'Aep2_5', 'Aep2_6')}

sampNum=length(libNames)

if (normSamp == 'Mrps17') {
normLibs=c('Mrps17_1', 'Mrps17_2')
normLib='Mrps17_2'}

if (normSamp == 'none') {
normLibs=''
normLib=''}


if (normSamp == 'Mrps17') {
normseries5pr <- scan(paste0(path,'bedFiles/', normLib, '_Mito_mRNA_lengths5p_P.bed'), list('',0,0,''))
normseries3pr <- scan(paste0(path,'bedFiles/', normLib, '_Mito_mRNA_lengths3p_P.bed'), list('',0,0,''))
}


# Get counts for the input/Mrps17 samples, summing replicates
# normdensities =c() # For library to normalize to
# 
# for (i in c(1:length(starts))) {
# positions=ranges[[i]]
# density=0 # reset read count for each gene
# 
# for (normLib in normLibs) {
# 
# x5prNorm=series5prNorm[[2]][series5prNorm[[2]] %in% positions]
# 
# density=length(x5prNorm)
# density=density+density
# }
# normdensities=c(normdensities, density)
# }
# 


##### First make a summed table and write to same directory as others
all5pr = data.table(NULL)
all3pr = data.table(NULL)
for (libName in libNames) {
	series5pr <- fread(paste0(path,'bedFiles/', libName, '_Mito_mRNA_lengths5p_P.bed'))
	series3pr <- fread(paste0(path,'bedFiles/', libName, '_Mito_mRNA_lengths3p_P.bed'))
	all5pr <- rbind(all5pr, series5pr)
	all3pr <- rbind(all3pr, series3pr)
}
all5pr <- all5pr[order(V2)]
all3pr <- all3pr[order(V2)]

write.table(all5pr, file=paste0(path, 'bedFiles/', sampName, '_Mito_mRNA_lengths5p_P.bed'), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(all3pr, file=paste0(path, 'bedFiles/', sampName, '_Mito_mRNA_lengths3p_P.bed'), row.names = FALSE, col.names = FALSE, quote = FALSE)
##### ##### ##### ##### ##### ##### 

# Add in the new summed sample to the list
libNames = c(libNames, sampName)





ptsize = .4


rphNorms = list() # list of the normalized counts for each gene for each library


for (libName in libNames) {

	series5pr <- scan(paste0(path,'bedFiles/', libName, '_Mito_mRNA_lengths5p_P.bed'), list('',0,0,''))
	series3pr <- scan(paste0(path,'bedFiles/', libName, '_Mito_mRNA_lengths3p_P.bed'), list('',0,0,''))
	# series5prMinus <- scan(paste0(path,'bedFiles/', libName, '_Mito_mRNA_lengths5p_M.bed'), list('',0,0,''))
	# series3prMinus <- scan(paste0(path,'bedFiles/', libName, '_Mito_mRNA_lengths3p_M.bed'), list('',0,0,''))

	densities=c() # For library of interest
	densitiesNorm =c() # Library of interest, normalized

	if (filetype == 'pdf') {
	if (zoom == 'start') {
	pdf(paste0(path, libName, '_5pr3pr_start.pdf'),width=5,height=2.2,pointsize=8)
	}
	if (zoom == 'stop') {
	pdf(paste0(path, libName, '_5pr3pr_stop.pdf'),width=3,height=2.2,pointsize=8)
	}
	if (zoom == 'none') {
	pdf(paste0(path, libName, '_5pr3pr_none.pdf'),width=12,height=2.2,pointsize=8)
	}
	if (zoom == 'widestart') {
	pdf(paste0(path, libName, '_5pr3pr_widestart.pdf'),width=9,height=2.2,pointsize=8)
	}
	if (zoom == '5end') {
	pdf(paste0(path, libName, '_5pr3pr_5end.pdf'),width=9,height=2.2,pointsize=8)
	}
	if (zoom == 'UTR') {
	pdf(paste0(path, libName, '_5pr3pr_UTR.pdf'),width=9,height=2.2,pointsize=8)
	}
	if (zoom == 'ProxUTR') {
	pdf(paste0(path, libName, '_5pr3pr_ProxUTR.pdf'),width=2.5,height=2.2,pointsize=8)
	}
	if (sampName == 'Di80') {
	pdf(paste0(path, libName, '_5pr3pr_5end.pdf'),width=9,height=3,pointsize=8)
	}


	}

	if (filetype == 'png') {
	if (zoom == 'start') {
	png(paste0(path, libName, '_5pr3pr_start.png'),width=3,height=16,pointsize=8,units = 'in', res = 300)
	}
	if (zoom == 'stop') {
	png(paste0(path, libName, '_5pr3pr_stop.png'),width=3,height=16,pointsize=8,units = 'in', res = 300)
	}
	if (zoom == 'none') {
	png(paste0(path, libName, '_5pr3pr_none.png'),width=6,height=8,pointsize=16,units = 'in', res = 300)
	}
	if (zoom == 'widestart') {
	png(paste0(path, libName, '_5pr3pr_widestart.png'),width=7,height=13,pointsize=8,units = 'in', res = 300)
	}
	if (zoom == '5end') {
	png(paste0(path, libName, '_5pr3pr_5end.png'),width=9,height=16,pointsize=8,units = 'in', res = 300)
	}
	if (zoom == 'UTR') {
	png(paste0(path, libName, '_5pr3pr_UTR.png'),width=9,height=16,pointsize=8,units = 'in', res = 300)
	if (zoom == 'ProxUTR') {
	png(paste0(path, libName, '_5pr3pr_ProxUTR.png'),width=2.5,height=2.2,pointsize=8,units = 'in', res = 300)
	}
	}

	par(mfrow=c(8,1))
	}

	# Set up lists for storing data for each gene
	for (i in c(1:length(starts))) {
		positions=ranges[[i]]

		x5pr=series5pr[[2]][series5pr[[2]] %in% positions | series3pr[[2]] %in% positions]
		y5pr=series5pr[[3]][series5pr[[2]] %in% positions | series3pr[[2]] %in% positions]
		x3pr=series3pr[[2]][series5pr[[2]] %in% positions | series3pr[[2]] %in% positions]
		y3pr=series3pr[[3]][series5pr[[2]] %in% positions | series3pr[[2]] %in% positions]

		# Save data for plotting summed values at the end
		x5prs[[length(x5prs) + 1 ]] <- x5pr
		y5prs[[length(y5prs) + 1 ]] <- y5pr
		x3prs[[length(x5prs) + 1 ]] <- x3pr
		y3prs[[length(y3prs) + 1 ]] <- y3pr

		counts = length(x5pr)

		if (downsample == 'yes' & counts > maxcount) {
		seed=sample(c(1:100))
		set.seed(seed)
		indexes = sample(c(1:length(x5pr)), maxcount)
		x5pr <- x5pr[indexes]
		y5pr <- y5pr[indexes]
		x3pr <- x3pr[indexes]
		y3pr <- y3pr[indexes]
		}

		if (normSamp != 'none') {
			x5prNorm=normseries5pr[[2]][normseries5pr[[2]] %in% positions]
			x3prNorm=normseries3pr[[2]][normseries3pr[[2]] %in% positions]
			if (downsample == 'yes' & length(counts) > maxcount) {
				x5prNorm <- x5prNorm[indexes]
				x3prNorm <- x3prNorm[indexes]
			}
		}

		# Number of mapppers in region
		density=length(x5pr)
		densities = c(densities, density)

		if (filetype == 'pdf') {
		transparency = 100/density #1000/density .5
		if (transparency > 1) { transparency = 1 }
		if (transparency < .05) { transparency = .05 }
		if (zoom == 'UTR') {
		transparency = 10000/density #1000/density 
		if (transparency > 1) { transparency = 1 }
		if (transparency < .1) { transparency = .1 }
		}
		}

		if (filetype == 'png') {
		transparency = 50/density #1000/density .5
		if (transparency > 1) { transparency = 1}
		if (transparency < .01) { transparency = .01 }
		}

		## Plot it, reordering rows so that densest points are plotted on top
		plot(x5pr,jitter(y5pr,2.5), col = alpha('black', transparency), ylim = ylimits, xaxt = 'n', xlab = genes[i], ylab = 'length', pch = 16, cex = ptsize, xlim = c(head(positions,1), tail(positions,1))) 

		points(x3pr,jitter(y3pr,2.5), col = alpha('firebrick3', transparency), xlab = genes[i], pch = 16, cex = ptsize)

		abline(v=starts[i], col = alpha('black',.4), lwd = .5)
		abline(v=starts[i]+38, col = alpha('black',.4), lwd = .5, lty=2)
		abline(v=starts[i]+22, col = alpha('black',.4), lwd = .5, lty=2)
		abline(v=starts[i]-16, col = alpha('black',.4), lwd = .5, lty=2)
		abline(v=starts[i]-38, col = alpha('black',.4), lwd = .5, lty=2)
		abline(v=starts[i]-94, col = alpha('black',.4), lwd = .5, lty=2)
		abline(v=starts[i]-200, col = alpha('black',.4), lwd = .5, lty=2)
		abline(v=starts[i]-350, col = alpha('black',.4), lwd = .5, lty=2)
		abline(v=starts[i]-950, col = alpha('black',.4), lwd = .5, lty=2)
		abline(v=ends[i], col = alpha('black',.4), lwd = .5)
		abline(h=38, lwd=.3, lty=2)

		axis(1, at=c(starts[i]-950,starts[i]-350, starts[i]-94, starts[i]-38, starts[i]-16, starts[i], starts[i]+22, starts[i]+38, ends[i]), labels=c('-950','-350','-94', '-38', '-16', 'start','+22','+38', 'stop'), cex.axis=1, las=2, tick=FALSE)

		legend('topleft', bty='n', text.col=c('black', 'red'), legend=c("5'", "3'"))

	}

	dev.off()





	# Plot pileups

	# pdf(paste0(path, libName, '_5pr3pr_widestart_lineplot.pdf'),width=3.5,height=2.2,pointsize=8)
	# 
	# x5pr_tab <- data.table(table(x5pr))
	# x3pr_tab <- data.table(table(x3pr))
	# x5prNorm_tab <- data.table(table(x5prNorm))
	# x3prNorm_tab <- data.table(table(x3prNorm))
	# 
	# x = as.numeric(x5pr_tab$x5pr)
	# y = x5pr_tab$N
	# z1 = max(y)
	# plot(x,y, col='black', type='l', xlim = c(positions[280], positions[length(positions)]-60), ylab = 'FP read counts', xlab = 'Position of mRNA from first nt in P site (nt)', xaxt = 'n', ylim = c(0, z1+500))
	# x = as.numeric(x3pr_tab$x3pr)
	# y = x3pr_tab$N
	# z2 = max(y)
	# lines(x,y, col='red', lty=1)

	# x = as.numeric(x5prNorm_tab$x5prNorm)
	# y = x5prNorm_tab$N
	# lines(x,y, col='indianred3', lty=1)
	# x = as.numeric(x3prNorm_tab$x3prNorm)
	# y = x3prNorm_tab$N
	# lines(x,y, col='indianred3', lty=2)

	# if (genes[i] == 'COX2') {
	# labels=c(-53, -16, 0)
	# labels2 = c(-48:-36, 11:31)
	# abline(v=starts[i] + labels2, lwd = 5, col = alpha('grey90', .5))

} # Finish going through each rep




# 
# 
# segments(c(starts[i] + labels[1], starts[i] + labels[2], starts[i] + labels[3]), c(0,0,0), c(starts[i] + labels[1], starts[i] + labels[2], starts[i] + labels[3]), rep(max(z1,z2), 3), col = 'black', lty = 2)
# 
# axis(1, at=seq((starts[i] - 60), (starts[i] + 60), by=20), labels=seq(-60, 60, by=20))
# text(c(starts[i] + labels[1], starts[i] + labels[2], starts[i]+ labels[3]), max(z1,z2)+200, labels=labels)
# 
# legend('topleft', bty='n', text.col=c('black', 'red'), legend=c("5'", "3'"))
# 
# 
# dev.off()
# }

# Plot counts

# Normalize counts within each library
# rph=densities/sum(densities)*100 # reads per hundred
# normrph=normdensities/sum(normdensities)*100
# 
# # Normalized by norm library
# rphNorm=rph/normrph
# rphNorms[[length(rphNorms) + 1 ]] <- rphNorm
# }
# 
# toPlot=matrix(unlist(rphNorms), nrow=sampNum, byrow=TRUE)
# 
# colors=rep(c('indianred', 'forestgreen', 'forestgreen', 'dodgerblue', 'forestgreen', 'purple', 'indianred','indianred'), each=sampNum)
# 
# pdf(paste0(path, sampName, '_', zoom, '_', u, 'upstream_',d,'downstream.pdf'),width=3,height=2.2,pointsize=8)
# 
# barplot(toPlot, beside=TRUE, names=genes,main=sampName, las=2, col=colors, ylab=paste0(normSamp, '-normalized reads near ', zoom))
# dev.off()

# graphics.off()







# source('/Users/Mary/Desktop/Data/SelectiveRP/Scripts/ScatterPlotsLengths.R')

