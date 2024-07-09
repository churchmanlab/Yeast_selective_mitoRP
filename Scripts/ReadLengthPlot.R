

library(data.table)
library(pheatmap)
library(inlmisc)
library(RColorBrewer)
library(stringr)
# display.brewer.all()


Folder = 'TAs'
Experiment = 'SelRP' # NoXlink SelRP Disome SelRP_2 SelRP_3 SelRP_Oxa1
Species = 'Aligned.Mito_mRNA.noDups_noSoft' # 
path = paste0('/Users/Mary/Desktop/Data/SelectiveRP/',Folder,'/MitoRP/LengthDistribution/')


# Length distributions 
DT <- data.table(read.table(paste0(path, Experiment,'_',Species, 'LengthDist.txt'), sep = '\t', header=TRUE, stringsAsFactors = FALSE, fill=TRUE))


# Replace NA values with 0
DT[is.na(DT)] <- 0

samples=colnames(DT)[2:length(DT)]

soi=samples
name='All' # Aep2 [1:12] Aep1 [13:18] Cbs1Cbs2 [19:26] Atp22Smt1Cbp1Pet309 [27:35]
n = length(soi) # number of samples, for reversing list later


# Normalize
normDT <- cbind(Length=DT$Length, DT[,..soi])
for (s in soi) {
normDT[, c(s) := DT[[s]]/sum(DT[[s]])*100]
}

# Remove very large values, which are all 0
# lower <- normDT[Length <=42]
# sampNormDT <- copy(lower)
# # Remove length column
# sampNormDT[, c('Length') := NULL]
# # Reverse order of columns
# sampList = c()
# for (i in seq(1,n)) {sampList = append(sampList,paste0('sample',i))}
# sampNames = sapply(sampList, get)
# revSampNames = rev(sampNames)
# # colors = c('dodgerblue','darkblue','tan1','tan3','tan1','tan3','tan1','tan3','tan1','tan3','tan1','tan3')
# # colors=rep(c('tan1','tan3'), times=(n/2+1))

sampNormDT <- normDT[,-1]

# setcolorder(sampNormDT, revSampNames)

# Heat Plot
pdf(paste0(path, Experiment,'_',Species, 'LengthDist_heatplot_',name,'.pdf'),
     width=4,
     height=5,
     )

   
pheatmap(t(as.matrix(sampNormDT, rownames=normDT$Length)), cluster_rows = FALSE, cluster_cols = FALSE, border_color = NA, col = GetColors(256, "YlOrBr"), na_col = "white",legend=TRUE, scale = "none", cex=.9) # cm.colors(256)  GetColors(256, "YlOrBr") brewer.pal(n=9,name='Blues')

dev.off()


# Line Plot

xlimits=c(20,50)
ylimits=c(0,21)
pdf(paste0(path, Experiment,'_',Species, 'LengthDist_lineplot_',name,'.pdf'),
     width=3,
     height=3.4,
     pointsize=10.5
     )
# colors = c(alpha('orange',.7), 'chocolate1')
colors=rep(c('orange', 'orange', 'indianred4', 'indianred4'),10)
# colors=colors()[1:n]

plot(spline(normDT$Length,normDT[[soi[1]]], method = 'n', n=150), type = 'l', lwd = 3, col = colors[1], ylim = ylimits, xlab = 'Length', ylab = '%', xlim = xlimits)
for (i in c(2:n)) {
lines(spline(normDT$Length,normDT[[soi[i]]], method = 'n', n=150), lwd = 3, col = colors[i], lty=1)
}

legend('topleft', legend=c(soi), col=c(colors), lty=1, lwd = 3, cex=0.7, bty='n') 


for (i in seq(1,length(soi),by=2)) {
plot(spline(normDT$Length,normDT[[soi[i]]], method = 'n', n=150), type = 'l', lwd = 2, col = colors[i], ylim = ylimits, xlab = 'Length', ylab = '%', xlim = xlimits)
# title(str_sub(sampNames[1],end=-3)) #, line=-2
lines(spline(normDT$Length,normDT[[soi[i+1]]], method = 'n', n=150), lwd = 2, col = colors[i+1], lty=1)
legend('topleft', legend=c(soi[i],soi[i+1]), col=c(colors[i], colors[i+1]), lty=1, cex=0.7, bty='n') 

# For plotting 3 to a plot
# lines(spline(normDT$Length,normDT[[soi[i+2]]], method = 'n', n=150), lwd = 2, col = colors[i+2], lty=1)
# legend('topleft', legend=c(soi[i],soi[i+1],soi[i+2]), col=c(colors[i], colors[i+1], colors[i+2]), lty=1, cex=0.7, bty='n') 

}

dev.off()

# source('/Users/Mary/Desktop/Data/SelectiveRP/Scripts/ReadLengthPlot.R')
