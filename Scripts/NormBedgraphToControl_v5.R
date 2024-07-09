

### Written 11/12/22 by M. Couvillion

### Normalize by taking the log2 fold-change of a sample bedGraph to a control bedGraph in a window across the mito genome
### Requires a bedGraph with every position filled in
### v3 adds a flag for positions where there isn't enough data
### v4 uses NA instead of 0 values to make it more clear where there isn't enough data
### v5 uses a sliding window instead of fixed

library(data.table)
library(stringr)
library(purrr)
library(rlist)
library('scales')
library(gtsummary)

args <- commandArgs(trailingOnly = TRUE)

normSamp <- args[1] #'Mrps17_1'
Lib = args[2]
window = as.numeric(args[3])
libsuffix = args[4]
normsuffix = args[5]
rawfilesuffix = args[6]
zerothreshold = as.numeric(args[7])


sampDT <- data.table(read.table(paste0(Lib,libsuffix, '.bedGraph'), row.names=3, skip=1))
normDT <- data.table(read.table(paste0(normSamp,normsuffix, '.bedGraph'), row.names=3, skip=1))
rawsampDT <- data.table(read.table(paste0('../',Lib,rawfilesuffix, '.bedGraph'), row.names=3, skip=1))
rawnormDT <- data.table(read.table(paste0('../',normSamp,rawfilesuffix, '.bedGraph'), row.names=3, skip=1))


DT <- cbind(sampDT, normDT[[3]], rawsampDT[[3]], rawnormDT[[3]]) # column of interest is 3 because the original third col is now the row names (see above)
names(DT) <- c('chr', 'pos0', 'samp', 'norm', 'rawsamp', 'rawnorm')

# Add column to aggregate on by window
# remainder=nrow(DT)%%window
# labvec = rep(1:((nrow(DT)+remainder)/window), each=window)
# DT[, labvec := labvec[1:nrow(DT)]]
DT <- DT[chr == 'chrM']



# Get values over sliding/rolling window
samp.Sum = frollsum(DT$samp, window, algo='exact', align='center', na.rm=TRUE)
norm.Sum = frollsum(DT$norm, window, algo='exact', align='center', na.rm=TRUE)
rawsamp.Sum = frollsum(DT$rawsamp, window, algo='exact', align='center', na.rm=TRUE)
rawnorm.Sum = frollsum(DT$rawnorm, window, algo='exact', align='center', na.rm=TRUE)
f = function(x) length(x[x==0])
sampzeroCount = frollapply(DT$norm, window, f, align='center')
normzeroCount = frollapply(DT$samp, window, f, align='center')

# Make new DT with info above
DT2=DT[,.(chr=chr, pos0=pos0, samp.Sum=samp.Sum, norm.Sum=norm.Sum, sampzeroCount=sampzeroCount,normzeroCount=normzeroCount, rawsamp.Sum.log1p=log(rawsamp.Sum+1), rawnorm.Sum.log1p=log(rawnorm.Sum+1))]


# If there is not enough coverage assign NA
maxzeros=zerothreshold # Can have up to zerothreshold zeros in window round(((window+1)*(zerothreshold/10)),0) # Can have up to threshold zero values
DT2[, samp.Norm := log(samp.Sum/norm.Sum, 2)][sampzeroCount > maxzeros, samp.Norm := NA][normzeroCount > maxzeros, samp.Norm := NA]
DT2[, samp.Norm.nolog := samp.Sum/norm.Sum][sampzeroCount > maxzeros, samp.Norm.nolog := NA][normzeroCount > maxzeros, samp.Norm.nolog := NA]
# DT2[, suff_data := 'T'][sampzeroCount > maxzeros, suff_data := 'F'][normzeroCount > maxzeros, suff_data := 'F']


DT2[DT2 == -Inf] <- NA 
DT2[DT2 == Inf] <- NA

# Make final bedGraph table
dt = data.table(chr = DT2$chr, pos0 = DT2$pos0, pos1 = DT2$pos0+1, val = DT2$samp.Norm, sampcount=DT2$rawsamp.Sum.log1p, normcount=DT2$rawnorm.Sum.log1p, val_nolog = DT2$samp.Norm.nolog)


filename = paste0(Lib,libsuffix,'_',normSamp,'_normRollWin',window,'_max',zerothreshold,'zeros.bedGraph')
cat("track type=bedGraph \n", file=filename)
write.table(dt, file=filename, sep=("\t"), quote=FALSE, row.names=FALSE, append=TRUE,col.names=FALSE)


# And make them with 0s instead of NA to load on IGV (both log-trans and no log)
DT2igv <- DT2

DT2igv$samp.Norm[is.na(DT2igv$samp.Norm)] <- 0
dt2 = data.table(chr = DT2igv$chr, pos0 = DT2igv$pos0, pos1 = DT2igv$pos0+1, val = DT2igv$samp.Norm, sampcount=DT2igv$rawsamp.Sum.log1p, normcount=DT2igv$rawnorm.Sum.log1p)

DT2igv$samp.Norm.nolog[is.na(DT2igv$samp.Norm.nolog)] <- 0
dt2nolog = data.table(chr = DT2igv$chr, pos0 = DT2igv$pos0, pos1 = DT2igv$pos0+1, val = DT2igv$samp.Norm.nolog, sampcount=DT2igv$rawsamp.Sum.log1p, normcount=DT2igv$rawnorm.Sum.log1p)

filename = paste0(Lib,libsuffix,'_',normSamp,'_normRollWin',window,'_max',zerothreshold,'zeros_forIGV.bedGraph')
cat("track type=bedGraph \n", file=filename)
write.table(dt2, file=filename, sep=("\t"), quote=FALSE, row.names=FALSE, append=TRUE,col.names=FALSE)

filename = paste0(Lib,libsuffix,'_',normSamp,'_normRollWin',window,'_max',zerothreshold,'zeros_forIGV_nolog.bedGraph')
cat("track type=bedGraph \n", file=filename)
write.table(dt2nolog, file=filename, sep=("\t"), quote=FALSE, row.names=FALSE, append=TRUE,col.names=FALSE)


