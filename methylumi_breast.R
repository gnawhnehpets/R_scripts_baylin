# gc()
# used (Mb) gc trigger (Mb) max used (Mb)
# Ncells 256905  6.9     407500 10.9   350000  9.4
# Vcells 222729  1.7     786432  6.0   786344  6.0

library(limma)
library(methylumi)
path <- "/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/geo/methylation/breast/"
setwd(path)
#save all idat files; char class
meta <- read.table("finalreport_cellinformation.txt", header=TRUE)
meta

#dat <- NULL
datchar <- NULL
nrow(meta)
for(i in 1:nrow(meta)){
#  dat <- rbind(dat, dir(pattern=paste(meta$barcode[i],meta$arrayposition[i],sep="_"), recursive=TRUE))
  datchar <- c(datchar, dir(pattern=paste(meta$barcode[i],meta$arrayposition[i],sep="_"), recursive=TRUE))
}
#allfiles <- dat
allfileschar <- datchar
length(datchar)
# ##############################################################################
# dim <- 1:40
# allfileschar <- datchar[dim]
# newdim <- unique(as.integer((dim+1)/2))
# metadat <- meta[newdim,]
# metadat
# gc()
# filename <- as.matrix(sub(".*/", "", allfileschar))
# dat <- methylumIDAT(barcodes = filename, idatPath="../filtered/")
# sampleNames(dat) <- paste("breast", metadat$cellline, metadat$treatment, metadat$day, sep="_")
# sampleNames(dat)
# gendat <- dat
# datname <- "1to40.txt" #1
# ##############################################################################
# dim2 <- 41:80
# allfileschar2 <- datchar[dim2]
# newdim2 <- unique(as.integer((dim2+1)/2))
# metadat2 <- meta[newdim2,]
# metadat2
# filename2 <- as.matrix(sub(".*/", "", allfileschar2))
# dat2 <- methylumIDAT(barcodes = filename2, idatPath="../filtered/")
# sampleNames(dat2) <- paste("breast", metadat2$cellline, metadat2$treatment, metadat2$day, sep="_")
# sampleNames(dat2)
# gendat <- dat2
# datname <- "41to80.txt" #2
# ##############################################################################
# dim3 <- 81:120
# allfileschar3 <- datchar[dim3]
# newdim3 <- unique(as.integer((dim3+1)/2))
# metadat3 <- meta[newdim3,]
# metadat3
# filename3 <- as.matrix(sub(".*/", "", allfileschar3))
# gc()
# dat3 <- methylumIDAT(barcodes = filename3, idatPath="../filtered/")
# sampleNames(dat3) <- paste("breast", metadat3$cellline, metadat3$treatment, metadat3$day, sep="_")
# sampleNames(dat3)
# gendat <- dat3
# datname <- "81to120.txt" #3
# ##############################################################################
# dim4 <- 121:160
# allfileschar4 <- datchar[dim4]
# newdim4 <- unique(as.integer((dim4+1)/2))
# metadat4 <- meta[newdim4,]
# metadat4
# filename4 <- as.matrix(sub(".*/", "", allfileschar4))
# dat4 <- methylumIDAT(barcodes = filename4, idatPath="../filtered/")
# sampleNames(dat4) <- paste("breast", metadat4$cellline, metadat4$treatment, metadat4$day, sep="_")
# sampleNames(dat4)
# 
# gendat <- dat4
# datname <- "121to160.txt" #4
# ##############################################################################
# dim5 <- 161:200
# allfileschar5 <- datchar[dim5]
# newdim5 <- unique(as.integer((dim5+1)/2))
# metadat5 <- meta[newdim5,]
# metadat5
# filename5 <- as.matrix(sub(".*/", "", allfileschar5))
# dat5 <- methylumIDAT(barcodes = filename5, idatPath="../filtered/")
# sampleNames(dat5) <- paste("breast", metadat5$cellline, metadat5$treatment, metadat5$day, sep="_")
# sampleNames(dat5)
# 
# gendat <- dat5
# datname <- "161to200.txt" #5
# ##############################################################################
# dim6 <- 201:240
# allfileschar6 <- datchar[dim6]
# newdim6 <- unique(as.integer((dim6+1)/2))
# metadat6 <- meta[newdim6,]
# metadat6
# filename6 <- as.matrix(sub(".*/", "", allfileschar6))
# dat6 <- methylumIDAT(barcodes = filename6, idatPath="../filtered/")
# sampleNames(dat6) <- paste("breast", metadat6$cellline, metadat6$treatment, metadat6$day, sep="_")
# sampleNames(dat6)
# gendat <- dat6
# datname <- "201to240.txt" #6
# ##############################################################################
# dim7 <- 241:280
# allfileschar7 <- datchar[dim7]
# newdim7 <- unique(as.integer((dim7+1)/2))
# metadat7 <- meta[newdim7,]
# metadat7
# filename7 <- as.matrix(sub(".*/", "", allfileschar7))
# dat7 <- methylumIDAT(barcodes = filename7, idatPath="../filtered/")
# sampleNames(dat7) <- paste("breast", metadat7$cellline, metadat7$treatment, metadat7$day, sep="_")
# sampleNames(dat7)
# gendat <- dat7
# datname <- "241to280.txt" #7
# ##############################################################################
# dim8 <- 281:320
# allfileschar8 <- datchar[dim8]
# newdim8 <- unique(as.integer((dim8+1)/2))
# metadat8 <- meta[newdim8,]
# metadat8
# filename8 <- as.matrix(sub(".*/", "", allfileschar8))
# dat8 <- methylumIDAT(barcodes = filename8, idatPath="../filtered/")
# sampleNames(dat8) <- paste("breast", metadat8$cellline, metadat8$treatment, metadat8$day, sep="_")
# sampleNames(dat8)
# gendat <- dat8
# datname <- "281to320.txt" #8
# ##############################################################################
# dim9 <- 321:360
# allfileschar9 <- datchar[dim9]
# newdim9 <- unique(as.integer((dim9+1)/2))
# metadat9 <- meta[newdim9,]
# metadat9
# filename9 <- as.matrix(sub(".*/", "", allfileschar9))
# dat9 <- methylumIDAT(barcodes = filename9, idatPath="../filtered/")
# sampleNames(dat9) <- paste("breast", metadat9$cellline, metadat9$treatment, metadat9$day, sep="_")
# sampleNames(dat9)
# gendat <- dat9
# datname <- "321to360.txt" #9
# ##############################################################################
# dim10 <- 361:400
# allfileschar10 <- datchar[dim10]
# newdim10 <- unique(as.integer((dim10+1)/2))
# metadat10 <- meta[newdim10,]
# metadat10
# filename10 <- as.matrix(sub(".*/", "", allfileschar10))
# dat10 <- methylumIDAT(barcodes = filename10, idatPath="../filtered/")
# sampleNames(dat10) <- paste("breast", metadat10$cellline, metadat10$treatment, metadat10$day, sep="_")
# sampleNames(dat10)
# gendat <- dat10
# datname <- "361to400.txt" #10
# ##############################################################################
# dim11 <- 401:420
# allfileschar11 <- datchar[dim11]
# newdim11 <- unique(as.integer((dim11+1)/2))
# metadat11 <- meta[newdim11,]
# metadat11
# filename11 <- as.matrix(sub(".*/", "", allfileschar11))
# dat11 <- methylumIDAT(barcodes = filename11, idatPath="../filtered/")
# sampleNames(dat11) <- paste("breast", metadat11$cellline, metadat11$treatment, metadat11$day, sep="_")
# sampleNames(dat11)
# gendat <- dat11
# datname <- "401to420.txt" #11


##################################################################################################################
##################################################################################################################
##################################################################################################################
tissue <- "breast"
betas <- betas(gendat)
pvals <- pvals(gendat)
methylated <- methylated(gendat)
unmethylated <- unmethylated(gendat)
colnames(betas)

dim(betas)
dim(pvals)
dim(methylated)
dim(unmethylated)

mproc <- NULL
msign <- NULL
for(i in 1:10){gc()}
for(i in 1:20){
#for(i in 1:25){
  #for `Matrix processed` worksheet
  temp <- cbind(betas[,i],pvals[,i])
  colnames(temp) <- c(colnames(gendat)[i], "Detection Pval")
  print(colnames(temp))
  mproc <- cbind(mproc, temp)
  
  #for `Matrix signal` worksheet
  temp2 <- cbind(unmethylated[,i],methylated[,i],pvals[,i])
  colnames(temp2) <- c(paste(colnames(gendat)[i], "Unmethylated Signal", sep=" "), paste(colnames(gendat)[i], "Methylated Signal", sep=" "), paste(colnames(gendat)[i], "Detection Pval", sep=" "))
  print(colnames(temp2))
  msign <- cbind(msign, temp2)
  #print(colnames(dat)[i])
}

rm(pvals)
rm(temp)
rm(temp2)
rm(methylated)
rm(unmethylated)
rm(betas)

gc()

#write.table(mproc, file = paste(tissue,"mproc", datname, sep="-"), quote = FALSE, row.names = TRUE, col.names=NA, sep="\t")
#write.table(msign, file = paste(tissue,"msign", datname, sep="-"), quote = FALSE, row.names = TRUE, col.names=NA, sep="\t")
write.table(mproc, file = paste(tissue,"mproc", datname, sep="-"), quote = FALSE, row.names = FALSE, col.names=TRUE, sep="\t")
write.table(msign, file = paste(tissue,"msign", datname, sep="-"), quote = FALSE, row.names = FALSE, col.names=TRUE, sep="\t")

# BEWARE, CLEARS ALL DATA SAVED IN RAM
gc()
rm(list=ls())
#
















# #############################################################################
# #############################################################################
# library(limma)
# library(methylumi)
# #setwd("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/geo/methylation/colon/7786923045/")
# path <- "/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/geo/methylation/breast/"
# setwd(path)
# meta <- read.table("finalreport_cellinformation.txt", header=TRUE)
# meta
# #save all idat files
# #files <- dir(pattern=".idat")
# #get all directories
# files <- dir()[file.info(dir())$isdir]
# files
# length(files)
# head(meta)
# 
# dir <- NULL
# dir
# for(i in 1:nrow(meta)){
#   dir <- c(dir, paste(meta$barcode[i], meta$arrayposition[i], sep="_"))
# }
# dir
# allfiles <- NULL
# for(i in 1:length(dir)){
#   num <- sub("_.*", "", dir[i])
#   pathtoidat <- paste(path, num, "/", sep="")
#   idatfile <- dir(path=pathtoidat, pattern=dir[i])  
#   print(idatfile)
#   allfiles <- c(allfiles, paste(path, num, "/", idatfile, sep=""))
# }
# allfiles
# length(allfiles)/2
# ?methylumIDAT
# dat <- methylumIDAT(barcodes = allfiles)
# dat
# betas <- betas(dat)
# pvals <- pvals(dat)
# methylated <- methylated(dat)
# unmethylated <- unmethylated(dat)
# dim(betas)
# dim(pvals)
# dim(methylated)
# dim(unmethylated)
