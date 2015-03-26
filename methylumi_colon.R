rm(list=ls())
library(limma)
library(methylumi)
path <- "/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/geo/methylation/colon/"
setwd(path)
#save all idat files; char class
meta <- read.table("finalreport_cellinformation.txt", header=TRUE)
meta

#dat <- NULL
datchar <- NULL
nrow(meta)
for(i in 1:nrow(meta)){
  datchar <- c(datchar, dir(pattern=paste(meta$barcode[i],meta$arrayposition[i],sep="_"), recursive=TRUE))
}
#allfiles <- dat
allfileschar <- datchar
length(datchar) # equal to max range of dim
# ##############################################################################
dim <- 1:40
allfileschar <- datchar[dim]
newdim <- unique(as.integer((dim+1)/2))
metadat <- meta[newdim,]
metadat
gc()
filename <- as.matrix(sub(".*/", "", allfileschar))
dat <- methylumIDAT(barcodes = filename, idatPath="../filteredcolon/")
sampleNames(dat) <- paste("colon", metadat$cellline, metadat$treatment, metadat$day, sep="_")
sampleNames(dat)
gendat <- dat
datname <- "1to40.txt" #1
# ##############################################################################
# dim2 <- 41:80
# allfileschar2 <- datchar[dim2]
# newdim2 <- unique(as.integer((dim2+1)/2))
# metadat2 <- meta[newdim2,]
# metadat2
# filename2 <- as.matrix(sub(".*/", "", allfileschar2))
# dat2 <- methylumIDAT(barcodes = filename2, idatPath="../filteredcolon/")
# sampleNames(dat2) <- paste("colon", metadat2$cellline, metadat2$treatment, metadat2$day, sep="_")
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
# dat3 <- methylumIDAT(barcodes = filename3, idatPath="../filteredcolon/")
# sampleNames(dat3) <- paste("colon", metadat3$cellline, metadat3$treatment, metadat3$day, sep="_")
# sampleNames(dat3)
# gendat <- dat3
# datname <- "81to120.txt" #3
# ##############################################################################
# dim4 <- 121:152
# allfileschar4 <- datchar[dim4]
# newdim4 <- unique(as.integer((dim4+1)/2))
# metadat4 <- meta[newdim4,]
# metadat4
# filename4 <- as.matrix(sub(".*/", "", allfileschar4))
# dat4 <- methylumIDAT(barcodes = filename4, idatPath="../filteredcolon/")
# sampleNames(dat4) <- paste("colon", metadat4$cellline, metadat4$treatment, metadat4$day, sep="_")
# sampleNames(dat4)
# gendat <- dat4
# datname <- "121to152.txt" #3

##################################################################################################################
##################################################################################################################
##################################################################################################################
tissue <- "colon"
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

# write.table(mproc, file = paste(tissue,"mproc", datname, sep="-"), quote = FALSE, row.names = TRUE, col.names=NA, sep="\t")
# write.table(msign, file = paste(tissue,"msign", datname, sep="-"), quote = FALSE, row.names = TRUE, col.names=NA, sep="\t")
write.table(mproc, file = paste(tissue,"mproc", datname, sep="-"), quote = FALSE, row.names = FALSE, col.names=TRUE, sep="\t")
write.table(msign, file = paste(tissue,"msign", datname, sep="-"), quote = FALSE, row.names = FALSE, col.names=TRUE, sep="\t")

# BEWARE, CLEARS ALL DATA SAVED IN RAM
for(i in 1:10){gc()}
rm(list=ls())
