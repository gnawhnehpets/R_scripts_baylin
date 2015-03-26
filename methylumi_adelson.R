library(limma)
library(methylumi)
path <- "/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/geo/methylation/adelson/"
setwd(path)
#save all idat files; char class
meta <- read.table("finalreport_cellinformation.txt", header=TRUE)
meta

#dat <- NULL
datchar <- NULL

# cellline treatment   day    barcode arrayposition
# 1      A2780        d3  Mock 7786915143        R01C01
# 2      CAOV3        d3 Mock2 7786915168        R03C02
# 3      EF027        d3  mock 8622007232        R01C01
# 4        ES2        d3  mock 8622007206        R01C01
# 5        Hey        d3  mock 8622007185        R03C02
# 6  Kuramochi        d3  mock 8622007232        R05C01
# 7      OAW28        d3  mock 8622007232        R03C02
# 8     OV2008        d3  mock 8622007206        R05C01
# 9     OVCAR3        d3  Mock 7786915168        R05C01
# 10    OVCAR5        d3  Mock 8622007206        R03C02
# 11    OVKATE        d3  mock 8622007185        R01C01
# 12     SKOV3        d3  mock 8795207151        R01C01
# 13     Tyknu        d3  mock 8622007185        R05C01


nrow(meta)
for(i in 1:nrow(meta)){
  datchar <- c(datchar, dir(path= "../filteredovarian", pattern=paste(meta$barcode[i],meta$arrayposition[i],sep="_"), recursive=TRUE))
}
#allfiles <- dat
allfileschar <- datchar
length(datchar)
# ##############################################################################
dim <- 1:34
allfileschar <- datchar[dim]
newdim <- unique(as.integer((dim+1)/2))
metadat <- meta[newdim,]
metadat
gc()
filename <- as.matrix(sub(".*/", "", allfileschar))
dat <- methylumIDAT(barcodes = filename, idatPath="../filteredovarian/")
sampleNames(dat) <- paste("ovarian", metadat$cellline, metadat$treatment, metadat$day, sep="_")
sampleNames(dat)
gendat <- dat
datname <- "1to34.txt" #1


##################################################################################################################
##################################################################################################################
##################################################################################################################
tissue <- "ovarian"
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

write.table(mproc, file = paste(tissue,"mproc", datname, sep="-"), quote = FALSE, row.names = TRUE, col.names=NA, sep="\t")
write.table(msign, file = paste(tissue,"msign", datname, sep="-"), quote = FALSE, row.names = TRUE, col.names=NA, sep="\t")
# write.table(mproc, file = paste(tissue,"mproc", datname, sep="-"), quote = FALSE, row.names = FALSE, col.names=TRUE, sep="\t")
# write.table(msign, file = paste(tissue,"msign", datname, sep="-"), quote = FALSE, row.names = FALSE, col.names=TRUE, sep="\t")

# BEWARE, CLEARS ALL DATA SAVED IN RAM
gc()
rm(list=ls())