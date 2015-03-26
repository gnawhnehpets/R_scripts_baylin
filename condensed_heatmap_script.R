rm(list=ls())
for(i in 1:20){gc()}
library(limma)
library(FactoMineR)
# setwd("/home/steve/Desktop/analysis/cohortoneandtwo//")
# target <- readTargets("/home/steve/Desktop/analysis/cohortoneandtwo//targets_full.txt")
# RG <- read.maimages(target, source="agilent", path="/home/steve/Desktop/analysis/cohortoneandtwo/")
setwd("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/m084b trial/")
target <- readTargets("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/cohortoneandtwo/targets_full.txt")
RG <- read.maimages(target, source="agilent", path="/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/cohortoneandtwo/")
# setwd("Y:/users/shwang26/m084b trial/")
# target <- readTargets("Y:/users/shwang26/cohortoneandtwo/targets_full.txt")
# RG <- read.maimages(target, source="agilent", path="Y:/users/shwang26/cohortoneandtwo/")

meta <- read.table("patientstatus_prepost_allsamples_jh.csv", header=TRUE)
meta
meta$sixmonthsoverall

meta[1:34,1]
meta[,1]
# setwd("Y:/users/shwang26/m084b trial/")
# target <- readTargets("Y:/users/shwang26/cohortoneandtwo/targets_full.txt")
# RG <- read.maimages(target, source="agilent", path="Y:/users/shwang26/cohortoneandtwo/")
cohorttwo<- c(target$Cy3_Sample, target$Cy5_sample)
cohorttwo
RG$weights <-  matrix(rep(RG$genes$ControlType,ncol(RG$R)),ncol=ncol(RG$R),byrow=F)
RG$weights[RG$genes$ControlType!=0,] <- 0
RG$weights[RG$genes$ControlType==0,] <- 1

as.data.frame(cbind(as.matrix(paste(target$Cy3, target$Cy3_Sample, sep="-")),as.matrix(paste(target$Cy5, target$Cy5_sample, sep="-"))))
as.data.frame(sub(".*_1_", "", RG$targets$FileName))

#remove outlier samples PH1824/1641
dat <- RG[,-18]
dat <- dat[,-6]
targets <- target[-18,]
targets <- targets[-6,]

as.matrix(sub(".*_1_", "", dat$targets$FileName)) 

samp <- data.frame(paste(targets$Cy3,targets$Cy3_Sample,sep="-"), paste(targets$Cy5, targets$Cy5_sample, sep="-"))
colnames(samp) <- c("cy3", "cy5")
samp
          cy3         cy5
#R1  post-JHH004  pre-JHH004
#M2   pre-PH1550  pre-PH1600
#R3   pre-JHH005 post-JHH005
#R4   pre-PH1604  pre-PH1606
#R5   pre-PH1622  pre-PH1631
6   pre-PH1623  pre-PH1632
#R7  post-PH1815  pre-PH1815
#R8   pre-PH1635  pre-PH1644
9   pre-PH1827 post-PH1827
#M10  pre-PH1811  pre-PH1816
#R11  pre-PH1636 post-PH1640
#R12  pre-PH1892  pre-PH1902
REDO 13 post-PH1612  pre-PH1612
#R14  pre-PH1868  pre-PH1887
#R15  pre-PH1616 post-PH1616
16 post-PH1910  pre-PH1910
#R17 post-PH1843  pre-PH1843
#R18  pre-PH1913  pre-PH1886
#R19 post-PH1844  pre-PH1844
REDO 20  pre-PH1861 post-PH1861
#R21 post-PH1869  pre-PH1869
#R22  pre-PH1871 post-PH1871
#R23 post-PH1876  pre-PH1876
24  pre-PH1641 post-PH1641
#O25 post-PH1824  pre-PH1824
x 26  pre-PH1545 post-PH1545
x 27 post-PH1565  pre-PH1565
x 28 post-PH1910  pre-PH1910
#O29  pre-PH1641 post-PH1640
x 30 post-PH1612  pre-PH1824
x 31  pre-PH1861 post-PH1900
> 

#1 - no changes
# ALL SAMPLEs
# pos <- c(1,2,7,8,11,14:16,18:22,24,25,3:6,9,10,12,13,17,23,26)
# PAIRS ONLY 
# pos <- c(1,2,7,8,11,14:16,18:22,24,25,10)
#2
#pos <- c(1,2,7,8,11,14:16,18:22,24,25)
#3
# pos <- c(8,14,16,19,20,22,24,25)
#4
#pos <- c(8,14,20,24,25)
#5
#pos <- c(8,14)
# LESLIE
pos <- c(11,14,17,20,25,26)

dat <- dat[,pos]
colnames(dat)
dim(dat)
targets <- targets[pos,]
dim(targets)
targets

as.data.frame(cbind(as.matrix(paste(targets$Cy3, targets$Cy3_Sample, sep="-")),as.matrix(paste(targets$Cy5, targets$Cy5_sample, sep="-"))))

names(dat)
# QA #Note: DO NOT REMOVE CONTROLS 
spottypes <- readSpotTypes()
dat$genes$Status <- controlStatus(spottypes, dat)
limma::plotMA3by2(dat, prefix="prenorm_MA_lesli")
colnames(dat$G)
names(dat)

boxplot(data.frame(log2(dat$Gb)),main="Green background", names = NULL)
dev.off()
boxplot(data.frame(log2(dat$Rb)),main="Red background", names = NULL)
dev.off()
boxplot(data.frame(log2(dat$G)),main="Green foreground", names = NULL)
dev.off()
boxplot(data.frame(log2(dat$R)),main="Red foreground", names = NULL)
dev.off()


# BG SUBTRACTION
# dat <- backgroundCorrect(dat, method="subtract")

#1,2,3,4
# normwithin <-normalizeWithinArrays(dat,method='loess',bc.method='normexp', offset=50)
# normwithin <-normalizeWithinArrays(dat,method='loess',bc.method='none')
# LESLIE
normwithin <-normalizeWithinArrays(dat,method='none',bc.method='none')
normbetween <-normalizeBetweenArrays(normwithin,method='Aquantile')

png(filename = "normalized - plotDensities(normwithin).png", width = 800, height = 800)  
plotDensities(normwithin, main = "RG densities - normwithin")
dev.off()
png(filename = "boxplot(normwithin$M).png", width = 1800, height = 1000)  
boxplot(normwithin$M, names=NULL)
dev.off()
# normbetween <-normalizeBetweenArrays(normwithin,method='Aquantile')
png(filename = "normalized - plotDensities(normbetween).png", width = 800, height = 800)  
plotDensities(normbetween, main = "RG densities - normbetween")
dev.off()

#Remove controls from normwithin/between
normwithin <- normwithin[normbetween$genes$ControlType==0,]
normbetween <- normbetween[normbetween$genes$ControlType==0,]

# Convert MA back to RG
RGb <- RG.MA(normbetween)

spottypes <- readSpotTypes()
RGb$genes$Status <- controlStatus(spottypes, RGb)
limma::plotMA3by2(RGb, prefix="MA_leslie")
head(RGb$G)
names(RGb)
png(filename = "green-postnorm", width = 1800, height = 1000)  
boxplot(data.frame(log2(RGb$G)),main="Green post normalization", names = NULL)
dev.off()
png(filename = "red-postnorm", width = 1800, height = 1000)  
boxplot(data.frame(log2(RGb$R)),main="Red post normalization", names = NULL)
dev.off()
boxplot(data.frame(log2(dat$Rb)),main="Red background", names = NULL)
boxplot(data.frame(log2(dat$Gb)),main="Green background", names = NULL)
boxplot(data.frame(log2(dat$Rb)),main="Red background", names = NULL)
boxplot(data.frame(log2(dat$G)),main="Green foreground", names = NULL)
boxplot(data.frame(log2(dat$R)),main="Red foreground", names = NULL)

# plotDensities(RGb)
# names(RGb)
# names(dat)
# # pre-normalization
# boxplot(data.frame(log2(dat$Gb)),main="Green background - pre-normalization", names=paste(targets$Cy3, targets$Cy3_Sample, sep="-"), las=2)
# boxplot(data.frame(log2(dat$Rb)),main="Red background - pre-normalization", names=paste(targets$Cy5, targets$Cy5_sample, sep="-"), las=2)
# 
# # post-normalization
# boxplot(data.frame(log2(RGb$G)),main="Green background - normalized", names=paste(targets$Cy3, targets$Cy3_Sample, sep="-"), las=2)
# boxplot(data.frame(log2(RGb$R)),main="Red background - normalized", names=paste(targets$Cy5, targets$Cy5_sample, sep="-"), las=2)

# cy3
cy3 = RGb$R
rownames(cy3) <- RGb$genes$GeneName
colnames(cy3) <- paste(targets$Cy3, targets$Cy3_Sample, sep="-")
colnames(cy3)
# cy5
cy5 = RGb$G
rownames(cy5) <- RGb$genes$GeneName
colnames(cy5) <- paste(targets$Cy5, targets$Cy5_sample, sep="-")
colnames(cy5)

npre <- colnames(cy3[,sub("-.*","",colnames(cy3))=="pre"])
npost <- colnames(cy3[,sub("-.*","",colnames(cy3))=="post"])
npre2 <- colnames(cy5[,sub("-.*","",colnames(cy5))=="pre"])
npost2 <- colnames(cy5[,sub("-.*","",colnames(cy5))=="post"])
c(npre, npre2)

c(npost2, npost)
ppre <- cy3[,sub("-.*","",colnames(cy3))=="pre"]
ppost <- cy3[,sub("-.*","",colnames(cy3))=="post"]
ppre2 <- cy5[,sub("-.*","",colnames(cy5))=="pre"]
ppost2 <- cy5[,sub("-.*","",colnames(cy5))=="post"]
pppre <- cbind(ppre, ppre2)
pppost <- cbind(ppost2, ppost)
colnames(pppost) == meta$patient
as.data.frame(colnames(pppre),colnames(pppost))
######################################################################################################################################
#r <- cbind(cy3, cy5)
# all pre/post (pairs)
allsamples <- cbind(pppre,pppost)
# just pre
allsamples.preonly <- pppre

# pre/post for pairs only
pairs <- cbind(pppre, pppost)
# pres for pairs only
pairs.preonly <- pppre

##########################
# DIAGNOSTIC
diag <- log2(allsamples)
datpca <- prcomp(t(diag), scale=TRUE)
png(filename = "pca_allsamples_scale.png", width = 800, height = 800)  
title <- "PCA of colon trial data\nall samples, nl/scale=T"
datpca <- prcomp(t(diag), scale=FALSE)
png(filename = "pca_allsamples_no-scale.png", width = 800, height = 800)  
title <- "PCA of colon trial data\nall samples, nl/scale=F"

plot(
  range(datpca$x[, 1]),range(datpca$x[, 2]), 
  type = "n", xlab = "pre", ylab = "post", 
  main = title
)
dimnames(datpca)
points(datpca$x[, 1][targets$Cy3 == "pre"], datpca$x[, 2][targets$Cy3 == "pre"], bg = "Red", pch = 21)
text(datpca$x[, 1][targets$Cy3 == "pre"], datpca$x[, 2][targets$Cy3 == "pre"], label = "pre", col = "red", cex = .8)
text(datpca$x[, 1][targets$Cy3 == "pre"], datpca$x[, 2][targets$Cy3 == "pre"], label = names(datpca$x[,1][targets$Cy3 == "pre"]), col = "red", cex = .8)
points(datpca$x[, 1][targets$Cy5 == "post"], datpca$x[, 2][targets$Cy5 == "post"], bg = "Blue", pch=21)
text(datpca$x[, 1][targets$Cy5 == "post"], datpca$x[, 2][targets$Cy5 == "post"], label = names(datpca$x[,1][targets$Cy3 == "post"]), col = "blue", cex = .8)
text(datpca$x[, 1][targets$Cy5 == "post"], datpca$x[, 2][targets$Cy5 == "post"], label = "post", col = "blue", cex = .8)
leg.names <- c("pre", "post")       
legend("topleft",leg.names,col=c("red","blue"),pch=15,cex=.7,horiz=F)

dev.off()
##########################


save(allsamples, file="allsamples")
save(allsamples.preonly, file="allsamples.preonly")
save(pairs, file="pairs")
save(pairs.preonly, file="pairs.preonly")

load("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/m084b trial/R objects/allsamples")
load("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/m084b trial/R objects/allsamples.preonly")
load("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/m084b trial/R objects/pairs")
load("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/m084b trial/R objects/pairs.preonly")

r <- allsamples
r <- allsamples.preonly
r <- pairs
r <- pairs.preonly

####################################################################################################################################################################
# all samples, lt/mt 2 months therapy
r <- allsamples
r <- allsamples.preonly

# NOTE: NEED TO INCLUDE ALL SAMPLES, NOT JUST PAIRS
missing <- "pre-PH1816"
x <- which(colnames(r)==missing, arr.ind=TRUE)
colnames(x)
m <- r[, -c(x)]

# # NOTE: for pairs only 
# r <- pairs
# r <- pairs.preonly
# m <- r

dim(r)
colnames(r)
dim(m)
colnames(m)
lt2mo = NULL
mt2mo = NULL
lt = NULL
mt = NULL

for(i in 1:length(colnames(m))){
  cname <- NULL
  cname <- colnames(m)[i]
  b <- NULL
  b <- cname == meta[,1]
  cname
  b
  if(meta$twomonthstherapy[b]==1){
    w <- NULL
    w <- which(colnames(m)==cname, arr.ind=TRUE)
    print(paste("MT",w,sep=":"))
    mt2mo = c(mt2mo, colnames(m)[i])
    mt <- c(mt, w)
  }else if(meta$twomonthstherapy[b]==0){
    v <- NULL
    v <- which(colnames(m)==cname, arr.ind=TRUE)
    print(paste("LT",v,sep=":"))
    lt2mo = c(lt2mo, colnames(m)[i])
    lt <- c(lt, v)
  }
}

more <- m[,mt]
less <- m[,lt]
dim(more)
dim(less)

pm <- paste("MORE-", colnames(m)[mt], sep="")
pl <- paste("LESS-", colnames(m)[lt], sep="")
colnames(more) <- pm
colnames(less) <- pl

# all samples
mo <- more[,c(1,10,2,11,3,12,5,14,6,15,4,7:9,13)]
le <- less[,c(1,26,2,27,3,28,4,29,14,30,15,31,16,32,17,33,18,34,19,35,23,36,5,6,7,8,9,10,11,12,13,20,21,22,24,25)]
# pair
mo <- more[,c(1,6,2,7,3,8,4,9,5,10)]
le <- less[,c(1,12,2,13,3,14,4,15,5,16,6,17,7,18,8,19,9,20,10,21,11,22)]
# pres only
# mo <- more[,c(1,10,2,11,3,12,5,14,6,15,4,7:9,13)]
# le <- less[,c(1,26,2,27,3,28,4,29,14,30,15,31,16,32,17,33,18,34,19,35,23,36,5,6,7,8,9,10,11,12,13,20,21,22,24,25)]
mo <- more
le <- less
data.frame(colnames(mo))
data.frame(colnames(le))
# mo <- mo[,c(1,3,5,7,9,11:14)]
# le <- le[,c(1,3,5,7,9,11,13,15,17,18,19,21,23:36)]


dim(mo)
dim(le)
r <- cbind(le, mo)
colnames(r)
####################################################################################################################################################################
####################################################################################################################################################################
# all samples, lt/mt 12 month overall survival
r <- allsamples
r <- allsamples.preonly

# NOTE: NEED TO INCLUDE ALL SAMPLES, NOT JUST PAIRS
missing <- "pre-PH1816"
x <- which(colnames(r)==missing, arr.ind=TRUE)
colnames(x)
m <- r[, -c(x)]

# # NOTE: for pairs only 
# r <- pairs
# r <- pairs.preonly
# m <- r

dim(r)
colnames(r)
dim(m)
colnames(m)
lt12mo = NULL
mt12mo = NULL
lt2 = NULL
mt2 = NULL
colnames(m)
for(i in 1:length(colnames(m))){
  cname <- NULL
  cname <- colnames(m)[i]
  b <- NULL
  b <- cname == meta[,1]
  cname
  b
  if(meta$twelvemonthsoverall[b]==1){
    w <- NULL
    w <- which(colnames(m)==cname, arr.ind=TRUE)
    print(paste("MT",w,sep=":"))
    mt12mo = c(mt12mo, colnames(m)[i])
    mt2 <- c(mt2, w)
  }else if(meta$twelvemonthsoverall[b]==0){
    v <- NULL
    v <- which(colnames(m)==cname, arr.ind=TRUE)
    print(paste("LT",v,sep=":"))
    lt12mo = c(lt12mo, colnames(m)[i])
    lt2 <- c(lt2, v)
  }
}

more2 <- m[,mt2]
less2 <- m[,lt2]
dim(more2)
dim(less2)

pm2 <- paste("MORE", colnames(m)[mt2], sep="-")
pl2 <- paste("LESS", colnames(m)[lt2], sep="-")
colnames(more2) <- pm2
colnames(less2) <- pl2
colnames(more2)

# all samples
mo2 <- more2[,c(1,9,2,10,4,12,5,13,6,14,3,7,8,11)]
le2 <- less2[,(c(1,27,2,28,3,29,4,30,5,31,15,32,16,33,17,34,18,35,19,36,23,37,6,7,8,9,10,11,12,13,14,20,21,22,24,25,26))]

# # pair
# mo2 <- more2[,c(1,6,2,7,3,8,4,9,5,10)]
# le2 <- less2[,c(1,12,2,13,3,14,4,15,5,16,6,17,7,18,8,19,9,20,10,21,11,22)]

# pres only
mo2 <- more2
le2 <- less2

dim(mo2)
dim(le2)
r <- cbind(le2, mo2)
colnames(r)
####################################################################################################################################################################
####################################################################################################################################################################
####################################################################################################################################################################
# all samples, lt/mt 9 month overall survival
r <- allsamples
r <- allsamples.preonly

# NOTE: NEED TO INCLUDE ALL SAMPLES, NOT JUST PAIRS
missing <- "pre-PH1816"
x <- which(colnames(r)==missing, arr.ind=TRUE)
colnames(x)
m <- r[, -c(x)]

# # NOTE: for pairs only 
# r <- pairs
# r <- pairs.preonly
# m <- r

dim(r)
colnames(r)
dim(m)
colnames(m)
lt9mo = NULL
mt9mo = NULL
lt3 = NULL
mt3 = NULL
meta$ninemonthsoverall[1]

for(i in 1:length(colnames(m))){
  cname <- NULL
  cname <- colnames(m)[i]
  b <- NULL
  b <- cname == meta[,1]
  cname
  b
  if(meta$ninemonthsoverall[b]==1){
    w <- NULL
    w <- which(colnames(m)==cname, arr.ind=TRUE)
    print(paste("MT",w,sep=":"))
    mt9mo = c(mt9mo, colnames(m)[i])
    mt3 <- c(mt3, w)
  }else if(meta$ninemonthsoverall[b]==0){
    v <- NULL
    v <- which(colnames(m)==cname, arr.ind=TRUE)
    print(paste("LT",v,sep=":"))
    lt9mo = c(lt9mo, colnames(m)[i])
    lt3 <- c(lt3, v)
  }
}

more3 <- m[,mt3]
less3 <- m[,lt3]
dim(more3)
dim(less3)

pm3 <- paste("MORE", colnames(m)[mt3], sep="-")
pl3 <- paste("LESS", colnames(m)[lt3], sep="-")
colnames(more3) <- pm3
colnames(less3) <- pl3
colnames(less3)
colnames(more3)

# all samples
mo3 <- more3[,c(1,11,2,12,3,13,6,15,7,16,8,17,4,5,9,10,14)]
le3 <- less3[,(c(1,25,2,26,3,27,4,28,13,29,14,30,15,31,16,32,17,33,21,34,5:12,18,19,20,22,23,24))]

# pair
mo3 <- more3[,c(1,7,2,8,3,9,4,10,5,11,6,12)]
le3 <- less3[,c(1,11,2,12,3,13,4,14,5,15,6,16,7,17,8,18,9,19,10,20)]

# pres only
mo3 <- more3
le3 <- less3

dim(mo3)
dim(le3)
r <- cbind(le3, mo3)
colnames(r)
####################################################################################################################################################################
# all samples, lt/mt 6 month overall survival
r <- allsamples
r <- allsamples.preonly

# NOTE: NEED TO INCLUDE ALL SAMPLES, NOT JUST PAIRS
missing <- "pre-PH1816"
x <- which(colnames(r)==missing, arr.ind=TRUE)
colnames(x)
m <- r[, -c(x)]

# # NOTE: for pairs only 
# r <- pairs
# r <- pairs.preonly
# m <- r

dim(r)
colnames(r)
dim(m)
colnames(m)
lt6mo = NULL
mt6mo = NULL
lt4 = NULL
mt4 = NULL
meta$sixmonthsoverall[1]

for(i in 1:length(colnames(m))){
  cname <- NULL
  cname <- colnames(m)[i]
  b <- NULL
  b <- cname == meta[,1]
  cname
  b
  if(meta$sixmonthsoverall[b]==1){
    w <- NULL
    w <- which(colnames(m)==cname, arr.ind=TRUE)
    print(paste("MT",w,sep=":"))
    mt6mo = c(mt6mo, colnames(m)[i])
    mt4 <- c(mt4, w)
  }else if(meta$sixmonthsoverall[b]==0){
    v <- NULL
    v <- which(colnames(m)==cname, arr.ind=TRUE)
    print(paste("LT",v,sep=":"))
    lt6mo = c(lt6mo, colnames(m)[i])
    lt4 <- c(lt4, v)
  }
}

more4 <- m[,mt4]
less4 <- m[,lt4]
dim(more4)
dim(less4)

pm4 <- paste("MORE", colnames(m)[mt4], sep="-")
pl4 <- paste("LESS", colnames(m)[lt4], sep="-")
colnames(more4) <- pm4
colnames(less4) <- pl4
colnames(less4)
colnames(more4)

# all samples
mo4 <- more4[,c(1,19,2,20,3,21,4,22,5,23,6,24,10,26,11,27,12,28,13,29,14,30,15,31,7,8,9,16,17,18,25)]
le4 <- less4[,c(1,17,9,18,10,19,14,20,2,3,4,5,6,7,8,11,12,13,15,16)]
# pairs only
mo4 <- more4[,c(1,13,2,14,3,15,4,16,5,17,6,18,7,19,8,20,9,21,10,22,11,23,12,24)]
le4 <- less4[,c(1,5,2,6,3,7,4,8)]

#pres only
mo4 <- more4
le4 <- less4

colnames(le4)
dim(mo4)
dim(le4)
r <- cbind(le4, mo4)
colnames(r)
dim(r)
####################################################################################################################################################################

r <- allsamples
r <- allsamples.preonly
r <- pairs
r <- pairs.preonly

new.dat <- apply(r,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
new.dat2mo <- apply(r,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
new.dat6mo <- apply(r,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
new.dat9mo <- apply(r,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
new.dat12mo <- apply(r,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 

new.datpair <- apply(r,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
new.dat2mopair <- apply(r,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
new.dat6mopair <- apply(r,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
new.dat9mopair <- apply(r,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
new.dat12mopair <- apply(r,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 

new.datpre <- apply(r,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
new.dat2mopre <- apply(r,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
new.dat6mopre <- apply(r,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
new.dat9mopre <- apply(r,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
new.dat12mopre <- apply(r,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 

new.datpairpre <- apply(r,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
# new.dat2mopairpre <- apply(r,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
# new.dat6mopairpre <- apply(r,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
# new.dat9mopairpre <- apply(r,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
# new.dat12mopairpre <- apply(r,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 

save(new.dat, file="new.dat")
save(new.dat2mo, file="new.dat2mo")
save(new.dat6mo, file="new.dat6mo")
save(new.dat9mo, file="new.dat9mo")
save(new.dat12mo, file="new.dat12mo")
save(new.datpair, file="new.datpair")
save(new.dat2mopair, file="new.dat2mopair")
save(new.dat6mopair, file="new.dat6mopair")
save(new.dat9mopair, file="new.dat9mopair")
save(new.dat12mopair, file="new.dat12mopair")
save(new.datpre, file="new.datpre")
save(new.dat2mopre, file="new.dat2mopre")
save(new.dat6mopre, file="new.dat6mopre")
save(new.dat9mopre, file="new.dat9mopre")
save(new.dat12mopre, file="new.dat12mopre")


load("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/m084b trial/R objects/new.dat")
load("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/m084b trial/R objects/new.dat2mo") 
load("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/m084b trial/R objects/new.dat6mo") 
load("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/m084b trial/R objects/new.dat9mo") 
load("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/m084b trial/R objects/new.dat12mo") 
load("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/m084b trial/R objects/new.datpair")
load("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/m084b trial/R objects/new.dat2mopair") 
load("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/m084b trial/R objects/new.dat6mopair") 
load("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/m084b trial/R objects/new.dat9mopair") 
load("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/m084b trial/R objects/new.dat12mopair") 
load("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/m084b trial/R objects/new.datpre")
load("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/m084b trial/R objects/new.dat2mopre") 
load("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/m084b trial/R objects/new.dat6mopre") 
load("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/m084b trial/R objects/new.dat9mopre") 
load("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/m084b trial/R objects/new.dat12mopre") 

####################################################################################################################################################################
####################################################################################################################################################################
####################################################################################################################################################################

##########################################
# AIM ####################################
##########################################
new.datfil <- new.dat
new.datfil <- r
dim(new.datfil)
new.dat["XCL2",]
r["KLRC2",]
one<-new.datfil["ADM",]
two<-new.datfil["ADORA2A",]
three<-new.datfil["ADORA2B",]
four<-new.datfil["ALOX5AP",]
five<-new.datfil["ANKRD1",]
six<-new.datfil["ANXA1",]
seven<-new.datfil["AOC3",]
eight<-new.datfil["AOX1",]
nine<-new.datfil["B2M",]
onezero<-new.datfil["BCL2",]
oneone<-new.datfil["C4BPB",]
onetwo<-new.datfil["CAMP",]
onethree<-new.datfil["CCL2",]
onefour<-new.datfil["CCL20",]
onefive<-new.datfil["CCL26",]
onesix<-new.datfil["CCL28",]
oneseven<-new.datfil["CCL3",]
oneeight<-new.datfil["CCL3L3",]
onenine<-new.datfil["CCL4",]
twozero<-new.datfil["CCL5",]
twoone<-new.datfil["CCR9",]
twotwo<-new.datfil["CD44",]
twothree<-new.datfil["CD81",]
twofour<-new.datfil["CRP",]
twofive<-new.datfil["CSF2",]
twosix<-new.datfil["CTGF",]
twoseven<-new.datfil["CTSS",]
twoeight<-new.datfil["CXCL1",]
twonine<-new.datfil["CXCL11",]
threezero<-new.datfil["CXCL12",]
threeone<-new.datfil["CXCL2",]
threetwo<-new.datfil["CXCL3",]
threethree<-new.datfil["CXCL6",]
threefour<-new.datfil["CXCL9",]
threefive<-new.datfil["CXCR4",]
threesix<-new.datfil["CXCR7",]
threeseven<-new.datfil["DCBLD2",]
threeeight<-new.datfil["DDX58",]
threenine<-new.datfil["DEFB103A",]
fourzero<-new.datfil["EGR1",]
fourone<-new.datfil["EIF2AK2",]
fourtwo<-new.datfil["EIF4E",]
fourthree<-new.datfil["EIF4G1",]
fourfour<-new.datfil["EREG",]
fourfive<-new.datfil["FOS",]
foursix<-new.datfil["FOSL1",]
fourseven<-new.datfil["GBP1",]
foureight<-new.datfil["GBP5",]
fournine<-new.datfil["GP9",]
fivezero<-new.datfil["HCP5",]
fiveone<-new.datfil["HERC5",]
fivetwo<-new.datfil["HLA-B",]
fivethree<-new.datfil["HLA-C",]
fivefour<-new.datfil["HLA-DMB",]
fivefive<-new.datfil["HLA-DRB1",]
fivesix<-new.datfil["HSP90AA1",]
fiveseven<-new.datfil["ICAM1",]
fiveeight<-new.datfil["IFI27",]
fivenine<-new.datfil["IFI6",]
sixzero<-new.datfil["IFIH1",]
sixone<-new.datfil["IFIT1",]
sixtwo<-new.datfil["IFIT2",]
sixthree<-new.datfil["IFIT3",]
sixfour<-new.datfil["IFITM1",]
sixfive<-new.datfil["IFNGR1",]
sixsix<-new.datfil["IL1A",]
sixseven<-new.datfil["IL1B",]
sixeight<-new.datfil["IL1R2",]
sixnine<-new.datfil["IL32",]
sevenzero<-new.datfil["IL6",]
sevenone<-new.datfil["IL6ST",]
seventwo<-new.datfil["IL8",]
seventhree<-new.datfil["INHBA",]
sevenfour<-new.datfil["IRF6",]
sevenfive<-new.datfil["IRF7",]
sevensix<-new.datfil["IRF8",]
sevenseven<-new.datfil["IRF9",]
seveneight<-new.datfil["ISG15",]
sevennine<-new.datfil["ISG20",]
eightzero<-new.datfil["ITGAV",]
eightone<-new.datfil["JAK2",]
eighttwo<-new.datfil["KCNN4",]
eightthree<-new.datfil["KLK8",]
eightfour<-new.datfil["KLRC2",]
eightfive<-new.datfil["LBP",]
eightsix<-new.datfil["LCK",]
eightseven<-new.datfil["LSP1",]
eighteight<-new.datfil["LY96",]
eightnine<-new.datfil["LYN",]
ninezero<-new.datfil["LYST",]
nineone<-new.datfil["MDK",]
ninetwo<-new.datfil["MRC2",]
ninethree<-new.datfil["MT2A",]
ninefour<-new.datfil["MX1",]
ninefive<-new.datfil["NCF1",]
ninesix<-new.datfil["NCF2",]
nineseven<-new.datfil["NLRP3",]
nineeight<-new.datfil["NLRX1",]
ninenine<-new.datfil["NOD2",]
onezerozero<-new.datfil["NUP35",]
onezeroone<-new.datfil["OAS1",]
onezerotwo<-new.datfil["OAS2",]
onezerothree<-new.datfil["OAS3",]
onezerofour<-new.datfil["OASL",]
onezerofive<-new.datfil["ORM1",]
onezerosix<-new.datfil["ORM2",]
onezeroseven<-new.datfil["PAGE1",]
onezeroeight<-new.datfil["PIK3R2",]
onezeronine<-new.datfil["PLA2G7",]
oneonezero<-new.datfil["PLAT",]
oneoneone<-new.datfil["PPBP",]
oneonetwo<-new.datfil["PROS1",]
oneonethree<-new.datfil["PSMA3",]
oneonefour<-new.datfil["PSMB8",]
oneonefive<-new.datfil["PSMB9",]
oneonesix<-new.datfil["PSMC6",]
oneoneseven<-new.datfil["PTX3",]
oneoneeight<-new.datfil["RPL26",]
oneonenine<-new.datfil["RPL38",]
onetwozero<-new.datfil["RSAD2",]
onetwoone<-new.datfil["S100A7",]
onetwotwo<-new.datfil["S100A8",]
onetwothree<-new.datfil["S100A9",]
onetwofour<-new.datfil["SERPINE1",]
onetwofive<-new.datfil["SPRR3",]
onetwosix<-new.datfil["STAT1",]
onetwoseven<-new.datfil["STAT5A",]
onetwoeight<-new.datfil["TAP1",]
onetwonine<-new.datfil["TFPI",]
onethreezero<-new.datfil["TGFB2",]
onethreeone<-new.datfil["THBD",]
onethreetwo<-new.datfil["TNFAIP3",]
onethreethree<-new.datfil["TPR",]
onethreefour<-new.datfil["TYROBP",]
onethreefive<-new.datfil["UBA7",]
onethreesix<-new.datfil["UBE2L6",]
onethreeseven<-new.datfil["USP18",]
onethreeeight<-new.datfil["VCAM1",]
onethreenine<-new.datfil["WAS",]
onefourzero<-new.datfil["XAF1",]
onefourone<-new.datfil["XPO1",]
onefourtwo<-new.datfil["AFAP1L2",]
onefourthree<-new.datfil["AIF1",]
onefourfour<-new.datfil["APOBEC3F",]
onefourfive<-new.datfil["APOBEC3G",]
onefoursix<-new.datfil["BNIP3",]
onefourseven<-new.datfil["CADM1",]
onefoureight<-new.datfil["CALR",]
onefournine<-new.datfil["CAMK2B",]
onefivezero<-new.datfil["CASP1",]
onefiveone<-new.datfil["CCR7",]
onefivetwo<-new.datfil["CD19",]
onefivethree<-new.datfil["CD83",]
onefivefour<-new.datfil["CDK1",]
onefivefive<-new.datfil["CEBPB",]
onefivesix<-new.datfil["CEBPG",]
onefiveseven<-new.datfil["CSF2RB",]
onefiveeight<-new.datfil["CXCL10",]
onefivenine<-new.datfil["CXCL16",]
onesixzero<-new.datfil["CXCR3",]
onesixone<-new.datfil["CYSLTR1",]
onesixtwo<-new.datfil["DEFB1",]
onesixthree<-new.datfil["EIF4A2",]
onesixfour<-new.datfil["F2",]
onesixfive<-new.datfil["F2R",]
onesixsix<-new.datfil["F5",]
onesixseven<-new.datfil["F7",]
onesixeight<-new.datfil["FAS",]
onesixnine<-new.datfil["FASLG",]
onesevenzero<-new.datfil["GAGE1",]
onesevenone<-new.datfil["GBP2",]
oneseventwo<-new.datfil["GTF2F2",]
oneseventhree<-new.datfil["GZMB",]
onesevenfour<-new.datfil["HLA-A",]
onesevenfive<-new.datfil["HLA-DMA",]
onesevensix<-new.datfil["HLA-DPB1",]
onesevenseven<-new.datfil["HLA-DRB3",]
oneseveneight<-new.datfil["HLA-E",]
onesevennine<-new.datfil["HLA-F",]
oneeightzero<-new.datfil["HP",]
oneeightone<-new.datfil["HRAS",]
oneeighttwo<-new.datfil["IFITM2",]
oneeightthree<-new.datfil["IFITM3",]
oneeightfour<-new.datfil["IL17RB",]
oneeightfive<-new.datfil["IL18",]
oneeightsix<-new.datfil["IL1R1",]
oneeightseven<-new.datfil["IL1RN",]
oneeighteight<-new.datfil["IL2RA",]
oneeightnine<-new.datfil["IL2RG",]
oneninezero<-new.datfil["IL6R",]
onenineone<-new.datfil["IL7R",]
oneninetwo<-new.datfil["INHBB",]
oneninethree<-new.datfil["IRAK1",]
oneninefour<-new.datfil["ITGB5",]
oneninefive<-new.datfil["KLRC3",]
oneninesix<-new.datfil["KLRC4",]
onenineseven<-new.datfil["KPNA2",]
onenineeight<-new.datfil["KPNA3",]
twozeronine<-new.datfil["LGALS3BP",]
twozerozero<-new.datfil["LY75",]
twozeroone<-new.datfil["LYZ",]
twozerotwo<-new.datfil["MAP2K4",]
twozerothree<-new.datfil["MGLL",]
twozerofour<-new.datfil["MIA3",]
twozerofive<-new.datfil["MICB",]
twozerosix<-new.datfil["MX2",]
twozeroseven<-new.datfil["NFATC4",]
twozeroeight<-new.datfil["NFKB2",]
twoonenine<-new.datfil["NMI",]
twoonezero<-new.datfil["NOS2",]
twooneone<-new.datfil["NRAS",]
twoonetwo<-new.datfil["NUP107",]
twoonethree<-new.datfil["NUP155",]
twoonefour<-new.datfil["NUP205",]
twoonefive<-new.datfil["NUP37",]
twoonesix<-new.datfil["NUP43",]
twooneseven<-new.datfil["NUP85",]
twooneeight<-new.datfil["NUP93",]
twotwonine<-new.datfil["OR2H2",]
twotwozero<-new.datfil["PELI3",]
twotwoone<-new.datfil["PF4",]
twotwotwo<-new.datfil["POLR2K",]
twotwothree<-new.datfil["POLR2L",]
twotwofour<-new.datfil["PRF1",]
twotwofive<-new.datfil["PRL",]
twotwosix<-new.datfil["PSG8",]
twotwoseven<-new.datfil["PSMA6",]
twotwoeight<-new.datfil["PSMB10",]
twothreenine<-new.datfil["PSMB3",]
twothreezero<-new.datfil["PSMB6",]
twothreeone<-new.datfil["PSMD1",]
twothreetwo<-new.datfil["PSMD10",]
twothreethree<-new.datfil["PSME2",]
twothreefour<-new.datfil["PTAFR",]
twothreefive<-new.datfil["PTPN1",]
twothreesix<-new.datfil["PYDC1",]
twothreeseven<-new.datfil["RBX1",]
twothreeeight<-new.datfil["RPL11",]
twofournine<-new.datfil["RPL12",]
twofourzero<-new.datfil["RPL14",]
twofourone<-new.datfil["RPL15",]
twofourtwo<-new.datfil["RPL37A",]
twofourthree<-new.datfil["RPL4",]
twofourfour<-new.datfil["RPL41",]
twofourfive<-new.datfil["RPLP1",]
twofoursix<-new.datfil["RPS11",]
twofourseven<-new.datfil["RPS14",]
twofoureight<-new.datfil["RPS18",]
twofivenine<-new.datfil["RPS23",]
twofivezero<-new.datfil["RPS27",]
twofiveone<-new.datfil["RPS28",]
twofivetwo<-new.datfil["RPS4Y1",]
twofivethree<-new.datfil["RPS6",]
twofivefour<-new.datfil["RPS8",]
twofivefive<-new.datfil["S100A12",]
twofivesix<-new.datfil["SCG2",]
twofiveseven<-new.datfil["SEC61B",]
twofiveeight<-new.datfil["SEC61G",]
twosixnine<-new.datfil["SEH1L",]
twosixzero<-new.datfil["SH2B1",]
twosixone<-new.datfil["SHC1",]
twosixtwo<-new.datfil["SOD1",]
twosixthree<-new.datfil["TCIRG1",]
twosixfour<-new.datfil["TFF3",]
twosixfive<-new.datfil["TLR3",]
twosixsix<-new.datfil["TPST1",]
twosixseven<-new.datfil["UBE2E1",]
twosixeight<-new.datfil["UBE2N",]
twosevennine<-new.datfil["UMOD",]
twosevenzero<-new.datfil["VWF",]
twosevenone<-new.datfil["APOL3",]
twoseventwo<-new.datfil["BNIP3L",]
twoseventhree<-new.datfil["C2",]
twosevenfour<-new.datfil["CCRL1",]
twosevenfive<-new.datfil["CD1D",]
twosevensix<-new.datfil["CD36",]
twosevenseven<-new.datfil["CD40",]
twoseveneight<-new.datfil["CFP",]
twoeightnine<-new.datfil["CHST2",]
twoeightzero<-new.datfil["COLEC12",]
twoeightone<-new.datfil["CSF2RA",]
twoeighttwo<-new.datfil["CSH1",]
twoeightthree<-new.datfil["CXCL5",]
twoeightfour<-new.datfil["CXCR6",]
twoeightfive<-new.datfil["DCDC2",]
twoeightsix<-new.datfil["DMBT1",]
twoeightseven<-new.datfil["ELF3",]
twoeighteight<-new.datfil["F12",]
twoninenine<-new.datfil["GBP4",]
twoninezero<-new.datfil["GH1",]
twonineone<-new.datfil["GPR68",]
twoninetwo<-new.datfil["HLA-DPA1",]
twoninethree<-new.datfil["HLA-G",]
twoninefour<-new.datfil["HOXB13",]
twoninefive<-new.datfil["IFI35",]
twoninesix<-new.datfil["IFNG",]
twonineseven<-new.datfil["IL29",]
twonineeight<-new.datfil["IL2RB",]
threezeronine<-new.datfil["KRT1",]
threezerozero<-new.datfil["LYVE1",]
threezeroone<-new.datfil["MAP3K8",]
threezerotwo<-new.datfil["MST1R",]
threezerothree<-new.datfil["NOX4",]
threezerofour<-new.datfil["PELI1",]
threezerofive<-new.datfil["PELI2",]
threezerosix<-new.datfil["PROC",]
threezeroseven<-new.datfil["PTPN6",]
threezeroeight<-new.datfil["RNASEL",]
threeonenine<-new.datfil["RPS12",]
threeonezero<-new.datfil["SP140",]
threeoneone<-new.datfil["STAB1",]
threeonetwo<-new.datfil["STAT2",]
threeonethree<-new.datfil["TNFAIP6",]
threeonefour<-new.datfil["TNIP1",]
threeonefive<-new.datfil["VAV1",]
threeonesix<-new.datfil["XCL1",]
threeoneseven<-new.datfil["XCL2",]

# aimgenes <- rbind(one,two,three,four,five,six,seven,eight,nine,onezero,oneone,onetwo,onethree,onefour,onefive,onesix,oneseven,oneeight,onenine,twozero,twoone,twotwo,twothree,twofour,twofive,twosix,twoseven,twoeight,twonine,threezero,threeone,threetwo,threethree,threefour,threefive,threesix,threeseven,threeeight,threenine,fourzero,fourone,fourtwo,fourthree,fourfour,fourfive,foursix,fourseven,foureight,fournine,fivezero,fiveone,fivetwo,fivethree,fivefour,fivefive,fivesix,fiveseven,fiveeight,fivenine,sixzero,sixone,sixtwo,sixthree,sixfour,sixfive,sixsix,sixseven,sixeight,sixnine,sevenzero,sevenone,seventwo,seventhree,sevenfour,sevenfive,sevensix,sevenseven,seveneight,sevennine,eightzero,eightone,eighttwo,eightthree,eightfour,eightfive,eightsix,eightseven,eighteight,eightnine,ninezero,nineone,ninetwo,ninethree,ninefour,ninefive,ninesix,nineseven,nineeight,ninenine,onezerozero,onezeroone,onezerotwo,onezerothree,onezerofour,onezerofive,onezerosix,onezeroseven,onezeroeight,onezeronine,oneonezero,oneoneone,oneonetwo,oneonethree,oneonefour,oneonefive,oneonesix,oneoneseven,oneoneeight,oneonenine,onetwozero,onetwoone,onetwotwo,onetwothree,onetwofour,onetwofive,onetwosix,onetwoseven,onetwoeight,onetwonine,onethreezero,onethreeone,onethreetwo,onethreethree,onethreefour,onethreefive,onethreesix,onethreeseven,onethreeeight,onethreenine,onefourzero,onefourone,onefourtwo,onefourthree,onefourfour,onefourfive,onefoursix,onefourseven,onefoureight,onefournine,onefivezero,onefiveone,onefivetwo,onefivethree,onefivefour,onefivefive,onefivesix,onefiveseven,onefiveeight,onefivenine,onesixzero,onesixone,onesixtwo,onesixthree,onesixfour,onesixfive,onesixsix,onesixseven,onesixeight,onesixnine,onesevenzero,onesevenone,oneseventwo,oneseventhree,onesevenfour,onesevenfive,onesevensix,onesevenseven,oneseveneight,onesevennine,oneeightzero,oneeightone,oneeighttwo,oneeightthree,oneeightfour,oneeightfive,oneeightsix,oneeightseven,oneeighteight,oneeightnine,oneninezero,onenineone,oneninetwo,oneninethree,oneninefour,oneninefive,oneninesix,onenineseven,onenineeight,twozeronine,twozerozero,twozeroone,twozerotwo,twozerothree,twozerofour,twozerofive,twozerosix,twozeroseven,twozeroeight,twoonenine,twoonezero,twooneone,twoonetwo,twoonethree,twoonefour,twoonefive,twoonesix,twooneseven,twooneeight,twotwonine,twotwozero,twotwoone,twotwotwo,twotwothree,twotwofour,twotwofive,twotwosix,twotwoseven,twotwoeight,twothreenine,twothreezero,twothreeone,twothreetwo,twothreethree,twothreefour,twothreefive,twothreesix,twothreeseven,twothreeeight,twofournine,twofourzero,twofourone,twofourtwo,twofourthree,twofourfour,twofourfive,twofoursix,twofourseven,twofoureight,twofivenine,twofivezero,twofiveone,twofivetwo,twofivethree,twofivefour,twofivefive,twofivesix,twofiveseven,twofiveeight,twosixnine,twosixzero,twosixone,twosixtwo,twosixthree,twosixfour,twosixfive,twosixsix,twosixseven,twosixeight,twosevennine,twosevenzero,twosevenone,twoseventwo,twoseventhree,twosevenfour,twosevenfive,twosevensix,twosevenseven,twoseveneight,twoeightnine,twoeightzero,twoeightone,twoeighttwo,twoeightthree,twoeightfour,twoeightfive,twoeightsix,twoeightseven,twoeighteight,twoninenine,twoninezero,twonineone,twoninetwo,twoninethree,twoninefour,twoninefive,twoninesix,twonineseven,twonineeight,threezeronine,threezerozero,threezeroone,threezerotwo,threezerothree,threezerofour,threezerofive,threezerosix,threezeroseven,threezeroeight,threeonenine,threeonezero,threeoneone,threeonetwo,threeonethree,threeonefour,threeonefive,threeonesix,threeoneseven)
# aimgenesnames <- c("ADM","ADORA2A","ADORA2B","ALOX5AP","ANKRD1","ANXA1","AOC3","AOX1","B2M","BCL2","C4BPB","CAMP","CCL2","CCL20","CCL26","CCL28","CCL3","CCL3L3","CCL4","CCL5","CCR9","CD44","CD81","CRP","CSF2","CTGF","CTSS","CXCL1","CXCL11","CXCL12","CXCL2","CXCL3","CXCL6","CXCL9","CXCR4","CXCR7","DCBLD2","DDX58","DEFB103A","EGR1","EIF2AK2","EIF4E","EIF4G1","EREG","FOS","FOSL1","GBP1","GBP5","GP9","HCP5","HERC5","HLA-B","HLA-C","HLA-DMB","HLA-DRB1","HSP90AA1","ICAM1","IFI27","IFI6","IFIH1","IFIT1","IFIT2","IFIT3","IFITM1","IFNGR1","IL1A","IL1B","IL1R2","IL32","IL6","IL6ST","IL8","INHBA","IRF6","IRF7","IRF8","IRF9","ISG15","ISG20","ITGAV","JAK2","KCNN4","KLK8","KLRC2","LBP","LCK","LSP1","LY96","LYN","LYST","MDK","MRC2","MT2A","MX1","NCF1","NCF2","NLRP3","NLRX1","NOD2","NUP35","OAS1","OAS2","OAS3","OASL","ORM1","ORM2","PAGE1","PIK3R2","PLA2G7","PLAT","PPBP","PROS1","PSMA3","PSMB8","PSMB9","PSMC6","PTX3","RPL26","RPL38","RSAD2","S100A7","S100A8","S100A9","SERPINE1","SPRR3","STAT1","STAT5A","TAP1","TFPI","TGFB2","THBD","TNFAIP3","TPR","TYROBP","UBA7","UBE2L6","USP18","VCAM1","WAS","XAF1","XPO1","AFAP1L2","AIF1","APOBEC3F","APOBEC3G","BNIP3","CADM1","CALR","CAMK2B","CASP1","CCR7","CD19","CD83","CDK1","CEBPB","CEBPG","CSF2RB","CXCL10","CXCL16","CXCR3","CYSLTR1","DEFB1","EIF4A2","F2","F2R","F5","F7","FAS","FASLG","GAGE1","GBP2","GTF2F2","GZMB","HLA-A","HLA-DMA","HLA-DPB1","HLA-DRB3","HLA-E","HLA-F","HP","HRAS","IFITM2","IFITM3","IL17RB","IL18","IL1R1","IL1RN","IL2RA","IL2RG","IL6R","IL7R","INHBB","IRAK1","ITGB5","KLRC3","KLRC4","KPNA2","KPNA3","LGALS3BP","LY75","LYZ","MAP2K4","MGLL","MIA3","MICB","MX2","NFATC4","NFKB2","NMI","NOS2","NRAS","NUP107","NUP155","NUP205","NUP37","NUP43","NUP85","NUP93","OR2H2","PELI3","PF4","POLR2K","POLR2L","PRF1","PRL","PSG8","PSMA6","PSMB10","PSMB3","PSMB6","PSMD1","PSMD10","PSME2","PTAFR","PTPN1","PYDC1","RBX1","RPL11","RPL12","RPL14","RPL15","RPL37A","RPL4","RPL41","RPLP1","RPS11","RPS14","RPS18","RPS23","RPS27","RPS28","RPS4Y1","RPS6","RPS8","S100A12","SCG2","SEC61B","SEC61G","SEH1L","SH2B1","SHC1","SOD1","TCIRG1","TFF3","TLR3","TPST1","UBE2E1","UBE2N","UMOD","VWF","APOL3","BNIP3L","C2","CCRL1","CD1D","CD36","CD40","CFP","CHST2","COLEC12","CSF2RA","CSH1","CXCL5","CXCR6","DCDC2","DMBT1","ELF3","F12","GBP4","GH1","GPR68","HLA-DPA1","HLA-G","HOXB13","IFI35","IFNG","IL29","IL2RB","KRT1","LYVE1","MAP3K8","MST1R","NOX4","PELI1","PELI2","PROC","PTPN6","RNASEL","RPS12","SP140","STAB1","STAT2","TNFAIP6","TNIP1","VAV1","XCL1","XCL2")
# 3 genes missing
aimgenes <- rbind(one,two,three,four,five,six,seven,eight,nine,onezero,oneone,onetwo,onethree,onefour,onefive,onesix,oneseven,oneeight,onenine,twozero,twoone,twotwo,twothree,twofour,twofive,twosix,twoseven,twoeight,twonine,threezero,threeone,threetwo,threethree,threefour,threefive,threesix,threeseven,threeeight,threenine,fourzero,fourone,fourtwo,fourthree,fourfour,fourfive,foursix,fourseven,foureight,fournine,fivezero,fiveone,fivetwo,fivethree,fivefour,fivefive,fivesix,fiveseven,fiveeight,fivenine,sixzero,sixone,sixtwo,sixthree,sixfour,sixfive,sixsix,sixseven,sixeight,sixnine,sevenzero,sevenone,seventwo,seventhree,sevenfour,sevenfive,sevensix,sevenseven,seveneight,sevennine,eightzero,eightone,eighttwo,eightthree,eightfive,eightsix,eightseven,eighteight,eightnine,ninezero,nineone,ninetwo,ninethree,ninefour,ninefive,ninesix,nineseven,nineeight,ninenine,onezerozero,onezeroone,onezerotwo,onezerothree,onezerofour,onezerofive,onezerosix,onezeroseven,onezeroeight,onezeronine,oneonezero,oneoneone,oneonetwo,oneonethree,oneonefour,oneonefive,oneonesix,oneoneseven,oneoneeight,oneonenine,onetwozero,onetwoone,onetwotwo,onetwothree,onetwofour,onetwofive,onetwosix,onetwoseven,onetwoeight,onetwonine,onethreezero,onethreeone,onethreetwo,onethreethree,onethreefour,onethreefive,onethreesix,onethreeseven,onethreeeight,onethreenine,onefourzero,onefourone,onefourtwo,onefourthree,onefourfour,onefourfive,onefoursix,onefourseven,onefoureight,onefournine,onefivezero,onefiveone,onefivetwo,onefivethree,onefivefour,onefivefive,onefivesix,onefiveseven,onefiveeight,onefivenine,onesixzero,onesixone,onesixtwo,onesixthree,onesixfour,onesixfive,onesixsix,onesixseven,onesixeight,onesixnine,onesevenzero,onesevenone,oneseventwo,oneseventhree,onesevenfour,onesevenfive,onesevensix,onesevenseven,oneseveneight,onesevennine,oneeightzero,oneeightone,oneeighttwo,oneeightthree,oneeightfour,oneeightfive,oneeightsix,oneeightseven,oneeighteight,oneeightnine,oneninezero,onenineone,oneninetwo,oneninethree,oneninefour,oneninefive,oneninesix,onenineseven,onenineeight,twozeronine,twozerozero,twozerotwo,twozerothree,twozerofour,twozerofive,twozerosix,twozeroseven,twozeroeight,twoonenine,twoonezero,twooneone,twoonetwo,twoonethree,twoonefour,twoonefive,twoonesix,twooneseven,twooneeight,twotwonine,twotwozero,twotwoone,twotwotwo,twotwothree,twotwofour,twotwofive,twotwosix,twotwoseven,twotwoeight,twothreenine,twothreezero,twothreeone,twothreetwo,twothreethree,twothreefour,twothreefive,twothreesix,twothreeseven,twothreeeight,twofournine,twofourzero,twofourone,twofourtwo,twofourthree,twofourfour,twofourfive,twofoursix,twofourseven,twofoureight,twofivenine,twofivezero,twofiveone,twofivetwo,twofivethree,twofivefour,twofivefive,twofivesix,twofiveseven,twofiveeight,twosixnine,twosixzero,twosixone,twosixtwo,twosixthree,twosixfour,twosixfive,twosixsix,twosixseven,twosixeight,twosevennine,twosevenzero,twosevenone,twoseventwo,twoseventhree,twosevenfour,twosevenfive,twosevensix,twosevenseven,twoseveneight,twoeightnine,twoeightzero,twoeightone,twoeighttwo,twoeightthree,twoeightfour,twoeightfive,twoeightsix,twoeightseven,twoeighteight,twoninenine,twoninezero,twonineone,twoninetwo,twoninethree,twoninefour,twoninefive,twoninesix,twonineseven,twonineeight,threezeronine,threezerozero,threezeroone,threezerotwo,threezerothree,threezerofour,threezerofive,threezerosix,threezeroseven,threezeroeight,threeonenine,threeonezero,threeoneone,threeonetwo,threeonethree,threeonefour,threeonefive,threeonesix)
aimgenesnames <- c("ADM","ADORA2A","ADORA2B","ALOX5AP","ANKRD1","ANXA1","AOC3","AOX1","B2M","BCL2","C4BPB","CAMP","CCL2","CCL20","CCL26","CCL28","CCL3","CCL3L3","CCL4","CCL5","CCR9","CD44","CD81","CRP","CSF2","CTGF","CTSS","CXCL1","CXCL11","CXCL12","CXCL2","CXCL3","CXCL6","CXCL9","CXCR4","CXCR7","DCBLD2","DDX58","DEFB103A","EGR1","EIF2AK2","EIF4E","EIF4G1","EREG","FOS","FOSL1","GBP1","GBP5","GP9","HCP5","HERC5","HLA-B","HLA-C","HLA-DMB","HLA-DRB1","HSP90AA1","ICAM1","IFI27","IFI6","IFIH1","IFIT1","IFIT2","IFIT3","IFITM1","IFNGR1","IL1A","IL1B","IL1R2","IL32","IL6","IL6ST","IL8","INHBA","IRF6","IRF7","IRF8","IRF9","ISG15","ISG20","ITGAV","JAK2","KCNN4","KLK8","LBP","LCK","LSP1","LY96","LYN","LYST","MDK","MRC2","MT2A","MX1","NCF1","NCF2","NLRP3","NLRX1","NOD2","NUP35","OAS1","OAS2","OAS3","OASL","ORM1","ORM2","PAGE1","PIK3R2","PLA2G7","PLAT","PPBP","PROS1","PSMA3","PSMB8","PSMB9","PSMC6","PTX3","RPL26","RPL38","RSAD2","S100A7","S100A8","S100A9","SERPINE1","SPRR3","STAT1","STAT5A","TAP1","TFPI","TGFB2","THBD","TNFAIP3","TPR","TYROBP","UBA7","UBE2L6","USP18","VCAM1","WAS","XAF1","XPO1","AFAP1L2","AIF1","APOBEC3F","APOBEC3G","BNIP3","CADM1","CALR","CAMK2B","CASP1","CCR7","CD19","CD83","CDK1","CEBPB","CEBPG","CSF2RB","CXCL10","CXCL16","CXCR3","CYSLTR1","DEFB1","EIF4A2","F2","F2R","F5","F7","FAS","FASLG","GAGE1","GBP2","GTF2F2","GZMB","HLA-A","HLA-DMA","HLA-DPB1","HLA-DRB3","HLA-E","HLA-F","HP","HRAS","IFITM2","IFITM3","IL17RB","IL18","IL1R1","IL1RN","IL2RA","IL2RG","IL6R","IL7R","INHBB","IRAK1","ITGB5","KLRC3","KLRC4","KPNA2","KPNA3","LGALS3BP","LY75","MAP2K4","MGLL","MIA3","MICB","MX2","NFATC4","NFKB2","NMI","NOS2","NRAS","NUP107","NUP155","NUP205","NUP37","NUP43","NUP85","NUP93","OR2H2","PELI3","PF4","POLR2K","POLR2L","PRF1","PRL","PSG8","PSMA6","PSMB10","PSMB3","PSMB6","PSMD1","PSMD10","PSME2","PTAFR","PTPN1","PYDC1","RBX1","RPL11","RPL12","RPL14","RPL15","RPL37A","RPL4","RPL41","RPLP1","RPS11","RPS14","RPS18","RPS23","RPS27","RPS28","RPS4Y1","RPS6","RPS8","S100A12","SCG2","SEC61B","SEC61G","SEH1L","SH2B1","SHC1","SOD1","TCIRG1","TFF3","TLR3","TPST1","UBE2E1","UBE2N","UMOD","VWF","APOL3","BNIP3L","C2","CCRL1","CD1D","CD36","CD40","CFP","CHST2","COLEC12","CSF2RA","CSH1","CXCL5","CXCR6","DCDC2","DMBT1","ELF3","F12","GBP4","GH1","GPR68","HLA-DPA1","HLA-G","HOXB13","IFI35","IFNG","IL29","IL2RB","KRT1","LYVE1","MAP3K8","MST1R","NOX4","PELI1","PELI2","PROC","PTPN6","RNASEL","RPS12","SP140","STAB1","STAT2","TNFAIP6","TNIP1","VAV1","XCL1")
length(aimgenesnames)
rownames(aimgenes) <- aimgenesnames
tempdat <- aimgenes
dim(tempdat)

# for(i in 1:nrow(aimgenes)){
#   colors <- rainbow(nrow(aimgenes))
#   genes <- rbind(aimgenes[i,])
#   genename <- as.character(rownames(aimgenes)[i])
#   rownames(genes) <- genename
#   color <- colors[i]
#   fn = as.character(paste("/home/steve/.gvfs//onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/m084b trial/aimgenes/","aim_",genename,".png", sep=""))
#   png(filename = fn, width = 1000, height = 1000)  
#   plot(c(1,ncol(r)),range(genes),type='n',main=paste("Profile plot of aimgene",rownames(aimgenes)[i], sep="_"),xlab="",ylab="Expression intensity",axes=F)
#   axis(side=1,at=c(1:ncol(genes)),labels=colnames(genes),cex.axis=0.8,las=2)
#   axis(side=2)
#   for(j in 1:nrow(genes)) {
#     ycoord <- as.numeric(genes[j,])
#     lines(c(1:ncol(genes)),ycoord,col=color,lwd=2)
#     points(c(1:ncol(genes)),ycoord,col=color,pch=19, lwd=1)
#   }
#   #legend for genes
#   legend("topright", legend=rownames(genes), fill=color, bg="white", cex = .8, ncol=1)
#   dev.off()
# }

library(gplots)
# ALL SAMPLES
newgenes = aimgenes

# PAIRS ONLY 
newgenes = aimgenes[,c(1,17,2,18,3,19,4,20,5,21,6,22,7,23,8,24,9,25,10,26,11,27,12,28,13,29,14,30,15,31,16,32)]
colnames(newgenes)
dim(newgenes)

# ALL PRES
plotgenes <- newgenes #or rnewgenes
title <- "pre only, AIM genes"
title <- "all samples, AIM genes"
title <- "all pairs, AIM genes"

colbreak = c(seq(-2,-.5,length=100),seq(-.5,.5,length=100),seq(.5,2,length=100))
colbreak = c(seq(-2,-.75,length=100),seq(-.75,.75,length=100),seq(.75,2,length=100))
colbreak = c(seq(-3,-1,length=100),seq(-1,1,length=100),seq(1,3,length=100))
mycolors <- colorRampPalette(c("red", "black", "green"))(n = 299)
par(cex.main=.8)

# pre-only #NEED TO SELECT 'pos' FOR ALL SAMPLES
ht <- heatmap.2(plotgenes[,1:35], main=title, col=mycolors, breaks=colbreak, density.info="none", trace="none", srtCol = 45, margins=c(7,10), cexCol= .7, offsetCol = -.6, dendrogram=c("none"), symm=F,symkey=F,symbreaks=T, scale="row", Colv = NA)
# all samples #NEED TO SELECT 'pos' FOR ALL SAMPLES
ht <- heatmap.2(plotgenes, main=title, col=mycolors, breaks=colbreak, density.info="none", trace="none", srtCol = 45, margins=c(7,10), cexCol= .7, offsetCol = -.6, dendrogram=c("none"), symm=F,symkey=F,symbreaks=T, scale="row", Colv = NA)

##################################################################################################################################################################################################################################################################################################################################
# BARPLOT
# tran <- as.matrix(t(newgenes))
# tran
# rownames(newgenes)[1:12]
# dim(newgenes)
# 
# for(i in 1:ncol(tran)){
#   g <- tran[,i]
#   genename <- rownames(newgenes)[i]
#   #dev.off()
#   fn = as.character(paste("/home/steve/.gvfs//onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/m084b trial/aimgenes/","aim_",genename,"_r.png", sep=""))
#   png(filename = fn, width = 1000, height = 1000)
#   bp <- barplot(g, xlab = "", main=paste("Profile plot of aim",rownames(newgenes)[i], sep="_"), ylab = "Expression intensity", axisnames=F, col=c(1,1,32,32))#rainbow(2))
#   #   labs <- rownames(tran)
#   #   text(cex=.7, x=x-1.25, labs, xpd=TRUE, srt=55)
#   axis(1, at=bp, labels=rownames(tran), cex.axis=.7, las=2)
#   axis(side=2)
#   dev.off()
# }

################################################
# ALLGENES #####################################
################################################
library(gplots)
new.datfil <- r
new.datfil <- new.dat
new.datfil <- new.datpairs
new.datfil <- new.dat2mo
new.datfil <- new.dat2mopres
new.datfil <- new.dat6mo
new.datfil <- new.dat9mo
new.datfil <- new.dat12mo
new.datfil <- new.dat2mopair
new.datfil <- new.dat6mopair
new.datfil <- new.dat9mopair
new.datfil <- new.dat12mopair

new.datfil <- log2(new.dat)
new.datfil <- log2(new.datpair)
new.datfil <- log2(new.dat2mo)
new.datfil <- log2(new.dat6mo)
new.datfil <- log2(new.dat9mo)
new.datfil <- log2(new.dat12mo)
new.datfil <- log2(new.dat2mopair)
new.datfil <- log2(new.dat6mopair)
new.datfil <- log2(new.dat9mopair)
new.datfil <- log2(new.dat12mopair)
new.datfil <- log2(new.datpre)
new.datfil <- log2(new.dat2mopre)
new.datfil <- log2(new.dat6mopre)
new.datfil <- log2(new.dat9mopre)
new.datfil <- log2(new.dat12mopre)
title <- paste("log pres only, more/less 6 mo survival\n", "limma-filtered genes, p<", pval, sep="")
title <- paste("log pres only, more/less 9 mo survival\n", "limma-filtered genes, p<", pval, sep="")
title <- paste("log pres only, more/less 12 mo survival\n", "limma-filtered genes, p<", pval, sep="")
title <- paste("log pres only, more/less 2mo therapy\n", "limma-filtered genes, p<", pval, sep="")
title <- paste("log pres only, all samples\n", "limma-filtered genes, p<", pval, sep="")
pval <- .01
title <- paste("log pres only, more/less 6 mo survival\n", "limma-filtered genes, p<", pval, sep="")

f <- factor(as.character(sub("-\\w+$", "", colnames(new.datfil)))) # for all samples, new.dat
f <- factor(as.character(sub("-\\w+-\\w+$", "", colnames(new.datfil))))

f
f2 <- factor(as.character(sub("-\\w+-\\w+$", "", colnames(new.dat2mopre))))
library(gplots)
design <- model.matrix(~f)
design
fit <- eBayes(lmFit(new.datfil, design))
selected  <- fit$p.value[, 2] < pval # or p<.01 #p.adjust(fit$p.value[, 2], method="BY")<.05
esetSel <- new.datfil[selected,]
dim(esetSel)
############
colnames(esetSel)
colnames(esetSel) <- sub("LESS-","",colnames(esetSel))
colnames(esetSel) <- sub("MORE-","*******",colnames(esetSel))
colnames(esetSel)
############

title
colbreak = c(seq(-2,-.75,length=100),seq(-.75,.75,length=100),seq(.75,2,length=100))
# colbreak = c(seq(-2,-.5,length=100),seq(-.5,.5,length=100),seq(.5,2,length=100))
mycolors <- colorRampPalette(c("red", "black", "green"))(n = 299)
par(cex.main=.8)
# pre-only #NEED TO SELECT 'pos' FOR ALL SAMPLES
# ht <- heatmap.2(esetSel[,1:35], main=title, col=mycolors, breaks=colbreak, density.info="none", trace="none", srtCol = 45, margins=c(7,10), cexCol= .7, offsetCol = -.6, dendrogram=c("none"), symm=F,symkey=F,symbreaks=T, scale="row", Colv = NA, cexRow=.5)
num_genes <- length(rownames(esetSel))
# all-samples #NEED TO SELECT 'pos' FOR ALL SAMPLES
filename <- gsub(", ", "_", title)
filename <- sub("p<.","p",filename)
filename <- sub("\n","_",filename)
filename <- sub("more/less","",filename)
filename <- gsub("\\s+", "", filename)
filename <- paste(filename, "_", num_genes, ".png", sep="")
title
filename

color.map <- function(classfactor) { if (classfactor=="LESS") "pink" else "lightblue" }
patientcolors <- unlist(lapply(f2, color.map))
patientcolors

dist.pear <- function(x) as.dist(1-cor(t(x)))
hclust.ave <- function(x) hclust(x, method="average")
##########
# dist.pear <- as.dist(1-cor(t(esetSel)))
# hclust.ave <- hclust(dist.pear, method="average")
# plot(hclust.ave, cex=.5)
# heatmap(esetSel, Rowv=as.dendrogram(hclust.ave), Colv=NA, scale="row")
# heatmap(log2(esetSel), Rowv=NA, Colv=as.dendrogram(hclust.ave))
# heatmap(log2(esetSel), Rowv=NA, Colv=as.dendrogram(hclust.ave))
##########

png(filename = filename, width=1000, height=800, units="px")
# ht <- heatmap.2(esetSel, main=title, col=mycolors, breaks=colbreak, density.info="none", trace="none", srtCol = 45, margins=c(7,10), offsetCol = -.6, dendrogram=c("none"), symm=F,symkey=F,symbreaks=T, scale="row", Colv = NA, cexRow=1, cexCol=1)
ht <- heatmap.2(esetSel, ColSideColors = patientcolors, main=title, col=mycolors, breaks=colbreak, density.info="none", trace="none", srtCol = 45, margins=c(7,10), offsetCol = -.6, dendrogram=c("none"), symm=F,symkey=F,symbreaks=T, scale="row", Colv = NA, cexRow=1, cexCol=1)
legend("left",c("less than 2mo therapy", "more than 2mo therapy"),col=c("pink","lightblue"),pch=10,cex=1,horiz=F)
dev.off()
dend_filename <- sub("_limma", "_dendrogram2_limma",filename)
dend_filename
png(filename = dend_filename, width=1000, height=800, units="px")
ht <- heatmap.2(esetSel, ColSideColors = patientcolors, main=title, col=mycolors, breaks=colbreak, density.info="none", trace="none", srtCol = 45, margins=c(7,10), offsetCol = -.6, dendrogram=c("column"), symm=F,symkey=F,symbreaks=T, scale="row", cexRow=1, cexCol=1)
heatmap.2(esetSel, main = title, ColSideColors = patientcolors, trace = "none", scale="row", col=mycolors, breaks=colbreak, density.info="none", srtCol = 45, margins=c(7,10), offsetCol = -.6, symm=F,symkey=F,symbreaks=T, cexRow=1, cexCol=1)
#heatmap.2(esetSel, main = title, ColSideColors = patientcolors, distfun=dist.pear, hclustfun=hclust.ave, trace = "none", scale="row", col=mycolors, breaks=colbreak, density.info="none", srtCol = 45, margins=c(7,10), offsetCol = -.6, symm=F,symkey=F,symbreaks=T, cexRow=1, cexCol=1)
legend("topright",c("less than 2mo therapy", "more than 2mo therapy"),col=c("pink","lightblue"),pch=15,cex=1,horiz=F)
dev.off()
heatmapgenes <- rownames(esetSel)[ht$rowInd]
heatmapfilename <- gsub(", ", "_", title)
heatmapfilename <- sub("p<.","p",heatmapfilename)
heatmapfilename <- sub("\n","_",heatmapfilename)
heatmapfilename <- sub("more/less","",heatmapfilename)
heatmapfilename <- gsub("\\s+", "", heatmapfilename)
heatmapfilename <- paste(heatmapfilename, "_", num_genes, "_genenames.txt", sep="")
heatmapfilename
write.table(heatmapgenes, file = heatmapfilename, quote=FALSE, row.names=FALSE)#paste(title,".csv",sep=""))



clust <- hclust(dist(esetSel))
plot(clust, cex = .5)
prop <- propexpr(esetSel)
#clustering
plotMDS(esetSel, cex=.7, labels = factor(as.character(sub("-\\w+$", "", colnames(new.datfil)))))#, col=c(rep("black",16), rep("red",16)), labels= c(rep("pre",16), rep("post",16)))
plotMDS(esetSel, main = title, cex=.7, col=c(rep("black",16), rep("red",16)), labels= c(rep("pre",16), rep("post",16)), xlab=NA, ylab=NA, )

# PCA
#dat.rmmis <- new.dat
dat.rmmis <- esetSel
dim(dat.rmmis)
datpca <- prcomp(t(dat.rmmis), scale=TRUE)
title <- "PCA of colon trial data\npre vs post, nl/scale=T"
datpca <- prcomp(t(dat.rmmis), scale=FALSE)
title <- "PCA of colon trial data (new.dat)\npre vs post, nl/scale=F"
# datpca <- prcomp(t(log2(new.dat)), scale=TRUE)
# title <- "PCA of colon trial data\npre vs post, l2/scale=T"
# datpca <- prcomp(t(log2(new.dat)), scale=FALSE)
# title <- "PCA of colon trial data\npre vs post, l2/scale=F"


plot(
  range(datpca$x[, 1]),range(datpca$x[, 2]), 
  type = "n", xlab = "pre", ylab = "post", 
  main = title
)
dimnames(datpca)
points(datpca$x[, 1][targets$Cy3 == "pre"], datpca$x[, 2][targets$Cy3 == "pre"], bg = "Red", pch = 21)
# text(datpca$x[, 1][targets$Cy3 == "pre"], datpca$x[, 2][targets$Cy3 == "pre"], label=names(datpca$x[,1][targets$Cy3 == "pre"]),pos=1,cex=0.5)
#text(datpca$x[, 1][targets$Cy3 == "pre"], datpca$x[, 2][targets$Cy3 == "pre"], label="pre",pos=1,cex=0.5)
points(datpca$x[, 1][targets$Cy5 == "post"], datpca$x[, 2][targets$Cy5 == "post"], bg = "Blue", pch=21)
# text(datpca$x[, 1][targets$Cy5 == "post"], datpca$x[, 2][targets$Cy5 == "post"], label=names(datpca$x[,1][targets$Cy5 == "post"]),pos=1,cex=0.5)
#text(datpca$x[, 1][targets$Cy5 == "post"], datpca$x[, 2][targets$Cy5 == "post"], label="post",pos=1,cex=0.5)
leg.names <- c("pre", "post")       
legend("topright",leg.names,col=c("red","blue"),pch=15,cex=.7,horiz=F)

##################################################################################################################################
#dat.rmmis <- r[, -c(which(colnames(new.dat)=="pre-PH1816"))]
# new.datfil <- new.dat12mo
# f <- factor(as.character(sub("-\\w+-\\w+$", "", colnames(new.datfil))))
# f
# design <- model.matrix(~f)
# design
# fit <- eBayes(lmFit(new.datfil, design))
# selected  <- fit$p.value[, 2] <.05 # or p<.01 #p.adjust(fit$p.value[, 2], method="BY")<.05
# limma <- new.datfil[selected,]
# dim(limma)

# dat.rmmis <- limma[, -c(which(colnames(new.dat)=="pre-PH1816"))]
# dat.rmmis <- new.dat[, -c(which(colnames(new.dat)=="pre-PH1816"))]
# dat.rmmis <- new.datfil
dat.rmmis <- new.dat12mopair
dim(dat.rmmis)
isscale <- TRUE
datpca <- prcomp(t(dat.rmmis), scale=isscale)
summary(datpca)
months <- 12
ispair <- "pairs only"
morelab <- paste(">", months, "mo", sep="")
morelab
lesslab <- paste("<", months, "mo", sep="")
title <- paste("PCA of colon trial data (", ispair, ")\n", morelab, "vs", lesslab, "survival, nl/scale =", isscale, sep = " ")
title
# datpca <- prcomp(t(log2(new.dat)), scale=TRUE)
# title <- "PCA of colon trial data\npre vs post, l2/scale=T"
# datpca <- prcomp(t(log2(new.dat)), scale=FALSE)
# title <- "PCA of colon trial data\npre vs post, l2/scale=F"

plot(
  range(datpca$x[, 1]),range(datpca$x[, 2]), 
  type = "n", xlab = morelab , ylab = lesslab, 
  main = title
)

whichmonth <- NULL
if(months == 12){
  whichmonth <- meta$twelvemonthsoverall
}else if(months == 6){
  whichmonth <- meta$sixmonthsoverall
}else if(months == 9){
  whichmonth <- meta$ninemonthsoverall
}else if(months == 2){
  whichmonth <- meta$twomonthstherapy
}
#points(datpca$x[, 1][targets$Cy3 == "pre"], datpca$x[, 2][targets$Cy3 == "pre"], bg = "Red", pch = 21)
points(datpca$x[, 1][whichmonth == 0], datpca$x[, 2][whichmonth == 0], bg = "Red", pch = 21)
points(datpca$x[, 1][whichmonth == 1], datpca$x[, 2][whichmonth == 1], bg = "Blue", pch = 21)
# data.frame(datpca$x[,2],datpca$x[,1], whichmonth)
length(whichmonth)
# # text(datpca$x[, 1][targets$Cy3 == "pre"], datpca$x[, 2][targets$Cy3 == "pre"], label=names(datpca$x[,1][targets$Cy3 == "pre"]),pos=1,cex=0.5)
# text(datpca$x[, 1][whichmonth == 0], datpca$x[, 2][whichmonth == 0], label=names(datpca$x[,1][whichmonth == 0]),pos=1,cex=0.5)#label="less",pos=1,cex=0.5)
# # text(datpca$x[, 1][targets$Cy5 == "post"], datpca$x[, 2][targets$Cy5 == "post"], label=names(datpca$x[,1][targets$Cy5 == "post"]),pos=1,cex=0.5)
# text(datpca$x[, 1][whichmonth == 1], datpca$x[, 2][whichmonth == 1], label=names(datpca$x[,1][whichmonth == 1]),pos=1,cex=0.5)#label="more",pos=1,cex=0.5)
leg.names <- c(lesslab, morelab)       
legend("topleft",leg.names,col=c("red","blue"),pch=15,cex=.7,horiz=F)

# #points(datpca$x[, 1][targets$Cy3 == "pre"], datpca$x[, 2][targets$Cy3 == "pre"], bg = "Red", pch = 21)
# points(datpca$x[, 1][meta$twelvemonthsoverall == 0], datpca$x[, 2][meta$twelvemonthsoverall == 0], bg = "Red", pch = 21)
# points(datpca$x[, 1][meta$twelvemonthsoverall == 1], datpca$x[, 2][meta$twelvemonthsoverall == 1], bg = "Blue", pch = 21)
# data.frame(datpca$x[,2],datpca$x[,1], meta$twelvemonthsoverall)
# length(meta$twelvemonthsoverall)
# # text(datpca$x[, 1][targets$Cy3 == "pre"], datpca$x[, 2][targets$Cy3 == "pre"], label=names(datpca$x[,1][targets$Cy3 == "pre"]),pos=1,cex=0.5)
# text(datpca$x[, 1][meta$twelvemonthsoverall == 0], datpca$x[, 2][meta$twelvemonthsoverall == 0], label=names(datpca$x[,1][meta$twelvemonthsoverall == 0]),pos=1,cex=0.5)#label="less",pos=1,cex=0.5)
# # text(datpca$x[, 1][targets$Cy5 == "post"], datpca$x[, 2][targets$Cy5 == "post"], label=names(datpca$x[,1][targets$Cy5 == "post"]),pos=1,cex=0.5)
# text(datpca$x[, 1][meta$twelvemonthsoverall == 1], datpca$x[, 2][meta$twelvemonthsoverall == 1], label=names(datpca$x[,1][meta$twelvemonthsoverall == 1]),pos=1,cex=0.5)#label="more",pos=1,cex=0.5)
# leg.names <- c(lesslab, morelab)       
# legend("topleft",leg.names,col=c("red","blue"),pch=15,cex=.7,horiz=F)