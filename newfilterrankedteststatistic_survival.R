# r <- 20
# c <- 5
# di_matrix <- matrix(round(runif(r * c, 1, 100)), ncol = c)
# #get median per gene
# genes_median <- apply(di_matrix, 1, median)
# #convert to 0 and 1
# di_matrix01 <- ifelse(di_matrix > genes_median, 1, 0)
# di_matrix
# genes_median
# t(di_matrix01)


# PRIMARY CODE FOR PRELIM2 ANALYSIS 

x2 1        post-PH1876        pre-PH1876
x2 2         pre-PH1900       post-PH1900
x1 3         pre-PH1811        pre-PH1816
x1 4         pre-PH1636       post-PH1640
x1 5         pre-PH1892        pre-PH1902
x1 6         pre-PH1622        pre-PH1631
x2 7         pre-JHH005       post-JHH005
8        post-PH1910        pre-PH1910
x9         pre-PH1604        pre-PH1606
x1 10       post-JHH004        pre-JHH004
x2 11       post-PH1612        pre-PH1612
x1 12        pre-PH1913        pre-PH1886
x1 13        pre-PH1635        pre-PH1644
14        pre-PH1861       post-PH1861
x2 15        pre-PH1616       post-PH1616
x3 16       post-PH1844        pre-PH1844    #NO DRY ICE
x1 17        pre-PH1623        pre-PH1632
x2 18        pre-PH1827       post-PH1827
x3 19       post-PH1815        pre-PH1815
20       post-PH1544        pre-PH1544
x2 21       post-PH1843        pre-PH1843
x3 22        pre-PH1871       post-PH1871
x1 23        pre-PH1868        pre-PH1887
24        pre-PH1545       post-PH1545
25       post-PH1869        pre-PH1869
x1 26        pre-PH1550        pre-PH1600
# 
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

#FINAL
rm(list=ls())
gc()
library(limma)
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
meta$twelvemonthsoverall

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
as.data.frame(sub(".*_1_", "", dat$targets$FileName))
# > data.frame(targets$Cy3_Sample, targets$Cy5_sample)
1              JHH004             JHH004
2              PH1550             PH1600
3              JHH005             JHH005
4              PH1604             PH1606
5              PH1622             PH1631
6              PH1623             PH1632
7              PH1815             PH1815
8              PH1635             PH1644
9              PH1827             PH1827
10             PH1811             PH1816
11             PH1636             PH1640
12             PH1892             PH1902
13             PH1612             PH1612
14             PH1868             PH1887
15             PH1616             PH1616
16             PH1910             PH1910
17             PH1843             PH1843
18             PH1913             PH1886
19             PH1844             PH1844
20             PH1861             PH1861
21             PH1869             PH1869
22             PH1871             PH1871
23             PH1876             PH1876
24             PH1641             PH1641
25             PH1824             PH1824
26             PH1545             PH1545
27             PH1565             PH1565
28             PH1910             PH1910
29             PH1641             PH1640
30             PH1612             PH1824
31             PH1861             PH1900



#1 - no changes
# ALL SAMPLEs
# pos <- c(1:26)
 pos <- c(1,2,7,8,11,14:16,18:22,24,25,3:6,9,10,12,13,17,23,26)
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

dat <- dat[,pos]
colnames(dat)
dim(dat)
targets <- targets[pos,]
dim(targets)
targets

as.data.frame(cbind(as.matrix(paste(targets$Cy3, targets$Cy3_Sample, sep="-")),as.matrix(paste(targets$Cy5, targets$Cy5_sample, sep="-"))))

# BG SUBTRACTION
# dat <- backgroundCorrect(dat, method="subtract")

#1,2,3,4
normwithin <-normalizeWithinArrays(dat,method='loess',bc.method='normexp', offset=50)
normbetween <-normalizeBetweenArrays(normwithin,method='Aquantile')

#Remove controls from normwithin/between
normwithin <- normwithin[normbetween$genes$ControlType==0,]
normbetween <- normbetween[normbetween$genes$ControlType==0,]

#Convert MA back to RG
RGb <- RG.MA(normbetween)

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

library(genefilter)
#rsd <- rowSds(dat.m)
#rsd <- rowSds(dat)
#dat <- cbind(cy3, cy5)

#1,2
dim(normwithin)
colnames(normwithin)
dat <- normwithin
tar2 <- targets
design <- modelMatrix(tar2, ref="pre")
#3
### [1] PH1876 post-Cy3_PH1876 pre-Cy5 - filtered
### [2] PH1900 pre-Cy3_PH1900 post-Cy5 - degraded, filtered
### [3] PH1811 pre-Cy3_PH1816 pre-Cy5 - non-pair
### [4] PH1636 pre-Cy3_PH1640 post-Cy5 - non-pair
### [5] PH1892 pre-Cy3_PH1902 pre-Cy5 - non-pair
### [6] PH1622 pre-Cy3_PH1631 pre-Cy5 - non-pair
### [7] JHH005 pre-Cy3_JHH005 post-Cy5 - degraded
# [8] PH1910 post-Cy3_PH1910 pre-Cy5
### [9] PH1604 pre-Cy3_PH1606 pre-Cy5 - non-pair
### [10] JHH004 post-Cy3_JHH004 pre-Cy5 - degraded sample
# [11] PH1612 post-Cy3_PH1612 pre-Cy5
### [12] PH1913 pre-Cy3_PH1886 pre-Cy5 - non-pair
### [13] PH1635 pre-Cy3_PH1644 pre-Cy5 - non-pair
# [14] PH1861 pre-Cy3_PH1861 post-Cy5
#2# [15] PH1616 pre-Cy3_PH1616 post-Cy5 - degraded
#2# [16] PH1844 post-Cy3_PH1844 pre-Cy5
### [17] PH1623 pre-Cy3_PH1632 pre-Cy5 - non-pair
# [18] PH1827 pre-Cy3_PH1827 post-Cy5
# [19] PH1815 post-Cy3_PH1815 pre-Cy5
### [20] PH1544 post-Cy3_PH1544 pre-Cy5 - filtered
#2# [21] PH1843 post-Cy3_PH1843 pre-Cy5 
#2# [22] PH1871 pre-Cy3_PH1871 post-Cy5
### [23] PH1868 pre-Cy3_PH1887 pre-Cy5 - non-pair
# [24] PH1545 pre-Cy3_PH1545 post-Cy5
### [25] PH1869 post-Cy3_PH1869 pre-Cy5 - pre sample is too degraded
### [26] PH1550 pre-Cy3_PH1600 pre-Cy5 - non-pair

# colnames(normwithin)
# #3
# filter1 <- c(1,2,3,4,5,6,9,10,12,13,17,20,23,25,26)
# #4
# filter2 <- c(1,2,3,4,5,6,9,10,12,13,15,16,17,20,21,22,23,25,26)
# dat <- normwithin[,-filter1]
# tar2 <- targets[-filter1,]
# 
# dat <- normwithin[,-filter2]
# tar2 <- targets[-filter2,]
# 
###############
#1,2,3,4
dim(dat)
colnames(dat)
design <- modelMatrix(tar2, ref="pre")

fit <- lmFit(dat, design)
fit <- eBayes(fit)
names(fit)
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# separate channel analysis

# convert targets frame to be channel (vs array oriented)
targets2 <- targetsA2C(tar2)
targets2

# create design matrix
u <- unique(targets2$Target)
u
f <- factor(targets2$Target, levels=u)
design <- model.matrix(~0+f)
colnames(design) <- u
design
corfit <- intraspotCorrelation(dat, design)
fit <- lmscFit(dat, design, correlation=corfit$consensus)

#making contrast between pre/post groups
cont.matrix <- makeContrasts("pre-post",levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, adjust="BH")
names(topTable(fit2, adjust="BH"))


edat <- topTable(fit2, adjust="BH", number = 50000)
dim(dat)
dim(edat)
colnames(edat)
position <- as.numeric(rownames(edat))
position
#rownames <- RGb$genes$GeneName[position]
rname <- as.matrix(RGb$genes$GeneName[position])
dim(rname)
head(rname)
fc <- as.matrix(edat$logFC)
dim(fc)
t <- as.matrix(edat$t)
dim(t)

#fitdat <- cbind(rname, fc)
#dim(fitdat)
#fitdat[,1]
#rownames(fitdat) <- rownames
#rownames(fitdat) <- RGb$genes$GeneName[position]
rownames(fc) <- RGb$genes$GeneName[position]
rownames(t) <- RGb$genes$GeneName[position]
#agg <- apply(fitdat,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
agg <- apply(fc,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
tagg <- apply(t,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
dim(agg)

colnames(fc)

createrankedlist <- function(x,y){
  genename <- as.data.frame(rownames(x))
  names(genename) <- "genename"
  new.dat <- cbind(genename,x)
  filename = y
  write.table(new.dat, file=y, quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")
}
#1 - all samples
#2 - 1861, 1910
#3 - paired only, filter1
#4 - paired only, filter2
#FC
#createrankedlist(agg, "bothcohorts_preranked_log2foldchange.rnk")
#createrankedlist(agg, "bothcohorts_preranked_log2foldchange2.rnk")
#createrankedlist(agg, "bothcohorts_preranked_log2foldchange3.rnk")
#createrankedlist(agg, "bothcohorts_preranked_log2foldchange4.rnk")
createrankedlist(agg, "bothcohorts_preranked_log2foldchange1b.rnk")
#t-statistic
#createrankedlist(tagg, "bothcohorts_preranked_t.rnk")
#createrankedlist(tagg, "bothcohorts_preranked_t2.rnk")
#createrankedlist(tagg, "bothcohorts_preranked_t3.rnk")
#createrankedlist(tagg, "bothcohorts_preranked_t4.rnk")
createrankedlist(tagg, "bothcohorts_preranked_t1b.rnk")

####################################################################################################################################################################
####################################################################################################################################################################
####################################################################################################################################################################

library(ggplot2)
#1,2
# normwithin <- normwithin
tar2 <- targets
normbetween <-normalizeBetweenArrays(normwithin,method='Aquantile')
normwithin <- normwithin[normbetween$genes$ControlType==0,]
normbetween <- normbetween[normbetween$genes$ControlType==0,]

#Convert MA back to RG
#1,2,3,4
RGb <- RG.MA(normbetween)
cy3 = RGb$R
rownames(cy3) <- RGb$genes$GeneName
colnames(cy3) <- paste(tar2$Cy3, tar2$Cy3_Sample, sep="-")
colnames(cy3)
# cy5
cy5 = RGb$G
rownames(cy5) <- RGb$genes$GeneName
colnames(cy5) <- paste(tar2$Cy5, tar2$Cy5_sample, sep="-")
colnames(cy5)

######################################################################################################################################
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
r <- cbind(pppre,pppost)
# just pre
# r <- pppre
colnames(r)
dim(r)

####################################################################################################################################################################
# all samples, lt/mt 2 month therapy
# NOTE: NEED TO INCLUDE ALL SAMPLES, NOT JUST PAIRS
missing <- "pre-PH1816"
x <- which(colnames(r)==missing, arr.ind=TRUE)
colnames(x)
m <- r[, -c(x)]
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
  #print((meta$twomonthstherapy[b]==1)==FALSE)
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

pm <- paste("MORE", colnames(m)[mt], sep="-")
pl <- paste("LESS", colnames(m)[lt], sep="-")
colnames(more) <- pm
colnames(less) <- pl
mo <- more[,c(1,10,2,11,3,12,5,14,6,15,4,7:9,13)]
colnames(less)
le <- less[,c(1,26,2,27,3,28,4,29,14,30,15,31,16,32,17,33,18,34,19,35,23,36,5,6,7,8,9,10,11,12,13,20,21,22,24,25)]
dim(mo)
dim(le)
r <- cbind(le, mo)
colnames(r)
####################################################################################################################################################################
####################################################################################################################################################################
# all samples, lt/mt 12 month overall survival
missing <- "pre-PH1816"
x <- which(colnames(r)==missing, arr.ind=TRUE)
#r[,-c(x)]
m <- r[, -c(x)]
dim(r)
colnames(r)
dim(m)
colnames(m)
lt12mo = NULL
mt12mo = NULL
lt2 = NULL
mt2 = NULL

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
mo2 <- more2[,c(1,9,2,10,4,12,5,13,6,14,3,7,8,11)]
#data.frame(colnames(less2))
#unique(colnames(less2[,c(1,27,2,28,3,29,4,30,5,31,15,32,16,33,17,34,18,35,19,36,23,37,6,7,8,9,10,11,12,13,14,20,21,22,24,25,26)]))
le2 <- less2[,(c(1,27,2,28,3,29,4,30,5,31,15,32,16,33,17,34,18,35,19,36,23,37,6,7,8,9,10,11,12,13,14,20,21,22,24,25,26))]
#le2 <- less[,c(1,26,2,27,3,28,4,29,14,30,15,31,16,32,17,33,18,34,19,35,23,36,5,6,7,8,9,10,11,12,13,20,21,22,24,25)]

dim(mo2)
dim(le2)
r <- cbind(le2, mo2)
colnames(r)
####################################################################################################################################################################
####################################################################################################################################################################
####################################################################################################################################################################
# all samples, lt/mt 9 month overall survival
missing <- "pre-PH1816"
x <- which(colnames(r)==missing, arr.ind=TRUE)
#r[,-c(x)]
m <- r[, -c(x)]
dim(r)
colnames(r)
dim(m)
colnames(m)
lt9mo = NULL
mt9mo = NULL
lt3 = NULL
mt3 = NULL

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
colnames(more3)
mo3 <- more3[,c(1,9,2,10,4,12,5,13,6,14,3,7,8,11)]
le3 <- less3[,(c(1,27,2,28,3,29,4,30,5,31,15,32,16,33,17,34,18,35,19,36,23,37,6,7,8,9,10,11,12,13,14,20,21,22,24,25,26))]
dim(mo3)
dim(le3)
r <- cbind(le3, mo3)
colnames(r)
####################################################################################################################################################################

new.dat <- apply(r,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
dim(new.dat)
colnames(new.dat)

####################################################################################################################################################################
####################################################################################################################################################################
####################################################################################################################################################################

#r <- combat
# Cluster dendogram
rt <- t(r)      #transpose dat
rt.dist <- dist(rt,method="euclidean")  # calculate distance
rt.clust <- hclust(rt.dist,method="single")  # calculate clusters
plot(rt.clust,labels=names(rt),cex=0.75)  # plot cluster tree

# CV v mean pl
dev.off()
r.mean <- apply(log2(r),2,mean)  	# calculate mean for each sample
r.sd <- sqrt(apply(log2(r),2,var))		# calculate st.deviation for each sample
r.cv <- r.sd/r.mean		

plot(r.mean,r.cv,main="colon trial dataset\nSample CV vs. Mean",xlab="Mean",ylab="CV",col='blue',cex=1.5,type="n")
points(r.mean,r.cv,bg="lightblue",col=1,pch=21)
text(r.mean,r.cv,label=dimnames(r)[[2]],pos=1,cex=0.6)
# with sedata, we would eliminate pre-PH1550, pre-PH1811(?), pre-PH1869 as outliers

# Avg correlation plot
dev.off()
r.cor <- cor(r)
r.avg <- apply(r.cor,1,mean)
par(oma=c(3,0.1,0.1,0.1))
plot(c(1,length(r.avg)),range(r.avg),type="n",xlab="",ylab="Avg r",main="Avg correlation of colon trial samples",axes=F)
points(r.avg,bg="red",col=1,pch=21,cex=1.25)
axis(1,at=c(1:length(r.avg)),labels=dimnames(r)[[2]],las=2,cex.lab=0.4,cex.axis=1)
axis(2)
abline(v=seq(0.5,62.5,1),col="grey")

# Pearson's correlation plot
library(gplots)
dev.off()
r.cor <- cor(r)
layout(matrix(c(1,1,1,1,1,1,1,1,2,2), 5, 2, byrow = TRUE))
par(oma=c(5,7,1,1))
cx <- rev(colorpanel(25,"yellow","black","blue"))
leg <- seq(min(r.cor,na.rm=T),max(r.cor,na.rm=T),length=10)
image(r.cor,main="Correlation plot colon biopsy data",axes=F,col=cx)
axis(1,at=seq(0,1,length=ncol(r.cor)),label=dimnames(r.cor)[[2]],cex.axis=0.9,las=2)
axis(2,at=seq(0,1,length=ncol(r.cor)),label=dimnames(r.cor)[[2]],cex.axis=0.9,las=2)
par(mar=rep(2,4))
image(as.matrix(leg),col=cx,axes=T)
tmp <- round(leg,2)
#axis(1,at=seq(0,1,length=length(leg)),labels=tmp,cex.axis=1)

# PCA biplot
# biplot(prcomp(t(log2(r))),cex=0.6,main="PCA biplot - genes and patient samples",col=c("black","grey"),expand=0.8)

dev.off()
datpca <- prcomp(t(r), cor=F)
plot(
  range(datpca$x[, 1]),range(datpca$x[, 2]), 
  type = "n", xlab = "pre", ylab = "post", 
  main = "PCA of colon trial data\npre vs post"
)
points(datpca$x[, 1][tar2$Cy3 == "pre"], datpca$x[, 2][tar2$Cy3 == "pre"], bg = "Red", pch = 21)
text(datpca$x[, 1][tar2$Cy3 == "pre"], datpca$x[, 2][tar2$Cy3 == "pre"], label=names(datpca$x[,1][tar2$Cy3 == "pre"]),pos=1,cex=.75)
points(datpca$x[, 1][tar2$Cy5 == "post"], datpca$x[, 2][tar2$Cy5 == "post"], bg = "Blue", pch=21)
text(datpca$x[, 1][tar2$Cy5 == "post"], datpca$x[, 2][tar2$Cy5 == "post"], label=names(datpca$x[,1][tar2$Cy5 == "post"]),pos=1,cex=.75)
#leg.names <- c("pre", "post")       
#legend(-30,12,leg.names,leg.col,pch=15,cex=.7,horiz=F)
summary(datpca)
dat.pca.var <- round(datpca$sdev^2 / sum(datpca$sdev^2)*100,2)
dev.off()
plot(c(1:length(dat.pca.var)),dat.pca.var,type="b",xlab="# components",ylab="% variance",pch=21,col=1,bg=3,cex=1.5)
title("Scree plot showing % variability explained by each eigenvalue\nKIU/OXF dataset")
#3
# first two components account for 86.5%of the variatino


#USE GENECARD FOR ALTERNATIVE GENENAMES

# # regex for gene name
# dimnames(new.dat)[[2]]
# #TMS1
# length(grep("tms", dimnames(new.dat)[[1]], ignore.case = TRUE))
# rownames(new.dat)[grep("tms", dimnames(new.dat)[[1]], ignore.case = TRUE)]
# #E-cad
# length(grep("cdh1", dimnames(new.dat)[[1]], ignore.case = TRUE))
# rownames(new.dat)[grep("cdh1", dimnames(new.dat)[[1]], ignore.case = TRUE)]
# #TIMP
# length(grep("timp", dimnames(new.dat)[[1]], ignore.case = TRUE))
# rownames(new.dat)[grep("timp", dimnames(new.dat)[[1]], ignore.case = TRUE)]
# #MLH1
# length(grep("mlh", dimnames(new.dat)[[1]], ignore.case = TRUE))
# rownames(new.dat)[grep("mlh", dimnames(new.dat)[[1]], ignore.case = TRUE)]
# #MGMT
# length(grep("mgmt", dimnames(new.dat)[[1]], ignore.case = TRUE))
# rownames(new.dat)[grep("mgmt", dimnames(new.dat)[[1]], ignore.case = TRUE)]
# #p16 AKA CDKN2A
# length(grep("p16", dimnames(new.dat)[[1]], ignore.case = TRUE))
# rownames(new.dat)[grep("p16", dimnames(new.dat)[[1]], ignore.case = TRUE)]
# #PD1
# length(grep("pd1", dimnames(new.dat)[[1]], ignore.case = TRUE))
# rownames(new.dat)[grep("pd1", dimnames(new.dat)[[1]], ignore.case = TRUE)]
# #PDL1 AKA CD274
# length(grep("cd274", dimnames(new.dat)[[1]], ignore.case = TRUE))
# rownames(new.dat)[grep("cd274", dimnames(new.dat)[[1]], ignore.case = TRUE)]
# #PD1
# length(grep("pdcd1", dimnames(new.dat)[[1]], ignore.case = TRUE))
# rownames(new.dat)[grep("pdcd1", dimnames(new.dat)[[1]], ignore.case = TRUE)]
# dim(r)
# length(grep("p16", dimnames(r)[[1]], ignore.case = TRUE))
# rownames(r)[grep("p16", dimnames(r)[[1]], ignore.case = TRUE)]

# as.data.frame(colnames(new.dat))
# #2
# pre_pos = c(1,4,5,7,10,11,13,14,17,19,20)
# #3
# pre_pos = c(1,4,5,7,9,10,13)
# #4
# pre_pos = c(2,3)

####################################################################################################################################################################
####################################################################################################################################################################
####################################################################################################################################################################

#1,2,3,4
new.datfil <- new.dat

as.matrix(colnames(new.dat))
colnames(new.datfil)
dim(new.datfil)
dim(new.dat)


##########################################
# #1 #####################################
##########################################

ATF1
FAM3C
HBD
HSPA13
UBE2B
C16orf87
MRPS30
C14orf142
TSGA14
TCTN3

##########################################
# #2 #####################################
##########################################

vgl <- new.datfil["VGLL1",]
ens <- new.datfil["ENST00000439198",]
tce <- new.datfil["TCEB2",]
lpl <- new.datfil["LPL",]
ppp <- new.datfil["PPP1CC",]
mmp <- new.datfil["MMP27",]
prd <- new.datfil["PRDX2",]
dek <- new.datfil["DEK",]
xm0 <- new.datfil["XM_002342506",]
glp <- new.datfil["GLP1R",]
ap4 <- new.datfil["AP4S1",]
a33 <- new.datfil["A_33_P3298830",]
c15 <- new.datfil["C15orf37",]
cop <- new.datfil["COPB1",]
c11 <- new.datfil["C11orf46",]
c4o <- new.datfil["C4orf33",]
chc <- new.datfil["CHCHD3",]
ens2 <- new.datfil["ENST00000512519",]
dph <- new.datfil["DPH2",]
mtp <- new.datfil["MTPN",]

allgenes <- rbind(vgl, ens, tce, lpl, ppp, mmp, prd, dek, xm0, glp, ap4, a33, c15, cop, c11, c4o, chc, ens2, dph, mtp)
allgenes
rownames(allgenes) <- as.matrix(c("VGLL1", "ENST00000439198", "TCEB2", "LPL", "PPP1CC", "MMP27", "PRDX2", "DEK", "XM_002342506", "GLP1R", "AP4S1", "A_33_P3298830", "C15orf37", "COPB1", "C11orf46" , "C4orf33", "CHCHD3", "ENST00000512519", "DPH2", "MTPN"))
genename <- "all genes"
colors <- colors()[c(130, 8, 461, 47, 84, 53, 554, 418, 640, 275, 143, 258, 345,123, 35,86, 92,29,70,84)]

for(i in 1:nrow(allgenes)){
  colors <- colors()[c(130, 8, 461, 47, 84, 53, 554, 418, 640, 275, 143, 258, 345,123, 35,86, 92,29,70,84)]
  print(i)
  genes <- rbind(allgenes[i,])
  genename <- as.character(rownames(allgenes)[i])
  rownames(genes) <- genename
  color <- colors[i]
  print(color)
  fn = as.character(paste("group2_",genename,".png", sep=""))
  png(filename = fn, width = 1000, height = 1000)
  
  plot(c(1,ncol(r)),range(genes),type='n',main=paste("Profile plot of group2",genename, sep="_"),xlab="",ylab="Expression intensity",axes=F)
  #1
  axis(side=1,at=c(1:ncol(genes)),labels=colnames(genes),cex.axis=0.8,las=2)
  #2
  #axis(side=1,at=c(1:22),labels=colnames(genes),cex.axis=0.4,las=2)
  axis(side=2)
  for(j in 1:nrow(genes)) {
    ycoord <- as.numeric(genes[j,])
    lines(c(1:ncol(genes)),ycoord,col=color,lwd=2)
    points(c(1:ncol(genes)),ycoord,col=color,pch=19, lwd=1)
  }
  #legend for genes
  legend("topright", legend=rownames(genes), fill=color, bg="white", cex = .8, ncol=1)
  dev.off()
}

##########################################
# #3 #####################################
##########################################
fam3 <- new.datfil["FAM3C",]
vgl <- new.datfil["VGLL1",]
soh <- new.datfil["SOHLH1",]
pof <- new.datfil["POFUT1",]
btg <- new.datfil["BTG3",]
kcn <- new.datfil["KCNJ5",]
rn2 <- new.datfil["RN28S1",]
prs <- new.datfil["PRSS45",]
pig <- new.datfil["PIGH",]
fam1 <- new.datfil["FAM170B",]
ens <- new.datfil["ENST00000409517",]
ens2 <- new.datfil["ENST00000340284",]
a33 <- new.datfil["A_33_P3230369",]
loc <- new.datfil["LOC100133224",]
dax <- new.datfil["DAXX",]
c13 <- new.datfil["C13orf31",]
krt <- new.datfil["KRT82",]
fbx <- new.datfil["FBXL21",]
med <- new.datfil["MED27",]
rgs <- new.datfil["RGS13",]


genes <- rbind(fam3, vgl, soh, pof, btg, kcn, rn2, prs, pig, fam1, ens, ens2, a33, loc, dax, c13, krt, fbx, med, rgs)
dim(genes)
allgenes <- genes
genename <- "all genes"
rownames(allgenes) <- as.matrix(c("FAM3C", "VGLL1","SOHLH1","POFUT1","BTG3","KCNJ5","RN28S1","PRSS45","PIGH","FAM170B","ENST00000409517","ENST00000340284","A_33_P3230369","LOC100133224","DAXX","C13orf31","KRT82","FBXL21","MED27","RGS13"))
rownames(allgenes)
colors <- colors()[c(130, 8, 461, 47, 84, 53, 554, 418, 640, 275, 143, 258, 345,123, 35,86, 92,29,63,84)]



for(i in 1:nrow(allgenes)){
  colors <- colors()[c(130, 8, 461, 47, 84, 53, 554, 418, 640, 275, 143, 258, 345,123, 35,86, 92,29,70,84)]
  print(i)
  genes <- rbind(allgenes[i,])
  genename <- as.character(rownames(allgenes)[i])
  rownames(genes) <- genename
  color <- colors[i]
  print(color)
  fn = as.character(paste("group3_",genename,".png", sep=""))
  png(filename = fn, width = 1000, height = 1000)
  
  plot(c(1,ncol(r)),range(genes),type='n',main=paste("Profile plot of group3",genename, sep="_"),xlab="",ylab="Expression intensity",axes=F)
  #1
  axis(side=1,at=c(1:ncol(genes)),labels=colnames(genes),cex.axis=0.8,las=2)
  #2
  #axis(side=1,at=c(1:22),labels=colnames(genes),cex.axis=0.4,las=2)
  axis(side=2)
  for(j in 1:nrow(genes)) {
    ycoord <- as.numeric(genes[j,])
    lines(c(1:ncol(genes)),ycoord,col=color,lwd=2)
    points(c(1:ncol(genes)),ycoord,col=color,pch=19, lwd=1)
  }
  #legend for genes
  legend("topright", legend=rownames(genes), fill=color, bg="white", cex = .8, ncol=1)
  dev.off()
}

##########################################
# #4 #####################################
##########################################
one <- new.datfil["POLE3",]
two <- new.datfil["MAGOHB",]
three <- new.datfil["AY927536",]
four <- new.datfil["GRIPAP1",]
five <- new.datfil["LOC100288842",]
six <- new.datfil["EI24",]
seven <- new.datfil["RBM18",]
eight <- new.datfil["CKS1B",]
nine <- new.datfil["ENST00000340284",]
ten <- new.datfil["MMGT1",]
eleven <- new.datfil["YWHAQ",]
twelve <- new.datfil["WDSUB1",]
thirteen <- new.datfil["A_33_P3351615",]
fourteen <- new.datfil["C10orf88",]
fifteen <- new.datfil["SNRPE",]
sixteen <- new.datfil["BRD7",]
seventeen <- new.datfil["FAM3C",]
eighteen <- new.datfil["GNAQ",]
nineteen <- new.datfil["C11orf73",]
twenty <- new.datfil["FAM170B",]
twentyone <- new.datfil["SNRNP27",]
twentytwo <- new.datfil["E2F5",]
twentythree <- new.datfil["RNF40",]
twentyfour <- new.datfil["DCTN5",]
twentyfive <- new.datfil["OR10H5",]
twentysix <- new.datfil["C2orf76",]
twentyseven <- new.datfil["KPNA4",]
twentyeight <- new.datfil["SRFBP1",]
twentynine <- new.datfil["PRDX3",]
thirty <- new.datfil["HIPK3",]
thirtyone <- new.datfil["ELAC1",]
thirtytwo <- new.datfil["RAD1",]
thirtythree <- new.datfil["FCHSD2",]
thirtyfour <- new.datfil["FNBP1L",]
thirtyfive <- new.datfil["COMMD10",]
thirtysix <- new.datfil["C16orf87",]
thirtyseven <- new.datfil["A_33_P3354574",]
thirtyeight <- new.datfil["CHAC2",]
thirtynine <- new.datfil["C2orf69",]
fourty <- new.datfil["ENST00000414544",]
fourtyone <- new.datfil["KLF11",]
fourtytwo <- new.datfil["TAF1A",]
fourtythree <- new.datfil["IQSEC3",]
fourtyfour <- new.datfil["POGZ",]
fourtyfive <- new.datfil["MFN1",]
fourtysix <- new.datfil["LOC390595",]
fourtyseven <- new.datfil["ATPBD4",]
fourtyeight <- new.datfil["C3orf26",]
fourtynine <- new.datfil["UAP1",]
fifty <- new.datfil["SOHLH1",]
fiftyone <- new.datfil["SHKBP1",]
fiftytwo <- new.datfil["ZNF230",]
fiftythree <- new.datfil["ATP6V0E2",]
fiftyfour <- new.datfil["C16orf80",]
fiftyfive <- new.datfil["RPS4XP21",]
fiftysix <- new.datfil["ZMAT5",]
fiftyseven <- new.datfil["PIGH",]
fiftyeight <- new.datfil["XM_002342506",]
fiftynine <- new.datfil["FAM161A",]
sixty <- new.datfil["A_33_P3370612",]
sixtyone <- new.datfil["GLRX3",]
sixtytwo <- new.datfil["RPRM",]
sixtythree <- new.datfil["CACYBP",]
sixtyfour <- new.datfil["PGM3",]
sixtyfive <- new.datfil["NAA50",]
sixtysix <- new.datfil["RCAN1",]
sixtyseven <- new.datfil["MED27",]
sixtyeight <- new.datfil["ENST00000391369",]
sixtynine <- new.datfil["DNAJB9",]
seventy <- new.datfil["HACE1",]
seventyone <- new.datfil["A_24_P169843",]
seventytwo <- new.datfil["ABCF2",]
seventythree <- new.datfil["FAM91A1",]
seventyfour <- new.datfil["ANKS3",]
seventyfive <- new.datfil["RBM7",]
seventysix <- new.datfil["C2orf29",]
seventyeight <- new.datfil["GPN1",]
seventynine <- new.datfil["GPN3",]
eighty <- new.datfil["ATG3",]
eightyone <- new.datfil["SDCBP",]
eightytwo <- new.datfil["LOC100131101",]
eightythree <- new.datfil["POFUT1",]
eightyfour <- new.datfil["L2HGDH",]
eightyfive <- new.datfil["MCM2",]
eightysix <- new.datfil["PAIP2",]
eightyseven <- new.datfil["WDR92",]
eightyeight <- new.datfil["EED",]
eightynine <- new.datfil["TSG1",]
ninety <- new.datfil["PTS",]
ninetyone <- new.datfil["TMEM126B",]
ninetytwo <- new.datfil["SYF2",]
ninetythree <- new.datfil["A_24_P67408",]
ninetyfour <- new.datfil["HARBI1",]
ninetyfive <- new.datfil["NUP93",]
ninetysix <- new.datfil["CFDP1",]
ninetyseven <- new.datfil["AA627135",]
ninetyeight <- new.datfil["LOC729291",]
ninetynine <- new.datfil["LOC643802",]
hundred <- new.datfil["DDX21",]
hundredone <- new.datfil["VSIG8",]
hundredtwo <- new.datfil["BTG3",]
hundredthree <- new.datfil["ZBTB40",]
hundredfour <- new.datfil["TAGLN2",]
hundredfive <- new.datfil["TSPAN13",]
hundredsix <- new.datfil["ENST00000434635",]
hundredseven <- new.datfil["C10orf84",]
hundredeight <- new.datfil["METAP2",]
hundrednine <- new.datfil["CSTF2T",]
oneone <- new.datfil["MRPL19",]
onetwo <- new.datfil["A_33_P3377714",]
onethree <- new.datfil["RPF2",]
onefour <- new.datfil["HBB",]
onefive <- new.datfil["SKP1",]
onesix <- new.datfil["LOC100129195",]
oneseven <- new.datfil["INTS3",]
oneeight <- new.datfil["VGLL1",]
onenine <- new.datfil["CAPZA2",]
twoone <- new.datfil["PSMC3",]
twotwo <- new.datfil["PIGM",]
twothree <- new.datfil["LRRC8B",]
twofour <- new.datfil["POLR2A",]
twofive <- new.datfil["HPSE",]
twosix <- new.datfil["VSTM1",]
twoseven <- new.datfil["OGFOD1",]
twoeight <- new.datfil["C13orf34",]
twonine <- new.datfil["ECT2",]
threeone <- new.datfil["C9orf80",]
threetwo <- new.datfil["CIRH1A",]
threethree <- new.datfil["C2orf60",]
threefour <- new.datfil["LARP4",]
threefive <- new.datfil["PCNP",]
threesix <- new.datfil["ZNF828",]
threeseven <- new.datfil["ENST00000356822",]
threeeight <- new.datfil["C11orf74",]

genes <- rbind(one,two,three,four,five,six,seven,eight,nine,ten,eleven,twelve,thirteen,fourteen,fifteen,sixteen,seventeen,eighteen,nineteen,twenty,twentyone,twentytwo,twentythree,twentyfour,twentyfive,twentysix,twentyseven,twentyeight,twentynine,thirty,thirtyone,thirtytwo,thirtythree,thirtyfour,thirtyfive,thirtysix,thirtyseven,thirtyeight,thirtynine,fourty,fourtyone,fourtytwo,fourtythree,fourtyfour,fourtyfive,fourtysix,fourtyseven,fourtyeight,fourtynine,fifty,fiftyone,fiftytwo,fiftythree,fiftyfour,fiftyfive,fiftysix,fiftyseven,fiftyeight,fiftynine,sixty,sixtyone,sixtytwo,sixtythree,sixtyfour,sixtyfive,sixtysix,sixtyseven,sixtyeight,sixtynine,seventy,seventyone,seventytwo,seventythree,seventyfour,seventyfive,seventysix,seventyeight,seventynine,eighty,eightyone,eightytwo,eightythree,eightyfour,eightyfive,eightysix,eightyseven,eightyeight,eightynine,ninety,ninetyone,ninetytwo,ninetythree,ninetyfour,ninetyfive,ninetysix,ninetyseven,ninetyeight,ninetynine,hundred,hundredone,hundredtwo,hundredthree,hundredfour,hundredfive,hundredsix,hundredseven,hundredeight,hundrednine,oneone,onetwo,onethree,onefour,onefive,onesix,oneseven,oneeight,onenine,twoone,twotwo,twothree,twofour,twofive,twosix,twoseven,twoeight,twonine,threeone,threetwo,threethree,threefour,threefive,threesix,threeseven,threeeight)
dim(genes)
allgenes <- genes
rownames(allgenes) <- as.matrix(c("POLE3","MAGOHB","AY927536","GRIPAP1","LOC100288842","EI24","RBM18","CKS1B","ENST00000340284","MMGT1","YWHAQ","WDSUB1","A_33_P3351615","C10orf88","SNRPE","BRD7","FAM3C","GNAQ","C11orf73","FAM170B","SNRNP27","E2F5","RNF40","DCTN5","OR10H5","C2orf76","KPNA4","SRFBP1","PRDX3","HIPK3","ELAC1","RAD1","FCHSD2","FNBP1L","COMMD10","C16orf87","A_33_P3354574","CHAC2","C2orf69","ENST00000414544","KLF11","TAF1A","IQSEC3","POGZ","MFN1","LOC390595","ATPBD4","C3orf26","UAP1","SOHLH1","SHKBP1","ZNF230","ATP6V0E2","C16orf80","RPS4XP21","ZMAT5","PIGH","XM_002342506","FAM161A","A_33_P3370612","GLRX3","RPRM","CACYBP","PGM3","NAA50","RCAN1","MED27","ENST00000391369","DNAJB9","HACE1","A_24_P169843","ABCF2","FAM91A1","ANKS3","RBM7","C2orf29","GPN1","GPN3","ATG3","SDCBP","LOC100131101","POFUT1","L2HGDH","MCM2","PAIP2","WDR92","EED","TGS1","PTS","TMEM126B","SYF2","A_24_P67408","HARBI1","NUP93","CFDP1","AA627135","LOC729291","LOC643802","DDX21","VSIG8","BTG3","ZBTB40","TAGLN2","TSPAN13","ENST00000434635","C10orf84","METAP2","CSTF2T","MRPL19","A_33_P3377714","RPF2","HBB","SKP1","LOC100129195","INTS3","VGLL1","CAPZA2","PSMC3","PIGM","LRRC8B","POLR2A","HPSE","VSTM1","OGFOD1","C13orf34","ECT2","C9orf80","CIRH1A","C2orf60","LARP4","PCNP","ZNF828","ENST00000356822","C11orf74"))
dim(allgenes)

genename <- "all genes"
as.data.frame(rownames(allgenes))
colors <- colors()[c(130, 8, 461, 47, 84, 53, 554, 418, 640, 275, 143, 258, 345,123, 35,86, 92,29,63,84, 130, 8, 461, 47, 84, 53, 554, 418, 640, 275, 143, 258, 345,123, 35,86, 92,29,63,84, 130, 8, 461, 47, 84, 53, 554, 418, 640, 275, 143, 258, 345,123, 35,86, 92,29,63,84, 130, 8, 461, 47, 84, 53, 554, 418, 640, 275, 143, 258, 345,123, 35,86, 92,29,63,84, 130, 8, 461, 47, 84, 53, 554, 418, 640, 275, 143, 258, 345,123, 35,86, 92,29,63,84, 130, 8, 461, 47, 84, 53, 554, 418, 640, 275, 143, 258, 345,123, 35,86, 92,29,63,84, 130, 8, 461, 47, 84, 53, 554, 418, 640, 275, 143, 258, 345,123)]

for(i in 1:nrow(allgenes)){
  colors <- colors()[c(130, 8, 461, 47, 84, 53, 554, 418, 640, 275, 143, 258, 345,123, 35,86, 92,29,63,84, 130, 8, 461, 47, 84, 53, 554, 418, 640, 275, 143, 258, 345,123, 35,86, 92,29,63,84, 130, 8, 461, 47, 84, 53, 554, 418, 640, 275, 143, 258, 345,123, 35,86, 92,29,63,84, 130, 8, 461, 47, 84, 53, 554, 418, 640, 275, 143, 258, 345,123, 35,86, 92,29,63,84, 130, 8, 461, 47, 84, 53, 554, 418, 640, 275, 143, 258, 345,123, 35,86, 92,29,63,84, 130, 8, 461, 47, 84, 53, 554, 418, 640, 275, 143, 258, 345,123, 35,86, 92,29,63,84, 130, 8, 461, 47, 84, 53, 554, 418, 640, 275, 143, 258, 345,123)]
  genes <- rbind(allgenes[i,])
  genename <- as.character(rownames(allgenes)[i])
  rownames(genes) <- genename
  color <- colors[134]
  fn = as.character(paste("group4_",genename,".png", sep=""))
  png(filename = fn, width = 1000, height = 1000)
  
  plot(c(1,ncol(r)),range(genes),type='n',main=paste("Profile plot of group4",genename, sep="_"),xlab="",ylab="Expression intensity",axes=F)
  axis(side=1,at=c(1:ncol(genes)),labels=colnames(genes),cex.axis=0.8,las=2)
  axis(side=2)
  for(j in 1:nrow(genes)) {
    ycoord <- as.numeric(genes[j,])
    lines(c(1:ncol(genes)),ycoord,col=color,lwd=2)
    points(c(1:ncol(genes)),ycoord,col=color,pch=19, lwd=1)
  }
  #legend for genes
  legend("topright", legend=rownames(genes), fill=color, bg="white", cex = .8, ncol=1)
  dev.off()
}

##########################################
# #5 #####################################
##########################################
mpg <- new.datfil["MPG",]
mtd <- new.datfil["MTDH",]
eif <- new.datfil["EIF1",]
c1q <- new.datfil["C1QBP",]
ppp <- new.datfil["PPP1CC",]
sna <- new.datfil["SNAPIN",]
pol <- new.datfil["POLR2K",]
tsg <- new.datfil["TSGA14",]
dek <- new.datfil["DEK",]
cyb <- new.datfil["CYB5B",]
ier <- new.datfil["IER3IP1",]
cni<- new.datfil["CNIH",]
mtp <- new.datfil["MTPN",]
stx <- new.datfil["STXBP3",]
lpl <- new.datfil["LPL",]
tce <- new.datfil["TCEB2",]
mcm <- new.datfil["MCM2",]
grn <- new.datfil["GRN",]
c14 <- new.datfil["C14orf142",]
eed <- new.datfil["EED",]
tct <- new.datfil["TCTN3",]
hmm <- new.datfil["HMMR",]
prd <- new.datfil["PRDX2",]
hbx <- new.datfil["HBXIP",]
a24 <- new.datfil["A_24_P273245",]
cop <- new.datfil["COPB1",]
ube <- new.datfil["UBE2Q2",]
prd3 <- new.datfil["PRDX3",]
chc <- new.datfil["CHCHD3",]
dph <- new.datfil["DPH2",]
usp <- new.datfil["USP1",]
mrp <- new.datfil["MRPL39",]
ple <- new.datfil["PLEKHG4",]
pja <- new.datfil["PJA1",]
nae <- new.datfil["NAE1",]


genes <- rbind(mpg,mtd,eif, c1q, ppp, sna, pol, tsg, dek, cyb, ier, cni, mtp, stx, lpl, tce, mcm, grn, c14, eed, tct, hmm, prd, hbx, a24, cop, ube, prd3, chc, dph, usp, mrp, ple, pja, nae)
dim(genes)
allgenes <- genes
genename <- "all genes"
rownames(allgenes) <- as.matrix(c("MPG","MTDH","EIF1","C1QBP","PPP1CC","SNAPIN","POLR2K","TSGA14","DEK","CYB5B","IER3IP1","CNIH","MTPN","STXBP3","LPL","TCEB2","MCM2","GRN","C14orf142","EED","TCTN3","HMMR","PRDX2","HBXIP","A_24_P273245","COPB1","UBE2Q2","PRDX3","CHCHD3","DPH2","USP1","MRPL39","PLEKHG4","PJA1","NAE1"))
rownames(allgenes)
colors <- colors()[c(130, 8, 461, 47, 84, 53, 554, 418, 640, 275, 143, 258, 345,123, 35,86, 92,29,70,84, 130, 8, 461, 47, 84, 53, 554, 418, 640, 275, 143, 258, 345,123, 35)]
dim(allgenes)
for(i in 1:nrow(allgenes)){
  colors <- colors()[c(130, 8, 461, 47, 84, 53, 554, 418, 640, 275, 143, 258, 345,123, 35,86, 92,29,70,84, 130, 8, 461, 47, 84, 53, 554, 418, 640, 275, 143, 258, 345,123, 35)]
  print(i)
  genes <- rbind(allgenes[i,])
  genename <- as.character(rownames(allgenes)[i])
  rownames(genes) <- genename
  color <- colors[i]
  print(color)
  fn = as.character(paste("group5_",genename,".png", sep=""))
  png(filename = fn, width = 1000, height = 1000)
  
  plot(c(1,ncol(r)),range(genes),type='n',main=paste("Profile plot of group5",genename, sep="_"),xlab="",ylab="Expression intensity",axes=F)
  #1
  axis(side=1,at=c(1:ncol(genes)),labels=colnames(genes),cex.axis=0.8,las=2)
  #2
  #axis(side=1,at=c(1:22),labels=colnames(genes),cex.axis=0.4,las=2)
  axis(side=2)
  for(j in 1:nrow(genes)) {
    ycoord <- as.numeric(genes[j,])
    lines(c(1:ncol(genes)),ycoord,col=color,lwd=2)
    points(c(1:ncol(genes)),ycoord,col=color,pch=19, lwd=1)
  }
  #legend for genes
  legend("topright", legend=rownames(genes), fill=color, bg="white", cex = .8, ncol=1)
  dev.off()
}


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
# rnewgenes <- newgenes[,-c(1,2)]
# colnames(rnewgenes)
plotgenes <- newgenes #or rnewgenes
title <- "all samples, AIM genes"
title <- "all pairs, AIM genes"
colbreak = c(seq(-2,-.25,length=100),seq(-.25,.25,length=100),seq(.25,2,length=100))

plotgenes <- rnewgenes
title <- "all pairs, without 1900"
colbreak = c(seq(-2,-.5,length=100),seq(-.5,.5,length=100),seq(.5,2,length=100))

mycolors <- colorRampPalette(c("red", "black", "green"))(n = 299)
par(cex.main=.8)

# heatmap.2(as.matrix(plotgenes), main=title, col=mycolors, breaks=colbreak, density.info="none", trace="none", dendrogram=c("row"), symm=F,symkey=F,symbreaks=T, scale="row", Colv = NA)
#heatmap.2(plotgenes, main=title, col=mycolors, breaks=colbreak, density.info="none", trace="none", dendrogram=c("row"), symm=F,symkey=F,symbreaks=T, scale="row", Colv = NA)
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
################################################################################################
#all samples
new.datfil <- new.dat[,c(1,17,2,18,3,19,4,20,5,21,6,22,7,23,8,24,9,25,10,26,11,27,12,28,13,29,14,30,15,31,16,32)]
f <- factor(as.character(sub("-\\w+-\\w+$", "", colnames(new.datfil))))

################################################################################################
#pairs only
new.datfil <- new.dat[,c(1,17,2,18,3,19,4,20,5,21,6,22,7,23,8,24,9,25,10,26,11,27,12,28,13,29,14,30,15,31,16,32)]
f <- factor(as.character(sub("-\\w+$", "", colnames(new.datfil))))
f
################################################################################################

dim(new.datfil)
colnames(new.datfil)


library(gplots)
f
design <- model.matrix(~f)
design
fit <- eBayes(lmFit(new.datfil, design))
selected  <- fit$p.value[, 2] <.05 # or p<.01 #p.adjust(fit$p.value[, 2], method="BY")<.05
esetSel <- new.datfil[selected,]
dim(esetSel)
# title <- "all samples, more/less 2mo. therapy, limma-filtered genes, p<.05"
title <- "all pairs, limma-filtered genes, p<.05"
colbreak = c(seq(-2,-.5,length=100),seq(-.5,.5,length=100),seq(.5,2,length=100))
mycolors <- colorRampPalette(c("red", "black", "green"))(n = 299)
par(cex.main=.8)
ht <- heatmap.2(esetSel, main=title, col=mycolors, breaks=colbreak, density.info="none", trace="none", srtCol = 45, margins=c(7,10), cexCol= .7, offsetCol = -.6, dendrogram=c("none"), symm=F,symkey=F,symbreaks=T, scale="row", Colv = NA)
#text(x = seq_along(labels), y=-.2, labels=colnames(esetSel), srt=45, xpd=TRUE, cex=.8)
clust <- hclust(dist(esetSel))
plot(clust, cex = .5)
prop <- propexpr(esetSel)
#clustering 
plotMDS(esetSel, cex=.7, labels = factor(as.character(sub("-\\w+$", "", colnames(new.datfil)))))#, col=c(rep("black",16), rep("red",16)), labels= c(rep("pre",16), rep("post",16)))
plotMDS(esetSel, cex=.7, col=c(rep("black",16), rep("red",16)), labels= c(rep("pre",16), rep("post",16)), xlab=NA, ylab=NA, )

corfit <- duplicateCorrelation(new.datfil, design)
# fittr <- lmFit(MA2, design, block = biolrep, cor = corfit$consensus) 
fittr <- lmFit(new.datfil, design, cor = corfit$consensus) 
fittr <- eBayes(fittr) 
tabtr<-toptable(fittr, number=nrow(new.datfil), sort.by="none") 
head(tabtr)
cols<-rep("black", nrow(tabtr))
cols[tabtr$logFC<=-1 & tabtr$adj.P.Val<=0.05]<-"red" 
cols[tabtr$logFC>=1 & tabtr$adj.P.Val<=0.05]<-"green" 
plot(x=tabtr$logFC, y=-log(tabtr$adj.P.Val), pch=16, cex=0.5, col=cols)
cols[tabtr$t<=-2 & tabtr$adj.P.Val<=0.05]<-"red" 
cols[tabtr$t>=2 & tabtr$adj.P.Val<=0.05]<-"green" 
plot(x=tabtr$t, y=tabtr$adj.P.Val, pch=16, cex=0.5, col=cols)

colnames(tabtr)
length(tabtr$logFC)
length(tabtr$adj.P.Val)
# sd <- 0.3*sqrt(4/rchisq(1000,df=4))
# x <- matrix(rnorm(1000*6,sd=sd),1000,6)
# rownames(x) <- paste("Gene",1:1000)
# x[1:50,4:6] <- x[1:50,4:6] + 2
# mds <- plotMDS(x,  col=c(rep("black",3), rep("red",3)) )
# plotMDS(mds,  col=c(rep("black",3), rep("red",3)), labels= c(rep("Grp1",3), rep("Grp2",3)))


for(i in 1:10){gc()}
newgenes <- new.datfil
# ALL PRES
# rnewgenes <- newgenes[,-c(1,2)]
# colnames(rnewgenes)
plotgenes <- newgenes #or rnewgenes
title <- "all pairs"
colbreak = c(seq(-2,-.25,length=100),seq(-.25,.25,length=100),seq(.25,2,length=100))

plotgenes <- rnewgenes
title <- "all pairs, without 1900"
colbreak = c(seq(-2,-.5,length=100),seq(-.5,.5,length=100),seq(.5,2,length=100))

mycolors <- colorRampPalette(c("red", "black", "green"))(n = 299)
par(cex.main=.8)

# heatmap.2(as.matrix(plotgenes), main=title, col=mycolors, breaks=colbreak, density.info="none", trace="none", dendrogram=c("row"), symm=F,symkey=F,symbreaks=T, scale="row", Colv = NA)
heatmap.2(plotgenes, main=title, col=mycolors, breaks=colbreak, density.info="none", trace="none", srtCol = 45, margins=c(7,10), cexCol= .7, offsetCol = -.6, dendrogram=c("row"), symm=F,symkey=F,symbreaks=T, scale="row", Colv = NA)

##########################################
# HC #####################################
##########################################

ecad <- new.datfil["CDH1",]
p16 <- new.datfil["CDKN2A",]
tms1 <- new.datfil["PYCARD",]
timp1 <- new.datfil["TIMP1",]
timp2 <- new.datfil["TIMP2",]
timp3 <- new.datfil["TIMP3",]
timp4 <- new.datfil["TIMP4",]
mlh1 <- new.datfil["MLH1",]
mgmt <- new.datfil["MGMT",]
pd1 <- new.datfil["PDCD1",]
pdl1 <- new.datfil["CD274",]
pdl2 <- new.datfil["PDCD1LG2",]
sfrp1 <- new.datfil["SFRP1",]
tfpi2 <- new.datfil["TFPI2",]
gata4 <- new.datfil["GATA4",]
gata5 <- new.datfil["GATA5",]
apc <- new.datfil["APC",]
chfr <- new.datfil["CHFR",]
rassf1 <- new.datfil["RASSF1",]
hin1 <- new.datfil["SCGB3A1",]

##########################################################################################################################################################################################
# CATEGORIZE VARIABLE
# di_matrix <- matrix(round(runif(r * c, 1, 100)), ncol = c)
# #get median per gene
# genes_median <- apply(di_matrix, 2, median)
# #convert to 0 and 1
# di_matrix01 <- ifelse(di_matrix > genes_median, 1, 0)
# di_matrix
# genes_median
# t(di_matrix01)
##########################################################################################################################################################################################
# for survival analysis
tempdat <- rbind(ecad, p16, tms1, timp1, timp2, timp3, timp4, mlh1, mgmt, pd1, pdl1, pdl2, sfrp1, tfpi2, gata4, gata5, apc, chfr, rassf1, hin1)
dim(tempdat)
tempdat
#get median per gene
genes_median <- apply(tempdat, 1, median)
genes_median
tempdat2 <- ifelse(tempdat > genes_median, 1, 0)
tempdat2
t(tempdat2)
#write.table(t(tempdat2), file="/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/m084b trial/patientstatus_genes.csv", sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)
write.table(t(tempdat2), file="/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/m084b trial/patientstatus_allgenes_pairsonly.csv", sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)
##########################################################################################################################################################################################
colnames(new.datfil)
dim(new.datfil)
genes <-rbind(ecad, p16, tms1, timp1, timp2, timp3, timp4, mlh1, mgmt, pd1, pdl1, pdl2, sfrp1, tfpi2, gata4, gata5, apc, chfr, rassf1, hin1)
allgenes <- genes
genename <- "all genes"
colors <- colors()[c(130, 8, 461, 47, 84, 53, 554, 418, 640, 275, 143, 258)]
#rownames(allgenes) <- as.matrix(c("ECAD", "P16", "TMS1", "TIMP1", "TIMP2", "TIMP3", "TIMP4", "MLH1", "MGMT", "PD1", "PDL1", "PDL2"))
rownames(allgenes)

for(i in 1:nrow(allgenes)){
  colors <- colors()[c(130, 8, 461, 47, 84, 53, 554, 418, 640, 275, 143, 258, 345,123, 35,86, 92,29,70,84)]
  print(i)
  genes <- rbind(allgenes[i,])
  genename <- as.character(rownames(allgenes)[i])
  rownames(genes) <- genename
  color <- colors[i]
  print(color)
  fn = as.character(paste("group5_",genename,".png", sep=""))
  png(filename = fn, width = 1000, height = 1000)
  
  plot(c(1,ncol(r)),range(genes),type='n',main=paste("Profile plot of group5",rownames(allgenes)[i], sep="_"),xlab="",ylab="Expression intensity",axes=F)
  #1
  axis(side=1,at=c(1:ncol(genes)),labels=colnames(genes),cex.axis=0.8,las=2)
  #2
  #axis(side=1,at=c(1:22),labels=colnames(genes),cex.axis=0.4,las=2)
  axis(side=2)
  for(j in 1:nrow(genes)) {
    ycoord <- as.numeric(genes[j,])
    lines(c(1:ncol(genes)),ycoord,col=color,lwd=2)
    points(c(1:ncol(genes)),ycoord,col=color,pch=19, lwd=1)
  }
  #legend for genes
  legend("topright", legend=rownames(genes), fill=color, bg="white", cex = .8, ncol=1)
  dev.off()
}

#head(table(new.dat))
# counts <- table(mtcars$vs, mtcars$gear)
# barplot(counts, main="Car Distribution by Gears and VS",
#         xlab="Number of Gears", col=c("darkblue","red"),
#         legend = rownames(counts), beside=TRUE)


dev.off()


dim(genes)
#>>>>>>
##################################################################################################################################################################################################################################################################################################################################
newgenes = genes[,c(1,17,2,18,3,19,4,20,5,21,6,22,7,23,8,24,9,25,10,26,11,27,12,28,13,29,14,30,15,31,16,32)]
tran <- as.matrix(t(newgenes))
tran
rownames(newgenes)[1]

# BARPLOT
for(i in 1:ncol(tran)){
  g <- tran[,i]
  genename <- rownames(newgenes)[i]
  fn = as.character(paste("/home/steve/.gvfs//onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/m084b trial/aimgenes/","hc_",genename,"_r.png", sep=""))
  png(filename = fn, width = 1000, height = 1000)
  bp <- barplot(g, xlab = "", main=paste("Profile plot of hc",rownames(newgenes)[i], sep="_"), ylab = "Expression intensity", axisnames=F, col=c(1,1,32,32))#rainbow(2))
  axis(1, at=bp, labels=rownames(tran), cex.axis=.7, las=2)
  axis(side=2)
}

##################################################################################################################################################################################################################################################################################################################################
dev.off()

#>>>>>>
plot(c(1,ncol(r)),range(genes),type='n',main=paste("Profile plot of",genename, sep=" "),xlab="Samples",ylab="Expression intensity",axes=F)
#1
axis(side=1,at=c(1:32),labels=colnames(genes),cex.axis=0.4,las=2)
#2
#axis(side=1,at=c(1:22),labels=colnames(genes),cex.axis=0.4,las=2)
axis(side=2)
for(i in 1:length(genes)) {
  r.y <- as.numeric(genes[i,])
  lines(c(1:ncol(genes)),r.y,col=colors[i],lwd=2)
}
#legend for genes
legend("topright", legend=rownames(genes), fill=colors, bg="white", cex = .7, ncol=2)
colors <- colors()[c(130, 8, 461, 47, 84, 53, 554, 418, 640, 275, 143, 258)]

# Hierarchical clustering - Pearson's correlation
#dat.cor <- cor(sedata, use = "pairwise.complete.obs", method = "pearson")
#par(mar = rep(4, 4))
#plot(dat.cor)
#text(dat.cor,label=dimnames(dat.cor)[[2]],pos=1,cex=0.5)

#hierarchical clustering plot
#dat.t <- t(sedata3)
# get pairwise distances
#dat.dist <- dist(dat.t)
# calculate and plot hierarchical tree
#par(mar = rep(4, 4))
#plot(hclust(dat.dist), labels=dimnames(sedata3)[[2]],main="HC Dendogram of Pre and Post Samples",cex.main=0.75)
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#1
ATF1
FAM3C
HBD
HSPA13
UBE2B
C16orf87
MRPS30
C14orf142
TSGA14
TCTN3
PIGH
LYRM1
A_24_P273245
INTS12
STEAP1
PLS1
NETO2
CCDC91
HBB
TIPIN
XRCC4
KLF3
MRPL39
ACN9
GPN3
#2
VGLL1
ENST00000439198
TCEB2
LPL
PPP1CC
MMP27
PRDX2
DEK
XM_002342506
GLP1R
AP4S1
A_33_P3298830
C15orf37
COPB1
C11orf46
C4orf33
CHCHD3
ENST00000512519
DPH2
MTPN
#3
FAM3C
VGLL1
SOHLH1
POFUT1
BTG3
KCNJ5
RN28S1
PRSS45
PIGH
FAM170B
ENST00000409517
ENST00000340284
A_33_P3230369
LOC100133224
DAXX
C13orf31
KRT82
FBXL21
MED27
RGS13
#4
POLE3
MAGOHB
AY927536
GRIPAP1
LOC100288842
EI24
RBM18
CKS1B
ENST00000340284
MMGT1
YWHAQ
A_33_P3351615
C10orf88
SNRPE
BRD7
FAM3C
GNAQ
C11orf73
FAM170B
SNRNP27
E2F5
RNF40
DCTN5
OR10H5
C2orf76
KPNA4
SRFBP1
PRDX3
HIPK3
ELAC1
RAD1
FCHSD2
FNBP1L
COMMD10
C16orf87
A_33_P3354574
CHAC2
C2orf69
ENST00000414544
KLF11
TAF1A
IQSEC3
POGZ
MFN1
LOC390595
ATPBD4
C3orf26
UAP1
SOHLH1
SHKBP1
ZNF230
ATP6V0E2
C16orf80
RPS4XP21
ZMAT5
PIGH
XM_002342506
FAM161A
A_33_P3370612
GLRX3
RPRM
CACYBP
PGM3
NAA50
RCAN1
MED27
ENST00000391369
DNAJB9
HACE1
A_24_P169843
ABCF2
FAM91A1
ANKS3
RBM7
C2orf29
GPN1
GPN3
ATG3
SDCBP
LOC100131101
POFUT1
L2HGDH
MCM2
PAIP2
WDR92
EED
TGS1
PTS
TMEM126B
SYF2
A_24_P67408
HARBI1
NUP93
CFDP1
AA627135
LOC729291
LOC643802
DDX21
VSIG8
BTG3
ZBTB40
TAGLN2
TSPAN13
ENST00000434635
C10orf84
METAP2
CSTF2T
MRPL19
A_33_P3377714
RPF2
HBB
SKP1
LOC100129195
INTS3
VGLL1
CAPZA2
PSMC3
PIGM
LRRC8B
POLR2A
HPSE
VSTM1
OGFOD1
C13orf34
ECT2
C9orf80
CIRH1A
C2orf60
LARP4
PCNP
ZNF828
ENST00000356822
C11orf74
WDSUB1
#5
MPG
MTDH
EIF1
C1QBP
PPP1CC
SNAPIN
POLR2K
TSGA14
DEK
CYB5B
IER3IP1
CNIH
MTPN
STXBP3
LPL
TCEB2
MCM2
GRN
C14orf142
EED
TCTN3
HMMR
PRDX2
HBXIP
A_24_P273245
COPB1
UBE2Q2
PRDX3
CHCHD3
DPH2
USP1
MRPL39
PLEKHG4
PJA1
NAE1
#AIM genes
ADM
ADORA2A
ADORA2B
ALOX5AP
ANKRD1
ANXA1
AOC3
AOX1
B2M
BCL2
C4BPB
CAMP
CCL2
CCL20
CCL26
CCL28
CCL3
CCL3L3
CCL4
CCL5
CCR9
CD44
CD81
CRP
CSF2
CTGF
CTSS
CXCL1
CXCL11
CXCL12
CXCL2
CXCL3
CXCL6
CXCL9
CXCR4
CXCR7
DCBLD2
DDX58
DEFB103A
EGR1
EIF2AK2
EIF4E
EIF4G1
EREG
FOS
FOSL1
GBP1
GBP5
GP9
HCP5
HERC5
HLA-B
HLA-C
HLA-DMB
HLA-DRB1
HSP90AA1
ICAM1
IFI27
IFI6
IFIH1
IFIT1
IFIT2
IFIT3
IFITM1
IFNGR1
IL1A
IL1B
IL1R2
IL32
IL6
IL6ST
IL8
INHBA
IRF6
IRF7
IRF8
IRF9
ISG15
ISG20
ITGAV
JAK2
KCNN4
KLK8
KLRC2
LBP
LCK
LSP1
LY96
LYN
LYST
MDK
MRC2
MT2A
MX1
NCF1
NCF2
NLRP3
NLRX1
NOD2
NUP35
OAS1
OAS2
OAS3
OASL
ORM1
ORM2
PAGE1
PIK3R2
PLA2G7
PLAT
PPBP
PROS1
PSMA3
PSMB8
PSMB9
PSMC6
PTX3
RPL26
RPL38
RSAD2
S100A7
S100A8
S100A9
SERPINE1
SPRR3
STAT1
STAT5A
TAP1
TFPI
TGFB2
THBD
TNFAIP3
TPR
TYROBP
UBA7
UBE2L6
USP18
VCAM1
WAS
XAF1
XPO1
AFAP1L2
AIF1
APOBEC3F
APOBEC3G
BNIP3
CADM1
CALR
CAMK2B
CASP1
CCR7
CD19
CD83
CDK1
CEBPB
CEBPG
CSF2RB
CXCL10
CXCL16
CXCR3
CYSLTR1
DEFB1
EIF4A2
F2
F2R
F5
F7
FAS
FASLG
GAGE1
GBP2
GTF2F2
GZMB
HLA-A
HLA-DMA
HLA-DPB1
HLA-DRB3
HLA-E
HLA-F
HP
HRAS
IFITM2
IFITM3
IL17RB
IL18
IL1R1
IL1RN
IL2RA
IL2RG
IL6R
IL7R
INHBB
IRAK1
ITGB5
KLRC3
KLRC4
KPNA2
KPNA3
LGALS3BP
LY75
LYZ
MAP2K4
MGLL
MIA3
MICB
MX2
NFATC4
NFKB2
NMI
NOS2
NRAS
NUP107
NUP155
NUP205
NUP37
NUP43
NUP85
NUP93
OR2H2
PELI3
PF4
POLR2K
POLR2L
PRF1
PRL
PSG8
PSMA6
PSMB10
PSMB3
PSMB6
PSMD1
PSMD10
PSME2
PTAFR
PTPN1
PYDC1
RBX1
RPL11
RPL12
RPL14
RPL15
RPL37A
RPL4
RPL41
RPLP1
RPS11
RPS14
RPS18
RPS23
RPS27
RPS28
RPS4Y1
RPS6
RPS8
S100A12
SCG2
SEC61B
SEC61G
SEH1L
SH2B1
SHC1
SOD1
TCIRG1
TFF3
TLR3
TPST1
UBE2E1
UBE2N
UMOD
VWF
APOL3
BNIP3L
C2
CCRL1
CD1D
CD36
CD40
CFP
CHST2
COLEC12
CSF2RA
CSH1
CXCL5
CXCR6
DCDC2
DMBT1
ELF3
F12
GBP4
GH1
GPR68
HLA-DPA1
HLA-G
HOXB13
IFI35
IFNG
IL29
IL2RB
KRT1
LYVE1
MAP3K8
MST1R
NOX4
PELI1
PELI2
PROC
PTPN6
RNASEL
RPS12
SP140
STAB1
STAT2
TNFAIP6
TNIP1
VAV1
XCL1
XCL2
