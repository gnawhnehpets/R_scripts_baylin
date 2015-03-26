# PRIMARY CODE FOR PRELIM2 ANALYSIS 
# FINAL
# library(limma)
# setwd("/home/steve/Desktop/analysis/cohortoneandtwo//")
# target <- readTargets("/home/steve/Desktop/analysis/cohortoneandtwo//targets_full.txt")
# RG <- read.maimages(target, source="agilent", path="/home/steve/Desktop/analysis/cohortoneandtwo/")
# cohorttwo<- c(target$Cy3_Sample, target$Cy5_sample)
# cohorttwo
# RG$weights <-  matrix(rep(RG$genes$ControlType,ncol(RG$R)),ncol=ncol(RG$R),byrow=F)
# RG$weights[RG$genes$ControlType!=0,] <- 0
# RG$weights[RG$genes$ControlType==0,] <- 1
# 
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
# as.data.frame(cbind(as.matrix(paste(target$Cy3, target$Cy3_Sample, sep="-")),as.matrix(paste(target$Cy5, target$Cy5_sample, sep="-"))))
# as.data.frame(sub(".*_1_", "", RG$targets$FileName))
# dat <- RG[,-18]
# dat <- dat[,-6]
# targets <- target[-18,]
# targets <- targets[-6,]
# 
# as.data.frame(sub(".*_1_", "", dat$targets$FileName))
# 
# #1
# pos <- c(1,2,7,8,11,14:16,18:22,24,25)
# #2
# pos <- c(8,14,16,19,20,22,24,25)
# #3
# pos <- c(8,14,20,24,25)
# 
# dat <- dat[,pos]
# colnames(dat)
# dim(dat)
# targets <- targets[pos,]
# dim(targets)
# targets
# 
# as.data.frame(cbind(as.matrix(paste(targets$Cy3, targets$Cy3_Sample, sep="-")),as.matrix(paste(targets$Cy5, targets$Cy5_sample, sep="-"))))
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

#FINAL
rm(list=ls())
library(limma)
setwd("/home/steve/Desktop/analysis/cohortoneandtwo//")
target <- readTargets("/home/steve/Desktop/analysis/cohortoneandtwo//targets_full.txt")
target
RG <- read.maimages(target, source="agilent", path="/home/steve/Desktop/analysis/cohortoneandtwo/")
cohorttwo<- c(target$Cy3_Sample, target$Cy5_sample)
cohorttwo
RG$weights <-  matrix(rep(RG$genes$ControlType,ncol(RG$R)),ncol=ncol(RG$R),byrow=F)
RG$weights[RG$genes$ControlType!=0,] <- 0
RG$weights[RG$genes$ControlType==0,] <- 1

# x2 1        post-PH1876        pre-PH1876
# x2 2         pre-PH1900       post-PH1900
# x1 3         pre-PH1811        pre-PH1816
# x1 4         pre-PH1636       post-PH1640
# x1 5         pre-PH1892        pre-PH1902
# x1 6         pre-PH1622        pre-PH1631
# x2 7         pre-JHH005       post-JHH005
# 8        post-PH1910        pre-PH1910
# x9         pre-PH1604        pre-PH1606
# x1 10       post-JHH004        pre-JHH004
# x2 11       post-PH1612        pre-PH1612
# x1 12        pre-PH1913        pre-PH1886
# x1 13        pre-PH1635        pre-PH1644
# 14        pre-PH1861       post-PH1861
# x2 15        pre-PH1616       post-PH1616
# x3 16       post-PH1844        pre-PH1844    #NO DRY ICE
# x1 17        pre-PH1623        pre-PH1632
# x2 18        pre-PH1827       post-PH1827
# x3 19       post-PH1815        pre-PH1815
# 20       post-PH1544        pre-PH1544
# x2 21       post-PH1843        pre-PH1843
# x3 22        pre-PH1871       post-PH1871
# x1 23        pre-PH1868        pre-PH1887
# 24        pre-PH1545       post-PH1545
# 25       post-PH1869        pre-PH1869
# x1 26        pre-PH1550        pre-PH1600

as.data.frame(cbind(as.matrix(paste(target$Cy3, target$Cy3_Sample, sep="-")),as.matrix(paste(target$Cy5, target$Cy5_sample, sep="-"))))
as.data.frame(sub(".*_1_", "", RG$targets$FileName))
#remove outlier samples PH1824/1641
dat <- RG[,-18]
dat <- dat[,-6]
targets <- target[-18,]
targets <- targets[-6,]

as.data.frame(sub(".*_1_", "", dat$targets$FileName))

#1 - no changes
#2
pos <- c(1,2,7,8,11,14:16,18:22,24,25)
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

##################################################################################
##################################################################################
##################################################################################

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

####################################################
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
as.data.frame(colnames(pppre),colnames(pppost))
####################################################
#r <- cbind(cy3, cy5)
r <- cbind(pppre,pppost)
colnames(r)
dim(r)

new.dat <- apply(r,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
dim(new.dat)
colnames(new.dat)

##################################################################################
##################################################################################
##################################################################################

#r <- combat
# Cluster dendogram
rt <- t(r)      #transpose dat
rt.dist <- dist(rt,method="euclidean")  # calculate distance
rt.clust <- hclust(rt.dist,method="single")  # calculate clusters
plot(rt.clust,labels=names(rt),cex=0.75)	# plot cluster tree

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

##################################################################################
##################################################################################
##################################################################################

#1,2,3,4
#new.datfil <- cbind(new.dat[,sub("-.*","",colnames(new.dat))=="pre"],new.dat[,sub("-.*","",colnames(new.dat))=="post"])
# all_pres <- new.dat[,sub("-.*","",colnames(new.dat))=="pre"]
# colnames(new.dat)
# all_posts <- new.dat[,sub("-.*","",colnames(new.dat))=="post"]
# sub("-.*","",colnames(new.dat))=="post"
# as.data.frame(colnames(all_pres),colnames(all_posts))
new.datfil <- new.dat




as.matrix(colnames(new.dat))
colnames(new.datfil)
dim(new.datfil)
dim(new.dat)
#new.datfil <- cbind(new.dat[,pre_pos], new.dat[,-pre_pos])

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
<- new.datfil["",]

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
# HC #####################################
##########################################
#e-cad
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
dim(new.datfil)
colnames(new.datfil)
dim(new.datfil)


genes <- rbind(ecad, p16, tms1, timp1, timp2, timp3, timp4, mlh1, mgmt, pd1, pdl1, pdl2)
allgenes <- genes
genename <- "all genes"
colors <- colors()[c(130, 8, 461, 47, 84, 53, 554, 418, 640, 275, 143, 258)]
rownames(allgenes) <- as.matrix(c("ECAD", "P16", "TMS1", "TIMP1", "TIMP2", "TIMP3", "TIMP4", "MLH1", "MGMT", "PD1", "PDL1", "PDL2"))
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




dev.off()
plot(c(1,ncol(r)),range(genes),type='n',main=paste("Profile plot of",genename, sep=" "),xlab="Samples",ylab="Expression intensity",axes=F)
#1
axis(side=1,at=c(1:30),labels=colnames(genes),cex.axis=0.4,las=2)
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