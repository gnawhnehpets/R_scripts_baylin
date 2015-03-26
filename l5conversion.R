#FINAL
library(limma)
setwd("/home/steve/Desktop/analysis/cohortoneandtwo//")
target <- readTargets("/home/steve/Desktop/analysis/cohortoneandtwo//targets_full.txt")
RG <- read.maimages(target, source="agilent", path="/home/steve/Desktop/analysis/cohortoneandtwo/")
cohorttwo<- c(target$Cy3_Sample, target$Cy5_sample)
cohorttwo
RG$weights <-  matrix(rep(RG$genes$ControlType,ncol(RG$R)),ncol=ncol(RG$R),byrow=F)
RG$weights[RG$genes$ControlType!=0,] <- 0
RG$weights[RG$genes$ControlType==0,] <- 1

as.data.frame(cbind(as.matrix(paste(target$Cy3, target$Cy3_Sample, sep="-")),as.matrix(paste(target$Cy5, target$Cy5_sample, sep="-"))))
#1
# Remove 1641, 1824 pairs. cols #6,18
dat <- RG[,-18]
dat <- dat[,-6]
targets <- target[-18,]
targets <- targets[-6,]

#2
# for 1910, 1861
dat <- RG[,c(9,15)]
dim(dat)
targets <- target[c(9,15),]
dim(targets)

as.data.frame(cbind(as.matrix(paste(targets$Cy3, targets$Cy3_Sample, sep="-")),as.matrix(paste(targets$Cy5, targets$Cy5_sample, sep="-"))))

#normalization
#normalization
#normwithin <-normalizeWithinArrays(RG,method='loess',bc.method='none')
normwithin <-normalizeWithinArrays(dat,method='loess',bc.method='normexp', offset=50)
normbetween <-normalizeBetweenArrays(normwithin,method='Aquantile')

#Remove controls from normwithin/between
normwithin <- normwithin[normbetween$genes$ControlType==0,]
normbetween <- normbetween[normbetween$genes$ControlType==0,]

#Convert MA back to RG
RGb <- RG.MA(normbetween)

plotDensities(RGb)
names(RGb)
names(dat)
# pre-normalization
boxplot(data.frame(log2(dat$Gb)),main="Green background - pre-normalization", names=paste(targets$Cy3, targets$Cy3_Sample, sep="-"), las=2)
boxplot(data.frame(log2(dat$Rb)),main="Red background - pre-normalization", names=paste(targets$Cy5, targets$Cy5_sample, sep="-"), las=2)

# post-normalization
boxplot(data.frame(log2(RGb$G)),main="Green background - normalized", names=paste(targets$Cy3, targets$Cy3_Sample, sep="-"), las=2)
boxplot(data.frame(log2(RGb$R)),main="Red background - normalized", names=paste(targets$Cy5, targets$Cy5_sample, sep="-"), las=2)

#1 Skip
#2
dat <- normwithin
dim(dat)
tar2 <- targets
tar2
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#3, 4
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
#4# [15] PH1616 pre-Cy3_PH1616 post-Cy5 - degraded
#4# [16] PH1844 post-Cy3_PH1844 pre-Cy5
### [17] PH1623 pre-Cy3_PH1632 pre-Cy5 - non-pair
# [18] PH1827 pre-Cy3_PH1827 post-Cy5
# [19] PH1815 post-Cy3_PH1815 pre-Cy5
### [20] PH1544 post-Cy3_PH1544 pre-Cy5 - filtered
#4# [21] PH1843 post-Cy3_PH1843 pre-Cy5 
#4# [22] PH1871 pre-Cy3_PH1871 post-Cy5
### [23] PH1868 pre-Cy3_PH1887 pre-Cy5 - non-pair
# [24] PH1545 pre-Cy3_PH1545 post-Cy5
### [25] PH1869 post-Cy3_PH1869 pre-Cy5 - pre sample is too degraded
### [26] PH1550 pre-Cy3_PH1600 pre-Cy5 - non-pair

#3
filter1 <- c(1,2,3,4,5,6,9,10,12,13,17,20,23,25,26)
dat <- normwithin[,-filter1]
tar2 <- targets[-filter1,]

#4
filter2 <- c(1,2,3,4,5,6,9,10,12,13,15,16,17,20,21,22,23,25,26)
dat <- normwithin[,-filter2]
tar2 <- targets[-filter2,]


normbetween <-normalizeBetweenArrays(dat,method='Aquantile')
#Remove controls from normwithin/between
dat <- dat[normbetween$genes$ControlType==0,]
normbetween <- normbetween[normbetween$genes$ControlType==0,]

#Convert MA back to RG
RGb <- RG.MA(normbetween)
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# cy3
cy3 = RGb$R
rownames(cy3) <- RGb$genes$GeneName
colnames(cy3) <- paste(tar2$Cy3, tar2$Cy3_Sample, sep="-")
colnames(cy3)
# cy5
cy5 = RGb$G
rownames(cy5) <- RGb$genes$GeneName
colnames(cy5) <- paste(tar2$Cy5, tar2$Cy5_sample, sep="-")
colnames(cy5)


library(genefilter)
#rsd <- rowSds(dat.m)
#rsd <- rowSds(dat)

dat <- cbind(cy3, cy5)
dat <- apply(dat,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
dim(dat)
fullname <- colnames(dat)
fullname
groupname <- sub("-.*", "", colnames(dat))
groupname
colnames(dat) <- fullname
colnames(dat)
pre <- as.matrix(dat[,sub("-.*","",colnames(dat))=="pre"])
dim(pre)
colnames(pre)
post <- as.matrix(dat[,sub("-.*","",colnames(dat))=="post"])
dim(post)
colnames(post)
prepost <- cbind(pre,post)
dim(prepost)
colnames(prepost)
dat.log <- log2(prepost)
dim(dat.log)
prename <- colnames(pre)
prename
postname <- colnames(post)
postname
t.test.all.genes <- function(x, d1, d2){
  x1 <- x[d1]
  x2 <- x[d2]
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  t.out <- t.test(x1, x2, alternative="two.sided", var.equal=T)
  out <- as.numeric(t.out$p.value)
  return(out)
}
prename
p.dat <- apply(dat.log, 1, t.test.all.genes, d1=prename, d2=postname)
length(p.dat)

# A histogram of the p-values and report how many probesets have a p<.05 and p<.01.  
# I divided alpha of 0.05 by the total number of probesets and report how many 
# probesets have a p-value less than this value. This is the Bonferroni correction
# step which is a very conservative p-value thresholding method to account for 
# multiple testing

length(p.dat[p.dat<.05]) #947 probesets have p < .05
length(p.dat[p.dat<.01]) #94 probesets have p < .01
# 7 genes in group #4filter2
length(p.dat)
b <- .05/length(p.dat)
b
length(p.dat[p.dat<b])

par(mfrow=c(1,2))
hist(p.dat,col="lightblue",xlab="p-values",main="P-value dist’n between\npre and post groups",cex.main=0.9)
abline(v=.05,col=2,lwd=2)
hist(-log10(p.dat),col="lightblue",xlab="log10(p-values)", main="-log10(p.dat) dist’n between\npre and groups",cex.main=0.9)
abline(v= -log10(.05),col=2,lwd=2)

# Calculate mean for each gene, fold change between groups
pre.m <- apply(dat.log[,prename], 1, mean, na.rm=T)
post.m <- apply(dat.log[,postname], 1, mean, na.rm=T)
fold <- pre.m-post.m
fold
fold.lin <- 2^fold
names(p.dat[p.dat<.05 & abs(fold.lin)>2])
#1, #3
#[1] "HBB" "HBD"    # both hemoglobin
names(p.dat[p.dat<.05 & abs(fold.lin)>1.5])
#1, #3
#[1] "HBB" "HBD"    # both hemoglobin   

#2
# C1QBP - complement component 1, q subcomponent binding protein
# CHCHD3 - coiled-coil-helix-coiled-coil-helix domain containing 3
# CNIH - cornichon family AMPA receptor auxiliary protein 1
# CYB5B - cytochrome b5 type B (outer mitochondrial membrane)
# EED - embryonic ectoderm development
# EIF1 - eukaryotic translation initiation factor 1
# IER3IP1 - immediate early response 3 interacting protein 1
# LPL - lipoprotein lipase
# MTDH - metadherin
# MTPN - myotrophin
# NAE1 - NEDD8 activating enzyme E1 subunit 1
# PJA1 - praja ring finger 1, E3 ubiquitin protein ligase
# PPP1CC - protein phosphatase 1, catalytic subunit, gamma isozyme
# PRDX3 - peroxiredoxin 3
# STXBP3 - syntaxin binding protein 3
# TCEB2 - transcription elongation factor B (SIII), polypeptide 2 (18kDa, elongin B)
# USP1 - ubiquitin specific peptidase 1

names(p.dat[p.dat<.05 & abs(fold.lin)>1.4])
#[1] "HBB" "HBD"    # both hemoglobin   #No filter, #3
#2
# A_24_P273245 - NR
# C14orf142 - chromosome 14 open reading frame 142
# C1QBP, CHCHD3, CNIH, CYB5B, EED, EIF1, IER3IP1, LPL, MTDH, MTPN, NAE1, PJA1, PPP1CC, PRDX3, USP1
# COPB1 - coatomer protein complex, subunit beta 1
# DEK - DEK oncogene
# DPH2 - DPH2 homolog (S. cerevisiae)
# GRN - granulin
# HBXIP - late endosomal/lysosomal adaptor, MAPK and MTOR activator 5
# HMMR - hyaluronan-mediated motility receptor (RHAMM)
# MCM2 - minichromosome maintenance complex component 2
# MPG - N-methylpurine-DNA glycosylase
# MRPL39 - mitochondrial ribosomal protein L39
# PLEKHG4 - pleckstrin homology domain containing, family G (with RhoGef domain) member 4
# POLR2K - polymerase (RNA) II (DNA directed) polypeptide K, 7.0kDa
# PRDX2 - peroxiredoxin 2
# SNAPIN - SNAP-associated protein
# STXBP3, TCEB2
# TCTN3 - tectonic family member 3
# TSGA14 - centrosomal protein 41kDa
# UBE2Q2 - ubiquitin-conjugating enzyme E2Q family member 2

names(p.dat[p.dat<.01 & abs(fold.lin)>1.4])
#2
# CHCHD3 - above
# COPB1 -
# DEK - 
# DPH2 -
# LPL - above
# MTPN - above
# PPP1CC - above
# PRDX2 - above
# TCEB2 - above

names(p.dat[p.dat<.01])
#4
# [1] "AK098126" "BICD1"    "C13orf26" "C9orf117" "CFP"      "CMTM2"    "GLUD2"
# AK098126 - mitochondrial ribosomal protein L42 pseudogene 1
# BICD1 - bicaudal D homolog 1 (Drosophila)
# C13orf26 - testis expressed 26
# C9orf117 - chromosome 9 open reading frame 117
# CFP - complement factor properdin
# CMTM2 - CKLF-like MARVEL transmembrane domain containing 2
# GLUD2 - glutamate dehydrogenase 2

names(p.dat[p.dat<.01 & abs(fold.lin)>2])


# Transform the p-value (-1*log10(p-value)) and create a volcano plot with the 
# p-value and fold change vectors (see the lecture notes).  Make sure to use a 
# log10 transformation for the p-value and a log2 (R function log2()) transformation 
# for the fold change.  Draw the horizontal lines at fold values of 2 and -2 (log2 value=1) 
# and the vertical p-value threshold line at p=.05 (remember that it is transformed in the plot).

#template
dev.off()
pval_cut = .05
fc_cut = 2
p.trans <- -1 * log10(p.dat)
plot(range(p.trans),range(fold),type='n',xlab='-1*log10(p-value)',ylab='fold change',main='Volcano Plot\npre and post group differences')
points(p.trans,fold,col='black',pch=21,bg=1)
points(p.trans[(p.trans> -log10(pval_cut)&fold>log2(fc_cut))],fold[(p.trans> -log10(pval_cut)&fold>log2(fc_cut))],col=1,bg=2,pch=21)
points(p.trans[(p.trans> -log10(pval_cut)&fold< -log2(fc_cut))],fold[(p.trans> -log10(pval_cut)&fold< -log2(fc_cut))],col=1,bg=3,pch=21)
abline(v= -log10(pval_cut))
abline(h= -log2(fc_cut))
abline(h=log2(fc_cut))

#2
dev.off()
pval_cut = .01
fc_cut = 2
p.trans <- -1 * log10(p.dat)
plot(range(p.trans),range(fold),type='n',xlab='-1*log10(p-value)',ylab='fold change',main='Volcano Plot\npre and post group differences')
points(p.trans,fold,col='black',pch=21,bg=1)
points(p.trans[(p.trans> -log10(pval_cut)&fold>log2(fc_cut))],fold[(p.trans> -log10(pval_cut)&fold>log2(fc_cut))],col=1,bg=2,pch=21)
points(p.trans[(p.trans> -log10(pval_cut)&fold< -log2(fc_cut))],fold[(p.trans> -log10(pval_cut)&fold< -log2(fc_cut))],col=1,bg=3,pch=21)
abline(v= -log10(pval_cut))
abline(h= -log2(fc_cut))
abline(h=log2(fc_cut))

#3
dev.off()
pval_cut = .05
fc_cut = 1.4
p.trans <- -1 * log10(p.dat)
plot(range(p.trans),range(fold),type='n',xlab='-1*log10(p-value)',ylab='fold change',main='Volcano Plot\npre and post group differences')
points(p.trans,fold,col='black',pch=21,bg=1)
points(p.trans[(p.trans> -log10(pval_cut)&fold>log2(fc_cut))],fold[(p.trans> -log10(pval_cut)&fold>log2(fc_cut))],col=1,bg=2,pch=21)
points(p.trans[(p.trans> -log10(pval_cut)&fold< -log2(fc_cut))],fold[(p.trans> -log10(pval_cut)&fold< -log2(fc_cut))],col=1,bg=3,pch=21)
abline(v= -log10(pval_cut))
abline(h= -log2(fc_cut))
abline(h=log2(fc_cut))


#4
dev.off()
p.trans <- -1 * log10(p.dat)
plot(range(p.trans),range(fold),type='n',xlab='-1*log10(p-value)',ylab='fold change',main='Volcano Plot\npre and postgroup differences')
points(p.trans,fold,col='black',pch=21,bg=1)
points(p.trans[(p.trans> -log10(.05))],fold[(p.trans> -log10(.05))],col=1,bg=2,pch=21)
points(p.trans[(p.trans> -log10(.05))],fold[(p.trans> -log10(.05))],col=1,bg=3,pch=21)
abline(v= -log10(.05))
abline(h= -log2(2))
abline(h=log2(2))
