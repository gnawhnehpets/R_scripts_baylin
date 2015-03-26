#PRIMARY CODE FOR PRELIM2 ANALYSIS 
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
dat <- RG[,-18]
dat <- dat[,-6]
targets <- target[-18,]
targets <- targets[-6,]

as.data.frame(sub(".*_1_", "", dat$targets$FileName))

#1 - no changes
#2
#pos <- c(1,2,7,8,11,14:16,18:22,24,25)
#3
#pos <- c(8,14,16,19,20,22,24,25)
#4
#pos <- c(8,14,20,24,25)
#5
pos <- c(8,14)

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

#1, 2
dat <- normwithin
tar2 <- targets

#Convert MA back to RG
RGb <- RG.MA(normbetween)


# 
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

#>>>> SKIP >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
# filter1 <- c(1,2,3,4,5,6,9,10,12,13,17,20,23,25,26)
# dat <- normwithin[,-filter1]
# tar2 <- targets[-filter1,]
# 
# #4
# filter2 <- c(1,2,3,4,5,6,9,10,12,13,15,16,17,20,21,22,23,25,26)
# dat <- normwithin[,-filter2]
# tar2 <- targets[-filter2,]
# #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 
# 
# normbetween <-normalizeBetweenArrays(dat,method='Aquantile')
# #Remove controls from normwithin/between
# dat <- dat[normbetween$genes$ControlType==0,]
# normbetween <- normbetween[normbetween$genes$ControlType==0,]
# 
# #Convert MA back to RG
# RGb <- RG.MA(normbetween)
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#1,2,3,4

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

#947 probesets have p < .05
length(p.dat[p.dat<.05]) 
#1 - 508
#2 - 142
#3 - 181
#4 - 531
#5 - 832

length(p.dat[p.dat<.01]) #94 probesets have p < .01
#1 - 52
#2 - 12
#3 - 17
#4 - 85
#5 - 158

# 7 genes in group #4filter2
length(p.dat)
b <- .05/length(p.dat)
b

length(p.dat[p.dat<b])
#1 - 0
#2 - 0
#3 - 0 
#4 - 0 
#5 - 0

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
sum(fold.lin>2)
#1, #2
# HBB
# HBD   

#3 - 0

#4
# HBB
# PAIP2 - poly(A) binding protein interacting protein 2

#5 - 0

names(p.dat[p.dat<.05 & abs(fold.lin)>1.5])
#1, #2
# HBB
# HBD

#3 - 0
#4
# A_24_P169843 - NR
# A_24_P67408 - NR
# A_33_P3351615 - NR
# ATG3 - autophagy related 3
# AY927536 - NR
# BTG3 - BTG family, member 3
# C11orf73 - chromosome 11 open reading frame 73
# C16orf80 - chromosome 16 open reading frame 80
# C3orf26 - cms1 ribosomal small subunit homolog (yeast)
# CAPZA2 - capping protein (actin filament) muscle Z-line, alpha 2
# CIRH1A - cirrhosis, autosomal recessive 1A (cirhin)
# CKS1B - CDC28 protein kinase regulatory subunit 1B
# DDX21 - DEAD (Asp-Glu-Ala-Asp) box helicase 21
# DNAJB9 - DnaJ (Hsp40) homolog, subfamily B, member 9
# EI24 - etoposide induced 2.4
# FAM3C - family with sequence similarity 3, member C
# GLRX3 - glutaredoxin 3
# GPN3 - GPN-loop GTPase 3
# HBB - hemoglobin, beta
# KPNA4 - karyopherin alpha 4 (importin alpha 3)
# METAP2 - methionyl aminopeptidase 2
# NAA50 - N(alpha)-acetyltransferase 50, NatE catalytic subunit
# PAIP2 - poly(A) binding protein interacting protein 2
# PCNP - PEST proteolytic signal containing nuclear protein
# PIGH - phosphatidylinositol glycan anchor biosynthesis, class H
# PRDX3 - peroxiredoxin 3
# PTS - 6-pyruvoyltetrahydropterin synthase
# RPF2 - ribosome production factor 2 homolog (S. cerevisiae)
# SDCBP - syndecan binding protein (syntenin)
# TMEM126B - transmembrane protein 126B
# TSPAN13 - tetraspanin 13

#5
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

#1
# ATF1 - activating transcription factor 1
# FAM3C - family with sequence similarity 3, member C
# GPN3 - GPN-loop GTPase 3
# HBB - hemoglobin, beta
# HBD - hemoglobin, delta
# KLF3 - Kruppel-like factor 3 (basic)
# NETO2 - neuropilin (NRP) and tolloid (TLL)-like 2
# PLS1 - plastin 1
# STEAP1 - six transmembrane epithelial antigen of the prostate 1
# TIPIN - TIMELESS interacting protein

names(p.dat[p.dat<.05 & abs(fold.lin)>1.4])
length(names(p.dat[p.dat<.05 & abs(fold.lin)>1.4]))
#1
# A_24_P273245 - NR
# ACN9 - ACN9 homolog (S. cerevisiae)
# ATF1 - activating transcription factor 1
# C14orf142 - chromosome 14 open reading frame 142
# C16orf87 - chromosome 16 open reading frame 87
# CCDC91 - coiled-coil domain containing 91
# FAM3C - family with sequence similarity 3, member C
# GPN3 - GPN-loop GTPase 3
# HBB
# HBD
# HSPA13 - heat shock protein 70kDa family, member 13
# INTS12 - integrator complex subunit 12
# KLF3 - Kruppel-like factor 3 (basic)
# LYRM1 - LYR motif containing 1
# MRPL39 - mitochondrial ribosomal protein L39
# MRPS30 - mitochondrial ribosomal protein S30
# NETO2 - neuropilin (NRP) and tolloid (TLL)-like 2
# PIGH - phosphatidylinositol glycan anchor biosynthesis, class H
# PLS1 - plastin 1
# STEAP1 - six transmembrane epithelial antigen of the prostate 1
# TCTN3 - tectonic family member 3
# TIPIN - TIMELESS interacting protein
# TSGA14 - centrosomal protein 41kDa
# UBE2B - ubiquitin-conjugating enzyme E2B
# XRCC4 - X-ray repair complementing defective repair in Chinese hamster cells 4

#2 - HBB, HBD

#3
# BTG3 - BTG family, member 3
# FAM3C - family with sequence similarity 3, member C
# PIGH - phosphatidylinositol glycan anchor biosynthesis, class H

#4
# A_24_P169843 - NR
# A_24_P67408 - NR
# A_33_P3351615 - NR
# AA627135 - NR
# ABCF2 - ATP-binding cassette, sub-family F (GCN20), member 2
# ATG3 - autophagy related 3
# AY927536 - ribosomal protein L10
# BRD7 - bromodomain containing 7
# BTG3 - BTG family, member 3
# C10orf88 - chromosome 10 open reading frame 88
# C11orf73 - chromosome 11 open reading frame 73
# C11orf74 - chromosome 11 open reading frame 74
# C16orf80 - chromosome 16 open reading frame 80
# C16orf87 - chromosome 16 open reading frame 87
# C2orf76 - chromosome 2 open reading frame 76
# C3orf26 - cms1 ribosomal small subunit homolog (yeast)
# CACYBP - calcyclin binding protein
# CAPZA2 - capping protein (actin filament) muscle Z-line, alpha 2
# CFDP1 - craniofacial development protein 1
# CIRH1A - cirrhosis, autosomal recessive 1A (cirhin)
# CKS1B - CDC28 protein kinase regulatory subunit 1B
# COMMD10 - COMM domain containing 10
# CSTF2T - cleavage stimulation factor, 3' pre-RNA, subunit 2, 64kDa, tau variant
# DCTN5 - dynactin 5 (p25)
# DDX21 - DEAD (Asp-Glu-Ala-Asp) box helicase 21
# DNAJB9 - DnaJ (Hsp40) homolog, subfamily B, member 9
# E2F5 - E2F transcription factor 5, p130-binding
# ECT2 - epithelial cell transforming sequence 2 oncogene
# EED - embryonic ectoderm development
# EI24 - etoposide induced 2.4
# FAM3C - family with sequence similarity 3, member C
# FNBP1L - formin binding protein 1-like
# GLRX3 - glutaredoxin 3
# GPN1 - GPN-loop GTPase 1
# GPN3 - GPN-loop GTPase 3
# HBB
# KPNA4 - karyopherin alpha 4 (importin alpha 3)
# MAGOHB - mago-nashi homolog B (Drosophila)
# MCM2 - minichromosome maintenance complex component 2
# METAP2 - methionyl aminopeptidase 2
# MMGT1 - membrane magnesium transporter 1
# MRPL19 - mitochondrial ribosomal protein L19
# NAA50 - N(alpha)-acetyltransferase 50, NatE catalytic subunit
# NUP93 - nucleoporin 93kDa
# PAIP2 - poly(A) binding protein interacting protein 2
# PCNP - PEST proteolytic signal containing nuclear protein
# PIGH - phosphatidylinositol glycan anchor biosynthesis, class H
# PRDX3 - peroxiredoxin 3
# PSMC3 - proteasome (prosome, macropain) 26S subunit, ATPase, 3
# PTS - 6-pyruvoyltetrahydropterin synthase
# RBM7 - RNA binding motif protein 7
# RPF2 - ibosome production factor 2 homolog (S. cerevisiae)
# SDCBP - syndecan binding protein (syntenin)
# SKP1 - S-phase kinase-associated protein 1
# SNRPE - small nuclear ribonucleoprotein polypeptide E
# TAGLN2 - transgelin 2
# TGS1 - trimethylguanosine synthase 1
# TMEM126B - transmembrane protein 126B
# TSPAN13 - tetraspanin 13
# UAP1 - UDP-N-acteylglucosamine pyrophosphorylase 1
# XM_002342506 - NR
# YWHAQ - tyrosine 3-monooxygenase/tryptophan 5-monooxygenase activation protein, theta polypeptide

#5 
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
#1 
# A_24_P273245

#2 HBB
#2
# CHCHD3 - above
# COPB1 - coatomer protein complex, subunit beta 1
# DEK - DEK oncogene
# DPH2 - DPH2 homolog (S. cerevisiae)
# LPL - above
# MTPN - above
# PPP1CC - above
# PRDX2 - above
# TCEB2 - above

#3 - 0

#4
# ABCF2 - ATP-binding cassette, sub-family F (GCN20), member 2
# C11orf74 - chromosome 11 open reading frame 74
# C2orf76 - chromosome 2 open reading frame 76
# DCTN5 - dynactin 5 (p25)
# EI24 - etoposide induced 2.4
# GPN1 - GPN-loop GTPase 1
# GPN3 - GPN-loop GTPase 3
# MAGOHB - mago-nashi homolog B (Drosophila)
# PIGH - phosphatidylinositol glycan anchor biosynthesis, class H
# PSMC3 - proteasome (prosome, macropain) 26S subunit, ATPase, 3
# PTS - 6-pyruvoyltetrahydropterin synthase
# SNRPE - small nuclear ribonucleoprotein polypeptide E

#5
# CHCHD3 - coiled-coil-helix-coiled-coil-helix domain containing 3
# COPB1 - coatomer protein complex, subunit beta 1
# DEK - DEK oncogene
# DPH2 - DPH2 homolog (S. cerevisiae)
# LPL - lipoprotein lipase
# MTPN - myotrophin
# PPP1CC - protein phosphatase 1, catalytic subunit, gamma isozyme
# PRDX2 - peroxiredoxin 2
# TCEB2 - transcription elongation factor B (SIII), polypeptide 2 (18kDa, elongin B) the protein is not connected to PTS

names(p.dat[p.dat<.01])
#2 
# A_33_P3298830 - NR
# AP4S1 - adaptor-related protein complex 4, sigma 1 subunit
# C11orf46 - ADP-ribosylation factor-like 14 effector protein
# C15orf37 - chromosome 15 open reading frame 37
# C4orf33 - chromosome 4 open reading frame 33
# ENST00000439198 - NR
# ENST00000512519 - NR
# GLP1R - glucagon-like peptide 1 receptor
# HBB
# MMP27 - matrix metallopeptidase 27
# VGLL1 - vestigial like 1 (Drosophila)
# XM_002342506 - NR

#3
# A_33_P3230369 - NR
# C13orf31 - laccase (multicopper oxidoreductase) domain containing 1
# DAXX - death-domain associated protein
# ENST00000340284 - NR
# ENST00000409517 - NR
# FAM170B - family with sequence similarity 170, member B
# FBXL21 - F-box and leucine-rich repeat protein 21 (gene/pseudogene)
# KCNJ5 - potassium inwardly-rectifying channel, subfamily J, member 5
# KRT82 - keratin 82
# LOC100133224 - NR
# MED27 - mediator complex subunit 27
# POFUT1 - mediator complex subunit 27
# PRSS45 - protease, serine, 45
# RGS13 - regulator of G-protein signaling 13
# RN28S1 - RNA, 28S ribosomal 5
# SOHLH1 - spermatogenesis and oogenesis specific basic helix-loop-helix 1
# VGLL1 - vestigial like 1 (Drosophila)

#4
# A_33_P3354574 - NR
# A_33_P3370612 - NR
# A_33_P3377714 - NR
# ABCF2 - ATP-binding cassette, sub-family F (GCN20), member 2
# ANKS3 - ankyrin repeat and sterile alpha motif domain containing 3
# ATP6V0E2 - ATPase, H+ transporting V0 subunit e2
# ATPBD4 - diphthamine biosynthesis 6
# C10orf84 - family with sequence similarity 204, member A
# C11orf74 - chromosome 11 open reading frame 74
# C13orf34 - bora, aurora kinase A activator
# C2orf29 - CCR4-NOT transcription complex, subunit 11
# C2orf60 - tRNA-yW synthesizing protein 5
# C2orf69 - chromosome 2 open reading frame 69
# C2orf76 - chromosome 2 open reading frame 76
# C9orf80 - INTS3 and NABP interacting protein
# CHAC2 - ChaC, cation transport regulator homolog 2 (E. coli)
# DCTN5 - dynactin 5 (p25)
# EI24 - etoposide induced 2.4
# ELAC1 - elaC ribonuclease Z 1
# ENST00000340284 - NR
# ENST00000356822 - NR
# ENST00000391369 - NR
# ENST00000414544 - GSN antisense RNA 1
# ENST00000434635 - NR 
# FAM161A - family with sequence similarity 161, member A
# FAM170B - family with sequence similarity 170, member B
# FAM91A1 - family with sequence similarity 91, member A1
# FCHSD2 - FCH and double SH3 domains 2
# GNAQ - guanine nucleotide binding protein (G protein), q polypeptide
# GPN1 - GPN-loop GTPase 1
# GPN3 - GPN-loop GTPase 3
# GRIPAP1 - GRIP1 associated protein 1
# HACE1 - HECT domain and ankyrin repeat containing E3 ubiquitin protein ligase 1
# HARBI1 - harbinger transposase derived 1
# HIPK3 - homeodomain interacting protein kinase 3
# HPSE - heparanase
# INTS3 - integrator complex subunit 3
# IQSEC3 - IQ motif and Sec7 domain 3
# KLF11 - Kruppel-like factor 11
# L2HGDH - L-2-hydroxyglutarate dehydrogenase
# LARP4 - La ribonucleoprotein domain family, member 4
# LOC100129195 - NR 
# LOC100131101 - NR 
# LOC100288842 - UDP-GlcNAc:betaGal beta-1,3-N-acetylglucosaminyltransferase 5 pseudogene
# LOC390595 - ubiquitin associated protein 1-like
# LOC643802 - u3 small nucleolar ribonucleoprotein protein MPP10-like
# LOC729291 - uncharacterized LOC729291
# LRRC8B - leucine rich repeat containing 8 family, member B
# MAGOHB - mago-nashi homolog B (Drosophila)
# MED27 - mediator complex subunit 27
# MFN1 - mitofusin 1
# NKX3-2 - NK3 homeobox 2
# OGFOD1 - 2-oxoglutarate and iron-dependent oxygenase domain containing 1
# OR10H5 - olfactory receptor, family 10, subfamily H, member 5
# PGM3 - phosphoglucomutase 3
# PIGH - phosphatidylinositol glycan anchor biosynthesis, class H
# PIGM - phosphatidylinositol glycan anchor biosynthesis, class M
# POFUT1 - protein O-fucosyltransferase 1
# POGZ - pogo transposable element with ZNF domain
# POLE3 - polymerase (DNA directed), epsilon 3, accessory subunit
# POLR2A - polymerase (RNA) II (DNA directed) polypeptide A, 220kDa
# PSMC3 - proteasome (prosome, macropain) 26S subunit, ATPase, 3
# PTS - 6-pyruvoyltetrahydropterin synthase
# RAD1 - RAD1 homolog (S. pombe)
# RBM18 - RNA binding motif protein 18
# RCAN1 - regulator of calcineurin 1
# RNF40 - ring finger protein 40, E3 ubiquitin protein ligase
# RPRM - reprimo, TP53 dependent G2 arrest mediator candidate
# RPS4XP21 - ribosomal protein S4X pseudogene 21
# SHKBP1 - SH3KBP1 binding protein 1
# SNRNP27 - small nuclear ribonucleoprotein 27kDa (U4/U6.U5)
# SNRPE - small nuclear ribonucleoprotein polypeptide E
# SOHLH1 - spermatogenesis and oogenesis specific basic helix-loop-helix 1
# SRFBP1 - serum response factor binding protein 1
# SYF2 - SYF2 pre-mRNA-splicing factor
# TAF1A - TATA box binding protein (TBP)-associated factor, RNA polymerase I, A, 48kDa
# VGLL1 - vestigial like 1 (Drosophila)
# VSIG8 - V-set and immunoglobulin domain containing 8
# VSTM1 - V-set and transmembrane domain containing 1
# WDR92 - WD repeat domain 92
# WDSUB1 - WD repeat, sterile alpha motif and U-box domain containing 1
# ZBTB40 - zinc finger and BTB domain containing 40
# ZMAT5 - zinc finger, matrin-type 5
# ZNF230 - zinc finger protein 230
# ZNF828 - chromosome alignment maintaining phosphoprotein 1

names(p.dat[p.dat<.01 & abs(fold.lin)>2])
#1 - 0
#2 - HBB
#3 - 0
#4 - 0
#5 - 0

# Transform the p-value (-1*log10(p-value)) and create a volcano plot with the 
# p-value and fold change vectors (see the lecture notes).  Make sure to use a 
# log10 transformation for the p-value and a log2 (R function log2()) transformation 
# for the fold change.  Draw the horizontal lines at fold values of 2 and -2 (log2 value=1) 
# and the vertical p-value threshold line at p=.05 (remember that it is transformed in the plot).

#template
dev.off()
op <- par(mfrow=c(2,2))
pval_cut = .05
fc_cut = 2
p.trans <- -1 * log10(p.dat)
plot(range(p.trans),range(fold),type='n',xlab='-1*log10(p-value)',ylab='fold change',main='Volcano Plot\npre and post group (p<.05, fc>abs(2))')
points(p.trans,fold,col='black',pch=21,bg=1)
points(p.trans[(p.trans> -log10(pval_cut)&fold>log2(fc_cut))],fold[(p.trans> -log10(pval_cut)&fold>log2(fc_cut))],col=1,bg=2,pch=21)
points(p.trans[(p.trans> -log10(pval_cut)&fold< -log2(fc_cut))],fold[(p.trans> -log10(pval_cut)&fold< -log2(fc_cut))],col=1,bg=3,pch=21)
abline(v= -log10(pval_cut))
abline(h= -log2(fc_cut))
abline(h=log2(fc_cut))

#2
pval_cut = .01
fc_cut = 2
p.trans <- -1 * log10(p.dat)
plot(range(p.trans),range(fold),type='n',xlab='-1*log10(p-value)',ylab='fold change',main='Volcano Plot\npre and post group (p<.01, fc>abs(2))')
points(p.trans,fold,col='black',pch=21,bg=1)
points(p.trans[(p.trans> -log10(pval_cut)&fold>log2(fc_cut))],fold[(p.trans> -log10(pval_cut)&fold>log2(fc_cut))],col=1,bg=2,pch=21)
points(p.trans[(p.trans> -log10(pval_cut)&fold< -log2(fc_cut))],fold[(p.trans> -log10(pval_cut)&fold< -log2(fc_cut))],col=1,bg=3,pch=21)
abline(v= -log10(pval_cut))
abline(h= -log2(fc_cut))
abline(h=log2(fc_cut))

#3
pval_cut = .05
fc_cut = 1.4
p.trans <- -1 * log10(p.dat)
plot(range(p.trans),range(fold),type='n',xlab='-1*log10(p-value)',ylab='fold change',main='Volcano Plot\npre and post group (p<.05, fc>abs(1.4))')
points(p.trans,fold,col='black',pch=21,bg=1)
points(p.trans[(p.trans> -log10(pval_cut)&fold>log2(fc_cut))],fold[(p.trans> -log10(pval_cut)&fold>log2(fc_cut))],col=1,bg=2,pch=21)
points(p.trans[(p.trans> -log10(pval_cut)&fold< -log2(fc_cut))],fold[(p.trans> -log10(pval_cut)&fold< -log2(fc_cut))],col=1,bg=3,pch=21)
abline(v= -log10(pval_cut))
abline(h= -log2(fc_cut))
abline(h=log2(fc_cut))

#4
pval_cut = .01
fc_cut = 1.4
p.trans <- -1 * log10(p.dat)
plot(range(p.trans),range(fold),type='n',xlab='-1*log10(p-value)',ylab='fold change',main='Volcano Plot\npre and post group (p<.01, fc>abs(1.4))')
points(p.trans,fold,col='black',pch=21,bg=1)
points(p.trans[(p.trans> -log10(pval_cut)&fold>log2(fc_cut))],fold[(p.trans> -log10(pval_cut)&fold>log2(fc_cut))],col=1,bg=2,pch=21)
points(p.trans[(p.trans> -log10(pval_cut)&fold< -log2(fc_cut))],fold[(p.trans> -log10(pval_cut)&fold< -log2(fc_cut))],col=1,bg=3,pch=21)
abline(v= -log10(pval_cut))
abline(h= -log2(fc_cut))
abline(h=log2(fc_cut))
par(op)

# #5
# dev.off()
# pval_cut = .1
# fc_cut = 2
# p.trans <- -1 * log10(p.dat)
# plot(range(p.trans),range(fold),type='n',xlab='-1*log10(p-value)',ylab='fold change',main='Volcano Plot\npre and post group differences')
# points(p.trans,fold,col='black',pch=21,bg=1)
# points(p.trans[(p.trans> -log10(pval_cut)&fold>log2(fc_cut))],fold[(p.trans> -log10(pval_cut)&fold>log2(fc_cut))],col=1,bg=2,pch=21)
# points(p.trans[(p.trans> -log10(pval_cut)&fold< -log2(fc_cut))],fold[(p.trans> -log10(pval_cut)&fold< -log2(fc_cut))],col=1,bg=3,pch=21)
# abline(v= -log10(pval_cut))
# abline(h= -log2(fc_cut))
# abline(h=log2(fc_cut))
