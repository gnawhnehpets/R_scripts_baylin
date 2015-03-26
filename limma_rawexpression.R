#FINAL
# this creates the .cls and .gct files for GSEA
library(limma)
setwd("/home/steve/Desktop/analysis/cohortoneandtwo//")
target <- readTargets("/home/steve/Desktop/analysis/cohortoneandtwo//targets_full.txt")
RG <- read.maimages(target, source="agilent", path="/home/steve/Desktop/analysis/cohortoneandtwo/")
cohorttwo<- c(target$Cy3_Sample, target$Cy5_sample)
cohorttwo
RG$weights <-  matrix(rep(RG$genes$ControlType,ncol(RG$R)),ncol=ncol(RG$R),byrow=F)
RG$weights[RG$genes$ControlType!=0,] <- 0
RG$weights[RG$genes$ControlType==0,] <- 1



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

as.data.frame(cbind(as.matrix(paste(target$Cy3, target$Cy3_Sample, sep="-")),as.matrix(paste(target$Cy5, target$Cy5_sample, sep="-"))))
as.data.frame(sub(".*_1_", "", RG$targets$FileName))
dat <- RG[,-18]
dat <- dat[,-6]
targets <- target[-18,]
targets <- targets[-6,]

as.data.frame(sub(".*_1_", "", dat$targets$FileName))

#1
pos <- c(1,2,7,8,11,14:16,18:22,24,25)
#2
pos <- c(8,14,16,19,20,22,24,25)
#3
pos <- c(8,14,20,24,25)

dat <- dat[,pos]
colnames(dat)
dim(dat)
targets <- targets[pos,]
dim(targets)
targets



match(dat$targets$FileName, targets$FileName)

###################################3
#normalization
#normalization
#normwithin <-normalizeWithinArrays(RG,method='loess',bc.method='none')
normwithin <- normalizeWithinArrays(dat,method='loess',bc.method='normexp', offset=50)
normbetween <- normalizeBetweenArrays(normwithin,method='Aquantile')

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

# Collapse rows where genenames are the same
cy3_agg <- apply(cy3,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
cy5_agg <- apply(cy5,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})})

# remove:
# 1 - pre-1869:   cy5_agg   col 25
# 2 - pre-1544:   cy5_agg   col 20
# 1 - post-1640:  cy5_agg   col 4
# 1 - pre-1876:   cy5_agg   col 1
# 1 - pre-1550:   cy3_agg   col 26
# 1 - pre-1622:   cy3_agg   col 6
# 2 - pre-1811:   cy3_agg   col 3
# 3 - pre-1632:   cy5_agg   col 17
# 3 - pre-1816:   cy5_agg   col 3
# 3 - pre-1868:   cy3_agg   col 23

cbind(as.data.frame(colnames(cy3_agg)),as.data.frame(colnames(cy5_agg)))
#remove2samp


# Double check
#j1 <- match(colnames(cy3_agg), colnames(cy3_agg[,delpos_cy3]))
#j1
#k1 <- which(is.na(j1))
#o1 <- colnames(cy3_agg[,k1])
#o1

#j2 <- match(colnames(cy5_agg), colnames(cy5_agg[,delpos_cy5]))
#j2
#k2 <- which(is.na(j2))
#o2 <- colnames(cy5_agg[,k2])
#o2

########################################################################
########################################################################

dat1 <- cy3_agg
dat2 <- cy5_agg
#dat1 <- dat1[,delpos_cy3]
#dat2 <- dat2[,delpos_cy5]
colnames(dat1)
colnames(dat2)
dim(dat1)
dim(dat2)
# QC: check if unique rownames in pre_agg is the same as unique rownames in pre
length(unique(rownames(dat1)))
length(unique(rownames(cy3)))
length(unique(rownames(dat2)))
length(unique(rownames(cy5)))

#all.dat <- cbind(cy3_agg, cy5_agg)
all.dat <- cbind(dat1,dat2)
colnames(all.dat)
all.cls <-as.matrix(sub("-.*","",colnames(all.dat)))
all.cls
all.col <- colnames(all.dat)
all.col

create_gct <- function(x,y,z){
  genename <- as.data.frame(rownames(x))
  gene_num <- length(rownames(x))
  sample_num <- length(colnames(x))
  description <- matrix(nrow=gene_num, ncol=1)
  new.dat <- cbind(genename,description,x)
  filename <- y
  firstline = "#1.2"
  secondline <- cbind(gene_num, sample_num)
  samplename <- z
  thirdline <- t(as.matrix(c("NAME", "Description", samplename)))
  write.table(firstline, file=y, quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")
  write.table(secondline, file=y, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")
  write.table(thirdline, file=y, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")
  write.table(new.dat, file=y, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")
}

#create_gct(all.dat, "expdat_group1.gct", all.col)
create_gct(all.dat, "expdat_group2.gct", all.col)
#create_gct(all.dat, "bothcohorts_expdat_allsamples_group2.gct", all.col)
#create_gct(all.dat, "bothcohorts_expdat_allsamples_group3.gct", all.col)

create_cls <- function(x,y){
  sample_num <- length(x)
  class_num <- length(unique(x))
  firstline <- paste(sample_num, class_num, "1", sep = " ")
  uniq <- c("#", unique(x))
  secondline <- paste(uniq, collapse = ' ')
  thirdline <- paste(x, collapse = ' ')
  write.table(firstline, file=y, quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")
  write.table(secondline, file=y, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")
  write.table(thirdline, file=y, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")
}

#create_cls(all.cls, "expdat_group1.cls")
create_cls(all.cls, "expdat_group2.cls")
#create_cls(all.cls, "bothcohorts_expdat_allsamples_group2.cls")
#create_cls(all.cls, "bothcohorts_expdat_allsamples_group3.cls")


###################################################################
# PAIRS ONLY ######################################################
###################################################################
# OPTIONAL
# Collapse rows where genenames are the same
#cy3_agg <- apply(cy3,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
#cy5_agg <- apply(cy5,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})})

# QC: check if unique rownames in pre_agg is the same as unique rownames in pre
length(unique(rownames(dat1)))
length(unique(rownames(cy3)))
length(unique(rownames(dat2)))
length(unique(rownames(cy5)))

colnames(dat1)
colnames(dat2)
targets$Cy3_Sample
targets$Cy5_sample
dim(dat1)
dim(dat2)
n1 <- sub(".*-","",colnames(dat1))
n2 <- sub(".*-","",colnames(dat2))
s <- intersect(n1,n2)
s
#m <- match(targets$Cy3_Sample, s)
#m
m1 <- match(n1, s)
m1
m2 <- match(n2, s)
m2


position1 <- which(!is.na(m1))
position1
position2 <- which(!is.na(m2))
position2
p1 <- dat1[,position1]
p2 <- dat2[,position2]
#targets$Cy3_Sample[position1]
colnames(dat1)[position1]
#targets$Cy5_sample[position2]
colnames(dat2)[position2]
as.data.frame(cbind(colnames(dat1)[position1],colnames(dat2)[position2]))
as.data.frame(cbind(colnames(p1),colnames(p2)))
colnames(p1)
colnames(p2)

pair.dat <- cbind(p1,p2)
colnames(pair.dat)
pair.cls <-as.matrix(sub("-.*","",colnames(pair.dat)))
pair.cls
pair.col <- colnames(pair.dat)
pair.col

create_gct <- function(x,y,z){
  genename <- as.data.frame(rownames(x))
  gene_num <- length(rownames(x))
  sample_num <- length(colnames(x))
  description <- matrix(nrow=gene_num, ncol=1)
  new.dat <- cbind(genename,description,x)
  filename <- y
  firstline = "#1.2"
  secondline <- cbind(gene_num, sample_num)
  samplename <- z
  thirdline <- t(as.matrix(c("NAME", "Description", samplename)))
  write.table(firstline, file=y, quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")
  write.table(secondline, file=y, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")
  write.table(thirdline, file=y, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")
  write.table(new.dat, file=y, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")
}

#create_gct(pair.dat, "bothcohorts_expdat_remove2samp_pairsonly.gct", pair.col)
#create_gct(pair.dat, "bothcohorts_expdat_pairsonly_group1.gct", pair.col)
#create_gct(pair.dat, "bothcohorts_expdat_pairsonly_group2.gct", pair.col)
create_gct(pair.dat, "bothcohorts_expdat_pairsonly_group3.gct", pair.col)

create_cls <- function(x,y){
  sample_num <- length(x)
  class_num <- length(unique(x))
  firstline <- paste(sample_num, class_num, "1", sep = " ")
  uniq <- c("#", unique(x))
  secondline <- paste(uniq, collapse = ' ')
  thirdline <- paste(x, collapse = ' ')
  write.table(firstline, file=y, quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")
  write.table(secondline, file=y, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")
  write.table(thirdline, file=y, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")
}

#create_cls(pair.cls, "bothcohorts_expdat_remove2samp_pairsonly.cls")
#create_cls(pair.cls, "bothcohorts_expdat_pairsonly_group1.cls")
#create_cls(pair.cls, "bothcohorts_expdat_pairsonly_group2.cls")
create_cls(pair.cls, "bothcohorts_expdat_pairsonly_group3.cls")
