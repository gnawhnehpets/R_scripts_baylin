# FOR REFERENCE
library(limma)
setwd("/home/steve/Desktop/analysis/cohortoneandtwo//")
target <- readTargets("/home/steve/Desktop/analysis/cohortoneandtwo//targets_full.txt")
RG <- read.maimages(target, source="agilent", path="/home/steve/Desktop/analysis/cohortoneandtwo/")
cohorttwo<- c(target$Cy3_Sample, target$Cy5_sample)
RG$weights <-  matrix(rep(RG$genes$ControlType,ncol(RG$R)),ncol=ncol(RG$R),byrow=F)
RG$weights[RG$genes$ControlType!=0,] <- 0
RG$weights[RG$genes$ControlType==0,] <- 1
# Remove 1641, 1824 pairs. cols #6,18
dat <- RG[,-18]
dat <- dat[,-6]
dim(dat)
dim(targets)
targets <- target[-18,]
targets <- targets[-6,]
dim(targets)

#normalization
#normwithin <-normalizeWithinArrays(RG,method='loess',bc.method='none')

#Includes 4 bad samples (2 bad pairs)
#normwithin <-normalizeWithinArrays(RG,method='loess',bc.method='normexp', offset=50)

normwithin <-normalizeWithinArrays(dat,method='loess',bc.method='normexp', offset=50)
normbetween <-normalizeBetweenArrays(normwithin,method='Aquantile')

#Remove controls from normwithin/between
normwithin <- normwithin[normbetween$genes$ControlType==0,]

plotDensities(normwithin)

#normbetween <- normbetween[normbetween$genes$ControlType==0,]

#Convert MA back to RG
#ma2rg <- RG.MA(nc_normbetween)
RGb <- RG.MA(normbetween)

plotDensities(RGb)

names(RGb)
names(RGb)

# pre-normalization
boxplot(data.frame(log2(dat$Gb)),main="Green background - pre-normalization", names=paste(targets$Cy3, targets$Cy3_Sample, sep="-"), las=2)
boxplot(data.frame(log2(dat$Rb)),main="Red background - pre-normalization", names=paste(targets$Cy5, targets$Cy5_sample, sep="-"), las=2)

# post-normalization
boxplot(data.frame(log2(RGb$G)),main="Green background - normalized", names=paste(targets$Cy3, targets$Cy3_Sample, sep="-"), las=2)
boxplot(data.frame(log2(RGb$R)),main="Red background - normalized", names=paste(targets$Cy5, targets$Cy5_sample, sep="-"), las=2)


# For all samples
cy3_num <- as.character(targets$numpre)
cy5_num <- as.character(targets$numpost)
length(cy3_num)
cy5_num

cy3temp <- strsplit(cy3_num, split="")
cy5temp <- strsplit(cy5_num, split="")
cy3_num <- matrix(unlist(cy3temp), ncol=3, byrow=TRUE)
cy5_num <- matrix(unlist(cy5temp), ncol=3, byrow=TRUE)
dim(cy3_num)
cy5_num

# get date/location/rin data for batch effects
cy3date <- as.matrix(cy3_num[,1])
cy3location <- as.matrix(cy3_num[,2])
cy3rin <- as.matrix(as.numeric(cy3_num[,3]))
head(cy3date)
dim(cy3location)
dim(cy3rin)
cy5date <- as.matrix(cy5_num[,1])
cy5location <- as.matrix(cy5_num[,2])
cy5rin <- as.matrix(as.numeric(cy5_num[,3]))
dim(cy5date)
dim(cy5location)
dim(cy5rin)

cy3names = as.matrix(paste(targets$Cy3, targets$Cy3_Sample, sep="-"))
cy5names = as.matrix(paste(targets$Cy5, targets$Cy5_sample, sep="-"))
cy3names
cy5names
# cy3
sva_cy3 = RGb$R
rownames(sva_cy3) <- RGb$genes$GeneName
colnames(sva_cy3) <- cy3names
colnames(sva_cy3)
# cy5
sva_cy5 = RGb$G
rownames(sva_cy5) <- RGb$genes$GeneName
colnames(sva_cy5) <- cy5names
colnames(sva_cy5)

# one matrix (edata)
# with "group-sample" naming format
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sva <- cbind(sva_cy3, sva_cy5)

# with "group" naming format
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
colnames(sva_cy3) <- targets$Cy3
colnames(sva_cy5) <- targets$Cy5
sva2 <- cbind(sva_cy3, sva_cy5)

sedata <- apply(sva,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
sedata2 <- apply(sva2,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 


# Exploratory
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
library(ggplot2)
r <- sedata
#r <- combat
# Cluster dendogram
rt <- t(r)      #transpose dat
rt.dist <- dist(rt,method="euclidean")  # calculate distance
rt.clust <- hclust(rt.dist,method="single")  # calculate clusters
plot(rt.clust,labels=names(rt),cex=0.75)	# plot cluster tree

# CV v mean plot
r.mean <- apply(log2(r),2,mean)  	# calculate mean for each sample
r.sd <- sqrt(apply(log2(r),2,var))		# calculate st.deviation for each sample
r.cv <- r.sd/r.mean		

plot(r.mean,r.cv,main="colon trial dataset\nSample CV vs. Mean",xlab="Mean",ylab="CV",col='blue',cex=1.5,type="n")
points(r.mean,r.cv,bg="lightblue",col=1,pch=21)
text(r.mean,r.cv,label=dimnames(r)[[2]],pos=1,cex=0.6)
# with sedata, we would eliminate pre-PH1550, pre-PH1811(?), pre-PH1869 as outliers

# Avg correlation plot
r.cor <- cor(r)
r.avg <- apply(r.cor,1,mean)
par(oma=c(3,0.1,0.1,0.1))
plot(c(1,length(r.avg)),range(r.avg),type="n",xlab="",ylab="Avg r",main="Avg correlation of colon trial samples",axes=F)
points(r.avg,bg="red",col=1,pch=21,cex=1.25)
axis(1,at=c(1:length(r.avg)),labels=dimnames(r)[[2]],las=2,cex.lab=0.4,cex.axis=0.6)
axis(2)
abline(v=seq(0.5,62.5,1),col="grey")

# Pearson's correlation plot
library(gplots)
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

datpca <- prcomp(t(sedata), cor=F)

plot(
  range(datpca$x[, 1]),range(datpca$x[, 2]), 
  type = "n", xlab = "pre", ylab = "post", 
  main = "PCA of colon trial data\npre vs post"
)
points(datpca$x[, 1][targets$Cy3 == "pre"], datpca$x[, 2][targets$Cy3 == "pre"], bg = "Red", pch = 21)
text(datpca$x[, 1][targets$Cy3 == "pre"], datpca$x[, 2][targets$Cy3 == "pre"], label=names(datpca$x[,1][targets$Cy3 == "pre"]),pos=1,cex=0.5)
points(datpca$x[, 1][targets$Cy5 == "post"], datpca$x[, 2][targets$Cy5 == "post"], bg = "Blue", pch=21)
text(datpca$x[, 1][targets$Cy5 == "post"], datpca$x[, 2][targets$Cy5 == "post"], label=names(datpca$x[,1][targets$Cy5 == "post"]),pos=1,cex=0.5)
leg.names <- c("pre", "post")       
legend(-30,12,leg.names,leg.col,pch=15,cex=.7,horiz=F)


# regex for gene name
dimnames(sedata)[[2]]
#TMS1
length(grep("tms", dimnames(sedata)[[1]], ignore.case = TRUE))
rownames(sedata)[grep("tms", dimnames(sedata)[[1]], ignore.case = TRUE)]
#E-cad
length(grep("cad", dimnames(sedata)[[1]], ignore.case = TRUE))
rownames(sedata)[grep("cad", dimnames(sedata)[[1]], ignore.case = TRUE)]
#TIMP
length(grep("timp", dimnames(sedata)[[1]], ignore.case = TRUE))
rownames(sedata)[grep("timp", dimnames(sedata)[[1]], ignore.case = TRUE)]
#MLH1
length(grep("mlh", dimnames(sedata)[[1]], ignore.case = TRUE))
rownames(sedata)[grep("mlh", dimnames(sedata)[[1]], ignore.case = TRUE)]
#MGMT
length(grep("mgmt", dimnames(sedata)[[1]], ignore.case = TRUE))
rownames(sedata)[grep("mgmt", dimnames(sedata)[[1]], ignore.case = TRUE)]
#p16
length(grep("p16", dimnames(sedata)[[1]], ignore.case = TRUE))
rownames(sedata)[grep("p16", dimnames(sedata)[[1]], ignore.case = TRUE)]

timp1 <- sedata["TIMP1",]
timp2 <- sedata["TIMP2",]
timp3 <- sedata["TIMP3",]
timp4 <- sedata["TIMP4",]
mlh1 <- sedata["MLH1",]
mgmt <- sedata["MGMT",]

genes <- rbind(timp1, timp2, timp3, timp4, mlh1, mgmt)
genes
colnames(genes)
dim(genes)

plot(c(1,ncol(r)),range(genes),type='n',main="Profile plot of genes",xlab="Samples",ylab="Expression intensity",axes=F)
axis(side=1,at=c(1:52),labels=colnames(genes),cex.axis=0.4,las=2)
axis(side=2)
for(i in 1:length(genes)) {
  r.y <- as.numeric(genes[i,])
  lines(c(1:ncol(genes)),r.y,col=i,lwd=2)
}

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

######################################
# Remove outliers ####################
######################################


# Remove outliers
as.data.frame(colnames(sedata))
dim(sedata)

# remove:
# 1 - pre-1869: col 51
# 2 - pre-1544: col 46
# 1 - pre-1640: col 30
# 1 - pre-1876: col 27
# 1 - pre-1550: col 26
# 1 - pre-1622: col 6
# 2 - pre-1811: col 3
# 3 - pre-1632: col 43
# 3 - pre-1816: col 29
# 3 - pre-1868: col 23

#1
#delpos <- c(51,30,27,26,6)
#2
#delpos <- c(51,46,30,27,26,6,3)
#3
delpos <- c(51,46,43,30,29,27,26,23,6,3)

new.dat <- sedata[,-delpos]
dimnames(new.dat)
length(grep("pre", dimnames(new.dat)[[2]], ignore.case = TRUE))
colnames(new.dat)[grep("pre", dimnames(new.dat)[[2]], ignore.case = TRUE)]
length(grep("post", dimnames(new.dat)[[2]], ignore.case = TRUE))
colnames(new.dat)[grep("post", dimnames(new.dat)[[2]], ignore.case = TRUE)]

######################################
# sva ################################
######################################

# BATCHES
# create single matrix for date/location/rin
date <- as.data.frame(rbind(cy3date, cy5date))
location <- as.data.frame(rbind(cy3location, cy5location))
rin <- as.data.frame(rbind(cy3rin, cy5rin))
dim(date)
dim(location)
dim(rin)
delpos
date <- as.data.frame(date[-delpos,])
location <- as.data.frame(location[-delpos,])
rin <- as.data.frame(rin[-delpos,])
colnames(date) <- "date"
colnames(location) <- "location"
colnames(rin) <- "rin"
dim(date)
dim(location)
dim(rin)

# For all samples
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
group_cy3 = as.data.frame(targets$Cy3)
names(group_cy3) <- "group"
group_cy5 = as.data.frame(targets$Cy5)
names(group_cy5) <- "group"
group <- rbind(group_cy3, group_cy5)
group <- as.data.frame(group[-delpos,])
group
colnames(group) <- "group"
#sample <- as.data.frame(rbind(cy3names, cy5names))
#sample <- as.data.frame(c(1:(length(targets$Cy5)*2)))
sample <- as.data.frame(c(1:length(rownames(group))))
sample <- as.data.frame(colnames(new.dat))
colnames(sample) <- "sample"
sample
# returns false/true statement based on RIN cutoff aka QC of samples

#rincutoff <- as.data.frame(factor(rin>4))
rincutoff <- as.data.frame(rin>4)
rincutoff[rincutoff=="FALSE"] <- "below"
rincutoff[rincutoff=="TRUE"] <- "above"
colnames(rincutoff) <- "rincutoff"
dim(sample)
dim(group)
dim(rincutoff)
#phen <- cbind(sample, group, date, location, rincutoff)
colnames(date) <- "batch"
phen <- cbind(sample, group, date, location)
phen
#rownames(phen) <- colnames(sedata)
rownames(phen) <- colnames(new.dat)
phen
library(sva)

# Create a model matrix for adjustment variables and variable of interest (rincutoff)
# http://www.bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf page 4
# mod = model.matrix(~as.factor(cancer), data=pheno)
# mod0 = model.matrix(~1,data=pheno)
#smod = model.matrix(~as.factor(rincutoff), data=phen)
mod = model.matrix(~as.factor(group), data=phen)
mod0 = model.matrix(~1,data=phen)
# This returns an expression matrix, with the same dimensions as original dataset. This new expression
# matrix is adjusted for batch.
combat = ComBat(dat=new.dat, batch=phen$batch, mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE)
pValuesComBat = f.pvalue(combat,smod,smod0)
qValuesComBat = p.adjust(pValuesComBat,method="BH")

dim(combat)

######################################
# combat = adjusted(new.dat) #########
######################################

#cy3_sva <- combat[,1:16]
num_of_samples <- length(colnames(combat))
num_of_samples

# POSITION OF LAST cy3 sample THAT WASN'T FILTERED
last_cy3 <- "post-PH1869"
colnames(combat)==last_cy3
which(colnames(combat)==last_cy3)
dim(combat)

cy3_sva <- combat[,1:which(colnames(combat)==last_cy3)]
colnames(cy3_sva) <- colnames(combat)[1:which(colnames(combat)==last_cy3)]
cy5_sva <- combat[,(which(colnames(combat)==last_cy3)+1):num_of_samples]
colnames(cy5_sva) <- colnames(combat)[(which(colnames(combat)==last_cy3)+1):num_of_samples]

########################################################################
########################################################################
#1
#is.na(match(colnames(cy3_sva),colnames(combat[,1:24])))
#is.na(match(colnames(cy5_sva),colnames(combat[,25:47])))
#2
#is.na(match(colnames(cy3_sva),colnames(combat[,1:23])))
#is.na(match(colnames(cy5_sva),colnames(combat[,24:45])))
#3
is.na(match(colnames(cy3_sva),colnames(combat[,1:22])))
is.na(match(colnames(cy5_sva),colnames(combat[,23:42])))
########################################################################
########################################################################

cbind(tail(colnames(cy3_sva)),head(colnames(cy5_sva)))
dim(cy3_sva)
dim(cy5_sva)
colnames(cy3_sva)
colnames(cy5_sva)
colnames(combat)
cy3names <- colnames(combat)[1:which(colnames(combat)==last_cy3)]
cy5names <- colnames(combat)[(which(colnames(combat)==last_cy3)+1):num_of_samples]
pre_names <- c(colnames(cy3_sva)[grep("pre", colnames(cy3_sva), ignore.case = TRUE)], colnames(cy5_sva)[grep("pre", colnames(cy5_sva), ignore.case = TRUE)])
post_names <- c(colnames(cy5_sva)[grep("post", colnames(cy5_sva), ignore.case = TRUE)], colnames(cy3_sva)[grep("post", colnames(cy3_sva), ignore.case = TRUE)])
pre_names
post_names
colnames(cy3_sva)<- sub("-.*","",colnames(cy3_sva))
colnames(cy5_sva)<- sub("-.*","",colnames(cy5_sva))
colnames(cy3_sva)
colnames(cy5_sva)
# NOTE: IMPORTANT! Cy3==PRE @position1 will be Cy5==POST
pre_sva <- cbind(cy3_sva[,colnames(cy3_sva)=="pre"], cy5_sva[,colnames(cy5_sva)=="pre"])
dim(pre_sva)
post_sva <- cbind(cy5_sva[,colnames(cy5_sva)=="post"], cy3_sva[,colnames(cy3_sva)=="post"])
dim(post_sva)
#pre_names <- c(targets$Cy3_Sample[colnames(cy3_sva)=="pre"], targets$Cy5_sample[colnames(cy5_sva)=="pre"])
#pre_names <- c(colnames(cy3_sva)[grep("pre", colnames(cy3_sva), ignore.case = TRUE)], colnames(cy5_sva)[grep("pre", colnames(cy5_sva), ignore.case = TRUE)])
pre_names
#post_names <- c(targets$Cy3_Sample[colnames(cy5_sva)=="post"], targets$Cy5_sample[colnames(cy3_sva)=="post"])
#post_names <- c(colnames(cy3_sva)[grep("post", colnames(cy3_sva), ignore.case = TRUE)], colnames(cy5_sva)[grep("post", colnames(cy5_sva), ignore.case = TRUE)])
post_names

#Make colnames "group-sample"
#colnames(pre_sva) <- paste(colnames(pre_sva),pre_names,sep="-")
colnames(pre_sva) <- pre_names
colnames(pre_sva)

#colnames(post_sva) <- paste(colnames(post_sva),post_names,sep="-")
colnames(post_sva) <- post_names
colnames(post_sva)

presva_mean <- as.data.frame(apply(pre_sva,1,mean))
rownames(presva_mean) <- rownames(pre_sva)
postsva_mean <- as.data.frame(apply(post_sva,1,mean))
rownames(postsva_mean) <- rownames(post_sva)
Msva <- log2(presva_mean/postsva_mean)

presva_median <- as.data.frame(apply(pre_sva,1,median))
rownames(presva_median) <- rownames(pre_sva)
postsva_median <- as.data.frame(apply(post_sva,1,median))
rownames(postsva_median) <- rownames(post_sva)
Medsva <- log2(presva_median/postsva_median)

#Create ranked list files, using M-values of comparison groups
createrankedlist <- function(x,y){
  genename <- as.data.frame(rownames(x))
  names(genename) <- "genename"
  new.dat <- cbind(genename,x)
  filename = y
  write.table(new.dat, file=y, quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")
}
#createrankedlist(Msva, "bothcohorts_preVpost-MEAN_SVA_allsamples.rnk")
#createrankedlist(Msva, "bothcohorts_preVpost-MEDIAN_SVA_allsamples.rnk")
#createrankedlist(Msva, "bothcohorts_preVpost-MEAN_SVA_allsamples2.rnk")
#createrankedlist(Msva, "bothcohorts_preVpost-MEDIAN_SVA_allsamples2.rnk")
createrankedlist(Msva, "bothcohorts_preVpost-MEAN_SVA_allsamples3.rnk")
createrankedlist(Msva, "bothcohorts_preVpost-MEDIAN_SVA_allsamples3.rnk")
>>>>>>>>>>>>>>>
  
  ###################################################################
# PAIRS ONLY ######################################################
###################################################################

#pre_sva <- cbind(cy3_sva[,colnames(cy3_sva)=="pre"], cy5_sva[,colnames(cy5_sva)=="pre"])
#post_sva <- cbind(cy5_sva[,colnames(cy5_sva)=="post"], cy3_sva[,colnames(cy3_sva)=="post"])
#colnames(pre_sva)
#colnames(post_sva)
#group.sample.pre <- colnames(pre_sva)
group.sample.pre <- pre_names
group.sample.pre
#group.sample.post <- colnames(post_sva)
group.sample.post <- post_names
group.sample.post
#colnames(pre_sva) <- sub(".*-","",colnames(pre_sva))
#colnames(post_sva) <- sub(".*-","",colnames(post_sva))
colnames(pre_sva) <- sub(".*-","", group.sample.pre)
colnames(post_sva) <- sub(".*-","", group.sample.post)
colnames(pre_sva)
colnames(post_sva)

#Only keep pairs
s <- intersect(colnames(pre_sva), colnames(post_sva))
#s <- intersect()
s
m_pre <- match(colnames(pre_sva), s)
m_pre
m_post <- match(colnames(post_sva), s)
m_post

position_pre <- which(!is.na(m_pre))
position_pre
position_post <- which(!is.na(m_post))
position_post

pre_pairs <- pre_sva[,position_pre]
dim(pre_pairs)
post_pairs <- post_sva[,position_post]
dim(post_pairs)

length(group.sample.pre)
length(group.sample.pre[position_pre])
group.sample.pre[position_pre]
group.sample.post[position_post]
colnames(pre_pairs) <- group.sample.pre[position_pre]
colnames(post_pairs) <- group.sample.post[position_post]
colnames(pre_pairs)
colnames(post_pairs)

# QC: check if unique rownames in pre_agg is the same as unique rownames in pre
length(unique(rownames(pre_sva)))
length(unique(rownames(pre_pairs)))
length(unique(rownames(post_sva)))
length(unique(rownames(post_pairs)))

# TAKING AVERAGE OF ALL PREs vs AVERAGE OF ALL POSTs
# Create single expression values from multiple samples (paired)
pre_mean_pair <- as.data.frame(apply(pre_pairs,1,mean))
rownames(pre_mean_pair) <- rownames(pre_pairs)
post_mean_pair <- as.data.frame(apply(post_pairs,1,mean))
rownames(post_mean_pair) <- rownames(post_pairs)
Mmeanpair <- log2(pre_mean_pair/post_mean_pair)

pre_median_pair <- as.data.frame(apply(pre_pairs,1,median))
rownames(pre_median_pair) <- rownames(pre_pairs)
post_median_pair <- as.data.frame(apply(post_pairs,1,median))
rownames(post_median_pair) <- rownames(post_pairs)
Mmedpair <- log2(pre_median_pair/post_median_pair)

#Create ranked list files, using M-values of comparison groups
createrankedlist <- function(x,y){
  genename <- as.data.frame(rownames(x))
  names(genename) <- "genename"
  new.dat <- cbind(genename,x)
  filename = y
  write.table(new.dat, file=y, quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")
}

#createrankedlist(Msva, "bothcohorts_preVpost-MEAN_SVA_pairs.rnk")
#createrankedlist(Msva, "bothcohorts_preVpost-MEDIAN_SVA_pairs.rnk")
#createrankedlist(Msva, "bothcohorts_preVpost-MEAN_SVA_pairs2.rnk")
#createrankedlist(Msva, "bothcohorts_preVpost-MEDIAN_SVA_pairs2.rnk")
createrankedlist(Msva, "bothcohorts_preVpost-MEAN_SVA_pairs3.rnk")
createrankedlist(Msva, "bothcohorts_preVpost-MEDIAN_SVA_pairs3.rnk")

# VS//
# TAKING SINGLE PRE vs SINGLE POST
dim(pre_pairs)
for(i in 1:length(colnames(pre_pairs))){
  preM <- as.matrix(pre_pairs[,i])
  postM <- as.matrix(post_pairs[,i])
  M <- log2(preM/postM)
  #M <- log2(pre_pairs[,i]/post_pairs[,i])
  sample <- colnames(pre_pairs)[i]
  filename <- paste("bothcohorts_preVpost_", sample, ".rnk")
  createrankedlist(M, filename)
}








# PCA plot comparison batch correction
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#sedata
#combat
datpca <- prcomp(t(sedata), cor=F)
datpca <- prcomp(t(combat), cor=F)

plot(
  range(datpca$x[, 1]),range(datpca$x[, 2]), 
  type = "n", xlab = "pre", ylab = "post", 
  main = "PCA of colon trial data\npre vs post"
)

names(datpca)
head(phen)
points(datpca$x[, 1][phen$group == "pre"], datpca$x[, 2][phen$group == "pre"], bg = "Red", pch = 21)
text(datpca$x[, 1][phen$group == "pre"], datpca$x[, 2][phen$group == "pre"], label=names(datpca$x[,1][phen$group == "pre"]),pos=1,cex=0.5)
points(datpca$x[, 1][phen$group == "post"], datpca$x[, 2][phen$group == "post"], bg = "Blue", pch=21)
text(datpca$x[, 1][phen$group == "post"], datpca$x[, 2][phen$group == "post"], label=names(datpca$x[,1][phen$group == "post"]),pos=1,cex=0.5)
leg.names <- c("pre", "post")       
legend(-30,12,leg.names,leg.col,pch=15,cex=.7,horiz=F)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>








# Outliers? 
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
library(ggplot2)
rt <- t(sedata3)    	#transpose dat
rt.dist <- dist(rt,method="euclidean")	# calculate distance
rt.clust <- hclust(rt.dist,method="single")	# calculate clusters
plot(rt.clust,labels=names(rt),cex=0.75)	# plot cluster tree

r.mean <- apply(log2(sedata3),2,mean)  	# calculate mean for each sample
r.sd <- sqrt(apply(log2(sedata3),2,var))		# calculate st.deviation for each sample
r.cv <- r.sd/r.mean		

plot(r.mean,r.cv,main="colon trial dataset\nSample CV vs. Mean",xlab="Mean",ylab="CV",col='blue',cex=1.5,type="n")
points(r.mean,r.cv,bg="lightblue",col=1,pch=21)
text(r.mean,r.cv,label=dimnames(sedata3)[[2]],pos=1,cex=0.5)

# Hierarchical clustering - Pearson's correlation
dat.cor <- cor(sedata3, use = "pairwise.complete.obs", method = "pearson")
par(mar = rep(4, 4))
plot(dat.cor)

#hierarchical clustering plot
dat.t <- t(sedata3)
# get pairwise distances
dat.dist <- dist(dat.t)
# calculate and plot hierarchical tree
par(mar = rep(4, 4))
plot(hclust(dat.dist), labels=dimnames(sedata3)[[2]],main="HC Dendogram of Pre and Post Samples",cex.main=0.75)
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>