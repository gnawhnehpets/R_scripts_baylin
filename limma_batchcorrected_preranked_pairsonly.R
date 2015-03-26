#FINAL
library(limma)
setwd("/home/steve/Desktop/analysis/cohortoneandtwo//")
target <- readTargets("/home/steve/Desktop/analysis/cohortoneandtwo//targets_full.txt")
cohorttwo<- c(targets$Cy3_Sample, targets$Cy5_sample)
RG <- read.maimages(targets, source="agilent", path="/home/steve/Desktop/analysis/cohortoneandtwo/")
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
#normbetween <- normbetween[normbetween$genes$ControlType==0,]

#Convert MA back to RG
#ma2rg <- RG.MA(nc_normbetween)
RGb <- RG.MA(normbetween)

######################################
# sva ################################
######################################

# For all samples
# cy3
cy3 = RGb$R
rownames(cy3) <- RGb$genes$GeneName
colnames(cy3) <- targets$Cy3_Sample
colnames(cy3)
# cy5
cy5 = RGb$G
rownames(cy5) <- RGb$genes$GeneName
colnames(cy5) <- targets$Cy5_sample
colnames(cy5)

s <- intersect(targets$Cy3_Sample, targets$Cy5_sample)
s
length(s)
m_cy3 <- match(colnames(cy3), s)
m_cy3
m_cy5 <- match(colnames(cy5), s)
m_cy5

position_cy3 <- which(!is.na(m_cy3))
position_cy3
position_cy5 <- which(!is.na(m_cy5))
position_cy5

cy3_pairs <- cy3[,position_cy3]
dim(cy3_pairs)
cy5_pairs <- cy5[,position_cy5]
dim(cy5_pairs)

colnames(cy3_pairs)
colnames(cy5_pairs)

cy3_num <- as.character(targets$numpre[position_cy3])
cy5_num <- as.character(targets$numpost[position_cy5])
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

cy3names = as.matrix(paste(targets$Cy3[position_cy3], targets$Cy3_Sample[position_cy3], sep="-"))
cy5names = as.matrix(paste(targets$Cy5[position_cy5], targets$Cy5_sample[position_cy5], sep="-"))
cy3names
cy5names

# one matrix (edata)
# with "group-sample" naming format
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sva_pairs <- cbind(cy3_pairs, cy5_pairs)
dim(sva_pairs)
colnames(sva_pairs)
# with "group" naming format
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
colnames(cy3_pairs) <- targets$Cy3[position_cy3]
colnames(cy5_pairs) <- targets$Cy5[position_cy5]
sva_pairs2 <- cbind(cy3_pairs, cy5_pairs)
dim(sva_pairs2)
colnames(sva_pairs2)

sedata <- apply(sva_pairs,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
sedata2 <- apply(sva_pairs2,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 

# BATCHES
# create single matrix for date/location/rin
date <- as.data.frame(rbind(cy3date, cy5date))
colnames(date) <- "date"
dim(date)
location <- as.data.frame(rbind(cy3location, cy5location))
colnames(location) <- "location"
dim(location)
rin <- as.data.frame(rbind(cy3rin, cy5rin))
colnames(rin) <- "rin"
dim(rin)

group_cy3 = as.data.frame(colnames(cy3_pairs))
names(group_cy3) <- "group"
group_cy5 = as.data.frame(colnames(cy5_pairs))
names(group_cy5) <- "group"
cbind(group_cy3, group_cy5) #each row should have both values
group <- rbind(group_cy3, group_cy5)
#sample <- as.data.frame(rbind(cy3names, cy5names))
sample <- as.data.frame(c(1:length(rownames(group))))
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
phen <- cbind(sample, group, location, rincutoff)
phen
phen[,4] <- factor(phen[,4])
phen
colnames(sedata) <- c(cy3names, cy5names)
rownames(phen) <- colnames(sedata)
phen

library(sva)

# Create a model matrix for adjustment variables and variable of interest (rincutoff)
# http://www.bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf page 4
# mod = model.matrix(~as.factor(cancer), data=pheno)
# mod0 = model.matrix(~1,data=pheno)
smod = model.matrix(~as.factor(group), data=phen)
smod0 = model.matrix(~1, data=phen)

# This returns an expression matrix, with the same dimensions as original dataset. This new expression
# matrix is adjusted for batch.
# combat_edata = ComBat(dat=edata, batch=batch, mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE)
combat = ComBat(dat=sedata, batch=rincutoff, mod=smod, numCovs=NULL, par.prior=TRUE, prior.plots=TRUE)




pValuesComBat = f.pvalue(combat,smod,smod0)
qValuesComBat = p.adjust(pValuesComBat,method="BH")

# combat = new dat

#cy3_sva <- combat[,1:16]
num_of_samples <- length(colnames(combat))
cy3_sva <- combat[,1:(num_of_samples/2)]
colnames(cy3_sva) <- targets$Cy3
#cy5_sva <- combat[,17:32]
cy5_sva <- combat[,((num_of_samples/2)+1):num_of_samples]
colnames(cy5_sva) <- targets$Cy5
targets$Cy3
targets$Cy5

# NOTE: IMPORTANT! Cy3==PRE @position1 will be Cy5==POST
pre_sva <- cbind(cy3_sva[,colnames(cy3_sva)=="pre"], cy5_sva[,colnames(cy5_sva)=="pre"])
dim(pre)
post_sva <- cbind(cy5_sva[,colnames(cy5_sva)=="post"], cy3_sva[,colnames(cy3_sva)=="post"])
dim(post)
pre_names <- c(targets$Cy3_Sample[colnames(cy3_sva)=="pre"], targets$Cy5_sample[colnames(cy5_sva)=="pre"])
pre_names
post_names <- c(targets$Cy3_Sample[colnames(cy5_sva)=="post"], targets$Cy5_sample[colnames(cy3_sva)=="post"])
post_names

#Make colnames "group-sample"
colnames(pre_sva) <- paste(colnames(pre_sva),pre_names,sep="-")
colnames(pre_sva)
colnames(post_sva) <- paste(colnames(post_sva),post_names,sep="-")
colnames(post_sva)

presva_mean <- as.data.frame(apply(pre_sva,1,mean))
rownames(presva_mean) <- rownames(pre_sva)
postsva_mean <- as.data.frame(apply(post_sva,1,mean))
rownames(postsva_mean) <- rownames(post_sva)
Msva <- log2(presva_mean/postsva_mean)


#Create ranked list files, using M-values of comparison groups
createrankedlist <- function(x,y){
  genename <- as.data.frame(rownames(x))
  names(genename) <- "genename"
  new.dat <- cbind(genename,x)
  filename = y
  write.table(new.dat, file=y, quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")
}

createrankedlist(Msva, "bothcohorts_preVpost_BATCH-RIN_allsamples.rnk")
>>>>>>>>>>>>>>>
  
  ###################################################################
# PAIRS ONLY ######################################################
###################################################################

#colnames(pre_sva) <- paste(colnames(pre_sva),pre_names,sep="-")
#colnames(post_sva) <- paste(colnames(post_sva),post_names,sep="-")
colnames(pre_sva) <- pre_names
colnames(post_sva) <- post_names
colnames(pre_sva)
colnames(post_sva)

#Only keep pairs
s <- intersect(pre_names, post_names)
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

colnames(pre_pairs)
colnames(post_pairs)

# QC: check if unique rownames in pre_agg is the same as unique rownames in pre
length(unique(rownames(pre_sva)))
length(unique(rownames(pre)))
length(unique(rownames(pre_pairs)))
length(unique(rownames(post_sva)))
length(unique(rownames(post)))
length(unique(rownames(post_pairs)))

# TAKING AVERAGE OF ALL PREs vs AVERAGE OF ALL POSTs
# Create single expression values from multiple samples (paired)
pre_mean <- as.data.frame(apply(pre_pairs,1,mean))
rownames(pre_mean) <- rownames(pre_pairs)
post_mean <- as.data.frame(apply(post_pairs,1,mean))
rownames(post_mean) <- rownames(post_pairs)
M <- log2(pre_mean/post_mean)

#Create ranked list files, using M-values of comparison groups
createrankedlist <- function(x,y){
  genename <- as.data.frame(rownames(x))
  names(genename) <- "genename"
  new.dat <- cbind(genename,x)
  filename = y
  write.table(new.dat, file=y, quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")
}

createrankedlist(Msva, "bothcohorts_preVpost_BATCH-RIN_pairsonly.rnk")

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
points(datpca$x[, 1][phen$group == "pre"], datpca$x[, 2][phen$group == "pre"], bg = "Red", pch = 21)
text(datpca$x[, 1][phen$group == "pre"], datpca$x[, 2][phen$group == "pre"], label=names(datpca$x[,1][phen$group == "pre"]),pos=1,cex=0.55)
points(datpca$x[, 1][phen$group == "post"], datpca$x[, 2][phen$group == "post"], bg = "Blue", pch=21)
text(datpca$x[, 1][phen$group == "post"], datpca$x[, 2][phen$group == "post"], label=names(datpca$x[,1][phen$group == "post"]),pos=1,cex=0.55)
#leg.names <- c("pre", "post")       
#legend(-50,20,leg.names,leg.col,pch=15,cex=.7,horiz=F)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>








# Outliers? 
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
library(ggplot2)
rt <- t(sedata3)      #transpose dat
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













#>>>>>>>>>>>>>>..
#http://rafalab.jhsph.edu/batch/

library(corpcor)
library(affy)
source("http://rafalab.jhsph.edu/batch/myplclust.R")

##read-in array data
tab=read.csv("bladdercels/bladdertab.csv",as.is=TRUE)
eset=justRMA(filenames=tab$filename,celfile.path="bladdercels/")


###definde classes
outcome=tab[,8]
bt=tab[,5]
Index=which(outcome=="sTCC")
Cplus=grep("CIS",bt[Index])
outcome[Index]="sTCC-CIS"
outcome[Index[Cplus]]="sTCC+CIS"
outcome[49:57]<-"Biopsy"

##get expression matrix
mat=exprs(eset)

###get date
dates=vector("character",ncol(mat))
for(i in seq(along=dates)){
  tmp=affyio::read.celfile.header(file.path("bladdercels",tab$filenam[i]),info="full")$DatHeader
  dates[i]=strsplit(tmp,"\ +")[[1]][8]
}
dates=as.Date(dates,"%m/%d/%Y")

##divide dates into batches
batch=dates-min(dates);
batch=as.numeric(cut(batch,c(-1,10,75,200,300,500)))

##compute distance between and perform clustering 
mydist=dist(t(mat))
hc=hclust(mydist)
##make cluster. show outcome in text, batch in color
myplclust(hc,lab=outcome,lab.col=batch)
##one can also use muli-dimensional scaling
cmd=cmdscale(mydist)
plot(cmd,type="n")
text(cmd,outcome,col=batch)
##note the normals separate by date


##obtain singular value decomp 
s=fast.svd(mat-rowMeans(mat))

##how much variability explained by each component
plot(s$d^2/sum(s$d^2))
abline(h=0.10)
##note first two component explain almost half variability
##what do the correlate with?

##note correlation with batch
boxplot(split(s$v[,1],batch))
boxplot(split(s$v[,2],batch))

###and confounding between batch and outcome
table(outcome,batch)

###this amount of confounding is hard to fix.
###But we can try using sva or combat.
##sva
##http://www.biostat.jhsph.edu/~jleek/sva/
###or combat 
##http://jlab.byu.edu//ComBat/Download.html
