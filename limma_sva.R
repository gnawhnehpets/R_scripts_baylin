library(limma)
setwd("/home/steve/Desktop/analysis/08092013 analysis/")
targets <- readTargets("/home/steve/Desktop/analysis/08092013 analysis/targets_full.txt")
RG <- read.maimages(targets, source="agilent", path="/home/steve/Desktop/analysis/08092013 analysis/")
RG$weights <-  matrix(rep(RG$genes$ControlType,ncol(RG$R)),ncol=ncol(RG$R),byrow=F)
RG$weights[RG$genes$ControlType!=0,] <- 0
RG$weights[RG$genes$ControlType==0,] <- 1
dat <- RG[,-18]
dat <- dat[,-6]

#normalization
#normalization
#normwithin <-normalizeWithinArrays(RG,method='loess',bc.method='none')
normwithin <-normalizeWithinArrays(RG,method='loess',bc.method='normexp', offset=50)
normbetween <-normalizeBetweenArrays(normwithin,method='Aquantile')

#Remove controls from normwithin/between
normwithin <- normwithin[normbetween$genes$ControlType==0,]
normbetween <- normbetween[normbetween$genes$ControlType==0,]

#Convert MA back to RG
#ma2rg <- RG.MA(nc_normbetween)
RGb <- RG.MA(normbetween)

# cy3
cy3 = RGb$R
rownames(cy3) <- RGb$genes$GeneName
colnames(cy3) <- targets$Cy3
#cy3lognorm =normalizeBetweenArrays(log2(cy3), method="quantile") #logscaled Cy3 data
#cy3norm =normalizeBetweenArrays(cy3, method="quantile") #nonlogscaled Cy3 data

# cy5
cy5 = RGb$G
rownames(cy5) <- RGb$genes$GeneName
colnames(cy5) <- targets$Cy5

# compartmentalize channels into pre- and post-treatment groups

# ALL SAMPLES
pre_all = cbind(cy3[,c(2:4,8:10,12,14,15)],cy5[,c(1,4:8,11,13,15,16)])
post_all = cbind(cy3[,c(1,5:7,11,13,16)],cy5[,c(2,3,9,10,12,14)])
#PAIRS ONLY
#> targets$Cy3_Sample[c(1,2,5:7,9:14,16)]
#[1] "PH1876" "PH1900" "PH1824" "PH1910" "PH1612" "PH1861" "PH1616" "PH1844" "PH1641" "PH1843" "PH1871" "PH1869"
#> targets$Cy5_sample[c(1,2,5:7,9:14,16)]
#[1] "PH1876" "PH1900" "PH1824" "PH1910" "PH1612" "PH1861" "PH1616" "PH1844" "PH1641" "PH1843" "PH1871" "PH1869"
#> targets$Cy3[c(1,2,5:7,9:14,16)]
#[1] "post" "pre"  "post" "post" "post" "pre"  "pre"  "post" "pre"  "post" "pre"  "post"
#> targets$Cy5[c(1,2,5:7,9:14,16)]
#[1] "pre"  "post" "pre"  "pre"  "pre"  "post" "post" "pre"  "post" "pre"  "post" "pre" 

targets$Cy3_Sample[c(1,2,5:7,9:14,16)]
targets$Cy5_sample[c(1,2,5:7,9:14,16)]
targets$Cy3[c(1,2,5:7,9:14,16)]
targets$Cy5[c(1,2,5:7,9:14,16)]
cy3_pairs <- RGb$R[,c(1,2,5:7,9:14,16)]
cy5_pairs <- RGb$G[,c(1,2,5:7,9:14,16)]
colnames(cy3_pairs) <- targets$Cy3[c(1,2,5:7,9:14,16)]
colnames(cy5_pairs) <- targets$Cy5[c(1,2,5:7,9:14,16)]
pre = as.matrix(cbind(cy3_pairs[,c(2,6,7,9,11)],cy5_pairs[,c(1,3:5,8,10,12)]))
post = as.matrix(cbind(cy3_pairs[,c(1,3:5,8,10,12)],cy5_pairs[,c(2,6,7,9,11)]))
colnames(pre) <- targets$Cy3_Sample[c(1,2,5:7,9:14,16)]
colnames(post) <- targets$Cy5_sample[c(1,2,5:7,9:14,16)]
rownames(pre) <- RGb$genes$GeneName
rownames(post) <- RGb$genes$GeneName
#

# Collapse rows where genenames are the same
preall_agg <- apply(pre_all,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
postall_agg <- apply(post_all,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})})
pre_agg <- apply(pre,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
post_agg <- apply(post,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})})

# QC: check if unique rownames in pre_agg is the same as unique rownames in pre
length(unique(rownames(pre_agg)))
length(unique(rownames(pre)))
length(unique(rownames(preall_agg)))
length(unique(rownames(pre_all)))

#> dim(pre_agg)
#[1] 27193    17
#> dim(post_agg)
#[1] 27193     5

# Create single expression values from multiple samples (ALL)
preall_mean <- as.data.frame(apply(preall_agg,1,mean))
rownames(preall_mean) <- rownames(preall_agg)
postall_mean <- as.data.frame(apply(postall_agg,1,mean))
rownames(postall_mean) <- rownames(postall_agg)
Mall <- log2(preall_mean/postall_mean)

# Create single expression values from multiple samples (paired)
pre_mean <- as.data.frame(apply(pre_agg,1,mean))
rownames(pre_mean) <- rownames(pre_agg)
post_mean <- as.data.frame(apply(post_agg,1,mean))
rownames(post_mean) <- rownames(post_agg)
M <- log2(pre_mean/post_mean)

#Create ranked list files, using M-values of comparison groups
createrankedlist <- function(x,y){
  genename <- as.data.frame(rownames(x))
  names(genename) <- "genename"
  new.dat <- cbind(genename,x)
  filename = y
  write.table(new.dat, file=y, quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")
}

createrankedlist(Mall, "preVSpost_colondata_allsamples.rnk")
createrankedlist(M, "preVSpost_colondata_pairs_only.rnk")
#createrankedlist(Mcc1,"controlVScycle1.rnk")
#createrankedlist(Mcc2,"controlVScycle2.rnk")
#createrankedlist(Mcc3b,"controlVScycle3b.rnk")
#createrankedlist(Mcc3a,"controlVScycle3a.rnk")
#createrankedlist(Mc2c3b,"cycle2VScycle3b.rnk")
#createrankedlist(Mc2c3a,"cycle2VScycle3a.rnk")

######################################
# sva ################################
######################################
cy3names = paste(targets$Cy3, targets$Cy3_Sample, sep="-")
cy5names = paste(targets$Cy5, targets$Cy5_sample, sep="-")
# cy3
sva_cy3 = RGb$R
rownames(sva_cy3) <- RGb$genes$GeneName
colnames(sva_cy3) <- cy3names

# cy5
sva_cy5 = RGb$G
rownames(sva_cy5) <- RGb$genes$GeneName
colnames(sva_cy5) <- cy5names

# one matrix (edata)
sva <- cbind(sva_cy3, sva_cy5)
sedata <- apply(sva,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 

# batches
date <- as.data.frame(c(1,1,2,2,1,1,2,2,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,2,2))
colnames(date) <- "date"
location <- as.data.frame(c(4,7,6,5,7,7,5,2,4,4,6,6,1,1,6,6,5,5,4,4,4,4,6,6,3,3,4,4,7,7,6,6))
colnames(location) <- "location"
group_cy3 = as.data.frame(targets$Cy3)
names(group_cy3) <- "group"
group_cy5 = as.data.frame(targets$Cy5)
names(group_cy5) <- "group"
group <- rbind(group_cy3, group_cy5)
colnames(group) <- "group"
sample <- as.data.frame(c(1:32))
colnames(sample) <- "sample"

phen <- cbind(sample, group, location)
rownames(phen) <- colnames(sedata)


# Create a model matrix for adjustment variables and variable of interest (location)
# mod = model.matrix(~as.factor(cancer), data=pheno)
# mod0 = model.matrix(~1,data=pheno)
smod = model.matrix(~as.factor(location), data=phen)
smod0 = model.matrix(~1, data=phen)

# n.sv = num.sv(edata,mod,method="leek")
#sn.sv = num.sv(sedata,smod,method="leek")
#sn.sv

#  svobj = sva(edata,mod,mod0,n.sv=n.sv)
#ssvobj = sva(sedata, smod, smod0, n.sv=sn.sv)

# This returns an expression matrix, with the same dimensions as original dataset. This new expression
# matrix is adjusted for batch.
# combat_edata = ComBat(dat=edata, batch=batch, mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE)
combat = ComBat(dat=sedata, batch=group, mod=smod, numCovs=NULL, par.prior=TRUE, prior.plots=TRUE)

pValuesComBat = f.pvalue(combat,smod,smod0)
qValuesComBat = p.adjust(pValuesComBat,method="BH")

cy3_sva <- combat[,1:16]
colnames(cy3_sva) <- targets$Cy3
cy5_sva <- combat[,17:32]
colnames(cy5_sva) <- targets$Cy5
pre_sva = cbind(cy3_sva[,c(2:4,8:10,12,14,15)],cy5_sva[,c(1,4:8,11,13,15,16)])
post_sva = cbind(cy3_sva[,c(1,5:7,11,13,16)],cy5_sva[,c(2,3,9,10,12,14)])
colnames(pre_sva)
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

createrankedlist(Msva, "location-date-sva_preVSpost_colondata_allsamples.rnk")
0