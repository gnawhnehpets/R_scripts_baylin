#COHORT1
library(limma)
setwd("/home/steve/Desktop/analysis/04232013 analysis/")
targets <- readTargets("/home/steve/Desktop/analysis/04232013 analysis/targets_full.txt")
cohortone<- c(targets$Cy3_Sample, targets$Cy5_sample)
RG <- read.maimages(targets, source="agilent", path="/home/steve/Desktop/analysis/04232013 analysis/")
RG$weights <-  matrix(rep(RG$genes$ControlType,ncol(RG$R)),ncol=ncol(RG$R),byrow=F)
RG$weights[RG$genes$ControlType!=0,] <- 0
RG$weights[RG$genes$ControlType==0,] <- 1

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

#PAIRS
targets$Cy3_Sample
targets$Cy5_sample
intersect(targets$Cy3_Sample, targets$Cy5_sample)
#Next two should be identical to 'intersect'
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

# Create single expression values from multiple samples (paired)
pre_mean <- as.data.frame(apply(pre_agg,1,mean))
rownames(pre_mean) <- rownames(pre_agg)
post_mean <- as.data.frame(apply(post_agg,1,mean))
rownames(post_mean) <- rownames(post_agg)

#Create ranked list files, using M-values of comparison groups
createrankedlist <- function(x,y){
  genename <- as.data.frame(rownames(x))
  names(genename) <- "genename"
  new.dat <- cbind(genename,x)
  filename = y
  write.table(new.dat, file=y, quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")
}

createrankedlist(pre_mean, "pre_only-exp_val.rnk")
createrankedlist(post_mean, "post_only-exp_val.rnk")
