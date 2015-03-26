#FINAL
library(limma)
setwd("/home/steve/Desktop/analysis/cohortoneandtwo/")
targets <- readTargets("/home/steve/Desktop/analysis/cohortoneandtwo/targets_full.txt")
cohorttwo<- c(targets$Cy3_Sample, targets$Cy5_sample)
RG <- read.maimages(targets, source="agilent", path="/home/steve/Desktop/analysis/cohortoneandtwo/")
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

#sort into pre/post
colnames(cy3)
colnames(cy5)
pre <- cbind(cy3[,colnames(cy3)=="pre"], cy5[,colnames(cy5)=="pre"])
dim(pre)
post <- cbind(cy3[,colnames(cy3)=="post"], cy5[,colnames(cy5)=="post"])
dim(post)
pre_names <- c(targets$Cy3_Sample[colnames(cy3)=="pre"], targets$Cy5_sample[colnames(cy5)=="pre"])
pre_names
post_names <- c(targets$Cy3_Sample[colnames(cy3)=="post"], targets$Cy5_sample[colnames(cy5)=="post"])
post_names

#Make colnames "pre/post - sample"
colnames(pre) <- paste(colnames(pre),pre_names,sep="-")
colnames(pre)
colnames(post) <- paste(colnames(post),post_names,sep="-")
colnames(post)

# Collapse rows where genenames are the same
pre_agg <- apply(pre,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
post_agg <- apply(post,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})})

# QC: check if unique rownames in pre_agg is the same as unique rownames in pre
length(unique(rownames(pre_agg)))
length(unique(rownames(pre)))
length(unique(rownames(post_agg)))
length(unique(rownames(post)))

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

createrankedlist(M, "bothcohorts_preVpost.rnk")

###################################################################
# PAIRS ONLY ######################################################
###################################################################

#sort into pre/post
colnames(cy3)
colnames(cy5)
# NOTE: IMPORTANT! Cy3==PRE @position1 will be Cy5==POST
pre <- cbind(cy3[,colnames(cy3)=="pre"], cy5[,colnames(cy5)=="pre"])
dim(pre)
post <- cbind(cy5[,colnames(cy5)=="post"], cy3[,colnames(cy3)=="post"])
dim(post)
pre_names <- c(targets$Cy3_Sample[colnames(cy3)=="pre"], targets$Cy5_sample[colnames(cy5)=="pre"])
pre_names
post_names <- c(targets$Cy3_Sample[colnames(cy5)=="post"], targets$Cy5_sample[colnames(cy3)=="post"])
post_names

#Make colnames = sample name
colnames(pre) <- pre_names
colnames(pre)
colnames(post) <- post_names
colnames(post)

# Collapse rows where genenames are the same
pre_agg <- apply(pre,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
post_agg <- apply(post,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})})

#Only keep pairs
s <- intersect(pre_names, post_names)
s
m_pre <- match(colnames(pre_agg), s)
m_pre
m_post <- match(colnames(post_agg), s)
m_post

position_pre <- which(!is.na(m_pre))
position_pre
position_post <- which(!is.na(m_post))
position_post

pre_pairs <- pre_agg[,position_pre]
dim(pre_pairs)
post_pairs <- post_agg[,position_post]
dim(post_pairs)

colnames(pre_pairs)
colnames(post_pairs)

# QC: check if unique rownames in pre_agg is the same as unique rownames in pre
length(unique(rownames(pre_agg)))
length(unique(rownames(pre)))
length(unique(rownames(pre_pairs)))
length(unique(rownames(post_agg)))
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

createrankedlist(M, "bothcohorts_preVpost_pairsonly.rnk")

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















#######################################################
#######################################################
#NOT DONE

# Separate channel analysis of two color-data
#file1: control vs c3a
#file2: control vs c1
#file3: control vs c2
#file4: control vs c3b

# convert targets frame to channel (vs array-oriented)
targets2 <- targetsA2C(targets)

# create design matrix
u <- unique(targets2$Target)
f <- factor(targets2$Target, levels=u)
design <- model.matrix(~0+f)
colnames(design) <- u
library(statmod)
corfit <- intraspotCorrelation(normbetween, design)
fit <- lmscFit(normbetween, design, correlation=corfit$consensus)

#logFC: log2 fold change between groups; negative = down-regulation in second group, positve = up-regulation in second group
#t: moderated t-statistic from empirical Bayes method
#B: log-odds that gene is differentially expressed
#P.value: raw p-value
#adj.P.val: p.value corrected for multiple comparisons using Benjamini and Hochberg's FDR

# control v cycle1
cont.matrix <- makeContrasts("control-c1",levels=design)
con_v_c1 <- contrasts.fit(fit, cont.matrix)
con_v_c1 <- eBayes(con_v_c1)
controlvc1 <- topTable(con_v_c1, coef=1.4, adjust="BH", n=50)

# control v cycle2
cont.matrix <- makeContrasts("control-c2",levels=design)
con_v_c2 <- contrasts.fit(fit, cont.matrix)
con_v_c2 <- eBayes(con_v_c2)
controlvc2 <- topTable(con_v_c2, coef=1.4, adjust="BH", n=50)

# control v cycle3b
cont.matrix <- makeContrasts("control-c3w",levels=design)
con_v_c3w <- contrasts.fit(fit, cont.matrix)
con_v_c3w <- eBayes(con_v_c3w)
controlvc3w <- topTable(con_v_c3w, coef=1.4, adjust="BH", n=50)

# control v cycle3a
cont.matrix <- makeContrasts("control-c3c",levels=design)
con_v_c3c <- contrasts.fit(fit, cont.matrix)
con_v_c3c <- eBayes(con_v_c3c)
controlvc3c <- topTable(con_v_c3c, coef=1.4, adjust="BH", n=50)

# cycle2 v cycle3b
cont.matrix <- makeContrasts("c2-c3w",levels=design)
c2_v_c3w <- contrasts.fit(fit, cont.matrix)
c2_v_c3w <- eBayes(c2_v_c3w)
c2vc3w <- topTable(c2_v_c3w, coef=1.4, adjust="BH", n=50)

# cycle2 v cycle3a
cont.matrix <- makeContrasts("c2-c3c",levels=design)
c2_v_c3c <- contrasts.fit(fit, cont.matrix)
c2_v_c3c <- eBayes(c2_v_c3c)
c2vc3c <- topTable(c2_v_c3c, coef=1.4, adjust="BH", n=43119)

#################################################################
#M-values derived via individual channel analysis

#Convert MA back to RG
#ma2rg <- RG.MA(normbetween)
RGb <- RG.MA(normbetween)

# IGNORE
#RGb = backgroundCorrect(ma2rg, method="normexp", offset=50)

# cy3 = controls
cy3 = RGb$R
rownames(cy3) <- RGb$genes$GeneName
colnames(cy3) <- targets$Cy3
#cy3lognorm =normalizeBetweenArrays(log2(cy3), method="quantile") #logscaled Cy3 data
#cy3norm =normalizeBetweenArrays(cy3, method="quantile") #nonlogscaled Cy3 data

# cy5 = experimental group (cycle1, cycle2, cycle3b, cycle3a)
cy5 = RGb$G
rownames(cy5) <- RGb$genes$GeneName
colnames(cy5) <- targets$Cy5
#cy5lognorm =normalizeBetweenArrays(log2(cy5), method="quantile") #logscaled Cy5 data
#cy5norm = normalizeBetweenArrays(cy5, method="quantile") #nonlogscaled Cy5 data

# M-values
# control vs cycle1 - 
Mcc1 = as.data.frame(log2(cy3norm[,2]/cy5norm[,2]))
#rownames(Mcc1) <- rownames(Mcc1)
colnames(Mcc1) <- "controlvs cycle1"
# control vs cycle2
Mcc2 = as.data.frame(log2(cy3norm[,3]/cy5norm[,3]))
#rownames(Mcc2) <- rownames(cy3norm)
colnames(Mcc2) <- "control vs cycle2"
# control vs cycle3b
Mcc3b =as.data.frame(log2(cy3norm[,4]/cy5norm[,4]))
#rownames(Mcc3b) <- rownames(cy3norm)
colnames(Mcc3b) <- "control vs cycle3b"
# control vs cycle3a
Mcc3a = as.data.frame(log2(cy3norm[,1]/cy5norm[,1]))
#rownames(Mcc3a) <- rownames(cy3norm)
colnames(Mcc3a) <- "control vs cycle3a"
# cycle2 vs cycle3b
Mc2c3b = as.data.frame(log2(cy5norm[,3]/cy5norm[,4]))
#rownames(Mc2c3b) <- rownames(cy3norm)
colnames(Mc2c3b) <- "cycle2 vs cycle3b"
# cycle2 vs cycle3a
Mc2c3a = as.data.frame(log2(cy5norm[,3]/cy5norm[,1]))
#rownames(Mc2c3a) <- rownames(cy3norm)
colnames(Mc2c3a) <- "cycle2 vs cycle3a"

#Create ranked list files, using M-values of comparison groups
createrankedlist <- function(x,y){
  genename <- as.data.frame(RGb$genes$GeneName)
  names(genename) <- "genename"
  new.dat <- cbind(genename,x)
  # collapse similar probes and get median
  agg <- aggregate(. ~ genename, data=new.dat, median)
  filename = y
  write.table(agg, file=y, quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")
}

createrankedlist(Mcc1,"controlVScycle1.rnk")
createrankedlist(Mcc2,"controlVScycle2.rnk")
createrankedlist(Mcc3b,"controlVScycle3b.rnk")
createrankedlist(Mcc3a,"controlVScycle3a.rnk")
createrankedlist(Mc2c3b,"cycle2VScycle3b.rnk")
createrankedlist(Mc2c3a,"cycle2VScycle3a.rnk")