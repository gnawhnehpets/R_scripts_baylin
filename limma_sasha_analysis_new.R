library(limma)
setwd("/home/steve/Desktop/analysis/01142014 Alexandra Borodovsky Agilent GE/")
targets <- readTargets("/home/steve/Desktop/analysis/01142014 Alexandra Borodovsky Agilent GE/targets_full.txt")
RG <- read.maimages(targets, source="agilent", path="/home/steve/Desktop/analysis/01142014 Alexandra Borodovsky Agilent GE/")


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

#M-values for conditions 1 to 4
controlc1M <- normwithin$M[,2]
controlc2M <- normwithin$M[,3]
controlc3bM <- normwithin$M[,4]
controlc3aM <- normwithin$M[,1]

#controlc1M <- normbetween$M[,2]
#controlc2M <- normbetween$M[,3]
#controlc3bM <- normbetween$M[,4]
#controlc3aM <- normbetween$M[,1]

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

createrankedlist(Mcc1,"controlVScycle3.rnk")

test(Mcc1,"data2.rnk")

test <- function(x,y){
  write.table(x, file=y)
}

createrankedlist <- function(x,y){
  genename <- as.data.frame(RGb$genes$GeneName)
  names(genename) <- "genename"
  new.dat <- cbind(genename,x)
  # collapse similar probes frand get median
  agg <- aggregate(. ~ genename, data=new.dat, median)
  filename = y
  write.table(agg[,1:2], file=y, quote = FALSE, row.names = FALSE, col.names = FALSE )
}
















pre = cbind(cy3norm[,c(1:6,9:11)],cy5norm[,c(1,3:5,7,8,10,11)])
post = cbind(cy3norm[,c(7,8)],cy5norm[,c(2,6,9)])
pre_mean <- as.data.frame(apply(pre,1,mean))
rownames(pre_mean) <- rownames(pre)
post_mean <- as.data.frame(apply(post,1,mean))
rownames(post_mean) <- rownames(post)
M <- log2(post_mean/pre_mean)

probeid <- as.data.frame(MA$genes$ProbeName)
names(probeid) <- "probeid"
dat <- MA$M
colnames(dat) <- paste("sample", MA$targets$SampleNumber, sep="")
genename <- as.data.frame(MA$genes$GeneName)
names(genename) <- "genename"
new.dat <- cbind(probeid,genename,dat)

# collapse similar probes frand get median
#dat.med <- aggregate(. ~ genes, data = MA, median)
#agg1 <- aggregate(. ~ probeid, data=new.dat, median)
agg <- aggregate(. ~ genename, data=new.dat, median)


MA = normalizeWithinArrays(RGb, method="loess", weights=RGb$weights)
logratio <- MA$M # ratio
#The median of the replicate probes:
Cy3.norm <- apply(cy3,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
Cy5.norm <- apply(cy5,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})})
LogRatio <- apply(logratio,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})})

#E2.avg <- avereps(normwithin, ID=MA$genes$ProbeName)
#E2.avg_RG <- RG.MA(E2.avg)

#probeid <- as.data.frame(MA$genes$ProbeName)
#names(probeid) <- "probeid"



probeid <- as.data.frame(MA$genes$ProbeName)
names(probeid) <- "probeid"
dat <- MA$M
colnames(dat) <- paste("sample", MA$targets$SampleNumber, sep="")
genename <- as.data.frame(MA$genes$GeneName)
names(genename) <- "genename"
new.dat <- cbind(probeid,genename,dat)

# collapse similar probes frand get median
#dat.med <- aggregate(. ~ genes, data = MA, median)
#agg1 <- aggregate(. ~ probeid, data=new.dat, median)
agg <- aggregate(. ~ genename, data=new.dat, median)
