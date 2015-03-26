library(limma)
setwd("/home/steve/Desktop/analysis/08092013 analysis/")
targets <- readTargets("/home/steve/Desktop/analysis/08092013 analysis/targets_full.txt")
cohorttwo<- c(targets$Cy3_Sample, targets$Cy5_sample)
RG <- read.maimages(targets, source="agilent", path="/home/steve/Desktop/analysis/08092013 analysis/")
RG$weights <-  matrix(rep(RG$genes$ControlType,ncol(RG$R)),ncol=ncol(RG$R),byrow=F)
RG$weights[RG$genes$ControlType!=0,] <- 0
RG$weights[RG$genes$ControlType==0,] <- 1

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

# Collapse rows where genenames are the same
cy3_agg <- apply(cy3,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})}) 
cy5_agg <- apply(cy5,2,function(v){tapply(v,names(v),function(x){median(x,na.rm=TRUE)})})

# QC: check if unique rownames in pre_agg is the same as unique rownames in pre
length(unique(rownames(cy3_agg)))
length(unique(rownames(cy3)))
length(unique(rownames(cy5_agg)))
length(unique(rownames(cy5)))

dat <- cbind(cy3_agg, cy5_agg)
colnames(dat)
cls <- as.matrix(colnames(dat))
sample <- c(targets$Cy3_Sample, targets$Cy5_sample)
col <- paste(cls, sample, sep = '-')
col

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

create_gct(dat, "cohort2_expdat.gct", col)

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

create_cls(cls, "cohort2_expdat.cls")

###################################################################
# PAIRS ONLY ######################################################
###################################################################
# QC: check if unique rownames in pre_agg is the same as unique rownames in pre
length(unique(rownames(cy3_agg)))
length(unique(rownames(cy3)))
length(unique(rownames(cy5_agg)))
length(unique(rownames(cy5)))

colnames(cy3_agg)
colnames(cy5_agg)
targets$Cy3_Sample
targets$Cy5_sample
#=d
s <- intersect(targets$Cy3_Sample, targets$Cy5_sample)
s
m <- match(targets$Cy3_Sample, s)
m
position <- which(!is.na(m))
position
p1 <- cy3_agg[,position]
p2 <- cy5_agg[,position]
targets$Cy3_Sample[position]
targets$Cy5_sample[position]
colnames(p1)
colnames(p2)

dat <- cbind(p1, p2)
colnames(dat)
cls <- as.matrix(colnames(dat))
sample <- c(targets$Cy3_Sample[position], targets$Cy5_sample[position])
col <- paste(cls, sample, sep = '-')
col

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

create_gct(dat, "cohort2_expdat_pairsonly.gct", col)

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

create_cls(cls, "cohort2_expdat_pairsonly.cls")
