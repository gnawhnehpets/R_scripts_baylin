library(limma)
setwd("/home/steve/Desktop/MA_core_download/ahuja/expression/08092013 analysis/")
#Just pre-post, post-pre data
#targets <- targets[c(2,6:9),]
targets <- readTargets("/home/steve/Desktop/MA_core_download/ahuja/expression/08092013 analysis/targets_full.txt")
targets_fil <- targets[c(1:3,5:7,9:14,16),]
RG <- read.maimages(targets, source="agilent", path="/home/steve/Desktop/MA_core_download/ahuja/expression/08092013 analysis/")
RG_fil <- read.maimages(targets_fil, source = "agilent", path ="/home/steve/Desktop/MA_core_download/ahuja/expression/08092013 analysis/")

#RG
RG$weights <-  matrix(rep(RG$genes$ControlType,ncol(RG$R)),ncol=ncol(RG$R),byrow=F)
RG$weights[RG$genes$ControlType!=0,] <- 0
RG$weights[RG$genes$ControlType==0,] <- 1

colnames(RG) <- paste(targets$Cy5,targets$Cy5_sample)
colnames(RG$G) <- paste(targets$Cy3,targets$Cy3_Sample)
colnames(RG$R) <- paste(targets$Cy5,targets$Cy5_sample)
rownames(RG$R) <- paste(RG$genes$unique)
rownames(RG$G) <- paste(RG$genes$unique)

#RG$genes$unique <- c(1:nrow(RG))
#rownames(RG$genes) <- paste(RG$genes$unique)
MA <- normalizeWithinArrays(RG, method="loess", "normexp", offset=50)
raw.RG <- normalizeWithinArrays(RG, method="none", "normexp", offset=50)
#Remove controls
sub.MA <- MA[MA$genes$ControlType==0,]
MA.bet <- normalizeBetweenArrays(sub.MA, method="Aquantile")
dat <- sub.MA$M

#MA.bet <- normalizeBetweenArrays(dat, method="Aquantile")
###########################33
library(genefilter)
#rsd <- rowSds(dat.m)
rsd <- rowSds(dat)

#Create a vector of logical values (T or F)
i <- rsd >= 2
# Use the logical vector to filter the data matrix
dat.d <- dat[i,]
colnames(dat.d)
dim(dat.d)

### Expression filter

# When filtering by expression, it is not sensible to assume that all arrays behave similarly. On some arrays, 
# the gene might be expressed but, for some reason, not expressed at all or expressed at very low levels on other 
# arrays. The expression filter, therefore, has to account this possible discrepancy. This is done by letting the
# gene pass the filter (and be included in the dataset) if the gene is expressed at the set level in at least some
# proportion of the samples. These kind of filters can be easily created using the functions `kOverA()` and 
# `pOverA()`. 
# `kOverA()` uses the absolute number of samples during filtering
# `pOverA()` uses the proportion of samples during filtering

ff <- pOverA(A=1, p=0.5)
i <- genefilter(dat, ff)
dat.fo <- dat[i,] #overexpressed genes
dim(dat.fo)

# We now have a set of over-expressed genes (dat.fo). If we also want to get under expressed genes, we can invert the
# matrix (i.e. make under-expressed values over-expressed and vice versa). The procedure is exactly the same as the one
# outlined above, but we add a minus sign in front of the name of the matrix, which inverts its values

ff <- pOverA(A=1, p=0.5)
i <- genefilter(-dat, ff)
dat.fu <- dat[i,] #underexpressed genes
dim(dat.fu)

# We can then combine these two matrices into one matrix, if we want to retain both under- and over-expressed genes in 
# the same dataset. This is accomplished via `cbind` to combine two matrices row-wise (i.e. add rows in 2nd matrix after
# rows of first matrix)
dat.f <- rbind(dat.fo, dat.fu)

###########################
design <- modelMatrix(targets, ref = "pre")
#fit <- lmFit(sub.MA, design)
#fit <- lmFit(dat.d, design)
fit <- lmFit(dat.f, design)
fit <- eBayes(fit)
topTable(fit, coef=1.4)

rn <- rn <- rownames(topTable(fit, coef=1.4, n=100))

#tt <- topTable(fit, coef=1.4, n=nrow(dat))
#rn <- rownames(tt)[tt$P.Value <= .005]

rn <- as.numeric(rn)

# Extract genes from original data. `dat.s` now stores the data for the genes that were selected to be differentially 
# expressed.
dat.s <- dat[rn,]



##########################################################################################################################
#Filtered MAs
##########################################################################################################################

#RG_fil
RG_fil$weights <-  matrix(rep(RG_fil$genes$ControlType,ncol(RG_fil$R)),ncol=ncol(RG_fil$R),byrow=F)
RG_fil$weights[RG_fil$genes$ControlType!=0,] <- 0
RG_fil$weights[RG_fil$genes$ControlType==0,] <- 1

colnames(RG_fil) <- paste(targets_fil$Cy5,targets_fil$Cy5_sample)
colnames(RG_fil$G) <- paste(targets_fil$Cy3,targets_fil$Cy3_Sample)
colnames(RG_fil$R) <- paste(targets_fil$Cy5,targets_fil$Cy5_sample)
rownames(RG_fil$R) <- paste(RG_fil$genes$unique)
rownames(RG_fil$G) <- paste(RG_fil$genes$unique)

MA_fil <- normalizeWithinArrays(RG_fil, method="loess", "normexp", offset=50)
raw.RG_fil <- normalizeWithinArrays(RG_fil, method="none", "normexp", offset=50)
#Remove controls
sub.MA_fil <- MA_fil[MA_fil$genes$ControlType==0,]
MA.bet_fil <- normalizeBetweenArrays(sub.MA_fil, method="Aquantile")
dat_fil <- sub.MA_fil$M

#MA.bet <- normalizeBetweenArrays(dat, method="Aquantile")
###########################33
library(genefilter)
#rsd <- rowSds(dat.m)
rsd_fil <- rowSds(dat_fil)

#Create a vector of logical values (T or F)
i <- rsd_fil >= 2
# Use the logical vector to filter the data matrix
dat.d_fil <- dat_fil[i,]
colnames(dat.d_fil)


### Expression filter

ff_fil <- pOverA(A=1, p=0.5)
i_fil <- genefilter(dat_fil, ff_fil)
dat.fo_fil <- dat_fil[i_fil,] #overexpressed genes

ff_fil <- pOverA(A=1, p=0.5)
i_fil <- genefilter(-dat_fil, ff_fil)
dat.fu_fil <- dat_fil[i_fil,] #underexpressed genes

dat.f_fil <- rbind(dat.fo_fil, dat.fu_fil)

###########################
design_fil <- modelMatrix(targets_fil, ref = "pre")
#fit <- lmFit(sub.MA, design)
#fit <- lmFit(dat.d, design)
fit_fil <- lmFit(dat.f_fil, design_fil)
fit_fil <- eBayes(fit_fil)
topTable(fit_fil, coef=1.4)

rn_fil <- rn_fil <- rownames(topTable(fit_fil, coef=1.4, n=100))

#tt <- topTable(fit, coef=1.4, n=nrow(dat))
#rn <- rownames(tt)[tt$P.Value <= .005]

rn_fil <- as.numeric(rn_fil)

# Extract genes from original data. `dat.s` now stores the data for the genes that were selected to be differentially 
# expressed.
dat.s_fil <- dat_fil[rn_fil,]
