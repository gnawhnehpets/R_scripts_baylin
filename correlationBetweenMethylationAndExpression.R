setwd("/home/jhmi/shwang/proj/methylation_ex/idats")
baseDir = "/home/jhmi/shwang/proj/methylation_ex/idats"
targets <- read.450k.sheet(baseDir, recursive=FALSE)

library(minfi)
library(limma)
#load expresion log2fc data, post/pre, for each gene
load("../../m084b/fil_pre_col")
load("../../m084b/fil_post_col")
load("../../m084b/fil_pre_med_col")
load("../../m084b/fil_post_med_col")
load("../../m084b/fil_genename") #genename

#create single value for each group
pre <- as.matrix(rowMedians(fil_pre_col))
post <- as.matrix(rowMedians(fil_post_col))
pre_med <- as.matrix(rowMedians(fil_pre_med_col))
post_med <- as.matrix(rowMedians(fil_post_med_col))
rownames(pre) <- genename
rownames(post) <- genename
rownames(pre_med) <- genename
rownames(post_med) <- genename

log2col <- log2(post/pre)
log2col1 <- log(post/pre)
log2medcol <- log2(post_med/pre_med)

#load delta beta values per gene
# load(file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/deltabeta")
# load(file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/deltabetaraw")

#correlate relationship between paired samples
########################################################################################
load(file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/Mset.swan")
load(file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/Mset.swan2")
load(file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/Mset.swan3")
# gb <- getBeta(Mset.swan)
# ga <- getAnnotation(Mset.swan)

#1: swan type (1-3)
#2: allsamples (1) or paired samples (2)
#3: choose delta: allprobes (1), significant (2), significant down (3), sig up (4)
gb <- getBeta(Mset.swan2)
ga <- getAnnotation(Mset.swan2)
colnames(gb) <- targets$Sample_Name

# all samples
gbpre <- gb[,grep("_pre",colnames(gb))]
gbpost <- gb[,grep("_post",colnames(gb))]
# paired samples only
gbpre <- gb[,c(2,4,6,8,14,16,25,27,30,32,38,40,43,46,48,50,55,59)]
gbpost <- gb[,c(1,3,5,7,13,15,24,26,29,31,37,39,42,45,47,49,54,58)]

premed <- as.matrix(rowMedians(gbpre))
postmed <- as.matrix(rowMedians(gbpost))
rownames(premed) <- rownames(gb)
rownames(postmed) <- rownames(gb)

# filter probes with starting beta > .5
index1 <- which(premed > .5)
pre_fil <- premed[index1,]
post_fil <- postmed[index1,]
# filter for probes with deltabeta >.2 ############
# deltabeta of probes with starting beta > .5
deltabeta <- post_fil - pre_fil
# deltabetaraw <- postmed - premed
# deltabeta of significant probes (startingbeta > .5, deltabeta > .2)
sigdeltabeta <- as.matrix(deltabeta[which(abs(deltabeta) > .2)])
# names of significant probes
sigprobes <- names(which(abs(deltabeta) > .2))
sigdeltabetagenes <- ga$UCSC_RefGene_Name[which(ga$Name %in% sigprobes)]
#index of significant probes from original dataset
index2 <- which(ga$Name %in% sigprobes)
# # filter probes 
# premed[sigprobes,]
# index2 <- which(rownames(premed) %in% sigprobes)
# premed[index2,]

########################################################################################
# load(file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/index")
# load(file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/ga")
deltabetagenes <- ga$UCSC_RefGene_Name[index1]
sigdeltabetagenes <- ga$UCSC_RefGene_Name[index2]
## CHOOSE ONE #########################################
#1. 
delta1 <- cbind(deltabeta, deltabetagenes)
index <- index1
delta <- delta1
#2.
delta2 <- cbind(sigdeltabeta, sigdeltabetagenes)
colnames(delta2) <- c("deltabeta", "deltabetagenes")
index <- index2
delta <- delta2

#3. deltabeta < -.2
downindex <- which(sigdeltabeta < -.2)
delta3 <- cbind(sigdeltabeta[downindex], sigdeltabetagenes[downindex])
colnames(delta2) <- c("deltabeta", "deltabetagenes")
delta <- delta3
index <- downindex

#4. deltabeta > .2
upindex <- which(sigdeltabeta > .2)
delta4 <- cbind(sigdeltabeta[upindex], sigdeltabetagenes[upindex])
colnames(delta2) <- c("deltabeta", "deltabetagenes")
delta <- delta4
index <- upindex
#######################################################
colnames(delta) <- c("deltabeta", "deltabetagenes")
#find genename based on indices of original 485k probes
# tmp <- setNames(strsplit(as.character(test.frame$amounts), split = ','), test.frame$name)
# data.frame(name = rep(names(tmp), sapply(tmp, length)), amounts = unlist(tmp), row.names = NULL)
# tmp <- setNames(strsplit(ga$UCSC_RefGene_Name, ";"), delta[,1])[index]
# tmp <- setNames(strsplit(ga$UCSC_RefGene_Name, ";")[index], delta[,1])
# tmp2 <- data.frame(probe = rep(names(tmp[index]), sapply(tmp, length)), genename = unlist(tmp), row.names = NULL)
# tmp2 <- as.matrix(tmp2)
# http://stackoverflow.com/questions/24249351/can-i-split-this-dataframe-up-so-that-theres-a-new-row-for-each-item-in-some-co
tmp <- setNames(strsplit(delta[,2], split = ';'), delta[,1])
tmp2 <- data.frame(probe = rep(names(tmp), sapply(tmp, length)), genename = unlist(tmp), row.names = NULL)


# #grab all instances where MMEL1 is in the genename
# delta[grep("TMEM52", delta[,2]),]

#uniqtmp2 <- as.matrix(unique(tmp2))
uniqtmp2 <- unique(tmp2)
#uniqtmp2[,1] <- as.numeric(uniqtmp2[,1])
tmpdat1 <- as.matrix(uniqtmp2[,1])
tmpdat2 <- as.matrix(uniqtmp2[,2])
tmpdat <- cbind(tmpdat1, tmpdat2)
colnames(tmpdat) <- c("probe", "genename")
# http://stackoverflow.com/questions/16777972/r-collapse-rows-and-sum-the-values-in-the-column
# aggregate(PSM ~ ID, data=x, FUN=sum)

## CHOOSE FUNCTION ####################################
#1.
tmpdat1 <- aggregate(as.numeric(tmpdat[,1]) ~ genename, data=tmpdat, FUN=median) #CHANGE FUNCTION
dim(tmpdat1) #9592
#######################################################
# intersect of expression data genes, infinium data genes
commongenes <- intersect(as.character(tmpdat1[,1]), rownames(log2col))
length(commongenes) #8892

## CHOOSE L2FC ########################################
#1.
l2fc <- as.matrix(log2col[commongenes,]) # dataframe #1
#2. 
l2fc <- as.matrix(log2medcol[commongenes,]) # dataframe #1
#######################################################

index4 <- which(tmpdat1[,1] %in% commongenes)
meth <- as.matrix(tmpdat1[index4,2]) # dataframe #2
rownames(meth) <- as.character(tmpdat1[index4,1])
length(intersect(rownames(meth), rownames(l2fc))) # checkpoint: should equal 5067

#alternative: indicates Ha; 'less' corresponds to negative association
ctest1 <- cor.test(meth, l2fc, method="spearman", conf.level=.95, alternative="less") 
ctest2 <- cor.test(meth, l2fc, method="spearman", conf.level=.95, alternative="greater") 
ctest3 <- cor.test(meth, l2fc, method="spearman", conf.level=.95, alternative="two.sided") 
print(paste(ctest1$p.value, ctest2$p.value, ctest3$p.value, sep="   "))

jpeg(filename="expressionVsMethylation_1-1-2.jpeg", height=600, width=1000)
plot(meth, l2fc, type="p", main="expression vs methylation", xlab="deltabeta", ylab="log2fc")
# Add fit lines
abline(lm(meth~l2fc), col="red") # regression line (y~x) 
lines(lowess(meth,l2fc), col="blue") # lowess line (x,y)
dev.off()

cor(cbind(meth,l2fc))
jpeg(filename="correlationBetweenMethylationAndExpression_1-3-2.jpeg", height=600, width=1000)
pairs(cbind(meth,l2fc))
dev.off()

# 
# aggregate function = mean
# Spearman's rank correlation rho
# 
# data:  meth and l2fc
# S = 22154356374, p-value = 0.06055
# alternative hypothesis: true rho is less than 0
# sample estimates:
#         rho 
# -0.02178069 

# aggregate function = median
# Spearman's rank correlation rho
# 
# data:  meth and l2fc
# S = 21952530118, p-value = 0.1874
# alternative hypothesis: true rho is less than 0
# sample estimates:
#         rho 
# -0.01247227 

# aggregate function = min