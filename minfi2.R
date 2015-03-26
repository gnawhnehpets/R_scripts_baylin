http://www.bioconductor.org/help/course-materials/2013/BioC2013/minfiLab.pdf
http://www.bioconductor.org/packages/release/bioc/vignettes/minfi/inst/doc/minfi.pdf

rm(list=ls())
for(i in 1:100){gc()}
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("minfiData")
library(minfi)
setwd("Y:/users/shwang26/methylation_ex/idats")
baseDir = "Y:/users/shwang26/methylation_ex/idats"
setwd("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/methylation_ex/idats")
baseDir = "/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/methylation_ex/idats"

#cluster 
setwd("/home/jhmi/shwang/proj/methylation_ex/idats")
baseDir = "/home/jhmi/shwang/proj/methylation_ex/idats"

targets <- read.450k.sheet(baseDir, recursive=FALSE)
targets$Basename <- paste(baseDir,"/",targets$Slide, "/", targets$Slide, "_", targets$Array, sep="")
pairtargets <- targets[c(1:8,13:16,24:27,29:32,37:40,42,43,45:50,54,55,58,59),]

RGset <- read.450k.exp(base = baseDir, targets = targets)
RGsetpairs <- read.450k.exp(base = baseDir, targets = pairtargets)
# save(RGset, file="Y:/users/shwang26/methylation_ex/idats/R objects/RGset")
save(RGset, file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/RGset")
# load("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/methylation_ex/idats/R objects/RGset")
# load("Y:/users/shwang26/methylation_ex/idats/R objects/RGset")
load("/home/jhmi/shwang/proj/methylation_ex/idats/R objects/RGset")
# RGset <- load(file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/RGsetpairs")
pd <- pData(RGset) #phenodat

# manifest <- getManifest(RGset)
# manifest
# probeInfo <- getProbeInfo(manifest)
# head(probeInfo)
# typeIProbes <- getProbeInfo(manifest, type = "I")$Name
# typeIIProbes <- getProbeInfo(manifest, type = "II")$Name

# QC
qcReport(RGset, sampNames = pd$Sample_Name, sampGroups = pd$Sample_Group, pdf = "qcReport_prenorm.pdf")

# Density plots of methylation Beta plots
jpeg(filename = "densityPlot.jpeg", height = 900, width = 1600)
densityPlot(RGset, sampGroups = pd$Sample_Group,main = "Beta", xlab = "Beta")
dev.off()

# Density bean plots of methylation Beta values
jpeg(filename = "densityBeanPlot.jpeg", height = 500, width = 500)
par(oma=c(2,10,1,1))
densityBeanPlot(RGset, sampGroups = pd$Sample_Group, sampNames = pd$Sample_Name)
dev.off()

# Preprocessing (normalization)
# takes input a RGChannelSet --> outputs MethylSet
# number of preprocessing options, each one implemented as preprocessXXX where XXX = name of method
# "Raw" method = convert Red and Green channel into Methylated and Unmethylated signal

# MSet.raw <- preprocessRaw(RGset)
# save(MSet.raw, file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/MSet.raw")
load(file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/MSet.raw")
# load("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/methylation_ex/idats/R objects/MSet.raw")
# load("Y:/users/shwang26/methylation_ex/idats/R objects/MSet.raw")

# qc <- getQC(MSet.raw)
# jpeg(filename = "plotQC.jpeg", height = 450, width = 800)
# plotQC(qc)
# dev.off()

# PRE-SWAN
# "Normal" method = background normalization, control normalization
# MSet.norm <- preprocessIllumina(RGset, bg.correct = TRUE, normalize = "controls", reference = 2)
# save(MSet.norm, file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/MSet.norm")
load(file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/MSet.norm")
# MSet.normbg <- bgcorrect.illumina(RGset)
#save(MSet.normbg, file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/MSet.normbg")
#load(file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/MSet.normbg")


jpeg(filename = "densityPlotMSetnorm.jpeg", height = 900, width = 1600)
densityPlot(MSet.norm, sampGroups = pd$Sample_Group,main = "Beta", xlab = "Beta")
dev.off()

jpeg(filename = "controlStripPlot.jpeg", height = 900, width = 1600)
controlStripPlot(RGset, controls="BISULFITE CONVERSION II", sampNames = pd$Sample_Name)
dev.off()

jpeg(filename = "mdsPlotMSetnorm.jpeg", height = 400, width = 600)
mdsPlot(MSet.norm, numPositions = 1000, sampGroups = pd$Sample_Group, sampNames = pd$Sample_Name)
dev.off()

MsetEx <- MSet.raw
save(MsetEx, file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/MsetEx")
load(file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/MsetEx")

RGsetEx <- RGset
save(RGsetEx, file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/RGsetEx")
load(file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/RGsetEx")

# http://www.bioconductor.org/help/course-materials/2013/BioC2013/minfiLab.pdf
gRatioSet.quantile <- preprocessQuantile(RGset, fixOutliers = TRUE,
                                         removeBadSamples = TRUE, badSampleCutoff = 10.5,
                                         quantileNormalize = TRUE, stratified = TRUE,
                                         mergeManifest = FALSE, sex = NULL)

save(gRatioSet.quantile, file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/gRatioSet.quantile")
load(file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/gRatioSet.quantile")

#Mset.swan <- preprocessSWAN(RGsetEx, MsetEx)
#Mset.swan2 <- preprocessSWAN(MSet.normbg)

#save(Mset.swan, file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/Mset.swan")
load(file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/Mset.swan")
#save(Mset.swan2, file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/Mset.swan2")
#load(file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/Mset.swan2")


jpeg(filename = "densityPlotSwan.jpeg", height = 900, width = 1600)
densityPlot(Mset.swan, sampGroups = pd$Sample_Group,main = "Beta", xlab = "Beta")
dev.off()

jpeg(filename = "mdsPlotMSetSwan.jpeg", height = 400, width = 600)
mdsPlot(Mset.swan, numPositions = 1000, sampGroups = pd$Sample_Group, sampNames = pd$Sample_Name)
dev.off()


#####################
# After preprocessing
#####################


# returns beta values betwen 0 and 1 where 1 is very high methylation; if type = "Illumina", then formula:
# beta = methylated / (methylated + unmeth + 100)
gb <- getBeta(Mset.swan)
colnames(gb) <- targets$Sample_Name
# gbmed <- as.matrix(rowMedians(gb))
# #gbmed <- as.matrix(rowMeans(gb))
# rownames(gbmed) <- rownames(gb)

#pairs
gbpre <- gb[,c(2,4,6,8,14,16,25,27,30,32,38,40,43,46,48,50,55,59)]
gbpost <- gb[,c(1,3,5,7,13,15,24,26,29,31,37,39,42,45,47,49,54,58)]
#all pres
gbpre <- gb[,grep("_pre",colnames(gb))]
gbpost <- gb[,grep("_post",colnames(gb))]

premed <- as.matrix(rowMedians(gbpre))
postmed <- as.matrix(rowMedians(gbpost))
premed2 <- as.matrix(rowMeans(gbpre))
postmed2 <- as.matrix(rowMeans(gbpost))
rownames(premed) <- rownames(gb)
rownames(postmed) <- rownames(gb)
index <- which(premed > .5) #all:266755, pairs: 269723
save(index, file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/index")

pre_fil <- premed[index,]
post_fil <- postmed[index,]
# rownames(premed)[index]
# rownames(postmed)[index]
deltabeta <- post_fil - pre_fil
deltabetaraw <- postmed - premed
#save(deltabeta, file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/deltabeta")
# load(file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/deltabeta")
#save(deltabetaraw, file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/deltabetaraw")
# load(file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/deltabetaraw")
sigprobes <- names(which(abs(deltabeta) > .2))
# starting beta of significant probes
sigstartingb <- pre_fil[sigprobes]
# delta beta of significant probes
sigdeltabeta <- deltabeta[sigprobes]
# > head(data.frame(cbind(sigstartingb, sigdeltabeta)), n=100)
# sigstartingb sigdeltabeta
# cg00650640    0.6533278   -0.3238181
# cg02285579    0.6077444   -0.2556652
# cg05064925    0.5649316   -0.2201263
# cg05135828    0.6199000   -0.2405886
# cg05747459    0.5124683   -0.2268940
# cg05788681    0.5408125   -0.2505591
# cg06650546    0.5036181   -0.2133168
# cg06779458    0.6101787   -0.2167926
# cg09196081    0.5461061   -0.2121622
# cg09687738    0.5019990   -0.2526744
# cg09725213    0.6104378   -0.2487565
# cg12105899    0.6891505   -0.2202986


index3 <- rownames(ga) %in% sigprobes
which(index3 == "TRUE")
#find probes at index
siggenes <- ga$UCSC_RefGene_Name[which(index3=="TRUE")] #genenames of significant probes

uniquesiggenes <- unique(siggenes)
uniquesig <- sort(unique(unlist(strsplit(uniquesiggenes, ";")))) #46 unique genes
length(uniquesig) #705

aim <- as.matrix(read.table("aimgenes.txt", header=FALSE))
intersect(aim, uniquesig)

print(intersect(uniquesig, aim_in_dmpgenes), quote=FALSE)
print(intersect(uniquesig, uniquegenes), quote=FALSE)


data.frame(colnames(gbpre), colnames(gbpost))
# Make betas.txt for COHCAP
# header <- t(as.matrix(c("Site ID", colnames(gb))))
# betamat <- cbind(rownames(gb), gb)
# head(betamat)
# beta <- rbind(header, betamat)
# head(col)
# write.table(beta, file="betas.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

colnames(Mset.swan) = pd$Sample_Name
colnames(RGset)
jpeg(filename="plotBetasbyType_comparison1.jpeg", height=900, width=1600)
par(mfrow=c(1,2))
# Distribution of betas for each sample
plotBetasByType(MsetEx[,1], main = "Raw")
plotBetasByType(Mset.swan[,1], main = "SWAN")
dev.off()

jpeg(filename="plotBetasbyType_swan1.jpeg", height=900, width=1600)
par(mfrow=c(3,4))
for(i in 1:12){
  plotBetasByType(Mset.swan[,i], main = colnames(Mset.swan)[i])
}
dev.off()

jpeg(filename="plotBetasbyType_swan2.jpeg", height=900, width=1600)
par(mfrow=c(3,4))
for(i in 13:24){
  plotBetasByType(Mset.swan[,i], main = colnames(Mset.swan)[i])
}
dev.off()

jpeg(filename="plotBetasbyType_swan3.jpeg", height=900, width=1600)
par(mfrow=c(3,4))
for(i in 25:36){
  plotBetasByType(Mset.swan[,i], main = colnames(Mset.swan)[i])
}
dev.off()

jpeg(filename="plotBetasbyType_swan4.jpeg", height=900, width=1600)
par(mfrow=c(3,4))
for(i in 37:48){
  plotBetasByType(Mset.swan[,i], main = colnames(Mset.swan)[i])
}
dev.off()

jpeg(filename="plotBetasbyType_swan5.jpeg", height=900, width=1600)
par(mfrow=c(3,4))
for(i in 49:60){
  plotBetasByType(Mset.swan[,i], main = colnames(Mset.swan)[i])
}
dev.off()


# plotBetasByType(Mset.swan[,1], main = colnames(Mset.swan)[1])
# plotBetasByType(Mset.swan[,2], main = colnames(Mset.swan)[2])
# colnames(Mset.swan)
# rownames(pd)
# data.frame(colnames(Mset.swan), rownames(pd))
# dev.off()

# Store beta values and/or M-values instead of the methylayed and unmethylated signals
# ratioSet <- ratioConvert(MSet.raw, what = "both", keepCN = TRUE)
ratioSetSwan <- ratioConvert(Mset.swan, what = "both", keepCN = TRUE)
save(ratioSetSwan, file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/ratioSetSwan")
load(file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/ratioSetSwan")
# Associate genomic coordinates tot he probes using mapToGenome: transform ratioSet to a GenomicRatioSet (hold M and/or Beta values together with associated genomic coordinates)
# gRatioSet <- mapToGenome(ratioSet, mergeManifest = TRUE)
gRatioSetSwan <- mapToGenome(ratioSetSwan, mergeManifest = TRUE)
save(gRatioSetSwan, file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/gRatioSetSwan")
load(file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/gRatioSetSwan")
showClass("GenomicRatioSet")
gbs <- getBeta(gRatioSetSwan)
# getM(gRatioSet)
# getCN(gRatioSet)
# sampleNames <- sampleNames(gRatioSetSwan)
# probeNames <- featureNames(gRatioSetSwan)
# pheno <- pData(gRatioSetSwan)
ga <- getAnnotation(gRatioSetSwan)
# gl <- getLocations(gRatioSetSwan)
# gsi <- getSnpInfo(gRatioSetSwan)
# gpt <- getProbeType(gRatioSetSwan)

# Island<-as.character(Anotation[which(Anotation$Relation_to_UCSC_CpG_Island=='Island'),1])
island<-as.character(ga$Name[which(ga$Relation_to_Island=='Island')])
head(island)
length(island)

# Annotation to sea/shore/island/shelf/etc
# gr$Relation_to_Island

########## onccbio #############################################################################
# CpG island probes
colnames(gbs) <- targets$Sample_Name

colo.beta.island <- gbs[island,]
jpeg(filename="coloBetaIsland_boxplot.jpeg", height=900, width=1600)
boxplot(colo.beta.island, las=2)
dev.off()

# CpG islands that are methylated
# breast.beta.Island.methylated<-breast.beta.Island[which(apply(breast.beta.Island,1,max)>=0.5),]
colo.beta.island.methylated <- colo.beta.island[which(apply(colo.beta.island,1,max)>=0.5),]
head(colo.beta.island.methylated)
dim(colo.beta.island.methylated)

# Promoter<-as.character(Anotation[grep('Promoter',Anotation[,32]),1])
names(ga)
unique(ga$Regulatory_Feature_Group)
promoter <-as.character(ga$Name[grep('Promoter',ga$Regulatory_Feature_Group)])
head(promoter)
length(promoter)
island.promoter<-intersect(island,promoter)
head(island.promoter)
length(island.promoter)

colo.beta.island.promoter<-gbs[island.promoter,]
dim(colo.beta.island.promoter)
jpeg(filename="coloBetaIslandPromoter_boxplot.jpeg", height=900, width=1600)
boxplot(colo.beta.island.promoter,las=2)
dev.off()

colo.beta.island.promoter.methylated<-colo.beta.island.promoter[which(apply(colo.beta.island.promoter,1,max)>=0.5),]
dim(colo.beta.island.promoter.methylated)
jpeg(filename="coloBetaIslandPromoterMethylated_boxplot.jpeg", height=900, width=1600)
boxplot(colo.beta.island.promoter.methylated, las=2)
dev.off()


nshore <- as.character(ga$Name[which(ga$Relation_to_Island=='N_Shore')])
colo.beta.nshore <- gbs[nshore,]
head(colo.beta.nshore)
dim(colo.beta.nshore)
jpeg(filename="coloBetaNShore_boxplot.jpeg", height=900, width=1600)
boxplot(colo.beta.nshore, las=2)
dev.off()

sshore <- as.character(ga$Name[which(ga$Relation_to_Island=='S_Shore')])
colo.beta.sshore <- gbs[sshore,]
head(colo.beta.sshore)
dim(colo.beta.sshore)
jpeg(filename="coloBetaSShore_boxplot.jpeg", height=900, width=1600)
boxplot(colo.beta.sshore, las=2)
dev.off()

sshelf <- as.character(ga$Name[which(ga$Relation_to_Island=='S_Shelf')])
colo.beta.sshelf <- gbs[sshelf,]
head(colo.beta.sshelf)
dim(colo.beta.sshelf)
jpeg(filename="coloBetaSShelf_boxplot.jpeg", height=900, width=1600)
boxplot(colo.beta.sshelf, las=2)
dev.off()

nshelf <- as.character(ga$Name[which(ga$Relation_to_Island=='N_Shelf')])
colo.beta.nshelf <- gbs[nshelf,]
head(colo.beta.nshelf)
dim(colo.beta.nshelf)
jpeg(filename="coloBetaNShelf_boxplot.jpeg", height=900, width=1600)
boxplot(colo.beta.nshelf, las=2)
dev.off()

# save probe name
originalprobelist <- rownames(ga)

# remove sex chromosome
rownames(ga) <- c(1:length(rownames(ga)))
sex.probe<-rownames(ga)[which(ga$chr=='chrX'|ga$chr=='chrY')]
length(sex.probe)
rownames(ga)


rownames(ga)[which(ga$UCSC_RefGene_Name=='ADM')]

adm <-gbs[rownames(ga)[which(ga$UCSC_RefGene_Name=='ADM')]]

################################################################################################
# To extract the genomic locations
gRanges <- rowData(gRatioSetSwan)
save(gRanges, file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/gRanges")
load(file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/gRanges")
head(gRanges, n= 3)

# Finding differentially methylated positions (DMPs)
# Categorical phenotypes - 'dmpFinder' uses F-test to identify positions that are differentially methylated between two (or more) groups

# Find differences between GroupA and GroupB
table(pd$Sample_Group)
mset <- Mset.swan
M <- getM(mset, type = "beta", betaThreshold = 0.001)
# M
# Returns a table of CpG positions sorted by differential methylation p-value
# Tests each genomic position for assocaition between methylation and a phenotype
dmp <- dmpFinder(M, pheno=pd$Sample_Group, type="categorical")
save(dmp, file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/dmp")
load(file="/home/jhmi/shwang/proj/methylation_ex/idats/R objects/dmp")
head(dmp)

# List of most differentially methylated probes
dmp_names <- rownames(dmp)[dmp$pval < .001]
length(dmp_names) #17395
# index of most dmp
index <- rownames(ga) %in% dmp_names
head(rownames(ga)[index])
# genename of most dmp
dmp_genename <- ga$UCSC_RefGene_Name[index]
uniquegenes <- unique(dmp_genename) #intersect with AIM gene list
uniquegenes <- unique(unlist(strsplit(uniquegenes, ";"))) #6159 unique genes


aim <- as.matrix(read.table("aimgenes.txt", header=FALSE))
# Are there dmp genes that are also AIM genes?
aim_in_dmpgenes <- intersect(uniquegenes, aim)
length(aim_in_dmpgenes)
#68 .001
#35 .0001
table <- cbind(rownames(ga), ga$UCSC_RefGene_Name)
tmp <- setNames(strsplit(ga$UCSC_RefGene_Name, ";"), table[,1])
# http://stackoverflow.com/questions/24920807/how-can-i-maintain-relationship-within-a-matrix-after-splitting-a-character-vect
table2 <- data.frame(probe = rep(names(tmp), sapply(tmp, length)), genename = unlist(tmp), row.names = NULL)
# head(table2, n=50)
index3 <- table2[,2] %in% aim_in_dmpgenes
table_of_aim_probes <- table2[index3,]
head(table_of_aim_probes, n=50)

intersect(uniquesig, aim_in_dmpgenes)
intersect(uniquesig, uniquegenes)


# plotCpG function - plot methylation levels at individual positions
jpeg(filename="plotCpGindividualPositions.jpeg", height=900, width=1900)
cpgs <- rownames(dmp)[1:75]
length(rownames(dmp))
par(mfrow=c(5,15))
plotCpg(mset, cpg=cpgs, pheno=pd$Sample_Group)
dev.off()



# # Continuous phenotypes - identify DMPs where mean methylation level varies with a continuous covariate using linear regression
# continuousPheno <- rnorm(nrow(pd))
# continuousPheno
# 
# dmp <- dmpFinder(mset, pheno=continuousPheno, type="continuous")
# # beta column gives change in mean phenotype for each unit increase of methylation. We now filter the DMP list to exclude positions with small effect size
# dmp[1:3,]
# dmp <- subset(dmp, abs(beta)>1)
# 
# # Visualize continuous DMPs
# jpeg(filename="plotCpGindividualPositions-continuous2.jpeg", height=900, width=1600)
# cpgs <- rownames(dmp)[1:4]
# par(mfrow=c(2,2))
# plotCpg(mset, cpg=cpgs, type="continuous", pheno=continuousPheno, xlab="Phenotype 1")
# dev.off()




# beta <- getBeta(gRatioSet.quantile)
# names(pData(gRatioSet.quantile))
# prepost <- pData(gRatioSet.quantile)$Sample_Group
# dmp <- dmpFinder(beta, pheno = prepost , type = "categorical")
# head(dmp)
# 
# pheno <- pData(gRatioSet.quantile)$Sample_Group
# designMatrix <- model.matrix(~ pheno)
# dmrs <- bumphunter(gRatioSet.quantile, design = designMatrix, cutoff = 0.05, B=1)

# CONVERT Mset.swan (methylset) to a genomicRatioSet
test <- Mset.swan
ratioSet <- ratioConvert(Mset.swan, what = "both", keepCN = TRUE)
gRatioSet <- mapToGenome(ratioSet, mergeManifest = TRUE)

beta <- getBeta(gRatioSet)
names(pData(gRatioSet))
prepost <- pData(gRatioSet)$Sample_Group
dmp <- dmpFinder(beta, pheno = prepost , type = "categorical")
head(dmp)

pheno <- pData(gRatioSet)$Sample_Group
designMatrix <- model.matrix(~ pheno)
dmrs <- bumphunter(gRatioSet, design = designMatrix, cutoff = 0.05, B=1)
# save(dmrs, file="./Robjects/dmrs")
load("./Robjects/dmrs")
dmrindex <- which(dmrs$table$p.value < .05)
dmrs2 <- dmrs$table
dmrindex2 <- which(dmrs2$p.value < .05)
dmrs3 <- dmrs2[dmrindex2,]




# colon biopsy: 