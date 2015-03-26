rm(list=ls())
for(i in 1:100){gc()}
library(minfi)
# cluster dir
setwd("/home/jhmi/shwang/proj/methylation_ex/idats")
baseDir = "/home/jhmi/shwang/proj/methylation_ex/idats"

# READ IN
targets <- read.450k.sheet(baseDir, recursive=FALSE)
targets$Basename <- paste(baseDir,"/",targets$Slide, "/", targets$Slide, "_", targets$Array, sep="")
RGset <- read.450k.exp(base = baseDir, targets = targets)
# load("/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/RGset")
# RGset <- load(file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/RGsetpairs")
load("/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/RGset")
pd <- pData(RGset) #phenodat

# PRE-SWAN
# MSet.raw <- preprocessRaw(RGset)
# MSet.norm <- preprocessIllumina(RGset, bg.correct = TRUE, normalize = "controls", reference = 2)
# MSet.bg <- bgcorrect.illumina(RGset)
# save(MSet.raw, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/MSet.raw")
# save(MSet.norm, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/MSet.norm")
# save(MSet.bg, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/MSet.bg")
# load(file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/MSet.raw")
# load(file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/MSet.norm")
# load(file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/MSet.bg")

# SWAN
# Mset.swan <- preprocessSWAN(RGset, MSet.raw) #original method
# Mset.swan2 <- preprocessSWAN(MSet.bg) #bg correction
# Mset.swan3 <- preprocessSWAN(RGset, MSet.norm) #ash's method
# save(Mset.swan, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/Mset.swan")
# save(Mset.swan2, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/Mset.swan2")
# save(Mset.swan3, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/Mset.swan3")
load(file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/Mset.swan")
# load(file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/Mset.swan2")
# load(file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/Mset.swan3")
# main points:
# 

#dat <- Mset.swan[,-c(1:8, 13:16, 24:27, 29:32,37:40, 42:43, 45:50,54:55,58:59)]
#dat2 <- cbind(Mset.swan[,c(1:8, 13:16, 24:27, 29:32,37:40, 42:43, 45:50,54:55,58:59)], Mset.swan[,-c(1:8, 13:16, 24:27, 29:32,37:40, 42:43, 45:50,54:55,58:59)])

# POST-SWAN - Method 1 #################################################################
# filter betas
# gb <- getBeta(Mset.swan)
# ga <- getAnnotation(Mset.swan)
# gb <- getBeta(Mset.swan)
# ga <- getAnnotation(Mset.swan)
# save(gb, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/gb")
# save(ga, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/ga")
load(file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/gb")
load(file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/ga")
colnames(gb) <- targets$Sample_Name
# all samples
gbpre <- gb[,grep("_pre",colnames(gb))]
gbpost <- gb[,grep("_post",colnames(gb))]
# paired samples only
gbpre <- gb[,c(2,4,6,8,14,16,25,27,30,32,38,40,43,46,48,50,55,59)]
gbpost <- gb[,c(1,3,5,7,13,15,24,26,29,31,37,39,42,45,47,49,54,58)]

# Individual diff methy probes/genes ####################################################################################################
ind <- 18
indexi <- which(gbpre[,ind] > .3)
deltarawi <- gbpost[,ind]-gbpre[, ind]
deltai <- gbpost[indexi,ind]-gbpre[indexi, ind]
sigprobesi <- names(which(deltai < -.2))
sigstartingbi <- gbpre[sigprobesi,ind]
index3i <- rownames(ga) %in% sigprobesi
which(index3i == "TRUE")
siggenesi <- ga$UCSC_RefGene_Name[which(index3i=="TRUE")] #annotation genenames of significant probes
uniquesiggenesi <- unique(siggenesi) #filter out copies
uniquesigi <- sort(unique(unlist(strsplit(uniquesiggenesi, ";"))))
length(uniquesigi)

#uniqsigpair1 <- uniquesigi
#uniqsigpair2 <- uniquesigi
#uniqsigpair3 <- uniquesigi
#uniqsigpair4 <- uniquesigi
#uniqsigpair5 <- uniquesigi
#uniqsigpair6 <- uniquesigi
#uniqsigpair7 <- uniquesigi
# uniqsigpair8 <- uniquesigi
# uniqsigpair9 <- uniquesigi
# uniqsigpair10 <- uniquesigi
# uniqsigpair11 <- uniquesigi
# uniqsigpair12 <- uniquesigi
# uniqsigpair13 <- uniquesigi
# uniqsigpair14 <- uniquesigi
# uniqsigpair15 <- uniquesigi
# uniqsigpair16 <- uniquesigi
# uniqsigpair17 <- uniquesigi
uniqsigpair18 <- uniquesigi

save(uniqsigpair1, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/uniqsigpair1")
save(uniqsigpair2, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/uniqsigpair2")
save(uniqsigpair3, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/uniqsigpair3")
save(uniqsigpair4, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/uniqsigpair4")
save(uniqsigpair5, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/uniqsigpair5")
save(uniqsigpair6, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/uniqsigpair6")
save(uniqsigpair7, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/uniqsigpair7")
save(uniqsigpair8, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/uniqsigpair8")
save(uniqsigpair9, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/uniqsigpair9")
save(uniqsigpair10, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/uniqsigpair10")
save(uniqsigpair11, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/uniqsigpair11")
save(uniqsigpair12, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/uniqsigpair12")
save(uniqsigpair13, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/uniqsigpair13")
save(uniqsigpair14, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/uniqsigpair14")
save(uniqsigpair15, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/uniqsigpair15")
save(uniqsigpair16, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/uniqsigpair16")
save(uniqsigpair17, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/uniqsigpair17")
save(uniqsigpair18, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/uniqsigpair1")


intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(uniqsigpair1,
          uniqsigpair2),
          uniqsigpair3),
          uniqsigpair4),
          uniqsigpair5),
          uniqsigpair6),
          uniqsigpair7),
          uniqsigpair8),
          uniqsigpair9),
          uniqsigpair10),
          uniqsigpair11),
          uniqsigpair12),
          uniqsigpair13),
          uniqsigpair14),
          uniqsigpair15),
          uniqsigpair16),
          uniqsigpair17),
          uniqsigpair18          
)


#########################################################################################################################################

premed <- as.matrix(rowMedians(gbpre))
postmed <- as.matrix(rowMedians(gbpost))
rownames(premed) <- rownames(gb)
rownames(postmed) <- rownames(gb)

# Finding which probes are linked to specific gene ####################################################################################################
colnames(Mset.swan) <- targets$Sample_Name
anndat <- cbind(ga$Name, ga$UCSC_RefGene_Name)
tmp <- setNames(strsplit(anndat[,2], split = ';'), anndat[,1])
tmp2 <- data.frame(probe = rep(names(tmp), sapply(tmp, length)), genename = unlist(tmp), row.names = NULL)
uniqtmp <- unique(tmp2)
uniqtmp
# index of dat associated with annotated genes
which(uniqtmp[,2] == "TET2") DNMT3A, ASXL1, RUNX1, NRAS
which(uniqtmp[,2] == "DNMT3A")
which(uniqtmp[,2] == "ASXL1")
which(uniqtmp[,2] == "RUNX1")
which(uniqtmp[,2] == "NRAS")
index <- which(uniqtmp[,2] %in% c("ARPC5", "CDKN2A", "TFPI2", "IGFBP3", "SFRP1", "GATA4", "GATA5", "SCGB3A1", "MGMT", "APC", "CHFR", "RASSF1", "MLH1"))
filteredtmp <- uniqtmp[index,]
geneprobes <- filteredtmp[,1]
geneprobes

index2 <- which(ga$Name %in% geneprobes)
# filteredcolonsamples <- colonsamples[geneprobes,]
# filteredcolonsamples <- Mset.swan[index2,]
mayofilteredprobes <- Mset.swan[index2,]
mayofilteredannotation <- getAnnotation(mayofilteredprobes)
mayofilteredbetas <- getBeta(mayofilteredprobes)
mayofilteredannotation$UCSC_RefGene_Name
# Loading required package: IlluminaHumanMethylation450kanno.ilmn12.hg19
save(mayofilteredbetas, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/mayofilteredbetas")
save(mayofilteredannotation, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/mayofilteredannotation")

#######################################################################################################################################################

# filter probes with starting beta > .5
index <- which(premed > .3) #317367
index <- which(premed > .4) #295616
index <- which(premed > .5) #267405
# length(index)

# save(index, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/index")
# load(file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/index")
pre_fil <- premed[index,]
post_fil <- postmed[index,]

# filter for probes with deltabeta >.2
deltabeta <- post_fil - pre_fil
deltabetaraw <- postmed - premed
sigprobes <- names(which(abs(deltabeta) > .2))
#sigprobes <- names(which(deltabeta < -.2))


# DATA
# starting beta of significant probes
sigstartingb <- pre_fil[sigprobes]
# delta beta of significant probes
sigdeltabeta <- deltabeta[sigprobes]
# find gene names of significant probes
index3 <- rownames(ga) %in% sigprobes
which(index3 == "TRUE")
siggenes <- ga$UCSC_RefGene_Name[which(index3=="TRUE")] #annotation genenames of significant probes
uniquesiggenes <- unique(siggenes) #filter out copies
uniquesig <- sort(unique(unlist(strsplit(uniquesiggenes, ";"))))
length(uniquesig)

print(uniquesig, quote=FALSE, row.names=FALSE)
# genes3 <- uniquesig #715    #895
# length(genes3)
# genes4 <- uniquesig #625    #843
# length(genes4)
# genes5 <- uniquesig #420      #670
# length(genes5)
# genes3_abs <- uniquesig #715    #895
# length(genes3_abs)
# genes4_abs <- uniquesig #625    #843
# length(genes4_abs)
genes5_abs <- uniquesig #420      #670
length(genes5_abs)


# save(genes3, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/genes3")
# save(genes4, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/genes4")
# save(genes5, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/genes5")
save(genes3_abs, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/genes3_abs")
save(genes4_abs, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/genes4_abs")
save(genes5_abs, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/genes5_abs")

load(file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/genes3")
load(file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/genes4")
load(file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/genes5")

which(genes3 %in% genes4 == FALSE)
genes.4 <- genes3[which(genes3 %in% genes4 == FALSE)]
print(genes.4, quote=FALSE, row.names=FALSE)

which(genes3 %in% genes5 == FALSE)
genes.5 <- genes3[which(genes3 %in% genes5 == FALSE)]
print(genes.5, quote=FALSE, row.names=FALSE)

which(genes4 %in% genes5 == FALSE)
genes.45 <- genes4[which(genes4 %in% genes5 == FALSE)]
print(genes.45, quote=FALSE, row.names=FALSE)


outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}
print(outersect(genes5, genes5_abs), quote=FALSE)
# sigprobes <- names(which(abs(deltabeta) > .2))
#.3 - 780
#.4 - 680
#.5 - 461
# sigprobes <- names(which(deltabeta < -.2))
#.3 - 724
#.4 - 626
#.5 - 423

getwd()
aim <- as.matrix(read.table("aimgenes.txt", header=FALSE))
cta <- as.matrix(read.table("CTAs.txt", header=FALSE))
hervs <- as.matrix(read.table("hervs.txt", header=FALSE))
intersect(aim, genes5)
intersect(aim, genes4)
intersect(aim, genes3)

# POST-SWAN - Method 2 #################################################################
# Finding differentially methylated positions (DMPs)
# Categorical phenotypes - 'dmpFinder' uses F-test to identify positions that are differentially methylated between two (or more) groups

# Find differences between GroupA and GroupB
table(pd$Sample_Group)



# minfi

mset <- Mset.swan
colnames(mset) <- targets$Sample_Name
M <- getM(mset, type = "beta", betaThreshold = 0.0001)
dmpnew <- dmpFinder(M, pheno=pd$Sample_Group, type="categorical")

save(dmp, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/dmp")
####################################################################
#individual comparison#
# WORK IN PROGRESS
# swanpre <- mset[,c(2,4,6,8,14,16,25,27,30,32,38,40,43,46,48,50,55,59)]
# swanpost <- mset[,c(1,3,5,7,13,15,24,26,29,31,37,39,42,45,47,49,54,58)]
line1 <- c(2,4,6,8,14,16,25,27,30,32,38,40,43,46,48,50,55,59)
line2 <- c(1,3,5,7,13,15,24,26,29,31,37,39,42,45,47,49,54,58)
line3 <- NULL
for(i in 1:18){
  line3 <- c(line3, line1[i], line2[i])
}
mset2 <- Mset.swan[,line3]
colnames(mset2) <- pd$Sample_Name[line3]
Mo <- getM(mset2, type = "beta", betaThreshold = 0.0001)
colnames(Mo)
ind <- c(1:4)
M <- Mo[,ind]
dmpnew <- dmpFinder(M, pheno=pd$Sample_Group[ind], type="categorical", shrinkVar=TRUE)

####################################################################
#M <- getM(mset, type = "beta", betaThreshold = 0.5)
# M
# Returns a table of CpG positions sorted by differential methylation p-value
# Tests each genomic position for assocaition between methylation and a phenotype
#load(file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/dmp")
dmp <- dmpnew
# List of most differentially methylated probes
dmp_names <- rownames(dmp)[dmp$qval < .05] # dmp_names <- rownames(dmp)[dmp$qval < .01] # dmp_names <- rownames(dmp)[dmp$qval < .0001]
length(dmp_names) #17395, 20773
# index of most dmp
index <- rownames(ga) %in% dmp_names
head(rownames(ga)[index]) #6
# genename of most dmp
dmp_genename <- ga$UCSC_RefGene_Name[index]
uniquegenes <- unique(dmp_genename)
uniquegenes <- unique(unlist(strsplit(uniquegenes, ";"))) #intersect with AIM gene list      #6159 unique genes, 7055
length(uniquegenes)


length(intersect(uniquegenes, gene5))
#print(sort(intersect(uniquegenes, genes3)), quote=FALSE)
print(sort(intersect(uniquegenes, genes5)), quote=FALSE)
data.frame(sort(intersect(uniquegenes, genes5))) #USE THIS OUTPUT WITH getHTMLsourceCode2.pl TO GET GENE ANNOTATION
save(uniquegenes, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/uniquegenes")
load(file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/uniquegenes")

# FILTER MINFI DMP FOR SB > .3, DB > .2
# filter only for dmp significant probes
pre_dmp <- as.matrix(premed[dmp_names,])
post_dmp <- as.matrix(postmed[dmp_names,])
# filter only for dmp significant probes that have starting beta of at least .3
index_dmp <- which(pre_dmp>.3)
fil_pre_dmp <- as.matrix(pre_dmp[index_dmp,])
fil_post_dmp <- as.matrix(post_dmp[index_dmp,])
length(fil_pre_dmp)
# filter only for dmp significant probes that have starting beta of at least .3 AND deltabeta of at least .2
deltabetadmp <- fil_post_dmp - fil_pre_dmp
sigprobes_dmp <- rownames(deltabetadmp)[which(abs(deltabetadmp) > .2)]
sigprobes_dmp <- rownames(deltabetadmp)[which(deltabetadmp < -.2)]
# find index in original annotation for significant probes
index_dmpga <- rownames(ga) %in% sigprobes_dmp
siggenes_dmp <- ga$UCSC_RefGene_Name[index_dmpga]
uniquesiggenes_dmp <- unique(siggenes_dmp) #filter out copies
uniquesig_dmp <- sort(unique(unlist(strsplit(uniquesiggenes_dmp, ";"))))
length(uniquesig_dmp)
save(uniquesig_dmp2, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/uniquesig_dmp2")

print(uniquesig_dmp, quote=FALSE, row.names=FALSE)



# CROSS WITH AIM GENE LIST
aim <- as.matrix(read.table("aimgenes.txt", header=FALSE)) # oncotarget AIM gene list
# Are there dmp genes that are also AIM genes?
aim_in_dmpgenes <- intersect(uniquegenes, aim)d
length(aim_in_dmpgenes)
#68 .001
#35 .0001
# http://stackoverflow.com/questions/24920807/how-can-i-maintain-relationship-within-a-matrix-after-splitting-a-character-vect
table <- cbind(rownames(ga), ga$UCSC_RefGene_Name)
tmp <- setNames(strsplit(ga$UCSC_RefGene_Name, ";"), table[,1])
table2 <- data.frame(probe = rep(names(tmp), sapply(tmp, length)), genename = unlist(tmp), row.names = NULL)
# head(table2, n=50)
index3 <- table2[,2] %in% aim_in_dmpgenes
table_of_aim_probes <- table2[index3,]

load(file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/Mset.swan")
ga <- getAnnotation(Mset.swan)

index <- which(table2[,1] %in% ga$Name)
column1 <- rep("NA", length(ga$Name))
# grab whatever information you want from annotation
probes <- ga$Name
genes <- ga$UCSC_RefGene_Name
accession <- ga$UCSC_RefGene_Accession
# create matrix of new annotation matrix
mat <- cbind(probes, genes, accession)

column2 <- rep("NA", length(table2[,1]))
column3 <- rep("NA", length(table2[,1]))
column4 <- rep("NA", length(table2[,1]))
newmat <- cbind(as.matrix(table2[,1]), column2, column3, column4)
colnames(newmat) <- c("probe", "probes", "genes", "accession")
head(newmat)

tmp <- NULL
for(i in 1:length(rownames(table2))){
#for(i in 1:5){
#  print(i)
  index <- which(mat[,1] %in% newmat[i,1])
#   newmat[i,2] <- mat[index, 1]
#   newmat[i,3] <- mat[index, 2]
#   newmat[i,4] <- mat[index, 3]
   newmat[i,2] <- ga$Name[i]
   newmat[i,3] <- ga$UCSC_RefGene_Name[i]
   newmat[i,4] <- ga$UCSC_RefGene_Accession[i]
}
new_ga <- ga[index,]
# index <- which(ga$Name %in% table2[,1])
# index2 <- duplicated(which(table2[,1] %in% ga$Name))

















###########################################################################################################

# If you doing a linear model fit, p-value corrected for multiple hyp can be used.
# Else deltabeta of >= 0.2 for probes that have starting beta (untreated) of > 0.5  is a safe start.

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
intersect(uniquesig, aim_in_dmpgenes)
intersect(uniquesig, uniquegenes)

table <- cbind(rownames(ga), ga$UCSC_RefGene_Name)
tmp <- setNames(strsplit(ga$UCSC_RefGene_Name, ";"), table[,1])
# http://stackoverflow.com/questions/24920807/how-can-i-maintain-relationship-within-a-matrix-after-splitting-a-character-vect
table2 <- data.frame(probe = rep(names(tmp), sapply(tmp, length)), genename = unlist(tmp), row.names = NULL)
# head(table2, n=50)
index3 <- table2[,2] %in% aim_in_dmpgenes
table_of_aim_probes <- table2[index3,]
head(table_of_aim_probes, n=50)

metsamp <- read.table("methylationsamplenames.txt", header=FALSE)
head(metsamp)
metsamp <- metsamp[,2]
metsamp <- gsub("_", "-", metsamp)
metsamp
meta$patient
which(meta$patient %in% metsamp)

##################################################################################################################

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
ratioSetSwan <- ratioConvert(dat2, what = "both", keepCN = TRUE)
# save(ratioSetSwan, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/ratioSetSwan")
# load(file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/ratioSetSwan")
# Associate genomic coordinates tot he probes using mapToGenome: transform ratioSet to a GenomicRatioSet (hold M and/or Beta values together with associated genomic coordinates)
# gRatioSet <- mapToGenome(ratioSet, mergeManifest = TRUE)
gRatioSetSwan <- mapToGenome(ratioSetSwan, mergeManifest = TRUE)
# save(gRatioSetSwan, file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/gRatioSetSwan")
# load(file="/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/gRatioSetSwan")
# showClass("GenomicRatioSet")
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




load("C:/Users/Steve/Downloads/wenbing.RData")
Mset.swan <- beta.tab

# BUMPHUNTER

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

Mset.swan
class(beta.tab)

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
