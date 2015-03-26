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
load("/home/jhmi/shwang/proj/methylation_ex/idats/Robjects/RGset")
pd <- pData(RGset) #phenodat

# PRE-SWAN
MSet.raw <- preprocessRaw(RGset)

# SWAN
Mset.swan <- preprocessSWAN(RGset, MSet.raw) #original method

# POST-SWAN - Method 1 #################################################################
# betas
gb <- getBeta(Mset.swan)
ga <- getAnnotation(Mset.swan)
# all samples
gbpre <- gb[,grep("_pre",colnames(gb))]
gbpost <- gb[,grep("_post",colnames(gb))]
# paired samples only
gbpre <- gb[,c(2,4,6,8,14,16,25,27,30,32,38,40,43,46,48,50,55,59)] #all pre samples
gbpost <- gb[,c(1,3,5,7,13,15,24,26,29,31,37,39,42,45,47,49,54,58)] #all post samples

premed <- as.matrix(rowMedians(gbpre))
postmed <- as.matrix(rowMedians(gbpost))
rownames(premed) <- rownames(gb)
rownames(postmed) <- rownames(gb)

# Finding which probes are linked to specific gene ####################
colnames(Mset.swan) <- targets$Sample_Name
anndat <- cbind(ga$Name, ga$UCSC_RefGene_Name)
tmp <- setNames(strsplit(anndat[,2], split = ';'), anndat[,1])
tmp2 <- data.frame(probe = rep(names(tmp), sapply(tmp, length)), genename = unlist(tmp), row.names = NULL)
uniqtmp <- unique(tmp2)
uniqtmp
# index of dat associated with annotated genes
which(uniqtmp[,2] == "TET2") # check to see if gene exists
index <- which(uniqtmp[,2] %in% c("ARPC5", "CDKN2A", "TFPI2"))
filteredtmp <- uniqtmp[index,]
geneprobes <- filteredtmp[,1]
geneprobes

index2 <- which(ga$Name %in% geneprobes)
mayofilteredprobes <- Mset.swan[index2,]
filtered_ga <- getAnnotation(mayofilteredprobes)
filtered_ga$UCSC_RefGene_Name
# Loading required package: IlluminaHumanMethylation450kanno.ilmn12.hg19

# Finding which probes are linked to specific gene ####################

# filter probes with starting beta > .5
index <- which(premed > .3) #317367
index <- which(premed > .4) #295616
index <- which(premed > .5) #267405
pre_fil <- premed[index,]
post_fil <- postmed[index,]

# filter for probes with deltabeta >.2
deltabeta <- post_fil - pre_fil
deltabetaraw <- postmed - premed
sigprobes <- names(which(deltabeta < -.2))

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
genes3 <- uniquesig #715    #895
genes4 <- uniquesig #625    #843
genes5 <- uniquesig #420      #670

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
mset <- Mset.swan
# colnames(mset) <- targets$Sample_Name
M <- getM(mset, type = "beta", betaThreshold = 0.0001)
dmp <- dmpFinder(M, pheno=pd$Sample_Group, type="categorical")
# M
# Returns a table of CpG positions sorted by differential methylation p-value
# Tests each genomic position for assocaition between methylation and a phenotype
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

# intersect minfi with hari's method (starting beta .5)
length(intersect(uniquegenes, gene5))
#print(sort(intersect(uniquegenes, genes3)), quote=FALSE)
print(sort(intersect(uniquegenes, genes5)), quote=FALSE)
data.frame(sort(intersect(uniquegenes, genes5))) #USE THIS OUTPUT WITH getHTMLsourceCode2.pl TO GET GENE ANNOTATION

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
table <- cbind(rownames(ga), ga$UCSC_RefGene_Name)
tmp <- setNames(strsplit(ga$UCSC_RefGene_Name, ";"), table[,1])
# http://stackoverflow.com/questions/24920807/how-can-i-maintain-relationship-within-a-matrix-after-splitting-a-character-vect
table2 <- data.frame(probe = rep(names(tmp), sapply(tmp, length)), genename = unlist(tmp), row.names = NULL)
# head(table2, n=50)
index3 <- table2[,2] %in% aim_in_dmpgenes
table_of_aim_probes <- table2[index3,]