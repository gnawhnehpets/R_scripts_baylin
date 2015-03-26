rm(list=ls())
for(i in 1:100){gc()}
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("minfiData")
library(minfi)

#mex dir
setwd("D:/mex_data/new_data")
baseDir = "D:/mex_data/new_data"

targets <- read.450k.sheet(baseDir, recursive=FALSE)
targets
targets$Basename <- paste(baseDir,"/",targets$Slide, "/", targets$Slide, "_", targets$Array, sep="")
#pairtargets <- targets[c(1:8,13:16,24:27,29:32,37:40,42,43,45:50,54,55,58,59),]

#RGset <- read.450k.exp(base = baseDir, targets = targets[1:30,])
RGset <- read.450k.exp(base = baseDir, targets = targets)
RGsetpairs <- read.450k.exp(base = baseDir, targets = pairtargets)
# save(RGset, file="D:/mex_data/new_data/RGset")
load("D:/mex_data/new_data/RGset")
pd <- pData(RGset) #phenodat
# MSet.raw <- preprocessRaw(RGset)
# save(MSet.raw, file="D:/mex_data/new_data/MSet.raw")
load(file="D:/mex_data/new_data/MSet.raw")

# PRE-SWAN
# "Normal" method = background normalization, control normalization
# MSet.norm <- preprocessIllumina(RGset, bg.correct = TRUE, normalize = "controls", reference = 2)
#save(MSet.norm, file="D:/mex_data/new_data/MSet.norm")
load(file="D:/mex_data/new_data/MSet.norm")
# MSet.normbg <- bgcorrect.illumina(RGset)
# save(MSet.normbg, file="D:/mex_data/new_data/MSet.normbg")
load(file="D:/mex_data/new_data/MSet.normbg")

# Mset.swan <- preprocessSWAN(RGset, MSet.norm) #ash's method
#save(Mset.swan, file="D:/mex_data/new_data/Mset.swan")
load(file="D:/mex_data/new_data/Mset.swan")

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

# save the result
write.table(gb,"test.txt",row.names=T,quote=F,sep="\t")
write.table(gb,"test2.txt",row.names=T,col.names=NA,quote=F,sep="\t")




######
######
######
#TEST#
######
######
######

processData = function(dirPathToCSV, nameOfBetaFile){
  library(minfi)

  #mex dir
  setwd(dirPathToCSV)
  baseDir = dirPathToCSV

  targets <- read.450k.sheet(baseDir, recursive=FALSE)
  targets
  targets$Basename <- paste(baseDir,"/",targets$Slide, "/", targets$Slide, "_", targets$Array, sep="")
  RGset <- read.450k.exp(base = baseDir, targets = targets)
  pd <- pData(RGset) #phenodat
  MSet.raw <- preprocessRaw(RGset)
  MSet.norm <- preprocessIllumina(RGset, bg.correct = TRUE, normalize = "controls", reference = 2)
  Mset.swan <- preprocessSWAN(RGset, MSet.norm) #ash's method
  
  gb <- getBeta(Mset.swan)
  colnames(gb) <- targets$Sample_Name
  # save the result
  #write.table(gb,"test.txt",row.names=T,quote=F,sep="\t")
  write.table(gb,nameOfBetaFile,row.names=T,col.names=NA,quote=F,sep="\t")
}

#path to directory where .CSV file is located
path <- "D:/mex_data/new_data"
#name of new file containing processed beta values
filename <- "testLungBetas.txt"
processData(path, filename)
getwd()
