http://www.bioconductor.org/help/course-materials/2013/BioC2013/minfiLab.pdf
http://www.bioconductor.org/packages/release/bioc/vignettes/minfi/inst/doc/minfi.pdf

rm(list=ls())
for(i in 1:100){gc()}
library(minfi)
baseDir = "/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files"
setwd("/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files")

targets <- read.450k.sheet(baseDir, recursive=FALSE)
targets$Basename <- paste(baseDir,"/",targets$Slide, "/", targets$Slide, "_", targets$Array, sep="")
targets$Basename <- paste(baseDir,"/",targets$Slide, "_", targets$Array, sep="")

RGset <- read.450k.exp(base = baseDir, targets = targets)
baseDir

# save(RGset, file="Y:/users/shwang26/methylation_ex/idats/R objects/RGset")
# save(RGset, file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/RGset")
# load("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/methylation_ex/idats/R objects/RGset")
# load("Y:/users/shwang26/methylation_ex/idats/R objects/RGset")
load("/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/RGset")
# RGset <- load(file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/RGsetpairs")
pd <- pData(RGset) #phenodat
pd$Cell <- gsub("_.*", "", pd$Sample_Name)
colnames(pd)[5] <- "Phenotype"
colnames(pd)[4] <- "Cell"
pd$Phenotype <- gsub("_\\d+", "", pd$Sample_Name)
pd$Sample_Plate <- "cell"

Mset <- mapToGenome(RGset, mergeManifest = FALSE)
qc=minfiQC(Mset,fixOutliers=FALSE,verbose=TRUE)$qc
# save(qc,file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/qc.rda")
load(file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/qc.rda")

ind=which((qc[,1]+qc[,2])/2>= 10.5)
Mset1 <- preprocessIllumina(RGset[,ind])
Mset1 <- mapToGenome(Mset1)
pd=pData(Mset1)
# save(pd,file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/pd.rda")
# save(Mset1,file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/mset1.rda")
load(file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/pd.rda")
load(file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/mset1.rda")
rm(Mset1);gc();gc()

tmp92 <- targets[,1:4]

# tmp92[,2] <- gsub("_\\d$", "", tmp92[,1])
# tmp92[,3] <- gsub("_.*", "", tmp92[,1])
 tmp92[,3] <- gsub("_\\d$", "", tmp92[,1])
 tmp92[,2] <- gsub("_.*", "", tmp92[,1])

# tmp92[,4] <- gsub("_\\d$", "", tmp92[,1])
tmp92[,4] <- tmp92[,2]
compTable <- tmp92
colnames(compTable)<-c("name","group","type2","type1")
rownames(compTable)=1:nrow(compTable)
# save(compTable,file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/compTable.rda")
load(file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/compTable.rda")

load(file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/mset1.rda")
cmset <- cpgCollapse(Mset1,type="Illumina")
# save(cmset,file="rdas/cmset1.rda")
save(cmset,file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/cmset1.rda")
#############################################################################################
if(TRUE){
  load("Robjects/compTable.rda")
  load("Robjects/mset1.rda")
  load("Robjects/pd.rda")
    
  ##navg is average of first level (usually normal but sometimes
  ##adenoma) and cavg is the average of the second level
  
  avg <- cavg <- inds <- ttests <- vector("list",nrow(compTable))
  names(navg)<-names(cavg)<-names(inds)<-names(ttests)<-rownames(compTable)
  for(i in 1:nrow(compTable)){
    print(i)
    #ind=which(pd$Tissue==compTable$tissue[i] & !is.na(pd$Phenotype) & (pd$Phenotype==compTable$type1[i] | pd$Phenotype==compTable$type2[i]))
    ind=which(pd$Cell==compTable$group[i] & !is.na(pd$Phenotype) & (pd$Phenotype==compTable$type1[i] | pd$Phenotype==compTable$type2[i]))
    ind=which(pd$Cell==compTable$group[i] & !is.na(pd$Phenotype) & (pd$Phenotype==compTable$type1[i] | pd$Phenotype==compTable$type2[i]))

    inds[[i]]<-ind
#     status=factor(pd$Phenotype[ind],levels=rev(compTable[i,3:4]))
    status=factor(pd$Phenotype[ind],levels=rev(compTable[i,3:4]))
    tmpbeta=getBeta(Mset1[,ind],type="Illumina")
#     navg[[i]]=rowMeans(tmpbeta[,status==levels(status)[1],drop=FALSE])
#     cavg[[i]]=rowMeans(tmpbeta[,status==levels(status)[2],drop=FALSE])
    navg[[i]]=rowMeans(tmpbeta[,status==levels(status)[1],drop=FALSE])
    cavg[[i]]=rowMeans(tmpbeta[,status==levels(status)[2],drop=FALSE])

    ttests[[i]]=genefilter::rowttests(-tmpbeta,status)
    naind=which(is.na(ttests[[i]]$p.value))

    if(length(naind)>0){
      cat("Fixing NAs for",compTable[i,1],"\n")
      ttests[[i]]$p.value[naind]<-1 ###for one of these we have all 0
    }
    ##the negative is for sign of effecs match between blocks and
    ##tests. 
    
  }
  save(status,cavg,navg,ttests,inds,file="rdas/probe-results.rda")
  rm(list=ls())
}
## NOT COMPLETE #############################################################################
### computing average values at probe level
## navg is average of first level (usually normal but sometimes adenoma)
## and cavg is the average of the second level

for(i in 1:nrow(compTable)){
#   print(i)
#   ind=which(pd$Tissue==compTable$tissue[i] & !is.na(pd$Phenotype) &
#               (pd$Phenotype==compTable$type1[i] | pd$Phenotype==compTable$type2[i]))
#   inds[[i]]<-ind
  
#   status=factor(pd$Phenotype[ind],levels=rev(compTable[i,3:4]))
  status=factor(pd$Cell[i],levels=rev(compTable[i,3:4]))
  # }
  
#   tmpbeta=getBeta(Mset1[,ind],type="Illumina")
  tmpbeta=getBeta(Mset1[,i],type="Illumina")
  navg[[i]]=rowMeans(tmpbeta[,status==levels(status)[2],drop=FALSE])
  cavg[[i]]=rowMeans(tmpbeta[,status==levels(status)[2],drop=FALSE])
  ttests[[i]]=genefilter::rowttests(-tmpbeta,status)
  naind=which(is.na(ttests[[i]]$p.value))
  if(length(naind)>0){
    cat("Fixing NAs for",compTable[i,1],"\n")
    ttests[[i]]$p.value[naind]<-1 ###for one of these we have all 0
  }
  ##the negative is for sign of effecs match between blocks and
  ##tests. 
}
############################################################################################

###computing just average values for collapsed data

  load("rdas/compTable.rda")
  load("rdas/cmset1.rda")
  cavg <- navg <- vector("list",nrow(compTable))
  names(cavg)<-names(navg)<-rownames(compTable)
  for(i in 1:nrow(compTable)){
    print(i)
    #ind=which(pd$Tissue==compTable$tissue[i] & !is.na(pd$Phenotype) & (pd$Phenotype==compTable$type1[i] | pd$Phenotype==compTable$type2[i]))
    # find indices where PD_TISSUE == compTABLE_TISSUE, AND
    #                   PD_PHENOTYPE is not empty, AND
    #                   PD_PHENOTYPE == compTABLE_COL4, OR PD_PHENOTYPE == compTABLE_COL3
    #ind=which(pd$Cell==compTable$cell[1])
    ind <- which(pd$Cell==compTable$cell[i] & !is.na(pd$Phenotype) & (pd$Phenotype==compTable$type2[i] | pd$Phenotype==compTable$type1[i]))
    status=factor(pd$Phenotype[ind],levels=rev(compTable[i,3:4]))
    
    tmpbeta=getBeta(cmset$object[,ind])
    navg[[i]]=rowMeans(tmpbeta[,status==levels(status)[1],drop=FALSE])
    cavg[[i]]=rowMeans(tmpbeta[,status==levels(status)[2],drop=FALSE])
  }
  save(cavg,navg,file="rdas/results2.rda")
  rm(list=ls())
}
######## Figure 3 ####################################################################
load("rdas/cmset1.rda")
load("rdas/compTable.rda")
load("rdas/results.rda")
pd=pData(cmset$object)
gr=granges(cmset$obj)
regionType<-factor(gr$type,c("Island","Shore","Shelf","OpenSea"))
Indexes <- split(seq(along=regionType),regionType)
tissues=c("breast","colon", "liver", "lung", "pancreas",  "thyroid")
#######################################################
###prepare for figues
### this is more than we need, but during testing we looked at other figs
#######################################################
n<-length(Indexes[[1]])
tissues=unique(compTable[,2])
# canceravg<-array(NA,dim=c(n,5,length(tissues)))
# normalavg<-array(NA,dim=c(n,5,length(tissues)))
bjavg <-array(NA,dim=c(n,5,length(tissues)))
nonbjavg <-array(NA,dim=c(n,5,length(tissues)))
bj_htertavg <-array(NA,dim=c(n,5,length(tissues)))
bj_sv40avg <-array(NA,dim=c(n,5,length(tissues)))
bj_hrasavg <-array(NA,dim=c(n,5,length(tissues)))
epavg <-array(NA,dim=c(n,5,length(tissues)))
senavg <-array(NA,dim=c(n,5,length(tissues)))

dists<-array(NA,dim=c(n,5))
dist2block <- array(NA,dim=c(n,length(tissues),6))
Ns<-array(NA,dim=c(n,length(tissues)))
N=sapply(cmset$blockInfo$indexes,length)

for(h in seq(along=tissues)){
#   thei=which(compTable[,2]==tissues[h] & compTable[,3]=="cancer" & compTable[,4]=="normal")
  thei=which(compTable[,2]==tissues[h])
#   tIndex<-which(pd$Tissue==tissues[h])
  tIndex<-which(pd$Cell==tissues[h])
  betas=getBeta(cmset$obj[,tIndex])
  
  ##THE 5 dims are shelf,shore,island,shore,shelf
# #   cancertmp=rowMeans(betas[,which(pd$Phenotype[tIndex]=="cancer")])
# #   normaltmp=rowMeans(betas[,which(pd$Phenotype[tIndex]=="normal")])
#   bjtmp <- rowMeans(betas[,which(pd$Cell[tIndex]=="BJ")])
#   bj_hterttmp <- rowMeans(betas[,which(pd$Cell[tIndex]=="BJ_hTERT")])
#   bj_sv40tmp <- rowMeans(betas[,which(pd$Cell[tIndex]=="BJ_SV40")])
#   bj_hrastmp <- rowMeans(betas[,which(pd$Cell[tIndex]=="BJ_Hras")])
#   eptmp <- rowMeans(betas[,which(pd$Cell[tIndex]=="EP")])
#   sentmp <- rowMeans(betas[,which(pd$Cell[tIndex]=="Sen")])
#   
# #   canceravg[,3,h]<-cancertmp[Indexes[[1]]]
# #   normalavg[,3,h]<-normaltmp[Indexes[[1]]]
#   bjavg[,3,h]<-bjtmp[Indexes[[1]]]
#   bj_hterttmp[,3,h]<-bj_hterttmp[Indexes[[1]]]
#   bj_sv40tmp[,3,h]<-bj_sv40tmp[Indexes[[1]]]
#   bj_hrastmp[,3,h]<-bj_hrastmp[Indexes[[1]]]
#   epavg[,3,h]<-eptmp[Indexes[[1]]]
#   senavg[,3,h]<-sentmp[Indexes[[1]]]
   bjtmp=rowMeans(betas[,which(pd$Phenotype[tIndex]=="BJ")])
   nonbjtmp=rowMeans(betas[,which(pd$Phenotype[tIndex]!="BJ")])
   bjavg[,3,h]<-bjtmp[Indexes[[1]]]
   nonbjavg[,3,h]<-nonbjtmp[Indexes[[1]]]


  for(i in 2:3){
    for(j in 1:2){
      
      k<-(4-i)+(j-1)*(i-1)*2  ##this gives us 2,4,1,5 in the order of for loop
      
      mapfunc=ifelse(j==1,"follow","precede")
      subj=gr[Indexes[[i]],] ##just shores
      map=do.call(mapfunc,list(x=gr[Indexes[[1]],],subject=subj,ignore.strand=TRUE))
      nona=which(!is.na(map))
      dd=rep(NA,length(map))
      if(j==1) dd[nona]=start(gr[Indexes[[1]][nona],])-end(subj[map[nona],]) else dd[nona]= -end(gr[Indexes[[1]][nona],])+start(subj[map[nona],])
      ind=Indexes[[i]][map[nona]]
#       canceravg[nona,k,h]<-cancertmp[ind]
#       normalavg[nona,k,h]<-normaltmp[ind]
      bjavg[nona,k,h]<-bjtmp[ind]
      nonbjavg[nona,k,h]<-nonbjtmp[ind]

#       bjavg[nona,k,h]<-bjtmp[ind]
#       bj_htertavg[nona,k,h]<-bj_hterttmp[ind]
#       bj_sv40avg[nona,k,h]<-bj_sv40tmp[ind]
#       bj_hrasavg[nona,k,h]<-bj_hrastmp[ind]
#       epavg[nona,k,h]<-eptmp[ind]
#       senavg[nona,k,h]<-sentmp[ind]
      dists[,k]<-dd
      Ns[nona,h]<-N[ind]
    }
    bl=bsseq::data.frame2GRanges(blocks[[thei]]$table,keepColumns=TRUE)
    qvs=qvalue::qvalue(bl$p.value)$q
    map<-distanceToNearest(gr[Indexes[[1]],],bl)
    dist2block[,h,1]<-values(map)$distance
    dist2block[,h,2]<-bl$fwer[subjectHits(map)]
    dist2block[,h,3]<-qvs[subjectHits(map)]
    dist2block[,h,4]<-bl$L[subjectHits(map)]
    dist2block[,h,5]<-bl$value[subjectHits(map)]
    dist2block[,h,6]<-subjectHits(map)
  }
}


######################################################################################


# manifest <- getManifest(RGset)
# manifest
# probeInfo <- getProbeInfo(manifest)
# head(probeInfo)
# typeIProbes <- getProbeInfo(manifest, type = "I")$Name
# typeIIProbes <- getProbeInfo(manifest, type = "II")$Name

# QC
qcReport(RGset, sampNames = pd$Sample_Name, sampGroups = pd$Cell, pdf = "qcReport_prenorm.pdf")

# Density plots of methylation Beta plots
png(filename = "densityPlot.jpeg", height = 900, width = 1600)
densityPlot(RGset, sampGroups = pd$Cell,main = "Beta", xlab = "Beta")
dev.off()

# Density bean plots of methylation Beta values
png(filename = "densityBeanPlot.jpeg", height = 500, width = 500)
par(oma=c(2,10,1,1))
densityBeanPlot(RGset, sampGroups = pd$Cell, sampNames = pd$Sample_Name)
dev.off()

# Preprocessing (normalization)
# takes input a RGChannelSet --> outputs MethylSet
# number of preprocessing options, each one implemented as preprocessXXX where XXX = name of method
# "Raw" method = convert Red and Green channel into Methylated and Unmethylated signal

#MSet.raw <- preprocessRaw(RGset)
#save(MSet.raw, file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/MSet.raw")
load(file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/MSet.raw")
# load("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/methylation_ex/idats/R objects/MSet.raw")
# load("Y:/users/shwang26/methylation_ex/idats/R objects/MSet.raw")

# qc <- getQC(MSet.raw)
# png(filename = "plotQC.jpeg", height = 450, width = 800)
# plotQC(qc)
# dev.off()

# PRE-SWAN
# "Normal" method = background normalization, control normalization
# MSet.norm <- preprocessIllumina(RGset, bg.correct = TRUE, normalize = "controls", reference = 2)
# save(MSet.norm, file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/MSet.norm")
load(file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/MSet.norm")

png(filename = "densityPlotMSetnorm.jpeg", height = 900, width = 1600)
densityPlot(MSet.norm, sampGroups = pd$Cell,main = "Beta", xlab = "Beta")
dev.off()

png(filename = "controlStripPlot.jpeg", height = 900, width = 1600)
controlStripPlot(RGset, controls="BISULFITE CONVERSION II", sampNames = pd$Sample_Name)
dev.off()

png(filename = "mdsPlotMSetnorm.jpeg", height = 400, width = 600)
mdsPlot(MSet.norm, numPositions = 1000, sampGroups = pd$Cell, sampNames = pd$Sample_Name)
dev.off()


png(filename = "densityPlotSwan.jpeg", height = 900, width = 1600)
densityPlot(Mset.swan, sampGroups = pd$Cell,main = "Beta", xlab = "Beta")
dev.off()

png(filename = "mdsPlotMSetSwan.jpeg", height = 400, width = 600)
mdsPlot(Mset.swan, numPositions = 1000, sampGroups = pd$Cell, sampNames = pd$Sample_Name)
dev.off()


#####################
# After preprocessing
#####################
dat <- MSet.norm
sampleNames(dat)
sampleNames(dat) <- targets$Cell
sampleNames(dat)

# returns beta values betwen 0 and 1 where 1 is very high methylation; if type = "Illumina", then formula:
# beta = methylated / (methylated + unmeth + 100)
gb <- getBeta(dat)
colnames(gb)
sampleNames(dat)
bj <- getBeta(dat[,c(1,7)])
bj_hert <- getBeta(dat[,c(2,8)])
bj_sv40 <- getBeta(dat[,c(3,9)])
bj_hras <- getBeta(dat[,c(4,10)])
ep <- getBeta(dat[,5])
sen <- getBeta(dat[,6])
head(bj_hert)
bj_m <- as.matrix(rowMeans(bj))
hert_m <- as.matrix(rowMeans(bj_hert))
sv40_m <- as.matrix(rowMeans(bj_sv40))
hras_m <- as.matrix(rowMeans(bj_hras))
head(bj)
(bj[2,1]+bj[2,2])/2
head(bj_m)
t.test(ep, bj_m)
t.test(ep, hert_m)
t.test(ep, sv40_m)
t.test(ep, hras_m)
t.test(ep, sen)

t.test(sen, bj_m)
t.test(sen, hert_m)
t.test(sen, sv40_m)
t.test(sen, hras_m)
t.test(bj[,1], bj[,2])
sleep
with(sleep, t.test(extra[group == 1], extra[group == 2]))
head(bj)
# gbmed <- as.matrix(rowMedians(gb))
# #gbmed <- as.matrix(rowMeans(gb))
# rownames(gbmed) <- rownames(gb)


# cpgCollapse(object, what = c("Beta", "M"), maxGap = 500,
#             blockMaxGap = 2.5 * 10^5, maxClusterWidth = 1500,
#             dataSummary = colMeans, na.rm = FALSE,
#             returnBlockInfo = TRUE, islandAnno = NULL, verbose = TRUE,
#             ...)
# Arguments
# objects             An object of class [Genomic]MethylSet.
# what                Should operation be performed on the M-scale or Beta-scale?
# maxGap              Maximum gap between CpGs in a cluster
# blockMaxGap         Maximum block gap
# maxClusterWidth     Maximum cluster width
# dataSummary         Function used to summarize methylation across CpGs in the cluster.
# na.rm               Should NAs be removed when summarizing? Passed on to the dataSummary function.
# returnBlockInfo     Should the block annotation table be returned in addition to the block table?
# islandAnno          Which Island annotation should be used. NULL indicates the default. This argument is only useful if the annotatio object contains more than one island annotation.
# verbose             Should the function be verbose?

ratioSet <- ratioConvert(dat2, what = "both", keepCN = TRUE)
dat <- mapToGenome(MSet.norm, mergeManifest = TRUE)
class(dat)

cc <- cpgCollapse(dat, what="M", returnBlockInfo=FALSE)
load(file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/cc")
cc

# Usage
# blockFinder(object, design, coef = 2, what = c("Beta", "M"),
#             cluster = NULL, cutoff = NULL,
#             pickCutoff = FALSE, pickCutoffQ = 0.99,
#             nullMethod = c("permutation","bootstrap"),
#             smooth = TRUE, smoothFunction = locfitByCluster,
#             B = ncol(permutations), permutations = NULL,
#             verbose = TRUE, bpSpan = 2.5*10^5,...)
# Arguments
# object          An object of class GenomicRatioSet.
# design          Design matrix with rows representing samples and columns representing covariates. Regression is applied to each row of mat.
# coef            An integer denoting the column of the design matrix containing the covariate of interest. The hunt for bumps will be only be done for the estimate of this coefficient.
# what            Should blockfinding be performed on M-values or Beta values?
# cluster         The clusters of locations that are to be analyzed together. In the case of microarrays, the clusters are many times supplied by the manufacturer. If not available the function clusterMaker can be used to cluster nearby locations.4 blockFinder
# cutoff          A numeric value. Values of the estimate of the genomic profile above the cutoff or below the negative of the cutoff will be used as candidate regions. It is possible to give two separate values (upper and lower bounds). If one value is given, the lower bound is minus the value.
# pickCutoff      Should a cutoff be picked automatically?
# pickCutoffQ     The quantile used for picking the cutoff using the permutation distribution.
# nullMethod      Method used to generate null candidate regions, must be one of ‘bootstrap’ or ‘permutation’ (defaults to ‘permutation’). However, if covariates in addition to the outcome of interest are included in the design matrix (ncol(design)>2), the ‘permutation’ approach is not recommended. See vignette and original paper for more information.
# smooth          A logical value. If TRUE the estimated profile will be smoothed with the smoother defined by smoothFunction
# smoothFunction  A function to be used for smoothing the estimate of the genomic profile. Two functions are provided by the package: loessByCluster and runmedByCluster.
# B               An integer denoting the number of resamples to use when computing null distributions. This defaults to 0. If permutations is supplied that defines the number of permutations/bootstraps and B is ignored. permutations is a matrix with columns providing indexes to be used to scramble the data and create a null distribution. If this matrix is not supplied and B>0 then these indexes created using the function sample. verbose Should the function be verbose?
# bpSpan          Smoothing span. Note that this defaults to a large value becuase we are searching for large scale changes.

targets[,"Cell"]
design <- model.matrix(~ 0 + factor(c(1,2,3,4,5,6,1,2,3,4)))
colnames(design) <- c("BJ", "BJ_hTERT", "SV40", "BJ_Hras", "EP", "Sen")
head(design)

bf <- blockFinder(cc, design=design, cutoff=.05)
bf <- blockFinder(cc, design=design, coef=6, what="M", cutoff=.05, smooth=FALSE)

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
save(index, file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/index")

pre_fil <- premed[index,]
post_fil <- postmed[index,]
# rownames(premed)[index]
# rownames(postmed)[index]
deltabeta <- post_fil - pre_fil
deltabetaraw <- postmed - premed
#save(deltabeta, file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/deltabeta")
# load(file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/deltabeta")
#save(deltabetaraw, file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/deltabetaraw")
# load(file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/deltabetaraw")
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
png(filename="plotBetasbyType_comparison1.jpeg", height=900, width=1600)
par(mfrow=c(1,2))
# Distribution of betas for each sample
plotBetasByType(MsetEx[,1], main = "Raw")
plotBetasByType(Mset.swan[,1], main = "SWAN")
dev.off()

png(filename="plotBetasbyType_swan1.jpeg", height=900, width=1600)
par(mfrow=c(3,4))
for(i in 1:12){
  plotBetasByType(Mset.swan[,i], main = colnames(Mset.swan)[i])
}
dev.off()

png(filename="plotBetasbyType_swan2.jpeg", height=900, width=1600)
par(mfrow=c(3,4))
for(i in 13:24){
  plotBetasByType(Mset.swan[,i], main = colnames(Mset.swan)[i])
}
dev.off()

png(filename="plotBetasbyType_swan3.jpeg", height=900, width=1600)
par(mfrow=c(3,4))
for(i in 25:36){
  plotBetasByType(Mset.swan[,i], main = colnames(Mset.swan)[i])
}
dev.off()

png(filename="plotBetasbyType_swan4.jpeg", height=900, width=1600)
par(mfrow=c(3,4))
for(i in 37:48){
  plotBetasByType(Mset.swan[,i], main = colnames(Mset.swan)[i])
}
dev.off()

png(filename="plotBetasbyType_swan5.jpeg", height=900, width=1600)
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
save(ratioSetSwan, file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/ratioSetSwan")
load(file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/ratioSetSwan")
# Associate genomic coordinates tot he probes using mapToGenome: transform ratioSet to a GenomicRatioSet (hold M and/or Beta values together with associated genomic coordinates)
# gRatioSet <- mapToGenome(ratioSet, mergeManifest = TRUE)
gRatioSetSwan <- mapToGenome(ratioSetSwan, mergeManifest = TRUE)
save(gRatioSetSwan, file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/gRatioSetSwan")
load(file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/gRatioSetSwan")
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
png(filename="coloBetaIsland_boxplot.jpeg", height=900, width=1600)
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
png(filename="coloBetaIslandPromoter_boxplot.jpeg", height=900, width=1600)
boxplot(colo.beta.island.promoter,las=2)
dev.off()

colo.beta.island.promoter.methylated<-colo.beta.island.promoter[which(apply(colo.beta.island.promoter,1,max)>=0.5),]
dim(colo.beta.island.promoter.methylated)
png(filename="coloBetaIslandPromoterMethylated_boxplot.jpeg", height=900, width=1600)
boxplot(colo.beta.island.promoter.methylated, las=2)
dev.off()


nshore <- as.character(ga$Name[which(ga$Relation_to_Island=='N_Shore')])
colo.beta.nshore <- gbs[nshore,]
head(colo.beta.nshore)
dim(colo.beta.nshore)
png(filename="coloBetaNShore_boxplot.jpeg", height=900, width=1600)
boxplot(colo.beta.nshore, las=2)
dev.off()

sshore <- as.character(ga$Name[which(ga$Relation_to_Island=='S_Shore')])
colo.beta.sshore <- gbs[sshore,]
head(colo.beta.sshore)
dim(colo.beta.sshore)
png(filename="coloBetaSShore_boxplot.jpeg", height=900, width=1600)
boxplot(colo.beta.sshore, las=2)
dev.off()

sshelf <- as.character(ga$Name[which(ga$Relation_to_Island=='S_Shelf')])
colo.beta.sshelf <- gbs[sshelf,]
head(colo.beta.sshelf)
dim(colo.beta.sshelf)
png(filename="coloBetaSShelf_boxplot.jpeg", height=900, width=1600)
boxplot(colo.beta.sshelf, las=2)
dev.off()

nshelf <- as.character(ga$Name[which(ga$Relation_to_Island=='N_Shelf')])
colo.beta.nshelf <- gbs[nshelf,]
head(colo.beta.nshelf)
dim(colo.beta.nshelf)
png(filename="coloBetaNShelf_boxplot.jpeg", height=900, width=1600)
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
save(gRanges, file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/gRanges")
load(file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/gRanges")
head(gRanges, n= 3)

# Finding differentially methylated positions (DMPs)
# Categorical phenotypes - 'dmpFinder' uses F-test to identify positions that are differentially methylated between two (or more) groups

# Find differences between GroupA and GroupB
table(pd$Cell)
mset <- Mset.swan
M <- getM(mset, type = "beta", betaThreshold = 0.001)
# M
# Returns a table of CpG positions sorted by differential methylation p-value
# Tests each genomic position for assocaition between methylation and a phenotype
dmp <- dmpFinder(M, pheno=pd$Cell, type="categorical")
save(dmp, file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/dmp")
load(file="/amber2/scratch/baylin/shwang/proj/Wenbing/BJ_system_Methylation/idat_files/Robjects/dmp")
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
png(filename="plotCpGindividualPositions.jpeg", height=900, width=1900)
cpgs <- rownames(dmp)[1:75]
length(rownames(dmp))
par(mfrow=c(5,15))
plotCpg(mset, cpg=cpgs, pheno=pd$Cell)
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
# png(filename="plotCpGindividualPositions-continuous2.jpeg", height=900, width=1600)
# cpgs <- rownames(dmp)[1:4]
# par(mfrow=c(2,2))
# plotCpg(mset, cpg=cpgs, type="continuous", pheno=continuousPheno, xlab="Phenotype 1")
# dev.off()




# beta <- getBeta(gRatioSet.quantile)
# names(pData(gRatioSet.quantile))
# prepost <- pData(gRatioSet.quantile)$Cell
# dmp <- dmpFinder(beta, pheno = prepost , type = "categorical")
# head(dmp)
# 
# pheno <- pData(gRatioSet.quantile)$Cell
# designMatrix <- model.matrix(~ pheno)
# dmrs <- bumphunter(gRatioSet.quantile, design = designMatrix, cutoff = 0.05, B=1)

# CONVERT Mset.swan (methylset) to a genomicRatioSet

test <- Mset.swan
ratioSet <- ratioConvert(Mset.swan, what = "both", keepCN = TRUE)
gRatioSet <- mapToGenome(ratioSet, mergeManifest = TRUE)

beta <- getBeta(gRatioSet)
names(pData(gRatioSet))
prepost <- pData(gRatioSet)$Cell
dmp <- dmpFinder(beta, pheno = prepost , type = "categorical")
head(dmp)

pheno <- pData(gRatioSet)$Cell
designMatrix <- model.matrix(~ pheno)
dmrs <- bumphunter(gRatioSet, design = designMatrix, cutoff = 0.05, B=1)
save(dmrs, file="./Robjects/dmrs")
load("./Robjects/dmrs")
dmrindex <- which(dmrs$table$p.value < .05)
dmrs2 <- dmrs$table
dmrindex2 <- which(dmrs2$p.value < .05)
dmrs3 <- dmrs2[dmrindex2,]
pheno



# colon biopsy: 