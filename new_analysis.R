library(limma)
setwd("/home/steve/Desktop/MA_core_download/ahuja/expression/08092013 analysis/")
#Just pre-post, post-pre data
#targets <- targets[c(2,6:9),]
targets <- readTargets("/home/steve/Desktop/MA_core_download/ahuja/expression/08092013 analysis/targets_full.txt")
targets_fil <- targets[c(1:3,5:7,9:14,16),]
RG <- read.maimages(targets, source="agilent", path="/home/steve/Desktop/MA_core_download/ahuja/expression/08092013 analysis/")
RG_fil <- read.maimages(targets_fil, source = "agilent", path ="/home/steve/Desktop/MA_core_download/ahuja/expression/08092013 analysis/")


RG$weights <-  matrix(rep(RG$genes$ControlType,ncol(RG$R)),ncol=ncol(RG$R),byrow=F)
RG$weights[RG$genes$ControlType!=0,] <- 0
RG$weights[RG$genes$ControlType==0,] <- 1

#normalization
#normwithin <-normalizeWithinArrays(RG,method='loess',bc.method='normexp', offset=50)
normwithin <-normalizeWithinArrays(RG,method='loess',bc.method='none')
normbetween <-normalizeBetweenArrays(normwithin,method='Aquantile')

#remove controls
nc_normbetween <-normbetween[normbetween$genes$ControlType==0,]

############################################################################################
############################################################################################

#Convert M values to RG intensities
ma2rg <- RG.MA(nc_normbetween)

#RGb = backgroundCorrect(RG, method="normexp", offset=50)
#RGb2 = backgroundCorrect(RG, method="none")
RGb = backgroundCorrect(ma2rg, method="none")

#All red intensities
cy3 = RGb$R
rownames(cy3) <- RGb$genes$ProbeName
colnames(cy3) <- paste(targets$Cy3, targets$Cy3_Sample)
#cy3lognorm =normalizeBetweenArrays(log2(cy3), method="quantile") #logscaled Cy3 data
cy3_log <- log2(cy3)

#All green intensities
cy5 = RGb$G
rownames(cy5) <- RGb$genes$ProbeName
colnames(cy5) <- paste(targets$Cy5, targets$Cy5_sample)
#cy5lognorm =normalizeBetweenArrays(log2(cy5), method="quantile") #logscaled Cy5 data
#cy5norm =normalizeBetweenArrays(cy5, method="quantile") #nonlogscaled Cy5 data
cy5_log <- log2(cy5)

########################################################
#ALTERNATIVE #1 - collapse probes into genes before subsetting into pre/post groups
#This seems to be the correct process
#> length(rownames(cy3))          #length of unique genenames in original MA
#[1] 43118
#> length(unique(rownames(cy3)))
#[1] 27193
#Using pre-post values instead of RGb
#Create dataframe
genes <- as.matrix(rownames(cy3))
#1st col=genename, 2nd col=cy3 values
three_ap <- cbind (genes, cy3)
five_ap <- cbind (genes, cy5)
colnames(three_ap)[1] <- "genename"
colnames(five_ap)[1] <- "genename"
three_ap <- as.data.frame(three_ap)
five_ap <- as.data.frame(five_ap)

#Collapse similar probes into their gene name
agg3 <- aggregate(. ~ genename, data=three_ap, median) #genes: 27193
agg5 <- aggregate(. ~ genename, data=five_ap, median) #genes: 27193

cy3_agg <- agg3[,c(2:17)]
cy5_agg <- agg5[,c(2:17)]
rownames(cy3_agg) <- agg3[,1]
rownames(cy5_agg) <- agg5[,1]
#Group pre's with pre's, post's with post's
agg_pre = cbind(cy3_agg[,c(2:4, 8:10, 12, 14:15)],cy5_agg[,c(1,4:8,11,13,15:16)])
agg_post = cbind(cy3_agg[,c(1,5:7,11,13,16)],cy5_agg[,c(2:3,9:10,12,14)])

### cy3_agg_log???

#############################################################
#####################################################
#SUCCESSFUL
#Create dataframe
#genes <- RGb$genes$GeneName
#genes <- as.matrix(RGb$genes$GeneName)
#1st col=genename, 2nd col=cy3 values
#three <- cbind (genes, RGb$G)
#five <- cbind (genes, RGb$R)
#colnames(three)[1] <- "genename"
#colnames(five)[1] <- "genename"
#three <- as.data.frame(three)
#five <- as.data.frame(five)

#Collapse similar probes into their gene name
#agg3 <- aggregate(. ~ genename, data=three, median)
#agg5 <- aggregate(. ~ genename, data=five, median)

######################################################
#ALTERNATIVE #2 - subset into pre/post group before collapsing probes into similar genes
#Group pre's with pre's, post's with post's
pre_pa = cbind(cy3[,c(2:4, 8:10, 12, 14:15)],cy5[,c(1,4:8,11,13,15:16)])
post_pa = cbind(cy3[,c(1,5:7,11,13,16)],cy5[,c(2:3,9:10,12,14)])

#pre_log = cbind(cy3_log[,c(2:4, 8:10, 12, 14:15)],cy5_log[,c(1,4:8,11,13,15:16)])
#post_log = cbind(cy3_log[,c(1,5:7,11,13,16)],cy5_log[,c(2:3,9:10,12,14)])

#Using pre-post values instead of RGb
#Create dataframe
genes_pa <- as.matrix(rownames(pre_pa))
#1st col=genename, 2nd col=cy3 values
three_pa <- cbind (genes_pa, pre_pa)
five_pa <- cbind (genes_pa, post_pa)
colnames(three_pa)[1] <- "genename"
colnames(five_pa)[1] <- "genename"
three_pa <- as.data.frame(three_pa)
five_pa <- as.data.frame(five_pa)

#Collapse similar probes into their gene name
pre_ag <- aggregate(. ~ genename, data=three_pa, median) #genes: 19158
post_ag <- aggregate(. ~ genename, data=five_pa, median) #genes: 15250
rownames(pre_ag) <- pre_ag[,1]
rownames(post_ag) <- post_ag[,1]
pre_agg <- pre_ag[,c(2:20)]
post_agg <- post_ag[,c(2:14)]

#####################################################


#foreach each channel, write out ranked list
#Need to adjust x in `paste("x_")`#Create files for `pre_agg` `post_agg` `agg_post` and `agg_pre`
writeout <- function(x){
    for(i in 1:length(colnames(x))){ 
    sample <- colnames(x)[i]
    filename <- paste("pre_agg_",i,"-",sample,".rnk")
    filename <- as.character(filename)
    filename <- gsub(" ","", filename)
    filename <- as.character(filename)
    write.table(x[,i], file=filename, row.names=row.names(x), sep="\t", col.names=FALSE, quote=FALSE)
  }
}
writeout(pre_agg)

1

######################################################
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

