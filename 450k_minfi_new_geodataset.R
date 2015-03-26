library(minfi)
# cluster dir
setwd("/home/jhmi/shwang/proj/mdanderson/idats")
baseDir = "/home/jhmi/shwang/proj/mdanderson/idats"

# READ IN
targets <- read.450k.sheet(baseDir, recursive=FALSE)
targets <- read.table(file="minfisample.txt", header=TRUE)

#targets$Basename <- paste(baseDir,"/",targets$Slide, "/", targets$Slide, "_", targets$Array, sep="")
RGset <- read.450k.exp(base = baseDir, targets = targets)
save(RGset, file="/home/jhmi/shwang/proj/mdanderson/idats/R objects/RGset")
load("/home/jhmi/shwang/proj/mdanderson/idats/R objects/RGset")

sampleNames(RGset)

gsub(".*/idats/\\d+/", "", targets$Basename)
which(sampleNames(RGset) %in% gsub(".*/idats/\\d+/", "", targets$Basename))
data.frame(targets$Sample_Name, gsub(".*/idats/\\d+/", "", targets$Basename))
meta <- data.frame(targets$Sample_Name, gsub(".*/idats/\\d+/", "", targets$Basename))
colnames(meta) <- c("SampleName", "BaseName")
#Convert filename to samplenames
sampleNames(RGset) <- targets$Sample_Name

# RGset <- load(file="/home/jhmi/shwang/proj/mdanderson/idats/R objects/RGsetpairs")
pd <- pData(RGset) #phenodat

# PRE-SWAN
MSet.raw <- preprocessRaw(RGset)
MSet.norm <- preprocessIllumina(RGset, bg.correct = TRUE, normalize = "controls", reference = 2)
MSet.bg <- bgcorrect.illumina(RGset)
save(MSet.raw, file="/home/jhmi/shwang/proj/mdanderson/idats/R objects/MSet.raw_geo")
save(MSet.norm, file="/home/jhmi/shwang/proj/mdanderson/idats/R objects/MSet.norm_geo")
save(MSet.bg, file="/home/jhmi/shwang/proj/mdanderson/idats/R objects/MSet.bg_geo")

# SWAN
Mset.swan <- preprocessSWAN(RGset, MSet.raw) #original method
Mset.swan2 <- preprocessSWAN(MSet.bg) #bg correction
Mset.swan3 <- preprocessSWAN(RGset, MSet.norm) #ash's method
save(Mset.swan, file="/home/jhmi/shwang/proj/mdanderson/idats/R objects/Mset.swan_geo")
save(Mset.swan2, file="/home/jhmi/shwang/proj/mdanderson/idats/R objects/Mset.swan2_geo")
save(Mset.swan3, file="/home/jhmi/shwang/proj/mdanderson/idats/R objects/Mset.swan3_geo")

Hello,

We met with the statistician who works with high dimensional data today.  Ive included Dr Yin on this email.   She is going to be taking over as primary statistician 
for the P2C very soon (replacing Dr Shi) and we are trying to get her up to speed on various projects in the program.

One of the issues on this study is that the primary endpoint was assessed on the 2nd cohort of patients (N=22) and after we changed the patient population due to excessive 
toxicity, enrolling a more diseased population than expected,  and thus adding 2 criteria for eligibility:

1)  No more than 30% of tumor burden in the liver
2)	Patients having more 3 or more prior therapies had to be of wild type KRAS status.

Having 2 different patient populations will impact the analysis.  This is problematic in terms of using all patients at one time.  First, the patient population was more 
diseased burdened than the final population.  Second, The RECIST criteria were modified when the study started, making patients show 30% increase in the sum of the lesions (vs 20%).
The 2nd group of patients were measured according to the 20% rule (ie, true RECST v1.1 requirements).  Given that we have the measurements, we could go back and identify the new 
event dates, but patients would have been treated a cycle or more longer than they should have.  Third, I do not think that we retrospectively determined volumetrics and KRAS 
status on the 24 patients initially enrolled.   We could potentially try to select out those patients in the first 24 who are technically the same as those enrolled in the 
latter group.   It would need some additional work, however.

With all of that in mind, here are the questions that we have:

1)	Normalization
	a)	How were the results normalized?   Was it by the Beta value or the M and U values?  Do you have a reference for the method used?

  The data was normalized by the SWAN method outlined in Bioconductor minfi package. http://www.bioconductor.org/packages/release/bioc/vignettes/minfi/inst/doc/minfi.pdf
This protocol reads in a set of .idat files containing the raw red/green channel intensities as a RGChannelSet object. This RG object is then processed into a MethylSet object
which contains normalized data and consists of two matrices containing the methylated and unmethylated evidence for each CpG. 

	b)	Were the samples normalized with all of the other tumor types you have mentioned before?  Given that youve sorted out only the colon patients, it leads us to believe 
that all of the tumor types (mentioned in Dec 2013) were ran at the same time and thus normalized together.  If so, does it matter that only a small number other probes are 
differentially methylated between the cancers?

  The colon samples were normalized only within the colon samples. 

c)	Should we normalize the 1st cohort separate from the 2nd cohort?  Is there any biologic reason to believe that the patients in the first group are different genetically 
from those in the 2nd group?

2)	What is the difference in the files between:
  a) 	Corrected and uncorrected?
b)	Methylated intensity and unmethylated intensity?   Are these the M & U needed to calculate the beta value?

Please disregard as these values were generated by Illuminas GenomeStudio software. 

3)	Is the filter pvalue the detection values?  If so, what threshold did you use?

In order to extract the values that you want, you can use the following minfi functions:
  For beta values, the getBeta() function: getBeta(Robject, type"", betaThreshold=0.001)
  For M values, the getM() function: getM(Robject, type"", betaThreshold=0.001)
  For Meth values, the getMeth() function: getMeth(Robject, type"", betaThreshold=0.001)
  For Unmeth values, the getUnmeth() function: getUnmeth(Robject, type"", betaThreshold=0.001)

The answers to these questions will help us on our end.   We really do need to work out the cohort that we will use and what types of analysis will be performed.  We are limited in terms of what we can say, as there are no responses in our study and the sample size is 24 and 22 for the 2 cohorts.  This also needs to be considered if this data is ever used again or merged with other samples and tumor types.

# only colon samples
colonsamples <- Mset.swan[,-c(61:76)]
save(colonsamples, file="/home/jhmi/shwang/proj/mdanderson/idats/colonsamples")
annotation <- getAnnotation(test)
annotation$UCSC_RefGene_Name
annotation$Name
#p16, SFRP 1, TFPI-2 GATA4, GATA5, HIN1, O6-MGMT, APC, CHFR, RASSF1A, and hMLH1
which(ga$UCSC_RefGene_Name == "ARPC5") #p16
which(ga$UCSC_RefGene_Name == "CDKN2A") #p16
which(ga$UCSC_RefGene_Name == "SFRP1")
which(ga$UCSC_RefGene_Name == "GATA4")
which(ga$UCSC_RefGene_Name == "GATA5")
which(ga$UCSC_RefGene_Name == "SCGB3A1") #HIN1
which(ga$UCSC_RefGene_Name == "MGMT") #O6
which(ga$UCSC_RefGene_Name == "APC")
which(ga$UCSC_RefGene_Name == "CHFR")
which(ga$UCSC_RefGene_Name == "RASSF1") #RASSF1A
which(ga$UCSC_RefGene_Name == "MLH1") #hMLH1

c("ARPC5", "CDKN2A", "SFRP1", "GATA4", "GATA5", "SCGB3A1", "MGMT", "APC", "CHFR", "RASSF1", "MLH1")

# find probe names associated with geneset
anndat <- cbind(ga$Name, ga$UCSC_RefGene_Name)
tmp <- setNames(strsplit(anndat[,2], split = ';'), anndat[,1])
tmp2 <- data.frame(probe = rep(names(tmp), sapply(tmp, length)), genename = unlist(tmp), row.names = NULL)
uniqtmp <- unique(tmp2)
uniqtmp
# index of dat associated with annotated genes
index <- which(uniqtmp[,2] %in% c("ARPC5", "CDKN2A", "SFRP1", "GATA4", "GATA5", "SCGB3A1", "MGMT", "APC", "CHFR", "RASSF1", "MLH1"))
filteredtmp <- uniqtmp[index,]
geneprobes <- filteredtmp[,1]
geneprobes

index2 <- which(ga$Name %in% geneprobes)
# filteredcolonsamples <- colonsamples[geneprobes,]
filteredcolonsamples <- colonsamples[index2,]
filteredcolonsamples <- Mset.swan[index2,]

save(filteredcolonsamples, file="/home/jhmi/shwang/proj/mdanderson/idats/filteredcolonsamples")


# POST-SWAN - Method 1 #################################################################
# filter betas
# gb <- getBeta(Mset.swan)
# ga <- getAnnotation(Mset.swan)
gb <- getBeta(Mset.swan)
ga <- getAnnotation(Mset.swan)

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
index <- which(premed > .3) #317367
# index <- which(premed > .4) #295616
# index <- which(premed > .5) #267405
length(index)

# save(index, file="/home/jhmi/shwang/proj/mdanderson/idats/R objects/index")
# load(file="/home/jhmi/shwang/proj/mdanderson/idats/R objects/index")
pre_fil <- premed[index,]
post_fil <- postmed[index,]

# filter for probes with deltabeta >.2
deltabeta <- post_fil - pre_fil
deltabetaraw <- postmed - premed
# sigprobes <- names(which(abs(deltabeta) > .2))
sigprobes <- names(which(deltabeta < -.2))
# sigprobes <- names(which(abs(deltabeta) > .2)) # .5 - 706

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
genes3 <- uniquesig #715
genes4 <- uniquesig #625
genes5 <- uniquesig #420
length(genes3)
which(genes3 %in% genes4 == FALSE)
genes.4 <- genes3[which(genes3 %in% genes4 == FALSE)]
print(genes.4, quote=FALSE, row.names=FALSE)

which(genes3 %in% genes5 == FALSE)
genes.5 <- genes3[which(genes3 %in% genes5 == FALSE)]
print(genes.5, quote=FALSE, row.names=FALSE)

which(genes4 %in% genes5 == FALSE)
genes.45 <- genes4[which(genes4 %in% genes5 == FALSE)]
print(genes.45, quote=FALSE, row.names=FALSE)

# sigprobes <- names(which(abs(deltabeta) > .2))
#.3 - 780
#.4 - 680
#.5 - 461
# sigprobes <- names(which(deltabeta < -.2))
#.3 - 724
#.4 - 626
#.5 - 423


aim <- as.matrix(read.table("aimgenes.txt", header=FALSE))
intersect(aim, uniquesig)


# POST-SWAN - Method 2 #################################################################
# Finding differentially methylated positions (DMPs)
# Categorical phenotypes - 'dmpFinder' uses F-test to identify positions that are differentially methylated between two (or more) groups

# Find differences between GroupA and GroupB
table(pd$Sample_Group)
mset <- Mset.swan
M <- getM(mset, type = "beta", betaThreshold = 0.001)
#M <- getM(mset, type = "beta", betaThreshold = 0.5)
# M
# Returns a table of CpG positions sorted by differential methylation p-value
# Tests each genomic position for assocaition between methylation and a phenotype
dmp <- dmpFinder(M, pheno=pd$Sample_Group, type="categorical")
#save(dmp, file="/home/jhmi/shwang/proj/mdanderson/idats/R objects/dmp")
#load(file="/home/jhmi/shwang/proj/mdanderson/idats/R objects/dmp")

# List of most differentially methylated probes
dmp_names <- rownames(dmp)[dmp$qval < .05]
dmp_names <- rownames(dmp)[dmp$qval < .01]
dmp_names <- rownames(dmp)[dmp$qval < .0001]
#dmp_names <- rownames(dmp)[dmp$pval < .000001]
length(dmp_names) #17395, 20773
# index of most dmp
index <- rownames(ga) %in% dmp_names
head(rownames(ga)[index]) #6
# genename of most dmp
dmp_genename <- ga$UCSC_RefGene_Name[index]
uniquegenes <- unique(dmp_genename)
uniquegenes <- unique(unlist(strsplit(uniquegenes, ";"))) #intersect with AIM gene list      #6159 unique genes, 7055
length(uniquegenes)
length(intersect(uniquegenes, genes3))
#print(sort(intersect(uniquegenes, genes3)), quote=FALSE)
print(sort(intersect(uniquegenes, genes5)), quote=FALSE)
data.frame(sort(intersect(uniquegenes, genes5))) #USE THIS OUTPUT WITH getHTMLsourceCode2.pl TO GET GENE ANNOTATION



# CROSS WITH AIM GENE LIST
aim <- as.matrix(read.table("aimgenes.txt", header=FALSE)) # oncotarget AIM gene list
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
