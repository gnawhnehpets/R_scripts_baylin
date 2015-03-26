rm(list=ls())
for(i in 1:20){gc()}
library(FDb.InfiniumMethylation.hg19)
m450 <- getPlatform(platform='HM450', genome='hg19')
show(m450)
class(m450)

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)

toggleProbes(IlluminaHumanMethylation450kanno.ilmn12.hg19, "all")

data(hg19.islands)
split(hg19.islands, seqnames(hg19.islands))

feat <- features(FDb.InfiniumMethylation.hg19)
feat
probes <- names(feat)
length(probes)

dat <- miscData(m450)








#biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(GenomicFeatures)
# https://www.ncbi.nlm.nih.gov/gene/?term=
trak2 <- "66008"
# https://www.ncbi.nlm.nih.gov/gene/?term=mlh1
mlh1 <- "4292"
txbygene <- transcriptsBy(txdb, by="gene")[trak2]
txbygene <- transcriptsBy(txdb, by="gene")[mlh1]
names(txbygene)
mget(names(txsByGene), org.Hs.egSYMBOL, ifnotfound=NA)
txbygene
# The transcript names corresponding to the trak2 gene will be used to subset the extracted intron and exon regions. The
# txbygene object is a GRangesList and the transcript names are a metadata column on the individual GRanges. To
# extract the names we must first 'flatten' or unlist txbygene.

txbygene
seqnames(txbygene)
chromosome <- seqnames(txbygene)
chromosome <- as.data.frame(chromosome)
chromosome

ranges(txbygene)
ranges <- ranges(txbygene)
ranges <- as.data.frame(ranges)
ranges

strand(txbygene)
strand <- strand(txbygene)
strand <- as.data.frame(strand)
strand

# names in this vector will give you ucsc entries
ucsc_id <- mcols(unlist(txbygene))$tx_name
ucsc_id
as.data.frame(ucsc_id)

cbind(ucsc_id, chromosome, ranges, strand)

source("http://bioconductor.org/biocLite.R")
biocLite("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other
names(anno)
probes <- rownames(anno)
probes
knowngenenames <- anno$UCSC_RefGene_Name
head(anno)
names(anno)
# [1] "Forward_Sequence"         "SourceSeq"                "Random_Loci"             
# [4] "Methyl27_Loci"            "UCSC_RefGene_Name"        "UCSC_RefGene_Accession"  
# [7] "UCSC_RefGene_Group"       "Phantom"                  "DMR"                     
# [10] "Enhancer"                 "HMM_Island"               "Regulatory_Feature_Name" 
# [13] "Regulatory_Feature_Group" "DHS"   

genename <- "MLH1"
indices <- grep(genename, anno$UCSC_RefGene_Name)
probename <- rownames(anno)[indices]
ucsc_group <- anno$UCSC_RefGene_Group[indices]
ucsc_gene <- anno$UCSC_RefGene_Name[indices]
genomic_location <- anno$Regulatory_Feature_Name[indices]
regulatory_feature_group <- anno$Regulatory_Feature_Group[indices]
sequence <- anno$SourceSeq[indices]
sequence
enhancer <- anno$Enhancer

dat <- as.data.frame(cbind(probename, genomic_location, sequence, ucsc_gene, ucsc_group, regulatory_feature_group))
dat

anno$HMM_Island

#EF-S 18-55mm IS STM Lens, EF-S 55-250mm IS Telephoto Lens, and Canon Camera Case




source("http://bioconductor.org/biocLite.R")
biocLite("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other
names(anno)
