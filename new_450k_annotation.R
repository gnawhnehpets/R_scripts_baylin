# source("http://bioconductor.org/biocLite.R")
# biocLite("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other
names(anno)
head(anno)
probes <- rownames(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other)
length(probes)
anno$UCSC_RefGene_Name
anno$UCSC_RefGene_Accession
anno$UCSC_RefGene_Group
length(anno$UCSC_RefGene_Name)
head(anno$UCSC_RefGene_Accession, n=100)
length(anno$UCSC_RefGene_Group)

dat <- data.frame(probes, anno$UCSC_RefGene_Accession)
colnames(dat) <- c("probes", "refgene")
head(dat)

# tmp <- setNames(strsplit(as.character(x$Gene), split=';'), x$Probe)
tmp <- setNames(strsplit(as.character(dat$refgene), split=';'), dat$probes)
head(tmp)
table <- data.frame(probename = rep(names(tmp), sapply(tmp, length)), refgenename = unlist(tmp), row.names = NULL)
final <- unique(table)
data.frame(head(dat, n=50), head(final, n=50))


source("http://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")
library("org.Hs.eg.db")

x <- org.Hs.egREFSEQ2EG
# Get the RefSeq identifier that are mapped to an entrez gene ID
mapped_seqs <- mappedkeys(x)
head(cbind(x, mapped_seqs))

# Convert to a list
xx <- as.list(x[mapped_seqs])
xx
if(length(xx) > 0) {
     # Get the entrez gene for the first five Refseqs
     xx[1:15]
     # Get the first one
     xx[[1]]
}
mapped_key <- as.matrix(unlist(xx))
mapped_key

final
index <- which(final[,2] %in% rownames(mapped_key))
index
entrez_id <- mapped_key[index]
length(entrez_id)
head(entrez_id)
dim(final)


mat <- matrix(, nrow=nrow(final), ncol=3)

mat[,1] <- as.character(final[,1])
mat[,2] <- as.character(final[,2])
mat[index,3]<- entrez_id
length(rownames(mapped_key))
length(final[,2])
head(mat, n=1000)


which(1:6 %in% 0:36)
which(0:36 %in% 1:6)

setwd("/home/steve/Downloads/")
dat <- read.delim("MCF7_450k_islands.txt")
dat <- dat[,c(1,2,4,5)]
x <- dat
head(x)
tmp <- setNames(strsplit(as.character(x$Gene), split=';'), x$Probe)
head(tmp)
table <- data.frame(probe = rep(names(tmp), sapply(tmp, length)), genename = unlist(tmp), MCF10a = rep(x[,3], sapply(tmp, length)), MCF7 = rep(x[,4], sapply(tmp, length)), row.names = NULL)
final <- unique(table)
data.frame(head(x, n=100), head(table, n=100), head(final, n=100))
data.frame(head(x, n=50), head(final, n=50))
