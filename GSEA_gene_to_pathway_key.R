# Read in annotation file
genetogenesets <- read.table("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/tcga/genetogroup4.txt", fill=TRUE)

# Create key
num_of_genes <- length(which(genetogenesets[,2]=="=>"))
num_of_entries <- nrow(genetogenesets) - num_of_genes
num_of_entries

matdat <- matrix(,nrow=num_of_entries, ncol=2)
gene <- vector()
counter=1
for(i in 1:nrow(genetogenesets)){
     if(genetogenesets[i,2]=="=>"){
          gene <- as.character(genetogenesets[i,1])
     }
     if(genetogenesets[i,2]==","){
          matdat[counter,1] <- gene
          matdat[counter,2] <- as.character(genetogenesets[i,1])
          counter=counter + 1
     }
}
dim(matdat)
dim(matdat[complete.cases(matdat),])
matdat <- matdat[complete.cases(matdat),]
head(matdat)
tail(matdat)

# matdat IS THE KEY
# Create key
pathway_key <- matdat
matdat
head(pathway_key)
tail(pathway_key)
# save(pathway_key, file="/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/tcga/pathway_key.rda")

load(file="/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/tcga/pathway_key.rda")
which(genetogenesets == ",")
num_of_entries

# topvar_genes = genes of interest
topvar_genes <- gsub("\\|\\d+", "", rownames(topvar))
length(topvar_genes)

pathkey_index <- which(pathway_key %in% topvar_genes)
result.mat <- pathway_key[pathkey_index,]
head(result.mat)
result.genesets <- result.mat[,2]

# print out top 25 pathways
head(rev(sort(table(result.genesets), )), n=30)
