### Creating own annotation package
source("http://www.bioconductor.org/biocLite.R")
biocLite("RH2")
biocLite("AnnotationForge")
biocLite("AnnotationDbi")
biocLite("RSQLite")
remove.packages("RSQLite")
install.packages("RSQLite")
setwd("/home/steve/Desktop/")

library(AnnotationDbi)
library(RSQLite)

source("http://www.bioconductor.org/biocLite.R")
biocLite("AnnotationForge")
library(AnnotationForge)

# agilent_ann <- "C:/Users/shwang26/Desktop/agilent_ann.txt"
agilent_ann <- "/home/steve/Downloads/agilent_ann.txt"
tmpout = tempdir()
makeDBPackage("HUMANCHIP_DB",
              affy=FALSE,
              prefix="ag4x44",
              fileName=agilent_ann,
              baseMapType="gb",
              outputDir = tmpout,
              version="1.0.0",
              manufacturer = "Agilent",
              chipName = "Human Human GE 4x44K v2",
              manufacturerUrl = "http://www.agilent.com")
# OUTPUT: Creating package in /tmp/RtmpB5Scvk/ag4x44.db 

install.packages("/tmp/RtmpB5Scvk/ag4x44.db/", repos = NULL, type="source")
library(ag4x44.db)


hcg110_IDs <- system.file("extdata",
            "hcg110_ID",
            package="AnnotationDbi")

tmpout = tempdir()
makeDBPackage("HUMANCHIP_DB",
              affy=FALSE,
              prefix="hcg110",
              fileName=hcg110_IDs,
              baseMapType="gb",
              outputDir = tmpout,
              version="1.0.0",
              manufacturer = "Affymetrix",
              chipName = "Human Cancer G110 Array",
              manufacturerUrl = "http://www.affymetrix.com")




biocLite("hcg110.db")



sessionInfo()




> sessionInfo()
R version 3.1.1 (2014-07-10)
Platform: x86_64-unknown-linux-gnu (64-bit)

locale:
     [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
[3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
[5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
[7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
[9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
     [1] stats4    parallel  stats     graphics  grDevices utils     datasets
[8] methods   base     

other attached packages:
     [1] AnnotationForge_1.8.2 org.Hs.eg.db_3.0.0    RSQLite_1.0.0        
[4] DBI_0.3.1             AnnotationDbi_1.28.1  GenomeInfoDb_1.2.4   
[7] IRanges_2.0.1         S4Vectors_0.4.0       Biobase_2.26.0       
[10] BiocGenerics_0.12.1  