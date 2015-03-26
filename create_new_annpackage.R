
library(AnnotationForge)
#library(AnnotationDbi)
#source("http://www.bioconductor.org/biocLite.R")
#biocLite("human.db0")


agilent_ann <- "/home/steve/Downloads/agilent_ann.txt"
t <- read.table(agilent_ann)
t
tmpout = tempdir()
makeDBPackage("HUMAN_DB",
              affy=FALSE,
              prefix="ag4x44",
              #fileName=agilent_ann,
              fileName="/home/steve/Downloads/agilent_ann.txt",
              baseMapType="gb",
              outputDir = tmpout,
              version="1.0.0",
              manufacturer = "Agilent",
              chipName = "Human Human GE 4x44K v2",
              manufacturerUrl = "http://www.agilent.com")

# OUTPUT: Creating package in /tmp/RtmpB5Scvk/ag4x44.db 

install.packages("/tmp/RtmpB5Scvk/ag4x44.db/", repos = NULL, type="source")
#install.packages("/tmp/Rtmp7erYfO/ag4x44.sqlite", repos=NULL, type="source")
library(ag4x44.db)
