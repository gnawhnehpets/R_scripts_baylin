# NOT WORKING

library(lumi)
library(methylumi)
library(IlluminaHumanMethylation450k.db)
setwd("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/zahnow/")
filenames <- dir(pattern=".txt")
filenames
file <- filenames[3]
file
dat <- read.delim(file, sep = "\t". nrows)


profile <- lumiMethyR(filename="FinalReport_ryen_samples methylation profile_05082012.txt", lib="IlluminaHumanMethylation450k.db")




setwd("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/geo/methylation/")
samps <- read.table(system.file("ovarian/FinalReport_ryen_samples table_05082012.txt", package = "methylumi"), sep = "\t", header = TRUE)
mldat <- methylumiR(system.file("extdata/exampledata.samples.txt", package = "methylumi"), qcfile = system.file("extdata/exampledata.controls.txt", package = "methylumi"), sampleDescriptions = samps)