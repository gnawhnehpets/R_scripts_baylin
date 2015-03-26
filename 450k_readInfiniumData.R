#
### BETA VALUES CALCULATION
#
# this R script extracts the beta values from user defined Infinium microarray data

# set the directory where the microarray data is stored:
setwd("D:/mex_data/new_data/")

# select the infinium data files, there are two options:
# enter the filenames manually
filenames = c("FinalReport_Sample Methylation Profile_RYen_11052012.txt")

# or: use an annotation file that contains the file names
#annotation = read.table("infiniumMicroarrayAnnotation.txt",header=TRUE,sep="\t")
#filenames = as.character(unique(annotation$microarrayFiles))
# WARNING: reading a lot of methylation files at the same time can take very long
# and will use quite some RAM memory
# it might be a good idea to make a selection of the files:
#filenames = filenames[8:10]

# read the methylation data
# create an empty list to store all the data
raw.meth.data = list()
for (t in c(1:length(filenames))){                                       
	raw.meth.data[[t]] = read.table(filenames[t], fill=T, skip=8, header=T, sep="\t", blank.lines.skip=TRUE, quote="\"", as.is=TRUE)
}                     

# combine the data in one table
# this combination is sometimes necessary as the time points for a single cell line
# can be found on seperate arrays
infinium = raw.meth.data[[1]]
if (length(filenames) > 1){
	for (t in 2:length(raw.meth.data)){
		infinium = cbind(infinium,raw.meth.data[[t]])
	}
}

# select the probe id, beta value and p value columns
infM = infinium[,c("TargetID",grep("AVG_Beta",colnames(infinium),value=TRUE),grep("Detection.Pval",colnames(infinium),value=TRUE))]

# select the probes which have a p value < 0.05
pValCol = grep("Detection.Pval",colnames(infM),value=TRUE)
for (t in 1:length(pValCol)){
	infM = infM[infM[,grep(pValCol[t],colnames(infM),value=TRUE)] < 0.05,]
}
infM = infM[,c("TargetID",grep("AVG_Beta",colnames(infM),value=TRUE))]

# clean up the column names
colnames(infM) = gsub("_\\d+.+\\.AVG_Beta$","",colnames(infM))
colnames(infM) = gsub("\\.","_",colnames(infM))
colnames(infM) = gsub("_+","_",colnames(infM))

# save the result
write.table(infM,"kate_extra_methylation_1.txt",row.names=F,quote=F,sep="\t")

