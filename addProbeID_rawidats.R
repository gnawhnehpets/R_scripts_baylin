library(minfi)
setwd("D:/mex_data/new_data")
##
## This script adds the MySQL database probe id (from the probeIds.txt file) to each Infinium probe in an Infinium microarray file.
## It also creates a SQL file to load the methylation data in the MySQL database.
##

addIds = function(probeIdFile,infDataFile,tableName){
	
	#pids = read.table(probeIdFile,header=FALSE,sep="\t",as.is=TRUE)
	#save(pids, file="D:/mex_data/new_data/pids")
	load(file="D:/mex_data/new_data/pids")
	colnames(pids) = c("TargetID","id")

	#infData = read.table(infDataFile,header=TRUE,sep="\t",as.is=TRUE)
#	infDataRaw = read.table(infDataFile,header=TRUE,row.names=1, sep="\t",as.is=TRUE)
  #save(infDataRaw, file="D:/mex_data/new_data/infDataRaw")
	load(file="D:/mex_data/new_data/infDataRaw")
  infData <- infDataRaw
  head(pids)
  dim(pids)
  dim(infData)
  length(rownames(infData))
  pidsIndex <- which(pids$TargetID %in% rownames(infData))
	length(which(pids$TargetID %in% rownames(infData)))
#	length(which(rownames(infData %in% pids$TargetID)))
  
	pids_fil <- pids[pidsIndex,]
  dim(pids_fil)
  infIndex <- which(rownames(infData) %in% pids$TargetID)
  infData_fil <- infData[infIndex,]
#  infData2 = merge(pids_fil,infData_fil)  
  class(pids_fil)
  class(infData_fil)
	#infData = infData[,c(2,1,3:ncol(infData))]
  head(pids_fil)
  infData_merged<- cbind(pids_fil[,c(2,1)],infData_fil)
  head(infData_merged)

  #name of filtered data with probe ID keys
  resultsFile = gsub(".txt","_withID.txt",infDataFile)
  #resultsFile
	write.table(infData_merged,resultsFile,row.names=FALSE,sep="\t",quote=FALSE)
  colnames(infData_merged)
  tableName = "tableName1"

	loadData = paste("DROP TABLE IF EXISTS ",tableName,";\nCREATE TABLE ",tableName,"(\nid INT(11),\nprobeID VARCHAR(15),\n",sep="")

  infData <- infData_merged
	for (t in 3:(ncol(infData)-1)){
		loadData = paste(loadData,paste(colnames(infData)[t]," DECIMAL(11,10),\n",sep=""))
	}

	loadData = paste(loadData,paste(colnames(infData)[ncol(infData)])," DECIMAL(11,10)\n",sep="")
	
	loadData = paste(loadData,");\nLOAD DATA LOCAL INFILE '",resultsFile,"' INTO TABLE ",tableName," IGNORE 1 LINES;\nALTER TABLE ",tableName," ADD PRIMARY KEY (id);",sep="")
	loadData
	# create a function to reverse strings
	# strReverse <- function(x) sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
	# use the funtion to reverse strings to remove the file name from a path
	# (I couldn't find a way to do it without reversing the string)
# 	wd = strReverse(resultsFile)
# 	wd = gsub("^txt.+?/","/",wd,perl=TRUE)
# 	wd = strReverse(wd)
  wd <- paste(getwd(),"/", sep="")

	cat(loadData,file=paste(wd,tableName,".sql",sep=""),sep="")
}

# select the file with the Infinium 450k probe ids:
probeIdFile = "D:/mex_data/Rscripts/probeIds.txt"
# select the file with the Infinium methylation data:
infDataFile = "D:/mex_data/new_data/test2.txt"

addIds(probeIdFile,infDataFile,"tableName1")






