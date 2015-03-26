# This script is for annotating Illumina methylation data

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
dat <- IlluminaHumanMethylation450kanno.ilmn12.hg19
class(dat)
dat

class(dat)



# ecad <- new.datfil["CDH1",]
# p16 <- new.datfil["CDKN2A",]
# tms1 <- new.datfil["PYCARD",]
# timp1 <- new.datfil["TIMP1",]
# timp2 <- new.datfil["TIMP2",]
# timp3 <- new.datfil["TIMP3",]
# timp4 <- new.datfil["TIMP4",]
# mlh1 <- new.datfil["MLH1",]
# mgmt <- new.datfil["MGMT",]
# pd1 <- new.datfil["PDCD1",]
# pdl1 <- new.datfil["CD274",]
# pdl2 <- new.datfil["PDCD1LG2",]
# sfrp1 <- new.datfil["SFRP1",]
# tfpi2 <- new.datfil["TFPI2",]
# gata4 <- new.datfil["GATA4",]
# gata5 <- new.datfil["GATA5",]
# apc <- new.datfil["APC",]
# chfr <- new.datfil["CHFR",]
# rassf1 <- new.datfil["RASSF1",]
# hin1 <- new.datfil["SCGB3A1",]

x <- IlluminaHumanMethylation450kCHRLOC
# Get the probe identifiers that are mapped to chromosome locations
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])
if(length(xx) > 0) {
  xx[[232]]
}


x <- IlluminaHumanMethylation450kCHR37
xx <- as.list(x)[1:10]
xx
xx

library(IlluminaHumanMethylation450k.db)
x <- IlluminaHumanMethylation450kSYMBOL
loc <- IlluminaHumanMethylation450kCHRLOC
loc <- IlluminaHumanMethylation450kCHR

x
mapped_probes <- mappedkeys(x)
mapped_loc <- mappedkeys(loc)
mapped_probes
xx <- as.list(x[mapped_probes])
yy <- as.list(loc[mapped_loc])
names(yy[[245345]])
names(yy[[2123]])
v <- yy[[2123]]
v
xx[[22653]]
gname <- "SFRP1"
find_location
if(length(xx) > 0) {
  for(i in 1:length(xx)){
    if(xx[[i]] == "SFRP1"){
      location <- yy[[i]]
      gene <- xx[[i]]
      #print(xx[[i]])
      print(paste(i,rownames(xx[[i]]),gene,names(location),location,sep=":"))
    }
  }
}


head(xx)
xx[[1]]
gname = "RASSF1"
gname = "CDKN2A"
names(yy[[14423]])

index <- NULL
for(i in 1:length(rownames(as.matrix(xx[xx[]==gname])))){
  index <- c(index, which(rownames(as.matrix(yy))==rownames(as.matrix(xx[xx[]==gname]))[i], arr.ind=TRUE))
}
index
length(index)
for(i in 1:length(index)){
#  p <- paste(i, names(yy[[i]]), sep=":")
#  p <- paste(i, names(yy[[i]]), sep=":")
  #yy[rownames(as.matrix(xx[xx[]==gname]))][1]
  num_loc <- length(p<- paste(i, as.matrix(xx[xx[]==gname])[i],rownames(as.matrix(xx[xx[]==gname]))[i], names(yy[[i]]), sep=":"))
  #p<- paste(i, as.matrix(xx[xx[]==gname])[i],rownames(as.matrix(xx[xx[]==gname]))[i], yy[rownames(as.matrix(xx[xx[]==gname]))], sep=":")
  if(num_loc>1){
    for(j in 1:num_loc){
      #print(paste(i, yy[rownames(as.matrix(xx[xx[]==gname]))][j], sep=":"))
      print(paste(as.matrix(xx[xx[]==gname])[i],rownames(as.matrix(xx[xx[]==gname]))[i], names(yy[[i]]), yy[rownames(as.matrix(xx[xx[]==gname]))][j][1], sep=":"))
      #print(paste(i, as.matrix(xx[xx[]==gname])[i],rownames(as.matrix(xx[xx[]==gname]))[i], names(yy[[i]]), yy[rownames(as.matrix(xx[xx[]==gname]))][j], sep=":"))
      
      #paste(as.matrix(paste(i, as.matrix(xx[xx[]==gname])[34],rownames(as.matrix(xx[xx[]==gname]))[34], names(yy[[34]]), sep=":")),,sep=":")
    }
  }else{
    p<- paste(i, as.matrix(xx[xx[]==gname])[i],rownames(as.matrix(xx[xx[]==gname]))[i], names(yy[[i]]), yy[rownames(as.matrix(xx[xx[]==gname]))][i],sep=":")
    print(p)
  }
}



paste(21,as.matrix(paste(i, as.matrix(xx[xx[]==gname])[34],rownames(as.matrix(xx[xx[]==gname]))[34], names(yy[[34]]), sep=":")),3213,sep=":")
yy[rownames(as.matrix(xx[xx[]==gname]))][1]
rownames(as.matrix(xx[xx[]==gname]))[1]
as.matrix(xx[xx[]==gname])[1]
which(rownames(as.matrix(yy))==rownames(as.matrix(xx[xx[]==gname]))[i], arr.ind=TRUE)


paste(as.matrix(xx[xx[]==gname]),rownames(as.matrix(xx[xx[]==gname])), yy[rownames(as.matrix(xx[xx[]==gname]))],sep=":")

names(yy[[14464]])

gname = "CDKN2A"
for(i in 1:length(rownames(as.matrix(xx[xx[]==gname])))){
  index <- which(rownames(as.matrix(yy))==rownames(as.matrix(xx[xx[]==gname]))[i], arr.ind=TRUE)
  chrom <- names(yy[[index]])
#  print(chrom)
  dat <- paste(as.matrix(xx[xx[]==gname]),rownames(as.matrix(xx[xx[]==gname])), chrom, yy[rownames(as.matrix(xx[xx[]==gname]))],sep=":")
  #dat <- paste(as.matrix(xx[xx[]==gname]),rownames(as.matrix(xx[xx[]==gname])), yy[rownames(as.matrix(xx[xx[]==gname]))],sep=":")
  print(dat)
}



yy[rownames(as.matrix(xx[xx[]=="CDKN2A"]))]
rownames(as.matrix(xx[xx[]=="CDKN2A"]))
which(rownames(yy[])==rownames(as.matrix(xx[xx[]=="CDKN2A"])), arr.ind=TRUE)
which(as.matrix(rownames(as.matrix(xx[xx[]=="CDKN2A"]))), arr.ind=TRUE)

yy$n
yy$cg04026675
head(yy)
head(rownames(yy[]))
paste(as.matrix(xx[xx[]=="CDKN2A"]),rownames(as.matrix(xx[xx[]=="CDKN2A"])), sep=":")
paste(as.matrix(xx[xx[]=="CDKN2A"]),rownames(as.matrix(xx[xx[]=="CDKN2A"])), yy[rownames(as.matrix(xx[xx[]=="CDKN2A"]))],sep=":")


# which(rownames(as.matrix(yy))==(names(yy[rownames(as.matrix(xx[xx[]=="CDKN2A"]))])[1]), arr.ind=TRUE)
# names(yy[rownames(as.matrix(xx[xx[]=="CDKN2A"]))])


head(as.matrix(xx))
rownames(as.matrix(xx[xx[]=="CDKN2A"]))
yy["cg04026675"]
as.matrix(xx[xx[]=="CDKN2A"])
xx[xx[]=="CDKN2A"]

#P16
rownames(as.matrix(xx[xx[]=="CDKN2A"]))

#SFRP1
rownames(as.matrix(xx[xx[]=="SFRP1"]))
#TFPI-2
rownames(as.matrix(xx[xx[]=="TFPI2"]))
#GATA4
rownames(as.matrix(xx[xx[]=="GATA4"]))
#GATA5
rownames(as.matrix(xx[xx[]=="GATA5"]))
#APC
rownames(as.matrix(xx[xx[]=="APC"]))
#CHFR
rownames(as.matrix(xx[xx[]=="CHFR"]))
#RASSF1A
rownames(as.matrix(xx[xx[]=="RASSF1"]))
#HIN1
rownames(as.matrix(xx[xx[]=="SCGB3A1"]))
#O6-MGMT
paste("(?^i:",rownames(as.matrix(xx[xx[]=="MGMT"])),")",sep="")
#MLH1
as.matrix(rownames(as.matrix(xx[xx[]=="MLH1"])))
(?^i:ah\D+a)

rownames(as.matrix(xx[xx[]=="CDKN2A"]))

subset <- as.matrix(c(
  #                 "#P16", paste("(?^i:",rownames(as.matrix(xx[xx[]=="CDKN2A"])),")",sep=""), 
  #                 "#SFRP1", paste("(?^i:",rownames(as.matrix(xx[xx[]=="SFRP1"])),")",sep=""),
  #                 "#TFP12", paste("(?^i:",rownames(as.matrix(xx[xx[]=="TFP12"])),")",sep=""),
  #                 "#GATA4", paste("(?^i:",rownames(as.matrix(xx[xx[]=="GATA4"])),")",sep=""),
  #                 "#GATA5", paste("(?^i:",rownames(as.matrix(xx[xx[]=="GATA5"])),")",sep=""),
  #                 "#APC", paste("(?^i:",rownames(as.matrix(xx[xx[]=="APC"])),")",sep=""),
  #                 "#CHFR", paste("(?^i:",rownames(as.matrix(xx[xx[]=="CHFR"])),")",sep=""),
  #                 "#RASSF1A", paste("(?^i:",rownames(as.matrix(xx[xx[]=="RASSF1"])),")",sep=""),
  #                 "#HIN1", paste("(?^i:",rownames(as.matrix(xx[xx[]=="SCGB3A1"])),")",sep=""),
  #                 "#O6-MGMT", paste("(?^i:",rownames(as.matrix(xx[xx[]=="MGMT"])),")",sep=""),
  #                 "#MLH1", paste("(?^i:",rownames(as.matrix(xx[xx[]=="MLH1"])),")",sep="")
                  "#P16", rownames(as.matrix(xx[xx[]=="CDKN2A"])), 
                  "#SFRP1", rownames(as.matrix(xx[xx[]=="SFRP1"])),
                  "#TFP12", rownames(as.matrix(xx[xx[]=="TFP12"])),
                  "#GATA4", rownames(as.matrix(xx[xx[]=="GATA4"])),
                  "#GATA5", rownames(as.matrix(xx[xx[]=="GATA5"])),
                  "#APC", rownames(as.matrix(xx[xx[]=="APC"])),
                  "#CHFR", rownames(as.matrix(xx[xx[]=="CHFR"])),
                  "#RASSF1A", rownames(as.matrix(xx[xx[]=="RASSF1"])),
                  "#HIN1", rownames(as.matrix(xx[xx[]=="SCGB3A1"])),
                  "#O6-MGMT", rownames(as.matrix(xx[xx[]=="MGMT"])),
                  "#MLH1", rownames(as.matrix(xx[xx[]=="MLH1"]))
  
))
subset
write.table(subset, file="/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/1197_SU2C HM450 Data/subset/450probe-genes.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
# p16, SFRP 1, TFPI-2 GATA4, GATA5, HIN1, O6-MGMT, APC, CHFR, RASSF1A, and hMLH1.”
# sfrp1 <- new.datfil["SFRP1",]
# tfpi2 <- new.datfil["TFPI2",]
# gata4 <- new.datfil["GATA4",]
# gata5 <- new.datfil["GATA5",]
# apc <- new.datfil["APC",]
# chfr <- new.datfil["CHFR",]
# rassf1 <- new.datfil["RASSF1",]
# hin1 <- new.datfil["SCGB3A1",]
# mgmt <- new.datfil["MGMT",]
#MLH1
