
head(Annot.probelevel)
Annot.probelevel$GeneName[grep('SLC22A1$',Annot.probelevel$GeneName)]

t<-rownames(Annot.probelevel)[which(Annot.probelevel$GeneName=='SLC22A1')]
t
SLC22A1.mock<-cbind(Annot.probelevel[t,c(2,3,5,6)],Brbeta.mock[t,],Cobeta.mock[t,],Ovbeta.mock[t,])
SLC22A1.aza<-cbind(Annot.probelevel[t,c(2,3,5,6)],Brdbeta[t,],Codbeta[t,],Ovdbeta[t,])

SLC22A1.annot<-Annot.probelevel[t,5:6]
SLC22A1.annot[SLC22A1.annot=='']=NA
unique(SLC22A1.annot[,1])
unique(SLC22A1.annot[,2])
row.annot<-gsub('Island','red',as.matrix(SLC22A1.annot))
colnames(row.annot)<-c('CPG_ISLAND','')
row.annot[row.annot!='red']<-'white'
row.annot[is.na(row.annot)]='white'
row.annot

png('./SLC22A1 in breast colon and ovarian.png',height=10,width=15,units="in",res=150)
par(oma=c(8,1,1,3))
dat.subset<-SLC22A1.mock[,5:ncol(SLC22A1.mock)]
p<-heatmap.3(as.matrix(dat.subset), 
hclust=function(x) hclust(x,method="ward"),
#Colv=NULL,
scale="none", 
RowSideColors=t(as.matrix(row.annot[rownames(dat.subset),])),
col=colorRampPalette(c("white","purple"))(89),
density.info="none",
keysize=0.8,
#symkey=F
)
dev.off()


png('./SLC22A1 aza delta beta in breast colon and ovarian.png',height=10,width=20,units="in",res=150)
par(oma=c(8,1,1,3))
dat.subset<-SLC22A1.aza[,5:ncol(SLC22A1.aza)]
q<-heatmap.3(as.matrix(dat.subset), 
hclust=function(x) hclust(x,method="ward"),
#Colv=NULL,
scale="none", 
RowSideColors=t(as.matrix(row.annot[rownames(dat.subset),])),
col=colorRampPalette(c('green',"white","purple"))(89),
density.info="none",
keysize=0.8,
symkey=T
)
dev.off()

---------------------------------------
t<-rownames(Annot.probelevel)[which(Annot.probelevel$GeneName=='SLC22A2')]
t
SLC22A2.mock<-cbind(Annot.probelevel[t,c(2,3,5,6)],Brbeta.mock[t,],Cobeta.mock[t,],Ovbeta.mock[t,])
SLC22A2.aza<-cbind(Annot.probelevel[t,c(2,3,5,6)],Brdbeta[t,],Codbeta[t,],Ovdbeta[t,])

SLC22A2.annot<-Annot.probelevel[t,5:6]
SLC22A2.annot[SLC22A2.annot=='']=NA
unique(SLC22A2.annot[,1])
unique(SLC22A2.annot[,2])
row.annot<-gsub('Island','red',as.matrix(SLC22A2.annot))
colnames(row.annot)<-c('CPG_ISLAND','')
row.annot[row.annot!='red']<-'white'
row.annot[is.na(row.annot)]='white'
row.annot


png('./SLC22A2 in breast colon and ovarian.png',height=10,width=15,units="in",res=150)
par(oma=c(8,1,1,3))
dat.subset<-SLC22A2.mock[,5:ncol(SLC22A2.mock)]
p<-heatmap.3(as.matrix(dat.subset), 
hclust=function(x) hclust(x,method="ward"),
#Colv=NULL,
scale="none", 
RowSideColors=t(as.matrix(row.annot[rownames(dat.subset),])),
col=colorRampPalette(c("white","purple"))(89),
density.info="none",
keysize=0.8,
#symkey=F
)
dev.off()


png('./SLC22A2 aza delta beta in breast colon and ovarian.png',height=10,width=20,units="in",res=150)
par(oma=c(8,1,1,3))
dat.subset<-SLC22A2.aza[,5:ncol(SLC22A2.aza)]
q<-heatmap.3(as.matrix(dat.subset), 
hclust=function(x) hclust(x,method="ward"),
#Colv=NULL,
scale="none", 
RowSideColors=t(as.matrix(row.annot[rownames(dat.subset),])),
col=colorRampPalette(c('green',"white","purple"))(89),
density.info="none",
keysize=0.8,
symkey=T
)
dev.off()

---------------------------------------
t<-rownames(Annot.probelevel)[which(Annot.probelevel$GeneName=='SLC22A3')]
t
SLC22A3.mock<-cbind(Annot.probelevel[t,c(2,3,5,6)],Brbeta.mock[t,],Cobeta.mock[t,],Ovbeta.mock[t,])
SLC22A3.aza<-cbind(Annot.probelevel[t,c(2,3,5,6)],Brdbeta[t,],Codbeta[t,],Ovdbeta[t,])


SLC22A3.annot<-Annot.probelevel[t,5:6]
SLC22A3.annot[SLC22A3.annot=='']=NA
unique(SLC22A3.annot[,1])
unique(SLC22A3.annot[,2])
row.annot<-gsub('Island','red',as.matrix(SLC22A3.annot))
colnames(row.annot)<-c('CPG_ISLAND','')
row.annot[row.annot!='red']<-'white'
row.annot[is.na(row.annot)]='white'
row.annot

png('./SLC22A3 in breast colon and ovarian.png',height=10,width=15,units="in",res=150)
par(oma=c(8,1,1,3))
dat.subset<-SLC22A3.mock[,5:ncol(SLC22A3.mock)]
p<-heatmap.3(as.matrix(dat.subset), 
hclust=function(x) hclust(x,method="ward"),
#Colv=NULL,
scale="none", 
RowSideColors=t(as.matrix(row.annot[rownames(dat.subset),])),
col=colorRampPalette(c("white","purple"))(89),
density.info="none",
keysize=0.8,
#symkey=F
)
dev.off()

png('./SLC22A3 aza delta beta in breast colon and ovarian.png',height=10,width=20,units="in",res=150)
par(oma=c(8,1,1,3))
dat.subset<-SLC22A3.aza[,5:ncol(SLC22A3.aza)]
q<-heatmap.3(as.matrix(dat.subset), 
hclust=function(x) hclust(x,method="ward"),
#Colv=NULL,
scale="none", 
RowSideColors=t(as.matrix(row.annot[rownames(dat.subset),])),
col=colorRampPalette(c('green',"white","purple"))(89),
density.info="none",
keysize=0.8,
symkey=T
)
dev.off()
----------------------------------
annotation<-Annot[c(rownames(SLC22A1.annot),rownames(SLC22A2.annot),rownames(SLC22A3.annot)),]
dim(annotation)
unique(rownames(annotation))
write.csv(annotation,file='SLC22A1 2 3 probe annotation.csv')
