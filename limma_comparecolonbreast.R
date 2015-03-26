rm(list=ls())
library(scales)
library(limma)
load("/home/steve/Desktop/analysis/cohortoneandtwo/att3/huili TN pre post not log2 scaled.rda")
load("/home/steve/Desktop/analysis/cohortoneandtwo/att3/huili mat.luminal.pre.post log2 scaled.rda")

colnames(TNtrial)               #no log2fc # this R object has JH016pre @ col#6,22   # this R object has JH016post @ col#4,12
colnames(mat.luminal.pre.post)  #log2      # this R object has JH016pre @ col#9      # this R object has JH016post @ col#10
x <- as.matrix(TNtrial[,4])
y <- as.matrix(TNtrial[,12])
# colnames(TNtrial)[6]
# colnames(TNtrial)[22]
# colnames(mat.luminal.pre.post)[9]
rownames(x) <- TNtrial[,2]


#y <- as.matrix(mat.luminal.pre.post[,9])
#y <- mat.luminal.pre.post[,17]
r.cor <- cor.test(log2(x),log2(y))
r.cor
dev.off()
plot(log2(x), log2(y),main="breast biopsies\npre v. pre comparison",xlab="TNtrial[,22]",ylab="TNtrial[,6]",col='blue',cex=1.5,type="n")
points(log2(x), log2(y), col=alpha("blue",0.25),pch=21)
#points(log2(x), log2(y), bg="lightblue",col=1,pch=21)
text(log2(x), log2(y), label=rownames(x),pos=1,cex=0.4)

#abline(lm( log2(TNtrial[,6]) ~ log2(TNtrial[,22])))
abline(lm( log2(y) ~ log2(x) ))
#abline(lm( log2(TNtrial[,22]) ~ log2(TNtrial[,4])))

length(r.cor$p.value)
plot(r.cor)

library(gplots)
#r.cor <- cor(r)
layout(matrix(c(1,1,1,1,1,1,1,1,2,2), 5, 2, byrow = TRUE))
par(oma=c(5,7,1,1))
cx <- rev(colorpanel(25,"yellow","black","blue"))
leg <- seq(min(r.cor,na.rm=T),max(r.cor,na.rm=T),length=10)
image(r.cor,main="Correlation plot colon biopsy data",axes=F,col=cx)
axis(1,at=seq(0,1,length=ncol(r.cor)),label=dimnames(r.cor)[[2]],cex.axis=0.9,las=2)
axis(2,at=seq(0,1,length=ncol(r.cor)),label=dimnames(r.cor)[[2]],cex.axis=0.9,las=2)
par(mar=rep(2,4))
image(as.matrix(leg),col=cx,axes=T)
tmp <- round(leg,2)

