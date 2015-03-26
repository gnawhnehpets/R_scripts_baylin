setwd("/home/steve/Downloads/")
dat <- read.delim("MCF7_450k_islands.txt")
dat <- dat[,c(1,2,4,5)]
x <- dat
head(x)
# Probe Gene MCF10a       MCF7
# 1 cg03630821 A1BG  0.939 0.94428220
# 2 cg14222245 A1BG  0.852 0.88098800
# 3 cg03123289 A1BG  0.129 0.34758480
# 4 cg11001216 A1BG  0.096 0.07003345
# 5 cg10734734 A1BG  0.019 0.01909087
# 6 cg22286978 A1BG  0.381 0.01209927
tmp <- setNames(strsplit(as.character(x$Gene), split=';'), x$Probe)
head(tmp)
table <- data.frame(probe = rep(names(tmp), sapply(tmp, length)), genename = unlist(tmp), MCF10a = rep(x[,3], sapply(tmp, length)), MCF7 = rep(x[,4], sapply(tmp, length)), row.names = NULL)
final <- unique(table)
data.frame(head(x, n=100), head(table, n=100), head(final, n=100))
data.frame(head(x, n=50), head(final, n=50))


head(x, n=20)
head(table2, n=50)

length(tmp)
dim(x)
dim(tmp)
head(table2)




xx <- lapply(1:nrow(x), FUN=function(i){
  splitValues <- strsplit(as.character(x$Gene[i]), ";")
  lengthSpliValues = length(splitValues[[1]])
  probe = rep(x$Probe[i], times=lengthSpliValues)
  betaMCF10a = rep(x$MCF10a[i], times=lengthSpliValues)
  betaMCF7 = rep(x$MCF7[i], times=lengthSpliValues)  
  probesAndGene <- cbind.data.frame(Probe=probe, Gene=unlist(splitValues), 
                                    MCF10a=betaMCF10a, MCF7=betaMCF7)
  return(probesAndGene)
})
rbindedDat <- do.call(rbind, xx)