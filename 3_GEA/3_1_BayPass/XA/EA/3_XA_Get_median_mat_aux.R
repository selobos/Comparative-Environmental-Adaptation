#Set directory
setwd("~Desktop/Comparative-Landscape-Genomics/3_GEA/3_1_BayPass/XA/EA")

###Calculate means for aux runs and output a csv file with only those above BF20
myList<- vector(mode = "list", length = 5)

for (i in 1:5) {
  myList[[i]]<-as.matrix(read.table(paste("xa_aux.",i,"_summary_betai.out",sep=""), header=T))
}

bfmedian <- sapply(1:ncol(myList[[1]]), function(j) {apply(do.call(cbind,lapply(myList,`[`,,j)), 1, median)})
colnames(bfmedian) <-colnames(myList[[1]])

bfmedian<-as.data.frame(bfmedian)
write.csv(bfmedian, file=paste0("xa_aux_median_summary_betai.out.csv"))

uniq_snpss=bfmedian[,"MRK"][which(bfmedian[,"BF.dB."]>20)] 