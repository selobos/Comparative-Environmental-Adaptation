#Set directory
setwd("~Desktop/Comparative-Landscape-Genomics/3_GEA/3_1_BayPass/GA/Core")

###Calculate means from initial baypass runs on core model and produce covariance matrix
myList<- vector(mode = "list", length = 5) 

for (i in 1:5) {
  myList[[i]]<-as.matrix(read.table(paste("ga_core.",i,"_mat_omega.out",sep="")))
}

#mean omega matrix
median <- sapply(1:ncol(myList[[1]]), function(j) {apply(do.call(cbind,lapply(myList,`[`,,j)), 1, median)})

write.table(median, file=paste0("ga_core.median_mat_omega.out"), row.names=FALSE, col.names=FALSE, sep="\t")
