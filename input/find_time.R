dirs <- list.dirs('.')
for(dir in dirs){
     fs <- list.files(dir,pattern='rv$',full.name=TRUE)
     for(f in fs){
     	   tab <- read.table(f,header=TRUE)
	   write.table(cbind(tab[,1],0),file=paste0(dir,'.tim'),row.names=FALSE,quote=FALSE,col.names=FALSE)
     }
}
