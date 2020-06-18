file <- commandArgs(trailingOnly=TRUE)
tab <- read.table(file,header=TRUE)
target <- gsub('.+/|_.+','',file)
paul <- read.table(paste0('../input/',target,'/',target,'_PFSraw.rv'),header=TRUE)
fout <- paste0('../results/',target,'_PFSpexo3.dat')
cat(fout,'\n')
###tauE as epoch
write.table(cbind(rowSums(tab[,3:4]),tab[,5],paul[,3]),file=fout,quote=FALSE,row.names=FALSE,col.names=c('BJD','RV','eRV'))
###BJDtdb as epoch
#write.table(cbind(rowSums(tab[,1:2]),tab[,3],paul[,3]),file=fout,quote=FALSE,row.names=FALSE,col.names=c('BJD','RV','eRV'))
