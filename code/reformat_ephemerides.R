de <- commandArgs(trailingOnly=TRUE)
if(length(de)==0){
    de <- '430t'
}
dir <- paste0('../data/de',de)
fs <- list.files(dir,pattern=paste0('ascp.+',de))
fh <- list.files(dir,pattern=paste0('header.',de),full.name=TRUE)
indH <- sort(fh,decreasing=TRUE,index.return=TRUE)$ix[1]
cat('read header:',fh[indH],'\n')
header <- readLines(fh[indH])
yy0 <- gsub(paste0('ascp|\\.',de),'',fs)
yy <- as.numeric(yy0)
dyr <- diff(yy)[1]
#ind1 <- which(yy>max(1950,min(yy)))[1]-1
ind1 <- which(yy>max(1650,min(yy)))[1]-1
ind2 <- which(yy>=min(max(yy),2050))[1]
yrs <- yy[ind1:ind2]
out <- c()
for(k in ind1:ind2){
    str <- readLines(paste0(dir,'/ascp',yy0[k],'.',de))
    Ncoeff <- as.integer(gsub('.+NCOEFF=','',header[1]))
    ind <- c(grep(paste0('  ',Ncoeff),str),length(str)+1)
    for(j in 1:(length(ind)-1)){
        s <- unlist(strsplit(gsub('D','e',paste(str[(ind[j]+1):(ind[j+1]-1)],collapse=' ')),' '))
        s <- as.numeric(s[s!=''])
        out <- rbind(out,s)
    }
}
fout <- paste0(dir,'/','DE',de,'_J',min(yrs),'_J',max(yrs)+dyr,'.dat')
cat(fout,'\n')
write.table(out,file=fout,quote=FALSE,row.names=FALSE,col.names=FALSE)
