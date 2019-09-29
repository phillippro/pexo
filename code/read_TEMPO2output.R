if(!exists('UNITS')){
uni <- 'tcb'
}else{
uni <- tolower(UNITS)
}
uni <- 'tdb'
##replace the following file by the tempo2 output file, polyco_new.dat, with the correct path
#fin <- paste0('~/Documents/tempo2/polyco_',uni,'0.dat')
fin <- paste0('~/Documents/tempo2/polyco_',uni,'.dat')
#fin <- paste0('~/Documents/tempo2/polyco_tdb.dat')
cat('input polyco emulation file:\n',fin,'\n')
tab <- readLines(fin)
N <- length(tab)
data <- c()
k <- 0
cs <- c('star','date','utc','mjd','dm','doppler','rms','rphase','val','sitename','tspan','ncoeff.obsFreq.Phase',paste0('coeff',1:12))
for(j in 1:N){
    if(grepl('^10700',tab[j])){
        if(j>1){
            data <- rbind(data,tmp)
        }
        tmp <- c()
        out <- unlist(strsplit(tab[j],' '))
        out <- out[out!='']
        tmp <- c(tmp,out)
    }else{
        out <- unlist(strsplit(tab[j],' '))
        out <- out[out!='']
        tmp <- c(tmp,out)
    }
}
data <- rbind(data,tmp)
colnames(data) <- cs
write.table(data,file='../data/tempo2_emulation.dat',row.names=FALSE,quote=FALSE)
#cat('doppler=',as.numeric(data[,'doppler'])*1e-4,'\n')
