options(digits=16)
tab <- readLines('HD128620_ALMA.dat')
tab  <- gsub('\t','',tab)
out <- c()
for(j in 1:length(tab)){
    tmp <- unlist(strsplit(tab[j],' '))
    tmp <- tmp[tmp!='']
    out <- rbind(out,as.numeric(tmp))
}
jd <- out[,6]+2400000.5
ra <- ((14+39/60+out[,7]/3600)/24)*360
dra <- out[,8]/3600/24*(360*3600)*1000#mas
dec <- -(60+49/60+out[,9]/3600)
ddec <- out[,10]*1e3#mas
lambda <- 299792458/(out[,2]*1e9)*1e6#micrometer
write.table(cbind(jd,ra,dra,dec,ddec,lambda),file='HD128620_ALMA.astro',quote=FALSE,row.names=FALSE,col.names=c('JD','ra','era','dec','edec','lambda'))

