options(digits=20)
stars <- c('HD128620','HD128621')
for(star in stars){
tab <- readLines(paste0(star,'_ALMA.dat'))
tab  <- gsub('\t','',tab)
decs <- out <- c()
for(j in 1:length(tab)){
    tmp <- unlist(strsplit(tab[j],' '))
    tmp <- tmp[tmp!='']
    out <- rbind(out,as.numeric(tmp))
}
jd <- out[,6]+2400000.5
epoch.outlier <- 2457007.9699005
ind <- which.min(abs(jd-epoch.outlier))
if(length(ind)>0){
    out <- out[-ind,]
    jd <- out[,6]+2400000.5
}
if(star==stars[1]){
    pp <- cbind(out[,c(6,7,8,9,10)])
}else{
    qq <- cbind(out[,c(6,7,8,9,10)])
}
ra <- ((14+39/60+out[,7]/3600)/24)*360#deg
dra <- out[,8]/3600/24*(360*3600)*1000#mas
dec <- -(60+49/60+out[,9]/3600)
decs <- cbind(decs,dec)
ddec <- out[,10]*1e3#mas
lambda <- 299792458/(out[,2]*1e9)*1e6#micrometer
tmp <- cbind(jd,ra,dra,dec,ddec,lambda)
#ind <- which.min(abs(jd-epoch.outlier))
#write.table(tmp[-ind,],file='../HD128621v1/HD128621_ALMA.abs',quote=FALSE,row.names=FALSE,col.names=c('JD','ra','era','dec','edec','lambda'))
write.table(tmp,file=paste0('../',star,'/',star,'_ALMA.abs'),quote=FALSE,row.names=FALSE,col.names=c('JD','ra','era','dec','edec','lambda'))
write.table(tmp,file=paste0('../',star,'v1/',star,'_ALMA.abs'),quote=FALSE,row.names=FALSE,col.names=c('JD','ra','era','dec','edec','lambda'))
}
jd <- pp[,1]+2400000.5
ddec <-(qq[,4]-pp[,4])*1e3
eddec <-sqrt(qq[,5]^2+pp[,5]^2)*1e3
dec <- rowSums(decs)/2*pi/180
dra <- (qq[,2]-pp[,2])/24*360*1e3*cos(dec)
edra <- sqrt(qq[,3]^2+pp[,3]^2)/24*360*1e3
star <- 'HD128621'
out <- cbind(jd,dra,edra,ddec,eddec)
write.table(out,file=paste0('../',star,'/',star,'_ALMA.rel'),quote=FALSE,row.names=FALSE,col.names=c('JD','ra','era','dec','edec'))
write.table(out,file=paste0('../',star,'v1/',star,'_ALMA.rel'),quote=FALSE,row.names=FALSE,col.names=c('JD','ra','era','dec','edec'))
