out <- list()
fout <- paste0('../results/',star,'_comparison.pdf')
cat(fout,'\n')
pdf(fout,16,16)
par(mfrow=c(2,2))
                                        #    tab <- read.table('~/Documents/projects/dwarfs/data/standard/HD128621/ustar.rv')
fin <- paste0('~/Documents/projects/dwarfs/data/standard/',star,'/',star,'_D12.dat')
tab <- read.table(fin,header=TRUE)
t1 <- tab[,1]+2400000
t2 <- rowSums(bjd.tdb)
yr1 <- time_Jd2yr(cbind(t1%/%1,t1%%1))
yr2 <- time_Jd2yr(bjd.tdb)
y1 <- tab[,2]
y2 <- RvFull$kt/1e3+0.322
orb <- approxfun(yr2,y2)
plot(yr1,y1,xlab='BJD',ylab='RV [km/s]')
lines(yr2,y2,col='red')
####relativistic redshift
plot(yr1,(y1-orb(yr1))*1e3,xlab='BJD',ylab='RV [m/s]',main='HARPS D12')
lines(yr2,-RvFull$gt,col='red')
lines(yr2,-RvFull$st,col='blue')

####another fit with TERRA reduction
data <- read.table('~/Documents/projects/dwarfs/data/aperture/HD128621_ew3_v2/HD128621_HARPS.dat',header=TRUE)
tmp <- adjust.data(data[,1],data[,2],yr2,y2,factor=1e-3)
yr1 <- tmp[,1]
y1 <- tmp[,2]
plot(yr1,y1,xlab='BJD',ylab='RV [km/s]',main='HARPS-TERRA')
lines(yr2,y2,col='red')
###residual
dy <- (y1-orb(yr1))*1e3
plot(yr1,dy,xlab='BJD',ylab='RV [m/s]')
lines(yr2,RvFull$gt,col='red')
lines(yr2,RvFull$st,col='blue')
out[['HARPS']] <- cbind(yr1,y1,dy,data[,3])

###another fit with UVES data
data <- read.table('~/Documents/projects/dwarfs/data/aperture/alphaCen/alphaCenB_chi.dat',header=TRUE)
tmp <- adjust.data(data[,1],data[,2],yr2,y2,factor=1e-3)
yr1 <- tmp[,1]
y1 <- tmp[,2]
plot(yr1,y1,xlab='BJD',ylab='RV [km/s]',main='CHIRON')
lines(yr2,y2,col='red')
###residual
dy <- (y1-orb(yr1))*1e3
plot(yr1,dy,xlab='BJD',ylab='RV [m/s]')
lines(yr2,RvFull$gt,col='red')
lines(yr2,RvFull$st,col='blue')
out[['CHIRON']] <- cbind(yr1,y1,dy,data[,3])
                                        #    abline(v=c(1988.89,2001.73),col='blue')

data <- read.table('~/Documents/projects/dwarfs/data/aperture/alphaCen/alphaCenB_es.dat',header=TRUE)
tmp <- adjust.data(data[,1],data[,2],yr2,y2,factor=1e-3)
yr1 <- tmp[,1]
y1 <- tmp[,2]
plot(yr1,y1,xlab='BJD',ylab='RV [km/s]',main='ES')
lines(yr2,y2,col='red')
###residual
dy <- (y1-orb(yr1))*1e3
plot(yr1,dy,xlab='BJD',ylab='RV [m/s]')
lines(yr2,RvFull$gt,col='red')
lines(yr2,RvFull$st,col='blue')
out[['ES']] <- cbind(yr1,y1,dy,data[,3])

data <- read.table('~/Documents/projects/dwarfs/data/aperture/HD128621/HD128621_UVES_unbinned.dat',header=TRUE)
tmp <- adjust.data(data[,1],data[,2],yr2,y2,factor=1e-3)
yr1 <- tmp[,1]
y1 <- tmp[,2]
plot(yr1,y1,xlab='BJD',ylab='RV [km/s]',main='UVES')
lines(yr2,y2,col='red')
###residual
dy <- (y1-orb(yr1))*1e3
plot(yr1,dy,xlab='BJD',ylab='RV [m/s]')
lines(yr2,RvFull$gt,col='red')
lines(yr2,RvFull$st,col='blue')
out[['UVES']] <- cbind(yr1,y1,dy,data[,3])

data <- read.table('~/Documents/projects/dwarfs/data/other_sets/alpha_cen_old/HD128621_aat.dat',header=TRUE)
tmp <- adjust.data(data[,1],data[,2],yr2,y2,factor=1e-3)
yr1 <- tmp[,1]
y1 <- tmp[,2]
ind <- which(yr1<2009)
yr1 <- yr1[ind]
y1 <- y1[ind]
data <- data[ind,]
plot(yr1,y1,xlab='BJD',ylab='RV [km/s]',main='AAT')
lines(yr2,y2,col='red')
###residual
dy <- (y1-orb(yr1))*1e3
plot(yr1,dy,xlab='BJD',ylab='RV [m/s]')
lines(yr2,RvFull$gt,col='red')
lines(yr2,RvFull$st,col='blue')
out[['AAT']] <- cbind(yr1,y1,dy,data[,3])

###
data <- read.table('~/Documents/projects/dwarfs/data/other_sets/alpha_cen_old/LCb.dat',header=TRUE)
tmp <- adjust.data(data[,1],data[,2],yr2,y2,factor=1e-3)
yr1 <- tmp[,1]
y1 <- tmp[,2]
plot(yr1,y1,xlab='BJD',ylab='RV [km/s]',main='LC')
lines(yr2,y2,col='red')
###residual
dy <- (y1-orb(yr1))*1e3
plot(yr1,dy,xlab='BJD',ylab='RV [m/s]')
lines(yr2,RvFull$gt,col='red')
lines(yr2,RvFull$st,col='blue')
out[['LC']] <- cbind(yr1,y1,dy,data[,3])

###
data <- read.table('~/Documents/projects/dwarfs/data/other_sets/alpha_cen_old/coralieb.dat',header=TRUE)
tmp <- adjust.data(data[,1],data[,2],yr2,y2,factor=1e-3)
yr1 <- tmp[,1]
y1 <- tmp[,2]
plot(yr1,y1,xlab='BJD',ylab='RV [km/s]',main='CORALIE')
lines(yr2,y2,col='red')
###residual
dy <- (y1-orb(yr1))*1e3
plot(yr1,dy,xlab='BJD',ylab='RV [m/s]')
out[['CORALIE']] <- cbind(yr1,y1,dy,data[,3])

###
data <- read.table('~/Documents/projects/dwarfs/data/other_sets/alpha_cen_old/VLCb.dat',header=TRUE)
tmp <- adjust.data(data[,1],data[,2],yr2,y2,factor=1e-3)
yr1 <- tmp[,1]
y1 <- tmp[,2]
plot(yr1,y1,xlab='BJD',ylab='RV [km/s]',main='VLC')
lines(yr2,y2,col='red')
###residual
dy <- (y1-orb(yr1))*1e3
plot(yr1,dy,xlab='BJD',ylab='RV [m/s]')
lines(yr2,RvFull$gt,col='red')
lines(yr2,RvFull$st,col='blue')
out[['VLC']] <- cbind(yr1,y1,dy,data[,3])
dev.off()
fout <- gsub('.pdf','.Robj',fout)
cat(fout,'\n')
zgt <- RvFull$gt
zst <- RvFull$st
save(out,yr2,y2,zgt,zst,file=fout)
