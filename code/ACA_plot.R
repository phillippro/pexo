out <- list()
fout <- paste0('../results/',star,'_comparison_',Par$RefType,'.pdf')
cat(fout,'\n')
pdf(fout,16,16)
par(mfrow=c(2,2))
fin <- '~/Documents/projects/dwarfs/data/other_sets/alpha_cen_old/HARPSa.dat'
tab <- read.table(fin,header=TRUE)
t1 <- tab[,1]
if(max(tab[,1])<2400000) t1 <- tab[,1]+2400000
yr1 <- time_Jd2yr(cbind(t1%/%1,t1%%1))
y1 <- tab[,2]

t2 <- rowSums(bjd.tdb)
yr2 <- time_Jd2yr(bjd.tdb)
y2 <- rv.out$kt/1e3+0.322

orb <- approxfun(yr2,y2)
plot(yr1,y1,xlab='BJD',ylab='RV [km/s]')
lines(yr2,y2,col='red')
####relativistic redshift
plot(yr1,(y1-orb(yr1))*1e3,xlab='BJD',ylab='RV [m/s]',main='HARPS D12')
lines(yr2,-rv.out$gt,col='red')
lines(yr2,-rv.out$st,col='blue')

####another fit with TERRA reduction
data <- read.table('~/Documents/projects/dwarfs/data/aperture/HD128620/HD128620_HARPS.dat',header=TRUE)
tmp <- adjust.data(data[,1],data[,2],yr2,y2,factor=1e-3)
yr1 <- tmp[,1]
y1 <- tmp[,2]
plot(yr1,y1,xlab='BJD',ylab='RV [km/s]',main='HARPS-TERRA')
lines(yr2,y2,col='red')
###residual
dy <- (y1-orb(yr1))*1e3
plot(yr1,dy,xlab='BJD',ylab='RV [m/s]')
lines(yr2,rv.out$gt,col='red')
lines(yr2,rv.out$st,col='blue')
out[['HARPS']] <- cbind(yr1,y1,dy,data[,3])

###another fit with UVES data
data <- read.table('~/Documents/projects/dwarfs/data/aperture/alphaCen/alphaCenA_chi.dat',header=TRUE)
tmp <- adjust.data(data[,1],data[,2],yr2,y2,factor=1e-3)
yr1 <- tmp[,1]
y1 <- tmp[,2]
plot(yr1,y1,xlab='BJD',ylab='RV [km/s]',main='CHIRON')
lines(yr2,y2,col='red')
###residual
dy <- (y1-orb(yr1))*1e3
plot(yr1,dy,xlab='BJD',ylab='RV [m/s]')
lines(yr2,rv.out$gt,col='red')
lines(yr2,rv.out$st,col='blue')
out[['CHIRON']] <- cbind(yr1,y1,dy,data[,3])
                                        #    abline(v=c(1988.89,2001.73),col='blue')

data <- read.table('~/Documents/projects/dwarfs/data/aperture/alphaCen/alphaCenA_es.dat',header=TRUE)
tmp <- adjust.data(data[,1],data[,2],yr2,y2,factor=1e-3)
yr1 <- tmp[,1]
y1 <- tmp[,2]
plot(yr1,y1,xlab='BJD',ylab='RV [km/s]',main='ES')
lines(yr2,y2,col='red')
###residual
dy <- (y1-orb(yr1))*1e3
plot(yr1,dy,xlab='BJD',ylab='RV [m/s]')
lines(yr2,rv.out$gt,col='red')
lines(yr2,rv.out$st,col='blue')
out[['ES']] <- cbind(yr1,y1,dy,data[,3])

out[['UVES']] <- cbind(NA,NA,NA,NA)

data <- read.table('~/Documents/projects/dwarfs/data/other_sets/alpha_cen_old/HD128620_aat.dat',header=TRUE)
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
lines(yr2,rv.out$gt,col='red')
lines(yr2,rv.out$st,col='blue')
out[['AAT']] <- cbind(yr1,y1,dy,data[,3])

###
data <- read.table('~/Documents/projects/dwarfs/data/other_sets/alpha_cen_old/LCa.dat',header=TRUE)
tmp <- adjust.data(data[,1],data[,2],yr2,y2,factor=1e-3)
yr1 <- tmp[,1]
y1 <- tmp[,2]
plot(yr1,y1,xlab='BJD',ylab='RV [km/s]',main='LC')
lines(yr2,y2,col='red')
###residual
dy <- (y1-orb(yr1))*1e3
plot(yr1,dy,xlab='BJD',ylab='RV [m/s]')
lines(yr2,rv.out$gt,col='red')
lines(yr2,rv.out$st,col='blue')
out[['LC']] <- cbind(yr1,y1,dy,data[,3])

###
data <- read.table('~/Documents/projects/dwarfs/data/other_sets/alpha_cen_old/coraliea.dat',header=TRUE)
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
data <- read.table('~/Documents/projects/dwarfs/data/other_sets/alpha_cen_old/VLCa.dat',header=TRUE)
tmp <- adjust.data(data[,1],data[,2],yr2,y2,factor=1e-3)
yr1 <- tmp[,1]
y1 <- tmp[,2]
plot(yr1,y1,xlab='BJD',ylab='RV [km/s]',main='VLC')
lines(yr2,y2,col='red')
###residual
dy <- (y1-orb(yr1))*1e3
plot(yr1,dy,xlab='BJD',ylab='RV [m/s]')
lines(yr2,rv.out$gt,col='red')
lines(yr2,rv.out$st,col='blue')
out[['VLC']] <- cbind(yr1,y1,dy,data[,3])
dev.off()
fout <- gsub('.pdf','.Robj',fout)
cat(fout,'\n')
zgt <- rv.out$gt
zst <- rv.out$st
save(out,yr2,y2,zgt,zst,file=fout)
