options(digits=20)
source('../../code/timing_function.R')
source('../../code/general_function.R')
source('../../code/constants.R')
saveData <- FALSE
aca <- read.table('astrometry-acena.txt',header=TRUE)
acb <- read.table('astrometry-acenb.txt',header=TRUE)
hipa <- c(ra=219.92041034,dec=-60.83514707,plx=742.12,pmra=-3678.19,pmdec=481.84)
hipb <- c(ra=219.91412833,dec=-60.83947139,plx=742.12,pmra=-3600.35,pmdec=952.11)#original hipparcos solution
yr <- aca[,'t']
jd <- time_Yr2jd(yr)
deca <- hipa['dec']+aca[,'obsY']/3.6e6
#raa <- hipa['ra']+aca[,'obsX']/3.6e6#assume that dalpha is alpha1-alpha2
raa <- hipa['ra']+aca[,'obsX']/3.6e6/cos(deca/180*pi)#assume that dalpha* is (alpha1-alpha2)*cos(delta)
eraa <- abs(aca[,'barStartX']-aca[,'barEndX'])/2#mas
edeca <- abs(aca[,'barStartY']-aca[,'barEndY'])/2#mas

decb <- hipb['dec']+acb[,'obsY']/3.6e6
#rab <- hipb['ra']+acb[,'obsX']/3.6e6
rab <- hipb['ra']+acb[,'obsX']/3.6e6/cos(decb/180*pi)#assume that dalpha* is (alpha1-alpha2)*cos(delta)
erab <- abs(acb[,'barStartX']-acb[,'barEndX'])/2#mas
edecb <- abs(acb[,'barStartY']-acb[,'barEndY'])/2#mas

out1 <- cbind(rowSums(jd),raa,eraa,deca,edeca)
out2 <- cbind(rowSums(jd),rab,erab,decb,edecb)
dra  <- (rab-raa)*cos((decb+deca)/2*pi/180)*3.6e6#mas
ddec <- (decb-deca)*3.6e6
#dra <- acb[,'obsX']-aca[,'obsX']
#ddec <- acb[,'obsY']-aca[,'obsY']
edra <- sqrt(eraa^2+erab^2)
eddec <- sqrt(edeca^2+edecb^2)
out3 <- cbind(rowSums(jd),dra,edra,ddec,eddec)

if(saveData){
write.table(cbind(out1,0.55),file='../HD128620/HD128620_hip.abs',quote=FALSE,row.names=FALSE,col.names=c('JD', 'ra', 'era', 'dec', 'edec','lambda'))
write.table(cbind(out1,0.55),file='../HD128620v1/HD128620_hip.abs',quote=FALSE,row.names=FALSE,col.names=c('JD', 'ra', 'era', 'dec', 'edec','lambda'))

write.table(cbind(out2,0.55),file='HD128621_hip.abs',quote=FALSE,row.names=FALSE,col.names=c('JD', 'ra', 'era', 'dec', 'edec','lambda'))
write.table(cbind(out2,0.55),file='../HD128621/HD128621_hip.abs',quote=FALSE,row.names=FALSE,col.names=c('JD', 'ra', 'era', 'dec', 'edec','lambda'))
write.table(cbind(out2,0.55),file='../HD128621v1/HD128621_hip.abs',quote=FALSE,row.names=FALSE,col.names=c('JD', 'ra', 'era', 'dec', 'edec','lambda'))
write.table(cbind(out2,0.55),file='../HD128621v2/HD128621_hip.abs',quote=FALSE,row.names=FALSE,col.names=c('JD', 'ra', 'era', 'dec', 'edec','lambda'))
write.table(cbind(out2,0.55),file='HD128621_hip.abs',quote=FALSE,row.names=FALSE,col.names=c('JD', 'ra', 'era', 'dec', 'edec','lambda'))

write.table(cbind(out3,0.55),file='../HD128621/HD128621.rel',quote=FALSE,row.names=FALSE,col.names=c('JD', 'dRA', 'edRA', 'dDE', 'edDE','lambda'))
write.table(cbind(out3,0.55),file='HD128621.rel',quote=FALSE,row.names=FALSE,col.names=c('JD', 'dRA', 'edRA', 'dDE', 'edDE','lambda'))
write.table(cbind(out3,0.55),file='../HD128621v1/HD128621_hip.rel',quote=FALSE,row.names=FALSE,col.names=c('JD', 'dRA', 'edRA', 'dDE', 'edDE','lambda'))
write.table(cbind(out3,0.55),file='../HD128621v2/HD128621_hip.rel',quote=FALSE,row.names=FALSE,col.names=c('JD', 'dRA', 'edRA', 'dDE', 'edDE','lambda'))
}

alma <- read.table('HD128621_ALMA.rel',header=TRUE)
wds <- read.table('HD128621_WDS.rel',header=TRUE)
k17 <- read.table('HD128621_K17.rel',header=TRUE)
pdf('test_hip.pdf',16,16)
par(mfrow=c(4,4))

plot(rowSums(jd),out3[,2],xlab='JD',ylab='dRA[mas]',main='hip')
plot(rowSums(jd),out3[,4],xlab='JD',ylab='dDE[mas]',main='hip')
plot(out3[,2],out3[,4],xlab='dra[mas]',ylab='ddec[mas]',main='hip')

plot(out1[,1],out1[,2],xlab='JD',ylab='dRA[mas]',main='absolute astrometry of alpha Cen A')
plot(out1[,1],out1[,4],xlab='JD',ylab='dDE[mas]',main='absolute astrometry of alpha Cen A')
plot(out1[,2],out1[,4],xlab='dra[mas]',ylab='ddec[mas]',main='absolute astrometry of alpha Cen A')

plot(out2[,1],out2[,2],xlab='JD',ylab='dRA[mas]',main='absolute astrometry of alpha Cen B')
plot(out2[,1],out2[,4],xlab='JD',ylab='dDE[mas]',main='absolute astrometry of alpha Cen B')
plot(out2[,2],out2[,4],xlab='dra[mas]',ylab='ddec[mas]',main='absolute astrometry of alpha Cen B')

plot(wds[,1],wds[,2],xlab='JD',ylab='dRA[mas]',main='wds')
plot(wds[,1],wds[,4],xlab='JD',ylab='dDE[mas]',main='wds')
plot(wds[,2],wds[,4],xlab='dra[mas]',ylab='ddec[mas]',main='wds')
points(wds[437,2],wds[437,4],col='yellow',pch=20)

plot(k17[,1],k17[,2],xlab='JD',ylab='dRA[mas]',main='K17')
plot(k17[,1],k17[,4],xlab='JD',ylab='dDE[mas]',main='K17')
plot(k17[,2],k17[,4],xlab='dra[mas]',ylab='ddec[mas]',main='K17')

if(exists('alma')){
   jd <- alma[,1]
   plot(jd,alma[,2],xlab='JD',ylab='dRA[mas]',main='alma')
   plot(jd,alma[,4],xlab='JD',ylab='dDE[mas]',main='alma')
   plot(alma[,2],alma[,4],xlab='dra[mas]',ylab='ddec[mas]',main='alma')
   tot <- rbind(as.matrix(out3[,1:5]),as.matrix(alma[,1:5]),as.matrix(wds[,1:5]),as.matrix(k17[,1:5]))
}else{
   tot <- rbind(as.matrix(out3[,1:5]),as.matrix(wds[,1:5]),as.matrix(k17[,1:5]))
}

plot(tot[,1],tot[,2],xlab='JD',ylab='dRA[mas]',main='tot')
plot(tot[,1],tot[,4],xlab='JD',ylab='dDE[mas]',main='tot')
plot(tot[,2],tot[,4],xlab='dra[mas]',ylab='ddec[mas]',main='tot',col='white')
points(out3[,2],out3[,4],col='red')
points(wds[,2],wds[,4],col='black')
points(k17[,2],k17[,4],col='green')
if(exists('alma')) points(alma[,2],alma[,4],col='blue')
points(0,0,pch='+')
points(wds[437,2],wds[437,4],col='yellow',pch=20)
dev.off()
