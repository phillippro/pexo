dir <- '../results/'
fname <- paste0(dir,star,'_astro_rv.pdf')
cat('output pdf:\n',fname,'\n')
pdf(fname,16,16)
par(mfrow=c(4,4))

zRT2 <- rowSums(uSB*SB[,4:6,drop=FALSE])/CKMPS#barycentric corrected redshift of the Target system barycenter (TSB)

###astrometry
uSB <- SB[,1:3]/sqrt(rowSums(SB[,1:3]^2))
dec.SB <- asin(uSB[,3])#rad
ra.SB <- xy2phi(uSB[,1],uSB[,2])
dra.plx <- (dr$OT[,2]-ra.SB)*pc2au#arcsec
ddec.plx <- (dr$OT[,1]-dec.SB)*pc2au#arcsec
dra.all <- (dr$all[,2]-dr$all[1,2])*pc2au
ddec.all <- (dr$all[,1]-dr$all[1,1])*pc2au
dra.OT <- (dr$OT[,2]-dr$OT[1,2])*pc2au
ddec.OT <- (dr$OT[,1]-dr$OT[1,1])*pc2au
dra.SB <- (ra.SB-ra.SB[1])*pc2au
ddec.SB <- (dec.SB-dec.SB[1])*pc2au
##OT
plot(dr$all[,2],dr$all[,1],xlab='ra [rad]',ylab='dec [rad]', main='absolute astrometry OT (with all effects)')
plot(dra.all,ddec.all,xlab='dra [as]',ylab='ddec [as]', main='relative astrometry OT (with all effects)')
plot(dra.OT,ddec.OT,xlab='dra [as]',ylab='ddec [as]', main='relative astrometry OT (geometry)')
##SB
plot(ra.SB,dec.SB,xlab='ra [rad]',ylab='dec [rad]', main='absolute astrometry SB')
plot(dra.SB,ddec.SB,xlab='dra [as]',ylab='ddec [as]', main='relative astrometry SB')
##parallax
plot(dra.plx,ddec.plx,xlab='dra [as]',ylab='ddec [as]', main='parallax')

##aberration
plot(dr$abe[,2],dr$abe[,1],xlab='dra [as]',ylab='ddec [as]',main='stellar aberration')
##lensing
plot(dr$len[,2],dr$len[,1],xlab='dra [as]',ylab='ddec [as]',main='lensing')

###RV
##absolute
plot(rowSums(bjd.tcb),RV,xlab='bjd',ylab='RV [m/s]',main='RV all')
plot(rowSums(bjd.tcb),RV+rv$ro,xlab='bjd',ylab='RV [m/s]',main='RV all+rvRO')
plot(rowSums(bjd.tcb),RV+rv$ro-rv$rt,xlab='bjd',ylab='RV [m/s]',main='RV all+rvRO-rvRT')
plot(rowSums(bjd.tcb),rv$bary,xlab='bjd',ylab='RV [m/s]',main='RV bary')
plot(rowSums(bjd.tcb),rv$obs,xlab='bjd',ylab='RV [m/s]',main='RV obs')
plot(rowSums(bjd.tcb),-rv$obs+rv$gvo-rv$ro-rv$so,xlab='bjd',ylab='RV [m/s]',main='-OBS+GVO-RO-SO')
plot(rowSums(bjd.tcb),RV+rv$ro+rv$gvo+rv$so,xlab='bjd',ylab='RV [m/s]',main='RV+GO+VO+RO+SO')
plot(rowSums(bjd.tcb),RV+rv$ro-rv$gt,xlab='bjd',ylab='RV [m/s]',main='RV all+rvRO-rvGT')
plot(rowSums(bjd.tcb),RV+rv$ro-rv$gt,xlab='bjd',ylab='RV [m/s]',main='RV all+rvRO-rvGT+rvVO')
plot(rowSums(bjd.tcb),RV-rv$rt,xlab='bjd',ylab='RV [m/s]',main='RV all-rvRT')
plot(rowSums(bjd.tcb),rv$gvo,xlab='bjd',ylab='RV [m/s]',main='rvGVO')
plot(rowSums(bjd.tcb),rv$gt,xlab='bjd',ylab='RV [m/s]',main='rvGT')
plot(rowSums(bjd.tcb),rv$vt,xlab='bjd',ylab='RV [m/s]',main='rvVT')
plot(rowSums(bjd.tcb),rv$ro,xlab='bjd',ylab='RV [m/s]',main='rvRO')
plot(rowSums(bjd.tcb),rv$rt,xlab='bjd',ylab='RV [m/s]',main='rvRT')
plot(rowSums(bjd.tcb),zRT2*CMPS,xlab='bjd',ylab='RV [m/s]',main='rvRT2')
plot(rowSums(bjd.tcb),rv$rt-zRT2*CMPS,xlab='bjd',ylab='RV [m/s]',main='rvRT-rvRT2')
plot(rowSums(bjd.tcb),-rv$so,xlab='bjd',ylab='RV [m/s]',main='-rvSO',type='l')
plot(rowSums(bjd.tcb),rv$st,xlab='bjd',ylab='RV [m/s]',main='rvST')

#plot(rowSums(bjd.tcb),dRV,xlab='bjd',ylab='dRV')
#plot(rowSums(bjd.tcb),RV-rv-rowSums(SO[,1:3]*uOT),xlab='bjd',ylab='RV-rv.star-rv.obs')
dev.off()
