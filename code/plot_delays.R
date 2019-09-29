fpdf <- '../results/compare_delays.pdf'
cat('output:\n',fpdf,'\n')
pdf(fpdf,16,16)
size <- 1.2
par(mfrow=c(4,4),cex.axis=size,cex.lab=size,cex=size)
jd.utc <- rowSums(tb$utc)
jd.tt <- rowSums(tb$tt)
##TT-UTC
dTT <- rowSums(tb$tt-tb$utc)*DAYSEC
plot(jd.utc,dTT,xlab='JD[UTC]',ylab='TT-UTC')

##Einstein delay in the solar system; TCB-TT
dES <- rowSums(tb$tcb-tb$tt)*DAYSEC
plot(jd.tt,dES,xlab='JD[UTC]',ylab='TCB-TT')

##Roemer delay in the solar system;
plot(jd.tt,outBJD$roemer,xlab='JD[UTC]',ylab='Solar Roemer delay [s]')

##Solar Shapiro delay
dt.ss <- outBJD$shapiro
plot(jd.tt,dt.ss,xlab='JD[UTC]',ylab='Solar Shapiro delay [s]')

##Einstein delay due to relative motion of SSB and TSB
plot(jd.tt,dt.eis,xlab='JD[UTC]',ylab='Interstellar Einstein delay [s]')

##Vacuum propagation delay
plot(jd.tt,dt.vis,xlab='JD[UTC]',ylab='Vacuum delay [s]')

##geometric delay
plot(jd.tt,dt.geo,xlab='JD[UTC]',ylab='Geometric delay [s]')

##tau.o-tau.e
Dt <- dES-dt.ss-dt.eis-dt.geo
plot(jd.tt,Dt,xlab='JD[UTC]',ylab=expression(tau[e]-tau[o]*' [s]'))

dev.off()
