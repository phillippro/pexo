####check times
fpdf <- '../results/clock_correction.pdf'
cat('output pdf:\n',fpdf,'\n')
JDtai <- time_Jd2cal(tb$tai)
JDtai <- JDtai[,1]+JDtai[,2]/12+(JDtai[,3]+JDtai[,4])/365.25
pdf(fpdf,6,4)
par(mar=c(5,5,1,8))
terms <- c('UTC-TAI','TT-TAI','TDB-TAI','TCB-TAI','TCG-TAI')
cols <- c('black','red','blue','green','orange')
plot(JDtai,(tb$utc[,2]-tb$tai[,2])*DAYSEC,xlab='JD[TAI]',ylab='dt [s]',type='l',ylim=c(-60,60),col=cols[1])
lines(JDtai,(tb$tt[,2]-tb$tai[,2])*DAYSEC,col=cols[2])
lines(JDtai,(tb$tdb[,2]-tb$tai[,2])*DAYSEC,col=cols[3])
lines(JDtai,(tb$tcb[,2]-tb$tai[,2])*DAYSEC,col=cols[4])
lines(JDtai,(tb$tcg[,2]-tb$tai[,2])*DAYSEC,col=cols[5])
abline(h=0,lty=2)
legend('topright',inset=c(-0.4,0),xpd=NA,legend=terms,col=cols,lty=rep(1,5),bty='n')
dev.off()
