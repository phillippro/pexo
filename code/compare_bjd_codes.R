options(digits=22)
library(ggplot2)
library(magicaxis)
library(latex2exp)
DJM0 <- 2400000.5
DAYSEC <- 86400.0
tb <- read.table('../data/10700tb.dat',header=TRUE)
E10 <- read.table('../data/10700E10.dat')
#T2 <- read.table('10700T2.dat',header=TRUE)
tmp <- readLines('../data/10700T2.dat')#deal with high precision JDtt data
T2 <- c()
for(j in 1:length(tmp)){
    s1 <- unlist(strsplit(tmp[j],split=' '))
    if(j==1){
        header <- c('BJDtdb0','BJDtdb1','JDtt0','JDtt1',s1[-(1:2)])
    }else{
        s2 <- unlist(strsplit(s1[2],split='\\.'))
        s2[2] <- paste0('0.',s2[2])
        s3 <- unlist(strsplit(s1[1],split='\\.'))
        s3[2] <- paste0('0.',s3[2])
        T2 <- rbind(T2,as.numeric(c(s3,s2,s1[-(1:2)])))
    }
}
colnames(T2) <- header
fname <- '../results/comparing_BJD_codes.pdf'
cat('output pdf:\n',fname,'\n')
pdf(fname,8,8)
par(mfrow=c(2,2))
jd.ref <- sum(tb[1,c('JDutc0','JDutc1')])
Nt <- nrow(tb)
utc.ref <- t(replicate(Nt,as.numeric(tb[1,c('JDutc0','JDutc1')])))
utc.rel <- rowSums(tb[,c('JDutc0','JDutc1')]-utc.ref)
utc.abs <- rowSums(tb[,c('JDutc0','JDutc1')])
bjd.tb <- (tb[,'BJDtdb0']-utc.ref[,1])+tb[,'BJDtdb1']-utc.ref[,2]
bjd.E10 <- E10[,1]-jd.ref
bjd.T2 <- (T2[,'BJDtdb0']+DJM0-utc.ref[,1])+T2[,'BJDtdb1']-utc.ref[,2]
dt.pe <- ((tb[,'BJDtdb1']-T2[,'BJDtdb1'])+(tb[,'BJDtdb0']-DJM0-T2[,'BJDtdb0']))*DAYSEC

##BJD-JD
plot(utc.rel,bjd.tb-utc.rel,xlab='JD[UTC]-JD[UTC0]',ylab='BJD[TDB]-JD[UTC0]',type='l')
lines(utc.rel,bjd.E10-utc.rel,col='red')
lines(utc.rel,bjd.T2-utc.rel,col='blue')

xlim=c(ceiling(min(utc.abs)),2456000)
ind <- which(utc.abs>=xlim[1] & utc.abs<=xlim[2])
##BJD[xx]-BJD[T2]
#y1 <- abs(bjd.tb-bjd.T2)*DAYSEC
#y1[y1==0] <- 1e-10
y1 <- abs(dt.pe)
y2 <- abs(bjd.E10-bjd.T2)*DAYSEC
plot(utc.abs[ind],y1[ind],xlab='JD[UTC]',ylab=expression('|'*Delta*'BJD[TDB]|'),type='l',log='y',ylim=range(y1,2*y2),xlim=c(ceiling(min(utc.abs)),2456000))
#magaxis(side=2)
lines(utc.abs[ind],y2[ind],col='blue')
legend('top',inset=c(0,-0.2),xpd=NA,legend=c('PEXO     ','E10'),col=c('black','blue'),horiz=TRUE,lty=1,bty='n')

##BJD[xx]-BJD[T2]
y1 <- abs(bjd.tb-bjd.T2)*DAYSEC
y1[y1==0] <- 1e-10
y2 <- abs(bjd.E10-bjd.T2)*DAYSEC
plot(utc.abs[ind],y1[ind],xlab='JD[UTC]',ylab='BJD w.r.t. TEMPO2 values',type='l',log='y',ylim=range(y1,y2),xlim=c(ceiling(min(utc.abs)),2456000))
lines(utc.abs,y2,col='blue')
legend('top',inset=c(0,-0.2),xpd=NA,legend=c('PEXO-TEMPO2','PEXO-TEMPO2'),col=c('black','blue'),horiz=TRUE,lty=1,bty='n')

##TT-UTC;
dtt.utc <- rowSums(tb[,c('JDtt0','JDtt1')]-tb[,c('JDutc0','JDutc1')])*DAYSEC-T2[,'tt.utc']
plot(utc.abs[ind],dtt.utc[ind],xlab='JD[UTC]',ylab=expression(Delta*'TT [s]'),type='l',xlim=c(ceiling(min(utc.abs)),2456000))

##TDB-TT;
dtdb.tt <- (tb[,'JDtdb1']-tb[,'JDtt1'])*DAYSEC-T2[,'tdb.tt']
plot(utc.abs[ind],dtdb.tt[ind],xlab='JD[UTC]',ylab='relative TDB-TT [s]',type='l',xlim=c(ceiling(min(utc.abs)),2456000))

##TT[tb]-TT[T2]
dtt <- (tb[,'JDtt1']-tb[,'JDutc1'])/DAYSEC-T2[,'tt.utc']
plot(utc.abs[ind],dtt[ind],xlab='JD[UTC]',ylab=expression(Delta*'TT [s]'),type='l',xlim=c(ceiling(min(utc.abs)),2456000))

##Roemer delay
Dt.roemer <- tb[,'roemer']-T2[,'roemer']
#plot(utc.abs[ind],Dt.roemer[ind],xlab='JD[UTC]',ylab=expression(Delta[R'\u0298']*' [s]'),type='l',xlim=c(ceiling(min(utc.abs)),2456000))
#plot(utc.abs[ind],Dt.roemer[ind],xlab='JD[UTC]',ylab=labs(y=expression(M['☉'])),type='l')
plot(utc.abs[ind],Dt.roemer[ind],xlab='JD[UTC]',ylab='Roemer delay difference [s]',type='l')

##Shapiro delay
Dt.shapiro <- tb[,'shapiro']-T2[,'shapiro']
plot(utc.abs[ind],Dt.shapiro[ind],xlab='JD[UTC]',ylab=expression('Shapiro delay difference [s]'),type='l')

##BJD[tb]-BJD[T2]
plot(utc.abs[ind],dt.pe[ind],xlab='JD[UTC]',ylab='BJD[tb]-BJD[T2] [s]',type='l')

##Roemer-Shapiro delay
dt.rs <- (tb[,'roemer']-tb[,'shapiro'])-(T2[,'roemer']-T2[,'shapiro'])
plot(utc.abs[ind],dt.rs[ind],xlab='JD[UTC]',ylab='Roemer-Shapiro w.r.t. T2 [s]',type='l',xlim=c(ceiling(min(utc.abs)),2456000))

##tropo
plot(utc.abs[ind],T2[ind,'tropo'],xlab='JD[UTC]',ylab='tropo delay [s]',type='l',xlim=c(ceiling(min(utc.abs)),2456000))

##BJD[tb]-(BJD[T2]+tropo)
dt.pet <- ((tb[,'BJDtdb1']-T2[,'BJDtdb1'])+(tb[,'BJDtdb0']-DJM0-T2[,'BJDtdb0']))*DAYSEC-T2[,'tropo']
plot(utc.abs[ind],dt.pet[ind],xlab='JD[UTC]',ylab='BJD[tb]-BJD[T2]-tropo [s]',type='l',xlim=c(ceiling(min(utc.abs)),2456000))

##dt.pet-dt.rs
plot(utc.abs[ind],dt.pet[ind]-dt.rs[ind],xlab='JD[UTC]',ylab='dt.pet-dt.rs',type='l',xlim=c(ceiling(min(utc.abs)),2456000))

plot(utc.abs[ind],T2[ind,'shapiroJ'],xlab='JD[UTC]',ylab='shapiroJ',type='l',xlim=c(ceiling(min(utc.abs)),2456000))

plot(utc.abs[ind],T2[ind,'shapiroS'],xlab='JD[UTC]',ylab='shapiroS',type='l',xlim=c(ceiling(min(utc.abs)),2456000))

plot(utc.abs[ind],T2[ind,'shapiroV'],xlab='JD[UTC]',ylab='shapiroV',type='l',xlim=c(ceiling(min(utc.abs)),2456000))

plot(utc.abs[ind],T2[ind,'shapiroU'],xlab='JD[UTC]',ylab='shapiroU',type='l',xlim=c(ceiling(min(utc.abs)),2456000))

plot(utc.abs[ind],T2[ind,'shapiroN'],xlab='JD[UTC]',ylab='shapiroN',type='l',xlim=c(ceiling(min(utc.abs)),2456000))

plot(utc.abs[ind],rowSums(T2[ind,c('shapiroJ','shapiroS','shapiroV','shapiroU','shapiroN')]),xlab='JD[UTC]',ylab='shapiroJSVUN',type='l',xlim=c(ceiling(min(utc.abs)),2456000))

dev.off()

