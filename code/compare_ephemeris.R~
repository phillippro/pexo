library(orthopolynom)
library(pracma)
#options(digits= 16)
source('constants.R')
source('sofa_function.R')
source('astrometry_function.R')
source('general_function.R')
source('timing_function.R')
source('rv_function.R')

Tstart <- proc.time()
####Read parameter files
des <- c(405,414,421,430,435,436,438)
source('read_input.R')
for(DE in des){
    Par$DE <- DE
    source('load_data.R')
####
    OutBary <- time_Utc2tb(utc,Par)
    OutTime <- time_Ta2te(OutBary,Par)
    if(DE==405) P405 <- OutBary
    if(DE==430) P430 <- OutBary
    if(DE==436) P436 <- OutBary
    if(DE==438) P438 <- OutBary
    if(DE==414) P414 <- OutBary
    if(DE==421) P421 <- OutBary
    if(DE==435) P435 <- OutBary

    if(DE==405) timing405 <- OutTime
    if(DE==430) timing430 <- OutTime
    if(DE==436) timing436 <- OutTime
    if(DE==438) timing438 <- OutTime
    if(DE==414) timing414 <- OutTime
    if(DE==421) timing421 <- OutTime
    if(DE==435) timing435 <- OutTime
}
TtTdbMethod <- Par$TtTdbMethod
Tstep <- diff(rowSums(utc))[1]
fout <- paste0('../input/Tstep',Tstep,'_',TtTdbMethod,'_DE.Robj')
save(P405,P414,P421,P430,P435,P436,P438,timing405,timing414,timing421,timing430,timing435,timing436,timing438,file=fout)

truncate <- TRUE
dir <- '../input/'
#load(paste0(dir,TtTdbMethod,'_DE.Robj'))
y405 <- -time_T2mT2(timing430$BJDtdb,timing405$BJDtdb)*1e9
y414 <- -time_T2mT2(timing430$BJDtdb,timing414$BJDtdb)*1e9
y421 <- -time_T2mT2(timing430$BJDtdb,timing421$BJDtdb)*1e9
y435 <- -time_T2mT2(timing430$BJDtdb,timing435$BJDtdb)*1e9
y436 <- -time_T2mT2(timing430$BJDtdb,timing436$BJDtdb)*1e9
y438 <- -time_T2mT2(timing430$BJDtdb,timing438$BJDtdb)*1e9

x405 <- -time_T2mT2(P430$JDtdb,P405$JDtdb)*1e9
x414 <- -time_T2mT2(P430$JDtdb,P414$JDtdb)*1e9
x421 <- -time_T2mT2(P430$JDtdb,P421$JDtdb)*1e9
x435 <- -time_T2mT2(P430$JDtdb,P435$JDtdb)*1e9
x436 <- -time_T2mT2(P430$JDtdb,P436$JDtdb)*1e9
x438 <- -time_T2mT2(P430$JDtdb,P438$JDtdb)*1e9

rho405 <- sqrt(rowSums((P405$SG[,1:3]-P430$SG[,1:3])^2))*1e3
rho414 <- sqrt(rowSums((P414$SG[,1:3]-P430$SG[,1:3])^2))*1e3
rho421 <- sqrt(rowSums((P421$SG[,1:3]-P430$SG[,1:3])^2))*1e3
rho435 <- sqrt(rowSums((P435$SG[,1:3]-P430$SG[,1:3])^2))*1e3
rho436 <- sqrt(rowSums((P436$SG[,1:3]-P430$SG[,1:3])^2))*1e3
rho438 <- sqrt(rowSums((P438$SG[,1:3]-P430$SG[,1:3])^2))*1e3

v405 <- sqrt(rowSums((P405$SG[,4:6]-P430$SG[,4:6])^2))*1e6
v414 <- sqrt(rowSums((P414$SG[,4:6]-P430$SG[,4:6])^2))*1e6
v421 <- sqrt(rowSums((P421$SG[,4:6]-P430$SG[,4:6])^2))*1e6
v435 <- sqrt(rowSums((P435$SG[,4:6]-P430$SG[,4:6])^2))*1e6
v436 <- sqrt(rowSums((P436$SG[,4:6]-P430$SG[,4:6])^2))*1e6
v438 <- sqrt(rowSums((P438$SG[,4:6]-P430$SG[,4:6])^2))*1e6

truncate <- FALSE
fout <- paste0('../results/ephemeris_comparison_BJDtdb_tttdb',TtTdbMethod,'_',truncate,'.pdf')
cat(fout,'\n')
pdf(fout,6,6)
size <- 1.2
#tmin <- min(rowSums(utc))
#tt <- rowSums(utc)-tmin
tmin <- min(rowSums(P430$JDtdb))
tt <- rowSums(P430$JDtdb)-tmin
par(mar=c(5,5,1,1),cex=size,cex.axis=size,cex.lab=size)
ylim <- range(y405,y414,y421,y435,y436,y438)
if(truncate) ylim <- range(y421,y435,y436,y438)
if(!truncate){
plot(tt,y405,xlab=expression(JD[TDB]-2442000.5),ylab=expression(Delta*BJD[TDB]*" [ns]"),ylim=ylim,type='l')
lines(tt,y414,col='blue')
lines(tt,y421,col='red')
}else{
plot(tt,y421,xlab=expression(JD[TDB]-2442000.5),ylab=expression(Delta*BJD[TDB]*" [ns]"),ylim=ylim,type='l',col='red')
}
lines(tt,y435,col='green')
lines(tt,y436,col='orange')
lines(tt,y438,col='steelblue')
#legend('topright',inset=c(-0.4,0),xpd=NA,legend=c('DE405','DE414','DE421','DE435','DE436','DE438'),lty=1,col=c('black','blue','red','green','orange','steelblue'),bty='n')
dev.off()

fout <- paste0('../results/ephemeris_comparison_JDtdb_tttdb',TtTdbMethod,'_',truncate,'.pdf')
cat(fout,'\n\n')
pdf(fout,6,6)
size <- 1.2
par(mar=c(5,5,1,1),cex=size,cex.axis=size,cex.lab=size)
ylim <- range(x405,x414,x421,x435,x436,x438)
if(truncate) ylim <- range(x421,x435,x436,x438)
if(!truncate){
plot(tt,x405,xlab=expression(JD[TDB]-2442000.5),ylab=expression(Delta*BJD[TDB]*" [ns]"),ylim=ylim,type='l')
lines(tt,x414,col='blue')
lines(tt,x421,col='red')
}else{
plot(tt,x421,xlab=expression(JD[TDB]-2442000.5),ylab=expression(Delta*BJD[TDB]*" [ns]"),ylim=ylim,type='l',col='red')
}
lines(tt,x435,col='green')
lines(tt,x436,col='orange')
lines(tt,x438,col='steelblue')
#legend('topright',inset=c(-0.4,0),xpd=NA,legend=c('DE405','DE414','DE421','DE435','DE436','DE438'),lty=1,col=c('black','blue','red','green','orange','steelblue'),bty='n')
dev.off()

fout <- paste0('../results/ephemeris_comparison_pos_tttdb',TtTdbMethod,'_',truncate,'.pdf')
cat(fout,'\n\n')
pdf(fout,6,6)
size <- 1.2
par(mar=c(5,5,1,1),cex=size,cex.axis=size,cex.lab=size)
ylim <- range(rho405,rho414,rho421,rho435,rho436,rho438)
if(truncate) ylim <- range(rho421,rho435,rho436,rho438)
if(!truncate){
plot(tt,rho405,xlab=expression(JD[TDB]-2442000.5),ylab=expression(Delta*r[SG]*" [m]"),ylim=ylim,type='l')
lines(tt,rho414,col='blue')
lines(tt,rho421,col='red')
}else{
plot(tt,rho421,xlab=expression(JD[TDB]-2442000.5),ylab=expression(Delta*r[SG]*" [m]"),ylim=ylim,type='l',col='red')
}
lines(tt,rho435,col='green')
lines(tt,rho436,col='orange')
lines(tt,rho438,col='steelblue')
#legend('topright',inset=c(-0.4,0),xpd=NA,legend=c('DE405','DE414','DE421','DE435','DE436','DE438'),lty=1,col=c('black','blue','red','green','orange','steelblue'),bty='n')
dev.off()

fout <- paste0('../results/ephemeris_comparison_vel_tttdb',TtTdbMethod,'_',truncate,'.pdf')
cat(fout,'\n\n')
pdf(fout,6,6)
size <- 1.2
par(mar=c(5,5,1,1),cex=size,cex.axis=size,cex.lab=size)
ylim <- range(v405,v414,v421,v435,v436,v438)
if(truncate) ylim <- range(v421,v435,v436,v438)
if(!truncate){
plot(tt,v405,xlab=expression(JD[TDB]-2442000.5),ylab=expression(Delta*v[SG]*" [mm/s]"),ylim=ylim,type='l')
lines(tt,v414,col='blue')
lines(tt,v421,col='red')
}else{
plot(tt,v421,xlab=expression(JD[TDB]-2442000.5),ylab=expression(Delta*v[SG]*" [mm/s]"),ylim=ylim,type='l',col='red')
}
lines(tt,v435,col='green')
lines(tt,v436,col='orange')
lines(tt,v438,col='steelblue')
#
dev.off()

pdf('../results/legend.pdf',6,0.5)
par(mar=c(0,0,0,0),cex=0.8)
plot(0,xaxt='n',yaxt='n',cex=0,bty='n',xlab='',ylab='')
legend('center',xpd=NA,legend=c('DE405','DE414','DE421','DE435','DE436','DE438'),lty=1,col=c('black','blue','red','green','orange','steelblue'),bty='n',horiz=TRUE)
dev.off()
