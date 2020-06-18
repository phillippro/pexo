library(optparse)
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
#des <- c(405,414,421,430,435,436,438)
des <- c(405,430)
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
#save(P405,P414,P421,P430,P435,P436,P438,timing405,timing414,timing421,timing430,timing435,timing436,timing438,file=fout)

truncate <- TRUE
dir <- '../input/'
#load(paste0(dir,TtTdbMethod,'_DE.Robj'))
ylist <- xlist <- rho.list <- vlist <- list()
colors <- c()
if(exists('timing405')){
    y405 <- -time_T2mT2(timing430$BJDtdb,timing405$BJDtdb)*1e9
    x405 <- -time_T2mT2(P430$JDtdb,P405$JDtdb)*1e9
    rho405 <- sqrt(rowSums((P405$SG[,1:3]-P430$SG[,1:3])^2))*1e3
    v405 <- sqrt(rowSums((P405$SG[,4:6]-P430$SG[,4:6])^2))*1e6
    ylist$DE405 <- y405
    xlist$DE405 <- x405
    vlist$DE405 <- v405
    rho.list$DE405 <- rho405
    colors <- c(colors,'black')
}
if(exists('timing414')){
    y414 <- -time_T2mT2(timing430$BJDtdb,timing414$BJDtdb)*1e9
    x414 <- -time_T2mT2(P430$JDtdb,P414$JDtdb)*1e9
    rho414 <- sqrt(rowSums((P414$SG[,1:3]-P430$SG[,1:3])^2))*1e3
    v414 <- sqrt(rowSums((P414$SG[,4:6]-P430$SG[,4:6])^2))*1e6
    ylist$DE414 <- y414
    xlist$DE414 <- x414
    vlist$DE414 <- v414
    rho.list$DE414 <- rho414
    colors <- c(colors,'blue')
}
if(exists('timing421')){
    y421 <- -time_T2mT2(timing430$BJDtdb,timing421$BJDtdb)*1e9
    x421 <- -time_T2mT2(P430$JDtdb,P421$JDtdb)*1e9
    rho421 <- sqrt(rowSums((P421$SG[,1:3]-P430$SG[,1:3])^2))*1e3
    v421 <- sqrt(rowSums((P421$SG[,4:6]-P430$SG[,4:6])^2))*1e6
    ylist$DE421 <- y421
    xlist$DE421 <- x421
    vlist$DE421 <- v421
    rho.list$DE421 <- rho421
    colors <- c(colors,'red')
}
if(exists('timing435')){
    y435 <- -time_T2mT2(timing430$BJDtdb,timing435$BJDtdb)*1e9
    x435 <- -time_T2mT2(P430$JDtdb,P435$JDtdb)*1e9
    rho435 <- sqrt(rowSums((P435$SG[,1:3]-P430$SG[,1:3])^2))*1e3
    v435 <- sqrt(rowSums((P435$SG[,4:6]-P430$SG[,4:6])^2))*1e6
    ylist$DE435 <- y435
    xlist$DE435 <- x435
    vlist$DE435 <- v435
    rho.list$DE435 <- rho435
    colors <- c(colors,'green')
}
if(exists('timing436')){
    y436 <- -time_T2mT2(timing430$BJDtdb,timing436$BJDtdb)*1e9
    x436 <- -time_T2mT2(P430$JDtdb,P436$JDtdb)*1e9
    rho436 <- sqrt(rowSums((P436$SG[,1:3]-P430$SG[,1:3])^2))*1e3
    v436 <- sqrt(rowSums((P436$SG[,4:6]-P430$SG[,4:6])^2))*1e6
    ylist$DE436 <- y436
    xlist$DE436 <- x436
    vlist$DE436 <- v436
    rho.list$DE436 <- rho436
    colors <- c(colors,'orange')
}

if(exists('timing438')){
    y438 <- -time_T2mT2(timing430$BJDtdb,timing438$BJDtdb)*1e9
    x438 <- -time_T2mT2(P430$JDtdb,P438$JDtdb)*1e9
    rho438 <- sqrt(rowSums((P438$SG[,1:3]-P430$SG[,1:3])^2))*1e3
    v438 <- sqrt(rowSums((P438$SG[,4:6]-P430$SG[,4:6])^2))*1e6
    ylist$DE438 <- y438
    xlist$DE438 <- x438
    vlist$DE438 <- v438
    rho.list$DE438 <- rho438
    colors <- c(colors,'steelblue')
}

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
ylim <- range(unlist(ylist))
plot(tt,unlist(ylist),xlab=expression(JD[TDB]-2442000.5),ylab=expression(Delta*BJD[TDB]*" [ns]"),ylim=ylim,type='l')
for(j in 1:length(ylist)){
    lines(tt,ylist[[j]],col=colors[j])
}
#legend('topright',inset=c(-0.4,0),xpd=NA,legend=names(ylist),lty=1,col=colors,bty='n')
dev.off()

fout <- paste0('../results/ephemeris_comparison_JDtdb_tttdb',TtTdbMethod,'_',truncate,'.pdf')
cat(fout,'\n\n')
pdf(fout,6,6)
size <- 1.2
par(mar=c(5,5,1,1),cex=size,cex.axis=size,cex.lab=size)
ylim <- range(unlist(xlist))
plot(tt,unlist(xlist),xlab=expression(JD[TDB]-2442000.5),ylab=expression(Delta*BJD[TDB]*" [ns]"),ylim=ylim,type='l',col='white')
for(j in 1:length(xlist)){
    lines(tt,xlist[[j]],col=colors[j])
}
#legend('topright',inset=c(-0.4,0),xpd=NA,legend=names(xlist),lty=1,col=colors,bty='n')
dev.off()

fout <- paste0('../results/ephemeris_comparison_pos_tttdb',TtTdbMethod,'_',truncate,'.pdf')
cat(fout,'\n\n')
pdf(fout,6,6)
size <- 1.2
par(mar=c(5,5,1,1),cex=size,cex.axis=size,cex.lab=size)
ylim <- range(unlist(rho.list))
plot(tt,unlist(rho.list),xlab=expression(JD[TDB]-2442000.5),ylab=expression(Delta*r[SG]*" [m]"),ylim=ylim,type='l',col='white')
for(j in 1:length(rho.list)){
    lines(tt,rho.list[[j]],col=colors[j])
}
#legend('topright',inset=c(-0.4,0),xpd=NA,legend=names(rho.list),lty=1,col=colors,bty='n')
dev.off()

fout <- paste0('../results/ephemeris_comparison_vel_tttdb',TtTdbMethod,'_',truncate,'.pdf')
cat(fout,'\n\n')
pdf(fout,6,6)
size <- 1.2
par(mar=c(5,5,1,1),cex=size,cex.axis=size,cex.lab=size)
ylim <- range(unlist(vlist))
plot(tt,unlist(vlist),xlab=expression(JD[TDB]-2442000.5),ylab=expression(Delta*v[SG]*" [mm/s]"),ylim=ylim,type='l',col='white')
for(j in 1:length(vlist)){
    lines(tt,vlist[[j]],col=colors[j])
}
#legend('topright',inset=c(-0.4,0),xpd=NA,legend=names(vlist),lty=1,col=colors,bty='n')
dev.off()

pdf('../results/legend.pdf',6,0.5)
par(mar=c(0,0,0,0),cex=0.8)
plot(0,xaxt='n',yaxt='n',cex=0,bty='n',xlab='',ylab='')
legend('center',xpd=NA,legend=names(ylist),lty=1,col=colors,bty='n',horiz=TRUE)
dev.off()
