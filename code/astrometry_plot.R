options(scipen=-1)
###catalog astrometry
EleLowLimit <- 5#deg
#if(grepl('S2',Par$star) & grepl('ref',Par$RefType)){
if(FALSE){
    index <- which(OutTime$elevation*180/pi>EleLowLimit)#only show epochs with elevation angle larger than 10 degree
    if(length(index)==0){
        index <- which(OutTime$elevation*180/pi>EleLowLimit/2)
    }
    if(length(index)==0){
        index <- which(OutTime$elevation*180/pi>EleLowLimit/4)
    }
    if(length(index)==0){
        index <- 1:Par$Nepoch
    }
}else{
    index <- 1:Par$Nepoch
}
rd.cat <- t(c(Par$ra,Par$dec)*180/pi)
utc1   <- jd.utc <- rowSums(utc)
uOT <- OutTime$uOT
uST <- OutTime$uST
uOB <- OutTime$uOB
uOC <-  OutTime$uOC
rSC <-  OutTime$rSC
rd.OT <- gen_Xyz2lb(uOT)*180/pi#deg
rd.OTobs <- OutAstro$DirObs*180/pi#deg
rd.ST <- gen_Xyz2lb(uST)*180/pi#deg
##astrometry difference relativisitc effects
dra.OT <- (rd.OTobs[,1]-rd.OT[,1])*3600
ddec.OT <- (rd.OTobs[,2]-rd.OT[,2])*3600
dra.ST2 <- (rd.ST[,1]-rd.cat[,1])*3600
ddec.ST2 <- (rd.ST[,2]-rd.cat[,2])*3600
dir.OT <- gen_Xyz2lb(uOT)
dir.OC <- gen_Xyz2lb(uOC)
dir.OB <- gen_Xyz2lb(uOB)

xlab <- expression(Delta*alpha*'* [as]')
ylab <- expression(Delta*delta*' [as]')
dir.out <- '../results/'
if(!file.exists(dir.out)) system(paste('mkdir',dir.out))
dt <- diff(rowSums(utc))[1]
Ntime <- nrow(utc)
fname <- paste0(Par$star,'_astrometry_',Par$BinaryModel,'_Ntime',Ntime,'_',Par$RefType,'.pdf')
fout <- paste0(dir.out,'absolute_',fname)
cat(fout,'\n')
pdf(fout,12,12)
size <- 1.2
par(mfrow=c(3,3),mar=c(5,5,2,1),cex.lab=size)
##in order of decreasing effects

###astrometry change due to the heliocentric motion of the barycenter
lbSB <- gen_Xyz2lb(OutTime$uSB)
#plot((lbSB[,1]-lbSB[1,1])*pc2au,(lbSB[,2]-lbSB[1,2])*pc2au,xlab=xlab,ylab=ylab,main=expression(alpha*' Centauri barycenter motion'),pch=20,cex=0.2)

###refraction effects
plot(OutAstro$OffRef[index,1],OutAstro$OffRef[index,2],xlab=xlab,ylab=ylab,main=expression('P1: Atmospheric refraction'),pch=20,cex=0.2)

##aberration
plot(-OutAstro$OffAbe[index,1],OutAstro$OffAbe[index,2],xlab=xlab,ylab=ylab,main=expression('P2: Stellar aberration'),pch=20,cex=0.2)

##Binary motion
####prepare markers
yy <- time_Jd2yr(utc)
##    ts <- seq(1980,2050,by=10)
ts <- seq(1980,2040,by=5)
if(grepl('HD209100',Par$star))    ts <- seq(1990,2030,by=2)
ts <- ts[ts<max(yy)]
inds <- sapply(ts, function(t) which.min(abs(yy-t)))

if(Par$Np>0){
    delta <- gen_Xyz2lb(OutTime$uOT)[,2]
    lbBT <- gen_CalOffset(OutTime$uST,OutTime$uSB,bref=delta)
    lbBC <- gen_CalOffset(gen_CalUnit(rSC),OutTime$uSB,bref=delta)
    lbTC <- gen_CalOffset(gen_CalUnit(rSC),OutTime$uST,bref=delta)
####plot binary orbit to be comparable with Fig. 1 in Pourbaix et al. 1999
                                        #x1 <- (-x*(m1+m2)/m2)
                                        #y1 <- (-y*(m1+m2)/m2)
    if(grepl('alpha|Cen|AC|HD209100',Par$star)){
        x <- lbTC[,1]
        y <- lbTC[,2]
        if(grepl('alpha|Cen|AC',Par$star)){
            xlim <- c(-30,15)
            ylim=c(15,-30)
        }else{
            xlim <- 1.1*rev(range(x[index]),0)
            ylim <- 1.1*range(y[index],0)
#            xlim <- c(3,-3)
#            ylim <- c(-3,3)
        }
        plot(x[index],y[index],xlab=xlab,ylab=ylab,main=expression('P3: Binary motion'),type='l',xlim=xlim,ylim=ylim)
	ind <- floor(length(index)/10)
        arrows(x[index[1]],y[index[1]],x[index[ind]],y[index[ind]],length=0.1,angle=30,code=2,col='red')
        points(x[inds],y[inds],pch='+',col='red')
        text(x[inds],y[inds],labels=ts,pos=c(3,rep(2,9),rep(4,3)),xpd=NA)
        points(0,0,pch='+')
###directions
        dx <- (max(x)-min(x))*0.2
        dy <- (max(y)-min(y))*0.2
        ymin <- xmin <- -30
        ymax <- xmax <- -20
        arrows(xmin,ymin,xmax,ymin,length=0.1,angle=30,code=2,col='darkgrey',lwd=2)
        text(x=xmax,y=ymin,labels='E',pos=4,cex=1.5)
        arrows(xmin,ymin,xmin,ymax,length=0.1,angle=30,code=2,col='darkgrey',lwd=2)
        text(xmin,ymax,,labels='N',pos=1,cex=1.5)
    }else{
        x <- -lbTC[,1]
        y <- -lbTC[,2]
        if(star=='HD42581'){
            rel <- read.table('/Users/ffeng/Documents/projects/dwarfs/data/combined/HD42581/HD42581_astrometry.rel',header=TRUE)
            xlim <- range(x[index],0,rel[,'dra']/1e3)
            ylim <- range(y[index],0,rel[,'ddec']/1e3)

        }else{
            xlim <- range(x[index],0)
            ylim <- range(y[index],0)
        }
        yy <- time_Jd2yr(OutTime$BJDtcb)
        trange <- range(yy)
        ts <- unique(round(seq(round(trange[1],0),round(trange[2],0),length.out=20)))
#        ts <- c(ts,2018.37965)
        inds <- sapply(ts, function(t) which.min(abs(yy-t)))
        if(Par$star=='S2' & FALSE){
#            xlim <- c(-0.08,0.05)
#            ylim <- c(-0.05,0.20)
            xlim <- c(-0.2,0.2)
            ylim <- c(-0.2,0.2)
        }
        plot(x[index],y[index],xlab=xlab,ylab=ylab,main=expression('P3: Binary motion'),pch='.',xlim=rev(xlim),ylim=ylim)
	ind <- max(2,floor(length(index)/10))
        arrows(x[index[ind-1]],y[index[ind-1]],x[index[ind]],y[index[ind]],length=0.1,angle=30,code=2,col='grey')
        points(x[inds],y[inds],pch='+',col='red')
        text(x[inds],y[inds],labels=ts,pos=c(3,rep(2,9),rep(4,3)),xpd=NA)
        points(0,0,pch='+')
        if(star=='HD42581'){
            rel <- read.table('/Users/ffeng/Documents/projects/dwarfs/data/combined/HD42581/HD42581_astrometry.rel',header=TRUE)
            jd.ref <- time_Cal2JD(cal=cbind(2019,11,10))
            ii <- which.min(abs(rowSums(utc)-sum(jd.ref)))
            arrows(0,0,x[ii],y[ii],length=0.1,angle=30,code=2,col='blue')
            points(rel[,'dra']/1e3,rel[,'ddec']/1e3,col='red')
        }
    }
}

###lensing effects
lenSun <- OutAstro$SolarDefList$Sun
plot(lenSun[index,1]*pc2au,lenSun[index,2]*pc2au,xlab=xlab,ylab=ylab,main=expression('P4: Sun lensing'),pch=20,cex=0.2)
if(Par$PlanetShapiro){
lenEarth <- OutAstro$SolarDefList$Earth
plot(lenEarth[index,1]*pc2au,lenEarth[index,2]*pc2au,xlab=xlab,ylab=ylab,main=expression('P5: Earth lensing'),pch=20,cex=0.2)
}
if(Par$star=='alphaCenA'){
plot(-OutAstro$OffLenT[,1],OutAstro$OffLenT[,2],xlab=xlab,ylab=ylab,main=expression('P6: '*alpha*' Centauri B lensing'),pch=20,cex=0.2)
}else{
plot(OutAstro$OffLenT[index,1],OutAstro$OffLenT[index,2],xlab=xlab,ylab=ylab,main=paste('P6:',Par$star,'companion lensing'),pch=20,cex=0.2)
}

###geometric direction
lb <- gen_Xyz2lb(uOT)
if(Par$star=='alphaCenA'){
plot(lb[index,1]*180/pi,lb[index,2]*180/pi,xlab=expression(alpha*' [deg]'),ylab=expression(delta*' [deg]'),main=expression('P7: Geometric '*alpha*' Centauri A position'),pch=20,cex=0.2)
}else{
plot(lb[index,1]*180/pi,lb[index,2]*180/pi,xlab=expression(alpha*' [deg]'),ylab=expression(delta*' [deg]'),main=paste('P7: Geometric',Par$star,'position'),pch=20,cex=0.2)
}

###observed direction
lb.woRef <- gen_Xyz2lb(uOT+OutAstro$dl.woRef)
if(Par$star=='alphaCenA'){
    plot(lb.woRef[index,1]/pi*180,lb.woRef[index,2]/pi*180,xlab=expression(alpha*' [deg]'),ylab=expression(delta*' [deg]'),main=expression('P8: Observed '*alpha*' Centauri A position without refraction'),pch=20,cex=0.2)
    plot(OutAstro$DirObs[index,1]*180/pi,OutAstro$DirObs[index,2]*180/pi,xlab=expression(alpha*' [deg]'),ylab=expression(delta*' [deg]'),main=expression('P9: Observed '*alpha*' Centauri A position with refraction'),pch=20,cex=0.2)
}else{
    plot(lb.woRef[index,1]/pi*180,lb.woRef[index,2]/pi*180,xlab=expression(alpha*' [deg]'),ylab=expression(delta*' [deg]'),main=paste('P8: Observed',Par$star,'position without refraction'),pch=20,cex=0.2)
    plot(OutAstro$DirObs[index,1]*180/pi,OutAstro$DirObs[index,2]*180/pi,xlab=expression(alpha*' [deg]'),ylab=expression(delta*' [deg]'),main=paste('P9: Observed',Par$star,'position with refraction'),pch=20,cex=0.2)
}
dev.off()

if(Par$binary){
ParNew <- fit_ChangePar(Par)
OutAstroC <- astro_FullModel(OutObs,OutTime,ParNew,Mlens=Par$mT,component='C')
####relative astrometry
if(Par$Np>0){
fout <- paste0(dir.out,'relative_',fname)
cat(fout,'\n')
pdf(fout,12,12)
size <- 1.2
par(mfrow=c(3,3),mar=c(5,5,2,1),cex.lab=size)
#cosd <- cos(dir.OB[,2])
cosd <- cos(dir.OT[,2])
#cosd <- cos(AstroPar['dec2']/180*pi)
##alpha* vs delta
plot((OutAstroC$OffRef[index,1]-OutAstro$OffRef[index,1])*cosd[index],OutAstroC$OffRef[index,2]-OutAstro$OffRef[index,2],xlab=expression(Delta*alpha*'* [as]'),ylab=expression(Delta*delta*' [as]'),main='P1: Differential refraction',pch=20,cex=0.2)
####refraction
dref <- sqrt((OutAstroC$OffRef[,1]-OutAstro$OffRef[,1])^2+(OutAstroC$OffRef[,2]-OutAstro$OffRef[,2])^2)
##jd vs rho
plot(jd.utc[index],dref[index],xlab=expression(JD[UTC]),ylab=expression(rho*' [as]'),main='P2: Differential refraction',pch=20,cex=0.2)
##ele vs rho
indG <- index[which(Par$ObsInfo[index,'ObsType']=='ground')]
if(length(indG)>0){
    plot(OutTime$elevation[indG]*180/pi,dref[indG],xlab=expression(Theta*' [deg]'),ylab=expression(rho*' [as]'),main='P3: Differential refraction',pch=20,cex=0.2)
    abline(h=60,lty=2)
    ##ele vs R
    plot(OutTime$elevation[indG]*180/pi,OutAstro$Ref[indG]*pc2au,xlab=expression(Theta*' [deg]'),ylab=expression(R*' [as]'),main='P4: Atmospheric refraction',pch=20,cex=0.2)
    abline(h=60,lty=2)
}

if(Par$star=='alphaCenA'){
    plot(-((OutAstroC$OffAbe[index,1]-OutAstro$OffAbe[index,1])*cosd[index]),(OutAstroC$OffAbe[index,2]-OutAstro$OffAbe[index,2]),xlab=expression(Delta*alpha*'* [deg]'),ylab=expression(Delta*delta*' [as]'),main='P5: Differential aberration',pch=20,cex=0.2)
}else{
    plot(((OutAstroC$OffAbe[index,1]-OutAstro$OffAbe[index,1])*cosd[index]),(OutAstroC$OffAbe[index,2]-OutAstro$OffAbe[index,2]),xlab=expression(Delta*alpha*'* [deg]'),ylab=expression(Delta*delta*' [as]'),main='P5: Differential aberration',pch=20,cex=0.2)
}
dabe <- sqrt((OutAstroC$OffAbe[,1]-OutAstro$OffAbe[,1])^2+(OutAstroC$OffAbe[,2]-OutAstro$OffAbe[,2])^2)
plot(jd.utc[index],dabe[index],xlab=expression(JD[UTC]),ylab=expression(rho*' [as]'),main='P6: Differential aberration',pch=20,cex=0.2)
#lvl <- cross(uOT,cross(SO[,4:6],uOT))
#plot(sqrt(rowSums(lvl^2)),dabe,xlab=expression(Theta*' [deg]'),ylab=expression(Delta*theta*' [as]'),main='Differential aberration',pch=20,cex=0.2)
####all differetial astrometry
#duOC <-  (OutAstroC$dir-dir.OC)
#duOT <- (OutAstro$dir-dir.OT)
#dtheta <- cbind((duOC[,1]-duOT[,1])*cosd[index],duOC[,2]-duOT[,2])*pc2au
#dtheta <-  (OutAstroC$all-OutAstro$all)*pc2au
#plot(dtheta[,1],dtheta[,2],xlab=expression(Delta*alpha*'* [as]'),ylab=expression(Delta*delta*' [as]'),main='Relative motion induced by non-geometric effects',pch=20,cex=0.2)
dlenS <- OutAstroC$SolarDefList$Sun*pc2au-OutAstro$SolarDefList$Sun*pc2au
if(Par$star=='alphaCenA'){
    plot(-dlenS[index,1],-dlenS[index,2],xlab=xlab,ylab=ylab,main='P7: Differential Solar lensing',pch=20,cex=0.2)
}else{
    plot(dlenS[index,1],dlenS[index,2],xlab=xlab,ylab=ylab,main='P7: Differential Solar lensing',pch=20,cex=0.2)
}

eta <- cbind((dir.OC[,1]-dir.OT[,1])*cosd,dir.OC[,2]-dir.OT[,2])*pc2au
if(Par$star=='S2' & FALSE){
    x <- eta[index,1]
    y <- eta[index,2]
    plot(x,y,xlab=expression(Delta*alpha*'* [as]'),ylab=expression(Delta*delta*' [as]'),main='P8: Geometric orbit',pch=20,cex=0.2,xlim=rev(range(x)))#,xlim=rev(),ylim=)
#abline(h=c(-0.01,0.18),xpd=TRUE)
#abline(v=c(c(-0.07,0.05)),xpd=TRUE)
}else if(Par$star=='HD197461'){
    eta <- -eta
    plot(eta[index,1],eta[index,2],xlab=expression(Delta*alpha*'* [as]'),ylab=expression(Delta*delta*' [as]'),main='P8: Geometric orbit',pch=20,cex=0.2,xlim=rev(range(eta[index,1])))
}else{
    plot(eta[index,1],eta[index,2],xlab=expression(Delta*alpha*'* [as]'),ylab=expression(Delta*delta*' [as]'),main='P8: Geometric orbit',pch=20,cex=0.2)
}
points(0,0,pch='*',cex=3)

eta.obs <- cbind((OutAstroC$DirObs[,1]-OutAstro$DirObs[,1])*cosd,OutAstroC$DirObs[,2]-OutAstro$DirObs[,2])*pc2au
eta.obs2 <- cbind(OutAstroC$DirObs[,1]-OutAstro$DirObs[,1],OutAstroC$DirObs[,2]-OutAstro$DirObs[,2])*pc2au
plot(eta.obs[index,1],eta.obs[index,2],xlab=expression(Delta*alpha*'* [as]'),ylab=expression(Delta*delta*' [as]'),main='P9: Observed orbit',pch=20,cex=0.2)
points(0,0,pch='*',cex=3)
}
#points(eta.obs2[,1],eta.obs2[,2],pch=20,cex=0.2,col='red')
dev.off()
}

