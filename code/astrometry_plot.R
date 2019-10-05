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
rd.OTobs <- OutAstroT$DirObs*180/pi#deg
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
fname <- paste0(Par$star,'_astrometry_',Par$BinaryModel,'_dt',gsub('\\.','',dt),'day_Ntime',Ntime,'_',Par$RefType,'.pdf')
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
plot(OutAstroT$OffRef[index,1],OutAstroT$OffRef[index,2],xlab=xlab,ylab=ylab,main=expression('P1: Atmospheric refraction'),pch=20,cex=0.2)

##aberration
plot(-OutAstroT$OffAbe[index,1],OutAstroT$OffAbe[index,2],xlab=xlab,ylab=ylab,main=expression('P2: Stellar aberration'),pch=20,cex=0.2)

##Binary motion
if(Par$binary){
    delta <- gen_Xyz2lb(OutTime$uOT)[,2]
    lbBT <- gen_CalOffset(OutTime$uST,OutTime$uSB,bref=delta)
    lbBC <- gen_CalOffset(gen_CalUnit(rSC),OutTime$uSB,bref=delta)
    lbTC <- gen_CalOffset(gen_CalUnit(rSC),OutTime$uST,bref=delta)
####prepare markers
    yy <- time_Jd2yr(utc)
#    ts <- seq(1980,2050,by=10)
    ts <- seq(1980,2040,by=5)
    ts <- ts[ts<max(yy)]
    inds <- sapply(ts, function(t) which.min(abs(yy-t)))
####plot binary orbit to be comparable with Fig. 1 in Pourbaix et al. 1999
                                        #x1 <- (-x*(m1+m2)/m2)
                                        #y1 <- (-y*(m1+m2)/m2)
    if(grepl('alpha|Cen|AC',Par$star)){
        x <- lbTC[,1]
        y <- lbTC[,2]
        plot(x[index],y[index],xlab=xlab,ylab=ylab,main=expression('P3: Binary motion'),type='l',xlim=c(-30,15),ylim=c(15,-30))
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
        xlim <- range(x[index])
        ylim <- range(y[index])
        yy <- time_Jd2yr(OutTime$BJDtcb)
        trange <- range(yy)
        ts <- unique(round(seq(round(trange[1],0),round(trange[2],0),length.out=30)))
        ts <- c(ts,2018.37965)
        inds <- sapply(ts, function(t) which.min(abs(yy-t)))
        if(Par$star=='S2' & FALSE){
#            xlim <- c(-0.08,0.05)
#            ylim <- c(-0.05,0.20)
            xlim <- c(-0.2,0.2)
            ylim <- c(-0.2,0.2)
        }
        plot(x[index],y[index],xlab=xlab,ylab=ylab,main=expression('P3: Binary motion'),pch='.',xlim=rev(xlim))
        if(Par$star=='S2' & FALSE){
            astro <- read.table('../data/S2_astrometry.dat')
            if(TRUE){
                x0 <- 0.99#mas
                y0 <- -0.85#mas
                vx0 <- -0.060#mas/yr
                vy0 <- 0.221#mas/yr
                dalpha0 <- -astro[,3]
                ddelta0 <- astro[,4]
                rho <- sqrt(dalpha0^2+ddelta0^2)
                psi <- atan2(astro[,4],-astro[,3])
                da <- -(x0+vx0*(astro[,2]+DJM0-Par$tpos)/DJY)*1e-3
                dd <- (y0+vy0*(astro[,2]+DJM0-Par$tpos)/DJY)*1e-3
                da <- 0
                dd <- 0
                dalpha <- dalpha0+da
                ddelta <- ddelta0+dd
            }
            points(dalpha,ddelta,pch=20,cex=0.5,col='blue')
            ind3 <- which.min(abs(astro[,2]+DJM0-Par$Tp))
            points(dalpha[ind3],ddelta[ind3],pch=20,cex=0.5,col='green')
            points(x[inds],y[inds],pch='+',col='red')
            text(x[inds],y[inds],labels=ts,pos=c(3,rep(2,9),rep(4,3)))
            points(0,0,pch='+')
            ind1 <- which.min(abs(OutTime$OutBT$U-Par$Omega))
            ind2 <- which.min(abs(OutTime$OutBT$U-(Par$Omega+pi)%%(2*pi)))
            lines(x[c(ind1,ind2)],y[c(ind1,ind2)],lty=2)
            plot(rho,psi)
            rho.model <- sqrt(x^2+y^2)
            psi.model <- atan2(y,x)
            points(rho.model,psi.model,col='blue',pch='.')
            tt <- time_Jd2yr(cbind(DJM0+astro[,2],0))
            tmodel <- time_Jd2yr(OutTime$tB)
            plot(tmodel,rho.model,xlab='tB',ylab='rho [as]')
            points(tt,rho,col='blue',pch='.')
            plot(tmodel,psi.model*180/pi,xlab='tB',ylab='psi [deg]')
            points(tt,psi*180/pi,col='blue',pch='.')
        }
    }
}

###lensing effects
lenSun <- OutAstroT$SolarDefList$Sun
lenEarth <- OutAstroT$SolarDefList$Earth
plot(lenSun[index,1]*pc2au,lenSun[index,2]*pc2au,xlab=xlab,ylab=ylab,main=expression('P4: Sun lensing'),pch=20,cex=0.2)
plot(lenEarth[index,1]*pc2au,lenEarth[index,2]*pc2au,xlab=xlab,ylab=ylab,main=expression('P5: Earth lensing'),pch=20,cex=0.2)
if(Par$star=='alphaCenA'){
plot(-OutAstroT$OffLenT[,1],OutAstroT$OffLenT[,2],xlab=xlab,ylab=ylab,main=expression('P6: '*alpha*' Centauri B lensing'),pch=20,cex=0.2)
}else{
plot(OutAstroT$OffLenT[index,1],OutAstroT$OffLenT[index,2],xlab=xlab,ylab=ylab,main=paste('P6:',Par$star,'companion lensing'),pch=20,cex=0.2)
}

###geometric direction
lb <- gen_Xyz2lb(uOT)
if(Par$star=='alphaCenA'){
plot(lb[index,1]*180/pi,lb[index,2]*180/pi,xlab=expression(alpha*' [deg]'),ylab=expression(delta*' [deg]'),main=expression('P7: Geometric '*alpha*' Centauri A position'),pch=20,cex=0.2)
}else{
plot(lb[index,1]*180/pi,lb[index,2]*180/pi,xlab=expression(alpha*' [deg]'),ylab=expression(delta*' [deg]'),main=paste('P7: Geometric',Par$star,'position'),pch=20,cex=0.2)
}

###observed direction
lb.woRef <- gen_Xyz2lb(uOT+OutAstroT$dl.woRef)
if(Par$star=='alphaCenA'){
    plot(lb.woRef[index,1]/pi*180,lb.woRef[index,2]/pi*180,xlab=expression(alpha*' [deg]'),ylab=expression(delta*' [deg]'),main=expression('P8: Observed '*alpha*' Centauri A position without refraction'),pch=20,cex=0.2)
    plot(OutAstroT$DirObs[index,1]*180/pi,OutAstroT$DirObs[index,2]*180/pi,xlab=expression(alpha*' [deg]'),ylab=expression(delta*' [deg]'),main=expression('P9: Observed '*alpha*' Centauri A position with refraction'),pch=20,cex=0.2)
}else{
    plot(lb.woRef[index,1]/pi*180,lb.woRef[index,2]/pi*180,xlab=expression(alpha*' [deg]'),ylab=expression(delta*' [deg]'),main=paste('P8: Observed',Par$star,'position without refraction'),pch=20,cex=0.2)
    plot(OutAstroT$DirObs[index,1]*180/pi,OutAstroT$DirObs[index,2]*180/pi,xlab=expression(alpha*' [deg]'),ylab=expression(delta*' [deg]'),main=paste('P9: Observed',Par$star,'position with refraction'),pch=20,cex=0.2)
}
dev.off()

####relative astrometry
if(Par$binary){
fout <- paste0(dir.out,'relative_',fname)
cat(fout,'\n')
pdf(fout,12,12)
size <- 1.2
par(mfrow=c(3,3),mar=c(5,5,2,1),cex.lab=size)
#cosd <- cos(dir.OB[,2])
cosd <- cos(dir.OT[,2])
#cosd <- cos(AstroPar['dec2']/180*pi)
##alpha* vs delta
plot((OutAstroC$OffRef[index,1]-OutAstroT$OffRef[index,1])*cosd[index],OutAstroC$OffRef[index,2]-OutAstroT$OffRef[index,2],xlab=expression(Delta*alpha*'* [as]'),ylab=expression(Delta*delta*' [as]'),main='P1: Differential refraction',pch=20,cex=0.2)
####refraction
dref <- sqrt((OutAstroC$OffRef[,1]-OutAstroT$OffRef[,1])^2+(OutAstroC$OffRef[,2]-OutAstroT$OffRef[,2])^2)
##jd vs rho
plot(jd.utc[index],dref[index],xlab=expression(JD[UTC]),ylab=expression(rho*' [as]'),main='P2: Differential refraction',pch=20,cex=0.2)
##ele vs rho
if(Par$ObsType=='ground'){
plot(OutTime$elevation[index]*180/pi,dref[index],xlab=expression(Theta*' [deg]'),ylab=expression(rho*' [as]'),main='P3: Differential refraction',pch=20,cex=0.2)
abline(h=60,lty=2)
##ele vs R
plot(OutTime$elevation[index]*180/pi,OutAstroT$Ref[index]*pc2au,xlab=expression(Theta*' [deg]'),ylab=expression(R*' [as]'),main='P4: Atmospheric refraction',pch=20,cex=0.2)
abline(h=60,lty=2)
}

if(Par$star=='alphaCenA'){
    plot(-((OutAstroC$OffAbe[index,1]-OutAstroT$OffAbe[index,1])*cosd[index]),(OutAstroC$OffAbe[index,2]-OutAstroT$OffAbe[index,2]),xlab=expression(Delta*alpha*'* [deg]'),ylab=expression(Delta*delta*' [as]'),main='P5: Differential aberration',pch=20,cex=0.2)
}else{
    plot(((OutAstroC$OffAbe[index,1]-OutAstroT$OffAbe[index,1])*cosd[index]),(OutAstroC$OffAbe[index,2]-OutAstroT$OffAbe[index,2]),xlab=expression(Delta*alpha*'* [deg]'),ylab=expression(Delta*delta*' [as]'),main='P5: Differential aberration',pch=20,cex=0.2)
}
dabe <- sqrt((OutAstroC$OffAbe[,1]-OutAstroT$OffAbe[,1])^2+(OutAstroC$OffAbe[,2]-OutAstroT$OffAbe[,2])^2)
plot(jd.utc[index],dabe[index],xlab=expression(JD[UTC]),ylab=expression(rho*' [as]'),main='P6: Differential aberration',pch=20,cex=0.2)
#lvl <- cross(uOT,cross(SO[,4:6],uOT))
#plot(sqrt(rowSums(lvl^2)),dabe,xlab=expression(Theta*' [deg]'),ylab=expression(Delta*theta*' [as]'),main='Differential aberration',pch=20,cex=0.2)
####all differetial astrometry
#duOC <-  (OutAstroC$dir-dir.OC)
#duOT <- (OutAstroT$dir-dir.OT)
#dtheta <- cbind((duOC[,1]-duOT[,1])*cosd[index],duOC[,2]-duOT[,2])*pc2au
#dtheta <-  (OutAstroC$all-OutAstroT$all)*pc2au
#plot(dtheta[,1],dtheta[,2],xlab=expression(Delta*alpha*'* [as]'),ylab=expression(Delta*delta*' [as]'),main='Relative motion induced by non-geometric effects',pch=20,cex=0.2)
dlenS <- OutAstroC$SolarDefList$Sun*pc2au-OutAstroT$SolarDefList$Sun*pc2au
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

eta.obs <- cbind((OutAstroC$DirObs[,1]-OutAstroT$DirObs[,1])*cosd,OutAstroC$DirObs[,2]-OutAstroT$DirObs[,2])*pc2au
eta.obs2 <- cbind(OutAstroC$DirObs[,1]-OutAstroT$DirObs[,1],OutAstroC$DirObs[,2]-OutAstroT$DirObs[,2])*pc2au
plot(eta.obs[index,1],eta.obs[index,2],xlab=expression(Delta*alpha*'* [as]'),ylab=expression(Delta*delta*' [as]'),main='P9: Observed orbit',pch=20,cex=0.2)
points(0,0,pch='*',cex=3)
#points(eta.obs2[,1],eta.obs2[,2],pch=20,cex=0.2,col='red')
dev.off()
}

####compare with data
if(FALSE){
fout <- paste0(dir.out,'comparison_',fname)
cat(fout,'\n')
pdf(fout,8,4)
size <- 1
par(mar=c(5,5,1,1),cex.lab=size)
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE), widths=c(1,1), heights=c(3,3))

tab <- read.table('~/Documents/projects/dwarfs/data/other_sets/alpha_cen_old/k17_astrometry.dat',header=TRUE)
##compare geometric binary orbit
rho <- sqrt(rowSums(eta^2))
theta <- gen_Xy2phi(eta[,1],eta[,2])*180/pi
mjds <- rowSums(time_Yr2jd(tab[,1]))-DJM0
mjd <- jd.utc-DJM0
ind <- which(mjd>44000 & mjd<59000)

mjd <- mjd[ind]
rho <- rho[ind]
theta <- theta[ind]
eta <- eta[ind,]

alpha1 <- eta[,1]
delta1 <- eta[,2]
alpha2 <- tab[,'rho']*cos(tab[,'theta']/180*pi)
delta2 <- tab[,'rho']*sin(tab[,'theta']/180*pi)

alpha.fun <- approxfun(mjd,alpha1)
delta.fun <- approxfun(mjd,delta1)
dalpha <- alpha2-alpha.fun(mjds)
ddelta <- delta2-delta.fun(mjds)


rho.fun <- approxfun(mjd,rho)
theta.fun <- approxfun(mjd,theta)
drho <- tab[,'rho']-rho.fun(mjds)
dtheta <- tab[,'theta']-theta.fun(mjds)
col.model <- 'red'
col.data <- 'black'
len <- 0.03
point.size <- 0.8

par(mar=c(0,5,1,1))
plot(mjd[index],rho[index],xlab='MJD',ylab='Angular separation [as]',pch=20,cex=0.2,xaxt='n',col=col.model)
points(mjds,tab[,'rho'],col=col.data,pch=20,cex=point.size)
lines(mjd,rho,col=col.model)
arrows(mjds,tab[,'rho']-tab[,'erho'],mjds,tab[,'rho']+tab[,'erho'],length=len,angle=90,code=3,col=col.data)

par(mar=c(0,5,1,1))
plot(mjd[index],theta[index],xlab='MJD',ylab='Position angle of B [deg]',pch=20,cex=0.2,xaxt='n',col=col.model)
lines(mjd,theta,col=col.model)
points(mjds,tab[,'theta'],col=col.data,pch=20,cex=point.size)
arrows(mjds,tab[,'theta']-tab[,'etheta'],mjds,tab[,'theta']+tab[,'etheta'],length=len,angle=90,code=3,col=col.data)

par(mar=c(5,5,0,1))
plot(mjds[index],drho[index],xlab='MJD',ylab='O-C [as]',col=col.data,pch=20,cex=point.size,ylim=range(1.5*min(drho),1.5*max(drho),0))
abline(h=0,lty=2,col='grey')
arrows(mjds,drho-tab[,'erho'],mjds,drho+tab[,'erho'],length=len,angle=90,code=3,col=col.data)
legend('topright',legend=paste0('RMS=',round(sd(drho),1),'as'),bty='n')

par(mar=c(5,5,0,1))
plot(mjds[index],dtheta[index],xlab='MJD',ylab='O-C [deg]',col=col.data,pch=20,cex=point.size,ylim=range(1.5*min(dtheta),3*max(dtheta),0))
abline(h=0,lty=2,col='grey')
arrows(mjds,dtheta-tab[,'etheta'],mjds,dtheta+tab[,'etheta'],length=len,angle=90,code=3,col=col.data)
legend('topright',legend=paste0('RMS=',round(sd(dtheta),1),'deg'),bty='n')

dev.off()
}
