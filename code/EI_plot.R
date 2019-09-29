pdf('EI_orbit.pdf',6,6)
size <- 1.2
par(mar=c(5,5,1,1),cex=size,cex.lab=size,cex.axis=size)

uSC <- calUnit(rSC)
dir.ST <- xyz2lb.vec(uST[,1:3])
dir.SC <- xyz2lb.vec(uSC[,1:3])
cosd <- cos(dir.SC[,2])
eta <- cbind((dir.ST[,1]-dir.SC[,1])*cosd,dir.ST[,2]-dir.SC[,2])*pc2au
ind0 <- which.min(abs(T0-rowSums(bjd.tdb)))#tmin for RV
ind1 <- which.min(abs(2457993.71-rowSums(bjd.tdb)))#tmax for RV
plot(eta[,1],eta[,2],xlab=expression(Delta*alpha*'* [as]'),ylab=expression(Delta*delta*' [as]'),main='',pch=20,cex=0.2,xlim=rev(range(eta[,1])))#,xlim=c(-3,3.2),ylim=c(-3,3.2))
x <- eta[,1]
y <- eta[,2]
dx <- (max(x)-min(x))*0.2
dy <- (max(y)-min(y))*0.2
ymin <- min(y)
xmin <- min(x)
ymax <- max(y)
xmax <- max(x)
arrows(xmax-dx,ymin,xmax,ymin,length=0.1,angle=30,code=2,col='darkgrey',lwd=2)
text(x=xmax,y=ymin,labels='E',pos=3,cex=1.5,col='darkgrey')
arrows(xmax-dx,ymin,xmax-dx,ymin+dy,length=0.1,angle=30,code=2,col='darkgrey',lwd=2)
text(xmax-dx,ymin+dy,,labels='N',pos=2,cex=1.5,col='darkgrey')

yrs <- seq(1990,2025,by=5)
jds <- rowSums(yr2jd(yrs))
bjd.gaia <- sum(tpos)
bjd.hip <- sum(tpos)
#yr.di <- c(2004.7,2008.5)
yr.di <- 2008.5
jd.di <- rowSums(yr2jd(yr.di))
inds <- sapply(jds,function(jd) which.min(abs(jd-rowSums(bjd.tdb))))
ind.di <- sapply(jd.di,function(jd) which.min(abs(jd-rowSums(bjd.tdb))))
points(eta[inds,1],eta[inds,2],pch='x',col='red')
points(0,0,pch='+',col='black')
#text(1.2*eta[inds,1],1.2*eta[inds,2],labels=yrs,col='red',pos=4)
text(eta[inds,1],eta[inds,2],labels=yrs,col='red',pos=c(4,4,4,2,2,4,1,4))
rhos <- sqrt(eta[inds,1]^2+eta[inds,2]^2)
cat('yrs=',yrs,'\n')
cat('rhos=',round(rhos,2),'\n')
points(eta[ind.di,1],eta[ind.di,2],pch='x',col='steelblue')
text(eta[ind.di,1],eta[ind.di,2],labels=yr.di,col='steelblue',pos=2)
rho <- sqrt(eta[ind.di,1]^2+eta[ind.di,2]^2)
cat('separation at Janson epoch:',round(rho,2),'\n')
dev.off()
