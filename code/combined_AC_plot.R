source('OrbitFunction.R')
source('orbit.R')
source('sofa.R')
fout <- paste0('../results/AC_combined_comparison.pdf')
cat(fout,'\n')
pdf(fout,12,5)
DJM0 <- 2400000.5
size <- 1.0
par(mfrow=c(1,3),mar=c(4.2,4.2,4,0.2),cex=size,cex.lab=size,cex.axis=size,oma=rep(1,4))#oma=c(0,0,0,0)
transp <- 50
model <- out.all <- list()
load('../results/ACA_comparison.Robj')
ind <- which(out$AAT[,1] < 2002)
out$AAT <- out$AAT[ind,]
out.all[['ACA']] <- out
model[['ACA']] <- list(yr2=yr2,y2=y2,vgt=zgt,vst=zst)
load('../results/ACB_comparison.Robj')
ind <- which(out$AAT[,1] < 2002)
out$AAT <- out$AAT[ind,]
out.all[['ACB']] <- out
model[['ACB']] <- list(yr2=yr2,y2=y2,vgt=zgt,vst=zst)

rv.min <- min(unlist(sapply(1:2,function(k) unlist(sapply(names(out.all[[k]]),function(i) min(out.all[[k]][[i]][,2]))))),na.rm=TRUE)
rv.max <- max(unlist(sapply(1:2,function(k) unlist(sapply(names(out.all[[k]]),function(i) max(out.all[[k]][[i]][,2]))))),na.rm=TRUE)
for(k in 1:2){
    yr <- unlist(sapply(names(out.all[[k]]),function(i) out.all[[k]][[i]][,1]))
    RV <- unlist(sapply(names(out.all[[k]]),function(i) out.all[[k]][[i]][,2]))
    dRV <- unlist(sapply(names(out.all[[k]]),function(i) out.all[[k]][[i]][,3]))
    cat('sd(dRV)=',sd(dRV,na.rm=TRUE),'\n')
    eRV <- unlist(sapply(names(out.all[[k]]),function(i) out.all[[k]][[i]][,4]))
    cols <- c('grey','blue','green','brown','cyan','steelblue','orange','purple')
    if(k==1){
        plot(yr,RV,xlab='Time [year]',ylab='RV [km/s]',col='white',ylim=c(rv.min,rv.max))
        ts <- seq(1995,2015,by=5)
        axis(side=3,at=ts,labels=rowSums(yr2jd(ts))-DJM0)
    }else{
        points(yr,RV,col='white')
    }
    for(j in 1:length(out.all[[k]])){
        n <- names(out.all[[k]])[j]
        t <- out.all[[k]][[j]][,1]
        y <- out.all[[k]][[j]][,2]
        ey <- out.all[[k]][[j]][,4]
        points(out.all[[k]][[j]][,1],out.all[[k]][[j]][,2],col=tcol(cols[j],transp))
    }
#arrows(t,y-ey*1e-3,t,y+ey*1e-3,length=0.05,angle=90,code=3,col=tcol(cols[j],transp))
    lines(model[[k]]$yr2,model[[k]]$y2,col='red')
}

for(k in 1:2){
    plot(yr,dRV,xlab='Time [year]',ylab='O-C [m/s]',col='white')
    for(j in 1:length(out.all[[k]])){
        n <- names(out.all[[k]])[j]
        t <- out.all[[k]][[j]][,1]
        y <- out.all[[k]][[j]][,3]
        ey <- out.all[[k]][[j]][,4]
        points(t,y,col=tcol(cols[j],transp))
        arrows(t,y-ey,t,y+ey,length=0.05,angle=90,code=3,col=tcol(cols[j],transp))
        abline(h=0,col='black',lty=2)
        ts <- seq(1995,2015,by=5)
    }
    axis(side=3,at=ts,labels=rowSums(yr2jd(ts))-DJM0)
    vr <- model[[k]]$vst+model[[k]]$vgt
    lines(model[[k]]$yr2,vr-mean(vr),col='blue')
    if(k==1){
        legend('top',inset=c(0,-0.3),xpd=NA,legend=names(out),col=cols[1:length(out.all[[k]])],pch=1,horiz=TRUE,bty='n')
    }
    if(k==1){
        legend('topright',xpd=NA,legend=expression(alpha*' Cen A'),bty='n')
    }else{
        legend('topright',xpd=NA,legend=expression(alpha*' Cen B'),bty='n')
    }
}
dev.off()
