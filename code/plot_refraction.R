fout <- paste0('refracor_',target,'.pdf')
cat(fout,'\n')
pdf(fout,16,16)
#data <- read.table('../../data/PFS/PFS04-27-19/HD10700_PFS.dat',header=TRUE)
fin <- paste0('../../data/PFS/PFS04-27-19/HD',target,'_PFS.dat')
if(!file.exists(fin)){
fin <- paste0('../../data/PFS/PFS04-27-19/',target,'_PFS.dat')
}
if(set=='PFS'){
data <- read.table(fin,header=TRUE)
}else{
data <- tab
}

par(mfrow=c(4,4))
###sorting
jd.utc <- rowSums(utc)#%%2400000
ind.sort <- sort(jd.utc,index.return=TRUE)$ix
jd.utc <- jd.utc[ind.sort]
Ratm <- Ratm[ind.sort]
RV.ref <- RV.ref[ind.sort]
zen <- zen[ind.sort]
elevation <- elevation[ind.sort]
delevation <- delevation[ind.sort]
zenith <- zenith[ind.sort,]

ind2 <- ind1 <- c()
if(length(jd.utc)!=nrow(data)){
    if(length(jd.utc)>nrow(data)){
        N <- nrow(data)
        tref <- data[,1]
    }else{
        N <- nrow(jd.utc)
        tref <- jd.utc
    }
    for(i in 1:N){
        ind1 <- c(ind1,which.min(abs(jd.utc-tref[i])))
        ind2 <- c(ind2,which.min(abs(data[,1]-tref[i])))
    }
    jd.utc <- jd.utc[ind1]
    Ratm <- Ratm[ind1]
    RV.ref <- RV.ref[ind1]
    zen <- zen[ind1]
    elevation <- elevation[ind1]
    delevation <- delevation[ind1]
    data <- data[ind2,]
}
###loading
#data[,1] <- data[,1]#%%2400000
if(target=='39091'){
    tmin <- 2458407.784-0.1
    tmax <- 2458407.80+0.1
}else{
    tmin <- c()
    tmax <- c()
}
dt <- max(sort(diff(jd.utc))[round(0.6*length(jd.utc))],0.05)
ind <- which(diff(jd.utc)>dt)
if(length(ind)>0){
    ii <- sort(diff(ind),decreasing=TRUE,index.return=TRUE)$ix[1:2]
    tmin <- c(tmin,jd.utc[ind[ii]+1]-0.005)
    tmax <- c(tmax,jd.utc[ind[ii+1]]+0.005)
}

###overall view
cols <- c('blue','green','orange','brown')
for(k in 1:(length(tmin)+1)){
    if(k==1){
        ind1 <- 1:nrow(data)
        ind2 <- 1:length(jd.utc)
        t1 <- data[ind1,1]
        t2 <- jd.utc[ind2]
        xlab <- 'JD'
    }else{
        ind1 <-  which(tmin[k-1]<= data[,1] & data[,1]<=tmax[k-1])
        ind2 <- which(tmin[k-1]<=jd.utc & jd.utc<=tmax[k-1])
        t1 <- (data[ind1,1]-tmin[k-1])*24#hr
        t2 <- (jd.utc[ind2]-tmin[k-1])*24
        xlab <- 'Dt [hour]'
    }
    RV <- data[ind1,2]
    dRV <- RV.ref[ind2]
    RV.new <- RV-dRV

    plot(t1,RV,xlab=xlab,ylab='RV[m/s]',main=paste0('RMS=',round(sd(RV),2)))
    if(k==1){
        abline(v=(tmin+tmax)/2,col=cols[1:length(tmin)])
    }

    plot(t2,dRV,xlab=xlab,ylab='dRV[m/s]',main=paste0('RMS=',round(sd(RV),2)))
    if(k==1){
        abline(v=(tmin+tmax)/2,col=cols[1:length(tmin)])
    }

    plot(dRV,RV,xlab='dRV[m/s]',ylab='RV[m/s]',main=paste0('r=',round(cor(data[,2],RV.ref),2)))

    plot(t1,RV.new,xlab=xlab,ylab='RV-dRV [m/s]',main=paste0('RMS=',round(sd(RV.new),2)))
    if(k==1){
        abline(v=(tmin+tmax)/2,col=cols[1:length(tmin)])
    }
}
dev.off()
