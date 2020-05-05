#aca <- read.table('../input/HD128620v1/HD128620_HARPS.dat',header=TRUE)
#acb <- read.table('../input/HD128621v1/HD128621_HARPS.dat',header=TRUE)
#acb <- read.table('/Users/ffeng/Documents/projects/dwarfs/data/aperture/alphaCenA/alphaCenA_HARPS.dat',header=TRUE)
#aca <- read.table('/Users/ffeng/Documents/projects/dwarfs/data/aperture/alphaCenB/alphaCenB_HARPS.dat',header=TRUE)
acb <- read.table('/Users/ffeng/Documents/projects/dwarfs/data/aperture/HD128621_ew3_v1/HD128621_TERRA_1AP1_erv.dat',header=TRUE)
aca <- read.table('/Users/ffeng/Documents/projects/dwarfs/data/aperture/HD128620/HD128620_TERRA_1AP1_erv.dat',header=TRUE)
ts <- c(aca[,1],acb[,1])+2400000
utc <- seq(min(ts),max(ts),length.out=1000)
utc <- time_ChangeBase(cbind(utc,0))
ParNew <- Par
if(TRUE){
ParIni <- ParOpt
ParIni['mC1'] <- 1.1
}else{
kepini <- Par$KepIni
names(kepini) <- paste0(names(Par$KepIni),'1')
ParIni <- c(Par$Ini,kepini)
}
for(n in names(ParIni)) ParNew[[n]] <- ParIni[n]
ParT <- ParC <- ParNew
ParC$Nepoch <- ParT$Nepoch <- nrow(utc)
OutBary <- time_Utc2tb(utc,ParT)
OutT <- time_Ta2te(OutBary,ParT)
RVT <- rv_FullModel(OutBary,OutT,ParT)
##
ParC$mC1 <- ParNew$mT
ParC$mT <- ParNew$mC1
ParC$omegaT1 <- Par$omegaT1-pi
OutC <- time_Ta2te(OutBary,ParC)
RVC <- rv_FullModel(OutBary,OutC,ParC)

##
jdT <- rowSums(OutT$BJDtdb)
jdC <- rowSums(OutC$BJDtdb)
rvT <- RVT$RvST-RVT$RvST[1]+acb[1,2]
rvC <- RVC$RvST-RVC$RvST[1]+aca[1,2]

pdf('ac_test.pdf',16,16)
par(mfrow=c(4,4))
plot(aca[,1]+2400000,aca[,2],xlab='BJD',ylab='rv',ylim=range(aca[,2],rvC))
points(jdC,rvC,col='red')

plot(acb[,1]+2400000,acb[,2],xlab='BJD',ylab='rv',ylim=range(acb[,2],rvT))
points(jdT,rvT,col='red')
dev.off()
