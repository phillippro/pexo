kepini <- Par$KepIni
names(kepini) <- paste0(names(Par$KepIni),'1')
ParFitB <- ParFitA <- ParIni <- c(Par$Ini,kepini)
ParA <- ParB <- Par
for(n in names(ParFitB)) ParB[[n]] <- ParFitB[n]

OutBary <- time_Utc2tb(utc,ParB)
OutTime <- time_Ta2te(OutBary,ParB)
tmp <- fit_LogLike(Data,RateBary,ParFitB,ParB)
modelB <- tmp$model

for(n in names(ParFitA)) ParA[[n]] <- ParFitA[n]
ParA$stars <- rev(Par$stars)
ParA$star <- Par$stars[2]
ParA$secondary <- Par$stars[1]
ParA$mT <- Par$mC1
ParFitA['mC1'] <- ParA$mC1 <- Par$mT
ParFitA['omegaT1'] <- ParA$omegaT1 <- (Par$omegaT1-pi)%%(2*pi)
OutTime <- time_Ta2te(OutBary,ParA)
tmp <- fit_LogLike(Data,RateBary,ParFitA,ParA)
modelA <- tmp$model

pdf('debug_figure.pdf',8,8)
par(mfrow=c(2,2))
for(star in Par$stars){
    ind <- which(Data$star==star)
    drv <- Data[ind,2]-modelB[ind,2]
    plot(Data[ind,1],drv,xlab='BJD',ylab='rv',main=paste(star,'for ',Par$stars[1],'as target;\n RMS=',round(sd(drv),2),'m/s'))

    drv <- Data[ind,2]-modelA[ind,2]
    plot(Data[ind,1],drv,xlab='BJD',ylab='rv',main=paste(star,'for ',Par$stars[2],'as target;\n RMS=',round(sd(drv),2),'m/s'))
}
dev.off()
