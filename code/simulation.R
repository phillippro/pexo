utc <- time_ChangeBase(cbind(Data[,1],0))
for(n in names(ParOpt)) Par[[n]] <- ParOpt[n]
ParT <- ParC <- Par
ParC$Nepoch <- ParT$Nepoch <- nrow(utc)
OutBary <- time_Utc2tb(utc,ParT)
OutTd <- time_Ta2te(OutBary,ParT)
RVTd <- rv_FullModel(OutBary,OutTd,ParT)

ParC$mC1 <- Par$mT
ParC$mT <- Par$mC1
ParC$omegaT1 <- Par$omegaT1-pi
OutCd <- time_Ta2te(OutBary,ParC)
RVCd <- rv_FullModel(OutBary,OutCd,ParC)

###simulated epochs
utc <- seq(min(Data[,1]),max(Data[,1]),length.out=1000)
utc <- time_ChangeBase(cbind(utc,0))
ParT <- ParC <- Par
ParC$Nepoch <- ParT$Nepoch <- nrow(utc)
OutBary <- time_Utc2tb(utc,ParT)
OutT <- time_Ta2te(OutBary,ParT)
RVT <- rv_FullModel(OutBary,OutT,ParT)

ParC$mC1 <- Par$mT
ParC$mT <- Par$mC1
ParC$omegaT1 <- Par$omegaT1-pi
OutC <- time_Ta2te(OutBary,ParC)
RVC <- rv_FullModel(OutBary,OutC,ParC)

source('plot_simulation.R')
