if(exists('OutBary')){
    OutBary$SO[,1:3] <- OutBary$SO[,1:3]/IFTE.K
    OutBary$SG[,1:3] <- OutBary$SG[,1:3]/IFTE.K
    OutBary$GO[,1:3] <- OutBary$GO[,1:3]/IFTE.K
    OutBary$MO[,1:3] <- OutBary$MO[,1:3]/IFTE.K
    OutBary$GM[,1:3] <- OutBary$GM[,1:3]/IFTE.K
    OutBary$TDBmTTgeo <- OutBary$TDBmTTgeo/IFTE.K
    OutBary$dzenith <- OutBary$dzenith*IFTE.K
}

if(exists('OutTime')){
    OutTime$RoemerTarget <- OutTime$RoemerTarget/IFTE.K
    OutTime$RoemerSolar <- OutTime$RoemerSolar/IFTE.K
    OutTime$RoemerSB <- OutTime$RoemerSB/IFTE.K
    OutTime$RoemerT2$all <- OutTime$RoemerT2$all/IFTE.K
    OutTime$EinsteinTarget <- OutTime$EinsteinTarget/IFTE.K
    OutTime$EinsteinIS <- OutTime$EinsteinIS/IFTE.K
    OutTime$VacuumIS <- OutTime$VacuumIS/IFTE.K
    OutTime$ShapiroTarget <- OutTime$EinsteinTarget/IFTE.K
    OutTime$ShapiroSolar <- OutTime$ShapiroSolar/IFTE.K
    OutTime$TropoDelay <- OutTime$TropoDelay/IFTE.K
    OutTime$TropoDelayT2 <- OutTime$TropoDelayT2/IFTE.K
    OutTime$TargetDelay <- OutTime$TargetDelay/IFTE.K
    OutTime$delevation <- OutTime$delevation*IFTE.K
    OutTime$delevationT2 <- OutTime$delevationT2*IFTE.K
    OutTime$rOT <- OutTime$rOT/IFTE.K
    OutTime$rOB <- OutTime$rOB/IFTE.K
    OutTime$rST <- OutTime$rST/IFTE.K
    OutTime$rBT <- OutTime$rBT/IFTE.K
    OutTime$SB[,1:3] <- OutTime$SB[,1:3]/IFTE.K
    for(n in names(OutTime$Eph)){
        OutTime$Eph[[n]][,1:3] <- OutTime$Eph[[n]][,1:3]/IFTE.K
    }
    OutTime$OutBT$BT[,1:3] <- OutTime$OutBT$BT[,1:3]/IFTE.K
    OutTime$OutBT$RBT <- OutTime$OutBT$RBT/IFTE.K
    for(n in names(OutTime$RoemerOrder)){
        OutTime$RoemerOrder[[n]] <- OutTime$RoemerOrder[[n]]/IFTE.K
    }
    for(n in names(OutTime$OL)){
        OutTime$OL[[n]][,c('x.au','y.au','z.au','ROL','SLx.km','SLy.km','SLz.km')] <- OutTime$OL[[n]][,c('x.au','y.au','z.au','ROL','SLx.km','SLy.km','SLz.km')]/IFTE.K
    }
    OutTime$rTC <- OutTime$rTC/IFTE.K
    OutTime$rOC <- OutTime$rOC/IFTE.K
    OutTime$rSC <- OutTime$rSC/IFTE.K
}
##There is no length, time and mass quantities in OutAstroT and OutRv list object and thus no convertion is made.

