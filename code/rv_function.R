###################################################################################
##This file includes functions for relativistic astrometry
################################################################################
rv_GravTarget <- function(M,r,g=1){
####################################
## Calculate gravitational redshift
##
## Input:
##   M - Lens mass (Msun)
##   r - Distance from the target to the companion (au)
##   g - gamma; g=1 for GR; g!=1 for other gravitational theories.
##
## Output:
##   Gravitational Doppler shift
####################################
    sr <- SRS*M*(1+g)/4
    sr/r
}

rv_GravSolar <- function(OL,g=1){
####################################
## Calculate gravitational redshift
##
## Input:
##   OL - Position and velocity vectors from the observer to the Solar System lenses
##   g - gamma; g=1 for GR; g!=1 for other gravitational theories.
##
## Output:
##   all - All gravitational Doppler shift
##   Zlist - Gravitational Doppler shift for individual sources
####################################
    ns <- names(OL)
    Zlist <- list()
    Zall <- 0
    for(n in ns){
        ROL <- sqrt(rowSums(OL[[n]][,1:3]^2))#au
        Z <- SRS*Mssp[n]*(1+g)/4/ROL
        Zall <- Z+Zall
        Zlist[[n]] <- Z
    }
    list(Zall=Zall,Zlist=Zlist)
}

rv_LenTarget <- function(ML,vSL,RLT,rOT,vSO,rOL,vST,g=1){
####################################
## Calculate lensing redshift in the Solar System
##
## Input:
##   ML - Lens mass (Msun)
##   vSL - Lens position relative the SSB
##   RLT - Distance from the lens to the target
##   rOT - Position vector from the observer to the target
##   vSO - Velocity of the observer with respect to the SSB
##   rOL - Position of the lens with respect to the observer
##   vST - Velocity of the the target with respect to the SSB
##   g - gamma; g=1 for GR; g!=1 for other gravitational theories.
##
## Output:
##   Lensing doppler shift
####################################
    ROT <- gen_CalLen(rOT)
    uOT <- gen_CalUnit(rOT)
    ROL <- sqrt(rowSums(rOL^2))
    lambda <- pracma::cross(uOT,pracma::cross(rOL,uOT))
    alpha <- (1+g)*SRS*ML*lambda/rowSums(lambda^2)#alpha
    dZ <- (vSL-RLT/ROT*vSO-ROL/ROT*vST)/Cauyr
    rowSums(dZ*alpha)
}

rv_LenSolar <- function(OL,rOT,vST,vSO,g=1,LenRVmethod='T2'){
####################################
## Calculate lensing redshift in the Solar System
## The PEXO method is given by eqs. 64 - 67 in Feng et al. 2019
## The T2 method is the time derivative of eqn. 34 in Feng et al. 2019 which is compatible with TEMPO2 method. The formula is as follows,
## Zshapiro = -(1+g)*G*m_i/c^3*(uOL*vOL-uOT*vOL)/(RSO*(1-cos(psi)))
##
## Input:
##   OL - Position and velocity vectors from the observer to the Solar System lenses
##   rOT - Position vectors from the observer to the target (pc)
##   vST - Velocity vector from the SSB to the target (au/yr)
##   vSO - Velocity vector from the SSB to the observer
##   g - gamma; g=1 for GR; g!=1 for other gravitational theories.
##
## Output:
##   all - All lensing Doppler shift
##   Zlist - Lensing Doppler shift for individual sources
####################################
    ns <- names(OL)
    Zlist <- list()
    Zall <- 0
    ROT <- gen_CalLen(rOT)*pc2au
    uOT <- gen_CalUnit(rOT)
    for(n in ns){
#        lambda <- pracma::cross(uOT,pracma::cross(OL[[n]][,c('x.au','y.au','z.au')],uOT))#au
        rOL <- OL[[n]][,c('x.au','y.au','z.au')]
        vOL <- OL[[n]][,c('vx.auyr','vy.auyr','vz.auyr')]
        cpsi <-  OL[[n]][,'cpsi']
        ROL <-  OL[[n]][,'ROL']
        if(LenRVmethod=='PEXO'){
            lambda <- rOL-rowSums(rOL*uOT)*uOT
            alpha <- (1+g)*SRS*Mssp[n]*lambda/rowSums(lambda^2)#alpha
            RLT <- gen_CalLen(rOT*pc2au-OL[[n]][,c('x.au','y.au','z.au')])#au
            vSL <- OL[[n]][,c('SLvx.kms','SLvy.kms','SLvz.kms')]/auyr2kms#au/yr
            dZ <- (vSL-RLT/ROT*vSO-OL[[n]][,'ROL']/ROT*vST)/Cauyr
            Z <- rowSums(dZ*alpha)
        }else{
#            ROL <- gen_CalLen(rOL)
            uOL <- gen_CalUnit(rOL)
            c1 <- rowSums(uOL*vOL)#au/yr
            c2 <- rowSums(uOT*vOL)#au/yr
            Amp <- 0.5*(1+g)*SRSLT*Mssp[n]/DJY/DAYSEC#shapiro delay amplitude
            Z <- -Amp*(c1-c2)/(1-cpsi)/ROL
        }
        Zall <- Zall+Z
        Zlist[[n]] <- Z
    }
    list(Zall=Zall,Zlist=Zlist)
}

rv_FullModel <- function(OutObs,OutTime,Par,component='T'){
####################################
## Full model of redshift/rv due to all effects
##
## Input:
##   OutObs - Output of time_Utc2tb
##   OutTime - Output of time_Ta2te
##   Par - Input parameters
##   component - component in the target system
##
## Output:
##   all - All RV effects
##   RvgsO - RV due to general and special relativity in the Solar System
##   RvgT - RV due to general relativity in the target system
##   RvlO - RV due to lensing by Solar System bodies
##   RvSO - RV due to motion of the SSB relative to the observer
##   RvsT - RV due to special relativity in the target system
##   RvST - RV due to the motion of the SSB relative to the target
##   RvSB - RV due to the motion of the SSB relative to the TSB
##   RvBT - RV due to the motion of the TSB relative to the target
##   RvSG - RV due to the motion of the SSB relative to the geocenter
##   RvGO - RV due to the motion of the geocenter relative to the observer
##   RvlT - RV due to lensing by target system bodies (Only companion lensing is considered assuming the lensing effect due to the target star is symmetric)
##   RvTropo - RV due to atmospheric/tropospheric refraction
##   RvRemote - Combined RV due to remote (or target system) effects
##   RvLocal - Combined RV due to local (or Solar System) effects
##   Zcomb - A list of all types of Doppler shifts (the corresponding meanings of names are simliar to that of rv terms)
####################################
    if(component!='T' & component!='C') stop('Error: the component argument is not T or C!')
    RT <- NA
    emrat <- OutObs$emrat
    dt <- ((OutTime$tS[,1]-Par$tpos)+OutTime$tS[,2])/DJY
    Nepoch <- length(dt)
    vSO <- OutObs$SO[,4:6]/auyr2kms#au/yr
    vSB <- OutTime$vSB
    vSG <- OutObs$SG[,4:6]/auyr2kms
    vGO <- OutObs$GO[,4:6]/auyr2kms
    Nt <- nrow(vGO)

    if(component=='T'){
        MT <- Par$mT
        MC <- Par$mC1
        if(any(names(Par)=='RT')) RT <- Par$RT#Target radii (Rsun)
        rST <- OutTime$rST
        rOT <- OutTime$rOT
        rOC <- OutTime$rOC
        ROT <- gen_CalLen(rOT)#pc
        uOT <- OutTime$uOT
        vST <- OutTime$vST#au/yr
        vSC <- OutTime$vSC
        if(Par$binary & Par$Np>0) rTC <- OutTime$rTC
        vBT <- OutTime$vBT#au/yr
    }else{
###exchange T and C if component=C
        MC <- Par$mT
        MT <- Par$mC1
        if(any(names(Par)=='RC')) RT <- Par$RC#Target radii (Rsun)
        rST <- OutTime$rSC
        rOT <- OutTime$rOC
        rOC <- OutTime$rOT
        ROT <- gen_CalLen(rOT)#pc
        uOT <- OutTime$uOC
        vST <- OutTime$vSC#au/yr
        vSC <- OutTime$vST#au/yr
        if(Par$binary & Par$Np>0) rTC <- -OutTime$rTC
        vBT <- OutTime$vSC-vSB
    }
    Ztropo <- OutTime$Ztropo
    dtt.dtcb <- 1/OutObs$dTCB.dTT
    VSO <- gen_CalLen(vSO)
    ZsO <- (VSO/Cauyr)^2/2#Doppler shift due to special relativity
    ZgsO.de <- -(dtt.dtcb-1)#positive for solar system; Gravitational Doppler shift in the solar system at the observer's location due to general (potential) and special relativistic (velocity) effect
    ZgSS <- rv_GravSolar(OutTime$OL)
    ZgO <- ZgSS$Zall
#    ZgsO <- ZgO+ZsO
    ZgsO <- ZgsO.de
    if(Par$Lensing){
        lensing <- rv_LenSolar(OutTime$OL,rOT,vST,vSO,g=1,LenRVmethod=Par$LenRVmethod)
        ZlO <- lensing$Zall
    }else{
        ZlO <- rep(0,Nt)
    }
    Zgc <- rep(0,Nt)#gravitational and convection Doppler shift of the target star
    if(component=='T'){
        star <- Par$star
        vn <- paste0('vOff.',star)
        if(any(names(Par)==vn)) Zgc <- Par[[vn]]/CMPS
    }else{
        star <- Par$stars[Par$stars!=Par$star]
        vn <- paste0('vOff.',star)
        if(any(names(Par)==vn)) Zgc <- Par[[vn]]/CMPS
    }
    ZgT <- Zgc
    if(Par$binary & Par$Np>0){
        if(Par$Einstein){
            RCT <- gen_CalLen(rTC)
            ZgT <- ZgT+rv_GravTarget(MC,RCT,g=Par$g)
            if(!is.na(RT)){
                Zself <- rv_GravTarget(MT,RT*Rsun.au,g=Par$g)#gravitational Doppler shift due to the target if the radius and mass of the target star is known
                                        #	    cat('Zself*CMPS=',Zself*CMPS,'\n')
                ZgT <- ZgT+Zself
            }
            if(length(ZgT)==1){
                ZgT <- rep(ZgT,nrow(uOT))
            }
            if(all(is.na(ZgT)))  ZgT <- rep(0,Nt)
        }else{
            ZgT <- rep(0,Nt)
        }
###lensing in the solar system at the observer's location; assuming the Sun is static in the SSB
###    ZLO <- lensing.shift(1,vLT=vSO,uOT=uOT,rOL=rOS)
###        ZLT <- lensing.shift(ML=MC,vLT=vCT,uOT=uOT,rOL=rOC*pc2au)#lensing in the target system; correspond to shapiro delay
        if(Par$Lensing){
            ZlT <- rv_LenTarget(ML=MC,vSL=vSC,RLT=RCT,rOT=rOT*pc2au,vSO=vSO,rOL=rOC*pc2au,vST=vST,g=1)
        }else{
            ZlT <- rep(0,Nt)
        }
    }else{
        ZlT <- ZgT <- rep(0,Nt)
    }
#the following is for test
#    if(TRUE) uOT <- OutTime$uOB
    RvST <- rowSums(uOT*vST)
    RvSB <- rowSums(uOT*vSB)
    RvBT <- rowSums(uOT*vBT)
    RvSG <- rowSums(uOT*vSG)
    RvGO <- rowSums(uOT*vGO)
    RvSO <- rowSums(uOT*vSO)

    ZST <- RvST/Cauyr#kinematic effects corresponding to Roemer delay
    ZSO <- RvSO/Cauyr
    VST <- sqrt(rowSums(vST^2))
    if(Par$Einstein){
        ZsT <- (VST/Cauyr)^2/2#Doppler shift due to special relativity
    }else{
        ZsT <- rep(0,Nt)
    }
    ZgsT <- ZgT+ZsT

####local and remote
    Zlocal <- (1-ZgsO)/(1+ZSO-ZlO-Ztropo)-1
    Zremote <- (1+ZST-ZlT)/(1-ZgsT)-1
    Z <- (1-ZgsO)/(1-ZgsT)*(1+ZST-ZlT)/(1+ZSO-ZlO-Ztropo)-1
###barycentric correction
    u0 <- t(replicate(Nepoch,Par$pqu[,3]))
    RvST0 <- rowSums(vST*u0)
    ZST0 <- RvST0/Cauyr
    ZLT <- (Par$rv/CKMPS)*(1000/Par$plx)*pc2au/Cauyr*(Par$pmra^2+Par$pmdec^2)*DMAS2R*DMAS2R*dt#light travel term
    if(!Par$binary){
##remove stellar motion and solar system effects
        ZB <- (1+ZSO-ZlO-Ztropo)/(1-ZgsO)*(1+ZST0)/(1+ZST)-1
    }else{
##use Z as the barycentric correction, remove all effects and only planetary perturbation is left
        ZB <- -Z
    }
    ZBwe <- (1+ZSO)/(1-ZgsO)*(1+ZST0)/(1+ZST)-1-ZlO-ZLT#Wright & Eastman 2014 version
    Zcomb <- list(Ztot=Z,Zlocal=Zlocal,Zremote=Zremote,ZgsT=ZgsT,ZgsO=ZgsO,ZST=ZST,ZST0=ZST0,ZlT=ZlT,ZgT=ZgT,ZsT=ZsT,ZSO=ZSO,ZlO=ZlO,Ztropo=Ztropo,ZgsO.de=ZgsO.de,ZgO=ZgO,ZsO=ZsO,ZB=ZB,ZBwe=ZBwe,ZlT=ZlT,ZST0=ZST0,ZgSS=ZgSS)
    if(Par$Lensing) Zcomb <- c(Zcomb,list(Zlensing=lensing$Zlist))

    RvTropo <- Ztropo*CMPS
    RvAll <- Z*CMPS#m/s; absolute RV
    RvRemote <- Zremote*CMPS
    RvLocal <- Zlocal*CMPS
    RvgsO <- ZgsO*CMPS#m/s
    RvgT <- ZgT*CMPS#
###    RvVO <- ZVO*CMPS#
    RvsT <- ZsT*CMPS#
    RvlO <- ZlO*CMPS#
    RvlT <- ZlT*CMPS#
    list(RvTot=RvAll,RvgsO=RvgsO,RvgT=RvgT,RvlO=RvlO,RvSO=RvSO*auyr2kms*1e3,RvsT=RvsT,RvST=RvST*auyr2kms*1e3,RvSB=RvSB*auyr2kms*1e3,RvBT=RvBT*auyr2kms*1e3,RvSG=RvSG*auyr2kms*1e3,RvGO=RvGO*auyr2kms*1e3,RvlT=RvlT,RvTropo=RvTropo,RvRemote=RvRemote,RvLocal=RvLocal,Zcomb=Zcomb)
}

rv_Numerical <- function(utc,Par,OutObs,OutTime){
####################################
## Full model of redshift/rv due to all effects; There are two approaches for numerical computation of RVs: calculate RVs for epochs+/-Tstep separatedly from for assigned epochs or calculate them simultaneously by combining the epochs and epochs+/-Tstep.
##
## Input:
##  utc - JD[utc]
##  Par - Input parameters
##
## Output:
##   all - All RV effects
##   RvgsO - RV due to general and special relativity in the Solar System
##   RvgT - RV due to general relativity in the target system
##   RvlO - RV due to lensing by Solar System bodies
##   RvSO - RV due to motion of the SSB relative to the observer
##   RvsT - RV due to special relativity in the target system
##   RvST - RV due to the motion of the SSB relative to the target
##   RvSB - RV due to the motion of the SSB relative to the TSB
##   RvBT - RV due to the motion of the TSB relative to the target
##   RvSG - RV due to the motion of the SSB relative to the geocenter
##   RvGO - RV due to the motion of the geocenter relative to the observer
##   RvlT - RV due to lensing by target system bodies (Only companion lensing is considered assuming the lensing effect due to the target star is symmetric)
##   RvTropo - RV due to atmospheric/tropospheric refraction
##   RvRemote - Combined RV due to remote (or target system) effects
##   RvLocal - Combined RV due to local (or Solar System) effects
##   Zcomb - A list of all types of Doppler shifts (the corresponding meanings of names are simliar to that of rv terms)
####################################
    Tstep <- Par$Tstep
    OutObs1 <- time_Utc2tb(cbind(utc[,1],utc[,2]+Tstep),Par)
    OutTime1 <- time_Ta2te(OutObs1,Par)
    OutObs0 <- time_Utc2tb(cbind(utc[,1],utc[,2]-Tstep),Par)
    OutTime0 <- time_Ta2te(OutObs0,Par)
    dTAI <- time_T2mT2(OutObs1$JDtai,OutObs0$JDtai)
    dTDB <- time_T2mT2(OutObs1$JDtai,OutObs0$JDtdb)
    dTCB <- time_T2mT2(OutObs1$JDtai,OutObs0$JDtcb)
    dBJDtdb <- time_T2mT2(OutObs1$BJDtdb,OutObs0$BJDtdb)
    dBJDtcb <- time_T2mT2(OutObs1$BJDtcb,OutObs0$BJDtcb)
    dTT <- time_T2mT2(OutObs1$JDtai,OutObs0$JDtt)
    dTCG <- time_T2mT2(OutObs1$JDtai,OutObs0$JDtcg)
    dUT1 <- time_T2mT2(OutObs1$JDut1,OutObs0$JDut1)
    dtB <-time_T2mT2(OutTime1$tB,OutTime0$tB)
    dtauE <-time_T2mT2(OutTime1$tauE,OutTime0$tauE)
    dRoemerTarget <- OutTime1$RoemerTarget-OutTime0$RoemerTarget
    dRoemerSolar <- OutTime1$RoemerSolar-OutTime0$RoemerSolar
    dEinsteinSolar <- OutTime1$EinsteinSolar-OutTime0$EinsteinSolar
    dEinsteinIS <- OutTime1$EinsteinIS-OutTime0$EinsteinIS
    dEinsteinTarget <- OutTime1$EinsteinTarget-OutTime0$EinsteinTarget
    dShapiroSolar <- OutTime1$ShapiroSolar-OutTime0$ShapiroSolar
    dShapiroTarget <- OutTime1$ShapiroTarget-OutTime0$ShapiroTarget
#    ZgsO <- -(OutObs$dTCB.dTT-1)
    ZgsO <- 1/OutObs$dTCB.dTT-1
    ZlS <- ZlO <- dShapiroSolar/dTCB
    ZkS <- ZSO <- rowSums(OutTime$uOT*OutObs$SO[,4:6])/CKMPS
    ZSG <- rowSums(OutTime$uOT*OutObs$SG[,4:6])/CKMPS
    ZGO <- rowSums(OutTime$uOT*OutObs$GO[,4:6])/CKMPS
    ZgsT <- dEinsteinTarget/dtB+gen_CalLen(OutTime$vST)^2/2/Cauyr^2
    ZlT <- dShapiroTarget/dtB
    ZkT <- ZST <- rowSums(OutTime$uOT*OutTime$vST)/Cauyr
    Ztropo <- OutTime$Ztropo
    Zremote <- (1+ZST-ZlT)/(1-ZgsT)-1
    Zlocal <- (1-ZgsO)/(1+ZSO-ZlS-Ztropo)-1
    Ztot <- (1-ZgsO)/(1-ZgsT)*(1+ZST-ZlT)/(1+ZSO-ZlO-Ztropo)-1
    RvTot <- Ztot*CMPS
    Zcomb <- list(Ztot=Ztot,Zlocal=Zlocal,Zremote=Zremote,ZgsT=ZgsT,ZgsO=ZgsO,ZST=ZST,ZlT=ZlT,ZSO=ZSO,ZlO=ZlO,Ztropo=Ztropo,ZlT=ZlT)
    list(RvTot=RvTot,RvgsO=ZgsO*CMPS,RvgsT=ZgsT*CMPS,RvST=ZST*CMPS,RvSG=ZSG*CMPS,RvGO=ZGO*CMPS,RvlO=ZlO*CMPS,RvSO=ZSO*CMPS,RvlT=ZlT*CMPS,RvTropo=Ztropo*CMPS,RvRemote=Zremote*CMPS,RvLocal=Zlocal*CMPS,Zcomb=Zcomb)
}
