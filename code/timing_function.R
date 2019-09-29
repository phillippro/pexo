time_CalHms2JD <- function(CalHms){
####################################
## Convet a calendar time in format of Year, Month, Day, Hour, Minute, Second to Julian Date
##
## Input:
##   CalHms - Time in the format of year, month, day, hour, miniute, and second
##
## Output:
##   JD - 2-part Julian Date
####################################
    cal2 <- calHms[,1:3,drop=FALSE]
    jd2 <- time_Cal2JD(cal2)
    fd <- cal[,4]/24+cal[,5]/(24*60)+cal[,6]/(24*3600)
    cbind(jd2[,1],jd2[,2]+fd)
}

time_Cal2JD <- function(cal){
####################################
## Gregorian Calendar date to Julian Date.
## ref: https://uk.mathworks.com/matlabcentral/fileexchange/15285-geodetic-toolbox?focused=6513448&tab=function
## CAL2JD  Converts calendar date to Julian date using algorithm
## from "Practical Ephemeris Calculations" by Oliver Montenbruck
##   (Springer-Verlag, 1989). Uses astronomical year for B.C. dates
##   (2 BC = -1 yr). Non-vectorized version. See also DOY2JD, GPS2JD,
##   JD2CAL, JD2DOW, JD2DOY, JD2GPS, JD2YR, YR2JD.
##
## Input:
##   cal - Calendar date with day fraction
##
## Output:
##   JD - 2-part or 2D Julian Date
####################################
    y <- yr <- cal[,1]
    m <- mn <- cal[,2]
    dy <- cal[,3]
    ind <- which(mn<=2)
    y[ind] <- y[ind]-1
    m[ind] <- m[ind]+12
    date1 <- 4.5+31*(10+12*1582)# Last day of Julian calendar (1582.10.04 Noon)
    date2 <- 15.5+31*(10+12*1582)# First day of Gregorian calendar (1582.10.15 Noon)
    date <- dy+31*(mn+12*yr)
    ind1 <- which(date<=date1)
    ind2 <- which(date>=date2)
    b <- y
    b[ind1] <- -2
    b[ind2] <- trunc(y/400) - trunc(y/100)
    if(length(ind1)==0 & length(ind2)==0){
        cat('Dates between October 5 & 15, 1582 do not exist!\n')
    }
    ind1 <- which(y>0)
    ind2 <- which(y<0)
    jd <- y
    jd[ind1] <- trunc(365.25*y[ind1]) + trunc(30.6001*(m[ind1]+1)) + b[ind1] + 1720996.5 + dy[ind1]
    jd[ind2] <- trunc(365.25*y[ind2]-0.75) + trunc(30.6001*(m[ind2]+1)) + b[ind2] + 1720996.5 + dy[ind2]
#    return(cbind(DJM0,jd-DJM0))
    return(cbind(jd%/%1,jd%%1))
}

time_T2mMjd <- function(t1,t2){
####################################
## The difference between a 2D time and a 1D MJD
## A 2-part time is composed by a base and a fraction
## The sum of a 2-part time is JD
##
## Input:
##   t1 -  2-part time
##   t2 -  1-part MJD
##
## Output:
##   Difference between t1 and t2 in units of second
####################################
    ((t1[,1]%%DJM0-t2)+t1[,2])*DAYSEC
}

time_T2mT2 <- function(t1,t2,toSecond=TRUE){
####################################
## The difference between two 2D times
## A 2D time is composed by a base and a fraction
## The sum of a 2D time is JD
##
## Input:
##   t1 -  first 2D time
##   t2 -  second 2D time
##
## Output:
##   Difference between t1 and t2 in units of second
####################################
    dt <- ((t1[,1]-t2[,1])+(t1[,2]-t2[,2]))
    if(toSecond){
        dt*DAYSEC
    }else{
        dt
    }
}
time_ChangeBase <- function(t,denom=1){
    t1 <- (rowSums(t)%/%denom)*denom
    t2 <- (t[,1,drop=FALSE]-t1)+t[,2,drop=FALSE]
    return(cbind(t1,t2))
}
time_Df2hms <- function(f){
####################################
## Day fraction to hour, minute and second
##
## Input:
##   f  Fraction of a day
##
## Output:
##   h hour
##   m minute
##   s second
####################################
    h <- floor(f*24)
    m <- floor(f*24*60-h*60)
    s <- f*24*3600-h*3600-m*60
    return(c(h,m,s))
}
time_Doy2jd <- function(yr,doy){
####################################
## Convert integer Year and nubmer of days since the begining of the year to 2-part Julian Date
## ref: https://www.mathworks.com/matlabcentral/fileexchange/15285-geodetic-toolbox
##
## Input:
##   yr - Integer Year
##   doy - Number of days since the begining of the year
##
## Output:
##   jd - 2-part Julian Date
####################################
    jd <- time_Cal2JD(cbind(yr,1,0))
    cbind(jd[,1]+doy,jd[,2])
}
time_Jd2yr <- function(jd){
####################################
## Convert Julian Date to Year
## ref: https://www.mathworks.com/matlabcentral/fileexchange/15285-geodetic-toolbox
##
## Input:
##   jd - Julian Date
##
## Output:
##   yr - Year
####################################
    cal <- time_Jd2cal(jd)
    jd0  <-  time_Cal2JD(cbind(cal[,1],1,1))
    jd1  <-  time_Cal2JD(cbind(cal[,1]+1,1,1))
    if(!is.null(jd0)){
        return(cal[,1] + rowSums(jd-jd0)/(rowSums(jd1-jd0)))
    }else{
        return(cal[,1] + (jd-jd0)/(jd1-jd0))
    }
}

time_Yr2jd <- function(yr){
####################################
## Convert Year to Julian Date
## ref: https://www.mathworks.com/matlabcentral/fileexchange/15285-geodetic-toolbox
##
## Input:
##   yr -  Year
##
## Output:
##   jd - Julian Date
####################################
    iyr <-  floor(yr)
    jd0  <-  time_Cal2JD(cbind(iyr,1,1))
    days  <-  rowSums(time_Cal2JD(cbind(iyr+1,1,1)) - jd0)
    doy  <-  (yr-iyr)*days + 1
    time_Doy2jd(iyr,doy)
}
time_Ut1mUtc <- function(mjd,tempo=TRUE){
####################################
## Given a epoch in MJD[UTC], calculate UT1-UTC.
## If TEMPO2 method is used, dut for early epochs is set to be the first entry of the EOP file.
## Duet for later epochs is fixed at the last entry of the EOP file.
## If TEMPO2 method is not used, dut=0 for epochs earlier or later than the EOP file start and end dates, respectively.
##
## Input:
##   mjd -  Modified Julian Date
##   tempo -  Whether the TEMPO2 method is used or not
##
## Output:
##   dut1 - UT1-UTC
####################################
    dut1 <- rep(NA,length(mjd))
    mjd.max <- max(eop[,'MJD'])
    ind1 <- which(mjd>mjd.max)
    if(Par$CompareT2){
        dut1[ind1] <- eop[nrow(eop),'UT1-UTC']
    }else{
        dut1[ind1] <- pred.dut1(mjd[ind1])
    }
    mjd.min <- min(eop[,'MJD'])
    ind2 <- which(mjd<mjd.min)
    if(Par$CompareT2){
        dut1[ind2] <- eop[1,'UT1-UTC']
    }else{
        dut1[ind2] <- 0
    }
    inds <- which(mjd>=mjd.min & mjd<=mjd.max)
    if(length(inds)>0){
        ind3 <- sapply(mjd[inds],function(x) which(x<eop[,'MJD'])[1]-1)
        dut1[inds] <- eop[ind3,'UT1-UTC']
    }
    dut1
}

time_Utc2tt <- function(utc,Par){
####################################
## Convert UTC to TT
##
## Input:
##   utc -  2-part UTC time
##   Par - Input parameters
##
## Output:
##   tt - 2-part TT time
##   tai - 2-part TAI time
##   leap - Leap seconds
####################################
    tmp <- sofa_Utctai(utc,Par$TaiType)
    tai <- tmp$tai
    tt <- sofa_Taitt(utc,tai,Par=Par)
    list(tt=tt,tai=tmp$tai,leap=tmp$leap)
}

time_GetEop <- function(utc,tai){
####################################
## Get Earth orientation parameters (EOP)
##
## Input:
##   utc -  2-part UTC time
##
## Output:
##   ut1 - UT1 time
##   tmp - a list of variables including
####################################
    mjd.utc <- (utc[,1]-DJM0)+utc[,2]
    mjd.max <- max(eop[,'MJD'])
    mjd.min <- min(eop[,'MJD'])
    if(any(mjd.utc<mjd.min)) cat('Warning: some epochs are earlier than the minium epoch of EOP file which is MJD',mjd.min,'!\n') #return(NULL)
    if(any(mjd.utc>mjd.max)) cat('Warning: some epochs are later than the maximum epoch of EOP file which is MJD',mjd.max,'!\n') #return(NULL)
    tmp <- time_Dut1dot(mjd.utc,TAImUTC=time_T2mT2(tai,utc))
#    dut1dot <- tmp$dut1dot
#    leap <- tmp$leap
    ut1 <- cbind(utc[,1],utc[,2]+tmp$dut1/DAYSEC)
    ut1 <- time_ChangeBase(ut1,1)
    c(list(ut1=ut1),tmp)
}

time_Itrs2icrsApp <- function(tt,utc,rpom,EopPar){
##Approximate procedures; from Wallance 2006
    Nt <- nrow(tt)
    tau <- (EopPar$ut1[,1]-DJ00)+EopPar$ut1[,2]
    era <- sofa_Era00(EopPar$ut1)
    Omega <- 2.182-9.242e-4*tau
    X <- 2.6603e-7*tau-33.2e-6*sin(Omega)
    Y <- -8.14e-14*tau^2+44.6e-6*cos(Omega)
    Mt2c <- Mc2t <- array(NA,dim=c(3,3,Nt))
    for(j in 1:nrow(tt)){
        Mcio <- rbind(c(1,0,X[j]),c(0,1,-Y[j]),c(X[j],Y[j],1))
        R <- R3(era[j])%*%Mcio
        val <- rpom[,,j]%*%R
        Mc2t[,,j] <- val
        Mt2c[,,j] <- t(val)
    }
    list(Mc2t=Mc2t,Mt2c=Mt2c)
}


time_Dut1dot <- function(mjd,TAImUTC,T2=TRUE){
####################################
## Calculate the time derivative of UT1-UTC
##
## Input:
##   mjd - MJD in UTC time standard
##   TAImUTC - TAI-UTC (s); for epochs earlier than 1962
##   T2 - whether or not use TEMPO2 approach; default is TRUE
##
## Output:
##   dut1dot - d(UT1-UTC)/dUTC
##   leap - Leap seconds
##   xp - xp parameter in units of rad in Earth rotation model
##   yp - yp parameter in units of rad in Earth rotation model
##   dut1 - UT1-UTC
####################################
##deal with earlier epochs
t1 <- proc.time()
    eop1 <- eop[,c('MJD','UT1-UTC','x','y')]
    mjd.max <- max(eop[,'MJD'])
    mjd.min <- min(eop[,'MJD'])
    mjd1900 <- time_Jd2mjd(time_Yr2jd(1900))
    indE <- which(mjd<mjd.min & mjd>=mjd1900)
    if(length(indE)>0 & exists('eop1900') ){
            UT1mUTC <- UT1mTAIfunc(mjd[indE])+TAImUTC[indE]
            eop1 <- rbind(as.matrix(cbind(mjd[indE],UT1mUTC,xfunc(mjd[indE]),yfunc(mjd[indE]))),as.matrix(eop1))
            colnames(eop1) <- c('MJD','UT1-UTC','x','y')
    }
    mjd.min <- min(eop1[,'MJD'])
    N <- nrow(eop1)
    ind.early <- which(mjd<mjd.min)
    ind.late <- which(mjd>mjd.max)
    if(T2){
        mjd[ind.early] <- mjd.min
        mjd[ind.late] <- mjd.max
        isamp <- sapply(mjd,function(x) which(x<eop1[,'MJD'])[1])
#        isamp <- match(x,eop1[,'MJD'])
        isamp[is.na(isamp)] <- nrow(eop1)
    }else{
        isamp <- sapply(mjd,function(x) which(x<eop1[,'MJD'])[1])

        mjd[is.na(isamp)] <- eop1[nrow(eop1),'MJD']
        isamp[is.na(isamp)] <- nrow(eop1)

        mjd[isamp==0] <- eop1[1,'MJD']
        isamp[isamp==0] <- 1
    }
    dut1s <- eop1[,'UT1-UTC']
#cope with leap second ... take it off second point as jump happens right AT second point
    dut11 <- dut1s[isamp]
    ind2 <-isamp-1
    dut10 <- dut1s[ind2]

    leap <- rep(0,length(mjd))
    leap[(dut11 - dut10) > 0.5] <- 1
    leap[(dut11 - dut10) < -0.5] <- -1
    dut1dot <- (dut11 - leap - dut10)/DAYSEC

    ##interpolate
    f <- (mjd - eop1[ind2,'MJD']) / (eop1[isamp,'MJD'] - eop1[ind2,'MJD'])
    xp <- eop1[ind2,'x'] + f*(eop1[isamp,'x'] - eop1[ind2,'x'])
    yp <- eop1[ind2,'y'] + f*(eop1[isamp,'y'] - eop1[ind2,'y'])
    dut1 <- dut10 + f*(dut11 - leap - dut10)
    return(list(dut1dot=dut1dot,leap=leap,xp=xp/DR2AS,yp=yp/DR2AS,dut1=dut1))
}

time_Itrs2icrs06 <- function(tt,utc,t,rpom,EopPar){
####################################
## Convert the ITRS to ICRS reference frame using the IAU2006 model.
## ref: Wallace & Capitaine 2006
## Input:
##   tt - TT
##   utc - UTC
##   t - Julian century since J2000.0
##   rpom - Polar-motion matrix
##   EopPar - Earth orientation parameters
##      dut1dot - Time derivative of dut1
##      leap - Leap seconds
##      xp - EOP parameter (rad)
##      yp - EOP parameter (rad)
##      dut1 - UT1-UTC
##
## Output:
##   Mc2t - rotation matrix from ICRS to ITRS
####################################
    Nt <- nrow(tt)
###Earth rotation angle
    era <- sofa_Era00(EopPar$ut1)
####all angles are in arcsec
####precision and nutation
    xi0 <- -0.016617
    eta0 <- -0.006819
    dalpha0 <- -0.0146
    eps0 <- 84381.406000000#P03
    psiA <- 5038.481507*t-1.0790069*t^2-1.14045e-3*t^3+0.132851e-3*t^4-0.0000951e-3*t^5
    omegaA <- eps0-25.754e-3*t+51.2623e-3*t^2-7.72503e-3*t^3-0.000467e-3*t^4+0.0003337e-3*t^5
    epsA <- eps0-46.836769*t-0.1831e-3*t^2+2.00340e-3*t^3-0.000576e-3*t^4-0.0000434e-3*t^5
    chiA <- 10.556403*t-2.3814292*t^2-1.21197e-3*t^3+0.170663e-3*t^4-0.0000560e-3*t^5
    tmp <- sofa_Nut00a(tt)
    dpsi00 <- tmp$dpsi/DAS2R#as
#    cat('dpsi00=',dpsi00,'\n')
    deps00 <- tmp$deps/DAS2R
    f <- -2.7774e-6*t
    dpsi <- dpsi00+(0.4697e-6+f)*dpsi00
    deps <- deps00+f*deps00
###B (bias):eta0,xi0,dalpha0
###P (precession): chiA, omegaA, psiA, epsA
###N (nutation): epsA, deps, dpsi, epsA
    B <- gen_R1(-eta0*DAS2R)%*%gen_R2(xi0*DAS2R)%*%gen_R3(dalpha0*DAS2R)
    Mt2c <- Mc2t <- array(NA,dim=c(3,3,Nt))
#Mcios
    for(j in 1:Nt){
        N <- gen_R1(-(epsA[j]+deps[j])*DAS2R)%*%gen_R3(-dpsi[j]*DAS2R)%*%gen_R1(epsA[j]*DAS2R)
        P <- gen_R3(chiA[j]*DAS2R)%*%gen_R1(-omegaA[j]*DAS2R)%*%gen_R3(-psiA[j]*DAS2R)%*%gen_R1(eps0*DAS2R)
        Mclass <- N%*%P%*%B
###calculate s
        x <- Mclass[3,1]
        y <- Mclass[3,2]
        s <- sofa_S06(tt[j,,drop=FALSE],x,y)
###CIO matrix (CIO stand for celestial intermediate pole)
        Mcio <- gen_CalMcio(s,x,y)
###GCRS-to-TIRS matrix
        R <- gen_R3(era[j])%*%Mcio
        val <- rpom[,,j]%*%R
        Mc2t[,,j] <- val
        Mt2c[,,j] <- t(val)
    }
    list(Mc2t=Mc2t,Mt2c=Mt2c)
}

time_BipmCorr <- function(mjd.utc,BIPM='2018'){
####################################
## Correction of TT-TAI
## The BIPM file is downloaded from ftp://ftp2.bipm.org/pub/tai/ttbipm/TTBIPM.17 and could be automatically updated
##  1st column : MJD at 0 h UTC
##  2nd   "    : TT(BIPM17) - EAL - 32.184 s, unit is one microsecond
##  3rd   "    : TT(BIPM17) - TAI - 32.184 s, unit is one microsecond
##
## Input:
##   mjd.utc - MJD in UTC standard
##   BIPM - BIPM version/year
##
## Output:
##  dt - Correction for TT-TAI in units of second
####################################
    mjd.min <- min(BIPMfile[,1])
    mjd.max <- max(BIPMfile[,1])
    Nt <- length(mjd.utc)
    dt <- rep(0,Nt)
    f <- approxfun(BIPMfile[,1],BIPMfile[,3])

    ind.match <- which(mjd.utc>=mjd.min & mjd.utc<=mjd.max)
    ind.future <- which(mjd.utc>mjd.max)
    ind.past <- which(mjd.utc<mjd.min)
    dt[ind.match] <- f(mjd.utc[ind.match])*1e-6
    if(BIPM=='2016'){
        dt[ind.future] <- 27679e-9-0.05*(mjd.utc[ind.future]-57749)*1e-9
    }
    if(BIPM=='2017'){
        dt[ind.future] <- 27661e-9
    }
    if(BIPM=='2018'){
        dt[ind.future] <- (27672.0 + 0.03*(mjd.utc[ind.future]-58479))*1e-9
    }
    if(length(dt)>length(ind.past) & length(ind.past)>0)  dt[ind.past] <- dt[-ind.past][1]
    return(dt)
}

time_TropoDelay <- function(mjd.utc,elevation,delevation,Par){
####################################
## Calculate hydrostatic tropospheric delay
## Ref: Edwards, Hobbs & Manchester 2006
## Default air pressure value is P=101.325 kPa.
##
## Input:
##   mjd.utc - MJD in UTC standard
##   elevation - Elevation angle of the target in radian
##   delevation - Time derivative of the elevation angle (rad/day)
##   Par - Input parameters
##
## Output:
##   dt - Troposheric delay in second
##   z - Redshift due to tropospheric delay
####################################
    pressure <- Par$pmb/10#kpa
    H <- Par$height*1e3#metre
    phi <- Par$phi#rad
###David et al., 1985; Elgered et al. 1991; https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/90JB00834
    dt.hz <-  0.02268*pressure/CMPS/(1-0.00266*cos(2*phi)-2.8e-7*H)
    dry <- astro_NMFhydro(mjd.utc,as.numeric(phi),as.numeric(H),as.numeric(elevation),as.numeric(delevation))
    list(dt=dt.hz*dry$map,z=dt.hz*dry$dmap)
}

time_Geo2obs <- function(utc,tai,tt,ObsPar,EopType='2006',n=1){
####################################
## Convert quantities from ITRS frame to ICRS frame
##
## Input:
##   tt - 2-part TT
##   utc - 2-part UTC
##   ObsPar - Obser
##   EopType - Type of EOP parameters
##   n - ellipsoid type
##   n=1     WGS84 -> default model in PEXO; more advanced than GRS80
##   n=2     GRS80 -> used by tempo2
##   n=3     WGS72
##
## Output:
##   rGO - Vector from geocenter to the observer in the ICRS frame in units of km
##   rGOitrs - Vector from geocenter to the observer in the ITRS frame in units of km
##   vGO - Geocentric velocity of the observer in ICRS frame in units of km/s
##   vGOitrs - Geocentric velocity of the observer in ITRS frame in units of km/s
##   ZenithIcrs - Zenith direction in ICRS frame
##   ZenithItrs - Zenith direction in ITRS frame
##   PoleIcrs - Pole direction in ICRS frame
##   PoleItrs - Pole direction in ITRS frame
##   EopPar - Earth orientation parameters
####################################
###ITRS position
    t <- ((tt[,1] - DJ00)+tt[,2])/DJC#Julian century since J2000.0
    Nt <- nrow(tt)
    robs.itrs <- as.numeric(c(Par$xtel,Par$ytel,Par$ztel))
    zenith.itrs <- c(cos(Par$elong)*cos(Par$phi),sin(Par$elong)*cos(Par$phi),sin(Par$phi))
    ##Get Earth orientation parameters
    EopPar <- time_GetEop(utc,tai)
    ##other ITRS quantities ...
    itrs <- gen_RvzItrs(robs.itrs,EopPar,Nt,sp=NULL)

    ##ITRS to ICRS; three methods
    if(EopType=='2006'){
        Mt2c <- time_Itrs2icrs06(tt,utc,t,itrs$rpom,EopPar)$Mt2c
    }else if(EopType=='2000B'){
        Mt2c <- sofa_C2t00b(tt,EopPar$ut1,EopPar$xp,EopPar$yp,itrs$rpom)$Mt2c
    }else if(EopType=='approx'){
        Mt2c <- time_Itrs2icrsApp(tt,utc,itrs$rpom,EopPar)$Mt2c
    }

    ##Multiply the itrs position/velocity vector by its transpose (=inverse) to transform trs->crs
    pole.icrs <- zenith.icrs <- robs.icrs <- vobs.icrs <- array(NA,dim=c(Par$Nepoch,3))
    for(j in 1:Par$Nepoch){
        robs.icrs[j,] <- Mt2c[,,j]%*%robs.itrs#position in CRS; km
        vobs.icrs[j,] <- Mt2c[,,j]%*%t(itrs$vobs[j,,drop=FALSE])#velocity in CRS; km/s
        zenith.icrs[j,] <- Mt2c[,,j]%*%zenith.itrs#zenith in CRS
        pole.icrs[j,] <- Mt2c[,,j]%*%itrs$pole[j,]#Pole in CRS
    }
    return(c(list(rGO=robs.icrs,rGOitrs=robs.itrs,vGO=vobs.icrs,vGOitrs=t(itrs$vobs[j,,drop=FALSE]),ZenithIcrs=zenith.icrs,ZenithItrs=zenith.itrs,PoleIcrs=pole.icrs,PoleItrs=itrs$pole),EopPar))
}

time_Jd2bjd <- function(OutBary,OutSB,OutEle,Par,CompareT2=FALSE){
####################################
## Convert from JD to BJD considering the Tropospheric delay, Roemer delay, Shapiro delay, Einstein delay;
## the input position vectors are in units of pc and the input velocity vectors are in units of
##
## Input:
##   OutBary - output of time_Utc2tb
##   OutSB - output of astro_CalSB
##   OutEle - output of astro_CalElevation
##   Par - Input parameters
##
## Output:
##   bjd.tcb - BJD in TCB standard
##   bjd.tdb - BJD in TDB standard
##   roemer - Roemer delay in units of second
##   ShapiroSolar - Shapiro delay in units of second
##   RoemerDelay - All types of Roemer delays in units of second
##   TropoDelay - Tropospheric delay in second
##   Ztropo - Tropospheric redshift
##   RoemerT2 - Roemer delay using the TEMPO2 method (equivalent to degraded PEXO called PEXOt)
####################################
###troposphere delay
    mjd.utc <- time_Jd2mjd(utc)
    if(Par$ObsType=='ground'){
        tropo <- time_TropoDelay(mjd.utc,OutEle$elevation,OutEle$delevation,Par)
        TropoDelay <- tropo$dt*OutBary$dTCB.dTT
        Ztropo <- tropo$z
    }else{
        TropoDelay <- Ztropo <- 0
    }
###Roemer delay in solar system
    roemer <- time_RoemerSolar(OutBary,OutSB,SB=FALSE)#day; in coordiante units
    if(CompareT2){
        roemerSB <- time_RoemerSolar(OutBary,OutSB,SB=TRUE)#day; in coordiante units
    }else{
        roemerSB <- NULL
    }
    DRS <- roemer$Roemer
###Shapiro delay in solar system
    ShapiroSolar <- time_ShapiroSolar(OutBary,OutSB,Par)
    DSS <- ShapiroSolar$dt.all#second
    DDT <- (DRS+DSS+TropoDelay)/DAYSEC#Sum of Roemer, Shapiroa and tropospheric delays in the solar system
  if(Par$Unit=='TCB'){
        bjd.tcb <- cbind(OutBary$JDtcb[,1],OutBary$JDtcb[,2]-DDT)
        bjd.tdb <- sofa_Tcbtdb(bjd.tcb)
    }else{
        bjd.tdb <- cbind(OutBary$JDtdb[,1],OutBary$JDtdb[,2]-DDT/IFTE.K)
        bjd.tcb <- sofa_Tdbtcb(bjd.tdb)
    }
    return(list(BJDtcb=bjd.tcb,BJDtdb=bjd.tdb,Roemer=DRS,ShapiroSolar=ShapiroSolar,RoemerOrder=list(Roemer1=roemer$Roemer1,Roemer2=roemer$Roemer2,Roemer3=roemer$Roemer3),TropoDelay=TropoDelay,Ztropo=Ztropo,RoemerT2=roemer$RoemerT2,RoemerSB=roemerSB))
}


time_TsTb1 <- function(tO,OS,Par){
####################################
## Estimate the first order delay with respect to the reference epoch tS=tB=tpos
##
## Input:
##   tO - Observer coordinate time, equvalent to TCB
##   OS - Position (km) and (km/s) velocity vectors from the observer to SSB
##   Par - Input parameters
##
## Output:
##   tSf - first order coordinate time at SSB
##   tBf - first order coordinate time at TSB
####################################
    tol <- 1e-10
    ra <- Par$ra
    dec <- Par$dec
    pmra <- Par$pmra#mas/yr
    pmdec <- Par$pmdec#mas/yr
    plx <- Par$plx#mas
    rv <- Par$rv
    dOS0 <- OS[,1]*Par$u[1]+OS[,2]*Par$u[2]+OS[,3]*Par$u[3]
    dt.rom <- dOS0*AULT/DAYSEC#day
    v2 <- (pmra^2+pmdec^2)/plx^2+(rv/auyr2kms)^2#au/yr
    v2.c2 <- v2/Cauyr^2
###initial guess of tS; only consider fractional time
    tSf  <- tO[,2]-dt.rom#fraction
    if(Par$SBscaling){
        beta <- rv/CKMPS
        dt <- (tO[,1]-Par$tpos)+tSf#tS-tpos
        dt.vp <- beta*dt#day
        dt.es <-2*beta^2
        tBf <- tSf-dt.vp-0.5*v2.c2*(dt-dt.vp)#day
    }else{
        tBf <- tSf
    }
    cbind(tSf,tBf)
}

time_Jd2mjd <- function(jd){
####################################
## Convert JD to MJD
##
## Input:
##   jd - 2-part Julian Date
##
## Output:
##   mjd - Modified Julian Date
####################################
    if(is.null(dim(jd))){
        jd <- t(jd)
    }
    (jd[,1]-DJM0)+jd[,2]
}

time_Tt2tbObs <- function(tcg,tt,TDBmTTgeo,TDBmTTgeo0,TDBmTTobs,Zgeo){
####################################
## Convert TT to TCB and TDB at the observer's location given multiple quantities derived from some TT-TDB methods (see time_Tt2tb for details)
## all times are in units of day
## ref: eq. 10.5 of https://www.iers.org/SharedDocs/Publikationen/EN/IERS/Publications/tn/TechnNote36/tn36_151.pdf?__blob=publicationFile&v=1
##
## Input:
##   tt - 2-part TT
##   TDBmTTgeo - TDB-TT
##   TDBmTTgeo0 - TDB-TT at the reference epoch JD[TT0]=2443144.5003725
##   TDBmTTobs - TDB-TT due to geocenter-observation offset
##   Zgeo - d(TDB-TT)/dTDB at the geocenter
##
## Output:
##   tcb - 2-part TCB
##   tdb - 2-part TDB
####################################
    dt <- (IFTE.LC*(rowSums(tt)-TT0)+TDBmTTgeo-TDBmTTgeo0)/(1-ELB)+TDBmTTobs
    tcb <- time_ChangeBase(cbind(tcg[,1],tcg[,2] + dt),1)
    tdb <- sofa_Tcbtdb(tcb)
    list(tcb=tcb,tdb=tdb)
}

time_Tt2tb <- function(utc,tai,tcg,tt,Par){
####################################
## Convert TT to TDB at the observer's location with multiple methods.
## Calculate TDB accounting for the relativistic effects at the observer's location.
## According to my tests, ephemeris and FB01 are the two most accurate methods.
## The returned quantities are in units compatible with TDB.
##
## Input:
##   utc - 2-part UTC
##   tai - 2-part TAI
##   tcg - 2-part TCG
##   tt - 2-part TT
##   Par - Input parameters defined in read_input.R
##
## Output:
##   tdb - 2-part TDB
##   tcb - 2-part TCB
##   geo - output of gen_GeoEph()
##   dTT.dTDB - dTT/dTDB
##   dTDB.dTT- dTDB/dTT
##   dTCB.dTT - dTCB/dTT
##   dzenith - dzenith/dt
##   emrat - Earth-Moon mass ratio
##   GM - Position and velocity vectors from the geocenter to the Moon
##   zTDBmTTgeo - d(TDB-TT)/dTDB at the geocenter
##   zTDBmTTobs - d(TDB-TT)/dTDB due to relativistic effect of the observer in the geocentric frame
##   zTDBmTTobsR - zobs due to geocenter-observer offset
##   zTDBmTTobsV - zobs due to Earth's rotation
##   TDBmTTgeo - TDB-TT at geocenter
##   TDBmTTobs - TDB-TT due to relativistic effect of the observer in the geocentric frame
##   OutGeo2obs - Output of time_Geo2obs()
####################################
    Ntry <- 10
    h <- tol <- 1e-16#day
##
    mjd.tt <- (tt[,1]-DJM0)+tt[,2]
#    t <- ((tt[,1] - DJ00)+tt[,2])/DJC#Julian century since J2000.0

    ##BCRS state of geocenter; solar system barycenter (SST) reference frame

    ##calculate initial observation term; geocentric reference frame

#    obs2 <- time_Geo2obs(cbind(tt[,1],tt[,2]+dt),cbind(utc[,1],utc[,2]+dt),ObsPar,EopType=EopType)

##Given TT and TDB(geo)-TT as a function of TDB(geo), TDB(geo) can be estimated through iterations
##ref: DE430 document by Folkner et al. 2014
##https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430_and_de431.pdf
##The following loop calculates TDB-TT at the geocenter using JPL ephemeris
    tt0 <- c(TT0%/%1,TT0%%1)
    jd <- rbind(tt0,tt)
    TDBmTTgeo1 <- 0
    TtTdbMethod <- Par$TtTdbMethod
    for(j in 1:Ntry){
##only consider first-order roemer delay from observer to geocenter
        if(TtTdbMethod=='eph' | Par$ObsType=='space'){
            cat('Calculate TDB-TT using JPL ephemeris!\n')
            tt2tdb.geo <- gen_CalEph(jd,body='TT-TDB (at geocenter)',DE=Par$DE)
            TDBmTTgeo <- -tt2tdb.geo[,1]#TDB-TT; s
            zTDBmTTgeo <- -tt2tdb.geo[-1,2]
            TDBmTTgeo0 <- TDBmTTgeo[1]
            TDBmTTgeo <- TDBmTTgeo[-1]
        }else if(TtTdbMethod=='FB01'){
            cat('Calculate Teph-TT as a function of Teph using FB01 method!\n')
            tmin <- min(rowSums(jd))
            tmax <- max(rowSums(jd))
            dt <- median(diff(rowSums(utc)))
            ddt <- 10#min(1,dt/10)
            maxday <- floor((tmax-tmin)/ddt)
            jd.ref <- tmin+(0:maxday)*ddt
            write.table(c(tmin,tmax,ddt,'F'),file='FB01/fb2001.in',quote=FALSE,row.names=FALSE,col.names=FALSE)
            system('cd FB01 ; ./fb2001 <fb2001.in >fb2001.out')
            fb <- read.table('FB01/fb2001.out',skip=3)
            inds <- sapply(1:nrow(jd),function(i) which.min(abs(jd.ref-sum(jd[i,]))))

            TDBmTTgeo <- as.numeric(gsub('D','e',as.character(fb[inds,2])))#day;(Teph-TT)(Teph)
            zTDBmTTgeo1 <- as.numeric(gsub('D','e',as.character(fb[inds,3])))#day
            zTDBmTTgeo <- zTDBmTTgeo1[-1]
            dTDBmTTgeo <- zTDBmTTgeo1*time_T2mT2(jd,time_ChangeBase(cbind(jd.ref[inds],0),1),toSecond=FALSE)#correction by shifting teph.ref to teph
            TDBmTTgeo <- (TDBmTTgeo+dTDBmTTgeo)*DAYSEC+IFTE.TEPH0#s; correct to Teph-TT
            TDBmTTgeo0 <- TDBmTTgeo[1]
            TDBmTTgeo <- TDBmTTgeo[-1]
        }else if(TtTdbMethod=='FBgeo'){
            cat('Calcuate TDB-TT using FBsofa method for geocenter!\n')
            tmp <- sofa_Dtdb(jd[,1],jd[,2],ut=rowSums(jd),elong=0,u=0,v=0)
            tmp1 <- sofa_Dtdb(jd[,1],jd[,2]+0.01,ut=rowSums(jd),elong=0,u=0,v=0)
            tmp2 <- sofa_Dtdb(jd[,1],jd[,2]-0.01,ut=rowSums(jd),elong=0,u=0,v=0)
            zTDBmTTgeo <- -(tmp2[-1,1]-tmp1[-1,1])/(0.02*DAYSEC)
#            zTDBmTTgeo <- tmp[-1,2]
            TDBmTTgeo <- tmp[,1]
            TDBmTTgeo0 <- TDBmTTgeo[1]
            TDBmTTgeo <- TDBmTTgeo[-1]
        }
        delT <- abs(TDBmTTgeo-TDBmTTgeo1)
        if(all(delT<tol)){
            break()
        }else{
            TDBmTTgeo1 <- TDBmTTgeo
        }
        jd <- rbind(tt0,cbind(tt[,1],tt[,2]+TDBmTTgeo/DAYSEC))
    }
    jd <- jd[-1,]
    if(j==Ntry){
        cat('Iteration exceed maximum!\n')
    }

####from terrestial to ephemeris coordinates
    dTT.dTDB.geo <- 1+zTDBmTTgeo
    dTDB.dTT.geo <- 1/dTT.dTDB.geo
    if(Par$ObsType=='ground'){
        OutGeo2obs <- time_Geo2obs(utc,tai,tt,Par)
    }else{
        OutGeo2obs <- list(rGO=SpaceObs[,1:3],vGO=SpaceObs[,4:6],ZenithICRS=NA,PoleICRS=NA,ut1=utc)
    }
    tdb0 <- tdb <- jd

###obs term for tt to tdb
    TDBmTTobs <- 0
    for(k in 1:Ntry){
        geo <- gen_GeoEph(tdb,DE=Par$DE)
        if(TtTdbMethod!='FBsofa'){
            TDBmTTobs <- rowSums(OutGeo2obs$rGO*dTDB.dTT.geo*geo$out[,6:8])/CKMPS^2#s
        }
        if((TtTdbMethod=='FBgeo' & Par$DE==405) | TtTdbMethod=='FB90'){
            tdb <- cbind(tt[,1],tt[,2]+(TDBmTTgeo+TDBmTTobs)/DAYSEC)#TDBmTTobs~2 mus
            tcb <-  sofa_Tdbtcb(tdb)
        }else{
            tmp <- time_Tt2tbObs(tcg,tt,TDBmTTgeo=TDBmTTgeo/DAYSEC,TDBmTTgeo0=TDBmTTgeo0/DAYSEC,TDBmTTobs=TDBmTTobs/DAYSEC,zTDBmTTgeo)
            tdb <-tmp$tdb
            tcb <- tmp$tcb
        }
        ddt <- abs(rowSums(tdb-tdb0))
        if(all(ddt<tol)){
            break()
        }else{
            tdb0 <- tdb
        }
    }

###ref: DE430 document equation 5: https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430_and_de431.pdf
###calculate zTDBmTTobs
    rGO <- OutGeo2obs$rGO*replicate(3,dTDB.dTT.geo)#km
    vGO <- OutGeo2obs$vGO#km/s
    rSG <- geo$out[,3:5]/au2km#au
    RSG <- sqrt(rowSums(rSG^2))
    uSG <- rSG/RSG
    ve <- geo$out[,6:8]#km/s
    ageo <- -uSG*4*pi^2/RSG^2
    TtTdbMethod <- Par$TtTdbMethod
    if(TtTdbMethod!='FBsofa'){
#    if(TRUE){
        zTDBmTTobsV <- rowSums((vGO/CKMPS)*(ve/CKMPS))
        zTDBmTTobsR <- rowSums(rGO/au2km*ageo)/Cauyr^2
        zTDBmTTobs <- zTDBmTTobsR+zTDBmTTobsV
    }else{
        zTDBmTTobs <- zTDBmTTobsR <- zTDBmTTobsV <- 0
    }
    Robs <- sqrt(Par$xtel^2+Par$ytel^2+Par$ztel^2)*IFTE.K#km; in TCB-compatible unit
    dzenith <- OutGeo2obs$vGO/Robs*DAYSEC#rad/day
    ##TDB at the observatory
    tcb <-  sofa_Tdbtcb(tdb)
    robust <- FALSE
    if(robust){
        ve.c <- ve/CKMPS
        rSO <- OutGeo2obs$rGO/au2km+rSG
        RSO <- sqrt(rowSums(rSO^2))
        RGO <- rGO/au2km
        w0E.c2 <- calZG(1,RSO)#only account for potential due to the Sun
        wLE.c2 <- 0#potential due to oblateness of external bodies; third order effect lead to fourth order dTT.dTDB effect
        var <- 1-(ELG-ELB)/(1-ELB)-(1-ELG)/(1-ELB)*(ve.c^2/2+w0E.c2+wLE.c2)
        dTT.dTDB <- 1-zTDBmTTgeo-zTDBmTTobs*((var+zTDBmTTobs)/var)
    }else{
        dTT.dTDB <- 1-zTDBmTTgeo-zTDBmTTobs
    }
    dTDB.dTT <- 1/dTT.dTDB
    dTCB.dTT <- IFTE.K*dTDB.dTT

###barrycorr approach
    return(c(list(JDtdb=tdb,JDtcb=tcb,geo=geo$out,dTT.dTDB=dTT.dTDB,dTDB.dTT=dTDB.dTT,dTCB.dTT=dTCB.dTT,dzenith=dzenith,emrat=geo$emrat,GM=geo$GM,TDBmTTgeo=TDBmTTgeo,TDBmTTobs=TDBmTTobs,zTDBmTTgeo=zTDBmTTgeo,zTDBmTTobs=zTDBmTTobs,zTDBmTTobsV=zTDBmTTobsV,zTDBmTTobsR=zTDBmTTobsR),OutGeo2obs))
}

time_Jd2cal <- function(jd){
####################################
## Julian Date to Gregorian year, month, day, fraction
## ref: https://uk.mathworks.com/matlabcentral/fileexchange/15285-geodetic-toolbox?focused=6513469&tab=function
## ref: Algorithm is from Fliegel and van Flandern (1968); also see julian.py
##
##  Input:
##     jd - Julian Date with day fraction
##
##  Output:
##     cal - Calendar date with day fraction
####################################

## Minimum and maximum allowed JD
   DJMIN = -68569.5
   DJMAX = 1e9

## Verify date is acceptable.
   if(any(rowSums(jd) < DJMIN) | any(rowSums(jd) > DJMAX)) return(NULL)

## Copy the date, big then small, and re-align to midnight. */
   jd1 <- jd[,1]+jd[,2]-0.5+1
#   f <- jd1-jd
   f <- (jd[,1]-as.integer(jd1))+jd[,2]-0.5+1
   jd <- as.integer(jd1)

## Express day in Gregorian calendar.
   l <- jd + 68569
   n <- (4 * l) %/% 146097
   l <- l-((146097 * n + 3) %/% 4)
   i <- (4000 * (l + 1)) %/% 1461001
   l <- l-(1461 * i) %/% 4 + 31
   k <- (80 * l) %/% 2447
   id <- l - (2447 * k) %/% 80
   l <- k %/% 11
   im <- k + 2 - 12* l
   iy <- 100 * (n - 49) + i + l
   fd <- f
   return(cbind(iy,im,id,fd))
}

time_Jd2calHms <- function(jd){
####################################
## Julian Date to Gregorian year, month, day, fraction
## ref: https://uk.mathworks.com/matlabcentral/fileexchange/15285-geodetic-toolbox?focused=6513469&tab=function
## ref: Algorithm is from Fliegel and van Flandern (1968); also see julian.py
##
##  Input:
##     jd - Julian Date with day fraction
##
##  Output:
##     cal - Calendar date with day fraction
####################################
    tmp <- time_Jd2cal(jd)
    fd <- tmp[,4]*DAYSEC
    h <- fd%/%3600
    m <- (fd-h*3600)%/%60
    s <- fd-h*3600-m*60
    return(cbind(tmp[,1:3,drop=FALSE],cbind(h,m,s)))
}

time_Utc2tb <- function(utc,Par){
####################################
## Convert UTC to TCB and TDB
##
## Input:
##   utc - 2-part JD[UTC] time
##
## Output:
##   JDtcb - 2-part TCB
##   JDtdb - 2-part TDB
##   JDtt - 2-part TT
##   JDtai - 2-part TAI
##   JDut1 - 2-part UT1
##   JDtcg - 2-part TCG
##   dTCB.dTT - dTCB/dTT
##   dTDB.dTT - dTDB/dTT
##   SO - Position and velocity vectors from SSB to observer (km; km/s)
##   SG - Position and velocity vectors from SSB to geocenter (km; km/s)
##   GO - Position and velocity vectors from geocenter to observer (km; km/s)
##   TDBmTTgeo - TDB-TT at geocenter
##   zTDBmTTgeo - d(TDB-TT)/dTDB at the geocenter
##   zTDBmTTobs - d(TDB-TT)/dTDB due to relativistic effect of the observer in the geocentric frame
##   zTDBmTTobsR - zTDBmTTobs due to geocenter-observer offset
##   zTDBmTTobsV - zTDBmTTobs due to Earth's rotation
##   leap - Leap seconds
##   xp - xp parameter in units of rad in Earth rotation model
##   yp - yp parameter in units of rad in Earth rotation model
##   zenith - Zenith direction in the ICRS frame
##   dzenith - dzenith/dt; not convert from TDB to TCB unit because it is always used to derive other observer-related quantities
##   emrat - Earth-Moon mass ratio
##   MO - Position and velocity vector from the Moon to observer
##   GM - Position and velocity vector from the geocenter to the Moon
####################################
###convert UTC to TT
    val <- time_Utc2tt(utc,Par)
    tt <- val$tt
    tai <- val$tai
    leap <- val$leap
    tcg <- sofa_Tttcg(tt)
    h <- 1e-12#day
    tmp <- time_Tt2tb(utc,tai,tcg,tt,Par)
    tdb <- tmp$JDtdb
    emrat <- tmp$emrat

###convert TDB quantities to TCB quantities; ref: https://syrte.obspm.fr/iauJD16/klioner.pdf
###TCB is recommended by TEMPO2, Edwards et al. 2006; however the 2006 IAU definition of TDB is probably as suitable as TCB, e.g. Eastman et al. 2010
    SG <- tmp$geo[,3:8,drop=FALSE]
    GM <- tmp$GM
    rGO <- tmp$rGO
###First convert all quantities from the ephemeris (TDB) unit to coordinate (TDB) unit; then provide output quantities according to the unit chosen by the user.
###    if(Par$Unit=='TCB'){
    SG[,1:3] <- SG[,1:3]*IFTE.K
#    rGO <- rGO*IFTE.K
    rGO <- rGO*tmp$dTCB.dTT
    GM[,1:3] <- GM[,1:3]*IFTE.K
###    }
    GO <- cbind(rGO,tmp$vGO)
    MO <- GO-GM#in coordinate unit
    SO <- SG+GO#in coordinate unit

###transfer into coordinate unit
    return(list(JDtcb=tmp$JDtcb,JDtdb=tdb,JDtt=tt,JDtai=tai,JDut1=tmp$ut1,JDtcg=tcg,dTCB.dTT=tmp$dTCB.dTT,dTDB.dTT=tmp$dTDB.dTT,SO=SO,SG=SG,GO=GO,vGO=tmp$vGO,zTDBmTTgeo=tmp$zTDBmTTgeo,zTDBmTTobs=tmp$zTDBmTTobs,TDBmTTgeo=tmp$TDBmTTgeo,zTDBmTTobsR=tmp$zTDBmTTobsR,zTDBmTTobsV=tmp$zTDBmTTobsV,leap=leap,xp=tmp$xp,yp=tmp$yp,zenith=tmp$ZenithIcrs,dzenith=tmp$dzenith,emrat=emrat,MO=MO,GM=GM))
}

time_ShapiroSolar <- function(OutBary,OutSB,Par){
####################################
## Shapiro delay in the solar system
## Although the time when the Sun deflect the light ray is different from the observed TDB by about 8 min
## The Sun's velocity with respect to the barycenter is about 11m/s and thus the Sun only moves by about 5.5 km and thus negligible.
##
## Input:
##   OurBary - Output of time_Utc2tb
##   OutSB - Output of astro_CalSB
##   Par - Input parameters
##
## Output:
##   dt.all - Shapiro delay due to all bodies in units of second
##   dt.list - Shapiro delay due to individual bodies in units of second
##   OL - Position and velocity vectors from the observer to lenses
####################################
    tdb <- OutBary$JDtdb
    ROB <- sqrt(rowSums(OutSB$rOB^2))
    dt.list <- list()
    dt.all <- 0
    ms <- Mssp
    ns <- c('Sun','Mercury','Venus','Earth','Moon','Mars','Jupiter','Saturn','Uranus','Neptune')
    Nl <- length(ns)
    if(!Par$PlanetShapiro) Nl <- 1
    Eph <- OL <- list()
    for(j in 1:Nl){
        if(ns[j]!='Earth' & ns[j]!='Moon'){
            eph <- gen_CalEph(tdb,body=ns[j],DE=Par$DE)
        }else if(ns[j]=='Earth'){
            eph <- OutBary$SG
        }else if(ns[j]=='Moon'){
            eph <- OutBary$SO-OutBary$GO
        }
        eph[,1:3] <- eph[,1:3]*IFTE.K#from ephemeris to coordinate length
        Eph[[ns[j]]] <- eph
        rSL <- eph[,1:3]/au2km
        rOL <- -OutBary$SO[,1:3]/au2km+rSL#au
        vOL <- (-OutBary$SO[,4:6]+eph[,4:6])/auyr2kms#au/yr
        VOL <- gen_CalLen(vOL)
        ROL <- gen_CalLen(rOL)
        ##cosine of angle between sun and observatory
        cos.psi <- rowSums(rOL*OutSB$rOB)/(ROL*ROB)
        dt <- -SRSLT*ms[j]*log(ROL*(1-cos.psi))#in second; TEMPO2 style
        dt.list[[ns[j]]] <- dt
        dt.all <- dt.all+dt
        tmp <- cbind(rOL,vOL,ROL,VOL,cos.psi,eph)
        tmp.tdb <- cbind(rOL/IFTE.K,vOL,ROL/IFTE.K,VOL,cos.psi,eph[,1:3]/IFTE.K,eph[,4:6])
        colnames(tmp) <- c('x.au','y.au','z.au','vx.auyr','vy.auyr','vz.auyr','ROL','VOL','cpsi','SLx.km','SLy.km','SLz.km','SLvx.kms','SLvy.kms','SLvz.kms')
        OL[[ns[j]]] <- tmp
    }
    return(list(dt.all=dt.all,dt.list=dt.list,OL=OL,Eph=Eph))
}

time_RoemerSolar <- function(OutBary,OutSB,SB=FALSE){
####################################
## Roemer delay (first order) and curvature delay (second order) in the solar system
##
## Input:
##   OutBary - Output of time_Utc2tb
##   OutSB - Output of astro_CalSB
##
## Output:
##   Roemer - Roemer delay due to all effects (day)
##   Roemer1 - First order Roemer delay (s)
##   Roemer2 - Second order Roemer delay (s)
##   Roemer3 - Third order Roemer delay (s)
##   RoemerT2 - Roemer delay using TEMPO2 method (s)
##
    rSO <- OutBary$SO[,1:3]/au2km#au
    rSB <- OutSB$SB[,1:3,drop=FALSE]
    RSB <- sqrt(rowSums(rSB^2))
    RSO <- sqrt(rowSums(rSO^2))
    ROB <- sqrt(rowSums(OutSB$rOB^2))
    Dt3 <- tt <- tt0 <- tt1 <- tt2 <- tt3  <- tt4 <- 0
    RSOpara <- rowSums(rSO*OutSB$uST)
###use uST as the reference unit vector
    if(SB){
        uR <- OutSB$uSB
        RR <- RSB
    }else{
        uR <- OutSB$uST
        RR <- OutSB$RST
    }
    rSOpara <- replicate(3,RSOpara)*uR
    rSOperp <- rSO-rSOpara
    Roemer2 <- 0.5*rowSums(rSOperp^2)/(RR*pc2au)*AULT#second
    Roemer1 <- -rowSums(rSO*uR)*AULT
    Rperp2 <- RSO^2-RSOpara^2
#    Dt3 <- 0.5*AULT*RSOpara*Rperp2/ROB^2/pc2au^2#less than 10^-10 s for a given moment, but is cummulative and thus is important
    OSpara <- -rowSums(rSO*uR)#au
    epsilon <- (RSO/RR/pc2au)^2+2*OSpara/RR/pc2au
    Roemer3 <- -1/16*epsilon^3*RSB*pc2au*AULT#second; third order
    if(Par$near & !SB){
        Roemer <- (OutSB$ROT-OutSB$RST)*pc2au*AULT#second
    }else{
        Roemer <- Roemer1+Roemer2+Roemer3
    }
    if(!is.null(OutSB$uSB.T2)){
        ##Rperp2 <- (RSO^2-RSOpara^2)*AULT#s
        ##In tempo2, r.perp is perpendicular to u0 rather than uOT or uOB or uSB; so the following will recover tempo2 prediction of parallax delay
        Rperp2 <- RSO^2 - (rSO[,1]*Par$u[1]+rSO[,2]*Par$u[2]+rSO[,3]*Par$u[3])^2
        tt3 <- +0.5*Rperp2*AULT*Par$plx*DMAS2R
        rlt <- -rSO*AULT
        tt0 <- rowSums(rlt*OutSB$uSB.T2$u0)
        tt1 <- rowSums(rlt*OutSB$uSB.T2$u1)
        tt2 <- rowSums(rlt*OutSB$uSB.T2$u2)
        tt4 <- rowSums(rlt*OutSB$uSB.T2$u3)
        tt <- tt0+tt1+tt2+tt3+tt4#rowSums(rlt*tempo$uSBt)
    }
    list(Roemer=Roemer,Roemer1=Roemer1,Roemer2=Roemer2,Roemer3=Roemer3,RoemerT2=list(all=tt,dt0=tt0,dt1=tt1,dt2=tt2,dt3=tt3,dt4=tt4))
}

time_ShapiroTarget <- function(m2,e,U,sini,omega){
    r <- (4*pi^2)*m2/Cauyr^3*DJY*DAYSEC#second
    brace <- 1-e*cos(U)-sini*(sin(omega)*(cos(U)-e)+sqrt(1-e^2)*cos(omega)*sin(U))
    dt <- -2*r*log(brace)#day
    return(dt)
}

time_VacuumDelay <- function(dRSB){
####################################
## Vacuum delay in interstellar space
##
## Input:
##   dRSB - distance between the SSB and the TSB (pc)
##
## Output:
##   vacuum delay (day)
####################################
    dRSB*pc2au*AULT#second
}

time_EinsteinTarget <- function(gamma,U){
####################################
## Einstein delay in the Target system
## ref: Blandford et al. 1975, eqs. (2.1)-(2.13)
##
## Input:
##   gamma - Einstein delay parameter
##   U - Eccentricity anormaly
##
## Output:
##  Einstein delay in the target system
####################################
    gamma*sin(U)
}

time_Ta2te <- function(OutBary,Par){
####################################
## Convert from arrival coordinate time at the telescope to the coordinate time of light emission from the surface of the target
##
## Input:
##   OutBary - The output of time_Utc2tb
##   Par - Input parameters
##
## Output:
##   tauE - Proper time of light emission (2-part JD)
##   tB - Coordinate arrival time of light at the TSB (2-part JD)
##   tS - Coordinate arrival time of light at the SSB (2-part JD)
##   BJDtcb - BJD[TCB]
##   BJDtdb - BJD[TDB]
##   TargetDelay - Total delay in the target system
##   RoemerTarget - Roemer delay in the target system
##   EinsteinTarget - Einstein delay in the target system
##   ShapiroTarget - Shapiro delay in the target system
##   RoemerSolar - Roemer delay in the Solar System using uST as the reference direction
##   RoemerSB - Roemer delay in the Solar System using uSB as the reference direction
##   ShapiroSolar - Shapiro delay in the Solar System
##   ShapiroPlanet - Shapiro delay due to individual Solar System bodies
##   OL - Position and velocity vectors from the observer to the Solar System bodies or lenses
##   ShapiroTarget - Shapiro delay in the target system
##   AbeTarget - aberration delay in the target system
##   VacuumIS - Vacuum delay
##   EinsteinIS - Einstein delay due to the motion of the TSB relative to the SSB
##   rOB - Position vector from the observer to the TSB (pc)
##   rBT - Position vector from the TSB to the target (au)
##   uOB - Direction from the observer to the TSB
##   vOB - Velocity vector from the observer to the TSB (au/yr)
##   RBT - Disance from the TSB to the target (au)
##   uBT - Direction from the TSB to the target
##   rOT - Position vector from the observer to the target (pc)
##   rST - Position vector from the SSB to the target (pc)
##   SB - Position and velocity vectors from the SSB to the TSB (pc,au/yr)
##   uOT - Direction vector from the observer to the target
##   vOT - Velocity vector from the observer to the target (au/yr)
##   uST - Direction vector from the SSB to the target
##   vST - Velocity vector from the SSB to the target (au/yr)
##   uSB - Direction vector from the SSB to the TSB
##   rSB - Position vector from the SSB to the TSB (pc)
##   vSB - Velocity vector from the SSB to the TSB (au/yr)
##   vBT - Velocity vector from the SSB to the target (au/yr)
##   U - Eccentricity anomaly
##   OutBT - Output of time_CalBT()
##   RoemerOrder - Different orders of Roemer delays in the Solar System
##   Ztropo - Redshift due to Tropospheric delay
##   TropoDelay - Tropospheric delay
##   TropoDelayT2 - Tropospheric delay using uSB calculated by TEMPO2 method
##   elevation - Elevation angle (rad)
##   delevation - delevation/dt
##   delevation1 - d(elevation1)/dt
##   ZenIn - Real zenith for incident light ray (scaler; rad)
##   uSB.T2 - uSB calculated using the TEMPO2 method (precision up to 2 orders Roemer delay)
##   RoemerT2 - Roemer delay calculated by TEMPO2 method (s)
####################################
##derive quantities from the output of time_Utc2tb
    tpos <- Par$tpos

###iteration to determine tS and tB
    Ntry <- 10
#    if(!Par$binary)   Ntry <- 1

####initial guess of tS and tB
    OS <- -cbind(OutBary$SO[,1:3]/au2km,OutBary$SO[,4:6]/auyr2kms)
    ts <- time_TsTb1(tO=OutBary$JDtcb,OS,Par)
    ts <- cbind(OutBary$JDtcb[,1],ts[,1])
    tS <- tB <- time_ChangeBase(cbind(OutBary$JDtcb[,1],ts[,2]),1)
###If not scaling, the full transformation from tSSB to tTSB (eqs. 43-48 in E06) will be adopted; However, since these transformations only change the scaling of various parameters, it is not detectable in exoplanet data. If post-Newtonian theores are tested using PEXO, this scaling factor could influence the result and thus should be included.
    tol <- 1e-12
    Ndata <- nrow(OutBary$JDtcb)
    tauE0 <- tB
    for(j in 1:Ntry){
#        cat('j=',j,'\n')
        ##calculate SB given tB
        OutSB <- astro_CalSB(tB,Par)#pc, au/yr
        rOB <- OS[,1:3,drop=FALSE]/pc2au+OutSB$SB[,1:3,drop=FALSE]#pc
        uOB <- gen_CalUnit(rOB)
        ROB <- gen_CalLen(rOB)
        if(j==1){
            uST <- OutSB$uSB
            rST <- OutSB$rSB
            uOT <- uOB
            ROT <- ROB
        }
        RST <- gen_CalLen(rST)
        vOB <- OS[,4:6,drop=FALSE]+OutSB$SB[,4:6,drop=FALSE]#au/yr
        OutSB <- c(OutSB,list(rOB=rOB,uOB=uOB,vOB=vOB,ROB=ROB,uOT=uOT,uST=uST,RST=RST,rST=rST))
        if(Par$near){
#for solar system objects, use (ROT-RST)/c
            OutSB <- c(OutSB,list(ROT=ROT))
        }
        if(Par$ObsType=='ground'){
            OutEle <- astro_CalElevation(OutBary$zenith,OutBary$dzenith,uOT=uOT)#uOB is only used for iteration and more precise elevation angle is calculated using uOT later
        }else{
            OutEle <- list(elevation=NA,delevation=NA,ZenIn=NA)
        }
        if(Par$CompareT2){
            OutEleT2 <- astro_CalElevation(OutBary$zenith,OutBary$dzenith,uOT=t(replicate(Par$Nepoch,Par$u)))#this is T2 method which ignores proper motion
        }else{
            OutEleT2 <- list(elevation=NA,delevation=NA,ZenIn=NA)
        }
        ##
        ##calculate tS given tO and SB
        OutBjd <- time_Jd2bjd(OutBary,OutSB,OutEle,Par=Par,CompareT2=Par$CompareT2)
        if(Par$CompareT2 & Par$ObsType=='ground'){
            OutBjdT2 <- time_Jd2bjd(OutBary,OutSB,OutEleT2,Par=Par,CompareT2=FALSE)
        }else{
            OutBjdT2 <- NULL
        }
        tS <- BJDtcb <- OutBjd$BJDtcb
        ##calculate
        if(Par$SBscaling){
            dRSB <- OutSB$dRSB
            VacuumIS <- time_VacuumDelay(dRSB)#second; always set zero?
            v2 <- (Par$pmra^2+Par$pmdec^2)/Par$plx^2+(Par$rv/auyr2kms)^2
            v2.c2 <- v2/Cauyr^2
            EinsteinIS <- 0.5*v2.c2*((tS[,1]-Par$tpos)+tS[,2]-VacuumIS/DAYSEC)*DAYSEC#second
###ref. E06 Eqn. 44-45
            tB1 <-cbind(tS[,1],tS[,2]-(VacuumIS+EinsteinIS)/DAYSEC)
            tB1-tB
            if(all(abs(time_T2mT2(tB1,tB))<tol)) break()
            tB <- tB1
        }else{
            tB <- tS
            VacuumIS <- EinsteinIS <- 0
        }

        BJDtdb <-OutBjd$BJDtdb
        RoemerSolar <- OutBjd$Roemer
        RoemerSB <- OutBjd$RoemerSB$Roemer
#        cat('RoemerSolar-RoemerSB=',head(RoemerSolar-RoemerSB),'s\n')
        ShapiroSolar <- OutBjd$ShapiroSolar$dt.all
        ShapiroPlanet <- OutBjd$ShapiroSolar$dt.list

        tB <- time_ChangeBase(tB)
                                        #tB <- cbind(tS[,1],tB1)
        if(Par$binary){
            OutBT <- gen_CalBT(tB=tB,Par)
            BT <- OutBT$BT
            U <- OutBT$U
            vBT <- BT[,4:6,drop=FALSE]#auy/r
            rBT <- BT[,1:3,drop=FALSE]#au
            RBT <-sqrt(rowSums(rBT^2))#au
            uBT <- rBT/RBT
        }else{
            BT <- array(0,dim=c(Par$Nepoch,6))
            vBT <- rBT <- array(0,dim=c(Par$Nepoch,3))
            RBT <- rep(0,Par$Nepoch)
            uBT <- array(0,dim=c(Par$Nepoch,3))
            U <- NA
            OutBT <- list(BT=BT,vBT=vBT,RBT=RBT,uBT=uBT,U=U)
        }
        rOT <- rOB+rBT[,1:3,drop=FALSE]/pc2au#pc
        ROT <- sqrt(rowSums(rOT^2))
        vOT <- vOB+vBT[,1:3,drop=FALSE]#auy/r
        uOT <- rOT/ROT
        rSB <- OutSB$SB[,1:3,drop=FALSE]
        vSB <- OutSB$SB[,4:6,drop=FALSE]
        rST <- rSB+rBT[,1:3,drop=FALSE]/pc2au#pc
        vST <- vSB+vBT[,1:3,drop=FALSE]#au/yr
        RST <- sqrt(rowSums(rST^2))
        uST <- rST/RST
        ROB <-sqrt(rowSums(rOB^2))
        dD <- (ROT-ROB)*pc2au#au
        RoemerTarget <- dD*AULT#second
        if(Par$binary){
            robust <- TRUE
            if(Par$BinaryModel=='kepler'){
                tmp <- gen_mass2dd(m=Par$mTC,m2=Par$mC,x=Par$DDGR$x0,ecc=Par$DDGR$ecc,an=2*pi/Par$DDGR$pb,BinaryUnit=Par$BinaryUnit)
                EinsteinTarget <- time_EinsteinTarget(tmp$gamma*DJY*DAYSEC,OutBT$U)#second
                ShapiroTarget <- time_ShapiroTarget(Par$mC,e=Par$e,U=OutBT$U,sini=sin(Par$I),omega=Par$omegaT)#second
                AbeTarget <- rep(0,Par$Nepoch)
            }else{
                EinsteinTarget <- time_EinsteinTarget(OutBT$gamma*DJY*DAYSEC,OutBT$U)#second
                ShapiroTarget <- OutBT$Ds#second
                AbeTarget <-OutBT$Da#second
#                RoemerEinsteinTarget <- OutBT$Dre
            }
        }else{
            AbeTarget <- ShapiroTarget <- RoemerTarget <- EinsteinTarget <- 0
        }
        TargetDelay <- RoemerTarget+EinsteinTarget+ShapiroTarget#total delay in the target system
        tauE <- cbind(tB[,1],tB[,2]-TargetDelay/DAYSEC)#light emission proper time (do not consider beaming); julian day
        dtauE <- abs(time_T2mT2(tauE,tauE0))
        if(all(dtauE<1e-9)) break()
        tauE0 <- tauE
    }

###derive more quantities for astrometry and RV modeling
    if(Par$binary){
        rBT <- BT[,1:3]#au/yr
        rTC <- -rBT*Par$mTC/Par$mC#au
        rOC <- rOT+rTC/pc2au#pc
        uOC <- gen_CalUnit(rOC)
        rTC <- -rBT*Par$mTC/Par$mC
        vTC <- -vBT*Par$mTC/Par$mC
        vSC <- vST+vTC#au/yr
        VTC <- sqrt(rowSums(vTC^2))#au/yr
        RTC <- sqrt(rowSums(rTC^2))#au
        rSC <- rST+rTC/pc2au#pc
    }else{
        rTC <- 0
        rOC <- rOT
        uOC <- uOT
        rSC <- rST
        vSC <- vST
    }

###return
    return(list(tauE=tauE,tB=tB,tS=tS,BJDtcb=BJDtcb,BJDtdb=BJDtdb,TargetDelay=TargetDelay,RoemerTarget=RoemerTarget,EinsteinTarget=EinsteinTarget,ShapiroTarget=ShapiroTarget,RoemerSolar=RoemerSolar,RoemerSB=RoemerSB,ShapiroSolar=ShapiroSolar,ShapiroPlanet=ShapiroPlanet,OL=OutBjd$ShapiroSolar$OL,Eph=OutBjd$ShapiroSolar$Eph,ShapiroTarget=ShapiroTarget,AbeTarget=AbeTarget,VacuumIS=VacuumIS,EinsteinIS=EinsteinIS,rOB=rOB,rBT=rBT,uOB=uOB,vOB=vOB,RBT=RBT,uBT=uBT,BT=BT,rOT=rOT,rST=rST,SB=OutSB$SB,uOT=uOT,vOT=vOT,uST=uST,vST=vST,uSB=OutSB$uSB,rSB=rSB,vSB=vSB,vBT=vBT,U=U,OutBT=OutBT,RoemerOrder=OutBjd$RoemerOrder,Ztropo=OutBjd$Ztropo,TropoDelay=OutBjd$TropoDelay,elevation=OutEle$elevation,delevation=OutEle$delevation,ZenIn=OutEle$ZenIn,uSB.T2=OutSB$uSB.T2,RoemerT2=OutBjd$RoemerT2,TropoDelayT2=OutBjdT2$TropoDelay,elevationT2=OutEleT2$elevation,delevationT2=OutEleT2$delevation,ZenInT2=OutEleT2$ZenIn,rTC=rTC,rOC=rOC,uOC=uOC,rSC=rSC,vSC=vSC))
}

time_ToJD <- function(epoch){
####################################
## Transform epochs into JD format
##
## Input:
##   epoch - 1-part or 2-part epochs
##
## Output:
##   epoch - 1-part or 2-part JD epochs
####################################
    if(is.null(dim(epoch))){
        ##1-part epochs
        if(all(epoch<2000000)) epoch <- epoch+DJM0
    }else{
        ##2-part epochs
        if(all(rowSums(epoch)<2000000)) epoch[,1] <- epoch[,1]+DJM0
    }
    epoch
}
