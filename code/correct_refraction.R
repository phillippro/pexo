###This routine is for the correction of atmospheric refraction effects in RV modeling
library(orthopolynom)
library(pracma)
source('constants.R')
####SOFA functions
source('nut00.R')
source('sofa.R')
####Other functions
source('atmosphere.R')
source('orbit.R')
source('deflect_redshift.R')
source('delay.R')
source('OrbitFunction.R')
source('tropo.R')
source('binary_models.R')

args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    target <- args[1]
    set <- args[2]
}else{
####settings
    ##target <- '39091'
#    target <- 'HD10700'
    target <- '10700'
    ##target <- 'WASP29'
    ##target <- '6793'
    ##target <- 'HD10700'
    ##target <- '39091'
#    set <- 'HARPS'
    set <- 'PFS'
}

###global parameters
TT.standard <- 'BIPM16'
UNITS <- 'TCB'
DE <- 430
eopType <- '06'

###load observation data
###load PFS fits header data for correction
#tab  <-  read.table('../../data/PFS/PFS_info.dat',header=TRUE,check.names=FALSE)
if(set=='PFS'){
    fin <- paste0('../../data/PFS/PFSfits_',target,'.dat')
}else{
    fin <- paste0('../../data/combined/',target,'/',target,'_',set,'.dat')
}
cat(fin,'\n')
tab  <-  read.table(fin,header=TRUE,check.names=FALSE)
if(nrow(tab)>1e3) tab <- tab[tab[,1]>(55525.63217-1) & tab[,1]<(57355.68011+1),]
#pfs <- pfs[1,,drop=FALSE]
if(set=='PFS'){
    N <- nrow(tab)
    cal <- array(NA,dim=c(N,6))
    ind.rm <- c()
    for(k in 1:N){
        tmp <- as.numeric(unlist(strsplit(as.character(tab[k,'photon-weighted-mid-time']),split=':')))
        if(any(is.na(tmp)) | length(tmp)<3){
            tmp <- as.numeric(unlist(strsplit(as.character(tab[k,'geometric-mid-time']),split=':')))
        }
        if(any(is.na(tmp)) | length(tmp)<3){
            ind.rm <- c(ind.rm,k)
        }else{
            jd <- as.integer(unlist(strsplit(as.character(tab[k,'UT-DATE']),split='-')))
            cal[k,] <- c(jd,tmp)
        }
    }
    if(length(ind.rm)>0){
        tab <- tab[-ind.rm,]
        cal <- cal[-ind.rm,]
    }

###convert to JD
    utc <- time_CalHms2JD(cal)
    utc <- change.base(utc,denom=1)

    ##position
    ra <- tab[,'RA-D']*pi/180
    dec <- tab[,'DEC-D']*pi/180
    ra[2:length(ra)] <- ra[1]
    dec[2:length(dec)] <- dec[1]
    tdk <- tab[,'AIRTOUT']+273.15#air temperature
    pmb <- tab[,'AIRPRESS']#air pressure
    rh <- 0.6#relative humidity; averaged value from the weather data: https://weather.lco.global/#/bpl
                                        #wl <- 0.8#micrometer
    wl <- 0.5#micrometer
    phi <- tab[,'SITELAT']*pi/180#latitute of observer in radian
    hm <- tab[,'SITEALT']#altitude of observer in meter
    airmass <- tab[,'AIRMASS']
    atm.par <- list(hm=hm,tdk=tdk,pmb=pmb,rh=rh,wl=wl,phi=phi,airmass=airmass)
###observatory site in PEXO formate
#    llh <- tab[1,c('SITELONG','SITELAT','SITEALT')]
###From Paul's obssite.pro code
    lat = -29.013983
    lon = -70.692633
    ht =   2407.92
    llh <- c(lon,lat,ht)
    llh[1:2] <- llh[1:2]*pi/180
}else if(set=='HARPS'){
    if(max(tab[,1])<2400000) tab[,1] <- tab[,1]+2400000#HARPS data uses this definition of MJD
    fithead <- read.table(paste0('../../data/',set,'/',target,'_HARPS_info.dat'),header=TRUE)
    inds <- match(round(tab[,1],5),round(bjd,5))
    tab <- cbind(tab[which(!is.na(inds)),1:3],fithead[inds[!is.na(inds)]])
    llh <- c(-70.7345*pi/180,-29.2584*pi/180,2400)
    hm <- llh[3]
    phi <- llh[2]
    wl <- 0.5
    rh <- tab[,'rhum']
    pmb <- tab[,'']
    tdk <- 273+tab[,'T']
}else{
#    if(max(tab[,1])<2400000) tab[,1] <- tab[,1]+DJM0
    if(max(tab[,1])<2400000) tab[,1] <- tab[,1]+2400000#HARPS data uses this definition of MJD
    utc <- cbind(tab[,1],0)
    utc <- time_ChangeBase(utc,denom=1)
    ra <- rep(026.00930287667,nrow(tab))*pi/180
    dec <- rep(-15.93379865094,nrow(tab))*pi/180
    llh <- c(-70.7345*pi/180,-29.2584*pi/180,2400)
    hm <- llh[3]
    phi <- llh[2]
    wl <- 0.5
    rh <- 0.6
    pmb <- 771#milibar
    tdk <- 273+14#K
}
Nobs <- nrow(utc)

ObsPar <- list(data=llh,type='llh')
if(ObsPar$type=='xyz'){
    ObsPar$robs <-ObsPar$data
    ObsPar$llh <- as.numeric(Gc2gd(1,ObsPar$robs))
}
if(ObsPar$type=='llh'){
    ObsPar$llh <- ObsPar$data
    ObsPar$robs <- as.numeric(Gd2gc(1,ObsPar$llh))
    ObsPar$type <- 'xyz'
}
ObsPar$R <- sqrt(sum(ObsPar$robs^2))#meter

####load data about the Earth roation model
source('load_data.R')

##calculate un-refrected zenith
obs <- time_Utc2tb(utc,ObsPar,DE=DE,type=eopType,TaiType='scale')
zenith <- obs$robs/(ObsPar$R*1e-3)
dzenith <- obs$dzenith#rad/day; omega
SO <- obs$state
SO[,1:3] <- SO[,1:3]/au2km
SO[,4:6] <- SO[,4:6]/auyr2kms
rSO <- SO[,1:3,drop=FALSE]
OS <- -SO
SO[,1:3] <- SO[,1:3,drop=FALSE]
ROS <- sqrt(rowSums(SO[,1:3,drop=FALSE]^2))
###
uOB <- cbind(cos(dec)*cos(ra),cos(dec)*sin(ra),sin(dec))
elevation <- asin(rowSums(zenith*uOB))
delevation <- 1/sqrt(1-(rowSums(zenith*uOB))^2)*(rowSums(dzenith*uOB))
ind.rm <- which(elevation<0)
Nrm <- length(ind.rm)
if(Nrm>0){
    cat('elevation angle < 0 for ',length(ind.rm),'UTC times and are removed!!!\n')
    elevation <- elevation[-ind.rm]
    delevation <- delevation[-ind.rm]
    if(exists('airmass')){
        if(length(airmass)>1){
            airmass <- airmass[-ind.rm]
        }
    }
    if(length(phi)>1){
        phi <- phi[-ind.rm]
    }
    if(length(hm)>1){
        hm <- hm[-ind.rm]
    }
    if(length(pmb)>1){
        pmb <- pmb[-ind.rm]
    }
    if(length(tdk)>1){
        tdk <- tdk[-ind.rm]
    }
    if(nrow(zenith)>1){
        zenith <- zenith[-ind.rm,]
    }
    uOB <- uOB[-ind.rm,]
    OS <- OS[-ind.rm,]
    utc <- utc[-ind.rm,]
    tab <- tab[-ind.rm,]
}
Nin <- Nobs-Nrm
if(length(hm)==1) hm <- rep(hm,Nin)
if(length(tdk)==1) tdk <- rep(tdk,Nin)
if(length(pmb)==1) pmb <- rep(pmb,Nin)
if(length(rh)==1) rh <- rep(rh,Nin)
if(length(wl)==1) wl <- rep(wl,Nin)
if(length(phi)==1) phi <- rep(phi,Nin)
atm.par <- list(hm=hm,tdk=tdk,pmb=pmb,rh=rh,wl=wl,phi=phi)

zen <- pi/2-elevation
dzen <- -delevation
#zen[zen==0] <- 1e-10

#tmp <- cal.refraction(zen,atm.par,type='refcoq')
tmp <- cal.refraction(zen,atm.par,type='refro')
Ratm <- tmp$R
tmp1 <- cal.refraction(zen+dzen*0.001,atm.par,type='refro')
dRatm <- (Ratm-tmp1$R)/0.001#rad/day
Ratm.vec <- (zenith-rowSums(zenith*uOB)*uOB)/sin(zen)*Ratm
RV.ref <- rowSums(Ratm.vec*OS[,4:6,drop=FALSE])*auyr2kms*1e3#m/s

if(FALSE){
    tmp0 <- cal.refraction(zen,atm.par,type='refcoq')
    Ratm0 <- tmp$R
    Ratm0.vec <- (zenith-rowSums(zenith*uOB)*uOB)/sin(zen)*Ratm
    RV0.ref <- rowSums(Ratm.vec*OS[,4:6,drop=FALSE])*auyr2kms*1e3#m/s

    tmp1 <- cal.refraction(zen,atm.par,type='B82+')
    Ratm1 <- tmp1$R
    Ratm1.vec <- (zenith-rowSums(zenith*uOB)*uOB)/sin(zen)*Ratm1
    RV1.ref <- rowSums(Ratm1.vec*OS[,4:6,drop=FALSE])*auyr2kms*1e3#m/s

    tmp2 <- cal.refraction(zen,atm.par,type='S86+')
    Ratm2 <- tmp2$R
    Ratm2.vec <- (zenith-rowSums(zenith*uOB)*uOB)/sin(zen)*Ratm2
    RV2.ref <- rowSums(Ratm2.vec*OS[,4:6,drop=FALSE])*auyr2kms*1e3#m/s

    tmp3 <- cal.refraction(zreal=zen,atm.par=atm.par,type='airmass')
    Ratm3 <- tmp1$R
    Ratm3.vec <- (zenith-rowSums(zenith*uOB)*uOB)/sin(zen)*Ratm1
    RV3.ref <- rowSums(Ratm.vec1*OS[,4:6,drop=FALSE])*auyr2kms*1e3#m/s
}

mjd.utc <- time_Jd2mjd(utc)
tropo <- time_TropoDelay(mjd.utc,elevation,delevation,atm.par)
Dt.tropo <- tropo$dt/DAYSEC
z.tropo <- tropo$z
RV.tropo <- z.tropo*CMPS

source('plot_refraction.R')
