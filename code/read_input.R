library("optparse")
option_list  <-  list(
  make_option(c("-m", "--mode"), type="character", default='emulate',
              help="PEXO mode: emulate or fit [optional; default=%default]", metavar="character"),
  make_option(c("-c", "--component"), type="character", default='TAR',
              help="PEXO model component: timing (T), astrometry (A), radial velocity (R) and their combinations [optional; default=%default]", metavar="character"),
	make_option(c("-t", "--time"), type="character", default=NULL,help='Timing file: epochs or times could be in 1-part or 2-part JD[UTC] format [mandatory if mode=emulate]'),
	make_option(c("-p", "--par"), type="character", default=NULL,help='Parameter file: parameters for models, observatory, for Keplerian/binary motion [mandatory]'),
	make_option(c("-v", "--var"), type="character", default=NULL,help='Output variables [optional; default=%default]'),
	make_option(c("-o", "--out"), type="character", default="out.txt",
              help="Output file name: relative or absolute path [optional; default= %default]", metavar="character"),
	make_option(c("-f", "--figure"), type="character", default=TRUE,
              help="Output figure and verbose: FALSE or TRUE [optional; default= %default]", metavar="character")
)

opt_parser  <-  OptionParser(option_list=option_list)
opt  <- parse_args(opt_parser)

if(!is.null(opt$time)){
    if(is.null(opt$time) & opt$mode=='emulate'){
        print_help(opt_parser)
        stop("Error: at least one argument must be supplied (input timing file).n", call.=FALSE)
    }

    if(is.null(opt$par) & opt$mode=='emulate'){
        print_help(opt_parser)
        stop("Error: at least one argument must be supplied (input parameter file).n", call.=FALSE)
    }

}else{
###gaia80yrby10day corresponds to Fig. 14 and 15 in the paper
#    opt$time <- '../input/time_gaia80yrby10day.txt'
#    opt$time <- '../input/time1.txt'
#    opt$time <- '../input/hubble.tim'
    opt$time <- '../input/TCpfs.tim'
###hip80yrby10day corresponds to Fig. 17 in the paper
#    opt$time <- '../input/hip80yrby10day.tim'
#    opt$time <- '../input/gaia80yrby10day.tim'
#    opt$time <- '../input/mjd42000to52000by10day.tim'
#    opt$time <- '2456640.5 2458462.5 0.5'
#    opt$time <- '48620 58850 100'
#    opt$time <- '38300 57600 100'
##low zenith (or high elevation) version of ../input/time_hip80yrby10day.txt
#    opt$time <- '../input/time_hip80yrby10day_LowZ.txt'
#    opt$time <- '../input/time_hip80yrby10day.txt'
#    opt$time <- '../input/time_PythonTest.txt'
#    opt$time <- '~/Documents/projects/ISO/borisov/time_horizon.txt'
#    opt$par <- '../input/ACAgaia.par'
#    opt$par <- '../input/TC_Fig11b.par'
#    opt$par <- '../input/TC_Fig11a.par'
#    opt$par <- '../input/PSR_J0740+6620.par'
#    opt$par <- '../input/S2_gravity.par'
#    opt$par <- '../input/S2_E03.par'
#    opt$par <- '../input/HD197461.par'
    opt$par <- '../input/TCpfs.par'
#    opt$var <- c('JDtai','JDut1','JDtt','JDtdb','JDtcb','BJDtdb','BJDtcb','tE','tS','RoemerSolar','RoemerTarget','EinsteinTarget','EinsteinIS','EinsteinTarget','elevation','TropoDelay','TargetDelay','VacuumIS','ShapiroSolar','ShapiroTarget')
    opt$var <- c('BJDtcb','BJDtdb','RvTot')
}
###Example:
##Rscript pexo.R -m emulate -t ../input/time1.txt -p ../input/model_ACAhip.par -v 'JDtai JDut1 JDtt JDtdb' -o out1.txt
##

###Read timing file
if(!file.exists(opt$time)){
    s <- unlist(strsplit(opt$time,split=' '))
    s <- s[s!='']
    StrError <- 'Error: the UTC timing data cannot be found or the UTC sequence cannot be generated from the argument after -t !\n'
    if(length(s)==3){
        ts <- try(seq(as.numeric(s[1]),as.numeric(s[2]),by=as.numeric(s[3])),TRUE)
        if(class(ts)=='try-error'){
            stop(StrError,call. = FALSE)
        }else{
            utc <- time_ChangeBase(cbind(time_ToJD(ts),0))
        }
    }else{
        stop(StrError,call. = FALSE)
    }
}else{
    utc <- read.table(opt$time)
    if(ncol(utc)==1){
        utc <- cbind(utc,0)
    }
    utc <- time_ToJD(utc)
    utc <- time_ChangeBase(utc,1)
}

###Read parameter file
try(if(!file.exists(opt$par)) stop('Error: file does not exist!'))
par0  <- read.table(opt$par)
pars <- as.character(par0[,-1])
cn <- names(pars) <- as.character(par0[,1])

######Get parameter values from the parameter file
Par <- list()

Par$component <- opt$component
Par$figure <- as.logical(opt$figure)
###Read global parameters
if(any(cn=='name')){
    Par$star <- gsub(' ','',pars['name'])
}else{
    Par$star <- substr(opt$par,1,5)
}

if(any(cn=='RefType')){
    if(pars['RefType']!='refro' & pars['RefType']!='refco' & pars['RefType']!='refcoq' & pars['RefType']!='none') stop(paste0('In',opt$par,', ref value is not valid, it should be refro, refco, or refcoq!'),call. = FALSE)
    Par$RefType <- pars['RefType']
}else{
    Par$RefType <- 'none'
}

##   EopType - Version of EOP parameter file
if(any(cn=='EopType')){
    if(pars['EopType']!='2000B' & pars['EopType']!='2006' & pars['EopType']!='approx') stop(paste0('In',opt$par,', EopType value is not valid, it should be 2000B, 2006, or approx!'),call. = FALSE)
    Par$EopType <- pars['EopType']
}else{
    Par$EopType <- '2006'
}

if(any(cn=='TaiType')){
    if(pars['TaiType']!='scale' & pars['TaiType']!='instant') stop(paste0('In',opt$par,', TaiType value is not valid, it should be instant or scale!'),call. = FALSE)
    Par$TaiType <- pars['TaiType']
}else{
    Par$TaiType <- 'instant'
}

if(any(cn=='TtType')){
    if(pars['TtType']!='BIPM' & pars['TtType']!='TAI') stop(paste0('In',opt$par,', TtType value is not valid, it should be BIPM or TAI!'),call. = FALSE)
    Par$TtType <- pars['TtType']
}else{
    Par$TtType <- 'BIPM'
}

if(any(cn=='unit')){
    if(pars['unit']!='TCB' & pars['unit']!='TDB') stop(paste0('In',opt$par,', unit value is not valid, it should be TCB or TDB!'),call. = FALSE)
    Par$Unit <- pars['unit']
}else{
    Par$Unit <- 'TCB'
}

##   DE - Version of JPL ephemeris
if(any(cn=='DE')){
##Remove t from the DE value and thus DE430t is the same as DE430. For both, the DE430t ephemerides would be used.
    Par$DE <- as.integer(gsub('t','',pars['DE']))
}else{
    Par$DE <- 430
}

if(any(cn=='TtTdbMethod')){
    if(pars['TtTdbMethod']!='eph' & pars['TtTdbMethod']!='FB01' & pars['TtTdbMethod']!='FBgeo' & pars['TtTdbMethod']!='FB90') stop(paste0('In',opt$par,', TtTdbMethod value is not valid, it should be eph or FB01 or FBsofa or FBgeo or FB90 or analytic!'),call. = FALSE)
    Par$TtTdbMethod <- pars['TtTdbMethod']
}else{
    Par$TtTdbMethod <- 'eph'
}

## SBscaling - Whether or not transform the coordinate time at the SSB (tS) to the coordinate time at the TSB (tB). Since it is a linear transformation, it could be absorbed into the fit of orbital period or other time-scale parameters. So we follow TEMPO2 not to do this transformation/scaling by using scaling=FALSE as default.
if(any(cn=='SBscaling')){
    if(pars['SBscaling']!='TRUE' & pars['SBscaling']!='FALSE') stop(paste0('In',opt$par,', SBscaling value is not valid, it should be TRUE or FALSE!'),call. = FALSE)
    Par$SBscaling <- as.logical(pars['SBscaling'])
}else{
    Par$SBscaling <- NULL
}

##   planet - Whether or not the planet Shapiro delay is calculated.
if(any(cn=='PlanetShapiro')){
    if(pars['PlanetShapiro']!='TRUE' & pars['PlanetShapiro']!='FALSE') stop(paste0('In',opt$par,', PlanetShapiro value is not valid, it should be TRUE or FALSE!'),call. = FALSE)
    Par$PlanetShapiro <- as.logical(pars['PlanetShapiro'])
}else{
    Par$PlanetShapiro <- NULL
}

###The method used for RV modeling, analytical is more efficient and thus is default.
if(any(cn=='RVmethod')){
    if(pars['RVmethod']!='numerical' & pars['RVmethod']!='analytical') stop(paste0('In',opt$par,', RVmethod value is not valid, it should be analytical or numerical!'),call. = FALSE)
    Par$RVmethod <- pars['RVmethod']
}else{
    Par$RVmethod <- 'analytical'
}

if(any(cn=='Tstep')){
    Par$Tstep <- as.numeric(pars['Tstep'])
}else{
    if(Par$RVmethod=='numerical') cat('warning: No time step is given for numerical modeling of RV and a step of 0.01 days would be used!\n')
    Par$Tstep <- 0.01
}

if(Par$RVmethod=='numerical'){
    ts <- c()
    for(j in 1:nrow(utc)){
        ts <- rbind(ts,rbind(c(utc[j,1],utc[j,2]-Par$Tstep),utc[j,],c(utc[j,1],utc[j,2]+Par$Tstep)))
    }
    utc <- ts
}

Par$Nepoch <- nrow(utc)
##   CompareT2 - whether or not to compare TEMPO2 Roemer delay (second order) and PEXO Roemer delay (>=3nd order)
if(any(cn=='CompareT2')){
    if(pars['CompareT2']!='TRUE' & pars['CompareT2']!='FALSE') stop(paste0('In',opt$par,', CompareT2 value is not valid, it should be TRUE or FALSE!'),call. = FALSE)
    Par$CompareT2 <- as.logical(pars['CompareT2'])
}else{
    Par$CompareT2 <- NULL
}


if(any(cn=='LenRVmethod')){
    if(pars['LenRVmethod']!='PEXO' & pars['LenRVmethod']!='T2') stop(paste0('In',opt$par,', LenRVmethod value is not valid, it should be PEXO or T2!'),call. = FALSE)
    Par$LenRVmethod <- pars['LenRVmethod']
}else{
    Par$LenRVmethod <- 'T2'
}

## Binary model type: classical kepler, DDGR (default), DD
Par$binary <- FALSE
Par$BinaryModel <- 'none'
if(any(cn=='BinaryModel')){
    if(pars['BinaryModel']=='DDGR' | pars['BinaryModel']=='kepler'){
        Par$binary <- TRUE
        Par$BinaryModel <- pars['BinaryModel']
    }
}

## PPN parameter g (or gamma)
if(any(cn=='g')){
    Par$g <- as.numeric(pars['g'])
    if(is.na(Par$g)) stop('Error: g value is not numerical!',call.=FALSE)
}else{
    Par$g <- 1
}

##   BinaryUnit - Unit used in the computation of byinary orbit
##   unit 1: au and year (auyr)
##   unit 2: natural unit (c=1,[L]=[M]=[T]=s)
if(any(cn=='BinaryUnit')){
    Par$BinaryUnit <- pars['BinaryUnit']
}else{
    Par$BinaryUnit <- 'auyr'
}

###Read observertory/telescope parameters
if(any(cn=='ellipsoid')){
    Par$ellipsoid <- pars['ellipsoid']
}else{
    Par$ellipsoid <- 'WGS84'
}
if(Par$ellipsoid=='WGS84'){
    Par$n <- 1
}else if(Par$ellipsoid=='GRS80'){
    Par$n <- 2
}else if(Par$ellipsoid=='WGS72'){
    Par$n <- 3
}else{
    stop('Error: reference ellipsoid for Earth rotation model is not given!')
}

GetObsInd <- function(target,ObsList){
    ind <- grep(tolower(target),tolower(ObsList[,'Name']))
    if(length(ind)>1){
        Ns <- c()
        for(k in ind){
            s1 <- unlist(strsplit(as.character(ObsList[k,'Name']),split=''))
            s1 <- tolower(s1[s1!=''])
            s2 <- unlist(strsplit(target,split=''))
            s2 <- tolower(s2[s2!=''])
            Ns <- c(Ns,length(which(!is.na(match(s1,s2)))))
        }
        ind <- ind[which.max(Ns)]
    }
 return(ind)
}

###If GPS location is specified, using the given ones instead of the ones from the observatory data file
if(any(cn=='xtel') & any(cn=='ytel') & any(cn=='ztel')){
    Par$xtel <- as.numeric(pars['xtel'])/1e3#km
    Par$ytel <- as.numeric(pars['ytel'])/1e3
    Par$ztel <- as.numeric(pars['ztel'])/1e3
    if(is.na(Par$xtel) | is.na(Par$ytel) | is.na(Par$ztel)) stop('Error: xtel, ytel or ztel value is not numerical!',call.=FALSE)
    eph <- sofa_Gc2gd(2,xyz=as.numeric(pars[c('xtel','ytel','ztel')]))
    Par$elong <- eph['elong']#rad
    Par$phi <- eph['phi']#rad
    Par$height <- eph['height']/1e3#km
}else if(any(cn=='phi') & any(cn=='elong') & any(cn=='height')){
    Par$height <- as.numeric(pars['height'])#km
    Par$elong <- as.numeric(pars['elong'])*pi/180#rad
    Par$phi <- as.numeric(pars['phi'])*pi/180#rad
    if(is.na(Par$height) | is.na(Par$elong) | is.na(Par$phi)) stop('Error: height, elong or phi value is not numerical!',call.=FALSE)
    xyz <- as.numeric(sofa_Gd2gc(Par$n,elong=Par$elong,phi=Par$phi,height=Par$height*1e3))
    Par$xtel <- xyz[1]/1e3#km
    Par$ytel <- xyz[2]/1e3#
    Par$ztel <- xyz[3]/1e3
}
Par$ObsType <- 'ground'#Observatory type: ground (default) or space

if(!any(names(Par)=='xtel') | !any(names(Par)=='ytel') | !any(names(Par)=='ztel') ){
    if(any(cn=='ObsCode')){
        Par$ObsCode <- as.character(pars['ObsCode'])
    }
    obs <- read.csv2('../observatories/observatory_MPC.csv',row.names=NULL,quote='')
    MatchObs <- FALSE
    if(any(cn=='observatory' | cn=='ObsCode')){
###looking for ground or space-based observatory ephemeris...
        if(!any(names(Par)=='ObsCode')){
            Par$observatory <- pars['observatory']#could be code or observatory name which might not be matched well
            ind <- which(obs[,'Code']==Par$observatory)
            if(length(ind)==0){
                ind <- GetObsInd(Par$observatory,obs)
            }
        }else{
            ind <- which(as.character(obs[,'Code'])==Par$ObsCode)
        }
                                        #
        if(!is.na(obs[ind,'x.m']) & !is.na(obs[ind,'y.m']) & !is.na(obs[ind,'z.m'])){
            if(!all(obs[ind,c('x.m','y.m','z.m')]==0)){
                MatchObs <- TRUE
                Par$ObsCode <- as.character(obs[ind,'Code'])
###Ground-based observatory
                Par$xtel <- obs[ind,'x.m']/1e3#km
                Par$ytel <- obs[ind,'y.m']/1e3
                Par$ztel <- obs[ind,'z.m']/1e3
                Par$elong <- obs[ind,'Elong.rad']
                Par$phi <- obs[ind,'PhiGd.rad']
                Par$height <- obs[ind,'height.m']/1e3#km
            }
        }
        if(!MatchObs){
###Space-based observatory
            Par$ObsType <- 'space'
###https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html
                                        #        obslist <- read.table('../observatories/spacecraft_code.txt',header=TRUE)
            if(!any(names(Par)=='ObsCode')){
                obslist <- read.table('../observatories/satellite_list.csv',header=TRUE,sep=',')
                ind <- GetObsInd(Par$observatory,obslist)
            }else{
                ind <- which(as.character(obslist[,'Code'])==ObsCode)
            }
            if(length(ind)>0){
                MatchObs <- TRUE
                Par$ObsCode <- ObjId <- as.numeric(obslist[ind,'Code'])
                jd.utc <- rowSums(utc)
                Nmax <- 500
                if(Par$Nepoch<=Nmax){
                    system(paste('python GetSpaceObsEph.py', ObjId, opt$time))
                    SpaceObs <- read.csv('../observatories/space.csv',header=TRUE)
                }else{
                    Nf <- ceiling(Par$Nepoch/Nmax)
                    SpaceObs <- c()
                    for(jj in 1:Nf){
                        ftest <- paste0('../observatories/test',jj,'.txt')
                        write.table(jd.utc[((jj-1)*Nmax+1):min(jj*Nmax,Par$Nepoch)],file=ftest,quote=FALSE,row.names=FALSE,col.names=FALSE)
                        system(paste('source GetEph.sh', ObjId, ftest))
                                        #                    system(paste('rm',ftest))
                        tmp <- read.csv('../observatories/space.csv',header=TRUE)
                        SpaceObs <- rbind(SpaceObs,tmp)
                    }
                }
                SpaceObs[,1:3] <- SpaceObs[,1:3]*au2km#km
                SpaceObs[,4:6] <- SpaceObs[,4:6]*DJY*auyr2kms#km/s
                R <- gen_CalLen(SpaceObs[,1:3])
                lb <- gen_Xyz2lb(SpaceObs[,1:3])
                equ <- gen_Ecl2equ(lb[,1],lb[,2])
                SpaceObs[,1:3] <- R*cbind(equ$x,equ$y,equ$z)
                V <- gen_CalLen(SpaceObs[,4:6])
                lb <- gen_Xyz2lb(SpaceObs[,4:6])
                equ <- gen_Ecl2equ(lb[,1],lb[,2])
                SpaceObs[,4:6] <- V*cbind(equ$x,equ$y,equ$z)
                SpaceObs <- as.matrix(SpaceObs)
            }
        }
###warning if no observatory name is found
        if(!MatchObs) cat('Warning: observatory name or code is not found in the MPC and JPL space observatory file. You may update ../observatories/spacecraft_code.txt or ../observatories/observatory_MPC.txt by adding new observatorie!\n')
    }
}
if(!any(names(Par)=='xtel') & !exists('SpaceObs')) stop('Error: Observatory data is not given!')

###Read atmospheric parameters
Par$hm <- Par$height*1e3#meter, used by routine refro()
if(any(cn=='tdk')){
    Par$tdk <- as.numeric(pars['tdk'])
    if(is.na(Par$tdk)) stop('Error: tdk value is not numerical!',call.=FALSE)
    if(all(is.na(Par$tdk))) Par$tdk <- as.numeric(read.table(pars['tdk']))
}else{
    Par$tdk <- 278
    cat('Warning: An ambient temperature at the observer (tdk) is not provided and tdk=278K is adopted!\n')
}

if(any(cn=='pmb')){
    Par$pmb <- as.numeric(pars['pmb'])#millibar
    if(is.na(Par$pmb)) stop('Error: pmb value is not numerical!',call.=FALSE)
    if(all(is.na(Par$pmb))) Par$pmb <- as.numeric(read.table(pars['pmb']))
}else{
    Par$pmb <- 1013.25
    cat('Warning: pressure at the observer (pmb) is not provided and tdk=1013.25 millibar is adopted!\n')
}

if(any(cn=='rh')){
    Par$rh <- as.numeric(pars['rh'])#
    if(is.na(Par$rh)) stop('Error: rh value is not numerical!',call.=FALSE)
    if(all(is.na(Par$rh))) Par$rh <- as.numeric(unlist(read.table(pars['rh'])))
}else{
    Par$rh <- 0.1
    cat('Warning: relative humidity at the observer (rh) is not provided and rh=0.1 is adopted!\n')
}

if(any(cn=='wl')){
    Par$wl <- as.numeric(pars['wl'])#micrometre
    if(is.na(Par$wl)) stop('Error: wl value is not numerical!',call.=FALSE)
}else{
    Par$wl <- 0.5
    cat('Warning: effective wavelength of the source (wl) is not provided and wl=0.5 micrometre is adopted!\n')
}

if(any(cn=='tlr')){
    Par$tlr <- as.numeric(pars['tlr'])#K/metre
    if(is.na(Par$tlr)) stop('Error: tlr value is not numerical!',call.=FALSE)
}else{
    Par$tlr <- 0.0065
}

###Read astrometry parameters
try(if(!any(cn=='epoch')) stop('Error: astrometry epoch is not given!'))
Par$tpos <- time_ToJD(as.numeric(pars['epoch']))#reference astrometry epoch
if(is.na(Par$tpos)) stop('Error: epoch value is not numerical!',call.=FALSE)

try(if(!any(cn=='ra')) stop('Error: ra is not given!'))
ra0 <- Par$ra <- as.numeric(pars['ra'])/180*pi#rad
if(is.na(Par$ra)) stop('Error: ra value is not numerical!',call.=FALSE)

try(if(!any(cn=='dec')) stop('Error: dec is not given!'))
dec0 <- Par$dec <- as.numeric(pars['dec'])/180*pi#rad
if(is.na(Par$dec)) stop('Error: dec value is not numerical!',call.=FALSE)

Par$p <- as.numeric(c(-sin(ra0),cos(ra0),0))#or hat{alpha} in tempo2
Par$q <- as.numeric(c(-sin(dec0)*cos(ra0),-sin(dec0)*sin(ra0),cos(dec0)))#or hat{delta} in tempo2
Par$u <- as.numeric(c(cos(dec0)*cos(ra0),cos(dec0)*sin(ra0),sin(dec0)))# or hat{R0} in tempo2

if(any(cn=='plx')){
    Par$plx <- as.numeric(pars['plx'])#mas
    if(is.na(Par$plx)) stop('Error: plx value is not numerical!',call.=FALSE)
}else if(any(cn=='dist')){
    Par$plx <- 1000/dist#dist in pc unit
}else{
    Par$plx <- 1e-6#equal to D=1 Mpc so that parallex can be ignored if not given.
}
Par$dist <- 1000/Par$plx
if(Par$dist>1){
    Par$near <- FALSE
}else{
    Par$near <- TRUE
}

if(any(cn=='pmra')){
    Par$pmra <- as.numeric(pars['pmra'])#mas/yr
    if(is.na(Par$pmra)) stop('Error: pmra value is not numerical!',call.=FALSE)
}else{
    Par$pmra <- 0
}

if(any(cn=='pmdec')){
    Par$pmdec <- as.numeric(pars['pmdec'])#mas/yr
    if(is.na(Par$pmdec)) stop('Error: pmdec value is not numerical!',call.=FALSE)
}else{
    Par$pmdec <- 0
}

if(any(cn=='rv')){
    Par$rv <- as.numeric(pars['rv'])#km/s
    if(is.na(Par$rv)) stop('Error: rv value is not numerical!',call.=FALSE)
}else{
    Par$rv <- 0
}

###Read Keplerian parameters, and if they are given the binary model should be used
if(any(cn=='mT')){
    Par$mT <- Par$mTC <- as.numeric(pars['mT'])#Msun
    if(is.na(Par$mT)) stop('Error: mT value is not numerical!',call.=FALSE)
}else{
    Par$mTC <- Par$mT <- 1
}
Par$aT <- Par$aTC <- Par$aC <- Par$P <- Par$e <- Par$I <- Par$omegaT <- Par$omegaC <- Par$Omega <- Par$Tp <- NA
Par$mC <- 0
if(Par$binary){
    try(if(!any(cn=='mC')) stop('Error: companion mass mC is not given!\n'))
    Par$mC <- as.numeric(pars['mC'])
    if(is.na(Par$mC)) stop('Error: mC value is not numerical!',call.=FALSE)
    Par$mTC <- Par$mTC+Par$mC
    if(any(cn=='P')){
        Par$P <- as.numeric(pars['P'])#year
        if(is.na(Par$P)) stop('Error: P value is not numerical!',call.=FALSE)
        Par$aTC <- (Par$P^2*Par$mTC)^(1/3)#au
        Par$aT <- Par$aTC*Par$mC/Par$mTC
        Par$aC <- Par$aTC*Par$mT/Par$mTC
    }else if(any(cn=='aC') | any(cn=='aT') | any(cn=='aTC')){
        if(any(cn=='aTC')){
            Par$aTC <- as.numeric(pars['aTC'])
            if(is.na(Par$aTC)) stop('Error: aTC value is not numerical!',call.=FALSE)
            Par$aT <- as.numeric(pars['aTC'])*Par$mC/Par$mTC
            Par$aC <- as.numeric(pars['aTC'])*Par$mT/Par$mTC
        }
        if(any(cn=='aT')){
            Par$aT <- as.numeric(pars['aT'])
            if(is.na(Par$aT)) stop('Error: aT value is not numerical!',call.=FALSE)
            Par$aC <- as.numeric(pars['aT'])/Par$mC*Par$mT
            Par$aTC <- as.numeric(pars['aT'])/Par$mC*Par$mTC
        }
        if(any(cn=='aC')){
            Par$aC <- as.numeric(pars['aC'])
            if(is.na(Par$aC)) stop('Error: aC value is not numerical!',call.=FALSE)
            Par$aT <- as.numeric(pars['aC'])/Par$mT*Par$mC
            Par$aTC <- as.numeric(pars['aC'])/Par$mT*Par$mTC
        }
        Par$P <- sqrt(Par$aTC^3/Par$mTC)#year
    }
    if(!any(names(Par)=='P')) stop('Error: orbital period P or semi-major axis is not given!\n')
    Par$Pd <- Par$P*DJY#day

    try(if(!any(cn=='e')) stop('Error: eccentricity e is not given!\n'))
    Par$e <- as.numeric(pars['e'])
    if(is.na(Par$e)) stop('Error: e value is not numerical!',call.=FALSE)


    try(if(!any(cn=='I')) stop('Error: Inclination I is not given!\n'))
    Par$I <- as.numeric(pars['I'])/180*pi#rad
    if(is.na(Par$I)) stop('Error: I value is not numerical!',call.=FALSE)

    try(if(!any(cn=='omegaC') & !any(cn=='omegaT')) stop('Error: argument of periastron omegaT or omegaC is not given!\n'))
    if(any(cn=='omegaT')){
        Par$omegaT <- as.numeric(pars['omegaT'])/180*pi#rad
        if(is.na(Par$omegaT)) stop('Error: omegaT value is not numerical!',call.=FALSE)
        Par$omegaC <- (as.numeric(pars['omegaT'])/180*pi+pi)%%(2*pi)
    }
    if(any(cn=='omegaC')){
        Par$omegaC <- as.numeric(pars['omegaC'])/180*pi#rad
        if(is.na(Par$omegaC)) stop('Error: omegaC value is not numerical!',call.=FALSE)
        Par$omegaT <- (as.numeric(pars['omegaC'])/180*pi+pi)%%(2*pi)
    }

    try(if(!any(cn=='Omega')) stop('Error: longitude of ascending node Omega is not given!\n'))
    Par$Omega <- as.numeric(pars['Omega'])/180*pi#rad
    if(is.na(Par$Omega)) stop('Error: Omega value is not numerical!',call.=FALSE)

    try(if(!any(cn=='Tp') & !any(cn=='Tc') & !any(cn=='M0')) stop('Error: mean or true anomaly or reference epoch at periastron or primary transit (i.e. M0 or Tp or Tc) is not given!\n'))

    if(any(cn=='T0')){
        Par$T0 <- time_ToJD(as.numeric(pars['T0']))#jd
        if(is.na(Par$T0)) stop('Error: T0 value is not numerical!',call.=FALSE)
    }else{
        Par$T0 <- rowSums(utc)[1]#first epoch of the utc time sequence
    }

    if(any(cn=='Tasc')){
        Par$Tasc <- as.numeric(pars['Tasc'])#jd; sometimes used in ELL1 model (Wex 1998)
        if(is.na(Par$Tasc)) stop('Error: Tasc value is not numerical!',call.=FALSE)
    }else{
        Par$Tasc <- NA#first epoch of the utc time sequence
    }

    if(any(cn=='Tc')){
        Par$Tc <- time_ToJD(as.numeric(pars['Tc']))
        if(is.na(Par$Tc)) stop('Error: Tc value is not numerical!',call.=FALSE)
        Par$Tp <- gen_Tc2tp(Par$e,Par$omegaT,Par$Pb,Par$Tc)
        Par$M0 <- (Par$Tp-Par$T0)/Par$Pb*2*pi
    }else if(any(cn=='Tp')){
        Par$Tp <- time_ToJD(as.numeric(pars['Tp']))
        if(is.na(Par$Tp)) stop('Error: Tp value is not numerical!',call.=FALSE)
        Par$Tc <- gen_CalTc(Par$Tp,Par$Pb,Par$e,Par$omegaT)
        Par$M0 <- (Par$Tp-Par$T0)/Par$Pb*2*pi
    }else if(any(cn=='M0')){
        Par$M0 <- as.numeric(pars['M0'])/180*pi#rad
        if(is.na(Par$M0)) stop('Error: M0 value is not numerical!',call.=FALSE)
        Par$Tp <- T0+M0/(2*pi)*Par$Pb#epoch at periastron
        Tc <- gen_CalTc(Par$Tp,Par$Pd,Par$e,Par$omegaT)#mid-transit epoch
    }
    Par$DD <- Par$DDGR <- list()
###Some default binary parameters for DD and DDGR
    Par$DDGR$xdot <- Par$DD$xdot <- 0
    Par$DDGR$edot <- Par$DD$edot <- 0
    Par$DDGR$omdot <- Par$DD$omdot <- 0
    Par$DDGR$pbdot <- Par$DD$pbdot <- 0
    Par$DDGR$xpbdot <- Par$DD$xpbdot <- 0
    Par$DDGR$wdot <- Par$DD$wdot <- 0
    Par$DD$dth <- 0
    Par$DD$gamma <- 0
    Par$DD$dr <- 0
    Par$DD$sini <- sin(Par$I)

###binary/companion parameters
    Par$DDGR$a1 <- as.numeric(Par$mC/Par$mTC*(Par$P^2*Par$mTC)^(1/3)*Par$DD$sini)#au
    Par$DDGR$x0 <- Par$DDGR$a1/YC#yr
    Par$DDGR$omz <- Par$omegaT
    Par$DDGR$pb <- Par$P
    Par$DD$ecc <- Par$DDGR$ecc <- Par$e
}

####load ephermies and EOP and other input data sets
source('load_data.R')
