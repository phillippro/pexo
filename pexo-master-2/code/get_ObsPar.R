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
                NepochMax <- 500
                if(Par$Nepoch<=NepochMax){
                    system(paste('python GetSpaceObsEph.py', ObjId, opt$time))
                    SpaceObs <- read.csv('../observatories/space.csv',header=TRUE)
                }else{
                    Nf <- ceiling(Par$Nepoch/NepochMax)
                    SpaceObs <- c()
                    for(jj in 1:Nf){
                        ftest <- paste0('../observatories/test',jj,'.txt')
                        write.table(jd.utc[((jj-1)*NepochMax+1):min(jj*NepochMax,Par$Nepoch)],file=ftest,quote=FALSE,row.names=FALSE,col.names=FALSE)
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
            }else if(grepl('hip|HIP',Par$observatory)){
                tab <- read.table('../observatories/hipparcos_state.txt')
                xfunc <- approxfun(tab[,1],tab[,2])
                yfunc <- approxfun(tab[,1],tab[,3])
                zfunc <- approxfun(tab[,1],tab[,4])
                vxfunc <- approxfun(tab[,1],tab[,5])
                vyfunc <- approxfun(tab[,1],tab[,6])
                vzfunc <- approxfun(tab[,1],tab[,7])
                SpaceObs <- cbind(xfunc(jd.utc),yfunc(jd.utc),zfunc(jd.utc))
            }
        }
###warning if no observatory name is found
        if(!MatchObs & opt$verbose) cat('Warning: observatory name or code is not found in the MPC and JPL space observatory file. You may update ../observatories/spacecraft_code.txt or ../observatories/observatory_MPC.txt by adding new observatorie!\n')
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
    if(opt$verbose) cat('Warning: An ambient temperature at the observer (tdk) is not provided and tdk=278K is adopted!\n')
}

if(any(cn=='pmb')){
    Par$pmb <- as.numeric(pars['pmb'])#millibar
    if(is.na(Par$pmb)) stop('Error: pmb value is not numerical!',call.=FALSE)
    if(all(is.na(Par$pmb))) Par$pmb <- as.numeric(read.table(pars['pmb']))
}else{
    Par$pmb <- 1013.25
    if(opt$verbose) cat('Warning: pressure at the observer (pmb) is not provided and tdk=1013.25 millibar is adopted!\n')
}

if(any(cn=='rh')){
    Par$rh <- as.numeric(pars['rh'])#
    if(is.na(Par$rh)) stop('Error: rh value is not numerical!',call.=FALSE)
    if(all(is.na(Par$rh))) Par$rh <- as.numeric(unlist(read.table(pars['rh'])))
}else{
    Par$rh <- 0.1
    if(opt$verbose) cat('Warning: relative humidity at the observer (rh) is not provided and rh=0.1 is adopted!\n')
}

if(any(cn=='wl')){
    Par$wl <- as.numeric(pars['wl'])#micrometre
    if(is.na(Par$wl)) stop('Error: wl value is not numerical!',call.=FALSE)
}else{
    Par$wl <- 0.5
    if(opt$verbose) cat('Warning: effective wavelength of the source (wl) is not provided and wl=0.5 micrometre is adopted!\n')
}

if(any(cn=='tlr')){
    Par$tlr <- as.numeric(pars['tlr'])#K/metre
    if(is.na(Par$tlr)) stop('Error: tlr value is not numerical!',call.=FALSE)
}else{
    Par$tlr <- 0.0065
}
