obscor <- c()
xn <- paste0('xtel.',ins)
yn <- paste0('ytel.',ins)
zn <- paste0('ztel.',ins)
en <- paste0('elong.',ins)
pn <- paste0('phi.',ins)
hn <- paste0('height.',ins)
tn <- paste0('ObsType.',ins)

obsn <- paste0('observatory.',ins)
coden <- paste0('ObsCode.',ins)
sn <- paste0('SpaceObs.',ins)

#brn <-paste0('bRv.',star,'.',ins)
#arn <-paste0('aRv.',star)
#bran <-paste0('bAstro.',star,'.ra',ins)
#bdecn <-paste0('bAstro.',star,'.dec',ins)

if(any(cn=='xtel')) Par[[xn]] <- Par$xtel
if(any(cn=='ytel')) Par[[yn]] <- Par$ytel
if(any(cn=='ztel')) Par[[zn]] <- Par$ztel
#if(any(cn=='bRv')) Par[[zn]] <- Par$ztel
if(any(cn=='elong')) Par[[en]] <- Par$elong
if(any(cn=='phi')) Par[[pn]] <- Par$phi
if(any(cn=='height')) Par[[hn]] <- Par$height
if(any(cn=='ObsType')) Par[[tn]] <- Par$ObsType
if(any(cn=='observatory')){
    Par[[obsn]] <- Par$observatory
}else{
    Par[[obsn]] <- ins
}
cn <- names(Par)
if(any(cn=='ObsCode')) Par[[coden]] <- Par$ObsCode
if(any(cn=='SpaceObs')) Par[[sn]] <- Par$SpaceObs
Par[[tn]] <- 'ground'#Observatory type: ground (default) or space
if(any(cn==xn) & any(cn==yn) & any(cn==zn)){
####customerized gound based observatory data
    eph <- sofa_Gc2gd(2,xyz=c(Par[[xn]]/1e3,Par[[yn]]/1e3,Par[[zn]]/1e3))
    Par[[en]] <- eph['elong']*180/pi#deg
    Par[[pn]] <- eph['phi']*180/pi#deg
    Par[[hn]] <- eph['height']/1e3#km
}else if(any(cn==en) & any(cn==pn) & any(cn==hn)){
####customerized gound based observatory data
    xyz <- as.numeric(sofa_Gd2gc(Par$n,elong=Par[[en]]*pi/180,phi=Par[[pn]]*pi/180,height=Par[[hn]]*1e3))
    Par[[xn]] <- xyz[1]/1e3#km
    Par[[yn]] <- xyz[2]/1e3#
    Par[[zn]] <- xyz[3]/1e3
}else{
####Automatically search for ground-based observatory data
    obs <- read.csv2('../observatories/observatory_MPC.csv',row.names=NULL,quote='')
    MatchObs <- FALSE
    if(any(cn==obsn | cn==coden)){
#    if(!grepl('hip|Hip',obsn) & !grepl('hip|Hip',coden)){
###looking for ground or space-based observatory ephemeris...
        ind <- c()
        if(any(cn==obsn)){
            Par$observatory <- Par[[obsn]]#could be code or observatory name which might not be matched well
            ind <- gen_GetObsInd(Par[[obsn]],obs)
        }else if(any(cn==coden)){
            ind <- which(as.character(obs[,'Code'])==Par[[coden]])
        }
        if(length(ind)>0){
            if(!is.na(obs[ind,'x.m']) & !is.na(obs[ind,'y.m']) & !is.na(obs[ind,'z.m'])){
                if(!all(obs[ind,c('x.m','y.m','z.m')]==0)){
                    MatchObs <- TRUE
                    Par[[coden]] <- as.character(obs[ind,'Code'])
                    Par[[obsn]] <- as.character(obs[ind,'Name'])
###Ground-based observatory
                    Par[[xn]] <- obs[ind,'x.m']/1e3#km
                    Par[[yn]] <- obs[ind,'y.m']/1e3
                    Par[[zn]] <- obs[ind,'z.m']/1e3
                    Par[[en]] <- obs[ind,'Elong.rad']
                    Par[[pn]] <- obs[ind,'PhiGd.rad']
                    Par[[hn]] <- obs[ind,'height.m']/1e3#km
                }
            }
        }

        if(!MatchObs){
####Automatically search for space-based observatory data
            Par[[tn]] <- 'space'
###https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html
                                        #        obslist <- read.table('../observatories/spacecraft_code.txt',header=TRUE)
            obslist <- read.table('../observatories/satellite_list.csv',header=TRUE,sep=',')
            if(any(cn==obsn)){
                ind <- gen_GetObsInd(Par[[obsn]],obslist)
            }else if(any(cn==coden)){
                ind <- which(as.character(obslist[,'Code'])==Par[[coden]])
            }
            if(length(ind)>0){
                MatchObs <- TRUE
                Par[[coden]] <- ObjId <- as.numeric(obslist[ind,'Code'])
                NepochMax <- 500
                if(Par$Nepoch<=NepochMax){
###single python process for small data
                    write.table(rowSums(utc),file='test.tim',quote=FALSE,row.names=FALSE,col.names=FALSE)
                    system(paste('python GetSpaceObsEph.py', ObjId, 'test.tim'))
                    SpaceObs <- read.csv('../observatories/space.csv',header=TRUE)
                    system('rm ../observatories/space.csv')
                    system('rm test.tim')
                }else{
###loop for big data
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
               Par[[sn]] <- SpaceObs#for space telescope
            }
        }
###warning if no observatory name is found
        if(!MatchObs) cat('Warning: observatory name or code is not found in the MPC and JPL space observatory file. You may update ../observatories/spacecraft_code.txt or ../observatories/observatory_MPC.txt by adding new observatorie!\n')
    }else{
        Par[[tn]] <- 'space'
        tab <- read.table('../observatories/hipparcos_state.txt')
        xfunc <- approxfun(tab[,1],tab[,2])
        yfunc <- approxfun(tab[,1],tab[,3])
        zfunc <- approxfun(tab[,1],tab[,4])
        vxfunc <- approxfun(tab[,1],tab[,5])
        vyfunc <- approxfun(tab[,1],tab[,6])
        vzfunc <- approxfun(tab[,1],tab[,7])
        jd.utc <- rowSums(utc)
        t <- jd.utc[index]
        SpaceObs <- cbind(xfunc(t),yfunc(t),zfunc(t),vxfunc(t),vyfunc(t),vzfunc(t))
        SpaceObs <- as.matrix(SpaceObs)
        Par[[sn]] <- SpaceObs
    }
}

if(!any(names(Par)==xn) & !any(sn==names(Par))) stop('Error: Observatory data is not given!')

if(Par[[tn]]=='ground'){
    xtel <- rep(Par[[xn]],length(index))
    ytel <- rep(Par[[yn]],length(index))
    ztel <- rep(Par[[zn]],length(index))
    elong <- rep(Par[[en]],length(index))
    phi <- rep(Par[[pn]],length(index))
    height <- rep(Par[[hn]],length(index))
    obscor <- cbind(xtel,ytel,ztel,elong,phi,height,Par[[tn]])#for ground telescope,
}else{
    obscor <- cbind(SpaceObs,Par[[tn]])
}

###Read atmospheric parameters
Par[[paste0('hm.',ins)]] <- Par[[hn]]*1e3#meter, used by routine refro()
ref <- paste0('RefType.',ins)
tdk <- paste0('tdk.',ins)
pmb <- paste0('pmb.',ins)
rh <- paste0('rh.',ins)
wl <- paste0('wl.',ins)
tlr <- paste0('tlr.',ins)
p <-  paste0('p.',ins)
q <-  paste0('q.',ins)
if(any(cn=='RefType')) Par[[ref]] <- Par$RefType
if(any(cn=='tdk')) Par[[tdk]] <- Par$tdk
if(any(cn=='pmb')) Par[[pmb]] <- Par$pmb
if(any(cn=='rh')) Par[[rh]] <- Par$rh
if(any(cn=='wl')) Par[[wl]] <- Par$wl
if(any(cn=='tlr')) Par[[tlr]] <- Par$tlr
if(any(cn=='p')) Par[[p]] <- Par$p
if(any(cn=='q')) Par[[q]] <- Par$q

nn <- c(ref,tdk,pmb,rh,wl,tlr,p,q)
dval <- c('none',278,1013.25,0.1,0.5,0.0065,0,0)
for(k in 1:length(nn)){
    n <-  nn[k]
    if(any(cn==n)){
        if(all(is.na(Par[[n]]))){
            if(!file.exists(Par[[n]])) stop('Error: ',n,' value is not valid!')
            Par[[n]] <- as.numeric(read.table(Par[[n]]))
        }
    }else{
        Par[[n]] <- dval[k]
        if(opt$verbose) cat('Warning:',n,'is not provided and ',dval[k],' is adopted!\n')
    }
    obscor <- cbind(obscor,Par[[n]])
}
