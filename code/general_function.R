########################################################################
###This file contains general-purpose functions which are called by
###timing, astrometry and RV functions.
########################################################################
gen_RvzItrs <- function(robs.itrs,EopPar,Nt,sp=NULL){
####################################
##Given the ITRS position of a telescope, calculate its ITRS velocity
##
## Input:
##   robs.itrs - ITRS position (km)
##   EopPar - Earth orientation parameters (EOP)
##   Nt - Number of epochs
##   sp - Polar motion parameter
##
## Output:
##   vobs.itrs - Observatory velocity in ITRS (km/s)
##   pole - ITRS unit vector of the Pole
##   rpom - Polar motion matrix
####################################
    if(is.null(sp)) sp <- rep(0,Nt)
#    rpom <- sofa_Pom00(EopPar$xp,EopPar$yp,sp=rep(0,length(EopPar$xp)))#polar motion matrix
    rpom <- sofa_Pom00(EopPar$xp,EopPar$yp,sp=sp)#polar motion matrix
    north <- t(t(c(0,0,1)))#Vector to +ve pole
    pole.itrs <- t(sapply(1:length(EopPar$xp),function(k) rpom[,,k]%*%north))#Spin pole in ITRS
    eradot <- 2.0*pi*1.00273781191135448*(1.0+EopPar$dut1dot)/DAYSEC#rad/s
    omega.itrs <- pole.itrs*eradot#Angular velocity in ITRS (rad/s)
    vobs.itrs <- t(sapply(1:Nt, function(k) cross(omega.itrs[k,],robs.itrs)))#km/s
    list(robs=robs.itrs,vobs=vobs.itrs,pole=pole.itrs,rpom=rpom)
}

gen_PMecl2equ <- function(lambda,beta,PMlambda,PMbeta,dt=1){
####################################
##Given proper motions (PM) in equatorial coordinates, find the PM in ecliptic coordinates
##
## Input:
##   lambda - ecliptic longitude (deg)
##   beta - ecliptic latitude (deg)
##   PMlambda - PM in lambda (or ecliptic longitude direction; mas/yr)
##   PMbeta - PM in beta (or ecliptic latitude direction; mas/yr)
##   dt - the time increasement used for numerical computation (yr)
##
## Output:
##   PMalpha - PM in alpha (or equatorial longitude direction; mas/yr)
##   PMdelta - PM in delta (or equatorial longitude direction; mas/yr)
##
## Note:
##   The PM in longtiude (PMx) is acturally PMx* following Hipparcos' and Gaia's conventions. i.e. PMx=PMalpha*cos(delta) or PMlambda*cos(beta); Thus dx=PMx/cos(y)
####################################
    lambda <- lambda/180*pi
    beta <- beta/180*pi
    dlambda <- PMlambda*dt*DMAS2R/cos(beta)
    lambda1 <- lambda+dlambda
    lambda0 <- lambda-dlambda
    dbeta <- PMbeta*dt*DMAS2R
    beta1 <- beta+dbeta
    beta0 <- beta-dbeta
    tmp1 <- gen_Ecl2equ(lambda1,beta1)
    tmp0 <- gen_Ecl2equ(lambda0,beta0)
    tmp <- gen_Ecl2equ(lambda,beta)
    alpha <- tmp$alpha*180/pi
    delta <- tmp$delta*180/pi
    PMalpha <- (tmp1$alpha-tmp0$alpha)/DMAS2R/(2*dt)*cos(tmp$delta)
    PMdelta <- (tmp1$delta-tmp0$delta)/DMAS2R/(2*dt)
    return(list(alpha=alpha,delta=delta,PMalpha=PMalpha,PMdelta=PMdelta))
}

gen_PMequ2ecl <- function(alpha,delta,PMalpha,PMdelta,dt=1){
####################################
##Given proper motions (PM) in equatorial coordinates, find the PM in ecliptic coordinates
##
## Input:
##   alpha - equatorial longitude (or ascending node; deg)
##   delta - equatorial latitude (deg)
##   PMalpha - PM in alpha (or equatorial longitude direction; mas/yr)
##   PMdelta - PM in delta (or equatorial longitude direction; mas/yr)
##   dt - the time increasement used for numerical computation (yr)
##
## Output:
##   PMlambda - PM in lambda (or ecliptic longitude direction; mas/yr)
##   PMbeta - PM in beta (or ecliptic latitude direction; mas/yr)
##
## Note:
##   The PM in longtiude (PMx) is acturally PMx* following Hipparcos' and Gaia's conventions. i.e. PMx=PMalpha*cos(delta) or PMlambda*cos(beta); Thus dx=PMx/cos(y)
####################################
    alpha <- alpha/180*pi
    delta <- delta/180*pi
    dalpha <- PMalpha*dt*DMAS2R/cos(delta)
    alpha1 <- alpha+dalpha
    alpha0 <- alpha-dalpha
    ddelta <- PMdelta*dt*DMAS2R
    delta1 <- delta+ddelta
    delta0 <- delta-ddelta
    tmp1 <- gen_Equ2ecl(alpha1,delta1)
    tmp0 <- gen_Equ2ecl(alpha0,delta0)
    tmp <- gen_Equ2ecl(alpha,delta)
    lambda <- tmp$lambda*180/pi
    beta <- tmp$beta*180/pi
    PMlambda <- (tmp1$lambda-tmp0$lambda)/DMAS2R/(2*dt)*cos(tmp$beta)
    PMbeta <- (tmp1$beta-tmp0$beta)/DMAS2R/(2*dt)
    return(list(lambda=lambda,beta=beta,PMlambda=PMlambda,PMbeta=PMbeta))
}

gen_Equ2ecl <- function(alpha, delta,tt=NULL){
####################################
## Convert the Equatorial coordiantes into ecliptic coordinates
##
## Input:
##   alpha -  RA (rad)
##   delta -  DEC (rad)
##   tt - TT (JD)
##
## Output:
##   x - x component of unit vector in the ecliptic coordinate system
##   y - y component of unit vector in the ecliptic coordinate system
##   z - z component of unit vector in the ecliptic coordinate system
##   lambda - Longitude in the ecliptic coordinate system (rad)
##   beta - Latitude in the ecliptic coordinate system (rad)
####################################
#    eps <- 23.439281/180*pi
#    eps <- 84381.448/206264.8062471
    if(!is.null(tt)){
        eps <- sofa_Obl80(tt)
    }else{
        eps <- 84381.448/pc2au
    }
    x <- cos(alpha)*cos(delta)
    y <- sin(eps)*sin(delta)+sin(alpha)*cos(delta)*cos(eps)
    z <- cos(eps)*sin(delta)-sin(alpha)*cos(delta)*sin(eps)
    lb <- gen_Xyz2lb(cbind(x,y,z))
    return(list(x=x,y=y,z=z,lambda=lb[,1],beta=lb[,2]))
}

gen_Ecl2equ <- function(lambda, beta, tt=NULL){
####################################
## Convert the Equatorial coordiantes into ecliptic coordinates
##
## Input:
##   lambda -  Longitude (rad)
##   delta -  Latitude (rad)
##   tt - TT (JD)
##
## Output:
##   x - x component of unit vector in the equatorial coordinate system
##   y - y component of unit vector in the equatorial coordinate system
##   z - z component of unit vector in the equatorial coordinate system
##   alpha - Longitude in the equatorial coordinate system (rad)
##   delta - Latitude in the equatorial coordinate system (rad)
####################################
#    eps <- 23.439281/180*pi
#    eps <- 84381.448/206265
    if(!is.null(tt)){
        eps <- sofa_Obl80(tt)
    }else{
        eps <- 84381.448/pc2au
    }
    x <- cos(lambda)*cos(beta)
    y <- cos(eps)*sin(lambda)*cos(beta)-sin(eps)*sin(beta)
    z <- sin(eps)*sin(lambda)*cos(beta)+cos(eps)*sin(beta)
    ad <- gen_Xyz2lb(cbind(x,y,z))
    return(list(x=x,y=y,z=z,alpha=ad[,1],delta=ad[,2]))
}

gen_CalLen <- function(x){
####################################
## Calculate the length or magnitude of a given vector
##
## Input:
##   x   A vector
##
## Output:
##   Unit vector of x
####################################
    sqrt(rowSums(x^2))
}
gen_CalOffset <- function(u1,u2,bref=NULL){
####################################
## Calculate the angular offset between two unit vectors
##
## Input:
##   u1 - Unit vector 1
##   u2 -  Unit vector 2
##   bref - the reference latitude; if null the mean of u1 and u2 delta would be used
##
## Output:
##   dl - Longitude offset (rad)
##   db - Latitude offset (rad)
####################################
    lb1 <- gen_Xyz2lb(u1)
    lb2 <- gen_Xyz2lb(u2)
    if(is.null(bref)){
        dl <-(lb1[,1]-lb2[,1])*cos((lb1[,2]+lb2[,2])/2)
    }else{
        dl <-(lb1[,1]-lb2[,1])*cos(bref)
    }
    db <- lb1[,2]-lb2[,2]
    cbind(dl=dl*pc2au,db=db*pc2au)#arcsec
}
gen_CalUnit <- function(x){
####################################
## Calculate the unit vector for a given vector
##
## Input:
##   x   A vector
##
## Output:
##   Unit vector of x
####################################
    x/gen_CalLen(x)
}
gen_Xy2phi <- function(x,y){
####################################
## Calculate the x,y coordinates into polar coordinates
##
## Input:
##   x - x coordinate
##   y - y coordinate
##
## Output:
##   phi - Counter-clockwise angle of (x,y) from the x axis
####################################
    phi <- atan(y/x)
    inds <- which(x<0)
    phi[inds] <- phi[inds]+pi
    inds <- which(x>=0 & y<0)
    phi[inds] <- phi[inds]+2*pi
    return(phi)
}
gen_CalAmp <- function(x){
####################################
## Calculate the peak-to-peak amplitude of the variation in a time series
##
## Input:
##   x   Time series
##
## Output:
##   Peak to peak magnitude of x
####################################
    max(x)-min(x)
}

gen_GetEMRAT <- function(DE){
####################################
## Get Earth-Moon mass ratio from a given JPL ephemeris
##
## Input:
##   DE   JPL ephemeris version
##
## Output:
##   EMRAT Earth-Moon mass ratio
####################################
    str <- DEheader
    ind <- grep('GROUP',str)
    Ngroup <- c()
    for(j in ind){
        tmp <- unlist(strsplit(str[j],' '))
        tmp <- tmp[tmp!='']
        Ngroup <- c(Ngroup,as.integer(tmp[2]))
    }
    ind1 <- ind[Ngroup==1040]
    ind2 <- ind[Ngroup==1041]
    tmp <- paste(str[(ind1+3):(ind2-2)],collapse=' ')
    tmp2 <- unlist(strsplit(tmp,' '))
    names <- tmp2[tmp2!='']
    ind.EMRAT <- which(names=='EMRAT')

    ind1 <- ind[Ngroup==1041]
    ind2 <- ind[Ngroup==1050]
    tmp <- paste(str[(ind1+3):(ind2-2)],sep=' ')
    tmp2 <- unlist(strsplit(tmp,' '))
    vals <- as.numeric(gsub('D','e',tmp2[tmp2!='']))
    vals[ind.EMRAT]
}

gen_CalEph <- function(tdb,body='Earth-Moon barycenter',DE=430){
####################################
## For each TDB, find the position and velocity as well as TT-TDB for a planet.
## Note that the geocenter's ephemeris is not available in JPL ephemerides.
## One need to get Earth-Moon barycenter and the Moon ephemerides to derive
## the geocenter's ephemeris. Note that DExxxt has the same position and velocity
## ephemerides as DExxx but with additional time ephemeris (or TT-TDB as a function of TDB).
## Since this function is only
## Input:
##   tdb    TDB time
##   body   A solar system object
##   DE     Version of JPL dynamical ephemerides (DE); the default is DE430t
##
## Output:
##   barycentric position and
####################################
    Nt <- nrow(tdb)
    quantities <- c('Mercury','Venus','Earth-Moon barycenter','Mars','Jupiter','Saturn','Uranus','Neptune','Pluto','Moon (geocentric)','Sun','Earth Nutations in longitude and obliquity (IAU 1980 model)','Lunar mantle libration','Lunar mantle angular velocity','TT-TDB (at geocenter)')
    ind.body <- which(quantities==body)
    Ndim <- 3
    if(body=='Earth Nutations in longitude and obliquity (IAU 1980 model)'){
        Ndim <- 2
    }else if(body=='TT-TDB (at geocenter)'){
        Ndim <- 1
    }
    str <- DEheader
###find the ephemeris of the Earth
##ref: https://space.stackexchange.com/questions/12506/what-is-the-exact-format-of-the-jpl-ephemeris-files
    ind <- grep('GROUP   1050',str)
    tmp3 <- c()
    for(i in 1:3){
        tmp <- str[ind+1+i]
        tmp2 <- unlist(strsplit(tmp,' '))
        tmp3 <- rbind(tmp3,as.integer(tmp2[tmp2!='']))
    }
    par <- tmp3[,ind.body]
    Nsub <- par[3]
    Ncheb <- par[2]
    ind.start <- par[1]

###barycenter of the Earth-moon center
    chebs <- chebyshev.t.polynomials(Ncheb-1)
    cheb.deriv <- polynomial.derivatives(chebs)
    cheb.func <- sapply(1:Ncheb, function(i) as.function(chebs[[i]]))
    dcheb.func <- sapply(1:Ncheb, function(i) as.function(cheb.deriv[[i]]))
    chebT <- chebX <- chebY <- chebZ <- array(NA,dim=c(Nt,Ncheb))
    t2 <- c()
    time1 <- proc.time()
    cheb <- array(NA,dim=c(Nt,Ncheb,Ndim))
    for(i in 1:Nt){
        jd.tdb <- sum(tdb[i,])
        ind <- which((DEfile[,1]-jd.tdb)>0)[1]-1
        dt <- (DEfile[ind,2]-DEfile[ind,1])/Nsub
        ts <- seq(DEfile[ind,1],DEfile[ind,2],by=dt)
        ind.sub <- which(ts[-1]>jd.tdb & ts[-length(ts)]<=jd.tdb)
        for(k in 1:Ndim){
            cheb[i,,k] <- as.numeric(DEfile[ind,ind.start-1+(ind.sub-1)*(Ncheb*Ndim)+(k-1)*Ncheb+(1:Ncheb)])
        }
###scaling the JD[TDB] so that it is in the range of [-1,1]
        t0 <- ts[ind.sub]
        t1 <- ts[ind.sub+1]
        t2 <- c(t2,((tdb[i,1]-(t0+t1)/2)+tdb[i,2])*2/(t1-t0))#scaled time
    }
    t2T <- 2/(t1-t0)/(3600*24)#per day to per second
    var <- dvar <- c()
    if(Nt==1){
        for(k in 1:Ndim){
            var <- cbind(var,sum(sapply(1:Ncheb, function(i) cheb.func[[i]](t2)*cheb[,i,k])))
            dvar <- cbind(dvar,sum(sapply(1:Ncheb, function(i) dcheb.func[[i]](t2)*cheb[,i,k]))*t2T)
        }
    }else{
        for(k in 1:Ndim){
            var <- cbind(var,rowSums(sapply(1:Ncheb, function(i) cheb.func[[i]](t2)*cheb[,i,k])))
            dvar <- cbind(dvar,rowSums(sapply(1:Ncheb, function(i) dcheb.func[[i]](t2)*cheb[,i,k]))*t2T)
        }
    }
    out <- cbind(var,dvar)#xyz[km], vxyz[km/s]
###barycenter of the Earth-moon center
    return(out)
}
gen_GeoEph <- function(tdb,DE=430){
####################################
## Ephemeris of Earth's geocenter in the SSB reference frame
##
## Input:
##   tdb  TDB time
##   DE   JPL Ephemeris version
##
## Output:
##   out  Including 2D TDB time, position in km, velocity in km/s, TT-TDB, d(TT-TDB)/dTDB
##   emrat Earth-Moon mass ratio
##   GM Lunar ephemeris in the geocentric refrence frame
####################################
    SEM <- gen_CalEph(tdb,body='Earth-Moon barycenter',DE=DE)
    GM <- gen_CalEph(tdb,body='Moon (geocentric)',DE=DE)
###
    EMRAT <- gen_GetEMRAT(DE)
###pos+vel of geocenter in the BCRS
    out <- SEM[,1:6,drop=FALSE]-GM[,1:6,drop=FALSE]/(1+EMRAT)

###TT-TDB (geocenter)
#    dt <- gen_CalEph(tdb,body='TT-TDB (at geocenter)',DE=DE)/DAYSEC#day
#    out <- cbind(tdb,out,dt)
    out <- cbind(tdb,out)
#    colnames(out) <- c('JD1[TDB]','JD2[TDB]','x.km','y.km','z.km','vx.kms','vy.kms','vz.kms','TT-TDB','d(TT-TDB)/dTDB')
    colnames(out) <- c('JD1[TDB]','JD2[TDB]','x.km','y.km','z.km','vx.kms','vy.kms','vz.kms')
    return(list(out=out,emrat=EMRAT,GM=GM))
}

gen_R1 <- function(theta){
####################################
## Rotation matrix for rotation around x axis
##
## Input:
##   theta -  Counter-clockwise rotation angle
##
## Output:
##   3D rotation matrix
####################################
    as.matrix(rbind(c(1,0,0),c(0,cos(theta),sin(theta)),c(0,-sin(theta),cos(theta))))
}

gen_R2 <- function(theta){
####################################
## Rotation matrix for rotation around y axis
##
## Input:
##   theta -  Counter-clockwise rotation angle
##
## Output:
##   3D rotation matrix
####################################
    as.matrix(rbind(c(cos(theta),0,-sin(theta)),c(0,1,0),c(sin(theta),0,cos(theta))))
}

gen_R3 <- function(theta){
####################################
## Rotation matrix for rotation around z axis
##
## Input:
##   theta -  Counter-clockwise rotation angle
##
## Output:
##   3D rotation matrix
####################################
    as.matrix(rbind(c(cos(theta),sin(theta),0),c(-sin(theta),cos(theta),0),c(0,0,1)))
}
gen_CalMcio <- function(s,x,y){
####################################
## Calculate the Mcio matrix
## CIO is the "celestial intermediate origin"
## CIP is the "celestial intermediate pole"
## ref: https://www.aanda.org/articles/aa/abs/2006/17/aa4550-05/aa4550-05.html
##
## Input:
##   s - CIO locator (a small angle)
##   x - x coordinate of CIP
##   y - y coordinate of CIP
##
## Output:
##   Mcio - geocentric celestial reference system (GCRS) to the celestial intermediate reference system (CIRS) matrix
####################################
    z <- sqrt(1-x^2-y^2)
    a <- 1/(1+z)
    rbind(
        c(cos(s)+a*x*(y*sin(s)-x*cos(s)), -sin(s)+a*y*(y*sin(s)-x*cos(s)), -(x*cos(s)-y*sin(s))),
        c(sin(s)-a*x*(y*cos(s)+x*sin(s)), cos(s)-a*y*(y*cos(s)+x*sin(s)), -(y*cos(s)+x*sin(s))),
        c(x,y,z)
    )
}

gen_refraction <- function(zenith,ZenIn,Par){
####################################
## Calculate atmospheric refraction
##
## Input:
##   zenith - Zenith vector
##   ZenIn - Real zenith angle for incident light ray
##   Par - Input parameters
##
## Output:
##   Ratm - Atmospheric refraction (rad)
##   ZenObs - Observed zenith direction
####################################
    Ntry <- 100
    tol <- 1/206265
    tol.rel <- 1e-5
    ZenObs <- ZenIn
    if(Par$RefType=='S86+'){
        s <- Par$pmb*0.1/101*283/Par$tdk
        elevation <- pi/2-ZenIn
        ref <- 1.02*cot(elevation+10.3/(elevation+5.11))/60/180*pi
        ZenObs <- ZenIn+ref
    }else{
        ref1 <- 0
        for(j in 1:Ntry){
            if(Par$RefType=='refro'){
                ref <- astro_refro(ZenObs,Par$hm,Par$tdk,Par$pmb,Par$rh,Par$wl,Par$phi,Par$tlr)
            }else if(Par$RefType=='refcoq' | Par$RefType=='refco'){
                zreal.deg <- ZenObs/180*pi
                ind0 <- which(zreal.deg>=5 & zreal.deg<=85)
                ind1 <- which(zreal.deg<5 | zreal.deg>85)
                ref <- rep(NA,Par$Nepoch)
                if(length(ind0)>0){
                    if(Par$RefType=='refcoq'){
                        refab <- astro_refcoq(Par$tdk,Par$pmb,Par$rh,Par$wl)
                    }else{
                        refab <- astro_refco(Par$hm,Par$tdk,Par$pmb,Par$rh,Par$wl,Par$phi)
                    }
                    A <- refab$refa
                    B <- refab$refb
                    ref[ind0] <- A*tan(ZenObs[ind0])+B*tan(ZenObs[ind0])^3
                }
                if(length(ind1)>0){
                    ref[ind1] <- astro_refro(ZenObs[ind1],Par$hm,Par$tdk,Par$pmb,Par$rh,Par$wl,Par$phi,Par$tlr)
                }
            }else if(Par$RefType=='B82+'){
                s <- Par$pmb*0.1/101*283/Par$tdk
                elevation <- pi/2-ZenObs
                ref <- s*cot(elevation+7.31/(elevation+4.4))/60/180*pi
            }
            dref <- ref-ref1
            ddref <- abs(dref/ref)
            if(all(ddref<tol.rel)){
                break()
            }else{
                ZenObs <- ZenIn-ref
                ref1 <- ref
            }
        }
    }
    list(Ratm=ref,ZenObs=ZenObs)
}

gen_KeplerSolver <- function(m,e){
####################################
## Solve the Kepler's Equation
##
## Input:
##   m - Mean anomaly
##   e - Eccentricity
##
## Output:
##   E1 - Eccentric anomaly
####################################
    tol = 1e-8
    E0 <- m
    Ntt <- 1e2
    for(k in 1:Ntt){
        E1 = E0-(E0-e*sin(E0)-m)/(sqrt((1-e*cos(E0))^2-(E0-e*sin(E0)-m)*(e*sin(E0))))
        if(all(abs(E1-E0)<tol)) break()
        if(k==Ntt){
            cat('Warning: Keplerian solver does not converge!\n')
            cat('Warning: length(which(abs(E1-E0)>tol))=',length(which(abs(E1-E0)>tol)),'\n')
        }
        E0 <- E1
    }
    return(E1)
}

gen_KeplerClassic <- function(tB,Par){
####################################
## Solve the Kepler's Equation
##
## Input:
##   tB - Coordinate time at the TSB
##   Par - Input parameters
##
## Output:
##   state - the coordinates of the barycentric position and velocity of the target in the sky-plane frame defined by the [p q u] triad; (au; au/yr)
##   U - eccentricity anomaly
##   xo - x coordiante of the target in the orbit-plane frame (au)
##   xy - y coordiante of the target in the orbit-plane frame (au)
##   A - One of the scaled Thiele Innes constants
##   B - One of the scaled Thiele Innes constants
##   F - One of the scaled Thiele Innes constants
##   G - One of the scaled Thiele Innes constants
##   C - One of the scaled Thiele Innes constants
##   H - One of the scaled Thiele Innes constants
####################################
    m1 <- Par$mT
    m <- Par$mTC
    m2 <- Par$mC
    mu <-m1*m2/m
    Py <- Par$P#year
    Pd <- Par$Pd#day
    e <- Par$e
    I <- Par$I
    Omega <- Par$Omega
    omega <- Par$omegaT
    dt <- (tB[,1]-Par$Tp)+tB[,2]
    M0 <- 2*pi*(dt%%Pd)/Pd#mean anomaly
    n <- 2*pi/Py#1/yr
    arr <- (m*Py^2)^{1/3}#au
    ar <- arr*m2/m
    E <- gen_KeplerSolver(M0,e)%%(2*pi)
####Keplerian motion in the orbital plane
    x <- ar*(cos(E)-e)#au
    y <- ar*(sqrt(1-e^2)*sin(E))
    vx <- -ar*n*sin(E)/(1-e*cos(E))#au/yr
    vy <- ar*n*sqrt(1-e^2)*cos(E)/(1-e*cos(E))
###Thiele Innes constants; ref. Wright et al. 2009
###Definition:(ex:North; ey:East; ez:target to observer)
#    ae <- gen_Xy2phi(cos(E)-e)/(1-e*cos(E),sqrt(1-e^2)*sin(E)/(1-e*cos(E)))
    A <- cos(Omega)*cos(omega)-sin(Omega)*sin(omega)*cos(I)
    B <- sin(Omega)*cos(omega)+cos(Omega)*sin(omega)*cos(I)
    F <- -cos(Omega)*sin(omega)-sin(Omega)*cos(omega)*cos(I)
    G <- -sin(Omega)*sin(omega)+cos(Omega)*cos(omega)*cos(I)
    C <- sin(omega)*sin(I)#this is negative in Catanzarite 2010
    H <- cos(omega)*sin(I)#this is negative in Catanzarite 2010
####sky plane coordinates [pb, qb, ub]
    Y <- A*x+F*y
    X <- B*x+G*y
    Z <- C*x+H*y#Z is from observer to target; opposite in Catanzarite 2010
    VY <- A*vx+F*vy
    VX <- B*vx+G*vy
    VZ <- C*vx+H*vy
    list(state=cbind(X,Y,Z,VX,VY,VZ),U=E,xo=x,yo=y,A=A,B=B,F=F,G=G,C=C,H=H)
}

gen_CalBT <- function(tB,Par){
####################################
## Calculate the barycentric motion of the target
## i.e. the position and velocity vectors from TSB to the target
##
## Input:
##   tB - Coordinate time at the TSB
##   Par - Input parameters
##
## Output:
##   BT - Positiona and velocity vectors from the SSB to the target
##   out - Output from gen_KeplerClassic() or gen_DDmodel() or gen_DDGRmodel()
####################################
    if(Par$BinaryModel=='kepler'){
##classical Keplerian orbit
        out <- gen_KeplerClassic(tB,Par)
        BT <- out$state
    }else if(Par$BinaryModel=='DD'){
        out <- gen_DDmodel(tB=tB,Par)
        BT <- cbind(out$rvec,out$vvec)
    }else if(Par$BinaryModel=='DDGR'){
##Einstein theory-based Keplerian orbit
        out <- gen_DDGRmodel(tB=tB,Par)#bbat!=te; change
        BT <- cbind(out$rvec,out$vvec)
    }
##transform BT from sky-plane coordinates, [p,q,u] to ICRS coordinates
    x <- BT[,1]%o%Par$p
    y <- BT[,2]%o%Par$q
    z <- BT[,3]%o%Par$u
    vx <- BT[,4]%o%Par$p
    vy <- BT[,5]%o%Par$q
    vz <- BT[,6]%o%Par$u
    rBT <- x+y+z
    vBT <- vx+vy+vz
#    if(Par$Unit=='TDB') rBT <- rBT/IFTE.K
    c(list(BT=cbind(rBT,vBT)),out)
}

gen_One2two <- function(var){
####################################
## Split a string into multiple parts
## It is used to read the output of TEMPO2 general2 extension
##
## Input:
##   var - variable or string
##
## Output:
##   Splited substrings
####################################
    var <- gsub('\\.',' 0\\.',var)
    var <- unlist(strsplit(var,' '))
    as.numeric(var[var!=''])
}

gen_GetPhase <- function(e,omega,type='primary'){
####################################
## Get the phase of the target in its barycentric motion
##
## Input:
##   e - eccentricity
##   omega - Argument of periastron
##   type -
##     primary - Primary transit
##     secondary - Secondary transit
##     periastron - Periastron
##     ascendingnode - Ascending node
##     descendingnode - Descending node
##     l4 - L4 point
##     l5 - L5 point
##
## Output:
##   Phase in fraction
####################################
    if(type=='primary') theta <- pi/2-omega
    if(type=='secondary') theta <- 3*pi/2-omega
    if(type=='periastron') theta <- 0
    if(type=='ascendingnode') theta <- -omega
    if(type=='descendingnode') theta <- pi-omega
    if(type=='l4') theta <- 5*pi/6-omega
    if(type=='l5') theta <- pi/6-omega
     E <- 2*atan(sqrt((1-e)/(1+e))*tan(theta/2))
     M <- E-e*sin(E)
     M/(2*pi)
}

gen_CalTc <- function(Tp,P,e,omega){
####################################
## Calculate the primary transit epoch
##
## Input:
##   Tp - Reference epoch at the periastron
##   P - orbital period
##   e - eccentricity
##   omega - Argument of periastron
##
## Output:
##   Tc - Reference epoch at the primary transit
####################################
    theta <- 0.5*pi-omega
    E <- 2*atan(sqrt((1-e)/(1+e))*tan(theta/2))
    M <- E-e*sin(E)
    Tc <- (M%%(2*pi))/(2*pi)*P+Tp
    Tc
}

gen_Tc2tp <- function(e,omega,P,Tc,type='primary'){
####################################
## Convert the primarty transit epoch to the perastron epoch
##
## Input:
##   e - eccentricity
##   omega - Argument of periastron
##   P - orbital period
##   Tc - Reference epoch at the primary transit
##   type - type of phase (see comments for gen_GetPhase() for details)
##
## Output:
##   Tp - Reference epoch at the periastron
####################################
    Tc-gen_GetPhase(e,omega)*P
}

gen_DDmodel <- function(tB,Par){
####################################
## Binary motion modeled by the DD binary model; ref: DD86
## In TEMPO2, natural units are used, i.e. c=g=1 which need to change into au, au/yr in PEXO.
##
## Input:
##   tB - Coordinate time of light arrival time at the TSB
##   Par -
##
## Output:
##   torb - Roemer delay due to the target barycentric motion
##   pb - orbital period
##   a1 - semi-major axis of the target barycentric motion multiplied by sin(I)
##   ecc - eccentricity
##   edot - eccentricity time derivative
##   omega - Argument of periastron
##   omdot - omega time derivative
##   pbdot - orbital period time derivative
##   sini - sin(I) where I is inclination
##   gamma - Timing model parameter
##   m2 - mass of the companion
##   a1dot - time derivative of the semi-major axis
##   rvec - the coordinates of the barycentric position of the target in the sky-plane frame defined by the [p q u] triad; (au)
##   vvec - the coordinates of the barycentric position of the target in the sky-plane frame defined by the [p q u] triad; (au/yr)
##   U - eccentrisic anomaly
##   Dre - Roemer and Einstein delay
##   Ds - Shapiro delay
##   Da - aberration time delay (eq. 9b in DD86)
##   x - Time-varing a1 equivalent: a1/c + xdot*Dt
##   xo - x coordiante of the target in the orbit-plane frame (au)
##   xy - y coordiante of the target in the orbit-plane frame (au)
####################################
    Nt <-nrow(tB)
    dr = 0.0##WHAT SHOULD THESE BE SET TO?
    dth = 0.0
    si <- Par$DDGR$sini
    am2 <- Par$mC
    am<- Par$mTC

    if(Par$BinaryUnit=='natural'){
        pb = Par$P*DJY*DAYSEC#s
    }else{
        pb = Par$P#year
    }
    an = 2.0*pi/pb#mean motion in unit of rad/second or rad/yr
    if(Par$BinaryUnit=='natural'){
        n <- rad2deg*365.25*86400.0*an#deg/yr
    }else{
        n <-rad2deg*an#deg/yr
    }
    omdot <-Par$DDGR$omdot#
    k  <-  omdot/n#[omdot]=deg/yr
    Omega <- Par$Omega

    if(Par$BinaryUnit=='natural'){
        m <- am*SUNMASS#in unit of second because g=c=1
        m2 = am2*SUNMASS
    }else{
        m <- am
        m2 <- am2
    }
    m1 <- m-m2
    ct = tB

    if(Par$BinaryUnit=='natural'){
        tt0  <-  ((ct[,1]-Par$Tp)+ct[,2])*DAYSEC#day to second
    }else{
        tt0  <-  ((ct[,1]-Par$Tp)+ct[,2])/DJY#day to year
    }

    gamma  <-  Par$DD$gamma#yr or s; magnitude of Einstein delay
    a0     <-  0.0## WHAT SHOULD THIS BE SET TO?
    b0     <-  0.0## WHAT SHOULD THIS BE SET TO?

    omz =Par$DD$omz

    xdot  <- Par$DD$xdot
    pbdot <-  Par$DD$pbdot
    edot <- Par$DD$edot
    xpbdot <- Par$DD$xpbdot

    ##a1==x0 in Eqn. 71 of E06; define x==a*sin(i)/c
    if(Par$BinaryUnit=='natural'){
        x = Par$DDGR$x0+xdot*tt0#[x]=s; [xdot]=s/s;
        ar <-x/si#s
    }else{
        x = Par$DDGR$x0+xdot*tt0#[x]=yr, [xdot]=yr/yr
        ar <- x*YC/si#au
    }
    arr <-ar/am2*am
    ecc = Par$e+edot*tt0
    er = ecc*(1.0+dr)
    eth = ecc*(1.0+dth)
    if(any(ecc < 0.0 |  ecc > 1.0)){
        cat("DDmodel: problem with eccentricity = ",ecc,"\n")
    }

    orbits = tt0/pb - 0.5*(pbdot+xpbdot)*(tt0/pb)*(tt0/pb)
    pb <- pb+(pbdot+xpbdot)*tt0
    norbits  <-  floor(orbits)

    phase=2.0*pi*(orbits-norbits)
    ## Compute eccentric anomaly u by iterating Kepler's equation.
    U <- gen_ComputeU(phase,ecc)

    ##  DD equations 17b, 17c, 29, and 46 through 52
    su=sin(U)
    cu=cos(U)
    onemecu=1.0-ecc*cu
    cae=(cu-ecc)/onemecu
    sae=sqrt(1.0-ecc^2)*su/onemecu
    ae=atan2(sae,cae)

####Eqn. 58 of E06 is wrong; Eqn. 4.11b of DD85 and Eqn. 17a DD86 gives ae=2*atan2((1+e)^2/(1-e)^2*tan(u/2))
    ae[ae<0]=ae[ae<0]+2.0*pi

    ae=2.0*pi*orbits + ae - phase
    omega=omz + k*ae
    sw=sin(omega)
    cw=cos(omega)
    alpha=x*sw#s or yr
    beta=x*sqrt(1-eth^2)*cw#s or yr
    bg=beta+gamma#gamma is from input pars
    dre <- alpha*(cu-er) + bg*su
    drep=-alpha*su + bg*cu
    drepp=-alpha*cu - bg*su
    anhat=an/onemecu

    ##DD equations 26, 27, 57:
    sqr1me2=sqrt(1-ecc^2)
    cume=cu-ecc
    brace=onemecu-si*(sw*cume+sqr1me2*cw*su)
    ##printf("GEORGE: si = %g, brace = %g\n",(double)si,(double)brace)
    dlogbr=log(brace)
    if(Par$BinaryUnit=='natural'){
        Ds <- -2*m2*dlogbr
    }else{
        r <-4*pi^2*m2/Cauyr^3#yr
        Ds <- -2*r*dlogbr#yr
    }
    #aberration delay
    if(Par$BinaryUnit=='natural'){
        Da=(a0*(sin(omega+ae) + ecc*sw) + b0*(cos(omega+ae) + ecc*cw))#s
    }else{
        Da=a0/Cauyr*(sin(omega+ae) + ecc*sw) + b0/Cauyr*(cos(omega+ae) + ecc*cw)#year
    }

    ##Now compute d2bar, the orbital time correction in DD equation 42. ##
    Dre <- dre*(1-anhat*drep+(anhat^2)*(drep^2 + 0.5*dre*drepp - 0.5*ecc*su*dre*drep/onemecu))
###change unit
    Dre <- Dre*DJY
    Ds <- Ds*DJY
    Da <- Da*DJY
###unit: day
    d2bar  <-  Dre+ Ds + Da
    if(Par$BinaryUnit=='natural'){
        torb <- -d2bar/DAYSEC
    }else{
        torb <- -d2bar*DJY
    }

    ##Now we need the partial derivatives. Use DD equations 62a - 62k. ##
    csigma=x*(-sw*su+sqr1me2*cw*cu)/onemecu#[csigma]=s
    ce=su*csigma-x*sw-ecc*x*cw*su/sqr1me2#[ce]=s
    cx=sw*cume+sqr1me2*cw*su#[cx]=1
    comega=x*(cw*cume-sqr1me2*sw*su)
    cgamma=su
    ##cdth=-ecc*ecc*x*cw*su/sqr1me2
    cm2=-2*dlogbr
    csi=2*m2*(sw*cume+sqr1me2*cw*su)/brace#
    a1dot <- Par$DDGR$xdot#a1dot equivalent to xdot
    if(FALSE){
        pb <- -csigma*an*DAYSEC*tt0/(pb*DAYSEC)
        omdot <- ae*comega/(an*360.0/(2.0*pi)*365.25*DAYSEC)
        t0 <- -csigma*an*DAYSEC
        pbdot <- 0.5*tt0*(-csigma*an*DAYSEC*tt0/(pb*DAYSEC))
        edot <- ce*tt0
        a1dot <- cx*tt0
    }

###derive vectors
    ci  <- cos(asin(si))
    sini <- rep(si,Nt)
    b <- ar*(1-er*cu)#s or au
#    p <-ar*(1-er^2)
    p <-ar*(1-er^2)
    theta <- omega+ae
    mu <-m1*m2/m
###DD85 B.9

    q <- mu/m
    if(Par$BinaryUnit=='natural'){
        s <-m^2/p^3
    }else{
        s <-(2*pi)^4*m^2/YC^2/p^3
    }
    Rw <-s*(1+er*cos(ae))^2*(3-q+3*er^2-3.5*q*er^2+(2-4*q)*er*cos(ae)+(-4+0.5*q)*er^2*(cos(ae))^2)
    Tw <-s*(1+er*cos(ae))^3*(4-2*q)*er*sin(ae)
    adot <- 2/(an*sqrt(1-er^2))*(Rw*er*sin(ae)+(1+er*cos(ae))*Tw)
    edot <- sqrt(1-er^2)/(an*ar)*(Rw*sin(ae)+(er+2*cos(ae)+er*(cos(ae))^2)/(1+er*cos(ae))*Tw)
    udot <-(an+edot*su)/(1-er*cu)
    bdot <-v <- adot*(1-er*cu)-ar*edot*cu+ar*er*su*udot
    theta.dot <- an*(1+er*cos(ae))^2/(1-er^2)^(3/2)
    if(FALSE){
        udot <- an/(1-ecc*cu)
        bdot <- a1dot*(1-ecc*cu)-ar*(cu*edot-su*udot)
        etilde <-(eth+er)/2
        theta.dot <- omdot*sqrt(1-eth^2)/(1-etilde*cu)^2#derived from Eqn.4.2 of DD85
    }
    ##rotation matrix
###default Definition:(ex:North; ey:East; ez:target to observer);
###Tempo2 definition: (ex: East, alpha increase; ey: North, delta increase; ez: obserer to target)
    ##Thiele Innes elements with omega=0
    xo <- b*cos(theta)
    yo <- b*sin(theta)
    vxo <- bdot*cos(theta)-b*sin(theta)*theta.dot
    vyo <- bdot*sin(theta)+b*cos(theta)*theta.dot

if(FALSE){
    x <- ar*(cos(E)-e)#au
    y <- ar*(sqrt(1-e^2)*sin(E))
    vx <- -ar*n*sin(E)/(1-ecc*cos(E))#au/yr
    vy <- ar*n*sqrt(1-e^2)*cos(E)/(1-ecc*cos(E))

    A <- cos(Omega)*cos(omega)-sin(Omega)*sin(omega)*cos(I)
    B <- sin(Omega)*cos(omega)+cos(Omega)*sin(omega)*cos(I)
    F <- -cos(Omega)*sin(omega)-sin(Omega)*cos(omega)*cos(I)
    G <- -sin(Omega)*sin(omega)+cos(Omega)*cos(omega)*cos(I)
    C <- sin(omega)*sin(I)#this is negative in Catanzarite 2010
    H <- cos(omega)*sin(I)#this is negative in Catanzarite 2010
####sky plane coordinates [pb, qb, ub]
}
    Y <- A*xo+F*yo
    X <- B*xo+G*yo
    Z <- C*x+H*y#Z is from observer to target; opposite in Catanzarite 2010
    VY <- A*vx+F*vy
    VX <- B*vx+G*vy
    VZ <- C*vx+H*vy
    Y <- A*xo+F*yo#ex is to the increasing direction of alpha
    X <- B*xo+G*yo#ey is to the increasing direction of delta
    Z <- C*xo+H*yo#ez is from observer to target; opposite in Catanzarite 2010
    VY <- A*vxo+F*vyo
    VX <- B*vxo+G*vyo
    VZ <- C*vxo+H*vyo
    rvec <- cbind(X,Y,Z)
    vvec <- cbind(VX,VY,VZ)
#        R1 <- rbind(c(-cos(Omega),sin(Omega),0),c(cos(Omega),sin(Omega),0),c(0,0,1))
#        R2 <- rbind(c(1,0,0),c(0,-cosi[j],-sini[j]),c(0,sini[j],-cosi[j]))
#        R <- R1%*%R2
#        r <- rbind(b[j]*cos(theta[j]),b[j]*sin(theta[j]),0)
        ##orbital plane
        ##from oribtal plane to projected plane
#        rdot <- rbind(bdot[j]*cos(theta[j])-b[j]*sin(theta[j])*theta.dot[j],bdot[j]*sin(theta[j])+b[j]*cos(theta[j])*theta.dot[j],0)
#        vvec <- rbind(vvec,t(R%*%rdot))
#    }
    list(torb=torb,pb=pb,a1=cx,ecc=ce,edot=edot,omega=omega,omdot=omdot,pbdot=pbdot,sini=csi,gamma=cgamma,m2=cm2*Tsun,a1dot=a1dot,rvec=rvec,vvec=vvec,U=U,Dre=Dre,Ds=Ds,Da=Da,x=x,xo=xo,yo=yo)
}

gen_ComputeU <- function(phase,ecc){
####################################
## Compute eccentric anomaly U by iterating Kepler's equation if
## eccentricity is set. The equation is solved using a Newton-Raphson
## technique and the S9 ##starting value in Odell & Gooding 1986
## CeMec 38 307 ##
##
## Input:
##   phase - Orbital phase
##   ecc - Eccentricity
##
## Output:
##   U - Eccentric anomaly
####################################
    U <- phase+ecc*sin(phase)*(1.0+ecc*cos(phase))
    for(k in 1:1000){
        fac  <-  1.0/(1.0-ecc*cos(U))#NOTE COULD BE WRONG IN DDmodel - SEE USE OF FAC!!!!
        dU <- (phase-(U-ecc*sin(U)))*fac
        U <- U+dU
        if(all(abs(dU)<1.0e-14)) break()
    }
    U
}

gen_DDGRmodel <- function(tB,Par){
####################################
## Binary motion modeled by the DDGR binary model (or DD model assume GR); ref: DD86
## In TEMPO2, natural units are used, i.e. c=g=1 which need to change into au, au/yr in PEXO.
## ref toDDGRmodel.C of tempo2
## Based on bnryddgr.f
## Damour & Deruelle model assuming the general theory of relativity
## Computes the pulsar orbit time, torb, at the time of observation, t =ct(n)-pepoch
## Pulsar proper time is TP = T + TORB
##
## Input:
##   tB - Coordinate time of light arrival time at the TSB
##   Par - Input parameters
##
## Output:
##   torb - Roemer delay due to the target barycentric motion
##   a1 - semi-major axis of the target barycentric motion multiplied by sin(I)
##   xpbdot - time derivative of xpb
##   sini - sin(I)
##   m2 - mass of the secondary
##   ar - semi-major axis of binary motion of the target with respect to the companion or vice verse
##   rvec - the coordinates of the barycentric position of the target in the sky-plane frame defined by the [p q u] triad; (au)
##   U - eccentrisic anomaly
##   Dre - Roemer and Einstein delay
##   Ds - Shapiro delay
##   Da - aberration time delay (eq. 9b in DD86)
##   pb - orbital period
##   dpb - variation of orbital period since the reference epoch
##   edot - eccentricity time derivative
##   omega - Argument of periastron
##   alpha -  x*sw (yr)
##   abdot - x*(cos(omega[1])-sin(omega[1]))*k*an#assume e=0
##   beta - x*sqrt(1-eth^2)*cw (yr)
##   gamma - Timing model parameter
##   dre - alpha*(cu-er) + bg*su
##   drep - -alpha*su + bg*cu#year
##   drepp - -alpha*cu - bg*su#
##   anhat - an/onemecu
##   Dre1 - first order Roemer+Einstein delay in the target system
##   Dre2 - first order Roemer+Einstein delay in the target system
##   x - Time-varing a1 equivalent: a1/c + xdot*Dt
##   dr - an^(2/3)*(3*m1^2+6*m1*m2+3*m2^2)/m^(4/3)*const^(2/3)
##   dth - an^(2/3)*(3.5*m1^2+6*m1*m2+2*m2^2)/m^(4/3)*const^(2/3)
##   an - angular frequency; 2.0*pi/pb
##   ecc - eccentricity
##   Tc - Primary transit epoch
##   Vc - Velocity of the target with respect to the companion
##   bdot - Time derivative of impact parameter; bdot <- v <- adot*(1-ecc*cu)-ar*edot*cu+ar*ecc*su*udot
##   theta.dot - variation of theta which is the angular distance from the target to the periastron
##   adot - time derivative of semi-major axis
##   xo - x coordiante of the target in the orbit-plane frame (au)
##   yo - y coordiante of the target in the orbit-plane frame (au)
##   out - output from gen_mass2dd()
####################################
    Nt <- nrow(tB)
    ct <-  tB#binary barycentric time
    if(Par$BinaryUnit=='natural'){
        tt0  <-  ((ct[,1]-Par$Tp)+ct[,2])*DAYSEC#day to second
    }else{
        tt0  <-  ((ct[,1]-Par$Tp)+ct[,2])/DJY#day to year
    }
    f0 <- 1#pulsar frequency; set to 1 for exoplanet cases

    xomdot <- afac <- 0#extra variation
    si <- Par$DD$sini
    am2 <-Par$mC
    am <-Par$mTC
    am1 <-Par$mT

    if(Par$BinaryUnit=='natural'){
        m <- am*SUNMASS#in unit of second because g=c=1
        m2 = am2*SUNMASS
    }else{
        m <- am
        m2 <- am2
    }
    m1 = m-m2

    if(Par$BinaryUnit=='natural'){
        pb <-  Par$DDGR$pb*DAYSEC
    }else{
        pb <-  Par$DDGR$pb
    }
    an  <-  2.0*pi/pb

####import parameters
    xdot <- Par$DD$a1dot
    omz  <-  Par$DDGR$omz
    xdot  <- Par$DDGR$xdot
    pbdot <-  Par$DDGR$pbdot
    edot <- Par$DDGR$edot
    xpbdot <- Par$DDGR$xpbdot
    Omega <- Par$Omega

    x = Par$DDGR$x0+xdot*tt0#in unit of time
    ar <- x*YC/si
    ecc <- Par$DD$ecc+edot*tt0

    out <- gen_mass2dd(m,m2,x,ecc,an)
    arr <- out$arr
    ar <- out$ar
    si <- out$si
    gamma <- out$gamma
    pbdot <- out$pbdot
    omdot <-out$omdot
    k <- out$xk

###DD Eqn. 36, 37
    if(Par$BinaryUnit=='natural'){
        dr  <-  (3.0*m1^2 + 6.0*m1*m2 + 2.0*m2^2)/(arr*m)
        dth  <-  (3.5*m1^2 + 6*m1*m2 + 2*m2*m2)/(arr*m)
    }else{
        dr  <-  4*pi^2/Cauyr^2*(3.0*m1^2 + 6.0*m1*m2 + 2.0*m2^2)/(arr*m)
        dth  <-  4*pi^2/Cauyr^2*(3.5*m1^2 + 6*m1*m2 + 2*m2*m2)/(arr*m)
    }
    er  <-  ecc*(1.0+dr)
    eth  <-  ecc*(1.0+dth)

    orbits  <- tt0/pb - 0.5*(pbdot+xpbdot)*(tt0/pb)^2
    dpb <- (pbdot+xpbdot)*tt0
    pb <- pb+dpb
    dpb <- SJY*dpb
    norbits  <- floor(orbits)
    phase <- 2.0*pi*(orbits-norbits)#mean anomaly
    #Compute eccentric anomaly u by iterating Kepler's equation.
    U <- gen_ComputeU(phase,ecc)

    ##DD equations 17a, 29
    ae  <-  2.0*atan(sqrt((1+ecc)/(1-ecc))*tan(0.5*U))
    ae[ae<0] <- ae[ae<0]+2.0*pi
    ae  <-  2.0*pi*orbits + ae-phase
    omega <- omz + (k+xomdot/(an*rad2deg*365.25*86400.0))*ae
    ##DD equations 46 through 52
    su <- sin(U)
    cu <- cos(U)
    sw <- sin(omega)
    cw <- cos(omega)
    alpha <- x*sw#yr
    beta <- x*sqrt(1-eth^2)*cw#yr
    abdot <-x*(cos(omega[1])-sin(omega[1]))*k*an#assume e=0
    bg <- beta+gamma
    dre <- alpha*(cu-er) + bg*su#year
#    dre.u0 <- alpha*(1-er)
#    dre.upi2 <- alpha*(-er) + bg
#    dre.upi <- alpha*(-1-er)
    drep <- -alpha*su + bg*cu#year
    drepp <- -alpha*cu - bg*su#year
    onemecu <- 1.0-ecc*cu
    anhat <- an/onemecu

###DD equations 26, 27, 57
###DD equations 26,27,57
    cume=cu-ecc
    sqr1me2=sqrt(1-ecc^2)
    brace=onemecu-si*(sw*cume+sqr1me2*cw*su)
    if(any(brace<=0))
    {
        cat("ERROR: In DDGR model, brace < 0\n")
    }
    dlogbr=log(brace)
    if(Par$BinaryUnit=='natural'){
        Ds <- -2*m2*dlogbr#s
    }else{
        r <-4*pi^2*m2/Cauyr^3
        Ds <- -2*r*dlogbr#yr
    }

    ##These will be different if spin axis not aligned -- IS THIS AN ASSUMPTION OF THE MODEL?
    a0aligned = an*ar/(2.0*pi*f0*si*sqr1me2)
    a0 = afac*a0aligned
    b0 = 0.0
    if(Par$BinaryUnit=='natural'){
        Da = a0*(sin(omega+ae)+ecc*sw) + b0*(cos(omega+ae) + ecc*cw)
    }else{
        Da <-a0/Cauyr*(sin(omega+ae)+ecc*sw) + b0/Cauyr*(cos(omega+ae) + ecc*cw)
    }

    ##Now compute d2bar, the orbital time correction in DD equation 42.
    Dre1 <- -anhat*dre*drep
    Dre2 <- anhat^2*dre*(drep^2 + 0.5*dre*drepp - 0.5*ecc*su*dre*drep/onemecu)
    Dre <- dre+Dre1+Dre2
#    Dre <- dre*(1-anhat*drep+(anhat^2)*(drep^2 + 0.5*dre*drepp - 0.5*ecc*su*dre*drep/onemecu))

###change unit
    Dre1 <- Dre1*DJY*DAYSEC
    Dre2 <- Dre2*DJY*DAYSEC
    Dre <- Dre*DJY*DAYSEC
    Ds <- Ds*DJY*DAYSEC
    Da <- Da*DJY*DAYSEC

###unit:day
    d2bar <- Dre + Ds + Da
    torb <- -d2bar#second

    ##Now get partial derivatives
    if(Par$BinaryUnit=='natural'){
        an0 <-  sqrt(m/arr^3)#1/s
    }else{
        an0 <- 2*pi*sqrt(m/arr^3)#1/yr
    }

    csigma=x*(-sw*su+sqr1me2*cw*cu)/onemecu
    ce=su*csigma-x*sw-ecc*x*cw*su/sqr1me2
    cx=sw*cume+sqr1me2*cw*su
    comega=x*(cw*cume-sqr1me2*sw*su)
    cgamma=su
    cm2=-2*dlogbr
    csi=2*m2*(sw*cume+sqr1me2*cw*su)/brace

    fact1=(m/(2*arr)) ##((m-m2)*m2/m^2 - 9)
    fact2=(3*m/(2*arr^4)) ##(1.0 + fact1)
    fact3=(m/2*arr) ##(m2/m^2-2*(m-m2)*m2/m^3)
    fact4=(1+fact1)*3*m/(2*arr^4*an0)
    fact5=an0*fact1/arr

    denumm=(1+fact1)/(2*arr^3*an0) + an0*(fact1/m+fact3)
    denomm=fact4+fact5
    darrdm=denumm/denomm
    dnum = an0*(m-2*m2)/(2*arr*m)
    denom = an0*fact1/arr + fact2/an0
    darrdm2 = dnum/denom
    dgmdm2 = ((m+2*m2)/arr - (m2*(m+m2)*darrdm2/arr^2))*ecc/(an*m)
    cdth=-ecc*ecc*x*cw*su/sqr1me2
    dthdm2 = -dth*darrdm2/arr - (m+m2)/(arr*m)
    dkdm = k/m - k*darrdm/arr
    dsidm2 = -(m*x/(arr*m2))*(1.0/m2+darrdm2/arr)
    ck = ae*comega
    dkdm2 = -k*darrdm2/arr
    cdr = -ecc*x*sw
    ddrdm2 = -dr*darrdm2/arr - 2*m2/(arr*m)
    dtdm2 = -2*dlogbr
    csini = 2*m2*(sw*cume+sqr1me2*cw*su)/brace
    dsidm=-(m*x/(arr*m2))*(-1.0/m+darrdm/arr)
    dpbdm = pbdot/(m-m2) - pbdot/(3*m)
    if(FALSE){
        cpbdot = -csigma*an*tt0^2/(2*pb)
        ddrdm = -dr/m - dr*darrdm/arr + 6/arr
        dpbdm2 = pbdot/m2 - pbdot/(m-m2)
        cm2 = dtdm2+cgamma*dgmdm2+csini*dsidm2+ck*dkdm2+cdr*ddrdm2+cdth*dthdm2+cpbdot*dpbdm2
        fact6=1.0/(arr*m)
        fact7=-(m+m2)/(arr*m^2)
        fact8=-(m+m2)*darrdm/(arr^2*m)
        dgamdm = (ecc*m2/an)*(fact6+fact7+fact8)

        dthdm=-dth/m - dth*darrdm/arr + (7*m-m2)/(arr*m)

        cm = ck*dkdm+cgamma*dgamdm+cdr*ddrdm+cdth*dthdm+cpbdot*dpbdm+csini*dsidm
        cpbdot <- 0.5*tt0*(-csigma*an*DAYSEC*tt0/(pb*DAYSEC))#?
        cpb=-csigma*an*DAYSEC*tt0/(pb*DAYSEC)#?
        ce=ce
        om=comega
        t0=-csigma*an*DAYSEC
        pbdot=pbdot
        xpbdot=pbdot
        sini=csi
        m2=cm2*Tsun
        mtot=cm*Tsun
        a1dot=cx*tt0
    }

###derive vectors
##inclination
    ci  <- cos(asin(si))
    sini <- rep(si,Nt)
    cosi  <- cos(asin(sini))
    b <- ar*(1-er*cu)#
    theta <- omega+ae
    mu <-m1*m2/m
###DD85 B.1-9
    if(Par$BinaryUnit=='natural'){
        adot <- Par$DDGR$xdot/Par$DDGR$sini
    }else{
        adot <- Par$DDGR$xdot/Par$DDGR$sini*Cauyr
    }
    udot <-(an+edot*su)/(1-ecc*cu)
    bdot <- v <- adot*(1-ecc*cu)-ar*edot*cu+ar*ecc*su*udot
    theta.dot <- an*(1+ecc*cos(ae))^2/(1-ecc^2)^(3/2)

##rotation matrix
    xo <- b*cos(theta)
    yo <- b*sin(theta)

    vxo <- bdot*cos(theta)-b*sin(theta)*theta.dot
    vyo <- bdot*sin(theta)+b*cos(theta)*theta.dot

    A <- cos(Omega)
    B <- sin(Omega)
    F <- -sin(Omega)*ci
    G <- cos(Omega)*ci
    C <- 0#this is negative in Catanzarite 2010
    H <- si#this is negative in Catanzarite 2010
####sky plane coordinates [pb, qb, ub]
    Y <- A*xo+F*yo#ex is to the increasing direction of alpha
    X <- B*xo+G*yo#ey is to the increasing direction of delta
    Z <- C*xo+H*yo#ez is from observer to target; opposite in Catanzarite 2010
    VY <- A*vxo+F*vyo
    VX <- B*vxo+G*vyo
    VZ <- C*vxo+H*vyo
    rvec <- cbind(X,Y,Z)
    vvec <- cbind(VX,VY,VZ)
###derive the primary transit time
    tt <- rowSums(tB)
    Pd <- Par$Pd
    n <- (tt-Par$Tc)%/%Pd
    dTc <- gen_CalTc((omega-omz)/(2*pi)*Pd,pb,ecc,Par$omega)
    Tc <- cbind(n*Pd+Par$Tp,dTc)
    f <- 0.5*pi-omega
    Vc <- 2*pi*sqrt(Par$mC^2/Par$mTC/arr/(1-ecc^2))*sqrt(1+2*ecc*cos(f)+ecc^2)
    c(list(torb=torb,xpbdot=xpbdot,sini=sini,m2=m2,ar=ar,rvec=rvec,vvec=vvec,U=U,Dre=Dre,Ds=Ds,Da=Da,pb=pb,dpb=dpb,omega=omega,alpha=alpha,abdot=abdot,beta=beta,gamma=gamma,dre=dre,drep=drep,drepp=drepp,anhat=anhat,Dre1=Dre1,Dre2=Dre2,x=x,xdot=xdot,dth=dth,sw=sw,cw=cw,ae=ae,er=er,eth=eth,dr=dr,dth=dth,an=an,ecc=ecc,Tc=Tc,Vc=Vc,bdot=bdot,theta.dot=theta.dot,xo=xo,yo=yo),out)
}

gen_mass2dd <- function(m,m2,x,ecc,an,BinaryUnit='auyr'){
####################################
## Ref to DDGRmodel.C of tempo2:
## Given system masses m,m2 and keplerian parameters x,ecc,an, calculate the values
## of arr,ar,si,gamma,pbdot under general relativity
## Convet a calendar time in format of Year, Month, Day, Hour, Minute, Second to Julian Date
##
## Input:
##   m - Total mass of primary and secondary
##   m2 - Mass of the companion
##   x - a1*sin(I)/c
##   ecc - eccentricity
##   an - angular frequency
##   BinaryUnit - whether or not Par$BinaryUnit units are used.
##
## Output:
##   arr - semi-major axis for companion-target binary motion
##   ar - semi-major axis for barycentric motion of the target
##   si - sin(I)
##   gamma - timing model parameter
##   pbdot - orbital period time derivative
##   xk - 3*m/(arr*(1-ecc^2))
##   omdot - omega time derivative
##   dr - an^(2/3)*(3*m1^2+6*m1*m2+3*m2^2)/m^(4/3)*const^(2/3)
##   dth - an^(2/3)*(3.5*m1^2+6*m1*m2+2*m2^2)/m^(4/3)*const^(2/3)
####################################
    m1 <- m-m2
    if(m<0) cat('ERROR: problem with target system mass <0!\n')
    arr <- gen_CalArr(m,m1,m2,an,BinaryUnit=Par$BinaryUnit)
    ar <- arr*m2/m#semi-major axis of m1
    si <- x*YC/ar
    if(any( si > 1.0 )){
        cat("SIN I > 1.0, setting to 1: should probably use DD model!\n")
        si[si>1]  <-  1.0
    }
    if(Par$BinaryUnit=='natural'){
        xk <-3*m/(arr*(1-ecc^2))
        gamma <- ecc*m2*(m1+2*m2)/(an*arr*m)#s
        const <-1
    }else{
        xk <-3*4*pi^2*m/(arr*(1-ecc^2))/Cauyr^2#
        gamma <-4*pi^2/Cauyr^2*ecc*m2*(m1+2*m2)/(an*arr*m)#yr
        const <- Tsun.yr
    }
    pbdot <- -(96*2*pi/5)*an^(5/3)*(1-ecc^2)^(-3.5)*(1+73*ecc^2/24+37/96*ecc^4)*m2*m1*m^(-1/3)*const^(5/3)
    omdot <- 3*(an)^(5/3)*m^(2/3)/(1-ecc^2)*const^(2/3)
    dr <- an^(2/3)*(3*m1^2+6*m1*m2+3*m2^2)/m^(4/3)*const^(2/3)
    dth <- an^(2/3)*(3.5*m1^2+6*m1*m2+2*m2^2)/m^(4/3)*const^(2/3)
    list(arr=arr,ar=ar,si=si,gamma=gamma,pbdot=pbdot,xk=xk,omdot=omdot,dr=dr,dth=dth)
}

gen_CalArr <- function(m,m1,m2,an,BinaryUnit='auyr'){
####################################
## Calculate the semi-major axis of a companion-target binary orbit
##
## Input:
##   m - Total mass of primary and secondary
##   m1 - Mass of the target
##   m2 - Mass of the companion
##   an - Angular frequency of the orbit
##   BinaryUnit - Unit used in the computation of byinary orbit
##
## Output:
##   arr - semi-major axis for companion-target binary motion
####################################
    ARRTOL <-1e-15
    if(Par$BinaryUnit=='natural'){
        arr0 <- (m/an^2)^(1/3)#s or au
    }else{
        arr0 <- (m*(2*pi/an)^2)^(1/3)
    }
    arr <- arr0#
    Ntry <- 1e3
    for(k in 1:Ntry){
        arr.old <- arr
        if(Par$BinaryUnit=='natural'){
            arr <- arr0*(1+(m2*m1/m^2-9)*0.5*m/arr)^(2/3)
        }else{
            arr <- arr0*(1+(m2*m1/m^2-9)*0.5*4*pi^2*m/arr/YC^2)^(2/3)
        }
        if(all(abs((arr-arr.old)/arr.old)<ARRTOL)) break()
    }
    if(Par$BinaryUnit=='natural'){
        arr <- arr0*(1+(m2*m1/m^2-9)*0.5*m/arr)^(2/3)
    }else{
        arr <- arr0*(1+(m2*m1/m^2-9)*0.5*4*pi^2*m/arr/YC^2)^(2/3)
    }
    arr
}

gen_Astro2cart <- function(state){
####################################
## Transform from astrometry quantities to Cartesian coordinates in a reference frame
##
## Input:
##   state - six dimentional astrometry data
##
## Output:
##   out - x, y, z, vx, vy, vz in the equatorial coordinate system
####################################
    auyr2kms <- 4.74047
    kpcmyr2auyr <- 1e3*206265/1e6
    kpcmyr2kms <- kpcmyr2auyr*auyr2kms
    if(is.null(dim(state))){
        state <- t(state)
    }
    ra <- state[,1]/180*pi; dec <- state[,2]/180*pi; plx <- state[,3]; pmra <- state[,4]; pmdec <- state[,5]; rv <- state[,6]
    d <- 1000/plx#pc
    xyz <- gen_Lb2xyz(ra,dec)*d#pc
###convert velocity
    vde <- pmdec/plx#au/yr
    vra <- pmra/plx#au/yr
    vp <- sqrt(vra^2+vde^2)#au/yr
    vr <- rv/auyr2kms#au/yr
    vx <- vr*cos(dec)*cos(ra)-vde*sin(dec)*cos(ra)-vra*sin(ra)##note: vr is positive if the star is moving away from the Sun
    vy <- vr*cos(dec)*sin(ra)-vde*sin(dec)*sin(ra)+vra*cos(ra)
    vz <- vr*sin(dec)+vde*cos(dec)
    out <- cbind(xyz,vx,vy,vz)
    colnames(out) <- c('x.pc','y.pc','z.pc','u.auyr','v.auyr','w.auyr')
    return(out)
}
gen_Lb2xyz <- function(l.rad,b.rad){
####################################
## Transform from angular pairs to  Cartesian quantitites
## Angular pairs could be (ra,dec) or (l,b)
##
## Input:
##   state - six dimentional astrometry data
##
## Output:
##   out - x, y, z, vx, vy, vz in the equatorial coordinate system
####################################
    x <- cos(b.rad)*cos(l.rad)
    y <- cos(b.rad)*sin(l.rad)
    z <- sin(b.rad)
    return(cbind(x,y,z))
}
gen_Xyz2lb <- function(xyz){
####################################
## Transform from Cartesian quantitites to angular pairs
## Angular pairs could be (ra,dec) or (l,b)
##
## Input:
##   vec - x,y,z
##
## Output:
##   l - pair 1 (rad)
##   b - pair 2 (rad)
####################################
    x <- xyz[,1]
    y <- xyz[,2]
    z <- xyz[,3]
    b <- atan(z/sqrt(x^2+y^2))
    l <- atan(y/x)
    inds <- which(x<0)
    l[inds] <- l[inds]+pi
    inds <- which(x>=0 & y<0)
    l[inds] <- l[inds]+2*pi
    return(cbind(l,b))
}
gen_Cart2astro <- function(state){
    auyr2kms <- 4.74047
    kpcmyr2auyr <- 1e3*206265/1e6
    kpcmyr2kms <- kpcmyr2auyr*auyr2kms
    x <- state[,1]; y <- state[,2]; z <- state[,3]; vx <- state[,4]; vy <- state[,5]; vz <- state[,6]
    N <- length(x)
    ad.rad <- gen_Xyz2lb(cbind(x,y,z))
    d <- sqrt(x^2+y^2+z^2)
    ra.rad <- ad.rad[,1]
    dec.rad <- ad.rad[,2]
    ra <- ad.rad[,1]*180/pi
    dec <- ad.rad[,2]*180/pi
###velocity to pm
    v <- sqrt(vx^2+vy^2+vz^2)
    vad <- gen_Xyz2lb(cbind(vx,vy,vz))
    vxe <- v*cos(vad[,2])*cos(vad[,1])
    vye <- v*cos(vad[,2])*sin(vad[,1])
    vze <- v*sin(vad[,2])
    vequ <- array(NA,dim=c(N,3))
    for(j in 1:N){
        rotz <- matrix(data=c(cos(ra.rad[j]),sin(ra.rad[j]),0.0,-sin(ra.rad[j]),cos(ra.rad[j]),0.0,0.0,0.0,1.0),nrow=3,ncol=3,byrow=TRUE)#o-xyz -> o-x'y'z'
        roty <- matrix(data=c(cos(dec.rad[j]),0.0,sin(dec.rad[j]),0.0,1.0,0.0,-sin(dec.rad[j]),0.0,cos(dec.rad[j])),nrow=3,ncol=3,byrow=TRUE)
        vequ[j,] <- roty%*%rotz%*%c(vxe[j],vye[j],vze[j])
    }
    pm.ra <- vequ[,2]/d*1000#mas/yr
    pm.dec <- vequ[,3]/d*1000#mas/yr
    rv <- vequ[,1]*auyr2kms#km/s
    out <- cbind(ra,dec,1000/d,pm.ra,pm.dec,rv)
    colnames(out) <- c('ra','dec','parallax','pmra','pmdec','radial_velocity')
    return(out)
}
