###this file include functions for astrometry modeling
astro_deg2hdms <- function(RAdeg,DEdeg){
####################################
## Convert from deg to h/d m s
##
## Input:
##   RAdeg - RA in degree
##   DECdeg - DEC in degree
##
## Output:
##   RAh - RA hour
##   RAm - RA minute
##   RAs - RA second
##   DECh - DEC hour
##   DECm - DEC minute
##   DECs - DEC second
####################################
    val <- RAdeg/360*24#hour
    RAh <- floor(val)
    RAm <- floor((val%%1)*60)
    RAs <- (((val%%1)*60)%%1)*60
    sig <- sign(DEdeg)
    DEdeg <- abs(DEdeg)
    DEd <- sig*floor(DEdeg)
    DEm <- floor((DEdeg%%1)*60)
    DEs <- (((DEdeg%%1)*60)%%1)*60
    return(cbind(RAh,RAm,RAs,DEd,DEm,DEs))
}

astro_hdms2deg <- function(RAh,RAm,RAs,DECd,DECm,DECs){
####################################
## Convert from deg to h/d m s
##
## Input:
##   RAh - RA hour
##   RAm - RA minute
##   RAs - RA second
##   DECh - DEC hour
##   DECm - DEC minute
##   DECs - DEC second
##
## Output:
##   RAdeg - RA in degree
##   DECdeg - DEC in degree
##   RArad - RA in units of radian
##   DECrad - DEC in units of radian
####################################
    RAdeg <- (RAh+RAm/60+RAs/3600)/24*360
    DECdeg <- sign(DECd)*(abs(DECd)+DECm/60+DECs/3600)
    RArad <- RAdeg/180*pi
    DECrad <- DECdeg/180*pi
    cbind(RAdeg,DECdeg,RArad,DECrad)
}

astro_refi <- function(DN,RDNDR){
####################################
## Calculate x/(x+y)
##
## Input:
##   DN - DN
##   RDNDR - RDNDR
##
## Output:
##   RDNDR/(DN+RDNDR)
####################################
    RDNDR/(DN+RDNDR)
}

astro_drange <- function(theta){
####################################
## Convert angular quantity to the range of [-pi,pi]
##
## Input:
##   theta - Angular quantity
##
## Output:
##   theta1 - Angular quantity in the range from -pi to pi
####################################
    theta1 <- theta%%(2*pi)
    theta1[theta1>=pi] <- theta1[theta1>=pi]-2*pi
    theta1[theta1<=-pi] <- theta1[theta1<=-pi]+2*pi
    theta1
}

astro_refro <- function(zobs,hm,tdk,pmb,rh,wl,phi,tlr=0.0065,eps=1e-8){
####################################
## This routine is adapted from slaRefro.3 in P.T.Wallace   Starlink   31 October 1993
## references:
## https://github.com/Starlink/starlink/blob/master/libraries/sla/refro.f
## https://www.fact-project.org/FACT++/palRefro_8c_a2739c3f077af03fc65baeed433ed23bc.html
## https://iopscience.iop.org/article/10.1086/316172/pdf
##
## Input:
##   zobs    -  Observed zenith distance of the source (radian)
##   hm      -  Height of the observer above sea level (metre)
##   tdk     -  Ambient temperature at the observer (K)
##   pmb     -  Pressure at the observer (millibar)
##   rh     -  Relative humidity at the observer (range 0-1)
##   wl      -  Effective wavelength of the source (micrometre)
##   phi     -  Latitude of the observer (radian, astronomical)
##   tlr     -  Temperature lapse rate in the troposphere (K/metre)
##   eps     -  Precision required to terminate iteration (radian)
##
## Output:
##   Ref - Refraction angle
####################################
    ##93 degrees in radians
    d93=1.623156204
    ##Universal gas constant
    gcr=8314.32
    ##Molecular weight of dry air
    dmd=28.9644
    ##Molecular weight of water vapour
    dmw=18.0152
    ##Mean Earth radius (metre)
    s=6378120
    ##Exponent of temperature dependence of water vapour pressure
    delta=18.36
    ##Height of tropopause (metre)
    ht=11000
    ##Upper limit for refractive effects (metre)
    hs=80000
    ##Numerical integration: maximum number of strips.
    ismax=16384
    ##observation epochs
    Nobs <- max(length(zobs),length(hm),length(tdk),length(pmb),length(rh),length(wl),length(phi),length(tlr))
    ##make sure all quantities having the same vector length
    if(length(zobs)<Nobs) zobs <- rep(zobs,Nobs)
    if(length(hm)<Nobs) hm <- rep(hm,Nobs)
    if(length(tdk)<Nobs) tdk <- rep(tdk,Nobs)
    if(length(pmb)<Nobs) pmb <- rep(pmb,Nobs)
    if(length(rh)<Nobs) rh <- rep(rh,Nobs)
    if(length(wl)<Nobs) wl <- rep(wl,Nobs)
    if(length(phi)<Nobs) phi <- rep(phi,Nobs)
    if(length(tlr)<Nobs) tlr <- rep(tlr,Nobs)

    ##  Transform ZOBS into the normal range.
    zobs1  <- astro_drange(zobs)
    zobs2 <- sapply(1:length(zobs1),function(i) min(abs(zobs1[i]),d93))

    ##  Keep other arguments within safe bounds.
    hmok <- hm
    hmok[hmok< -1e-3] <- -1e-3
    hmok[hmok> hs] <- hs

    tdkok <- tdk
    tdkok[tdkok<100] <- 100
    tdkok[tdkok>500] <- 500

    pmbok <- pmb
    pmbok[pmbok<0] <- 0
    pmbok[pmbok>1e4] <- 1e4

    rhok <- rh
    rhok[rhok<0] <- 0
    rhok[rhok>1] <- 1

    wlok <- wl
    wlok[wlok<0.1] <- 0.1

    alpha <- tlr
    alpha[alpha<0.001] <- 0.001
    alpha[alpha>0.01] <- 0.01

    ##Tolerance for iteration
    tol <- min(max(abs(eps),1e-12),0.1)/2

    ##Decide whether optical/IR or radial case - switch at 100 microns
    ind.optic <- which(wlok <= 100)
    ind.radio <- which(wlok >100)

    ##  Set up model atmosphere parameters defined at the observer.
    wlsq <- wlok^2
    gb <- 9.784*(1-0.0026*cos(2*phi)-2.8e-7*hmok)
    a <- rep(77.689e-6,length(wlok))
    if(length(ind.optic)>0){
        a[ind.optic] <- (287.6155+(1.62887+0.01360/wlsq[ind.optic])/wlsq[ind.optic])*273.15e-6/1013.25
    }
    gamal <- gb*dmd/gcr
    gamma <- gamal/alpha

    gamm2 <- gamma-2
    delm2 <- delta-2
    tdc <- tdkok-273.15
    psat <- 10^((0.7859+0.03477*tdc)/(1+0.00412*tdc))*(1+pmbok*(4.5e-6+6e-10*tdc^2))
    pwo <- rep(0,length(rhok))
    ind <- which(pmbok>0)
    if(length(ind)>0){
        pwo[ind] = rhok[ind]*psat[ind]/(1.0-(1.0-rhok[ind])*psat[ind]/pmbok[ind])
    }
    w = pwo*(1.0-dmw/dmd)*gamma/(delta-gamma)
    c1 = a*(pmbok+w)/tdkok
    c2 <- rep(NA,length(pwo))
    if(length(ind.optic)>0)  c2[ind.optic] <-  (a[ind.optic]*w[ind.optic]+11.2684e-6*pwo[ind.optic])/tdkok[ind.optic]
    if(length(ind.radio)>0)  c2[ind.radio] <-  (a[ind.radio]*w[ind.radio]+6.3938e-6*pwo[ind.radio])/tdkok[ind.radio]
    c3 = (gamma-1.0)*alpha*c1/tdkok
    c4 = (delta-1.0)*alpha*c2/tdkok
    c6 <- c5 <- rep(0,length(pwo))
    if(length(ind.radio)>0){
        c5[ind.radio] = 375463e-6*pwo[ind.radio]/tdkok[ind.radio]
        c6[ind.radio] = c5[ind.radio]*delm2[ind.radio]*alpha[ind.radio]/tdkok[ind.radio]^2
    }
    ##conditions at the observer.
    r0 = s+hmok
    tmp <- astro_atmt(r0,tdkok,alpha,gamm2,delm2,c1,c2,c3,c4,c5,c6,r0)
    tempo <- tmp$T
    dn0 <- tmp$dn
    rdndr0 <- tmp$rdndr
    sk0 = dn0*r0*sin(zobs2)
    f0 = astro_refi(dn0,rdndr0)

    ##conditions in the troposphere at the tropopause.
    rt = s+sapply(1:length(hmok),function(i) max(ht,hmok[i]))
    tmp <- astro_atmt(r0,tdkok,alpha,gamm2,delm2,c1,c2,c3,c4,c5,c6,rt)
    tt <- tmp$T
    dnt <- tmp$dn
    rdndrt <- tmp$rdndr
    sine = sk0/(rt*dnt)

    pp <- sapply(1:length(sine),function(i) max(1.0-sine[i]^2,0.0))
    zt = atan2(sine,sqrt(pp))
    ft = astro_refi(dnt,rdndrt)

    ##conditions in the stratosphere at the tropopause.
    tmp <- astro_atms(rt,tt,dnt,gamal,rt)
    dnts <- tmp$dn
    rdndrp <- tmp$rdndr
    sine = sk0/(rt*dnts)
    pp <- sapply(1:length(sine),function(i) max(1.0-sine[i]^2,0.0))
    zts = atan2(sine,sqrt(pp))
    fts = astro_refi(dnts,rdndrp)

    ##conditions at the stratosphere limit.
    rs = s+hs
    tmp <- astro_atms(rt,tt,dnt,gamal,rs)
    dns <- tmp$dn
    rdndrs <- tmp$rdndr
    sine = sk0/(rs*dns)
    pp <- sapply(1:length(sine),function(i) max(1.0-sine[i]^2,0.0))
    zs = atan2(sine,sqrt(pp))
    fs = astro_refi(dns,rdndrs)

    ##variable initialization to avoid compiler warning.
    reft = 0.0

    ##integrate the refraction integral in two parts;  first in the
    ##troposphere (k=1), then in the stratosphere (k=2).

    for(k in 1:2){
        ##initialize previous refraction to ensure at least two iterations.
        refold = 1.0

        ##start off with 8 strips.
        IS = 8

        ## start z, z range, and start and end values.
        if(k==1){
            z0 = zobs2
            zrange = zt-z0
            fb = f0
            ff = ft
        }else{
            z0 = zts
            zrange = zs-z0
            fb = fts
            ff = fs
        }

        ##sums of odd and even values.
        fo = 0.0
        fe = 0.0

        ##first time through the loop we have to do every point.
        n = 1

        ##start of iteration loop (terminates at specified precision).
        loop = 1
        k3 <- 0
        while(loop){
            k3 <- k3+1
            ##strip width.
            h = zrange/IS

            ##initialize distance from earth centre for quadrature pass.
            if (k == 1) {
                r = r0
            } else {
                r = rt
            }

            ##one pass (no need to compute evens after first time).
            i <- 1
#            for(i in 1:(is-1)) {
            while(i<IS){
                i <- i+n
                ##           sine of observed zenith distance.
                sz  <- sin(z0+h*i)

                ##         find r (to the nearest metre, maximum four iterations).
                if(all(sz > 1e-20)){
                    w = sk0/sz
                    rg = r
                    dr = 1.0e6
                    j = 0
                    while(all(abs(dr) > 1.0) & j < 4){
                        j <- j+1
                        if(k==1){
                            tmp <- astro_atmt(r0,tdkok,alpha,gamm2,delm2,c1,c2,c3,c4,c5,c6,rg)
                            tg <- tmp$tg
                            dn <- tmp$dn
                            rdndr <- tmp$rdndr
                        }else{
                            tmp <- astro_atms(rt,tt,dnt,gamal,rg)
                            dn <- tmp$dn
                            rdndr <- tmp$rdndr
                        }
                        dr = (rg*dn-w)/(dn+rdndr)
                        rg = rg-dr
                    }
                    r = rg
                }

                ##           find the refractive index and integrand at r.
                if(k==1) {
                    tmp <- astro_atmt(r0,tdkok,alpha,gamm2,delm2,c1,c2,c3,c4,c5,c6,r)
                    t <- tmp$T
                    dn <- tmp$dn
                    rdndr <- tmp$rdndr
                } else {
                    tmp <- astro_atms(rt,tt,dnt,gamal,r)
                    dn <- tmp$dn
                    rdndr <- tmp$rdndr
                }
                f = astro_refi(dn,rdndr)

                ##      accumulate odd and (first time only) even values.
                if(n==1 & i%%2 == 0){
                    fe  <- fe+ f
                } else {
                    fo  <- fo+ f
                }
            }

            ##      evaluate the integrand using simpson's rule.
            refp = h*(fb+4.0*fo+2.0*fe+ff)/3.0

            ##        has the required precision been achieved (or can't be)? */
            if(all(abs(refp-refold) > tol) & IS < ismax){
                ##           no: prepare for next iteration.

                ##           save current value for convergence test.
                refold = refp

                ##           double the number of strips.
                IS  <-  2*IS

                ##           sum of all current values = sum of next pass's even values.
                fe  <-  fe + fo

                ##           prepare for new odd values.
                fo = 0.0

                ##      skip even values next time.
                n = 2
            }else{
                ## yes: save troposphere component and terminate the loop.
                if (k==1) reft = refp
                loop = 0
            }
        }
    }
    ## result.
    ref = reft+refp
    ref[zobs1<0] <- -ref[zobs1<0]
    return(ref)
}

astro_atmt <- function( r0,  T0,  alpha,  gamm2, delm2,  c1,  c2,  c3,  c4, c5,  c6,  r){
####################################
## Internal routine used by refro
## Refractive index and derivative with respect to height for the
## troposphere. The routine and comments are adapted from sla__ATMT
## https://github.com/scottransom/pyslalib/blob/master/atmt.f
##
## Input:
##   r0      -    Height of observer from centre of the Earth (metre)
##   T0      -    Temperature at the observer (K)
##   alpha   -    alpha
##   gamm2   -    gamma minus 2  ( see HMNAO paper)
##   delm2   -    delta minus 2
##   c1      -    useful term
##   c2      -    useful term
##   c3      -    useful term  (see source)
##   c4      -    useful term of sla_REFRO
##   c5      -    useful term
##   c6      -    useful term
##   r       -    current distance from the centre of the Earth (metre)
##
## Output:
##   T       -    Temperature at r (K)
##   dn      -    Refractive index at r
##   rdndr   -    r * rate the refractive index is changing at R
##
##  Note that in the optical case c5 and c6 are zero.
####################################
    T = max( min( T0 - alpha*(r-r0), 320.0), 100.0 )
    tt0 = T / T0
    tt0gm2 = tt0^gamm2
    tt0dm2 = tt0^delm2
    dn = 1.0 + ( c1 * tt0gm2 - ( c2 - c5 / T ) * tt0dm2 ) * tt0
    rdndr = r * ( -c3 * tt0gm2 + ( c4 - c6 / tt0 ) * tt0dm2 )
    rdndr1 <- -r *c3 * tt0gm2
    rdndr2 <- r *( c4 - c6 / tt0 ) * tt0dm2
    list(T=T,dn=dn,rdndr=rdndr)
}

astro_atms <-  function(rt, tt, dnt, gamal, r){
####################################
## Internal routine used by refro
## Refractive index and derivative with respect to height for the
## stratosphere.
##
## Input:
##   rt      -    height of tropopause from centre of the Earth (metre)
##   tt      -    temperature at the tropopause (K)
##   dnt     -    refractive index at the tropopause
##   gamal   -    constant of the atmospheric model = G*MD/R
##   r      -    current distance from the centre of the Earth (metre)
##
## Output:
##   dn      -    refractive index at R
##   rdndr   -    R * rate the refractive index is changing at RInternal routine used by refro
####################################
  b = gamal / tt
  w = (dnt - 1.0) * exp( -b * (r-rt) )
  dn = 1.0 + w
  rdndr = -r * b * w
  list(dn=dn,rdndr=rdndr)
}

astro_refco <- function(  hm,  tdk,  pmb,  rh, wl,  phi,  tlr=0.0065,eps=1e-8){
####################################
## Determine the constants A and B in the atmospheric refraction
## model dZ = A tan Z + B tan**3 Z.
##  Z is the "observed" zenith distance (i.e. affected by refraction)
##  and dZ is what to add to Z to give the "topocentric" (i.e. in vacuo)
##  zenith distance.
##  ref: https://iopscience.iop.org/article/10.1086/316172/pdf
##
##  Input:
##    hm      -     height of the observer above sea level (metre)
##    tdk     -     ambient temperature at the observer (K)
##    pmb     -     pressure at the observer (millibar)
##    rh      -     relative humidity at the observer (range 0-1)
##    wl      -     effective wavelength of the source (micrometre)
##    phi     -     latitude of the observer (radian, astronomical)
##    tlr     -     temperature lapse rate in the troposphere (K/metre)
##    eps     -     precision required to terminate iteration (radian)
##
##  Output:
##    refa    -     tan Z coefficient (radian)
##    refb    -     tan**3 Z coefficient (radian)
##
##  Called:  refro
####################################
    ##  Sample zenith distances: arctan(1) and arctan(4)
    ATN1 <-  0.7853981633974483
    ATN4 <-  1.325817663668033

    ##Determine refraction for the two sample zenith distances
    r1 <- astro_refro(ATN1,hm,tdk,pmb,rh,wl,phi,tlr,eps)
    r2 <- astro_refro(ATN4,hm,tdk,pmb,rh,wl,phi,tlr,eps)

    ##Solve for refraction constants
    refa = (64.0*r1-r2)/60.0
    refb = (r2-4.0*r1)/60.0
    list(refa=refa,refb=refb)
}

astro_refcoq <- function(tdk,pmb,rh,wl){
####################################
## Determine the constants A and B in the atmospheric refraction
## model dZ = A tan Z + B tan**3 Z.  This is a fast alternative
## to the sla_REFCO routine - see notes.
## Z is the "observed" zenith distance (i.e. affected by refraction)
## and dZ is what to add to Z to give the "topocentric" (i.e. in vacuo)
## zenith distance.
## ref: https://github.com/Starlink/starlink/blob/master/libraries/sla/refcoq.f
##
## Input:
##   tdk      -      ambient temperature at the observer (K)
##   pmb      -      pressure at the observer (millibar)
##   rh       -      relative humidity at the observer (range 0-1)
##   wl       -      effective wavelength of the source (micrometre)
##
## Output:
##   refa     -      tan Z coefficient (radian)
##   refb     -      tan**3 Z coefficient (radian)
##
## Note:
##   The radio refraction is chosen by specifying WL > 100 micrometres.
####################################
    Nobs <- max(length(tdk),length(pmb),length(rh),length(wl))
    ##make sure all quantities having the same vector length
    if(length(tdk)<Nobs) tdk <- rep(tdk,Nobs)
    if(length(pmb)<Nobs) pmb <- rep(pmb,Nobs)
    if(length(rh)<Nobs) rh <- rep(rh,Nobs)
    if(length(wl)<Nobs) wl <- rep(wl,Nobs)

    ##  Decide whether optical/IR or radio case:  switch at 100 microns.
    ind.optic <- which(wl <= 100)
    ind.radio <- which(wl > 100)

    ## Restrict parameters to safe values.
    T <- tdk
    T[T<100] <-100
    T[T>500] <-500
    p  <- pmb
    p[p<0] <- 0
    p[p>1e4] <- 1e4
    r <- rh
    r[r<0] <- 0
    r[r>1] <- 1
    w <- wl
    w[w<0.1] <- 0.1
    w[w>1e6] <- 1e6

    ##  Water vapour pressure at the observer.
    pw <- rep(0,length(p))
    ind <- which(p>0)
    tdc = T[ind]-273.15
    ps = 10^((0.7859+0.03477*tdc)/(1+0.00412*tdc))*(1+p[ind]*(4.5e-6+6e-10*tdc^2))
    pw[ind] <- r[ind]*ps/(1-(1-r[ind])*ps/p[ind])

    ##Refractive index minus 1 at the observer.
    wlsq  <-  w^2
    gamma = ((77.53484e-6+(4.39108e-7+3.666e-9/wlsq)/wlsq)*p -11.2684e-6*pw)/T
    if(length(ind.radio)>0) gamma[ind.radio] <- (77.6890e-6*p[ind.radio]-(6.3938e-6-0.375463e0/T[ind.radio])*pw[ind.radio])/T[ind.radio]

    ##Formula for beta adapted from Stone, with empirical adjustments.
    beta <- 4.4474e-6*T
    if(length(ind.radio)>0){
        beta[ind.radio] <- beta[ind.radio]-0.0074*pw[ind.radio]*beta[ind.radio]
    }
    ##Refraction constants from Green.
    refa = gamma*(1-beta)
    refb = -gamma*(beta-gamma/2)

    list(refa=refa,refb=refb)
}


astro_NMFhydro <- function(mjd.utc,phi,H,elevation,delevation){
####################################
## Determine the NMF map for the calculation of hydrostatic delay
##
## Input:
##   mjd.utc      -      2-part MJD[UTC]
##   phi      - Observation site latitude
##   H       - Observation site height
##   elevation       - Observation site elevation
##   delevation       - Time derivative of observation site elevation
##
## Output:
##   map     -      Mapping value for dry (or hydrostatic) component
##   dmap     -     Time derivative of map
##
## Note:
##   defaul is to ignore wet delays
##   Zenith delay from Davies et al (Radio Sci 20 1593)
##   tempo2 typo: only use phi rather than 2*phi in cos()
##   zenith.delay.hydrostatic  <-  0.02268 * pressure/(SPEED.LIGHT*(1.0-0.00266*cos(obs->latitude.grs80)-2.8e-7*obs->height.grs80))
####################################
    avgs.a <- c(1.2769934e-3, 1.2683230e-3, 1.2465397e-3, 1.2196049e-3, 1.2045996e-3)
    avgs.b <- c(2.9153695e-3,
                2.9152299e-3,
                2.9288445e-3,
                2.9022565e-3,
                2.9024912e-3)
    avgs.c <- c(62.610505e-3,
                62.837393e-3,
                63.721774e-3,
                63.824265e-3,
                64.258455e-3)
    amps.a <- c(0.0,
                1.2709626e-5,
                2.6523662e-5,
                3.4000452e-5,
                4.1202191e-5)
    amps.b <- c(0.0,
                2.1414979e-5,
                3.0160779e-5,
                7.2562722e-5,
                11.723375e-5)
    amps.c <- c(0.0,
                0.0128400e-5,
                4.3497037e-5,
                84.795348e-5,
                170.37206e-5)
    ##Height correction coefficients
    ah <- 2.53e-5
    bh <- 5.49e-3
    ch <- 1.14e-3

###work out which latitutde to use
  abs.lat <- abs(phi)*180/pi
  ilat1  <-  floor(abs.lat / 15.0)-1
  frac  <-  abs.lat/15.0 - 1.0 - ilat1
  ilat2 <- ilat1+1

  ind1 <- which(ilat1<0)
  if(length(ind1)>0){
      ilat1[ind1] <- ilat2[ind1]  <-  0
  }

  ind2 <- which(ind1>=4)
  if(length(ind2)>0){
      ilat1[ind2]  <-  ilat2[ind2]  <-  4
  }
  ilat1 <- ilat1+1#index adapted from C to R
  ilat2 <- ilat2+1#index adapted from C to R

###work out the phase wrt DOY 28 of 2005 (or any other year), mjd 53398
    cos.phase = cos((mjd.utc - 53398.0)*2*pi/365.25)*sign(phi)

    frac <- 0
    sin.elevation <- sin(elevation)
    dsin.elevation <- cos(elevation)*delevation
###get interpolated values for a, b, c
    a  <- frac*avgs.a[ilat2] + (1.0-frac)*avgs.a[ilat1] - (frac*amps.a[ilat2] + (1.0-frac)*amps.a[ilat1]) * cos.phase
    b  <-  frac*avgs.b[ilat2] + (1.0-frac)*avgs.b[ilat1] - (frac*amps.b[ilat2] + (1.0-frac)*amps.b[ilat1]) * cos.phase
    c  <-  frac*avgs.c[ilat2] + (1.0-frac)*avgs.c[ilat1] - (frac*amps.c[ilat2] + (1.0-frac)*amps.c[ilat1]) * cos.phase

    ind <- which(ilat1==ilat2)
    if(length(ind)>0){
        a[ind] <-  avgs.a[ilat1[ind]] - amps.a[ilat1[ind]]*cos.phase[ind]
        b[ind] <-  avgs.b[ilat1[ind]] - amps.b[ilat1[ind]]*cos.phase[ind]
        c[ind] <-  avgs.c[ilat1[ind]] - amps.c[ilat1[ind]]*cos.phase[ind]
    }

###get the basic mapping function value (no height correction)
    basic <- (1.0+a/(1.0+b/(1.0+c)))/(sin.elevation+a/(sin.elevation+b/(sin.elevation+c)))
    dt <- 0.001#day
    sin.elevation1 <- sin.elevation+dt*dsin.elevation
    dbasic <- ((1.0+a/(1.0+b/(1.0+c)))/(sin.elevation1+a/(sin.elevation1+b/(sin.elevation1+c))) - basic)/(dt*DAYSEC)

###
    height.correction <- H*1.0e-3*(1.0/sin.elevation -(1.0+ah/(1.0+bh/(1.0+ch)))/(sin.elevation+ah/(sin.elevation+bh/(sin.elevation+ch))))
    dheight.correction <- (H*1.0e-3*(1.0/sin.elevation1 -(1.0+ah/(1.0+bh/(1.0+ch)))/(sin.elevation1+ah/(sin.elevation1+bh/(sin.elevation1+ch))))-height.correction)/(dt*DAYSEC)
    list(map=basic+height.correction,dmap=dbasic+dheight.correction)
}

NMF.wet <- function(phi,elevation){
####################################
## Determine the NMF map for the calculation of hydrostatic delay
##
## Input:
##   phi     -  Observation site latitude
##   elevation   - Observation site elevation
##
## Output:
##   map     -  Mapping function for wet component
####################################
  avgs.a <- c(5.8021897e-4,
        5.6794847e-4,
        5.8118019e-4,
        5.9727542e-4,
        6.1641693e-4)
    avgs.b <- c(1.4275268e-3,
        1.5138625e-3,
        1.4572752e-3,
        1.5007428e-3,
        1.7599082e-3)
    avgs.c <- c(4.3472961e-2,
        4.6729510e-2,
        4.3908931e-2,
        4.4626982e-2,
        5.4736038e-2)
  frac <- 0
  sin.elevation <- sin(elevation)

  abs.lat <- abs(phi)*180/pi
  ilat1  <-  floor(abs.lat / 15.0)-1
  frac  <-  abs.lat/15.0 - 1.0 - ilat1
  ilat2 <- ilat1+1

  ind1 <- which(ilat1<0)
  if(length(ind1)>0){
      ilat1[ind1] <- ilat2[ind1]  <-  0
  }

  ind2 <- which(ind1>=4)
  if(length(ind2)>0){
      ilat1[ind2]  <-  ilat2[ind2]  <-  4
  }
  ilat1 <- ilat1+1#adapt from C to R
  ilat2 <- ilat2+1#adapt from C to R

  if(ilat1==ilat2){
        a = as[ilat1]
        b = bs[ilat1]
        c = cs[ilat1]
  }else{
        a = frac*as[ilat2] + (1.0-frac)*as[ilat1]
        b = frac*bs[ilat2] + (1.0-frac)*bs[ilat1]
        c = frac*bs[ilat2] + (1.0-frac)*cs[ilat1]
    }
  (1+a/(1+b/(1+c)))/(sin.elevation+a/(sin.elevation+b/(sin.elevation+c)))
}


astro_aberration <- function(lo,vSO,RSO,MA=1,g=1){
####################################
## Calculate the stellar aberration
## Aberration effect
## ref:KLIONER 2003
##
## Input:
##   lo - Direction of light ray before observed
##   vSO - Velocity of the observer relative to the SSB (au/yr)
##   RSO - Distance from the observer relative to the SSB (au)
##   MA - Mass of the body which impose gravitationl force onto the observer (Msolar)
##   g - gamma; g=1 for GR; g!=1 for other gravitational theories.
##
## Output:
##   uo - Direction of the light source from the observer's perspective
##   dl1 - first order aberration vector
##   dl2 - second order aberration vector
##   dl3 - third order aberration vector
##   dl - total aberration vector
####################################
    lvl <- pracma::cross(lo,pracma::cross(vSO,lo))
    vlv <- pracma::cross(vSO,pracma::cross(lo,vSO))
    lv <- rowSums(lo*vSO)
    w <- 4*pi^2*MA/RSO#au/yr
    p1 <- lvl
    p2 <- lv*lvl+0.5*vlv
    p3 <- (lv^2+(1+g)*w)*lvl+0.5*lv*vlv
    dl1 <- p1/Cauyr
    dl2 <- p2/Cauyr^2
    dl3 <- p3/Cauyr^3
    dl <- dl1+dl2+dl3
    uo <- -lo+dl
    list(uo=uo/gen_CalLen(uo),dl1=dl1,dl2=dl2,dl3=dl3,dl=dl)
}

astro_CalSB <- function(tB,Par){
####################################
## Solve the Kepler's Equation
##
## Input:
##   tB - Coordinate time at the TSB
##   Par - Input parameters
##
## Output:
##   SB - Position and velocity vectors from the SSB to TSB (pc; au/yr)
##   RSB - Distance from SSB to TSB (pc)
##   dRSB - Difference between the RSB at time t and the RSB at time tpos (pc)
##   uSB.T2 - uSB calculated using the TEMPO2 method
####################################
    d0 <- 1/Par$plx#kpc
    pqu <- Par$pqu
#    cat('p=',pqu[,1],'\n')
#    cat('q=',pqu[,2],'\n')
#    cat('u=',pqu[,3],'\n')
#    cat('Par$pmra=',Par$pmra,'\n')
#    cat('Par$pmdec=',Par$pmdec,'\n')
#    cat('Par$rv=',Par$rv,'\n')
    v1 <- c(d0*Par$pmra,d0*Par$pmdec,Par$rv/auyr2kms)#au/yr in p,q,u system, assuming zero galactic acceleration, vSB(t)=vSB(t0)
    vSB <- as.numeric(pqu[,1]*v1[1]+pqu[,2]*v1[2]+pqu[,3]*v1[3])#au/yr
    dt <- as.numeric((tB[,1]-Par$tpos)+(tB[,2]))/DJY#year
#    if(UNITS=='TDB'){
#        dt <- dt/IFTE.K
#    }
    drSB <- cbind(dt*vSB[1],dt*vSB[2],dt*vSB[3])/pc2au#pc
    rSB0 <- d0*1e3*c(cos(Par$dec)*cos(Par$ra),cos(Par$dec)*sin(Par$ra),sin(Par$dec))#pc
    rSB <- cbind(drSB[,1]+rSB0[1],drSB[,2]+rSB0[2],drSB[,3]+rSB0[3])
#    if(UNITS=='TDB'){
#        rSB <- rSB/IFTE.K
#        rSB0 <- rSB0/IFTE.K
#    }
    RSB <- sqrt(rowSums(rSB^2))
    RSB0 <- sqrt(sum(rSB0^2))
#    SB <- cbind(asNumeric(rSB),t(replicate(nrow(t),vSB)))
    SB <- cbind(rSB,t(replicate(nrow(tB),vSB)))
###compared with the TEMPO2 approach
    if(Par$CompareT2){
        mu.perp <- Par$pmra*Par$pqu[,1]+Par$pmdec*Par$pqu[,2]/IFTE.K#mas/yr; from ephemeris to coordinate
        mu.para <- Par$rv/auyr2kms/d0/IFTE.K#mas/yr
        u1 <- outer(dt,mu.perp*DMAS2R,FUN='*')#rad
        u2 <- -outer(dt^2*DMAS2R^2,0.5*sum(mu.perp^2)*Par$pqu[,3],'*')#rad
        u3 <- -outer(dt^2*DMAS2R^2,mu.para*mu.perp,'*')#rad
        u0 <- t(replicate(nrow(tB),Par$pqu[,3]))
        du1 <- t(replicate(length(dt),mu.perp*DMAS2R))#rad/yr
        du2 <- -outer(2*dt*DMAS2R^2,0.5*sum(mu.perp^2)*Par$pqu[,3],'*')#rad/yr
        du3 <- -outer(2*dt*DMAS2R^2,mu.para*mu.perp,'*')#rad/yr
        uSBt <- u0+u1+u2+u3
        uSB.T2 <- list(uSBt=uSBt,u0=u0,u1=u1,u2=u2,u3=u3,du1=du1,du2=du2,du3=du3)
    }else{
        uSB.T2 <- NULL
    }
    ##uSBt <- uSBt/sqrt(rowSums(uSB1^2))
    list(SB=SB,uSB=gen_CalUnit(SB[,1:3]),rSB=rSB,RSB=RSB,dRSB=RSB-RSB0,uSB.T2=uSB.T2)
}
astro_LenSolar <- function(l1,OL,g=1){
####################################
## Calculate lensing by Solar System bodies
##
## Input:
##   l1 - Direction (unit vector) of the light ray before entering the Solar System
##   OL - Position and velocity vectors from the observer to the lens
##   g - gamma; g=1 for GR; g!=1 for other gravitational theories.
##
## Output:
##   l - Light ray direction after gravitational lensing
##   Def - Deflection vector due to lensing
##   DefList - Deflection angles of individual lenses
####################################
    ns <- names(OL)
    DefList <- list()
    Def <- 0
    for(n in ns){
        sr <- SRS*Mssp[n]*(1+g)/2#au
        dL <- pracma::cross(l1,pracma::cross(OL[[n]][,c('x.au','y.au','z.au')],l1))
        DL <- sqrt(rowSums(dL^2))
        dl <- -(sr*dL/(DL^2))*(1+OL[[n]][,'cpsi'])#rad
        Def <- Def+dl
        DefList[[n]] <- dl
    }
    l2 <- l1+Def
    list(l=l2,Def=Def,DefList=DefList)
}
astro_LenTarget <- function(l1,MA,rOT,rOA,rTA,g=1){
####################################
## Calculate lensing by Solar System bodies
##
## Input:
##   l1 - Direction of the light ray before entering the Solar System
##   OL - Position and velocity vectors from the observer to the lens
##   g - gamma; g=1 for GR; g!=1 for other gravitational theories.
##
## Output:
##   l - Light ray direction after gravitational lensing
##   dl - Deflection angle due to lensing
####################################
    pp <- (1+g)*4*pi^2*MA/Cauyr^2
    RTA <- sqrt(rowSums(rTA^2))
    ROT <- sqrt(rowSums(rOT^2))
    ROA <- sqrt(rowSums(rOA^2))
    oto <- pracma::cross(as.matrix(rOT),pracma::cross(as.matrix(rTA),as.matrix(rOA)))
    oaa <- ROT*ROA*(RTA*ROA+rowSums(rOA*rTA))
    dl <- -pp*oto/oaa#vector
    list(l=l1+dl,dl=dl)
}
astro_Relative <- function(rBT,rSO,tS,Par){
####################################
## Calculate relative astrometry
##
## Input:
##   BT - Position and velocity vectors from the TSB to the target
##   SO - Position and velocity vectors from the SSB to the observer
##   tS - Coordinate time at the SSB
##   Par - Input parameters
##
## Output:
##   xi1 - 1st order offset in the direction of increasing right ascension or x direction
##   eta1 - 1st order offset in the direction of increasing decension or y direction
##   xi2 - 2nd order offset in +x
##   eta2 - 2nd order offset in +y
##   xi3 - 3rd order offset in +x
##   eta3 - 3rd order offset in +y
##   xi - Full offset in +x
##   eta - Full offset in +y
####################################
    dt <- ((tS[,1]-Par$tpos)+tS[,2])/DJY
    pmra <- Par$pmra
    pmdec <- Par$pmdec
    plx <- Par$plx
    pmrv <- Par$rv/auyr2kms*plx#mas/yr
    R <- (rBT-rSO)*plx
    p <- Par$pqu[,1]
    q <- Par$pqu[,2]
    u <- Par$pqu[,3]
    pR <- R[,1]*p[1]+R[,2]*p[2]+R[,3]*p[3]
    qR <- R[,1]*q[1]+R[,2]*q[2]+R[,3]*q[3]
    uR <- R[,1]*u[1]+R[,2]*u[2]+R[,3]*u[3]
    murat <-pmra*dt
    mudect <-pmdec*dt
    murart2 <- pmra*pmrv*dt^2*DMAS2R
    mudecrt2 <- pmdec*pmrv*dt^2*DMAS2R
    upR <- (pmra*uR+pmrv*pR)*dt*DMAS2R
    uqR <- (pmdec*uR+pmrv*qR)*dt*DMAS2R
    puR <- pR*uR*DMAS2R
    quR <- qR*uR*DMAS2R
    xi1 <- murat
    eta1 <- mudect
    xi2 <- pR-murart2-upR-puR
    eta2 <- qR-mudecrt2-uqR-quR
    xi12 <- xi1+xi2
    eta12 <- eta1+eta2
    xi3 <- (xi12)^3/6*DMAS2R^2+pmra*pmrv^2*dt^3*DMAS2R^2
    eta3 <- (eta12)^3/6*DMAS2R^2+pmdec*pmrv^2*dt^3*DMAS2R^2
    return(list(xi1=xi1,eta1=eta1,xi2=xi2,eta2=eta2,xi3=xi3,eta3=eta3,xi=xi1+xi2+xi3,eta=eta1+eta2+eta3))
}

astro_FullModel <- function(OutBary,OutTime,Par,Mlens=1,component='T'){
####################################
## Calculate relative astrometry
##
## Input:
##   OutBary - Output of time_Utc2tb
##   OutTime - Output of time_Ta2te
##   Input - Input parameter
##   Mlens - Lense mass
##   component - target star (T) or companion star (C)
##
## Output:
#    list(dir=dir,uo=uo,lo=lo,li=li,ll=ll,all=dth.all,abe=dth.abe,abe1=dth.abe1,abe2=dth.abe2,lenT=dth.lenT,ref=dth.ref,R=Ref,Rvec=Ref.vec,lensing=lensing)
##   DirObs - Observed direction of the target in RA and DEC in ICRS
##   uo - Observed direction of the target in cartesian coordinates in the ICRS frame
##   lo - Light ray direction before aberration and after atmospheric refraction
##   li - Light ray before atmospheric refraction and after Solar System lensing
##   ll - Light ray before Solar System lensing and after leaving the target system
##   OffAll - Offsets due to all effects
##   OffAbe - Offsets due to stellar aberration
##   OffAbe1 - Offsets due to 1st order stellar aberration
##   OffAbe2 - Offsets due to 2nd order stellar aberration
##   OffLenT - Offsets due to target system lensing
##   OffLenS - Offsets due to lensing effects in the Solar System
##   DefList - Lensing deflection angle due to various objects in the solar system
##   OffRef - Offsets due to atmospheric refraction
##   Ref - Scalar atmospheric refraction
##   ref - Vector atmospheric refraction
##   OffLensing - Offsets due to gravitational lensing by Solar System bodies
####################################
    uOT <- OutTime$uOT
    rOT <- OutTime$rOT
    vST <- OutTime$vST#au/yr
    rST <- OutTime$rST#pc
    rBT <- OutTime$BT[,1:3]
    rTC <- OutTime$rTC
    rOC <- OutTime$rOC
    uOC <- OutTime$uOC
    rSC <- OutTime$rSC
    vSC <- OutTime$vSC
    SO <- OutBary$SO
    zenith <- OutBary$zenith

###direction of emission light
    if(component=='T'){
        le <- -uOT#approximation of sigma
    }else{
        le <- -uOC
    }

##lensing in the target system
###PPN parameter; GR (g=1)
    if(is.null(Mlens)) Mlens <- 0
    if(Mlens>0 & Par$binary & Par$Np>0){
        if(component=='T'){
            tmp <- astro_LenTarget(l1=le,MA=Mlens,rOT=rOT*pc2au,rOA=rOC*pc2au,rTA=rTC,g=Par$g)
        }else{
            tmp <- astro_LenTarget(l1=le,MA=Mlens,rOT=rOC*pc2au,rOA=rOT*pc2au,rTA=-rTC,g=Par$g)
        }
        ll <-tmp$l
        llmle <- tmp$dl
    }else{
        ll <- le
        llmle <- 0
    }

###lensing in the solar system
    LenSolar <- astro_LenSolar(l1=ll,OL=OutTime$OL, g=Par$g)
    li <- LenSolar$l
    limll <- LenSolar$Def

####atmospheric refraction
                                        #    elevation[elevation<0] <- -elevation[elevation<0]
    elevation <- rep(NA,Par$Nepoch)
    Ref <- rep(0,Par$Nepoch)
    ref <- cbind(Ref,Ref,Ref)
    lo <- li
    indG <- which(Par$ObsInfo[,'ObsType']=='ground' & Par$ObsInfo[,'RefType']!='none')
    if(length(indG)>0){
        if(component=='T'){
            elevation <- asin(rowSums(zenith*uOT))
        }else{
            elevation <- asin(rowSums(zenith*uOC))
        }
        zen <- pi/2-elevation#angle
        tmp <- gen_refraction(zenith=OutBary$zenith,ZenIn=zen,Par)
        Ref <- tmp$R
        if(component=='T'){
            ref <- (zenith-rowSums(zenith*uOT)*uOT)/sin(zen)*Ref
        }else{
            ref <- (zenith-rowSums(zenith*uOC)*uOC)/sin(zen)*Ref
        }
        lo <- li-ref
        lo <- lo/sqrt(rowSums(lo^2))
    }
    lomli <- lo-li

###observed direction after aberration
    rSO <- SO[,1:3,drop=FALSE]
    RSO <- sqrt(rowSums(rSO^2))
    abe <- astro_aberration(lo,vSO=SO[,4:6,drop=FALSE]/auyr2kms,RSO=RSO,MA=1,g=Par$g)#observed direction
    uo <- abe$uo
    uommlo <- abe$dl
    uommlo1 <- abe$dl1
    uommlo2 <- abe$dl2
    uommlo3 <- abe$dl3
    dir <- gen_Xyz2lb(uo)

##Offsets; Note that the u and -u directions have opposite right ascension but have the same declination. So (dra(du),ddec(du)) would be different from -(dra(-du),ddec(-du)); we need to convert all light directions to their opposite directions (i.e. from observer to the source) in order to get the correct offsets
##astrometry modeling sequence: uOT -> le -> ll -> li -> lo -> uo (typically uOT=-le)
##backward modeling: uo -> -lo -> -li -> -ll -> -le -> uOT (typically uOT=-le)
    delta <- gen_Xyz2lb(uOT)[,2]
    dl.all <- uommlo-lomli-limll-llmle#uo-uOT
    dl.woRef <- uommlo-limll-llmle
##angle-based method
    OffAll <- gen_CalOffset(uo,uOT,bref=delta)
    OffAbe <- gen_CalOffset(uo,-lo,bref=delta)
    OffAbe1 <- gen_CalOffset(-lo+uommlo1,-lo,bref=delta)
    OffAbe2 <- gen_CalOffset(-lo+uommlo2,-lo,bref=delta)
    OffAbe3 <- gen_CalOffset(-lo+uommlo3,-lo,bref=delta)
    OffRef <- gen_CalOffset(-lo,-li,bref=delta)
    OffLenS <- gen_CalOffset(-li,-ll,bref=delta)
    OffLenT <- gen_CalOffset(-ll,uOT,bref=delta)
##
    list(DirObs=dir,uo=uo,lo=lo,li=li,ll=ll,OffAll=OffAll,OffAbe=OffAbe,OffAbe1=OffAbe1,OffAbe2=OffAbe2,OffLenT=OffLenT,OffLenS=OffLenS,SolarDef=LenSolar$Def,SolarDefList=LenSolar$DefList,OffRef=OffRef,Ref=Ref,ref=ref,uommlo=uommlo,uommlo1=uommlo1,uommlo2=uommlo2,uommlo3=uommlo3,lomli=lomli,limll=limll,llmle=llmle,dl.all=dl.all,dl.woRef=dl.woRef)
}

astro_CalElevation <- function(zenith,dzenith,uOT,verbose=FALSE){
####################################
## Calculate elevation angle of the target
##
## Input:
##   zenith - Unit vector along the zenith direction
##   dzenith - Time derivative of dzenith
##   uOT - Unit vector from the observer to the target
##
## Output:
##   elevation - Elevation angle (rad)
##   delevation - Time derivative of the elevation angle (rad)
##   ZenIn - Zenith angle of the incident light
####################################
    elevation <- asin(rowSums(zenith*uOT))
    delevation <- 1/sqrt(1-(rowSums(zenith*uOT))^2)*(rowSums(dzenith*uOT))
    ind <- which(elevation<0)
    if(length(ind)>0){
        if(verbose) cat('elevation angle < 0 for ',length(ind),'UTC epochs!\n')
    }
#    elevation <- abs(elevation)
    ZenIn <- pi/2-elevation
    ZenIn[ZenIn==0] <- 1e-10
    list(elevation=elevation,delevation=delevation,ZenIn=ZenIn)
}

astro_CalRefraction <- function(zenith,ZenIn,uOT,Par){
####################################
## Calculate atmospheric refraction
##
## Input:
##   zenith - Unit vector from geocenter to the observer
##   ZenIn - Zenith angle of the incident light
##   uOT - Unit vector from the observer to the target
##   Par - Input parameter
##   RefType - Calculation method for atmospheric refraction; 'refro' or 'refroq' or 'NA'(no atmospheric refraction is calculated); "refro" is the default
##
## Output:
##   Ref - Scalar refraction
##   ref - Vector refraction
##   RvRef - Radial velocity variation due to refraction
####################################
    if(Par$RefType!='' & Par$ObsType=='ground'){
        OutRef <- gen_refraction(zenith,ZenIn,Par)
        Ref <- OutRef$R
        ref <- (zenith-rowSums(zenith*uOT)*uOT)/sin(ZenIn)*Ratm
        RvRef <- rowSums(Ratm.vec*OS1[,4:6,drop=FALSE])*auyr2kms*1e3#m/s
    }else{
        Ref <- ref <- RvRef <- 0
    }
    list(Ref=Ref,ref=ref,RvRef=RvRef)
}

astro_RaDec2pqu <- function(ra,dec){
####################################
## Convert from RA and DEC to [p q u] triad
##
## Input:
##   ra - Right ascension (rad)
##   dec - Declination (rad)
##
## Output:
##  p - Direction of increasing RA
##  q - Direction of increasing DEC
##  u - Direction of increasing distance
####################################
    p <- as.numeric(c(-sin(ra),cos(ra),0))
    q <- as.numeric(c(-sin(dec)*cos(ra),-sin(dec)*sin(ra),cos(dec)))
    u <- as.numeric(c(cos(dec)*cos(ra),cos(dec)*sin(ra),sin(dec)))
    cbind(p,q,u)
}
