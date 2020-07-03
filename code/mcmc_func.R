###In this new file, I have chagned the calcuation memthod of mu1 and correct a bug in the expression of Vr.ma;
####add celerite SHO kernel
library(fields)
library(MASS)
library(foreach)
library(doMC)
library(Matrix)
library(doParallel)
library(parallel)
library(kernlab)
library(glasso)
library(JPEN)
library(matrixcalc)
library(mvtnorm)
library(MASS)
library(rootSolve)
tol1 <- 1e-16
###modeL Of Keplerian Motions Of planets
###According to http://w.astro.berkeley.edu/~kclubb/pdf/RV_Derivation.pdf and Ford 2006
###The mass is from Berger, D. H.; et al. (2006). "First Results from the CHARA Array. IV. The Interferometric Radii of Low-Mass Stars".
###Ms=0.404
###assign names to the parameter vectors of a RV model
#Naming system: Keplerian parameters: {{per},{K},{e},{omega},{Mo}}; trend pars: {a, b}; noise pars (white noise: {s}, Gaussian process red noise: {sigma.white,l}, ARMA red noise:{{phi},alpha,{w},beta}>)
##kepler solver
kep.mv <- function(m,e){
    tol = 1e-6
    m <- m%%(2*pi)
    E0 <- m+e*sin(m)+e^2*sin(m)*cos(m)+0.5*e^3*sin(m)*(3*cos(m)^2-1)#initial value
    Ntt <- 1e4
    for(k in 1:Ntt){
        t1 <- cos(E0)
        t2 <- -1+e*t1
        t3 <- sin(E0)
        t4 <- e*t3
        t5 <- -E0+t4+m
        t6 <- t5/(0.5*t5*t4/t2+t2)
        E1 = E0-t5/((1/2*t3 - 1/6*t1*t6)*e*t6+t2)
        if(all(abs(E1-E0)<tol)) break()
        if(k==Ntt) cat('Keplerian solver does not converge!\n')
        E0 <- E1
    }
    return(E1%%(2*pi))
}
deg2hdms <- function(RAdeg,DEdeg){
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

hdms2deg <- function(rah,ram,ras,ded,dem,des){
    RA.deg <- (rah+ram/60+ras/3600)/24*360
    DEC.deg <- sign(ded)*(abs(ded)+dem/60+des/3600)
    RA.rad <- RA.deg/180*pi
    DEC.rad <- DEC.deg/180*pi
    cbind(RA.deg,DEC.deg,RA.rad,DEC.rad)
}
kep.mt <- function(m,e){
    tol = 1e-8
    E <- rep(NA,length(m))
    Ntt <- 1000
    for(j in 1:length(m)){
        E0 <- m[j]
        for(k in 1:Ntt){
            E1 = E0-(E0-e*sin(E0)-m[j])/(sqrt((1-e*cos(E0))^2-(E0-e*sin(E0)-m[j])*(e*sin(E0))))
            if(abs(E1-E0)<tol) break()
            if(j==Ntt) cat('Keplerian solver does not converge!\n')
            E0 <- E1
        }
        E[j] <- E1%%(2*pi)
    }
    return(E)
}
getphase <- function(e,omega,type='primary'){
##e: eccentricity
##omega: argument of ascending node
##T: true anomaly
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
getM0 <- function(e,omega,P,T,T0,type='primary'){
    Tp <- T-getphase(e,omega)*P
    ((T0-Tp)%%P)*2*pi/P
}
kep.mt2 <- function(m,e){
    tol = 1e-8
    E0 <- m
    Ntt <- 1e3
    for(k in 1:Ntt){
        E1 = E0-(E0-e*sin(E0)-m)/(sqrt((1-e*cos(E0))^2-(E0-e*sin(E0)-m)*(e*sin(E0))))
        if(all(abs(E1-E0)<tol)) break()
#        if(k==Ntt) cat('Keplerian solver does not converge:',e,m,E0,E1,'!\n')
        E0 <- E1
    }
    if(k==Ntt){
        cat('Keplerian solver does not converge!\n')
        cat('length(which(abs(E1-E0)>tol))=',length(which(abs(E1-E0)>tol)),'\n')
    }
    return(E1)
}
##refer to Murison A Practical Method for Solving the Kepler Equation
eps3 <- function(m,e,x){
    t1 <- cos(x)
    t2 <- -1+e*t1
    t3 <- sin(x)
    t4 <- e*t3
    t5 <- -x+t4+m
    t6 <- t5/(1/2*t5*t4/t2+t2)
    return(t5/((1/2*t3-1/6*t1*t6)*e*t6+t2))
}
KeplerStart3 <- function(m,e){
    t34 <- e^2
    t35 <- e*t34
    t33 <- cos(m)
    return(m+(-1/2*t35+e+(t34+3/2*t33*t35)*t33)*sin(m))
}
kep.murison2 <- function(m,e,tol=1e-6){
    Mnorm <- m%%(2*pi)
    E0 <- KeplerStart3(Mnorm,e)
    Ntry <- 1000
    for(k in 1:Ntry){
        E1 <- E0-eps3(Mnorm,e,E0)
        if(k==Ntry) 'Kepler solver failed to converge!\n'
        if(all(abs(E1-E0)<tol)) break()
        E0 <- E1
    }
    return(E1)
}
##R package: uniroot
kep.R <- function(m,e){
    E <- rep(NA,length(m))
    for(j in 1:length(m)){
        y <- function(x) x-e*sin(x)-m[j]%%(2*pi)
        E[j] <- uniroot(y,interval=c(0,2*pi))$root
    }
    return(E)
}
divide.pars <- function(pars){
    #select Keplerian pars
    if(Np>0){
        ind.kep <- 1:(Nkeppar*Np)
    }else{
        ind.kep <- c()
    }
    #select noise parameters for different aperture
#    ind.noise <- c()
    ind.trend <- sort((1:Npar)[grepl('\\da',names(pars)) | grepl('\\db',names(pars)) | grepl('\\ds',names(pars))])
    ind.noise <- (1:Npar)[-c(ind.kep,ind.trend)]
    return(list(kep=ind.kep,trend=ind.trend,noise=ind.noise))
}

assign.names <- function(par.vector,Np,p=0,q=0,n=0,basis='natural'){
    nams <- c()
    #names for keplerian part
    if(Np>0){
        if(basis=='natural'){
            if(prior.type=='e0'){
                nams.kep <- c('per','K','omega')
            }else{
                nams.kep <- c('per','K','e','omega','Mo')
            }
        }else if(basis=='linear1'){
            if(prior.type=='e0'){
                nams.kep <- c('per','lnK','Tc')
            }else{
                nams.kep <- c('per','lnK','sqresinw','sqrecosw','Tc')
            }
        }
        for(j in 1:Np){
            nams <- c(nams,paste(nams.kep,j,sep=''))
        }
    }
    for(i1 in 1:length(p)){
        if(i1==1 & par.global=='trend'){
            nams <- c(nams,paste0('a',i1,1:Npoly))
        }
        nams <- c(nams,paste0('b',i1))
        nams <- c(nams,paste0('s',i1))
        if(p[i1]>0){
            nams.ar <- 'phi'
            for(j in 1:p[i1]){
                nams <- c(nams,paste0(nams.ar,i1,j))
            }
            nams <- c(nams,paste0('alpha',i1))
        }
        if(q[i1]>0){
                nams.ma <- 'w'
                for(j in 1:q[i1]){
                    nams <- c(nams,paste0(nams.ma,i1,j))
                }
                nams <- c(nams,paste0('beta',i1))
        }
        if(ns[i1]>0){
            nams <- c(nams,paste0('c',i1,1:ns[i1]))
        }
    }
    if(FALSE){
#    if(TRUE){
        cat('length(nams)=',length(nams),'\n')
        cat('nams=',nams,'\n')
        cat('length(par.vector)=',length(par.vector),'\n')
    }
    names(par.vector) <- nams
    return(par.vector)
}
plot.labels.simple <- function(par.name){
    pats <- c('per','K','e','omega','Mo','a','b','beta','alpha','w','m','tau','s','c','lnK','sqrecosw','sqresinw','Tc')
    new.name <- par.name
    for(pat in pats){
        ind <- grep(paste0('^',pat),par.name)
        if(length(ind)>0) new.name[ind] <- paste0('bquote(',pat,'[',gsub(paste0('^',pat),'',par.name[ind]),'])')
    }
    sapply(new.name, function(x) eval(parse(text=x)))
#    nams <- vector('expression',length(new.name))
#    for(k in 1:length(nams)){
#        labs[k] <- substitute(lab.name,list(lab.name=nams[k]))
#    }
#    labs
}

plot.labels <- function(Nset,Np){
    nams <- c()
    if(Np>0){
        nam.kep <- c()
        if(prior.type!='e0'){
            nam.kep <- c(nam.kep,c('period[day]','period[day]',expression(nu*"[day"^{-1}*']'),'K[m/s]','e',expression(omega*'[rad]'),expression(M[0]*"[rad]")))
        }else{
            nam.kep <- c(nam.kep,c('period[day]','period[day]',expression(nu*"[day"^{-1}*']'),'K[m/s]',expression(omega*'[rad]')))
        }
        nams <- c(nams,rep(nam.kep,Np))
    }
    for(i3 in 1:Nset){
        if(i3==1 | par.global!='trend'){
            if(Npoly>0){
                if(any(grepl('^a',ep.par))){
                    nams <- c(nams,paste0(j1,'a',1:Npoly,'_',i3,'[m/s/year]'))
                }
            }
        }
        if(any(grepl('^b',ep.par))){
            nams <- c(nams,paste0('b',1:Nb,'_',i3,'[m/s]'))
        }
	    if(any(grepl('^s',ep.par))){
	        nams <- c(nams,paste0(paste0('s',1:Njitter),i3,'[m/s]'))
	    }
            if((noise.model=='TJ' | noise.model=='ARMATJ' | noise.model=='TJAR' | noise.model=='ARMATJAR') & !all(Inds==0)){
                nams <- c(nams,'ss[m/s/[Sindex]]')
                if(ntj==3){
                    nams <- c(nams,c('sa[m/s/[BIS]]','sb[m/s/[FWHM]]'))
                }
            }
            if((noise.model=='GP' | noise.model=='GPR') & !(gp.single & i3>1)){
                if(gp.type=='abs' | gp.type=='sq'){
                    nams <- c(nams,c(expression(sigma[red]*'[m/s]'),'l[day]'))
                }
                if(gp.type=='qp'){
                    nams <- c(nams,c(expression(sigma[red]*'[m/s]'),'l[day]',expression(P[gp]*'[day]'),expression(l[gp])))
                }
                if(gp.type=='sho'){
                    nams <- c(nams,c(expression(sigma[red]*'[m/s]'),expression('ln('*tau[l]*')[day]'),expression('ln('*P[rot]*')[day]')))
                }
            }
            if(noise.model=='ARMA' | noise.model=='ARMATJ' | noise.model=='ARMATJAR' | noise.model=='ARMAPSID'){
                                        #AR(p)
                if(p>0){
                    for(i2 in 1:p){
                        if(any(grepl('ar',ep.par))){
                            nams <- c(nams,paste0(expression(phi),i2))
                        }
                    }
                    nams <- c(nams,c(paste0(expression("log("*eta*"[day])"),i3)))
                }
                                        #MA(q)
                if(q>0){
                    if(any(grepl('ma',ep.par))){
                        for(i2 in 1:q){
                            nams <- c(nams,c(paste0(expression(w),i2)))
                        }
                        nams <- c(nams,c(expression("log("*tau*"[day])")))
                    }
                }
###indices
            if(!all(Inds==0)){
                cname <- proxy.names
                if(length(nepoch)==1){
                    nams <- c(nams,paste0(paste0('c',cname[1:length(Inds)],'_'),i3,'[m/s]'))
                }else if(any(grepl('^c',ep.par))){
                    int <- as.integer(gsub('c','',ep.par[grep('^c',ep.par)]))
                    nams <- c(nams,paste0(paste0('c',cname[int]),i3,'[m/s]'))
                }
            }

        }
    }
    labs <- vector('expression',length(nams))
    for(k in 1:length(nams)){
        labs[k] <- substitute(lab.name,list(lab.name=nams[k]))
    }
     return(labs)
}
cal.trend <- function(a,b,t){
    trend <- rep(b,length(t))
#cat('a=',a,'\n')
#cat('b=',b,'\n')
#cat('head(t)=',head(t),'\n')
    for(j3 in 1:length(a)){
        trend <- trend+a[j3]*t^j3#t[,j3]
    }
    return(trend)
}
####functions for mcmc algorithms for detection of keplerian signal in RV data
####Keplerian model with 1 planet
RV.kepler <- function(pars.kep,tt=NA,prior.kep=prior.type,period.kep=period.par,injection=FALSE,kep.only=FALSE,noise.only=FALSE){
    if(all(is.na(tt))){
        sim.kep <- FALSE
        Nrep <- length(ins)
    }else{
        sim.kep <- TRUE
        Nrep <- 1
    }
    if(!sim.kep){
        tt <- trv.all
    }
    Np.kep <- length(grep('per',names(pars.kep)))
    tt <- tt-tmin
    if(Np.kep>0 & !noise.only){
        Ps = exp(pars.kep[grep('^per([[:digit:]]{1})',names(pars.kep))])
        ##        pers[ind.transit] <- log(pers[ind.transit])
        if(basis=='natural'){
            Ks = pars.kep[grep('^K([[:digit:]]{1})',names(pars.kep))]
        }else if(basis=='linear1'){
            Ks = exp(pars.kep[grep('^lnK([[:digit:]]{1})',names(pars.kep))])
        }
        if(prior.kep!='e0'){
            es = pars.kep[grep('^e([[:digit:]]{1})',names(pars.kep))]
            Mos = pars.kep[grep('^Mo([[:digit:]]{1})',names(pars.kep))]
        }
        omegas = pars.kep[grep('^omega([[:digit:]]{1})',names(pars.kep))]
        sqresin = pars.kep[grep('^sqresinw([[:digit:]]{1})',names(pars.kep))]
        sqrecos = pars.kep[grep('^sqrecosw([[:digit:]]{1})',names(pars.kep))]
        Tc = pars.kep[grep('^Tc([[:digit:]]{1})',names(pars.kep))]

        dVr.p <- rep(0,length(tt))
        for(h in 1:Np.kep){
            K <- Ks[h]
            P <- Ps[h]
            if(basis=='natural'){
                omega <- omegas[h]
                if(prior.kep!='e0'){
                    e <- es[h]
                    Mo <- Mos[h]
                    m <- (Mo+2*pi*tt/P)%%(2*pi)#mean anomaly
                    E <- kep.mt2(m,e)
                    T <- 2*atan(sqrt((1+e)/(1-e))*tan(E/2))
                    dVr.p <- dVr.p+K*(cos(omega+T)+e*cos(omega))
                }else{
                    Mo <- Mos[h]
                    m <- (Mo+2*pi*tt/P)%%(2*pi)#mean anomaly
                    dVr.p <- dVr.p+K*cos(m)
                }
            }else if(basis=='linear1'){
                if(prior.kep!='e0'){
                    e <- sqresin[h]^2+sqrecos[h]^2
                    if(sqrecos[h]!=0){
                        omega <- atan(sqresin[h]/sqrecos[h])
                    }else{
                        omega <- atan(sqresin[h]/1e-3)
                    }
                    Mo <- getM0(e=e,omega=omega,P=P,T=Tc[h],T0=tmin,type='primary')
                    m <- (Mo+2*pi*tt/P)%%(2*pi)#mean anomaly
                    E <- kep.mt2(m,e)
                    T <- 2*atan(sqrt((1+e)/(1-e))*tan(E/2))#true anomaly
                    dVr.p <- dVr.p+K*(cos(omega+T)+e*cos(omega))
                }else{
                    Mo <- (tmin-Tc)*2*pi/P
                    m <- (Mo+2*pi*tt/P)%%(2*pi)#mean anomaly
                    dVr.p <- dVr.p+K*cos(m)
                }
            }
        }
    }else{
        dVr.p <- rep(0,length(tt))
    }
    ind.na <- which(is.na(tt))
    if(length(ind.na)>0){
        dVr.p[ind.na] <- NA
    }

###calculate noise component
    dVr.kep <- list()
    for(j2 in 1:Nrep){
###assign parameters and data
        nqp <- out[[ins[j2]]]$noise$nqp
        if(!sim.kep){
            dVr.kep[[ins[j2]]] <- dVr.p[out[[ins[j2]]]$index]
            tt <- out[[ins[j2]]]$RV[,1]-tmin
        }else{
            dVr.kep <- dVr.p
        }
        if(!kep.only){
            ind <- grep('^a1',names(pars.kep))
            if(length(ind)>0){
                at <- pars.kep[ind]
            }else{
                at <- 0
            }
#            cat('at=',at,'\n')
#            cat('b=',pars.kep[paste0('b',j2)],'\n')
            trend <- cal.trend(a=at,b=pars.kep[paste0('b',j2)],t=tt/time.unit)#
                                        #        cat('range(trend)=',range(trend),'\n')
#            cat('range(trend)=',range(trend),'\n')
            dVr.kep[[ins[j2]]] <- dVr.kep[[ins[j2]]]+trend
            if(nqp[1]>0){
                ind <- grep(paste0('c'),1:nqp[1])
                if(length(ind)>0){
                    cs <- par.kep[ind]
                    proxies <-out[[j2]]$noise$proxy.opt
                    if(length(cs)>1){
                        tmp <- rowSums(t(cs*t(proxies)))
                    }else{
                        tmp <- cs*proxies
                    }
                    dVr.kep[[ins[j2]]] <- dVr.kep[[ins[j2]]] + tmp
                }
            }
        }
    }
    return(dVr.kep)
}
####red noise kernel
ker.red <- function(x,y,sigma.red,l,Pgp=1,lp=1,type='abs'){
    if(type=='abs'){
        return(sigma.red^2*exp(-abs(x-y)/l))
    }
    if(type=='sq'){
        return(sigma.red^2*exp(-(x-y)^2/(2*l^2)))
    }
    if(type=='qp'){
#        return(sigma.red^2*exp(-(x-y)^2/(2*l^2)-(sin(pi*(x-y)/Pgp))^2/(2*lp^2)))
        val <- sigma.red^2*exp(-(x-y)^2/(2*l^2)-2*(sin(pi*(x-y)/Pgp))^2/(lp^2))#Haywood
        return(val)
    }
}
ker.sho <- function(x,y,S0,logQ,logProt){
    tau <- abs(x-y)
    Q <- exp(logQ)
    w0 <- 2*pi/exp(logProt)
    eta <- abs(1-(4*Q^2)^(-1))^0.5
    A <- S0*w0*Q*exp(-w0*tau/(2*Q))
    if(Q>0 & Q<1/2){
        A*(cosh(eta*w0*tau)+sinh(eta*w0*tau)/(2*eta*Q))
    }else if(Q==1/2){
        A*(2*(1+w0*tau))
    }else if(Q>1/2){
        A*(cos(eta*w0*tau)+sin(eta*w0*tau)/(2*eta*Q))
    }
}
ker.celerite <- function(x,y,term){
    a <- term[,1]
    b <- term[,2]
    c <- term[,3]
    d <- term[,4]
    R <- length(a)
    J <- R/2
    N <- length(x)
    tau <- abs(x-y)
    a*exp(-c*tau)*cos(d*tau)+b*exp(-c*tau)*sin(d*tau)
}
sho.term <- function(S0,logQ,logProt){
#    w0 <- 2*pi/Prot
    Q <- exp(logQ)
    w0 <- 2*pi/exp(logProt)
    a <- b <- c <- d<- 0
    if(Q<0.5){
        f <- sqrt(1-4*Q^2)
        a <- 0.5*S0*w0*Q*c(1+1/f,1-1/f)
        c <- 0.5*w0/Q*c(1-f,1+f)
    }else{
        f <- sqrt(4.0*Q^2-1)
        a <- S0*w0*Q
        b <- S0*w0*Q/f
        c <- 0.5*w0/Q
        d <- 0.5*w0*f/Q
    }
    return(cbind(a,b,c,d))
}
##refer to cholesky.h in celerite python module
celerite <- function(t,dy,y,s,term){
    a <- term[,1]
    b <- term[,2]
    c <- term[,3]
    d <- term[,4]
###a, b, c, d are for complex terms
###ar, cr are for real terms
#    Jc <- length(a)
#    Jr <- length(ar)
#    J <- Jc+Jr
#    R <- 2*J-Jr
    R <- 2*length(a)
    J <- length(a)
    N <- length(t)
    phi <- Wp <- Up <- Vp <- array(NA,dim=c(N,R))
    A <- D <- array(0,dim=c(N,N))
    S <- array(NA,dim=c(N,R,R))
#    diag(D) <- s^2+sum(a)+sum(ar)
    js <- 1:J
    ns <- 1:N
    dt <- outer(d,t,'*')
    ct <- outer(c,t,'*')
    cdt <- cos(dt)
    sdt <- sin(dt)
    nct <- exp(-ct)
    pct <- exp(ct)
    diag(A) <- dy^2+s^2+sum(a)
####pre-conditioned variables
    Up[,2*js-1] <- a%*%cdt+b%*%sdt
    Up[,2*js] <- a%*%sdt-b%*%cdt
    Wp[,2*js-1] <- cdt
    Wp[,2*js] <- sdt
    ect <- outer(c,t[2:N]-t[1:(N-1)],'*')
    phi[2:N,2*js-1] <- phi[2:N,2*js] <- exp(-ect)
    phi[1,] <- 0
####calculate S, D and W
    S[1,,] <- 0
    D[1,1] <- A[1,1]
    Wp[1,] <- 1/D[1,1]*Wp[1,]
    for(n in 2:N){
        S[n,,] <- outer(phi[n-1,],phi[n-1,],'*')*(S[n-1,,]+D[n-1,n-1]*outer(Wp[n-1,],Wp[n-1,],'*'))
        D[n,n] <- A[n,n]-Up[n,]%*%(S[n,,]%*%Up[n,])
        Wp[n,] <- 1/D[n,n]*(Wp[n,]-Up[n,]%*%S[n,,])
    }
    lndetK <- sum(log(diag(D)))
    f <- array(NA,dim=c(N,R))
    z <- rep(NA,N)
    z[1] <- y[1]
    f[1,] <- 0
    yky <- y[1]^2/D[1,1]
    for(n in 2:N){
        f[n,] <- phi[n,]*(f[n-1,]+Wp[n-1,]*z[n-1])
        z[n] <- y[n]-sum(Up[n,]*f[n,])
        yky <- yky+z[n]^2/D[n,n]
    }
    return(list(lndetK=lndetK,yky=yky))
}
inv.cel <- function(t,dy,s,a,b,c,d,ar,cr,y){
    R <- length(a)
    N <- length(t)
    tmp <- mat.cel(t,dy,s,a,b,c,d,ar,cr)
    W <- tmp$W
    U <- tmp$U
    D <- tmp$D
    phi <- tmp$phi
    f <- array(NA,dim=c(N,R))
    z <- rep(NA,N)
    z[1] <- y[1]
    f[1,] <- 0
    out <- y[1]^2/D[1,1]
    for(n in 2:N){
        f[n,2:R] <- phi[n,2:R]*(f[n-1,]+W[n-1,]*z[n-1])
        z[n] <- y[n]-sum(U[n,]*f[n,])
        out <- out+z[n]^2/D[n,n]
    }
    return(out)
}
ker.rot <- function(B, C, L, tau, Prot){
    B/(2+C)*exp(-tau/L)*(cos(2*pi*tau/Prot)+(1+C))
}
sho.chol <- function(tau,S0,Q,w0){
    K <- ker.sho(tau,S0,Q,w0)
    Sn
}
###red noise model
rednoise.cov <- function(sigma.red,l,tt=trv,tol=1e-12,Pgp=1,lp=1,type='abs'){
#method 1: easy but not so many choices of kernels
#    ker <- laplacedot(1/l)
#    cov.red <- kernelMatrix(ker,tt)#+diag(eRV^2)#to guarantee that it is positive definite
#method 2: more time-consuming
#    cov.red <- outer(1:length(RV),1:length(RV),Vectorize(function(u,v) ker.red(sigma.red,l,tt[u],tt[v],Pgp=Pgp,lp=lp,type=gp.type)))
#method 3:
    SE <- function(x,y) ker.red(x,y,sigma.red=sigma.red,l=l,Pgp=Pgp,lp=lp,type=gp.type)
    cov.red <- outer(tt, tt, SE)
    return(cov.red)
}
##arma noise model
arma <- function(t,ymodel,ydata,pars,p,q,ind.set=1,ARMAtype='abs'){
    dv.arma <- 0
    ##AR(p.ar) model
    if(p>0){
        dVr.ar <- 0
        phi <- pars[paste0('phi',ind.set,1:p)]
        alpha <- pars[paste0('alpha',ind.set)]
        for(j in 1:p){
            dVr.ar <- dVr.ar + c(rep(0,j),phi[j]*exp((t[-(length(t)+1-(1:j))]-t[-(1:j)])/exp(alpha))*ydata[-(length(t)+1-(1:j))])
        }
        dv.arma <- dv.arma + dVr.ar
    }
    ##MA(q) model
    if(q>0){
        dVr.ma <- 0
        w <- pars[paste0('w',ind.set,1:q)]
        beta <- pars[paste0('beta',ind.set)]
        res <- ydata-ymodel
        for(j in 1:q){
            dVr.ma <- dVr.ma + c(rep(0,j),w[j]*exp(-abs(t[-(length(t)+1-(1:j))]-t[-(1:j)])/exp(beta))*res[-(length(t)+1-(1:j))])
        }
        dv.arma <- dv.arma + dVr.ma
    }
    return(dv.arma)
}
# ----- Define a function for plotting a matrix ----- #
myImagePlot <- function(x, ...){
     min <- min(x)
     max <- max(x)
     yLabels <- rownames(x)
     xLabels <- colnames(x)
     title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
       min <- Lst$zlim[1]
       max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
       yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
       xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
       title <- Lst$title
    }
  }
# check for null values
if( is.null(xLabels) ){
   xLabels <- c(1:ncol(x))
}
if( is.null(yLabels) ){
   yLabels <- c(1:nrow(x))
}

layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))

 # Red and green range from 0 to 1 while Blue ranges from 1 to 0
 ColorRamp <- rgb( seq(0,1,length=256),  # Red
                   seq(0,1,length=256),  # Green
                   seq(1,0,length=256))  # Blue
 ColorLevels <- seq(min, max, length=length(ColorRamp))

 # Reverse Y axis
 reverse <- nrow(x) : 1
 yLabels <- yLabels[reverse]
 x <- x[reverse,]

 # Data Map
 par(mar = c(3,5,2.5,2))
 image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
 ylab="", axes=FALSE, zlim=c(min,max))
 if( !is.null(title) ){
    title(main=title)
 }
axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
 axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
 cex.axis=0.7)

 # Color Scale
 par(mar = c(3,2.5,2.5,2))
 image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=ColorRamp,
      xlab="",ylab="",
      xaxt="n")

 layout(1)
}
matrix.power <- function(mat,power){
    for(j1 in 1:ncol(mat)){
        ind <- which(mat[,j1]>1e-3)
        mat[ind,j1] <- mat[ind,j1]^power
    }
    return(mat)
}
#likelihood
loglikelihood <- function(pars){
    RV.kep  <-  try(RV.kepler(pars.kep=pars),TRUE)
#    if(class(RV.kep)=='try-error') RV.kep <- RV.kepler(pars.kep=startvalue)
    logLike <- 0
    for(k3 in 1:length(ins)){
        trv <-out[[ins[k3]]]$RV[,1]
        rv.data <- out[[ins[k3]]]$RV[,2]
        erv <- out[[ins[k3]]]$RV[,3]
        rv.model <- rv.kep <- RV.kep[[ins[k3]]]
        nqp <- out[[ins[k3]]]$noise$nqp
        s <- pars[paste0('s',k3)]
#        cat('nqp=',nqp,'\n')
#        cat('range(rv.data)=',range(rv.data),'\n')
#        cat('range(rv.kep)=',range(rv.kep),'\n')
        if(nqp[2]>0 | nqp[3]>0){
            rv.arma <- arma(t=trv,ymodel=rv.kep,ydata=rv.data,pars=pars,ind.set=k3,p=nqp[3],q=nqp[2])
#            cat('range(rv.arma)=',range(rv.arma),'\n')
            rv.model <- rv.model +rv.arma
        }
        ll <- sum(dnorm(rv.data,mean=rv.model,sd=sqrt(erv^2+s^2),log=T))
#        cat('k3=',k3,'\n')
#        cat('ll=',ll,'\n')
        logLike <-  logLike +ll
    }
    return(logLike)
}

cal.residual <- function(pars){
    RV.all = RV.kepler(pars.kep=pars)
    RV.kep = RV.kepler(pars.kep=pars,kep.only=TRUE)
    trend <- arma <- all <- noise <- signal <- residual.all <- residual.noise <- residual.sig <- list()
    for(k3 in 1:length(ins)){
        trv <-out[[ins[k3]]]$RV[,1]
        rv.data <- out[[ins[k3]]]$RV[,2]
        erv <- out[[ins[k3]]]$RV[,3]
        rv.all <- RV.all[[ins[k3]]]
        rv.sig <-RV.kep[[ins[k3]]]
        rv.trend <- rv.all-rv.sig
        nqp <- out[[ins[k3]]]$noise$nqp
        if((nqp[2]>0 | nqp[3]>0)){
            rv.arma <- arma(t=trv,ymodel=rv.all,ydata=rv.data,pars=pars,ind.set=k3,p=nqp[3],q=nqp[2])
#            cat('range(rv.arma)=',range(rv.arma),'\n')
            rv.all <- rv.all +rv.arma
        }else{
            rv.arma <-rep(0,length(trv))
        }
        res.all <- rv.data-rv.all
        res.sig <- rv.data-rv.sig
        rv.noise <- rv.all-rv.sig
        res.noise <- rv.data-rv.noise
        residual.all[[ins[k3]]] <- res.all
        residual.sig[[ins[k3]]] <- res.sig
        residual.noise[[ins[k3]]] <- res.noise
        signal[[ins[k3]]] <- rv.sig
        noise[[ins[k3]]] <- rv.noise
        all[[ins[k3]]] <- rv.all
        arma[[ins[k3]]] <- rv.arma
        trend[[ins[k3]]] <- rv.trend
    }
    list(res.sig=residual.sig,res.all=residual.all,res.noise=residual.noise,rv.signal=signal,rv.noise=noise,rv.all=all,rv.arma=arma,rv.trend=trend)
}

tjar <- function(t=trv,x=Sindex[,1],phi.sab,alpha.sab,psi.sab,par.tj=1,symmetry=TJAR.sym){
    rv.tjar <- 0
    for(i0 in 1:length(phi.sab)){
        if(symmetry=='sym'){
            rv.tjar <- rv.tjar + par.tj*c(rep(0,i0),phi.sab[i0]*exp(alpha.sab/time.unit*(t[-(length(t)+1-(1:i0))]-t[-(1:i0)]))*x[-(length(t)+1-(1:i0))]) + par.tj*c(psi.sab[i0]*exp(-alpha.sab/time.unit*(t[-(1:i0)]-t[-(length(t)+1-(1:i0))]))*x[-(1:i0)],rep(0,i0))
        }else if(symmetry=='past'){
            rv.tjar <- rv.tjar + par.tj*c(rep(0,i0),phi.sab[i0]*exp(alpha.sab/time.unit*(t[-(length(t)+1-(1:i0))]-t[-(1:i0)]))*x[-(length(t)+1-(1:i0))])
        }else if(symmetry=='future'){
            rv.tjar <- rv.tjar + par.tj*c(phi.sab[i0]*exp(-alpha.sab/time.unit*(t[-(1:i0)]-t[-(length(t)+1-(1:i0))]))*x[-(1:i0)],rep(0,i0))
        }
    }
    return(rv.tjar)
}
prior.func <- function(pars){
    logprior <- 0
    Np <- length(grep('per',names(pars)))
    if(Np>0){
        pers  <-  pars[grep('per([[:digit:]]{1})',names(pars))]
        if(period.par=='logP'){
            Ps <- exp(pers)
        }else if(period.par=='P'){
            Ps <- pers
        }else if(period.par=='nu'){
            Ps <- 1/pers
        }
        for(j4 in 1:Np){
            ii <- which(ind.transit==j4)
            if(length(ii)>0){
                per.logprior <- dnorm(Ps[j4],mean=Ptransit[ii],sd=ePtransit[ii],log=TRUE)
                if(per.logprior< -100 & Niter0<1e5) cat('per.logprior=',per.logprior,'\n')
                logprior <- logprior+per.logprior
            }
        }
        if(basis=='natural'){
            Ks <-  pars[grep('K([[:digit:]]{1})',names(pars))]
            omegas <-  pars[grep('omega([[:digit:]]{1})',names(pars))]
            Omegas <-  pars[grep('Omega([[:digit:]]{1})',names(pars))]
            Mos <-  pars[grep('Mo([[:digit:]]{1})',names(pars))]
            if(prior.type!='e0'){
                es = pars[grep('e([[:digit:]]{1})',names(pars))]
            }else{
                es <-rep(0,Np)
            }
            for(k in 1:Np){
                if(prior.type=='jeffrey'){
                    Kprior = (Ks[k]+K0)^-1*(log(1+Kmax/K0*(Pmin/Ps[k])^(1/3)))^-1#according to Table 1 of Ford & Gregory 2007
                }else{
                    Kprior = 1/(Kmax-Kmin)
                }
                logprior <- logprior+log(Kprior)
                if(prior.type!='e0' & Esd!=1){
                    eprior = 2*dnorm(es[k],mean=0,sd=Esd)#normalized semi-Gaussian distribution
                }else{
                    eprior = 1
                }
                logprior <- logprior+log(eprior)
                logprior <- logprior+2*log(1/(2*pi))#omega and Omega
                Mo.prior <- log(1/(2*pi))
####if type=natural, then there is no prior on Tc
                if(Ntransit>0 & !is.na(Tc)){
                    ii <- which(k==ind.transit)
                    if(length(ii)>0 & FALSE){
                        Mo0 <- getM0(e=es[k],omega=omegas[k],P=Ps[k],T=Tc[ii],T0=tmin,type='primary')
                        Mo1 <- getM0(e=es[k],omega=omegas[k],P=Ps[k],T=Tc[ii]-eTc[ii],T0=tmin,type='primary')
                        Mo2 <- getM0(e=es[k],omega=omegas[k],P=Ps[k],T=Tc[ii]+eTc[ii],T0=tmin,type='primary')
                        dMo <- abs((Mo2-Mo1)/2)
                        cat('Mo0',Mo0,'\n')
                        cat('Mos[k]',Mos[k],'\n')
                        cat('dMo',dMo,'\n')
                        Mo.prior <- dnorm(Mos[k],mean=Mo0,sd=dMo,log=TRUE)
                        cat('Mo.prior=',Mo.prior,'\n')
                    }
                }
                logprior <- logprior+Mo.prior
            }
        }else if(basis=='linear1'){
            sqresinw <- pars[grep('sqresinw([[:digit:]]{1})',names(pars))]
            sqrecosw <- pars[grep('sqrecosw([[:digit:]]{1})',names(pars))]
            Tcs <-  pars[grep('Tc([[:digit:]]{1})',names(pars))]
###assume uniform prior for eccentricity if basis=linear1
            if(prior.type!='e0' & Esd!=1 & FALSE){
                eprior = 2*dnorm(sqresinw^2+sqrecosw^2,mean=0,sd=Esd)#normalized semi-Gaussian distribution
            }else{
                eprior = 1
            }
            logprior <- logprior+sum(log(eprior))
            if(!is.na(Tc) & Ntransit>0){
                ttc <- Tcs[ind.transit[ind.transit!=0]]
                logprior.Tc <- dnorm(ttc,Tc,eTc,log=TRUE)
                if(logprior.Tc< -100 & Niter0<1e5)  cat('logprior.Tc=',logprior.Tc,'\n')
                logprior <- logprior+sum(logprior.Tc)
            }
        }
    }
####noise model
    for(i1 in 1:length(ins)){
        nqp <- out[[ins[i1]]]$noise$nqp
        n <- nqp[1];q <- nqp[2];p <- nqp[3]
        phiprior = wprior = 1/(phi.max-phi.min)
        if(p>0){
            alpha.prior <- 1/(alpha.max-alpha.min)
            logprior <- logprior +p*log(phiprior)+log(alpha.prior)
        }
        if(q>0){
            beta.prior <- 1/(beta.max-beta.min)
            logprior <- logprior +q*log(wprior)+log(beta.prior)
        }
    }
    return(logprior)
}

#posterior distribution
posterior <- function(param,tem=1){
    llike <- loglikelihood(param)
    pr <- prior.func(param)
    post <- llike*tem+pr
    return(list(loglike=llike,logprior=pr,post=post))
}

proposalfunction.simple <- function(param,cov.adapt){
    Ntt <- 1e4
    for(k in 1:Ntt){
        cov.all <- cov.adapt
                                        #        param.new <- try(mvrnorm(n=1,mu=param,Sigma=cov.all,tol=tol1),TRUE)#this could cause some non-zero values for fix value, but it is too small to be accounted for.
        param.new <- try(mvrnorm(n=1,mu=param,Sigma=cov.all,tol=tol1),TRUE)
                                        #        cat('diag(cov.all)=',diag(cov.all),'\n')
                                        #        param.new <- mvrnorm(n=1,mu=param,Sigma=cov.all,tol=tol1)
        if(class(param.new)=='try-error'){
            param.new <- try(mvrnorm(n=1,mu=param,Sigma=nearPD(cov.all)$mat,tol=tol1),TRUE)
        }
        ind1 <- grep('Mo\\d',names(param.new))
        ind2 <- grep('omega\\d',names(param.new))
        ind3 <- grep('Tc\\d',names(param.new))
        if(length(ind1)>0) param.new[ind1] <- param.new[ind1]%%(2*pi)
        if(length(ind2)>0) param.new[ind2] <- param.new[ind2]%%(2*pi)
####don't mod the transit epoch by period
        if(length(ind3)>0 & all(is.na(Tc))){
                                        #cat('param.new[ind3]=',param.new[ind3],'\n')
            if(prior.type=='e0'){
                param.new[ind3] <- tmin+(param.new[ind3]-tmin)%%param.new[ind3-2]
            }else{
                param.new[ind3] <- tmin+(param.new[ind3]-tmin)%%param.new[ind3-4]
            }
            ##cat('2param.new[ind3]=',param.new[ind3],'\n')
        }
        ind1 <- grep('sqresinw',names(param.new))
        ind2 <- grep('sqrecosw',names(param.new))
        if(length(ind1)>0){
            ind <- c(ind1,ind2)
            if(all(param.new[-ind]>par.min[-ind] & param.new[-ind]<par.max[-ind]) & all((param.new[ind1]^2+param.new[ind2]^2)<1)) break()
        }else{
            if(all(param.new>par.min & param.new<par.max)) break()
        }
    }
    if(k==Ntt){
        cat('k=',k,'\n')
        cat('par.max-par.min=',par.max-par.min,'\n')
        cat('param.new[1]=',param.new[1],'\n')
        cat('startvalue[1]=',startvalue[1],'\n')
        cat('par.min[1]=',par.min[1],'\n')
        cat('par.max[1]=',par.max[1],'\n')
        cat('name(param>par.max)=',names(param)[which(param.new>par.max)],'\n')
        cat('name(param<par.min)=',names(param)[which(param.new<par.min)],'\n')
        break()
                                        #        cat('tough var:',names(param.new)[param.new<par.min | param.new>par.max],'\n')
                                        #        cat('diag(cov.all)=',diag(cov.all),'\n')
    }
    return(param.new)
}

#Metropolis algorithm####
#####Initial proposal function
proposalfunction <- function(param,cov.adapt){
    Ntt <- 1e4
    for(k in 1:Ntt){
        cov.all <- cov.adapt
        param.new <- try(mvrnorm(n=1,mu=param,Sigma=cov.all,tol=1e-10,empirical=FALSE))#this could cause some non-zero values for fix value, but it is too small to be accounted for.
     if(class(param.new)=='try-error'){
            param.new <- try(mvrnorm(n=1,mu=param,Sigma=nearPD(cov.all),tol=tol1),TRUE)
        }
        if(Np>0){
            M0.sim <- param.new[grep('Mo([[:digit:]]{1})',names(param.new))]
            param.new[grep('Mo([[:digit:]]{1})',names(param.new))] <- M0.sim%%(2*pi)
            omega.sim <- param.new[grep('omega([[:digit:]]{1})',names(param.new))]
            param.new[grep('omega([[:digit:]]{1})',names(param.new))] <- omega.sim%%(2*pi)
        }
        logplast <- param.new[grep(paste0('per',Np),names(param.new))]
        if(!quantify & Np>0){
            logplast <- param.new[grep(paste0('per',Np),names(param.new))]
            if(any(logplast<logPmaxs & logplast>logPmins)){
                logic.per <- TRUE
            }else{
                logic.per <- FALSE
            }
            if(all(param.new>par.min & param.new<par.max) & all(logic.per)) break()
        }else{
            if(all(param.new>par.min & param.new<par.max)) break()
        }
    }
#    cat('k=',k,'\n')
    if(k==Ntt & FALSE){
	cat('param.new=',param.new,'\n')
	cat('par.min=',par.min,'\n')
	cat('par.max=',par.max,'\n')
        cat('The times of generating the proposed parameters reach the maximum value!\n')
    }
    names(param.new) <- names(startvalue)
    return(param.new)
}

#####decide whether to accept or adjust the new period values
newper <- function(pars0,pars,pers.low,pers.up){
    ind.last <- which(names(pars)==paste0('per',Np))
    if(Np>1){
        inds = grep('per([[:digit:]]{1})',names(pars))
        ind.former <- inds[inds!=ind.last]
        per.former <- pars[ind.former]
        per0.former <- pars0[ind.former]
    }
    per.last <- pars[ind.last]
    tmp <- par.sec(per.last,pers.low,pers.up)
    logic.all <- logic1 <- tmp$logic
    pars[ind.last] <- tmp$par.new
###former period parameters
    if(length(pers.up)>1 & !fixP){
        pf.low <- pers.up[-length(pers.up)]
        pf.up <- pers.low[-1]
        if(Np>1){
            logic2 <- c()
            for(k1 in 1:length(per.former)){
                ind <- which(per0.former[k1]>pf.low & per0.former[k1]<pf.up)
                if(per.former[k1]>pf.low & per.former[k1]<pf.up){
                    logic2 <- c(logic2,TRUE)
                }else{
                    logic2 <- c(logic2,FALSE)
                }
            }
            logic.all <- all(logic1,logic2)
        }
    }
    return(list(pars=pars,logic=logic.all))
}
####
par.sec <- function(par,pars.low,pars.up){
    if(all(par>pars.low & par<pars.up)){
        log1 <- TRUE
        par.new <- par
    }else if(par< max(pars.up) & par>min(pars.low)){
        dpar.low <- min(abs(par-pars.low))
        dpar.up <- min(abs(par-pars.up))
        if(dpar.low<dpar.up){
            ind.min <- which.min(abs(par-pars.low))
            par.new <- min(pars.low[ind.min]+dpar.low,pars.up[ind.min])
        }else{
            ind.min <- which.min(abs(par-pars.up))
            par.new <- max(pars.up[ind.min]-dpar.up,pars.low[ind.min])
        }
        log1 <- TRUE
    }else{
        log1 <- FALSE
        par.new <- par
    }
    return(list(par.new=par.new,logic=log1))
}
####adaptive proposal function
###calculate covariance matrix
covariance.n0 <- function(mat,eps=1e-6){
    cov.par <- Sd*cov(mat)+Sd*eps*diag(ncol(mat))
    return(cov.par)
}

covariance.rep <- function(parms,cov1,mu1,n,Nupdate=1,eps=1e-6){
    if(n%%Nupdate==0){
        Sd <- 2.4^2/Npar
        mu2 <- (mu1*(n-1)+parms)/n#mu1 and mu2 are the mean vectors of parameters for n-1 and n iterations
        N <- n-1
        cov2 <- (N-1)*cov1/N+Sd*(N*mu1%*%t(mu1)-(N+1)*mu2%*%t(mu2)+parms%*%t(parms)+eps*diag(length(parms)))/N
    }else{
        cov2 <- cov1
    }
    return(cov2)
}
###run MH algorithm
run.metropolis.MCMC <- function(startvalue,cov.start,iterations,n0=2,verbose=FALSE,tem=1){
    Npar <- length(startvalue)
    chain  <-  array(dim=c(iterations+1,Npar))
    logpost = loglike = rep(NA,iterations+1)
    colnames(chain) <- names(startvalue)
    chain[1,]<- startvalue
    mu1 <- chain[1,]
    logpost.out = posterior(chain[1,],tem=tem)
    logpost[1] <- logpost.pre <- logpost.out$post
    loglike[1] <- loglike.pre <- logpost.out$loglike
    logprior.pre <- logpost.out$logprior
    dt0 <- 0
    cov.adapt <- array(data=0,dim=c(Npar,Npar))
    t.start <- proc.time()
    for(i in 1:iterations){
        if(i == n0){
            cov.adapt <- covariance.n0(mat=chain[1:n0,])
        }else if(i > n0){
            cov.adapt <- covariance.rep(chain[i,],cov.adapt,mu1,i)
        }else{
            cov.adapt <- cov.start
        }
        proposal = proposalfunction.simple(chain[i,],cov.adapt)
        proprop <- posterior(proposal,tem=tem)
        logpost.prop <- proprop$post
	logprior.prop <- proprop$logprior
        loglike.prop <- proprop$loglike

        logpost.cur <- logpost.pre
        loglike.cur <- loglike.pre
	logprior.cur <- logprior.pre
        if(is.na(logpost.prop)){
            probab = 0
        }else{
            probab = exp(logpost.prop-logpost.cur)
        }
        if(runif(1)<probab){
            chain[i+1,]=proposal
###values for next proposal estimation
            logpost.pre <- logpost.prop
            loglike.pre <- loglike.prop
	    logprior.pre <- logprior.prop
        }else{
            chain[i+1,]=chain[i,]
            logpost.pre <- logpost.cur
            loglike.pre <- loglike.cur
            logprior.pre <- logprior.cur
        }
##save values
        logpost[i+1]<- logprior.pre + loglike.pre
        loglike[i+1] <- loglike.pre
        mu1 <- (mu1*i+chain[i+1,])/(i+1)
##moniter parameters
        if((i%%100==0 | i==1) & verbose){
#            cat('i=',i,';tem=',tem,';diag(cov.adapt)=',diag(cov.adapt),';')
            cat('i=',i,';tem=',tem,';')
            cat('acceptance percentage:',100*(1-mean(duplicated(chain[1:(i+1),1:Npar]))))
            if(Np>=1){
                pars.now <- chain[i+1,]
                pers = pars.now[grep('per([[:digit:]]{1})',names(pars.now))]
                Ks = pars.now[grep('K([[:digit:]]{1})',names(pars.now))]
                es = pars.now[grep('e([[:digit:]]{1})',names(pars.now))]
                if(period.par=='nu'){
                    P.cur <- 1/pers
                }else if(period.par=='P'){
                    P.cur <- pers
                }else if(period.par=='logP'){
                    P.cur <- exp(pers)
                }else{
                    cat('The period parameter is not recognized!\n')
                }
                cat('; period= ',P.cur,'; K=', Ks)
#                cat('; period= ',P.cur,'; Pmin=',Pmin,'; Pmax=',Pmax,'K=', Kmul)
                if(prior.type!='e0'){
                    cat('; e=',es,'\n')
                }
            }
            cat('; max(logpost)=',max(logpost[1:(i+1)]),'; maximum likelihood:',max(loglike[1:(i+1)]),'\n')
            t.start <- proc.time()
        }
        conv <- Rhat <- NULL
        if(i==iterations){
#            cat('i=',i,'\n')
###check convergence of chains
            Nsub <- 5
            mcmc <- chain[,!grepl('Mo|omega|logpost|loglike',colnames(chain))]
            chain.len <- nrow(mcmc)
            subchain.len <- floor(chain.len/Nsub)
            var.sub <- array(data=NA,dim=c(Nsub,ncol(mcmc)))#s
            mean.sub <- array(data=NA,dim=c(Nsub,ncol(mcmc)))#theta.b
            meanvar <- rep(NA,ncol(mcmc))#W
            mean.all <- rep(NA,ncol(mcmc))#theta.bb
            var.single.est <- rep(NA,ncol(mcmc))#B
            for(j1 in 1:ncol(mcmc)){
                for(i1 in 1:Nsub){
                    var.sub[i1,j1] <- var(mcmc[((i1-1)*subchain.len+1):(i1*subchain.len),j1])
                    mean.sub[i1,j1] <- mean(mcmc[((i1-1)*subchain.len+1):(i1*subchain.len),j1])
                }
                meanvar[j1] <- mean(var.sub[,j1])
                mean.all[j1] <- mean(mean.sub[,j1])
                var.single.est[j1] <- subchain.len/(Nsub-1)*sum((mean.sub[,j1]-mean.all[j1])^2)
            }
            var.est <- (subchain.len-1)/subchain.len*meanvar+1/subchain.len*var.single.est#Var.hat
            Rhat <- sqrt(var.est/meanvar)
            acceptance <- 100*(1-mean(duplicated(chain[1:i,])))
            if(any(is.na(Rhat))) Rhat[is.na(Rhat)] <- 1
            if(any(Rhat>1.1)){
                conv <- FALSE
            }else{
                conv <- TRUE
#                chain <- chain[1:i,]
#                if(i>1e5){
#                    break()
#                }
            }
        }
    }
    val <- list(out=cbind(chain,logpost,loglike),Rhat=Rhat,conv=conv,acc=acceptance)
#    cat('boundary of period=',exp(par.min[1]),exp(par.max[1]),'d\n')
#    cat('range of period=',exp(range(chain[,(Np-1)*Nkeppar+1])),'d\n')
    return(val)
}
bin.simple <- function(data,Nbin){
    x <- data[,1]
    y <- data[,2]
    dy <- data[,3]
    xmid <- seq(min(x),max(x),length.out=Nbin+1)
    x1 <- xmid[-1]
    x2 <- xmid[-length(xmid)]
#    xbin <- (x1+x2)/2
    xbin <- sapply(1:Nbin,function(j) sum(x[x<x1[j] & x>x2[j]]/dy[x<x1[j] & x>x2[j]]^2)/sum(1/dy[x<x1[j] & x>x2[j]]^2))
    ybin <- sapply(1:Nbin,function(j) sum(y[x<x1[j] & x>x2[j]]/dy[x<x1[j] & x>x2[j]]^2)/sum(1/dy[x<x1[j] & x>x2[j]]^2))
    dybin <- sapply(1:Nbin,function(j) 1/sqrt(sum(1/dy[x<x1[j] & x>x2[j]]^2)))
    tmp <- cbind(xbin,ybin,dybin)
    tmp[which(!is.na(tmp[,1])),,drop=FALSE]
}

binning.post <- function(par.val,post.val,like.val){
    par.min<- min(par.val)
    par.max <- max(par.val)
    bin <- (par.max-par.min)/Nbins
    post.bin <- c()
    like.bin <- c()
    par.bin <- c()
    ind.na <- c()
    for(i in 1:Nbins){
        ind <- which((par.val>=par.min+(i-1)*bin) & (par.val<par.min+i*bin))
        if(length(ind)==0){
            post.bin <- c(post.bin,min(post.val))
	    like.bin <- c(like.bin,min(like.val))
            ind.na <- c(ind.na,i)
        }else{
#            post.bin <- c(post.bin,max(post.val[ind]))
#            like.bin <- c(like.bin,max(like.val[ind]))
            post.bin <- c(post.bin,NA)
            like.bin <- c(like.bin,NA)
        }
        par.bin <- c(par.bin,par.min+(2*i-1)*bin/2)
    }
    return(list(likes=like.bin,posts=post.bin,pars=par.bin,ind.na=ind.na))
}
###a more efficient way to binning parameters and log posteriors/likelihoods
binning.post2 <- function(par.val,post.val,like.val,Nbins){
    ind <- sort(par.val,index.return=TRUE)$ix
    par.sort <- par.val[ind]
    post.sort <- post.val[ind]
    like.sort <- like.val[ind]
    p1 <- hist(par.sort,breaks=seq(min(par.sort),max(par.sort),length.out=Nbins+1),plot=FALSE)
    index <- c(0,cumsum(p1$counts))
    post.max <- rep(NA,length(p1$mids))
    like.max <- rep(NA,length(p1$mids))
    ind.na <- c()
    for(j in 1:(length(index)-1)){
        if(index[j+1]>index[j]){
            post.max[j] <- max(post.sort[(index[j]+1):index[j+1]])
            like.max[j] <- max(like.sort[(index[j]+1):index[j+1]])
        }
    }
    ind.na <- which(is.na(post.max))
    return(list(likes=like.max,posts=post.max,pars=p1$mids,ind.na=ind.na))
}
kepler.tm <- function(Kp,wp,ep,Ea){
    Kp*sqrt(1-ep^2)*(cos(wp)*sqrt(1-ep^2)*cos(Ea)-sin(wp)*sin(Ea)/(1-ep*cos(Ea)))
}

kepler.ford <- function(Kp,wp,ep,Ea){
    Tp <- 2*atan(sqrt((1+ep)/(1-ep))*tan(Ea/2))
    rv <- Kp*(cos(wp+Tp)+ep*cos(wp))
    return(rv)
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

tcol <- function(color, percent = 50, name = NULL) {
    rgb.val <- col2rgb(color)
    t.col <- rgb(rgb.val[1,], rgb.val[2,], rgb.val[3,],
                 max = 255,
                 alpha = (100-percent)*255/100,
                 names = name)
    invisible(t.col)
}

calTc <- function(Tp,P,e,omega){
    theta <- 0.5*pi-omega
    E <- 2*atan(sqrt((1-e)/(1+e))*tan(theta/2))
    M <- E-e*sin(E)
    Tc <- (M%%(2*pi))/(2*pi)*P+Tp
    Tc
}

M02Tp <- function(M0,T0,P){
    T0-(M0%%(2*pi))*P/(2*pi)
}

Tc2M0.circular <- function(Tc,T0,P,omega){
    Tp <- Tc-(((0.5*pi-omega)/(2*pi))%%1)*P
    M0 <- (((T0-Tp)/P)%%1)*2*pi
    M0
}

data.distr <- function(x,lp=NULL,xlab,ylab,main='',oneside=FALSE,plotf=TRUE){
    xs <- seq(min(x),max(x),length.out=1e3)
    fitnorm <- fitdistr(x,"normal")
    p <- hist(x,plot=FALSE)
    xfit <- length(x)*mean(diff(p$mids))*dnorm(xs,fitnorm$estimate[1],fitnorm$estimate[2])
    ylim <- range(xfit,p$counts)
    if(plotf){
        plot(p,xlab=xlab,ylab=ylab,main=main,ylim=ylim)
        lines(xs,xfit,col='red')
    }
    xopt <- NA
    if(!is.null(lp)){
        xopt <- x[which.max(lp)]
    }
    x1=Mode(x)
    x2=mean(x)
    x3=sd(x)
    x4=skewness(x)
    x5=kurtosis(x)
    xs = sort(x)
    x1per = max(min(xs),xs[floor(length(xs)*0.01)])
    x99per = min(xs[ceiling(length(xs)*0.99)],max(xs))
    x10per = max(min(xs),xs[floor(length(xs)*0.1)])
    x90per = min(xs[ceiling(length(xs)*0.9)],max(xs))
    xminus.1sig = max(min(xs),xs[floor(length(xs)*0.15865)])
    xplus.1sig = min(xs[ceiling(length(xs)*(1-0.15865))],max(xs))
#    abline(v=c(x1per,x99per),col='blue')
    if(plotf){
        if(!oneside){
            legend('topleft',legend=c(as.expression(bquote('mode ='~.(format(x1,digit=3)))),as.expression(bquote(mu~'='~.(format(x2,digit=3)))),as.expression(bquote(sigma~'='~.(format(x3,digit=3))))),bty='n')
#            legend('topright',legend=c(as.expression(bquote('MAP ='~.(format(xopt,digit=3)))),as.expression(bquote(mu^3~'='~.(format(x4,digit=3)))),as.expression(bquote(mu^4~'='~.(format(x5,digit=3))))),bty='n')
            legend('topright',legend=c(as.expression(bquote('MAP ='~.(format(xopt,digit=3)))),as.expression(bquote(q[10]~'='~.(format(x10per,digit=3)))),as.expression(bquote(q[90]~'='~.(format(x90per,digit=3))))),bty='n')
        }else{
            legend('topleft',legend=c(as.expression(bquote('mode ='~.(format(x1,digit=3)))),as.expression(bquote(mu~'='~.(format(x2,digit=3)))),as.expression(bquote(sigma~'='~.(format(x3,digit=3)))),as.expression(bquote(mu^3~'='~.(format(x4,digit=3)))),as.expression(bquote(mu^4~'='~.(format(x5,digit=3))))),bty='n')
        }
    }
    tmp <- c(xopt=xopt,x1per=x1per,x99per=x99per,x10per=x10per,x90per=x90per,xminus.1sig=xminus.1sig,xplus.1sig=xplus.1sig,mode=x1,mean=x2,sd=x3,skewness=x4,kurtosis=x5)
    return(tmp)
}
matrix.distr <- function(y,lp=NULL,xlab,ylab,main='',oneside=FALSE,plotf=TRUE){
    tmp <- c()
    for(kk in 1:ncol(y)){
        x <- y[,kk]
        xs <- seq(min(x),max(x),length.out=1e3)
        fitnorm <- fitdistr(x,"normal")
        p <- hist(x,plot=FALSE)
        xfit <- length(x)*mean(diff(p$mids))*dnorm(xs,fitnorm$estimate[1],fitnorm$estimate[2])
        ylim <- range(xfit,p$counts)
        if(plotf){
            plot(p,xlab=xlab,ylab=ylab,main=main,ylim=ylim)
            lines(xs,xfit,col='red')
        }
        xopt <- NA
        if(!is.null(lp)){
            xopt <- x[which.max(lp)]
        }
        x1=Mode(x)
        x2=mean(x)
        x3=sd(x)
        x4=skewness(x)
        x5=kurtosis(x)
        xs = sort(x)
        x1per = max(min(xs),xs[floor(length(xs)*0.01)])
        x99per = min(xs[ceiling(length(xs)*0.99)],max(xs))
        x10per = max(min(xs),xs[floor(length(xs)*0.1)])
        x90per = min(xs[ceiling(length(xs)*0.9)],max(xs))
        xminus.1sig = max(min(xs),xs[floor(length(xs)*0.15865)])
        xplus.1sig = min(xs[ceiling(length(xs)*(1-0.15865))],max(xs))
                                        #    abline(v=c(x1per,x99per),col='blue')
        if(plotf){
            if(!oneside){
                legend('topleft',legend=c(as.expression(bquote('mode ='~.(format(x1,digit=3)))),as.expression(bquote(mu~'='~.(format(x2,digit=3)))),as.expression(bquote(sigma~'='~.(format(x3,digit=3))))),bty='n')
                                        #            legend('topright',legend=c(as.expression(bquote('MAP ='~.(format(xopt,digit=3)))),as.expression(bquote(mu^3~'='~.(format(x4,digit=3)))),as.expression(bquote(mu^4~'='~.(format(x5,digit=3))))),bty='n')
                legend('topright',legend=c(as.expression(bquote('MAP ='~.(format(xopt,digit=3)))),as.expression(bquote(q[10]~'='~.(format(x10per,digit=3)))),as.expression(bquote(q[90]~'='~.(format(x90per,digit=3))))),bty='n')
            }else{
                legend('topleft',legend=c(as.expression(bquote('mode ='~.(format(x1,digit=3)))),as.expression(bquote(mu~'='~.(format(x2,digit=3)))),as.expression(bquote(sigma~'='~.(format(x3,digit=3)))),as.expression(bquote(mu^3~'='~.(format(x4,digit=3)))),as.expression(bquote(mu^4~'='~.(format(x5,digit=3))))),bty='n')
            }
        }
        tmp <- cbind(tmp,c(xopt=xopt,x1per=x1per,x99per=x99per,x10per=x10per,x90per=x90per,xminus.1sig=xminus.1sig,xplus.1sig=xplus.1sig,mode=x1,mean=x2,sd=x3,skewness=x4,kurtosis=x5))
    }
    colnames(tmp) <- colnames(y)
    return(tmp)
}
####function to select output files from MCMC chain
select.file <- function(sub.out){
   fout <- NA
   if(length(sub.out)!=0){
      ind1 <- which(as.logical(sub.out[,5])==TRUE & as.numeric(sub.out[,3])>10 & as.numeric(sub.out[,3])<35)
      if(length(ind1)!=0){
        out2 <- sub.out[ind1,]
        if(length(ind1)==1){
           tmp <- out2
        }else{
#          ind2 <- which.max(as.numeric(out2[,4]))
           ind2 <- which.min(as.numeric(out2[,3]))
           tmp <- out2[ind2,]
        }
        fout <- as.character(tmp[1])
      }
  }
  if(!is.na(fout)){
    cat('Choose this mcmc results for further investigations: ', fout,'\n')
  }else{
    cat('No qualified mcmc chain for further investigations!\n')
  }
   return(fout)
}

###a function to reset -Inf element in an evidences
inf.rm <- function(E,B,flag=NA){
    ind.m <- which(!is.na(E) & E!=Inf & E!=-Inf,arr.ind=T)
    Emin <- min(E[ind.m])
    ind.b <- which(!is.na(B) & B!=Inf & B!=-Inf,arr.ind=T)
    Bmin <- min(B[ind.b])
    Edim <- dim(E)
    cat('Emin=',Emin,'\n')
    for(i in 1:Edim[3]){
        for(j in 1:Edim[2]){
            if(any(E[,j,i]==-Inf) | any(E[,j,i]==Inf) | any(is.na(E[,j,i])) | B[j,i]==Inf | B[j,i]==-Inf | is.na(B[j,i])){
                if(flag!=NA){
                    E[,j,i] <- Emin
                    B[j,i] <- 0
                }else{
                    E[,j,i] <- NA
                    B[j,i] <- NA
                }
            }
        }
    }
    return(list(E=E,B=B))
}
stemPlot <- function(x,y,pch=16,linecol=1,clinecol=1,add=FALSE,pair=FALSE,...){
    if(!add){
        plot(x,y,pch=pch,...)
    }else{
        points(x,y,pch=pch)
    }
    if(!pair){
        for (i in 1:length(x)){
            lines(c(x[i],x[i]), c(0,y[i]),col=linecol)
        }
    }else{
        ind.pair <- which((1:as.integer(length(x)/2))%%2==1)
        ind.solid <- c(2*ind.pair-1,2*ind.pair)
        ind.dashed <- c(1:length(x))[-ind.solid]
        for (i in 1:length(x)){
            if(any(i==ind.solid)){
                lines(c(x[i],x[i]), c(0,y[i]),col=linecol)
            }else{
                lines(c(x[i],x[i]), c(0,y[i]),col=linecol,lty=2)
            }
        }
    }
 #   lines(c(x[1]-2,x[length(x)]+2), c(0,0),col=clinecol,lty=3)
    abline(h=0,lty=2)
}
###functions
tpm <- function(priors, likes, lambda, h){
    prod.pre <- c(rep(0,h),likes[-((length(likes)-h+1):length(likes))])
    Ptpm <- sum(priors*likes/((1-lambda)*likes*priors+lambda*prod.pre))/sum(priors/((1-lambda)*likes*priors+lambda*prod.pre))
    return(log(Ptpm))
}
calc.alpha <- function(param,cov.adapt,post0,order='01'){
#        par.prop <- rbind(par.prop,par.tmp)
        post.prop <- posterior(param,tem=tem)$post
        if(order=='01'){
            if(is.na(post.prop)){
                probab = 0
            }else{
                probab = exp(post.prop-post0)
            }
        }else{
            if(is.na(post.prop)){
                probab = 1
            }else{
                probab = exp(post0-post.prop)
            }
        }
        min(1,probab)#since q1=q0, probab*q0/q1=probab; ref Chib and Jeliazkov 2001
}
####parallel computing
AMH <- function(nburn,Ns){
    if(nburn==0) nburn <- 1
    out.amh <- foreach(n = 1:Ncores, .combine='comb',.multicombine=TRUE) %dopar% {
        tmp <- run.metropolis.MCMC(startvalue,cov.start,Ns)
        list(out=cbind(tmp$chain[-(1:nburn),],tmp$post[-(1:nburn)],tmp$like[-(1:nburn)]),cov=tmp$cov)
#        cbind(tmp$chain[-(1:nburn),],tmp$post[-(1:nburn)],tmp$like[-(1:nburn)])
    }
#    return(out.amh)
    return(list(out=out.amh$out,cov=out.amh$cov))
}
####parallel chain returning covariance
AMH2 <- function(nburn,Ns){
    if(nburn==0) nburn <- 1
    out.amh <- foreach(n = 1:Ncores, .combine='comb',.multicombine=TRUE) %dopar% {
        tmp <- run.metropolis.MCMC(startvalue,cov.start,Ns)
        list(out=cbind(tmp$chain[-(1:nburn),],tmp$post[-(1:nburn)],tmp$like[-(1:nburn)]),cov=tmp$cov)
    }
    return(list(out=out.amh$out,cov=out.amh$cov))
}
####parallel chain with different temperture
AMH3 <- function(nburn,Ns){
    soft <- TRUE
    if(!soft){
        zeta <- tem.min
        Nlow <- floor(Ncores/2)+1
        Nup <- ceiling(Ncores/2)
        eta <- log(tem/tem.min)/log(Nlow)
        tems <- tempering(1:Nlow,zeta,eta)
        zeta <- tem
        eta <- log(1/tem)/log(Nup)
        tems <- c(tems,sort(tempering(Nup:1,zeta,eta))[-1])
    }else{
        zeta <- tem.min
        eta <- log(tem/tem.min)/log(Ncores)
        tems <- tempering(1:Ncores,zeta,eta)
    }
    tems <- sort(tems,decreasing=TRUE)
    if(nburn==0) nburn <- 1
    out.amh <- foreach(n = 1:Ncores, .combine='comb',.multicombine=TRUE) %dopar% {
        tem <- tems[n]
#        if(n==1){
#        if(n<=Ncores/2){
#        if(n<=3*Ncores/4){
#        if(n<=length(ps) & n<=Ncores/2){
        par.tmp <- pars
#        par.tmp <- par.end
        if(n<=(length(par.tmp)/Npar) & n<=(Ncores/2)){
            if(length(par.tmp)==Npar){
                startvalue <- par.tmp
            }else{
                startvalue <- par.tmp[n,]
            }
#            indP <- (Np-1)*Nkeppar+1
#            if(cov.start[indP,indP]<1e-3){
#                cov.start[indP,indP] <- startvalue[indP]*1e-3
#            }

            tem <- tems0[Ntrace[n]]
#            tem <- 1
        }else{
#            tem <- tems[1]
#            tem <- tems0[1]
            tem <- 1
            source('prepare_par.R',local=TRUE)
        }
#        source('mcmc_func.R',local=TRUE)
        tmp <- run.metropolis.MCMC(startvalue,cov.start,Ns)
        list(out=cbind(tmp$chain[-(1:nburn),],tmp$post[-(1:nburn)],tmp$like[-(1:nburn)]),cov=tmp$cov)
    }
    return(list(out=out.amh$out,cov=out.amh$cov))
}
###tempering function; power function
tempering <- function(x,a,b){
    a*x^b
}
####MCMC stuck steps calculation
stuck <- function(arr,post,Nstuck.min){
#    index <- which(duplicated(arr[,((Np-1)*Nkeppar+1):Npar]))
    index <- which(duplicated(arr[,(Np-1)*Nkeppar+1]))
    Ns <- 1
    Nstuck <- c()
    ind <- c()
    ind.stuck <- c()
    for(k in 2:length(index)){
        if((index[k]-index[k-1])==1){
            Ns <- Ns+1
            ind <- c(ind,index[k])
        }else if(length(ind)>0){
            Nstuck <- c(Nstuck,Ns)
            ind.stuck <- c(ind.stuck,ind[length(ind)])
            Ns <- 1
            ind <- c()
        }else{
	    Ns <- 1
            ind <- c()
	}
    }
    index <- which(Nstuck>Nstuck.min)
    if(length(index)>0){
        nstuck <- Nstuck[index]
        indstuck <- ind.stuck[index]
        index1 <- which.max(post[indstuck])
        indstuck <- indstuck[index1]
        nstuck <- nstuck[index1]
    }else{
        index1 <- which.max(Nstuck)
        indstuck <- ind.stuck[index1]
        nstuck <- Nstuck[index1]
    }
    ind.max <- which.max(post)
    per.stuck <- arr[indstuck,(Np-1)*Nkeppar+1]
    per.max <- arr[ind.max,(Np-1)*Nkeppar+1]
    #shift the to the local maxima
    if(abs(per.stuck-per.max)<(0.01*per.stuck)){
        indstuck <- ind.max
    }
#else{
#        indstuck <- 1
#        nstuck <- 1
#    }
    return(list(Nstuck=nstuck,ind.stuck=indstuck))
}
# Print the hostname for each cluster member
GetClusterInfo <- function() {
  info <- Sys.info()[c("nodename", "machine")]
  paste("Node:", info[1], "with CPU type", info[2])
}

###new Plows and Pups for divided period space
period.division <- function(plows,pups,Pmax,Pmin){
    plows.new <- c()
    pups.new <- c()
    if(any(c(plows,pups)>Pmax) | any(c(plows,pups)<Pmin)){
        for(i0 in 1:length(plows)){
            if((pups[i0]<Pmin & plows[i0]<Pmin) | (plows[i0]>Pmax & pups[i0]>Pmax)){
                plows.new <- plows.new
                pups.new <- pups.new
            }else if(plows[i0]<=Pmin & pups[i0]>=Pmin & pups[i0]<=Pmax){
                plows.new <- c(plows.new,Pmin)
                pups.new <- c(pups.new,pups[i0])
            }else if(plows[i0]<=Pmin & pups[i0]>Pmax){
                plows.new <- c(plows.new,Pmin)
                pups.new <- c(pups.new,Pmax)
            }else if(plows[i0]>=Pmin & pups[i0]<=Pmax){
                plows.new <- c(plows.new,plows[i0])
                pups.new <- c(pups.new,pups[i0])
            }else if(plows[i0]>=Pmin & plows[i0]<=Pmax & pups[i0]>Pmax){
                plows.new <- c(plows.new,plows[i0])
                pups.new <- c(pups.new,Pmax)
            }
        }
    }else{
        plows.new <- plows
        pups.new <- pups
    }
    return(list(plow=plows.new,pup=pups.new))
}

###tell foreach how to combine output
comb <- function(x, ...) {
      mapply(rbind,x,...,SIMPLIFY=FALSE)
}
####calculate timie
time.calc <- function(t1){
    dur <- as.numeric(proc.time()[3]-t1[3])
    time.consumed <- paste(floor(dur/3600),'h',floor((dur%%3600)/60),'m',dur%%60,'s',sep='')
    cat('Time consumed: ',time.consumed,'\n')
}
####weighted time binning
wtb <- function(t,x,ex,dt=1,sj=0){
#    t <- sort(t)
    t0 <- t[1]#the start time point
###note: weight w=1/(ex^2+sj^2)
    ts <- t[1]
    xs <- x[1]
    exs <- ex[1]
    tnew <- c()
    xnew <- c()
    exnew <- c()
    for(i1 in 2:length(t)){
        if((t[i1]-t0)<dt & i1<length(t)){
            ts <- c(ts,t[i1])
            xs <- c(xs,x[i1])
            exs <- c(exs,ex[i1])
        }else{
            tnew <- c(tnew,mean(ts))
            xnew <- c(xnew,xs%*%(1/(exs^2+sj^2))/sum(1/(exs^2+sj^2)))
#            xnew <- c(xnew,mean(xs))
            exnew <- c(exnew,sqrt(1/sum(1/(exs^2+sj^2))))
            ts <- t[i1]
            xs <- x[i1]
            exs <- ex[i1]
            t0 <- t[i1]
        }
    }
    return(cbind(tnew,xnew,exnew))
}
wtb.simple <- function(t,x,ex,dt=1,sj=0,N=NULL){
#    t <- sort(t)
    t0 <- t[1]#the start time point
###note: weight w=1/(ex^2+sj^2)
    tmin <- min(t)
    tmax <- max(t)
    if(is.null(N)){
        ts <- seq(tmin,tmax,by=dt)
    }else{
        ts <- seq(tmin,tmax,length.out=N+1)
    }
    out <- c()
    for(j in 1:(length(ts)-1)){
        if(j==1){
            inds <- which(t>=ts[j] & t<=ts[j+1])
        }else{
            inds <- which(t>ts[j] & t<=ts[j+1])
        }
        if(length(inds)>0){
            tmean <- mean(t[inds])
            xmean <- x[inds]%*%(1/(ex[inds]^2+sj^2))/sum(1/(ex[inds]^2+sj^2))
            exmean <- sqrt(1/sum(1/(ex[inds]^2+sj^2)))
            out <- rbind(out,c(tmean,xmean,exmean))
        }
    }
    return(out)
}

#####extract variable values from file names
extract.Nsamp <- function(f){
    Ns <- c()
    for(i in 1:length(f)){
        f1 <- gsub('.*Nsamp','',f[i])
        Ns <- c(Ns,as.integer(gsub('_.+','',f1)))
    }
    return(Ns)
}
extract.tem <- function(f){
    f1 <- gsub('.*tem','',f)
    return(as.numeric(gsub('_acc\\d+','',f1)))
}
extract.Ndata <- function(f){
    if(grepl('Ndata',f)){
        f1 <- gsub('.*Ndata','',f)
    }else{
        f1 <- gsub('.*w0_N','',f)
    }
    return(as.integer(gsub('_.+','',f1)))
}
extract.Lmax <- function(f){
    f1 <- gsub('.+negLmax','',f)
    f2 <- gsub('-.+','',f1)
    f3 <- gsub('.pdf','',f2)
    return(-as.numeric(f3))
}
extract.acc <- function(f){
    f1 <- gsub('.*acc','',f)
    return(as.numeric(f1))
}
combine.index <- function(ind1,ind2,norm=FALSE){
    if(length(ind1)!=0){
        if(!is.matrix(ind1)) ind1 <- matrix(ind1,ncol=1)
        if(!is.matrix(ind2)) ind2 <- matrix(ind2,ncol=1)
        if(nrow(ind2)<nrow(ind1)){
            ind <- cbind(ind1,c(ind2,rep(NA,nrow(ind1)-nrow(ind2))))
        }else if(nrow(ind2)>nrow(ind1)){
            ind1 <- rbind(ind1,matrix(NA,nrow=nrow(ind2)-nrow(ind1),ncol=ncol(ind1)))
            ind <- cbind(ind1,ind2)
        }else{
            ind <- cbind(ind1,ind2)
        }
    }else{
        ind <- ind2
    }
    if(norm){
        ind <- scale(ind)
    }
    if(!is.matrix(ind)){
        ind <- matrix(ind,ncol=1)
    }
    return(ind)
}
####period par transformation
period.transformation <- function(pers,period.type=period.par){
    if(period.type=='logP'){
        Pers <- exp(pers)
    }else if(period.type=='nu'){
        Pers <- 1/pers
    }else if(period.type=='P'){
        Pers <- pers
    }else{
        Pers <- NA
    }
    return(Pers)
}
period.trans2 <- function(ps,period.type=period.par){
    if(period.type=='logP'){
        Pers <- log(ps)
    }else if(period.type=='nu'){
        Pers <- 1/ps
    }else if(period.type=='P'){
        Pers <- ps
    }else{
        Pers <- NA
    }
    return(Pers)
}

K2msini <- function(K,P,e,Ms){
    G <- 4*pi^2*(1.496e11)^3/(365.25*24*3600)^2#m^3*s^-2*Msun^-1
    Me2s <- 3.003e-6#Earth mass in solar unit
    Mj2s <- 1/1048
    Mpj <- K*Ms/Mj2s*(2*pi*G*Ms/(P*24*3600))^(-1/3)*(1-e^2)^0.5#
    Mpe <- Mpj*Mj2s/Me2s
    ap <- ((P/yr2d)^2*Ms)^(1/3)#au
#    as <- K*P*24*3600/(2*pi)/1.496e11#m
#    Mp2 <- Ms*(as/ap)*1048#Mj; for circular orbit
#    Mp3 <- K*(P/365.25)^(1/3)*Ms^(2/3)/28.4#Mj; for circular orbit
#    return(c(Mpj,Mpe,ap))
    return(list(mj=Mpj,me=Mpe,a=ap))
}

fix <- function(x,par){
    if(length(par)>0){
        for(j in 1:length(par)){
            ind <- grep(par[j],names(startvalue))
            if(is.matrix(x)){
                x[ind,] <- 0
                x[,ind] <- 0
            }else{
#                if(grepl('per',par[j])){
                    x[ind] <- startvalue[ind]
#                }else{
#                    x[ind] <- 0
#                }
            }
        }
    }
    return(x)
}
lag.sta <- function(Ms,m1,a1,e1,m2,e2){
    mu1 <- m1/Ms
    mu2 <- m2/Ms
    alpha <- mu1+mu2
    gamma1 <- (1-e1^2)^0.5
    gamma2 <- (1-e2^2)^0.5
    f <- function(delta){
        alpha^-3*(mu1+mu2/delta^2)*(mu1*gamma1+mu2*gamma2*delta)^2-1-3^(4/3)*mu1*mu2/alpha^(4/3)
    }
    d <- uniroot.all(f,c(0,10))
    return(d^2*a1)
}
show.peaks <- function(ps,powers,levels){
    ind <- which(powers==max(powers) | (powers>(max(powers)-log(100)) & powers>levels[3]))
    pmax <- ps[ind]
    ppmax <- powers[ind]
    j0 <- 1
    p0 <- pmax[1]
    pp0 <- ppmax[1]
    pms <- p0
    pos <- pp0
    if(length(pmax)>1){
        for(j in 2:length(pmax)){
            if(abs(pmax[j]-p0) < 0.1*p0){
                if(ppmax[j]>pp0){
                    j0 <- j
                    p0 <- pmax[j]
                    pp0 <- ppmax[j0]
                    pms[length(pms)] <- p0
                    pos[length(pos)] <- pp0
                }
		    }else{
                j0 <- j
                p0 <- pmax[j]
                pp0 <- ppmax[j0]
                pms <- c(pms,p0)
                pos <- c(pos,pp0)
            }
        }
    }else{
        pms <- pmax
        pos <- ppmax
    }
    return(cbind(pms,pos))
}
check.window <- function(t,Dt,Nbin){
    n <- Nbin-1
    dt <- (max(t)-min(t)-Dt)/n
    tstart <- min(t)+(0:n)*dt
    tend <- min(t)+(0:n)*dt+Dt
    Ns <- c()
    for(j in 1:n){
        ind <- which(t>tstart[j] & t<tend[j])
        Ns <- c(Ns,length(ind))
    }
    return(Ns)
}

combine.list <- function(sets){
    N <- length(sets)
    data <- c()
    for(j in 1:N){
        data <- c(data,sets[[j]])
    }
    return(data)
}
