library(RColorBrewer)
library(plotly)
#library(MathJaxR)
library(astro)
api <- FALSE
#source('mcmc_func1.R')
source('mcmc_func.R')
###read table
tab <- read.csv('exoplanet_eu.csv',sep=',')
#####calculate omega.dot
p <- as.numeric(tab[,'orbital_period'])
a <- as.numeric(tab[,'semi_major_axis'])
e <- as.numeric(tab[,'eccentricity'])
e[is.na(e)] <- 0
ms <- as.numeric(tab[,'star_mass'])
mp <- as.numeric(tab[,'mass'])
eomega.min <- as.numeric(tab[,'omega_error_min'])
eomega.max <- as.numeric(tab[,'omega_error_max'])
####
Mj2s <- 9.543e-4#Jupiter mass in unit of solar mass
Mj2e <- 317.8#Jupiter mass
###find barycorr biased targets
out <- cbind(p,a,e,ms,mp)
colnames(out) <- c('P[day]','a[au]','e','m[Msol]','m[Mjup]')
for(j in 1:nrow(out)){
    ind <- which(is.na(out[j,]))
    if(length(ind)>0){
        if(any(ind==1)){
            out[j,1] <- sqrt(a[j]^3/ms[j])*365.25
        }
        if(any(ind==2)){
            out[j,2] <- ((p[j]/365.25)^2*ms[j])^{1/3}
        }
    }
}
p <- out[,1]#day
a <- out[,2]
e <- out[,3]
ms <- out[,4]
mp <- out[,5]*Mj2s
as <- a*mp/(mp+ms)
omega.dot <- 7.78/(1-e^2)*ms/(a/0.05)/p
#ind.pre <- which(omega.dot>1 & eomega.max>0 & eomega.max<1 & tab[,'detection_type']=='Radial Velocity')
ind.pre <- which(omega.dot>1 & eomega.max>0 & eomega.max<1)
cbind(omega.dot[ind.pre],tab[ind.pre,1:30],tab[ind.pre,'detection_type'])
