library(RColorBrewer)
library(plotly)
#library(MathJaxR)
library(astro)
api <- FALSE
#source('mcmc_func1.R')
source('mcmc_func.R')
Sys.setenv('MAPBOX_TOKEN' = 'pk.eyJ1IjoiZmFibyIsImEiOiJjam84c2YzZGYwMXpuM3JsYzczemdrMGloIn0.TVLQwOZsCuRRO6M9eqKCpQ')
if(api){
    Sys.setenv("plotly_username"="fabo")
    Sys.setenv("plotly_api_key"="tZIFLORwVdD9kTztrBbu")
}

###settings
compute_bins <- function(x, n) {
    list(
        start = min(x),
        end = max(x),
        size = (max(x) - min(x)) / n
    )
}
Nbin <- 20
threshold <- 1
####
Mj2s <- 9.543e-4#Jupiter mass in unit of solar mass
Mj2e <- 317.8#Jupiter mass in unit of Earth mass
yy <- ms <- exp(seq(log(1e-3*Mj2s),log(1e3*Mj2s),by=0.1))
Ms <- ms+1#total mass
#as <- exp(seq(log(1e-3),log(1e3),by=0.1))
Ps <-  exp(seq(log(1e-3),log(1e10),by=0.1))#day
a <- (Ps/365.25)^(2/3)#au
es <- c(0, 0.5, 0.9)
pm <- list()
omega.dot.th <- 10#c(1,10,100)
for(k in 1:length(es)){
    tmp <- 7.78/(1-es[k]^2)/(a/0.05)/Ps
    mts <- omega.dot.th/tmp#solar mass
    pm[[k]] <- cbind(Ps,mts)
}
###read table
tab <- read.csv('planets_nasa.csv')
###classification
ind.transit <- which(tab[,'pl_discmethod']=='Transit')
ind.rv <- which(tab[,'pl_discmethod']=='Radial Velocity')
ind.imaging <- which(tab[,'pl_discmethod']=='Imaging')
ind.etv <- which(tab[,'pl_discmethod']=='Eclipse Timing Variations')
ind.astrometry <- which(tab[,'pl_discmethod']=='Astrometry')
ind.obm <- which(tab[,'pl_discmethod']=='Orbital Brightness Modulation')
ind.ptv <- which(tab[,'pl_discmethod']=='Pulsation Timing Variations')
ind.microlensing <- which(tab[,'pl_discmethod']=='Microlensing')
ind.ttv <- which(tab[,'pl_discmethod']=='Transit Timing Variations')
ind.pt <- which(tab[,'pl_discmethod']=='Pulsar Timing')
#ind.select <- c(ind.transit,ind.rv)
ind.select <- 1:nrow(tab)

#####calculate omega.dot
p <- as.numeric(tab[ind.select,'pl_orbper'])
a <- as.numeric(tab[ind.select,'pl_orbsmax'])
e <- as.numeric(tab[ind.select,'pl_orbeccen'])
e[is.na(e)] <- 0
ms <- as.numeric(tab[ind.select,'st_mass'])
mp <- as.numeric(tab[ind.select,'pl_bmassj'])
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
ind.gr <- ind.select[which(omega.dot>omega.dot.th)]

tab[ind.select,'pl_orbper'] <- out[,1]
tab[ind.select,'pl_orbeccen'] <- out[,3]
ind.select <- ind.select[tab[ind.select,'pl_orbper']<1e2 & tab[ind.select,'st_mass']>0.1]
#ind.gr <- ind.select[which(omega.dot>100)]

ep <- c(0.2056,0.0068,0.0167,0.0934,0.0484,0.0542,0.0472,0.0086)#0.2488
Pp <- c(87.969,224.701,365.256,686.98,11.862*365.256,29.457*365.256,84.011*365.256,164.79*365.256)#247.68*365.256
pname <- c('Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune')
mp <- c(0.0553,0.815,1.000,0.107,317.83,95.159,14.536,17.147)/317.83 #Mj#0.0021
###add solar system planets
if(TRUE){
    tmp <- array(NA,dim=c(8,ncol(tab)))
    colnames(tmp) <- colnames(tab)
    tmp[,'pl_name'] <- pname
    tmp[,'pl_hostname'] <- 'Sun'
    tmp[,'pl_orbper'] <- Pp
    tmp[,'pl_bmassj'] <- mp
    tmp[,'pl_orbeccen'] <- ep
    tmp[,'pl_discmethod'] <- 'Solar System Planets'
    tmp[,'st_mass'] <- 1
    tmp[,'st_teff'] <- 5772
    tab <- rbind(tab,tmp)
}
tmp <- rep('NSS',nrow(tab))
tmp[length(tmp)-(7:0)] <- 'SS'
ind.exo <- 1:(nrow(tab)-8)
ind.mass <- which(!is.na(tab[,'pl_bmassj']))
ind.sun <- nrow(tab)-(7:0)
tab <- cbind(tab,type=tmp)#,pname=c(rep('',length(ind.exo)),pname))
####methods
methods <- as.character(unique(tab[,'pl_discmethod']))
ppos <- c('middle right','bottom right','middle right','middle right','middle right','middle right','bottom right','middle right')
q <- plot_ly(type = 'scatter', mode = 'markers')%>%add_trace(x = as.numeric(tab[ind.select,'pl_orbper']),y=as.numeric(tab[ind.select,'st_mass']),color=~as.character(tab[ind.select,'pl_discmethod']))%>%add_trace(x = as.numeric(tab[ind.gr,'pl_orbper']),y=as.numeric(tab[ind.gr,'st_mass']),name='Planets with precession >10 degree/century', marker = list(symbol='circle-open',size = 10), color = I('blue'))
xlim <- range(as.numeric(tab[ind.select,'pl_orbper']),na.rm =TRUE)
ylim <- range(as.numeric(tab[ind.select,'st_mass']),na.rm =TRUE)
dashes <- c('','dash','dot')
ns <- paste('e=',es)
for(j in 1:length(pm)){
    ind <- which(pm[[j]][,1]>xlim[1] & pm[[j]][,1]<xlim[2] & pm[[j]][,2]>ylim[1] & pm[[j]][,2]<ylim[2])
    q <- add_trace(q,x = pm[[j]][ind,1], y=pm[[j]][ind,2],color=I('black'),name = paste0('10 degree/century precession for e=',es[j]),line=list(dash=dashes[j]),mode = 'lines')
}

q <- layout(q,xaxis = list(type = "log",exponentformat='e',title = "Orbital Period [day]",showexponent='all',dtick=1),yaxis = list(type = "log",exponentformat='e',title = "Stellar Mass [solar mass]",showexponent='all',dtick=1),traceorder='grouped')

##combine plots
#s <- subplot(p1,p2,p3,p4,nrows = 2, heights = c(0.2, 0.8), widths = c(0.8, 0.2), margin = 0,shareX = TRUE, shareY = TRUE,titleY = TRUE,titleX = TRUE)
#p <- layout(s, showlegend = FALSE)
# Create a shareable link to your chart
# Set up API credentials: https://plot.ly/r/getting-startedj
if(api){
    chart_link <- api_create(p, filename="exoplanet_relativity_test")
}
#orca(p, "exo_dist.jpg")

#orca(s, fpair)

# Create a shareable link to your chart
# Set up API credentials: https://plot.ly/r/getting-started
if(FALSE){
chart_link <- api_create(p, filename = "exoplanet_barycorr_bias")
chart_link
}


