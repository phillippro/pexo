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
xx <- as <- exp(seq(log(1e-3),log(1e3),by=0.1))
ds <- c(1, 10, 100)
plxs <- 1000/ds
vtot <- 5e4#m/s
pm <- list()
pm2 <- list()
pm3 <- list()
pm4 <- list()
for(k in 1:length(ds)){
    plx <- 1000/ds[k]#mas
    vs.sqa <- 2*pi*ms/sqrt(Ms)
    tmp <- 1e-3/vtot*1e3*206265/plx/vs.sqa
    as <- 1/tmp^2#1 cm/s threshold
    tmp2 <- 2e-3/vtot*1e3*206265/plx/vs.sqa
    as2 <- 1/tmp2^2#1 cm/s threshold
    as3 <- 1e-2/vtot*1e3*206265/plx*Ms/ms#1 cm/s threshold
    as4 <- 2e-3/vtot*1e3*206265/plx*Ms/ms#1 cm/s threshold

    Ps <- sqrt(as^3/Ms)*365.25#day
    Ps2 <- sqrt(as2^3/Ms)*365.25#day
    Ps3 <- sqrt(as3^3/Ms)*365.25#day
    Ps4 <- sqrt(as4^3/Ms)*365.25#day

    ms2 <- ms/Mj2s#Mjup
    pm[[k]] <- cbind(Ps,ms2)
    pm2[[k]] <- cbind(Ps2,ms2)
    pm3[[k]] <- cbind(Ps3,ms2)
    pm4[[k]] <- cbind(Ps4,ms2)
}
###read table
tab <- read.csv('planets_nasa.csv')
p <- as.numeric(tab[,'pl_orbper'])
a <- as.numeric(tab[,'pl_orbsmax'])
ms <- as.numeric(tab[,'st_mass'])
mp <- as.numeric(tab[,'pl_bmassj'])
###find barycorr biased targets
out <- cbind(p,a,ms,mp)
colnames(out) <- c('P[day]','a[au]','m[Msol]','m[Mjup]')
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
p <- out[,1]/365.25
a <- out[,2]
ms <- out[,3]
mp <- out[,4]*Mj2s
as <- a*mp/(mp+ms)
dvr.period <- as*vtot/1e3/206265*1000/tab[,'st_dist']
vs <- 2*pi*as/p
dpm <- vs*plx
dvr.trend <- vtot*dpm/1e3/206265
ind.per <- which(dvr.period>0.01)
ind.trend <- which(dvr.trend>0.001)
tab[,'pl_orbper'] <- out[,1]

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
p <- plot_ly(type = 'scatter', mode = 'markers')%>%add_trace(x = as.numeric(tab[ind.exo,'pl_orbper']),y=as.numeric(tab[ind.exo,'pl_bmassj']),color=~as.character(tab[ind.exo,'pl_discmethod']))%>%add_trace(x = as.numeric(tab[ind.trend,'pl_orbper']),y=as.numeric(tab[ind.trend,'pl_bmassj']),name='Planets with >1 mm/s/yr trend bias',marker = list(symbol='circle-open',size = 10,color='blue'))%>%add_trace(x = as.numeric(tab[ind.per,'pl_orbper']),y=as.numeric(tab[ind.per,'pl_bmassj']),name='Planets with >1 cm/s periodic bias', marker = list(symbol='circle-open',size = 10, color = 'orange'))
ylim <- range(as.numeric(tab[ind.exo,'pl_bmassj']),na.rm =TRUE)
xlim <- range(as.numeric(tab[ind.exo,'pl_orbper']),na.rm =TRUE)
dashes <- c('','dash','dot')
ns <- paste(ds,'pc')
for(j in 1:length(pm)){
    ind <- which(pm[[j]][,1]>xlim[1] & pm[[j]][,1]<xlim[2] & pm[[j]][,2]>ylim[1] & pm[[j]][,2]<ylim[2])
    p <- add_trace(p,x = pm[[j]][ind,1], y=pm[[j]][ind,2],color=I('black'),name = paste0('1 mm/s/yr trend bias for d=',ns[j]),line=list(dash=dashes[j]),mode = 'lines')
}
if(threshold>1){
for(j in 1:length(pm)){
    ind2 <- which(pm2[[j]][,1]>xlim[1] & pm2[[j]][,1]<xlim[2] & pm2[[j]][,2]>ylim[1] & pm2[[j]][,2]<ylim[2])
    p <- add_trace(p,x = pm2[[j]][ind2,1], y=pm2[[j]][ind2,2],color=I('darkgrey'),name = paste0('2 mm/s/yr trend bias for d=',ns[j]),line=list(dash=dashes[j]),mode = 'lines')
}
}
for(j in 1:length(pm)){
    ind <- which(pm3[[j]][,1]>xlim[1] & pm3[[j]][,1]<xlim[2] & pm3[[j]][,2]>ylim[1] & pm3[[j]][,2]<ylim[2])
    p <- add_trace(p,x = pm3[[j]][ind,1], y=pm3[[j]][ind,2],color=I('darkgrey'),name = paste0('1 cm/s periodic bias for d=',ns[j]),line=list(dash=dashes[j]),mode = 'lines')
}
if(threshold>1){
for(j in 1:length(pm)){
    ind2 <- which(pm4[[j]][,1]>xlim[1] & pm4[[j]][,1]<xlim[2] & pm4[[j]][,2]>ylim[1] & pm4[[j]][,2]<ylim[2])
    p <- add_trace(p,x = pm4[[j]][ind2,1], y=pm4[[j]][ind2,2],color=I('pink'),name = paste0('2 mm/s threshold for d=',ns[j]),line=list(dash=dashes[j]),mode = 'lines')
}
}

p <- layout(p,xaxis = list(type = "log",exponentformat='e',title = "Orbital Period [day]",showexponent='all',dtick=1),yaxis = list(type = "log",exponentformat='e',title = "Mass [M<sub>Jup</sub>]",showexponent='all',dtick=1),traceorder='grouped')

##combine plots
#s <- subplot(p1,p2,p3,p4,nrows = 2, heights = c(0.2, 0.8), widths = c(0.8, 0.2), margin = 0,shareX = TRUE, shareY = TRUE,titleY = TRUE,titleX = TRUE)
#p <- layout(s, showlegend = FALSE)
# Create a shareable link to your chart
# Set up API credentials: https://plot.ly/r/getting-startedj
if(api){
    chart_link <- api_create(p, filename="exoplanet_barycorr_bias")
}
#orca(p, "exo_dist.jpg")

#orca(s, fpair)

# Create a shareable link to your chart
# Set up API credentials: https://plot.ly/r/getting-started
if(FALSE){
chart_link <- api_create(p, filename = "exoplanet_barycorr_bias")
chart_link
}


