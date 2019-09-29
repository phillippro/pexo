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
###
mT <- c(0.1,1,10)
mC <- exp(seq(log(1e-3),log(150),by=0.1))
#es <- c(0.1,0.2,0.8)
es <- rep(0.1,3)
mtot <- Ps.gr <- Ps.couple <- Ps.min <- array(NA,dim=c(length(mC),length(es)))
vr0 <- 0.01#m/s
vSB <- 50#km/s
for(j in 1:length(es)){
    mtot[,j] <- mC+mT[j]
    tmp1 <- 5.92*es[j]/(1-es[j]^2)*mtot[,j]^{-4/3}*mC^2
    tmp2 <- 9.94*es[j]/sqrt(1-es[j]^2)*mtot[,j]^{-2/3}*mC*vSB/50
    tmp3 <- 5.92*es[j]/(1-es[j]^2)*mtot[,j]^(-1/3)*mC
    Ps.min[,j] <- (tmp1/vr0)^1.5
    Ps.couple[,j] <- (tmp2/vr0)^3
    Ps.gr[,j] <- (tmp3/vr0)^1.5
}
tab <- read.csv('binary_malkov12_tab2.csv',sep='|')
e <- eccentricity <- tab[,'e']
mtot1 <- tab[,'Mdyn']
mC1 <- mT1 <- 0.5*mtot1
P1 <- tab[,'Per']#yr
dvmin <-5.92*e/(1-e^2)*P1^(-2/3)*mtot1^(-4/3)*mC1^2
dvmax <- dvmin+9.94*e/sqrt(1-e^2)^P1^(-1/3)*mtot1^(-2/3)*mC1*vSB/50
dvgr <- 5.92*e/(1-e^2)*P1^(-2/3)*mtot1^(-1/3)*mC1
ind.gr <- which(dvgr>1 & P1<10)
ind.min <- which(dvmin>1 & P1<10)
ind.max <- which(dvmax>1 & P1<10)

####methods
q <- plot_ly(type = 'scatter', mode = 'markers')%>%add_trace(x = tab[,'Per'],y=tab[,'Mdyn'],color=~eccentricity,name='Binary data',marker=list(opacity=0.5))
xlim <- range(tab[,'Per'],na.rm =TRUE)
ylim <- range(tab[,'Mdyn'],na.rm =TRUE)

#q <- add_trace(q,x = tab[ind.min,'Per'],y=tab[ind.min,'Mdyn'],name='Binaries with δv<sub>srT,min</sub>=1 cm/s', marker = list(symbol='circle-open',size = 10), color = I('blue'))

q <- add_trace(q,x = tab[ind.gr,'Per'],y=tab[ind.gr,'Mdyn'],name='δv<sub>grT</sub>>1 m/s and P<10 years', marker = list(symbol='circle-open',size = 10), color = I('blue'))

for(j in 1:length(es)){
    ind <- which(Ps.min[,j]>xlim[1] & Ps.min[,j]<xlim[2] & mtot[,j]>ylim[1] & mtot[,j]<ylim[2])
    q <- add_trace(q,x = Ps.min[ind,j],y=mtot[ind,j],name=paste0('m<sub>T</sub>=',round(mT[j],1),' M<sub>⊙</sub>, e=', round(es[j],2),', δv<sub>srT,min</sub>=1 cm/s'),mode = 'lines')
}
for(j in 1:length(es)){
    ind <- which(Ps.couple[,j]>xlim[1] & Ps.couple[,j]<xlim[2] & mtot[,j]>ylim[1] & mtot[,j]<ylim[2])
    q <- add_trace(q,x = Ps.couple[ind,j],y=mtot[ind,j],name=paste0('m<sub>T</sub>=',round(mT[j],1),' M<sub>⊙</sub>, e=', round(es[j],2),', δv<sub>srT,couple</sub>=1 cm/s'),mode = 'lines',line=list(dash='dash'))
}
for(j in 1:length(es)){
    ind <- which(Ps.gr[,j]>xlim[1] & Ps.gr[,j]<xlim[2] & mtot[,j]>ylim[1] & mtot[,j]<ylim[2])
    q <- add_trace(q,x = Ps.gr[ind,j],y=mtot[ind,j],name=paste0('m<sub>T</sub>=',round(mT[j],1),' M<sub>⊙</sub>, e=', round(es[j],2),', δv<sub>grT</sub>=1 cm/s'),mode = 'lines',line=list(dash='dot'))
}

q <- layout(q,xaxis = list(type = "log",exponentformat='e',title = "P [year]",showexponent='all',dtick=1),yaxis = list(type = "log",exponentformat='e',title = "m<sub>T</sub> + m<sub>C</sub> [M<sub>⊙</sub>]",showexponent='all',dtick=1),traceorder='grouped')
#


