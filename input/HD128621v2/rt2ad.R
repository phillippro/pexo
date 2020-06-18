options(digits=20)
source('../../code/timing_function.R')
source('../../code/general_function.R')
tab <- read.table('k17_astrometry.dat',header=TRUE)
yr <- tab[,1]
instr <- c('NACO',paste0('ALMA.B',1:10),'An15','An14','An12','An06','An11','An08','WDS','Hip','J94','P91','T85','F86','J95','SOFI','2MASS')
ab <- c(7,4,3,2,1.6,1.3,1.0,0.7,0.45,0.35)
wl <- c(2.17,ab*1e3,rep(0.55,13),1.7,2.159)#micrometer
ind <- match(tab[,'ref'],instr)
ind <- ind[!is.na(ind)]
jd <- time_Yr2jd(yr)
rho <- tab[,'rho']*1e3#mas
erho <- tab[,'erho']*1e3#mas
theta <- tab[,'theta']#deg
etheta <- tab[,'etheta']#deg
tmp <- gen_Rt2ad(rho,erho,theta,etheta)
write.table(cbind(rowSums(jd),tmp,wl[ind]),file='HD128621_K17.rel',quote=FALSE,row.names=FALSE,col.names=c('JD', 'dRA', 'edRA', 'dDE', 'edDE','wl'))
