####RV decoupling bias
mC <- AstroPar['m2']
mT <- AstroPar['m1']
if(mC==0){
    a <- 1e6#au; because there is no companion in
    plx <- 1e-6#mas
}else{
    a <- par.comp['a.au']
    plx <- AstroPar['plx1']#mas
}
#Dt <- gen_CalAmp(rowSums(utc))/DJY#yr
Dt <- 10#yr
dvr.trend <- 1.52*mC*(mC+mT)^(-0.5)*a^(-0.5)*plx*Dt*1e-3#m/s
dvr.period <- 0.24*a*plx*1e-3#m/s
dt.trend <- 15.17*mC*(mC+mT)^(-0.5)*a^(-0.5)*plx*Dt*1e-6#second
dt.period <- 2.4*a*plx*1e-6#second
du.trend <- 6.27*mC*(mC+mT)^(-0.5)*a^(-0.5)*plx*Dt*1e-3#as
du.period <- 0.99**a*plx*1e-3#as
cat('Timing trend bias',Dt,'yr:',dt.trend,'s\n')
cat('Timing period bias:',dt.period,'s\n')
cat('Astrometry trend bias',Dt,'yr:',du.trend,'as\n')
cat('Astrometry period bias:',du.period,'as\n')
cat('RV trend bias',Dt,'yr:',dvr.trend,'m/s\n')
cat('RV period bias:',dvr.period,'m/s\n')

