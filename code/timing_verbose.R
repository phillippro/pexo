optional <- FALSE
cat('\nTiming effects:\n')
if(Par$CompareT2){
    cat('Difference in Roemer delay using uSB and uST as the reference direction:',gen_CalAmp(OutTime$RoemerSB-OutTime$RoemerSolar),'second\n')
}
cat('Instant second order roemer delay in the solar system:',gen_CalAmp(OutTime$RoemerOrder$Roemer2),'second\n')
if(!all(OutTime$RoemerT2$all==0)){
    cat('Cummulative third order roemer delay in the solar system:',gen_CalAmp(OutTime$RoemerSB-OutTime$RoemerT2$all),'second\n')
}
dt.einsteinS <- rowSums(OutObs$JDtdb-OutObs$JDtt)*DAYSEC
cat('Einstein delay in the solar system:',max(abs(dt.einsteinS)),'second\n')
cat('Shapiro delay due to all solar system objects:',gen_CalAmp(OutTime$ShapiroSolar),'second\n')

if(Par$PlanetShapiro){
    ns <- c('Sun','Mercury','Venus','Earth','Moon','Mars','Jupiter','Saturn','Uranus','Neptune')
}else{
    ns <- 'Sun'
}
for(n in ns){
    cat('Shapiro delay due to',n,':',gen_CalAmp(OutTime$ShapiroPlanet[[n]]),'second\n')
}
if(Par$PlanetShapiro) cat('Shapiro delay due to all planets:',gen_CalAmp(OutTime$ShapiroSolar-OutTime$ShapiroPlanet[[1]]),'second\n')
Roemer0 <- (OutObs$SO[,1]/au2km*Par$pqu[1,3]+OutObs$SO[,2]/au2km*Par$pqu[2,3]+OutObs$SO[,3]/au2km*Par$pqu[3,3])*AULT#second
cat('Delay due to proper motion:',gen_CalAmp(OutTime$RoemerSolar-Roemer0),'second\n')
if(Par$binary){
    cat('Roemer delay in the target system:',gen_CalAmp(OutTime$RoemerTarget),'second\n')
    cat('Einstein delay in the target system:',gen_CalAmp(OutTime$EinsteinTarget),'second\n')
    cat('Shapiro delay in TS:',gen_CalAmp(OutTime$ShapiroTarget),'second\n')
}else{
    cat('Roemer delay in the target system:',0,'second\n')
    cat('Einstein delay in the target system:',0,'second\n')
    cat('Shapiro delay in TS:',0,'second\n')
}

ind <- which(abs(OutTime$elevation*180/pi)>10)
if(length(ind)>0){
    cat('Delay due to atmospheric effect:',gen_CalAmp(OutTime$TropoDelay[ind]),'second\n')
}
