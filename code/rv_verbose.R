rel.pos <- astro_Relative(rBT=OutTime$BT[,1:3],rSO=OutObs$SO[,1:3]/au2km,tS=OutTime$tS,Par)
rel1 <-sqrt(rel.pos$xi1^2+rel.pos$eta1^2)/1e3
rel2 <-sqrt(rel.pos$xi2^2+rel.pos$eta2^2)/1e3
rel3 <-sqrt(rel.pos$xi3^2+rel.pos$eta3^2)/1e3
uOT <- OutTime$uOT
uST <- OutTime$uST
du <- uOT-uST#real difference
rSO <- OutObs$SO[,1:3]/au2km#au
rSO.perp <- rSO-rowSums(rSO*uST)*uST#ref. Wright & Eastman 2014, page 6
RST <- gen_CalLen(OutTime$rST)
du1 <- -rSO.perp/(RST*pc2au)
du2 <- rSO.perp^2/(RST*pc2au)/2/Cauyr
ddu <- (du-du1)
drv2 <- OutRv$RvTot*sqrt(rowSums(ddu^2))
cat('\nRV effects:\n')
cat('Relativistic effects in solar system:',gen_CalAmp(OutRv$RvgsO),'m/s\n')
if(Par$Lensing){
    if(Par$PlanetShapiro){
        ns <- c('Sun','Mercury','Venus','Earth','Moon','Mars','Jupiter','Saturn','Uranus','Neptune')
    }else{
        ns <- 'Sun'
    }
    for(n in ns){
        cat('Lensing shift by',n,':',gen_CalAmp(OutRv$Zcomb$Zlensing[[n]])*CMPS,'m/s\n')
    }
    cat('Lensing by all SS planets and the Moon:',gen_CalAmp(OutRv$RvlO -OutRv$Zcomb$Zlensing$Sun*CMPS),'m/s\n')
}
cat('SR effect in TS:',gen_CalAmp(OutRv$RvsT),'m/s\n')
cat('GR effect in TS:',gen_CalAmp(OutRv$RvgT),'m/s\n')
cat('lensing effect in target system:',gen_CalAmp(OutRv$RvlT),'m/s\n')
cat('Second order geometric effects in TS:',gen_CalAmp(rel2*OutRv$RvTot/pc2au),'m/s\n')
cat('Third order geometric effects in TS:',gen_CalAmp(rel3*OutRv$RvTot/pc2au),'m/s\n')
ind10 <- which(abs(OutTime$elevation*180/pi)>10)
ind5 <- which(abs(OutTime$elevation*180/pi)>5)
cat('Tropospheric RV:',gen_CalAmp(OutRv$RvTropo[ind10]),'m/s for elevation > 10 deg\n')
cat('Tropospheric RV:',gen_CalAmp(OutRv$RvTropo[ind5]),'m/s for elevation > 5 deg\n')
cat('Tropospheric RV:',gen_CalAmp(OutRv$RvTropo),'m/s for all elevation angles\n')
