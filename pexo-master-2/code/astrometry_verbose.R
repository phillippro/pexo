rel.pos <- astro_Relative(rBT=OutTime$BT[,1:3],rSO=OutObs$SO[,1:3]/au2km,tS=OutTime$tS,Par)
rel1 <-sqrt(rel.pos$xi1^2+rel.pos$eta1^2)/1e3
rel2 <-sqrt(rel.pos$xi2^2+rel.pos$eta2^2)/1e3
rel3 <-sqrt(rel.pos$xi3^2+rel.pos$eta3^2)/1e3
ind <- which(abs(OutTime$elevation*180/pi)>10)
cat('\nAstrometry effects:\n')
cat('Second order aberration:',gen_CalAmp(sqrt(rowSums(OutAstro$OffAbe1-OutAstro$OffAbe2)^2)),'as\n')
cat('Third order aberration:',gen_CalAmp(sqrt(rowSums((OutAstro$OffAbe-OutAstro$OffAbe2)^2))),'as\n')
if(Par$PlanetShapiro){
ns <- c('Sun','Mercury','Venus','Earth','Moon','Mars','Jupiter','Saturn','Uranus','Neptune')
for(n in ns){
    cat('Lensing by',n,':',gen_CalAmp(gen_CalLen(OutAstro$SolarDefList[[n]]))*pc2au,'as\n')
}
cat('Lensing by all SS planets and the Moon:',gen_CalAmp(gen_CalLen(OutAstro$SolarDef-OutAstro$SolarDefList$Sun))*pc2au,'as\n')
}
cat('Target system lensing:',gen_CalAmp(sqrt(rowSums((OutAstro$OffLenT)^2))),'as\n')
cat('Second order geometric effect:',gen_CalAmp(rel2),'as\n')
cat('Third order geometric effect:',gen_CalAmp(rel3),'as\n')
cat('Atmospheric refraction:',gen_CalAmp(OutAstro$R*pc2au),'as\n')
