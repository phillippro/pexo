EpochHip <- 2448349.0625
EpochGaia <- 2457206.375
strs <- unique(Data[,c('star','instrument','type')])
ss <- c()
for(j in 1:nrow(strs)){
    if(strs[j,3]=='rv'){
        ss <- c(ss,paste0('jitterRv.',strs[j,1],'.',strs[j,2],' 0 -1e6 1e6 U'))
        ss <- c(ss,paste0('bRv.',strs[j,1],'.',strs[j,2],' 0 -1e6 1e6 U'))
    }else if(strs[j,3]=='rel' | strs[j,3]=='abs'){
        ss <- c(ss,paste0('logjitterAstro.',strs[j,1],'.',strs[j,2],' 0 -10 10 U'))
#        ss <- c(ss,paste0('xtel.',strs[j,2],' 0'))
#        ss <- c(ss,paste0('ytel.',strs[j,2],' 0'))
#        ss <- c(ss,paste0('ztel.',strs[j,2],' 0'))
    }

    if(strs[j,2]=='AAT'){
        ss <- c(ss,'phi.AAT -31.27704')
        ss <- c(ss,'elong.AAT 149.0661')
        ss <- c(ss,'height.AAT 1.164')
    }
    if(strs[j,2]=='PFS'){
        ss <- c(ss,'phi.PFS -29.013983')
        ss <- c(ss,'elong.PFS -70.692633')
        ss <- c(ss,'height.PFS 2.40792')
    }
    if(strs[j,2]=='MIKE'){
        ss <- c(ss,'phi.MIKE -29.013983')
        ss <- c(ss,'elong.MIKE -70.692633')
        ss <- c(ss,'height.MIKE 2.40792')
    }
    if(strs[j,2]=='KECK'){
        ss <- c(ss,'phi.KECK 19.82636')
        ss <- c(ss,'elong.KECK -155.47501')
        ss <- c(ss,'height.KECK 4.145')
    }
    if(strs[j,2]=='SOPHIE' | strs[j,2]=='ELODIE'){
        ss <- c(ss,paste0('phi.',strs[j,2],' 43.92944'))
        ss <- c(ss,paste0('elong.',strs[j,2],' 5.7125'))
        ss <- c(ss,paste0('height.',strs[j,2],' 0.65'))
    }
    if(grepl('CARM',strs[j,2])){
        ss <- c(ss,paste0('phi.',strs[j,2],' 37.220791'))
        ss <- c(ss,paste0('elong.',strs[j,2],' -2.546847'))
        ss <- c(ss,paste0('height.',strs[j,2],' 2.168'))
    }
    if(strs[j,2]=='UVES'){
        ss <- c(ss,'phi.UVES -24.625407')
        ss <- c(ss,'elong.UVES -70.403015')
        ss <- c(ss,'height.UVES 2.648')
    }
    if(strs[j,2]=='HARPN'){
        ss <- c(ss,'phi.HARPN 28.754')
        ss <- c(ss,'elong.HARPN -17.88814')
        ss <- c(ss,'height.HARPN 2.37')
    }
    if(strs[j,2]=='APF'){
        ss <- c(ss,'phi.APF 37.3425')
        ss <- c(ss,'elong.APF -121.63825')
        ss <- c(ss,'height.APF 1.274')
    }
    if(grepl('HARPS',strs[j,2])){
        ss <- c(ss,paste0('phi.',strs[j,2],' -29.2584'))
        ss <- c(ss,paste0('elong.',strs[j,2],' -70.7345'))
        ss <- c(ss,paste0('height.',strs[j,2],' 2.4'))
    }
    if(strs[j,2]=='WDS'){
        ss <- c(ss,'phi.WDS 0')
        ss <- c(ss,'elong.WDS 0')
        ss <- c(ss,'height.WDS 0')
    }
}
####add astrometric data
astro  <- read.csv('../data/RVstarInfo.csv')
star.name <- c(gsub(' ','',Par$star),gsub('GJ','GL',Par$star),gsub('-','',Par$star))
an <- read.table('../data/alias_name.txt')
for(s in star.name){
    ind <- which(s==gsub(' ','', astro[,'target']) | s==gsub(' ','', astro[,'ID']) | s==gsub(' ','', astro[,'StarKnown']))
    index <- match(s,an[,1])
    index <- index[!is.na(index)]
    if(length(ind)==0 & length(index)>0){
       s1 <- an[index,2]
       cat(s1,'\n')
       ind <- which(s1==gsub(' ','', astro[,'target']) | s1==gsub(' ','', astro[,'ID']) | s1==gsub(' ','', astro[,'StarKnown']))
    }
    if(length(ind)>0) break()
}
cname <- c('ra','dec','pmra','pmdec','parallax','radial_velocity')
star <- Par$star
if(length(ind)==0){
    cat('No information found in ../data/RVstarInfo.csv, please add other names in ../data/alias_name.txt in order to use Gaia astrometry; otherwise Hipparcos astrometry will be used!')
    cat('case0\n')
    source('simbad_cor.R')
}else{
    cat('case1\n')
    ra <- astro[ind[1],'RA']
    dec <- astro[ind[1],'DEC']
    if(is.na(ra))  source('simbad_cor.R')
}
if(!exists('StarAstro')){
if(is.na(ra)) stop(paste0('No position information found for ',star,'!'))
if(!is.na(astro[ind[1],'pmra'])){
    cat('case2\n')
    StarAstro <- as.numeric(astro[ind[1],cname])
    if(is.na(StarAstro[6])) StarAstro[6] <- as.numeric(as.character(astro[ind[1],'RVsimbad']))
    if(is.na(StarAstro[6])) StarAstro[6] <- 0
    epoch <- EpochGaia
    if(!is.na(astro[ind[1],'Ms'])) ss <- c(ss,paste0('mT ',as.character(astro[ind[1],'Ms'])))
}else{
    cat('case3\n')
   cat('No Gaia astrometry info found in ../data/RVstarInfo.csv and looking for astrometry in revised Hipparcos catalog!\n')
    astro2 <- read.csv('../data/hip2.csv',sep='|',check.names=FALSE)
   ang <- sqrt((ra-astro2[,'_RAJ2000'])^2+(dec-astro2[,'_DEJ2000'])^2)
    ind <- which(ang<1/60)
   if(length(ind)>0){
       index <- ind[which.min(ang[ind])]
       StarAstro <- c(as.numeric(astro2[index,c('RArad','DErad','pmRA','pmDE','Plx')]),0)
       epoch <- 2448348.750
   }else{
       stop(paste0('No astrometry found for ',star,' in Gaia and Hipparcos!'))
   }
}
}
if(Par$star=='HD128620' | Par$star=='HD128621'){
StarAstro <- c(219.917756905051,-60.83698524600996,-3641.862228256955,701.5804096412021,740.8767738648643,-20)
epoch <- 2448348.750
}
names(StarAstro) <- cname
s <- paste(c('ra','dec','pmra','pmdec','plx','rv'),StarAstro)
ss <- c(ss,paste0('epoch ',epoch))
ss <- c(ss,s)

###add planetary parameters if there is not any
##find signal
if(Par$Nmax>0){
    ss <- c(ss,paste('logmC',opt$orbit['logmC'],'-10 10 U'))
    ss <- c(ss,paste('logP',opt$orbit['logP'],'-10 20 U'))
    ss <- c(ss,paste('e',opt$orbit['e'],'0 1 U'))
    ss <- c(ss,paste('I',opt$orbit['I'],'0 180 U'))
    ss <- c(ss,paste('omegaT',opt$orbit['omegaT'],'0 360 U'))
    if(!grepl('R',opt$component)){
       if(opt$orbit['Omega']<180) ss <- c(ss,paste('Omega',opt$orbit['Omega'],'0 180 U'))
       if(opt$orbit['Omega']>=180) ss <- c(ss,paste('Omega',opt$orbit['Omega'],'180 360 U'))
    }else{
       ss <- c(ss,paste('Omega',opt$orbit['Omega'],'0 360 U'))
    }
    ss <- c(ss,paste('Tp',opt$orbit['Tp'],opt$orbit['Tp']-1e5,opt$orbit['Tp']+1e5,'U'))
    ss <- c(ss,'BinaryModel kepler')
}else{
    ss <- c(ss,'BinaryModel none')
}
tmp <- c(tmp,ss)
