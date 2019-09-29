source('Dtdb.R')
source('Dtdbv1.R')
jd <- tb$tt
cat('Calcuate TDB-TT using JPL ephemeris!\n')
tt2tdb.geo <- gen_CalEph(jd,body='TT-TDB (at geocenter)',DE=DE)
Dt.geo <- -tt2tdb.geo[,1]#TDB-TT; s
zgeo <- -tt2tdb.geo[,2]

##FB01
cat('Calcuate Teph-TT as a function of Teph using FB01 method!\n')
tmin <- min(rowSums(jd))
tmax <- max(rowSums(jd))
ddt <- min(1,Tstep/10)
maxday <- floor((tmax-tmin)/ddt)
jd.ref <- tmin+(0:maxday)*ddt
write.table(c(tmin,tmax,ddt,'F'),file='FB01/fb2001.in',quote=FALSE,row.names=FALSE,col.names=FALSE)
system('cd FB01 ; ./fb2001 <fb2001.in >fb2001.out')
fb <- read.table('FB01/fb2001.out',skip=3)
inds <- sapply(1:nrow(jd),function(i) which.min(abs(jd.ref-sum(jd[i,]))))
Dt.geo0 <- as.numeric(gsub('D','e',as.character(fb[inds,2])))#day;(Teph-TT)(Teph)
zgeo0 <- as.numeric(gsub('D','e',as.character(fb[inds,3])))
dDt.geo <- zgeo0*time_T2mT2(jd,time_ChangeBase(cbind(jd.ref[inds],0),1),toSecond=FALSE)#correction by shifting teph.ref to teph
Dt.geo0 <- (Dt.geo0+dDt.geo)*DAYSEC+IFTE.TEPH0#s; correct to Teph-TT

cat('Calcuate TDB-TT using FBsofa method for geocenter!\n')
obs <- time_Geo2obs(jd,utc,ObsPar,eop.type=eop.type,n=n)
tmp <- Dtdb(jd[,1],jd[,2],ut=rowSums(obs$ut1),elong=0,u=0,v=0)
tmp1 <- Dtdb(jd[,1],jd[,2]+1e-3,ut=rowSums(obs$ut1),elong=0,u=0,v=0)
tmp2 <- Dtdb(jd[,1],jd[,2]-1e-3,ut=rowSums(obs$ut1),elong=0,u=0,v=0)
zgeo1 <- -(tmp2[,1]-tmp1[,1])/(2e-3*DAYSEC)
zgeo2 <- tmp[,2]
Dt.geo1 <- Dt.geo2 <- tmp[,1]

source('Dtdbv1.R')
cat('Calcuate TDB-TT using FBsofa method for geocenter!\n')
tmp <- Dtdb(jd)
Dt.geo3   <- tmp[,1]
zgeo3 <- tmp[,2]

cat('head(Dt.geo)=',head(Dt.geo),'s\n')
cat('head(Dt.geo0)=',head(Dt.geo0),'s\n')
cat('head(Dt.geo1)=',head(Dt.geo1),'s\n')
cat('head(Dt.geo2)=',head(Dt.geo2),'s\n')
cat('head(Dt.geo3)=',head(Dt.geo3),'s\n\n')

cat('head(zgeo)=',head(zgeo),'\n')
cat('head(zgeo0)=',head(zgeo0),'\n')
cat('head(zgeo1)=',head(zgeo1),'\n')
cat('head(zgeo2)=',head(zgeo2),'\n')
cat('head(zgeo3)=',head(zgeo3),'\n')
