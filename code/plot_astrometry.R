dir <- '../results/'
fname <- paste0(dir,star,'_df',df,'_tbase',tbase,'_type',type,'.pdf')
cat('output pdf:\n',fname,'\n')
pdf(fname,8,8)
par(mfrow=c(2,2),mar=c(5,5,1,1))
Ntime <- nrow(data)
Ncol <- ncol(data)
JDutc <-rowSums(utc)
pos.tempo <- data[,Ncol-(2:0)]
plot(pos.tempo[,1],pos.tempo[,2],xlab='x',ylab='y',main='TEMPO2')
plot(pos.tempo[,1],pos.tempo[,3],xlab='x',ylab='z',main='TEMPO2')
plot(uOT[,1],uOT[,2],xlab='x',ylab='y',main='PEXO')
plot(uOT[,1],uOT[,3],xlab='x',ylab='z',main='PEXO')
dx <- uOT[,1]-pos.tempo[,1]
dy <- uOT[,2]-pos.tempo[,2]
dz <- uOT[,3]-pos.tempo[,3]
plot(dx,dy,xlab='dx',ylab='dy')
plot(dx,dz,xlab='dx',ylab='dz')
du <- sqrt(dx^2+dy^2+dz^2)
plot(JDutc,du,xlab='JD[UTC]',ylab=expression(Delta*theta))
dev.off()
