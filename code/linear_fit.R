###Fit linear function to each season of HARPS data
tab <- read.table('../../dwarfs/bary/data/HD128621_HARPScons.rv',header=TRUE)
inds <- which(diff(tab[,1])>100)
ind2 <- c(inds,nrow(tab))
ind1 <- c(1,inds+1)
fout <- '../../dwarfs/bary/data/HD128621_linear.pdf'
cat(fout,'\n')
pdf(fout,16,16)
par(mfrow=c(4,4))
for(j in 1:length(inds)){
#      cat('ind1=',ind1[j],'\n')
#      cat('ind2=',ind2[j],'\n')
      y <- tab[ind1[j]:ind2[j],2]
      t <- tab[ind1[j]:ind2[j],1]
      plot(t,y,xlab='tauE',ylab='RV [m/s]',main=paste(j,'season'))
      fit <- lm(y~t)
      res <- fit$residuals
      tab[ind1[j]:ind2[j],2] <- res
}
plot(tab[,1],tab[,2],xlab='tauE',ylab='RV [m/s]',main='all')
abline(v=0.5*(tab[ind1[-1],1]+tab[ind2[-length(ind2)],1]),col='red')
dev.off()
fout1 <- '../../dwarfs/bary/data/HD128621_HARPSlin.rv'
write.table(tab,file=fout1,quote=FALSE,row.names=FALSE)
out <- wtb(tab[,1],tab[,2],tab[,3],dt=15/24/60)
fout2 <- '../../dwarfs/bary/data/HD128621_HARPSbin15min.rv'
cat(fout1,'\n')
cat(fout2,'\n')
write.table(out,file=fout2,quote=FALSE,row.names=FALSE)
