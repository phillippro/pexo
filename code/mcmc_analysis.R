Npar <- ncol(mc)-2
pdf('mcmc_trace.pdf',16,16)
par(mfrow=c(4,4))
cn <- colnames(mc)[1:Npar]
for(j in 1:Npar){
      plot(mc[,j],xlab='iteration',ylab=cn[j])
      hist(mc[,j],xlab=cn[j],ylab='Freq.')
}
dev.off()