Npar <- ncol(mc)-2
###add onto previous plots rather than generating additional pdfs
#ff <- paste0(fname,'_mcmc.pdf')
#cat(ff,'\n')
#pdf(ff,16,16)
#par(mfrow=c(4,4))
cn <- colnames(mc)[1:Npar]
if(Par$Niter>=1e5){
    inds <- sort(sample(1:nrow(mc),1e4))
}else{
    inds <- 1:nrow(mc)
}
for(j in 1:Npar){
      plot(inds,mc[inds,j],xlab='Iteration',ylab=cn[j])
      plot(inds,mc[inds,'loglike'],xlab='Iteration',ylab='Log Likelihood')
      plot(mc[inds,j],mc[inds,'loglike'],xlab='Iteration',ylab='Log Likelihood')
      hist(mc[inds,j],xlab=cn[j],ylab='Freq.')
}
###pair plot
pairs(mc[inds,1:Npar])
#dev.off()
