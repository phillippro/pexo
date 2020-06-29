args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
   f <- args[1]
}else{
   f <- '../../dwarfs/bary/data/HD128621_HARPScons.rv'
}
tab <- read.table(f,header=TRUE)
fit <- lm(tab[,2]~poly(tab[,1],3))
write.table(cbind(tab[,1],fit$residuals,tab[,3]),file=f,quote=FALSE,row.names=FALSE,col.names=c('BJD','RV','eRV'))