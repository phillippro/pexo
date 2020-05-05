res <- Data[,2]-mcmc$res[[paste0('sig',Nsig)]]$model[,2]
tauE <- mcmc$res[[paste0('sig',Nsig)]]$model[,1]
for(star in Par$stars){
    ind <- which(Data$star==star)
    out <- cbind(tauE[ind],res[ind],Data[ind,3])
    fout <- paste0('../results/',star,'_pexo.dat')
    cat(fout,'\n')
    write.table(out,file=fout,quote=FALSE,row.names=FALSE,col.names=c('tE','RV','eRV'))
}
