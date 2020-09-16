if(file.exists('mcmc_func.R')){
    source('mcmc_func.R')
}else{
    source('../code/mcmc_func.R')
}
if(!exists('outf')) outf <- TRUE
#if(!exists('outf')) outf <- FALSE
fit <- fit_LogLike(Data,OutObs,RateObs=RateObs,ParFit=ParOpt,Par=Par,OutTime0=OutTime0,TimeUpdate=TRUE)
model <- fit$model
llmax <- round(fit$llike)

dir.out <- '../results/'
####plot
fname <- paste0(dir.out,star,'_',opt$component,'_Nmax',Par$Nmax,'_llmax',llmax,'_N',Par$Niter,'_einstein',Par$Einstein)
fpdf <- paste0(fname,'.pdf')
cat(fpdf,'\n')
pdf(fpdf,16,16)
par(mfrow=c(4,4))
for(star in Par$stars){
    index <- which(Data$type=='rv' & Data$star==star)
    if(length(index)>0){
        tauE <- model[index,1]
        rv <- Data[index,2]-model[index,2]#residual
        erv <- Data[index,3]
        ins <- Data[index,'instrument']
	ind <- sort(tauE,index.return=TRUE)$ix
        if(outf){
            out <- cbind(tauE,rv,erv,ins)[ind,]
            fout <- paste0(dir.out,star,'_comb.rv')
            cat(fout,'\n')
            write.table(out,file=fout,quote=FALSE,row.names=FALSE,col.names=c('tauE','RV','eRV'))
            if(star=='HD128620'){
                ind <- which((out[,2]-mean(out[,2]))<3*sd(out[,2]) & out[,2]> -15)
            }else{
                ind <- which((out[,2]-mean(out[,2]))<3*sd(out[,2]))
            }
        }
###save residual RVs
        inss <- unique(Data$instrument[index])
        for(instr in inss){
            ind <- index[Data[index,'instrument']==instr]
            plot(Data[ind,1],Data[ind,2],xlab='jd',ylab='rv',main=paste('rv for',instr,star),xlim=range(model[ind,1],Data[ind,1]),ylim=range(model[ind,2],Data[ind,2]))
                                        #                points(model[ind,1],model[ind,2],col='red')
            points(Data[ind,1],model[ind,2],col='red')
            plot(Data[ind,1],Data[ind,2]-model[ind,2],xlab='jd',ylab='rv',main=paste('rv residual for',instr,star,';sd(res)=',round(sd(Data[ind,2]-model[ind,2]),3)))
            tauE <- model[ind,1]
            rv <- Data[ind,2]-model[ind,2]#residual
            erv <- Data[ind,3]
            if(outf){
                out <- cbind(tauE,rv,erv)
###             fout <- paste0('../../dwarfs/bary/data/',star,'_',instr,'.rv')
                fout <- paste0(dir.out,star,'_',instr,'.rv')
                cat(fout,'\n')
                write.table(out,file=fout,quote=FALSE,row.names=FALSE,col.names=c('tauE','RV','eRV'))
            }
            ind <- which((out[,2]-median(out[,2]))<5*sd(out[,2]))
            out1 <- out[ind,]
        }
    }

    inds <- which(Data$type=='abs' & Data$star==star)
    if(length(inds)>0){
        inss <- unique(Data[inds,'instrument'])
        for(instr in inss){
            index <- inds[Data[inds,'instrument']==instr]
            plot(Data[index,2]*180/pi,Data[index,4]*180/pi,xlab='ra[deg]',ylab='dec[deg]',main=paste(instr,'absolute astrometry for',star),xlim=range(model[index,2]*180/pi,Data[index,2]*180/pi),ylim=range(model[index,4]*180/pi,Data[index,4]*180/pi))
            points(model[index,2]*180/pi,model[index,4]*180/pi,col='red')
            plot((Data[index,2]-model[index,2])*206264.8,(Data[index,4]-model[index,4])*206264.8,xlab='ra[as]',ylab='dec[as]',main=paste(instr,'astrometry residual for',star))
            plot(Data[index,1],(Data[index,2]-model[index,2])*206264.8,xlab='JD',ylab='ra[as]',main=paste(instr,'astrometry residual for',star))
            plot(Data[index,1],(Data[index,4]-model[index,4])*206264.8,xlab='JD',ylab='dec[as]',main=paste(instr,'astrometry residual for',star))
            if(outf){
                ##save residual RVs
                dra <- (Data[index,2]-model[index,2])*cos(model[index,4])/DMAS2R
                ddec <- (Data[index,4]-model[index,4])/DMAS2R
                tauE <- model[index,1]
                out <- cbind(tauE,dra,Data[index,3],ddec,Data[index,5])
                fout <- paste0(dir.out,star,'_',instr,'.abs')
                cat(fout,'\n')
                write.table(out,file=fout,quote=FALSE,row.names=FALSE,col.names=c('tauE','dra','era','ddec','edec'))
                fout1 <- paste0(dir.out,star,'_',instr,'_dRA.abs')
                cat(fout1,'\n')
                write.table(out[,1:3],file=fout1,quote=FALSE,row.names=FALSE,col.names=c('tauE','dra','era'))
                fout2 <- paste0(dir.out,star,'_',instr,'_dDEC.abs')
                cat(fout2,'\n')
                write.table(out[,c(1,4:5)],file=fout2,quote=FALSE,row.names=FALSE,col.names=c('tauE','ddec','edec'))
            }
        }
    }
}
####trace and posterior plots
if(any(Par$types=='rv')) source('mcmc_analysis.R')

inds <- which(Data$type=='rel')
if(length(inds)>0 & outf){
    inss <- unique(Data[inds,'instrument'])
    for(instr in inss){
        index <- inds[Data[inds,'instrument']==instr]
        plot(Data[index,2],Data[index,4],xlab='ra*',ylab='dec',main=paste(instr,'relative astrometry'),xlim=range(model[index,2],Data[index,2]),ylim=range(model[index,4],Data[index,4]))
        points(model[index,2],model[index,4],col='red')
        plot(Data[index,2]-model[index,2],Data[index,4]-model[index,4],xlab='ra*',ylab='dec',main=paste(instr,'relative astrometry residual'))
        dra <- Data[index,2]-model[index,2]
        ddec <- Data[index,4]-model[index,4]
        tauE <- model[index,1]
        out <- cbind(tauE,dra,Data[index,3],ddec,Data[index,5])
        fout <- paste0(dir.out,star,'_',instr,'.rel')
        cat(fout,'\n')
        write.table(out,file=fout,quote=FALSE,row.names=FALSE,col.names=c('tauE','dra','era','ddec','edec'))
        fout1 <- paste0(dir.out,star,'_',instr,'_dRA.rel')
        cat(fout1,'\n')
        write.table(out[,1:3],file=fout1,quote=FALSE,row.names=FALSE,col.names=c('tauE','dra','era'))
        fout2 <- paste0(dir.out,star,'_',instr,'_dDEC.rel')
        cat(fout2,'\n')
        write.table(out[,c(1,4:5)],file=fout2,quote=FALSE,row.names=FALSE,col.names=c('tauE','ddec','edec'))
    }
}
dev.off()

###save all outputs as R objects
fobj <- paste0(fname,'.Robj')
if(grepl('Robj',opt$out)){
    fobj <- opt$out
}
cat('\nOutput Robj file:\n')
cat(fobj,'\n')
save(list=ls(all=TRUE),file=fobj)


