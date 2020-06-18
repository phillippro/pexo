###Read timing or data files
Data <- c()
lambda <- list()
if(opt$mode=='emulate'){
    if(!file.exists(opt$time)){
        s <- unlist(strsplit(opt$time,split=' '))
        s <- s[s!='']
        StrError <- 'Error: the UTC timing data cannot be found or the UTC sequence cannot be generated from the argument after -t !\n'
        if(length(s)>=3){
            ts <- try(seq(as.numeric(s[1]),as.numeric(s[2]),by=as.numeric(s[3])),TRUE)
            if(class(ts)=='try-error'){
                stop(StrError,call. = FALSE)
            }else{
                utc <- time_ChangeBase(cbind(time_ToJD(ts),0))
            }
        }else{
            stop(StrError,call. = FALSE)
        }
    }else{
        utc <- read.table(opt$time)
#        instrument <- gsub('.+_|\\..+','',opt$time)
        if(ncol(utc)==1){
            utc <- cbind(utc,0)
        }
        utc <- time_ToJD(utc)
        utc <- time_ChangeBase(utc,1)
    }
    Data <- cbind(rowSums(utc),NA,NA,NA,NA,Par$star,'rv',instrument,0.5)
    Par$ins[[Par$star]]$rv <- opt$ins
}else{
    Nt <- 0
    for(kk in 1:length(stars)){
        star <-  stars[kk]
        fin <- paste0('../input/',stars[kk])
        fs <- list.files(path=fin,full.name=TRUE)
        ind.delay <- grep('delay$',fs)
        Ndelay <- length(ind.delay)
        Par[[star]]$NsetDelay <- length(ind.delay)
        inss <- c()
        if(Ndelay>0){
            cat(Ndelay,'.delay file found in ',fin,'\n')
            for(j in ind.delay){
                ins <- gsub('.delay|.+_','',fs[j])
#                if(any(ins==opt$ins)){
                if(TRUE){
                    inss <- c(inss,ins)
                    tab <- read.table(fs[j])
                    if(class(tab[2,1])!='numeric'){
                        tab <- read.table(fs[j],header=TRUE)
                    }
###                tab[,2] <- tab[,2]-mean(tab[,2])#automatically remove offset
                    tmp <- cbind(time_ToJD(tab[,1]),tab[,2:3],NA,NA,as.character(star),'delay',ins,NA)
                    Data <- rbind(Data,as.matrix(tmp))
                }else{
                    if(opt$verbose) cat('No .delay file found in ',fin,'!\n')
                }
            }
        }else{
            if(opt$verbose) cat('No .delay file found in ',fin,'!\n')
        }
        Par$ins[[star]]$delay <- inss

        ind.rv <- grep('rv$',fs)
        Nrv <- length(ind.rv)
        Par[[star]]$NsetRv <- length(ind.rv)
        inss <- c()
        if(length(ind.rv)>0 & grepl('R',Par$component)){
            cat(Nrv,'.rv file found in ',fin,'!\n')
            for(j in ind.rv){
                ins <- gsub('.rv|.+_','',fs[j])
#                if(any(ins==opt$ins)){
                if(TRUE){
                    inss <- c(inss,ins)
                    tab <- read.table(fs[j])
                    if(class(tab[2,1])!='numeric'){
                        tab <- read.table(fs[j],header=TRUE)
                    }
                    rvbar <- mean(tab[,2])
                                        #                tab[,2] <- tab[,2]-rvbar#automatically remove offset
                    tmp <- cbind(time_ToJD(tab[,1]),tab[,2:3],NA,NA,as.character(star),'rv',ins,NA)
                    Data <- rbind(Data,as.matrix(tmp))
                }else{
                    if(opt$verbose) cat('No .rv file found in ',fin,'!\n')
                }
            }
        }else{
            if(opt$verbose) cat('No .rv file found in ',fin,'!\n')
        }
        Par$ins[[star]]$rv <- inss

        ind.abs <- grep('abs$',fs)
        Nabs <- length(ind.abs)
        Par[[star]]$NsetAbs <- length(ind.abs)
        inss <- c()
        if(length(ind.abs)>0 & grepl('A',Par$component)){
            cat(Nabs,'.abs file found in ',fin,'!\n')
            for(j in ind.abs){
                ins <- gsub('.abs|.+_','',fs[j])
#                if(any(ins==opt$ins)){
                if(TRUE){
                    inss <- c(inss,ins)
                    tab <- read.table(fs[j])
                    if(class(tab[2,1])!='numeric'){
                        tab <- read.table(fs[j],header=TRUE)
                    }
                    tab[,c(2,4)] <- tab[,c(2,4)]/180*pi#from deg to rad
                    Nd <- nrow(tab)
                    if(ncol(tab)>5){
                        tmp <- cbind(time_ToJD(tab[,1]),tab[,2:5],as.character(star),'abs',ins,tab[,6])
                    }else{
                        tmp <- cbind(time_ToJD(tab[,1]),tab[,2:5],as.character(star),'abs',ins,NA)
                    }
                    Data <- rbind(Data,as.matrix(tmp))
                }else{
                    if(opt$verbose) cat('No .abs file found in ',fin,'!\n')
                }
            }
        }else{
            if(opt$verbose) cat('No .abs file found in ',fin,'!\n')
        }
        Par$ins[[star]]$abs <- inss

        ind.rel <- grep('rel$',fs)
        Nrel <- length(ind.rel)
        if(Nrel==0 & opt$verbose) cat('No .rel file found in ',fin,'!\n')
        Par[[star]]$NsetRel <- length(ind.rel)
        inss <- c()
        if(length(ind.rel)>0 & grepl('A',Par$component)){
            cat(Nrel,'.rel file found in ',fin,'!\n')
            for(j in ind.rel){
                ins <- gsub('.rel|.+_','',fs[j])
                inss <- c(inss,ins)
                tab <- read.table(fs[j])
                if(class(tab[2,1])!='numeric'){
                    tab <- read.table(fs[j],header=TRUE)
                }
                Nd <- nrow(tab)
                if(ncol(tab)>5){
                    tmp <- cbind(time_ToJD(tab[,1]),tab[,2:5],as.character(star),'rel',ins,tab[,6])
                }else{
                    tmp <- cbind(time_ToJD(tab[,1]),tab[,2:5],as.character(star),'rel',ins,NA)
                }
                Data <- rbind(Data,as.matrix(tmp))
            }
        }else{
            if(opt$verbose) cat('No .rel file found in ',fin,'!\n')
        }
        Par$ins[[star]]$rel <- inss
#        if(Nrv==0 & Nabs==0 & Nrel==0 & Ndelay==0) stop('No data file found in',fin,'\n')
    }
}
colnames(Data) <- c('utc','V1','eV1','V2','eV2','star','type','instrument','wavelength')
Data <- fit_changeType(Data,c(rep('n',5),rep('c',3),'n'))
utc <- time_ChangeBase(cbind(Data[,1],0),1)
Par$Dtypes <- unique(Data$type)
Par$tmin <- min(rowSums(utc))
Par$Nepoch <- nrow(utc)

###get index for different types and sets of data
Par$targets <- unique(Data[,'star'])
Par$types <- unique(Data[,'type'])
Par$instruments <- unique(Data[,'instrument'])
index <- list()
for(star in Par$targets){
    index[[star]] <- list()
    index[[star]]$all <- which(Data[,'star']==star)
    index[[star]]$types <- unique(Data[Data[,'star']==star,'type'])
    for(type in Par$types){
        index[[star]][[type]] <- list()
        index[[star]][[type]]$instruments <- unique(Data[Data[,'star']==star & Data[,'type']==type,'instrument'])
        index[[star]][[type]]$all <- which(Data[,'star']==star & Data[,'type']==type)
        for(ins in Par$instruments){
            index[[star]][[type]][[ins]] <- which(Data[,'star']==star & Data[,'instrument']==ins & Data[,'type']==type)
        }
    }
}
Par$index <- index
