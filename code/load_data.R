#dir <- '../data/'
dir <- '../../pexo/data/'
###read data
#eop <- read.table(paste0(dir,'dut1_orientation.txt'),header=TRUE,check.names=FALSE)
EopFile <- paste0(dir,'eopc04_IAU2000.62-now')
tmp <- readLines(EopFile)
cname <- tmp[11]
cname <- gsub(' Err','_Err',cname)
cname <- gsub('Date','Y M D',cname)
cname <- unlist(strsplit(cname,' '))
cname <- cname[cname!='']
eop <- read.table(EopFile,skip=14)
f1900 <- '../data/38_EOP_C01.1900-NOW_V2013_0138.txt'
if(!file.exists(f1900)){
    cat('Warn: EOP parameter file',f1900,' for epochs from 1900 to 1962 cannot be found\n')
}else{
    eop1900 <- read.table(f1900,header=TRUE)
    mjd1900 <- time_Jd2mjd(time_Yr2jd(eop1900[,'an']))
    xfunc <- approxfun(mjd1900,eop1900[,2])
    yfunc <- approxfun(mjd1900,eop1900[,4])
    UT1mTAIfunc <- approxfun(mjd1900,eop1900[,6])#second
}
                                        #eop1846 <- read.table('../data/186_EOP_C01_2000.1846_NOW_V2013_01186.txt',header=TRUE)
colnames(eop) <- cname

####other files
if(Par$TtType=='BIPM'){
     source('update_bipm.R')
}
##Dates and Delta(AT)s
delat <- read.csv(paste0(dir,'tai_utc.csv'))
dirDE <- paste0(dir,'de',Par$DE)
dirt <- paste0(dirDE,'t')
if(!file.exists(dirDE) & !file.exists(dirt)){
    cat('\n',dirDE,'or',dirt,'is not found and try to download from JPL website!\n')
    system(paste0('source download_ephemerides.sh ',Par$DE,'t'),intern=TRUE)
    if(!file.exists(dirt)){
        system(paste('source download_ephemerides.sh',Par$DE),intern=TRUE)
        if(file.exists(dirDE)) system(paste('Rscript reformat_ephemerides.R',Par$DE))
    }else{
        system(paste0('Rscript reformat_ephemerides.R ',Par$DE,'t'),intern=TRUE)
    }
}
if(!file.exists(dirt) & !file.exists(dirDE)){
#    dirs <- list.dirs(path=dir)
#    dirs <- gsub('\\..\\/data|\\/','',dirs)
#    cat('dirs=',dirs,'\n')
    stop(paste0('In',opt$par,', DE is not a valid JPL version of ephemerides'))
}

if(!file.exists(dirt)){
    cat(paste0('Warn: DE',Par$DE,'t is not found and the FB01 method is used to calculate TDB-TT!\n'))
    Par$TtTdbMethod <- 'FB01'
    dirt <- dirDE
}
fs <- list.files(dirt,pattern='DE.+dat',full.name=TRUE)
f1 <- fs[1]
DEfile <- read.table(f1)
fh <- list.files(dirt,pattern='header',full.name=TRUE)
indH <- sort(fh,decreasing=TRUE,index.return=TRUE)$ix[1]
cat('read header:',fh[indH],'\n')
DEheader <- readLines(fh[indH])
#DEheader <- readLines(paste0(dirt,'header.',Par$DE,'t'))
##for nutation calculation
xls <- read.csv(paste0(dir,'coeff_nut00a1.csv'))
xpl <- read.csv(paste0(dir,'coeff_nut00a2.csv'))
##
#fb <- read.table('../data/TDB.1950.2050')


