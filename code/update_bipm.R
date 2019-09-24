#dir <- '../data/'
fs <- list.files(dir,'TTBIPM\\.')
yr <- gsub('.+\\.','',fs)
yr.now <- as.integer(gsub('-.+','',Sys.Date()))
if(!any(yr!=yr.now)){
    try(system(paste0('wget ftp://ftp2.bipm.org/pub/tai/ttbipm/TTBIPM.',yr.now,' -P ',dir)),TRUE)
}
fs <- list.files(dir,'TTBIPM\\.',full.name=TRUE)
yr <- gsub('.+\\.','',fs)
fnow <- paste0(dir,'TTBIPM.',yr.now)
if(length(fs)==fnow){
    fbipm <- fnow
}else{
    fbipm <- fs[which.max(yr)]
}
Par$BIPM <- max(yr)
str <- readLines(fbipm)
ind <- grep('^42589',str)
BIPMfile <- c()
for(j in ind:length(str)){
    s <- unlist(strsplit(str[j],' '))
    BIPMfile <- rbind(BIPMfile,as.numeric(s[s!='']))
}
colnames(BIPMfile) <- c('MJD','TT-EAL-32.184(microsecond)', 'TT-TAI-32.184(microsecond)')

