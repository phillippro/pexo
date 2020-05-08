source('../code/sofa_function.R')
source('../code/timing_function.R')
tab <- readLines('hipparcos-tle.txt')
ind1 <- grep('^1',tab)
jds <- c()
for(k in 1:length(ind1)){
    str <- unlist(strsplit(tab[ind1[k]],split=' '))
    str <- str[str!='']
    tmp <- as.numeric(str[5])
    yr <- floor(tmp/1e3)
    day <- tmp-1e3*yr
    yr <- yr+1900+day/365.25
    jds <- c(jds,rowSums(time_Yr2jd(yr)))
}
write.table(jds,file='hipparcos_tle_jd.txt',quote=FALSE,row.names=FALSE,col.names=FALSE)
