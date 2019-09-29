if(dt>0){
tab <- readLines(paste0('tempo2/g2vary',CLK,'e-',-log10(dt),'rv',df,'type',type,'.out'))
}else{
tab <- readLines(paste0('tempo2/g2vary',CLK,'dt0rv',df,'type',type,'.out'))
}
ind1 <- grep('^Starting',tab)
ind2 <- grep('^Finished',tab)
data <- c()#jd.tt (2 part),bjd.tcb (2 part)
for(i in (ind1+1):(ind2-1)){
    tmp <- unlist(strsplit(tab[i],split=' '))
    tmp <- tmp[tmp!='']
    utcs <- gen_One2two(tmp[1])
    tcb <- gen_One2two(tmp[2])
    tt <- gen_One2two(tmp[3])
    s <- as.numeric(tmp[4])
    sJ <- as.numeric(tmp[5])
    sS <- as.numeric(tmp[6])
    sV <- as.numeric(tmp[7])
    sU <- as.numeric(tmp[8])
    sN <- as.numeric(tmp[9])
    if(dt>0){
        data <- rbind(data,c(tt,tcb,utcs,s,sJ,sS,sV,sU,sN))
    }else{
        pos <- as.numeric(tmp[10:12])
        data <- rbind(data,c(tt,tcb,utcs,s,sJ,sS,sV,sU,sN,pos))
    }
}
data[,c(1,3,5)] <- data[,c(1,3,5)]+DJM0

jd.tts0 <- floor(rowSums(data[,1:2]))
jd.tts1 <- (data[,1]-jd.tts0)+data[,2]
jd.tts <- cbind(jd.tts0,jd.tts1)

bjd.tcbs0 <- floor(rowSums(data[,3:4]))
#bjd.tcbs1 <- (data[,3]-bjd.tcbs0)+data[,4]
bjd.tcbs1 <- (data[,3]-jd.tts0)+data[,4]
bjd.tcbs <- cbind(jd.tts0,bjd.tcbs1)

jd.utcs0 <- floor(rowSums(data[,5:6]))
jd.utcs1 <- (data[,5]-jd.utcs0)+data[,6]
jd.utcs <- cbind(jd.utcs0,jd.utcs1)

if(FALSE){
    dutc <- jd.utcs[,2]-utc[,2]
    jd.tts[,2] <- jd.tts[,2]-dutc
    bjd.tcbs[,2] <- bjd.tcbs[,2]-dutc
}

#zg <- numDeriv.2part(bjd.tcbs,jd.tts,N=3)-1

