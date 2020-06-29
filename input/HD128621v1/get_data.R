source('../../code/timing_function.R')
source('../../code/sofa_function.R')
tab <- readLines('wds14396-6050.txt')
out <- c()
for(j in 1:length(tab)){
tab[j] <- gsub(':','',tab[j])
tmp <- unlist(strsplit(tab[j],' '))
tmp <- tmp[tmp!='' & tmp!='\\:']
ind <- which(tmp=='.')
if(length(ind)>0) tmp[ind] <- 'NA'
out <- rbind(out,tmp[1:5])
}
ind <- which(out[,1]==1991.25)
write.table(out,file='HD128621_WDS.dat',quote=FALSE,row.names=FALSE,col.names=c('date','theta','terr','rho','rerr'))
tab <- read.table('HD128621_WDS.dat',header=TRUE)

pdf('test.pdf',6,6)
rho <- tab[,'rho']
erho <- tab[,'rerr']
th <- tab[,'theta']/180*pi
eth <- tab[,'terr']/180*pi
t <- tab[,'date']
t <- rowSums(time_Yr2jd(t))
thip <- time_Yr2jd(1991.25)
x <- rho*sin(th)#alpha*
y <- rho*cos(th)#delta
ex <- sqrt((rho*cos(th)*eth)^2+(erho*sin(th))^2)
ey <- sqrt((rho*sin(th)*eth)^2+(erho*cos(th))^2)
ind <- which(is.na(x) | is.na(y))
if(length(ind)>0){
    t <- t[-ind]
    x <- x[-ind]
    y <- y[-ind]
    ex <- ex[-ind]
    ey <- ey[-ind]
}
ex[is.na(ex)] <- 1
ey[is.na(ey)] <- 1

plot(x,y,xlab='alpha*',ylab='delta',xlim=c(-30,20),ylim=c(-20,30))
dev.off()

ind <- which.min(abs(sum(thip)-t))
cat('ind=',ind,'\n')
write.table(cbind(t,x*1e3,ex*1e3,y*1e3,ey*1e3),file='HD128621_WDS.rel',quote=FALSE,row.names=FALSE,col.names=c('JD','alpha','ealpha','delta','edelta'))
