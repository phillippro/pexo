###This file contains plot functions
options(scipen=999)
options(digits=5)
library(gridExtra)
library(grid)
library(cowplot)
library(ggplot2)
library(lattice)

plot_OC1D <- function(raw,mp,res,xlab='JD [day]',ylab1='RV [m/s]',ylab2='O-C [m/s]',Nbin=10,FoldType='phase',alpha.data=0.2,alpha.fit=0.5,data.bin=NULL,res.bin=NULL,FitType='point',bsize=2,title=NULL){
########################################################################
## 1D scatter plot combined with residual plot
##
## input:
##    raw - raw data
##    mp - model prediction for best-fit parameters
##    res - residual
##    data.bin - binned data
##    res.bin - binned residual
##
## output:
##    ggplot
########################################################################
###raw data
###phase plot
    sets <- unique(raw[,'instrument'])
    Nset <- length(sets)
    colnames(mp)[1:3] <- colnames(res)[1:3] <- colnames(raw)[1:3] <- c('t','y','dy')
    tmin <- min(raw[,1])
    mp[,1] <- mp[,1]-tmin
    res[,1] <- raw[,1] <- raw[,1]-tmin
    xlab <- gsub('JD',paste0('JD-',round(tmin,2)),xlab)
    g1 <- ggplot(raw, aes(x=t, y=y))+geom_point(alpha=alpha.data)+geom_errorbar(aes(ymin=y-dy,ymax=y+dy),width=0.2,alpha = alpha.data)+theme_gray(base_size=bsize)
    if(Nset==1) g1 <- g1+theme(legend.position = "none")
    if(FitType=='line'){
        g2 <- g1+geom_line(data=mp,mapping=aes(x=t,y=y), color = "red",size=1,alpha=alpha.fit)+xlab(NULL)+ylab(ylab1)+ggtitle(title)
    }else{
        g2 <- g1+geom_point(data=mp,mapping=aes(x=t,y=y), color = "red",size=1,alpha=alpha.fit)+xlab(NULL)+ylab(ylab1)+ggtitle(title)
    }
    if(!is.null(data.bin)){
        FigRaw <- g2+geom_point(data=data.bin,aes(x=t,y=y,colour=instrument),size=3)+geom_errorbar(data=data.bin,aes(ymin=y-dy,ymax=y+dy),width=1)+coord_cartesian(xlim=range(raw[,1]), ylim = range(mp[,2],data.bin[,2],mean(raw[,2])+2*sd(raw[,2]),mean(raw[,2])-2*sd(raw[,2])))
    }else{
        FigRaw <- g2+coord_cartesian(xlim=range(raw[,1]),ylim = range(mp[,2],raw[,2]))
    }
###residual plot
    FigRes <- ggplot(res, aes(t, y))+geom_point(alpha =alpha.data)+xlab(xlab)+ylab(ylab2)+geom_errorbar(aes(ymin=y-dy,ymax=y+dy),alpha = alpha.data)+theme_gray(base_size=bsize)
    if(!is.null(res.bin)){
        ylim <- range(res.bin[,2],mean(res[,2])+2*sd(res[,2]),mean(res[,2])-2*sd(res[,2]))
        FigRes <- FigRes+geom_point(data=res.bin,aes(x=t,y=y),size=3)+geom_errorbar(data=res.bin,aes(ymin=y-dy,ymax=y+dy),width=1)+coord_cartesian(ylim = ylim)+ annotate("text", x = 0.9*max(res[,1]), y = ylim[2], label = paste0("RMS=",round(sd(res[,2]),2),' m/s'))
    }
###combined plot
    PlotComb <- plot_grid(
        plot_grid(
            FigRaw+ theme(legend.position = "none",plot.title = element_text(hjust = 0.5)),
            FigRes+theme(legend.position = "none"),
            align='v',ncol=1,axis = "lr",rel_heights = c(1,0.5)
        ), plot_grid(
            get_legend(FigRaw),
             ggplot()+theme_minimal()+geom_blank(),
             ncol =1,rel_heights = c(1,0.5)
           )
, rel_widths = c(9.5,0.5)
    )
    plot(PlotComb)
#    plot(p)
}

plot_OC2D <- function(raw,mp,res,xlab1='RA [deg]',xlab2='RA residual [deg]',ylab1='DEC [as]',ylab2='DEC residual [as]',Nbin=10,FoldType='phase',alpha.data=0.2,alpha.fit=0.5,data.bin=NULL,res.bin=NULL,FitType='point',bsize=1,title=NULL){
########################################################################
## 2D scatter plot combined with residual plot
##
## input:
##    raw - raw data
##    mp - model prediction for best-fit parameters
##    res - residual
##    data.bin - binned data
##    res.bin - binned residual
##
## output:
##    ggplot
########################################################################
###raw data
###phase plot
    sets <- unique(raw[,'instrument'])
    Nset <- length(sets)
    colnames(mp)[1:5] <- colnames(res)[1:5] <- colnames(raw)[1:5] <- c('t','x','dx','y','dy')
    g1 <- ggplot(raw, aes(x=x, y=y))+geom_point(alpha=alpha.data)+theme_gray(base_size=bsize)#+geom_errorbar(aes(ymin=y-dy,ymax=y+dy),width=0.2,alpha = alpha)+geom_errorbarh(aes(xmin=x-dx,xmax=x+dx),width=0.2,alpha = alpha)
    if(Nset==1) g1 <- g1+theme(legend.position = "none")

    if(FitType=='line'){
        g2 <- g1+geom_line(data=mp,mapping=aes(x=x,y=y), color = "red",size=1,alpha=alpha.fit)+xlab(xlab1)+ylab(ylab1)+ggtitle(title)
    }else{
        g2 <- g1+geom_point(data=mp,mapping=aes(x=x,y=y), color = "red",size=1,alpha=alpha.fit)+xlab(xlab1)+ylab(ylab1)+ggtitle(title)
    }
    if(!is.null(data.bin)){
        FigRaw <- g2+geom_point(data=data.bin,aes(x=x,y=y,colour=instrument),size=3)+coord_cartesian(ylim = range(mp[,2],data.bin[,2],mean(raw[,2])+2*sd(raw[,2]),mean(raw[,2])-2*sd(raw[,2])))+geom_errorbar(aes(ymin=y-dy,ymax=y+dy),width=0.2,alpha = alpha.data)+geom_errorbarh(aes(xmin=x-dx,xmax=x+dx),width=0.2,alpha = alpha.data)
    }else{
        FigRaw <- g2+coord_cartesian(ylim = range(mp[,4],raw[,4]),xlim=range(mp[,2],raw[,2]))
    }
###residual plot
    FigRes <- ggplot(res, aes(x, y))+geom_point(alpha =alpha.data)+xlab(xlab2)+ylab(ylab2)+geom_errorbar(aes(ymin=y-dy,ymax=y+dy),width=0.2,alpha = alpha.data)+geom_errorbarh(aes(xmin=x-dx,xmax=x+dx),width=0.2,alpha = alpha.data)+theme_gray(base_size=bsize)#+ annotate("text", x = 0.9*max(res[,1]), y = 0.9*max(res[,2]), label = paste0("RMS=",round(sd(res[,2]),2),' m/s'))
    if(!is.null(res.bin)){
        ylim <- range(res.bin[,2],mean(res[,2])+2*sd(res[,2]),mean(res[,2])-2*sd(res[,2]))
        FigRes <- FigRes+geom_point(data=res.bin,aes(x=x,y=y),size=3)+coord_cartesian(ylim = ylim)+ geom_errorbar(aes(ymin=y-dy,ymax=y+dy),width=0.2,alpha = alpha.data)+geom_errorbarh(aes(xmin=x-dx,xmax=x+dx),width=0.2,alpha = alpha.data)
    }
###combined plot
    PlotComb <- plot_grid(
        plot_grid(
            FigRaw+ theme(legend.position = "none",plot.title = element_text(hjust = 0.5)),
            FigRes+theme(legend.position = "none"),
            align='v',ncol=1,axis = "lr",rel_heights = c(1,0.5)
        ), plot_grid(
            get_legend(FigRaw),
             ggplot()+theme_minimal()+geom_blank(),
             ncol =1,rel_heights = c(1,0.5)
           )
, rel_widths = c(9.5,0.5)
    )
    plot(PlotComb)
#    plot(p)
}
