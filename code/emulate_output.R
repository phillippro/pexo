if(any(grepl(',',opt$var))){
    vars <- unlist(strsplit(opt$var,split=','))
}else{
    vars <- unlist(strsplit(opt$var,split=' '))
}
vars <- vars[vars!='']

OutAll <- c(OutObs,OutTime,JDutc=utc)
if(grepl('A',Par$component)){
    OutAll <- c(OutAll,OutAstro)
}
if(grepl('R',Par$component)){
    OutAll <- c(OutAll,OutRv)
}
###expand the shapiro list
tmp <- OutAll$ShapiroPlanet
names(tmp) <- paste0('Shapiro',names(tmp))
OutAll <- c(OutAll,tmp)
###expand the ephemerides list
tmp <- OutAll$Eph
names(tmp) <- paste0('Eph',names(tmp))
OutAll <- c(OutAll,tmp)
###expand the Roemer delay list
OutAll <- c(OutAll,OutAll$RoemerOrder)
###expand deflection angle list for solar system objects
if(grepl('A',Par$component)){
    tmp <- OutAll$SolarDefList
    names(tmp) <- paste0('Def',names(tmp))
    OutAll <- c(OutAll,tmp)
###expand deflection angle list for solar system objects
    tmp <- OutAll$SolarDefList
    names(tmp) <- paste0('Def',names(tmp))
    OutAll <- c(OutAll,tmp)
}
###expand the redshift list
if(grepl('R',Par$component)){
    OutAll <- c(OutAll,OutAll$Zcomb)
}
###
cn <- names(OutAll)
##get the name of list object
ns <- c()
out <- c()
for(var in vars){
    if(any(cn==var)){
        if(!is.null(dim(OutAll[[var]]))){
            ns <- c(ns,paste0(var,1:ncol(OutAll[[var]])))
        }else{
            ns <- c(ns,var)
        }
        out <- cbind(out,OutAll[[var]])
    }else{
        cat('Warning:',var,'cannot be found in the output variable list!\n')
    }
}
colnames(out) <- ns
cat('\nsave data to',opt$out,'\n')
write.table(out,file=opt$out,quote=FALSE,row.names=FALSE)

####plot
if(!file.exists('../input')) system('mkdir ../input')
if(grepl('T',Par$component)){
    source('timing_verbose.R')
    source('timing_plot.R')
}
if(grepl('A',Par$component)){
    source('astrometry_verbose.R')
    source('astrometry_plot.R')
}
if(grepl('R',Par$component)){
    source('rv_verbose.R')
    source('rv_plot.R')
}
