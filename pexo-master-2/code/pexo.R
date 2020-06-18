#library("Rmpfr")
#library("magicaxis")
library(optparse)
library(orthopolynom)
library(pracma)
#options(digits= 16)
source('constants.R')
source('sofa_function.R')
source('astrometry_function.R')
source('general_function.R')
source('timing_function.R')
source('rv_function.R')

Tstart <- proc.time()
####Read parameter files
source('read_input.R')
if(grepl('T',Par$component)){
    OutBary <- time_Utc2tb(utc,Par)
    OutTime <- time_Ta2te(OutBary,Par)
}

if(grepl('A',Par$component)){
    OutAstroT <- astro_FullModel(OutBary,OutTime,Par,Mlens=Par$mC,component='T')
    OutAstroC <- astro_FullModel(OutBary,OutTime,Par,Mlens=Par$mT,component='C')
}

if(grepl('R',Par$component)){
        OutRv <- rv_FullModel(OutBary,OutTime,Par)
#        OutRv <- rv_Numerical(utc,Par,OutBary,OutTime)
}
cat('ok5\n')
if(Par$Unit=='TDB'){
    source('output_tcb2tdb.R')
}
cat('ok6\n')

if(Par$figure){
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
}

Tend <- proc.time()
cat('duration of the process: ',round((Tend-Tstart)[3],1),'second\n')

if(!is.null(opt$var)){
    source('save_output.R')
}
