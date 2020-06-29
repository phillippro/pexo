 for(star in Par$stars){
        ind1 <- which(Data[,'star']==star)
        if(length(ind1)>0){
            tps <- unique(Data[ind1,'type'])
            for(tp in tps){
                ind2 <- which(Data[,'star']==star & Data[,'type']==tp)
                if(length(ind2)>0){
                    inss <- unique(Data[ind2,'instrument'])
                    for(ins in inss){
                        cat('star:',star,'type:',tp,';ins=',ins,'\n')
                        ind3 <- which(Data[,'star']==star & Data[,'type']==tp & Data[,'instrument']==ins)
                        if(length(ind3)>0){
                            if(tp=='rv'){
                                n <- paste0('bRv.',star,'.',ins)
                                if(any(names(ParIni)==n)) ParIni[n] <- ParIni[n]+mean(Data[ind3,2]-model[ind3,2])
                            }else{
                                for(coord in c('RA','DEC')){
                                    if(coord=='RA') colum <- 2
                                    if(coord=='DEC') colum <- 4
                                    if(tp=='abs'){
#                                        n <- paste0('b',toupper(coord),'.',star,'.',ins)
                                        n <- paste0(tolower(coord),'Off')
                                    }else{
                                        n <- paste0('b',toupper(coord),'.rel.',ins)
                                    }
                                    cat('n=',n,'\n')
                                    cat('any(names(ParIni)==n)',any(names(ParIni)==n),'\n')
                                    if(any(names(ParIni)==n)){
                                        if(tp=='abs'){
                                            res <- (Data[ind3,colum]-model[ind3,colum])*cos(Data[ind3,4])/DMAS2R
#                                            res <- (Data[ind3,colum]-model[ind3,colum])/DMAS2R
#                                            cat('res=',range(res),'mas\n')
                                            cat('n=',n,'\n')
                                            cat('value=',mean(res),'\n')
                                            ParIni[n] <- mean(res)
                                            n1 <- paste0('pm',tolower(coord),'Off')
                                            if(any(names(ParIni)==n1) & FALSE){
                                                res1 <- res*cos(Data[ind3,4])
                                                tmp <- lm(res1~Data[ind3,1])
                                                ParIni[n1] <- tmp$coefficients[2]*DJY
                                                cat(star,';',ins,';',n1,'=',ParIni[n1],'\n')
                                            }
                                        }else{
                                            ParIni[n] <- mean(Data[ind3,colum]-model[ind3,colum])
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
