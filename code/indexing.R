###This is to index the data of different types and instruments
Par$targets <- unique(Data[,'star'])
Par$types <- unique(Data[,'type'])
Par$instruments <- unique(Data[,'instrument'])
index <- list()
for(star in Par$targets){
    index[[star]] <- list()
    index[[star]]$all <- which(Data[,'star']==star)
    index[[star]]$types <- unique(Data[Data[,'star']==star,'type'])
    index.astro <- c()
    ins.astro <- c()
    for(type in Par$types){
        index[[star]][[type]] <- list()
        instr <- index[[star]][[type]]$instruments <- unique(Data[Data[,'star']==star & Data[,'type']==type,'instrument'])
        index[[star]][[type]]$all <- which(Data[,'star']==star & Data[,'type']==type)
        if(type=='abs' | type=='rel') index.astro <- c(index.astro,index[[star]][[type]]$all)
        if(type=='abs' | type=='rel') ins.astro <- c(ins.astro,index[[star]][[type]]$instruments)
        for(ins in instr){
            index[[star]][[type]][[ins]]$ind0 <- which(Data[,'star']==star & Data[,'instrument']==ins & Data[,'type']==type)
            index[[star]][[type]][[ins]]$ind1 <- which(Data[index[[star]][[type]]$all,'instrument']==ins)
        }
    }
    index[[star]]$astro <- list()
    index[[star]]$astro$all <- index.astro
    index[[star]]$astro$instruments <- ins.astro
    if(length(index.astro)>0){
        for(type in Par$types[Par$types=='abs' | Par$types=='rel']){
            instr <- unique(Data[Data[,'star']==star & Data[,'type']==type,'instrument'])
            index[[star]][[type]]$ind2 <- which(Data[index.astro,'type']==type)
            for(ins in instr){
                index[[star]][[type]][[ins]]$ind2 <- which(Data[index.astro,'instrument']==ins & Data[index.astro,'type']==type)
            }
        }
    }
}
Par$index <- index
