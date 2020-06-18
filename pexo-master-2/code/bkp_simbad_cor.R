if(!exists('star')) star <- commandArgs(trailingOnly=TRUE)
cmd <- paste('python get_astrometry_for_star.py',star)
py.path <- '/Users/ffeng/miniconda2/bin/python'
if(file.exists(py.path)) cmd <- paste(py.path,'get_cor_from_simbad.py',star)
cat(cmd,'\n')
system(cmd)
if(file.exists('test.txt')){
    astro <- read.table('test.txt',header=TRUE)
    ra <- astro[1,1]
    dec <- astro[1,2]
    plx <- astro[1,3]
    pmra <- astro[1,4]
    pmdec <- astro[1,5]
    rv <- astro[1,6]
    StarAstro <- c(ra,dec,pmra,pmdec,plx,rv)
    Par$epoch <- epoch <- 2457206.375
    cat('star=',star,'\n')
    cat('ra=',ra,'\n')
    cat('dec=',dec,'\n')
    cat('plx=',plx,'\n')
    cat('pmra=',pmra,'\n')
    cat('pmdec=',pmdec,'\n')
    cat('rv=',rv,'\n')
    system('rm test.txt')
}
