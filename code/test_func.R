foo <- function(a,b){
    x=a*b
    y=a+b
    list(x=x,y=y)
}

ff <- function(z){
    list2env(foo(1,2),envir=environment())
    cat('z+x=',z+x,'\n')
    cat('z+y=',z+y,'\n')
}
