
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 2){
    stop("Two arguments required, source .Robj and target .txt", call.=FALSE)
}

source <- args[1]
target <- args[2]

load(source)
write.table(ParStat, file=target, quote=FALSE)
