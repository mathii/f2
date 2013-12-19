## convert the AA map to macs format. Don't make a cut.map, because we're just using this
## to simulate, not to make inference. Be a bit careful with this one. It's
## just for testing. Assuming the aa map just has two columns, position and total length

args <- commandArgs(TRUE)
aa.map <- args[1]
macs.map <- args[2]

aa.map <- read.table(aa.map, header=TRUE, as.is=TRUE)
lo <- min(aa.map[,1])
hi <- max(aa.map[,1])

N <- NROW(aa.map)

from <- aa.map[1:(N-1),1]
to <- aa.map[2:N,1]

map.start <- aa.map[1:(N-1),2]
map.end <- aa.map[2:N,2]

rate <- (10^6)*(map.end-map.start)/(to-from)

from <- (from-lo)/(hi-lo)
to <- (to-lo)/(hi-lo)

new.map <- cbind(from, to, rate)

write.table(new.map, macs.map, sep="\t", row.names=FALSE, col.names=FALSE)
