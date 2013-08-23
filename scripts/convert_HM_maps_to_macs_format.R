## Really dumb, jsut convert the hapmap recombination maps to macs format

args <- commandArgs(TRUE)
hm.map <- args[1]
macs.map <- args[2]
cut.map <- args[3]

hm.map <- read.table(hm.map, header=TRUE)
lo <- min(hm.map$Position.bp.)
hi <- max(hm.map$Position.bp.)

pos <- (hm.map$Position.bp.-lo)/(hi-lo)
new.map <- cbind( pos[1:(length(pos)-1)], pos[2:length(pos)], hm.map$Rate.cM.Mb.[1:(length(pos)-1)] )
write.table(new.map, macs.map, sep="\t", row.names=FALSE, col.names=FALSE)

hm.map.cut <- hm.map
hm.map.cut[,2] <- hm.map[,2]-min(hm.map[,2])
write.table(hm.map.cut, cut.map, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
