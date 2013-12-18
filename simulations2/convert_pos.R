## Read in a list and divide by a number - must be an easier way...

args <- commandArgs(TRUE)
data <- scan(args[1])
data <- data/as.numeric(args[2])
write.table(args[3], matrix(data, nrow=1), sep=" ", row.names=F, col.names=F)
