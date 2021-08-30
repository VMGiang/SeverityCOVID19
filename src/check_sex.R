# Read in file
args <- commandArgs(trailingOnly = T)
tar_dir <- args[1]
pop_tar <- args[2]

setwd(tar_dir)

valid <- read.table(paste(pop_tar,".valid.sample",sep=''), header=T)
dat <- read.table(paste(pop_tar,".QC.sexcheck",sep=''), header=T)
valid <- subset(dat, STATUS=="OK" & FID %in% valid$FID)
write.table(valid[,c("FID", "IID")], paste(pop_tar,".QC.valid",sep =""), row.names=F, col.names=F, sep="\t", quote=F) 
