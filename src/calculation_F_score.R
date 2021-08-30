args <- commandArgs(trailingOnly = T)
tar_dir <- args[1]
tar_pop <- args[2]

# tar_dir <- '/mnt/ThamHoang/gt_v1/data/cad/GSE90073/output/bfile/'
# tar_pop <- 'GSE90073'
setwd(tar_dir)

dat <- read.table(paste(tar_pop,".QC.het",sep=''), header=T) # Read in the target.het file, specify it has header
m <- mean(dat$F) # Calculate the mean  
s <- sd(dat$F) # Calculate the SD
valid <- subset(dat, F <= m+3*s & F >= m-3*s) # Get any samples with F coefficient within 3 SD of the population mean
write.table(valid[,c(1,2)], paste(tar_pop,".valid.sample",sep = ''), quote=F, row.names=F) # print FID and IID for valid samples
