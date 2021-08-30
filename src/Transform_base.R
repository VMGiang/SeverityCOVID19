library(data.table)
library(dplyr)
args <- commandArgs(trailingOnly = T)
base_dir <- args[1]
basepop <- args[2]

dat <- read.table(gzfile(paste0(base_dir,basepop,".QC.gz")),
                  header=T)
dat$BETA <- log(dat$OR)
write.table(dat, paste0(base_dir,basepop,".QC.Transformed"),
            quote=F, row.names=F)
