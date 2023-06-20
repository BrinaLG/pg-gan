# Project 1: Admixture and GANs

path = "\\Users\\brina\\Dropbox\\Brown University\\Research\\Projekt 1\\3 Data\\"
df <- read.table(paste(path,"igsr_samples.tsv",sep=""), sep = '\t', header = TRUE)

POP_list <- unique(df$Population.code)

#for POP un POP_list


for (POP in POP_list) {
  sample.names <- matrix(df[which(df[,"Population.code"]==POP),"Sample.name"])
  write.table(sample.names,paste(path,paste(POP,"_gan.txt",sep=""),sep=""),quote = FALSE, row.names = FALSE,col.names = FALSE)
}


