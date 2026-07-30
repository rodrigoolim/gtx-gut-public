#library(MicrobeR)
library(qiime2R)
#library(RColorBrewer)
#library(ggplot2)
#library(stringi)

args = commandArgs(trailingOnly=TRUE)
rawdata  <- args[1]
outtable <- args[2]

otus <- read_qza( rawdata )
toplot <- otus$data

write.csv( toplot, outtable ) ##### WRITE FULL DATA