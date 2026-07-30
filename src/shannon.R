# library(vegan)
# library(qiime2R)

# args = commandArgs(trailingOnly=TRUE)
# #abundance <- read.csv(args[1], sep=',')

# abundance <- read_qza( args[1] )
# abundance <- t(abundance$data)

# #patientabundance <- read.csv("DUBAI_Specie_heatmap_abundance.csv", row.names=1)
# output_path <- args[2]
# #abundance <- abundance[abundance$level != 'p',]
# #abundance <- abundance[abundance$level != 'o',]
# #abundance$level <- NULL
# #abundance$taxon <- NULL
# #head(abundance)
# gutdiversity2 <- diversity(abundance, index="shannon", MARGIN=2)
# write.csv(data.frame(shannon_diversity=gutdiversity2), output_path)

# png(filename=gsub(".csv", ".png", output_path), height=1000, width=1000)
# my_variable <- gutdiversity2
# # Layout to split the screen
# layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))
 
# # Draw the boxplot and the histogram 
# par(mar=c(0, 3.1, 1.1, 2.1))
# boxplot(my_variable , horizontal=TRUE , xaxt="n" , col=rgb(0.8,0.8,0,0.5) , frame=F)
# par(mar=c(4, 3.1, 1.1, 2.1))
# hist(my_variable , breaks=40 , col=rgb(0.2,0.8,0.5,0.5) , border=F , main="" , xlab="Shannon Diversity Index")
# dev.off()

library(vegan)

args <- commandArgs(trailingOnly=TRUE)
in_path <- args[1]
output_path <- args[2]

if (grepl("\\.qza$", in_path, ignore.case=TRUE)) {
  library(qiime2R)
  abundance <- read_qza(in_path)
  abundance <- t(abundance$data)
} else {
  # CSV no formato: taxon | sample1 | sample2 | ...
  df <- read.csv(in_path, check.names=FALSE)

  # remove 1a coluna (taxon)
  df <- df[, -1, drop=FALSE]

  abundance <- as.matrix(df)
  mode(abundance) <- "numeric"
  abundance[is.na(abundance)] <- 0
}

gutdiversity2 <- diversity(abundance, index="shannon", MARGIN=2)
write.csv(data.frame(shannon_diversity=gutdiversity2), output_path)

png(filename=gsub(".csv", ".png", output_path), height=1000, width=1000)
my_variable <- gutdiversity2

layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))

par(mar=c(0, 3.1, 1.1, 2.1))
boxplot(my_variable , horizontal=TRUE , xaxt="n" , col=rgb(0.8,0.8,0,0.5) , frame=F)

par(mar=c(4, 3.1, 1.1, 2.1))
hist(my_variable , breaks=40 , col=rgb(0.2,0.8,0.5,0.5) , border=F , main="" , xlab="Shannon Diversity Index")
dev.off()