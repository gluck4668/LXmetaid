
library(openxlsx)

gene_data <- read.xlsx("gene_data.xlsx")
meta_data <- read.xlsx("meta_data.xlsx")

usethis::use_data(gene_data,overwrite = T)
usethis::use_data(meta_data,overwrite = T)

rm(list=ls())

data(gene_data)
data(meta_data)

