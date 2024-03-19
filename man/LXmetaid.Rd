\name{LXmetaid}
\alias{LXmetaid}
\title{LXmetaid}
\usage{
LXmetaid(gene_data,meta_data,species)
}
\description{ LXmetaid can be used to convert HMDB ID to KEGG ID.

}
\examples{

#----代谢物ID转换-----------

install.packages("devtools")
library(devtools)

install_github("gluck4668/LXmetaid")

library(LXmetaid)

??LXmetaid
#---------------------
data(gene_data)
data(meta_data)
#--------------------

rm(list=ls())

gene_data <- "gene_data.xlsx"
meta_data <- "meta_data.xlsx"

species <- "rat" # "human", "rat", or "mouse"

LXmetaid(gene_data,meta_data,species)

}
