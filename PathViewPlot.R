# pathView plot
library(pathview)
#
library(dplyr)
source("DE_GSE38376_GEO2R_code.R")
top <- DE_analysis()
tT <- top

pathView.dat = tT$logFC
names(pathView.dat) = tT$Gene.symbol
pathview(gene.data = pathView.dat, 
         pathway.id = "04350",
         species = "hsa",
         kegg.native = T,
         gene.idtype = "SYMBOL",
         limit = list(-6,6)
)
