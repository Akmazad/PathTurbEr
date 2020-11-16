# This script is inspired by following website
# http://www.gettinggeneticsdone.com/2012/03/pathway-analysis-for-high-throughput.html

# These are bioconductor packages. See http://www.bioconductor.org/install/ for installation instructions	

 
# significant genes is a vector of fold changes where the names
# are ENTREZ gene IDs. The background set is a vector of all the 
# genes represented on the platform.
library(dplyr)

source("DE_GSE38376_GEO2R_code.R")
# top <- DE_analysis()
top <- tT
top.filt <- subset(top, adj.P.Val<0.01) %>% 
  dplyr::select(c("Gene.ID","logFC")) %>%
  dplyr::filter(Gene.ID != "" & 
                  !base::grepl("///", Gene.ID, fixed = T) & 
                  !duplicated(Gene.ID)) %>%
  as.data.frame() 

sig_genes <- top.filt$logFC
names(sig_genes) <- top.filt$Gene.ID

res <- spia(de=sig_genes, all = top.filt$Gene.ID, plots = T, verbose = T)
#two-way evidence plot -- doesn't work
# res$pG = combfunc(res$pNDE, res$pPERT, combine = "norminv")
# res$pGFdr=p.adjust(res$pG,"fdr")
# res$pGFWER=p.adjust(res$pG,"bonferroni")
# plotP(res,threshold=0.05)

# # make custom SPIA data from only signaling pathways (manually downloaded)
# makeSPIAdata(kgml.path="C:\\Users\\Az\\Google Drive\\Aberrant Signaling Pathway Analyses\\New Experiment\\GSE16179\\Dataset\\KEGG_XML\\",organism="hsa",out.path="./")
# 
# # run SPIA.
# res<-spia(de=sig_genes, all=all_genes, organism="hsa",data.dir="./",plots=TRUE)
# 
# head(res)
# plotP(res, threshold=0.05)
# 
# # save the result
# write.csv(res[,-12],"C:\\Users\\Az\\Google Drive\\Aberrant Signaling Pathway Analyses\\New Experiment\\GSE16179\\SPIA_results_nP_GSE16179_45_new_2.csv")
