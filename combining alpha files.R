library(data.table)
library(dplyr)

# ns = fread("../Project1/Project1/bin/Debug/TGF-BETA_SIGNALING_PATHWAY_NS_net.csv", header = F) %>% as.data.frame()
# har = fread("../Project1/Project1/bin/Debug/TGF-BETA_SIGNALING_PATHWAY_HAR_net.csv", header = F) %>% as.data.frame()
# mh = fread("../Project1/Project1/bin/Debug/TGF-BETA_SIGNALING_PATHWAY_MH_net.csv", header = F) %>% as.data.frame()
# 
# full = rbind(ns,har,mh)
# full.names = c(full$V1 %>% as.character(), full$V2 %>% as.character()) %>% unique()

fc.dat = fread("data/GSE38376_DE_genes_limma.csv")

ns.alpha = fread("../PathTurbEr/Preprocessing_MCMC_sampling/Preprocessing_MCMC_sampling/bin/Debug/JAGS_output/TGF-BETA_SIGNALING_PATHWAY_alpha_NS.csv") %>% 
  dplyr::select(c("geneID","alpha_mean")) %>%
  as.data.frame()
har.alpha = fread("../PathTurbEr/Preprocessing_MCMC_sampling/Preprocessing_MCMC_sampling/bin/Debug/JAGS_output/TGF-BETA_SIGNALING_PATHWAY_alpha_HAR.csv") %>% 
  dplyr::select(c("geneID","alpha_mean")) %>%
  as.data.frame()
mh.alpha = fread("../PathTurbEr/Preprocessing_MCMC_sampling/Preprocessing_MCMC_sampling/bin/Debug/JAGS_output/TGF-BETA_SIGNALING_PATHWAY_alpha_MH.csv") %>% 
  dplyr::select(c("geneID","alpha_mean")) %>%
  as.data.frame()

all.genes = data.frame(geneID = c(ns.alpha$geneID, har.alpha$geneID, mh.alpha$geneID) %>% unique())
all = dplyr::left_join(all.genes, ns.alpha, by="geneID") %>%
  dplyr::left_join(., har.alpha, by = "geneID") %>% 
  dplyr::left_join(., mh.alpha, by = "geneID") %>%
  dplyr::rename("From NS sampling" = "alpha_mean.x",
                "From HAR sampling" = "alpha_mean.y",
                "From MH sampling" = "alpha_mean")
all$LogFC = fc.dat[base::match(all$geneID, fc.dat$Gene.symbol), "logFC"] %>% unlist(use.names = F)
all$Gene.title = fc.dat[base::match(all$geneID, fc.dat$Gene.symbol), "Gene.title"] %>% unlist(use.names = F)

colnames(all)[5] = "Log Fold-Change"
all = all %>% dplyr::filter(`From NS sampling` > 0 | 
                              `From HAR sampling` > 0 |
                              `From MH sampling` > 0)
all$`From NS sampling` = ifelse(all$`From NS sampling` > 0, all$`From NS sampling`, "--")
all$`From HAR sampling` = ifelse(all$`From HAR sampling` > 0, all$`From HAR sampling`, "--")
all$`From MH sampling` = ifelse(all$`From MH sampling` > 0, all$`From MH sampling`, "--")
fwrite(all, file="data/Table1_rev.csv")

ns.alpha$shape = ifelse(ns.alpha$alpha_mean > 0, "trianle", "circle")
har.alpha$shape = ifelse(har.alpha$alpha_mean > 0, "trianle", "circle")
mh.alpha$shape = ifelse(mh.alpha$alpha_mean > 0, "trianle", "circle")


ns.alpha$logFC = fc.dat[base::match(ns.alpha$geneID, fc.dat$Gene.symbol), "logFC"]
har.alpha$logFC = fc.dat[base::match(har.alpha$geneID, fc.dat$Gene.symbol), "logFC"]
mh.alpha$logFC = fc.dat[base::match(mh.alpha$geneID, fc.dat$Gene.symbol), "logFC"]

fwrite(ns.alpha, file = "data/NS_cyto_node_attribute_rev.csv")
fwrite(har.alpha, file = "data/HAR_cyto_node_attribute_rev.csv")
fwrite(mh.alpha, file = "data/MH_cyto_node_attribute_rev.csv")

# Helper function to display Venn diagram: DataNovia (https://www.datanovia.com/en/blog/venn-diagram-with-r-or-rstudio-a-million-ways/)
# display_venn <- function(x, ...){
#   library(VennDiagram)
#   grid.newpage()
#   venn_object <- venn.diagram(x, filename = NULL, ...)
#   grid.draw(venn_object)
# }

# display_venn(
#   ,
#   category.names = c("NS" , "HAR", "MH"),
#   fill = c("#E69F00", "#56B4E9", "#009E73")
# )
