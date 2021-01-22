
disease.gene.ORA.venn <- function(disease.genes, driver.genes, title, hyperPvalue = ""){
  require(org.Hs.eg.db)
  require(RVenn)
  require(ggplot2)
  dat <- list(disease.genes, driver.genes)
  ven.dat = Venn(dat)
  univ = org.Hs.egSYMBOL %>% as.character() # length(univ) = 61217
  er = enrichment_test(ven.dat, 1, 2, univ = univ)
  paste0(er$Significance)

  q <- qplot(er$Overlap_Counts, geom = "blank") +
    geom_histogram(fill = "orange", bins = 8, color = "black", alpha=0.2) +
    geom_vline(xintercept = length(overlap(ven.dat, c(1, 2))), color = "firebrick2",
               size = 2, linetype = "dashed", alpha = 0.7) +
    geom_text(aes(x = length(overlap(ven.dat, c(1, 2)))-0.05, y = 3000, 
                  label = paste0("p-value (empirical): ", er$Significance)), 
              color = "firebrick2", size = 5, angle=90, vjust = -1) +
    geom_text(aes(x = length(overlap(ven.dat, c(1, 2)))+0.1, y = 4400, 
                  label = paste0("p-value (hypergeometric test): ", hyperPvalue)), 
              color = "blue", size = 5, angle=90, vjust = 1) + 
    ggtitle(title) +
    scale_x_continuous(name = "Overlap Counts") +
    scale_y_continuous(name = "Frequency") + theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20, face = "bold"))
    # theme_classic(plot.title = element_text(hjust = 0.5), text = element_text(size = 20, face = "bold"))
  return(q)
}

library(data.table)
library(dplyr)

all = fread("data/Table1_rev.csv")

ns.drivers = all %>% dplyr::filter(`From NS sampling` != '--') %>% dplyr::select(geneID) %>% unlist(use.names = F)
har.drivers = all %>% dplyr::filter(`From HAR sampling` != '--') %>% dplyr::select(geneID) %>% unlist(use.names = F)
mh.drivers = all %>% dplyr::filter(`From MH sampling` != '--') %>% dplyr::select(geneID) %>% unlist(use.names = F)

disease.genes = fread("data/Breast Cancer Genes [PMID 32101536].txt", header = F, encoding = "UTF-8") %>% unlist(use.names = F)
driver.genes = ns.drivers %>% as.character()
ns_plot <- disease.gene.ORA.venn(disease.genes, driver.genes, title = "NS", hyperPvalue = "1.42e-19")
ggsave("data/ns_plot_ORA.pdf", ns_plot)

driver.genes = har.drivers %>% as.character()
har_plot <- disease.gene.ORA.venn(disease.genes, driver.genes, title = "HAR", hyperPvalue = "1.42e-19")
ggsave("data/har_plot_ORA.pdf", har_plot)


driver.genes = mh.drivers %>% as.character()
mh_plot <- disease.gene.ORA.venn(disease.genes, driver.genes, title = "MH", hyperPvalue = "6.95e-22")
ggsave("data/mh_plot_ORA.pdf", mh_plot)




# driver.genes = union(ns.drivers, har.drivers) %>% union(mh.drivers)  %>% unique() %>% as.character()
# disease.gene.ORA.venn(disease.genes, driver.genes, title = "NS U HAR U MH")
# 
# driver.genes = intersect(ns.drivers, har.drivers) %>% intersect(mh.drivers)  %>% unique() %>% as.character()
# disease.gene.ORA.venn(disease.genes, driver.genes, title = "NS and HAR and MH")


