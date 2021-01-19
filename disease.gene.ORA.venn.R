
disease.gene.ORA.venn <- function(disease.genes, driver.genes, title){
  require(org.Hs.eg.db)
  dat <- list(disease.genes, driver.genes)
  ven.dat = Venn(dat)
  univ = org.Hs.egSYMBOL %>% as.character() # length(univ) = 61217
  er = enrichment_test(ven.dat, 1, 2, univ = univ)
  paste0(er$Significance)

  q <- qplot(er$Overlap_Counts, geom = "blank") +
    geom_histogram(fill = "lemonchiffon4", bins = 8, color = "black") +
    geom_vline(xintercept = length(overlap(ven.dat, c(1, 2))), color = "firebrick2",
               size = 2, linetype = "dashed", alpha = 0.7) +
    geom_text(aes(x = length(overlap(ven.dat, c(1, 2))), y = 1000, 
                  label = paste0("p-value:", er$Significance)), 
              color = "firebrick2", angle=90, vjust = -1) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(name = "Overlap Counts") +
    scale_y_continuous(name = "Frequency")
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
disease.gene.ORA.venn(disease.genes, driver.genes, title = "NS")

driver.genes = har.drivers %>% as.character()
disease.gene.ORA.venn(disease.genes, driver.genes, title = "HAR")

driver.genes = mh.drivers %>% as.character()
disease.gene.ORA.venn(disease.genes, driver.genes, title = "MH")

driver.genes = union(ns.drivers, har.drivers) %>% union(mh.drivers)  %>% unique() %>% as.character()
disease.gene.ORA.venn(disease.genes, driver.genes, title = "NS U HAR U MH")

driver.genes = intersect(ns.drivers, har.drivers) %>% intersect(mh.drivers)  %>% unique() %>% as.character()
disease.gene.ORA.venn(disease.genes, driver.genes, title = "NS and HAR and MH")


