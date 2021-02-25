#' Title
#'
#' @param disease.genes e.g. Breast cancer genes
#' @param driver.genes 
#' @param population 
#' @param p.adjust.methods 
#'
#' @return
#' @export
#'
#' @examples
disease.gene.ORA <- function(disease.genes, driver.genes){
  require(org.Hs.eg.db)
  disease.genes = disease.genes[,1] %>% unlist(use.names = F)
  M = disease.genes %>% length()
  K = driver.genes %>% length()
  x = intersect(disease.genes,driver.genes) %>% length()
  # N = union(disease.genes,driver.genes) %>% length()
  N = org.Hs.egSYMBOL %>% as.character() %>% length()
  
  res = phyper(x, M, N-M, K, lower.tail = FALSE)
  return(res)
}

# library(data.table)
# library(dplyr)
# 
# all = fread("data/Table1_rev.csv")
# 
# ns.drivers = all %>% dplyr::filter(`From NS sampling` != '--') %>% dplyr::select(geneID) %>% unlist(use.names = F)
# har.drivers = all %>% dplyr::filter(`From HAR sampling` != '--') %>% dplyr::select(geneID) %>% unlist(use.names = F)
# mh.drivers = all %>% dplyr::filter(`From MH sampling` != '--') %>% dplyr::select(geneID) %>% unlist(use.names = F)
# 
# disease.gene.ORA(disease.genes = fread("data/Breast Cancer Genes [PMID 32101536].txt", header = F, encoding = "UTF-8"),
#                  driver.genes = ns.drivers)
# disease.gene.ORA(disease.genes = fread("data/Breast Cancer Genes [PMID 32101536].txt", header = F, encoding = "UTF-8"),
#                  driver.genes = har.drivers)
# disease.gene.ORA(disease.genes = fread("data/Breast Cancer Genes [PMID 32101536].txt", header = F, encoding = "UTF-8"),
#                  driver.genes = mh.drivers)
# disease.gene.ORA(disease.genes = fread("data/Breast Cancer Genes [PMID 32101536].txt", header = F, encoding = "UTF-8"),
#                  driver.genes = base::intersect(ns.drivers, har.drivers) %>% base::intersect(mh.drivers))
# 
# # can try qplot for better visualization the enrichment test
# link: https://cran.r-project.org/web/packages/RVenn/vignettes/vignette.html
