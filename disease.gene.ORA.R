#' Title
#'
#' @param disease.genes 
#' @param driver.genes 
#' @param population 
#' @param p.adjust.methods 
#'
#' @return
#' @export
#'
#' @examples
disease.gene.ORA <- function(disease.genes, driver.genes){
  disease.genes = disease.genes[,1] %>% unlist(use.names = F)
  M = disease.genes %>% length()
  K = driver.genes %>% length()
  n = intersect(disease.genes,driver.genes) %>% length()
  U = union(disease.genes,driver.genes) %>% length()
  res = phyper(x, M, N-M, K, lower.tail = FALSE)
  return(res)
}

library(data.table)
library(dplyr)

all = fread("data/Table1_rev.csv")

ns.drivers = all %>% dplyr::filter(`From NS sampling` != '--') %>% dplyr::select(geneID) %>% unlist(use.names = F)
har.drivers = all %>% dplyr::filter(`From HAR sampling` != '--') %>% dplyr::select(geneID) %>% unlist(use.names = F)
mh.drivers = all %>% dplyr::filter(`From MH sampling` != '--') %>% dplyr::select(geneID) %>% unlist(use.names = F)

disease.gene.ORA(disease.genes = fread("data/Breast Cancer Genes [PMID 32101536].txt", header = F, encoding = "UTF-8"),
                 driver.genes = ns.drivers)
disease.gene.ORA(disease.genes = fread("data/Breast Cancer Genes [PMID 32101536].txt", header = F, encoding = "UTF-8"),
                 driver.genes = har.drivers)
disease.gene.ORA(disease.genes = fread("data/Breast Cancer Genes [PMID 32101536].txt", header = F, encoding = "UTF-8"),
                 driver.genes = mh.drivers)
