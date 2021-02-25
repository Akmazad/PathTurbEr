

get_d <- function(dist_matrix){
  d <- mean(c(apply(dist_matrix,1,min), apply(dist_matrix,2,min)))

  return(d)
}
getRandD <- function(a_degree, net){
  return(base::sample(which(degree(net) == a_degree),1, replace = T))
}
permuteTest <- function(net, t, d, d_td, N){
  r <- c()
  for (i in 1:N) {
    t_rand <- sapply(as.numeric(degree(net, v = t)), getRandD,net)
    d_rand <- sapply(as.numeric(degree(net, v = d)), getRandD,net)
    d_td_rand <- get_d(shortest.paths(net, v = t_rand %>% unique(), to=d_rand %>% unique()))
    r<-c(r,d_td_rand)
  }
  r[!is.finite(r)] <- NA
  m <- mean(r, na.rm = T)
  s <- sd(r, na.rm = T)
  z <- (d_td - m)/s
  return(z)
}


#' Title
#'
#' @param disease.genes 
#' @param driver.genes 
#'
#' @return
#' @export
#'
#' @examples
disease.proximity <- function(net, disease.genes, driver.genes){
  disease.genes = disease.genes[,1] %>% unlist(use.names = F)
  d <- disease.genes
  keep <- which(d %in% V(net)$name)
  d <- unique(d[keep])
  
  t <- driver.genes
  keep <- which(t %in% V(net)$name)
  t <- unique(t[keep])
  
  d_td <- get_d(shortest.paths(net, v = t, to=d))  
  
  z <- permuteTest(net, t, d, d_td, 1000)
  p <- pnorm(-abs(z))
  
  return(list(d_td = d_td, z = z, p = p))
}

# library(data.table)
# library(igraph)
# human.ppi <- fread("data/i2d.human.anno.ppi.Genes.csv") %>% as.data.frame() # can look for only the signaling net
# ppiNet <- graph_from_data_frame(human.ppi[,c("symbol1", "symbol2")], directed = F) %>% igraph::simplify()
# 
# all = fread("data/Table1_rev.csv")
# 
# ns.drivers = all %>% dplyr::filter(`From NS sampling` != '--') %>% dplyr::select(geneID) %>% unlist(use.names = F)
# har.drivers = all %>% dplyr::filter(`From HAR sampling` != '--') %>% dplyr::select(geneID) %>% unlist(use.names = F)
# mh.drivers = all %>% dplyr::filter(`From MH sampling` != '--') %>% dplyr::select(geneID) %>% unlist(use.names = F)
# 
# 
# # disease.proximity(net = ppiNet, disease.genes = fread("data/Breast Cancer Genes [PMID 32101536].txt", header = F, encoding = "UTF-8"),
# #                   driver.genes = ns.drivers)
# # # [1] 0.001274818
# # disease.proximity(net = ppiNet, disease.genes = fread("data/Breast Cancer Genes [PMID 32101536].txt", header = F, encoding = "UTF-8"),
# #                   driver.genes = har.drivers)
# # # [1] 5.734396e-06
# # disease.proximity(net = ppiNet, disease.genes = fread("data/Breast Cancer Genes [PMID 32101536].txt", header = F, encoding = "UTF-8"),
# #                   driver.genes = mh.drivers)
# # # [1] 4.035954e-06
# # disease.proximity(net = ppiNet, disease.genes = fread("data/Breast Cancer Genes [PMID 32101536].txt", header = F, encoding = "UTF-8"),
# #                   driver.genes = base::intersect(ns.drivers, har.drivers) %>% base::intersect(mh.drivers))
# # # [1] 2.769058e-06
# 
# disease.proximity(net = ppiNet, disease.genes = fread("data/Breast Cancer Genes [PMID 32101536].txt", header = F, encoding = "UTF-8"),
#                                      driver.genes = ns.drivers)
# # $d_td
# # [1] 1.207627
# #
# # $z
# # [1] -2.874407
# #
# # $p
# # [1] 0.002023933
# 
# disease.proximity(net = ppiNet, disease.genes = fread("data/Breast Cancer Genes [PMID 32101536].txt", header = F, encoding = "UTF-8"),
#                                        driver.genes = har.drivers)
# # $d_td
# # [1] 1.148305
# #
# # $z
# # [1] -4.58863
# #
# # $p
# # [1] 2.230821e-06
# 
# disease.proximity(net = ppiNet, disease.genes = fread("data/Breast Cancer Genes [PMID 32101536].txt", header = F, encoding = "UTF-8"),
#                                         driver.genes = mh.drivers)
# # $d_td
# # [1] 1.253333
# #
# # $z
# # [1] -4.228998
# #
# # $p
# # [1] 1.173672e-05
# 
# disease.proximity(net = ppiNet, disease.genes = fread("data/Breast Cancer Genes [PMID 32101536].txt", header = F, encoding = "UTF-8"),
#                                          driver.genes = base::intersect(ns.drivers, har.drivers) %>% base::intersect(mh.drivers))
# $d_td
# [1] 1.538462
#
# $z
# [1] -4.483789
#
# $p
# [1] 3.666455e-06

