#' Title
#'
#' @param driver.genes 
#' @param disease.genes 
#' @param pop.filepath 
#'
#' @return
#' @export
#'
#' @examples
sankey.plot.visual <- function(driver.genes, disease.genes, pop.filepath){
  library(networkD3)
  # require(org.Hs.eg.db)
  require(stringi)
  # require(AnnotationDbi)
  require(dplyr)
  require(data.table)
  # require(clusterProfiler)
  
  
  num.col = max(count.fields(pop.filepath, sep = "\t", blank.lines.skip = TRUE), na.rm = TRUE)
  pathway.data = read.table(pop.filepath, sep = "\t", header = FALSE, stringsAsFactors = FALSE, fill = TRUE, col.names = 1:num.col, blank.lines.skip = TRUE, quote = "")
  # newGO_terms = GO_terms[-which(rowSums(GO_terms != "", na.rm = T)-1 <= 15),]
  # newGO_terms = newGO_terms[-which(rowSums(newGO_terms != "", na.rm = T)-1 >= 100),]
  pathway.data = dplyr::select(pathway.data, -2) # throw away the the column with all NA (2nd-Column)
  
  # all.GO.terms = newGO_terms
  row.names(pathway.data) = gsub(")","",lapply(lapply(stri_split_fixed(stri_reverse(pathway.data[,1]),"(",n = 2), FUN = stri_reverse), "[[",1))
  nTerms = nrow(pathway.data)
  
  nodes = data.frame(name = c(driver.genes, pathway.data[,1], disease.genes) %>% unique())
  edges = data.frame(source = c(),
                     target = c())
  links = data.frame(source = c(),
                     target = c())
  
  for (i in 1:nTerms) {
    print(i)
    path.name = pathway.data[i,1]
    pathway.genes.symbol = as.character(pathway.data[i,which(pathway.data[i,] != '')])[-1]
    dr.overlap = intersect(pathway.genes.symbol, driver.genes)
    if(length(dr.overlap) > 0){
      source = dr.overlap
      target = path.name
      edges = rbind(edges, cbind(source, target) %>% as.data.frame()) 
      
      # source = base::match(dr.overlap, nodes[,1]) - 1
      # target = base::match(path.name, nodes[,1]) - 1
      # links = rbind(links, cbind(source, target) %>% as.data.frame()) 
    }
    dis.overlap = intersect(pathway.genes.symbol, disease.genes)
    if(length(dis.overlap) > 0){
      source = path.name
      target = dis.overlap
      edges = rbind(edges, cbind(source, target) %>% as.data.frame()) 
      
      # source = base::match(path.name, nodes[,1]) - 1
      # target = base::match(dis.overlap, nodes[,1]) - 1
      # links = rbind(links, cbind(source, target) %>% as.data.frame()) 
    }
    
  }
  
  nodes = dplyr::inner_join(nodes, data.frame(name = c(edges[,1] %>% as.character(), edges[,2]  %>% as.character()) %>% unique()), by = "name")
  
  links = data.frame(source = base::match(edges[,1] %>% as.character(), nodes[,1]) - 1,
                     target = base::match(edges[,2] %>% as.character(), nodes[,1]) - 1)
  
  sankeyNetwork(Links = links, Nodes = nodes, Source = 'source',
                Target = 'target', Value = 1, NodeID = 'name',
                fontSize = 12, nodeWidth = 30)
}

library(data.table)
library(dplyr)
all = fread("data/Table1_rev.csv")

ns.drivers = all %>% dplyr::filter(`From NS sampling` != '--') %>% dplyr::select(geneID) %>% unlist(use.names = F)
har.drivers = all %>% dplyr::filter(`From HAR sampling` != '--') %>% dplyr::select(geneID) %>% unlist(use.names = F)
mh.drivers = all %>% dplyr::filter(`From MH sampling` != '--') %>% dplyr::select(geneID) %>% unlist(use.names = F)
all.drivers = c(ns.drivers, har.drivers, mh.drivers) %>% unique()

filename = "data/KEGG_45_SIGNALING.csv"

disease.genes = fread("data/Breast Cancer Genes [PMID 32101536].txt", header = F, encoding = "UTF-8")

sankey.plot.visual(driver.genes = all.drivers, disease.genes = disease.genes[,1] %>% unlist(use.names = F), pop.filepath = filename)