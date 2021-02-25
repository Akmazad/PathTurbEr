#' Title
#'
#' @param driver.genes 
#' @param disease.genes 
#' @param pop.filepath 
#' @param combined 
#' @param preambleEdges 
#'
#' @return
#' @export
#'
#' @examples
sankey.plot.visual <- function(driver.genes, disease.genes, pop.filepath, combined=F, preambleEdges=NA){
  library(networkD3)
  require(stringi)
  require(dplyr)
  require(data.table)
  
  pathway.data = fread(pop.filepath, header = T) %>% dplyr::select(1,2)
  nTerms = nrow(pathway.data)
  
  disease.genes = setdiff(disease.genes, driver.genes)
  
  driverNodes = data.frame(name = driver.genes, group = "drivers")
  pathwayNodes = data.frame(name = pathway.data[,1] %>% 
                              unlist(use.names = F) %>% 
                              gsub("_", " ", ., fixed = T) %>% 
                              gsub(" PATHWAY","",., fixed = T), 
                            group = "pathways")
  diseaseNodes = data.frame(name = disease.genes, group = "disease")
  # nodes = data.frame(name = c(driver.genes, 
  #                             pathway.data[,1] %>% unlist(use.names = F), 
  #                             disease.genes) %>% unique())
  nodes = bind_rows(driverNodes, pathwayNodes, diseaseNodes)
  edges = data.frame(source = c(),
                     target = c(),
                     group = c())
  links = data.frame(source = c(),
                     target = c(),
                     group = c())
  
  if(combined){
    if(is.na(preambleEdges)){
      stop("for combined set, the preamble edge list shouldn't be empty!!")
    }
    nodes = data.frame(name = c("NS", "HAR", "MH"), group = "algorithms") %>% bind_rows(.,nodes)
    edges = preambleEdges
  }
  
  for (i in 1:nTerms) {
    print(i)
    path.name = pathway.data[i,1] %>% as.character() %>% gsub("_", " ", ., fixed = T) %>% gsub(" PATHWAY","",., fixed = T)
    pathway.genes.symbol = pathway.data[i,2] %>% as.character() %>% strsplit(.,split = "/",fixed = T) %>% unlist(use.names = F)
    dr.overlap = intersect(pathway.genes.symbol, driver.genes)
    if(length(dr.overlap) > 0){
      source = dr.overlap
      target = path.name
      group = "dr_pathway"
      edges = rbind(edges, cbind(source, target, group) %>% as.data.frame()) 
    }
    dis.overlap = intersect(pathway.genes.symbol, disease.genes)
    if(length(dis.overlap) > 0){
      source = path.name
      target = dis.overlap
      group = "pathway_dis"
      edges = rbind(edges, cbind(source, target, group) %>% as.data.frame()) 
    }
    
  }
  
  nodes = dplyr::inner_join(nodes, data.frame(name = c(edges[,1] %>% as.character(), 
                                                       edges[,2]  %>% as.character()) %>% unique()
                                              ), 
                            by = "name")
  
  links = data.frame(source = base::match(edges[,1] %>% as.character(), nodes[,1]) - 1,
                     target = base::match(edges[,2] %>% as.character(), nodes[,1]) - 1,
                     group = base::match(edges[,3] %>% as.character(), nodes[,1]) - 1)
  
  fwrite(edges, file = "shankey_edge_data.csv")
 
  sankeyNetwork(Links = links, Nodes = nodes, Source = 'source', LinkGroup = 'group', fontFamily = "Helvetica", 
                Target = 'target', Value = 1, NodeID = 'name', NodeGroup = "group", sinksRight = T,
                
                fontSize = 13, nodeWidth = 30, iterations = 0)
}

library(data.table)
library(dplyr)
all = fread("data/Table1_rev.csv")

ns.drivers = all %>% dplyr::filter(`From NS sampling` != '--') %>% dplyr::select(geneID) %>% unlist(use.names = F)
har.drivers = all %>% dplyr::filter(`From HAR sampling` != '--') %>% dplyr::select(geneID) %>% unlist(use.names = F)
mh.drivers = all %>% dplyr::filter(`From MH sampling` != '--') %>% dplyr::select(geneID) %>% unlist(use.names = F)
all.drivers = c(ns.drivers, har.drivers, mh.drivers) %>% unique()
common.drivers = intersect(ns.drivers, har.drivers) %>% intersect(mh.drivers)

filename = "data/pathEnrichResult_KEGG_45_signaling_rev.csv"

disease.genes = fread("data/Breast Cancer Genes [PMID 32101536].txt", header = F, encoding = "UTF-8")

preambleEdges = bind_rows(
  data.frame(source = "NS", target = ns.drivers, group = "ns_driver"),
  data.frame(source = "HAR", target = har.drivers, group = "har_driver"),
  data.frame(source = "MH", target = mh.drivers, group = "mh_driver")
)

sankey.plot.visual(driver.genes = all.drivers, 
                   disease.genes = disease.genes[,1] %>% unlist(use.names = F), 
                   pop.filepath = filename, combined = T, preambleEdges = preambleEdges)

sankey.plot.visual(driver.genes = mh.drivers, 
                   disease.genes = disease.genes[,1] %>% unlist(use.names = F), 
                   pop.filepath = filename)
