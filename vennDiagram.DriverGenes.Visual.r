
#' This function draws venn diagramm of overlapping driver genes, derived from Neighbourhoood Sampler, Hit-and-Run sampler and Metropolish-Hasting sampler 
#'
#' @param x is dataframe, where the first, second and third columns belongs to alpha-values of drivers from Neighbourhoood Sampler, Hit-and-Run sampler and Metropolish-Hasting sampler
#'
#' @return a list; firt element is venn diagramm image, and the second one is setMap (with heatmap)
#' @export
#'
#' @examples
vennDiagram.DriverGenes.Visual <- function(x){
  library(purrr)
  library(RVenn)
  library(ggplot2)
  
  venn.dat = list("NS" = all %>% dplyr::filter(`From NS sampling` != '--') %>% dplyr::select(geneID) %>% unlist(use.names = F),
                  "HAR" = all %>% dplyr::filter(`From HAR sampling` != '--') %>% dplyr::select(geneID)  %>% unlist(use.names = F),
                  "MH" = all %>% dplyr::filter(`From MH sampling` != '--') %>% dplyr::select(geneID)  %>% unlist(use.names = F))
  
  # construct the venn object
  venn.obj = Venn(venn.dat)
  
  p <-list(ggvenn(venn.obj, slice = c(1, 2, 3)), 
           setmap(venn.obj))
  
  return(p)
}

vennDiagram.DriverGenes.Visual(x = all)[[1]]
