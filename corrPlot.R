#' Title
#'
#' @return
#' @export
#'
#' @examples
corrPlot <- function(){
  require(data.table)
  require(dplyr)
  require(igraph)
  require(infotheo)
  human.sigNet = fread("../Data/HuamnSignalingNet_WangLab.csv") %>% 
    dplyr::select(c(2,4)) %>% 
    as.data.frame() %>% 
    graph_from_data_frame(directed = F)
  
  filelist <- list.files(path = "C:\\Users\\Azad\\OneDrive - UNSW\\own files\\Fund application\\Salem\\Project-1\\Part 1\\Identifying key players in aberrant pathways\\Project1\\Project1\\bin\\Debug\\JAGS_output\\", pattern = "*.csv", full.names = T)
  
  res <- do.call(rbind, lapply(filelist, FUN = function(x){
    print(x %>% basename())
    alpha = fread(x) %>% 
      as.data.frame() %>%
      dplyr::select(c("geneID", "alpha_mean")) %>%
      dplyr::filter(alpha_mean > 0.0) %>% 
      as.data.frame()
    
    keep =  intersect (alpha$geneID, V(human.sigNet)$name) 
    deg = do.call("rbind", lapply(keep, FUN = function(y){
      # print(y)
      return(c(y, igraph::degree(human.sigNet, v = y, mode = "all", loops = F)))
    })) %>% as.data.frame()
    
    colnames(deg) = c("geneID", "deg")
    
    m = dplyr::inner_join(alpha, deg, by="geneID") %>% 
      dplyr::select(c(2,3)) %>%  
      as.data.frame() #%>% 
    # discretize()
    
    corr_val = cor(x = m$deg %>% as.numeric(), y=m$alpha_mean, method = "pearson") %>% abs() 
    return(c(x %>% basename(), corr_val))
    # mutinformation(m$deg , m$alpha_mean, method = "emp") %>% print()
  })) 
  colnames(res) <- c("PathwayName","Corr_deg_alpha")
  fwrite(res, file = "corr_deg_alpha_P53_KEGG_signalling.csv")
}

corrPlot()

library(ggplot2)
library(dplyr)
library(data.table)
library(plotly)
pathway = "TGF-BETA_SIGNALING_PATHWAY"
NS <- fread(paste0("../Project1/Project1/bin/Debug/BNMCMC output/LogFx_NS_", pathway, ".csv")) %>% as.data.frame()
HAR <- fread(paste0("../Project1/Project1/bin/Debug/BNMCMC output/LogFx_HAR_", pathway, ".csv")) %>% as.data.frame()
MH <- fread(paste0("../Project1/Project1/bin/Debug/BNMCMC output/LogFx_MH_", pathway, ".csv")) %>% as.data.frame()
plot.dat = cbind(1:nrow(NS), NS[,1],HAR[,1], MH[,1]) %>% as.data.frame()
colnames(plot.dat) = c("iteration", "NS", "HAR", "MH")

fig <- plot_ly(plot.dat, x = ~iteration, y = ~NS, name = 'NS', type = 'scatter', mode = 'lines') 
fig <- fig %>% add_trace(y = ~HAR, name = 'HAR', mode = 'lines') #%>%
  # add_segments(x = 3000,
  #              y = min(c(NS$V1,HAR$V1,MH$V1)),
  #              xend = 5000,
  #              yend = max(c(NS$V1,HAR$V1,MH$V1)),
  #              color = I("yellow"),
  #              showlegend = FALSE, inherit = T
  # )
fig <- fig %>% add_trace(y = ~MH, name = 'MH', mode = 'lines') %>%
  layout(title = "Optimal STP inference: Bayesian network sampling with MCMC methods",
         paper_bgcolor='rgb(255,255,255)', 
         xaxis = list(title = "Iteration", showline = T, zeroline = T, showgrid = F),
         yaxis = list (title = "Log Posteriori", showline = T, zeroline = T, showgrid = F)) %>%
  config(displayModeBar = F)
fig
# draw standardized expression in case study compared to control
std.dat <- fread("../Project1/Project1/bin/Debug/GSE38376_standardizedOnly.csv") %>% as.data.frame()

fig <- plot_ly(type = 'box') %>%
  add_boxplot(y = std.dat$V2, name = "Only Whiskers", boxpoints = T,
              marker = list(color = 'rgb(9,56,125)'),
              line = list(color = 'rgb(9,56,125)')) %>%
fig


# pathView plot
library(pathview)
pathView.dat = tT$logFC
names(pathView.dat) = tT$Gene.symbol
pathview(gene.data = pathView.dat, 
         pathway.id = "04350",
         species = "hsa",
         kegg.native = T,
         gene.idtype = "SYMBOL",
         limit = list(-6,6)
         )
#
library(dplyr)
source("DE_GSE38376_GEO2R_code.R")
# top <- DE_analysis()
top <- tT

NS_alpha = fread("data/P53 signaling/P53_SIGNALING_PATHWAY_alpha_NS.csv") %>% 
  dplyr::filter(alpha_mean > 0) %>% 
  # dplyr::left_join(., top, by=c("geneID" = "Gene.symbol")) %>%
  as.data.frame()
HAR_alpha = fread("data/P53 signaling/P53_SIGNALING_PATHWAY_alpha_HAR.csv") %>% 
  dplyr::filter(alpha_mean > 0) %>% 
  # dplyr::left_join(., top, by=c("geneID" = "Gene.symbol")) %>%
  as.data.frame()
MH_alpha = fread("data/P53 signaling/P53_SIGNALING_PATHWAY_alpha_MH.csv") %>%  
  dplyr::filter(alpha_mean > 0) %>% 
  # dplyr::left_join(., top, by=c("geneID" = "Gene.symbol")) %>%
  as.data.frame()

# cor.test(NS_alpha$alpha_mean, top[base::match(NS_alpha$geneID, top$Gene.symbol),"logFC"])
# cor.test(HAR_alpha$alpha_mean, top[base::match(HAR_alpha$geneID, top$Gene.symbol),"logFC"])
# cor.test(MH_alpha$alpha_mean, top[base::match(MH_alpha$geneID, top$Gene.symbol),"logFC"])
NS_alpha = dplyr::inner_join(NS_alpha, top, by=c("geneID" = "Gene.symbol"))
cor.test(NS_alpha$alpha_mean, top[base::match(NS_alpha$geneID, top$Gene.symbol),"logFC"])
cor.test(HAR_alpha$alpha_mean, top[base::match(HAR_alpha$geneID, top$Gene.symbol),"logFC"])
cor.test(MH_alpha$alpha_mean, top[base::match(MH_alpha$geneID, top$Gene.symbol),"logFC"])
