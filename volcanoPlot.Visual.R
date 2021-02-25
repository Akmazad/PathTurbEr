volcanoPlot.Visual <- function(fc.dat, title="VolcanoPlot", subtitle="Enhance Plot"){
  require(EnhancedVolcano)
  
  EnhancedVolcano(fc.dat,
                  lab = fc.dat$Gene.symbol,
                  x = 'logFC',
                  y = 'P.Value',
                  pCutoff = 10e-14,
                  title = title,
                  titleLabSize = 28,
                  subtitleLabSize = 24,
                  subtitle = subtitle,
                  FCcutoff = 0.5,
                  pointSize = 4.0,
                  labSize = 6.0,
                  labCol = 'black',
                  labFace = 'bold',
                  boxedLabels = TRUE,
                  colAlpha = 4/5,
                  legendPosition = 'right',
                  legendLabSize = 20,
                  legendIconSize = 4.0,
                  drawConnectors = TRUE, 
                  widthConnectors = 1.0,
                  colConnectors = 'black') %>% 
    return()
}
library(dplyr)
library(data.table)
library(ggplot2)
fc.dat = fread("data/GSE38376_DE_genes_limma.csv") %>% 
  # subset(., adj.P.Val<0.00001) %>% 
  dplyr::select(c("Gene.ID", "Gene.symbol","logFC", "P.Value")) %>%
  dplyr::filter(Gene.ID != "" & 
                  !base::grepl("///", Gene.ID, fixed = T) & 
                  !duplicated(Gene.ID)) %>%
  as.data.frame() 

p <- volcanoPlot.Visual(fc.dat = fc.dat,
                        title = "Differential Expression of GSE38376",
                        subtitle = "Lapatinib sensitive-VS-resitance conditions in breast cancer cell-line (SKBR3)")
p
ggsave(
  paste0("volcanoPlot_DE_limma.pdf"),
  p,
  path = "data/" ,
  device = "pdf",
  limitsize = F,
  width = 20,
  height = 20
)
