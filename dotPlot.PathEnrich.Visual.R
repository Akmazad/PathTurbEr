dotPlot.PathEnrich.Visual <- function(enrichmentResultFile = myEnrichmentDF, x = "GeneRatio", color = "pvalue",
                                      showCategory=20, size=NULL, split = NULL,
                                      font.size=20, title = "", orderBy="GeneRatio", decreasing=TRUE){
  library(enrichplot)
  library(ggplot2)
  library(DOSE)
  library(data.table)
  
  colorBy <- match.arg(color, c("pvalue", "p.adjust"))
  
  # ------------ This code-snippet can be used for further flexibilty
  if (x == "geneRatio" || x == "GeneRatio") {
    x <- "GeneRatio"
    if (is.null(size))
      size <- "Count"
  } else if (x == "count" || x == "Count") {
    x <- "Count"
    if (is.null(size))
      size <- "GeneRatio"
  } else if (is(x, "formula")) {
    x <- as.character(x)[2]
    if (is.null(size))
      size <- "Count"
  } else {
    ## message("invalid x, setting to 'GeneRatio' by default")
    ## x <- "GeneRatio"
    ## size <- "Count"
    if (is.null(size))
      size  <- "Count"
  }
  
  # Read the data
  df =  fread(enrichmentResultFile) %>% as.data.frame()
  # library(tidyr)
  # aa = lapply(df$Pathway.geneSymbol, FUN = function(x){
  #   strsplit(x, split = "/", fixed = T) %>% length() %>% return()
  # })
  # df <- dplyr::mutate(df, GeneRatio = eval(parse(text=GeneRatio)))
  
  # print(paste0(df$GeneRatio))
  
  idx <- order(df[[orderBy]], decreasing = decreasing)
  df$Description <- factor(df$Description, levels=rev(unique(df$Description[idx])))
  # plot the graph
  p <- ggplot(df, aes_string(x=x, y="Description", size=size, color=colorBy)) +
    geom_point() +
    scale_color_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE)) +
    # scale_color_gradientn(name = color, colors=sig_palette, guide=guide_colorbar(reverse=TRUE)) +
    ylab(NULL) + ggtitle(title) + theme_dose(font.size) + theme(legend.text=element_text(size=16)) + scale_size(range=c(3, 8))
  
  # 
  # ggsave(
  #   paste0("cnetPlot_", title, ".pdf"),
  #   p,
  #   path = "Data/PathEnrichResults/" ,
  #   device = "pdf",
  #   limitsize = F,
  #   width = 20,
  #   height = 20
  # )
  
  return(p)
}

p = dotPlot.PathEnrich.Visual(enrichmentResultFile = "data/pathEnrichResult_KEGG_45_signaling.csv")

p

ggsave(
  paste0("dotPlot_DE_PathEnrich.pdf"),
  p,
  path = "data/" ,
  device = "pdf",
  limitsize = F,
  width = 20,
  height = 20
)
