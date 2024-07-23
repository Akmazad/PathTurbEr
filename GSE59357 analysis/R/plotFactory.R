volcanoPlot.Visual <- function(fc.dat, title="VolcanoPlot", subtitle="Enhance Plot"){
  require(EnhancedVolcano)
  
  EnhancedVolcano(fc.dat,
                  lab = fc.dat$Gene.symbol,
                  x = 'logFC',
                  y = 'adj.P.Val',
                  pCutoff = 10e-2,
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
  df$Description = df$Description %>% 
    gsub(pattern = "_PATHWAY", replacement = "", fixed = T) %>% 
    gsub(pattern = "_", replacement = " ", fixed = T)
  
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

lollipopPlot.PathEnrich.Visual <- function(enrich.df = NULL){
  if(is.null(enrich.df)){
    stop("Enrichment data file must be provided!!!")
  }
  library(ggplot2)
  library(forcats)
  library(enrichplot)
  library(dplyr)
  library(hrbrthemes)
  library(extrafont)
  library(extrafontdb)
  
  
  enrich.df %>%
    dplyr::mutate(qscore = -log10(pvalue)) %>%
    dplyr::mutate(Description = Description %>% 
                    gsub(" Homo sapiens .*","",.))  %>% # for Reactome
    dplyr::mutate(Description = Description %>% 
                    gsub(" \\(GO:.*","",.))  %>%   # for GO terms
    dplyr::mutate(Description = Description %>%
                    gsub(" WP.*","",.))  %>% # for GO terms
    ggplot( 
      aes(qscore, fct_reorder(Description, qscore))) + 
    geom_segment(aes(xend=0, yend = Description)) +
    geom_point(aes(color=qscore, size = Count)) +
    scale_color_viridis_c(guide=guide_colorbar(title="-log10(adj.P)"), direction = -1) +
    # scale_color_viridis_c(option = "B", limits = c(0, 1)) +
    # scale_size_continuous(range=c(2, 10)) +
    theme_minimal() +
    xlab(paste0("-log10(pvalue)")) +
    ylab(NULL) + 
    theme(axis.line = element_line()) +
    # theme_ipsum() +
    ggtitle("Enriched Pathways")
}
cnetPlot.PathEnrich.Visual <-
  function(enrich.df,
           showCategory = 100,
           fc.dat   = NULL,
           layout = "kk",
           colorEdge = F,
           circular = T,
           analysis = "common",
           node_label = "all") {
    require(enrichplot)
    require(ggplot2)
    require(data.table)
    require(qdapTools)
    require(ggplot2)
    require(cowplot)
    require(tidyr)
    require(dplyr)
    require(stringr)
    require(DOSE)
    require(ggraph)
    require(igraph)
    require(hrbrthemes)
    
    
    # throw away homo sapiens part, GO, and Wiki part
    enrich.df = tidyr::separate(enrich.df,
                                1,
                                into = c("Description", "Junk"),
                                sep = " Homo sapiens "
    )[, -2] # for Reactome
    enrich.df = tidyr::separate(enrich.df,
                                1,
                                into = c("Description", "Junk"),
                                sep = " \\(GO:"
    )[, -2]  # for GO terms
    enrich.df = tidyr::separate(enrich.df,
                                1,
                                into = c("Description", "Junk"),
                                sep = " WP"
    )[, -2]  # for WikiPathways
    
    n <- update_n(enrich.df, showCategory)
    geneSets <- extract_geneSets(enrich.df, n)
    # print(geneSets)
    g <- list2graph(geneSets)
    
    size <- sapply(geneSets, length)
    V(g)$size <- min(size) / 2
    
    n <- length(geneSets)
    n <- length(V(g)) - n
    V(g)$size[(n + 1):length(V(g))] <- size
    # print(paste0(V(g)$name[(n + 1):length(V(g))]))
    # print(V(g)$name)
    
    node_label <-
      match.arg(node_label, c("category", "gene", "all", "none"))
    if (circular) {
      # geom_edge <- geom_edge_bend
      # geom_edge <- geom_edge_bend
      geom_edge <- geom_edge_bundle_minimal
    }else{
      geom_edge <- geom_edge_bundle_minimal
    }
    
    if (colorEdge) {
      E(g)$category <- rep(names(geneSets), sapply(geneSets, length))
      edge_layer <- geom_edge(aes_(color = ~ category), alpha = .3)
    }else {
      edge_layer <- geom_edge(aes(alpha = stat(index)),
                              strength=0.5,
                              colour = 'darkgrey',
                              show.legend=F)
    }
    
    # load Foldchange
    if (!is.null(fc.dat)) {
      colnames(fc.dat) = c("SYMBOL", "FC")
      # print(head(fc.dat))
      # V(g)$color <- "#009999"
      V(g)$color <- NA
      temp = dplyr::inner_join(data.frame(SYMBOL = V(g)$name), fc.dat, by =
                                 "SYMBOL")
      # palette <- fc_palette(min(fc.dat$FC, na.rm = T):max(fc.dat$FC, na.rm=T))
      palette <- fc_palette(-6:6)
      # print(head(V(g)$name))
      print(min(fc.dat$FC, na.rm = T))
      print(max(fc.dat$FC, na.rm=T))
      
      # tkplot(g)
      V(g)$color[1:n] <- temp$FC
      
      # degree
      V(g)$degree = degree(g)
      
      
      ggraph(g, layout = layout, circular = circular) +
        edge_layer +
        geom_node_point(aes_(
          fill = ~as.numeric(color),
          size = ~size), color = "black", shape = 21) +
        scale_fill_gradient2(
          name = "Log2(Fold Change)",
          low = scales::muted("green"),
          mid = "white",
          high = scales::muted("red"),
          midpoint = 0,
          space = "Lab",
          na.value = "#e8ebb2",
          guide = "colourbar"
          # aesthetics = "colour"
        ) +
        scale_size(range = c(3, 10), breaks = unique(round(seq(
          min(size), max(size), length.out = 4)))) +
        theme_void() + 
        geom_node_text(aes_(label = ~name), repel = TRUE)
    }
  }



update_n <- function(df, showCategory) {
  if (!is.numeric(showCategory)) {
    return(showCategory)
  }
  n <- showCategory
  if (nrow(df) < n) {
    n <- nrow(df)
  }
  
  return(n)
}

extract_geneSets <- function(df, n) {
  n <- update_n(df, n)
  df <- df[1:n, ]
  geneSets <- str_split(df[, 3], "/")
  names(geneSets) = df[, 1]
  return(geneSets) ## if n is a vector of Description
}

list2graph <- function(inputList) {
  x <- list2df(inputList)
  # print(x)
  g <- graph_from_data_frame(x, directed = FALSE)
  return(g)
}

fc_readable <- function(x, foldChange = NULL) {
  if (is.null(foldChange))
    return(NULL)
  
  
  gid <- names(foldChange)
  if (is(x, 'gseaResult')) {
    ii <- gid %in% names(x@geneList)
  } else {
    ii <- gid %in% x@gene
  }
  gid[ii] <- x@gene2Symbol[gid[ii]]
  names(foldChange) <- gid
  return(foldChange)
}

fc_palette <- function(fc) {
  if (all(fc > 0, na.rm = TRUE)) {
    palette <- enrichplot::color_palette(c("gray", "red"))
  } else if (all(fc < 0, na.rm = TRUE)) {
    palette <- enrichplot::color_palette(c("darkgreen", "gray"))
  } else {
    palette <-
      enrichplot::color_palette(c("darkgreen", "#0AFF34", "#B3B3B3", "#FF6347", "red"))
    # enrichplot::color_palette(c("darkgreen", "gray", "red"))
    
  }
  return(palette)
}

bnmcmc_convergence.Plot <- function(pathwayName, data.dir = paste0("MCMC sampling/bin/Debug/BNMCMC_output/")){
  # read data
  set.seed(123)
  library(data.table)
  library(viridis)
  library(ggplot2)

  
  df = fread(file = paste0(data.dir, "LogFx_MH_", pathwayName, ".csv"))
  iterations <- df$V2
  
  log_fx_MH <- df %>%  dplyr::select(V1) %>% unlist(use.names = F)
  log_fx_NS <- fread(file = paste0(data.dir, "LogFx_NS_", pathwayName, ".csv")) %>%
    dplyr::select(V1) %>% unlist(use.names = F)
  log_fx_HAR <- fread(file = paste0(data.dir, "LogFx_HAR_", pathwayName, ".csv")) %>%
    dplyr::select(V1) %>% unlist(use.names = F)
  
  plot.df <- data.frame(
    iteration = rep(iterations, 3),
    log_fx = c(log_fx_MH, log_fx_NS, log_fx_HAR),
    algorithm = rep(c("Metropolis-Hastings", "Neighbourhood Sampler", "Hit-and-Run Sampler"), each = length(iterations))
  )
  
  # Create the line plot
  # plot.df %>%
  #   ggplot(aes(x = iteration, y = log_fx, color = algorithm, linetype = algorithm)) +
  #   geom_line() +
  #   scale_color_viridis_c() +
  #   labs(title = "Convergence of Log-Likelihood Over Iterations",
  #        x = "Iterations",
  #        y = "Log-Likelihood",
  #        color = "Algorithm",
  #        linetype = "Algorithm") +
  #   theme_minimal() +
  #   theme(
  #     plot.title = element_text(hjust = 0.5),
  #     leglegend.position = "bottom"
  #   )
  # 
  # # faceted line plot
  # ggplot(plot.df, aes(x = iteration, y = log_fx, color = algorithm, linetype = algorithm)) +
  #   geom_line() +
  #   facet_wrap(~ algorithm, scales = "free_y") +
  #   scale_color_viridis_d() +
  #   labs(title = "Convergence of Log-Likelihood Over Iterations",
  #        x = "Iterations",
  #        y = "Log-Likelihood",
  #        color = "Algorithm",
  #        linetype = "Algorithm") +
  #   theme_minimal() +
  #   theme(
  #     plot.title = element_text(hjust = 0.5),
  #     legend.position = "bottom"
  #   )
  
  # Smoothed line plot
  ggplot(plot.df, aes(x = iteration, y = log_fx, color = algorithm, linetype = algorithm)) +
    geom_line(alpha = 0.3) +  # Adding transparency to the original lines
    geom_smooth(se = FALSE, method = "loess", span = 0.1) +
    scale_color_viridis(discrete = T, option = "H") +
    labs(title = "Smoothed Convergence of Log-Likelihood Over Iterations",
         x = "Iterations",
         y = "Log-Likelihood",
         color = "Algorithm",
         linetype = "Algorithm") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
  
}

alpha_heatmap.Plot <- function(plot.df){
  require(dplyr)
  require(ComplexHeatmap)
  require(viridis)
  require(imputeTS)
  require(circlize)

  plot.df <- plot.df %>%
    na_replace(.,0) 
  plot.df %>%
    ComplexHeatmap::Heatmap(name = "Perturbation Score (alpha)", 
                            cluster_rows = T,
                            cluster_columns = T,
                            show_column_dend = F,
                            border = T,
                            cell_fun = function(j, i, x, y, width, height, fill) {
                              grid.text(sprintf("%.3f", plot.df[i, j]), x, y, gp = gpar(fontsize = 7))
                            },
                            # border_gp = gpar(col = "black"),
                            # col = viridis(100, direction = -1),
                            col = colorRamp2(c(0, base::max(plot.df, na.rm = T)), 
                                             c("white", "indianred1")),
                            show_row_names = T,
                            show_column_names = T, row_title = NULL, 
                            show_row_dend = F, row_names_gp = gpar(fontsize = 6))
  
}

