Pathway_geneSetEnrichmentAnalysis <- function(givenSet = givenSet, 
                                              pop.filepath = pop.filepath,
                                              background_genes,
                                              cutoff_col = "p.adjust",
                                              cutoff_th = 0.05,
                                              outFilePath = outFilePath){
  require(org.Hs.eg.db)
  require(stringi)
  require(AnnotationDbi)
  require(dplyr)
  require(data.table)
  require(clusterProfiler)
  
  
  num.col = max(count.fields(pop.filepath, 
                             sep = "\t", 
                             blank.lines.skip = TRUE), 
                na.rm = TRUE)
  pathway.data = read.table(pop.filepath, 
                            sep = "\t", header = FALSE, 
                            stringsAsFactors = FALSE, fill = TRUE, 
                            col.names = 1:num.col, 
                            blank.lines.skip = TRUE, quote = "")
  # newGO_terms = GO_terms[-which(rowSums(GO_terms != "", na.rm = T)-1 <= 15),]
  # newGO_terms = newGO_terms[-which(rowSums(newGO_terms != "", na.rm = T)-1 >= 100),]
  pathway.data = dplyr::select(pathway.data, -2) # throw away the the column with all NA (2nd-Column)
  
  # all.GO.terms = newGO_terms
  # row.names(pathway.data) = gsub(")","",lapply(lapply(stri_split_fixed(stri_reverse(pathway.data[,1]),"(",n = 2), FUN = stri_reverse), "[[",1))
  row.names(pathway.data) = pathway.data[, 1]
  nTerms = nrow(pathway.data)
  
  # pop.genes = unique(c(as.matrix(pathway.data[,-1])))
  # temp <- org.Hs.egSYMBOL
  # mapped_genes <- mappedkeys(temp)
  # temp <- as.list(temp[mapped_genes])
  
  # make the population set the union of background set and the pathway genes
  
  pop.genes <- union(background_genes,
                     pathway.data %>% dplyr::select(-c("X1")) %>% 
                       unlist(use.names = F) %>% unique())
  
  N = length(pop.genes)
  
  # enrichment test using HyperGeometric test
  info.gene.overrep = data.frame(matrix(Inf, nrow = nrow(pathway.data), ncol = 7))
  
  K = length(givenSet)
  
  
  for (i in 1:nTerms) {
    
    pathway.genes.symbol = as.character(pathway.data[i,which(pathway.data[i,] != '')])[-1]
    # pathway.genes.df = bitr(pathway.genes.symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    # pathway.genes.entrezID = as.character(pathway.genes.df$ENTREZID)
    
    M = length(pathway.genes.symbol)
    print(M)
    
    x.overlap.genes.symbol = intersect(pathway.genes.symbol, givenSet)
    
    print(i)
    x = length(x.overlap.genes.symbol)
    if(x > 2){
      info.gene.overrep[i,1] = pathway.data[i,1]
      info.gene.overrep[i,2] = paste0(pathway.genes.symbol, collapse = "/")
      info.gene.overrep[i,3] = paste0(x.overlap.genes.symbol, collapse = "/")
      
      x.overlap.genes.entrezID = bitr(x.overlap.genes.symbol, 
                                      fromType = "SYMBOL", 
                                      toType = "ENTREZID", 
                                      OrgDb = org.Hs.eg.db) %>%
        dplyr::select(ENTREZID) %>% unlist(use.names = F)
      
      
      info.gene.overrep[i,4] = paste0(x.overlap.genes.entrezID, collapse = "/")
      info.gene.overrep[i,5] = x
      info.gene.overrep[i,6] = paste0(x, "/", K)
      
      info.gene.overrep[i,7] = phyper(x, M, N-M, K, lower.tail = FALSE) #He wrote : overlap.genes-1
      #insert FDR val
      # info.gene.overrep[i,8] = p.adjust(info.gene.overrep[i,7], method = "fdr", n = nTerms)
    }
  }
  colnames(info.gene.overrep) = c("Description", 
                                  "Pathway.geneSymbol",
                                  "Overlapping.geneSymbol", 
                                  "Overlapping.geneID", 
                                  "Count", "GeneRatio", 
                                  "pvalue")
  info.gene.overrep = info.gene.overrep %>% 
    dplyr::filter(pvalue != Inf) %>% 
    dplyr::mutate("p.adjust" = p.adjust(pvalue, 
                                        method = "fdr", 
                                        n = nrow(.)))
  res <- info.gene.overrep %>%
    dplyr::filter(!!dplyr::sym(cutoff_col) < cutoff_th) %>%
    dplyr::arrange(!!dplyr::sym(cutoff_col))
     
  # fwrite(res, file = outFilePath, sep = ",", col.names = T, quote = "auto")
  return(res)
}

get_geneSymbol_to_ENTREZID <- function(geneSymbol_list = NULL) {
  if (is.null(geneSymbol_list)) {
    stop("No geneSymbol list provided")
  }
  
  require(org.Hs.eg.db)
  hs <- org.Hs.eg.db
  
  AnnotationDbi::select(hs,
                        keys = geneSymbol_list,
                        columns = c("ENTREZID", "SYMBOL"),
                        keytype = "SYMBOL") %>%
    dplyr::rename("ENTREZ" = ENTREZID) %>%
    dplyr::select(ENTREZ) %>% unlist(use.names = F) %>% return()
}

perturber_analysis <- function(pathwayName,
                               data.dir = paste0("MCMC sampling/bin/Debug/JAGS_output/")){
  require(data.table)
  require(tibble)
  
  alpha_MH <- fread(file = paste0(data.dir, pathwayName, "_alpha_MH.csv")) %>%
    dplyr::filter(alpha_mean > 0) %>% 
    dplyr::select(c(geneID, "alpha_MH" = alpha_mean))
  alpha_NS <- fread(file = paste0(data.dir, pathwayName, "_alpha_NS.csv")) %>%
    dplyr::filter(alpha_mean > 0) %>% 
    dplyr::select(c(geneID, "alpha_NS" = alpha_mean))
  alpha_HAR <- fread(file = paste0(data.dir, pathwayName, "_alpha_HAR.csv")) %>%
    dplyr::filter(alpha_mean > 0) %>% 
    dplyr::select(c(geneID, "alpha_HAR" = alpha_mean))
  
  plot.df <- alpha_MH %>% 
    dplyr::full_join(alpha_NS, by="geneID") %>%
    dplyr::full_join(alpha_HAR, by="geneID") %>%
    as_tibble() %>%
    tibble::column_to_rownames("geneID") %>%
    as.matrix()
  
  plot.df %>% 
    alpha_heatmap.Plot()
  
  
  
  
}
