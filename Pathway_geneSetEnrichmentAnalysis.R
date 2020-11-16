Pathway_geneSetEnrichmentAnalysis <- function(givenSet = givenSet, pop.filepath = pop.filepath, outFilePath = outFilePath){
  require(org.Hs.eg.db)
  require(stringi)
  require(AnnotationDbi)
  require(dplyr)
  require(data.table)
  require(clusterProfiler)
  

  num.col = max(count.fields(pop.filepath, sep = "\t", blank.lines.skip = TRUE), na.rm = TRUE)
  pathway.data = read.table(pop.filepath, sep = "\t", header = FALSE, stringsAsFactors = FALSE, fill = TRUE, col.names = 1:num.col, blank.lines.skip = TRUE, quote = "")
  # newGO_terms = GO_terms[-which(rowSums(GO_terms != "", na.rm = T)-1 <= 15),]
  # newGO_terms = newGO_terms[-which(rowSums(newGO_terms != "", na.rm = T)-1 >= 100),]
  pathway.data = dplyr::select(pathway.data, -2) # throw away the the column with all NA (2nd-Column)
  
  # all.GO.terms = newGO_terms
  row.names(pathway.data) = gsub(")","",lapply(lapply(stri_split_fixed(stri_reverse(pathway.data[,1]),"(",n = 2), FUN = stri_reverse), "[[",1))
  nTerms = nrow(pathway.data)
  
  # pop.genes = unique(c(as.matrix(pathway.data[,-1])))
  temp <- org.Hs.egSYMBOL
  mapped_genes <- mappedkeys(temp)
  temp <- as.list(temp[mapped_genes])
  pop.genes <- unlist(temp, use.names = F)
  
  N = length(pop.genes)
  
  # enrichment test using HyperGeometric test
  info.gene.overrep = data.frame(matrix(Inf, nrow = nrow(pathway.data), ncol = 8))
  
  K = length(givenSet)
  
  
  for (i in 1:nTerms) {
    
    pathway.genes.symbol = as.character(pathway.data[i,which(pathway.data[i,] != '')])[-1]
    # pathway.genes.df = bitr(pathway.genes.symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    # pathway.genes.entrezID = as.character(pathway.genes.df$ENTREZID)
    
    M = length(pathway.genes.symbol)
    
    x.overlap.genes.symbol = intersect(pathway.genes.symbol, givenSet)
    
    print(i)
    x = length(x.overlap.genes.symbol)
    info.gene.overrep[i,1] = pathway.data[i,1]
    if(x > 0){
      info.gene.overrep[i,2] = paste0(pathway.genes.symbol, collapse = "/")
      info.gene.overrep[i,3] = paste0(x.overlap.genes.symbol, collapse = "/")
      
      x.overlap.genes.df = bitr(x.overlap.genes.symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
      x.overlap.genes.entrezID = as.character(x.overlap.genes.df$ENTREZID)
      
      info.gene.overrep[i,4] = paste0(x.overlap.genes.entrezID, collapse = "/")
      info.gene.overrep[i,5] = x
      info.gene.overrep[i,6] = paste0(x, "/", K)
      
      info.gene.overrep[i,7] = phyper(x, M, N-M, K, lower.tail = FALSE) #He wrote : overlap.genes-1
      #insert FDR val
      info.gene.overrep[i,8] = p.adjust(info.gene.overrep[i,7], method = "fdr", n = nTerms)
    }
  }
  colnames(info.gene.overrep) = c("Description", "Pathway.geneSymbol","Overlapping.geneSymbol", "Overlapping.geneID", "Count", "GeneRatio", "pvalue", "p.adjust")
  res <- info.gene.overrep[which(info.gene.overrep$pvalue != Inf),]
  fwrite(res, file = outFilePath, sep = ",", col.names = T, quote = "auto")
  # return(res)
}

library(dplyr)

source("DE_GSE38376_GEO2R_code.R")
# top <- DE_analysis()
top <- tT
top.filt <- subset(top, adj.P.Val<0.00001) %>% 
  dplyr::select(c("Gene.ID", "Gene.symbol","logFC")) %>%
  dplyr::filter(Gene.ID != "" & 
                  !base::grepl("///", Gene.ID, fixed = T) & 
                  !duplicated(Gene.ID)) %>%
  as.data.frame() 

givenSet = top.filt$Gene.symbol %>% unique()
filename = "data/KEGG_45_SIGNALING.csv"
outFilePath = paste0("data/pathEnrichResult_KEGG_45_signaling.csv")
Pathway_geneSetEnrichmentAnalysis(givenSet = givenSet, pop.filepath = filename, outFilePath)
