# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Fri Nov 13 21:20:09 EST 2020

################################################################
#   Differential expression analysis with limma
DE_analysis <- function(threshold = 0.01,
                        outdir = "data/"){
  # Required libraries
  library(GEOquery)
  library(limma)
  library(data.table)
  library(dplyr)
  
  # Load series and platform data from GEO
  geo_accession = "GSE59357"
  gset <- getGEO(geo_accession, GSEMatrix = TRUE, AnnotGPL = TRUE)
  if (length(gset) > 1) idx <- grep("GPL10558", attr(gset, "names")) else idx <- 1
  gset <- gset[[idx]]
  
  # Make proper column names to match toptable 
  fvarLabels(gset) <- make.names(fvarLabels(gset))
  
  # Group membership for all samples
  gsms <- "000000000111111111"
  sml <- strsplit(gsms, split = "")[[1]]
  
  # Log2 transformation
  ex <- exprs(gset)
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
  LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
  if (LogC) { 
    ex[which(ex <= 0)] <- NaN
    exprs(gset) <- log2(ex) 
  }
  exprs(gset) <- normalizeBetweenArrays(exprs(gset)) # Normalize data
  
  # Assign samples to groups and set up design matrix
  gs <- factor(sml)
  groups <- make.names(c("sensitive", "resistant"))
  levels(gs) <- groups
  gset$group <- gs
  design <- model.matrix(~group + 0, gset)
  colnames(design) <- levels(gs)
  
  gset <- gset[complete.cases(exprs(gset)), ] # Skip missing values
  
  # split and case and control data matrix seperately
  extract_case_vs_control_data_matrix(gset = gset, outdir = outdir)
  
  # calculate precision weights and show plot of mean-variance trend
  v <- vooma(gset, design, plot=T)
  # OR weights by group
  # v <- voomaByGroup(gset, group=groups, design, plot=T, cex=0.1, pch=".", col=1:nlevels(gs))
  v$genes <- fData(gset) # attach gene annotations
  
  # fit linear model
  fit  <- lmFit(v)
  
  # set up contrasts of interest and recalculate model coefficients
  cts <- paste(groups[1], groups[2], sep="-")
  cont.matrix <- makeContrasts(contrasts=cts, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  
  # compute statistics and table of top significant genes
  fit2 <- eBayes(fit2, 0.01)
  tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(fit2))
  
  # tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","GB_ACC","Gene.symbol","Gene.title"))
  tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.ID","Gene.symbol", "Gene.title"))
  # write.table(tT, file=stdout(), row.names=F, sep="\t")
  
  
  fwrite(tT, file=paste0(outdir, "DE_genes_limma.csv"))
  return(tT)
}

extract_case_vs_control_data_matrix <- function(gset, outdir){
  require(data.table)
  require(tidyverse)
  require(Biobase)
  
  
  # Convert probe-level data to gene-level data
  annotation <- fData(gset)
  probe2gene <- annotation[, c("ID", "Gene.symbol")]
  probe2gene <- probe2gene[probe2gene$Gene.symbol != "", ]
  
  # Aggregate probe-level data by taking the mean for each gene
  ex_df <- as.data.frame(ex)
  ex_df$ProbeID <- rownames(ex_df)
  ex_df <- inner_join(ex_df, probe2gene, by = c("ProbeID" = "ID"))
  gene_ex <- ex_df %>%
    group_by(Gene.symbol) %>%
    summarize(across(starts_with("GSM"), mean, na.rm = TRUE))
  
  # Convert gene_ex back to a matrix with Gene.symbol as row names
  gene_ex <- as.data.frame(gene_ex)
  rownames(gene_ex) <- gene_ex$Gene.symbol
  gene_ex$Gene.symbol <- NULL
  gene_ex <- as.matrix(gene_ex)
  
  
  # Extract sample names for each group
  sensitive_samples <- colnames(exprs(gset))[gset$group == "sensitive"]
  resistant_samples <- colnames(exprs(gset))[gset$group == "resistant"]
  
  
  # Subset expression data matrix for each group
  sensitive_matrix <- gene_ex[, sensitive_samples, drop = FALSE]
  resistant_matrix <- gene_ex[, resistant_samples, drop = FALSE]
  
  
  
  # Save each subset matrix to a file
  sensitive_matrix %>% as.data.frame() %>%
    fwrite(file = paste0(outdir, "sensitive_group_matrix.csv"), 
         col.names = F,  row.names = TRUE)
  resistant_matrix %>% as.data.frame() %>%
    fwrite(file = paste0(outdir, "resistant_group_matrix.csv"), 
         col.names = F,  row.names = TRUE)
  
}