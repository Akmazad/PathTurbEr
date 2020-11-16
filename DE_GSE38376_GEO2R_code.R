# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Fri Nov 13 21:20:09 EST 2020

################################################################
#   Differential expression analysis with limma
DE_analysis <- function(threshold = 0.01){
  library(Biobase)
  library(GEOquery)
  library(limma)
  
  # load series and platform data from GEO
  
  gset <- getGEO("GSE38376", GSEMatrix =TRUE, AnnotGPL=TRUE)
  if (length(gset) > 1) idx <- grep("GPL6947", attr(gset, "names")) else idx <- 1
  gset <- gset[[idx]]
  
  # make proper column names to match toptable 
  fvarLabels(gset) <- make.names(fvarLabels(gset))
  
  # group names for all samples
  gsms <- "111111111000000000"
  sml <- c()
  for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
  
  # log2 transform
  ex <- exprs(gset)
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }
  
  # set up the data and proceed with analysis
  sml <- paste("G", sml, sep="")    # set group names
  fl <- as.factor(sml)
  gset$description <- fl
  design <- model.matrix(~ description + 0, gset)
  colnames(design) <- levels(fl)
  fit <- lmFit(gset, design)
  cont.matrix <- makeContrasts(G1-G0, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2, 0.01)
  tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(fit2))
  
  tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol", "Gene.ID","Gene.title"))
  library(data.table)
  fwrite(tT, file="GSE38376_DE_genes_limma.csv")
  return(tT)
}

################################################################
# #   Boxplot for selected GEO samples
# library(Biobase)
# library(GEOquery)
# 
# # load series and platform data from GEO

# gset <- getGEO("GSE38376", GSEMatrix =TRUE, getGPL=FALSE)
# if (length(gset) > 1) idx <- grep("GPL6947", attr(gset, "names")) else idx <- 1
# gset <- gset[[idx]]
# 
# # group names for all samples in a series
# gsms <- "111111111000000000"
# sml <- c()
# for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
# sml <- paste("G", sml, sep="")  # set group names
# 
# # order samples by group
# ex <- exprs(gset)[ , order(sml)]
# sml <- sml[order(sml)]
# fl <- as.factor(sml)
# labels <- c("case_resistance","control_sensitive")
# 
# # set parameters and draw the plot
# palette(c("#f4dfdf","#dfeaf4", "#AABBCC"))
# dev.new(width=4+dim(gset)[[2]]/5, height=6)
# par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
# title <- paste ("GSE38376", '/', annotation(gset), " selected samples", sep ='')
# boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
# legend("topleft", labels, fill=palette(), bty="n")
