set.seed(1234)
library(dplyr)
library(data.table)
setwd("GSE59357 analysis/")
source("R/DE_GSE59357_GEO2R_code.R")
source("R/downstream_analysis.R")
source("R/plotFactory.R")

geo_accession = "GSE59357"

# --------- EXPLORATORY PERTURBATION ANALYSIS -------
# find DEGs ------
top = DE_analysis()

library(dplyr)
library(data.table)
library(ggplot2)
fc.dat = top %>% 
  # subset(., adj.P.Val<0.00001) %>% 
  dplyr::select(c("Gene.ID", "Gene.symbol","logFC", "P.Value", "adj.P.Val")) %>%
  dplyr::filter(Gene.ID != "" & 
                  !base::grepl("///", Gene.ID, fixed = T) & 
                  !duplicated(Gene.ID)) %>%
  as.data.frame() 

p <- volcanoPlot.Visual(fc.dat = fc.dat,
                        title = paste0("Differential Expression of ", geo_accession),
                        subtitle = "Dasatinib sensitive-VS-resitance conditions in pancreatic cancer cell-line")
ggsave(
  paste0("volcanoPlot_DE_limma.pdf"),
  p,
  path = "data/" ,
  device = "pdf",
  limitsize = F,
  width = 20,
  height = 20
)


# ---- downstream analysis ----
# pathway enrichment -----
top.filt <- top %>% 
  dplyr::filter(adj.P.Val < 0.01) %>% 
  dplyr::select(c("Gene.ID", "Gene.symbol","logFC")) %>%
  dplyr::filter(Gene.ID != "" & 
                  !base::grepl("///", Gene.ID, fixed = T) & 
                  !duplicated(Gene.ID)) %>%
  as.data.frame() 
givenSet = top.filt$Gene.symbol %>% unique()
filename = "data/KEGG_selected_pathways.csv"
outFilePath = paste0("data/pathEnrichResult_KEGG_45_signaling_rev.csv")
background_genes = top$Gene.symbol %>% unique()
res <- Pathway_geneSetEnrichmentAnalysis(givenSet = givenSet, 
                                  pop.filepath = filename, cutoff_col = "pvalue",
                                  background_genes = background_genes,
                                  outFilePath = outFilePath)

# enrichment visuals -----
dotPlot.PathEnrich.Visual(
  enrichmentResultFile = "data/pathEnrichResult_KEGG_45_signaling_rev.csv") %>%
ggsave(
  filename = paste0("dotPlot_DE_PathEnrich.pdf"),
  path = "data/" ,
  device = "pdf",
  limitsize = F,
  width = 20,
  height = 20
)
lollipopPlot.PathEnrich.Visual(enrich.df = res) %>%
  ggsave(filename = paste0("data/",
                           # "[",deType,"] ",
                           "lollipopPlot_plot",
                           # ora.popName[pop_i],
                           ".pdf"), device = cairo_pdf,
         limitsize = F, width = 9, height = 7)
cnt.plot <- cnetPlot.PathEnrich.Visual(enrich.df = res,
                           circular = F,
                           layout = "kk",
                           # layout = "kk",
                           fc.dat = top %>% 
                             dplyr::select(c("Gene.symbol","logFC"))
) 
cnt.plot %>% 
  ggsave(filename = paste0("data/", 
                           # "[",deType,"] ",
                           "cnetPlot_plot", 
                           # ora.popName[pop_i],
                           ".pdf"),
         limitsize = F, width = 12, height = 12)

# Pathway activation analysis [SPIA] -----
require(SPIA)
sig_genes <- top.filt$logFC
names(sig_genes) <- top.filt$Gene.ID

# with default KEGG pathways ----
res <- spia(de=sig_genes, all = top$Gene.ID, plots = T, verbose = T)
fwrite(res, file="SPIA_res_default.csv")
plotP(res, threshold = 0.05)
#two-way evidence plot -- doesn't work
res$pG = combfunc(res$pNDE, res$pPERT, combine = "norminv")
res$pGFdr=p.adjust(res$pG,"fdr")
res$pGFWER=p.adjust(res$pG,"bonferroni")
pdf(file = "data/two_way_evidence_plot.pdf")
plotP(res,threshold=0.05)
dev.off()
# with KEGG signalling pathways ----
makeSPIAdata(kgml.path="data/KEGG_XML",organism="hsa",out.path="./")
spia_res<-spia(de=sig_genes, all=top$Gene.ID, organism="hsa",data.dir="./",plots=TRUE)
spia_res %>% fwrite(file = paste0(outdir, "SPIA_res_signaling_only.csv"))
plotP(spia_res, threshold = 0.05)
#two-way evidence plot -- 
spia_res <- spia_res %>%
  dplyr::mutate(pG = combfunc(pNDE, pPERT, combine = "norminv"))
spia_res$pGFdr=p.adjust(spia_res$pG,"fdr")
spia_res$pGFWER=p.adjust(spia_res$pG,"bonferroni")
pdf(file = "data/two_way_evidence_plot.pdf")
plotP(spia_res,threshold=0.05)
dev.off()
# pathView plot of a specific pathway -------
require(pathview)
pathView.dat = top$logFC
names(pathView.dat) = top$Gene.symbol
pathview(gene.data = pathView.dat, 
         pathway.id = "04012", # Erbb pathway
         species = "hsa",
         kegg.native = T,
         gene.idtype = "SYMBOL",
         limit = list(-6,6)
)
pathview(gene.data = pathView.dat, 
          pathway.id = "04912", # Pancreatic cancer pathway
          species = "hsa",
          kegg.native = T,
          gene.idtype = "SYMBOL",
          limit = list(-6,6)
)
pathview(gene.data = pathView.dat, 
          pathway.id = "04912", # GnRH signaling pathway
          species = "hsa",
          kegg.native = T,
          gene.idtype = "SYMBOL",
          limit = list(-6,6)
)


# --------- CAUSAL PERTURBATION NETWORK MODELLING -------
causal_path_bn = "ErbB_signaling_pathway"
p <- bnmcmc_convergence.Plot(pathwayName = causal_path_bn)
p %>%
  ggsave(
    paste0(causal_path_bn, "_bnmcmc_convergence.pdf"),
    p,
    path = "data/" ,
    device = "pdf",
    limitsize = F,
    width = 6,
    height = 6
  )
p <- perturber_analysis(pathwayName = causal_path_bn)
pdf(file = "data/PerturBer_alpha_scores.pdf", height = 8, width = 6)
p
dev.off()
# do BNMCMC modeling using NS, HAR and MH algorithms for BN structure learning
# --------- STATISTICAL MODELLING PERTURBED NETWORK -------
# --------- PATHWAY PERTURBATION: VALIDATION -------