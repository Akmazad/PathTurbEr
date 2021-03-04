# PathTurbEr
Pathway Perturbation Driver identification 

## Installation
----
- Visual Studio 
- JAGS: Can be downloaded from [```here```](https://sourceforge.net/projects/mcmc-jags). After installing JAGS, the binary file path (e.g. ```C:\Program Files\JAGS\JAGS-4.3.0\x64\bin```) should be used within the code to run. Moreover, the SharpJags.dll file should be added as a reference into the visual studio project.

## Parameter list
- Overall parameters: bnmcmc_method (i.e. NS, HAR, MH, or, All), and more
- BNMCMC parameters: maxParent, maxChil, nSamplingIter, nBurnIn
- JAGS parameters: gamma_prior_a, gamma_prior_b, nSamplingIter, nBurnIn

## Running a demo
----
### Dataset
- Case data: [GSE38376] Laptinib-resistant gene expression of SKBR3 Breast cancer cell-line formatted as gene-level data, [```which can be found here```](https://github.com/Akmazad/PathTurbEr/blob/master/data/R_GE_data_GSE38376.csv).
- Control data: [GSE38376] Laptinib-sensitive gene expression of SKBR3 Breast cancer cell-line formatted as gene-level data, [```which can be found here```](https://github.com/Akmazad/PathTurbEr/blob/master/data/nR_GE_data_GSE38376.csv).
- Pathway data: Signalling pathways collected from KEGG database, [```which can be found here```](https://github.com/Akmazad/PathTurbEr/blob/master/data/KEGG_45_SIGNALING.csv)

### DE analysis
Run [```DE_GSE38376_GEO2R_code.R```](https://github.com/Akmazad/PathTurbEr/blob/master/DE_GSE38376_GEO2R_code.R) in R to conduct differential expression analysis. The output of this step should later be filtered based on desired threshold values of parameter, e.g. logFC,  p-value, adj.pval, etc. Note, this file independently collect the same gene expression data as above. 

### Pathway Enrichment of DEGs
Run [```Pathway_geneSetEnrichmentAnalysis.R```](https://github.com/Akmazad/PathTurbEr/blob/master/Pathway_geneSetEnrichmentAnalysis.R) for finding enrichment of the DEGs, found from previous step.

### MCMC sampling for optimal BN structure learning for a particular STP perturbation
Run [```BaseModule.cs```](https://github.com/Akmazad/PathTurbEr/blob/master/MCMC%20sampling/BaseModule.cs) within a visual studio project (e.g. MCMC sampling) for generating 


### MCMC sampling for statistical modeling of perturbation driver characterization

### Post-processing and analysis result

# Structure of this repository [TBD]
- create a logfile when running
- create command-line interface for MCMC sampling
- call that from R using systems function
- Make an R package 
