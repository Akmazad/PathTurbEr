# PathTurbEr
Pathway Perturbation Driver identification 

## Acknowledgement
If you find PathTeurbEr as useful for your research, please cite our work by including the following citation:
- <b>Discovering novel cancer bio-markers in acquired lapatinib resistance using Bayesian methods.</b> <i>Briefings in Bioinformatics</i> [```Link to the paper```](https://doi.org/10.1093/bib/bbab137)
- Citation:
```
@article{10.1093/bib/bbab137,
    author = {Azad, A K M and Alyami, Salem A},
    title = "{Discovering novel cancer bio-markers in acquired lapatinib resistance using Bayesian methods}",
    journal = {Briefings in Bioinformatics},
    year = {2021},
    month = {04},
    issn = {1477-4054},
    doi = {10.1093/bib/bbab137},
    url = {https://doi.org/10.1093/bib/bbab137},
    note = {bbab137},
    eprint = {https://academic.oup.com/bib/advance-article-pdf/doi/10.1093/bib/bbab137/37074977/bbab137.pdf},
}
```


## Installation
----
- Visual Studio 
- JAGS: Can be downloaded from [```here```](https://sourceforge.net/projects/mcmc-jags). After installing JAGS, the binary file path (e.g. ```C:\Program Files\JAGS\JAGS-4.3.0\x64\bin```) should be used within the code to run. Moreover, the SharpJags.dll file should be added as a reference into the visual studio project.

## Parameter list
- Overall parameters: bnmcmc_method (i.e. NS, HAR, MH, or, All), and more
- BNMCMC parameters: maxParent, maxChild, nSamplingIter, nBurnIn
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

### MCMC sampling (a Visual Studio project)
### Optimal BN structure learning for a particular STP perturbation
Run [```BaseModule.cs```](https://github.com/Akmazad/PathTurbEr/blob/master/MCMC%20sampling/BaseModule.cs) for generating optimal STP perturbation BN from Neighbourhood sampler, Hit-and-Run sampler, and Metropolis-Hasting sampler. The files will be saved within the bin\Debug\BNMCMC output\ directory of the project under the default settings.

#### Statistical modeling of perturbation driver characterization
Run [```bugsSampling.cs```](https://github.com/Akmazad/PathTurbEr/blob/master/MCMC%20sampling/JAGS/bugsSampling.cs) for genetring ```alpha``` values each of the inferred BN structures. The files will be saved within the bin\Debug\JAGS_output\ directory of the project under the default settings.

# Structure of this repository [TBD]
- create a logfile when running
- create command-line interface for MCMC sampling
- call that from R using systems function
- Make an R package 
