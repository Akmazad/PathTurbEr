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

### MCMC sampling for optimal BN structure learning for a particular STP perturbation

# Structure of this repository [TBD]
- create a logfile when running
- create command-line interface for MCMC sampling
- call that from R using systems function
- Make an R package 
