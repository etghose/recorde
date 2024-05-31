# ReCorDE
Recurrent Correlation of Drugs with Enrichment

Reference: August J. John, Emily T. Ghose, Huanyao Gao, Meagan Luck, Dabin Jeong, Krishna R. Kalari, Liewei Wang. (2024) ReCorDE a framework for identifying drug classes targeting shared vulnerabilities with applications to synergistic drug discovery

ReCorDE identifies drug classes targeting shared vulnerabilities using large-scale, publicly available monotherapy cytotoxicity data from cancer cell lines; this information can aid efforts to identify novel synergistic drug combinations. ReCorDe identifies drug pairs with recurrently correlated response phenotypes, in multiple datasets; it then identifies drug class combinations over-represented in this set of recurrently correlated drug pairs. These enriched drug class combinations represent drug classes targeting shared vulnerabilities. 

This repository contains code used to implement ReCorDe as shown in John et. al 2024.

## File Descriptions
### Input
- **CleanMonotherapyDrugResponse.RDS** --- a list of data frames. Each list element is a data frame containing clean drug response data from a different publicly available drug response data set. Data frames are given as drug x cell line with cell values corresponding to the AUC.
- **NameStandardization.RDS** --- a list of data frames. Each list element is a data frame containing necessary drug name standardization information.
- **ATC_map.RDS** -- a data frame. Maps standardized drug names to ATC code drug classes.
- **ATC_dictionary.RDS** -- a data frame. Gives the definitions of each ATC code drug class; also contains information about whether the ATC code is canonical or not.

### Scripts
- **configReCorDE.R** --- configuration script. Always sourced by Correlations.R and Enrichment.R. Installs dependencies and initializes ReCorDE parameters. Basic ReCorDE parameters can easily be changed using this script (e.g. correlation coefficient cutoff, classes to exclude post-FDR, etc.). Please see configReCorDE.R comments for further details. Parameters are set as presented in John et al. 2024 by default. 
- **Correlations.R** -- correlation script. Runs the correlation step of ReCorDE. Takes CleanMonotherapyDrugResponse.RDS and NameStandardization.RDS as input. Outputs the unmerged, unpruned correlations results (UnprunedCorrelation.RDS) and the pruned correlation table (DCSCAD.RDS).
- **Enrichment.R** -- enrichment script. Runs the enrichment step of ReCorDE. Takes DCSCAD.RDS, ATC_map.RDS, and ATC_dictionary.RDS as input. 

### Other
- **sessionInfo.txt** -- session information from analysis run for John et al. 2024.

## Running ReCorDE
```r
source("Correlations.R"")
source("Enrichment.R")
```
