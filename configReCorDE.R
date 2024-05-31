############################################################
###############Install Dependencies and Load Pkgs################
############################################################
#Install dependencies if need be
dependencies = c("plyr", "dplyr", "stringr", "pbapply")
dependencies = lapply(dependencies, function(x){
  if(eval(parse(text=paste("require(",x,")")))){
    op = NA
  } else{
    op = x
  }
  return(op)
}) 
dependencies = unlist(dependencies)
dependencies = na.omit(dependencies)
if(length(dependencies) > 0){
  install.packages(pkgs = dependencies)
}

#load libraries
library(stringr)
library(dplyr)
library(pbapply)

rm(dependencies)

############################################
###############Make Output Directory################
############################################
dir.create("output", showWarnings = FALSE)

############################################
###############Parameters################
############################################
######General Parameters
#seed value
seedval = 51

#######Correlation (DCSCAD) parameters
#minds_dcscad -- minimum number of datasets that drug combination must appear in to be tested during correlation
minds_dcsad = 2
#min.cl -- minimum number of cell lines that drugs must share to be tested during correlation
min.cl = 10
#cut_dcscad -- rho cutoff for the dcscad
cut_dcscad = 0.25
#operator_dcscad -- boolean operator used for dcscad cutoff (e.g. "<" would mean rho must be less than some cutoff)
operator_dcscad = ">=" 
#abstf_dcscad  -- whether we consider magnitude of rho or just rho when judging the cutoff
abstf_dcscad = T
#signcut_dcscad -- maximum cross database direction disagreement for dcscad. Given as fraction of minds_dscad.
signcut_dcscad = 0.5

#########Enrichment parameters
#postadjexcl -- strings used to exclude class combinations from final enrichment table after FDR correction. By default, "miscellaneous"-type classes 
postadjexcl = c("D01AE", "L01AX", "L01CX", "L01DC", "L01EX", "L01XX")
#adj.type -- type of multiple testing correction to use. Can be "bonferroni" or "BH". By default "BH"
adj.type = "BH"
