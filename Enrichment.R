###########################################################################
##########################------Run Configuration---------#########################
###########################################################################
source("configReCorDE.R")
set.seed(seedval)

##################################################################################################
################################Load Data########################################
##################################################################################################
ATC_map = readRDS("input/ATC_map.RDS")
ATC_dictionary = readRDS("input/ATC_dictionary.RDS")

UnprunedCorrelation <- readRDS("output/UnprunedCorrelation.RDS")
DCSCAD = readRDS("output/DCSCAD.RDS")
##################################################################################################
################################Functions########################################
##################################################################################################

#function takes two columns ("DrugA" and "DrugB") from a df. Columns DrugA and DrugB should be drugs
#maps DrugA and DrugB to classes, and then put them together in an alphabetical combo (returns a vector of class combinations)
#map is a df that has a column NAME (containing drugs) and CLASS (with corresponding ATC code)
drug2class = function(df, map){
  cpdA = df$DrugA
  cpdB = df$DrugB
  cpdA2class = plyr::mapvalues(cpdA, from=map$NAME, to=map$CLASS, warn_missing = F)
  cpdB2class = plyr::mapvalues(cpdB, from=map$NAME, to=map$CLASS, warn_missing = F)
  #standardize combos (alphabetical order before pasting together)
  df.tmp = cbind.data.frame(cpdA2class, cpdB2class)
  output = apply(df.tmp, MARGIN=1, function(x){
    vec.sort = sort(c(x[1], x[2]))
    class.comb = paste(vec.sort[1], vec.sort[2], sep="-")
    return(class.comb)
  })
  return(output)
}
#enrichment function
#given a target and universe frequency table, computes enrichment scores for all 
#entries in a frequency table
enrich = function(target, universe, correction = NULL){
  #target is df with two columns based off annotations from signif hits: combo and freq
  #universe is same format but for all possibilities in universe
  output = apply(target, MARGIN=1, FUN = function(x){
    combo = x[1]
    #N is total number of combinations (annotations) considered (combination universe)
    #n is number of signif combos (annotations)
    #m is number of combos (annotations) in universe that match the combo of interest we are testing for enrichment
    #k is number of combo of interest in the signif subset
    N = sum(universe[,"freq"])
    n = sum(target[,"freq"])
    m = universe[universe$ClassComb == combo, "freq"]
    k = as.numeric(x[2])
    #initialize fisher table
    fisher_table = data.frame(matrix(data=NA, nrow=2, ncol=2))
    rownames(fisher_table) = c('term_of_interest', "not_term_of_interest")
    colnames(fisher_table) = c("signif", "not_signif")
    #fill in fisher table
    fisher_table[1,1] = k
    fisher_table[1,2] = m-k
    fisher_table[2,1] = n-k
    fisher_table[2,2] = N + k - n - m
    if (0 %in% colSums(fisher_table)){print(c(combo, "colsum"))}
    if (0 %in% rowSums(fisher_table)){print(c(combo, "rowsum"))} 
    #do Hypergeometric test
    enrichment.test = fisher.test(fisher_table, alternative="greater")
    #will return vector of p-value, odds ratio estimate, lower bound CI, upper bound CI
    return(c(combo,enrichment.test[["p.value"]], enrichment.test[["estimate"]][["odds ratio"]], enrichment.test[["conf.int"]][1], enrichment.test[["conf.int"]][2]))
  })
  output = as.data.frame(t(output))
  colnames(output) = c("ClassCombination", "p_value", "odds_ratio", "CI_lower_bound", "CI_upper_bound")
  output$p_value = as.numeric(output$p_value)
  output$odds_ratio = as.numeric(output$odds_ratio)
  output$CI_lower_bound = as.numeric(output$CI_lower_bound)
  output$CI_upper_bound = as.numeric(output$CI_upper_bound)
  if(is.null(correction) == T){
    return(output)
  }
  if(correction == "BH"){
    padj = p.adjust(p = output$p_value, method="BH")
    output = cbind.data.frame(output, padj)
    output$padj = as.numeric(output$padj)
    return(output)
  }
  if(correction == "bonferroni"){
    padj = p.adjust(p = output$p_value, method="bonferroni")
    output = cbind.data.frame(output, padj)
    output$padj = as.numeric(output$padj)
    return(output)
  }
}


##################################################################################################
################################Background Dist########################################
##################################################################################################
background_dist = lapply(UnprunedCorrelation, function(df.tmp){
  op = drug2class(df.tmp[,c("DrugA", "DrugB")], map = ATC_map)
}) %>% unlist() %>% table() %>% data.frame() 
colnames(background_dist) = c("ClassComb", "freq")
background_dist = background_dist %>% mutate(ClassComb = as.character(ClassComb))

##################################################################################################
################################Target Dist########################################
##################################################################################################
dcscad.drugclass = drug2class(DCSCAD[,c("DrugA", "DrugB")], map = ATC_map)
dcscad.weights = apply(DCSCAD, MARGIN = 1, function(x){
  druglevelW = x[4:length(x)]
  druglevelW = druglevelW[!(is.na(druglevelW))] %>% length()
  druglevelW = druglevelW/2
})
target_dist = lapply(1:length(dcscad.drugclass), function(i){
  vec = rep(dcscad.drugclass[i], dcscad.weights[i])
}) %>% unlist() %>% table() %>% data.frame()
colnames(target_dist) = c("ClassComb", "freq")
target_dist = target_dist %>% mutate(ClassComb = as.character(ClassComb))


##################################################################################################
################################Enrich########################################
##################################################################################################
enrichment.df = enrich(target = target_dist, universe = background_dist, correction = paste0(adj.type))

saveRDS(enrichment.df, file = "output/RawEnrichment.RDS")
##################################################################################################
################################PostProcessing########################################
##################################################################################################
#exclude like and readjust
sameclass = paste0(ATC_map$CLASS, "-", ATC_map$CLASS)
enrichment.df = subset(enrichment.df, !(enrichment.df$ClassCombination %in% sameclass))
enrichment.df$padj = p.adjust(p = enrichment.df$p_value, method = paste0(adj.type))

#exclude other classes from final res
enrichment.df = subset(enrichment.df, !grepl(pattern = paste0(postadjexcl, collapse = "|"), x = enrichment.df$ClassCombination))
enrichment.df = subset(enrichment.df, enrichment.df$padj < 0.05)
saveRDS(enrichment.df, file = "output/FinalResults.RDS")
