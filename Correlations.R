###########################################################################
##########################------Run Configuration---------#########################
###########################################################################
print("Configuring")
source("configReCorDE.R")
set.seed(seedval)
###########################################################################
##########################------Load-Data---------#########################
###########################################################################
drug_resp = readRDS("input/CleanMonotherapyDrugResponse.RDS")
name_std = readRDS("input/NameStandardization.RDS")


###########################################################################
#################------Standardize names---------#########################
###########################################################################
drug_respW = lapply(1:length(drug_resp), function(i){
  drug.tmp = drug_resp[[i]]
  name.tmp = name_std[[i]]
  rownames(drug.tmp) = plyr::mapvalues(rownames(drug.tmp), from = name.tmp$IF.THIS, to = name.tmp$THEN.THIS, warn_missing = F) %>% toupper()
  rownames(drug.tmp) = gsub(pattern = "-", replacement = "_", x = rownames(drug.tmp))
  return(drug.tmp)
})
names(drug_respW) = names(drug_resp)
###########################################################################
#################------Correlation Prep---------#########################
###########################################################################
#identify drugs in 2+ datasets
all.drugs = lapply(drug_respW, rownames) %>% unlist()
all.drugs = table(all.drugs) %>% data.frame() %>% filter(Freq >= minds_dcsad)
all.drugs = all.drugs$all.drugs %>% as.character()

#subset each drug response set down to those drugs
drug_respW = lapply(drug_respW, function(df.tmp){
  op = df.tmp[rownames(df.tmp) %in% all.drugs,]
})

print("InitializingCorrelations")
#grab drug combinations tested in more than min.cl of the same cell lines
comb2test = lapply(drug_respW, function(df.tmp){
  #grab combo and drug response matrix, format if necessary
  comb.tmp = combn(x = rownames(df.tmp), m = 2) %>% data.frame()
  rstor = rownames(df.tmp)
  df.tmp = df.tmp %>% t() %>% data.frame()
  colnames(df.tmp) = rstor
  #grab combinations with > min.cl same cell lines tested
  tf = apply(comb.tmp, MARGIN = 2, function(tmp){
    n.overlap = df.tmp[,colnames(df.tmp) %in% tmp] %>% na.omit() %>% nrow()
    n.overlap[is.null(n.overlap)] = 0
    if(n.overlap >= min.cl){
      return(T)
    } else{return(F)}
  })
  op = comb.tmp[,tf]
  #sort
  op = apply(op, MARGIN = 2, sort) %>% data.frame()
})

#grab combinations tested in minds datasets
preprunecomb = lapply(comb2test, function(df.tmp){
  conc.tmp = paste0(df.tmp[1,], "-", df.tmp[2,])
}) %>% unlist()
preprunecomb = table(preprunecomb) %>% data.frame() %>% filter(Freq >= minds_dcsad)
preprunecomb = preprunecomb$preprunecomb %>% as.character()

#update comb2test
comb2test = lapply(comb2test, function(df.tmp){
  conc.tmp = paste0(df.tmp[1,], "-", df.tmp[2,])
  op = df.tmp[,conc.tmp %in% preprunecomb]
})

###########################################################################
#################------Run Correlations---------#########################
###########################################################################
cor.list = lapply(1:length(comb2test), function(i){
  #grab combos 2 test and the drug response matrices
  comb2test.tmp = comb2test[[i]]
  drug_res_tmp =  drug_respW[[i]]
  print(paste0("Performing correlations for dataset", " ", names(comb2test)[i]))
  #format 
  rstor = rownames(drug_res_tmp)
  drug_res_tmp = drug_res_tmp %>% t() %>% data.frame()
  colnames(drug_res_tmp) = rstor
  op = pbapply(comb2test.tmp, MARGIN = 2, function(comb.tmp){
    df.tmp = drug_res_tmp[,colnames(drug_res_tmp) %in% comb.tmp] %>% na.omit()
    ncl = nrow(df.tmp)
    d.a = comb.tmp[1]
    d.b = comb.tmp[2]
    set.seed(seedval)
    tst = cor.test(x = as.numeric(df.tmp[,d.a]), y = as.numeric(df.tmp[,d.b]), method = "spearman") %>% suppressWarnings()
    tst.op = c("DrugA" = d.a, "DrugB" = d.b, "DrugAB" = paste(d.a, d.b, sep = "-"), "coef" = tst$estimate %>% unname(), "p_raw" = tst$p.value, "ncl" = ncl)
  })
  op = op %>% t() %>% data.frame() %>% type.convert(as.is=T)
})
names(cor.list) = names(comb2test)

#save
saveRDS(object = cor.list, file = "output/UnprunedCorrelation.RDS")

###########################################################################
#################------Prune DCSCAD---------#########################
###########################################################################
#subset correlation results down to combinations meeting pruning criteria
cor.list.pruned = lapply(cor.list, function(df.tmp){
  df.tmpW = df.tmp
  #considers if we should look at absolute value or not for coef pruning
  if(abstf_dcscad == T){
    df.tmpW$coef = df.tmpW$coef %>% abs()
  }
  if(abstf_dcscad == F){
    df.tmpW$coef = df.tmpW$coef
  }
  #identify entries meeting coef pruning criteria
  tf = paste0(df.tmpW$coef, " ", operator_dcscad, " ", cut_dcscad)
  tf = sapply(tf, function(x){
    op = parse(text = x) %>% eval()
  })
  df.tmp = df.tmp[tf,]
  return(df.tmp)
})

#identify drug combinations still present in multiple datasets after pruning
postprunecomb = lapply(cor.list.pruned, function(x){return(x$DrugAB)}) %>% unlist()
postprunecomb = table(postprunecomb) %>% data.frame() %>% filter(Freq >= minds_dcsad)
postprunecomb = postprunecomb$postprunecomb %>% as.character()

cor.list.pruned = lapply(cor.list.pruned, function(df.tmp){
  df.tmp = subset(df.tmp, df.tmp$DrugAB %in% postprunecomb)
})

#shove everything together
cor.list.pruned.prep = lapply(1:length(cor.list.pruned), function(i){
  df.tmp = cor.list.pruned[[i]]
  df.tmp = df.tmp %>% select(-c(ncl)) %>% dplyr::rename(Coef = coef, P = p_raw)
  colnames(df.tmp)[colnames(df.tmp) %in% c("Coef", "P")] = paste0(colnames(df.tmp)[colnames(df.tmp) %in% c("Coef", "P")], "_", names(cor.list.pruned)[i])
  return(df.tmp)
})
names(cor.list.pruned.prep) = names(cor.list.pruned)

for (i in 1:length(cor.list.pruned.prep)){
  df.tmp = cor.list.pruned.prep[[i]]
  if(i == 1){
    #initialize df
    cor_pp = df.tmp %>% select(DrugA, DrugB, DrugAB)
  }
  cor_pp = full_join(cor_pp, df.tmp, by = c("DrugA", "DrugB", "DrugAB"))
  rm(df.tmp, i)
}

#identify and remove with > fraction unmatching signs
sign_tf = cor_pp[,grepl(pattern="Coef", x = colnames(cor_pp))] %>% sign()
sign_tf = apply(sign_tf, MARGIN =1, function(x){
  npos = length(x[x==1] %>% na.omit())
  nneg = length(x[x==-1] %>% na.omit())
  chk.tmp = c("pos" = npos, "neg" = nneg)
  frxn.agree = max(chk.tmp)/sum(chk.tmp)
  if(frxn.agree > signcut_dcscad){
    return(T)
  } else{
    return(F)
  }
})

cor_pp = cor_pp[sign_tf,]

saveRDS(object = cor_pp, file = "output/DCSCAD.RDS")

print("Correlations done")



