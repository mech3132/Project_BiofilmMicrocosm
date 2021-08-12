#!bin/bash/Rscript

library(tidyverse)
### Inhibitory Sampling ####
inhibDat <- read.delim("../../IsolateInfo/all_isolate_info_combined.txt")

inhibDat_filt <- inhibDat %>% filter(BiofilmFormer, !is.na(Inhibitory))

date_start <- c("2020-07-27")
njars <- 15
jarNames <- 1:15
RichLevels <- c(1,3,10)
totalJarsPerRep <- (length(RichLevels)*2+1)*2
nIsolatesPerRep <- max(RichLevels)
nreps <- njars%/%totalJarsPerRep

set.seed(280720201)
inhibDat_inhib <- inhibDat %>% filter(Inhibitory, BiofilmFormer) %>%
  mutate(randomOrder = sample(1:nrow(.), replace=FALSE)) %>%
  arrange(randomOrder) %>% 
  filter(randomOrder <= nIsolatesPerRep) %>%
  mutate(L0 = NA,L1 = c(rep(RichLevels[1], RichLevels[1]), rep(NA, n()-RichLevels[1])), L2 = c(rep(RichLevels[2], RichLevels[2]), rep(NA, n()-RichLevels[2])), L3=c(rep(RichLevels[3], RichLevels[3]), rep(NA, n()-RichLevels[3]))
         , Inhibitory="Inhib", Rep="Rep1", DateStart=date_start) 
# %>%
#   pivot_longer(cols=c(N1,N3,N10))
  # 
  # mutate(RichLevel = unlist(sapply(1:length(RichLevel), FUN=function(x) {rep(x=paste0("N",RichLevel)[x], times=RichLevel[x])}))) %>%
  # mutate(Inhibitory="Inhib", Rep="Rep1", DateStart=date_start) %>%
  # unite(Inhibitory,DateStart,RichLevel,Rep, col=SampleID, remove=FALSE)


set.seed(280720200)
inhibDat_notinhib <- inhibDat %>% filter(!Inhibitory, BiofilmFormer) %>%
  mutate(randomOrder = sample(1:nrow(.), replace=FALSE)) %>%
  arrange(randomOrder) %>% 
  filter(randomOrder <= nIsolatesPerRep) %>%
  mutate(L0= NA, L1 = c(rep(RichLevels[1], RichLevels[1]), rep(NA, n()-RichLevels[1])), L2 = c(rep(RichLevels[2], RichLevels[2]), rep(NA, n()-RichLevels[2])), L3=c(rep(RichLevels[3], RichLevels[3]), rep(NA, n()-RichLevels[3]))
         , Inhibitory="Non", Rep="Rep1", DateStart=date_start) 
# 
#   mutate(RichLevel = unlist(sapply(1:length(RichLevel), FUN=function(x) {rep(x=paste0("N",RichLevel)[x], times=RichLevel[x])}))) %>%
#   mutate(Inhibitory="Non", Rep="Rep1", DateStart=date_start) %>%
#   unite(Inhibitory,DateStart,RichLevel, Rep, col=SampleID, remove=FALSE)

## To print for chosen isolates
inhibDat_inhib %>% full_join(inhibDat_notinhib) %>%
  pivot_longer(cols=c(L1,L2,L3), names_to="RichNames", values_to="RichLevel") %>%
  arrange(IsolateID) %>%
  write.table(file="print_chosen_isolates_long.txt", row.names=FALSE, quote=FALSE, sep="\t")

inhibDat_inhib %>% full_join(inhibDat_notinhib) %>%
  arrange(IsolateID) %>%
  select(IsolateID, g_sp,Inhibitory, L1, L2, L3) %>%
  write.table(file="print_chosen_isolates_short.txt", row.names=FALSE, quote=FALSE, sep="\t")
  
## To print for downstream exp analyses
fullSet <- inhibDat_inhib %>% full_join(inhibDat_notinhib) %>%
  pivot_longer(cols=c(L0,L1,L2,L3), names_to="RichNames", values_to="RichLevel") %>% 
  filter(RichNames =="L0" | !is.na(RichLevel)) %>% 
  mutate(RichLevel = ifelse(is.na(RichLevel), 0, RichLevel)
         , IsolateID = ifelse(RichLevel==0, NA, IsolateID)
         , Inhibitory = ifelse(RichLevel==0, "Con",Inhibitory)) %>%
  select(IsolateID, Inhibitory, Rep, DateStart, RichLevel) %>%
  group_by(Inhibitory, Rep, DateStart,RichLevel) %>%
  mutate(subIsolateID=paste0("Isolate",seq(1,n()))) %>%
  ungroup() %>%
  spread(key=subIsolateID, value=IsolateID) %>%
  unite(Inhibitory,DateStart,RichLevel, Rep, col=SampleID, remove=FALSE) 
View(fullSet)

fullSet <- inhibDat_inhib %>% full_join(inhibDat_notinhib) %>%
  pivot_longer(cols=c(L0,L1,L2,L3), names_to="RichNames", values_to="RichLevel") %>% 
  filter(RichNames =="L0" | !is.na(RichLevel)) %>% 
  mutate(RichLevel = ifelse(is.na(RichLevel), 0, RichLevel)
         , IsolateID = ifelse(RichLevel==0, NA, IsolateID)
         , Inhibitory = ifelse(RichLevel==0, "Con",Inhibitory)) %>%
  select(IsolateID, Inhibitory, Rep, DateStart, RichLevel) %>%
  # group_by(Inhibitory, Rep, DateStart,RichLevel) %>%
  # mutate(subIsolateID=paste0("Isolate",seq(1,n()))) %>% 
  # ungroup() %>%
  mutate(Value=TRUE) %>% 
  # spread(key=IsolateID, value=Value) %>% View()
  pivot_wider(names_from = IsolateID, values_from=Value, values_fn={unique}) %>%
  select(-`NA`) %>%
  unite(Inhibitory,DateStart,RichLevel, Rep, col=SampleID, remove=FALSE) 

CV_set <- fullSet %>% mutate(Type="CV")
Bd_set <- fullSet %>% mutate(Type="Bd")
set.seed(280720203)
metadata <- CV_set %>% full_join(Bd_set) %>% ungroup() %>%
  mutate(JarID = sample(jarNames[1:n()], replace=FALSE)) %>%
  arrange(JarID) %>%
  select(SampleID, JarID, RichLevel, Rep, DateStart, everything()) %>% arrange(SampleID)
metadata %>%
  write.table(file="jar_metadata_rep1.txt", sep="\t", quote=FALSE, row.names=FALSE)

forCultureMixing <- metadata %>% arrange(JarID) %>%mutate(JarID = paste0("Jar_",JarID)) %>%
  select(-c(SampleID, RichLevel, Rep, DateStart, Inhibitory, Type)) %>%
  column_to_rownames(var="JarID") 

weights <- 1/rowSums(forCultureMixing,na.rm = TRUE)
weights[is.infinite(weights)] <- 0

forCultureMixing <- forCultureMixing * weights
forCultureMixing[is.na(forCultureMixing)] <- 0
forCultureMixing %>% rownames_to_column(var="JarID") %>%
  write.table(file="forCultureMixing.txt", sep="\t", quote=FALSE, row.names = FALSE)



