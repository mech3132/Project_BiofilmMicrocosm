#!bin/bash/Rscript

library(tidyverse)
setwd("./Rep12_31Aug2020/")
### Inhibitory Sampling ####
inhibDat <- read.delim("../../IsolateInfo/all_isolate_info_combined.txt")
# usedIso <- as.vector(read.delim("../Rep1b_25Aug2020/isolateListUSED.txt", header = FALSE)$V1)
## Usually, we'll use "remaining Isos in set" not used ones
# remainIso 

inhibDat_filt <- inhibDat %>% filter(BiofilmFormer, !is.na(Inhibitory))

date_start <- c("2020-08-31")
njars <- 30
jarNames <- 1:28
RichLevels <- c(1,3,10)
totalJarsPerRep <- (length(RichLevels)*2+1)*2
nIsolatesPerRep <- max(RichLevels)
nreps <- njars%/%totalJarsPerRep
# nIsolatesTotal <- nIsolatesPerRep*nreps
# 
# ## Check if we need to re-cycle some isolates
# inhibDat_removedlastRound <- inhibDat_filt %>% 
#   filter(!(IsolateID %in% usedIso))

##### INHIB #####
set.seed(9823472)
inhibDat_inhib1 <- inhibDat_filt %>% filter(Inhibitory, BiofilmFormer) %>%
  mutate(randomOrder = sample(1:nrow(.), replace=FALSE)) %>%
  arrange(randomOrder) %>% 
  filter(randomOrder <= nIsolatesPerRep) %>%
  mutate(randomOrder = 1:nrow(.)) %>%
  mutate(L0 = NA,L1 = c(rep(RichLevels[1], RichLevels[1]), rep(NA, n()-RichLevels[1])), L2 = c(rep(RichLevels[2], RichLevels[2]), rep(NA, n()-RichLevels[2])), L3=c(rep(RichLevels[3], RichLevels[3]), rep(NA, n()-RichLevels[3]))
         , Inhibitory="Inhib", Rep="Rep1", DateStart=date_start) 
inhibDat_inhib2 <- inhibDat_filt %>% filter(Inhibitory, BiofilmFormer) %>%
  mutate(randomOrder = sample(1:nrow(.), replace=FALSE)) %>%
  arrange(randomOrder) %>% 
  filter(randomOrder <= nIsolatesPerRep) %>%
  mutate(randomOrder = 1:nrow(.)) %>%
  mutate(L0 = NA,L1 = c(rep(RichLevels[1], RichLevels[1]), rep(NA, n()-RichLevels[1])), L2 = c(rep(RichLevels[2], RichLevels[2]), rep(NA, n()-RichLevels[2])), L3=c(rep(RichLevels[3], RichLevels[3]), rep(NA, n()-RichLevels[3]))
         , Inhibitory="Inhib", Rep="Rep2", DateStart=date_start) 
inhibDat_inhib <- rbind(inhibDat_inhib1, inhibDat_inhib2)

##### NOT INHIB #####
inhibDat_notinhib1 <- inhibDat_filt %>% filter(!Inhibitory, BiofilmFormer) %>%
  mutate(randomOrder = sample(1:nrow(.), replace=FALSE)) %>%
  arrange(randomOrder) %>% 
  filter(randomOrder <= nIsolatesPerRep) %>%
  mutate(randomOrder = 1:nrow(.)) %>%
  mutate(L0 = NA,L1 = c(rep(RichLevels[1], RichLevels[1]), rep(NA, n()-RichLevels[1])), L2 = c(rep(RichLevels[2], RichLevels[2]), rep(NA, n()-RichLevels[2])), L3=c(rep(RichLevels[3], RichLevels[3]), rep(NA, n()-RichLevels[3]))
         , Inhibitory="Non", Rep="Rep1", DateStart=date_start) 
inhibDat_notinhib2 <- inhibDat_filt %>% filter(!Inhibitory, BiofilmFormer) %>%
  mutate(randomOrder = sample(1:nrow(.), replace=FALSE)) %>%
  arrange(randomOrder) %>% 
  filter(randomOrder <= nIsolatesPerRep) %>%
  mutate(randomOrder = 1:nrow(.)) %>%
  mutate(L0 = NA,L1 = c(rep(RichLevels[1], RichLevels[1]), rep(NA, n()-RichLevels[1])), L2 = c(rep(RichLevels[2], RichLevels[2]), rep(NA, n()-RichLevels[2])), L3=c(rep(RichLevels[3], RichLevels[3]), rep(NA, n()-RichLevels[3]))
         , Inhibitory="Non", Rep="Rep2", DateStart=date_start) 
inhibDat_notinhib <- rbind(inhibDat_notinhib1, inhibDat_notinhib2)


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
  group_by(Inhibitory,DateStart,RichLevel,Rep) %>%
  mutate(subIsolateID=paste0("Isolate",seq(1,n()))) %>% 
  ungroup() %>%
  spread(key=subIsolateID, value=IsolateID) %>%
  unite(Inhibitory,DateStart,RichLevel,Rep, col=SampleID, remove=FALSE) 
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

metadata <- CV_set %>% full_join(Bd_set) %>% ungroup() %>%
  mutate(JarID = sample(jarNames[1:n()], replace=FALSE)) %>%
  arrange(JarID) %>%
  select(SampleID, JarID, RichLevel, Rep, Type, DateStart, everything()) %>% arrange(SampleID) %>%
  mutate(CV_raw=NA, CV=NA, Bd_n=NA)

#### Check that we haven't done this combination before.

agg_dat <- read.delim("../AggregatedData/metadata_allReps.txt")
colnames(agg_dat) <- gsub("^X","",colnames(agg_dat))

agg_dat_filt <- agg_dat %>% filter(!(Rep %in% c("Rep1","Rep2")))
agg_all_noCon <- agg_dat_filt %>%
  full_join(metadata) %>% filter(Type=="Bd", Inhibitory !="Con")
whichDup <- duplicated(agg_all_noCon[,which(colnames(agg_all_noCon) %in% inhibDat_filt$IsolateID)])
agg_all_noCon[whichDup,]

## Write out
metadata %>%
  write.table(file="jar_metadata_rep12.txt", sep="\t", quote=FALSE, row.names=FALSE)

forCultureMixing <- metadata %>% arrange(JarID) %>%mutate(JarID = paste0("Jar_",JarID)) %>%
  select(-c(SampleID, RichLevel, Rep, DateStart, Inhibitory, Type, CV_raw, CV, Bd_n)) %>%
  column_to_rownames(var="JarID") 

weights <- 1/rowSums(forCultureMixing,na.rm = TRUE)
weights[is.infinite(weights)] <- 0

forCultureMixing <- forCultureMixing * weights
forCultureMixing[is.na(forCultureMixing)] <- 0
forCultureMixing %>% rownames_to_column(var="JarID") %>%
  write.table(file="forCultureMixing.txt", sep="\t", quote=FALSE, row.names = FALSE)

## For this set, I actually did it in a different place.
ODplateOrder <- c(sample(colnames(forCultureMixing)), sample(colnames(forCultureMixing)))

remaining <- sample(c(ODplateOrder, rep("TRYP",3)))
plateMap <- data.frame(BLANK=toupper(letters[1:8]))
col_name <- 1
currentCol <- 1
while (currentCol<=12) {
  until <- ifelse(length(remaining)>=8, 8, length(remaining))
  if (until==0) {
    newcol <- rep('',8)
  }else if (until<8) {
    newcol <- c(remaining[1:until], rep('', 8-length(remaining)))
  } else {
    newcol <- remaining[1:until]
  }
  plateMap <- cbind(plateMap, newcol)
  remaining <- remaining[-c(1:until)]
  currentCol <- currentCol+1
}
colnames(plateMap) = 0:12

write.table(plateMap, quote=FALSE, row.names=F, sep="\t",file="platemap_OD_reading.txt")

isoList <- sort(colnames(forCultureMixing))
write.table(isoList, row.names=FALSE, col.names = FALSE, quote=FALSE, file="isoList.txt")

setwd("..")
