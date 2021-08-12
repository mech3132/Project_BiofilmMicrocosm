#!bin/bash

library(tidyverse)
library(lubridate)
#### DNA extraction organization ####
setwd("AggregatedData/DNA_extraction/")
## Read in sample mapping file
dat <- read.delim("../metadata_allReps.txt") %>%
  filter(!is.na(SampleID)) %>%
  mutate(JJarID = paste0("J",JarID))
colnames(dat) <- gsub("^X","",colnames(dat))

#### Preinoculation samples ####
dat_preinoc <- dat %>% select(-c(plate_bdinhib, CV_var, CV, Bd_n, cClean, JarID)) %>%
  spread(key=Type, value=JJarID) %>%
  rename(original_sampleID = SampleID, Bd_jar = Bd, CV_jar = CV) %>%
  unite(col=info_sampleID, original_sampleID, Bd_jar, CV_jar, sep="_", remove = FALSE) %>%
  mutate(sampling_method = "pellet") %>% 
  mutate(DateSampled = as.Date(DateStart) + days(4), sampleType = "Inoc") 
# Add in controls
dat_preinoc_con <- dat_preinoc %>% select(DateStart, DateSampled) %>% distinct() %>%
  mutate(info_sampleID = paste0("springH2O-for-dil_Inoc_",DateStart)
        , sampling_method = "pellet"
        , sampleType = "sample_control")

#### Post biofilm samples ####
dat_biofilm <- dat %>% filter(Type == "CV") %>% 
  mutate(sampling_method = "whole membrane"
         , sampleType = "PreBd") %>%
  rename(original_sampleID = SampleID) %>%
  unite(col=info_sampleID, original_sampleID, JJarID, sampleType, sep="_", remove=FALSE) %>% 
  select(-c(plate_bdinhib, Type, volPellet, Bd_n, JJarID))  %>%
  mutate(DateSampled = as.Date(DateStart) + days(8))
# Add in controls
dat_biofilm_con <- dat_biofilm %>% select(DateStart, DateSampled) %>% distinct() %>%
  mutate(info_sampleID = paste0("molecH2O-for-rinse_PreBd_",DateStart)
         , sampling_method = "water"
         , sampleType = "sample_control")

#### Post Bd exposure ####
dat_postbd <- dat %>% filter(Type == "Bd") %>%
  mutate(sampling_method = "swab"
         , sampleType = "Bd") %>%
  rename(original_sampleID = SampleID) %>%
  unite(col=info_sampleID, original_sampleID, JJarID, sampleType, sep="_", remove=FALSE) %>%
  select(-c(plate_bdinhib, Type, volPellet, CV, CV_var, JJarID)) %>%
  mutate(DateSampled = as.Date(DateStart) + days(12))
dat_postbd_mem <- dat %>% filter(Type == "Bd") %>%
  mutate(sampling_method = "swab"
         , sampleType = "Bdmem") %>%
  rename(original_sampleID = SampleID) %>%
  unite(col=info_sampleID, original_sampleID, JJarID, sampleType, sep="_", remove=FALSE) %>%
  select(-c(plate_bdinhib, Type, volPellet, CV, CV_var, JJarID)) %>%
  mutate(DateSampled = as.Date(DateStart) + days(12))
# Add in controls
dat_postbd_con1 <- dat_postbd %>% select(DateStart, DateSampled) %>% distinct() %>%
  mutate(info_sampleID = paste0("molecH2O-for-rinse_PostBd_",DateStart)
         , sampling_method = "water"
         , sampleType = "sample_control")
dat_postbd_con2 <- dat_postbd %>% select(DateStart, DateSampled) %>% distinct() %>%
  mutate(info_sampleID = paste0("surfaceSwab_PostBd_",DateStart)
         , sampling_method = "swab"
         , sampleType = "sample_control")


#### Merge all DNA samples together
allSamples_to_extract <- full_join(dat_preinoc,dat_preinoc_con) %>%
  full_join(dat_biofilm) %>%
  full_join(dat_biofilm_con) %>%
  full_join(dat_postbd) %>%
  full_join(dat_postbd_mem) %>%
  full_join(dat_postbd_con1) %>%
  full_join(dat_postbd_con2)

# Check row number is preserved
nrow(allSamples_to_extract)
nrow(dat_preinoc) + nrow(dat_preinoc_con) + nrow(dat_biofilm) + nrow(dat_biofilm_con) + 
  nrow(dat_postbd) + nrow(dat_postbd_mem) + nrow(dat_postbd_con1) + nrow(dat_postbd_con2)

#### Separate out into extraction "groups" ####
set.seed(8245345)

allSamples_to_extract_remaining <- allSamples_to_extract
# Set one for Andrew; do 7 membranes from Rep 3
A1_randomized <- allSamples_to_extract_remaining %>% filter(sampleType == "PreBd", Rep == "Rep3") %>% 
  mutate(randomOrder = sample(1:nrow(.))) %>%
  arrange(randomOrder) %>%
  mutate(extraction_set = c(rep("A1",4), rep("A2",3)))

allSamples_to_extract_remaining <- allSamples_to_extract %>% filter(!c(info_sampleID %in% A1_randomized$info_sampleID))

# Set two for Andrew; do 7 membranes from Rep 4
A2_randomized <- allSamples_to_extract_remaining %>% filter(sampleType == "PreBd", Rep == "Rep4") %>%
  mutate(randomOrder = sample(1:nrow(.))) %>%
  arrange(randomOrder) %>%
  mutate(extraction_set = c(rep("A2",4), rep("A1",3)))


allSamples_to_extract_remaining <- allSamples_to_extract_remaining %>% filter(!c(info_sampleID %in% A2_randomized$info_sampleID))

# Set three for Andrew; swabs and water from Pre-bd R3/4
temp_con <- allSamples_to_extract_remaining %>% filter(sampleType == "sample_control", sampling_method == "water",DateStart == "2020-09-07", DateSampled == "2020-09-15")
A3_randomized <- allSamples_to_extract_remaining %>% filter((sampleType == "Bd" & Rep %in% c("Rep4", "Rep3")) ) %>%
  full_join(temp_con) %>%
  mutate(randomOrder = sample(1:nrow(.))
         , extraction_set = "A3") %>%
  arrange(randomOrder) 
allSamples_to_extract_remaining <- allSamples_to_extract_remaining %>% filter(!c(info_sampleID %in% A3_randomized$info_sampleID))

# For Melissa
# Practise extraction
allSamples_to_extract_remaining %>% filter(sampleType == "Bdmem", Rep == "Rep3")


preBd <- allSamples_to_extract_remaining %>% filter(sampleType == "PreBd") %>% mutate(randomOrder = sample(1:nrow(.))) %>% arrange(randomOrder)
postBd <- allSamples_to_extract_remaining %>% filter(sampleType == "Bd") %>% mutate(randomOrder = sample(1:nrow(.))) %>% arrange(randomOrder)
watercon <- allSamples_to_extract_remaining %>% filter(sampleType == "sample_control", sampling_method == "water") %>% mutate(randomOrder = sample(1:nrow(.))) %>% arrange(randomOrder)
swabcon <- allSamples_to_extract_remaining %>% filter(sampleType == "sample_control", sampling_method == "swab") %>% mutate(randomOrder = sample(1:nrow(.))) %>% arrange(randomOrder)

M1_randomized <- full_join(preBd[1:10,], postBd[1:10,]) %>%
  full_join(watercon[1,]) %>% full_join( swabcon[1:2,]) %>% mutate(extraction_set = "M1")
M2_randomized <- full_join(preBd[11:20,], postBd[11:20,]) %>%
  full_join(watercon[2,]) %>% full_join( swabcon[3:4,]) %>% mutate(extraction_set = "M2")
M3_randomized <- full_join(preBd[21:30,], postBd[21:30,]) %>%
  full_join(watercon[3,]) %>% full_join( swabcon[5:6,]) %>% mutate(extraction_set = "M3")
M4_randomized <- full_join(preBd[31:40,], postBd[31:40,]) %>%
  full_join(watercon[4:6,]) %>% mutate(extraction_set = "M4")
M5_randomized <- full_join(preBd[41:50,], postBd[41:50,]) %>%
  full_join(watercon[7:9,]) %>% mutate(extraction_set = "M5")
M6_randomized <- full_join(preBd[51:61,], postBd[51:61,]) %>%
  full_join(watercon[10,]) %>% mutate(extraction_set = "M6")
M7_randomized <- full_join(preBd[62:70,], postBd[62:71,]) %>%
  full_join(watercon[11,]) %>% mutate(extraction_set = "M7")

# num_per_set <- 19
# allSamples_to_extract_randomized <- allSamples_to_extract %>% 
#   filter(!(sampleType %in% c("Bdmem", "Inoc")), !(sampleType == "sample_control"  & sampling_method == "swab")) %>%
#   mutate(random_extraction_order = sample(1:nrow(.), replace=FALSE)) %>%
#   rowwise() %>%
#   mutate(numset = which(random_extraction_order <= seq(0, to = ceiling(nrow(.)/num_per_set)*num_per_set, by=num_per_set)[-1])[1]
#          , extr_set = LETTERS[numset]) %>% 
#   arrange(random_extraction_order) 
# 
# ## add in extraction controls
# allSamples_to_extract_randomized_withcon <- allSamples_to_extract_randomized %>% select(extr_set) %>% distinct() %>%
#   mutate(info_sampleID = paste0("extr_con_",extr_set), sampleType = "extr_con") %>%
#   full_join(allSamples_to_extract_randomized) %>%
#   arrange(Rep, JarID) %>%
#   select(info_sampleID, original_sampleID, everything())


## Number of rows if we subtract the samples above
# 365-nrow(dat_postbd_con2)- nrow(dat_postbd_mem)-nrow(dat_preinoc)


allSamples_to_extract_randomized_withcon <- A1_randomized %>%
  full_join(A2_randomized) %>% full_join(A3_randomized) %>%
  full_join(M1_randomized) %>% full_join(M2_randomized) %>% full_join(M3_randomized) %>%
  full_join(M4_randomized) %>% full_join(M5_randomized) %>% full_join(M6_randomized) %>%
  full_join(M7_randomized)

# nrow(allSamples_to_extract_randomized_withcon)
#### For printing ####
extraction_plan <- allSamples_to_extract_randomized_withcon %>%
  mutate(randomOrder = sample(1:nrow(.))) %>%
  arrange(extraction_set) %>% group_by(extraction_set) %>%
  mutate(rank_within = rank(randomOrder)) %>% ungroup() %>%
  select(info_sampleID, extraction_set, rank_within) %>% 
  pivot_wider(names_from = extraction_set, values_from = info_sampleID) %>%
  arrange(rank_within)

write.table(extraction_plan, file = "./extraction_plan.txt",quote = FALSE, sep = "\t", row.names = FALSE)
write.table(allSamples_to_extract_randomized_withcon, file = "./extraction_metadata.txt",quote = FALSE, sep = "\t", row.names = FALSE)
setwd("../../")
