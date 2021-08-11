#!bin/bash

##### NOTE: THIS WAS USED FOR FIRST REP WHEN ALL ISOLATE RESULTS HADN'T COME IN YET ####
### Combining info for all isolates ###
library(tidyverse)

inhib <- read.delim("00_isolate_info/Alex_inhib_assay/inhibitory_PA_edited.txt", header=TRUE)
biofilm <- read.delim("00_pilots/Biofilm_assay/biofilm_assay_allisolates.txt")
biofilm_former <- read.delim("00_pilots/Biofilm_assay/list_se_sig_biofilmformers.txt")

inhib_collapsed <- inhib %>%
  select("Sample_ID", "Rep1_Measureday3", "Rep2_Measureday3") %>%
  rename(Rep1 = Rep1_Measureday3, Rep2 = Rep2_Measureday3, IsolateID = Sample_ID) %>%
  gather(-c(IsolateID), key=Replicate, value=InhibZonemm) %>%
  group_by(IsolateID) %>% summarize(aveInhibZonemm = mean(InhibZonemm, na.rm = TRUE))

CV_collapsed <- biofilm %>% 
  select(IsolateID, g_sp, CV) %>%
  group_by(IsolateID, g_sp) %>% summarize(aveCV = mean(CV, na.rm=TRUE))

bfFormer <- biofilm_former %>%
  select(IsolateID) %>% mutate(BiofilmFormer=TRUE)

growth_collapsed <- biofilm %>% 
  select(IsolateID, growth_midpoint1, growth_midpoint2, growth_midpoint3) %>% distinct() %>%
  gather(-IsolateID, key=Replicate, value=GrowthMidpoint) %>%
  group_by(IsolateID) %>% summarize(medGrowthMidpoint = median(GrowthMidpoint, na.rm=TRUE))

## Need to add biofilm former Y/N metric


allIsolateInfo <- CV_collapsed %>% full_join(inhib_collapsed) %>% full_join(growth_collapsed) %>% full_join(bfFormer) %>%
  filter(IsolateID !="#N/A") %>%
  mutate(BiofilmFormer = ifelse(is.na(BiofilmFormer), FALSE, BiofilmFormer)
        , Inhibitory = ifelse(aveInhibZonemm>=1, TRUE, ifelse(aveInhibZonemm==0, FALSE, NA))) 

# write.table(allIsolateInfo, sep="\t", row.names = FALSE, quote=FALSE, file="all_isolate_info_combined.txt")

write.table(allIsolateInfo, sep="\t", row.names = FALSE, quote=FALSE, file="all_isolate_info_limited_for_rep1.txt")



