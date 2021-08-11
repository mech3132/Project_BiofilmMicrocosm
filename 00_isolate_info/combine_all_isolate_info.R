#!bin/bash

### Combining info for all isolates ###
library(tidyverse)

inhib <- read.delim("00_isolate_info/Alex_inhib_assay/Master_bd_inhibition_data_filtered31aug2020.txt", header=TRUE)
biofilm <- read.delim("00_pilots/Biofilm_assay/biofilm_assay_allisolates.txt")
biofilm_former <- read.delim("00_pilots/Biofilm_assay/list_se_sig_biofilmformers.txt")

inhib_collapsed <- inhib %>%
  select("Isolate_ID", "Max_inhib_zone") %>%
  rename(aveInhibZonemm=Max_inhib_zone
         , IsolateID = Isolate_ID)

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
        , Inhibitory = ifelse(aveInhibZonemm>=1, TRUE, ifelse(aveInhibZonemm==0, FALSE, NA))) %>% 
  mutate(medGrowthMidpoint = as.numeric(medGrowthMidpoint))
write.table(allIsolateInfo, sep="\t", row.names = FALSE, quote=FALSE, file="00_isolate_info/all_isolate_info_combined.txt")

