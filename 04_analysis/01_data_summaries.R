#!bin/bash
##########  Adjusting files and filtering out things ############

library(tidyverse)
setwd("04_analysis/")
dir.create("01_data_summaries")
##### Loading ####
# dat_raw <- read.delim("../03_datacleaning/downstream/metadata_allReps_collapsed.txt")
# dat_split_raw <- read.delim("../03_datacleaning/downstream/metadata_allReps.txt")
dat_raw <- read.delim("../03_datacleaning/05_final_seq_filtering/downstream/asv_and_meta_merged_collapsed.txt")
dat_split_raw <- read.delim("../03_datacleaning/05_final_seq_filtering/downstream/asv_and_meta_merged.txt")
dat_inhibmeta <- read.delim("../03_datacleaning/05_final_seq_filtering/downstream/metadata_inhibitory_summary.txt")
isolate_info <- read.delim("../03_datacleaning/05_final_seq_filtering/downstream/isolate_info_filtered.txt")

#### Renaming; calculating all required columns ####
dat <- dat_raw %>%
  # mutate(CV_lwrlog10 = log10(CV-CV_sd*2 + 0.01), CV_uprlog10 = log10(CV+CV_sd*2 + 0.01)) %>%
  # mutate(qPCR_bact_lwrlog10 = log10(10^qPCR_log10_bact -qPCR_bact_sd*2 + 1), qPCR_bact_uprlog10 = log10(10^qPCR_log10_bact +qPCR_bact_sd*2 + 1)) %>%
  mutate(CV_log10 = log10(CV+0.01), Bd_micro_log10 = log10(Micro_Bd_density+1) , qPCR_bact_log10_inhib = log10(inhibReads+1)) %>%
  rename(Bd_qPCR_log10 = qPCR_log10_Bd, qPCR_bact_log10 = qPCR_log10_bact, Richness = B_obsRich) %>%
  rename(inhibRich = B_inhibRich, inhibRichFrac = B_inhibRichFrac, shannon=B_shannon, pielou = B_pielou) %>%
  # mutate(Bd_prod_log10 = log10((10^qPCR_log10_Bd)/(10^Micro_Bd_density_log10))) %>%
  mutate(DateStart = factor(DateStart))  %>%
  mutate(Inhibitory = ifelse(Inhibitory=="Inhib",1,0)) %>%
  mutate(Inhib = ifelse(inhibRich == 0, 0, 1)) %>%
  select(SampleID, Rep, DateStart, RichLevel,Inhibitory,Inhib, CV_log10, qPCR_bact_log10, qPCR_bact_log10_inhib, Richness, Bd_micro_log10, Bd_qPCR_log10, inhibRich, inhibRichFrac, inhibReads, inhibProp, shannon)

#### Creating output file ####
LOG = "./01_data_summaries/log.txt"
cat("SUMMARY", file = LOG,sep="\n\n")

#### Removing NAs and deciding on final dataset ####
cat(paste0("Number of hypothetical replicates: ", 12*7), file=LOG ,sep="\n\n", append = TRUE)
nrow(dat)
cat(paste0("Number of partial replicates: ",nrow(dat)), file=LOG ,sep="\n\n", append = TRUE)
# Filter out incomplete cases
dat_filt<- dat %>% drop_na()
nrow(dat_filt)
cat(paste0("Number of complete replicates: ",nrow(dat_filt)), file=LOG ,sep="\n\n", append = TRUE)
# Number of each type of treatment
sink(LOG, append = TRUE)
cat("\n\n")
dat_filt %>% select(RichLevel, Inhibitory) %>% table()
sink()
# Ranges of each variable
sink(LOG, append=TRUE)
cat("\n\n")
dat_filt %>% 
  mutate(CV = 10^CV_log10-0.01, qPCR_bact = 10^qPCR_bact_log10-1, Bd_micro = 10^Bd_micro_log10-1, Bd_qPCR = 10^Bd_qPCR_log10-1)%>%
  select(CV, qPCR_bact, Richness, Bd_micro, Bd_qPCR, inhibRich, inhibProp, shannon) %>%
  pivot_longer(cols=c(CV, qPCR_bact, Richness, Bd_micro, Bd_qPCR, inhibRich, inhibProp, shannon), names_to = "Variable", values_to = "Value") %>%
  group_by(Variable) %>%
  summarize(Mean = mean(Value), SD = sd(Value), Low=min(Value), High=max(Value)) %>%
  as.data.frame()
sink()

## Getting ASV counts
sink(LOG, append = TRUE)
cat("\n\n")
cat("Number of ASVs in final dataset:")
dat_split_raw %>% select(starts_with("ISO")) %>% ncol()
sink()


### number of inhibitory vs not
sink(LOG, append = TRUE)
cat("\n\n")
cat("Range of inhibitory values:")
apply(dat_inhibmeta[,-1], 2, function(x) c(mean=mean(x), min=min(x), max=max(x)))
sink()

### ASV counts by inhib/noninhib
sink(LOG, append = TRUE)
cat("\n\n")
cat("Number of inhibitory vs non inhibitory :")
isolate_info %>% select(Inhibitory) %>% table()
sink()

#### Make full version #####
dat_full <- dat_split_raw %>% filter(Type=="CV") %>% select(SampleID, starts_with("ISO_")) %>%
  right_join(dat_filt)

# dat_full %>%  select(SampleID, starts_with("ISO"), Richness) %>% View()
#### Centering variables ####
dir.create("01_data_summaries/downstream")
# dat_numeric <- dat_filt %>% select(Bd_micro_log10, Bd_qPCR_log10, qPCR_bact_log10, CV_log10, Richness)
dat_numeric <- dat_filt %>% select(qPCR_bact_log10, CV_log10, Richness, inhibRich, inhibRichFrac, inhibProp, inhibReads, shannon, qPCR_bact_log10_inhib)
dat_numeric_scaled <- scale(dat_numeric, scale=FALSE)
dat_scaled <- as.data.frame(dat_numeric_scaled)
dat_centers <- attributes(dat_numeric_scaled)$`scaled:center`
write.table(dat_centers, "01_data_summaries/downstream/dat_centers.txt", row.names = FALSE, quote = FALSE, sep = "\t" )
# Add back other variables, add centers back
names(dat_centers) <- paste0("Centre_", names(dat_centers))
t(as.data.frame(dat_centers))
dat_scaled$SampleID <- dat_filt$SampleID
dat_scaled$Rep <- dat_filt$Rep
dat_scaled$DateStart <- dat_filt$DateStart
dat_scaled$RichLevel <- dat_filt$RichLevel
dat_scaled$Inhibitory <- dat_filt$Inhibitory
dat_scaled$Inhib <- dat_filt$Inhib
dat_scaled$inhibProp_raw <- dat_filt$inhibProp
dat_scaled$inhibRichFrac_raw <- dat_filt$inhibRichFrac
dat_scaled$Bd_micro_log10 <- dat_filt$Bd_micro_log10
dat_scaled$Bd_qPCR_log10 <- dat_filt$Bd_qPCR_log10
dat_scaled_final <- dat_scaled %>% select(SampleID, Rep, DateStart, RichLevel, Inhibitory, everything()) %>%
  cbind(t(as.data.frame(dat_centers)))

write.table(dat_scaled_final, file="01_data_summaries/downstream/dat_scaled_final.txt", row.names = FALSE, quote = FALSE, sep = "\t" )
write.table(dat, file="01_data_summaries/downstream/dat_nonscaled_final.txt", row.names = FALSE, quote = FALSE, sep = "\t" )
write.table(dat_full, file="01_data_summaries/downstream/dat_full.txt", row.names = FALSE, quote = FALSE, sep = "\t" )


setwd("..")

