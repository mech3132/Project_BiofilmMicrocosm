#!bin/bash
# library(MASS) # fitting normal distributions
library(tidyverse)
#### This aggregates all data ####
setwd("03_datacleaning")

allReps <- "../01_experiment_metadata/downstream/aggregated_metadata.txt"

CV <- "../02_externally_generated_data/CV/downstream/CV_adj.txt"
Bd_count <- "../02_externally_generated_data/Bd_counting/downstream/dat_count_agg.txt" 
Bd_measure <- "../02_externally_generated_data//Bd_counting/downstream/dat_meas_agg.txt"
qPCR <- "../02_externally_generated_data/qPCRs/downstream/qPCR_results_combined_edited.txt"
tax <- "01_seq_data_cleaning/downstream/asv_metadata.txt"
otu <- "01_seq_data_cleaning/downstream/metadata_with_sequence.txt"


allinhib <- "../00_isolate_info/manual_chosen_isolates/Inhib_to_include.txt"
allnoninhib <- "../00_isolate_info/manual_chosen_isolates/nonInhib_to_include.txt"

blastmatchPWD <- "../02_externally_generated_data/16S_sequencing/downstream/allInhibData.txt"
dir.create("02_merge_all_data/")
#### Load files #####

allReps_dat <- read.delim(allReps)
CV_dat <- read.delim(file=CV)
qPCR_dat <- read.delim(file=qPCR)
Bd_count_dat <- read.delim(file=Bd_count)
Bd_count_dat <- Bd_count_dat %>%rename(sampleID_qPCR = sampleID) # Rename to match others
Bd_measure_dat <- read.delim(file=Bd_measure)
Bd_measure_dat <- Bd_measure_dat %>% rename(sampleID_qPCR = sampleID) # rename to match others
seq_dat <- read.delim(otu)
taxa <- read.delim(tax)
allIso <- c(read.delim(allinhib, header=FALSE)$V1, read.delim(allnoninhib, header=FALSE)$V1)
blastmatch <- read.delim(blastmatchPWD)
#### Edit ####
colnames(allReps_dat) <- gsub("^X","",colnames(allReps_dat))
allReps_dat_edited <- allReps_dat %>% 
  mutate(volPellet = ifelse(is.na(volPellet), 1, volPellet)) %>%
  select(SampleID, sampleID_qPCR, JarID, DateStart, Rep, Inhibitory, RichLevel, introRich, Type, one_of(allIso), volPellet, cClean)

#Rename isolate IDs so you can just grep them
tempColnames <- colnames(allReps_dat_edited[which(colnames(allReps_dat_edited) %in% allIso)])
colnames(allReps_dat_edited)[which(colnames(allReps_dat_edited) %in% allIso)] <- paste0("ISO_",tempColnames)
#### Alter datasets to keep naming convention consistent and succinct ####

## qPCR ##
qPCR_dat_edited <- qPCR_dat %>% select(sampleID, ave_log10copies, sd, lwr_95, upr_95) %>%
  rename(sampleID_qPCR = sampleID
         , qPCR_log10 = ave_log10copies
         , qPCR_sd = sd
         , qPCR_lwr_95_total = lwr_95
         , qPCR_upr_95_total = upr_95)

## CV ##
CV_dat_edited <- CV_dat %>%
  rename(CV=CV_adj, CV_sd = sd_CV_adj)

## Bd microscopy ##
Bd_count_dat_edited <- Bd_count_dat %>%
  rename(Micro_Bd_density=zoosp_dens_micro, Micro_Bd_density_sd=zoosp_count_sd_micro, Micro_Bd_truezero=Bd_true_zero) %>%
  select(sampleID_qPCR, Rep, JarID, Micro_Bd_density,Micro_Bd_density_sd,Micro_Bd_truezero, clouds_present, pinpricks_present)
  
Bd_measure_dat_edited <- Bd_measure_dat %>% 
  rename(Micro_Bd_meansize = ave_zoosp_size_micro, Micro_Bd_meansize_sd = sd_zoosp_size_micro,  Micro_Bd_medsize = med_zoosp_size_micro) %>%
  select(sampleID_qPCR, Rep, JarID, Micro_Bd_meansize, Micro_Bd_meansize_sd, Micro_Bd_medsize)

### Merge all data

allReps_merged <- allReps_dat_edited %>%
  left_join(CV_dat_edited) %>%
  left_join(qPCR_dat_edited) %>%
  left_join(Bd_count_dat_edited) %>%
  left_join(Bd_measure_dat_edited) %>% 
  mutate(Micro_Bd_truezero = ifelse(!is.na(Micro_Bd_medsize)|clouds_present|pinpricks_present, FALSE, Micro_Bd_truezero)) %>%
  left_join(seq_dat) 

# %>%
#   select(-clouds_present, -pinpricks_present)
# 
# allReps_withqPCR <- allReps %>%
#   # mutate(introRich = introRich) %>%
#   left_join(qPCR_dat_edited) %>%
#   left_join(Bd_count_dat) %>%
#   left_join(Bd_measure_dat) %>% 
#   # left_join(seq_dat) %>%
#   mutate(Bd_true_zero = ifelse(!is.na(med_zoosp_size_micro)|clouds_present|pinpricks_present, FALSE, Bd_true_zero))

## Quick peak at qPCR data
allReps_merged %>% 
  ggplot() + geom_point(aes(x=RichLevel, y=qPCR_log10, col=Type))
# Looks good, but noisy.

## Quick peak at richness data
allReps_merged %>% filter(Type == "CV") %>%
  ggplot() + geom_jitter(aes(x=RichLevel, y=obsRich), width=0.1, height=0)

##### Clean data and collapse into "replicate" ######
# First, is there anything wrong with cClean = FALSE?
sum(allReps_merged$cClean==FALSE)
allReps_merged %>% filter(is.na(qPCR_log10))
allReps_merged %>% 
  ggplot() + geom_jitter(aes(x=RichLevel, y=qPCR_log10, col=cClean), width=0.2, height=0, alpha=0.3) +
  facet_grid(Type~.)
allReps_merged %>% filter(Type=="CV") %>%
  ggplot() + geom_jitter(aes(x=RichLevel, y=CV, col=cClean), width=0.2, height=0, alpha=0.3) +
  facet_grid(Type~.)

# Does volPellet have an effect on anything?
allReps_merged %>% 
  ggplot() + geom_jitter(aes(x=RichLevel, y=qPCR_log10, col=volPellet), width=0.2, height=0, alpha=0.3) +
  facet_grid(Type~.)
allReps_merged %>% filter(Type=="CV") %>%
  ggplot() + geom_jitter(aes(x=RichLevel, y=CV, col=volPellet), width=0.2, height=0, alpha=0.3) +
  facet_grid(Type~.)
# Nope
allReps_merged %>% 
  ggplot() + geom_jitter(aes(x=RichLevel, y=qPCR_log10, col=factor(REMOVE_SEQCONTAM)), width=0.2, height=0, alpha=0.3) +
  facet_grid(Type~.)
allReps_merged %>% filter(Type=="CV") %>%
  ggplot() + geom_jitter(aes(x=RichLevel, y=CV, col=factor(REMOVE_SEQCONTAM)), width=0.2, height=0, alpha=0.3) +
  facet_grid(Type~.)
allReps_merged %>% filter(Type=="Bd") %>%
  ggplot() + geom_jitter(aes(x=RichLevel, y=qPCR_log10, col=factor(REMOVE_SEQCONTAM)), width=0.2, height=0, alpha=0.3) +
  facet_grid(Type~.)

####  Rules for removing outliers: ####
# The individual ASV plots (byRep_ASV_thresholds) will dictate which samples are discarded.
# A sample will be discarded if:
# (1) The threshold algorithm appears to do a poor job of identifying the "real" ASVs
# (2) The separation between "real" and "noise" ASVs looks particularily diffult to parse out
# (3) The sample looks entirely contaminated with bacteria.

### This a manually curated list of jars that we should remove because it looks contaminated
toRemoveSamples <-c("R1_J27_Pre", "R1_J22_Pre" # Rich 10 Inhib. ASVs look very different between pre and post
                    ,"R3_J7_Pre", "R3_J1_Bd" # Rich 1, Inhib, contaminated
                    ,"R3_J24_Pre" ,"R3_J18_Bd" # Rich 3, Inhib, contaminated
                    ,"R3_J4_Pre", "R3_J12_Bd" # Rich 10, Inhib, contaminated
                    , "R4_J17_Bd" # Rich 10, Inhib, contaminated. Can't tell low from high
                    # , "R6_J7_Bd" # Keep; looks okay
                    , "R6_J24_Bd"
                    , "R7_J15_Bd" # Paenibacillus contam
                    # , "R7_J21_Pre" # can't tell what's real and not
                    , "R8_J7_Pre" # Rich 3, Non, just looks weird. Could also get rid of entier R8 non
                    , "R10_J1_Bd"# Rich1, Non, looks nothing like Pre and pre matches rest
                    , "R12_J18_Pre", "R12_J14_Bd" # Rich1, Non, Looks nothing like rest of samples
)
allReps_merged_filt <- allReps_merged %>% filter(!(sampleID_qPCR %in% toRemoveSamples))

##### Matching between pre and post ######
## Filter out ASVs that are not present in the dataset; compare ASVs that don't match pre and post
tempASVs <- names(which(colSums(allReps_merged_filt %>% select(starts_with("ASV_")))==0))
allReps_merged_filt2 <- allReps_merged_filt %>% select(-one_of(tempASVs))
actualASVsKept <- allReps_merged_filt2 %>% select(starts_with("ASV")) %>% colnames()

# Remove ASVs that are no longer present in filtered dataset
tempASV_CV <- names(which(colSums(allReps_merged_filt2 %>% filter(Type=="CV") %>% select(starts_with("ASV_")))>0))
tempASV_Bd <- names(which(colSums(allReps_merged_filt2 %>% filter(Type=="Bd") %>% select(starts_with("ASV_")))>0))

# Compare ones that differ
inCV_not_Bd <- tempASV_CV[!(tempASV_CV %in% tempASV_Bd)]
inBd_not_CV <- tempASV_Bd[!(tempASV_Bd %in% tempASV_CV)]
# Let's look at the samples where this happens
temp_nonSharedASVs <- allReps_merged_filt2 %>% #select(SampleID, sampleID_qPCR, Rep, Inhibitory, RichLevel, Type) %>%
  mutate(total = rowSums(across(.cols=all_of(c(inCV_not_Bd, inBd_not_CV))))) %>% filter(total>0)
relaventASVs <- names(which(temp_nonSharedASVs %>% select(starts_with("ASV")) %>% colSums() >0))
temp_nonSharedASVs %>% select(SampleID, sampleID_qPCR, Rep, Inhibitory, RichLevel, Type, one_of(relaventASVs)) %>%
  pivot_longer(starts_with("ASV"), names_to= "ASVID", values_to="Reads") %>%
  left_join(taxa %>% select(ASVID, ScientificName_unique)) %>%
  mutate(Reads = ifelse(Reads==0, NA, Reads)) %>%
  rowwise() %>% mutate(notShared = ASVID %in% c(inCV_not_Bd, inBd_not_CV)) %>% ungroup() %>%
  ggplot() + geom_bar(aes(y=Reads, x=ScientificName_unique, fill=notShared), stat = "identity") +
  facet_wrap(Type~sampleID_qPCR , scales = "free")+ theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
### MANUAL NOTES:
# All three of these can probably be kept, since they are either inconsequential, or could be "real" and just
# disspeared in one of the other samples due to filering. 


### separate out; and then merge
allReps_CV <- allReps_merged_filt2 %>% filter(Type=="CV") %>% 
  rename(qPCR_log10_bact = qPCR_log10, qPCR_bact_sd = qPCR_sd, qPCR_bact_upr_95_total = qPCR_upr_95_total, qPCR_bact_lwr_95_total = qPCR_lwr_95_total
         , B_JarID = JarID, B_obsRich=obsRich) %>%
  select(-c(sampleID_qPCR, cClean, volPellet, Type, REMOVE_SEQCONTAM, starts_with("Micro"), starts_with("ASV_"), starts_with("ISO_"), clouds_present, pinpricks_present))

allReps_Bd <- allReps_merged_filt2 %>% filter(Type=="Bd") %>% 
  rename(qPCR_log10_Bd = qPCR_log10, qPCR_Bd_sd = qPCR_sd, qPCR_Bd_upr_95_total = qPCR_upr_95_total, qPCR_Bd_lwr_95_total = qPCR_lwr_95_total
         , A_JarID = JarID, A_obsRich=obsRich) %>%
  select(-c(sampleID_qPCR, cClean, volPellet, Type, REMOVE_SEQCONTAM, starts_with("CV"), starts_with("ASV_"), starts_with("ISO_")))


allReps_collapsed <- full_join(allReps_CV, allReps_Bd)

write.table(allReps_merged_filt2, quote=FALSE, sep="\t", row.names = FALSE, file="02_merge_all_data/metadata_allReps.txt")
write.table(allReps_collapsed, quote=FALSE, sep="\t", row.names = FALSE, file="02_merge_all_data/metadata_allReps_collapsed.txt")

#### Alex blast results ####
remainingASVs <- allReps_merged_filt %>% select(starts_with("ASV")) %>% colnames()
remainingASVs <- gsub("ASV_","",remainingASVs)
blast_matches <- blastmatch %>% filter(ASVID %in% remainingASVs) %>% select(ASVID, avePercIdent_blast_alex, avePercIdent_blast_woodhams, Matched_AlexID)


write.table(blast_matches, quote=FALSE, sep = "\t", row.names = FALSE, file="02_merge_all_data/blast_matches.txt")
# number of alex matches
sum(!is.na(blast_matches$Matched_AlexID))
# Number of unique alex matches
length(unique(blast_matches$Matched_AlexID[!is.na(blast_matches$Matched_AlexID)]))

setwd("..")

