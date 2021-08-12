#!bin/bash

###### Filtering out contaminants ###########
library(tidyverse)
library(ape)
library(phytools)

setwd("03_datacleaning/")
dir.create("01_seq_data_cleaning")
####### Load files ################

otu <- read.delim("../02_externally_generated_data/16S_sequencing/downstream/otu_table.txt", row.names = 1)
inhibMeta <- read.delim("../02_externally_generated_data/16S_sequencing/downstream/allInhibData.txt")
tree <- read.tree("../02_externally_generated_data/16S_sequencing/downstream/filtered_tree.tre")
taxa <- read.delim("../02_externally_generated_data/16S_sequencing/downstream/taxonomy.txt")
qPCR <- read.delim("../02_externally_generated_data/qPCRs/downstream/qPCR_results_combined_edited.txt")
# alexIsolates <- read.delim("../02_externally_generated_data/16S_sequencing/intermediate_files/allIsolates.txt", header = FALSE)
alex_isos <- read.delim("../02_externally_generated_data/16S_sequencing/intermediate_files/imported_databases/alexKeepSamesunique.txt")
extraction_meta <- read.delim("../02_externally_generated_data/DNA_extraction/DNA_platemap_forJess.txt")

meta <- read.delim("../01_experiment_metadata/downstream/aggregated_metadata.txt") 
colnames(meta) <- gsub("^X","", colnames(meta))

########## OTU table adjustments #############

### Parse out alex ID's
alex_isos <- alex_isos %>% rename(ASVID = X.SampleID) %>% separate(ASVID, into = c("IsolateID"), sep = "-", remove=FALSE, extra = "drop")
## combine metadata for qpcr and experimental
qPCR_wide <- qPCR %>% 
  mutate(qPCR_copies = 10^(ave_log10copies)) %>% select(sampleID, primers, qPCR_copies, lwr_95, upr_95) %>% 
  pivot_wider(names_from = primers, values_from = qPCR_copies) %>%
  rename(qPCR_bact = `515F/806R`, qPCR_Bd = `ITS1-3/5.8SChytr`) %>%
  rename(sampleID_qPCR= sampleID)
full_meta <- meta %>% left_join(qPCR_wide)

# Get all OTUs
allOTUs <- rownames(otu)
# Rename OTU samplenames
colnames(otu) <- gsub(".","_", gsub("..","-", colnames(otu), fixed = TRUE), fixed=TRUE)

# Make relative abundance, PA
otu_relabund <- apply(otu, MARGIN = 2, FUN = function(x) x/sum(x))
otu_presenceabsence <- otu
otu_presenceabsence[otu_presenceabsence>0] <- 1
#Transpose so we can combine with qPCR
otu_with_qPCR <- t(otu_relabund) %>% as.data.frame() %>% rownames_to_column(var = "sampleID_qPCR") %>% as_tibble() %>%left_join(full_meta) %>%
  filter(!is.na(qPCR_bact)) %>% select(sampleID_qPCR, one_of(allOTUs), qPCR_bact) %>% as.data.frame()

rownames(otu_with_qPCR) <- otu_with_qPCR$sampleID_qPCR

otu_real_abund <- apply(otu_with_qPCR[,allOTUs], MARGIN = 2, FUN = function(x) x*otu_with_qPCR$qPCR_bact)
# Check I did that math right; they should be equal below
# round(rowSums(otu_real_abund)) == round(otu_with_qPCR$qPCR_bact)
# Now round each individually
otu_real_abund_round <- round(otu_real_abund)



##### Need to do some sanity checks, to make sure all data is lining up properly #####
# Make an OTU table with raw values
otu_raw <- t(otu) %>% as.data.frame() %>% rownames_to_column(var = "sampleID_qPCR")
# Do a version with adjust values
otu_prop <- t(otu_relabund) %>% as.data.frame() %>% rownames_to_column(var = "sampleID_qPCR")
# OTUs with qPCR-adjusted values
otu_freq <- otu_real_abund_round %>% as.data.frame() %>% rownames_to_column(var = "sampleID_qPCR")
# OTU presence/absence table
otu_pa <- t(otu_presenceabsence) %>% as.data.frame() %>% rownames_to_column(var = "sampleID_qPCR")

### First, check that sample qPCRs make sense
full_meta %>% filter(!is.na(qPCR_bact)) %>% ggplot() + geom_point(aes(x=introRich, y = qPCR_bact)) +
  scale_y_log10() + geom_hline(aes(yintercept = 1000), lty=2, col="red")
# It looks like things that are less than 1000 are generally just noise-- and this looks like it could all make sense.

## Next, look at raw data
otu_pa %>% mutate(obsRich = rowSums(across(one_of(allOTUs)))) %>% select(sampleID_qPCR, obsRich) %>%
  left_join(full_meta) %>% select(sampleID_qPCR, obsRich, introRich)

# Let's look at the composition of introRIch==0 and see what it looks like
otu_prop %>% left_join(full_meta) %>% 
  mutate(Type = ifelse(is.na(SampleID), "Control",Type)) %>% 
  filter(Type %in% c("CV", "Control")) %>%
  select(SampleID, sampleID_qPCR, introRich, RichLevel, Rep, Inhibitory, one_of(allOTUs)) %>% 
  mutate(Rep = factor(Rep, levels = c("Rep1","Rep2","Rep3","Rep4","Rep5","Rep6","Rep7","Rep8","Rep9","Rep10","Rep11","Rep12","Rep13"))) %>%
  pivot_longer(cols = one_of(allOTUs), names_to="ASVID", values_to = "propReads") %>%
  ggplot() + geom_bar(aes(x=sampleID_qPCR, y = propReads, fill=ASVID), stat="identity", show.legend = FALSE) + 
  facet_grid(. ~ RichLevel, drop = TRUE, scales="free_x", space = "free_x") +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))

# Look at raw OTU richness in all samples
# otu_pa[,-match(allOTUs,colnames(otu_pa))] # Checking that allOTUs is, in fact, "all OTUs"
data.frame(sampleID_qPCR=otu_pa$sampleID_qPCR, observed_otus= rowSums(otu_pa[,allOTUs])) %>% left_join(full_meta) %>% 
  mutate(Type = ifelse(is.na(SampleID), "Control",Type)) %>% 
  filter(Type %in% c("CV", "Control")) %>%
  select(SampleID, sampleID_qPCR, introRich, RichLevel, observed_otus, Rep, Inhibitory) %>% 
  mutate(Rep = factor(Rep, levels = c("Rep1","Rep2","Rep3","Rep4","Rep5","Rep6","Rep7","Rep8","Rep9","Rep10","Rep11","Rep12","Rep13"))) %>%
  # pivot_longer(cols = one_of(allOTUs), names_to="ASVID", values_to = "propReads") %>%
  ggplot() + geom_bar(aes(x=sampleID_qPCR, y = observed_otus), stat="identity") + 
  facet_grid(. ~ RichLevel, drop = TRUE, scales="free_x", space = "free_x") +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))

# Look at raw read amounts-- I know this is flawed because of proportion issues, but I just want to see if low biomass samples ACTUALLY have less reads
otu_raw %>% left_join(full_meta) %>% 
  mutate(Type = ifelse(is.na(SampleID), "Control",Type)) %>% 
  filter(Type %in% c("CV", "Control")) %>%
  select(SampleID, sampleID_qPCR, introRich, RichLevel, Rep, Inhibitory, one_of(allOTUs)) %>% 
  mutate(Rep = factor(Rep, levels = c("Rep1","Rep2","Rep3","Rep4","Rep5","Rep6","Rep7","Rep8","Rep9","Rep10","Rep11","Rep12","Rep13"))) %>%
  pivot_longer(cols = one_of(allOTUs), names_to="ASVID", values_to = "rawReads") %>%
  ggplot() + geom_bar(aes(x=sampleID_qPCR, y = rawReads, fill=ASVID), stat="identity", show.legend = FALSE) + 
  facet_grid(. ~ RichLevel, drop = TRUE, scales="free_x", space = "free_x") +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))
otu_raw %>% left_join(full_meta) %>% 
  mutate(Type = ifelse(is.na(SampleID), "Control",Type)) %>% 
  filter(Type %in% c("Bd", "Control")) %>%
  select(SampleID, sampleID_qPCR, introRich, RichLevel, Rep, Inhibitory, one_of(allOTUs)) %>% 
  mutate(Rep = factor(Rep, levels = c("Rep1","Rep2","Rep3","Rep4","Rep5","Rep6","Rep7","Rep8","Rep9","Rep10","Rep11","Rep12","Rep13"))) %>%
  pivot_longer(cols = one_of(allOTUs), names_to="ASVID", values_to = "rawReads") %>%
  ggplot() + geom_bar(aes(x=sampleID_qPCR, y = rawReads, fill=ASVID), stat="identity", show.legend = FALSE) + 
  facet_grid(. ~ RichLevel, drop = TRUE, scales="free_x", space = "free_x") +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))
# And finally, let's look at "true" read amount
otu_freq %>% left_join(full_meta) %>%
  pivot_longer(cols = one_of(allOTUs), names_to = "ASVID", values_to = "Reads") %>%
  ggplot(aes(x = sampleID_qPCR, y = Reads, fill = ASVID)) + geom_bar(stat = "identity", show.legend = FALSE) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  facet_grid (.~ RichLevel, drop = TRUE, scales = "free_x", space =)

# ggsave("01_seq_data_cleaning/rawReads_agains_sampleID.png", height=5, width=13
#        ,otu_raw %>% left_join(full_meta) %>% 
#          mutate(Type = ifelse(is.na(SampleID), "Control",Type)) %>% 
#          filter(Type %in% c("CV", "Control"), !(sampleID_qPCR %in% c("R1_J1_Pre","R6_J6_Pre"))) %>%
#          select(SampleID, sampleID_qPCR, introRich, RichLevel, Rep, Inhibitory, one_of(allOTUs)) %>% 
#          mutate(Rep = factor(Rep, levels = c("Rep1","Rep2","Rep3","Rep4","Rep5","Rep6","Rep7","Rep8","Rep9","Rep10","Rep11","Rep12","Rep13"))) %>%
#          pivot_longer(cols = one_of(allOTUs), names_to="ASVID", values_to = "rawReads") %>%
#          ggplot() + geom_bar(aes(x=sampleID_qPCR, y = rawReads, fill=ASVID), stat="identity", show.legend = FALSE) + 
#          facet_grid(. ~ RichLevel, drop = TRUE, scales="free_x", space = "free_x") +
#          theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))
#        )

##### Now, let's look at abundance of OTUs #####

# First, plot qPCR against introduced richness
full_meta %>%
  ggplot(aes(x=introRich, y=qPCR_bact, col = Inhibitory)) + geom_point() +scale_y_log10()

# Combine with basic metadata
full_asv_and_meta <- otu_real_abund_round %>% as.data.frame() %>% rownames_to_column(var = "sampleID_qPCR") %>%
  mutate(sums = rowSums(across(where(is.numeric)))) %>% filter(!is.na(sums), sums>0) %>%
  select(-sums) %>% left_join(full_meta) %>%   mutate(Control = ifelse(Inhibitory=="Con", TRUE, FALSE)) 

## Plotting histograms of reads by richness of treatment; there is clearly a lump of possible contaminant DNA.
gg_histreads <- full_asv_and_meta %>%
  mutate(Control=ifelse(Control, "Control", "Real Sample")) %>%
  mutate(introRich = factor(introRich)) %>%
  pivot_longer(cols= all_of(allOTUs), names_to = "ASV", values_to = "reads") %>%
  ggplot(aes(x=reads, fill = introRich)) + geom_histogram(bins = 70) + scale_y_log10() + scale_x_log10() + facet_grid(Control~.) +
  xlab("qPCR-adjusted number of reads") + ylab("Count") + labs(fill="Number of isolates\nintroduced")
gg_histreads
ggsave(filename = "01_seq_data_cleaning/histogram_of_all_reads.png"
       , height=4, width=6, gg_histreads)


# Identify singletons
singletonDF <- data.frame(sampleCount = colSums(full_asv_and_meta[,allOTUs]>0)) %>% rownames_to_column(var = "ASVID")
gg_dotsingles <- full_asv_and_meta %>%
  pivot_longer(cols= all_of(allOTUs), names_to = "ASVID", values_to = "reads") %>%
  left_join(singletonDF) %>% mutate(singleton = ifelse(sampleCount==1, TRUE, FALSE)) %>%
  ggplot(aes(x=reads, y = factor(introRich), col = singleton)) + geom_jitter(height=0.2, width=0) + scale_x_log10()  +geom_vline(aes(xintercept = 1000)) +
  labs(col="Is singleton?") +xlab("qPCR_adjusted number of reads") + ylab("Number of introduced isolates")
gg_dotsingles
ggsave(filename = "01_seq_data_cleaning/dotplot_singletons.png", height=4, width=6, gg_dotsingles)
## Could choose 1000 to be a pretty good cutoff for what's considered "real" bacteria; and what is likely contamination
# However, I'm actually par down the OTU list by setting a global cutoff; and then looking at individual OTUs

###### 2-WAY THRESHOLD/KMEANS TEST: for Pre #####

# Use k-means thresholding to divide "real" against "fake"
allReps <- unique(full_asv_and_meta$Rep)
dir.create("01_seq_data_cleaning/byRep_ASV_thresholds")
toKeepBySampleAll <- tibble()
toKeepBySampleAll2 <-tibble()
for ( r in allReps) {
  # r="Rep10"
  subMeta <- full_asv_and_meta%>% filter(Rep == r)
  tempSamples <- subMeta %>% filter(Inhibitory !="Con") %>% pull(sampleID_qPCR)
  toKeepByRep <- tibble()
  toKeepByRep2 <- tibble()
  #If any treatment has a lower maxread than the control, then trust kmeans and not threshold
  conMaxReads <- max(subMeta %>% filter(Inhibitory == "Con") %>% select(one_of(allOTUs)))
  for (s in tempSamples) {
    # s = "R10_J21_Pre"
    # s= tempSamples[1]
    # Filter out zeros
    reads_temp <- subMeta  %>% filter(sampleID_qPCR == s) %>%
      pivot_longer(cols = one_of(allOTUs), names_to = "ASVID", values_to = "Reads") %>% filter(Reads > 0) %>%
      mutate(logReads = log10(Reads)) %>% select(sampleID_qPCR, ASVID, logReads, Reads)
    tempRich <- subMeta  %>% filter(sampleID_qPCR == s) %>% select(introRich)
    # If less than two; take the maximum one
    if (nrow(reads_temp)<3) {
      toKeep <- which.max(reads_temp$logReads)
    } else {
      # Cluster
      kmeans_temp <- kmeans(reads_temp$logReads, centers = 2)
      # Get the larger group; compare to threshold of 1000
      toKeep_kmeans <- c(which(kmeans_temp$cluster ==  which.max(kmeans_temp$centers)))
      toKeep_1000 <- which(reads_temp$Reads>1000)

      # For each dot within a range, find out if it's close to a farther one or not.
      allVal <- sort(reads_temp$Reads)
      # Get points inside and outside of test range
      if ( any(allVal<500)) {
        below <- max(allVal[allVal<500])
      } else { below = min(allVal) }
      if ( any(allVal>2000) ) {
        above <- min(allVal[allVal > 2000])
      } else {
        above = max(allVal)
      }
      between <- allVal[allVal>500 & allVal <2000]
      
      if (length(between)==0) {
        cutoff <- 1000
      } else {
        for (i in between) {
          direction <- apply(cbind(c(below,between,above)-below, above - c(below,between,above)), MARGIN=1, which.min)
          cutoff <- mean(c(c(below,between,above)[which(diff(direction)==1)], c(below,between,above)[which(diff(direction)==1)+1]))
        }
        # toKeep_cutoff1 <- which(reads_temp$Reads>cutoff)
        # 
        # # Test against a straight cutoff to see which is better
        # diff_cutoff <- tempRich - length(toKeep_cutoff1)
        # diff_1000 <- tempRich - length(toKeep_1000)
        # if (abs(diff_1000) == abs(diff_cutoff)) {
        #   ind_best_cutoff <- which.max(c(diff_1000, diff_cutoff))
        # } else {
        #   ind_best_cutoff <- which.min(c(abs(diff_1000),abs(diff_cutoff)))
        # }
        # cutoff <- c(1000, cutoff)[ind_best_cutoff]
      }
      toKeep_cutoff <- which(reads_temp$Reads>cutoff)
    
      # which is closest to target?
      diff_kmeans <- tempRich - length(toKeep_kmeans)
      diff_cutoff <- tempRich - length(toKeep_cutoff)
      # If one of them is -1 and the other is +1; then use the +1
      ind_best <- NA
      if (abs(diff_kmeans) == abs(diff_cutoff)) {
        ind_best <- which.max(c(diff_kmeans, diff_cutoff))
      } else {
        ind_best <- which.min(c(abs(diff_kmeans),abs(diff_cutoff)))
      }
      toKeep <- unlist(list(toKeep_kmeans,toKeep_cutoff)[ind_best])
    }
    # If the maxread of this sample was less than the control, then use kmeans. Probably qPCR messed up.
    if (conMaxReads > max(reads_temp$Reads)) {
      toKeep <- toKeep_kmeans
    }
    # add into other sample data tibble
    toKeepByRep <- rbind(toKeepByRep, reads_temp[toKeep,-3])
    toKeepByRep2 <- rbind(toKeepByRep2, reads_temp[toKeep_kmeans,-3])
    
  }
  # For plotting; highlighting
  toKeepByRep$KEEP <- TRUE
  toKeepByRep2$KEEP <- TRUE
  toKeepByRep2 <- toKeepByRep2 %>% filter(!is.na(ASVID))
  
  ggsave(paste0("01_seq_data_cleaning/byRep_ASV_thresholds/",r,"_adaptivethresholds.png"), width = 5, height = 3,
         subMeta %>%
           pivot_longer(cols = one_of(allOTUs), names_to = "ASVID", values_to = "Reads") %>% 
           select(sampleID_qPCR, introRich, Inhibitory, ASVID, Reads) %>%
           left_join(toKeepByRep) %>% mutate(KEEP = ifelse(is.na(KEEP), FALSE, KEEP)) %>%
           # ggplot(aes(y = factor(RichLevel), x = Reads, col = Inhibitory)) +
           ggplot(aes(y = factor(introRich), x = Reads, col = sampleID_qPCR, pch = KEEP)) + facet_grid(Inhibitory~., scales = "free_y", drop = TRUE) +
           geom_jitter(width = 0, height=0.2) + geom_vline(aes(xintercept=500)) +geom_vline(aes(xintercept=2000)) +  scale_x_log10() +
           scale_shape_manual(values =c(21,19))
  )
  ggsave(paste0("01_seq_data_cleaning/byRep_ASV_thresholds/",r,"_kmeansthresholds.png"), width = 5, height = 3,
         subMeta %>%
           pivot_longer(cols = one_of(allOTUs), names_to = "ASVID", values_to = "Reads") %>% 
           select(sampleID_qPCR, introRich, Inhibitory, ASVID, Reads) %>%
           left_join(toKeepByRep2) %>% mutate(KEEP = ifelse(is.na(KEEP), FALSE, KEEP)) %>%
           # ggplot(aes(y = factor(RichLevel), x = Reads, col = Inhibitory)) +
           ggplot(aes(y = factor(introRich), x = Reads, col = sampleID_qPCR, pch = KEEP)) + facet_grid(Inhibitory~., scales = "free_y", drop = TRUE) +
           geom_jitter(width = 0, height=0.2) + geom_vline(aes(xintercept=500)) + geom_vline(aes(xintercept=2000))+scale_x_log10() +
           scale_shape_manual(values =c(21,19))
  )
  # Merge with rest of dataset
  toKeepBySampleAll <- rbind(toKeepBySampleAll, toKeepByRep)
  toKeepBySampleAll2 <- rbind(toKeepBySampleAll2, toKeepByRep2)
  
}
### I'M GOING TO USE K-MEANS THRESHOLDING HERE; MIGHT GO BACK AND CHANGE LATER ###
# reduce duplicates to get list of OTUs
filtOTUs <- unique(toKeepBySampleAll$ASVID)
# How many reads did we actually lose?
1-sum(otu_raw[,filtOTUs])/sum(otu_raw[,allOTUs])
# 1.8% of reads; not too bad.


# Make new table
asv_table_filt <- toKeepBySampleAll %>% select(-KEEP) %>% 
  pivot_wider(names_from = ASVID, values_from = Reads) %>%
  full_join(full_meta %>% filter(Type == "CV") %>%select(sampleID_qPCR)) %>% select(sampleID_qPCR, everything())
asv_table_filt[is.na(asv_table_filt)] <- 0
full_asv_and_meta_filt1 <-  full_meta %>% filter(Type =="CV") %>%
  full_join(asv_table_filt) %>% 
  mutate(obsRich = rowSums(across(one_of(filtOTUs))>0)) %>%
  select(sampleID_qPCR, Inhibitory, introRich, obsRich, Rep, one_of(filtOTUs), everything()) 

#### Are remaining ASVs singletons/rare? ####
## Here, we update otu tables and calcuate the prevalence and rareity of these ASVs before filtering
# If they were singletons or extremely rare in the "pre" dataset, that suggests we might have a problem with contamination
# Okay, now let's look for additional contaminants
otu_raw_filt <- otu_raw %>% select(sampleID_qPCR, one_of(filtOTUs)) %>% 
  pivot_longer(cols=one_of(filtOTUs), names_to="ASVID", values_to = "rawReads") %>%
  right_join(toKeepBySampleAll) %>% select(sampleID_qPCR, ASVID, rawReads) %>%
  pivot_wider(names_from = ASVID, values_from=rawReads) %>% as.data.frame()
otu_raw_filt[is.na(otu_raw_filt)] <- 0

otu_prop_filt <- otu_prop %>% select(sampleID_qPCR, one_of(filtOTUs)) %>% 
  pivot_longer(cols=one_of(filtOTUs), names_to="ASVID", values_to = "rawReads") %>%
  right_join(toKeepBySampleAll) %>% select(sampleID_qPCR, ASVID, rawReads) %>%
  pivot_wider(names_from = ASVID, values_from=rawReads) %>% as.data.frame()
otu_prop_filt[is.na(otu_prop_filt)] <- 0

otu_freq_filt <- otu_freq %>% select(sampleID_qPCR, one_of(filtOTUs)) %>% 
  pivot_longer(cols=one_of(filtOTUs), names_to="ASVID", values_to = "rawReads") %>%
  right_join(toKeepBySampleAll) %>% select(sampleID_qPCR, ASVID, rawReads) %>%
  pivot_wider(names_from = ASVID, values_from=rawReads) %>% as.data.frame()
otu_freq_filt[is.na(otu_freq_filt)] <- 0

otu_pa_filt <- otu_pa %>% select(sampleID_qPCR, one_of(filtOTUs)) %>% 
  pivot_longer(cols=one_of(filtOTUs), names_to="ASVID", values_to = "rawReads") %>%
  right_join(toKeepBySampleAll) %>% select(sampleID_qPCR, ASVID, rawReads) %>%
  pivot_wider(names_from = ASVID, values_from=rawReads) %>% as.data.frame()
otu_pa_filt[is.na(otu_pa_filt)] <- 0
# Any ASVs that have a count less than 100
names(which(colSums(otu_raw_filt[,filtOTUs])<100))
# Any remaining singletons?
any(otu_raw_filt[,filtOTUs] <2 & otu_raw_filt[,filtOTUs] !=0)
# What about less than 10 reads originally?
sum(otu_raw_filt[,filtOTUs] < 10 &otu_raw_filt[,filtOTUs] !=0)

## Flagging singletons and low-abundance ASVs
# What about less than 1% originally?
asv_percL1 <- names(which(colSums(otu_prop_filt[,filtOTUs] < 0.01 & otu_prop_filt[,filtOTUs,] !=0)>0))
# Occurance less than one?
asv_prevL1 <- names(which(colSums(otu_pa_filt[,filtOTUs]) ==1))

# Samples that this occurs in
# What about less than 1% originally?
sample_percL1 <- otu_prop_filt[which(rowSums(otu_prop_filt[,asv_percL1] < 0.01 & otu_prop_filt[,asv_percL1,] !=0)>0),"sampleID_qPCR"]
# Occurance less than one?
sample_prevL1 <- otu_pa_filt[which(rowSums(otu_pa_filt[,asv_prevL1])>0),"sampleID_qPCR"]

### Looking at presence of low-proportion ASVs
full_asv_and_meta_filt1 %>%
  filter(sampleID_qPCR %in% sample_percL1) %>%
  pivot_longer(cols = one_of(filtOTUs), names_to = "ASVID", values_to = "Reads") %>% 
  select(sampleID_qPCR, introRich, Inhibitory, ASVID, Reads) %>%
  rowwise() %>% mutate(prevone = ifelse(ASVID %in% asv_percL1, TRUE, FALSE)) %>%
  ungroup() %>% 
  # ggplot(aes(y = factor(RichLevel), x = Reads, col = Inhibitory)) +
  ggplot(aes(y = sampleID_qPCR, x = Reads, col = factor(introRich), pch = prevone)) + facet_grid(Inhibitory~., scales = "free_y", drop = TRUE) +
  geom_jitter(width = 0, height=0.2) + geom_vline(aes(xintercept=1000))+scale_x_log10() +
  scale_shape_manual(values =c(21,19))
## Thelow prop ASVs all look like they SHOULD belong; they're all in the correct clusters.


### Looking at presence of low-prevalence ASVs
full_asv_and_meta_filt1 %>%
  filter(sampleID_qPCR %in% sample_prevL1) %>%
  pivot_longer(cols = one_of(filtOTUs), names_to = "ASVID", values_to = "Reads") %>% 
  select(sampleID_qPCR, introRich, Inhibitory, ASVID, Reads) %>%
  rowwise() %>% mutate(prevone = ifelse(ASVID %in% asv_prevL1, TRUE, FALSE)) %>%
  ungroup() %>% 
  # ggplot(aes(y = factor(RichLevel), x = Reads, col = Inhibitory)) +
  ggplot(aes(y = sampleID_qPCR, x = Reads, col = factor(introRich), pch = prevone)) + facet_grid(Inhibitory~., scales = "free_y", drop = TRUE) +
  geom_jitter(width = 0, height=0.2) + geom_vline(aes(xintercept=1000))+scale_x_log10() +
  scale_shape_manual(values =c(21,19))
## The single prevalence ASVs all look like they SHOULD belong; they're all in the correct clusters.

#### Finally, let's look at ASVs that are in both
# Things that are rare AND low-prevalence; let's cross-reference them with controls to see if it's
# PCR contamination or extraction contamination of some sort
inBoth <- full_asv_and_meta_filt1 %>% pivot_longer(cols = one_of(filtOTUs), names_to = "ASVID", values_to = "Reads") %>%
  filter(Reads>0) %>% select(ASVID, Inhibitory) %>% distinct() %>% arrange(ASVID) %>%
  mutate(dup = duplicated(ASVID)) %>% filter(dup) %>% pull(ASVID)

# Cross-ref this with OTUs that are abundant in controls
# Get OTUs that were in controls
exp_control_only <- otu[,grep("H2O|SWAB", colnames(otu))] %>%  mutate(total = rowSums(.)) %>% filter(total>0) %>% select(-total) %>%
  apply(MARGIN = 1, max) %>% as.data.frame() %>% rename(maxReads_control = paste0(".")) %>% rownames_to_column(var = "ASVID")
extr_control_only <- otu[,grep("CON", colnames(otu))] %>%  mutate(total = rowSums(.)) %>% filter(total>0) %>% select(-total) %>%
  apply(MARGIN = 1, max) %>% as.data.frame() %>% rename(maxReads_control = paste0(".")) %>% rownames_to_column(var = "ASVID")
pcr_control_only <- otu[,grep("NTC", colnames(otu))] %>%  mutate(total = rowSums(.)) %>% filter(total>0) %>% select(-total) %>%
  apply(MARGIN = 1, max) %>% as.data.frame() %>% rename(maxReads_control = paste0(".")) %>% rownames_to_column(var = "ASVID")
# 
# toRemove_abundantInControls <- otu[,-grep("CON|H2O|SWAB|NTC", colnames(otu))] %>% 
#   apply(MARGIN = 1, max) %>% as.data.frame() %>% rename(maxReads_samples = paste0(".")) %>% rownames_to_column(var = "ASVID") %>%
#   right_join(control_only) %>%
#   mutate(orderMag_sample = log10(maxReads_samples+1), orderMag_control = log10(maxReads_control+1), toRemove = (orderMag_sample - orderMag_control) < 1) %>%
#   filter(toRemove) %>% pull(ASVID)

exp_control_only %>% rowwise() %>% mutate(inBoth = ifelse(ASVID %in% inBoth, ASVID, "Other") ) %>%
  ggplot() + geom_jitter(aes(x=maxReads_control, y=1,col =inBoth), width = 0) + scale_x_log10()
# Looks like these are just as abundant as any of the other "control" contaminants-- so likely not GROWTH contamination

extr_control_only %>% rowwise() %>% mutate(inBoth = ifelse(ASVID %in% inBoth, ASVID, "Other") ) %>%
  ggplot() + geom_jitter(aes(x=maxReads_control, y=1,col =inBoth), width = 0) + scale_x_log10()
# Looks like these are just as abundant as any other "controL' contaminants-- so likely not EXTRACTION contamination

pcr_control_only %>% rowwise() %>% mutate(inBoth = ifelse(ASVID %in% inBoth, ASVID, "Other") ) %>%
  ggplot() + geom_jitter(aes(x=maxReads_control, y=1,col =inBoth), width = 0) + scale_x_log10()
# Same as above; nothing looks like a signficant outlier.

# Overlaps between inhibitory and controls
dir.create("01_seq_data_cleaning/overlapping_asvs_inhib_non")
toRemove_from <- tibble()
universal_contam <- tibble()
for ( b in inBoth) {
  # b = "5648dccee530d68ceb3e4d7d22cf8756"
  ratio_temp <- full_asv_and_meta_filt1 %>% rowwise() %>%
    mutate(maxCount = max(across(one_of(filtOTUs)))) %>% ungroup() %>%
    select(sampleID_qPCR, Inhibitory, introRich, Rep, paste0(b), maxCount) %>%
    filter(maxCount>0, get(b)>0) %>%
    mutate(ratio = get(b)/maxCount) 
  maxRatioInhib <- ratio_temp %>% filter(Inhibitory == "Inhib") %>% pull(ratio) %>%max()
  maxRatioNon <- ratio_temp %>% filter(Inhibitory == "Non") %>% pull(ratio) %>%max()
  
  if ( !any(table(ratio_temp[,c("Inhibitory")])<=1) | (maxRatioInhib>0.5 & maxRatioNon>0.5) ) { # check if one side too short to do t-test
    if ( (maxRatioInhib>0.5 & maxRatioNon>0.5)) {
      p_temp <- 1
    } else {
      if (!any(table(ratio_temp[,c("Inhibitory")])<=1) ) {
        p_temp <- t.test(ratio ~ Inhibitory, data=ratio_temp)$p.value
      } else {
        p_temp <- 1
      }
    }
   
  } else {
    p_temp <- 0
  }
  # Get info about isolate
  tax_temp <- taxa %>% filter(ASVID == b) %>% select(ASVID, ScientificName_unique)
  inhibornot_temp <- inhibMeta %>% filter(ASVID == b) %>% select(ASVID, inhibitory_probability_alex, inhibitory_probability_woodhams, inhibitory_probability_tree)
  
  if (p_temp < 0.05) {
    inhibMean <- ratio_temp %>% filter(Inhibitory == "Inhib") %>% select(ratio) %>% colMeans()
    nonMean <- ratio_temp %>% filter(Inhibitory == "Non") %>% select(ratio) %>% colMeans()
    toFilter <- c("Inhib","Non")[which.min(c(inhibMean, nonMean))]
    temp_meta <- ratio_temp %>% filter(Inhibitory == toFilter) %>% select(sampleID_qPCR, paste0(b)) %>% 
      pivot_longer(cols = one_of(b), names_to = "ASVID", values_to = "Reads") %>% left_join(tax_temp) %>% left_join(inhibornot_temp)
    toRemove_from <- rbind(toRemove_from, temp_meta)
  } else {
    universal_contam <- rbind(universal_contam, tax_temp %>% full_join(inhibornot_temp))
  }
  ggsave(paste0("01_seq_data_cleaning/overlapping_asvs_inhib_non/", b,".png"), height=3, width = 8
    ,full_asv_and_meta_filt1 %>% rowwise() %>%
      mutate(maxCount = max(across(one_of(filtOTUs)))) %>% ungroup() %>%
      select(sampleID_qPCR, Inhibitory, introRich, Rep, paste0(b), maxCount) %>%
      pivot_longer(cols = one_of(c(b,"maxCount")), names_to = "ASVID", values_to = "Reads") %>%
      filter(Reads>0) %>%
      ggplot() + geom_bar(aes(x=introRich, y= Reads, fill=ASVID), stat = "identity", position="dodge2", width = 2) +
      facet_grid(Inhibitory~Rep) + theme(axis.text.x = element_text(angle=90))
  )
  
}

# Find IDs of universal contams
universal_contam
# toRemove_from %>% View()



## OLD (for adaptive thresholds) Manual override notes and decisions:
# ASVID 5648dccee530d68ceb3e4d7d22cf8756 is a Pseudomonas that is predicted to be inhibitory, but doesn't match anything from the databases.
# It is abundant in many different samples, so it is "real"-- I will need to keep this in all samples.

# ASVID 63afe8e6aac58bf0d670a82ca5bc574c is a Yersiniaceae that matches isolates in Alex and Woodhams. It is predicted to be NON inhibitory
# It is also mostly seen in non-samples; all Inhib
# Consider removing samples Rep5-0 and Rep 5-1 entirely due to contamination?

# ASVID 82dece6e35540738ba450a0c3a90b5a0 is Serratia, and should be kept in all "Non", remove from all inhib

# ASVID dcba105f35d8ebc9e22269c7491ad3a7 is a Stenotrophomonas that matches things in Alex and Woodhams. It is predicted to be inhibitory.
# Strangely, it's almost only found in non samples. InhibRep 1,5,7 is contaminated with it. Keep in all samples; create flag.

# ASVID 38c27ceaed634984c1225a82648cf571 is a Pseudomonas that matches only Woodhams; predicted to be inhibitroy.
# It's in R8_J7 but actually pretty low; so we're going to remove from this sample only.

# ASVID fb67b286b0f781b0de13d50179318995 is also a Stenotrophonmonas that matches Alex and Woodhams. It is predicted to be inhibitory.
# It's actually pretty low in abundance in Non; low in Rep 8. Also, it is ONLY found in Rep 7/8. 
# It is in both Rich 1 and 3 in Rep 7 and 8; could be contaminant so I'll remove it from Rep 7

### I will need to think more carefully about how to deal with these contaminants 
## Will consider them after assigning specific identities to each ASV.

############# Plotting some contaminants in extraction controls ###########
dir.create("./01_seq_data_cleaning/contaminant_plots_extractionControls")
for ( asv in extr_control_only$ASVID) {
  tempTaxa <-taxa %>% filter(ASVID == asv) %>% pull(ScientificName_unique)
    temp_Data <- otu_prop %>% select(sampleID_qPCR, one_of(asv)) %>% filter(!is.na(paste0(asv))) %>% left_join(full_asv_and_meta%>% select(sampleID_qPCR, qPCR_bact))  %>%
      rename(Reads=paste0(asv))
  # temp_Data <- full_asv_and_meta %>% select(sampleID_qPCR, one_of(asv), qPCR_bact) %>% rename(Reads = paste0(asv)) 
  if ( sum(temp_Data$Reads) > 0 ) {
    tempP <- coef(summary(lm(Reads ~ qPCR_bact, data=temp_Data)))[2,4] # Get pvalue
    ggsave(filename = paste0("01_seq_data_cleaning/contaminant_plots_extractionControls/", asv, "_",tempTaxa, ".png"),
           ggplot(temp_Data, aes(x=qPCR_bact, y=Reads)) + geom_point() + ggtitle(paste0("p = ", tempP))
    )
  }
}
dir.create("./01_seq_data_cleaning/contaminant_plots_pcrControls")
for ( asv in pcr_control_only$ASVID) {
  tempTaxa <-taxa %>% filter(ASVID == asv) %>% pull(ScientificName_unique)
     temp_Data <- otu_prop %>% select(sampleID_qPCR, one_of(asv)) %>% filter(!is.na(paste0(asv))) %>% left_join(full_asv_and_meta%>% select(sampleID_qPCR, qPCR_bact))  %>%
      rename(Reads=paste0(asv))
  # temp_Data <- full_asv_and_meta %>% select(sampleID_qPCR, one_of(asv), qPCR_bact) %>% rename(Reads = paste0(asv)) 
  if ( sum(temp_Data$Reads) > 0 ) {
    tempP <- coef(summary(lm(Reads ~ qPCR_bact, data=temp_Data)))[2,4] # Get pvalue
    ggsave(filename = paste0("01_seq_data_cleaning/contaminant_plots_pcrControls/", asv, "_",tempTaxa, ".png"),
           ggplot(temp_Data, aes(x=qPCR_bact, y=Reads)) + geom_point() + ggtitle(paste0("p = ", tempP))
    )
  }
}
# Things to remove
cutoff_remove <- allOTUs[!(allOTUs %in% filtOTUs)]
contaminants <- c()
allToRemove_unique <- unique(c(cutoff_remove, contaminants))
# get all things to remove; unique
allToRemove_unique_toSave <- data.frame(ASVID=allToRemove_unique) %>% left_join(taxa) %>% select(ASVID, Taxon, ScientificName_unique)
write.table(allToRemove_unique_toSave, file = "01_seq_data_cleaning/removed_ASVs_from_dataset.txt", quote=FALSE, row.names=FALSE, sep = "\t")

### Calculate proportion of sequences lost
1-sum(otu_raw[,filtOTUs])/sum(otu_raw[,allOTUs])
# 1.84%! Nice.
length(filtOTUs) # 48 sequences remaining
length(allToRemove_unique) # 849 sequences removed due to low abundance


##### Plotting Alex's trees #######
## Check if we can cluster any OTUs together, because they are highly correlated or inversely so
tree.filt1 <- keep.tip(tree, filtOTUs)

# TLet's see what their identities are, and if they cluster together with other OTUs?
tree.filt1_relablled <- tree.filt1
tree.filt1_relablled$tip.label <- taxa[match(tree.filt1_relablled$tip.label, taxa$ASVID), "ScientificName_unique"]

png(filename = "01_seq_data_cleaning/tree_myisolates_filtered.png",  height=600, width=480)
plot(tree.filt1_relablled)
dev.off()

# Collapse tips that are monophyletic; and whose taxaIDs are of the same taxonomic group.
# If there are multiple levels in one monophyletic cluster (e.g. g, g, o) then collapse at o, as long as all tips belonging to
# that group are within the order (for example). 

taxa.filt1 <- taxa %>% filter(ASVID %in% tree.filt1$tip.label)

# Sort tips from outer to inner edges
rootNode <-   length(tree.filt1$tip.label) + 1
allTips <- tree.filt1$tip.label
node_path_length <- sapply(1:length(allTips), FUN = function(x) length(nodepath(tree.filt1, from = rootNode, to = x)))
# Order by this path length
remainingTips <- allTips[order(node_path_length)]

finalTips <- list()
while (length(remainingTips)>0) {
  currTipNodes <- match(remainingTips[1], tree.filt1$tip.label)
  currPath <- rev(nodepath(tree.filt1, from =  rootNode, to = currTipNodes))
  deepest_monophyl_found = FALSE
  currNodeTestn = 2 # Start at second node; first node is lone node
  while ( !deepest_monophyl_found ) {
    currNodeTest = currPath[currNodeTestn]
    tempdesc <- getDescendants(tree.filt1, currNodeTest)
    # Get rid of nodes; keep just tips
    tempdesc <- tempdesc[-which(tempdesc>rootNode)] 
    # get IDs of tips
    currTaxa <- taxa.filt1 %>% filter(ASVID %in% tree.filt1$tip.label[tempdesc])
    found_monophyly_or_not <- FALSE
    lvls <- c("Species","Genus","Family","Order","Class","Phylum","Domain")
    lvl <- 1
    while ( !found_monophyly_or_not ) {
      currlvl <- lvls[lvl]
      allUnique <- length(unique(currTaxa[,currlvl])[!is.na(unique(currTaxa[,currlvl]))])<=1
      if (allUnique) {
        allPresent <- !any(is.na(currTaxa[,currlvl]))
        if (allPresent) {
          #### MONOPHYLETIC GROUP FOUND ####
          found_monophyly_or_not <- TRUE
          currNodeTestn <- currNodeTestn + 1
        } else {
          lvl= lvl + 1
        }
      } else {
        ##### NOT MONOPHYLETIC GROUP ####
        found_monophyly_or_not <- TRUE
        deepest_monophyl_found <- TRUE
        # finalTips[[length(finalTips) + 1]] <- remainingTips[1]
        # # remove this tip
        # remainingTips <- remainingTips[-1]
      }
    }
  }
  monophyletic_group <- getDescendants(tree.filt1, currPath[currNodeTestn-1]) # Get the previous node, since all others would have ended up not mono
  monophyletic_tips <- tree.filt1$tip.label[monophyletic_group[monophyletic_group<rootNode]]
  monophyletic_tips_noNAs <- match(monophyletic_tips, remainingTips)[!is.na(match(monophyletic_tips, remainingTips))]
  
  nested_tips <- unlist(lapply(finalTips, FUN = function(x) any(!is.na(match(monophyletic_tips, x)))) )
  toReplace <- ifelse(length(nested_tips)>0, which(nested_tips), NA)
  if (length(toReplace)>0 & any(!is.na(toReplace))) {
    finalTips[[toReplace]] <- monophyletic_tips
  } else {
    finalTips[[length(finalTips)+1]] <- remainingTips[monophyletic_tips_noNAs]
    
  }

  remainingTips <- remainingTips[-monophyletic_tips_noNAs]
  
}
##### See whether these samples that are grouped correlated with each other?

### Maybe check Alex's isolates to see if they match up to anything

alexTips <- alex_isos$ASVID[alex_isos$ASVID %in% tree$tip.label]

tree.filt_walex <- keep.tip(tree, c(filtOTUs, alexTips))
# Turn into color
colorTips <- rep("black", length(tree.filt_walex$tip.label))
alexIsos_inhib <- inhibMeta  %>% filter(dataset == "alex", knownInhibitory == "Inhib") %>% pull(ASVID)
alexIsos_non <- inhibMeta  %>% filter(dataset == "alex", knownInhibitory == "Non") %>% pull(ASVID)
# View(inhibMeta)
colorTips[ match(alexIsos_inhib, tree.filt_walex$tip.label) ] <- "red"
colorTips[ match(alexIsos_non, tree.filt_walex$tip.label) ] <- "blue"

# colorTips[grep(".ab1", tree.filt_walex$tip.label)] <- "red"
taxaNames <- taxa[match(tree.filt_walex$tip.label, taxa$ASVID), "ScientificName_unique"]
combined_names <- paste(taxa[match(tree.filt_walex$tip.label, taxa$ASVID),"ScientificName_unique"]
                         , taxa[match(tree.filt_walex$tip.label, taxa$ASVID),"ASVID"], sep ="--")
taxaNames[grep(".ab1", tree.filt_walex$tip.label)] <- combined_names[grep(".ab1", tree.filt_walex$tip.label)]

tree.filt_walex$tip.label <- taxaNames

png(file = "01_seq_data_cleaning/tree_filtered_ASVS_withAlexIsolates.png",height=1000, width=800)
plot(tree.filt_walex, tip.color = colorTips)
dev.off()

# Get inhibitory data
inhibColRAmp <- colorRampPalette(colors = c("blue","grey","red"))
colRange <- inhibColRAmp(n = 101)

colorTips <- inhibMeta %>% filter(ASVID %in% filtOTUs) %>%
  mutate(inhibitory = ifelse(!is.na(inhibitory_probability_alex), inhibitory_probability_alex, ifelse(!is.na(inhibitory_probability_woodhams), inhibitory_probability_woodhams, inhibitory_probability_tree))) %>%
  select(ASVID, inhibitory) %>% arrange(tree.filt1$tip.label) %>%
  mutate(colorTip = colRange[round(as.numeric(inhibitory)*100)+1]
         , colorTip = ifelse(is.na(colorTip), "black", colorTip)) %>%
  pull(colorTip)

png("01_seq_data_cleaning/myIsolates_predicted_inhibitory.png", height=600)
plot(tree.filt1_relablled, tip.color = colorTips)
dev.off()

########## 2-WAY THRESHOLD/KMEANS TEST: for Post #########

# Get full asv and meta for post-bd
full_asv_and_meta_POST <- t(otu_relabund) %>% as.data.frame() %>%rownames_to_column(var="sampleID_qPCR") %>% 
  mutate(PropSums = rowSums(across(where(is.numeric)))) %>% filter(!is.na(PropSums), PropSums>0) %>% select(-PropSums) %>%
  left_join(full_meta) %>%   mutate(Control = ifelse(Inhibitory=="Con", TRUE, FALSE)) %>%
  filter(Type == "Bd") %>% select(-qPCR_bact)

## Basic cutoff

## Plotting histograms of reads by richness of treatment; there is clearly a lump of possible contaminant DNA.
gg_histreads_post <- full_asv_and_meta_POST %>%
  mutate(Control = ifelse(Control, "Control","Real sample")) %>%
  pivot_longer(cols= all_of(allOTUs), names_to = "ASV", values_to = "reads") %>%
  ggplot(aes(x=reads, fill = factor(introRich))) + geom_histogram(bins = 70) + 
  scale_y_log10() + scale_x_log10() + facet_grid(Control~.) +
  geom_vline(aes(xintercept = 0.03))+
  xlab("Relative abundance of reads") + ylab("Count") + labs(fill="Number of introduced\nisolates")
gg_histreads_post
# Seems like 0.3 is a good sweet spot for a cutoff
ggsave(filename = "01_seq_data_cleaning/histogram_of_all_reads_post.png", height=4, width=6, 
       gg_histreads_post)

full_asv_and_meta_POST %>%
  pivot_longer(cols= all_of(allOTUs), names_to = "ASVID", values_to = "reads") %>%
  ggplot(aes(x=reads, y = factor(introRich))) + geom_jitter(height=0.2, width=0) + scale_x_log10()  +geom_vline(aes(xintercept = 0.03)) 

# Verify that controls look okay still
full_asv_and_meta_POST_raw <- t(otu) %>% as.data.frame() %>%rownames_to_column(var="sampleID_qPCR") %>% 
  mutate(PropSums = rowSums(across(where(is.numeric)))) %>% filter(!is.na(PropSums), PropSums>0) %>% select(-PropSums) %>%
  left_join(full_meta) %>%   mutate(Control = ifelse(Inhibitory=="Con", TRUE, FALSE)) %>%
  filter(Type == "Bd") %>% select(-qPCR_bact)
full_asv_and_meta_POST_raw %>% select(sampleID_qPCR, Inhibitory, Rep, RichLevel, one_of(allOTUs)) %>%
  pivot_longer(cols = one_of(allOTUs), names_to="ASVID", values_to = "reads") %>%
  ggplot() + geom_bar(aes(x=RichLevel, y=reads, fill = ASVID), stat="identity", show.legend = FALSE) +
  facet_grid(Inhibitory~Rep, scales = "free_x", space = "free_x")
# Whater are those things in the controls?
controlSums <- full_asv_and_meta_POST_raw %>% filter(Control) %>% select(one_of(allOTUs)) %>% colSums()
full_asv_and_meta_POST_raw %>% filter(Control) %>% select(sampleID_qPCR, Inhibitory, Rep, RichLevel, one_of(allOTUs)) %>%
  pivot_longer(cols = one_of(allOTUs), names_to="ASVID", values_to = "reads") %>%
  filter(reads>1000) %>% left_join(taxa) %>%
  ggplot() + geom_bar(aes(x=RichLevel, y=reads, fill = ScientificName_unique), stat="identity", show.legend = TRUE) +
  facet_grid(~Rep, scales = "free_x", space = "free_x")
ggsave("01_seq_data_cleaning/experiment_controls_postbd_taxasummary.png", height=5, width=8,
       full_asv_and_meta_POST_raw %>% filter(Control) %>% select(sampleID_qPCR, Inhibitory, Rep, RichLevel, one_of(allOTUs)) %>%
         pivot_longer(cols = one_of(allOTUs), names_to="ASVID", values_to = "reads") %>%
         filter(reads>1000) %>% left_join(taxa) %>%
         ggplot() + geom_bar(aes(x=RichLevel, y=reads, fill = ScientificName_unique), stat="identity", show.legend = TRUE) +
         facet_grid(~Rep, scales = "free_x", space = "free_x"))
full_asv_and_meta_POST_raw %>% select(sampleID_qPCR, Inhibitory, Rep, RichLevel, one_of(allOTUs)) %>%
  pivot_longer(cols = one_of(allOTUs), names_to="ASVID", values_to = "reads") %>%
  filter(reads>10000) %>% left_join(taxa) %>%
  ggplot() + geom_bar(aes(x=RichLevel, y=reads, fill = ScientificName_unique), stat="identity") +
  facet_grid(Inhibitory~Rep, scales = "free_x", space = "free_x")
taxa %>% filter(ASVID %in% names(which(controlSums>10000)))

# First, manually go through each sample and see if default 0.03 cutoff looks okay; adjust in edge cases
allReps <- unique(full_asv_and_meta_POST$Rep)
dir.create("01_seq_data_cleaning/byRep_ASV_thresholds_post")
toKeepBySampleAll_post <- tibble()
toKeepBySampleAll_post2 <-tibble()
for ( r in allReps) {
  # r="Rep1"
  subMeta <- full_asv_and_meta_POST%>% filter(Rep == r)
  tempSamples <- subMeta %>% filter(Inhibitory !="Con") %>% pull(sampleID_qPCR)
  toKeepByRep <- tibble()
  toKeepByRep2 <- tibble()
  for (s in tempSamples) {
    # s = "R1_J1_Pre"
    # s= tempSamples[1]
    # Filter out zeros
    reads_temp <- subMeta  %>% filter(sampleID_qPCR == s) %>%
      pivot_longer(cols = one_of(allOTUs), names_to = "ASVID", values_to = "Reads") %>% filter(Reads > 0) %>%
      mutate(logReads = log10(Reads)) %>% select(sampleID_qPCR, ASVID, logReads, Reads)
    tempRich <- subMeta  %>% filter(sampleID_qPCR == s) %>% select(introRich)
    # If less than two; take the one over 1000
    if (nrow(reads_temp)<3) {
      toKeep <- which.max(reads_temp$logReads)
    } else {
      # Cluster
      kmeans_temp <- kmeans(reads_temp$logReads, centers = 2)
      # Get the larger group; compare to threshold of 1000
      toKeep_kmeans <- c(which(kmeans_temp$cluster ==  which.max(kmeans_temp$centers)))
      toKeep_1000 <- which(reads_temp$Reads>0.03)
      
      # For each dot within a range, find out if it's close to a farther one or not.
      allVal <- sort(reads_temp$Reads)
      # Get points inside and outside of test range
      if ( any(allVal<0.015)) {
        below <- max(allVal[allVal<0.015])
      } else { below = min(allVal) }
      if ( any(allVal>0.06) ) {
        above <- min(allVal[allVal > 0.06])
      } else {
        above = max(allVal)
      }
      between <- allVal[allVal>0.015 & allVal <0.06]
      
      if (length(between)==0) {
        cutoff <- 0.03
      } else {
        for (i in between) {
          direction <- apply(cbind(c(below,between,above)-below, above - c(below,between,above)), MARGIN=1, which.min)
          cutoff <- mean(c(c(below,between,above)[which(diff(direction)==1)], c(below,between,above)[which(diff(direction)==1)+1]))
        }
        # toKeep_cutoff1 <- which(reads_temp$Reads>cutoff)
        # 
        # # Test against a straight cutoff to see which is better
        # diff_cutoff <- tempRich - length(toKeep_cutoff1)
        # diff_1000 <- tempRich - length(toKeep_1000)
        # if (abs(diff_1000) == abs(diff_cutoff)) {
        #   ind_best_cutoff <- which.max(c(diff_1000, diff_cutoff))
        # } else {
        #   ind_best_cutoff <- which.min(c(abs(diff_1000),abs(diff_cutoff)))
        # }
        # cutoff <- c(1000, cutoff)[ind_best_cutoff]
      }
      toKeep_cutoff <- which(reads_temp$Reads>cutoff)
      
      # which is closest to target?
      diff_kmeans <- tempRich - length(toKeep_kmeans)
      diff_cutoff <- tempRich - length(toKeep_cutoff)
      # If one of them is -1 and the other is +1; then use the +1
      ind_best <- NA
      if (abs(diff_kmeans) == abs(diff_cutoff)) {
        ind_best <- which.max(c(diff_kmeans, diff_cutoff))
      } else {
        ind_best <- which.min(c(abs(diff_kmeans),abs(diff_cutoff)))
      }
      toKeep <- unlist(list(toKeep_kmeans,toKeep_cutoff)[ind_best])
    }
    # add into other sample data tibble
    toKeepByRep <- rbind(toKeepByRep, reads_temp[toKeep,-3])
    toKeepByRep2 <- rbind(toKeepByRep2, reads_temp[toKeep_kmeans,-3])
    
  }
  # For plotting; highlighting
  toKeepByRep$KEEP <- TRUE
  toKeepByRep2$KEEP <- TRUE
  
  ggsave(paste0("01_seq_data_cleaning/byRep_ASV_thresholds_post/",r,"_adaptivethresholds.png"), width = 5, height = 3,
         subMeta %>%
           pivot_longer(cols = one_of(allOTUs), names_to = "ASVID", values_to = "Reads") %>% 
           select(sampleID_qPCR, introRich, Inhibitory, ASVID, Reads) %>%
           left_join(toKeepByRep) %>% mutate(KEEP = ifelse(is.na(KEEP), FALSE, KEEP)) %>%
           # ggplot(aes(y = factor(RichLevel), x = Reads, col = Inhibitory)) +
           ggplot(aes(y = factor(introRich), x = Reads, col = sampleID_qPCR, pch = KEEP)) + facet_grid(Inhibitory~., scales = "free_y", drop = TRUE) +
           geom_jitter(width = 0, height=0.2) + geom_vline(aes(xintercept=0.015)) +geom_vline(aes(xintercept=0.06)) +  scale_x_log10() +
           scale_shape_manual(values =c(21,19))
  )
  ggsave(paste0("01_seq_data_cleaning/byRep_ASV_thresholds_post/",r,"_kmeansthresholds.png"), width = 5, height = 3,
         subMeta %>%
           pivot_longer(cols = one_of(allOTUs), names_to = "ASVID", values_to = "Reads") %>% 
           select(sampleID_qPCR, introRich, Inhibitory, ASVID, Reads) %>%
           left_join(toKeepByRep2) %>% mutate(KEEP = ifelse(is.na(KEEP), FALSE, KEEP)) %>%
           # ggplot(aes(y = factor(RichLevel), x = Reads, col = Inhibitory)) +
           ggplot(aes(y = factor(introRich), x = Reads, col = sampleID_qPCR, pch = KEEP)) + facet_grid(Inhibitory~., scales = "free_y", drop = TRUE) +
           geom_jitter(width = 0, height=0.2) + geom_vline(aes(xintercept=0.015)) + geom_vline(aes(xintercept=0.06))+scale_x_log10() +
           scale_shape_manual(values =c(21,19))
  )
  # Merge with rest of dataset
  toKeepBySampleAll_post <- rbind(toKeepBySampleAll_post, toKeepByRep)
  toKeepBySampleAll_post2 <- rbind(toKeepBySampleAll_post2, toKeepByRep2)
  
}
# reduce duplicates to get list of OTUs
filtOTUs_post <- unique(toKeepBySampleAll_post$ASVID)
# How many reads did we actually lose?
1-sum(t(otu)[,filtOTUs_post])/sum(t(otu)[,allOTUs])
# 1.3%
## Final counts
sum(t(otu)[,allOTUs])-sum(t(otu)[,filtOTUs_post])

1-sum(t(otu)[,unique(c(filtOTUs_post, filtOTUs))])/sum(t(otu)[,allOTUs])
# LOST 0.01064868 reads total

# Number of features lost:
length(unique(allOTUs)) - length(unique(c(filtOTUs_post, filtOTUs)))
# 839

### Filter out OTUs 

# Make new table
asv_table_filt_POST <- toKeepBySampleAll_post %>% select(-KEEP) %>% 
  pivot_wider(names_from = ASVID, values_from = Reads) %>%
  full_join(full_meta %>% filter(Type == "Bd") %>%select(sampleID_qPCR)) %>% select(sampleID_qPCR, everything())
asv_table_filt_POST[is.na(asv_table_filt_POST)] <- 0
full_asv_and_meta_filt1_post <-  full_meta %>% filter(Type =="Bd") %>%
  full_join(asv_table_filt_POST) %>% 
  mutate(obsRich = rowSums(across(one_of(filtOTUs_post))>0)) %>%
  select(sampleID_qPCR, Inhibitory, introRich, obsRich, Rep, one_of(filtOTUs_post), everything()) 


#### Look at data by extraction batch and seq prep batch; are there any biases? ####
# Get list of ASVs present in controls
control_present_asvs <- otu %>% t() %>% as.data.frame() %>% rownames_to_column(var = "sampleID_qPCR") %>%
  filter(sampleID_qPCR %in% c("NTC1","NTC2")) %>%
  pivot_longer(one_of(allOTUs), names_to = "ASVID", values_to = "Reads" ) %>%
  filter(Reads>0) %>% left_join(taxa %>% select(ASVID, ScientificName_unique)) %>%
  select(ASVID, ScientificName_unique, sampleID_qPCR) %>% distinct() %>%
  rename(PCR_batch = sampleID_qPCR)
  
otu_table_by_extraction_btch <- otu %>% t() %>% as.data.frame() %>% rownames_to_column(var = "sampleID_qPCR") %>%
  pivot_longer(one_of(allOTUs), names_to = "ASVID", values_to = "Reads" ) %>%
  # mutate(Reads = ifelse(is.na(Reads), 0, Reads)) %>%
  group_by(sampleID_qPCR) %>% mutate(totalReads = sum(Reads)) %>% ungroup() %>%
  mutate(RelAbund = Reads/totalReads) %>% select(-Reads, -totalReads) %>%
  left_join(extraction_meta %>% rename(sampleID_qPCR = sampleID) %>% select(sampleID_qPCR, extractionID)) %>%
  filter(!is.na(extractionID)) %>%
  separate(extractionID, into = c("extraction_batch","extractionSubID"), remove=FALSE ) %>%
  group_by(extraction_batch, ASVID) %>% mutate(prev = sum(RelAbund>0, na.rm = TRUE)/n()) %>% ungroup() %>%
  left_join(taxa %>% select(ASVID, ScientificName_unique)) %>%
  filter(RelAbund>0)

## Decay curve of ASV prevalence
otu_table_by_extraction_btch %>%
  mutate(`0`=prev>0.0, `0.1` = prev>0.1, `0.2` = prev>0.2, `0.3` = prev>0.3, `0.4` = prev>0.4, `0.6` = prev>0.6, `0.8`=prev>0.8) %>%
  select(extraction_batch, ScientificName_unique, starts_with("0")) %>% distinct() %>%
  pivot_longer(cols = starts_with("0"), names_to = "PrevCutoff", values_to = "TF") %>%
  group_by(extraction_batch, PrevCutoff) %>% summarize(nASV = sum(TF)) %>%
  ungroup() %>% mutate(PrevCutoff = as.numeric(PrevCutoff)) %>%
  ggplot() + geom_line(aes(x=PrevCutoff, y = nASV, col = extraction_batch)) +
  xlab("Prevalence Cutoff") + ylab("Number ASV remaining")

## Plot out high-prevalence ASVs to see if they should be removed
otu_table_by_extraction_btch %>%
  filter( prev>0.7, ASVID %in% c(filtOTUs, filtOTUs_post)) %>%
  left_join(rbind(toKeepBySampleAll,toKeepBySampleAll_post)) %>% 
  filter(!is.na(Reads)) %>%
  # filter(extractionID == asv_table_full_withextrID$extractionID[1]) %>%
  ggplot() + geom_jitter(aes(x=sampleID_qPCR, y=RelAbund, col = ScientificName_unique), show.legend = TRUE, width=0.2, height=0) +
  # geom_line(aes(x=sampleID_qPCR, y=RelAbund, col = ASVID), show.legend = FALSE) +
  facet_wrap(extraction_batch~., drop = TRUE, scale='free') +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

high_prevalence_by_batch <- otu_table_by_extraction_btch %>%
  filter( prev>0.7, ASVID %in% c(filtOTUs, filtOTUs_post)) %>%
  left_join(rbind(toKeepBySampleAll,toKeepBySampleAll_post)) %>% 
  filter(!is.na(Reads)) %>% select(ASVID, ScientificName_unique, extraction_batch) %>% distinct()
all_possible_contam_ASVs <- high_prevalence_by_batch %>% full_join(control_present_asvs)

write.table(all_possible_contam_ASVs, file = "01_seq_data_cleaning/high_prev_ASVs_possible_contam.txt"
            , quote=FALSE, row.names = FALSE, sep="\t")

######## Filtering out ASVs that are possibly contaminants ######
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

temp_meta <- full_join(full_asv_and_meta_filt1, full_asv_and_meta_filt1_post) 
temp_meta[is.na(temp_meta)] <- 0

dir.create("01_seq_data_cleaning/pre_post_ASV_prevalence")
allReps <- temp_meta %>% filter(Rep !="Rep13") %>% pull(Rep) %>% unique()
nExtrBatch <- length(unique(all_possible_contam_ASVs$extraction_batch))
for ( r in allReps ) {
  ggsave(paste0("01_seq_data_cleaning/pre_post_ASV_prevalence/",r,".png"), width = 8, height=6
         ,temp_meta %>% 
           filter(Rep == r, Inhibitory %in% c("Inhib","Non")) %>%
           mutate(Dataset = factor(ifelse(Type=="CV", "Pre","Post"), levels=c("Pre","Post"))) %>%
           select(-one_of(alex_isos$IsolateID), SampleID, sampleID_qPCR, Rep, Inhibitory, RichLevel, introRich, obsRich, Dataset, one_of(c(filtOTUs, filtOTUs_post))) %>%
           pivot_longer(cols = one_of(c(filtOTUs, filtOTUs_post)), names_to = "ASVID", values_to = "Abundance" ) %>%
           group_by(ASVID) %>% mutate(ASVTOTAL = sum(Abundance)) %>% ungroup() %>% filter(ASVTOTAL>0) %>% select(-ASVTOTAL) %>%
           group_by(sampleID_qPCR) %>% mutate(totalAbund = sum(Abundance)) %>% ungroup() %>%
           mutate(RelAbund = Abundance/totalAbund) %>% mutate(RelAbund = ifelse(RelAbund==0, NA, RelAbund)) %>%
           left_join(taxa %>% select(ASVID, ScientificName_unique)) %>% left_join(all_possible_contam_ASVs) %>%
           mutate(high_prev_in_extraction_batch = ifelse(is.na(extraction_batch), "None", extraction_batch)
                  , found_in_PCR_batch = ifelse(is.na(PCR_batch), "None", PCR_batch)) %>%
           # mutate(toRemove = ifelse(sampleID_qPCR %in% toRemove, TRUE, FALSE)) %>%
           ggplot() + geom_point(aes(x=Dataset, y= ScientificName_unique, cex = RelAbund, pch = found_in_PCR_batch, col = high_prev_in_extraction_batch)) +
           facet_grid(Inhibitory ~RichLevel) + ylab("ASV ID") + xlab("Sample Type (Pre/post Bd exposure)") +
           scale_color_manual(values = c(gg_color_hue(nExtrBatch-1), None="black"))
  ) 
}

### manual adjust
# "R5_J15_Pre" # Rich3, Inhib, Remove everything less than 1000 count; just noise here (Bacillus_unknown.40)
# R8 Non Pre (R8_J28_Pre, R8_J7_Pre, R8_J18_Pre) # For all, consider increasing threshold to 17500
# R9 Non Pre (R9_J19_Pre, R9_J28_Pre, R9_J2_Pre) # For all, consider increasing threshold to 10000
# "R1_J22_Bd" # Remove everything less than 0.06
# "R4_J2_Bd" # Remove everything less than 0.03; just noise here.
# "R11_J26_Bd" # Remove everything less than 0.1 reads
# "R12_J19_Bd" # Remove everything less than 0.03;

### Make a new table for pre, with new thresholding requirements
toKeepBySampleAll_manual <- toKeepBySampleAll %>% 
  filter(!(sampleID_qPCR %in% c("R5_J15_Pre") & Reads<1000)) %>%
  filter(!(sampleID_qPCR %in% c("R8_J28_Pre", "R8_J7_Pre", "R8_J18_Pre") & Reads < 17500)) %>%
  filter(!(sampleID_qPCR %in% c("R9_J19_Pre", "R9_J28_Pre", "R9_J2_Pre") & Reads < 10000))
filtOTUs_manual <- unique(toKeepBySampleAll_manual$ASVID)
asv_table_filt_manual <- toKeepBySampleAll_manual %>% select(-KEEP) %>% 
  pivot_wider(names_from = ASVID, values_from = Reads) %>%
  full_join(full_meta %>% filter(Type == "CV") %>%select(sampleID_qPCR)) %>% select(sampleID_qPCR, everything())
asv_table_filt_manual[is.na(asv_table_filt_manual)] <- 0
full_asv_and_meta_filt1_manual <-  full_meta %>% filter(Type =="CV") %>%
  full_join(asv_table_filt_manual) %>% 
  mutate(obsRich = rowSums(across(one_of(filtOTUs_manual))>0)) %>%
  select(sampleID_qPCR, Inhibitory, introRich, obsRich, Rep, one_of(filtOTUs_manual), everything()) 

### Make a new table for pre, with new thresholding requirements
toKeepBySampleAll_post_manual <- toKeepBySampleAll_post %>% 
  filter(!(sampleID_qPCR %in% c("R1_J22_Bd") & Reads<0.06)) %>%
  filter(!(sampleID_qPCR %in% c("R4_J2_Bd") & Reads<0.03)) %>%
  filter(!(sampleID_qPCR %in% c("R11_J26_Bd") & Reads < 0.1)) %>%
  filter(!(sampleID_qPCR %in% c("R12_J19_Bd") & Reads < 0.03))
filtOTUs_post_manual <- unique(toKeepBySampleAll_post_manual$ASVID)
asv_table_filt_POST_manual <- toKeepBySampleAll_post_manual %>% select(-KEEP) %>% 
  pivot_wider(names_from = ASVID, values_from = Reads) %>%
  full_join(full_meta %>% filter(Type == "Bd") %>%select(sampleID_qPCR)) %>% select(sampleID_qPCR, everything())
asv_table_filt_POST_manual[is.na(asv_table_filt_POST_manual)] <- 0
full_asv_and_meta_filt1_post_manual <-  full_meta %>% filter(Type =="Bd") %>%
  full_join(asv_table_filt_POST_manual) %>% 
  mutate(obsRich = rowSums(across(one_of(filtOTUs_post_manual))>0)) %>%
  select(sampleID_qPCR, Inhibitory, introRich, obsRich, Rep, one_of(filtOTUs_post_manual), everything()) 

###### Obs vs Intr Rich ########
# Data
full_asv_and_meta_filt1_manual %>% group_by(sampleID_qPCR) %>%
  ggplot() + geom_point(aes(x=introRich, y = obsRich))
full_asv_and_meta_filt1_post_manual %>% group_by(sampleID_qPCR) %>%
  ggplot() + geom_point(aes(x=introRich, y = obsRich))

## Compare the two
before <- full_asv_and_meta_filt1_manual %>% select(Rep, Inhibitory, RichLevel, obsRich) %>% rename(obsRich_before=obsRich)
after <- full_asv_and_meta_filt1_post_manual %>% select(Rep, Inhibitory, RichLevel, obsRich) %>% rename(obsRich_after=obsRich)
full_join(before,after) %>% ggplot() + geom_jitter(aes(x=obsRich_before, y=obsRich_after, col = Rep), width=0.1, height=0.1) +
  geom_abline(aes(intercept=0,slope=1))

# Remove the samples that are extreme outliers -- flag these
full_asv_and_meta_filt2 <- full_asv_and_meta_filt1_manual %>% mutate(REMOVE_SEQCONTAM = !(obsRich<introRich+5)) %>%
  select(SampleID, sampleID_qPCR, introRich, obsRich,one_of(filtOTUs_manual),REMOVE_SEQCONTAM) 

full_asv_and_meta_filt2_post <- full_asv_and_meta_filt1_post_manual  %>% left_join(before) %>% mutate(REMOVE_SEQCONTAM = !(obsRich<obsRich_before+5)) %>%
  select(SampleID, sampleID_qPCR, introRich, obsRich, one_of(filtOTUs_post_manual),REMOVE_SEQCONTAM) 

final_meta <- full_join(full_asv_and_meta_filt2, full_asv_and_meta_filt2_post) 
final_meta[is.na(final_meta)] <- 0

# Change ASV names so they start with ASV-
colnames(final_meta)[which(colnames(final_meta) %in% c(filtOTUs_manual, filtOTUs_post_manual))] <- paste0("ASV_", colnames(final_meta)[which(colnames(final_meta) %in% c(filtOTUs_manual, filtOTUs_post_manual))])

# right now, the asv_table_filt is "true" reads, whereas the POST is "sequences" reads. 
# Need to make a diff version of the "true" reads one so they match
long_otu <- t(otu) %>% as.data.frame() %>%rownames_to_column(var="sampleID_qPCR") %>%
  pivot_longer(cols=one_of(allOTUs), names_to = "ASVID", values_to = "rawReads")
asv_table_final <- full_join(toKeepBySampleAll_manual, toKeepBySampleAll_post_manual) %>% left_join(long_otu) %>%
  select(-c(Reads,KEEP)) %>% 
  pivot_wider(names_from=ASVID, values_from=rawReads)
asv_table_final[is.na(asv_table_final)] <- 0
colnames(asv_table_final)[which(colnames(asv_table_final) %in% c(filtOTUs_manual, filtOTUs_post_manual))] <- paste0("ASV_",colnames(asv_table_final)[which(colnames(asv_table_final) %in% c(filtOTUs_manual, filtOTUs_post_manual))])
ASV_metadata <- inhibMeta %>% filter(ASVID %in% c(filtOTUs_manual, filtOTUs_post_manual)) %>% select(-c("knownInhibitory","Inhibitory01", "dataset", "IsolateID","Isolate","tipColInhib","tipColDataset")) %>%
  rowwise() %>% mutate(inhibitory_final = weighted.mean(c(inhibitory_probability_alex, inhibitory_probability_woodhams, inhibitory_probability_tree), c(0.6, 0.3, 0.1), na.rm = TRUE)) %>%
  left_join(taxa) %>% select(ASVID, ScientificName_unique, inhibitory_final, inhibitory_probability_alex,inhibitory_probability_woodhams, everything()) %>%
  mutate(temp = "ASV", ASVID_original=ASVID) %>% unite(temp, ASVID, col="ASVID", sep="_")


# Get filtered metadata on inhibitory and taxonomy
dir.create("01_seq_data_cleaning/downstream")

write.table(final_meta, file = "01_seq_data_cleaning/downstream/metadata_with_sequence.txt", quote=FALSE, row.names = FALSE, sep = "\t")
write.table(asv_table_final, file = "01_seq_data_cleaning/downstream/asv_table_filtered.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(ASV_metadata, file = "01_seq_data_cleaning/downstream/asv_metadata.txt", quote=FALSE, row.names = FALSE, sep = "\t")

# Make version of raw table for future reference
pre_otu_real <- otu_real_abund_round %>% as.data.frame() %>%
  rownames_to_column(var="sampleID_qPCR")
write.table(pre_otu_real, file = "01_seq_data_cleaning/downstream/pre_otu_real_abundance.txt", quote=FALSE, row.names=FALSE, sep="\t")


setwd("..")
