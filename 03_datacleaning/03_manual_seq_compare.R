#!bin/bash

###### Filtering out contaminants ###########
library(tidyverse)
library(ape)
library(phytools)

setwd("03_datacleaning/")
dir.create("03_manual_seq_compare")
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

#### Downstream from sequence filtering
# asv_filtered <- read.delim("01_seq_data_cleaning/downstream/asv_table_filtered.txt")
meta_asv_filtered <- read.delim("02_merge_all_data/metadata_allReps.txt")
ASV_fromFiltered <- colnames(meta_asv_filtered)[grep("ASV_",colnames(meta_asv_filtered))]
ASV_absent_pre <- names(which(colSums(meta_asv_filtered[,ASV_fromFiltered])==0))
ASV_fromFiltered_present <- ASV_fromFiltered[!(ASV_fromFiltered %in% ASV_absent_pre)]

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
allASVs_ASV <- paste0("ASV_",allOTUs)
# Rename OTU samplenames
colnames(otu) <- gsub(".","_", gsub("..","-", colnames(otu), fixed = TRUE), fixed=TRUE)

##### Combine ASV with meta #####
# Bind with qPCR and find out "real" abundances
# Make relative abundance, PA
otu_relabund <- apply(otu, MARGIN = 2, FUN = function(x) x/sum(x))

#Transpose so we can combine with qPCR
otu_with_qPCR <- t(otu_relabund) %>% as.data.frame() %>% 
  rename_with(.fn=~paste0("ASV_",.x)) %>%
  rownames_to_column(var = "sampleID_qPCR") %>% as_tibble() %>%left_join(full_meta) %>%
  filter(!is.na(qPCR_bact)) %>% select(sampleID_qPCR, one_of(allASVs_ASV), qPCR_bact) %>% as.data.frame()
rownames(otu_with_qPCR) <- otu_with_qPCR$sampleID_qPCR
otu_real_abund <- apply(otu_with_qPCR[,allASVs_ASV], MARGIN = 2, FUN = function(x) x*otu_with_qPCR$qPCR_bact)
# Round
otu_real_abund_round <- round(otu_real_abund)

# Combine with basic metadata
full_asv_and_meta <- otu_real_abund_round %>% as.data.frame() %>% rownames_to_column(var = "sampleID_qPCR") %>%
  mutate(sums = rowSums(across(where(is.numeric)))) %>% filter(!is.na(sums), sums>0) %>%
  select(-sums) %>% left_join(full_meta) %>%   mutate(Control = ifelse(Inhibitory=="Con", TRUE, FALSE)) 

full_asvfilt_meta <- full_meta %>% select(SampleID, sampleID_qPCR, RichLevel, introRich, Rep, Type, Inhibitory, one_of(alex_isos$IsolateID)) %>%
  right_join(meta_asv_filtered)

addASVBack <- allASVs_ASV[!(allASVs_ASV %in% colnames(full_asvfilt_meta))]
blankASVs <- as.data.frame(matrix(0, ncol=length(addASVBack), nrow=nrow(full_asvfilt_meta), dimnames = list(NULL, addASVBack)))
full_asvfilt_meta <- cbind(full_asvfilt_meta, blankASVs)
#### Look at individual ASV distributions

dat_CVonly <- full_asv_and_meta %>% filter(Type == "CV")
dat_CVonly_FILT <- full_asvfilt_meta %>% filter(Type=="CV")

#### By ISOLATE
allIsolatesKept <- alex_isos$IsolateID[alex_isos$IsolateID %in% colnames(dat_CVonly)]

isolateQuery_potentialMatches <- data.frame()
isolateQuery_potentialMatches_FILT <- data.frame()

for (iso in allIsolatesKept) {
  # for ALL ASVs
  datCurrentIso <- dat_CVonly %>%
    filter(get(iso))%>% 
    select(sampleID_qPCR, introRich, starts_with("ASV")) 
  ASVToKeep <- datCurrentIso %>% select(starts_with("ASV")) %>% colSums() %>% as.data.frame() %>%
    filter(.>0) %>% rename(Reads=".") %>% rownames()
  summary_potentialMatches <- datCurrentIso %>% 
    pivot_longer(-c(sampleID_qPCR, introRich), names_to="ASVID", values_to="Reads") %>%
    mutate(PA = Reads>0, approxCells = 1/introRich) %>%
    mutate(IsolateID=paste0(iso)) %>% group_by(ASVID) %>% 
    summarise(IsolateID=paste0(iso), percentMatch_isolateQuery = sum(PA)/length(PA), propMatch_isolateQuery = sum(PA*approxCells)/sum(approxCells)) %>%
    select(IsolateID, ASVID, percentMatch_isolateQuery,propMatch_isolateQuery)
  isolateQuery_potentialMatches <- rbind(isolateQuery_potentialMatches, summary_potentialMatches)
  
  # for ALL ASVs
  datCurrentIso_FILT <- dat_CVonly_FILT %>%
    filter(get(iso))%>% 
    select(sampleID_qPCR, introRich, starts_with("ASV")) 
  ASVToKeep_FILT <- datCurrentIso_FILT %>% select(starts_with("ASV")) %>% colSums() %>% as.data.frame() %>%
    filter(.>0) %>% rename(Reads=".") %>% rownames()
  summary_potentialMatches_FILT <- datCurrentIso_FILT %>% 
    pivot_longer(-c(sampleID_qPCR, introRich), names_to="ASVID", values_to="Reads") %>%
    mutate(PA = Reads>0, approxCells = 1/introRich) %>%
    mutate(IsolateID=paste0(iso)) %>% group_by(ASVID) %>% 
    summarise(IsolateID=paste0(iso), percentMatch_isolateQuery = sum(PA)/length(PA), propMatch_isolateQuery = sum(PA*approxCells)/sum(approxCells)) %>%
    select(IsolateID, ASVID, percentMatch_isolateQuery,propMatch_isolateQuery)
  isolateQuery_potentialMatches_FILT <- rbind(isolateQuery_potentialMatches_FILT, summary_potentialMatches_FILT)
}

# Look at distributions of matches
isolateQuery_potentialMatches %>%
  ggplot() + geom_histogram(aes(x=propMatch_isolateQuery), bins=15) +
  facet_wrap(.~IsolateID)
isolateQuery_potentialMatches_FILT %>%
  ggplot() + geom_histogram(aes(x=propMatch_isolateQuery), bins=15) +
  facet_wrap(.~IsolateID)

## Compare two metrics
isolateQuery_potentialMatches %>%
  ggplot() + geom_point(aes(x=percentMatch_isolateQuery, y=propMatch_isolateQuery)) +
  facet_wrap(.~IsolateID)
isolateQuery_potentialMatches_FILT %>%
  ggplot() + geom_point(aes(x=percentMatch_isolateQuery, y=propMatch_isolateQuery)) +
  facet_wrap(.~IsolateID)
# # Filter out matches that are "Maximum" match
# isolateQuery_potentialMatches_FILT <- isolateQuery_potentialMatches %>% group_by(IsolateQuery) %>%
#   mutate(maxMatch=max(percentMatch)) %>% filter(percentMatch==maxMatch) %>% ungroup()

### Now, do reverse match search by ASVID
allASVsKept <- names(which(dat_CVonly %>% select(starts_with("ASV")) %>% colSums() >0))
asvQuery_potentialMatches <- data.frame()
for ( a in allASVsKept) {
  datCurrentASV <- dat_CVonly %>% filter(Inhibitory !="Con") %>%
    filter(get(a)>0)%>% 
    select(sampleID_qPCR, one_of(a), one_of(allIsolatesKept)) 
  summary_potentialMatches <- datCurrentASV %>% 
    pivot_longer(-c(sampleID_qPCR, one_of(a)), names_to="IsolateID", values_to="PA") %>%
    mutate(PA = ifelse(is.na(PA), FALSE, PA)) %>%
    group_by(IsolateID) %>%
    summarize(propMatch_ASVQuery = sum(PA*get(a))/sum(get(a)), percentMatch_ASVQuery = sum(PA)/n())  %>%
    mutate(ASVID=a)
  
  asvQuery_potentialMatches <- rbind(asvQuery_potentialMatches, summary_potentialMatches)  
}

#Now with filtered
allASVsKept_FILT <- names(which(dat_CVonly_FILT %>% select(starts_with("ASV")) %>% colSums() >0))
asvQuery_potentialMatches_FILT <- data.frame()
for ( a in allASVsKept_FILT) {
  datCurrentASV_FILT <- dat_CVonly_FILT %>% filter(Inhibitory !="Con") %>%
    filter(get(a)>0)%>% 
    select(sampleID_qPCR, one_of(a), one_of(allIsolatesKept)) 
  summary_potentialMatches_FILT <- datCurrentASV_FILT %>% 
    pivot_longer(-c(sampleID_qPCR, one_of(a)), names_to="IsolateID", values_to="PA") %>%
    mutate(PA = ifelse(is.na(PA), FALSE, PA)) %>%
    group_by(IsolateID) %>%
    summarize(propMatch_ASVQuery = sum(PA*get(a))/sum(get(a)), percentMatch_ASVQuery = sum(PA)/n())  %>%
    mutate(ASVID=a)
  
  asvQuery_potentialMatches_FILT <- rbind(asvQuery_potentialMatches_FILT, summary_potentialMatches_FILT)  
}


# Look at distributions of matches
asvQuery_potentialMatches %>%
  ggplot() + geom_histogram(aes(x=propMatch_ASVQuery), bins=15) +
  facet_wrap(.~IsolateID)
asvQuery_potentialMatches_FILT %>%
  ggplot() + geom_histogram(aes(x=propMatch_ASVQuery), bins=15) +
  facet_wrap(.~IsolateID)

## Compare two metrics
asvQuery_potentialMatches %>%
  ggplot() + geom_point(aes(x=percentMatch_ASVQuery, y=propMatch_ASVQuery)) +
  facet_wrap(.~IsolateID)
asvQuery_potentialMatches_FILT %>%
  ggplot() + geom_point(aes(x=percentMatch_ASVQuery, y=propMatch_ASVQuery)) +
  facet_wrap(.~IsolateID)

### Plot by individual isolate

twoWayQuery_matches <- isolateQuery_potentialMatches %>% full_join(asvQuery_potentialMatches)
twoWayQuery_matches_FILT <- isolateQuery_potentialMatches_FILT %>% full_join(asvQuery_potentialMatches_FILT)

### Filter ASVs by whether they match >50%
twoWayQuery_matches_filt <- twoWayQuery_matches %>% 
  group_by(IsolateID) %>% mutate(isoRank = rank(-propMatch_isolateQuery, ties.method = "min")) %>% 
  mutate(asvRank = rank(-propMatch_ASVQuery, ties.method = "min")) %>% ungroup() %>% 
  # filter(propMatch_isolateQuery>=isoThresh, propMatch_ASVQuery>=asvThresh)
  filter(((asvRank<5 & propMatch_isolateQuery>0.25 )|( propMatch_isolateQuery>0.5 & propMatch_ASVQuery>0.5)))

twoWayQuery_matches_FILT_filt <- twoWayQuery_matches_FILT %>% 
  group_by(IsolateID) %>% mutate(isoRank = rank(-propMatch_isolateQuery, ties.method = "min")) %>% 
  mutate(asvRank = rank(-propMatch_ASVQuery, ties.method = "min")) %>% ungroup() %>% 
  # filter(propMatch_isolateQuery>=isoThresh, propMatch_ASVQuery>=asvThresh)
  filter(( propMatch_isolateQuery>0.25 & propMatch_ASVQuery>0.25))
# Give unique ASV references so the plots look nicer
# ASVREF_FILT <- data.frame(ASVID = unique(twoWayQuery_matches_FILT_filt$ASVID), ASVRef = paste0("REF_",seq(1,length(unique(twoWayQuery_matches_FILT_filt$ASVID)))))
# Give unique ASV references so the plots look nicer
ASVREF <- data.frame(ASVID = unique(c(twoWayQuery_matches_filt$ASVID, twoWayQuery_matches_FILT_filt$ASVID)), ASVRef = paste0("REF_",seq(1,length(unique(c(twoWayQuery_matches_filt$ASVID, twoWayQuery_matches_FILT_filt$ASVID))))))
twoWayQuery_matches_filt <- twoWayQuery_matches_filt %>% left_join(ASVREF)
twoWayQuery_matches_FILT_filt <- twoWayQuery_matches_FILT_filt %>% left_join(ASVREF)
# 
# 
# currN <- which(allIsolatesKept=="2E")
# # currN <- 1
# temp_asv <- twoWayQuery_matches_filt %>% 
#   filter(IsolateID==allIsolatesKept[currN]) %>%
#   pull(ASVID)
# temp_asv_FILT <- twoWayQuery_matches_FILT_filt %>% 
#   filter(IsolateID==allIsolatesKept[currN]) %>%
#   pull(ASVID)
# 
# twoWayQuery_matches %>%
#   filter(IsolateID==allIsolatesKept[currN]) %>% 
#   rowwise() %>% mutate(Kept=ifelse(ASVID %in% temp_asv, TRUE, FALSE))  %>% ungroup() %>%
#   ggplot() + geom_point(aes(x=propMatch_isolateQuery, y=propMatch_ASVQuery, col=Kept)) +
#   labs(main=allIsolatesKept[currN])
# 
# twoWayQuery_matches_FILT %>%
#   filter(IsolateID==allIsolatesKept[currN]) %>% 
#   rowwise() %>% mutate(Kept=ifelse(ASVID %in% temp_asv, TRUE, FALSE))  %>% ungroup() %>%
#   ggplot() + geom_point(aes(x=propMatch_isolateQuery, y=propMatch_ASVQuery, col=Kept)) +
#   labs(main=allIsolatesKept[currN])

### Plot all data by ASV and isolate to see if it matches up
# twoWayQuery_matches_filt %>% group_by(IsolateID) %>%
#   summarise(nASVs = length(unique(ASVID))) %>% arrange(-nASVs)
# twoWayQuery_matches_FILT_filt %>% group_by(IsolateID) %>%
#   summarise(nASVs = length(unique(ASVID))) %>% arrange(-nASVs)
# Trying doing multiple ASVs per isolate
dir.create("03_manual_seq_compare/byIsolate_asvoccurance")
for ( tempIsolate in allIsolatesKept) {
  # tempIsolate <- allIsolatesKept[32]
  tempCompare <- twoWayQuery_matches_filt %>% filter(IsolateID==tempIsolate)
  tempASV <- tempCompare$ASVID
  # FOr filt
  tempCompare_FILT <- twoWayQuery_matches_FILT_filt %>% filter(IsolateID==tempIsolate)
  tempASV_FILT <- tempCompare_FILT$ASVID
  if (length(c(tempASV,tempASV_FILT))==0) {
    dat_temp <- full_asv_and_meta %>% filter(Inhibitory!="Con", Type=="CV") %>%
      select(SampleID, sampleID_qPCR, RichLevel, Rep, Inhibitory, one_of(c(tempASV,tempASV_FILT)), one_of(tempIsolate)) %>%
      mutate(Rep = factor(Rep, levels=c(paste0("Rep",seq(1,12))))) %>%
      rename_at(vars(all_of(tempIsolate)), ~"IsolateIntroduced") %>%
      mutate(IsolateIntroduced=ifelse(is.na(IsolateIntroduced), FALSE, TRUE)) %>%
      mutate(ASVID="NoASVMatches", Reads = NA)
    tempASV <- "NoASVMatches"
  } else {
    dat_temp <- full_asv_and_meta %>% filter(Inhibitory!="Con", Type=="CV") %>%
      select(SampleID, sampleID_qPCR, RichLevel, Rep, Inhibitory, one_of(c(tempASV,tempASV_FILT)), one_of(tempIsolate)) %>%
      mutate(Rep = factor(Rep, levels=c(paste0("Rep",seq(1,12))))) %>%
      rename_at(vars(all_of(tempIsolate)), ~"IsolateIntroduced") %>%
      mutate(IsolateIntroduced=ifelse(is.na(IsolateIntroduced), FALSE, TRUE)) %>%
      pivot_longer(one_of(c(tempASV,tempASV_FILT)), names_to = "ASVID", values_to="Reads")
  }
  if (length(tempASV_FILT)==0 && tempASV == "NoASVMatches") {
    dat_temp_FILT <- full_asvfilt_meta %>% filter(Inhibitory!="Con", Type=="CV") %>%
      select(SampleID, sampleID_qPCR, RichLevel, Rep, Inhibitory, one_of(c(tempASV,tempASV_FILT)), one_of(tempIsolate)) %>%
      mutate(Rep = factor(Rep, levels=c(paste0("Rep",seq(1,12))))) %>%
      rename_at(vars(all_of(tempIsolate)), ~"IsolateIntroduced") %>%
      mutate(IsolateIntroduced=ifelse(is.na(IsolateIntroduced), FALSE, TRUE)) %>%
      mutate(ASVID="NoASVMatches", Reads = NA)
    tempASV_FILT <- "NoASVMatches"
  } else {
    dat_temp_FILT <- full_asvfilt_meta %>% filter(Inhibitory!="Con", Type=="CV") %>%
      select(SampleID, sampleID_qPCR, RichLevel, Rep, Inhibitory, one_of(c(tempASV,tempASV_FILT)), one_of(tempIsolate)) %>%
      mutate(Rep = factor(Rep, levels=c(paste0("Rep",seq(1,12))))) %>%
      rename_at(vars(all_of(tempIsolate)), ~"IsolateIntroduced") %>%
      mutate(IsolateIntroduced=ifelse(is.na(IsolateIntroduced), FALSE, TRUE)) %>%
      pivot_longer(one_of(c(tempASV,tempASV_FILT)), names_to = "ASVID", values_to="Reads")
  }
  ### Plotting
  gg_isotemp <- dat_temp %>%
    left_join(ASVREF) %>%
    # rename_at(vars(all_of(tempASV)), ~"Reads") %>%
    mutate(Reads = ifelse(Reads==0, NA, Reads)) %>%
    ggplot(aes(x=factor(RichLevel), y=Rep)) + geom_point(aes(pch=IsolateIntroduced, fill=log(Reads+1), cex=IsolateIntroduced), col="black") +
    scale_shape_manual(values=c(21,23)) +
    scale_size_manual(values=c(2,5)) +
    scale_fill_gradient(low="white", high="red", na.value = "grey") +
    facet_grid(Inhibitory ~ ASVRef) +xlab("Richness (hypothetical)") 
  gg_isotemp
  gg_isotemp_FILT <- dat_temp_FILT %>%
    left_join(ASVREF) %>%
    # rename_at(vars(all_of(tempASV)), ~"Reads") %>%
    mutate(Reads = ifelse(Reads==0, NA, Reads)) %>%
    ggplot(aes(x=factor(RichLevel), y=Rep)) + geom_point(aes(pch=IsolateIntroduced, fill=log(Reads+1), cex=IsolateIntroduced), col="black") +
    scale_shape_manual(values=c(21,23)) +
    scale_size_manual(values=c(2,5)) +
    scale_fill_gradient(low="white", high="red", na.value = "grey") +
    facet_grid(Inhibitory ~ ASVRef) +xlab("Richness (hypothetical)") 
  gg_isotemp_FILT
  ggsave(filename = paste0("03_manual_seq_compare/byIsolate_asvoccurance/", tempIsolate, "_gr50ASV.png"), height=7, width=2+1.5*length(c(tempASV,tempASV_FILT)), limitsize=FALSE
         , gg_isotemp)
  ggsave(filename = paste0("03_manual_seq_compare/byIsolate_asvoccurance/", tempIsolate, "_gr50ASV_FILT.png"), height=7, width=2+1.5*length(c(tempASV,tempASV_FILT)), limitsize=FALSE
         , gg_isotemp_FILT)
}


#### Save a spreadsheet to note which ones are possible and which ones to eliminate
taxa2 <- taxa %>% unite(Class, Order, Family, Genus, Species, col="taxastring")
# alex_isos %>% left_join(taxa2) %>% select(IsolateID, taxastring)
alex_isos_IDs <- alex_isos %>% left_join(taxa2) %>% select(IsolateID, taxastring) %>% rename(Taxon_Isolate = taxastring) %>%
  arrange(IsolateID) %>% mutate(dup = duplicated(IsolateID)) %>% filter(!dup) %>% select(-dup)
blank_spreadsheet <- twoWayQuery_matches_filt %>% full_join(twoWayQuery_matches_FILT_filt) %>% left_join(taxa2 %>% mutate(ASVID=paste0("ASV_",ASVID))) %>%
  mutate(Possible = "") %>% select(IsolateID, ASVID, ASVRef, Possible, taxastring) %>% distinct() %>% left_join(alex_isos_IDs) %>%
  arrange(IsolateID, ASVRef)
write.table(blank_spreadsheet, file = "03_manual_seq_compare/blank_spreadsheet_for_possible_matches2.txt", sep="\t", quote=FALSE, row.names = FALSE)


##### Re-upload MANUAL assignments ######
seq_compare <- read.csv("03_manual_seq_compare/blank_spreadsheet_for_possible_matches_CHECKED.csv", sep="\t")
# table(seq_compare$Possible)
### First, get all matches that we are confident in
seq_compare %>% filter(Possible=="MATCH") %>% mutate(SameTaxa=taxastring==Taxon_Isolate) %>% pull(SameTaxa)
exactMatches <- seq_compare %>% filter(Possible=="MATCH") %>% select(IsolateID, ASVRef, ASVID)

### Then, get all non-confident matches, but filter out confident ones
seq_compare_remaining <- seq_compare %>% filter(!(ASVRef %in% exactMatches$ASVRef)) %>% 
  filter(!(IsolateID %in% exactMatches$IsolateID)) %>%
  filter(Possible%in%c("Potential")) 

### Confirm potentials that only match one ASV; and ASVs that only match one Isolate; and are also in filtered set
seq_compare_remaining %>%
  mutate(SameTaxa = taxastring==Taxon_Isolate) %>%
  mutate(IsoDup = duplicated(IsolateID), ASVdup = duplicated(ASVRef)) %>%
  group_by(IsolateID) %>% mutate(OnlyASVHit_byiso = !any(IsoDup), OnlyIsolateHit_byASV = !any(ASVdup)) %>% 
  filter(SameTaxa,OnlyASVHit_byiso,OnlyIsolateHit_byASV)
# 2F Confirmed as Ref_71
# 35E confirmed as REF_49 (maybe merged with REF_50, REF_51)
# 45E confirmed as REF_28
exactMatches <- seq_compare %>% filter(IsolateID=="2F"&ASVRef=="REF_71" |
                                         IsolateID=="35E"&ASVRef%in%c("REF_49") |
                                         IsolateID=="45E"&ASVRef=="REF_28"  ) %>% 
  select(IsolateID, ASVRef, ASVID) %>%
  full_join(exactMatches)
# Update remaining
seq_compare_remaining2  <- seq_compare_remaining %>%
  filter(!(IsolateID %in% exactMatches$IsolateID)) %>%
  filter(!(ASVRef %in% exactMatches$ASVRef))


### Relax requirement for SameTaxa being true
seq_compare_remaining2 %>%
  mutate(SameTaxa = taxastring==Taxon_Isolate) %>%
  mutate(IsoDup = duplicated(IsolateID), ASVdup = duplicated(ASVRef)) %>%
  group_by(IsolateID) %>% mutate(OnlyASVHit_byiso = !any(IsoDup), OnlyIsolateHit_byASV = !any(ASVdup)) %>% 
  filter(OnlyASVHit_byiso,OnlyIsolateHit_byASV) %>%
  select(IsolateID, ASVRef)
# 13C confirmed REF_51
# 24L could be REF_100, 102, 103 (103 is closest taxonomic match, but not abundant)- NOT 101
# 25A confirmed REF_98; could be combined with REF_97 and REF_99
# 36G confirmed REF_43; could be merged with REF_41,42, 44 (44 overlap with 36A, but 36A better fit elsewhere)
# 40A confirmed REF_3
# 50B confirmed REF_21; could be merged with REF19, 20, 22
exactMatches <- seq_compare %>% 
  filter(IsolateID=="13C"&ASVRef%in%c("REF_51") |
          IsolateID=="24L"& ASVRef%in% c("REF_100", "REF_102", "REF_103") |
           IsolateID=="25A"&ASVRef%in%c("REF_98","REF_97","REF_99") |
           IsolateID=="36G"&ASVRef%in% c("REF_41","REF_42")  |
           IsolateID=="40A"&ASVRef%in% c("REF_3")  |
           IsolateID=="50B"&ASVRef%in% c("REF_21","REF_19","REF_20","REF_22")  
  ) %>% 
  select(IsolateID, ASVRef, ASVID) %>%
  full_join(exactMatches)

# Update remaining
seq_compare_remaining3  <- seq_compare_remaining2 %>%
  filter(!(IsolateID %in% exactMatches$IsolateID)) %>%
  filter(!(ASVRef %in% exactMatches$ASVRef))


### Look at Isolates with ONLY one potential hit
seq_compare_remaining3 %>%
  mutate(SameTaxa = taxastring==Taxon_Isolate) %>%
  mutate(IsoDup = duplicated(IsolateID), ASVdup = duplicated(ASVRef)) %>%
  group_by(IsolateID) %>% mutate(OnlyASVHit_byiso = !any(IsoDup), OnlyIsolateHit_byASV = !any(ASVdup)) %>% 
  filter(OnlyASVHit_byiso) %>%
  select(IsolateID, ASVRef)
# 33D is likely REF_59 (maybe merged with REF_58)
# 40B COULD be REF34, but this matches 35B better
# 50A could be REF_25, but actually could also be REF_11-- and REF_11 doesn't match anything else.
# 50A is REF_11
exactMatches <- seq_compare %>% 
  filter(IsolateID=="33D"& ASVRef%in% c("REF_59", "REF_58") |
           IsolateID=="50A"&ASVRef%in%c("REF_11") 
  ) %>% 
  select(IsolateID, ASVRef, ASVID) %>%
  full_join(exactMatches)

# Update remaining
seq_compare_remaining4  <- seq_compare_remaining3 %>%
  filter(!(IsolateID %in% exactMatches$IsolateID)) %>%
  filter(!(ASVRef %in% exactMatches$ASVRef))


### Look at ASVs with ONLY one potential hit
seq_compare_remaining4 %>%
  mutate(SameTaxa = taxastring==Taxon_Isolate) %>%
  mutate(IsoDup = duplicated(IsolateID), ASVdup = duplicated(ASVRef)) %>%
  group_by(IsolateID) %>% mutate(OnlyASVHit_byiso = !any(IsoDup), OnlyIsolateHit_byASV = !any(ASVdup)) %>% 
  filter(OnlyIsolateHit_byASV) %>%
  select(IsolateID, ASVRef, taxastring, Taxon_Isolate)
# 35B is REF_34
# 36A is 44; could be 48 too but not likely
# 36E is REF_45 and REF_46, probably merged
# 37H is REF_39 and REF_40, probably merged
exactMatches <- seq_compare %>% 
  filter(IsolateID=="35B"& ASVRef%in% c("REF_34") |
           IsolateID=="36A"&ASVRef%in%c("REF_44")  |
           IsolateID=="36E"&ASVRef%in%c("REF_45","REF_46")  |
           IsolateID=="37H"&ASVRef%in%c("REF_39","REF_40")  
           
  ) %>% 
  select(IsolateID, ASVRef, ASVID) %>%
  full_join(exactMatches)

# Update remaining
seq_compare_remaining5  <- seq_compare_remaining4 %>%
  filter(!(IsolateID %in% exactMatches$IsolateID)) %>%
  filter(!(ASVRef %in% exactMatches$ASVRef))

#### Manually assign remaining 4 Isolates
seq_compare_remaining5 %>%   select(IsolateID, ASVRef, taxastring, Taxon_Isolate) %>% arrange(ASVRef)
# REF_25 belongs to 2H
# REF_14 belongs to 34B
# REF_81 belongs to 2E
# REF_38 belongs to 32G, but might also be a contaminant.
exactMatches <- seq_compare %>% 
  filter(IsolateID=="2H"& ASVRef%in% c("REF_25") |
           IsolateID=="34B"&ASVRef%in%c("REF_14")  |
           IsolateID=="2E"&ASVRef%in%c("REF_81")  |
           IsolateID=="32G"&ASVRef%in%c("REF_38")  
         
  ) %>% 
  select(IsolateID, ASVRef, ASVID) %>%
  full_join(exactMatches)

### Check out merged ones
seq_compare %>% filter(!(ASVID %in% exactMatches$ASVID)) %>%
  filter(!(IsolateID %in% exactMatches$IsolateID)) %>%
  filter(Possible=="Merge")
# 34F is REF 54, 55, 56, 57, 69; maybe all merged? 
exactMatches <- seq_compare %>% 
  filter(IsolateID=="34F"& ASVRef%in% c("REF_54", "REF_55", "REF_56", "REF_57", "REF_69") 
  ) %>% 
  select(IsolateID, ASVRef, ASVID) %>%
  full_join(exactMatches)



### Check which ASVs have not been assigned an isolate
importantASVs <- ASV_fromFiltered_present
unmatched_ASVs <- importantASVs[!(importantASVs %in% exactMatches$ASVID)]
unmatched_isolates <- allIsolatesKept[!allIsolatesKept %in% exactMatches$IsolateID]

#### UNMATCHED ############

#### Get plots of unmatched ASVs
# meta_asv_filtered$sampleID_qPCR
# meta_asv_filtered %>% select(sampleID_qPCR, RichLevel, Inhibitory, one_of(unmatched_ASVs)) %>% View()
########### All these remaining ASVs aren't actually found in any of the "Pre" ones!
# I will actually manually remove this because I don't think it actually means anything? Let's check.
taxa %>% mutate(ASVID = paste0("ASV_", ASVID)) %>% filter(ASVID==unmatched_ASVs)


dir.create("03_manual_seq_compare/byIsolate_asvoccurance_unmatched")
for ( tempIsolate in unmatched_isolates) {
  # tempIsolate <- allIsolatesKept[32]
  tempCompare <- twoWayQuery_matches %>% filter(IsolateID==tempIsolate, ASVID %in% unmatched_ASVs)
  tempASV <- tempCompare$ASVID
  # FOr filt
  tempCompare_FILT <- twoWayQuery_matches_FILT %>% filter(IsolateID==tempIsolate, ASVID %in% unmatched_ASVs)
  tempASV_FILT <- tempCompare_FILT$ASVID
  if (length(c(tempASV,tempASV_FILT))==0) {
    dat_temp <- full_asv_and_meta %>% filter(Inhibitory!="Con", Type=="CV") %>%
      select(SampleID, sampleID_qPCR, RichLevel, Rep, Inhibitory, one_of(c(tempASV,tempASV_FILT)), one_of(tempIsolate)) %>%
      mutate(Rep = factor(Rep, levels=c(paste0("Rep",seq(1,12))))) %>%
      rename_at(vars(all_of(tempIsolate)), ~"IsolateIntroduced") %>%
      mutate(IsolateIntroduced=ifelse(is.na(IsolateIntroduced), FALSE, TRUE)) %>%
      mutate(ASVID="NoASVMatches", Reads = NA)
    tempASV <- "NoASVMatches"
  } else {
    dat_temp <- full_asv_and_meta %>% filter(Inhibitory!="Con", Type=="CV") %>%
      select(SampleID, sampleID_qPCR, RichLevel, Rep, Inhibitory, one_of(c(tempASV,tempASV_FILT)), one_of(tempIsolate)) %>%
      mutate(Rep = factor(Rep, levels=c(paste0("Rep",seq(1,12))))) %>%
      rename_at(vars(all_of(tempIsolate)), ~"IsolateIntroduced") %>%
      mutate(IsolateIntroduced=ifelse(is.na(IsolateIntroduced), FALSE, TRUE)) %>%
      pivot_longer(one_of(c(tempASV,tempASV_FILT)), names_to = "ASVID", values_to="Reads")
  }
  if (length(tempASV_FILT)==0 && tempASV == "NoASVMatches") {
    dat_temp_FILT <- full_asvfilt_meta %>% filter(Inhibitory!="Con", Type=="CV") %>%
      select(SampleID, sampleID_qPCR, RichLevel, Rep, Inhibitory, one_of(c(tempASV,tempASV_FILT)), one_of(tempIsolate)) %>%
      mutate(Rep = factor(Rep, levels=c(paste0("Rep",seq(1,12))))) %>%
      rename_at(vars(all_of(tempIsolate)), ~"IsolateIntroduced") %>%
      mutate(IsolateIntroduced=ifelse(is.na(IsolateIntroduced), FALSE, TRUE)) %>%
      mutate(ASVID="NoASVMatches", Reads = NA)
    tempASV_FILT <- "NoASVMatches"
  } else {
    dat_temp_FILT <- full_asvfilt_meta %>% filter(Inhibitory!="Con", Type=="CV") %>%
      select(SampleID, sampleID_qPCR, RichLevel, Rep, Inhibitory, one_of(c(tempASV,tempASV_FILT)), one_of(tempIsolate)) %>%
      mutate(Rep = factor(Rep, levels=c(paste0("Rep",seq(1,12))))) %>%
      rename_at(vars(all_of(tempIsolate)), ~"IsolateIntroduced") %>%
      mutate(IsolateIntroduced=ifelse(is.na(IsolateIntroduced), FALSE, TRUE)) %>%
      pivot_longer(one_of(c(tempASV,tempASV_FILT)), names_to = "ASVID", values_to="Reads")
  }
  
  ### Plotting
  gg_isotemp <- dat_temp %>%
    left_join(ASVREF) %>%
    rowwise() %>%
    mutate(ASVRef = ifelse(is.na(ASVRef), substr(ASVID, nchar(ASVID)-8, nchar(ASVID)),ASVRef)) %>% ungroup() %>%
    # rename_at(vars(all_of(tempASV)), ~"Reads") %>%
    mutate(Reads = ifelse(Reads==0, NA, Reads)) %>%
    ggplot(aes(x=factor(RichLevel), y=Rep)) + geom_point(aes(pch=IsolateIntroduced, fill=log(Reads+1), cex=IsolateIntroduced), col="black") +
    scale_shape_manual(values=c(21,23)) +
    scale_size_manual(values=c(2,5)) +
    scale_fill_gradient(low="white", high="red", na.value = "grey") +
    facet_grid(Inhibitory ~ ASVRef) +xlab("Richness (hypothetical)") 
  # gg_isotemp
  gg_isotemp_FILT <- dat_temp_FILT %>%
    left_join(ASVREF) %>%
    # rename_at(vars(all_of(tempASV)), ~"Reads") %>%
    mutate(Reads = ifelse(Reads==0, NA, Reads)) %>%
    ggplot(aes(x=factor(RichLevel), y=Rep)) + geom_point(aes(pch=IsolateIntroduced, fill=log(Reads+1), cex=IsolateIntroduced), col="black") +
    scale_shape_manual(values=c(21,23)) +
    scale_size_manual(values=c(2,5)) +
    scale_fill_gradient(low="white", high="red", na.value = "grey") +
    facet_grid(Inhibitory ~ ASVRef) +xlab("Richness (hypothetical)") 
  gg_isotemp_FILT
  ggsave(filename = paste0("03_manual_seq_compare/byIsolate_asvoccurance_unmatched/", tempIsolate, "_leftover.png"), height=7, width=2+1.5*length(c(tempASV,tempASV_FILT)), limitsize=FALSE
         , gg_isotemp)
  ggsave(filename = paste0("03_manual_seq_compare/byIsolate_asvoccurance_unmatched/", tempIsolate, "_leftover_FILT.png"), height=7, width=2+1.5*length(c(tempASV,tempASV_FILT)), limitsize=FALSE
         , gg_isotemp_FILT)
}

### Most likely 37E!!!! Just not a very good grower, apparently. 
exactMatches <- seq_compare %>% 
  filter(IsolateID=="37E"& ASVRef%in% c("REF_70") 
  ) %>% 
  select(IsolateID, ASVRef, ASVID) %>%
  full_join(exactMatches)

######## Save exact matches and ASVIDs so we can built fasttree with the sequences
subsetASVs_for_tree <- exactMatches %>% select(ASVID, IsolateID, ASVRef) %>%
  full_join(alex_isos)

dir.create("03_manual_seq_compare/downstream")
write.table(subsetASVs_for_tree, file="03_manual_seq_compare/downstream/subsetASVs_for_tree_df.txt", sep="\t", quote=FALSE, row.names = FALSE)

# For downstream filtering; need to remove "ASV_"
asvs_for_downstream <- gsub("ASV_","",subsetASVs_for_tree$ASVID)
write.table(asvs_for_downstream, file="03_manual_seq_compare/downstream/subsetASVs_for_tree_list.txt", sep="\t", quote=FALSE, row.names = FALSE, col.names = "#SampleID")


setwd("..")
