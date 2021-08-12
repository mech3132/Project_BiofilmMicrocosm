#!bin/bash Rscript

library(ape)
library(vegan)
library(tidyverse)

setwd("03_datacleaning/")
dir.create("05_final_seq_filtering")

#### Load files ########
asv_and_meta <- read.delim("02_merge_all_data/metadata_allReps.txt")
# asv_raw <- read.delim("01_seq_data_cleaning/downstream/raw_otu_relative_abundance.txt")
# asv_meta <- read.delim("01_seq_data_cleaning/downstream/asv_metadata.txt")
filtSeq_meta <- read.delim("03_manual_seq_compare/downstream/subsetASVs_for_tree_df.txt")
treeFilt <- read.tree("04_tree_mapping_for_filtered_seqs/raxml_tree_filteredseqs.tre/tree.nwk")
allASV_inhibData <- read.delim("../00_isolate_info/all_isolate_info_combined.txt", sep = ",")
allASV_taxonomy <- read.delim("../02_externally_generated_data/16S_sequencing/downstream/taxonomy.txt")
otu_raw <- read.delim("../02_externally_generated_data/16S_sequencing/downstream/otu_table.txt", row.names = 1)
#### Filter to include only relavent
treeTips <- treeFilt$tip.label

# Fix header names for OTU
colnames(otu_raw) <- gsub(".","_", gsub("..","-",colnames(otu_raw), fixed=TRUE), fixed=TRUE)
otu_all <- otu_raw %>% t() %>% as.data.frame() %>%
  rownames_to_column(var="sampleID_qPCR") 

# Get filtered files
taxa_filt <- allASV_taxonomy %>% filter(ASVID %in% treeTips) %>% full_join(data.frame(ASVID=treeTips)) %>%
  rowwise() %>% mutate(dataset=ifelse(length(grep(".ab1", ASVID, fixed=TRUE))>0, "alex","mine")) %>% ungroup() %>%
  rowwise() %>% mutate(IsolateID = ifelse(dataset=="alex", strsplit(ASVID, split="-")[[1]], NA)) %>%
  ungroup()
allIsolates <- unique(taxa_filt$IsolateID[!is.na(taxa_filt$IsolateID)])
inhibData_filt <- allASV_inhibData %>% select(IsolateID, Inhibitory) %>% filter(IsolateID %in% allIsolates) %>%
  mutate(Inhibitory01 = ifelse(Inhibitory,1,0))

### Make metadata of filered sequences
ASV_meta_full <- taxa_filt %>% full_join(inhibData_filt) %>% filter(ASVID %in% treeTips) %>%
  mutate(dataset = ifelse(is.na(dataset), "alex", dataset)) %>% rowwise() %>%
  mutate(IsolateID = ifelse(dataset=="alex", strsplit(ASVID, split="-")[[1]], NA)) %>% ungroup() %>% 
  mutate(treeTipLabel = ifelse(dataset=="alex", IsolateID, ScientificName_unique)) %>% 
  select(ASVID, IsolateID, Taxon, Domain, Phylum, Class, Order, Family, Genus, Species, ScientificName_unique, dataset, Inhibitory01, treeTipLabel) %>%
  mutate(treePlotColor = ifelse(dataset=="mine", "black", ifelse(Inhibitory01==1, "red", "blue")))
# ASV_meta_full %>% View()

########### Try mapping ONLY important ASVs and then cross-checking ######
importantASVs <- gsub("ASV_","",asv_and_meta %>% select(starts_with("ASV")) %>% colnames())
alexASVs <- treeFilt$tip.label[grep(".ab1",treeFilt$tip.label, fixed=TRUE)]
nonImportantTips <- treeFilt$tip.label[!(treeFilt$tip.label %in% c(importantASVs,alexASVs))]
treeFilt_FILT_tipsChanged <- drop.tip(treeFilt, tip = nonImportantTips)

treeFilt_FILT_tipsChanged$tip.label <- pull(ASV_meta_full[match(treeFilt_FILT_tipsChanged$tip.label, ASV_meta_full$ASVID),"treeTipLabel"])
pdf("05_final_seq_filtering/raxml_tree_of_filtered_seqs_ONLYIMPORTANT.pdf", height=15, width=10)
plot(treeFilt_FILT_tipsChanged, tip.color = pull(ASV_meta_full[match(treeFilt_FILT_tipsChanged$tip.label, ASV_meta_full$treeTipLabel),"treePlotColor"]))
dev.off()

# And also print the supposed matching
toWrite_filt <- filtSeq_meta %>% left_join(ASV_meta_full %>% mutate(ASVID = paste0("ASV_", ASVID)) %>% select(ASVID, ScientificName_unique, dataset)) %>%
  filter(dataset=="mine") %>% select(-dataset) %>% arrange(IsolateID) %>%
  left_join(filtSeq_meta) %>%
  filter(ASVID %in% paste0("ASV_",importantASVs))
write.table(toWrite_filt, file="05_final_seq_filtering/table_to_crosscheck_with_tree_IMPORTANTONLY.txt", sep="\t", quote=FALSE, row.names = FALSE)



#### Upload after manual edits #####
treeMatches <- read.delim("05_final_seq_filtering/table_to_crosscheck_with_tree_IMPORTANTONLY_ANNOTATED.csv")
final_assignments <- treeMatches %>% filter(perfectmatch)
unassignedASVs <- importantASVs[!(paste0("ASV_",importantASVs) %in% final_assignments$ASVID)]
unassignedIsolates <- allIsolates[!(allIsolates %in% final_assignments$IsolateID)]
unassignedAlex <-ASV_meta_full %>% filter(IsolateID %in% unassignedIsolates) %>% pull(ASVID)
tipsToDrop <- treeFilt$tip.label[!(treeFilt$tip.label %in% c(unassignedASVs, unassignedAlex))]

##### Look at tree
treeFilt_tipsChanged_unassigned <- drop.tip(treeFilt, tip = tipsToDrop)
# NOTE: the one that I can't find is that extra in Bd set that "appeared" later on!

treeFilt_tipsChanged_unassigned$tip.label <- pull(ASV_meta_full[match(treeFilt_tipsChanged_unassigned$tip.label, ASV_meta_full$ASVID),"treeTipLabel"])
pdf("05_final_seq_filtering/raxml_tree_of_filtered_seqs_LEFTOVERS.pdf", height=10, width=15)
plot(treeFilt_tipsChanged_unassigned, tip.color = pull(ASV_meta_full[match(treeFilt_tipsChanged_unassigned$tip.label, ASV_meta_full$treeTipLabel),"treePlotColor"]))
dev.off()

# And also print the supposed matching
newList_nonassignedASVs <- ASV_meta_full %>% filter(dataset=="mine") %>% select(ASVID, ScientificName_unique) %>%
  filter(ASVID %in% unassignedASVs) %>%
  mutate(ASVID = paste0("ASV_",ASVID)) %>%
  left_join(filtSeq_meta) %>%
  rowwise() %>%
  mutate(IsolateID = ifelse(IsolateID %in% unassignedIsolates, IsolateID, NA))%>% ungroup()
write.table(newList_nonassignedASVs, file="05_final_seq_filtering/table_to_crosscheck_with_tree_UNASSIGNED.txt", sep="\t", quote=FALSE, row.names = FALSE)

#### Generate images of observed vs expected sampling schemes to compare
dir.create("05_final_seq_filtering/unassignedASV_maps")
for ( asv in unassignedASVs) {
  tempName <- ASV_meta_full %>% filter(ASVID==asv) %>% pull(ScientificName_unique)
  gg_tempasv <- asv_and_meta %>%
    # filter(Type=="CV", Inhibitory !="Con") %>% 
    filter( Inhibitory !="Con") %>%
    mutate(Type = factor(Type, levels=c("CV","Bd"))) %>%
    select(sampleID_qPCR, Rep, Inhibitory, RichLevel, Type) %>%
    left_join(otu_all) %>%
    select(Rep, Inhibitory, RichLevel, Type, one_of(asv)) %>%
    mutate(Rep = factor(Rep, levels=c(paste0("Rep",seq(1,12))))) %>%
    rename_at(vars(all_of(asv)), ~"Reads") %>%
    group_by(Type) %>%
    mutate(Reads = Reads/sum(Reads)) %>%
    mutate(Reads = ifelse(Reads==0, NA, Reads)) %>%
    ggplot() + geom_point(aes(x=factor(RichLevel), y=Rep, col=Reads), cex=4) +
    facet_grid(Inhibitory~Type) +xlab("RichLevel") + labs(title=paste0(tempName, ": ", asv))+
    scale_color_gradient(low="black", high="red", na.value = "lightgrey")
  ggsave(filename = paste0("05_final_seq_filtering/unassignedASV_maps/",asv,".png"), 
         height=6, width=6
         ,gg_tempasv)
  gg_tempasv2 <- asv_and_meta %>%
    # filter(Type=="CV", Inhibitory !="Con") %>% 
    filter( Inhibitory !="Con") %>%
    mutate(Type = factor(Type, levels=c("CV","Bd"))) %>%
    select(Rep, Inhibitory, RichLevel, Type, one_of(paste0("ASV_",asv))) %>%
    mutate(Rep = factor(Rep, levels=c(paste0("Rep",seq(1,12))))) %>%
    rename_at(vars(all_of(paste0("ASV_",asv))), ~"Reads") %>%
    group_by(Type) %>%
    mutate(Reads = Reads/sum(Reads)) %>%
    mutate(Reads = ifelse(Reads==0, NA, Reads)) %>%
    ggplot() + geom_point(aes(x=factor(RichLevel), y=Rep, col=Reads), cex=4) +
    facet_grid(Inhibitory~Type) +xlab("RichLevel") + labs(title=paste0(tempName, ": ", asv))+
    scale_color_gradient(low="black", high="red", na.value = "lightgrey")
  ggsave(filename = paste0("05_final_seq_filtering/unassignedASV_maps/",asv,"_FILT.png"), 
         height=6, width=6
         ,gg_tempasv2)
}

dir.create("05_final_seq_filtering/unassignedIso_maps")
actualIsolates <- gsub("ISO_","",asv_and_meta %>% select(starts_with("ISO")) %>% colnames())
unassignedIsolates2 <- unassignedIsolates[unassignedIsolates %in% actualIsolates]
for ( iso in unassignedIsolates2) {
  # tempName <- ASV_meta_full %>% filter(IsolateID==iso) %>% pull(ScientificName_unique)
  gg_tempiso <- asv_and_meta %>% filter(Type=="CV", Inhibitory !="Con") %>% 
    select(Rep, Inhibitory, RichLevel, one_of(paste0("ISO_",iso))) %>%
    mutate(Rep = factor(Rep, levels=c(paste0("Rep",seq(1,12))))) %>%
    rename_at(vars(all_of(paste0("ISO_",iso))), ~"Present") %>%
    mutate(Reads = ifelse(Present, TRUE, NA)) %>%
    ggplot() + geom_point(aes(x=factor(RichLevel), y=Rep, col=Reads), cex=5) +
    facet_wrap(Inhibitory~., ncol=1) +xlab("RichLevel") + labs(title=paste0(iso))
  ggsave(filename = paste0("05_final_seq_filtering/unassignedIso_maps/",iso,".png"), 
         height=6, width=4,gg_tempiso)
}


#### Investigate Bd settlement for strange cases
# ASV_meta_full %>%
#   filter(IsolateID%in%c("37E","34B","32E","2F","BTB_59","40A","2E","2H","50A","34F")) %>%
#   select(IsolateID, ScientificName_unique)
# 
# ASV_meta_full %>%
#   filter(IsolateID%in%allIsolates) %>%
#   select(IsolateID, ScientificName_unique) %>% View()

## Load manual notes for leftovers ####

# Note: manually added the pseudomonas back in that was only in Bd. Won't do any harm anyway, and pretty sure it's the same as all the other pseudomonases
annotate_leftovers <- read.csv("05_final_seq_filtering/table_to_crosscheck_with_tree_UNASSIGNED_ANNOTATED.csv")

final_assignments_all <- annotate_leftovers %>%
  mutate(IsolateID=proposed, perfectmatch=FALSE
         , phylogeneticmatch=ifelse(is.na(phylogeneticmatch), FALSE, TRUE)) %>%
  select(ASVID, IsolateID, ASVRef, ScientificName_unique, perfectmatch, phylogeneticmatch) %>%
  full_join(final_assignments)

#### Plots of final match-ups ####

#### Generate images of observed vs expected sampling schemes to compare
dir.create("05_final_seq_filtering/final_matchups")
# importantASVs_filt <- importantASVs[importantASVs!="14c6f745ca3896169903ea6ea03d48ac"]
for ( asv in importantASVs) {
  tempName <- ASV_meta_full %>% filter(ASVID==asv) %>% pull(ScientificName_unique)
  tempIso <- final_assignments_all %>% filter(ASVID==paste0("ASV_",asv)) %>% pull(IsolateID)
  tempPlan <- asv_and_meta %>% filter(Type=="CV", Inhibitory !="Con") %>% 
    select(Rep, Inhibitory, RichLevel, one_of(paste0("ISO_",iso))) %>%
    mutate(Rep = factor(Rep, levels=c(paste0("Rep",seq(1,12))))) %>%
    rename_at(vars(all_of(paste0("ISO_",iso))), ~"Present") %>%
    mutate(Reads = ifelse(Present, TRUE, NA))
  
  gg_tempmatch <- asv_and_meta %>%
    # filter(Type=="CV", Inhibitory !="Con") %>% 
    filter( Inhibitory !="Con") %>%
    mutate(Type = factor(Type, levels=c("CV","Bd"))) %>% 
    select(Rep, Inhibitory, RichLevel, Type, one_of(paste0("ASV_",asv)), one_of(paste0("ISO_",tempIso))) %>%
    mutate(Rep = factor(Rep, levels=c(paste0("Rep",seq(1,12))))) %>%
    rename_at(vars(all_of(paste0("ASV_",asv))), ~"Reads") %>%
    rename_at(vars(all_of(paste0("ISO_",tempIso))), ~"Present") %>%
    mutate(Present = ifelse(is.na(Present), FALSE, Present)) %>%
    # mutate(Present = ifelse(Present, 24, 21)) %>%
    group_by(Type) %>%
    mutate(Reads = Reads/sum(Reads)) %>%
    mutate(Reads = ifelse(Reads==0, NA, Reads)) %>%
    ggplot() + geom_point(aes(x=factor(RichLevel), y=Rep, fill=Reads, pch=Present), cex=4) +
    facet_grid(Inhibitory~Type) +xlab("RichLevel") + labs(title=paste0(tempName, ": ", asv))+
    scale_fill_gradient(low="black", high="red", na.value = "lightgrey")+
    scale_shape_manual(values=c(21, 24))
  ggsave(filename = paste0("05_final_seq_filtering/final_matchups/",tempIso,"_",asv,".png"), 
         height=6, width=6
         ,gg_tempmatch)
}


## Summarize findings:
# How many had perfect, single, matches?
PerfectSingleMatches_withphylogeny <- final_assignments_all %>% group_by(IsolateID) %>% mutate(number=n()) %>% ungroup() %>%
  filter(number==1, perfectmatch==TRUE, phylogeneticmatch==TRUE)
PerfectSingleMatches_nophylogeny <- final_assignments_all %>% group_by(IsolateID) %>% mutate(number=n()) %>% ungroup() %>%
  filter(number==1, perfectmatch==TRUE) %>%
  filter(!(IsolateID %in% PerfectSingleMatches_withphylogeny$IsolateID))
SingleMatches_withphylogeny <- final_assignments_all %>% group_by(IsolateID) %>% mutate(number=n()) %>% ungroup() %>%
  filter(number==1, phylogeneticmatch==TRUE) %>%
  filter(!(IsolateID %in% c(PerfectSingleMatches_withphylogeny$IsolateID, PerfectSingleMatches_nophylogeny$IsolateID)))
MultipleMatches_withphylogeny <- final_assignments_all %>% group_by(IsolateID) %>% mutate(number=n()) %>% ungroup() %>%
  filter(phylogeneticmatch==TRUE) %>%
  filter(!(IsolateID %in% c(PerfectSingleMatches_withphylogeny$IsolateID, PerfectSingleMatches_nophylogeny$IsolateID, SingleMatches_withphylogeny$IsolateID)))
Remaining <- final_assignments_all %>%
  filter(!(IsolateID %in% c(PerfectSingleMatches_withphylogeny$IsolateID, PerfectSingleMatches_nophylogeny$IsolateID, SingleMatches_withphylogeny$IsolateID, MultipleMatches_withphylogeny$IsolateID)))


##### Create final version of metadata #####
meta_only <- asv_and_meta %>% select(-starts_with("ASV"), -starts_with("ISO"), -obsRich, -REMOVE_SEQCONTAM)
ISO_only <- asv_and_meta %>% select( sampleID_qPCR, starts_with("ISO"))
ASV_only <- asv_and_meta %>% select( sampleID_qPCR, starts_with("ASV"))

# Keep ASVs as they are
ASV_raw <- ASV_only %>% rowwise() %>%
  mutate(obsRich_rawASV = sum(across(starts_with("ASV"))>0)
         , shannon_rawASV = diversity(across(starts_with("ASV")),index = "shannon")
         , pielou_rawASV = shannon_rawASV/log(obsRich_rawASV)) %>% 
  mutate(shannon_rawASV = ifelse(length(grep("Bd",sampleID_qPCR))>0, NA, shannon_rawASV)
         ,pielou_rawASV = ifelse(length(grep("Bd",sampleID_qPCR))>0, NA, pielou_rawASV)) %>% ungroup() 

# Merge ASVs into isolates and then get richness
ASV_merged_long <- ASV_only %>%
  pivot_longer(-sampleID_qPCR, names_to = "ASVID", values_to = "Reads") %>%
  left_join(final_assignments_all %>% select(ASVID, IsolateID)) %>%
  group_by(sampleID_qPCR, IsolateID) %>%
  summarise(sumReads = sum(Reads)) %>%
  ungroup() %>% mutate(IsolateID = paste0("ISO_", IsolateID))

# Get inhibitory Stats

meta_inhibitory <- inhibData_filt %>% select(IsolateID, Inhibitory01) %>% distinct() %>%
  mutate(IsolateID = paste0("ISO_",IsolateID)) %>% 
  right_join(ASV_merged_long) %>% 
  mutate(inhibRich_single = ifelse(sumReads>0, Inhibitory01, 0), inhibRich_prop = sumReads*Inhibitory01, PA = sumReads>0) %>%
  group_by(sampleID_qPCR) %>% 
  summarise(inhibRich = sum(inhibRich_single), totalRich = sum(PA), inhibReads = sum(inhibRich_prop), inhibProp = sum(inhibRich_prop)/sum(sumReads)) %>%
  mutate(inhibRichFrac = inhibRich/totalRich, inhibProp = ifelse(totalRich==0, 0, inhibProp), inhibRichFrac = ifelse(totalRich==0, 0, inhibRichFrac) ) %>%
  select(-totalRich)

ASV_merged <- ASV_merged_long %>%
  pivot_wider(names_from=IsolateID, values_from = sumReads) %>%
  rowwise() %>%
  mutate(obsRich_mergedASV = sum(across(starts_with("ISO"))>0)
         , shannon_mergedASV = diversity(across(starts_with("ISO")),index = "shannon")
         , pielou_mergedASV = shannon_mergedASV/log(obsRich_mergedASV)) %>% 
  mutate(shannon_mergedASV = ifelse(length(grep("Bd",sampleID_qPCR))>0, NA, shannon_mergedASV)
         ,pielou_mergedASV = ifelse(length(grep("Bd",sampleID_qPCR))>0, NA, pielou_mergedASV)) %>% ungroup() 

# Merge ASVs into isolates, remove the samples that should be zero, and then get richness
ASV_merged_filt <- ISO_only %>% pivot_longer(-sampleID_qPCR, names_to = "IsolateID", values_to = "Expected") %>%
  mutate(Expected = ifelse(is.na(Expected), FALSE, Expected)) %>%
  left_join(ASV_merged_long) %>%
  mutate(sumReads = ifelse(is.na(sumReads), 0, sumReads)) %>%
  mutate(finalReads = Expected*sumReads) %>%
  select(sampleID_qPCR, IsolateID, finalReads) %>%
  pivot_wider(names_from=IsolateID, values_from = finalReads)%>% rowwise() %>%
  mutate(obsRich_mergedASVfilt = sum(across(starts_with("ISO"))>0)
         , shannon_mergedASVfilt = diversity(across(starts_with("ISO")),index = "shannon")
         , pielou_mergedASVfilt = shannon_mergedASVfilt/log(obsRich_mergedASVfilt)) %>% 
  mutate(shannon_mergedASVfilt = ifelse(length(grep("Bd",sampleID_qPCR))>0, NA, shannon_mergedASVfilt)
         ,pielou_mergedASVfilt = ifelse(length(grep("Bd",sampleID_qPCR))>0, NA, pielou_mergedASVfilt)) %>% ungroup() 
  

###### Save for downstream #######
dir.create("05_final_seq_filtering/downstream")
write.table(ASV_raw, file="05_final_seq_filtering/downstream/ASV_only_raw.txt", sep="\t", quote=FALSE, row.names = FALSE)
write.table(ASV_merged, file="05_final_seq_filtering/downstream/ASV_only_merged_by_iso.txt", sep="\t", quote=FALSE, row.names = FALSE)
write.table(ASV_merged_filt, file="05_final_seq_filtering/downstream/ASV_only_merged_by_iso_filt.txt", sep="\t", quote=FALSE, row.names = FALSE)
write.table(meta_inhibitory, file="05_final_seq_filtering/downstream/metadata_inhibitory_summary.txt", sep="\t", quote=FALSE, row.names = FALSE)
allData_merged <- meta_only %>% full_join(ASV_merged) %>% full_join(meta_inhibitory)
allData_merged_filt <- meta_only %>% full_join(ASV_merged_filt) 

write.table(allData_merged, file="05_final_seq_filtering/downstream/asv_and_meta_merged.txt", sep="\t", quote=FALSE, row.names = FALSE)
write.table(allData_merged_filt, file="05_final_seq_filtering/downstream/asv_and_meta_merged_and_filt.txt", sep="\t", quote=FALSE, row.names = FALSE)

## Merge CV and Bd and collapse

### separate out; and then merge
allData_merged_CV <- allData_merged %>% filter(Type=="CV") %>% 
  rename(qPCR_log10_bact = qPCR_log10, qPCR_bact_sd = qPCR_sd, qPCR_bact_upr_95_total = qPCR_upr_95_total, qPCR_bact_lwr_95_total = qPCR_lwr_95_total
         , B_JarID = JarID, B_obsRich=obsRich_mergedASV, B_shannon = shannon_mergedASV, B_pielou = pielou_mergedASV, B_inhibRich = inhibRich, B_inhibRichFrac = inhibRichFrac, B_inhibProp = inhibProp, B_inhibReads = inhibReads) %>%
  select(-c(sampleID_qPCR, cClean, volPellet, Type, starts_with("Micro"), starts_with("ISO_"), clouds_present, pinpricks_present))

allData_merged_Bd <- allData_merged %>% filter(Type=="Bd") %>% 
  rename(qPCR_log10_Bd = qPCR_log10, qPCR_Bd_sd = qPCR_sd, qPCR_Bd_upr_95_total = qPCR_upr_95_total, qPCR_Bd_lwr_95_total = qPCR_lwr_95_total
         , A_JarID = JarID, A_obsRich=obsRich_mergedASV, A_inhibRich = inhibRich, A_inhibRichFrac = inhibRichFrac) %>%
  select(-c(sampleID_qPCR, cClean, volPellet, Type, starts_with("CV"), starts_with("ISO_")))


allData_merged_collapsed <- full_join(allData_merged_CV, allData_merged_Bd)
write.table(allData_merged_collapsed, quote=FALSE, sep="\t", row.names = FALSE, file="05_final_seq_filtering/downstream/asv_and_meta_merged_collapsed.txt")

### Save metadata for isolates
isolates_kept <- gsub("ISO_","",allData_merged %>% select(starts_with("ISO")) %>% colnames())
isolates_remaining <- allASV_inhibData %>% filter(IsolateID %in% isolates_kept) %>% distinct() 
write.table(isolates_remaining, quote=FALSE, sep="\t", row.names = FALSE, file="05_final_seq_filtering/downstream/isolate_info_filtered.txt")

###### Quick EDA to check things

allData_merged_collapsed %>% 
  ggplot(aes(x=B_obsRich, y=A_obsRich, col=Inhibitory)) +
  geom_jitter(width=0.2, height=0.2)

allData_merged_collapsed %>% 
  ggplot(aes(x=B_inhibRich, y=A_inhibRich, col=Inhibitory)) +
  geom_jitter(width=0.2, height=0.2)

allData_merged_collapsed %>% 
  ggplot() + geom_point(aes(x=B_obsRich, y=CV))

allData_merged_collapsed %>% 
  ggplot() + geom_jitter(aes(x=B_obsRich, y=qPCR_log10_bact, col=log(CV+0.01)), cex=3, width=0.2, height=0)

allData_merged_collapsed %>% 
  ggplot() + geom_jitter(aes(x=B_obsRich, y=B_shannon), cex=3, width=0.2, height=0)
allData_merged_collapsed %>% 
  ggplot() + geom_jitter(aes(x=B_obsRich, y=B_pielou), cex=3, width=0.2, height=0)

allData_merged_collapsed %>%
  ggplot() + geom_point(aes(x=B_obsRich, y=qPCR_log10_Bd, col=B_inhibRichFrac), cex=3)
allData_merged_collapsed %>%
  ggplot() + geom_point(aes(x=B_obsRich, y=qPCR_log10_Bd, col=log10(B_inhibReads+1)), cex=3)
allData_merged_collapsed %>%
  ggplot() + geom_point(aes(x=B_obsRich, y=qPCR_log10_Bd, col=B_inhibProp), cex=3)

allData_merged_collapsed %>%
  ggplot() + geom_point(aes(x=B_inhibRich, y=qPCR_log10_Bd), cex=3)
allData_merged_collapsed %>%
  ggplot() + geom_point(aes(x=B_inhibRichFrac, y=qPCR_log10_Bd), cex=3)
allData_merged_collapsed %>%
  ggplot() + geom_point(aes(x=B_inhibReads, y=qPCR_log10_Bd), cex=3)
allData_merged_collapsed %>%
  ggplot() + geom_point(aes(x=B_inhibProp, y=qPCR_log10_Bd), cex=3)

### TO BE CONTINUED...
allData_merged_collapsed %>% filter(RichLevel!=0) %>%
  ggplot(aes(x=log(CV+0.01), y=qPCR_log10_Bd, col=Inhibitory)) + geom_point(cex=3)+
  geom_smooth(method="lm")
allData_merged_collapsed %>% filter(RichLevel!=0) %>%
  ggplot() + geom_point(aes(col=log(CV+0.01), y=qPCR_log10_Bd, x=B_inhibRich), cex=3)
setwd("..")



