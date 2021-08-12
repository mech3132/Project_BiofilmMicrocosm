#!bin/Rscript
##### Goals: Bind raw OTU data with metadata
library(tidyverse)
library(phytools)

setwd("02_externally_generated_data/16S_sequencing/")
source("../../code/inferring_traits_from_tree.R")
dir.create("downstream")
####### Load files ################
otu <- read.delim("intermediate_files/filtering_table/01_otu_table_SEPPfiltered_SILVA.txt", skip=1)
tree <- read.tree("intermediate_files/phylogeny/SEPP_inserted_tree_rooted_SILVA/tree.nwk")

woodhams <- read.delim("raw_data/woodhams_database/Amphibian-skin_bacteria_metadata.txt")
inhibIso <- read.delim("../../00_isolate_info/manual_chosen_isolates/Inhib_to_include.txt", header=FALSE) %>% rename(IsolateID=V1) %>% mutate(Inhibitory = "Inhib")
nonIso <- read.delim("../../00_isolate_info/manual_chosen_isolates/nonInhib_to_include.txt", header=FALSE)%>% rename(IsolateID=V1) %>% mutate(Inhibitory = "Non")
alexIsos <- rbind(inhibIso, nonIso)
alex_ids <- read.delim("intermediate_files/imported_databases/alexKeepSamesunique.txt") %>% 
  separate(col=X.SampleID, into = c("IsolateID","other"), sep = "-", remove=FALSE) %>%
  group_by(IsolateID) %>% mutate(dup = n()) %>% rowwise() %>%
  mutate(redo = ifelse(length(grep("_R_", other))>0, TRUE, FALSE)
         , keep = ifelse(dup==1, TRUE, ifelse(redo, TRUE, FALSE))) %>% ungroup() %>% filter(keep) %>%
  left_join(alexIsos) %>%
  rename(Isolate=IsolateID, IsolateID = X.SampleID) %>% select(IsolateID, Isolate, Inhibitory) %>%
  filter(!is.na(Inhibitory))
vsearch_taxa_SILVA <- read.delim("intermediate_files/bayes_classifier/SILVA_classifications_biofilm_alex_woodhams_merged/taxonomy.tsv")

blast_woodhams <- read.delim("intermediate_files/blast/blast_against_Woodhams_wheader.txt")
blast_alex <- read.delim("intermediate_files/blast/blast_against_Alex_wheader.txt")

# Create metadatafile of what should and shouldn't be inhibitory
allIsoInhib <- woodhams %>% rename(IsolateID = X.SampleID, Inhibitory = Bd_inhibition) %>% select(IsolateID, Inhibitory) %>% 
  mutate(Inhibitory = ifelse(Inhibitory == "inhibitory","Inhib", ifelse(Inhibitory =="ns", "Non", Inhibitory))) %>% filter(Inhibitory != "unknown") %>%
  full_join(alex_ids) 
# Get all sequences I want to keep
keepSeqs <- c(otu$X.OTU.ID, allIsoInhib$IsolateID)
# Remove those not in tree
keepSeqs <- keepSeqs[keepSeqs %in% tree$tip.label]
# Filter tree
tree.filt <- keep.tip(tree, keepSeqs)

# Get tip colors to match up with inhib/not inhib
ASV_metadata <- allIsoInhib[match(tree.filt$tip.label,allIsoInhib$IsolateID),] %>%
  rename(knownInhibitory = Inhibitory) %>%
  mutate(ASVID = tree.filt$tip.label) %>%mutate(tipColInhib = ifelse(is.na(knownInhibitory), "grey", ifelse(knownInhibitory == "Inhib", "purple",ifelse(knownInhibitory == "enhancing", "green","yellow")))) %>%
  mutate(dataset = ifelse(!is.na(Isolate), "alex", ifelse(!is.na(IsolateID), "woodhams","mine"))
         ,tipColDataset = ifelse(dataset=="alex","blue", ifelse(dataset == "woodhams","brown", "red"))) %>%
  mutate(Inhibitory01 = ifelse(knownInhibitory == "Inhib",1,0)) %>% select(ASVID, knownInhibitory, Inhibitory01, dataset, everything())
# View(ASV_metadata)
# Get vectors of tips from each dataset
allAlexTips <- ASV_metadata %>% filter(dataset=="alex") %>% pull(ASVID)
allWoodhamsTips <- ASV_metadata %>% filter(dataset=="woodhams") %>% pull(ASVID)
allMineTips <- ASV_metadata %>% filter(dataset=="mine") %>% pull(ASVID)

# Clean up Vsearch taxonomy
# Note: I turn all uncultured and metagenome entries into NA.
vsearch_taxa_SILVA_edited <- vsearch_taxa_SILVA %>% separate(Taxon, into = c("Domain","Phylum","Class","Order","Family","Genus","Species"), sep = "; ", remove=FALSE, fill="right") %>%
  rowwise() %>% 
  mutate(Species = ifelse(length(grep("uncultured|metagenome", Species))>0, "s__", Species)
         , Genus = gsub("uncultured","",Genus)
         , Family = gsub("uncultured","",Family)
         , Order = gsub("uncultured","",Order)
         , Class = gsub("uncultured","",Class)
         , Phylum = gsub("uncultured","",Phylum)) %>%
  mutate(temp_Domain = ifelse(Domain=="d__","Unknown_isolate",Domain)
         , temp_Phylum = ifelse(Phylum=="p__", ifelse(length(grep("isolate",temp_Domain))==0,paste(temp_Domain,"isolate", sep="_"), temp_Domain),Phylum)
         , temp_Class = ifelse(Class=="c__", ifelse(length(grep("isolate",temp_Phylum))==0,paste(temp_Phylum,"isolate", sep="_"), temp_Phylum),Class)
         , temp_Order = ifelse(Order=="o__", ifelse(length(grep("isolate",temp_Class))==0, paste(temp_Class,"isolate", sep="_"),temp_Class),Order)
         , temp_Family = ifelse(Family=="f__",ifelse(length(grep("isolate",temp_Order))==0, paste(temp_Order,"isolate", sep="_"),temp_Order),Family)
         , temp_Genus = ifelse(Genus=="g__", ifelse(length(grep("isolate",temp_Family))==0, paste(temp_Family,"isolate", sep="_"),temp_Family),Genus)
         , temp_Species = ifelse(Species=="s__", ifelse(length(grep("isolate",temp_Genus))==0, paste(temp_Genus,"isolate", sep="_"),temp_Genus),Species)) %>%
  mutate(temp_Domain = ifelse(is.na(temp_Domain), "Unassigned", temp_Domain)
         , temp_Phylum = ifelse(is.na(temp_Phylum), temp_Domain, temp_Phylum)
         , temp_Class = ifelse(is.na(temp_Class), temp_Phylum, temp_Class)
         , temp_Order = ifelse(is.na(temp_Order), temp_Class, temp_Order)
         , temp_Family = ifelse(is.na(temp_Family), temp_Order, temp_Family)
         , temp_Genus = ifelse(is.na(temp_Genus), temp_Family, temp_Genus)
         , ScientificName = ifelse(is.na(temp_Species), paste(temp_Genus,"unknown", sep = "_"), temp_Species)
  ) %>% 
  mutate(Domain = ifelse(Domain == "d__",NA, Domain)
         , Phylum = ifelse(Phylum == "p__",NA, Phylum)
         , Class = ifelse(Class == "c__",NA, Class)
         , Order = ifelse(Order == "o__",NA, Order)
         , Family = ifelse(Family == "f__",NA, Family)
         , Genus = ifelse(Genus == "g__",NA, Genus)
         , Species = ifelse(Species == "s__",NA, Species)
  )  %>%
  select(-starts_with("temp_")) %>% ungroup() %>%
  mutate(ScientificName_unique = make.unique(ScientificName))  %>% 
  rename(ASVID = Feature.ID)

# Identify top BLAST hits
blast_against_woodhams_summary <- blast_woodhams %>% filter(Length>=150) %>% 
  group_by(QueryID) %>% mutate(minevalue = min(evalue)) %>% ungroup() %>%
  filter(minevalue == evalue) %>% group_by(QueryID) %>% mutate(maxPercIdent = max(PercentIdentity)) %>% ungroup() %>%
  filter(maxPercIdent == PercentIdentity) %>%
  rename(X.SampleID=MatchedSequenceID) %>% left_join(woodhams) %>% rename(Matched_WoodhamsID = X.SampleID) %>%
  select(QueryID, Matched_WoodhamsID, PercentIdentity, Length, Bd_inhibition) %>% 
  filter(Bd_inhibition != "unknown") %>% mutate(Bd_inhibition = ifelse(Bd_inhibition == "inhibitory", 1, 0)) %>%
  group_by(QueryID) %>% summarize(avePercIdent_blast_woodhams = mean(PercentIdentity), inhibitory_probability_woodhams = mean(Bd_inhibition))

blast_against_alex_summary <- blast_alex %>% filter(Length>=150) %>% 
  group_by(QueryID) %>% mutate(minevalue = min(evalue)) %>% ungroup() %>%
  filter(minevalue == evalue) %>% group_by(QueryID) %>% mutate(maxPercIdent = max(PercentIdentity)) %>% ungroup() %>%
  filter(maxPercIdent == PercentIdentity) %>%
  rename(Matched_AlexID=MatchedSequenceID) %>% 
  separate(Matched_AlexID, into = c("IsolateID"), remove=FALSE, extra="drop") %>% 
  left_join(alexIsos) %>% 
  select(QueryID, Matched_AlexID, PercentIdentity, Length, Inhibitory) %>% mutate(Inhibitory = ifelse(Inhibitory == "Inhib", 1, 0)) %>%
  group_by(QueryID) %>% summarize(avePercIdent_blast_alex = mean(PercentIdentity), inhibitory_probability_alex = mean(Inhibitory), Matched_AlexID)
# Sanity check; are Alex's isolates making sense when blasted against each other?
blast_against_alex_summary %>% filter(QueryID %in% allAlexTips) %>% mutate(matched = c(QueryID!=Matched_AlexID)) %>% pull(matched) %>% any()
# Alex's blast always turns up her same sequences good.

## Infer all tips from reference tips
inhibitory_inference <- infer_binary_trait(test_isos = c(allMineTips, allAlexTips, allWoodhamsTips)
                                             , reference_isos = c(allAlexTips, allWoodhamsTips)
                                             , full_tree = tree.filt
                                             , metadata = ASV_metadata
                                             , sampleid_colname = "ASVID"
                                             , trait = "Inhibitory01"
                                             , max_cophen_dist = 0.5)

allInhibData <- inhibitory_inference %>% rename(QueryID = test_isos, best_matches_from_tree = best_matches, tree_cophen_dist = match_cophen_dist, inhibitory_probability_tree = probability) %>%
  full_join(blast_against_woodhams_summary) %>% full_join(blast_against_alex_summary) %>%
  rename(ASVID = QueryID) %>%
  select(ASVID, inhibitory_probability_alex, inhibitory_probability_woodhams, inhibitory_probability_tree, tree_cophen_dist, everything())  %>%
  full_join(ASV_metadata) %>%
  select(-best_matches_from_tree) # I think this is causing issues?
## Sanity check on known isolates
allInhibData %>% filter(dataset != "mine") %>% select(ASVID, Inhibitory01, dataset, inhibitory_probability_alex, inhibitory_probability_woodhams, inhibitory_probability_tree) %>%
  pivot_longer(cols= c("inhibitory_probability_alex", "inhibitory_probability_woodhams", "inhibitory_probability_tree"), names_to = "Metric", values_to = "PredictedInhibitory") %>%
  mutate(Error = abs(Inhibitory01 - PredictedInhibitory)) %>%
  ggplot() + geom_boxplot(aes(x=Metric, y=Error, col=dataset)) 
# Phylogenetic tree is NOT very good at identifying whether or not something is inhibitory

##### Prepare sequencing data for future merging of other datasets #####
otu_edit <- otu %>% rename(ASVID = X.OTU.ID)

write.tree(tree.filt, file = "downstream/filtered_tree.tre")
write.table(allInhibData, file = "downstream/allInhibData.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(otu_edit, file = "downstream/otu_table.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(vsearch_taxa_SILVA_edited, file = "downstream/taxonomy.txt", quote = FALSE, row.names = FALSE, sep = "\t")
setwd("../../")

