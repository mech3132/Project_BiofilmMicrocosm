#!bin/bash
library(tidyverse)
setwd("01_experiment_metadata/")
########## Aggregate experimental metadata #########

#### This aggregates all data ####
# Data fp
# rep1 <- "../Rep1b_25Aug2020/jar_metadata_rep1.txt"
rep12 <- "./Rep12_31Aug2020/jar_metadata_rep12_used.txt"
rep34 <- "./Rep34_7Sept2020/jar_metadata_rep34_used.txt"
rep56 <- "./Rep56_14Sept2020/jar_metadata_rep56_used.txt"
rep78 <- "./Rep78_19Sept2020/jar_metadata_rep78_used.txt"
rep910 <- "./Rep910_28Sept2020/jar_metadata_rep910_used.txt"
rep1112 <- "./Rep1112_5Oct2020/jar_metadata_rep1112_used.txt"


allinhib <- "../00_isolate_info/manual_chosen_isolates/Inhib_to_include.txt"
allnoninhib <- "../00_isolate_info/manual_chosen_isolates/nonInhib_to_include.txt"

#### Load files #####
# rep1_dat <- read.delim(file=rep1)
rep12_dat <- read.delim(file=rep12)
rep34_dat <- read.delim(file=rep34)
rep56_dat <- read.delim(file=rep56)
rep78_dat <- read.delim(file=rep78)
rep910_dat <- read.delim(file=rep910)
rep1112_dat <- read.delim(file=rep1112)
allIso <- c(read.delim(allinhib, header=FALSE)$V1, read.delim(allnoninhib, header=FALSE)$V1)

#### Aggregate ####
allReps <- rep12_dat %>% 
  full_join(rep34_dat) %>% 
  full_join(rep56_dat) %>% 
  full_join(rep78_dat) %>% 
  full_join(rep910_dat) %>%
  full_join(rep1112_dat)

#### Adjust column names ####
#### Edit ####
colnames(allReps) <- gsub("^X","",colnames(allReps))
allReps <- allReps %>%
  select(SampleID, JarID, RichLevel, Rep, DateStart, Inhibitory, Type, everything()) %>%
  arrange(Rep, Inhibitory, RichLevel, Type) %>%
  mutate(suffix=ifelse(Type=="Bd", "Bd", ifelse(Type=="CV","Pre",NA))
         ,Rep_n=gsub("^Rep","R",Rep)
         ,J=paste0("J",JarID)) %>%
  unite(Rep_n, J, suffix,col=sampleID_qPCR)

### Check for duplicates
allReps_noCon <- allReps %>%filter(Type=="Bd", Inhibitory !="Con")

whichDup <- duplicated(allReps_noCon[,which(colnames(allReps_noCon) %in% allIso)])
allReps_noCon[whichDup,]


# Get real introduced richness
introRich <- allReps %>% select(one_of(allIso)) %>% rowSums(na.rm=TRUE)
allReps <- allReps %>% mutate(introRich=introRich)

# Remove obsolete columns by including all the ones below:
allReps <- allReps %>% select(SampleID, JarID, RichLevel, introRich, Rep, DateStart, Inhibitory, Type, one_of(allIso), volPellet, cClean, sampleID_qPCR)

### Save aggregated metadata
dir.create("downstream")
write.table(allReps, file="downstream/aggregated_metadata.txt", quote=FALSE, row.names = FALSE, sep="\t")

setwd("..")

