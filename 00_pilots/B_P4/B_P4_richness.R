#### Biofilm 4- richness ####
library(tidyverse)
dat <- read.delim("B_P4_planning_richness.txt")
bf <- read.delim("../Biofilm_assay/biofilm_assay_allisolates.txt")

# change headers
colnames(dat) <- gsub("X","",colnames(dat))

biofilmthickness <- bf[match(colnames(dat), bf$IsolateID),"CV"]
BFT <- data.frame(Isolate=colnames(dat), biofilmthickness=biofilmthickness) %>%
  filter(!is.na(biofilmthickness)) 

# Multiply
for ( i in BFT$Isolate) {
  dat[,i] <- dat[,i]*BFT[which(BFT$Isolate==i),"biofilmthickness"]
}

dat$estBFT <- dat %>%
  select(one_of(as.character(BFT$Isolate))) %>%
  rowSums()


####### Plotting #######

dat %>%
  ggplot() +geom_point(aes(x=PreCount, y=PostCount))
dat %>%
  ggplot() +geom_point(aes(x=PlatedR, y=PostR))
dat %>%
  ggplot() +geom_point(aes(x=Richness, y=PostR))
dat %>%
  ggplot() +geom_point(aes(x=Richness, y=PlatedR))


dat %>%
  ggplot() +geom_point(aes(x=PostR, y=CV))
dat %>%
  ggplot() +geom_point(aes(x=Richness, y=CV))
dat %>%
  ggplot() +geom_point(aes(x=estBFT, y=CV))


  
