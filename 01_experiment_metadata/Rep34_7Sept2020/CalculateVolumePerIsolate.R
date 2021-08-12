#!bin/bash

library(tidyverse)
setwd("Rep34_7Sept2020/")
#### For calculating the volume per tube ####

metadata <- read.delim("forCultureMixing.txt", row.names = 1)
metadata2 <- read.delim("jar_metadata_rep34_used.txt")
colnames(metadata) <- gsub("X","",colnames(metadata))
## Using 2F_WRONG 

### Subtract the tryptone controls

PR_raw <- read.delim("OD_isolates_rep34.txt", stringsAsFactors = FALSE)
# remove 20D
PR_raw <- PR_raw %>% filter(IsolateID!="20D")

trypt <- PR_raw %>% filter(IsolateID=="TRYP") %>% pull(OD) %>% mean()

## Summarize by rep
PR <- PR_raw %>% group_by(IsolateID) %>% 
  filter(IsolateID!="TRYP") %>%summarize(aveOD=mean(OD)-trypt, diff = abs(diff(OD)))

### Check to make sure readings are consistent between reps
PR %>% mutate(concerning = ifelse(diff > 0.2, TRUE, FALSE)) %>%
  mutate(too_dil = ifelse(aveOD<=0, TRUE, FALSE)) %>%
  select(IsolateID, too_dil, concerning)%>% View()

volume <- PR %>% mutate(volPerMl = 0.05/aveOD) %>%
  mutate(volPerMl = ifelse(volPerMl<0, 1, volPerMl)) %>%
  select(IsolateID, aveOD, volPerMl) 
write.table(volume, file="OD_summary.txt", row.names=FALSE, quote=FALSE, sep="\t")

volume <- volume %>% select(IsolateID, volPerMl)
# Collapse by jar
extrameta <- metadata2 %>% select(SampleID, JarID, Type) %>% spread(key=Type, value=JarID) %>%
  mutate(Bd = paste0("Jar_",Bd), CV = paste0("Jar_",CV))
  
# volumes <- metadata %>% t() %>% data.frame() %>% rownames_to_column(var="IsolateID") %>%
#   full_join(volume) %>% mutate_each(funs(.*volPerMl*2), starts_with("Jar_"))

volumes <- metadata %>% t() %>% data.frame() %>% rownames_to_column(var="IsolateID") %>%
  full_join(volume) %>% mutate_each(~.*volPerMl*2, starts_with("Jar_"))

combinedVolumes <- data.frame(IsolateID = volumes$IsolateID, volPerMl = volumes$volPerMl)
for ( r in 1:nrow(extrameta)) {
  twoJars <- as.vector(unlist(extrameta[r,c("Bd","CV")]))
  newColName <- paste0(twoJars[1], twoJars[2])
  newCol <- data.frame(rowSums(volumes[,c(twoJars)]))
  colnames(newCol) <- newColName
  combinedVolumes <- cbind(newCol, combinedVolumes)
}
combinedVolumes <- combinedVolumes %>% select(IsolateID, everything())
# 
# volumes <- metadata %>% t() %>% data.frame() %>% rownames_to_column(var="IsolateID") %>%
#   full_join(volume) %>% mutate_each(funs(.*volPerMl*2), starts_with("Jar_"))

combinedVolumes
water_vol <- 4-colSums(combinedVolumes[,-which(colnames(combinedVolumes)%in%c("IsolateID","volPerMl"))], na.rm = TRUE)
finalVolumes <- rbind(combinedVolumes, c(IsolateID="Water", water_vol, volPerMl = NA))

write.table(finalVolumes, file="mixingVolumes.txt", sep="\t", row.names = FALSE, quote = FALSE)

# 
# finalVolumes[,-which(colnames(finalVolumes)%in%c("IsolateID","volPerMl"))] %>%
#   mutate_all(funs(as.numeric(.))) %>% colSums


setwd("..")

