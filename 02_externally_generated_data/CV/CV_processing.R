#!bin/bash

### Process CV results
library(tidyverse)
library(lme4)
setwd("02_externally_generated_data/CV/")
cv <- read.delim("raw_data/CV_all_16Oct2020.txt")

cv_cleaned <- cv %>% select(OriginalCV_platecode, JarID, Rep, Row, Column, OD_rep, CV) %>%
  mutate(blank = ifelse(Rep == "BLANK", 1,0)
         , well = paste0(Row,Column)
         , OD_rep = paste0("ODrep",OD_rep)
          , OriginalCV_platecode = ifelse(Rep=="BLANK", paste0(Rep,Row,Column), OriginalCV_platecode))

cv_agg <- cv_cleaned %>% group_by(OriginalCV_platecode, JarID, Rep, Row, Column) %>%
  summarize(aveOD590 = mean(CV, na.rm=TRUE)
            , nRep = length(CV)
            , sd_CV = sd(CV, na.rm=TRUE)) %>%
  mutate(blank = ifelse(Rep == "BLANK", 1,0)
         , OriginalCV_platecode = ifelse(Rep=="BLANK", paste0(Rep,Row,Column), OriginalCV_platecode))

blank_meansd <- cv_agg %>% filter(Rep == "BLANK") %>% ungroup() %>%
  group_by(Rep) %>% summarize(aveOD590 = mean(aveOD590)
                              , sd = sqrt(mean(sd_CV^2))) %>% ungroup() 

cv_agg_adj <- cv_agg %>% filter(Rep !="BLANK") %>%
  mutate(blank_mean = as.numeric(blank_meansd[,"aveOD590"]), blank_sd = as.numeric(blank_meansd[,"sd"])) %>%
  mutate(CV_adj = aveOD590-blank_mean
         , sd_CV_adj = sqrt(sd_CV^2 + blank_sd^2)) %>% ungroup() %>%
  select(JarID, Rep, CV_adj, sd_CV_adj) %>%
  mutate(CV_adj = ifelse(CV_adj<0, 0,CV_adj))

dir.create("downstream")
write.table(cv_agg_adj, file="downstream/CV_adj.txt", quote=FALSE, row.names=FALSE, sep="\t")

setwd("../..")
