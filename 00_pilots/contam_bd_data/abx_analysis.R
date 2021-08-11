#!bin/bash

# setwd("Documents/PhD/Project_biofilm/contam_bd_data/")
library(tidyverse)


dat <- read.delim("plate_layout_abxtest.txt", na.strings = "")
dat_adj <- dat %>%
  mutate(Content=ifelse(is.na(Content), "0BL0",Content)) %>%
  separate(Content, into=c("bact","treat","rep"), sep=c(1,3)) %>%
  filter(treat!="BL") %>%
  mutate(Bacteria=ifelse(bact==0, "No Bacteria", ifelse(bact==1, "13C: Hafnia alvei", ifelse(bact==2, "20A: B. wiedmannii", "25B: P. weihenstephanensis")))
         , Treatment=ifelse(treat=="AB","Full ABX",ifelse(treat=="DA","Diluted ABX", ifelse(treat=="DI","0.5% T",ifelse(treat=="PW","Plate-washed","1% T"))))) 
dat_adj2 <- dat_adj %>%
  filter(bact==0) %>%
  group_by(Bacteria, Treatment) %>%
  summarize(Con_average=mean(OD600)) %>% ungroup() %>%
  select(Treatment, Con_average) %>%
  right_join(dat_adj) %>% 
  mutate(OD600_adj = OD600-Con_average) %>%
  mutate(Treatment = factor(Treatment, levels=c("Full ABX"
                                                ,"Diluted ABX"
                                                ,"Plate-washed"
                                                , "0.5% T"
                                                , "1% T"))) 

dat_bars <- dat_adj2%>%
  group_by(Bacteria, Treatment) %>%
  summarize(mean_OD = mean(OD600_adj), se_OD = sd(OD600_adj)/sqrt(n())) %>% ungroup()
  
ggsave(filename = "ABX_results.pdf", height=8, width=4
       ,ggplot() +
         geom_col(data=dat_bars, aes(x=Treatment, y=mean_OD))+
         geom_segment(data=dat_bars,aes(x=Treatment, xend=Treatment, y=mean_OD-se_OD, yend=mean_OD+se_OD)) +
         geom_jitter(data=dat_adj2,aes(x=Treatment, y=OD600_adj),alpha=0.5, width=0.1, height=0, col="red") +
         facet_grid(Bacteria~.) + theme(axis.text.x = element_text(angle=90))+
         ylab("OD600 (standardized by controls)"))

