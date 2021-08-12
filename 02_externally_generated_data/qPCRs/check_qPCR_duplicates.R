#!bin/bash

# Checking qPCR consistency between reps

library(tidyverse)
setwd("02_externally_generated_data//qPCRs")
dir.create("EDA_qPCR")

qdat <- read.delim("raw_data/qPCR_results_combined.txt") 
mdat <- read.delim("raw_data/platemaps/qPCR_platemaps_filled.txt")
# diss1 <- read.delim("qPCR_data/dissociation_curve_Bact_Col3-12_19Nov2020 - Instrument Data - Text Format 1.txt")
dat <- left_join(qdat, mdat) %>% filter(Run != "Bact_1") %>% #unite(Run, sampleID, col = same_run_duplicates, remove=FALSE) %>%
  mutate(RunRep = ifelse(Run == "Bact_1-2", "Run1", ifelse(Run=="Bact_3-12_Rep1","Run2",ifelse(Run=="Bact_3-12_Rep2","Run3",ifelse(Run=="Bd_lastrow","Run1",ifelse(Run=="Bd_1-11_Rep1","Run2",ifelse(Run=="Bd_1-11_Rep2","Run3",ifelse(Run=="Bd_trips", "Run4", ifelse(Run=="Bact_trips", "Run4", NA))))))))) %>%
  unite(RunRep, replicate_within_plate, col=fullReps, remove=FALSE)
## Need to remove curves that look wrong
## Look at prelimiary standard curves; remove those that don't look correct.
dat %>% filter(well_designation == "gblock_standard" ) %>%
  ggplot() + geom_point(aes(x=Ct_dR, y = log2(Quantity_copies), col=sampleID)) +
  facet_wrap(.~Run)
# Bact_1-2; Gblock-9
# Bact 3-12_Rep1; Gblock-9
# Bact 3-12_Rep2; Gblock-9
# Bact_trips; Gblock-8, Gblock-9
# Bd_1-11_Rep1; Gblock-10
# Bd_lastrow; Gblock-1
## make dataframe of values to renove
removeStandards <- data.frame(Run = c("Bact_1-2"
                                      , "Bact_3-12_Rep1"
                                      , "Bact_3-12_Rep2"
                                      , "Bact_trips"
                                      , "Bact_trips"
                                      , "Bd_1-11_Rep1"
                                      , "Bd_lastrow")
                             , sampleID = c(rep("Gblock e-9", 4)
                                           , "Gblock e-8"
                                           , "Gblock e-10"
                                           , "Gblock e-1")
                             , remove=TRUE)

dat_filt <- dat %>% left_join(removeStandards) %>% 
  mutate(remove = ifelse(is.na(remove), FALSE, remove)) %>%
  filter(!remove) %>% select(-remove)

## Manually calculate G-block standard adjustments for each plate. 
best_fit_standards <- dat_filt %>% filter(well_designation == "gblock_standard" ) %>%
  group_by(Run) %>% 
  separate(sampleID, into= c("A","temp"),sep="-", remove=FALSE) %>%
  mutate(temp = as.numeric(temp)
         , Quantity_copies_logged = log10(1.95*10^(10-temp))) %>%
  select(Run, Ct_dR, Quantity_copies_logged, Slope_dR, Threshold_dR) %>%
  mutate(slope = lm(Ct_dR ~ Quantity_copies_logged)$coefficients[2]
         , intercept = lm(Ct_dR ~ Quantity_copies_logged)$coefficients[1]
         , R2 = summary(lm(Ct_dR ~ Quantity_copies_logged))$adj.r.squared) %>% ungroup() %>%
  select(Run, slope, intercept, R2) %>%distinct()

#### Adjust Bd_lastrow to use G-block standard ####
dat_corrected <- dat_filt %>% left_join(best_fit_standards) %>% rowwise(.) %>%
  mutate(log10_CopyNumber = (Ct_dR-intercept)/slope
         ,CopyNumber= 10^log10_CopyNumber) %>%
  rename(QC_MxPro = Quantity_copies) 
  


#### Check NTC to see if there are any consistent cases of contamination
dat_corrected %>% filter(well_designation == "NTC") %>% 
  ggplot() + geom_point(aes(x=Run, y=CopyNumber)) + theme(axis.text.x = element_text(angle=90))
dat_corrected %>% 
  ggplot() + geom_point(aes(x=well_designation, y=log10_CopyNumber))
# Looks like there's reagent contamination-- subtract from values.
contam <- dat_corrected %>% filter(well_designation == "NTC", log10_CopyNumber>0) %>%
  group_by(Run) %>% summarize(baseline_contam= mean(log10_CopyNumber)) 
  

dat_corrected <- dat_corrected %>% left_join(contam) %>% 
  mutate(baseline_contam = ifelse(is.na(baseline_contam), 0, baseline_contam)) %>%
  mutate(log10_CopyNumber_adj = log10_CopyNumber-baseline_contam)

## Check that corrected values correspond to original
dat_corrected %>%
  ggplot() + geom_point(aes(x=log10(QC_MxPro), y=log10_CopyNumber))
dat_corrected %>%
  ggplot() + geom_point(aes(x=log10(QC_MxPro), y=log10_CopyNumber_adj))
# Not sure whether I should adjust for baseline NTC contam...  I think I will NOT use the contam-adj values.

### Look at lowest standard curve threshold
dat_corrected %>%
  # filter(well_designation=="gblock_standard") %>%
  group_by(Run, well_designation) %>%
  summarize(lowest = min(log10_CopyNumber, na.rm=TRUE)) %>%
  pivot_wider(names_from=well_designation, values_from=lowest)
dat_corrected %>%
  group_by(primers, well_designation) %>%
  summarize(lowest = min(log10_CopyNumber, na.rm=TRUE)) %>%
  pivot_wider(names_from=well_designation, values_from=lowest)

#### Plotting #####

#### Standard curves check ####
ggsave(filename = "EDA_qPCR/standardCurves.pdf", width=6, height=4
  ,dat_corrected %>% filter(well_designation == "gblock_standard") %>% 
  # mutate(RunRep = ifelse(Run == "Bact_1-2", "Run1", ifelse(Run=="Bact_3-12_Rep1","Run2",ifelse(Run=="Bact_3-12_Rep2","Run3",ifelse(Run=="Bd_lastrow","Run1",ifelse(Run=="Bd_1-11_Rep1","Run2",ifelse(Run=="Bd_1-11_Rep2","Run3",NA))))))) %>%
  ggplot(aes(x=log10_CopyNumber, y=Ct_dR, col=RunRep)) +
  geom_point() +
  geom_smooth(method="lm") +
  facet_grid(.~primers, drop = TRUE, scales = "free") + theme_bw()+
    xlab("Log 10 Copy number")+ylab("Cycle at threshold")+
    labs(title="Gblock standard curves"))
ggsave(filename = "EDA_qPCR/standardCurves_adj.pdf", width=6, height=4
       ,dat_corrected %>% filter(well_designation == "gblock_standard") %>% 
         # mutate(RunRep = ifelse(Run == "Bact_1-2", "Run1", ifelse(Run=="Bact_3-12_Rep1","Run2",ifelse(Run=="Bact_3-12_Rep2","Run3",ifelse(Run=="Bd_lastrow","Run1",ifelse(Run=="Bd_1-11_Rep1","Run2",ifelse(Run=="Bd_1-11_Rep2","Run3",NA))))))) %>%
         ggplot(aes(x=log10_CopyNumber_adj, y=Ct_dR, col=RunRep)) +
         geom_point() +
         geom_smooth(method="lm") +
         facet_grid(.~primers, drop = TRUE, scales = "free") + theme_bw()+
         xlab("Log 10 Copy number")+ylab("Cycle at threshold")+
         labs(title="Gblock standard curves"))

# dat_corrected %>% filter(well_designation == "gblock_standard") %>% 
#   # mutate(RunRep = ifelse(Run == "Bact_1-2", "Run1", ifelse(Run=="Bact_3-12_Rep1","Run2",ifelse(Run=="Bact_3-12_Rep2","Run3",ifelse(Run=="Bd_lastrow","Run1",ifelse(Run=="Bd_1-11_Rep1","Run2",ifelse(Run=="Bd_1-11_Rep2","Run3",NA))))))) %>%
#   ggplot(aes(x=log10_CopyNumber_adj, y=Ct_dR, col=RunRep)) +
#   geom_point() +
#   geom_smooth(method="lm") +
#   facet_grid(.~primers, drop = TRUE, scales = "free") + theme_bw()+
#   xlab("Log 10 Copy number")+ylab("Cycle at threshold")+
#   labs(title="Gblock standard curves")

ggsave(filename = "EDA_qPCR/Rsq_stCurves.pdf", height=3, width=4
       , dat_corrected %>% select(primers, RunRep, R2) %>% distinct() %>%
  ggplot() +   geom_boxplot(aes(x=primers, y=R2)) + 
  geom_point(aes(x=primers, y=R2, col=RunRep)) +
  ylab("RSq of Standard Curve") + xlab("Primer set"))

dat_corrected %>% select(primers, Run,slope) %>% distinct() %>%
  ggplot() + 
  geom_boxplot(aes(x=primers, y=slope)) +
  geom_point(aes(x=primers, y=slope, col=Run))

dat_corrected %>% select(primers, Run,Efficiency_percent) %>% distinct() %>%
  ggplot() + geom_point(aes(x=primers, y=Efficiency_percent, col=Run)) +
  geom_boxplot(aes(x=primers, y=Efficiency_percent))

dat_corrected %>% filter(is.na(Ct_dR), well_designation=="sample") %>% as.data.frame()
# See if these samples were low-biomass to begin with, and see their dup/triplicates
dat_corrected[dat_corrected$sampleID %in% c(dat_corrected %>% filter(is.na(Ct_dR), well_designation=="sample")  %>% pull(sampleID)),] %>%
  arrange(sampleID) %>% as.data.frame()

## Quick and dirty check
dat_variance <- dat_corrected %>% filter(well_designation == "sample" ) %>% 
  mutate(edge = (Row %in% c("A","H") | Column %in% c(1,12))) %>% # indicate if it is an edge case
  group_by(sampleID) %>% mutate(anyedge = any(edge), anyevap = any(evaporation)) %>% ungroup() %>% # are any duplicates evaporated?
  select(sampleID, log10_CopyNumber_adj, fullReps, anyedge, anyevap, primers) %>% 
  pivot_wider(names_from = fullReps, values_from=log10_CopyNumber_adj) %>%
  # mutate(Extra = ifelse(is.na(R3), FALSE, TRUE)) %>%
  rowwise(.) %>% mutate(sd = sd((c(Run1_1,Run1_2,Run2_1,Run2_2, Run3_1, Run3_2, Run4_1, Run4_2)), na.rm=TRUE)
                        , ave_log10copies = mean((c(Run1_1,Run1_2,Run2_1,Run2_2, Run3_1, Run3_2, Run4_1, Run4_2)), na.rm=TRUE)
                        , nRep = sum(!is.na(c(Run1_1,Run1_2,Run2_1,Run2_2, Run3_1, Run3_2, Run4_1, Run4_2))))
 
# First, check ones that have >2 replicates and see if they are consistent with each other.
dat_variance %>% 
  filter(nRep>2) %>% select(-anyedge,-anyevap, -primers)
# Omit outliers by summarizing within-sample variability, and then omitting those that go beyond 2SD if >2 replicates
sumsq_n <- dat_variance %>% 
  mutate(sumsq = sum(((c(Run1_1,Run1_2,Run2_1,Run2_2, Run3_1, Run3_2, Run4_1, Run4_2))-ave_log10copies)^2, na.rm=TRUE)
         ,nRep = sum(!is.na(c(Run1_1,Run1_2,Run2_1,Run2_2, Run3_1, Run3_2, Run4_1, Run4_2)))) %>%
  select(sumsq, nRep) %>%
  colSums()
totalSD <- sqrt(sumsq_n[1]/(sumsq_n[2]-1))
# Now, use this to calculate lower and upper variance limits.
dat_outliers <- dat_variance %>%
  mutate(totalStDev = totalSD
         , lwr_95 = ave_log10copies - 2*totalStDev
         , upr_95 = ave_log10copies + 2*totalStDev) %>% 
  pivot_longer(cols = c(Run1_1,Run1_2,Run2_1,Run2_2, Run3_1, Run3_2, Run4_1, Run4_2), names_to = "fullReps", values_to = "log10_CopyNumber_adj") %>% 
  mutate(outlier = (log10_CopyNumber_adj<lwr_95 | log10_CopyNumber_adj>upr_95)) %>%
    select(sampleID, ave_log10copies, lwr_95, upr_95, fullReps, log10_CopyNumber_adj, outlier)%>%
  filter(!is.na(log10_CopyNumber_adj))
# Look at all the replicates, and see if there enough to consider not re-doing them
outlier_meta <- dat_outliers %>% filter(outlier)
outlier_samples <- unique(pull(dat_outliers[which(dat_outliers$outlier),"sampleID"]))
# Get lwr and uppr bounds
lwr_uppr <- dat_outliers %>% select(sampleID, lwr_95, upr_95) %>% distinct()

dat_variance %>%left_join(lwr_uppr) %>% filter(sampleID %in% outlier_samples|nRep==1) %>% select(-c(primers, anyedge))
#### This includes the re-do triplicates (and the singles that lost their duplicate)-- now we consider whether
# we cant to remove outliers or pool them all together. 

# Omit Run 4 for R12_J22 and R6_J7-- those seem to be outliers in and of themselves. 
# Omit Run 2 for R2_J17_Pre, R2_J6_Bd; that seems to be the outliers.
# Average all 4 for R12_J3-- I think it is just very low abundance.

dat_variance_adj <- dat_variance %>%
  mutate(Run4_1 = ifelse(sampleID %in% c("R12_J22_Pre", "R6_J7_Bd"), NA, Run4_1)
         , Run4_2 = ifelse(sampleID %in% c("R12_J22_Pre", "R6_J7_Bd"), NA, Run4_2)
         , Run2_1 = ifelse(sampleID %in% c("R2_J17_Pre", "R2_J6_Bd"), NA, Run2_1)) %>%
  mutate(Run1_ave = mean(c(Run1_1, Run1_2), na.rm=TRUE)
         , Run2_ave = mean(c(Run2_1, Run2_2), na.rm=TRUE)
         , Run3_ave = mean(c(Run3_1, Run3_2), na.rm=TRUE)
         , Run4_ave = mean(c(Run4_1, Run4_2), na.rm=TRUE))  %>%
  mutate(ave_over_same_plate = sum(!is.na(c(Run1_ave, Run2_ave, Run3_ave, Run4_ave)))==1
    , sd = ifelse(sum(!is.na(c(Run1_ave, Run2_ave, Run3_ave, Run4_ave)))==1, sd(c(Run1_1, Run1_2, Run2_1, Run2_2, Run3_1, Run3_2, Run4_1, Run4_2), na.rm=TRUE)
                     , sd(c(Run1_ave, Run2_ave, Run3_ave, Run4_ave), na.rm=TRUE))  
         , ave_log10copies = ifelse(sum(!is.na(c(Run1_ave, Run2_ave, Run3_ave, Run4_ave)))==1, mean(c(Run1_1, Run1_2, Run2_1, Run2_2, Run3_1, Run3_2, Run4_1, Run4_2), na.rm=TRUE)
                                  , mean(c(Run1_ave, Run2_ave, Run3_ave, Run4_ave), na.rm=TRUE))) %>%
  # Finally, choose the 2 replicates to compare consistency
  mutate(R1 = ifelse(ave_over_same_plate, c(Run1_1,Run1_2,Run2_1,Run2_2, Run3_1, Run3_2, Run4_1, Run4_2)[-which(is.na(c(Run1_1,Run1_2,Run2_1,Run2_2, Run3_1, Run3_2, Run4_1, Run4_2)))][1]
                   , c(Run1_ave,Run2_ave, Run3_ave, Run4_ave)[-which(is.na(c(Run1_ave,Run2_ave, Run3_ave, Run4_ave)))][1])
         ,R2 = ifelse(ave_over_same_plate, c(Run1_1,Run1_2,Run2_1,Run2_2, Run3_1, Run3_2, Run4_1, Run4_2)[-which(is.na(c(Run1_1,Run1_2,Run2_1,Run2_2, Run3_1, Run3_2, Run4_1, Run4_2)))][2]
                      , c(Run1_ave,Run2_ave, Run3_ave, Run4_ave)[-which(is.na(c(Run1_ave,Run2_ave, Run3_ave, Run4_ave)))][2])
         ,R3 = ifelse(ave_over_same_plate, c(Run1_1,Run1_2,Run2_1,Run2_2, Run3_1, Run3_2, Run4_1, Run4_2)[-which(is.na(c(Run1_1,Run1_2,Run2_1,Run2_2, Run3_1, Run3_2, Run4_1, Run4_2)))][3]
                      , c(Run1_ave,Run2_ave, Run3_ave, Run4_ave)[-which(is.na(c(Run1_ave,Run2_ave, Run3_ave, Run4_ave)))][3])
         ,R4 = ifelse(ave_over_same_plate, c(Run1_1,Run1_2,Run2_1,Run2_2, Run3_1, Run3_2, Run4_1, Run4_2)[-which(is.na(c(Run1_1,Run1_2,Run2_1,Run2_2, Run3_1, Run3_2, Run4_1, Run4_2)))][4]
                       , c(Run1_ave,Run2_ave, Run3_ave, Run4_ave)[-which(is.na(c(Run1_ave,Run2_ave, Run3_ave, Run4_ave)))][4]))

# Identify which ones are outliers
ggsave(filename = "EDA_qPCR/Replicate_consistency.pdf", width=4, height=3.5
       ,dat_variance_adj %>%
  mutate(Gr_than_2SD_variation = sd>2) %>%
  ggplot() + geom_point(aes(x=(R1), y=(R2), col=Gr_than_2SD_variation)) +
  geom_abline(aes(slope=1, intercept=0)) +
  scale_color_discrete(name = ">2 St.Dev") + xlab("Rep 1 Copy estimate") +
  ylab("Rep 2 Copy estimate"))

ggsave(filename = "EDA_qPCR/Replicate_consistency_primersplit.pdf", width=5, height=3
       ,dat_variance_adj %>%
         mutate(Gr_than_2SD_variation = sd>2) %>%
         ggplot() + geom_point(aes(x=(R1), y=(R2), col=Gr_than_2SD_variation)) +
         geom_abline(aes(slope=1, intercept=0)) +
         scale_color_discrete(name = ">2 St.Dev") + xlab("Rep 1 Copy estimate (log10)") +
         ylab("Rep 2 Copy estimate (log10)") + 
         facet_grid(.~primers))

dat_variance_adj %>%
  ggplot() + geom_histogram(aes(x=log10(sd), col=ave_over_same_plate)) +
  facet_grid(.~primers)

# Standardized variance plots
dat_variance_adj %>%
  ggplot() + geom_histogram(aes(x=log10(sd), col=ave_over_same_plate))
dat_variance_adj %>%
  ggplot() + geom_histogram(aes(x=log10(sd), col=anyedge))
dat_variance_adj %>%
  ggplot() + geom_histogram(aes(x=log10(sd), col=primers))

# Is between-run vs within-run variance any different from each other?
dat_variance_adj %>% filter(nRep == 2) %>% 
  mutate(diff = abs((R1)-(R2))) %>% select(sampleID, R1, R2, diff, ave_over_same_plate, primers) %>%
  ggplot() + # geom_boxplot(aes(x=dup_on_same_run, y=diff))
  geom_histogram(aes(x=log10(diff), col=ave_over_same_plate))+
  facet_grid(.~primers)
dat_variance_adj %>% filter(nRep == 2) %>% 
  mutate(diff = abs((R1)-(R2))) %>% select(sampleID, R1, R2, diff, ave_over_same_plate, primers) %>%
  ggplot() +  geom_boxplot(aes(x=ave_over_same_plate, y=log10(diff)))+
  # geom_histogram(aes(x=log(diff), col=dup_on_same_run))+
  facet_grid(.~primers)

### Save output ####
#*** Before saving, make sure to recalculate total SD with new, filtered data.
# Omit outliers by summarizing within-sample variability, and then omitting those that go beyond 2SD if >2 replicates
sumsq_n <- dat_variance_adj %>% 
  mutate(sumsq = sum(((c(Run1_1,Run1_2,Run2_1,Run2_2, Run3_1, Run3_2, Run4_1, Run4_2))-ave_log10copies)^2, na.rm=TRUE)
         ,nRep = sum(!is.na(c(Run1_1,Run1_2,Run2_1,Run2_2, Run3_1, Run3_2, Run4_1, Run4_2)))) %>%
  select(sumsq, nRep) %>%
  colSums()
totalSD <- sqrt(sumsq_n[1]/(sumsq_n[2]-1))
# Now, use this to calculate lower and upper variance limits; save as final.

dat_variance_final <- dat_variance_adj %>%
  mutate(totalStDev = totalSD
         , lwr_95 = ave_log10copies - 2*totalStDev
         , upr_95 = ave_log10copies + 2*totalStDev) 

dat_variance_final %>%
  ggplot() + geom_point(aes(x=R1, y=R2)) +
  geom_line(aes(x=ave_log10copies, y=lwr_95)) +
  geom_line(aes(x=ave_log10copies, y=upr_95))
9/nrow(dat_variance_final) # about a 5% "outside" rate for 95% range. Acceptable.

#### Test baseline Bd DNA in non-sporangial ones.
dir.create("downstream")
write.table(dat_variance_final, file="downstream/qPCR_results_combined_edited.txt", quote=FALSE, row.names = FALSE, sep="\t")

setwd("../../")

