#!bin/bash
library(bayestestR)
library(ggbiplot)
library(tidyverse)
##### Plotting for manuscript ####
setwd("04_analysis/")
dir.create("03_figure_generation")

######## Loading ##########
dat <- read.delim("01_data_summaries/downstream/dat_scaled_final.txt")
samps_bayes <- read.delim("02_statistical_analyses/samplesAdjusted_bayes_mixedmodel.txt")
samps_bayes_noInhib <- read.delim("02_statistical_analyses/samplesAdjusted_bayes_mixedmodel_noInhib.txt")
samps_bayes_qpcrbymicr <- read.delim("02_statistical_analyses/samplesAdjust_bayes_qpcrbymicro.txt")
samps_bayes_qpcrbymicr_noInhib <- read.delim("02_statistical_analyses/samplesAdjust_bayes_qpcrbymicro_noInhib.txt")
summary_bayes <- read.delim("02_statistical_analyses/summary_bayes_mixedmodel.txt")
summary_bayes_noInhib <- read.delim("02_statistical_analyses/summary_bayes_mixedmodel_noinhib.txt")
summary_bayes_qpcrbymicr <- read.delim("02_statistical_analyses/summary_bayes_pcrbymicro.txt")
summary_bayes_qpcrbymicr_noInhib <- read.delim("02_statistical_analyses/summary_bayes_pcrbymicro_noInhib.txt")
isolate_info <- read.delim("../00_isolate_info/all_isolate_info_combined.txt", sep =",")
original_dat_withASV <- read.delim("../01_experiment_metadata/downstream/aggregated_metadata.txt")
# dat_centres <- read.delim("01_data_summaries/downstream/dat_centers.txt")
dat_full_with_iso_summary <- read.delim("02_statistical_analyses/dat_full_with_iso_summary.txt")
## For ASV-level 
isolate_indiv_effects_res <- read.delim("02_statistical_analyses/isolate_indiv_effects_res.txt")
isolate_indiv_effects_summary<- read.delim("02_statistical_analyses/isolate_indiv_effects_summary.txt")
challenge_iso_summary<- read.delim("02_statistical_analyses/challenge_iso_summary.txt")
isolate_indiv_effects_together<- read.delim("02_statistical_analyses/isolate_indiv_effects_together.txt")

# Re-centre data
dat_c <- dat %>% 
  mutate(Richness = Richness + Centre_Richness, CV_log10 = CV_log10+ Centre_CV_log10, qPCR_bact_log10 = qPCR_bact_log10 + Centre_qPCR_bact_log10)
#### /----------Data Summaries--------/ ####
gg_predhistogram <- dat_c %>% 
  # mutate(Richness = Richness + Centre_Richness, CV_log10 = CV_log10+ Centre_CV_log10, qPCR_bact_log10 = qPCR_bact_log10 + Centre_qPCR_bact_log10) %>%
  select(SampleID, Richness, CV_log10, qPCR_bact_log10, inhibRichFrac_raw) %>%
  pivot_longer(-SampleID, names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = ifelse(Metric == "Richness", "Bacterial richness",ifelse(Metric == "CV_log10", "Biofilm thickness \n(CV log10+0.01)", ifelse(Metric=="qPCR_bact_log10","Bacterial cell density \n(log10+1)", "Fraction of ASVs\nthat were inhibitory")))) %>%
  mutate(Metric = factor(Metric, levels = c("Bacterial richness", "Bacterial cell density \n(log10+1)", "Biofilm thickness \n(CV log10+0.01)", "Fraction of ASVs\nthat were inhibitory"))) %>%
  ggplot() + geom_histogram(aes(x=Value)) +
  facet_wrap(.~Metric, scales="free") + ylab("Count")
  
ggsave(filename = "03_figure_generation/histogram_predictors.png", height=6, width=6
       ,gg_predhistogram)

gg_resphistogram <- dat_c %>% 
  select(SampleID, Bd_qPCR_log10, Bd_micro_log10) %>%
  pivot_longer(-SampleID, names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = ifelse(Metric == "Bd_qPCR_log10", "Bd gene copies \n(qPCR log10+1)","Bd microscopy counts \n(log10 +1 per mm2)")) %>%
  ggplot() + geom_histogram(aes(x=Value)) +
  facet_wrap(.~Metric, scales="free") + ylab("Count")

ggsave(filename = "03_figure_generation/histogram_response.png", height=3, width=5
       ,gg_resphistogram)

#### Raw data plots ####
dat_c %>%
  # mutate(Inhibitory_designation = ifelse(Inhibitory==1, "Inhibitory", "Non-inhibitory")) %>%
  ggplot(aes(x=CV_log10, y=Bd_micro_log10, col=inhibRichFrac_raw)) + geom_point()

##### Bayes ######
# Random effects
gg_bayes_randef <- summary_bayes%>% filter(Predictor == "Set") %>%
  mutate(response = ifelse(response=="Bd_micro_log10", "Bd density (log10 +1, /mm2)", "Bd qPCR copies (log10)")) %>%
  ggplot() + geom_segment(aes(x=lower95twoside, xend=upper95twoside, y=GroupOrType, yend=GroupOrType)) +
  geom_point(aes(x=Median, y=GroupOrType)) +
  geom_vline(aes(xintercept=0),col="black", lty=3) +
  facet_grid(.~response, scales = "free") +
  xlab("Coefficient estimate") +ylab("Starting date of set") 
ggsave(filename = "03_figure_generation/bayes_samp_raneff.png", width = 9, height=5
       ,gg_bayes_randef)

# Look at individual histograms to check for good sampling 

gg_bayes_fixef_micro <- summary_bayes %>% filter(Predictor != "Set", response == "Bd_micro_log10") %>%
  mutate(GroupOrType = ifelse(GroupOrType == "Non",  "Non-inhibitory \nbiofilms", ifelse(GroupOrType == "Inhib", "Inhibitory \nbiofilms", "Diff between \n Inhib/Non-inhib"))) %>%
  mutate(GroupOrType = factor(GroupOrType, levels = rev(c("Non-inhibitory \nbiofilms","Inhibitory \nbiofilms","Diff between \n Inhib/Non-inhib")))) %>%
  mutate(Predictor = ifelse(Predictor == "CV_log10", "Biofilm thickness \n(log10 CV+0.01)",ifelse(Predictor == "qPCR_bact_log10","Bacterial density \n(log10 qPCR)", ifelse(Predictor == "Richness", "Bacterial richness", Predictor)) )) %>%
  mutate(Predictor = factor(Predictor, levels = c("Intercept", "Bacterial richness", "Bacterial density \n(log10 qPCR)", "Biofilm thickness \n(log10 CV+0.01)"))) %>%
  mutate(Significant = PD<0.025) %>%
  filter(GroupOrType != "Diff between \n Inhib/Non-inhib") %>%
  ggplot() + geom_segment(aes(x=lower95twoside, xend=upper95twoside, y=GroupOrType, yend=GroupOrType, col=Significant)) +
  geom_point(aes(x=Median, y=GroupOrType, col=Significant)) +
  geom_vline(aes(xintercept=0),col="black", lty=3) +
  facet_grid(.~Predictor, scales = "free") +
  scale_color_manual(values=c("grey","darkred")) +ylab("Biofilm type") +xlab("Effect on Bd microscopy counts (coefficient estimate)") +
  labs(col="Does 95% CI\n include zero?")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

ggsave(filename = "03_figure_generation/bayes_samp_fixef_micro.png", width = 7, height=3
       ,gg_bayes_fixef_micro)

gg_bayes_fixef_qpcr <- summary_bayes %>% filter(Predictor != "Set", response == "Bd_qPCR_log10") %>%
  mutate(GroupOrType = ifelse(GroupOrType == "Non",  "Non-inhibitory \nbiofilms", ifelse(GroupOrType == "Inhib", "Inhibitory \nbiofilms", "Diff between \n Inhib/Non-inhib"))) %>%
  mutate(GroupOrType = factor(GroupOrType, levels = rev(c("Non-inhibitory \nbiofilms","Inhibitory \nbiofilms","Diff between \n Inhib/Non-inhib")))) %>%
  mutate(Predictor = ifelse(Predictor == "CV_log10", "Biofilm thickness \n(log10 CV+0.01)",ifelse(Predictor == "qPCR_bact_log10","Bacterial density \n(log10 qPCR)", ifelse(Predictor == "Richness", "Bacterial richness", Predictor)) )) %>%
  mutate(Predictor = factor(Predictor, levels = c("Intercept", "Bacterial richness", "Bacterial density \n(log10 qPCR)", "Biofilm thickness \n(log10 CV+0.01)"))) %>%
  mutate(Significant = PD<0.025) %>%
  filter(GroupOrType != "Diff between \n Inhib/Non-inhib") %>%
  ggplot() + geom_segment(aes(x=lower95twoside, xend=upper95twoside, y=GroupOrType, yend=GroupOrType, col=Significant)) +
  geom_point(aes(x=Median, y=GroupOrType, col=Significant)) +
  geom_vline(aes(xintercept=0),col="black", lty=3) +
  facet_grid(.~Predictor, scales = "free") +
  scale_color_manual(values=c("grey","darkred")) +ylab("Biofilm type") +xlab("Effect on Bd qPCR copies (coefficient estimate)") +
  labs(col="Does 95% CI\n include zero?")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))


ggsave(filename = "03_figure_generation/bayes_samp_fixef_qpcr.png", width = 7, height=3
       ,gg_bayes_fixef_qpcr)



gg_bayes_fixef_qpcrbymicro <- summary_bayes_qpcrbymicr %>% filter(Predictor != "Set") %>%
  mutate(GroupOrType = ifelse(GroupOrType == "Non",  "Non-inhibitory \nbiofilms", ifelse(GroupOrType == "Inhib", "Inhibitory \nbiofilms", "Diff between \n Inhib/Non-inhib"))) %>%
  mutate(GroupOrType = factor(GroupOrType, levels = rev(c("Non-inhibitory \nbiofilms","Inhibitory \nbiofilms","Diff between \n Inhib/Non-inhib")))) %>%
  mutate(Predictor = ifelse(Predictor == "CV_log10", "Biofilm thickness \n(log10 CV+0.01)",ifelse(Predictor == "qPCR_bact_log10","Bacterial density \n(log10 qPCR)", ifelse(Predictor == "Richness", "Bacterial richness", ifelse(Predictor == "Bd_micro_log10","Bd zoosporangia \n(log10+1 /mm2)", Predictor))) )) %>%
  mutate(Predictor = factor(Predictor, levels = c("Intercept", "Bd zoosporangia \n(log10+1 /mm2)","Bacterial richness", "Bacterial density \n(log10 qPCR)", "Biofilm thickness \n(log10 CV+0.01)"))) %>%
  mutate(Significant = PD<0.025) %>%
  filter(GroupOrType != "Diff between \n Inhib/Non-inhib", Predictor != "Bd zoosporangia \n(log10+1 /mm2)") %>%
  ggplot() + geom_segment(aes(x=lower95twoside, xend=upper95twoside, y=GroupOrType, yend=GroupOrType, col=Significant)) +
  geom_point(aes(x=Median, y=GroupOrType, col=Significant)) +
  geom_vline(aes(xintercept=0),col="black", lty=3) +
  facet_grid(.~Predictor, scales = "free") +
  scale_color_manual(values=c("grey","darkred")) +xlab("Effect on residual Bd qPCR copies\nafter accounting for Bd microscopy counts") +ylab("Biofilm type") +
  labs(col="Does 95% CI\n include zero?")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

ggsave(filename = "03_figure_generation/bayes_samp_fixef_qpcrbymicro.png", width = 7, height=3
       ,gg_bayes_fixef_qpcrbymicro)

##### PCA #####
dat_prcomp <- dat %>% select(Bd_micro_log10, Bd_qPCR_log10, CV_log10, qPCR_bact_log10, Richness, inhibRichFrac_raw) %>%
  dplyr::rename("Bd microscopy" = Bd_micro_log10, "Bd qPCR"= Bd_qPCR_log10, "Biofilm thickness (CV)"= CV_log10, "Bact density (qPCR)"=qPCR_bact_log10, "Bact Richness"=Richness, "Fraction of ASVs\ninhibitory"=inhibRichFrac_raw)
  # mutate(Predictor = ifelse(Predictor == "CV_log10", "Biofilm thickness \n(log10 CV+0.01)",ifelse(Predictor == "qPCR_bact_log10","Bacterial density \n(log10 qPCR)", ifelse(Predictor == "Richness", "Bacterial richness", ifelse(Predictor == "Bd_micro_log10","Bd zoosporangia \n(log10+1 /mm2)", Predictor))) )) 
prcomp_all <- prcomp(dat_prcomp, scale. = TRUE)
# Flip PCA1
# prcomp_all$rotation <- prcomp_all$rotation*-1
ggsave("03_figure_generation/PCA.png", width=6.5, height=4,
       ggbiplot(prcomp_all, cex = 0.5, groups = dat_prcomp$`Bd microscopy`, choices = c(1,2), varname.adjust = 1.1) +
         xlim(-3,3.5) + ylim(-3,3) +
         labs(col = "Bd microscopy \n(log10 + 1 /mm2)")
)
# NO inhibitory; for presentations
dat_prcomp_noinhib <- dat %>% select(Bd_micro_log10, Bd_qPCR_log10, CV_log10, qPCR_bact_log10, Richness) %>%
  dplyr::rename("Bd microscopy" = Bd_micro_log10, "Bd qPCR"= Bd_qPCR_log10, "Biofilm thickness(CV)"= CV_log10, "Bact density (qPCR)"=qPCR_bact_log10, "Bact Richness"=Richness)
# mutate(Predictor = ifelse(Predictor == "CV_log10", "Biofilm thickness \n(log10 CV+0.01)",ifelse(Predictor == "qPCR_bact_log10","Bacterial density \n(log10 qPCR)", ifelse(Predictor == "Richness", "Bacterial richness", ifelse(Predictor == "Bd_micro_log10","Bd zoosporangia \n(log10+1 /mm2)", Predictor))) )) 
prcomp_noinhib <- prcomp(dat_prcomp_noinhib, scale. = TRUE)
prcomp_noinhib$rotation <- -1*prcomp_noinhib$rotation
prcomp_noinhib$x <- -1*prcomp_noinhib$x
# Flip PCA1
# prcomp_all$rotation <- prcomp_all$rotation*-1
ggsave("03_figure_generation/PCA_noinhib.png", width=6.5, height=4,
       ggbiplot(prcomp_noinhib, cex = 0.5, groups = dat_prcomp$`Bd microscopy`, choices = c(1,2), varname.adjust = 1.1) +
         ylim(-3,3) + xlim(-3.5,3) + 
         labs(col = "Bd microscopy \n(log10 + 1 /mm2)") 
)

#### /--------- Predicting CV using community -------/####

### Full comunity ####
# dat_full_with_iso_summary <- dat_iso_summary %>% full_join(dat_full) %>%
#   mutate(CV = 10^CV_log10-0.01) 

dat_CV_assessment <- dat_full_with_iso_summary %>%
  filter(maxCV!=0) %>%
  rowwise() %>% mutate(propDom = max(across(starts_with("ISO")))/sum(across(starts_with("ISO")))) %>% ungroup() %>%
  mutate(propDom = ifelse(Richness==1, NA, propDom)) %>%
  # select(CV, sumCV, maxCV, adjCV, propDom,domIsoCV, Richness) %>% 
  select(CV, maxCV, adjCV, sumCV, propDom, Richness) %>%
  pivot_longer(-c(propDom, Richness, CV), values_to = "CV_predict", names_to = "CV_type")  %>%
  mutate(CV_type = ifelse(CV_type=="adjCV", "PROPORTIONAL:\nWeighted average of crystal violet\n across members",
                          ifelse(CV_type=="maxCV", "COMPLEMENTARY:\nMaximum crystal violet intensity\n out of all members", 
                                 ifelse(CV_type=="sumCV", "ADDITIVE:\nSummed crystal violet intensity\nacross all members", 
                                        ifelse(CV_type=="adjCV_raw", "Read-weighted sum of crystal violet\nacross members",NA))))) %>%
  mutate(Rich1 = ifelse(Richness==1, "1", "Greater than 1")) %>%
  mutate(CV_type = factor(CV_type, levels=c("ADDITIVE:\nSummed crystal violet intensity\nacross all members", "COMPLEMENTARY:\nMaximum crystal violet intensity\n out of all members","PROPORTIONAL:\nWeighted average of crystal violet\n across members","Read-weighted sum of crystal violet\nacross members")))
gg_CVpredict <- dat_CV_assessment %>% ggplot() + 
  geom_point(aes(x=log10(CV_predict+0.01), y=log10(CV+0.01), pch=factor(Rich1), fill=propDom), cex=3, col="black") +
  geom_smooth(data=dat_CV_assessment %>% filter(Richness==1), aes(x=log10(CV_predict+0.01), y=log10(CV+0.01)), col="black", method="lm", fullrange=TRUE) +
  # geom_point(aes(x=log10(sumCV+0.01), y=log10(CV+0.01), col=factor(RichLevel)), pch=21) +
  # geom_point(aes(x=log10(maxCV+0.01), y=log10(CV+0.01), col=factor(RichLevel)), pch=8) +
  xlab("Predicted crystal violet intensities\nbased on community composition and 96-well assay")+
  ylab("Observed Crystal violet intensity\nof microcosm community") +
  labs(fill="Proportion of community\nthat is isolate with\nhighest CV", pch="Richness")  +
  facet_grid(.~CV_type, scale="free") +
  scale_fill_gradient(low="darkgrey", high="purple", na.value="black") +
  scale_shape_manual(values=c(19, 24))+theme_bw()
gg_CVpredict
ggsave(filename = "03_figure_generation/predict_CV_with_CV.png", width = 9, height=4
       ,gg_CVpredict)
### Singles only ####
load("02_statistical_analyses/dat_singlesOnly.RData")
gg_singlestrends <- dat_singlesOnly %>% 
  mutate(Predictor = ifelse(Predictor == "CV_log10", "Biofilm thickness\n(log10 CV + 0.01)", "Biofilm cell density\n(log10 qPCR copies+1)")) %>%
  mutate(Bd = ifelse(Bd == "Bd_micro_log10", "Bd microscopy counts", "Bd qPCR copies")) %>%
  ggplot(aes(x= predictor_value, y=response_value, col=Isolate_designation)) + geom_point() +
  # geom_smooth(method="lm", se = FALSE) +
  geom_abline(aes(slope=m, intercept=b, col=Isolate_designation, lty=PD)) +
  facet_grid(Bd ~ Predictor, scales="free", switch="both") +xlab("Predictor value") + ylab("Bd metric") +
  theme_bw()
gg_singlestrends
ggsave(filename = "03_figure_generation/singles_trends_biofilm.png", height=4, width=6,
       gg_singlestrends)

#### /------- Predicting Bd settlement using inhibitory -------/####
gg_InhibPredict <- dat_full_with_iso_summary %>% filter(Richness>0, inhibRich>0) %>%
  select(Bd_micro_log10, Bd_qPCR_log10, maxinhib, suminhib, adjinhib) %>%
  pivot_longer(-c(Bd_micro_log10, Bd_qPCR_log10), values_to="Value", names_to="Metric") %>%
  pivot_longer(-c(Value, Metric), values_to="Bd", names_to="Response") %>%
  mutate(Metric = ifelse(Metric=="maxinhib", "COMPLEMENTARY:\nMaximum inhibitory zone of\nmost effective inhibiting isolate", 
                         ifelse(Metric == "suminhib", "ADDITIVE:\nSum of inhibitory zones of\nall inhibiting isolates"
                                , ifelse(Metric == "adjinhib", "AVERAGE:\nAverage zone of inhibition\nacross inhibitory isolates", 
                                         ifelse(Metric == "inhibRich", "RICHNESS:\nNumber of inhibitory isolates", NA))))) %>%
  mutate(Metric = factor(Metric, levels=c("ADDITIVE:\nSum of inhibitory zones of\nall inhibiting isolates",
                                          "COMPLEMENTARY:\nMaximum inhibitory zone of\nmost effective inhibiting isolate", 
                                          "AVERAGE:\nAverage zone of inhibition\nacross inhibitory isolates", 
                                          "RICHNESS:\nNumber of inhibitory isolates"))) %>%
  mutate(Response = ifelse(Response=="Bd_micro_log10", "Bd microscopy counts (log10)", "Bd qPCR copies (log10)"))  %>%
  ggplot() + geom_jitter(aes(x=Value, y=Bd), cex=3, width=0.3, height=0) +
  geom_smooth(aes(x=Value, y=Bd), method="lm") +
  facet_grid(Response~Metric, scales = "free", switch = "both") + ylab("Bd metric")+
  theme_bw()
gg_InhibPredict
ggsave(filename = "03_figure_generation/predict_Bd_with_inhib.png", width = 9, height=5
       ,gg_InhibPredict)


#### /----Correlation plots----/ ####

#### Rich vs biomass ####
gg_rich_biomass <- dat_c %>%
  ggplot() + 
  geom_jitter(aes(x=Richness, y= qPCR_bact_log10, col = CV_log10), cex = 2, width=0.2, height=0) +
  scale_color_continuous(low = "darkgrey", high="purple") +
  theme_bw() +ylab("Bacterial cell density \n(log10 qPCR)") + xlab("Bacterial community richness") +
  labs(col = "Biofilm thickness\n(log10 CV)")
ggsave(filename = "03_figure_generation/corr_rich_biomass.png", width = 5, height=3
       ,gg_rich_biomass)

#### CV vs qPCR ####
gg_cv_qpcr <- dat_c %>%
  ggplot() + geom_point(aes(x=CV_log10, y=qPCR_bact_log10), cex =2) +
  theme_bw() + ylab("Bacterial cell density\n (log10 qPCR)") + xlab("Biofilm thickness\n(log10 CV)")
ggsave(filename = "03_figure_generation/corr_cv_qpcr.png", width = 4, height=3
       ,gg_cv_qpcr)

#### CV vs qPCR ####
# gg_cv_inhibitory <- dat_c %>%
#   mutate(Inhibitory = ifelse(Inhibitory == 1, "Inhibitory", "Non-inhibitory")) %>%
#   ggplot(aes(x=Inhibitory, y=CV_log10)) + geom_violin() + geom_jitter(cex =2, width = 0.3, height=0) +
#   theme_bw() + xlab("Inhibitory designation for biofilm community") + ylab("Biofilm thickness\n(log10 CV)")
# ggsave(filename = "03_figure_generation/corr_cv_inhibitory.png", width = 4, height=3
#        ,gg_cv_inhibitory)

gg_cv_inhibitory <- dat_c %>%
  # mutate(Inhibitory = ifelse(Inhibitory == 1, "Inhibitory", "Non-inhibitory")) %>%
  ggplot(aes(x=inhibRichFrac_raw, y=CV_log10)) +geom_point() +
  geom_smooth(method="lm") +
  theme_bw() + xlab("Fraction of ASVs that were inhibitory") + ylab("Biofilm thickness\n(log10 CV)")
ggsave(filename = "03_figure_generation/corr_cv_inhibitory.png", width = 4, height=3
       ,gg_cv_inhibitory)

#### Bd qPCR vs Bd micro ####
gg_bd_withinhib <- dat_c %>%
  # mutate(Inhibitory = ifelse(Inhibitory == 1, "Inhibitory", "Non-Inhibitory")) %>%
  ggplot() + geom_point(aes(x=Bd_micro_log10, y=Bd_qPCR_log10, col = inhibRichFrac_raw)) +
  geom_smooth(aes(x=Bd_micro_log10, y=Bd_qPCR_log10), col="black", method="lm") +
  scale_color_continuous(low="darkgrey", high="purple") +
  xlab("Bd zoosporangia density \n(log10+1 /mm2)") +ylab("Bd copy abundance \n(log10+1 qPCR)") +
  labs(col = "Fraction of ASVs\nthat were inhibitory") +
  theme_bw()
ggsave(filename = "03_figure_generation/corr_bd_withinhib.png",width=4.5, height=3
       ,gg_bd_withinhib)

#### CV not monotomic ####
gg_qPCRmonotonic <- dat_c %>%
  ggplot() + 
  geom_jitter(aes(x=Richness, y= qPCR_bact_log10, col = CV_log10), cex = 2, width=0.2, height=0) +
  geom_line(aes(x=Richness, y= qPCR_bact_log10, group=Rep), col="lightgrey") +
  scale_color_continuous(low = "darkgrey", high="purple") +
  theme_bw() +ylab("Biofilm thickness\n(log10 CV)") + xlab("Bacterial community richness") +
  labs(col = "Bacterial cell density \n(log10 qPCR)")
gg_qPCRmonotonic

gg_CVmonotonic <- dat_c %>%
  ggplot() + 
  geom_point(aes(x=Richness, y= CV_log10, col = qPCR_bact_log10), cex = 2) +
  geom_line(aes(x=Richness, y= CV_log10, group=Rep)) +
  scale_color_continuous(low = "darkgrey", high="blue") +
  theme_bw() +ylab("Biofilm thickness\n(log10 CV)") + xlab("Bacterial community richness") +
  labs(col = "Bacterial cell density \n(log10 qPCR)")
gg_CVmonotonic
ggsave(filename = "03_figure_generation/corr_CV_notmonotomic.png", width = 5, height=3
       ,gg_CVmonotonic)

# #### /------------ Plotting marginal effects ------------/ ####
#### With inhibitory marginals ####

B_micro_mat <- samps_bayes %>% filter(response == "Bd_micro_log10", GroupOrType !="Inhib") %>%
  select(Coefficient, Estimate) %>% group_by(Coefficient) %>% mutate(iter = seq(1, n())) %>%
  pivot_wider(names_from="Coefficient", values_from = "Estimate") %>% select(-iter) %>% t()

B_qpcr_mat <- samps_bayes %>% filter(response == "Bd_qPCR_log10", GroupOrType !="Inhib") %>%
  select(Coefficient, Estimate) %>% group_by(Coefficient) %>% mutate(iter = seq(1, n())) %>%
  pivot_wider(names_from="Coefficient", values_from = "Estimate") %>% select(-iter) %>% t()

B_qpcrbymicro_mat <- samps_bayes_qpcrbymicr %>% filter(response == "qPCRbyMicro", GroupOrType !="Inhib") %>%
  select(Coefficient, Estimate) %>% group_by(Coefficient) %>% mutate(iter = seq(1, n())) %>%
  pivot_wider(names_from="Coefficient", values_from = "Estimate") %>% select(-iter) %>% t()

### RECENTRE-ING
dat_centres <- dat %>%
  select(starts_with("Centre")) %>% distinct() %>% select(-Centre_qPCR_bact_log10_inhib) %>%
  rename_at(vars(everything()), function(x) gsub("Centre_","",x)) %>%
  t() %>% as.data.frame() %>% rownames_to_column(var="Predictor") %>% rename("Centre"=V1)

######## TRY NON1-0 INHIBITORY DES

X_mat_micro  <- dat %>%
  mutate(D = 1) %>% pivot_wider(names_from = DateStart, values_from = D, values_fill=0) %>%
  mutate(InterceptXNon = 1, RichnessXNon = Richness, CV_log10XNon = CV_log10, qPCR_bact_log10XNon = qPCR_bact_log10
         , InterceptXInhibEff = inhibRichFrac_raw, RichnessXInhibEff = Richness*inhibRichFrac_raw, CV_log10XInhibEff = CV_log10*inhibRichFrac_raw, qPCR_bact_log10XInhibEff = qPCR_bact_log10*inhibRichFrac_raw) %>%
  select(one_of(rownames(B_micro_mat)))
X_mat_qpcr  <- dat %>%
  mutate(D = 1) %>% pivot_wider(names_from = DateStart, values_from = D, values_fill=0) %>%
  mutate(InterceptXNon = 1, RichnessXNon = Richness, CV_log10XNon = CV_log10, qPCR_bact_log10XNon = qPCR_bact_log10
         , InterceptXInhibEff = inhibRichFrac_raw, RichnessXInhibEff = Richness*inhibRichFrac_raw, CV_log10XInhibEff = CV_log10*inhibRichFrac_raw, qPCR_bact_log10XInhibEff = qPCR_bact_log10*inhibRichFrac_raw) %>%
  select(one_of(rownames(B_qpcr_mat)))
X_mat_qpcrbymicro  <- dat %>%
  mutate(D = 1) %>% pivot_wider(names_from = DateStart, values_from = D, values_fill=0) %>%
  mutate(InterceptXNon = 1, RichnessXNon = Richness, CV_log10XNon = CV_log10, qPCR_bact_log10XNon = qPCR_bact_log10
         , InterceptXInhibEff = inhibRichFrac_raw, RichnessXInhibEff = Richness*inhibRichFrac_raw, CV_log10XInhibEff = CV_log10*inhibRichFrac_raw, qPCR_bact_log10XInhibEff = qPCR_bact_log10*inhibRichFrac_raw
         , Bd_micro_log10XNon = Bd_micro_log10) %>%
  select(one_of(rownames(B_qpcrbymicro_mat)))

# dat %>% select(starts_with("Centre")) %>% distinct() %>% mutate(Intercept = 0, InhibEff = 0)

y_mat_micro <- dat %>% select(Bd_micro_log10)
y_mat_qpcr <- dat %>% select(Bd_qPCR_log10)
y_mat_qpcrbymicro <- dat %>% select(Bd_qPCR_log10) %>% rename("qPCRbyMicro"=Bd_qPCR_log10)
############ ADD IN INTERCEPT???

allMarginalForPlotting <- data.frame()
allCredForPlotting <- data.frame()
allSlopesForPlotting <- data.frame()
#### Intercepts NOT included in residual bd
for ( x in c("Richness","CV_log10","qPCR_bact_log10")) {
  # x="Richness"
  # x="qPCR_bact_log10"
  for ( y in c("Bd_micro_log10", "Bd_qPCR_log10", "qPCRbyMicro")) {
    # y="Bd_micro_log10"
    # Calling specific matrices; setting y labels
    if ( y == "Bd_micro_log10") {
      # B_bdtemp <- B_micro
      y_mat_temp <- y_mat_micro
      y_lab <- "Residual variation \nin Bd (counts)"
      X_mat <- X_mat_micro
      B_bdtemp <- B_micro_mat[match(colnames(X_mat), rownames(B_micro_mat)),]
      
    } else if (y=="Bd_qPCR_log10") {
      # B_bdtemp <- B_qpcr
      y_mat_temp <- y_mat_qpcr
      y_lab <- "Residual variation \nin Bd (qPCR)"
      X_mat <- X_mat_qpcr
      B_bdtemp <- B_qpcr_mat[match(colnames(X_mat), rownames(B_qpcr_mat)),]
      
    } else {
      y_mat_temp <- y_mat_qpcrbymicro
      y_lab <- "Residual variation \nin Bd (residual qPCR)"
      X_mat <- X_mat_qpcrbymicro
      B_bdtemp <- B_qpcrbymicro_mat[match(colnames(X_mat), rownames(B_qpcrbymicro_mat)),]
    }
    # aesthetics for plotting x axes
    if (x == "qPCR_bact_log10") {
      x_lab <- "Bacterial density \n(log10 qPCR)"
    } else if (x == "CV_log10") {
      x_lab <- "Biofilm thickness \n(log10 CV)"
    } else if (x == "Richness") {
      x_lab <- "Bacterial community \n Richness"
    }
    # Which to exclude to plot for marginal
    pos <- c(grep("Intercept", colnames(X_mat)), grep(x, colnames(X_mat)))
    # Get the marginal predictor values we're actually looking for
    X_marg <- as.matrix(X_mat)[,c(pos)]
    B_marg <- as.matrix(B_bdtemp[c(pos),])
    # Get the marginal effects by getting pred y
    X_temp <- as.matrix(X_mat)[,-c(pos)]
    B_temp <- as.matrix(B_bdtemp[-c(pos),])
    # colnames(X_temp) == rownames(B_temp)
    y_marg <- X_temp %*% B_temp
    y_res <- y_mat_temp - apply(y_marg,1,median)
    # Get slope of the marginal effects
    samps_non_slope <- B_marg[grep(paste0(x,"XNon"), rownames(B_marg)),]
    samps_inhib_slope <- colSums(B_marg[grep(paste0(x), rownames(B_marg)),])
    # Get intercept of overall
    samps_non_inter <-  B_marg[grep(paste0("InterceptXNon"), rownames(B_marg)),] - samps_non_slope*dat_centres[which(dat_centres$Predictor==x), "Centre"]
    samps_inhib_inter <- colSums(B_marg[grep(paste0("Intercept"), rownames(B_marg)),])- samps_inhib_slope*dat_centres[which(dat_centres$Predictor==x), "Centre"]
    
    # Adjust X_marg to be "real" values
    X_marg_adj <- X_marg + dat[,grep(paste0("Centre_",x,"$"), colnames(dat))]
    X_all_adjmarg <- cbind(X_temp, X_marg_adj, Inhibitory_fraction = X_mat$InterceptXInhibEff) %>%
      as_tibble() %>%
      mutate(Inhibitory_designation = ifelse(Inhibitory_fraction==0, 0, 1))
    # mutate(Inhibitory_designation = ifelse(Inhibitory_fraction==0, "Non-inhibitory", "Inhibitory"))
    
    ### Summary sheets ###
    ## Make datasheet
    tempNewDat <- data.frame(y_res, X_all_adjmarg ) %>%
      rename(X = paste0(x,"XNon"), Y = paste0(y))  %>%
      select(X, Y, Inhibitory_fraction,Inhibitory_designation) %>% mutate(Response = paste0(y), Predictor = paste0(x))
    ## Make slopes and intercept
    tempmb <- t(apply(B_marg,1,median)) %>% as.data.frame() %>%
      mutate(InterceptXInhib = InterceptXNon + InterceptXInhibEff, temp = get(paste0(x,"XNon")) + get(paste0(x,"XInhibEff"))) %>% #rename(paste0(x,"XInhib") = temp) %>%
      select(InterceptXNon, InterceptXInhib, paste0(x, "XNon"), temp) %>%
      pivot_longer(everything() , names_to = "TEMP", values_to = "Estimate") %>%
      mutate(TEMP = ifelse(TEMP == "temp", paste0(x,"XInhib"), TEMP)) %>%
      separate(TEMP, into = c("Predictor","Inhibitory_designation"), sep = "X") %>% 
      # mutate(Inhibitory_designation = ifelse(Inhibitory_designation == "Inhib", "Inhibitory", "Non-inhibitory")) %>%
      mutate(Inhibitory_designation = ifelse(Inhibitory_designation == "Inhib", 1, 0)) %>%
      rownames_to_column(var="row") %>%
      pivot_wider(-c("row"), names_from = Predictor, values_from = Estimate) %>%
      rename(slope = paste0(x)) %>%
      mutate(Response = paste0(y), Predictor = paste0(x)) %>%
      left_join(dat_centres) %>% mutate(intercept = Intercept-slope*Centre) %>% select(-Centre, -Intercept)
    
    ### Simulation ###
    xsim <- seq(min(X_marg_adj[,paste0(x,"XNon")]), max(X_marg_adj[,paste0(x,"XNon")]), length.out = 100)
    # Non-inhibitory
    non_sim_values <- matrix(sapply(1:length(samps_non_slope), function(x) samps_non_slope[x]*xsim + samps_non_inter[x]),ncol=length(samps_non_slope))
    CredInterNon <- data.frame(lower95 = apply(non_sim_values, 1, function(x) ci(x, method="HDI")$CI_low), upper95 = apply(non_sim_values, 1, function(x) ci(x, method="HDI")$CI_high),xsim=xsim) %>%
      mutate(InterceptXInhibEff = 0)
    # mutate(InterceptXInhibEff = "Non-inhibitory")
    # Inhibitory
    inhib_sim_values <- matrix(sapply(1:length(samps_inhib_slope), function(x) samps_inhib_slope[x]*xsim + samps_inhib_inter[x]),ncol=length(samps_inhib_slope))
    CredInterInhib <- data.frame(lower95 = apply(inhib_sim_values, 1, function(x) ci(x, method="HDI")$CI_low), upper95 = apply(inhib_sim_values, 1, function(x) ci(x, method="HDI")$CI_high),xsim=xsim) %>%
      mutate(InterceptXInhibEff = 1)
    # mutate(InterceptXInhibEff = "Inhibitory")
    CredInterAll <- rbind(CredInterInhib, CredInterNon)
    ## Make Cred datasheet
    tempNewCred <- CredInterAll %>%
      mutate(Response = paste0(y), Predictor = paste0(x))
    
    # Look at marginal effect
    #data.frame(y_res, X_all_adjmarg ) %>%
    #   mutate(InterceptXInhibEff = ifelse(InterceptXInhibEff==1, "Inhibitory","Non-inhibitory"))
    
    
    gg_temp <-   tempNewDat %>% #ggplot(aes(x=get(paste0(x,"XNon")), y=get(y), col=factor(InterceptXInhibEff))) +
      mutate(Inhibitory_designation = as.numeric(Inhibitory_designation))  %>%
      ggplot(aes(x=X, y=Y)) +
      geom_point(aes(col=Inhibitory_fraction)) +
      scale_color_gradient(low="grey", high="salmon") +
      geom_ribbon(data = tempNewCred, aes(x=xsim, ymin=lower95, ymax=upper95, fill=factor(InterceptXInhibEff)), alpha=0.1, inherit.aes = FALSE, show.legend = FALSE) +
      scale_fill_manual(values=c("grey","salmon")) +
      geom_abline(data=tempmb, aes(intercept=intercept, slope = slope, col = Inhibitory_designation)) +
      labs(col="Biofilm type") +xlab(x_lab) +ylab(y_lab) 
    # tempNewDat %>% #ggplot(aes(x=get(paste0(x,"XNon")), y=get(y), col=factor(InterceptXInhibEff))) +
    #   ggplot(aes(x=X, y=Y, col=Inhibitory_designation)) +
    #   geom_point() +
    #   geom_ribbon(data = tempNewCred, aes(x=xsim, ymin=lower95, ymax=upper95, fill=InterceptXInhibEff), alpha=0.1, inherit.aes = FALSE, show.legend = FALSE) +
    #   geom_abline(data=tempmb, aes(intercept=intercept, slope = slope, col = Inhibitory_designation)) +
    #   labs(col="Biofilm type") +xlab(x_lab) +ylab(y_lab)
    # gg_temp
    ggsave(filename = paste0("03_figure_generation/marginal_", y,"_",x,".png"), width=4, height=3
           ,gg_temp)
    # Add onto data frame
    allMarginalForPlotting <- rbind(allMarginalForPlotting, tempNewDat)
    allCredForPlotting <- rbind(allCredForPlotting, tempNewCred)
    allSlopesForPlotting <- rbind(allSlopesForPlotting, tempmb)
  }
}



# 
# X_mat_micro  <- dat %>%
#   mutate(D = 1) %>% pivot_wider(names_from = DateStart, values_from = D, values_fill=0) %>%
#   mutate(InterceptXNon = 1, RichnessXNon = Richness, CV_log10XNon = CV_log10, qPCR_bact_log10XNon = qPCR_bact_log10
#          , InterceptXInhibEff = Inhibitory, RichnessXInhibEff = Richness*Inhibitory, CV_log10XInhibEff = CV_log10*Inhibitory, qPCR_bact_log10XInhibEff = qPCR_bact_log10*Inhibitory) %>%
#   select(one_of(rownames(B_micro_mat)))
# X_mat_qpcr  <- dat %>%
#   mutate(D = 1) %>% pivot_wider(names_from = DateStart, values_from = D, values_fill=0) %>%
#   mutate(InterceptXNon = 1, RichnessXNon = Richness, CV_log10XNon = CV_log10, qPCR_bact_log10XNon = qPCR_bact_log10
#          , InterceptXInhibEff = Inhibitory, RichnessXInhibEff = Richness*Inhibitory, CV_log10XInhibEff = CV_log10*Inhibitory, qPCR_bact_log10XInhibEff = qPCR_bact_log10*Inhibitory) %>%
#   select(one_of(rownames(B_qpcr_mat)))
# X_mat_qpcrbymicro  <- dat %>%
#   mutate(D = 1) %>% pivot_wider(names_from = DateStart, values_from = D, values_fill=0) %>%
#   mutate(InterceptXNon = 1, RichnessXNon = Richness, CV_log10XNon = CV_log10, qPCR_bact_log10XNon = qPCR_bact_log10
#          , InterceptXInhibEff = Inhibitory, RichnessXInhibEff = Richness*Inhibitory, CV_log10XInhibEff = CV_log10*Inhibitory, qPCR_bact_log10XInhibEff = qPCR_bact_log10*Inhibitory
#          , Bd_micro_log10XNon = Bd_micro_log10) %>%
#   select(one_of(rownames(B_qpcrbymicro_mat)))
# 
# # dat %>% select(starts_with("Centre")) %>% distinct() %>% mutate(Intercept = 0, InhibEff = 0)
# 
# y_mat_micro <- dat %>% select(Bd_micro_log10)
# y_mat_qpcr <- dat %>% select(Bd_qPCR_log10)
# y_mat_qpcrbymicro <- dat %>% select(Bd_qPCR_log10) %>% rename("qPCRbyMicro"=Bd_qPCR_log10)
# ############ ADD IN INTERCEPT???
# 
# allMarginalForPlotting <- data.frame()
# allCredForPlotting <- data.frame()
# allSlopesForPlotting <- data.frame()
# #### Intercepts NOT included in residual bd
# for ( x in c("Richness","CV_log10","qPCR_bact_log10")) {
#   # x="Richness"
#   # x="qPCR_bact_log10"
#   for ( y in c("Bd_micro_log10", "Bd_qPCR_log10", "qPCRbyMicro")) {
#     # y="Bd_micro_log10"
#     # Calling specific matrices; setting y labels
#     if ( y == "Bd_micro_log10") {
#       # B_bdtemp <- B_micro
#       y_mat_temp <- y_mat_micro
#       y_lab <- "Residual variation \nin Bd (counts)"
#       X_mat <- X_mat_micro
#       B_bdtemp <- B_micro_mat[match(colnames(X_mat), rownames(B_micro_mat)),]
# 
#     } else if (y=="Bd_qPCR_log10") {
#       # B_bdtemp <- B_qpcr
#       y_mat_temp <- y_mat_qpcr
#       y_lab <- "Residual variation \nin Bd (qPCR)"
#       X_mat <- X_mat_qpcr
#       B_bdtemp <- B_qpcr_mat[match(colnames(X_mat), rownames(B_qpcr_mat)),]
# 
#     } else {
#       y_mat_temp <- y_mat_qpcrbymicro
#       y_lab <- "Residual variation \nin Bd (residual qPCR)"
#       X_mat <- X_mat_qpcrbymicro
#       B_bdtemp <- B_qpcrbymicro_mat[match(colnames(X_mat), rownames(B_qpcrbymicro_mat)),]
#     }
#     # aesthetics for plotting x axes
#     if (x == "qPCR_bact_log10") {
#       x_lab <- "Bacterial density \n(log10 qPCR)"
#     } else if (x == "CV_log10") {
#       x_lab <- "Biofilm thickness \n(log10 CV)"
#     } else if (x == "Richness") {
#       x_lab <- "Bacterial community \n Richness"
#     }
#     # Which to exclude to plot for marginal
#     pos <- c(grep("Intercept", colnames(X_mat)), grep(x, colnames(X_mat)))
#     # Get the marginal predictor values we're actually looking for
#     X_marg <- as.matrix(X_mat)[,c(pos)]
#     B_marg <- as.matrix(B_bdtemp[c(pos),])
#     # Get the marginal effects by getting pred y
#     X_temp <- as.matrix(X_mat)[,-c(pos)]
#     B_temp <- as.matrix(B_bdtemp[-c(pos),])
#     # colnames(X_temp) == rownames(B_temp)
#     y_marg <- X_temp %*% B_temp
#     y_res <- y_mat_temp - apply(y_marg,1,median)
#     # Get slope of the marginal effects
#     samps_non_slope <- B_marg[grep(paste0(x,"XNon"), rownames(B_marg)),]
#     samps_inhib_slope <- colSums(B_marg[grep(paste0(x), rownames(B_marg)),])
#     # Get intercept of overall
#     samps_non_inter <-  B_marg[grep(paste0("InterceptXNon"), rownames(B_marg)),] - samps_non_slope*dat_centres[which(dat_centres$Predictor==x), "Centre"]
#     samps_inhib_inter <- colSums(B_marg[grep(paste0("Intercept"), rownames(B_marg)),])- samps_inhib_slope*dat_centres[which(dat_centres$Predictor==x), "Centre"]
# 
#     # Adjust X_marg to be "real" values
#     X_marg_adj <- X_marg + dat[,grep(paste0("Centre_",x,"$"), colnames(dat))]
#     X_all_adjmarg <- cbind(X_temp, X_marg_adj, Inhibitory_designation = X_mat$InterceptXInhibEff) %>%
#       as_tibble() %>%
#       mutate(Inhibitory_designation = ifelse(Inhibitory_designation==1, "Inhibitory", "Non-inhibitory"))
# 
#     ### Summary sheets ###
#     ## Make datasheet
#     tempNewDat <- data.frame(y_res, X_all_adjmarg ) %>%
#       rename(X = paste0(x,"XNon"), Y = paste0(y))  %>%
#       select(X, Y, Inhibitory_designation) %>% mutate(Response = paste0(y), Predictor = paste0(x))
#     ## Make slopes and intercept
#     tempmb <- t(apply(B_marg,1,median)) %>% as.data.frame() %>%
#       mutate(InterceptXInhib = InterceptXNon + InterceptXInhibEff, temp = get(paste0(x,"XNon")) + get(paste0(x,"XInhibEff"))) %>% #rename(paste0(x,"XInhib") = temp) %>%
#       select(InterceptXNon, InterceptXInhib, paste0(x, "XNon"), temp) %>%
#       pivot_longer(everything() , names_to = "TEMP", values_to = "Estimate") %>%
#       mutate(TEMP = ifelse(TEMP == "temp", paste0(x,"XInhib"), TEMP)) %>%
#       separate(TEMP, into = c("Predictor","Inhibitory_designation"), sep = "X") %>% mutate(Inhibitory_designation = ifelse(Inhibitory_designation == "Inhib", "Inhibitory", "Non-inhibitory")) %>%
#       rownames_to_column(var="row") %>%
#       pivot_wider(-c("row"), names_from = Predictor, values_from = Estimate) %>%
#       rename(slope = paste0(x)) %>%
#       mutate(Response = paste0(y), Predictor = paste0(x)) %>%
#       left_join(dat_centres) %>% mutate(intercept = Intercept-slope*Centre) %>% select(-Centre, -Intercept)
# 
#     ### Simulation ###
#     xsim <- seq(min(X_marg_adj[,paste0(x,"XNon")]), max(X_marg_adj[,paste0(x,"XNon")]), length.out = 100)
#     # Non-inhibitory
#     non_sim_values <- matrix(sapply(1:length(samps_non_slope), function(x) samps_non_slope[x]*xsim + samps_non_inter[x]),ncol=length(samps_non_slope))
#     CredInterNon <- data.frame(lower95 = apply(non_sim_values, 1, function(x) ci(x, method="HDI")$CI_low), upper95 = apply(non_sim_values, 1, function(x) ci(x, method="HDI")$CI_high),xsim=xsim) %>%
#       mutate(InterceptXInhibEff = "Non-inhibitory")
#     # Inhibitory
#     inhib_sim_values <- matrix(sapply(1:length(samps_inhib_slope), function(x) samps_inhib_slope[x]*xsim + samps_inhib_inter[x]),ncol=length(samps_inhib_slope))
#     CredInterInhib <- data.frame(lower95 = apply(inhib_sim_values, 1, function(x) ci(x, method="HDI")$CI_low), upper95 = apply(inhib_sim_values, 1, function(x) ci(x, method="HDI")$CI_high),xsim=xsim) %>%
#       mutate(InterceptXInhibEff = "Inhibitory")
#     CredInterAll <- rbind(CredInterInhib, CredInterNon)
#     ## Make Cred datasheet
#     tempNewCred <- CredInterAll %>%
#       mutate(Response = paste0(y), Predictor = paste0(x))
# 
#     # Look at marginal effect
#     #data.frame(y_res, X_all_adjmarg ) %>%
#     #   mutate(InterceptXInhibEff = ifelse(InterceptXInhibEff==1, "Inhibitory","Non-inhibitory"))
#     gg_temp <-  tempNewDat %>% #ggplot(aes(x=get(paste0(x,"XNon")), y=get(y), col=factor(InterceptXInhibEff))) +
#       ggplot(aes(x=X, y=Y, col=Inhibitory_designation)) +
#       geom_point() +
#       geom_ribbon(data = tempNewCred, aes(x=xsim, ymin=lower95, ymax=upper95, fill=InterceptXInhibEff), alpha=0.1, inherit.aes = FALSE, show.legend = FALSE) +
#       geom_abline(data=tempmb, aes(intercept=intercept, slope = slope, col = Inhibitory_designation)) +
#       labs(col="Biofilm type") +xlab(x_lab) +ylab(y_lab)
#     # gg_temp
#     ggsave(filename = paste0("03_figure_generation/marginal_", y,"_",x,".png"), width=4, height=3
#            ,gg_temp)
#     # Add onto data frame
#     allMarginalForPlotting <- rbind(allMarginalForPlotting, tempNewDat)
#     allCredForPlotting <- rbind(allCredForPlotting, tempNewCred)
#     allSlopesForPlotting <- rbind(allSlopesForPlotting, tempmb)
#   }
# }
# 

#### Without inhibitory marginals #####
B_micro_noInhib_mat <- samps_bayes_noInhib %>% filter(response == "Bd_micro_log10", GroupOrType !="Inhib") %>%
  select(Coefficient, Estimate) %>% group_by(Coefficient) %>% mutate(iter = seq(1, n())) %>%
  pivot_wider(names_from="Coefficient", values_from = "Estimate") %>% select(-iter) %>% t()
B_qpcr_noInhib_mat <- samps_bayes_noInhib %>% filter(response == "Bd_qPCR_log10", GroupOrType !="Inhib") %>%
  select(Coefficient, Estimate) %>% group_by(Coefficient) %>% mutate(iter = seq(1, n())) %>%
  pivot_wider(names_from="Coefficient", values_from = "Estimate") %>% select(-iter) %>% t()
B_qpcrbymicro_noInhib_mat <- samps_bayes_qpcrbymicr_noInhib %>% filter(GroupOrType !="Inhib") %>%
  select(Coefficient, Estimate) %>% group_by(Coefficient) %>% mutate(iter = seq(1, n())) %>%
  pivot_wider(names_from="Coefficient", values_from = "Estimate") %>% select(-iter) %>% t()
X_mat_micro_noInhib  <- dat %>%
  mutate(D = 1) %>% pivot_wider(names_from = DateStart, values_from = D, values_fill=0) %>%
  mutate(Intercept = 1, Richness = Richness, CV_log10 = CV_log10, qPCR_bact_log10 = qPCR_bact_log10) %>%
  select(one_of(rownames(B_qpcr_noInhib_mat)))
X_mat_qpcr_noInhib  <- dat %>%
  mutate(D = 1) %>% pivot_wider(names_from = DateStart, values_from = D, values_fill=0) %>%
  mutate(Intercept = 1, Richness = Richness, CV_log10 = CV_log10, qPCR_bact_log10 = qPCR_bact_log10) %>%
  select(one_of(rownames(B_qpcr_noInhib_mat)))
X_mat_qpcrbymicro_noInhib  <- dat %>%
  mutate(D = 1) %>% pivot_wider(names_from = DateStart, values_from = D, values_fill=0) %>%
  mutate(Intercept = 1, Richness = Richness, CV_log10 = CV_log10, qPCR_bact_log10 = qPCR_bact_log10) %>%
  select(one_of(rownames(B_qpcrbymicro_noInhib_mat)))
X_Inhibitory <- dat %>%
  mutate(D = 1) %>% pivot_wider(names_from = DateStart, values_from = D, values_fill=0) %>%
  mutate(Inhibitory_designation = ifelse(Inhibitory==1, "Inhibitory", "Non-inhibitory")) %>% select(Inhibitory_designation)

allMarginalForPlotting_noInhib <- data.frame()
allCredForPlotting_noInhib <- data.frame()
allSlopesForPlotting_noInhib <- data.frame()
#### Intercepts NOT included in residual bd
for ( x in c("Richness","CV_log10","qPCR_bact_log10")) {
  # x="Richness"
  for ( y in c("Bd_micro_log10", "Bd_qPCR_log10", "qPCRbyMicro")) {
  # for ( y in c("Bd_micro_log10", "Bd_qPCR_log10")) {
    # y="Bd_qPCR_log10"
    # Calling specific matrices; setting y labels
    if ( y == "Bd_micro_log10") {
      # B_bdtemp <- B_micro
      y_mat_temp <- y_mat_micro
      y_lab <- "Residual variation \nin Bd (counts)"
      X_mat <- X_mat_micro_noInhib
      B_bdtemp <- B_micro_noInhib_mat[match(colnames(X_mat), rownames(B_micro_noInhib_mat)),]

    } else if (y=="Bd_qPCR_log10") {
      # B_bdtemp <- B_qpcr
      y_mat_temp <- y_mat_qpcr
      y_lab <- "Residual variation \nin Bd (qPCR)"
      X_mat <- X_mat_qpcr_noInhib
      B_bdtemp <- B_qpcr_noInhib_mat[match(colnames(X_mat), rownames(B_qpcr_noInhib_mat)),]

    } else {
      y_mat_temp <- y_mat_qpcrbymicro
      y_lab <- "Residual variation \nin Bd (residual qPCR)"
      X_mat <- X_mat_qpcrbymicro_noInhib
      B_bdtemp <- B_qpcrbymicro_noInhib_mat[match(colnames(X_mat), rownames(B_qpcrbymicro_noInhib_mat)),]
    }
    # aesthetics for plotting x axes
    if (x == "qPCR_bact_log10") {
      x_lab <- "Bacterial density \n(log10 qPCR)"
    } else if (x == "CV_log10") {
      x_lab <- "Biofilm thickness \n(log10 CV)"
    } else if (x == "Richness") {
      x_lab <- "Bacterial community \n Richness"
    }
    # Which to exclude to plot for marginal
    pos <- c(grep("Intercept", colnames(X_mat)), grep(x, colnames(X_mat)))
    # Get the marginal predictor values we're actually looking for
    X_marg <- as.matrix(X_mat)[,c(pos)]
    B_marg <- as.matrix(B_bdtemp[c(pos),])
    # Get the marginal effects by getting pred y
    X_temp <- as.matrix(X_mat)[,-c(pos)]
    B_temp <- as.matrix(B_bdtemp[-c(pos),])
    # colnames(X_temp) == rownames(B_temp)
    y_marg <- X_temp %*% B_temp
    y_res <- y_mat_temp - apply(y_marg,1,median)
    # Get slope of the marginal effects
    samps_slope <- B_marg[grep(paste0(x), rownames(B_marg)),]
    # Get intercept of overall
    samps_inter <-  B_marg[grep(paste0("Intercept"), rownames(B_marg)),] - samps_slope*dat_centres[which(dat_centres$Predictor==x), "Centre"]

    # Adjust X_marg to be "real" values
    X_marg_adj <- X_marg + dat[,grep(paste0("Centre_",x,"$"), colnames(dat))]
    X_all_adjmarg <- cbind(X_temp, X_marg_adj) %>%
      as_tibble()
    ### Summary sheets ###
    ## Make datasheet
    tempNewDat <- data.frame(y_res, X_all_adjmarg,X_Inhibitory) %>%
      rename(X = paste0(x), Y = paste0(y))  %>%
      select(X, Y, Inhibitory_designation) %>% mutate(Response = paste0(y), Predictor = paste0(x))
    ## Make slopes and intercept
    tempmb <- t(apply(B_marg,1,median)) %>% as.data.frame() %>%
      # mutate(InterceptXInhib = InterceptXNon + InterceptXInhibEff, temp = get(paste0(x,"XNon")) + get(paste0(x,"XInhibEff"))) %>% #rename(paste0(x,"XInhib") = temp) %>%
      # select(InterceptXNon, InterceptXInhib, paste0(x, "XNon"), temp) %>%
      # pivot_longer(everything() , names_to = "TEMP", values_to = "Estimate") %>%
      # mutate(TEMP = ifelse(TEMP == "temp", paste0(x,"XInhib"), TEMP)) %>%
      # separate(TEMP, into = c("Predictor","Inhibitory_designation"), sep = "X") %>% mutate(Inhibitory_designation = ifelse(Inhibitory_designation == "Inhib", "Inhibitory", "Non-inhibitory")) %>%
      # rownames_to_column(var="row") %>%
      # pivot_wider(-c("row"), names_from = Predictor, values_from = Estimate) %>%
      rename(slope = paste0(x)) %>%
      mutate(Response = paste0(y), Predictor = paste0(x)) %>%
      left_join(dat_centres) %>% mutate(intercept = Intercept-slope*Centre) %>% select(-Centre, -Intercept)

    ### Simulation ###
    xsim <- seq(min(X_marg_adj[,paste0(x)]), max(X_marg_adj[,paste0(x)]), length.out = 100)
    # Non-inhibitory
    sim_values <- matrix(sapply(1:length(samps_slope), function(x) samps_slope[x]*xsim + samps_inter[x]),ncol=length(samps_slope))
    CredInter <- data.frame(lower95 = apply(sim_values, 1, function(x) ci(x, method="HDI")$CI_low), upper95 = apply(sim_values, 1, function(x) ci(x, method="HDI")$CI_high),xsim=xsim)
    # # Inhibitory
    # inhib_sim_values <- matrix(sapply(1:length(samps_inhib_slope), function(x) samps_inhib_slope[x]*xsim + samps_inhib_inter[x]),ncol=length(samps_inhib_slope))
    # CredInterInhib <- data.frame(lower95 = apply(inhib_sim_values, 1, function(x) ci(x, method="HDI")$CI_low), upper95 = apply(inhib_sim_values, 1, function(x) ci(x, method="HDI")$CI_high),xsim=xsim) %>%
    #   mutate(InterceptXInhibEff = "Inhibitory")
    # CredInterAll <- rbind(CredInterInhib, CredInterNon)
    ## Make Cred datasheet
    tempNewCred <- CredInter %>%
      mutate(Response = paste0(y), Predictor = paste0(x))

    # Look at marginal effect
    #data.frame(y_res, X_all_adjmarg ) %>%
    #   mutate(InterceptXInhibEff = ifelse(InterceptXInhibEff==1, "Inhibitory","Non-inhibitory"))
    gg_temp <-  tempNewDat %>% #ggplot(aes(x=get(paste0(x,"XNon")), y=get(y), col=factor(InterceptXInhibEff))) +
      ggplot(aes(x=X, y=Y)) +
      geom_point() +
      geom_ribbon(data = tempNewCred, aes(x=xsim, ymin=lower95, ymax=upper95), alpha=0.1, inherit.aes = FALSE, show.legend = FALSE) +
      geom_abline(data=tempmb, aes(intercept=intercept, slope = slope)) +
      labs(col="Biofilm type") +xlab(x_lab) +ylab(y_lab)
    # gg_temp
    ggsave(filename = paste0("03_figure_generation/marginal_", y,"_",x,"_noInhib.png"), width=4, height=3
           ,gg_temp)
    # Add onto data frame
    allMarginalForPlotting_noInhib <- rbind(allMarginalForPlotting_noInhib, tempNewDat)
    allCredForPlotting_noInhib <- rbind(allCredForPlotting_noInhib, tempNewCred)
    allSlopesForPlotting_noInhib <- rbind(allSlopesForPlotting_noInhib, tempmb)
  }
}

#
# #### Intercepts already incorporated in residual bd
# allMarginalForPlotting <- data.frame()
# allCredForPlotting <- data.frame()
# allSlopesForPlotting <- data.frame()
# for ( x in c("Richness","CV_log10","qPCR_bact_log10")) {
#   # x="Richness"
#   for ( y in c("Bd_micro_log10", "Bd_qPCR_log10", "qPCRbyMicro")) {
#     # y="Bd_micro_log10"
#     # Calling specific matrices; setting y labels
#     if ( y == "Bd_micro_log10") {
#       # B_bdtemp <- B_micro
#       y_mat_temp <- y_mat_micro
#       y_lab <- "Residual variation \nin Bd (counts)"
#       X_mat <- X_mat_micro
#       B_bdtemp <- B_micro_mat[match(colnames(X_mat), rownames(B_micro_mat)),]
#       # y_title <- y
#
#     } else if (y=="Bd_qPCR_log10") {
#       # B_bdtemp <- B_qpcr
#       y_mat_temp <- y_mat_qpcr
#       y_lab <- "Residual variation \nin Bd (qPCR)"
#       X_mat <- X_mat_qpcr
#       B_bdtemp <- B_qpcr_mat[match(colnames(X_mat), rownames(B_qpcr_mat)),]
#       # y_title <- y
#
#     } else {
#       y_mat_temp <- y_mat_qpcrbymicro
#       y_lab <- "Residual variation \nin Bd (residual qPCR)"
#       X_mat <- X_mat_qpcrbymicro
#       B_bdtemp <- B_qpcrbymicro_mat[match(colnames(X_mat), rownames(B_qpcrbymicro_mat)),]
#       # y_title <- y
#       # y <- "Bd_qPCR_log10"
#     }
#     # aesthetics for plotting x axes
#     if (x == "qPCR_bact_log10") {
#       x_lab <- "Bacterial density \n(log10 qPCR)"
#     } else if (x == "CV_log10") {
#       x_lab <- "Biofilm thickness \n(log10 CV)"
#     } else if (x == "Richness") {
#       x_lab <- "Bacterial community \n Richness"
#     }
#     # Which to exclude to plot for marginal
#     pos <- grep(x, colnames(X_mat))
#     # Get the marginal predictor values we're actually looking for
#     X_marg <- as.matrix(X_mat)[,c(pos)]
#     B_marg <- as.matrix(B_bdtemp[c(pos),])
#     # Get the marginal effects by getting pred y
#     X_temp <- as.matrix(X_mat)[,-c(pos)]
#     B_temp <- as.matrix(B_bdtemp[-c(pos),])
#     # colnames(X_temp) == rownames(B_temp)
#     y_marg <- X_temp %*% B_temp
#     y_res <- y_mat_temp - apply(y_marg,1,median)
#     # Get slope of the marginal effects
#     samps_non_slope <- B_marg[grep("Non", rownames(B_marg)),]
#     samps_inhib_slope <- colSums(B_marg)
#     # Adjust X_marg to be "real" values
#     X_marg_adj <- X_marg + dat[,grep(paste0("Centre_",x), colnames(dat))]
#     X_all_adjmarg <- cbind(X_temp, X_marg_adj)
#
#     ### Summary sheets ###
#     ## Make datasheet
#     tempNewDat <- data.frame(y_res, X_all_adjmarg ) %>%
#     mutate(InterceptXInhibEff = ifelse(InterceptXInhibEff==1, "Inhibitory","Non-inhibitory")) %>%
#       rename(Inhibitory_designation = InterceptXInhibEff, X = paste0(x,"XNon"), Y = paste0(y))  %>%
#       select(X, Y, Inhibitory_designation) %>% mutate(Response = paste0(y), Predictor = paste0(x))
#     # tempNewDat <- data.frame(y_res, X_mat ) %>%
#     #   mutate(InterceptXInhibEff = ifelse(InterceptXInhibEff==1, "Inhibitory","Non-inhibitory")) %>%
#     #   rename(Inhibitory_designation = InterceptXInhibEff, X = paste0(x,"XNon"), Y = paste0(y))  %>%
#     #   select(X, Y, Inhibitory_designation) %>% mutate(Response = paste0(y), Predictor = paste0(x))
#     ## Make slopes and intercept
#     tempmb <- t(apply(B_marg,1,median)) %>% as.data.frame() %>%
#       mutate(Inhibitory = get(paste0(x,"XNon")) + get(paste0(x,"XInhibEff"))) %>% rename("Non-inhibitory" = paste0(x,"XNon")) %>%
#       select("Non-inhibitory", "Inhibitory") %>%
#       mutate(Response = paste0(y), Predictor = paste0(x)) %>%
#       pivot_longer(-c(Response, Predictor), names_to = "Inhibitory_designation", values_to = "slope") %>%
#       left_join(dat_centres) %>% mutate(intercept = 0-slope*Centre) %>% select(-Centre)
#
#     ### Simulation ###
#     xsim <- seq(min(X_marg_adj[,paste0(x,"XNon")]), max(X_marg_adj[,paste0(x,"XNon")]), length.out = 100)
#     intercept_non <- tempmb %>% filter(Inhibitory_designation == "Non-inhibitory") %>% pull(intercept)
#     intercept_inhib <- tempmb %>% filter(Inhibitory_designation == "Inhibitory") %>% pull(intercept)
#     # Non-inhibitory
#     non_sim_values <- matrix(sapply(1:length(samps_non_slope), function(x) samps_non_slope[x]*xsim + intercept_non),ncol=length(samps_non_slope))
#     CredInterNon <- data.frame(lower95 = apply(non_sim_values, 1, function(x) ci(x, method="HDI")$CI_low), upper95 = apply(non_sim_values, 1, function(x) ci(x, method="HDI")$CI_high),xsim=xsim) %>%
#       mutate(InterceptXInhibEff = "Non-inhibitory")
#     # Inhibitory
#     inhib_sim_values <- matrix(sapply(1:length(samps_inhib_slope), function(x) samps_inhib_slope[x]*xsim + intercept_inhib),ncol=length(samps_inhib_slope))
#     CredInterInhib <- data.frame(lower95 = apply(inhib_sim_values, 1, function(x) ci(x, method="HDI")$CI_low), upper95 = apply(inhib_sim_values, 1, function(x) ci(x, method="HDI")$CI_high),xsim=xsim) %>%
#       mutate(InterceptXInhibEff = "Inhibitory")
#     CredInterAll <- rbind(CredInterInhib, CredInterNon)
#     ## Make Cred datasheet
#     tempNewCred <- CredInterAll %>%
#       mutate(Response = paste0(y), Predictor = paste0(x))
#
#     # Look at marginal effect
#     #data.frame(y_res, X_all_adjmarg ) %>%
#     #   mutate(InterceptXInhibEff = ifelse(InterceptXInhibEff==1, "Inhibitory","Non-inhibitory"))
#     gg_temp <-  tempNewDat %>% #ggplot(aes(x=get(paste0(x,"XNon")), y=get(y), col=factor(InterceptXInhibEff))) +
#       ggplot(aes(x=X, y=Y, col=Inhibitory_designation)) +
#       geom_point() +
#       # geom_smooth(method = "lm") +
#       geom_ribbon(data = tempNewCred, aes(x=xsim, ymin=lower95, ymax=upper95, fill=InterceptXInhibEff), alpha=0.1, inherit.aes = FALSE, show.legend = FALSE) +
#       geom_abline(data=tempmb, aes(intercept=intercept, slope = slope, col = Inhibitory_designation)) +
#       geom_abline(data=tempmb,aes(intercept=intercept, slope = slope, col = Inhibitory_designation)) +
#       labs(col="Biofilm type") +xlab(x_lab) +ylab(y_lab)
#    # gg_temp
#    ggsave(filename = paste0("03_figure_generation/marginal_", y,"_",x,".png"), width=4, height=3
#      ,gg_temp)
#    # Add onto data frame
#    allMarginalForPlotting <- rbind(allMarginalForPlotting, tempNewDat)
#    allCredForPlotting <- rbind(allCredForPlotting, tempNewCred)
#    allSlopesForPlotting <- rbind(allSlopesForPlotting, tempmb)
#   }
# }

#### Multi-panel with Inhib ####

# Relabelling facets, transforming x back into un-centred
pred.labs <- c(Richness ="Bacterial richness \n(# ASVs)"
               , qPCR_bact_log10 ="Bacterial cell density \n(log10 qPCR)"
               , CV_log10="Biofilm thickness \n(log10 CV)")
resp.labs <- c(Bd_qPCR_log10= "Bd qPCR copies (log10)"
               , Bd_micro_log10 = "Bd microscopy counts (log10)"
               , qPCRbyMicro = "Bd qPCR copies (log10)\n Microscopy count-adjusted")
# Change factor level order
allMarginalForPlotting <- allMarginalForPlotting %>%
  mutate(Predictor = factor(Predictor, levels = c("Richness","qPCR_bact_log10", "CV_log10")))
allCredForPlotting <- allCredForPlotting%>%
  mutate(Predictor = factor(Predictor, levels = c("Richness","qPCR_bact_log10", "CV_log10")))
allSlopesForPlotting <- allSlopesForPlotting%>%
  mutate(Predictor = factor(Predictor, levels = c("Richness","qPCR_bact_log10", "CV_log10")))

# Make credible interval for lines; first, get samps for intercepts
# samps_bayes %>% filter(GroupOrType != "InhibEff")
#
# summary_bayes_qpcrbymicr
# allSlopesForPlotting_wsign <- full_join(summary_bayes,summary_bayes_qpcrbymicr) %>% rename(Response=response) %>%filter(GroupOrType %in% c("Non","Inhib")) %>%
#   mutate(significant = PD<0.025) %>%
#   mutate(posProb = factor(ifelse(PD<0.025, ">95%", ifelse(PD<0.05, ">90%", "<=89%")), levels = c(">95%",">90%", "<=89%"))) %>%
#   mutate(significance = factor(ifelse(PD<0.025, "95% CredI", ifelse(PD<0.05, "90% CredI", "<90%")), levels=c("95% CredI", "90% CredI", "<90%"))) %>%
#   select(Coefficient, Response, significant,significance,posProb) %>% separate(Coefficient, into=c("Predictor","Inhibitory_designation"), sep="X") %>%
#   mutate(Inhibitory_designation = ifelse(Inhibitory_designation=="Inhib", "Inhibitory", "Non-inhibitory")) %>%
#   right_join(allSlopesForPlotting)%>% mutate(Predictor = factor(Predictor, levels = c("Richness","qPCR_bact_log10", "CV_log10")))

allSlopesForPlotting_wsign <- full_join(summary_bayes,summary_bayes_qpcrbymicr) %>% rename(Response=response) %>%filter(GroupOrType %in% c("Non","Inhib")) %>%
  mutate(significant = PD<0.025) %>%
  mutate(posProb = factor(ifelse(PD<0.025, ">95%", ifelse(PD<0.05, ">90%", "<=89%")), levels = c(">95%",">90%", "<=89%"))) %>%
  mutate(significance = factor(ifelse(PD<0.025, "95% CredI", ifelse(PD<0.05, "90% CredI", "<90%")), levels=c("95% CredI", "90% CredI", "<90%"))) %>%
  select(Coefficient, Response, significant,significance,posProb) %>% separate(Coefficient, into=c("Predictor","Inhibitory_designation"), sep="X") %>%
  mutate(Inhibitory_designation = ifelse(Inhibitory_designation=="Inhib", 1, 0)) %>%
  right_join(allSlopesForPlotting)%>% mutate(Predictor = factor(Predictor, levels = c("Richness","qPCR_bact_log10", "CV_log10")))

##### all marginal effects from all models
gg_allmarg <- allMarginalForPlotting %>%
  ggplot(aes(x=X, y=Y)) +
  geom_point(aes( col = Inhibitory_fraction)) +
  #geom_smooth(method="lm") +
  ylab("Residual variation in Bd") +
  xlab("Predictor variable being tested for marginal effects") +
  geom_ribbon(data = allCredForPlotting, aes(x=xsim, ymin=lower95, ymax=upper95, fill=factor(InterceptXInhibEff)), alpha=0.1, inherit.aes = FALSE, show.legend = FALSE) +
  scale_fill_manual(values=c("darkgrey","salmon"))+
  geom_abline(data = allSlopesForPlotting_wsign, aes(intercept=intercept, slope = slope, col = Inhibitory_designation, lty = posProb)) +
  scale_color_gradient(low="darkgrey",high="salmon") +
  # geom_abline(data=allSlopesForPlotting, aes(intercept=0, slope = B_marg[1] + B_marg[2], col = "Inhibitory")) +
  facet_grid(Response ~ Predictor, scales = "free", switch = "both", labeller=labeller(Predictor=pred.labs, Response=resp.labs)) +
  scale_linetype_manual(values = c(1,2,3))  +
  # facet_grid(factor(Response, levels=c("Bd_micro_log10", "Bd_qPCR_log10"))~factor(Predictor, levels=c("Richness","qPCR_bact_log10","CV_log10")), scales = "free", switch = "both", labeller=labeller(Predictor=pred.labs, Response=resp.labs)) +
  labs(col="Biofilm\nfraction inhibitory", lty="Posterior probability\n (Probability that \nslope is > 0)")+
  theme_bw()
gg_allmarg
ggsave(filename = "03_figure_generation/marginal_ALL.png", height=7, width=8
       ,gg_allmarg)

##### Two main responses: micro an qpcr
gg_allmarg_twomain <- allMarginalForPlotting %>%
  filter(Response != "qPCRbyMicro") %>%
  ggplot(aes(x=X, y=Y, col =Inhibitory_fraction)) +
  geom_point() +
  #geom_smooth(method="lm") +
  ylab("Residual variation in Bd") +
  xlab("Predictor variable being tested for marginal effects") +
  geom_ribbon(data = allCredForPlotting %>%filter(Response != "qPCRbyMicro"), aes(x=xsim, ymin=lower95, ymax=upper95, fill=factor(InterceptXInhibEff)), alpha=0.1, inherit.aes = FALSE, show.legend = FALSE) +
  scale_fill_manual(values=c("darkgrey","salmon")) +
  geom_abline(data = allSlopesForPlotting_wsign%>%filter(Response != "qPCRbyMicro"), aes(intercept=intercept, slope = slope, col = Inhibitory_designation, lty = posProb)) +
  scale_color_gradient(low="darkgrey", high="salmon") +
  # geom_abline(data=allSlopesForPlotting, aes(intercept=0, slope = B_marg[1] + B_marg[2], col = "Inhibitory")) +
  facet_grid(Response ~ Predictor, scales = "free", switch = "both", labeller=labeller(Predictor=pred.labs, Response=resp.labs)) +
  scale_linetype_manual(values = c(1,2,3))  +
  # facet_grid(factor(Response, levels=c("Bd_micro_log10", "Bd_qPCR_log10"))~factor(Predictor, levels=c("Richness","qPCR_bact_log10","CV_log10")), scales = "free", switch = "both", labeller=labeller(Predictor=pred.labs, Response=resp.labs)) +
  labs(col="Biofilm\nfraction inhibitory", lty="Posterior probability\n (Probability that \nslope is > 0)") +
  theme_bw()
gg_allmarg_twomain

ggsave(filename = "03_figure_generation/marginal_ALL_micro_qPCR.png", height=5, width=8
       ,gg_allmarg_twomain)


  
##### Micro and qPCR adjust only

gg_allmarg_micro_and_qpcradj <- allMarginalForPlotting %>%
  filter(Response != "Bd_qPCR_log10") %>%
  ggplot(aes(x=X, y=Y, col =Inhibitory_fraction)) +
  geom_point() +
  #geom_smooth(method="lm") +
  ylab("Residual variation in Bd") +
  xlab("Predictor variable being tested for marginal effects") +
  geom_ribbon(data = allCredForPlotting %>% filter(Response != "Bd_qPCR_log10"), aes(x=xsim, ymin=lower95, ymax=upper95, fill=factor(InterceptXInhibEff)), alpha=0.1, inherit.aes = FALSE, show.legend = FALSE) +
  geom_abline(data = allSlopesForPlotting_wsign %>% filter(Response != "Bd_qPCR_log10"), aes(intercept=intercept, slope = slope, col = Inhibitory_designation, lty = posProb)) +
  # geom_abline(data=allSlopesForPlotting, aes(intercept=0, slope = B_marg[1] + B_marg[2], col = "Inhibitory")) +
  facet_grid(Response ~ Predictor, scales = "free", switch = "both", labeller=labeller(Predictor=pred.labs, Response=resp.labs)) +
  scale_linetype_manual(values = c(1,2,3))  +
  scale_fill_manual(values=c("darkgrey","salmon")) +
  scale_color_gradient(low="darkgrey", high="salmon") +
  # facet_grid(factor(Response, levels=c("Bd_micro_log10", "Bd_qPCR_log10"))~factor(Predictor, levels=c("Richness","qPCR_bact_log10","CV_log10")), scales = "free", switch = "both", labeller=labeller(Predictor=pred.labs, Response=resp.labs)) +
  labs(col="Biofilm designation", lty="Posterior probability\n(Probability that \nslope is > 0)")+
  theme_bw()

ggsave(filename = "03_figure_generation/marginal_ALL_micro_qPCRadj.png", height=5, width=8
       ,gg_allmarg_micro_and_qpcradj)



#### Multi-panel, NO Inhib ####

# Relabelling facets, transforming x back into un-centred
pred.labs <- c(Richness ="Bacterial richness \n(# ASVs)"
               , qPCR_bact_log10 ="Bacterial cell density \n(log10 qPCR)"
               , CV_log10="Biofilm thickness \n(log10 CV)")
resp.labs <- c(Bd_qPCR_log10= "Bd qPCR copies (log10)"
               , Bd_micro_log10 = "Bd microscopy counts (log10)"
               , qPCRbyMicro = "Bd qPCR copies (log10)\n Microscopy count-adjusted")
# Change factor level order
allMarginalForPlotting_noInhib <- allMarginalForPlotting_noInhib %>%
  mutate(Predictor = factor(Predictor, levels = c("Richness","qPCR_bact_log10", "CV_log10")))
allCredForPlotting_noInhib <- allCredForPlotting_noInhib%>%
  mutate(Predictor = factor(Predictor, levels = c("Richness","qPCR_bact_log10", "CV_log10")))
allSlopesForPlotting_noInhib <- allSlopesForPlotting_noInhib%>%
  mutate(Predictor = factor(Predictor, levels = c("Richness","qPCR_bact_log10", "CV_log10")))

# Make credible interval for lines; first, get samps for intercepts
# samps_bayes %>% filter(GroupOrType != "InhibEff")
#
# summary_bayes_qpcrbymicr
allSlopesForPlotting_wsign_noInhib <- full_join(summary_bayes_noInhib,summary_bayes_qpcrbymicr_noInhib)  %>%
  rename(Response=response) %>% #filter(GroupOrType %in% c("Non","Inhib")) %>%
  mutate(significant = PD<0.025) %>%
  mutate(posProb = factor(ifelse(PD<0.025, ">95%", ifelse(PD<0.05, ">90%", "<=89%")), levels = c(">95%",">90%", "<=89%"))) %>%
  mutate(significance = factor(ifelse(PD<0.025, "95% CredI", ifelse(PD<0.05, "90% CredI", "<90%")), levels=c("95% CredI", "90% CredI", "<90%"))) %>%
  select(Predictor, Response, significant,significance,posProb) %>% #separate(Coefficient, into=c("Predictor","Inhibitory_designation"), sep="X") %>%
  # filter(Predictor %in% c("Richness","CV_log10","qPCR_bact_log10")) %>%
#mutate(Inhibitory_designation = ifelse(Inhibitory_designation=="Inhib", "Inhibitory", "Non-inhibitory")) %>%
  right_join(allSlopesForPlotting_noInhib)%>% mutate(Predictor = factor(Predictor, levels = c("Richness","qPCR_bact_log10", "CV_log10")))

##### all marginal effects from all models
gg_allmarg_noInhib_nocol <- allMarginalForPlotting_noInhib %>%
  ggplot(aes(x=X, y=Y)) +
  geom_point() +
  #geom_smooth(method="lm") +
  ylab("Residual variation in Bd") +
  xlab("Predictor variable being tested for marginal effects") +
  geom_ribbon(data = allCredForPlotting_noInhib, aes(x=xsim, ymin=lower95, ymax=upper95), alpha=0.1, inherit.aes = FALSE, show.legend = FALSE) +
  geom_abline(data = allSlopesForPlotting_wsign_noInhib, aes(intercept=intercept, slope = slope, lty = posProb)) +
  # geom_abline(data=allSlopesForPlotting, aes(intercept=0, slope = B_marg[1] + B_marg[2], col = "Inhibitory")) +
  facet_grid(Response ~ Predictor, scales = "free", switch = "both", labeller=labeller(Predictor=pred.labs, Response=resp.labs)) +
  scale_linetype_manual(values = c(1,2,3))  +
  labs(col="Biofilm designation", lty="Posterior probability\n (Probability that \nslope is > 0)")

ggsave(filename = "03_figure_generation/marginal_ALL_noInhib_nocolour.png", height=5, width=8
       ,gg_allmarg_noInhib_nocol)

gg_allmarg_noInhib_withcol <- allMarginalForPlotting_noInhib %>%
  ggplot(aes(x=X, y=Y, col =Inhibitory_designation)) +
  geom_point() +
  #geom_smooth(method="lm") +
  ylab("Residual variation in Bd") +
  xlab("Predictor variable being tested for marginal effects") +
  # geom_ribbon(data = allCredForPlotting_noInhib, aes(x=xsim, ymin=lower95, ymax=upper95), alpha=0.1, inherit.aes = FALSE, show.legend = FALSE) +
  geom_abline(data = allSlopesForPlotting_wsign_noInhib, aes(intercept=intercept, slope = slope, lty = posProb)) +
  geom_smooth(method = "lm") +
  # geom_abline(data=allSlopesForPlotting, aes(intercept=0, slope = B_marg[1] + B_marg[2], col = "Inhibitory")) +
  facet_grid(Response ~ Predictor, scales = "free", switch = "both", labeller=labeller(Predictor=pred.labs, Response=resp.labs)) +
  scale_linetype_manual(values = c(1,2,3))  +
  # facet_grid(factor(Response, levels=c("Bd_micro_log10", "Bd_qPCR_log10"))~factor(Predictor, levels=c("Richness","qPCR_bact_log10","CV_log10")), scales = "free", switch = "both", labeller=labeller(Predictor=pred.labs, Response=resp.labs)) +
  labs(col="Biofilm designation", lty="Posterior probability\n (Probability that \nslope is > 0)")

ggsave(filename = "03_figure_generation/marginal_ALL_noInhib_withcolour.png", height=5, width=8
       ,gg_allmarg_noInhib_withcol)


##### Two main responses: micro an qpcr
gg_allmarg_twomain_noInhib_nocol <- allMarginalForPlotting_noInhib %>%
  filter(Response != "qPCRbyMicro") %>%
  ggplot(aes(x=X, y=Y)) +
  geom_point() +
  #geom_smooth(method="lm") +
  ylab("Residual variation in Bd") +
  xlab("Predictor variable being tested for marginal effects") +
  geom_ribbon(data = allCredForPlotting_noInhib %>%filter(Response != "qPCRbyMicro"), aes(x=xsim, ymin=lower95, ymax=upper95), alpha=0.1, inherit.aes = FALSE, show.legend = FALSE) +
  geom_abline(data = allSlopesForPlotting_wsign_noInhib%>%filter(Response != "qPCRbyMicro"), aes(intercept=intercept, slope = slope, lty = posProb)) +
  # geom_abline(data=allSlopesForPlotting, aes(intercept=0, slope = B_marg[1] + B_marg[2], col = "Inhibitory")) +
  facet_grid(Response ~ Predictor, scales = "free", switch = "both", labeller=labeller(Predictor=pred.labs, Response=resp.labs)) +
  scale_linetype_manual(values = c(1,2,3))  +
  # facet_grid(factor(Response, levels=c("Bd_micro_log10", "Bd_qPCR_log10"))~factor(Predictor, levels=c("Richness","qPCR_bact_log10","CV_log10")), scales = "free", switch = "both", labeller=labeller(Predictor=pred.labs, Response=resp.labs)) +
  labs(col="Biofilm designation", lty="Posterior probability\n (Probability that \nslope is > 0)")

ggsave(filename = "03_figure_generation/marginal_ALL_micro_qPCR_noInhib_noCol.png", height=5, width=8
       ,gg_allmarg_twomain_noInhib_nocol)

gg_allmarg_twomain_noInhib_withcol <- allMarginalForPlotting_noInhib %>%
  filter(Response != "qPCRbyMicro") %>%
  ggplot(aes(x=X, y=Y, col=Inhibitory_designation)) +
  geom_point() +
  geom_smooth(method="lm") +
  ylab("Residual variation in Bd") +
  xlab("Predictor variable being tested for marginal effects") +
  # geom_ribbon(data = allCredForPlotting_noInhib %>%filter(Response != "qPCRbyMicro"), aes(x=xsim, ymin=lower95, ymax=upper95), alpha=0.1, inherit.aes = FALSE, show.legend = FALSE) +
  geom_abline(data = allSlopesForPlotting_wsign_noInhib%>%filter(Response != "qPCRbyMicro"), aes(intercept=intercept, slope = slope, lty = posProb)) +
  # geom_abline(data=allSlopesForPlotting, aes(intercept=0, slope = B_marg[1] + B_marg[2], col = "Inhibitory")) +
  facet_grid(Response ~ Predictor, scales = "free", switch = "both", labeller=labeller(Predictor=pred.labs, Response=resp.labs)) +
  scale_linetype_manual(values = c(1,2,3))  +
  # facet_grid(factor(Response, levels=c("Bd_micro_log10", "Bd_qPCR_log10"))~factor(Predictor, levels=c("Richness","qPCR_bact_log10","CV_log10")), scales = "free", switch = "both", labeller=labeller(Predictor=pred.labs, Response=resp.labs)) +
  labs(col="Biofilm designation", lty="Posterior probability\n (Probability that \nslope is > 0)")

ggsave(filename = "03_figure_generation/marginal_ALL_micro_qPCR_noInhib_withCol.png", height=5, width=8
       ,gg_allmarg_twomain_noInhib_withcol)

##### Micro and qPCR adjust only

gg_allmarg_micro_and_qpcradj_noInhib_nocol <- allMarginalForPlotting_noInhib %>%
  filter(Response != "Bd_qPCR_log10") %>%
  ggplot(aes(x=X, y=Y)) +
  geom_point() +
  #geom_smooth(method="lm") +
  ylab("Residual variation in Bd") +
  xlab("Predictor variable being tested for marginal effects") +
  geom_ribbon(data = allCredForPlotting_noInhib %>% filter(Response != "Bd_qPCR_log10"), aes(x=xsim, ymin=lower95, ymax=upper95), alpha=0.1, inherit.aes = FALSE, show.legend = FALSE) +
  geom_abline(data = allSlopesForPlotting_wsign_noInhib %>% filter(Response != "Bd_qPCR_log10"), aes(intercept=intercept, slope = slope, lty = posProb)) +
  # geom_abline(data=allSlopesForPlotting, aes(intercept=0, slope = B_marg[1] + B_marg[2], col = "Inhibitory")) +
  facet_grid(Response ~ Predictor, scales = "free", switch = "both", labeller=labeller(Predictor=pred.labs, Response=resp.labs)) +
  scale_linetype_manual(values = c(1,2,3))  +
  # facet_grid(factor(Response, levels=c("Bd_micro_log10", "Bd_qPCR_log10"))~factor(Predictor, levels=c("Richness","qPCR_bact_log10","CV_log10")), scales = "free", switch = "both", labeller=labeller(Predictor=pred.labs, Response=resp.labs)) +
  labs(col="Biofilm designation", lty="Posterior probability\n(Probability that \nslope is > 0)")

ggsave(filename = "03_figure_generation/marginal_ALL_micro_qPCRadj_noInhib_nocol.png", height=5, width=8
       ,gg_allmarg_micro_and_qpcradj_noInhib_nocol)



gg_allmarg_micro_and_qpcradj_noInhib_withcol <- allMarginalForPlotting_noInhib %>%
  filter(Response != "Bd_qPCR_log10") %>%
  ggplot(aes(x=X, y=Y, col =Inhibitory_designation)) +
  geom_point() +
  geom_smooth(method="lm") +
  ylab("Residual variation in Bd") +
  xlab("Predictor variable being tested for marginal effects") +
  # geom_ribbon(data = allCredForPlotting_noInhib %>% filter(Response != "Bd_qPCR_log10"), aes(x=xsim, ymin=lower95, ymax=upper95), alpha=0.1, inherit.aes = FALSE, show.legend = FALSE) +
  geom_abline(data = allSlopesForPlotting_wsign_noInhib %>% filter(Response != "Bd_qPCR_log10"), aes(intercept=intercept, slope = slope, lty = posProb)) +
  # geom_abline(data=allSlopesForPlotting, aes(intercept=0, slope = B_marg[1] + B_marg[2], col = "Inhibitory")) +
  facet_grid(Response ~ Predictor, scales = "free", switch = "both", labeller=labeller(Predictor=pred.labs, Response=resp.labs)) +
  scale_linetype_manual(values = c(1,2,3))  +
  # facet_grid(factor(Response, levels=c("Bd_micro_log10", "Bd_qPCR_log10"))~factor(Predictor, levels=c("Richness","qPCR_bact_log10","CV_log10")), scales = "free", switch = "both", labeller=labeller(Predictor=pred.labs, Response=resp.labs)) +
  labs(col="Biofilm designation", lty="Posterior probability\n(Probability that \nslope is > 0)")

ggsave(filename = "03_figure_generation/marginal_ALL_micro_qPCRadj_noInhib_withcol.png", height=5, width=8
       ,gg_allmarg_micro_and_qpcradj_noInhib_withcol)

#### FOR PRESENTATIONS: blank figure and simple figure

gg_BLANK <- allMarginalForPlotting_noInhib %>%
  filter(Response != "Bd_qPCR_log10") %>%
  ggplot(aes(x=X, y=Y)) +
  geom_point(alpha=0) +
  # geom_smooth(method="lm", alpha=0, se = FALSE) +
  ylab("Residual variation in Bd") +
  xlab("Predictor variable being tested for marginal effects") +
  # geom_ribbon(data = allCredForPlotting_noInhib %>% filter(Response != "Bd_qPCR_log10"), aes(x=xsim, ymin=lower95, ymax=upper95), alpha=0.1, inherit.aes = FALSE, show.legend = FALSE) +
  geom_abline(data = allSlopesForPlotting_wsign_noInhib %>% filter(Response != "Bd_qPCR_log10"), aes(intercept=intercept, slope = slope, lty = posProb), alpha=0) +
  # geom_abline(data=allSlopesForPlotting, aes(intercept=0, slope = B_marg[1] + B_marg[2], col = "Inhibitory")) +
  facet_grid(Response ~ Predictor, scales = "free", switch = "both", labeller=labeller(Predictor=pred.labs, Response=resp.labs)) +
  scale_linetype_manual(values = c(1,2,3))  +
  # facet_grid(factor(Response, levels=c("Bd_micro_log10", "Bd_qPCR_log10"))~factor(Predictor, levels=c("Richness","qPCR_bact_log10","CV_log10")), scales = "free", switch = "both", labeller=labeller(Predictor=pred.labs, Response=resp.labs)) +
  labs(col="Biofilm designation", lty="Posterior probability\n(Probability that \nslope is > 0)")
gg_BLANK


ggsave(filename = "03_figure_generation/marginal_ALL_micro_qPCRadj_noInhib_BLANK.png", height=5, width=8
       ,gg_BLANK)


gg_allmarg_micro_and_qpcradj_SIMPLE <- allMarginalForPlotting %>%
  filter(Response != "Bd_qPCR_log10") %>%
  ggplot(aes(x=X, y=Y, col =Inhibitory_designation)) +
  geom_point(alpha=0) +
  #geom_smooth(method="lm") +
  ylab("Residual variation in Bd") +
  xlab("Predictor variable being tested for marginal effects") +
  geom_ribbon(data = allCredForPlotting %>% filter(Response != "Bd_qPCR_log10"), aes(x=xsim, ymin=lower95, ymax=upper95, fill=InterceptXInhibEff), alpha=0, inherit.aes = FALSE, show.legend = FALSE) +
  geom_abline(data = allSlopesForPlotting_wsign %>% filter(Response != "Bd_qPCR_log10"), aes(intercept=intercept, slope = slope, col = Inhibitory_designation, lty = posProb)) +
  # geom_abline(data=allSlopesForPlotting, aes(intercept=0, slope = B_marg[1] + B_marg[2], col = "Inhibitory")) +
  facet_grid(Response ~ Predictor, scales = "free", switch = "both", labeller=labeller(Predictor=pred.labs, Response=resp.labs)) +
  scale_linetype_manual(values = c(1,2,3))  +
  # facet_grid(factor(Response, levels=c("Bd_micro_log10", "Bd_qPCR_log10"))~factor(Predictor, levels=c("Richness","qPCR_bact_log10","CV_log10")), scales = "free", switch = "both", labeller=labeller(Predictor=pred.labs, Response=resp.labs)) +
  labs(col="Biofilm designation", lty="Posterior probability\n(Probability that \nslope is > 0)")

ggsave(filename = "03_figure_generation/marginal_ALL_micro_qPCRadj_SIMPLE.png", height=5, width=8
       ,gg_allmarg_micro_and_qpcradj_SIMPLE)


#### /------------ Raw data plots ------------/ ####
# Relabelling facets, transforming x back into un-centred
pred.labs <- c(Richness ="Bacterial richness \n(# ASVs)"
               , qPCR_bact_log10 ="Bacterial cell density \n(log10 qPCR)"
               , CV_log10="Biofilm thickness \n(log10 CV)")
resp.labs <- c(Bd_qPCR_log10= "Bd qPCR copies (log10)"
               , Bd_micro_log10 = "Bd microscopy counts (log10)"
               , qPCRbyMicro = "Bd qPCR copies (log10)\n Microscopy count-adjusted")

gg_rawdat <- dat %>% 
  mutate(Richness = Richness + Centre_Richness, qPCR_bact_log10 = qPCR_bact_log10 + Centre_qPCR_bact_log10, CV_log10 = CV_log10 + Centre_CV_log10, inhibRichFrac = inhibRichFrac_raw) %>%
  select(Bd_micro_log10, Bd_qPCR_log10, Inhibitory, Richness, qPCR_bact_log10, CV_log10,inhibRichFrac) %>%
  pivot_longer(starts_with("Bd_"), names_to = "Response", values_to = "Y") %>%
  pivot_longer(c(Richness, qPCR_bact_log10, CV_log10), names_to = "Predictors", values_to = "X") %>%
  mutate(Inhibitory_designation = ifelse(inhibRichFrac > 0, "At least one\ninhibitory ASV", "Non-inhibitory ASVs")) %>%
  mutate(Predictors = factor(Predictors, levels = c("Richness","qPCR_bact_log10", "CV_log10")))  %>%
  ggplot(aes(x=X, y = Y, col = Inhibitory_designation)) + 
  geom_point() +
  geom_smooth(method = "lm", se=FALSE) +
  facet_grid(Response ~ Predictors, scales = "free", switch = "both",labeller=labeller(Predictors=pred.labs, Response=resp.labs)) +
  xlab("Predictors") + ylab("Response metrics") + labs(col = "Biofilm type")

ggsave("03_figure_generation/raw_dat_plot.png", height=5, width=8,
       gg_rawdat)

gg_rawdat_nozoer <- dat %>% 
  mutate(Richness = Richness + Centre_Richness, qPCR_bact_log10 = qPCR_bact_log10 + Centre_qPCR_bact_log10, CV_log10 = CV_log10 + Centre_CV_log10, inhibRichFrac = inhibRichFrac_raw) %>%
  filter(Richness>0) %>%
  select(Bd_micro_log10, Bd_qPCR_log10, Inhibitory, Richness, qPCR_bact_log10, CV_log10,inhibRichFrac) %>%
  pivot_longer(starts_with("Bd_"), names_to = "Response", values_to = "Y") %>%
  pivot_longer(c(Richness, qPCR_bact_log10, CV_log10), names_to = "Predictors", values_to = "X") %>%
  mutate(Inhibitory_designation = ifelse(inhibRichFrac > 0, "At least one\ninhibitory ASV", "Non-inhibitory ASVs")) %>%
  mutate(Predictors = factor(Predictors, levels = c("Richness","qPCR_bact_log10", "CV_log10")))  %>%
  ggplot(aes(x=X, y = Y, col = Inhibitory_designation)) + 
  geom_point() +
  geom_smooth(method = "lm", se=FALSE) +
  facet_grid(Response ~ Predictors, scales = "free", switch = "both",labeller=labeller(Predictors=pred.labs, Response=resp.labs)) +
  xlab("Predictors") + ylab("Response metrics") + labs(col = "Biofilm type")
ggsave("03_figure_generation/raw_dat_nozeros_plot.png", height=5, width=8,
       gg_rawdat_nozoer)

### more subtle

gg_rawdat_adj <- dat %>% 
  mutate(Richness = Richness + Centre_Richness, qPCR_bact_log10 = qPCR_bact_log10 + Centre_qPCR_bact_log10, CV_log10 = CV_log10 + Centre_CV_log10, inhibRichFrac = inhibRichFrac_raw) %>%
  select(Bd_micro_log10, Bd_qPCR_log10, Inhibitory, Richness, qPCR_bact_log10, CV_log10,inhibRichFrac) %>%
  pivot_longer(starts_with("Bd_"), names_to = "Response", values_to = "Y") %>%
  pivot_longer(c(Richness, qPCR_bact_log10, CV_log10), names_to = "Predictors", values_to = "X") %>%
  mutate(Inhibitory_designation = ifelse(inhibRichFrac > 0, "At least one\ninhibitory ASV", "Non-inhibitory ASVs")) %>%
  # mutate(inhibRich_adj = ifelse(inhibRichFrac==0, NA, inhibRichFrac)) %>%
  mutate(Predictors = factor(Predictors, levels = c("Richness","qPCR_bact_log10", "CV_log10")))  %>%
  ggplot(aes(x=X, y = Y)) + 
  # geom_point(aes(fill=inhibRich_adj), col=rgb(0,0,0,0),pch=21) +
  geom_point(aes(col=inhibRichFrac)) +
  geom_smooth(aes(group=Inhibitory, col=Inhibitory),method = "lm", se=FALSE) +
  facet_grid(Response ~ Predictors, scales = "free", switch = "both",labeller=labeller(Predictors=pred.labs, Response=resp.labs)) +
  xlab("Predictors") + ylab("Response metrics") + labs(col="Fraction of ASVs\nthat were inhibitory")+
  # scale_fill_gradient(low="darkgrey",high="purple", na.value = "black")
scale_color_gradient(low="darkgrey",high="salmon", na.value = "black")
gg_rawdat_adj

ggsave("03_figure_generation/raw_dat_plot_alt.png", height=5, width=8,
       gg_rawdat_adj)


#### /------------ ASV-level plots ------------/ ####
#### ASV competition ####

gg_competition <- challenge_iso_summary %>%
  rowwise() %>% mutate(Isolate = paste0(g_sp, " (",IsolateID,")")) %>% ungroup() %>%
  mutate(Isolate = factor(Isolate, levels=unique(Isolate))) %>%
  ggplot() + geom_linerange(aes(x=Isolate, ymin=fractionBeat-seFractionBeat, ymax=fractionBeat+seFractionBeat)) +
  geom_point(aes(x=Isolate, y=fractionBeat, col=aveCV, pch=Inhibitory), cex=4) +
  theme_bw() +theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ylab("Fraction of pair-wise comparisons\n where isolate was more abundant") +
  labs(col="Biofilm forming\nability (96-well plate)\n(CV)", pch="Does isolate create\nBd zone of inhibition?") +
  scale_color_gradient(low="darkgrey", high="purple") 
gg_competition
ggsave(filename = "03_figure_generation/ASV_competition.png"
       ,height=6, width=6, gg_competition)

#### ASV effects on bd  ####
gg_alleffects_together <- isolate_indiv_effects_together%>%
  rowwise() %>% mutate(Isolate = paste0(g_sp, " (",IsolateID,")")) %>% ungroup() %>%
  arrange(median) %>% mutate(Isolate = factor(Isolate, levels=unique(Isolate))) %>%
  ggplot () + geom_point(aes(x=median, y=Isolate, pch=Inhibitory, col=PD<0.025), cex=3) +
  geom_segment(aes(x=lower95, xend=upper95, y=Isolate, yend=Isolate,lty=PD<0.025, col=PD<0.025)) +
  geom_vline(aes(xintercept=0), col="red") +
  xlab("Coefficient estimate") +
  scale_color_manual(values=c("darkgrey","darkred"), na.value="black") +
  scale_linetype_manual(values=c(3,1)) +
  theme_bw() +
  scale_shape_manual(values=c(19,17))+
  labs(pch="Isolate inhibited\nBd on plates") +
  facet_grid(.~Response, labeller = labeller(Response=c(Bd_micro = "Effect of isolate on\nBd microscopy counts", Bd_qpcr = "Effect of isolate on\nBd qPCR copies")))

gg_alleffects_together
ggsave(filename = "03_figure_generation/ASV_allASVstogether.png"
       ,height=5, width=7, gg_alleffects_together)

gg_individualAddEffects <- isolate_indiv_effects_summary%>%
  rowwise() %>% mutate(Isolate = paste0(g_sp, " (",IsolateID,")")) %>% ungroup() %>%
  arrange(median) %>% mutate(Isolate = factor(Isolate, levels=unique(Isolate))) %>%
  ggplot () + geom_point(aes(x=median, y=Isolate, col=Inhibitory)) +
  geom_segment(aes(x=lower95, xend=upper95, y=Isolate, yend=Isolate, col=Inhibitory, lty=PD<0.025)) +
  geom_vline(aes(xintercept=0), col="red") +
  xlab("Effect of isolate on Bd microscopy counts") +
  scale_color_manual(values=c("darkgrey","salmon"), na.value="black") +
  scale_linetype_manual(values=c(2,1)) +theme_bw() +
  labs(col="Isolate inhibited\nBd on plates")
gg_individualAddEffects
ggsave(filename = "03_figure_generation/ASV_allASVswithotherpredictors.png"
       ,height=5, width=6, gg_individualAddEffects)

gg_residualonly <- isolate_indiv_effects_res %>%
  rowwise() %>% mutate(Isolate = paste0(g_sp, " (",IsolateID,")")) %>% ungroup() %>%
  arrang(median) %>% mutate(Isolate = factor(Isolate, levels=unique(Isolate))) %>%
  ggplot () + geom_point(aes(x=median, y=Isolate)) +
  geom_segment(aes(x=lower95, xend=upper95, y=Isolate, yend=Isolate)) +
  geom_vline(aes(xintercept=0), col="red") +
  xlab("Effect of isolate on Bd microscopy counts") +theme_bw() 
gg_residualonly
ggsave(filename = "03_figure_generation/ASV_allASVs_onresiduals.png"
       ,height=5, width=6, gg_residualonly)

#### ASV characteristics #####
# Get colnames
colnames(original_dat_withASV) <- gsub("^X","",colnames(original_dat_withASV))
# make filtered isoalte table
allIsolates <- colnames(original_dat_withASV[,which(colnames(original_dat_withASV) %in% isolate_info$IsolateID)])
# Add in isolates that were retained vs not
isolateInfo_long <- isolate_indiv_effects_summary %>% mutate(PRESENT=TRUE) %>% select(IsolateID, PRESENT) %>%
  full_join(isolate_info %>% filter(IsolateID %in% allIsolates)) %>% 
  mutate(PRESENT = ifelse(is.na(PRESENT), 0, 1), medGrowthMidpoint = as.numeric(medGrowthMidpoint)) %>%
  rowwise() %>% mutate(Isolate = paste0(g_sp, " (",IsolateID,")")) %>% ungroup() %>%
  # rename(`Biofilm thickness\n(CV)` = aveCV, `Time to reach\n1/2 carrying capacity (h)` = medGrowthMidpoint, `Persisted in experiment` = PRESENT, `Diameter of inhibition\n(mm)`=aveInhibZonemm) %>%
  select(-c(g_sp, BiofilmFormer)) %>%
  pivot_longer(-c(Isolate,IsolateID, Inhibitory), names_to = "Isolate trait", values_to="Value") %>%
  arrange(Inhibitory, Isolate) %>%
  mutate(Isolate = factor(Isolate, levels=unique(Isolate)), IsolateID = factor(IsolateID, levels=unique(IsolateID))) %>%
  select(-Inhibitory)


gg_present <- isolateInfo_long %>%
  filter(`Isolate trait` == "PRESENT") %>%
  mutate(Value=factor(ifelse(Value==1, "Yes","No"),levels=c("No","Yes")))%>%
  ggplot() + geom_point(aes(x="", y=Isolate, pch=factor(Value)), cex=4, show.legend = FALSE) +
  xlab("Detected\nin experiment")+
  ylab("Full list of introduced isolates") +
  scale_shape_manual(values=c(4,19)) +
  theme(legend.position = c(0.8,0.8), axis.text.y = element_text(hjust=0)) +labs(pch="Isolate\npresent?")
gg_present


gg_growth <- isolateInfo_long %>%
  filter(`Isolate trait` == "medGrowthMidpoint") %>%
  ggplot() + geom_point(aes(x=Value, y=IsolateID)) +
  xlab("Time to reach\n1/2 carrying capacity(h)")+ylab("")
gg_growth

gg_cv <- isolateInfo_long %>%
  filter(`Isolate trait` == "aveCV") %>%
  ggplot() + geom_bar(aes(x=IsolateID, y=Value), stat="identity", fill="purple") +
  coord_flip() +
  ylab("Biofilm thickness\n(CV OD @ 590nm)") +xlab("") 
gg_cv

gg_zone <- isolateInfo_long %>%
  filter(`Isolate trait` == "aveInhibZonemm") %>%
  mutate(Value=ifelse(Value==0, NA,Value)) %>%
  ggplot() + geom_point(aes(y=IsolateID, x="Approximate zone size", cex=Value),col="salmon", show.legend = FALSE) +
  xlab("Relative diameter\n of inhibition zone") + ylab("") +labs(cex="Diameter of\nzone of inhibition\n(mm)")
gg_zone

gg_allisolate_info <- gridExtra::grid.arrange(gg_present, gg_growth, gg_cv, gg_zone, layout_matrix = rbind(c(1,1,1,2,2,3,3,4,4)))

ggsave(filename="03_figure_generation/isolate_summary.png"
       ,gg_allisolate_info, height=6, width=10)
 
### Add "Isolates added" ###
### In main, only include individual isolates, include full version in supplement 

setwd("..")

