#!bin/bash

library(rstanarm)
library(bayestestR)
library(lme4)
library(car)
library(ggbiplot)
library(tidyverse)


######### Statistical analyses ###########
setwd("04_analysis/")
dir.create("02_statistical_analyses")
### Load ###
dat <- read.delim("01_data_summaries/downstream/dat_scaled_final.txt")
dat_unscaled <- read.delim("01_data_summaries/downstream/dat_nonscaled_final.txt")
dat_full <- read.delim("01_data_summaries/downstream/dat_full.txt")
isolate_info <- read.delim("../03_datacleaning/05_final_seq_filtering/downstream/isolate_info_filtered.txt")

##### Co-correlations between predictor variables ####
cor_predictors <- matrix(nrow = 0, ncol = 6, dimnames = list(NULL, c("X","Y", "Correlation","t","df","pvalue")))
relavent_predictors <- c("CV_log10", "qPCR_bact_log10", "Richness", "inhibRichFrac_raw")
for ( a in 1:(length(relavent_predictors)-1)) {
  for (b in (a+1):length(relavent_predictors)) {
    # a=1
    # b=2
    A = relavent_predictors[a]
    B = relavent_predictors[b]
    corTemp <- cor.test(dat[,A],dat[,B])
    cor_predictors <- rbind(cor_predictors
                            , data.frame(X=A, Y=B, Correlation=signif(corTemp$estimate,3), t=signif(corTemp$statistic,3), df=corTemp$parameter, pvalue=signif(corTemp$p.value,3)))
  }
}
write.table(cor_predictors, file="02_statistical_analyses/predictor_correlations.txt", quote = FALSE, row.names = FALSE, sep="\t")

# Response correlation
sink(file = "02_statistical_analyses/response_correlations.txt")
cor.test(dat$Bd_micro_log10, dat$Bd_qPCR_log10)
sink()

# Follow-up CV vs Inhibitory
sink(file = "02_statistical_analyses/corr_CV_vs_inhib.txt")
t.test(dat$CV_log10~as.character(dat$Inhibitory))
sink() 

# Follow-up CV vs inhibProp
sink(file = "02_statistical_analyses/corr_CV_vs_inhibProp.txt")
summary(lm(CV_log10 ~ inhibProp, data=dat))
summary(lm(CV_log10 ~ inhibRichFrac_raw, data=dat))
sink() 

# Follow-up CV vs inhibProp
corr_samps <- rstan::extract(stan_glm(CV_log10 ~ inhibRichFrac_raw, data=dat, iter=5000)$stanfit)
sink(file = "02_statistical_analyses/corr_CV_vs_inhibfrac_bayes.txt")
"Median: "
median(corr_samps$beta)
"\nCI: "
c(ci(corr_samps$beta,method="HDI")$CI_low, ci(corr_samps$beta,method="HDI")$CI_high)
"\nPD: "
min(c(sum(corr_samps$beta>0), sum(corr_samps$beta<0)))/length(corr_samps$beta)
sink() 

#### Variance Inflation Factor ####
lmer_model_qpcr <- lmer(Bd_qPCR_log10 ~ inhibRichFrac_raw + inhibRichFrac_raw:Richness + Richness + CV_log10 + qPCR_bact_log10 + CV_log10:inhibRichFrac_raw + qPCR_bact_log10:inhibRichFrac_raw+ (1|DateStart), data = dat)
lmer_model_micro <- lmer(Bd_micro_log10 ~ inhibRichFrac_raw + inhibRichFrac_raw:Richness + Richness + CV_log10 + qPCR_bact_log10 + CV_log10:inhibRichFrac_raw + qPCR_bact_log10:inhibRichFrac_raw+(1|DateStart), data = dat)
lmer_model_qpcrbymicro <- lmer(Bd_qPCR_log10 ~ Bd_micro_log10 + inhibRichFrac_raw + inhibRichFrac_raw:Richness + Richness + CV_log10 + qPCR_bact_log10 + CV_log10:inhibRichFrac_raw + qPCR_bact_log10:inhibRichFrac_raw+(1|DateStart), data = dat)
vif_all <- data.frame(Bd_qpcr_model = vif(lmer_model_qpcr)) %>% rownames_to_column(var="Predictor") %>%
  full_join(data.frame(Bd_micro_model = vif(lmer_model_micro))%>% rownames_to_column(var="Predictor")) %>%
  full_join(data.frame(Bd_qpcrbymicro_model = vif(lmer_model_qpcrbymicro))%>% rownames_to_column(var="Predictor"))
vif_all

# 
# lmer_model_qpcr <- lmer(Bd_qPCR_log10 ~ inhibRichFrac_raw + inhibRichFrac_raw:Richness + Richness + CV_log10 + qPCR_bact_log10 + CV_log10:inhibProp_raw + qPCR_bact_log10:inhibProp_raw+ (1|DateStart), data = dat)
# lmer_model_micro <- lmer(Bd_micro_log10 ~ inhibRichFrac_raw + inhibRichFrac_raw:Richness + Richness + CV_log10 + qPCR_bact_log10 + CV_log10:inhibProp_raw + qPCR_bact_log10:inhibProp_raw+(1|DateStart), data = dat)
# lmer_model_qpcrbymicro <- lmer(Bd_qPCR_log10 ~ Bd_micro_log10 + inhibRichFrac_raw + inhibRichFrac_raw:Richness + Richness + CV_log10 + qPCR_bact_log10 + CV_log10:inhibProp_raw + qPCR_bact_log10:inhibProp_raw+(1|DateStart), data = dat)
# vif_all <- data.frame(Bd_qpcr_model = vif(lmer_model_qpcr)) %>% rownames_to_column(var="Predictor") %>%
#   full_join(data.frame(Bd_micro_model = vif(lmer_model_micro))%>% rownames_to_column(var="Predictor")) %>%
#   full_join(data.frame(Bd_qpcrbymicro_model = vif(lmer_model_qpcrbymicro))%>% rownames_to_column(var="Predictor"))
# vif_all
write.table(vif_all, sep="\t", quote=FALSE, row.names=FALSE, file="02_statistical_analyses/VIF_table_predictors.txt")
sink("02_statistical_analyses/VIF_of_scaled_predictors.txt")
cat("############ Variance Inflation Factor of scaled predictors\n")
vif_all
cat("\n\n")

cat("############ qPCR model\n")
summary(lmer_model_qpcr)
cat("\n\n")
cat("############ Micro model\n")
summary(lmer_model_micro)
sink()

#### Bayesian models- Without inhibitory ####
dir.create("02_statistical_analyses/rstanarm_samps")

if (!file.exists("02_statistical_analyses/rstanarm_samps/allSamps_noinhib.RData")) {
  allSamps_noinhib <- list()
  for ( response in c("Bd_micro_log10", "Bd_qPCR_log10")) {
    frml <- paste0(response, "~ Richness + qPCR_bact_log10 + CV_log10 + (1|DateStart)")
    stanlmer_temp <- stan_lmer(formula(frml)
                               , data = dat
                               , prior = normal(0,5,autoscale=TRUE)
                               # , prior_intercept = exponential(autoscale=TRUE)
                               , adapt_delta = 0.99
                               , iter = 5000
    )
    ## Extract samples
    stanlmer_samps_temp <- rstan::extract(stanlmer_temp$stanfit)
    colnames(stanlmer_samps_temp$beta) <- names(fixef(stanlmer_temp))[-1]
    colnames(stanlmer_samps_temp$alpha) <- "Intercept"
    colnames(stanlmer_samps_temp$b) <- c(rownames(ranef(stanlmer_temp)$DateStart),"OTHER")
    allSamps_noinhib[[paste0(response)]] <- stanlmer_samps_temp
  }
  save(allSamps_noinhib, file = "02_statistical_analyses/rstanarm_samps/allSamps_noinhib.RData")
} else {
  load("02_statistical_analyses/rstanarm_samps/allSamps_noinhib.RData")
}


allSamps_noinhib_adj <- data.frame()
summary_bayes_noinhib <- data.frame()
for ( r in c("Bd_micro_log10", "Bd_qPCR_log10") ) {
  sampsTemp <- allSamps_noinhib[[r]]
  adj_samps <- data.frame(sampsTemp$alpha, sampsTemp$beta) %>% 
    cbind(sampsTemp$b[,1:6]) %>% as_tibble() %>%
    pivot_longer(everything(), names_to = "Coefficient", values_to = "Estimate") %>%
    # separate(Coefficient, into = c("Predictor", "GroupOrType"), sep = "X", remove=FALSE) %>%
    # mutate(GroupOrType = ifelse(is.na(GroupOrType), Predictor, GroupOrType), Predictor = ifelse(GroupOrType==Predictor, "Set", Predictor ),  response=r) %>%
    mutate(GroupOrType = ifelse(!Coefficient %in% c("Intercept", "Richness","CV_log10","qPCR_bact_log10"), Coefficient, "Predictor")) %>%
    mutate(Predictor = ifelse(Coefficient %in% c("Intercept", "Richness","CV_log10","qPCR_bact_log10"), Coefficient, "Set") ) %>%
    mutate(Predictor = factor(Predictor, levels = c("Intercept","Richness","CV_log10","qPCR_bact_log10","Set"))) %>%
    mutate(response = r) %>%
    select(Coefficient, Predictor, GroupOrType, Estimate, response)
  # Get centres for intercept
  # centres_temp <- dat %>% select(starts_with("Centre")) %>% distinct()
  # intercept_adj <- centres_temp[, colnames(centres_temp)[grep(r, colnames(centres_temp))]]
  summary_samps <- adj_samps %>%
    # mutate(intercept_adj = intercept_adj) %>%
    # mutate(Estimate = ifelse(Predictor == "Intercept", Estimate + intercept_adj, Estimate)) %>% select(-intercept_adj) %>%
    group_by(Coefficient, Predictor, GroupOrType, response) %>% 
    summarize(Median = signif(median(Estimate),3),lower95twoside=signif(ci(Estimate, method="HDI")$CI_low,3), lower95oneside = signif(quantile(Estimate, 0.05),3), upper95twoside = signif(ci(Estimate, method="HDI")$CI_high,3), upper95oneside = signif(quantile(Estimate, 0.95)), PD = signif(min(c(sum(Estimate>0), sum(Estimate<0)))/n(),3),minP = 1/n()) %>% 
    ungroup() %>%
    select(Coefficient, Predictor, GroupOrType, response, Median, lower95twoside, upper95twoside, lower95oneside, upper95oneside, PD)
  
  allSamps_noinhib_adj <- rbind(allSamps_noinhib_adj , adj_samps)
  summary_bayes_noinhib <- rbind(summary_bayes_noinhib, summary_samps)
}

## Look at histograms
# Random effects
allSamps_noinhib_adj %>% filter(Predictor == "Set") %>%
  ggplot(aes(x=Estimate, fill = GroupOrType)) +
  geom_histogram(bins=100, col = "black") + facet_grid(response ~ GroupOrType) +
  geom_vline(aes(xintercept=0), col="red", lty=2)
# Look at individual histograms to check for good sampling 
allSamps_noinhib_adj %>% 
  left_join(summary_bayes_noinhib) %>%
  filter(Predictor != "Set", response == "Bd_micro_log10") %>%
  mutate(GroupOrType = factor(GroupOrType, levels = c("Non","Inhib","InhibEff"))) %>%
  mutate(Significant = PD<0.025) %>%
  ggplot(aes(x=Estimate, fill = Significant)) +
  geom_histogram(bins=100) + facet_grid(. ~ Predictor) +
  geom_vline(aes(xintercept=0), col="red", lty=2) +
  scale_fill_manual(values=c("grey","red"))
allSamps_noinhib_adj %>% 
  left_join(summary_bayes_noinhib) %>%
  filter(Predictor != "Set", response == "Bd_qPCR_log10") %>%
  mutate(GroupOrType = factor(GroupOrType, levels = c("Non","Inhib","InhibEff"))) %>%
  mutate(Significant = PD<0.025) %>%
  ggplot(aes(x=Estimate, fill = Significant)) +
  geom_histogram(bins=100) + facet_grid(. ~ Predictor) +
  geom_vline(aes(xintercept=0), col="red", lty=2) +
  scale_fill_manual(values=c("grey","red"))

write.table(summary_bayes_noinhib, file="02_statistical_analyses/summary_bayes_mixedmodel_noinhib.txt", quote=FALSE, row.names = FALSE, sep = "\t")
write.table(allSamps_noinhib_adj, file="02_statistical_analyses/samplesAdjusted_bayes_mixedmodel_noInhib.txt", quote=FALSE, row.names = FALSE, sep = "\t")

#### Bayesian models- Basic ####
if (!file.exists("02_statistical_analyses/rstanarm_samps/allSamps.RData")) {
  allModel_basic <- list()
  allSamps <- list()
  for ( response in c("Bd_micro_log10", "Bd_qPCR_log10")) {
    # frml <- paste0(response, "~ Richness + qPCR_bact_log10 + CV_log10 +  inhibRich + qPCR_bact_log10_inhib + (1|DateStart)")
    frml <- paste0(response, "~ Richness + qPCR_bact_log10 + CV_log10  +inhibRichFrac_raw +  Richness:inhibRichFrac_raw + qPCR_bact_log10:inhibRichFrac_raw + CV_log10:inhibRichFrac_raw + (1|DateStart)")
    stanlmer_temp <- stan_lmer(formula(frml)
                               , data = dat
                               , prior = normal(0,5,autoscale=TRUE)
                               # , prior_intercept = exponential(autoscale=TRUE)
                               , adapt_delta = 0.99
                               , iter = 5000
    )
    plot(stanlmer_temp)
    ## Extract samples
    stanlmer_samps_temp <- rstan::extract(stanlmer_temp$stanfit)
    colnames(stanlmer_samps_temp$beta) <- names(fixef(stanlmer_temp))[-1]
    colnames(stanlmer_samps_temp$alpha) <- "Intercept"
    colnames(stanlmer_samps_temp$b) <- c(rownames(ranef(stanlmer_temp)$DateStart),"OTHER")
    allSamps[[paste0(response)]] <- stanlmer_samps_temp
    allModel_basic[[paste0(response)]] <- stanlmer_temp
  }
  save(allSamps, file = "02_statistical_analyses/rstanarm_samps/allSamps.RData")
  save(allModel_basic, file = "02_statistical_analyses/rstanarm_samps/allModel_basic.RData")
} else {
  load("02_statistical_analyses/rstanarm_samps/allSamps.RData")
  load("02_statistical_analyses/rstanarm_samps/allModel_basic.RData")
}
plot(allModel_basic$Bd_micro_log10)
plot(allModel_basic$Bd_qPCR_log10)

allSamps_adj <- data.frame()
summary_bayes <- data.frame()
for ( r in c("Bd_micro_log10", "Bd_qPCR_log10") ) {
  sampsTemp <- allSamps[[r]]
  adj_samps <- data.frame(sampsTemp$alpha, sampsTemp$beta) %>% rownames_to_column(var="iter") %>%
    pivot_longer(-iter, names_to = "Predictor", values_to = "Estimate") %>%
    mutate(Predictor = ifelse(Predictor == "inhibRichFrac_raw", "Intercept.inhibRichFrac_raw", Predictor)) %>%
    separate(Predictor, into = c("Main","Interaction"), sep = "\\.") %>%
    mutate(Interaction = ifelse(is.na(Interaction), "Non", "InhibEff")) %>%
    pivot_wider(names_from = Interaction, values_from = Estimate, values_fill = 0) %>%
    mutate(Inhib = Non + InhibEff) %>% select(iter, Main, Non,InhibEff, Inhib) %>%
    # mutate(Inhib = ifelse(Non==Inhib, NA, Inhib), InhibEff = ifelse(InhibEff==0, NA, InhibEff)) %>%
    pivot_longer(c(Non, InhibEff, Inhib), names_to = "Interaction", values_to = "Estimate") %>% filter(!is.na(Estimate)) %>%
    unite(Main, Interaction, col = "Predictor", sep = "X") %>%
    pivot_wider(names_from = Predictor, values_from = Estimate) %>%
    select(-c(iter))  %>%
    cbind(sampsTemp$b[,1:6]) %>% as_tibble() %>%
    pivot_longer(everything(), names_to = "Coefficient", values_to = "Estimate") %>%
    separate(Coefficient, into = c("Predictor", "GroupOrType"), sep = "X", remove=FALSE) %>% 
    mutate(GroupOrType = ifelse(is.na(GroupOrType), Predictor, GroupOrType), Predictor = ifelse(GroupOrType==Predictor, "Set", Predictor ),  response=r) %>%
    mutate(Predictor = factor(Predictor, levels = c("Intercept","Richness","CV_log10","qPCR_bact_log10","Set")))
  # Get centres for intercept
  # centres_temp <- dat %>% select(starts_with("Centre")) %>% distinct()
  # intercept_adj <- centres_temp[, colnames(centres_temp)[grep(r, colnames(centres_temp))]]
  summary_samps <- adj_samps %>%
    # mutate(intercept_adj = intercept_adj) %>%
    # mutate(Estimate = ifelse(Predictor == "Intercept", Estimate + intercept_adj, Estimate)) %>% select(-intercept_adj) %>%
    group_by(Coefficient, Predictor, GroupOrType, response) %>%
    summarize(Median = signif(median(Estimate),3),lower95twoside=signif(ci(Estimate, method="HDI")$CI_low,3), lower95oneside = signif(quantile(Estimate, 0.05),3), upper95twoside = signif(ci(Estimate, method="HDI")$CI_high,3), upper95oneside = signif(quantile(Estimate, 0.95)), PD = signif(min(c(sum(Estimate>0), sum(Estimate<0)))/n(),3), ,minP = 1/n()) %>%
    ungroup() %>%
    select(Coefficient, Predictor, GroupOrType, response, Median, lower95twoside, upper95twoside, lower95oneside, upper95oneside, PD)

    allSamps_adj <- rbind(allSamps_adj , adj_samps)
  summary_bayes <- rbind(summary_bayes, summary_samps)
}

## Look at histograms
# Random effects
allSamps_adj %>% filter(Predictor == "Set") %>%
  ggplot(aes(x=Estimate, fill = GroupOrType)) +
  geom_histogram(bins=100, col = "black") + facet_grid(response ~ GroupOrType) +
  geom_vline(aes(xintercept=0), col="red", lty=2)
# Look at individual histograms to check for good sampling 
allSamps_adj %>% 
  left_join(summary_bayes) %>%
  filter(Predictor != "Set", response == "Bd_micro_log10") %>%
  mutate(GroupOrType = factor(GroupOrType, levels = c("Non","Inhib","InhibEff"))) %>%
  mutate(Significant = PD<0.025) %>%
  ggplot(aes(x=Estimate, fill = Significant)) +
  geom_histogram(col = "black", bins=100) + facet_grid(GroupOrType ~ Predictor) +
  geom_vline(aes(xintercept=0), col="red", lty=2) +
  scale_fill_manual(values=c("grey","red"))
allSamps_adj %>% 
  left_join(summary_bayes) %>%
  filter(Predictor != "Set", response == "Bd_qPCR_log10") %>%
  mutate(GroupOrType = factor(GroupOrType, levels = c("Non","Inhib","InhibEff"))) %>%
  mutate(Significant = PD<0.025) %>%
  ggplot(aes(x=Estimate, fill = Significant)) +
  geom_histogram(col = "black", bins=100) + facet_grid(GroupOrType ~ Predictor) +
  geom_vline(aes(xintercept=0), col="red", lty=2) +
  scale_fill_manual(values=c("grey","red"))
# unique(allSamps_adj$response)

write.table(summary_bayes, file="02_statistical_analyses/summary_bayes_mixedmodel.txt", quote=FALSE, row.names = FALSE, sep = "\t")
write.table(allSamps_adj, file="02_statistical_analyses/samplesAdjusted_bayes_mixedmodel.txt", quote=FALSE, row.names = FALSE, sep = "\t")


##### Bayesian qPCR adjusted by microscopy- without inhibitory ####
if (!file.exists("02_statistical_analyses/rstanarm_samps/qpcr_by_micro_samps_noInhib.RData")) {
  stanlmer_qpcrbymicro_noInhib <- stan_lmer(Bd_qPCR_log10 ~ Bd_micro_log10 + Richness + qPCR_bact_log10 + CV_log10  + (1|DateStart)
                                    , data = dat
                                    , prior = normal(0,5,autoscale=TRUE)
                                    # , prior_intercept = exponential(autoscale=TRUE)
                                    , adapt_delta = 0.99
                                    , iter = 5000
  )
  ## Extract samples
  stanlmer_samps_qpcrbymicro_noInhib <- rstan::extract(stanlmer_qpcrbymicro_noInhib$stanfit)
  colnames(stanlmer_samps_qpcrbymicro_noInhib$beta) <- names(fixef(stanlmer_qpcrbymicro_noInhib))[-1]
  colnames(stanlmer_samps_qpcrbymicro_noInhib$alpha) <- "Intercept"
  colnames(stanlmer_samps_qpcrbymicro_noInhib$b) <- c(rownames(ranef(stanlmer_qpcrbymicro_noInhib)$DateStart),"OTHER")
  save(stanlmer_samps_qpcrbymicro_noInhib, file = "02_statistical_analyses/rstanarm_samps/qpcr_by_micro_samps.RData")
} else {
  load("02_statistical_analyses/rstanarm_samps/qpcr_by_micro_samps.RData")
}

adj_samps_qpcrbymicro_noInhib <- data.frame(stanlmer_samps_qpcrbymicro_noInhib$alpha, stanlmer_samps_qpcrbymicro_noInhib$beta) %>% 
  cbind(sampsTemp$b[,1:6]) %>% as_tibble() %>%
  pivot_longer(everything(), names_to = "Coefficient", values_to = "Estimate") %>%
  mutate(GroupOrType = ifelse(Coefficient %in% c("Intercept","Bd_micro_log10","Richness","qPCR_bact_log10", "CV_log10"), "Predictor", Coefficient)) %>%
  mutate(Predictor = ifelse(Coefficient %in% c("Intercept","Bd_micro_log10","Richness","qPCR_bact_log10", "CV_log10"), Coefficient, "Set")) %>%
  mutate(response = "qPCRbyMicro") %>%
  mutate(Predictor = factor(Predictor, levels = c("Intercept","Bd_micro_log10", "Richness","CV_log10","qPCR_bact_log10","Set")))

# Get centres for intercept
# centres_temp <- dat %>% select(starts_with("Centre")) %>% distinct()
# intercept_adj <- centres_temp[, colnames(centres_temp)[grep(r, colnames(centres_temp))]]
summary_samps_qpcrbymicro_noInhib <- adj_samps_qpcrbymicro_noInhib %>%
  # mutate(intercept_adj = intercept_adj) %>%
  # mutate(Estimate = ifelse(Predictor == "Intercept", Estimate + intercept_adj, Estimate)) %>% select(-intercept_adj) %>%
  group_by(Coefficient, Predictor, GroupOrType, response) %>% 
  summarize(Median = signif(median(Estimate),3),lower95twoside=signif(ci(Estimate, method="HDI")$CI_low,3), lower95oneside = signif(quantile(Estimate, 0.05),3), upper95twoside = signif(ci(Estimate, method="HDI")$CI_high,3), upper95oneside = signif(quantile(Estimate, 0.95)), PD = signif(min(c(sum(Estimate>0),sum(Estimate<0)))/n(),3), minP = 1/n()) %>% 
  ungroup() %>%
  select(Coefficient, Predictor, GroupOrType, response, Median, lower95twoside, upper95twoside, lower95oneside, upper95oneside, PD)

## Plot
# Random effects
adj_samps_qpcrbymicro_noInhib %>% filter(Predictor == "Set") %>%
  ggplot(aes(x=Estimate, fill = GroupOrType)) +
  geom_histogram(bins=100, col = "black") + facet_grid(response ~ GroupOrType) +
  geom_vline(aes(xintercept=0), col="red", lty=2)
# Look at individual histograms to check for good sampling 
adj_samps_qpcrbymicro_noInhib %>% 
  left_join(summary_samps_qpcrbymicro_noInhib) %>%
  filter(Predictor != "Set") %>%
  # mutate(GroupOrType = factor(GroupOrType, levels = c("Non","Inhib","InhibEff"))) %>%
  mutate(Significant = PD<0.025) %>%
  ggplot(aes(x=Estimate, fill = Significant)) +
  geom_histogram(col = "black", bins=100) + facet_grid(GroupOrType ~ Predictor) +
  geom_vline(aes(xintercept=0), col="red", lty=2) +
  scale_fill_manual(values=c("grey","red"))

write.table(summary_samps_qpcrbymicro_noInhib, file = "02_statistical_analyses/summary_bayes_pcrbymicro_noInhib.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(adj_samps_qpcrbymicro_noInhib, file = "02_statistical_analyses/samplesAdjust_bayes_qpcrbymicro_noInhib.txt", row.names = FALSE, quote = FALSE, sep = "\t")

stan_test <- stan_glm(Bd_qPCR_log10 ~ Bd_micro_log10, data=dat)
dat$Bd_residuals <- stan_test$residuals

##### Bayesian qPCR adjusted by microscopy- basic ####
if (!file.exists("02_statistical_analyses/rstanarm_samps/qpcr_by_micro_samps_basic.RData")) {
  # stanlmer_simple <- stan_glm(Bd_qPCR_log10 ~ Bd_micro_log10, data =dat)
  # sampsTemp <- rstan::extract(stanlmer_simple$stanfit)
  # median_prior <- median(sampsTemp$beta)
  # sd_prior <- sd(sampsTemp$beta)
  stanlmer_qpcrbymicro <- stan_lmer(Bd_qPCR_log10 ~ Bd_micro_log10 + Richness + qPCR_bact_log10 + CV_log10  + inhibRichFrac_raw + Richness:inhibRichFrac_raw + qPCR_bact_log10:inhibRichFrac_raw + CV_log10:inhibRichFrac_raw + (1|DateStart)
                            , data = dat
                             , prior = normal(0,5, autoscale = TRUE)
                             # , prior_intercept = exponential(autoscale=TRUE)
                             , adapt_delta = 0.99
                             , iter = 5000
  )
  # 
  # stan_test <- stan_glm(Bd_qPCR_log10 ~ Bd_micro_log10, data=dat)
  # dat$Bd_residuals <- stan_test$residuals
  # stanlmer_qpcrbymicro <- stan_lmer(Bd_residuals ~ Richness + qPCR_bact_log10 + CV_log10  + inhibRichFrac_raw + Richness:inhibRichFrac_raw + qPCR_bact_log10:inhibRichFrac_raw + CV_log10:inhibRichFrac_raw + (1|DateStart)
  #                                   , data = dat
  #                                   , prior = normal(0,5, autoscale = TRUE)
  #                                   # , prior_intercept = exponential(autoscale=TRUE)
  #                                   , adapt_delta = 0.99
  #                                   , iter = 5000
  # )
  # plot(stanlmer_qpcrbymicro, pars = names(fixef(stanlmer_qpcrbymicro)))
  ## Extract samples
  stanlmer_samps_qpcrbymicro <- rstan::extract(stanlmer_qpcrbymicro$stanfit)
  colnames(stanlmer_samps_qpcrbymicro$beta) <- names(fixef(stanlmer_qpcrbymicro))[-1]
  colnames(stanlmer_samps_qpcrbymicro$alpha) <- "Intercept"
  colnames(stanlmer_samps_qpcrbymicro$b) <- c(rownames(ranef(stanlmer_qpcrbymicro)$DateStart),"OTHER")
  save(stanlmer_samps_qpcrbymicro, file = "02_statistical_analyses/rstanarm_samps/qpcr_by_micro_samps_basic.RData")
  save(stanlmer_qpcrbymicro, file = "02_statistical_analyses/rstanarm_samps/stanlmer_qpcrbymicro.RData")
} else {
  load("02_statistical_analyses/rstanarm_samps/qpcr_by_micro_samps_basic.RData")
  load("02_statistical_analyses/rstanarm_samps/stanlmer_qpcrbymicro.RData")
}
# sampsTemp <- rstan::extract(stanlmer_qpcrbymicro$stanfit)
adj_samps_qpcrbymicro <- data.frame(stanlmer_samps_qpcrbymicro$alpha, stanlmer_samps_qpcrbymicro$beta) %>% rownames_to_column(var="iter") %>%
  pivot_longer(-iter, names_to = "Predictor", values_to = "Estimate") %>%
  mutate(Predictor = ifelse(Predictor == "inhibRichFrac_raw", "Intercept.inhibRichFrac_raw", Predictor)) %>%
  separate(Predictor, into = c("Main","Interaction"), sep = "\\.") %>%
  mutate(Interaction = ifelse(is.na(Interaction), "Non", "InhibEff")) %>%
  pivot_wider(names_from = Interaction, values_from = Estimate, values_fill = 0) %>%
  mutate(Inhib = Non + InhibEff) %>% select(iter, Main, Non,InhibEff, Inhib) %>%
  mutate(Inhib = ifelse(Non==Inhib, NA, Inhib), InhibEff = ifelse(InhibEff==0, NA, InhibEff)) %>%
  pivot_longer(c(Non, InhibEff, Inhib), names_to = "Interaction", values_to = "Estimate") %>% filter(!is.na(Estimate)) %>%
  unite(Main, Interaction, col = "Predictor", sep = "X") %>%
  pivot_wider(names_from = Predictor, values_from = Estimate) %>%
  select(-c(iter))  %>%
  cbind(stanlmer_samps_qpcrbymicro$b[,1:6]) %>% as_tibble() %>%
  pivot_longer(everything(), names_to = "Coefficient", values_to = "Estimate") %>%
  separate(Coefficient, into = c("Predictor", "GroupOrType"), sep = "X", remove=FALSE) %>%
  mutate(GroupOrType = ifelse(is.na(GroupOrType), Predictor, GroupOrType), Predictor = ifelse(GroupOrType==Predictor, "Set", Predictor ),  response="qPCRbyMicro") %>%
  mutate(Predictor = factor(Predictor, levels = c("Intercept","Bd_micro_log10", "Richness","CV_log10","qPCR_bact_log10","Set")))

# Get centres for intercept
# centres_temp <- dat %>% select(starts_with("Centre")) %>% distinct()
# intercept_adj <- centres_temp[, colnames(centres_temp)[grep(r, colnames(centres_temp))]]
summary_samps_qpcrbymicro <- adj_samps_qpcrbymicro %>%
  # mutate(intercept_adj = intercept_adj) %>%
  # mutate(Estimate = ifelse(Predictor == "Intercept", Estimate + intercept_adj, Estimate)) %>% select(-intercept_adj) %>%
  group_by(Coefficient, Predictor, GroupOrType, response) %>% 
  summarize(Median = signif(median(Estimate),3),lower95twoside=signif(ci(Estimate, method="HDI")$CI_low,3), lower95oneside = signif(quantile(Estimate, 0.05),3), upper95twoside = signif(ci(Estimate, method="HDI")$CI_high,3), upper95oneside = signif(quantile(Estimate, 0.95)), PD = signif(min(c(sum(Estimate>0),sum(Estimate<0)))/n(),3), minP = 1/n()) %>% 
  ungroup() %>%
  select(Coefficient, Predictor, GroupOrType, response, Median, lower95twoside, upper95twoside, lower95oneside, upper95oneside, PD)

## Plot
# Random effects
adj_samps_qpcrbymicro %>% filter(Predictor == "Set") %>%
  ggplot(aes(x=Estimate, fill = GroupOrType)) +
  geom_histogram(bins=100, col = "black") + facet_grid(response ~ GroupOrType) +
  geom_vline(aes(xintercept=0), col="red", lty=2)
# Look at individual histograms to check for good sampling 
adj_samps_qpcrbymicro %>% 
  left_join(summary_samps_qpcrbymicro) %>%
  filter(Predictor != "Set") %>%
  mutate(GroupOrType = factor(GroupOrType, levels = c("Non","Inhib","InhibEff"))) %>%
  mutate(Significant = PD<0.025) %>%
  ggplot(aes(x=Estimate, fill = Significant, col=Significant)) +
  geom_histogram(bins=100) + facet_grid(GroupOrType ~ Predictor) +
  geom_vline(aes(xintercept=0), col="red", lty=2) +
  scale_fill_manual(values=c("grey","red"))+
  scale_color_manual(values=c("grey","red"))

# summary_samps_qpcrbymicro <- summary_samps_qpcrbymicro %>% mutate(PosteriorProbability = 1-pvalue)

write.table(summary_samps_qpcrbymicro, file = "02_statistical_analyses/summary_bayes_pcrbymicro.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(adj_samps_qpcrbymicro, file = "02_statistical_analyses/samplesAdjust_bayes_qpcrbymicro.txt", row.names = FALSE, quote = FALSE, sep = "\t")


#### @~~~~~~~~~~~~~ASV-level effects~~~~~~~~~~~~~@ ####
#### PCA of ASVs ####
# Look at individual ASVs
dat_iso <- dat_full %>% select(starts_with("ISO"))
dat_iso_log <- log10(dat_iso+1)
dat_iso_log_bd <- cbind(dat_iso_log, dat_full[,c("Bd_micro_log10", "Bd_qPCR_log10")])
pca_iso <- prcomp(dat_iso_log_bd)
ggbiplot(pca_iso)


### Summing isolate qualities; compare ####

#### Predicting CV #####
### See if proportion within sample lines up with median growth and ave CV
allIso <- dat_full %>% select(starts_with("ISO")) %>% colnames()
### If you sum up individual biofilm formers, do you get what you expect?
dat_iso_summary <- cbind(SampleID=dat_full$SampleID, dat_iso) %>%
  mutate(maxinhib=NA, suminhib=NA, adjinhib=NA, adjinhib_raw = NA, maxCV=NA, sumCV=NA, adjCV=NA, adjCV_raw=NA, domIsoCV=NA, aveMidpoint=NA, maxMidpoint=NA)
for ( r in 1:nrow(dat_iso_summary)) {
  # r <- 1
  currentIso <- dat_iso_summary[r,] %>% select(starts_with("ISO")) %>%
    pivot_longer(c(allIso),names_to="Isolate", values_to="Reads") %>% 
    mutate(IsolateID=gsub("ISO_","",Isolate)) %>%
    filter(Reads>0) %>%
    mutate(PropAbund = Reads/sum(Reads)) 
  if ( nrow(currentIso)==0) {
    tempSummary <- c(0,0,0,0,0,0,0,0,0,0,0)
  } else {
    tempSummary <- currentIso %>% left_join(isolate_info) %>%
      mutate(maxReads = max(Reads)) %>%
      rowwise() %>%
      mutate(propinhib=aveInhibZonemm*PropAbund, propCV = aveCV*PropAbund
             , rawinhib = aveInhibZonemm*(Reads), rawCV=aveCV*(Reads)
             , domIso = Reads==maxReads
             ) %>% 
      group_by() %>%
      summarise(maxinhib=max(aveInhibZonemm), suminhib=sum(aveInhibZonemm), adjinhib=sum(propinhib), adjinhib_raw = sum(rawinhib)
                , maxCV = max(aveCV), sumCV = sum(aveCV, na.rm=TRUE), adjCV = sum(propCV), adjCV_raw = sum(rawCV), domIsoCV = sum(aveCV*domIso), aveMidpoint = mean(medGrowthMidpoint), maxMidpoint = max(medGrowthMidpoint))
    
  }

  dat_iso_summary[r,c("maxinhib","suminhib","adjinhib","adjinhib_raw","maxCV","sumCV","adjCV", "adjCV_raw","domIsoCV", "aveMidpoint", "maxMidpoint")] <- tempSummary
  
}


### Maximum CV predicts total CV the best-- thickness is saturating
dat_full_with_iso_summary <- dat_iso_summary %>% full_join(dat_full) %>%
  mutate(CV = 10^CV_log10-0.01) 
write.table(dat_full_with_iso_summary, file="02_statistical_analyses/dat_full_with_iso_summary.txt", quote=FALSE, row.names = FALSE, sep="\t")

dat_CV_assessment <- dat_full_with_iso_summary %>%
  filter(maxCV!=0) %>%
  rowwise() %>% mutate(propDom = max(across(starts_with("ISO")))/sum(across(starts_with("ISO")))) %>% ungroup() %>%
  mutate(propDom = ifelse(Richness==1, NA, propDom)) %>%
  # select(CV, sumCV, maxCV, adjCV, propDom,domIsoCV, Richness) %>% 
  select(CV, maxCV, adjCV, sumCV, propDom, Richness) %>%
  pivot_longer(-c(propDom, Richness, CV), values_to = "CV_predict", names_to = "CV_type")  %>%
  mutate(CV_type = ifelse(CV_type=="adjCV", "Weighted average of crystal violet\n across members",
                          ifelse(CV_type=="maxCV", "Maximum crystal violet intensity\n out of all members", 
                                 ifelse(CV_type=="sumCV", "Summed crystal violet intensity\nacross all members", 
                                        ifelse(CV_type=="adjCV_raw", "Read-weighted sum of crystal violet\nacross members",NA))))) %>%
  mutate(Rich1 = ifelse(Richness==1, "1", "Greater than 1")) %>%
  mutate(CV_type = factor(CV_type, levels=c("Weighted average of crystal violet\n across members", "Maximum crystal violet intensity\n out of all members", "Summed crystal violet intensity\nacross all members","Read-weighted sum of crystal violet\nacross members")))
dat_CV_assessment %>% ggplot() + 
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

#### Mean squared error from model
dat_CV_Rich1 <- dat_iso_summary %>%left_join(dat_full) %>% filter(Richness==1)
set.seed(2005)
stan_CV_R1 <- stan_lmer(CV_log10 ~ log10(maxCV+0.01) + (1|DateStart), data=dat_CV_Rich1)
samps_CV_R1 <- rstan::extract(stan_CV_R1$stanfit)
CV_R1_fit <- c(median(samps_CV_R1$alpha), median(samps_CV_R1$beta))

### Now calculate the predicted value from each
dat_CV_witherrors <- dat_CV_assessment %>% 
  mutate(CV_R1model = CV_R1_fit[1] + CV_R1_fit[2]*log10(CV_predict+0.01)) %>%
  mutate(Error = CV_R1model-log10(CV+0.01), SqE = Error^2) 
dat_CV_witherrors %>% 
  # filter(CV_type == "Summed crystal violet intensity\nacross all members") %>%
  ggplot() + geom_point(aes(x=log10(CV_predict+0.01), y=Error, col=CV_type)) +
  geom_hline(aes(yintercept = 0), lty=2, col="black")

PredictCV_stats <- dat_CV_witherrors %>% filter(Richness !=1) %>% select(Error, SqE, CV_type) %>%  
  group_by(CV_type) %>% summarise(MSE = mean(SqE), MBE = mean(Error), MAE = mean(abs(Error)), RMSE = sqrt(MSE))
write.table(PredictCV_stats, file="02_statistical_analyses/PredictCV_stats.txt", quote=FALSE, row.names = FALSE, sep="\t")  

#### Predicting Bd, Rich 1 ####
## Linear models
dat_singlesModel <-  dat_full_with_iso_summary %>% filter(Richness == 1) %>%
  select(Bd_micro_log10, Bd_qPCR_log10, CV_log10, qPCR_bact_log10, inhibRich, DateStart)
set.seed(23408)
# stan_singlesOnly_micro <- stan_lmer(Bd_micro_log10 ~ CV_log10 + qPCR_bact_log10 +inhibRich  + CV_log10:inhibRich+qPCR_bact_log10:inhibRich(1|DateStart), data=dat_singlesModel, iter = 5000)
# stan_singlesOnly_qpcr <- stan_lmer(Bd_qPCR_log10 ~ CV_log10 + qPCR_bact_log10 +inhibRich + CV_log10:inhibRich +qPCR_bact_log10:inhibRich+ (1|DateStart), data=dat_singlesModel, iter = 5000)
stan_singlesOnly_micro_CV <- stan_lmer(Bd_micro_log10 ~ CV_log10*inhibRich + (1|DateStart), data=dat_singlesModel, iter = 5000, adapt_delta = 0.99)
stan_singlesOnly_micro_qpcr <- stan_lmer(Bd_micro_log10 ~ qPCR_bact_log10*inhibRich+ (1|DateStart), data=dat_singlesModel, iter = 5000, adapt_delta = 0.99)
stan_singlesOnly_qpcr_CV <- stan_lmer(Bd_qPCR_log10 ~  CV_log10*inhibRich + (1|DateStart), data=dat_singlesModel, iter = 5000, adapt_delta = 0.99)
stan_singlesOnly_qpcr_qpcr <- stan_lmer(Bd_qPCR_log10 ~ qPCR_bact_log10*inhibRich+ (1|DateStart), data=dat_singlesModel, iter = 5000, adapt_delta = 0.99)
plot(stan_singlesOnly_micro_CV)
plot(stan_singlesOnly_micro_qpcr)
plot(stan_singlesOnly_qpcr_CV)
plot(stan_singlesOnly_qpcr_qpcr)
# Get samps
samps_singles_micro_CV <- rstan::extract(stan_singlesOnly_micro_CV$stanfit)
samps_singles_qPCR_CV <- rstan::extract(stan_singlesOnly_qpcr_CV$stanfit)
summary_singles_CV <- rbind(cbind(samps_singles_micro_CV$alpha, samps_singles_micro_CV$beta, response="micro"),cbind(samps_singles_qPCR_CV$alpha, samps_singles_qPCR_CV$beta, response="qpcr")) %>% as.data.frame() %>% 
  rename_at(vars(paste0("V",seq(1,4))), ~names(fixef(stan_singlesOnly_micro_CV))) %>%
  rownames_to_column(var="iter") %>%
  pivot_longer(-c(iter,response), values_to = "samp", names_to = "Predictor") %>%
  mutate(Predictor = ifelse(Predictor == "inhibRich", "(Intercept):inhibRich", Predictor)) %>%
  separate(Predictor, into = c("Main", "Inhibitory"), sep=":", remove=FALSE) %>%
  mutate(Inhibitory = ifelse(is.na(Inhibitory), "Non", "EffInhib"))  %>%
  select(-Predictor) %>%pivot_wider(names_from=Inhibitory, values_from=samp) %>%
  mutate(Non = as.numeric(Non), EffInhib = as.numeric(EffInhib)) %>%
  mutate(Inhib = as.numeric(Non) + as.numeric(EffInhib)) %>%
  pivot_longer(-c(iter,Main, response), values_to = "samps", names_to = "Inhibitory") %>%
  mutate(predictor="CV") %>%
  group_by(Main, Inhibitory, response, predictor) %>%
  summarize(median = median(samps), lower95 = ci(samps, method="HDI")$CI_low, upper95 = ci(samps, method="HDI")$CI_high, PD = min(c(sum(samps<0), sum(samps>0))/n())) %>%
  arrange(response, Main, Inhibitory) 

samps_singles_micro_qpcr <- rstan::extract(stan_singlesOnly_micro_qpcr$stanfit)
samps_singles_qPCR_qpcr <- rstan::extract(stan_singlesOnly_qpcr_qpcr$stanfit)
summary_singles_qpcr <- rbind(cbind(samps_singles_micro_qpcr$alpha, samps_singles_micro_qpcr$beta, response="micro"),cbind(samps_singles_qPCR_qpcr$alpha, samps_singles_qPCR_qpcr$beta, response="qpcr")) %>% as.data.frame() %>% 
  rename_at(vars(paste0("V",seq(1,4))), ~names(fixef(stan_singlesOnly_micro_qpcr))) %>%
  rownames_to_column(var="iter") %>%
  pivot_longer(-c(iter,response), values_to = "samp", names_to = "Predictor") %>%
  mutate(Predictor = ifelse(Predictor == "inhibRich", "(Intercept):inhibRich", Predictor)) %>%
  separate(Predictor, into = c("Main", "Inhibitory"), sep=":", remove=FALSE) %>%
  mutate(Inhibitory = ifelse(is.na(Inhibitory), "Non", "EffInhib"))  %>%
  select(-Predictor) %>%pivot_wider(names_from=Inhibitory, values_from=samp) %>%
  mutate(Non = as.numeric(Non), EffInhib = as.numeric(EffInhib)) %>%
  mutate(Inhib = as.numeric(Non) + as.numeric(EffInhib)) %>%
  pivot_longer(-c(iter,Main, response), values_to = "samps", names_to = "Inhibitory") %>%
  mutate(predictor = "qpcr") %>%
  group_by(Main, Inhibitory, response, predictor) %>%
  summarize(median = median(samps), lower95 = ci(samps, method="HDI")$CI_low, upper95 = ci(samps, method="HDI")$CI_high, PD = min(c(sum(samps<0), sum(samps>0))/n())) %>%
  arrange(response, Main, Inhibitory)

write.table(summary_singles_CV, file="02_statistical_analyses/summary_bayes_singles_CV.txt", quote=FALSE, row.names = FALSE, sep="\t")
write.table(summary_singles_qpcr, file="02_statistical_analyses/summary_bayes_singles_qpcr.txt", quote=FALSE, row.names = FALSE, sep="\t")

slopePDs <- rbind(summary_singles_CV, summary_singles_qpcr) %>% ungroup() %>%
  filter(Main!="(Intercept)") %>%
  select(Inhibitory, response, predictor, PD) %>%
  mutate(PD = ifelse(PD<0.025, "PD<0.025", "PD>=0.025")) %>%
  distinct()

bestfit_models_singles <- rbind(summary_singles_CV, summary_singles_qpcr) %>% ungroup() %>%
  filter(Inhibitory !="EffInhib") %>% select(Inhibitory, Main, response, predictor, median) %>%
  left_join(slopePDs) %>%
  mutate(Main = ifelse(Main=="(Intercept)", "b","m")) %>%
  pivot_wider(names_from = Main, values_from = median) %>%
  rename(Isolate_designation = Inhibitory, Predictor = predictor, Bd = response) %>% 
  mutate(Isolate_designation = ifelse(Isolate_designation=="Non", "Not inhibitory", "Inhibitory")) %>%
  mutate(Bd = ifelse(Bd=="qpcr", "Bd_qPCR_log10","Bd_micro_log10")) %>%
  mutate(Predictor = ifelse(Predictor == "CV", "CV_log10","qPCR_bact_log10"))

## Plot
dat_singlesOnly <- dat_full_with_iso_summary %>% filter(Richness == 1) %>%
  select(Bd_micro_log10, Bd_qPCR_log10, CV_log10, qPCR_bact_log10, inhibRich) %>%
  pivot_longer(-c(Bd_micro_log10, Bd_qPCR_log10, inhibRich), values_to = "predictor_value", names_to = "Predictor") %>%
  pivot_longer(-c(predictor_value, Predictor, inhibRich), values_to = "response_value", names_to = "Bd") %>%
  mutate(inhibRich = ifelse(inhibRich==1, "Inhibitory", "Not inhibitory")) %>%
  rename(Isolate_designation = inhibRich) %>%
  left_join(bestfit_models_singles)
dat_singlesOnly %>% 
  mutate(Predictor = ifelse(Predictor == "CV_log10", "Biofilm thickness\n(log10 CV + 0.01)", "Biofilm cell density\n(log10 qPCR copies+1)")) %>%
  mutate(Bd = ifelse(Bd == "Bd_micro_log10", "Bd microscopy counts", "Bd qPCR copies")) %>%
  ggplot(aes(x= predictor_value, y=response_value, col=Isolate_designation)) + geom_point() +
  # geom_smooth(method="lm", se = FALSE) +
  geom_abline(aes(slope=m, intercept=b, col=Isolate_designation, lty=PD)) +
  facet_grid(Bd ~ Predictor, scales="free", switch="both") +xlab("Predictor value") + ylab("Bd metric")
save(dat_singlesOnly, file="02_statistical_analyses/dat_singlesOnly.RData")


#### Predicting cell density with growth rate ####
dat_full_with_iso_summary %>% 
  filter(Richness>0) %>%
  mutate(Rich1 = ifelse(Richness==1, "1", "Greater than 1")) %>%
  select(qPCR_bact_log10, Rich1, aveMidpoint, maxMidpoint) %>%
  pivot_longer(-c(qPCR_bact_log10, Rich1), names_to = "Metric", values_to = "Growthrate") %>%
  ggplot() + geom_point(aes(x=log(Growthrate), y=qPCR_bact_log10, col=Rich1)) +
  facet_grid(.~Metric)

### No interesting enough to report. 

#### Predicting Bd with inhibitory traits ####

dat_iso_summary %>% full_join(dat_full) %>% filter(Richness>0) %>% 
  mutate(Rich1 = ifelse(Richness==1, "1", "Greater than 1")) %>%
  select(Bd_micro_log10, Bd_qPCR_log10, maxinhib, suminhib, adjinhib, inhibRich, Rich1) %>%
  pivot_longer(-c(Bd_micro_log10, Bd_qPCR_log10, Rich1), values_to="Value", names_to="Metric") %>%
  pivot_longer(-c(Value, Metric, Rich1), values_to="Bd", names_to="Response") %>%
  mutate(Metric = ifelse(Metric=="maxinhib", "COMPLEMENTARY:\nMaximum inhibitory zone of\nmost effective inhibiting isolate", 
                         ifelse(Metric == "suminhib", "ADDITIVE:\nSum of inhibitory zones of\nall inhibiting isolates"
                                , ifelse(Metric == "adjinhib", "AVERAGE:\nAverage zone of inhibition\nacross inhibitory isolates", 
                                         ifelse(Metric == "inhibRich", "RICHNESS:\nNumber of inhibitory isolates", NA))))) %>%
  mutate(Metric = factor(Metric, levels=c("ADDITIVE:\nSum of inhibitory zones of\nall inhibiting isolates",
                                          "COMPLEMENTARY:\nMaximum inhibitory zone of\nmost effective inhibiting isolate", 
                                          "AVERAGE:\nAverage zone of inhibition\nacross inhibitory isolates", 
                                          "RICHNESS:\nNumber of inhibitory isolates"))) %>%
  mutate(Response = ifelse(Response=="Bd_micro_log10", "Bd microscopy counts (log10)", "Bd qPCR copies (log10"))  %>%
  ggplot() + geom_jitter(aes(x=Value, y=Bd, col=Rich1), cex=3, width=0.3, height=0) +
  geom_smooth(aes(x=Value, y=Bd), method="lm") +
  facet_grid(Response~Metric, scales = "free")


## Only inhibitory

dat_iso_summary %>% full_join(dat_full) %>% filter(Richness>0) %>% filter(inhibRich>1) %>%
  mutate(Rich1 = ifelse(Richness==1, "1", "Greater than 1")) %>%
  select(Bd_micro_log10, Bd_qPCR_log10, maxinhib, suminhib, adjinhib, inhibRich, Rich1) %>%
  pivot_longer(-c(Bd_micro_log10, Bd_qPCR_log10, Rich1), values_to="Value", names_to="Metric") %>%
  pivot_longer(-c(Value, Metric, Rich1), values_to="Bd", names_to="Response") %>%
  mutate(Metric = ifelse(Metric=="maxinhib", "COMPLEMENTARY:\nMaximum inhibitory zone of\nmost effective inhibiting isolate", 
                         ifelse(Metric == "suminhib", "ADDITIVE:\nSum of inhibitory zones of\nall inhibiting isolates"
                                , ifelse(Metric == "adjinhib", "AVERAGE:\nAverage zone of inhibition\nacross inhibitory isolates", 
                                         ifelse(Metric == "inhibRich", "RICHNESS:\nNumber of inhibitory isolates", NA))))) %>%
  mutate(Metric = factor(Metric, levels=c("ADDITIVE:\nSum of inhibitory zones of\nall inhibiting isolates",
                                          "COMPLEMENTARY:\nMaximum inhibitory zone of\nmost effective inhibiting isolate", 
                                          "AVERAGE:\nAverage zone of inhibition\nacross inhibitory isolates", 
                                          "RICHNESS:\nNumber of inhibitory isolates"))) %>%
  mutate(Response = ifelse(Response=="Bd_micro_log10", "Bd microscopy counts (log10)", "Bd qPCR copies (log10"))  %>%
  ggplot() + geom_jitter(aes(x=Value, y=Bd, col=Rich1), cex=3, width=0.3, height=0) +
  geom_smooth(aes(x=Value, y=Bd), method="lm") +
  facet_grid(Response~Metric, scales = "free")


dat_inhibRichOnly <- dat_iso_summary %>% full_join(dat_full) %>% filter(Richness>0, inhibRich>1)

stan_predBdmicro_suminhib <- stan_lmer(Bd_micro_log10 ~ suminhib + (1|DateStart), data=dat_inhibRichOnly, iter=5000, adapt_delta=0.99)
stan_predBdmicro_maxinhib <- stan_lmer(Bd_micro_log10 ~ maxinhib + (1|DateStart), data=dat_inhibRichOnly, iter=5000, adapt_delta=0.99)
stan_predBdmicro_adjinhib <- stan_lmer(Bd_micro_log10 ~ adjinhib + (1|DateStart), data=dat_inhibRichOnly, iter=5000, adapt_delta=0.99)

stan_predBdqpcr_suminhib <- stan_lmer(Bd_qPCR_log10 ~ suminhib + (1|DateStart), data=dat_inhibRichOnly, iter=5000, adapt_delta=0.99)
stan_predBdqpcr_maxinhib <- stan_lmer(Bd_qPCR_log10 ~ maxinhib + (1|DateStart), data=dat_inhibRichOnly, iter=5000, adapt_delta=0.99)
stan_predBdqpcr_adjinhib <- stan_lmer(Bd_qPCR_log10 ~ adjinhib + (1|DateStart), data=dat_inhibRichOnly, iter=5000, adapt_delta=0.99)

plot(stan_predBdmicro_suminhib)
plot(stan_predBdmicro_maxinhib)
plot(stan_predBdmicro_adjinhib)
plot(stan_predBdqpcr_suminhib)
plot(stan_predBdqpcr_maxinhib)
plot(stan_predBdqpcr_adjinhib)

###### WORKING using bayes to find fit of each?#####

dat_iso_summary %>% full_join(dat_full) %>% filter(Richness>0) %>%
  select(Bd_micro_log10, Bd_qPCR_log10, inhibRich, inhibRichFrac, inhibProp) %>%
  pivot_longer(-c(Bd_micro_log10, Bd_qPCR_log10), values_to="Value", names_to="Metric") %>%
  mutate(Metric = ifelse(Metric=="inhibProp", "Proportion of isolates\nthat are inhibitory", 
                         ifelse(Metric == "inhibRich", "Number of isolates\nthat are inhibitory"
                                , ifelse(Metric == "inhibRichFrac", "Proportion of community\nthat is inhibitory", NA)))) %>%
  mutate(Metric = factor(Metric, levels=c("Number of isolates\nthat are inhibitory", "Proportion of isolates\nthat are inhibitory", "Proportion of community\nthat is inhibitory"))) %>%
  ggplot() + geom_point(aes(x=Value, y=Bd_micro_log10), cex=5) +
  facet_grid(.~Metric, scales = "free")

dat_iso_summary %>% full_join(dat_full) %>% filter(Richness>0) %>%
  ggplot() + geom_point(aes(x=log10(adjinhib_raw), y=Bd_micro_log10, col=CV_log10), cex=5)
dat_iso_summary %>% full_join(dat_full) %>% filter(Richness>0) %>%
  ggplot() + geom_point(aes(x=log10(adjinhib_raw), y=Bd_qPCR_log10, col=CV_log10), cex=5)
dat_iso_summary %>% full_join(dat_full) %>% filter(Richness>0) %>%
  ggplot() + geom_point(aes(x=adjCV_raw, y=Bd_micro_log10, col=adjinhib_raw), cex=5)

dat_iso_summary %>% full_join(dat_full) %>% filter(Richness>0) %>%
  ggplot() + geom_point(aes(x=log10(adjinhib_raw), y=Bd_qPCR_log10, col=CV_log10), cex=5)
dat_iso_summary %>% full_join(dat_full) %>% filter(Richness>0) %>%
  ggplot() + geom_point(aes(x=adjCV_raw, y=Bd_qPCR_log10, col=log10(adjinhib_raw)), cex=5)

######## Individual isolates stan_lmer ######
allIso_filt <- allIso[which(allIso!="ISO_24D")]

if ( !file.exists("02_statistical_analyses/rstanarm_samps/indiv_iso_samps.RData") ) {
  isoEffectSamps <- data.frame(iter=1:10000)
  for ( iso in allIso_filt) {
    # iso <- "ISO_BTB_59"
    ### Include as predictor in model
    dat_temp <- dat_full %>% rename_at(vars(all_of(iso)), ~"CurrentIso") %>%
      mutate(Present = CurrentIso>0)%>% 
      select(SampleID, CurrentIso, Present) %>%
      left_join(dat)
    lmer_temp <- stan_lmer(Bd_micro_log10 ~ Richness + CV_log10 + qPCR_bact_log10 + inhibRichFrac_raw + inhibRichFrac_raw:Richness + CV_log10:inhibRichFrac_raw + qPCR_bact_log10*inhibRichFrac_raw +Present + (1|DateStart), data=dat_temp
                           , iter=5000
                           , adapt_delta = 0.99)
    currentSamps <- as.data.frame(rstan::extract(lmer_temp$stanfit)$beta)
    colnames(currentSamps) <- names(fixef(lmer_temp))[-1]
    tempSamps <- data.frame(temp=currentSamps$PresentTRUE)
    colnames(tempSamps) <- iso
    isoEffectSamps <- cbind(isoEffectSamps, tempSamps)
    # # Bd
    # dat_full %>% rename_at(vars(all_of(iso)), ~"CurrentIso") %>%
    #   mutate(Present = CurrentIso>0) %>%
    #   ggplot() + geom_jitter(aes(x=log10(CurrentIso+1), y=Bd_micro_log10))
    # dat_full %>% rename_at(vars(all_of(iso)), ~"CurrentIso") %>%
    #   mutate(Present = CurrentIso>0) %>%
    #   ggplot() + geom_jitter(aes(x=log10(CurrentIso+1), y=Bd_qPCR_log10))
    # # CV/qpcr
    # dat_full %>% rename_at(vars(all_of(iso)), ~"CurrentIso") %>%
    #   mutate(Present = CurrentIso>0) %>%
    #   ggplot() + geom_point(aes(x=Present, y=CV_log10))
    # dat_full %>% rename_at(vars(all_of(iso)), ~"CurrentIso") %>%
    #   mutate(Present = CurrentIso>0) %>%
    #   ggplot() + geom_point(aes(x=Present, y=qPCR_bact_log10))
    # 

  }
  save(isoEffectSamps, file = "02_statistical_analyses/rstanarm_samps/indiv_iso_samps.RData")
  write.table(isoEffectSamps, file="02_statistical_analyses/isoEffectSamps.txt", quote=FALSE, row.names = FALSE, sep="\t")
  
} else {
  load("02_statistical_analyses/rstanarm_samps/indiv_iso_samps.RData")
}

isolate_indiv_effects_summary <- isoEffectSamps %>% pivot_longer(one_of(allIso_filt), names_to="Isolate", values_to="Effect") %>%
  group_by(Isolate) %>%
  summarise(median = median(Effect), lower95 = ci(Effect, method="HDI")$CI_low, upper95 = ci(Effect, method="HDI")$CI_high
            , PD = min(c(sum(Effect>0), sum(Effect<0)))/n()) %>%
  ungroup() %>% 
  mutate(IsolateID = gsub("ISO_","",Isolate), sig = PD<0.05) %>%
  left_join(isolate_info) %>%
  arrange(median) %>%
  mutate(IsolateID = factor(IsolateID, levels=unique(IsolateID)))

isolate_indiv_effects_summary%>%
  ggplot () + geom_point(aes(x=median, y=IsolateID, col=Inhibitory)) +
  geom_segment(aes(x=lower95, xend=upper95, y=IsolateID, yend=IsolateID, col=Inhibitory, lty=PD<0.025)) +
  geom_vline(aes(xintercept=0), col="red") +
  xlab("Effect of isolate on Bd microscopy counts") +
  scale_color_manual(values=c("darkgrey","salmon"), na.value="black") +
  scale_linetype_manual(values=c(2,1)) +theme_bw() +
  labs(col="Isolate inhibited\nBd on plates")
save(isolate_indiv_effects_summary, file="02_statistical_analyses/isolate_indiv_effects_summary.RData")
write.table(isolate_indiv_effects_summary, file="02_statistical_analyses/isolate_indiv_effects_summary.txt", quote=FALSE, row.names = FALSE, sep="\t")

#### Model with EVERYTHING in it #####
dat_everything <- dat_full %>% select(SampleID, starts_with("ISO")) %>% select(-ISO_24D) %>%
  left_join(dat)

iso_frml <- paste(allIso_filt, collapse = "+")
base_frml <- paste0("Bd_micro_log10 ~",iso_frml," +Richness + CV_log10 + qPCR_bact_log10 + inhibRichFrac_raw + inhibRichFrac_raw:Richness + CV_log10:inhibRichFrac_raw + qPCR_bact_log10*inhibRichFrac_raw + (1|DateStart)")
set.seed(3492)
stan_everything <- stan_glmer(as.formula(base_frml),  data=dat_everything, iter=10000)

iso_everything <- as.data.frame(rstan::extract(stan_everything$stanfit)$beta)
colnames(iso_everything) <- names(fixef(stan_everything))[-1]
isolate_everything_summary <- iso_everything %>% pivot_longer(one_of(colnames(iso_everything)), names_to="Predictor", values_to="samp") %>%
  group_by(Predictor) %>%
  summarise(median = median(samp), lower95 = ci(samp, method="HDI")$CI_low, upper95 = ci(samp, method="HDI")$CI_high
            , PD = min(c(sum(samp>0), sum(samp<0)))/n()) %>%
  ungroup() %>% 
  mutate(IsolateID = gsub("ISO_","",Predictor), sig = PD<0.05) %>%
  mutate(IsolateID = ifelse(Predictor=="Richness", "Richness",IsolateID)) %>%
  left_join(isolate_info) %>%
  arrange(median) %>%
  mutate(IsolateID = factor(IsolateID, levels=unique(IsolateID)))

isolate_everything_summary%>%
  filter(IsolateID %in% gsub("ISO_","",allIso_filt)) %>%
  rowwise() %>% mutate(Isolate = paste0(g_sp, " (",IsolateID,")")) %>% ungroup() %>%
  mutate(Isolate = factor(Isolate, levels=unique(Isolate))) %>%
  ggplot () + geom_point(aes(x=median, y=Isolate, col=Inhibitory)) +
  geom_segment(aes(x=lower95, xend=upper95, y=Isolate, yend=Isolate, col=Inhibitory, lty=PD<0.025)) +
  geom_vline(aes(xintercept=0), col="red") +
  xlab("Effect of isolate on Bd microscopy counts") +
  scale_color_manual(values=c("darkgrey","salmon"), na.value="black") +
  scale_linetype_manual(values=c(2,1)) +theme_bw() +
  labs(col="Isolate inhibited\nBd on plates")

isolate_everything_summary%>%
  filter(!IsolateID %in% gsub("ISO_","",allIso_filt)) %>%
  ggplot () + geom_point(aes(x=median, y=Predictor, col=Inhibitory)) +
  geom_segment(aes(x=lower95, xend=upper95, y=Predictor, yend=Predictor, col=Inhibitory, lty=PD<0.025)) +
  geom_vline(aes(xintercept=0), col="red") +
  xlab("Effect of isolate on Bd microscopy counts") +
  scale_color_manual(values=c("darkgrey","salmon"), na.value="black") +
  scale_linetype_manual(values=c(2,1)) +theme_bw() +
  labs(col="Isolate inhibited\nBd on plates")

save(isolate_everything_summary, file="02_statistical_analyses/isolate_everything_summary.RData")
write.table(isolate_everything_summary, file="02_statistical_analyses/isolate_everything_summary.RData", quote=FALSE, row.names = FALSE, sep="\t")


#### Include a model with JUST isolates ####
dat_iso_log_micro <- dat_iso_log_bd %>% select(-Bd_qPCR_log10, -ISO_24D) #%>% mutate(Richness = dat_full$Richness) 
dat_iso_log_qPCR <- dat_iso_log_bd %>% select(-Bd_micro_log10, -ISO_24D) #%>% mutate(Richness = dat_full$Richness) 

## Micro
set.seed(3492)
stan_isoonly_micro <- stan_glm(Bd_micro_log10 ~ ., data=dat_iso_log_micro, iter=5000)
iso_effects_micro <- as.data.frame(rstan::extract(stan_isoonly_micro$stanfit)$beta)
colnames(iso_effects_micro) <- names(fixef(stan_isoonly_micro))[-1]
isolate_indiv_effects_together_micro <- iso_effects_micro %>% pivot_longer(one_of(colnames(iso_effects_micro)), names_to="Predictor", values_to="samp") %>%
  group_by(Predictor) %>%
  summarise(median = median(samp), lower95 = ci(samp, method="HDI")$CI_low, upper95 = ci(samp, method="HDI")$CI_high
            , PD = min(c(sum(samp>0), sum(samp<0)))/n()) %>%
  ungroup() %>% 
  mutate(IsolateID = gsub("ISO_","",Predictor), sig = PD<0.05) %>%
  mutate(IsolateID = ifelse(Predictor=="Richness", "Richness",IsolateID)) %>%
  left_join(isolate_info) %>%
  arrange(median) %>%
  mutate(IsolateID = factor(IsolateID, levels=unique(IsolateID)), Response = "Bd_micro")

## qPCR
set.seed(3492)
stan_isoonly_qpcr <- stan_glm(Bd_qPCR_log10 ~ ., data=dat_iso_log_qPCR, iter=5000)
iso_effects_qpcr <- as.data.frame(rstan::extract(stan_isoonly_qpcr$stanfit)$beta)
colnames(iso_effects_qpcr) <- names(fixef(stan_isoonly_qpcr))[-1]
isolate_indiv_effects_together_qpcr <- iso_effects_qpcr %>% pivot_longer(one_of(colnames(iso_effects_qpcr)), names_to="Predictor", values_to="samp") %>%
  group_by(Predictor) %>%
  summarise(median = median(samp), lower95 = ci(samp, method="HDI")$CI_low, upper95 = ci(samp, method="HDI")$CI_high
            , PD = min(c(sum(samp>0), sum(samp<0)))/n()) %>%
  ungroup() %>% 
  mutate(IsolateID = gsub("ISO_","",Predictor), sig = PD<0.05) %>%
  mutate(IsolateID = ifelse(Predictor=="Richness", "Richness",IsolateID)) %>%
  left_join(isolate_info) %>%
  arrange(median) %>%
  mutate(IsolateID = factor(IsolateID, levels=unique(IsolateID)), Response = "Bd_qpcr")

isolate_indiv_effects_together <- rbind(isolate_indiv_effects_together_micro, isolate_indiv_effects_together_qpcr)

isolate_indiv_effects_together%>%
  ggplot () + geom_point(aes(x=median, y=IsolateID, col=Inhibitory)) +
  geom_segment(aes(x=lower95, xend=upper95, y=IsolateID, yend=IsolateID, col=Inhibitory, lty=PD<0.025)) +
  geom_vline(aes(xintercept=0), col="red") +
  xlab("Coefficient estimate") +
  scale_color_manual(values=c("darkgrey","salmon"), na.value="black") +
  scale_linetype_manual(values=c(2,1)) +theme_bw() +
  labs(col="Isolate inhibited\nBd on plates") +
  # facet_grid(.~Response)
  facet_grid(.~Response, labeller = labeller(Response=c(Bd_micro = "Effect of isolate on\nBd microscopy counts", Bd_qpcr = "Effect of isolate on\nBd qPCR copies")))

  
save(isolate_indiv_effects_together, file="02_statistical_analyses/isolate_indiv_effects_together.RData")
write.table(isolate_indiv_effects_together, file="02_statistical_analyses/isolate_indiv_effects_together.txt", quote=FALSE, row.names = FALSE, sep="\t")

#### Presence/absence instead ####
# dat_iso_PA <- dat_iso
# dat_iso_PA[dat_iso_PA>0] <- 1
# dat_iso_log_micro_PA <- dat_iso_PA %>% select(-ISO_24D) %>% mutate(Bd_micro_log10 = dat_full$Bd_micro_log10) 
# setseed(236)
# stan_isoonly_PA <- stan_glm(Bd_micro_log10 ~ ., data=dat_iso_log_micro_PA, iter=5000)
# iso_effects_PA <- as.data.frame(rstan::extract(stan_isoonly_PA$stanfit)$beta)
# colnames(iso_effects_PA) <- names(fixef(stan_isoonly_PA))[-1]
# isolate_indiv_effects_PA_together <- iso_effects_PA %>% pivot_longer(one_of(colnames(iso_effects_PA)), names_to="Predictor", values_to="samp") %>%
#   group_by(Predictor) %>%
#   summarise(median = median(samp), lower95 = ci(samp, method="HDI")$CI_low, upper95 = ci(samp, method="HDI")$CI_high
#             , PD = min(c(sum(samp>0), sum(samp<0)))/n()) %>%
#   ungroup() %>% 
#   mutate(IsolateID = gsub("ISO_","",Predictor), sig = PD<0.05) %>%
#   mutate(IsolateID = ifelse(Predictor=="Richness", "Richness",IsolateID)) %>%
#   left_join(isolate_info) %>%
#   arrange(median) %>%
#   mutate(IsolateID = factor(IsolateID, levels=unique(IsolateID)))
# 
# isolate_indiv_effects_PA_together%>%
#   ggplot () + geom_point(aes(x=median, y=IsolateID, col=Inhibitory)) +
#   geom_segment(aes(x=lower95, xend=upper95, y=IsolateID, yend=IsolateID, col=Inhibitory, lty=PD<0.025)) +
#   geom_vline(aes(xintercept=0), col="red") +
#   xlab("Effect of isolate on Bd microscopy counts") +
#   scale_color_manual(values=c("darkgrey","salmon"), na.value="black") +
#   scale_linetype_manual(values=c(2,1)) +theme_bw() +
#   labs(col="Isolate inhibited\nBd on plates")


### Effect on residual instead ####
stan_micro_model <- allModel_basic[["Bd_micro_log10"]]
Bd_residuals <- data.frame(Bd_residuals=stan_micro_model$residuals, SampleID=stan_micro_model$data %>% select(SampleID))
dat_bdres_iso <- dat_full %>% select(SampleID, starts_with("ISO")) %>% select(-ISO_24D) %>%
  left_join(Bd_residuals) %>% select(-SampleID)
setseed(8765)
stan_isoonly_res <- stan_glm(Bd_residuals ~ ., data=dat_bdres_iso, iter=10000)
iso_effects_res <- as.data.frame(rstan::extract(stan_isoonly_res$stanfit)$beta)
colnames(iso_effects_res) <- names(fixef(stan_isoonly_res))[-1]
isolate_indiv_effects_res <- iso_effects_res %>% pivot_longer(one_of(colnames(iso_effects_res)), names_to="Predictor", values_to="samp") %>%
  group_by(Predictor) %>%
  summarise(median = median(samp), lower95 = ci(samp, method="HDI")$CI_low, upper95 = ci(samp, method="HDI")$CI_high
            , PD = min(c(sum(samp>0), sum(samp<0)))/n()) %>%
  ungroup() %>%
  mutate(IsolateID = gsub("ISO_","",Predictor), sig = PD<0.05) %>%
  mutate(IsolateID = ifelse(Predictor=="Richness", "Richness",IsolateID)) %>%
  left_join(isolate_info) %>%
  arrange(median) %>%
  mutate(IsolateID = factor(IsolateID, levels=unique(IsolateID)))

isolate_indiv_effects_res%>%
  ggplot () + geom_point(aes(x=median, y=IsolateID, col=aveCV)) +
  geom_segment(aes(x=lower95, xend=upper95, y=IsolateID, yend=IsolateID, col=aveCV, lty=PD<0.025)) +
  geom_vline(aes(xintercept=0), col="red") +
  xlab("Effect of isolate on Bd microscopy counts") +
  # scale_color_manual(values=c("darkgrey","salmon"), na.value="black") +
  scale_linetype_manual(values=c(2,1)) +theme_bw() #+
  # labs(col="Isolate inhibited\nBd on plates")

save(isolate_indiv_effects_res, file = "02_statistical_analyses/isolate_indiv_effects_res.RData")
write.table(isolate_indiv_effects_res, file = "02_statistical_analyses/isolate_indiv_effects_res.txt", quote=FALSE, row.names = FALSE, sep="\t")

#### Isolates outcompeting one another ####
challenge_mat <- matrix(ncol=length(allIso_filt), nrow=length(allIso_filt), dimnames = list(allIso_filt, allIso_filt))
for ( iso1 in allIso_filt) {
  for (iso2 in allIso_filt[which(allIso_filt!=iso1)]) {
    challenge_mat[iso1, iso2] <- dat_iso %>% select(one_of(c(iso1, iso2))) %>%
      filter(get(iso1)>0 & get(iso2)>0) %>%
      mutate(win=sign(get(iso1)-get(iso2))>0) %>%
      summarise(compete = sum(win)/n()) %>% as.numeric()
  }
}
challenge_iso_adj <- challenge_mat %>% as.data.frame() %>%rownames_to_column(var="Isolate_testing") %>%
  pivot_longer(-Isolate_testing, names_to="Isolate_comparing", values_to = "Result") %>%
  filter(!is.na(Result)) %>%
  group_by(Isolate_testing) %>%
  mutate(fraction = mean(Result, na.rm=TRUE)) %>% ungroup() %>%
  arrange(fraction) %>% mutate(Isolate_testing = factor(Isolate_testing, levels=unique(Isolate_testing))) %>%
  mutate(IsolateID = gsub("ISO_","",Isolate_testing)) %>%
  left_join(isolate_info) 
challenge_iso_summary <- challenge_iso_adj %>%
  group_by(Isolate_testing, IsolateID) %>%
  summarise(fractionBeat = mean(Result), seFractionBeat = sd(Result)/sqrt(n())) %>%
  left_join(isolate_info)
challenge_iso_summary %>%
  rowwise() %>% mutate(Isolate = paste0(g_sp, " (",IsolateID,")")) %>% ungroup() %>%
  mutate(Isolate = factor(Isolate, levels=unique(Isolate))) %>%
  ggplot() + geom_linerange(aes(x=Isolate, ymin=fractionBeat-seFractionBeat, ymax=fractionBeat+seFractionBeat)) +
  geom_point(aes(x=Isolate, y=fractionBeat, col=aveCV), cex=4) +
  theme_bw() +theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ylab("Fraction of pair-wise comparisons\n where isolate was more abundant") +
  labs(col="Biofilm forming\nability (CV)") +
  scale_color_gradient(low="darkgrey", high="purple") 

save(challenge_iso_summary, file="02_statistical_analyses/challenge_iso_summary.RData")
write.table(challenge_iso_summary, file="02_statistical_analyses/challenge_iso_summary.txt", quote=FALSE, row.names = FALSE, sep="\t")

#### Saturating effects of CV, qPCR
isolate_info
# get isolate PA table
iso_table <- dat_full %>% select(SampleID, starts_with("ISO"))
iso_table[,-1] <- ifelse(iso_table[,-1]>0,1,0)
iso_table_mat <- iso_table[,-1]
rownames(iso_table_mat) <- iso_table[,1]
# 
# iso_table_rev <- t(iso_table[,-1])
# colnames(iso_table_rev) <- iso_table[,1]

# get coefficients
iso_cv_coef <- isolate_info %>%
  rowwise() %>% mutate(Isolate = paste0("ISO_", IsolateID)) %>% ungroup() %>% 
  select(Isolate, aveCV)
iso_cv_mat <- iso_cv_coef[match(colnames(iso_table_mat), iso_cv_coef$Isolate),2]
iso_mat_CV <- apply(as.matrix(iso_table_mat), MARGIN=1, FUN=function(x) x*as.matrix(iso_cv_mat))
# Aggregate data
summary_CV_by_sample <- t(apply(iso_mat_CV, 2, function(x) c(maxCV=max(x), sumCV=sum(x)))) %>%
  as.data.frame() %>%rownames_to_column(var="SampleID")

dat_unscaled %>% left_join(summary_CV_by_sample) %>%
  mutate(CV = 10^CV_log10-0.001) %>%
  select(SampleID, CV, maxCV, sumCV) %>%
  pivot_longer(-c(SampleID, CV), names_to = "SummaryType", values_to = "CV_summary") %>%
  mutate(SummaryType = ifelse(SummaryType == "maxCV", "Maximum CV of\nall isolates in biofilm", "Sum of all individual\nisolate CV values")) %>%
  ggplot() + geom_point(aes(x=CV_summary, y=CV)) +
  facet_grid(.~SummaryType, scales="free") +
  xlab("Expected CV based on biofilm composition") + ylab("Observed CV in mixed biofilms")

## How about if we look at richness ==1?
rowSums(iso_table_mat) %>% as.data.frame()
dat_unscaled %>% filter(Richness==1)


dat_unscaled %>% left_join(summary_CV_by_sample) %>%
  ggplot() + geom_point(aes(x=sumCV, y=CV_log10))

#### Single isolate test


left_join(dat_everything) %>%
  ggplot() + geom_point(aes(x=CV_log10, y=Bd_micro_log10, col=factor(Inhibitory), pch=factor(RichLevel)))
dat_unscaled %>% #filter(Richness==1) %>%
  left_join(dat_everything) %>%
  ggplot() + geom_point(aes(x=CV_log10, y=Bd_micro_log10, col=factor(Inhibitory), pch=factor(RichLevel)))
##

# Testing
# dat_full %>%
#   ggplot() + geom_point(aes(x=RichLevel, y=Bd_micro_log10, col=CV_log10), cex=3) +
#   facet_grid(.~Inhibitory) +
#   scale_color_gradient(low='grey', high='purple') +
#   theme_bw()
# dat_full %>%
#   ggplot() + geom_point(aes(x=RichLevel, y=Bd_micro_log10, col=qPCR_bact_log10), cex=3) +
#   facet_grid(.~Inhibitory) +
#   scale_color_gradient(low='grey', high='purple') +
#   theme_bw()

setwd("..")

