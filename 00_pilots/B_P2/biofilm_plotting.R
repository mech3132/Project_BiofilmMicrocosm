#!bin/bash
library(tidyverse)
library(gridExtra)
library(car)
setwd("~/Documents/PhD/Project_biofilm/")
dat <- read.delim("biofilm_results.txt", header=TRUE, stringsAsFactors = FALSE)
dat <- dat %>%
  mutate(log_CFUfilter = log(CFU_estimate_filter), log_CFUwater = log(CFU_estimate_water)) %>%
  mutate(Media_source=factor(Media_source, levels=c("IN","OUT","CON")))
dat %>%
  ggplot(aes(x=Time_day, y=CV_adj, col=Media_source)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  scale_color_manual(values=c("red","blue","darkgrey"))

dat %>%
  ggplot(aes(x=Time_day, y=log_CFUwater, col=Media_source)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE)+
  scale_color_manual(values=c("red","blue","darkgrey"))
dat %>%
  ggplot(aes(x=Time_day, y=log_CFUfilter, col=Media_source)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE)+
  scale_color_manual(values=c("red","blue","darkgrey"))

g_legend <- function(a.gplot){ # from stack overflow
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

dat %>%
  filter(Richness!=0) %>%
  mutate(Time_day=factor(Time_day)) %>%
  ggplot() +
  geom_boxplot(aes(x=Time_day, y=log_CFUfilter, col=Media_source)) +
  scale_color_manual(values=c("red","blue"))

# Boxplots
pdf("biofilm_plots.pdf", height=4, width=4)
grid.arrange(
dat %>%
  mutate(Time_day=factor(Time_day)) %>%
  ggplot() +
  geom_boxplot(aes(x=Time_day, y=CV_adj, col=Media_source), show.legend=FALSE) +
  scale_color_manual(values=c("red","blue","darkgrey")) +
  ylab("CV Intensity") + xlab("Day")
, dat %>%
  mutate(Time_day = factor(Time_day)) %>%
  ggplot() +
  geom_boxplot(aes(x=Time_day, y=log_CFUfilter, col=Media_source), show.legend=FALSE)+
  scale_color_manual(values=c("red","blue","darkgrey")) +
  ylab("CFU (Filter)") + xlab("Day")
, g_legend(dat %>%
             mutate(Time_day=factor(Time_day)) %>%
             ggplot() +
             geom_boxplot(aes(x=Time_day, y=CV_adj, col=Media_source)) +
             scale_color_manual(values=c("red","blue","darkgrey")))
, dat %>%
  mutate(Time_day = factor(Time_day)) %>%
  ggplot() +
  geom_boxplot(aes(x=Time_day, y=log_CFUwater, col=Media_source), show.legend=FALSE)+
  scale_color_manual(values=c("red","blue","darkgrey"))+
  ylab("CFU (Water)") + xlab("Day")
, nrow=2)
dev.off()

pdf("richness.pdf", height=4, width=4)
dat %>%
  filter(Richness!=0) %>%
  mutate(Time_day=factor(Time_day)) %>%
  ggplot() +
  geom_boxplot(aes(x=Time_day, y=Richness, col=Media_source)) +
  scale_color_manual(values=c("red","blue"))
dev.off()

### Stats
anova(lm(CV_adj ~ Media_source*Time_day, data=dat))

dat.temp <- dat %>%
  filter(!is.na(log_CFUfilter), is.finite(log_CFUfilter))

anova(lm(CFU_estimate_filter ~ Media_source*Time_day, data=dat.temp))

dat.3 <- dat %>%
  filter(Time_day==3)
Anova(lm(Richness ~ Media_source, data=dat.3))

### Sterile nutrient tubes?

pdf("contaminated_tubes.pdf", width=5, height=3)
dat %>%
  unite(Rep, Media_source,col=ReplicateTreatment, remove=FALSE) %>%
  mutate(Contaminated=as.factor(CFU_estimate_tube)
         , ReplicateTreatment=factor(ReplicateTreatment, levels=c("1_CON","1_IN","2_IN","3_IN","1_OUT","2_OUT","3_OUT"))) %>%
  ggplot() +
  geom_point(aes(x=Time_day, y=ReplicateTreatment, col=Contaminated))
dev.off()


