##### plotting cv biofilm assay ######
library(tidyverse)
dat <- read.delim("biofilm_assay_allisolates.txt", stringsAsFactors = FALSE)

dat_filt <- dat %>%
  filter(TubeID!="EDGE")

dat_filt <- dat_filt %>%
  group_by(TubeID) %>%
  summarize(ave_cv = mean(CV)) %>%
  right_join(dat_filt)

# Make genus and species ID unique; get first of every tube
dat_filt <- dat_filt %>%
  select(TubeID, g_sp) %>%
  mutate(TubeID_first = duplicated(TubeID)) %>%
  filter(!TubeID_first) %>%
  mutate(g_sp_unique = make.unique(as.character(g_sp))) %>%
  select(TubeID, g_sp_unique) %>%
  full_join(dat_filt)

# Edges only-- may be some contaminated ones.
dat %>%
  filter(TubeID=="EDGE") %>%
  ggplot() +
  geom_histogram(aes(x=CV))

# Plotting boxplot of ordered
dat_filt %>%
  arrange(ave_cv) %>%
  mutate(TubeID=factor(TubeID, levels=unique(TubeID))) %>%
  ggplot(aes(x=factor(TubeID), y=CV)) +
  geom_boxplot()+
  geom_point(aes(col=factor(PlateID)))

CON_only <- dat_filt %>%
  filter(TubeID=="CON") 

# Plotting histogram

dat_filt %>%
  ggplot() +
  geom_histogram(aes(x=CV, y=..density..), bins=40) +
  geom_histogram(data=CON_only, aes(x=CV, y=..density..), bins=40, col="red")


# Fitting normal distribution to controls
con.mean <-  mean(CON_only$CV)
con.sd <- sd(CON_only$CV)
con.se <- sd(CON_only$CV)/sqrt(length(CON_only$CV))
sig_thresh <- qnorm(c(0.975), mean=con.mean, sd=con.sd)
sig_thresh_se <- qnorm(c(0.975), mean=con.mean, sd=con.se)


dat_biofilm_members <- dat_filt %>%
  group_by(TubeID) %>%
  summarize(mean=mean(CV), sd=sd(CV)) %>%
  mutate(sig_lwr = qnorm(0.025, mean=mean, sd=sd)) %>%
  mutate(sig_lwr_se = qnorm(0.025, mean=mean, sd=sd/sqrt(3))) %>%
  mutate(biofilm_former_sd = ifelse(sig_lwr>sig_thresh, T,F)) %>%
  mutate(biofilm_former_se = ifelse(sig_lwr_se>sig_thresh_se, T,F))

dat_filt_withbiofilm <- dat_filt %>%
  left_join(dat_biofilm_members) %>%
  arrange(ave_cv) %>%
  mutate(g_sp_unique=factor(g_sp_unique, levels=unique(g_sp_unique))) %>%
  mutate(CAUTION=ifelse(cloudy<2, ave_cv, NA))

ggsave("biofilm_former_sd.pdf", width=9, height=5
       , dat_filt_withbiofilm %>%
         ggplot(aes(x=factor(g_sp_unique), y=CV)) +
         geom_boxplot(aes(col=biofilm_former_sd), show.legend = FALSE)+
         # geom_point(aes(y=CAUTION), col="darkred",alpha=0.5) +
         xlab("IsolateID")+
         ylab("CV intensity (OD@590nm)") +
         theme(axis.text.x=element_text(angle=90, size=6, hjust = 1, vjust=0.5))
       )
ggsave("biofilm_former_sd_crosscheck.pdf", width=9, height=5
       , dat_filt_withbiofilm %>%
         ggplot(aes(x=factor(g_sp_unique), y=CV)) +
         geom_boxplot(aes(col=biofilm_former_sd), show.legend = FALSE)+
         geom_point(aes(y=CAUTION), col="darkred",alpha=0.5) +
         xlab("IsolateID")+
         ylab("CV intensity (OD@590nm)") +
         theme(axis.text.x=element_text(angle=90, size=6, hjust = 1, vjust=0.5))
)

ggsave("biofilm_former_se.pdf", width=9, height=5
       , dat_filt_withbiofilm %>%
         ggplot(aes(x=factor(g_sp_unique), y=CV)) +
         geom_boxplot(aes(col=biofilm_former_se), show.legend = FALSE)+
         # geom_point(aes(y=CAUTION), col="darkred",alpha=0.5) +
         xlab("Isolate#")+
         ylab("CV intensity (OD@590nm)") +
         theme(axis.text.x=element_text(angle=90, size=6, hjust = 1, vjust=0.5))
)

ggsave("biofilm_former_se_withlegend.pdf", width=12, height=5
       , dat_filt_withbiofilm %>%
         mutate(Biofilm_former = biofilm_former_se) %>%
         ggplot(aes(x=factor(g_sp_unique), y=CV)) +
         geom_boxplot(aes(col=biofilm_former_se), show.legend = TRUE)+
         # geom_point(aes(y=CAUTION), col="darkred",alpha=0.5) +
         xlab("Isolate#")+
         ylab("CV intensity (OD@590nm)") +
         theme(axis.text.x=element_text(angle=90, size=6, hjust = 1, vjust=0.5))
)

ggsave("biofilm_former_se_crosscheck.pdf", width=9, height=5
       , dat_filt_withbiofilm %>%
         ggplot(aes(x=factor(g_sp_unique), y=CV)) +
         geom_boxplot(aes(col=biofilm_former_se), show.legend = FALSE)+
         geom_point(aes(y=CAUTION), col="darkred",alpha=0.5) +
         xlab("Isolate#")+
         ylab("CV intensity (OD@590nm)") +
         theme(axis.text.x=element_text(angle=90, size=6, hjust = 1, vjust=0.5))
)

ggsave("biofilm_data.pdf", width=9, height=5
       , dat_filt_withbiofilm %>%
         ggplot(aes(x=factor(g_sp_unique), y=CV)) +
         geom_boxplot()+
         # geom_point(col="grey",alpha=0.5) +
         xlab("Isolate#")+
         ylab("CV intensity (OD@590nm)") +
         theme(axis.text.x=element_text(angle=90, size=6, hjust = 1, vjust=0.5))
)

dat_filt_withbiofilm %>%
  filter(biofilm_former_se) %>%
  group_by(TubeID, g_sp_unique, IsolateID) %>%
  summarize(ave_cv = mean(ave_cv)) %>%
  arrange(-ave_cv) %>%
  write.table(file="list_se_sig_biofilmformers.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")


#### Is biofilm forming abilities a function of growth speed?
# dat_filt_withbiofilm %>%
  # gather(growth_midpoint1, growth_midpoint2, growth_midpoint3, key=Rep, value=growth) %>%
  # ggplot() + geom_histogram(aes(x=as.numeric(growth)))
  
ggsave("growthrate_vs_biofilmformer.pdf"
       , dat_filt_withbiofilm %>%
         mutate(growth_midpoint1=as.numeric(growth_midpoint1), growth_midpoint2 = as.numeric(growth_midpoint2), growth_midpoint3= as.numeric(growth_midpoint3)) %>%
         mutate(growth_midpoint1 = ifelse(growth_midpoint1>30 | growth_midpoint1<=0, NA, growth_midpoint1)
                , growth_midpoint2 = ifelse(growth_midpoint2>30| growth_midpoint2<=0, NA, growth_midpoint2)
                , growth_midpoint3 = ifelse(growth_midpoint3>30| growth_midpoint3<=0, NA, growth_midpoint3)) %>%
         group_by(TubeID, g_sp_unique, IsolateID, growth_midpoint1, growth_midpoint2, growth_midpoint3) %>%
         summarize(mean_growth=mean(c(growth_midpoint1,growth_midpoint2,growth_midpoint3), na.rm = TRUE)
                   , sd_growth=sd(c(growth_midpoint1,growth_midpoint2,growth_midpoint3), na.rm = TRUE)
                   , mean_cv = mean(CV), sd_cv = sd(CV)) %>%
         ungroup() %>%
         ggplot() +
         geom_point(aes(x=mean_growth, y=mean_cv)) +
         geom_segment(aes(x = mean_growth-sd_growth, xend=mean_growth+sd_growth, y=mean_cv, yend=mean_cv), col="red", alpha=0.2)+
         geom_segment(aes(x = mean_growth, xend=mean_growth, y=mean_cv-sd_cv, yend=mean_cv+sd_cv), col="red", alpha=0.2)+
         xlab("Growth rate (1/2kmax, hours)")+
         ylab("CV intensity (OD@590)")+
         theme(axis.title = element_text(size=16))
       )
  
