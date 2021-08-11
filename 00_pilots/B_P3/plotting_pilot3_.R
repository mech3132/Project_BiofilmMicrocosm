### Plotting cv for pilot 3 ###
library(tidyverse)
dat <- read.csv("B_P3_planning.csv")

dat <- dat %>%
  mutate(treatment_name = factor(treatment_name, levels=c("control","burk","rhodo","chryseo","pseudo","mixed"))) %>%
  rename(Number_species_on_platecount = growth_plate_count) 
# CV
dat %>%
  ggplot() +
  geom_boxplot(aes(x=treatment_name, y=log(Adj_cv_membrane1_methanol))) +
  geom_point(aes(x=treatment_name, y=log(Adj_cv_membrane1_methanol), col=Number_species_on_platecount), position = position_jitter(width=0.2, height=0)) +
  xlab("Treatment")+
  ylab("log CV (normalized)")
  
ggsave(filename = "cv_by_treatment_raw.pdf", width=7, height=3
       , dat %>%
  ggplot() +
  geom_boxplot(aes(x=treatment_name, y=CV_membrane1_methanol)) +
  geom_point(aes(x=treatment_name, y=CV_membrane1_methanol, col=Number_species_on_platecount), position = position_jitter(width=0.2, height=0))+
  xlab("Treatment")+
  ylab("CV intensity (OD @ 590nm)")+
  ylim(c(0,0.9)) 
)

ggsave(filename = "cv_by_treatment_raw_polished.pdf", width=4, height=2.5
       , dat %>%
         mutate(treatment_name=as.character(treatment_name)) %>%
         mutate(treatment_name = ifelse(treatment_name=="control","No_bact"
                                        , ifelse(treatment_name=="mixed","All_bact"
                                                 , ifelse(treatment_name=="burk","Bact_1"
                                                          , ifelse(treatment_name=="rhodo", "Bact_2"
                                                                   , ifelse(treatment_name=="chryseo","Bact_3"
                                                                            , "Bact_4")))))) %>%
         mutate(treatment_name = factor(treatment_name, levels=c("No_bact","Bact_1","Bact_2","Bact_3","Bact_4","All_bact"))) %>%
         ggplot() +
         geom_boxplot(aes(x=treatment_name, y=CV_membrane1_methanol)) +
         # geom_point(aes(x=treatment_name, y=CV_membrane1_methanol, col=Number_species_on_platecount), position = position_jitter(width=0.2, height=0))+
         xlab("Bacteria in biofilm")+
         ylab("Biofilm thickness (CV intensity)")+
         ylim(c(0,0.9)) 
)
# dat %>%
#   ggplot() +
#   geom_boxplot(aes(x=treatment_name, y=log(CV_membrane1_methanol))) +
#   geom_point(aes(x=treatment_name, y=log(CV_membrane1_methanol), col=Number_species_on_platecount), position = position_jitter(width=0.2, height=0))+
#   xlab("Treatment")+
#   ylab("log CV intensity (OD @ 590nm)")

# two sides of CV
ggsave("cv_membrane_comparison.pdf", width=7, height=4
       , dat %>%
         # select(JarID, treatment_type, treatment_name, feeding, Number_species_on_platecount, CV_membrane1_methanol, CV_membrane2) %>%
         # gather(CV_membrane1_methanol, CV_membrane2, key=membrane, value=CV) %>%
         ggplot() +
         geom_point(aes(x=log(CV_membrane1_methanol), y=log(CV_membrane2), col=Number_species_on_platecount, pch=treatment_name)) +
         geom_abline(aes(intercept=0, slope=1)) +
         xlab("Membrane 1 (fixed with methanol)") +ylab("Membrane 2 (not fixed)")
       
       )

    

