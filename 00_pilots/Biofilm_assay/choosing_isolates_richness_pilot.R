### Choosing random isolate combos for richness pilot ####
library(tidyverse)
dat <- read.delim("list_se_sig_biofilmformers.txt")
set.seed(24893)
dat_edited <- dat %>%
  filter(ave_cv>0.18) %>%
  mutate(randomOrder=sample(1:n())) 
# %>%
#   mutate(r1a=NA, r1b=NA
#          , r2=NA
#          , r5=NA
#          , r10=NA
#          , r15=NA
#          , r20=NA
#          , r30=NA
#          , r50=NA) %>%

# First, choose total 50 to include
dat_edited[,"r50"] <- TRUE


