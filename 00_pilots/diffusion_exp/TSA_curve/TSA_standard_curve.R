library("tidyverse")
library("ggplot2")

dat <- read.csv("myc_tsa_st_23sept2019.csv")

dat %>%
    # filter(Protein.Conc. < 1) %>%
    ggplot(aes(x=g_per_L, y=Protein.Conc.)) + 
    geom_point() +
    geom_smooth(method="lm")

# Fit linear model
dat.lm <- lm(Protein.Conc. ~ g_per_L, data=dat)
# Get 95% prediction intervals
newdat <- data.frame(g_per_L=exp(seq(log(0.000001),log(60), length.out = 500)))
dat.pred <- predict(dat.lm, newdata = newdat, interval = "predict" )
# Get minimum g/L we could detect
dat.pred.0 <- dat.pred %>%
    as_tibble() %>%
    filter(lwr>0) 
min_gpL <- dat.pred.0$fit[1]
# What is that in dilution?
min_gpL/30
#factors of 2
log(min_gpL/30, base=2)

# What is minimum reliable protein reading?
predict(dat.lm, newdata=data.frame(g_per_L=min_gpL/30))

# 1 mL dissolved in how much liquid?
1/(min_gpL/30)

# jars aremax 250mL, so we should be able to detect protein even if 1Ml is dissolved in ALL liquid


###### Now, create new line, excluding the "zeros" ######
dat.filt <- dat %>%
    filter(g_per_L>min_gpL)

dat.filt  %>%
    # filter(Protein.Conc. < 1) %>%
    ggplot(aes(x=g_per_L, y=Protein.Conc.)) + 
    geom_point() +
    geom_smooth(method="lm")

# Fit linear model
dat.filt.lm <- lm(Protein.Conc. ~ g_per_L, data=dat.filt)
summary(dat.filt.lm)
# Get 95% prediction intervals
newdat <- data.frame(g_per_L=exp(seq(log(0.000001),log(60), length.out = 500)))
dat.pred <- predict(dat.filt.lm, newdata = newdat, interval = "predict" )
# Get minimum g/L we could detect
dat.pred.0 <- dat.pred %>%
    as_tibble() %>%
    filter(lwr>0) 
min_gpL <- dat.pred.0$fit[1]
# What is that in dilution?
min_gpL/30
#factors of 2
log(min_gpL/30, base=2)

# What is minimum reliable protein reading?
predict(dat.filt.lm, newdata=data.frame(g_per_L=min_gpL/30))

# 1 mL dissolved in how much liquid?
1/(min_gpL/30)
