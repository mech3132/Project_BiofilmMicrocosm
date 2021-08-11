library("tidyverse")
library("ggplot2")

diff <- read.delim("diffusion_exp.txt")
diff <- diff %>%
    mutate(adj_conc = ifelse(A280 < 0.02, 0, A280))

diff %>%
    # filter(i.o = I) %>%
    ggplot(aes(x=T_elapsed, y= (adj_conc), col=Rep)) + 
    # geom_point() +
    geom_line() +
    facet_wrap(~i.o, scales = "free_y")

pdf("diffusion_both.pdf", width=5, height=3)
diff %>%
    mutate(i.o=ifelse(i.o=="I","Inside","Outside"), Days=T_elapsed/60/24) %>%
    ggplot(aes(x=Days, y= (adj_conc), col=Rep)) + 
    # geom_point() +
    geom_line() +
    facet_wrap(~i.o, scales = "free_y") +
    ylab("Absorbance (280nm)")
dev.off()

pdf("diffusion_both_log.pdf", width=5, height=3)
diff %>%
    mutate(i.o=ifelse(i.o=="I","Inside","Outside"), Days=T_elapsed/60/24) %>%
    rename(Inside_or_Outside=i.o) %>%
    ggplot(aes(x=Days, y= log(adj_conc), col=Rep, lty=Inside_or_Outside)) + 
    # geom_point() +
    geom_line() +
    ylab("log Absorbance (280nm)")
dev.off()

ave.time <- diff %>%
    group_by(TP) %>%
    summarize(Time  = mean(T_elapsed))
diff_diffC <- diff %>%
    select(TP, adj_conc, Rep, i.o) %>%
    left_join(ave.time) %>%
    spread(key=i.o, value=adj_conc) %>%
    mutate(diff_C = I-O)

pdf("diffusion_difference.pdf", width=4, height=3)
diff_diffC %>%
    mutate(Days=Time/60/24) %>%
    ggplot(aes(x=Days, y=diff_C, col=Rep)) +
    geom_line() +
    ylab("Difference in absorbance (In-Out)")
dev.off() 


decay_func <- function(k) {1.25*exp(1)^(k*seq(0,8000))}

fit <- data.frame(Time=seq(0,8000), diff_C = decay_func(-0.000804719), Rep = as.character("FIT"))

diff_diffC %>%
    ggplot(aes(x=Time, y=diff_C, col=Rep)) +
    geom_line() + 
    geom_line(data=fit, aes(x=Time, y=diff_C), col="black")
