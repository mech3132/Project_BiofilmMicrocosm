#!bin/bash

### Processing Bd microscopy counts and measurements ####
library(tidyverse)
setwd("02_externally_generated_data/Bd_counting/")
dat_count <- read.delim("raw_data/Bd_count_data.txt")
dat_meas <- read.delim("raw_data/Bd_measurement_data.txt")
dat_zeros <- read.delim("raw_data/Bd_count_true_zeros.txt")

## Count data processing ##
colnames(dat_count)
# Look at problematic ones first
dat_count[which(dat_count$too_blurry_to_count==TRUE),c("Rep","Jar","Picture","Count","Notes")] # %>% View() # View makes it pop up in a nice window
# After reviewing all of them, I think we should just omit these counts. They are too difficult to count.
# I filter these images out below.


# Look at true zeros-- basically, went through and manually indicated whether I thought the Jar/Rep and a TRUE zero, or if it was a "I can't tell if pinpricks are sporangia" zero
dat_zeros_clean <- dat_zeros %>%
  mutate(sampleID = paste0("R",Rep, "_J", Jar,"_","Bd")) %>%
  mutate(Rep = paste0("Rep",Rep)) %>%
  rename(JarID = Jar) 

### QUADRANTS note
# 1=TopLeft, 2=TopRight, 3=BottomLeft, 4=BottomRIght, 5=Top, 6=Bottom 7=Left, 8=Right, 9=Diagonla, 10=FULL
dat_count_cleaned <- dat_count %>% 
  mutate(too_blurry_to_count=ifelse(is.na(too_blurry_to_count), FALSE, too_blurry_to_count)) %>% #### Blurry counts ####
  filter(!too_blurry_to_count) %>% # Filter out all "NOT" "too_blurry_to_count"
  mutate(pinpricks_present = ifelse(pinpricks_present%in% c("TRUE","TRUE?") , TRUE, FALSE)) %>%  #### Pinpricks ####
  mutate(clouds_present = ifelse(is.na(clouds_present), FALSE, clouds_present))%>%   #### Clouds ####
  mutate(very_small = ifelse(very_small %in% c("TRUE","Possibly","TRUE?"), TRUE, FALSE)) %>%   #### very small ####
  mutate(Count = ifelse(Count=="??", NA, as.numeric(Count))) %>%   #### very small ####
  rowwise() %>%   #### Fixing quadrants ####
  mutate(Full = !any(!is.na(c(TopLeft, TopRight, BottomLeft, BottomRight)))
         , Top = sum(c(TopLeft, TopRight, is.na(c(BottomLeft, BottomRight))))==4
         , Bottom = sum(c(BottomLeft, BottomRight, is.na(c(TopLeft, TopRight ))))==4
         , Left = sum(c(BottomLeft, TopLeft, is.na(c(BottomRight, TopRight ))))==4
         , Right = sum(c(BottomRight, TopRight, is.na(c(BottomLeft, TopLeft ))))==4
         , OnlyTopLeft = sum(c(TopLeft, is.na(c(BottomLeft, TopRight, BottomRight ))))==4
         , OnlyTopRight = sum(c(TopRight, is.na(c(BottomLeft, TopLeft, BottomRight ))))==4
         , OnlyBottomLeft = sum(c(BottomLeft, is.na(c(TopLeft, TopRight, BottomRight ))))==4
         , OnlyBottomRight = sum(c(BottomRight, is.na(c(BottomLeft, TopRight, TopLeft ))))==4
         , Diagonal = (sum(c(BottomRight, TopLeft, is.na(c(TopRight, BottomLeft))))==4 |sum(c(TopRight, BottomLeft, is.na(c(BottomRight, TopLeft))))==4)
         )  %>%
  mutate(Top = ifelse(is.na(Top), FALSE, Top)
         , Bottom = ifelse(is.na(Bottom), FALSE, Bottom)
         , Left = ifelse(is.na(Left), FALSE, Left)
         , Right = ifelse(is.na(Right), FALSE, Right)
         , OnlyTopLeft = ifelse(is.na(OnlyTopLeft), FALSE, OnlyTopLeft)
         , OnlyTopRight = ifelse(is.na(OnlyTopRight), FALSE, OnlyTopRight)
         , OnlyBottomLeft = ifelse(is.na(OnlyBottomLeft), FALSE, OnlyBottomLeft)
         , OnlyBottomRight = ifelse(is.na(OnlyBottomRight), FALSE, OnlyBottomRight)
         , Diagonal = ifelse(is.na(Diagonal), FALSE, Diagonal)
  ) %>%
  mutate(image_quadrant = ifelse(Full, 10, ifelse(OnlyTopLeft,1,ifelse(OnlyTopRight,2,ifelse(OnlyBottomLeft,3,ifelse(OnlyBottomRight,4,ifelse(Top,5,ifelse(Bottom,6,ifelse(Left,7,ifelse(Right,8,ifelse(Diagonal,9,NA))))))))))) %>%   #### Manual area in put ####
  mutate(Area_whole = 0.16651143, Area_unit = "mm2"
         , Area_counted = ifelse(image_quadrant %in% c(1,2,3,4), Area_whole/4, ifelse(image_quadrant %in% c(5,6,7,8,9), Area_whole/2, Area_whole))
         ,sampleID = paste0("R",Rep, "_J", Jar,"_","Bd")) %>%
  select(sampleID, Rep, Jar, Picture, Count, very_small, pinpricks_present, clouds_present, image_quadrant, Area_counted, Area_unit, Notes) %>% #### remove unwanted cols####
  filter(!is.na(Rep)) %>%
  mutate(zoosp_dens_perimage = Count/Area_counted
    , Rep = paste0("Rep",Rep))%>%
  rename(JarID = Jar)

# Aggregate counts by Sample
dat_count_agg <- dat_count_cleaned %>%
  group_by(sampleID, Rep, JarID) %>%
  summarize(zoosp_dens_micro = sum(Count, na.rm=TRUE)/sum(Area_counted)
            , zoosp_count_sd_micro = sd(zoosp_dens_perimage, na.rm=TRUE)
            , max_zoosp_count_micro = max(Count, na.rm=TRUE)
            , min_zoosp_count_micro = min(Count, na.rm=TRUE)
            , pinpricks_present = any(pinpricks_present)
            , clouds_present = any(clouds_present)) %>% left_join(dat_zeros_clean)
  
## Plotting variability of counts
dat_count_cleaned %>%
  ggplot()+geom_point(aes(x=as.factor(JarID), y=zoosp_dens_perimage, col=pinpricks_present)) +
  facet_wrap(.~Rep, nrow=3)

##### Now to measurement data #####
dat_meas_cleaned <- dat_meas %>% 
  filter(!is.na(Diameter)) %>% # remove unsures
  mutate(conv_pix_in_0.1mm = 1031.2316
         , diameter_mm = Diameter/1031.2316*0.1
         , diameter_um = diameter_mm*1000) %>%
  mutate(sampleID = paste0("R",Rep, "_J", Jar,"_","Bd")
         ,Rep = paste0("Rep",Rep)) %>%
  rename(JarID = Jar)

# Aggergate data
dat_meas_agg <- dat_meas_cleaned %>%
  group_by(sampleID, Rep, JarID) %>%
  summarize(ave_zoosp_size_micro = mean(diameter_um, na.rm=TRUE)
            , med_zoosp_size_micro = median(diameter_um, na.rm=TRUE)
            , sd_zoosp_size_micro = sd(diameter_um, na.rm=TRUE)
            , min_zoosp_size_micro = min(diameter_um, na.rm=TRUE)
            , max_zoosp_size_micro = max(diameter_um, na.rm=TRUE)
            , n_zoosp_size_micro = length(.)) 

## Plotting variability of size
dat_meas_cleaned %>%
  ggplot() + geom_violin(aes(x=as.factor(JarID), y=diameter_um)) +
  facet_wrap(.~Rep)

### Write out for downstream

write.table(dat_count_cleaned, "downstream/dat_count_clean.txt", quote=FALSE, row.names = FALSE,sep="\t")
write.table(dat_count_agg, "downstream/dat_count_agg.txt", quote=FALSE, row.names = FALSE,sep="\t")

write.table(dat_meas_cleaned, "downstream/dat_meas_clean.txt", quote=FALSE, row.names = FALSE,sep="\t")
write.table(dat_meas_agg, "downstream/dat_meas_agg.txt", quote=FALSE, row.names = FALSE,sep="\t")

setwd("../..")

