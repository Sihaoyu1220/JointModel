library(tidyverse)
library(ggplot2)
library(ggpubr)

# Import the raw data
data <- read.table(here::here("Chapter3","Data","full2dat.txt"), header = TRUE)

# Viral load detection limit is 2
# Define the viral rebound as two consecutive viral rise
data <- data %>% 
  mutate(`CD4` = ifelse(cd4 > 200, 1, 0)) %>% 
  mutate(censor = ifelse(censor == -1, 1, 0)) %>% 
  mutate(limit =0) %>% 
  group_by(patid) %>% 
  mutate(prev = lag(rna)) %>% 
  mutate(diff = rna - prev) %>% 
  mutate(nextdiff = lead(diff)) %>% 
  mutate(rebound = ifelse(nextdiff > 0 & diff > 0, 1, 0)) %>% 
  mutate(rebound = ifelse(is.na(rebound), 0, rebound)) %>% 
  mutate(reboundtime = ifelse(rebound == 1, day, 0)) %>%
  mutate(reboundtime = ifelse(sum(reboundtime) != 0, max(reboundtime), max(day)+1)) %>% 
  mutate(reboundnew = ifelse(day >= reboundtime, 1, 0)) %>% 
  dplyr::select(1:7, 13) %>% 
  mutate(decay = ifelse(reboundnew == 1, 0, 1))

# Manually define some decay and rebound
data <- data %>% 
  mutate(decay = ifelse(patid == 610470, 1, decay)) %>% 
  mutate(reboundnew = ifelse(decay == 1, 0, 1))
data2 <- data
data2$reboundnew[c(31, 53, 67, 68, 76, 112, 131, 148, 172, 189, 217, 314, 361)] <- 1
data2 <- data2 %>% 
  mutate(decay = ifelse(reboundnew == 1, 0, 1))

# Pivot long
data.long <- pivot_longer(data2, cols=c(rna, CD4), names_to = "name", values_to = "observation")
data.long <- data.long %>% 
  mutate(observation = ifelse(is.na(observation), ".", observation)) %>% 
  mutate(censor = ifelse(name == "CD4", 0, censor))

# New decay/rebound definition
data3 <- data.long %>% group_by(patid) %>% mutate(reboundnew = ifelse(sum(reboundnew)!=0 & decay==1, lead(reboundnew, 2), reboundnew))
data3 <- data3 %>% mutate(decay = ifelse(reboundnew == 1, 0, 1))

# Define Ti
data3 <- data3 %>% 
  group_by(patid) %>% 
  mutate(Ti = case_when(decay == 1 ~ day)) %>% 
  mutate(Ti = max(Ti, na.rm = TRUE))
write.csv(data3, here::here("Chapter3","Data","Cleaned_data","full2dat_cleaned.csv"), row.names = FALSE)
