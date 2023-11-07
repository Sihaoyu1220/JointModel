library(tidyverse)
library(ggplot2)
library(gridExtra)

raw_data<-read.csv("/Users/gsiii/Library/CloudStorage/OneDrive-UBC/PhD/Unified/Code/raw_data.csv",header = T)
# If RNA_V = -1, indicate the observation as censored and replace RNA_V with RNA_L.
# Take log10 to RNA_V. Convert time from days to months. 
# Set lower boundary (limit) for the censoring interval to 0.
# Remove rows with unobserved RNA_V and time
raw_data <- raw_data %>% 
  mutate(threshold = 40) %>% 
  mutate(censor = ifelse(RNA_V < threshold, 1, 0)) %>% 
  mutate(RNA_V = ifelse(RNA_V < threshold, threshold, RNA_V)) %>% 
  mutate(log10 = log10(RNA_V)) %>% 
  mutate(time = days_from_seroco/30) %>% 
  mutate(limit = 0) %>% 
  filter(!is.na(RNA_V) & !is.na(time))
# Remove outliers if there is any
# if(!is.null(outliers)){
#   raw_data <- raw_data %>% 
#     filter(!(PATIENT %in% outliers))
# }
raw_data <- raw_data %>% 
  group_by(PATIENT) %>% 
  mutate(ind = case_when(treatment == 1 ~ time)) %>%  
  mutate(ind = max(ind, na.rm=TRUE)) 
dup <- raw_data %>% 
  filter(time == ind) %>% 
  mutate(treatment == 0)
dataplot <- rbind(raw_data, dup)

# Split the dataset into before_ART (treatment=1) and after_ART(treatment=0 after treatment=1). 
# Note: for after_ART, ignore all "treatment=0" before "treatment=1".
# Change time after ART to the time since ART interruption 
# For example: if before_ART time is c(1,2,3), after_ART time is c(3.5,4,5), then 
#             the time since ART interruption for after_ART is c(0.5,1,2).
before_ART <- raw_data %>% 
  filter(treatment==1) 
after_ART <- raw_data %>% 
  group_by(PATIENT) %>% 
  mutate(ind = case_when(treatment == 1 ~ time)) %>%  
  mutate(ind = max(ind, na.rm=TRUE)) %>% 
  filter(time > ind) %>% 
  mutate(time = time-ind) %>% 
  dplyr::select(-ind)
# plot raw_data (do not plot any "treatment=0" before "treatment=1")
combined <- rbind(before_ART,after_ART)
combined<-combined%>% 
  arrange(PATIENT,days_from_seroco)
combined$PATIENT<-as.factor(combined$PATIENT)
combined$censor<-as.factor(combined$censor)
combined$treatment<-as.factor(combined$treatment)
ggplot(combined, aes(x=days_from_seroco, y=log10)) + 
  geom_point(aes(shape=factor(censor)),size=2, stroke=0) +
  geom_line(aes(group=PATIENT,color=treatment)) +
  scale_x_continuous("Day") + 
  scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)"))+
  scale_shape_manual(values=c(16, 17),labels = c("No", "Yes"))+
  scale_color_manual(values=c("steelblue","black"),labels = c("Following ART interruption", "During ART"))+
  labs(color = "Phase", shape = "Censored")+ggtitle("Plot for all observations")+
  theme_classic()+theme(text=element_text(size=14),panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"))

id <- c(1,2,6)
subpop <- combined %>% filter(PATIENT %in% id)
ggplot(subpop, aes(x=days_from_seroco, y=log10)) + 
  geom_point(aes(shape=factor(censor)),size=2, stroke=0) +
  geom_line(aes(group=PATIENT,color=treatment)) +
  scale_x_continuous("Day") + 
  scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)"))+
  scale_shape_manual(values=c(16, 17),labels = c("No", "Yes"))+
  scale_color_manual(values=c("steelblue","black"),labels = c("Following ART interruption", "During ART"))+
  labs(color = "Phase", shape = "Censored")+ggtitle("Plot for all observations")+
  theme_classic()+theme(text=element_text(size=14),panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"))
