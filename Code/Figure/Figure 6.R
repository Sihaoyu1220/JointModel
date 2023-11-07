library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ggpubr)
setwd("/Users/gsiii/Library/CloudStorage/OneDrive-UBC/PhD/Unified/Code")

######################## main function: covariates on beta1 in rebound model =======
# filename is a string to indicate the raw data file name. 
# Raw data includes "PATIENT","CD4_V","RNA_V","RNA_L","treatment", and "days_from_seroco"
# num_bootstrap is the number of bootstrap want to repeat
# outliers is a vector of PATIENT ID which are outliers, e.g., outliers = c(1,5,10)
# num_months is a number to indicate the time interval for viral rebound analysis. 
filename = "raw_data.csv"
# outliers=c(4,6,12,13,14,15,27,37,60,61)
outliers = c(4,20,27)
# num_months=9

######################### Clean raw data ================
# Read raw data
raw_data<-read.csv(filename,header = T)
# If RNA_V = -1, indicate the observation as censored and replace RNA_V with RNA_L.
# Take log10 to RNA_V. Convert time from days to months. 
# Set lower boundary (limit) for the censoring interval to 0.
# Remove rows with unobserved RNA_V and time
raw_data <- raw_data %>% 
  mutate(censor = ifelse(RNA_V==-1, 1, 0)) %>% 
  mutate(DL = log10(40)) %>% 
  mutate(RNA_V = ifelse(RNA_V==-1, DL, RNA_V)) %>% 
  mutate(log10 = log10(RNA_V)) %>% 
  mutate(time = days_from_seroco/30) %>% 
  mutate(limit = 0) %>% 
  filter(!is.na(RNA_V) & !is.na(time))
# Remove outliers if there is any
if(!is.null(outliers)){
  raw_data <- raw_data %>% 
    filter(!(PATIENT %in% outliers))
}
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

# lead the treatment time for a better plot
data_replicate <- combined %>% group_by(PATIENT, treatment) %>% slice(c(n())) %>% filter(treatment == 1) %>% mutate(treatment = 0)
data_plot <- rbind(combined, data_replicate) %>% arrange(PATIENT, days_from_seroco)
plot_all <- ggplot(data_plot, aes(x=days_from_seroco, y=log10, group=interaction(PATIENT, treatment))) + 
  geom_point(aes(shape=factor(censor)),size=3) +
  geom_line(aes(linetype=treatment)) +
  scale_x_continuous("Day") + 
  scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)"))+
  scale_linetype_manual(name = "Phase", values = c("1" = "solid", "0" = "dashed"), labels=c('rebound','decay'))+
  scale_shape_manual(name = "Censored", values = c(20, 17),labels=c('No', 'Yes'))+
  ggtitle("(A) Viral load trajectories for all patients")+
  theme_classic()+theme(text=element_text(size=14),panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"),legend.position = "none")

# Plot viral load trajectories for randomly selected 4 patients
ind5 <- c(7, 12, 38, 46, 50)
pat5 <- unique(data_plot$PATIENT)[ind5]
pat <- data_plot %>% 
  filter(PATIENT %in% pat5)
plot_random5 <- ggplot(pat, aes(x=days_from_seroco, y=log10, group=interaction(PATIENT, treatment))) + 
  geom_point(aes(shape=factor(censor)),size=3) +
  geom_line(aes(linetype=treatment)) +
  scale_x_continuous("Day") + 
  scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)"))+
  scale_linetype_manual(name = "Phase", values = c("1" = "solid", "0" = "dashed"),labels=c('rebound','decay'))+
  scale_shape_manual(name = "Censored", values = c(20, 17),labels=c('No', 'Yes'))+
  ggtitle("(B) Viral load trajectories for five randomly selected patients")+
  theme_classic()+theme(text=element_text(size=14),panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"))

ggarrange(plot_all, plot_random5, widths = c(0.7, 1))

# Plot decay and rebound separately
ind3 <- c(2,25,70)
pat <- data_plot %>% 
  filter(PATIENT %in% ind3)
pat <- pat %>% 
  group_by(PATIENT) %>% 
  mutate(Ti = case_when(treatment == 1 ~ days_from_seroco)) %>%  
  mutate(Ti = max(Ti, na.rm=TRUE)) %>% 
  mutate(trebound = days_from_seroco-Ti)
before_ART <- pat %>% 
  filter(treatment == 1)
after_ART <- pat %>% 
  filter(treatment == 0)
plot_before <- ggplot(before_ART, aes(x=days_from_seroco, y=log10)) + 
  geom_point(aes(shape=factor(censor)),size=3) +
  geom_line(aes(linetype=PATIENT)) +
  scale_x_continuous("Day") + 
  scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)"))+
  scale_shape_manual(name = "Censored", values = c(20, 17),labels=c('No', 'Yes'))+
  ggtitle("(A) Viral decay trajectories")+
  guides(linetype = FALSE) +
  theme_classic()+theme(text=element_text(size=13),panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"))
plot_after <- ggplot(after_ART, aes(x=trebound, y=log10)) + 
  geom_point(aes(shape=factor(censor)),size=3) +
  geom_line(aes(linetype=PATIENT)) +
  scale_x_continuous("Day") + 
  scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)"))+
  scale_shape_manual(name = "Censored", values = c(20, 17),labels=c('No', 'Yes'))+
  ggtitle("(B) Viral rebound trajectories")+
  guides(linetype = FALSE) +
  theme_classic()+theme(text=element_text(size=13),panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"))
ggarrange(plot_before, plot_after, common.legend = TRUE, widths = c(0.45,0.45), legend = "right")

## summary statistics
raw_data<-read.csv(filename,header = T)
length(unique(raw_data$PATIENT))
t <- raw_data %>% group_by(PATIENT) %>% summarise(max(days_from_seroco))
summary(t)
ni <- raw_data %>% group_by(PATIENT) %>% summarise(n = n())
summary(ni)
summary(raw_data$censor)
duringART <- raw_data %>% 
  filter(treatment == 1)
t <- duringART %>% group_by(PATIENT) %>% summarise(max(days_from_seroco))
summary(t)
ni <- duringART %>% group_by(PATIENT) %>% summarise(n = n())
summary(ni)
summary(as.factor(duringART$RNA_V))
414/695
afterART <- raw_data %>% 
  filter(treatment == 0)
t <- afterART %>% group_by(PATIENT) %>% summarise(max(days_from_seroco-ati_start))
summary(t)
ni <- afterART %>% group_by(PATIENT) %>% summarise(n = n())
summary(ni)
summary(as.factor(afterART$RNA_V))
122/1199

