library(tidyverse)
library(ggplot2)
library(ggpubr)
setwd(here::here())

# Import the cleaned data
data <- read.csv("Chapter3/Data/Cleaned_data/full2dat_cleaned.csv")

# Update observation to the real CD4 value for CD4
# (not the binary value)
data_plot <- data %>% 
  mutate(observation = ifelse(name == "CD4", cd4, observation))

# Random select 4 patients
ind <- c(13,3,7, 37)
pat <- unique(data_plot$patid)[ind]

# Plot for patient 1
d1 <- data_plot %>% filter(patid %in% pat[1]) %>% 
  mutate(observation = ifelse(name == "CD4", observation/50, observation))
p1 <- ggplot(d1, aes(x=day, y = observation)) +
  geom_point(aes(shape=factor(censor)), size=3)+
  geom_line(aes(linetype = name)) +
  scale_shape_manual(name = "Censored", values = c(20, 17),labels=c('No', 'Yes'))+
  scale_linetype_manual(name = "Observation", values = c("rna" = "solid", "CD4" = "dashed"),labels=c('CD4', 'Viral load'))+
  scale_x_continuous("Day") +
  scale_y_continuous(name = bquote("Viral load (in" ~ log[10]~"-scale)"),
                     sec.axis = sec_axis(~.*50, name=bquote("CD4 cell count (in cells /" ~ mm^3 ~ ")") ))+
  theme_classic()+theme(text=element_text(size=14),panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"))

# Plot for patient 2
d2 <- data_plot %>% filter(patid %in% pat[2]) %>% 
  mutate(observation = ifelse(name == "CD4", observation/20, observation))
p2 <- ggplot(d2, aes(x=day, y = observation)) +
  geom_point(aes(shape=factor(censor)), size=3)+
  geom_line(aes(linetype = name)) +
  scale_shape_manual(name = "Censored", values = c(20, 17),labels=c('No', 'Yes'))+
  scale_linetype_manual(name = "Observation", values = c("rna" = "solid", "CD4" = "dashed"),labels=c('CD4', 'Viral'))+
  scale_x_continuous("Day") +
  scale_y_continuous(name = bquote("Viral load (in" ~ log[10]~"-scale)"),
                     sec.axis = sec_axis(~.*20, name=bquote("CD4 cell count (in cells /" ~ mm^3 ~ ")") ))+
  theme_classic()+theme(text=element_text(size=14),panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"),legend.position = "none")

# Plot for patient 3
d3 <- data_plot %>% filter(patid %in% pat[3]) %>% 
  mutate(observation = ifelse(name == "CD4", observation/50, observation))
p3 <- ggplot(d3, aes(x=day, y = observation)) +
  geom_point(aes(shape=factor(censor)), size=3)+
  geom_line(aes(linetype = name)) +
  scale_shape_manual(name = "Censored", values = c(20, 17),labels=c('No', 'Yes'))+
  scale_linetype_manual(name = "Observation", values = c("rna" = "solid", "CD4" = "dashed"),labels=c('CD4', 'Viral'))+
  scale_x_continuous("Day") +
  scale_y_continuous(name = bquote("Viral load (in" ~ log[10]~"-scale)"),
                     sec.axis = sec_axis(~.*50, name=bquote("CD4 cell count (in cells /" ~ mm^3 ~ ")") ))+
  theme_classic()+theme(text=element_text(size=14),panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"),legend.position = "none")

# Plot for patient 4
d4 <- data_plot %>% filter(patid %in% pat[4]) %>% 
  mutate(observation = ifelse(name == "CD4", observation/80, observation))
p4 <- ggplot(d4, aes(x=day, y = observation)) +
  geom_point(aes(shape=factor(censor)), size=3)+
  geom_line(aes(linetype = name)) +
  scale_shape_manual(name = "Censored", values = c(20, 17),labels=c('No', 'Yes'))+
  scale_linetype_manual(name = "Observation", values = c("rna" = "solid", "CD4" = "dashed"),labels=c('CD4', 'Viral'))+
  scale_x_continuous("Day") +
  scale_y_continuous(name = bquote("Viral load (in" ~ log[10]~"-scale)"),
                     sec.axis = sec_axis(~.*80, name=bquote("CD4 cell count (in cells /" ~ mm^3 ~ ")") ))+
  theme_classic()+theme(text=element_text(size=14),panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"),legend.position = "none")

# Save Figure 1
png(file="Figure/Results/Figure 1.png", width=1000, height=500)
ggarrange(p1, p2, p3, p4, nrow = 2, ncol=2, common.legend = TRUE, legend = "right")
dev.off()
