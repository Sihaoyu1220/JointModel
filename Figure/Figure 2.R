library(tidyverse)
library(ggplot2)
library(ggpubr)
setwd(here::here())

# Import the cleaned data
data <- read.csv("Chapter3/Data/Cleaned_data/full2dat_cleaned.csv")

# Select the viral load
viral <- data %>% filter(name == "rna")
viral$reboundnew <- as.factor(viral$reboundnew)
viral$censor <- as.factor(viral$censor)
# Plot for viral load
p1 <- ggplot(viral, aes(x=day, y=observation, group=interaction(patid, reboundnew))) +
  geom_point(aes(shape=censor),size=3) +
  geom_line(aes(linetype = reboundnew)) +
  scale_x_continuous("Day") +
  scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)"))+
  scale_linetype_manual(name = "Phase", values = c("0" = "solid", "1" = "dashed"),labels=c('decay', 'rebound'))+
  scale_shape_manual(name = "Censored", values = c(20, 17),labels=c('No', 'Yes'))+
  ggtitle("(B) Viral load trajectories for all patients")+
  theme_classic()+theme(text=element_text(size=14),panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"))

# Select for CD4 
CD4 <- data %>% filter(name == "CD4")
# Plot for.CD4
p3 <- ggplot(CD4, aes(x=day, y=cd4, group=(patid))) +
  geom_point(size=2) +
  geom_line() + 
  scale_x_continuous("Day") +
  scale_y_continuous(bquote("CD4 cell count (in cells /" ~ mm^3 ~ ")") )+
  ggtitle("(A) CD4 cell count trajectories for all patients")+
  theme_classic()+theme(text=element_text(size=14),panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"),legend.position = "none")

# Save Figure 2
png(file="Figure/Results/Figure 2.png", width=1200, height=600)
ggarrange(p3, p1, widths = c(0.85, 1))
dev.off()