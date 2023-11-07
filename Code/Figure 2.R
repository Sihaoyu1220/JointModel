library(tidyverse)
library(ggplot2)
library(ggpubr)

# plot
data <- read.csv("Data/Cleaned_data/full2dat_cleaned.csv")
data2 <- data %>% filter(name == "rna")
data2$reboundnew <- as.factor(data2$reboundnew)
data2$censor <- as.factor(data2$censor)
p1 <- ggplot(data2, aes(x=day, y=observation, group=interaction(patid, reboundnew))) +
  geom_point(aes(shape=censor),size=3) +
  geom_line(aes(linetype = reboundnew)) +
  scale_x_continuous("Day") +
  scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)"))+
  scale_linetype_manual(name = "Phase", values = c("0" = "solid", "1" = "dashed"),labels=c('decay', 'rebound'))+
  scale_shape_manual(name = "Censored", values = c(20, 17),labels=c('No', 'Yes'))+
  ggtitle("(B) Viral load trajectories for all patients")+
  theme_classic()+theme(text=element_text(size=14),panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"))

data4 <- data %>% filter(name == "CD4")
original_data <- read.table("Data/full2dat.txt", header = TRUE)
original_data <- original_data %>% select(patid, day, cd4)
data5 <- merge(data4, original_data, by = c("patid", "day"))
data5$decay <- as.factor(data5$decay)
data5 <- data5 %>% arrange(patid, day) %>% 
  group_by(patid) %>% 
  mutate(reboundnext = lead(reboundnew)) %>% 
  mutate(dup = ifelse(reboundnext != reboundnew, 1, 0)) 
dup <- data5 %>% filter(dup==1) %>% mutate(reboundnew = 1)
data6 <- rbind(data5, dup) %>% arrange(patid, day)
data6$reboundnew <- as.factor(data6$reboundnew)
p3 <- ggplot(data6, aes(x=day, y=cd4, group=interaction(patid, reboundnew))) +
  geom_point(size=2) +
  geom_line() + 
  #  geom_line(aes(linetype = reboundnew)) +
  scale_x_continuous("Day") +
  scale_y_continuous(bquote("CD4 cell count (in cells /" ~ mm^3 ~ ")") )+
  ggtitle("(A) CD4 cell count trajectories for all patients")+
  #geom_hline(yintercept=200, linetype="dashed", size=1)+
  theme_classic()+theme(text=element_text(size=14),panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"),legend.position = "none")
# ggarrange(p1,p3, common.legend = TRUE,legend="right")
ggarrange(p3, p1, widths = c(0.85, 1))
