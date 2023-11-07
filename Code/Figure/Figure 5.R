library(ggpubr)
library(tidyverse)
library(ggplot2)
library(splines)

indparest <- read.csv(here::here("Code/Chapter3/indest.csv"), header = TRUE)
data <- read.csv(here::here("Data/Cleaned_data/full2dat_cleaned.csv"), header = TRUE)
# patid 610282 
pat1 <- data %>% filter(patid == 610282 & name == "rna") 
indparest1 <- indparest %>% filter(id == 610282)
P1 <- indparest1$p1
lambda1 <- indparest1$b1
P2 <- indparest1$p2
beta1 <- indparest1$beta1
beta2 <- indparest1$beta2
beta3 <- indparest1$beta3
beta4 <- indparest1$beta4
decaylast <- unique(pat1$Ti)
rebound1st <- min(pat1[which(pat1$decay == 0), ]$day)
datamid <- data.frame(day = seq(decaylast, rebound1st, 0.1))
datamid <- datamid %>% 
  mutate(ydecay = log10( exp(P1-lambda1*day) + exp(P2))) %>% 
  mutate(yrebound = beta1 * (day-decaylast) / ((day-decaylast)  + exp(beta2-beta3*(day-decaylast))) + beta4) %>% 
  group_by(day) %>% 
  mutate(observation = max(ydecay, yrebound))

fun1 <- function(x){
  ifelse(x<=decaylast, log10( exp(P1-lambda1*x) + exp(P2)), NA)
}
fun2 <- function(x){
  ifelse(x>=rebound1st, beta1 * (x-decaylast) / ((x-decaylast)  + exp(beta2-beta3*(x-decaylast))) + beta4, NA)
}

(p1 <- ggplot(pat1, aes(x = day, y = observation)) + geom_point(size=3) + 
    stat_function(fun = fun1) + stat_function(fun = fun2) + geom_smooth(data=datamid, se=FALSE, color = "black",size=0.5) +
    xlab("Time") + scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)")) + theme_classic() + theme(text = element_text(size = 20)))


# patid 610434
pat2 <- data %>% filter(patid == 610434 & name == "rna") 
indparest2 <- indparest %>% filter(id == 610434)
P1 <- indparest2$p1
lambda1 <- indparest2$b1
P2 <- indparest2$p2
beta1 <- indparest2$beta1
beta2 <- indparest2$beta2
beta3 <- indparest2$beta3
beta4 <- indparest2$beta4

decaylast <- unique(pat2$Ti)
rebound1st <- min(pat2[which(pat2$decay == 0), ]$day)
datamid2 <- data.frame(day = seq(decaylast-0.3, rebound1st, 0.1))
datamid2 <- datamid2 %>% 
  mutate(ydecay = log10( exp(P1-lambda1*day) + exp(P2))) %>% 
  mutate(yrebound = beta1 * (day-decaylast) / ((day-decaylast)  + exp(beta2-beta3*(day-decaylast))) + beta4) %>% 
  group_by(day) %>% 
  mutate(observation = max(ydecay, yrebound))

fun3 <- function(x){
  ifelse(x<=decaylast, log10( exp(P1-lambda1*x) + exp(P2)), NA)
}
fun4 <- function(x){
  ifelse(x>=rebound1st, beta1 * (x-decaylast) / ((x-decaylast)  + exp(beta2-beta3*(x-decaylast))) + beta4, NA)
}

(p2 <- ggplot(pat2, aes(x = day, y = observation)) + geom_point(size=3) + 
    stat_function(fun = fun3) + stat_function(fun = fun4) + geom_smooth(data=datamid2, se=FALSE, color = "black",size=0.5, span=0.3) +
    xlab("Time") + scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)")) + theme_classic() + theme(text = element_text(size = 20)))

# patid 271775
pat3 <- data %>% filter(patid == 271775 & name == "rna") 
indparest3 <- indparest %>% filter(id == 271775)
P1 <- indparest3$p1
lambda1 <- indparest3$b1
P2 <- indparest3$p2
beta1 <- indparest3$beta1
beta2 <- indparest3$beta2
beta3 <- indparest3$beta3
beta4 <- indparest3$beta4
decaylast <- unique(pat3$Ti)
rebound1st <- min(pat3[which(pat3$decay == 0), ]$day)
datamid <- data.frame(day = seq(decaylast-3, rebound1st+3, 0.1))
datamid <- datamid %>% 
  mutate(ydecay = log10( exp(P1-lambda1*day) + exp(P2))) %>% 
  mutate(yrebound = beta1 * (day-decaylast) / ((day-decaylast)  + exp(beta2-beta3*(day-decaylast))) + beta4) %>% 
  mutate(observation = ifelse(day <= 0.5*(decaylast+rebound1st), ydecay, yrebound))

fun1 <- function(x){
  ifelse(x<=decaylast, log10( exp(P1-lambda1*x) + exp(P2)), NA)
}
fun2 <- function(x){
  ifelse(x>=rebound1st, beta1 * (x-decaylast) / ((x-decaylast)  + exp(beta2-beta3*(x-decaylast))) + beta4, NA)
}

(p3 <- ggplot(pat3, aes(x = day, y = observation)) + 
    geom_point(aes(shape=factor(censor)),size=3,  stroke=0) + 
    scale_shape_discrete(labels=c("Observed data", "Censored data"),name="Data")+
    stat_function(fun = fun1) + stat_function(fun = fun2) + geom_smooth(data=datamid, se=FALSE, color = "black",size=0.5, span = 1) +
    xlab("Time") + scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)")) + theme_classic()  + theme(text = element_text(size = 20)))

# patid 610429
pat4 <- data %>% filter(patid == 610429 & name == "rna") 
indparest4 <- indparest %>% filter(id == 610429)
P1 <- indparest4$p1
lambda1 <- indparest4$b1
P2 <- indparest4$p2
beta1 <- indparest4$beta1
beta2 <- indparest4$beta2
beta3 <- indparest4$beta3
beta4 <- indparest4$beta4


decaylast <- unique(pat4$Ti)
rebound1st <- min(pat4[which(pat4$decay == 0), ]$day)
datamid <- data.frame(day = seq(decaylast-1, rebound1st+1, 0.1))
datamid <- datamid %>% 
  mutate(ydecay = log10( exp(P1-lambda1*day) + exp(P2))) %>% 
  mutate(yrebound = beta1 * (day-decaylast) / ((day-decaylast)  + exp(beta2-beta3*(day-decaylast))) + beta4) %>% 
  mutate(observation = ifelse(day <= 0.5*(decaylast+rebound1st), ydecay, yrebound))

fun1 <- function(x){
  ifelse(x<=decaylast, log10( exp(P1-lambda1*x) + exp(P2)), NA)
}
fun2 <- function(x){
  ifelse(x>=rebound1st, beta1 * (x-decaylast) / ((x-decaylast)  + exp(beta2-beta3*(x-decaylast))) + beta4, NA)
}

(p4 <- ggplot(pat4, aes(x = day, y = observation)) + 
    geom_point(aes(shape=factor(censor)),size=3,  stroke=0) + 
    scale_shape_discrete(labels=c("Observed data", "Censored data"),name="Data")+
    stat_function(fun = fun1) + stat_function(fun = fun2) + geom_smooth(data=datamid, se=FALSE, color = "black",size=0.5, span = 1) +
    xlab("Time") + scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)")) + theme_classic()  + theme(text = element_text(size = 20)))

