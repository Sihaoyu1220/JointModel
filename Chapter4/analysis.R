library(lixoftConnectors)
library(mvtnorm)
library(nlme)
library(MASS)
library(saemix)
library(tidyverse)
library(testit)
library(here)
library(docopt)
library(ggplot2)
library(lmtest)
library(gridExtra)
library(splines)
library(stats)
library(xtable)
library(GLDEX)
library(ggalt)
initializeLixoftConnectors(software="monolix")
setwd("/Users/gsiii/Library/CloudStorage/OneDrive-UBC/PhD/Unified/Code")
#setwd(here("Desktop","PhD","Unified","Code"))

######################## main function: covariates on beta1 in rebound model =======
# filename is a string to indicate the raw data file name. 
# Raw data includes "PATIENT","CD4_V","RNA_V","RNA_L","treatment", and "days_from_seroco"
# num_bootstrap is the number of bootstrap want to repeat
# outliers is a vector of PATIENT ID which are outliers, e.g., outliers = c(1,5,10)
# num_months is a number to indicate the time interval for viral rebound analysis. 
filename = "raw_data.csv"
# outliers=c(4,6,12,13,14,15,27,37,60,61)
outliers = c(4,20,27)

######################### Clean raw data ================
# Read raw data
raw_data<-read.csv(filename,header = T)
# Set the detection limit to be 40.
raw_data <- raw_data %>%
  mutate(DL = 40) 
# Find lowest detection limit for each patient
# raw_data <- raw_data %>%
#   group_by(PATIENT) %>%
#   mutate(RNA_L = min(RNA_L,na.rm = TRUE)) 
# If RNA_V = -1, indicate the observation as censored and replace RNA_V with RNA_L.
# Take log10 to RNA_V. Convert time from days to months. 
# Set lower boundary (limit) for the censoring interval to 0.
# Remove rows with unobserved RNA_V and time. 
raw_data <- raw_data %>% 
  mutate(censor = ifelse(RNA_V<=40, 1, 0)) %>% 
  mutate(RNA_V = ifelse(censor==1, DL, RNA_V)) %>% 
  mutate(log10 = log10(RNA_V)) %>% 
  mutate(time_raw = days_from_seroco/30) %>% 
  mutate(limit = 0) %>% 
  filter(!is.na(RNA_V) & !is.na(time_raw))
# Remove outliers if there is any
if(!is.null(outliers)){
  raw_data <- raw_data %>% 
    filter(!(PATIENT %in% outliers))
}
# Find interruption time Ti for each patient. 
# Create a new variable "phase" to indicate viral decay and viral rebound.
# Change time for rebound data to the time since ART interruption 
# For example: if during ART time is c(1,2,3), viral rebound time is c(3.5,4,5), then 
#             the time since ART interruption for after_ART is c(0.5,1,2).
# Select useful variables.
data <- raw_data %>% 
  group_by(PATIENT) %>% 
  mutate(Ti = case_when(treatment == 1 ~ time_raw)) %>%  
  mutate(Ti = max(Ti, na.rm=TRUE)) %>% 
  mutate(phase = ifelse(time_raw <= Ti, "decay", "rebound")) %>% 
  mutate(time = ifelse(phase=="decay", time_raw, time_raw-Ti)) %>% 
  mutate(decay = ifelse(phase == "decay", 1, 0)) %>% 
  mutate(rebound = ifelse(phase == "rebound", 1, 0))
# Remove outliers
data <- data %>% 
  filter(!(PATIENT == 4 & phase == "rebound")) %>% 
  filter(!(PATIENT == 30 & phase == "rebound")) %>% 
  filter(!(PATIENT == 5 & days_from_seroco == 1089)) %>% 
  filter(!(PATIENT == 11 & days_from_seroco == 712)) %>% 
  filter(!(PATIENT == 11 & days_from_seroco == 748)) %>% 
  filter(!(PATIENT == 12 & days_from_seroco == 989)) %>% 
  filter(!(PATIENT == 17 & days_from_seroco == 740)) %>% 
  filter(!(PATIENT == 27)) %>% 
  filter(!(PATIENT == 35 & days_from_seroco == 1181)) %>% 
  filter(!(PATIENT == 42 & days_from_seroco == 1190)) %>% 
  filter(!(PATIENT == 32 & days_from_seroco ==  742)) %>% 
  filter(!(PATIENT == 32 & days_from_seroco ==  810))
# Split before_ART and after_ART based on phase.
before_ART <- data %>% filter(phase=="decay") 
after_ART <- data %>% filter(phase=="rebound")
# Save the two datasets locally. This is required when run Monolix.
write.csv(before_ART,"before_ART.csv",row.names = FALSE)
write.csv(after_ART,"after_ART.csv",row.names = FALSE)
write.csv(data,"cleaned_data.csv",row.names = FALSE)
# Read cleaned dataset
before_ART <- read.csv("before_ART.csv")
after_ART <- read.csv("after_ART.csv")
data <- read.csv("cleaned_data.csv")

######################## Nonparametric ========
# AICs and BICs for CD4_V and log_10 of the viral loads.
lme1.1 <- lme(log10~1, data=after_ART, random=~1|PATIENT, method="ML")

basis2 <- ns(after_ART$time, df=1, intercept=F)
after_ART2 <- data.frame(after_ART, basis2)
names(after_ART2)[ncol(after_ART2)] <- "z1"
lme2.1 <- lme(log10~z1, data=after_ART2, random=~1|PATIENT, method="ML")
lme2.2 <- lme(log10~z1, data=after_ART2, random=~z1|PATIENT, method="ML")

basis3 <- ns(after_ART$time, df=2, intercept=F)
after_ART3 <- data.frame(after_ART, basis3)
names(after_ART3)[c(ncol(after_ART3)-1,ncol(after_ART3))] <- c("z1","z2")
lme3.1 <- lme(log10~z1+z2, data=after_ART3, random=~1|PATIENT, method="ML")
lme3.2 <- lme(log10~z1+z2, data=after_ART3, random=~z1|PATIENT, method="ML")
lme3.3 <- lme(log10~z1+z2, data=after_ART3, random=~z1+z2|PATIENT, method="ML",control =list(msMaxIter = 1000, msMaxEval = 1000))

basis4 <- ns(after_ART$time, df=3, intercept=F)
after_ART4 <- data.frame(after_ART, basis4)
names(after_ART4)[c(ncol(after_ART4)-2,ncol(after_ART4)-1,ncol(after_ART4))] <- c("z1","z2","z3")
lme4.1 <- lme(log10~z1+z2+z3, data=after_ART4, random=~1|PATIENT, method="ML")
lme4.2 <- lme(log10~z1+z2+z3, data=after_ART4, random=~z1|PATIENT, method="ML")
lme4.3 <- lme(log10~z1+z2+z3, data=after_ART4, random=~z1+z2|PATIENT, method="ML")
lme4.4 <- lme(log10~z1+z2+z3, data=after_ART4, random=~z1+z2+z3|PATIENT, method="ML",control =list(msMaxIter = 1000, msMaxEval = 1000))

anova(lme1.1, lme2.1,lme2.2,lme3.1,lme3.2,lme3.3,lme4.1,lme4.2,lme4.3,lme4.4)
# Select model lme3.3 for log_10 of the viral loads.
pred <- predict(lme3.3)
# 
# anova(cd4.lme1.1, cd4.lme2.1,cd4.lme2.2,cd4.lme3.1,cd4.lme3.2,cd4.lme3.3,cd4.lme4.1,cd4.lme4.2,
#       cd4.lme4.3,cd4.lme4.4)
# # Select model lme3.3 for CD4_V of the viral loads.
# cd4.pred <- predict(cd4.lme3.3)
# 
# ind <- which.na(after_ART3$CD4_V)
# for(i in 0:(length(ind)-1)) 
#   cd4.pred <- append(cd4.pred, NA, after=(ind[i+1]-1))
# cd4preddata <- cbind(after_ART3, cd4.pred)

# par(mfrow=c(3,3))
# ID <- unique(after_ART3$PATIENT)
# for(i in 1:4){
#   subdata <- cd4preddata[cd4preddata$PATIENT==ID[i],]
#   completedata <- subset(subdata, !is.na(cd4.pred))
#   plot(completedata$time, completedata$cd4.pred, type="l", 
#        ylim=c(0, max(subdata$CD4_V, na.rm = TRUE)), 
#        xlab="Time after treatment (re-scaled)", 
#        ylab="CD4 count")
#   points(subdata$time, subdata$CD4_V, pch=1)
# }

# 
basis3_before <- ns(before_ART$time, df=2, intercept=F)
before_ART3 <- data.frame(before_ART, basis3_before)
names(before_ART3)[c(ncol(before_ART3)-1,ncol(before_ART3))] <- c("z1","z2")

fulldata <- rbind(before_ART3, after_ART3)

fulldata[is.na(fulldata)] <- '.'
fulldata$CD4_V <- as.numeric(fulldata$CD4_V)

# threeobs <- fulldata %>% 
#   pivot_longer(c(CD4_V, log10), values_to = "observation", names_to = "name") %>% 
#   mutate(name = ifelse(name == "log10", phase, name)) %>% 
#   mutate(censor = ifelse(name == "CD4_V", 0, censor)) %>% 
#   mutate(observation = ifelse(name=="CD4_V"&phase=="decay", '.', observation))
# write.csv(threeobs,"threeobs.csv",row.names = FALSE)

# threeobslog <- fulldata %>% 
#   mutate(CD4_V = log10(CD4_V)) %>% 
#   pivot_longer(c(CD4_V, log10), values_to = "observation", names_to = "name") %>% 
#   mutate(censor = ifelse(name == "CD4_V", 0, censor)) %>% 
#   mutate(observation = ifelse(name=="CD4_V"&phase=="decay", '.', observation))
# write.csv(threeobslog,"threeobslog.csv",row.names = FALSE)

threeobssqrt <- fulldata %>% 
  mutate(CD4_V = sqrt(CD4_V)) %>% 
  pivot_longer(c(CD4_V, log10), values_to = "observation", names_to = "name") %>% 
  mutate(name = ifelse(phase=="decay"&name=="log10", "decay", ifelse(phase=="rebound"&name=="log10", "rebound", name))) %>% 
  mutate(censor = ifelse(name == "CD4_V", 0, censor)) %>% 
  filter(!(name=="CD4_V"& phase=="decay")) %>% 
  filter(!(observation == ".")) %>% 
  dplyr::select(PATIENT, censor, time_raw, Ti, limit, z1, z2, name, observation)
write.csv(threeobssqrt,"threeobssqrt.csv",row.names = FALSE)


######################## Find population parameter estimates before ART interruption in Monolix ============
# model1: y_ij = log10(exp(P_1i-lambda_1i*t_ij)+exp(P_2i-lambda_2i*t_ij))+e_ij
# P_1i = P_1 + b_1i, P_2i = P_2 + b_2i, lambda_1i = lambda_1+b_3i, lambda_2i = lambda_2+b_4i
# b_i ~ N(0,B), e_ij ~ N(0, sigma^2)
threeobssqrt_basis3 <- read.csv("/Users/gsiii/Library/CloudStorage/OneDrive-UBC/PhD/Unified/Code/threeobssqrt.csv")
data = list(dataFile = paste0('/Users/gsiii/Library/CloudStorage/OneDrive-UBC/PhD/Unified/Code/threeobssqrt.csv'),
            headerTypes =c("id","cens", "time","regressor","limit","regressor","regressor","obsid","observation"),
            observationTypes = list(CD4_V = "continuous", decay = "continuous", rebound = "continuous"))
#mapping = list("CD4_V"= "yCD4_V", "viral"="yviral"))
modelFile = paste0('/Users/gsiii/Library/CloudStorage/OneDrive-UBC/PhD/Unified/Code/threeobsmodel5_reg.txt')
# create a new project by setting a data set and a structural model
newProject(data = data, modelFile =modelFile)
# getObservationInformation()
setErrorModel(list(yCD4_V="constant",  ydecay = "constant",  yrebound = "constant"))
# covariates on beta50
# getCovariateInformation()
# set tasks in scenario
scenario <- getScenario()
scenario$tasks = c(populationParameterEstimation = T,
                   conditionalModeEstimation = T,
                   conditionalDistributionSampling = T,
                   standardErrorEstimation=T,
                   logLikelihoodEstimation=T)
scenario$linearization = FALSE
# setPopulationParameterEstimationSettings(nbexploratoryiterations=1000, nbsmoothingiterations=500)
# setStandardErrorEstimationSettings(miniterations = 200)
setIndividualParameterVariability(beta50 = FALSE)
setIndividualParameterDistribution(p1="normal",b1="normal",p2="normal",b2="normal", beta1="normal", beta2="normal", beta3="normal",beta4="normal",beta50="normal", beta51="normal",beta52="normal",beta53="normal", alpha0="normal",alpha1="normal")
# getVariabilityLevels()

setCorrelationBlocks(ID = list(c("b1","b2", "beta1", "beta3")))
setPopulationParameterInformation(p1_pop = list(initialValue = 17.8), 
                                  b1_pop = list(initialValue = 4.3),
                                  p2_pop = list(initialValue = 3.2),
                                  b2_pop = list(initialValue = 0.07),
                                  beta1_pop = list(initialValue = 3.5),
                                  beta2_pop = list(initialValue = 8),
                                  beta3_pop = list(initialValue = 3.45),
                                  beta4_pop = list(initialValue = 1.4),
                                  beta50_pop = list(initialValue = -0.14),
                                  beta51_pop = list(initialValue =3.5),
                                  beta52_pop = list(initialValue = -6.4),
                                  beta53_pop = list(initialValue = -6.9),
                                  alpha0_pop = list(initialValue = 26.6),
                                  alpha1_pop = list(initialValue = -0.22))
# run the estimation
setScenario(scenario)
runScenario()
# store the estimates
popestimates <- getEstimatedPopulationParameters()
# store the random effects 
# bi <- as.data.frame(getEstimatedRandomEffects()$conditionalMean)
# store the s.e.
popse <- as.numeric(getEstimatedStandardErrors()[["stochasticApproximation"]][["se"]][1:14])
# # Find G
omega_viral <- popestimates[15:27]
G <- matrix(diag(omega_viral)^2,ncol=13)
G[4,3] <- G[3,4] <- popestimates[28] * popestimates[16]^2 * popestimates[18]^2
G[5,3] <- G[3,5] <- popestimates[29] * popestimates[16]^2 * popestimates[19]^2
G[7,3] <- G[3,7] <- popestimates[30] * popestimates[16]^2 * popestimates[21]^2
G[4,5] <- G[5,4] <- popestimates[31] * popestimates[18]^2 * popestimates[19]^2
G[4,7] <- G[7,4] <- popestimates[32] * popestimates[18]^2 * popestimates[21]^2
G[5,7] <- G[7,5] <- popestimates[33] * popestimates[19]^2 * popestimates[21]^2

# Find sigma
sigma_CD4 <- popestimates[34]
sigma_decay <- popestimates[35]
sigma_rebound <- popestimates[36]
# Summary of the results
names <- c("$P_1$","$\\lambda_1$","$P_2$","$\\lambda_2$", "$\\beta_1$", "$\\beta_{2}$", "$\\beta_3$","$\\beta_4$", "\\beta_{50}", "$\\beta_{51}$", "$\\beta_{52}$", "$\\beta_{53}$","$\\alpha_0$","$\\alpha_1$")
z <- popestimates[1:14]/popse
p <- pnorm(abs(z), lower.tail = FALSE)*2
summary <- as.data.frame(cbind(names, round(popestimates[1:14], 3), round(popse, 3), round(z, 3), round(p, 3)))
rownames(summary) <- NULL
colnames(summary) <- c("Parameter", "Estimates", "SE", "z-value", "p-value")
summary
popestimates


latex_table<-xtable(summary, type = "latex",align=c("cccccc"))
digits(latex_table)<-c(0,0,3,3,3,3)
print(latex_table, file = "latex_table_Ch4_est2.tex",include.rownames=FALSE,sanitize.text.function = function(x){x},hline.after = c(-1,0,13))


# Latex Table for correlation 
cortable_names <- c("Random Effects","$b_{3i}$","$b_{4i}$", "$\\tau_{1i}$","$\\tau_{3i}$")
cortable <- data.frame(matrix(nrow=4, ncol=length(cortable_names)))
colnames(cortable) <- cortable_names
cortable$`Random Effects` <- c("$b_{3i}$","$b_{4i}$", "$\\tau_{1i}$","$\\tau_{3i}$")
cortable[1,2] <- cortable[2,3]<- cortable[3,4]<- cortable[4,5] <- 1
cortable[2,2] <- cortable[1,3] <- popestimates[28]
cortable[3,2] <- cortable[1,4] <-  popestimates[29]
cortable[4,2] <- cortable[1,5] <- popestimates[30]
cortable[3,3] <- cortable[2,4] <-  popestimates[31]
cortable[4,3] <- cortable[2,5] <-  popestimates[32]
cortable[4,4] <- cortable[3,5] <- popestimates[33]

latex_table<-xtable(cortable, type = "latex",align=c("cccccc"))
digits(latex_table)<-c(0,0,3,3,3,3)
print(latex_table, file = "latex_table_Ch4_cor.tex",include.rownames=FALSE,sanitize.text.function = function(x){x},hline.after = c(-1,0,4))


####### Individual fits ========
data <- read.csv("/Users/gsiii/Library/CloudStorage/OneDrive-UBC/PhD/Unified/Code/threeobssqrt.csv")
# get individual parameters
indpar <- getEstimatedIndividualParameters()$saem

# ID = 1
indfit <- function(ID){
  pat1 <- data %>% filter(PATIENT == ID) 
  decay <- pat1 %>% filter(name == "decay") %>% mutate(observation = ifelse(censor == 1, NA, observation))
  rebound <- pat1 %>% filter(name == "rebound")%>% mutate(observation = ifelse(censor == 1, NA, observation))
  CD4 <- pat1 %>% filter(name == "CD4_V") %>% 
    mutate(time=time_raw - Ti)
  
  indparest1 <- indpar %>% filter(id == ID)
  P1 <- indparest1$p1
  lambda1 <- indparest1$b1
  P2 <- indparest1$p2
  lambda2 <- indparest1$b2
  beta1 <- indparest1$beta1
  beta2 <- indparest1$beta2
  beta3 <- indparest1$beta3
  beta4 <- indparest1$beta4
  beta50 <- indparest1$beta50
  beta51 <- indparest1$beta51
  beta52 <- indparest1$beta52
  beta53 <- indparest1$beta53
  alpha0 <- indparest1$alpha0
  alpha1 <- indparest1$alpha1
  
  decayfun <- function(x){
    log10( exp(P1-lambda1*x) + exp(P2-lambda2 * x))
  }
  cd4fun <- function(x){
    alpha0 + alpha1 * x
  }
  
  decay_censor <- decay[which(decay$censor==1),]$time
  line_df <- data.frame(
    x = decay_censor,
    y = 0,
    xend = decay_censor,
    yend = log10(40)
  )
  (p1 <- ggplot(decay, aes(x = time_raw, y = observation)) + geom_point(size=2) + 
      stat_function(fun = decayfun) + geom_segment(data = line_df, 
                                                   mapping = aes(x=x, y=y, xend=xend, yend=yend), 
                                                   inherit.aes = FALSE, linewidth = 1, color = "grey")+
      xlab("Time (in months)") +  scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)")) + theme_classic() + theme(text = element_text(size = 14)))
  
  rebound <- rebound %>%
    mutate(CD4 = alpha0 + alpha1 * (time_raw-Ti)) %>%
    mutate(beta5nonpar = beta50* CD4 + beta51 + beta52 * z1 + beta53 * z2 ) %>%
    mutate(fitted = beta1 * (time_raw-Ti) / ( (time_raw-Ti) + exp( beta2 -beta3 * (time_raw-Ti) )) + beta4/(1+exp(beta5nonpar*(time_raw-Ti)))) %>%
    mutate(time=time_raw - Ti)
  # newdata <- data.frame(time = seq(min(rebound$time), max(rebound$time),0.1))
  # basis <- ns(newdata$time, df=2, intercept=F)
  # newdata <- data.frame(newdata, basis)
  # names(newdata)[c(ncol(newdata)-1,ncol(newdata))] <- c("z1","z2")
  # newdata <- newdata %>% 
  #   mutate(CD4 = alpha0 + alpha1 * time) %>% 
  #   mutate(beta5nonpar = beta50* CD4 + beta51 + beta52 * z1 + beta53 * z2 ) %>% 
  #   mutate(fitted = beta1 * time / (time + exp( beta2 -beta3 * time )) + beta4/(1+exp(beta5nonpar*(time_raw-Ti)))) 

  rebound_censor <- rebound[which(rebound$censor==1),]$time
  line_df <- data.frame(
    x = ifelse(length(rebound_censor)==0, 0, rebound_censor),
    y = ifelse(length(rebound_censor)==0, min(rebound$observation), 0),
    xend = ifelse(length(rebound_censor)==0, 0, rebound_censor),
    yend = ifelse(length(rebound_censor)==0, min(rebound$observation), log10(40))
  )
  
  (p2 <- ggplot(rebound, aes(x = time, y = observation)) + geom_point(size=2) + 
      geom_segment(data = line_df, 
                   mapping = aes(x=x, y=y, xend=xend, yend=yend), 
                   inherit.aes = FALSE, linewidth = 1, color = "grey")+
      geom_xspline(aes(x =time, y=fitted),  spline_shape = 1, color = "black", linewidth=0.5)+
      xlab("Time (in months)") + scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)")) + theme_classic() + theme(text = element_text(size = 14)))
  (p3 <- ggplot(CD4, aes(x = time, y = observation)) + geom_point(size=2) + 
      stat_function(fun = cd4fun) + 
      xlab("Time (in months)") + scale_y_continuous(bquote("CD4 cell count (in cells /" ~ mm^3 ~")")) + theme_classic() + theme(text = element_text(size = 14)))
  (p4 <- ggplot(rebound, aes(x = time, y = observation)) + geom_point(size=2) +
      geom_segment(data = line_df, 
                   mapping = aes(x=x, y=y, xend=xend, yend=yend), 
                   inherit.aes = FALSE, linewidth = 1, color = "grey")+
      geom_line(aes(x =time, y=fitted),se=FALSE, color = "black", linewidth=0.5)+
      xlab("Time (in months)") + scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)")) + theme_classic() + theme(text = element_text(size = 14)))
  
  #grid.arrange(p1, p2, p3,p4, ncol=2)
  grid.arrange(p1, p2, p3, ncol=3)
}
indfit(3)
indfit(8)
indfit(19)
indfit(2)



