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
library(GLDEX)
library(lme4)
library(xtable)
library(berryFunctions)
initializeLixoftConnectors(software="monolix")
setwd("/Users/gsiii/Library/CloudStorage/OneDrive-UBC/PhD/Joint")
# setwd("/Users/visstar/Library/CloudStorage/OneDrive-UBC/PhD/Joint")


num_sim <- 100  # number of simulation
num_patient <- 200  # number of patients
censor_value <- log10(25)  # censor value 
num_boots = 100


#Set true values for p1, b1, p2, sigma1, and covariance matrix D for viral loads before ART interruption
#Follow model 1: y_ij=log10(exp(p_1i-lambda_1i*t_ij)+exp(p_2i-lambda_2i*t_ij))+e_ij
#Random effects: p_1i=p1+b_1i,p_2i=p1+b_2i,lambda_1i=lambda1+b_3i,lambda_2i=lambda2+b_4i,
#where e_ij~N(0,sigma_v^2) and b_i~N(0,D)
P1=11
lambda1 = 0.3
P2 = 2.8
beta1 = 1.8
beta2=8
beta3 = 0.4
beta4 = 1.3

gamma0 <- 8
gamma1 <- -2.2

sigma1 <- 0.1

A <- diag(c(2, 0.01, 0.05))
B <- diag(c(0.05, 0, 0.01, 0))
sigma2 <- 0.3

parnames <- c("gamma0","gamma1")
parnames_viral <- c("P1", "lambda1", "P2", "beta1", "beta2", "beta3", "beta4")

# Choose time from true_value dataset. 
time0 <- c(0, 8, 14, 28, 56, 65, 84, 87, 90 , 93, 96,  99, 102, 104, 108, 111, 116, 120, 125,140) 


#Create a data frame for simulation estimates.
ncol <- length(parnames)
ncol_viral <- length(parnames_viral)

estimates_naive <- estimates_empirical <- estimates_new <- data.frame(matrix(nrow=num_sim,ncol=ncol))
colnames(estimates_naive) <- colnames(estimates_empirical) <- colnames(estimates_new) <- parnames
estimates_new_viral <- data.frame(matrix(nrow=num_sim,ncol=ncol_viral))
colnames(estimates_new_viral) <- parnames_viral

SE_new <- SE_empirical <- SE_naive <- data.frame(matrix(nrow=num_sim,ncol=ncol))
colnames(SE_naive) <- colnames(SE_empirical) <- colnames(SE_new) <- parnames
SE_new_viral <- data.frame(matrix(nrow=num_sim,ncol=ncol_viral))
colnames(SE_new_viral) <- parnames_viral


bootstrap_nonpar_new<- bootstrap_par_new <- data.frame(matrix(nrow=num_boots,ncol=ncol))
colnames(bootstrap_par_new) <- colnames(bootstrap_nonpar_new) <- parnames
bootstrap_par_new_viral <- bootstrap_nonpar_new_viral <- data.frame(matrix(nrow=num_boots,ncol =ncol_viral))
colnames(bootstrap_par_new_viral) <- colnames(bootstrap_par_new_viral) <- parnames_viral


bootstrap_par_SE_new <- bootstrap_nonpar_SE_new <- data.frame(matrix(nrow = num_sim,ncol = ncol))
colnames(bootstrap_par_SE_new) <- colnames(bootstrap_nonpar_SE_new) <- parnames
bootstrap_par_SE_new_viral<- bootstrap_nonpar_SE_new_viral<- data.frame(matrix(nrow=num_sim,ncol = ncol_viral))
colnames(bootstrap_par_SE_new_viral) <- colnames(bootstrap_nonpar_SE_new_viral) <- parnames_viral


true_value <- c(gamma0, gamma1)
true_value_viral <- c(P1, lambda1, P2, beta1, beta2, beta3, beta4)

i=1
for (i in 1:10){
  print(glue::glue("i = ",i))
  
  # Create dataset for viral load: 1/2*num_patient patients with time1 and 1/2*num_patient patients with time2
  PATIENT <- rep(1:num_patient, each = length(time0))
  time <- rep(time0, num_patient)
  decay <- rep(c(rep(1, 7), rep(0, 13)),num_patient)
  rebound <- rep(c(rep(0, 7), rep(1, 13)),num_patient)
  Ti <- rep(84, num_patient*20)
  data0 <- as.data.frame(cbind(PATIENT, time, decay, rebound, Ti))
  
  # simulate a_i, b_i, and e_ij
  ai_sim <- mvrnorm(num_patient, rep(0,dim(A)[1]), A)
  colnames(ai_sim) <- c("a1_sim","a2_sim","a3_sim")
  bi_sim <- mvrnorm(num_patient, rep(0,dim(B)[1]), B)
  colnames(bi_sim) <- c("b1_sim", "b2_sim" ,"b3_sim","b4_sim")
  g0_sim = rnorm(num_patient, 0, sigma2)
  ranef_sim <- cbind(PATIENT=c(1:num_patient), ai_sim, bi_sim, g0_sim)
  data <- merge(data0, ranef_sim, by="PATIENT", all = TRUE)
  data <- data %>% 
    mutate(e = rnorm(nrow(data), 0, sigma1)) %>% 
    mutate(w_star = log10(exp(P1+a1_sim-(lambda1+a2_sim)*time)+exp(P2+a3_sim))*decay + ((beta1+b1_sim) * (time-Ti) / ((time-Ti)  + exp(beta2 + b2_sim - (beta3+b3_sim) * (time-Ti))) + (beta4 + b4_sim))*rebound) %>% 
    mutate(log10rna = w_star + e) %>% 
    mutate(censor = ifelse(log10rna < censor_value, 1, 0)) %>% 
    mutate(log10rna = ifelse(log10rna < censor_value, censor_value, log10rna)) %>% 
    mutate(limit = 0) %>% 
    mutate(xb = gamma0 + g0_sim + gamma1 * w_star) %>% 
    mutate(p = 1/(1+exp(-xb))) %>% 
    mutate(CD4 = rbinom(n = nrow(data), size = 1, p = p)) 
  data <- data %>% 
    dplyr::select(c(PATIENT, time, decay, rebound, Ti, censor, limit, log10rna, CD4))
  data_new <- pivot_longer(data, cols = c(8:9), names_to = "name", values_to = "observation")
  data_new <- data_new %>% 
    mutate(censor = ifelse(name == "CD4", 0, censor))
  data_empirical <- data_new %>% 
    mutate(observation = ifelse(censor == 1 & name=="log10rna", 0.5*observation, observation))
  data_naive <- data %>% 
    mutate(log10rna = ifelse(censor == 1, 0.5*log10rna, log10rna))
  summary(as.factor(data$censor))
  summary(as.factor(data$CD4))
  ggplot(data[PATIENT%in%c(1:5),], aes(x=time, y=log10rna)) +
    geom_point(aes(fill=factor(censor)),size=2, shape=21, stroke=0) +
    geom_line(aes(group=PATIENT)) +
    scale_x_continuous("Day") +
    scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)"))+
    scale_fill_manual(values=c("black","red"),labels = c("Observed data", "Censored data"))+
    scale_color_manual(values=c("steelblue","black"),labels = c("Following ART interruption", "During ART"))+
    labs(color = "ART status", fill = "Data type")+ggtitle("Plot for all observations")+
    theme_classic()+theme(text=element_text(size=14),panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"))
  
  # Save the simulated data locally
  write.table(data_new, "Code/Simulation/sim_data_new.txt", sep = ",", quote = FALSE, row.names = FALSE)
  write.table(data_naive, "Code/Simulation/sim_data_naive.txt", sep = ",", quote = FALSE, row.names = FALSE)
  write.table(data_empirical, "Code/Simulation/sim_data_empirical.txt", sep = ",", quote = FALSE, row.names = FALSE)
  
  # Our model 
  # Read in data for Monolix
  data_joint = list(dataFile = paste0('Code/Simulation/sim_data_new.txt'),
                    headerTypes =c("id","time","regressor", "regressor", "regressor","cens","limit","obsid","observation"),
                    observationTypes =list(CD4 = "categorical", log10rna = "continuous"))
  modelFile = paste0('Code/model_joint_full6.txt')
  # create a new project by setting a data set and a structural model
  newProject(data = data_joint, modelFile = modelFile)
  
  #getObservationInformation()
  setErrorModel(ylog10rna = "constant")
  setObservationDistribution(ylog10rna= "normal")
  # covariates on beta50
  # getCovariateInformation()
  # setCovariateModel(gamma0 = c(trtarm=TRUE))
  # set tasks in scenario
  scenario <- getScenario()
  scenario$tasks = c(populationParameterEstimation = T,
                     conditionalModeEstimation = T,
                     conditionalDistributionSampling = T,
                     standardErrorEstimation=T,
                     logLikelihoodEstimation=T)
  scenario$linearization = FALSE
  setIndividualParameterVariability(p1 = TRUE, b1 = TRUE, p2 = TRUE, beta1 = TRUE, beta2 = FALSE, beta3 = TRUE, beta4 = FALSE, gamma0 = TRUE, gamma1 = FALSE)
  setPopulationParameterEstimationSettings(nbexploratoryiterations=1000, nbsmoothingiterations=500)
  setIndividualParameterDistribution(p1="normal",b1="normal",p2="normal",beta1="normal",beta2="normal",beta3="normal",beta4="normal", gamma0="normal", gamma1="normal")
  setPopulationParameterInformation(p1_pop = list(initialValue = 11), 
                                    b1_pop = list(initialValue = 0.3),
                                    p2_pop = list(initialValue = 2.8),
                                    beta1_pop = list(initialValue = 1.8),
                                    beta2_pop = list(initialValue = 8),
                                    beta3_pop = list(initialValue = 0.4),
                                    beta4_pop = list(initialValue = 1.3),
                                    gamma0_pop = list(initialValue = 8), 
                                    gamma1_pop  = list(initialValue = -2.5))
  # getVariabilityLevels()
  # setCorrelationBlocks(ID = list(c("b1", "p1", "b2", "p2","beta1","beta2","beta3")))
  # run the estimation
  setScenario(scenario)
  runScenario()
  # store the estimates
  estimates_new[i,] <- getEstimatedPopulationParameters()[8:9]
  SE_new[i,] <- as.numeric(getEstimatedStandardErrors()[["stochasticApproximation"]][["se"]][8:9])
  estimates_new_viral[i,] <- getEstimatedPopulationParameters()[1:7]
  SE_new_viral[i,] <- as.numeric(getEstimatedStandardErrors()[["stochasticApproximation"]][["se"]][1:7])
  
  # Naive model 
  data_naive <- read.table("Code/Simulation/sim_data_naive.txt", sep=",", header = TRUE)
  model_naive <- glmmPQL(CD4 ~ log10rna, random = ~ 1|PATIENT,family=binomial, data=data_naive)
  estimates_naive[i,] <- summary(model_naive)$coefficients$fixed
  SE_naive[i,] <- summary(model_naive)$tTable[,2]
  
  # Empirical model
  # Read in data for Monolix
  data_empirical = list(dataFile = paste0('Code/Simulation/sim_data_empirical.txt'),
                        headerTypes =c("id","time","ignore","ignore","ignore","ignore","ignore","obsid","observation"),
                        observationTypes =list(CD4 = "categorical", log10rna = "continuous"))
  modelFile = paste0('Code/model_joint_empirical.txt')
  # create a new project by setting a data set and a structural model
  newProject(data = data_empirical, modelFile = modelFile)
  
  #getObservationInformation()
  setErrorModel(list(ylog10rna = "constant"))
  setObservationDistribution(ylog10rna= "normal")
  # covariates on beta50
  # getCovariateInformation()
  # setCovariateModel(gamma0 = c(trtarm=TRUE))
  # set tasks in scenario
  scenario <- getScenario()
  scenario$tasks = c(populationParameterEstimation = T,
                     conditionalModeEstimation = T,
                     conditionalDistributionSampling = T,
                     standardErrorEstimation=T,
                     logLikelihoodEstimation=T)
  scenario$linearization = FALSE
  setIndividualParameterVariability(alpha0 = TRUE, alpha1 = TRUE, alpha2 = TRUE, gamma0 = TRUE, gamma1 = FALSE)
  setPopulationParameterEstimationSettings(nbexploratoryiterations=1000)
  setIndividualParameterDistribution(alpha0="normal",alpha1="normal",alpha2="normal", gamma0="normal", gamma1="normal")
  setPopulationParameterInformation(gamma0_pop = list(initialValue = 8), 
                                    gamma1_pop  = list(initialValue = -2.5))
  # getVariabilityLevels()
  # setCorrelationBlocks(ID = list(c("b1", "p1", "b2", "p2","beta1","beta2","beta3")))
  # run the estimation
  setScenario(scenario)
  runScenario()
  # store the estimates
  estimates_empirical[i,] <- getEstimatedPopulationParameters()[4:5]
  SE_empirical[i,] <- as.numeric(getEstimatedStandardErrors()[["stochasticApproximation"]][["se"]][4:5])
  
  while(any(is.na(c(SE_empirical[i, ], SE_new[i, ], SE_naive[i, ], SE_new_viral[i, ])))){
    # simulate a_i, b_i, and e_ij
    ai_sim <- mvrnorm(num_patient, rep(0,dim(A)[1]), A)
    colnames(ai_sim) <- c("a1_sim","a2_sim","a3_sim")
    bi_sim <- mvrnorm(num_patient, rep(0,dim(B)[1]), B)
    colnames(bi_sim) <- c("b1_sim", "b2_sim" ,"b3_sim","b4_sim")
    g0_sim = rnorm(num_patient, 0, sigma2)
    ranef_sim <- cbind(PATIENT=c(1:num_patient), ai_sim, bi_sim, g0_sim)
    data <- merge(data0, ranef_sim, by="PATIENT", all = TRUE)
    data <- data %>% 
      mutate(e = rnorm(nrow(data), 0, sigma1)) %>% 
      mutate(w_star = log10(exp(P1+a1_sim-(lambda1+a2_sim)*time)+exp(P2+a3_sim))*decay + ((beta1+b1_sim) * (time-Ti) / ((time-Ti)  + exp(beta2 + b2_sim - (beta3+b3_sim) * (time-Ti))) + beta4 + b4_sim)*rebound) %>% 
      mutate(log10rna = w_star + e) %>% 
      mutate(censor = ifelse(log10rna < censor_value, 1, 0)) %>% 
      mutate(log10rna = ifelse(log10rna < censor_value, censor_value, log10rna)) %>% 
      mutate(limit = 0) %>% 
      mutate(xb = gamma0 + g0_sim + gamma1 * w_star) %>% 
      mutate(p = 1/(1+exp(-xb))) %>% 
      mutate(CD4 = rbinom(n = nrow(data), size = 1, p = p)) 
    data <- data %>% 
      dplyr::select(c(PATIENT, time, decay, rebound, Ti, censor, limit, log10rna, CD4))
    data_new <- pivot_longer(data, cols = c(8:9), names_to = "name", values_to = "observation")
    data_new <- data_new %>% 
      mutate(censor = ifelse(name == "CD4", 0, censor))
    data_empirical <- data_new %>% 
      mutate(observation = ifelse(censor == 1 & name=="log10rna", 0.5*observation, observation))
    data_naive <- data %>% 
      mutate(log10rna = ifelse(censor == 1, 0.5*log10rna, log10rna))

    # Save the simulated data locally
    write.table(data_new, "Code/Simulation/sim_data_new.txt", sep = ",", quote = FALSE, row.names = FALSE)
    write.table(data_naive, "Code/Simulation/sim_data_naive.txt", sep = ",", quote = FALSE, row.names = FALSE)
    write.table(data_empirical, "Code/Simulation/sim_data_empirical.txt", sep = ",", quote = FALSE, row.names = FALSE)
    
    # Our model 
    # Read in data for Monolix
    data_joint = list(dataFile = paste0('Code/Simulation/sim_data_new.txt'),
                      headerTypes =c("id","time","regressor", "regressor", "regressor","cens","limit","obsid","observation"),
                      observationTypes =list(CD4 = "categorical", log10rna = "continuous"))
    modelFile = paste0('Code/model_joint_full6.txt')
    # create a new project by setting a data set and a structural model
    newProject(data = data_joint, modelFile = modelFile)
    
    #getObservationInformation()
    setErrorModel(ylog10rna = "constant")
    setObservationDistribution(ylog10rna= "normal")
    # covariates on beta50
    # getCovariateInformation()
    # setCovariateModel(gamma0 = c(trtarm=TRUE))
    # set tasks in scenario
    scenario <- getScenario()
    scenario$tasks = c(populationParameterEstimation = T,
                       conditionalModeEstimation = T,
                       conditionalDistributionSampling = T,
                       standardErrorEstimation=T,
                       logLikelihoodEstimation=T)
    scenario$linearization = FALSE
    setIndividualParameterVariability(p1 = TRUE, b1 = TRUE, p2 = TRUE, beta1 = TRUE, beta2 = FALSE, beta3 = TRUE, beta4 = FALSE, gamma0 = TRUE, gamma1 = FALSE)
    setPopulationParameterEstimationSettings(nbexploratoryiterations=1000, nbsmoothingiterations=500)
    setIndividualParameterDistribution(p1="normal",b1="normal",p2="normal",beta1="normal",beta2="normal",beta3="normal",beta4="normal", gamma0="normal", gamma1="normal")
    setPopulationParameterInformation(p1_pop = list(initialValue = 11), 
                                      b1_pop = list(initialValue = 0.3),
                                      p2_pop = list(initialValue = 2.8),
                                      beta1_pop = list(initialValue = 1.8),
                                      beta2_pop = list(initialValue = 8),
                                      beta3_pop = list(initialValue = 0.4),
                                      beta4_pop = list(initialValue = 1.3),
                                      gamma0_pop = list(initialValue = 8), 
                                      gamma1_pop  = list(initialValue = -2.5))
    # getVariabilityLevels()
    # setCorrelationBlocks(ID = list(c("b1", "p1", "b2", "p2","beta1","beta2","beta3")))
    # run the estimation
    setScenario(scenario)
    runScenario()
    # store the estimates
    estimates_new[i,] <- getEstimatedPopulationParameters()[8:9]
    SE_new[i,] <- as.numeric(getEstimatedStandardErrors()[["stochasticApproximation"]][["se"]][8:9])
    estimates_new_viral[i,] <- getEstimatedPopulationParameters()[1:7]
    SE_new_viral[i,] <- as.numeric(getEstimatedStandardErrors()[["stochasticApproximation"]][["se"]][1:7])
    
    # Naive model 
    data_naive <- read.table("Code/Simulation/sim_data_naive.txt", sep=",", header = TRUE)
    model_naive <- glmmPQL(CD4 ~ log10rna, random = ~ 1|PATIENT,family=binomial, data=data_naive)
    estimates_naive[i,] <- summary(model_naive)$coefficients$fixed
    SE_naive[i,] <- summary(model_naive)$tTable[,2]
    
    # Empirical model
    # Read in data for Monolix
    data_empirical = list(dataFile = paste0('Code/Simulation/sim_data_empirical.txt'),
                          headerTypes =c("id","time","ignore","ignore","ignore","ignore","ignore","obsid","observation"),
                          observationTypes =list(CD4 = "categorical", log10rna = "continuous"))
    modelFile = paste0('Code/model_joint_empirical.txt')
    # create a new project by setting a data set and a structural model
    newProject(data = data_empirical, modelFile = modelFile)
    
    #getObservationInformation()
    setErrorModel(list(ylog10rna = "constant"))
    setObservationDistribution(ylog10rna= "normal")
    # covariates on beta50
    # getCovariateInformation()
    # setCovariateModel(gamma0 = c(trtarm=TRUE))
    # set tasks in scenario
    scenario <- getScenario()
    scenario$tasks = c(populationParameterEstimation = T,
                       conditionalModeEstimation = T,
                       conditionalDistributionSampling = T,
                       standardErrorEstimation=T,
                       logLikelihoodEstimation=T)
    scenario$linearization = FALSE
    setIndividualParameterVariability(alpha0 = TRUE, alpha1 = TRUE, alpha2 = TRUE, gamma0 = TRUE, gamma1 = FALSE)
    setPopulationParameterEstimationSettings(nbexploratoryiterations=1000)
    setIndividualParameterDistribution(alpha0="normal",alpha1="normal",alpha2="normal", gamma0="normal", gamma1="normal")
    setPopulationParameterInformation(gamma0_pop = list(initialValue = 8), 
                                      gamma1_pop  = list(initialValue = -2.5))
    # getVariabilityLevels()
    # setCorrelationBlocks(ID = list(c("b1", "p1", "b2", "p2","beta1","beta2","beta3")))
    # run the estimation
    setScenario(scenario)
    runScenario()
    # store the estimates
    estimates_empirical[i,] <- getEstimatedPopulationParameters()[4:5]
    SE_empirical[i,] <- as.numeric(getEstimatedStandardErrors()[["stochasticApproximation"]][["se"]][4:5])
  }
  
  for (j in 1:num_boots){
    print(glue::glue("i = ",i," j = ",j))
    
    # Parametric Bootstrap
    # simulate a_i, b_i, and e_ij
    ai_sim <- mvrnorm(num_patient, rep(0,dim(A)[1]), A)
    colnames(ai_sim) <- c("a1_sim","a2_sim","a3_sim")
    bi_sim <- mvrnorm(num_patient, rep(0,dim(B)[1]), B)
    colnames(bi_sim) <- c("b1_sim", "b2_sim" ,"b3_sim","b4_sim")
    g0_sim = rnorm(num_patient, 0, sigma2)
    ranef_sim <- cbind(PATIENT=c(1:num_patient), ai_sim, bi_sim, g0_sim)
    data <- merge(data0, ranef_sim, by="PATIENT", all = TRUE)
    data <- data %>% 
      mutate(e = rnorm(nrow(data), 0, sigma1)) %>% 
      mutate(w_star = log10(exp(estimates_new_viral[i, 1]+a1_sim-(estimates_new_viral[i, 2]+a2_sim)*time)+exp(estimates_new_viral[i, 3]+a3_sim))*decay + ((estimates_new_viral[i, 4]+b1_sim) * (time-Ti) / ((time-Ti)  + exp(estimates_new_viral[i, 5] + b2_sim - (estimates_new_viral[i, 6]+b3_sim) * (time-Ti))) + estimates_new_viral[i, 7] + b4_sim)*rebound) %>% 
      mutate(log10rna = w_star + e) %>% 
      mutate(censor = ifelse(log10rna < censor_value, 1, 0)) %>% 
      mutate(log10rna = ifelse(log10rna < censor_value, censor_value, log10rna)) %>% 
      mutate(limit = 0) %>% 
      mutate(xb = estimates_new[i, 1] + g0_sim + estimates_new[i, 2] * w_star) %>% 
      mutate(p = 1/(1+exp(-xb))) %>% 
      mutate(CD4 = rbinom(n = nrow(data), size = 1, p = p)) 
    data <- data %>% 
      dplyr::select(c(PATIENT, time, decay, rebound, Ti, censor, limit, log10rna, CD4))
    data_new_boot <- pivot_longer(data, cols = c(8:9), names_to = "name", values_to = "observation")
    data_new_boot <- data_new_boot %>% 
      mutate(censor = ifelse(name == "CD4", 0, censor))
    ggplot(data[PATIENT%in%c(1:5),], aes(x=time, y=log10rna)) +
      geom_point(aes(fill=factor(censor)),size=2, shape=21, stroke=0) +
      geom_line(aes(group=PATIENT)) +
      scale_x_continuous("Day") +
      scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)"))+
      scale_fill_manual(values=c("black","red"),labels = c("Observed data", "Censored data"))+
      scale_color_manual(values=c("steelblue","black"),labels = c("Following ART interruption", "During ART"))+
      labs(color = "ART status", fill = "Data type")+ggtitle("Plot for all observations")+
      theme_classic()+theme(text=element_text(size=14),panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"))
    
    # Save the simulated data locally
    write.table(data_new_boot, "Code/Simulation/sim_data_new_boot.txt", sep = ",", quote = FALSE, row.names = FALSE)
    
    # Our model
    # Read in data for Monolix
    data_new_boot = list(dataFile = paste0('Code/Simulation/sim_data_new_boot.txt'),
                         headerTypes =c("id","time","regressor", "regressor", "regressor","cens","limit","obsid","observation"),
                         observationTypes =list(CD4 = "categorical", log10rna = "continuous"))
    modelFile = paste0('Code/model_joint_full6.txt')
    # create a new project by setting a data set and a structural model
    newProject(data = data_new_boot, modelFile = modelFile)
    
    #getObservationInformation()
    setErrorModel(ylog10rna = "constant")
    setObservationDistribution(ylog10rna= "normal")
    # covariates on beta50
    # getCovariateInformation()
    # setCovariateModel(gamma0 = c(trtarm=TRUE))
    # set tasks in scenario
    scenario <- getScenario()
    scenario$tasks = c(populationParameterEstimation = T,
                       conditionalModeEstimation = T,
                       conditionalDistributionSampling = T,
                       standardErrorEstimation=T,
                       logLikelihoodEstimation=T)
    scenario$linearization = FALSE
    setIndividualParameterVariability(p1 = TRUE, b1 = TRUE, p2 = TRUE, beta1 = TRUE, beta2 = FALSE, beta3 = TRUE, beta4 = FALSE, gamma0 = TRUE, gamma1 = FALSE)
    setPopulationParameterEstimationSettings(nbexploratoryiterations=1000, nbsmoothingiterations=500)
    setIndividualParameterDistribution(p1="normal",b1="normal",p2="normal",beta1="normal",beta2="normal",beta3="normal",beta4="normal", gamma0="normal", gamma1="normal")
    setPopulationParameterInformation(p1_pop = list(initialValue = 11),
                                      b1_pop = list(initialValue = 0.3),
                                      p2_pop = list(initialValue = 2.8),
                                      beta1_pop = list(initialValue = 1.8),
                                      beta2_pop = list(initialValue = 8),
                                      beta3_pop = list(initialValue = 0.4),
                                      beta4_pop = list(initialValue = 1.3),
                                      gamma0_pop = list(initialValue = 8),
                                      gamma1_pop  = list(initialValue = -2.5))
    # getVariabilityLevels()
    # setCorrelationBlocks(ID = list(c("b1", "p1", "b2", "p2","beta1","beta2","beta3")))
    # run the estimation
    setScenario(scenario)
    runScenario()
    # store the estimates
    bootstrap_par_new[j,] <- getEstimatedPopulationParameters()[8:9]
    bootstrap_par_new_viral[j,] <- getEstimatedPopulationParameters()[1:7]
    
    # Nonparametric bootstrap
    pat.select <- sample(c(1:num_patient), size = num_patient, replace = TRUE)
    data_new_boot <- data_new[which(data_new$PATIENT == pat.select[1]), ]
    for (k in 2:num_patient){
      data_new_boot <- rbind(data_new_boot, data_new[which(data_new$PATIENT == pat.select[k]), ])
    }
    # Save the simulated data locally
    write.table(data_new_boot, "Code/Simulation/sim_data_new_boot.txt", sep = ",", quote = FALSE, row.names = FALSE)
    
    # Our model
    # Read in data for Monolix
    data_new_boot = list(dataFile = paste0('Code/Simulation/sim_data_new_boot.txt'),
                         headerTypes =c("id","time","regressor", "regressor", "regressor","cens","limit","obsid","observation"),
                         observationTypes =list(CD4 = "categorical", log10rna = "continuous"))
    modelFile = paste0('Code/model_joint_full6.txt')
    # create a new project by setting a data set and a structural model
    newProject(data = data_new_boot, modelFile = modelFile)
    
    #getObservationInformation()
    setErrorModel(ylog10rna = "constant")
    setObservationDistribution(ylog10rna= "normal")
    # covariates on beta50
    # getCovariateInformation()
    # setCovariateModel(gamma0 = c(trtarm=TRUE))
    # set tasks in scenario
    scenario <- getScenario()
    scenario$tasks = c(populationParameterEstimation = T,
                       conditionalModeEstimation = T,
                       conditionalDistributionSampling = T,
                       standardErrorEstimation=T,
                       logLikelihoodEstimation=T)
    scenario$linearization = FALSE
    setIndividualParameterVariability(p1 = TRUE, b1 = TRUE, p2 = TRUE, beta1 = TRUE, beta2 = FALSE, beta3 = TRUE, beta4 = FALSE, gamma0 = TRUE, gamma1 = FALSE)
    setPopulationParameterEstimationSettings(nbexploratoryiterations=1000, nbsmoothingiterations=500)
    setIndividualParameterDistribution(p1="normal",b1="normal",p2="normal",beta1="normal",beta2="normal",beta3="normal",beta4="normal", gamma0="normal", gamma1="normal")
    setPopulationParameterInformation(p1_pop = list(initialValue = 11),
                                      b1_pop = list(initialValue = 0.3),
                                      p2_pop = list(initialValue = 2.8),
                                      beta1_pop = list(initialValue = 1.8),
                                      beta2_pop = list(initialValue = 8),
                                      beta3_pop = list(initialValue = 0.4),
                                      beta4_pop = list(initialValue = 1.3),
                                      gamma0_pop = list(initialValue = 8),
                                      gamma1_pop  = list(initialValue = -2.5))
    # getVariabilityLevels()
    # setCorrelationBlocks(ID = list(c("b1", "p1", "b2", "p2","beta1","beta2","beta3")))
    # run the estimation
    setScenario(scenario)
    runScenario()
    # store the estimates
    bootstrap_nonpar_new[j,] <- getEstimatedPopulationParameters()[8:9]
    bootstrap_nonpar_new_viral[j,] <- getEstimatedPopulationParameters()[1:7]
  }
  bootstrap_par_SE_new[i,] <- sapply(bootstrap_par_new,sd)
  bootstrap_par_SE_new_viral[i,] <- sapply(bootstrap_par_new_viral,sd)
  bootstrap_nonpar_SE_new[i,] <- sapply(bootstrap_nonpar_new,sd)
  bootstrap_nonpar_SE_new_viral[i,] <- sapply(bootstrap_nonpar_new_viral,sd)
  
  saveRDS(estimates_new, "Code/Simulation/estimates_new_boot.RDS")
  saveRDS(estimates_new_viral, "Code/Simulation/estimates_new_viral_boot.RDS")
  saveRDS(estimates_naive, "Code/Simulation/estimates_naive_boot.RDS")
  saveRDS(estimates_empirical, "Code/Simulation/estimates_empirical_boot.RDS")
  saveRDS(SE_new, "Code/Simulation/SE_new_boot.RDS")
  saveRDS(SE_new_viral, "Code/Simulation/SE_new_viral_boot.RDS")
  saveRDS(SE_naive, "Code/Simulation/SE_naive_boot.RDS")
  saveRDS(SE_empirical, "Code/Simulation/SE_empirical_boot.RDS")
  saveRDS(bootstrap_par_SE_new,"Code/Simulation/bootstrap_par_SE_new_boot.RDS")
  saveRDS(bootstrap_par_SE_new_viral,"Code/Simulation/bootstrap_par_SE_new_viral_boot.RDS")
  saveRDS(bootstrap_nonpar_SE_new,"Code/Simulation/bootstrap_nonpar_SE_new_boot.RDS")
  saveRDS(bootstrap_nonpar_SE_new_viral,"Code/Simulation/bootstrap_nonpar_SE_new_viral_boot.RDS")
}  



estimates_new <- readRDS("Code/Simulation/estimates_new_boot.RDS")
SE_new <- readRDS("Code/Simulation/SE_new_boot.RDS")
bootstrap_par_SE_new  <- readRDS("Code/Simulation/bootstrap_par_SE_new_boot.RDS")
bootstrap_nonpar_SE_new  <- readRDS("Code/Simulation/bootstrap_nonpar_SE_new_boot.RDS")
estimates_new_viral <- readRDS("Code/Simulation/estimates_new_viral_boot.RDS")
SE_new_viral <- readRDS("Code/Simulation/SE_new_viral_boot.RDS")
bootstrap_par_SE_new_viral  <- readRDS("Code/Simulation/bootstrap_par_SE_new_viral_boot.RDS")
bootstrap_nonpar_SE_new_viral  <- readRDS("Code/Simulation/bootstrap_nonpar_SE_new_viral_boot.RDS")
estimates_naive <- readRDS("Code/Simulation/estimates_naive_boot.RDS")
SE_naive <- readRDS("Code/Simulation/SE_naive_boot.RDS")
estimates_empirical <- readRDS("Code/Simulation/estimates_empirical_boot.RDS")
SE_empirical <- readRDS("Code/Simulation/SE_empirical_boot.RDS")

