# Install libraries if needed
if (!require("data.table", quietly = TRUE)) {install.packages("data.table")}
if (!require("INLA", quietly = TRUE)) {install.packages("INLA")}
if (!require("spdep", quietly = TRUE)) {install.packages("spdep")}
if (!require("knitr", quietly = TRUE)) {install.packages("knitr")}

# Load libraries
library("data.table")
library("INLA")
library("spdep")
library("knitr")
library("kableExtra")

inla.ShANOVA.2x2 <- function(obs, exp, gr, fac.names = NULL, lev.fac1 = NULL, lev.fac2 = NULL) {
  
  ## DISCLAIMER: Observed and expected values have to be given in an specific order. 
    # Consider n the number of areas, 
      # the first n values (1:n) of the obs (exp) should be the ones belonging to the
        # FIRST level (the first position of the lev.fac1 vector argument) of the
        # FIRST factor (the first position of the fac.names vector argument) and
        # to the FIRST level (the first position of the lev.fac2 vector argument) 
        # of the SECOND factor (the second position of the fac.names vector argument).
      # The n following values ((n+1):2n) of the obs (exp) should be the ones 
        # belonging to the FIRST level (the first position of the lev.fac1 vector 
        # argument) of the FIRST factor (the first position of the fac.names 
        # vector argument) and to the SECOND level (the second position of the 
        # lev.fac2 vector argument) of the SECOND factor (the second position of
        # the fac.names vector argument).
      # The n following values ((2n+1):3n) of the obs (exp) should be the ones 
        # belonging to the SECOND level (the second position of the lev.fac1 vector 
        # argument) of the FIRST factor (the first position of the fac.names 
        # vector argument) and to the FIRST level (the first position of the 
        # lev.fac2 vector argument) of the SECOND factor (the second position of
        # the fac.names vector argument).
      # The n following values ((3n+1):4n) of the obs (exp) should be the ones 
        # belonging to the SECOND level (the second position of the lev.fac1 vector 
        # argument) of the FIRST factor (the first position of the fac.names 
        # vector argument) and to the SECOND level (the second position of the 
        # lev.fac2 vector argument) of the SECOND factor (the second position of
        # the fac.names vector argument).
      # The first n values are O1, the second n values are O2, the third n values
        # are O3 and the last n values are O4
  
  ## Print warnings (warnings are printed as they occur)
  options(warn=1)
  
  ## Define uniform priors for the standard deviation
  prior.prec <- list(prior = "expression: log_dens = 0 - log(2) - theta / 2; return(log_dens);", initial = 0)
  
  ## Define number of levels and number of areas
  n.levels <- c(2, 2)
  n.areas <- gr$n
  
  ## Get maximum number of groups/diseases
  n.groups <- prod(n.levels)
  
  ## Check if obs and exp are the same
  if (sum(obs) != as.integer(sum(exp))){
    stop("ERROR: Number of total observations and expected values are not the same.")
  }
  
  ## Check if data has the proper length
  if (length(obs) != n.areas*n.groups){
    stop("ERROR: The length of the observations does not correspond with the number of groups and areas.")
  }
  
  if (length(exp) != n.areas*n.groups){
    stop("ERROR: The length of the expected values does not correspond with the number of groups and areas.")
  }
  
  ## Empty list
  data.INLA <- list(OBS_f1l1_f2l1 = matrix(NA, nrow = n.areas*n.groups, ncol = n.groups), EXP_f1l1_f2l1 = matrix(NA, nrow = n.areas),
                    OBS_f1l2_f2l1 = matrix(NA, nrow = n.areas*n.groups, ncol = n.groups), EXP_f1l2_f2l1 = matrix(NA, nrow = n.areas),
                    OBS_f1l1_f2l2 = matrix(NA, nrow = n.areas*n.groups, ncol = n.groups), EXP_f1l1_f2l2 = matrix(NA, nrow = n.areas),
                    OBS_f1l2_f2l2 = matrix(NA, nrow = n.areas*n.groups, ncol = n.groups), EXP_f1l2_f2l2 = matrix(NA, nrow = n.areas),
                    OBS_f2l1_f1l1 = matrix(NA, nrow = n.areas*n.groups, ncol = n.groups), EXP_f2l1_f1l1 = matrix(NA, nrow = n.areas),
                    OBS_f2l2_f1l1 = matrix(NA, nrow = n.areas*n.groups, ncol = n.groups), EXP_f2l2_f1l1 = matrix(NA, nrow = n.areas),
                    OBS_f2l1_f1l2 = matrix(NA, nrow = n.areas*n.groups, ncol = n.groups), EXP_f2l1_f1l2 = matrix(NA, nrow = n.areas),
                    OBS_f2l2_f1l2 = matrix(NA, nrow = n.areas*n.groups, ncol = n.groups), EXP_f2l2_f1l2 = matrix(NA, nrow = n.areas))
  
  ## Here we order the Observed and Expected values following the default order: ------------------
    # O1, O2, O3 and O4
      data.INLA$OBS_f1l1_f2l1[1:n.areas, 1] <- obs[1:n.areas] # O1 
      data.INLA$OBS_f1l1_f2l1[n.areas + 1:n.areas, 2] <- obs[n.areas + 1:n.areas] # O2 
      data.INLA$OBS_f1l1_f2l1[2*n.areas + 1:n.areas, 3] <- obs[2*n.areas + 1:n.areas] # O3 
      data.INLA$OBS_f1l1_f2l1[3*n.areas + 1:n.areas, 4] <- obs[3*n.areas + 1:n.areas] # O4 
    # E1, E2, E3 and E4
      data.INLA$EXP_f1l1_f2l1[1:n.areas] <- exp[1:n.areas] # E1
      data.INLA$EXP_f1l1_f2l1[n.areas + 1:n.areas] <- exp[n.areas + 1:n.areas] # E2
      data.INLA$EXP_f1l1_f2l1[2*n.areas + 1:n.areas] <- exp[2*n.areas + 1:n.areas] # E3
      data.INLA$EXP_f1l1_f2l1[3*n.areas + 1:n.areas] <- exp[3*n.areas + 1:n.areas] # E4
      
  ## Here we order the Observed and Expected values following this order: ------------------
    # O3, O4, O1 and O2
      data.INLA$OBS_f1l2_f2l1[1:n.areas, 1] <- obs[2*n.areas + 1:n.areas] # O3
      data.INLA$OBS_f1l2_f2l1[n.areas + 1:n.areas, 2] <- obs[3*n.areas + 1:n.areas] # O4
      data.INLA$OBS_f1l2_f2l1[2*n.areas + 1:n.areas, 3] <- obs[1:n.areas] # O1
      data.INLA$OBS_f1l2_f2l1[3*n.areas + 1:n.areas, 4] <- obs[n.areas + 1:n.areas] # O2
    # E3, E4, E1 and E2
      data.INLA$EXP_f1l2_f2l1[1:n.areas] <- exp[2*n.areas + 1:n.areas] # E3
      data.INLA$EXP_f1l2_f2l1[n.areas + 1:n.areas] <- exp[3*n.areas + 1:n.areas] # E4
      data.INLA$EXP_f1l2_f2l1[2*n.areas + 1:n.areas] <- exp[1:n.areas] # E1
      data.INLA$EXP_f1l2_f2l1[3*n.areas + 1:n.areas] <- exp[n.areas + 1:n.areas] # E2
  
  ## Here we order the Observed and Expected values following this order: ------------------
    # O2, O1, O4 and O3
      data.INLA$OBS_f1l1_f2l2[1:n.areas, 1] <- obs[n.areas + 1:n.areas]  # O2
      data.INLA$OBS_f1l1_f2l2[n.areas + 1:n.areas, 2] <- obs[1:n.areas] # O1
      data.INLA$OBS_f1l1_f2l2[2*n.areas + 1:n.areas, 3] <- obs[3*n.areas + 1:n.areas] # O4
      data.INLA$OBS_f1l1_f2l2[3*n.areas + 1:n.areas, 4] <- obs[2*n.areas + 1:n.areas] # O3
    # E2, E1, E4 and E3
      data.INLA$EXP_f1l1_f2l2[1:n.areas, 1] <- exp[n.areas + 1:n.areas]  # E2
      data.INLA$EXP_f1l1_f2l2[n.areas + 1:n.areas] <- exp[1:n.areas] # E1
      data.INLA$EXP_f1l1_f2l2[2*n.areas + 1:n.areas] <- exp[3*n.areas + 1:n.areas] # E4
      data.INLA$EXP_f1l1_f2l2[3*n.areas + 1:n.areas] <- exp[2*n.areas + 1:n.areas] # E3
  
  ## Here we order the Observed and Expected values following this order: ------------------
    # O4, O3, O2 and O1
      data.INLA$OBS_f1l2_f2l2[1:n.areas, 1] <- obs[3*n.areas + 1:n.areas] # E4
      data.INLA$OBS_f1l2_f2l2[n.areas + 1:n.areas, 2] <- obs[2*n.areas + 1:n.areas] # E3
      data.INLA$OBS_f1l2_f2l2[2*n.areas + 1:n.areas, 3] <- obs[n.areas + 1:n.areas] # E2
      data.INLA$OBS_f1l2_f2l2[3*n.areas + 1:n.areas, 4] <- obs[1:n.areas] # E1
    # E4, E3, E2 and E1
      data.INLA$EXP_f1l2_f2l2[1:n.areas, 1] <- exp[3*n.areas + 1:n.areas] # E4
      data.INLA$EXP_f1l2_f2l2[n.areas + 1:n.areas] <- exp[2*n.areas + 1:n.areas] # E3
      data.INLA$EXP_f1l2_f2l2[2*n.areas + 1:n.areas] <- exp[n.areas + 1:n.areas] # E2
      data.INLA$EXP_f1l2_f2l2[3*n.areas + 1:n.areas] <- exp[1:n.areas] # E1
  
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
  
  ## Here the observed are ordered following this order: ------------------
    # O1, O3, O2 and O4
      data.INLA$OBS_f2l1_f1l1[1:n.areas, 1] <- obs[1:n.areas] # O1
      data.INLA$OBS_f2l1_f1l1[n.areas + 1:n.areas, 2] <- obs[2*n.areas + 1:n.areas] # O3
      data.INLA$OBS_f2l1_f1l1[2*n.areas + 1:n.areas, 3] <- obs[n.areas + 1:n.areas] # O2
      data.INLA$OBS_f2l1_f1l1[3*n.areas + 1:n.areas, 4] <- obs[3*n.areas + 1:n.areas]  # O4
    # E1, E3, E2 and E4
      data.INLA$EXP_f2l1_f1l1[1:n.areas] <- exp[1:n.areas] # O1
      data.INLA$EXP_f2l1_f1l1[n.areas + 1:n.areas] <- exp[2*n.areas + 1:n.areas] # O3
      data.INLA$EXP_f2l1_f1l1[2*n.areas + 1:n.areas] <- exp[n.areas + 1:n.areas] # O2
      data.INLA$EXP_f2l1_f1l1[3*n.areas + 1:n.areas] <- exp[3*n.areas + 1:n.areas]  # O4
  
  ## Here the observed are ordered following this order: ------------------
    # O2, O4, O1 and O3
      data.INLA$OBS_f2l2_f1l1[1:n.areas, 1] <- obs[n.areas + 1:n.areas] # O2
      data.INLA$OBS_f2l2_f1l1[n.areas + 1:n.areas, 2] <- obs[3*n.areas + 1:n.areas] # O4
      data.INLA$OBS_f2l2_f1l1[2*n.areas + 1:n.areas, 3] <- obs[1:n.areas] # O1
      data.INLA$OBS_f2l2_f1l1[3*n.areas + 1:n.areas, 4] <- obs[2*n.areas + 1:n.areas] # O3
    # E2, E4, E1 and E3 
      data.INLA$EXP_f2l2_f1l1[1:n.areas] <- exp[n.areas + 1:n.areas] # E2
      data.INLA$EXP_f2l2_f1l1[n.areas + 1:n.areas] <- exp[3*n.areas + 1:n.areas] # E4
      data.INLA$EXP_f2l2_f1l1[2*n.areas + 1:n.areas] <- exp[1:n.areas] # E1
      data.INLA$EXP_f2l2_f1l1[3*n.areas + 1:n.areas] <- exp[2*n.areas + 1:n.areas] # E3
  
  ## Here the observed are ordered following this order: ------------------
    # O3, O1, O4 and O2
      data.INLA$OBS_f2l1_f1l2[1:n.areas, 1] <- obs[2*n.areas + 1:n.areas] # O3
      data.INLA$OBS_f2l1_f1l2[n.areas + 1:n.areas, 2] <- obs[1:n.areas] # O1
      data.INLA$OBS_f2l1_f1l2[2*n.areas + 1:n.areas, 3] <- obs[3*n.areas + 1:n.areas] # O4
      data.INLA$OBS_f2l1_f1l2[3*n.areas + 1:n.areas, 4] <- obs[n.areas + 1:n.areas] # O2
    # E3, E1, E4 and E2
      data.INLA$EXP_f2l1_f1l2[1:n.areas] <- exp[2*n.areas + 1:n.areas] # E3
      data.INLA$EXP_f2l1_f1l2[n.areas + 1:n.areas] <- exp[1:n.areas] # E1
      data.INLA$EXP_f2l1_f1l2[2*n.areas + 1:n.areas] <- exp[3*n.areas + 1:n.areas] # E4
      data.INLA$EXP_f2l1_f1l2[3*n.areas + 1:n.areas] <- exp[n.areas + 1:n.areas] # E2
  
  ## Here the observed are ordered following this order: ------------------
    # O4, O2, O3 and O1
      data.INLA$OBS_f2l2_f1l2[1:n.areas, 1] <- obs[3*n.areas + 1:n.areas] # O4
      data.INLA$OBS_f2l2_f1l2[n.areas + 1:n.areas, 2] <- obs[n.areas + 1:n.areas] # O2
      data.INLA$OBS_f2l2_f1l2[2*n.areas + 1:n.areas, 3] <- obs[2*n.areas + 1:n.areas] # O3
      data.INLA$OBS_f2l2_f1l2[3*n.areas + 1:n.areas, 4] <- obs[1:n.areas] # O1
    # E4, E2, E3 and E1
      data.INLA$EXP_f2l2_f1l2[1:n.areas] <- exp[3*n.areas + 1:n.areas] # E4
      data.INLA$EXP_f2l2_f1l2[n.areas + 1:n.areas] <- exp[n.areas + 1:n.areas] # E2
      data.INLA$EXP_f2l2_f1l2[2*n.areas + 1:n.areas] <- exp[2*n.areas + 1:n.areas] # E3
      data.INLA$EXP_f2l2_f1l2[3*n.areas + 1:n.areas] <- exp[1:n.areas] # E1
  
  ## Intercepts for each group ------------------
    # Here the sequence is: "f1l1-f2l1" "f1l1-f2l2" "f1l2-f2l1" "f1l2-f2l2"
    data.INLA$alpha_f1l1_f2l1 <- rep(c(paste0(lev.fac1[1], "-", lev.fac2[1]), paste0(lev.fac1[1], "-", lev.fac2[2]), 
                                  paste0(lev.fac1[2], "-", lev.fac2[1]), paste0(lev.fac1[2], "-", lev.fac2[2])), 
                                each = n.areas)
    data.INLA$alpha_f1l1_f2l1 <- as.factor(data.INLA$alpha_f1l1_f2l1)
  
    # Here the sequence is: "f1l2-f2l1" "f1l2-f2l2" "f1l1-f2l1" "f1l1-f2l2"
    data.INLA$alpha_f1l2_f2l1 <- rep(c(paste0(lev.fac1[2], "-", lev.fac2[1]), paste0(lev.fac1[2], "-", lev.fac2[2]), 
                                  paste0(lev.fac1[1], "-", lev.fac2[1]), paste0(lev.fac1[1], "-", lev.fac2[2])), 
                                each = n.areas)
    data.INLA$alpha_f1l2_f2l1 <- as.factor(data.INLA$alpha_f1l2_f2l1)
  
    # Here the sequence is: "f1l1-f2l2" "f1l1-f2l1" "f1l2-f2l2" "f1l2-f2l1"
    data.INLA$alpha_f1l1_f2l2 <- rep(c(paste0(lev.fac1[1], "-", lev.fac2[2]), paste0(lev.fac1[1], "-", lev.fac2[1]), 
                                  paste0(lev.fac1[2], "-", lev.fac2[2]), paste0(lev.fac1[2], "-", lev.fac2[1])), 
                                each = n.areas)
    data.INLA$alpha_f1l1_f2l2 <- as.factor(data.INLA$alpha_f1l1_f2l2)
  
    # Here the sequence is: "f1l2-f2l2" "f1l2-f2l1" "f1l1-f2l2" "f1l1-f2l1"
    data.INLA$alpha_f1l2_f2l2 <- rep(c(paste0(lev.fac1[2], "-", lev.fac2[2]), paste0(lev.fac1[2], "-", lev.fac2[1]), 
                                  paste0(lev.fac1[1], "-", lev.fac2[2]), paste0(lev.fac1[1], "-", lev.fac2[1])), 
                                each = n.areas)
    data.INLA$alpha_f1l2_f2l2 <- as.factor(data.INLA$alpha_f1l2_f2l2)
  
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
  
    # Here the sequence is: "f2l1-f1l1" "f2l1-f1l2" "f2l2-f1l1" "f2l2-f1l2"
    data.INLA$alpha_f2l1_f1l1 <- rep(c(paste0(lev.fac2[1], "-", lev.fac1[1]), paste0(lev.fac2[1], "-", lev.fac1[2]), 
                                       paste0(lev.fac2[2], "-", lev.fac1[1]), paste0(lev.fac2[2], "-", lev.fac1[2])), 
                                     each = n.areas)
    data.INLA$alpha_f2l1_f1l1 <- as.factor(data.INLA$alpha_f2l1_f1l1)
  
    # Here the sequence is: "f2l2-f1l1" "f2l2-f1l2" "f2l1-f1l1" "f2l1-f1l2"
    data.INLA$alpha_f2l2_f1l1 <- rep(c(paste0(lev.fac2[2], "-", lev.fac1[1]), paste0(lev.fac2[2], "-", lev.fac1[2]), 
                                       paste0(lev.fac2[1], "-", lev.fac1[1]), paste0(lev.fac2[1], "-", lev.fac1[2])), 
                                     each = n.areas)
    data.INLA$alpha_f2l2_f1l1 <- as.factor(data.INLA$alpha_f2l2_f1l1)
  
    # Here the sequence is: "f2l1-f1l2" "f2l1-f1l1" "f2l2-f1l2" "f2l2-f1l1"
    data.INLA$alpha_f2l1_f1l2 <- rep(c(paste0(lev.fac2[1], "-", lev.fac1[2]), paste0(lev.fac2[1], "-", lev.fac1[1]), 
                                       paste0(lev.fac2[2], "-", lev.fac1[2]), paste0(lev.fac2[2], "-", lev.fac1[1])), 
                                     each = n.areas)
    data.INLA$alpha_f2l1_f1l2 <- as.factor(data.INLA$alpha_f2l1_f1l2)
  
    # Here the sequence is: "f2l2-f1l2" "f2l2-f1l1" "f2l1-f1l2" "f2l1-f1l1"
    data.INLA$alpha_f2l2_f1l2 <- rep(c(paste0(lev.fac2[2], "-", lev.fac1[2]), paste0(lev.fac2[2], "-", lev.fac1[1]), 
                                       paste0(lev.fac2[1], "-", lev.fac1[2]), paste0(lev.fac2[1], "-", lev.fac1[1])), 
                                     each = n.areas)
    data.INLA$alpha_f2l2_f1l2 <- as.factor(data.INLA$alpha_f2l2_f1l2)
    
    
  ## Create IDs for each area and group ------------------
  data.INLA$AREAID <- rep(1:n.areas, n.groups)
    
  ## Create IDs for shared spatial effects - phi ------------------
  data.INLA$ID_g1 <- data.INLA$AREAID
  data.INLA$ID_g1[-c(1:n.areas)] <- NA 
  data.INLA$ID_g2 <- data.INLA$AREAID
  data.INLA$ID_g2[-(n.areas + 1:n.areas)] <- NA 
  data.INLA$ID_g3 <- data.INLA$AREAID
  data.INLA$ID_g3[-(2 * n.areas + 1:n.areas)] <- NA 
  data.INLA$ID_g4 <- data.INLA$AREAID
  data.INLA$ID_g4[-(3 * n.areas + 1:n.areas)] <- NA 
  
  ## Create IDs for IID Effect - omega_j ------------------ 
  data.INLA$omega_j <- 1:(542*4)
  
  ####################################################### MODEL 0 ####################################################### 
  # This scenario considers a different intercept in the linear predictor of each group
  # and one individual IID Effect. (M0 in the article)
  
  # Formulas for the model 
  formula <- OBS_f1l1_f2l1 ~ 0 + 
    
    # Intercept
    alpha_f1l1_f2l1 +
    
    # IID Effect - omega_j
    f(omega_j, model = "iid") 
  
  # Run model in INLA
  M0 <- inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l1_f2l1, 
               control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE), 
               control.inla = list(tolerance.step = 1e-8))
  
  ####################################################### MODEL 1 ####################################################### 
  # This scenario considers only a different intercept and a different spatial 
  # effect in the linear predictor of each group. (M1 in the article)
  
  # Individual spatial effects - phi_g
  data.INLA$phi_1 <- data.INLA$ID_g1
  data.INLA$phi_2 <- data.INLA$ID_g2
  data.INLA$phi_3 <- data.INLA$ID_g3
  data.INLA$phi_4 <- data.INLA$ID_g4
  
  # Formulas for the model 
  formula <- OBS_f1l1_f2l1 ~ 0 + 
    
    # Intercept
    alpha_f1l1_f2l1 +
    
    # Individual spatial effects - phi_g
    f(phi_1, model = "besag", graph = gr, 
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=TRUE, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) + 
    f(phi_2, model = "besag", graph = gr, 
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=TRUE, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) + 
    f(phi_3, model = "besag", graph = gr, 
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=TRUE, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) + 
    f(phi_4, model = "besag", graph = gr, 
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=TRUE, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    
    # IID Effect - omega_j
    f(omega_j, model = "iid") 
  
  # Run Model
  M1 <- inla(formula = formula,data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l1_f2l1, 
               control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE), 
               control.inla = list(tolerance.step = 1e-8))
  
  ####################################################### MODEL 2 ####################################################### 
  # This scenario considers a different intercept and the same shared spatial 
  # effect in the linear predictor of each group. (M2 in the article)
  # Four different combinations have to be fit where each group is considered as 
  # a reference group in each of the four models. 
  
  # Define IDs for overall shared spatial effect - phi_11
  data.INLA$phi_11 <- data.INLA$ID_g1
  data.INLA$phi_11_g2 <- data.INLA$ID_g2
  data.INLA$phi_11_g3 <- data.INLA$ID_g3
  data.INLA$phi_11_g4 <- data.INLA$ID_g4
  
  # IID Effect - omega_j
  data.INLA$omega_j <- 1:(542*4)
  
  ### F1L1-F2L1 ---------------------------------------------------------------------------------------------------------------------
  
    ## Formulas for the model 
    formula <- OBS_f1l1_f2l1 ~ 0 + 
      
      # Intercept
      alpha_f1l1_f2l1 +
      
      # Overall shared spatial effect - phi_11
      f(phi_11, model = "besag", graph = gr, 
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
        ) + 
      f(phi_11_g2, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g3, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g4, copy = "phi_11", fixed = TRUE) +
      
      # IID Effect - omega_j
      f(omega_j, model = "iid") 
    
    ## Run model in INLA
    M2_f1l1_f2l1 <- inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l1_f2l1, 
                           control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE), 
                           control.inla = list(tolerance.step = 1e-5))
  
  ### F1L2-F2L1 ---------------------------------------------------------------------------------------------------------------------
    
    ## Formulas for the model 
    formula <- OBS_f1l2_f2l1 ~ 0 + 
      
      # Intercept
      alpha_f1l2_f2l1 +
      
      # Overall shared spatial effect - phi_11
      f(phi_11, model = "besag", graph = gr, 
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      f(phi_11_g2, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g3, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g4, copy = "phi_11", fixed = TRUE) +
      
      # IID Effect - omega_j
      f(omega_j, model = "iid") 
    
    ## Run model in INLA
    M2_f1l2_f2l1 <- inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l2_f2l1, 
                             control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE), 
                             control.inla = list(tolerance.step = 1e-5))
    
  ### F1L1-F2L2 ---------------------------------------------------------------------------------------------------------------------
  
    ## Formulas for the model 
    formula <- OBS_f1l1_f2l2 ~ 0 + 
      
      # Intercept
      alpha_f1l1_f2l2 +
      
      # Overall shared spatial effect - phi_11
      f(phi_11, model = "besag", graph = gr, 
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      f(phi_11_g2, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g3, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g4, copy = "phi_11", fixed = TRUE) +
      
      # IID Effect - omega_j
      f(omega_j, model = "iid") 
    
    ## Run model in INLA
    M2_f1l1_f2l2 <- inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l1_f2l2, 
                           control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE), 
                           control.inla = list(tolerance.step = 1e-5))
    
  ### F1L2-F2L2 ---------------------------------------------------------------------------------------------------------------------
  
    ## Formulas for the model 
    formula <- OBS_f1l2_f2l2 ~ 0 + 
      
      # Intercept
      alpha_f1l2_f2l2 +
      
      # Overall shared spatial effect - phi_11
      f(phi_11, model = "besag", graph = gr, 
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      f(phi_11_g2, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g3, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g4, copy = "phi_11", fixed = TRUE) +
      
      # IID Effect - omega_j
      f(omega_j, model = "iid") 
    
    ## Run model in INLA
    M2_f1l2_f2l2 <- inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l2_f2l2, 
                           control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE), 
                           control.inla = list(tolerance.step = 1e-5))
  
  ####################################################### MODEL 3 ####################################################### 
  # This scenario considers a different intercept, the same shared spatial 
  # effect in the linear predictor of each group and a specific efect for the 
  # second category of the first factor (M3 in the article). Two different
  # combinations have to be fit where each of the two levels is considered
  # as a reference category of the first factor. 
  
  # Overall shared spatial effect - phi_11
  data.INLA$phi_11 <- data.INLA$ID_g1
  data.INLA$phi_11_g2 <- data.INLA$ID_g2
  data.INLA$phi_11_g3 <- data.INLA$ID_g3
  data.INLA$phi_11_g4 <- data.INLA$ID_g4
  
  # Define IDs for 2-group shared effect - phi_12
  data.INLA$phi_12  <- data.INLA$ID_g3
  data.INLA$phi_12_g4 <- data.INLA$ID_g4
  
  ### F1L1-F2L1 ---------------------------------------------------------------------------------------------------------------------
  
    ## Formulas for the model 
    formula <- OBS_f1l1_f2l1 ~ 0 + 
      
      # Intercept
      alpha_f1l1_f2l1 +
      
      # Overall shared spatial effect - phi_11
      f(phi_11, model = "besag", graph = gr, 
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      f(phi_11_g2, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g3, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g4, copy = "phi_11", fixed = TRUE) +
      
      # 2 groups shared effect - phi_12
      f(phi_12, model = "besag", graph = gr, 
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      f(phi_12_g4, copy = "phi_12", fixed = TRUE) +
      
      # IID Effect - omega_j
      f(omega_j, model = "iid") 
    
    ## Run model in INLA
    M3_f1l1_f2l1 <- inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l1_f2l1, 
                 control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE), 
                 control.inla = list(tolerance.step = 1e-5))
  
  ### F1L2-F2L1 ---------------------------------------------------------------------------------------------------------------------
  
    ## Formulas for the model 
    formula <- OBS_f1l2_f2l1 ~ 0 + 
      
      # Intercept
      alpha_f1l2_f2l1 +
      
      # Overall shared spatial effect - phi_11
      f(phi_11, model = "besag", graph = gr, 
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      f(phi_11_g2, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g3, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g4, copy = "phi_11", fixed = TRUE) +
      
      # 2 groups shared effect - phi_12
      f(phi_12, model = "besag", graph = gr, 
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      f(phi_12_g4, copy = "phi_12", fixed = TRUE) +
      
      # IID Effect - omega_j
      f(omega_j, model = "iid") 
    
    ## Run model in INLA
    M3_f1l2_f2l1 <- inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l2_f2l1, 
                             control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE), 
                             control.inla = list(tolerance.step = 1e-5))
  
  ####################################################### MODEL 4 ####################################################### 
  # This scenario considers a different intercept, the same shared spatial 
  # effect in the linear predictor of each group and a specific efect for the 
  # second category of the second factor (M4 in the article). Two different
  # combinations have to be fit where each of the two levels is considered
  # as a reference category of the second factor. 
  
  # Define IDs for 2-group shared effect - phi_21
  data.INLA$phi_21  <- data.INLA$ID_g2
  data.INLA$phi_21_g4 <- data.INLA$ID_g4
  
  ### F1L1-F2L1 ---------------------------------------------------------------------------------------------------------------------
  
    ## Formulas for the model 
    formula <- OBS_f1l1_f2l1 ~ 0 + 
      
      # Intercept
      alpha_f1l1_f2l1 +
      
      # Overall shared spatial effect - phi_11
      f(phi_11, model = "besag", graph = gr, 
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      f(phi_11_g2, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g3, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g4, copy = "phi_11", fixed = TRUE) +
      
      # 2 groups shared effect - phi_12
      f(phi_21, model = "besag", graph = gr, 
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      f(phi_21_g4, copy = "phi_21", fixed = TRUE) +
      
      # IID Effect - omega_j
      f(omega_j, model = "iid") 
    
    ## Run model in INLA
    M4_f1l1_f2l1 <- inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l1_f2l1, 
                             control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE), 
                             control.inla = list(tolerance.step = 1e-5))
  
  ### F1L1-F2L2 ---------------------------------------------------------------------------------------------------------------------
  
    ## Formulas for the model 
    formula <- OBS_f1l1_f2l2 ~ 0 + 
      
      # Intercept
      alpha_f1l1_f2l2 +
      
      # Overall shared spatial effect - phi_11
      f(phi_11, model = "besag", graph = gr, 
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      f(phi_11_g2, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g3, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g4, copy = "phi_11", fixed = TRUE) +
      
      # 2 groups shared effect - phi_12
      f(phi_21, model = "besag", graph = gr, 
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      f(phi_21_g4, copy = "phi_21", fixed = TRUE) +
      
      # IID Effect - omega_j
      f(omega_j, model = "iid") 
    
    ## Run model in INLA
    M4_f1l1_f2l2 <- inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l1_f2l2, 
                             control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE), 
                             control.inla = list(tolerance.step = 1e-5))
  
  ####################################################### MODEL 5 ####################################################### 
  # This scenario considers a different intercept, the same shared spatial 
  # effect in the linear predictor of each group, a specific spatial effect for
  # the second category of the first factor and a diferent specific spatial 
  # effect the second category of the second factor. (M5 in the article)
  # Four different combinations have to be fit.
  
  # Define IDs for overall shared spatial effect - phi_11
  data.INLA$phi_11 <- data.INLA$ID_g1
  data.INLA$phi_11_g2 <- data.INLA$ID_g2
  data.INLA$phi_11_g3 <- data.INLA$ID_g3
  data.INLA$phi_11_g4 <- data.INLA$ID_g4
  
  # Define IDs for 2 groups shared effect - phi_12
  data.INLA$phi_12  <- data.INLA$ID_g3
  data.INLA$phi_12_g4 <- data.INLA$ID_g4
  
  # Define IDs for 2 groups shared effect - phi_21
  data.INLA$phi_21 <- data.INLA$ID_g2
  data.INLA$phi_21_g4 <- data.INLA$ID_g4
  
  ### F1L1-F2L1 ---------------------------------------------------------------------------------------------------------------------
  
  ## Formulas for the model 
  formula <- OBS_f1l1_f2l1 ~ 0 + 
    # Intercept
    alpha_f1l1_f2l1 +
    
    # Overall shared spatial effect - phi_11
    f(phi_11, model = "besag", graph = gr, 
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=TRUE, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) + 
    f(phi_11_g2, copy = "phi_11", fixed = TRUE) + 
    f(phi_11_g3, copy = "phi_11", fixed = TRUE) + 
    f(phi_11_g4, copy = "phi_11", fixed = TRUE) +
    
    # 2 groups shared effect - phi_12
    f(phi_12, model = "besag", graph = gr, 
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=TRUE, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) + 
    f(phi_12_g4, copy = "phi_12", fixed = TRUE) + 
    
    # 2 groups shared effect - phi_21
    f(phi_21, model = "besag", graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=TRUE, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) + 
    f(phi_21_g4, copy = "phi_21", fixed = TRUE) +
    
    # IID Effect - omega_j
    f(omega_j, model = "iid") 
  
  ## Run model in INLA
  M5_f1l1_f2l1 <- inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l1_f2l1, 
                 control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE), 
                 control.inla = list(tolerance.step = 1e-5))
  
  ### F1L2-F2L1 ---------------------------------------------------------------------------------------------------------------------
  
  ## Formulas for the model 
  formula <- OBS_f1l2_f2l1 ~ 0 + 
    # Intercept
    alpha_f1l2_f2l1 +
    
    # Overall shared spatial effect - phi_11
    f(phi_11, model = "besag", graph = gr, 
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=TRUE, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) + 
    f(phi_11_g2, copy = "phi_11", fixed = TRUE) + 
    f(phi_11_g3, copy = "phi_11", fixed = TRUE) + 
    f(phi_11_g4, copy = "phi_11", fixed = TRUE) +
    
    # 2 groups shared effect - phi_12
    f(phi_12, model = "besag", graph = gr, 
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=TRUE, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) + 
    f(phi_12_g4, copy = "phi_12", fixed = TRUE) + 
    
    # 2 groups shared effect - phi_21
    f(phi_21, model = "besag", graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=TRUE, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) + 
    f(phi_21_g4, copy = "phi_21", fixed = TRUE) +
    
    # IID Effect - omega_j
    f(omega_j, model = "iid") 
  
  ## Run model in INLA
  M5_f1l2_f2l1 <- inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l2_f2l1, 
                           control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE), 
                           control.inla = list(tolerance.step = 1e-5))
  
  ### F1L1-F2L2 ---------------------------------------------------------------------------------------------------------------------
  
  ## Formulas for the model 
  formula <- OBS_f1l1_f2l2 ~ 0 + 
    # Intercept
    alpha_f1l1_f2l2 +
    
    # Overall shared spatial effect - phi_11
    f(phi_11, model = "besag", graph = gr, 
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=TRUE, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) + 
    f(phi_11_g2, copy = "phi_11", fixed = TRUE) + 
    f(phi_11_g3, copy = "phi_11", fixed = TRUE) + 
    f(phi_11_g4, copy = "phi_11", fixed = TRUE) +
    
    # 2 groups shared effect - phi_12
    f(phi_12, model = "besag", graph = gr, 
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=TRUE, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) + 
    f(phi_12_g4, copy = "phi_12", fixed = TRUE) + 
    
    # 2 groups shared effect - phi_21
    f(phi_21, model = "besag", graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=TRUE, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) + 
    f(phi_21_g4, copy = "phi_21", fixed = TRUE) +
    
    # IID Effect - omega_j
    f(omega_j, model = "iid") 
  
  ## Run model in INLA
  M5_f1l1_f2l2 <- inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l1_f2l2, 
                           control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE), 
                           control.inla = list(tolerance.step = 1e-5))
  
  ### F1L2-F2L2 ---------------------------------------------------------------------------------------------------------------------
  
  ## Formulas for the model 
  formula <- OBS_f1l2_f2l2 ~ 0 + 
    # Intercept
    alpha_f1l2_f2l2 +
    
    # Overall shared spatial effect - phi_11
    f(phi_11, model = "besag", graph = gr, 
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=TRUE, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) + 
    f(phi_11_g2, copy = "phi_11", fixed = TRUE) + 
    f(phi_11_g3, copy = "phi_11", fixed = TRUE) + 
    f(phi_11_g4, copy = "phi_11", fixed = TRUE) +
    
    # 2 groups shared effect - phi_12
    f(phi_12, model = "besag", graph = gr, 
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=TRUE, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) + 
    f(phi_12_g4, copy = "phi_12", fixed = TRUE) + 
    
    # 2 groups shared effect - phi_21
    f(phi_21, model = "besag", graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=TRUE, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) + 
    f(phi_21_g4, copy = "phi_21", fixed = TRUE) +
    
    # IID Effect - omega_j
    f(omega_j, model = "iid") 
  
  ## Run model in INLA
  M5_f1l2_f2l2 <- inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l2_f2l2, 
                           control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE), 
                           control.inla = list(tolerance.step = 1e-5))
  
  ####################################################### MODEL 6 ####################################################### 
  # This scenario considers a different intercept, the same shared spatial 
  # effect in the linear predictor of each group, a specific spatial effect for
  # the second category of the first factor,  a diferent specific spatial 
  # effect the second category of the second factor (just in one group), and
  # another specific spatial effect for the interaction between two factors. 
  # (M6 in the article)
  # Eight different combinations have to be fit.
  
  # Overall shared spatial effect - phi_11
  data.INLA$phi_11 <- data.INLA$ID_g1
  data.INLA$phi_11_g2 <- data.INLA$ID_g2
  data.INLA$phi_11_g3 <- data.INLA$ID_g3
  data.INLA$phi_11_g4 <- data.INLA$ID_g4
  
  # # Define IDs for 2 groups shared effect - phi_12
  data.INLA$phi_12  <- data.INLA$ID_g3
  data.INLA$phi_12_g4 <- data.INLA$ID_g4
  
  # # Define IDs for individual effect 1 - phi_21
  data.INLA$phi_21 <- data.INLA$ID_g2
  
  # # Define IDs for individual effect 2 - phi_22
  data.INLA$phi_22 <- data.INLA$ID_g4
  
  ### F1L1-F2L1 ---------------------------------------------------------------------------------------------------------------------
  
    ## Formulas for the model 
    formula <- OBS_f1l1_f2l1 ~ 0 + 
      # Intercept
      alpha_f1l1_f2l1 +
      
      # Overall shared spatial effect - phi_11
      f(phi_11, model = "besag", graph = gr, 
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      f(phi_11_g2, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g3, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g4, copy = "phi_11", fixed = TRUE) +
      
      # 2 groups shared effect - phi_12
      f(phi_12, model = "besag", graph = gr, 
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      f(phi_12_g4, copy = "phi_12", fixed = TRUE) + 
      
      # Individual effect 1 - phi_21
      f(phi_21, model = "besag", graph = gr,
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      
      # Individual effect 1 - phi_22
      f(phi_22, model = "besag", graph = gr,
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) +
      
      # IID Effect - omega_j
      f(omega_j, model = "iid") 
    
    ## Run model in INLA
    M6_f1l1_f2l1 <- inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l1_f2l1, 
                   control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE), 
                   control.inla = list(tolerance.step = 1e-5))
  
  ### F1L2-F2L1 ---------------------------------------------------------------------------------------------------------------------
  
    ## Formulas for the model 
    formula <- OBS_f1l2_f2l1 ~ 0 + 
      # Intercept
      alpha_f1l2_f2l1 +
      
      # Overall shared spatial effect - phi_11
      f(phi_11, model = "besag", graph = gr, 
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      f(phi_11_g2, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g3, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g4, copy = "phi_11", fixed = TRUE) +
      
      # 2 groups shared effect - phi_12
      f(phi_12, model = "besag", graph = gr, 
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      f(phi_12_g4, copy = "phi_12", fixed = TRUE) + 
      
      # Individual effect 1 - phi_21
      f(phi_21, model = "besag", graph = gr,
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      
      # Individual effect 1 - phi_22
      f(phi_22, model = "besag", graph = gr,
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) +
      
      # IID Effect - omega_j
      f(omega_j, model = "iid") 
    
    ## Run model in INLA
    M6_f1l2_f2l1 <- inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l2_f2l1, 
                             control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE), 
                             control.inla = list(tolerance.step = 1e-5))
  
  ### F1L1-F2L2 ---------------------------------------------------------------------------------------------------------------------
  
    ## Formulas for the model 
    formula <- OBS_f1l1_f2l2 ~ 0 + 
      # Intercept
      alpha_f1l1_f2l2 +
      
      # Overall shared spatial effect - phi_11
      f(phi_11, model = "besag", graph = gr, 
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      f(phi_11_g2, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g3, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g4, copy = "phi_11", fixed = TRUE) +
      
      # 2 groups shared effect - phi_12
      f(phi_12, model = "besag", graph = gr, 
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      f(phi_12_g4, copy = "phi_12", fixed = TRUE) + 
      
      # Individual effect 1 - phi_21
      f(phi_21, model = "besag", graph = gr,
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      
      # Individual effect 1 - phi_22
      f(phi_22, model = "besag", graph = gr,
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) +
      
      # IID Effect - omega_j
      f(omega_j, model = "iid") 
    
    ## Run model in INLA
    M6_f1l1_f2l2 <- inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l1_f2l2, 
                             control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE), 
                             control.inla = list(tolerance.step = 1e-5))
  
  ### F1L2-F2L2 ---------------------------------------------------------------------------------------------------------------------
  
    ## Formulas for the model 
    formula <- OBS_f1l2_f2l2 ~ 0 + 
      # Intercept
      alpha_f1l2_f2l2 +
      
      # Overall shared spatial effect - phi_11
      f(phi_11, model = "besag", graph = gr, 
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      f(phi_11_g2, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g3, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g4, copy = "phi_11", fixed = TRUE) +
      
      # 2 groups shared effect - phi_12
      f(phi_12, model = "besag", graph = gr, 
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      f(phi_12_g4, copy = "phi_12", fixed = TRUE) + 
      
      # Individual effect 1 - phi_21
      f(phi_21, model = "besag", graph = gr,
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      
      # Individual effect 1 - phi_22
      f(phi_22, model = "besag", graph = gr,
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) +
      
      # IID Effect - omega_j
      f(omega_j, model = "iid") 
    
    ## Run model in INLA
    M6_f1l2_f2l2 <- inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l2_f2l2, 
                             control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE), 
                             control.inla = list(tolerance.step = 1e-5))
  
  ### F2L1-F1L1 ---------------------------------------------------------------------------------------------------------------------
  
    ## Formulas for the model 
    formula <- OBS_f2l1_f1l1 ~ 0 + 
      # Intercept
      alpha_f2l1_f1l1 +
      
      # Overall shared spatial effect - phi_11
      f(phi_11, model = "besag", graph = gr, 
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      f(phi_11_g2, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g3, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g4, copy = "phi_11", fixed = TRUE) +
      
      # 2 groups shared effect - phi_12
      f(phi_12, model = "besag", graph = gr, 
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      f(phi_12_g4, copy = "phi_12", fixed = TRUE) + 
      
      # Individual effect 1 - phi_21
      f(phi_21, model = "besag", graph = gr,
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      
      # Individual effect 1 - phi_22
      f(phi_22, model = "besag", graph = gr,
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) +
      
      # IID Effect - omega_j
      f(omega_j, model = "iid") 
    
    ## Run model in INLA
    M6_f2l1_f1l1 <- inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f2l1_f1l1, 
                             control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE), 
                             control.inla = list(tolerance.step = 1e-5))
  
  ### F2L2-F1L1 ---------------------------------------------------------------------------------------------------------------------
  
    ## Formulas for the model 
    formula <- OBS_f2l2_f1l1 ~ 0 + 
      # Intercept
      alpha_f2l2_f1l1 +
      
      # Overall shared spatial effect - phi_11
      f(phi_11, model = "besag", graph = gr, 
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      f(phi_11_g2, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g3, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g4, copy = "phi_11", fixed = TRUE) +
      
      # 2 groups shared effect - phi_12
      f(phi_12, model = "besag", graph = gr, 
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      f(phi_12_g4, copy = "phi_12", fixed = TRUE) + 
      
      # Individual effect 1 - phi_21
      f(phi_21, model = "besag", graph = gr,
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      
      # Individual effect 1 - phi_22
      f(phi_22, model = "besag", graph = gr,
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) +
      
      # IID Effect - omega_j
      f(omega_j, model = "iid") 
    
    ## Run model in INLA
    M6_f2l2_f1l1 <- inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f2l2_f1l1, 
                             control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE), 
                             control.inla = list(tolerance.step = 1e-5))
  
  ### F2L1-F1L2 ---------------------------------------------------------------------------------------------------------------------
  
    ## Formulas for the model 
    formula <- OBS_f2l1_f1l2 ~ 0 + 
      # Intercept
      alpha_f2l1_f1l2 +
      
      # Overall shared spatial effect - phi_11
      f(phi_11, model = "besag", graph = gr, 
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      f(phi_11_g2, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g3, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g4, copy = "phi_11", fixed = TRUE) +
      
      # 2 groups shared effect - phi_12
      f(phi_12, model = "besag", graph = gr, 
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      f(phi_12_g4, copy = "phi_12", fixed = TRUE) + 
      
      # Individual effect 1 - phi_21
      f(phi_21, model = "besag", graph = gr,
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      
      # Individual effect 1 - phi_22
      f(phi_22, model = "besag", graph = gr,
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) +
      
      # IID Effect - omega_j
      f(omega_j, model = "iid") 
    
    ## Run model in INLA
    M6_f2l1_f1l2 <- inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f2l1_f1l2, 
                             control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE), 
                             control.inla = list(tolerance.step = 1e-5))
  
  ### F2L2-F2L2 ---------------------------------------------------------------------------------------------------------------------
  
    ## Formulas for the model 
    formula <- OBS_f2l2_f1l2 ~ 0 + 
      # Intercept
      alpha_f2l2_f1l2 +
      
      # Overall shared spatial effect - phi_11
      f(phi_11, model = "besag", graph = gr, 
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      f(phi_11_g2, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g3, copy = "phi_11", fixed = TRUE) + 
      f(phi_11_g4, copy = "phi_11", fixed = TRUE) +
      
      # 2 groups shared effect - phi_12
      f(phi_12, model = "besag", graph = gr, 
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      f(phi_12_g4, copy = "phi_12", fixed = TRUE) + 
      
      # Individual effect 1 - phi_21
      f(phi_21, model = "besag", graph = gr,
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) + 
      
      # Individual effect 1 - phi_22
      f(phi_22, model = "besag", graph = gr,
        constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
        scale.model=TRUE, # Scale spatial effects
        hyper = list(prec = prior.prec) # Prior Distributions
      ) +
      
      # IID Effect - omega_j
      f(omega_j, model = "iid") 
    
    ## Run model in INLA
    M6_f2l2_f1l2 <- inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f2l2_f1l2, 
                             control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE), 
                             control.inla = list(tolerance.step = 1e-5))
  
  ####################################################### RESULTS ####################################################### 
    
  ## Add the combination of the groups in the final order to each model for easier interpretation -------------------------------------
    # M0
    M0$groups <- c(paste0(lev.fac1[1], "-", lev.fac2[1], " group"), paste0(lev.fac1[1], "-", lev.fac2[2], " group"), 
                     paste0(lev.fac1[2], "-", lev.fac2[1], " group"), paste0(lev.fac1[2], "-", lev.fac2[2], " group"))
    # M1
    M1$groups <- c(paste0(lev.fac1[1], "-", lev.fac2[1], " group"), paste0(lev.fac1[1], "-", lev.fac2[2], " group"), 
                     paste0(lev.fac1[2], "-", lev.fac2[1], " group"), paste0(lev.fac1[2], "-", lev.fac2[2], " group"))
    # M2
    M2_f1l1_f2l1$groups <- c(paste0(lev.fac1[1], "-", lev.fac2[1], " group"), paste0(lev.fac1[1], "-", lev.fac2[2], " group"), 
                               paste0(lev.fac1[2], "-", lev.fac2[1], " group"), paste0(lev.fac1[2], "-", lev.fac2[2], " group"))
    
    M2_f1l1_f2l2$groups <- c(paste0(lev.fac1[1], "-", lev.fac2[2], " group"), paste0(lev.fac1[1], "-", lev.fac2[1], " group"), 
                               paste0(lev.fac1[2], "-", lev.fac2[2], " group"), paste0(lev.fac1[2], "-", lev.fac2[1], " group"))
    
    M2_f1l2_f2l1$groups <- c(paste0(lev.fac1[2], "-", lev.fac2[1], " group"), paste0(lev.fac1[2], "-", lev.fac2[2], " group"), 
                               paste0(lev.fac1[1], "-", lev.fac2[1], " group"), paste0(lev.fac1[1], "-", lev.fac2[2], " group"))
    
    M2_f1l2_f2l2$groups <- c(paste0(lev.fac1[2], "-", lev.fac2[2], " group"), paste0(lev.fac1[2], "-", lev.fac2[1], " group"), 
                               paste0(lev.fac1[1], "-", lev.fac2[2], " group"), paste0(lev.fac1[1], "-", lev.fac2[1], " group"))
    # M3
    M3_f1l1_f2l1$groups <- c(paste0(lev.fac1[1], "-", lev.fac2[1], " group"), paste0(lev.fac1[1], "-", lev.fac2[2], " group"), 
                               paste0(lev.fac1[2], "-", lev.fac2[1], " group"), paste0(lev.fac1[2], "-", lev.fac2[2], " group"))
    
    M3_f1l2_f2l1$groups <- c(paste0(lev.fac1[2], "-", lev.fac2[1], " group"), paste0(lev.fac1[2], "-", lev.fac2[2], " group"), 
                               paste0(lev.fac1[1], "-", lev.fac2[1], " group"), paste0(lev.fac1[1], "-", lev.fac2[2], " group"))
    # M4
    M4_f1l1_f2l1$groups <- c(paste0(lev.fac1[1], "-", lev.fac2[1], " group"), paste0(lev.fac1[1], "-", lev.fac2[2], " group"), 
                               paste0(lev.fac1[2], "-", lev.fac2[1], " group"), paste0(lev.fac1[2], "-", lev.fac2[2], " group"))
    
    M4_f1l1_f2l2$groups <- c(paste0(lev.fac1[1], "-", lev.fac2[2], " group"), paste0(lev.fac1[1], "-", lev.fac2[1], " group"), 
                               paste0(lev.fac1[2], "-", lev.fac2[2], " group"), paste0(lev.fac1[2], "-", lev.fac2[1], " group"))
    # M5
    M5_f1l1_f2l1$groups <- c(paste0(lev.fac1[1], "-", lev.fac2[1], " group"), paste0(lev.fac1[1], "-", lev.fac2[2], " group"), 
                               paste0(lev.fac1[2], "-", lev.fac2[1], " group"), paste0(lev.fac1[2], "-", lev.fac2[2], " group"))
    
    M5_f1l1_f2l2$groups <- c(paste0(lev.fac1[1], "-", lev.fac2[2], " group"), paste0(lev.fac1[1], "-", lev.fac2[1], " group"), 
                               paste0(lev.fac1[2], "-", lev.fac2[2], " group"), paste0(lev.fac1[2], "-", lev.fac2[1], " group"))
    
    M5_f1l2_f2l1$groups <- c(paste0(lev.fac1[2], "-", lev.fac2[1], " group"), paste0(lev.fac1[2], "-", lev.fac2[2], " group"), 
                               paste0(lev.fac1[1], "-", lev.fac2[1], " group"), paste0(lev.fac1[1], "-", lev.fac2[2], " group"))
    
    M5_f1l2_f2l2$groups <- c(paste0(lev.fac1[2], "-", lev.fac2[2], " group"), paste0(lev.fac1[2], "-", lev.fac2[1], " group"), 
                               paste0(lev.fac1[1], "-", lev.fac2[2], " group"), paste0(lev.fac1[1], "-", lev.fac2[1], " group"))
    # M6
    M6_f1l1_f2l1$groups <- c(paste0(lev.fac1[1], "-", lev.fac2[1], " group"), paste0(lev.fac1[1], "-", lev.fac2[2], " group"), 
                               paste0(lev.fac1[2], "-", lev.fac2[1], " group"), paste0(lev.fac1[2], "-", lev.fac2[2], " group"))
    
    M6_f1l1_f2l2$groups <- c(paste0(lev.fac1[1], "-", lev.fac2[2], " group"), paste0(lev.fac1[1], "-", lev.fac2[1], " group"), 
                               paste0(lev.fac1[2], "-", lev.fac2[2], " group"), paste0(lev.fac1[2], "-", lev.fac2[1], " group"))
    
    M6_f1l2_f2l1$groups <- c(paste0(lev.fac1[2], "-", lev.fac2[1], " group"), paste0(lev.fac1[2], "-", lev.fac2[2], " group"), 
                               paste0(lev.fac1[1], "-", lev.fac2[1], " group"), paste0(lev.fac1[1], "-", lev.fac2[2], " group"))
    
    M6_f1l2_f2l2$groups <- c(paste0(lev.fac1[2], "-", lev.fac2[2], " group"), paste0(lev.fac1[2], "-", lev.fac2[1], " group"), 
                               paste0(lev.fac1[1], "-", lev.fac2[2], " group"), paste0(lev.fac1[1], "-", lev.fac2[1], " group"))
    
    M6_f2l1_f1l1$groups <- c(paste0(lev.fac2[1], "-", lev.fac1[1], " group"), paste0(lev.fac2[1], "-", lev.fac1[2], " group"), 
                               paste0(lev.fac2[2], "-", lev.fac1[1], " group"), paste0(lev.fac2[2], "-", lev.fac1[2], " group"))
    
    M6_f2l1_f1l2$groups <- c(paste0(lev.fac2[1], "-", lev.fac1[2], " group"), paste0(lev.fac2[1], "-", lev.fac1[1], " group"), 
                               paste0(lev.fac2[2], "-", lev.fac1[2], " group"), paste0(lev.fac2[2], "-", lev.fac1[1], " group"))
    
    M6_f2l2_f1l1$groups <- c(paste0(lev.fac2[2], "-", lev.fac1[1], " group"), paste0(lev.fac2[2], "-", lev.fac1[2], " group"), 
                               paste0(lev.fac2[1], "-", lev.fac1[1], " group"), paste0(lev.fac2[1], "-", lev.fac1[2], " group"))
    
    M6_f2l2_f1l2$groups <- c(paste0(lev.fac2[2], "-", lev.fac1[2], " group"), paste0(lev.fac2[2], "-", lev.fac1[1], " group"), 
                               paste0(lev.fac2[1], "-", lev.fac1[2], " group"), paste0(lev.fac2[1], "-", lev.fac1[1], " group"))
  
  ## Add name of the level being adjusted in each shared spatial effect ---------------------------------------------------------------
    
    # M1
    M1$sp_effects <- c(paste0("phi_1(", lev.fac1[1], "-", lev.fac2[1], ")"), paste0("phi_2(", lev.fac1[1], "-", lev.fac2[2], ")"), 
                       paste0("phi_3(", lev.fac1[2], "-", lev.fac2[1], ")"), paste0("phi_4(", lev.fac1[2], "-", lev.fac2[2], ")"))
  
    # M2
    M2_f1l1_f2l1$sp_effects <- c(paste0("phi_11(General[", lev.fac1[1], "-", lev.fac2[1], "]"))
    M2_f1l1_f2l2$sp_effects <- c(paste0("phi_11(General[", lev.fac1[1], "-", lev.fac2[2], "]"))
    M2_f1l2_f2l1$sp_effects <- c(paste0("phi_11(General[", lev.fac1[2], "-", lev.fac2[1], "]"))
    M2_f1l2_f2l2$sp_effects <- c(paste0("phi_11(General[", lev.fac1[2], "-", lev.fac2[2], "]"))
    
    # M3
    M3_f1l1_f2l1$sp_effects <- c(paste0("phi_11(General[", lev.fac1[1], "-", lev.fac2[1], "])"), paste0("phi_12(", lev.fac1[2], " effect", ")"))
    
    M3_f1l2_f2l1$sp_effects <- c(paste0("phi_11(General[", lev.fac1[2], "-", lev.fac2[1], "])"), paste0("phi_12(", lev.fac1[1], " effect", ")"))
    
    # M4
    M4_f1l1_f2l1$sp_effects <- c(paste0("phi_11(General[", lev.fac1[1], "-", lev.fac2[1], "])"), paste0("phi_21(", lev.fac2[2], " effect", ")"))
    
    M4_f1l1_f2l2$sp_effects <- c(paste0("phi_11(General[", lev.fac1[1], "-", lev.fac2[2], "])"), paste0("phi_21(", lev.fac2[1], " effect", ")"))
    
    # M5
    M5_f1l1_f2l1$sp_effects <- c(paste0("phi_11(General[", lev.fac1[1], "-", lev.fac2[1], "])"), paste0("phi_12(", lev.fac1[2], " effect", ")"), 
                                 paste0("phi_21(", lev.fac2[2], " effect", ")"))
    
    M5_f1l1_f2l2$sp_effects <- c(paste0("phi_11(General[", lev.fac1[1], "-", lev.fac2[2], "])"), paste0("phi_12(", lev.fac1[2], " effect", ")"), 
                                 paste0("phi_21(", lev.fac2[1], " effect", ")"))
    
    M5_f1l2_f2l1$sp_effects <- c(paste0("phi_11(General[", lev.fac1[2], "-", lev.fac2[1], "])"), paste0("phi_12(", lev.fac1[1], " effect", ")"), 
                                 paste0("phi_21(", lev.fac2[2], " effect", ")"))
    
    M5_f1l2_f2l2$sp_effects <- c(paste0("phi_11(General[", lev.fac1[2], "-", lev.fac2[2], "])"), paste0("phi_12(", lev.fac1[1], " effect", ")"), 
                                 paste0("phi_21(", lev.fac2[1], " effect", ")"))
    
    # M6
    M6_f1l1_f2l1$sp_effects <- c(paste0("phi_11(General[", lev.fac1[1], "-", lev.fac2[1], "])"), paste0("phi_12(", lev.fac1[2], " effect", ")"), 
                                 paste0("phi_21(", lev.fac1[1], "*", lev.fac2[2], " effect", ")"), paste0("phi_22(", lev.fac1[2], "*", lev.fac2[2], " effect", ")"))
    
    M6_f1l1_f2l2$sp_effects <- c(paste0("phi_11(General[", lev.fac1[1], "-", lev.fac2[2], "])"), paste0("phi_12(", lev.fac1[2], " effect", ")"), 
                                 paste0("phi_21(", lev.fac1[1], "*", lev.fac2[1], " effect", ")"), paste0("phi_22(", lev.fac1[2], "*", lev.fac2[1], " effect", ")"))
    
    M6_f1l2_f2l1$sp_effects <- c(paste0("phi_11(General[", lev.fac1[2], "-", lev.fac2[1], "])"), paste0("phi_12(", lev.fac1[1], " effect", ")"), 
                                 paste0("phi_21(", lev.fac1[2], "*", lev.fac2[2], " effect", ")"), paste0("phi_22(", lev.fac1[1], "*", lev.fac2[2], " effect", ")"))
    
    M6_f1l2_f2l2$sp_effects <- c(paste0("phi_11(General[", lev.fac1[2], "-", lev.fac2[2], "])"), paste0("phi_12(", lev.fac1[1], " effect", ")"), 
                                 paste0("phi_21(", lev.fac1[2], "*", lev.fac2[1], " effect", ")"), paste0("phi_22(", lev.fac1[1], "*", lev.fac2[1], " effect", ")"))
    
    
    M6_f2l1_f1l1$sp_effects <- c(paste0("phi_11(General[", lev.fac2[1], "-", lev.fac1[1], "])"), paste0("phi_12(", lev.fac2[2], " effect", ")"), 
                                 paste0("phi_21(", lev.fac2[1], "*", lev.fac1[2], " effect", ")"), paste0("phi_22(", lev.fac2[2], "*", lev.fac1[2], " effect", ")"))
    
    M6_f2l1_f1l2$sp_effects <- c(paste0("phi_11(General[", lev.fac2[1], "-", lev.fac1[2], "])"), paste0("phi_12(", lev.fac2[2], " effect", ")"), 
                                 paste0("phi_21(", lev.fac2[1], "*", lev.fac1[1], " effect", ")"), paste0("phi_22(", lev.fac2[2], "*", lev.fac1[1], " effect", ")"))
    
    M6_f2l2_f1l1$sp_effects <- c(paste0("phi_11(General[", lev.fac2[2], "-", lev.fac1[1], "])"), paste0("phi_12(", lev.fac2[1], " effect", ")"), 
                                 paste0("phi_21(", lev.fac2[2], "*", lev.fac1[2], " effect", ")"), paste0("phi_22(", lev.fac2[1], "*", lev.fac1[2], " effect", ")"))
    
    M6_f2l2_f1l2$sp_effects <- c(paste0("phi_11(General[", lev.fac2[2], "-", lev.fac1[2], "])"), paste0("phi_12(", lev.fac2[1], " effect", ")"), 
                                 paste0("phi_21(", lev.fac2[2], "*", lev.fac1[1], " effect", ")"), paste0("phi_22(", lev.fac2[1], "*", lev.fac1[1], " effect", ")"))
  
  ## Save models into list
  data.models <- list(M0, M1, # M0 & M1
                      M2_f1l1_f2l1, M2_f1l2_f2l1, M2_f1l1_f2l2, M2_f1l2_f2l2, # M2
                      M3_f1l1_f2l1, M3_f1l2_f2l1, # M3
                      M4_f1l1_f2l1, M4_f1l1_f2l2, # M4
                      M5_f1l1_f2l1, M5_f1l2_f2l1, M5_f1l1_f2l2, M5_f1l2_f2l2, # M5
                      M6_f1l1_f2l1, M6_f1l2_f2l1, M6_f1l1_f2l2, M6_f1l2_f2l2, # M6
                      M6_f2l1_f1l1, M6_f2l2_f1l1, M6_f2l1_f1l2, M6_f2l2_f1l2) # M6
  
  ## Add name to each model
  names(data.models)[1] <- "M0" 
  names(data.models)[2] <- "M1" 
  
  names(data.models)[3] <- paste0("M2.(", lev.fac1[1],")") 
  names(data.models)[4] <- paste0("M2.(", lev.fac1[2],")") 
  names(data.models)[5] <- paste0("M2.(", lev.fac2[1],")") 
  names(data.models)[6] <- paste0("M2.(", lev.fac2[2],")") 
  
  names(data.models)[7] <- paste0("M3.", fac.names[1], "(", lev.fac1[1],")") 
  names(data.models)[8] <- paste0("M3.", fac.names[1], "(", lev.fac1[2],")") 
  
  names(data.models)[9] <- paste0("M4.", fac.names[2], "(", lev.fac2[1],")") 
  names(data.models)[10] <- paste0("M4.", fac.names[2], "(", lev.fac2[2],")") 
  
  names(data.models)[11] <- paste0("M5.", fac.names[1], "(", lev.fac1[1],")", "+", fac.names[2], "(", lev.fac2[1],")")
  names(data.models)[12] <- paste0("M5.", fac.names[1], "(", lev.fac1[2],")", "+", fac.names[2], "(", lev.fac2[1],")")
  names(data.models)[13] <- paste0("M5.", fac.names[1], "(", lev.fac1[1],")", "+", fac.names[2], "(", lev.fac2[2],")")
  names(data.models)[14] <- paste0("M5.", fac.names[1], "(", lev.fac1[2],")", "+", fac.names[2], "(", lev.fac2[2],")")
  
  names(data.models)[15] <- paste0("M6.", fac.names[1], "(", lev.fac1[1],")", "*", fac.names[2], "(", lev.fac2[1],")")
  names(data.models)[16] <- paste0("M6.", fac.names[1], "(", lev.fac1[2],")", "*", fac.names[2], "(", lev.fac2[1],")")
  names(data.models)[17] <- paste0("M6.", fac.names[1], "(", lev.fac1[1],")", "*", fac.names[2], "(", lev.fac2[2],")")
  names(data.models)[18] <- paste0("M6.", fac.names[1], "(", lev.fac1[2],")", "*", fac.names[2], "(", lev.fac2[2],")")
  
  names(data.models)[19] <- paste0("M6.", fac.names[2], "(", lev.fac2[1],")", "*", fac.names[1], "(", lev.fac1[1],")")
  names(data.models)[20] <- paste0("M6.", fac.names[2], "(", lev.fac2[2],")", "*", fac.names[1], "(", lev.fac1[1],")")
  names(data.models)[21] <- paste0("M6.", fac.names[2], "(", lev.fac2[1],")", "*", fac.names[1], "(", lev.fac1[2],")")
  names(data.models)[22] <- paste0("M6.", fac.names[2], "(", lev.fac2[2],")", "*", fac.names[1], "(", lev.fac1[2],")")
  
  data_summary <- data.frame(
    "MODEL"=c("M0", "M1", 
              paste0("M2-ind(", lev.fac1[1], "-", lev.fac2[1], ")"), paste0("M2-ind(", lev.fac1[2], "-", lev.fac2[1], ")"), 
              paste0("M2-ind(", lev.fac2[1], "-", lev.fac1[2], ")"), paste0("M2-ind(", lev.fac2[2], "-", lev.fac1[2], ")"),
              
              paste0("M3-", fac.names[1], ".(", lev.fac1[1],")"), paste0("M3-", fac.names[1], ".(", lev.fac1[2],")"),
              
              paste0("M4-", fac.names[2], ".(", lev.fac2[1],")"), paste0("M4-", fac.names[2], ".(", lev.fac2[2],")"), 
              
              paste0("M5-", fac.names[1], ".(", lev.fac1[1],")", "+", fac.names[2], ".(", lev.fac2[1],")"),
              paste0("M5-", fac.names[1], ".(", lev.fac1[2],")", "+", fac.names[2], ".(", lev.fac2[1],")"), 
              paste0("M5-", fac.names[1], ".(", lev.fac1[1],")", "+", fac.names[2], ".(", lev.fac2[2],")"),
              paste0("M5-", fac.names[1], ".(", lev.fac1[2],")", "+", fac.names[2], ".(", lev.fac2[2],")"),
              
              paste0("M6-", fac.names[1], ".(", lev.fac1[1],")", "*", fac.names[2], ".(", lev.fac2[1],")"),
              paste0("M6-", fac.names[1], ".(", lev.fac1[2],")", "*", fac.names[2], ".(", lev.fac2[1],")"),
              paste0("M6-", fac.names[1], ".(", lev.fac1[1],")", "*", fac.names[2], ".(", lev.fac2[2],")"),
              paste0("M6-", fac.names[1], ".(", lev.fac1[2],")", "*", fac.names[2], ".(", lev.fac2[2],")"),
              
              paste0("M6-", fac.names[2], ".(", lev.fac2[1],")", "*", fac.names[1], ".(", lev.fac1[1],")"),
              paste0("M6-", fac.names[2], ".(", lev.fac2[2],")", "*", fac.names[1], ".(", lev.fac1[1],")"),
              paste0("M6-", fac.names[2], ".(", lev.fac2[1],")", "*", fac.names[1], ".(", lev.fac1[2],")"),
              paste0("M6-", fac.names[2], ".(", lev.fac2[2],")", "*", fac.names[1], ".(", lev.fac1[2],")")
             ), 
    
    "DIC"=c(round(M0$dic$dic, 1), # M0
            round(M1$dic$dic, 1), # M1
            round(M2_f1l1_f2l1$dic$dic, 1), round(M2_f1l2_f2l1$dic$dic, 1), round(M2_f1l1_f2l2$dic$dic, 1), round(M2_f1l2_f2l2$dic$dic, 1), # M2
            round(M3_f1l1_f2l1$dic$dic, 1), round(M3_f1l2_f2l1$dic$dic, 1), # M3
            round(M4_f1l1_f2l1$dic$dic, 1), round(M4_f1l1_f2l2$dic$dic, 1), # M4
            round(M5_f1l1_f2l1$dic$dic, 1), round(M5_f1l2_f2l1$dic$dic, 1), round(M5_f1l1_f2l2$dic$dic, 1), round(M5_f1l2_f2l2$dic$dic, 1), # M5
            round(M6_f1l1_f2l1$dic$dic, 1), round(M6_f1l2_f2l1$dic$dic, 1), round(M6_f1l1_f2l2$dic$dic, 1), round(M6_f1l2_f2l2$dic$dic, 1), # M6
            round(M6_f2l1_f1l1$dic$dic, 1), round(M6_f2l2_f1l1$dic$dic, 1), round(M6_f2l1_f1l2$dic$dic, 1), round(M6_f2l2_f1l2$dic$dic, 1) # M6
            ), 
    
    "WAIC"=c(round(M0$waic$waic, 1), # M0
             round(M1$waic$waic, 1), # M1
             round(M2_f1l1_f2l1$waic$waic, 1), round(M2_f1l2_f2l1$waic$waic, 1), round(M2_f1l1_f2l2$waic$waic, 1), round(M2_f1l2_f2l2$waic$waic, 1), # M2
             round(M3_f1l1_f2l1$waic$waic, 1), round(M3_f1l2_f2l1$waic$waic, 1), # M3
             round(M4_f1l1_f2l1$waic$waic, 1), round(M4_f1l1_f2l2$waic$waic, 1), # M4
             round(M5_f1l1_f2l1$waic$waic, 1), round(M5_f1l2_f2l1$waic$waic, 1), round(M5_f1l1_f2l2$waic$waic, 1), round(M5_f1l2_f2l2$waic$waic, 1), # M5
             round(M6_f1l1_f2l1$waic$waic, 1), round(M6_f1l2_f2l1$waic$waic, 1), round(M6_f1l1_f2l2$waic$waic, 1), round(M6_f1l2_f2l2$waic$waic, 1), # M6
             round(M6_f2l1_f1l1$waic$waic, 1), round(M6_f2l2_f1l1$waic$waic, 1), round(M6_f2l1_f1l2$waic$waic, 1), round(M6_f2l2_f1l2$waic$waic, 1) # M6
             ), 
    "CPU"=c(M0$cpu.used[4], # M0
            M1$cpu.used[4], # M1
            M2_f1l1_f2l1$cpu.used[4], M2_f1l2_f2l1$cpu.used[4], M2_f1l1_f2l2$cpu.used[4], M2_f1l2_f2l2$cpu.used[4], # M2
            M3_f1l1_f2l1$cpu.used[4], M3_f1l2_f2l1$cpu.used[4], # M3
            M4_f1l1_f2l1$cpu.used[4], M4_f1l1_f2l2$cpu.used[4], # M4
            M5_f1l1_f2l1$cpu.used[4], M5_f1l2_f2l1$cpu.used[4], M5_f1l1_f2l2$cpu.used[4], M5_f1l2_f2l2$cpu.used[4], # M5
            M6_f1l1_f2l1$cpu.used[4], M6_f1l2_f2l1$cpu.used[4], M6_f1l1_f2l2$cpu.used[4], M6_f1l2_f2l2$cpu.used[4], # M6
            M6_f2l1_f1l1$cpu.used[4], M6_f2l2_f1l1$cpu.used[4], M6_f2l1_f1l2$cpu.used[4], M6_f2l2_f1l2$cpu.used[4] # M6
            )
    )
  
  data.models[[23]] <- data_summary
  names(data.models)[23] <- "Summary"
  print(data_summary)
  
  # Devolvemos lista con todos los modelos
  return(data.models)
}
