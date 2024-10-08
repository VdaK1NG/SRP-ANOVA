# Load libraries
library("INLA")
library("spdep")
library("ggplot2")
library("data.table")
library("RColorBrewer")
library("broom")
library("tidyverse")
library("santoku")
library("gridExtra")
library("kableExtra")

# Customise ggplot
PSIC <- "#133BF2"
DIF <- c("#133BF2", "#7189F7", "#FFFFFF", "#FF867A", "#FF2F1B")
CScale_dif <- colorRampPalette(DIF)

tema_VDK <- theme_update(plot.title = element_text(size=18, face= "bold", colour= "grey43", hjust = 0.5), 
                         axis.title.x = element_text(size=14, face="bold", colour = "black"), 
                         axis.title.y = element_text(size=14, face="bold", colour = "black"), 
                         panel.border = element_rect(colour = "black", fill=NA, size=0.3), 
                         legend.title = element_text(size=10, face="bold", colour = "black"),
                         legend.text = element_text(size=8, colour = "black"),
                         strip.background = element_rect(color="black", fill=PSIC, size=0.5, linetype="solid"),
                         strip.text.x = element_text(color = "white", face = "bold"))

# Create function to fix map labels
newLabels <- function(x, dig.lab = 4){
  lev <- levels(x)
  pattern <- paste0("^[\\(\\[][-]*\\d*\\.\\d{", dig.lab, "}|,[-]*\\d*\\.\\d{", dig.lab, "}")
  m <- gregexpr(pattern = pattern, levels(x))
  y <- regmatches(lev, m)
  y <- sapply(y, paste, collapse = "")
  y <- paste0(y, substring(lev, nchar(lev)))
  y
}

# Load function
source("functions.R")

# Load data
load("./Datos/data_paper.Rdata")
load("dat_112_CV.Rdata")
mun_cv_df <- st_as_sf(mun_cv)

# Define values of the dataset
n_areas <- 542 # number of municipalities

# Run Models

  ## Period VS Caller 
  period.vs.caller <- inla.ShANOVA.2x2(obs = c(vic_pre, test_pre, vic_covid, test_covid), exp = EXP_112_call_covid, gr = g_ccaa_cv, 
                                      fac.names = c("Period", "Caller"), lev.fac1 = c("preCOVID", "postCOVID"), 
                                      lev.fac2 = c("Per. In Cri.", "Byst."))
  ## Caller VS Gender
  caller.vs.gender <- inla.ShANOVA.2x2(obs = c(h_vic, m_vic, h_test, m_test), exp = EXP_112_sex_call, gr = g_ccaa_cv, 
                                       fac.names = c("Caller", "Gender"), lev.fac1 = c("Per. In Cri.", "Byst."), 
                                       lev.fac2 = c("Male", "Female"))
  
  ## Period vs Gender
  period.vs.gender <- inla.ShANOVA.2x2(obs = c(h_pre, m_pre, h_covid, m_covid), exp = EXP_112_sex_covid, gr = g_ccaa_cv, 
                                      fac.names = c("Period", "Gender"), lev.fac1 = c("preCOVID", "postCOVID"), 
                                      lev.fac2 = c("Male", "Female"))
  
# Prepare data for figures
mods <- c("period.vs.caller", "period.vs.gender", "caller.vs.gender") # names of the objects
  
  ## Prepare Spatial Effects

    ### M1 -/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-
    for (i in 1:length(mods)) {
      
      data.temp <- get(mods[i])
      data.mod <- data.temp[[2]]
      data_g1 <- data.mod$summary.random$phi_1$mean
      data_g2 <- data.mod$summary.random$phi_2$mean
      data_g3 <- data.mod$summary.random$phi_3$mean
      data_g4 <-data.mod$summary.random$phi_4$mean
      
      mun_cv_df <- cbind(mun_cv_df, data_g1)
      mun_cv_df <- cbind(mun_cv_df, data_g2)
      mun_cv_df <- cbind(mun_cv_df, data_g3)
      mun_cv_df <- cbind(mun_cv_df, data_g4)
      
      names(mun_cv_df)[names(mun_cv_df) == "data_g1"] <- paste0(mods[i], "_mod", 2, ".sp.1")
      names(mun_cv_df)[names(mun_cv_df) == "data_g2"] <- paste0(mods[i], "_mod", 2, ".sp.2")
      names(mun_cv_df)[names(mun_cv_df) == "data_g3"] <- paste0(mods[i], "_mod", 2, ".sp.3")
      names(mun_cv_df)[names(mun_cv_df) == "data_g4"] <- paste0(mods[i], "_mod", 2, ".sp.4")
      
    }
  
    ### M2 -/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-
    for (i in 1:length(mods)) {
      for (j in 3:6) {
        data.temp <- get(mods[i])
        data.mod <- data.temp[[j]]
        data_g1 <- data.mod$summary.random$phi_11$mean
        data_g2 <- NA
        data_g3 <- NA
        data_g4 <- NA
        
        mun_cv_df <- cbind(mun_cv_df, data_g1)
        mun_cv_df <- cbind(mun_cv_df, data_g2)
        mun_cv_df <- cbind(mun_cv_df, data_g3)
        mun_cv_df <- cbind(mun_cv_df, data_g4)
        
        names(mun_cv_df)[names(mun_cv_df) == "data_g1"] <- paste0(mods[i], "_mod", j, ".sp.1")
        names(mun_cv_df)[names(mun_cv_df) == "data_g2"] <- paste0(mods[i], "_mod", j, ".sp.2")
        names(mun_cv_df)[names(mun_cv_df) == "data_g3"] <- paste0(mods[i], "_mod", j, ".sp.3")
        names(mun_cv_df)[names(mun_cv_df) == "data_g4"] <- paste0(mods[i], "_mod", j, ".sp.4")
      }
    }
  
    ### M3 -/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-
    for (i in 1:length(mods)) {
      for (j in 7:8) {
        data.temp <- get(mods[i])
        data.mod <- data.temp[[j]]
        data_g1 <- data.mod$summary.random$phi_11$mean
        data_g2 <- data.mod$summary.random$phi_12$mean
        data_g3 <- NA
        data_g4 <- NA
        
        mun_cv_df <- cbind(mun_cv_df, data_g1)
        mun_cv_df <- cbind(mun_cv_df, data_g2)
        mun_cv_df <- cbind(mun_cv_df, data_g3)
        mun_cv_df <- cbind(mun_cv_df, data_g4)
        
        names(mun_cv_df)[names(mun_cv_df) == "data_g1"] <- paste0(mods[i], "_mod", j, ".sp.1")
        names(mun_cv_df)[names(mun_cv_df) == "data_g2"] <- paste0(mods[i], "_mod", j, ".sp.2")
        names(mun_cv_df)[names(mun_cv_df) == "data_g3"] <- paste0(mods[i], "_mod", j, ".sp.3")
        names(mun_cv_df)[names(mun_cv_df) == "data_g4"] <- paste0(mods[i], "_mod", j, ".sp.4")
      }
    }
  
    ### M4 -/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-
    for (i in 1:length(mods)) {
      for (j in 9:10) {
        data.temp <- get(mods[i])
        data.mod <- data.temp[[j]]
        data_g1 <- data.mod$summary.random$phi_11$mean
        data_g2 <- data.mod$summary.random$phi_21$mean
        data_g3 <- NA
        data_g4 <- NA
        
        mun_cv_df <- cbind(mun_cv_df, data_g1)
        mun_cv_df <- cbind(mun_cv_df, data_g2)
        mun_cv_df <- cbind(mun_cv_df, data_g3)
        mun_cv_df <- cbind(mun_cv_df, data_g4)
        
        names(mun_cv_df)[names(mun_cv_df) == "data_g1"] <- paste0(mods[i], "_mod", j, ".sp.1")
        names(mun_cv_df)[names(mun_cv_df) == "data_g2"] <- paste0(mods[i], "_mod", j, ".sp.2")
        names(mun_cv_df)[names(mun_cv_df) == "data_g3"] <- paste0(mods[i], "_mod", j, ".sp.3")
        names(mun_cv_df)[names(mun_cv_df) == "data_g4"] <- paste0(mods[i], "_mod", j, ".sp.4")
      }
    }
  
    ### M5 -/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-
    for (i in 1:length(mods)) {
      for (j in 11:14) {
        data.temp <- get(mods[i])
        data.mod <- data.temp[[j]]
        data_g1 <- data.mod$summary.random$phi_11$mean
        data_g2 <- data.mod$summary.random$phi_12$mean
        data_g3 <- data.mod$summary.random$phi_21$mean
        data_g4 <- NA
        
        mun_cv_df <- cbind(mun_cv_df, data_g1)
        mun_cv_df <- cbind(mun_cv_df, data_g2)
        mun_cv_df <- cbind(mun_cv_df, data_g3)
        mun_cv_df <- cbind(mun_cv_df, data_g4)
        
        names(mun_cv_df)[names(mun_cv_df) == "data_g1"] <- paste0(mods[i], "_mod", j, ".sp.1")
        names(mun_cv_df)[names(mun_cv_df) == "data_g2"] <- paste0(mods[i], "_mod", j, ".sp.2")
        names(mun_cv_df)[names(mun_cv_df) == "data_g3"] <- paste0(mods[i], "_mod", j, ".sp.3")
        names(mun_cv_df)[names(mun_cv_df) == "data_g4"] <- paste0(mods[i], "_mod", j, ".sp.4")
      }
    }
  
    ### M6 -/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-
    for (i in 1:length(mods)) {
      for (j in 15:22) {
        data.temp <- get(mods[i])
        data.mod <- data.temp[[j]]
        data_g1 <- data.mod$summary.random$phi_11$mean
        data_g2 <- data.mod$summary.random$phi_12$mean
        data_g3 <- data.mod$summary.random$phi_21$mean
        data_g4 <- data.mod$summary.random$phi_22$mean
        
        mun_cv_df <- cbind(mun_cv_df, data_g1)
        mun_cv_df <- cbind(mun_cv_df, data_g2)
        mun_cv_df <- cbind(mun_cv_df, data_g3)
        mun_cv_df <- cbind(mun_cv_df, data_g4)
        
        names(mun_cv_df)[names(mun_cv_df) == "data_g1"] <- paste0(mods[i], "_mod", j, ".sp.1")
        names(mun_cv_df)[names(mun_cv_df) == "data_g2"] <- paste0(mods[i], "_mod", j, ".sp.2")
        names(mun_cv_df)[names(mun_cv_df) == "data_g3"] <- paste0(mods[i], "_mod", j, ".sp.3")
        names(mun_cv_df)[names(mun_cv_df) == "data_g4"] <- paste0(mods[i], "_mod", j, ".sp.4")
      }
    }
  
  ## Prepare Relative Risks Maps
  for (i in 1:length(mods)) {
    
    for (j in 1:22) {
      data.temp <- get(mods[i])
      data.mod <- data.temp[[j]]
      data.mod <- data.mod$summary.fitted.values$mean
      
      data_g1 <- data.mod[1:n_areas]
      data_g2 <- data.mod[(n_areas+1):(n_areas*2)]
      data_g3 <- data.mod[((n_areas*2)+1):(n_areas*3)]
      data_g4 <- data.mod[((n_areas*3)+1):(n_areas*4)]
      
      mun_cv_df <- cbind(mun_cv_df, data_g1)
      mun_cv_df <- cbind(mun_cv_df, data_g2)
      mun_cv_df <- cbind(mun_cv_df, data_g3)
      mun_cv_df <- cbind(mun_cv_df, data_g4)
      
      names(mun_cv_df)[names(mun_cv_df) == "data_g1"] <- paste0(mods[i], "_mod", j, ".RME.G1")
      names(mun_cv_df)[names(mun_cv_df) == "data_g2"] <- paste0(mods[i], "_mod", j, ".RME.G2")
      names(mun_cv_df)[names(mun_cv_df) == "data_g3"] <- paste0(mods[i], "_mod", j, ".RME.G3")
      names(mun_cv_df)[names(mun_cv_df) == "data_g4"] <- paste0(mods[i], "_mod", j, ".RME.G4")
      
    }
    
  }
  
  ## Add names to groups
  fig.names <- data.frame("mods" = rep(mods, each=22), "n_mod" = rep(1:22, length(mods)), "group1" = NA, "group2" = NA, "group3" = NA, "group4" = NA, 
                          "sp_ef_1"="Not Included", "sp_ef_2"="Not Included", "sp_ef_3"="Not Included", "sp_ef_4"="Not Included")
  
  for (i in 1:nrow(fig.names)) {
    
    data.temp <- get(fig.names$mods[i])
    data.mod <- data.temp[[fig.names$n_mod[i]]]
    fig.names$group1[i] <- data.mod$groups[1]
    fig.names$group2[i] <- data.mod$groups[2]
    fig.names$group3[i] <- data.mod$groups[3]
    fig.names$group4[i] <- data.mod$groups[4]
    
  }
  
  ## Add names to effects
  for (i in 1:nrow(fig.names)) {
    
    if(fig.names$n_mod[i] %in% c(3:6)){
      
      data.temp <- get(fig.names$mods[i])
      data.mod <- data.temp[[fig.names$n_mod[i]]]
      fig.names$sp_ef_1[i] <- data.mod$sp_effects[1]
      
    }else if(fig.names$n_mod[i] %in% c(7:10)){
      
      data.temp <- get(fig.names$mods[i])
      data.mod <- data.temp[[fig.names$n_mod[i]]]
      fig.names$sp_ef_1[i] <- data.mod$sp_effects[1]
      fig.names$sp_ef_2[i] <- data.mod$sp_effects[2]
      
    }else if(fig.names$n_mod[i] %in% c(11:14)){
      
      data.temp <- get(fig.names$mods[i])
      data.mod <- data.temp[[fig.names$n_mod[i]]]
      fig.names$sp_ef_1[i] <- data.mod$sp_effects[1]
      fig.names$sp_ef_2[i] <- data.mod$sp_effects[2]
      fig.names$sp_ef_3[i] <- data.mod$sp_effects[3]
      
    }else if(fig.names$n_mod[i] %in% c(2, 15:22)){
      
      data.temp <- get(fig.names$mods[i])
      data.mod <- data.temp[[fig.names$n_mod[i]]]
      fig.names$sp_ef_1[i] <- data.mod$sp_effects[1]
      fig.names$sp_ef_2[i] <- data.mod$sp_effects[2]
      fig.names$sp_ef_3[i] <- data.mod$sp_effects[3]
      fig.names$sp_ef_4[i] <- data.mod$sp_effects[4]
    }
  }
  
# Figures - Spatial Effects
  
  ## Period vs Gender -/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-
  
  ### Define model to paint
  input_factors <- "period.vs.gender"
  input_nmod <- 11
  
  ### Select model to paint
  col_numg1 <- which(names(mun_cv_df) == paste0(input_factors, "_mod", input_nmod, ".sp.1"))
  col_numg2 <- which(names(mun_cv_df) == paste0(input_factors, "_mod", input_nmod, ".sp.2"))
  col_numg3 <- which(names(mun_cv_df) == paste0(input_factors, "_mod", input_nmod, ".sp.3"))
  col_numg4 <- which(names(mun_cv_df) == paste0(input_factors, "_mod", input_nmod, ".sp.4"))
  
  title_sp_plot1 <- fig.names %>% filter(mods == input_factors, n_mod == input_nmod) %>% select(sp_ef_1) %>% pull()
  title_sp_plot2 <- fig.names %>% filter(mods == input_factors, n_mod == input_nmod) %>% select(sp_ef_2) %>% pull()
  title_sp_plot3 <- fig.names %>% filter(mods == input_factors, n_mod == input_nmod) %>% select(sp_ef_3) %>% pull()
  title_sp_plot4 <- fig.names %>% filter(mods == input_factors, n_mod == input_nmod) %>% select(sp_ef_4) %>% pull()
  
  # Fill values for group 1
  mun_cv_df$fill_val1 <- cut_number(as.data.frame(mun_cv_df[,col_numg1])[,1], 5)
  mun_cv_df$fill_val1 <- factor(mun_cv_df$fill_val1, labels = newLabels(mun_cv_df$fill_val1, dig.lab = 2))
  
  # Fill values for group 2
  mun_cv_df$fill_val2 <- cut_number(as.data.frame(mun_cv_df[,col_numg2])[,1], 5)
  mun_cv_df$fill_val2 <- factor(mun_cv_df$fill_val2, labels = newLabels(mun_cv_df$fill_val2, dig.lab = 2))
  
  # Fill values for group 3
  mun_cv_df$fill_val3 <- cut_number(as.data.frame(mun_cv_df[,col_numg3])[,1], 5)
  mun_cv_df$fill_val3 <- factor(mun_cv_df$fill_val3, labels = newLabels(mun_cv_df$fill_val3, dig.lab = 2))
  
  fig1 <- ggplot() + 
    geom_sf(color = "black", linewidth = 0.25, data  = mun_cv_df, aes(fill = fill_val1)) + 
    scale_fill_manual("", values = CScale_dif(5), na.value="darkgray", 
                      guide = guide_legend(direction = "horizontal", nrow = 3, keywidth = 1)) +
    theme(axis.text = element_blank(), axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank(), 
          legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=18), 
          strip.text.x = element_text(size = 12, face = "bold.italic")) +  scale_y_continuous(expand = c(0, 0)) + 
    ggtitle(title_sp_plot1) + xlab("") + ylab("")
  
  fig2 <- ggplot() + 
    geom_sf(color = "black", linewidth = 0.25, data  = mun_cv_df, aes(fill = fill_val2)) + 
    scale_fill_manual("", values = CScale_dif(5), na.value="darkgray", 
                      guide = guide_legend(direction = "horizontal", nrow = 3, keywidth = 1)) +
    theme(axis.text = element_blank(), axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank(), 
          legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=18), 
          strip.text.x = element_text(size = 12, face = "bold.italic")) +  scale_y_continuous(expand = c(0, 0)) + 
    ggtitle(title_sp_plot2) + xlab("") + ylab("")
  
  fig3 <- ggplot() + 
    geom_sf(color = "black", linewidth = 0.25, data  = mun_cv_df, aes(fill = fill_val3)) + 
    scale_fill_manual("", values = CScale_dif(5), na.value="darkgray", 
                      guide = guide_legend(direction = "horizontal", nrow = 3, keywidth = 1)) +
    theme(axis.text = element_blank(), axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank(), 
          legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=18), 
          strip.text.x = element_text(size = 12, face = "bold.italic")) +  scale_y_continuous(expand = c(0, 0)) + 
    ggtitle(title_sp_plot3) + xlab("") + ylab("")
  
  fig4 <- ggplot() + 
    geom_sf(color = "black", linewidth = 0.25, data  = mun_cv_df, aes(fill = NA)) + 
    scale_fill_manual("", values = CScale_dif(5), na.value="darkgray", 
                      guide = guide_legend(direction = "horizontal", nrow = 3, keywidth = 1)) +
    theme(axis.text = element_blank(), axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank(), 
          legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=18), 
          strip.text.x = element_text(size = 12, face = "bold.italic")) +  scale_y_continuous(expand = c(0, 0)) + 
    ggtitle(title_sp_plot4) + xlab("") + ylab("")
  
  grid.arrange(fig1, fig2, fig3, fig4, ncol=4)
  
  dev.print(device = png, filename = "fig1.png", width = 1400, height = 700)
  
  ## Caller VS Period -/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-
  
    ### Define model to paint
    input_factors <- "period.vs.caller"
    input_nmod <- 22
    
    ### Select model to paint
    col_numg1 <- which(names(mun_cv_df) == paste0(input_factors, "_mod", input_nmod, ".sp.1"))
    col_numg2 <- which(names(mun_cv_df) == paste0(input_factors, "_mod", input_nmod, ".sp.2"))
    col_numg3 <- which(names(mun_cv_df) == paste0(input_factors, "_mod", input_nmod, ".sp.3"))
    col_numg4 <- which(names(mun_cv_df) == paste0(input_factors, "_mod", input_nmod, ".sp.4"))
    
    title_sp_plot1 <- fig.names %>% filter(mods == input_factors, n_mod == input_nmod) %>% select(sp_ef_1) %>% pull()
    title_sp_plot2 <- fig.names %>% filter(mods == input_factors, n_mod == input_nmod) %>% select(sp_ef_2) %>% pull()
    title_sp_plot3 <- fig.names %>% filter(mods == input_factors, n_mod == input_nmod) %>% select(sp_ef_3) %>% pull()
    title_sp_plot4 <- fig.names %>% filter(mods == input_factors, n_mod == input_nmod) %>% select(sp_ef_4) %>% pull()
    
    ### Fill values for group 1
    mun_cv_df$fill_val1 <- cut_number(as.data.frame(mun_cv_df[,col_numg1])[,1], 5)
    mun_cv_df$fill_val1 <- factor(mun_cv_df$fill_val1, labels = newLabels(mun_cv_df$fill_val1, dig.lab = 2))
    
    ### Fill values for group 2
    mun_cv_df$fill_val2 <- cut_number(as.data.frame(mun_cv_df[,col_numg2])[,1], 5)
    mun_cv_df$fill_val2 <- factor(mun_cv_df$fill_val2, labels = newLabels(mun_cv_df$fill_val2, dig.lab = 2))
    
    ### Fill values for group 3
    mun_cv_df$fill_val3 <- cut_number(as.data.frame(mun_cv_df[,col_numg3])[,1], 5)
    mun_cv_df$fill_val3 <- factor(mun_cv_df$fill_val3, labels = newLabels(mun_cv_df$fill_val3, dig.lab = 2))
    
    ### Fill values for group 4
    mun_cv_df$fill_val4 <- cut_number(as.data.frame(mun_cv_df[,col_numg4])[,1], 5)
    mun_cv_df$fill_val4 <- factor(mun_cv_df$fill_val4, labels = newLabels(mun_cv_df$fill_val4, dig.lab = 2))
    
    fig1 <- ggplot() + 
      geom_sf(color = "black", linewidth = 0.25, data  = mun_cv_df, aes(fill = fill_val1)) + 
      scale_fill_manual("", values = CScale_dif(5), na.value="darkgray", 
                        guide = guide_legend(direction = "horizontal", nrow = 3, keywidth = 1)) +
      theme(axis.text = element_blank(), axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank(), 
            legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=18), 
            strip.text.x = element_text(size = 12, face = "bold.italic")) +  scale_y_continuous(expand = c(0, 0)) + 
      ggtitle(title_sp_plot1) + xlab("") + ylab("")
    
    fig2 <- ggplot() + 
      geom_sf(color = "black", linewidth = 0.25, data  = mun_cv_df, aes(fill = fill_val2)) + 
      scale_fill_manual("", values = CScale_dif(5), na.value="darkgray", 
                        guide = guide_legend(direction = "horizontal", nrow = 3, keywidth = 1)) +
      theme(axis.text = element_blank(), axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank(), 
            legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=18), 
            strip.text.x = element_text(size = 12, face = "bold.italic")) +  scale_y_continuous(expand = c(0, 0)) + 
      ggtitle(title_sp_plot2) + xlab("") + ylab("")
    
    fig3 <- ggplot() + 
      geom_sf(color = "black", linewidth = 0.25, data  = mun_cv_df, aes(fill = fill_val3)) + 
      scale_fill_manual("", values = CScale_dif(5), na.value="darkgray", 
                        guide = guide_legend(direction = "horizontal", nrow = 3, keywidth = 1)) +
      theme(axis.text = element_blank(), axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank(), 
            legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=18), 
            strip.text.x = element_text(size = 12, face = "bold.italic")) +  scale_y_continuous(expand = c(0, 0)) + 
      ggtitle(title_sp_plot3) + xlab("") + ylab("")
    
    fig4 <- ggplot() + 
      geom_sf(color = "black", linewidth = 0.25, data  = mun_cv_df, aes(fill = fill_val4)) + 
      scale_fill_manual("", values = CScale_dif(5), na.value="darkgray", 
                        guide = guide_legend(direction = "horizontal", nrow = 3, keywidth = 1)) +
      theme(axis.text = element_blank(), axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank(), 
            legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=18), 
            strip.text.x = element_text(size = 12, face = "bold.italic")) +  scale_y_continuous(expand = c(0, 0)) + 
      ggtitle(title_sp_plot4) + xlab("") + ylab("")
    
    grid.arrange(fig1, fig2, fig3, fig4, ncol=4)
    
    dev.print(device = png, filename = "fig3.png", width = 1400, height = 700)
  
  ## Caller VS Gender -/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-
    
    ### Define model to paint
    input_factors <- "caller.vs.gender"
    input_nmod <- 21
    
    ### Select model to paint
    col_numg1 <- which(names(mun_cv_df) == paste0(input_factors, "_mod", input_nmod, ".sp.1"))
    col_numg2 <- which(names(mun_cv_df) == paste0(input_factors, "_mod", input_nmod, ".sp.2"))
    col_numg3 <- which(names(mun_cv_df) == paste0(input_factors, "_mod", input_nmod, ".sp.3"))
    col_numg4 <- which(names(mun_cv_df) == paste0(input_factors, "_mod", input_nmod, ".sp.4"))
    
    title_sp_plot1 <- fig.names %>% filter(mods == input_factors, n_mod == input_nmod) %>% select(sp_ef_1) %>% pull()
    title_sp_plot2 <- fig.names %>% filter(mods == input_factors, n_mod == input_nmod) %>% select(sp_ef_2) %>% pull()
    title_sp_plot3 <- fig.names %>% filter(mods == input_factors, n_mod == input_nmod) %>% select(sp_ef_3) %>% pull()
    title_sp_plot4 <- fig.names %>% filter(mods == input_factors, n_mod == input_nmod) %>% select(sp_ef_4) %>% pull()
    
    ### Fill values for group 1
    mun_cv_df$fill_val1 <- cut_number(as.data.frame(mun_cv_df[,col_numg1])[,1], 5)
    mun_cv_df$fill_val1 <- factor(mun_cv_df$fill_val1, labels = newLabels(mun_cv_df$fill_val1, dig.lab = 2))
    
    ### Fill values for group 2
    mun_cv_df$fill_val2 <- cut_number(as.data.frame(mun_cv_df[,col_numg2])[,1], 5)
    mun_cv_df$fill_val2 <- factor(mun_cv_df$fill_val2, labels = newLabels(mun_cv_df$fill_val2, dig.lab = 2))
    
    ### Fill values for group 3
    mun_cv_df$fill_val3 <- cut_number(as.data.frame(mun_cv_df[,col_numg3])[,1], 5)
    mun_cv_df$fill_val3 <- factor(mun_cv_df$fill_val3, labels = newLabels(mun_cv_df$fill_val3, dig.lab = 2))
    
    ### Fill values for group 4
    mun_cv_df$fill_val4 <- cut_number(as.data.frame(mun_cv_df[,col_numg4])[,1], 5)
    mun_cv_df$fill_val4 <- factor(mun_cv_df$fill_val4, labels = newLabels(mun_cv_df$fill_val4, dig.lab = 2))
    
    fig1 <- ggplot() + 
      geom_sf(color = "black", linewidth = 0.25, data  = mun_cv_df, aes(fill = fill_val1)) + 
      scale_fill_manual("", values = CScale_dif(5), na.value="darkgray", 
                        guide = guide_legend(direction = "horizontal", nrow = 3, keywidth = 1)) +
      theme(axis.text = element_blank(), axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank(), 
            legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=18), 
            strip.text.x = element_text(size = 12, face = "bold.italic")) +  scale_y_continuous(expand = c(0, 0)) + 
      ggtitle(title_sp_plot1) + xlab("") + ylab("")
    
    fig2 <- ggplot() + 
      geom_sf(color = "black", linewidth = 0.25, data  = mun_cv_df, aes(fill = fill_val2)) + 
      scale_fill_manual("", values = CScale_dif(5), na.value="darkgray", 
                        guide = guide_legend(direction = "horizontal", nrow = 3, keywidth = 1)) +
      theme(axis.text = element_blank(), axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank(), 
            legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=18), 
            strip.text.x = element_text(size = 12, face = "bold.italic")) +  scale_y_continuous(expand = c(0, 0)) + 
      ggtitle(title_sp_plot2) + xlab("") + ylab("")
    
    fig3 <- ggplot() + 
      geom_sf(color = "black", linewidth = 0.25, data  = mun_cv_df, aes(fill = fill_val3)) + 
      scale_fill_manual("", values = CScale_dif(5), na.value="darkgray", 
                        guide = guide_legend(direction = "horizontal", nrow = 3, keywidth = 1)) +
      theme(axis.text = element_blank(), axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank(), 
            legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=18), 
            strip.text.x = element_text(size = 12, face = "bold.italic")) +  scale_y_continuous(expand = c(0, 0)) + 
      ggtitle(title_sp_plot3) + xlab("") + ylab("")
    
    fig4 <- ggplot() + 
      geom_sf(color = "black", linewidth = 0.25, data  = mun_cv_df, aes(fill = fill_val4)) + 
      scale_fill_manual("", values = CScale_dif(5), na.value="darkgray", 
                        guide = guide_legend(direction = "horizontal", nrow = 3, keywidth = 1)) +
      theme(axis.text = element_blank(), axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank(), 
            legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=18), 
            strip.text.x = element_text(size = 12, face = "bold.italic")) +  scale_y_continuous(expand = c(0, 0)) + 
      ggtitle(title_sp_plot4) + xlab("") + ylab("")
    
    grid.arrange(fig1, fig2, fig3, fig4, ncol=4)
    
    dev.print(device = png, filename = "fig5.png", width = 1400, height = 700)
  
# Figures - Relative Risks
    
    ## Period vs Gender
    
    ### Define model to paint
    input_factors <- "period.vs.gender"
    input_nmod <- 11
    
    ### Select mod to paint
    col_numg1 <- which(names(mun_cv_df) == paste0(input_factors, "_mod", input_nmod, ".RME.G1"))
    col_numg2 <- which(names(mun_cv_df) == paste0(input_factors, "_mod", input_nmod, ".RME.G2"))
    col_numg3 <- which(names(mun_cv_df) == paste0(input_factors, "_mod", input_nmod, ".RME.G3"))
    col_numg4 <- which(names(mun_cv_df) == paste0(input_factors, "_mod", input_nmod, ".RME.G4"))
    
    title_rme_plot1 <- fig.names %>% filter(mods == input_factors, n_mod == input_nmod) %>% select(group1) %>% pull()
    title_rme_plot2 <- fig.names %>% filter(mods == input_factors, n_mod == input_nmod) %>% select(group2) %>% pull()
    title_rme_plot3 <- fig.names %>% filter(mods == input_factors, n_mod == input_nmod) %>% select(group3) %>% pull()
    title_rme_plot4 <- fig.names %>% filter(mods == input_factors, n_mod == input_nmod) %>% select(group4) %>% pull()
    
    # Fill values for group 1
    mun_cv_df$fill_rme1 <- cut(as.data.frame(mun_cv_df[,col_numg1])[,1], breaks = c(0, 0.5, 0.9, 1.1, 2, Inf), 
                               labels = c("0, 0.5", "0.5, 0.9", "0.9, 1.1", "1.1, 2", "> 2"), include.lowest = TRUE)
    
    # Fill values for group 2
    mun_cv_df$fill_rme2 <- cut(as.data.frame(mun_cv_df[,col_numg2])[,1], breaks = c(0, 0.5, 0.9, 1.1, 2, Inf), 
                               labels = c("0, 0.5", "0.5, 0.9", "0.9, 1.1", "1.1, 2", "> 2"), include.lowest = TRUE)
    
    # Fill values for group 3
    mun_cv_df$fill_rme3 <- cut(as.data.frame(mun_cv_df[,col_numg3])[,1], breaks = c(0, 0.5, 0.9, 1.1, 2, Inf), 
                               labels = c("0, 0.5", "0.5, 0.9", "0.9, 1.1", "1.1, 2", "> 2"), include.lowest = TRUE)
    
    # Fill values for group 4
    mun_cv_df$fill_rme4 <- cut(as.data.frame(mun_cv_df[,col_numg4])[,1], breaks = c(0, 0.5, 0.9, 1.1, 2, Inf), 
                               labels = c("0, 0.5", "0.5, 0.9", "0.9, 1.1", "1.1, 2", "> 2"), include.lowest = TRUE)
    
    fig1 <- ggplot() + 
      geom_sf( color = "black", linewidth = 0.25, data  = mun_cv_df, aes(fill = fill_rme1)) + 
      scale_fill_manual("", values = CScale_dif(5), na.value="red", 
                        guide = guide_legend(direction = "horizontal", nrow = 3, keywidth = 1)) +
      theme(axis.text = element_blank(), axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank(), 
            legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=18), 
            strip.text.x = element_text(size = 12, face = "bold.italic")) +  scale_y_continuous(expand = c(0, 0)) + 
      ggtitle(title_rme_plot1) + xlab("") + ylab("")
    
    fig2 <- ggplot() + 
      geom_sf( color = "black", linewidth = 0.25, data  = mun_cv_df, aes(fill = fill_rme2)) + 
      scale_fill_manual("", values = CScale_dif(5), na.value="red", 
                        guide = guide_legend(direction = "horizontal", nrow = 3, keywidth = 1)) +
      theme(axis.text = element_blank(), axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank(), 
            legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=18), 
            strip.text.x = element_text(size = 12, face = "bold.italic")) +  scale_y_continuous(expand = c(0, 0)) + 
      ggtitle(title_rme_plot2) + xlab("") + ylab("")
    
    fig3 <- ggplot() + 
      geom_sf( color = "black", linewidth = 0.25, data  = mun_cv_df, aes(fill = fill_rme3)) + 
      scale_fill_manual("", values = CScale_dif(5), na.value="red", 
                        guide = guide_legend(direction = "horizontal", nrow = 3, keywidth = 1)) +
      theme(axis.text = element_blank(), axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank(), 
            legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=18), 
            strip.text.x = element_text(size = 12, face = "bold.italic")) +  scale_y_continuous(expand = c(0, 0)) + 
      ggtitle(title_rme_plot3) + xlab("") + ylab("")
    
    fig4 <- ggplot() + 
      geom_sf( color = "black", linewidth = 0.25, data  = mun_cv_df, aes(fill = fill_rme4)) + 
      scale_fill_manual("", values = CScale_dif(5), na.value="red", 
                        guide = guide_legend(direction = "horizontal", nrow = 3, keywidth = 1)) +
      theme(axis.text = element_blank(), axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank(), 
            legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=18), 
            strip.text.x = element_text(size = 12, face = "bold.italic")) +  scale_y_continuous(expand = c(0, 0)) + 
      ggtitle(title_rme_plot4) + xlab("") + ylab("")
    
    grid.arrange(fig1, fig2, fig3, fig4, ncol=4)
    
    dev.print(device = png, filename = "fig2.png", width = 1400, height = 700)
  
  ## Caller VS Period 
    
    ### Define model to paint
    input_factors <- "period.vs.caller"
    input_nmod <- 22
  
    ### Select mod to paint
    col_numg1 <- which(names(mun_cv_df) == paste0(input_factors, "_mod", input_nmod, ".RME.G1"))
    col_numg2 <- which(names(mun_cv_df) == paste0(input_factors, "_mod", input_nmod, ".RME.G2"))
    col_numg3 <- which(names(mun_cv_df) == paste0(input_factors, "_mod", input_nmod, ".RME.G3"))
    col_numg4 <- which(names(mun_cv_df) == paste0(input_factors, "_mod", input_nmod, ".RME.G4"))
    
    title_rme_plot1 <- fig.names %>% filter(mods == input_factors, n_mod == input_nmod) %>% select(group1) %>% pull()
    title_rme_plot2 <- fig.names %>% filter(mods == input_factors, n_mod == input_nmod) %>% select(group2) %>% pull()
    title_rme_plot3 <- fig.names %>% filter(mods == input_factors, n_mod == input_nmod) %>% select(group3) %>% pull()
    title_rme_plot4 <- fig.names %>% filter(mods == input_factors, n_mod == input_nmod) %>% select(group4) %>% pull()
    
    # Fill values for group 1
    mun_cv_df$fill_rme1 <- cut(as.data.frame(mun_cv_df[,col_numg1])[,1], breaks = c(0, 0.5, 0.9, 1.1, 2, Inf), 
                               labels = c("0, 0.5", "0.5, 0.9", "0.9, 1.1", "1.1, 2", "> 2"), include.lowest = TRUE)
    
    # Fill values for group 2
    mun_cv_df$fill_rme2 <- cut(as.data.frame(mun_cv_df[,col_numg2])[,1], breaks = c(0, 0.5, 0.9, 1.1, 2, Inf), 
                               labels = c("0, 0.5", "0.5, 0.9", "0.9, 1.1", "1.1, 2", "> 2"), include.lowest = TRUE)
    
    # Fill values for group 3
    mun_cv_df$fill_rme3 <- cut(as.data.frame(mun_cv_df[,col_numg3])[,1], breaks = c(0, 0.5, 0.9, 1.1, 2, Inf), 
                               labels = c("0, 0.5", "0.5, 0.9", "0.9, 1.1", "1.1, 2", "> 2"), include.lowest = TRUE)
    
    # Fill values for group 4
    mun_cv_df$fill_rme4 <- cut(as.data.frame(mun_cv_df[,col_numg4])[,1], breaks = c(0, 0.5, 0.9, 1.1, 2, Inf), 
                               labels = c("0, 0.5", "0.5, 0.9", "0.9, 1.1", "1.1, 2", "> 2"), include.lowest = TRUE)
    
    fig1 <- ggplot() + 
      geom_sf( color = "black", linewidth = 0.25, data  = mun_cv_df, aes(fill = fill_rme1)) + 
      scale_fill_manual("", values = CScale_dif(5), na.value="red", 
                        guide = guide_legend(direction = "horizontal", nrow = 3, keywidth = 1)) +
      theme(axis.text = element_blank(), axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank(), 
            legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=18), 
            strip.text.x = element_text(size = 12, face = "bold.italic")) +  scale_y_continuous(expand = c(0, 0)) + 
      ggtitle(title_rme_plot1) + xlab("") + ylab("")
    
    fig2 <- ggplot() + 
      geom_sf( color = "black", linewidth = 0.25, data  = mun_cv_df, aes(fill = fill_rme2)) + 
      scale_fill_manual("", values = CScale_dif(5), na.value="red", 
                        guide = guide_legend(direction = "horizontal", nrow = 3, keywidth = 1)) +
      theme(axis.text = element_blank(), axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank(), 
            legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=18), 
            strip.text.x = element_text(size = 12, face = "bold.italic")) +  scale_y_continuous(expand = c(0, 0)) + 
      ggtitle(title_rme_plot2) + xlab("") + ylab("")
    
    fig3 <- ggplot() + 
      geom_sf( color = "black", linewidth = 0.25, data  = mun_cv_df, aes(fill = fill_rme3)) + 
      scale_fill_manual("", values = CScale_dif(5), na.value="red", 
                        guide = guide_legend(direction = "horizontal", nrow = 3, keywidth = 1)) +
      theme(axis.text = element_blank(), axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank(), 
            legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=18), 
            strip.text.x = element_text(size = 12, face = "bold.italic")) +  scale_y_continuous(expand = c(0, 0)) + 
      ggtitle(title_rme_plot3) + xlab("") + ylab("")
    
    fig4 <- ggplot() + 
      geom_sf( color = "black", linewidth = 0.25, data  = mun_cv_df, aes(fill = fill_rme4)) + 
      scale_fill_manual("", values = CScale_dif(5), na.value="red", 
                        guide = guide_legend(direction = "horizontal", nrow = 3, keywidth = 1)) +
      theme(axis.text = element_blank(), axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank(), 
            legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=18), 
            strip.text.x = element_text(size = 12, face = "bold.italic")) +  scale_y_continuous(expand = c(0, 0)) + 
      ggtitle(title_rme_plot4) + xlab("") + ylab("")
    
    grid.arrange(fig1, fig2, fig3, fig4, ncol=4)
    
    dev.print(device = png, filename = "fig4.png", width = 1400, height = 700)
  
  ## Caller VS Gender
    
    ### Define model to paint
    input_factors <- "caller.vs.gender"
    input_nmod <- 21
    
    ### Select mod to paint
    col_numg1 <- which(names(mun_cv_df) == paste0(input_factors, "_mod", input_nmod, ".RME.G1"))
    col_numg2 <- which(names(mun_cv_df) == paste0(input_factors, "_mod", input_nmod, ".RME.G2"))
    col_numg3 <- which(names(mun_cv_df) == paste0(input_factors, "_mod", input_nmod, ".RME.G3"))
    col_numg4 <- which(names(mun_cv_df) == paste0(input_factors, "_mod", input_nmod, ".RME.G4"))
    
    title_rme_plot1 <- fig.names %>% filter(mods == input_factors, n_mod == input_nmod) %>% select(group1) %>% pull()
    title_rme_plot2 <- fig.names %>% filter(mods == input_factors, n_mod == input_nmod) %>% select(group2) %>% pull()
    title_rme_plot3 <- fig.names %>% filter(mods == input_factors, n_mod == input_nmod) %>% select(group3) %>% pull()
    title_rme_plot4 <- fig.names %>% filter(mods == input_factors, n_mod == input_nmod) %>% select(group4) %>% pull()
    
    # Fill values for group 1
    mun_cv_df$fill_rme1 <- cut(as.data.frame(mun_cv_df[,col_numg1])[,1], breaks = c(0, 0.5, 0.9, 1.1, 2, Inf), 
                               labels = c("0, 0.5", "0.5, 0.9", "0.9, 1.1", "1.1, 2", "> 2"), include.lowest = TRUE)
    
    # Fill values for group 2
    mun_cv_df$fill_rme2 <- cut(as.data.frame(mun_cv_df[,col_numg2])[,1], breaks = c(0, 0.5, 0.9, 1.1, 2, Inf), 
                               labels = c("0, 0.5", "0.5, 0.9", "0.9, 1.1", "1.1, 2", "> 2"), include.lowest = TRUE)
    
    # Fill values for group 3
    mun_cv_df$fill_rme3 <- cut(as.data.frame(mun_cv_df[,col_numg3])[,1], breaks = c(0, 0.5, 0.9, 1.1, 2, Inf), 
                               labels = c("0, 0.5", "0.5, 0.9", "0.9, 1.1", "1.1, 2", "> 2"), include.lowest = TRUE)
    
    # Fill values for group 4
    mun_cv_df$fill_rme4 <- cut(as.data.frame(mun_cv_df[,col_numg4])[,1], breaks = c(0, 0.5, 0.9, 1.1, 2, Inf), 
                               labels = c("0, 0.5", "0.5, 0.9", "0.9, 1.1", "1.1, 2", "> 2"), include.lowest = TRUE)
    
    fig1 <- ggplot() + 
      geom_sf( color = "black", linewidth = 0.25, data  = mun_cv_df, aes(fill = fill_rme1)) + 
      scale_fill_manual("", values = CScale_dif(5), na.value="red", 
                        guide = guide_legend(direction = "horizontal", nrow = 3, keywidth = 1)) +
      theme(axis.text = element_blank(), axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank(), 
            legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=18), 
            strip.text.x = element_text(size = 12, face = "bold.italic")) +  scale_y_continuous(expand = c(0, 0)) + 
      ggtitle(title_rme_plot1) + xlab("") + ylab("")
    
    fig2 <- ggplot() + 
      geom_sf( color = "black", linewidth = 0.25, data  = mun_cv_df, aes(fill = fill_rme2)) + 
      scale_fill_manual("", values = CScale_dif(5), na.value="red", 
                        guide = guide_legend(direction = "horizontal", nrow = 3, keywidth = 1)) +
      theme(axis.text = element_blank(), axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank(), 
            legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=18), 
            strip.text.x = element_text(size = 12, face = "bold.italic")) +  scale_y_continuous(expand = c(0, 0)) + 
      ggtitle(title_rme_plot2) + xlab("") + ylab("")
    
    fig3 <- ggplot() + 
      geom_sf( color = "black", linewidth = 0.25, data  = mun_cv_df, aes(fill = fill_rme3)) + 
      scale_fill_manual("", values = CScale_dif(5), na.value="red", 
                        guide = guide_legend(direction = "horizontal", nrow = 3, keywidth = 1)) +
      theme(axis.text = element_blank(), axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank(), 
            legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=18), 
            strip.text.x = element_text(size = 12, face = "bold.italic")) +  scale_y_continuous(expand = c(0, 0)) + 
      ggtitle(title_rme_plot3) + xlab("") + ylab("")
    
    fig4 <- ggplot() + 
      geom_sf( color = "black", linewidth = 0.25, data  = mun_cv_df, aes(fill = fill_rme4)) + 
      scale_fill_manual("", values = CScale_dif(5), na.value="red", 
                        guide = guide_legend(direction = "horizontal", nrow = 3, keywidth = 1)) +
      theme(axis.text = element_blank(), axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank(), 
            legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=18), 
            strip.text.x = element_text(size = 12, face = "bold.italic")) +  scale_y_continuous(expand = c(0, 0)) + 
      ggtitle(title_rme_plot4) + xlab("") + ylab("")
    
    grid.arrange(fig1, fig2, fig3, fig4, ncol=4)
    
    dev.print(device = png, filename = "fig6.png", width = 1400, height = 700)
  
  