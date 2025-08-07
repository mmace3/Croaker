#-------------------------------------------------------------------------------
# INFO FOR CODE IN THIS FILE
#-------------------------------------------------------------------------------

# Created: 06 August 2025 by MMMIII

# This file contains code for running index standardization models and saving
# output for use in a summary document. The index standardization is part of
# the current ASMFC Atlantic Croaker stock assessment.

#-------------------------------------------------------------------------------
# Load packages, set options
#-------------------------------------------------------------------------------
library(DHARMa)
library(tidyverse)
library(cowplot)
library(glmmTMB)

source("utils.R") # get_data() function

#-------------------------------------------------------------------------------
# Read in Data
#-------------------------------------------------------------------------------


croaker_catch <- read.csv("ChesMMAP_AtlanticCroaker_Catch_2002_2024_new_adjusted.csv",
                          header = TRUE) %>%
  mutate(time = strptime(time, "%Y-%m-%dT%H:%M:%S", tz="UTC")) %>%
  mutate(date = as.Date(time, format = "%Y/%m/%d")) %>%
  mutate(year = format(date, "%Y")) %>%
  mutate(month_int = as.numeric(format(date, "%m"))) %>%
  mutate(month_name = format(date, "%b")) %>%
  mutate(month_name_factor = factor(month_name, 
                                    levels = c("Mar", "Apr", "May", "Jun", "Jul", "Aug", 
                                               "Sep", "Oct", "Nov"), 
                                    ordered = TRUE)) %>%
  mutate(region_dstrat = paste(region, dstrat, sep = "_")) %>%
  mutate(new_dstrat = if_else(depth < 12.2, 1, 2)) %>%
  mutate(new_region_dstrat = paste(region, new_dstrat, sep = "_")) %>%
  filter(!is.na(depth)) %>%
  mutate(new_region = case_when(region %in% c("1", "A", "2") ~ 1,
                                region %in% c("3", "B") ~ 2,
                                region %in% c("4", "C") ~ 3,
                                region %in% c("5", "D") ~ 4)) %>%
  mutate(pos = if_else(count == 0, 0, 1)) %>%
  mutate(CPUE = count/areasw) %>%
  mutate(CPUE_km = count/(areasw/1000))%>%
  mutate(new_station = gsub('[ ]', '', station)) %>%
  mutate(new_region_name_factor = as.factor(new_region)) %>%
  mutate(year_factor = as.factor(year)) %>%
  mutate(areasw_km = areasw/1000)


#-------------------------------------------------------------------------------
# Run models in loop and save output to be used in markdown document
#-------------------------------------------------------------------------------

# create data frame with models

models_glmmTMB <-
  data.frame(model = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 
                       16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29),
             season = c("summer",
                        "summer",
                        "summer",
                        "summer",
                        "summer",
                        "summer",
                        "summer",
                        "summer",
                        "summer",
                        "summer",
                        "summer",
                        "summer",
                        "summer",
                        "summer",
                        "summer",
                        
                        "fall",
                        "fall",
                        "fall",
                        "fall",
                        "fall",
                        "fall",
                        "fall",
                        "fall",
                        "fall",
                        "fall",
                        "fall",
                        "fall",
                        "fall",
                        "fall"),
             season_model = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                              1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
             
             lhs = c("count", 
                     "count",
                     "count",
                     "count",
                     "count",
                     "count",
                     "count",
                     "count",
                     "count",
                     "count",
                     "count",
                     "count",
                     "count",
                     "count",
                     "count",
                     
                     "count",
                     "count",
                     "count",
                     "count",
                     "count",
                     "count",
                     "count",
                     "count",
                     "count",
                     "count",
                     "count",
                     "count",
                     "count",
                     "count"),
             
             mean = c("year_factor + offset(log(areasw))",
                      "year_factor + depth + offset(log(areasw))",
                      "year_factor + DO + offset(log(areasw))",
                      "year_factor + SA + offset(log(areasw))",
                      "year_factor + WT + offset(log(areasw))",
                      "year_factor + depth + DO + SA + WT + offset(log(areasw))",
                      "year_factor + offset(log(areasw))",
                      "year_factor + new_region_name_factor + offset(log(areasw))",
                      "year_factor + new_region_name_factor + offset(log(areasw))",
                      "year_factor + new_region_name_factor + depth + offset(log(areasw))",
                      "year_factor + new_region_name_factor + SA + offset(log(areasw))",
                      "year_factor + new_region_name_factor + DO + offset(log(areasw))",
                      "year_factor + new_region_name_factor + WT + offset(log(areasw))",
                      "year_factor + new_region_name_factor + depth + DO + SA + WT + offset(log(areasw))",
                      "year_factor + new_region_name_factor + month_name_factor + offset(log(areasw))",
                      
                      "year_factor + offset(log(areasw))",
                      "year_factor + new_region_name_factor + offset(log(areasw))",
                      "year_factor + new_region_name_factor + depth + offset(log(areasw))",
                      "year_factor*new_region_name_factor + offset(log(areasw))",
                      "year_factor*new_region_name_factor + depth + offset(log(areasw))",
                      "year_factor + new_region_name_factor + month_name_factor + offset(log(areasw))",
                      "year_factor + new_region_name_factor + month_name_factor + offset(log(areasw))",
                      "year_factor + new_region_name_factor + month_name_factor + depth + offset(log(areasw))",
                      "year_factor*new_region_name_factor + month_name_factor + offset(log(areasw))",
                      "year_factor*new_region_name_factor + month_name_factor + depth + offset(log(areasw))",
                      "year_factor*new_region_name_factor + month_name_factor + offset(log(areasw))",
                      "year_factor*new_region_name_factor + month_name_factor + offset(log(areasw))",
                      "year_factor*new_region_name_factor + month_name_factor + offset(log(areasw))",
                      "year_factor*new_region_name_factor + month_name_factor + depth + offset(log(areasw))"
             ),
             dispformula = c("1",
                             "1",
                             "1",
                             "1",
                             "1",
                             "1",
                             "year_factor",
                             "year_factor",
                             "year_factor + new_region_name_factor",
                             "year_factor + new_region_name_factor",
                             "year_factor + new_region_name_factor",
                             "year_factor + new_region_name_factor",
                             "year_factor + new_region_name_factor",
                             "year_factor + new_region_name_factor",
                             "year_factor + new_region_name_factor",
                             
                             "1",
                             "1",
                             "1",
                             "1",
                             "1",
                             "1",
                             "1",
                             "1",
                             "1",
                             "1",
                             "1",
                             "1",
                             "1",
                             "1"
             ),
             ziformula = c("0",
                           "0",
                           "0",
                           "0",
                           "0",
                           "0",
                           "0",
                           "0",
                           "0",
                           "0",
                           "0",
                           "0",
                           "0",
                           "0",
                           "0",
                           
                           "0",
                           "0",
                           "0",
                           "0",
                           "0",
                           "0",
                           "new_region_name_factor + month_name_factor",
                           "0",
                           "0",
                           "0",
                           "year_factor*new_region_name_factor + month_name_factor",
                           "new_region_name_factor + month_name_factor",
                           "new_region_name_factor*month_name_factor",
                           "new_region_name_factor*month_name_factor"
             ),
             family = c("nbinom2",
                        "nbinom2",
                        "nbinom2",
                        "nbinom2",
                        "nbinom2",
                        "nbinom2",
                        "nbinom2",
                        "nbinom2",
                        "nbinom2",
                        "nbinom2",
                        "nbinom2",
                        "nbinom2",
                        "nbinom2",
                        "nbinom2",
                        "nbinom2",
                        
                        "nbinom2",
                        "nbinom2",
                        "nbinom2",
                        "nbinom2",
                        "nbinom2",
                        "nbinom2",
                        "nbinom2",
                        "nbinom2",
                        "nbinom2",
                        "nbinom2",
                        "nbinom2",
                        "nbinom2",
                        "nbinom2",
                        "nbinom2")
  )

# function used below to run models in loop. Calls glmmTMB.

try_glmmTMB<- function(LHS, 
                       RHS, 
                       disp = "~1", 
                       family = "gaussian()", 
                       zinf = "~0", 
                       data, 
                       iter = NA_real_)
{
  
  
  mean <- as.formula(paste0(LHS, " ~ ", RHS))
  ziformula <- as.formula(paste0("~", zinf))
  dispformula <- as.formula(paste0("~", disp))
  # print(mean)
  temp_model <- glmmTMB(formula = mean,
                        ziformula = ziformula,
                        dispformula = dispformula,
                        family = family,
                        data = data)
  #         data = data)
  # tryCatch(expr = {glmmTMB(formula = mean,
  #                         ziformula = ziformula,
  #                         dispformula = dispformula,
  #                         family = family,
  #                         data = data)},
  #                 
  #          warning = function(w)
  #          {
  #            if(!is.na(iter)) message("There was a warning on iteration ", iter)
  #            message("The was a warning: ", w)
  #            message("The model is: \n", as.formula(paste0(LHS, " ~ ", RHS)))
  #            return(list(1, w))
  #            
  #          },
  #          error = function(e)
  #          {
  #            if(!is.na(iter)) message("There was a error on iteration ", iter)
  #            message("The was an error: \n", e)
  #            message("The model is: \n", paste0(LHS, " ~ ", RHS))
  #            return(list(1, e))
  #          }
  #          
  #          
  # )
  # 
  return(temp_model)
}

# create output directory for results

if(!dir.exists("output/glmmTMB"))
{
  dir.create("output/glmmTMB") 
}


# run models and save output
# shouldn't run into errors during evaluation through loop so no code to
# try to run things and then see if there is an error or not.

for(j in unique(models_glmmTMB$season)) {
  
  if(j == "summer") temp_months <- c("Jul", "Jun")
  if(j == "fall") temp_months <- c("Sep", "Oct", "Nov")
  if(!(j %in% c("summer", "fall"))) temp_months <- -99      
  
  temp_my_data <- subset(croaker_catch, month_name %in% temp_months & !is.na(WT) & !is.na(SA) & !is.na(DO)) %>%
    droplevels()
  
  temp_models_glmmTMB <- subset(models_glmmTMB, season == j)
  
  n_temp <- nrow(temp_models_glmmTMB)
  
  print(paste0("This is season ", j))

  
  
  for(i in c(1:n_temp)) {
    
    print(paste0("This is model ", i))
    print(paste0("The mean for model ", i, " in season ", j, " is ", temp_models_glmmTMB$mean[i]))
    print(paste0("The dispersion for model ", i , " in season ", j, " is ", temp_models_glmmTMB$dispformula[i]))
    print(paste0("The zi for model ", i , " in season ", j, " is ", temp_models_glmmTMB$ziformula[i]))

    
    temp_results <- try_glmmTMB(temp_models_glmmTMB$lhs[i],
                                temp_models_glmmTMB$mean[i],
                                disp = temp_models_glmmTMB$dispformula[i],
                                zinf = temp_models_glmmTMB$ziformula[i],
                                iter = i,
                                family = temp_models_glmmTMB$family[i],
                                data = temp_my_data)

    temp_dir <- paste0("output/glmmTMB/model_", j, "_", i)
    
    if(!dir.exists(temp_dir))
    {
      dir.create(temp_dir) 
    }
    print(paste0("temp_dir for model ", i, " in season ", j, " is ", temp_dir))
    saveRDS(temp_results,
            file = paste0(temp_dir, "/model_", j, "_", i, ".rds"))
    
    print(paste0("Done with season ", j, " model ", i))
  }
}


# Now get model results and create and save graphs for later use in document.

models_glmmTMB_results <- data.frame(season = character(),
                                     season_model = integer(),
                                     converged = integer(),
                                     pdh = integer(),
                                     AIC = double(),
                                     dev_exp = double()
)


for(j in unique(models_glmmTMB$season)) {
  
  if(j == "summer") temp_months <- c("Jul", "Jun")
  if(j == "fall") temp_months <- c("Sep", "Oct", "Nov")
  if(!(j %in% c("summer", "fall"))) temp_months <- -99      
  
  temp_my_data <- subset(croaker_catch, month_name %in% temp_months & !is.na(WT) & !is.na(SA) & !is.na(DO)) %>%
    droplevels()
  
  
  null_model_temp <- glmmTMB(count ~ 1,
                             ziformula = ~0,
                             dispformula = ~1,
                             data = temp_my_data, 
                             family = nbinom2)
  
  null_dev_temp <- summary(null_model_temp)[["logLik"]][1] * -2
  
  
  n_temp <- nrow(subset(models_glmmTMB, season == j))
  print(paste0("This is season ", j))
  
  
  for(i in c(1:n_temp))
  {
    print(paste0("This is model ", i))
    temp_dir <- paste0("output/glmmTMB/model_", j, "_", i)
    
    temp_model <- readRDS(paste0(temp_dir, "/model_", j, "_", i, ".rds"))
    temp_model_number <- i
    
    temp_converged <- ifelse(temp_model[["fit"]][["convergence"]] == 0, 1, 0)
    temp_pdh <- ifelse(temp_model[["sdr"]][["pdHess"]] == TRUE, 1, 0)
    
    temp_row <- which(models_glmmTMB$season_model == i & models_glmmTMB$season == j)
    
    if(temp_converged == TRUE & temp_pdh == TRUE){
      
      temp_AIC <- AIC(temp_model)
      
      
      summary_temp_model <- summary(temp_model)
      res_dev <- -2*summary_temp_model[["logLik"]][1]
      
      temp_dev_exp <- (1 - (res_dev/null_dev_temp))*100
      
      # make plots
      dharma_fit <- simulateResiduals(fittedModel = temp_model)
      
      # qq plot for residuals
      png(filename = paste0(temp_dir, "/model_", j, "_", i, "_qq_plot.png"))
      plotQQunif(dharma_fit)
      dev.off()
      
      # overall model
      png(filename = paste0(temp_dir, "/model_", j, "_", i, "_model_fit.png"))
      plotResiduals(dharma_fit, rank = TRUE, smoothScatter = TRUE,  quantreg=T)
      dev.off()
      
      
      temp_names_list <- strsplit(models_glmmTMB$mean[temp_row], " + ", fixed = TRUE)
      if(grepl("*", temp_names_list, fixed=TRUE)) temp_names_list <- strsplit(temp_names_list[[1]], "*", fixed = TRUE)
      temp_names <- unlist(temp_names_list)
      temp_names <- temp_names[which(temp_names != "offset(log(areasw))")]
      
      # plots for each factor in mean model
      for(k in temp_names)
      {
        temp_variable <- temp_my_data[,k]
        
        png(file = paste0(temp_dir,"/model_", j, "_", i, "_model_fit_", k, ".png"))
        plotResiduals(dharma_fit, 
                      rank = TRUE, 
                      smoothScatter = TRUE, 
                      form = temp_variable,
                      quantreg = TRUE,
                      xlab = k)
        dev.off()
        
        
        
        # rootgram
        fit <- fitted(temp_model)
        obs <- temp_my_data$count
        all <- get_data(fit = fit, obs = obs)
        
        rootgram_plot <-
          ggplot(all, aes(x = mids, y = sqrt_fit)) + 
          geom_rect(aes(ymin = sqrt_fit - sqrt_obs,
                        ymax = sqrt_fit,
                        xmin = low_brks,
                        xmax = high_brks),
                    fill = "grey",
                    alpha = 1) +
          geom_line(aes(group = 1), color = "red") +
          geom_point(color = "red") +
          geom_hline(yintercept = 0, linetype = "dashed") +
          labs(x = "interval mid-points",
               y = "sqrt(Frequency)") +
          scale_x_continuous(breaks = all$mids,
                             labels = all$mids) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
        
        ggsave(rootgram_plot,
               file = paste0(temp_dir,"/model_", j, "_", i, "_model_fit_rootgram.png"),
               width = 4,
               height = 4,
               units = "in")
        
        
      }
      
    } else {
      
      temp_AIC <- NA_real_
      temp_dev_exp <- NA_real_
      
    }

    
    
    models_glmmTMB_results[temp_row, "season_model"] <- temp_model_number
    models_glmmTMB_results[temp_row, "season"] <- j
    models_glmmTMB_results[temp_row, "converged"] <- temp_converged
    models_glmmTMB_results[temp_row, "pdh"] <- temp_pdh
    models_glmmTMB_results[temp_row, "AIC"] <- temp_AIC
    models_glmmTMB_results[temp_row, "dev_exp"] <- temp_dev_exp

    print(paste0("Done with season ", j, " model ", i))
    
  } # model i loop closed
  
} # season j loop closed


# Save some results to table for model comparisons
models_glmmTMB_table <-
  models_glmmTMB_results %>%
  mutate(season_model = as.numeric(season_model)) %>%
  full_join(models_glmmTMB, by = c("season", "season_model"))

write.csv(models_glmmTMB_table,
          file = "output/glmmTMB/glmmTMB_results_table.csv",
          row.names = FALSE)

