##----------------------------------------------------------------------------------------------- ##
#                             PREDICTIG DEPRESSION TRAJECTORIES
##----------------------------------------------------------------------------------------------- ##
#
# Description:  This script uses the same csv file as baselineDifferences_trajectories.R and  
#               tries to predict the trajectories with baseline psychological, bio-
#               chemical and genetic data.
#               
#
# Authors:      Plas
# Date:         November 2024
# Version:      1.0
# R.version:    4.2.2 (2022-10-31)
# Rstudio:      2023.12.1+402
#
## ---------------------------------------------------------------------------------------------- ##


# ------------------------------------------------------------------------------------------------ #
#                                      Settings & Dependencies
# ------------------------------------------------------------------------------------------------ #

# Define path of the Rproject to get and save files
save_loc = paste("/Users/aplas2/surfdrive - Plas-2, A. (Xandra)@surfdrive.surf.nl/Documents/PhD/p_PRISMO/")


# Import libraries
# # --------------------------------------------------- #

# read files
library(readxl)
library(writexl)

# data manipulation
library(tidyverse)
library(dplyr)

# statistics/models
library(smotefamily)
library(fastDummies)
library(glmnet)
library(caret)
library(e1071)
library(pROC)
library(MLmetrics)

# Other
library(doMC)

# Functions
#--------------------------------#
source(paste(save_loc, "c_Trajectories_Depression/functions.R", sep=""))




# ------------------------------------------------------------------------------------------------ #
#                                          Data Collection
# ------------------------------------------------------------------------------------------------ #


df_train <- read_csv(paste(save_loc, "a_PRISMO_overall/Data/df_total_A_train.csv", sep = ""), show_col_types = FALSE)
df_test <- read_csv(paste(save_loc, "a_PRISMO_overall/Data/df_total_A_test.csv", sep = ""), show_col_types = FALSE)

# how many people per trajectory? (132 were removed because they missed >80% of the baseline data)
table(df_test$outcome)
prop.table(table(df_test$outcome)) * 100

# rename outcome column and remove demographics
df_train <- df_train %>%
  dplyr::rename(trajectory = outcome) %>%
  dplyr::select(-c("age_cat", "A_CIS", "A_ETI", "A_SCL", "A_SCL_dep", "A_SRIP", "PRS_PTSD", "PRS_dep")) # remove all questionnaire total scores

# rename outcome column and remove demographics
df_test <- df_test %>%
  dplyr::rename(trajectory = outcome) %>%
  dplyr::select(-c("age_cat", "A_CIS", "A_ETI", "A_SCL", "A_SCL_dep", "A_SRIP")) # remove all questionnaire total scores

df_train$trajectory <- as.factor(df_train$trajectory)
df_test$trajectory <- as.factor(df_test$trajectory)




# get variable names instead of R column names for readability in paper
df_varNames <- read_excel(paste(save_loc, "c_Trajectories_Depression/Data/variableName_vs_RcodeName.xlsx", sep = ""))




# ------------------------------------------------------------------------------------------------ #
#                                  Data Preparation
# ------------------------------------------------------------------------------------------------ #

# dummies of categorical variables and scale continuous variables
xy_train <- to_dummy_and_scaled(df_train) #use function from functions.R
xy_test <- to_dummy_and_scaled(df_test) #use function from functions.R

# create balanced training data set
# xy_train_balanced <- SCUT(xy_train, "trajectory")
# xy_train_balanced <- SMOTE(xy_train %>% select(-trajectory), xy_train$trajectory, dup_size = 3)[["data"]]
# xy_train_balanced <- SMOTE(xy_train_balanced %>% select(-class), xy_train_balanced$class, dup_size = 1)[["data"]]
# xy_train_balanced <- xy_train_balanced %>%
#   rename(trajectory = class)

# ------------------------------------------------------------------------------------------------ #
#                                           Create model
# ------------------------------------------------------------------------------------------------ #

x_train <- xy_train %>% select(-trajectory)
y_train <- factor(xy_train$trajectory, levels = c(1, 2, 3, 4), labels = make.names(c("R", "IS", "C", "LOI")))


#Initialize parallel processing
registerDoMC(detectCores()-2)
getDoParWorkers()

set.seed(1234)

# Define parameter grid
param_grid <- expand.grid(alpha = seq(0, 1, by = 0.1), 
                          lambda = seq(0.0001, 0.1, by = 0.001))

control <- trainControl(method = "repeatedcv", 
                        number = 10, 
                        repeats = 10, 
                        search = "grid",
                        classProbs = TRUE,
                        summaryFunction = multiClassSummary,
                        savePredictions = "all") 


enet_model <- train(x = x_train,
                    y = y_train,
                    method = "glmnet",
                    family = "multinomial",
                    type.multinomial = "grouped",
                    metric = "AUC",
                    tuneLength = 50,
                    tuneGrid = param_grid,
                    trControl = control)


saveRDS(enet_model, paste(save_loc, "/c_Trajectories_Depression/Results/enet_model.RDS", sep=""))
enet_model <- readRDS(paste(save_loc, "/c_Trajectories_Depression/Results/enet_model.RDS", sep=""))


enet_y_predicted <- predict(enet_model, x_train)
cm <- confusionMatrix(enet_y_predicted, y_train)
print(cm)

# ------------------------------------------------------------------------------------------------ #
#                                    Predict and save results
# ------------------------------------------------------------------------------------------------ #

x_test <- xy_test %>% select(-trajectory)
y_test <- xy_test$trajectory
y_test <- factor(y_test, levels = c(1, 2, 3, 4), labels = c("R", "IS", "C", "LOI"))


# predictions on test set
enet_y_predicted <- predict(enet_model, x_test)
cm <- confusionMatrix(enet_y_predicted, y_test)
print(cm)


# Get odds ratios for each class
odds_ratios <- calculate_odds_ratios(enet_model$finalModel, enet_model$bestTune$lambda)

# Access each class's data frame
odds_ratios_R <- odds_ratios[["R"]]     # Example: Data frame for class "R"
odds_ratios_IS <- odds_ratios[["IS"]]   # Example: Data frame for class "IS"
odds_ratios_C <- odds_ratios[["C"]]   # Example: Data frame for class "C"
odds_ratios_LOI <- odds_ratios[["LOI"]]   # Example: Data frame for class "LOI"

write.csv(odds_ratios_R, paste(save_loc, "c_Trajectories_Depression/Results/oddsRatios_R.csv", sep = ""), row.names = FALSE)
write.csv(odds_ratios_IS, paste(save_loc, "c_Trajectories_Depression/Results/oddsRatios_IS.csv", sep = ""), row.names = FALSE)
write.csv(odds_ratios_C, paste(save_loc, "c_Trajectories_Depression/Results/oddsRatios_C.csv", sep = ""), row.names = FALSE)
write.csv(odds_ratios_LOI, paste(save_loc, "c_Trajectories_Depression/Results/oddsRatios_LOI.csv", sep = ""), row.names = FALSE)





# ------------------------------------------------------------------------------------------------ #
#                                           Visualizations
# ------------------------------------------------------------------------------------------------ #

# Get AUC/ROC
enet_probs <- predict(enet_model, x_test, type = "prob")

# Initialize list to store ROC curves
roc_curves <- list()

# Calculate ROC curve for each class
for (class in levels(y_test)) {
  binary_labels <- as.numeric(y_test == class)
  roc_curve <- roc(binary_labels, enet_probs[, class])
  roc_curves[[class]] <- roc_curve
}
names(roc_curves) <- c("Resilient", "Intermediate-stable", "Symptomatic-chronic", "Late-onset-increasing")


# Get coefficients
coefs <- coef(enet_model$finalModel, enet_model$bestTune$lambda)
names(coefs) <- c("Resilient", "Intermediate-stable", "Symptomatic-chronic", "Late-onset-increasing")


plot_list <- list()
i = 1
j = 2

# Create figures
for (trajectory in names(coefs)) {
  print(trajectory)
  trajectory_coefs <- as.data.frame(as.matrix(coefs[[trajectory]]))
  
  #---------------------------#
  ## plot ROC curves
  #---------------------------#
  roc_plot <- ggplot(data.table("Specificity" = roc_curves[[trajectory]][["specificities"]], "Sensitivity" = roc_curves[[trajectory]][["sensitivities"]]), aes(Specificity, Sensitivity)) +
    geom_path(color = "#e1700e") +  
    scale_x_reverse() +
    geom_abline(intercept = 1, slope = 1, colour = "grey") +
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
    theme_classic() +
    labs(title = paste0("AUROC = ", signif(roc_curves[[trajectory]][["auc"]], 3), sep = ""))
  
  file_name <- paste(save_loc, "c_Trajectories_Depression/Figures/ROC_", trajectory, ".svg", sep = "")
  ggsave(file_name, roc_plot, device = "svg", width = 360, height = 260, units = "mm")
  
  
  #---------------------------#
  ## plot coefficients
  #---------------------------#
  # make dataframe with correct variable names for coefs
  trajectory_coefs$r_name <- rownames(trajectory_coefs)
  trajectory_coefs$r_name <- gsub("^A_", "", trajectory_coefs$r_name)
  trajectory_coefs$r_name <- ifelse(trajectory_coefs$r_name %in% c("NPY", "AVP", "OX", "TES", "SHBG", "DHEA", "GABA"), 
                                    paste0("BIO_", trajectory_coefs$r_name), 
                                    trajectory_coefs$r_name)
  trajectory_coefs <- merge(trajectory_coefs, df_varNames, by="r_name", all.x = TRUE)
  
  df_trajectory_coefs <- data.frame(r_name = trajectory_coefs$r_name,
                                    variable = trajectory_coefs$variable,
                                    coefficient = trajectory_coefs[, "1"])
  df_trajectory_coefs <- df_trajectory_coefs[df_trajectory_coefs$r_name != "(Intercept)", ] # remove intercept
  
  
  coefs_plot <- ggplot(df_trajectory_coefs, aes(x = reorder(variable, coefficient), y = coefficient)) +
    geom_bar(stat = "identity", fill = "#e1700e", color = "white") +
    coord_flip() +
    labs(x = "Predictors", y = "Coefficient Value") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8)) +
    ggtitle(paste("Trajectory: ", trajectory, " - Coefficients of the best performing model", sep=""))
  
  
  file_name <- paste(save_loc, "c_Trajectories_Depression/Figures/coefs_", trajectory, ".svg", sep = "")
  ggsave(file_name, coefs_plot, device = "svg", width = 360, height = 260, units = "mm")

  
  # TOP 10
  # Prepare data: Select the top 10 coefficients by absolute value
  df_trajectory_coefs <- df_trajectory_coefs %>%
    mutate(abs_coefficient = abs(coefficient)) %>%          # Add a column for absolute values
    arrange(desc(abs_coefficient)) %>%                     # Sort by absolute value
    slice(1:10) %>%                                        # Select the top 10 rows
    arrange(desc(coefficient))                             # Re-sort by actual coefficient value for plotting
  
  # Plot the top 10 coefficients
  coefs_plot_top10 <- ggplot(df_trajectory_coefs, aes(x = reorder(variable, coefficient), y = coefficient)) +
    geom_bar(stat = "identity", fill = "#e1700e", color = "white") +
    coord_flip() +
    labs(x = "Predictors", y = "Coefficient Value") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8)) +
    theme(
      axis.text.y = element_text(size = 8),
      plot.title = element_text(hjust = 0.7, vjust = 1),  # Adjust title position
      plot.margin = margin(t = 20, r = 10, b = 10, l = 10)  # Increase the top margin to avoid overlap
    ) +
    ggtitle(paste("Trajectory: ", trajectory, " - Top 10 coefficients of the best performing model",  sep=""))
  
  
  file_name <- paste(save_loc, "c_Trajectories_Depression/Figures/coefs_top10_", trajectory, ".svg", sep = "")
  ggsave(file_name, coefs_plot, device = "svg", width = 360, height = 260, units = "mm")
  
  
  #---------------------------#
  ## Save Plots to List
  #---------------------------#
  plot_list[[i]] <- coefs_plot_top10  # Add Coefficients plot
  plot_list[[j]] <- roc_plot  # Add ROC plot
  i <- i + 2
  j <- j + 2
  
  
}

final_plot <- egg::ggarrange(plots = plot_list, ncol = 2, widths = c(2,1))
file_name <- paste(save_loc, "c_Trajectories_Depression/Figures/coefs_roc_grid.svg", sep = "")
ggsave(file_name, final_plot, device = "svg", width = 260, height = 360, units = "mm")





























# ------------------------------------------------------------------------------------------------ #
#                         Create model - significant baseline differences
# ------------------------------------------------------------------------------------------------ #

xy_train_balanced_sigBaseline <- xy_train_balanced %>%
  select("A_SCL_ago", "A_SCL_ang", "A_SCL_som", "A_SCL_ins", "A_SCL_sen", "A_SCL_hos", "A_SCL_sla", 
         "A_SRIP_re_ex","A_SRIP_avo", "A_SRIP_hyp","A_CIS_fat_sev", "A_CIS_conc", "A_CIS_mot", "A_CIS_act",
         "A_UBOS_EE", "A_UBOS_DP", "A_UBOS_C", "A_TCI_HA", "A_TCI_SD", "A_TCI_ST", "A_ETI_emo", "trajectory")

x_train_sigBaseline <- xy_train_balanced_sigBaseline %>% select(-trajectory)
y_train_sigBaseline <- factor(xy_train_balanced_sigBaseline$trajectory, levels = c(1, 2, 3, 4), labels = make.names(c("R", "IS", "C", "LOI")))

set.seed(1234)

enet_model_sigBaseline <- train(x = x_train_sigBaseline,
                                y = y_train_sigBaseline,
                                method = "glmnet",
                                family = "multinomial",
                                type.multinomial = "grouped",
                                tuneLength = 50,
                                tuneGrid = param_grid,
                                trControl = control)


saveRDS(enet_model_sigBaseline, paste(save_loc, "/c_Trajectories_Depression/Results/enet_model_sigBaseline.RDS", sep=""))
enet_model_sigBaseline <- readRDS(paste(save_loc, "/c_Trajectories_Depression/Results/enet_model_sigBaseline.RDS", sep=""))



# ------------------------------------------------------------------------------------------------ #
#                                    Predict and save results
# ------------------------------------------------------------------------------------------------ #

x_test_sigBaseline <- xy_test %>% select(-trajectory)
y_test_sigBaseline <- xy_test$trajectory
y_test_sigBaseline <- factor(y_test_sigBaseline, levels = c(1, 2, 3, 4), labels = c("R", "IS", "C", "LOI"))


x_test_sigBaseline <- x_test_sigBaseline %>%
  select("A_SCL_ago", "A_SCL_ang", "A_SCL_som", "A_SCL_ins", "A_SCL_sen", "A_SCL_hos", "A_SCL_sla", 
         "A_SRIP_re_ex","A_SRIP_avo", "A_SRIP_hyp","A_CIS_fat_sev", "A_CIS_conc", "A_CIS_mot", "A_CIS_act",
         "A_UBOS_EE", "A_UBOS_DP", "A_UBOS_C", "A_TCI_HA", "A_TCI_SD", "A_TCI_ST", "A_ETI_emo")


## predictions on test set:
enet_y_predicted_sigBaseline <- predict(enet_model_sigBaseline, x_test_sigBaseline)
cm <- confusionMatrix(enet_y_predicted_sigBaseline, y_test_sigBaseline)
print(cm)




# ------------------------------------------------------------------------------------------------ #
#                                           Visualizations
# ------------------------------------------------------------------------------------------------ #

# Get coefficients
coefs_sigBaseline <- coef(enet_model_sigBaseline$finalModel, enet_model_sigBaseline$bestTune$lambda)
names(coefs_sigBaseline) <- c("Resilient", "Intermediate-stable", "Symptomatic-chronic", "Late-onset-increasing")

for (trajectory in names(coefs_sigBaseline)) {
  print(trajectory)
  trajectory_coefs <- as.data.frame(as.matrix(coefs_sigBaseline[[trajectory]]))
  
  # make dataframe with correct variable names
  trajectory_coefs$r_name <- rownames(trajectory_coefs)
  trajectory_coefs$r_name <- gsub("^A_", "", trajectory_coefs$r_name)
  trajectory_coefs$r_name <- ifelse(trajectory_coefs$r_name %in% c("NPY", "AVP", "OX", "TES", "SHBG", "DHEA", "GABA"), 
                                    paste0("BIO_", trajectory_coefs$r_name), 
                                    trajectory_coefs$r_name)
  trajectory_coefs <- merge(trajectory_coefs, df_varNames, by="r_name", all.x = TRUE)
  
  df_trajectory_coefs <- data.frame(r_name = trajectory_coefs$r_name,
                                    variable = trajectory_coefs$variable,
                                    coefficient = trajectory_coefs[, "1"])
  df_trajectory_coefs <- df_trajectory_coefs[df_trajectory_coefs$r_name != "(Intercept)", ] # remove intercept
  
  
  coefs_plot <- ggplot(df_trajectory_coefs, aes(x = reorder(variable, coefficient), y = coefficient)) +
    geom_bar(stat = "identity", fill = "#e1700e", color = "white") +
    coord_flip() +
    labs(x = "Predictors", y = "Coefficient Value") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8)) +
    ggtitle(paste("Coefficients of the Best Performing Model - Trajectory: ", trajectory, sep=""))
  
  
  file_name <- paste(save_loc, "c_Trajectories_Depression/Figures/coefs_sigBaseline_", trajectory, ".svg", sep = "")
  ggsave(file_name, coefs_plot, device = "svg", width = 360, height = 260, units = "mm")
  
}


