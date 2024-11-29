##----------------------------------------------------------------------------------------------- ##
#                   DIFFERENCES IN BASELINE CHARACTERISTICS BETWEEN TRAJECTORIES
##----------------------------------------------------------------------------------------------- ##
#
# Description:  This script reads a csv file with psychological, biochemical and genetic data,
#               and investigates the differences in baseline characteristics between
#               different trajectories of depressive symptoms (resilient (65%), 
#               intermediate-stable (20%), symptomatic-chronic (9%), and 
#               late-onset-increasing (6%)) using a nonparametric multivariate analysis from the
#               npvm package. Post hoc tests are performed with univariate Kruskal-Wallis tests
#               for which p values are corrected with Bonferroni. For significant KW-tests, 
#               pairwise comparison with Dunn's test was performed in which p values were also
#               corrected for multiple comparisons with the Bonferroni correction.
#               
#               Boxplots are created to visualize the differences between trajectories and
#               effect sizes (Cliff's Delta) of significant pairwise comparisons are plotted.
#
# Authors:      Plas
# Date:         August 2024
# Version:      1.0
# R.version:    4.2.2 (2022-10-31)
# Rstudio:      2023.12.1+402
#
## ---------------------------------------------------------------------------------------------- ##


# ------------------------------------------------------------------------------------------------ #
#                                      Settings & Dependencies
# ------------------------------------------------------------------------------------------------ #

# numbers of external volumes
home <- "home-4"
heronderzoek <- "heronderzoek-5"

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

# statistics
library(mvnormtest)
library(npmv)
library(rstatix)
library(rcompanion)

# visualizations
library(ggplot2)
library(ggpattern)
library(ggbreak)
library(gridExtra)
library(ggpubr)
library(egg)
library(svglite)


# ------------------------------------------------------------------------------------------------ #
#                                          Data Collection
# ------------------------------------------------------------------------------------------------ #

df_total <- read_csv(paste(save_loc, "a_PRISMO_overall/Data/df_total_A.csv", sep = ""), show_col_types = FALSE)


# how many people per trajectory? (132 were removed because they missed >80% of the baseline data)
table(df_total$outcome)
prop.table(table(df_total$outcome)) * 100

# rename outcome column and remove demographics
df_total <- df_total %>%
  dplyr::rename(trajectory = outcome) %>%
  dplyr::select(-c("gender", "age", "age_cat", "rank_cat", "education_cat", "work_function", 
         "yr_deployment_cat", "Prev_deployment_dummy")) %>% # remove demographics (already investigated in first paper)
  dplyr::select(-c("A_CIS", "A_ETI", "A_SCL", "A_SCL_dep", "A_SRIP")) # remove all questionnaire total scores

df_total$trajectory <- as.factor(df_total$trajectory)


# get variable names instead of R column names for readability in paper
df_varNames <- read_csv(paste(save_loc, "c_Trajectories_Depression/Data/variableName_vs_RcodeName.csv", sep = ""), show_col_types = FALSE)



# ------------------------------------------------------------------------------------------------ #
#                                  Data Exploration And Preparation
# ------------------------------------------------------------------------------------------------ #

# summary(df_total)
# 
# 
# # Visualizations
# # --------------------------------------------------- #
# # histogram
# ggplot(df_total[df_total$AVP < 10,], aes(x=AVP)) + 
#   geom_histogram()
# 
# # boxplots
# df_total[df_total$AVP < 55,] %>%
#   ggplot( aes(x=trajectory, y=AVP, fill=trajectory)) +
#   geom_boxplot() +
#   theme(legend.position="none",
#         plot.title = element_text(size=11)) +
#   ggtitle("Basic boxplot") +
#   xlab("")
# 
# 
# 
# 
# # Assumptions
# # --------------------------------------------------- #
# # variance-covariance matrices
# var_cov_matrix <- by(df_total %>% select(-trajectory), df_total$trajectory, cov)
# # variance and covariances are largest in trajectory 3 and 4 (if var highest in smaller group,
# # ANOVA is liberal)
# 
# 
# # multivariate normality
# traj1 <- t(df_total[df_total["trajectory"] == 1, 8:36])
# traj2 <- t(df_total[df_total["trajectory"] == 2, 8:36])
# traj3 <- t(df_total[df_total["trajectory"] == 3, 8:36])
# traj4 <- t(df_total[df_total["trajectory"] == 4, 8:36])
# 
# mshapiro.test(traj1)
# mshapiro.test(traj2)
# mshapiro.test(traj3)
# mshapiro.test(traj4)
# # there is no multivariate normality





# ------------------------------------------------------------------------------------------------ #
#                               Statistical analysis: multivariate test
# ------------------------------------------------------------------------------------------------ #

# prepare formula's for psychological, biochemical, and genetic variables separately

#psychological variables
formula_quest <- as.formula(paste("A_CIS_fat_sev | A_CIS_conc | A_CIS_mot | A_CIS_act |", 
                                  "A_ETI_gen | A_ETI_phys | A_ETI_emo | A_ETI_sex |",
                                  "A_SCL_ago | A_SCL_ang | A_SCL_som | A_SCL_ins |",
                                  "A_SCL_sen | A_SCL_hos | A_SCL_sla |",
                                  "A_UBOS_EE | A_UBOS_DP | A_UBOS_C |",
                                  "A_TCI_NS | A_TCI_HA | A_TCI_RD | A_TCI_P |",
                                  "A_TCI_SD | A_TCI_C | A_TCI_ST |",
                                  "A_SRIP_re_ex | A_SRIP_avo | A_SRIP_hyp ~ trajectory"))
print(formula_quest)

#biochemical variables
formula_bio <- as.formula("NPY | AVP | OX | TES | SHBG | DHEA | GABA ~ trajectory")
print(formula_bio)

#genetic variables
formula_PRS <- as.formula("PRS_dep | PRS_PTSD ~ trajectory")
print(formula_PRS)



# Multivariate Test
# --------------------------------------------------- #
res_psych <- nonpartest(formula_quest, data = df_total, permreps = 1000) #imputed data
res_bio <- nonpartest(formula_bio, data = df_total, permreps = 1000) #imputed data
res_PRS <- nonpartest(formula_PRS, data = df_total, permreps = 1000) #not imputed

write.csv(res_psych[["results"]], paste(save_loc, "c_Trajectories_Depression/Results/compTraj_multiTest_psych.csv", sep = ""), row.names = FALSE)
#results for psychological variables are significant so continue with post hoc tests

write.csv(res_bio[["results"]], paste(save_loc, "c_Trajectories_Depression/Results/compTraj_multiTest_bio.csv", sep = ""), row.names = FALSE)
write.csv(res_PRS[["results"]], paste(save_loc, "c_Trajectories_Depression/Results/compTraj_multiTest_PRS.csv", sep = ""), row.names = FALSE)
#results for biochemicals and PRS are not significant so do not continue with post hoc tests



# ------------------------------------------------------------------------------------------------ #
#                    Statistical analysis: Post-Hoc tests for questionnaires
# ------------------------------------------------------------------------------------------------ #


# Preparation
# --------------------------------------------------- #
df_total <- df_total %>%
  rename_with(~sub("^A_", "", .), -trajectory)


# Post hoc: univariate kruskal-wallis
# --------------------------------------------------- #

df_p_values_kw <- data.frame(cat = character(), 
                          var = character(), 
                          kw = numeric(),
                          df = numeric(),
                          p_value = numeric(),
                          stringsAsFactors = FALSE)

categories <- c("SCL", "SRIP", "CIS", "UBOS", "TCI", "ETI")
for (cat in categories){
  vars <- names(df_total)[grep(cat, names(df_total))]

  
  plots_list <- list()
  for (var in vars){
    
    res_kw <- kruskal.test(df_total[[var]] ~ trajectory, data = df_total)
    
    print(paste("Processing:", var, "p_value:", res_kw[["p.value"]]))
    
    df_p_values_kw <- rbind(df_p_values_kw, data.frame(cat = cat, var = var, 
                                                       kw = res_kw$statistic,
                                                       df = res_kw$parameter,
                                                       p_value = res_kw$p.value))
  } 
}  

# correct for multiple comparisons
df_p_values_kw$adjusted_p <- p.adjust(df_p_values_kw$p_value, method = "bonferroni")


# customize and save as table for paper
p_value_table_kw <- merge(df_varNames, df_p_values_kw, 
                        by.x = "r_name", by.y = "var", all.y = TRUE, sort = FALSE)
p_value_table_kw[, c("kw", "p_value", "adjusted_p")] <- 
  round(p_value_table_kw[,c("kw", "p_value", "adjusted_p")], 2)
p_value_table_kw <- p_value_table_kw[, c("variable", "kw", "df", "p_value", "adjusted_p")] 

p_value_table_kw$p_value <- ifelse(p_value_table_kw$p_value == 0.00, "<0.01", p_value_table_kw$p_value)
p_value_table_kw$adjusted_p <- ifelse(p_value_table_kw$adjusted_p == 0.00, "<0.01", p_value_table_kw$adjusted_p)
p_value_table_kw$p_value <- gsub("^(?!<)", "  ", p_value_table_kw$p_value, perl = TRUE)
p_value_table_kw$adjusted_p <- gsub("^(?!<)", "  ", p_value_table_kw$adjusted_p, perl = TRUE)

write_xlsx(p_value_table_kw, paste(save_loc, "c_Trajectories_Depression/Results/p_value_table_kw.xlsx", sep = "")) 





# Save means per variable, per trajectory
# --------------------------------------------------- #

summary_df <- df_total
colnames(summary_df) <- gsub("_", "", colnames(summary_df))

summary_df <- summary_df %>%
  group_by(trajectory) %>%
  summarise(across(everything(), list(median = median, min = min, max = max), .names = "{.col}_{.fn}")) %>%
  pivot_longer(cols = -trajectory, names_to = "variable", values_to = "value") %>%
  separate(variable, into = c("variable", "stat"), sep = "_", extra = "merge") %>%
  pivot_wider(names_from = "stat", values_from = "value") %>%
  select(variable, trajectory, median, min, max)

order_index <- match(summary_df$variable, unique(summary_df$variable))
summary_df <- summary_df[order(order_index), ]

trajectory_labels <- c("Resilient", "Intermediate-stable", "Symptomatic-chronic", "Late-onset-increasing")
summary_df$trajectory <- factor(trajectory_labels[as.numeric(summary_df$trajectory)], levels = trajectory_labels)

summary_df[, c("median", "min", "max")] <- 
  round(summary_df[,c("median", "min", "max")], 2)

write_xlsx(summary_df, paste(save_loc, "c_Trajectories_Depression/Results/descriptive_statistics_trajectory_differences.xlsx", sep = ""))









# post hoc: paiwise comparison with Dunn's test
# --------------------------------------------------- #


df_p_values_dunn <- data.frame(var = character(), 
                             trajA = character(),
                             trajB = character(),
                             dunn = numeric(),
                             p_value = numeric(),
                             adjusted_p = numeric(),
                             cd = numeric(),
                             ci_lower = numeric(),
                             ci_upper = numeric(),
                             stringsAsFactors = FALSE)

for (i in 1:nrow(df_p_values_kw)){
  
  p_value <- df_p_values_kw[i, "adjusted_p"]
  cat <- df_p_values_kw[i, "cat"]
  var <- df_p_values_kw[i, "var"]
  
  print(paste("Processing:", var, "p_value:", p_value))
  
  if (p_value < .05){
    
    # Dunn's test
    res_pw_dunn <-  dunn_test(as.formula(paste(var, "~ trajectory")),
                              data = df_total, p.adjust.method = "bonferroni")
    
    # add results to df_p_values_dunn
    df <- data.frame(var = character(6), 
                     trajA = character(6),
                     trajB = character(6),
                     dunn = numeric(6),
                     p_value = numeric(6),
                     adjusted_p = numeric(6),
                     cd = numeric(6),
                     ci_lower = numeric(6),
                     ci_upper = numeric(6),
                     stringsAsFactors = FALSE)
    
    df$var <- res_pw_dunn$.y.
    df$trajA <- res_pw_dunn$group1
    df$trajB <- res_pw_dunn$group2
    df$dunn <- res_pw_dunn$statistic
    df$p_value <- res_pw_dunn$p
    df$adjusted_p <- res_pw_dunn$p.adj
    
    # add effect size (Vargha and Delaney's A)
    for (j  in 1:nrow(df)){
      trajA <- as.numeric(df[j, "trajA"][[1]])
      trajB <- as.numeric(df[j, "trajB"][[1]])
      
      # calculate Cliff's Delta
      res_pw_cd <- cliffDelta(as.formula(paste(var, "~ trajectory")),
                        data = df_total[df_total$trajectory == trajA | df_total$trajectory == trajB,], 
                        ci = TRUE)
      df[df$trajA == trajA & df$trajB == trajB, "cd"] <- res_pw_cd$Cliff.delta
      df[df$trajA == trajA & df$trajB == trajB, "ci_lower"] <- res_pw_cd$lower.ci
      df[df$trajA == trajA & df$trajB == trajB, "ci_upper"] <- res_pw_cd$upper.ci
  
    }
    df_p_values_dunn <- rbind(df_p_values_dunn, df)
  }
}
rm(trajA); rm(trajB)


df_p_values_dunn$significance <- ifelse(df_p_values_dunn$adjusted_p < 0.01, "**",
                                      ifelse(df_p_values_dunn$adjusted_p < 0.05, "*", NA))

write.csv(df_p_values_dunn, paste(save_loc, "c_Trajectories_Depression/Results/df_p_values_dunn.csv", sep = ""), row.names = FALSE)


# customize and save as table for paper
df_p_values_dunn <- merge(df_varNames, df_p_values_dunn, 
                          by.x = "r_name", by.y = "var", all.x = TRUE, sort = FALSE)
df_p_values_dunn <- df_p_values_dunn[!is.na(df_p_values_dunn$trajA), ]
df_p_values_dunn <- df_p_values_dunn %>%
  rename(var = r_name)

p_value_table_dunn <- df_p_values_dunn %>%
  mutate(comparison = paste(trajA, " vs. ", trajB, sep = ""))
p_value_table_dunn <- p_value_table_dunn[, c("variable", "comparison", "dunn", "p_value", "adjusted_p", "significance")]
p_value_table_dunn[, c("dunn", "p_value", "adjusted_p")] <- 
  round(p_value_table_dunn[,c("dunn", "p_value", "adjusted_p")], 2)

p_value_table_dunn$p_value <- ifelse(p_value_table_dunn$p_value == 0.00, "<0.01", p_value_table_dunn$p_value)
p_value_table_dunn$adjusted_p <- ifelse(p_value_table_dunn$adjusted_p == 0.00, "<0.01", p_value_table_dunn$adjusted_p)
p_value_table_dunn$p_value <- gsub("^(?!<)", "  ", p_value_table_dunn$p_value, perl = TRUE)
p_value_table_dunn$adjusted_p <- gsub("^(?!<)", "  ", p_value_table_dunn$adjusted_p, perl = TRUE)

write_xlsx(p_value_table_dunn, paste(save_loc, "c_Trajectories_Depression/Results/p_value_table_dunn.xlsx", sep = "")) 




# ------------------------------------------------------------------------------------------------ #
#                                 Visualize Results: boxplots
# ------------------------------------------------------------------------------------------------ #
  


# Preparation
# --------------------------------------------------- #

# save ranges per subscale
scl_range <- data.frame(lower = c(7, 10, 12, 9, 18, 6, 3), upper = c(35, 50, 60, 45, 90, 30, 15))
rownames(scl_range) <- c("SCL_ago", "SCL_ang", "SCL_som", "SCL_ins", "SCL_sen", "SCL_hos", "SCL_sla")
srip_range <- data.frame(lower = c(6, 9, 7), upper = c(24, 36, 63))
rownames(srip_range) <- c("SRIP_re_ex", "SRIP_avo", "SRIP_hyp")
cis_range <- data.frame(lower = c(0, 0, 0, 0), upper = c(56, 35, 28, 21))
rownames(cis_range) <- c("CIS_fat_sev", "CIS_conc", "CIS_mot", "CIS_act")
ubos_range <- data.frame(lower = c(0, 0, 0), upper = c(6, 6, 6))
rownames(ubos_range) <- c("UBOS_EE", "UBOS_DP", "UBOS_C")
tci_range <- data.frame(lower = c(0, 0, 0, 0, 0, 0, 0), upper = c(15, 15, 15, 15, 15, 15, 15))
rownames(tci_range) <- c("TCI_NS", "TCI_HA", "TCI_RD", "TCI_P", "TCI_SD", "TCI_C", "TCI_ST")
eti_range <- data.frame(lower = c(0, 0, 0, 0), upper = c(11, 5, 5, 6))
rownames(eti_range) <- c("ETI_gen", "ETI_phys", "ETI_emo", "ETI_sex")

range_list <- list(SCL = scl_range, SRIP = srip_range, CIS = cis_range, UBOS = ubos_range,
                   TCI = tci_range, ETI = eti_range)

    
# Add trajectory labels
trajectory_labels <- c("R", "IS", "SC", "LOI")
df_total$trajectory_label <- factor(trajectory_labels[as.numeric(df_total$trajectory)], levels = trajectory_labels)



# Create plot
# --------------------------------------------------- #

plots_list <- list()
for (var in df_p_values_dunn$var){
  cat <- df_p_values_kw[df_p_values_kw$var == var, "cat"]
  
  # get data
  df_fig <- df_total[, c(var, "trajectory_label")]
  colnames(df_fig)[1] <- "score"
  
  # plot
  fig <- df_fig %>%
    ggplot( aes(x=trajectory_label, y=score, fill=trajectory_label)) +
    geom_boxplot(width=0.95, outlier.size = 0.1) 
  
  # customize y-axis
  if (cat %in% names(range_list)){
    fig <- fig + 
      ylim(c(range_list[[cat]][[var, "lower"]], 
             (range_list[[cat]][[var, "upper"]]  + ((range_list[[cat]][[var, "upper"]]/20)) * 3)))
  }
  
  # customize color and theme
  title <- df_p_values_dunn[df_p_values_dunn$var == var, "variable"][1]
  wrapped_title <- strwrap(title, width = 20)
  wrapped_title <- paste(wrapped_title, collapse = "\n")
  
  fig <- fig + 
    scale_fill_manual(values = c("#7FC524", "#e1700e", "#2494c5", "#0E4668"), 
                      guide = "none") +
    ggtitle(paste(" \n \n", wrapped_title, sep = "")) +
    theme_classic() +
    theme(legend.position = "none", 
          plot.title = element_text(size = 12, hjust = 0.5),
          axis.title = element_blank(),
          axis.line.x = element_line(size = 1), axis.line.y = element_line(size = 2),
          axis.text.x = element_text(size = 8), 
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 8),
          axis.line.y.right = element_blank(),
          axis.ticks.y.right = element_blank(), 
          axis.text.y.right = element_blank())
                            
  
  # add significance
  p_values <- df_p_values_dunn[df_p_values_dunn$var == var,]
  p_values$y_values <- c(range_list[[cat]][[var, "upper"]] - (range_list[[cat]][[var, "upper"]]/10),
                         range_list[[cat]][[var, "upper"]] - (range_list[[cat]][[var, "upper"]]/20),
                         range_list[[cat]][[var, "upper"]],
                         range_list[[cat]][[var, "upper"]] - ((range_list[[cat]][[var, "upper"]]/20)) * 3,
                         range_list[[cat]][[var, "upper"]] - (range_list[[cat]][[var, "upper"]]/10),
                         range_list[[cat]][[var, "upper"]] - (range_list[[cat]][[var, "upper"]]/20))
  p_values$y_values <- c(range_list[[cat]][[var, "upper"]] + (range_list[[cat]][[var, "upper"]]/30),
                         range_list[[cat]][[var, "upper"]] + (range_list[[cat]][[var, "upper"]]/10),
                         range_list[[cat]][[var, "upper"]] + ((range_list[[cat]][[var, "upper"]]/20)) * 3,
                         range_list[[cat]][[var, "upper"]],
                         range_list[[cat]][[var, "upper"]] + (range_list[[cat]][[var, "upper"]]/20),
                         range_list[[cat]][[var, "upper"]] + (range_list[[cat]][[var, "upper"]]/10))
  p_values <- p_values[complete.cases(p_values), ]
  
  if (nrow(p_values) > 0){
    fig <- fig +
      ylab(paste("Scale range: ", range_list[[cat]][[var, "lower"]], "-", range_list[[cat]][[var, "upper"]], sep = "")) +
      ggpubr::geom_signif(xmin = as.numeric(p_values$trajA) + 0.05,
                          xmax = as.numeric(p_values$trajB) - 0.05,
                          y_position = p_values$y_values, vjust = 0.5,
                          annotation = p_values$significance)
  }
  
  plots_list[[var]] <- fig
}

# re-order plots for paper
plots_list <- plots_list[c("SCL_ago", "SCL_ang", "SCL_som", "SCL_ins", "SCL_sen", "SCL_hos", "SCL_sla", 
                           "SRIP_re_ex","SRIP_avo", "SRIP_hyp","CIS_fat_sev", "CIS_conc", "CIS_mot", "CIS_act",
                           "UBOS_EE", "UBOS_DP", "UBOS_C", "TCI_HA", "TCI_SD", "TCI_ST", "ETI_emo")]
final_plot <- egg::ggarrange(plots = plots_list, ncol = 7)
file_name <- paste(save_loc, "c_Trajectories_Depression/Figures/sig_boxplot_grid.svg", sep="")
ggsave(file_name, final_plot, device = "svg", width = 360, height = 260, units = "mm")





# ------------------------------------------------------------------------------------------------ #
#                                 Visualize Results: effect sizes
# ------------------------------------------------------------------------------------------------ #



# Preparation
# --------------------------------------------------- #
# read data
df_fig <- read_csv(paste(save_loc, "c_Trajectories_Depression/Results/df_p_values_dunn.csv", sep = ""), show_col_types = FALSE)
df_fig <- merge(df_varNames, df_fig, 
                by.x = "r_name", by.y = "var", all.x = TRUE, sort = FALSE)
df_fig <- df_fig %>%
  rename(var = r_name)

# remove all non significant comparisons
df_fig <- df_fig[complete.cases(df_fig$significance), ]

# prepare data
df_fig <- df_fig %>%
  mutate(combined_var_traj = paste(variable, ": ", trajA, " vs. ", trajB, sep = ""))
df_fig$combined_var_traj <- factor(df_fig$combined_var_traj, levels = df_fig$combined_var_traj)

df_fig[, c("dunn", "p_value", "adjusted_p")] <- 
  round(df_fig[,c("dunn", "p_value", "adjusted_p")], 2)




# plot 1 (middle): effect sizes figure
# --------------------------------------------------- #
p <- df_fig %>%
  ggplot(aes(y = fct_rev(combined_var_traj))) + 
  xlim(c(-1, 1)) +
  theme_classic()

p <- p +
  geom_linerange(aes(xmin=ci_lower, xmax=ci_upper)) +
  geom_point(aes(x=cd), shape=16, size=3, col = "#e1700e")

# add grey/white background
p <- p + annotate("rect",
                  xmin = rep(-Inf, times = 12),
                  xmax = rep(Inf, times = 12),
                  ymin = c(4.5, 8.5, 13.5, 19.5, 23.5, 28.5, 36.5, 44.5, 52.5, 61.5, 71.5, 80.5),
                  ymax = c(6.5, 12.5, 17.5, 20.5, 27.5, 32.5, 40.5, 48.5, 56.5, 66.5, 75.5, 84.5),
                  fill = "grey80", alpha = 0.5)

# add lines
p <- p +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_vline(xintercept = -0.43, linetype="dashed", colour = "#e1700e") +
  geom_vline(xintercept = 0.43, linetype="dashed", colour = "#e1700e") +
  labs(x="Cliff's Delta", y="")

# add data
p <- p +
  geom_linerange(aes(xmin=ci_lower, xmax=ci_upper)) +
  geom_point(aes(x=cd), shape=16, size=3, col = "#e1700e")

# add orange rectangle
p <- p + annotate("rect",
                  xmin = -0.43,
                  xmax = 0.43,
                  ymin = 0,
                  ymax = nrow(df_fig)+1,
                  fill = "#e1700e", alpha = 0.5)

# customize
p <- p + 
  coord_cartesian(ylim=c(1,nrow(df_fig)+1))
  theme(axis.text.y = element_text(size = 8))
p

# save file
file_name <- paste(save_loc, "c_Trajectories_Depression/Figures/eff_sizes.svg", sep="")
ggsave(file_name, p, device = "svg", width = 250, height = 350, units = "mm")

# remove axis 
p_mid <- p + 
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank())
p_mid


# plot 2 (left): effect sizes text
# --------------------------------------------------- #

# prepare
# combine effect size and ci intervals
df_fig2 <- df_fig %>%
  mutate(cd_ci = paste(sprintf('%.2f', cd), " (", sprintf('%.2f', ci_lower), "-", sprintf('%.2f', ci_upper), ")", sep = ""))
df_fig2$cd_ci <- gsub("^(?!-)", " ", df_fig2$cd_ci, perl = TRUE)

# change to character
df_fig2 <- df_fig2 %>%
  mutate(across(c("adjusted_p", "ci_lower", "ci_upper"), ~sprintf("%.2f", .)))

top_row <- data.frame(
  var = "",
  variable = "",
  trajA = "",
  trajB = "",
  dunn = "",
  p_value = "",
  adjusted_p = "p-value",
  cd = "",
  ci_lower = "",
  ci_upper = "",
  significance = "",
  combined_var_traj = "Variable", 
  cd_ci = "Cliff's Delta (95% CI)")
df_fig2 <- rbind(top_row, df_fig2)

df_fig2$combined_var_traj <- factor(df_fig2$combined_var_traj, levels = df_fig2$combined_var_traj)


# create p_left
p_left <- df_fig2 %>%
  ggplot(aes(y = fct_rev(combined_var_traj))) +
  theme_classic()

# add text
p_left <- p_left +
  geom_text(aes(x = 0, label = combined_var_traj), 
            hjust = 0, size = 4, fontface = "bold")

p_left <- p_left +
  geom_text(aes(x = 1, label = cd_ci),
            hjust = 0, size = 4,
            fontface = ifelse(df_fig2$cd_ci == "Cliff's Delta (95% CI)", "bold", "plain"))

# customize
p_left <- p_left +
  theme_void() +
  coord_cartesian(xlim = c(0, 1.4))

p_left


# plot 3 (right): p-values text
# --------------------------------------------------- #

df_fig2$adjusted_p <- ifelse(df_fig2$adjusted_p == "0.00", "<0.01", df_fig2$adjusted_p)
df_fig2$adjusted_p <- gsub("^(?!<)", "  ", df_fig2$adjusted_p, perl = TRUE)
df_fig2$adjusted_p <- gsub("^\\s+(p-value)", "\\1", df_fig2$adjusted_p)

p_right <- df_fig2 %>%
  ggplot() +
  geom_text(aes(x = 0, y = fct_rev(combined_var_traj), label = adjusted_p),
            hjust = 0, size = 4,
            fontface = ifelse(df_fig2$adjusted_p == "p-value", "bold", "plain")) +
  theme_void() 

p_right


# Combine and save
# --------------------------------------------------- #
final_plot <- ggarrange(p_left, p_mid, p_right, nrow = 1, widths = c(1, 0.6, 0.2))

file_name <- paste(save_loc, "c_Trajectories_Depression/Figures/eff_sizes_plus.svg", sep="")
ggsave(file_name, final_plot, device = "svg", width = 250, height = 350, units = "mm")



