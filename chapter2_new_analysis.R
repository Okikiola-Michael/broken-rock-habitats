setwd("C:/Users/alegb/Documents/WSU/Phd/Research/chapter_2/chapter_2_data_analysis/chapter_2_obj2")

#============================
# Packages
# ============================
# install.packages(c("dplyr","tidyr","purrr","ggplot2","lme4","mgcv","broom","yardstick","ranger","xgboost","glmnet","kernlab"))
# Optional:
# install.packages(c("randomForest","merf","catboost"))

library(dplyr)
library(tidymodels)
library(Boruta)
library(caret)
library(ranger)
library(spdep)
library(sf)
library(SAEforest)
library(maps)
# #library(xgboost)
# library(glmnet)
# library(kernlab)
# library(readr)         
# library(randomForest)  
# library(pdp)           
# library(ggcorrplot)
# library(tidyr)
# library(purrr)
# library(ggplot2)
# library(lme4)
# library(mgcv)
# library(broom)
# library(yardstick)

# ============================
# 
# ============================
# Compute using Arjan's function
source("C:/Users/alegb/Documents/WSU/Phd/programming/r_functions/model.stats2_arjan.R")

rawdata <- read.csv("2024_2025_chap2_obj2_data_rf.csv")
names(rawdata)
View(rawdata)

#rawdata <- rawdata[-c(146, 44), ]

#shapefile
#df <- st_read("2025_chapter2_data_all_sites_rf.shp")


# ============================
# 2) 80/20 split by Random sampling
# ============================

set.seed(50)
response_col <- "d50"
id_col <- "grid_id"
train_idx <- createDataPartition(rawdata[[response_col]], p = 0.7, list = FALSE)
train_df  <- rawdata[train_idx, ]
test_df   <- rawdata[-train_idx, ]
View(train_df)

#write.csv(train_df, "training_data_rf_5m.csv")

# Check the number of rows per site
# split_summary <- train_df %>%
#   count(!!sym(id_col)) %>%
#   tidyr::pivot_wider(names_from = plot, values_from = n, values_fill = 0)
# 
# print(split_summary)




#==================================================================================
#feature selection

# data_subset <- rawdata[,c("d84", "B_std", "B_mean", "BR_mean", "BR_std", "CI_mean", "CI_std", "EXG_mean", "EXG_std", "GR_std", "GR_mean","G_std", "G_mean", "NBR_mean", "NBR_std", "NRG_std", "NRG_mean", "R_Std", "NRB_mean", "NRB_std", "RG_std", "RG_mean", "RGB_mean", "RGB_std", "RR_mean", "RR_mean", "VARI_mean", "VARI_std",  "Aspect_std", "Aspect_mean", "CurvPlan_std", "CurvPlan_mean", "CurvPr_std", "CurvPr_mean", "CurvT_std", "CurvT_mean","Dem_std", "Dem_mean","Flow_Dir_std", "Flow_Dir_mean","MinmaxTRI_mean", "MinmaxTRI_std", "SdTRI_mean", "SdTRI_std", "Slope_std", "Slope_mean", "TPI_std", "TPI_mean", "TWI_mean", "TWI_std", "VRM_mean","VRM_std")]

# ============ 
# 12 for circularity, 
# 10 for rock volume,
# 6 for d50, 7 for d84,  
# 8 and 9 for a50, and a84.
# 10:51: predictor variables 
# 13:55 predictors for a84/A84)

train_df <- read.csv("training_data_rf_5m.csv")
test_df <- read.csv("testing_data_rf_5m.csv")

names(train_df)
data_subset <- train_df[, c(8, 13:55)]
View(data_subset)

# ========  Recursive Feature Selection Analysis

ctrl <- rfeControl(functions = rfFuncs, # Using Random Forest functions
                   method = "cv",
                   number = 10,
                   repeats = 5,
                   verbose = FALSE)


rfe_profile <- rfe(x = data_subset[,c(2:43)],
                   y = data_subset$a84,
                   sizes = c(1:ncol(data_subset[,c(2:43)])),
                   rfeControl = ctrl)

print(rfe_profile)
# The selected features
predictors(rfe_profile)
# Visualize performance across feature subsets
plot(rfe_profile, type = c("g", "o"))

###### Data selection based on the best predictors

train_df1 <- subset(train_df,
                    select = c('a50', 
                               predictors(rfe_profile)))
                    
test_df1 <- subset(test_df, 
                   select = c('a50', 
                              predictors(rfe_profile)))


# 
train_df1 <- subset(train_df,
                    select = c('a50', "sd_aspect", "sd_elev", "mean_aspect", "sd_eco_prof" ))
# 
test_df1 <- subset(test_df, select = c('a50', "sd_aspect", "sd_elev", "mean_aspect", "sd_eco_prof" ))

#==================================================================================




# ========  Boruta Analysis

# set.seed(150)
# boruta <- Boruta(d84 ~., data = data_subset, maxRuns = 500)
# print(boruta)
# plot(boruta, las = 2, cex.axis = 0.5)
# 
# #attStats(boruta)
# 
# final_vars <- getSelectedAttributes(boruta, withTentative = FALSE)
# cat("Boruta selected:", paste(final_vars, collapse=", "), "\n")
# 
# train_df1 <- subset(train_df,
#                     select = c('d84', final_vars))
# test_df1 <- subset(test_df, select = c('d84',  final_vars))


# ===========================================================================
#                   Hyperparameter Tuning
# ===========================================================================

# Model and Workflow

# Random forest model with tunable hyperparameters
rf_spec <- rand_forest(
  mode  = "regression",   # continuous outcome
  mtry  = tune(),         # number of predictors sampled at each split
  trees = tune(),         # number of trees
  min_n = tune()          # minimum node size
) %>%
  set_engine("ranger")     # ranger backend (fast + stable)

#Build a workflow for the model

wf <- workflow() %>%
  add_model(rf_spec) %>%
  add_formula(a50 ~ .)

# --------------------------
# Hyperparameter Tuning
# --------------------------

# 5-fold cross-validation on training data
folds <- vfold_cv(train_df1, v = 5)

# Define parameter ranges to explore
param_space <- tune::extract_parameter_set_dials(wf)

param_space <- param_space %>%
  update(
    trees = trees(c(300, 1500)),     # search tree counts from 300 to 1500
    min_n = min_n(c(2, 20))         # search min node sizes from 2 to 20
  ) %>%
  finalize(train_df1)                    # auto-set mtry range from data

# Generate 50 candidate parameter sets
grid <- grid_space_filling (param_space, size = 30)


# Tune the model via CV
tuned <- tune_grid(
  wf,
  resamples = folds,
  grid = grid,
  metrics = metric_set(rsq)
)

# Select best hyperparameters (lowest RMSE)
best <- select_best(tuned, metric = "rsq")

# Insert best parameters back into workflow
final_wf <- finalize_workflow(wf, best)

# --------------------------
# Final Model Fit + Test Evaluation
# --------------------------

final_fit <- fit(final_wf, train_df1)

# Predict on test data
test_pred <- predict(final_fit, test_df1)
test_pred$a50 <- test_df1$a50

metrics_Arjan <- model.stats2(test_pred$a50, test_pred$.pred)
metrics_Arjan

# Compute RMSE + R²
# test_metrics <- metric_set(rmse, rsq, mae)(test_pred, truth = d50, estimate = .pred)
# print(test_metrics)




# ============================
# Plotting the predictor variables for the variables of interest
# ============================

# install.packages(c("ComplexUpset", "dplyr", "tidyr", "ggplot2"))
# 
# library(ComplexUpset)
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# 
# predictors <- list(
#   D50 = c("Aspect (Std)", "Total Curvature (Std)", "Aspect (Mean)", "Profile Curvature (Std)", 
#            "Elevation (Std)", "Elevation (Mean)", "Flow Direction (Std)", "Plan Curvature (Std)", 
#            "Green-Red Ratio (Std)", "Red-RGB Ratio (Mean)"),
#   
#   D84 = c("Total Curvature (Std)", "Profile Curvature (Std)", "Aspect (Mean)", 
#           "Elevation (Mean)", "Plan Curvature (Std)", "Flow Direction (Std)",
#           "Greenness Index (Std)"),
#   
#   A50 = c("Aspect (Std)", "Elevation (Std)","Aspect (Mean)", "Profile Curvature (Std)"),
#   
#   A84 = c("Aspect (Std)", "Elevation (Std)","Aspect (Mean)", "Profile Curvature (Std)"),
#   
#   Volume = c("Profile Curvature (Std)", "Total Curvature (Std)", "Plan Curvature (Std)", 
#              "Elevation (Mean)", "Elevation (Std)", "Blue-Red Ratio (Mean)", "Aspect (Std)",
#              "Blue-Red Ratio (Std)", "Red-Red Ratio (Std)", "Red-Red Ratio (Mean)"),
#   
#   Circularity = c("TPI (Std)", "Plan Curvature (Std)", "Roughness (Std)", "Slope (Mean)", "Elevation (Mean)",
#                   "Total Curvature (Std)","Slope (Std)", "Roughness (Mean)", "Flow Direction (Std)", 
#                   "Green-Red Ratio (Std)", "Aspect (Mean)","Blue-Red Ratio (Std)", "Green-Red Ratio (Mean)")
# )
# 
# 
# long_df <- stack(predictors)
# colnames(long_df) <- c("Predictor", "Variable")
# long_df$present <- 1
# 
# long_df <- long_df |>
#   distinct(Predictor, Variable, .keep_all = TRUE)
# 
# 
# wide_df <- long_df |>
#   complete(Predictor, Variable, fill = list(present = 0)) |>
#   spread(Variable, present)
# 
# View(wide_df)
# ggplot2::theme_set(NULL)
# 
# 
# upset(
#   wide_df,
#   intersect = colnames(wide_df)[-1],
#   base_annotations = list("Intersection size" = intersection_size(counts = TRUE)),
#   min_size = 1,
#   themes= upset_themes
#   # width_ratio = 0.25
# ) +
#   theme(axis.text.x = element_text( "X", angle = 45, hjust = 1),
#     axis.text.y = element_text(size = 10)
#   )
# 
# 


# ============================
# 2) 80/20 split by stratum at ID level
# ============================

# site_col <- "plot"
# # # # Perform 80/20 split per site
# df_split <- rawdata %>%
#   group_by(.data[[site_col]]) %>%
#   mutate(
#     split = if_else(row_number() %in%
#                       sample(row_number(), size = floor(0.7 * n())),
#                     "train", "test")
#   ) %>%
#   ungroup()
# 
# # Separate training and testing data
# train_df <- df_split %>% filter(split == "train")
# test_df  <- df_split %>% filter(split == "test")
# 
# plot(train_df$d50)
# plot(test_df$d50)
# 
# # Check the number of rows per site in each split
# split_summary <- df_split %>%
#   count(!!sym(site_col), split) %>%
#   tidyr::pivot_wider(names_from = split, values_from = n, values_fill = 0)
# 
# print(split_summary)

#removed the outliers

#train_df <- train_df[-c(98, 29), ]


# write.csv(train_df, "all_training_data_jan_2_2026.csv")
# write.csv(test_df, "all_testing_data_jan_2_2026.csv")

# world_map <- map_data("usa")
# plot(world_map)
# 
# ggplot() +
#   geom_map(data = world_map, map = world_map,
#            aes(long, lat, map_id = region),
#            color = "black", fill = "lightgray", size = 0.1) +
#   geom_point(data = rawdata, aes(x = lon, y = lat),
#              color = "red", size = 3, alpha = 0.8) +
#   labs(title = "Geographic Coordinates on a World Map") 
#coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE)

# # Prediction
# use_random_effect <- id_col %in% names(rawdata)   # set TRUE if random effect column exists
# 



# ============================
# Moran's I correlation

# Create spatial neighbors (k-nearest or distance)
# coords <- st_coordinates(df)
# nb <- knn2nb(knearneigh(coords, k = 10))       #nearest neighbors
# lw <- nb2listw(nb, style = "W")                #spatial weights
# 
# # Moran's I test
# moran <- moran.test(df$d84, lw)
# moran
# 
# #===============================
# install.packages("ncf")
# library(ncf)
# 
# coords <- st_coordinates(df)
# correlog <- correlog(x = coords[,1],
#                      y = coords[,2],
#                      z = df$d84,
#                      increment = 10,
#                      resamp = 1000)
# plot(correlog)



# ========================================================


# ggplot() +
#   geom_point(aes(x = test_pred$d50, y = test_pred$.pred), alpha = 0.7) +
#   geom_abline(slope = 1, intercept = 0, linetype = 2) +
#   labs(#title = "Predicted vs Observed (Test Set)", 
#     x = "Observed d50", y = "Predicted d50") +
#   theme_minimal(base_size = 10) + 
#   theme(legend.position = 'top') + theme(legend.title = element_blank()) 
# 


# # 
# write.csv(okwe_train_df, "okwe_train_df.csv")
# write.csv(okwe_test_df, "okwe_test_df.csv")

#====================================================================================
# Correlation Reduction
# correl <- cor(train_df1)
# highcorr <- findCorrelation(correl, cutoff = 0.75)
# #abc <- train_df1[, -highcorr]
# 
# train_df1 <- train_df1[, -highcorr]
# test_df1 <- test_df1[, -highcorr]

#====================================================================================

# -------------------------------
# 2) Define Predictors
# -------------------------------
# train_df1 <- subset(train_df, 
#                     select = c('d50',final_vars))
# test_df1 <- subset(test_df, select = c('d50',  final_vars))

# ============================================================
# A) RANDOM FOREST
# ============================================================
# id_col       <- "site"
# response_col <- "d50"
# pred_cols <- setdiff(names(train_df1), c(response_col, id_col))
# 
# rf_formula <- as.formula(paste(response_col, "~", paste(pred_cols, collapse = " + ")))
# rf_fit <- ranger(
#   formula = rf_formula,
#   data = train_df1,
#   num.trees = 500,
#   importance = "impurity"
#   # mtry = max(1, floor(sqrt(length(pred_cols)))),
#   #seed = 42
# )
# barchart(sort(rf_fit$variable.importance), decreasing = 'TRUE')
# 
# rf_metrics_train <- postResample(rf_fit$predictions, obs = train_df1[[response_col]])
# cat(" Testing RMSE: ", rf_metrics_train["RMSE"], "\n Testing R-squared: ",
#     rf_metrics_train["Rsquared"], "\n")
# 
# 
# obs_pred_train <- tibble(Observed  = train_df1[[response_col]],
#                          Predicted = rf_fit$predictions)
# cor_train <- cor(obs_pred_train)
# cor_train
# 
# # Prediction using test data
# 
# rf_pred <- predict(rf_fit, data = test_df1)$predictions
# # rf_metrics <- metrics_reg(test_df1[[response_col]], rf_pred) %>% mutate(Model = "RandomForest")
# rf_metrics <- postResample(rf_pred, obs = test_df1[[response_col]])
# cat(" Testing RMSE: ", rf_metrics["RMSE"], "\n Testing R-squared: ", rf_metrics["Rsquared"], "\n")
# 
# obs_pred_test <- tibble(Observed  = test_df1[[response_col]],
#                         Predicted = rf_pred)
# cor_testrf <- cor(obs_pred_test)
# cor_testrf
# 
# # Compute using Arjan's function
# source("C:/Users/alegb/Documents/WSU/Phd/programming/r_functions/model.stats2_arjan.R")
# 
# metrics_Arjan <- model.stats2(obs_pred_test$Observed, obs_pred_test$Predicted)
# metrics_Arjan
# 


# ============================================================
# B) Mixed Effect Random Forest
# ============================================================
# Y = as.data.frame(train_df1[, 1])
# X = as.data.frame(train_df1[, c(3:15)])
# test_df1_rep_merf <- as.data.frame(test_df1[, 1])
# test_df1_pred_merf <- as.data.frame(test_df1[, c(1:15)])

# Define the variables i=ot include "site" as random effects

# train_df2 <- subset(train_df, select = c('d50', 'site',  "Dem_mean",  "CurvT_std"))
# 
# test_df2 <- subset(test_df, select = c('d84','site',  "Dem_mean",  "CurvT_std"))
# 
# mrf_fit <- MERFranger(Y = train_df2[, 1], X = train_df2[, c(3, 4)] , random = "(1|site)",  data = train_df2)
# 
# mrf_pred <- predict(mrf_fit, test_df2)

# =======  Metrics by site for merf

# test_df2['d84_pred'] <- mrf_pred
# site_result <- subset(test_df2, select = c("d84", "site", "d84_pred"))
# merf_sitemetrics <- site_result %>%
#   group_by(site) %>%
#   reframe(postResample(d84, obs = d84_pred)["RMSE"])

# ========================================================

# metrics_Arjan2 <- model.stats2(test_df2[, 1], mrf_pred)
# metrics_Arjan2
# 
# 
# mrf_metrics <- postResample(mrf_pred, obs = test_df2[, 1])
# cat(" Testing RMSE: ", mrf_metrics["RMSE"], "\n Testing R-squared: ", mrf_metrics["Rsquared"], "\n")
# 
# obs_pred_testmrf <- tibble(Observed  = test_df2[, 1],
#                         Predicted = mrf_pred)
# cor_testmrf <- cor(obs_pred_testmrf)
# cor_testmrf


# Random Effects

# merf_ref <- as.data.frame(mrf_fit$RandomEffects)
# names(merf_ref)
# 
# ggplot(merf_ref, aes( x = grp, y = condval, fill = grp)) +
#   geom_col(alpha = 0.7) +
#   labs(#title = "Predicted vs Observed (Test Set)",
#      y = "Intercept") +
#   theme(legend.title = element_blank()) + theme(axis.ticks.x =  element_blank()) +
#   theme(axis.text.x  = element_blank()) + theme(axis.title.x =  element_blank()) 
# 

# ============================================================
# C) LINEAR MIXED EFFECTS (random intercept by id)
# ============================================================
# lme_metrics <- NULL
# lme_pred <- NULL
# lme_formula <- as.formula(paste0(response_col, " ~ ", paste(pred_cols, collapse = " + "), " + (1|", id_col, ")"))
# lme_fit <- lmer(lme_formula, data = train_df1, REML = TRUE)
# lme_pred <- predict(lme_fit, newdata = test_df1, allow.new.levels = TRUE)
# lme_metrics <- postResample(lme_pred, obs = test_df1[[response_col]])
# cat(" Testing RMSE: ", lme_metrics["RMSE"], "\n Testing R-squared: ", lme_metrics["Rsquared"], "\n")
# 
# obs_pred_test_lme <- tibble(Observed  = test_df1[[response_col]],
#                         Predicted = lme_pred)
# cor_test_lme <- cor(obs_pred_test_lme)
# cor_test_lme






# ============================================================
# D) Combine and Visualize
# ============================================================
all_metrics <- bind_rows(rf_metrics, lme_metrics) %>% relocate(Model)
print(all_metrics)

metrics_long <- all_metrics %>%
  pivot_longer(cols = c(RMSE, MAE, R2), names_to = "Metric", values_to = "Value")

ggplot(metrics_long, aes(x = Model, y = Value)) +
  geom_col() +
  facet_wrap(~ Metric, scales = "free_y") +
  labs(title = "Model Performance: RF vs XGB vs LME", x = NULL, y = NULL) +
  theme_minimal(base_size = 12)

# Predicted vs Observed plot
pred_df <- bind_rows(
  tibble(Model = "Random Forest", y_true = test_df1[[response_col]], y_pred = rf_pred),
  tibble(Model = "Mixed Effect RF", y_true = test_df1[[response_col]], y_pred = mrf_pred),
  # tibble(Model = "LME (random intercept)", y_true = test_df1[[response_col]], y_pred = lme_pred)
)

ggplot() +
  geom_point(aes(x = test_pred$d84, y = test_pred$.pred), alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs(#title = "Predicted vs Observed (Test Set)", 
    x = "Observed D84", y = "Predicted D84") +
  theme_minimal(base_size = 10) + 
  theme(legend.position = 'top') + theme(legend.title = element_blank()) 


ggplot(rawdata, aes( y = d84, fill = site)) +
  geom_boxplot(alpha = 0.7) +
  # geom_abline(slope = 1, intercept = 0, linetype = 2) +
  # labs(#title = "Predicted vs Observed (Test Set)", 
  #   x = "Observed D84", y = "Predicted D84") +
  theme_minimal(base_size = 12) +
  theme(axis.line.x = element_line()) + theme(axis.text.x = element_blank())


ggplot(rawdata, aes( y = d50, fill = site)) +
  geom_boxplot(alpha = 0.7) +
  # geom_abline(slope = 1, intercept = 0, linetype = 2) +
  # labs(#title = "Predicted vs Observed (Test Set)", 
  #   x = "Observed D84", y = "Predicted D84") +
  theme_minimal(base_size = 12) + theme(axis.text.x = element_blank())


ggplot(rawdata, aes(x = site,  y = d84, fill = site)) +
  geom_boxplot(alpha = 0.7) +
  # geom_abline(slope = 1, intercept = 0, linetype = 2) +
  # labs(#title = "Predicted vs Observed (Test Set)", 
  #   x = "Observed D84", y = "Predicted D84") +
  theme_minimal(base_size = 12) +
  theme(axis.line.x = element_line()) + theme(axis.text.x = element_blank())

# ============================================================
# (Optional) Save
# ============================================================
# write.csv(all_metrics, "model_metrics.csv", row.names = FALSE)
# saveRDS(list(rf = rf_fit, xgb = xgb_fit, lme = if (exists("lme_fit")) lme_fit else NULL),
#         file = "models.rds")




