rm(list = ls())
library(tidyverse)
library(readr)
library(janitor)
library(scales)
library(ggpubr)
library(minpack.lm)
options(scipen =999)

# this reads in all the functions needed for fitting
source("code/F_standard_curve_fitting_functions.R")

# STEP 1) ensure data are in the right format and visualise

example_data <- read_csv("data-inputs/dummy_data.csv")

# ensure data has these column names (it can have more columns)
# needs to be in LONG format i.e. cant have a seperate mfi columns for each antigen
head(example_data)
table(example_data$antigen)

# plot the data to make sure it looks correct

ggplot(example_data, aes(y = mfi, x = antigen_conc_iu_m_l , color = as.factor(antigen))) +
  geom_line(linewidth = 0.8) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()


# STEP 2) fit the standard curve model for each parameter

hrp2_params <- fit_params_function(dat = example_data, antigen_name = "hrp2")

hrp2_params

# STEP 3) predict mfi over a sequence of antigen concentration values to
#         check the model fit 

hrp2_standard_curve <- make_model_predictions(data = example_data,
                                              fitted_params = hrp2_params)

p1 = ggplot() +
  geom_line(data = hrp2_standard_curve, aes(x = antigen_conc,
                                            y = predicted_mfi, 
                                            color = "Standard Curve"),
            linewidth = 0.8) +
  geom_point(data =dplyr::filter(example_data, antigen == "hrp2"),
             aes(x = antigen_conc_iu_m_l, y = mfi, 
                 color = "Data")) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual("", values = c("black", "deeppink3")) +
  ggtitle("Fitted Standard Curve for HRP2") +
  labs(x = "Antigen Concentration (IU/ml)", y = "MFI") +
  theme_bw()


# STEP 4) Predict antigen concentrations from a given MFI value

model_params = hrp2_params


preds <- pred_antigen_conc_from_mfi(mfi_values = c(100, 1000),
                           model_params = hrp2_params)


# can double check that these point do indeed fall on the SC 
p1 + geom_point(data = preds, aes(x = antigen_conc_pred,  y = mfi),
                size = 2, color = "dodgerblue1")







