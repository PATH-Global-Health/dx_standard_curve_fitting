# dx_standard_curve_fitting

Code to fit a 5-parameter standard curve to Luminex data.

## Setup

Clear the environment and load required libraries:
```r
rm(list = ls())

library(tidyverse)
library(readr)
library(janitor)
library(scales)
library(ggpubr)
library(minpack.lm)

options(scipen = 999)
```

Load the functions needed for fitting:
```r
source("F_standard_curve_fitting_functions.R")
```

## Step 1 — Load data and visualise

Read in data and check format. Data must be in **long format** with a single `mfi` column (not a separate column per antigen) and include at minimum the columns shown below:
```r
example_data <- read_csv("dummy_data.csv")

head(example_data)
table(example_data$antigen)
```

Plot the raw data to check it looks correct:
```r
ggplot(example_data, aes(y = mfi, x = antigen_conc_iu_m_l, color = as.factor(antigen))) +
  geom_line(linewidth = 0.8) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()
```

## Step 2 — Fit the standard curve model

Fit the 5-parameter model for each antigen:
```r
hrp2_params <- fit_params_function(dat = example_data, antigen_name = "hrp2")
hrp2_params
```

## Step 3 — Predict MFI over a concentration sequence

Generate model predictions across a sequence of antigen concentration values to check the fit:
```r
hrp2_standard_curve <- make_model_predictions(data = example_data,
                                              fitted_params = hrp2_params)

p1 <- ggplot() +
  geom_line(data = hrp2_standard_curve, aes(x = antigen_conc,
                                            y = predicted_mfi,
                                            color = "Standard Curve"),
            linewidth = 0.8) +
  geom_point(data = dplyr::filter(example_data, antigen == "hrp2"),
             aes(x = antigen_conc_iu_m_l, y = mfi,
                 color = "Data")) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual("", values = c("black", "deeppink3")) +
  ggtitle("Fitted standard curve for HRP2") +
  labs(x = "Antigen concentration (IU/ml)", y = "MFI") +
  theme_bw()
```

## Step 4 — Predict antigen concentration from MFI

Back-calculate antigen concentrations from given MFI values:
```r
predicted_antigen_conc <- pred_antigen_conc_from_mfi(mfi_values = c(100, 1000),
                                    model_params = hrp2_params)
```

Overlay predictions onto the standard curve plot to verify they fall on the curve:
```r
p1 + geom_point(data = predicted_antigen_conc, aes(x = antigen_conc_pred, y = mfi),
                size = 2, color = "dodgerblue1")
```


