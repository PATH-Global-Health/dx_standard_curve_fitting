#### FITTING STANDARD CURVES ####
# standard approach is a five parameter sigmoid curve f(x)=A+(D−A)/(1+exp(B(C−x)))^S
# https://www.r-bloggers.com/2019/11/five-parameters-logistic-regression/

getInitial5PL <- function(x, y){
  s <- getInitial(y ~ SSfpl(x, A, D, xmid, inverseB),
                  data = data.frame(x = x, y = y))
  c(A = s[["A"]], 
    B = 1/s[["inverseB"]], 
    xmid = s[["xmid"]], 
    D = s[["D"]], 
    S = 1)
}

fit5pl <- function(x, y){
  startingValues <- getInitial5PL(x, y)
  fit <- tryCatch({
    nlsLM(
      y ~ A + (D-A)/(1 + exp(log(2^(1/S)-1) + B*(xmid-x)))^S,
      data = data.frame(x = x, y = y),
      start = startingValues,
      lower = c(-Inf, 0, -Inf, -Inf, 0),
      control = nls.lm.control(maxiter = 1024, maxfev=10000))
  }, error = function(e){
    paste0("Failure of model fitting: ", e$message)
  })
  if(class(fit) == "nls" && fit[["convInfo"]][["isConv"]]){
    fit
  }else if(class(fit) == "nls" && !fit[["convInfo"]][["isConv"]]){
    "Convergence not achieved"
  }else{ # in this case, 'fit' is the error message
    fit
  }
}

# fit the model and extract the parameters
fit_params_function <- function(dat, antigen_name){
  
  # select antigen to fit
  dd = dat %>% filter(antigen == antigen_name)
  # fit the 5param sigmoid
  fit_mod = fit5pl(x = log10(dd$antigen_conc_iu_m_l),
                   y = log10(dd$mfi))
  
  # predict over a sequence of hrp2 concentratins to check fit
  mod_params = coef(fit_mod) %>% t() %>% data.frame()
  return(mod_params) 
  
}


# make the fitted standard curve

make_model_predictions <- function(data, fitted_params, antigen_conc_min, antigen_conc_max){
  
  antigen_conc_min = max(c(min(data$antigen_conc_iu_m_l), 0.001))
  antigen_conc_max = max(data$antigen_conc_iu_m_l)
  
  x_log10 = seq(log10(antigen_conc_min), log10(antigen_conc_max), by=0.05)
  preds_mod = data.frame(x_log10 = x_log10,
                         x = 10^x_log10,
                         y_log10 = fitted_params$A + (fitted_params$D-fitted_params$A)/(1 + exp(log(2^(1/fitted_params$S)-1) + fitted_params$B*(fitted_params$xmid-x_log10)))^fitted_params$S)
  
  preds_mod$y = 10^preds_mod$y_log10 
  pred_output <- preds_mod %>% 
    dplyr::select(antigen_conc = x,
                  predicted_mfi = y)
  
  return(pred_output)
}
  

# predict antigen concentration from a given MFI value

pred_antigen_conc_from_mfi <- function(mfi_values, model_params){
  y = log10(mfi_values)
  A = model_params$A
  D = model_params$D
  B = model_params$B
  xmid = model_params$xmid
  S = model_params$S
  
  antigen_conc_pred_log10 <- xmid - (1 / B) * log(((D - A) / (y - A))^(1 / S) - 1) +
    (1 / B) * log(2^(1 / S) - 1)
  
  return(data.frame(mfi = mfi_values, 
                    antigen_conc_pred = 10^antigen_conc_pred_log10))
}

