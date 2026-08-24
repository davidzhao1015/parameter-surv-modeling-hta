# ========================================================================
# This example script shows how to choose, fit, and validate parametric
# survival models obeying recommended framework by NICE TSD14 guide.
# ========================================================================


# === Setup ===
# Install required packages
if (!requireNamespace("survival", quietly = TRUE)) {
  install.packages("survival")
}
if (!requireNamespace("flexsurv", quietly = TRUE)) {
  install.packages("flexsurv")
}
# if (!requireNamespace("survminer", quietly = TRUE)) {
#   install.packages("survminer")
# }
if (!requireNamespace("GofCens", quietly = TRUE)) {
  install.packages("GofCens")
}

# Print versions of all the packages
packageVersion("survival") # v3.8.9
packageVersion("flexsurv") # v2.3.2
# packageVersion("survminer") # v0.5.2
packageVersion("GofCens") # 1.5


# Load required libraries
# Load required libraries
library("survival")
library("flexsurv")
# library("survminer")
library("GofCens")

library("styler") # For code formatting


# === Load data ===

# Control arm
ipd_ctrl <- read.csv("IPD_ctrl_gp.csv")

# Treatment arm
ipd_trt <- read.csv("IPD_treat_gp.csv")

head(ipd_ctrl)
head(ipd_trt)

# Merge datasets of two arms
ipd_merged <- rbind(ipd_ctrl, ipd_trt)

head(ipd_merged)



# === Log cumulative hazard plot ===

surv_km <- survival::survfit(Surv(time, status) ~ treat,
  data = ipd_merged,
  type = "kaplan-meier"
)

# Plot cumulative hazard (for Exponential distribution assumption)
# Show legend
old_par <- par(no.readonly = TRUE)

par(mfrow = c(1, 2)) # Set up a 1x2 plotting area

plot(surv_km,
  fun = "cumhaz",
  xlab = "Time (months)",
  ylab = "Cumulative hazard",
  col = c("black", "green"),
  main = "Cumulative hazard plot"
)

legend("topleft",
  legend = c("Control", "Treatment"),
  col = c("black", "green"),
  lty = 1,
  lwd = 1,
  bty = "n"
)

# Log both axis (for Weibull distribution assumption)
plot(surv_km,
  fun = "cumhaz", 
  log = "xy",
  xlab = "Log, Time (months)",
  ylab = "Log, Cumulative hazard",
  col = c("black", "green"),
  main = "Log cumulative hazard plot"
)

legend("topleft",
  legend = c("Control", "Treatment"),
  col = c("black", "green"),
  lty = 1,
  lwd = 1,
  bty = "n"
)
par(old_par) # Restore original plotting parameters






# === Overlap KM vs parameteric models (2 arms * 6 matrix) ===

data_list  <- list(
  ipd_ctrl = ipd_ctrl,
  ipd_trt = ipd_trt
)

distr_list <- c(
  "exponential",
  "weibull",
  "gompertz",
  "lognormal",
  "loglogistic"
)

pdf(
  file = "overlap_KM_vs_parametric_models.pdf",
  width = 10,
  height = 7,
  onefile = TRUE
)

for (data_name in names(data_list)) {
  data <- data_list[[data_name]]

  for (distr in distr_list) {
    cat("Fitting distribution:", distr, "for dataset:", data_name, "\n")

    fit <- flexsurv::flexsurvreg(Surv(time, status) ~ 1,
      data = data,
      dist = distr
    )

    plot(
      fit,
      type = "survival",
      ci = TRUE,
      main = paste(distr, "vs KM"),
      xlab = "Time (months)",
      ylab = "Survival probability"
    )
  }
}

dev.off()





# === Q-Q plot for right censored data ===

# Note: Q-Q plot is a graphical method for comparing two
# probability distributions by plotting their quantiles against
# each other. In this case, we will use it to assess if right censored data
# has a similar distribution to the fitted parametric model.


pdf(
  file = "survival_distribution_diagnostics.pdf",
  width = 10,
  height = 7,
  onefile = TRUE
)

for (distr in distr_list) {
  cat("Q-Q plot for distribution:", distr, "\n")

  tryCatch(
    {
      GofCens::probPlot(
        ipd_ctrl$time,
        ipd_ctrl$status,
        distr = distr,
        plots = c("PP", "QQ", "SP", "ER"),
        main = paste(distr, "diagnostic plots")
      )
    },
    error = function(error_message) {
      warning(
        "Distribution '", distr,
        "' is not available for GofCens::probPlot(): ",
        conditionMessage(error_message)
      )
    }
  )
}

dev.off()

# Todo: gompertz is not supported by GofCens::probPlot() function, so
# we need to find an alternative way to generate a Q-Q plot for the
# Gompertz distribution.







# === Residuals plot (Cox-snell) ===

cs <- flexsurv::coxsnell_flexsurvreg(fit_gompertz_ctrl)

surv <- survfit(Surv(cs$est, ipd_ctrl$status) ~ 1)
plot(surv, fun="cumhaz")
abline(0, 1, col="red")


# Todo: Add residuals plot for treatment arm and other distributions.








# ==== Fit six common parametric survival models with treatment arm ====

model_distr_list <- c(
  "exponential",
  "weibull",
  "gompertz",
  "lognormal",
  "llogis",
  "gamma"
)

fitted_models <- lapply(model_distr_list, function(distr) {
  cat("Fitting distribution:", distr, "\n")

  flexsurv::flexsurvreg(
    Surv(time, status) ~ treat,
    data = ipd_merged,
    dist = distr
  )
})

names(fitted_models) <- model_distr_list



# ---- Evaluate AIC of fitted models ----

aic_values <- sapply(fitted_models, function(model) {
  model$AIC
})

# Create a data frame to display AIC values for each distribution
aic_df <- data.frame(
  Distribution = names(aic_values),
  AIC = aic_values
)

# Drop row name
rownames(aic_df) <- NULL

# Sort the data frame by AIC values in ascending order
aic_df <- aic_df[order(aic_df$AIC), ]

aic_df # gompertz has the lowest AIC value




# --- Overlap KM vs best fitted model (Gompertz) ---

pdf(
  file = "overlap_KM_vs_best_fitted_models.pdf",
  width = 10,
  height = 7,
  onefile = TRUE
)

for (mod in names(fitted_models)) {

  cat("Fitted model:", mod, "\n")

  par(mfrow = c(1, 2))

  # Control arm

  km_ctrl <- survfit(Surv(time, status) ~ 1, 
  data = ipd_merged[ipd_merged$treat == 0, ])

  plot(
    km_ctrl,
    main = paste("Control arm: KM vs", mod),
    xlab = "Time",
    ylab = "Survival probability",
    lty = 1,
    col = "black",
    mark.time = FALSE
  )

  lines(
    fitted_models[[mod]],
    newdata = data.frame(treat = 0),
    type = "survival",
    ci = FALSE,
    lty = 1,
    col = "blue"
  )

  # Treatment arm

  km_trt <- survfit(Surv(time, status) ~ 1, 
  data = ipd_merged[ipd_merged$treat == 1, ])

  plot(
    km_trt,
    main = paste("Treatment arm: KM vs", mod),
    xlab = "Time",
    ylab = "Survival probability",
    lty = 2,
    mark.time = FALSE,
    col = "black"
  )

  lines(
    fitted_models[[mod]],
    newdata = data.frame(treat = 1),
    type = "survival",
    ci = FALSE,
    lty = 2,
    col = "blue"
  )

  legend(
    "bottomleft",
    legend = c(
      "Control - KM",
      "Control - Fitted",
      "Treatment - KM",
      "Treatment - Fitted"
    ),
    col = c("black", "blue", "black", "blue"),
    lty = c(1, 1, 2, 2),
    lwd = 1,
    bty = "n"
  )

}    

par(mfrow = c(1, 1)) # Reset plotting area to 1x1

dev.off()
