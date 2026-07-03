#================A project in comparing KM estimator and CIF===================
#=======================while a competing risk exists==========================

require(readr)
require(tidyr)
require(dplyr)
require(ggplot2)
require(survival)
require(survminer)
require(broom)
require(purrr)
require(cmprsk)

#====The event of interest is relapse while having a serious illness===========
#====Someone could face death while in the survey==============================
#====It's a major problem considering death as censoring=======================

#====CASE 1: No competing risk================================================

#A function to simulate data

simulate_cr_data<-function(n, lambda1, lambda2){
  t_relapse<-rexp(n, rate = lambda1)
  
  if (lambda2==0){
    t_death<- rep(Inf, n)
  } else {
    t_death<-rexp(n, rate=lambda2)
  }
  
  time<-pmin(t_relapse, t_death)
  
  status<-ifelse(t_relapse<t_death, 1,2)
  
  data.frame(
    time = time,
    status = status
  )
  
}

#Administrative cutoff at 5 years

apply_admin_censoring<-function(dat, tau=5){
  dat$time_obs<-pmin(dat$time, tau)
  dat$status_obs<-ifelse(dat$time<=tau, dat$status, 0)
  dat
}

#Getting some data for n=1000 patients (without competing risk)

dat<-simulate_cr_data(
  n=1000,
  lambda1=0.08,
  lambda2=0
)

head(dat)
#Admin censoring
dat<-apply_admin_censoring(dat)
head(dat)
summary(dat)
#Kaplan-Meier for these data
km_1<-survfit(Surv(time_obs,status_obs) ~ 1, data=dat)

summary(km_1)$table

ggsurvplot(km_1,
           conf.int = TRUE,
           xlab="Survival Time (in years)",
           ylab="Survival Probability",
           title="Kaplan-Meier Curve for no competing risk")

#Getting a prediction for the risk of relapse in 5 years
t0<-5

km_1_sum<-summary(km_1, times=5, extend=TRUE)
km_1_t0<-km_1_sum$surv
km_risk_t0<-1-km_1_t0

km_risk_t0

#Naive K-M estiamtes 5-year relapse risk around 33.4.


#CI for relapse risk

km_lower_risk<-1 - km_1_sum$upper
km_upper_risk<-1 - km_1_sum$lower

c(
  estimate = km_risk_t0,
  lower = km_lower_risk,
  upper = km_upper_risk
)

# 95$ CI for relapse risk in 5 years is [0.3, 0.36]


#----Computing CIF----------------------------------
require(cmprsk)

cif_1<-cuminc(
  ftime = dat$time_obs,
  fstatus = dat$status_obs,
  cencode = 0
)

names(cif_1)

cif_1_relapse<-cif_1[["1 1"]]
cif_1_relapse

#Getting CIF at 5 years

get_cif_at<-function(cif_component, t0){
  i <- findInterval(t0, cif_component$time)
  
  if (i==0){
    return(c(estimate=0, se=0))
  }
  
  est<-cif_component$est[i]
  se<-sqrt(cif_component$var[i])
  
  c(estimate = est, se = se)
}

cif_1_t0<-get_cif_at(cif_1_relapse, t0=5)

cif_1_t0

comparison<-data.frame(
  method = c("Naive Kaplan-Meier", "CIF"),
  estimate_5y = c(km_risk_t0, cif_1_t0['estimate'])
)
comparison

#=========It was expected that both K-M and CIF methods agree to each other.


#====CASE 2: Existence of competing risk (death) with λ_2=0.08=λ_1
require(purrr)


dat_2<-simulate_cr_data(n=1000, lambda1 = 0.08, lambda2 = 0.08)
head(dat_2)
dat_2<-apply_admin_censoring(dat_2)
head(dat_2)
summary(dat_2)
relapse_event<-ifelse(dat_2$status_obs==1, 1,0)
dat_2$relapse_event<-relapse_event
#Apply naive K-M estimate
km_2<-survfit(Surv(time_obs, relapse_event)~1, data=dat_2)

ggsurvplot(km_2,
           conf.int = TRUE,
           xlab='Time (in Years)',
           ylab="Survival Prob.",
           title="K-M where death is trated as censoring")

km_2_sum<-summary(km_2, times=t0, extend=TRUE)
km_2_sum$surv
km_2_risk<-1 - km_2_sum$surv
km_2_risk


#CIF estimate

estimate_cif_risk<-function(dat, t0, cause=1){
  
  cif_fit<-cuminc(
    ftime=dat$time_obs,
    fstatus=dat$status_obs,
    cencode = 0
  )
  
  cif_names<-names(cif_fit)
  cif_names<-cif_names[cif_names != "Tests"]
  
  cause_name<-cif_names[grepl(paste0(" ", cause, "$"), cif_names)][1]
  
  cif_component<-cif_fit[[cause_name]]
  
  i<-findInterval(t0, cif_component$time)
  
  if (i==0){
    return(0)
  }
  
  as.numeric(cif_component$est[i])
}


true_cif<-function(t0, lambda1, lambda2){
  
  if (lambda2==0){
    return(1-exp(-lambda1*t0))
  }
  
  lambda1/(lambda1+lambda2)*(1-exp(-(lambda1+lambda2)*t0))
}

table(dat_2$status_obs)
round(prop.table(table(dat_2$status_obs)),2)

cif_est<-estimate_cif_risk(dat_2, t0=5, cause=1)

truth<-true_cif(t0=5, lambda1=0.08, lambda2=0.08)

data.frame(
  method=c("Naive K-M", "CIF", "True CIF"),
  estimate=c(km_2_risk, cif_est, truth)
)


#----Plot both K-M and CIF in one dataset

plot_one_dataset <- function(dat, t0 = 5) {
  
  # Naive KM
  relapse_event <- ifelse(dat$status_obs == 1, 1, 0)
  
  km_fit <- survfit(
    Surv(time_obs, relapse_event) ~ 1,
    data = dat
  )
  
  km_df <- data.frame(
    time = km_fit$time,
    estimate = 1 - km_fit$surv,
    method = "Naive Kaplan-Meier"
  )
  
  # CIF
  cif_fit <- cuminc(
    ftime = dat$time_obs,
    fstatus = dat$status_obs,
    cencode = 0
  )
  
  cif_names <- names(cif_fit)
  cif_names <- cif_names[cif_names != "Tests"]
  cause1_name <- cif_names[grepl(" 1$", cif_names)][1]
  
  cif_relapse <- cif_fit[[cause1_name]]
  
  cif_df <- data.frame(
    time = cif_relapse$time,
    estimate = cif_relapse$est,
    method = "CIF"
  )
  
  plot_df <- rbind(km_df, cif_df)
  
  ggplot(plot_df, aes(x = time, y = estimate, linetype = method)) +
    geom_step(linewidth = 1) +
    coord_cartesian(xlim = c(0, t0), ylim = c(0, NA)) +
    labs(
      x = "Time",
      y = "Estimated relapse probability",
      linetype = "Method",
      title = "Naive Kaplan-Meier vs CIF",
      subtitle = "Competing event: death before relapse"
    ) +
    theme_minimal()
}

plot_one_dataset(dat_2, t0=5)


#=====Let's summarize and create a model that does everything done above.

#1. A function to simulate data with competing risks

simulate_cr_data<-function(n, lambda1, lambda2, tau=5){
  
  if (lambda1<=0) stop("lambda1 must be positive.")
  if (lambda2<0) stop("lambda2 must be non=negative.")
  
  t_relapse<-rexp(n, rate=lambda1)
  
  
  if (lambda2==0){
    t_death<-rep(Inf, n)
  } else {
    t_death<-rexp(n, rate=lambda2)
  }
  
  time_true<-pmin(t_relapse, t_death)
  
  status_true<-ifelse(t_relapse < t_death, 1,2)
  
  time_obs<-pmin(time_true, tau)
  status_obs<-ifelse(time_true<=tau, status_true, 0)
  
  data.frame(
    time_obs = time_obs,
    status_obs = status_obs
  )
}

#2. Naive Kaplan-Meier. Death as censoring

estimate_km_risk<-function(dat, t0=5){
  
  relapse_event<-ifelse(dat$status_obs ==1, 1,0)
  
  km_fit<-survfit(
    Surv(time_obs, relapse_event) ~ 1, 
    data = dat
  )
  
  km_sum<-summary(km_fit, times=t0, extend=TRUE)
  
  risk<-1- km_sum$surv
  
  as.numeric(risk)
}


#3. CIF Estimate

estimate_cif_risk<-function(dat, t0=5, cause =1 ){
  
  cif_fit<- cuminc(
    ftime = dat$time_obs,
    fstatus = dat$status_obs,
    cencode = 0
  )
  
  cif_names<-names(cif_fit)
  cif_names<-cif_names[cif_names != "Tests"]
  
  pattern<-paste0("(^|\\s)", cause, "$")
  cause_name<-cif_names[grepl(pattern, cif_names)][1]
  
  if (is.na(cause_name)) {
    return(NA_real_)
  }
  
  cif_component<-cif_fit[[cause_name]]
  
  i<-findInterval(t0, cif_component$time)
  
  if (i==0) {
    return(0)
  }
  
  as.numeric(cif_component$est[i])
  
  
}

#4. Analytical true CIF

true_cif<-function(t0, lambda1, lambda2){
  
  if(lambda2==0){
    return(1-exp(-lambda1*t0))
  }
  
  lambda1/(lambda1 + lambda2) * (1-exp(-(lambda1+lambda2)*t0))
}

#6. A simulation repetition
run_one_replication<-function(n, lambda1, lambda2, t0=5){
  
  dat<-simulate_cr_data(
    n=n,
    lambda1 = lambda1,
    lambda2 = lambda2,
    tau = t0
  )
  
  km_est<-estimate_km_risk(dat, t0=t0)
  cif_est<-estimate_cif_risk(dat, t0=t0, cause=1)
  truth<-true_cif(t0, lambda1, lambda2)
  
  data.frame(
    n = n,
    t0 = t0,
    lambda1 = lambda1,
    lambda2 = lambda2,
    truth = truth,
    km_est = km_est,
    cif_est = cif_est,
    km_bias = km_est - truth,
    cif_bias = cif_est - truth
  )
}


#7. 1000 simulations for a scenario

run_scenario<-function(B = 1000, n, lambda1, lambda2, t0=5){
  
  map_dfr(seq_len(B), function(b){
    
    out<-run_one_replication(
      n=n,
      lambda1 = lambda1,
      lambda2 = lambda2,
      t0 =t0
    )
    
    out$replication<-b
    
    out
  })
}

#---An example---------
set.seed(123)

res_008<-run_scenario(
  B=1000,
  n=1000,
  lambda1=0.08,
  lambda2 = 0.08,
  t0=5
)
head(res_008)


#8. 1000 simulations for many lambda2 values

run_simulation_grid<-function(
    B=1000,
    n_values=1000,
    lambda1_values = 0.08,
    lambda2_values=c(0, 0.04, 0.08, 0.16),
    t0_values=5,
    seed=NULL
){
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  grid<-expand.grid(
    n= n_values,
    lambda1 = lambda1_values,
    lambda2 = lambda2_values,
    t0 = t0_values
  )
  
  
  results<-pmap_dfr(
    grid,
    function(n, lambda1, lambda2, t0){
      
      run_scenario(
        B= B,
        n=n,
        lambda1=lambda1,
        lambda2=lambda2,
        t0=t0
      )
      
    }
  )
  
  results
}


results<-run_simulation_grid(
  B=1000,
  n_values = 1000,
  lambda1_values = 0.08,
  lambda2_values=c(0, 0.04, 0.08, 0.16),
  t0_values=5,
  seed=123
)

head(results)


#9. Summary table

summarise_simulation_results<-function(results){
  
  results %>%
    group_by(n, t0, lambda1, lambda2)%>%
    summarise(
      B=n(),
      
      truth=first(truth),
      
      mean_km=mean(km_est, na.rm=TRUE),
      mean_cif=mean(cif_est, na.rm=TRUE),
      
      bias_km=mean(km_bias, na.rm=TRUE),
      bias_cif=mean(cif_bias, na.rm=TRUE),
      
      relative_bias_km=bias_km/truth,
      relative_bias_cif=bias_cif/truth,
      
      sd_km=sd(km_est, na.rm=TRUE),
      sd_cif=sd(cif_est, na.rm=TRUE),
      
      mse_km=mean((km_est-truth)^2, na.rm=TRUE),
      mse_cif=mean((cif_est-truth)^2, na.rm=TRUE),
      
      .groups="drop"
    )
}

summary_results<-summarise_simulation_results(results)

summary_results

summary_results %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

#10. Everything you need in one place
results <- run_simulation_grid(
  B = 1000,
  n_values = 1000,
  lambda1_values = 0.08,
  lambda2_values = c(0, 0.04, 0.08, 0.16),
  t0_values = 5,
  seed = 123
)

summary_results <- summarise_simulation_results(results)

summary_results %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

#11. Plot estimates

plot_estimates_vs_lambda2<-function(summary_results){
  
  plot_df<-summary_results %>%
    select(lambda2, truth, mean_km, mean_cif)%>%
    pivot_longer(
      cols=c(truth, mean_km, mean_cif),
      names_to="method",
      values_to="estimate"
    )%>%
    mutate(
      method = recode(
        method,
        truth="True CIF",
        mean_km="Naive Kaplan-Meier",
        mean_cif='CIF'
      )
    )
  
  ggplot(plot_df, aes(x=lambda2, y=estimate, linetype=method))+
    geom_line(linewidth=1)+
    geom_point(size=2)+
    labs(
      title = "Estimated 5-year relapse risk",
      subtitle = "Naive Kaplan-Meier overestimates risk as competing risk increases",
      x=expression(lambda[2] ~ "(death before relapse rate)"),
      y="Estimated relapse probability",
      linetype="Method"
    )+
    theme_minimal()
}

plot_estimates_vs_lambda2(summary_results)

#12. Plot bias vs λ_2

plot_bias_vs_lambda2<-function(summary_results){
  
  plot_df<-summary_results %>%
    select(lambda2, bias_km, bias_cif)%>%
    pivot_longer(
      cols=c(bias_km, bias_cif),
      names_to="method",
      values_to="bias"
    ) %>%
    mutate(
      method=recode(
        method,
        bias_km="Naive Kaplan-Meier",
        bias_cif="CIF"
      )
    )
  
  ggplot(plot_df, aes(x=lambda2, y=bias, linetype=method))+
    geom_line(linewidth=1)+
    geom_point(size=2)+
    geom_hline(yintercept = 0)+
    labs(
      title="Bias in estimated 5-year relapse risk",
      subtitle="Bias of naive Kaplan-Meier increases with the competing risk event rate",
      x=expression(lambda[2]~"(death before relapse rate)"),
      y="Bias",
      linetype="Method"
    )+
    theme_minimal()
}

plot_bias_vs_lambda2(summary_results)


#=====Extensions: Let's work with some more lambda2 values

results_more<-run_simulation_grid(
  B=1000,
  n_values = 1000,
  lambda1_values = 0.08,
  lambda2_values = seq(0, 0.3, by=0.02),
  t0_values=5,
  seed=123
)

summary_more<-summarise_simulation_results(results_more)

plot_estimates_vs_lambda2(summary_more)
plot_bias_vs_lambda2(summary_more)


#======Testing for different sample sizes

results_n<-run_simulation_grid(
  B=1000,
  n_values=c(250, 500, 1000, 2000),
  lambda1_values = 0.08,
  lambda2_values = c(0,0.04,0.08, 0.16),
  t0_values=5,
  seed=123
)

summary_n<-summarise_simulation_results(results_n)
plot_estimates_vs_lambda2(summary_n)
plot_bias_vs_lambda2(summary_n)

sample_size_table <- summary_n |>
  transmute(
    n,
    truth,
    mean_km,
    mean_cif,
    bias_km,
    bias_cif,
    sd_km,
    sd_cif
  ) |>
  mutate(across(where(is.numeric), ~ round(.x, 4)))

sample_size_table

kable(
  sample_size_table_report,
  caption = "Sample size sensitivity analysis.",
  col.names = c(
    "n",
    "True CIF",
    "Naive KM",
    "CIF",
    "KM bias",
    "CIF bias",
    "SD KM",
    "SD CIF"
  )
)

#====Now that we've created the model, let's get some results

results<-run_simulation_grid(
  B=1000,
  n_values = 1000,
  lambda1_values = 0.08,
  lambda2_values = c(0,0.04, 0.08, 0.16),
  t0_values = 5,
  seed=123
)

summary_results<-summarise_simulation_results(results)

summary_results_rounded<-summary_results %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

#Diagnostic table
diagnostic_table<-summary_results %>%
  transmute(
    lambda2,
    true_5y_risk=truth,
    km_5y_risk=mean_km,
    cif_5y_risk=mean_cif,
    km_absolute_overestimation=mean_km-truth,
    km_relative_overestimation=100*(mean_km-truth)/truth
  )%>%
  mutate(across(where(is.numeric), ~round(.x, 4)))

diagnostic_table$km_relative_overestimation

(estimate_plot<-plot_estimates_vs_lambda2(summary_results))
(bias_plot<-plot_bias_vs_lambda2(summary_results))


ggsave(
  filename = "estimate_vs_lambda2.png",
  plot = estimate_plot,
  width=7,
  height=5,
  dpi=300
)

ggsave(
  filename = "bias_vs_lambda2.png",
  plot= bias_plot,
  width=7,
  height=5,
  dpi=300
)

write.csv(
  summary_results_rounded,
  file = "summary_results.csv",
  row.names = FALSE
)

write.csv(
  diagnostic_table,
  file = "diagnostic_table.csv",
  row.names=FALSE
)

write.csv(
  sample_size_table,
  file = "sample_size_table.csv",
  row.names=FALSE
)


#A report ready table
report_table <- summary_results %>%
  transmute(
    lambda2,
    true_5y_risk = truth,
    naive_km_5y_risk = mean_km,
    cif_5y_risk = mean_cif,
    km_abs_bias = bias_km,
    cif_abs_bias = bias_cif,
    km_relative_bias_percent = 100 * relative_bias_km,
    cif_relative_bias_percent = 100 * relative_bias_cif
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

report_table

write.csv(
  report_table,
  file ="report_table.csv",
  row.names=FALSE
)














