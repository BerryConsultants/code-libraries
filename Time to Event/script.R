
source('utils_tte.R')
key_arms = c('Control', 'Dose 1', 'Dose 2', 'Dose 3', 'Dose 4')
key_palette = c('black', brewer.pal(9, 'Set1')[1:(length(key_arms) - 1)])

data_seed = 12345

s_j_true = c(0, 5, 10, 15, 20) # cutpoints
lambda_true = c(0.06, 0.06, 0.03, 0.03) # baseline hazard rates
theta_true = log(c(1.5, 2, 2.5, 1.75)) # hazard ratios

N = 250 # sample size 
D = length(theta_true) # number of doses 

rnd_ratio = c(4, 2, 2, 2, 2) # randomization ratio
blk_size = sum(rnd_ratio) # block size

# Generate the random allocations
set.seed(data_seed)
dose_label = tibble(BLOCK = c(rep(1:floor(N/blk_size), each = blk_size), 
                              rep(floor(N/blk_size) + 1, N - length(rep(1:floor(N/blk_size), 
                                                                        each = blk_size))))) %>% 
  group_by(BLOCK) %>% 
  mutate(TRTNUM = sample(rep(0:D, times = rnd_ratio), length(BLOCK))) %>% 
  ungroup() %>% 
  pull(TRTNUM)

# Put the allocations into the dose matrix format
X = matrix(0, N, D)
for (s in 1:D){
  X[which(dose_label == s),s] = 1
}

# Show X in output
X %>% 
  colnames_inplace(sprintf('Dose %s', 1:D)) %>% 
  head(n = 6) %>% 
  knitr::kable(booktabs = TRUE, 
               caption = 'First rows of the dose allocation matrix X.') %>% 
  kableExtra::kable_styling(bootstrap_options = c("striped", "hover"))

accr_profile = NULL
accr_profile$peak_rate = 3 # in weeks
accr_profile$ramp_up = 0 # in weeks

J = length(lambda_true)
max_followup = 20
t_grid = seq(0, 20, length.out = 500)

# Now we use the inverse cdf method to sample from the likelihood model
u = c(0.952, 0.657, 0.537, 0.444, 0.142)

# Check which time interval the observation should fall into
s_j_temp = s_j_true
s_j_temp[J+1] = +Inf
H_temp = c(0, cumsum(lambda_true * diff(s_j_temp)))
time = event = interval_idx = array(NA, dim = length(u))
for (i in 1:length(u)) {
  interval_idx[i] = head(which(log(u[i]) >= - H_temp), 1) - 1
  
  if (length(interval_idx[i]) == 0){ # censored at the observation window endpoint
    time[i] = max_followup
    event[i] = 0
  }
  else { # sampled from the model in that interval
    event[i] = 1
    time[i] = s_j_true[interval_idx[i]] - 1/lambda_true[interval_idx[i]] * (log(u[i]) + H_temp[interval_idx[i]])
    if (time[i] > max_followup){ # it still could be censored in that interval
      time[i] = max_followup
      event[i] = 0
    }
  }
}

i = 1
ggplot() +
  geom_line(aes(x = t_grid, y = St(t_grid, s_j_true, lambda_true)), size = 1) +
  geom_hline(yintercept = u[i], size = 1, lty = 2, col = 'red') +
  geom_vline(xintercept = time[i], size = 1, col = 'red', lty = 2) +
  geom_point(aes(x = time[i], y = u[i]), col = 'red', size = 4) + 
  geom_vline(xintercept = s_j_true, lty = 3) +
  scale_y_continuous(name = TeX("$S(t)$"), limits = c(0, 1), expand = c(0, 0), 
                     breaks = c(0, 0.25, 0.50, 0.75, 1.00, u[i]), labels = c('0.00', '0.25', '0.50', '0.75', '1.00', 'u')) + 
  scale_x_continuous(name = 'time (weeks)', limits = c(0, 30), expand = c(0, 0), 
                     breaks = c(0, 10, 20, 30, time[i]), labels = c('0', '10', '20', '30', 't')) + 
  theme(text = element_text(size = 16))
i = 2
ggplot() +
  geom_line(aes(x = t_grid, y = St(t_grid, s_j_true, lambda_true)), size = 1) +
  geom_hline(yintercept = u[i], size = 1, lty = 2, col = 'red') +
  geom_vline(xintercept = time[i], size = 1, col = 'red', lty = 2) +
  geom_point(aes(x = time[i], y = u[i]), col = 'red', size = 4) + 
  geom_vline(xintercept = s_j_true, lty = 3) +
  scale_y_continuous(name = TeX("$S(t)$"), limits = c(0, 1), expand = c(0, 0), 
                     breaks = c(0, 0.25, 0.50, 0.75, 1.00, u[i]), labels = c('0.00', '0.25', '0.50', '0.75', '1.00', 'u')) + 
  scale_x_continuous(name = 'time (weeks)', limits = c(0, 30), expand = c(0, 0), 
                     breaks = c(0, 10, 20, 30, time[i]), labels = c('0', '10', '20', '30', 't')) + 
  theme(text = element_text(size = 16))
i = 3
ggplot() +
  geom_line(aes(x = t_grid, y = St(t_grid, s_j_true, lambda_true)), size = 1) +
  geom_hline(yintercept = u[i], size = 1, lty = 2, col = 'red') +
  geom_vline(xintercept = time[i], size = 1, col = 'red', lty = 2) +
  geom_point(aes(x = time[i], y = u[i]), col = 'red', size = 4) + 
  geom_vline(xintercept = s_j_true, lty = 3) +
  scale_y_continuous(name = TeX("$S(t)$"), limits = c(0, 1), expand = c(0, 0), 
                     breaks = c(0, 0.25, 0.50, 0.75, 1.00, u[i]), labels = c('0.00', '0.25', '0.50', '0.75', '1.00', 'u')) + 
  scale_x_continuous(name = 'time (weeks)', limits = c(0, 30), expand = c(0, 0), 
                     breaks = c(0, 10, 20, 30, time[i]), labels = c('0', '10', '20', '30', 't')) + 
  theme(text = element_text(size = 16))
i = 4
ggplot() +
  geom_line(aes(x = t_grid, y = St(t_grid, s_j_true, lambda_true)), size = 1) +
  geom_hline(yintercept = u[i], size = 1, lty = 2, col = 'red') +
  geom_vline(xintercept = time[i], size = 1, col = 'red', lty = 2) +
  geom_point(aes(x = time[i], y = u[i]), col = 'red', size = 4) + 
  geom_vline(xintercept = s_j_true, lty = 3) +
  scale_y_continuous(name = TeX("$S(t)$"), limits = c(0, 1), expand = c(0, 0), 
                     breaks = c(0, 0.25, 0.50, 0.75, 1.00, u[i]), labels = c('0.00', '0.25', '0.50', '0.75', '1.00', 'u')) + 
  scale_x_continuous(name = 'time (weeks)', limits = c(0, 30), expand = c(0, 0), 
                     breaks = c(0, 10, 20, 30, time[i]), labels = c('0', '10', '20', '30', 't')) + 
  theme(text = element_text(size = 16))
i = 5
ggplot() +
  geom_line(aes(x = t_grid, y = St(t_grid, s_j_true, lambda_true)), size = 1) +
  geom_hline(yintercept = u[i], size = 1, lty = 2, col = 'red') +
  geom_vline(xintercept = time[i], size = 1, col = 'red', lty = 2) +
  geom_point(aes(x = time[i], y = u[i]), col = 'red', size = 4) + 
  geom_vline(xintercept = s_j_true, lty = 3) +
  scale_y_continuous(name = TeX("$S(t)$"), limits = c(0, 1), expand = c(0, 0), 
                     breaks = c(0, 0.25, 0.50, 0.75, 1.00, u[i]), labels = c('0.00', '0.25', '0.50', '0.75', '1.00', 'u')) + 
  scale_x_continuous(name = 'time (weeks)', limits = c(0, 30), expand = c(0, 0), 
                     breaks = c(0, 10, 20, 30, time[i]), labels = c('0', '10', '20', '30', 't')) + 
  theme(text = element_text(size = 16))

set.seed(data_seed)
data = sample_pwc_data(n = N, # sample size
                       X = X, # dose allocation matrix
                       lambda = lambda_true, # baseline hazard rates
                       theta = theta_true, # log hazard ratios
                       s_j = s_j_true, # cutpoints 
                       T_max = 0, # follow-up after last randomized patient
                       accr_profile = accr_profile, # accrual profile
                       t0 = 0 # current week (0 for beginning of time)
)

# Censor the observations at the maximum follow-up
data$event[which(data$time > 20)] = 0
data$time[which(data$time > 20)] = 20

y = as.numeric(data$time) # either censoring time (if event == 0) or death time (if event == 1)
nu = as.numeric(data$event) # event indicator

# True survival curves
tibble(X = rep(t_grid, length(key_arms)), 
       y = c(St(t_grid, s_j_true, lambda_true), 
             St(t_grid, s_j_true, exp(theta_true[1]) * lambda_true), 
             St(t_grid, s_j_true, exp(theta_true[2]) * lambda_true), 
             St(t_grid, s_j_true, exp(theta_true[3]) * lambda_true), 
             St(t_grid, s_j_true, exp(theta_true[4]) * lambda_true)), 
       TRTPN = factor(rep(key_arms, each = length(t_grid)), levels = key_arms)) %>% 
  ggplot() + 
  geom_line(aes(x = X, y = y, col = TRTPN), size = 1) + 
  geom_vline(xintercept = s_j_true, lty = 3, size = 1) + 
  scale_x_continuous(name = 'Follow-up time (weeks)', expand = c(0, 0)) + 
  scale_y_continuous(name = TeX("$S(t)$"), limits = c(0, 1), expand = c(0, 0)) +
  scale_color_manual(name = 'Arm', values = key_palette) +
  theme(legend.position = 'top', text = element_text(size = 16))

# KM Plot
dat_KM = tibble(Duration = y,
                Outcome = nu,
                Dose = factor(dose_label))
fit_KM <- survfit(Surv(Duration, Outcome) ~ Dose, data = dat_KM)
p_km = ggsurvplot(fit_KM,
                  data = dat_KM,
                  size = 1, 
                  legend.title = "Arm",
                  legend.labs = key_arms,
                  palette = key_palette,
                  axes.offset = F,
                  ggtheme = theme_bw(base_size = 16)
)
p_km$plot +
  geom_vline(xintercept = s_j_true, size = 1, lty = 3) +
  labs(x = 'Follow-up time (weeks)', y = TeX("$S(t)$"))

tibble(TIMESINCERAND = data$dates) %>%   
  mutate(Subject = row_number()) %>% 
  mutate(EXPACCR = Lambda(t = TIMESINCERAND, peak_rate = accr_profile$peak_rate, 
                          ramp_up = accr_profile$ramp_up)) %>% 
  ggplot() + 
  geom_step(aes(x = TIMESINCERAND, y = Subject, col = 'observed'), size = 1) + 
  geom_line(aes(x = TIMESINCERAND, y = EXPACCR, 
                col = paste0(accr_profile$peak_rate, ' per week')), 
            size = 1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'Weeks since start of accrual', 
       y = 'Cumulative number of subjects', 
       title = 'Observed vs simulated accrual') + 
  scale_color_brewer(name = '', palette = 'Dark2') + 
  theme(legend.position = "top")

model_seed = 11111

# Define interval breaks (in weeks)
eps = 1E-6
s_j = c(0, eps + c(10, 20)) 
J = length(s_j) - 1

lambda_grid = seq(0, max(rgamma(1000, 1, 1 / round(- log (1 - 0.4) / 20, 4))), length.out = 200)
ggplot() +
  geom_line(aes(x = lambda_grid, y = dgamma(lambda_grid, 1, 1 / round(- log (1 - 0.4) / 20, 4))), col = key_palette[1], size = 2) +
  geom_vline(xintercept = round(- log (1 - 0.4) / 20, 4), col = 'red', lty = 3, size = 1) +
  scale_x_continuous(name = TeX("$\\lambda_{j}$"), expand = c(0, 0)) +
  scale_y_continuous(name = 'density', expand = c(0, 0)) +
  theme(text = element_text(size = 16))

ggplot() + 
  geom_line(aes(x = t_grid, y = St(t_grid, c(0, 10, 20), rep(round(- log (1 - 0.4) / 20, 4), 2))), col = key_palette[1], size = 2) +
  geom_point(aes(x = 20, y = 0.6), col = 'red', size = 3) +
  geom_vline(xintercept = c(0, 10, 20), lty = 3, size = 1) + 
  geom_hline(yintercept = 0.60, col = 'red', lty = 3, size = 1) +
  scale_x_continuous(name = 'time (weeks)', limits = c(0, 30), expand = c(0, 0)) + 
  scale_y_continuous(name = TeX("$S(t)$"), limits = c(0, 1), expand = c(0, 0), 
                     breaks = c(0, 0.25, 0.5, 0.6, 0.75, 1)) +
  theme(text = element_text(size = 16))

# Prior hazard hyperparameters
a_j = rep(1, J)
b_j = rep(1/round(- log (1 - 0.4) / 20, 4), J)

alpha = 1
beta = 1

# Calculate summary statistics
summ_stats = get_summary_stats(y, nu, s_j)

# Put data in list for STAN
stan_data <- list(N = N,
                  X = X,
                  D = D,
                  J = J,
                  delta = summ_stats$delta_i,
                  H = summ_stats$Hij,
                  a_j = a_j,
                  b_j = b_j, 
                  alpha = alpha, 
                  beta = beta
)


# Run STAN model
n_iters = 5000
burnin = 2500
n_chains = 2
set.seed(model_seed)
fit <- stan(file = './TTE_pwc.stan',
            data = stan_data,
            chains = n_chains,
            warmup = burnin,
            iter = n_iters,
            seed = model_seed) %>%
  extract()

# Extract model output
samples = fit$lambda %>%
  cbind(fit$theta) %>%
  cbind(fit$tau2) %>%
  colnames_inplace(c(sprintf('lambda[%s]', 1:J), sprintf('theta[%s]', 1:D), 'tau^2'))

# Traceplots
mcmc_trace(samples, facet_args = list(labeller = ggplot2::label_parsed, 
                                      ncol = 2))

# Effective sample size
print(round(coda::effectiveSize(samples)))

samples %>% 
  as_tibble %>% 
  pivot_longer(cols = everything(), names_to = 'Parameter') %>% 
  group_by(Parameter) %>% 
  summarize(Mean = mean(value), 
            Sd = sd(value),
            Median = median(value), 
            Low = quantile(value, probs = 0.025), 
            Upp = quantile(value, probs = 0.975), 
            ESS = round(coda::effectiveSize(value))) %>% 
  mutate(Mean = sprintf("$%0.3f$", Mean), 
         Sd = sprintf("$%0.3f$", Sd),
         MCI = sprintf("$%0.3f$ ($%0.3f$, $%0.3f$)", Median, Low, Upp), 
         ESS = sprintf("$%s$", ESS)) %>% 
  select(Parameter, Mean, `Std. Dev.` = Sd, `Median (95% CI)` = MCI, 
         `Effective Sample Size` = ESS) %>% 
  mutate(Parameter = str_replace_all(Parameter, '\\[', '\\_'), 
         Parameter = str_remove_all(Parameter, '\\]'), 
         Parameter = paste0("$\\", Parameter, "$")) %>% 
  knitr::kable(booktabs = TRUE, 
               caption = 'Posterior estimates of the model parameters.') %>% 
  kableExtra::kable_styling(bootstrap_options = c("striped", "hover"))

# Survival curves for each treatment arm
post_surv_ctr = St(t = t_grid,
                   s_j = s_j,
                   lambda = samples[,1:J])
post_surv_trt1 = St(t = t_grid,
                    s_j = s_j,
                    lambda = exp(samples[,J+1]) * samples[,1:J])
post_surv_trt2 = St(t = t_grid,
                    s_j = s_j,
                    lambda = exp(samples[,J+2]) * samples[,1:J])
post_surv_trt3 = St(t = t_grid,
                    s_j = s_j,
                    lambda = exp(samples[,J+3]) * samples[,1:J])
post_surv_trt4 = St(t = t_grid,
                    s_j = s_j,
                    lambda = exp(samples[,J+4]) * samples[,1:J])
post_survival = post_surv_ctr %>%
  reshape2::melt(varnames = c('t', 'iteration')) %>%
  mutate(t = t_grid[t],
         TRTPN = key_arms[1]) %>%
  bind_rows(post_surv_trt1 %>%
              reshape2::melt(varnames = c('t', 'iteration')) %>%
              mutate(t = t_grid[t],
                     TRTPN = key_arms[2])) %>%
  bind_rows(post_surv_trt2 %>%
              reshape2::melt(varnames = c('t', 'iteration')) %>%
              mutate(t = t_grid[t],
                     TRTPN = key_arms[3])) %>%
  bind_rows(post_surv_trt3 %>%
              reshape2::melt(varnames = c('t', 'iteration')) %>%
              mutate(t = t_grid[t],
                     TRTPN = key_arms[4])) %>%
  bind_rows(post_surv_trt4 %>%
              reshape2::melt(varnames = c('t', 'iteration')) %>%
              mutate(t = t_grid[t],
                     TRTPN = key_arms[5])) %>%
  group_by(t, TRTPN) %>%
  summarize(mean = mean(value),
            low = quantile(value, probs = 0.025),
            upp = quantile(value, probs = 0.975)) %>%
  ungroup() %>%
  mutate(TRTPN = factor(TRTPN, levels = key_arms))


# KM Estimate
dat_KM = tibble(ExpTime = y, PRIM = nu, Dose = factor(key_arms[dose_label + 1]))
fit_KM <- survfit(Surv(ExpTime, PRIM) ~ Dose, data = dat_KM)
p_km = ggsurvplot(fit_KM,
                  data = dat_KM,
                  palette = key_palette,
                  legend = "top",
                  legend.title = "Arm",
                  legend.labs = levels(dat_KM$Dose),
                  axes.offset = FALSE,
                  break.time.by = 10,
                  ggtheme = custom_theme()
)

# Posterior survival curves and KM plot
p_km = p_km$plot +
  geom_line(data = post_survival, aes(x = t, y = mean, col = TRTPN), size = 1, lty = 4) +
  # geom_ribbon(data = post_survival, 
  #             aes(x = t, ymin = low, ymax = upp, fill = TRTPN), 
  #             alpha = 0.25, inherit.aes = FALSE) +
  geom_vline(xintercept = s_j, lty = 3) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  labs(x = 'Follow-up time (weeks)', y = 'Survival probability') +
  theme(legend.position = 'top')
p_km

# Dose response curve (HR density for each dose)
samples %>%
  as_tibble() %>%
  select(contains('theta')) %>%
  add_column('theta[0]' = 0, .before = 'theta[1]') %>%
  pivot_longer(cols = everything(), names_to = 'TRTPN', values_to = 'HR',
               names_prefix = 'theta') %>%
  mutate(TRTPN = as.numeric(str_remove_all(TRTPN, "\\[|\\]")), 
         HR = exp(HR),
         TRTPN = factor(key_arms[as.numeric(TRTPN) + 1], levels = key_arms)) %>%
  group_by(TRTPN) %>%
  summarize(low = quantile(HR, probs = 0.025),
            upp = quantile(HR, probs = 0.975),
            mean = mean(HR),
            median = median(HR)) %>%
  ggplot() +
  geom_point(aes(x = TRTPN, y = mean), pch = 10, size = 4, stroke = 1) +
  geom_line(aes(x = TRTPN, y = mean, group = 1), size = 1) +
  geom_line(aes(x = TRTPN, y = low, group = 1), lty = 3) +
  geom_line(aes(x = TRTPN, y = upp, group = 1), lty = 3) +
  geom_point(aes(x = 1:5, y = c(1, exp(theta_true))), col = 'red') + 
  geom_ribbon(aes(x = TRTPN, ymin = low, ymax = upp, group = 1), alpha = 0.25) +
  labs(x = '', y = 'Hazard ratios')

# Baseline hazard function
reshape2::melt(ht(t_grid, s_j, samples[,1:J]), varnames = c('Iteration', 'X')) %>%
  mutate(X = t_grid[X]) %>%
  group_by(X) %>%
  dplyr::summarize(ybar = mean(value),
                   UI = quantile(value, probs = 0.975),
                   LI = quantile(value, probs = 0.025)) %>%
  ggplot() +
  geom_line(aes(x = X, y = ybar), size = 1) +
  geom_ribbon(aes(x = X, ymin = LI, ymax = UI), alpha = 0.3) +
  geom_step(data = data.frame(s_j_true = s_j_true, 
                              lambda_true = c(lambda_true, tail(lambda_true, 1))),
            aes(x = s_j_true, y = lambda_true), size = 1, col = 'red') +
  geom_vline(xintercept = s_j, lty = 3) +
  scale_x_continuous(name = 'time (weeks)', expand = c(0, 0)) + 
  labs(y = TeX("$\\lambda_{j}$"))

# Define interval breaks (in weeks)
s_j = c(0, eps + c(10, 20)) 
J = length(s_j) - 1

# Prior hazard hyperparameters
a_j = rep(1, J)
b_j = rep(1/round(- log (1 - 0.4) / 20, 4), J)

alpha = 10
beta = 0.1

# Calculate summary statistics
summ_stats = get_summary_stats(y, nu, s_j)

# Put data in list for STAN
stan_data <- list(N = N,
                  X = X,
                  D = D,
                  J = J,
                  delta = summ_stats$delta_i,
                  H = summ_stats$Hij,
                  a_j = a_j,
                  b_j = b_j, 
                  alpha = alpha, 
                  beta = beta
)


# Run STAN model
set.seed(model_seed)
fit <- stan(file = './TTE_pwc.stan',
            data = stan_data,
            chains = n_chains,
            warmup = burnin,
            iter = n_iters,
            seed = model_seed) %>%
  extract()

# Extract model output
samples = fit$lambda %>%
  cbind(fit$theta) %>%
  cbind(fit$tau2) %>%
  colnames_inplace(c(sprintf('lambda[%s]', 1:J), sprintf('theta[%s]', 1:D), 'tau^2'))


# Dose response curve (HR density for each dose)
samples %>%
  as_tibble() %>%
  select(contains('theta')) %>%
  add_column('theta[0]' = 0, .before = 'theta[1]') %>%
  pivot_longer(cols = everything(), names_to = 'TRTPN', values_to = 'HR',
               names_prefix = 'theta') %>%
  mutate(TRTPN = as.numeric(str_remove_all(TRTPN, "\\[|\\]")),
         HR = exp(HR),
         TRTPN = factor(key_arms[as.numeric(TRTPN) + 1], levels = key_arms)) %>%
  group_by(TRTPN) %>%
  summarize(low = quantile(HR, probs = 0.025),
            upp = quantile(HR, probs = 0.975),
            mean = mean(HR),
            median = median(HR)) %>%
  ggplot() +
  geom_point(aes(x = TRTPN, y = mean), pch = 10, size = 4, stroke = 1) +
  geom_line(aes(x = TRTPN, y = mean, group = 1), size = 1) +
  geom_line(aes(x = TRTPN, y = low, group = 1), lty = 3) +
  geom_line(aes(x = TRTPN, y = upp, group = 1), lty = 3) +
  geom_point(aes(x = 1:5, y = c(1, exp(theta_true))), col = 'red') +
  geom_ribbon(aes(x = TRTPN, ymin = low, ymax = upp, group = 1), alpha = 0.25) +
  labs(x = '', y = 'Hazard ratios')

# Define interval breaks (in weeks)
s_j = c(0, eps + seq(2.5, 20, by = 2.5))
J = length(s_j) - 1

# Prior hazard hyperparameters
a_j = rep(1, J)
b_j = rep(1/round(- log (1 - 0.4) / 20, 4), J)

alpha = 1
beta = 1

# Calculate summary statistics
summ_stats = get_summary_stats(y, nu, s_j)

# Put data in list for STAN
stan_data <- list(N = N,
                  X = X,
                  D = D,
                  J = J,
                  delta = summ_stats$delta_i,
                  H = summ_stats$Hij,
                  a_j = a_j,
                  b_j = b_j, 
                  alpha = alpha, 
                  beta = beta
)


# Run STAN model
set.seed(model_seed)
fit <- stan(file = './TTE_pwc.stan',
            data = stan_data,
            chains = n_chains,
            warmup = burnin,
            iter = n_iters,
            seed = model_seed) %>%
  extract()


# Extract model output
samples = fit$lambda %>%
  cbind(fit$theta) %>%
  cbind(fit$tau2) %>%
  colnames_inplace(c(sprintf('lambda[%s]', 1:J), sprintf('theta[%s]', 1:D), 'tau^2'))


# Survival curves for each treatment arm
post_surv_ctr = St(t = t_grid,
                   s_j = s_j,
                   lambda = samples[,1:J])
post_surv_trt1 = St(t = t_grid,
                    s_j = s_j,
                    lambda = exp(samples[,J+1]) * samples[,1:J])
post_surv_trt2 = St(t = t_grid,
                    s_j = s_j,
                    lambda = exp(samples[,J+2]) * samples[,1:J])
post_surv_trt3 = St(t = t_grid,
                    s_j = s_j,
                    lambda = exp(samples[,J+3]) * samples[,1:J])
post_surv_trt4 = St(t = t_grid,
                    s_j = s_j,
                    lambda = exp(samples[,J+4]) * samples[,1:J])
post_survival = post_surv_ctr %>%
  reshape2::melt(varnames = c('t', 'iteration')) %>%
  mutate(t = t_grid[t],
         TRTPN = key_arms[1]) %>%
  bind_rows(post_surv_trt1 %>%
              reshape2::melt(varnames = c('t', 'iteration')) %>%
              mutate(t = t_grid[t],
                     TRTPN = key_arms[2])) %>%
  bind_rows(post_surv_trt2 %>%
              reshape2::melt(varnames = c('t', 'iteration')) %>%
              mutate(t = t_grid[t],
                     TRTPN = key_arms[3])) %>%
  bind_rows(post_surv_trt3 %>%
              reshape2::melt(varnames = c('t', 'iteration')) %>%
              mutate(t = t_grid[t],
                     TRTPN = key_arms[4])) %>%
  bind_rows(post_surv_trt4 %>%
              reshape2::melt(varnames = c('t', 'iteration')) %>%
              mutate(t = t_grid[t],
                     TRTPN = key_arms[5])) %>%
  group_by(t, TRTPN) %>%
  summarize(mean = mean(value),
            low = quantile(value, probs = 0.025),
            upp = quantile(value, probs = 0.975)) %>%
  ungroup() %>%
  mutate(TRTPN = factor(TRTPN, levels = key_arms))


# KM Estimate
dat_KM = tibble(ExpTime = y, PRIM = nu, Dose = factor(key_arms[dose_label + 1]))
fit_KM <- survfit(Surv(ExpTime, PRIM) ~ Dose, data = dat_KM)
p_km = ggsurvplot(fit_KM,
                  data = dat_KM,
                  palette = key_palette,
                  legend = "top",
                  legend.title = "Arm",
                  legend.labs = levels(dat_KM$Dose),
                  axes.offset = FALSE,
                  break.time.by = 2.5,
                  ggtheme = custom_theme()
)

# Posterior survival curves and KM plot
p_km = p_km$plot +
  geom_line(data = post_survival, aes(x = t, y = mean, col = TRTPN), size = 1, lty = 4) +
  # geom_ribbon(data = post_survival,
  #             aes(x = t, ymin = low, ymax = upp, fill = TRTPN),
  #             alpha = 0.25, inherit.aes = FALSE) +
  geom_vline(xintercept = s_j, lty = 3) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  labs(x = 'Follow-up time (weeks)', y = 'Survival probability') +
  theme(legend.position = 'top')
p_km


# Dose response curve (HR density for each dose)
samples %>%
  as_tibble() %>%
  select(contains('theta')) %>%
  add_column('theta[0]' = 0, .before = 'theta[1]') %>%
  pivot_longer(cols = everything(), names_to = 'TRTPN', values_to = 'HR',
               names_prefix = 'theta') %>%
  mutate(TRTPN = as.numeric(str_remove_all(TRTPN, "\\[|\\]")),
         HR = exp(HR),
         TRTPN = factor(key_arms[as.numeric(TRTPN) + 1], levels = key_arms)) %>%
  group_by(TRTPN) %>%
  summarize(low = quantile(HR, probs = 0.025),
            upp = quantile(HR, probs = 0.975),
            mean = mean(HR),
            median = median(HR)) %>%
  ggplot() +
  geom_point(aes(x = TRTPN, y = mean), pch = 10, size = 4, stroke = 1) +
  geom_line(aes(x = TRTPN, y = mean, group = 1), size = 1) +
  geom_line(aes(x = TRTPN, y = low, group = 1), lty = 3) +
  geom_line(aes(x = TRTPN, y = upp, group = 1), lty = 3) +
  geom_point(aes(x = 1:5, y = c(1, exp(theta_true))), col = 'red') +
  geom_ribbon(aes(x = TRTPN, ymin = low, ymax = upp, group = 1), alpha = 0.25) +
  labs(x = '', y = 'Hazard ratios')


# Baseline hazard function
reshape2::melt(ht(t_grid, s_j, samples[,1:J]), varnames = c('Iteration', 'X')) %>%
  mutate(X = t_grid[X]) %>%
  group_by(X) %>%
  dplyr::summarize(ybar = mean(value),
                   UI = quantile(value, probs = 0.975),
                   LI = quantile(value, probs = 0.025)) %>%
  ggplot() +
  geom_line(aes(x = X, y = ybar), size = 1) +
  geom_ribbon(aes(x = X, ymin = LI, ymax = UI), alpha = 0.3) +
  geom_step(data = data.frame(s_j_true = s_j_true,
                              lambda_true = c(lambda_true, tail(lambda_true, 1))),
            aes(x = s_j_true, y = lambda_true), size = 1, col = 'red') +
  geom_vline(xintercept = s_j, lty = 3) +
  scale_x_continuous(name = 'time (weeks)', expand = c(0, 0)) + 
  labs(y = TeX("$\\lambda_{j}$"))
