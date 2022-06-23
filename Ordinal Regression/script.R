

rm(list=ls())

source('./utils_OR.R')


cols = brewer.pal(9, 'Set1')

data_seed = 4321
weeks_to_months = 6/26
key_trt = c('Control', 'Treatment')


# We use the same setup of the design
OR_true = 2.2
p_ctr_true = c(0.07, 0.2, 0.28, 0.2, 0.15, 0.1)
p_trt_true = OR_transform(p_ctr_true, OR_true)


# Generate the data
set.seed(data_seed)
N = 300
dose_label = sample(c(rep(0, N/2), rep(1, N/2)))
TRTPN = factor(key_trt[dose_label + 1], levels = key_trt)
K = length(p_ctr_true)

# Accrual profile
accr_profile = NULL
accr_profile$peak_rate = 12 * weeks_to_months # in weeks
accr_profile$ramp_up = 24 / weeks_to_months # in weeks

# Transition probabilities
P_tr = matrix(c(1, 0, 0, 0, 0, 0, 
                0, 0.2, 0.4, 0.4, 0, 0, 
                0, 0, 0.34, 0.66, 0, 0, 
                0, 0, 0.22, 0.67, 0.11, 0, 
                0, 0, 0, 0.71, 0.29, 0, 
                0, 0, 0, 0, 0, 1), 
              K, K, byrow = T)


set.seed(data_seed)
data = sample_long_ordinal_data(N = N, 
                                TRTPN, 
                                p_ctr_true = p_ctr_true, 
                                p_trt_true = p_trt_true, 
                                P_tr, 
                                T_max = 0, 
                                accr_profile = accr_profile, 
                                t0 = 0)



# Set prior counts to 0 so that only the empirical transitions are 
# counted
pr_counts = list(beta = rep(0, K),
                 alpha = matrix(c(0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0),
                                K, K, byrow = T))

# Updated parameters for the longitudinal model
imp_pars = data %>%
  select(TRTPN,
         OTC30FL, RESP30,
         OTC90FL, RESP90,
         OTC180FL, RESP180) %>%
  get_imputation_pars(pr_counts, key_trt)

# Plot the number of transitions under control group
plt_ctr = plot_Ptrans(list(p30 = imp_pars$p30_ctr,
                           p90 = imp_pars$p90_ctr,
                           p180 = imp_pars$p180_ctr), 
                      digits = 0, 
                      palette = 'Greens')

# Plot the number of transitions under treatment group
plt_trt = plot_Ptrans(list(p30 = imp_pars$p30_trt,
                           p90 = imp_pars$p90_trt,
                           p180 = imp_pars$p180_trt), 
                      digits = 0)



# Descriptive statistics
data %>% 
  mutate(TIMESINCERAND = as.numeric(difftime(RANDDT, min(RANDDT), units = 'weeks'))) %>% 
  arrange(TIMESINCERAND) %>% 
  mutate(Subject = row_number(), 
         EXPACCR = Lambda(t = TIMESINCERAND, peak_rate = accr_profile$peak_rate, 
                          ramp_up = accr_profile$ramp_up)) %>% 
  select(TIMESINCERAND, Subject, EXPACCR) %>% 
  ggplot() + 
  geom_step(aes(x = TIMESINCERAND, y = Subject, col = 'observed'), size = 1) + 
  geom_line(aes(x = TIMESINCERAND, y = EXPACCR, 
                col = sprintf('%.2f per week', accr_profile$peak_rate)), 
            size = 1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'Weeks since start of accrual', 
       y = 'Cumulative number of subjects', 
       title = 'Observed vs simulated accrual') + 
  scale_color_brewer(name = '', palette = 'Dark2') + 
  theme(legend.position = "top")
# ggsave(paste0(dir_figures, 'accrual_plot.pdf'), width = 10, height = 8, units = 'in')


data %>%
  filter(OTC180FL == 1) %>% 
  group_by(TRTPN, RESP180) %>% 
  summarise(n = n()) %>% 
  mutate(freq = n / sum(n)) %>% 
  right_join(list(
    TRTPN = levels(TRTPN),
    RESP180 = 0:(K-1)) %>%
      cross_df(), 
    by = c('TRTPN', 'RESP180')) %>%
  mutate(TRTPN = factor(TRTPN), 
         n = ifelse(is.na(n), 0, n), 
         freq = ifelse(is.na(freq), 0, freq)) %>% 
  arrange(TRTPN, RESP180) %>% 
  group_by(TRTPN) %>% 
  mutate(freq_sum = cumsum(freq)) %>% 
  ggplot() +
  geom_bar(aes(x = freq, y = TRTPN, fill = factor(RESP180, levels = 0:(K-1))), position = position_stack(reverse = TRUE), stat = 'identity') + 
  labs(x = 'Proportion of subjects', 
       y = '', 
       title = 'Frequency of the response at 180 days') + 
  scale_fill_brewer(name = 'Response', palette = "RdBu", direction=-1) +
  theme(legend.position = "top") + 
  guides(fill = guide_legend(nrow = 1))
# ggsave(paste0(dir_figures, 'stacked_responses.pdf'), width = 10, height = 8, units = 'in')


data %>%
  filter(OTC180FL == 1) %>% 
  group_by(TRTPN, RESP180) %>% 
  summarise(n = n()) %>% 
  mutate(freq = n / sum(n)) %>% 
  right_join(list(
    TRTPN = levels(TRTPN),
    RESP180 = 0:(K-1)) %>%
      cross_df(), 
    by = c('TRTPN','RESP180')) %>%
  mutate(TRTPN = factor(TRTPN), 
         n = ifelse(is.na(n), 0, n), 
         freq = ifelse(is.na(freq), 0, freq)) %>% 
  arrange(TRTPN, RESP180) %>% 
  group_by(TRTPN) %>% 
  mutate(freq_sum = cumsum(freq), 
         imp_month_plot = ifelse(TRTPN == key_trt[1], RESP180, RESP180)) %>% 
  ggplot() + 
  geom_step(aes(x = imp_month_plot, y = freq_sum, col = TRTPN), size = 1) +
  geom_bar(aes(fill = TRTPN, col = TRTPN, y = freq, x = RESP180),
           width = 0.6, position = position_dodge(width = 0.7), stat = "identity") +
  labs(x = 'Proportion of subjects', 
       y = '', 
       title = 'Frequency of the response at 180 days') + 
  scale_x_continuous(breaks = 0:(K-1)) + 
  scale_colour_brewer(name = '', palette = 'Set2') + 
  scale_fill_brewer(name = '', palette = 'Set2') + 
  theme(legend.position = "top")
# ggsave(paste0(dir_figures, 'bars_responses.pdf'), width = 10, height = 8, units = 'in')



data %>%
  group_by(TRTPN) %>% 
  count(RESP30, RESP90, RESP180) %>% 
  mutate(id = row_number()) %>% 
  pivot_longer(cols = contains('RESP')) %>% 
  mutate(value = factor(value, levels = c(NA, as.character(0:(K-1))), exclude = NULL), 
         name = factor(name, levels = c('RESP30', 'RESP90', 'RESP180'))) %>% 
  group_by(name) %>% 
  ggplot(aes(x = name, y = n, stratum = value, fill = value, alluvium = id)) +
  ggalluvial::geom_stratum(alpha = 1) +
  ggalluvial::geom_flow() +
  facet_wrap(~ TRTPN) + 
  scale_fill_brewer(name = 'Response', breaks = c(NA, as.character(0:(K-1))), labels = c('Missing', as.character(0:(K-1))), palette = "RdBu", direction = -1, na.value = 'grey75') + 
  scale_x_discrete(name = '', labels = c('30 Days', '90 Days', '180 Days'), expand = c(0.1, 0.1)) +
  scale_y_continuous(name = 'Proportion of patients', expand = c(0.01, 0.01))


# Longitudinal model ------------------------------------------------------

model_seed = 123

# Priors for the longitudinal model
long_prior = list(beta = rep(1, K), # priors for p30 
                  alpha = matrix(rep(1/K, K*K), 
                                 K, K, byrow = T))

imp_pars = get_imputation_pars(data, long_prior, key_trt)


# Plot the posteriors under treatment and control 
plt_ctr = plot_Ptrans(list(p30 = imp_pars$p30_ctr, 
                           p90 = imp_pars$p90_ctr, 
                           p180 = imp_pars$p180_ctr), 
                      palette = 'Greens')

plt_trt = plot_Ptrans(list(p30 = imp_pars$p30_trt, 
                           p90 = imp_pars$p90_trt, 
                           p180 = imp_pars$p180_trt))


# Fit the ordinal regression model ----------------------------------------


# Now fit the frequentist version
fit_freq = MASS::polr(factor(RESP180) ~ TRTPN, data = data, Hess = TRUE)
fit_summ = summary(fit_freq)

log(OR_true); - fit_summ$coefficients[1,1]
logit(head(cumsum(p_ctr_true), -1)); as.numeric(fit_summ$coefficients[-1,1])

obs_OR = exp(- fit_summ$coefficients[1,1])
t_val = - fit_summ$coefficients[1,1] / fit_summ$coefficients[1,2]
dof = fit_freq$df.residual

pval = 1 - pt(q = t_val, df = dof)
pval



# Bayesian version of the model -------------------------------------------


# STAN data
stan_data <- list(N = nrow(data %>% filter(OTC180FL == 1)),
                  D = K, 
                  X = as.numeric(data %>% filter(OTC180FL == 1) %>% pull(TRTPN)) - 1, 
                  y = data %>% filter(OTC180FL == 1) %>% pull(RESP180) + 1
)

# Run STAN model
n_iters = 5000
burnin = 2500
n_thin = 2
n_chains = 4
samp_size = n_chains * (n_iters - burnin) / n_thin
set.seed(model_seed)
fit <- stan(file = 'ord_logistic.stan',
            data = stan_data,
            chains = n_chains,
            warmup = burnin,
            iter = n_iters,
            thin = n_thin,
            seed = model_seed) %>%
  extract()

# Extract model output
samples = fit$theta %>%
  cbind(fit$c) %>%
  colnames_inplace(c('theta', sprintf('c_%s', 1:(K - 1))))


log(OR_true); mean(samples[,1])
logit(head(cumsum(p_ctr_true), -1)); as.numeric(colMeans(samples[,2:K]))

# Show the posterior distributions of the cutpoints 
samples %>% 
  as_tibble %>% 
  select(contains('c_')) %>% 
  pivot_longer(cols = everything(), names_to = 'cutpoint') %>% 
  ggplot() + 
  geom_histogram(aes(x = value, y = ..density.., fill = cutpoint), 
                 col = 'white', alpha = 0.4) +
  scale_x_continuous(name = TeX('$c_{j}$')) + 
  scale_fill_brewer(name = '', palette = 'Dark2', labels = sprintf('j = %i', 1:(K-1))) + 
  theme(legend.position = 'top') + 
  guides(fill = guide_legend(nrow = 1))


# Show the posterior distribution of the class probabilities under the control group
p_ctr_post = cbind(rep(0, samp_size), inv_logit(samples[,2:K]), rep(1, samp_size)) %>% 
  apply(1, diff) %>% 
  t() %>% 
  colnames_inplace(sprintf('p_%s', 0:(K - 1)))
p_ctr_est = p_ctr_post %>% 
  as_tibble %>% 
  pivot_longer(cols = everything(), names_to = 'param') %>% 
  group_by(param) %>% 
  summarize(mean = mean(value), 
            low = quantile(value, 0.025), 
            upp = quantile(value, 0.975)) %>% 
  mutate(true = p_ctr_true, 
         param = as.numeric(str_remove(param, 'p_')))

p_ctr_est %>% 
  rbind(tail(p_ctr_est, 1) %>% mutate(param = param + 1)) %>% 
  ggplot() + 
  geom_step(aes(x = param - 0.5, y = mean), col = brewer.pal(8, 'Dark2')[1], size = 1) + 
  pammtools::geom_stepribbon(aes(x = param - 0.5, ymin = low, ymax = upp), fill = brewer.pal(8, 'Dark2')[1], alpha = 0.3) + 
  geom_step(aes(x = param - 0.5, y = true)) + 
  scale_y_continuous(name = '', expand = c(0, 0), limits = c(0, max(p_ctr_est$upp) + 0.1)) + 
  scale_x_continuous(name = '', breaks = 0:(K-1), labels = TeX(sprintf("$p_{%s}$", 0:(K-1))), expand = c(0, 0))


# Show the posterior distribution of the class probabilities under the treatment group
p_trt_post = cbind(rep(0, samp_size), inv_logit(samples[,2:K] + samples[,1]), rep(1, samp_size)) %>% 
  apply(1, diff) %>% 
  t() %>% 
  colnames_inplace(sprintf('p_%s', 0:(K - 1)))
p_trt_est = p_trt_post %>% 
  as_tibble %>% 
  pivot_longer(cols = everything(), names_to = 'param') %>% 
  group_by(param) %>% 
  summarize(mean = mean(value), 
            low = quantile(value, 0.025), 
            upp = quantile(value, 0.975)) %>% 
  mutate(true = p_trt_true, 
         param = as.numeric(str_remove(param, 'p_')))

p_trt_est %>% 
  rbind(tail(p_trt_est, 1) %>% mutate(param = param + 1)) %>% 
  ggplot() + 
  geom_step(aes(x = param - 0.5, y = mean), col = brewer.pal(8, 'Dark2')[2], size = 1) + 
  pammtools::geom_stepribbon(aes(x = param - 0.5, ymin = low, ymax = upp), fill = brewer.pal(8, 'Dark2')[2], alpha = 0.3) + 
  geom_step(aes(x = param - 0.5, y = true)) + 
  scale_y_continuous(name = '', expand = c(0, 0), limits = c(0, max(p_trt_est$upp) + 0.1)) + 
  scale_x_continuous(name = '', breaks = 0:(K-1), labels = TeX(sprintf("$p_{%s}$", 0:(K-1))), expand = c(0, 0))



# Calculate predictive probabilities --------------------------------------

fit_pp = fit_model(df = data, long_prior = long_prior, n_max = 500, 
                   key_trt = key_trt, iters = 1000)

cat("The probability of success at the current sample size is", round(mean(fit_pp$PPn), 2), "\n")
cat("The probability of success at the maxiumum sample size is", round(mean(fit_pp$PPmax), 2), "\n")

