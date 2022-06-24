

rm(list=ls())

source('./utils_plots.R')

weeks_to_months = 6/26
cols = brewer.pal(9, 'Set1')


data = read_csv(file = 'data_ordinal.csv') %>% 
  mutate(TRTPN = factor(TRTPN))
key_trt = levels(data$TRTPN)

digits = 0
n_digits = if_else(digits == 0, "%s", paste0("%0.", digits, "f"))
K = data %>% 
  filter(!is.na(RESP180)) %>% 
  select(RESP180) %>% 
  n_distinct



# Priors for the longitudinal model
long_prior = list(beta = rep(0, K), # priors for p30 
                  alpha = matrix(rep(0, K*K), 
                                 K, K, byrow = T))

imp_pars = get_imputation_pars(data, long_prior, key_trt)


imp_pars$p90_ctr %>%
  as_tibble() %>% 
  colnames_inplace(as.character(0:(K - 1))) %>% 
  mutate(RESP30 = as.character(0:(K - 1)), .before = 1) %>% 
  pivot_longer(
    cols = -RESP30,
    names_to = "RESP90",
    values_to = "transitions"
  ) %>% 
  group_by(RESP30) %>% 
  mutate(p_transitions = transitions / sum(transitions), 
         p_transitions = if_else(is.na(p_transitions), 0, p_transitions)) %>% 
  ungroup() %>% 
  mutate(label = sprintf(n_digits, transitions)) %>%
  plot_ly(
    x = ~RESP90, 
    y = ~RESP30,
    z = ~p_transitions, 
    text = ~label, 
    type = "heatmap", 
    hovertemplate = paste0('State <b>%{y}</b> -> State <b>%{x}</b>',
                           '<br><b>%{text}</b><extra></extra>'),
    colors = "Greens", 
    showscale = FALSE
  ) %>% 
  add_annotations(showarrow = FALSE, 
                  font = list(size = 16)) %>%
  layout(title = '', 
         xaxis = list(title = 'Response at 90 days', 
                      showgrid = FALSE, 
                      zeroline = FALSE), 
         yaxis = list(title = 'Response at 30 days', 
                      autorange = "reversed", 
                      # scaleanchor = "x", 
                      # scaleratio = 1, 
                      showgrid = FALSE, 
                      zeroline = FALSE)
  )





# Accrual profile
accr_profile = NULL
accr_profile$peak_rate = 12 * weeks_to_months # in weeks
accr_profile$ramp_up = 24 / weeks_to_months # in weeks


data %>% 
  mutate(TIMESINCERAND = as.numeric(difftime(RANDDT, min(RANDDT), units = 'weeks'))) %>% 
  arrange(TIMESINCERAND) %>% 
  mutate(Subject = row_number(), 
         EXPACCR = Lambda(t = TIMESINCERAND, peak_rate = accr_profile$peak_rate, 
                          ramp_up = accr_profile$ramp_up)) %>% 
  select(RANDDT, Subject, EXPACCR) %>% 
  plot_ly(
    type = 'scatter', 
    mode = 'lines', 
    line = list(color = cols), 
    hoverinfo = 'x+y') %>% 
  add_trace(x = ~RANDDT, 
            y = ~Subject, 
            name = 'observed', 
            line = list(shape = 'hv')) %>% 
  add_trace(x = ~RANDDT, 
            y = ~EXPACCR, 
            name = sprintf('%.2f per week', accr_profile$peak_rate)) %>%
  layout(title = 'Observed vs simulated accrual',
         xaxis = list(title = 'Randomization date'),
         yaxis = list(title = 'Cumulative number of subjects'), 
         legend = list(orientation = "h",   # show entries horizontally
                       xanchor = "center",  # use center of legend as anchor
                       x = 0.5, 
                       y = 1)
  )


data %>%
  filter(OTC180FL == 1) %>% 
  group_by(TRTPN, RESP180) %>% 
  summarise(n = n()) %>% 
  mutate(freq = n / sum(n)) %>% 
  right_join(list(
    TRTPN = key_trt,
    RESP180 = 0:(K-1)) %>%
      cross_df(), 
    by = c('TRTPN', 'RESP180')) %>%
  mutate(n = ifelse(is.na(n), 0, n), 
         freq = ifelse(is.na(freq), 0, freq)) %>% 
  arrange(TRTPN, RESP180) %>% 
  group_by(TRTPN) %>% 
  mutate(freq_sum = cumsum(freq)) %>% 
  plot_ly(x = ~freq, 
          y = ~TRTPN, 
          color = ~factor(RESP180, levels = 1:K-1),
          text = ~as.character(RESP180),
          colors = brewer.pal(K, 'RdBu')[K:1],
          type = 'bar', 
          hovertemplate = paste0('<b>%{y}</b><br>', 
                                 'Proportion with response = %{text}:<br>',
                                 '%{x:.2f}<extra></extra>')
  ) %>% 
  layout(title = 'Frequency of the response at 180 days',
         xaxis = list(title = 'Proportion of subjects'),
         yaxis = list(title = ''), 
         legend = list(title = list(text='<b>Response</b>'), 
                       traceorder = 'normal', 
                       orientation = "h",  # show entries horizontally
                       xanchor = "center", # use center of legend as anchor
                       x = 0.5, 
                       y = 1.02), 
         barmode = 'stack'
  )



data %>%
  filter(OTC180FL == 1) %>% 
  group_by(TRTPN, RESP180) %>% 
  summarise(n = n()) %>% 
  mutate(freq = n / sum(n)) %>% 
  right_join(list(
    TRTPN = key_trt,
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
  plot_ly(x = ~RESP180,
          y = ~freq, 
          color = ~TRTPN,
          text = ~as.character(TRTPN),
          colors = 'Dark2',
          type = 'bar',
          hovertemplate = paste0('<b>%{text}</b><br>',
                                 'Proportion with response = %{x}:<br>',
                                 '%{y:.2f}<extra></extra>')
  ) %>% 
  add_trace(x = ~imp_month_plot, 
            y = ~freq_sum, 
            color = ~TRTPN,
            text = ~as.character(TRTPN),
            type = 'scatter', 
            mode = 'lines', 
            line = list(shape = 'hv'), 
            showlegend = F, 
            hovertemplate = paste0('<b>%{text}</b><br>',
                                   'Cumulative proportion with response <= %{x}:<br>',
                                   '%{y:.2f}<extra></extra>')) %>% 
  layout(title = 'Frequency of the response at 180 days',
         xaxis = list(title = 'Proportion of subjects'),
         yaxis = list(title = ''), 
         legend = list(orientation = "h",  # show entries horizontally
                       xanchor = "center", # use center of legend as anchor
                       x = 0.5, 
                       y = 1), 
         barmode = 'group', 
         bargap = 0.3, 
         bargroupgap = 0.05
  )


# STAN data
stan_data <- list(N = nrow(data %>% filter(OTC180FL == 1)),
                  D = K, 
                  X = as.numeric(data %>% filter(OTC180FL == 1) %>% pull(TRTPN)) - 1, 
                  y = data %>% filter(OTC180FL == 1) %>% pull(RESP180) + 1
)


# Run STAN model
model_seed = 123
n_iters = 5000
burnin = 2500
n_thin = 2
n_chains = 4
samp_size = n_chains * (n_iters - burnin) / n_thin
set.seed(model_seed)
fit <- stan(file = '../Ordinal Regression/ord_logistic.stan',
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
  colnames_inplace(c('theta', sprintf('c%s', 1:(K - 1))))


samples %>% 
  as_tibble %>% 
  select(contains('c')) %>% 
  pivot_longer(cols = everything(), names_to = 'cutpoint') %>% 
  plot_ly(x = ~value, color = ~factor(cutpoint), 
          alpha = 0.6, type = "histogram", 
          histnorm = "probability") %>% 
  layout(barmode = "overlay", 
         xaxis = list(title = 'Cutpoints'), 
         yaxis = list(title = 'density'), 
         legend = list(orientation = "h",  # show entries horizontally
                       xanchor = "center", # use center of legend as anchor
                       x = 0.5, 
                       y = 1))



# Show the posterior distribution of the class probabilities under the control group
p_ctr_true = c(0.07, 0.2, 0.28, 0.2, 0.15, 0.1)
p_ctr_post = cbind(rep(0, samp_size), inv_logit(samples[,2:K]), rep(1, samp_size)) %>% 
  apply(1, diff) %>% 
  t() %>% 
  colnames_inplace(sprintf('p%s', 0:(K - 1)))
p_ctr_est = p_ctr_post %>% 
  as_tibble %>% 
  pivot_longer(cols = everything(), names_to = 'param') %>% 
  group_by(param) %>% 
  summarize(mean = mean(value), 
            low95 = quantile(value, 0.025), 
            upp95 = quantile(value, 0.975), 
            low50 = quantile(value, 0.25), 
            upp50 = quantile(value, 0.75)) %>% 
  mutate(true = p_ctr_true, 
         param = as.numeric(str_remove(param, 'p')))

p_ctr_est %>% 
  rbind(tail(p_ctr_est, 1) %>% mutate(param = param + 1)) %>% 
  plot_ly(x = ~param, y = ~upp95, 
          type = 'scatter', mode = 'lines',
          line = list(color = 'transparent', shape = 'hv'),
          showlegend = FALSE, name = '95% CI') %>% 
  add_trace(y = ~low95, 
            fill = 'tonexty', fillcolor='rgba(0,100,80,0.3)', 
            name = '95% CI') %>% 
  add_trace(y = ~upp50, 
            name = '50% CI') %>% 
  add_trace(y = ~low50, 
            fill = 'tonexty', fillcolor='rgba(0,100,80,0.4)', 
            name = '50% CI') %>% 
  add_trace(y = ~mean, 
            line = list(color='rgb(0,100,80)'),
            name = 'Mean') %>% 
  add_trace(y = ~true, 
            name = 'True', 
            line = list(color='black')) %>% 
  layout(title = "",
         xaxis = list(title = "", 
                      range = c(0, K),
                      ticktext = list("p0", "p1", "p2", "p3", "p4", "p5"), 
                      tickvals = list(0.5, 1.5, 2.5, 3.5, 4.5, 5.5)),
         yaxis = list(title = "", 
                      range = c(0, 0.5)))

