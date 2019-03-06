source("r/header.R")

# Prepare data ----

stomata <- read_csv("data/limonium-stomata.csv") %>%
  group_by(spp, treatment, rep) %>%
  summarize(
    sd = sd[surface == "aB"] + sd[surface == "aD"],
    sp = (sp[surface == "aB"] + sp[surface == "aD"]) / 2
  ) %>%
  ungroup()

lags <- read_csv("data/limonium-lags.csv")

stomata %<>% full_join(select(lags, spp, treatment, rep, gs), 
                       by = c("spp", "treatment", "rep"))

# Calculate gsmax ----

# Based on Monteith and Unsworth Table A.3 (pg. 376)
D_wv <- set_units(24.9e-6, m ^ 2 / s)
v <- set_units(2.24e-2, m ^ 3/ mol)
stomata$b <- biophysical_constant(D_wv, v)
stomata$m <- set_units(morphological_constant(0.5, 0.5, 0.5))
stomata$sd %<>% set_units(1 / mm ^ 2)
stomata$sp %<>% set_units(um)
stomata %<>% 
  mutate(
    s = 0.5 * (2 * sp) ^ 2,
    gsmax = b * m * sd * sqrt(s)
  )
stomata$gsmax %<>% set_units(mol/m^2/s)

stomata$treatment %<>% factor(levels = c("WW", "WD"))

# Fit models ----

# m2a <- brm(log(gs) ~ gsmax * treatment + (1 | spp), 
#            data = mutate_if(stomata, ~ inherits(.x, "units"), drop_units),
#            cores = 4)
# m2b <- brm(log(gs) ~ treatment + gsmax + (1 | spp), 
#            data = mutate_if(stomata, ~ inherits(.x, "units"), drop_units),
#            cores = 4)
# write_rds(m2a, "objects/m2a.rds")
# write_rds(m2b, "objects/m2b.rds")

m2a <- read_rds("objects/m2a.rds")
m2b <- read_rds("objects/m2b.rds")

# Parameter summary ----

pars2 <- c("b_Intercept", "b_treatmentWD", "b_gsmax")

m2_tab <- m2b$fit %>%
  tidy(pars = pars2, rhat = TRUE, ess = TRUE, conf.int = TRUE) %>%
  left_join(pars2 %>%
              as.data.frame(m2b$fit, pars = .) %>%
              summarize_all(get_p) %>%
              gather(term, p), by = "term")

write_rds(m2_tab, "objects/m2_tab.rds")


# Figure 4 ----

r2 <- bayes_R2(m2b)

fig4 <- stomata %>%
  filter(!is.na(gs), !is.na(gsmax)) %>%
  select(spp, Treatment = treatment, gs, gsmax) %>%
  mutate_if(function(.x) inherits(.x, "units"), drop_units)

ex <- m2b %>%
  as.data.frame(pars = c("b_Intercept", "b_treatmentWD", "b_gsmax")) %>%
  mutate(iter = 1:nrow(.)) %>%
  crossing(x = seq(min(fig4$gsmax), max(fig4$gsmax), length.out = 1e2),
           drought = c(0, 1)) %>%
  mutate(
    mu_gs = exp(b_Intercept + b_gsmax * x + b_treatmentWD * drought)
  ) %>%
  group_by(x, drought) %>%
  summarize(
    est_gs = median(mu_gs),
    lci_gs = hdi(mu_gs)[1, 1],
    uci_gs = hdi(mu_gs)[1, 2]
  ) %>%
  mutate(Treatment = ifelse(drought == 0, "WW", "WD"))

f4 <- ggplot(fig4, aes(x = gsmax, y = gs, fill = Treatment)) +  
  scale_x_continuous(limits = c(1, 3)) +
  scale_y_log10() +
  scale_fill_manual(values = c("black", "white")) +
  annotate(
    "text", 3, 0.05, hjust = 1, vjust = 0, parse = TRUE,
    label = glue("atop(paste(italic(R) ^ 2, \" = {r2}\"), paste(italic(P), \"{p}\"))",
                 r2 = round(r2["R2", "Estimate"], 2), 
                 p = ifelse(
                   m2_tab$p[m2_tab$term == "b_gsmax"] < 0.001,
                   " < 0.001",
                   str_c(" = ", round(m2_tab$p[m2_tab$term == "b_gsmax"], 3)
                   )))
  ) +
  geom_ribbon(data = filter(ex, Treatment == "WW"), 
              mapping = aes(x = x, ymin = lci_gs, ymax = uci_gs), 
              inherit.aes = FALSE, fill = "grey50", color = "grey25") +
  geom_ribbon(data = filter(ex, Treatment == "WD"), 
              mapping = aes(x = x, ymin = lci_gs, ymax = uci_gs), 
              inherit.aes = FALSE, fill = "grey50", color = "grey25") +
  geom_line(data = filter(ex, Treatment == "WW"), 
            mapping = aes(x = x, y = est_gs), 
            inherit.aes = FALSE, size = 1.2) +
  geom_line(data = filter(ex, Treatment == "WD"), 
            mapping = aes(x = x, y = est_gs), 
            size = 1.2, inherit.aes = FALSE, linetype = "dashed") +
  geom_point(size = 2, shape = 21) +
  xlab(expression(paste("Anatomical ", italic(g)[smax], " [mol ", H[2], "O ", m^-2~s^-1, "]"))) +
  ylab(expression(paste(italic(g)[s], " [mol ", H[2], "O ", m^-2~s^-1, "]"))) +
  theme_cdm() +
  theme(
    legend.position = "top",
    plot.margin = unit(c(0, 0, 0, 0), units = "cm")
  ) +
  NULL

ggsave("figures/fig4.eps", width = 4, height = 4)
