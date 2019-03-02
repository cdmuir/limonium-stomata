source("r/header.R")

stomata <- read_csv("data/limonium-stomata.csv") %>%
  group_by(spp, treatment, rep) %>%
  summarize(
    sd = sd[surface == "aB"] + sd[surface == "aD"],
    sp = (sp[surface == "aB"] + sp[surface == "aD"]) / 2
  ) %>% 
  ungroup()

rubisco <- read_csv("data/limonium-kinetics.csv")

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

stomata %<>% mutate_if(function(.x) inherits(.x, "units"), drop_units)

# Prepare for Stan ----
m1_stan <- stomata %>%
  select(spp_stomata = spp, treatment, sd1 = sd, sp1 = sp, b, m) %>%
  mutate(drought = as.numeric(treatment == "WD")) %>%
  filter(!is.na(sd1), !is.na(sp1)) %>%
  as.list()

m1_stan$spp_stomata %<>% 
  as.factor() %>% 
  as.integer()

m1_stan$n_stomata <- length(m1_stan$spp_stomata)
m1_stan$n_spp <- max(m1_stan$spp_stomata)

m1_stan$b %<>% first()
m1_stan$m %<>% first()

m1_stan$kcatc <- filter(rubisco, trait == "kcatc") %>%
  pull(value)

m1_stan$spp_kcatc <- filter(rubisco, trait == "kcatc") %>%
  pull(spp) %>%
  as.factor() %>%
  as.numeric()

m1_stan$n_kcatc <- length(m1_stan$spp_kcatc)

m1_stan$sco <- filter(rubisco, trait == "sco") %>%
  pull(value)

m1_stan$spp_sco <- filter(rubisco, trait == "sco") %>%
  pull(spp) %>%
  as.factor() %>%
  as.numeric()

m1_stan$n_sco <- length(m1_stan$spp_sco)

write_rds(m1_stan, "objects/m1_stan.rds")

# m1 <- stan("m1.stan", data = m1_stan, cores = 4, iter = 2e4, thin = 1e1)
# write_rds(m1, "objects/m1.rds")

m1 <- read_rds("objects/m1.rds")

# Parameter summary ----

pars1 <- c("b0_sd", "b_sd_drought", "b0_sp", "b_sp_drought", 
           "sigma_sd", "sigma_sp", "cor_1", "rescor", "Omega_rubisco[1,2]")

m1_tab <- tidy(m1, pars = pars1, rhat = TRUE, ess = TRUE, conf.int = TRUE) %>%
  left_join(pars1 %>%
              as.data.frame(m1, pars = .) %>%
              summarize_all(get_p) %>%
              gather(term, p), by = "term")

write_rds(m1_tab, "objects/m1_tab.rds")

# Figure 2 ----

fig2 <- stomata %>%
  filter(!is.na(sd), !is.na(s)) %>%
  select(spp, treatment, sd, s, gsmax) %>%
  group_by(spp, treatment) %>%
  summarize(
    sd_sd = sd(sd),
    sd_s = sd(s),
    sd_gsmax = sd(gsmax),
    sd = mean(sd),
    s = mean(s),
    gsmax = mean(gsmax),
    n = n()
  ) %>%
  mutate(
    se_sd = sd_sd / sqrt(n),
    se_s = sd_s / sqrt(n),
    se_gsmax = sd_gsmax / sqrt(n)
  ) %>%
  mutate_if(function(.x) inherits(.x, "units"), drop_units)

## Panel A ----

f2a <- ggplot(fig2, aes(
  x = sd, y = s,
  xmin = sd - se_sd, xmax = sd + se_sd,
  ymin = s - se_s, ymax = s + se_s,
  fill = treatment)
) +
  scale_x_continuous(limits = c(70, 180), breaks = seq(75, 175, 25),
                     position = "top") +
  scale_y_continuous(limits = c(900, 2100), breaks = c(1000, 1500, 2000)) +
  scale_fill_manual("Treatment", values = c("black", "white")) +
  geom_errorbar() +
  geom_errorbarh() +
  geom_point(size = 3, shape = 21) +
  xlab(expression(paste("Stomatal density [m", m^-2, "]"))) +
  ylab(expression(paste("Stomatal size [", mu, m^2, "]"))) +
  theme_cdm() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(color = rgb(0, 0, 0, alpha = 0)),
    plot.margin = unit(c(0, 0, 0, 0), units = "cm")
  ) +
  NULL

## Panel B ----

f2b <- ggplot(fig2, aes(
  x = sd, y = gsmax,
  xmin = sd - se_sd, xmax = sd + se_sd,
  ymin = gsmax - se_gsmax, ymax = gsmax + se_gsmax,
  fill = treatment)
) +
  scale_x_continuous(limits = c(70, 180), breaks = seq(75, 175, 25)) +
  scale_y_continuous(limits = c(1.2, 2.8), breaks = seq(1.25, 2.75, 0.5)) +
  scale_fill_manual(values = c("black", "white")) +
  geom_errorbar() +
  geom_errorbarh() +
  geom_point(size = 3, shape = 21) +
  xlab(expression(paste("Stomatal density [m", m^-2, "]"))) +
  ylab(expression(paste("Anatomical ", italic(g)[smax], " [mol ", H[2], "O ", m^-2~s^-1, "]"))) +
  theme_cdm() +
  theme(
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), units = "cm")
  ) +
  NULL

## Panel C ----

f2c <- ggplot(fig2, aes(
  x = s, y = gsmax,
  xmin = s - se_s, xmax = s + se_s,
  ymin = gsmax - se_gsmax, ymax = gsmax + se_gsmax,
  fill = treatment)
) +  
  scale_x_continuous(limits = c(900, 2100), breaks = c(1000, 1500, 2000)) +
  scale_y_continuous(limits = c(1.2, 2.8), breaks = seq(1.25, 2.75, 0.5),
                     position = "right") +
  scale_fill_manual(values = c("black", "white")) +
  geom_errorbar() +
  geom_errorbarh() +
  geom_point(size = 3, shape = 21) +
  xlab(expression(paste("Stomatal size [", mu, m^2, "]"))) +
  ylab(expression(paste("Anatomical ", italic(g)[smax], " [mol ", H[2], "O ", m^-2~s^-1, "]"))) +
  theme_cdm() +
  theme(
    axis.title.y = element_text(color = rgb(0, 0, 0, alpha = 0)),
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), units = "cm")
  ) +
  NULL

fguide <- get_legend(f2a + theme(legend.position = "left"))

left_col <- plot_grid(f2a, f2b, labels = c("A", "B"), ncol = 1, align = "v", 
                      axis = "r", hjust = -5, vjust = c(5, 1.5))
rght_col <- plot_grid(fguide, f2c, labels = c("", "C"), ncol = 1, align = "hv", 
                      axis = "l")#, label_x = 1, hjust = 1.5)
plot_grid(left_col, rght_col, ncol = 2, align = "hv")

ggsave("figures/fig2.eps", w = 6.5, h = 6.175)

# Figure 3 ----

ex <- m1 %>%
  as.data.frame(
    pars = c("gsmax_spp", "kcatc_spp", "sco_spp")
  ) %>%
  mutate(iter = 1:nrow(.))

fig3 <- m1 %>%
  tidy(pars = c("gsmax_spp", "kcatc_spp", "sco_spp")) %>%
  mutate(
    spp1 = str_replace(term, "^[a-z]+_spp\\[([0-9]+)\\]$", "\\1"),
    trait = str_replace(term, "^*([a-z]+)_spp\\[[0-9]+\\]$", "\\1")
  ) %>%
  select(-term) %>%
  rename(mu = estimate, se = std.error) %>%
  gather(key, value, -trait, -spp1) %>%
  mutate(term = str_c(key, "_", trait)) %>%
  select(term, spp1, value) %>%
  gather(key, value, -spp1, -term) %>%
  select(-key) %>%
  spread(term, value)

b1 <- ex %>%
  gather(par, value, -iter) %>%
  mutate(
    spp = str_replace(par, "^[[:print:]]+_spp\\[([0-9]+)\\]$", "\\1"),
    trait = str_replace(par, "^([[:print:]]+)_spp\\[([0-9]+)\\]$", "\\1")
  ) %>%
  select(-par) %>%
  group_by(iter) %>%
  spread(trait, value) %>%
  summarise(
    b0_kcatc = coef(lm(kcatc ~ gsmax))[1],
    b_kcatc_gsmax = coef(lm(kcatc ~ gsmax))[2],
    b0_sco = coef(lm(sco ~ gsmax))[1],
    b_sco_gsmax = coef(lm(sco ~ gsmax))[2]
  ) %>%
  select(-iter)

m1_tab %<>% 
  bind_rows(
    b1 %>%
      tidyMCMC(estimate.method = "median", conf.int = TRUE, 
               conf.method = "HPDinterval") %>%
      left_join(b1 %>%
                  summarize_all(get_p) %>%
                  gather(term, p), by = "term"))

b1 %<>% 
  crossing(x = seq(min(fig3$mu_gsmax), 
                   max(fig3$mu_gsmax), length.out = 1e2)) %>%
  mutate(
    mu_kcatc = b0_kcatc + b_kcatc_gsmax * x,
    mu_sco = b0_sco + b_sco_gsmax * x
  ) %>%
  group_by(x) %>%
  summarize(
    est_kcatc = median(mu_kcatc),
    lci_kcatc = hdi(mu_kcatc)[1, 1],
    uci_kcatc = hdi(mu_kcatc)[1, 2],
    est_sco = median(mu_sco),
    lci_sco = hdi(mu_sco)[1, 1],
    uci_sco = hdi(mu_sco)[1, 2]
  )

r2 <- as.data.frame(m1, c("gsmax_spp", "kcatc_spp", "sco_spp")) %>%
  mutate(iter = 1:nrow(.)) %>%
  gather(par, value, -iter) %>%
  mutate(
    spp = str_replace(par, "^[[:print:]]+_spp\\[([0-9]+)\\]$", "\\1"),
    trait = str_replace(par, "^([[:print:]]+)_spp\\[([0-9]+)\\]$", "\\1")
  ) %>%
  select(-par) %>%
  group_by(iter) %>%
  spread(trait, value) %>%
  summarize(
    r_kcatc = cor(gsmax, kcatc),
    r_sco = cor(gsmax, sco),
    r2_kcatc = cor(gsmax, kcatc) ^ 2,
    r2_sco = cor(gsmax, sco) ^ 2
  ) %>%
  ungroup() %>%
  select(-iter) %>%
  tidyMCMC(estimate.method = "median", conf.int = TRUE, conf.method = "HPDinterval")

## Panel A ----

f3a <- ggplot(fig3, aes(
  x = mu_gsmax, y = mu_kcatc,
  xmin = mu_gsmax - se_gsmax, xmax = mu_gsmax + se_gsmax,
  ymin = mu_kcatc - se_kcatc, ymax = mu_kcatc + se_kcatc
)) +  
  scale_x_continuous(limits = c(1.2, 2.8), breaks = seq(1.25, 2.75, 0.5),
                     position = "top") +
  geom_ribbon(data = b1, 
              mapping = aes(x = x, ymin = lci_kcatc, ymax = uci_kcatc), 
              inherit.aes = FALSE, fill = "grey50") +
  geom_line(data = b1, mapping = aes(x = x, y = est_kcatc), 
            size = 1.2, inherit.aes = FALSE) +
  geom_errorbar() +
  geom_errorbarh() +
  geom_point(size = 3) +
  annotate(
    "text", 1.25, 3.75, hjust = 0, vjust = 1, parse = TRUE,
    label = glue("atop(paste(italic(R) ^ 2, \" = {r2}\"), paste(italic(P), \"{p}\"))",
                 r2 = round(r2$estimate[r2$term == "r2_kcatc"], 2), 
                 p = ifelse(
                   m1_tab$p[m1_tab$term == "b_kcatc_gsmax"] < 0.001,
                   " < 0.001",
                   str_c(" = ", round(m1_tab$p[m1_tab$term == "b_kcatc_gsmax"], 3)
                 )))
    ) +
      ylab(expression(paste(italic(k)[textstyle(cat)]^textstyle(c), " [", s ^ -1, "]"))) +
  xlab(expression(paste("Anatomical ", italic(g)[smax], " [mol ", H[2], "O ", m^-2~s^-1, "]"))) +
  theme_cdm() +
  theme(
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), units = "cm"),
    axis.title.x = element_text(color = rgb(0, 0, 0, alpha = 0))
  ) +
  NULL

## Panel B ----

f3b <- ggplot(fig3, aes(
  x = mu_gsmax, y = mu_sco,
  xmin = mu_gsmax - se_gsmax, xmax = mu_gsmax + se_gsmax,
  ymin = mu_sco - se_sco, ymax = mu_sco + se_sco
)) +  
  scale_x_continuous(limits = c(1.2, 2.8), breaks = seq(1.25, 2.75, 0.5)) +
  geom_ribbon(data = b1, 
              mapping = aes(x = x, ymin = lci_sco, ymax = uci_sco), 
              inherit.aes = FALSE, fill = "gray50") +
  geom_line(data = b1, mapping = aes(x = x, y = est_sco), 
            size = 1.2, inherit.aes = FALSE) +
  geom_errorbar() +
  geom_errorbarh() +
  geom_point(size = 3) +
  annotate(
    "text", 2.75, 120, hjust = 1, vjust = 1, parse = TRUE,
    label = glue("atop(paste(italic(R) ^ 2, \" = {r2}\"), paste(italic(P), \" = {p}\"))",
                 r2 = round(r2$estimate[r2$term == "r2_sco"], 2), 
                 p = round(m1_tab$p[m1_tab$term == "b_sco_gsmax"], 3))
  ) +
  ylab(expression(paste(italic(S)[C/O], " [mol ", mol ^ -1, "]"))) +
  xlab(expression(paste("Anatomical ", italic(g)[smax], " [mol ", H[2], "O ", m^-2~s^-1, "]"))) +
  theme_cdm() +
  theme(
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), units = "cm")
  ) +
  NULL

plot_grid(f3a, f3b, ncol = 1, align = "hv", labels = c("A", "B"), axis = "t")

ggsave("figures/fig3.eps", w = 3.25, h = 6.5)

