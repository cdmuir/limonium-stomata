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

# There is evidence for species x drought interactions on stomatal traits
m1_df <- stomata %>%
  select(spp_stomata = spp, treatment, sd1 = sd, sp1 = sp, b, m) %>%
  mutate(drought = as.numeric(treatment == "WD")) %>%
  filter(!is.na(sd1), !is.na(sp1))

if (run) {
  
  ma <- brm(mvbind(sd1, sp1) ~ drought + (1 | spp_stomata), data = m1_df,
            cores = 4, chains = 4, seed = 55902492) %>%
    add_criterion("loo", reloo = TRUE)
  mb <- brm(mvbind(sd1, sp1) ~ drought + (drought | spp_stomata), data = m1_df,
            cores = 4, chains = 4, seed = 167740183) %>%
    add_criterion("loo", reloo = TRUE)
  
  write_rds(ma, "objects/ma_gsmaxdiff.rds")
  write_rds(mb, "objects/mb_gsmaxdiff.rds")
  
} else {
  
  ma <- read_rds("objects/ma_gsmaxdiff.rds")
  mb <- read_rds("objects/mb_gsmaxdiff.rds")
  
}

bayes_R2(ma)
bayes_R2(mb)
loo_compare(ma, mb, criterion = "loo")

if (run) {
  
  # m1a: no species x drought interaction
  m1a <- stan("m1a.stan", data = m1_stan, cores = 4, iter = 2e4, thin = 1e1,
              seed = 353601735)
  write_rds(m1a, "objects/m1a.rds")

  # m1b: including species x drought interaction
  m1b <- stan("m1b.stan", data = m1_stan, cores = 4, iter = 2e4, thin = 1e1,
              seed = 605663820)
  write_rds(m1b, "objects/m1b.rds")

} else {
  
  m1a <- read_rds("objects/m1a.rds")
  m1b <- read_rds("objects/m1b.rds")

}

# Parameter summary ----

pars1 <- c("b0_sd", "b_sd_drought", "b0_sp", "b_sp_drought", 
           "sigma_sd", "sigma_sp", "cor_1", "cor_2", "rescor", 
           "Omega_rubisco[1,2]")

m1_tab <- tidy(m1b, pars = pars1, rhat = TRUE, ess = TRUE, conf.int = TRUE) %>%
  left_join(pars1 %>%
              as.data.frame(m1b, pars = .) %>%
              summarize_all(get_p) %>%
              gather(term, p), by = "term")

write_rds(m1_tab, "objects/m1_tab.rds")

# Figure S2 ----

figS2 <- stomata %>%
  filter(!is.na(sd), !is.na(sp)) %>%
  select(spp, treatment, sd, sp, gsmax) %>%
  group_by(spp, treatment) %>%
  summarize(
    sd_sd = sd(sd),
    sd_sp = sd(sp),
    sd_gsmax = sd(gsmax),
    sd = mean(sd),
    sp = mean(sp),
    gsmax = mean(gsmax),
    n = n()
  ) %>%
  mutate(
    se_sd = sd_sd / sqrt(n),
    se_sp = sd_sp / sqrt(n),
    se_gsmax = sd_gsmax / sqrt(n)
  ) %>%
  mutate_if(function(.x) inherits(.x, "units"), drop_units)

## Panel A ----

fS2a <- ggplot(figS2, aes(
  x = sd, y = sp,
  xmin = sd - se_sd, xmax = sd + se_sd,
  ymin = sp - se_sp, ymax = sp + se_sp,
  fill = treatment)
) +
  scale_x_continuous(limits = c(70, 180), breaks = seq(75, 175, 25),
                     position = "top") +
  scale_y_continuous(limits = c(20, 32), breaks = seq(21, 30, 3)) +
  scale_fill_manual("Treatment", values = c("black", "white")) +
  geom_errorbar() +
  geom_errorbarh() +
  geom_point(size = 3, shape = 21) +
  xlab(expression(paste("Stomatal density [m", m^-2, "]"))) +
  ylab(expression(paste("Stomatal pore lengt", h^phantom(1), " [", mu, "m]"))) +
  theme_cdm() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(color = rgb(0, 0, 0, alpha = 0)),
    plot.margin = unit(c(0, 0, 0, 0), units = "cm")
  ) +
  NULL

## Panel B ----

fS2b <- ggplot(figS2, aes(
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

fS2c <- ggplot(figS2, aes(
  x = sp, y = gsmax,
  xmin = sp - se_sp, xmax = sp + se_sp,
  ymin = gsmax - se_gsmax, ymax = gsmax + se_gsmax,
  fill = treatment)
) +  
  scale_x_continuous(limits = c(20, 32), breaks = seq(21, 30, 3)) +
  scale_y_continuous(limits = c(1.2, 2.8), breaks = seq(1.25, 2.75, 0.5),
                     position = "right") +
  scale_fill_manual(values = c("black", "white")) +
  geom_errorbar() +
  geom_errorbarh() +
  geom_point(size = 3, shape = 21) +
  xlab(expression(paste("Stomatal pore lengt", h^phantom(1), " [", mu, "m]"))) +
  ylab(expression(paste("Anatomical ", italic(g)[smax], " [mol ", H[2], "O ", m^-2~s^-1, "]"))) +
  theme_cdm() +
  theme(
    axis.title.y = element_text(color = rgb(0, 0, 0, alpha = 0)),
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), units = "cm")
  ) +
  NULL

fguide <- get_legend(fS2a + theme(legend.position = "left"))

left_col <- plot_grid(fS2a, fS2b, labels = c("A", "B"), ncol = 1, align = "v", 
                      axis = "r", hjust = -5, vjust = c(5, 1.5))
rght_col <- plot_grid(fguide, fS2c, labels = c("", "C"), ncol = 1, 
                      align = "hv", axis = "l")
plot_grid(left_col, rght_col, ncol = 2, align = "hv")

ggsave("figures/figS2.pdf", w = 6.5, h = 6.175)

# Figure 2 ----

ex <- m1b %>%
  as.data.frame(
    pars = c("gsmax_ww_spp", "gsmax_wd_spp", "kcatc_spp", "sco_spp")
  ) %>%
  mutate(iter = 1:nrow(.))

fig2 <- m1b %>%
  tidy(pars = c("gsmax_ww_spp", "gsmax_wd_spp", "kcatc_spp", "sco_spp")) %>%
  mutate(
    spp1 = str_replace(term, "^[a-z]+_[a-z]*_*spp\\[([0-9]+)\\]$", "\\1"),
    trait = str_replace(term, "^*([a-z]+)_[a-z]*_*spp\\[[0-9]+\\]$", "\\1"),
    treatment = str_replace(term, "^[a-z]+_([a-z]*)_*spp\\[[0-9]+\\]$", "\\1")
  ) %>%
  select(-term) %>%
  rename(mu = estimate, se = std.error) %>%
  gather(key, value, -trait, -treatment, -spp1) %>%
  mutate(treatment = str_replace(treatment, "^$", "ww")) %>%
  bind_rows(
    filter(. , trait != "gsmax") %>%
      mutate(treatment = "wd")
  ) %>%
  group_by(spp1, trait, treatment) %>%
  spread(key, value)

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
    b0_kcatc_wd = coef(lm(kcatc ~ gsmax_wd))[1],
    b_kcatc_gsmax_wd = coef(lm(kcatc ~ gsmax_wd))[2],
    b0_sco_wd = coef(lm(sco ~ gsmax_wd))[1],
    b_sco_gsmax_wd = coef(lm(sco ~ gsmax_wd))[2],
    b0_kcatc_ww = coef(lm(kcatc ~ gsmax_ww))[1],
    b_kcatc_gsmax_ww = coef(lm(kcatc ~ gsmax_ww))[2],
    b0_sco_ww = coef(lm(sco ~ gsmax_ww))[1],
    b_sco_gsmax_ww = coef(lm(sco ~ gsmax_ww))[2]
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
  crossing(x = seq(min(filter(fig2, trait == "gsmax") %>%
                         transmute(lower = mu - se) %>%
                         pull(lower)),
                   max(filter(fig2, trait == "gsmax") %>%
                         transmute(upper = mu + se) %>%
                         pull(upper)),
                   length.out = 1e2)) %>%
  mutate(
    mu_kcatc_wd = b0_kcatc_wd + b_kcatc_gsmax_wd * x,
    mu_sco_wd = b0_sco_wd + b_sco_gsmax_wd * x,
    mu_kcatc_ww = b0_kcatc_ww + b_kcatc_gsmax_ww * x,
    mu_sco_ww = b0_sco_ww + b_sco_gsmax_ww * x
  ) %>%
  select_at(vars(x, mu_kcatc_wd, mu_sco_wd, mu_kcatc_ww, mu_sco_ww)) %>%
  group_by(x) %>%
  point_interval(.point = median, .interval = hdi) %>%
  gather(key, value, -x, -.width, -.point, -.interval) %>%
  mutate(treatment = str_replace(key, "^mu_[:alpha:]+_(w[w|d])[\\.]*[(low|upp)er]*$", "\\1")) %>%
  group_by(x, .width, .point, .interval, treatment) %>%
  mutate(key = str_replace(key, "^(mu_[:alpha:]+)(_w[w|d])([\\.]*[(low|upp)er]*)$", "\\1\\3")) %>%
  spread(key, value) %>%
  rename(
    mu_gsmax = x, 
    lower_kcatc = mu_kcatc.lower,
    upper_kcatc = mu_kcatc.upper,
    lower_sco = mu_sco.lower,
    upper_sco = mu_sco.upper
  ) %>%
  mutate(
    lower_gsmax = 0,
    upper_gsmax = 0
  ) %>%
  mutate(Treatment = factor(toupper(treatment), levels = c("WW", "WD")))

r2 <- as.data.frame(m1b, c("gsmax_wd_spp", "gsmax_ww_spp", "kcatc_spp",
                           "sco_spp")) %>%
  mutate(iter = 1:nrow(.)) %>%
  gather(par, value, -iter) %>%
  mutate(
    spp = str_replace(par, "^[[:print:]]+_[a-z]*_*spp\\[([0-9]+)\\]$", "\\1"),
    trait = str_replace(par, "^([[:print:]]+)_spp\\[([0-9]+)\\]$", "\\1")
  ) %>%
  select(-par) %>%
  group_by(iter) %>%
  spread(trait, value) %>%
  summarize(
    r_kcatc_wd = cor(gsmax_wd, kcatc),
    r_sco_wd = cor(gsmax_wd, sco),
    r2_kcatc_wd = cor(gsmax_wd, kcatc) ^ 2,
    r2_sco_wd = cor(gsmax_wd, sco) ^ 2,
    r_kcatc_ww = cor(gsmax_ww, kcatc),
    r_sco_ww = cor(gsmax_ww, sco),
    r2_kcatc_ww = cor(gsmax_ww, kcatc) ^ 2,
    r2_sco_ww = cor(gsmax_ww, sco) ^ 2
  ) %>%
  ungroup() %>%
  select(-iter) %>%
  point_interval(.point = median, .interval = hdi)

fig2 %<>%
  ungroup() %>%
  mutate(lower = mu - se, upper = mu + se) %>%
  select(-se) %>%
  unite("trait_treatment", trait, treatment) %>%
  gather(term, value, -spp1, -trait_treatment) %>%
  separate(trait_treatment, c("trait", "treatment")) %>%
  group_by(treatment) %>%
  unite("term_trait", term, trait) %>%
  spread(term_trait, value) %>%
  mutate(Treatment = factor(toupper(treatment), levels = c("WW", "WD")))

l1 <- full_join(
  r2 %>%
    select_at(vars(matches("^r2_(kcatc|sco)_w(w|d)$"))) %>%
    gather(term, r2) %>%
    mutate(
      response = str_replace(term, "^r2_(kcatc|sco)_w(w|d)$", "\\1"),
      treatment = str_replace(term, "^r2_(kcatc|sco)_w(w|d)$", "w\\2"),
      r2 = round(r2, 2)
    ) %>%
    select(r2, response, treatment),
  m1_tab %>%
    filter(str_detect(term, "^b_(kcatc|sco)_gsmax_w[w|d]{1}$")) %>%
    mutate(
      response = str_replace(term, "^b_(kcatc|sco)_gsmax_w(w|d)$", "\\1"),
      treatment = str_replace(term, "^b_(kcatc|sco)_gsmax_w(w|d)$", "w\\2"),
      p = ifelse(p < 0.001, " < 0.001", str_c(" = ", round(p, 3)))
    ) %>%
    select(p, response, treatment),
  by = c("response", "treatment")
) %>% 
  mutate(
    label = glue("atop(paste(italic(R) ^ 2, \" = {r2}\"), paste(italic(P), \"{p}\"))", r2 = r2, p = p),
    mu_gsmax = ifelse(response == "kcatc", 1.25, 2.75),
    mu_kcatc = ifelse(response == "kcatc", 3.75, NA),
    mu_sco = ifelse(response == "kcatc", NA, 120),
    lower_gsmax = 0,
    upper_gsmax = 0,
    lower_kcatc = 0,
    upper_kcatc = 0,
    lower_sco = 0,
    upper_sco = 0
  ) %>%
  mutate(Treatment = factor(toupper(treatment), levels = c("WW", "WD")))

## Panel A ----

f2a <- ggplot(fig2, aes(
  x = mu_gsmax, y = mu_kcatc,
  xmin = lower_gsmax, xmax = upper_gsmax,
  ymin = lower_kcatc, ymax = upper_kcatc,
  fill = Treatment
)) +  
  facet_wrap(. ~ Treatment) +
  scale_x_continuous(limits = c(1.2, 2.8), breaks = seq(1.25, 2.75, 0.5)) +
  scale_y_continuous(limits = c(1.5, 3.75), breaks = seq(1.75, 3.75, 0.5)) +
  scale_fill_manual("Treatment", values = c("black", "white")) +
  geom_ribbon(data = b1, fill = "grey50") +
  geom_line(data = b1, size = 1.2, lineend = "round") +
  geom_errorbar() +
  geom_errorbarh() +
  geom_point(size = 3, shape = 21) +
  geom_text(data = filter(l1, response == "kcatc"), 
            aes(label = label), parse = TRUE, hjust = 0, vjust = 1) +
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

f2b <- ggplot(fig2, aes(
  x = mu_gsmax, y = mu_sco,
  xmin = lower_gsmax, xmax = upper_gsmax,
  ymin = lower_sco, ymax = upper_sco,
  fill = Treatment
)) +  
  facet_wrap(. ~ Treatment) +
  scale_x_continuous(limits = c(1.2, 2.8), breaks = seq(1.25, 2.75, 0.5)) +
  scale_y_continuous(limits = c(105, 120), breaks = seq(105, 120, 5)) +
  scale_fill_manual("Treatment", values = c("black", "white")) +
  geom_ribbon(data = b1, fill = "grey50") +
  geom_line(data = b1, size = 1.2, lineend = "round") +
  geom_errorbar() +
  geom_errorbarh() +
  geom_point(size = 3, shape = 21) +
  geom_text(data = filter(l1, response == "sco"), 
            aes(label = label), parse = TRUE, hjust = 1, vjust = 1) +
  ylab(expression(paste(italic(S)[C/O], " [mol ", mol ^ -1, "]"))) +
  xlab(expression(paste("Anatomical ", italic(g)[smax], " [mol ", H[2], "O ", m^-2~s^-1, "]"))) +
  theme_cdm() +
  theme(
    legend.position = "top",
    plot.margin = unit(c(0, 0, 0, 0), units = "cm")
  ) +
  NULL

plot_grid(f2a, f2b, ncol = 1, align = "h", labels = c("A", "B"), axis = "t")

ggsave("figures/fig2.pdf", w = 4.5, h = 6)

# Figure S3 ----
# Effect of drought treatment on gsmax by species

stomata <- read_csv("data/limonium-stomata.csv") %>%
  mutate(spp1 = as.numeric(as.factor(spp))) %>%
  group_by(spp) %>%
  summarize(spp1 = as.character(first(spp1)))

figS3 <- as.data.frame(m1b, pars = c("gsmax_ww_spp", "gsmax_wd_spp"))  %>%
  mutate(iter = 1:nrow(.)) %>%
  gather(key, value, -iter) %>%
  mutate(
    spp1 = str_replace(key, "^gsmax_w(w|d)_spp\\[([0-9]{1,2})\\]$", "\\2"),
    treatment = str_replace(key, "^gsmax_w(w|d)_spp\\[([0-9]{1,2})\\]$", 
                            "gsmax_w\\1")
  ) %>%
  select(-key) %>%
  group_by(iter, spp1) %>%
  spread(treatment, value) %>%
  ungroup() %>%
  mutate(gsmax_diff = gsmax_wd / gsmax_ww) %>%
  select(spp1, gsmax_diff) %>%
  group_by(spp1) %>%
  point_interval(.width = 0.95, .point = median, .interval = hdi) %>%
  arrange(gsmax_diff) %>%
  full_join(stomata, by = "spp1") %>%
  mutate(spp = factor(spp, levels = spp))

ggplot(figS3, aes(spp, gsmax_diff)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_pointinterval() +
  xlab("Species") +
  ylab(expression(frac(italic(g)[paste("smax",",WW")], italic(g)[paste("smax",",WD")]))) +
  theme_cdm() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    axis.title.y = element_text(angle = 0, vjust = 0.5)
  )

ggsave("figures/figS3.pdf", width = 4, height = 4)
