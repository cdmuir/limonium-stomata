source("r/header.R")

stomata <- read_csv("data/limonium-stomata.csv")
stomata$rep %<>% as.character()
stomata$treatment %<>% factor(levels = c("WW", "WD"))

# 1. Correlation between ab- and adaxial stomatal density ----

s1 <- stomata %>%
  select(spp, treatment, rep, surface, sd) %>%
  filter(!is.na(sd)) %>%
  spread(surface, sd) %>%
  mutate_if(is.numeric, log) %>%
  rename_if(is.numeric, str_replace, pattern = "([[:alpha:]]+)", 
            replacement = "log_\\1")

if (run) {
  
  mS1a <- brm(mvbind(log_aB, log_aD) ~ treatment + (treatment | p | spp), 
             data = s1, chains = 4, cores = 4, seed = 799241667)
  write_rds(mS1a, "objects/mS1a.rds")
  
} else {
  
  mS1a <- read_rds("objects/mS1a.rds")

}

# Parameter summary ----

parsS1a <- c("b_logaB_Intercept", "b_logaD_Intercept", 
             "b_logaB_treatmentWD", "b_logaD_treatmentWD", 
             "cor_spp__logaB_Intercept__logaD_Intercept", 
             "rescor__logaB__logaD")

mS1a_tab <- mS1a$fit %>%
  tidy(pars = parsS1a, rhat = TRUE, ess = TRUE, conf.int = TRUE) %>%
  left_join(parsS1a %>%
              as.data.frame(mS1a$fit, pars = .) %>%
              summarize_all(get_p) %>%
              gather(term, p), by = "term")

write_rds(mS1a_tab, "objects/mS1a_tab.rds")

# Figure S1a ----

new1a <- crossing(spp = unique(s1$spp), treatment = c("WW", "WD"))

figS1a <- predict(mS1a, new1a) %>%
  set_colnames(c("estimate", "se", "lowerci", "upperci")) %>%
  as.data.frame() %>%
  bind_cols(new1a)

figS1a$treatment %<>% factor(levels = c("WW", "WD"))

fS1a <- ggplot(figS1a, aes(x = estimate.logaB, y = estimate.logaD,
                          xmin = estimate.logaB - se.logaB,
                          xmax = estimate.logaB + se.logaB,
                          ymin = estimate.logaD - se.logaD,
                          ymax = estimate.logaD + se.logaD,
                          fill = treatment, colour = treatment)) +
  scale_x_continuous(trans = "exp",
                     limits = c(log(30), log(85)),
                     breaks = log(seq(30, 80, 10)),
                     labels = seq(30, 80, 10)) +
  scale_y_continuous(trans = "exp",
                     limits = c(log(30), log(85)),
                     breaks = log(seq(30, 80, 10)),
                     labels = seq(30, 80, 10)) +
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed") +
  geom_errorbar(colour = "black") +
  geom_errorbarh(colour = "black") +
  geom_point(size = 3, shape = 21) +
  scale_fill_manual("Treatment", values = c("black", "white")) +
  scale_colour_manual("Treatment", values = c("white", "black")) +
  xlab(expression(paste("Abaxial stomatal density [m", m^-2, "]"))) +
  ylab(expression(paste("Adaxial stomatal density [m", m^-2, "]"))) +
  theme_cdm() +
  theme(
    legend.position = "top",
    plot.margin = unit(c(0, 0, 0, 0), units = "cm")
  ) +
  NULL

# 2. Correlation between ab- and adaxial stomatal pore length ----

s2 <- stomata %>%
  select(spp, treatment, rep, surface, sp) %>%
  filter(!is.na(sp)) %>%
  spread(surface, sp) %>%
  mutate_if(is.numeric, log) %>%
  rename_if(is.numeric, str_replace, pattern = "([[:alpha:]]+)", 
            replacement = "log_\\1")

if (run) {
  
  mS1b <- brm(mvbind(log_aB, log_aD) ~ treatment + (treatment | p | spp), 
             data = s2, chains = 4, cores = 4, seed = 100623159)
  write_rds(mS1b, "objects/mS1b.rds")

} else {
  
  mS1b <- read_rds("objects/mS1b.rds")

}

# Parameter summary ----

parsS1b <- c("b_logaB_Intercept", "b_logaD_Intercept", 
             "b_logaB_treatmentWD", "b_logaD_treatmentWD", 
             "cor_spp__logaB_Intercept__logaD_Intercept", 
             "rescor__logaB__logaD")

mS1b_tab <- mS1b$fit %>%
  tidy(pars = parsS1b, rhat = TRUE, ess = TRUE, conf.int = TRUE) %>%
  left_join(parsS1b %>%
              as.data.frame(mS1b$fit, pars = .) %>%
              summarize_all(get_p) %>%
              gather(term, p), by = "term")

write_rds(mS1b_tab, "objects/mS1b_tab.rds")

# Figure S1b ----

new1b <- crossing(spp = unique(s2$spp), treatment = c("WW", "WD"))

figS1b <- predict(mS1b, new1b) %>%
  set_colnames(c("estimate", "se", "lowerci", "upperci")) %>%
  as.data.frame() %>%
  bind_cols(new1b)

figS1b$treatment %<>% factor(levels = c("WW", "WD"))

fS1b <- ggplot(figS1b, aes(x = estimate.logaB, y = estimate.logaD,
                           xmin = estimate.logaB - se.logaB,
                           xmax = estimate.logaB + se.logaB,
                           ymin = estimate.logaD - se.logaD,
                           ymax = estimate.logaD + se.logaD,
                           fill = treatment, colour = treatment)) +
  scale_x_continuous(trans = "exp",
                     limits = c(log(19), log(35)),
                     breaks = log(seq(21, 33, 3)),
                     labels = seq(21, 33, 3)) +
  scale_y_continuous(trans = "exp",
                     limits = c(log(19), log(35)),
                     breaks = log(seq(21, 33, 3)),
                     labels = seq(21, 33, 3)) +
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed") +
  geom_errorbar(colour = "black") +
  geom_errorbarh(colour = "black") +
  geom_point(size = 3, shape = 21) +
  scale_fill_manual("Treatment", values = c("black", "white")) +
  scale_colour_manual("Treatment", values = c("white", "black")) +
  xlab(expression(paste("Abaxial stomatal pore length [", mu, "m]"))) +
  ylab(expression(paste("Adaxial stomatal pore length [", mu, "m]"))) +
  theme_cdm() +
  theme(
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), units = "cm")
  ) +
  NULL

plot_grid(fS1a, fS1b, nrow = 1, align = "hv", labels = c("A", "B"), axis = "t")

ggsave("figures/figS1.pdf", w = 6.5, h = 3.25)
