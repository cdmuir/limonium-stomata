source("r/header.R")

# Prepare data ----

stomata <- read_csv("data/limonium-stomata.csv") %>%
  group_by(spp, treatment, rep) %>%
  summarize(
    sd = sd[surface == "aB"] + sd[surface == "aD"],
    sp = (sp[surface == "aB"] + sp[surface == "aD"]) / 2
  )

lags <- read_csv("data/limonium-lags.csv")

stomata %<>% full_join(select(lags, spp, treatment, rep, la = LA1), 
                       by = c("spp", "treatment", "rep"))

stomata$treatment %<>% factor(levels = c("WW", "WD"))

# Fit model ----

if (run) {
  
  m3a <- brm(mvbind(sd, sp) ~ la * treatment + (1 | spp), 
             data = stomata, cores = 4, seed = 812890523) %>%
    add_criterion("loo", reloo = TRUE)
  m3b <- brm(mvbind(sd, sp) ~ la + la : treatment + (1 | spp), 
             data = stomata, cores = 4, seed = 393331214) %>%
    add_criterion("loo", reloo = TRUE)
  m3c <- brm(mvbind(sd, sp) ~ la + (1 | spp), 
             data = stomata, cores = 4, seed = 215291214) %>%
    add_criterion("loo", reloo = TRUE)

  write_rds(m3a, "objects/m3a.rds")
  write_rds(m3b, "objects/m3b.rds")
  write_rds(m3c, "objects/m3c.rds")

} else {
  
  m3a <- read_rds("objects/m3a.rds")
  m3b <- read_rds("objects/m3b.rds")
  m3c <- read_rds("objects/m3c.rds")

}

loo_compare(m3a, m3b, m3c)

# Parameter summary ----

pars3 <- c("b_sd_Intercept", "b_sp_Intercept", "b_sd_la", "b_sp_la")

m3_tab <- m3c$fit %>%
  tidy(pars = pars3, rhat = TRUE, ess = TRUE, conf.int = TRUE) %>%
  left_join(pars3 %>%
              as.data.frame(m3c$fit, pars = .) %>%
              summarize_all(get_p) %>%
              gather(term, p), by = "term")


write_rds(m3_tab, "objects/m3_tab.rds")

# Figure 4 ----

ex <- m3c %>%
  as.data.frame(pars = pars3) %>%
  mutate(iter = 1:nrow(.)) %>% 
  crossing(x = seq(min(stomata$la, na.rm = TRUE), 
                   max(stomata$la, na.rm = TRUE), length.out = 1e2)) %>%
  mutate(
    mu_sd = b_sd_Intercept + b_sd_la * x,
    mu_sp = b_sp_Intercept + b_sp_la * x
  ) %>%
  group_by(x) %>%
  summarize(
    est_sd = median(mu_sd),
    lci_sd = hdi(mu_sd)[1, 1],
    uci_sd = hdi(mu_sd)[1, 2],
    est_sp = median(mu_sp),
    lci_sp = hdi(mu_sp)[1, 1],
    uci_sp = hdi(mu_sp)[1, 2]
  )

fig4 <- stomata %>%
  filter(!is.na(la), !is.na(sd), !is.na(sp)) %>%
  select(spp, treatment, la, sd, sp) %>%
  group_by(spp, treatment) %>%
  summarize(
    sd_la = sd(la, na.rm = TRUE),
    sd_sd = sd(sd, na.rm = TRUE),
    sd_sp = sd(sp, na.rm = TRUE),
    la = mean(la, na.rm = TRUE),
    sd = mean(sd, na.rm = TRUE),
    sp = mean(sp, na.rm = TRUE),
    n = n()
  ) %>%
  mutate(
    se_la = sd_la / sqrt(n),
    se_sd = sd_sd / sqrt(n),
    se_sp = sd_sp / sqrt(n)
  ) %>%
  ungroup() %>%
  mutate_if(function(.x) inherits(.x, "units"), drop_units)

r2 <- bayes_R2(m3c)

## Panel A ----

f4a <- ggplot(fig4, aes(
  x = la, y = sd,
  xmin = la - se_la, xmax = la + se_la,
  ymin = sd - se_sd, ymax = sd + se_sd,
  fill = treatment
)) +  
  scale_fill_manual("Treatment", values = c("black", "white")) +
  geom_ribbon(data = ex, 
              mapping = aes(x = x, ymin = lci_sd, ymax = uci_sd), 
              inherit.aes = FALSE, fill = "grey50") +
  geom_line(data = ex, mapping = aes(x = x, y = est_sd), 
            size = 1.2, inherit.aes = FALSE) +
  geom_errorbar() +
  geom_errorbarh() +
  geom_point(size = 3, shape =21) +
  annotate(
    "text", max(fig4$la, na.rm = TRUE), max(fig4$sd, na.rm = TRUE), 
    hjust = 0.25, vjust = 0, parse = TRUE,
    label = glue("atop(paste(italic(R) ^ 2, \" = {r2}\"), paste(italic(P), \" = {p}\"))",
                 r2 = round(r2["R2sd", "Estimate"], 2), 
                 p = signif(m3_tab$p[m3_tab$term == "b_sd_la"], 1))
  ) +
  xlab(expression(paste("Leaf area [m", m ^ 2, "]"))) +
  ylab(expression(paste("Stomatal density [m", m^-2, "]"))) +
  theme_cdm() +
  theme(
    legend.position = "top",
    plot.margin = unit(c(0, 0, 0, 0), units = "cm")
  ) +
  NULL

## Panel B ----

f4b <- ggplot(fig4, aes(
  x = la, y = sp,
  xmin = la - se_la, xmax = la + se_la,
  ymin = sp - se_sp, ymax = sp + se_sp,
  fill = treatment
)) +  
  scale_y_continuous(limits = c(20, 35)) +
  scale_fill_manual("Treatment", values = c("black", "white")) +
  geom_ribbon(data = ex, 
              mapping = aes(x = x, ymin = lci_sp, ymax = uci_sp), 
              inherit.aes = FALSE, fill = "grey50") +
  geom_line(data = ex, mapping = aes(x = x, y = est_sp), 
            size = 1.2, inherit.aes = FALSE) +
  geom_errorbar() +
  geom_errorbarh() +
  geom_point(size = 3, shape =21  ) +
  annotate(
    "text", max(fig4$la, na.rm = TRUE), 32, 
    hjust = 0.25, vjust = 0.15, parse = TRUE,
    label = glue("atop(paste(italic(R) ^ 2, \" = {r2}\"), paste(italic(P), \" = {p}\"))",
                 r2 = round(r2["R2sp", "Estimate"], 2), 
                 p = signif(m3_tab$p[m3_tab$term == "b_sp_la"], 3))
  ) +
  xlab(expression(paste("Leaf area [m", m ^ 2, "]"))) +
  ylab(expression(paste("Stomatal pore length [", mu, "m]"))) +  
  theme_cdm() +
  theme(
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), units = "cm")
  ) +
  NULL

plot_grid(f4a, f4b, ncol = 1, align = "hv", labels = c("A", "B"), axis = "t")

ggsave("figures/fig4.pdf", w = 3.25, h = 6.5)
