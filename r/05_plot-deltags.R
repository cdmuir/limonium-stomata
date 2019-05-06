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

# Fit models ----

if (run) {
  
  ma <- brm(log(gs) ~ treatment + (treatment | spp), 
             data = mutate_if(stomata, ~ inherits(.x, "units"), drop_units),
             cores = 4, seed = 461457053) %>%
    add_criterion("loo", reloo = TRUE)
  mb <- brm(log(gs) ~ treatment + (1 | spp), 
             data = mutate_if(stomata, ~ inherits(.x, "units"), drop_units),
             cores = 4, seed = 407140810) %>%
    add_criterion("loo", reloo = TRUE)
  
  write_rds(ma, "objects/ma_gsdiff.rds")
  write_rds(mb, "objects/mb_gsdiff.rds")

} else {
  
  ma <- read_rds("objects/ma_gsdiff.rds")
  mb <- read_rds("objects/mb_gsdiff.rds")

}

loo_compare(ma, mb, criterion = "loo")

gs_diff <- as.data.frame(ma$fit) %>%
  select_at(vars(matches("^(r|b)_"))) %>%
  mutate(iter = 1:nrow(.)) %>%
  mutate_at(vars(matches("^r_spp\\[([A-Z]{3}),Intercept\\]$")),
            ~ .x + b_Intercept) %>%
  mutate_at(vars(matches("^r_spp\\[([A-Z]{3}),treatmentWW\\]$")),
            ~ .x + b_treatmentWW) %>%
  select(-b_Intercept, -b_treatmentWW) %>%
  gather(key, value, -iter) %>%
  mutate(
    spp = str_replace(key, "^r_spp\\[([A-Z]{3}),([:alpha:]+)\\]$", "\\1"),
    treatment = str_replace(key, "^r_spp\\[([A-Z]{3}),([:alpha:]+)\\]$", "\\2")
  ) %>%
  select(-key) %>%
  group_by(iter, spp) %>%
  spread(treatment, value) %>%
  mutate(gs_diff = exp(treatmentWW)) %>%
  ungroup() %>%
  select(iter, spp, gs_diff) 

# Figure S4 ----

figS4 <- gs_diff %>%
  select(-iter) %>%
  group_by(spp) %>%
  point_interval(.width = 0.95, .point = median, .interval = hdci) %>%
  arrange(gs_diff) %>%
  mutate(spp = factor(spp, levels = spp))

ggplot(figS4, aes(spp, gs_diff)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_pointinterval() +
  xlab("Species") +
  ylab(expression(frac(italic(g)[paste("s",",WW")], italic(g)[paste("s",",WD")]))) +
  theme_cdm() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    axis.title.y = element_text(angle = 0, vjust = 0.5)
  )

ggsave("figures/figS4.pdf", width = 4, height = 4)

# Correlation between anatomical and physiological response

m1b <- read_rds("objects/m1b.rds")

gsmax_diff <- as.data.frame(m1b, pars = c("gsmax_ww_spp", "gsmax_wd_spp"))  %>%
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
  select(iter, spp1, gsmax_diff)

gs_diff %<>% mutate(spp1 = as.character(as.numeric(as.factor(spp))))

full_join(gsmax_diff, gs_diff, by = c("iter", "spp1")) %>%
  group_by(iter) %>%
  summarize(r = cor(gsmax_diff, gs_diff)) %>%
  select(-iter) %>%
  # pull(r) %>%
  # get_p()
  point_interval(.width = 0.95, .point = median, .interval = hdi)
