source("r/header.R")

stom_dens <- read_excel("data/Limo_mac_V3_Tofol.xlsx", 
                        sheet = "stomatal density") %>%
  filter(!is.na(sd)) %>%
  group_by(spp, treatment, rep, surface) %>%
  summarize(sd = mean(sd), si = mean(si))

stom_pore <- read_excel("data/Limo_mac_V3_Tofol.xlsx", 
                        sheet = "stomatal pore length") %>%
  filter(!is.na(sp)) %>%
  group_by(spp, treatment, rep, surface) %>%
  summarize(sp = mean(sp))

stom_all <- full_join(stom_dens, stom_pore, by = c("spp", "treatment", "rep", "surface"))

write_csv(stom_all, "data/limonium-stomata.csv")

lags <- read_excel("data/Limo_mac_V3_Tofol.xlsx", 
                   sheet = "leaf area and gs") %>%
  select(spp = Sp...3, treatment, rep = plant, LA1, gs) %>%
  filter(!(is.na(LA1) & is.na(gs)))

write_csv(lags, "data/limonium-lags.csv")
