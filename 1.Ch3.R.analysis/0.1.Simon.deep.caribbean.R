library(tidyverse)
#library(brms)
library(modelr)
#library(tidybayes)
library(fishflux)
library(fishualize)
#library(rfishbase)
library(vegan)
library(ggrepel)
library(clustsig)
library(dendextend)
library(ggdendro)
library(tinter)
library(mgcv)
library(purrr)
#library(patchwork)
library(FD)
#library(AICcmodavg)

#library(devtools) #to enable the download of an archived version of clustsig
#install_github("cran/clustsig")
library(clustsig)

carib.pal <- fish(4, option = "Bodianus_rufus")

# read in raw data files
bonaire <- read.csv(file = "raw fish counts/Bonaire_Raw_2020.csv") %>%
  mutate(hab = "deep")
curacao <- read.csv(file = "raw fish counts/Curacao_Raw_2020.csv") %>%
  mutate(hab = "deep")
roatan_deep <- read.csv(file = "raw fish counts/Roatan_Raw_2020_onlyDeepReefs.csv") %>%
  mutate(hab = "deep")
roatan <- read.csv(file = "raw fish counts/Roatan_Raw_2020.csv") %>%
  mutate(hab = "shallow")
statia <- read.csv(file = "raw fish counts/Statia_raw_2020.csv") %>%
  mutate(hab = "deep")

# read in correction factors
bon.corr <- read.csv("CorrectionFactors_SamplingEfforts/BonCoeff_2020.csv") %>%
  separate(depth, into = c("dband", "m")) %>%
  mutate(dband = as.numeric(dband)) %>%
  mutate(location = "Bonaire")
cur.corr <- read.csv("CorrectionFactors_SamplingEfforts/CurCoeff_2020.csv") %>%
  separate(depth, into = c("dband", "m")) %>%
  mutate(dband = as.numeric(dband)) %>%
  mutate(location = "Curacao")
roa.corr <- read.csv("CorrectionFactors_SamplingEfforts/RoaCoeff_2020.csv") %>%
  separate(depth, into = c("dband", "m")) %>%
  mutate(dband = as.numeric(dband)) %>%
  mutate(location = "Roatan")
sta.corr <- read.csv("CorrectionFactors_SamplingEfforts/StaCoeff_2020.csv") %>%
  separate(X, into = c("dband", "m")) %>%
  rename(COEF = "coeff") %>%
  mutate(dband = as.numeric(dband)) %>%
  mutate(location = "St. Eustatius")

coefs <- bind_rows(bon.corr, cur.corr, roa.corr, sta.corr)

carib.deep <- bind_rows(bonaire, curacao, roatan, roatan_deep, statia) %>%
  mutate(dband = case_when(meters >= 0 & meters < 10 ~ 0,
                           meters >= 10 & meters < 20 ~ 10,
                           meters >= 20 & meters < 30 ~ 20,
                           meters >= 30 & meters < 40 ~ 30,
                           meters >= 40 & meters < 50 ~ 40,
                           meters >= 50 & meters < 60 ~ 50,
                           meters >= 60 & meters < 70 ~ 60,
                           meters >= 70 & meters < 80 ~ 70,
                           meters >= 80 & meters < 90 ~ 80,
                           meters >= 90 & meters < 100 ~ 90,
                           meters >= 100 & meters < 110 ~ 100,
                           meters >= 110 & meters < 120 ~ 110,
                           meters >= 120 & meters < 130 ~ 120,
                           meters >= 130 & meters < 140 ~ 130,
                           meters >= 140 & meters < 150 ~ 140,
                           meters >= 150 & meters < 160 ~ 150,
                           meters >= 160 & meters < 170 ~ 160,
                           meters >= 170 & meters < 180 ~ 170,
                           meters >= 180 & meters < 190 ~ 180,
                           meters >= 190 & meters < 200 ~ 190,
                           meters >= 200 & meters < 210 ~ 200,
                           meters >= 210 & meters < 220 ~ 210,
                           meters >= 220 & meters < 230 ~ 220,
                           meters >= 230 & meters < 240 ~ 230,
                           meters >= 240 & meters < 250 ~ 240,
                           meters >= 250 & meters < 260 ~ 250,
                           meters >= 260 & meters < 270 ~ 260,
                           meters >= 270 & meters < 280 ~ 270,
                           meters >= 280 & meters < 290 ~ 280,
                           meters >= 290 & meters < 300 ~ 290,
                           meters >= 300 & meters < 310 ~ 300,
                           meters >= 310 & meters < 320 ~ 310,
                           meters >= 320 & meters < 330 ~ 320,
                           meters >= 330 & meters < 340 ~ 330,
                           meters >= 340 & meters < 350 ~ 340,
                           meters >= 350 & meters < 360 ~ 350,
                           meters >= 360 & meters < 370 ~ 360,
                           meters >= 370 & meters < 380 ~ 370,
                           meters >= 380 & meters < 390 ~ 380,
                           meters >= 390 & meters < 400 ~ 390,
                           meters >= 400 & meters < 410 ~ 400,
                           meters >= 410 & meters < 420 ~ 410,
                           meters >= 420 & meters < 430 ~ 420,
                           meters >= 430 & meters < 440 ~ 430,
                           meters >= 440 & meters < 450 ~ 440,
                           meters >= 450 & meters < 460 ~ 450,
                           meters >= 460 & meters < 470 ~ 460,
                           meters >= 470 & meters < 480 ~ 470,
                           meters >= 480 & meters < 490 ~ 480,
                           meters >= 490 & meters < 500 ~ 490)) %>%
  rename(species = "name") %>%
  group_by(species, dband, location, hab) %>%
  summarize(abundance = n()) %>%
  left_join(coefs) %>%
  mutate(abu.corr = abundance*COEF)
  
summary(carib.deep)

# look for false names -- need to replace a bunch with FB comopatible names
# choosing close congeners
# mute for quick analysis

tax.names <- data.frame("species" = c("Bathyanthias sp.","Bembrops sp.","Cantherhines sp.","Carangoides ruber","Caulolatilus plumieri",
                                      "Coryphopterus curasub","Decapterus sp.","Decodon shallow","Epigonidae","Epinephalus striatus",
                                      "Fistularia sp.","Flammeo marianus","Haemulon vittata","Lipogramma idabeli","Lipogramma regium",
                                      "Myctoperca bonaci","Neoniphon coruscum","Palatogobius incendius","Phenacoscorpius nebrius",
                                      "Pinnichthys aimoriensis","Plectranthias sp.","Polylepion sp.","Psilotris laurae","Rypticus saponaceous",
                                      "Scarus viride","Serranus tabacarrius","Stegastes viride","Synagrops egretta","Varicus adamsi",
                                      "Varicus cephalocellatus","Varicus decorum","Varicus veliguttatus")) %>%
  mutate(species.fb = case_when(species == "Bathyanthias sp." ~ "Bathyanthias mexicanus",
                                species == "Bembrops sp." ~ "Bembrops gobioides",
                                species == "Cantherhines sp." ~ "Cantherhines pullus",
                                species == "Carangoides ruber" ~ "Caranx ruber",
                                species == "Caulolatilus plumieri" ~ "Malacanthus plumieri",
                                species == "Coryphopterus curasub" ~ "Coryphopterus glaucofraenum",
                                species == "Decapterus sp." ~ "Decapterus macarellus",
                                species == "Decodon shallow" ~ "Decodon puellaris",
                                species == "Epigonidae" ~ "Epigonus denticulatus",
                                species == "Epinephalus striatus" ~ "Epinephelus striatus",
                                species == "Fistularia sp." ~ "Fistularia commersonii",
                                species == "Flammeo marianus" ~ "Neoniphon marianus",
                                species == "Haemulon vittata" ~ "Haemulon vittatum",
                                species == "Lipogramma idabeli" ~ "Lipogramma evides",
                                species == "Lipogramma regium" ~ "Lipogramma evides",
                                species == "Myctoperca bonaci" ~ "Mycteroperca bonaci",
                                species == "Neoniphon coruscum" ~ "Sargocentron coruscum",
                                species == "Palatogobius incendius" ~ "Palatogobius grandoculus",
                                species == "Phenacoscorpius nebrius" ~ "Phenacoscorpius nebris",
                                species == "Pinnichthys aimoriensis" ~ "Chriolepis zebra",
                                species == "Plectranthias sp." ~ "Plectranthias garrupellus",
                                species == "Polylepion sp." ~ "Polylepion russelli",
                                species == "Psilotris laurae" ~ "Psilotris boehlkei",
                                species == "Rypticus saponaceous" ~ "Rypticus saponaceus",
                                species == "Scarus viride" ~ "Sparisoma viride",
                                species == "Serranus tabacarrius" ~ "Serranus tabacarius",
                                species == "Stegastes viride" ~ "Sparisoma viride",
                                species == "Synagrops egretta" ~ "Synagrops bellus",
                                species == "Varicus adamsi" ~ "Varicus bucca",
                                species == "Varicus cephalocellatus" ~ "Varicus bucca",
                                species == "Varicus decorum" ~ "Varicus bucca",
                                species == "Varicus veliguttatus" ~ "Varicus bucca"))
carib.deep.fishbase <- carib.deep %>%
  left_join(tax.names) %>%
  distinct() %>%
  mutate(species.fb = case_when(is.na(species.fb) ~ species,
                                TRUE ~ species.fb))
# 

# sp.list <- data.frame("species.fb" = carib.deep.fishbase$species.fb) %>%
#   distinct()
# 
# lwl <- purrr::map(as.list(sp.list$species.fb), possibly(fishflux::find_lw, otherwise =  NA_real_))
# lwl2 <- lapply(lwl, function(x){
#   if (is.na(x)){
#     res <- data.frame(
#       species = NA, lwa_m = NA, lwa_sd = NA, lwb_m = NA, lwb_sd = NA
#     )
#     return(res)
#   } else{
#     return(x)
#   }
# })
# 
# lw_data <- plyr::ldply(lwl2) %>%
#   rename(species.fb = "species")

# get size data from master list
# actino.all <- read.csv("actino_all.csv") %>%
#   drop_na(Length) %>%
#   rename(species.fb = "sciname")
# 
# 
# actino.gen <- actino.all %>%
#   group_by(Genus) %>%
#   summarize(maxlgen = mean(Length), maxdepthgen = max(DepthRangeDeep)) %>%
#   rename(genus = "Genus")
# 
# 
# lmax.fb <- data.frame("species.fb" = sp.list$species.fb) %>%
#   separate(species.fb, into = c("genus", "species"), remove = F) %>%
#   left_join(actino.all) %>%
#   select(species.fb, genus, species, DepthRangeDeep, Length) %>%
#   left_join(actino.gen) %>%
#   mutate(maxTL = case_when(is.na(Length) ~ maxlgen,
#                            TRUE ~ Length),
#          maxDepth = case_when(is.na(DepthRangeDeep) ~ maxdepthgen,
#                               TRUE ~ DepthRangeDeep)) %>%
#   select(species.fb, maxTL, maxDepth)
# 
# carib.deep.traits <- carib.deep.fishbase %>%
#   left_join(lw_data) %>%
#   left_join(lmax.fb) %>%
#   mutate(maxTL = case_when(species == "Vomerogobius flavus" ~ 3,
#                            species == "Sphyraenops bairdianus" ~ 17,
#                            species == "Choranthias tenuis" ~ 12.5,
#                            TRUE ~ maxTL)) %>%
#   mutate(length = maxTL*0.75) %>%
#   mutate(biomass = abu.corr*(lwa_m*length^lwb_m))

# 
# write.csv(carib.deep.traits, file = "carib.deep.traits.csv")

carib.deep.traits <- read.csv(file = "carib.deep.traits.csv") %>%
  select(-X)

carib.deep.traits.species <- carib.deep.traits %>%
  select(species, maxTL, maxDepth) %>%
  distinct()
# 
# write.csv(carib.deep.traits.species, file = "carib.deep.traits.species.csv")

carib.fish.sum <- carib.deep.traits %>%
  group_by(dband, location, hab, species) %>%
  summarize(tot.abu = sum(abu.corr), u.abu = sum(abundance), tot.biom = sum(biomass)) %>%
  ungroup() %>%
  group_by(dband, location, hab) %>%
  summarize(abundance = sum(tot.abu), u.abundance = sum(u.abu), biomass = sum(tot.biom), spric = n()) 

#Abundances across locations

cur.abu <- carib.fish.sum %>%
  filter(location == "Curacao")

cur.abu.gam <- gam(log10(abundance) ~ s(dband), data = cur.abu)
summary(cur.abu.gam)
par(mfrow = c(2,2))
gam.check(cur.abu.gam)
cur.abu.new <- data.frame(dband = as.numeric(seq(min(cur.abu$dband), 
                                                 max(cur.abu$dband), 
                                                 length.out = 500)))
cur.abu.pred <- as.data.frame(predict(cur.abu.gam, data.frame(cur.abu.new), se.fit = T)) %>%
  rename(abu_predicted = 1) %>%
  mutate(abu_predicted_raw = exp(abu_predicted)) %>%
  mutate(lci = abu_predicted-1.96*se.fit,
         uci = abu_predicted+1.96*se.fit) %>%
  bind_cols(cur.abu.new)

cur.abu.plot <- ggplot(cur.abu.pred, aes(x = dband, y = abu_predicted)) +
  geom_line(lty = 1, color = carib.pal[1]) +
  geom_ribbon(aes(x = dband, ymin = lci, ymax = uci), alpha = 0.2, fill = carib.pal[1]) +
  geom_point(data = cur.abu, aes(x = dband, y = log10(abundance)), shape = 21, alpha = 0.75, color = "grey23", fill = carib.pal[1]) +
  geom_vline(xintercept = 100, lty = 2, color = "grey69") +
  geom_vline(xintercept = 150, lty = 2, color = "grey69") +
  geom_vline(xintercept = 230, lty = 2, color = "grey69") +
  geom_vline(xintercept = 290, lty = 2, color = "grey69") +
  theme_bw() +
  ylab("log10 abundance per depth band") +
  xlab("Depth (m)") +
  scale_x_continuous(limits = c(40,300), breaks = seq(40,300,10))
cur.abu.plot


### Bonaire

bon.abu <- carib.fish.sum %>%
  filter(location == "Bonaire")

bon.abu.gam <- gam(log10(abundance) ~ s(dband), data = bon.abu)
summary(bon.abu.gam)
par(mfrow = c(2,2))
gam.check(bon.abu.gam)
bon.abu.new <- data.frame(dband = as.numeric(seq(min(bon.abu$dband), 
                                                 max(bon.abu$dband), 
                                                 length.out = 500)))
bon.abu.pred <- as.data.frame(predict(bon.abu.gam, data.frame(bon.abu.new), se.fit = T)) %>%
  rename(abu_predicted = 1) %>%
  mutate(abu_predicted_raw = exp(abu_predicted)) %>%
  mutate(lci = abu_predicted-1.96*se.fit,
         uci = abu_predicted+1.96*se.fit) %>%
  bind_cols(bon.abu.new)

bon.abu.plot <- ggplot(bon.abu.pred, aes(x = dband, y = abu_predicted)) +
  geom_line(lty = 1, color = carib.pal[2]) +
  geom_ribbon(aes(x = dband, ymin = lci, ymax = uci), alpha = 0.2, fill = carib.pal[2]) +
  geom_point(data = bon.abu, aes(x = dband, y = log10(abundance)), shape = 22, alpha = 0.75, color = "grey23", fill = carib.pal[2]) +
  geom_vline(xintercept = 100, lty = 2, color = "grey69") +
  geom_vline(xintercept = 160, lty = 2, color = "grey69") +
  geom_vline(xintercept = 230, lty = 2, color = "grey69") +
  geom_vline(xintercept = 280, lty = 2, color = "grey69") +
  theme_bw() +
  ylab("log10 abundance per depth band") +
  xlab("Depth (m)") +
  scale_x_continuous(limits = c(40,300), breaks = seq(40,300,10))
bon.abu.plot


### Roatan

roa.abu <- carib.fish.sum %>%
  filter(location == "Roatan") %>%
  # filter(dband < 300) %>%
  filter(hab == "deep")

roa.abu.gam <- gam(log10(abundance) ~ s(dband), data = roa.abu, method = "REML")
summary(roa.abu.gam)
par(mfrow = c(2,2))
gam.check(roa.abu.gam)
roa.abu.new <- data.frame(dband = as.numeric(seq(min(roa.abu$dband), 
                                                 max(roa.abu$dband), 
                                                 length.out = 500)))
roa.abu.pred <- as.data.frame(predict(roa.abu.gam, data.frame(roa.abu.new), se.fit = T)) %>%
  rename(abu_predicted = 1) %>%
  mutate(abu_predicted_raw = exp(abu_predicted)) %>%
  mutate(lci = abu_predicted-1.96*se.fit,
         uci = abu_predicted+1.96*se.fit) %>%
  bind_cols(roa.abu.new)

roa.abu.plot <- ggplot(roa.abu.pred, aes(x = dband, y = abu_predicted)) +
  geom_line(lty = 1, color = carib.pal[3]) +
  geom_ribbon(aes(x = dband, ymin = lci, ymax = uci), alpha = 0.2, fill = carib.pal[3]) +
  geom_point(data = roa.abu, aes(x = dband, y = log10(abundance)), shape = 23, alpha = 0.75, color = "grey23", fill = carib.pal[3]) +
  geom_vline(xintercept = 90, lty = 2, color = "grey69") +
  geom_vline(xintercept = 140, lty = 2, color = "grey69") +
  geom_vline(xintercept = 200, lty = 2, color = "grey69") +
  geom_vline(xintercept = 290, lty = 2, color = "grey69") +
  theme_bw() +
  ylab("log10 abundance per depth band") +
  xlab("Depth (m)") 
  # scale_x_continuous(limits = c(40,max(roa.abu$dband)), breaks = seq(40,max(roa.abu$dband),10))
roa.abu.plot


### Statia

sta.abu <- carib.fish.sum %>%
  filter(location == "St. Eustatius") %>%
  filter(dband < 300)

sta.abu.gam <- gam(log10(abundance) ~ s(dband), data = sta.abu, method = "REML")
summary(sta.abu.gam)
par()
gam.check(sta.abu.gam)
sta.abu.new <- data.frame(dband = as.numeric(seq(min(sta.abu$dband), 
                                                 max(sta.abu$dband), 
                                                 length.out = 500)))
sta.abu.pred <- as.data.frame(predict(sta.abu.gam, data.frame(sta.abu.new), se.fit = T)) %>%
  rename(abu_predicted = 1) %>%
  mutate(abu_predicted_raw = exp(abu_predicted)) %>%
  mutate(lci = abu_predicted-1.96*se.fit,
         uci = abu_predicted+1.96*se.fit) %>%
  bind_cols(sta.abu.new)

sta.abu.plot <- ggplot(sta.abu.pred, aes(x = dband, y = abu_predicted)) +
  geom_line(lty = 1, color = carib.pal[4]) +
  geom_ribbon(aes(x = dband, ymin = lci, ymax = uci), alpha = 0.2, fill = carib.pal[4]) +
  geom_point(data = sta.abu, aes(x = dband, y = log10(abundance)), shape = 24, alpha = 0.75, color = "grey23", fill = carib.pal[4]) +
  geom_vline(xintercept = 90, lty = 2, color = "grey69") +
  geom_vline(xintercept = 140, lty = 2, color = "grey69") +
  geom_vline(xintercept = 200, lty = 2, color = "grey69") +
  geom_vline(xintercept = 270, lty = 2, color = "grey69") +
  theme_bw() +
  ylab("log10 abundance per depth band") +
  xlab("Depth (m)") +
  scale_x_continuous(limits = c(40,max(sta.abu$dband)), breaks = seq(40,max(sta.abu$dband),10))
sta.abu.plot


abundance.plots <- (cur.abu.plot + bon.abu.plot + roa.abu.plot + sta.abu.plot) + plot_annotation(tag_levels = "A")
ggsave(abundance.plots, file = "abundance.plots.lin.png", width = 14, height = 8)


biomass.plot <- ggplot(carib.fish.sum, aes(x = dband, y = log10(biomass), color = location, fill = location)) +
  geom_point(aes(shape = location), alpha = 0.75) +
  geom_smooth(method = "gam") +
  theme_bw() +
  scale_fill_fish_d(option = "Bodianus_rufus") +
  scale_color_fish_d(option = "Bodianus_rufus") +
  facet_wrap(.~location, scales = "free")

ggsave(biomass.plot, file = "biomass.png", width = 8, height = 6)


spric.plot <- ggplot(carib.fish.sum, aes(x = dband, y = log10(spric), color = location, fill = location)) +
  geom_point(aes(shape = location), alpha = 0.75) +
  geom_smooth(method = "gam") +
  theme_bw() +
  scale_fill_fish_d(option = "Bodianus_rufus") +
  scale_color_fish_d(option = "Bodianus_rufus") +
  facet_wrap(.~location)

spric.plot



### size structure
carib.size <- carib.deep.traits %>%
  select(species, species.fb, dband, location, hab, abundance, COEF, abu.corr, maxTL, length, biomass, lwa_m, lwb_m) %>%
  as_tibble() %>%
  mutate(obs = map(abundance, ~rep_len(1, .x))) %>%
  unnest() %>%
  select(-abundance, -obs)

carib.size %>%
  filter(species == "Abudefduf saxatilis")




#### community composition ####
carib.deep.fishbase <- read.csv("file.fish.carib.grouped.csv")
carib.deep.com <- carib.deep.fishbase %>%
  select(species, dband, location, hab, abundance, abu.corr) %>%
  filter(hab == "deep") %>%
  filter(dband <= 300) %>%
  # mutate(dband.coarse = case_when(dband >= 0 & dband < 50 ~ "0-50",
  #                                 dband >= 50 & dband < 100 ~ "50-100",
  #                                 dband >= 100 & dband < 150 ~ "100-150",
  #                                 dband >= 150 & dband < 200 ~ "150-200",
  #                                 dband >= 200 & dband < 300 ~ "200-300",
  #                                 TRUE ~ ">300")) %>%
  group_by(species, dband, location) %>%
  summarize(abundance = round(sqrt(sqrt(sum(abu.corr))))) %>%
  spread(species, abundance, fill =  0)

carib.deep.pca <- prcomp(carib.deep.com[-c(1:2)], center = TRUE, scale = TRUE)
summary(carib.deep.pca)

carib.deep.points <- 
  # first convert the pca results to a tibble
  as_tibble(carib.deep.pca$x) %>% 
  bind_cols(carib.deep.com) %>%
  mutate(zone = case_when(dband >= 40 & dband < 80 ~ "Upper Mesophotic",
                          dband >= 80 & dband < 130 ~ "Lower Mesophotic",
                          dband >= 130 & dband < 190 ~ "Upper Rariphotic",
                          dband >= 190 & dband < 310 ~ "Lower Rariphotic",
                          TRUE ~ "Abyss"))

carib.deep.hull <- 
  carib.deep.points %>% 
  group_by(location, zone) %>% 
  slice(chull(PC1, PC2))

carib.deep.load.PC1pos <- 
  as_tibble(carib.deep.pca$rotation, rownames = 'variable') %>%
  top_n(3, PC1) %>%
  rename(species = "variable") %>%
  select(species, PC1, PC2)

carib.deep.load.PC1 <- 
  as_tibble(carib.deep.pca$rotation, rownames = 'variable') %>%
  top_n(-3, PC1) %>%
  rename(species = "variable") %>%
  select(species, PC1, PC2) %>%
  bind_rows(carib.deep.load.PC1pos)

carib.deep.load.PC2pos <- 
  as_tibble(carib.deep.pca$rotation, rownames = 'variable') %>%
  top_n(3, PC2) %>%
  rename(species = "variable") %>%
  select(species, PC1, PC2)

carib.deep.load.PC <- 
  as_tibble(carib.deep.pca$rotation, rownames = 'variable') %>%
  top_n(-3, PC2) %>%
  rename(species = "variable") %>%
  select(species, PC1, PC2) %>%
  bind_rows(carib.deep.load.PC2pos) %>%
  full_join(carib.deep.load.PC1)



pca.plot <- ggplot(carib.deep.points) +
  geom_jitter(data = carib.deep.points, aes(x = PC1, y = PC2, 
                                            color = location, 
                                            fill = location, 
                                            shape = zone), alpha = 0.5, color = "black") +
  geom_polygon(data = carib.deep.hull,
               aes(x = PC1, y = PC2, fill = location,
                   colour = location,
                   group = zone),
               alpha = 0.3,
               show.legend = TRUE) +
  geom_segment(data = carib.deep.load.PC, 
               aes(x = 0, y = 0, 
                   xend = PC1 * 50,
                   yend = PC2 * 50), lwd = 0.1, 
               arrow = arrow(length = unit(1/2, 'picas'), type = "open")) +
  geom_text_repel(data = carib.deep.load.PC, aes(x = PC1*52, 
                                          y = PC2*52, label = species), size = 2, segment.color = "grey69") +
  theme_bw() +
  scale_fill_fish_d(option = "Bodianus_rufus") +
  scale_color_fish_d(option = "Bodianus_rufus") +
  scale_shape_manual(values = c(21:25)) +
  facet_wrap(.~location)

pca.plot



#### full cluster ####
carib.deep.com <- as.data.frame(carib.deep.com)
rownames(carib.deep.com) <- interaction(carib.deep.com$dband, carib.deep.com$location)
carib.com.cl <- carib.deep.com[-c(1:2)]

carib.clust <- simprof(data = carib.com.cl, method.cluster = "ward.D", method.distance = "czekanowski")
carib.clust
carib.clust.veg <- hclust(vegdist(carib.com.cl, method = "bray"), method = "ward.D")
plot(carib.clust.veg)
carib.dend <- as.dendrogram(carib.clust.veg)
carib.dend.data <- dendro_data(carib.dend, type = "rectangle")

carib.deep.com.h <- carib.deep.com %>%
  unite(label, dband, location, remove = F, sep = ".")

carib.dend.helper <- carib.dend.data$labels %>%
  inner_join(carib.deep.com.h) 


carib.dend.data$labels <- carib.dend.helper

carib.dend.plot <- ggplot(carib.dend.data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), color = "grey46")+
  geom_text(data = carib.dend.data$labels, aes(x, y, label = dband, color = location, alpha = dband),
            hjust = 1, angle = 90, size = 3, fontface = "bold") +
  scale_color_fish(option = "Bodianus_rufus", discrete = T) +
  theme_bw() +
  theme(legend.position = "top",
        line = element_blank(),
        title = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
carib.dend.plot

ggsave(carib.dend.plot, file =  "carib.dend.plot.png", width = 10, height = 8)

### split up PCAs by location and cluster analysis

# Curacao
cur.com <- carib.deep.fishbase %>%
  filter(location == "Curacao") %>%
  select(species, dband, location, hab, abundance, abu.corr) %>%
  filter(hab == "deep") %>%
  # mutate(dband.coarse = case_when(dband >= 0 & dband < 50 ~ "0-50",
  #                                 dband >= 50 & dband < 100 ~ "50-100",
  #                                 dband >= 100 & dband < 150 ~ "100-150",
  #                                 dband >= 150 & dband < 200 ~ "150-200",
  #                                 dband >= 200 & dband < 300 ~ "200-300",
  #                                 TRUE ~ ">300")) %>%
  group_by(species, dband) %>%
  summarize(abundance = sqrt(sum(abu.corr))) %>%
  spread(species, abundance, fill =  0) %>%
  as.data.frame()

rownames(cur.com) <- cur.com$dband
cur.com.cl <- cur.com[-1]

cur.clust <- simprof(data = cur.com.cl, method.cluster = "complete", method.distance="braycurtis",
                     method.transform="identity", alpha=0.0000001,
                     sample.orientation="row", const=0,
                     silent=TRUE, increment=100,
                     undef.zero=TRUE, warn.braycurtis=TRUE)
cur.clust

# Curacao: 8 groups, from shallow to deep
cur.res.tib <- tibble(cl1 = list(cur.clust$significantclusters[[1]]),
                      cl2 = list(cur.clust$significantclusters[[2]]),
                      cl3 = list(cur.clust$significantclusters[[3]]),
                      cl4 = list(cur.clust$significantclusters[[4]]),
                      cl5 = list(cur.clust$significantclusters[[5]]),
                      cl6 = list(cur.clust$significantclusters[[6]]),
                      cl7 = list(cur.clust$significantclusters[[7]])) %>%
  unnest_longer(cl1) %>%
  unnest_longer(cl2) %>%
  unnest_longer(cl3) %>%
  unnest_longer(cl4) %>%
  unnest_longer(cl5) %>%
  unnest_longer(cl6) %>%
  unnest_longer(cl7) %>%
  gather(key = "cluster", value = "dband") %>%
  distinct() %>%
  mutate(dband = as.numeric(dband)) %>%
  full_join(cur.com) %>%
  select(cluster, dband) %>%
  mutate(label = as.factor(dband)) %>%
  arrange(dband) %>%
  mutate(run_num = 1:27) %>%
  group_by(cluster) %>%
  mutate(ord_num = min(run_num))

# hclust 
cur.clust.veg <- hclust(vegdist(cur.com.cl, method = "bray"), method = "complete")
cur.dend <- as.dendrogram(cur.clust.veg)
cur.dend.data <- dendro_data(cur.dend, type = "rectangle")

cur.dend.helper <- cur.dend.data$labels %>%
  inner_join(cur.res.tib) 


cur.dend.data$labels <- cur.dend.helper

# get colors for each level using fishualize and tinter package
carib.pal <- fish(4, option = "Bodianus_rufus")
cur.pal <- tinter(carib.pal[1], steps = 5)

cur.plot <- ggplot(cur.dend.data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), color = "grey46")+
  geom_text(data = cur.dend.data$labels, aes(x, y, label = label, color = as.factor(ord_num)),
            hjust = 1, angle = 90, size = 3, fontface = "bold") +
  scale_color_manual(values = cur.pal[-1]) +
  theme_bw() +
  theme(legend.position = "none",
        line = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  ggtitle("Curaçao")
cur.plot


### pca plot
cur.pca <- prcomp(cur.com[-1], center = TRUE, scale = TRUE)
summary(cur.pca)

cur.points <- 
  # first convert the pca results to a tibble
  as_tibble(cur.pca$x) %>% 
  bind_cols(cur.com) %>%
  inner_join(cur.dend.helper)

cur.hull <- 
  cur.points %>% 
  group_by(cluster) %>% 
  slice(chull(PC1, PC2))

cur.simper <- with(cur.points, simper(cur.com[-1], cluster))
cur.simsum <- summary(cur.simper)

cur.simspecies <- data_frame("species" = as.character())

for (i in 1:length(cur.simsum)){
  spec <- as.data.frame(cur.simsum[i]) %>%
    mutate(species = rownames(.)) %>%
    select(1,7) %>%
    filter(. > 0.05) %>%
    select(species)
  
  cur.simspecies <- bind_rows(cur.simspecies, spec)
  
}

cur.species <- cur.simspecies %>%
  distinct()

cur.loads <- 
  as_tibble(cur.pca$rotation, rownames = 'species') %>%
  inner_join(cur.species) %>%
  select(species, PC1, PC2)


cur.pca.plot <- ggplot(cur.points) +
  geom_jitter(data = cur.points, aes(x = PC1, y = PC2, 
                                            color = as.factor(ord_num), 
                                            fill = as.factor(ord_num)), shape = 21, alpha = 0.5, color = "black") +
  geom_polygon(data = cur.hull,
               aes(x = PC1, y = PC2, fill = as.factor(ord_num),
                   colour = as.factor(ord_num),
                   group = as.factor(ord_num)),
               alpha = 0.3,
               show.legend = TRUE) +
  geom_segment(data = cur.loads, 
               aes(x = 0, y = 0, 
                   xend = PC1*30,
                   yend = PC2*30), lwd = 0.1, 
               arrow = arrow(length = unit(1/2, 'picas'), type = "open")) +
  geom_text_repel(data = cur.loads, aes(x = PC1*30, 
                                                 y = PC2*30, label = species), size = 2, segment.color = "grey69") +
  theme_bw() +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_blank()) +
  scale_fill_manual(values = cur.pal[-1]) +
  scale_color_manual(values = cur.pal[-1]) +
  guides(colour = guide_legend(nrow = 1),
         fill = guide_legend(nrow = 1))
cur.pca.plot



# Bonaire
bon.com <- carib.deep.fishbase %>%
  filter(location == "Bonaire") %>%
  select(species, dband, location, hab, abundance, abu.corr) %>%
  mutate(dband.grouped = case_when(dband >= 190 & dband <= 200 ~ 190,
                                  dband >= 220 & dband <= 250 ~ 220,
                                  dband >= 280 & dband <= 300 ~ 280,
                                  TRUE ~ dband)) %>%
  group_by(species, dband.grouped) %>%
  summarize(abundance = sqrt(sum(abundance))) %>%
  spread(species, abundance, fill =  0) %>%
  as.data.frame()


rownames(bon.com) <- bon.com$dband.grouped
bon.com.cl <- bon.com[-1]

bon.clust <- simprof(data = bon.com.cl, method.cluster = "complete", method.distance="braycurtis",
                     method.transform="identity", alpha=0.0000001,
                     sample.orientation="row", const=0,
                     silent=TRUE, increment=100,
                     undef.zero=TRUE, warn.braycurtis=TRUE)
bon.clust
# Bonaire: 7 groups, from shallow to deep
bon.res.tib <- tibble(cl1 = list(bon.clust$significantclusters[[1]]),
                      cl2 = list(bon.clust$significantclusters[[2]]),
                      cl3 = list(bon.clust$significantclusters[[3]]),
                      cl4 = list(bon.clust$significantclusters[[4]]),
                      cl5 = list(bon.clust$significantclusters[[5]]),
                      cl6 = list(bon.clust$significantclusters[[6]]),
                      cl7 = list(bon.clust$significantclusters[[7]])) %>%
  unnest_longer(cl1) %>%
  unnest_longer(cl2) %>%
  unnest_longer(cl3) %>%
  unnest_longer(cl4) %>%
  unnest_longer(cl5) %>%
  unnest_longer(cl6) %>%
  unnest_longer(cl7) %>%
  gather(key = "cluster", value = "dband.grouped") %>%
  distinct() %>%
  mutate(dband.grouped = as.numeric(dband.grouped)) %>%
  full_join(bon.com) %>%
  select(cluster, dband.grouped) %>%
  mutate(label = as.factor(dband.grouped)) %>%
  arrange(dband.grouped) %>%
  mutate(run_num = 1:21) %>%
  group_by(cluster) %>%
  mutate(ord_num = min(run_num))

# hclust 
bon.clust.veg <- hclust(vegdist(bon.com.cl, method = "bray"), method = "complete")
bon.dend <- as.dendrogram(bon.clust.veg)
bon.dend.data <- dendro_data(bon.dend, type = "rectangle")

bon.dend.helper <- bon.dend.data$labels %>%
  inner_join(bon.res.tib) 


bon.dend.data$labels <- bon.dend.helper

# get colors for each level using fishualize and tinter package
bon.pal <- tinter(carib.pal[2], steps = 5)

bon.plot <- ggplot(bon.dend.data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), color = "grey46")+
  geom_text(data = bon.dend.data$labels, aes(x, y, label = label, color = as.factor(ord_num)),
            hjust = 1, angle = 90, size = 3, fontface = "bold") +
  scale_color_manual(values = bon.pal[-1]) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(legend.position = "none",
        line = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  ggtitle("Bonaire")
bon.plot

### Bonaire pca plot
bon.pca <- prcomp(bon.com[-1], center = TRUE, scale = TRUE)
summary(bon.pca)

bon.points <- 
  # first convert the pca results to a tibble
  as_tibble(bon.pca$x) %>% 
  bind_cols(bon.com) %>%
  inner_join(bon.dend.helper)

bon.hull <- 
  bon.points %>% 
  group_by(cluster) %>% 
  slice(chull(PC1, PC2))

bon.simper <- with(bon.points, simper(bon.com[-1], cluster))
bon.simsum <- summary(bon.simper)

bon.simspecies <- data_frame("species" = as.character())

for (i in 1:length(bon.simsum)){
  spec <- as.data.frame(bon.simsum[i]) %>%
    mutate(species = rownames(.)) %>%
    select(1,7) %>%
    filter(. > 0.05) %>%
    select(species)
  
  bon.simspecies <- bind_rows(bon.simspecies, spec)
  
}

bon.species <- bon.simspecies %>%
  distinct()

bon.loads <- 
  as_tibble(bon.pca$rotation, rownames = 'species') %>%
  inner_join(bon.species) %>%
  select(species, PC1, PC2)


bon.pca.plot <- ggplot(bon.points) +
  geom_jitter(data = bon.points, aes(x = PC1, y = PC2, 
                                     color = as.factor(ord_num), 
                                     fill = as.factor(ord_num)), shape = 21, alpha = 0.5, color = "black") +
  geom_polygon(data = bon.hull,
               aes(x = PC1, y = PC2, fill = as.factor(ord_num),
                   colour = as.factor(ord_num),
                   group = as.factor(ord_num)),
               alpha = 0.3,
               show.legend = TRUE) +
  geom_segment(data = bon.loads, 
               aes(x = 0, y = 0, 
                   xend = PC1*30,
                   yend = PC2*30), lwd = 0.1, 
               arrow = arrow(length = unit(1/2, 'picas'), type = "open")) +
  geom_text_repel(data = bon.loads, aes(x = PC1*30, 
                                        y = PC2*30, label = species), size = 2, segment.color = "grey69") +
  theme_bw() +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_blank()) +
  scale_fill_manual(values = bon.pal) +
  scale_color_manual(values = bon.pal) +
  guides(colour = guide_legend(nrow = 1),
         fill = guide_legend(nrow = 1))
bon.pca.plot




# Roatan
roa.com <- carib.deep.fishbase %>%
  filter(location == "Roatan") %>%
  filter(dband >= 40 & dband <= 300) %>%
  select(species, dband, location, hab, abundance, abu.corr) %>%
  mutate(dband.grouped = case_when(dband >= 200 & dband <= 210 ~ 200,
                                   dband >= 180 & dband <= 190 ~ 180,
                                   dband >= 230 & dband <= 240 ~ 230,
                                   dband >= 260 & dband <= 290 ~ 260,
                                                            TRUE ~ dband)) %>%
  group_by(species, dband.grouped) %>%
  summarize(abundance = sqrt(sum(abu.corr))) %>%
  spread(species, abundance, fill =  0) %>%
  as.data.frame()
  
rownames(roa.com) <- roa.com$dband.grouped
roa.com.cl <- roa.com[-1]

roa.clust <- simprof(data = roa.com.cl, method.cluster = "complete", method.distance="braycurtis",
                     method.transform="identity", alpha=0.0000001,
                     sample.orientation="row", const=0,
                     silent=TRUE, increment=100,
                     undef.zero=TRUE, warn.braycurtis=TRUE)
roa.clust

# Roatan: 5 groups
roa.res.tib <- tibble(cl1 = list(roa.clust$significantclusters[[1]]),
                      cl2 = list(roa.clust$significantclusters[[2]]),
                      cl3 = list(roa.clust$significantclusters[[3]]),
                      cl4 = list(roa.clust$significantclusters[[4]]),
                      cl5 = list(roa.clust$significantclusters[[5]])) %>%
  unnest_longer(cl1) %>%
  unnest_longer(cl2) %>%
  unnest_longer(cl3) %>%
  unnest_longer(cl4) %>%
  unnest_longer(cl5) %>%
  gather(key = "cluster", value = "dband.grouped") %>%
  distinct() %>%
  mutate(dband.grouped = as.numeric(dband.grouped))%>%
  full_join(roa.com) %>%
  select(cluster, dband.grouped) %>%
  mutate(label = as.numeric(dband.grouped)) %>%
  arrange(dband.grouped) %>%
  mutate(run_num = 1:21) %>%
  group_by(cluster) %>%
  mutate(ord_num = min(run_num))

# hclust 
roa.clust.veg <- hclust(vegdist(roa.com.cl, method = "bray"), method = "complete")
roa.dend <- as.dendrogram(roa.clust.veg)
roa.dend.data <- dendro_data(roa.dend, type = "rectangle")

roa.dend.helper <- roa.dend.data$labels %>%
  mutate(label = as.numeric(label))%>%
  inner_join(roa.res.tib)
  

roa.dend.data$labels <- roa.dend.helper

# get colors for each level using fishualize and tinter package
carib.pal <- fish(4, option = "Bodianus_rufus")
roa.pal <- tinter(carib.pal[3], steps = 4)

roa.plot <- ggplot(roa.dend.data$segments, horiz = TRUE) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), color = "grey46")+
  geom_text(data = roa.dend.data$labels, aes(x, y, label = label, color = as.factor(ord_num)),
            hjust = 1, angle = 90, size = 3, fontface = "bold") +
  scale_color_manual(values = roa.pal[-1]) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(legend.position = "none",
        line = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  ggtitle("Roatan")
roa.plot


### Roatan pca plot
roa.pca <- prcomp(roa.com[-1], center = TRUE, scale = TRUE)
summary(roa.pca)

roa.points <- 
  # first convert the pca results to a tibble
  as_tibble(roa.pca$x) %>% 
  bind_cols(roa.com) %>%
  inner_join(roa.dend.helper)

roa.hull <- 
  roa.points %>% 
  group_by(cluster) %>% 
  slice(chull(PC1, PC2))

roa.simper <- with(roa.points, simper(roa.com[-1], cluster))
roa.simsum <- summary(roa.simper)

roa.simspecies <- data_frame("species" = as.character())

for (i in 1:length(roa.simsum)){
  spec <- as.data.frame(roa.simsum[i]) %>%
    mutate(species = rownames(.)) %>%
    select(1,7) %>%
    filter(. > 0.05) %>%
    select(species)
  
  roa.simspecies <- bind_rows(roa.simspecies, spec)
  
}

roa.species <- roa.simspecies %>%
  distinct()

roa.loads <- 
  as_tibble(roa.pca$rotation, rownames = 'species') %>%
  inner_join(roa.species) %>%
  select(species, PC1, PC2)


roa.pca.plot <- ggplot(roa.points) +
  geom_jitter(data = roa.points, aes(x = PC1, y = PC2, 
                                     color = as.factor(ord_num), 
                                     fill = as.factor(ord_num)), shape = 21, alpha = 0.5, color = "black") +
  geom_polygon(data = roa.hull,
               aes(x = PC1, y = PC2, fill = as.factor(ord_num),
                   colour = as.factor(ord_num),
                   group = as.factor(ord_num)),
               alpha = 0.3,
               show.legend = TRUE) +
  geom_segment(data = roa.loads, 
               aes(x = 0, y = 0, 
                   xend = PC1*30,
                   yend = PC2*30), lwd = 0.1, 
               arrow = arrow(length = unit(1/2, 'picas'), type = "open")) +
  geom_text_repel(data = roa.loads, aes(x = PC1*30, 
                                        y = PC2*30, label = species), size = 2, segment.color = "grey69") +
  theme_bw() +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_blank()) +
  scale_fill_manual(values = roa.pal) +
  scale_color_manual(values = roa.pal) +
  guides(colour = guide_legend(nrow = 1),
         fill = guide_legend(nrow = 1))
roa.pca.plot



# St. Eustatius
sta.com <- carib.deep.fishbase %>%
  filter(location == "St. Eustatius") %>%
  select(species, dband, location, hab, abundance, abu.corr) %>%
  filter(hab == "deep") %>%
  filter(dband <= 300) %>%
  mutate(dband.grouped = case_when(dband >= 240 & dband <= 270 ~ 240,
                                   TRUE ~ dband)) %>%
  group_by(species, dband.grouped) %>%
  summarize(abundance = sqrt(sum(abu.corr))) %>%
  spread(species, abundance, fill =  0) %>%
  as.data.frame()

rownames(sta.com) <- sta.com$dband.grouped
sta.com.cl <- sta.com[-1]

sta.clust <- simprof(data = sta.com.cl, method.cluster = "complete", method.distance="braycurtis",
                     method.transform="identity", alpha=0.0000001,
                     sample.orientation="row", const=0,
                     silent=TRUE, increment=100,
                     undef.zero=TRUE, warn.braycurtis=TRUE)
sta.clust

# 8 groups
sta.res.tib <- tibble(cl1 = list(sta.clust$significantclusters[[1]]),
                      cl2 = list(sta.clust$significantclusters[[2]]),
                      cl3 = list(sta.clust$significantclusters[[3]]),
                      cl4 = list(sta.clust$significantclusters[[4]])) %>%
  unnest_longer(cl1) %>%
  unnest_longer(cl2) %>%
  unnest_longer(cl3) %>%
  unnest_longer(cl4) %>%
  gather(key = "cluster", value = "dband.grouped") %>%
  distinct() %>%
  mutate(dband.grouped = as.numeric(dband.grouped)) %>%
  full_join(sta.com) %>%
  select(cluster, dband.grouped) %>%
  mutate(label = as.numeric(dband.grouped)) %>%
  arrange(dband.grouped) %>%
  mutate(run_num = 1:21) %>%
  group_by(cluster) %>%
  mutate(ord_num = min(run_num))

# hclust 
sta.clust.veg <- hclust(vegdist(sta.com.cl, method = "bray"), method = "complete")
sta.dend <- as.dendrogram(sta.clust.veg)
sta.dend.data <- dendro_data(sta.dend, type = "rectangle")

sta.dend.helper <- sta.dend.data$labels %>%
  mutate(label = as.numeric(label)) %>%
  inner_join(sta.res.tib)

sta.dend.data$labels <- sta.dend.helper


# get colors for each level using fishualize and tinter package
carib.pal <- fish(4, option = "Bodianus_rufus")
sta.pal <- tinter(carib.pal[4], steps = 3)


sta.plot <- ggplot(sta.dend.data$segments, horiz = TRUE) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), color = "grey46")+
  geom_text(data = sta.dend.data$labels, aes(x, y, label = label, color = as.factor(ord_num)),
            hjust = 1, angle = 90, size = 3, fontface = "bold") +
  scale_color_manual(values = sta.pal[-1]) +
  theme_bw() +
  theme(legend.position = "none",
        line = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  ggtitle("St. Eustatius")
sta.plot

### Statia pca plot
sta.pca <- prcomp(sta.com[-1], center = TRUE, scale = TRUE)
summary(sta.pca)

sta.points <- 
  # first convert the pca results to a tibble
  as_tibble(sta.pca$x) %>% 
  bind_cols(sta.com) %>%
  inner_join(sta.dend.helper)

sta.hull <- 
  sta.points %>% 
  group_by(cluster) %>% 
  slice(chull(PC1, PC2))

sta.simper <- with(sta.points, simper(sta.com[-1], cluster))
sta.simsum <- summary(sta.simper)

sta.simspecies <- data_frame("species" = as.character())

for (i in 1:length(sta.simsum)){
  spec <- as.data.frame(sta.simsum[i]) %>%
    mutate(species = rownames(.)) %>%
    select(1,7) %>%
    filter(. > 0.05) %>%
    select(species)
  
  sta.simspecies <- bind_rows(sta.simspecies, spec)
  
}

sta.species <- sta.simspecies %>%
  distinct()

sta.loads <- 
  as_tibble(sta.pca$rotation, rownames = 'species') %>%
  inner_join(sta.species) %>%
  select(species, PC1, PC2)


sta.pca.plot <- ggplot(sta.points) +
  geom_jitter(data = sta.points, aes(x = PC1, y = PC2, 
                                     color = as.factor(ord_num), 
                                     fill = as.factor(ord_num)), shape = 21, alpha = 0.5, color = "black") +
  geom_polygon(data = sta.hull,
               aes(x = PC1, y = PC2, fill = as.factor(ord_num),
                   colour = as.factor(ord_num),
                   group = as.factor(ord_num)),
               alpha = 0.3,
               show.legend = TRUE) +
  geom_segment(data = sta.loads,
               aes(x = 0, y = 0,
                   xend = PC1*30,
                   yend = PC2*30), lwd = 0.1,
               arrow = arrow(length = unit(1/2, 'picas'), type = "open")) +
  geom_text_repel(data = sta.loads, aes(x = PC1*30, 
                                        y = PC2*30, label = species), size = 2, segment.color = "grey69") +
  theme_bw() +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_blank()) +
  scale_fill_manual(values = sta.pal) +
  scale_color_manual(values = sta.pal) +
  guides(colour = guide_legend(nrow = 1),
         fill = guide_legend(nrow = 1))
sta.pca.plot

cluster.plots <- (cur.plot/cur.pca.plot | bon.plot/bon.pca.plot | roa.plot/roa.pca.plot | sta.plot/sta.pca.plot)
ggsave(cluster.plots, file = "cluster.plots.png", width = 18, height = 11)


### check Ross' data
sfdata <- read.csv(file = "SFGC.csv") %>%
  filter(dband >= 40) %>%
  filter(dband < 300)%>%
  select(which(!colSums(sfdata, na.rm=TRUE) %in% 0))

rownames(sfdata) <- sfdata$dband
sf.com.cl <- sfdata[-1]

sf.clust <- simprof(data = sf.com.cl, method.cluster = "complete", method.distance = "binary")
simprof.plot(sf.clust) 


### check Semmler's data
semmlerdata <- read.csv(file = "BioGoMX_r.csv") %>%
  filter(dband >= 40) %>%
  filter(dband < 300)%>%
  select(which(!colSums(sfdata, na.rm=TRUE) %in% 0))

rownames(semmlerdata) <- semmlerdata$dband
semmlerdata.cl <- semmlerdata[-1]

semmler.clust <- simprof(data = semmlerdata.cl, method.cluster = "complete", method.distance = "binary")
simprof.plot(semmler.clust) 



# 5 groups
sf.res.tib <- tibble(cl1 = list(sta.clust$significantclusters[[5]]),
                      cl2 = list(sta.clust$significantclusters[[6]]),
                      cl3 = list(sta.clust$significantclusters[[4]]),
                      cl4 = list(sta.clust$significantclusters[[3]]),
                      cl5 = list(sta.clust$significantclusters[[1]]),
                      cl6 = list(sta.clust$significantclusters[[2]])) %>%
  unnest_longer(cl1) %>%
  unnest_longer(cl2) %>%
  unnest_longer(cl3) %>%
  unnest_longer(cl4) %>%
  unnest_longer(cl5) %>%
  unnest_longer(cl6) %>%
  gather(key = "cluster", value = "dband") %>%
  distinct() %>%
  mutate(dband = as.numeric(dband)) %>%
  full_join(sta.com) %>%
  select(cluster, dband) %>%
  mutate(label = as.factor(dband))

# hclust 
sf.clust.veg <- hclust(vegdist(sf.com.cl, method = "jaccard"), method = "ward.D")
sf.dend <- as.dendrogram(sf.clust.veg)
sf.dend.data <- dendro_data(sf.dend, type = "rectangle")

sf.dend.data$labels$label <- as.numeric(sf.dend.data$labels$label)


# get colors for each level using fishualize and tinter package
carib.pal <- fish(4, option = "Bodianus_rufus")
sta.pal <- tinter(carib.pal[4], steps = 4)


sf.plot <- ggplot(sf.dend.data$segments, horiz = TRUE) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), color = "grey46")+
  geom_text(data = sf.dend.data$labels, aes(x, y, label = label, color = label),
            hjust = 1, angle = 90, size = 3, fontface = "bold") +
  theme_bw() +
  theme(legend.position = "none",
        line = element_blank(),
        title = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_color_fish(option = "Prionace_glauca", direction = -1)
sf.plot
ggsave(sf.plot, file = "caribbean.clustering.png", width = 8, height = 6)

### try some comps


%>%
  untangle(method = "step1side") %>% # Find the best alignment layout
  tanglegram()        
cor.dendlist(dend_list, method = "cophenetic")
cor_bakers_gamma(cur.dend, roa.dend)

sta.com$dband
bon.com$dband
roa.com$dband
cur.com.corr$dband

cur.com.corr <- cur.com %>%
  filter(dband < 280) %>%
  filter(dband %not_in% c(240, 250))
c.cur <- hclust(vegdist(cur.com.corr[-1], method = "bray"), method = "ward.D")
cur.dend.comp <- as.dendrogram(c.cur)

bon.com.corr <- bon.com %>%
  filter(dband < 280) %>%
  filter(dband %not_in% c(240, 250))
c.bon <- hclust(vegdist(bon.com.corr[-1], method = "bray"), method = "ward.D")
bon.dend.comp <- as.dendrogram(c.bon)

roa.com.corr <- roa.com %>%
  filter(dband < 280) %>%
  filter(dband %not_in% c(240, 250))
c.roa <- hclust(vegdist(roa.com.corr[-1], method = "bray"), method = "ward.D")
roa.dend.comp <- as.dendrogram(c.roa)

sta.com.corr <- sta.com %>%
  filter(dband < 280) %>%
  filter(dband %not_in% c(240, 250))
c.sta <- hclust(vegdist(sta.com.corr[-1], method = "bray"), method = "ward.D")
sta.dend.comp <- as.dendrogram(c.sta)

sf.com.corr <- sfdata %>%
  filter(dband < 280) %>%
  filter(dband %not_in% c(240, 250))
rownames(sf.com.corr) <- sf.com.corr$dband

c.sf <- hclust(vegdist(sf.com.corr[-1], method = "jaccard"), method = "ward.D")
sf.dend.comp <- as.dendrogram(c.sf)

dend_list <- dendlist("Curacao" = cur.dend.comp, 
                      "Bonaire" = bon.dend.comp, 
                      "Roatan" = roa.dend.comp,
                      "Statia" = sta.dend.comp,
                      "All" = sf.dend.comp)
den.corrs <- cor.dendlist(dend_list, method = "cophenetic") %>%
  reshape::melt() %>%
  rename(Loc1 = "X1", Loc2 = "X2", Correlation = "value") %>%
  distinct()

correlation.plot <- ggplot(data = den.corrs, aes(x=Loc1, y=Loc2, fill=Correlation)) + 
  geom_tile(alpha = 0.75, color = "grey23") +
  geom_text(aes(label = round(Correlation, 2)), fontface = "italic") +
  theme_bw() +
  scale_fill_gradientn(limits = c(0,1),
                       colours=c(fish(256, option = "Hypsypops_rubicundus", end = 0.8))) +
  theme(axis.title = element_blank(),
        axis.text = element_text(color = "black"))
correlation.plot
ggsave(correlation.plot, file = "dendrogram.correlations.png", width = 8, height = 7)

dendlist("Curacao" = cur.dend.comp,
                      "Statia" =sta.dend.comp) %>%
  untangle(method = "step1side") %>% 
  tanglegram(
    highlight_distinct_edges = F, # Turn-off dashed lines
    common_subtrees_color_lines = F, # Turn-off line colors
    common_subtrees_color_branches = TRUE)


cur.bon.comp <- cur.com %>%
  filter(dband < 300) %>%
  filter(dband %not_in% c(240))
c.cur.bon <- hclust(vegdist(cur.bon.comp[-1], method = "bray"), method = "ward.D")
cur.dend.bon <- as.dendrogram(c.cur.bon)

bon.cur.comp <- bon.com %>%
  filter(dband < 300) %>%
  filter(dband %not_in% c(240))
c.bon.cur <- hclust(vegdist(bon.cur.comp[-1], method = "bray"), method = "ward.D")
bon.dend.cur <- as.dendrogram(c.bon.cur)

dend.list.cur.bon <- dendlist("Curacao" = cur.dend.bon, 
                      "Statia" = cur.dend.cur) %>%
  untangle(method = "step1side") %>% 
  tanglegram(
    highlight_distinct_edges = F, # Turn-off dashed lines
    common_subtrees_color_lines = F, # Turn-off line colors
    common_subtrees_color_branches = TRUE)


cur.bon <- dendlist(dend_list$dend1, dend_list$dend2) %>%
  untangle(method = "step1side") %>% 
  tanglegram(
    highlight_distinct_edges = FALSE, # Turn-off dashed lines
    common_subtrees_color_lines = T, # Turn-off line colors
    common_subtrees_color_branches = TRUE)



#################################
#################################
#### trait based analysis #######

# get species names metadata 
traits.species = carib.deep.traits %>%
  select(species, species.fb) %>%
  distinct() %>%
  mutate(species.traits = case_when(species.fb == "Centropristis fuscula" ~ "Centropristis philadelphica",
                                    species.fb == "Chriolepis zebra" ~ "Pinnichthys aimoriensis",
                                    species.fb == "Chromis enchrysura" ~ "Chromis enchrysurus",
                                    species.fb == "Fistularia commersonii" ~ "Fistularia tabacaria",
                                    species.fb == "Haemulon chrysargyreum" ~ "Brachygenys chrysargyrea",
                                    species.fb == "Lipogramma haberi" ~ "Lipogramma haberorum",
                                    species.fb == "Polylepion russelli" ~ "Polylepion species A",
                                    species.fb == "Sargocentron coruscum" ~ "Neoniphon coruscum",
                                    species.fb == "Stegastes variabilis" ~ "Stegastes xanthurus",
                                    TRUE ~ species.fb))



t.trophic <- read.csv("traits_trophic.csv") %>%
  rename_all( ~ paste0("trophic_", .)) %>%
  rename(species.traits = 1) %>%
  mutate(across(2:ncol(.), ~as.numeric(recode(., "x" = 1, .default = 0)))) %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>%
  filter(sum > 0) %>%
  # mutate(across(2:ncol(.), ~./sum)) %>%
  select(-sum) %>%
  right_join(traits.species)

t.diet <- read.csv("traits_diet.csv") %>%
  rename_all( ~ paste0("diet_", .)) %>%
  rename(species.traits = 1) %>%
  mutate(diet_MobileBenthicInverts = case_when(diet_Mobile.benthic.crustacea..shrimps.crabs. == "x" ~ "x",
                              diet_Mobile.benthic.gastropods.bivalves == "x" ~ "x",
                              diet_Mobile.benthic.worms == "x" ~ "x",
                              diet_Sea.stars.cucumbers.urchins == "x" ~ "x",
                              diet_Octopus.squid.cuttlefish == "x" ~ "x",
                                        TRUE ~ ""),
         diet_Zooplankton = case_when(diet_Pelagic.crustacea == "x" ~ "x",
                                      diet_Pelagic.fish.eggs == "x" ~ "x",
                                      diet_Pelagic.fish.larvae == "x" ~ "x",
                                      diet_Pelagic.jellyfish.ctenophores == "x" ~ "x",
                                      diet_Zooplankton == "x" ~ "x",
                                      TRUE ~ ""),
         diet_SessileInverts = case_when(diet_Sessile.crustacea == "x" ~ "x",
                                         diet_Sessile.molluscs == "x" ~ "x",
                                         diet_Sessile.worms == "x" ~ "x",
                                         diet_Soft.corals.hydroids == "x" ~ "x",
                                         diet_Sponges.seasquirts.bryozoa == "x" ~ "x",
                                      TRUE ~ ""),
         diet_Megafauna = case_when(diet_Sharks.rays == "x" ~ "x",
                                    diet_Sea.snakes.mammals.turtles.birds == "x" ~ "x",
                                    TRUE ~ "")) %>%
  select(-diet_Mobile.benthic.crustacea..shrimps.crabs.,
         -diet_Mobile.benthic.gastropods.bivalves,
         -diet_Mobile.benthic.worms,
         -diet_Sea.stars.cucumbers.urchins,
         -diet_Octopus.squid.cuttlefish,
         -diet_Pelagic.crustacea,
         -diet_Pelagic.fish.eggs,
         -diet_Pelagic.fish.larvae,
         -diet_Pelagic.jellyfish.ctenophores,
         -diet_Sessile.crustacea,
         -diet_Sessile.molluscs,
         -diet_Sessile.worms,
         -diet_Soft.corals.hydroids,
         -diet_Sponges.seasquirts.bryozoa,
         -diet_Sharks.rays,
         -diet_Sea.snakes.mammals.turtles.birds) %>%
  mutate(across(2:ncol(.), ~as.numeric(recode(., "x" = 1, .default = 0)))) %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>%
  filter(sum > 0) %>%
  # mutate(across(2:ncol(.), ~./sum)) %>%
  select(-sum) %>%
  right_join(traits.species)

t.position <- read.csv("traits_position.csv") %>%
  rename_all( ~ paste0("position_", .)) %>%
  rename(species.traits = 1) %>%
  mutate(across(2:ncol(.), ~as.numeric(recode(., "x" = 1, .default = 0)))) %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>%
  filter(sum > 0) %>%
  # mutate(across(2:ncol(.), ~./sum)) %>%
  select(-sum) %>%
  right_join(traits.species)

t.eggs <- read.csv("traits_eggs.csv") %>%
  rename_all( ~ paste0("eggs_", .)) %>%
  rename(species.traits = 1) %>%
  mutate(across(2:ncol(.), ~as.numeric(recode(., "x" = 1, .default = 0)))) %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>%
  filter(sum > 0) %>%
  # mutate(across(2:ncol(.), ~./sum)) %>%
  select(-sum) %>%
  right_join(traits.species)

t.repro <- read.csv("traits_eggs.csv") %>%
  rename_all( ~ paste0("repro_", .)) %>%
  rename(species.traits = 1) %>%
  mutate(repro_ParentalCare = case_when(repro_Brooded == "x" ~ "x",
                                        repro_Live.birth == "x" ~ "x",
                                        TRUE ~ "")) %>%
  select(-repro_Brooded, -repro_Live.birth) %>%
  mutate(across(2:ncol(.), ~as.numeric(recode(., "x" = 1, .default = 0)))) %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>%
  filter(sum > 0) %>%
  # mutate(across(2:ncol(.), ~./sum)) %>%
  select(-sum) %>%
  right_join(traits.species)

t.size <- read.csv("traits_size.csv") %>%
  rename_all( ~ paste0("size_", .)) %>%
  rename(species.traits = 1) %>%
  drop_na(.) %>%
  right_join(traits.species)

t.size.cat <- read.csv("traits_size.csv") %>%
  rename_all( ~ paste0("size_", .)) %>%
  rename(species.traits = 1) %>%
  drop_na(.) %>%
  mutate(size_xs = case_when(size_Size.cm <= 5 ~ 1,
                        TRUE ~ 0),
         size_s = case_when(size_Size.cm > 5 & size_Size.cm <= 10 ~ 1,
                       TRUE ~ 0),
         size_m = case_when(size_Size.cm > 10 & size_Size.cm <= 20 ~ 1,
                       TRUE ~ 0),
         size_l = case_when(size_Size.cm > 20 & size_Size.cm <= 40 ~ 1,
                       TRUE ~ 0),
         size_xl = case_when(size_Size.cm > 40 ~ 1,
                       TRUE ~ 0)) %>%
  select(-size_Size.cm) %>%
  right_join(traits.species)

traits.all <- t.size.cat %>%
  inner_join(t.trophic) %>%
  inner_join(t.diet) %>%
  inner_join(t.position) %>%
  inner_join(t.repro) %>%
  relocate(species.fb)


trait.weights <- c(rep((1/(ncol(t.size.cat)-3)),ncol(t.size.cat)-3), 
                   rep((1/(ncol(t.trophic)-3)),ncol(t.trophic)-3),
                   rep((1/(ncol(t.diet)-3)),ncol(t.diet)-3),
                   rep((1/(ncol(t.position)-3)),ncol(t.position)-3),
                   rep((1/(ncol(t.repro)-3)),ncol(t.repro)-3))

### FD Curacao
cur.names <- data.frame("species" = colnames(cur.com.cl)) 
cur.traits <- cur.names %>%
  left_join(traits.all) %>%
  select(-species.fb, -species.traits) %>%
  filter(species != "Sphyraenops bairdianus")

rownames(cur.traits) = cur.traits$species
cur.traits.rn <- cur.traits[-1]

cur.com.fd <- cur.com.cl %>%
  select(-`Sphyraenops bairdianus`)

cur.namecheck <- data.frame("colnames" = colnames(cur.com.fd), "rownames" = cur.traits$species) %>%
  mutate(namecheck = case_when(colnames == rownames ~ "TRUE",
                               TRUE ~ "FALSE"))

fd.cur <- dbFD(cur.traits.rn, cur.com.fd, w = trait.weights, w.abun = TRUE, m = 2, CWM.type = "all", corr = "cailliez", stand.FRic = TRUE)
fd.cur
cur.fd.res <- data.frame("dband" = as.numeric(rownames(cur.com.fd)),
                         "location" = "Curacao",
                         "spric" = fd.cur$nbsp,
                         "fric" = fd.cur$FRic,
                         "feve" = fd.cur$FEve,
                         "raoq" = fd.cur$RaoQ)

### FD Bonaire
bon.names <- data.frame("species" = colnames(bon.com.cl)) 
bon.traits <- bon.names %>%
  left_join(traits.all) %>%
  select(-species.fb, -species.traits) %>%
  filter(species != "Carangoides bartholomaei")

rownames(bon.traits) = bon.traits$species
bon.traits.rn <- bon.traits[-1]

bon.com.fd <- bon.com.cl %>%
  select(-`Carangoides bartholomaei`)

bon.namecheck <- data.frame("colnames" = colnames(bon.com.fd), "rownames" = bon.traits$species) %>%
  mutate(namecheck = case_when(colnames == rownames ~ "TRUE",
                               TRUE ~ "FALSE"))

fd.bon <- dbFD(bon.traits.rn, bon.com.fd, w = trait.weights, w.abun = TRUE, CWM.type = "all", m = 2, corr = "cailliez", stand.FRic = TRUE)
fd.bon

bon.fd.res <- data.frame("dband" = as.numeric(rownames(bon.com.fd)),
                         "location" = "Bonaire",
                         "spric" = fd.bon$nbsp,
                         "fric" = fd.bon$FRic,
                         "feve" = fd.bon$FEve,
                         "raoq" = fd.bon$RaoQ)


### FD Roatan
roa.names <- data.frame("species" = colnames(roa.com.cl)) 
roa.traits <- roa.names %>%
  left_join(traits.all) %>%
  select(-species.fb, -species.traits) %>%
  filter(!species %in% c("Epigonidae", "Neoscopelus macrolepidotus"))

rownames(roa.traits) = roa.traits$species
roa.traits.rn <- roa.traits[-1]

roa.com.fd <- roa.com.cl %>%
  select(-Epigonidae) %>%
  filter(rowSums(.) != 0)

rowSums(roa.com.fd)
roa.namecheck <- data.frame("colnames" = colnames(roa.com.fd), "rownames" = roa.traits$species) %>%
  mutate(namecheck = case_when(colnames == rownames ~ "TRUE",
                               TRUE ~ "FALSE"))

fd.roa <- dbFD(roa.traits.rn, roa.com.fd, w = trait.weights, w.abun = TRUE, m = 3, CWM.type = "all", corr = "cailliez", stand.FRic = TRUE)
fd.roa

roa.fd.res <- data.frame("dband" = as.numeric(rownames(roa.com.fd)),
                         "location" = "Roatan",
                         "spric" = fd.roa$nbsp,
                         "fric" = fd.roa$FRic,
                         "feve" = fd.roa$FEve,
                         "raoq" = fd.roa$RaoQ)

### FD Statia
sta.names <- data.frame("species" = colnames(sta.com.cl)) 
sta.traits <- sta.names %>%
  left_join(traits.all) %>%
  select(-species.fb, -species.traits) 

rownames(sta.traits) = sta.traits$species
sta.traits.rn <- sta.traits[-1]

sta.com.fd <- sta.com.cl 

sta.namecheck <- data.frame("colnames" = colnames(sta.com.fd), "rownames" = sta.traits$species) %>%
  mutate(namecheck = case_when(colnames == rownames ~ "TRUE",
                               TRUE ~ "FALSE"))

fd.sta <- dbFD(sta.traits.rn, sta.com.fd, w = trait.weights, w.abun = TRUE, m = 3, CWM.type = "all", corr = "cailliez", 
               stand.FRic = TRUE)
fd.sta

sta.fd.res <- data.frame("dband" = as.numeric(rownames(sta.com.fd)),
                         "location" = "Statia",
                         "spric" = fd.sta$nbsp,
                         "fric" = fd.sta$FRic,
                         "feve" = fd.sta$FEve,
                         "raoq" = fd.sta$RaoQ)


fd.res.all <- bind_rows(cur.fd.res, bon.fd.res, roa.fd.res, sta.fd.res) %>%
  drop_na(fric) %>%
  mutate(location = as.factor(location))
  
GGally::ggpairs(fd.res.all)

### SPRIC GAM
spric.gam = gam(spric ~ s(dband) + location, family = poisson, data = fd.res.all)
summary(spric.gam)

spric.newdat <- fd.res.all %>%
  data_grid(dband = seq(40,300,1),
            location = fd.res.all$location)
spric.pred.gam <- data.frame(predict(spric.gam, spric.newdat, se.fit = TRUE, type = "link")) %>%
  add_column(dband = spric.newdat$dband) %>%
  add_column(location = spric.newdat$location) %>%
  mutate(pred.spric = exp(fit), pred.spric.se = exp(se.fit)) %>%
  mutate(lci = pred.spric-1.96*pred.spric.se,
         uci = pred.spric+1.96*pred.spric.se)

spric.plot <- ggplot(spric.pred.gam, aes(x = dband, y = pred.spric, color = location)) +
  geom_line(aes(group = location, color = location)) +
  geom_ribbon(aes(ymin = lci, ymax = uci, group = location, fill = location), alpha = 0.1) +
  geom_point(data = fd.res.all, aes(x = dband, y = spric, fill = location, shape = location), color = "grey69") +
  theme_bw() +
  scale_color_fish_d(option = "Bodianus_rufus") +
  scale_fill_fish_d(option = "Bodianus_rufus") +
  scale_shape_manual(values = c(21:24)) +
  theme(legend.position = "right",
        legend.title = element_blank())
spric.plot     

### Raos'sQ GAM
rao.gam <- gam(raoq ~ s(dband) + location, data = fd.res.all)
summary(rao.gam)


rao.newdat <- fd.res.all %>%
  data_grid(dband = seq(40,300,1),
            location = fd.res.all$location)
rao.pred.gam <- data.frame(predict(rao.gam, rao.newdat, se.fit = TRUE, type = "link")) %>%
  add_column(dband = rao.newdat$dband) %>%
  add_column(location = rao.newdat$location) %>%
  mutate(pred.rao = fit, pred.rao.se = se.fit) %>%
  mutate(lci = pred.rao-1.96*pred.rao.se,
         uci = pred.rao+1.96*pred.rao.se) %>%
  unite(gv, dband, location, sep = "_", remove = F)

rao.plot <- ggplot(rao.pred.gam, aes(x = dband, y = pred.rao, color = location)) +
  geom_line(aes(group = location, color = location)) +
  geom_ribbon(aes(ymin = lci, ymax = uci, group = location, fill = location), alpha = 0.1) +
  geom_point(data = fd.res.all, aes(x = dband, y = raoq, fill = location, shape = location), color = "grey69") +
  theme_bw() +
  scale_color_fish_d(option = "Bodianus_rufus") +
  scale_fill_fish_d(option = "Bodianus_rufus") +
  scale_shape_manual(values = c(21:24)) +
  theme(legend.position = "none",
        legend.title = element_blank())
rao.plot                   

div.depth.plots <- (spric.plot / rao.plot) + theme(legend.position = "right") +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect")
  
div.depth.plots
ggsave(div.depth.plots, file = "divdepthplot.png", width = 6, height = 10)
### FRIC GAM
# fric.gam <- gam(fric ~ s(dband) + location, data = fd.res.all)
# summary(fric.gam)
# 
# 
# fric.newdat <- fd.res.all %>%
#   data_grid(dband = seq(40,300,1),
#             location = fd.res.all$location)
# fric.pred.gam <- data.frame(predict(fric.gam, fric.newdat, se.fit = TRUE, type = "link")) %>%
#   add_column(dband = fric.newdat$dband) %>%
#   add_column(location = fric.newdat$location) %>%
#   mutate(pred.fric = fit, pred.fric.se = se.fit) %>%
#   mutate(lci = pred.fric-1.96*pred.fric.se,
#          uci = pred.fric+1.96*pred.fric.se)
# 
# fric.plot <- ggplot(fric.pred.gam, aes(x = dband, y = pred.fric, color = location)) +
#   geom_line(aes(group = location, color = location)) +
#   geom_ribbon(aes(ymin = lci, ymax = uci, group = location, fill = location), alpha = 0.1) +
#   geom_point(data = fd.res.all, aes(x = dband, y = fric, fill = location, shape = location), color = "grey69") +
#   theme_bw() +
#   scale_color_fish_d(option = "Bodianus_rufus") +
#   scale_fill_fish_d(option = "Bodianus_rufus") +
#   scale_shape_manual(values = c(21:24))
# fric.plot       

# s = ggplot(fd.res.all, aes(x = spric, y = raoq, color = dband)) +
#   geom_point() +
#   geom_smooth(method = "loess")
# 
# data.loess <- loess(dband ~ spric * raoq, data = fd.res.all)
# xgrid <-  seq(min(fd.res.all$spric), max(fd.res.all$spric), 0.1)
# ygrid <-  seq(min(fd.res.all$raoq), max(fd.res.all$raoq), 0.1)
# # Generate a dataframe with every possible combination of wt and hp
# data.fit <-  expand.grid(spric = xgrid, raoq = ygrid)
# # Feed the dataframe into the loess model and receive a matrix output with estimates of
# # acceleration for each combination of wt and hp
# pred.dat <-  predict(data.loess, newdata = data.fit)
# mtrx.melt <- reshape2::melt(pred.dat, id.vars = c("spric", "raoq"), measure.vars = "dband")
# names(mtrx.melt) <- c("spric", "raoq", "dband")
# mtrx.melt$spric <- as.numeric(str_sub(mtrx.melt$spric, str_locate(mtrx.melt$spric, "=")[1,1] + 1))
# mtrx.melt$raoq <- as.numeric(str_sub(mtrx.melt$raoq, str_locate(mtrx.melt$raoq, "=")[1,1] + 1))
# 
# mtrx.filt <- mtrx.melt %>%
#   filter(dband >= 0)
# 
# plot1 <- ggplot(mtrx.filt, aes(x = spric, y = raoq, z = dband)) +
#   stat_contour(geom = "polygon", aes(fill = ..level..)) +
#   geom_tile(aes(fill = dband), alpha = 0.3) +
#   stat_contour(bins = 10) +
#   theme_bw()
# plot1
# 
# mtrx.filt$equalSpace <- cut(mtrx.filt$dband, 10)
# # Sort the segments in ascending order
# breaks <- levels(unique(mtrx.filt$equalSpace))
# 
# 
# plot2 <- ggplot(mtrx.filt, aes(x = spric, y = raoq, z = dband)) +
#   geom_tile(aes(fill = equalSpace)) +
#   geom_contour(color = "white", alpha = 0.5) +
#   theme_bw() +
#   theme(legend.position = "top",
#         legend.title = element_blank()) +
#   scale_fill_fish_d(option = "Prionace_glauca", direction = -1) +
#   geom_jitter(data = fd.res.all, aes(x = spric, y = raoq, shape = location), width = 0.25, height = 0.25, size = 2.5, fill = "white", color = "black") +
#   geom_smooth(data = fd.res.all, aes(x = spric, y = raoq), method = "lm", formula = y ~ x + I(x^2), color = "black") +
#   scale_shape_manual(values = c(21:24)) +
#   scale_x_reverse()
# plot2

# 
# cur.cwm.size <- data.frame("dband" = as.numeric(rownames(cur.com.fd)),
#                        "location" = "Curacao",
#                        "size" = fd.cur$CWM$size_Size.cm)
# bon.cwm.size <- data.frame("dband" = as.numeric(rownames(bon.com.fd)),
#                            "location" = "Bonaire",
#                            "spric" = fd.bon$nbsp,
#                            "size" = fd.bon$CWM$size_Size.cm)
# roa.cwm.size <- data.frame("dband" = as.numeric(rownames(roa.com.fd)),
#                            "location" = "Roatan",
#                            "spric" = fd.roa$nbsp,
#                            "size" = fd.roa$CWM$size_Size.cm)
# sta.cwm.size <- data.frame("dband" = as.numeric(rownames(sta.com.fd)),
#                            "location" = "Statia",
#                            "spric" = fd.sta$nbsp,
#                            "size" = fd.sta$CWM$size_Size.cm)
# all.cwm.size <- bind_rows(cur.cwm.size, bon.cwm.size, roa.cwm.size, sta.cwm.size) %>%
#   filter(spric > 1) %>%
#   filter(dband <= 300)
# 
# size <- ggplot(all.cwm.size, aes(x = dband, y = size, group = location)) +
#   geom_point(aes(shape = location, fill = location), color = "black", alpha = 0.5) +
#   stat_smooth(aes(color = location, fill = location), alpha = 0.25, method = "lm", formula = y ~ x + I(x^2)) +
#   theme_bw() +
#   scale_shape_manual(values = c(21:24)) +
#   scale_fill_fish_d(option = "Bodianus_rufus") +
#   scale_color_fish_d(option = "Bodianus_rufus")
# size
# 
# 
# 
# size.m <- lme(size ~ dband + I(dband^2), random = ~1|location, data = all.cwm.size)
# summary(size.m)
# par(mfrow = c(2,2))
# plot(size.m)
# size.new <- data.frame(dband = as.numeric(rep(seq(min(all.cwm.size$dband), 
#                                                  max(all.cwm.size$dband), 
#                                                  length.out = 200)),4),
#                        location = rep(c("Curacao", "Bonaire", "Roatan", "Statia"), each = 200))
# 
# predictions <- predictSE.lme(size.m, data.frame(size.new),se.fit = TRUE, level = 0)
# 
# size.pred <- size.new %>%
#   add_column("pred.size" = predictions$fit) %>%
#   add_column("se.fit" = predictions$se.fit) %>%
#   mutate(lci = pred.size-1.96*se.fit,
#          uci = pred.size+1.96*se.fit)
#   
# size.plot <- ggplot(size.pred, aes(x = dband, y = pred.size)) +
#   geom_line(color = "black") +
#   geom_ribbon(aes(x = dband, ymin = lci, ymax = uci), alpha = 0.2) +
#   geom_point(data = all.cwm.size, aes(x = dband, y = size, shape = location, fill = location), alpha = 0.75, color = "grey23") +
#   theme_bw() +
#   ylab("Community-weighted mean size") +
#   xlab("Depth (m)") +
#   scale_shape_manual(values = c(21:24)) +
#   scale_fill_fish_d(option = "Bodianus_rufus") +
#   scale_color_fish_d(option = "Bodianus_rufus")
# size.plot
# 



### try pca plots for CWM
cur.cwms <- fd.cur$CWM %>%
  select_if(stringr::str_detect(names(.), "_1"))
cur.cwm.size <- cur.cwms %>%
  select(contains("size")) %>%
  mutate(across(where(is.character), as.numeric)) %>%
  mutate("dband" = as.numeric(rownames(cur.com.fd)),
         "location" = "Curacao") %>%
  relocate(dband, location) 
cur.cwm.trophic <- cur.cwms %>%
  select(contains("trophic")) %>%
  mutate(across(where(is.character), as.numeric)) %>%
  mutate("dband" = as.numeric(rownames(cur.com.fd)),
  "location" = "Curacao") %>%
  relocate(dband, location)
cur.cwm.diet <- cur.cwms %>%
  select(contains("diet")) %>%
  mutate(across(where(is.character), as.numeric)) %>%
  mutate("dband" = as.numeric(rownames(cur.com.fd)),
         "location" = "Curacao") %>%
  relocate(dband, location)
cur.cwm.position <- cur.cwms %>%
  select(contains("position")) %>%
  mutate(across(where(is.character), as.numeric)) %>%
  mutate("dband" = as.numeric(rownames(cur.com.fd)),
         "location" = "Curacao") %>%
  relocate(dband, location)
cur.cwm.repro <- cur.cwms %>%
  select(contains("repro")) %>%
  mutate(across(where(is.character), as.numeric)) %>%
  mutate("dband" = as.numeric(rownames(cur.com.fd)),
         "location" = "Curacao") %>%
  relocate(dband, location)

bon.cwms <- fd.bon$CWM %>%
  select_if(stringr::str_detect(names(.), "_1"))

bon.cwm.size <- bon.cwms %>%
  select(contains("size")) %>%
  mutate(across(where(is.character), as.numeric)) %>%
  mutate("dband" = as.numeric(rownames(bon.com.fd)),
         "location" = "Bonaire") %>%
  relocate(dband, location)
bon.cwm.trophic <- bon.cwms %>%
  select(contains("trophic")) %>%
  mutate(across(where(is.character), as.numeric)) %>%
  mutate("dband" = as.numeric(rownames(bon.com.fd)),
         "location" = "Bonaire") %>%
  relocate(dband, location)
bon.cwm.diet <- bon.cwms %>%
  select(contains("diet")) %>%
  mutate(across(where(is.character), as.numeric)) %>%
  mutate("dband" = as.numeric(rownames(bon.com.fd)),
         "location" = "Bonaire") %>%
  relocate(dband, location)
bon.cwm.position <- bon.cwms %>%
  select(contains("position")) %>%
  mutate(across(where(is.character), as.numeric)) %>%
  mutate("dband" = as.numeric(rownames(bon.com.fd)),
         "location" = "Bonaire") %>%
  relocate(dband, location)
bon.cwm.repro <- bon.cwms %>%
  select(contains("repro")) %>%
  mutate(across(where(is.character), as.numeric)) %>%
  mutate("dband" = as.numeric(rownames(bon.com.fd)),
         "location" = "Bonaire") %>%
  relocate(dband, location)

### roatan

roa.cwms <- fd.roa$CWM %>%
  select_if(stringr::str_detect(names(.), "_1"))
roa.cwm.size <- roa.cwms %>%
  select(contains("size")) %>%
  mutate(across(where(is.character), as.numeric)) %>%
  mutate("dband" = as.numeric(rownames(roa.com.fd)),
         "location" = "Roatan") %>%
  relocate(dband, location)
roa.cwm.trophic <- roa.cwms %>%
  select(contains("trophic")) %>%
  mutate(across(where(is.character), as.numeric)) %>%
  mutate("dband" = as.numeric(rownames(roa.com.fd)),
         "location" = "Roatan") %>%
  relocate(dband, location)
roa.cwm.diet <- roa.cwms %>%
  select(contains("diet")) %>%
  mutate(across(where(is.character), as.numeric)) %>%
  mutate("dband" = as.numeric(rownames(roa.com.fd)),
         "location" = "Roatan") %>%
  relocate(dband, location)
roa.cwm.position <- roa.cwms %>%
  select(contains("position")) %>%
  mutate(across(where(is.character), as.numeric)) %>%
  mutate("dband" = as.numeric(rownames(roa.com.fd)),
         "location" = "Roatan") %>%
  relocate(dband, location)
roa.cwm.repro <- roa.cwms %>%
  select(contains("repro")) %>%
  mutate(across(where(is.character), as.numeric)) %>%
  mutate("dband" = as.numeric(rownames(roa.com.fd)),
         "location" = "Roatan") %>%
  relocate(dband, location)

### Statia
sta.cwms <- fd.sta$CWM %>%
  select_if(stringr::str_detect(names(.), "_1"))

sta.cwm.size <- sta.cwms %>%
  select(contains("size")) %>%
  mutate(across(where(is.character), as.numeric)) %>%
  mutate("dband" = as.numeric(rownames(sta.com.fd)),
         "location" = "Statia") %>%
  relocate(dband, location)


sta.cwm.trophic <- sta.cwms %>%
  select(contains("trophic")) %>%
  mutate(across(where(is.character), as.numeric)) %>%
  mutate("dband" = as.numeric(rownames(sta.com.fd)),
         "location" = "Statia") %>%
  relocate(dband, location)

sta.cwm.diet <- sta.cwms %>%
  select(contains("diet")) %>%
  mutate(across(where(is.character), as.numeric)) %>%
  mutate("dband" = as.numeric(rownames(sta.com.fd)),
         "location" = "Statia") %>%
  relocate(dband, location)
sta.cwm.position <- sta.cwms %>%
  select(contains("position")) %>%
  mutate(across(where(is.character), as.numeric)) %>%
  mutate("dband" = as.numeric(rownames(sta.com.fd)),
         "location" = "Statia") %>%
  relocate(dband, location)
sta.cwm.repro <- sta.cwms %>%
  select(contains("repro")) %>%
  mutate(across(where(is.character), as.numeric)) %>%
  mutate("dband" = as.numeric(rownames(sta.com.fd)),
         "location" = "Statia") %>%
  relocate(dband, location)


cwm.trophic <- bind_rows(cur.cwm.trophic, bon.cwm.trophic, roa.cwm.trophic, sta.cwm.trophic) %>%
  select_if(negate(function(col) is.numeric(col) && sum(col) == 0)) %>%
  filter(dband > 40 & dband <= 300)
cwm.env <- cwm.trophic %>%
  select(dband, location) 

cwm.trophic.rda <- rda(cwm.trophic[-c(1:2)] ~ dband + location, data = cwm.env)
summary(cwm.trophic.rda)
anova(cwm.trophic.rda, by="term")
RsquareAdj(cwm.trophic.rda)

cwm.trophic.scores <- as.tibble(cwm.trophic.rda$CCA$wa) %>%
  bind_cols(cwm.trophic)
cwm.trophic.constraints <- data.frame(cwm.trophic.rda$CCA$biplot) %>%
  mutate(constraint = rownames(.))
cwm.trophic.vectors <- data.frame(cwm.trophic.rda$CCA$v) %>%
  mutate(vector = rownames(.)) %>%
  mutate(trophic = str_remove(vector, "trophic_"))

cwm.trophic.plot <- ggplot(cwm.trophic.scores, aes(x = RDA1, y = RDA2)) +
  geom_jitter(aes(fill = dband, shape = location), color = "grey69", size = 2, width = 0.01, height = 0.01) +
  geom_segment(data = cwm.trophic.constraints, 
               aes(x = 0, y = 0, 
                   xend = RDA1,
                   yend = RDA2)) +
  annotate('text', x = (cwm.trophic.constraints$RDA1), y = (cwm.trophic.constraints$RDA2),
           label = cwm.trophic.constraints$constraint,
           size = 2) +
  geom_point(data = cwm.trophic.vectors, aes(x = RDA1, y = RDA2), fill = "firebrick", color = "black", shape = 22, size = 2.5) +
  geom_text_repel(data = cwm.trophic.vectors, aes(x = RDA1, y = RDA2, label = trophic), size = 2, color = "firebrick") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_fish(option = "Prionace_glauca", direction = -1) +
  scale_shape_manual(values = c(21:24))
cwm.trophic.plot  


### diet
cwm.diet<- bind_rows(cur.cwm.diet, bon.cwm.diet, roa.cwm.diet, sta.cwm.diet) %>%
  select_if(negate(function(col) is.numeric(col) && sum(col) == 0)) %>%
  filter(dband > 40 & dband <= 300)
cwm.env <- cwm.diet %>%
  select(dband, location) 

cwm.diet.rda <- rda(cwm.diet[-c(1:2)] ~ dband + location, data = cwm.env)
summary(cwm.diet.rda)
anova(cwm.diet.rda, by="term")
RsquareAdj(cwm.diet.rda)

cwm.diet.scores <- as.tibble(cwm.diet.rda$CCA$wa) %>%
  bind_cols(cwm.diet)
cwm.diet.constraints <- data.frame(cwm.diet.rda$CCA$biplot) %>%
  mutate(constraint = rownames(.))
cwm.diet.vectors <- data.frame(cwm.diet.rda$CCA$v) %>%
  mutate(vector = rownames(.)) %>%
  mutate(diet = str_remove(vector, "diet_"))

cwm.diet.plot <- ggplot(cwm.diet.scores, aes(x = RDA1, y = RDA2)) +
  geom_jitter(aes(fill = dband, shape = location), color = "grey69", size = 2, width = 0.01, height = 0.01) +
  geom_segment(data = cwm.diet.constraints, 
               aes(x = 0, y = 0, 
                   xend = RDA1,
                   yend = RDA2)) +
  annotate('text', x = (cwm.diet.constraints$RDA1), y = (cwm.diet.constraints$RDA2),
           label = cwm.diet.constraints$constraint,
           size = 2) +
  geom_point(data = cwm.diet.vectors, aes(x = RDA1, y = RDA2), fill = "firebrick", color = "black", shape = 22, size = 2.5) +
  geom_text_repel(data = cwm.diet.vectors, aes(x = RDA1, y = RDA2, label = diet), size = 2, color = "firebrick") +
  theme_bw() +
  scale_fill_fish(option = "Prionace_glauca", direction = -1) +
  scale_shape_manual(values = c(21:24))
cwm.diet.plot  

#### position
cwm.position <- bind_rows(cur.cwm.position, bon.cwm.position, roa.cwm.position, sta.cwm.position) %>%
  select_if(negate(function(col) is.numeric(col) && sum(col) == 0)) %>%
  filter(dband > 40 & dband <= 300)
cwm.env <- cwm.position %>%
  select(dband, location) 

cwm.position.rda <- rda(cwm.position[-c(1:2)] ~ dband + location, data = cwm.env)
summary(cwm.position.rda)
anova(cwm.position.rda, by="term")
RsquareAdj(cwm.position.rda)

cwm.position.scores <- as.tibble(cwm.position.rda$CCA$wa) %>%
  bind_cols(cwm.position)
cwm.position.constraints <- data.frame(cwm.position.rda$CCA$biplot) %>%
  mutate(constraint = rownames(.))
cwm.position.vectors <- data.frame(cwm.position.rda$CCA$v) %>%
  mutate(vector = rownames(.)) %>%
  mutate(position = str_remove(vector, "position_"))

cwm.position.plot <- ggplot(cwm.position.scores, aes(x = RDA1, y = RDA2)) +
  geom_jitter(aes(fill = dband, shape = location), color = "grey69", size = 2, width = 0.01, height = 0.01) +
  geom_segment(data = cwm.position.constraints, 
               aes(x = 0, y = 0, 
                   xend = RDA1,
                   yend = RDA2)) +
  annotate('text', x = (cwm.position.constraints$RDA1), y = (cwm.position.constraints$RDA2),
           label = cwm.position.constraints$constraint,
           size = 2) +
  geom_point(data = cwm.position.vectors, aes(x = RDA1, y = RDA2), fill = "firebrick", color = "black", shape = 22, size = 2.5) +
  geom_text_repel(data = cwm.position.vectors, aes(x = RDA1, y = RDA2, label = position), size = 2, color = "firebrick") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_fish(option = "Prionace_glauca", direction = -1) +
  scale_shape_manual(values = c(21:24))
cwm.position.plot  



#### repro
cwm.repro <- bind_rows(cur.cwm.repro, bon.cwm.repro, roa.cwm.repro, sta.cwm.repro) %>%
  select_if(negate(function(col) is.numeric(col) && sum(col) == 0)) %>%
  filter(dband > 40 & dband <= 300)
cwm.env <- cwm.repro %>%
  select(dband, location) 

cwm.repro.rda <- rda(cwm.repro[-c(1:2)] ~ dband + location, data = cwm.env)
summary(cwm.repro.rda)
anova(cwm.repro.rda, by="term")
RsquareAdj(cwm.repro.rda)

cwm.repro.scores <- as.tibble(cwm.repro.rda$CCA$wa) %>%
  bind_cols(cwm.repro)
cwm.repro.constraints <- data.frame(cwm.repro.rda$CCA$biplot) %>%
  mutate(constraint = rownames(.))
cwm.repro.vectors <- data.frame(cwm.repro.rda$CCA$v) %>%
  mutate(vector = rownames(.)) %>%
  mutate(repro = str_remove(vector, "repro_"))

cwm.repro.plot <- ggplot(cwm.repro.scores, aes(x = RDA1, y = RDA2)) +
  geom_jitter(aes(fill = dband, shape = location), color = "grey69", size = 2, width = 0.01, height = 0.01) +
  geom_segment(data = cwm.repro.constraints, 
               aes(x = 0, y = 0, 
                   xend = RDA1,
                   yend = RDA2)) +
  annotate('text', x = (cwm.repro.constraints$RDA1), y = (cwm.repro.constraints$RDA2),
           label = cwm.repro.constraints$constraint,
           size = 2) +
  geom_point(data = cwm.repro.vectors, aes(x = RDA1, y = RDA2), fill = "firebrick", color = "black", shape = 22, size = 2.5) +
  geom_text_repel(data = cwm.repro.vectors, aes(x = RDA1, y = RDA2, label = repro), size = 2, color = "firebrick") +
  theme_bw() +
  scale_fill_fish(option = "Prionace_glauca", direction = -1) +
  scale_shape_manual(values = c(21:24))
cwm.repro.plot  

cwm.plots <- (cwm.trophic.plot + cwm.diet.plot + cwm.position.plot + cwm.repro.plot & theme(legend.position = "right")) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect")
cwm.plots
ggsave(cwm.plots, file = "cwm.plots.png", width = 12, height = 10)



#### CWM beta models
cwm.size <- bind_rows(cur.cwm.size, bon.cwm.size, roa.cwm.size, sta.cwm.size) %>%
  select_if(negate(function(col) is.numeric(col) && sum(col) == 0)) 

cwm.size.beta <- cwm.size %>%
  gather(3:7, key = "sizeclass", value = "cwm") %>%
  mutate(sizeclass = as.factor(str_remove(sizeclass, "_1")))

cwm.size.gam <- gam(cwm ~ s(dband, by = sizeclass), family = betar, data = cwm.size.beta, method = "REML")
summary(cwm.size.gam)

pred.beta.new <- cwm.size.beta %>%
  data_grid(dband = seq(40,300,1),
            sizeclass = cur.beta$sizeclass)

size.pred.beta <- data.frame(predict(cwm.size.gam, pred.beta.new, se.fit = TRUE)) %>%
  add_column(dband = pred.beta.new$dband) %>%
  add_column(sizeclass = pred.beta.new$sizeclass) %>%
  mutate(lci = fit-1.96*se.fit,
         uci = fit+1.96*se.fit) %>%
  mutate(pred.cwm = boot::inv.logit(fit), pred.lci = boot::inv.logit(lci), pred.uci = boot::inv.logit(uci)) %>%
  mutate(sizeclass = str_remove(sizeclass, "size_")) %>%
  mutate(sizeclass.reord = recode(sizeclass, l = 4,
                                  m = 3,
                                  s = 2,
                                  xl = 5,
                                  xs = 1))
  
size.plot <- ggplot(size.pred.beta, aes(x = dband, y = pred.cwm)) +
  geom_line(aes(group = forcats::fct_reorder(sizeclass, sizeclass.reord), color = forcats::fct_reorder(sizeclass, sizeclass.reord)), lty = 1, lwd = 0.5) +
  geom_ribbon(aes(ymin = pred.lci, ymax = pred.uci, group = forcats::fct_reorder(sizeclass, sizeclass.reord), fill = forcats::fct_reorder(sizeclass, sizeclass.reord)), alpha = 0.1) +
  theme_bw() +
  scale_color_fish_d(option = "Centropyge_loricula") +
  scale_fill_fish_d(option = "Centropyge_loricula") +
  theme(legend.title = element_blank())
size.plot                   


###cwm trophic
#### CWM beta models
cwm.trophic <- bind_rows(cur.cwm.trophic, bon.cwm.trophic, roa.cwm.trophic, sta.cwm.trophic) %>%
  select_if(negate(function(col) is.numeric(col) && sum(col) == 0)) %>%
  filter(dband > 40 & dband <= 300)

cwm.trophic.beta <- cwm.trophic %>%
  gather(3:9, key = "trophic", value = "cwm") %>%
  mutate(trophic = as.factor(str_remove(trophic, "trophic_")))%>%
  mutate(trophic = as.factor(str_remove(trophic, "_1")))


cwm.trophic.gam <- gam(cwm ~ s(dband, by = trophic), bs = "cs", family = betar, data = cwm.trophic.beta, method = "REML")
summary(cwm.trophic.gam)

trophic.pred.beta.new <- cwm.trophic.beta %>%
  data_grid(dband = seq(40,300,1),
            trophic = cwm.trophic.beta$trophic)

trophic.pred.beta <- data.frame(predict(cwm.trophic.gam, trophic.pred.beta.new, se.fit = TRUE)) %>%
  add_column(dband = trophic.pred.beta.new$dband) %>%
  add_column(trophic = trophic.pred.beta.new$trophic) %>%
  mutate(lci = fit-1.96*se.fit,
         uci = fit+1.96*se.fit) %>%
  mutate(pred.cwm = boot::inv.logit(fit), pred.lci = boot::inv.logit(lci), pred.uci = boot::inv.logit(uci))

trophic.plot <- ggplot(trophic.pred.beta, aes(x = dband, y = pred.cwm)) +
  geom_line(aes(group = trophic, color = trophic), lty = 1, lwd = 0.5) +
  geom_ribbon(aes(ymin = pred.lci, ymax = pred.uci, group = trophic, fill = trophic), alpha = 0.1) +
  theme_bw() +
  scale_color_fish_d(option = "Centropyge_loricula") +
  scale_fill_fish_d(option = "Centropyge_loricula") +
  theme(legend.title = element_blank())
trophic.plot    


#### CWM diet models
cwm.diet <- bind_rows(cur.cwm.diet, bon.cwm.diet, roa.cwm.diet, sta.cwm.diet) %>%
  select_if(negate(function(col) is.numeric(col) && sum(col) == 0)) %>%
  filter(dband > 40 & dband <= 300)

cwm.diet.beta <- cwm.diet %>%
  gather(3:11, key = "diet", value = "cwm") %>%
  mutate(diet = as.factor(str_remove(diet, "diet_")))%>%
  mutate(diet = as.factor(str_remove(diet, "_1")))


cwm.diet.gam <- gam(cwm ~ s(dband, by = diet), bs = "cs", family = betar, data = cwm.diet.beta, method = "REML")
summary(cwm.diet.gam)

diet.pred.beta.new <- cwm.diet.beta %>%
  data_grid(dband = seq(40,300,1),
            diet = cwm.diet.beta$diet)

diet.pred.beta <- data.frame(predict(cwm.diet.gam, diet.pred.beta.new, se.fit = TRUE)) %>%
  add_column(dband = diet.pred.beta.new$dband) %>%
  add_column(diet = diet.pred.beta.new$diet) %>%
  mutate(lci = fit-1.96*se.fit,
         uci = fit+1.96*se.fit) %>%
  mutate(pred.cwm = boot::inv.logit(fit), pred.lci = boot::inv.logit(lci), pred.uci = boot::inv.logit(uci))

diet.plot <- ggplot(diet.pred.beta, aes(x = dband, y = pred.cwm)) +
  geom_line(aes(group = diet, color = diet), lty = 1, lwd = 0.5) +
  geom_ribbon(aes(ymin = pred.lci, ymax = pred.uci, group = diet, fill = diet), alpha = 0.1) +
  theme_bw() +
  scale_color_fish_d(option = "Centropyge_loricula") +
  scale_fill_fish_d(option = "Centropyge_loricula") +
  theme(legend.title = element_blank())
diet.plot  

#### CWM position models
cwm.position <- bind_rows(cur.cwm.position, bon.cwm.position, roa.cwm.position, sta.cwm.position) %>%
  select_if(negate(function(col) is.numeric(col) && sum(col) == 0)) %>%
  filter(dband > 40 & dband <= 300)

cwm.position.beta <- cwm.position %>%
  gather(3:7, key = "position", value = "cwm") %>%
  mutate(position = as.factor(str_remove(position, "position_")))%>%
  mutate(position = as.factor(str_remove(position, "_1")))


cwm.position.gam <- gam(cwm ~ s(dband, by = position), bs = "cs", family = betar, data = cwm.position.beta, method = "REML")
summary(cwm.position.gam)

position.pred.beta.new <- cwm.position.beta %>%
  data_grid(dband = seq(40,300,1),
            position = cwm.position.beta$position)

position.pred.beta <- data.frame(predict(cwm.position.gam, position.pred.beta.new, se.fit = TRUE)) %>%
  add_column(dband = position.pred.beta.new$dband) %>%
  add_column(position = position.pred.beta.new$position) %>%
  mutate(lci = fit-1.96*se.fit,
         uci = fit+1.96*se.fit) %>%
  mutate(pred.cwm = boot::inv.logit(fit), pred.lci = boot::inv.logit(lci), pred.uci = boot::inv.logit(uci))

position.plot <- ggplot(position.pred.beta, aes(x = dband, y = pred.cwm)) +
  geom_line(aes(group = position, color = position), lty = 1, lwd = 0.5) +
  geom_ribbon(aes(ymin = pred.lci, ymax = pred.uci, group = position, fill = position), alpha = 0.1) +
  theme_bw() +
  scale_color_fish_d(option = "Centropyge_loricula") +
  scale_fill_fish_d(option = "Centropyge_loricula") +
  theme(legend.title = element_blank())
position.plot  

###
#### CWM diet models
cwm.repro <- bind_rows(cur.cwm.repro, bon.cwm.repro, roa.cwm.repro, sta.cwm.repro) %>%
  select_if(negate(function(col) is.numeric(col) && sum(col) == 0)) %>%
  filter(dband > 40 & dband <= 300)

cwm.repro.beta <- cwm.repro %>%
  gather(3:7, key = "repro", value = "cwm") %>%
  mutate(repro = as.factor(str_remove(repro, "repro_")))%>%
  mutate(repro = as.factor(str_remove(repro, "_1")))


cwm.repro.gam <- gam(cwm ~ s(dband, by = repro), bs = "cs", family = betar, data = cwm.repro.beta, method = "REML")
summary(cwm.repro.gam)

repro.pred.beta.new <- cwm.repro.beta %>%
  data_grid(dband = seq(40,300,1),
            repro = cwm.repro.beta$repro)

repro.pred.beta <- data.frame(predict(cwm.repro.gam, repro.pred.beta.new, se.fit = TRUE)) %>%
  add_column(dband = repro.pred.beta.new$dband) %>%
  add_column(repro = repro.pred.beta.new$repro) %>%
  mutate(lci = fit-1.96*se.fit,
         uci = fit+1.96*se.fit) %>%
  mutate(pred.cwm = boot::inv.logit(fit), pred.lci = boot::inv.logit(lci), pred.uci = boot::inv.logit(uci))

repro.plot <- ggplot(repro.pred.beta, aes(x = dband, y = pred.cwm)) +
  geom_line(aes(group = repro, color = repro), lty = 1, lwd = 0.5) +
  geom_ribbon(aes(ymin = pred.lci, ymax = pred.uci, group = repro, fill = repro), alpha = 0.1) +
  theme_bw() +
  scale_color_fish_d(option = "Centropyge_loricula") +
  scale_fill_fish_d(option = "Centropyge_loricula") +
  theme(legend.title = element_blank())
repro.plot  

cwm.plots.beta <- size.plot / trophic.plot / diet.plot / position.plot / repro.plot +
  plot_annotation(tag_levels = 'A')
ggsave(cwm.plots.beta, file = "cwm.plots.beta.png", width = 8, height = 15)
