library(tidyverse) #
library(modelr) #part of the tidyverse ecosystem and provides tools for creating and working with data grids, which are often used in modeling and data analysis tasks.
library(fishualize) # provides color scales based on fish color
library(ggrepel)# extension of the popular ggplot2 package and provides additional functionality for controlling the placement and positioning of text labels in a ggplot2 plot.
library(tinter) # generates palettes / tints for graphs
library(mgcv) # mixed generalized additive models 
library(GGally) #functions and features that complement the functionality of the popular ggplot2 package.
library(hillR) #to calculate richness, evenness and dominance 
library(betapart) #to calculate beta diversity 


#### 0.1 Data prep - fish observations ####
# read in raw data files
bonaire <- read.csv(file = "raw.data/raw fish counts/Bonaire_Raw_2020.csv") 
curacao <- read.csv(file = "raw.data/raw fish counts/Curacao_Raw_2020.csv")
roatan <- read.csv(file = "raw.data/raw fish counts/Roatan_Raw_2020.csv") 
statia <- read.csv(file = "raw.data/raw fish counts/Statia_raw_2020.csv") 

#### 0.1.1. Correction factors (to normalize abundance by time spent) ####
bon.corr <- read.csv("raw.data/CorrectionFactors_SamplingEfforts/BonCoeff_2020.csv") %>%
  separate(depth, into = c("dband", "m")) %>%
  mutate(dband = as.numeric(dband)) %>%
  mutate(location = "Bonaire")
cur.corr <- read.csv("raw.data/CorrectionFactors_SamplingEfforts/CurCoeff_2020.csv") %>%
  separate(depth, into = c("dband", "m")) %>%
  mutate(dband = as.numeric(dband)) %>%
  mutate(location = "Curacao")
roa.corr <- read.csv("raw.data/CorrectionFactors_SamplingEfforts/RoaCoeff_2020.csv") %>%
  separate(depth, into = c("dband", "m")) %>%
  mutate(dband = as.numeric(dband)) %>%
  mutate(location = "Roatan")
sta.corr <- read.csv("raw.data/CorrectionFactors_SamplingEfforts/StaCoeff_2020.csv") %>%
  separate(X, into = c("dband", "m")) %>%
  rename(COEF = "coeff") %>%
  mutate(dband = as.numeric(dband)) %>%
  mutate(location = "St. Eustatius")

coefs <- bind_rows(bon.corr, cur.corr, roa.corr, sta.corr)

#write.csv(coefs,"raw.data/CorrectionFactors_SamplingEfforts/file.correction.coeffs.grouped.csv")


#### 0.1.2. Check species name and assign depth bands ####
carib.deep <- bind_rows(bonaire, curacao, roatan, statia)%>%
  rename(species = "name")

# look for false names -- need to replace with FB compatible names
# choosing close congeners

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
  #distinct() %>% # selects distinct rows from data frame 
  mutate(species.fb = case_when(is.na(species.fb) ~ species,
                                TRUE ~ species.fb)) # attributes the original sp. name when no correction is needed
#write.csv(carib.deep.fishbase,"2.Ch3.R.analysis/2.file.fish.grouped.raw.csv")

coefs <- read.csv("raw.data/CorrectionFactors_SamplingEfforts/file.correction.coeffs.grouped.csv")
carib.deep.fishbase <- carib.deep.fishbase %>%
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
  group_by(species, dband, location) %>%
  summarize(abundance = n()) %>% # sum of the number of observation per depth, location, species
  left_join(coefs) %>%
  mutate(abu.corr = abundance*COEF) ## abundance corrected by sampling effort coeff

#write.csv(carib.deep.fishbase,"2.Ch3.R.analysis/2.file.fish.carib.grouped.csv")

#### 1.1 RAW abundance ####
#### 1.1.0. data prep ####
carib.deep.traits <- read.csv(file = "1.Ch3.R.analysis/2.carib.deep.traits.csv")

## total abundance 
total_abundance <- read.csv(file = "1.Ch3.R.analysis/2.carib.deep.traits.csv")%>%
  filter(dband>30 & dband<310) %>%
  summarize(tot.abu=sum(abundance),tot.species=length(unique(species)))
total_abundance

## list of species observed + their max depth and total length
carib.deep.traits.species <- carib.deep.traits %>%
  select(species, maxTL, maxDepth) %>%
  distinct()

# write.csv(carib.deep.traits.species, file = "2.carib.deep.traits.species.csv")

## total abundance and biomass (across species) per location per depth
carib.fish.sum <- carib.deep.traits %>%
  group_by(dband, location, species) %>%
  summarize(tot.abu = sum(abu.corr), u.abu = sum(abundance), tot.biom = sum(biomass)) %>%
  ungroup() %>%
  group_by(dband, location) %>%
  summarize(abundance = sum(tot.abu), u.abundance = sum(u.abu), biomass = sum(tot.biom), spric = n()) %>%
  bind_rows(data_frame(dband=c(240,250),location=c("Bonaire","St. Eustatius"),abundance=0,u.abundance=0,
                       biomass=0,spric=0)) 
  


#### 1.1.1. Curacao ####
cur.abu <- carib.fish.sum %>%
  filter(location == "Curacao")

cur.abu.gam <- gam(log10(abundance) ~ s(dband), data = cur.abu)
summary(cur.abu.gam) # sign correlation (<0.001) between abundance and depth bands
par(mfrow = c(2,2))
gam.check(cur.abu.gam)

## modelling continuous variation of abundance across depth based on results
# sub dividing the depth range into 500 values to plot a continuous prediction line
cur.abu.new <- data.frame(dband = as.numeric(seq(min(cur.abu$dband), 
                                                 max(cur.abu$dband), 
                                                 length.out = 500)))

## calculating predicted values of abundance for the 500 depth values based on the gam model
## includes SE associated with predictions
cur.abu.pred <- as.data.frame(predict(cur.abu.gam, data.frame(cur.abu.new), se.fit = T)) %>%
  rename(abu_predicted = 1) %>%
  mutate(abu_predicted_raw = exp(abu_predicted)) %>% # log values back to raw
  mutate(lci = abu_predicted-1.96*se.fit,
         uci = abu_predicted+1.96*se.fit) %>%  # upper and lower 95CI
  bind_cols(cur.abu.new) # adds depth values

## abundance of fish in Curacao
cur.abu.plot <- ggplot(cur.abu.pred, aes(x = dband, y = abu_predicted)) +
  theme_classic()+
  geom_line(lty = 1, color = carib.pal[1]) + # predicted values 
  geom_ribbon(aes(x = dband, ymin = lci, ymax = uci), alpha = 0.2, fill = carib.pal[1]) + # 95 CI
  geom_point(data = cur.abu, aes(x = dband, y = log10(abundance)), shape = 21,
             alpha = 0.75, color = "grey23", fill = carib.pal[1]) + # original data
  geom_vline(xintercept = 100, lty = 2, color = "grey69") +
  geom_vline(xintercept = 150, lty = 2, color = "grey69") +
  geom_vline(xintercept = 230, lty = 2, color = "grey69") +
  geom_vline(xintercept = 290, lty = 2, color = "grey69") +
  ylab("log10 abundance") +
  xlab("depth (m)") +
  scale_x_continuous(limits = c(40,300), breaks = seq(40,300,20))+
  ylim(0.3,4)
cur.abu.plot


#### 1.1.2. Bonaire ####
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

ggplot(bon.abu.pred, aes(x = dband, y = abu_predicted)) +
  geom_line(lty = 1, color = carib.pal[2]) +
  geom_ribbon(aes(x = dband, ymin = lci, ymax = uci), alpha = 0.2, fill = carib.pal[2]) +
  geom_point(data = bon.abu, aes(x = dband, y = log10(abundance)), shape = 22, alpha = 0.75, color = "grey23", fill = carib.pal[2]) +
  geom_vline(xintercept = 100, lty = 2, color = "grey69") +
  geom_vline(xintercept = 160, lty = 2, color = "grey69") +
  geom_vline(xintercept = 230, lty = 2, color = "grey69") +
  geom_vline(xintercept = 280, lty = 2, color = "grey69") +
  theme_classic() +
  ylab("log10 abundance") +
  xlab("depth (m)") +
  scale_x_continuous(limits = c(40,300), breaks = seq(40,300,20))+
  ylim(-0.1,4)



#### 1.1.3. Roatan ####
roa.abu <- carib.fish.sum %>%
  filter(location == "Roatan")
  # filter(dband < 300) 

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
  geom_line(lty = 1, color = carib.pal[4]) +
  geom_ribbon(aes(x = dband, ymin = lci, ymax = uci), alpha = 0.2, fill = carib.pal[4]) +
  geom_point(data = roa.abu, aes(x = dband, y = log10(abundance)), shape = 23, 
             alpha = 0.75, color = "grey23", fill = carib.pal[4]) +
  geom_vline(xintercept = 90, lty = 2, color = "grey69") +
  geom_vline(xintercept = 140, lty = 2, color = "grey69") +
  geom_vline(xintercept = 200, lty = 2, color = "grey69") +
  geom_vline(xintercept = 300, lty = 2, color = "grey69") +
  theme_classic() +
  ylab("log10 abundance") +
  xlab("depth (m)")+
  ylim(0.3,4)+
  xlim(0,500)
# scale_x_continuous(limits = c(40,max(roa.abu$dband)), breaks = seq(40,max(roa.abu$dband),10))
roa.abu.plot


#### 1.1.4. Statia ####
sta.abu <- carib.fish.sum %>%
  filter(location == "St. Eustatius") 

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
  theme_classic() +
  ylab("log10 abundance") +
  xlab("depth (m)") +
  scale_x_continuous(limits = c(40,max(sta.abu$dband)), breaks = seq(40,max(sta.abu$dband),20))+
  ylim(0.3,4)
sta.abu.plot


abundance.plots <- (cur.abu.plot + bon.abu.plot + roa.abu.plot + sta.abu.plot) + plot_annotation(tag_levels = "A")
#ggsave(abundance.plots, file = "abundance.plots.lin.png", width = 14, height = 8)

#### 1.1.5. all locations in same plot ####
FD.values <- read.csv("file_FD.all.islands.csv") ## used to build prediction table
carib.fish.sum$location <- as.character(carib.fish.sum$location)
carib.fish.sum$location[carib.fish.sum$location=="St. Eustatius"] <- "Statia"

abu.gam = gam(abundance ~ s(dband) + location, family = poisson, data = carib.fish.sum)
summary(spric.gam)

# empty table with 1 line per (meter depth x site) to predict abundance
abu.newdat <- FD.values %>%
  # data_grid expands a table using all variables provided
  data_grid(dband = seq(40,300,1), 
            location = FD.values$location)%>%
  mutate(location = factor(location, levels = c("Curacao", "Bonaire","Statia", "Roatan")))

abu.pred.gam <- data.frame(predict(abu.gam, abu.newdat, se.fit = TRUE, type = "link")) %>%
  add_column(dband = abu.newdat$dband) %>%
  add_column(location = abu.newdat$location) %>%
  mutate(pred.abu = exp(fit), pred.abu.se = exp(se.fit)) %>%
  mutate(lci = pred.abu-1.96*pred.abu.se,
         uci = pred.abu+1.96*pred.abu.se)%>%
  mutate(location = factor(location, levels = c("Curacao", "Bonaire","Statia", "Roatan")))

carib.fish.sum <- carib.fish.sum %>%
  mutate(location = factor(location, levels = c("Curacao", "Bonaire","Statia", "Roatan")))


ggplot(abu.pred.gam, aes(x = dband, y = log10(pred.abu), color = location)) +
  geom_line(aes(group = location, color = location)) +
  geom_ribbon(aes(ymin = log10(lci), ymax = log10(uci), group = location, fill = location), alpha = 0.1) +
  geom_point(data = carib.fish.sum, aes(x = dband, y = log10(abundance), fill = location, shape = location), color = "grey69") +
  theme_bw() +
  scale_color_fish_d(option = "Bodianus_rufus") +
  scale_fill_fish_d(option = "Bodianus_rufus") +
  scale_shape_manual(values = c(21:24)) +
  theme(legend.position = "right",
        legend.title = element_blank())+
  ylab("abundance (log10)")+
  xlab("depth")+
  xlim(0,310)


#### 1.2. NORM abundance ####
#### 1.2.0. data prep ####
carib.deep.traits <- read.csv(file = "carib.deep.traits.csv")

## list of species observed + their max depth and total length
carib.deep.traits.species <- carib.deep.traits %>%
  select(species, maxTL, maxDepth) %>%
  distinct()

# write.csv(carib.deep.traits.species, file = "2.carib.deep.traits.species.csv")

## total abundance and biomass (across species) per location per depth
carib.fish.sum <- carib.deep.traits %>%
  group_by(dband, location, species) %>%
  summarize(tot.abu = sum(abu.corr), u.abu = sum(abundance), tot.biom = sum(biomass)) %>%
  ungroup() %>%
  group_by(dband, location) %>%
  summarize(abundance = sum(tot.abu), u.abundance = sum(u.abu), biomass = sum(tot.biom), spric = n()) %>%
  bind_rows(data_frame(dband=c(240,250),location=c("Bonaire","St. Eustatius"),abundance=0,u.abundance=0,
                       biomass=0,spric=0)) 

carib.fish.sum$location <- as.character(carib.fish.sum$location)
carib.fish.sum$location[carib.fish.sum$location=="St. Eustatius"] <- "Statia" 
carib.fish.sum <- carib.fish.sum %>%
  mutate(location = factor(location, levels = c("Curacao", "Bonaire","Statia", "Roatan")))

#### 1.2.1 facet wrap ####
carib.fish.sum <- carib.fish.sum %>%
  filter(dband>30 & dband<310)

ggplot(carib.fish.sum, aes(x = dband, y = abu.norm, color = location, fill = location)) +
  geom_point(aes(shape = location), alpha = 0.75) +
  geom_smooth(method = "gam") +
  theme_bw() +
  theme(legend.position = "none") +  # Remove legends
  ylab("abundance") +
  xlab("depth (m)")+
  scale_fill_fish_d(option = "Bodianus_rufus") +
  scale_color_fish_d(option = "Bodianus_rufus") +
  facet_wrap(.~location)



#### 2.0. Data prep (run before) ####
carib.fish.sum <- read.csv(file = "2.Ch3.R.analysis/2.file.fish.carib.grouped.csv") %>%
  #filter(dband>30 & dband <310) %>%
  group_by(dband, location, species) %>%
  summarize(abu.corr = sum(abu.corr), abu.raw = sum(abundance)) %>%
  ungroup() %>%
  group_by(dband, location) %>%
  summarize(abu.corr = sum(abu.corr), abu.raw = sum(abu.raw), spric = n()) %>%
  mutate(location = factor(location, levels = c("Curacao", "Bonaire","St. Eustatius", "Roatan")))

## normalizing from 0 to 1 using the maximal abundance value at each location
for (k in c("Curacao", "Bonaire","St. Eustatius", "Roatan")){
  max.abu <- max(carib.fish.sum[carib.fish.sum$location==k & carib.fish.sum$dband>30,"abu.corr"])
  print(max.abu)
  carib.fish.sum$abu.norm[carib.fish.sum$location==k] <- carib.fish.sum$abu.corr[carib.fish.sum$location==k]/max.abu
}
carib.fish.sum$abu.norm <- unlist(carib.fish.sum$abu.norm)*100

#### 2.1 abundance ####
## raw abundance
ggplot(carib.fish.sum, aes(x = dband, y = abu.raw , color = location, fill = location)) +
  geom_point(aes(shape = location), alpha = 0.75) +
  geom_smooth(method = "gam") +
  theme_bw() +
  theme(legend.position = "none") +  # Remove legends
  ylab("abundance") +
  xlab("depth (m)")+
  scale_fill_fish_d(option = "Bodianus_rufus") +
  scale_color_fish_d(option = "Bodianus_rufus") +
  facet_wrap(.~location)

## corrected abundance - facet wrap across locations
ggplot(carib.fish.sum, aes(x = dband, y = abu.corr, color = location, fill = location)) +
  geom_point(aes(shape = location), alpha = 0.75) +
  geom_smooth(method = "gam") +
  ylab("abundance") +
  xlab("depth (m)")+
  theme_bw()+
  theme(legend.position = "none")+  # Remove legends
  scale_fill_fish_d(option = "Bodianus_rufus") +
  scale_color_fish_d(option = "Bodianus_rufus") +
  scale_x_continuous(limits=c(40,300),breaks=c(40,100,200,300))+
  #ylim(-50000,350000)+
  facet_wrap(.~location,scales = "free") ## to have all scales adapted to each graph, use scales = "free"

#### FIGURE 2 - LEFT  corrected and normalized abundance - all on same plot ####
abu.gam = gam(abu.norm ~ s(dband) + location, family = poisson, data = carib.fish.sum)
summary(abu.gam)

# empty table with 1 line per (meter depth x site) to predict abundance
abu.newdat <- read.csv("2.Ch3.R.analysis/2.file.fish.carib.grouped.csv") %>%
  # data_grid expands a table using all variables provided
  data_grid(dband = seq(10,480,1), 
            location = c("Curacao", "Bonaire","St. Eustatius", "Roatan"))%>%
  mutate(location = factor(location, levels = c("Curacao", "Bonaire","St. Eustatius", "Roatan")))%>%
  filter((dband<300 & dband>40 & location!="Roatan") | location=="Roatan") # only keeps Roatan for 1-40 and 300-480m

abu.pred.gam <- data.frame(predict(abu.gam, abu.newdat, se.fit = TRUE, type = "link")) %>%
  add_column(dband = abu.newdat$dband) %>%
  add_column(location = abu.newdat$location) %>%
  mutate(pred.abu = exp(fit), pred.abu.se = exp(se.fit)) %>%
  mutate(lci = pred.abu-1.96*pred.abu.se,
         uci = pred.abu+1.96*pred.abu.se)%>%
  mutate(location = factor(location, levels = c("Curacao", "Bonaire","St. Eustatius", "Roatan")))

new.pal=c("#ffa600","#ef5675","#21918c","#003f5c")
#plot
ggplot(abu.pred.gam, aes(x = dband, y = pred.abu, color = location)) +
  geom_point(data = carib.fish.sum, aes(x = dband, y = abu.norm, fill = location,
                                        shape = location), color = "white") +
  geom_line(aes(group = location, color = location)) +
  geom_ribbon(aes(ymin = lci, ymax = uci, group = location, fill = location), 
              alpha = 0.2,color=NA) +
   theme_classic() +
  scale_color_manual(values=new.pal)+
  scale_fill_manual(values=new.pal)+
  #scale_color_fish_d(option = "Bodianus_rufus") +
  #scale_fill_fish_d(option = "Bodianus_rufus") +
  scale_shape_manual(values = c(21:24)) +
  #theme(legend.position = "none")+
  ylab("relative abundance")+
  xlab("depth")+
  xlim(10,480)+
  ylim(-7,100)

#### 2.3. species richness ####
spric.plot <- ggplot(carib.fish.sum, aes(x = dband, y = spric, color = location, fill = location)) +
  geom_point(aes(shape = location), alpha = 0.75) +
  geom_smooth(method = "gam") +
  theme_bw() +
  theme(legend.position = "none") +  # Remove legends
  ylab("species richness") +
  xlab("depth (m)")+
  scale_fill_fish_d(option = "Bodianus_rufus") +
  scale_color_fish_d(option = "Bodianus_rufus") +
  facet_wrap(.~location)
spric.plot
#ggsave(spric.plot, file = "Figures/spric.png", width = 8, height = 6)


#### FIG 2 RIGHT: species richness - all in same plot ####
spric.gam = gam(spric ~ s(dband) + location, family = poisson, data = carib.fish.sum)
summary(spric.gam)

spric.newdat <- read.csv("2.Ch3.R.analysis/2.file.fish.carib.grouped.csv") %>%
  # data_grid expands a table using all variables provided
  data_grid(dband = seq(10,480,1), 
            location = c("Curacao", "Bonaire","St. Eustatius", "Roatan"))%>%
  mutate(location = factor(location, levels = c("Curacao", "Bonaire","St. Eustatius", "Roatan")))%>%
  filter((dband<300 & dband>40 & location!="Roatan") | location=="Roatan") 

spric.pred.gam <- data.frame(predict(spric.gam, spric.newdat, se.fit = TRUE, type = "link")) %>%
  add_column(dband = spric.newdat$dband) %>%
  add_column(location = spric.newdat$location) %>%
  mutate(pred.spric = exp(fit), pred.spric.se = exp(se.fit)) %>%
  mutate(lci = pred.spric-1.96*pred.spric.se,
         uci = pred.spric+1.96*pred.spric.se)%>%
  mutate(location = factor(location, levels = c("Curacao", "Bonaire","St. Eustatius", "Roatan")))

new.pal=c("#ffa600","#ef5675","#21918c","#003f5c")

ggplot(spric.pred.gam, aes(x = dband, y = pred.spric, color = location)) +
  geom_point(data = carib.fish.sum, aes(x = dband, y = spric, fill = location, 
                                        shape = location),color = "white") +
  geom_ribbon(aes(ymin = lci, ymax = uci, group = location, fill = location), 
              alpha = 0.2,colour = NA) +
  geom_line(aes(group = location, color = location)) +
  theme_classic() +
  scale_color_manual(values=new.pal)+
  scale_fill_manual(values=new.pal)+
  #scale_color_fish_d(option = "Trimma_lantana") +
  #scale_fill_fish_d(option = "Trimma_lantana") +
  scale_shape_manual(values = c(21:24)) +
  #theme(legend.position = "none")+
  ylab("species richness")+
  xlab("depth")+
  xlim(10,480)

#### 3. Taxonomic diversity metrics using Hill numbers ####
fish.data <- read.csv(file = "2.Ch3.R.analysis/2.file.fish.carib.grouped.csv")
#### 3.1. Basic plotting ####
## data prep
fish.curacao <- fish.data %>%
  filter(location == "Curacao") %>%
  select(species,dband,abu.corr) 

hill.curacao <- fish.curacao%>%
  spread(species, abu.corr, fill =  0) 

### taxonomic diversity 
taxa.hill0 <- hill_taxa(hill.curacao,q=0)
taxa.hill1 <- hill_taxa(hill.curacao,q=1)
taxa.hill2 <- hill_taxa(hill.curacao,q=2)

div.curacao <- data.frame(depth=seq(40,300,10), richness=taxa.hill0, evenness=taxa.hill1,
                          dominance=taxa.hill2)

### total species richness
ncol(hill.curacao)

# richness raw
ggplot(div.curacao,aes(x=depth,y=richness))+
  geom_point()
# richness log10
ggplot(div.curacao,aes(x=depth,y=log10(richness)))+
  geom_point()

# evenness - Shannon
ggplot(div.curacao,aes(x=depth,y=evenness))+
  geom_point()

# dominance- Simpson
ggplot(div.curacao,aes(x=depth,y=dominance))+
  geom_point()

#### 3.2. Plots with prediction and CI ####
#### 3.2.1 Richness ####
## richness
cur.div.gam <- gam(richness ~ s(depth), data = div.curacao) ## here taking the "smooth" of dband.
summary(cur.div.gam) # sign correlation (<0.001) between richness and depth
par(mfrow = c(2,2))
gam.check(cur.div.gam)

# sub dividing the depth range into 500 values to plot a continuous prediction line
cur.div.new <- data.frame(depth = as.numeric(seq(40,300,length.out = 500)))

## calculating predicted values of abundance for the 500 depth values based on the gam model
## includes SE associated with predictions
cur.div.pred <- as.data.frame(predict(cur.div.gam, data.frame(cur.div.new), se.fit = T)) %>%
  rename(div_predicted = 1) %>% #renames column 1
  mutate(lci = div_predicted-1.96*se.fit,
         uci = div_predicted+1.96*se.fit) %>%  # upper and lower 95CI
  bind_cols(cur.div.new) # adds depth values

## figure
carib.pal <- fish(4, option = "Bodianus_rufus")

ggplot(cur.div.pred, aes(x = depth, y = div_predicted)) +
  theme_classic()+
  geom_line(lty = 1, color = carib.pal[1]) + # predicted values 
  geom_ribbon(aes(x = depth, ymin = lci, ymax = uci), alpha = 0.2, fill = carib.pal[1]) + # 95 CI
  geom_point(data = div.curacao, aes(x = depth, y = richness), shape = 21,
             alpha = 0.75, color = "grey23", fill = carib.pal[1]) + # original data
  geom_vline(xintercept = 100, lty = 2, color = "grey69") +
  geom_vline(xintercept = 150, lty = 2, color = "grey69") +
  geom_vline(xintercept = 230, lty = 2, color = "grey69") +
  geom_vline(xintercept = 290, lty = 2, color = "grey69") +
  ylab("richness") +
  xlab("depth (m)") +
  scale_x_continuous(limits = c(40,300), breaks = seq(40,300,20))


#### 3.2.2 Evenness ####
cur.div.gam <- gam(evenness ~ s(depth), data = div.curacao) ## here taking the "smooth" of dband.
summary(cur.div.gam) # sign correlation (<0.001) between evenness and depth
par(mfrow = c(2,2))
gam.check(cur.div.gam)

# sub dividing the depth range into 500 values to plot a continuous prediction line
cur.div.new <- data.frame(depth = as.numeric(seq(40,300,length.out = 500)))

## calculating predicted values of abundance for the 500 depth values based on the gam model
## includes SE associated with predictions
cur.div.pred <- as.data.frame(predict(cur.div.gam, data.frame(cur.div.new), se.fit = T)) %>%
  rename(div_predicted = 1) %>% #renames column 1
  mutate(lci = div_predicted-1.96*se.fit,
         uci = div_predicted+1.96*se.fit) %>%  # upper and lower 95CI
  bind_cols(cur.div.new) # adds depth values

## figure
carib.pal <- fish(4, option = "Bodianus_rufus")

ggplot(cur.div.pred, aes(x = depth, y = div_predicted)) +
  theme_classic()+
  geom_line(lty = 1, color = carib.pal[1]) + # predicted values 
  geom_ribbon(aes(x = depth, ymin = lci, ymax = uci), alpha = 0.2, fill = carib.pal[1]) + # 95 CI
  geom_point(data = div.curacao, aes(x = depth, y = evenness), shape = 21,
             alpha = 0.75, color = "grey23", fill = carib.pal[1]) + # original data
  geom_vline(xintercept = 100, lty = 2, color = "grey69") +
  geom_vline(xintercept = 150, lty = 2, color = "grey69") +
  geom_vline(xintercept = 230, lty = 2, color = "grey69") +
  geom_vline(xintercept = 290, lty = 2, color = "grey69") +
  ylab("evenness") +
  xlab("depth (m)") +
  scale_x_continuous(limits = c(40,300), breaks = seq(40,300,20))

#### 3.2.3 Dominance ####
cur.div.gam <- gam(dominance ~ s(depth), data = div.curacao) ## here taking the "smooth" of dband.
summary(cur.div.gam) # sign correlation (<0.001) between evenness and depth
par(mfrow = c(2,2))
gam.check(cur.div.gam)

# sub dividing the depth range into 500 values to plot a continuous prediction line
cur.div.new <- data.frame(depth = as.numeric(seq(40,300,length.out = 500)))

## calculating predicted values of abundance for the 500 depth values based on the gam model
## includes SE associated with predictions
cur.div.pred <- as.data.frame(predict(cur.div.gam, data.frame(cur.div.new), se.fit = T)) %>%
  rename(div_predicted = 1) %>% #renames column 1
  mutate(lci = div_predicted-1.96*se.fit,
         uci = div_predicted+1.96*se.fit) %>%  # upper and lower 95CI
  bind_cols(cur.div.new) # adds depth values

## figure
carib.pal <- fish(4, option = "Bodianus_rufus")

ggplot(cur.div.pred, aes(x = depth, y = div_predicted)) +
  theme_classic()+
  geom_line(lty = 1, color = carib.pal[1]) + # predicted values 
  geom_ribbon(aes(x = depth, ymin = lci, ymax = uci), alpha = 0.2, fill = carib.pal[1]) + # 95 CI
  geom_point(data = div.curacao, aes(x = depth, y = dominance), shape = 21,
             alpha = 0.75, color = "grey23", fill = carib.pal[1]) + # original data
  geom_vline(xintercept = 100, lty = 2, color = "grey69") +
  geom_vline(xintercept = 150, lty = 2, color = "grey69") +
  geom_vline(xintercept = 230, lty = 2, color = "grey69") +
  geom_vline(xintercept = 290, lty = 2, color = "grey69") +
  ylab("dominance") +
  xlab("depth (m)") +
  scale_x_continuous(limits = c(40,300), breaks = seq(40,300,20))


#### 3.2.4. Dissimilarity ####
## transform the abundance matrix to presence/absence matrix 
betadiv.curacao <- hill.curacao %>%
  select(-dband) %>%
  mutate_all(~ case_when(.>0 ~ 1, TRUE~0)) #when abundance>0, abundance becomes 1

## pair-wise dissimilarity between each depth bin + nestdeness and turnover components
beta.pair <- beta.pair(betadiv.curacao,index.family="sorensen")

## run analyses using community clusters to compare instead of individual 10m depth
betadiv.curacao.cluster <- fish.curacao %>%
  mutate(cluster=case_when(dband<=70 ~ "upper mesophotic",
                       dband<=120~ "lower mesophotic",
                       dband<=180~ "upper rariphotic",
                       dband>180~ "lower rariphotic"))%>%
  mutate(cluster=factor(cluster,levels=c("upper mesophotic","lower mesophotic","upper rariphotic","lower rariphotic")))%>%
  group_by(cluster,species) %>%
  summarize(abu.corr=sum(abu.corr)) %>%
  ungroup()%>% #actions will no longer be performed grouped by clusters
  spread(species, abu.corr, fill =  0) %>%
  select(-cluster) %>%
  mutate_all(~ case_when(.>0 ~ 1, TRUE~0))

## beta diversity across clusters
beta.cluster.curacao <- beta.pair(betadiv.curacao.cluster,index.family="sorensen")

## graph: upper mesophotic vs. other depths
taxa.hill0 <- hill_taxa(betadiv.curacao.cluster,q=0)

beta.cluster.curacao.graph <- data.frame(cluster=rep(c("upper mesophotic","lower mesophotic","upper rariphotic","lower rariphotic"),3),
                                         type=rep(c("beta","nestedness","turnover"),each=4),
                                         betadiv=c(0,beta.cluster.curacao$beta.sor[1:3],
                                                     0,beta.cluster.curacao$beta.sne[1:3],
                                                     0,beta.cluster.curacao$beta.sim[1:3]),
                                         richness=taxa.hill0)  %>%
 mutate(cluster=factor(cluster,levels=c("upper mesophotic","lower mesophotic","upper rariphotic","lower rariphotic")))

plot.betadiv.cur <- ggplot(beta.cluster.curacao.graph)+
  geom_point(aes(x=cluster,y=betadiv,shape=type,color=type),position=position_dodge(width=0.3),
             size=5)+
  scale_shape_manual(values=c(16,22,24))+
  scale_color_viridis_d(end=0.7)+
  ggtitle("Curaçao")+
  #theme_classic()+
  theme(panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    text = element_text(color = "black"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10,),
    axis.ticks = element_line(color = "black"),
  plot.margin = unit(c(1,1,1,1),"lines"))+
  coord_cartesian(ylim=c(0,1),clip = "off")+
  xlab("")+
  scale_x_discrete(labels=str_wrap(beta.cluster.curacao.graph$cluster,6))+
  ylab("dissimilarity with upper mesophotic")
  #annotate(x=seq(1,4,1),y=-1,label=taxa.hill0,geom="text")
plot.betadiv.cur

## comparison between all depth bins 
beta.cluster.curacao.graph.2 <- data.frame(cluster=c("lower mesophotic","upper rariphotic","lower rariphotic",
                                                     "upper rariphotic ","lower rariphotic ","lower rariphotic  "),
                                         type=rep(c("beta","nestedness","turnover"),each=6),
                                         betadiv=c(beta.cluster.curacao$beta.sor,
                                                   beta.cluster.curacao$beta.sne,
                                                   beta.cluster.curacao$beta.sim)) %>%
  mutate(cluster=factor(cluster,levels=c("lower mesophotic","upper rariphotic","lower rariphotic",
                                         "upper rariphotic ","lower rariphotic ","lower rariphotic  ")))

plot.betadiv.cur <- ggplot(beta.cluster.curacao.graph.2)+
  geom_point(aes(x=cluster,y=betadiv,shape=type,color=type),position=position_dodge(width=0.3),
             size=5)+
  scale_shape_manual(values=c(16,22,24))+
  scale_color_viridis_d(end=0.7)+
  #theme_classic()+
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        text = element_text(color = "black"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10,),
        axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(4,1,1,1),"lines"))+
  xlab("")+
  geom_vline(xintercept=c(3.5,5.5),linetype="dashed",color="grey")+
  scale_x_discrete(labels=str_wrap(beta.cluster.curacao.graph.2$cluster,6))+
  ylab("dissimilarity")+
annotate("text",x=c(2,4.5,6.2),y=1.1,
         label=str_wrap(c("vs. upper mesophotic","vs. lower mesophotic",
                          "vs. upper rariphotic"),10),size=4)+
  annotate("text",x=0,y=1.2,label="Curaçao")+
  coord_cartesian(ylim=c(0,1),clip ="off")
plot.betadiv.cur

#### 3.3. All repeated for Bonaire ####
## data prep
fish.bonaire <- fish.data %>%
  filter(location == "Bonaire") %>%
  select(species,dband,abu.corr) 

hill.bonaire <- fish.bonaire%>%
  spread(species, abu.corr, fill =  0) 

# add manually the lack of any observations for the 240 m depth bin.
hill.bonaire <- rbind(hill.bonaire,c(240,rep(0,120)))

### taxonomic diversity 
taxa.hill0 <- hill_taxa(hill.bonaire,q=0)
taxa.hill1 <- hill_taxa(hill.bonaire,q=1)
taxa.hill2 <- hill_taxa(hill.bonaire,q=2)

div.bonaire <- data.frame(depth=c(seq(40,230,10),seq(250,300,10),240), richness=taxa.hill0, evenness=taxa.hill1,
                          dominance=taxa.hill2)

### total species richness
ncol(hill.bonaire)

# richness raw
ggplot(div.bonaire,aes(x=depth,y=richness))+
  geom_point()

# evenness - Shannon
ggplot(div.bonaire,aes(x=depth,y=evenness))+
  geom_point()

# dominance- Simpson
ggplot(div.bonaire,aes(x=depth,y=dominance))+
  geom_point()

### Plots with prediction and CI ###
###  Richness ###
bon.div.gam <- gam(richness ~ s(depth), data = div.bonaire) ## here taking the "smooth" of dband.
summary(bon.div.gam) # sign correlation (<0.001) between richness and depth
par(mfrow = c(2,2))
gam.check(bon.div.gam)

# sub dividing the depth range into 500 values to plot a continuous prediction line
bon.div.new <- data.frame(depth = as.numeric(seq(40,300,length.out = 500)))

## calculating predicted values of abundance for the 500 depth values based on the gam model
## includes SE associated with predictions
bon.div.pred <- as.data.frame(predict(bon.div.gam, data.frame(bon.div.new), se.fit = T)) %>%
  rename(div_predicted = 1) %>% #renames column 1
  mutate(lci = div_predicted-1.96*se.fit,
         uci = div_predicted+1.96*se.fit) %>%  # upper and lower 95CI
  bind_cols(bon.div.new) # adds depth values

## figure
carib.pal <- fish(4, option = "Bodianus_rufus")

ggplot(bon.div.pred, aes(x = depth, y = div_predicted)) +
  theme_classic()+
  geom_line(lty = 1, color = carib.pal[1]) + # predicted values 
  geom_ribbon(aes(x = depth, ymin = lci, ymax = uci), alpha = 0.2, fill = carib.pal[1]) + # 95 CI
  geom_point(data = div.bonaire, aes(x = depth, y = richness), shape = 21,
             alpha = 0.75, color = "grey23", fill = carib.pal[1]) + # original data
  geom_vline(xintercept = 100, lty = 2, color = "grey69") +
  geom_vline(xintercept = 150, lty = 2, color = "grey69") +
  geom_vline(xintercept = 230, lty = 2, color = "grey69") +
  geom_vline(xintercept = 290, lty = 2, color = "grey69") +
  ylab("richness") +
  xlab("depth (m)") +
  scale_x_continuous(limits = c(40,300), breaks = seq(40,300,20))

### 3.2.2 Evenness ###
bon.div.gam <- gam(evenness ~ s(depth), data = div.bonaire) ## here taking the "smooth" of dband.
summary(bon.div.gam) # sign correlation (<0.001) between evenness and depth
par(mfrow = c(2,2))
gam.check(bon.div.gam)

# sub dividing the depth range into 500 values to plot a continuous prediction line
bon.div.new <- data.frame(depth = as.numeric(seq(40,300,length.out = 500)))

## calculating predicted values of abundance for the 500 depth values based on the gam model
## includes SE associated with predictions
bon.div.pred <- as.data.frame(predict(bon.div.gam, data.frame(bon.div.new), se.fit = T)) %>%
  rename(div_predicted = 1) %>% #renames column 1
  mutate(lci = div_predicted-1.96*se.fit,
         uci = div_predicted+1.96*se.fit) %>%  # upper and lower 95CI
  bind_cols(bon.div.new) # adds depth values

## figure
carib.pal <- fish(4, option = "Bodianus_rufus")

ggplot(bon.div.pred, aes(x = depth, y = div_predicted)) +
  theme_classic()+
  geom_line(lty = 1, color = carib.pal[1]) + # predicted values 
  geom_ribbon(aes(x = depth, ymin = lci, ymax = uci), alpha = 0.2, fill = carib.pal[1]) + # 95 CI
  geom_point(data = div.bonaire, aes(x = depth, y = evenness), shape = 21,
             alpha = 0.75, color = "grey23", fill = carib.pal[1]) + # original data
  geom_vline(xintercept = 100, lty = 2, color = "grey69") +
  geom_vline(xintercept = 150, lty = 2, color = "grey69") +
  geom_vline(xintercept = 230, lty = 2, color = "grey69") +
  geom_vline(xintercept = 290, lty = 2, color = "grey69") +
  ylab("evenness") +
  xlab("depth (m)") +
  scale_x_continuous(limits = c(40,300), breaks = seq(40,300,20))

### 3.2.3 Dominance ###
bon.div.gam <- gam(dominance ~ s(depth), data = div.bonaire) ## here taking the "smooth" of dband.
summary(bon.div.gam) # sign correlation (<0.001) between evenness and depth
par(mfrow = c(2,2))
gam.check(bon.div.gam)

# sub dividing the depth range into 500 values to plot a continuous prediction line
bon.div.new <- data.frame(depth = as.numeric(seq(40,300,length.out = 500)))

## calculating predicted values of abundance for the 500 depth values based on the gam model
## includes SE associated with predictions
bon.div.pred <- as.data.frame(predict(bon.div.gam, data.frame(bon.div.new), se.fit = T)) %>%
  rename(div_predicted = 1) %>% #renames column 1
  mutate(lci = div_predicted-1.96*se.fit,
         uci = div_predicted+1.96*se.fit) %>%  # upper and lower 95CI
  bind_cols(bon.div.new) # adds depth values

## figure
carib.pal <- fish(4, option = "Bodianus_rufus")

ggplot(bon.div.pred, aes(x = depth, y = div_predicted)) +
  theme_classic()+
  geom_line(lty = 1, color = carib.pal[1]) + # predicted values 
  geom_ribbon(aes(x = depth, ymin = lci, ymax = uci), alpha = 0.2, fill = carib.pal[1]) + # 95 CI
  geom_point(data = div.bonaire, aes(x = depth, y = dominance), shape = 21,
             alpha = 0.75, color = "grey23", fill = carib.pal[1]) + # original data
  geom_vline(xintercept = 100, lty = 2, color = "grey69") +
  geom_vline(xintercept = 150, lty = 2, color = "grey69") +
  geom_vline(xintercept = 230, lty = 2, color = "grey69") +
  geom_vline(xintercept = 290, lty = 2, color = "grey69") +
  ylab("dominance") +
  xlab("depth (m)") +
  scale_x_continuous(limits = c(40,300), breaks = seq(40,300,20))

#### Dissimilarity 

betadiv.bonaire.cluster <- fish.bonaire %>%
  mutate(cluster=case_when(dband<70 ~ "upper mesophotic",
                       dband<120~ "lower mesophotic",
                       dband<170~ "upper rariphotic",
                       dband>=170~ "lower rariphotic")) %>%
  mutate(cluster=factor(cluster,levels=c("upper mesophotic","lower mesophotic","upper rariphotic","lower rariphotic")))%>%
  group_by(cluster,species) %>%
  summarize(abu.corr=sum(abu.corr)) %>%
  ungroup()%>% #actions will no longer be performed grouped by clusters
  spread(species, abu.corr, fill =  0) %>%
  select(-cluster) %>%
  mutate_all(~ case_when(.>0 ~ 1, TRUE~0))

beta.cluster.bonaire <- beta.pair(betadiv.bonaire.cluster,index.family="sorensen")

## species richness per cluster
taxa.hill0 <- hill_taxa(betadiv.bonaire.cluster,q=0)

beta.cluster.bonaire.graph <- data.frame(cluster=c("lower mesophotic","upper rariphotic","lower rariphotic",
                                                     "upper rariphotic ","lower rariphotic ","lower rariphotic  "),
                                           type=rep(c("beta","nestedness","turnover"),each=6),
                                           betadiv=c(beta.cluster.bonaire$beta.sor,
                                                     beta.cluster.bonaire$beta.sne,
                                                     beta.cluster.bonaire$beta.sim)) %>%
  mutate(cluster=factor(cluster,levels=c("lower mesophotic","upper rariphotic","lower rariphotic",
                                         "upper rariphotic ","lower rariphotic ","lower rariphotic  ")))

plot.betadiv.bon <- ggplot(beta.cluster.bonaire.graph)+
  geom_point(aes(x=cluster,y=betadiv,shape=type,color=type),position=position_dodge(width=0.3),
             size=5)+
  scale_shape_manual(values=c(16,22,24))+
  scale_color_viridis_d(end=0.7)+
  #theme_classic()+
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        text = element_text(color = "black"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10,),
        axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(4,1,1,1),"lines"))+
  xlab("")+
  geom_vline(xintercept=c(3.5,5.5),linetype="dashed",color="grey")+
  scale_x_discrete(labels=str_wrap(beta.cluster.bonaire.graph$cluster,6))+
  ylab("dissimilarity")+
  annotate("text",x=c(2,4.5,6.2),y=1.1,
           label=str_wrap(c("vs. upper mesophotic","vs. lower mesophotic",
                            "vs. upper rariphotic"),10),size=4)+
  annotate("text",x=0,y=1.2,label="Bonaire")+
  coord_cartesian(ylim=c(0,1),clip ="off")
plot.betadiv.bon


#### 3.4. All repeated for Statia ####
## data prep
fish.statia <- fish.data %>%
  filter(location == "St. Eustatius") %>%
  select(species,dband,abu.corr) 

hill.statia <- fish.statia%>%
  spread(species, abu.corr, fill =  0) 

# add manually the lack of any observations for the 250 m depth bin.
hill.statia <- rbind(hill.statia,c(250,rep(0,120)))

### taxonomic diversity 
taxa.hill0 <- hill_taxa(hill.statia,q=0)
taxa.hill1 <- hill_taxa(hill.statia,q=1)
taxa.hill2 <- hill_taxa(hill.statia,q=2)

## have to add three lines where spric = 0
div.statia <- data.frame(depth=c(seq(40,240,10),250,260,270), richness=taxa.hill0, evenness=taxa.hill1,
                          dominance=taxa.hill2)

## total species richness
ncol(hill.statia)

# richness raw
ggplot(div.statia,aes(x=depth,y=richness))+
  geom_point()

# evenness - Shannon
ggplot(div.statia,aes(x=depth,y=evenness))+
  geom_point()

# dominance- Simpson
ggplot(div.statia,aes(x=depth,y=dominance))+
  geom_point()

### Plots with prediction and CI ###
###  Richness ###
sta.div.gam <- gam(richness ~ s(depth), data = div.statia) ## here taking the "smooth" of dband.
summary(sta.div.gam) # sign correlation (<0.001) between richness and depth
par(mfrow = c(2,2))
gam.check(sta.div.gam)

# sub dividing the depth range into 500 values to plot a continuous prediction line
sta.div.new <- data.frame(depth = as.numeric(seq(40,300,length.out = 500)))

## calculating predicted values of abundance for the 500 depth values based on the gam model
## includes SE associated with predictions
sta.div.pred <- as.data.frame(predict(sta.div.gam, data.frame(sta.div.new), se.fit = T)) %>%
  rename(div_predicted = 1) %>% #renames column 1
  mutate(lci = div_predicted-1.96*se.fit,
         uci = div_predicted+1.96*se.fit) %>%  # upper and lower 95CI
  bind_cols(sta.div.new) # adds depth values

## figure
carib.pal <- fish(4, option = "Bodianus_rufus")

ggplot(sta.div.pred, aes(x = depth, y = div_predicted)) +
  theme_classic()+
  geom_line(lty = 1, color = carib.pal[3]) + # predicted values 
  geom_ribbon(aes(x = depth, ymin = lci, ymax = uci), alpha = 0.2, fill = carib.pal[3]) + # 95 CI
  geom_point(data = div.statia, aes(x = depth, y = richness), shape = 21,
             alpha = 0.75, color = "grey23", fill = carib.pal[3]) + # original data
  geom_vline(xintercept = 100, lty = 2, color = "grey69") +
  geom_vline(xintercept = 150, lty = 2, color = "grey69") +
  geom_vline(xintercept = 230, lty = 2, color = "grey69") +
  geom_vline(xintercept = 290, lty = 2, color = "grey69") +
  ylab("richness") +
  xlab("depth (m)") +
  scale_x_continuous(limits = c(40,300), breaks = seq(40,300,20))

### Evenness ###
sta.div.gam <- gam(evenness ~ s(depth), data = div.statia) ## here taking the "smooth" of dband.
summary(sta.div.gam) # sign correlation (<0.001) between evenness and depth
par(mfrow = c(2,2))
gam.check(sta.div.gam)

# sub dividing the depth range into 500 values to plot a continuous prediction line
sta.div.new <- data.frame(depth = as.numeric(seq(40,300,length.out = 500)))

## calculating predicted values of abundance for the 500 depth values based on the gam model
## includes SE associated with predictions
sta.div.pred <- as.data.frame(predict(sta.div.gam, data.frame(sta.div.new), se.fit = T)) %>%
  rename(div_predicted = 1) %>% #renames column 1
  mutate(lci = div_predicted-1.96*se.fit,
         uci = div_predicted+1.96*se.fit) %>%  # upper and lower 95CI
  bind_cols(sta.div.new) # adds depth values

## figure
ggplot(sta.div.pred, aes(x = depth, y = div_predicted)) +
  theme_classic()+
  geom_line(lty = 1, color = carib.pal[3]) + # predicted values 
  geom_ribbon(aes(x = depth, ymin = lci, ymax = uci), alpha = 0.2, fill = carib.pal[3]) + # 95 CI
  geom_point(data = div.statia, aes(x = depth, y = evenness), shape = 21,
             alpha = 0.75, color = "grey23", fill = carib.pal[3]) + # original data
  geom_vline(xintercept = 100, lty = 2, color = "grey69") +
  geom_vline(xintercept = 150, lty = 2, color = "grey69") +
  geom_vline(xintercept = 230, lty = 2, color = "grey69") +
  geom_vline(xintercept = 290, lty = 2, color = "grey69") +
  ylab("evenness") +
  xlab("depth (m)") +
  scale_x_continuous(limits = c(40,300), breaks = seq(40,300,20))

### Dominance ###
sta.div.gam <- gam(dominance ~ s(depth), data = div.statia) ## here taking the "smooth" of dband.
summary(sta.div.gam) # sign correlation (<0.001) between evenness and depth
par(mfrow = c(2,2))
gam.check(sta.div.gam)

# sub dividing the depth range into 500 values to plot a continuous prediction line
sta.div.new <- data.frame(depth = as.numeric(seq(40,300,length.out = 500)))

## calculating predicted values of abundance for the 500 depth values based on the gam model
## includes SE associated with predictions
sta.div.pred <- as.data.frame(predict(sta.div.gam, data.frame(sta.div.new), se.fit = T)) %>%
  rename(div_predicted = 1) %>% #renames column 1
  mutate(lci = div_predicted-1.96*se.fit,
         uci = div_predicted+1.96*se.fit) %>%  # upper and lower 95CI
  bind_cols(sta.div.new) # adds depth values

## figure
carib.pal <- fish(4, option = "Bodianus_rufus")

ggplot(sta.div.pred, aes(x = depth, y = div_predicted)) +
  theme_classic()+
  geom_line(lty = 1, color = carib.pal[3]) + # predicted values 
  geom_ribbon(aes(x = depth, ymin = lci, ymax = uci), alpha = 0.2, fill = carib.pal[3]) + # 95 CI
  geom_point(data = div.statia, aes(x = depth, y = dominance), shape = 21,
             alpha = 0.75, color = "grey23", fill = carib.pal[3]) + # original data
  geom_vline(xintercept = 100, lty = 2, color = "grey69") +
  geom_vline(xintercept = 150, lty = 2, color = "grey69") +
  geom_vline(xintercept = 230, lty = 2, color = "grey69") +
  geom_vline(xintercept = 290, lty = 2, color = "grey69") +
  ylab("dominance") +
  xlab("depth (m)") +
  scale_x_continuous(limits = c(40,300), breaks = seq(40,300,20))

#### Dissimilarity 
betadiv.statia.cluster <- fish.statia %>%
  mutate(cluster=case_when(dband<100 ~ "upper mesophotic",
                       dband<140 ~ "lower mesophotic",
                       dband<190~ "upper rariphotic",
                       dband>=190~ "lower rariphotic")) %>%
  mutate(cluster=factor(cluster,levels=c("upper mesophotic","lower mesophotic",
                                         "upper rariphotic","lower rariphotic"))) %>%
  group_by(cluster,species) %>%
  summarize(abu.corr=sum(abu.corr)) %>%
  ungroup()%>% #actions will no longer be performed grouped by clusters
  spread(species, abu.corr, fill =  0) %>%
  select(-cluster) %>%
  mutate_all(~ case_when(.>0 ~ 1, TRUE~0))

beta.cluster.statia <- beta.pair(betadiv.statia.cluster,index.family="sorensen")
## species richness per cluster
beta.cluster.statia.graph <- data.frame(cluster=c("lower mesophotic","upper rariphotic","lower rariphotic",
                                                   "upper rariphotic ","lower rariphotic ","lower rariphotic  "),
                                         type=rep(c("beta","nestedness","turnover"),each=6),
                                         betadiv=c(beta.cluster.statia$beta.sor,
                                                   beta.cluster.statia$beta.sne,
                                                   beta.cluster.statia$beta.sim)) %>%
  mutate(cluster=factor(cluster,levels=c("lower mesophotic","upper rariphotic","lower rariphotic",
                                         "upper rariphotic ","lower rariphotic ","lower rariphotic  ")))

plot.betadiv.sta <- ggplot(beta.cluster.statia.graph)+
  geom_point(aes(x=cluster,y=betadiv,shape=type,color=type),position=position_dodge(width=0.3),
             size=5)+
  scale_shape_manual(values=c(16,22,24))+
  scale_color_viridis_d(end=0.7)+
  #theme_classic()+
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        text = element_text(color = "black"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10,),
        axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(4,1,1,1),"lines"))+
  xlab("")+
  geom_vline(xintercept=c(3.5,5.5),linetype="dashed",color="grey")+
  scale_x_discrete(labels=str_wrap(beta.cluster.statia.graph$cluster,6))+
  ylab("dissimilarity")+
  annotate("text",x=c(2,4.5,6.2),y=1.1,
           label=str_wrap(c("vs. upper mesophotic","vs. lower mesophotic",
                            "vs. upper rariphotic"),10),size=4)+
  annotate("text",x=0,y=1.2,label="Statia")+
  coord_cartesian(ylim=c(0,1),clip ="off")
plot.betadiv.sta

#### 3.5 All repeated for Roatan ####
## data prep
fish.roatan <- fish.data %>%
  filter(location == "Roatan") %>%
  filter(dband > 30 & dband < 310) %>%
  select(species,dband,abu.corr)

hill.roatan <- fish.roatan %>%
  spread(species, abu.corr, fill =  0) 


### taxonomic diversity 
taxa.hill0 <- hill_taxa(hill.roatan,q=0)
taxa.hill1 <- hill_taxa(hill.roatan,q=1)
taxa.hill2 <- hill_taxa(hill.roatan,q=2)


div.roatan <- data.frame(depth=sort(unique(fish.roatan$dband)), richness=taxa.hill0, evenness=taxa.hill1,
                         dominance=taxa.hill2)

## total species richness
ncol(hill.roatan)

# richness raw
ggplot(div.roatan,aes(x=depth,y=richness))+
  geom_point()

# evenness - Shannon
ggplot(div.roatan,aes(x=depth,y=evenness))+
  geom_point()

# dominance- Simpson
ggplot(div.roatan,aes(x=depth,y=dominance))+
  geom_point()

### Plots with prediction and CI ###
###  Richness ###
roa.div.gam <- gam(richness ~ s(depth), data = div.roatan) ## here taking the "smooth" of dband.
summary(roa.div.gam) # sign correlation (<0.001) between richness and depth
par(mfrow = c(2,2))
gam.check(roa.div.gam)

# sub dividing the depth range into 500 values to plot a continuous prediction line
roa.div.new <- data.frame(depth = as.numeric(seq(10,500,length.out = 800)))

## calculating predicted values of abundance for the 500 depth values based on the gam model
## includes SE associated with predictions
roa.div.pred <- as.data.frame(predict(roa.div.gam, data.frame(roa.div.new), se.fit = T)) %>%
  rename(div_predicted = 1) %>% #renames column 1
  mutate(lci = div_predicted-1.96*se.fit,
         uci = div_predicted+1.96*se.fit) %>%  # upper and lower 95CI
  bind_cols(roa.div.new) # adds depth values

## figure
carib.pal <- fish(4, option = "Bodianus_rufus")

ggplot(roa.div.pred, aes(x = depth, y = div_predicted)) +
  theme_classic()+
  geom_line(lty = 1, color = carib.pal[4]) + # predicted values 
  geom_ribbon(aes(x = depth, ymin = lci, ymax = uci), alpha = 0.2, fill = carib.pal[4]) + # 95 CI
  geom_point(data = div.roatan, aes(x = depth, y = richness), shape = 21,
             alpha = 0.75, color = "grey23", fill = carib.pal[4]) + # original data
  geom_vline(xintercept = 100, lty = 2, color = "grey69") +
  geom_vline(xintercept = 150, lty = 2, color = "grey69") +
  geom_vline(xintercept = 230, lty = 2, color = "grey69") +
  geom_vline(xintercept = 290, lty = 2, color = "grey69") +
  ylab("richness") +
  xlab("depth (m)") +
  scale_x_continuous(limits = c(10,500), breaks = seq(10,500,50))

### Evenness ###
roa.div.gam <- gam(evenness ~ s(depth), data = div.roatan) ## here taking the "smooth" of dband.
summary(roa.div.gam) # sign correlation (<0.001) between richness and depth
par(mfrow = c(2,2))
gam.check(roa.div.gam)

# sub dividing the depth range into 500 values to plot a continuous prediction line
roa.div.new <- data.frame(depth = as.numeric(seq(10,500,length.out = 800)))

## calculating predicted values of abundance for the 500 depth values based on the gam model
## includes SE associated with predictions
roa.div.pred <- as.data.frame(predict(roa.div.gam, data.frame(roa.div.new), se.fit = T)) %>%
  rename(div_predicted = 1) %>% #renames column 1
  mutate(lci = div_predicted-1.96*se.fit,
         uci = div_predicted+1.96*se.fit) %>%  # upper and lower 95CI
  bind_cols(roa.div.new) # adds depth values

## figure
carib.pal <- fish(4, option = "Bodianus_rufus")

ggplot(roa.div.pred, aes(x = depth, y = div_predicted)) +
  theme_classic()+
  geom_line(lty = 1, color = carib.pal[4]) + # predicted values 
  geom_ribbon(aes(x = depth, ymin = lci, ymax = uci), alpha = 0.2, fill = carib.pal[4]) + # 95 CI
  geom_point(data = div.roatan, aes(x = depth, y = evenness), shape = 21,
             alpha = 0.75, color = "grey23", fill = carib.pal[4]) + # original data
  geom_vline(xintercept = 100, lty = 2, color = "grey69") +
  geom_vline(xintercept = 150, lty = 2, color = "grey69") +
  geom_vline(xintercept = 230, lty = 2, color = "grey69") +
  geom_vline(xintercept = 290, lty = 2, color = "grey69") +
  ylab("evenness") +
  xlab("depth (m)") +
  scale_x_continuous(limits = c(10,500), breaks = seq(10,500,50))

### Dominance ###
roa.div.gam <- gam(dominance ~ s(depth), data = div.roatan) ## here taking the "smooth" of dband.
summary(roa.div.gam) # sign correlation (<0.001) between richness and depth
par(mfrow = c(2,2))
gam.check(roa.div.gam)

# sub dividing the depth range into 500 values to plot a continuous prediction line
roa.div.new <- data.frame(depth = as.numeric(seq(10,500,length.out = 800)))

## calculating predicted values of abundance for the 500 depth values based on the gam model
## includes SE associated with predictions
roa.div.pred <- as.data.frame(predict(roa.div.gam, data.frame(roa.div.new), se.fit = T)) %>%
  rename(div_predicted = 1) %>% #renames column 1
  mutate(lci = div_predicted-1.96*se.fit,
         uci = div_predicted+1.96*se.fit) %>%  # upper and lower 95CI
  bind_cols(roa.div.new) # adds depth values

## figure
carib.pal <- fish(4, option = "Bodianus_rufus")

ggplot(roa.div.pred, aes(x = depth, y = div_predicted)) +
  theme_classic()+
  geom_line(lty = 1, color = carib.pal[4]) + # predicted values 
  geom_ribbon(aes(x = depth, ymin = lci, ymax = uci), alpha = 0.2, fill = carib.pal[4]) + # 95 CI
  geom_point(data = div.roatan, aes(x = depth, y = dominance), shape = 21,
             alpha = 0.75, color = "grey23", fill = carib.pal[4]) + # original data
  geom_vline(xintercept = 100, lty = 2, color = "grey69") +
  geom_vline(xintercept = 150, lty = 2, color = "grey69") +
  geom_vline(xintercept = 230, lty = 2, color = "grey69") +
  geom_vline(xintercept = 290, lty = 2, color = "grey69") +
  ylab("dominance") +
  xlab("depth (m)") +
  scale_x_continuous(limits = c(10,500), breaks = seq(10,500,50))

#### Dissimilarity 

betadiv.roatan.cluster <- fish.roatan %>%
  filter(dband<= 300)%>%
  mutate(cluster=case_when(dband<100 ~ "upper mesophotic",
                       dband<150 ~ "lower mesophotic",
                       dband<210~ "upper rariphotic",
                       dband>=210~ "lower rariphotic")) %>%
  mutate(cluster=factor(cluster,levels=c("upper mesophotic","lower mesophotic",
                                         "upper rariphotic","lower rariphotic"))) %>%
  group_by(cluster,species) %>%
  summarize(abu.corr=sum(abu.corr)) %>%
  ungroup()%>% #actions will no longer be performed grouped by clusters
  spread(species, abu.corr, fill =  0) %>%
  select(-cluster) %>%
  mutate_all(~ case_when(.>0 ~ 1, TRUE~0))

beta.cluster.roatan <- beta.pair(betadiv.roatan.cluster,index.family="sorensen")

beta.cluster.roatan.graph <- data.frame(cluster=c("lower mesophotic","upper rariphotic","lower rariphotic",
                                                  "upper rariphotic ","lower rariphotic ","lower rariphotic  "),
                                        type=rep(c("beta","nestedness","turnover"),each=6),
                                        betadiv=c(beta.cluster.roatan$beta.sor,
                                                  beta.cluster.roatan$beta.sne,
                                                  beta.cluster.roatan$beta.sim)) %>%
  mutate(cluster=factor(cluster,levels=c("lower mesophotic","upper rariphotic","lower rariphotic",
                                         "upper rariphotic ","lower rariphotic ","lower rariphotic  ")))

plot.betadiv.roa <- ggplot(beta.cluster.roatan.graph)+
  geom_point(aes(x=cluster,y=betadiv,shape=type,color=type),position=position_dodge(width=0.3),
             size=5)+
  scale_shape_manual(values=c(16,22,24))+
  scale_color_viridis_d(end=0.7)+
  #theme_classic()+
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        text = element_text(color = "black"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10,),
        axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(4,1,1,1),"lines"))+
  xlab("")+
  geom_vline(xintercept=c(3.5,5.5),linetype="dashed",color="grey")+
  scale_x_discrete(labels=str_wrap(beta.cluster.roatan.graph$cluster,6))+
  ylab("dissimilarity")+
  annotate("text",x=c(2,4.5,6.2),y=1.1,
           label=str_wrap(c("vs. upper mesophotic","vs. lower mesophotic",
                            "vs. upper rariphotic"),10),size=4)+
  annotate("text",x=0,y=1.2,label="Roatan")+
  coord_cartesian(ylim=c(0,1),clip ="off")
plot.betadiv.roa

#### 3.6. Averages across sites ####
betadiv.roatan.fullrange <- read.csv("2.Ch3.R.analysis/betadiv.roatan.fullrange.csv") %>%
  select(-X)%>%
  mutate(se=c(0))
colnames(betadiv.roatan.fullrange) <- c("cluster","type", "avg","se")

beta.cluster.avg <- beta.cluster.roatan.graph %>%
  left_join(beta.cluster.bonaire.graph, by=c("cluster","type")) %>%
  left_join(beta.cluster.curacao.graph.2, by=c("cluster","type")) %>%
  left_join(beta.cluster.statia.graph, by=c("cluster","type")) %>%
  rowwise() %>%
  mutate(avg=mean(c(betadiv.x,betadiv.y,betadiv.x.x,betadiv.y.y)), 
         se=sd(c(betadiv.x,betadiv.y,betadiv.x.x,betadiv.y.y))/2) %>%
  select(cluster,type,avg,se) %>%
  rbind(betadiv.roatan.fullrange[c(1,15,16,30,31,45),]) 
beta.cluster.avg$cluster <- as.character(beta.cluster.avg$cluster)
beta.cluster.avg$cluster <- factor(beta.cluster.avg$cluster, levels= c("upper mesophotic","lower mesophotic","upper rariphotic","lower rariphotic",
        "upper rariphotic ","lower rariphotic ","lower rariphotic  ", "deep sea"))

labels=c("upper mesophotic","lower mesophotic","upper rariphotic","lower rariphotic",
        "upper rariphotic ","lower rariphotic ","lower rariphotic  ", "deep sea")

ggplot(beta.cluster.avg)+
  geom_point(data=beta.cluster.avg,aes(x=cluster,y=avg,shape=type,color=type),position=position_dodge(width=0),
             size=5.5)+
  geom_errorbar(aes(x=cluster,y=avg,ymin=avg-se, ymax=avg+se), width=.05,
                position=position_dodge(.9),color="black") +
  scale_shape_manual(values=c(15,16,17))+
  #scale_color_viridis_d(end=0.7)+
  scale_color_manual(values=c("#B26969","#FFDEBC","#E7A72B"))+
  #theme_classic()+
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        text = element_text(color = "black"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10,),
        axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(1,1,5,1),"lines"))+
  xlab("")+
  geom_vline(xintercept=c(1.5,4.5,6.5,7.5),linetype="dashed",color="grey")+
  scale_x_discrete(labels=str_wrap(labels,6))+
  ylab("dissimilarity")+
  annotate("text",x=c(1,3,5.5,7,8),y=-0.2,
          label=str_wrap(c("vs. altiphotic","vs. upper mesophotic","vs. lower mesophotic",
                           "vs. upper rariphotic","vs. lower rariphotic"),10),size=3.5)+
coord_cartesian(ylim=c(0,1),clip ="off")

#### Figure 5: Turnover and nestdesness between depth ####
# Need to run all section 3
beta.cluster.avg <- beta.cluster.avg %>%
  filter(type!="beta")

ggplot(beta.cluster.avg)+
  geom_point(data=beta.cluster.avg,aes(x=cluster,y=avg,shape=type,color=type),position=position_dodge(width=0),
             size=5.5)+
  geom_errorbar(aes(x=cluster,y=avg,ymin=avg-se, ymax=avg+se), width=.05,
                position=position_dodge(.9),color="black") +
  scale_shape_manual(values=c(15,16,17))+
  #scale_color_viridis_d(end=0.7)+
  scale_color_manual(values=c("#B26969","#FFDEBC","#E7A72B"))+
  #theme_classic()+
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        text = element_text(color = "black"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10,),
        axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(1,1,5,1),"lines"))+
  xlab("")+
  geom_vline(xintercept=c(1.5,4.5,6.5,7.5),linetype="dashed",color="grey")+
  scale_x_discrete(labels=str_wrap(labels,6))+
  ylab("dissimilarity")+
  annotate("text",x=c(1,3,5.5,7,8),y=-0.2,
           label=str_wrap(c("vs. altiphotic","vs. upper mesophotic","vs. lower mesophotic",
                            "vs. upper rariphotic","vs. lower rariphotic"),10),size=3.5)+
  coord_cartesian(ylim=c(0,1),clip ="off")



#### 4. Between site beta diversity ####
## Need to run all the 3. sections to run this part ##
betadiv.clusters <- rbind(betadiv.statia.cluster,betadiv.curacao.cluster,
                          betadiv.bonaire.cluster,betadiv.roatan.cluster)

betadiv.clusters <- fish.data%>%
  select(species,dband,location,abu.corr)%>%
  filter(dband>30 & dband<310)%>%
  mutate(cluster=case_when(dband<=70 & location=="Curacao"~ "upper mesophotic",
                           dband<=120 & location=="Curacao"~ "lower mesophotic",
                           dband<=180 & location=="Curacao"~ "upper rariphotic",
                           dband>180 & location=="Curacao"~ "lower rariphotic",
                      dband<100 & location=="Roatan"~ "upper mesophotic",
                     dband<150 & location=="Roatan"~ "lower mesophotic",
                     dband<210 & location=="Roatan"~ "upper rariphotic",
                     dband>=210 & location=="Roatan"~ "lower rariphotic",
                     dband<100 & location=="St. Eustatius"~ "upper mesophotic",
                     dband<140 & location=="St. Eustatius"~ "lower mesophotic",
                     dband<190 & location=="St. Eustatius"~ "upper rariphotic",
                     dband>=190 & location=="St. Eustatius"~ "lower rariphotic",
                     dband<70 & location=="Bonaire"~ "upper mesophotic",
                     dband<120 & location=="Bonaire"~ "lower mesophotic",
                     dband<170 & location=="Bonaire"~ "upper rariphotic",
                     dband>=170 & location=="Bonaire"~ "lower rariphotic")) %>%
  #mutate(cluster=factor(cluster,levels=c("upper mesophotic","lower mesophotic","upper rariphotic","lower rariphotic")))%>%
  group_by(cluster,species,location) %>%
  summarize(abu.corr=sum(abu.corr))

#### upper mesophotic 
betadiv.upMesoP <- betadiv.clusters %>%
  ungroup() %>%
  filter(cluster=="upper mesophotic") %>%
  spread(species, abu.corr, fill =  0) %>%
  mutate(location=factor(location,levels=c("Curacao","Bonaire","St. Eustatius","Roatan"))) %>%
  arrange(location) %>%
  select(-cluster,-location) %>%
  mutate_all(~ case_when(.>0 ~ 1, TRUE~0))

beta.cluster.mesoP <- beta.pair(betadiv.upMesoP,index.family="sorensen")

plot.mesoP <- data_frame(sites1=c(rep("Curacao",3),rep("Bonaire",2),"St. Eustatius"),
                        sites2=c("Bonaire","St. Eustatius","Roatan","St. Eustatius","Roatan","Roatan"),
                        betadiv=beta.cluster.mesoP$beta.sor)%>%
  mutate(sites1=factor(sites1,levels=c("Curacao","Bonaire","St. Eustatius")),
         sites2=factor(sites2,levels=rev(c("Bonaire","St. Eustatius","Roatan"))))
  

plot.beta.upMeso <- ggplot(data=plot.mesoP, aes(x=sites2, y=sites1,fill=betadiv))+
  geom_tile(color="grey")+
  scale_fill_viridis(option="magma",begin=1,end=0.55,limits=c(0.4,0.85))+
  coord_flip()+
  scale_y_discrete(position="right")+
  theme(panel.grid.major = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))+
  labs(fill="dissimilarity",y="",x="")+
  ggtitle("upper mesophotic")

#### lower mesophotic 
betadiv.lowMesoP <- betadiv.clusters %>%
  ungroup() %>%
  filter(cluster=="lower mesophotic") %>%
  spread(species, abu.corr, fill =  0) %>%
  mutate(location=factor(location,levels=c("Curacao","Bonaire","St. Eustatius","Roatan"))) %>%
  arrange(location) %>%
  select(-cluster,-location) %>%
  mutate_all(~ case_when(.>0 ~ 1, TRUE~0))

beta.cluster.mesoP <- beta.pair(betadiv.lowMesoP,index.family="sorensen")

plot.mesoP <- data_frame(sites1=c(rep("Curacao",3),rep("Bonaire",2),"St. Eustatius"),
                         sites2=c("Bonaire","St. Eustatius","Roatan","St. Eustatius","Roatan","Roatan"),
                         betadiv=beta.cluster.mesoP$beta.sor)%>%
  mutate(sites1=factor(sites1,levels=c("Curacao","Bonaire","St. Eustatius")),
         sites2=factor(sites2,levels=rev(c("Bonaire","St. Eustatius","Roatan"))))



plot.beta.lowMeso <-ggplot(data=plot.mesoP, aes(x=sites2, y=sites1,fill=betadiv))+
  geom_tile(color="grey")+
  scale_fill_viridis(option="magma",begin=1,end=0.55,limits=c(0.4,0.85))+
  coord_flip()+
  scale_y_discrete(position="right")+
  theme(panel.grid.major = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))+
  labs(fill="dissimilarity",y="",x="")+
  ggtitle("lower mesophotic")

#### upper rariphotic
betadiv.upRariP <- betadiv.clusters %>%
  ungroup() %>%
  filter(cluster=="upper rariphotic") %>%
  spread(species, abu.corr, fill =  0) %>%
  mutate(location=factor(location,levels=c("Curacao","Bonaire","St. Eustatius","Roatan"))) %>%
  arrange(location) %>%
  select(-cluster,-location) %>%
  mutate_all(~ case_when(.>0 ~ 1, TRUE~0))

beta.cluster.rariP <- beta.pair(betadiv.upRariP,index.family="sorensen")

plot.rariP <- data_frame(sites1=c(rep("Curacao",3),rep("Bonaire",2),"St. Eustatius"),
                         sites2=c("Bonaire","St. Eustatius","Roatan","St. Eustatius","Roatan","Roatan"),
                         betadiv=beta.cluster.rariP$beta.sor)%>%
  mutate(sites1=factor(sites1,levels=c("Curacao","Bonaire","St. Eustatius")),
         sites2=factor(sites2,levels=rev(c("Bonaire","St. Eustatius","Roatan"))))



plot.beta.upRari <- ggplot(data=plot.rariP, aes(x=sites2, y=sites1,fill=betadiv))+
  geom_tile(color="grey")+
  scale_fill_viridis(option="magma",begin=1,end=0.55,limits=c(0.4,0.85))+
  coord_flip()+
  scale_y_discrete(position="right")+
  theme(panel.grid.major = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))+
  labs(fill="dissimilarity",y="",x="")+
  ggtitle("upper rariphotic")

#### lower rariphotic
betadiv.lowRariP <- betadiv.clusters %>%
  ungroup() %>%
  filter(cluster=="lower rariphotic") %>%
  spread(species, abu.corr, fill =  0) %>%
  mutate(location=factor(location,levels=c("Curacao","Bonaire","St. Eustatius","Roatan"))) %>%
  arrange(location) %>%
  select(-cluster,-location) %>%
  mutate_all(~ case_when(.>0 ~ 1, TRUE~0))

beta.cluster.rariP <- beta.pair(betadiv.lowRariP,index.family="sorensen")

plot.rariP <- data_frame(sites1=c(rep("Curacao",3),rep("Bonaire",2),"St. Eustatius"),
                         sites2=c("Bonaire","St. Eustatius","Roatan","St. Eustatius","Roatan","Roatan"),
                         betadiv=beta.cluster.rariP$beta.sor)%>%
  mutate(sites1=factor(sites1,levels=c("Curacao","Bonaire","St. Eustatius")),
         sites2=factor(sites2,levels=rev(c("Bonaire","St. Eustatius","Roatan"))))



plot.beta.lowRari <-ggplot(data=plot.rariP, aes(x=sites2, y=sites1,fill=betadiv))+
  geom_tile(color="grey")+
  scale_fill_viridis(option="magma",begin=1,end=0.55,limits=c(0.4,0.85))+
  coord_flip()+
  scale_y_discrete(position="right")+
  theme(panel.grid.major = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))+
  labs(fill="dissimilarity",y="",x="")+
  ggtitle("lower rariphotic")

transsite.beta.plots <- ggarrange(plot.beta.upMeso, plot.beta.lowMeso+ rremove("ylab"),plot.beta.upRari+ rremove("ylab"),plot.beta.lowRari+ rremove("ylab"),
                           ncol = 2, nrow=2, common.legend=TRUE,legend="bottom")
transsite.beta.plots


#### 5. Rarefaction curves ####
#### Curacao ####
rarefaction.cur <- read.csv(file = "raw fish counts/Curacao_Raw_2020.csv") 
rarefaction.cur <- rarefaction.cur[sample(nrow(rarefaction.cur)), ] #shuffles the ordre of lines

total.species <- vector(length=nrow(rarefaction.cur))
for (k in 1:nrow(rarefaction.cur)){
  tot=length(unique(rarefaction.cur$name[1:k]))
  total.species[k] <- tot
}

ggplot()+
  geom_point(aes(x=1:nrow(rarefaction.cur),y=total.species))+
  ylab("total species richness")+
  xlab("total number of observations")+
  ggtitle("Curaçao - all depths")

## only for deep species
rarefaction.cur.deep <-rarefaction.cur %>%
  filter(meters<150)

total.species.deep <- vector(length=nrow(rarefaction.cur.deep))
for (k in 1:nrow(rarefaction.cur.deep)){
  tot=length(unique(rarefaction.cur.deep$name[1:k]))
  total.species.deep[k] <- tot
}

ggplot()+
  geom_point(aes(x=1:nrow(rarefaction.cur.deep),y=total.species.deep))+
  ylab("total species richness")+
  xlab("total number of observations")+
  ggtitle("Curaçao - deep")

#### Bonaire ####
rarefaction.bon <- read.csv(file = "raw fish counts/Bonaire_Raw_2020.csv") 
rarefaction.bon <- rarefaction.bon[sample(nrow(rarefaction.bon)), ] #shuffles the ordre of lines

total.species <- vector(length=nrow(rarefaction.bon))
for (k in 1:nrow(rarefaction.bon)){
  tot=length(unique(rarefaction.bon$name[1:k]))
  total.species[k] <- tot
}

ggplot()+
  geom_point(aes(x=1:nrow(rarefaction.bon),y=total.species))+
  ylab("total species richness")+
  xlab("total number of observations")+
  ggtitle("Bonaire - all depths")

#### Statia ####
rarefaction.sta <- read.csv(file = "raw fish counts/Statia_Raw_2020.csv") 
rarefaction.sta <- rarefaction.sta[sample(nrow(rarefaction.sta)), ] #shuffles the ordre of lines

total.species <- vector(length=nrow(rarefaction.sta))
for (k in 1:nrow(rarefaction.sta)){
  tot=length(unique(rarefaction.sta$name[1:k]))
  total.species[k] <- tot
}

ggplot()+
  geom_point(aes(x=1:nrow(rarefaction.sta),y=total.species))+
  ylab("total species richness")+
  xlab("total number of observations")+
  ggtitle("Statia - all depths")

#### Roatan ####
rarefaction.roa <- read.csv(file = "raw fish counts/Roatan_Raw_2020.csv") %>%
  filter(meters>40 & meters<310)
rarefaction.roa <- rarefaction.roa[sample(nrow(rarefaction.roa)), ] #shuffles the ordre of lines

total.species <- vector(length=nrow(rarefaction.roa))
for (k in 1:nrow(rarefaction.roa)){
  tot=length(unique(rarefaction.roa$name[1:k]))
  total.species[k] <- tot
}

ggplot()+
  geom_point(aes(x=1:nrow(rarefaction.roa),y=total.species))+
  ylab("total species richness")+
  xlab("total number of observations")+
  ggtitle("Roatan - all depths")

#### 6. Species list at each site ####
carib.species <- read.csv("2.Ch3.R.analysis/2.file.fish.carib.grouped.csv") %>%
  filter(dband>30 & dband<310) %>%
  select(species)%>%
  unique()

bonaire.species <- read.csv("2.Ch3.R.analysis/2.file.fish.carib.grouped.csv") %>%
  filter(location=="Bonaire") %>%
  select(species) %>%
  unique()%>%
  mutate(bonaire=c(1))

curacao.species <- read.csv("2.Ch3.R.analysis/2.file.fish.carib.grouped.csv") %>%
  filter(location=="Curacao") %>%
  select(species) %>%
  unique()%>%
  mutate(curacao=c(1))

statia.species <- read.csv("2.Ch3.R.analysis/2.file.fish.carib.grouped.csv") %>%
  filter(location=="St. Eustatius") %>%
  select(species) %>%
  unique()%>%
  mutate(statia=c(1))

roatan.species <- read.csv("2.Ch3.R.analysis/2.file.fish.carib.grouped.csv") %>%
  filter(location=="Roatan" & dband>30 & dband<310) %>%
  select(species) %>%
  unique()%>%
  mutate(roatan=c(1))

## presence/absence of species at each site 
carib.species.list <- carib.species %>%
  left_join(bonaire.species) %>%
  left_join(curacao.species) %>%
  left_join(statia.species) %>%
  left_join(roatan.species)
#write.csv(carib.species.list,"2.Ch3.R.analysis/2.file.carib.species.list.csv")

## species only observed at one site 
carib.species.list <- read.csv("2.Ch3.R.analysis/2.file.carib.species.list.csv")

endemics <- carib.species.list %>%
  select(-1)%>%
  filter(rowSums(. == 1,na.rm=TRUE) == 1) %>% #selects lines where "1" is only present once. 
  mutate(endemic="yes") %>%
  select(species,endemic)

endemic.depth <- endemics %>%
  right_join(read.csv("2.Ch3.R.analysis/2.file.species.depth.csv")) %>%
  mutate(predominance=case_when(is.na(predominance)~rating,TRUE~predominance),
         endemic=case_when(is.na(endemic)~"no",TRUE~endemic)) %>%
  select(-X) 

#write.csv(endemic.depth,"endemic.depth.csv")

endemic.category <- endemic.depth %>%
  group_by(predominance,endemic) %>%
  summarize(total=n())
 

