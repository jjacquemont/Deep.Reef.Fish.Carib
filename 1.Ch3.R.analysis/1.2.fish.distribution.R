library(tidyverse)
library(ggplot2)
library(fishualize)
library(viridis)
library(ggpubr)

palette <- viridis(5,begin=0.9,end=0)

#### 1. Depth attributions ####
#### 1.1 Curacao ####
carib.deep.fishbase <- read.csv("2.file.fish.carib.grouped.csv")
carib.species <- read.csv("2.carib.deep.traits.species.csv")

cur.com <- carib.deep.fishbase %>%
  filter(location == "Curacao") %>%
  select(species, dband, location, abundance, abu.corr) %>%
  group_by(species, dband) %>%
  summarize(abu.corr = sum(abu.corr)) %>%
  inner_join(carib.species)

species.distribution.cur <- data.frame(species=c(),u.meso=c(),l.meso=c(),
                                       u.rari=c(),l.rari=c(),tot.abu=c())
for (k in unique(cur.com$species)){
  species <- cur.com[cur.com$species==k,]
  u.meso <- sum(species$abu.corr[species$dband<80])/sum(species$abu.corr)
  l.meso <- sum(species$abu.corr[species$dband<130 & species$dband>=80])/sum(species$abu.corr)
  u.rari <- sum(species$abu.corr[species$dband<190 & species$dband>=130])/sum(species$abu.corr)
  l.rari <- sum(species$abu.corr[species$dband>=190])/sum(species$abu.corr)
  tot.abu <- sum(species$abu.corr)
  species.distribution.cur <- rbind(species.distribution.cur,
                                    c(k,u.meso,l.meso,u.rari,l.rari,tot.abu))
}

colnames(species.distribution.cur) <- c("species","u.meso","l.meso","u.rari","l.rari","tot.abu")
species.distribution.cur <- species.distribution.cur %>%
  mutate(u.meso=as.numeric(u.meso),l.meso=as.numeric(l.meso),u.rari=as.numeric(u.rari),
         l.rari=as.numeric(l.rari),tot.abu=as.numeric(tot.abu)) %>%
  #add number so that species can be ordered by depth predominance
  arrange(u.meso,desc(l.rari),l.meso,u.rari) %>%
  mutate(species=factor(species,levels=species)) %>%
  #stacking columns for ggplot
  gather(key="depthzone", value="abundance",u.meso,l.meso,u.rari,l.rari,tot.abu)  %>%
  mutate(location="Curacao",
    depthzone=factor(depthzone,levels=c("u.meso","l.meso","u.rari","l.rari","tot.abu"))) 

ggplot(species.distribution.cur, aes(y=species,x=depthzone))+
  geom_point(aes(color=abundance))+
  scale_color_viridis(begin=1, end=0)+
  theme_classic()

rel.abundance <- species.distribution.cur %>%
  filter(depthzone!="tot.abu")

ggplot(rel.abundance, aes(y=species,x=depthzone))+
  geom_tile(aes(fill=abundance))+
  xlab("relative abundance per depth zone")+
  scale_x_discrete(position="top")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  ylab("")

#### 1.2 Roatan ####
## Roatan comes second because we have info on altiphotic distribution that
# is used to inform Statia and Bonaire
carib.deep.fishbase <- read.csv("2.file.fish.carib.grouped.csv")

#### 1.2.1. proportion of observation per species per depth zone ####
roa.com <- carib.deep.fishbase %>%
  filter(location == "Roatan")

species.distribution.roa <- data.frame(species=c(),alti= c(),u.meso=c(),l.meso=c(),
                                       u.rari=c(),l.rari=c(),tot.abu=c())
for (k in unique(roa.com$species)){
  species <- roa.com[roa.com$species==k,]
  alti <- sum(species$abu.corr[species$dband<40])/sum(species$abu.corr)
  u.meso <- sum(species$abu.corr[species$dband<100 & species$dband>=40])/sum(species$abu.corr)
  l.meso <- sum(species$abu.corr[species$dband<150 & species$dband>=100])/sum(species$abu.corr)
  u.rari <- sum(species$abu.corr[species$dband<210 & species$dband>=150])/sum(species$abu.corr)
  l.rari <- sum(species$abu.corr[species$dband>=210])/sum(species$abu.corr)
  tot.abu <- sum(species$abu.corr)
  species.distribution.roa <- rbind(species.distribution.roa,
                                    c(k,alti,u.meso,l.meso,u.rari,l.rari,tot.abu))
}

colnames(species.distribution.roa) <- c("species","alti","u.meso","l.meso","u.rari","l.rari","tot.abu")
species.distribution.roa <- species.distribution.roa %>%
  mutate(alti=as.numeric(alti),u.meso=as.numeric(u.meso),l.meso=as.numeric(l.meso),u.rari=as.numeric(u.rari),
         l.rari=as.numeric(l.rari),tot.abu=as.numeric(tot.abu)) %>%
  #add number so that species can be ordered by depth predominance
  arrange(alti,u.meso,desc(l.rari),l.meso,u.rari) %>%
  mutate(species=factor(species,levels=species)) %>%
  #stacking columns for ggplot
  gather(key="depthzone", value="abundance",alti, u.meso,l.meso,u.rari,l.rari,tot.abu)  %>%
  mutate(location="Roatan",
    depthzone=factor(depthzone,levels=c("alti","u.meso","l.meso","u.rari","l.rari","tot.abu"))) 

#ggplot(species.distribution.roa, aes(y=species,x=depthzone))+
 # geom_point(aes(color=abundance))+
  #scale_color_viridis(begin=1, end=0)+
  #theme_classic()

rel.abundance <- species.distribution.roa %>%
  filter(depthzone!="tot.abu")

ggplot(rel.abundance, aes(y=species,x=depthzone))+
  geom_tile(aes(fill=abundance))+
  xlab("relative abundance per depth zone")+
  scale_x_discrete(position="top")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),axis.text = element_text(size = 5))+
  ylab("")

## same but for deep only
roa.com.deep <- carib.deep.fishbase %>%
  filter(location == "Roatan",dband>30) 

species.distribution.roa.deep <- data.frame(species=c(),u.meso=c(),l.meso=c(),
                                       u.rari=c(),l.rari=c(),tot.abu=c())
for (k in unique(roa.com.deep$species)){
  species <- roa.com[roa.com$species==k,]
  u.meso <- sum(species$abu.corr[species$dband<100])/sum(species$abu.corr)
  l.meso <- sum(species$abu.corr[species$dband<150 & species$dband>=100])/sum(species$abu.corr)
  u.rari <- sum(species$abu.corr[species$dband<210 & species$dband>=150])/sum(species$abu.corr)
  l.rari <- sum(species$abu.corr[species$dband>=210])/sum(species$abu.corr)
  tot.abu <- sum(species$abu.corr)
  species.distribution.roa.deep <- rbind(species.distribution.roa.deep,
                                    c(k,u.meso,l.meso,u.rari,l.rari,tot.abu))
}

colnames(species.distribution.roa.deep) <- c("species","u.meso","l.meso","u.rari","l.rari","tot.abu")
species.distribution.roa.deep <- species.distribution.roa.deep %>%
  mutate(u.meso=as.numeric(u.meso),l.meso=as.numeric(l.meso),u.rari=as.numeric(u.rari),
         l.rari=as.numeric(l.rari),tot.abu=as.numeric(tot.abu)) %>%
  #add number so that species can be ordered by depth predominance
  arrange(u.meso,desc(l.rari),l.meso,u.rari) %>%
  mutate(species=factor(species,levels=species)) %>%
  #stacking columns for ggplot
  gather(key="depthzone", value="abundance", u.meso,l.meso,u.rari,l.rari,tot.abu)  %>%
  mutate(location="Roatan",
         depthzone=factor(depthzone,levels=c("u.meso","l.meso","u.rari","l.rari","tot.abu"))) 


#### 1.2.2. distribution across depths ####
roa.com <- carib.deep.fishbase %>%
  filter(location == "Roatan") %>%
  mutate(cluster=case_when(dband<30 ~"altiphotic",
                           dband<150 ~"mesophotic",
                           TRUE~"rariphotic"))%>%
  select(cluster,dband,species,abu.corr)

combinations <- expand.grid(
  species = unique(roa.com$species),
  cluster = unique(roa.com$cluster))

roa.cluster <- combinations %>%
  left_join(roa.com)%>%
  mutate(abu.corr = coalesce(abu.corr, 0)) %>%  # Replace NAs with 0
  group_by(species, cluster) %>%
  summarize(abu.corr = sum(abu.corr))

roa.fish <- unique(roa.cluster$species)
fish.distrib <- vector(length=length(roa.fish))

for (k in 1:length(roa.fish)){
  fish <- roa.fish[k]
  distrib <- roa.cluster[roa.cluster$species==fish,]
  tot.abu <- sum(distrib$abu.corr)
  if (distrib[distrib$cluster=="altiphotic","abu.corr"]/tot.abu>0.2){
    fish.distrib[k] <- "AM"}
  else if (distrib[distrib$cluster=="mesophotic","abu.corr"]/tot.abu>0.75){
    fish.distrib[k] <- "M"}
  else if (distrib[distrib$cluster=="rariphotic","abu.corr"]/tot.abu>0.75){
    fish.distrib[k] <- "R"}
  else{fish.distrib[k] <- "RM"}
  #print(distrib)}
}

roa.fish.distrib <- data_frame(species=roa.fish,distrib.roa=fish.distrib)

## manual change for deep sea species 
roa.fish.distrib$distrib.roa[roa.fish.distrib$species=="Synagrops bellus"] <- "DS"
roa.fish.distrib$distrib.roa[roa.fish.distrib$species=="Synagrops egretta"] <- "DS"
roa.fish.distrib$distrib.roa[roa.fish.distrib$species=="Hollardia hollardi"] <- "DS"

## adding this info to the file
carib.species <- read.csv("2.carib.deep.traits.species.csv") %>%
  select(1:4)%>%
  left_join(roa.fish.distrib)%>%
  mutate(distrib.roa=case_when(distrib.cur=="AM" ~ "AM",
                               distrib.cur=="DS"~ "DS",
                               TRUE~as.character(distrib.roa)))


write_csv(carib.species,"2.carib.deep.traits.species.csv")

#### 1.3. Bonaire ####
carib.deep.fishbase <- read.csv("2.file.fish.carib.grouped.csv")

#### 1.3.1. proportion of observation per species per depth zone ####
bon.com <- carib.deep.fishbase %>%
  filter(location == "Bonaire")

species.distribution.bon <- data.frame(species=c(),u.meso=c(),l.meso=c(),
                                       u.rari=c(),l.rari=c(),tot.abu=c())
for (k in unique(bon.com$species)){
  species <- bon.com[bon.com$species==k,]
  u.meso <- sum(species$abu.corr[species$dband<70])/sum(species$abu.corr)
  l.meso <- sum(species$abu.corr[species$dband<120 & species$dband>=70])/sum(species$abu.corr)
  u.rari <- sum(species$abu.corr[species$dband<170 & species$dband>=120])/sum(species$abu.corr)
  l.rari <- sum(species$abu.corr[species$dband>=170])/sum(species$abu.corr)
  tot.abu <- sum(species$abu.corr)
  species.distribution.bon <- rbind(species.distribution.bon,
                                    c(k,u.meso,l.meso,u.rari,l.rari,tot.abu))
}

colnames(species.distribution.bon) <- c("species","u.meso","l.meso","u.rari","l.rari","tot.abu")
species.distribution.bon <- species.distribution.bon %>%
  mutate(u.meso=as.numeric(u.meso),l.meso=as.numeric(l.meso),u.rari=as.numeric(u.rari),
         l.rari=as.numeric(l.rari),tot.abu=as.numeric(tot.abu)) %>%
  #add number so that species can be ordered by depth predominance
  arrange(u.meso,desc(l.rari),l.meso,u.rari) %>%
  mutate(species=factor(species,levels=species)) %>%
  #stacking columns for ggplot
  gather(key="depthzone", value="abundance",u.meso,l.meso,u.rari,l.rari,tot.abu)  %>%
  mutate(location="Bonaire",
    depthzone=factor(depthzone,levels=c("u.meso","l.meso","u.rari","l.rari","tot.abu"))) 

ggplot(species.distribution.bon, aes(y=species,x=depthzone))+
  geom_point(aes(color=abundance))+
  scale_color_viridis(begin=1, end=0)+
  theme_classic()

rel.abundance <- species.distribution.bon %>%
  filter(depthzone!="tot.abu")

ggplot(rel.abundance, aes(y=species,x=depthzone))+
  geom_tile(aes(fill=abundance))+
  xlab("relative abundance per depth zone")+
  scale_x_discrete(position="top")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),axis.text = element_text(size = 6))+
  ylab("")

#### 1.3.2. distribution across depths ####
bon.com <- carib.deep.fishbase %>%
  filter(location == "Bonaire") %>%
  mutate(cluster=case_when(dband<120 ~"mesophotic",
                           TRUE~"rariphotic"))%>%
  select(cluster,dband,species,abu.corr)

combinations <- expand.grid(
  species = unique(bon.com$species),
  cluster = unique(bon.com$cluster))

bon.cluster <- combinations %>%
  left_join(bon.com)%>%
  mutate(abu.corr = coalesce(abu.corr, 0)) %>%  # Replace NAs with 0
  group_by(species, cluster) %>%
  summarize(abu.corr = sum(abu.corr))

bon.fish <- unique(bon.cluster$species)
fish.distrib <- vector(length=length(bon.fish))

for (k in 1:length(bon.fish)){
  fish <- bon.fish[k]
  distrib <- bon.cluster[bon.cluster$species==fish,]
  tot.abu <- sum(distrib$abu.corr)
  if (distrib[distrib$cluster=="mesophotic","abu.corr"]/tot.abu>0.75){
    print("M")
    fish.distrib[k] <- "M"}
  else if (distrib[distrib$cluster=="rariphotic","abu.corr"]/tot.abu>0.75){
    fish.distrib[k] <- "R"}
  else{print("RM")
    fish.distrib[k] <- "RM"}
}

bon.fish.distrib <- data_frame(species=bon.fish,distrib.bon=fish.distrib)


## adding this info to the file
carib.species <- read.csv("2.carib.deep.traits.species.csv") %>%
  select(1:5)%>%
  left_join(bon.fish.distrib)%>%
  mutate(distrib.bon=case_when(distrib.cur=="AM" | distrib.roa=="AM"~ "AM",
                               distrib.cur=="DS"~ "DS",
                               TRUE~as.character(distrib.bon)))


write_csv(carib.species,"2.carib.deep.traits.species.csv")

#### 1.4. Statia ####
#### 1.4.1. proportion of observation per species per depth zone ####
sta.com <- read.csv("2.file.fish.carib.grouped.csv") %>%
  filter(location == "St. Eustatius")

species.distribution.sta <- data.frame(species=c(),u.meso=c(),l.meso=c(),
                                       u.rari=c(),l.rari=c(),tot.abu=c())
for (k in unique(sta.com$species)){
  species <- sta.com[sta.com$species==k,]
  u.meso <- sum(species$abu.corr[species$dband<100])/sum(species$abu.corr)
  l.meso <- sum(species$abu.corr[species$dband<140 & species$dband>=100])/sum(species$abu.corr)
  u.rari <- sum(species$abu.corr[species$dband<190 & species$dband>=140])/sum(species$abu.corr)
  l.rari <- sum(species$abu.corr[species$dband>=190])/sum(species$abu.corr)
  tot.abu <- sum(species$abu.corr)
  species.distribution.sta <- rbind(species.distribution.sta,
                                    c(k,u.meso,l.meso,u.rari,l.rari,tot.abu))
}

colnames(species.distribution.sta) <- c("species","u.meso","l.meso","u.rari","l.rari","tot.abu")
species.distribution.sta <- species.distribution.sta %>%
  mutate(u.meso=as.numeric(u.meso),l.meso=as.numeric(l.meso),u.rari=as.numeric(u.rari),
         l.rari=as.numeric(l.rari),tot.abu=as.numeric(tot.abu)) %>%
  #add number so that species can be ordered by depth predominance
  arrange(u.meso,desc(l.rari),l.meso,u.rari) %>%
  mutate(species=factor(species,levels=species)) %>%
  #stacking columns for ggplot
  gather(key="depthzone", value="abundance",u.meso,l.meso,u.rari,l.rari,tot.abu)  %>%
  mutate(location="Statia",
    depthzone=factor(depthzone,levels=c("u.meso","l.meso","u.rari","l.rari","tot.abu"))) 

ggplot(species.distribution.sta, aes(y=species,x=depthzone))+
  geom_point(aes(color=abundance))+
  scale_color_viridis(begin=1, end=0)+
  theme_classic()

rel.abundance <- species.distribution.sta %>%
  filter(depthzone!="tot.abu")

ggplot(rel.abundance, aes(y=species,x=depthzone))+
  geom_tile(aes(fill=abundance))+
  xlab("relative abundance per depth zone")+
  scale_x_discrete(position="top")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),axis.text = element_text(size = 6))+
  ylab("")

#### 1.4.2. distribution across depths ####
sta.com <- read.csv("2.file.fish.carib.grouped.csv") %>%
  filter(location == "St. Eustatius") %>%
  mutate(cluster=case_when(dband<140 ~"mesophotic",
                           TRUE~"rariphotic"))%>%
  select(cluster,dband,species,abu.corr)

combinations <- expand.grid(
  species = unique(sta.com$species),
  cluster = unique(sta.com$cluster))

sta.cluster <- combinations %>%
  left_join(sta.com)%>%
  mutate(abu.corr = coalesce(abu.corr, 0)) %>%  # Replace NAs with 0
  group_by(species, cluster) %>%
  summarize(abu.corr = sum(abu.corr))

sta.fish <- unique(sta.cluster$species)
fish.distrib <- vector(length=length(sta.fish))

for (k in 1:length(sta.fish)){
  fish <- sta.fish[k]
  distrib <- sta.cluster[sta.cluster$species==fish,]
  tot.abu <- sum(distrib$abu.corr)
  if (distrib[distrib$cluster=="mesophotic","abu.corr"]/tot.abu>0.75){
    print("M")
    fish.distrib[k] <- "M"}
  else if (distrib[distrib$cluster=="rariphotic","abu.corr"]/tot.abu>0.75){
    fish.distrib[k] <- "R"}
  else{print("RM")
    fish.distrib[k] <- "RM"}
}

sta.fish.distrib <- data_frame(species=sta.fish,distrib.sta=fish.distrib)

## adding this info to the file
carib.species <- read.csv("2.carib.deep.traits.species.csv") %>%
  select(1:6)%>%
  left_join(sta.fish.distrib)%>%
  mutate(distrib.sta=case_when(distrib.cur=="AM" | distrib.roa=="AM"~ "AM",
                               distrib.cur=="DS"~ "DS",
                               TRUE~as.character(distrib.sta)))

## adding expert knowledge 
LT.rating <- read.csv("LT_rating.csv") 

carib.species <- carib.species %>%
  left_join(LT.rating) %>%
  mutate(across(starts_with("distrib"), ~case_when(!is.na(.) & !is.na(Luke.rating)~Luke.rating,
                                                   TRUE~.)))
         
#write_csv(depth.rating,"2.carib.deep.traits.species.csv")

#### 1.5. Depth rating across 4 sites ####
species.distribution <- rbind(species.distribution.cur,species.distribution.bon,
                              species.distribution.sta,species.distribution.roa.deep)

species.distribution.overall <- data.frame(species=c(),u.meso=c(),l.meso=c(),
                                       u.rari=c(),l.rari=c(),tot.abu=c())
for (k in unique(species.distribution$species)){
  species <- species.distribution[species.distribution$species==k,]
  #alti <- mean(species$abundance[species$depthzone=="alti"])
  u.meso <- mean(species$abundance[species$depthzone=="u.meso"])
  l.meso <- mean(species$abundance[species$depthzone=="l.meso"])
  u.rari <- mean(species$abundance[species$depthzone=="u.rari"])
  l.rari <- mean(species$abundance[species$depthzone=="l.rari"])
  tot.abu <- sum(species$abundance)
  species.distribution.overall <- rbind(species.distribution.overall,
                                    c(k,u.meso,l.meso,u.rari,l.rari,tot.abu))
}

colnames(species.distribution.overall) <- c("species","u.meso","l.meso","u.rari","l.rari","tot.abu")
species.distribution.overall <- species.distribution.overall %>%
  mutate(u.meso=as.numeric(u.meso),l.meso=as.numeric(l.meso),u.rari=as.numeric(u.rari),
         l.rari=as.numeric(l.rari),tot.abu=as.numeric(tot.abu)) %>%
  #add number so that species can be ordered by depth predominance
  arrange(u.meso,l.meso,desc(l.rari),u.rari) %>%
  mutate(species=factor(species,levels=species)) %>%
  #overallcking columns for ggplot
  gather(key="depthzone", value="abundance",u.meso,l.meso,u.rari,l.rari,tot.abu)  %>%
  mutate(location="overall",
         depthzone=factor(depthzone,levels=c("u.meso","l.meso","u.rari","l.rari","tot.abu"))) 

rel.abundance <- species.distribution.overall %>%
  filter(depthzone!="tot.abu") %>%
  arrange(species)
rel.abundance1 <- rel.abundance[1:nrow(rel.abundance)/2,]
rel.abundance2 <- rel.abundance[(nrow(rel.abundance)/2):nrow(rel.abundance),]

ggplot(rel.abundance1, aes(y=species,x=depthzone))+
  geom_tile(aes(fill=abundance))+
  xlab("relative abundance per depth zone")+
  scale_x_discrete(position="top",labels=str_wrap(c("upper mesophotic","lower mesophotic",
                                           "upper rariphotic","lower rariphotic"),7))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),axis.text = element_text(size = 6))+
  ylab("")+
  coord_cartesian(clip = "off")+
 # annotate(geom="text",y=1:length(unique(species.distribution.overall$species)),
   #        x=4.6,label=round(species.distribution.overall$abundance[species.distribution.overall$depthzone=="tot.abu"],0),
  #         size=2)
ggplot(rel.abundance2, aes(y=species,x=depthzone))+
  geom_tile(aes(fill=abundance))+
  xlab("relative abundance per depth zone")+
  scale_x_discrete(position="top",labels=str_wrap(c("upper mesophotic","lower mesophotic",
                                                    "upper rariphotic","lower rariphotic"),7))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),axis.text = element_text(size = 6))+
  ylab("")+
  coord_cartesian(clip = "off")
  
#### attributing a depth predominance to each species based on their distribution
# across the four sites, and assessing the contribution of depth predominance to
# the diversity and abundance across depth at each site 
rel.abundance.spread <- rel.abundance %>%
  spread(depthzone,abundance) %>%
  mutate(predominance=case_when(
    u.meso+l.meso>0.75 ~ "M",
    u.rari+l.rari>0.75 ~ "R",
    TRUE ~ "MR")) %>%
  left_join(read_csv("expert.ratings.csv")) %>%
  mutate(predominance=case_when(
   !is.na(rating)~rating,
    TRUE~predominance)) %>%
  select(predominance,species)  %>%
  left_join(carib.deep.fishbase) %>%
  group_by(predominance,location,dband) %>%
  summarize(abu.corr=sum(abu.corr),spric=n()) %>%
  mutate(predominance=factor(predominance, levels=c("AM","M","MR","R","DS")),
         location=factor(location, levels=c("Curacao","Bonaire","St. Eustatius", "Roatan")))

# to obtain dashed line for site specific community breaks 
dashed_lines <- data_frame(location=rep(c("Curacao","Bonaire","St. Eustatius", "Roatan"),each=3),
                           breaks=c(70,120,180,60,110,160,90,130,180,90,140,180)) %>%
  mutate(location=factor(location, levels=c("Curacao","Bonaire","St. Eustatius", "Roatan")))

## plot for species richness
blue_palette <- c("#BAE0F3","#87CEEB","#6CBDE9","#50ABE7", "#4895EF", "#4361EE", "#2835AF", "#12086F")
mesoP_palette <- c("#DCEEF3","#C2E2EA","#A7D5E1","#8DC8D8","#72BBCE")
rariP_palette <- c("#DDC5EB","#BC90DB","#831EB6")
pal <- c("AM"="#d0dc36","M"="#50ABE7","MR"="#BC90DB","R"="#831EB6","DS"="#12086F")
ggplot(data=rel.abundance.spread,aes(x=dband,y=spric,color=predominance))+
  geom_point(aes(shape = predominance),size=3) +
  scale_shape_manual(values=c(18,18,18,18,12))+
  xlab("depth")+
  ylab("species richness")+
  scale_color_manual(values=pal,
                     labels=c("alti/mesophotic","mesophotic","meso/rariphotic",
                              "rariphotic","deep sea"),name="predominant depth")+
  guides(shape = "none")+
  xlim(300,40)+
  theme_classic()+
  geom_vline(data = dashed_lines, aes(xintercept = breaks), linetype = "dashed",
             color = "grey69",size=0.3)+
  coord_flip()+
  facet_wrap(~location, ncol=4)

## bar plots
rel.abundance.realms <- rel.abundance %>%
  spread(depthzone,abundance) %>%
  mutate(predominance=case_when(
    u.meso+l.meso>0.75 ~ "M",
    u.rari+l.rari>0.75 ~ "R",
    TRUE ~ "MR")) %>%
  left_join(read_csv("expert.ratings.csv")) %>%
  mutate(predominance=case_when(
    !is.na(rating)~rating,
    TRUE~predominance)) %>%
  select(predominance,species)  %>%
  left_join(carib.deep.fishbase) %>%
  group_by(species,location,dband) %>%
  mutate(depth.realm=case_when(
    location=="Curacao" & dband<=70 ~ "upper mesophotic",
    location=="Curacao" & dband<=120 ~ "lower mesophotic",
    location=="Curacao" & dband<=180 ~ "upper rariphotic",
    location=="Curacao" & dband >180 ~ "lower rariphotic",
    location=="Bonaire" & dband<=60 ~ "upper mesophotic",
    location=="Bonaire" & dband<=110 ~ "lower mesophotic",
    location=="Bonaire" & dband<=160 ~ "upper rariphotic",
    location=="Bonaire" & dband >160 ~ "lower rariphotic",
    location=="St. Eustatius" & dband<=90 ~ "upper mesophotic",
    location=="St. Eustatius" & dband<=130 ~ "lower mesophotic",
    location=="St. Eustatius" & dband<=180 ~ "upper rariphotic",
    location=="St. Eustatius" & dband >180 ~ "lower rariphotic",
    location=="Roatan" & dband<=90 ~ "upper mesophotic",
    location=="Roatan" & dband<=140 ~ "lower mesophotic",
    location=="Roatan" & dband<=180 ~ "upper rariphotic",
    location=="Roatan" & dband >180 ~ "lower rariphotic")) %>%
  group_by(depth.realm,species,predominance) %>%
  summarize(abu.corr=sum(abu.corr))%>%
  group_by(depth.realm,predominance) %>%
  summarize(spric=n(),abu.corr=sum(abu.corr))%>%
 mutate(predominance=factor(predominance,levels=c("AM","M","MR","R","DS")),
  depth.realm=factor(depth.realm, levels=rev(c("upper mesophotic","lower mesophotic",
                                     "upper rariphotic","lower rariphotic"))))

## bar plot richness
ggplot(data=rel.abundance.realms, aes(x=depth.realm,y=spric,fill=predominance))+
  geom_bar(position="stack",stat="identity",width=0.4) +

## bar plot abundance
ggplot(data=rel.abundance.realms, aes(x=depth.realm,y=abu.corr,fill=predominance))+
  geom_bar(position="stack",stat="identity",width=0.4) +
  scale_fill_manual(values=pal,
                    labels=c("alti/mesophotic","mesophotic","meso/rariphotic",
                             "rariphotic","deep sea"),
                    name=str_wrap("species depth affinity"),7)+
  scale_x_discrete(label=str_wrap(levels(rel.abundance.realms$depth.realm),7))+
  scale_y_continuous(breaks=c(0,50000))+
  theme(axis.line.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  xlab("")+
  ylab("relative abundance")+
  coord_flip()


## plot for abundance
ggplot(data=rel.abundance.spread,aes(x=dband,y=log(abu.corr),color=predominance))+
  geom_point(aes(shape = predominance),size=3) +
  xlab("depth")+
  ylab("species richness")+
  scale_color_viridis_d(labels=c("alti/mesophotic","mesophotic","meso/rariphotic",
                                 "rariphotic","deep sea"),name="predominant depth",
                        begin=0.9,end=0)+
  #scale_color_discrete(labels=c("alti/mesophotic","mesophotic","meso/rariphotic",
  #                    "rariphotic","deep sea"),name="predominant depth")+
  guides(shape = "none")+
  xlim(300,40)+
  theme_classic()+
  geom_vline(xintercept = 80, lty = 2, color = "grey69",size=0.3) +
  geom_vline(xintercept = 130, lty = 2, color = "grey69",size=0.3) +
  geom_vline(xintercept = 190, lty = 2, color = "grey69",size=0.3) +
  coord_flip()+
  #scale_x_reverse()+
  facet_wrap(~location, ncol=4)


#### 2. Contribution of depth predominance across depth ####
#### 2.1 Curacao ####
abu.com <- cur.com %>%
  select(distrib.cur, dband, abu.corr) %>%
  group_by(distrib.cur, dband) %>%
  summarize(abu.corr = sum(abu.corr))

abu.com <- abu.com %>%
  mutate(distrib.cur=factor(distrib.cur,levels=c("AM","M","RM","R","DS")))

blue_palette <- c("#BAE0F3","#87CEEB","#6CBDE9","#50ABE7", "#4895EF", "#4361EE", "#2835AF", "#12086F")
mesoP_palette <- c("#DCEEF3","#C2E2EA","#A7D5E1","#8DC8D8","#72BBCE")
rariP_palette <- c("#DDC5EB","#BC90DB","#831EB6")
final_palette <- c("#72BBCE","#4895EF","#DDC5EB","#BC90DB","#12086F")

ggplot(data=abu.com,aes(x=dband,y=abu.corr,color=distrib.cur))+
  geom_point(aes(shape = distrib.cur),size=3) +
  xlab("depth")+
  ylab("abundance")+
  ggtitle("Curaçao")+
  scale_color_viridis_d(labels=c("alti/mesophotic","mesophotic","meso/rariphotic",
                                 "rariphotic","deep sea"),name="predominant depth",
                        begin=0.9,end=0)+
  #scale_color_discrete(labels=c("alti/mesophotic","mesophotic","meso/rariphotic",
  #                    "rariphotic","deep sea"),name="predominant depth")+
  guides(shape = "none")+
  theme_classic()+
  geom_vline(xintercept = 80, lty = 2, color = "grey69",size=0.3) +
  geom_vline(xintercept = 130, lty = 2, color = "grey69",size=0.3) +
  geom_vline(xintercept = 190, lty = 2, color = "grey69",size=0.3) +
  coord_flip()+
  scale_x_reverse()

#### 5.1.2.
abu.clus <- cur.com %>%
  mutate(cluster=case_when(dband<=70 ~ "upper mesophotic",
                           dband<=120~ "lower mesophotic",
                           dband<=180~ "upper rariphotic",
                           dband>180~ "lower rariphotic")) %>%
  group_by(cluster,distrib.cur) %>%
  summarize(abu.corr=sum(abu.corr))%>%
  mutate(distrib.cur=factor(distrib.cur,levels=c("AM","M","RM","R","DS")))%>%
  mutate(cluster=factor(cluster,levels=rev(c("upper mesophotic","lower mesophotic","upper rariphotic","lower rariphotic"))))

##plot
ggplot(data=abu.clus, aes(x=cluster,y=abu.corr,fill=distrib.cur))+
  geom_bar(position="stack",stat="identity",width=0.2) +
  scale_fill_viridis_d(labels=c("alti/mesophotic","mesophotic","meso/rariphotic",
                                "rariphotic","deep sea"),name="predominant depth",
                       begin=0.9,end=0)+
  scale_x_discrete(label=str_wrap(levels(abu.clus$cluster),7))+
  scale_y_continuous(breaks=c(0,2000))+
  theme(axis.line.x = element_blank(),
        #axis.text.x = element_blank(),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+
  ylab("total abundance")+
  coord_flip()

#### 2.1.2. Contribution of depth predominance categories to richness
spric.com <- cur.com %>%
  group_by(distrib.cur, dband) %>%
  summarize(spric = n()) %>%
  mutate(distrib.cur=factor(distrib.cur,levels=c("AM","M","RM","R","DS")))

cur.spric.points <- ggplot(data=spric.com,aes(x=dband,y=spric,color=distrib.cur))+
  geom_point(aes(shape = distrib.cur),size=3) +
  xlab("depth")+
  ylab("species richness")+
  ggtitle("Curacao")+
  scale_color_viridis_d(labels=c("alti/mesophotic","mesophotic","meso/rariphotic",
                                 "rariphotic","deep sea"),name="predominant depth",
                        begin=0.9,end=0)+
  #scale_color_discrete(labels=c("alti/mesophotic","mesophotic","meso/rariphotic",
  #                    "rariphotic","deep sea"),name="predominant depth")+
  guides(shape = "none")+
  theme_classic()+
  geom_vline(xintercept = 80, lty = 2, color = "grey69",size=0.3) +
  geom_vline(xintercept = 130, lty = 2, color = "grey69",size=0.3) +
  geom_vline(xintercept = 190, lty = 2, color = "grey69",size=0.3) +
  coord_flip()+
  ggtitle("Curaçao")+
  scale_x_reverse()

#### 5.1.2.
spric.clus <- cur.com %>%
  mutate(cluster=case_when(dband<=70 ~ "upper mesophotic",
                           dband<=120~ "lower mesophotic",
                           dband<=180~ "upper rariphotic",
                           dband>180~ "lower rariphotic")) %>%
  group_by(cluster,species,distrib.cur) %>%
  summarize(abundance=sum(abu.corr))%>%
  group_by(cluster,distrib.cur) %>%
  summarize(spric=n())%>%
  mutate(distrib.cur=factor(distrib.cur,levels=c("AM","M","RM","R","DS")))%>%
  mutate(cluster=factor(cluster,levels=rev(c("upper mesophotic","lower mesophotic","upper rariphotic","lower rariphotic"))))

##plot
cur.spric.barplot <- ggplot(data=spric.clus, aes(x=cluster,y=spric,fill=distrib.cur))+
  geom_bar(position="stack",stat="identity",width=0.4) +
  scale_fill_viridis_d(labels=c("alti/mesophotic","mesophotic","meso/rariphotic",
                                "rariphotic","deep sea"),name="predominant depth",
                       begin=0.9,end=0)+
  scale_x_discrete(label=str_wrap(levels(abu.clus$cluster),7))+
  theme(axis.line.x = element_blank(),
       # axis.text.y=element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  xlab("")+
  ggtitle("Curaçao")+
  ylab("total richness")+
  guides(fill="none")+
  coord_flip()


#### 2.2. Roatan ####
####  abundance
roa.abu.com <- roa.com %>%
  inner_join(carib.species) %>%
  group_by(distrib.roa, dband) %>%
  summarize(abu.corr = sum(abu.corr)) %>%
  mutate(distrib.roa=factor(distrib.roa,levels=c("AM","M","RM","R","DS")))%>%
  filter(dband<310 & dband>30)

palette <- viridis(5,begin=0.9,end=0)

ggplot(data=roa.abu.com,aes(x=dband,y=abu.corr,color=distrib.roa))+
  geom_point(aes(shape = distrib.roa),size=3) +
  xlab("depth")+
  ylab("abundance")+
  ggtitle("Roatan")+
  scale_color_manual(labels=c("alti/mesophotic","mesophotic","meso/rariphotic",
                              "rariphotic","deep sea"),name="predominant depth",
                     values = palette[1:5])+
  #scale_color_discrete(labels=c("alti/mesophotic","mesophotic","meso/rariphotic",
  #                    "rariphotic","deep sea"),name="predominant depth")+
  guides(shape = "none")+
  theme_classic()+
  geom_vline(xintercept = 100, lty = 2, color = "grey69",size=0.3) +
  geom_vline(xintercept = 150, lty = 2, color = "grey69",size=0.3) +
  geom_vline(xintercept = 210, lty = 2, color = "grey69",size=0.3) +
  coord_flip()+
  scale_x_reverse()

#### Barplot
roa.abu.clus <- roa.abu.com %>%
  mutate(cluster=case_when(dband<100 ~ "upper mesophotic",
                           dband<150 ~ "lower mesophotic",
                           dband<210~ "upper rariphotic",
                           dband>=210~ "lower rariphotic")) %>%
  group_by(cluster,distrib.roa) %>%
  summarize(abu.corr=sum(abu.corr))%>%
  mutate(distrib.roa=factor(distrib.roa,levels=c("AM","M","RM","R","DS")))%>%
  mutate(cluster=factor(cluster,levels=rev(c("upper mesophotic","lower mesophotic","upper rariphotic","lower rariphotic"))))

##plot
ggplot(data=roa.abu.clus, aes(x=cluster,y=abu.corr,fill=distrib.roa))+
  geom_bar(position="stack",stat="identity",width=0.3) +
  scale_fill_manual(labels=c("alti/mesophotic","mesophotic","meso/rariphotic",
                             "rariphotic","deep sea"),name="predominant depth",
                    values = palette[1:5])+
  scale_x_discrete(label=str_wrap(levels(roa.abu.clus$cluster),7))+
  scale_y_continuous(breaks=c(0,2000))+
  theme(axis.line.x = element_blank(),
        #axis.text.x = element_blank(),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank())+
  ggtitle("Roatan")+
  xlab("")+
  ylab("total abundance")+
  coord_flip()

#### 2.2.2. richness
roa.spric.com <- roa.com %>%
  inner_join(carib.species) %>%
  group_by(distrib.roa, dband) %>%
  summarize(spric = n())%>%
  mutate(distrib.roa=factor(distrib.roa,levels=c("AM","M","RM","R","DS"))) %>%
  filter(dband<310 & dband>30)


roa.spric.points <- ggplot(data=roa.spric.com,aes(x=dband,y=spric,color=distrib.roa))+
  geom_point(aes(shape = distrib.roa),size=3) +
  xlab("depth")+
  ylab("species richness")+
  ggtitle("Roatan")+
  scale_color_manual(labels=c("alti/mesophotic","mesophotic","meso/rariphotic",
                              "rariphotic","deep sea"),name="predominant depth",
                     values = palette[1:5])+
  guides(shape = "none")+
  theme_classic()+
  geom_vline(xintercept = 40, lty = 2, color = "grey69",size=0.3) +
  geom_vline(xintercept = 100, lty = 2, color = "grey69",size=0.3) +
  geom_vline(xintercept = 150, lty = 2, color = "grey69",size=0.3) +
  geom_vline(xintercept = 210, lty = 2, color = "grey69",size=0.3) +
  coord_flip()+
  ggtitle("Roatan")+
  scale_x_reverse()

#### 5.1.2.
roa.spric.clus <- roa.com %>%
  inner_join(carib.species) %>%
  mutate(cluster=case_when(dband<100 ~ "upper mesophotic",
                           dband<150 ~ "lower mesophotic",
                           dband<210~ "upper rariphotic",
                           dband>=210~ "lower rariphotic")) %>%
  group_by(cluster,species,distrib.roa) %>%
  summarize(abundance=sum(abu.corr))%>%
  group_by(cluster,distrib.roa) %>%
  summarize(spric=n())%>%
  mutate(distrib.roa=factor(distrib.roa,levels=c("AM","M","RM","R","DS")))%>%
  mutate(cluster=factor(cluster,levels=rev(c("upper mesophotic","lower mesophotic","upper rariphotic","lower rariphotic"))))


##plot
roa.spric.barplot <- ggplot(data=roa.spric.clus, aes(x=cluster,y=spric,fill=distrib.roa))+
  geom_bar(position="stack",stat="identity",width=0.4) +
  scale_fill_manual(labels=c("alti/mesophotic","mesophotic","meso/rariphotic",
                             "rariphotic","deep sea"),name=str_wrap("fish predominant depth",6),
                    values = palette[1:5])+
  scale_x_discrete(label=str_wrap(levels(roa.spric.clus$cluster),7))+
  theme(axis.line.x = element_blank(),
        axis.text.y=element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  xlab("")+
  ylab("total richness")+
  ggtitle("Roatan")+
  #guides(fill="none")+
  coord_flip()

#### 2.3. Bonaire ####
#### 3.1 Contribution of depth predominance categories to abundance
bon.abu.com <- bon.com %>%
  inner_join(carib.species) %>%
  group_by(distrib.bon, dband) %>%
  summarize(abu.corr = sum(abu.corr)) %>%
  mutate(distrib.bon=factor(distrib.bon,levels=c("AM","M","RM","R","DS")))



ggplot(data=bon.abu.com,aes(x=dband,y=abu.corr,color=distrib.bon))+
  geom_point(aes(shape = distrib.bon),size=3) +
  xlab("depth")+
  ylab("abundance")+
  ggtitle("Bonaire")+
  scale_color_manual(labels=c("alti/mesophotic","mesophotic","meso/rariphotic",
                                 "rariphotic","deep sea"),name="predominant depth",
                        values = palette[1:5])+
  #scale_color_discrete(labels=c("alti/mesophotic","mesophotic","meso/rariphotic",
  #                    "rariphotic","deep sea"),name="predominant depth")+
  guides(shape = "none")+
  theme_classic()+
  geom_vline(xintercept = 70, lty = 2, color = "grey69",size=0.3) +
  geom_vline(xintercept = 120, lty = 2, color = "grey69",size=0.3) +
  geom_vline(xintercept = 170, lty = 2, color = "grey69",size=0.3) +
  coord_flip()+
  scale_x_reverse()

#### Barplot
bon.abu.clus <- bon.abu.com %>%
  mutate(cluster=case_when(dband<70 ~ "upper mesophotic",
                           dband<120~ "lower mesophotic",
                           dband<170~ "upper rariphotic",
                           dband>=170~ "lower rariphotic")) %>%
  group_by(cluster,distrib.bon) %>%
  summarize(abu.corr=sum(abu.corr))%>%
  mutate(distrib.bon=factor(distrib.bon,levels=c("AM","M","RM","R","DS")))%>%
  mutate(cluster=factor(cluster,levels=rev(c("upper mesophotic","lower mesophotic","upper rariphotic","lower rariphotic"))))

##plot
ggplot(data=bon.abu.clus, aes(x=cluster,y=abu.corr,fill=distrib.bon))+
  geom_bar(position="stack",stat="identity",width=0.3) +
  scale_fill_manual(labels=c("alti/mesophotic","mesophotic","meso/rariphotic",
                              "rariphotic","deep sea"),name="predominant depth",
                     values = palette[1:5])+
  scale_x_discrete(label=str_wrap(levels(bon.abu.clus$cluster),7))+
  scale_y_continuous(breaks=c(0,2000))+
  theme(axis.line.x = element_blank(),
        #axis.text.x = element_blank(),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+
  ylab("total abundance")+
  coord_flip()

#### 1.3. Contribution of depth predominance categories to richness
bon.spric.com <- bon.com %>%
  inner_join(carib.species) %>%
  group_by(distrib.bon, dband) %>%
  summarize(spric = n())%>%
  mutate(distrib.bon=factor(distrib.bon,levels=c("AM","M","RM","R","DS"))) 


bon.spric.points <- ggplot(data=bon.spric.com,aes(x=dband,y=spric,color=distrib.bon))+
  geom_point(aes(shape = distrib.bon),size=3) +
  xlab("depth")+
  ylab("species richness")+
  ggtitle("Bonaire")+
  scale_color_manual(labels=c("alti/mesophotic","mesophotic","meso/rariphotic",
                             "rariphotic","deep sea"),name="predominant depth",
                    values = palette[1:5])+
  guides(shape = "none")+
  theme_classic()+
  geom_vline(xintercept = 70, lty = 2, color = "grey69",size=0.3) +
  geom_vline(xintercept = 120, lty = 2, color = "grey69",size=0.3) +
  geom_vline(xintercept = 170, lty = 2, color = "grey69",size=0.3) +
  coord_flip()+
  ggtitle("Bonaire")+
  scale_x_reverse()

#### 5.1.2.
bon.spric.clus <- bon.com %>%
  inner_join(carib.species) %>%
  mutate(cluster=case_when(dband<70 ~ "upper mesophotic",
                           dband<120~ "lower mesophotic",
                           dband<170~ "upper rariphotic",
                           dband>=170~ "lower rariphotic")) %>%
  group_by(cluster,species,distrib.bon) %>%
  summarize(abundance=sum(abu.corr))%>%
  group_by(cluster,distrib.bon) %>%
  summarize(spric=n())%>%
  mutate(distrib.bon=factor(distrib.bon,levels=c("AM","M","RM","R","DS")))%>%
  mutate(cluster=factor(cluster,levels=rev(c("upper mesophotic","lower mesophotic","upper rariphotic","lower rariphotic"))))


##plot
bon.spric.barplot <- ggplot(data=bon.spric.clus, aes(x=cluster,y=spric,fill=distrib.bon))+
  geom_bar(position="stack",stat="identity",width=0.4) +
  scale_fill_manual(labels=c("alti/mesophotic","mesophotic","meso/rariphotic",
                              "rariphotic","deep sea"),name="predominant depth",
                     values = palette[1:5])+
  scale_x_discrete(label=str_wrap(levels(bon.abu.clus$cluster),7))+
  theme(axis.line.x = element_blank(),
        axis.text.y=element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  xlab("")+
  ggtitle("Bonaire")+
  ylab("total richness")+
  guides(fill="none")+
  coord_flip()


#### 2.4. Statia ####
#### 4.1. Contribution of depth predominance categories to abundance
sta.abu.com <- sta.com %>%
  inner_join(carib.species) %>%
  group_by(distrib.sta, dband) %>%
  summarize(abu.corr = sum(abu.corr)) %>%
  mutate(distrib.sta=factor(distrib.sta,levels=c("AM","M","RM","R","DS")))

palette <- viridis(5,begin=0.9,end=0)

ggplot(data=sta.abu.com,aes(x=dband,y=abu.corr,color=distrib.sta))+
  geom_point(aes(shape = distrib.sta),size=3) +
  xlab("depth")+
  ylab("abundance")+
  ggtitle("Statia")+
  scale_color_manual(labels=c("alti/mesophotic","mesophotic","meso/rariphotic",
                              "rariphotic","deep sea"),name="predominant depth",
                     values = palette[1:4])+
  #scale_color_discrete(labels=c("alti/mesophotic","mesophotic","meso/rariphotic",
  #                    "rariphotic","deep sea"),name="predominant depth")+
  guides(shape = "none")+
  theme_classic()+
  geom_vline(xintercept = 100, lty = 2, color = "grey69",size=0.3) +
  geom_vline(xintercept = 190, lty = 2, color = "grey69",size=0.3) +
  #geom_vline(xintercept = 170, lty = 2, color = "grey69",size=0.3) +
  coord_flip()+
  scale_x_reverse()

#### Barplot
sta.abu.clus <- sta.abu.com %>%
  mutate(cluster=case_when(dband<100 ~ "mesophotic",
                           dband<190~ "upper rariphotic",
                           dband>=190~ "lower rariphotic")) %>%
  group_by(cluster,distrib.sta) %>%
  summarize(abu.corr=sum(abu.corr))%>%
  mutate(distrib.sta=factor(distrib.sta,levels=c("AM","M","RM","R","DS")))%>%
  mutate(cluster=factor(cluster,levels=rev(c("mesophotic","upper rariphotic","lower rariphotic"))))

##plot
ggplot(data=sta.abu.clus, aes(x=cluster,y=abu.corr,fill=distrib.sta))+
  geom_bar(position = "stack", stat="identity",width=0.3) +
  scale_fill_manual(labels=c("alti/mesophotic","mesophotic","meso/rariphotic",
                             "rariphotic","deep sea"),name="predominant depth",
                    values = palette[1:4])+
  scale_x_discrete(label=str_wrap(levels(sta.abu.clus$cluster),7))+
  scale_y_continuous(breaks=c(0,2000))+
  theme(axis.line.x = element_blank(),
        #axis.text.x = element_blank(),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+
  ylab("total abundance")+
  coord_flip()

#### 2.3. Contribution of depth predominance categories to richness
sta.spric.com <- sta.com %>%
  inner_join(carib.species) %>%
  group_by(distrib.sta, dband) %>%
  summarize(spric = n())%>%
  mutate(distrib.sta=factor(distrib.sta,levels=c("AM","M","RM","R","DS"))) 


sta.spric.points <- ggplot(data=sta.spric.com,aes(x=dband,y=spric,color=distrib.sta))+
  geom_point(aes(shape = distrib.sta),size=3) +
  xlab("depth")+
  ylab("species richness")+
  ggtitle("Statia")+
  scale_color_manual(labels=c("alti/mesophotic","mesophotic","meso/rariphotic",
                              "rariphotic","deep sea"),name="predominant depth",
                     values = palette[1:4])+
  guides(shape = "none")+
  theme_classic()+
  geom_vline(xintercept = 100, lty = 2, color = "grey69",size=0.3) +
  geom_vline(xintercept = 140, lty = 2, color = "grey69",size=0.3) +
  geom_vline(xintercept = 190, lty = 2, color = "grey69",size=0.3) +
  coord_flip()+
  ggtitle("Statia")+
  scale_x_reverse()

#### 5.1.2.
sta.spric.clus <- sta.com %>%
  inner_join(carib.species) %>%
  mutate(cluster=case_when(dband<100 ~ "upper mesophotic",
                           dband<140~ "lower mesophotic",
                           dband<190~ "upper rariphotic",
                           dband>=190~ "lower rariphotic")) %>%
  group_by(cluster,species,distrib.sta) %>%
  summarize(abundance=sum(abu.corr))%>%
  group_by(cluster,distrib.sta) %>%
  summarize(spric=n())%>%
  mutate(distrib.sta=factor(distrib.sta,levels=c("AM","M","RM","R","DS")))%>%
  mutate(cluster=factor(cluster,levels=rev(c("upper mesophotic","lower mesophotic","upper rariphotic","lower rariphotic"))))


##plot
sta.spric.barplot <- ggplot(data=sta.spric.clus, aes(x=cluster,y=spric,fill=distrib.sta))+
  geom_bar(position="stack",stat="identity",width=0.4) +
  scale_fill_manual(labels=c("alti/mesophotic","mesophotic","meso/rariphotic",
                             "rariphotic","deep sea"),name="predominant depth",
                    values = palette[1:4])+
  scale_x_discrete(label=str_wrap(levels(sta.spric.clus$cluster),7))+
  theme(axis.line.x = element_blank(),
        axis.text.y=element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  xlab("")+
  ggtitle("Statia")+
  ylab("total richness")+
  guides(fill="none")+
  coord_flip()

#### 2.5. combined plot ####
library(ggpubr)
ggarrange(cur.spric.barplot, bon.spric.barplot,
                                        sta.spric.barplot, roa.spric.barplot, ncol = 4,
                                        common.legend=TRUE,legend="bottom")

ggarrange(cur.spric.points, bon.spric.points + rremove("ylab"),
          sta.spric.points+ rremove("ylab"), roa.spric.points+ rremove("ylab"), ncol = 4,
          common.legend=TRUE,legend="bottom")

#### 3. Top abundant species in each cluster  ####
carib.deep.fishbase <- read.csv("2.file.fish.carib.grouped.csv")

## Curacao
abu.sp <- carib.deep.fishbase %>%
  filter(location=="Curacao",dband>30,dband<310)%>%
  select(species, dband,  abundance, abu.corr) %>%
  group_by(species, dband) %>%
  summarize(abu.corr = sum(abu.corr)) %>%
  mutate(cluster=case_when(dband<=70 ~ "upper mesophotic",
                           dband<=120~ "lower mesophotic",
                           dband<=180~ "upper rariphotic",
                           dband>180~ "lower rariphotic")) %>%
  group_by(cluster,species) %>%
  summarize(abu.corr=sum(abu.corr)) 

abundant.species <- data.frame(cluster=c(),species=c(),abu.corr=c(),prop=c(),
                                   totprop=c(),location=c())

for (k in unique(abu.sp$cluster)){
  cluster <- abu.sp[abu.sp$cluster==k,]
  abu.cluster <- cluster[order(cluster$abu.corr, decreasing = TRUE)[1:5],]
  abu.cluster$prop <- abu.cluster$abu.corr/sum(cluster$abu.corr)
  abu.cluster$totprop <- sum(abu.cluster$prop)
  abu.cluster$location <- "Curacao"
  abundant.species <- bind_rows(abundant.species,abu.cluster)
}

## Bonaire
abu.sp <- carib.deep.fishbase %>%
  filter(location=="Bonaire",dband>30,dband<310)%>%
  select(species, dband,  abundance, abu.corr) %>%
  group_by(species, dband) %>%
  summarize(abu.corr = sum(abu.corr)) %>%
  mutate(cluster=case_when(dband<=60 ~ "upper mesophotic",
                           dband<=110~ "lower mesophotic",
                           dband<=160~ "upper rariphotic",
                           dband>160~ "lower rariphotic")) %>%
  group_by(cluster,species) %>%
  summarize(abu.corr=sum(abu.corr)) 

abundant.species <- data.frame(cluster=c(),species=c(),abu.corr=c(),prop=c(),
                               totprop=c(),location=c())

for (k in unique(abu.sp$cluster)){
  cluster <- abu.sp[abu.sp$cluster==k,]
  abu.cluster <- cluster[order(cluster$abu.corr, decreasing = TRUE)[1:5],]
  abu.cluster$prop <- abu.cluster$abu.corr/sum(cluster$abu.corr)
  abu.cluster$totprop <- sum(abu.cluster$prop)
  abu.cluster$location <- "Bonaire"
  abundant.species <- bind_rows(abundant.species,abu.cluster)
}

## Statia
abu.sp <- carib.deep.fishbase %>%
  filter(location=="St. Eustatius",dband>30,dband<310)%>%
  select(species, dband,  abundance, abu.corr) %>%
  group_by(species, dband) %>%
  summarize(abu.corr = sum(abu.corr)) %>%
  mutate(cluster=case_when(dband<=90 ~ "upper mesophotic",
                           dband<=140~ "lower mesophotic",
                           dband<=200~ "upper rariphotic",
                           dband>200~ "lower rariphotic")) %>%
  group_by(cluster,species) %>%
  summarize(abu.corr=sum(abu.corr)) 

abundant.species <- data.frame(cluster=c(),species=c(),abu.corr=c(),prop=c(),
                               totprop=c(),location=c())

for (k in unique(abu.sp$cluster)){
  cluster <- abu.sp[abu.sp$cluster==k,]
  abu.cluster <- cluster[order(cluster$abu.corr, decreasing = TRUE)[1:5],]
  abu.cluster$prop <- abu.cluster$abu.corr/sum(cluster$abu.corr)
  abu.cluster$totprop <- sum(abu.cluster$prop)
  abu.cluster$location <- "Statia"
  abundant.species <- bind_rows(abundant.species,abu.cluster)
}


## Roatan
abu.sp <- carib.deep.fishbase %>%
  filter(location=="Roatan",dband>30,dband<310)%>%
  select(species, dband, abundance, abu.corr) %>%
  group_by(species, dband) %>%
  summarize(abu.corr = sum(abu.corr)) %>%
  mutate(cluster=case_when(dband<=100 ~ "upper mesophotic",
                           dband<=150~ "lower mesophotic",
                           dband<=210~ "upper rariphotic",
                           dband<= 300~ "lower rariphotic")) %>%
  group_by(cluster,species) %>%
  summarize(abu.corr=sum(abu.corr)) 

abundant.species <- data.frame(cluster=c(),species=c(),abu.corr=c(),prop=c(),
                               totprop=c())

for (k in unique(abu.sp$cluster)){
    cluster <- abu.sp[abu.sp$cluster==k,]
    abu.cluster <- cluster[order(cluster$abu.corr, decreasing = TRUE)[1:10],]
    abu.cluster$prop <- abu.cluster$abu.corr/sum(cluster$abu.corr)
    abu.cluster$totprop <- sum(abu.cluster$prop)
    abundant.species <- bind_rows(abundant.species,abu.cluster)
  }

#### 4. Looking at the distribution of species across sites #### 
#### 4.1. Bullisichtys caribbeaus ####
b.carib <- read.csv("2.file.fish.carib.grouped.csv") %>%
  filter(species=="Bullisichthys caribbaeus") %>%
  group_by(location,dband)%>%
  summarize(abu.corr=sum(abu.corr))%>%
  mutate(location=factor(location,levels=c("Curacao","Bonaire","St. Eustatius","Roatan")))

carib.pal <- fish(4, option = "Bodianus_rufus")

ggplot(data=b.carib,aes(x=dband, y=log(abu.corr),color=location))+
  geom_rect(xmin=130,xmax=80,ymin=7.5, ymax=8,fill=carib.pal[1],color="white")+
  geom_rect(xmin=70,xmax=120,ymin=8.5, ymax=9,fill=carib.pal[2],color="white")+
  geom_rect(xmin=100,xmax=140,ymin=9.5, ymax=10,fill=carib.pal[3],color="white")+
  geom_rect(xmin=100,xmax=150,ymin=10.5, ymax=11,fill=carib.pal[4],color="white")+
  geom_point(aes(shape = location),size=3) +
  scale_color_manual(values=carib.pal)+
  #coord_flip()+
  xlim(40,300)+
  theme_classic()+
  ylim(0,11)


#### 5.2. Pronotogrammus martinicensis ####
p.marti <- read.csv("2.file.fish.carib.grouped.csv") %>%
  filter(species=="Pronotogrammus martinicensis") %>%
  group_by(location,dband)%>%
  summarize(abu.corr=sum(abu.corr))%>%
  mutate(location=factor(location,levels=c("Curacao","Bonaire","St. Eustatius","Roatan")))

carib.pal <- fish(4, option = "Bodianus_rufus")

ggplot(data=p.marti,aes(x=dband, y=log(abu.corr),color=location))+
  geom_rect(xmin=130,xmax=190,ymin=7.5, ymax=8,fill=carib.pal[1],color="white")+
  geom_rect(xmin=170,xmax=120,ymin=8.5, ymax=9,fill=carib.pal[2],color="white")+
  geom_rect(xmin=190,xmax=140,ymin=9.5, ymax=10,fill=carib.pal[3],color="white")+
  #geom_rect(xmin=210,xmax=150,ymin=10.5, ymax=11,fill=carib.pal[4],color="white")+
  geom_point(aes(shape = location),size=3) +
  scale_color_manual(values=carib.pal)+
  #coord_flip()+
  xlim(40,300)+
  ylim(0,11)+
  theme_classic()

#### 5.3. Chrionema squamentum ####
c.squam <- read.csv("2.file.fish.carib.grouped.csv") %>%
  filter(species=="Chrionema squamentum") %>%
  group_by(location,dband)%>%
  summarize(abu.corr=sum(abu.corr))%>%
  mutate(location=factor(location,levels=c("Curacao","Bonaire","St. Eustatius","Roatan")))

carib.pal <- fish(4, option = "Bodianus_rufus")
rect.y <- log(max(c.squam$abu.corr))+0.5
  
ggplot(data=c.squam,aes(x=dband, y=log(abu.corr),color=location))+
  geom_rect(xmin=300,xmax=190,ymin=rect.y, ymax=rect.y+0.3,fill=carib.pal[1],color="white")+
  geom_rect(xmin=170,xmax=300,ymin=rect.y+0.6, ymax=rect.y+0.9,fill=carib.pal[2],color="white")+
  geom_rect(xmin=190,xmax=300,ymin=rect.y+1.2, ymax=rect.y+1.5,fill=carib.pal[3],color="white")+
  geom_rect(xmin=210,xmax=300,ymin=rect.y+1.8, ymax=rect.y+2.1,fill=carib.pal[4],color="white")+
  geom_point(aes(shape = location),size=3) +
  scale_color_manual(values=carib.pal)+
  #coord_flip()+
  xlim(40,400)+
  theme_classic()+
  ylim(0,8)

#### 5.4. Ostichthys trachypoma ####
o.trachy <- read.csv("2.file.fish.carib.grouped.csv") %>%
  filter(species=="Ostichthys trachypoma") %>%
  group_by(location,dband)%>%
  summarize(abu.corr=sum(abu.corr))%>%
  mutate(location=factor(location,levels=c("Curacao","Bonaire","St. Eustatius","Roatan")))

carib.pal <- fish(4, option = "Bodianus_rufus")
rect.y <- log(max(o.trachy$abu.corr))+0.5

ggplot(data=o.trachy,aes(x=dband, y=log(abu.corr),color=location))+
  geom_rect(xmin=300,xmax=190,ymin=rect.y, ymax=rect.y+0.3,fill=carib.pal[1],color="white")+
  geom_rect(xmin=170,xmax=300,ymin=rect.y+0.6, ymax=rect.y+0.9,fill=carib.pal[2],color="white")+
  geom_rect(xmin=190,xmax=300,ymin=rect.y+1.2, ymax=rect.y+1.5,fill=carib.pal[3],color="white")+
  geom_rect(xmin=210,xmax=300,ymin=rect.y+1.8, ymax=rect.y+2.1,fill=carib.pal[4],color="white")+
  geom_point(aes(shape = location),size=3) +
  scale_color_manual(values=carib.pal)+
  #coord_flip()+
  xlim(40,400)+
  theme_classic()+
  ylim(0,7)

#### 5.5. Symphysanodon berryi ####
s.berr <- read.csv("2.file.fish.carib.grouped.csv") %>%
  filter(species=="Symphysanodon berryi") %>%
  group_by(location,dband)%>%
  summarize(abu.corr=sum(abu.corr))%>%
  mutate(location=factor(location,levels=c("Curacao","Bonaire","St. Eustatius","Roatan")))

carib.pal <- fish(4, option = "Bodianus_rufus")
rect.y <- log(max(s.berr$abu.corr))+0.5

ggplot(data=s.berr,aes(x=dband, y=log(abu.corr),color=location))+
  geom_rect(xmin=300,xmax=190,ymin=rect.y, ymax=rect.y+0.3,fill=carib.pal[1],color="white")+
  #geom_rect(xmin=170,xmax=300,ymin=rect.y+0.6, ymax=rect.y+0.9,fill=carib.pal[2],color="white")+
  geom_rect(xmin=190,xmax=300,ymin=rect.y+1.2, ymax=rect.y+1.5,fill=carib.pal[3],color="white")+
  geom_rect(xmin=210,xmax=300,ymin=rect.y+1.8, ymax=rect.y+2.1,fill=carib.pal[4],color="white")+
  geom_point(aes(shape = location),size=3) +
  scale_color_manual(values=carib.pal[c(1,3,4)])+
  scale_shape_manual(values = c("Curacao" = 16,"St. Eustatius"=15, "Roatan" = 3))+
  #coord_flip()+
  xlim(40,400)+
  theme_classic()+
  ylim(0,7)

#### 5.6. Gephyroberyx darwinii ####
g.darw <- read.csv("2.file.fish.carib.grouped.csv") %>%
  filter(species=="Gephyroberyx darwinii") %>%
  group_by(location,dband)%>%
  summarize(abu.corr=sum(abu.corr))%>%
  mutate(location=factor(location,levels=c("Curacao","Bonaire","St. Eustatius","Roatan")))

carib.pal <- fish(4, option = "Bodianus_rufus")
rect.y <- log(max(g.darw$abu.corr))+0.5

ggplot(data=g.darw,aes(x=dband, y=log(abu.corr),color=location))+
  geom_rect(xmin=300,xmax=190,ymin=rect.y, ymax=rect.y+0.3,fill=carib.pal[1],color="white")+
  #geom_rect(xmin=170,xmax=300,ymin=rect.y+0.6, ymax=rect.y+0.9,fill=carib.pal[2],color="white")+
  #geom_rect(xmin=190,xmax=300,ymin=rect.y+1.2, ymax=rect.y+1.5,fill=carib.pal[3],color="white")+
  geom_rect(xmin=210,xmax=300,ymin=rect.y+0.6, ymax=rect.y+0.9,fill=carib.pal[4],color="white")+
  geom_point(aes(shape = location),size=3) +
  scale_color_manual(values=carib.pal[c(1,4)])+
  scale_shape_manual(values = c("Curacao" = 16, "Roatan" = 3))+
  #coord_flip()+
  xlim(40,500)+
  theme_classic()+
  ylim(0,4.5)

#### 5.6. Palatogobius incendius ####
p.inc <- read.csv("2.file.fish.carib.grouped.csv") %>%
  filter(species=="Palatogobius incendius") %>%
  group_by(location,dband)%>%
  summarize(abu.corr=sum(abu.corr))%>%
  mutate(location=factor(location,levels=c("Curacao","Bonaire","St. Eustatius","Roatan")))

carib.pal <- fish(4, option = "Bodianus_rufus")
rect.y <- log(max(p.inc$abu.corr))+0.5

ggplot(data=p.inc,aes(x=dband, y=log(abu.corr),color=location))+
  geom_rect(xmin=190,xmax=80,ymin=7.5, ymax=8,fill=carib.pal[1],color="white")+
  geom_rect(xmin=70,xmax=170,ymin=8.5, ymax=9,fill=carib.pal[2],color="white")+
  geom_rect(xmin=190,xmax=100,ymin=9.5, ymax=10,fill=carib.pal[3],color="white")+
  geom_rect(xmin=100,xmax=210,ymin=10.5, ymax=11,fill=carib.pal[4],color="white")+
  geom_point(aes(shape = location),size=3) +
  scale_color_manual(values=carib.pal)+
  #coord_flip()+
  xlim(40,300)+
  theme_classic()+
  ylim(0,10.5)

#### 5.7. Chromis insolata ####
c.inso <- read.csv("2.file.fish.carib.grouped.csv") %>%
  filter(species=="Chromis insolata") %>%
  group_by(location,dband)%>%
  summarize(abu.corr=sum(abu.corr))%>%
  mutate(location=factor(location,levels=c("Curacao","Bonaire","St. Eustatius","Roatan")))

carib.pal <- fish(4, option = "Bodianus_rufus")
rect.y <- log(max(c.inso$abu.corr))+0.5

ggplot(data=c.inso,aes(x=dband, y=log(abu.corr),color=location))+
  geom_rect(xmin=-40,xmax=-80,ymin=-2, ymax=-2.3,fill=carib.pal[1],color="white")+
  geom_rect(xmin=-40,xmax=-70,ymin=-2.5, ymax=-2.8,fill=carib.pal[2],color="white")+
  geom_rect(xmin=-40,xmax=-100,ymin=-3, ymax=-3.3,fill=carib.pal[3],color="white")+
  geom_rect(xmin=-100,xmax=-40,ymin=-3.5, ymax=-3.8,fill=carib.pal[4],color="white")+
  geom_point(aes(shape = location),size=3) +
  scale_color_manual(values=carib.pal)+
  coord_flip(clip = "off")+
  theme(axis.text.y=element_text(colour="black"),
        axis.line.y=element_line(colour="black"),
        axis.text.x=element_text(colour="black"),
        panel.background = element_blank(),
        plot.margin = unit(c(1,1,1,3),"lines"))+
  #ylim(0,8)+
  #coord_cartesian(clip = "off")+
  scale_x_reverse(limits = c(150,0), name ="")+
  scale_y_continuous(position="right",name = "abundance (log)",)

#### 5.8. Chromis cyanea ####
c.cyan <- read.csv("2.file.fish.carib.grouped.csv") %>%
  filter(species=="Chromis cyanea") %>%
  group_by(location,dband)%>%
  summarize(abu.corr=sum(abu.corr))%>%
  mutate(location=factor(location,levels=c("Curacao","Bonaire","St. Eustatius","Roatan")))

carib.pal <- fish(4, option = "Bodianus_rufus")
rect.y <- log(max(c.cyan$abu.corr))+0.5

ggplot(data=c.cyan,aes(x=dband, y=log(abu.corr),color=location))+
  geom_rect(xmin=-40,xmax=-80,ymin=-2, ymax=-2.3,fill=carib.pal[1],color="white")+
  geom_rect(xmin=-40,xmax=-70,ymin=-2.5, ymax=-2.8,fill=carib.pal[2],color="white")+
  geom_rect(xmin=-40,xmax=-100,ymin=-3, ymax=-3.3,fill=carib.pal[3],color="white")+
  geom_rect(xmin=-100,xmax=-40,ymin=-3.5, ymax=-3.8,fill=carib.pal[4],color="white")+
  geom_point(aes(shape = location),size=3) +
  scale_color_manual(values=carib.pal)+
  coord_flip(clip = "off")+
  theme(axis.text.y=element_text(colour="black"),
        axis.line.y=element_line(colour="black"),
        axis.text.x=element_text(colour="black"),
        panel.background = element_blank(),
        plot.margin = unit(c(1,1,1,7),"lines"))+
  #ylim(0,8)+
  #coord_cartesian(clip = "off")+
  scale_x_reverse(limits = c(100,0), name ="")+
  scale_y_continuous(position="right",name = "abundance (log)",)

#### 5.9. Stegastes partitus ####
s.parti <- read.csv("2.file.fish.carib.grouped.csv") %>%
  filter(species=="Stegastes partitus") %>%
  group_by(location,dband)%>%
  summarize(abu.corr=sum(abu.corr))%>%
  mutate(location=factor(location,levels=c("Curacao","Bonaire","St. Eustatius","Roatan")))

carib.pal <- fish(4, option = "Bodianus_rufus")
rect.y <- log(max(s.parti$abu.corr))+0.5

ggplot(data=s.parti,aes(x=dband, y=log(abu.corr),color=location))+
  geom_rect(xmin=-40,xmax=-80,ymin=-2, ymax=-2.3,fill=carib.pal[1],color="white")+
  geom_rect(xmin=-40,xmax=-70,ymin=-2.5, ymax=-2.8,fill=carib.pal[2],color="white")+
  geom_rect(xmin=-40,xmax=-100,ymin=-3, ymax=-3.3,fill=carib.pal[3],color="white")+
  geom_rect(xmin=-100,xmax=-40,ymin=-3.5, ymax=-3.8,fill=carib.pal[4],color="white")+
  geom_point(aes(shape = location),size=3) +
  scale_color_manual(values=carib.pal)+
  coord_flip(clip = "off")+
  theme(axis.text.y=element_text(colour="black"),
        axis.line.y=element_line(colour="black"),
        axis.text.x=element_text(colour="black"),
        panel.background = element_blank(),
        plot.margin = unit(c(1,1,1,3),"lines"))+
  #ylim(0,8)+
  #coord_cartesian(clip = "off")+
  scale_x_reverse(limits = c(100,0), name ="")+
  scale_y_continuous(position="right",name = "abundance (log)",)

#### 5.10. Antilligobius nikkiae ####
a.nik <- read.csv("2.file.fish.carib.grouped.csv") %>%
  filter(species=="Antilligobius nikkiae") %>%
  group_by(location,dband)%>%
  summarize(abu.corr=sum(abu.corr))%>%
  mutate(location=factor(location,levels=c("Curacao","Bonaire","St. Eustatius","Roatan")))

carib.pal <- fish(4, option = "Bodianus_rufus")
rect.y <- log(max(a.nik$abu.corr))+0.5

ggplot(data=a.nik,aes(x=dband, y=log(abu.corr),color=location))+
  geom_rect(xmin=130,xmax=80,ymin=7.5, ymax=8,fill=carib.pal[1],color="white")+
  geom_rect(xmin=70,xmax=120,ymin=8.5, ymax=9,fill=carib.pal[2],color="white")+
  geom_rect(xmin=100,xmax=140,ymin=9.5, ymax=10,fill=carib.pal[3],color="white")+
  geom_rect(xmin=100,xmax=150,ymin=10.5, ymax=11,fill=carib.pal[4],color="white")+
  geom_point(aes(shape = location),size=3) +
  scale_color_manual(values=carib.pal)+
  #coord_flip()+
  xlim(40,300)+
  theme_classic()+
  ylim(0,10.5)


#### 6. Sisters species - one species per graph ####
#### 6.1. Lipogramma  ####
lipogram <- read.csv("2.file.fish.carib.grouped.csv") %>%
  #filter(species=="Lipogramma haberi" | species=="Lipogramma evides"| species== "Lipogramma levinsoni") %>%
  filter(grepl("Lipogramma", species))%>%
  group_by(location,dband, species)%>%
  summarize(abu.corr=sum(abu.corr))%>%
  mutate(location=factor(location,levels=c("Curacao","Bonaire","St. Eustatius","Roatan")))

carib.pal <- fish(4, option = "Bodianus_rufus")
#rect.y <- log(max(a.nik$abu.corr))+0.5

ggplot(data=lipogram,aes(x=dband, y=log(abu.corr),color=location))+
  #geom_rect(xmin=130,xmax=80,ymin=7.5, ymax=8,fill=carib.pal[1],color="white")+
  #geom_rect(xmin=70,xmax=120,ymin=8.5, ymax=9,fill=carib.pal[2],color="white")+
  #geom_rect(xmin=100,xmax=140,ymin=9.5, ymax=10,fill=carib.pal[3],color="white")+
  #geom_rect(xmin=100,xmax=150,ymin=10.5, ymax=11,fill=carib.pal[4],color="white")+
  geom_point(size=3) +
  scale_color_manual(values=carib.pal)+
  facet_wrap(.~species,ncol=1)+
  #coord_flip()+
  xlim(40,300)+
  theme_classic()
  #ylim(0,4)

## only most common species 
lipogram <- read.csv("2.file.fish.carib.grouped.csv") %>%
  filter(species=="Lipogramma klayi" | species=="Lipogramma evides"| species== "Lipogramma levinsoni") %>%
  group_by(location,dband, species)%>%
  summarize(abu.corr=sum(abu.corr))%>%
  mutate(location=factor(location,levels=c("Curacao","Bonaire","St. Eustatius","Roatan")),
         species=factor(species,levels=c("Lipogramma klayi","Lipogramma levinsoni","Lipogramma evides")))

ggplot(data=lipogram,aes(x=dband, y=log(abu.corr),color=location))+
  #geom_rect(xmin=130,xmax=80,ymin=7.5, ymax=8,fill=carib.pal[1],color="white")+
  #geom_rect(xmin=70,xmax=120,ymin=8.5, ymax=9,fill=carib.pal[2],color="white")+
  #geom_rect(xmin=100,xmax=140,ymin=9.5, ymax=10,fill=carib.pal[3],color="white")+
  #geom_rect(xmin=100,xmax=150,ymin=10.5, ymax=11,fill=carib.pal[4],color="white")+
  geom_point(size=3) +
  scale_color_manual(values=carib.pal)+
  facet_wrap(.~species,ncol=1)+
  #coord_flip()+
  xlim(40,300)+
  theme_classic()



#### 6.2. Liopropoma ####
lioprop <- read.csv("2.file.fish.carib.grouped.csv") %>%
  filter(grepl("Liopropoma", species))%>%
  group_by(location,dband, species)%>%
  summarize(abu.corr=sum(abu.corr))%>%
  mutate(location=factor(location,levels=c("Curacao","Bonaire","St. Eustatius","Roatan")))

carib.pal <- fish(4, option = "Bodianus_rufus")
#rect.y <- log(max(a.nik$abu.corr))+0.5

ggplot(data=lioprop,aes(x=dband, y=log(abu.corr),color=location))+
  #geom_rect(xmin=130,xmax=80,ymin=7.5, ymax=8,fill=carib.pal[1],color="white")+
  #geom_rect(xmin=70,xmax=120,ymin=8.5, ymax=9,fill=carib.pal[2],color="white")+
  #geom_rect(xmin=100,xmax=140,ymin=9.5, ymax=10,fill=carib.pal[3],color="white")+
  #geom_rect(xmin=100,xmax=150,ymin=10.5, ymax=11,fill=carib.pal[4],color="white")+
  geom_point(size=3) +
  scale_color_manual(values=carib.pal)+
  facet_wrap(.~species,ncol=1)+
  #coord_flip()+
  xlim(40,300)+
  theme_classic()


## only most common species 
lioprop <- read.csv("2.file.fish.carib.grouped.csv") %>%
  filter(species=="Liopropoma aberrans" | species=="Liopropoma mowbrayi"| species== "Liopropoma olneyi") %>%
  group_by(location,dband, species)%>%
  summarize(abu.corr=sum(abu.corr))%>%
  mutate(location=factor(location,levels=c("Curacao","Bonaire","St. Eustatius","Roatan")),
         species=factor(species,levels=c("Liopropoma mowbrayi","Liopropoma aberrans","Liopropoma olneyi")))

ggplot(data=lioprop,aes(x=dband, y=log(abu.corr),color=location))+
  #geom_rect(xmin=130,xmax=80,ymin=7.5, ymax=8,fill=carib.pal[1],color="white")+
  #geom_rect(xmin=70,xmax=120,ymin=8.5, ymax=9,fill=carib.pal[2],color="white")+
  #geom_rect(xmin=100,xmax=140,ymin=9.5, ymax=10,fill=carib.pal[3],color="white")+
  #geom_rect(xmin=100,xmax=150,ymin=10.5, ymax=11,fill=carib.pal[4],color="white")+
  geom_point(size=3) +
  scale_color_manual(values=carib.pal)+
  facet_wrap(.~species,ncol=1)+
  #coord_flip()+
  xlim(40,300)+
  theme_classic()


#### 6.3. Prognathodes ####
prognathodes <- read.csv("2.file.fish.carib.grouped.csv") %>%
  filter(grepl("Prognathodes", species))%>%
  group_by(location,dband, species)%>%
  summarize(abu.corr=sum(abu.corr))%>%
  mutate(location=factor(location,levels=c("Curacao","Bonaire","St. Eustatius","Roatan")))

carib.pal <- fish(4, option = "Bodianus_rufus")
#rect.y <- log(max(a.nik$abu.corr))+0.5

ggplot(data=prognathodes,aes(x=dband, y=log(abu.corr),color=location))+
   geom_point(size=3) +
  scale_color_manual(values=carib.pal)+
  facet_wrap(.~species,ncol=1)+
  #coord_flip()+
  xlim(40,300)+
  theme_classic()

#### 6.4. Chromis ####
chromis <- read.csv("2.file.fish.carib.grouped.csv") %>%
  filter(grepl("Chromis", species))%>%
  group_by(location,dband, species)%>%
  summarize(abu.corr=sum(abu.corr))%>%
  mutate(location=factor(location,levels=c("Curacao","Bonaire","St. Eustatius","Roatan")),
         species=factor(species,levels=c("Chromis cyanea","Chromis insolata","Chromis multilineata",
                                         "Chromis scotti","Chromis enchrysura")))

carib.pal <- fish(4, option = "Bodianus_rufus")
#rect.y <- log(max(a.nik$abu.corr))+0.5

ggplot(data=chromis,aes(x=dband, y=log(abu.corr),color=location))+
  geom_point(size=3) +
  scale_color_manual(values=carib.pal)+
  facet_wrap(.~species,ncol=1)+
  #coord_flip()+
  xlim(0,300)+
  theme_classic()

#### 6.5. Palatogobius ####
gobius <- read.csv("2.file.fish.carib.grouped.csv") %>%
  filter(species=="Palatogobius incendius"| species=="Palatogobius grandoculus"|
         species=="Antilligobius nikkiae")%>%
  group_by(location,dband, species)%>%
  summarize(abu.corr=sum(abu.corr))%>%
  mutate(location=factor(location,levels=c("Curacao","Bonaire","St. Eustatius","Roatan")),
         species=factor(species,levels=c("Antilligobius nikkiae","Palatogobius incendius",
                                         "Palatogobius grandoculus")))

carib.pal <- fish(4, option = "Bodianus_rufus")
#rect.y <- log(max(a.nik$abu.corr))+0.5

ggplot(data=gobius,aes(x=dband, y=log(abu.corr),color=location))+
  geom_point(size=3) +
  scale_color_manual(values=carib.pal)+
  facet_wrap(.~species,ncol=1)+
  #coord_flip()+
  xlim(0,300)+
  theme_classic()


#### 7. Sister species - one genus per graph ####
#### 7.1. Lipogramma ####
lipogram <- read.csv("2.file.fish.carib.grouped.csv") %>%
  filter(species=="Lipogramma klayi" | species=="Lipogramma evides"| species== "Lipogramma levinsoni") %>%
  group_by(dband, species)%>%
  summarize(abu.corr=sum(abu.corr),abu.raw=sum(abundance))

lipogram.means <- lipogram %>%
  mutate(weighted.depth=dband*abu.corr) %>%
  group_by(species) %>%
  summarize(tot.depth=sum(weighted.depth),abu.corr=sum(abu.corr), 
            abu.raw=sum(abu.raw))%>%
  mutate(mean.depth=tot.depth/abu.corr)

## stat test
lipogram.stat <- read.csv("2.file.fish.carib.grouped.csv") %>%
  filter(species=="Lipogramma klayi" | species=="Lipogramma evides"| species== "Lipogramma levinsoni")

lipogram.aov <- aov(data=lipogram, dband~species)
summary(lipogram.aov) ## p <0.001
TukeyHSD(lipogram.aov)

## plot
carib.pal <- fish(3, option = "Bodianus_rufus")
lipogram.plot <- ggplot(data=lipogram,aes(x=dband, y=log(abu.corr),color=species))+
  geom_vline(xintercept=lipogram.means$mean.depth,color=carib.pal)+
  geom_point(size=2) +
  scale_color_manual(values=carib.pal)+
  coord_flip()+
  theme_classic()+
  scale_x_reverse(limits=c(300,0))+
  xlab("depth")+
  ylab("abundance (log)")+
  ggtitle("Lipogramma")

lipogram.plot

#### 7.2. Liopropoma ####
## L.carmabi (n=5) kept to add to collected speciments
## L. santi (n=4) so removed 
lioprop <- read.csv("2.file.fish.carib.grouped.csv") %>%
filter(grepl("Liopropoma", species))%>%
 filter(species!="Liopropoma santi")%>%
  group_by(dband, species)%>%
  summarize(abu.corr=sum(abu.corr),abu.raw=sum(abundance))

lioprop.means <- lioprop %>%
  mutate(weighted.depth=dband*abu.corr) %>%
  group_by(species) %>%
  summarize(tot.depth=sum(weighted.depth),abu.corr=sum(abu.corr), 
                          abu.raw=sum(abu.raw))%>%
  mutate(mean.depth=tot.depth/abu.corr)

## stat test
lioprop.aov <- aov(data=lioprop, dband~species)
summary(lioprop.aov) ## p <0.001
TukeyHSD(lioprop.aov)

## plot
carib.pal <- fish(4, option = "Bodianus_rufus")
lioprop.plot <- ggplot(data=lioprop,aes(x=dband, y=log(abu.corr),color=species))+
  geom_vline(xintercept=lioprop.means$mean.depth,color=carib.pal)+
  geom_point(size=2) +
  scale_color_manual(values=carib.pal)+
  coord_flip()+
  theme_classic()+
  scale_x_reverse(limits=c(300,0))+
  xlab("depth")+
  ylab("abundance (log)")+
  ggtitle("Liopropoma")

lioprop.plot

#### 7.3. Prognathodes ####
prognathodes <- read.csv("2.file.fish.carib.grouped.csv") %>%
  filter(grepl("Prognathodes", species))%>%
  group_by(dband, species)%>%
  summarize(abu.corr=sum(abu.corr),abu.raw=sum(abundance))

prognathodes.means <- prognathodes %>%
  mutate(weighted.depth=dband*abu.corr) %>%
  group_by(species) %>%
  summarize(tot.depth=sum(weighted.depth),abu.corr=sum(abu.corr), 
            abu.raw=sum(abu.raw))%>%
  mutate(mean.depth=tot.depth/abu.corr)

## stat test
t.test(data=prognathodes, dband~species) ## p =0.002

## plot
carib.pal <- fish(2, option = "Bodianus_rufus")
prognathodes.plot <- ggplot(data=prognathodes,aes(x=dband, y=log(abu.corr),color=species))+
  geom_vline(xintercept=prognathodes.means$mean.depth,color=carib.pal)+
  geom_point(size=2) +
  scale_color_manual(values=carib.pal)+
  coord_flip()+
  theme_classic()+
  scale_x_reverse(limits=c(300,0))+
  xlab("depth")+
  ylab("abundance (log)")+
  ggtitle("Prognathodes")

prognathodes.plot


#### 7.4. Serranus ####
## 4 species with n>10: luciopercanus, nostospilus, phoebe and tortugarum
serranus <- read.csv("2.file.fish.carib.grouped.csv") %>%
  filter(species=="Serranus luciopercanus"| species=="Serranus notospilus" |
           species=="Serranus phoebe"|species=="Serranus tortugarum")%>%
  group_by(dband, species)%>%
  summarize(abu.corr=sum(abu.corr),abu.raw=sum(abundance))

serranus.means <- serranus %>%
  mutate(weighted.depth=dband*abu.corr) %>%
  group_by(species) %>%
  summarize(tot.depth=sum(weighted.depth),abu.corr=sum(abu.corr), 
            abu.raw=sum(abu.raw))%>%
  mutate(mean.depth=tot.depth/abu.corr)

## stat test
serranus.aov <- aov(data=serranus, dband~species)
summary(serranus.aov) ## p <0.001
TukeyHSD(serranus.aov)

## plot
carib.pal <- fish(4, option = "Bodianus_rufus")
serranus.plot <- ggplot(data=serranus,aes(x=dband, y=log(abu.corr),color=species))+
  geom_vline(xintercept=serranus.means$mean.depth,color=carib.pal)+
  geom_point(size=2) +
  scale_color_manual(values=carib.pal)+
  coord_flip()+
  theme_classic()+
  scale_x_reverse(limits=c(300,0))+
  xlab("depth")+
  ylab("abundance (log)")+
  ggtitle("Serranus")
serranus.plot


#### 7.5. Palatogobius ####
gobius <- read.csv("1.Ch3.R.analysis/1.file.fish.carib.grouped.csv") %>%
  filter(species=="Palatogobius incendius"| species=="Palatogobius grandoculus")%>%
  group_by(dband, species)%>%
  summarize(abu.corr=sum(abu.corr),abu.raw=sum(abundance))

gobius.means <- gobius %>%
  mutate(weighted.depth=dband*abu.corr) %>%
  group_by(species) %>%
  summarize(tot.depth=sum(weighted.depth),abu.corr=sum(abu.corr), 
            abu.raw=sum(abu.raw))%>%
  mutate(mean.depth=tot.depth/abu.corr)

## stat test
t.test(data=gobius, dband~species)

## plot
carib.pal <- fish(2, option = "Bodianus_rufus")
gobius.plot <- ggplot(data=gobius,aes(x=dband, y=log(abu.corr),color=species))+
  geom_vline(xintercept=gobius.means$mean.depth,color=carib.pal)+
  geom_point(size=2) +
  scale_color_manual(values=carib.pal)+
  coord_flip()+
  theme_classic()+
  scale_x_reverse(limits=c(300,0))+
  xlab("depth")+
  ylab("abundance (log)")+
  ggtitle("Palatogobius")

gobius.plot

#### 7.6 Symphysanodon ####
symphy <- read.csv("2.file.fish.carib.grouped.csv") %>%
  filter(grepl("Symphysanodon", species))%>%
  group_by(dband, species)%>%
  summarize(abu.corr=sum(abu.corr),abu.raw=sum(abundance))

symphy.means <- symphy %>%
  mutate(weighted.depth=dband*abu.corr) %>%
  group_by(species) %>%
  summarize(tot.depth=sum(weighted.depth),abu.corr=sum(abu.corr), 
            abu.raw=sum(abu.raw))%>%
  mutate(mean.depth=tot.depth/abu.corr)

## stat test
t.test(data=symphy, dband~species)

## plot
carib.pal <- fish(2, option = "Bodianus_rufus")
symphy.plot <- ggplot(data=symphy,aes(x=dband, y=log(abu.corr),color=species))+
  geom_vline(xintercept=symphy.means$mean.depth,color=carib.pal)+
  geom_point(size=2) +
  scale_color_manual(values=carib.pal)+
  coord_flip()+
  theme_classic()+
  scale_x_reverse(limits=c(300,0))+
  xlab("depth")+
  ylab("abundance (log)")+
  ggtitle("Symphysanodon")

symphy.plot

#### 7.7. Chromis ####
## C. scotti removed because of possibly wrong IDs
chromis <- read.csv("2.file.fish.carib.grouped.csv") %>%
  filter(grepl("Chromis", species))%>%
  filter(species!="Chromis scotti")%>%
  group_by(dband, species)%>%
  summarize(abu.corr=sum(abu.corr),abu.raw=sum(abundance))

chromis.means <- chromis %>%
  mutate(weighted.depth=dband*abu.corr) %>%
  group_by(species) %>%
  summarize(tot.depth=sum(weighted.depth),abu.corr=sum(abu.corr), 
            abu.raw=sum(abu.raw))%>%
  mutate(mean.depth=tot.depth/abu.corr)

## stat test
chromis.aov <- aov(data=chromis, dband~species)
summary(chromis.aov) ## p <0.001
TukeyHSD(chromis.aov)

## plot
carib.pal <- fish(4, option = "Bodianus_rufus")
chromis.plot <- ggplot(data=chromis,aes(x=dband, y=log(abu.corr),color=species))+
  geom_vline(xintercept=chromis.means$mean.depth,color=carib.pal)+
  geom_point(size=2) +
  scale_color_manual(values=carib.pal)+
  coord_flip()+
  theme_classic()+
  #scale_x_reverse()+
  scale_x_reverse(limits=c(300,0))+
  xlab("depth")+
  ylab("abundance (log)")+
  ggtitle("Chromis")

chromis.plot

#### 7.8. Haemulon ####
haemulon <- read.csv("2.file.fish.carib.grouped.csv") %>%
  filter(grepl("Haemulon", species))%>%
  filter(species!="Haemulon bonariense" & species!="Haemulon plumierii" &
           species!="Haemulon sciurus")%>% #less than 10 observations
  group_by(dband, species)%>%
  summarize(abu.corr=sum(abu.corr),abu.raw=sum(abundance))

haemulon.means <- haemulon %>%
  mutate(weighted.depth=dband*abu.corr) %>%
  group_by(species) %>%
  summarize(tot.depth=sum(weighted.depth),abu.corr=sum(abu.corr), 
            abu.raw=sum(abu.raw))%>%
  mutate(mean.depth=tot.depth/abu.corr)

carib.pal <- fish(4, option = "Bodianus_rufus")

haemulon.plot <- ggplot(data=haemulon,aes(x=dband, y=log(abu.corr),color=species))+
  geom_vline(xintercept=haemulon.means$mean.depth,color=carib.pal)+
  geom_point(size=2) +
  scale_color_manual(values=carib.pal,name=NULL,
                     labels = c("H. chrysargyreum","H. flavolineatym",
                                "H. striatum","H. vittata"))+
  coord_flip()+
  theme_classic()+
  #scale_x_reverse()+
  scale_x_reverse(limits=c(300,0))+
  xlab("depth")+
  ylab("abundance (log)")+
  ggtitle("Haemulon")

haemulon.plot

#### 7.9. Plots all together ####
genus.plots <- ggarrange(chromis.plot, haemulon.plot+ rremove("ylab"),lioprop.plot+ rremove("ylab"),
                           prognathodes.plot+ rremove("ylab"),serranus.plot,
                         lipogram.plot+ rremove("ylab"),gobius.plot+ rremove("ylab"),symphy.plot+ rremove("ylab"),
                          ncol = 4, nrow=2, legend="none")
genus.plots
