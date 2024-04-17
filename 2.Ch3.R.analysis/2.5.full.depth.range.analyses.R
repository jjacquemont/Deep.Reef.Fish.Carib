library(tidyverse) #
library(modelr) #part of the tidyverse ecosystem and provides tools for creating and working with data grids, which are often used in modeling and data analysis tasks.
library(fishualize) # provides color scales based on fish color
library(vegan) # ecological and environmental data analysis
library(ggrepel)# extension of the popular ggplot2 package and provides additional functionality for controlling the placement and positioning of text labels in a ggplot2 plot.
library(cluster) # alternative package to clustsig
library(ggdendro) # to perform dendograms
library(tinter) # generates palettes / tints for graphs
library(mgcv) # mixed generalized additive models 
library(FD) # functional diversity package
library(GGally) #functions and features that complement the functionality of the popular ggplot2 package.
#library(devtools) #to enable the download of an archived version of clustsig
#install_github("cran/clustsig")
library(clustsig)
library(hillR) #to calculate richness, evenness and dominance 
library(betapart) #to calculate beta diversity 

## This script analyzes the community structure for Roatan and Curacao for which 
# data is available above 30 m and below 300 m in the case of Roatan.

#### 1.1. Roatan dendrogram ####
roa.com.full <- read.csv("2.Ch3.R.analysis/2.file.fish.binned.csv") %>%
  filter(location == "Roatan") %>%
  select(species, dband, location, abundance, abu.corr) %>%
  mutate(dband.grouped = case_when(dband >= 310 & dband <= 330 ~ 310,
                                   dband >= 350 & dband <= 360 ~ 350,
                                   dband >= 390 & dband <= 400 ~ 390,
                                   dband >= 410 & dband <= 430 ~ 410,
                                   dband >= 440 & dband <= 460 ~ 440,
                                   dband >= 470 & dband <= 480 ~ 470,
                                   TRUE ~ as.numeric(dband))) %>%
  group_by(species, dband.grouped) %>%
  summarize(abundance = sqrt(sum(abu.corr))) %>%
  spread(species, abundance, fill =  0) %>%
  as.data.frame()

rownames(roa.com.full) <- roa.com.full$dband.grouped
roa.com.cl.full <- roa.com.full[-1]

roa.clust.full <- simprof(data = roa.com.cl.full, method.cluster = "complete", method.distance="braycurtis",
                          method.transform="identity", alpha=0.0000001,
                          sample.orientation="row", const=0,
                          silent=TRUE, increment=100,
                          undef.zero=TRUE, warn.braycurtis=TRUE)
roa.clust.full

# Roatan: 6 groups
roa.res.tib.full <- tibble(cl1 = list(roa.clust.full$significantclusters[[1]]),
                           cl2 = list(roa.clust.full$significantclusters[[2]]),
                           cl3 = list(roa.clust.full$significantclusters[[3]]),
                           cl4 = list(roa.clust.full$significantclusters[[4]]),
                           cl5 = list(roa.clust.full$significantclusters[[5]]),
                           cl6 = list(roa.clust.full$significantclusters[[6]])) %>%
  unnest_longer(cl1) %>%
  unnest_longer(cl2) %>%
  unnest_longer(cl3) %>%
  unnest_longer(cl4) %>%
  unnest_longer(cl5) %>%
  unnest_longer(cl6) %>%
  gather(key = "cluster", value = "dband.grouped") %>%
  distinct() %>%
  mutate(dband.grouped = as.numeric(dband.grouped))%>%
  full_join(roa.com.full) %>%
  select(cluster, dband.grouped) %>%
  mutate(label = as.numeric(dband.grouped)) %>%
  arrange(dband.grouped) %>%
  mutate(run_num = 1:36) %>%
  group_by(cluster) %>%
  mutate(ord_num = min(run_num))

# hclust 
roa.clust.veg.full <- hclust(vegdist(roa.com.cl.full, method = "bray"), method = "complete")
roa.dend.full <- as.dendrogram(roa.clust.veg.full)
roa.dend.data.full <- dendro_data(roa.dend.full, type = "rectangle")

roa.dend.helper.full <- roa.dend.data.full$labels %>%
  mutate(label = as.numeric(label))%>%
  inner_join(roa.res.tib.full)


roa.dend.data.full$labels <- roa.dend.helper.full %>%
  mutate(label = case_when(dband.grouped==310 ~ "310-330",
                           dband.grouped == 350  ~ "350-360",
                           dband.grouped == 390  ~ "390-400",
                           dband.grouped == 410~ "410-430",
                           dband.grouped == 440 ~ "440-460",
                           dband.grouped == 470 ~ "470-480",
                           TRUE ~ as.character(label)))


# get colors for each level using fishualize and tinter package
carib.pal <- fish(4, option = "Bodianus_rufus")
roa.pal <- tinter(carib.pal[4], steps = 5)
blue_palette <- c("#BAE0F3","#87CEEB","#6CBDE9","#50ABE7", "#4895EF", "#4361EE", "#2835AF", "#12086F")
mesoP_palette <- c("#DCEEF3","#C2E2EA","#A7D5E1","#8DC8D8","#72BBCE")
rariP_palette <- c("#DDC5EB","#BC90DB","#831EB6")
pal <- c(blue_palette[3],blue_palette[3],rariP_palette[1],rariP_palette[2])

roa.plot.full <- ggplot(roa.dend.data.full$segments, horiz = TRUE) + 
  geom_rect(xmin=26.7,xmax=31.3,ymin=-0.25, ymax=0.95,fill=blue_palette[3])+
  geom_rect(xmin=31.7,xmax=36.3,ymin=-0.25, ymax=0.95,fill=rariP_palette[1])+
  geom_rect(xmin=21.7,xmax=26.3,ymin=-0.25, ymax=0.95,fill=blue_palette[1])+
  geom_rect(xmin=17.7,xmax=20.3,ymin=-0.25, ymax=0.95,fill=blue_palette[1])+
  geom_rect(xmin=20.7,xmax=21.3,ymin=-0.25, ymax=0.95,fill="#d0dc36")+
  geom_rect(xmin=7.7,xmax=17.3,ymin=-0.25, ymax=0.95,fill=rariP_palette[2])+
  geom_rect(xmin=0.3,xmax=7.3,ymin=-0.25, ymax=0.95,fill="#831EB6")+
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), color = "black")+
  geom_text(data = roa.dend.data.full$labels, aes(x, y, label = label),color = "black",
            hjust = 1, angle = 90, size = 3, fontface = "bold") +
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
  geom_segment(aes(x=0.3,xend=7.3,y=-0.25,yend=-0.25),size=1)+
  geom_segment(aes(x=7.8,xend=17.3,y=-0.25,yend=-0.25),size=1)+
  geom_segment(aes(x=17.8,xend=20.3,y=-0.25,yend=-0.25),size=1)+
  geom_segment(aes(x=20.8,xend=21.3,y=-0.25,yend=-0.25),size=1)+
  geom_segment(aes(x=21.8,xend=26.3,y=-0.25,yend=-0.25),size=1)+
  geom_segment(aes(x=26.8,xend=36.3,y=-0.25,yend=-0.25),size=1)+
  ggtitle("Roatan")+
  ylim(-0.5,1)
roa.plot.full

### 1.2. Roatan pca plot ####
roa.pca.full <- prcomp(roa.com.full[-1], center = TRUE, scale = TRUE)
summary(roa.pca.full)

roa.points <- 
  # first convert the pca results to a tibble
  as_tibble(roa.pca.full$x) %>% 
  bind_cols(roa.com.full) %>%
  inner_join(roa.dend.helper.full)  %>%
  mutate(cat=case_when(dband.grouped<20~ "altiphotic",
                        dband.grouped<100 ~ "upper mesophotic",
                       dband.grouped<150 ~ "lower mesophotic",
                       dband.grouped<210~ "upper rariphotic",
                       dband.grouped<300~ "lower rariphotic",
                       TRUE~ "deep-sea")) %>%
  mutate(cat=factor(cat,levels=c("altiphotic","upper mesophotic","lower mesophotic",
                                 "upper rariphotic","lower rariphotic","deep-sea")))


roa.hull <- 
  roa.points %>% 
  group_by(cluster) %>% 
  slice(chull(PC1, PC2))%>%
  mutate(cat=case_when(dband.grouped<20~ "altiphotic",
                       dband.grouped<100 ~ "upper mesophotic",
                       dband.grouped<150 ~ "lower mesophotic",
                       dband.grouped<210~ "upper rariphotic",
                       dband.grouped<300~ "lower rariphotic",
                       TRUE~ "deep-sea")) %>%
  mutate(cat=factor(cat,levels=c("altiphotic","upper mesophotic","lower mesophotic",
                                 "upper rariphotic","lower rariphotic","deep-sea")))


roa.simper <- with(roa.points, simper(roa.com.full[-1], cluster))
roa.simsum <- summary(roa.simper)

#### SIMPER between depth zones ####
roa.simper.zones <- simper(roa.com.full, roa.points$cat) 
roa.simsum.zones <- summary(roa.simper.zones)

# species driving differences between lower rariphotic and below
simper.deep <- roa.simsum.zones[["lower rariphotic_deep-sea"]]


roa.simspecies.zones <- data_frame("species" = as.character())
for (i in 1:length(roa.simsum.zones)){ # iterates for each element of list
  spec <- as.data.frame(roa.simsum.zones[i]) %>%
    mutate(species = rownames(.)) %>%
    select(c(1,7,8)) %>% # "average","p-value", "species"
    filter(.[[2]] < 0.05) %>%
    arrange(desc(.[[1]]))
  roa.simspecies.zones <- bind_rows(roa.simspecies.zones, spec)
}
roa.simspecies.zones$location="Roatan"


roa.species <- roa.simspecies %>%
  distinct()

roa.loads <- 
  as_tibble(roa.pca.full$rotation, rownames = 'species') %>%
  inner_join(roa.species) %>%
  select(species, PC1, PC2)


pal.roa=c("#DCEEF3",blue_palette[1],blue_palette[3],rariP_palette[1],"#831EB6","#12086F")
roa.pca.plot <- ggplot(roa.points) +
  geom_jitter(data = roa.points, aes(x = PC1, y = PC2, fill = cat), 
              shape = 21, alpha = 0.5, color = "black") +
  geom_polygon(data = roa.hull,
               aes(x = PC1, y = PC2, fill = cat,
                   colour = cat,
                   group = as.factor(ord_num)),
               alpha = 0.3,
               show.legend = TRUE) +
  geom_segment(data = roa.loads[c(1:7),], 
               aes(x = 0, y = 0, 
                   xend = PC1*30,
                   yend = PC2*30), lwd = 0.1, 
               arrow = arrow(length = unit(1/2, 'picas'), type = "open")) +
  geom_text_repel(data = roa.loads[c(1:7),], aes(x = PC1*30, 
                                        y = PC2*30, label = species), size = 2, segment.color = "grey69") +
  theme_classic() +
  scale_fill_manual(values = pal.roa) +
  scale_color_manual(values = pal.roa)
#guides(colour = guide_legend(nrow = 1), fill = guide_legend(nrow = 1))
roa.pca.plot

#### 1.3. Zoom on Roatan PCA for rariphotic and deep sea ####
roa.com.deep <- read.csv("2.Ch3.R.analysis/2.file.fish.carib.grouped.csv") %>%
  filter(location == "Roatan" & dband>150) %>%
  select(species, dband, location, abundance, abu.corr) %>%
  mutate(dband.grouped = case_when(dband >= 310 & dband <= 330 ~ 310,
                                   dband >= 350 & dband <= 360 ~ 350,
                                   dband >= 390 & dband <= 400 ~ 390,
                                   dband >= 410 & dband <= 430 ~ 410,
                                   dband >= 440 & dband <= 460 ~ 440,
                                   dband >= 470 & dband <= 480 ~ 470,
                                   TRUE ~ as.numeric(dband))) %>%
  group_by(species, dband.grouped) %>%
  summarize(abundance = sqrt(sum(abu.corr))) %>%
  spread(species, abundance, fill =  0) %>%
  as.data.frame()

rownames(roa.com.deep) <- roa.com.deep$dband.grouped
roa.com.cl.deep <- roa.com.deep[-1]
roa.com.deep <- roa.com.deep%>%
  filter(dband.grouped>150)

roa.pca.deep <- prcomp(roa.com.deep[-1], center = TRUE, scale = TRUE)
summary(roa.pca.deep)

roa.points <- 
  # first convert the pca results to a tibble
  as_tibble(roa.pca.deep$x) %>% 
  bind_cols(roa.com.deep) %>%
  inner_join(roa.dend.helper.full)  %>%
  filter(dband.grouped>150) %>%
  mutate(cat=case_when(dband.grouped<210~ "upper rariphotic",
                                dband.grouped<300~ "lower rariphotic",
                                TRUE~ "deep-sea")) %>%
           mutate(cat=factor(cat,levels=c("upper rariphotic","lower rariphotic","deep-sea")))


roa.hull <- 
  roa.points %>% 
  group_by(cluster) %>% 
  slice(chull(PC1, PC2))%>%
  mutate(cat=case_when(dband.grouped<210~ "upper rariphotic",
                       dband.grouped<300~ "lower rariphotic",
                       TRUE~ "deep-sea")) %>%
  mutate(cat=factor(cat,levels=c("upper rariphotic","lower rariphotic","deep-sea")))


roa.simper <- with(roa.points, simper(roa.com.deep[-1], cluster))
roa.simsum <- summary(roa.simper)

roa.simspecies <- data_frame("species" = as.character())

for (i in 1:length(roa.simsum)){
  spec <- as.data.frame(roa.simsum[i]) %>%
    mutate(species = rownames(.)) %>%
    select(1,8) %>%
    filter(.[[1]] > 0.05) %>% 
    select(species)
  
  roa.simspecies <- bind_rows(roa.simspecies, spec)
  
}

roa.species <- roa.simspecies %>%
  distinct()

roa.loads <- 
  as_tibble(roa.pca.deep$rotation, rownames = 'species') %>%
  inner_join(roa.species) %>%
  select(species, PC1, PC2)


pal.roa=c(rariP_palette[1],"#831EB6","#12086F")
roa.pca.plot <- ggplot(roa.points) +
  geom_jitter(data = roa.points, aes(x = PC1, y = PC2, fill = cat), 
              shape = 21, alpha = 0.5, color = "black") +
  geom_polygon(data = roa.hull,
               aes(x = PC1, y = PC2, fill = cat,
                   colour = cat,
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
  theme_classic() +
  scale_fill_manual(values = pal.roa) +
  scale_color_manual(values = pal.roa)
#guides(colour = guide_legend(nrow = 1), fill = guide_legend(nrow = 1))
roa.pca.plot

#### 1.4. MDS plot ####
####  MDS
mds_roa <- cmdscale(vegdist(roa.com.cl.full, method = "bray")) %>%
  as_tibble() %>%
  mutate(group=roa.points$cat,dband=roa.points$dband.grouped)
colnames(mds_roa) <- c("Dim.1", "Dim.2","cluster","dband")

pal <- c("#d0dc36","#BAE0F3","#6CBDE9","#DDC5EB","#BC90DB","#831EB6")
roa.mds.plot <- ggscatter(mds_roa, x = "Dim.1", y = "Dim.2", 
                          palette = pal,
                          size = 1.5,
                          label=mds_roa$dband,
                          fill="cluster",
                          color = "black",
                          ellipse.alpha = 0.7,
                          ellipse = TRUE,
                          ellipse.type = "convex",
                          repel = TRUE)
roa.mds.plot

#### 2. Contribution of depth affinity groups ####
## for 0-500 m
roa.com.all <- read.csv("2.Ch3.R.analysis/2.file.fish.carib.grouped.csv") %>%
  filter(location == "Roatan")

species.distribution.roa.all <- data.frame(species=c(),alti= c(),u.meso=c(),l.meso=c(),
                                           u.rari=c(),l.rari=c(),deep.sea=c(),tot.abu=c())
for (k in unique(roa.com.all$species)){
  species <- roa.com.all[roa.com.all$species==k,]
  alti <- sum(species$abu.corr[species$dband<40])/sum(species$abu.corr)
  u.meso <- sum(species$abu.corr[species$dband<100 & species$dband>=40])/sum(species$abu.corr)
  l.meso <- sum(species$abu.corr[species$dband<150 & species$dband>=100])/sum(species$abu.corr)
  u.rari <- sum(species$abu.corr[species$dband<210 & species$dband>=150])/sum(species$abu.corr)
  l.rari <- sum(species$abu.corr[species$dband<300 & species$dband>=210])/sum(species$abu.corr)
  deep.sea <- sum(species$abu.corr[species$dband>300])/sum(species$abu.corr)
  
  tot.abu <- sum(species$abu.corr)
  species.distribution.roa.all <- rbind(species.distribution.roa.all,
                                        c(k,alti,u.meso,l.meso,u.rari,l.rari,deep.sea,tot.abu))
}

colnames(species.distribution.roa.all) <- c("species","alti","u.meso","l.meso","u.rari","l.rari","deep.sea","tot.abu")
species.distribution.roa.all <- species.distribution.roa.all %>%
  mutate(alti=as.numeric(alti),u.meso=as.numeric(u.meso),l.meso=as.numeric(l.meso),u.rari=as.numeric(u.rari),
         l.rari=as.numeric(l.rari),deep.sea=as.numeric(deep.sea),tot.abu=as.numeric(tot.abu)) %>%
  #add number so that species can be ordered by depth predominance
  arrange(desc(deep.sea),alti,u.meso,l.meso,u.rari,l.rari) %>%
  mutate(species=factor(species,levels=species)) %>%
  #stacking columns for ggplot
  gather(key="depthzone", value="abundance",alti, u.meso,l.meso,u.rari,l.rari,deep.sea,tot.abu)  %>%
  mutate(location="Roatan",
         depthzone=factor(depthzone,levels=c("alti","u.meso","l.meso","u.rari","l.rari","deep.sea","tot.abu"))) 

#### plot 
rel.abundance.roa.all <- species.distribution.roa.all %>%
  filter(depthzone!="tot.abu")

ggplot(rel.abundance.roa.all, aes(y=species,x=depthzone))+
  geom_tile(aes(fill=abundance))+
  xlab("relative abundance per depth zone")+
  scale_x_discrete(position="top")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),axis.text = element_text(size = 5))+
  ylab("")

#### 2.2. Species present below the rariphotic ####
roa.deep <- read.csv("2.Ch3.R.analysis/2.file.fish.carib.grouped.csv") %>%
  filter(location == "Roatan" & dband>300) %>%
  group_by(species) %>%
  summarize(abu.corr = sum(abu.corr),abu.raw=sum(abundance)) 

#### 2.3. Contribution of depth affinities ####

####  2.3.1. abundance ####
roa.abu.com <- roa.com.all[-1] %>%
  inner_join(read.csv("2.Ch3.R.analysis/2.file.species.depth.csv"),by="species") %>%
  group_by(distrib.roa, dband) %>%
  summarize(abu.corr = sum(abu.corr)) %>%
  mutate(distrib.roa=factor(distrib.roa,levels=c("AM","M","RM","R","DS")))

palette <- viridis(5,begin=0.9,end=0)

ggplot(data=roa.abu.com,aes(x=dband,y=abu.corr,color=distrib.roa))+
  geom_point(aes(shape = distrib.roa),size=3) +
  xlab("depth")+
  ylab("abundance")+
  ggtitle("Roatan")+
  scale_color_manual(values = palette[1:5])+
  #scale_color_discrete(labels=c("alti/mesophotic","mesophotic","meso/rariphotic",
  #                    "rariphotic","deep sea"),name="predominant depth")+
  #guides(shape = "none")+
  theme_classic()+
  geom_vline(xintercept = 100, lty = 2, color = "grey69",size=0.3) +
  geom_vline(xintercept = 150, lty = 2, color = "grey69",size=0.3) +
  geom_vline(xintercept = 210, lty = 2, color = "grey69",size=0.3) +
  coord_flip()+
  scale_x_reverse()

#### Barplot
roa.abu.clus <- roa.abu.com %>%
  mutate(cluster=case_when(dband<40 ~ "altiphotic",
                           dband<100 ~ "upper mesophotic",
                           dband<150 ~ "lower mesophotic",
                           dband<210~ "upper rariphotic",
                           dband<300~ "lower rariphotic",
                           TRUE ~ "deep sea")) %>%
  group_by(cluster,distrib.roa) %>%
  summarize(abu.corr=sum(abu.corr))%>%
  mutate(distrib.roa=factor(distrib.roa,levels=c("AM","M","RM","R","DS")))%>%
  mutate(cluster=factor(cluster,levels=rev(c("altiphotic","upper mesophotic","lower mesophotic",
                                             "upper rariphotic","lower rariphotic","deep sea"))))

##bar plot
pal <- c("AM"="#d0dc36","M"="#50ABE7","RM"="#BC90DB","R"="#831EB6","DS"="#12086F")

ggplot(data=roa.abu.clus, aes(x=cluster,y=sqrt(abu.corr),fill=distrib.roa))+
  geom_bar(position="stack",stat="identity",width=0.3) +
  scale_fill_manual(values=pal,
                     labels=c("alti/mesophotic","mesophotic","meso/rariphotic",
                              "rariphotic","deep sea"),name="predominant depth")+
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
  ylab("total abundance (sqrt)")+
  coord_flip()

#### 2.3.2. richness ####
roa.spric.com <- roa.com.all[-1] %>%
  inner_join(read.csv("2.Ch3.R.analysis/2.file.species.depth.csv"),by="species") %>%
  group_by(distrib.roa, dband) %>%
  summarize(spric = n())%>%
  mutate(distrib.roa=factor(distrib.roa,levels=c("AM","M","RM","R","DS"))) 

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
roa.spric.points

#### bar plot
roa.spric.clus <- roa.com.all[-1] %>%
  inner_join(read.csv("2.Ch3.R.analysis/2.file.species.depth.csv"),by="species") %>%
  mutate(cluster=case_when(dband<40 ~ "altiphotic",
                           dband<100 ~ "upper mesophotic",
                           dband<150 ~ "lower mesophotic",
                           dband<210~ "upper rariphotic",
                           dband<300~ "lower rariphotic",
                           TRUE ~ "deep sea")) %>%
  group_by(cluster,species,distrib.roa) %>%
  summarize(abundance=sum(abu.corr))%>%
  group_by(cluster,distrib.roa) %>%
  summarize(spric=n())%>%
  mutate(distrib.roa=factor(distrib.roa,levels=c("AM","M","RM","R","DS")))%>%
  mutate(cluster=factor(cluster,levels=rev(c("altiphotic","upper mesophotic","lower mesophotic",
                                             "upper rariphotic","lower rariphotic","deep sea"))))

##bar plot
pal <- c("AM"="#d0dc36","M"="#50ABE7","RM"="#BC90DB","R"="#831EB6","DS"="#12086F")

roa.spric.barplot <- ggplot(data=roa.spric.clus, aes(x=cluster,y=spric,fill=distrib.roa))+
  geom_bar(position="stack",stat="identity",width=0.4) +
  scale_fill_manual(values=pal,
                     labels=c("alti/mesophotic","mesophotic","meso/rariphotic",
                              "rariphotic","deep sea"),name="predominant depth")+
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
roa.spric.barplot

#### 1.7. Diversity ####
#### 1.7.1. Species richness ####
## data prep
fish.roatan <- read.csv("2.Ch3.R.analysis/2.file.fish.carib.grouped.csv") %>%
  filter(location == "Roatan") %>%
  select(species,dband,abu.corr)

hill.roatan <- fish.roatan %>%
  spread(species, abu.corr, fill =  0) 


### taxonomic diversity 
taxa.hill0 <- hill_taxa(hill.roatan,q=0)
taxa.hill1 <- hill_taxa(hill.roatan,q=1)
taxa.hill2 <- hill_taxa(hill.roatan,q=2)


div.roatan <- data.frame(depth=sort(unique(fish.roatan$dband)), richness=taxa.hill0, evenness=taxa.hill1,
                         dominance=taxa.hill2)

### Plots with prediction and CI ###
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
  geom_vline(xintercept = 90, lty = 2, color = "grey69") +
  geom_vline(xintercept = 140, lty = 2, color = "grey69") +
  geom_vline(xintercept = 200, lty = 2, color = "grey69") +
  geom_vline(xintercept = 300, lty = 2, color = "grey69") +
  ylab("richness") +
  xlab("depth (m)") 
  #coord_flip()+
  #scale_x_reverse()


#### 1.7.2. beta-diversity ####
betadiv.roatan.cluster <- read.csv("2.Ch3.R.analysis/2.file.fish.carib.grouped.csv") %>%
  filter(location == "Roatan") %>%
  select(species,dband,abu.corr)  %>%
  mutate(cluster=case_when(dband<40 ~ "altiphotic",
                            dband<100 ~ "upper mesophotic",
                           dband<150 ~ "lower mesophotic",
                           dband<210~ "upper rariphotic",
                           dband<300~ "lower rariphotic",
                           TRUE~"deep sea")) %>%
  mutate(cluster=factor(cluster,levels=c("altiphotic","upper mesophotic","lower mesophotic",
                                         "upper rariphotic","lower rariphotic","deep sea"))) %>%
  group_by(cluster,species) %>%
  summarize(abu.corr=sum(abu.corr)) %>%
  ungroup()%>% #actions will no longer be performed grouped by clusters
  spread(species, abu.corr, fill =  0) %>%
  select(-cluster) %>%
  mutate_all(~ case_when(.>0 ~ 1, TRUE~0))

beta.cluster.roatan <- beta.pair(betadiv.roatan.cluster,index.family="sorensen")

beta.cluster.roatan.graph <- data.frame(cluster=c("upper mesophotic","lower mesophotic","upper rariphotic","lower rariphotic","deep sea",
                                                  "lower mesophotic ","upper rariphotic ","lower rariphotic ", "deep sea ",
                                                  "upper rariphotic  ","lower rariphotic  ", "deep sea  ",
                                                  "lower rariphotic   ", "deep sea   ",
                                                  "deep sea    "),
                                        type=rep(c("beta","nestedness","turnover"),each=15),
                                        betadiv=c(beta.cluster.roatan$beta.sor,
                                                  beta.cluster.roatan$beta.sne,
                                                  beta.cluster.roatan$beta.sim)) %>%
  mutate(cluster=factor(cluster,levels=c("upper mesophotic","lower mesophotic","upper rariphotic","lower rariphotic","deep sea",
                                        "lower mesophotic ","upper rariphotic ","lower rariphotic ", "deep sea ",
                                        "upper rariphotic  ","lower rariphotic  ", "deep sea  ",
                                        "lower rariphotic   ", "deep sea   ",
                                        "deep sea    ")))

#write.csv(beta.cluster.roatan.graph,"2.Ch3.R.analysis/betadiv.roatan.fullrange.csv")

plot.betadiv.roa <- ggplot(beta.cluster.roatan.graph)+
  geom_point(aes(x=cluster,y=betadiv,shape=type,color=type),position=position_dodge(width=0.3),
             size=5)+
  scale_shape_manual(values=c(16,22,24))+
  scale_color_viridis_d(begin=0.9,end=0)+
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
  geom_vline(xintercept=c(5.5,9.5,12.5,14.5),linetype="dashed",color="grey")+
  scale_x_discrete(labels=str_wrap(beta.cluster.roatan.graph$cluster,6))+
  ylab("dissimilarity")+
  annotate("text",x=c(3,7.5,11,13.5,15),y=1.1,
           label=str_wrap(c("vs. altiphotic","vs. upper mesophotic","vs. lower mesophotic",
                            "vs. upper rariphotic","vs. lower rariphotic"),10),size=4)+
  annotate("text",x=0,y=1.2,label="Roatan")+
  coord_cartesian(ylim=c(0,1),clip ="off")
plot.betadiv.roa



#### 2. Trans-site dendrogram ####
carib.all <- read.csv("2.Ch3.R.analysis/2.file.fish.binned.csv") %>%
  select(species, dband, location, abundance, abu.corr) %>%
  mutate(dband.grouped = case_when(location=="Roatan" & dband >= 310 & dband <= 330 ~ 310,
                                   location=="Roatan" & dband >= 350 & dband <= 360 ~ 350,
                                   location=="Roatan" & dband >= 390 & dband <= 400 ~ 390,
                                   location=="Roatan" & dband >= 410 & dband <= 430 ~ 410,
                                   location=="Roatan" & dband >= 440 & dband <= 460 ~ 440,
                                   location=="Roatan" & dband >= 470 & dband <= 480 ~ 470,
                                   location=="St. Eustatius" & dband >= 240 & dband <= 270 ~ 240,
                                   location=="Bonaire" & dband >= 190 & dband <= 200 ~ 190,
                                   location=="Bonaire" & dband >= 220 & dband <= 250 ~ 220,
                                   location=="Bonaire" & dband >= 280 & dband <= 300 ~ 280,
                                   TRUE ~ as.numeric(dband))) %>%
  group_by(species, dband.grouped, location) %>%
  summarize(abundance = round(sqrt(sum(abu.corr)))) %>%
  # pivots the table from having variable as columns to species as columns
  spread(species, abundance, fill =  0)

carib.all <- as.data.frame(carib.all)
rownames(carib.all) <- interaction(carib.all$dband.grouped, carib.all$location)
carib.com.all <- carib.all[-c(1:2)] # removes location and depth columns

## SIMPROF
#ward.d clustering
carib.clust.all <- simprof(data = carib.com.all, method.cluster = "ward.D", method.distance="braycurtis",
                           method.transform="identity", alpha=0.0000001,
                           sample.orientation="row", const=0,
                           silent=TRUE, increment=100,
                           undef.zero=TRUE, warn.braycurtis=TRUE)
#saveRDS(carib.clust2,"2.Ch3.R.analysis/carib.clust.wardD.RData")

## "vegdist" calculates ecological similarity and diversity analyses based on matrices
## "hclust" takes a distance matrix and performs agglomerative hierarchical clustering
#carib.clust2 <- readRDS("2.Ch3.R.analysis/carib.clust.wardD.RData")
carib.clust.veg.all <- hclust(vegdist(carib.com.all, method = "bray"), method = "ward.D")

carib.dend.all <- as.dendrogram(carib.clust.veg.all)
carib.dend.data.all <- dendro_data(carib.dend.all, type = "rectangle")

carib.deep.com.h <- carib.all %>%
  unite(label, dband.grouped, location, remove = F, sep = ".")

carib.dend.helper <- carib.dend.data.all$labels %>%
  inner_join(carib.deep.com.h) 

carib.dend.data.all$labels <- carib.dend.helper %>%
  mutate(label = case_when(location=="Roatan" & dband.grouped >= 310 & dband.grouped <= 330 ~ "310-330",
                           location=="Roatan" & dband.grouped >= 350 & dband.grouped <= 360 ~ "350-360",
                           location=="Roatan" & dband.grouped >= 390 & dband.grouped <= 400 ~ "390-400",
                           location=="Roatan" & dband.grouped >= 410 & dband.grouped <= 430 ~ "410-430",
                           location=="Roatan" & dband.grouped >= 440 & dband.grouped <= 460 ~ "440-460",
                           location=="Roatan" & dband.grouped >= 470 & dband.grouped <= 480 ~ "470-480",
                           location=="St. Eustatius" & dband.grouped >= 240 & dband.grouped <= 270 ~ "240-270",
                           location=="Bonaire" & dband.grouped >= 190 & dband.grouped <= 200 ~ "190-200",
                           location=="Bonaire" & dband.grouped >= 220 & dband.grouped <= 250 ~ "220-250",
                           location=="Bonaire" & dband.grouped >= 280 & dband.grouped <= 300 ~ "280-300",
                           TRUE ~ as.character(dband.grouped))) %>%
  mutate(location=factor(location,levels=c("Curacao","Bonaire","St. Eustatius","Roatan")))

new.pal=c("#ffa600","#ef5675","#3DAC78","#003f5c")
carib.dend.plot <- ggplot(carib.dend.data.all$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), color = "grey46")+
  geom_text(data = carib.dend.data.all$labels, aes(x, y, label = label, color = location),
            hjust = 0,  size = 2.5, fontface = "bold") +
  #geom_label(data = carib.dend.data2$labels, aes(x, y, label = label, fill = location),
  #          hjust = 0,  size = 2.5, fontface = "bold")+
  scale_color_manual(values=new.pal)+
  theme_bw() +
  #scale_x_continuous(labels=c(1:96),breaks=c(1:96))+
  theme(legend.position = "top",
        line = element_blank(),
        title = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(0,0,0,1),"lines"))+
  geom_segment(aes(x=0.8,xend=5.2,y=-1.5,yend=-1.5),size=1)+
  geom_segment(aes(x=5.8,xend=12.2,y=-1.5,yend=-1.5),size=1)+
  geom_segment(aes(x=12.8,xend=25.2,y=-1.5,yend=-1.5),size=1)+
  geom_segment(aes(x=25.8,xend=38.2,y=-1.5,yend=-1.5),size=1)+
  geom_segment(aes(x=38.8,xend=48.2,y=-1.5,yend=-1.5),size=1)+
  geom_segment(aes(x=48.8,xend=55.2,y=-1.5,yend=-1.5),size=1)+
  geom_segment(aes(x=55.8,xend=58.2,y=-1.5,yend=-1.5),size=1)+
  geom_segment(aes(x=58.8,xend=63.2,y=-1.5,yend=-1.5),size=1)+
  geom_segment(aes(x=63.8,xend=69.2,y=-1.5,yend=-1.5),size=1)+
  geom_segment(aes(x=69.8,xend=76.2,y=-1.5,yend=-1.5),size=1)+
  geom_segment(aes(x=76.8,xend=80.2,y=-1.5,yend=-1.5),size=1)+
  geom_segment(aes(x=80.8,xend=84.2,y=-1.5,yend=-1.5),size=1)+
  geom_segment(aes(x=84.8,xend=89.2,y=-1.5,yend=-1.5),size=1)+
  geom_segment(aes(x=89.8,xend=98.2,y=-1.5,yend=-1.5),size=1)+
  geom_segment(aes(x=98.8,xend=105.2,y=-1.5,yend=-1.5),size=1)+
  #ylim(-3,10)+ 
  coord_flip()+
  scale_y_reverse(limits=c(10,-3))
carib.dend.plot

