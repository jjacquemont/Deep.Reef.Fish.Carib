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

#### 1. Size structure and community composition - PCA ####
carib.deep.traits <- read.csv(file = "2.Ch3.R.analysis/2.carib.deep.traits.csv") %>%
  select(-X)
carib.deep.fishbase <- read.csv("2.Ch3.R.analysis/2.file.fish.carib.grouped.csv")
## from one line per species per location to one line per individual per location.
carib.size <- carib.deep.traits %>%
  select(species, species.fb, dband, location, hab, abundance, COEF, abu.corr, maxTL, length, biomass, lwa_m, lwb_m) %>%
  as_tibble() %>%
  # mutate: adds or modifies columns in a data frame (tydiverse)
  # map function applies a given function to each element of a list/vector.
  # rep_len function to create a vector of ones with a length specified by .x
  # creates a new column "obs", where each value in the abundance column is replaced by a vector of ones with a length specified by that value.
  mutate(obs = map(abundance, ~rep_len(1, .x))) %>% 
  unnest(cols = c(obs)) %>% # expands the obs column, creating separate rows for each element in the lists in the obs column while duplicating the other columns as necessary.
  select(-abundance, -obs)

### community composition
carib.deep.com <- carib.deep.fishbase %>%
  select(species, dband, location, abundance, abu.corr) %>%
  filter(dband <= 300 & dband >= 40) %>%
  group_by(species, dband, location) %>%
  summarize(abundance = round(sqrt(sqrt(sum(abu.corr))))) %>%
  # pivots the table from having variable as columns to species as columns
  spread(species, abundance, fill =  0)

#add.data <- data.frame(c("Bonaire",240,rep(0,212), "St. Eustatius",
#                                          250, rep(0,212)),rows=2)
#colnames(add.data) <- colnames(carib.deep.com)
#carib.deep.com <- rbind(carib.deep.com, add.data)

##principal component analysis
        # center: a logical value indicating whether the variables should be shifted to be zero centered.
        # scale: a logical value indicating whether the variables should be scaled to have unit variance before the analysis takes place.
carib.deep.pca <- prcomp(carib.deep.com[-c(1:2)], center = TRUE, scale = TRUE)
summary(carib.deep.pca)

carib.deep.points <- 
  # first convert the pca results to a tibble
  as_tibble(carib.deep.pca$x) %>% 
  # add fish species column
  bind_cols(carib.deep.com) %>%
  mutate(zone = case_when(dband >= 40 & dband < 80 ~ "Upper Mesophotic",
                          dband >= 80 & dband < 130 ~ "Lower Mesophotic",
                          dband >= 130 & dband < 190 ~ "Upper Rariphotic",
                          dband >= 190 & dband < 310 ~ "Lower Rariphotic",
                          TRUE ~ "Abyss"))

## taking a set of points defined by PC1 and PC2, calculating the convex hull of those points, 
## and then returning a subset of those points using the indices computed by chull. 
#The result is a selection of points from the original data that form the convex hull.
carib.deep.hull <- 
  carib.deep.points %>% 
  group_by(location, zone) %>% 
  # chull (PC1, PC2) returns a vector of indices representing the points that form the convex hull.
  #  slice function subsets a data frame or a tibble based on row indices.
  slice(chull(PC1, PC2))

carib.deep.load.PC1pos <- 
  as_tibble(carib.deep.pca$rotation, rownames = 'variable') %>%
  top_n(3, PC1) %>%
  rename(species = "variable") %>%
  select(species, PC1, PC2)

carib.deep.load.PC1 <- 
  as_tibble(carib.deep.pca$rotation, rownames = 'variable') %>%
  # select the 3 lines with the highest values in PC1
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



#### 2. Trans-site dendogram ####
carib.deep.com <- carib.deep.fishbase %>%
  select(species, dband, location, abundance, abu.corr) %>%
  filter(dband <= 300 & dband >= 40) %>%
  mutate(dband.grouped = case_when(location=="St. Eustatius" & dband >= 240 & dband <= 270 ~ 240,
                                   location=="Bonaire" & dband >= 190 & dband <= 200 ~ 190,
                                   location=="Bonaire" & dband >= 220 & dband <= 250 ~ 220,
                                   location=="Bonaire" & dband >= 280 & dband <= 300 ~ 280,
                                   TRUE ~ as.numeric(dband))) %>%
  group_by(species, dband.grouped, location) %>%
  summarize(abundance = round(sqrt(sum(abu.corr)))) %>%
  # pivots the table from having variable as columns to species as columns
  spread(species, abundance, fill =  0)

carib.deep.com <- as.data.frame(carib.deep.com)
rownames(carib.deep.com) <- interaction(carib.deep.com$dband.grouped, carib.deep.com$location)
carib.com.cl <- carib.deep.com[-c(1:2)] # removes location and depth columns

## SIMPROF
#determining the number of significant clusters produced using hclust() with the assumption of no a priori groups.
carib.clust <- simprof(data = carib.com.cl, method.cluster = "complete", method.distance="braycurtis",
                     method.transform="identity", alpha=0.0000001,
                     sample.orientation="row", const=0,
                     silent=TRUE, increment=100,
                     undef.zero=TRUE, warn.braycurtis=TRUE)
#ward.d clustering
carib.clust2 <- simprof(data = carib.com.cl, method.cluster = "ward.D", method.distance="braycurtis",
                       method.transform="identity", alpha=0.0000001,
                       sample.orientation="row", const=0,
                       silent=TRUE, increment=100,
                       undef.zero=TRUE, warn.braycurtis=TRUE)
#saveRDS(carib.clust2,"2.Ch3.R.analysis/carib.clust.wardD.RData")

## "vegdist" calculates ecological similarity and diversity analyses based on matrices
## "hclust" takes a distance matrix and performs agglomerative hierarchical clustering
carib.clust2 <- readRDS("2.Ch3.R.analysis/carib.clust.wardD.RData")
carib.clust.veg <- hclust(vegdist(carib.com.cl, method = "bray"), method = "complete")
carib.clust.veg2 <- hclust(vegdist(carib.com.cl, method = "bray"), method = "ward.D")

carib.dend <- as.dendrogram(carib.clust.veg)
carib.dend2 <- as.dendrogram(carib.clust.veg2)
carib.dend.data <- dendro_data(carib.dend, type = "rectangle")
carib.dend.data2 <- dendro_data(carib.dend2, type = "rectangle")

carib.deep.com.h <- carib.deep.com %>%
  unite(label, dband.grouped, location, remove = F, sep = ".")

carib.dend.helper <- carib.dend.data$labels %>%
  inner_join(carib.deep.com.h) 
carib.dend.helper2 <- carib.dend.data2$labels %>%
  inner_join(carib.deep.com.h) 

carib.dend.data$labels <- carib.dend.helper %>%
  mutate(label = case_when(location=="St. Eustatius" & dband.grouped >= 240 & dband.grouped <= 270 ~ "240-270",
                                   location=="Bonaire" & dband.grouped >= 190 & dband.grouped <= 200 ~ "190-200",
                                   location=="Bonaire" & dband.grouped >= 220 & dband.grouped <= 250 ~ "220-250",
                                   location=="Bonaire" & dband.grouped >= 280 & dband.grouped <= 300 ~ "280-300",
                                   TRUE ~ as.character(dband.grouped))) %>%
  mutate(location=factor(location,levels=c("Curacao","Bonaire","St. Eustatius","Roatan")))

carib.dend.data2$labels <- carib.dend.helper2 %>%
  mutate(label = case_when(location=="St. Eustatius" & dband.grouped >= 240 & dband.grouped <= 270 ~ "240-270",
                           location=="Bonaire" & dband.grouped >= 190 & dband.grouped <= 200 ~ "190-200",
                           location=="Bonaire" & dband.grouped >= 220 & dband.grouped <= 250 ~ "220-250",
                           location=="Bonaire" & dband.grouped >= 280 & dband.grouped <= 300 ~ "280-300",
                           TRUE ~ as.character(dband.grouped))) %>%
  mutate(location=factor(location,levels=c("Curacao","Bonaire","St. Eustatius","Roatan")))

carib.dend.plot <- ggplot(carib.dend.data2$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), color = "grey46")+
  geom_text(data = carib.dend.data2$labels, aes(x, y, label = label, color = location),
            hjust = 0,  size = 2.5, fontface = "bold") +
  #geom_label(data = carib.dend.data2$labels, aes(x, y, label = label, fill = location),
   #          hjust = 0,  size = 2.5, fontface = "bold")+
  scale_color_fish(option = "Bodianus_rufus", discrete = T) +
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
  geom_segment(aes(x=5.8,xend=14.2,y=-1.5,yend=-1.5),size=1)+
  geom_segment(aes(x=14.8,xend=21.2,y=-1.5,yend=-1.5),size=1)+
  geom_segment(aes(x=21.8,xend=28.2,y=-1.5,yend=-1.5),size=1)+
  geom_segment(aes(x=28.8,xend=31.2,y=-1.5,yend=-1.5),size=1)+
  geom_segment(aes(x=31.8,xend=34.2,y=-1.5,yend=-1.5),size=1)+
  geom_segment(aes(x=34.8,xend=39.2,y=-1.5,yend=-1.5),size=1)+
  geom_segment(aes(x=39.8,xend=46.2,y=-1.5,yend=-1.5),size=1)+
  geom_segment(aes(x=46.8,xend=50.2,y=-1.5,yend=-1.5),size=1)+
  geom_segment(aes(x=50.8,xend=54.2,y=-1.5,yend=-1.5),size=1)+
  geom_segment(aes(x=54.8,xend=57.2,y=-1.5,yend=-1.5),size=1)+
  geom_segment(aes(x=57.8,xend=71.2,y=-1.5,yend=-1.5),size=1)+
  geom_segment(aes(x=71.8,xend=76.2,y=-1.5,yend=-1.5),size=1)+
  geom_segment(aes(x=76.8,xend=83.2,y=-1.5,yend=-1.5),size=1)+
  geom_segment(aes(x=83.8,xend=96.2,y=-1.5,yend=-1.5),size=1)+
  #ylim(-3,10)+ 
  coord_flip()+
  scale_y_reverse(limits=c(10,-3))
carib.dend.plot


#### 3. PCA and SIMPROF by location ####
#### 3.1. Curacao ####

cur.com <- read.csv("2.Ch3.R.analysis/2.file.fish.carib.grouped.csv") %>%
  filter(location == "Curacao") %>%
  select(species, dband, location, abundance, abu.corr) %>%
  group_by(species, dband) %>%
  summarize(abundance = sqrt(sum(abu.corr))) %>%
  spread(species, abundance, fill =  0) %>%
  as.data.frame()

rownames(cur.com) <- cur.com$dband
cur.com.cl <- cur.com[-1]
#write.csv(cur.com.cl,"file_Curacao_community_cluster.csv")

#### 3.1.0. SIMPROF - saved as file ####
  #determining the number of significant clusters produced using hclust() with the assumption of no a priori groups.
cur.clust <- simprof(data = cur.com.cl, method.cluster = "complete", method.distance="braycurtis",
                     method.transform="identity", alpha=0.0000001,
                     sample.orientation="row", const=0,
                     silent=TRUE, increment=100,
                     undef.zero=TRUE, warn.braycurtis=TRUE)

cur.clust.ward <- simprof(data = cur.com.cl, method.cluster = "ward.D", method.distance="czekanowski",
                     method.transform="identity", alpha=0.0000001,
                     sample.orientation="row", const=0,
                     silent=TRUE, increment=100,
                     undef.zero=TRUE, warn.braycurtis=TRUE)
#saveRDS(cur.clust,"2.Ch3.R.analysis/2.cur.clust.RData")


#### 3.1.1 dendogram ####
# Curacao: 7 groups, from shallow to deep
cur.clust <- readRDS("2.Ch3.R.analysis/2.cur.clust.RData")
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
  gather(key = "cluster", value = "dband") %>% # transform wide to long format 
  distinct() %>% # table with cluster and associated dbands
  mutate(dband = as.numeric(dband)) %>%
  # full_join(cur.com) %>%
  #  select(cluster, dband) %>%
  mutate(label = as.factor(dband)) %>%
  arrange(dband) %>%
  mutate(run_num = 1:27) %>%
  group_by(cluster) %>%
  mutate(ord_num = min(run_num))

# vegdist() : distance between each depth bin given relative abundance of species
# hlcust(): hierarchical clustering. 
  #takes a dissimilarity table as produced by vegdist()
  #returns the order of observations by similarity ($order), and the order of merges ($merge)
cur.clust.veg <- hclust(vegdist(cur.com.cl, method = "bray"), method = "complete")
cur.clust.veg.ward <- hclust(vegdist(cur.com.cl, method = "bray"), method = "ward.D")

# changes class of object from hclust to dendogram 
cur.dend <- as.dendrogram(cur.clust.veg) 

# length and coordinates of each branch in the dendogram ($segments)
cur.dend.data <- dendro_data(cur.dend, type = "rectangle")

# joining with values of dbands
cur.dend.helper <- cur.dend.data$labels %>%
  inner_join(cur.res.tib) 
cur.dend.data$labels <- cur.dend.helper

# get colors for each level using fishualize and tinter package
carib.pal <- fish(4, option = "Bodianus_rufus")
cur.pal <- tinter(carib.pal[1], steps = 5)
blue_palette <- c("#BAE0F3","#87CEEB","#6CBDE9","#50ABE7", "#4895EF", "#4361EE", "#2835AF", "#12086F")
mesoP_palette <- c("#DCEEF3","#C2E2EA","#A7D5E1","#8DC8D8","#72BBCE")
rariP_palette <- c("#DDC5EB","#BC90DB","#831EB6")
pal <- c(blue_palette[3],blue_palette[3],rariP_palette[1],rariP_palette[2])

pdf("2.Ch3.R.analysis/graphs/cur.dendro.pdf")
cur.plot <- ggplot(cur.dend.data$segments) + 
  geom_rect(xmin=22.7,xmax=27.3,ymin=-0.13, ymax=0.95,fill=blue_palette[3])+
  geom_rect(xmin=18.7,xmax=22.3,ymin=-0.13, ymax=0.95,fill=blue_palette[1])+
  geom_rect(xmin=12.7,xmax=18.3,ymin=-0.13, ymax=0.95,fill=rariP_palette[1])+
  geom_rect(xmin=0.7,xmax=12.3,ymin=-0.13, ymax=0.95,fill=rariP_palette[2])+
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), color = "black")+
  geom_text(data = cur.dend.data$labels, aes(x, y, label = label), color = "black", #below for yellow gradient
                                             #color = as.factor(ord_num)), 
            hjust = 1, angle = 90, size = 3, fontface = "bold") +
  scale_color_manual(values = cur.pal[-1]) +
  #theme of plot
  theme(legend.position = "none",
        line = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  # clusters highlight
  geom_segment(aes(x=0.8,xend=3.2,y=-0.13,yend=-0.13),size=1)+
  geom_segment(aes(x=3.8,xend=8.2,y=-0.13,yend=-0.13),size=1)+
  geom_segment(aes(x=8.8,xend=12.2,y=-0.13,yend=-0.13),size=1)+
  geom_segment(aes(x=12.8,xend=18.2,y=-0.13,yend=-0.13),size=1)+
  geom_segment(aes(x=18.8,xend=22.2,y=-0.13,yend=-0.13),size=1)+
  geom_segment(aes(x=22.8,xend=24.2,y=-0.13,yend=-0.13),size=1)+
  geom_segment(aes(x=24.8,xend=27.2,y=-0.13,yend=-0.13),size=1)+
  ggtitle("Curaçao")+
  ylab("")+
  ylim(-0.2,1)
cur.plot

dev.off()


#### 3.1.2. PCA plot ####
cur.pca <- prcomp(cur.com.cl, center = TRUE, scale = TRUE)
summary(cur.pca)

## Species most contributing to PC
cur.pca.speciesPC1 <- as.tibble(cur.pca$rotation[,c(1:4)]) %>% # PC1-PC4
  mutate(species=rownames(cur.pca$rotation[,c(1:4)]), location="Curacao",PC="PC1") %>%
  arrange(desc(abs(PC1))) %>%
  slice_head(n=5)
cur.pca.speciesPC1

cur.pca.speciesPC2 <- as.tibble(cur.pca$rotation[,c(1:4)]) %>% # PC1-PC4
  mutate(species=rownames(cur.pca$rotation[,c(1:4)]), location="Curacao",PC="PC2") %>%
  arrange(desc(abs(PC2))) %>%
  slice_head(n=5)
cur.pca.speciesPC2

cur.pca.species <- distinct(bind_rows(cur.pca.speciesPC1,cur.pca.speciesPC2))
  
## linking the clusters from h. clustering with coordinates in PCA
cur.points <- 
  as_tibble(cur.pca$x) %>% # coordinates of the depth band points in the 27 projected dimensions
  bind_cols(cur.com) %>%
  inner_join(cur.dend.helper)%>%
  mutate(cat=case_when(dband<=70 ~ "upper mesophotic",
                                    dband<=120~ "lower mesophotic",
                                    dband<=180~ "upper rariphotic",
                                    dband>180~ "lower rariphotic"))%>%
  mutate(cat=factor(cat,levels=c("upper mesophotic","lower mesophotic","upper rariphotic","lower rariphotic")))



### SIMPER: contribution of species to the pair-wise differences between clusters 
cur.simper <- simper(cur.com.cl, cur.points$cluster) # length is 21 (nb of different pairs among 7 clusters)
cur.simsum <- summary(cur.simper)

# SIMPER but for depth zones and not cluster 
cur.simper.zones <- simper(cur.com.cl, cur.points$cat) 
cur.simsum.zones <- summary(cur.simper.zones)

## obtaining species name that contribute to > 0.05 dissimilarity between clusters 
  # and their coordinates 
cur.simspecies <- data_frame("species" = as.character())
for (i in 1:length(cur.simsum)){  #length is 21, nb of pairs among 7 clusters
  spec <- as.data.frame(cur.simsum[i]) %>%
    mutate(species = rownames(.)) %>%
    select(1,7,8) %>% #"average" (contribution to dissimilarity) and "species"
    filter(.[[2]] < 0.05 & .[[1]] > 0.05) #pvalue<0.05, avg contribution>0.05
  cur.simspecies <- bind_rows(cur.simspecies, spec)
}
# removing duplicates (species contributing to distinctiveness in several pairs)
cur.species <- cur.simspecies %>%
  distinct(species)

## same, but for depth zones instead of clusters
cur.simspecies.zones <- data_frame("species" = as.character())
for (i in 1:length(cur.simsum.zones)){  #length is 21, nb of pairs among 7 clusters
  spec <- as.data.frame(cur.simsum.zones[i]) %>%
    mutate(species = rownames(.)) %>%
    select(1,7,8) %>% #"average" (contribution to dissimilarity) and "species"
    filter(.[[2]] < 0.05 & .[[1]] > 0.055) #pvalue<0.05, avg contribution>0.05
  cur.simspecies.zones <- bind_rows(cur.simspecies.zones, spec)
}
cur.simspecies.zones$location="Curaçao"

## obtaining coordinates of species (SIMPER)
cur.loads <- 
  # $rotation contains the coordinates of each fish species in the 27 projected dimensions
  as_tibble(cur.pca$rotation, rownames = 'species') %>%  
  inner_join(cur.species) %>% # selects only species contributing to dissimilarity
  select(species, PC1, PC2) # selects only coordinates in the 2 main directions

## obtaining coordinates of species (PC drivers)
cur.loads.PC <- 
  # $rotation contains the coordinates of each fish species in the 27 projected dimensions
  as_tibble(cur.pca$rotation, rownames = 'species') %>%  
  inner_join(cur.pca.species) %>% # selects only species contributing to dissimilarity
  select(species, PC1, PC2)

### coordinates of depth bands forming the convex hull of each cluster in the PC1, PC2 dimensions
  # used to make the polygons representing each cluster in the PCA plot
  # represents all depth bands except depth bands "contained" in a polygon cluster
cur.hull <- cur.points %>% 
  group_by(cluster) %>% # using cluster from hierarchial clustering hclust()
  # chull: coordinates of points which constitute the outside border (=convex hull) of a cluster of points
  slice(chull(PC1, PC2)) %>%
  mutate(cat=case_when(dband<=80 ~ "upper mesophotic",
                       dband<=120~ "lower mesophotic",
                       dband<=180~ "upper rariphotic",
                       dband>180~ "lower rariphotic")) %>%
  mutate(cat=factor(cat,levels=c("upper mesophotic","lower mesophotic","upper rariphotic","lower rariphotic")))


#### plot
pal <- c(blue_palette[1],blue_palette[3],rariP_palette[1],rariP_palette[2])

pdf("2.Ch3.R.analysis/graphs/cur.pca.pdf")
cur.pca.plot <- ggplot(cur.points) +
  geom_jitter(data = cur.points, aes(x = PC1, y = PC2, 
              fill = cat), shape = 21, alpha = 0.5, color = "black") +
  # polygon shape of each 7 cluster
  geom_polygon(data = cur.hull,
               aes(x = PC1, y = PC2, fill = cat,
                   colour = cat, group = as.factor(ord_num)),
               alpha = 0.3,  show.legend = TRUE) +
  # significant species from SIMPER analyses
  geom_segment(data = cur.loads, aes(x = 0, y = 0, xend = PC1*30, yend = PC2*30), 
               lwd = 0.1, 
               arrow = arrow(length = unit(1/2, 'picas'), type = "open")) +
  geom_text_repel(data = cur.loads[c(1:6),], aes(x = PC1*30,  y = PC2*30, label = species), 
                  size = 4, segment.color = "grey69")+
  theme_classic() +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) 
 # guides(colour = guide_legend(nrow = 1),
  #       fill = guide_legend(nrow = 1))
cur.pca.plot
dev.off()

## MDS trial
mds_cur <- cmdscale(vegdist(cur.com.cl, method = "bray"))
plot(mds_cur,col=cur.hull$cat, pch = 19, xlab = "MDS1", ylab = "MDS2")

#### 3.2. Bonaire ####
bon.com <- read.csv("2.Ch3.R.analysis/2.file.fish.carib.grouped.csv") %>%
  filter(location == "Bonaire") %>%
  select(species, dband, location, abundance, abu.corr) %>%
  #mutate(dband.grouped = case_when(dband >= 190 & dband <= 200 ~ 190,
   #                                dband >= 220 & dband <= 250 ~ 220,
    #                               dband >= 280 & dband <= 300 ~ 280,
     #                              TRUE ~ as.numeric(dband))) %>% 
  group_by(species, dband) %>% 
  summarize(abundance = sqrt(sum(abu.corr))) %>%
  spread(species, abundance, fill =  0) %>%
  as.data.frame()


rownames(bon.com) <- bon.com$dband
bon.com.cl <- bon.com[-1]
#write.csv(bon.com.cl,"file_bon_com_cl.csv")

bon.clust <- simprof(data = bon.com.cl, method.cluster = "complete", method.distance="braycurtis",
                     method.transform="identity", alpha=0.0000001,
                     sample.orientation="row", const=0,
                     silent=TRUE, increment=100,
                     undef.zero=TRUE, warn.braycurtis=TRUE)
bon.clust
# Bonaire: 8 groups, from shallow to deep
bon.res.tib <- tibble(cl1 = list(bon.clust$significantclusters[[1]]),
                      cl2 = list(bon.clust$significantclusters[[2]]),
                      cl3 = list(bon.clust$significantclusters[[3]]),
                      cl4 = list(bon.clust$significantclusters[[4]]),
                      cl5 = list(bon.clust$significantclusters[[5]]),
                      cl6 = list(bon.clust$significantclusters[[6]]),
                      cl7 = list(bon.clust$significantclusters[[7]]),
                      cl8 = list(bon.clust$significantclusters[[8]])) %>%
  unnest_longer(cl1) %>%
  unnest_longer(cl2) %>%
  unnest_longer(cl3) %>%
  unnest_longer(cl4) %>%
  unnest_longer(cl5) %>%
  unnest_longer(cl6) %>%
  unnest_longer(cl7) %>%
  unnest_longer(cl8) %>%
  gather(key = "cluster", value = "dband") %>%
  distinct() %>%
  mutate(dband = as.numeric(dband)) %>%
  full_join(bon.com) %>%
  select(cluster, dband) %>%
  mutate(label = as.factor(dband)) %>%
  arrange(dband) %>%
  mutate(run_num = 1:26) %>%
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
carib.pal <- fish(4, option = "Bodianus_rufus")
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

bon.simper <- simper(bon.com[-1], bon.points$cluster)
bon.simsum <- summary(bon.simper)

bon.simspecies <- data_frame("species" = as.character())

for (i in 1:length(bon.simsum)){
    spec <- as.data.frame(cur.simsum[i]) %>%
      mutate(species = rownames(.)) %>%
      select(1,7,8) %>% #"average" (contribution to dissimilarity) and "species"
      filter(.[[2]] < 0.05 & .[[1]] > 0.05) #pvalue<0.05, avg contribution>0.05
  bon.simspecies <- bind_rows(bon.simspecies, spec)
  
}

bon.species <- bon.simspecies %>%
  select(species)%>%
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





#### 3.3. St. Eustatius ####
sta.com <- read.csv("2.Ch3.R.analysis/2.file.fish.carib.grouped.csv")e %>%
  filter(location == "St. Eustatius") %>%
  select(species, dband, location, abundance, abu.corr) %>%
  filter(dband <= 300) %>%
  #mutate(dband.grouped = case_when(dband >= 240 & dband <= 270 ~ 240,
  #                                TRUE ~ dband)) %>%
  group_by(species, dband) %>%
  summarize(abundance = sqrt(sum(abu.corr))) %>%
  spread(species, abundance, fill =  0) %>%
  as.data.frame()

rownames(sta.com) <- sta.com$dband
sta.com.cl <- sta.com[-1]

#write.csv(sta.com.cl,"file_sta.com.cl.csv")

sta.clust <- simprof(data = sta.com.cl, method.cluster = "complete", method.distance="braycurtis",
                     method.transform="identity", alpha=0.0000001,
                     sample.orientation="row", const=0,
                     silent=TRUE, increment=100,
                     undef.zero=TRUE, warn.braycurtis=TRUE)
sta.clust

# 5 groups
sta.res.tib <- tibble(cl1 = list(sta.clust$significantclusters[[1]]),
                      cl2 = list(sta.clust$significantclusters[[2]]),
                      cl3 = list(sta.clust$significantclusters[[3]]),
                      cl4 = list(sta.clust$significantclusters[[4]]),
                      cl5 = list(sta.clust$significantclusters[[5]])) %>%
  unnest_longer(cl1) %>%
  unnest_longer(cl2) %>%
  unnest_longer(cl3) %>%
  unnest_longer(cl4) %>%
  unnest_longer(cl5) %>%
  gather(key = "cluster", value = "dband") %>%
  distinct() %>%
  mutate(dband = as.numeric(dband)) %>%
  full_join(sta.com) %>%
  select(cluster, dband) %>%
  mutate(label = as.numeric(dband)) %>%
  arrange(dband) %>%
  mutate(run_num = 1:23) %>%
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
sta.pal <- tinter(carib.pal[3], steps = 4)


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
    select(1,8) %>%
    filter(.[[1]] > 0.05) %>% 
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

#### 3.4. Roatan ####
roa.com <- read.csv("2.Ch3.R.analysis/2.file.fish.carib.grouped.csv") %>%
  filter(location == "Roatan") %>%
  filter(dband >= 40 & dband <= 300) %>%
  select(species, dband, location, abundance, abu.corr) %>%
  group_by(species, dband) %>%
  summarize(abundance = sqrt(sum(abu.corr))) %>%
  spread(species, abundance, fill =  0) %>%
  as.data.frame()

rownames(roa.com) <- roa.com$dband
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
  gather(key = "cluster", value = "dband") %>%
  distinct() %>%
  mutate(dband = as.numeric(dband))%>%
  full_join(roa.com) %>%
  select(cluster, dband) %>%
  mutate(label = as.numeric(dband)) %>%
  arrange(dband) %>%
  mutate(run_num = 1:27) %>%
  group_by(cluster) %>%
  mutate(ord_num = min(run_num))
#write.csv(roa.res.tib,"2.Ch3.R.analysis/2.roa.res.tib.csv")

# hclust 
roa.res.tib <- read.csv("2.Ch3.R.analysis/2.roa.res.tib.csv")
roa.clust.veg <- hclust(vegdist(roa.com.cl, method = "bray"), method = "complete")
roa.dend <- as.dendrogram(roa.clust.veg)
roa.dend.data <- dendro_data(roa.dend, type = "rectangle")

roa.dend.helper <- roa.dend.data$labels %>%
  mutate(label = as.numeric(label))%>%
  inner_join(roa.res.tib)


roa.dend.data$labels <- roa.dend.helper

# get colors for each level using fishualize and tinter package
carib.pal <- fish(4, option = "Bodianus_rufus")
roa.pal <- tinter(carib.pal[4], steps = 5)

pdf("2.Ch3.R.analysis/graphs/roa.dendro.pdf")
roa.plot <- ggplot(roa.dend.data$segments, horiz = TRUE) + 
  geom_rect(xmin=16.7,xmax=21.3,ymin=-0.13, ymax=0.95,fill=blue_palette[3])+
  geom_rect(xmin=11.7,xmax=16.3,ymin=-0.13, ymax=0.95,fill=rariP_palette[1])+
  geom_rect(xmin=21.7,xmax=27.3,ymin=-0.13, ymax=0.95,fill=blue_palette[1])+
  geom_rect(xmin=0.7,xmax=11.3,ymin=-0.13, ymax=0.95,fill=rariP_palette[2])+
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), color = "black")+
  geom_text(data = roa.dend.data$labels, aes(x, y, label = label),color = "black",
            hjust = 1, angle = 60, size = 3, fontface = "bold") +
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
  geom_segment(aes(x=0.8,xend=11.2,y=-0.13,yend=-0.13),size=1)+
  geom_segment(aes(x=11.8,xend=16.2,y=-0.13,yend=-0.13),size=1)+
  geom_segment(aes(x=16.8,xend=21.2,y=-0.13,yend=-0.13),size=1)+
  geom_segment(aes(x=21.8,xend=24.2,y=-0.13,yend=-0.13),size=1)+
  geom_segment(aes(x=24.8,xend=27.2,y=-0.13,yend=-0.13),size=1)+
  ggtitle("Roatan")+
  ylim(-0.15,1)
roa.plot
dev.off()

### Roatan pca plot
roa.pca <- prcomp(roa.com[-1], center = TRUE, scale = TRUE)
summary(roa.pca)

## Species most contributing to PC
roa.pca.speciesPC1 <- as.tibble(roa.pca$rotation[,c(1:4)]) %>% # PC1-PC4
  mutate(species=rownames(roa.pca$rotation[,c(1:4)]), location="Roatan",PC="PC1") %>%
  arrange(desc(abs(PC1))) %>%
  slice_head(n=5)
roa.pca.speciesPC1

roa.pca.speciesPC2 <- as.tibble(roa.pca$rotation[,c(1:4)]) %>% # PC1-PC4
  mutate(species=rownames(roa.pca$rotation[,c(1:4)]), location="Roatan",PC="PC2") %>%
  arrange(desc(abs(PC2))) %>%
  slice_head(n=5)
roa.pca.speciesPC2

roa.points <- 
  # first convert the pca results to a tibble
  as_tibble(roa.pca$x) %>% 
  bind_cols(roa.com) %>%
  inner_join(roa.dend.helper)  %>%
  mutate(cat=case_when(dband<100 ~ "upper mesophotic",
                       dband<150 ~ "lower mesophotic",
                       dband<210~ "upper rariphotic",
                       dband>=210~ "lower rariphotic")) %>%
  mutate(cat=factor(cat,levels=c("upper mesophotic","lower mesophotic","upper rariphotic","lower rariphotic")))


roa.hull <- 
  roa.points %>% 
  group_by(cluster) %>% 
  slice(chull(PC1, PC2))%>%
  mutate(cat=case_when(dband<100 ~ "upper mesophotic",
                       dband<150 ~ "lower mesophotic",
                       dband<210~ "upper rariphotic",
                       dband>=210~ "lower rariphotic")) %>%
  mutate(cat=factor(cat,levels=c("upper mesophotic","lower mesophotic","upper rariphotic","lower rariphotic")))


roa.simper <- with(roa.points, simper(roa.com[-1], cluster))
roa.simsum <- summary(roa.simper)

roa.simspecies <- data_frame("species" = as.character())

for (i in 1:length(roa.simsum)){
  spec <- as.data.frame(roa.simsum[i]) %>%
    mutate(species = rownames(.)) %>%
    select(1,8) %>%
    filter(.[[1]] > 0.05) %>% 
    select(species)
  
  roa.simspecies <- bind_rows(roa.simspecies, spec)}

roa.species <- roa.simspecies %>%
  distinct()

# SIMPER but for depth zones and not cluster 
roa.simper.zones <- simper(roa.com[-1], roa.points$cat) 
roa.simsum.zones <- summary(roa.simper.zones)

roa.simspecies.zones <- data_frame("species" = as.character())
for (i in 1:length(roa.simsum.zones)){  #length is 21, nb of pairs among 7 clusters
  spec <- as.data.frame(roa.simsum.zones[i]) %>%
    mutate(species = rownames(.)) %>%
    select(1,7,8) %>% #"average" (contribution to dissimilarity) and "species"
    filter(.[[2]] < 0.05 & .[[1]] > 0.05) #pvalue<0.05, avg contribution>0.05
  roa.simspecies.zones <- bind_rows(roa.simspecies.zones, spec)
}
roa.simspecies.zones$location="Roatan"

roa.loads <- 
  as_tibble(roa.pca$rotation, rownames = 'species') %>%
  inner_join(roa.species) %>%
  select(species, PC1, PC2)

pdf("2.Ch3.R.analysis/graphs/roa.pca.pdf")
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
               aes(x = 0, y = 0, xend = PC1*30,yend = PC2*30), 
               lwd = 0.1, arrow = arrow(length = unit(1/2, 'picas'), type = "open")) +
  geom_text_repel(data = roa.loads[c(1:5),], aes(x = PC1*30,y = PC2*30, label = species), 
                  size = 4) +
  theme_classic() +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal)
  #guides(colour = guide_legend(nrow = 1), fill = guide_legend(nrow = 1))
roa.pca.plot
dev.off()

#### 4. PCA and SIMPROF by location - pooled  ####
### pooling depth band < 5 individuals
#### 4.1. Bonaire ####
bon.com <- read.csv("2.Ch3.R.analysis/2.file.fish.carib.grouped.csv") %>%
  filter(location == "Bonaire") %>%
  select(species, dband, location, abundance, abu.corr) %>%
  mutate(dband.grouped = case_when(dband >= 190 & dband <= 200 ~ 190,
                                  dband >= 220 & dband <= 250 ~ 220,
                                 dband >= 280 & dband <= 300 ~ 280,
                                TRUE ~ as.numeric(dband))) %>% 
  group_by(species, dband.grouped) %>% 
  summarize(abundance = sqrt(sum(abu.corr))) %>%
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
#saveRDS(bon.clust,"2.Ch3.R.analysis/2.cur.clust.RData")

# Bonaire: 6 groups, from shallow to deep
bon.clust <- readRDS("2.Ch3.R.analysis/2.cur.clust.RData")
bon.res.tib <- tibble(cl1 = list(bon.clust$significantclusters[[1]]),
                      cl2 = list(bon.clust$significantclusters[[2]]),
                      cl3 = list(bon.clust$significantclusters[[3]]),
                      cl4 = list(bon.clust$significantclusters[[4]]),
                    cl5 = list(bon.clust$significantclusters[[5]]),
                    cl6 = list(bon.clust$significantclusters[[6]])) %>%
  unnest_longer(cl1) %>%
  unnest_longer(cl2) %>%
  unnest_longer(cl3) %>%
  unnest_longer(cl4) %>%
 unnest_longer(cl5) %>%
unnest_longer(cl6) %>%
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


bon.clust.veg <- hclust(vegdist(bon.com.cl, method = "bray"), method = "complete")
bon.dend <- as.dendrogram(bon.clust.veg)
bon.dend.data <- dendro_data(bon.dend, type = "rectangle")

bon.dend.helper <- bon.dend.data$labels %>%
  inner_join(bon.res.tib) 


bon.dend.data$labels <- bon.dend.helper %>%
  mutate(label= case_when(label == "190"  ~ "190-200",
                          label == "220"  ~ "220-250",
                          label == "280" ~ "280-300",
                          TRUE ~ as.character(label)))

# get colors for each level using fishualize and tinter package
blue_palette <- c("#BAE0F3","#87CEEB","#6CBDE9","#50ABE7", "#4895EF", "#4361EE", "#2835AF", "#12086F")
mesoP_palette <- c("#DCEEF3","#C2E2EA","#A7D5E1","#8DC8D8","#72BBCE")
rariP_palette <- c("#DDC5EB","#BC90DB","#831EB6")

pdf("bonaire.dendro.pdf")
bon.plot.pooled <- ggplot(bon.dend.data$segments) + 
  geom_rect(xmin=16.7,xmax=21.3,ymin=-0.2, ymax=0.95,fill=blue_palette[3])+
  geom_rect(xmin=11.7,xmax=16.3,ymin=-0.2, ymax=0.95,fill=rariP_palette[1])+
  geom_rect(xmin=8.7,xmax=11.3,ymin=-0.2, ymax=0.95,fill=blue_palette[1])+
  geom_rect(xmin=0.7,xmax=8.3,ymin=-0.2, ymax=0.95,fill=rariP_palette[2])+
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), color = "black")+
  geom_text(data = bon.dend.data$labels,aes(x, y=-0.01,label = str_wrap(label,3)) , color="black",
                                             #color = as.factor(ord_num)),
            hjust = 1, angle = 60, size = 3, fontface = "bold") +
  theme(legend.position = "none",
        line = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_segment(aes(x=0.8,xend=3.2,y=-0.2,yend=-0.2),size=1)+
  geom_segment(aes(x=3.8,xend=6.2,y=-0.2,yend=-0.2),size=1)+
  geom_segment(aes(x=6.8,xend=8.2,y=-0.2,yend=-0.2),size=1)+
  geom_segment(aes(x=8.8,xend=11.2,y=-0.2,yend=-0.2),size=1)+
  geom_segment(aes(x=11.8,xend=16.2,y=-0.2,yend=-0.2),size=1)+
  geom_segment(aes(x=16.8,xend=21.2,y=-0.2,yend=-0.2),size=1)+
  ggtitle("Bonaire")+
  ylab("")+
  ylim(-0.2,1)
bon.plot.pooled
dev.off()

### Bonaire pca plot
bon.pca <- prcomp(bon.com[-1], center = TRUE, scale = TRUE)
summary(bon.pca)

## Species most contributing to PC
bon.pca.speciesPC1 <- as.tibble(bon.pca$rotation[,c(1:4)]) %>% # PC1-PC4
  mutate(species=rownames(bon.pca$rotation[,c(1:4)]), location="Bonaire",PC="PC1") %>%
  arrange(desc(abs(PC1))) %>%
  slice_head(n=5)
bon.pca.speciesPC1

bon.pca.speciesPC2 <- as.tibble(bon.pca$rotation[,c(1:4)]) %>% # PC1-PC4
  mutate(species=rownames(bon.pca$rotation[,c(1:4)]), location="Bonaire",PC="PC2") %>%
  arrange(desc(abs(PC2))) %>%
  slice_head(n=5)
bon.pca.speciesPC2

## points 
bon.points <- 
  # first convert the pca results to a tibble
  as_tibble(bon.pca$x) %>% 
  bind_cols(bon.com) %>%
  inner_join(bon.dend.helper)%>% 
  mutate(cat=case_when(dband.grouped<70 ~ "upper mesophotic",
                       dband.grouped<120~ "lower mesophotic",
                       dband.grouped<170~ "upper rariphotic",
                       dband.grouped>=170~ "lower rariphotic")) %>%
  mutate(cat=factor(cat,levels=c("upper mesophotic","lower mesophotic","upper rariphotic","lower rariphotic")))


bon.hull <- 
  bon.points %>% 
  group_by(cluster) %>% 
  slice(chull(PC1, PC2))%>% 
 mutate(cat=case_when(dband.grouped<70 ~ "upper mesophotic",
                         dband.grouped<120~ "lower mesophotic",
                         dband.grouped<170~ "upper rariphotic",
                         dband.grouped>=170~ "lower rariphotic")) %>%
  mutate(cat=factor(cat,levels=c("upper mesophotic","lower mesophotic","upper rariphotic","lower rariphotic")))


bon.simper <- simper(bon.com[-1], bon.points$cluster)
bon.simsum <- summary(bon.simper)

bon.simspecies <- data_frame("species" = as.character())

for (i in 1:length(bon.simsum)){
  spec <- as.data.frame(bon.simsum[i]) %>%
    mutate(species = rownames(.)) %>%
    select(1,7,8) %>% #"average" (contribution to dissimilarity) and "species"
    filter(.[[2]] < 0.05 & .[[1]] > 0.07) %>%  #pvalue<0.05, avg contribution>0.05
  select(species)
  bon.simspecies <- bind_rows(bon.simspecies, spec)
  
}

bon.species <- bon.simspecies %>%
  distinct()

# SIMPER but for depth zones and not cluster 
bon.simper.zones <- simper(bon.com[-1], bon.points$cat) 
bon.simsum.zones <- summary(bon.simper.zones)

bon.simspecies.zones <- data_frame("species" = as.character())
for (i in 1:length(bon.simsum.zones)){  #length is 21, nb of pairs among 7 clusters
  spec <- as.data.frame(bon.simsum.zones[i]) %>%
    mutate(species = rownames(.)) %>%
    select(1,7,8) %>% #"average" (contribution to dissimilarity) and "species"
    filter(.[[2]] < 0.05 & .[[1]] > 0.05) #pvalue<0.05, avg contribution>0.05
  bon.simspecies.zones <- bind_rows(bon.simspecies.zones, spec)
}
bon.simspecies.zones$location="Bonaire"

bon.loads <- 
  as_tibble(bon.pca$rotation, rownames = 'species') %>%
  inner_join(bon.species) %>%
  select(species, PC1, PC2)

pdf("2.Ch3.R.analysis/graphs/bonaire.pca.pdf") 
bon.pca.plot <- ggplot(bon.points) +
  geom_jitter(data = bon.points, aes(x = PC1, y = PC2, fill = cat),
              shape = 21, alpha = 0.5, color = "black") +
  geom_polygon(data = bon.hull,
               aes(x = PC1, y = PC2, fill = cat,
                   colour = cat,
                   group = as.factor(ord_num)),
               alpha = 0.3,
               show.legend = TRUE) +
  geom_segment(data = bon.loads[c(1:5),], 
               aes(x = 0, y = 0, 
                   xend = PC1*30,
                   yend = PC2*30), lwd = 0.1, 
               arrow = arrow(length = unit(1/2, 'picas'), type = "open")) +
  geom_text_repel(data = bon.loads[c(1:5),], aes(x = PC1*30,  y = PC2*30, label = species), 
                  size = 4, segment.color = "NA")+
  theme_classic() +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) 
bon.pca.plot
dev.off()

#### 4.2. St. Eustatius ####
sta.com <- read.csv("2.Ch3.R.analysis/2.file.fish.carib.grouped.csv") %>%
  filter(location == "St. Eustatius") %>%
  select(species, dband, location, abundance, abu.corr) %>%
  mutate(dband.grouped = case_when(dband >= 240 & dband <= 270 ~ 240,
                                  TRUE ~ as.numeric(dband))) %>%
  group_by(species, dband.grouped) %>%
  summarize(abundance = sqrt(sum(abu.corr))) %>%
  spread(species, abundance, fill =  0) %>%
  as.data.frame()

rownames(sta.com) <- sta.com$dband.grouped
sta.com.cl <- sta.com[-1]

#### no need to run - file saved ####

sta.clust <- simprof(data = sta.com.cl, method.cluster = "complete", method.distance="braycurtis",
                     method.transform="identity", alpha=0.0000001,
                     sample.orientation="row", const=0,
                     silent=TRUE, increment=100,
                     undef.zero=TRUE, warn.braycurtis=TRUE)
sta.clust

# 3 groups
sta.res.tib <- tibble(cl1 = list(sta.clust$significantclusters[[1]]),
                      cl2 = list(sta.clust$significantclusters[[2]]),
                      cl3 = list(sta.clust$significantclusters[[3]])) %>%
  unnest_longer(cl1) %>%
  unnest_longer(cl2) %>%
  unnest_longer(cl3) %>%
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

#write_csv(sta.res.tib,"2.Ch3.R.analysis/2.sta.res.tib.csv")

#### Run again starting here ####
sta.res.tib <- read.csv("2.Ch3.R.analysis/2.sta.res.tib.csv")

sta.clust.veg <- hclust(vegdist(sta.com.cl, method = "bray"), method = "complete")
sta.dend <- as.dendrogram(sta.clust.veg)
sta.dend.data <- dendro_data(sta.dend, type = "rectangle")

sta.dend.helper <- sta.dend.data$labels %>%
  mutate(label = as.numeric(label)) %>%
  inner_join(sta.res.tib)

sta.dend.data$labels <- sta.dend.helper


# get colors for each level using fishualize and tinter package
sta.dend.data$labels <- sta.dend.helper %>%
  mutate(label= case_when(label == "240"  ~ "240-270",
                          TRUE ~ as.character(label)))

pdf("2.Ch3.R.analysis/graphs/statia.dendro.pdf") 
sta.plot.pooled <- ggplot(sta.dend.data$segments, horiz = TRUE) + 
  geom_rect(xmin=12.7,xmax=17.8,ymin=-0.22, ymax=0.95,fill=blue_palette[3])+
  geom_rect(xmin=17.6,xmax=21.3,ymin=-0.22, ymax=0.95,fill=rariP_palette[1])+
  geom_rect(xmin=6.7,xmax=12.3,ymin=-0.22, ymax=0.95,fill=blue_palette[1])+
  geom_rect(xmin=0.7,xmax=6.3,ymin=-0.22, ymax=0.95,fill=rariP_palette[2])+
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), color = "black")+
  geom_text(data = sta.dend.data$labels, aes(x, y=-0.01, label = label), color="black",
            hjust = 1, angle = 60, size = 3, fontface = "bold") +
  theme(legend.position = "none",
        line = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_segment(aes(x=0.8,xend=6.2,y=-0.22,yend=-0.22),size=1)+
  geom_segment(aes(x=6.8,xend=12.2,y=-0.22,yend=-0.22),size=1)+
  geom_segment(aes(x=12.8,xend=21.2,y=-0.22,yend=-0.22),size=1)+
  ylab("")+
  ylim(-0.3,1)+
  ggtitle("St. Eustatius")
sta.plot.pooled
dev.off()

### Statia pca plot
sta.pca <- prcomp(sta.com[-1], center = TRUE, scale = TRUE)
summary(sta.pca)

## Species most contributing to PC
sta.pca.speciesPC1 <- as.tibble(sta.pca$rotation[,c(1:4)]) %>% # PC1-PC4
  mutate(species=rownames(sta.pca$rotation[,c(1:4)]), location="Statia",PC="PC1") %>%
  arrange(desc(abs(PC1))) %>%
  slice_head(n=5)
sta.pca.speciesPC1

sta.pca.speciesPC2 <- as.tibble(sta.pca$rotation[,c(1:4)]) %>% # PC1-PC4
  mutate(species=rownames(sta.pca$rotation[,c(1:4)]), location="Statia",PC="PC2") %>%
  arrange(desc(abs(PC2))) %>%
  slice_head(n=5)
sta.pca.speciesPC2

## points
sta.points <- as_tibble(sta.pca$x) %>% # convert the pca results to a tibble
  bind_cols(sta.com) %>%
  inner_join(sta.dend.helper)%>%
  mutate(cat=case_when(dband.grouped<100 ~ "upper mesophotic",
                           dband.grouped<140~ "lower mesophotic",
                           dband.grouped<190~ "upper rariphotic",
                           dband.grouped>=190~ "lower rariphotic")) %>%
  mutate(cat=factor(cat,levels=c("upper mesophotic","lower mesophotic","upper rariphotic","lower rariphotic")))


sta.hull <- 
  sta.points %>% 
  group_by(cluster) %>% 
  slice(chull(PC1, PC2))%>%
  mutate(cat=case_when(dband.grouped<100 ~ "upper mesophotic",
                       dband.grouped<140~ "lower mesophotic",
                       dband.grouped<190~ "upper rariphotic",
                       dband.grouped>=190~ "lower rariphotic")) %>%
  mutate(cat=factor(cat,levels=c("upper mesophotic","lower mesophotic","upper rariphotic","lower rariphotic")))


sta.simper <- with(sta.points, simper(sta.com[-1], cluster))
sta.simsum <- summary(sta.simper)

sta.simspecies <- data_frame("species" = as.character())

for (i in 1:length(sta.simsum)){
  spec <- as.data.frame(sta.simsum[i]) %>%
    mutate(species = rownames(.)) %>%
    select(1,7,8) %>% #"average" (contribution to dissimilarity) and "species"
    filter(.[[2]] < 0.05 & .[[1]] > 0.05) %>% #pvalue<0.05, avg contribution>0.05
    select(species)
  sta.simspecies <- bind_rows(sta.simspecies, spec)
}

sta.species <- sta.simspecies %>%
  distinct()

# SIMPER but for depth zones and not cluster 
sta.simper.zones <- simper(sta.com[-1], sta.points$cat) 
sta.simsum.zones <- summary(sta.simper.zones)

sta.simspecies.zones <- data_frame("species" = as.character())
for (i in 1:length(sta.simsum.zones)){  #length is 21, nb of pairs among 7 clusters
  spec <- as.data.frame(sta.simsum.zones[i]) %>%
    mutate(species = rownames(.)) %>%
    select(1,7,8) %>% #"average" (contribution to dissimilarity) and "species"
    filter(.[[2]] < 0.05 & .[[1]] > 0.05) #pvalue<0.05, avg contribution>0.05
  sta.simspecies.zones <- bind_rows(sta.simspecies.zones, spec)
}
sta.simspecies.zones$location="Statia"

sta.loads <- 
  as_tibble(sta.pca$rotation, rownames = 'species') %>%
  inner_join(sta.species) %>%
  select(species, PC1, PC2)


pdf("2.Ch3.R.analysis/graphs/statia.pca.pdf") 
sta.pca.pooled.plot <- ggplot(sta.points) +
  geom_jitter(data = sta.points, aes(x = PC1, y = PC2, fill = cat),
              shape = 21, alpha = 0.5, color = "black") +
  geom_polygon(data = sta.hull,
               aes(x = PC1, y = PC2, fill = cat,
                   colour = cat,
                   group = as.factor(ord_num)),
               alpha = 0.3,
               show.legend = TRUE) +
  geom_segment(data = sta.loads,
               aes(x = 0, y = 0, xend = PC1*30, yend = PC2*30), lwd = 0.1,
               arrow = arrow(length = unit(1/2, 'picas'), type = "open")) +
  geom_text_repel(data = sta.loads, aes(x = PC1*30, 
                                        y = PC2*30, label = species), size = 4, segment.color = "grey69") +
  theme_classic() +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) 
 # guides(colour = guide_legend(nrow = 1), fill = guide_legend(nrow = 1))
sta.pca.pooled.plot
dev.off()

#### 4.3. Plots all together ####

## arrangement 1 
cluster.plots <- ggarrange(cur.plot, bon.plot.pooled+ rremove("ylab"),
                           cur.pca.plot,bon.pca.plot+ rremove("ylab"),
                           sta.plot.pooled+ rremove("ylab"),roa.plot+ rremove("ylab"),
                           sta.pca.pooled.plot+ rremove("ylab"), roa.pca.plot+ rremove("ylab"),
                           ncol = 2, nrow=4, common.legend=TRUE,legend="bottom")
cluster.plots

ggsave(cluster.plots, file = "plot.dendo.pca.png", width = 7, height = 12)

## arrangement 2 
cluster.plots <- ggarrange(cur.plot, cur.pca.plot + rremove("xlab"),
                           bon.plot.pooled,bon.pca.plot+ rremove("xlab"),
                           sta.plot.pooled,sta.pca.pooled.plot+ rremove("xlab"), 
                           roa.plot,roa.pca.plot,
                           ncol = 2, nrow=4, legend="none")
cluster.plots

ggsave(cluster.plots, file = "2.Ch3.R.analysis/graphs/combined.dendro.pca.png", width = 10, height = 14)

## top 5 species for PC1 and PC2 
pca.species.PC <- rbind(cur.pca.speciesPC1,cur.pca.speciesPC2,
                        bon.pca.speciesPC1,bon.pca.speciesPC2,
                        sta.pca.speciesPC1,sta.pca.speciesPC2,
                        roa.pca.speciesPC1,roa.pca.speciesPC2)
inner_join(cur.pca.speciesPC1,bon.pca.speciesPC1,by="species")

write.csv(pca.species.PC,"stats.main.pca.species.csv")

## significant species SIMPER for depth zones
simper.depthzones <- rbind(cur.simspecies.zones,bon.simspecies.zones,
                           sta.simspecies.zones,roa.simspecies.zones)
write.csv(simper.depthzones,"stat.simper.depthzones.csv")

#### 5. PERMANOVA to test depth zones ####
library(vegan)
## Curacao 
cur.permanova <- cur.com %>%
  mutate(zone=case_when(dband<=70 ~ "upper mesophotic",
                       dband<=120~ "lower mesophotic",
                       dband<=180~ "upper rariphotic",
                       dband>180~ "lower rariphotic")) %>%
  select(zone)
  
adonis2(cur.com[,-1]~cur.permanova$zone)
  

## Bonaire 
bon.permanova <- bon.com %>% 
  mutate(zone=case_when(dband.grouped<70 ~ "upper mesophotic",
                       dband.grouped<120~ "lower mesophotic",
                       dband.grouped<170~ "upper rariphotic",
                       dband.grouped>=170~ "lower rariphotic"))%>%
  select(zone)

adonis2(bon.com[,-1]~bon.permanova$zone)

## Statia
sta.permanova <- sta.com %>%
  mutate(zone=case_when(dband.grouped<100 ~ "upper mesophotic",
                       dband.grouped<140~ "lower mesophotic",
                       dband.grouped<190~ "upper rariphotic",
                       dband.grouped>=190~ "lower rariphotic")) %>%
  select(zone)

adonis2(sta.com[,-1]~sta.permanova$zone)

## Roatan
roa.permanova <- roa.com %>%
  mutate(zone=case_when(dband<100 ~ "upper mesophotic",
                       dband<150 ~ "lower mesophotic",
                       dband<210~ "upper rariphotic",
                       dband>=210~ "lower rariphotic")) %>%
  select(zone)

adonis2(roa.com[,-1]~roa.permanova$zone)