#---------------------
# Taxonomic Bias
# Analyses and Plots
#
# I. Alpha diversity by site (Table 1)
# II. Beta diversity by site (for Figure 2d-f) 
# III. Ratio of Occurrences that at iNat vs Specimens (Figure 2a-c)
# IV. Species richness by order by site (Figure 2d-f)
#---------------------
require(dplyr)
require(ggplot2)
require(tidyverse)
require(ggExtra)
require(cowplot)
require(grid)
require(gridExtra)

# Set working directory
setwd("")

#-----------------------
# LOAD IN DATA
#
dat <- read.csv("full_gbif_for_analysis.csv")
table(dat$place_name)

dat$place_name[dat$place_name %in% "Marble/Salmon Mountains"] <- "Marble Mountains"
dat$place_name[dat$place_name %in% "One Tam"] <- "Mount Tamalpais"

# Filter non-verified species and duplicates out 
dat %>% filter(verifiedTaxon=="verified") %>% filter(keepWhenDeduplicating=="keep") -> dats

#-----------------------
# I. Alpha diversity table: number of unique species per list
#

# alpha diversity by site (both lists combined)
dats %>% select(place_name, species) %>% distinct() %>% group_by(place_name) %>% count() %>%
  rename(alphaCombined=n) -> alphaT

# alpha diversity by data type  
dats %>% select(place_name, basisofrecord, species) %>% distinct() %>% group_by(place_name, basisofrecord) %>% 
  count() %>% spread(., basisofrecord, n) %>%
  rename(alpha_iNat = HUMAN_OBSERVATION, alpha_herb=PRESERVED_SPECIMEN) -> alphaEach

# Alpha diversity table: number of unique species per list
dats %>% select(place_name, Origin, basisofrecord, species) %>% distinct() %>%
  group_by(place_name, Origin, basisofrecord, species) %>% count() %>% ungroup() %>% 
  spread(., basisofrecord, n) %>% 
  filter(is.na(HUMAN_OBSERVATION) | is.na(PRESERVED_SPECIMEN)) %>% 
  mutate(unique.who = if_else(is.na(HUMAN_OBSERVATION), "specimen", "inat")) %>% 
  group_by(place_name,unique.who,Origin) %>% count() %>% ungroup() %>%
  mutate(type=paste("unique", unique.who, Origin, sep = "_")) %>% select(-unique.who, - Origin) %>%
  spread(., type, n) -> alpha_unique  

# Combine and reorignize to make table 
full_join(alphaT, alphaEach, by=c("place_name")) %>%
  full_join(alpha_unique, by=c("place_name")) %>% 
  select(place_name, alphaCombined, alpha_herb, unique_specimen_native, unique_specimen_naturalized, 
         alpha_iNat, unique_inat_native, unique_inat_naturalized) -> TableOne

write.csv(TableOne, "Table1.csv")

#---------------------------------------- 
# II. Beta diversity - Jaccard dissimilarity 
# Calculate Beta diversity for each site 

require(vegan)

places <- unique(dats$place_name)

beta.out <- data.frame()

for(i in 1:length(places)){
  #  i=1
  a.list <- dats %>% filter(place_name==places[i])
  
  a.list %>% group_by(basisofrecord, species) %>% count() %>% ungroup() %>% 
    spread(., basisofrecord, n) %>% 
    # Change non zero values to ones to make it a presence absence table
    mutate_if(is.numeric, ~1 * (. > 0)) %>%
    mutate_if(is.numeric, ~replace_na(.,0)) -> a_wide
  
  nrow(a_wide) # 
  ncol(a_wide) # 
  
  a <- data.frame(t(a_wide))
  #View(a)
  # Remove the species name from the top row
  a <- a[-1,]
  a <- mutate_if(a, is.character, as.numeric)
  
  a <- as.matrix(a)
  
  dist.j <- vegdist(a, method="jaccard", binary = TRUE) %>% as.matrix()
  
  data.frame(beta=dist.j[1,2], place_name=places[i]) %>% rbind(beta.out) -> beta.out
  
}

# Values inserted into Figure 2
beta.out %>% mutate(beta = round(beta, 2))

#-----------
# III. Ratio of Occurrences that at iNat vs Specimens (Figure 2a-c)

# Place level sample size for ratios
mm.dat <- dats %>% filter(place_name=="Marble Mountains") 

# MM
dats %>% select(place_name, basisofrecord,  order, species) %>%
  group_by(place_name, basisofrecord,  order, species) %>% count() %>% 
  spread(., basisofrecord, n) %>% filter(place_name=="Marble Mountains") %>%
  mutate(HUMAN_OBSERVATION=replace_na(HUMAN_OBSERVATION, 0), PRESERVED_SPECIMEN=replace_na(PRESERVED_SPECIMEN, 0)) %>%
  mutate(sp.n = HUMAN_OBSERVATION + PRESERVED_SPECIMEN) %>% 
  mutate(HO.prop = HUMAN_OBSERVATION/sp.n, PS.prop=PRESERVED_SPECIMEN/sp.n) %>% 
  select(place_name, order, species, HO.prop, PS.prop) %>% 
  # Order by frequency of Human obs
  arrange(desc(HO.prop)) %>% 
  tibble::rowid_to_column(var="rank") %>%
  # Percent of observations by species dominated by each data source 
  #filter(HO.prop > PS.prop) %>% nrow() /(length(unique(mm.dat$species))) # 0.116
  #filter(HO.prop < PS.prop) %>% nrow() / (length(unique(mm.dat$species))) # 0.862
  #filter(HO.prop == PS.prop) %>% nrow() / (length(unique(mm.dat$species))) # 0.021  
  pivot_longer(cols = c("HO.prop", "PS.prop"), names_to = "basisofrecord",
               names_prefix = "wk", values_to = "count",values_drop_na = TRUE) %>% 
  ggplot(aes(fill=basisofrecord, y=count, x=reorder(species, -rank))) + 
  geom_bar(position="fill", stat="identity", width=1) + theme(legend.position="none") + 
  scale_fill_manual(values=c("HO.prop"="darkgoldenrod1", "PS.prop"="darkslategray")) +
  geom_hline(yintercept = 0.5, color="gray99", linetype="dashed") +
  ylab("Ratio of iNat to herbarium observations") + xlab("") +
  theme(axis.title.x = element_blank(), axis.text.x  =element_blank(),axis.ticks.x=element_blank(),
        axis.title.y = element_blank(), axis.text.y  =element_text(size=6))   -> mm.rank

# TAM
tam.dat <- dats %>% filter(place_name=="Mount Tamalpais") 

dats %>% select(place_name, basisofrecord,  order, species) %>%
  group_by(place_name, basisofrecord, order, species) %>% count() %>% 
  spread(., basisofrecord, n) %>% filter(place_name=="Mount Tamalpais") %>%
  mutate(HUMAN_OBSERVATION=replace_na(HUMAN_OBSERVATION, 0), PRESERVED_SPECIMEN=replace_na(PRESERVED_SPECIMEN, 0)) %>%
  mutate(sp.n = HUMAN_OBSERVATION + PRESERVED_SPECIMEN) %>% 
  mutate(HO.prop = HUMAN_OBSERVATION/sp.n, PS.prop=PRESERVED_SPECIMEN/sp.n) %>% 
  select(place_name, order, species, HO.prop, PS.prop) %>% 
  # Order by frequency of Human obs
  arrange(desc(HO.prop)) %>% tibble::rowid_to_column(var="rank") %>%
  # Percent of observations by species dominated by each data source 
  #filter(HO.prop > PS.prop) %>% nrow() / (length(unique(tam.dat$species))) # 0.710
  #filter(HO.prop < PS.prop) %>% nrow() / (length(unique(tam.dat$species))) # 0.2361
  #filter(HO.prop == PS.prop) %>% nrow() / (length(unique(tam.dat$species))) # 0.0545
  pivot_longer(cols = c("HO.prop", "PS.prop"), names_to = "basisofrecord",
               names_prefix = "wk", values_to = "count",values_drop_na = TRUE) %>% 
  ggplot(aes(fill=basisofrecord, y=count, x=reorder(species, -rank))) + 
  geom_bar(position="fill", stat="identity", width=1) + theme(legend.position="none") + 
  scale_fill_manual(values=c("HO.prop"="darkgoldenrod1", "PS.prop"="darkslategray")) +
  geom_hline(yintercept = 0.5, color="gray99", linetype="dashed") +
  ylab("Ratio of iNat to herbarium observations") + xlab("") +
  theme(axis.title.x = element_blank(), axis.text.x  =element_blank(),axis.ticks.x=element_blank(),
        axis.title.y = element_blank(), axis.text.y  =element_text(size=6))  -> tam.rank

# SMM
smm.dat <- dats %>% filter(place_name=="Santa Monica Mountains") 

dats %>% select(place_name, basisofrecord,  order, species) %>%
  group_by(place_name, basisofrecord,  order, species) %>% count() %>% 
  spread(., basisofrecord, n) %>% filter(place_name=="Santa Monica Mountains") %>%
  mutate(HUMAN_OBSERVATION=replace_na(HUMAN_OBSERVATION, 0), PRESERVED_SPECIMEN=replace_na(PRESERVED_SPECIMEN, 0)) %>%
  mutate(sp.n = HUMAN_OBSERVATION + PRESERVED_SPECIMEN) %>% 
  mutate(HO.prop = HUMAN_OBSERVATION/sp.n, PS.prop=PRESERVED_SPECIMEN/sp.n) %>% 
  select(place_name,  order, species, HO.prop, PS.prop) %>% 
  # Order by frequency of Human obs
  arrange(desc(HO.prop)) %>% 
  tibble::rowid_to_column(var="rank") %>%
  # Percent of observations by species dominated by each data source 
  #filter(HO.prop > PS.prop) %>% nrow() / (length(unique(smm.dat$species))) # 0.516
  #filter(HO.prop < PS.prop) %>% nrow() / (length(unique(smm.dat$species))) # 0.454
  #filter(HO.prop == PS.prop) %>% nrow() / (length(unique(smm.dat$species))) # 0.029
  pivot_longer(cols = c("HO.prop", "PS.prop"), names_to = "basisofrecord",
               names_prefix = "wk", values_to = "count",values_drop_na = TRUE) %>% 
  ggplot(aes(fill=basisofrecord, y=count, x=reorder(species, -rank))) + 
  geom_bar(position="fill", stat="identity", width=1) + theme(legend.position="none") + 
  scale_fill_manual(values=c("HO.prop"="darkgoldenrod1", "PS.prop"="darkslategray")) +
  geom_hline(yintercept = 0.5, color="gray99", linetype="dashed") +
  ylab("Ratio of iNat to herbarium observations") + xlab("") +
  theme(axis.title.x = element_blank(), axis.text.x  =element_blank(),axis.ticks.x=element_blank(),
        axis.title.y = element_blank(), axis.text.y  =element_text(size=6))   -> smm.rank

#------ Combine
a5 <- ggplotGrob(mm.rank) 
a6 <- ggplotGrob(tam.rank) 
a7 <- ggplotGrob(smm.rank) 

combo <- cowplot::plot_grid(a5, a6, a7, ncol=3,  
                            labels=c('(a)', '(b)', '(c)'),
                            label_fontface = "plain",
                            label_size = 8, hjust = 0, label_x = 0.01, align = "h")

y.grob <- textGrob("Ratio of herbarium (green) to iNaturalist (gold) records", gp=gpar(fontface="plain", fontsize=8), rot=90)

#add to plot
Figure2a <- grid.arrange(arrangeGrob(combo, left = y.grob))

ggsave("Figure2_ObsRatio_HerbInat.pdf", Figure2a, dpi = 600, width=6, height=2, units="in")


#------------------------
# IV. Species richness by order by site (Figure 2d-f)
# Summary showing the unique species for each basis of record

##### MT
dats %>% filter(place_name=="Mount Tamalpais") %>% select(basisofrecord, order, species) %>%
  distinct() %>% mutate(occ=1) %>%
  spread(., basisofrecord, occ) %>%  
  mutate(HUMAN_OBSERVATION = replace_na(HUMAN_OBSERVATION, 0), PRESERVED_SPECIMEN = replace_na(PRESERVED_SPECIMEN, 0)) %>%
  mutate(uniqueStat = if_else(HUMAN_OBSERVATION==1 & PRESERVED_SPECIMEN==1, "both",
                              if_else(HUMAN_OBSERVATION==1 & PRESERVED_SPECIMEN==0, "iNatOnly", 
                                      if_else(PRESERVED_SPECIMEN==1 & HUMAN_OBSERVATION==0, "HerbOnly", "check")))) %>%
  #group_by(uniqueStat) %>% count()
  select(species, uniqueStat) %>% full_join(tam.dat, by=c("species")) %>%
  select(basisofrecord, order, uniqueStat, species) %>% distinct() %>% 
  group_by(basisofrecord, order, uniqueStat) %>% count() %>% arrange(order) %>% ungroup() %>% 
  mutate(uniqueStat = if_else(uniqueStat=="both", paste(basisofrecord, "both", sep="_"), uniqueStat)) %>% 
  mutate(uniqueStat = factor(uniqueStat, levels=c("HerbOnly", "iNatOnly", "HUMAN_OBSERVATION_both", 
                                                  "PRESERVED_SPECIMEN_both"))) %>%
  mutate(basisofrecord = factor(basisofrecord, levels=c("PRESERVED_SPECIMEN", "HUMAN_OBSERVATION"))) %>%
  ggplot() +
  geom_bar(aes(x = basisofrecord, y = n, fill = uniqueStat),  position = "stack", stat = "identity") +
  facet_grid(~ order, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) +
  scale_fill_manual(values=c("HUMAN_OBSERVATION_both"="darkgoldenrod1", "iNatOnly"="#DC863B",
                             "HerbOnly"="#A2A475", "PRESERVED_SPECIMEN_both"="darkslategray")) +
  ylab("Species richness by order") +  theme(legend.position="none") + 
  theme(axis.title.x = element_blank(), axis.text.x  =element_blank(),axis.ticks.x=element_blank(),
        axis.title.y = element_blank(), axis.text.y  =element_text(size=6),
        strip.text.x.bottom = element_text(angle=90)) -> tam.order.plot

##### SMM
dats %>% filter(place_name=="Santa Monica Mountains") %>% select(basisofrecord, order, species) %>%
  distinct() %>% mutate(occ=1) %>%
  spread(., basisofrecord, occ) %>%  
  mutate(HUMAN_OBSERVATION = replace_na(HUMAN_OBSERVATION, 0), PRESERVED_SPECIMEN = replace_na(PRESERVED_SPECIMEN, 0)) %>%
  mutate(uniqueStat = if_else(HUMAN_OBSERVATION==1 & PRESERVED_SPECIMEN==1, "both",
                              if_else(HUMAN_OBSERVATION==1 & PRESERVED_SPECIMEN==0, "iNatOnly", 
                                      if_else(PRESERVED_SPECIMEN==1 & HUMAN_OBSERVATION==0, "HerbOnly", "check")))) %>%
  #group_by(uniqueStat) %>% count()
  select(species, uniqueStat) %>% full_join(smm.dat, by=c("species")) %>%
  select(basisofrecord, order, uniqueStat, species) %>% distinct() %>% 
  group_by(basisofrecord, order, uniqueStat) %>% count() %>% arrange(order) %>% ungroup() %>% 
  mutate(uniqueStat = if_else(uniqueStat=="both", paste(basisofrecord, "both", sep="_"), uniqueStat)) %>% 
  mutate(uniqueStat = factor(uniqueStat, levels=c("HerbOnly", "iNatOnly", "HUMAN_OBSERVATION_both", 
                                                  "PRESERVED_SPECIMEN_both"))) %>%
  mutate(basisofrecord = factor(basisofrecord, levels=c("PRESERVED_SPECIMEN", "HUMAN_OBSERVATION"))) %>%
  ggplot() +
  geom_bar(aes(x = basisofrecord, y = n, fill = uniqueStat),  position = "stack", stat = "identity") +
  facet_grid(~ order, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) +
  scale_fill_manual(values=c("HUMAN_OBSERVATION_both"="darkgoldenrod1", "iNatOnly"="#DC863B",
                             "HerbOnly"="#A2A475", "PRESERVED_SPECIMEN_both"="darkslategray")) +
  ylab("Species richness by order") +  theme(legend.position="none") + 
  theme(axis.title.x = element_blank(), axis.text.x  =element_blank(),axis.ticks.x=element_blank(),
        axis.title.y = element_blank(), axis.text.y  =element_text(size=6),
        strip.text.x.bottom = element_text(angle=90)) -> smm.order.plot

##### MM
dats %>% filter(place_name=="Marble Mountains") %>% select(basisofrecord, order, species) %>%
  distinct() %>% mutate(occ=1) %>%
  spread(., basisofrecord, occ) %>%  
  mutate(HUMAN_OBSERVATION = replace_na(HUMAN_OBSERVATION, 0), PRESERVED_SPECIMEN = replace_na(PRESERVED_SPECIMEN, 0)) %>%
  mutate(uniqueStat = if_else(HUMAN_OBSERVATION==1 & PRESERVED_SPECIMEN==1, "both",
                              if_else(HUMAN_OBSERVATION==1 & PRESERVED_SPECIMEN==0, "iNatOnly", 
                                      if_else(PRESERVED_SPECIMEN==1 & HUMAN_OBSERVATION==0, "HerbOnly", "check")))) %>%
  #group_by(uniqueStat) %>% count()
  select(species, uniqueStat) %>% full_join(mm.dat, by=c("species")) %>% 
  select(basisofrecord, order, uniqueStat, species) %>% distinct() %>% 
  group_by(basisofrecord, order, uniqueStat) %>% count() %>% arrange(order) %>% ungroup() %>% 
  mutate(uniqueStat = if_else(uniqueStat=="both", paste(basisofrecord, "both", sep="_"), uniqueStat)) %>% 
  mutate(uniqueStat = factor(uniqueStat, levels=c("HerbOnly", "iNatOnly", "HUMAN_OBSERVATION_both", 
                                                  "PRESERVED_SPECIMEN_both"))) %>%
  mutate(basisofrecord = factor(basisofrecord, levels=c("PRESERVED_SPECIMEN", "HUMAN_OBSERVATION"))) %>%
  ggplot() +
  geom_bar(aes(x = basisofrecord, y = n, fill = uniqueStat),  position = "stack", stat = "identity") +
  facet_grid(~ order, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) +
  scale_fill_manual(values=c("HUMAN_OBSERVATION_both"="darkgoldenrod1", "iNatOnly"="#DC863B",
                             "HerbOnly"="#A2A475", "PRESERVED_SPECIMEN_both"="darkslategray")) +
  ylab("Species richness by order") +  theme(legend.position="none") + 
  theme(axis.title.x = element_blank(), axis.text.x  =element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_blank(), axis.text.y =element_text(size=6), 
        strip.text.x.bottom = element_text(angle=90))  -> mm.order.plot

#------ Combine
a1 <- ggplotGrob(mm.order.plot) 
a2 <- ggplotGrob(tam.order.plot) 
a3 <- ggplotGrob(smm.order.plot) 

combo <- cowplot::plot_grid(a1, a2, a3, ncol=1,  
                            labels=c('(d)', '(e)', '(f)'),
                            label_fontface = "plain",
                            label_size = 8, hjust = 0, label_x = 0.01, align = "h")

y.grob <- textGrob("Species richness by order", gp=gpar(fontface="plain", fontsize=8), rot=90)

#add to plot
Figure2d <- grid.arrange(arrangeGrob(combo, left = y.grob)) 

ggsave("Figure2_OrderRichness_HerbInat.pdf", Figure2d, dpi = 600, width=7, height=8, units="in")


