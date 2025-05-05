#---------------------
#
# Phylogenetic diversity
# I. Calculate PD and Pagel's lambda 
# II. Make circle charts (Figure 3)
#---------------------
require(dplyr)
require(ggplot2)
require(tidyverse)
library(canaper)
require(ape)
require(picante)
require(janitor)
require(phytools)
library(ggtree)
library(treeio)
library(patchwork)
#BiocManager::install("ggtreeExtra")
library(ggtreeExtra)

#------
# Set Working directory
setwd()

#-----------------------
# LOAD IN DATA

# Read in the tree from Thornhill et al. 2017. Sequence matrix and tree files for Spatial phylogenetics of the native California flora (Thronhill et al. SMC Biology)
# https://datadryad.org/dataset/doi:10.6078/D1VD4P
ca_tree <- read.nexus("Californian_clades_tree_final.nex")

# Records
dat <- read.csv("full_gbif_for_analysis.csv")

dat$place_name[dat$place_name %in% "Marble/Salmon Mountains"] <- "Marble Mountains"
dat$place_name[dat$place_name %in% "One Tam"] <- "Mount Tamalpais"

# Filter to verified species and remove duplicates 
dat %>% filter(verifiedTaxon=="verified") %>% filter(keepWhenDeduplicating=="keep") -> dats

# Do all the native species shave an OTU?
dats %>% filter(Origin=="native") %>% filter(is.na(OTU))
# Only one doesn't: Allium crenulatum - NA-technically this has not yet been added to the Jepson, although native, so don't include in analysis

#-----------------------
# I. Calculate PD and Pagel's lambda 

#----------------- LOOP IT
place <- unique(dats$place_name)
pd.out <- data.frame()
plam.out <- data.frame()

for(i in 1:length(place)){
  
 # i=1
  a <- dats[dats$place_name==place[i],]

  # PD based on basis of record
  a %>% select(basisofrecord, species, OTU, Origin) %>% 
    filter(!is.na(OTU)) %>%
    # Who doesn't have an OTU?
    #filter(is.na(OTU)) %>% filter(Origin=="native") %>% distinct() # no mismatches
    #group_by(basisofrecord, OTU) %>% count() %>% rename(occ=n) %>%
    select(basisofrecord, OTU) %>% mutate(occ=1) %>%  distinct() %>%
    spread(., basisofrecord, occ) %>% 
    mutate_if(is.numeric, ~replace_na(.,0)) -> site_wide

  site_matrix <- data.frame(t(site_wide)) %>% row_to_names(row_number = 1) %>% mutate_if(is.character, as.numeric)
  
  pd(site_matrix, ca_tree) %>% data.frame() %>% rownames_to_column(var="basisofrecord") %>% 
    mutate(place_name=place[i]) %>% rbind(pd.out) -> pd.out
    
  # PD data sources together
  a %>% select(place_name, species, OTU, Origin) %>% 
    # Who doesn't have an OTU?
    #filter(is.na(OTU)) %>% filter(Origin=="native") %>% distinct() # all the nonnative species so looks good
    filter(!is.na(OTU)) %>% 
    select(place_name, OTU) %>% mutate(occ=1) %>%  distinct() %>%
    spread(., place_name, occ) %>% 
    mutate_if(is.numeric, ~replace_na(.,0)) -> together_wide
  
  together_matrix <- data.frame(t(together_wide)) %>% row_to_names(row_number = 1) %>% mutate_if(is.character, as.numeric)
  
  pd(together_matrix, ca_tree) %>% data.frame() %>% rownames_to_column(var="basisofrecord") %>% 
    mutate(place_name=place[i]) %>% rbind(pd.out) -> pd.out

  # Pagel's lambda
  # iNaturalist observations
  a %>% select(basisofrecord, species, OTU, Origin) %>% 
    filter(!is.na(OTU)) %>% filter(basisofrecord=="HUMAN_OBSERVATION") %>% 
    group_by(OTU) %>% count() %>% 
    column_to_rownames(var="OTU") -> a.nat
  
  Obs.nat <-setNames(log(a.nat$n), rownames(a.nat))
  
  pnat <- phylosig(ca_tree,Obs.nat,method="lambda",test=TRUE)
  
  data.frame(plambda = round(pnat$lambda,2), pval= round(pnat$P,3), place_name=place[i], basisofrecord="HUMAN_OBSERVATION") %>% 
    rbind(plam.out) -> plam.out
  
  # Specimen records 
  a %>% select(basisofrecord, species, OTU, Origin) %>% 
    filter(!is.na(OTU)) %>% filter(basisofrecord=="PRESERVED_SPECIMEN") %>% 
    group_by(OTU) %>% count() %>% 
    column_to_rownames(var="OTU") -> a.herb
  
  Obs.herb <-setNames(log(a.herb$n), rownames(a.herb))
  
  pherb <- phylosig(ca_tree, Obs.herb, method="lambda",test=TRUE)
  
  data.frame(plambda = round(pherb$lambda,2), pval= round(pherb$P,3), place_name=place[i], basisofrecord="PRESERVED_SPECIMEN") %>% 
    rbind(plam.out) -> plam.out
  
  # All records together
  a %>% select(basisofrecord, species, OTU, Origin) %>% 
    filter(!is.na(OTU)) %>% group_by(OTU) %>% count() %>% 
    column_to_rownames(var="OTU") -> a.both
  
  Obs.both <-setNames(log(a.both$n), rownames(a.both))
  
  pboth <- phylosig(ca_tree, Obs.both, method="lambda",test=TRUE)
  
  data.frame(plambda = round(pboth$lambda,2), pval= round(pboth$P,3), place_name=place[i], basisofrecord="combined") %>% 
    rbind(plam.out) -> plam.out
  
  }


# Summary Table, for combined PD values 
pd.out %>%  mutate(pdval = paste(round(PD, 2)," (", SR , ")", sep="")) %>%
  select(place_name, basisofrecord, pdval) %>%
  spread(., basisofrecord, pdval) 

plam.out %>% mutate(plamval = paste(plambda," (", pval , ")", sep="")) %>%
  select(place_name, basisofrecord, plamval) %>%
  spread(., basisofrecord, plamval) -> plam.summary


#-----------------------
# II. Make circle charts


#---------------- Marble Mountains
mm.dat <- dats %>% filter(place_name=="Marble Mountains")
# mm tree
mm_tips2include <- as.vector(na.omit(unique(mm.dat$OTU))) # 438 unique species on MM
mmTree_temp <- drop.tip(ca_tree, ca_tree$tip.label[-match(mm_tips2include, ca_tree$tip.label)])

# base tree
mmp <- ggtree(mmTree_temp, layout = "circular", size = 0.25, open.angle = 10) + xlim(-1,4.5) + ylim(0,450)
mmTree <- mmp %<+% mm.dat

# counts for each dataset and combined
mm.dat.obs <- mm.dat %>% group_by(basisofrecord) %>% filter(basisofrecord == "HUMAN_OBSERVATION") %>% group_by(OTU) %>% drop_na(OTU) %>% count()
mm.dat.ps <- mm.dat %>% group_by(basisofrecord) %>% filter(basisofrecord == "PRESERVED_SPECIMEN") %>% group_by(OTU) %>% drop_na(OTU) %>% count()
mm.dat.comb <- mm.dat %>% group_by(OTU) %>% drop_na(OTU) %>% count()

# iNat observations on mm
mmTree_obs <- mmTree + geom_fruit(data = mm.dat.obs,
                                  geom = geom_col,
                                  mapping = aes(y = OTU, x = sqrt(n)),
                                  fill = 'darkgoldenrod1',
                                  size = 0.2,
                                  pwidth = 1.45
) +
  theme(plot.margin = margin(.1,.1,.1,.1)) +
  annotate("text", x = -0.55, y = 150, label = paste0("PD==", round(as.numeric(pd.out[pd.out$basisofrecord=="HUMAN_OBSERVATION" & pd.out$place_name=="Marble Mountains","PD"]),2)), parse = TRUE, size = 2.5) +
  annotate("text", x = -0.55, y = 200, label = paste0("OTU==", as.numeric(pd.out[pd.out$basisofrecord=="HUMAN_OBSERVATION" & pd.out$place_name=="Marble Mountains","SR"])), parse = TRUE, size = 2.5) +
  annotate("text", x = -0.55, y = 250, label = paste0("lambda==", as.numeric(plam.out[plam.out$basisofrecord=="HUMAN_OBSERVATION" & pd.out$place_name=="Marble Mountains","plambda"])), parse = TRUE, size = 2.5) +
  annotate("text", x = -0.55, y = 300, label = paste0("p==", as.numeric(plam.out[plam.out$basisofrecord=="HUMAN_OBSERVATION" & pd.out$place_name=="Marble Mountains","pval"])), parse = TRUE, size = 2.5) 



# Herbarium records on mm
mmTree_ps <- mmTree + geom_fruit(data = mm.dat.ps,
                                 geom = geom_col,
                                 mapping = aes(y = OTU, x = sqrt(n)),
                                 fill = 'darkslategray',
                                 pwidth = 1.45) +
  theme(plot.margin = margin(.1,.1,.1,.1)) +
  annotate("text", x = -0.55, y = 150, label = paste0("PD==", round(as.numeric(pd.out[pd.out$basisofrecord=="PRESERVED_SPECIMEN" & pd.out$place_name=="Marble Mountains","PD"]),2)), parse = TRUE, size = 2.5) +
  annotate("text", x = -0.55, y = 200, label = paste0("OTU==", as.numeric(pd.out[pd.out$basisofrecord=="PRESERVED_SPECIMEN" & pd.out$place_name=="Marble Mountains","SR"])), parse = TRUE, size = 2.5) +
  annotate("text", x = -0.55, y = 250, label = paste0("lambda==", as.numeric(plam.out[plam.out$basisofrecord=="PRESERVED_SPECIMEN" & pd.out$place_name=="Marble Mountains","plambda"])), parse = TRUE, size = 2.5) +
  annotate("text", x = -0.55, y = 300, label = paste0("p==", as.numeric(plam.out[plam.out$basisofrecord=="PRESERVED_SPECIMEN" & pd.out$place_name=="Marble Mountains","pval"])), parse = TRUE, size = 2.5) 



# iNat + Herbarium records on mm
mmTree_comb <- mmTree + geom_fruit(data = mm.dat.comb,
                                   geom = geom_col,
                                   mapping = aes(y = OTU, x = sqrt(n)),
                                   fill = 'mediumpurple4',
                                   pwidth = 1.45,
) +
  theme(plot.margin = margin(.1,.1,.1,.1)) +
  annotate("text", x = -0.55, y = 150, label = paste0("PD==", round(as.numeric(pd.out[pd.out$basisofrecord=="Marble Mountains" & pd.out$place_name=="Marble Mountains","PD"]),2)), parse = TRUE, size = 2.5) +
  annotate("text", x = -0.55, y = 200, label = paste0("OTU==", as.numeric(pd.out[pd.out$basisofrecord=="Marble Mountains" & pd.out$place_name=="Marble Mountains","SR"])), parse = TRUE, size = 2.5) +
  annotate("text", x = -0.55, y = 250, label = paste0("lambda==", as.numeric(plam.out[plam.out$basisofrecord=="combined" & pd.out$place_name=="Marble Mountains","plambda"])), parse = TRUE, size = 2.5) +
  annotate("text", x = -0.55, y = 300, label = paste0("p==", as.numeric(plam.out[plam.out$basisofrecord=="combined" & pd.out$place_name=="Marble Mountains","pval"])), parse = TRUE, size = 2.5) 



#--- Mount Tamalpais
tam.dat <- dats %>% filter(place_name=="Mount Tamalpais")

tam_tips2include <- as.vector(na.omit(unique(tam.dat$OTU))) # 403 unique species on MT
tamTree_temp <- drop.tip(ca_tree, ca_tree$tip.label[-match(tam_tips2include, ca_tree$tip.label)])

# base tree
tamp <- ggtree(tamTree_temp, layout = "circular", size = 0.25, open.angle = 10) + xlim(-1,4.5) + ylim(0,450)
tamTree <- tamp %<+% tam.dat

# counts for each dataset and combined
tam.dat.obs <- tam.dat %>% group_by(basisofrecord) %>% filter(basisofrecord == "HUMAN_OBSERVATION") %>% group_by(OTU) %>% drop_na(OTU) %>% count()
tam.dat.ps <- tam.dat %>% group_by(basisofrecord) %>% filter(basisofrecord == "PRESERVED_SPECIMEN") %>% group_by(OTU) %>% drop_na(OTU) %>% count()
tam.dat.comb <- tam.dat %>% group_by(OTU) %>% drop_na(OTU) %>% count()

# iNat observations on Tam
tamTree_obs <- tamTree + geom_fruit(data = tam.dat.obs,
                                    geom = geom_col,
                                    mapping = aes(y = OTU, x = sqrt(n)),
                                    fill = 'darkgoldenrod1',
                                    size = 0.2,
                                    pwidth = 1.45
) +
  theme(plot.margin = margin(.1,.1,.1,.1)) +
  annotate("text", x = -0.55, y = 150, label = paste0("PD==", round(as.numeric(pd.out[pd.out$basisofrecord=="HUMAN_OBSERVATION" & pd.out$place_name=="Mount Tamalpais","PD"]),2)), parse = TRUE, size = 2.5) +
  annotate("text", x = -0.55, y = 200, label = paste0("OTU==", as.numeric(pd.out[pd.out$basisofrecord=="HUMAN_OBSERVATION" & pd.out$place_name=="Mount Tamalpais","SR"])), parse = TRUE, size = 2.5) +
  annotate("text", x = -0.55, y = 250, label = paste0("lambda==", as.numeric(plam.out[plam.out$basisofrecord=="HUMAN_OBSERVATION" & pd.out$place_name=="Mount Tamalpais","plambda"])), parse = TRUE, size = 2.5) +
  annotate("text", x = -0.55, y = 300, label = paste0("p==", as.numeric(plam.out[plam.out$basisofrecord=="HUMAN_OBSERVATION" & pd.out$place_name=="Mount Tamalpais","pval"])), parse = TRUE, size = 2.5) 



# Herbarium records on Tam
tamTree_ps <- tamTree + geom_fruit(data = tam.dat.ps,
                                   geom = geom_col,
                                   mapping = aes(y = OTU, x = sqrt(n)),
                                   fill = 'darkslategray',
                                   pwidth = 1.45) +
  theme(plot.margin = margin(.1,.1,.1,.1)) +
  annotate("text", x = -0.55, y = 150, label = paste0("PD==", round(as.numeric(pd.out[pd.out$basisofrecord=="PRESERVED_SPECIMEN" & pd.out$place_name=="Mount Tamalpais","PD"]),2)), parse = TRUE, size = 2.5) +
  annotate("text", x = -0.55, y = 200, label = paste0("OTU==", as.numeric(pd.out[pd.out$basisofrecord=="PRESERVED_SPECIMEN" & pd.out$place_name=="Mount Tamalpais","SR"])), parse = TRUE, size = 2.5) +
  annotate("text", x = -0.55, y = 250, label = paste0("lambda==", as.numeric(plam.out[plam.out$basisofrecord=="PRESERVED_SPECIMEN" & pd.out$place_name=="Mount Tamalpais","plambda"])), parse = TRUE, size = 2.5) +
  annotate("text", x = -0.55, y = 300, label = paste0("p==", as.numeric(plam.out[plam.out$basisofrecord=="PRESERVED_SPECIMEN" & pd.out$place_name=="Mount Tamalpais","pval"])), parse = TRUE, size = 2.5) 



# iNat + Herbarium records on Tam
tamTree_comb <- tamTree + geom_fruit(data = tam.dat.comb,
                                     geom = geom_col,
                                     mapping = aes(y = OTU, x = sqrt(n)),
                                     fill = 'mediumpurple4',
                                     pwidth = 1.45,
) +
  theme(plot.margin = margin(.1,.1,.1,.1)) +
  annotate("text", x = -0.55, y = 150, label = paste0("PD==", round(as.numeric(pd.out[pd.out$basisofrecord=="Mount Tamalpais" & pd.out$place_name=="Mount Tamalpais","PD"]),2)), parse = TRUE, size = 2.5) +
  annotate("text", x = -0.55, y = 200, label = paste0("OTU==", as.numeric(pd.out[pd.out$basisofrecord=="Mount Tamalpais" & pd.out$place_name=="Mount Tamalpais","SR"])), parse = TRUE, size = 2.5) +
  annotate("text", x = -0.55, y = 250, label = paste0("lambda==", as.numeric(plam.out[plam.out$basisofrecord=="combined" & pd.out$place_name=="Mount Tamalpais","plambda"])), parse = TRUE, size = 2.5) +
  annotate("text", x = -0.55, y = 300, label = paste0("p==", as.numeric(plam.out[plam.out$basisofrecord=="combined" & pd.out$place_name=="Mount Tamalpais","pval"])), parse = TRUE, size = 2.5) 


#----------------- Santa Monica Mountains
smm.dat <- dats %>% filter(place_name=="Santa Monica Mountains")

smm_tips2include <- as.vector(na.omit(unique(smm.dat$OTU))) # 435 unique species on MT
smmTree_temp <- drop.tip(ca_tree, ca_tree$tip.label[-match(smm_tips2include, ca_tree$tip.label)])

# base tree
smmp <- ggtree(smmTree_temp, layout = "circular", size = 0.25, open.angle = 10) + xlim(-1,4.5) + ylim(0,450)
smmTree <- smmp %<+% smm.dat

# counts for each dataset and combined
smm.dat.obs <- smm.dat %>% group_by(basisofrecord) %>% filter(basisofrecord == "HUMAN_OBSERVATION") %>% group_by(OTU) %>% drop_na(OTU) %>% count()
smm.dat.ps <- smm.dat %>% group_by(basisofrecord) %>% filter(basisofrecord == "PRESERVED_SPECIMEN") %>% group_by(OTU) %>% drop_na(OTU) %>% count()
smm.dat.comb <- smm.dat %>% group_by(OTU) %>% drop_na(OTU) %>% count()

# iNat observations on smm
smmTree_obs <- smmTree + geom_fruit(data = smm.dat.obs,
                                    geom = geom_col,
                                    mapping = aes(y = OTU, x = sqrt(n)),
                                    fill = 'darkgoldenrod1',
                                    size = 0.2,
                                    pwidth = 1.45
) +
  theme(plot.margin = margin(.1,.1,.1,.1)) +
  annotate("text", x = -0.55, y = 150, label = paste0("PD==", round(as.numeric(pd.out[pd.out$basisofrecord=="HUMAN_OBSERVATION" & pd.out$place_name=="Santa Monica Mountains","PD"]),2)), parse = TRUE, size = 2.5) +
  annotate("text", x = -0.55, y = 200, label = paste0("OTU==", as.numeric(pd.out[pd.out$basisofrecord=="HUMAN_OBSERVATION" & pd.out$place_name=="Santa Monica Mountains","SR"])), parse = TRUE, size = 2.5) +
  annotate("text", x = -0.55, y = 250, label = paste0("lambda==", as.numeric(plam.out[plam.out$basisofrecord=="HUMAN_OBSERVATION" & pd.out$place_name=="Santa Monica Mountains","plambda"])), parse = TRUE, size = 2.5) +
  annotate("text", x = -0.55, y = 300, label = paste0("p==", as.numeric(plam.out[plam.out$basisofrecord=="HUMAN_OBSERVATION" & pd.out$place_name=="Santa Monica Mountains","pval"])), parse = TRUE, size = 2.5) 


# Herbarium records on smm
smmTree_ps <- smmTree + geom_fruit(data = smm.dat.ps,
                                   geom = geom_col,
                                   mapping = aes(y = OTU, x = sqrt(n)),
                                   fill = 'darkslategray',
                                   pwidth = 1.45) +
  theme(plot.margin = margin(.1,.1,.1,.1)) +
  annotate("text", x = -0.55, y = 150, label = paste0("PD==", round(as.numeric(pd.out[pd.out$basisofrecord=="PRESERVED_SPECIMEN" & pd.out$place_name=="Santa Monica Mountains","PD"]),2)), parse = TRUE, size = 2.5) +
  annotate("text", x = -0.55, y = 200, label = paste0("OTU==", as.numeric(pd.out[pd.out$basisofrecord=="PRESERVED_SPECIMEN" & pd.out$place_name=="Santa Monica Mountains","SR"])), parse = TRUE, size = 2.5) +
  annotate("text", x = -0.55, y = 250, label = paste0("lambda==", as.numeric(plam.out[plam.out$basisofrecord=="PRESERVED_SPECIMEN" & pd.out$place_name=="Santa Monica Mountains","plambda"])), parse = TRUE, size = 2.5) +
  annotate("text", x = -0.55, y = 300, label = paste0("p==", as.numeric(plam.out[plam.out$basisofrecord=="PRESERVED_SPECIMEN" & pd.out$place_name=="Santa Monica Mountains","pval"])), parse = TRUE, size = 2.5) 


# iNat + Herbarium records on smm
smmTree_comb <- smmTree + geom_fruit(data = smm.dat.comb,
                                     geom = geom_col,
                                     mapping = aes(y = OTU, x = sqrt(n)),
                                     fill = 'mediumpurple4',
                                     pwidth = 1.45,
) +
  theme(plot.margin = margin(.1,.1,.1,.1)) +
  annotate("text", x = -0.55, y = 150, label = paste0("PD==", round(as.numeric(pd.out[pd.out$basisofrecord=="Santa Monica Mountains" & pd.out$place_name=="Santa Monica Mountains","PD"]),2)), parse = TRUE, size = 2.5) +
  annotate("text", x = -0.55, y = 200, label = paste0("OTU==", as.numeric(pd.out[pd.out$basisofrecord=="Santa Monica Mountains" & pd.out$place_name=="Santa Monica Mountains","SR"])), parse = TRUE, size = 2.5) +
  annotate("text", x = -0.55, y = 250, label = paste0("lambda==", as.numeric(plam.out[plam.out$basisofrecord=="combined" & pd.out$place_name=="Santa Monica Mountains","plambda"])), parse = TRUE, size = 2.5) +
  annotate("text", x = -0.55, y = 300, label = paste0("p==", as.numeric(plam.out[plam.out$basisofrecord=="combined" & pd.out$place_name=="Santa Monica Mountains","pval"])), parse = TRUE, size = 2.5) 


# Combine the trees
require(cowplot)

#------ Combine
m1 <- ggplotGrob(mmTree_ps) 
m2 <- ggplotGrob(mmTree_obs) 
m3 <- ggplotGrob(mmTree_comb)

t1 <- ggplotGrob(tamTree_ps) 
t2 <- ggplotGrob(tamTree_obs) 
t3 <- ggplotGrob(tamTree_comb)

s1 <- ggplotGrob(smmTree_ps) 
s2 <- ggplotGrob(smmTree_obs) 
s3 <- ggplotGrob(smmTree_comb)

combo <- cowplot::plot_grid(m1, t1, s1, 
                            m2, t2, s2, m3, t3, s3, ncol=3,  
                            labels=c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)'),
                            label_fontface = "plain",
                            label_size = 8, hjust = 0, label_x = 0.01, align = "hv")


ggsave("Figure3.pdf", combo, dpi = 600, width=15, height=15, units="in")





