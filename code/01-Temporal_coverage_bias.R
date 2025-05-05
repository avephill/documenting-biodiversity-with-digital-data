#---------------------
# TEMMPORAL BIAS
# Analyses and Plots
#
# I. Species accumulation through time (Figure 4a-c)
# II. Native and non native species accumulation (Figure 4d-f)
# III. Species richness through time (Figure 4g-i)
# IV. Change in the number of records through time (Figure 5a-f)
# V. Change in the number of collectors/observers through time (Figure 5g-l)
#---------------------
require(dplyr)
require(ggplot2)
require(tidyverse)
require(ggExtra)
require(cowplot)
require(grid)
require(gridExtra)
require(mgcv)

# Set working directory
setwd()

#-----------------------
# LOAD IN DATA
dat <- read.csv("full_gbif_for_analysis.csv")
table(dat$place_name)

dat$place_name[dat$place_name %in% "Marble/Salmon Mountains"] <- "Marble Mountains"
dat$place_name[dat$place_name %in% "One Tam"] <- "Mount Tamalpais"

dat %>% filter(verifiedTaxon=="verified") %>% filter(keepWhenDeduplicating=="keep") -> dats

#------------------------------
# TEMMPORAL BIAS
# I. Species accumulation through time (Figure 4a-c)

# How many samples don't have a year associated with them?
dats %>% filter(is.na(year)) %>% nrow() /nrow(dats) 
# Remove from temporal analysis

# I. Species accumulation through time (Figure 4a-c)
# What year was the first records for each place and sample type  
dats %>% group_by(place_name, basisofrecord) %>% summarise(min(year, na.rm = TRUE)) %>%
  rename(firstYear= `min(year, na.rm = TRUE)`) %>% ungroup() %>%
  spread(., basisofrecord, firstYear) -> firstYear

# For each data type at each location, add the start and end years of data collection
# to account for years without data  
data.frame(place_name="Mount Tamalpais", basisofrecord="HUMAN_OBSERVATION", year = seq(as.numeric(firstYear[firstYear$place_name=="Mount Tamalpais","HUMAN_OBSERVATION"]), 2024, by=1)) %>% 
  rbind(data.frame(place_name="Mount Tamalpais", basisofrecord="PRESERVED_SPECIMEN", year = seq(as.numeric(firstYear[firstYear$place_name=="Mount Tamalpais","PRESERVED_SPECIMEN"]), 2024, by=1))) %>%
  rbind(data.frame(place_name="Marble Mountains", basisofrecord="HUMAN_OBSERVATION", year = seq(as.numeric(firstYear[firstYear$place_name=="Marble Mountains","HUMAN_OBSERVATION"]), 2024, by=1))) %>% 
  rbind(data.frame(place_name="Marble Mountains", basisofrecord="PRESERVED_SPECIMEN", year = seq(as.numeric(firstYear[firstYear$place_name=="Marble Mountains","PRESERVED_SPECIMEN"]), 2024, by=1))) %>%
  rbind(data.frame(place_name="Santa Monica Mountains", basisofrecord="HUMAN_OBSERVATION", year = seq(as.numeric(firstYear[firstYear$place_name=="Santa Monica Mountains","HUMAN_OBSERVATION"]), 2024, by=1))) %>% 
  rbind(data.frame(place_name="Santa Monica Mountains", basisofrecord="PRESERVED_SPECIMEN", year = seq(as.numeric(firstYear[firstYear$place_name=="Santa Monica Mountains","PRESERVED_SPECIMEN"]), 2024, by=1))) %>%
  rename(first.obs=year) -> syears

# Calculate the cumulative species diversity -- for both data types through time
dats %>% select(year, place_name, species) %>%
  filter(!is.na(year)) %>%
  group_by(place_name,species) %>% summarise(first.obs=min(year)) %>% 
  group_by(place_name, first.obs) %>% count() %>% ungroup() %>% 
  arrange(first.obs) %>% 
  group_by(place_name)  %>% mutate(cum_n = cumsum(n)) %>% 
  mutate(basisofrecord="combined") %>%
  select(place_name, basisofrecord, first.obs, n, cum_n) -> togetherSUM

# Create plots with species accumulation
dats %>% select(year, place_name, basisofrecord, species) %>%
  filter(!is.na(year)) %>%
  group_by(place_name, basisofrecord, species) %>% summarise(first.obs=min(year)) %>%
  group_by(place_name, basisofrecord, first.obs) %>% count() %>% ungroup() %>% 
  full_join(syears, by=c("place_name", "basisofrecord", "first.obs")) %>% 
  arrange(first.obs) %>% mutate(n =replace_na(n, 0)) %>%
  group_by(place_name, basisofrecord) %>% mutate(cum_n = cumsum(n)) %>% 
  rbind(togetherSUM) %>% select(-n) %>% 
  # When the diversity of combined becomes greater than the diversity of Preserved specimens alone 
  #spread(., basisofrecord, cum_n) %>%  filter(!is.na(combined)) %>% filter(combined==PRESERVED_SPECIMEN) %>%
  #ungroup() %>% group_by(place_name) %>% summarize(max(first.obs))
  # Adding in Cumulative records for each location // the year when the cumulative SR > herbarium SR
  mutate(cum_n = if_else(place_name=="Marble Mountains" & basisofrecord=="combined" & first.obs<=2011, NA, cum_n)) %>%
  mutate(cum_n = if_else(place_name=="Mount Tamalpais" & basisofrecord=="combined" & first.obs<=1983, NA, cum_n)) %>%
  mutate(cum_n = if_else(place_name=="Santa Monica Mountains" & basisofrecord=="combined" & first.obs<=1987, NA, cum_n)) %>%
  mutate(place_name = factor(place_name, levels=c("Marble Mountains", "Mount Tamalpais",  "Santa Monica Mountains"))) %>%
  mutate(basisofrecord = factor(basisofrecord, 
                                levels=c("HUMAN_OBSERVATION", "PRESERVED_SPECIMEN", "combined"))) %>%
  ggplot(aes(x=first.obs, y=cum_n, color=basisofrecord)) + #geom_point(shape=21, size=1.5, color="white") + 
  geom_point() + xlab("Year") + 
  ylab("Species richness (accumulation)") + theme_bw() + theme(legend.position="none") + 
  scale_fill_manual(values=c("combined"="mediumpurple1", "HUMAN_OBSERVATION"="darkgoldenrod1", "PRESERVED_SPECIMEN"="darkslategray")) +
  scale_color_manual(values=c("combined"="mediumpurple1", "HUMAN_OBSERVATION"="darkgoldenrod1", "PRESERVED_SPECIMEN"="darkslategray")) +
  facet_wrap(~place_name) -> site.specAcc
# Warnings the 281 are all of the NAs for the combined records column 


#------------------- 
# II. Native and non native species accumulation (Figure 4d-f)

# Based on the start of each data type, create a dumby data frame to fill in years when there may not have been data collected
dats %>% group_by(place_name, basisofrecord, Origin) %>% summarise(min(year, na.rm = TRUE)) %>%
  rename(firstYear= `min(year, na.rm = TRUE)`) %>% ungroup() %>% 
  mutate(type=paste(basisofrecord, Origin, sep = "_")) %>% select(-Origin, -basisofrecord) %>%
  spread(., type, firstYear) -> firstYearOrigin

data.frame(place_name="Mount Tamalpais", basisofrecord="HUMAN_OBSERVATION", Origin="native", year = seq(as.numeric(firstYearOrigin[firstYearOrigin$place_name=="Mount Tamalpais","HUMAN_OBSERVATION_native"]), 2024, by=1)) %>% 
  rbind(data.frame(place_name="Mount Tamalpais", basisofrecord="HUMAN_OBSERVATION", Origin="naturalized", year = seq(as.numeric(firstYearOrigin[firstYearOrigin$place_name=="Mount Tamalpais","HUMAN_OBSERVATION_naturalized"]), 2024, by=1))) %>% 
  rbind(data.frame(place_name="Mount Tamalpais", basisofrecord="PRESERVED_SPECIMEN", Origin="native", year = seq(as.numeric(firstYearOrigin[firstYearOrigin$place_name=="Mount Tamalpais","PRESERVED_SPECIMEN_native"]), 2024, by=1))) %>%
  rbind(data.frame(place_name="Mount Tamalpais", basisofrecord="PRESERVED_SPECIMEN", Origin="naturalized", year = seq(as.numeric(firstYearOrigin[firstYearOrigin$place_name=="Mount Tamalpais","PRESERVED_SPECIMEN_naturalized"]), 2024, by=1))) %>%
  rbind(data.frame(place_name="Marble Mountains", basisofrecord="HUMAN_OBSERVATION", Origin="native", year = seq(as.numeric(firstYearOrigin[firstYearOrigin$place_name=="Marble Mountains","HUMAN_OBSERVATION_native"]), 2024, by=1))) %>% 
  rbind(data.frame(place_name="Marble Mountains", basisofrecord="HUMAN_OBSERVATION", Origin="naturalized", year = seq(as.numeric(firstYearOrigin[firstYearOrigin$place_name=="Marble Mountains","HUMAN_OBSERVATION_naturalized"]), 2024, by=1))) %>% 
  rbind(data.frame(place_name="Marble Mountains", basisofrecord="PRESERVED_SPECIMEN", Origin="native", year = seq(as.numeric(firstYearOrigin[firstYearOrigin$place_name=="Marble Mountains","PRESERVED_SPECIMEN_native"]), 2024, by=1))) %>%
  rbind(data.frame(place_name="Marble Mountains", basisofrecord="PRESERVED_SPECIMEN", Origin="naturalized", year = seq(as.numeric(firstYearOrigin[firstYearOrigin$place_name=="Marble Mountains","PRESERVED_SPECIMEN_naturalized"]), 2024, by=1))) %>%
  rbind(data.frame(place_name="Santa Monica Mountains", basisofrecord="HUMAN_OBSERVATION", Origin="native", year = seq(as.numeric(firstYearOrigin[firstYearOrigin$place_name=="Santa Monica Mountains","HUMAN_OBSERVATION_native"]), 2024, by=1))) %>% 
  rbind(data.frame(place_name="Santa Monica Mountains", basisofrecord="HUMAN_OBSERVATION", Origin="naturalized", year = seq(as.numeric(firstYearOrigin[firstYearOrigin$place_name=="Santa Monica Mountains","HUMAN_OBSERVATION_naturalized"]), 2024, by=1))) %>% 
  rbind(data.frame(place_name="Santa Monica Mountains", basisofrecord="PRESERVED_SPECIMEN", Origin="native", year = seq(as.numeric(firstYearOrigin[firstYearOrigin$place_name=="Santa Monica Mountains","PRESERVED_SPECIMEN_native"]), 2024, by=1))) %>%
  rbind(data.frame(place_name="Santa Monica Mountains", basisofrecord="PRESERVED_SPECIMEN", Origin="naturalized", year = seq(as.numeric(firstYearOrigin[firstYearOrigin$place_name=="Santa Monica Mountains","PRESERVED_SPECIMEN_naturalized"]), 2024, by=1))) %>%
  rename(first.obs=year) -> syears.nn


# Plot
dats %>% select(year, place_name, Origin, basisofrecord, species) %>%
  filter(!is.na(year)) %>%
  group_by(place_name, basisofrecord, Origin, species) %>% summarise(first.obs=min(year)) %>% 
  group_by(place_name, basisofrecord, Origin, first.obs) %>% count() %>% ungroup() %>% 
  full_join(syears.nn, by=c("place_name", "basisofrecord", "Origin", "first.obs")) %>% 
  mutate(originBasis = paste(basisofrecord, Origin, sep="_")) %>% 
  arrange(first.obs) %>% mutate(n =replace_na(n, 0)) %>%
  group_by(place_name, originBasis)  %>% mutate(cum_n = cumsum(n)) %>% 
  mutate(place_name = factor(place_name, levels=c( "Marble Mountains", "Mount Tamalpais", "Santa Monica Mountains"))) %>% 
  ggplot(aes(x=first.obs, y=cum_n, color=originBasis)) + geom_point(size=1) + xlab("Year") + 
  ylab("Species richness (accumulation)") +  theme_bw() + theme(legend.position="none") + 
  scale_color_manual(values=c("HUMAN_OBSERVATION_native"="darkgoldenrod1", "HUMAN_OBSERVATION_naturalized"="#DC863B",
                             # "combined_native"="mediumpurple4", "combined_naturalized"="mediumpurple1",
                              "PRESERVED_SPECIMEN_native"="darkslategray",  "PRESERVED_SPECIMEN_naturalized"="#A2A475")) +
  facet_wrap(~place_name) -> native.non.specAcc

#------------
# III. Species richness through time (Figure 4g-i)

dats %>% filter(!is.na(year)) %>%
  # Removed 2024 because we didn't have complete sampling effort for that year
  filter(year<2024) %>%
  select(place_name, basisofrecord, year, species) %>% distinct() %>%
  group_by(place_name, basisofrecord, year) %>% count() %>% 
  mutate(place_name = factor(place_name, levels=c( "Marble Mountains", "Mount Tamalpais", "Santa Monica Mountains"))) %>% 
  # when does iNat richness surpass herbarium?
  #  spread(., basisofrecord, n) %>% filter(!is.na(HUMAN_OBSERVATION)) %>% filter(!is.na(PRESERVED_SPECIMEN)) %>%
  # filter(HUMAN_OBSERVATION>PRESERVED_SPECIMEN) %>% View() 
  ggplot() + 
  geom_point(aes(x=year, y=n, color=basisofrecord), alpha=0.7 ) + #, alpha=0.7 
  geom_line(aes(x=year, y=n, color=basisofrecord), alpha=0.7 ) + 
  ylab("alpha diversity") + theme_bw() + theme(legend.position="none") + 
  scale_color_manual(values=c("HUMAN_OBSERVATION"="darkgoldenrod1", "PRESERVED_SPECIMEN"="darkslategray")) +
  facet_wrap( ~place_name) -> alpha.time


#------ Combine
a1 <- ggplotGrob(site.specAcc) 
a2 <- ggplotGrob(native.non.specAcc) 
# From below
a3 <- ggplotGrob(alpha.time)

combo <- cowplot::plot_grid(a1, a2, a3, ncol=1,  
                            labels=c('(a)', '(b)', '(c)', '(d)'),
                            label_fontface = "plain",
                            label_size = 8, hjust = 0, label_x = 0.01, align = "h")

ggsave("Figure_4.pdf", combo, dpi = 600, width=8, height=7, units="in")

#---------------
# IV. Change in the number of records through time (Figure 5a-f)

#------------------------
# TAM
t.tm <- dats %>% filter(!is.na(year)) %>% group_by(place_name, basisofrecord, year) %>% count() %>% ungroup() %>% 
  filter(year<2024) %>% filter(place_name=="Mount Tamalpais") %>% mutate(log.n = log(n), fbasisofrecord=as.factor(basisofrecord))
summary(t.tm)

# Model fitting - Count data
m.tam.p <- gam(n ~ s(year, by=fbasisofrecord, bs = "cr"), family=poisson, data=t.tm)

# Checking for overdisispersion
E1 <- resid(m.tam.p, type = "pearson")
Overdispersion <- sum(E1^2) / (m.tam.p$df.res)
Overdispersion # Very high

m.tam.nb <- gam(n ~ s(year, by=fbasisofrecord, bs = "cr"), family=nb(), data=t.tm)
qqnorm(resid(m.tam.nb), pch=16, main="Negative binomial GAM"); qqline(resid(m.tam.nb)) 

E1 <- resid(m.tam.nb, type = "pearson")
Overdispersion <- sum(E1^2) / (m.tam.nb$df.res)
Overdispersion # good

# Extract e.d.f. and p-value
summary(m.tam.nb)

#--- Predict
newdat.tam.r <- with(t.tm, expand.grid(year=seq(min(year), max(year), by=1),
                                     fbasisofrecord = unique(fbasisofrecord)))
newdat.tam.r$fit <- predict(m.tam.nb, newdata = newdat.tam.r, type="link", se = TRUE)$fit
newdat.tam.r$se.fit <- predict(m.tam.nb, newdata = newdat.tam.r, type="link", se = TRUE)$se.fit
newdat.tam.r$place_name <-"Mount Tamalpais"
# 
min.ho <- t.tm %>% filter(fbasisofrecord=="HUMAN_OBSERVATION") %>% summarise(min(year))
newdat.tam.r <- newdat.tam.r %>% filter(!(fbasisofrecord=="HUMAN_OBSERVATION" & year < (min.ho$`min(year)`)))

#------------------------
# Santa Monica mountains
t.sm <- dats %>% filter(!is.na(year)) %>% group_by(place_name, basisofrecord, year) %>% count() %>% ungroup() %>% 
  filter(year<2024) %>% filter(place_name=="Santa Monica Mountains") %>% mutate(log.n = log(n), fbasisofrecord=as.factor(basisofrecord))
summary(t.sm)

m.sm.nb <- gam(n ~ s(year, by=fbasisofrecord, bs = "cr"), family=nb(), data=t.sm)
qqnorm(resid(m.sm.nb), pch=16, main="Negative binomial GAM"); qqline(resid(m.sm.nb)) 

# Extract e.d.f. and p-value
summary(m.sm.nb)

#--- Predict
newdat.sm.r <- with(t.sm, expand.grid(year=seq(min(year), max(year), by=1),
                                    fbasisofrecord = unique(fbasisofrecord)))
newdat.sm.r$fit <- predict(m.sm.nb, newdata = newdat.sm.r, type="link",  se = TRUE)$fit
newdat.sm.r$se.fit <- predict(m.sm.nb, newdata = newdat.sm.r, type="link",  se = TRUE)$se.fit
newdat.sm.r$place_name <-"Santa Monica Mountains"
min.ho <- t.sm %>% filter(fbasisofrecord=="HUMAN_OBSERVATION") %>% summarise(min(year))
newdat.sm.r <- newdat.sm.r %>% filter(!(fbasisofrecord=="HUMAN_OBSERVATION" & year < (min.ho$`min(year)`)))

#------------------------
# Marble Mountains
t.mm <- dats %>% filter(!is.na(year)) %>% group_by(place_name, basisofrecord, year) %>% count() %>% ungroup() %>% 
  filter(year<2024) %>% filter(place_name=="Marble Mountains") %>% mutate(log.n = log(n), fbasisofrecord=as.factor(basisofrecord))
summary(t.mm)

m.mm.nb <- gam(n ~ s(year, by=fbasisofrecord, bs = "cr"), family=nb(), data=t.mm)
qqnorm(resid(m.mm.nb), pch=16, main="NB GLM"); qqline(resid(m.mm.nb)) 

# Extract e.d.f. and p-value
summary(m.mm.nb)

#--- Predict
newdat.mm.r <- with(t.mm, expand.grid(year=seq(min(year), max(year), by=1),
                                    fbasisofrecord = unique(fbasisofrecord)))
newdat.mm.r$fit <- predict(m.mm.nb, newdata = newdat.mm.r, type="link", se = TRUE)$fit
newdat.mm.r$se.fit <- predict(m.mm.nb, newdata = newdat.mm.r, type="link", se = TRUE)$se.fit
newdat.mm.r$place_name <-"Marble Mountains"
min.ho <- t.mm %>% filter(fbasisofrecord=="HUMAN_OBSERVATION") %>% summarise(min(year))
newdat.mm.r <- newdat.mm.r %>% filter(!(fbasisofrecord=="HUMAN_OBSERVATION" & year < (min.ho$`min(year)`)))

#--- Combine data
newdat.REC <- rbind(newdat.mm.r, newdat.sm.r, newdat.tam.r) %>%
  mutate(est = exp(fit), low=exp(fit - 1.96*se.fit), high=exp(fit+1.96*se.fit), 
         basisofrecord=as.character(fbasisofrecord))

summary(newdat.REC)

#------ Plot the number of records through time
dats %>% filter(!is.na(year)) %>% group_by(place_name, basisofrecord, year) %>% count() %>% ungroup() %>% 
  mutate(place_name = factor(place_name, levels=c("Marble Mountains", "Mount Tamalpais",  "Santa Monica Mountains"))) %>%
  ggplot(aes(group=basisofrecord, y=n, x=year)) + 
  geom_point( size=1, color="gray75") + theme_bw() +
  geom_ribbon(data=newdat.REC, aes(y=year, ymin=low, ymax=high, fill=basisofrecord),linetype = 0,alpha=0.75) +
  geom_line(data=newdat.REC, aes(y=est, x=year, group=basisofrecord, color=basisofrecord)) +
  ylab("Count of records") +  theme(legend.position="none") + 
  scale_color_manual(values=c("HUMAN_OBSERVATION"="darkgoldenrod1", "PRESERVED_SPECIMEN"="darkslategray")) +
  scale_fill_manual(values=c("HUMAN_OBSERVATION"="gold1", "PRESERVED_SPECIMEN"="darkslategray4")) +
  facet_wrap(basisofrecord~place_name, scales = "free_y") -> effortTime.Records.plot

#--------------------------------
# V. Change in the number of collectors/observers through time (Figure 5g-l)

#------------------------
# TAM
o.tm <- dats %>% filter(!is.na(year)) %>% select(place_name, basisofrecord, recordedby, year) %>% distinct() %>%
  group_by(place_name, basisofrecord, year) %>% count() %>% ungroup() %>% 
  filter(year<2024) %>% filter(place_name=="Mount Tamalpais") %>% mutate(log.n = log(n), fbasisofrecord=as.factor(basisofrecord))

# Model fitting 
m.tam.obs.p <- gam(n ~ s(year, by=fbasisofrecord, bs = "cr"), family=poisson, data=o.tm)

# Checking for overdisispersion
E1 <- resid(m.tam.obs.p, type = "pearson")
Overdispersion <- sum(E1^2) / (m.tam.obs.p$df.res)
Overdispersion # High

m.tam.obs.nb <- gam(n ~ s(year, by=fbasisofrecord, bs = "cr"), family=nb(), data=o.tm)
qqnorm(resid(m.tam.obs.nb), pch=16, main="Negative binomial GAM"); qqline(resid(m.tam.obs.nb)) 

E1 <- resid(m.tam.obs.nb, type = "pearson")
Overdispersion <- sum(E1^2) / (m.tam.obs.nb$df.res)
Overdispersion # Better

# Extract e.d.f. and p-value
summary(m.tam.obs.nb)

#--- Predict
newdat.tam <- with(o.tm, expand.grid(year=seq(min(year), max(year), by=1),
                                     fbasisofrecord = unique(fbasisofrecord)))
newdat.tam$fit <- predict(m.tam.obs.nb, newdata = newdat.tam, type="link", se = TRUE)$fit
newdat.tam$se.fit <- predict(m.tam.obs.nb, newdata = newdat.tam, type="link", se = TRUE)$se.fit
newdat.tam$place_name <-"Mount Tamalpais"
min.ho <- t.tm %>% filter(fbasisofrecord=="HUMAN_OBSERVATION") %>% summarise(min(year))
newdat.tam <- newdat.tam %>% filter(!(fbasisofrecord=="HUMAN_OBSERVATION" & year < (min.ho$`min(year)`)))

#------------------------
# Santa Monica mountains
o.sm <- dats %>% filter(!is.na(year)) %>% select(place_name, basisofrecord, recordedby, year) %>% distinct() %>%
  group_by(place_name, basisofrecord, year) %>% count() %>% ungroup() %>%   
  filter(year<2024) %>% filter(place_name=="Santa Monica Mountains") %>% mutate(log.n = log(n), fbasisofrecord=as.factor(basisofrecord))
summary(t.sm)

m.sm.obs.nb <- gam(n ~ s(year, by=fbasisofrecord, bs = "cr"), family=nb(), data=o.sm)
qqnorm(resid(m.sm.obs.nb), pch=16, main="Gaussian GAM"); qqline(resid(m.sm.obs.nb)) 

# Extract e.d.f. and p-value
summary(m.sm.obs.nb)

#--- Predict
newdat.sm <- with(o.sm, expand.grid(year=seq(min(year), max(year), by=1),
                                    fbasisofrecord = unique(fbasisofrecord)))
newdat.sm$fit <- predict(m.sm.obs.nb, newdata = newdat.sm, type="link",  se = TRUE)$fit
newdat.sm$se.fit <- predict(m.sm.obs.nb, newdata = newdat.sm, type="link",  se = TRUE)$se.fit
newdat.sm$place_name <-"Santa Monica Mountains"
min.ho <- t.sm %>% filter(fbasisofrecord=="HUMAN_OBSERVATION") %>% summarise(min(year))
newdat.sm <- newdat.sm %>% filter(!(fbasisofrecord=="HUMAN_OBSERVATION" & year < (min.ho$`min(year)`)))

#------------------------
# Marble Mountains
o.mm <- dats %>% filter(!is.na(year)) %>% select(place_name, basisofrecord, recordedby, year) %>% distinct() %>%
  group_by(place_name, basisofrecord, year) %>% count() %>% ungroup() %>%   
  filter(year<2024) %>% filter(place_name=="Marble Mountains") %>% mutate(log.n = log(n), fbasisofrecord=as.factor(basisofrecord))
summary(o.mm)

m.mm.obs.nb <- gam(n ~ s(year, by=fbasisofrecord, bs = "cr"), family=nb(), data=o.mm)
qqnorm(resid(m.mm.obs.nb), pch=16, main="NB GLM"); qqline(resid(m.mm.obs.nb)) 

# Extract e.d.f. and p-value
summary(m.mm.obs.nb)

#--- Predict
newdat.mm <- with(o.mm, expand.grid(year=seq(min(year), max(year), by=1),
                                    fbasisofrecord = unique(fbasisofrecord)))
newdat.mm$fit <- predict(m.mm.obs.nb, newdata = newdat.mm, type="link", se = TRUE)$fit
newdat.mm$se.fit <- predict(m.mm.obs.nb, newdata = newdat.mm, type="link", se = TRUE)$se.fit
newdat.mm$place_name <-"Marble Mountains"
min.ho <- t.mm %>% filter(fbasisofrecord=="HUMAN_OBSERVATION") %>% summarise(min(year))
newdat.mm <- newdat.mm %>% filter(!(fbasisofrecord=="HUMAN_OBSERVATION" & year < (min.ho$`min(year)`)))

#----- Combine data
newdat.obs.REC <- rbind(newdat.mm, newdat.sm, newdat.tam) %>%
  mutate(est = exp(fit), low=exp(fit - 1.96*se.fit), high=exp(fit+1.96*se.fit), 
         basisofrecord=as.character(fbasisofrecord))

summary(newdat.obs.REC)

#---------------------
dats %>% filter(!is.na(year)) %>% 
  select(place_name, basisofrecord, year, recordedby) %>% distinct() %>%
  group_by(place_name, basisofrecord, year) %>% count() %>% ungroup() %>% 
  mutate(place_name = factor(place_name, levels=c("Marble Mountains", "Mount Tamalpais",  "Santa Monica Mountains"))) %>%
  ggplot(aes(group=basisofrecord, y=n, x=year)) + 
  geom_point( size=1, color="gray75") + theme_bw() +
  geom_ribbon(data=newdat.obs.REC, aes(y=year, ymin=low, ymax=high, fill=basisofrecord),linetype = 0,alpha=0.75) +
  geom_line(data=newdat.obs.REC, aes(y=est, x=year, group=basisofrecord, color=basisofrecord)) +
  ylab("Number of observers") +  theme(legend.position="none") + 
  scale_color_manual(values=c("HUMAN_OBSERVATION"="darkgoldenrod1", "PRESERVED_SPECIMEN"="darkslategray")) +
  scale_fill_manual(values=c("HUMAN_OBSERVATION"="gold1", "PRESERVED_SPECIMEN"="darkslategray4")) +
  facet_wrap(basisofrecord~place_name, scales = "free_y") -> ObsTime.plot


#------ Combine
a4 <- ggplotGrob(effortTime.Records.plot) 
a5 <- ggplotGrob(ObsTime.plot) 


combo <- cowplot::plot_grid(a4, a5, ncol=1,  
                            labels=c('(a)', '(b)', '(c)', '(d)'),
                            label_fontface = "plain",
                            label_size = 8, hjust = 0, label_x = 0.01, align = "h")


ggsave("Figure5.pdf", combo, dpi = 600, width=8, height=8, units="in")


