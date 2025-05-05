# This script performs statistical analysis on the distances between species
# occurrences and their nearest roads, using Gamma regression models. The analysis
# compares distances between iNaturalist observations and herbarium specimens for each
# study area, including model diagnostics and confidence intervals.

require(dplyr)
require(ggplot2)
require(tidyverse)
require(ggExtra)
require(cowplot)
require(grid)
require(gridExtra)


# Load data
rd <- read.csv("data/occ_nearest_trail_2025-04-22-2.csv")

head(rd)
table(rd$place_name)

#----------- Tam
tm.rd <- rd %>% filter(place_name == "One Tam")
hist(tm.rd$distance)
summary(tm.rd)

mod.1 <- glm(distance ~ basisofrecord, data = tm.rd, family = Gamma(link = "log"))
summary(mod.1)

# The second value "basisofrecordPRESERVED_SPECIMEN" is the parameter estimate
mod.1$coefficients
# And this gives the 95% CI around the parameter estimate
confint(mod.1)
qqnorm(resid(mod.1), pch = 16, main = "Gamma link=log")
qqline(resid(mod.1))

# Check the other ones, Gamma log link looks best
mod.gaus <- lm(distance ~ basisofrecord, data = tm.rd)
qqnorm(resid(mod.gaus), pch = 16, main = "Gaussian")
qqline(resid(mod.gaus))
tm.rd$log_dist <- log(tm.rd$distance)
mod.log <- lm(log_dist ~ basisofrecord, data = tm.rd)
qqnorm(resid(mod.log), pch = 16, main = "Gaussian log transformed")
qqline(resid(mod.log))

#---- MM
mm.rd <- rd %>% filter(place_name == "Marble/Salmon Mountains")
hist(mm.rd$distance)
summary(mm.rd)

mod.2 <- glm(distance ~ basisofrecord, data = mm.rd, family = Gamma(link = "log"))
summary(mod.2)

# The second value "basisofrecordPRESERVED_SPECIMEN" is the parameter estimate
mod.2$coefficients
# And this gives the 95% CI around the parameter estimate
confint(mod.2)
qqnorm(resid(mod.2), pch = 16, main = "Gamma link=log")
qqline(resid(mod.2))


#---- Santa MM
smm.rd <- rd %>% filter(place_name == "Santa Monica Mountains")
hist(smm.rd$distance)
summary(smm.rd)

mod.3 <- glm(distance ~ basisofrecord, data = smm.rd, family = Gamma(link = "log"))
summary(mod.3)

# The second value "basisofrecordPRESERVED_SPECIMEN" is the parameter estimate
mod.3$coefficients
# And this gives the 95% CI around the parameter estimate
confint(mod.3)
qqnorm(resid(mod.3), pch = 16, main = "Gamma link=log")
qqline(resid(mod.3))
