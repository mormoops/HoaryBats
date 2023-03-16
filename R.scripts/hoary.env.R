############################################################################################################
## Hoary Bat niche dynamics analysis - ecospat
## Aim: test for differences in niche overlap, equivalency, & similarity
## Author: Angelo Soto-Centeno
## Paper: phenotypic convergence of Hoary Bats - Soto-Centeno and Simmons 2022 Scientific Reports
############################################################################################################

# library
library(ecospat) 
library(ade4)

# Load data for the niche dynamics analysis (i.e colonized range vs native range)
haw <- read.csv(file = "PATH/TO/FILE/haw.csv", header = T) # e.g. hawaii range
nam <- read.csv(file = "PATH/TO/FILE/nam.csv", header = T) # e.g. N Am range
sam <- read.csv(file = "PATH/TO/FILE/sam.csv", header = T) # e.g. S Am range

# Climate variables chosen
  # see below on Pearson Correlation, Regressions, and VIFs
# "bio4"        "bio5"        "bio15"       "bio16"       "bio17"      "bioalt"


# Niche dynamics quantification and comparison
### PCA-ENVIRONMENT
### --- v each section is repetitive for each pairwise comparison --- ###
  # the pca is calibrated on all the sites of the study area
  # Calibrating the PCA in the whole study area, including native and colonized ranges (same as PCAenv in Broenniman et al. 2012)
  # this points out to only the environmental variables [, 3:9] in the order: native range, colonize range
### test order:
  # 1. South America (native) vs North America (colonized)
  # 2. North America (native) vs Hawaii (colonized)
  # 3. South America (native) vs Hawaii (colonized)


#####################################################
## 1. S America (native) vs N America (colonized)
#####################################################

# a. create PCA-environment
pca.env.sn <- dudi.pca(rbind(sam, nam)[, 3:8], center = T, scale = T, scannf = F, nf = 2)
    # plot variable contribution
ecospat.plot.contrib(contrib = pca.env.sn$co, eigen = pca.env.sn$eig)

# b. predict PCA scores
  # PCA scores for the whole study area
scores.globclim.sn <- pca.env.sn$li
  # PCA scores for the species native distribution # 10 indicates the species_occ column; 3:9 is the range of environmental variables
scores.sp.sam <- suprow(pca.env.sn, sam[which(sam[, 9] == 1), 3:8])$li  
  # PCA scores for the species colonized distribution # 10 indicates the species_occ column; 3:9 is the range of environmental variables
scores.sp.nam <- suprow(pca.env.sn, nam[which(nam[, 9] == 1), 3:8])$li  
  # PCA scores for the whole native study area
scores.clim.sam <- suprow(pca.env.sn, sam[, 3:8])$li  
  # PCA scores for the whole colonized study area
scores.clim.nam <- suprow(pca.env.sn, nam[, 3:8])$li  

# c. calculate the occurrence densities
  # grid for a species in the native range (e.g. South America)
grid.clim.sam <- ecospat.grid.clim.dyn(glob = scores.globclim.sn, glob1 = scores.clim.sam, 
                                       sp = scores.sp.sam, R = 100, th.sp = 0)
# grid for a species in the colonized range (e.g. North America)
grid.clim.nam <- ecospat.grid.clim.dyn(glob = scores.globclim.sn, glob1 = scores.clim.nam, 
                                       sp = scores.sp.nam, R = 100, th.sp = 0)

# d. calculate the index of niche overlap
  # compute Schoener's D
ecospat.niche.overlap(grid.clim.sam, grid.clim.nam, cor = T)$D
    # D = 0.4450455

# e. perform the Niche Equivalency Test with ecospat.niche.equivalency.test()
  # NOTE: 1000 replicates recommended
eq.test.sn <- ecospat.niche.equivalency.test(grid.clim.sam, grid.clim.nam, rep = 1000,
                                             intersection = 0.1, 
                                             overlap.alternative = "higher",
                                             expansion.alternative = "lower",
                                             stability.alternative = "higher",
                                             unfilling.alternative = "lower")
    # plot equivalency test
ecospat.plot.overlap.test(eq.test.sn, "D", "Equivalency: S Am vs N Am")

# f. peroform the Niche Similarity Test with ecospat.niche.similarity.test()
sim.test.sn <- ecospat.niche.similarity.test(grid.clim.sam, grid.clim.nam, rep = 1000, 
                                             overlap.alternative = "higher",
                                             expansion.alternative = "lower",
                                             stability.alternative = "higher",
                                             unfilling.alternative = "lower",
                                             intersection = 0.1, rand.type = 2)
    # Plot Similarity test
ecospat.plot.overlap.test(sim.test.sn, "D", "Similarity S Am vs N Am")

# g. delimiting niche categories and quantifying niche dynamics in analogue climates with ecospat.niche.dyn.index()
ecospat.niche.dyn.index(grid.clim.sam, grid.clim.nam, intersection = 0.1)$dynamic.index.w
  # print
niche.dyn.sn$dynamic.index.w
# expansion stability unfilling S Am vs N Am
# 0.3004385 0.6995615 0.2583219

# h. Visualizing niche categories, niche dynamics and climate analogy between ranges with ecospat.plot.niche.dyn()
ecospat.plot.niche.dyn(grid.clim.sam, grid.clim.nam, quant = 0.1, interest = 2,
                       title = "Niche Overlap S. America vs N. America", name.axis1 = "PC1", name.axis2 = "PC2")
  # add the niche centroids & arrows  
    # NOTE: only plot vector for one centroid - Native to Colonized range
ecospat.shift.centroids(scores.clim.sam, scores.clim.nam, col = "white")


# i. Plot Similarity test for niche expansion, stability and unfilling
ecospat.plot.overlap.test(sim.test.sn, "expansion", "Niche Similarity S. America vs N. America") # test for niche expansion
ecospat.plot.overlap.test(sim.test.sn, "stability", "Niche Similarity S. America vs N. America") # test for niche expansion
ecospat.plot.overlap.test(sim.test.sn, "unfilling", "Niche Similarity S. America vs N. America") # test for niche expansion

#         expansion stability unfilling S Am vs N Am
#         0.3004385 0.6995615 0.2583219
# pvalue  0.16284   0.16284   0.05195

# j. make grid of native and colonized niches for a single variable
grid.clim.t.sam <- ecospat.grid.clim.dyn(glob = as.data.frame(rbind(sam, nam)[ , 8]),
                                         glob1 = as.data.frame(sam[ , 8]),
                                         sp = as.data.frame(sam[which(sam[ , 9] == 1), 8]),
                                         R = 1000, th.sp = 0)
  # gridding the colonized niche
grid.clim.t.nam <- ecospat.grid.clim.dyn(glob = as.data.frame(rbind(sam, nam)[ , 8]),
                                         glob1 = as.data.frame(nam[ , 8]),
                                         sp = as.data.frame(nam[which(nam[ , 9] == 1), 8]),
                                         R = 1000, th.sp = 0)
  # create the niche dynamic index for that single variable
t.dyn <- ecospat.niche.dyn.index(grid.clim.t.sam, grid.clim.t.nam,
                                intersection = 0.5)
  # plot the niche dynamic index for altitude
ecospat.plot.niche.dyn(grid.clim.t.sam, grid.clim.t.nam, quant = 0,
                       interest = 2, title = "Niche Overlap",
                       name.axis1 = "Elevation (bioalt)")

#####################################################
## 2. N America (native) vs Hawaii (colonized)
#####################################################

# repeat code as in the first comparison


#####################################################
## 3. S America (native) vs Hawaii (colonized)
#####################################################

# repeat code as in the first comparison



### --- clear up your Global Environment --- ###
rm(list = ls())
### --- clear up plots --- ###
dev.off()



####################################################################
# Run multiple linear regression to test the relationship of phenotype (PC1) vs environment
# Examine variable co-linearity using Variance Inflation Factors
# NOTE: explore environmental variable co-linearity
# DATE: 30 Nov 2022
####################################################################


# load libraries
library(raster)
library(tidyverse)
library(car)


# upload mcl.data.csv to run the regression analyses (pre-wrangled data after correlation)
data <- read.csv(file = "PATH/TO/DIRECTORY/mcl.data.csv", header = T)

# Environmental variables chosen
# "bio4"        "bio5"        "bio15"       "bio16"       "bio17"      "bioalt"


## Full multiple linear regression with all variables
full.model <- lm(PC1 ~ bio4 + bio5 + bio15 + bio16 + bio17 + bioalt + lon + lat, data = ddata)
summary(full.model)

# calculate the VIF for each predictor variable in the model
vif(full.model)

# make vector of VIF values
vif_values <- vif(full.model)

# plot bar chart
barplot(vif_values, main = "VIF values of full model", horiz = T, col = "steelblue", xlim = c(0, 6))

# add vertical line at 5
abline(v = 5, lwd = 3, lty = 2)

## NOTE: variables VIF > 5 are considered highly correlated and must be removed from the model.

# Reduced model excluding those variables with VIF > 5
red.model <- lm(PC1 ~ bio4 + bio5 + bio16 + bioalt + lon + lat, data = ddata)
summary(red.model)

# calculate the VIF for each predictor variable in the model
vif(red.model)

# make vector of VIF values
vif_values.red <- vif(red.model)

# plot bar chart
barplot(vif_values.red, main = "VIF values of reduced model", horiz = T, col = "steelblue", xlim = c(0, 6))

# add vertical line at 5
abline(v = 5, lwd = 3, lty = 2)

## regression plots
# plot the regression PC1 vs lat
  # create color scale
    # colors for cinereus/semotus/villosissimus
colrs.csv <- c("#000000", "#4292c6", "#993404") # black / blue / brown
  # set the plot theme
theme_set(theme_light())

  # PC1 vs latitude
Fig.lat <- ggplot(cdata, aes(x = lat, y = PC1)) +
  stat_smooth(method = "lm", formula = y ~ x, se = F) +
  geom_point(aes(colour = Ssp), size = 5, position = position_jitter(width = 0.2), alpha = 0.5) +
  labs(x = "Latitude", y = "Phenotypic Variation (PC1)") +
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18),
        legend.position = "none") # change legend.position to "right" so it can be displayed
# add the custom color palete
Fig.lat + scale_color_manual(values = colrs.csv) # or you can set colors manually by using c()

## repeat the code above for each comparison
  # PC1 vs longitude
  # PC1 vs elevation
  # PC1 vs bio16
  # PC1 vs bio4
  # PC1 vs bio4

#### ---- FIN ---- ####