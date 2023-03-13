### --- Script to analyze phenotypic data for Hoary bat species --- ###
### --- Author: J.A. Soto-Centeno
### --- Last date modified: 1 Nov 2022 --- ###
### --- IN: Soto-Centeno & Simmons 2022 Scientific Reports

## This script runs analyses for species delimitation using morphological data collected by
## Soto-Centeno & Simmons 2022. The script implements data imputation based on the 
## Multivariate Imputation by Chained Equations (VanBurren et al. 2011), custom Linear
## Discriminant Analysis and Principal Component Analysis

## Main script sections:
## 1. Packages and dependencies
## 2. Preliminary data management & imputation (mice)
## 3. Data wrangling to create multiple partitions for LDA & PCA


# load packages
library(tidyverse)
library(caret)
library(MASS) # does the lda
library(car) # does scatterplot matrix
library(gridExtra)
library(kernlab)
# data imputation 
library(mice)
library(ade4)

## 1. load clean data & create
data_all <- read.csv(file = "PATH/TO/*.csv", header = T)


## 2. wrangle a phenotypic partition
  # A. morphological comparison of the current species of hoary bats: cinereus vs semotus vs villosissimus
data <- data_all %>% filter(Ssp == "cinereus" | Ssp == "semotus" | Ssp == "villosissimus")

  # B. morphological comparison of the recognized sub-species of Lasiurus villosissimus
    # L.v. villosissimus - Peru, Bolivia, Paraguay, Brazil, Uruguay, Argentina
    # L.v. brasiliensis - SE Brazil -- Minas Gerais & Sao Paulo
    # L.v. grayi - Chile
    # L.v. palescens - Colombia and Venezuela
    # L.v. "galapagos" - unknown subspecies
vil <- subset(data, Country == "Argentina" | Country == "Uruguay" | Country == "Chile" | Country == "Peru" |
              Country == "Colombia" | Country == "Paraguay" | Country == "Venezuela" | Country == "Galapagos")

  # C. morphological comparison of Lasiurus cinereus based on geography
    # L. cinereus from Canada, Mexico, and USA
    # start by subsetting by Country, and then by excluding Hawaii from the US data
cin <- subset(data, Country == "Canada" | Country == "Mexico" | Country == "USA")
    # now subset by excluding Hawaii
cin1 <- subset(cin, Locality1 != "Hawaii")


## 3. wrangle the matrix
  # extract species names
colnames(data)
sp.names <- data[8]  # change number depending labels: e.g. for Gal, SAm, Hisp is 30
  # log transform the data
sp.ln <- log(data[, c(11:26)])
  # combine names + ln data
sp <- cbind(sp.names, sp.ln)

  # check variable names & fix the Ssp to sp
colnames(sp)
colnames(sp)[1] <- paste("sp") # change number depending labels


## 4. create train/test splits
set.seed(123)
training <- sp$sp %>% createDataPartition(p = 0.75, list = F)
# training set
train.sp <- sp[training, ]
# testing set
test.sp <- sp[-training, ]

# visualize the data splits
# scatterplots show if data is Gausian
# only look at the training set to prevent bias
scatterplotMatrix(train.sp[2:9])
scatterplotMatrix(train.sp[10:17])

# boxplots help you see the distribution of means
# only look at the training set to prevent bias
boxplot(train.sp[, c(2:17)], main = "Raw Data")

# preprocess to center and scaele the data (i.e. normalize)
prep.train.sp <- preProcess(train.sp[, c(2:17)], method = c("center", "scale"))
prep.train.sp.d <- predict(prep.train.sp, train.sp[, c(2:17)])
# boxplot to verify normalization 
boxplot(prep.train.sp.d, main = "Normalized Data")

## 5. fit the LDA model to training dataset
model.fit <- train(sp ~., data = train.sp, preProcess = c("center", "scale"), method = "lda")
# check results
model.fit$finalModel

# make the model predicitons to testing dataset
predictions <- predict(model.fit, newdata = test.sp)
# make a table to ensure both data & reference are factors
t <- table(factor(predictions), factor(test.sp$sp))
# calculate a confusion matrix
confusionMatrix(t)
# if both prediction & reference are factors, then use "confusionMatrix(predictions, test.sp$sp)"

# make cross validation set
kfoldcv <- trainControl(method = "cv", number = 5) 
performance.metric <- "Accuracy"


### 6. final LDA models ###
set.seed(123)
### run linear discriminant analysis LDA
hoary.lda <- lda(sp ~ ., data = sp) 
hoary.lda

# confusion matrix of LDA machine learning classification
hoary.lda.predict <- train(sp ~ .,  data = sp, method = "lda", metric = performance.metric,
                            trControl = kfoldcv, preProcess = c("center", "scale"))

# calculate a confusion matrix to see the accuracy of classification
confusionMatrix(as.factor(sp$sp), predict(hoary.lda.predict, sp))


## 7. plotting
# create stacked histograms of the LDA values for csv subspecies
hoary.lda.values <- predict(hoary.lda)
ldahist(hoary.lda.values$x[ ,1], g = sp$sp, 
        type = "histogram", col = "light gray") # [ ,1] indicates LD1, change number for the LD you want to plot


# create a scatterplot of LD to visualize discrimination
# convert LD data into a data.frame for cinereus/semotus/villosissimus
hoary_LD <- data.frame(type = sp[,1], lda = hoary.lda.values$x)


### step-by-step ggplot for LDA with all features
# set the plot theme
theme_set(theme_light())
  # colors for cinereus/semotus/villosissimus
colrs.csv <- c("#000000", "#4292c6", "#993404") # black / blue / brown

  # OR use grayscale colors
colrs.grsc <- c("#737373", "#000000", "#bdbdbd") # gray / black / darkgray


F1 <- ggplot(hoary_LD, aes(lda.LD1, lda.LD2, color = type)) + 
  geom_point(size = 3, alpha = 0.7) + stat_ellipse(level = 0.68)

  ## Optional: create non-overlapping labels for each point to ID individuals
library(ggrepel)
F.text <- ggplot(hoary_LD, aes(lda.LD1, lda.LD2, color = type)) + 
  geom_point(size = 3) + stat_ellipse(level = 0.68) +
  geom_text_repel(aes(label = data$Country))
  ## here used Locality1 to distinguish the state from where each specimen came fron
  ## this identifies that L semotus resembles L cinereous from MX and L villosissimus from Galapagos & Colombia

# continue main plot - Island, N America, S America
F1a <- F1 + scale_color_manual(values = colrs.sagh,
                               aesthetics = c("colour", "fill"),
                               labels = c("Galápagos", "Hispaniola", "South America")) +
  labs(x = "LD1 (XXX%)", y = "LD2 (XXX%)") +
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18),
        legend.position = "none")


# 8. save LD values per specimen to plot 
write.csv(hoary_LD, file = "PATH/TO/SAVE/*.csv", row.names = F)

### get the xy values for the Canada+MX+USA comparison of LDAs ###
colnames(data)
xy <- data[, 28:29]
# now cbind
hoary.LD.xy <- cbind(hoary_LD, xy)
head(hoary.LD.xy)

## run a linear regression of LD1 against Latitude to examine size variation
cin.lat.lm <- lm(lda.LD1 ~ lat, data = hoary.LD.xy)
#example <- lm(PC1 ~ bio4 + bio16 + bio17 + lon + lat, data = ddata)
summary(cin.lat.lm)

# LD1 vs latitude
Fig.lat <- ggplot(hoary.LD.xy, aes(x = lat, y = lda.LD1)) +
  stat_smooth(method = "lm", formula = y ~ x, se = F) +
  geom_point(aes(colour = type), size = 5, position = position_jitter(width = 0.2), alpha = 0.5) +
  labs(x = "Latitude", y = "Phenotypic Variation (LD1)") +
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18),
        legend.position = "none")
# add the custom color palete
Fig.lat + scale_color_manual(values = colrs.sagh)

### use mean values for the most important cranial characters to evaluate size ###
# most important characters for LD1 = GLS_24, ROST.LEN_26, MASTOID_34, LOW.TR_49
  # NOTE: measurements with the highest coefficients for LD1 are all related to length
# a. subset each country
can <- subset(data, Country == "Canada")
mex <- subset(data, Country == "Mexico")
usa <- subset(data, Country == "USA")
# b. get mean GSL
mean(can$GSL_24)
mean(mex$GSL_24)
mean(usa$GSL_24)
# c. get mean ROST.LEN_26
mean(can$ROST.LEN_26)
mean(mex$ROST.LEN_26)
mean(usa$ROST.LEN_26)
# d. get mean MASTOID_34
mean(can$MASTOID_34)
mean(mex$MASTOID_34)
mean(usa$MASTOID_34)
# e. get mean LOW.TR_49
mean(can$LOW.TR_49)
mean(mex$LOW.TR_49)
mean(usa$LOW.TR_49)

## create a box plot panel for L cinereus characters: GLS_24, ROST.LEN_26, MASTOID_34, LOW.TR_49
# source: http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/77-facilitating-exploratory-data-visualization-application-to-tcga-genomic-data/
# first extract the charaacter columns of the main dataset
colnames(data)
lcin <- data[, c(8,11,12,15,21)]
# change the order of Country so that it shows not alphabetically
level_ord <- c('Canada', 'USA', 'Mexico')

# create character boxplots
box.GSL <- lcin %>%
  ggplot(aes(x = Country, y = GSL_24)) +
  geom_boxplot(aes(x = factor(Country, level = level_ord)), alpha = 0.5, fill = c("#01665e", "#993404", "#053061")) +
  geom_point(aes(colour = Country), size = 3, position = position_jitter(width = 0.1), alpha = 0.6) +
  labs(x = "Country", y = "Greatest skull length") +
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        legend.position = "none") # change legend.position to: "right" then add < , legend.text = element_text(size = 12) >  so it can be displayed
# add custom color palette
gsl.f <- box.GSL + scale_color_manual(values = c("#01665e", "#053061", "#993404"))

box.RL <- lcin %>%
  ggplot(aes(x = Country, y = ROST.LEN_26)) +
  geom_boxplot(aes(x = factor(Country, level = level_ord)), alpha = 0.5, fill = c("#01665e", "#993404", "#053061")) +
  geom_point(aes(colour = Country), size = 3, position = position_jitter(width = 0.1), alpha = 0.6) +
  labs(x = "Country", y = "Rostral length") +
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), 
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        legend.position = "none") # change legend.position to: "right" then add < , legend.text = element_text(size = 12) >  so it can be displayed
# add custom color palette
rl.f <- box.RL + scale_color_manual(values = c("#01665e", "#053061", "#993404"))

box.MAS <- lcin %>%
  ggplot(aes(x = Country, y = MASTOID_34)) +
  geom_boxplot(aes(x = factor(Country, level = level_ord)), alpha = 0.5, fill = c("#01665e", "#993404", "#053061")) +
  geom_point(aes(colour = Country), size = 3, position = position_jitter(width = 0.2), alpha = 0.6) +
  labs(x = "Country", y = "Breadth at mastoid") +
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), 
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        legend.position = "none") # change legend.position to: "right" then add < , legend.text = element_text(size = 12) >  so it can be displayed
# add custom color palette
mas.f <- box.MAS + scale_color_manual(values = c("#01665e", "#053061", "#993404"))

box.TR <- lcin %>%
  ggplot(aes(x = Country, y = LOW.TR_49)) +
  geom_boxplot(aes(x = factor(Country, level = level_ord)), alpha = 0.5, fill = c("#01665e", "#993404", "#053061")) +
  geom_point(aes(colour = Country), size = 3, position = position_jitter(width = 0.2), alpha = 0.6) +
  labs(x = "Country", y = "Lower toothrow") +
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), 
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        legend.position = "none") # change legend.position to: "right" then add < , legend.text = element_text(size = 12) >  so it can be displayed
# add custom color palette
tr.f <- box.TR + scale_color_manual(values = c("#01665e", "#053061", "#993404"))

# arrange the four boxplots
grid.arrange(gsl.f, rl.f, mas.f, tr.f,
             ncol = 2, nrow = 2, heights = c(4, 4))

### --- clear up your Global Environment --- ###
rm(list = ls())
### --- clear up plots --- ###
dev.off()



## for the villosissimus species comparison ##
data <- vil # replace the data.frame
colnames(data)
fix(data) # fix a column with the subspecies names
sp.names <- data[30]  # change number: this one uses "Country" to define subspecies of villosissimus
# log transform the data
sp.ln <- log(data[, c(11:26)])
# combine names + ln data
sp <- cbind(sp.names, sp.ln)

# check variable names & fix the Ssp to sp
colnames(sp)
colnames(sp)[1] <- paste("sp")  # change name: will be different for each dataset depending on the labels you want


# colors for L villosissimus subspecies
# Galapagos / grayi / pasescens / villosissimus
colrs.csv <- c("#FF9400", "#F9E231", "#047101", "#00A1FE", "#1500FF")

library(ggbiplot) # for pca plotting

# run PCA, examine PC contributions & rotation
vil_pca <- prcomp(sp.ln, center = T, scale. = T)
summary(vil_pca)
vil_pca$rotation

#create scree plot of pca proportion of variance
ggscreeplot(vil_pca)

# extract the PCs into a data.frame for plotting
# vil_pca$x <- as.data.frame(vil_pca)

# create pca biplot
# set the theme of the plot
theme_set(theme_bw())
# set raw plot
vil_45 <- ggbiplot(vil_pca, choices = c(4,5), obs.scale = 1, var.scale = 1, groups = data$vssp, ellipse = T, 
                    circle = F, var.axes = F, size = 2.5, alpha = 0.75)
## change var.axes = T to see the character vectors
## ellipses represent 68% Gausian data ellipses ~ standard deviation & size shows the variance
## used choices = c(x,y) to test different combinations of PCs in multiple plots
# set final plot
vil_45 + scale_color_manual(values = colrs.csv,
                             labels = c("Galápagos", "L.v. grayi", "L.v. pallescens", "L.v. villosissimus")) +
  labs(x = "PC4 (8.4%)", y = "PC5 (7.2%)") +
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), legend.position = "none") +
  xlim(-6,6) + ylim(-4,4)
       # legend.title = element_blank(), legend.text = element_text(face = "italic", size = 15))
          # replace the line of text above for legend.position = "none" to get the text legend
  # NOTE: used fixed x and y axis limits to make the plots size comparable for multi-plot figure


