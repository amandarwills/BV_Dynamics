################################################################################
#Title: Multivariate analysis of immune markers and bacterial abundance data
#Project: BV Dynamics
#Author: Amanda Williams 
#Date Last Modified: 20231115
#See Readme file for details
################################################################################

# Clear workspace and restart R
rm(list=ls())  
.rs.restartR()

# Install new vegan
library(remotes)
remotes::install_github("vegandevs/vegan")

# Load libraries
library(devtools) 
library(magrittr)
library(dplyr)
library(vegan)
library(mixOmics)
library(ggplot2)
library(ggbiplot)
library(ggforce)
library(scales)
library(ggfortify)
library(ggpubr)
library(gridExtra)

# Set seed
set.seed(54321)


# Set working directory
workingDir <- "/Users/amandawilliams/Desktop/BVDynamics/DATA"
setwd(workingDir)

# Set data files names
sample.info.file <- "BV_dynamics_Complete_10Nov2023.csv"
sample.class.file <- "Sample_class_within_person.csv"
inf.sample.class.file <- "Sample_class_infection.csv"

# Set plot titles
MET.pca.plot.title <- "PCA of patients with BV and treated with MET immune data"
AZI.pca.plot.title <- "PCA of patients with BV and treated with AZI immune data"
MET.plsda.plot.title <- "sPLS-DA of patients with BV and treated with MET immune data"

# Define plot settings
options(repr.plot.width = 2, repr.plot.height = 3)
label.decimal.rounding <- 1
plot.text.size  <- 12
plot.title.size <- 10

status.shapes<-c(16,17,15,3,7, 8, 9, 10)

STATUS.Colors <- c("AFTER.AZI" = "#00ebf7", "DURING.AZI" = "#01ff00", 
                  "PRIOR.TO.AZI" = "#006600", "PRIOR.TO.SBV" = "#090c47", 
                  "BV.DX" = "#660099", "AFTER.MET" = "#fad700", "DURING.MET" = "#ff7c00", 
                  "PRIOR.TO.MET" = "#ff0000")

Immune.Colors <- c("sEcad" = "#8DD3C7", "IFNa2a" = "#FFFFB3", "IL17A" = "#BEBADA", 
                   "IL1a" = "#FB8072", "IL1b" = "#80B1D3", "IL6" = "#FDB462", 
                   "IL8" = "#B3DE69", "IP10" = "#FCCDE5", "MIG" = "#D9D9D9", 
                   "MIP1b" = "#BC80BD", "MIP3a" = "#CCEBC5", "MMP9" = "#FFED6F")


# Load files
data <- read.table(sample.info.file, header=T, row.names = NULL, check.names=FALSE, sep=',')
patient.sample.class <- read.table(sample.class.file, header=F, sep=',')
infection.sample.class <- read.table(inf.sample.class.file, header=F, sep=',')


### Format data ###

#Remove unecessary columns
colnames(data)
drop <- c("visit", "visit.bv", "visit.metro",	"visit.azi", "visit.bv.details", 
          "visit.metro.details", "visit.azi.details", "SID", "Metagenomic_total", 
          "Metabolomic_all", "sampleID_16S", "barcode", "available", "aliquoted", 
          "aliquot.status", "diarySymbv", "MG_2023", "CK_2023")  
master <- data[,!(names(data) %in% drop)] 

# Reorder columns for easier analysis
colnames(master)
master <- master[, c(15, 1, 16, 2, 19, 17, 18, 22:24, 3:14, 20, 21, 25)] 

# Remove first instance of duplicated rows (bv.dx)
master <- master %>% 
              group_by(UID) %>% 
                              filter(duplicated(UID) | n()==1) 
master <- data.frame(master)

################################## MET STUDY ###################################

### log transform immune marker abundances ###

# Filter for only MET study
MET.master <- master %>% 
                        group_by(study) %>% 
                                          filter(study == "Met_Study") 
MET.master <- data.frame(MET.master)
rownames(MET.master) <- MET.master[[1]]
MET.master[[1]] <- NULL

# Pull immune marker data to new df
colnames(MET.master)
MET.immune.data <- MET.master[, c(10:21)] 

# Remove NAs until NAs fixed
MET.immune.data[is.na(MET.immune.data)] <- 0
sum(is.na(MET.immune.data))
which(is.na(MET.immune.data))

# Perform log10 transformation
log.MET.immune.data <- log10(MET.immune.data+1)

# Rejoin log transformed immune data to complete df
#data.no.immune <- data[, c(14, 15, 1, 18, 16, 17, 21:23, 19, 20, 24)] 
#analysis.data <- data.no.immune %>%
                                  #left_join(data.log.cyto, by = "UID") 


### PCA ###

# Rename df
#t.log.MET.immune.data <- t(log.MET.immune.data)
#MET.pca.in <- data.frame(t.log.MET.immune.data)
MET.pca.in <- data.frame(log.MET.immune.data)

# Create PCA object (base R)
MET.pca.out <- prcomp(MET.pca.in,
                      center = FALSE,
                      scale. = FALSE) # Use scale = TRUE if your variables are on different scales (e.g. for abiotic variables)

summary(MET.pca.out)
MET.pca.out.summary <- summary(MET.pca.out)


# Create plotting df by pulling out axes and joining to metadata
MET.PC1.df <- data.frame(MET.pca.out$x[ , 1])
names(MET.PC1.df)[1] <- 'PC1'

MET.PC2.df <- data.frame(MET.pca.out$x[,2])
names(MET.PC2.df)[1] <- 'PC2'

# Remove immune data from master
MET.metadata <- MET.master[, c(1:9, 22:24)]

# Create function to merge multiple df together by rownames
MyMerge <- function(x, y){
  df            <- merge(x, y, by = 0, all.x= F, all.y= F)
  rownames(df)  <- df$Row.names
  df$Row.names  <- NULL
  return(df)}

plot.pca <- Reduce(MyMerge, list(MET.metadata, MET.PC1.df, MET.PC2.df))

# Calculate the percent of variance explained by first two axes and format for plot labels
MET.PCA.axis1.var <- MET.pca.out.summary$importance[2,1]
MET.PCA.axis2.var <- MET.pca.out.summary$importance[2,2]
MET.PCA.axis1.var.label <- paste("PC1 (", round(MET.PCA.axis1.var * 100, label.decimal.rounding), "%)", sep = '')
MET.PCA.axis2.var.label <- paste("PC2 (", round(MET.PCA.axis2.var * 100, label.decimal.rounding), "%)", sep = '')

# PCA plot
MET.pca.plot <- ggplot(plot.pca, aes(x = PC1, y = PC2)) +
                        geom_point(aes(color = STATUS, shape = ID)) +
                        scale_color_manual(values = STATUS.Colors) +
                        scale_shape_manual(values = status.shapes) +
                        geom_mark_ellipse(aes(color = STATUS),
                                expand = unit(0.5, "mm")) +
                        labs(x = MET.PCA.axis1.var.label,
                             y = MET.PCA.axis2.var.label) +
                        ggtitle(MET.pca.plot.title) +
                        theme(text = element_text(family = "Helvetica", size = plot.text.size),
                        plot.title = element_text(family = "Helvetica", size = plot.title.size, hjust = 0.5)) +
                        theme(aspect.ratio = 1)

MET.pca.plot


### PERMANOVA ###

# Use square root or proportions to minimize influence of most abundant groups
MET.data.mat <- sqrt(MET.pca.in)
sum(is.na(MET.data.mat))
which(is.na(MET.data.mat), arr.ind = TRUE) #where NA's are located, if present

# Create a dissimilarity matrix (vegan)
MET.data.dist <- vegdist(MET.data.mat, method = 'bray')

# Format master df to remove immune data & NAs
MET.perm.master <- MET.master[ , -c(10:21)] 
MET.perm.master[is.na(MET.perm.master)] <- 0
MET.perm.master$STATUS[MET.perm.master$STATUS == 'PRIOR.TO.SBV'] <- 'PRIOR.TO.MET'
MET.perm.master$STATUS[MET.perm.master$STATUS == 'BV.DX'] <- 'PRIOR.TO.MET'

# Run perMANOVA (vegan)

# Set seed
set.seed(69)

(h <- how(within = Within(type = "series"), plots = Plots(strata = MET.perm.master$ID, type = "none")))

MET.mod.data <- adonis2(MET.data.dist ~ CST * ID * STATUS, data = MET.perm.master, permutations = how(nperm = 5040, minperm = 9999), method = "bray", strata = MET.perm.master$ID)
MET.mod.data

MET.mod.data2 <- adonis2(MET.data.dist ~ ID * CST * STATUS, data = MET.perm.master, permutations = how(nperm = 5040, minperm = 9999), method = "bray", strata = MET.perm.master$ID)
MET.mod.data2


with(master, adonis2(MET.data.dist ~ STATUS * CST, data = MET.perm.master, permutations = how(nperm = 5040, minperm = 9999), strata = MET.perm.master$ID, method = 'bray'))



################################## AZI STUDY ###################################

### log transform immune marker abundances ###

# Filter for only AZI study
AZI.master <- master %>% 
                      group_by(study) %>% 
                                        filter(study == "Azi_Study") 
AZI.master <- data.frame(AZI.master)
rownames(AZI.master) <- AZI.master[[1]]
AZI.master[[1]] <- NULL

# Pull immune marker data to new df
colnames(AZI.master)
AZI.immune.data <- AZI.master[, c(10:21)] 

# Remove NAs until NAs fixed
AZI.immune.data[is.na(AZI.immune.data)] <- 0
sum(is.na(AZI.immune.data))
which(is.na(AZI.immune.data))

# Perform log10 transformation
AZI.log.immune.data <- log10(AZI.immune.data+1)


### PCA ###

# Rename df
#t.AZI.log.immune.data <- t(AZI.log.immune.data)
#AZI.pca.in <- data.frame(t.AZI.log.immune.data)
AZI.pca.in <- data.frame(AZI.log.immune.data)

# Create PCA object (base R)
AZI.pca.out <- prcomp(AZI.pca.in,
                      center = FALSE,
                      scale. = FALSE) # Use scale = TRUE if your variables are on different scales (e.g. for abiotic variables)

summary(AZI.pca.out)
AZI.pca.out.summary <- summary(AZI.pca.out)


# Create plotting df by pulling out axes and joining to metadata
AZI.PC1.df <- data.frame(AZI.pca.out$x[,1])
names(AZI.PC1.df)[1] <- 'PC1'

AZI.PC2.df <- data.frame(AZI.pca.out$x[,2])
names(AZI.PC2.df)[1] <- 'PC2'

# Remove immune data from master
AZI.metadata <- AZI.master[, c(1:9, 22:24)]

# Create function to merge multiple df together by rownames
MyMerge <- function(x, y){
  df            <- merge(x, y, by = 0, all.x= F, all.y= F)
  rownames(df)  <- df$Row.names
  df$Row.names  <- NULL
  return(df)}

AZI.plot.pca <- Reduce(MyMerge, list(AZI.metadata, AZI.PC1.df, AZI.PC2.df))

# Calculate the percent of variance explained by first two axes and format for plot labels
AZI.PCA.axis1.var <- AZI.pca.out.summary$importance[2,1]
AZI.PCA.axis2.var <- AZI.pca.out.summary$importance[2,2]
AZI.PCA.axis1.var.label <- paste("PC1 (", round(AZI.PCA.axis1.var * 100, label.decimal.rounding), "%)", sep = '')
AZI.PCA.axis2.var.label <- paste("PC2 (", round(AZI.PCA.axis2.var * 100, label.decimal.rounding), "%)", sep = '')

# PCA plot
AZI.pca.plot <- ggplot(AZI.plot.pca, aes(x = PC1, y = PC2)) +
                  geom_point(aes(color = STATUS, shape = ID)) +
                  scale_color_manual(values = STATUS.Colors) +
                  scale_shape_manual(values = status.shapes) +
                  geom_mark_ellipse(aes(color = STATUS),
                                    expand = unit(0.5, "mm")) +
                  labs(x = AZI.PCA.axis1.var.label,
                       y = AZI.PCA.axis2.var.label) +
                  ggtitle(AZI.pca.plot.title) +
                  theme(text = element_text(family = "Helvetica", size = plot.text.size),
                        plot.title = element_text(family = "Helvetica", size = plot.title.size, hjust = 0.5)) +
                  theme(aspect.ratio = 1)

AZI.pca.plot


### PERMANOVA ###

# Use square root or proportions to minimize influence of most abundant groups
AZI.data.mat <- sqrt(AZI.pca.in)
sum(is.na(AZI.data.mat))
which(is.na(AZI.data.mat), arr.ind = TRUE) #where NA's are located, if present

# Create a dissimilarity matrix (vegan)
AZI.data.dist <- vegdist(AZI.data.mat, method = 'bray')

# Format master df to remove immune data & NAs
AZI.perm.master <- AZI.master[ , -c(10:21)] 
AZI.perm.master[is.na(AZI.perm.master)] <- 0

# Run perMANOVA (vegan)

# Set seed
set.seed(69)

(h <- how(within = Within(type = "series"), plots = Plots(strata = AZI.perm.master$ID, type = "none")))

AZI.mod.data <- adonis2(AZI.data.dist ~ STATUS * ID, data = AZI.perm.master, permutations = how(nperm = 5040, minperm = 9999), method = "bray", strata = AZI.perm.master$ID)
AZI.mod.data

with(master, adonis2(AZI.data.dist ~ STATUS * CST, data = AZI.perm.master, permutations = how(nperm = 5040, minperm = 9999), strata = AZI.perm.master$ID, method = 'bray'))

