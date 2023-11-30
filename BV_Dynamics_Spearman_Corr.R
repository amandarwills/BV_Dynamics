################################################################################
#Title: Spearman correlation of immune markers and bacterial abundances
#Project: BV Dynamics
#Author: Amanda Williams 
#Date Last Modified: 20231127
#See Readme file for details
################################################################################

# Clear workspace and restart R
rm(list=ls())  
.rs.restartR()

# Load libraries
library(devtools) 
library(magrittr)
require(readxl)
library(tidyverse)
library(dplyr)
library(matrixStats)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

# Set seed
set.seed(54321)


# Set working directory
workingDir <- "/Users/amandawilliams/Desktop/BVDynamics/DATA/"
setwd(workingDir)

# Load files
hmp.initial <- as.data.frame(read_excel(path = paste(workingDir,"AW_HMP_metadata_updated_taxonomy_StR_CSTs_151123.xlsx", sep = ""), sheet = 1, col_names = TRUE, na = "NA"))
immune.initial <- read.table("BV_dynamics_Complete_10Nov2023.csv", header = T, row.names = NULL, check.names = FALSE, sep = ',')

### Format immune data ###

# Remove unecessary columns
colnames(immune.initial)
drop <- c("visit", "visit.bv", "visit.metro",	"visit.azi", "visit.bv.details", 
          "visit.metro.details", "visit.azi.details", "SID", "Metagenomic_total", 
          "Metabolomic_all", "sampleID_16S", "barcode", "available", "aliquoted", 
          "aliquot.status", "diarySymbv", "MG_2023", "CK_2023")  
immune <- immune.initial[,!(names(immune.initial) %in% drop)] 

# Reorder columns for easier analysis
colnames(immune)
immune <- immune[ , c(15, 1, 16, 2, 19, 17, 18, 22:24, 3:14, 20, 21, 25)] 

# Remove first instance of duplicated rows (bv.dx)
immune <- immune %>% 
                  group_by(UID) %>% 
                                  filter(duplicated(UID) | n() == 1) 
immune <- data.frame(immune)

# Filter for only MET study
immune <- immune %>% 
                    group_by(study) %>% 
                                      filter(study == "Met_Study") 
immune$STATUS[immune$STATUS == 'PRIOR.TO.SBV'] <- 'PRIOR.TO.MET'
immune$STATUS[immune$STATUS == 'BV.DX'] <- 'PRIOR.TO.MET'


### Pull out rows from hmp.initial with patients of interest ###
patients <- c("UAB003","UAB005","UAB017","UAB053", "UAB127","UAB128","UAB130","UAB135")
hmp <- hmp.initial[hmp.initial$SID %in% patients, ]

# Join immune to hmp, keeping all of immune
data <- immune %>%
                  left_join(hmp, by = "UID", suffix = c("",".y")) %>%
                                                                    select(-ends_with(".y"))

### Perform Spearman Corr (taxa-cytokine) for each patient ###
for (i in patients) {
  data2 <- data %>% 
                  filter(ID == i) 
  
  ## Pull out immune data and format for analysis
  immune.4.analysis <- data2[ , c(1, 11:22)]
  immune.4.analysis <- data.frame(immune.4.analysis)
  rownames(immune.4.analysis) <- immune.4.analysis$UID
  immune.4.analysis[is.na(immune.4.analysis)] <- 0
  immune.4.analysis <- immune.4.analysis[ , -1]
  
  # Create list of filtered variables (cytokines)
  t.immune.4.analysis <- t(immune.4.analysis)
  t.immune.4.analysis <- data.frame(t.immune.4.analysis)
  cytokines <- rownames_to_column(t.immune.4.analysis, var = "Cytokines")
  cytokines <- data.frame(cytokines$Cytokines)%>%
                                              rename(Cytokines = 1)
  
  
  ## Pull out amplicon data and format for analysis
  amp.4.analysis <- data2[ , c(1, 58:256)]
  amp.4.analysis <- data.frame(amp.4.analysis)
  rownames(amp.4.analysis) <- amp.4.analysis$UID
  amp.4.analysis <- amp.4.analysis[ , -1]
  
  # Filter amplicon data for species with >X reads
  amp.4.analysis[is.na(amp.4.analysis)] <- 0
  mat_fun <- function(m){
                          m2 <- apply(m,  2,  function(x) as.numeric(paste(x)))
                          colnames(m2) <- colnames(m)
                          rownames(m2) <- rownames(m)
                          m2 <- m2[ , colSums(m2) > 200]
                          m2 <- data.frame(m2)
                          return(m2)
  }
  
  amp.4.analysis <- mat_fun(amp.4.analysis)
  
  # Create list of filtered variables (species)
  t.amp.4.analysis <- t(amp.4.analysis)
  t.amp.4.analysis <- data.frame(t.amp.4.analysis)
  species <- rownames_to_column(t.amp.4.analysis, var = "Species")
  species <- data.frame(species$Species) %>%
                                            rename(Species = 1)
  
  ## Join data frames for analysis
  joined.data <- amp.4.analysis %>% merge(immune.4.analysis, by = 0)
  rownames(joined.data) <- joined.data$Row.names
  joined.data$Row.names <- NULL
  
  ## Create data frame of every combination of metabolite and class with leading data frame name
  combinations <- tidyr::crossing(species, cytokines)
  combinations.merge <- combinations %>% rename(Species = 1,
                                                Cytokines = 2)
  combinations <- t(combinations)
  combinations <- data.frame(combinations)
  names(combinations) = str_sub(names(combinations), 2)
  
  ## Run the Spearman correlation test
  spearman.output <- lapply(combinations, function(x) {
                                        (cor.test(joined.data[ ,x[1]], joined.data[ ,x[2]], 
                                                          method = "spearman",
                                                          exact = FALSE))
  })
  
  ## Format spearman output list into data frame
  spearman.output.df <- data.frame(matrix(unlist(spearman.output), nrow = length(spearman.output), byrow = TRUE))
  
  spearman.output.df <- spearman.output.df[, c(1:5)]
  spearman.output.df <- spearman.output.df %>%
                                        rename(S.statistic = 1,
                                               pavalue = 2,
                                               rho = 3,
                                               null.value = 4,
                                               type = 5)
  
  spearman.output.df <- spearman.output.df %>% merge(combinations.merge, by = 0)
  spearman.output.df <- spearman.output.df[ , -1]
  spearman.output.df$rho <- as.numeric(spearman.output.df$rho)
  
  ## Analysis of species of interest 
  
  # Filter results for high significance correlations
  output.analysis <- spearman.output.df %>%
                                        arrange(desc(rho))
  
  
  gard.data <- output.analysis %>% 
                                  dplyr::filter(Species == "Gardnerella_vaginalis")
  lacto.data <- output.analysis %>% 
                              filter(if_any(everything(), ~ grepl('Lactobacillus',.)), .by_group = TRUE) %>%
                                                                                                          arrange(Cytokines)
  
  lacto.gard.data <- lacto.data %>% full_join (gard.data)
  
  
  ## Unmelt ambient table and format for plotting
  corr.results <- spearman.output.df[ , c(6, 7, 3)] 
  corr.results <- corr.results[order(corr.results$Species, corr.results$Cytokines),]
  
  
  corr.results <- dcast(corr.results, Species ~ Cytokines)
  corr.results <- corr.results %>% 
                                  arrange(Species)
  rownames(corr.results) <- corr.results$Species
  corr.results <- corr.results[ , -1]
  
  # Transform df into matrix
  corr.results.mat <- as.matrix(corr.results)
  mode(corr.results.mat)
  
  
  ## Plot heatmap 
  col.scheme <- colorRamp2(seq(min(corr.results.mat), max(corr.results.mat), length = 3), c("blue", "#EEEEEE", "red"), space = "RGB")
  
  all.heatmap.plot <- Heatmap(corr.results.mat, 
                              col = col.scheme,
                              heatmap_legend_param = list(title = "Correlation coefficient", title_gp = gpar(fontsize = 8), 
                                                          labels_gp = gpar(fontsize = 8)),
                              column_title_gp = gpar(fontsize = 10),
                              row_names_max_width = unit(4, "cm"),
                              row_names_gp = gpar(fontsize = 6),
                              column_names_max_height = unit(4, "cm"),
                              column_names_gp = gpar(fontsize = 6),
                              rect_gp = gpar(col = "white", lwd = 1),
                              show_row_dend = TRUE,
                              row_dend_side = "right",
                              show_column_dend = TRUE,
                              column_dend_side = "bottom",
                              width = unit(8, "cm"), height = unit(20, "cm"),
                              column_title = (paste("Spearman correlation of amplicon and cytokine data - ", i, sep = "")))
  
  all.heatmap.plot
  
  ## Plot species of interest (Lacto, Gard)
  
  # Subset species of interest for analysis from data frame
  t.corr.results <- t(corr.results) %>%
                                      data.frame()
  lacto.gard <- t.corr.results %>% 
                                  dplyr:: select(starts_with("Lactobacillus") | starts_with("Gardnerella"))
  lacto.gard.mat <- t(lacto.gard)

  ## Plot heatmap 
  col.scheme2 <- colorRamp2(seq(min(lacto.gard.mat), max(lacto.gard.mat), length = 3), c("blue", "#EEEEEE", "red"), space = "RGB")
  
  lac.gar.heatmap.plot <- Heatmap(lacto.gard.mat, 
                              col = col.scheme2,
                              heatmap_legend_param = list(title = "Correlation coefficient", title_gp = gpar(fontsize = 8), 
                                                          labels_gp = gpar(fontsize = 8)),
                              column_title_gp = gpar(fontsize = 10),
                              row_names_max_width = unit(4, "cm"),
                              row_names_gp = gpar(fontsize = 6),
                              column_names_max_height = unit(6, "cm"),
                              column_names_gp = gpar(fontsize = 6),
                              rect_gp = gpar(col = "white", lwd = 1),
                              show_row_dend = TRUE,
                              row_dend_side = "right",
                              show_column_dend = TRUE,
                              column_dend_side = "bottom",
                              width = unit(6, "cm"), height = unit(8, "cm"),
                              column_title = (paste("Spearman correlation of amplicon and cytokine data - ", i, sep = "")))
  
  lac.gar.heatmap.plot
  
  ## Save plots
  pdf(file = paste(i, "_Spearman_Plots.pdf", sep = ""), paper="letter", width = 8, height = 25, onefile = TRUE)
    par(mai = c(1, 0.5, 0.5, 0.25), pin = c(3, 3))
    print(all.heatmap.plot)
    print(lac.gar.heatmap.plot)
  dev.off()
  
}
