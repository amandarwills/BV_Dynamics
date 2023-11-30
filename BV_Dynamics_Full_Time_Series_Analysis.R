################################################################################
#Title: Analysis of full time series for MET treated patients
#Project: BV Dynamics
#Author: Amanda Williams 
#Date Last Modified: 20231122
#See Readme file for details
################################################################################

# Clear workspace and restart R
rm(list=ls())  
.rs.restartR()

# Load libraries
library(devtools) 
library(magrittr)
require(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpmisc)
library(vcd)
library(ggpubr)

# Set seed
set.seed(54321)


# Set working directory
workingDir <- "/Users/amandawilliams/Desktop/BVDynamics/DATA/"
setwd(workingDir)

# Load files
hmp.initial <- as.data.frame(read_excel(path = paste(workingDir,"AW_HMP_metadata_updated_taxonomy_StR_CSTs_151123.xlsx", sep=""), sheet = 1, col_names = TRUE, na = "NA"))
immune.initial <- read.table("BV_dynamics_Complete_10Nov2023.csv", header=T, row.names = NULL, check.names=FALSE, sep=',')

### Format immune data ###

#Remove unecessary columns
colnames(immune.initial)
drop <- c("visit", "visit.bv", "visit.metro",	"visit.azi", "visit.bv.details", 
          "visit.metro.details", "visit.azi.details", "SID", "Metagenomic_total", 
          "Metabolomic_all", "sampleID_16S", "barcode", "available", "aliquoted", 
          "aliquot.status", "diarySymbv", "MG_2023", "CK_2023")  
immune <- immune.initial[,!(names(immune.initial) %in% drop)] 

# Reorder columns for easier analysis
colnames(immune)
immune <- immune[, c(15, 1, 16, 2, 19, 17, 18, 22:24, 3:14, 20, 21, 25)] 

# Remove first instance of duplicated rows (bv.dx)
immune <- immune %>% 
                    group_by(UID) %>% 
                                      filter(duplicated(UID) | n()==1) 
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

# Join immune to hmp, keeping all of hmp
data <- hmp %>%
              left_join(immune, by = "UID", suffix = c("",".y")) %>%
                                                                    select(-ends_with(".y"))

# Remove NAs in certain columns
data$MENSTRUATION[is.na(data$MENSTRUATION)] <- 0

# Code Nugent Class and CST columns for plotting 
data <- data.frame(append(data, list(Nugent_Code = data$NUGENT_CLASS), after = 8))
data <- data.frame(append(data, list(CST_Code = data$CST), after = 10))

data$Nugent_Code[data$Nugent_Code == 'NO_BV'] <- 1
data$Nugent_Code[data$Nugent_Code == 'INTER_BV'] <- 1.5
data$Nugent_Code[data$Nugent_Code == 'BV'] <- 2

unique(data$CST)
data$CST_Code[data$CST_Code == 'I'] <- 1
data$CST_Code[data$CST_Code == 'III'] <- 1.25
data$CST_Code[data$CST_Code == 'IV-A'] <- 1.5
data$CST_Code[data$CST_Code == 'IV-B'] <- 1.75
data$CST_Code[data$CST_Code == 'IV-C'] <- 2

data <- data %>% 
              mutate_at(c(9, 11), as.numeric)
data$WEEK <- as.character(data$WEEK)


### Fill in STATUS column based on dates ### 
data$ID <- seq.int(nrow(data))

data <- data %>% 
                relocate(STATUS, .after = DAY) %>%
                                                  relocate(ID)

data <- data %>% 
                  mutate(TREATMENT = case_when(ID <= 8  ~ "Prior to MET",
                                               ID >= 9 & ID <= 15  ~ "During MET",
                                               ID >= 52 & ID <= 59  ~ "Prior to MET",
                                               ID >= 60 & ID <= 66  ~ "During MET",
                                               ID >= 71 & ID <= 77  ~ "Prior to MET",
                                               ID >= 78 & ID <= 84  ~ "During MET",
                                               ID >= 141 & ID <= 186  ~ "Prior to MET",
                                               ID >= 187 & ID <= 193  ~ "During MET",
                                               ID >= 213 & ID <= 247  ~ "Prior to MET",
                                               ID >= 248 & ID <= 254  ~ "During MET",
                                               ID >= 283 & ID <= 310  ~ "Prior to MET",
                                               ID >= 311 & ID <= 317  ~ "During MET",
                                               ID >= 353 & ID <= 360  ~ "Prior to MET",
                                               ID >= 361 & ID <= 367  ~ "During MET",
                                               ID >= 383 & ID <= 392  ~ "Prior to MET",
                                               ID >= 393 & ID <= 399  ~ "During MET",
                                               ID >= 423 & ID <= 459  ~ "Prior to MET",
                                               ID >= 460 & ID <= 466  ~ "During MET",
                                               ID >= 493 & ID <= 532  ~ "Prior to MET",
                                               ID >= 533 & ID <= 539  ~ "During MET",
                                               TRUE ~ "After MET")) %>%
                                                                        relocate(TREATMENT, .after = STATUS)       
                  

data <- data[!is.na(data$NUGENT_CLASS), ]

### Plots of Nugent Class and CST correlation ### 

# Define plot settings
theme_set(theme_bw())

Colors <- c("Prior to MET" = "#238B45", "During MET" = "#9ECAE1", "After MET" = "#00441B", 
            "NO_BV" = "#CCECE6", "INTER_BV" = "#66C2A4", "BV" = "#00441B", "I" = "#FE0308", 
            "I-A" = "#FE0308", "I-B" = "#F6D3DA", "II" = "#86C61A", "III" = "#FF7200", 
            "III-A" = "#FF7200", "III-B" = "#F8A40E", "IV" = "#221886", "IV-A" = "#448A73", 
            "IV-B" = "#221886", "IV-C" = "#C0ACD3", "IV-C0" = "#989898", "IV-C1" = "#EF53A7", 
            "IV-C2" = "#A7DDDC", "IV-C3" = "#98C999", "IV-C4" = "#7F0B7C", "V" = "#FAE50D", 
            "NA" = "#FFFFFF")

## Produce plots for each patient (CST and Nugent day series, Jitter, Association) 
for (i in patients) {
  data2 <- data %>% 
                  filter(SID == i) 

  # Create CST and Nugent day series scatter plots 
  lower <- (min(data2$SERIAL, na.rm = TRUE)) + 1
  upper <- (max(data2$SERIAL, na.rm = TRUE)) - 1
  
  met.rows <- which(grepl("During MET", data2$TREATMENT))
  x <- str_count(data2$TREATMENT, "During MET")  
  x.sum <- sum(x)
  
  if (x.sum <= 7) {
    met.date1 <- (min(met.rows, na.rm = TRUE)) - 0.5
    met.date2 <- (max(met.rows, na.rm = TRUE)) + 0.5
    vertical.lines <- c(met.date1, met.date2) %>%
      sort()
    text1 <- met.date1 + 3.5
    
    ## Plot ##
    n <- ggplot(data2, aes(x = SERIAL, y = Nugent_Code)) +
                geom_point(size = 3, aes(colour = CST, shape = TREATMENT)) +
                scale_color_manual(values = Colors) +
                guides(colour = guide_legend(title = "CST Group")) +
                scale_x_continuous(name = "Day") +
                xlim(lower, upper) +
                scale_y_continuous(name = "Nugent Score", 
                                   limits = c(1, 2.1), breaks = seq(1, 2, 0.5),
                                   labels = c('No BV', 'Inter BV', 'BV')) +
                theme(axis.text.y = element_text(angle = 45, vjust = 1, hjust = 1)) +
                geom_vline(xintercept = vertical.lines) +
                annotate("text", x = text1, y = 2.1, label = "MET Treatment", size = 3) +
                ggtitle(paste("Correlation between Nugent Score and Day - ", i, sep = "")) +
                labs(shape = "Treatment Status") 
    
    c <- ggplot(data2, aes(x = SERIAL, y = CST_Code)) +
                geom_point(size = 3, aes(colour = NUGENT_CLASS, shape = TREATMENT)) +
                scale_color_manual(values = Colors) +
                guides(colour = guide_legend(title = "Nugent Score")) +
                scale_x_continuous(name = "Day") +
                xlim(lower, upper) +
                xlim(lower, upper) +
                scale_y_continuous(name = "CST Group ID", 
                                   limits = c(1, 2.1), breaks = seq(1, 2, 0.25),
                                   labels = c('I', 'III', 'IV-A', 'IV-B', 'IV-C')) +
                geom_vline(xintercept = vertical.lines) +
                annotate("text", x = text1, y = 2.1, label = "MET Treatment", size = 3) +
                ggtitle(paste("Correlation between CST and Day - ", i, sep = "")) +
                labs(shape = "Treatment Status") 
    
  } else {
    met.date1 <- (min(met.rows, na.rm = TRUE)) - 0.5
    met.date2 <- met.date1 + 7
    met.date4 <- (max(met.rows, na.rm = TRUE)) + 0.5
    met.date3 <- met.date4 - 7
    vertical.lines <- c(met.date1, met.date2, met.date3, met.date4) 
    text1 <- met.date1 + 3.5
    text2 <- met.date3 + 3.5
    
    ## Plot ##
    n <- ggplot(data2, aes(x = SERIAL, y = Nugent_Code)) +
                geom_point(size = 3, aes(colour = CST, shape = TREATMENT)) +
                scale_color_manual(values = Colors) +
                guides(colour = guide_legend(title = "CST Group")) +
                scale_x_continuous(name = "Day") +
                xlim(lower, upper) +
                scale_y_continuous(name = "Nugent Score", 
                                   limits = c(1, 2.1), breaks = seq(1, 2, 0.5),
                                   labels = c('No BV', 'Inter BV', 'BV')) +
                theme(axis.text.y = element_text(angle = 45, vjust = 1, hjust = 1)) +
                geom_vline(xintercept = vertical.lines) +
                annotate("text", x = text1, y = 2.1, label = "MET Treatment", size = 3) +
                annotate("text", x = text2, y = 2.1,  label = "MET Treatment", size = 3) +
                ggtitle(paste("Correlation between Nugent Score and Day - ", i, sep = "")) +
                labs(shape = "Treatment Status") 
    
    c <- ggplot(data2, aes(x = SERIAL, y = CST_Code)) +
                geom_point(size = 3, aes(colour = NUGENT_CLASS, shape = TREATMENT)) +
                scale_color_manual(values = Colors) +
                guides(colour = guide_legend(title = "Nugent Score")) +
                scale_x_continuous(name = "Day") +
                xlim(lower, upper) +
                scale_y_continuous(name = "CST Group ID", 
                                   limits = c(1, 2.1), breaks = seq(1, 2, 0.25),
                                   labels = c('I', 'III', 'IV-A', 'IV-B', 'IV-C')) +
                geom_vline(xintercept = vertical.lines) +
                annotate("text", x = text1, y = 2.1, label = "MET Treatment", size = 3) +
                annotate("text", x = text2, y = 2.1,  label = "MET Treatment", size = 3) +
                ggtitle(paste("Correlation between CST and Day - ", i, sep = "")) +
                labs(shape = "Treatment Status") 
  }
  
  
  ## Create jitter plot
  j <- ggplot(data2, aes(x = Nugent_Code, y = CST_Code)) +
              geom_jitter(size = 2, width = 0.1, aes(colour = CST)) +
              scale_color_manual(values = Colors) +
              guides(colour = guide_legend(title = "CST")) +
              scale_x_continuous(name = "Nugent Score",
                                 limits = c(0.8, 2.2), breaks = seq(1, 2, 0.5),
                                 labels = c('No BV', 'Intermediate BV', 'BV')) +
              scale_y_continuous(name = "CST",
                                 limits = c(0.8, 2.2), breaks = seq(1, 2, 0.25),
                                 labels = c('I', 'III', 'IV-A', 'IV-B', 'IV-C')) +
              ggtitle(paste("Correlation between CST and Day - ", i, sep = ""))
  
  
  ## Create a contingency table from cross-classifying factors
  ass <- xtabs(~CST + NUGENT_CLASS, data = data2)
  ass.mat <- as.matrix(ass)
  
  # Create association plot
  a <- assocplot(ass.mat, col = c("lightblue", "tomato"),
                 xlab = "CST", ylab = "Infection Status", 
                 main = (paste("Correlation between CST and infection status - ", i, sep = "")))
  
  ## Save plots
  pdf(file = paste( i, "_Corr_Plots", ".pdf", sep = ""), width = 6, height = 6, onefile = TRUE)
    print(j)
    print(a)
  dev.off()
  
  # Combine day series plots
  Figs <- ggarrange(n, c,
                    ncol = 1, nrow = ,
                    legend = "right",
                    common.legend = FALSE)
  
  ggsave(Figs, filename = paste( i, "_NugentScore_CST_Day_Series", ".pdf", sep = ""), 
         width = 12, height = 8, dpi = 600)

}





