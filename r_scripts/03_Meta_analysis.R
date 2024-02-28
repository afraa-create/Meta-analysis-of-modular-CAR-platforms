# This is Systematic Review data analysis
# Date:14.01.2024

# part(3) data analysis

# Load and clean data ---------------------------

devtools::install_github("MathiasHarrer/dmetar")

library(readr)

tumor_burden <- read_csv("01_tidy_data/Tumor_Burden.csv", na = "NA")
tumor_volume <- read_csv("01_tidy_data/Tumor_Volume.csv", na = "NA")

library(dmetar)

library(metafor)

library(tidyverse)

library(dplyr)
library(meta)

tumor_burden <- janitor::clean_names(tumor_burden)
tumor_volume <- janitor::clean_names(tumor_volume)

tumor_burden %>% 
  rename( sd_ex = sd_x )

tumor_burden <-
  tumor_burden  %>%
  mutate(study_id = as.factor(study_id),
         car_design = as.factor(car_design),
         car_hinge = as.factor(car_hinge),
         car_tm = as.factor(car_tm),
         car_costim_1 = as.factor(car_costim_1),
         soluble_module = as.factor(soluble_module),
         target_antigen_1 = as.factor(target_antigen_1),
         target_antigen_2 = as.factor(target_antigen_2),
         frequency_car = as.factor(frequency_car),
         frequency_soluble_module_dose_s = as.factor(frequency_soluble_module_dose_s),
         cell_line = as.factor(cell_line)
  )


tumor_burden <- tumor_burden[tumor_burden$sd_x != 0, ]
tumor_burden[2, "study_id"] <- "Cho et al., 2018"
tumor_burden_2 <- tumor_burden[tumor_burden$effect_pc != 0, ]

tumor_volume <-
  tumor_volume  %>%
  mutate(study_id = as.factor(study_id),
         car_design = as.factor(car_design),
         car_hinge = as.factor(car_hinge),
         car_tm = as.factor(car_tm),
         car_costim_1 = as.factor(car_costim_1),
         soluble_module = as.factor(soluble_module),
         target_antigen_1 = as.factor(target_antigen_1),
         target_antigen_2 = as.factor(target_antigen_2),
         frequency_car = as.factor(frequency_car),
         frequency_soluble_module_dose_s = as.factor(frequency_soluble_module_dose_s),
         cell_line = as.factor(cell_line)
  )

tumor_volume <- tumor_volume[complete.cases(tumor_volume$sd_ex), ]
levels(tumor_volume$study_id) <- c(levels(tumor_volume$study_id), "Stock et al., 2022")
tumor_volume[16, "study_id"] <- "Stock et al., 2022"

tumor_volume_2 <- tumor_volume[!is.na(tumor_volume$effect_pc), ]

# Meta analysis - tumor burden and volume---------------------------

library(esc)
m.cont <- metacont(n.e = no_ex,
                   mean.e = mean_ex,
                   sd.e = sd_x,
                   n.c = no_nc,
                   mean.c = mean_nc,
                   sd.c = sd_nc,
                   studlab = study_id,
                   data = tumor_burden,
                   sm = "SMD",
                   method.smd = "Hedges",
                   fixed = FALSE,
                   random = TRUE,
                   method.tau = "REML",
                   hakn = TRUE,
                   title = "Tumor Burden")

png(file = "03_plots/forestplot1_CAR name.png", width = 5600, height = 1600, res = 300)

forest(m.cont, 
       sortvar = TE,
       leftcols = c("studlab","car_design", "n.e", "mean.e", "sd.e", "n.c", "mean.c", "sd.c","effect", "ci", "w.random"),
       leftlabs = c("Author, year", "Modular CAR Platform Name", NA, NA, NA, NA, NA, NA,NA,NA,NA),layout = "RevMan5",
              digits.sd = 2,
       digits.tau2 = 2,  colgap = "0.8cm",
       colgap.forest = "2cm",col.by = "black",
       col.square = "#7AA6DCFF",
       col.inside = "black",
       col.square.lines = "#7AA6DCFF",
       label.left = "Favors Experimental",
       label.right = "Favors Control",
       lab.e = "Experimental",
       lab.c = "Negative Control",
       overall = TRUE,
       overall.hetstat = TRUE,
       prediction = TRUE, 
       print.tau2 = TRUE)
dev.off()

# subgroup = car design 
m2.cont <- metacont(n.e = no_ex,
                   mean.e = mean_ex,
                   sd.e = sd_x,
                   n.c = no_nc,
                   mean.c = mean_nc,
                   sd.c = sd_nc,
                   studlab = study_id,
                   data = tumor_burden,
                   sm = "SMD",
                   method.smd = "Hedges",
                   fixed = FALSE,
                   random = TRUE,
                   method.tau = "REML",
                   hakn = TRUE,
                   byvar = car_design,
                   print.byvar = FALSE,
                   title = "Tumor Burden")

png(file = "03_plots/forestplot2.png", width = 5200, height = 3200, res = 300)

forest(m2.cont, layout = "RevMan5",prediction = TRUE, 
       digits.sd = 2,
       digits.tau2 = 2,  colgap = "0.8cm",
       colgap.forest = "2cm",col.by = "black",
       col.square = "#7AA6DCFF",
       col.inside = "black",
       col.square.lines = "#7AA6DCFF",
       test.effect.subgroup.random = FALSE,
       overall = TRUE,
       overall.hetstat = TRUE,
       test.subgroup.random = FALSE,
       random.subgroup =FALSE,
       prediction.subgroup = FALSE,
       print.Q.subgroup =FALSE,
       subgroup.hetstat = FALSE,
       sep.subgroup = "10cm",
       label.left = "Favors Experimental",
       label.right = "Favors Control",
       lab.e = "Experimental",
       lab.c = "Negative Control", pooled.totals =FALSE)

dev.off()
png(file = "03_plots/forestplot2.png", width = 5200, height = 3200, res = 300)

forest(m2.cont, layout = "RevMan5",prediction = TRUE, 
       digits.sd = 2,
       digits.tau2 = 2,  colgap = "0.8cm",
       colgap.forest = "2cm",col.by = "black",
       col.square = "#7AA6DCFF",
       col.inside = "black",
       col.square.lines = "#7AA6DCFF",
       test.effect.subgroup.random = FALSE,
       overall = TRUE,
       overall.hetstat = TRUE,
       test.subgroup.random = FALSE,
       random.subgroup =FALSE,
       prediction.subgroup = FALSE,
       print.Q.subgroup =FALSE,
       subgroup.hetstat = FALSE,
       sep.subgroup = "10cm",
       label.left = "Favors Experimental",
       label.right = "Favors Control",
       lab.e = "Experimental",
       lab.c = "Negative Control", pooled.totals =FALSE)

dev.off()


#---------------------------------------------------

m.cont1 <- metacont(n.e = no_ex,
                   mean.e = mean_ex,
                   sd.e = sd_x,
                   n.c = no_pc,
                   mean.c = effect_pc,
                   sd.c = error_pc,
                   studlab = study_id,
                   data = tumor_burden_2,
                   sm = "SMD",
                   method.smd = "Hedges",
                   fixed = FALSE,
                   random = TRUE,
                   method.tau = "REML",
                   hakn = TRUE,
                   title = "Tumor Burden")

png(file = "03_plots/forestplot3_car name.png", width = 5500, height = 1600, res = 300)

forest(m.cont1, 
       sortvar = TE,
       leftcols = c("studlab","car_design", "n.e", "mean.e", "sd.e", "n.c", "mean.c", "sd.c","effect", "ci", "w.random"),
       leftlabs = c("Author, year", "Modular CAR Platform Name", NA, NA, NA, NA, NA, NA,NA,NA,NA),layout = "RevMan5",
       digits.sd = 2,
       digits.tau2 = 2,  colgap = "0.8cm",
       colgap.forest = "2cm",col.by = "black",
       col.square = "#7AA6DCFF",
       col.inside = "black",
       col.square.lines = "#7AA6DCFF",
       label.left = "Favors Experimental",
       label.right = "Favors Control",
       lab.e = "Experimental",
       lab.c = "Positive Control",
       overall = TRUE,
       overall.hetstat = TRUE,
       prediction = TRUE, 
       print.tau2 = TRUE)
dev.off()

# subgroup = car design 
m2.cont1 <- metacont(n.e = no_ex,
                     mean.e = mean_ex,
                     sd.e = sd_x,
                     n.c = no_pc,
                     mean.c = effect_pc,
                     sd.c = error_pc,
                     studlab = study_id,
                     data = tumor_burden_2,
                     sm = "SMD",
                     method.smd = "Hedges",
                     fixed = FALSE,
                     random = TRUE,
                     method.tau = "REML",
                     hakn = TRUE,
                    byvar = car_design,
                    print.byvar = FALSE,
                    test.subgroup = FALSE,
                    title = "Tumor Burden")

png(file = "03_plots/forestplot4-1.png", width = 5200, height = 3200, res = 300)

forest(m2.cont1, layout = "RevMan5",prediction = TRUE, 
       digits.sd = 2,
       digits.tau2 = 2,  colgap = "0.8cm",
       colgap.forest = "2cm",col.by = "black",
       col.square = "#7AA6DCFF",
       col.inside = "black",
       col.square.lines = "#7AA6DCFF",
       test.effect.subgroup.random = FALSE,
       overall = TRUE,
       overall.hetstat = TRUE,
       test.subgroup.random = FALSE,
       random.subgroup =FALSE,
       prediction.subgroup = FALSE,
       print.Q.subgroup =FALSE,
       subgroup.hetstat = FALSE,
       sep.subgroup = "10cm",
       label.left = "Favors Experimental",
       label.right = "Favors Control",
       lab.e = "Experimental",
       lab.c = "Positive Control", pooled.totals =FALSE)

dev.off()
png(file = "03_plots/forestplot4.png", width = 5200, height = 3200, res = 300)


forest(m2.cont1, layout = "RevMan5",prediction = TRUE, 
       digits.sd = 2,
       digits.tau2 = 2,  colgap = "0.8cm",
       colgap.forest = "2cm",col.by = "black",
       col.square = "#7AA6DCFF",
       col.inside = "black",
       col.square.lines = "#7AA6DCFF",
       test.subgroup.random = FALSE,
       random.subgroup =FALSE,
       prediction.subgroup = FALSE,
       print.Q.subgroup =FALSE,
       subgroup.hetstat = FALSE,
       overall = TRUE,
       overall.hetstat = TRUE,
       label.left = "Favors Experimental ",
       label.right = "Favors Control",
       lab.e = "Experimental",
       lab.c = "Positive Control")

dev.off()


# Volume----------------------------------------------------------------------------------

m.cont2 <- metacont(n.e = no_ex,
                   mean.e = mean_ex,
                   sd.e = sd_ex,
                   n.c = no_nc,
                   mean.c = mean_nc,
                   sd.c = sd_nc,
                   studlab = study_id,
                   data = tumor_volume,
                   sm = "SMD",
                   method.smd = "Hedges",
                   fixed = FALSE,
                   random = TRUE,
                   method.tau = "REML",
                   hakn = TRUE,
                   title = "Tumor Volume")

png(file = "03_plots/forestplot5_car name.png", width = 5500, height = 2000, res = 300)

forest(m.cont2, 
         sortvar = TE,
       leftcols = c("studlab","car_design", "n.e", "mean.e", "sd.e", "n.c", "mean.c", "sd.c","effect", "ci", "w.random"),
       leftlabs = c("Author, year", "Modular CAR Platform Name", NA, NA, NA, NA, NA, NA,NA,NA,NA),layout = "RevMan5",
       digits.sd = 2,
       digits.tau2 = 2,  colgap = "0.8cm",
       colgap.forest = "2cm",col.by = "black",
       col.square = "#EEA236FF",
       col.inside = "black",
       col.square.lines = "#EEA236FF",
       label.left = "Favors Experimental",
       label.right = "Favors Control",
       lab.e = "Experimental",
       lab.c = "Negative Control",overall = TRUE,
       overall.hetstat = TRUE,
       prediction = TRUE, 
       print.tau2 = TRUE)
dev.off()


m2.cont2 <- metacont(n.e = no_ex,
                     mean.e = mean_ex,
                     sd.e = sd_ex,
                     n.c = no_nc,
                     mean.c = mean_nc,
                     sd.c = sd_nc,
                     studlab = study_id,
                     data = tumor_volume,
                     sm = "SMD",
                     method.smd = "Hedges",
                     fixed = FALSE,
                     random = TRUE,
                     method.tau = "REML",
                     hakn = TRUE,
                     byvar = car_design,
                     print.byvar = FALSE,
                     title = "Tumor Volume")

png(file = "03_plots/forestplot6-1.png", width = 5200, height = 3200, res = 300)

forest(m2.cont2, layout = "RevMan5",prediction = TRUE, 
       digits.sd = 2,
       digits.tau2 = 2,  colgap = "0.8cm",
       colgap.forest = "2cm",col.by = "black",
       col.square = "#EEA236FF",
       col.inside = "black",
       col.square.lines = "#EEA236FF",
       test.effect.subgroup.random = FALSE,
       overall = TRUE,
       overall.hetstat = TRUE,
       test.subgroup.random = FALSE,
       random.subgroup =FALSE,
       prediction.subgroup = FALSE,
       print.Q.subgroup =FALSE,
       subgroup.hetstat = FALSE,
       label.left = "Favors Experimental",
       label.right = "Favors Control",
       lab.e = "Experimental",
       lab.c = "Negative Control", pooled.totals =FALSE)

dev.off()


png(file = "03_plots/forestplot6.png", width = 5200, height = 3200, res = 300)


forest(m2.cont2, layout = "RevMan5",prediction = TRUE, 
       digits.sd = 2,
       digits.tau2 = 2,  colgap = "0.8cm",
       colgap.forest = "2cm",col.by = "black",
       col.square = "#7AA6DCFF",
       col.inside = "black",
       col.square.lines = "#7AA6DCFF",
       test.effect.subgroup.random = TRUE,
       overall = TRUE,
       overall.hetstat = TRUE,
       label.left = "Tumor Volume Reduction",
       label.right = "Tumor Volume Increase",
       lab.e = "Experimental",
       lab.c = "Negative Control")

dev.off()


m.cont3 <- metacont(n.e = no_ex,
                    mean.e = mean_ex,
                    sd.e = sd_ex,
                    n.c = no_pc,
                    mean.c = effect_pc,
                    sd.c = error_pc,
                    studlab = study_id,
                    data = tumor_volume_2,
                    sm = "SMD",
                    method.smd = "Hedges",
                    fixed = FALSE,
                    random = TRUE,
                    method.tau = "REML",
                    hakn = TRUE,
                    title = "Tumor Burden")

png(file = "03_plots/forestplot7_car name.png", width = 5200, height = 2000, res = 300)

forest(m.cont3, 
       sortvar = TE,
       leftcols = c("studlab","car_design", "n.e", "mean.e", "sd.e", "n.c", "mean.c", "sd.c","effect", "ci", "w.random"),
       leftlabs = c("Author, year", "Modular CAR Platform Name", NA, NA, NA, NA, NA, NA,NA,NA,NA),layout = "RevMan5",
       digits.sd = 2,
       digits.tau2 = 2,  colgap = "0.8cm",
       colgap.forest = "2cm",col.by = "black",
       col.square = "#EEA236FF",
       col.inside = "black",
       col.square.lines = "#EEA236FF",
       label.left = "Favors Experimental",
       label.right = "Favors Control",
       lab.e = "Experimental",
       lab.c = "Positive Control",overall = TRUE,
       overall.hetstat = TRUE,
       prediction = TRUE, 
       print.tau2 = TRUE)
dev.off()




# Heterogeneity -------------------------------------------

# add prediction interval
m.cont4 <- update.meta(m.cont2, prediction = TRUE)

# Basic outlier removal
library(dmetar)

find.outliers(m.cont)
find.outliers(m.cont1)
find.outliers(m.cont2)
find.outliers(m.cont3)

# Influence analysis

m.cont.inf <- InfluenceAnalysis(m.cont, random = TRUE)
m.cont.inf2 <- InfluenceAnalysis(m.cont2, random = TRUE)
 
 # Leave-One-Out Meta-Analysis Results
 png(file = "03_plots/inf-es.png", width = 5200, height = 2000, res = 300)
plot(m.cont.inf, "es")
dev.off()
png(file = "03_plots/inf-i2.png", width = 5400, height = 2000, res = 300)
 plot(m.cont.inf, "i2")
 dev.off()
 png(file = "03_plots/inf-es2.png", width = 5200, height = 2000, res = 300)
 plot(m.cont.inf2, "es")
 dev.off()
 png(file = "03_plots/inf-i22.png", width = 5400, height = 2000, res = 300)
 plot(m.cont.inf2, "i2")
 dev.off()
# GOSH Plot Analysis (effect size−heterogeneity pattern)
 m.rma <- rma(yi = m.cont$TE,
              sei = m.cont$seTE,
              method = m.cont$method.tau,
              test = "knha")
 
 res.gosh <- gosh(m.rma) 
 
 res.gosh.diag <- gosh.diagnostics(res.gosh, 
                                   km.params = list(centers = 2),
                                   db.params = list(eps = 0.08, 
                                                    MinPts = 50))
 res.gosh.diag
 
 m.cont.new <- update(m.cont, exclude = c(1,6,9)) 
   summary( m.cont.new)
   summary(m.cont)
  
   m.rma2 <- rma(yi = m.cont2$TE,
                sei = m.cont2$seTE,
                method = m.cont2$method.tau,
                test = "knha")
   
   res.gosh1 <- gosh(m.rma2) 
   
   res.gosh.diag1 <- gosh.diagnostics(res.gosh1, 
                                     km.params = list(centers = 2),
                                     db.params = list(eps = 0.08, 
                                                      MinPts = 50))
   res.gosh.diag1
   
   m.cont.new1 <- update(m.cont1, exclude = c(1,5,8)) 
   summary( m.cont.new1)
   
# Publication Bias --------------------------------------------------------
   
   # Produce funnel plot
   funnel(m.cont,
               studlab = TRUE)
   
   # Define fill colors for contour
   col.contour = c("gray75", "gray85", "gray95")
   
   # Generate funnel plot (we do not include study labels here)
   png(file = "03_plots/funnel-tumor-burden1.png", width = 3000, height = 2000, res = 300)  
    funnel(m.cont,xlim = c(-15, 15),
               contour = c(0.9, 0.95, 0.99),
               col.contour = col.contour)
   
   # Add a legend
   legend(x = 11, y = 0.01, 
          legend = c("p < 0.1", "p < 0.05", "p < 0.01"),
          fill = col.contour)
   title("Contour-Enhanced Funnel Plot (Tumor Burden vs negative control)")
   
   dev.off()
   
   png(file = "03_plots/funnel-tumor-burden2.png", width = 3000, height = 2000, res = 300)  
   funnel(m.cont1,xlim = c(-5, 5),
          contour = c(0.9, 0.95, 0.99),
          col.contour = col.contour)
   
   # Add a legend
   legend(x = 3.8, y = 0.01, 
          legend = c("p < 0.1", "p < 0.05", "p < 0.01"),
          fill = col.contour)
   title("Contour-Enhanced Funnel Plot (Tumor Burden vs positive control)")
   
   dev.off()
   
   png(file = "03_plots/funnel-tumor-volume1.png", width = 3000, height = 2000, res = 300)  
   funnel(m.cont2,xlim = c(-10, 10),
          contour = c(0.9, 0.95, 0.99),
          col.contour = col.contour)
   
   # Add a legend
   legend(x = 7.8, y = 0.01, 
          legend = c("p < 0.1", "p < 0.05", "p < 0.01"),
          fill = col.contour)
   title("Contour-Enhanced Funnel Plot (Tumor Volume vs negative control)")
   
   dev.off()
   
 
   
# Egger’s Regression Test
   
  eggers.test(m.cont)
  metabias(m.cont, method.bias = "linreg")
   