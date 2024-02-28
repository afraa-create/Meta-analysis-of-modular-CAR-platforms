# This is Systematic Review data analysis
# Date:29.01.2024

# part(4) data analysis-median survival
  
# load libraries
  library(readr)
  library(dmetar)
  library(metafor)
  library(tidyverse)
  library(dplyr)

# Load and clean data
  median_survival <- read_csv("01_tidy_data/Median_survival.csv", na = "NA") 
  median_survival <- janitor::clean_names(median_survival)
  median_survival <- median_survival %>% 
    mutate(car_design = as.factor(car_design),
           study_id = as.factor(study_id))
  median_survival_1 <- median_survival[median_survival$size_NC != 0, ]
    # table contains median values of each group in studies with the sample size
  
# calculate hazard rate from median survival for each group
  median_survival$hazard_rate_ex_1 <- log(2) / median_survival$effect_ex_1
  median_survival$hazard_rate_nc_1 <- log(2) / median_survival$effect_nc_1
  median_survival$hazard_rate_nc_2 <- log(2) / median_survival$effect_nc_2
  median_survival$hazard_rate_nc_3 <- log(2) / median_survival$effect_nc_3
  median_survival$hazard_rate_pc_1 <- log(2) / median_survival$effect_pc_1

# combine hazard rates within one study
  median_survival$Hazard_rate_PC <- median_survival$hazard_rate_pc_1
  median_survival$Hazard_rate_EX <- median_survival$hazard_rate_ex_1
  median_survival$Nxnc_1 <- median_survival$no_nc_1 * median_survival$hazard_rate_nc_1
  median_survival$Nxnc_2 <- median_survival$no_nc_2 * median_survival$hazard_rate_nc_2
  median_survival$Nxnc_3 <- median_survival$no_nc_1 * median_survival$hazard_rate_nc_3
  
  sum_NC <- c("Nxnc_1", "Nxnc_2", "Nxnc_3")
  sum_size <- c("no_nc_1", "no_nc_2", "no_nc_3")
  
  median_survival$Hazard_rate_NC <- rowSums(median_survival[sum_NC], na.rm = TRUE) / rowSums(median_survival[sum_size], na.rm = TRUE)
  
  median_survival$size_NC <- rowSums(median_survival[sum_size], na.rm = TRUE)
  median_survival$size_PC <- median_survival$no_pc_1
  median_survival$size_EX <- median_survival$no_ex_1
  
  median_survival$MSR <- median_survival$Hazard_rate_NC / median_survival$Hazard_rate_EX
  
  sum_weight <- c("size_EX", "size_NC")
  median_survival$weight <- rowSums(median_survival[sum_weight], na.rm = TRUE)
  
  median_survival$yi <- log(median_survival$MSR)
  median_survival$vi <- sqrt(1/(median_survival$weight))
  median_survival$Lower <-median_survival$yi - 1.96*sqrt(median_survival$vi)
  median_survival$Upper <-median_survival$yi + 1.96*sqrt(median_survival$vi)
  

  
  median_survival_1$yi <- log(median_survival_1$MSR)
  median_survival_1$vi <- sqrt(1/(median_survival_1$weight))
  median_survival_1$vi2 <- 1/(median_survival_1$weight^2)
  median_survival_1$Lower <-median_survival_1$yi - 1.96*sqrt(median_survival_1$vi2)
  median_survival_1$Upper <-median_survival_1$yi + 1.96*sqrt(median_survival_1$vi2)
  
# Perform a meta-analysis
  
  meta_result <- rma(yi = yi, vi  = vi2, method = "REML", data = median_survival_1)
  
  meta_result1 <- rma(yi = yi, vi  = vi2, method = "REML", data = median_survival_1, mods = ~ car_design)
  
  summary_result <- summary(meta_result)
  
  forest(meta_result, slab = median_survival_1$study_id, header= c("Author(s) and Year", "Log[MSR]   95%CI"))
  
  m.gen <- metagen(TE = yi,
                   lower = Lower,
                   upper = Upper,
                   n.e = size_EX,
                   n.c = size_NC,
                   studlab = study_id,
                   data = median_survival_1,
                   sm = "HR",
                   fixed = FALSE,
                   random = TRUE,
                   method.tau = "REML",
                   hakn = TRUE,
                   prediction = TRUE,)
  png(file = "03_plots/forestplot10.png", width = 5200, height = 3200, res = 300)
  forest(m.gen, sortvar = TE,
         prediction = TRUE, 
         print.tau2 = TRUE,
         layout = "RevMan5",
         leftcols = c("studlab","car_design", "TE", "seTE", "n.e", "n.c", "effect", "ci","w.random"),
                  digits.sd = 2,
         digits.tau2 = 2,
         colgap = "0.8cm",
         colgap.forest = "2cm",
         col.by = "black",
         col.square = "#CD534CFF",
         col.inside = "black",
         col.square.lines = "#CD534CFF",
         label.left = "Favors Control",
         label.right = "Favors Experimental",
         leftlabs = c("Author, year","Modular CAR Platform Name", "logMSR", NA, "Experimental Total", "Control Total", "MSR", NA, NA),
         smlab = "Median Survival Ratio 
IV. Random, 95% CI",
         )
  dev.off() 
#positive ctrl---------------------------------------------------------------  

  median_survival_2 <-median_survival[, c("study_id","car_design", "cell_line", "Hazard_rate_PC", "Hazard_rate_EX",
                                            "size_PC", "size_EX")]
  median_survival_2 <- median_survival[!is.na(median_survival_2$Hazard_rate_PC), ]
  
  median_survival_2$MSR <- median_survival_2$Hazard_rate_PC / median_survival_2$Hazard_rate_EX
  sum_weight <- c("size_EX", "size_PC")
  median_survival_2$weight <- rowSums(median_survival_2[sum_weight], na.rm = TRUE)  
  median_survival_2$yi <- log(median_survival_2$MSR)
  median_survival_2$vi <- sqrt(1/(median_survival_2$weight))
  median_survival_2$vi2 <- 1/(median_survival_2$weight^2)
  median_survival_2$Lower <-median_survival_2$yi - 1.96*sqrt(median_survival_2$vi2)
  median_survival_2$Upper <-median_survival_2$yi + 1.96*sqrt(median_survival_2$vi2)
  
  m.gen_1 <- metagen(TE = yi,
                   lower = Lower,
                   upper = Upper,
                   n.e = size_EX,
                   n.c = size_PC,
                   studlab = study_id,
                   data = median_survival_2,
                   sm = "HR",
                   fixed = FALSE,
                   random = TRUE,
                   method.tau = "REML",
                   hakn = TRUE,
                   prediction = TRUE,)
  png(file = "03_plots/forestplot11.png", width = 5200, height = 3200, res = 300)
  forest(m.gen_1, sortvar = TE,
         prediction = TRUE, 
         print.tau2 = TRUE,
         layout = "RevMan5",
         digits.sd = 2,
         digits.tau2 = 2,
         colgap = "0.8cm",
         colgap.forest = "2cm",
         col.by = "black",
         col.square = "#CD534CFF",
         col.inside = "black",
         col.square.lines = "#CD534CFF",
         label.left = "Worse Survival",
         label.right = "Better Survival",
         leftlabs = c(NA, "logMSR", NA, NA, NA, NA , "Median Survival Ratio 
IV. Random, 95% CI"),
         smlab = "Median Survival Ratio 
IV. Random, 95% CI")
  dev.off() 
 # subgroup analysis ----------------------------------------------
  m.gen_2 <- metagen(TE = yi,
                   lower = Lower,
                   upper = Upper,
                   n.e = size_EX,
                   n.c = size_NC,
                   studlab = study_id,
                   data = median_survival_1,
                   sm = "HR",
                   fixed = FALSE,
                   random = TRUE,
                   method.tau = "REML",
                   hakn = TRUE,
                   prediction = TRUE,
                   byvar = car_design,
                   print.byvar = FALSE)
  png(file = "03_plots/forestplot12-2.png", width = 4500, height = 3600, res = 300)
  forest(m.gen_2, sortvar = TE,
         prediction = TRUE, 
         print.tau2 = TRUE,
         layout = "RevMan5",
         digits.sd = 2,
         digits.tau2 = 2,
         colgap = "0.8cm",
         colgap.forest = "2cm",
         col.by = "black",
         col.square = "#CD534CFF",
         col.inside = "black",
         col.square.lines = "#CD534CFF",
         leftlabs = c(NA, "logMSR", NA, NA, NA, NA , "Median Survival Ratio 
IV. Random, 95% CI"),
         smlab = "Median Survival Ratio 
IV. Random, 95% CI",
         test.effect.subgroup.random = FALSE,
         overall = TRUE,
         overall.hetstat = TRUE,
         test.subgroup.random = FALSE,
         random.subgroup =FALSE,
         prediction.subgroup = FALSE,
         print.Q.subgroup =FALSE,
         subgroup.hetstat = FALSE,
         label.left = "Favors Control",
         label.right = "Favors Experimental",
         lab.e = "Experimental",
         lab.c = "Control", pooled.totals =FALSE)
  dev.off() 
  
  m.gen_3 <- metagen(TE = yi,
                     lower = Lower,
                     upper = Upper,
                     n.e = size_EX,
                     n.c = size_PC,
                     studlab = study_id,
                     data = median_survival_2,
                     sm = "HR",
                     fixed = FALSE,
                     random = TRUE,
                     method.tau = "REML",
                     hakn = TRUE,
                     prediction = TRUE,byvar = car_design,
                     print.byvar = FALSE)
  png(file = "03_plots/forestplot13-2.png", width = 5200, height = 3200, res = 300)
  forest(m.gen_3, sortvar = TE,
         prediction = TRUE, 
         print.tau2 = TRUE,
         layout = "RevMan5",
         digits.sd = 2,
         digits.tau2 = 2,
         colgap = "0.8cm",
         colgap.forest = "2cm",
         col.by = "black",
         col.square = "#CD534CFF",
         col.inside = "black",
         col.square.lines = "#CD534CFF",
         leftlabs = c(NA, "logMSR", NA, NA, NA, NA , "Median Survival Ratio 
IV. Random, 95% CI"),
         smlab = "Median Survival Ratio 
IV. Random, 95% CI",
  test.effect.subgroup.random = FALSE,
  overall = TRUE,
  overall.hetstat = TRUE,
  test.subgroup.random = FALSE,
  random.subgroup =FALSE,
  prediction.subgroup = FALSE,
  print.Q.subgroup =FALSE,
  subgroup.hetstat = FALSE,
  label.left = "Favors Control",
  label.right = "Favors Experimental",
  lab.e = "Experimental",
  lab.c = "Control", pooled.totals =FALSE)
  dev.off() 
  
  # Heterogeneity -------------------------------------------
  
  # add prediction interval
  m.cont4 <- update.meta(m.cont2, prediction = TRUE)
  
  # Basic outlier removal
  library(dmetar)
  
  find.outliers(m.gen)
  find.outliers(m.gen_1)
  
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
  g.rma <- rma(yi = m.gen$TE,
               sei = m.gen$seTE,
               method = m.gen$method.tau,
               test = "knha")
  
  g.gosh <- gosh(g.rma) 
  
  g.gosh.diag <- gosh.diagnostics(g.gosh, 
                                    km.params = list(centers = 2),
                                    db.params = list(eps = 0.08, 
                                                     MinPts = 50))
  res.gosh.diag
  
  m.gen.new <- update(m.gen, exclude = c(1,6,9)) 
  summary( m.gen.new)
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
  
  m.cont.new1 <- update(m.cont2, exclude = c(1,5,8))

  summary( m.cont.new1)
 
  
  # Publication Bias --------------------------------------------------------
  
  # Produce funnel plot
  funnel(m.cont,
         studlab = TRUE)
  
  # Define fill colors for contour
  col.contour = c("gray75", "gray85", "gray95")
  
  # Generate funnel plot (we do not include study labels here)
  funnel(m.cont,
         contour = c(0.9, 0.95, 0.99),
         col.contour = col.contour)
  
  # Add a legend
  legend(x = 0.6, y = 0.01, 
         legend = c("p < 0.1", "p < 0.05", "p < 0.01"),
         fill = col.contour)
  # Egger’s Regression Test
  
  metabias(m.cont, method.bias = "linreg")  