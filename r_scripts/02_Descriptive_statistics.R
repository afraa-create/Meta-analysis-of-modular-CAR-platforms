# This is Systematic Review data analysis
# Date:14.01.2024
# part(1) data wrangling

# Load data ---------------------------

library(readr)
study_characteristics <- read_csv("01_tidy_data/study_characteristics.csv",na = "NA")

# Load libraries ---------------------------  

library(r4np)

library(dplyr)

library(naniar)

library(forcats)

library(ggplot2)

library (tidyr)
# Wrangle data ---------------------------

glimpse(study_characteristics) # view each variable in the data set and their data types

skimr::skim(study_characteristics) # to skim the data set

study_characteristics <- janitor::clean_names(study_characteristics) # clean column names

study_characteristics %>% 
  rename( secondary_study_location = collaborators )
study_characteristics$secodary_study_location <- study_characteristics$collaborators # rename column

# change column type from character to factor
study_characteristics <-
  study_characteristics %>%
  mutate(study_id = as.factor(study_id),
         journal = as.factor(journal),
         primary_study_location = as.factor(primary_study_location),
         secodary_study_location = as.factor(secodary_study_location),
         modular_car_platform = as.factor(modular_car_platform),
         car_generation = as.factor(car_generation),
         car_ecd = as.factor(car_ecd),
         car_hinge = as.factor(car_hinge),
         car_tm = as.factor(car_tm),
         car_costim_1 = as.factor(car_costim_1),
         car_costim_2 = as.factor(car_costim_2),
         soluble_module = as.factor(soluble_module),
         target_antigen_1 = as.factor(target_antigen_1),
         target_antigen_2 = as.factor(target_antigen_2),
         target_antigen_3 = as.factor(target_antigen_3),
         tumor_type = as.factor(tumor_type),
         frequency_car = as.factor(frequency_car),
         frequency_soluble_module_dose_s = as.factor(frequency_soluble_module_dose_s),
         mice_strain = as.factor(mice_strain),
         mice_provider = as.factor(mice_provider),
         mice_sex = as.factor(mice_sex),
         xenograft_type = as.factor(xenograft_type),
         route_administration = as.factor(route_administration),
         selection_bias_sequence_generation = as.factor(selection_bias_sequence_generation),
         selection_bias_baseline_characteristics = as.factor(selection_bias_baseline_characteristics),
         selection_bias_allocation_concealment = as.factor(selection_bias_allocation_concealment),
         performance_bias_random_housing = as.factor(performance_bias_random_housing),
         performance_bias_blinding = as.factor(performance_bias_blinding),
         detection_bias_random_outcome_assessment = as.factor(detection_bias_random_outcome_assessment),
         detection_bias_blinding = as.factor(detection_bias_blinding),
         attrition_bias_incomplete_outcome_data = as.factor(attrition_bias_incomplete_outcome_data),
         reporting_bias_selective_outcome_reporting = as.factor(reporting_bias_selective_outcome_reporting),
         other_other_sources_of_bias = as.factor(other_other_sources_of_bias)
  )

# change values in column
study_characteristics <-
  study_characteristics %>% 
  mutate(mice_sex = fct_recode(mice_sex, "Female"="F", "Male"="M", "Male and Female"="Both"),
         mice_sex = replace_na(as.character(mice_sex), "Not Reported"),
         mice_sex = as_factor(mice_sex )
  )
# part(2) descriptive statistics-----------------------------------------------------

library(ggplot2)
library(ggpubr)
library(ggsci)
library(scales)
# Pie charts: Tumor type and mice sex--------------------------------------------------------------------- 

df_1 <- study_characteristics %>%
  group_by(mice_sex) %>%
  summarise(counts = n())
df_1 <- df_1 %>%
  arrange(desc(mice_sex)) %>%
  mutate(prop = round(counts*100/sum(counts), 1),
         lab.ypos = cumsum(prop) - 0.5*prop)

p_1 <- ggplot(df_1, aes(x = "", y = prop, fill = mice_sex)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  geom_text(aes(y = lab.ypos, label = prop), color = "white")+
  coord_polar("y", start = 0)+
  ggpubr::fill_palette("jco")+
  theme_void() + labs(fill = "Mice Gender")

# save plot high resolution

ggsave(filename = "E:/Review/Modular-CAR-Meta-Analysis/mice_gender.jpg", plot = p_1, width = 8, height = 6, dpi = 300)

df_2 <- study_characteristics %>%
  group_by(tumor_type) %>%
  summarise(counts = n())
df_2 <- df_2 %>%
  arrange(desc(tumor_type)) %>%
  mutate(prop = round(counts*100/sum(counts), 1),
         lab.ypos = cumsum(prop) - 0.5*prop)

p_2 <- ggplot(df_2, aes(x = "", y = prop, fill = tumor_type)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  geom_text(aes(y = lab.ypos, label = prop), color = "white")+
  coord_polar("y", start = 0)+
  ggpubr::fill_palette("jco")+
  theme_void() + labs(fill = "Tumor Type")

# save plot high resolution

ggsave(filename = "E:/Review/Modular-CAR-Meta-Analysis/tumor_type.jpg", plot = p_2, width = 8, height = 6, dpi = 300) 


# Plot outcome measures stacked bar chart---------------------------------------------------------------------

study_characteristics$Outcome1_basic_car_activity <- study_characteristics$basic_car_activity
study_characteristics$Outcome2_tumor_burden <- study_characteristics$tumor_burden
study_characteristics$Outcome3_tumor_volume <- study_characteristics$tumor_volume
study_characteristics$Outcome4_percent_servival <- study_characteristics$percent_servival
study_characteristics$Outcome5_percent_tumor_free <- study_characteristics$ percent_tumor_free
study_characteristics$Outcome6_percent_car_t_pb <- study_characteristics$percent_car_t_pb
study_characteristics$Outcome7_body_weight <- study_characteristics$body_weight
study_characteristics$Outcome8_serum_levels_soluble_module <- study_characteristics$serum_levels_soluble_module
study_characteristics$Outcome9_pb_t_phenotype <- study_characteristics$pb_t_phenotype
study_characteristics$Outcome10_human_cytokines_pb<- study_characteristics$human_cytokines_pb
study_characteristics$Outcome11_metastases <- study_characteristics$metastases
study_characteristics$Outcome12_pharmacokinetics <- study_characteristics$pharmacokinetics

library(tidyverse)

df_3 <- study_characteristics[, c("Outcome1_basic_car_activity", "Outcome2_tumor_burden",
                                  "Outcome3_tumor_volume", "Outcome4_percent_servival", "Outcome5_percent_tumor_free",
                                  "Outcome6_percent_car_t_pb", "Outcome7_body_weight", "Outcome8_serum_levels_soluble_module",
                                  "Outcome9_pb_t_phenotype", "Outcome10_human_cytokines_pb", "Outcome11_metastases", "Outcome12_pharmacokinetics")]

df_3 <-
  df_3 %>% 
    mutate(Outcome1_basic_car_activity = replace_na(as.character(Outcome1_basic_car_activity), "Not Applicable"),
           Outcome1_basic_car_activity = as_factor(Outcome1_basic_car_activity ),
           Outcome2_tumor_burden = as.factor(Outcome2_tumor_burden),
           Outcome3_tumor_volume = as.factor(Outcome3_tumor_volume),
           Outcome4_percent_servival = as.factor(Outcome4_percent_servival),
           Outcome5_percent_tumor_free = as.factor(Outcome5_percent_tumor_free),
           Outcome6_percent_car_t_pb = as.factor(Outcome6_percent_car_t_pb),
           Outcome7_body_weight = as.factor(Outcome7_body_weight),
           Outcome8_serum_levels_soluble_module = as.factor(Outcome8_serum_levels_soluble_module),
           Outcome9_pb_t_phenotype = as.factor(Outcome9_pb_t_phenotype),
           Outcome10_human_cytokines_pb = as.factor(Outcome10_human_cytokines_pb),
           Outcome11_metastases = as.factor(Outcome11_metastases),
           Outcome12_pharmacokinetics = as.factor(Outcome12_pharmacokinetics) )                                


# Function to calculate percentage of each value for each column
calc_percentage <- function(column) {
  column <- factor(column, levels = c("TRUE", "FALSE", "Not Applicable"))
  prop.table(table(column)) * 100
}

percentage_df_3 <- data.frame(lapply(df_3, calc_percentage))
percentage_df_31 <- data.frame(lapply(df_3, calc_percentage))
percentage_df_31 %>% select(where(is.double))
                            
df_44 <- as.data.frame(t(percentage_df_31)) %>%
  rownames_to_column("Outcomes") %>%
  rename_with(~str_replace(., "Outcome", ""), starts_with("Outcome")) %>%
  pivot_longer(cols = -Outcomes, names_to = "status", values_to = "percentage") %>%
  mutate(Outcomes = factor(Outcomes))

df_4 <-
   percentage_df_3 %>% select(where(is.double))
df_4 <- df_4 %>% 
  rename( Outcome_1 = Outcome1_basic_car_activity.Freq,
          Outcome_2 = Outcome2_tumor_burden.Freq,
          Outcome_3 = Outcome3_tumor_volume.Freq,
            Outcome_4 = Outcome4_percent_servival.Freq, 
            Outcome_5 = Outcome5_percent_tumor_free.Freq,
            Outcome_6 = Outcome6_percent_car_t_pb.Freq,
            Outcome_7 = Outcome7_body_weight.Freq,
            Outcome_8 = Outcome8_serum_levels_soluble_module.Freq,
            Outcome_9 = Outcome9_pb_t_phenotype.Freq,
            Outcome_10 = Outcome10_human_cytokines_pb.Freq,
            Outcome_11 = Outcome11_metastases.Freq,
            Outcome_12 = Outcome12_pharmacokinetics.Freq)

df_4<- t(df_4)
df_4<- as.data.frame(df_4)
library(tibble)
df_4<- rownames_to_column(df_4, var = "Outcomes")
df_4 <- df_4 %>% 
  rename(Reported = V1,
         Not_Reported = V2,
         Not_Applicable = V3)
df_4 <- df_4 %>%
  rename (Not_applicable = NOt_Applicable)
library(tidyr)
df_long <- pivot_longer(df_4, cols = -Outcomes, names_to = "status", values_to = "percentage")
df_long <- df_long %>% 
  mutate(Outcomes= as.factor(Outcomes))

new_outcome_names <- c("Baseline CAR Activity", "Human Cytokines in the Peripheral Blood","Metastases Formation","Pharmacokinetics Properties of the Soluble Module"," Tumor Burden",
                       "Tumor Volume", "Median Survival", "Percent Tumor Free", "Peripheral Blood CAR-T Cells Quantification", "Body Weight",
                       "Serum Levels of Soluble Module Post-implant", "Peripheral Blood T-cells Phenotype")
status_colors <- c("Reported" = "#7AA6DCFF",
                  "Not_Reported" = "#CD534CFF",
                  "Not_applicable" = "#868686FF")
# Create stacked bar chart
p_3 <-ggplot(df_long, aes(x = Outcomes, y = percentage, fill = status)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(labels = new_outcome_names) + 
  scale_fill_manual(values = status_colors,name = "Status", labels= c("Not Applicable", "Not Reported","Reported") )+
  labs(x = NULL, y = NULL) +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9, color = "black")
        ) +
  coord_flip()

ggsave(filename = "E:/Review/Modular-CAR-Meta-Analysis/Outcome_measures.jpg", plot = p_3, width = 8, height = 6, dpi = 300)
  
# Plot ROB stacked bar chart---------------------------------------------------------------------

df_5 <- study_characteristics[, c("selection_bias_sequence_generation", "selection_bias_baseline_characteristics",
                                  "selection_bias_allocation_concealment", "performance_bias_random_housing", "performance_bias_blinding",
                                  "detection_bias_random_outcome_assessment", "detection_bias_blinding", "attrition_bias_incomplete_outcome_data",
                                  "reporting_bias_selective_outcome_reporting", "other_other_sources_of_bias")]


# Function to calculate percentage of each value for each column
calc_percentage_1 <- function(column) {
  column <- factor(column, levels = c("High", "Low", "Unclear"))
  prop.table(table(column)) * 100
}

percentage_df_4 <- data.frame(lapply(df_5, calc_percentage_1))

df_6 <-
  percentage_df_4 %>% select(where(is.double))

df_6<- t(df_6)

df_6<- as.data.frame(df_6)

df_6<- rownames_to_column(df_6, var = "ROB")

df_6 <- df_6 %>% 
  rename(High = V1,
         Low = V2,
         Unclear = V3)

df_long_1 <- pivot_longer(df_6, cols = -ROB, names_to = "ROB2", values_to = "percentage")
df_long_1 <- df_long_1 %>% 
  mutate(ROB2= as.factor(ROB2))

new_ROB_names <- c("Incomplete Outcome Data", "Detection Bias - Blinding", "Random Outcome Assessment",
                   "Other Sources of Bias", " Performance Bias - Blinding", "Random Housing", "Selective Outcome Reporting",
                   "Allocation Concealment", "Baseline Characteristics", "Sequence Generation")
status_colors_2 <- c("High" = "#CD534CFF",
                   "Low" = "#5CB85CFF",
                   "Unclear" = "#EEA236FF")

# Create stacked bar chart
p_4 <-ggplot(df_long_1, aes(x = ROB, y = percentage, fill = ROB2)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(labels = new_ROB_names) +
  scale_fill_manual(values = status_colors_2,name = "Risk of Bias", labels= c("High", "Low", "Unclear") )+
  labs(x = NULL, y = NULL) +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9, color = "black")) +
  coord_flip()

ggsave(filename = "E:/Review/Modular-CAR-Meta-Analysis/ROB.jpg", plot = p_4, width = 8, height = 6, dpi = 300)


p_5 <- ggarrange(p_1, p_2, p_3, p_4, 
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)

ggsave(filename = "E:/Review/Modular-CAR-Meta-Analysis/all_characteristics.jpg", plot = p_5, width = 14, height = 7, dpi = 300)


# Statistics---------------------------------------------------------------------------------------------------
study_characteristics %>% 
  count(tumor_type)
mean(study_characteristics$mice_age,na.rm = TRUE)
mean(study_characteristics$mice_age,na.rm = TRUE)
median(study_characteristics$mice_age,na.rm = TRUE)
mean(study_characteristics$mice_number_total,na.rm = TRUE)
median(study_characteristics$mice_number_total,na.rm = TRUE)
sum(study_characteristics$mice_number_total,na.rm = TRUE)
summary(study_characteristics$immunogenicity_factors_report)
summary(study_characteristics$tumor_burden)

summary_table <- study_characteristics %>%
  summarise_if(is.logical, list(
    sum_true = ~sum(., na.rm = TRUE),
    sum_false = ~sum(!., na.rm = TRUE),
    num_na = ~sum(is.na(.))
  ))
summary_table <- t(summary_table)
 
# Other outcome measures (Body weight)---------------------------------------------------------------------------------------
body_weight <- read_csv("00_raw_data/07_Body_weight.csv",na = "NA")

write.csv(body_weight1, "01_tidy_data/Body_weight1.csv", row.names=FALSE)
per_body_weight <- read_csv("01_tidy_data/Body_weight1.csv",na = "NA")

# plot data with ready means and sd
weight_smd <- metacont(n.e = ex_size,
                       mean.e = ex_mean,
                       sd.e = ex_sd,
                       n.c = ctrl_n,
                       mean.c = ctrl_mean,
                       sd.c = ctrl_sd,
                       studlab = study_id,
                       data = body_weight,
                       sm = "SMD",
                       method.smd = "Hedges",
                       fixed = FALSE,
                       random = TRUE,
                       method.tau = "REML",
                       hakn = TRUE,
                       title = "Body weight")

body_weight$SMD <- weight_smd$TE
body_weight$ci_low <- weight_smd$lower
body_weight$ci_high <- weight_smd$upper
body_weight$index <- c(1:8)
body_weight$weight <- weight_smd$w.random

ploty1 <-ggplot(data=body_weight, aes(y=index, x=SMD, xmin=ci_low, xmax=ci_high)) +
  geom_errorbarh(height=0, size =10, alpha = 0.9 , show.legend = FALSE,  color = "#0073C2FF" ) +
  geom_segment(aes(x = SMD, xend = SMD+0.01, yend = index), size = 10) +
  scale_y_continuous(name = "", breaks=1:nrow(body_weight), labels=body_weight$study_id) +
  labs( x='Std Mean Difference, 95% CI') +
  theme_classic2() +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme(axis.text.y = element_text(size = 10, color = "black",face="bold"),
        axis.text.x = element_text(size = 10, color = "black", face="bold"),
        axis.title.x = element_text(size=12, face="bold", colour = "black"))

ggsave(filename = "E:/Review/Modular-CAR-Meta-Analysis/03_plots/body_weight.jpg",
       plot = ploty1, width = 8, height = 6, dpi = 300)

write.csv(per_body_weight, "01_tidy_data/Body_weight.csv", row.names=FALSE)
# Other outcome measures (Cytokine)---------------------------------------------------------------------------------------
cytokine <- read_csv("00_raw_data/06_Cytokine.csv",na = "NA")
cytokine_smd_n <- metacont(n.e = ex_n,
                       mean.e = ex_effect,
                       sd.e = ex_sd,
                       n.c = ctrl_n,
                       mean.c = ctrl_effect,
                       sd.c = ctrl_sd,
                       studlab = study_id,
                       data = cytokine,
                       sm = "SMD",
                       method.smd = "Hedges",
                       fixed = FALSE,
                       random = TRUE,
                       method.tau = "REML",
                       hakn = TRUE,
                       title = "Cytokine")

cytokine_smd_p <- metacont(n.e = ex_n,
                           mean.e = ex_effect,
                           sd.e = ex_sd,
                           n.c = ctrl_p_n,
                           mean.c = ctrl_p_effect,
                           sd.c = ctrl_p_sd,
                           studlab = study_id,
                           data = cytokine,
                           sm = "SMD",
                           method.smd = "Hedges",
                           fixed = FALSE,
                           random = TRUE,
                           method.tau = "REML",
                           hakn = TRUE,
                           title = "Cytokine")

cytokine$SMD_n <- cytokine_smd_n$TE
cytokine$ci_low_n <- cytokine_smd_n$lower
cytokine$ci_high_n <- cytokine_smd_n$upper
cytokine$index <- c(1:17)
cytokine$weight_n <-cytokine_smd_n$w.random
cytokine$SMD_p <- cytokine_smd_p$TE
cytokine$ci_low_p <- cytokine_smd_p$lower
cytokine$ci_high_p <- cytokine_smd_p$upper
cytokine$weight_p <-cytokine_smd_p$w.random
cytokine$Cytokine <- as.factor(cytokine$Cytokine)
cytokine[17, "study_id"] <- "Lu et al., 2019"

  

ploty2 <- ggplot(data=cytokine, aes(y=index, x=SMD_n, xmin=ci_low_n, xmax=ci_high_n )) +
  scale_color_manual(values = c("IFN-y" = "#CD534CFF", "IL-2" = "#5CB85CFF", "TNF-a" = "#EFC000FF" , 
                                "IL-18" = "#868686FF", "IL-4" = "#8F7700FF", "IL-6" ="#003C67FF")) +
  geom_errorbarh(height=0, size =10, alpha = 0.9 ,  aes(colour = Cytokine), show.legend = FALSE ) +
  geom_segment(aes(x = SMD_n, xend = SMD_n+0.01, yend = index), size = 10) +
  facet_grid2(study_id ~ Cytokine, scales = "free", switch = "y", independent = "y", shrink = TRUE)+
  labs( x='Std. Mean Difference, 95% CI') +
  theme_classic2() +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme (axis.title.y = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.y=element_blank(),
         axis.text.x = element_text(size = 10, color = "black"),
         axis.title.x = element_text(size=12, face="bold", colour = "black"),
         strip.placement = "outside",                      
         strip.background = element_blank(),
         strip.text.y = element_text(size = 10, colour = "black", face = "bold"),
         strip.text.x = element_text(size = 12, colour = "black", face = "bold"),
         strip.text.y.left = element_text(hjust = 1, angle = 0),
         panel.spacing = unit(0, "lines"))

ggsave(filename = "E:/Review/Modular-CAR-Meta-Analysis/03_plots/cytokine_negative.jpg",
       plot = ploty2, width = 14, height = 6, dpi = 300) 

write.csv(cytokine, "01_tidy_data/Cytokine.csv", row.names=FALSE)

cytokine_p <- cytokine[!is.na(cytokine$SMD_p), ]
cytokine_p$index <- c(1:10)
ploty3 <- ggplot(data=cytokine_p, aes(y=index, x=SMD_p, xmin=ci_low_p, xmax=ci_high_p)) +
  scale_color_manual(values = c("IFN-y" = "#CD534CFF", "IL-2" = "#5CB85CFF", "TNF-a" = "#EFC000FF" , 
                                "IL-18" = "#868686FF", "IL-4" = "#8F7700FF", "IL-6" ="#003C67FF")) +
  geom_errorbarh(height=0, size =10, alpha = 0.9 ,  aes(colour = Cytokine), show.legend = FALSE ) +
  geom_segment(aes(x = SMD_p, xend = SMD_p+0.01, yend = index), size = 10) +
  facet_grid2(study_id ~ Cytokine, scales = "free", switch = "y", independent = "y", shrink = TRUE)+
  labs( x='Std. Mean Difference, 95% CI') +
  theme_classic2() +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme (axis.title.y = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.y=element_blank(),
         axis.text.x = element_text(size = 10, color = "black"),
         axis.title.x = element_text(size=12, face="bold", colour = "black"),
         strip.placement = "outside",                      
         strip.background = element_blank(),
         strip.text.y = element_text(size = 10, colour = "black", face = "bold"),
         strip.text.x = element_text(size = 12, colour = "black", face = "bold"),
         strip.text.y.left = element_text(hjust = 1, angle = 0),
         panel.spacing = unit(0, "lines"))

ggsave(filename = "E:/Review/Modular-CAR-Meta-Analysis/03_plots/cytokine_positive.jpg",
       plot = ploty3, width = 10, height = 6, dpi = 300)


ploty4 <- ggarrange(ggarrange(ploty1, ploty3, ncol = 2, labels = c("A", "B")), 
                    ggarrange(ploty2, labels = "C"), 
                    nrow = 2)

ggsave(filename = "E:/Review/Modular-CAR-Meta-Analysis/03_plots/other_outcomes.jpg",
       plot = ploty4, width = 16, height = 12, dpi = 300)



outcome_columns <- grep("^Outcome", names(study_characteristics), value = TRUE)
df_33 <- study_characteristics %>%
  select(all_of(outcome_columns)) %>%
  mutate(across(everything(), ~replace_na(as.character(.), "Not Applicable"))) %>%
  mutate(across(everything(), as.factor))