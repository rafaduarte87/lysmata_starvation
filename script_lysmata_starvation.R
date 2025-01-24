# Loading the packages ---------------------------------------------------------

library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(performance)
library(emmeans)
library(patchwork)
library(vegan)
library(ggfortify)
library(ggcorrplot)
library(ggtext)
library(pairwiseAdonis)

# Importing the datasets -------------------------------------------------------

# Experiment 1 - biomarkers

data_experiment_1 <-
  read_excel("lysmata_starvation_dataset.xlsx", sheet = 1,
             range = "A1:M151", na = "NA") |> 
  mutate(food = factor(food),
         time = factor(time))

data_experiment_1$food_time <-
  as.factor(paste(data_experiment_1$food,
            data_experiment_1$time,
            sep = "_"))

# let's remove the data for female 4 from the treatment S24 and for female 5 from the treatment F24 - values are all strange

data_experiment_1_edited <-
  data_experiment_1 |> 
  filter(!sample %in% c(66:70, 96:100))

# Experiment 2 - biomarkers

data_experiment_2 <-
  read_excel("lysmata_starvation_dataset.xlsx", sheet = 2,
             range = "A1:M181", na = "NA") |> 
  mutate(stage = factor(stage,
                        levels = c("early_pl", "late_pl", "z9")),
         food = factor(food))

data_experiment_2$food_stage <-
  as.factor(paste(data_experiment_2$food,
                  data_experiment_2$stage,
                  sep = "_"))

# Pooled data - biomarkers

data_nmds <-
  read_excel("lysmata_starvation_dataset.xlsx", sheet = 3,
             range = "A1:M331", na = "NA") |> 
  mutate(stage_1 = factor(stage_1,
                        levels = c("z2", "z9", "early_pl", "late_pl")),
         stage_2 = factor(stage_2,
                          levels = c("z2", "z9", "pl")),
         food = factor(food),
         time = factor(time)) |> 
  drop_na()

data_nmds_edited <-
  data_nmds |> 
  filter(!sample %in% c(66:70, 96:100))

# Experiment 2 - size and time of metamorphosis

data_size_time <- 
  read_excel("lysmata_starvation_dataset.xlsx", sheet = 4,
              range = "A1:F181", na = "NA") |> 
  mutate(stage = factor(stage,
                        levels = c("early_pl", "late_pl", "z9")),
         food = factor(food)) 

# BIOMARKERS -------------------------------------------------------------------

## Non-metric Multidimensional Scaling (nMDS) ----------------------------------

# How do the overall biomarkers vary between larvae that faced different feeding regimes?

### Running the nMDS -----------------------------------------------------------

nmds_biomarkers <-
  metaMDS(comm = data_nmds_edited |> 
                 dplyr::select(total_protein_g:TAC),
          autotransform = FALSE,
          distance = "euclidean",
          k = 2,
          maxit = 300,
          try = 40,
          trymax = 100)

# evaluating the model fitting 

nmds_biomarkers$stress

plot(nmds_biomarkers$diss, nmds_biomarkers$dist)

stressplot(object = nmds_biomarkers, lwd = 5)

### Creating the graphic -------------------------------------------------------

# adding the nMDS coordinates to the dataset

nmds_dataset <- 
  bind_cols(data_nmds_edited,
            as_tibble(nmds_biomarkers$points))

# changing the names of the graph facets

stage_labs <- 
  c("Zoea 2", "Zoea 9", "Post-larvae")

names(stage_labs) <- 
  c("z2", "z9", "pl")

# plotting

graph_nmds_biomarkers_a <-
  ggplot(nmds_dataset,
         aes(x = MDS1, y = MDS2)) +
  geom_point(aes(fill = food, shape = time), size = 4, alpha = 0.8) +
  facet_wrap(~ stage_2, ncol = 1,
             labeller = labeller(stage_2 = stage_labs)) 

graph_nmds_biomarkers_b <-
  graph_nmds_biomarkers_a +
  scale_x_continuous(name = "",
                     limits = c(-210, 110)) +
  scale_y_continuous(name = "",
                     limits = c(-100, 100)) +
  scale_fill_manual(name = "Feeding treatment",
                    values = c("#336699", "#FF0000" )) +
  scale_shape_manual(name = "Time",
                     values = c(21, 22, 24)) 

graph_nmds_biomarkers_c <-
  graph_nmds_biomarkers_b +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    legend.box.spacing = unit(0.05, "cm"),
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.9),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),
    legend.spacing.x = unit(0.25, "cm"),
    legend.key.size = unit(0.25, "cm"),
    legend.background = element_rect(fill = "transparent"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(size = 12, face = "bold"),
    strip.background = element_rect(colour = "black", linewidth = 0.5),
    panel.spacing.y = unit(0.25, "cm")) 

graph_nmds_biomarkers_d <-
  graph_nmds_biomarkers_c +
  geom_text(label = "2D Stress = 0.018", x = -190, y = 95, 
            size = 3, fontface = "bold") +
  guides(shape = guide_legend(
                 ncol =  3,
                 override.aes = list(size = 2.5)),
         fill = guide_legend(
                ncol = 2,
                override.aes = list(size = 2.5, shape = 21)))

# saving the graph

ggsave("Figure 2.png", graph_nmds_biomarkers_d,
       height = 22.5, width = 15, units = "cm", dpi = 600)

## PERMANOVA -------------------------------------------------------------------

# Testing differences in all biomarkers for the Experiment 1 and Experiment 2

### Experiment 1 - Zoea 2 ------------------------------------------------------

# PERMANOVA did not allow NAs so we will need to remove them

# separating the numerical data and the factors in different matrices

data_biomarkers_1 <-
  data_experiment_1_edited |>
  drop_na() |>
  dplyr::select(total_protein_g : TAC)

factors_1 <-
  data_experiment_1_edited |> 
  drop_na() |> 
  dplyr::select(female : time)

# female is a "blocking factor" since larvae hatching from the same female were more similar to each other - we need to add the strata to the model

# running the PERMANOVA

with(factors_1, 
     adonis2(decostand(data_biomarkers_1, method = "normalize") ~ food * time,
             method = "euc", by = "terms",
             data = factors_1, strata = female))

data_experiment_1_pairwise <-
  data_experiment_1_edited |> 
  drop_na() |> 
  as.data.frame()

# pairwise comparison tests

# time 12h

data_pairwise_experiment_1_12 <-
  data_experiment_1_edited |> 
  drop_na() |> 
  filter(time == "12") |> 
  as.data.frame()

data_biomarkers_12 <-
  data_pairwise_experiment_1_12 |> 
  dplyr::select(total_protein_g : TAC)

data_distance_12 <-
  vegdist(decostand(data_biomarkers_12, method = "normalize"), 
          method = "euclidean")

pairwise.adonis2(
  data_distance_12 ~ food, strata = 'female',
  data = data_pairwise_experiment_1_12)

# time 24h

data_pairwise_experiment_1_24 <-
  data_experiment_1_edited |> 
  drop_na() |> 
  filter(time == "24") |> 
  as.data.frame()

data_biomarkers_24 <-
  data_pairwise_experiment_1_24 |> 
  dplyr::select(total_protein_g : TAC)

data_distance_24 <-
  vegdist(decostand(data_biomarkers_24, method = "normalize"), 
          method = "euclidean")

pairwise.adonis2(
  data_distance_24 ~ food, strata = 'female',
  data = data_pairwise_experiment_1_24)

# time 48h

data_pairwise_experiment_1_48 <-
  data_experiment_1_edited |> 
  drop_na() |> 
  filter(time == "48") |> 
  as.data.frame()

data_biomarkers_48 <-
  data_pairwise_experiment_1_48 |> 
  dplyr::select(total_protein_g : TAC)

data_distance_48 <-
  vegdist(decostand(data_biomarkers_48, method = "normalize"), 
          method = "euclidean")

pairwise.adonis2(
  data_distance_48 ~ food, strata = 'female',
  data = data_pairwise_experiment_1_48)

# using SIMPER to calculate what biomarker contribute mostly to differentiate the groups 

sim_food_1 <- 
  with(factors_1, 
       simper(decostand(data_biomarkers_1, method = "normalize"),
              food))

summary(sim_food_1) 

# SOD, GST and total protein were the most important biomarkers discriminating larvae that received food or were starved

sim_time_1 <- 
  with(factors_1, 
       simper(decostand(data_biomarkers_1, method = "normalize"), 
              time))

summary(sim_time_1)

# SOD and GST were the most important biomarkers discriminating larvae across the different feeding times

### Experiment 2 - zoea 9 x early and late post-larvae -------------------------

# separating the numerical data and the factors in different matrices

data_biomarkers_2 <-
  data_experiment_2 |>
  drop_na() |> 
  dplyr::select(total_protein_g : TAC)

factors_2 <-
  data_experiment_2 |> 
  drop_na() |> 
  dplyr::select(stage : female)

# female is a "blocking factor" since larvae from the same female were more similar to each other - we need to add the strata to the model

# running the PERMANOVA

with(factors_2, 
     adonis2(decostand(data_biomarkers_2, method = "normalize") ~ food * stage,
             method = "euc", by = "terms",
             data = factors_2, strata = female))

# pairwise comparison tests

# stage

data_pairwise_experiment_2 <-
  data_experiment_2 |> 
  drop_na() |> 
  as.data.frame()

data_distance_experiment_2 <-
  vegdist(data_biomarkers_2, method = "euclidean")

pairwise.adonis2(
  data_distance_experiment_2 ~ stage, strata = 'female', 
  data = data_pairwise_experiment_2)

# using SIMPER to calculate what biomarker contribute mostly to differentiate the stages

sim_stages_2 <- 
  with(factors_2, 
       simper(decostand(data_biomarkers_2, method = "normalize"), 
              stage))

summary(sim_stages_2) 

# total protein, SOD and GST were the most important biomarkers discriminating early and late post-larvae settlers. CAT and TAC were the most important biomarkers discriminating zoea 9 and early post-larvae settlers

## Linear mixed-effects models -------------------------------------------------

# Testing differences in each biomarker separately for the Experiment 1 and Experiment 2

### Experiment 1 - Zoea 2 ------------------------------------------------------

# total protein 

model_protein_1 <-
  lmer(sqrt(total_protein_g) ~ food * time + (1|female),
       na.action = na.omit,
       data = data_experiment_1_edited)

check_model(model_protein_1, check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_protein_1)

anova(model_protein_1)

pairs(emmeans(model_protein_1, "time"))

pairs(emmeans(model_protein_1, "food"))

data_experiment_1_edited |> 
  group_by(food) |> 
  summarise(mean_tp = mean(total_protein_g))

# GST 

model_gst_1 <-
  lmer(sqrt(GST) ~ food * time + (1|female),
       na.action = na.omit,
       data = data_experiment_1_edited)

check_model(model_gst_1, check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_gst_1)

anova(model_gst_1)

pairs(emmeans(model_gst_1, "time"))

# SOD

model_sod_1 <-
  lmer(sqrt(SOD) ~ food * time + (1|female),
       na.action = na.omit,
       data = data_experiment_1_edited)

check_model(model_sod_1, check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_sod_1)

anova(model_sod_1)

pairs(emmeans(model_sod_1, "food", by = "time"))

pairs(emmeans(model_sod_1, "time", by = "food"))

# CAT

model_cat_1 <-
  lmer(sqrt(CAT + 1) ~ food * time + (1|female),
       na.action = na.omit,
       data = data_experiment_1_edited)

check_model(model_cat_1, check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_cat_1)

anova(model_cat_1)

pairs(emmeans(model_cat_1, "time"))

pairs(emmeans(model_cat_1, "food"))

data_experiment_1_edited |> 
  group_by(food) |> 
  summarise(mean_tp = mean(CAT, na.rm = TRUE))

# LPO

model_lpo_1 <-
  lmer(log(LPO) ~ food * time + (1|female),
       na.action = na.omit,
       data = data_experiment_1_edited)

check_model(model_lpo_1, check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_lpo_1)

anova(model_lpo_1)

pairs(emmeans(model_lpo_1, "food", by = "time"))

pairs(emmeans(model_lpo_1, "time", by = "food"))

# TAC

model_tac_1 <-
  lmer(log(TAC) ~ food * time + (1|female),
       na.action = na.omit,
       data = data_experiment_1_edited)

check_model(model_tac_1, check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_tac_1)

anova(model_tac_1)

pairs(emmeans(model_tac_1, "food", by = "time"))

pairs(emmeans(model_tac_1, "time", by = "food"))

### Experiment 2 - Zoea 9 x early and late post-larvae -------------------------

# total protein 

model_protein_2 <-
  lmer(sqrt(total_protein_g) ~ food * stage + (1|female),
       na.action = na.omit,
       data = data_experiment_2)

check_model(model_protein_2, check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_protein_2)

anova(model_protein_2)

pairs(emmeans(model_protein_2, "stage"))

# GST 

model_gst_2 <-
  lmer(sqrt(GST) ~ food * stage + (1|female),
       na.action = na.omit,
       data = data_experiment_2)

check_model(model_gst_2, check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_gst_2)

anova(model_gst_2)

pairs(emmeans(model_gst_2, "stage"))

# SOD

model_sod_2 <-
  lmer(sqrt(SOD + 7.5) ~ food * stage + (1|female),
       na.action = na.omit,
       data = data_experiment_2)

check_model(model_sod_2, check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_sod_2)

anova(model_sod_2)

pairs(emmeans(model_sod_2, "food", by = "stage"))

# CAT

model_cat_2 <-
  lmer(sqrt(CAT + 1.1) ~ food * stage + (1|female),
       na.action = na.omit,
       data = data_experiment_2)

check_model(model_cat_2, check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_cat_2)

anova(model_cat_2)

pairs(emmeans(model_cat_2, "stage"))

pairs(emmeans(model_cat_2, "food"))

# LPO

model_lpo_2 <-
  lmer(log(LPO) ~ food * stage + (1|female),
       na.action = na.omit,
       data = data_experiment_2)

check_model(model_lpo_2, check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_lpo_2)

anova(model_lpo_2)

pairs(emmeans(model_lpo_2, "stage"))

# TAC

model_tac_2 <-
  lmer(TAC ~ food * stage + (1|female),
       na.action = na.omit,
       data = data_experiment_2)

check_model(model_tac_2, check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_tac_2)

anova(model_tac_2)

pairs(emmeans(model_tac_2, "stage"))

## Graphics --------------------------------------------------------------------

### Experiment 1 - Zoea 2 ------------------------------------------------------

#### TP (Total Protein) --------------------------------------------------------

graph_tp_1_a <-
  ggplot(data_experiment_1_edited,
         aes(x = food, y = total_protein_g, fill = food_time)) +
  geom_violin(width = 0.8, alpha = 0.5, show.legend = FALSE,
              position = position_dodge(0.75)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.75, width = 0.2,
               position = position_dodge(0.75)) +
  labs(x = "**Treatment**", 
       y = "**TP (mg of total protein / g wet weight)**",
       fill = "**Time**")

graph_tp_1_b <-
  graph_tp_1_a +
  scale_x_discrete(labels = c("Fed", "Starved")) +
  scale_y_continuous(breaks = seq(8, 24, by = 4),
                     limits = c(8, 27)) +
  scale_fill_manual(values = c("#9999FF", "#0066FF", "#000099",
                               "#FF9999", "#FF0000", "#990000"),
                    labels = c("F12", "F24", "F48", "S12", "S24", "S48"))

graph_tp_1_c <-
  graph_tp_1_b +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 10, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.x = element_markdown(size = 12, 
                       margin = margin(t = 10, r = 20, b = 10, l = 20)),
        axis.title.y = element_markdown(size = 12, 
                       margin = margin(t = 10, r = 10, b = 10, l = 10)),
        legend.position = "inside",
        legend.position.inside = c(0.45, 0.89),
        legend.title = element_markdown(size = 10),
        legend.text = element_text(size = 8),
        legend.spacing.x = unit(0.2, "cm"),
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_rect(fill = "transparent"),
        panel.border = element_rect(linewidth = 1, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(
    byrow = TRUE, nrow = 1))

#### GST (Glutationa-S-Transferase) --------------------------------------------

graph_gst_1_a <-
  ggplot(data_experiment_1_edited,
         aes(x = food, y = GST, fill = food_time)) +
  geom_violin(width = 0.8, alpha = 0.5, show.legend = FALSE,
              position = position_dodge(0.75)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.75, width = 0.2,
               position = position_dodge(0.75)) +
  labs(x = "**Treatment**", 
       y = "**GST (nmol / min / mg of TP)**",
       fill = "**Time**")

graph_gst_1_b <-
  graph_gst_1_a +
  scale_x_discrete(labels = c("Fed", "Starved")) +
  scale_y_continuous(breaks = seq(50, 250, by = 50),
                     limits = c(50, 250)) +
  scale_fill_manual(values = c("#9999FF", "#0066FF", "#000099",
                               "#FF9999", "#FF0000", "#990000"),
                    labels = rep(c("12", "24", "48"), 2))

graph_gst_1_c <-
  graph_gst_1_b +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 10, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.x = element_markdown(size = 12, 
                       margin = margin(t = 10, r = 20, b = 10, l = 20)),
        axis.title.y = element_markdown(size = 12, 
                      margin = margin(t = 10, r = 10, b = 10, l = 10)),
        legend.position = "none",
        panel.border = element_rect(linewidth = 1, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#### SOD (Superoxide Dismutase) ------------------------------------------------

graph_sod_1_a <-
  ggplot(data_experiment_1_edited,
         aes(x = food, y = SOD, fill = food_time)) +
  geom_violin(width = 0.8, alpha = 0.5, show.legend = FALSE,
              position = position_dodge(0.75)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.75, width = 0.2,
               position = position_dodge(0.75)) +
  labs(x = "**Treatment**", 
       y = "**SOD (% inhibition / mg of TP)**",
       fill = "**Time**")

graph_sod_1_b <-
  graph_sod_1_a +
  scale_x_discrete(labels = c("Fed", "Starved")) +
  scale_y_continuous(breaks = seq(50, 250, by = 50),
                     limits = c(50, 250)) +
  scale_fill_manual(values = c("#9999FF", "#0066FF", "#000099",
                               "#FF9999", "#FF0000", "#990000"),
                    labels = rep(c("12", "24", "48"), 2))

graph_sod_1_c <-
  graph_sod_1_b +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 10, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.x = element_markdown(size = 12, 
                       margin = margin(t = 10, r = 20, b = 10, l = 20)),
        axis.title.y = element_markdown(size = 12, 
                       margin = margin(t = 10, r = 10, b = 10, l = 10)),
        legend.position = "none",
        panel.border = element_rect(linewidth = 1, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#### CAT (Catalase) ------------------------------------------------------------

graph_cat_1_a <-
  ggplot(data_experiment_1_edited |> 
         drop_na(CAT),
         aes(x = food, y = CAT, fill = food_time)) +
  geom_violin(width = 0.8, alpha = 0.5, show.legend = FALSE,
              position = position_dodge(0.75)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.75, width = 0.2,
               position = position_dodge(0.75)) +
  labs(x = "**Treatment**", 
       y = expression(bold("CAT"~("\u03bc"*M~"/"~"min"~"/"~"mg of TP"))),
       fill = "**Time**")

graph_cat_1_b <-
  graph_cat_1_a +
  scale_x_discrete(labels = c("Fed", "Starved")) +
  scale_y_continuous(breaks = seq(0, 10, by = 2.5),
                     limits = c(-0.25, 10.25)) +
  scale_fill_manual(values = c("#9999FF", "#0066FF", "#000099",
                               "#FF9999", "#FF0000", "#990000"),
                    labels = rep(c("12", "24", "48"), 2))

graph_cat_1_c <-
  graph_cat_1_b +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 10, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.x = element_markdown(size = 12, 
                       margin = margin(t = 10, r = 20, b = 10, l = 20)),
        axis.title.y = element_text(size = 12,
                       margin = margin(t = 10, r = 10, b = 10, l = 10)),
        legend.position = "none",
        panel.border = element_rect(linewidth = 1, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#### LPO (Lipid Peroxidation) --------------------------------------------------

graph_lpo_1_a <-
  ggplot(data_experiment_1_edited |> 
         drop_na(LPO),
         aes(x = food, y = LPO, fill = food_time)) +
  geom_violin(width = 0.8, alpha = 0.5, show.legend = FALSE,
              position = position_dodge(0.75)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.75, width = 0.2,
               position = position_dodge(0.75)) +
  labs(x = "**Treatment**", 
       y = "**LPO (nmol / mg of TP)**",
       fill = "**Time**")

graph_lpo_1_b <-
  graph_lpo_1_a +
  scale_x_discrete(labels = c("Fed", "Starved")) +
  scale_y_continuous(breaks = seq(0, 0.08, by = 0.02),
                     limits = c(0, 0.08)) +
  scale_fill_manual(values = c("#9999FF", "#0066FF", "#000099",
                               "#FF9999", "#FF0000", "#990000"),
                    labels = rep(c("12", "24", "48"), 2))

graph_lpo_1_c <-
  graph_lpo_1_b +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 10, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.x = element_markdown(size = 12, 
                       margin = margin(t = 10, r = 20, b = 10, l = 20)),
        axis.title.y = element_markdown(size = 12,
                       margin = margin(t = 10, r = 10, b = 10, l = 10)),
        legend.position = "none",
        panel.border = element_rect(linewidth = 1, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#### TAC (Total Antioxidant Capacity) ------------------------------------------

graph_tac_1_a <-
  ggplot(data_experiment_1_edited |> 
         drop_na(TAC),
         aes(x = food, y = TAC, fill = food_time)) +
  geom_violin(width = 0.8, alpha = 0.5, show.legend = FALSE,
              position = position_dodge(0.75)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.75, width = 0.2,
               position = position_dodge(0.75)) +
  labs(x = "**Treatment**", 
       y = "**TAC (mmol / mg of TP)**",
       fill = "**Time**")

graph_tac_1_b <-
  graph_tac_1_a +
  scale_x_discrete(labels = c("Fed", "Starved")) +
  scale_y_continuous(breaks = seq(0, 0.8, by = 0.2),
                     limits = c(0, 0.8)) +
  scale_fill_manual(values = c("#9999FF", "#0066FF", "#000099",
                               "#FF9999", "#FF0000", "#990000"),
                    labels = rep(c("12", "24", "48"), 2))

graph_tac_1_c <-
  graph_tac_1_b +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 10, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.x = element_markdown(size = 12, 
                       margin = margin(t = 10, r = 20, b = 10, l = 20)),
        axis.title.y = element_markdown(size = 12,
                       margin = margin(t = 10, r = 10, b = 10, l = 10)),
        legend.position = "none",
        panel.border = element_rect(linewidth = 1, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

### Creating the final figure --------------------------------------------------

graph_biomarkers_experiment_1 <-
  graph_tp_1_c +
  graph_gst_1_c +
  graph_cat_1_c +
  graph_sod_1_c +
  graph_lpo_1_c +
  graph_tac_1_c +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

# remove x title and text for plot 1

graph_biomarkers_experiment_1[[1]] <-
  graph_biomarkers_experiment_1[[1]] +
  theme(
    plot.margin = margin(0.5, 0.25, 0.25, 0.5, "cm"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank())

# remove x title and text for plot 2

graph_biomarkers_experiment_1[[2]] <-
  graph_biomarkers_experiment_1[[2]] +
  theme(
    plot.margin = margin(0.5, 0.5, 0, 0, "cm"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank())

# remove x title and text for plot 3

graph_biomarkers_experiment_1[[3]] <-
  graph_biomarkers_experiment_1[[3]] +
  theme(
    plot.margin = margin(0, 0.25, 0.25, 0.5, "cm"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank())

# remove x title and text for plot 4

graph_biomarkers_experiment_1[[4]] <-
  graph_biomarkers_experiment_1[[4]] +
  theme(
    plot.margin = margin(0, 0.5, 0, 0, "cm"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank())

# change margins for plots 5

graph_biomarkers_experiment_1[[5]] <-
  graph_biomarkers_experiment_1[[5]] +
  theme(plot.margin = margin(0, 0.25, 0.25, 0.5, "cm"))

# change margins for plot 6

graph_biomarkers_experiment_1[[6]] <-
  graph_biomarkers_experiment_1[[6]] +
  theme(
    plot.margin = margin(0, 0.5, 0.5, 0, "cm"))

ggsave("Figure 3.png",
       graph_biomarkers_experiment_1,
       height = 25, width = 25, units = "cm", dpi = 600)

### Experiment 2 - Zoea 9 and Postlarvae ---------------------------------------

#### TP (Total Protein) ---------------------------------------------------------

graph_tp_2_a <-
  ggplot(data_experiment_2 |> 
         drop_na(total_protein_g),
         aes(x = stage, y = total_protein_g, fill = food)) +
  geom_violin(width = 0.8, alpha = 0.5, show.legend = FALSE,
              position = position_dodge(0.75)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.75, width = 0.2,
               position = position_dodge(0.75)) +
  labs(x = "**Stage**", 
       y = "**TP (mg of total protein / g wet weight)**",
       fill = "**Treatment**")

graph_tp_2_b <-
  graph_tp_2_a +
  scale_x_discrete(labels = c("Early PL", "Late PL", "Z9")) +
  scale_y_continuous(breaks = seq(15, 55, by = 10),
                     limits = c(15, 56)) +
  scale_fill_manual(labels = c("Fed", "Starved"),
                    values = c("#336699", "#FF0000"))

graph_tp_2_c <-
  graph_tp_2_b +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 10, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.x = element_markdown(size = 12, 
                       margin = margin(t = 10, r = 20, b = 10, l = 20)),
        axis.title.y = element_markdown(size = 12, 
                       margin = margin(t = 10, r = 10, b = 10, l = 10)),
        legend.position = "inside",
        legend.position.inside = c(0.22, 0.89),
        legend.title = element_markdown(size = 10),
        legend.text = element_text(size = 8),
        legend.spacing.x = unit(0.25, "cm"),
        legend.key.size = unit(0.5, "cm"),
        legend.background = element_rect(fill = "transparent"),
        panel.border = element_rect(linewidth = 1, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(
    byrow = TRUE, nrow = 1))

#### GST (Glutationa-S-Transferase) --------------------------------------------

graph_gst_2_a <-
  ggplot(data_experiment_2,
         aes(x = stage, y = GST, fill = food)) +
  geom_violin(width = 0.8, alpha = 0.5, show.legend = FALSE,
              position = position_dodge(0.75)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.75, width = 0.2,
               position = position_dodge(0.75)) +
  labs(x = "**Stage**", 
       y = "**GST (nmol / min / mg of TP)**",
       fill = "**Treatment**")

graph_gst_2_b <-
  graph_gst_2_a +
  scale_x_discrete(labels = c("Early PL", "Late PL", "Z9")) +
  scale_y_continuous(breaks = seq(40, 200, by = 40),
                     limits = c(39, 200)) +
  scale_fill_manual(labels = c("Fed", "Starved"),
                    values = c("#336699", "#FF0000"))

graph_gst_2_c <-
  graph_gst_2_b +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 10, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.x = element_markdown(size = 12, 
                       margin = margin(t = 10, r = 20, b = 10, l = 20)),
        axis.title.y = element_markdown(size = 12, 
                       margin = margin(t = 10, r = 10, b = 10, l = 10)),
        legend.position = "none",
        panel.border = element_rect(linewidth = 1, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#### SOD (Superoxide Dismutase) ------------------------------------------------

graph_sod_2_a <-
  ggplot(data_experiment_2 |> 
         drop_na(SOD),
         aes(x = stage, y = SOD, fill = food)) +
  geom_violin(width = 0.8, alpha = 0.5, show.legend = FALSE,
              position = position_dodge(0.75)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.75, width = 0.2,
               position = position_dodge(0.75)) +
  labs(x = "**Stage**", 
       y = "**SOD (% inhibition / mg of TP)**",
       fill = "**Treatment**")

graph_sod_2_b <-
  graph_sod_2_a +
  scale_x_discrete(labels = c("Early PL", "Late PL", "Z9")) +
  scale_y_continuous(breaks = seq(-10, 90, by = 25),
                     limits = c(-10, 90)) +
  scale_fill_manual(labels = c("Fed", "Starved"),
                    values = c("#336699", "#FF0000"))

graph_sod_2_c <-
  graph_sod_2_b +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 10, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.x = element_markdown(size = 12, 
                       margin = margin(t = 10, r = 20, b = 10, l = 20)),
        axis.title.y = element_markdown(size = 12, 
                       margin = margin(t = 10, r = 10, b = 10, l = 10)),
        legend.position = "none",
        panel.border = element_rect(linewidth = 1, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#### CAT (Catalase) ------------------------------------------------------------

graph_cat_2_a <-
  ggplot(data_experiment_2 |> 
         drop_na(CAT),
         aes(x = stage, y = CAT, fill = food)) +
  geom_violin(width = 0.8, alpha = 0.5, show.legend = FALSE,
              position = position_dodge(0.75)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.75, width = 0.2,
               position = position_dodge(0.75)) +
  labs(x = "**Stage**", 
       y = expression(bold("CAT"~("\u03bc"*M~"/"~"min"~"/"~"mg of TP"))),
       fill = "**Treatment**")

graph_cat_2_b <-
  graph_cat_2_a +
  scale_x_discrete(labels = c("Early PL", "Late PL", "Z9")) +
  scale_y_continuous(breaks = seq(-2, 6, by = 2),
                     limits = c(-2, 6)) +
  scale_fill_manual(labels = c("Fed", "Starved"),
                    values = c("#336699", "#FF0000"))

graph_cat_2_c <-
  graph_cat_2_b +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 10, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.x = element_markdown(size = 12, 
                       margin = margin(t = 10, r = 20, b = 10, l = 20)),
        axis.title.y = element_text(size = 12,
                       margin = margin(t = 10, r = 10, b = 10, l = 10)),
        legend.position = "none",
        panel.border = element_rect(linewidth = 1, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#### LPO (Lipid Peroxidation) --------------------------------------------------

graph_lpo_2_a <-
  ggplot(data_experiment_2,
         aes(x = stage, y = LPO, fill = food)) +
  geom_violin(width = 0.8, alpha = 0.5, show.legend = FALSE,
              position = position_dodge(0.75)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.75, width = 0.2,
               position = position_dodge(0.75)) +
  labs(x = "**Stage**", 
       y = "**LPO (nmol / mg of TP)**",
       fill = "**Treatment**")

graph_lpo_2_b <-
  graph_lpo_2_a +
  scale_x_discrete(labels = c("Early PL", "Late PL", "Z9")) +
  scale_y_continuous(breaks = seq(0, 0.06, by = 0.015),
                     limits = c(0, 0.06)) +
  scale_fill_manual(labels = c("Fed", "Starved"),
                    values = c("#336699", "#FF0000"))

graph_lpo_2_c <-
  graph_lpo_2_b +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 10, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.x = element_markdown(size = 12, 
                       margin = margin(t = 10, r = 20, b = 10, l = 20)),
        axis.title.y = element_markdown(size = 12,
                       margin = margin(t = 10, r = 10, b = 10, l = 10)),
        legend.position = "none",
        panel.border = element_rect(linewidth = 1, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#### TAC (Total Antioxidant Capacity) ------------------------------------------

graph_tac_2_a <-
  ggplot(data_experiment_2,
         aes(x = stage, y = TAC, fill = food)) +
  geom_violin(width = 0.8, alpha = 0.5, show.legend = FALSE,
              position = position_dodge(0.75)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.75, width = 0.2,
               position = position_dodge(0.75)) +
  labs(x = "**Stage**", 
       y = "**TAC (mmol / mg of TP)**",
       fill = "**Treatment**")

graph_tac_2_b <-
  graph_tac_2_a +
  scale_x_discrete(labels = c("Early PL", "Late PL", "Z9")) +
  scale_y_continuous(breaks = seq(0, 0.8, by = 0.2),
                     limits = c(-0.05, 0.85)) +
  scale_fill_manual(labels = c("Fed", "Starved"),
                    values = c("#336699", "#FF0000"))

graph_tac_2_c <-
  graph_tac_2_b +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 10, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.x = element_markdown(size = 12, 
                       margin = margin(t = 10, r = 20, b = 10, l = 20)),
        axis.title.y = element_markdown(size = 12,
                       margin = margin(t = 10, r = 10, b = 10, l = 10)),
        legend.position = "none",
        panel.border = element_rect(linewidth = 1, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#### Creating the final figure --------------------------------------------------

graph_biomarkers_experiment_2 <-
  graph_tp_2_c +
  graph_gst_2_c +
  graph_cat_2_c +
  graph_sod_2_c +
  graph_lpo_2_c +
  graph_tac_2_c +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

# remove x title and text for plot 1

graph_biomarkers_experiment_2[[1]] <-
  graph_biomarkers_experiment_2[[1]] +
  theme(
    plot.margin = margin(0.5, 0.25, 0.25, 0.5, "cm"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank())

# remove x title and text for plot 2

graph_biomarkers_experiment_2[[2]] <-
  graph_biomarkers_experiment_2[[2]] +
  theme(
    plot.margin = margin(0.5, 0.5, 0, 0, "cm"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank())

# remove x title and text for plot 3

graph_biomarkers_experiment_2[[3]] <-
  graph_biomarkers_experiment_2[[3]] +
  theme(
    plot.margin = margin(0, 0.25, 0.25, 0.5, "cm"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank())

# remove x title and text for plot 4

graph_biomarkers_experiment_2[[4]] <-
  graph_biomarkers_experiment_2[[4]] +
  theme(
    plot.margin = margin(0, 0.5, 0, 0, "cm"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank())

# change margins for plots 5

graph_biomarkers_experiment_2[[5]] <-
  graph_biomarkers_experiment_2[[5]] +
  theme(plot.margin = margin(0, 0.25, 0.25, 0.5, "cm"))

# change margins for plot 6

graph_biomarkers_experiment_2[[6]] <-
  graph_biomarkers_experiment_2[[6]] +
  theme(
    plot.margin = margin(0, 0.5, 0.5, 0, "cm"))

ggsave("Figure 4.png",
       graph_biomarkers_experiment_2,
       height = 25, width = 25, units = "cm", dpi = 600)

# SIZE AND TIME OF METAMORPHOSIS -----------------------------------------------

## Linear models ---------------------------------------------------------------

# Let's test whether the size of the larvae and the time of metamorphosis differed between the larval stages (early and late post-larvae and zoea 9) and the feeding treatments

model_size <-
  lmer(size ~ stage * food + (1|female),
       data = data_size_time)

check_model(model_size, check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_size)

anova(model_size)

pairs(emmeans(model_size, "stage"))

pairs(emmeans(model_size, "food"))

pairs(emmeans(model_size, "food", by = "stage"))

data_size_time |> 
  group_by(stage) |> 
  summarise(mean_size = mean(size))

data_size_time |> 
  group_by(food) |> 
  summarise(mean_size = mean(size))

model_time <-
  lmer(time ~ stage * food + (1|female),
       data = data_size_time)

check_model(model_time, check = c("qq", "normality", "homogeneity"))

check_heteroscedasticity(model_time)

anova(model_time)

pairs(emmeans(model_time, "stage"))

pairs(emmeans(model_time, "food"))

data_size_time |> 
  group_by(stage) |> 
  summarise(mean_time = mean(time))

data_size_time |> 
  filter(!stage == "z9") |> 
  group_by(food) |> 
  summarise(mean_time = mean(time))

## Graphic ---------------------------------------------------------------------

data_size_time_summary <-
  data_size_time |> 
  group_by(stage, food) |> 
  summarise(sample_size = n(),
            mean_time = mean(time),
            mean_size = mean(size),
            sd_time = sd(time),
            sd_size = sd(size)) |> 
  mutate(se_time = sd_time / sqrt(sample_size),
         se_size = sd_size / sqrt(sample_size))

graph_size_time_a <-
  ggplot(data = data_size_time_summary,
         aes(x = mean_time, y = mean_size, 
             fill = food, shape = stage)) + 
  geom_errorbar(aes(ymin = mean_size - sd_size,
                    ymax = mean_size + sd_size), 
                width = 0.5,
                position = position_dodge(0)) + 
  geom_errorbarh(aes(xmin = mean_time - sd_time,
                     xmax = mean_time + sd_time),
                 height = 0.2) + 
  geom_point(size = 5) + 
  labs(x = "**Time (days)**", 
       y = "**Size (mm)**",
       fill = "**Treatment**",
       shape = "**Stage**")

graph_size_time_b <-
  graph_size_time_a +
  scale_x_continuous(breaks = seq(25, 45, by = 5),
                     limits = c(25, 45.2)) +
  scale_y_continuous(breaks = seq(6, 10, by = 1),
                     limits = c(5.8, 10)) +
  scale_fill_manual(labels = c("Fed", "Starved"),
                    values = c("#336699", "#FF0000")) +
  scale_shape_manual(labels = c("Early PL", "Late PL", "Z9"),
                     values = c(21, 22, 24))

graph_size_time_c <-
  graph_size_time_b +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 10, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.x = element_markdown(size = 12, 
                       margin = margin(t = 10, r = 20, b = 10, l = 20)),
        axis.title.y = element_markdown(size = 12, 
                       margin = margin(t = 10, r = 10, b = 10, l = 10)),
        legend.position = "inside",
        legend.position.inside = c(0.2, 0.85),
        legend.title = element_markdown(size = 10),
        legend.text = element_text(size = 8),
        legend.spacing.x = unit(0.25, "cm"),
        legend.key.size = unit(0.5, "cm"),
        legend.background = element_rect(fill = "transparent"),
        panel.border = element_rect(linewidth = 1, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(byrow = TRUE, nrow = 1, 
                             override.aes = list(shape = 21)),
         shape = guide_legend(byrow = TRUE, nrow = 1))

ggsave("Figure 1.png",
       graph_size_time_c,
       height = 15, width = 20, units = "cm", dpi = 600)
