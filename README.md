# Codes and data of the manuscript "Post-hatching starvation increases oxidative stress responses and prolongs larval life in the marine ornamental shrimp *Lysmata seticaudata*" submitted to Aquaculture.

**Authors**: Rafael Campos Duarte; Daniel Acebes; Diana Madeira; Ruben X. G. Silva; Joana F. Fernandes; Fernando Ricardo; Ricardo Calado

**Contact**: Rafael Campos Duarte, rafaduarte87@gmail.com / rafael.duarte@ua.pt

**Abstract**: We experimentally assessed the effects of different feeding regimes on newly hatched larvae of the shrimp *Lysmata seticaudata* using a range of biomarkers of oxidative stress and antioxidant response. We showed that starving newly hatched larvae promoted short-term stress responses at the subsequent larval stage, particularly by increasing the activity of oxidative stress enzymes and the total antioxidant capacity of larvae. Such a brief starvation period at the beginning of larval development had significant long-term effects, delaying the metamorphosis to post-larvae and resulting in the metamorphosis of larger shrimp. 

**Responsible for running the experiments and collecting data**: Daniel Acebes

**Responsible for writing the codes**: Rafael Campos Duarte

## List of files
* **dataset**: lysmata_starvation_dataset.xlsx
* **script**: script_lysmata_starvation.R

The script contains the code for the analyses and figures of the different activities conducted in the study. 

## Session info
* R version 4.4.2 (2024-10-31 ucrt) -- Platform: x86_64-w64-mingw32/x64 (64-bit)

## Packages installed and their versions in R 4.4.2
* emmeans v.1.10.6
* ggtext v.0.1.2
* lme4 v.1.1-36
* lmerTest v.3.1-3
* pairwiseAdonis v.0.4.1
* patchwork v.1.3.0
* performance v.0.13.0
* readxl v.1.4.3
* tidyverse v.2.0.0
* vegan v.2.6-8

**Workflow**: To replicate the figures and statistical analyses from the manuscript, run the codes available on the file script_lysmata_starvation.R

**General information about the scripts and datasets**: The script and datasets refer to the analyses of: (i) **Experiment 1** - the total protein content and the activity of five oxidative stress biomarkers of zoea 2 larvae of the shrimp *L. seticaudata* which were either fed or starved during the first 12, 24 or 48 hours after hatching; (ii) **Experiment 2** - the total protein content and the activity of five oxidative stress biomarkers of *L. seticaudata* post-larvae that metamorphosed in the first (as early post-larvae) or last (as late post-larvae) 10 days after the first metamorphosed individual, as well as zoea 9 larvae that did not undergo metamorphosis, all which were either fed or starved during the first 48 hours after hatching as zoea 1 larvae; (iii) **Experiment 2** - the time of metamorphosis and the size of both *L. seticaudata* post-larvae and zoea 9. The repository consists of two files: one code script and one Excel file with four different spreadsheets.

## Description of the dataset and script

### DATASET

* #### lysmata_starvation_dataset.xlsx 

This dataset contains four spreadsheets. The first spreadsheet shows the total protein content and the activity of five oxidative stress biomarkers of zoea 2 larvae that were either fed or starved over different periods after hatching (Experiment 1). The second spreadsheet shows the total protein content and the activity of five oxidative stress biomarkers of early and late post-larvae as well as zoea 9 that did not undergo a metamorphosis, all of which were either fed or starved during the first 48 hours after hatching (Experiment 2). The third spreadsheet shows the joined data of Experiments 1 and 2 for multivariate analyses. The fourth spreadsheet shows the time of metamorphosis and the size of both post-larvae and zoea 9 (Experiment 2).

**Column names and description**:

*experiment_1*

* sample: the sample identification 
* female: the identification of the shrimp female from which the larvae hatched (1 to 5)
* food: the feeding treatment of zoea 1 larvae (fed - F or starved - S)
* time: the duration of the feeding treatment (12, 24 or 48 hours)
* tank: the identification of the water bath tank where the 24-well plates were deployed
* total_protein_ml: the amount of total protein (TP) per sample volume (mg/ml)
* total_protein_microL: the amount of total protein (TP) after homogenization in buffer solution
* total_protein_g: the amount of total protein (TP) corrected by the sample weight (mg/g)
* GST: glutathione-S-transferase activity (nmol/minute/mg of TP)
* SOD: superoxide dismutase inhibition (% of inhibition/mg of TP)
* CAT: catalase activity (microM/minute/mg of TP)
* LPO: lipid peroxidation (nmol/mg of TP)
* TAC: total antioxidant capacity (nmol/mg of TP)

*experiment_2*

* sample: the sample identification
* stage: the shrimp developmental stage (zoea 9 - z9, early post-larvae - early_pl or late post-larvae - late_pl)
* food: the feeding treatment during the zoea 1 stage (fed - F or starved - S)
* tank: the identification of the larval cultivation tank
* female: the identification of the shrimp female from which the larvae hatched (1 to 3)
* total_protein_ml: the amount of total protein (TP) per sample volume (mg/ml)
* total_protein_microL: the amount of total protein (TP) after homogenization in buffer solution
* total_protein_g: the amount of total protein (TP) corrected by the sample weight (mg/g)
* GST: glutathione-S-transferase activity (nmol/minute/mg of TP)
* SOD: superoxide dismutase inhibition (% of inhibition/mg of TP)
* CAT: catalase activity (microM/minute/mg of TP)
* LPO: lipid peroxidation (nmol/mg of TP)
* TAC: total antioxidant capacity (nmol/mg of TP)

*multivariate*

* sample: the sample identification
* stage_1: the shrimp developmental stage (zoea 2 - z2, zoea 9 - z9, early post-larvae - early_pl or late post-larvae - late_pl)
* stage_2: the shrimp developmental stage (zoea 2 - z2, zoea 9 - z9, post-larvae - pl)
* food: the feeding treatment of zoea 1 larvae (fed - F or starved - S)
* time: the duration of the feeding treatment (12, 24 or 48 hours)
* total_protein_ml: the amount of total protein (TP) per sample volume (mg/ml)
* total_protein_microL: the amount of total protein (TP) after homogenization in buffer solution
* total_protein_g: the amount of total protein (TP) corrected by the sample weight (mg/g)
* GST: glutathione-S-transferase activity (nmol/minute/mg of TP)
* SOD: superoxide dismutase inhibition (% of inhibition/mg of TP)
* CAT: catalase activity (microM/minute/mg of TP)
* LPO: lipid peroxidation (nmol/mg of TP)
* TAC: total antioxidant capacity (nmol/mg of TP)

*time_size*

* stage: the shrimp developmental stage (zoea 9 - z9, early post-larvae - early_pl or late post-larvae - late_pl)
* food: the feeding treatment during the zoea 1 stage (fed - F or starved - S)
* tank: the identification of the larval cultivation tank
* female: the identification of the shrimp female from which the larvae hatched (1 to 3)
* size: the size of the larva or post-larva (in mm)
* time: the number of days after hatching for the metamorphosis to post-larvae or for the maintenance as zoea 9

### SCRIPTS

* #### script_lysmata_starvation.R

Contains the statistical analyses and the codes to create: 

 * **Figure 1**: variation in size (in mm) and time (in days) of metamorphosis of the shrimp *Lysmata seticaudata* according to feeding treatment (larvae fed or starved in the first 48 hours after hatching) and larval duration of post-larvae (classified as early and late settlers), in comparison to larvae that did not undergo metamorphosis, remaining as zoea 9 at the end of the experiment.

 * **Figure 2**: non-metric multidimensional scaling (nMDS) plot illustrating the overall variation in total protein content and in the profile of five different oxidative stress biomarkers throughout the larval development of the shrimp *Lysmata seticaudata* according to different feeding treatments during zoea 1 stage (larvae fed – F or starved – S over different time periods).
   
 * **Figure 3**: levels of total protein content and profiles of five different oxidative stress biomarkers measured on the zoea 2 larval stage of the shrimp *Lysmata seticaudata*, which were fed or starved during 12, 24 or 48 hours after hatching (Experiment 1).

 * **Figure 4**: levels of total protein content and profiles of five different oxidative stress biomarkers measured on early and late settling post-larvae and zoea 9 that did not undergo metamorphosis of the shrimp *Lysmata seticaudata*, which were fed or starved during 48 hours after hatching (Experiment 2).

## Funding Sources

* Fundação de Ciência e Tecnologia (FCT)
