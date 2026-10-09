library(ggplot2)
library(gridExtra)
library(dplyr)
library(fuzzyjoin)
library(export)
library(ggsignif)
library(data.table)
library(patchwork)
library(tableone)
library(tidyr)
library(multcomp)
library(mgcv)
library(mclust)

source("./corrections_ABCDS.R")
source("./code_to_plot.R")
source("./Fig1_otherfuncs.R")
source("./Revision_functions.R")
source("./slidingWindowFuncs.R")
source("./Fig1_updatedFunction.R")

################################################################################
##ABCDS
################################################################################

latency <- read.csv("../ABCDS_LONI_20240702/Age_at_Event_and_Latency_11Jul2024.csv")
demogs <- read.csv("../ABCDS_LONI_20250528/Participant_Demographics_28May2025.csv")
demogs <- demogs[, c("subject_label", "event_sequence", "de_gender", "ds_vs_control_flag")]
demogs <- demogs[!duplicated(demogs),]

latency <- latency %>%
  group_by(subject_label) %>%
  mutate(
    baseline_age = min(age_at_visit, na.rm = TRUE),
    is_min_row = age_at_visit == baseline_age,
    
    clinical_AgefromBaseline = if_else(
      is_min_row, 
      baseline_age,
      baseline_age + (clinical_latency_in_days / 365.25)
    ),
    amy_AgefromBaseline = if_else(
      is_min_row, 
      baseline_age,
      baseline_age + (amy_latency_in_days / 365.25)
    ),
    tau_AgefromBaseline = if_else(
      is_min_row, 
      baseline_age,
      baseline_age + (tau_latency_in_days / 365.25)
    ),
    fdg_AgefromBaseline = if_else(
      is_min_row, 
      baseline_age,
      baseline_age + (fdg_latency_in_days / 365.25)
    ),
    mri_AgefromBaseline = if_else(
      is_min_row, 
      baseline_age,
      baseline_age + (mri_latency_in_days / 365.25)
    ),
    csf_AgefromBaseline = if_else(
      is_min_row, 
      baseline_age,
      baseline_age + (csf_latency_in_days / 365.25)
    )
  ) %>%
  ungroup() %>%
  dplyr::select(-is_min_row)

mcrt <- read.csv(".././ABCDS_LONI_20250528/Cued_Recall_28May2025.csv") #Krasny says 28 for mci and 23 for dementia

ABCDS_corr <- read.csv(".././ABCDS_LONI_20250528/freeze_5_corrections_2025.07.31.csv")


corr_sel <- ABCDS_corr %>%
  dplyr::select(
    SUBJECT_LABEL_BAD,
    subject_label_new = subject_label,
    event_sequence_new = event_sequence,
    CaseControl
  )
corr_sel$CaseControl <- as.factor(corr_sel$CaseControl)
levels(corr_sel$CaseControl) <- c("DS", "Control")

demogs$ds_vs_control_flag <- as.factor(demogs$ds_vs_control_flag)
demogs <- implement_corrections(demogs,
                                subject_label_col = "subject_label",
                                ds_flag_col = "ds_vs_control_flag",
                                corr_sel = corr_sel)

ABCDS_apoe <- read.csv(".././ABCDS_LONI_20250528/ApoE_Genotyping_Results_28May2025.csv")
ABCDS_apoe$APOE_grouped <- ifelse(ABCDS_apoe$allele_combo %in% c("E2/E4", "E4/E2", "E3/E4", "E4/E4", "E4/E3"), "APOE4",
                                  ifelse(ABCDS_apoe$allele_combo %in% c("E3/E3", "E3/3E"), "APOE3", 
                                         ifelse(ABCDS_apoe$allele_combo %in% c("E2/E2", "E2/E3", "E3/E2"), "APOE2", NA)))
ABCDS <- read.csv(".././ABCDS_LONI_20250528/abcds_biospecimen_results_28May2025.csv")
dx <- read.csv(".././ABCDS_LONI_20250528/Consensus_28May2025.csv")
dx <- dx[, c("subject_label", "event_sequence", "consensus_dx")]
dx <- merge(dx, ABCDS_apoe[, c("subject_label", "APOE_grouped", "allele_combo")], by = "subject_label", all = FALSE)
dx <- merge(dx, demogs, by = "subject_label", all.x = TRUE, all.y = FALSE)
GFAP <- ABCDS[ABCDS$TESTNAME %in% "GFAP",]
GFAP <- GFAP[!is.na(GFAP$SUBJECT_LABEL),]

NFL <- ABCDS[ABCDS$TESTNAME == "NF-light",]
ABCDS <- ABCDS[ABCDS$TESTNAME %in% c("P-tau217"),]
ABCDS <- ABCDS[ABCDS$UNITS %in% c("Normalised Mean (pg/ml)"),]
ABCDS <- ABCDS[!is.na(ABCDS$SUBJECT_LABEL),]
colnames(ABCDS)[6] <- c( "pTau217")
ABCDS$pTau217 <- as.numeric(ABCDS$pTau217)
ABCDS <- merge(ABCDS, demogs, by.x = "SUBJECT_LABEL", by.y = "subject_label", all.x = TRUE, all.y = FALSE)
ABCDS <- merge(ABCDS, ABCDS_apoe[, c("subject_label", "APOE_grouped")], by.x = "SUBJECT_LABEL", by.y = "subject_label", all = FALSE)
ABCDS <- merge(ABCDS, mcrt[, c("subject_label", "event_sequence", "trs")], by.x = c("SUBJECT_LABEL", "EVENT_SEQUENCE"),
               by.y = c("subject_label", "event_sequence"), all.x = TRUE, all.y = FALSE)
pTau217_thresh_ABCDS <- mean(ABCDS[ABCDS$age_at_visit <= 35 & !duplicated(ABCDS$SUBJECT_LABEL),]$pTau217, na.rm = TRUE) + 
  1.96 *sd(ABCDS[ABCDS$age_at_visit <= 35 & !duplicated(ABCDS$SUBJECT_LABEL),]$pTau217, na.rm = TRUE) 

ABCDS_CL <- read.csv(".././ABCDS_DF3_AmyloidOnly/all_amy_120924.csv")
ABCDS_CL <- ABCDS_CL[ABCDS_CL$PET.TC.QC.Status == "Passed", c("Subject", "VISIT", "WUSTLcentiloid")]
ABCDS_CL$event_sequence <- as.numeric(substr(ABCDS_CL$VISIT, start = nchar(ABCDS_CL$VISIT), stop = nchar(ABCDS_CL$VISIT)))
ABCDS_CL <- merge(ABCDS_CL, demogs[!duplicated(demogs$subject_label), 
                                   c("subject_label", "de_gender", "ds_vs_control_flag")], by.x = "Subject", by.y = "subject_label", all.x = TRUE, all.y = FALSE)
ABCDS_CL <- ABCDS_CL[!is.na(ABCDS_CL$ds_vs_control_flag),]
ABCDS_CL <- merge(ABCDS_CL, ABCDS_apoe[, c("subject_label", "APOE_grouped")], by.x = "Subject", by.y = "subject_label", all = FALSE)
ABCDS_CL <- ABCDS_CL[!is.na(ABCDS_CL$APOE_grouped),]
ABCDS_CL <- ABCDS_CL[order(ABCDS_CL$Subject, ABCDS_CL$event_sequence),]
ABCDS_CL <- ABCDS_CL[!duplicated(ABCDS_CL$Subject),]
ABCDS_CL <- merge(ABCDS_CL, latency[, c("subject_label", "event_sequence", "amy_AgefromBaseline")],
                  by.x = c("Subject", "event_sequence"), by.y = c("subject_label", "event_sequence"), all = FALSE)


CL_thresh_ABCDS <- mean(ABCDS_CL[ABCDS_CL$age_at_visit <= 35 & !duplicated(ABCDS_CL$Subject),]$WUSTLcentiloid, na.rm = TRUE) + 
  1.96 *sd(ABCDS_CL[ABCDS_CL$age_at_visit <= 35 & !duplicated(ABCDS_CL$Subject),]$WUSTLcentiloid, na.rm = TRUE) 
ABCDS_Hipp <- readxl::read_xlsx(".././ABCDS_DF3_AmyloidOnly/ABCDS_amy_volumes.xlsx")
ABCDS_Hipp <- ABCDS_Hipp[, c("MRI", "Right-Hippocampus", "Left-Hippocampus", "EstimatedTotalIntraCranialVol")]
IDlist <- read.csv(".././ABCDS_DF3_AmyloidOnly/LONI_masterlist.csv")

ABCDS_Hipp$event_sequence <- substr(ABCDS_Hipp$MRI, start = nchar(ABCDS_Hipp$MRI), stop = nchar(ABCDS_Hipp$MRI))
ABCDS_Hipp$U19_ID <-  gsub("_.*", "", ABCDS_Hipp$MRI)
ABCDS_Hipp <- merge(ABCDS_Hipp, IDlist[, c("Subject.ID", "U19_ID")], by = "U19_ID", all.x = TRUE, all.y = FALSE)
ABCDS_Hipp <- ABCDS_Hipp[!duplicated(ABCDS_Hipp),]

ABCDS_Hipp$Subject.ID <- ifelse(is.na(ABCDS_Hipp$Subject.ID), ABCDS_Hipp$U19_ID, ABCDS_Hipp$Subject.ID)
ABCDS_Hipp$event_sequence <- as.factor(ABCDS_Hipp$event_sequence)
levels(ABCDS_Hipp$event_sequence) <- c("1", "2", "3", "1")
ABCDS_Hipp$hippo <- ABCDS_Hipp$`Right-Hippocampus` + ABCDS_Hipp$`Left-Hippocampus`
ABCDS_Hipp <- ABCDS_Hipp[order(ABCDS_Hipp$Subject.ID, ABCDS_Hipp$event_sequence),]
ABCDS_Hipp <- ABCDS_Hipp[!duplicated(ABCDS_Hipp$Subject.ID),]
ABCDS_Hipp <- merge(ABCDS_Hipp, ABCDS_apoe[, c("subject_label", "APOE_grouped")], by.x = "Subject.ID", by.y = "subject_label", all = FALSE)
ABCDS_Hipp <- merge(ABCDS_Hipp, latency[, c("subject_label", "event_sequence", "mri_AgefromBaseline")],
                    by.x = c("Subject.ID", "event_sequence"),
                    by.y = c("subject_label", "event_sequence"), all = FALSE)

mod <- lm(hippo ~ EstimatedTotalIntraCranialVol, data = ABCDS_Hipp, model = TRUE)
ABCDS_Hipp$hippo_resid <- NA
# Step 4: Insert residuals using the "na.action" attribute
used_rows <- as.numeric(rownames(model.frame(mod)))  # TRUE numeric row indices

ABCDS_Hipp$hippo_resid[used_rows] <- resid(mod)

ABCDS_Hipp$hippo_Z <- (ABCDS_Hipp$hippo_resid - mean(ABCDS_Hipp[ABCDS_Hipp$mri_AgefromBaseline < 35,]$hippo_resid)) /
  sd(ABCDS_Hipp[ABCDS_Hipp$mri_AgefromBaseline < 35,]$hippo_resid)

prefunc <- read.csv(".././ABCDS_LONI_20250528/Premorbid_Functioning_Level_28May2025.csv")
prefunc <- prefunc[, c("subject_label", "prefunclevel")]
prefunc <- prefunc[!is.na(prefunc$prefunclevel),]
prefunc <- prefunc[!duplicated(prefunc$subject_label),]

prefunc$Intellectual_Disability <- as.factor(prefunc$prefunclevel)
levels(prefunc$prefunclevel) <- c("Mild", "Moderate", "Severe")

ABCDS <- merge(ABCDS, prefunc, by.x = "SUBJECT_LABEL", by.y = "subject_label", all = FALSE)
ABCDS_CL <- merge(ABCDS_CL, prefunc, by.x = "Subject", by.y = "subject_label", all = FALSE)
GFAP <- merge(GFAP, prefunc, by.x = "SUBJECT_LABEL", by.y = "subject_label", all = FALSE)
NFL <- merge(NFL, prefunc, by.x = "SUBJECT_LABEL", by.y = "subject_label", all = FALSE)
ABCDS_Hipp <- merge(ABCDS_Hipp, prefunc, by.x = "Subject.ID", by.y = "subject_label", all = FALSE)
mcrt <- merge(mcrt, prefunc, by.x = "subject_label", by.y = "subject_label", all = FALSE)


source("./DABNI_maximizingAPOE.R")
DABNI_ID <- readxl::read_xlsx("./Data_202512/NPS_20251127.xlsx")
DABNI_ID <- DABNI_ID[, c("NHC", "Intellectual_Disability")]
DABNI_ID <- DABNI_ID[!duplicated(DABNI_ID$NHC),]
DABNI_CSF <- readxl::read_xlsx("./Data_202512/Biofluids_20251210.xlsx")
DABNI_CSF$CSFAB42_AB40 <- DABNI_CSF$CSF_AB142_A0010 / DABNI_CSF$CSF_AB140_A0010
DABNI_CSF <- DABNI_CSF[, c("NHC", "APOE", "Sex", "Age_at_CSF", "CSFAB42_AB40")]
# # DABNI <- DABNI[!DABNI$APOE == "24",]
DABNI$APOE_grouped <- recode(as.factor(DABNI$APOE),
                             "22" = "APOE2",
                             "23" = "APOE2",
                             "24" = "APOE2",
                             "33" = "APOE3",
                             "34" = "APOE4",
                             "43" = "APOE4",
                             "44" = "APOE4")
DABNI$APOE_grouped <- factor(DABNI$APOE_grouped, levels = c("APOE2", "APOE3", "APOE4"))

DABNI_CSF$APOE_grouped <- recode(as.factor(DABNI_CSF$APOE),
                             "22" = "APOE2",
                             "23" = "APOE2",
                             "24" = "APOE2",
                             "33" = "APOE3",
                             "34" = "APOE4",
                             "43" = "APOE4",
                             "44" = "APOE4")
DABNI_CSF$APOE_grouped <- factor(DABNI_CSF$APOE_grouped, levels = c("APOE2", "APOE3", "APOE4"))
# 
# 
# 
# DABNI_Hipp <- readxl::read_xlsx("./Data_202512/Hipp_Volume_T1_20251127.xlsx")
# DABNI_Hipp <- DABNI_Hipp[, c("NHC", "Sex", "Age_at_MRI", "FS_Left-Hippocampus", "FS_Right-Hippocampus", "FS_EstimatedTotalIntraCranialVol")]
DABNI_Hipp <- DABNI_Hipp[order(DABNI_Hipp$NHC, DABNI_Hipp$Age_at_MRI),]

# 
# DABNI_Hipp <- merge(DABNI_Hipp, DABNI[!duplicated(DABNI$NHC), c("NHC", "APOE_grouped")],
#                     by.x = "NHC", by.y = "NHC", all = FALSE)

DABNI_Hipp$cohort <- "DABNI"
DABNI_Hipp$hippo <- DABNI_Hipp$`FS_Left-Hippocampus` + DABNI_Hipp$`FS_Right-Hippocampus`
mod <- lm(hippo ~ FS_EstimatedTotalIntraCranialVol, data = DABNI_Hipp, model = TRUE)
DABNI_Hipp$hippo_resid <- NA
# Step 4: Insert residuals using the "na.action" attribute
used_rows <- as.numeric(rownames(model.frame(mod)))  # TRUE numeric row indices

DABNI_Hipp$hippo_resid[used_rows] <- resid(mod)
DABNI_Hipp$hippo_Z <- (DABNI_Hipp$hippo_resid - mean(DABNI_Hipp[DABNI_Hipp$Age_at_MRI < 35,]$hippo_resid)) /
  sd(DABNI_Hipp[DABNI_Hipp$Age_at_MRI < 35,]$hippo_resid)

ABCDS_Hipp$cohort <- "ABCDS"

colnames(DABNI_Hipp) <- c("Subject.ID", "mri_AgefromBaseline",
                          "Hippo_left", "hippo_right", "ICV", "APOE_grouped",
                          "de_gender", 
                          "cohort", "hippo", "hippo_resid", "hippo_Z")
DABNI_ID$Intellectual_Disability <- as.factor(DABNI_ID$Intellectual_Disability)
DABNI <- merge(DABNI, DABNI_ID, by = "NHC", all = FALSE)
DABNI_CL <- merge(DABNI_CL, DABNI_ID, by = "NHC", all = FALSE)
DABNI_CSF <- merge(DABNI_CSF, DABNI_ID, by = "NHC", all = FALSE)
DABNI_dx <- merge(DABNI_dx, DABNI_ID, by = "NHC", all = FALSE)
DABNI_Hipp <- merge(DABNI_Hipp, DABNI_ID, by.x = "Subject.ID", by.y = "NHC", all = FALSE)
################################################################################
################################################################################

#AB42/AB40
DABNI_CSF[!duplicated(DABNI_CSF$NHC),]
DABNI_CSF$Age_at_CSF <- as.numeric(DABNI_CSF$Age_at_CSF)
DABNI_CSF$NHC <- as.factor(DABNI_CSF$NHC)
DABNI_CSF$Sex <- as.factor(DABNI_CSF$Sex)
DABNI_CSF$APOE_grouped <- relevel(factor(DABNI_CSF$APOE_grouped), ref = "APOE3")

diagnostics_CSF_sex <- evaluate_gam_k_selection(DABNI_CSF, y_col="CSFAB42_AB40", age_col="Age_at_CSF", 
                                                apoe_col="APOE_grouped", sex_col = "Sex", intellectual_disability_col <- "Intellectual_Disability")


res_CSF_sex <- run_apoe_model_cohort(
  DABNI_CSF,
  apoe_col = "APOE_grouped",
  y_col = "CSFAB42_AB40",
  age_col = "Age_at_CSF",
  intellectual_disability_col = "Intellectual_Disability",
  cohort = NULL,
  sex_col = "Sex",
  reference_sex = "M",                                  
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  k_val = 4,
  markSig = FALSE
) 


export_gam_summary_to_excel(
  gam_summary = res_CSF_sex$model_obj, 
  file_path   = "./SuppMaterials/Supplemental_Table_CSF_DABNI_control_for_sex.xlsx", 
  model_label = "Supplemental Table X: GAM Trajectories of CSF AB42/AB40 by APOE and Sex and ID in DABNI"
)


CL <- readxl::read_xlsx("./Data_202512/CL_20251127.xlsx")
# CL <- CL[!CL$APOE == "24",]
CL$APOE_grouped <- recode(as.factor(CL$APOE),
                          "22" = "APOE2",
                          "23" = "APOE2",
                          "24" = "APOE2",
                          "33" = "APOE3",
                          "34" = "APOE4",
                          "43" = "APOE4",
                          "44" = "APOE4")
CL <- CL[!is.na(CL$CL_SPM12),]
CL <- CL[order(CL$NHC, CL$Age_at_PET),]
CL <- CL[!duplicated(CL$NHC),]

colnames(CL) <- c("Subject", "Session", "CL_drop", "WUSTLcentiloid",
                  "Tracer", "APOE", "amy_AgefromBaseline", "de_gender", "APOE_grouped")
CL$cohort <- "DABNI"
ABCDS_CL$de_gender <- as.factor(ABCDS_CL$de_gender)
levels(ABCDS_CL$de_gender) <- c("M", "F")
ABCDS_CL$cohort <- "ABCDS"
CL <- merge(CL, DABNI_ID, by.x = "Subject", by.y = "NHC", all = FALSE)
levels(ABCDS_CL$Intellectual_Disability) <- c("Mild", "Moderate", "Severe")
CL <- rbind(CL[, c("Subject","WUSTLcentiloid",
                   "amy_AgefromBaseline", "APOE_grouped", "cohort", "de_gender", "Intellectual_Disability")],
            ABCDS_CL[, c("Subject","WUSTLcentiloid",
                         "amy_AgefromBaseline", "APOE_grouped", "cohort", "de_gender", "Intellectual_Disability")])

CL$Subject <- as.factor(CL$Subject)
CL <- CL[!is.na(CL$amy_AgefromBaseline),]
CL$cohort <- factor(CL$cohort, levels = c("DABNI", "ABCDS"))
CL$APOE_grouped <- relevel(factor(CL$APOE_grouped), ref = "APOE3")

ggplot(CL, aes(x = WUSTLcentiloid, colour = cohort, group = cohort, fill = cohort)) + geom_histogram(alpha = 0.5, position = "identity") + theme_bw() +
  theme(legend.position = "bottom")


gam_select_CL_sex <- evaluate_gam_k_selection(df = CL, y_col = "WUSTLcentiloid", 
                                          age_col = "amy_AgefromBaseline", 
                                          apoe_col = "APOE_grouped", 
                                          intellectual_disability_col = "Intellectual_Disability",
                                          
                                          sex_col = "de_gender",
                                          cohort = "cohort") 


CL$WUSTLcentiloid_shift <- CL$WUSTLcentiloid - min(CL$WUSTLcentiloid, na.rm = TRUE) + 0.1


res_CL_cohort_sex <- run_apoe_model_cohort(
  df = CL,
  apoe_col = "APOE_grouped",
  age_col = "amy_AgefromBaseline",
  y_col = "WUSTLcentiloid_shift",
  intellectual_disability_col = "Intellectual_Disability",
  cohort = "cohort",
  reference_apoe = "APOE3",
  reference_sex = "M",  
  family_choice = Gamma(link = "log"),
  y_shift = min(CL$WUSTLcentiloid),
  
  sex = "de_gender",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c( "#3399FF", "#33CC33","#FF4D4D"),
  k_val = 4,
  markSig = TRUE
)
summary(res_CL_cohort_sex$model_obj)
p2 <- res_CL_cohort_sex[[1]] + xlab("Age") + ylab("Cortical Amyloid Burden [CL]") + scale_shape_manual(values = c(1, 4, 4))


export_gam_summary_to_excel(
  gam_summary = res_CL_cohort_sex$model_obj, 
  file_path   = "./SuppMaterials/Supplemental_Table_CL_cohort_control_for_sex.xlsx", 
  model_label = "Supplemental Table X: GAM Trajectories of Amyloid PET by APOE and Sex"
)

###############################################################################
##pTau217
###############################################################################

gam_select_ptau217_DABNI_sex <- evaluate_gam_k_selection(  df = DABNI,
                                                           apoe_col = "APOE_grouped",
                                                           age_col = "Age_at_CSF",
                                                           intellectual_disability_col = "Intellectual_Disability",
                                                           
                                                           y_col = "PLASMA_PTAU217_S0217",
                                                           sex_col = "Sex", k_range = 3:10)


res_ptau217_DABNI_sex <- run_apoe_model_cohort(
  df = DABNI,
  apoe_col = "APOE_grouped",
  age_col = "Age_at_CSF",
  y_col = "PLASMA_PTAU217_S0217",
  intellectual_disability_col = "Intellectual_Disability",
  reference_sex = "M",                                  
  family_choice = Gamma(link = "log"),
  sex_col = "Sex",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c( "#3399FF", "#33CC33","#FF4D4D"),
  markSig = TRUE,
  k_val = 5
)
# Set layout to a 2x2 grid
par(mfrow = c(2, 2))
# Plot all diagnostic charts
gam.check(res_ptau217_DABNI_sex$model_obj)
# Reset plotting layout back to standard 1x1
par(mfrow = c(1, 1))
summary(res_ptau217_DABNI_sex$model_obj)

export_gam_summary_to_excel(
  gam_summary = res_ptau217_DABNI_sex$model_obj, 
  file_path   = "./SuppMaterials/Supplemental_Table_pTau217_DABNI_control_for_sex.xlsx", 
  model_label = "Supplemental Table X: GAM Trajectories of Plasma pTau217 by APOE and Sex in DABNI"
)



pTau217 <- merge(ABCDS[!duplicated(ABCDS$SUBJECT_LABEL), c("SUBJECT_LABEL", "EVENT_SEQUENCE",
                                                           "pTau217", "APOE_grouped", "de_gender", "Intellectual_Disability")],
                 latency[, c("subject_label", "event_sequence", "clinical_AgefromBaseline")],
                 by.x = c("SUBJECT_LABEL", "EVENT_SEQUENCE"),
                 by.y = c("subject_label", "event_sequence"), all = FALSE)

gam_select_ptau217_ABCDS_sex <- evaluate_gam_k_selection(  df = pTau217,
                                                           apoe_col = "APOE_grouped",
                                                           age_col = "clinical_AgefromBaseline",
                                                           intellectual_disability_col = "Intellectual_Disability",
                                                           
                                                           sex_col = "de_gender",
                                                           y_col = "pTau217", k_range = 3:10)

pTau217$de_gender <- as.factor(pTau217$de_gender)
levels(pTau217$de_gender) <- c("M", "F")

res_ptau217_ABCDS_sex <- run_apoe_model_cohort(
  df = pTau217[!is.na(pTau217$de_gender),],
  apoe_col = "APOE_grouped",
  age_col = "clinical_AgefromBaseline",
  y_col = "pTau217",
  sex_col = "de_gender",
  intellectual_disability_col = "Intellectual_Disability",
  reference_sex = "M",                                  
  family_choice = Gamma(link = "log"),
  
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c( "#3399FF", "#33CC33","#FF4D4D"),
  markSig = TRUE,
  k_val = 4
)

# Set layout to a 2x2 grid
par(mfrow = c(2, 2))
# Plot all diagnostic charts
gam.check(res_ptau217_ABCDS_sex$model_obj)
# Reset plotting layout back to standard 1x1
par(mfrow = c(1, 1))

export_gam_summary_to_excel(
  gam_summary = res_ptau217_ABCDS_sex$model_obj, 
  file_path   = "./SuppMaterials/Supplemental_Table_pTau217_ABCDS_control_for_sex.xlsx", 
  model_label = "Supplemental Table X: GAM Trajectories of Plasma pTau217 by APOE and Sex in ABCDS"
)

#GFAP
DABNI_GFAP <- DABNI[!is.na(DABNI$PLASMA_GFAP_S0122),c("NHC", "Age_at_CSF", "Sex", "PLASMA_GFAP_S0122", "APOE_grouped", "Intellectual_Disability")]
levels(GFAP$Intellectual_Disability) <- c("Mild", "Moderate", "Severe")
GFAP <- merge(GFAP, ABCDS_apoe[, c("subject_label", "APOE_grouped")],
              by.x = "SUBJECT_LABEL", by.y = "subject_label", all = FALSE)
GFAP <- merge(GFAP, latency[, c("subject_label", "event_sequence", "clinical_AgefromBaseline")],
              by.x = c("SUBJECT_LABEL", "EVENT_SEQUENCE"),
              by.y = c("subject_label", "event_sequence"), all = FALSE)
GFAP <- merge(GFAP, demogs[!duplicated(demogs$subject_label), c("subject_label", "de_gender")],
              by.x = c("SUBJECT_LABEL"), by.y = c("subject_label"), all = FALSE)

colnames(DABNI_GFAP) <- c("SUBJECT_LABEL", "clinical_AgefromBaseline", "de_gender", "TESTVALUE", "APOE_grouped", "Intellectual_Disability")
DABNI_GFAP$cohort <- "DABNI"
GFAP$cohort <- "ABCDS"

colnames(GFAP)[c(14:15)] <- c("APOE_grouped","clinical_AgefromBaseline")

GFAP <- rbind(DABNI_GFAP, GFAP[!duplicated(GFAP$SUBJECT_LABEL), 
                               c("SUBJECT_LABEL", "clinical_AgefromBaseline", "de_gender", 
                                 "TESTVALUE", "APOE_grouped", "cohort", "Intellectual_Disability")])
GFAP$TESTVALUE <- as.numeric(GFAP$TESTVALUE)

GFAP$cohort <- factor(GFAP$cohort, levels = c("DABNI", "ABCDS"))
GFAP$de_gender <- as.factor(GFAP$de_gender)
levels(GFAP$de_gender) <- c("M", "F", "F", "M")



gam_select_GFAP_sex <- evaluate_gam_k_selection(df = GFAP, y_col = "TESTVALUE", 
                                            age_col = "clinical_AgefromBaseline", 
                                            intellectual_disability_col = "Intellectual_Disability",
                                            apoe_col = "APOE_grouped", sex_col = "de_gender",
                                            cohort = "cohort") 

res_GFAP_cohort_sex <- run_apoe_model_cohort(
  df = GFAP,
  apoe_col = "APOE_grouped",
  age_col = "clinical_AgefromBaseline",
  intellectual_disability_col = "Intellectual_Disability",
  y_col = "TESTVALUE",
  cohort = "cohort",
  reference_sex = "M",                                  
  family_choice = Gamma(link = "log"),
  
  sex_col = "de_gender",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = TRUE,
  
  k_val = 5
)


# Set layout to a 2x2 grid
par(mfrow = c(2, 2))
# Plot all diagnostic charts
gam.check(res_GFAP_cohort_sex$model_obj)
# Reset plotting layout back to standard 1x1
par(mfrow = c(1, 1))


summary(res_GFAP_cohort_sex$model_obj)

export_gam_summary_to_excel(
  gam_summary = res_GFAP_cohort_sex$model_obj, 
  file_path   = "./SuppMaterials/Supplemental_Table_GFAP_cohort_control_for_sex.xlsx", 
  model_label = "Supplemental Table X: GAM Trajectories of Plasma GFAP by APOE and Sex"
)


#NFL
DABNI_NFL <- DABNI[!is.na(DABNI$PLASMA_NFLIGHT_S0100),c("NHC", "Age_at_CSF", "PLASMA_NFLIGHT_S0100", "APOE_grouped", "Sex", "Intellectual_Disability")]
NFL <- merge(NFL[, c("SUBJECT_LABEL", "EVENT_SEQUENCE", "TESTVALUE")], ABCDS_apoe[, c("subject_label", "APOE_grouped")],
              by.x = "SUBJECT_LABEL", by.y = "subject_label", all = FALSE)
NFL <- merge(NFL, latency[, c("subject_label", "event_sequence", "clinical_AgefromBaseline")],
              by.x = c("SUBJECT_LABEL", "EVENT_SEQUENCE"),
              by.y = c("subject_label", "event_sequence"), all = FALSE)
NFL <- merge(NFL, demogs[!duplicated(demogs$subject_label),c("subject_label", "de_gender")],
             by.x = "SUBJECT_LABEL", by.y = "subject_label", all = FALSE)
NFL <- merge(NFL, prefunc, by.x = "SUBJECT_LABEL", by.y = "subject_label", all = FALSE)
levels(NFL$Intellectual_Disability) <- c("Mild", "Moderate", "Severe") 
colnames(DABNI_NFL) <- c("SUBJECT_LABEL", "clinical_AgefromBaseline", "TESTVALUE", "APOE_grouped", "de_gender", "Intellectual_Disability")
DABNI_NFL$cohort <- "DABNI"
NFL$cohort <- "ABCDS"

NFL <- rbind(DABNI_NFL, NFL[!duplicated(NFL$SUBJECT_LABEL), 
                               c("SUBJECT_LABEL", "clinical_AgefromBaseline", "TESTVALUE", "APOE_grouped", "de_gender", "cohort", "Intellectual_Disability")])
NFL$TESTVALUE <- as.numeric(NFL$TESTVALUE)
NFL$de_gender <- as.factor(NFL$de_gender)
levels(NFL$de_gender) <- c("M", "F", "F", "M")

NFL$TESTVALUE <- log(NFL$TESTVALUE)
NFL$cohort <- factor(NFL$cohort, levels = c("DABNI", "ABCDS"))


gam_select_NFL_sex <- evaluate_gam_k_selection(df = NFL, y_col = "TESTVALUE", 
                                           age_col = "clinical_AgefromBaseline", 
                                           apoe_col = "APOE_grouped", 
                                           intellectual_disability_col = "Intellectual_Disability",
                                           
                                           sex_col = "de_gender",
                                           cohort = "cohort", k_range = 3:10) 


res_NFL_cohort_sex <- run_apoe_model_cohort(
  df = NFL,
  apoe_col = "APOE_grouped",
  age_col = "clinical_AgefromBaseline",
  sex_col = "de_gender",
  y_col = "TESTVALUE",
  intellectual_disability_col = "Intellectual_Disability",
  reference_sex = "M",                                  
  
  cohort = "cohort",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF",  "#33CC33","#FF4D4D"),
  k_val = 4,
  markSig = TRUE,
  
)

# Set layout to a 2x2 grid
par(mfrow = c(2, 2))
# Plot all diagnostic charts
gam.check(res_NFL_cohort_sex$model_obj)
# Reset plotting layout back to standard 1x1
par(mfrow = c(1, 1))



summary(res_NFL_cohort_sex$model_obj)


export_gam_summary_to_excel(
  gam_summary = res_NFL_cohort_sex$model_obj, 
  file_path   = "./SuppMaterials/Supplemental_Table_NFL_cohort_control_for_sex.xlsx", 
  model_label = "Supplemental Table X: GAM Trajectories of Plasma NFL by APOE and Sex"
)


#Hippocampal volume
ABCDS_Hipp <- merge(ABCDS_Hipp, demogs[!duplicated(demogs$subject_label), c("subject_label", "de_gender")],
                    by.x = "Subject.ID", by.y = "subject_label", all = FALSE)
levels(ABCDS_Hipp$Intellectual_Disability) <- c("Mild", "Moderate", "Severe")
hippo <- rbind(DABNI_Hipp[, c("Subject.ID", "mri_AgefromBaseline",
                              "APOE_grouped", "de_gender",
                              "cohort",  "hippo_Z", "Intellectual_Disability")],
               ABCDS_Hipp[, c("Subject.ID", "mri_AgefromBaseline",
                              "APOE_grouped", "de_gender",
                              "cohort",  "hippo_Z", "Intellectual_Disability")])
hippo$cohort <- factor(hippo$cohort, levels = c("DABNI", "ABCDS"))
hippo$de_gender <- as.factor(hippo$de_gender)
levels(hippo$de_gender) <- c("M", "F", "F", "M")
hippo$APOE_grouped <- as.factor(hippo$APOE_grouped)
levels(hippo$APOE_grouped) <- c("APOE2", "APOE2", "APOE4",
                                "APOE3", "APOE4", "APOE4", 
                                "APOE4", "APOE3", "APOE2")

gam_select_hippo_sex <- evaluate_gam_k_selection(df = hippo, y_col = "hippo_Z", 
                                             age_col = "mri_AgefromBaseline", 
                                             apoe_col = "APOE_grouped", 
                                             intellectual_disability_col = "Intellectual_Disability",
                                             
                                             sex_col = "de_gender",
                                             cohort = "cohort", k_range = 3:10) 


res_hippo_cohort_sex <- run_apoe_model_cohort(
  df = hippo,
  apoe_col = "APOE_grouped",
  age_col = "mri_AgefromBaseline",
  y_col = "hippo_Z",
  cohort = "cohort",
  sex_col = "de_gender",
  reference_sex = "M",                                  
  
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  k_val = 4,
  markSig = TRUE
)

# Set layout to a 2x2 grid
par(mfrow = c(2, 2))
# Plot all diagnostic charts
gam.check(res_hippo_cohort_sex$model_obj)
# Reset plotting layout back to standard 1x1
par(mfrow = c(1, 1))

summary(res_hippo_cohort_sex$model_obj)


export_gam_summary_to_excel(
  gam_summary = res_hippo_cohort_sex$model_obj, 
  file_path   = "./SuppMaterials/Supplemental_Table_hippo_cohort_control_for_sex.xlsx", 
  model_label = "Supplemental Table X: GAM Trajectories of Hippocampal Volume by APOE and Sex"
)


#mcrt
mcrt <- mcrt[order(mcrt$subject_label, mcrt$event_sequence), c("subject_label",
                                                               "event_sequence", "trs")]
mcrt <- mcrt[!duplicated(mcrt$subject_label),]
mcrt <- merge(mcrt, ABCDS_apoe[, c("subject_label", "APOE_grouped")],
              by = "subject_label", all = FALSE)
mcrt <- merge(mcrt, demogs[!duplicated(demogs$subject_label), c("subject_label", "de_gender")],
              by.x = "subject_label", by.y = "subject_label", all = FALSE)

mcrt <- merge(mcrt, latency[, c("subject_label", "event_sequence", "clinical_AgefromBaseline")],
              by = c("subject_label", "event_sequence"), all = FALSE)
mcrt$cohort <- "ABCDS"
mcrt <- merge(mcrt, prefunc, by = "subject_label", all = FALSE)
levels(mcrt$Intellectual_Disability) <- c("Mild", "Moderate", "Severe")

DABNI_dx <- readxl::read_xlsx("./Data_202512/NPS_20251127.xlsx")
DABNI_dx <- DABNI_dx[, c("NHC", "Sex", "Age_NPS", "Diag", "crt_rli_rfi_total", "Intellectual_Disability")]
DABNI_dx <- merge(DABNI_dx, DABNI[, c("NHC", "APOE_grouped")], by = "NHC", all = FALSE)
DABNI_dx <- DABNI_dx[order(DABNI_dx$NHC, DABNI_dx$Age_NPS),]
DABNI_dx <- DABNI_dx[!duplicated(DABNI_dx$NHC),]
DABNI_dx$cohort <- "DABNI"
colnames(DABNI_dx) <- c("subject_label", "de_gender", "clinical_AgefromBaseline",
                        "diag", "trs", "Intellectual_Disability", "APOE_grouped", "cohort")


mcrt <- rbind(DABNI_dx[, c("subject_label",  "clinical_AgefromBaseline",
                           "trs", "APOE_grouped", "de_gender", "cohort", "Intellectual_Disability")],
              mcrt[, c("subject_label",  "clinical_AgefromBaseline",
                               "trs", "APOE_grouped", "de_gender", "cohort", "Intellectual_Disability")])
mcrt$trs <- as.numeric(mcrt$trs)
mcrt$cohort <- factor(mcrt$cohort, levels = c("DABNI", "ABCDS"))
mcrt$de_gender <- as.factor(mcrt$de_gender)
levels(mcrt$de_gender) <- c("M", "F", "F", "M")

gam_select_mcrt <- evaluate_gam_k_selection(df = mcrt, y_col = "trs", 
                                            age_col = "clinical_AgefromBaseline", 
                                            apoe_col = "APOE_grouped", 
                                            intellectual_disability_col = "Intellectual_Disability",
                                            cohort = "cohort", k_range = 3:10) 
gam_select_mcrt_sex <- evaluate_gam_k_selection(df = mcrt, y_col = "trs", 
                                            age_col = "clinical_AgefromBaseline", 
                                            apoe_col = "APOE_grouped", 
                                            intellectual_disability_col = "Intellectual_Disability",
                                            
                                            cohort = "cohort", 
                                            sex_col = "de_gender", k_range = 3:10)

# 1. Convert score to work with distribution
mcrt$score_prop <- 36 - mcrt$trs



res_mcrt_cohort_sex <- run_apoe_model_cohort(
  df = mcrt,
  apoe_col = "APOE_grouped",
  age_col = "clinical_AgefromBaseline",
  y_col = "score_prop",
  family = tw(link = "log"),
  cohort = "cohort",
  reference_sex = "M",                                  
  
  sex_col = "de_gender",
  intellectual_disability_col = "Intellectual_Disability",
  
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c( "#3399FF","#33CC33", "#FF4D4D"),
  YLIM= c(0, 90),
  k_val = 5,
  markSig = TRUE
)


# Extract model predictions
mcrt_pred <- res_mcrt_cohort_sex$predictions

# Transform predictions back to original TRS scale
mcrt_pred$trs_fit <- 36 - mcrt_pred$forecast
mcrt_pred$trs_min <- 36 - mcrt_pred$max
mcrt_pred$trs_max <- 36 - mcrt_pred$min

# Plot original data + model predictions on original TRS scale
p_mcrt_original <- ggplot2::ggplot(
  mcrt,
  ggplot2::aes(
    x = clinical_AgefromBaseline,
    y = trs,
    colour = APOE_grouped,
    group = APOE_grouped,
    shape = cohort
  )
) +
  ggplot2::geom_point(alpha = 0.3) +
  ggplot2::geom_ribbon(
    data = mcrt_pred,
    ggplot2::aes(
      x = clinical_AgefromBaseline,
      y = trs_fit,
      ymin = trs_min,
      ymax = trs_max,
      fill = APOE_grouped
    ),
    alpha = 0.4,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_segment(
    data = res_mcrt_cohort_sex$significance_summary,
    ggplot2::aes(
      x = start_age,
      xend = end_age,
      y = y - 55,
      yend = y - 55
    ),
    inherit.aes = FALSE,
    linewidth = 0.8,
    colour = "black"
  ) +
  ggplot2::geom_line(
    data = mcrt_pred,
    ggplot2::aes(
      x = clinical_AgefromBaseline,
      y = trs_fit,
      colour = APOE_grouped,
      group = APOE_grouped,
      shape = cohort
    ),
    inherit.aes = FALSE
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::scale_fill_manual(
    values = c("#3399FF", "#33CC33", "#FF4D4D")
  ) +
  ggplot2::scale_colour_manual(
    values = c("#3399FF", "#33CC33", "#FF4D4D")
  ) + ggplot2::scale_shape_manual(values = c(1, 4)) +
  ggplot2::scale_x_continuous(
    breaks = c(20, 30, 40, 50, 60, 70)
  ) +
  ggplot2::coord_cartesian(
    xlim = c(18, 72),
    ylim = c(0, 48)
  ) +
  ggplot2::xlab("Age") +
  ggplot2::ylab("modified Cued Recall Test (mCRT)")

p_mcrt_original



# Set layout to a 2x2 grid
par(mfrow = c(2, 2))
# Plot all diagnostic charts
gam.check(res_mcrt_cohort_sex$model_obj)
# Reset plotting layout back to standard 1x1
par(mfrow = c(1, 1))


export_gam_summary_to_excel(
  gam_summary = res_mcrt_cohort_sex$model_obj, 
  file_path   = "./SuppMaterials/Supplemental_Table_mcrt_cohort_control_for_sex_and_ID.xlsx", 
  model_label = "Supplemental Table X: GAM Trajectories of mCRT by APOE and Sex"
)



lemon::grid_arrange_shared_legend(res_CSF_sex$plot + aes(shape = APOE_grouped) + 
                                    scale_shape_manual(values = c(1, 1, 1), guide = "none") + ggtitle("A."), 
                                  res_CL_cohort_sex$plot + ggtitle("B."), 
                                  res_ptau217_DABNI_sex$plot + aes(shape = APOE_grouped) + 
                                    scale_shape_manual(values = c(1, 1, 1), guide = "none") +  ggtitle("C."),
                                  res_ptau217_ABCDS_sex$plot + 
                                    aes(shape = APOE_grouped) + 
                                    scale_shape_manual(values = c(4, 4, 4), guide = "none") + 
                                    ggtitle("D."),
                                  res_GFAP_cohort_sex$plot + ggtitle("E."), res_NFL_cohort_sex$plot + ggtitle("F."), 
                                  res_hippo_cohort_sex$plot + ggtitle("G."), p_mcrt_original + ggtitle("H."), nrow = 4, ncol = 2)

graph2ppt(file = "./figures/Fig1.pptx", width = 6, height = 11)

# res_mcrt_cohort_censoredModel$plot
# graph2ppt(file = "./figures/SuppMaterials_BetaRegressionformCRT.pptx", width = 6, height = 11)

grid.arrange(diagnostics_CSF_sex$plots$aic_plot + ggtitle("A. CSF AB42/AB40, DABNI"), 
                                  gam_select_CL_sex$plots$aic_plot + ggtitle("B. Amyloid PET"), 
                                  gam_select_ptau217_DABNI_sex$plots$aic_plot +  ggtitle("C. Plasma pTau217, DABNI"),
                                  gam_select_ptau217_ABCDS_sex$plots$aic_plot + 
                                    ggtitle("D. Plasma pTau217, ABCDS"),
                                  gam_select_GFAP_sex$plots$aic_plot + ggtitle("E. Plasma GFAP"), gam_select_NFL_sex$plots$aic_plot + ggtitle("F. Plasma NFL"), 
                                  gam_select_hippo_sex$plots$aic_plot + ggtitle("G. Hippocampal Volume"), gam_select_mcrt_sex$plots$aic_plot + 
                                    ggtitle("H. modified Cued Recall Test (mCRT)"), nrow = 4, ncol = 2)

graph2ppt(file = "./SuppMaterials/k_selection_plots.pptx", width = 6, height = 11)


###############################################################################
##Stratified by sex
###############################################################################


pTau217$SUBJECT_LABEL <- as.character(pTau217$SUBJECT_LABEL)
pTau217$EVENT_SEQUENCE <- as.numeric(pTau217$EVENT_SEQUENCE)
pTau217 <- merge(pTau217, demogs, by.x = c("SUBJECT_LABEL", "EVENT_SEQUENCE"),
                 by.y = c("subject_label", "event_sequence"), all = FALSE)
pTau217$Sex <- as.factor(pTau217$de_gender.y)
levels(pTau217$Sex) <- c("M", "F")


res_CSF_by_sex <- compare_sexes_by_apoe(
  df = DABNI_CSF,
  apoe_col = "APOE_grouped",
  sex_col = "Sex",                                 # Changed to primary stratifying variable
  y_col = "CSFAB42_AB40",
  age_col = "Age_at_CSF",

  # cohort = "cohort",
  fill_values = c("F" = "#FF4D4D", "M" = "#3399FF"),   # Default sex-stratified palette
  colour_values = c("F" = "#FF4D4D", "M" = "#3399FF"),
  sig_bar_colors = c("sig_APOE2" = "purple", "sig_APOE3" = "orange", "sig_APOE4" = "darkred"),
  markSig = TRUE,
  k_val = 4,
  reference_cohort = NULL,
  XLIM = c(18, 72),
  YLIM = FALSE
)


res_CL_by_sex <- compare_sexes_by_apoe(
  df = CL,
  apoe_col = "APOE_grouped",
  sex_col = "de_gender",                                 # Changed to primary stratifying variable
  y_col = "WUSTLcentiloid_shift",
  family_choice = Gamma(link = "log"),
  age_col = "amy_AgefromBaseline",
  cohort = "cohort",
  fill_values = c("F" = "#FF4D4D", "M" = "#3399FF"),   # Default sex-stratified palette
  colour_values = c("F" = "#FF4D4D", "M" = "#3399FF"),
  sig_bar_colors = c("sig_APOE2" = "purple", "sig_APOE3" = "orange", "sig_APOE4" = "darkred"),
  markSig = TRUE,
  k_val = 4,
  reference_cohort = NULL,
  XLIM = c(18, 72),
  YLIM = FALSE
)

res_pTau217_DABNI_by_sex <- compare_sexes_by_apoe(
  df = DABNI,
  apoe_col = "APOE_grouped",
  sex_col = "Sex",                                 # Changed to primary stratifying variable
  y_col = "PLASMA_PTAU217_S0217",
  age_col = "Age_at_CSF",
  family_choice = Gamma(link = "log"),
  fill_values = c("F" = "#FF4D4D", "M" = "#3399FF"),   # Default sex-stratified palette
  colour_values = c("F" = "#FF4D4D", "M" = "#3399FF"),
  sig_bar_colors = c("sig_APOE2" = "purple", "sig_APOE3" = "orange", "sig_APOE4" = "darkred"),
  markSig = TRUE,
  k_val = 4,
  reference_cohort = NULL,
  XLIM = c(18, 72),
  YLIM = FALSE
)

res_pTau217_ABCDS_by_sex <- compare_sexes_by_apoe(
  df = pTau217,
  apoe_col = "APOE_grouped",
  sex_col = "Sex",                                 # Changed to primary stratifying variable
  y_col = "pTau217",
  age_col = "clinical_AgefromBaseline",
  fill_values = c("F" = "#FF4D4D", "M" = "#3399FF"),   # Default sex-stratified palette
  colour_values = c("F" = "#FF4D4D", "M" = "#3399FF"),
  sig_bar_colors = c("sig_APOE2" = "purple", "sig_APOE3" = "orange", "sig_APOE4" = "darkred"),
  family_choice = Gamma(link = "log"),
  markSig = TRUE,
  k_val = 4,
  reference_cohort = NULL,
  XLIM = c(18, 72),
  YLIM = FALSE
)


res_GFAP_by_sex <- compare_sexes_by_apoe(
  df = GFAP,
  apoe_col = "APOE_grouped",
  sex_col = "de_gender",                                 # Changed to primary stratifying variable
  y_col = "TESTVALUE",
  age_col = "clinical_AgefromBaseline",
  cohort = "cohort",
  fill_values = c("F" = "#FF4D4D", "M" = "#3399FF"),   # Default sex-stratified palette
  colour_values = c("F" = "#FF4D4D", "M" = "#3399FF"),
  sig_bar_colors = c("sig_APOE2" = "purple", "sig_APOE3" = "orange", "sig_APOE4" = "darkred"),
  markSig = TRUE,
  family_choice = Gamma(link = "log"),
  k_val = 5,
  reference_cohort = NULL,
  XLIM = c(18, 72),
  YLIM = FALSE
)


res_NFL_by_sex <- compare_sexes_by_apoe(
  df = NFL,
  apoe_col = "APOE_grouped",
  sex_col = "de_gender",                                 # Changed to primary stratifying variable
  y_col = "TESTVALUE",
  age_col = "clinical_AgefromBaseline",
  cohort = "cohort",
  fill_values = c("F" = "#FF4D4D", "M" = "#3399FF"),   # Default sex-stratified palette
  colour_values = c("F" = "#FF4D4D", "M" = "#3399FF"),
  sig_bar_colors = c("sig_APOE2" = "purple", "sig_APOE3" = "orange", "sig_APOE4" = "darkred"),
  markSig = TRUE,
  k_val = 4,
  reference_cohort = NULL,
  XLIM = c(18, 72),
  YLIM = FALSE
)

res_hippo_by_sex <- compare_sexes_by_apoe(
  df = hippo,
  apoe_col = "APOE_grouped",
  sex_col = "de_gender",                                 # Changed to primary stratifying variable
  y_col = "hippo_Z",
  age_col = "mri_AgefromBaseline",
  cohort = "cohort",
  fill_values = c("F" = "#FF4D4D", "M" = "#3399FF"),   # Default sex-stratified palette
  colour_values = c("F" = "#FF4D4D", "M" = "#3399FF"),
  sig_bar_colors = c("sig_APOE2" = "purple", "sig_APOE3" = "orange", "sig_APOE4" = "darkred"),
  markSig = TRUE,
  k_val = 3,
  reference_cohort = NULL,
  XLIM = c(18, 72),
  YLIM = FALSE
)

res_mcrt_by_sex <- compare_sexes_by_apoe(
  df = mcrt,
  apoe_col = "APOE_grouped",
  sex_col = "de_gender",                                 # Changed to primary stratifying variable
  y_col = "score_prop",
  age_col = "clinical_AgefromBaseline",
  cohort = "cohort",
  fill_values = c("F" = "#FF4D4D", "M" = "#3399FF"),   # Default sex-stratified palette
  colour_values = c("F" = "#FF4D4D", "M" = "#3399FF"),
  sig_bar_colors = c("sig_APOE2" = "purple", "sig_APOE3" = "orange", "sig_APOE4" = "darkred"),
  markSig = TRUE,
  k_val = 7,
  family_choice = tw(link = "log"),
  reference_cohort = NULL,
  XLIM = c(18, 72),
  YLIM = c(0, 46)
)



# Extract model predictions
mcrt_pred_by_sex <- res_mcrt_by_sex$predictions

# Transform predictions back to original TRS scale
# Transform predictions back to raw 0-36 mCRT scale
mcrt_pred_by_sex$trs_fit <- 36 - mcrt_pred_by_sex$forecast
mcrt_pred_by_sex$trs_min <- 36 - mcrt_pred_by_sex$max  # max deficit gives min score
mcrt_pred_by_sex$trs_max <- 36 - mcrt_pred_by_sex$min  # min deficit gives max score

# Prepare significance summary for plotting above the raw data
sig_summary_by_sex <- res_mcrt_by_sex$significance_summary
if (!is.null(sig_summary_by_sex) && nrow(sig_summary_by_sex) > 0) {
  # Stack bars above 36 (e.g. 38, 40, 42)
  sig_summary_by_sex$y_plot <- 37 + (sig_summary_by_sex$grp * 2)
}

p_mcrt_original_by_sex <- ggplot2::ggplot(
  mcrt,
  ggplot2::aes(
    x = clinical_AgefromBaseline,
    y = trs,
    colour = de_gender,
    group = interaction(APOE_grouped, de_gender),
    shape = cohort
  )
) +
  # Points for raw observations
  ggplot2::geom_point(alpha = 0.3) +
  
  # Fitted confidence ribbons back on raw mCRT scale
  ggplot2::geom_ribbon(
    data = mcrt_pred_by_sex,
    ggplot2::aes(
      x = clinical_AgefromBaseline,
      y = trs_fit,
      ymin = trs_min,
      ymax = trs_max,
      fill = de_gender,
      group = interaction(APOE_grouped, de_gender)
    ),
    alpha = 0.25,
    inherit.aes = FALSE
  ) +
  
  # GAM trajectory lines
  ggplot2::geom_line(
    data = mcrt_pred_by_sex,
    ggplot2::aes(
      x = clinical_AgefromBaseline,
      y = trs_fit,
      colour = de_gender,
      group = interaction(APOE_grouped, de_gender)
    ),
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  
  # Significance bars positioned neatly above y = 36
  ggplot2::geom_segment(
    data = sig_summary_by_sex,
    ggplot2::aes(
      x = start_age,
      xend = end_age,
      y = y_plot,
      yend = y_plot
    ),
    inherit.aes = FALSE,
    linewidth = 0.8,
    colour = "black"
  ) +
  
  # Layout and formatting
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  ggplot2::scale_fill_manual(values = c("F" = "#FF4D4D", "M" = "#3399FF")) +
  ggplot2::scale_colour_manual(values = c("F" = "#FF4D4D", "M" = "#3399FF")) +
  ggplot2::scale_shape_manual(values = c(1, 4)) +
  ggplot2::scale_x_continuous(breaks = c(20, 30, 40, 50, 60, 70)) +
  
  # Facet by APOE status so female vs male trajectories are clear per genotype
  ggplot2::facet_wrap(~ APOE_grouped) +
  
  # Set viewport range (0 to 46 gives space for the significance bars above 36)
  ggplot2::coord_cartesian(
    xlim = c(18, 72),
    ylim = c(0, 46)
  ) +
  ggplot2::xlab("Age") +
  ggplot2::ylab("modified Cued Recall Test (mCRT)")

print(p_mcrt_original_by_sex)

lemon::grid_arrange_shared_legend(res_CSF_by_sex$plot  + aes(shape = APOE_grouped) + 
                                    scale_shape_manual(values = c(1, 1, 1), guide = "none") + ggtitle("A."),
                                  res_CL_by_sex$plot   + ggtitle("B."),
                                  res_pTau217_DABNI_by_sex$plot +  aes(shape = APOE_grouped) + 
                                    scale_shape_manual(values = c(1, 1, 1), guide = "none") + ggtitle("C."),
                                  res_pTau217_ABCDS_by_sex$plot +  aes(shape = APOE_grouped) + 
                                    scale_shape_manual(values = c(4, 4, 4), guide = "none") + ggtitle("D."),
                                  res_GFAP_by_sex$plot + ggtitle("E."),
                                  res_NFL_by_sex$plot + ggtitle("F."),
                                  res_hippo_by_sex$plot + ggtitle("G."),
                                  p_mcrt_original_by_sex + ggtitle("H."), nrow = 4, ncol = 2)
graph2ppt(file = "./figures/Fig2_SexStratified_FacetPlots_withBars.pptx", width = 6, height = 11)


DABNI$Sex <- as.factor(DABNI$Sex)
CL$de_gender <- as.factor(CL$de_gender)
DABNI_CSF$Sex <- as.factor(DABNI_CSF$Sex)

apoe_sex_CSF <- run_apoe_sex_model_cohort(
    df = DABNI_CSF,
    apoe_col = "APOE_grouped",
    sex_col = "Sex",                                 # Changed to primary stratifying variable
    y_col = "CSFAB42_AB40",
    age_col = "Age_at_CSF",
    reference_sex = "M",                                  

    # cohort = "cohort",
    colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),                                        # 3 colors mapped to APOE2, APOE3, APOE4 names
    k_val = 4,
    reference_cohort = NULL,
    XLIM = c(18, 72),                                     
    YLIM = FALSE
)

apoe_sex_CSF$significance_summary$y <-0.01 *  seq_len(nrow(apoe_sex_CSF$significance_summary)) + 0.1
apoe_sex_CSF$plot <- apoe_sex_CSF$plot +
  geom_segment(data = data.frame(apoe_sex_CSF$significance_summary),
               aes(x = start_age, xend = end_age,
                   y = y, yend = y),
               linewidth = 1, inherit.aes = FALSE)

export_gam_summary_to_excel(
  gam_summary = apoe_sex_CSF$model_obj, 
  file_path   = "./SuppMaterials/Supplemental_Table_CSF_by_APOE_and_sex.xlsx", 
  model_label = "Supplemental Table X: GAM Trajectories of CSF AB42/AB40 by APOE and Sex"
)

# Set layout to a 2x2 grid
par(mfrow = c(2, 2))
# Plot all diagnostic charts
gam.check(apoe_sex_CSF$model_obj)
# Reset plotting layout back to standard 1x1
par(mfrow = c(1, 1))

apoe_sex_CL <- run_apoe_sex_model_cohort(
  df = CL,
  apoe_col = "APOE_grouped",
  sex_col = "de_gender",                                 # Changed to primary stratifying variable
  y_col = "WUSTLcentiloid_shift",
  age_col = "amy_AgefromBaseline",
  cohort = "cohort",
  reference_sex = "M",                                  
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),                                        # 3 colors mapped to APOE2, APOE3, APOE4 names
  k_val = 4,
  family_choice = Gamma(link = "log"),
  y_shift = min(CL$WUSTLcentiloid),
  reference_cohort = NULL,
  XLIM = c(18, 72),                                     
  YLIM = FALSE
)
export_gam_summary_to_excel(
  gam_summary = apoe_sex_CL$model_obj, 
  file_path   = "./SuppMaterials/Supplemental_Table_CL_by_APOE_and_sex.xlsx", 
  model_label = "Supplemental Table X: GAM Trajectories of Amyloid PET by APOE and Sex"
)


# apoe_sex_CL$significance_summary$y <-2 *  seq_len(nrow(apoe_sex_CL$significance_summary)) + 80
# apoe_sex_CL$plot <- apoe_sex_CL$plot +
#   geom_segment(data = data.frame(apoe_sex_CL$significance_summary),
#                aes(x = start_age, xend = end_age,
#                    y = y, yend = y),
#                linewidth = 1, inherit.aes = FALSE)

# Set layout to a 2x2 grid
par(mfrow = c(2, 2))
# Plot all diagnostic charts
gam.check(apoe_sex_CL$model_obj)
# Reset plotting layout back to standard 1x1
par(mfrow = c(1, 1))

apoe_sex_pTau217_DABNI <- run_apoe_sex_model_cohort(
  df = DABNI,
  apoe_col = "APOE_grouped",
  sex_col = "Sex",                                 # Changed to primary stratifying variable
  y_col = "PLASMA_PTAU217_S0217",
  age_col = "Age_at_CSF",
  reference_sex = "M",                                  
  family_choice = Gamma(link = "log"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),                                        # 3 colors mapped to APOE2, APOE3, APOE4 names
  k_val = 4,
  reference_cohort = NULL,
  XLIM = c(18, 72),                                     
  YLIM = FALSE
)
export_gam_summary_to_excel(
  gam_summary = apoe_sex_pTau217_DABNI$model_obj, 
  file_path   = "./SuppMaterials/Supplemental_Table_PlasmapTau217_DABNI_by_APOE_and_sex.xlsx", 
  model_label = "Supplemental Table X: GAM Trajectories of Plasma pTau217 in DABNI by APOE and Sex"
)




# Set layout to a 2x2 grid
par(mfrow = c(2, 2))
# Plot all diagnostic charts
gam.check(apoe_sex_pTau217_DABNI$model_obj)
# Reset plotting layout back to standard 1x1
par(mfrow = c(1, 1))


apoe_sex_pTau217_ABCDS <- run_apoe_sex_model_cohort(
  df = pTau217,
  apoe_col = "APOE_grouped",
  sex_col = "Sex",                                 # Changed to primary stratifying variable
  y_col = "pTau217",
  reference_sex = "M",                                  
  age_col = "clinical_AgefromBaseline",
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),                                        # 3 colors mapped to APOE2, APOE3, APOE4 names
  k_val = 4,
  reference_cohort = NULL,
  XLIM = c(18, 72),                                     
  YLIM = FALSE
)
export_gam_summary_to_excel(
  gam_summary = apoe_sex_pTau217_ABCDS$model_obj, 
  file_path   = "./SuppMaterials/Supplemental_Table_PlasmapTau217_ABCDS_by_APOE_and_sex.xlsx", 
  model_label = "Supplemental Table X: GAM Trajectories of Plasma pTau217 in ABCDS by APOE and Sex"
)

apoe_sex_pTau217_ABCDS$significance_summary$y <-0.5 *  seq_len(nrow(apoe_sex_pTau217_ABCDS$significance_summary)) + 2
apoe_sex_pTau217_ABCDS$plot <- apoe_sex_pTau217_ABCDS$plot +
  geom_segment(data = data.frame(apoe_sex_pTau217_ABCDS$significance_summary),
               aes(x = start_age, xend = end_age,
                   y = y, yend = y),
               linewidth = 1, inherit.aes = FALSE)

# Set layout to a 2x2 grid
par(mfrow = c(2, 2))
# Plot all diagnostic charts
gam.check(apoe_sex_pTau217_ABCDS$model_obj)
# Reset plotting layout back to standard 1x1
par(mfrow = c(1, 1))
apoe_sex_GFAP <- run_apoe_sex_model_cohort(
  df = GFAP,
  apoe_col = "APOE_grouped",
  sex_col = "de_gender",                                 # Changed to primary stratifying variable
  y_col = "TESTVALUE",
  reference_sex = "M",                                  
  family_choice = Gamma(link = "log"),
  age_col = "clinical_AgefromBaseline",
  cohort = "cohort",
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),                                        # 3 colors mapped to APOE2, APOE3, APOE4 names
  k_val = 4,
  reference_cohort = NULL,
  XLIM = c(18, 72),                                     
  YLIM = c(0, 1000)
)
export_gam_summary_to_excel(
  gam_summary = apoe_sex_GFAP$model_obj, 
  file_path   = "./SuppMaterials/Supplemental_Table_PlasmaGFAP_by_APOE_and_sex.xlsx", 
  model_label = "Supplemental Table X: GAM Trajectories of Plasma GFAP by APOE and Sex"
)

apoe_sex_GFAP$significance_summary$y <-50 *  seq_len(nrow(apoe_sex_GFAP$significance_summary)) + 500
apoe_sex_GFAP$plot <- apoe_sex_GFAP$plot +
  geom_segment(data = data.frame(apoe_sex_GFAP$significance_summary),
               aes(x = start_age, xend = end_age,
                   y = y, yend = y),
               linewidth = 1, inherit.aes = FALSE)


# Set layout to a 2x2 grid
par(mfrow = c(2, 2))
# Plot all diagnostic charts
gam.check(apoe_sex_GFAP$model_obj)
# Reset plotting layout back to standard 1x1
par(mfrow = c(1, 1))
apoe_sex_NFL <- run_apoe_sex_model_cohort(
  df = NFL,
  apoe_col = "APOE_grouped",
  sex_col = "de_gender",                                 # Changed to primary stratifying variable
  y_col = "TESTVALUE",
  age_col = "clinical_AgefromBaseline",
  cohort = "cohort",
  reference_sex = "M",                                  
  
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),                                        # 3 colors mapped to APOE2, APOE3, APOE4 names
  k_val = 4,
  reference_cohort = NULL,
  XLIM = c(18, 72),                                     
  YLIM = FALSE
)
export_gam_summary_to_excel(
  gam_summary = apoe_sex_NFL$model_obj, 
  file_path   = "./SuppMaterials/Supplemental_Table_PlasmaNFL_by_APOE_and_sex.xlsx", 
  model_label = "Supplemental Table X: GAM Trajectories of Plasma NFL by APOE and Sex"
)

apoe_sex_NFL$significance_summary$y <-0.5 *  seq_len(nrow(apoe_sex_NFL$significance_summary)) + 4
apoe_sex_NFL$plot <- apoe_sex_NFL$plot +
  geom_segment(data = data.frame(apoe_sex_NFL$significance_summary),
               aes(x = start_age, xend = end_age,
                   y = y, yend = y),
               linewidth = 1, inherit.aes = FALSE)

# Set layout to a 2x2 grid
par(mfrow = c(2, 2))
# Plot all diagnostic charts
gam.check(apoe_sex_NFL$model_obj)
# Reset plotting layout back to standard 1x1
par(mfrow = c(1, 1))

apoe_sex_hippo <- run_apoe_sex_model_cohort(
  df = hippo,
  apoe_col = "APOE_grouped",
  sex_col = "de_gender",                                 # Changed to primary stratifying variable
  y_col = "hippo_Z",
  age_col = "mri_AgefromBaseline",
  cohort = "cohort",
  reference_sex = "M",                                  
  
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),                                        # 3 colors mapped to APOE2, APOE3, APOE4 names
  k_val = 3,
  reference_cohort = NULL,
  XLIM = c(18, 72),                                     
  YLIM = FALSE
)
export_gam_summary_to_excel(
  gam_summary = apoe_sex_hippo$model_obj, 
  file_path   = "./SuppMaterials/Supplemental_Table_hippo_by_APOE_and_sex.xlsx", 
  model_label = "Supplemental Table X: GAM Trajectories of Hippocampal Volume by APOE and Sex"
)

# Set layout to a 2x2 grid
par(mfrow = c(2, 2))
# Plot all diagnostic charts
gam.check(apoe_sex_hippo$model_obj)
# Reset plotting layout back to standard 1x1
par(mfrow = c(1, 1))

apoe_sex_hippo$significance_summary$y <-0.25 *  seq_len(nrow(apoe_sex_hippo$significance_summary)) + 0.5
apoe_sex_hippo$plot <- apoe_sex_hippo$plot +
  geom_segment(data = data.frame(apoe_sex_hippo$significance_summary),
               aes(x = start_age, xend = end_age,
                   y = y, yend = y),
               linewidth = 1, inherit.aes = FALSE)



apoe_sex_mcrt <- run_apoe_sex_model_cohort(
  df = mcrt,
  apoe_col = "APOE_grouped",
  sex_col = "de_gender",                                 # Changed to primary stratifying variable
  y_col = "score_prop",
  age_col = "clinical_AgefromBaseline",
  cohort = "cohort",
  reference_sex = "M",                                  
  family_choice = tw(link = "log"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),                                        # 3 colors mapped to APOE2, APOE3, APOE4 names
  k_val = 7,
  reference_cohort = NULL,
  XLIM = c(18, 72),
  YLIM = c(0, 5)
)

export_gam_summary_to_excel(
  gam_summary = apoe_sex_mcrt$model_obj, 
  file_path   = "./SuppMaterials/Supplemental_Table_mcrt_by_APOE_and_sex.xlsx", 
  model_label = "Supplemental Table X: GAM Trajectories of mCRT by APOE and Sex"
)


# Extract model predictions
mcrt_pred_sex <- apoe_sex_mcrt$predictions

# Transform predictions back to original TRS scale
mcrt_pred_sex$trs_fit <- 36 - mcrt_pred_sex$forecast
mcrt_pred_sex$trs_min <- 36 - mcrt_pred_sex$max
mcrt_pred_sex$trs_max <- 36 - mcrt_pred_sex$min
apoe_sex_mcrt$significance_summary$y <- seq_len(nrow(apoe_sex_mcrt$significance_summary)) + 40

# Plot original data + model predictions on original TRS scale
p_mcrt_original <- ggplot(mcrt,
  aes(x = clinical_AgefromBaseline, y = trs, colour = APOE_grouped,
    group = APOE_grouped, shape = cohort)) +
  geom_point(alpha = 0.3) +
  geom_segment(data = apoe_sex_mcrt$significance_summary,
      aes(x = start_age, xend = end_age,
      y = y,  yend = y),
    inherit.aes = FALSE, linewidth = 1.2,
    colour = "black") +
  geom_line(data = mcrt_pred_sex, aes(
      x = clinical_AgefromBaseline, y = trs_fit,
      colour = APOE_grouped, group = apoe_sex,
        linetype = de_gender
    ),
    inherit.aes = FALSE
  ) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_fill_manual(
    values = c("#3399FF", "#33CC33", "#FF4D4D")
  ) +
  scale_colour_manual(
    values = c("#3399FF", "#33CC33", "#FF4D4D")
  ) + scale_shape_manual(values = c(1, 4)) +
  scale_x_continuous(breaks = c(20, 30, 40, 50, 60, 70)) +
  coord_cartesian( xlim = c(18, 72), ylim = c(0, 48)) +
  xlab("Age") +
  ylab("modified Cued Recall Test (mCRT)")





# Set layout to a 2x2 grid
par(mfrow = c(2, 2))
# Plot all diagnostic charts
gam.check(apoe_sex_mcrt$model_obj)
# Reset plotting layout back to standard 1x1
par(mfrow = c(1, 1))
lemon::grid_arrange_shared_legend(apoe_sex_CSF$plot, apoe_sex_CL$plot,
             apoe_sex_pTau217_DABNI$plot, apoe_sex_pTau217_ABCDS$plot,
             apoe_sex_GFAP$plot, apoe_sex_NFL$plot,
             apoe_sex_hippo$plot, p_mcrt_original, nrow = 4, ncol = 2)
graph2ppt(file = "./figures/Fig2_SexStratified.pptx", width = 6, height = 11)
