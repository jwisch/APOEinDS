library(ggplot2)
library(gridExtra)
library(dplyr)
library(fuzzyjoin)
library(export)
library(ggsignif)
library(data.table)
library(patchwork)

#https://alz-journals.onlinelibrary.wiley.com/doi/10.1002/alz.13859
source("./slidingWindowFuncs.R")
source("./DABNI_maximizingAPOE.R")
#Core 1, Amyloid: Amyloid PET
#Core 1, T1: p-Tau217
#Core 2, T2: Tau PET, ABCDS Only
#I: GFAP
#N: Anatomic MRI - Hippocampal Volume
visit_latency <- read.csv("../ABCDS_LONI_20240702/Age_at_Event_and_Latency_11Jul2024.csv")
demogs <- read.csv("../ABCDS_LONI_20250528/Participant_Demographics_28May2025.csv")
demogs <- demogs[, c("subject_label", "de_gender")]
demogs <- demogs[!duplicated(demogs),]

mcrt <- read.csv(".././ABCDS_LONI_20250528/Cued_Recall_28May2025.csv") #Krasny says 28 for mci and 23 for dementia

ABCDS_corr <- read.csv(".././ABCDS_LONI_20250528/freeze_5_corrections_2025.07.31.csv")
ABCDS_apoe <- read.csv(".././ABCDS_LONI_20250528/ApoE_Genotyping_Results_28May2025.csv")
ABCDS_apoe$APOE_grouped <- ifelse(ABCDS_apoe$allele_combo %in% c("E2/E4", "E4/E2", "E3/E4", "E4/E4", "E4/E3"), "APOE4",
                                  ifelse(ABCDS_apoe$allele_combo %in% c("E3/E3", "E3/3E"), "APOE3", 
                                         ifelse(ABCDS_apoe$allele_combo %in% c("E2/E2", "E2/E3", "E3/E2"), "APOE2", NA)))
df_ABCDS <- read.csv(".././ABCDS_LONI_20250528/abcds_biospecimen_results_28May2025.csv")
dx <- read.csv(".././ABCDS_LONI_20250528/Consensus_28May2025.csv")
dx <- dx[, c("subject_label", "event_sequence", "consensus_dx")]
dx <- merge(dx, ABCDS_apoe[, c("subject_label", "APOE_grouped")], by = "subject_label", all = FALSE)
dx <- merge(dx, demogs, by = "subject_label", all.x = TRUE, all.y = FALSE)

# DABNI_dx <- readxl::read_xlsx("./Data_202512/NPS_20251127.xlsx")
# DABNI_dx <- DABNI_dx[, c("NHC", "Sex", "Age_NPS", "Diag", "crt_rli_rfi_total")]

ABCDS <- df_ABCDS[df_ABCDS$TESTNAME %in% c("P-tau217"),]
ABCDS <- ABCDS[ABCDS$UNITS %in% c("Normalised Mean (pg/ml)"),]
ABCDS <- ABCDS[!is.na(ABCDS$SUBJECT_LABEL),]
colnames(ABCDS)[6] <- c( "pTau217")
ABCDS <- merge(ABCDS, visit_latency[, c("subject_label", "event_sequence", "age_at_visit")],
               by.x = c("SUBJECT_LABEL", "EVENT_SEQUENCE"), by.y = c("subject_label", "event_sequence"), all = FALSE)
ABCDS$pTau217 <- as.numeric(ABCDS$pTau217)
ABCDS <- merge(ABCDS, demogs, by.x = "SUBJECT_LABEL", by.y = "subject_label", all.x = TRUE, all.y = FALSE)
ABCDS <- merge(ABCDS, ABCDS_apoe[, c("subject_label", "APOE_grouped")], by.x = "SUBJECT_LABEL", by.y = "subject_label", all = FALSE)
ABCDS <- merge(ABCDS, mcrt[, c("subject_label", "event_sequence", "trs")], by.x = c("SUBJECT_LABEL", "EVENT_SEQUENCE"),
               by.y = c("subject_label", "event_sequence"), all.x = TRUE, all.y = FALSE)

# DABNI <- readxl::read_xlsx("./Data_202512/Biofluids_20251125.xlsx")
# DABNI <- DABNI[, c("NHC", "APOE", "Sex", "Age_at_CSF", "PLASMA_PTAU217_S0217")]

pTau217_thresh_DABNI <- mean(DABNI[DABNI$Age_at_CSF <= 35 & !duplicated(DABNI$NHC),]$PLASMA_PTAU217_S0217, na.rm = TRUE) + 
  1.96 *sd(DABNI[DABNI$Age_at_CSF <= 35 & !duplicated(DABNI$NHC),]$PLASMA_PTAU217_S0217, na.rm = TRUE) 

pTau217_thresh_ABCDS <- mean(ABCDS[ABCDS$age_at_visit <= 35 & !duplicated(ABCDS$SUBJECT_LABEL),]$pTau217, na.rm = TRUE) + 
  1.96 *sd(ABCDS[ABCDS$age_at_visit <= 35 & !duplicated(ABCDS$SUBJECT_LABEL),]$pTau217, na.rm = TRUE) 


# DABNI_CL <- readxl::read_xlsx("./Data_202512/CL_20251127.xlsx")

CL_thresh_DABNI <- mean(DABNI_CL[DABNI_CL$Age_at_PET <= 35 & !duplicated(DABNI_CL$NHC),]$CL_SPM12, na.rm = TRUE) + 
  1.96 *sd(DABNI_CL[DABNI_CL$Age_at_PET <= 35 & !duplicated(DABNI_CL$NHC),]$CL_SPM12, na.rm = TRUE) 

ABCDS_CL <- read.csv(".././ABCDS_DF3_AmyloidOnly/all_amy_120924.csv")
ABCDS_CL <- ABCDS_CL[ABCDS_CL$PET.TC.QC.Status == "Passed", c("Subject", "VISIT", "WUSTLcentiloid")]
ABCDS_CL$event_sequence <- as.numeric(substr(ABCDS_CL$VISIT, start = nchar(ABCDS_CL$VISIT), stop = nchar(ABCDS_CL$VISIT)))
ABCDS_CL <- merge(ABCDS_CL, visit_latency[, c("subject_label", "event_sequence", "age_at_visit")],
                  by.x = c("Subject", "event_sequence"), by.y = c("subject_label", "event_sequence"), all = FALSE)
ABCDS_CL <- merge(ABCDS_CL, demogs, by.x = "Subject", by.y = "subject_label", all.x = TRUE, all.y = FALSE)
ABCDS_CL <- merge(ABCDS_CL, ABCDS_apoe, by.x = "Subject", by.y = "subject_label", all = FALSE)

CL_thresh_ABCDS <- mean(ABCDS_CL[ABCDS_CL$age_at_visit <= 35 & !duplicated(ABCDS_CL$Subject),]$WUSTLcentiloid, na.rm = TRUE) + 
  1.96 *sd(ABCDS_CL[ABCDS_CL$age_at_visit <= 35 & !duplicated(ABCDS_CL$Subject),]$WUSTLcentiloid, na.rm = TRUE) 

# DABNI_Hipp <- readxl::read_xlsx("./Data_202512/Hipp_Volume_T1_20251127.xlsx")
# DABNI_Hipp <- DABNI_Hipp[, c("NHC", "Sex", "Age_at_MRI", "FS_Left-Hippocampus", "FS_Right-Hippocampus", "FS_EstimatedTotalIntraCranialVol")]


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
hip_mod_abcds <- lm(hippo ~ EstimatedTotalIntraCranialVol, data = ABCDS_Hipp)

ABCDS_Hipp$hippo_ICVresid <- resid(hip_mod_abcds)
visit_latency$event_sequence <- as.factor(visit_latency$event_sequence)
ABCDS_Hipp <- merge(ABCDS_Hipp, visit_latency[, c("subject_label", "event_sequence", "age_at_visit")],
                    by.x = c("Subject.ID", "event_sequence"), by.y = c("subject_label", "event_sequence"), all = FALSE)
ABCDS_Hipp <- merge(ABCDS_Hipp, demogs, by.x = "Subject.ID", by.y = "subject_label", all.x = TRUE, all.y = FALSE)
ABCDS_Hipp <- merge(ABCDS_Hipp, ABCDS_apoe[, c("subject_label", "APOE_grouped")], by.x = "Subject.ID", by.y = "subject_label", all = FALSE)
ABCDS_Hipp$hippo_Z <- (ABCDS_Hipp$hippo_ICVresid - mean(ABCDS_Hipp[ABCDS_Hipp$age_at_visit <= 35,]$hippo_ICVresid, na.rm = TRUE)) /
  sd(ABCDS_Hipp[ABCDS_Hipp$age_at_visit <= 35,]$hippo_ICVresid, na.rm = TRUE)

hippo_thresh_ABCDS <- mean(ABCDS_Hipp[ABCDS_Hipp$age_at_visit <= 35 & !duplicated(ABCDS_Hipp$Subject.ID),]$hippo_Z, na.rm = TRUE) - 
  1.96 * sd(ABCDS_Hipp[ABCDS_Hipp$age_at_visit <= 35 & !duplicated(ABCDS_Hipp$Subject.ID),]$hippo_Z, na.rm = TRUE)


#Cleaning up hippocampal volume
DABNI_Hipp$hippo <- DABNI_Hipp$`FS_Left-Hippocampus` + DABNI_Hipp$`FS_Right-Hippocampus`
hip_mod <- lm(hippo ~ FS_EstimatedTotalIntraCranialVol, data = DABNI_Hipp)
DABNI_Hipp$hippo_ICVresid <- resid(hip_mod)
DABNI_Hipp$hippo_Z <- (DABNI_Hipp$hippo_ICVresid - mean(DABNI_Hipp[DABNI_Hipp$Age_at_MRI <= 35,]$hippo_ICVresid, na.rm = TRUE)) /
  sd(DABNI_Hipp[DABNI_Hipp$Age_at_MRI <= 35,]$hippo_ICVresid, na.rm = TRUE)


ggplot() + geom_histogram(data = ABCDS_Hipp, aes(x = hippo_Z), fill = "red", alpha = 0.4) + 
  geom_histogram(data = DABNI_Hipp, aes(x = hippo_Z), fill = "blue", alpha = 0.4) + 
  theme_bw() 

hist(ABCDS_Hipp$hippo_Z)
hist(DABNI_Hipp$hippo_Z)

hippo_thresh_DABNI <- mean(DABNI_Hipp[DABNI_Hipp$Age_at_MRI <= 35 & !duplicated(DABNI_Hipp$NHC),]$hippo_Z, na.rm = TRUE) - 
  1.96 * sd(DABNI_Hipp[DABNI_Hipp$Age_at_MRI <= 35 & !duplicated(DABNI_Hipp$NHC),]$hippo_Z, na.rm = TRUE)


DABNI$pTau217_positive <- ifelse(DABNI$PLASMA_PTAU217_S0217 > pTau217_thresh_DABNI, 1, 0)
DABNI_CL$CL_positive <- ifelse(DABNI_CL$CL_SPM12 > CL_thresh_DABNI, 1, 0)
DABNI_Hipp$hippo_positive <- ifelse(DABNI_Hipp$hippo_Z < hippo_thresh_DABNI, 1, 0)


DABNI$APOE_grouped <- ifelse(DABNI$APOE %in% c("22", "23"), "APOE2",
                             ifelse(DABNI$APOE %in% c("33"), "APOE3",
                                    ifelse(DABNI$APOE %in% c("24", "34", "44"), "APOE4", NA)))
DABNI_CL$APOE_grouped <- ifelse(DABNI_CL$APOE %in% c("22", "23"), "APOE2",
                                ifelse(DABNI_CL$APOE %in% c("33"), "APOE3",
                                       ifelse(DABNI_CL$APOE %in% c("24", "34", "44"), "APOE4", NA)))
# dabni_apoe <- rbind(DABNI[, c("NHC", "APOE_grouped")],
#                     DABNI_CL[, c("NHC", "APOE_grouped")])
# dabni_apoe <- dabni_apoe[!duplicated(dabni_apoe),]
# DABNI_Hipp <- merge(DABNI_Hipp, dabni_apoe, by = "NHC", all = FALSE)
# DABNI_dx <- merge(DABNI_dx, dabni_apoe, by = "NHC", all = FALSE)
DABNI_Hipp$APOE_grouped <- ifelse(DABNI_Hipp$APOE %in% c("22", "23"), "APOE2",
                                  ifelse(DABNI_Hipp$APOE %in% c("33"), "APOE3",
                                         ifelse(DABNI_Hipp$APOE %in% c("24", "34", "44"), "APOE4", NA)))
DABNI_dx$APOE_grouped <- ifelse(DABNI_dx$APOE %in% c("22", "23"), "APOE2",
                                ifelse(DABNI_dx$APOE %in% c("33"), "APOE3",
                                       ifelse(DABNI_dx$APOE %in% c("24", "34", "44"), "APOE4", NA)))



ABCDS$pTau217_positive <- ifelse(ABCDS$pTau217 > pTau217_thresh_ABCDS, 1, 0)
ABCDS_CL$CL_positive <- ifelse(ABCDS_CL$WUSTLcentiloid > CL_thresh_ABCDS, 1, 0)
ABCDS_Hipp$hippo_positive <- ifelse(ABCDS_Hipp$hippo_Z < hippo_thresh_ABCDS, 1, 0)

DABNI_dx$mcrt_dementia <- ifelse(DABNI_dx$crt_rli_rfi_total < 23, 1, 0)
DABNI_dx$mcrt_MCI <- ifelse(DABNI_dx$crt_rli_rfi_total < 28 & DABNI_dx$crt_rli_rfi_total > 22, 1, 0)
DABNI_dx$mcrt_impaired <- ifelse(DABNI_dx$crt_rli_rfi_total < 28, 1, 0)

ABCDS$mcrt_dementia <- ifelse(ABCDS$trs < 23, 1, 0)
ABCDS$mcrt_MCI <- ifelse(ABCDS$trs < 28 & ABCDS$trs > 22, 1, 0)
ABCDS$mcrt_impaired <- ifelse(ABCDS$trs < 28, 1, 0)

DABNI_CSF <- DABNI[!is.na(DABNI$CSF_AB142_A0010),]
DABNI_CSF <- DABNI_CSF %>% rename(ID = NHC, Age = Age_at_CSF, APOE_grouped = APOE_grouped,
                                  SEX = Sex)
DABNI_CSF$SEX <- as.factor(DABNI_CSF$SEX)
levels(DABNI_CSF$SEX) <- c("FEMALE", "MALE")
DABNI_CSF$Group <- paste0(DABNI_CSF$SEX, DABNI_CSF$APOE_grouped)

results_CSF_grouped_bs <- sliding_proportion_bs(DABNI_CSF[!duplicated(DABNI_CSF$ID) & !is.na(DABNI_CSF$APOE_grouped),], age_col = "Age", group_col = "Group",
                                                outcome_col = "pTau217_positive", window_size = 10)

results_CSF_grouped_bs <- results_CSF_grouped_bs %>%
  group_by(Group) %>%
  mutate(
    ci_lower_smooth = predict(
      loess(ci_lower ~ window_center,
            span = 0.2,
            na.action = na.exclude),
      newdata = data.frame(window_center = window_center)
    ),
    
    ci_upper_smooth = predict(
      loess(ci_upper ~ window_center,
            span = 0.2,
            na.action = na.exclude),
      newdata = data.frame(window_center = window_center)
    )
  ) %>%
  ungroup()

get_hist_of_counts(results_CSF_grouped_bs[results_CSF_grouped_bs$Group %in% c("FEMALEAPOE2", "FEMALEAPOE4",
                                                                                    "MALEAPOE2", "MALEAPOE4"),]) + ggtitle("A. CSF AB42/AB40") 