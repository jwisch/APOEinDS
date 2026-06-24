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
mcrt_thresh <- 28
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
DABNI_dx$mcrt_MCI <- ifelse(DABNI_dx$crt_rli_rfi_total < mcrt_thresh & DABNI_dx$crt_rli_rfi_total > 22, 1, 0)
DABNI_dx$mcrt_impaired <- ifelse(DABNI_dx$crt_rli_rfi_total < mcrt_thresh, 1, 0)

ABCDS$mcrt_dementia <- ifelse(ABCDS$trs < 23, 1, 0)
ABCDS$mcrt_MCI <- ifelse(ABCDS$trs < mcrt_thresh & ABCDS$trs > 22, 1, 0)
ABCDS$mcrt_impaired <- ifelse(ABCDS$trs < mcrt_thresh, 1, 0)

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

DABNI_CL <- DABNI_CL %>% rename(ID = NHC, Age = Age_at_PET, APOE_grouped = APOE_grouped,
                                SEX = Sex, Apos = CL_positive)
ABCDS_CL <- ABCDS_CL %>% rename(ID = Subject, Age = age_at_visit, APOE_grouped = APOE_grouped,
                                SEX = de_gender, Apos = CL_positive)


CL <- rbind(DABNI_CL[, c("ID", "Age", "APOE_grouped", "SEX", "Apos")],
            ABCDS_CL[, c("ID", "Age", "APOE_grouped", "SEX", "Apos")])
CL$SEX <- as.factor(CL$SEX)
levels(CL$SEX) <- c("MALE", "FEMALE", "FEMALE", "MALE")

CL$Group <- paste0(CL$SEX, CL$APOE_grouped)

results_CL <- sliding_proportion(CL[!duplicated(CL$ID),], age_col = "Age", group_col = "APOE_grouped",
                                 outcome_col = "Apos", window_size = 10)

results_CL_bs <- sliding_proportion_bs(CL[!duplicated(CL$ID),], age_col = "Age", group_col = "APOE_grouped",
                                 outcome_col = "Apos", window_size = 10)

results_CL_grouped_bs <- sliding_proportion_bs(CL[!duplicated(CL$ID) & !is.na(CL$APOE_grouped),], age_col = "Age", group_col = "Group",
                                       outcome_col = "Apos", window_size = 10)

results_CL_bs <- results_CL_bs[!is.na(results_CL_bs$APOE_grouped),] %>%
  group_by(APOE_grouped) %>%
  mutate(
    ci_lower_smooth = predict(loess(ci_lower ~ window_center, span = 0.2)),
    ci_upper_smooth = predict(loess(ci_upper ~ window_center, span = 0.2))
  ) %>%
  ungroup()

results_CL_grouped_bs <- results_CL_grouped_bs %>%
  group_by(Group) %>%
  mutate(
    ci_lower_smooth = predict(loess(ci_lower ~ window_center, span = 0.2)),
    ci_upper_smooth = predict(loess(ci_upper ~ window_center, span = 0.2))
  ) %>%
  ungroup()



# Grouped application
df_age50_ci <- results_CL_bs %>%
  group_by(APOE_grouped) %>%
  summarise(estimate_age50_ci(cur_data_all()), .groups = "drop")

df_age50_ci_grouped <- results_CL_grouped_bs %>%
  group_by(Group) %>%
  summarise(estimate_age50_ci(cur_data_all()), .groups = "drop")

# Step 1: Get APOE4 crossing age
apoe4 <- results_CL_bs %>% filter(APOE_grouped == "APOE4") %>% arrange(window_center)
age_apoe4_50 <- interp_cross(apoe4$window_center, apoe4$prop_pos, cutoff = 0.5)
get_prop_at_age(results_CL_bs[results_CL_bs$APOE_grouped == "APOE2",], age_apoe4_50)
get_prop_at_age(results_CL_bs[results_CL_bs$APOE_grouped == "APOE3",], age_apoe4_50)


df_age50_ci


# Step 1: Get APOE4 crossing age
apoe4 <- results_CL_grouped_bs %>% filter(Group == "FEMALEAPOE4") %>% arrange(window_center)
age_apoe4_50 <- interp_cross(apoe4$window_center, apoe4$prop_pos, cutoff = 0.5)
get_prop_at_age(results_CL_grouped_bs[results_CL_grouped_bs$Group == "MALEAPOE2",], age_apoe4_50)
get_prop_at_age(results_CL_grouped_bs[results_CL_grouped_bs$Group == "MALEAPOE4",], age_apoe4_50)



DABNI_pTau217 <- DABNI %>% rename(ID = NHC, Age = Age_at_CSF, APOE_grouped = APOE_grouped,
                                SEX = Sex, Tpos = pTau217_positive)
ABCDS_pTau217 <- ABCDS %>% rename(ID = SUBJECT_LABEL, Age = age_at_visit, APOE_grouped = APOE_grouped,
                                SEX = de_gender, Tpos = pTau217_positive)


pTau <- rbind(DABNI_pTau217[, c("ID", "Age", "APOE_grouped", "SEX", "Tpos")],
            ABCDS_pTau217[, c("ID", "Age", "APOE_grouped", "SEX", "Tpos")])
pTau$SEX <- as.factor(pTau$SEX)
levels(pTau$SEX) <- c("MALE", "FEMALE", "FEMALE", "MALE")
pTau$Group <- paste0(pTau$SEX, pTau$APOE_grouped)

results_pTau <- sliding_proportion(pTau[!duplicated(pTau$ID),], age_col = "Age", group_col = "APOE_grouped",
                                 outcome_col = "Tpos", window_size = 10)
results_pTau_bs <- sliding_proportion_bs(pTau[!duplicated(pTau$ID),], age_col = "Age", group_col = "APOE_grouped",
                                   outcome_col = "Tpos", window_size = 10)

results_pTau_grouped_bs <- sliding_proportion_bs(pTau[!duplicated(pTau$ID),], age_col = "Age", group_col = "Group",
                                               outcome_col = "Tpos", window_size = 10)

ABCDS_pTau217$SEX <- as.factor(ABCDS_pTau217$SEX)
levels(ABCDS_pTau217$SEX) <- c("MALE", "FEMALE")
ABCDS_pTau217$Group <- paste0(ABCDS_pTau217$SEX, ABCDS_pTau217$APOE_grouped)
results_pTau_ABCDS_bs <- sliding_proportion_bs(ABCDS_pTau217[!duplicated(ABCDS_pTau217$ID),], age_col = "Age", group_col = "Group",
                                                 outcome_col = "Tpos", window_size = 10)
DABNI_pTau217$SEX <- as.factor(DABNI_pTau217$SEX)
levels(DABNI_pTau217$SEX) <- c("FEMALE", "MALE")
DABNI_pTau217$Group <- paste0(DABNI_pTau217$SEX, DABNI_pTau217$APOE_grouped)
results_pTau_DABNI_bs <- sliding_proportion_bs(DABNI_pTau217[!duplicated(DABNI_pTau217$ID),], age_col = "Age", group_col = "Group",
                                               outcome_col = "Tpos", window_size = 10)

results_pTau_bs <- results_pTau_bs %>%
  group_by(APOE_grouped) %>%
  mutate(
    ci_lower_smooth = predict(loess(ci_lower ~ window_center, span = 0.2)),
    ci_upper_smooth = predict(loess(ci_upper ~ window_center, span = 0.2))
  ) %>%
  ungroup()

results_pTau_grouped_bs <- results_pTau_grouped_bs %>%
  group_by(Group) %>%
  mutate(
    ci_lower_smooth = predict(loess(ci_lower ~ window_center, span = 0.2),
                              newdata = data.frame(window_center = window_center)),
    ci_upper_smooth = predict(loess(ci_upper ~ window_center, span = 0.2),
                              newdata = data.frame(window_center = window_center))
  ) %>%
  ungroup()


results_pTau_bs %>%
  group_by(APOE_grouped) %>%
  summarise(estimate_age50_ci(cur_data_all()), .groups = "drop")


results_pTau_ABCDS_bs <- results_pTau_ABCDS_bs %>%
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

results_pTau_DABNI_bs <- results_pTau_DABNI_bs %>%
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

#Need to do GFAP & NFL######################################################
DABNI_GFAP <- DABNI[, c("NHC", "APOE_grouped", "Sex", "Age_at_CSF", "PLASMA_GFAP_S0122")]
ABCDS_GFAP <- df_ABCDS[df_ABCDS$TESTNAME == "GFAP", c("SUBJECT_LABEL", "EVENT_SEQUENCE", "TESTVALUE")]
ABCDS_GFAP <- merge(ABCDS_GFAP, ABCDS_apoe[, c("subject_label", "APOE_grouped")],
                    by.x = "SUBJECT_LABEL", by.y = "subject_label", all = FALSE)
ABCDS_GFAP <- merge(ABCDS_GFAP, demogs, by.x = "SUBJECT_LABEL", by.y = "subject_label", all = FALSE)
ABCDS_GFAP$de_gender <- as.factor(ABCDS_GFAP$de_gender)
levels(ABCDS_GFAP$de_gender) <- c("MALE", "FEMALE")
ABCDS_GFAP <- ABCDS_GFAP[!is.na(ABCDS_GFAP$de_gender),]
ABCDS_GFAP <- ABCDS_GFAP[!duplicated(ABCDS_GFAP$SUBJECT_LABEL),]
ABCDS_GFAP <- merge(ABCDS_GFAP, visit_latency[, c("subject_label", "event_sequence", "age_at_visit")],
                    by.x = c("SUBJECT_LABEL", "EVENT_SEQUENCE"), by.y = c("subject_label", "event_sequence"), all = FALSE)
ABCDS_GFAP <- ABCDS_GFAP %>% rename(ID = SUBJECT_LABEL, Age = age_at_visit, APOE_grouped = APOE_grouped,
                                    SEX = de_gender)
ABCDS_GFAP$Pos_dummy <- 1
DABNI_GFAP <- DABNI_GFAP %>% rename(ID = NHC, Age = Age_at_CSF, APOE_grouped = APOE_grouped,
                                    SEX = Sex)
DABNI_GFAP$SEX <- as.factor(DABNI_GFAP$SEX)
levels(DABNI_GFAP$SEX) <- c("FEMALE", "MALE")
DABNI_GFAP$Pos_dummy <- 0

GFAP <- rbind(DABNI_GFAP[, c("ID", "Age", "APOE_grouped", "SEX", "Pos_dummy")], 
              ABCDS_GFAP[, c("ID", "Age", "APOE_grouped", "SEX", "Pos_dummy")])
GFAP$Group <- paste0(GFAP$SEX, GFAP$APOE_grouped)
results_GFAP_grouped_bs <- sliding_proportion_bs(GFAP[!duplicated(GFAP$ID),], age_col = "Age", group_col = "Group",
                                                 outcome_col = "Pos_dummy", window_size = 10)
results_GFAP_grouped_bs <- results_GFAP_grouped_bs %>%
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



DABNI_NFL <- DABNI[, c("NHC", "APOE_grouped", "Sex", "Age_at_CSF", "PLASMA_NFLIGHT_S0100")]
ABCDS_NFL <- df_ABCDS[df_ABCDS$TESTNAME == "NF-light", c("SUBJECT_LABEL", "EVENT_SEQUENCE", "TESTVALUE")]
ABCDS_NFL <- merge(ABCDS_NFL, ABCDS_apoe[, c("subject_label", "APOE_grouped")],
                    by.x = "SUBJECT_LABEL", by.y = "subject_label", all = FALSE)
ABCDS_NFL <- merge(ABCDS_NFL, demogs, by.x = "SUBJECT_LABEL", by.y = "subject_label", all = FALSE)
ABCDS_NFL$de_gender <- as.factor(ABCDS_NFL$de_gender)
levels(ABCDS_NFL$de_gender) <- c("MALE", "FEMALE")
ABCDS_NFL <- ABCDS_NFL[!is.na(ABCDS_NFL$de_gender),]
ABCDS_NFL <- ABCDS_NFL[!duplicated(ABCDS_NFL$SUBJECT_LABEL),]
ABCDS_NFL <- merge(ABCDS_NFL, visit_latency[, c("subject_label", "event_sequence", "age_at_visit")],
                    by.x = c("SUBJECT_LABEL", "EVENT_SEQUENCE"), by.y = c("subject_label", "event_sequence"), all = FALSE)
ABCDS_NFL <- ABCDS_NFL %>% rename(ID = SUBJECT_LABEL, Age = age_at_visit, APOE_grouped = APOE_grouped,
                                    SEX = de_gender)
ABCDS_NFL$Pos_dummy <- 1
DABNI_NFL <- DABNI_NFL %>% rename(ID = NHC, Age = Age_at_CSF, APOE_grouped = APOE_grouped,
                                    SEX = Sex)
DABNI_NFL$SEX <- as.factor(DABNI_NFL$SEX)
levels(DABNI_NFL$SEX) <- c("FEMALE", "MALE")
DABNI_NFL$Pos_dummy <- 0

NFL <- rbind(DABNI_NFL[, c("ID", "Age", "APOE_grouped", "SEX", "Pos_dummy")], 
              ABCDS_NFL[, c("ID", "Age", "APOE_grouped", "SEX", "Pos_dummy")])
NFL$Group <- paste0(NFL$SEX, NFL$APOE_grouped)
results_NFL_grouped_bs <- sliding_proportion_bs(NFL[!duplicated(NFL$ID),], age_col = "Age", group_col = "Group",
                                                 outcome_col = "Pos_dummy", window_size = 10)
results_NFL_grouped_bs <- results_NFL_grouped_bs %>%
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
###############################################################################

apoe4 <- results_pTau_bs %>% filter(APOE_grouped == "APOE4") %>% arrange(window_center)
age_apoe4_50 <- interp_cross(apoe4$window_center, apoe4$prop_pos, cutoff = 0.5)
get_prop_at_age(results_pTau_bs[results_pTau_bs$APOE_grouped == "APOE2",], age_apoe4_50)
get_prop_at_age(results_pTau_bs[results_pTau_bs$APOE_grouped == "APOE3",], age_apoe4_50)
interp_cross(apoe4$window_center, apoe4$ci_lower_smooth, cutoff = 0.5)
interp_cross(apoe4$window_center, apoe4$ci_upper_smooth, cutoff = 0.5)

DABNI_Hipp <- DABNI_Hipp %>% rename(ID = NHC, Age = Age_at_MRI, APOE_grouped = APOE_grouped,
                                  SEX = Sex, Npos = hippo_positive)
ABCDS_Hipp <- ABCDS_Hipp %>% rename(ID = Subject.ID, Age = age_at_visit, APOE_grouped = APOE_grouped,
                                  SEX = de_gender, Npos = hippo_positive)
Hippo <- rbind(DABNI_Hipp[, c("ID", "Age", "APOE_grouped", "SEX", "Npos")],
              ABCDS_Hipp[, c("ID", "Age", "APOE_grouped", "SEX", "Npos")])
Hippo$SEX <- as.factor(Hippo$SEX)
levels(Hippo$SEX) <- c("MALE", "FEMALE", "FEMALE", "MALE")

Hippo$Group <- paste0(Hippo$SEX, Hippo$APOE_grouped)

results_Hippo <- sliding_proportion(Hippo[!duplicated(Hippo$ID),], age_col = "Age", group_col = "APOE_grouped",
                                   outcome_col = "Npos", window_size = 10)
results_Hippo_bs <- sliding_proportion_bs(Hippo[!duplicated(Hippo$ID),], age_col = "Age", group_col = "APOE_grouped",
                                    outcome_col = "Npos", window_size = 10)

results_Hippo_grouped_bs <- sliding_proportion_bs(Hippo[!duplicated(Hippo$ID),], age_col = "Age", group_col = "Group",
                                          outcome_col = "Npos", window_size = 10)
results_Hippo_bs <- results_Hippo_bs %>%
  group_by(APOE_grouped) %>%
  mutate(
    ci_lower_smooth = predict(loess(ci_lower ~ window_center, span = 0.2)),
    ci_upper_smooth = predict(loess(ci_upper ~ window_center, span = 0.2))
  ) %>%
  ungroup()

results_Hippo_grouped_bs <- results_Hippo_grouped_bs %>%
  group_by(Group) %>%
  mutate(
    ci_lower_smooth = predict(loess(ci_lower ~ window_center, span = 0.2)),
    ci_upper_smooth = predict(loess(ci_upper ~ window_center, span = 0.2))
  ) %>%
  ungroup()



results_Hippo_bs %>%
  group_by(APOE_grouped) %>%
  summarise(estimate_age50_ci(cur_data_all()), .groups = "drop")

apoe4 <- results_Hippo_bs %>% filter(APOE_grouped == "APOE4") %>% arrange(window_center)
age_apoe4_50 <- interp_cross(apoe4$window_center, apoe4$prop_pos, cutoff = 0.5)
get_prop_at_age(results_Hippo_bs[results_Hippo_bs$APOE_grouped == "APOE2",], age_apoe4_50)
get_prop_at_age(results_Hippo_bs[results_Hippo_bs$APOE_grouped == "APOE3",], age_apoe4_50)

DABNI$Impaired <- ifelse(DABNI$Diag %in% c("aDS"), "asymptomatic",
                            ifelse(DABNI$Diag %in% c("dDS", "pDS"), "symptomatic", "other"))

DABNI_dx <- DABNI_dx %>% rename(ID = NHC, Age = Age_NPS, APOE_grouped = APOE_grouped,
                                    SEX = Sex)
DABNI <- DABNI %>% rename(ID = NHC, Age = Age_at_CSF, APOE_grouped = APOE_grouped,
                                SEX = Sex)

dx$Impaired <- ifelse(dx$consensus_dx == 0, "asymptomatic",
                      ifelse(dx$consensus_dx %in% c(1, 2), "symptomatic", "other"))
dx <- merge(dx, visit_latency[, c("subject_label", "event_sequence", "age_at_visit")], 
            by = c("subject_label", "event_sequence"), all = FALSE)

dx <- dx %>% rename(ID = subject_label, Age = age_at_visit, APOE_grouped = APOE_grouped,
                                    SEX = de_gender)
dx <- rbind(dx[, c("ID", "Age", "APOE_grouped", "SEX", "Impaired")],
               DABNI[, c("ID", "Age", "APOE_grouped", "SEX", "Impaired")])
dx$SEX <- as.factor(dx$SEX)
levels(dx$SEX) <- c("MALE", "FEMALE", "FEMALE", "MALE")
dx$Group <- paste0(dx$SEX, dx$APOE_grouped)

dx$Impaired <- as.factor(dx$Impaired)
levels(dx$Impaired) <- c(0, NA, 1)
dx <- dx[!is.na(dx$Impaired),]
results_dx <- sliding_proportion(dx[!duplicated(dx$ID),], age_col = "Age", group_col = "APOE_grouped",
                                    outcome_col = "Impaired", window_size = 10)
results_dx_bs <- sliding_proportion_bs(dx[!duplicated(dx$ID),], age_col = "Age", group_col = "APOE_grouped",
                                 outcome_col = "Impaired", window_size = 10)
results_dx_grouped_bs <- sliding_proportion_bs(dx[!duplicated(dx$ID),], age_col = "Age", group_col = "Group",
                                       outcome_col = "Impaired", window_size = 10)

results_dx_bs <- results_dx_bs %>%
  group_by(APOE_grouped) %>%
  mutate(
    ci_lower_smooth = predict(loess(ci_lower ~ window_center, span = 0.5)),
    ci_upper_smooth = predict(loess(ci_upper ~ window_center, span = 0.5))
  ) %>%
  ungroup()

results_dx_grouped_bs <- results_dx_grouped_bs %>%
  group_by(Group) %>%
  mutate(
    ci_lower_smooth = predict(loess(ci_lower ~ window_center, span = 0.5)),
    ci_upper_smooth = predict(loess(ci_upper ~ window_center, span = 0.5))
  ) %>%
  ungroup()

results_dx_bs %>%
  group_by(APOE_grouped) %>%
  summarise(estimate_age50_ci(cur_data_all()), .groups = "drop")


apoe4 <- results_dx_bs %>% filter(APOE_grouped == "APOE4") %>% arrange(window_center)
age_apoe4_50 <- interp_cross(apoe4$window_center, apoe4$prop_pos, cutoff = 0.5)
get_prop_at_age(results_dx_bs[results_dx_bs$APOE_grouped == "APOE2",], age_apoe4_50)
get_prop_at_age(results_dx_bs[results_dx_bs$APOE_grouped == "APOE3",], age_apoe4_50)

##########MCRT
ABCDS <- ABCDS %>% rename(ID = SUBJECT_LABEL, Age = age_at_visit, APOE_grouped = APOE_grouped,
                    SEX = de_gender)
mcrt <- rbind(ABCDS[, c("ID", "SEX", "Age", "APOE_grouped", "mcrt_dementia", "mcrt_MCI", "mcrt_impaired")],
              DABNI_dx[, c("ID", "SEX", "Age", "APOE_grouped", "mcrt_dementia", "mcrt_MCI", "mcrt_impaired")])
mcrt <- rbind(DABNI_dx[, c("ID", "SEX", "Age", "APOE_grouped", "mcrt_dementia", "mcrt_MCI", "mcrt_impaired")])
mcrt$SEX <- as.factor(mcrt$SEX)
levels(mcrt$SEX) <- c("MALE", "FEMALE", "FEMALE", "MALE")
mcrt$Group <- paste0(mcrt$SEX, mcrt$APOE_grouped)


results_mcrt <- sliding_proportion(mcrt[!duplicated(mcrt$ID) & !is.na(mcrt$APOE_grouped),], age_col = "Age", group_col = "APOE_grouped",
                                 outcome_col = "mcrt_dementia", window_size = 10)
results_mcrt_bs <- sliding_proportion_bs(mcrt[!duplicated(mcrt$ID) & !is.na(mcrt$APOE_grouped),], age_col = "Age", group_col = "APOE_grouped",
                                       outcome_col = "mcrt_dementia", window_size = 10)
results_mcrt_grouped_bs <- sliding_proportion_bs(mcrt[!duplicated(mcrt$ID) & !is.na(mcrt$APOE_grouped),], age_col = "Age", group_col = "Group",
                                               outcome_col = "mcrt_dementia", window_size = 10)

results_mcrt_bs <- results_mcrt_bs %>%
  group_by(APOE_grouped) %>%
  mutate(
    ci_lower_smooth = predict(
      loess(ci_lower ~ window_center, span = 0.5),
      newdata = tibble(window_center = window_center)
    ),
    ci_upper_smooth = predict(
      loess(ci_upper ~ window_center, span = 0.5),
      newdata = tibble(window_center = window_center)
    )
  ) %>%
  ungroup()

results_mcrt_grouped_bs <- results_mcrt_grouped_bs %>%
  group_by(Group) %>%
  mutate(
    ci_lower_smooth = predict(
      loess(ci_lower ~ window_center, span = 0.5),
      newdata = tibble(window_center = window_center)
    ),
    ci_upper_smooth = predict(
      loess(ci_upper ~ window_center, span = 0.5),
      newdata = tibble(window_center = window_center)
    )
  ) %>%
  ungroup()


apoe4 <- results_mcrt_bs %>% filter(APOE_grouped == "APOE4") %>% arrange(window_center)
age_apoe4_50 <- interp_cross(apoe4$window_center, apoe4$prop_pos, cutoff = 0.5)
get_prop_at_age(results_mcrt_bs[results_mcrt_bs$APOE_grouped == "APOE2",], age_apoe4_50)
get_prop_at_age(results_mcrt_bs[results_mcrt_bs$APOE_grouped == "APOE3",], age_apoe4_50)


results_mcrt_MCI <- sliding_proportion(mcrt[!duplicated(mcrt$ID),], age_col = "Age", group_col = "APOE_grouped",
                                   outcome_col = "mcrt_MCI", window_size = 10)
results_mcrt_MCI_bs <- sliding_proportion_bs(mcrt[!duplicated(mcrt$ID),], age_col = "Age", group_col = "APOE_grouped",
                                         outcome_col = "mcrt_MCI", window_size = 10)
results_mcrt_MCI_grouped_bs <- sliding_proportion_bs(mcrt[!duplicated(mcrt$ID),], age_col = "Age", group_col = "Group",
                                                 outcome_col = "mcrt_MCI", window_size = 10)

results_mcrt_MCI_bs <- results_mcrt_MCI_bs %>%
  group_by(APOE_grouped) %>%
  mutate(
    ci_lower_smooth = predict(
      loess(ci_lower ~ window_center, span = 0.5),
      newdata = tibble(window_center = window_center)
    ),
    ci_upper_smooth = predict(
      loess(ci_upper ~ window_center, span = 0.5),
      newdata = tibble(window_center = window_center)
    )
  ) %>%
  ungroup()

results_mcrt_MCI_grouped_bs <- results_mcrt_MCI_grouped_bs %>%
  group_by(Group) %>%
  mutate(
    ci_lower_smooth = predict(
      loess(ci_lower ~ window_center, span = 0.5),
      newdata = tibble(window_center = window_center)
    ),
    ci_upper_smooth = predict(
      loess(ci_upper ~ window_center, span = 0.5),
      newdata = tibble(window_center = window_center)
    )
  ) %>%
  ungroup()



results_mcrt_MCI_bs %>%
  group_by(APOE_grouped) %>%
  summarise(estimate_age50_ci(cur_data_all()), .groups = "drop")


apoe4 <- results_mcrt_MCI_bs %>% filter(APOE_grouped == "APOE4") %>% arrange(window_center)
age_apoe4_50 <- interp_cross(apoe4$window_center, apoe4$prop_pos, cutoff = 0.5)
get_prop_at_age(results_mcrt_MCI_bs[results_mcrt_MCI_bs$APOE_grouped == "APOE2",], age_apoe4_50)
get_prop_at_age(results_mcrt_MCI_bs[results_mcrt_MCI_bs$APOE_grouped == "APOE3",], age_apoe4_50)



results_mcrt_impaired <- sliding_proportion(mcrt[!duplicated(mcrt$ID),], age_col = "Age", group_col = "APOE_grouped",
                                       outcome_col = "mcrt_impaired", window_size = 10)
results_mcrt_impaired_bs <- sliding_proportion_bs(mcrt[!duplicated(mcrt$ID),], age_col = "Age", group_col = "APOE_grouped",
                                             outcome_col = "mcrt_impaired", window_size = 10)
results_mcrt_impaired_grouped_bs <- sliding_proportion_bs(mcrt[!duplicated(mcrt$ID),], age_col = "Age", group_col = "Group",
                                                     outcome_col = "mcrt_impaired", window_size = 10)


apoe4 <- results_mcrt_impaired_bs %>% filter(APOE_grouped == "APOE4") %>% arrange(window_center)
age_apoe4_50 <- interp_cross(apoe4$window_center, apoe4$prop_pos, cutoff = 0.5)
get_prop_at_age(results_mcrt_MCI_bs[results_mcrt_MCI_bs$APOE_grouped == "APOE2",], age_apoe4_50)
get_prop_at_age(results_mcrt_MCI_bs[results_mcrt_MCI_bs$APOE_grouped == "APOE3",], age_apoe4_50)

results_mcrt_impaired_bs <- results_mcrt_impaired_bs %>%
  group_by(APOE_grouped) %>%
  mutate(
    ci_lower_smooth = {
      ok <- !is.na(ci_lower) & !is.na(window_center)
      if (sum(ok) > 3) {
        fit <- loess(ci_lower ~ window_center, data = cur_data()[ok, ], span = 0.5)
        preds <- predict(fit, newdata = data.frame(window_center = window_center))
      } else {
        preds <- rep(NA_real_, n())
      }
      preds
    },
    ci_upper_smooth = {
      ok <- !is.na(ci_upper) & !is.na(window_center)
      if (sum(ok) > 3) {
        fit <- loess(ci_upper ~ window_center, data = cur_data()[ok, ], span = 0.5)
        preds <- predict(fit, newdata = data.frame(window_center = window_center))
      } else {
        preds <- rep(NA_real_, n())
      }
      preds
    }
  ) %>%
  ungroup()

results_mcrt_impaired_grouped_bs <- results_mcrt_impaired_grouped_bs %>%
  group_by(Group) %>%
  mutate(
    ci_lower_smooth = {
      ok <- !is.na(ci_lower) & !is.na(window_center)
      if (sum(ok) > 3) {
        fit <- loess(ci_lower ~ window_center, data = cur_data()[ok, ], span = 0.5)
        preds <- predict(fit, newdata = data.frame(window_center = window_center))
      } else {
        preds <- rep(NA_real_, n())
      }
      preds
    },
    ci_upper_smooth = {
      ok <- !is.na(ci_upper) & !is.na(window_center)
      if (sum(ok) > 3) {
        fit <- loess(ci_upper ~ window_center, data = cur_data()[ok, ], span = 0.5)
        preds <- predict(fit, newdata = data.frame(window_center = window_center))
      } else {
        preds <- rep(NA_real_, n())
      }
      preds
    }
  ) %>%
  ungroup()

results_mcrt_impaired_bs %>%
  group_by(APOE_grouped) %>%
  summarise(estimate_age50_ci(cur_data_all()), .groups = "drop")


apoe4 <- results_mcrt_impaired_bs %>% filter(APOE_grouped == "APOE4") %>% arrange(window_center)
age_apoe4_50 <- interp_cross(apoe4$window_center, apoe4$prop_pos, cutoff = 0.5)
get_prop_at_age(results_mcrt_impaired_bs[results_mcrt_impaired_bs$APOE_grouped == "APOE2",], age_apoe4_50)
get_prop_at_age(results_mcrt_impaired_bs[results_mcrt_impaired_bs$APOE_grouped == "APOE3",], age_apoe4_50)


interp_cross(apoe4$window_center, apoe4$ci_upper_smooth, cutoff = 0.5)
interp_cross(apoe4$window_center, apoe4$ci_lower_smooth, cutoff = 0.5)

results_CL_bs$APOE_grouped <- factor(results_CL_bs$APOE_grouped, levels = c("APOE4", "APOE3", "APOE2"))
results_Hippo_bs$APOE_grouped <- factor(results_Hippo_bs$APOE_grouped, levels = c("APOE4", "APOE3", "APOE2"))
results_mcrt_impaired_bs$APOE_grouped <- factor(results_mcrt_impaired_bs$APOE_grouped, levels = c("APOE4", "APOE3", "APOE2"))
results_dx_bs$APOE_grouped <- factor(results_dx_bs$APOE_grouped, levels = c("APOE4", "APOE3", "APOE2"))


p1 <- get_bs_plot(results_CL_bs[!is.na(results_CL_bs$APOE_grouped),])+ ggtitle("A. Amyloid PET")
p2 <- get_bs_plot(results_Hippo_bs[!is.na(results_Hippo_bs$APOE_grouped),])+ ggtitle("B. Hippocampal Volume")
p3 <- get_bs_plot(results_mcrt_impaired_bs[!is.na(results_mcrt_impaired_bs$APOE_grouped),])+ ggtitle(" MCRT Total Score < 28")
p4 <- get_bs_plot(results_dx_bs[!is.na(results_dx_bs$APOE_grouped),])+ ggtitle("D. Cognitive Impairment")

final_plot <- p1 | p2 | p3 | p4
final_plot

graph2ppt(file = "./figures/Prevalence_by_APOE.pptx", width = 10, height = 5.78)

results_CSF_grouped_bs$Group <- factor(results_CSF_grouped_bs$Group, levels = c("FEMALEAPOE4", "FEMALEAPOE3",
                                                                              "FEMALEAPOE2","MALEAPOE4",
                                                                              "MALEAPOE3",
                                                                              "MALEAPOE2"))
results_CL_grouped_bs$Group <- factor(results_CL_grouped_bs$Group, levels = c("FEMALEAPOE4", "FEMALEAPOE3",
                                                                                  "FEMALEAPOE2","MALEAPOE4",
                                                                                  "MALEAPOE3",
                                                                                  "MALEAPOE2"))
results_mcrt_grouped_bs$Group <- factor(results_mcrt_grouped_bs$Group, levels = c("FEMALEAPOE4", "FEMALEAPOE3",
                                                                              "FEMALEAPOE2","MALEAPOE4",
                                                                              "MALEAPOE3",
                                                                              "MALEAPOE2"))
results_Hippo_grouped_bs$Group <- factor(results_Hippo_grouped_bs$Group, levels = c("FEMALEAPOE4", "FEMALEAPOE3",
                                                                              "FEMALEAPOE2","MALEAPOE4",
                                                                              "MALEAPOE3",
                                                                              "MALEAPOE2"))
results_dx_grouped_bs$Group <- factor(results_dx_grouped_bs$Group, levels = c("FEMALEAPOE4", "FEMALEAPOE3",
                                                                              "FEMALEAPOE2","MALEAPOE4",
                                                                              "MALEAPOE3",
                                                                              "MALEAPOE2"))

p1 <- get_bs_plot_sex(results_CL_grouped_bs[results_CL_grouped_bs$Group %in% c("FEMALEAPOE2", "FEMALEAPOE3", "FEMALEAPOE4",
                                                                         "MALEAPOE2", "MALEAPOE3", "MALEAPOE4"),]) + ggtitle("A. Amyloid PET") 
p3 <- get_bs_plot_sex(results_mcrt_grouped_bs[results_mcrt_grouped_bs$Group %in% c("FEMALEAPOE2", "FEMALEAPOE3", "FEMALEAPOE4",
                                                                         "MALEAPOE2", "MALEAPOE3", "MALEAPOE4"),])+ ggtitle("C. MCRT Total Score < 28")

p2 <- get_bs_plot_sex(results_Hippo_grouped_bs[results_Hippo_grouped_bs$Group %in% c("FEMALEAPOE2", "FEMALEAPOE3", "FEMALEAPOE4",
                                                                             "MALEAPOE2", "MALEAPOE3", "MALEAPOE4"),])+ ggtitle("B. Hippocampal Volume")
p4 <- get_bs_plot_sex(results_dx_grouped_bs[results_dx_grouped_bs$Group %in% c("FEMALEAPOE2", "FEMALEAPOE3", "FEMALEAPOE4",
                                                                             "MALEAPOE2", "MALEAPOE3", "MALEAPOE4"),])+ ggtitle("D. Cognitive Impairment")
final_plot_sex <- p1 | p2 | p3 | p4
final_plot_sex

graph2ppt(file = "./figures/Prevalence_by_APOE_SEX.pptx", width = 10, height = 5.78)

p1 <- get_hist_of_counts(results_CSF_grouped_bs[results_CSF_grouped_bs$Group %in% c("FEMALEAPOE2", "FEMALEAPOE3", "FEMALEAPOE4",
                                                                         "MALEAPOE2", "MALEAPOE3", "MALEAPOE4"),]) + ggtitle("A. CSF AB42/AB40") 
p2 <- get_hist_of_counts(results_CL_grouped_bs[results_CL_grouped_bs$Group %in% c("FEMALEAPOE2", "FEMALEAPOE3", "FEMALEAPOE4",
                                                                         "MALEAPOE2", "MALEAPOE3", "MALEAPOE4"),]) + ggtitle("B. Amyloid PET") 
p3 <- get_hist_of_counts(results_pTau_DABNI_bs[results_pTau_DABNI_bs$Group %in% c("FEMALEAPOE2", "FEMALEAPOE3", "FEMALEAPOE4",
                                                                         "MALEAPOE2", "MALEAPOE3", "MALEAPOE4"),]) + ggtitle("C. Plasma pTau217 - AlzPATH") 
p4 <- get_hist_of_counts(results_pTau_ABCDS_bs[results_pTau_ABCDS_bs$Group %in% c("FEMALEAPOE2", "FEMALEAPOE3", "FEMALEAPOE4",
                                                                         "MALEAPOE2", "MALEAPOE3", "MALEAPOE4"),]) + ggtitle("D. Plasma pTau217 - Lilly") 
p5 <- get_hist_of_counts(results_GFAP_grouped_bs[results_GFAP_grouped_bs$Group %in% c("FEMALEAPOE2", "FEMALEAPOE3", "FEMALEAPOE4",
                                                                         "MALEAPOE2", "MALEAPOE3", "MALEAPOE4"),]) + ggtitle("E. Plasma GFAP") 
p6 <- get_hist_of_counts(results_NFL_grouped_bs[results_NFL_grouped_bs$Group %in% c("FEMALEAPOE2", "FEMALEAPOE3", "FEMALEAPOE4",
                                                                             "MALEAPOE2", "MALEAPOE3", "MALEAPOE4"),]) + ggtitle("F. Plasma NFL") 
p7 <- get_hist_of_counts(results_Hippo_grouped_bs[results_Hippo_grouped_bs$Group %in% c("FEMALEAPOE2", "FEMALEAPOE3", "FEMALEAPOE4",
                                                                               "MALEAPOE2", "MALEAPOE3", "MALEAPOE4"),])+ ggtitle("G. Hippocampal Volume")
p8 <- get_hist_of_counts(results_mcrt_grouped_bs[results_mcrt_grouped_bs$Group %in% c("FEMALEAPOE2", "FEMALEAPOE3", "FEMALEAPOE4",
                                                                             "MALEAPOE2", "MALEAPOE3", "MALEAPOE4"),])+ ggtitle("H. mCRT")
lemon::grid_arrange_shared_legend(p1 ,p2, p3, p4, p5, p6, p7, p8, nrow = 4, ncol = 2)
graph2ppt(file = "./figures/Prevalence_by_APOE_SEX_allBiomarkers_forSupplement.pptx", width = 5, height = 5.78)

p1 <- get_hist_of_counts(results_CSF_grouped_bs[results_CSF_grouped_bs$Group %in% c("FEMALEAPOE2", "FEMALEAPOE4",
                                                                                    "MALEAPOE2", "MALEAPOE4"),]) + ggtitle("A. CSF AB42/AB40") 
p2 <- get_hist_of_counts(results_CL_grouped_bs[results_CL_grouped_bs$Group %in% c("FEMALEAPOE2", "FEMALEAPOE4",
                                                                                  "MALEAPOE2",  "MALEAPOE4"),]) + ggtitle("B. Amyloid PET") 
p3 <- get_hist_of_counts(results_pTau_DABNI_bs[results_pTau_DABNI_bs$Group %in% c("FEMALEAPOE2", "FEMALEAPOE4",
                                                                                  "MALEAPOE2", "MALEAPOE4"),]) + ggtitle("C. Plasma pTau217 - AlzPATH") 
p4 <- get_hist_of_counts(results_pTau_ABCDS_bs[results_pTau_ABCDS_bs$Group %in% c("FEMALEAPOE2",  "FEMALEAPOE4",
                                                                                  "MALEAPOE2",  "MALEAPOE4"),]) + ggtitle("D. Plasma pTau217 - Lilly") 
p5 <- get_hist_of_counts(results_GFAP_grouped_bs[results_GFAP_grouped_bs$Group %in% c("FEMALEAPOE2",  "FEMALEAPOE4",
                                                                                      "MALEAPOE2", "MALEAPOE4"),]) + ggtitle("E. Plasma GFAP") 
p6 <- get_hist_of_counts(results_NFL_grouped_bs[results_NFL_grouped_bs$Group %in% c("FEMALEAPOE2", "FEMALEAPOE4",
                                                                                    "MALEAPOE2", "MALEAPOE4"),]) + ggtitle("F. Plasma NFL") 
p7 <- get_hist_of_counts(results_Hippo_grouped_bs[results_Hippo_grouped_bs$Group %in% c("FEMALEAPOE2","FEMALEAPOE4",
                                                                                        "MALEAPOE2", "MALEAPOE4"),])+ ggtitle("G. Hippocampal Volume")
p8 <- get_hist_of_counts(results_mcrt_grouped_bs[results_mcrt_grouped_bs$Group %in% c("FEMALEAPOE2",  "FEMALEAPOE4",
                                                                                      "MALEAPOE2", "MALEAPOE4"),])+ ggtitle("H. mCRT")
lemon::grid_arrange_shared_legend(p1 ,p2, p3, p4, p5, p6, p7, p8, nrow = 4, ncol = 2)
graph2ppt(file = "./figures/Prevalence_by_APOE_SEX_allBiomarkers_noAPOE3_forSupplement.pptx", width = 5, height = 5.78)


get_bs_plot_sex(results_mcrt_grouped_bs[results_mcrt_grouped_bs$Group %in% c("MALEAPOE2"),])+ ggtitle("C. MCRT Total Score < 28")

results_pTau_bs$APOE_grouped <- factor(results_pTau_bs$APOE_grouped, levels = c("APOE4", "APOE3", "APOE2"))
results_pTau_grouped_bs$Group <- factor(results_pTau_grouped_bs$Group, levels = c("FEMALEAPOE4", "FEMALEAPOE3",
                                                                                  "FEMALEAPOE2","MALEAPOE4",
                                                                                   "MALEAPOE3",
                                                                                   "MALEAPOE2"))

p1 <- get_bs_plot(results_pTau_bs[!is.na(results_pTau_bs$APOE_grouped),])+ ggtitle("A. Plasma pTau\nStratified by APOE")
p2 <- get_bs_plot_sex(results_pTau_grouped_bs[results_pTau_grouped_bs$Group %in% c("FEMALEAPOE2", "FEMALEAPOE3", "FEMALEAPOE4",
                                                                                     "MALEAPOE2", "MALEAPOE3", "MALEAPOE4"),])+ ggtitle("B. Plasma pTau\nStratified by APOE and Sex")
final_plot_ptau <- p1 | p2 
final_plot_ptau

graph2ppt(file = "./figures/Supp_Prevalence_pTau217.pptx", width = 10, height = 5.78)

