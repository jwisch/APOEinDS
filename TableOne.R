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
library(stringr)
source("./corrections_ABCDS.R")
source("./code_to_plot.R")
source("./DABNI_maximizingAPOE.R")
#https://alz-journals.onlinelibrary.wiley.com/doi/10.1002/alz.13859
#Core 1, Amyloid: Amyloid PET
#Core 1, T1: p-Tau217
#Core 2, T2: Tau PET, ABCDS Only
#I: GFAP
#N: Anatomic MRI - Hippocampal Volume
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

# DABNI_dx <- readxl::read_xlsx("./Data_202512/NPS_20251127.xlsx")
DABNI_dx$crt_rli_rfi_total <- as.numeric(DABNI_dx$crt_rli_rfi_total)
DABNI_dx <- DABNI_dx[!is.na(DABNI_dx$crt_rli_rfi_total),]

tmp <- setDT(DABNI_dx)[, .N, by = "NHC"]
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


pTau217_thresh_DABNI <- mean(DABNI[DABNI$Age_at_CSF <= 35 & !duplicated(DABNI$NHC),]$PLASMA_PTAU217_S0217, na.rm = TRUE) + 
  1.96 *sd(DABNI[DABNI$Age_at_CSF <= 35 & !duplicated(DABNI$NHC),]$PLASMA_PTAU217_S0217, na.rm = TRUE) 

pTau217_thresh_ABCDS <- mean(ABCDS[ABCDS$age_at_visit <= 35 & !duplicated(ABCDS$SUBJECT_LABEL),]$pTau217, na.rm = TRUE) + 
  1.96 *sd(ABCDS[ABCDS$age_at_visit <= 35 & !duplicated(ABCDS$SUBJECT_LABEL),]$pTau217, na.rm = TRUE) 


CL_thresh_DABNI <- mean(DABNI_CL[DABNI_CL$Age_at_PET <= 35 & !duplicated(DABNI_CL$NHC),]$CL_SPM12, na.rm = TRUE) + 
  1.96 *sd(DABNI_CL[DABNI_CL$Age_at_PET <= 35 & !duplicated(DABNI_CL$NHC),]$CL_SPM12, na.rm = TRUE) 

ABCDS_CL <- read.csv(".././ABCDS_DF3_AmyloidOnly/all_amy_120924.csv")
ABCDS_CL <- ABCDS_CL[ABCDS_CL$PET.TC.QC.Status == "Passed", c("Subject", "VISIT", "WUSTLcentiloid")]
ABCDS_CL$event_sequence <- as.numeric(substr(ABCDS_CL$VISIT, start = nchar(ABCDS_CL$VISIT), stop = nchar(ABCDS_CL$VISIT)))
ABCDS_CL <- merge(ABCDS_CL, demogs, by.x = "Subject", by.y = "subject_label", all.x = TRUE, all.y = FALSE)
ABCDS_CL <- merge(ABCDS_CL, ABCDS_apoe, by.x = "Subject", by.y = "subject_label", all = FALSE)

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

ABCDS <- ABCDS[!duplicated(ABCDS$SUBJECT_LABEL),]

ABCDS_apoe <- merge(ABCDS_apoe[, c("subject_label", "allele_combo", "APOE_grouped"),],
                    ABCDS[, c("SUBJECT_LABEL", "EVENT_SEQUENCE", "pTau217")],
                    by.x = "subject_label", by.y = "SUBJECT_LABEL", all.x = TRUE, all.y = FALSE)
colnames(ABCDS_apoe)[4] <- c("Event_Sequence_pTau217")
rm(ABCDS)
ABCDS_CL <- ABCDS_CL[order(ABCDS_CL$Subject, ABCDS_CL$event_sequence.x),]
ABCDS_apoe <- merge(ABCDS_apoe,
                    ABCDS_CL[!duplicated(ABCDS_CL$Subject), c("Subject", "event_sequence.x", "WUSTLcentiloid")],
                    by.x = "subject_label", by.y = "Subject", all.x = TRUE, all.y = FALSE)
colnames(ABCDS_apoe)[6] <- c("Event_Sequence_CL")
rm(ABCDS_CL)

rm(ABCDS_corr)

ABCDS_apoe <- merge(ABCDS_apoe, ABCDS_Hipp[!duplicated(ABCDS_Hipp$Subject.ID), c("Subject.ID", "event_sequence", "hippo")],
                    by.x = "subject_label", by.y = "Subject.ID", all.x = TRUE, all.y = FALSE)
mcrt <- mcrt[order(mcrt$subject_label, mcrt$event_sequence),]
ABCDS_apoe <- merge(ABCDS_apoe, mcrt[!duplicated(mcrt$subject_label), c("subject_label", "event_sequence", "trs")],
                    by = "subject_label", all.x = TRUE, all.y = FALSE)
dx <- dx[order(dx$subject_label, dx$event_sequence.y),]
ABCDS_apoe <- merge(ABCDS_apoe, dx[!duplicated(dx$subject_label), c("subject_label", "event_sequence.y", "consensus_dx")],
                    by = "subject_label", all.x = TRUE, all.y = FALSE)

rm(mcrt, dx, ABCDS_Hipp)
colnames(ABCDS_apoe)[c(8, 10, 12)] <- c("Event_Sequence_Hippo", "Event_Sequence_trs", "Event_Sequence_dx")
ABCDS <- merge(ABCDS_apoe, demogs[!duplicated(demogs$subject_label), c("subject_label", "de_gender", "ds_vs_control_flag")],
               by = "subject_label", all.x = TRUE, all.y = FALSE)
ABCDS <- merge(ABCDS, latency[, c("subject_label", "event_sequence", "clinical_AgefromBaseline")],
               by.x = c("subject_label", "Event_Sequence_dx"), by.y = c("subject_label", "event_sequence"), all.x = TRUE, all.y = FALSE)
ABCDS <- merge(ABCDS, GFAP[, c("SUBJECT_LABEL", "EVENT_SEQUENCE", "TESTVALUE")], 
               by.x = "subject_label", by.y = "SUBJECT_LABEL", all.x = TRUE, all.y = FALSE)
colnames(ABCDS)[17:18] <- c("Event_Sequence_GFAP", "GFAP")
ABCDS <- merge(ABCDS, NFL[, c("SUBJECT_LABEL", "EVENT_SEQUENCE", "TESTVALUE")],
               by.x = "subject_label", by.y = "SUBJECT_LABEL", all.x = TRUE, all.y = FALSE)
colnames(ABCDS)[19:20] <- c("Event_Sequence_NFL", "NFL")
rm(ABCDS_apoe, corr_sel, IDlist, demogs, GFAP, NFL)
ABCDS$cohort <- "ABCDS"


DABNI <- merge(DABNI_dx, DABNI, by = "NHC", all = TRUE)
DABNI <- merge(DABNI_CL, DABNI, by = "NHC", all = TRUE)

DABNI <- merge(DABNI[, c("NHC", "CL_SPM12", "Age_at_PET", "Age_NPS", "Diag", "crt_rli_rfi_total",
                         "Age_at_CSF", "PLASMA_PTAU217_S0217", "PLASMA_GFAP_S0122",
                         "PLASMA_NFLIGHT_S0100", "CSF_AB142_A0010")],
               df, by = "NHC", all = FALSE)
DABNI <- DABNI[!is.na(DABNI$NHC),]
DABNI <- DABNI[!duplicated(DABNI$NHC),]
# rm(DABNI_dx, DABNI_CL)
# DABNI_Hipp <- DABNI_Hipp[order(DABNI_Hipp$NHC, DABNI_Hipp$Age_at_MRI),]
# 
DABNI <- merge(DABNI, DABNI_Hipp[, c("NHC", "Age_at_MRI",
                                     "FS_Left-Hippocampus",
                                     "FS_Right-Hippocampus",
                                     "FS_EstimatedTotalIntraCranialVol")], by = "NHC", all = TRUE)
DABNI <- DABNI[!duplicated(DABNI$NHC),]
# row_mode <- function(x) {
#   x <- x[!is.na(x)]
#   if (length(x) == 0) return(NA)
#   ux <- unique(x)
#   ux[which.max(tabulate(match(x, ux)))]
# }
# # column names
# cn <- colnames(df)
# 
# # extract base names for .x / .y variables
# base_vars <- cn[grepl("\\.(x|y)$", cn)] %>%
#   str_replace("\\.(x|y)$", "") %>%
#   unique()
# 
# for (v in base_vars) {
#   
#   cols <- grep(paste0("^", v, "\\.(x|y)$"), names(df), value = TRUE)
#   cols <- unique(cols)   # <-- THIS FIXES THE ERROR
#   
#   df[cols] <- lapply(df[cols], as.character)
#   
#   df[[v]] <- apply(df[cols], 1, row_mode)
# }
# 
# 
# df <- df[, c("NHC", "APOE", "Sex", "Diag", "CL_SPM12", "crt_rli_rfi_total",
#              "Age_at_CSF", "PLASMA_PTAU217_S0217", "PLASMA_GFAP_S0122",               
#              "PLASMA_NFLIGHT_S0100", "CSF_AB142_A0010", "Age_at_MRI", 
#              "FS_Left-Hippocampus", "FS_Right-Hippocampus",
#              "FS_EstimatedTotalIntraCranialVol")]

# DABNI <- df

names(DABNI) <- c("subject_label","WUSTLcentiloid", "Age_at_PET",
                  "clinical_AgefromBaseline","consensus_dx",  "trs",
                  "Age_at_CSF", "pTau217", "GFAP",
                  "NFL", "CSFAB42AB40","APOE", 
                   "de_gender",  "Age_at_MRI",
                  "Hippo_left", "Hippo_right", "ICV")
                  
DABNI$cohort <- "DABNI"
DABNI$hippo <- DABNI$Hippo_left + DABNI$Hippo_right
mod <- lm(hippo ~ ICV, data = DABNI, model = TRUE)
DABNI$hippo_resid <- NA
# Step 4: Insert residuals using the "na.action" attribute
used_rows <- as.numeric(rownames(model.frame(mod)))  # TRUE numeric row indices

DABNI$hippo_resid[used_rows] <- resid(mod)
DABNI$APOE <- as.factor(DABNI$APOE)
DABNI$APOE_grouped <- recode(DABNI$APOE,
                             "22" = "APOE2",
                             "23" = "APOE2",
                             "24" = "APOE4",
                             "33" = "APOE3",
                             "34" = "APOE4",
                             "43" = "APOE4",
                             "44" = "APOE4")
ABCDS$CSFAB42AB40 <- NA

DABNI$APOE <- recode(DABNI$APOE,
                             "22" = "22",
                             "23" = "23",
                             "24" = "24",
                             "33" = "33",
                             "34" = "34",
                             "43" = "34",
                             "44" = "44")

ABCDS$APOE <- recode(ABCDS$allele_combo,
                     "E2/E2" = "22", 
                     "E2/E3" = "23",
                     "E2/E4" = "24",
                     "E3/3E" = "33",
                     "E3/E2" = "23",
                     "E3/E3" = "33",
                     "E3/E4" = "34",
                     "E4/E2" = "24",
                     "E4/E3" = "34",
                     "E4/E4" = "44")

ABCDS <- ABCDS[!duplicated(ABCDS$subject_label),]

df <- rbind(DABNI[, c("subject_label", "de_gender", "APOE", "APOE_grouped", "clinical_AgefromBaseline", 
                      "WUSTLcentiloid", "CSFAB42AB40", "pTau217", "GFAP", "NFL", "trs", "consensus_dx", "hippo", "cohort") ],
            ABCDS[ABCDS$ds_vs_control_flag == "DS", c("subject_label", "de_gender", "APOE", "APOE_grouped", "clinical_AgefromBaseline", 
            "WUSTLcentiloid", "CSFAB42AB40", "pTau217", "GFAP", "NFL", "trs", "consensus_dx", "hippo", "cohort")])


df <- df[!is.na(df$APOE),]
# df <- df[!duplicated(df$subject_label),]

df$hasMRI <- ifelse(!is.na(df$hippo), 1, 0)
df$hasCL <- ifelse(!is.na(df$WUSTLcentiloid), 1, 0)
df$hasCSF <- ifelse(!is.na(df$CSFAB42AB40), 1, 0)
df$haspTau217 <- ifelse(!is.na(df$pTau217), 1, 0)
df$hasGFAP <- ifelse(!is.na(df$GFAP), 1, 0)
df$hasNFL <- ifelse(!is.na(df$NFL), 1, 0)
df$hasmCRT <- ifelse(!is.na(df$trs), 1, 0)
df$hasDX <- ifelse(!is.na(df$consensus_dx), 1, 0)
df$dx_grouped <- recode(df$consensus_dx,
                             "0" = "asymptomatic",
                             "1" = "symptomatic",
                             "2" = "symptomatic",
                             "3" = "no consensus",
                             "aDS" = "asymptomatic",
                             "dDS" = "symptomatic",
                             "pDS" = "symptomatic",
                        "uDS" = "no consensus")
df$de_gender <- recode(df$de_gender,
                        "1" = "M",
                       "2" = "F",
                       "F" = "F", 
                       "M" = "M")
df$hasMRI <- ifelse(!is.na(df$hippo), 1, 0)

vars <- c("de_gender", "APOE", "clinical_AgefromBaseline", "dx_grouped",
          "hasCL", "hasCSF", "haspTau217", "hasGFAP", "hasNFL", "hasMRI", "hasmCRT", "hasDX")
catVars <- c("de_gender", "APOE", "dx_grouped",
             "hasCL", "hasCSF", "haspTau217", "hasGFAP", "hasNFL", "hasMRI", "hasmCRT", "hasDX")
CreateTableOne(vars = vars, factorVars = catVars, strata = "cohort", data = df)


df$gender_APOE <- as.factor(paste0(df$de_gender, df$APOE_grouped))
vars <- c( "clinical_AgefromBaseline", "dx_grouped",
          "hasCL", "hasCSF", "haspTau217", "hasGFAP", "hasNFL", "hasMRI", "hasmCRT", "hasDX", "cohort")
catVars <- c("dx_grouped",
             "hasCL", "hasCSF", "haspTau217", "hasGFAP", "hasNFL", "hasMRI", "hasmCRT", "hasDX", "cohort")
CreateTableOne(vars = vars, factorVars = catVars, strata = "gender_APOE", data = df)

CreateTableOne(vars = vars, factorVars = catVars, strata = "APOE_grouped", data = df)



library(ComplexUpset)

cols <- c("hasCL", "hasCSF", "haspTau217", "hasGFAP",
          "hasNFL", "hasMRI", "hasmCRT", "hasDX")

df_upset <- df %>%
  mutate(across(all_of(cols), ~ as.integer(.x == 1)),
         cohort = factor(cohort))

ComplexUpset::upset(
  df_upset,
  cols,
  name = "Modalities",
  min_size = 10,
  base_annotations = list(
    'Intersection size' = intersection_size(
      counts = TRUE,
      mapping = aes(fill = cohort)
    )
  )  
) + scale_fill_manual(values = c("#F8766D",   # blue
                                   "#619CFF")) +
  scale_colour_manual(values = c("#F8766D",   # blue
                               "#619CFF")) + theme_bw()

graph2ppt(file = "./figures/UPSETR.pptx", width = 10, height = 5.78)

UpSetR::upset(
  df_upset[, cols],
  ) + theme_bw()


df$age_cut <- cut(df$clinical_AgefromBaseline, c(0, 35, 40, 45, 50, 55, 75))

values = c(
  "APOE2.asymptomatic" = "#CCE5FF",
  "APOE2.symptomatic"   = "#3399FF",
  "APOE3.asymptomatic"     = "#D6EFD6",
  "APOE3.symptomatic"      = "#33CC33",
  "APOE4.asymptomatic" = "#FFD6D6",
  "APOE4.symptomatic"   = "#FF4D4D"
)
labels = c(
  "APOE2 — Asymptomatic",
  "APOE2 — Symptomatic",
  "APOE3 — Asymptomatic",
  "APOE3 — Symptomatic",
  "APOE4 — Asymptomatic",
  "APOE4 — Symptomatic"
)


group_levels <- c("APOE2", "APOE3", "APOE4")
obj <- plot_group_proportions(df,"APOE_grouped", group_levels, VALS = values, LABS = labels)
obj
graph2ppt(file = "./figures/Fig3.pptx", width = 8, height = 8)

plot_group_proportions(df[df$de_gender == "F",],"APOE_grouped", group_levels, VALS = values, LABS = labels)
graph2ppt(file = "./figures/Fig3_Analog_forSupp_Females.pptx", width = 8, height = 8)

plot_group_proportions(df[df$de_gender == "M",],"APOE_grouped", group_levels, VALS = values, LABS = labels)
graph2ppt(file = "./figures/Fig3_Analog_forSupp_Males.pptx", width = 8, height = 8)

p1 <- plot_group_proportions(df[df$APOE_grouped == "APOE2",], "de_gender", group_levels = c("M", "F"),
                       VALS = c("M.asymptomatic" = "#CCE5FF", "M.symptomatic"= "#3399FF",
                                "F.asymptomatic" = "#CCE5FF", "F.symptomatic"= "#3399FF"), 
                       LABS = c("Male - Asymptomatic", "Male - Symptomatic",
                                "Female - Asymptomatic", "Female - Symptomatic"))

p2 <- plot_group_proportions(df[df$APOE_grouped == "APOE3",], "de_gender", group_levels = c("M", "F"),
                       VALS = c("M.asymptomatic" = "#D6EFD6", "M.symptomatic"= "#33CC33",
                                "F.asymptomatic" = "#D6EFD6", "F.symptomatic"= "#33CC33"), 
                       LABS = c("Male - Asymptomatic", "Male - Symptomatic",
                                "Female - Asymptomatic", "Female - Symptomatic"))

p3 <- plot_group_proportions(df[df$APOE_grouped == "APOE4",], "de_gender", group_levels = c("M", "F"),
                       VALS = c("M.asymptomatic" = "#FFD6D6", "M.symptomatic"= "#FF4D4D",
                                "F.asymptomatic" = "#FFD6D6", "F.symptomatic"= "#FF4D4D"), 
                       LABS = c("Male - Asymptomatic", "Male - Symptomatic",
                                "Female - Asymptomatic", "Female - Symptomatic"))
library(gridExtra)

grid.arrange(p1[[1]], p2[[1]], p3[[1]], nrow = 1)
graph2ppt(file = "./figures/Fig3_Analog_forSupp_Everyone.pptx", width = 11, height = 8)

plot_group_proportions(df[df$cohort == "ABCDS",],"APOE_grouped", group_levels, VALS = values, LABS = labels)
graph2ppt(file = "./figures/Fig3_Analog_forSupp_ABCDSonly.pptx", width = 8, height = 8)
plot_group_proportions(df[df$cohort == "DABNI",],"APOE_grouped", group_levels, VALS = values, LABS = labels)
graph2ppt(file = "./figures/Fig3_Analog_forSupp_DABNIonly.pptx", width = 8, height = 8)


library(mgcv)
df$APOE_grouped <- factor(df$APOE_grouped, levels = c("APOE2", "APOE3", "APOE4"))
res_CSF <- run_apoe_model(
  df = df,
  apoe_col = "APOE_grouped",
  y_col = "CSFAB42AB40",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D")
)
res_CL <- run_apoe_model(
  df = df,
  apoe_col = "APOE_grouped",
  y_col = "WUSTLcentiloid",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D")
)

res_CL$p
res_CL$summary_table
# 
# 
# res_GFAP <- run_apoe_model(
#   df = df,
#   apoe_col = "APOE_grouped",
#   y_col = "GFAP",
#   fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
#   colour_values = c("#3399FF", "#33CC33", "#FF4D4D")
# )
# 
# res_NFL <- run_apoe_model(
#   df = df,
#   apoe_col = "APOE_grouped",
#   y_col = "NFL",
#   fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
#   colour_values = c("#3399FF", "#33CC33", "#FF4D4D")
# )
# 
# res_hippo <- run_apoe_model(
#   df = df,
#   apoe_col = "APOE_grouped",
#   y_col = "hippo",
#   fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
#   colour_values = c("#3399FF", "#33CC33", "#FF4D4D")
# )
# 
# res_mcrt<- run_apoe_model(
#   df = df,
#   apoe_col = "APOE_grouped",
#   y_col = "trs",
#   fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
#   colour_values = c("#3399FF", "#33CC33", "#FF4D4D")
# )
