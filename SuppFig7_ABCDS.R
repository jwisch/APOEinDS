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


ABCDS_Hipp$cohort <- "ABCDS"

################################################################################
################################################################################
ABCDS_CL$de_gender <- as.factor(ABCDS_CL$de_gender)
levels(ABCDS_CL$de_gender) <- c("M", "F")
ABCDS_CL$cohort <- "ABCDS"
levels(ABCDS_CL$Intellectual_Disability) <- c("Mild", "Moderate", "Severe")
CL <- ABCDS_CL[, c("Subject","WUSTLcentiloid",
                         "amy_AgefromBaseline", "APOE_grouped", "cohort", "de_gender", "Intellectual_Disability")]

CL$Subject <- as.factor(CL$Subject)
CL <- CL[!is.na(CL$amy_AgefromBaseline),]
CL$cohort <- factor(CL$cohort, levels = c("DABNI", "ABCDS"))
CL$APOE_grouped <- relevel(factor(CL$APOE_grouped), ref = "APOE3")


CL$WUSTLcentiloid_shift <- CL$WUSTLcentiloid - min(CL$WUSTLcentiloid, na.rm = TRUE) + 0.1


res_CL_cohort_sex <- run_apoe_model_cohort(
  df = CL,
  apoe_col = "APOE_grouped",
  age_col = "amy_AgefromBaseline",
  y_col = "WUSTLcentiloid_shift",
  intellectual_disability_col = "Intellectual_Disability",
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

###############################################################################
##pTau217
###############################################################################

pTau217 <- merge(ABCDS[!duplicated(ABCDS$SUBJECT_LABEL), c("SUBJECT_LABEL", "EVENT_SEQUENCE",
                                                           "pTau217", "APOE_grouped", "de_gender", "Intellectual_Disability")],
                 latency[, c("subject_label", "event_sequence", "clinical_AgefromBaseline")],
                 by.x = c("SUBJECT_LABEL", "EVENT_SEQUENCE"),
                 by.y = c("subject_label", "event_sequence"), all = FALSE)


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

#GFAP
levels(GFAP$Intellectual_Disability) <- c("Mild", "Moderate", "Severe")
GFAP <- merge(GFAP, ABCDS_apoe[, c("subject_label", "APOE_grouped")],
              by.x = "SUBJECT_LABEL", by.y = "subject_label", all = FALSE)
GFAP <- merge(GFAP, latency[, c("subject_label", "event_sequence", "clinical_AgefromBaseline")],
              by.x = c("SUBJECT_LABEL", "EVENT_SEQUENCE"),
              by.y = c("subject_label", "event_sequence"), all = FALSE)
GFAP <- merge(GFAP, demogs[!duplicated(demogs$subject_label), c("subject_label", "de_gender")],
              by.x = c("SUBJECT_LABEL"), by.y = c("subject_label"), all = FALSE)

colnames(GFAP)[c(14:15)] <- c("APOE_grouped","clinical_AgefromBaseline")

GFAP$TESTVALUE <- as.numeric(GFAP$TESTVALUE)

GFAP$cohort <- factor(GFAP$cohort, levels = c("DABNI", "ABCDS"))
GFAP$de_gender <- as.factor(GFAP$de_gender)


res_GFAP_cohort_sex <- run_apoe_model_cohort(
  df = GFAP,
  apoe_col = "APOE_grouped",
  age_col = "clinical_AgefromBaseline",
  intellectual_disability_col = "Intellectual_Disability",
  y_col = "TESTVALUE",
  reference_sex = "M",                                  
  family_choice = Gamma(link = "log"),
  
  sex_col = "de_gender",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = TRUE,
  
  k_val = 5
)



#NFL
NFL <- merge(NFL[, c("SUBJECT_LABEL", "EVENT_SEQUENCE", "TESTVALUE")], ABCDS_apoe[, c("subject_label", "APOE_grouped")],
             by.x = "SUBJECT_LABEL", by.y = "subject_label", all = FALSE)
NFL <- merge(NFL, latency[, c("subject_label", "event_sequence", "clinical_AgefromBaseline")],
             by.x = c("SUBJECT_LABEL", "EVENT_SEQUENCE"),
             by.y = c("subject_label", "event_sequence"), all = FALSE)
NFL <- merge(NFL, demogs[!duplicated(demogs$subject_label),c("subject_label", "de_gender")],
             by.x = "SUBJECT_LABEL", by.y = "subject_label", all = FALSE)
NFL <- merge(NFL, prefunc, by.x = "SUBJECT_LABEL", by.y = "subject_label", all = FALSE)
levels(NFL$Intellectual_Disability) <- c("Mild", "Moderate", "Severe") 
NFL$TESTVALUE <- as.numeric(NFL$TESTVALUE)
NFL$de_gender <- as.factor(NFL$de_gender)

NFL$TESTVALUE <- log(NFL$TESTVALUE)



res_NFL_cohort_sex <- run_apoe_model_cohort(
  df = NFL,
  apoe_col = "APOE_grouped",
  age_col = "clinical_AgefromBaseline",
  sex_col = "de_gender",
  y_col = "TESTVALUE",
  intellectual_disability_col = "Intellectual_Disability",
  reference_sex = "M",                                  
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF",  "#33CC33","#FF4D4D"),
  k_val = 4,
  markSig = TRUE,
  
)


#Hippocampal volume
ABCDS_Hipp <- merge(ABCDS_Hipp, demogs[!duplicated(demogs$subject_label), c("subject_label", "de_gender")],
                    by.x = "Subject.ID", by.y = "subject_label", all = FALSE)
levels(ABCDS_Hipp$Intellectual_Disability) <- c("Mild", "Moderate", "Severe")
hippo <- ABCDS_Hipp[, c("Subject.ID", "mri_AgefromBaseline",
                              "APOE_grouped", "de_gender",
                              "cohort",  "hippo_Z", "Intellectual_Disability")]
hippo$de_gender <- as.factor(hippo$de_gender)
hippo$APOE_grouped <- as.factor(hippo$APOE_grouped)
levels(hippo$APOE_grouped) <- c("APOE2", "APOE2", "APOE4",
                                "APOE3", "APOE4", "APOE4", 
                                "APOE4", "APOE3", "APOE2")


res_hippo_cohort_sex <- run_apoe_model_cohort(
  df = hippo,
  apoe_col = "APOE_grouped",
  age_col = "mri_AgefromBaseline",
  y_col = "hippo_Z",
  sex_col = "de_gender",
  reference_sex = "M",                                  
  
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  k_val = 4,
  markSig = TRUE
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
mcrt <- merge(mcrt, prefunc, by = "subject_label", all = FALSE)
levels(mcrt$Intellectual_Disability) <- c("Mild", "Moderate", "Severe")

mcrt <-  mcrt[, c("subject_label",  "clinical_AgefromBaseline",
                       "trs", "APOE_grouped", "de_gender", "cohort", "Intellectual_Disability")]
mcrt$trs <- as.numeric(mcrt$trs)
mcrt$de_gender <- as.factor(mcrt$de_gender)


# 1. Convert score to work with distribution
mcrt$score_prop <- 36 - mcrt$trs



res_mcrt_cohort_sex <- run_apoe_model_cohort(
  df = mcrt,
  apoe_col = "APOE_grouped",
  age_col = "clinical_AgefromBaseline",
  y_col = "score_prop",
  family = tw(link = "log"),
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
    group = APOE_grouped
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
      group = APOE_grouped
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




lemon::grid_arrange_shared_legend(res_CL_cohort_sex$plot + ggtitle("A."), 
                                  res_ptau217_ABCDS_sex$plot + ggtitle("B."),
                                  res_GFAP_cohort_sex$plot + ggtitle("C."), res_NFL_cohort_sex$plot + ggtitle("D."), 
                                  res_hippo_cohort_sex$plot + ggtitle("E."), p_mcrt_original + ggtitle("F."), nrow = 3, ncol = 2)

graph2ppt(file = "./figures/Fig1_ABCDSonly.pptx", width = 6, height = 11)