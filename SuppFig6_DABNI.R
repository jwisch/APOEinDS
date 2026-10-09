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


#GFAP
GFAP <- DABNI[!is.na(DABNI$PLASMA_GFAP_S0122),c("NHC", "Age_at_CSF", "Sex", "PLASMA_GFAP_S0122", "APOE_grouped", "Intellectual_Disability")]
colnames(GFAP) <- c("SUBJECT_LABEL", "clinical_AgefromBaseline", "de_gender", "TESTVALUE", "APOE_grouped", "Intellectual_Disability")
GFAP$cohort <- "DABNI"
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
NFL <- DABNI[!is.na(DABNI$PLASMA_NFLIGHT_S0100),c("NHC", "Age_at_CSF", "PLASMA_NFLIGHT_S0100", "APOE_grouped", "Sex", "Intellectual_Disability")]
colnames(NFL) <- c("SUBJECT_LABEL", "clinical_AgefromBaseline", "TESTVALUE", "APOE_grouped", "de_gender", "Intellectual_Disability")
NFL$cohort <- "DABNI"
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
hippo <- DABNI_Hipp[, c("Subject.ID", "mri_AgefromBaseline",
                              "APOE_grouped", "de_gender",
                              "hippo_Z", "Intellectual_Disability")]
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


DABNI_dx <- readxl::read_xlsx("./Data_202512/NPS_20251127.xlsx")
DABNI_dx <- DABNI_dx[, c("NHC", "Sex", "Age_NPS", "Diag", "crt_rli_rfi_total", "Intellectual_Disability")]
DABNI_dx <- merge(DABNI_dx, DABNI[, c("NHC", "APOE_grouped")], by = "NHC", all = FALSE)
DABNI_dx <- DABNI_dx[order(DABNI_dx$NHC, DABNI_dx$Age_NPS),]
DABNI_dx <- DABNI_dx[!duplicated(DABNI_dx$NHC),]
DABNI_dx$cohort <- "DABNI"
colnames(DABNI_dx) <- c("subject_label", "de_gender", "clinical_AgefromBaseline",
                        "diag", "trs", "Intellectual_Disability", "APOE_grouped", "cohort")


mcrt <- DABNI_dx[, c("subject_label",  "clinical_AgefromBaseline",
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
      group = APOE_grouped    ),
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
lay <- rbind(c(1, 2), c(3, 3), c(4, 5), c(6, 7))
lemon::grid_arrange_shared_legend(res_CSF_sex$plot + aes(shape = APOE_grouped) + 
                                    scale_shape_manual(values = c(1, 1, 1), guide = "none") + ggtitle("A."), 
                                  res_CL_cohort_sex$plot + ggtitle("B."), 
                                  res_ptau217_DABNI_sex$plot + aes(shape = APOE_grouped) + 
                                    scale_shape_manual(values = c(1, 1, 1), guide = "none") +  ggtitle("C."),
                                  res_GFAP_cohort_sex$plot + ggtitle("D."), res_NFL_cohort_sex$plot + ggtitle("E."), 
                                  res_hippo_cohort_sex$plot + ggtitle("F."), p_mcrt_original + ggtitle("G."), layout_matrix = lay)

graph2ppt(file = "./figures/Fig1_DABNIonly.pptx", width = 6, height = 11)
