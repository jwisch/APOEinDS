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
colnames(DABNI_CL) <- c("ID", "CL", "Age", "APOE", "SEX", "Apos", "APOE_grouped")

DABNI_Hipp$APOE_grouped <- ifelse(DABNI_Hipp$APOE %in% c("22", "23"), "APOE2",
                                  ifelse(DABNI_Hipp$APOE %in% c("33"), "APOE3",
                                         ifelse(DABNI_Hipp$APOE %in% c("24", "34", "44"), "APOE4", NA)))
DABNI_dx$APOE_grouped <- ifelse(DABNI_dx$APOE %in% c("22", "23"), "APOE2",
                                ifelse(DABNI_dx$APOE %in% c("33"), "APOE3",
                                       ifelse(DABNI_dx$APOE %in% c("24", "34", "44"), "APOE4", NA)))



ABCDS$pTau217_positive <- ifelse(ABCDS$pTau217 > pTau217_thresh_ABCDS, 1, 0)
ABCDS_CL$CL_positive <- ifelse(ABCDS_CL$WUSTLcentiloid > CL_thresh_ABCDS, 1, 0)
ABCDS_Hipp$hippo_positive <- ifelse(ABCDS_Hipp$hippo_Z < hippo_thresh_ABCDS, 1, 0)
colnames(ABCDS_CL) <- c("ID", "event_seq", "VISIT", "CL", "Age", "SEX", "APOE", "updated", "APOE_grouped", "Apos")

DABNI_dx$mcrt_dementia <- ifelse(DABNI_dx$crt_rli_rfi_total < 23, 1, 0)
DABNI_dx$mcrt_MCI <- ifelse(DABNI_dx$crt_rli_rfi_total < 28 & DABNI_dx$crt_rli_rfi_total > 22, 1, 0)
DABNI_dx$mcrt_impaired <- ifelse(DABNI_dx$crt_rli_rfi_total < 28, 1, 0)

ABCDS$mcrt_dementia <- ifelse(ABCDS$trs < 23, 1, 0)
ABCDS$mcrt_MCI <- ifelse(ABCDS$trs < 28 & ABCDS$trs > 22, 1, 0)
ABCDS$mcrt_impaired <- ifelse(ABCDS$trs < 28, 1, 0)





get_hist_of_counts <- function(DF, YLIM = NULL, nvals = 3){
  
  p <- ggplot(DF, aes(x = window_center, y = n, fill = Group)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(y = "Count") +
    xlim(c(25, 55)) 
  

  # optionally add y-axis limits
  if (!is.null(YLIM)) {
    p <- p + ylim(YLIM)
  }
  
  return(p)
}


plot_grouped_histogram <- function(
    df,
    value_col,
    group_col = "Group",
    fill_colors = NULL,
    binwidth = NULL,
    bins = 30,
    title = NULL,
    x_label = NULL,
    nvals = 3,
    XLAB = NULL
) {

  # Remove NAs for clean plotting
  plot_df <- df %>%
    filter(!is.na(.data[[value_col]]), !is.na(.data[[group_col]])) %>%
    mutate(!!sym(group_col) := factor(.data[[group_col]]))
  
  # Build base plot
  p <- ggplot(plot_df, aes(x = .data[[value_col]], fill = .data[[group_col]])) +
    geom_histogram(
      aes(y = ..count..),
      binwidth = binwidth,
      bins = if (is.null(binwidth)) bins else NULL,
      color = "black",
      alpha = 0.7,
      position = "stack"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = "none", axis.title.y = element_blank()
    ) + labs(x = XLAB)
  
  # Apply custom colors if supplied
  if(nvals == 3){
    p <- p +
      scale_fill_manual(values = c(
        "#619CFF",
        "#00BA38",
        "#F8766D"
      ))
  }
  if(nvals == 2){
    p <- p +
      scale_fill_manual(values = c(
        "#619CFF",
        "#F8766D"
      ))
  }
  
  return(p)
}


DABNI_CSF <- DABNI[!is.na(DABNI$CSF_AB142_A0010),]
DABNI_CSF <- DABNI_CSF %>% rename(ID = NHC, Age = Age_at_CSF, APOE_grouped = APOE_grouped,
                                  SEX = Sex)
DABNI_CSF$SEX <- as.factor(DABNI_CSF$SEX)
levels(DABNI_CSF$SEX) <- c("FEMALE", "MALE")
DABNI_CSF$Group <- paste0(DABNI_CSF$APOE_grouped)


hist_CSF <- plot_grouped_histogram(DABNI_CSF[DABNI_CSF$Group %in% c("APOE2", "APOE3", "APOE4"),], value_col = "CSF_AB142_A0010", XLAB = "CSF AB42/AB40")
hist_CSF_no3 <- plot_grouped_histogram(DABNI_CSF[DABNI_CSF$Group %in% c("APOE2",  "APOE4"),], value_col = "CSF_AB142_A0010", XLAB = "CSF AB42/AB40", nvals = 2)




CL <- rbind(DABNI_CL[, c("ID", "Age", "APOE_grouped", "SEX", "Apos", "CL")],
            ABCDS_CL[, c("ID", "Age", "APOE_grouped", "SEX", "Apos", "CL")])
CL$SEX <- as.factor(CL$SEX)
levels(CL$SEX) <- c("MALE", "FEMALE", "FEMALE", "MALE")

CL$Group <- paste0(CL$APOE_grouped)
DABNI_CL$Group <- DABNI_CL$APOE_grouped
ABCDS_CL$Group <- ABCDS_CL$APOE_grouped


hist_CL <- plot_grouped_histogram(CL[CL$Group %in% c("APOE2", "APOE3", "APOE4"),], value_col = "CL", XLAB = "Cortical Amyloid Burden (CL)")
hist_CL_DABNI <- plot_grouped_histogram(DABNI_CL[DABNI_CL$Group %in% c("APOE2", "APOE3", "APOE4"),], value_col = "CL", XLAB = "Cortical Amyloid Burden (CL)")
hist_CL_ABCDS <- plot_grouped_histogram(ABCDS_CL[ABCDS_CL$Group %in% c("APOE2", "APOE3", "APOE4"),], value_col = "CL", XLAB = "Cortical Amyloid Burden (CL)")

hist_CL_no3s <- plot_grouped_histogram(CL[CL$Group %in% c("APOE2",  "APOE4"),], value_col = "CL", nvals = 2, XLAB = "Cortical Amyloid Burden (CL)")
hist_CL_DABNI_no3s <- plot_grouped_histogram(DABNI_CL[DABNI_CL$Group %in% c("APOE2",  "APOE4"),], value_col = "CL",  nvals = 2,XLAB = "Cortical Amyloid Burden (CL)")
hist_CL_ABCDS_no3s <- plot_grouped_histogram(ABCDS_CL[ABCDS_CL$Group %in% c("APOE2", "APOE4"),], value_col = "CL",  nvals = 2,XLAB = "Cortical Amyloid Burden (CL)")


DABNI_pTau217 <- DABNI %>% rename(ID = NHC, Age = Age_at_CSF, APOE_grouped = APOE_grouped,
                                  SEX = Sex, Tpos = pTau217_positive)
ABCDS_pTau217 <- ABCDS %>% rename(ID = SUBJECT_LABEL, Age = age_at_visit, APOE_grouped = APOE_grouped,
                                  SEX = de_gender, Tpos = pTau217_positive)
DABNI_pTau217$Group <- DABNI_pTau217$APOE_grouped
ABCDS_pTau217$Group <- ABCDS_pTau217$APOE_grouped

hist_pTau_DABNI <- plot_grouped_histogram(DABNI_pTau217[DABNI_pTau217$Group %in% c("APOE2", "APOE3", "APOE4"),], value_col = "PLASMA_PTAU217_S0217", XLAB = "Plasma pTau217 [AlzPATH]")
hist_pTau_ABCDS <- plot_grouped_histogram(ABCDS_pTau217[!duplicated(ABCDS_pTau217$ID) & ABCDS_pTau217$Group %in% c("APOE2", "APOE3", "APOE4"),], value_col = "pTau217", XLAB = "Plasma pTau217 [Lilly]")

hist_pTau_DABNI_no3s <- plot_grouped_histogram(DABNI_pTau217[DABNI_pTau217$Group %in% c("APOE2", "APOE4"),], value_col = "PLASMA_PTAU217_S0217", XLAB = "Plasma pTau217 [AlzPATH]", nvals = 2)
hist_pTau_ABCDS_no3s <- plot_grouped_histogram(ABCDS_pTau217[!duplicated(ABCDS_pTau217$ID) & ABCDS_pTau217$Group %in% c("APOE2",  "APOE4"),], value_col = "pTau217", XLAB = "Plasma pTau217 [Lilly]", nvals = 2)


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

colnames(DABNI_GFAP)[5] <- "GFAP"
colnames(ABCDS_GFAP)[3] <- "GFAP"
GFAP <- rbind(DABNI_GFAP[, c("ID", "Age", "APOE_grouped", "SEX", "GFAP")], 
              ABCDS_GFAP[, c("ID", "Age", "APOE_grouped", "SEX", "GFAP")])
GFAP$GFAP <- as.numeric(GFAP$GFAP)
ABCDS_GFAP$GFAP <- as.numeric(ABCDS_GFAP$GFAP)

GFAP$Group <- paste0( GFAP$APOE_grouped)
DABNI_GFAP$Group <- paste0( DABNI_GFAP$APOE_grouped)
ABCDS_GFAP$Group <- paste0( ABCDS_GFAP$APOE_grouped)

hist_GFAP <- plot_grouped_histogram(GFAP[GFAP$Group %in% c("APOE2", "APOE3", "APOE4"),], value_col = "GFAP", XLAB = "Plasma GFAP")
hist_GFAP_DABNI <- plot_grouped_histogram(DABNI_GFAP[DABNI_GFAP$Group %in% c("APOE2", "APOE3", "APOE4"),], value_col = "GFAP", XLAB = "Plasma GFAP")
hist_GFAP_ABCDS <- plot_grouped_histogram(ABCDS_GFAP[ABCDS_GFAP$Group %in% c("APOE2", "APOE3", "APOE4"),], value_col = "GFAP", XLAB = "Plasma GFAP")

hist_GFAP_no3s <- plot_grouped_histogram(GFAP[GFAP$Group %in% c("APOE2", "APOE4"),], value_col = "GFAP", XLAB = "Plasma GFAP", nvals = 2)
hist_GFAP_DABNI_no3s <- plot_grouped_histogram(DABNI_GFAP[DABNI_GFAP$Group %in% c("APOE2", "APOE4"),], value_col = "GFAP", XLAB = "Plasma GFAP", nvals = 2)
hist_GFAP_ABCDS_no3s <- plot_grouped_histogram(ABCDS_GFAP[ABCDS_GFAP$Group %in% c("APOE2",  "APOE4"),], value_col = "GFAP", XLAB = "Plasma GFAP", nvals = 2)



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
                                  SEX = de_gender, NFL = TESTVALUE)
DABNI_NFL <- DABNI_NFL %>% rename(ID = NHC, Age = Age_at_CSF, APOE_grouped = APOE_grouped,
                                  SEX = Sex, NFL = PLASMA_NFLIGHT_S0100)
DABNI_NFL$SEX <- as.factor(DABNI_NFL$SEX)

DABNI_NFL$Group <- paste0(DABNI_NFL$APOE_grouped)
ABCDS_NFL$Group <- paste0(ABCDS_NFL$APOE_grouped)

NFL <- rbind(DABNI_NFL[, c("ID", "Age", "APOE_grouped", "SEX", "NFL", "Group")], 
             ABCDS_NFL[, c("ID", "Age", "APOE_grouped", "SEX", "NFL", "Group")])
NFL$NFL <- log(as.numeric(NFL$NFL))
ABCDS_NFL$NFL <- log(as.numeric(ABCDS_NFL$NFL))
DABNI_NFL$NFL <- log(DABNI_NFL$NFL)
hist_NFL <- plot_grouped_histogram(NFL[NFL$Group %in% c("APOE2", "APOE3", "APOE4"),], value_col = "NFL", XLAB = "log(Plasma NFL)")
hist_NFL_DABNI <- plot_grouped_histogram(DABNI_NFL[DABNI_NFL$Group %in% c("APOE2", "APOE3", "APOE4"),], value_col = "NFL", XLAB = "log(Plasma NFL)")
hist_NFL_ABCDS <- plot_grouped_histogram(ABCDS_NFL[ABCDS_NFL$Group %in% c("APOE2", "APOE3", "APOE4"),], value_col = "NFL", XLAB = "log(Plasma NFL)")


DABNI_Hipp <- DABNI_Hipp %>% rename(ID = NHC, Age = Age_at_MRI, APOE_grouped = APOE_grouped,
                                    SEX = Sex, Npos = hippo_positive)
ABCDS_Hipp <- ABCDS_Hipp %>% rename(ID = Subject.ID, Age = age_at_visit, APOE_grouped = APOE_grouped,
                                    SEX = de_gender, Npos = hippo_positive)
Hippo <- rbind(DABNI_Hipp[, c("ID", "Age", "APOE_grouped", "SEX", "hippo_Z")],
               ABCDS_Hipp[, c("ID", "Age", "APOE_grouped", "SEX", "hippo_Z")])

Hippo$Group <- paste0(Hippo$APOE_grouped)
DABNI_Hipp$Group <- paste0(DABNI_Hipp$APOE_grouped)
ABCDS_Hipp$Group <- paste0(ABCDS_Hipp$APOE_grouped)


hist_Hippo <- plot_grouped_histogram(Hippo[Hippo$Group %in% c("APOE2", "APOE3", "APOE4"),], value_col = "hippo_Z", XLAB = "Hippocampal Vol [Z]")
hist_DABNI_Hipp <- plot_grouped_histogram(DABNI_Hipp[DABNI_Hipp$Group %in% c("APOE2", "APOE3", "APOE4"),], value_col = "hippo_Z", XLAB = "Hippocampal Vol [Z]")
hist_ABCDS_Hipp <- plot_grouped_histogram(ABCDS_Hipp[ABCDS_Hipp$Group %in% c("APOE2", "APOE3", "APOE4"),], value_col = "hippo_Z", XLAB = "Hippocampal Vol [Z]")





#TODO: ADD MCRT

ABCDS_mcrt <- ABCDS[!duplicated(ABCDS$SUBJECT_LABEL), c("SUBJECT_LABEL", "trs", "APOE_grouped")]


#ABCDS - trs
DABNI_NPS <- readxl::read_xlsx("./Data_202512/NPS_20251127.xlsx")
DABNI_NPS <- DABNI_NPS[, c("NHC", "APOE", "crt_rli_rfi_total")]
DABNI_NPS$APOE <- as.factor(DABNI_NPS$APOE)

DABNI_NPS$APOE_grouped <- dplyr::recode(
  DABNI_NPS$APOE,
  "22" = "APOE2",
  "23" = "APOE2",
  "24" = "APOE4",
  "33" = "APOE3",
  "34" = "APOE4",
  "43" = "APOE4",
  "44" = "APOE4")


colnames(ABCDS_mcrt) <- c("ID", "mcrt", "Group")
colnames(DABNI_NPS) <- c("ID", "APOE", "mcrt", "Group")

mcrt <- rbind(ABCDS_mcrt[, c("ID", "mcrt", "Group")],
              DABNI_NPS[, c("ID", "mcrt", "Group")])
mcrt$mcrt <- as.numeric(mcrt$mcrt)
ABCDS_mcrt$mcrt <- as.numeric(ABCDS_mcrt$mcrt)
DABNI_NPS$mcrt <- as.numeric(DABNI_NPS$mcrt)

hist_mcrt <- plot_grouped_histogram(mcrt[mcrt$Group %in% c("APOE2", "APOE3", "APOE4"),], value_col = "mcrt", XLAB = "modified Cued Recall Test (mCRT)")
hist_DABNI_mcrt <-plot_grouped_histogram(DABNI_NPS[DABNI_NPS$Group %in% c("APOE2", "APOE3", "APOE4"),], value_col = "mcrt", XLAB = "modified Cued Recall Test (mCRT)")
hist_ABCDS_mcrt <- plot_grouped_histogram(ABCDS_mcrt[ABCDS_mcrt$Group %in% c("APOE2", "APOE3", "APOE4"),], value_col = "mcrt", XLAB = "modified Cued Recall Test (mCRT)")




layout_matrx <- rbind(c(1, 1, 1, 2, 2, 3),
                      c(1, 1, 1, 2, 2, 4),
                      c(5, 5, 5, 6, 6, 6),
                      c(5, 5, 5, 6, 6, 6),
                      c(7, 7, 8, 10, 10, 11),
                      c(7, 7, 9, 10, 10, 12),
                      c(13, 13, 14, 16, 16, 17),
                      c(13, 13, 15, 16, 16, 18))
grid.arrange(hist_CSF, hist_CL, hist_CL_ABCDS,
             hist_CL_DABNI, hist_pTau_ABCDS, hist_pTau_DABNI,
             hist_GFAP, hist_GFAP_ABCDS, hist_GFAP_DABNI,
             hist_NFL, hist_NFL_ABCDS, hist_NFL_DABNI,
             hist_Hippo, hist_ABCDS_Hipp, hist_DABNI_Hipp,
             hist_mcrt, hist_ABCDS_mcrt, hist_DABNI_mcrt, layout_matrix = layout_matrx)
graph2ppt(file = "./figures/SuppFig4_DataAvailHistograms.pptx", width = 8, height = 9)

grid.arrange(hist_Hippo, hist_ABCDS_Hipp, hist_DABNI_Hipp)
