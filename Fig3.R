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
dx <- read.csv(".././ABCDS_LONI_20250528/Consensus_28May2025.csv")
dx <- dx[, c("subject_label", "event_sequence", "consensus_dx")]
dx <- merge(dx, ABCDS_apoe[, c("subject_label", "APOE_grouped")], by = "subject_label", all = FALSE)
dx <- merge(dx, demogs, by = "subject_label", all.x = TRUE, all.y = FALSE)
dx <- dx[!is.na(dx$de_gender),]

dx <- dx[order(dx$subject_label, dx$event_sequence),]
dx <- dx[!duplicated(dx$subject_label),]

dx <- merge(dx, visit_latency[, c("subject_label", "event_sequence", "age_at_visit")],
               by.x = c("subject_label", "event_sequence"), by.y = c("subject_label", "event_sequence"), all = FALSE)
colnames(dx) <- c("ID", "visit", "dx", "APOE", "sex", "age")

DABNI_dx <- readxl::read_xlsx("./Data_202512/NPS_20251127.xlsx")
DABNI_dx <- DABNI_dx[, c("NHC", "Sex", "Age_NPS", "Diag", "APOE")]
DABNI_dx <- DABNI_dx[order(DABNI_dx$NHC, DABNI_dx$Age_NPS),]
DABNI_dx <- DABNI_dx[!duplicated(DABNI_dx$NHC),]


colnames(DABNI_dx) <- c("ID", "sex",  "age", "dx", "APOE")


dx <- rbind(dx[, c("ID", "sex",  "age", "dx", "APOE")], DABNI_dx)
dx$dx <- as.factor(dx$dx)
levels(dx$dx) <- c("Asymptomatic", "Symptomatic", "Symptomatic", "noConsensus", "Asymptomatic",
                   "Symptomatic", "Symptomatic", "noConsensus")
dx$APOE <- as.factor(dx$APOE)
levels(dx$APOE) <- c("APOE2", "APOE2", "APOE4", "APOE3", "APOE4", "APOE4", "APOE4", "APOE2", "APOE3", "APOE4")

dx$age_group <- cut(dx$age, breaks = c(0, 35, 40, 45, 50, 55, 100))
levels(dx$age_group) <- c("< 35 Years Old", "35 – 40 Years Old", "40 – 45 Years Old",
                         "45 – 50 Years Old", "50 – 55 Years Old", "> 55 Years Old")
dx <- dx[!(dx$dx == "noConsensus"),]
library(ggplot2)
library(dplyr)
library(ggpubr)
library(rstatix)

# ------------------------------------------------------------------------------
# 1. Data Preparation
# Convert dx to numeric (0 = Asymptomatic, 1 = Symptomatic) for ANOVA
# ------------------------------------------------------------------------------
plot_df <- dx %>%
  mutate(
    APOE = factor(APOE, levels = c("APOE2", "APOE3", "APOE4")),
    age_group = factor(age_group, levels = c(
      "< 35 Years Old", "35 – 40 Years Old", "40 – 45 Years Old",
      "45 – 50 Years Old", "50 – 55 Years Old", "> 55 Years Old"
    )),
    dx = factor(dx, levels = c("Asymptomatic", "Symptomatic")),
    
    # Numeric binary outcome for ANOVA model
    dx_num = ifelse(dx == "Symptomatic", 1, 0),
    
    # Custom legend categories
    fill_cat = factor(
      paste(APOE, "—", dx),
      levels = c(
        "APOE2 — Asymptomatic", "APOE3 — Asymptomatic", "APOE4 — Asymptomatic",
        "APOE2 — Symptomatic",  "APOE3 — Symptomatic",  "APOE4 — Symptomatic"
      )
    )
  )

# Sample counts per APOE & age group
n_counts <- plot_df %>%
  group_by(age_group, APOE) %>%
  summarize(N = n(), .groups = "drop") %>%
  mutate(label = paste0("N = ", N))

# 1. Ensure dx_num exists in plot_df

# ------------------------------------------------------------------------------
# 2. Automated ANOVA + Post-Hoc Pairwise Comparisons
# ------------------------------------------------------------------------------
p_value_df <- plot_df %>%
  group_by(age_group) %>%
  group_modify(~ pairwise_fisher_test(table(.x$dx, .x$APOE), p.adjust.method = "bonferroni")) %>%
  mutate(
    p.adj.signif = case_when(
      p.adj < 0.001 ~ "p<0.001",
      TRUE ~ paste0("p=", round(p.adj, 2))
    ),
    y.position = case_when(
      (group1 == "APOE2" & group2 == "APOE3") | (group1 == "APOE3" & group2 == "APOE2") ~ 1.03,
      (group1 == "APOE2" & group2 == "APOE4") | (group1 == "APOE4" & group2 == "APOE2") ~ 1.22,
      (group1 == "APOE3" & group2 == "APOE4") | (group1 == "APOE4" & group2 == "APOE3") ~ 1.12,
      TRUE ~ 1.10
    )
  )

# Color Palette
custom_colors <- c(
  "APOE2 — Asymptomatic" = "#d0e5ff", "APOE2 — Symptomatic" = "#2b90f5",
  "APOE3 — Asymptomatic" = "#dbf0db", "APOE3 — Symptomatic" = "#32cd32",
  "APOE4 — Asymptomatic" = "#ffd6d6", "APOE4 — Symptomatic" = "#ff3b30"
)

# ------------------------------------------------------------------------------
# 3. Plotting
# ------------------------------------------------------------------------------
ggplot(plot_df, aes(x = APOE, fill = fill_cat)) +
  geom_bar(position = "fill", color = "black", width = 0.8) +
  
  # Sample size text (N = xx)
  geom_text(
    data = n_counts,
    aes(x = APOE, y = 0.08, label = label),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 3.8
  ) +
  
  # Dynamic ANOVA post-hoc p-value brackets
  stat_pvalue_manual(
    p_value_df,
    label = "p.adj.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    group.by = "age_group",
    tip.length = 0,
    label.size = 3.5
  ) +
  
  facet_wrap(~ age_group, ncol = 3, scales = "free_x") +
  scale_fill_manual(values = custom_colors, name = NULL) +
  scale_y_continuous(
    labels = scales::percent_format(),
    breaks = c(0, 0.5, 1.0),
    limits = c(0, 1.35),
    expand = c(0, 0)
  ) +
  labs(x = "Group", y = "Proportion") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey92"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey85", color = "black"),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    legend.box = "horizontal",
    axis.text = element_text(color = "black")
  ) +
  guides(fill = guide_legend(ncol = 3, byrow = TRUE))
