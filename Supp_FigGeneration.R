#DABNI ONLY

res_CSF <- run_apoe_model(
  DABNI_CSF,
  apoe_col = "APOE_grouped",
  y_col = "CSFAB42_AB40",
  age_col = "Age_at_CSF",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = TRUE
) 

res_CL <- run_apoe_model_cohort(
  df = CL[CL$cohort == "DABNI",],
  apoe_col = "APOE_grouped",
  age_col = "amy_AgefromBaseline",
  y_col = "WUSTLcentiloid",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = TRUE
)

res_ptau217_DABNI <- run_apoe_model_cohort(
  df = DABNI,
  apoe_col = "APOE_grouped",
  age_col = "Age_at_CSF",
  y_col = "PLASMA_PTAU217_S0217",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = TRUE
)

# res_ptau217_ABCDS <- run_apoe_model_cohort(
#   df = pTau217,
#   apoe_col = "APOE_grouped",
#   age_col = "clinical_AgefromBaseline",
#   y_col = "pTau217",
#   fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
#   colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
#   markSig = TRUE
# )

res_GFAP <- run_apoe_model_cohort(
  df = GFAP[GFAP$cohort == "DABNI",],
  apoe_col = "APOE_grouped",
  age_col = "clinical_AgefromBaseline",
  y_col = "TESTVALUE",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = TRUE,
  k = "5"
)
res_NFL <- run_apoe_model_cohort(
  df = NFL[NFL$cohort == "DABNI",],
  apoe_col = "APOE_grouped",
  age_col = "clinical_AgefromBaseline",
  y_col = "TESTVALUE",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = TRUE
)

sig_bar_colors = c(
  "sig_APOE2" = "#7570B3", # Muted Deep Purple
  "sig_APOE3" = "#D95F02", # Burnt Copper Orange
  "sig_APOE4" = "#666666"  # Charcoal Slate Gray
)
res_hippo <- run_apoe_model_cohort(
  df = hippo[hippo$cohort == "DABNI",],
  apoe_col = "APOE_grouped",
  age_col = "mri_AgefromBaseline",
  y_col = "hippo_Z",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = TRUE,
  k = "3"
)
res_mcrt <- run_apoe_model_cohort(
  df = mcrt[mcrt$cohort == "DABNI",],
  apoe_col = "APOE_grouped",
  age_col = "clinical_AgefromBaseline",
  y_col = "trs",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = TRUE,
  k = "7",
  YLIM = c(0, 46)
)


p1 <- res_CSF[[1]]
p2 <- res_CL[[1]]
p3 <- res_GFAP[[1]]
p4 <- res_NFL[[1]]
p5 <- res_hippo[[1]]
p6 <- res_mcrt[[1]]
layout_matrix <- rbind(c(1, 2), c(3, 3), c(4, 5), c(6, 7))
lemon::grid_arrange_shared_legend(p1 + ggtitle("A."), p2 + ggtitle("B.") + scale_shape_manual(values = c(1, 4, 4)), 
                                  res_ptau217_DABNI[[1]] + scale_shape_manual(values = c(1)) + ggtitle("C."), 
                                  
                                  p3 + ggtitle("D.") + scale_shape_manual(values = c(1, 4, 4)), p4 + ggtitle("E.") + scale_shape_manual(values = c(1, 4, 4)), 
                                  p5 + ggtitle("F.") + scale_shape_manual(values = c(1, 4, 4)), p6 + ggtitle("G.") + scale_shape_manual(values = c(1, 4, 4)), layout_matrix = layout_matrix)

graph2ppt(file = "./SuppMaterials/Fig1Equivalent_DABNI_bars.pptx", width = 6, height = 11)




#DABNI ONLY

res_CSF <- run_apoe_model(
  DABNI_CSF,
  apoe_col = "APOE_grouped",
  y_col = "CSFAB42_AB40",
  age_col = "Age_at_CSF",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = FALSE
) 

res_CL <- run_apoe_model_cohort(
  df = CL[CL$cohort == "DABNI",],
  apoe_col = "APOE_grouped",
  age_col = "amy_AgefromBaseline",
  y_col = "WUSTLcentiloid",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = FALSE
)

res_ptau217_DABNI <- run_apoe_model_cohort(
  df = DABNI,
  apoe_col = "APOE_grouped",
  age_col = "Age_at_CSF",
  y_col = "PLASMA_PTAU217_S0217",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = FALSE
)

# res_ptau217_ABCDS <- run_apoe_model_cohort(
#   df = pTau217,
#   apoe_col = "APOE_grouped",
#   age_col = "clinical_AgefromBaseline",
#   y_col = "pTau217",
#   fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
#   colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
#   markSig = TRUE
# )

res_GFAP <- run_apoe_model_cohort(
  df = GFAP[GFAP$cohort == "DABNI",],
  apoe_col = "APOE_grouped",
  age_col = "clinical_AgefromBaseline",
  y_col = "TESTVALUE",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = FALSE,
  k = "5"
)
res_NFL <- run_apoe_model_cohort(
  df = NFL[NFL$cohort == "DABNI",],
  apoe_col = "APOE_grouped",
  age_col = "clinical_AgefromBaseline",
  y_col = "TESTVALUE",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = FALSE
)

sig_bar_colors = c(
  "sig_APOE2" = "#7570B3", # Muted Deep Purple
  "sig_APOE3" = "#D95F02", # Burnt Copper Orange
  "sig_APOE4" = "#666666"  # Charcoal Slate Gray
)
res_hippo <- run_apoe_model_cohort(
  df = hippo[hippo$cohort == "DABNI",],
  apoe_col = "APOE_grouped",
  age_col = "mri_AgefromBaseline",
  y_col = "hippo_Z",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = FALSE,
  k = "3"
)
res_mcrt <- run_apoe_model_cohort(
  df = mcrt[mcrt$cohort == "DABNI",],
  apoe_col = "APOE_grouped",
  age_col = "clinical_AgefromBaseline",
  y_col = "trs",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = FALSE,
  k = "7",
  YLIM = c(0, 46)
)


p1 <- res_CSF[[1]]
p2 <- res_CL[[1]]
p3 <- res_GFAP[[1]]
p4 <- res_NFL[[1]]
p5 <- res_hippo[[1]]
p6 <- res_mcrt[[1]]
layout_matrix <- rbind(c(1, 2), c(3, 3), c(4, 5), c(6, 7))
lemon::grid_arrange_shared_legend(p1 + ggtitle("A."), p2 + ggtitle("B.") + scale_shape_manual(values = c(1, 4, 4)), 
                                  res_ptau217_DABNI[[1]] + scale_shape_manual(values = c(1)) + ggtitle("C."), 
                                  
                                  p3 + ggtitle("D.") + scale_shape_manual(values = c(1, 4, 4)), p4 + ggtitle("E.") + scale_shape_manual(values = c(1, 4, 4)), 
                                  p5 + ggtitle("F.") + scale_shape_manual(values = c(1, 4, 4)), p6 + ggtitle("G.") + scale_shape_manual(values = c(1, 4, 4)), layout_matrix = layout_matrix)

graph2ppt(file = "./SuppMaterials/Fig1Equivalent_DABNI_NObars.pptx", width = 6, height = 11)


#DABNI ONLY

res_CSF <- run_apoe_model(
  DABNI_CSF,
  apoe_col = "APOE_grouped",
  y_col = "CSFAB42_AB40",
  age_col = "Age_at_CSF",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = TRUE
) 

res_CL <- run_apoe_model_cohort(
  df = CL[CL$cohort == "DABNI",],
  apoe_col = "APOE_grouped",
  age_col = "amy_AgefromBaseline",
  y_col = "WUSTLcentiloid",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = TRUE
)

res_ptau217_DABNI <- run_apoe_model_cohort(
  df = DABNI,
  apoe_col = "APOE_grouped",
  age_col = "Age_at_CSF",
  y_col = "PLASMA_PTAU217_S0217",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = TRUE
)

# res_ptau217_ABCDS <- run_apoe_model_cohort(
#   df = pTau217,
#   apoe_col = "APOE_grouped",
#   age_col = "clinical_AgefromBaseline",
#   y_col = "pTau217",
#   fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
#   colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
#   markSig = TRUE
# )

res_GFAP <- run_apoe_model_cohort(
  df = GFAP[GFAP$cohort == "DABNI",],
  apoe_col = "APOE_grouped",
  age_col = "clinical_AgefromBaseline",
  y_col = "TESTVALUE",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = TRUE,
  k = "5"
)
res_NFL <- run_apoe_model_cohort(
  df = NFL[NFL$cohort == "DABNI",],
  apoe_col = "APOE_grouped",
  age_col = "clinical_AgefromBaseline",
  y_col = "TESTVALUE",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = TRUE
)

sig_bar_colors = c(
  "sig_APOE2" = "#7570B3", # Muted Deep Purple
  "sig_APOE3" = "#D95F02", # Burnt Copper Orange
  "sig_APOE4" = "#666666"  # Charcoal Slate Gray
)
res_hippo <- run_apoe_model_cohort(
  df = hippo[hippo$cohort == "DABNI",],
  apoe_col = "APOE_grouped",
  age_col = "mri_AgefromBaseline",
  y_col = "hippo_Z",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = TRUE,
  k = "3"
)
res_mcrt <- run_apoe_model_cohort(
  df = mcrt[mcrt$cohort == "DABNI",],
  apoe_col = "APOE_grouped",
  age_col = "clinical_AgefromBaseline",
  y_col = "trs",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = TRUE,
  k = "7",
  YLIM = c(0, 46)
)


p1 <- res_CSF[[1]]
p2 <- res_CL[[1]]
p3 <- res_GFAP[[1]]
p4 <- res_NFL[[1]]
p5 <- res_hippo[[1]]
p6 <- res_mcrt[[1]]
layout_matrix <- rbind(c(1, 2), c(3, 3), c(4, 5), c(6, 7))
lemon::grid_arrange_shared_legend(p1 + ggtitle("A."), p2 + ggtitle("B.") + scale_shape_manual(values = c(1, 4, 4)), 
                                  res_ptau217_DABNI[[1]] + scale_shape_manual(values = c(1)) + ggtitle("C."), 
                                  
                                  p3 + ggtitle("D.") + scale_shape_manual(values = c(1, 4, 4)), p4 + ggtitle("E.") + scale_shape_manual(values = c(1, 4, 4)), 
                                  p5 + ggtitle("F.") + scale_shape_manual(values = c(1, 4, 4)), p6 + ggtitle("G.") + scale_shape_manual(values = c(1, 4, 4)), layout_matrix = layout_matrix)

graph2ppt(file = "./SuppMaterials/Fig1Equivalent_DABNI_bars.pptx", width = 6, height = 11)




#ABCDS ONLY


res_CL <- run_apoe_model_cohort(
  df = CL[CL$cohort == "ABCDS",],
  apoe_col = "APOE_grouped",
  age_col = "amy_AgefromBaseline",
  y_col = "WUSTLcentiloid",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = FALSE
)


res_ptau217_ABCDS <- run_apoe_model_cohort(
  df = pTau217,
  apoe_col = "APOE_grouped",
  age_col = "clinical_AgefromBaseline",
  y_col = "pTau217",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = FALSE
)

res_GFAP <- run_apoe_model_cohort(
  df = GFAP[GFAP$cohort == "ABCDS",],
  apoe_col = "APOE_grouped",
  age_col = "clinical_AgefromBaseline",
  y_col = "TESTVALUE",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = FALSE,
  k = "5"
)
res_NFL <- run_apoe_model_cohort(
  df = NFL[NFL$cohort == "ABCDS",],
  apoe_col = "APOE_grouped",
  age_col = "clinical_AgefromBaseline",
  y_col = "TESTVALUE",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = FALSE
)

sig_bar_colors = c(
  "sig_APOE2" = "#7570B3", # Muted Deep Purple
  "sig_APOE3" = "#D95F02", # Burnt Copper Orange
  "sig_APOE4" = "#666666"  # Charcoal Slate Gray
)
res_hippo <- run_apoe_model_cohort(
  df = hippo[hippo$cohort == "ABCDS",],
  apoe_col = "APOE_grouped",
  age_col = "mri_AgefromBaseline",
  y_col = "hippo_Z",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = FALSE,
  k = "3"
)
res_mcrt <- run_apoe_model_cohort(
  df = mcrt[mcrt$cohort == "ABCDS",],
  apoe_col = "APOE_grouped",
  age_col = "clinical_AgefromBaseline",
  y_col = "trs",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = FALSE,
  k = "7",
  YLIM = c(0, 46)
)


p2 <- res_CL[[1]]
p3 <- res_GFAP[[1]]
p4 <- res_NFL[[1]]
p5 <- res_hippo[[1]]
p6 <- res_mcrt[[1]]
lemon::grid_arrange_shared_legend( p2 + ggtitle("A.") + scale_shape_manual(values = c(1)), 
                                  res_ptau217_ABCDS[[1]] + scale_shape_manual(values = c(4)) + ggtitle("B."), 
                                  
                                  p3 + ggtitle("C.") + scale_shape_manual(values = c( 4)), p4 + ggtitle("D.") + scale_shape_manual(values = c( 4)), 
                                  p5 + ggtitle("E.") + scale_shape_manual(values = c(4)), p6 + ggtitle("F.") + scale_shape_manual(values = c( 4)), 
                                  nrow = 3, ncol = 2)

graph2ppt(file = "./SuppMaterials/Fig1Equivalent_ABCDS_NObars.pptx", width = 6, height = 11)




res_CL <- run_apoe_model_cohort(
  df = CL[CL$cohort == "ABCDS",],
  apoe_col = "APOE_grouped",
  age_col = "amy_AgefromBaseline",
  y_col = "WUSTLcentiloid",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = TRUE
)


res_ptau217_ABCDS <- run_apoe_model_cohort(
  df = pTau217,
  apoe_col = "APOE_grouped",
  age_col = "clinical_AgefromBaseline",
  y_col = "pTau217",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig =TRUE
)

res_GFAP <- run_apoe_model_cohort(
  df = GFAP[GFAP$cohort == "ABCDS",],
  apoe_col = "APOE_grouped",
  age_col = "clinical_AgefromBaseline",
  y_col = "TESTVALUE",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig =TRUE,
  k = "5"
)
res_NFL <- run_apoe_model_cohort(
  df = NFL[NFL$cohort == "ABCDS",],
  apoe_col = "APOE_grouped",
  age_col = "clinical_AgefromBaseline",
  y_col = "TESTVALUE",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig =TRUE
)

sig_bar_colors = c(
  "sig_APOE2" = "#7570B3", # Muted Deep Purple
  "sig_APOE3" = "#D95F02", # Burnt Copper Orange
  "sig_APOE4" = "#666666"  # Charcoal Slate Gray
)
res_hippo <- run_apoe_model_cohort(
  df = hippo[hippo$cohort == "ABCDS",],
  apoe_col = "APOE_grouped",
  age_col = "mri_AgefromBaseline",
  y_col = "hippo_Z",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = TRUE,
  k = "3"
)
res_mcrt <- run_apoe_model_cohort(
  df = mcrt[mcrt$cohort == "ABCDS",],
  apoe_col = "APOE_grouped",
  age_col = "clinical_AgefromBaseline",
  y_col = "trs",
  fill_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  colour_values = c("#3399FF", "#33CC33", "#FF4D4D"),
  markSig = TRUE,
  k = "7",
  YLIM = c(0, 46)
)


p2 <- res_CL[[1]]
p3 <- res_GFAP[[1]]
p4 <- res_NFL[[1]]
p5 <- res_hippo[[1]]
p6 <- res_mcrt[[1]]
lemon::grid_arrange_shared_legend( p2 + ggtitle("A.") + scale_shape_manual(values = c(1)), 
                                   res_ptau217_ABCDS[[1]] + scale_shape_manual(values = c(4)) + ggtitle("B."), 
                                   
                                   p3 + ggtitle("C.") + scale_shape_manual(values = c( 4)), p4 + ggtitle("D.") + scale_shape_manual(values = c( 4)), 
                                   p5 + ggtitle("E.") + scale_shape_manual(values = c(4)), p6 + ggtitle("F.") + scale_shape_manual(values = c( 4)), 
                                   nrow = 3, ncol = 2)

graph2ppt(file = "./SuppMaterials/Fig1Equivalent_ABCDS_bars.pptx", width = 6, height = 11)