

evaluate_gam_k_selection <- function(df, y_col, age_col, apoe_col, cohort = NULL, sex_col = NULL, intellectual_disability_col = NULL, k_range = 3:10) {
  
  # Ensure core factors/numerics are set
  df[[apoe_col]] <- as.factor(df[[apoe_col]])
  df[[age_col]]   <- as.numeric(df[[age_col]])
  
  # Conditionally convert covariates to factors if provided
  if (!is.null(cohort) && cohort %in% colnames(df)) {
    df[[cohort]] <- as.factor(df[[cohort]])
  }
  if (!is.null(sex_col) && sex_col %in% colnames(df)) {
    df[[sex_col]] <- as.factor(df[[sex_col]])
  }
  
  # Clean ID column factor mapping
  if (!is.null(intellectual_disability_col) && intellectual_disability_col %in% colnames(df)) {
    df[[intellectual_disability_col]] <- factor(df[[intellectual_disability_col]])

  }
  
  
  results_list <- list()
  
  for (k in k_range) {
    # 1. Build the formula string dynamically
    formula_str <- paste0(
      y_col, " ~ ", apoe_col, " + s(", age_col, ", k = ",k,")",
      # APOE-specific deviation
      " + s(",age_col, ", by = ", apoe_col,", k = ",
      k, ")"
    )    
    if (!is.null(cohort) && cohort %in% colnames(df)) {
      formula_str <- paste0(formula_str, " + ", cohort)
    }
    if (!is.null(sex_col) && sex_col %in% colnames(df)) {
      formula_str <- paste0(formula_str, " + ", sex_col)
    }
    
    if (!is.null(intellectual_disability_col) && intellectual_disability_col %in% colnames(df)) {
      formula_str <- paste0(formula_str, " + ", intellectual_disability_col)
    }
    
    formula_gam <- as.formula(formula_str)
    
    # Fit the model for the current k
    mod <- gam(formula_gam, method = "REML", data = df)
    
    # 2. Extract k-check diagnostics
    # k.check returns a matrix with columns: k', edf, k-index, and p-value
    kc <- k.check(mod)
    kc_df <- as.data.frame(kc)
    
    # Add metadata
    kc_df$k_tested <- k
    kc_df$AIC <- AIC(mod)
    kc_df$smooth_term <- rownames(kc)
    
    results_list[[as.character(k)]] <- kc_df
  }
  
  # Combine results
  full_diagnostics <- bind_rows(results_list)
  
  # 3. Summarize to find the "Optimal" point
  # We look for where the p-value stops being significant (p > 0.05)
  # and where AIC stabilizes.
  summary_stats <- full_diagnostics %>%
    group_by(k_tested) %>%
    summarise(
      mean_k_index = mean(`k-index`, na.rm = TRUE),
      min_p_value = min(`p-value`, na.rm = TRUE),
      total_edf = sum(edf, na.rm = TRUE),
      AIC = first(AIC),
      .groups = 'drop'
    )
  
  # 4. Generate Diagnostic Plots
  p1 <- ggplot(summary_stats, aes(x = k_tested, y = min_p_value)) +
    geom_line() + geom_point() +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    labs(title = "K-Check p-values (Min across groups)", 
         subtitle = "Values above red line suggest k is sufficient",
         y = "p-value", x = "k value") +
    theme_bw()
  
  p2 <- ggplot(summary_stats, aes(x = k_tested, y = AIC)) +
    geom_line() + geom_point() +
    labs(title = "AIC by Basis Dimension (k)", 
         y = "AIC", x = "k value") +
    theme_bw()
  
  return(list(
    full_data = full_diagnostics,
    summary = summary_stats,
    plots = list(p_val_plot = p1, aic_plot = p2)
  ))
}






#' Export GAM Summary to Publication-Ready Excel Table
#'
#' @param gam_summary A summary.gam object or a fitted gam model object.
#' @param file_path Character string indicating where to save the .xlsx file.
#' @param model_label Optional title or description for the table header.
#'
export_gam_summary_to_excel <- function(gam_summary, file_path, model_label = "Supplemental Table: GAM Model Statistics") {
  
  # Dependencies
  if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Please install 'openxlsx' package.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install 'dplyr' package.")
  
  # Handle cases where user passed the raw model object instead of the summary
  if (inherits(gam_summary, "gam")) {
    gam_summary <- summary(gam_summary)
  }
  
  # ---------------------------------------------------------
  # 1. Parse Parametric (Fixed) Effects
  # ---------------------------------------------------------
  p_df <- as.data.frame(gam_summary$p.table)
  p_df$Term <- rownames(p_df)
  p_df <- p_df[, c("Term", "Estimate", "Std. Error", "t value", "Pr(>|t|)")]
  colnames(p_df) <- c("Variable / Predictor", "Estimate", "Std. Error", "t-statistic", "p-value")
  
  # Format metrics to standard clinical journal precisions
  p_df$Estimate     <- sprintf("%.3f", p_df$Estimate)
  p_df$`Std. Error` <- sprintf("%.3f", p_df$`Std. Error`)
  p_df$`t-statistic`<- sprintf("%.2f", p_df$`t-statistic`)
  p_df$`p-value`    <- ifelse(as.numeric(p_df$`p-value`) < 0.001, "<0.001", sprintf("%.3f", as.numeric(p_df$`p-value`)))
  
  # ---------------------------------------------------------
  # 2. Parse Smooth Non-Linear Effects
  # ---------------------------------------------------------
  s_df <- as.data.frame(gam_summary$s.table)
  s_df$Term <- rownames(s_df)
  s_df <- s_df[, c("Term", "edf", "Ref.df", "F", "p-value")]
  colnames(s_df) <- c("Smooth Term Spline", "Estimated DF (edf)", "Reference DF", "F-statistic", "p-value")
  
  s_df$`Estimated DF (edf)` <- sprintf("%.2f", s_df$`Estimated DF (edf)`)
  s_df$`Reference DF`       <- sprintf("%.2f", s_df$`Reference DF`)
  s_df$`F-statistic`        <- sprintf("%.2f", s_df$`F-statistic`)
  s_df$`p-value`            <- ifelse(as.numeric(s_df$`p-value`) < 0.001, "<0.001", sprintf("%.3f", as.numeric(s_df$`p-value`)))
  
  # ---------------------------------------------------------
  # 3. Parse Global Model Diagnostics
  # ---------------------------------------------------------
  formula_string <- paste(deparse(gam_summary$formula), collapse = " ")
  
  global_df <- data.frame(
    Diagnostic = c("Sample Size (N)", "Adjusted R-squared", "Deviance Explained", "Smoothing Selection Method"),
    Value = c(
      as.character(gam_summary$n),
      sprintf("%.3f", gam_summary$r.sq),
      sprintf("%.1f%%", gam_summary$dev.expl * 100),
      gsub("-", "", as.character(gam_summary$method))
    ),
    stringsAsFactors = FALSE
  )
  
  # ---------------------------------------------------------
  # 4. Build Workbook & Professional Layout Styling
  # ---------------------------------------------------------
  wb <- openxlsx::createWorkbook()
  sheet_name <- "GAM_Summary_Table"
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::showGridLines(wb, sheet_name, show = TRUE)
  
  # Define styles
  title_style   <- openxlsx::createStyle(fontName = "Arial", fontSize = 14, textDecoration = "bold")
  formula_style <- openxlsx::createStyle(fontName = "Arial", fontSize = 10, textDecoration = "italic", fontColour = "#555555")
  section_style <- openxlsx::createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold", fgFill = "#F2F2F2")
  header_style  <- openxlsx::createStyle(fontName = "Arial", fontSize = 10, textDecoration = "bold", 
                                         border = c("top", "bottom"), borderColour = "#000000", borderStyle = "thin",
                                         halign = "left")
  data_style    <- openxlsx::createStyle(fontName = "Arial", fontSize = 10, halign = "left")
  bold_lbl_style<- openxlsx::createStyle(fontName = "Arial", fontSize = 10, textDecoration = "bold")
  bottom_style  <- openxlsx::createStyle(border = "bottom", borderColour = "#000000", borderStyle = "thin")
  
  # Write elements sequentially with explicit spacing control
  curr_row <- 1
  openxlsx::writeData(wb, sheet_name, model_label, startRow = curr_row, startCol = 1)
  openxlsx::addStyle(wb, sheet_name, title_style, rows = curr_row, cols = 1)
  
  curr_row <- curr_row + 1
  openxlsx::writeData(wb, sheet_name, paste("Model Specification Formula:", formula_string), startRow = curr_row, startCol = 1)
  openxlsx::addStyle(wb, sheet_name, formula_style, rows = curr_row, cols = 1)
  
  # Block 1: Parametric
  curr_row <- curr_row + 2
  openxlsx::writeData(wb, sheet_name, "1. Parametric Coefficients (Fixed Linear Effects)", startRow = curr_row, startCol = 1)
  openxlsx::addStyle(wb, sheet_name, section_style, rows = curr_row, cols = 1:5, gridExpand = TRUE)
  
  curr_row <- curr_row + 1
  p_start <- curr_row
  openxlsx::writeData(wb, sheet_name, p_df, startRow = curr_row, startCol = 1, headerStyle = header_style)
  curr_row <- curr_row + nrow(p_df) + 1
  openxlsx::addStyle(wb, sheet_name, data_style, rows = (p_start+1):(curr_row-1), cols = 1:5, gridExpand = TRUE)
  openxlsx::addStyle(wb, sheet_name, bottom_style, rows = curr_row - 1, cols = 1:5)
  
  # Block 2: Smooth terms
  curr_row <- curr_row + 1
  openxlsx::writeData(wb, sheet_name, "2. Smooth Terms (Non-linear Spline Fit Assessment)", startRow = curr_row, startCol = 1)
  openxlsx::addStyle(wb, sheet_name, section_style, rows = curr_row, cols = 1:5, gridExpand = TRUE)
  
  curr_row <- curr_row + 1
  s_start <- curr_row
  openxlsx::writeData(wb, sheet_name, s_df, startRow = curr_row, startCol = 1, headerStyle = header_style)
  curr_row <- curr_row + nrow(s_df) + 1
  openxlsx::addStyle(wb, sheet_name, data_style, rows = (s_start+1):(curr_row-1), cols = 1:5, gridExpand = TRUE)
  openxlsx::addStyle(wb, sheet_name, bottom_style, rows = curr_row - 1, cols = 1:5)
  
  # Block 3: Model performance diagnostics
  curr_row <- curr_row + 1
  openxlsx::writeData(wb, sheet_name, "3. Global Fit & Diagnostic Performance Metrics", startRow = curr_row, startCol = 1)
  openxlsx::addStyle(wb, sheet_name, section_style, rows = curr_row, cols = 1:2, gridExpand = TRUE)
  
  curr_row <- curr_row + 1
  g_start <- curr_row
  openxlsx::writeData(wb, sheet_name, global_df, startRow = curr_row, startCol = 1, colNames = FALSE)
  curr_row <- curr_row + nrow(global_df)
  openxlsx::addStyle(wb, sheet_name, bold_lbl_style, rows = g_start:(curr_row-1), cols = 1)
  openxlsx::addStyle(wb, sheet_name, data_style, rows = g_start:(curr_row-1), cols = 2)
  openxlsx::addStyle(wb, sheet_name, bottom_style, rows = curr_row - 1, cols = 1:2)
  
  # Auto-adjust column widths safely
  openxlsx::setColWidths(wb, sheet_name, cols = 1:5, widths = "auto")
  
  # Save out workbook
  openxlsx::saveWorkbook(wb, file_path, overwrite = TRUE)
  message(paste("Successfully saved styled GAM summary to:", file_path))
}



compare_sexes_by_apoe <- function(
    df,
    apoe_col = "APOE_grouped",
    sex_col = "de_gender",                             
    y_col = "WUSTLcentiloid",
    age_col = "clinical_AgefromBaseline",
    cohort = "cohort",
    fill_values = c("Female" = "#E41A1C", "Male" = "#377EB8"),   
    colour_values = c("Female" = "#E41A1C", "Male" = "#377EB8"),
    sig_bar_colors = c("sig_APOE2" = "purple", "sig_APOE3" = "orange", "sig_APOE4" = "darkred"),
    markSig = FALSE,
    k_val = 4,
    reference_cohort = NULL,
    XLIM = c(18, 72),
    YLIM = FALSE
) {
  library(mgcv)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(ggplot2)
  
  # ---- Fix column types ----
  df[[apoe_col]]  <- factor(df[[apoe_col]], levels = c("APOE2", "APOE3", "APOE4"))
  df[[sex_col]]   <- factor(df[[sex_col]])
  df[[age_col]]   <- as.numeric(df[[age_col]])
  
  # Generate explicit interaction tracking variable for the smooth split
  df$apoe_sex <- interaction(df[[apoe_col]], df[[sex_col]], sep = "_")
  
  if (!is.null(cohort) && cohort %in% colnames(df)) {
    df[[cohort]]    <- factor(df[[cohort]])
    if (is.null(reference_cohort)) {
      reference_cohort <- levels(df[[cohort]])[1]
    }
  }
  
  # ---- Determine Group-Specific Age Bounds (N >= 5 at limits per Sex/APOE group) ----
  bounds_df <- df %>%
    group_by(!!sym(apoe_col), !!sym(sex_col)) %>%
    summarise(
      min_valid_age = if(n() >= 5) sort(!!sym(age_col))[5] else min(!!sym(age_col), na.rm=TRUE),
      max_valid_age = if(n() >= 5) sort(!!sym(age_col), decreasing = TRUE)[5] else max(!!sym(age_col), na.rm=TRUE),
      .groups = 'drop'
    )
  
  # ---- Fit GAM with Sex-by-APOE Interaction Trajectories ----
  formula_str <- paste0(y_col, " ~ apoe_sex + s(", age_col, ", by = apoe_sex, k = ", k_val, ")",
                        " + s(", age_col, ", k = ",k_val,")")
  if (!is.null(cohort) && cohort %in% colnames(df)) {
    formula_str <- paste0(formula_str, " + ", cohort)
  }
  
  formula_gam <- as.formula(formula_str)
  mod <- gam(formula_gam, method = "REML", data = df)
  
  # ---- Prediction Grid Setup ----
  grid_list <- list(
    seq(min(df[[age_col]], na.rm = TRUE), max(df[[age_col]], na.rm = TRUE), by = 0.1),
    levels(df[[apoe_col]]),
    levels(df[[sex_col]])
  )
  grid_names <- c(age_col, apoe_col, sex_col)
  
  if (!is.null(cohort) && cohort %in% colnames(df)) {
    grid_list <- c(grid_list, list(reference_cohort))
    grid_names <- c(grid_names, cohort)
  }
  
  mod_df <- expand.grid(grid_list, stringsAsFactors = FALSE)
  colnames(mod_df) <- grid_names
  mod_df$apoe_sex <- interaction(mod_df[[apoe_col]], mod_df[[sex_col]], sep = "_")
  
  if (!is.null(cohort) && cohort %in% colnames(df)) {
    mod_df[[cohort]] <- factor(mod_df[[cohort]], levels = levels(df[[cohort]]))
  }
  
  # ---- Generate Estimates & Compute Bonferroni Critical Value ----
  tmp <- predict(mod, newdata = mod_df, se.fit = TRUE)
  mod_df$forecast <- tmp$fit
  mod_df$se        <- tmp$se.fit
  
  # Bonferroni adjustment: 3 planned tests (M vs F within APOE2, APOE3, APOE4)
  n_comparisons <- 3
  alpha_adj     <- 0.05 / n_comparisons
  z_crit        <- qnorm(1 - alpha_adj / 2) 
  
  mod_df$min      <- mod_df$forecast - z_crit * mod_df$se
  mod_df$max      <- mod_df$forecast + z_crit * mod_df$se
  
  # ---- Apply Constraints ----
  mod_df <- mod_df %>%
    left_join(bounds_df, by = c(apoe_col, sex_col)) %>%
    filter(!!sym(age_col) >= min_valid_age & !!sym(age_col) <= max_valid_age)
  
  if (is.numeric(YLIM) && length(YLIM) == 2) {
    mod_df$min <- pmax(mod_df$min, YLIM[1])
    mod_df$max <- pmin(mod_df$max, YLIM[2])
  }
  
  # ---- Pivot Wider to Compare Sexes within Genotypes Safely ----
  mod_df_wide <- mod_df %>%
    dplyr::select(all_of(c(age_col, apoe_col, sex_col, "forecast", "se"))) %>%
    pivot_wider(names_from = all_of(sex_col), values_from = c("forecast", "se"))
  
  sex_levels <- levels(df[[sex_col]])
  f_level <- sex_levels[1] 
  m_level <- sex_levels[2] 
  
  f_fore_col <- paste0("forecast_", f_level)
  m_fore_col <- paste0("forecast_", m_level)
  f_se_col   <- paste0("se_", f_level)
  m_se_col   <- paste0("se_", m_level)
  
  if (f_fore_col %in% names(mod_df_wide) && m_fore_col %in% names(mod_df_wide)) {
    mod_df_wide <- mod_df_wide %>%
      mutate(
        diff_M_F = .data[[m_fore_col]] - .data[[f_fore_col]],
        se_M_F   = sqrt(.data[[m_se_col]]^2 + .data[[f_se_col]]^2),
        diff_low = diff_M_F - z_crit * se_M_F,
        diff_hi  = diff_M_F + z_crit * se_M_F,
        sig_diff = diff_low > 0 | diff_hi < 0
      )
    
    # Extract continuous significant age windows per individual genotype panel safely
    significance_list <- list()
    for (apoe_val in levels(df[[apoe_col]])) {
      sub_df <- mod_df_wide %>% filter(.data[[apoe_col]] == apoe_val)
      if (nrow(sub_df) > 0) {
        ranges <- get_sig_ranges(as.data.table(sub_df), "sig_diff", age_col)
        if (!is.null(ranges) && nrow(ranges) > 0) {
          ranges[, contrast := paste0("sig_", apoe_val)]
          ranges[, (apoe_col) := apoe_val] 
          significance_list[[apoe_val]] <- ranges
        }
      }
    }
    significance_summary <- rbindlist(significance_list, fill = TRUE)
  } else {
    significance_summary <- data.table()
  }
  
  # ---- Build Stratified Facet Plot ----
  plot_aes <- aes(
    x = !!sym(age_col),
    y = !!sym(y_col),
    colour = !!sym(sex_col),
    group = !!sym(sex_col)
  )
  
  if (!is.null(cohort) && cohort %in% colnames(df)) {
    plot_aes$shape <- sym(cohort)
  }
  
  p <- ggplot(df, plot_aes) +
    geom_point(alpha = 0.3) +
    geom_ribbon(
      data = mod_df,
      aes(
        x = !!sym(age_col),
        y = forecast,
        ymin = min,
        ymax = max,
        fill = !!sym(sex_col)
      ),
      alpha = 0.25,
      inherit.aes = FALSE
    ) +
    geom_line(
      data = mod_df,
      aes(
        x = !!sym(age_col),
        y = forecast,
        group = !!sym(sex_col)
      ),
      linewidth = 1
    ) +
    facet_wrap(vars(!!sym(apoe_col))) +                            
    theme_bw() +
    theme(legend.position = "bottom", strip.background = element_blank()) +
    scale_fill_manual(values = fill_values) +
    scale_x_continuous(breaks = c(30, 50, 70))
  
  if (!is.null(cohort) && cohort %in% colnames(df)) {
    p <- p + scale_shape_manual(values = c(1, 4))
  }
  
  # ---- Add Correctly Faceted Significance Bars ----
  if (nrow(significance_summary) > 0 & markSig == TRUE) {
    y_range_max <- if(is.numeric(YLIM)) YLIM[2] else max(df[[y_col]], na.rm = TRUE)
    y_range_min <- if(is.numeric(YLIM)) YLIM[1] else min(df[[y_col]], na.rm = TRUE)
    y_buffer    <- (y_range_max - y_range_min) * 0.05 
    
    significance_summary[, y := y_range_max + y_buffer]
    
    p <- p + geom_segment(
      data = significance_summary,
      aes(
        x = start_age,
        xend = end_age,
        y = y,
        yend = y,
        colour = contrast
      ),
      inherit.aes = FALSE,
      linewidth = 2
    ) +
      scale_colour_manual(
        values = c(colour_values, 
                   sig_APOE2 = sig_bar_colors["sig_APOE2"],
                   sig_APOE3 = sig_bar_colors["sig_APOE3"],
                   sig_APOE4 = sig_bar_colors["sig_APOE4"]),
        guide = guide_legend(override.aes = list(shape = NA))
      )
  } else {
    p <- p + scale_colour_manual(values = colour_values) # Fixed typo here (<- instead of :=)
  }
  
  # ---- Unified Viewport Windowing ----
  final_xlim <- if (is.numeric(XLIM) && length(XLIM) == 2) XLIM else c(min(df[[age_col]], na.rm = TRUE), max(df[[age_col]], na.rm = TRUE))
  final_ylim <- if (is.numeric(YLIM) && length(YLIM) == 2) YLIM else NULL
  
  p <- p + coord_cartesian(xlim = final_xlim, ylim = final_ylim)
  
  return(list(
    plot = p,
    predictions = mod_df,
    significance_summary = significance_summary,
    model_obj = mod,
    z_value_used = z_crit
  ))
}




run_apoe_sex_model_cohort <- function(
    df,
    apoe_col = "APOE_grouped",
    sex_col = "de_gender",                                 # Made mandatory to support the 6-group interaction
    y_col = "WUSTLcentiloid",
    age_col = "clinical_AgefromBaseline",
    intellectual_disability_col = "Intellectual_Disability",
    cohort = "cohort",
    colour_values,                                        # 3 colors mapped to APOE2, APOE3, APOE4 names
    k_val = 4,
    reference_cohort = NULL,
    reference_sex = NULL,
    reference_id = NULL,
    XLIM = c(18, 72),                                     
    YLIM = FALSE
) {
  
  library(mgcv)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(ggplot2)
  
  # ---- Fix column types ----
  df[[apoe_col]]  <- factor(df[[apoe_col]], levels = c("APOE2", "APOE3", "APOE4"))
  df[[sex_col]]   <- factor(df[[sex_col]]) # Assumes levels are ordered such that Female/Male or F/M map cleanly
  df[[age_col]]   <- as.numeric(df[[age_col]])
  
  # Generate explicit interaction variable to fit individual smooths safely
  df$apoe_sex <- interaction(df[[apoe_col]], df[[sex_col]], sep = "_")
  
  if (!is.null(cohort) && cohort %in% colnames(df)) {
    df[[cohort]]    <- factor(df[[cohort]])
    if (is.null(reference_cohort)) {
      reference_cohort <- levels(df[[cohort]])[1]
    }
  }
  if (!is.null(sex_col) && sex_col %in% colnames(df)) {
    df[[sex_col]]    <- factor(df[[sex_col]])
    if (is.null(reference_sex)) {
      reference_sex <- levels(df[[sex_col]])[1]
    }
  }
  
  # UNCOMMENTED AND FIXED: Prepare Intellectual Disability column and baseline choice
  if (!is.null(intellectual_disability_col) && intellectual_disability_col %in% colnames(df)) {
    df[[intellectual_disability_col]] <- factor(df[[intellectual_disability_col]])
    if (is.null(reference_id)) {
      reference_id <- levels(df[[intellectual_disability_col]])[1]
    }
  }
  
  # ---- Determine Group-Specific Age Bounds (N >= 5 at limits per 6 unique cells) ----
  bounds_df <- df %>%
    group_by(!!sym(apoe_col), !!sym(sex_col)) %>%
    summarise(
      min_valid_age = if(n() >= 5) sort(!!sym(age_col))[5] else min(!!sym(age_col), na.rm=TRUE),
      max_valid_age = if(n() >= 5) sort(!!sym(age_col), decreasing = TRUE)[5] else max(!!sym(age_col), na.rm=TRUE),
      .groups = 'drop'
    )
  
  # ---- Fit GAM (Evaluating association as a function of both APOE and Sex) ----
  formula_str <- paste0(y_col, " ~ apoe_sex + s(", age_col, ", by = apoe_sex, k = ", k_val, ")",
                        " + s(", age_col, ", k = ", k_val, ")")
  
  if (!is.null(cohort) && cohort %in% colnames(df)) {
    formula_str <- paste0(formula_str, " + ", cohort)
  }
  
  if (!is.null(intellectual_disability_col) && intellectual_disability_col %in% colnames(df)) {
    formula_str <- paste0(formula_str, " + ", intellectual_disability_col)
  }
  
  formula_gam <- as.formula(formula_str)
  mod <- gam(formula_gam, method = "REML", data = df)
  
  # ---- Prediction Grid (Balanced at reference parameters) ----
  grid_list <- list(
    seq(min(df[[age_col]], na.rm = TRUE), max(df[[age_col]], na.rm = TRUE), by = 0.1),
    levels(df[[apoe_col]]),
    levels(df[[sex_col]])
  )
  grid_names <- c(age_col, apoe_col, sex_col)
  
  if (!is.null(cohort) && cohort %in% colnames(df)) {
    grid_list <- c(grid_list, list(reference_cohort))
    grid_names <- c(grid_names, cohort)
  }
  
  # ADDED: Include Intellectual Disability reference level in the prediction grid matrix
  if (!is.null(intellectual_disability_col) && intellectual_disability_col %in% colnames(df)) {
    grid_list <- c(grid_list, list(reference_id))
    grid_names <- c(grid_names, intellectual_disability_col)
  }
  
  mod_df <- expand.grid(grid_list, stringsAsFactors = FALSE)
  colnames(mod_df) <- grid_names
  mod_df$apoe_sex <- interaction(mod_df[[apoe_col]], mod_df[[sex_col]], sep = "_")
  
  if (!is.null(cohort) && cohort %in% colnames(df)) {
    mod_df[[cohort]] <- factor(mod_df[[cohort]], levels = levels(df[[cohort]]))
  }
  
  # ADDED: Enforce factor attributes on the grid column to match your model structure
  if (!is.null(intellectual_disability_col) && intellectual_disability_col %in% colnames(df)) {
    mod_df[[intellectual_disability_col]] <- factor(mod_df[[intellectual_disability_col]], levels = levels(df[[intellectual_disability_col]]))
  }
  
  # ---- Generate Model Predictions ----
  tmp <- predict(mod, newdata = mod_df, se.fit = TRUE)
  mod_df$forecast <- tmp$fit
  
  # ---- Apply Sample Density Constraints ----
  mod_df <- mod_df %>%
    left_join(bounds_df, by = c(apoe_col, sex_col)) %>%
    filter(!!sym(age_col) >= min_valid_age & !!sym(age_col) <= max_valid_age)
  
  # ---- Plot (6 Subgroups, Shared Colors, Conditional Linetypes, No Ribbons) ----
  linetype_vector <- c(
    "M"      = "solid", 
    "Male"   = "solid", 
    "male"   = "solid",
    "F"      = "dashed", 
    "Female" = "dashed", 
    "female" = "dashed"
  )
  
  plot_aes <- aes(
    x = !!sym(age_col),
    y = !!sym(y_col),
    colour = !!sym(apoe_col)
  )
  
  if (!is.null(cohort) && cohort %in% colnames(df)) {
    plot_aes$shape <- sym(cohort)
  }
  
  p <- ggplot(df, plot_aes) +
    geom_point(alpha = 0.2) + 
    geom_line(
      data = mod_df,
      aes(
        x = !!sym(age_col),
        y = forecast,
        group = apoe_sex,
        linetype = !!sym(sex_col) 
      ),
      linewidth = 1
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal"
    ) +
    scale_colour_manual(values = colour_values) +
    scale_linetype_manual(values = linetype_vector) + 
    scale_x_continuous(breaks = c(20, 30, 40, 50, 60, 70)) +
    xlab("Age")
  
  if (!is.null(cohort) && cohort %in% colnames(df)) {
    p <- p + scale_shape_manual(values = c(1, 4))
  }
  
  # ---- Coordinate Viewport Windows ----
  final_xlim <- if (is.numeric(XLIM) && length(XLIM) == 2) XLIM else c(min(df[[age_col]], na.rm = TRUE), max(df[[age_col]], na.rm = TRUE))
  final_ylim <- if (is.numeric(YLIM) && length(YLIM) == 2) YLIM else NULL
  
  p <- p + coord_cartesian(xlim = final_xlim, ylim = final_ylim)
  
  return(list(
    plot = p,
    predictions = mod_df,
    model_obj = mod
  ))
}


run_apoe_model_cohort <- function(
    df,
    apoe_col = "APOE_grouped",
    y_col = "WUSTLcentiloid",
    age_col = "clinical_AgefromBaseline",
    intellectual_disability_col = "Intellectual_Disability",
    cohort = NULL,                                         # Changed default to NULL to match your call safely
    sex_col = NULL,                                        
    fill_values,
    colour_values,
    sig_bar_colors = c("sig_2_3" = "red", "sig_2_4" = "blue", "sig_3_4" = "green"),
    markSig = FALSE,
    k_val = 4,
    reference_apoe = "APOE3",
    reference_cohort = NULL,
    reference_sex = NULL,                                  
    reference_id = NULL,                                  
    XLIM = c(18, 72),                                     
    YLIM = FALSE
) {
  
  # ---- Fix column types ----
  df[[apoe_col]]  <- factor(df[[apoe_col]], levels = c("APOE2", "APOE3", "APOE4"))
  other_levels <- setdiff(c("APOE2", "APOE3", "APOE4"), reference_apoe)
  df[[apoe_col]]  <- factor(df[[apoe_col]], levels = c(reference_apoe, other_levels))
  df[[age_col]]   <- as.numeric(df[[age_col]])
  
  if (!is.null(cohort) && cohort %in% colnames(df)) {
    df[[cohort]]    <- factor(df[[cohort]])
    if (is.null(reference_cohort)) {
      reference_cohort <- levels(df[[cohort]])[1]
    }
  }
  
  if (!is.null(sex_col) && sex_col %in% colnames(df)) {
    df[[sex_col]]    <- factor(df[[sex_col]])
    if (is.null(reference_sex)) {
      reference_sex <- levels(df[[sex_col]])[1]
    }
  }
  
  if (!is.null(intellectual_disability_col) && intellectual_disability_col %in% colnames(df)) {
    df[[intellectual_disability_col]] <- factor(df[[intellectual_disability_col]])
    if (is.null(reference_id)) {
      reference_id <- levels(df[[intellectual_disability_col]])[1]
    }
  }
  
  # ---- Determine Group-Specific Age Bounds (NA Safe) ----
  bounds_df <- df %>%
    filter(!is.na(!!sym(age_col))) %>%                     # Strip NAs out first
    group_by(!!sym(apoe_col)) %>%
    summarise(
      # Safe fallback: if fewer than 5 obs, grab the absolute min/max instead of [5]
      min_valid_age = if(n() >= 5) sort(!!sym(age_col))[5] else min(!!sym(age_col)),
      max_valid_age = if(n() >= 5) sort(!!sym(age_col), decreasing = TRUE)[5] else max(!!sym(age_col)),
      .groups = 'drop'
    )
  
  # ---- Fit GAM (cohort-adjusted) ----
  formula_str <- paste0(
    y_col, " ~ ", apoe_col, " + s(", age_col, ", k = ", k_val, ")",
    " + s(", age_col, ", by = ", apoe_col, ", k = ", k_val, ")"
  )
  
  if (!is.null(cohort) && cohort %in% colnames(df)) {
    formula_str <- paste0(formula_str, " + ", cohort)
  }
  if (!is.null(sex_col) && sex_col %in% colnames(df)) {
    formula_str <- paste0(formula_str, " + ", sex_col)
  }
  if (!is.null(intellectual_disability_col) && intellectual_disability_col %in% colnames(df)) {
    formula_str <- paste0(formula_str, " + ", intellectual_disability_col)
  }
  
  formula_gam <- as.formula(formula_str)
  mod <- gam(formula_gam, method = "REML", data = df)
  
  # ---- Prediction Grid Setup ----
  grid_list <- list(
    seq(min(df[[age_col]], na.rm = TRUE), max(df[[age_col]], na.rm = TRUE), by = 0.1),
    levels(df[[apoe_col]])
  )
  grid_names <- c(age_col, apoe_col)
  
  if (!is.null(cohort) && cohort %in% colnames(df)) {
    grid_list <- c(grid_list, list(reference_cohort))
    grid_names <- c(grid_names, cohort)
  }
  if (!is.null(sex_col) && sex_col %in% colnames(df)) {
    grid_list <- c(grid_list, list(reference_sex))
    grid_names <- c(grid_names, sex_col)
  }
  if (!is.null(intellectual_disability_col) && intellectual_disability_col %in% colnames(df)) {
    grid_list <- c(grid_list, list(reference_id))
    grid_names <- c(grid_names, intellectual_disability_col)
  }
  
  mod_df <- expand.grid(grid_list, stringsAsFactors = FALSE)
  colnames(mod_df) <- grid_names
  
  if (!is.null(cohort) && cohort %in% colnames(df)) {
    mod_df[[cohort]] <- factor(mod_df[[cohort]], levels = levels(df[[cohort]]))
  }
  if (!is.null(sex_col) && sex_col %in% colnames(df)) {
    mod_df[[sex_col]] <- factor(mod_df[[sex_col]], levels = levels(df[[sex_col]]))
  }
  if (!is.null(intellectual_disability_col) && intellectual_disability_col %in% colnames(df)) {
    mod_df[[intellectual_disability_col]] <- factor(mod_df[[intellectual_disability_col]], levels = levels(df[[intellectual_disability_col]]))
  }
  
  tmp <- predict(mod, newdata = mod_df, se.fit = TRUE)
  mod_df$forecast <- tmp$fit
  mod_df$se        <- tmp$se.fit
  mod_df$min      <- mod_df$forecast - 1.96 * mod_df$se
  mod_df$max      <- mod_df$forecast + 1.96 * mod_df$se
  
  # ---- Apply Constraints ----
  mod_df <- mod_df %>%
    left_join(bounds_df, by = apoe_col) %>%
    filter(!!sym(age_col) >= min_valid_age & !!sym(age_col) <= max_valid_age)
  
  if (is.numeric(YLIM) && length(YLIM) == 2) {
    mod_df$min <- pmax(mod_df$min, YLIM[1])
    mod_df$max <- pmin(mod_df$max, YLIM[2])
  }
  
  # ---- Pivot wider for differences ----
  df_wide <- mod_df %>%
    dplyr::select(all_of(c(age_col, apoe_col, "forecast", "se"))) %>%
    pivot_wider(names_from = all_of(apoe_col), values_from = c("forecast", "se"))
  
  DT <- as.data.table(df_wide)
  
  req_cols <- c("forecast_APOE2", "forecast_APOE3", "forecast_APOE4")
  if (all(req_cols %in% names(DT))) {
    DT[, `:=`(
      diff_2_3 = forecast_APOE2 - forecast_APOE3,
      diff_2_4 = forecast_APOE2 - forecast_APOE4,
      diff_3_4 = forecast_APOE3 - forecast_APOE4,
      se_2_3 = sqrt(se_APOE2^2 + se_APOE3^2),
      se_2_4 = sqrt(se_APOE2^2 + se_APOE4^2),
      se_3_4 = sqrt(se_APOE3^2 + se_APOE4^2)
    )]
    
    DT[, `:=`(
      diff_2_3_low = diff_2_3 - 1.96 * se_2_3, diff_2_3_hi  = diff_2_3 + 1.96 * se_2_3,
      diff_2_4_low = diff_2_4 - 1.96 * se_2_4, diff_2_4_hi  = diff_2_4 + 1.96 * se_2_4,
      diff_3_4_low = diff_3_4 - 1.96 * se_3_4, diff_3_4_hi  = diff_3_4 + 1.96 * se_3_4
    )]
    
    DT[, `:=`(
      sig_2_3 = diff_2_3_low > 0 | diff_2_3_hi < 0,
      sig_2_4 = diff_2_4_low > 0 | diff_2_4_hi < 0,
      sig_3_4 = diff_3_4_low > 0 | diff_3_4_hi < 0
    )]
    
    # Simple check to make sure your custom function 'get_sig_ranges' exists in your global environment
    if (exists("get_sig_ranges")) {
      significance_summary <- rbindlist(list(
        get_sig_ranges(DT, "sig_2_3", age_col),
        get_sig_ranges(DT, "sig_2_4", age_col),
        get_sig_ranges(DT, "sig_3_4", age_col)
      ), fill = TRUE)
    } else {
      significance_summary <- data.table()
    }
  } else {
    significance_summary <- data.table()
  }
  
  # ---- Plot ----
  plot_aes <- aes(x = !!sym(age_col), y = !!sym(y_col), colour = !!sym(apoe_col), group = !!sym(apoe_col))
  if (!is.null(cohort) && cohort %in% colnames(df)) { plot_aes$shape <- sym(cohort) }
  
  p <- ggplot(df, plot_aes) +
    geom_point(alpha = 0.3) +
    geom_ribbon(data = mod_df, aes(x = !!sym(age_col), y = forecast, ymin = min, ymax = max, fill = !!sym(apoe_col)), alpha = 0.4, inherit.aes = FALSE) +
    geom_line(data = mod_df, aes(x = !!sym(age_col), y = forecast, group = !!sym(apoe_col))) +
    theme_bw() + theme(legend.position = "bottom") +
    scale_fill_manual(values = fill_values) + scale_colour_manual(values = colour_values)  +
    scale_x_continuous(breaks = c(20, 30, 40, 50, 60, 70)) + xlab("Age")
  
  if (!is.null(cohort) && cohort %in% colnames(df)) { p <- p + scale_shape_manual(values = c(1, 4)) }
  
  if (nrow(significance_summary) > 0 && markSig == TRUE) {
    y_range_max <- if(is.numeric(YLIM)) YLIM[2] else max(df[[y_col]], na.rm = TRUE)
    y_range_min <- if(is.numeric(YLIM)) YLIM[1] else min(df[[y_col]], na.rm = TRUE)
    y_buffer <- (y_range_max - y_range_min) * 0.05 
    contrast_levels <- unique(significance_summary$contrast)
    y_positions <- setNames(seq(from = y_range_max + y_buffer, by = y_buffer, length.out = length(contrast_levels)), contrast_levels)
    significance_summary[, y := y_positions[contrast]]
    
    p <- p + geom_segment(data = significance_summary, aes(x = start_age, xend = end_age, y = y, yend = y, colour = contrast), inherit.aes = FALSE, linewidth = 2) +
      scale_colour_manual(values = c(colour_values, sig_2_3 = sig_bar_colors["sig_2_3"], sig_2_4 = sig_bar_colors["sig_2_4"], sig_3_4 = sig_bar_colors["sig_3_4"]), guide = guide_legend(override.aes = list(shape = NA)))
  }
  
  final_xlim <- if (is.numeric(XLIM) && length(XLIM) == 2) XLIM else c(min(df[[age_col]], na.rm = TRUE), max(df[[age_col]], na.rm = TRUE))
  final_ylim <- if (is.numeric(YLIM) && length(YLIM) == 2) YLIM else NULL
  p <- p + coord_cartesian(xlim = final_xlim, ylim = final_ylim)
  
  return(list(plot = p, predictions = mod_df, significance_summary = significance_summary, model_obj = mod))
}

