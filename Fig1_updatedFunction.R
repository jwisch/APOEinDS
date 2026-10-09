run_apoe_model_cohort <- function(
    df,
    apoe_col = "APOE_grouped",
    y_col = "WUSTLcentiloid",
    age_col = "clinical_AgefromBaseline",
    intellectual_disability_col = "Intellectual_Disability",
    cohort = NULL,
    sex_col = NULL,
    fill_values,
    colour_values,
    family_choice = gaussian(),
    sig_bar_colors = c(
      "sig_2_3" = "#BABBB1",
      "sig_2_4" = "#D67B44",
      "sig_3_4" = "#514289"
    ),
    markSig = FALSE,
    k_val = 4,
    y_shift = 0,
    reference_apoe = "APOE3",
    reference_cohort = NULL,
    reference_sex = NULL,
    reference_id = NULL,
    XLIM = c(18, 72),
    YLIM = FALSE
) {
  
  # ============================================================
  # 1. Fix column types
  # ============================================================
  
  df[[apoe_col]] <- factor(
    df[[apoe_col]],
    levels = c("APOE2", "APOE3", "APOE4")
  )
  
  other_levels <- setdiff(
    c("APOE2", "APOE3", "APOE4"),
    reference_apoe
  )
  
  df[[apoe_col]] <- factor(
    df[[apoe_col]],
    levels = c(reference_apoe, other_levels)
  )
  
  df[[age_col]] <- as.numeric(df[[age_col]])
  
  if (!is.null(cohort) && cohort %in% colnames(df)) {
    df[[cohort]] <- factor(df[[cohort]])
    
    if (is.null(reference_cohort)) {
      reference_cohort <- levels(df[[cohort]])[1]
    }
  }
  
  if (!is.null(sex_col) && sex_col %in% colnames(df)) {
    df[[sex_col]] <- factor(df[[sex_col]])
    
    if (is.null(reference_sex)) {
      reference_sex <- levels(df[[sex_col]])[1]
    }
  }
  
  if (!is.null(intellectual_disability_col) &&
      intellectual_disability_col %in% colnames(df)) {
    
    df[[intellectual_disability_col]] <-
      factor(df[[intellectual_disability_col]])
    
    if (is.null(reference_id)) {
      reference_id <- levels(df[[intellectual_disability_col]])[1]
    }
  }
  
  
  # ============================================================
  # 2. Determine group-specific age bounds
  # ============================================================
  
  bounds_df <- df %>%
    dplyr::filter(!is.na(!!rlang::sym(age_col))) %>%
    dplyr::group_by(!!rlang::sym(apoe_col)) %>%
    dplyr::summarise(
      min_valid_age =
        if (dplyr::n() >= 5) {
          sort(!!rlang::sym(age_col))[5]
        } else {
          min(!!rlang::sym(age_col))
        },
      
      max_valid_age =
        if (dplyr::n() >= 5) {
          sort(
            !!rlang::sym(age_col),
            decreasing = TRUE
          )[5]
        } else {
          max(!!rlang::sym(age_col))
        },
      
      .groups = "drop"
    )
  
  
  # ============================================================
  # 3. Construct GAM formula
  # ============================================================
  
  formula_str <- paste0(
    y_col,
    " ~ ",
    apoe_col,
    " + s(",
    age_col,
    ", k = ",
    k_val,
    ")",
    " + s(",
    age_col,
    ", by = ",
    apoe_col,
    ", k = ",
    k_val,
    ")"
  )
  
  if (!is.null(cohort) && cohort %in% colnames(df)) {
    formula_str <- paste0(
      formula_str,
      " + ",
      cohort
    )
  }
  
  if (!is.null(sex_col) && sex_col %in% colnames(df)) {
    formula_str <- paste0(
      formula_str,
      " + ",
      sex_col
    )
  }
  
  if (!is.null(intellectual_disability_col) &&
      intellectual_disability_col %in% colnames(df)) {
    
    formula_str <- paste0(
      formula_str,
      " + ",
      intellectual_disability_col
    )
  }
  
  formula_gam <- as.formula(formula_str)
  
  
  # ============================================================
  # 4. Fit GAM
  # ============================================================
  
  mod <- mgcv::gam(
    formula_gam,
    method = "REML",
    data = df,
    family = family_choice
  )
  
  
  # ============================================================
  # 5. Capture gam.check output
  # ============================================================
  
  pdf(NULL)
  
  check_output <- capture.output(
    mgcv::gam.check(mod)
  )
  
  dev.off()
  
  
  # ============================================================
  # 6. Construct prediction grid
  # ============================================================
  
  age_grid <- seq(
    min(df[[age_col]], na.rm = TRUE),
    max(df[[age_col]], na.rm = TRUE),
    by = 0.1
  )
  
  grid_list <- list(
    age_grid,
    levels(df[[apoe_col]])
  )
  
  grid_names <- c(
    age_col,
    apoe_col
  )
  
  if (!is.null(cohort) &&
      cohort %in% colnames(df)) {
    
    grid_list <- c(
      grid_list,
      list(reference_cohort)
    )
    
    grid_names <- c(
      grid_names,
      cohort
    )
  }
  
  if (!is.null(sex_col) &&
      sex_col %in% colnames(df)) {
    
    grid_list <- c(
      grid_list,
      list(reference_sex)
    )
    
    grid_names <- c(
      grid_names,
      sex_col
    )
  }
  
  if (!is.null(intellectual_disability_col) &&
      intellectual_disability_col %in% colnames(df)) {
    
    grid_list <- c(
      grid_list,
      list(reference_id)
    )
    
    grid_names <- c(
      grid_names,
      intellectual_disability_col
    )
  }
  
  mod_df <- expand.grid(
    grid_list,
    stringsAsFactors = FALSE
  )
  
  colnames(mod_df) <- grid_names
  
  
  # Restore factor levels for prediction
  
  if (!is.null(cohort) &&
      cohort %in% colnames(df)) {
    
    mod_df[[cohort]] <- factor(
      mod_df[[cohort]],
      levels = levels(df[[cohort]])
    )
  }
  
  if (!is.null(sex_col) &&
      sex_col %in% colnames(df)) {
    
    mod_df[[sex_col]] <- factor(
      mod_df[[sex_col]],
      levels = levels(df[[sex_col]])
    )
  }
  
  if (!is.null(intellectual_disability_col) &&
      intellectual_disability_col %in% colnames(df)) {
    
    mod_df[[intellectual_disability_col]] <- factor(
      mod_df[[intellectual_disability_col]],
      levels = levels(
        df[[intellectual_disability_col]]
      )
    )
  }
  
  
  # ------------------------------------------------------------
  # 7. Link-scale prediction
  # ------------------------------------------------------------
  
  pred_link <- predict(
    mod,
    newdata = mod_df,
    type = "link",
    se.fit = TRUE
  )
  
  z <- qnorm(0.975)
  
  mod_df$eta <- as.numeric(pred_link$fit)
  mod_df$se_eta <- as.numeric(pred_link$se.fit)
  
  # Response-scale transformation
  linkinv <- mod$family$linkinv
  mu_eta <- mod$family$mu.eta(mod_df$eta)
  
  mod_df$forecast <- linkinv(mod_df$eta)
  
  # Delta-method SE on response scale
  mod_df$se <- abs(mu_eta) * mod_df$se_eta
  
  # Response-scale confidence interval
  mod_df$min <- linkinv(
    mod_df$eta - z * mod_df$se_eta
  )
  
  mod_df$max <- linkinv(
    mod_df$eta + z * mod_df$se_eta
  )
  
  
  # ------------------------------------------------------------
  # 8. Age constraints
  # ------------------------------------------------------------
  
  mod_df <- mod_df %>%
    dplyr::left_join(
      bounds_df,
      by = apoe_col
    ) %>%
    dplyr::filter(
      !!rlang::sym(age_col) >= min_valid_age &
        !!rlang::sym(age_col) <= max_valid_age
    )
  
  
  if (nrow(mod_df) == 0) {
    stop("No prediction points remain after applying the APOE-specific age bounds.")
  }
  
  
  # ------------------------------------------------------------
  # 9. Recalculate prediction covariance AFTER filtering
  # ------------------------------------------------------------
  
  Xp <- predict(
    mod,
    newdata = mod_df,
    type = "lpmatrix"
  )
  
  Vp <- vcov(mod)
  
  pred_cov_link <- Xp %*% Vp %*% t(Xp)
  
  
  # ------------------------------------------------------------
  # 10. Plot CI clipping
  # ------------------------------------------------------------
  
  if (is.numeric(YLIM) && length(YLIM) == 2) {
    
    mod_df$min_plot <- pmax(
      mod_df$min,
      YLIM[1]
    )
    
    mod_df$max_plot <- pmin(
      mod_df$max,
      YLIM[2]
    )
    
  } else {
    
    mod_df$min_plot <- mod_df$min
    mod_df$max_plot <- mod_df$max
    
  }
  
  
  # ------------------------------------------------------------
  # 11. APOE contrasts
  # ------------------------------------------------------------
  
  if (
    all(
      c("APOE2", "APOE3", "APOE4") %in%
      levels(df[[apoe_col]])
    )
  ) {
    
    # Use the already-filtered prediction data
    pred_long <- mod_df %>%
      dplyr::select(
        dplyr::all_of(
          c(
            age_col,
            apoe_col,
            "eta",
            "forecast"
          )
        )
      )
    
    pred_2 <- pred_long %>%
      dplyr::filter(
        !!rlang::sym(apoe_col) == "APOE2"
      ) %>%
      dplyr::rename(
        forecast_APOE2 = forecast,
        eta_APOE2 = eta
      )
    
    pred_3 <- pred_long %>%
      dplyr::filter(
        !!rlang::sym(apoe_col) == "APOE3"
      ) %>%
      dplyr::rename(
        forecast_APOE3 = forecast,
        eta_APOE3 = eta
      )
    
    pred_4 <- pred_long %>%
      dplyr::filter(
        !!rlang::sym(apoe_col) == "APOE4"
      ) %>%
      dplyr::rename(
        forecast_APOE4 = forecast,
        eta_APOE4 = eta
      )
    
    
    # Only compare ages where all three APOE groups
    # have a valid prediction
    df_wide <- pred_2 %>%
      dplyr::inner_join(
        pred_3,
        by = age_col
      ) %>%
      dplyr::inner_join(
        pred_4,
        by = age_col
      )
    
    
    contrast_results <- lapply(
      seq_len(nrow(df_wide)),
      function(i) {
        
        age_i <- df_wide[[age_col]][i]
        
        r2 <- which(
          mod_df[[apoe_col]] == "APOE2" &
            mod_df[[age_col]] == age_i
        )
        
        r3 <- which(
          mod_df[[apoe_col]] == "APOE3" &
            mod_df[[age_col]] == age_i
        )
        
        r4 <- which(
          mod_df[[apoe_col]] == "APOE4" &
            mod_df[[age_col]] == age_i
        )
        
        if (
          length(r2) != 1 ||
          length(r3) != 1 ||
          length(r4) != 1
        ) {
          return(NULL)
        }
        
        r <- c(r2, r3, r4)
        
        
        # ------------------------------------------------------
        # Predictions at this age
        # ------------------------------------------------------
        
        eta_i <- mod_df$eta[r]
        mu_i <- mod_df$forecast[r]
        
        
        # ------------------------------------------------------
        # Covariance of the three predictions
        # ------------------------------------------------------
        
        Sigma_eta <- pred_cov_link[
          r,
          r,
          drop = FALSE
        ]
        
        
        # ------------------------------------------------------
        # Transform covariance from link to response scale
        # ------------------------------------------------------
        
        deriv_i <- mod$family$mu.eta(eta_i)
        
        J <- diag(
          as.numeric(deriv_i)
        )
        
        Sigma_mu <- J %*%
          Sigma_eta %*%
          J
        
        
        # ------------------------------------------------------
        # Contrasts
        # ------------------------------------------------------
        
        c23 <- c(1, -1, 0)
        c24 <- c(1, 0, -1)
        c34 <- c(0, 1, -1)
        
        
        diff_23 <- sum(c23 * mu_i)
        diff_24 <- sum(c24 * mu_i)
        diff_34 <- sum(c34 * mu_i)
        
        
        se_23 <- sqrt(
          max(
            0,
            as.numeric(
              t(c23) %*%
                Sigma_mu %*%
                c23
            )
          )
        )
        
        se_24 <- sqrt(
          max(
            0,
            as.numeric(
              t(c24) %*%
                Sigma_mu %*%
                c24
            )
          )
        )
        
        se_34 <- sqrt(
          max(
            0,
            as.numeric(
              t(c34) %*%
                Sigma_mu %*%
                c34
            )
          )
        )
        
        
        data.frame(
          age = age_i,
          
          diff_2_3 = diff_23,
          diff_2_4 = diff_24,
          diff_3_4 = diff_34,
          
          se_2_3 = se_23,
          se_2_4 = se_24,
          se_3_4 = se_34
        )
      }
    )
    
    
    contrast_results <- dplyr::bind_rows(
      contrast_results
    )
    
    names(contrast_results)[
      names(contrast_results) == "age"
    ] <- age_col
    
    
    # ----------------------------------------------------------
    # Confidence intervals for contrasts
    # ----------------------------------------------------------
    
    DT <- data.table::as.data.table(
      contrast_results
    )
    
    DT[, `:=`(
      
      diff_2_3_low =
        diff_2_3 - z * se_2_3,
      
      diff_2_3_hi =
        diff_2_3 + z * se_2_3,
      
      diff_2_4_low =
        diff_2_4 - z * se_2_4,
      
      diff_2_4_hi =
        diff_2_4 + z * se_2_4,
      
      diff_3_4_low =
        diff_3_4 - z * se_3_4,
      
      diff_3_4_hi =
        diff_3_4 + z * se_3_4
      
    )]
    
    
    # ----------------------------------------------------------
    # Significance
    # ----------------------------------------------------------
    
    DT[, `:=`(
      
      sig_2_3 =
        diff_2_3_low > 0 |
        diff_2_3_hi < 0,
      
      sig_2_4 =
        diff_2_4_low > 0 |
        diff_2_4_hi < 0,
      
      sig_3_4 =
        diff_3_4_low > 0 |
        diff_3_4_hi < 0
      
    )]
    
    
    # ----------------------------------------------------------
    # Significant age ranges
    # ----------------------------------------------------------
    
    if (exists("get_sig_ranges")) {
      
      significance_summary <- data.table::rbindlist(
        
        list(
          
          get_sig_ranges(
            DT,
            "sig_2_3",
            age_col
          ),
          
          get_sig_ranges(
            DT,
            "sig_2_4",
            age_col
          ),
          
          get_sig_ranges(
            DT,
            "sig_3_4",
            age_col
          )
          
        ),
        
        fill = TRUE
      )
      
    } else {
      
      significance_summary <-
        data.table::data.table()
      
    }
    
  } else {
    
    DT <- data.table::data.table()
    
    significance_summary <-
      data.table::data.table()
    
  }
  
  
  # 13. plot
  
  plot_aes <- ggplot2::aes(
    x = !!rlang::sym(age_col),
    y = !!rlang::sym(y_col) + y_shift,
    colour = !!rlang::sym(apoe_col),
    group = !!rlang::sym(apoe_col)
  )
  
  if (!is.null(cohort) && cohort %in% colnames(df)) {
    plot_aes$shape <- rlang::sym(cohort)
  }
  
  p <- ggplot2::ggplot(df, plot_aes) +
    
    # Raw data
    ggplot2::geom_point(alpha = 0.3) +
    
    # Confidence intervals
    ggplot2::geom_ribbon(
      data = mod_df,
      ggplot2::aes(
        x = !!rlang::sym(age_col),
        y = forecast + y_shift ,
        ymin = min_plot + y_shift,
        ymax = max_plot + y_shift,
        fill = !!rlang::sym(apoe_col)
      ),
      alpha = 0.4,
      inherit.aes = FALSE
    ) +
    
    # Model fits
    ggplot2::geom_line(
      data = mod_df,
      ggplot2::aes(
        x = !!rlang::sym(age_col),
        y = forecast + y_shift,
        colour = !!rlang::sym(apoe_col),
        group = !!rlang::sym(apoe_col)
      ),
      inherit.aes = FALSE
    ) +
    
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom") +
    
    # Fill colors for confidence intervals
    ggplot2::scale_fill_manual(
      values = fill_values
    ) +
    
    # APOE colors
    ggplot2::scale_colour_manual(
      values = c(
        "APOE2" = colour_values[1],
        "APOE3" = colour_values[2],
        "APOE4" = colour_values[3]
      )
    ) +
    
    ggplot2::scale_x_continuous(
      breaks = c(20, 30, 40, 50, 60, 70)
    ) +
    
    ggplot2::xlab("Age")
  
  
  # Shape scale for cohort
  if (!is.null(cohort) && cohort %in% colnames(df)) {
    p <- p +
      ggplot2::scale_shape_manual(
        values = c(1, 4)
      )
  }
  
  
  # ------------------------------------------------------------
  # Significance bars
  # ------------------------------------------------------------
  
  if (nrow(significance_summary) > 0 && markSig == TRUE) {
    
    y_range_max <- if (is.numeric(YLIM)) {
      YLIM[2]
    } else {
      max(df[[y_col]], na.rm = TRUE)
    }
    
    y_range_min <- if (is.numeric(YLIM)) {
      YLIM[1]
    } else {
      min(df[[y_col]], na.rm = TRUE)
    }
    
    y_buffer <- (y_range_max - y_range_min) * 0.05
    
    contrast_levels <- unique(significance_summary$contrast)
    
    y_positions <- setNames(
      seq(
        from = y_range_max + y_buffer,
        by = y_buffer,
        length.out = length(contrast_levels)
      ),
      contrast_levels
    )
    
    significance_summary[, y := y_positions[contrast]]
    
    
    # Add significance bars
    p <- p +
      ggplot2::geom_segment(
        data = significance_summary,
        ggplot2::aes(
          x = start_age,
          xend = end_age,
          y = y + y_shift,
          yend = y + y_shift,
          colour = contrast,
          group = contrast
        ),
        inherit.aes = FALSE,
        linewidth = 2
      )
    
    
    # IMPORTANT:
    # Rebuild ONE combined colour scale containing both
    # APOE colors and significance-bar colors.
    combined_colour_values <- c(
      "APOE2" = colour_values[1],
      "APOE3" = colour_values[2],
      "APOE4" = colour_values[3],
      "sig_2_3" = sig_bar_colors["sig_2_3"],
      "sig_2_4" = sig_bar_colors["sig_2_4"],
      "sig_3_4" = sig_bar_colors["sig_3_4"]
    )
    
    p <- p +
      ggplot2::scale_colour_manual(
        values = combined_colour_values,
        guide = ggplot2::guide_legend(
          override.aes = list(shape = NA)
        )
      )
  }
  
  
  # ------------------------------------------------------------
  # Axis limits
  # ------------------------------------------------------------
  
  final_xlim <- if (
    is.numeric(XLIM) && length(XLIM) == 2
  ) {
    XLIM
  } else {
    c(
      min(df[[age_col]], na.rm = TRUE),
      max(df[[age_col]], na.rm = TRUE)
    )
  }
  
  final_ylim <- if (
    is.numeric(YLIM) && length(YLIM) == 2
  ) {
    YLIM
  } else {
    NULL
  }
  
  p <- p +
    ggplot2::coord_cartesian(
      xlim = final_xlim,
      ylim = final_ylim
    )
  
  
  return(
    list(
      plot = p,
      predictions = mod_df,
      significance_summary = significance_summary,
      model_obj = mod,
      gam_check = check_output
    )
  )
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
    sig_bar_colors = c(
      "sig_APOE2" = "purple",
      "sig_APOE3" = "orange",
      "sig_APOE4" = "darkred"
    ),
    markSig = FALSE,
    k_val = 4,
    family_choice = gaussian(),
    reference_cohort = NULL,
    XLIM = c(18, 72),
    YLIM = FALSE
) {
  
  # ============================================================
  # 1. Packages
  # ============================================================
  
  library(mgcv)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(ggplot2)
  
  
  # ============================================================
  # 2. Fix column types
  # ============================================================
  
  df[[apoe_col]] <- factor(
    df[[apoe_col]],
    levels = c("APOE2", "APOE3", "APOE4")
  )
  
  df[[sex_col]] <- factor(df[[sex_col]])
  
  df[[age_col]] <- as.numeric(df[[age_col]])
  
  
  # Remove observations with missing variables needed for the model
  required_cols <- c(
    y_col,
    age_col,
    apoe_col,
    sex_col
  )
  
  if (!is.null(cohort) && cohort %in% colnames(df)) {
    required_cols <- c(required_cols, cohort)
  }
  
  df <- df[
    complete.cases(df[, required_cols, drop = FALSE]),
    ,
    drop = FALSE
  ]
  
  
  # ============================================================
  # 3. Validate sex levels
  # ============================================================
  
  sex_levels <- levels(df[[sex_col]])
  
  if (length(sex_levels) != 2) {
    stop(
      "The function requires exactly two sex levels. ",
      "Found: ",
      paste(sex_levels, collapse = ", ")
    )
  }
  
  f_level <- sex_levels[1]
  m_level <- sex_levels[2]
  
  
  # ============================================================
  # 4. Create APOE × sex interaction
  # ============================================================
  
  df$apoe_sex <- interaction(
    df[[apoe_col]],
    df[[sex_col]],
    sep = "_",
    drop = TRUE
  )
  
  
  # ============================================================
  # 5. Cohort
  # ============================================================
  
  if (!is.null(cohort) && cohort %in% colnames(df)) {
    
    df[[cohort]] <- factor(df[[cohort]])
    
    if (is.null(reference_cohort)) {
      reference_cohort <- levels(df[[cohort]])[1]
    }
    
    if (!reference_cohort %in% levels(df[[cohort]])) {
      stop(
        "reference_cohort = '",
        reference_cohort,
        "' is not a valid level of ",
        cohort,
        ". Valid levels are: ",
        paste(levels(df[[cohort]]), collapse = ", ")
      )
    }
  }
  
  
  # ============================================================
  # 6. APOE/sex-specific age bounds
  # ============================================================
  
  bounds_df <- df %>%
    dplyr::group_by(
      !!rlang::sym(apoe_col),
      !!rlang::sym(sex_col)
    ) %>%
    dplyr::summarise(
      
      min_valid_age =
        if (dplyr::n() >= 5) {
          sort(!!rlang::sym(age_col))[5]
        } else {
          min(
            !!rlang::sym(age_col),
            na.rm = TRUE
          )
        },
      
      max_valid_age =
        if (dplyr::n() >= 5) {
          sort(
            !!rlang::sym(age_col),
            decreasing = TRUE
          )[5]
        } else {
          max(
            !!rlang::sym(age_col),
            na.rm = TRUE
          )
        },
      
      .groups = "drop"
    )
  
  
  # ============================================================
  # 7. Fit GAM
  # ============================================================
  
  formula_str <- paste0(
    y_col,
    " ~ apoe_sex",
    " + s(",
    age_col,
    ", by = apoe_sex, k = ",
    k_val,
    ")",
    " + s(",
    age_col,
    ", k = ",
    k_val,
    ")"
  )
  
  if (!is.null(cohort) && cohort %in% colnames(df)) {
    
    formula_str <- paste0(
      formula_str,
      " + ",
      cohort
    )
  }
  
  formula_gam <- as.formula(formula_str)
  
  mod <- mgcv::gam(
    formula_gam,
    method = "REML",
    data = df,
    family = family_choice
  )
  
  
  # ============================================================
  # 8. gam.check
  # ============================================================
  
  pdf(NULL)
  
  check_output <- capture.output(
    mgcv::gam.check(mod)
  )
  
  dev.off()
  
  
  # ============================================================
  # 9. Prediction grid
  # ============================================================
  
  age_grid <- seq(
    min(df[[age_col]], na.rm = TRUE),
    max(df[[age_col]], na.rm = TRUE),
    by = 0.1
  )
  
  grid_list <- list(
    age_grid,
    levels(df[[apoe_col]]),
    levels(df[[sex_col]])
  )
  
  grid_names <- c(
    age_col,
    apoe_col,
    sex_col
  )
  
  if (!is.null(cohort) && cohort %in% colnames(df)) {
    
    grid_list <- c(
      grid_list,
      list(reference_cohort)
    )
    
    grid_names <- c(
      grid_names,
      cohort
    )
  }
  
  mod_df <- expand.grid(
    grid_list,
    stringsAsFactors = FALSE
  )
  
  colnames(mod_df) <- grid_names
  
  
  # Restore factor types
  mod_df[[apoe_col]] <- factor(
    mod_df[[apoe_col]],
    levels = levels(df[[apoe_col]])
  )
  
  mod_df[[sex_col]] <- factor(
    mod_df[[sex_col]],
    levels = levels(df[[sex_col]])
  )
  
  if (!is.null(cohort) && cohort %in% colnames(df)) {
    
    mod_df[[cohort]] <- factor(
      mod_df[[cohort]],
      levels = levels(df[[cohort]])
    )
  }
  
  
  # Create same interaction variable used by model
  mod_df$apoe_sex <- interaction(
    mod_df[[apoe_col]],
    mod_df[[sex_col]],
    sep = "_",
    drop = TRUE
  )
  
  
  # ============================================================
  # 10. Apply APOE × sex age constraints BEFORE prediction
  # ============================================================
  
  mod_df <- mod_df %>%
    dplyr::left_join(
      bounds_df,
      by = c(apoe_col, sex_col)
    ) %>%
    dplyr::filter(
      !!rlang::sym(age_col) >= min_valid_age,
      !!rlang::sym(age_col) <= max_valid_age
    )
  
  
  if (nrow(mod_df) == 0) {
    stop(
      "No prediction points remain after applying ",
      "the APOE × sex-specific age bounds."
    )
  }
  
  
  # ============================================================
  # 11. Link-scale prediction
  # ============================================================
  
  pred_link <- predict(
    mod,
    newdata = mod_df,
    type = "link",
    se.fit = TRUE
  )
  
  z <- qnorm(0.975)
  
  mod_df$eta <- as.numeric(
    pred_link$fit
  )
  
  mod_df$se_eta <- as.numeric(
    pred_link$se.fit
  )
  
  
  # ============================================================
  # 12. Transform prediction to response scale
  # ============================================================
  
  linkinv <- mod$family$linkinv
  
  mu_eta <- mod$family$mu.eta(
    mod_df$eta
  )
  
  mod_df$forecast <- linkinv(
    mod_df$eta
  )
  
  # Delta-method SE on response scale
  mod_df$se <- abs(mu_eta) *
    mod_df$se_eta
  
  
  # Response-scale confidence intervals
  #
  # IMPORTANT:
  # Transform the endpoints rather than doing
  # forecast +/- z*SE when the link is nonlinear.
  
  mod_df$min <- linkinv(
    mod_df$eta - z * mod_df$se_eta
  )
  
  mod_df$max <- linkinv(
    mod_df$eta + z * mod_df$se_eta
  )
  
  
  # ============================================================
  # 13. Prediction covariance matrix
  # ============================================================
  
  Xp <- predict(
    mod,
    newdata = mod_df,
    type = "lpmatrix"
  )
  
  Vp <- vcov(mod)
  
  pred_cov_link <- Xp %*%
    Vp %*%
    t(Xp)
  
  
  # ============================================================
  # 14. Plot CI clipping
  # ============================================================
  
  if (
    is.numeric(YLIM) &&
    length(YLIM) == 2
  ) {
    
    mod_df$min_plot <- pmax(
      mod_df$min,
      YLIM[1]
    )
    
    mod_df$max_plot <- pmin(
      mod_df$max,
      YLIM[2]
    )
    
  } else {
    
    mod_df$min_plot <- mod_df$min
    mod_df$max_plot <- mod_df$max
  }
  
  
  # ============================================================
  # 15. Sex contrasts within each APOE group
  # ============================================================
  
  significance_list <- list()
  
  
  for (apoe_val in levels(df[[apoe_col]])) {
    
    # Rows corresponding to this APOE group
    # in the prediction grid
    apoe_rows <- which(
      mod_df[[apoe_col]] == apoe_val
    )
    
    if (length(apoe_rows) == 0) {
      next
    }
    
    
    # Ages for which both sexes have predictions
    female_rows <- apoe_rows[
      mod_df[[sex_col]][apoe_rows] == f_level
    ]
    
    male_rows <- apoe_rows[
      mod_df[[sex_col]][apoe_rows] == m_level
    ]
    
    
    common_ages <- intersect(
      mod_df[[age_col]][female_rows],
      mod_df[[age_col]][male_rows]
    )
    
    
    if (length(common_ages) == 0) {
      next
    }
    
    
    contrast_results <- lapply(
      common_ages,
      function(age_i) {
        
        r_f <- female_rows[
          mod_df[[age_col]][female_rows] == age_i
        ]
        
        r_m <- male_rows[
          mod_df[[age_col]][male_rows] == age_i
        ]
        
        
        if (
          length(r_f) != 1 ||
          length(r_m) != 1
        ) {
          return(NULL)
        }
        
        
        r <- c(
          r_f,
          r_m
        )
        
        
        # ------------------------------------------------------
        # Link-scale predictions
        # ------------------------------------------------------
        
        eta_i <- mod_df$eta[r]
        
        
        # ------------------------------------------------------
        # Response-scale predictions
        # ------------------------------------------------------
        
        mu_i <- mod_df$forecast[r]
        
        
        # ------------------------------------------------------
        # Covariance of male/female predictions
        # ------------------------------------------------------
        
        Sigma_eta <- pred_cov_link[
          r,
          r,
          drop = FALSE
        ]
        
        
        # ------------------------------------------------------
        # Transform covariance to response scale
        # ------------------------------------------------------
        
        deriv_i <- mod$family$mu.eta(
          eta_i
        )
        
        J <- diag(
          as.numeric(deriv_i)
        )
        
        Sigma_mu <- J %*%
          Sigma_eta %*%
          J
        
        
        # ------------------------------------------------------
        # Male - Female contrast
        # ------------------------------------------------------
        
        contrast <- c(
          -1,
          1
        )
        
        diff_M_F <- sum(
          contrast * mu_i
        )
        
        
        se_M_F <- sqrt(
          max(
            0,
            as.numeric(
              t(contrast) %*%
                Sigma_mu %*%
                contrast
            )
          )
        )
        
        
        data.frame(
          
          age = age_i,
          
          diff_M_F = diff_M_F,
          
          se_M_F = se_M_F
          
        )
      }
    )
    
    
    contrast_results <- dplyr::bind_rows(
      contrast_results
    )
    
    
    if (nrow(contrast_results) == 0) {
      next
    }
    
    
    names(contrast_results)[
      names(contrast_results) == "age"
    ] <- age_col
    
    
    DT_tmp <- data.table::as.data.table(
      contrast_results
    )
    
    
    # ----------------------------------------------------------
    # Bonferroni correction for 3 planned APOE-specific tests
    # ----------------------------------------------------------
    
    n_comparisons <- 3
    
    alpha_adj <- 0.05 /
      n_comparisons
    
    z_crit <- qnorm(
      1 - alpha_adj / 2
    )
    
    
    DT_tmp[, `:=`(
      
      diff_low =
        diff_M_F -
        z_crit * se_M_F,
      
      diff_hi =
        diff_M_F +
        z_crit * se_M_F
      
    )]
    
    
    DT_tmp[, sig_diff := (
      diff_low > 0 |
        diff_hi < 0
    )]
    
    
    # ----------------------------------------------------------
    # Extract continuous significant age ranges
    # ----------------------------------------------------------
    
    ranges <- get_sig_ranges(
      DT_tmp,
      "sig_diff",
      age_col
    )
    
    
    if (
      !is.null(ranges) &&
      nrow(ranges) > 0
    ) {
      
      ranges[, contrast :=
               paste0(
                 "sig_",
                 apoe_val
               )]
      
      ranges[, (apoe_col) :=
               apoe_val]
      
      significance_list[[apoe_val]] <-
        ranges
    }
  }
  
  
  # ============================================================
  # 16. Combine significance ranges
  # ============================================================
  
  if (length(significance_list) > 0) {
    
    significance_summary <-
      data.table::rbindlist(
        significance_list,
        fill = TRUE
      )
    
  } else {
    
    significance_summary <-
      data.table::data.table()
  }
  
  
  # ============================================================
  # 17. Plot
  # ============================================================
  
  plot_aes <- ggplot2::aes(
    x = !!rlang::sym(age_col),
    y = !!rlang::sym(y_col),
    colour = !!rlang::sym(sex_col),
    group = !!rlang::sym(sex_col)
  )
  
  if (
    !is.null(cohort) &&
    cohort %in% colnames(df)
  ) {
    
    plot_aes$shape <-
      rlang::sym(cohort)
  }
  
  
  p <- ggplot2::ggplot(
    df,
    plot_aes
  ) +
    
    ggplot2::geom_point(
      alpha = 0.3
    ) +
    
    ggplot2::geom_ribbon(
      data = mod_df,
      
      ggplot2::aes(
        x = !!rlang::sym(age_col),
        y = forecast,
        ymin = min_plot,
        ymax = max_plot,
        fill = !!rlang::sym(sex_col),
        group = !!rlang::sym(sex_col)
      ),
      
      alpha = 0.25,
      
      inherit.aes = FALSE
    ) +
    
    ggplot2::geom_line(
      data = mod_df,
      
      ggplot2::aes(
        x = !!rlang::sym(age_col),
        y = forecast,
        colour = !!rlang::sym(sex_col),
        group = !!rlang::sym(sex_col)
      ),
      
      linewidth = 1,
      
      inherit.aes = FALSE
    ) +
    
    ggplot2::facet_wrap(
      ggplot2::vars(
        !!rlang::sym(apoe_col)
      )
    ) +
    
    ggplot2::theme_bw() +
    
    ggplot2::theme(
      legend.position = "bottom",
      strip.background = ggplot2::element_blank()
    ) +
    
    ggplot2::scale_fill_manual(
      values = fill_values
    ) +
    
    ggplot2::scale_colour_manual(
      values = colour_values
    ) +
    
    ggplot2::scale_x_continuous(
      breaks = c(30, 50, 70)
    )
  
  
  if (
    !is.null(cohort) &&
    cohort %in% colnames(df)
  ) {
    
    p <- p +
      ggplot2::scale_shape_manual(
        values = c(1, 4)
      )
  }
  
  
  # ============================================================
  # 18. Significance bars
  # ============================================================
  
  if (
    nrow(significance_summary) > 0 &&
    markSig == TRUE
  ) {
    
    y_range_max <- if (
      is.numeric(YLIM)
    ) {
      YLIM[2]
    } else {
      max(
        df[[y_col]],
        na.rm = TRUE
      )
    }
    
    
    y_range_min <- if (
      is.numeric(YLIM)
    ) {
      YLIM[1]
    } else {
      min(
        df[[y_col]],
        na.rm = TRUE
      )
    }
    
    
    y_buffer <- (
      y_range_max -
        y_range_min
    ) * 0.05
    
    
    # One bar per APOE panel
    significance_summary[
      ,
      y := y_range_max + y_buffer
    ]
    
    
    p <- p +
      
      ggplot2::geom_segment(
        data = significance_summary,
        
        ggplot2::aes(
          x = start_age,
          xend = end_age,
          y = y,
          yend = y,
          colour = contrast,
          group = contrast
        ),
        
        inherit.aes = FALSE,
        
        linewidth = 2
      )
    
    
    # IMPORTANT:
    # One combined colour scale for both
    # sex trajectories AND significance bars.
    
    combined_colour_values <- c(
      colour_values,
      sig_APOE2 =
        sig_bar_colors["sig_APOE2"],
      sig_APOE3 =
        sig_bar_colors["sig_APOE3"],
      sig_APOE4 =
        sig_bar_colors["sig_APOE4"]
    )
    
    
    p <- p +
      ggplot2::scale_colour_manual(
        values = combined_colour_values,
        guide = ggplot2::guide_legend(
          override.aes = list(
            shape = NA
          )
        )
      )
  }
  
  
  # ============================================================
  # 19. Axis limits
  # ============================================================
  
  final_xlim <- if (
    is.numeric(XLIM) &&
    length(XLIM) == 2
  ) {
    
    XLIM
    
  } else {
    
    c(
      min(
        df[[age_col]],
        na.rm = TRUE
      ),
      max(
        df[[age_col]],
        na.rm = TRUE
      )
    )
  }
  
  
  final_ylim <- if (
    is.numeric(YLIM) &&
    length(YLIM) == 2
  ) {
    
    YLIM
    
  } else {
    
    NULL
  }
  
  
  p <- p +
    ggplot2::coord_cartesian(
      xlim = final_xlim,
      ylim = final_ylim
    )
  
  
  # ============================================================
  # 20. Return
  # ============================================================
  
  return(
    list(
      plot = p,
      predictions = mod_df,
      significance_summary = significance_summary,
      model_obj = mod,
      z_value_used = z_crit,
      gam_check = check_output
    )
  )
}

run_apoe_sex_model_cohort <- function(
    df,
    apoe_col = "APOE_grouped",
    sex_col = "de_gender",
    y_col = "WUSTLcentiloid",
    age_col = "clinical_AgefromBaseline",
    intellectual_disability_col = "Intellectual_Disability",
    cohort = "cohort",
    family_choice = gaussian(),
    
    colour_values = c(
      "APOE2" = "#E41A1C",
      "APOE3" = "#377EB8",
      "APOE4" = "#4DAF4A"
    ),
    
    k_val = 4,
    reference_cohort = NULL,
    reference_sex = NULL,
    reference_id = NULL,
    
    XLIM = c(18, 72),
    YLIM = FALSE,
    
    # Shift displayed y values without changing the fitted model
    y_shift = 0,
    
    # Whether to display response-scale 95% CI ribbons
    show_ribbon = FALSE,
    ribbon_alpha = 0.20
) {
  
  library(mgcv)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(ggplot2)
  
  # ============================================================
  # 1. Check required columns
  # ============================================================
  
  required_cols <- c(
    y_col,
    age_col,
    apoe_col,
    sex_col
  )
  
  if (!is.null(cohort) && cohort %in% colnames(df)) {
    required_cols <- c(
      required_cols,
      cohort
    )
  }
  
  if (
    !is.null(intellectual_disability_col) &&
    intellectual_disability_col %in% colnames(df)
  ) {
    required_cols <- c(
      required_cols,
      intellectual_disability_col
    )
  }
  
  missing_cols <- setdiff(
    required_cols,
    colnames(df)
  )
  
  if (length(missing_cols) > 0) {
    stop(
      "The following required columns are missing from df: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  # ============================================================
  # 2. Fix column types
  # ============================================================
  
  df[[apoe_col]] <- factor(
    df[[apoe_col]],
    levels = c(
      "APOE2",
      "APOE3",
      "APOE4"
    )
  )
  
  df[[sex_col]] <- factor(
    df[[sex_col]]
  )
  
  df[[age_col]] <- as.numeric(
    df[[age_col]]
  )
  
  df[[y_col]] <- as.numeric(
    df[[y_col]]
  )
  
  # ============================================================
  # 3. Remove incomplete observations
  # ============================================================
  
  df <- df[
    complete.cases(
      df[, required_cols, drop = FALSE]
    ),
    ,
    drop = FALSE
  ]
  
  # ============================================================
  # 4. Validate APOE levels
  # ============================================================
  
  if (
    !all(
      c(
        "APOE2",
        "APOE3",
        "APOE4"
      ) %in%
      unique(
        as.character(
          df[[apoe_col]]
        )
      )
    )
  ) {
    stop(
      "APOE column must contain APOE2, APOE3, and APOE4."
    )
  }
  
  df[[apoe_col]] <- factor(
    df[[apoe_col]],
    levels = c(
      "APOE2",
      "APOE3",
      "APOE4"
    )
  )
  
  # ============================================================
  # 5. Validate sex levels
  # ============================================================
  
  sex_levels <- levels(
    droplevels(
      df[[sex_col]]
    )
  )
  
  if (length(sex_levels) != 2) {
    stop(
      "The function requires exactly two sex levels. Found: ",
      paste(
        sex_levels,
        collapse = ", "
      )
    )
  }
  
  df[[sex_col]] <- droplevels(
    df[[sex_col]]
  )
  
  # ============================================================
  # 6. Prepare cohort
  # ============================================================
  
  if (
    !is.null(cohort) &&
    cohort %in% colnames(df)
  ) {
    
    df[[cohort]] <- factor(
      df[[cohort]]
    )
    
    if (is.null(reference_cohort)) {
      reference_cohort <-
        levels(
          df[[cohort]]
        )[1]
    }
    
    if (
      !reference_cohort %in%
      levels(
        df[[cohort]]
      )
    ) {
      stop(
        "reference_cohort = '",
        reference_cohort,
        "' is not a valid level of ",
        cohort,
        ". Valid levels are: ",
        paste(
          levels(
            df[[cohort]]
          ),
          collapse = ", "
        )
      )
    }
  }
  
  # ============================================================
  # 7. Prepare reference sex
  # ============================================================
  
  if (is.null(reference_sex)) {
    reference_sex <- sex_levels[1]
  }
  
  if (
    !reference_sex %in%
    sex_levels
  ) {
    stop(
      "reference_sex = '",
      reference_sex,
      "' is not a valid level of ",
      sex_col,
      ". Valid levels are: ",
      paste(
        sex_levels,
        collapse = ", "
      )
    )
  }
  
  # ============================================================
  # 8. Prepare intellectual disability
  # ============================================================
  
  if (
    !is.null(intellectual_disability_col) &&
    intellectual_disability_col %in% colnames(df)
  ) {
    
    df[[intellectual_disability_col]] <-
      factor(
        df[[intellectual_disability_col]]
      )
    
    if (is.null(reference_id)) {
      reference_id <-
        levels(
          df[[intellectual_disability_col]]
        )[1]
    }
    
    if (
      !reference_id %in%
      levels(
        df[[intellectual_disability_col]]
      )
    ) {
      stop(
        "reference_id = '",
        reference_id,
        "' is not a valid level of ",
        intellectual_disability_col,
        ". Valid levels are: ",
        paste(
          levels(
            df[[intellectual_disability_col]]
          ),
          collapse = ", "
        )
      )
    }
  }
  
  # ============================================================
  # 9. Explicit APOE × sex interaction
  # ============================================================
  
  df$apoe_sex <- interaction(
    df[[apoe_col]],
    df[[sex_col]],
    sep = "_",
    drop = TRUE
  )
  
  # ============================================================
  # 10. Determine group-specific age bounds
  # ============================================================
  
  bounds_df <- df %>%
    dplyr::group_by(
      !!rlang::sym(apoe_col),
      !!rlang::sym(sex_col)
    ) %>%
    dplyr::summarise(
      
      min_valid_age =
        if (dplyr::n() >= 5) {
          sort(
            !!rlang::sym(age_col)
          )[5]
        } else {
          min(
            !!rlang::sym(age_col),
            na.rm = TRUE
          )
        },
      
      max_valid_age =
        if (dplyr::n() >= 5) {
          sort(
            !!rlang::sym(age_col),
            decreasing = TRUE
          )[5]
        } else {
          max(
            !!rlang::sym(age_col),
            na.rm = TRUE
          )
        },
      
      .groups = "drop"
    )
  
  # ============================================================
  # 11. Construct GAM formula
  # ============================================================
  
  formula_str <- paste0(
    y_col,
    " ~ apoe_sex",
    " + s(",
    age_col,
    ", by = apoe_sex, k = ",
    k_val,
    ")",
    " + s(",
    age_col,
    ", k = ",
    k_val,
    ")"
  )
  
  if (
    !is.null(cohort) &&
    cohort %in% colnames(df)
  ) {
    
    formula_str <- paste0(
      formula_str,
      " + ",
      cohort
    )
  }
  
  if (
    !is.null(intellectual_disability_col) &&
    intellectual_disability_col %in% colnames(df)
  ) {
    
    formula_str <- paste0(
      formula_str,
      " + ",
      intellectual_disability_col
    )
  }
  
  formula_gam <- as.formula(
    formula_str
  )
  
  # ============================================================
  # 12. Fit GAM
  # ============================================================
  
  mod <- mgcv::gam(
    formula_gam,
    method = "REML",
    data = df,
    family = family_choice
  )
  
  # ============================================================
  # 13. Capture gam.check output
  # ============================================================
  
  check_output <- capture.output(
    mgcv::gam.check(
      mod,
      plot = FALSE
    )
  )
  
  # ============================================================
  # 14. Build prediction grid
  # ============================================================
  
  age_grid <- seq(
    min(
      df[[age_col]],
      na.rm = TRUE
    ),
    max(
      df[[age_col]],
      na.rm = TRUE
    ),
    by = 0.1
  )
  
  grid_list <- list(
    age_grid,
    levels(
      df[[apoe_col]]
    ),
    levels(
      df[[sex_col]]
    )
  )
  
  grid_names <- c(
    age_col,
    apoe_col,
    sex_col
  )
  
  if (
    !is.null(cohort) &&
    cohort %in% colnames(df)
  ) {
    
    grid_list <- c(
      grid_list,
      list(reference_cohort)
    )
    
    grid_names <- c(
      grid_names,
      cohort
    )
  }
  
  if (
    !is.null(intellectual_disability_col) &&
    intellectual_disability_col %in% colnames(df)
  ) {
    
    grid_list <- c(
      grid_list,
      list(reference_id)
    )
    
    grid_names <- c(
      grid_names,
      intellectual_disability_col
    )
  }
  
  mod_df <- expand.grid(
    grid_list,
    stringsAsFactors = FALSE
  )
  
  colnames(mod_df) <- grid_names
  
  # Restore factor levels
  mod_df[[apoe_col]] <- factor(
    mod_df[[apoe_col]],
    levels = levels(
      df[[apoe_col]]
    )
  )
  
  mod_df[[sex_col]] <- factor(
    mod_df[[sex_col]],
    levels = levels(
      df[[sex_col]]
    )
  )
  
  if (
    !is.null(cohort) &&
    cohort %in% colnames(df)
  ) {
    
    mod_df[[cohort]] <- factor(
      mod_df[[cohort]],
      levels = levels(
        df[[cohort]]
      )
    )
  }
  
  if (
    !is.null(intellectual_disability_col) &&
    intellectual_disability_col %in% colnames(df)
  ) {
    
    mod_df[[intellectual_disability_col]] <-
      factor(
        mod_df[[intellectual_disability_col]],
        levels =
          levels(
            df[[intellectual_disability_col]]
          )
      )
  }
  
  mod_df$apoe_sex <- interaction(
    mod_df[[apoe_col]],
    mod_df[[sex_col]],
    sep = "_",
    drop = TRUE
  )
  
  # ============================================================
  # 15. Apply age-density constraints BEFORE prediction
  # ============================================================
  
  mod_df <- mod_df %>%
    dplyr::left_join(
      bounds_df,
      by = c(
        apoe_col,
        sex_col
      )
    ) %>%
    dplyr::filter(
      !!rlang::sym(age_col) >= min_valid_age,
      !!rlang::sym(age_col) <= max_valid_age
    )
  
  if (nrow(mod_df) == 0) {
    stop(
      "No prediction points remain after applying ",
      "the APOE × sex-specific age bounds."
    )
  }
  
  # ============================================================
  # 16. Link-scale prediction
  # ============================================================
  
  pred_link <- predict(
    mod,
    newdata = mod_df,
    type = "link",
    se.fit = TRUE
  )
  
  # Pointwise 95% CI
  z <- qnorm(
    0.975
  )
  
  mod_df$eta <- as.numeric(
    pred_link$fit
  )
  
  mod_df$se_eta <- as.numeric(
    pred_link$se.fit
  )
  
  # ============================================================
  # 17. Convert predictions to response scale
  # ============================================================
  
  linkinv <- mod$family$linkinv
  
  mu_eta <- mod$family$mu.eta(
    mod_df$eta
  )
  
  mod_df$forecast <- linkinv(
    mod_df$eta
  )
  
  # Delta-method response-scale SE
  mod_df$se <- abs(
    mu_eta
  ) * mod_df$se_eta
  
  # Correct response-scale CI for nonlinear links
  mod_df$min <- linkinv(
    mod_df$eta -
      z * mod_df$se_eta
  )
  
  mod_df$max <- linkinv(
    mod_df$eta +
      z * mod_df$se_eta
  )
  
  # ============================================================
  # 18. Full prediction covariance
  # ============================================================
  
  Xp <- predict(
    mod,
    newdata = mod_df,
    type = "lpmatrix"
  )
  
  Vp <- vcov(
    mod
  )
  
  pred_cov_link <- Xp %*%
    Vp %*%
    t(Xp)
  
  # ============================================================
  # 19. Significant age ranges
  #
  # Tests:
  #
  #   M vs F within APOE2
  #   M vs F within APOE3
  #   M vs F within APOE4
  #
  # = 3 planned comparisons total.
  #
  # SEs use the FULL prediction covariance matrix.
  # ============================================================
  
  significance_list <- list()
  
  n_comparisons <- 3
  
  alpha_adj <- 0.05 / n_comparisons
  
  z_crit <- qnorm(
    1 - alpha_adj / 2
  )
  
  # ------------------------------------------------------------
  # Determine which sex is the reference and which is the
  # comparison. The contrast is:
  #
  #   reference_sex - other_sex
  #
  # If reference_sex = "M", this is M - F.
  # ------------------------------------------------------------
  
  comparison_sex <- setdiff(
    sex_levels,
    reference_sex
  )
  
  if (length(comparison_sex) != 1) {
    stop(
      "Could not determine comparison sex."
    )
  }
  
  comparison_sex <- comparison_sex[1]
  
  # ------------------------------------------------------------
  # Test each APOE group separately
  # ------------------------------------------------------------
  
  for (
    apoe_val in c(
      "APOE2",
      "APOE3",
      "APOE4"
    )
  ) {
    
    # Rows for this APOE group and both sexes
    apoe_rows <- which(
      mod_df[[apoe_col]] == apoe_val
    )
    
    if (length(apoe_rows) == 0) {
      next
    }
    
    # Ages where both sexes have predictions
    ages_ref <- unique(
      mod_df[[age_col]][
        apoe_rows &
          mod_df[[sex_col]] == reference_sex
      ]
    )
    
    ages_comp <- unique(
      mod_df[[age_col]][
        apoe_rows &
          mod_df[[sex_col]] == comparison_sex
      ]
    )
    
    common_ages <- intersect(
      ages_ref,
      ages_comp
    )
    
    if (length(common_ages) == 0) {
      next
    }
    
    # ----------------------------------------------------------
    # Calculate M vs F contrast at each age
    # ----------------------------------------------------------
    
    contrast_results <- lapply(
      common_ages,
      function(age_i) {
        
        ref_row <- which(
          mod_df[[apoe_col]] == apoe_val &
            mod_df[[sex_col]] == reference_sex &
            mod_df[[age_col]] == age_i
        )
        
        comp_row <- which(
          mod_df[[apoe_col]] == apoe_val &
            mod_df[[sex_col]] == comparison_sex &
            mod_df[[age_col]] == age_i
        )
        
        if (
          length(ref_row) != 1 ||
          length(comp_row) != 1
        ) {
          return(NULL)
        }
        
        r <- c(
          ref_row,
          comp_row
        )
        
        # ------------------------------------------------------
        # Response-scale predictions
        # ------------------------------------------------------
        
        mu_i <- mod_df$forecast[r]
        
        # ------------------------------------------------------
        # Link-scale covariance
        # ------------------------------------------------------
        
        Sigma_eta <-
          pred_cov_link[
            r,
            r,
            drop = FALSE
          ]
        
        # ------------------------------------------------------
        # Transform covariance to response scale
        # ------------------------------------------------------
        
        eta_i <- mod_df$eta[r]
        
        deriv_i <-
          mod$family$mu.eta(
            eta_i
          )
        
        J <- diag(
          as.numeric(
            deriv_i
          )
        )
        
        Sigma_mu <- J %*%
          Sigma_eta %*%
          J
        
        # ------------------------------------------------------
        # Contrast:
        #
        # reference sex - comparison sex
        #
        # e.g. M - F
        # ------------------------------------------------------
        
        contrast_vector <- c(
          1,
          -1
        )
        
        diff <- sum(
          contrast_vector *
            mu_i
        )
        
        se_diff <- sqrt(
          max(
            0,
            as.numeric(
              t(
                contrast_vector
              ) %*%
                Sigma_mu %*%
                contrast_vector
            )
          )
        )
        
        data.frame(
          age = age_i,
          diff = diff,
          se = se_diff
        )
      }
    )
    
    contrast_results <-
      dplyr::bind_rows(
        contrast_results
      )
    
    if (
      nrow(
        contrast_results
      ) == 0
    ) {
      next
    }
    
    names(
      contrast_results
    )[
      names(
        contrast_results
      ) == "age"
    ] <- age_col
    
    DT_tmp <-
      data.table::as.data.table(
        contrast_results
      )
    
    # ----------------------------------------------------------
    # Bonferroni-adjusted CI
    # ----------------------------------------------------------
    
    DT_tmp[
      ,
      `:=`(
        diff_low =
          diff -
          z_crit * se,
        
        diff_hi =
          diff +
          z_crit * se
      )
    ]
    
    # ----------------------------------------------------------
    # Significant if adjusted CI excludes zero
    # ----------------------------------------------------------
    
    DT_tmp[
      ,
      sig_diff :=
        diff_low > 0 |
        diff_hi < 0
    ]
    
    # ----------------------------------------------------------
    # Convert significant points into age ranges
    # ----------------------------------------------------------
    
    ranges <- get_sig_ranges(
      DT_tmp,
      "sig_diff",
      age_col
    )
    
    if (
      !is.null(ranges) &&
      nrow(ranges) > 0
    ) {
      
      ranges <- data.table::as.data.table(
        ranges
      )
      
      ranges[,
        contrast :=
          paste0(
            reference_sex,
            "_vs_",
            comparison_sex
          )]
      
      ranges[,
        (apoe_col) :=
          apoe_val]
      
      ranges[,
        (sex_col) :=
          paste0(
            reference_sex,
            "_vs_",
            comparison_sex
          )]
      
      significance_list[[length(significance_list) + 1]] <- ranges
    }
  }
  
  # ------------------------------------------------------------
  # Combine all significant ranges
  # ------------------------------------------------------------
  
  if (
    length(
      significance_list
    ) > 0
  ) {
    
    significance_summary <-
      data.table::rbindlist(
        significance_list,
        fill = TRUE
      )
    
  } else {
    
    significance_summary <-
      data.table::data.table()
  }
  
  # ============================================================
  # 20. Plot-scale shifted predictions
  # ============================================================
  
  mod_df$forecast_plot <-
    mod_df$forecast +
    y_shift
  
  mod_df$min_plot <-
    mod_df$min +
    y_shift
  
  mod_df$max_plot <-
    mod_df$max +
    y_shift
  
  df_plot <- df
  
  df_plot$y_plot <-
    df_plot[[y_col]] +
    y_shift
  
  # ============================================================
  # 21. Plot
  # ============================================================
  
  linetype_vector <- c(
    "M" = "solid",
    "Male" = "solid",
    "male" = "solid",
    "F" = "dashed",
    "Female" = "dashed",
    "female" = "dashed"
  )
  
  linetype_vector <-
    linetype_vector[
      names(linetype_vector) %in%
        levels(
          df[[sex_col]]
        )
    ]
  
  plot_aes <- ggplot2::aes(
    x = !!rlang::sym(age_col),
    y = y_plot,
    colour = !!rlang::sym(apoe_col)
  )
  
  if (
    !is.null(cohort) &&
    cohort %in% colnames(df)
  ) {
    
    plot_aes$shape <-
      rlang::sym(cohort)
  }
  
  p <- ggplot2::ggplot(
    df_plot,
    plot_aes
  ) +
    
    ggplot2::geom_point(
      alpha = 0.2
    )
  
  # ============================================================
  # 22. Optional confidence ribbons
  # ============================================================
  
  if (show_ribbon) {
    
    p <- p +
      ggplot2::geom_ribbon(
        data = mod_df,
        ggplot2::aes(
          x = !!rlang::sym(age_col),
          ymin = min_plot,
          ymax = max_plot,
          fill = !!rlang::sym(apoe_col),
          group = apoe_sex
        ),
        alpha = ribbon_alpha,
        inherit.aes = FALSE
      )
  }
  
  # ============================================================
  # 23. Fitted trajectories
  # ============================================================
  
  p <- p +
    
    ggplot2::geom_line(
      data = mod_df,
      ggplot2::aes(
        x = !!rlang::sym(age_col),
        y = forecast_plot,
        group = apoe_sex,
        colour = !!rlang::sym(apoe_col),
        linetype = !!rlang::sym(sex_col)
      ),
      linewidth = 1,
      inherit.aes = FALSE
    ) +
    
    ggplot2::theme_bw() +
    
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "horizontal"
    ) +
    
    ggplot2::scale_colour_manual(
      values = colour_values
    ) +
    
    ggplot2::scale_linetype_manual(
      values = linetype_vector
    ) +
    
    ggplot2::scale_x_continuous(
      breaks = c(
        20,
        30,
        40,
        50,
        60,
        70
      )
    ) +
    
    ggplot2::xlab(
      "Age"
    )
  
  if (show_ribbon) {
    
    p <- p +
      ggplot2::scale_fill_manual(
        values = colour_values,
        guide = "none"
      )
  }
  
  if (
    !is.null(cohort) &&
    cohort %in% colnames(df)
  ) {
    
    p <- p +
      ggplot2::scale_shape_manual(
        values = c(
          1,
          4
        )
      )
  }
  
  # ============================================================
  # 24. Coordinate viewport
  #
  # YLIM is applied to the DISPLAYED scale after y_shift.
  # ============================================================
  
  final_xlim <- if (
    is.numeric(XLIM) &&
    length(XLIM) == 2
  ) {
    
    XLIM
    
  } else {
    
    c(
      min(
        df[[age_col]],
        na.rm = TRUE
      ),
      max(
        df[[age_col]],
        na.rm = TRUE
      )
    )
  }
  
  final_ylim <- if (
    is.numeric(YLIM) &&
    length(YLIM) == 2
  ) {
    
    YLIM
    
  } else {
    
    NULL
  }
  
  p <- p +
    ggplot2::coord_cartesian(
      xlim = final_xlim,
      ylim = final_ylim
    )
  
  # ============================================================
  # 25. Return
  # ============================================================
  
  return(
    list(
      plot = p,
      
      # Response-scale predictions
      predictions = mod_df,
      
      # Original data plus shifted plotting outcome
      plotting_data = df_plot,
      
      # Fitted model
      model_obj = mod,
      
      # Full link-scale prediction covariance
      prediction_covariance_link =
        pred_cov_link,
      
      # Significant age ranges
      significance_summary =
        significance_summary,
      
      # Pointwise 95% CI z value
      z_value_used = z,
      
      # Bonferroni-adjusted z value used for
      # significance testing
      z_value_significance =
        z_crit,
      
      # Model diagnostics
      gam_check = check_output
    )
  )
}