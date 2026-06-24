library(mgcv)
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
# ---- Identify continuous ranges of significance ----
get_sig_ranges <- function(DT, col, age_col) {
  DT2 <- DT[get(col) == TRUE, .(age_val = get(age_col))]
  if (nrow(DT2) == 0) return(NULL)
  
  diffs <- diff(DT2$age_val)
  step <- as.numeric(names(sort(table(round(diffs, 6)), decreasing = TRUE)[1]))
  if (is.na(step) || step <= 0) step <- median(diffs)
  
  tol <- step * 0.01
  DT2[, lag_age := shift(age_val)]
  DT2[, grp := cumsum(abs(age_val - lag_age - step) > tol | is.na(lag_age))]
  
  DT2[, .(
    contrast = col,
    start_age = min(age_val),
    end_age   = max(age_val),
    n_points  = .N
  ), by = grp]
}


run_apoe_model <- function(
    df,
    apoe_col = "APOE_grouped",
    y_col = "WUSTLcentiloid",
    age_col = "clinical_AgefromBaseline",
    cohort = "cohort",
    fill_values,
    colour_values,
    sig_bar_colors = c("sig_2_3" = "red", "sig_2_4" = "blue", "sig_3_4" = "green"),
    markSig = FALSE,
    k_val = 4, 
    scatter = TRUE,
    YLIM = FALSE
) {
  # ---- Fix column types ----
  df[[apoe_col]] <- factor(df[[apoe_col]], levels = c("APOE2", "APOE3", "APOE4"))
  df[[age_col]] <- as.numeric(df[[age_col]]) 
  
  # ---- Determine Group-Specific Age Bounds (N >= 5 at limits) ----
  # Calculate the age limits for each group where at least 5 points exist 
  # above and below the boundary.
  bounds_df <- df %>%
    group_by(!!sym(apoe_col)) %>%
    summarise(
      min_valid_age = sort(!!sym(age_col))[5],
      max_valid_age = sort(!!sym(age_col), decreasing = TRUE)[5],
      .groups = 'drop'
    )
  
  # ---- Fit GAM ----
  formula_gam <- as.formula(
    paste0(y_col, " ~ s(", age_col, ", by = ", apoe_col, ", k = ", k_val, ")")
  )
  mod <- gam(formula_gam, method = "REML", data = df)
  
  # ---- Prediction Grid ----
  mod_df <- expand.grid(
    seq(min(df[[age_col]], na.rm = TRUE), max(df[[age_col]], na.rm = TRUE), by = 0.1),
    levels(df[[apoe_col]])
  )
  colnames(mod_df) <- c(age_col, apoe_col)
  
  tmp <- predict(mod, newdata = mod_df, se.fit = TRUE)
  mod_df$forecast <- tmp$fit
  mod_df$se        <- tmp$se.fit
  mod_df$min      <- mod_df$forecast - 1.96 * mod_df$se
  mod_df$max      <- mod_df$forecast + 1.96 * mod_df$se
  
  # ---- Apply Constraints ----
  
  # 1. Truncate model fits based on the N=5 data density rule
  mod_df <- mod_df %>%
    left_join(bounds_df, by = apoe_col) %>%
    filter(!!sym(age_col) >= min_valid_age & !!sym(age_col) <= max_valid_age)
  
  # 2. Truncate confidence intervals to YLIM to prevent clipping/removal
  if (is.numeric(YLIM) && length(YLIM) == 2) {
    mod_df$min <- pmax(mod_df$min, YLIM[1])
    mod_df$max <- pmin(mod_df$max, YLIM[2])
  }
  
  # ---- Pivot wider for differences ----
  df_wide <- mod_df %>%
    dplyr::select(all_of(c(age_col, apoe_col, "forecast", "se"))) %>%
    pivot_wider(names_from = all_of(apoe_col),
                values_from = c("forecast", "se"))
  
  DT <- as.data.table(df_wide)
  
  # ---- Compute pairwise differences and significance ----
  #Justification for not using CI overlap but instd computing difference:
  # You cannot determine significance by comparing whether two 95% CIs overlap.
  # You must compute the CI for the difference between groups.
  # Cumming & Finch (2005) — CI overlap is not a correct significance test
  # 
  # Cumming, G., & Finch, S. (2005). "Inference by eye: Confidence intervals and how to read pictures of data." American Psychologist, 60(2), 170–180.
  # https://doi.org/10.1037/0003-066X.60.2.170
  # 
  # Relevant points from the paper:
  #   
  #   “Checking whether two independent 95% CIs overlap is a conservative test of the hypothesis that the two means differ.”
  # 
  # “Two 95% CIs can overlap by as much as ~25% and still correspond to a statistically significant difference (p < .05).”
  # 
  # “For accurate inference, construct a CI for the difference between the means.”
  # 
  # This is directly relevant to comparing two models’ MAE estimates.
  # 
  # **2. Efron & Tibshirani (1994) — Bootstrap the difference
  # 
  # (Essential reference for bootstrap hypothesis testing)**
  #   Efron, B., & Tibshirani, R. (1994). An Introduction to the Bootstrap. CRC Press.
  # 
  # Key idea from the book:
  #   
  #   When comparing two estimators , the bootstrap distribution of the difference 
  # is the correct basis for inference.
  # 
  # Quoting: “The bootstrap estimate of the standard error of a difference in parameters is obtained by resampling the paired differences.”
  # 
  # This directly supports using paired bootstrap differences of MAE rather than CI overlap.
  # 
  # 3. Altman & Bland (2011) — CI overlap is not a hypothesis test
  # 
  # Altman, D. G., & Bland, J. M. (2011). "How to compare two means." BMJ, 343, d2304.
  # https://doi.org/10.1136/bmj.d2304
  # 
  # They explain:
  #   
  #   “The common practice of examining whether the confidence intervals overlap is conservative and can be misleading.”
  # 
  # “The proper comparison is a confidence interval for the difference between means, not separate intervals.”
  # 
  # This applies directly to your situation with two MAE curves.
  # 
  # 4. Gelman & Hill (2006) — CI overlap is not equivalent to a test of differences
  # 
  # Gelman, A., & Hill, J. (2006). Data Analysis Using Regression and Multilevel/Hierarchical Models. Cambridge University Press.
  # 
  # They emphasize:
  #   
  #   “Confidence intervals for group-specific estimates are not a substitute for intervals on their difference.”
  
  # Ensure necessary columns exist before DT operations
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
      diff_2_3_low = diff_2_3 - 1.96 * se_2_3,
      diff_2_3_hi  = diff_2_3 + 1.96 * se_2_3,
      
      diff_2_4_low = diff_2_4 - 1.96 * se_2_4,
      diff_2_4_hi  = diff_2_4 + 1.96 * se_2_4,
      
      diff_3_4_low = diff_3_4 - 1.96 * se_3_4,
      diff_3_4_hi  = diff_3_4 + 1.96 * se_3_4
    )]
    
    DT[, `:=`(
      sig_2_3 = diff_2_3_low > 0 | diff_2_3_hi < 0,
      sig_2_4 = diff_2_4_low > 0 | diff_2_4_hi < 0,
      sig_3_4 = diff_3_4_low > 0 | diff_3_4_hi < 0
    )]
    
    significance_summary <- rbindlist(list(
      get_sig_ranges(DT, "sig_2_3", age_col),
      get_sig_ranges(DT, "sig_2_4", age_col),
      get_sig_ranges(DT, "sig_3_4", age_col)
    ), fill = TRUE)
  } else {
    significance_summary <- data.table()
  }
  
  # ---- Plot ----
  p <- ggplot(df, aes(
    x = !!sym(age_col),
    y = !!sym(y_col),
    colour = !!sym(apoe_col),
    shape = !!sym(cohort),
    group = !!sym(apoe_col)
  )) +
    geom_ribbon(
      data = mod_df,
      aes(
        x = !!sym(age_col),
        y = forecast,
        ymin = min,
        ymax = max,
        fill = !!sym(apoe_col)
      ),
      alpha = 0.4,
      inherit.aes = FALSE
    ) +
    geom_line(
      data = mod_df,
      aes(
        x = !!sym(age_col),
        y = forecast,
        group = !!sym(apoe_col)
      )
    ) +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_fill_manual(values = fill_values) +
    scale_colour_manual(values = colour_values) +
    scale_shape_manual(values = c(1, 3)) +
    coord_cartesian(xlim = c(min(df[[age_col]], na.rm = TRUE), max(df[[age_col]], na.rm = TRUE)))
  
  if(scatter == TRUE){
    p <- p + geom_point(data = df, aes(
      x = !!sym(age_col),
      y = !!sym(y_col),
      colour = !!sym(apoe_col),
      shape = !!sym(cohort),
      group = !!sym(apoe_col)
    ), alpha = 0.3) 
  }
  
  if(is.numeric(YLIM) && length(YLIM) == 2){
    p <- p + coord_cartesian(ylim = YLIM)
  }
  
  # ---- Add colored bars for significance ----
  if (nrow(significance_summary) > 0 & markSig == TRUE) {
    
    # Use specified YLIM or data range for placement
    y_range_max <- if(is.numeric(YLIM)) YLIM[2] else max(df[[y_col]], na.rm = TRUE)
    y_range_min <- if(is.numeric(YLIM)) YLIM[1] else min(df[[y_col]], na.rm = TRUE)
    
    y_buffer <- (y_range_max - y_range_min) * 0.05 
    contrast_levels <- unique(significance_summary$contrast)
    y_positions <- setNames(seq(from = y_range_max + y_buffer, 
                                by = y_buffer, 
                                length.out = length(contrast_levels)), 
                            contrast_levels)
    
    significance_summary[, y := y_positions[contrast]]
    
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
                   sig_2_3 = sig_bar_colors["sig_2_3"],
                   sig_2_4 = sig_bar_colors["sig_2_4"],
                   sig_3_4 = sig_bar_colors["sig_3_4"]),
        guide = guide_legend(override.aes = list(shape = NA))
      )
  }
  
  return(list(
    plot = p,
    predictions = mod_df,
    significance_summary = significance_summary
  ))
}

run_apoe_model_24s <- function(
    df,
    apoe_col = "APOE_grouped",
    y_col = "WUSTLcentiloid",
    age_col = "clinical_AgefromBaseline",
    cohort = "cohort",
    fill_values,
    colour_values,
    sig_bar_colors = c(
      "sig_2_3"  = "red",
      "sig_2_24" = "orange",
      "sig_2_4"  = "blue",
      "sig_24_3" = "purple",
      "sig_24_4" = "brown",
      "sig_3_4"  = "green"
    ),
    markSig = FALSE,
    k_val = 4
) {
  
  ## ---- Fix column types ----
  apoe_levels <- c("APOE2", "APOE24", "APOE3", "APOE4")
  df[[apoe_col]] <- factor(df[[apoe_col]], levels = apoe_levels)
  df[[age_col]]  <- as.numeric(df[[age_col]])
  
  ## ---- Fit GAM ----
  formula_gam <- as.formula(
    paste0(y_col, " ~ s(", age_col, ", by = ", apoe_col, ", k = ", k_val, ")")
  )
  mod <- gam(formula_gam, method = "REML", data = df)
  
  ## ---- Prediction Grid ----
  mod_df <- expand.grid(
    seq(min(df[[age_col]]), max(df[[age_col]]), by = 0.1),
    apoe_levels
  )
  colnames(mod_df) <- c(age_col, apoe_col)
  
  tmp <- predict(mod, newdata = mod_df, se.fit = TRUE)
  mod_df$forecast <- tmp$fit
  mod_df$se       <- tmp$se.fit
  mod_df$min      <- mod_df$forecast - 1.96 * mod_df$se
  mod_df$max      <- mod_df$forecast + 1.96 * mod_df$se
  
  ## ---- Pivot wider ----
  df_wide <- mod_df %>%
    dplyr::select(all_of(c(age_col, apoe_col, "forecast", "se"))) %>%
    tidyr::pivot_wider(
      names_from  = all_of(apoe_col),
      values_from = c("forecast", "se")
    )
  
  DT <- data.table::as.data.table(df_wide)
  
  ## ---- Define contrasts ----
  contrasts <- list(
    c("APOE2",  "APOE3"),
    c("APOE2",  "APOE24"),
    c("APOE2",  "APOE4"),
    c("APOE24", "APOE3"),
    c("APOE24", "APOE4"),
    c("APOE3",  "APOE4")
  )
  
  ## ---- Compute differences + significance ----
  for (ct in contrasts) {
    g1 <- ct[1]
    g2 <- ct[2]
    tag <- paste0(
      sub("APOE", "", g1), "_",
      sub("APOE", "", g2)
    )
    
    DT[, paste0("diff_", tag) :=
         get(paste0("forecast_", g1)) -
         get(paste0("forecast_", g2))]
    
    DT[, paste0("se_", tag) :=
         sqrt(get(paste0("se_", g1))^2 +
                get(paste0("se_", g2))^2)]
    
    DT[, paste0("low_", tag) :=
         get(paste0("diff_", tag)) -
         1.96 * get(paste0("se_", tag))]
    
    DT[, paste0("hi_", tag) :=
         get(paste0("diff_", tag)) +
         1.96 * get(paste0("se_", tag))]
    
    DT[, paste0("sig_", tag) :=
         get(paste0("low_", tag)) > 0 |
         get(paste0("hi_", tag))  < 0]
  }
  
  ## ---- Collect significance ranges ----
  significance_summary <- rbindlist(
    lapply(
      paste0("sig_", sapply(contrasts, function(x)
        paste0(sub("APOE","",x[1]), "_", sub("APOE","",x[2])))),
      function(s) get_sig_ranges(DT, s, age_col)
    ),
    fill = TRUE
  )
  
  ## ---- Plot ----
  p <- ggplot(df, aes(
    x = !!sym(age_col),
    y = !!sym(y_col),
    colour = !!sym(apoe_col),
    shape = !!sym(cohort),
    group  = !!sym(apoe_col)
  )) +
    geom_point(alpha = 0.3) +
    geom_ribbon(
      data = mod_df,
      aes(
        x = !!sym(age_col),
        y = forecast,
        ymin = min,
        ymax = max,
        fill = !!sym(apoe_col)
      ),
      alpha = 0.4
    ) +
    geom_line(
      data = mod_df,
      aes(x = !!sym(age_col), y = forecast)
    ) +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_fill_manual(values = fill_values) +
    scale_colour_manual(values = colour_values) +
    scale_shape_manual(values = c(1, 3))
  
  ## ---- Significance bars ----
  if (markSig && nrow(significance_summary) > 0) {
    
    y_max <- max(df[[y_col]], na.rm = TRUE)
    y_buffer <- diff(range(df[[y_col]], na.rm = TRUE)) * 0.05
    
    contrast_levels <- unique(significance_summary$contrast)
    y_positions <- setNames(
      seq(y_max + y_buffer,
          by = y_buffer,
          length.out = length(contrast_levels)),
      contrast_levels
    )
    
    significance_summary[, y := y_positions[contrast]]
    
    p <- p +
      geom_segment(
        data = significance_summary,
        aes(
          x = start_age,
          xend = end_age,
          y = y,
          yend = y,
          colour = contrast
        ),
        inherit.aes = FALSE,
        size = 2
      ) +
      scale_colour_manual(
        values = c(colour_values, sig_bar_colors),
        guide  = guide_legend(override.aes = list(shape = NA))
      )
  }
  
  return(list(
    plot = p,
    predictions = mod_df,
    significance_summary = significance_summary
  ))
}



run_group_model <- function(
    df,
    group_col = "group_col",      # now 6 groups instead of 3
    y_col = "WUSTLcentiloid",
    age_col = "clinical_AgefromBaseline",
    cohort = "cohort",
    fill_values = c("#FDAA9F", "#66D98E", "#A6C9FF",
                    "#F8766D", "#00BA38", "#619CFF"),
    colour_values = c("#FDAA9F", "#66D98E", "#A6C9FF",
                      "#F8766D", "#00BA38", "#619CFF"),
    linetype_values = c("solid", "solid", "solid",
                        "dashed", "dashed", "dashed"),
    k_val = 4
) {
  
  # ---- Fix column types ----
  df[[group_col]] <- factor(df[[group_col]])
  df[[age_col]]   <- as.numeric(df[[age_col]])
  
  # ---- Fit GAM ----
  formula_gam <- as.formula(
    paste0(y_col, " ~ s(", age_col, ", by = ", group_col, ", k = ", k_val, ")")
  )
  mod <- mgcv::gam(formula_gam, method = "REML", data = df)
  
  # ---- Prediction Grid ----
  mod_df <- expand.grid(
    age = seq(min(df[[age_col]]), max(df[[age_col]]), by = 0.1),
    grp = levels(df[[group_col]])
  )
  names(mod_df) <- c(age_col, group_col)
  
  tmp <- predict(mod, newdata = mod_df, se.fit = TRUE)
  mod_df$forecast <- tmp$fit
  mod_df$se       <- tmp$se.fit
  mod_df$min      <- mod_df$forecast - 1.96 * mod_df$se
  mod_df$max      <- mod_df$forecast + 1.96 * mod_df$se
  
  # ---- Plot ----
  p <- ggplot(df, aes(
    x = !!sym(age_col),
    y = !!sym(y_col),
    colour = !!sym(group_col),
    shape = !!sym(cohort),
    linetype = !!sym(group_col),
    group = !!sym(group_col)
  )) +
    geom_point(alpha = 0.3) +
    # geom_ribbon(
    #   data = mod_df,
    #   aes(
    #     x = !!sym(age_col),
    #     y = forecast,
    #     ymin = min,
    #     ymax = max,
    #     fill = !!sym(group_col)
    #   ),
    #   alpha = 0.35
    # ) +
    geom_line(
      data = mod_df,
      aes(
        x = !!sym(age_col),
        y = forecast,
        group = !!sym(group_col)
      ),
      linewidth = 1
    ) +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_fill_manual(values = fill_values) +
    scale_colour_manual(values = colour_values) +
    scale_linetype_manual(values = linetype_values) +
    scale_shape_manual(values = c(1, 4, 4)) +
    xlim(c(min(df[[age_col]]), max(df[[age_col]])))
  
  return(p)
}

run_group_sex_model <- function(
    df,
    group_col = "group_col",
    y_col = "WUSTLcentiloid",
    age_col = "clinical_AgefromBaseline",
    cohort = "cohort",   # set to NULL to disable cohort control
    fill_values = c("#FDAA9F", "#66D98E", "#A6C9FF",
                    "#F8766D", "#00BA38", "#619CFF"),
    colour_values = c("#FDAA9F", "#66D98E", "#A6C9FF",
                      "#F8766D", "#00BA38", "#619CFF"),
    linetype_values = c("solid","solid","solid",
                        "dashed","dashed","dashed"),
    k_val = 4,
    alpha = 0.05,
    YLIM = FALSE
) {
  
  library(data.table)
  library(mgcv)
  library(ggplot2)
  
  # ---- Prep ----
  df[[group_col]] <- factor(
    df[[group_col]],
    levels = c("APOE2F","APOE2M",
               "APOE3F","APOE3M",
               "APOE4F","APOE4M")
  )
  
  df[[age_col]] <- as.numeric(df[[age_col]])
  
  # ---- Optional cohort handling ----
  use_cohort <- !is.null(cohort)
  
  if (use_cohort) {
    df[[cohort]] <- droplevels(factor(df[[cohort]]))
    has_cohort <- nlevels(df[[cohort]]) > 1
  } else {
    has_cohort <- FALSE
  }
  
  # ---- Build GAM formula ----
  if (has_cohort) {
    formula_gam <- as.formula(
      paste0(
        y_col,
        " ~ ",
        cohort,
        " + s(", age_col, ", by = ", group_col, ", k = ", k_val, ")"
      )
    )
  } else {
    formula_gam <- as.formula(
      paste0(
        y_col,
        " ~ s(", age_col, ", by = ", group_col, ", k = ", k_val, ")"
      )
    )
  }

  
  # ---- Fit GAM ----
  mod <- gam(
    formula_gam,
    method = "REML",
    data = df
  )
  
  # ---- Prediction grid ----
  if (has_cohort) {
    ref_cohort <- levels(df[[cohort]])[1]
    
    mod_df <- expand.grid(
      age = seq(min(df[[age_col]]), max(df[[age_col]]), by = 0.1),
      grp = levels(df[[group_col]]),
      cohort = ref_cohort
    )
    names(mod_df) <- c(age_col, group_col, cohort)
  } else {
    mod_df <- expand.grid(
      age = seq(min(df[[age_col]]), max(df[[age_col]]), by = 0.1),
      grp = levels(df[[group_col]])
    )
    names(mod_df) <- c(age_col, group_col)
  }
  
  pred <- predict(mod, mod_df, se.fit = TRUE)
  mod_df$fit <- pred$fit
  mod_df$se  <- pred$se.fit
  
  # ---- Wide format for contrasts ----
  DT <- as.data.table(mod_df)
  form <- as.formula(paste(age_col, "~", group_col))
  
  W <- dcast(
    DT,
    formula = form,
    value.var = c("fit", "se")
  )
  
  # ---- Define contrasts ----
  contrasts <- list(
    sig_APOE2 = c("APOE2F","APOE2M"),
    sig_APOE3 = c("APOE3F","APOE3M"),
    sig_APOE4 = c("APOE4F","APOE4M")
  )
  
  # ---- Compute differences, SEs, p-values ----
  for (nm in names(contrasts)) {
    gF <- contrasts[[nm]][1]
    gM <- contrasts[[nm]][2]
    
    W[, paste0("diff_", nm) :=
        get(paste0("fit_", gF)) - get(paste0("fit_", gM))]
    
    W[, paste0("se_", nm) :=
        sqrt(get(paste0("se_", gF))^2 +
               get(paste0("se_", gM))^2)]
    
    W[, paste0("z_", nm) :=
        get(paste0("diff_", nm)) /
        get(paste0("se_", nm))]
    
    W[, paste0("p_", nm) :=
        2 * pnorm(-abs(get(paste0("z_", nm))))]
  }
  
  # ---- Holm correction ----
  pcols <- paste0("p_", names(contrasts))
  W[, (pcols) := lapply(.SD, p.adjust, method = "holm"),
    .SDcols = pcols]
  
  for (nm in names(contrasts)) {
    W[, (nm) := get(paste0("p_", nm)) < alpha]
  }
  
  # ---- Extract continuous age ranges ----
  get_ranges <- function(DT, flag, age_col) {
    X <- DT[get(flag) == TRUE, .(age = get(age_col))]
    if (nrow(X) == 0) return(NULL)
    
    step <- median(diff(X$age))
    X[, grp := cumsum(c(TRUE, diff(age) > step * 1.01))]
    
    X[, .(
      contrast = flag,
      start_age = min(age),
      end_age   = max(age)
    ), by = grp]
  }
  
  sig_ranges <- rbindlist(
    lapply(names(contrasts), get_ranges, DT = W, age_col = age_col),
    fill = TRUE
  )
  
  # ---- Plot ----
  p <- ggplot(df, aes(
    x = !!sym(age_col),
    y = !!sym(y_col),
    colour = !!sym(group_col),
    linetype = !!sym(group_col),
    group = !!sym(group_col)
  )) +
    geom_point(alpha = 0.3) +
    geom_line(data = mod_df, aes(y = fit), linewidth = 1) +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_colour_manual(values = colour_values) +
    scale_linetype_manual(values = linetype_values)
  
  if (use_cohort) {
    p <- p + aes(shape = !!sym(cohort))
  }
  
  if (!is.null(sig_ranges) && nrow(sig_ranges) > 0) {
    sig_ranges[, y := Inf - as.numeric(factor(contrast)) * 0.03]
    
    p <- p +
      geom_segment(
        data = sig_ranges,
        aes(x = start_age, xend = end_age,
            y = y, yend = y,
            colour = contrast),
        inherit.aes = FALSE,
        linewidth = 2
      )
  }
  
  if(length(YLIM) > 1){
    p <- p + ylim(YLIM)
  }
  return(list(
    plot = p,
    significance_ranges = sig_ranges,
    stats = W
  ))
}
