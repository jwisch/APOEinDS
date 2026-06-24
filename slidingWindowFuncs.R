# Function for sliding window proportions
sliding_proportion <- function(df, age_col = "Age", group_col = "APOE_grouped",
                               outcome_col = "Apos", window_size = 5) {
  
  # Get age range
  min_age <- floor(min(df[[age_col]], na.rm = TRUE))
  max_age <- ceiling(max(df[[age_col]], na.rm = TRUE))
  
  # Define window centers (you could also use seq(min_age, max_age, by=1))
  window_centers <- seq(min_age, max_age, by = 1)
  
  # Loop over windows
  result <- lapply(window_centers, function(center) {
    lower <- center - window_size/2
    upper <- center + window_size/2
    
    df %>%
      filter(.data[[age_col]] >= lower & .data[[age_col]] < upper) %>%
      group_by(.data[[group_col]]) %>%
      summarise(
        n = n(),
        n_pos = sum(.data[[outcome_col]] == 1, na.rm = TRUE),
        prop_pos = ifelse(n > 0, n_pos / n, NA_real_),
        .groups = "drop"
      ) %>%
      mutate(window_center = center,
             window_lower = lower,
             window_upper = upper)
  }) %>%
    bind_rows()
  
  return(result)
}


sliding_proportion_bs <- function(df,
                                  age_col = "Age",
                                  group_col = "APOE_grouped",
                                  outcome_col = "Apos",
                                  window_size = 5,
                                  n_boot = 1000,
                                  ci = c(0.025, 0.975)) {
  # Convert to data.table
  setDT(df)
  
  # Age range
  min_age <- floor(min(df[[age_col]], na.rm = TRUE))
  max_age <- ceiling(max(df[[age_col]], na.rm = TRUE))
  window_centers <- seq(min_age, max_age, by = 1)
  
  res_list <- vector("list", length(window_centers))
  
  for (i in seq_along(window_centers)) {
    center <- window_centers[i]
    lower <- center - window_size/2
    upper <- center + window_size/2
    
    df_window <- df[get(age_col) >= lower & get(age_col) < upper]
    if (nrow(df_window) == 0L) next
    
    # Observed summary
    observed <- df_window[,
                          .(n = .N,
                            n_pos = sum(get(outcome_col) == 1, na.rm = TRUE)),
                          by = group_col
    ][, prop_pos := fifelse(n > 0, n_pos / n, NA_real_)]
    
    # Bootstrap
    boot_res <- df_window[,
                          {
                            vals <- get(outcome_col)
                            n <- .N
                            if (n == 0) {
                              list(ci_lower = NA_real_, ci_upper = NA_real_)
                            } else {
                              boot_props <- replicate(
                                n_boot,
                                {
                                  samp <- sample(vals, size = ceiling(0.8 * n), replace = TRUE)
                                  mean(samp == 1, na.rm = TRUE)
                                }
                              )
                              list(ci_lower = quantile(boot_props, probs = ci[1], na.rm = TRUE),
                                   ci_upper = quantile(boot_props, probs = ci[2], na.rm = TRUE))
                            }
                          },
                          by = group_col
    ]
    
    out <- merge(observed, boot_res, by = group_col, all.x = TRUE)
    out[, `:=`(window_center = center,
               window_lower = lower,
               window_upper = upper)]
    
    res_list[[i]] <- out
  }
  
  result <- rbindlist(res_list, use.names = TRUE, fill = TRUE)
  return(result)
}





get_bs_plot <- function(DF){
  # Base plot
  # Histogram on top
  hist_top <- ggplot(DF, aes(x = window_center, y = n, fill = APOE_grouped)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),  # remove x-axis title
      axis.text.x = element_blank(),   # remove x-axis labels
      axis.ticks.x = element_blank(),
      legend.position = "none"         # hide legend to avoid duplication
    ) +
    labs(y = "Count") + xlim(c(25, 55))     + scale_fill_manual(values =  c("#F8766D",   # blue
                                                                            "#00BA38",   # green
                                                                            "#619CFF"))# red
  
  # Main plot with smoothed lines and ribbon
  main_plot <- ggplot(DF, aes(x = window_center, y = prop_pos, group = APOE_grouped, colour = APOE_grouped)) +
    geom_smooth(method = "gam", se = FALSE) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = APOE_grouped), alpha = 0.2, colour = NA) +
    theme_bw() +
    theme(legend.position = "bottom") +
    xlim(c(25, 55)) +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    labs(y = "Proportion Positive", x = "Age") +
    scale_fill_manual(values =  c("#F8766D",   # blue
                                  "#00BA38",   # green
                                  "#619CFF")# red
    ) +
    scale_colour_manual(values =  c("#F8766D",   # blue
                                    "#00BA38",   # green
                                    "#619CFF"))# red
  
  # Stack histogram above main plot
  p <- hist_top / main_plot + plot_layout(heights = c(1, 3))
  return(p)}


# General interpolation helper
interp_cross <- function(x, y, cutoff = 0.5) {
  idx <- which(y >= cutoff)[1]
  if (is.na(idx) || idx == 1) return(NA_real_)
  
  x0 <- x[idx - 1]; x1 <- x[idx]
  y0 <- y[idx - 1]; y1 <- y[idx]
  
  x0 + (cutoff - y0) * (x1 - x0) / (y1 - y0)
}

# # Apply to prop_pos, ci_lower, ci_upper
# estimate_age50_ci <- function(df, cutoff = 0.5) {
#   df <- arrange(df, window_center)
#   
#   age50       <- interp_cross(df$window_center, df$prop_pos, cutoff)
#   age50_lower <- interp_cross(df$window_center, df$ci_upper, cutoff) # earliest
#   age50_upper <- interp_cross(df$window_center, df$ci_lower, cutoff) # latest
#   
#   tibble(
#     age50 = age50,
#     age50_lower = age50_lower,
#     age50_upper = age50_upper
#   )
# }


# Helper to interpolate y at a given x
interp_y <- function(x, y, x_target) {
  idx <- which(x >= x_target)[1]
  if (is.na(idx) || idx == 1) return(NA_real_)
  
  x0 <- x[idx - 1]; x1 <- x[idx]
  y0 <- y[idx - 1]; y1 <- y[idx]
  
  y0 + (x_target - x0) * (y1 - y0) / (x1 - x0)
}



# Step 2: Get prop_pos for APOE2 and APOE3 at that age
get_prop_at_age <- function(df, age) {
  df <- arrange(df, window_center)
  interp_y(df$window_center, df$prop_pos, age)
}


get_hist_of_counts <- function(DF, YLIM = NULL){
  
  p <- ggplot(DF, aes(x = window_center, y = n, fill = Group)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(y = "Count") +
    xlim(c(25, 55)) +
    scale_fill_manual(values = c(
      "#FDAA9F",
      "#66D98E",
      "#A6C9FF",
      "#F8766D",
      "#00BA38",
      "#619CFF"
    ))
  
  # optionally add y-axis limits
  if (!is.null(YLIM)) {
    p <- p + ylim(YLIM)
  }
  
  return(p)
}
get_bs_plot_sex <- function(DF, YLIM = NULL){
  # Base plot
  # Histogram on top
  hist_top <- get_hist_of_counts(DF, YLIM)  + theme(      axis.title.x = element_blank(),  # remove x-axis title
                                                          axis.text.x = element_blank(),   # remove x-axis labels
                                                          axis.ticks.x = element_blank(),
                                                          legend.position = "none" )
  
  
  # Main plot with smoothed lines and ribbon
  main_plot <- ggplot(DF, aes(x = window_center, y = prop_pos, group = Group, colour = Group, linetype = Group)) +
    geom_smooth(method = "gam", se = FALSE) +
    # geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = Group), alpha = 0.2, colour = NA) +
    theme_bw() +
    theme(legend.position = "bottom") +
    xlim(c(25, 55)) +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    labs(y = "Proportion Positive", x = "Age") +
    scale_fill_manual(values =  c("#FDAA9F",   # blue
                                  "#66D98E",   # green
                                  "#A6C9FF",
                                  "#F8766D",   # blue
                                  "#00BA38",   # green
                                  "#619CFF")# red
    ) +
    scale_colour_manual(values =  c("#FDAA9F",   # blue
                                    "#66D98E",   # green
                                    "#A6C9FF",
                                    "#F8766D",   # blue
                                    "#00BA38",   # green
                                    "#619CFF")) +# red
    scale_linetype_manual(values = c("solid", "solid", "solid", "dashed", "dashed", "dashed"))
  
  # Stack histogram above main plot
  p <- hist_top / main_plot + plot_layout(heights = c(1, 3))
  return(p)}

safe_interp_cross <- function(x, y, cutoff = 0.5) {
  # coerce to numeric
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  # remove NA pairs
  ok <- !(is.na(x) | is.na(y))
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n == 0) return(NA_real_)
  
  # --- CHANGE MADE HERE: use latest x instead of earliest ---
  eq_idx <- which(y == cutoff)
  if (length(eq_idx) >= 1) {
    return(min(x[eq_idx], na.rm = TRUE))   # <–– updated
  }
  
  # Look for bracketing consecutive points
  d <- y - cutoff
  if (length(d) < 2) return(NA_real_)
  
  # find first sign-change index
  sign_change_idx <- which(d[-length(d)] * d[-1] < 0)
  if (length(sign_change_idx) == 0) {
    return(NA_real_)
  }
  
  i <- sign_change_idx[1]
  x1 <- x[i]; x2 <- x[i+1]
  y1 <- y[i]; y2 <- y[i+1]
  
  # linear interpolation
  if (y2 == y1) return(mean(c(x1, x2)))
  x_at <- x1 + (cutoff - y1) * (x2 - x1) / (y2 - y1)
  return(as.numeric(x_at))
}


estimate_age50_ci <- function(df, cutoff = 0.5) {
  df <- arrange(df, window_center)
  
  age50       <- safe_interp_cross(df$window_center, df$prop_pos, cutoff)
  age50_lower <- safe_interp_cross(df$window_center, df$ci_upper, cutoff) # earliest
  age50_upper <- safe_interp_cross(df$window_center, df$ci_lower, cutoff) # latest
  
  tibble(
    age50 = age50,
    age50_lower = age50_lower,
    age50_upper = age50_upper
  )
}
