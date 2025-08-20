setwd("~/Imperial/Thesis/Data")

library(readr)
library(dplyr)
library(metafor)
library(gt)

data <- read_csv("ponds_data_extraction.csv")

# Filter out rows with missing essential data
data_c <- data %>%
  filter(!is.na(i_bd_mean) & !is.na(c_bd_mean) &
           !is.na(i_sample_size) & !is.na(c_sample_size))

# Convert percentages to number of ponds
data_c <- data_c %>%
  mutate(
    i_bd_mean = ifelse(
      bd_metric_unit == "percentage of ponds (%)",
      (i_bd_mean / 100) * i_sample_size,
      i_bd_mean
    ),
    c_bd_mean = ifelse(
      bd_metric_unit == "percentage of ponds (%)",
      (c_bd_mean / 100) * c_sample_size,
      c_bd_mean
    ),
    bd_metric_unit = ifelse(
      bd_metric_unit == "percentage of ponds (%)",
      "number of ponds",
      bd_metric_unit
    )
  )

# Create counts of ponds with and without amphibians
data_c <- data_c %>%
  mutate(
    i_with = round(i_bd_mean),  
    i_without = i_sample_size - i_with,
    c_with = round(c_bd_mean),  
    c_without = c_sample_size - c_with
  )

# Calculate log odds ratio using escalc
data_clean <- escalc(
  measure = "OR",    
  ai = i_with,       
  bi = i_without,    
  ci = c_with,       
  di = c_without,    
  data = data_c
)

# Multilevel meta-analysis
meta_all <- rma.mv(yi, vi, random = ~ 1 | study_ID, data = data_clean, method = "REML")
print(meta_all)

# Study-level summary for forest plot
study_summary <- data_clean %>%
  group_by(label) %>%
  summarise(
    n_effects = n(),
    mean_lor = mean(yi, na.rm = TRUE),
    pooled_var = sum(vi, na.rm = TRUE) / n()^2,
    .groups = 'drop'
  ) %>%
  mutate(
    pooled_se = sqrt(pooled_var),
    ci_lower = mean_lor - 1.96 * pooled_se,
    ci_upper = mean_lor + 1.96 * pooled_se,
    odds_ratio = exp(mean_lor),
    ci_lower_or = exp(ci_lower),
    ci_upper_or = exp(ci_upper)
  )

# Meta-analysis for study-level data
meta_combined <- rma(yi = mean_lor, vi = pooled_var, data = study_summary, method = "REML")
print(meta_combined)

# Create forest plot for all effect sizes - metafor
metafor::forest(
  meta_all,
  slab = data_clean$label,
  xlim = c(-15, 20),
  cex = 0.75,
  header = TRUE,
  mlab = "Overall Effect Size",
  xlab = "Log Response Ratio (± 95% CI)"
)

# Create table for forest plot
plot_data <- study_summary %>%
  mutate(
    Study = paste(label),
    Estimate = round(mean_lor, 2),   # now using mean_lor instead of mean_lrr
    LC = round(ci_lower, 2),
    UC = round(ci_upper, 2),
    NES = n_effects
  ) %>%
  select(Study, Estimate, LC, UC, NES) %>%
  rename(mean = Estimate, lower = LC, upper = UC)

# Summary 
summary_row <- data.frame(
  Study = "Overall Effect Size",
  mean = as.numeric(meta_combined$b),         
  lower = as.numeric(meta_combined$ci.lb),
  upper = as.numeric(meta_combined$ci.ub),
  NES = sum(study_summary$n_effects)
)

# Combine study-level and overall rows
forest_data <- rbind(plot_data, summary_row)

# Format table text for forest plot
table_text <- cbind(
  Study = forest_data$Study,
  `Confidence Interval` = "",
  Estimate = sprintf("%.2f", forest_data$mean),
  `CI Lower` = sprintf("%.2f", forest_data$lower),
  `CI Upper` = sprintf("%.2f", forest_data$upper),
  `#ES` = forest_data$NES
)

# Forest plot (log odds ratio)
forest_plot <- forestploter::forest(
  table_text,
  est = forest_data$mean,
  lower = forest_data$lower,
  upper = forest_data$upper,
  ci_column = 2,
  xlim = c(-3, 6),                       
  ticks_at = c(-2, -1, 0, 1, 2, 3, 4, 5, 6),
  ref_line = 0,                          
  sizes = 0.5,
  col_widths = c(2.5, 100, 2, 1.2, 1.2, 1.2),
  new_page = TRUE
)
print(forest_plot)

# Funnel plots
# For all individual effect sizes
funnel_all <- funnel(meta_all,
                     xlab = "Log Odds Ratio",
                     ylab = "Standard Error",
                     main = "Funnel Plot - Publication Bias Assessment")

# For study-level data
funnel_comb <- funnel(meta_combined,
                      xlab = "Log Odds Ratio",
                      ylab = "Standard Error",
                      main = "Funnel Plot - Publication Bias Assessment")

# Publication bias tests
# convert column names
data_c1 <- data_clean %>%
  rename(
    mean_lrr = yi,    # effect size column
    lrr_var = vi    # variance column
  )

# simple meta-analysis for eggers test
meta_simple <- rma(yi = mean_lrr,
                   vi = lrr_var,
                   data = data_c1,
                   method = "REML")

# Run Egger's regression test
egger_test <- regtest(meta_simple, model = "lm")  
egger_test

# Create summary table 
# Create categories 
summary_table <- data %>%
  mutate(
    # Taxon categories
    taxon_group = case_when(
      grepl("Newt|Triturus|Lissotriton|Ichthyosaura", taxon_scientific, ignore.case = TRUE) ~ "Newts",
      grepl("Toad|Bufo|Bombina|Pelobates|Pseudepidalea", taxon_scientific, ignore.case = TRUE) ~ "Toads",
      grepl("Frog|Rana|Pelophylax|Hyla", taxon_scientific, ignore.case = TRUE) ~ "Frogs",
      grepl("Salamander|Salamandra", taxon_scientific, ignore.case = TRUE) ~ "Salamanders",
      TRUE ~ "Other"
    ),
    
    # Biodiversity metric
    metric_type = case_when(
      bd_metric == "occurrence" | bd_metric == "occurrence frequency" ~ "Occurrence",
      grepl("abundance", bd_metric, ignore.case = TRUE) ~ "Abundance",
      TRUE ~ "Other"
    ),
    
    # Climate categories
    climate_cat = case_when(
      grepl("temperate", climate, ignore.case = TRUE) ~ "Temperate",
      grepl("continental", climate, ignore.case = TRUE) ~ "Continental",
      TRUE ~ "Other"
    ),
    
    # Scale categories
    scale_cat = case_when(
      scale_of_intervention == "very small (1)" ~ "Very small-scale",
      scale_of_intervention == "small (2)" ~ "Small-scale",
      scale_of_intervention == "medium (3)"~ "Medium-scale",
      scale_of_intervention == "large (4)" ~ "Large-scale",
      TRUE ~ "Other"
    )
  )

# Create summary statistics table
create_summary_table <- function(data, grouping_var, group_name) {
  data %>%
    group_by(!!sym(grouping_var)) %>%
    summarise(
      moderator_level = first(!!sym(grouping_var)),
      number_of_effect_sizes = n(),
      .groups = 'drop'
    ) %>%
    filter(!is.na(moderator_level) & moderator_level != "Other" & moderator_level != "Unknown") %>%
    select(moderator_level, number_of_effect_sizes) %>%
    mutate(category = group_name)

}

# Create individual summary tables
taxon_summary <- create_summary_table(summary_table, "taxon_other", "Taxon")
metric_summary <- create_summary_table(summary_table, "metric_type", "Metric")
climate_summary <- create_summary_table(summary_table, "climate_cat", "Climate")
scale_summary <- create_summary_table(summary_table, "scale_cat", "Scale")

# Combine all summaries
final_summary_table <- bind_rows(
  taxon_summary,
  metric_summary,
  climate_summary,
  scale_summary
) %>%
  select(moderator_level, number_of_effect_sizes)

final_summary_table <- final_summary_table %>%
  rename(
    `Moderator level` = moderator_level,
    `Number of effect sizes` = number_of_effect_sizes
  )

print(final_summary_table, row.names = FALSE)
gt(final_summary_table)

final_summary_table %>%
  gt() %>%
  tab_header(
    title = "Summary of Effect Sizes by Moderator Level"
  ) %>%
  gtsave("summary_table.png")