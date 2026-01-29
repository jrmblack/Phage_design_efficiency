library(pROC)
library(broom)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(ggtree)
library(treeio)
library(ape)
library(tidytree)
library(patchwork)
library(Biostrings)
library(cowplot)
library(scales)
library(extrafont)
#font_import()  # This takes a while
loadfonts()

library(extrafont)
loadfonts(device = "pdf")

df_final <- readr::read_tsv("~/Downloads/df_final.tsv")

theme_pub <- theme_bw(base_size = 10, base_family = "Arial") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#444444")
    strip.text = element_text(color = "white", face = "bold", family = "Arial"),
    legend.text = element_text(size = 9, family = "Arial"),
    legend.title = element_text(face = "bold", family = "Arial"),
    axis.text = element_text(color = "black", family = "Arial"),
    axis.title = element_text(face = "bold", family = "Arial"),
    legend.background = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_blank()
  )

numeric_cols_force <- c(
  "evo2_7b_likelihood", "evo1_131k_base_likelihood", "king_et_al_score", 
  "NT250_pseudollh", # Added NT
  "median_patristic_dist", "median_blast_identity",
  "reference_genome_percent_identity", "average_protein_percent_identity",
  "tropism_protein_mmseqs_percent_identity", "architecture_similarity_score",
  "num_syntenic_genes", "gc_content", "max_nt_homopolymer_length", 
  "total_num_genes", "genome_length", 
  "prodigal_orf_count", "prodigal_coding_density", 
  "gibson_overlap2_tm"
)

analysis_df <- df_final %>%
  mutate(validated = as.logical(validated)) %>%
  mutate(across(any_of(numeric_cols_force), as.numeric)) %>%
  filter(!is.na(evo2_7b_likelihood)) %>%
  
  # Transformations
  mutate(
    gc_dist_from_optimum = abs(gc_content - 44.8),
    len_dist_from_optimum = abs(genome_length - 5386),
    median_blast_distance = 100 - median_blast_identity,
    inv_patristic_dist  = -median_patristic_dist,
    inv_blast_dist      = -median_blast_distance,
    inv_gc_dist         = -gc_dist_from_optimum,
    inv_len_dist        = -len_dist_from_optimum,
    inv_homopolymer     = -max_nt_homopolymer_length
  )

calculate_precision_curve <- function(data, score_col, label_base) {
  sub_data <- data %>% filter(!is.na(.data[[score_col]]))
  auc_val <- roc(sub_data$validated, sub_data[[score_col]], quiet=TRUE)$auc
  label_full <- paste0(label_base, " [AUC: ", sprintf("%.2f", auc_val), "]")
  
  sub_data %>%
    arrange(desc(.data[[score_col]])) %>%
    mutate(
      rank = row_number(),
      is_validated = ifelse(validated, 1, 0),
      cumulative_hits = cumsum(is_validated),
      success_rate = cumulative_hits / rank,
      model_label = label_full
    ) %>%
    select(rank, success_rate, model_label)
}

curve_evo2   <- calculate_precision_curve(analysis_df, "evo2_7b_likelihood", "Evo2 (7b)")
curve_evo1   <- calculate_precision_curve(analysis_df, "evo1_131k_base_likelihood", "Evo1 (131k)")
curve_phage  <- calculate_precision_curve(analysis_df, "king_et_al_score", "Evo2 (Bacteriophage)")
curve_nt     <- calculate_precision_curve(analysis_df, "NT250_pseudollh", "NT-250m (Baseline)")

all_curves <- bind_rows(curve_evo2, curve_evo1, curve_phage, curve_nt)
baseline_rate <- sum(analysis_df$validated) / nrow(analysis_df)

# --- Colors ---
labels <- unique(all_curves$model_label)
color_map_base <- c("Evo2 (7b)"="#FC4E07", "Evo1 (131k)"="#2E9FDF", 
                    "Evo2 (Bacteriophage)"="#E7B800", "NT-250m (Baseline)"="purple")
final_colors <- setNames(character(length(labels)), labels)
for(lbl in labels) { for(base in names(color_map_base)) { if(grepl(base, lbl, fixed=TRUE)) final_colors[lbl] <- color_map_base[base] } }

# --- Plot 1c ---
p1c_arial <- ggplot(all_curves, aes(x = rank, y = success_rate, color = model_label)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = baseline_rate, linetype = "dashed", color = "black") +

  scale_color_manual(values = final_colors) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +

  labs(x = NULL, y = NULL, color = NULL) + 
  
  theme_pub +
  theme(legend.position = "right")

print(p1c_arial)


# Plot 1c
p1c_integrated <- ggplot(all_curves, aes(x = rank, y = success_rate, color = model_label)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = baseline_rate, linetype = "dashed") +
  scale_color_manual(values = final_colors) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  theme_pub + theme(legend.position = "right") +
  labs(x = "Ranked Sequences", y = "Cumulative Viability Rate")

print(p1c_integrated)
# ==============================================================================
# 4. SAVE
# ==============================================================================
ggsave("~/Documents/PhD/Projects/Biosecurity/bacteriophage_paper/figure_1c_topN_clean.png", p1c_clean, width = 6, height = 3,dpi = 600, device = "png")

###################################################
###################################################
###################################################
########################   now 1d       ###########
###################################################
###################################################
###################################################
###################################################
work_dir <- "~/Downloads"
tree_file <- file.path(work_dir, "combined_analysis_aligned.fasta.contree")
id_file   <- file.path(work_dir, "natural_ids.txt")

tree <- read.tree(tree_file)
try({ library(phytools); tree <- midpoint.root(tree) }, silent=TRUE)

natural_ids <- readLines(id_file)
validated_ids <- df_final$prompt_id[df_final$validated == TRUE]

# Define Greek Label
natural_label <- expression(paste("Natural ", Phi, "X174 Genomes"))

tree_data <- data.frame(label = tree$tip.label) %>%
  mutate(
    type = ifelse(label %in% natural_ids, "Natural Reference", "Synthetic Genome"),
    is_validated = label %in% validated_ids,
    
    display_group = case_when(
      type == "Natural Reference" ~ "Natural",
      is_validated ~ "Viable",
      TRUE ~ "Non-Viable"
    ),
    
    display_group = factor(display_group, levels = c("Non-Viable", "Natural", "Viable"))
  ) %>%
  arrange(display_group)

# ==============================================================================
# 2. GENERATE LEGEND (Vertical, Ordered, Smaller Dots)
# ==============================================================================
p_legend_source <- ggtree(tree, layout = "rectangular") %<+% tree_data +
  geom_tippoint(aes(color = display_group, alpha = display_group), size = 2.5) + 
  
  scale_color_manual(
    values = c("Natural" = "black", "Viable" = "#FC4E07", "Non-Viable" = "grey70"),
    breaks = c("Natural", "Viable", "Non-Viable"),
    labels = c("Natural" = natural_label, 
               "Viable" = "Synthetic (Viable)", 
               "Non-Viable" = "Synthetic (Non-Viable)")
  ) +
  
  scale_alpha_manual(
    values = c("Natural" = 1.0, "Viable" = 1.0, "Non-Viable" = 0.6),
    breaks = c("Natural", "Viable", "Non-Viable"),
    labels = c("Natural" = natural_label, 
               "Viable" = "Synthetic (Viable)", 
               "Non-Viable" = "Synthetic (Non-Viable)")
  ) +
  
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 10, family = "sans"), # 'sans' fixes the font error
    legend.key = element_blank(),
    legend.spacing.y = unit(0.5, "cm")
  ) +
  
  guides(
    color = guide_legend(ncol = 1, override.aes = list(size = 2.5)),
    alpha = guide_legend(ncol = 1, override.aes = list(size = 2.5))
  )

# Extract the legend
legend_vertical <- get_legend(p_legend_source)

# ==============================================================================
# 3. GENERATE TREE 
# ==============================================================================
p_tree_rect <- ggtree(tree, layout = "rectangular", size = 0.3) %<+% tree_data +
  
  geom_tippoint(aes(color = display_group, alpha = display_group), size = 2.5) +
  
  scale_color_manual(values = c("Natural" = "black", "Viable" = "#FC4E07", "Non-Viable" = "grey70")) +
  scale_alpha_manual(values = c("Natural" = 1.0, "Viable" = 1.0, "Non-Viable" = 0.6)) +
  coord_cartesian(clip = "off") + 
  geom_treescale(x = 0, y = 5, width = 0.1, fontsize = 3) +
  
  theme(legend.position = "none",
        plot.margin = margin(10, 10, 10, 10))

# ==============================================================================
# 4. SAVE 
# ==============================================================================
# Save Tree
ggsave("~/Documents/PhD/Projects/Biosecurity/bacteriophage_paper/figure_1d_tree_clean.png", 
       p_tree_rect, width = 3.8, height = 6.5, dpi = 600, device = "png")

# Save Legend
ggsave("~/Documents/PhD/Projects/Biosecurity/bacteriophage_paper/figure_1d_legend.png", 
       legend_vertical, width = 3, height = 2, dpi = 600, device = "png")
###################################################
###################################################
###################################################
############next up 1e ############################
###################################################
###################################################
###################################################

plot_data <- df_final %>%
  mutate(across(c(median_patristic_dist, evo2_7b_likelihood), as.numeric)) %>%
  filter(!is.na(evo2_7b_likelihood), !is.na(median_patristic_dist))

# ==============================================================================
# 2. PLOT
# ==============================================================================
p_scatter_clean <- ggplot(plot_data, aes(x = median_patristic_dist, y = evo2_7b_likelihood)) +
  
  geom_point(color = "grey30", alpha = 0.6, size = 2.5) +
  
  geom_smooth(method = "lm", formula = y ~ x, color = "black", 
              linetype = "dashed", size = 0.8, se = TRUE, fill = "grey85") +
  labs(x = NULL,
       y = NULL) +
  theme_bw(base_size = 10, base_family = "sans") +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black", size = 10)
  )

print(p_scatter_clean)
# ==============================================================================
# 3. SAVE
# ==============================================================================

ggsave("~/Documents/PhD/Projects/Biosecurity/bacteriophage_paper/figure_1e_scatter_clean.png", 
       p_scatter_clean, width = 3.5, height = 3, dpi = 600, device = "png")

###################################################
###################################################
############next up 1f ############################
###################################################
###################################################
###################################################

numeric_cols_force <- c(
  "evo2_7b_likelihood", "evo1_131k_base_likelihood", "king_et_al_score", 
  "NT250_pseudollh", "evo2_1b_likelihood",
  "median_patristic_dist", "median_blast_identity",
  "reference_genome_percent_identity", "average_protein_percent_identity",
  "tropism_protein_mmseqs_percent_identity", "architecture_similarity_score",
  "num_syntenic_genes", "gc_content", "max_nt_homopolymer_length", 
  "total_num_genes", "genome_length", 
  "prodigal_orf_count", "prodigal_coding_density", 
  "gibson_overlap2_tm"
)

analysis_df <- df_final %>%
  mutate(validated = as.logical(validated)) %>%
  mutate(across(any_of(numeric_cols_force), as.numeric)) %>%
  filter(!is.na(evo2_7b_likelihood)) %>%
  
  # Transformations
  mutate(
    gc_dist_from_optimum = abs(gc_content - 44.8),
    len_dist_from_optimum = abs(genome_length - 5386),
    median_blast_distance = 100 - median_blast_identity,
    
    # Inversions (Positive = Success)
    inv_patristic_dist  = -median_patristic_dist,
    inv_blast_dist      = -median_blast_distance,
    inv_gc_dist         = -gc_dist_from_optimum,
    inv_len_dist        = -len_dist_from_optimum,
    inv_homopolymer     = -max_nt_homopolymer_length
  )

# ==============================================================================
# 2. DEFINE GROUPS
# ==============================================================================
vars_model <- c("evo2_7b_likelihood") 

vars_bio <- c(
  "inv_gc_dist", "inv_homopolymer", "inv_len_dist",
  "total_num_genes", "prodigal_orf_count", "prodigal_coding_density",
  "gibson_overlap2_tm"
)

vars_evo <- c(
  "inv_patristic_dist", "inv_blast_dist",
  "reference_genome_percent_identity", "average_protein_percent_identity",
  "tropism_protein_mmseqs_percent_identity",
  "architecture_similarity_score", "num_syntenic_genes"
)

all_vars_list <- list(
  " " = vars_model,
  "Genome Properties" = vars_bio,
  "Homology/Evolutionary Similarity" = vars_evo
)

# ==============================================================================
# 3. RUN MODELS
# ==============================================================================
results_list <- list()

for (group_name in names(all_vars_list)) {
  vars <- all_vars_list[[group_name]]
  vars <- vars[vars %in% names(analysis_df)]
  
  for (var in vars) {
    temp_data <- analysis_df %>%
      select(validated, x_var = all_of(var)) %>%
      drop_na() %>% filter(var(x_var) > 0) %>% 
      mutate(x_scaled = scale(x_var))
    
    if(nrow(temp_data) < 10) next
    
    try({
      model <- glm(validated ~ x_scaled, data = temp_data, family = "binomial")
      res <- tidy(model, conf.int = TRUE) %>%
        filter(term == "x_scaled") %>%
        mutate(variable = var, group = group_name)
      results_list[[length(results_list) + 1]] <- res
    }, silent = TRUE)
  }
}

# ==============================================================================
# 4. PROCESS LABELS
# ==============================================================================
all_res <- bind_rows(results_list) %>%
  mutate(
    # Simple Binary Significance
    is_significant = p.value < 0.05,
    
    label = recode(variable,
                   "evo2_7b_likelihood" = "Evo2 (7b) Likelihood",
                   "inv_patristic_dist" = "Patristic Dist. (-)",
                   "inv_blast_dist" = "BLAST Dist. (-)",
                   "reference_genome_percent_identity" = "Ref. Genome Identity (%)",
                   "average_protein_percent_identity" = "Avg. Protein Identity (%)",
                   "tropism_protein_mmseqs_percent_identity" = "Spike Protein Identity (%)",
                   "architecture_similarity_score" = "Architecture Score",
                   "num_syntenic_genes" = "# Syntenic Genes",
                   "inv_gc_dist" = "GC Dist. (from WT)",
                   "inv_homopolymer" = "Max Homopolymer",
                   "inv_len_dist" = "Length Dist. (from WT)",
                   "total_num_genes" = "Total # Genes",
                   "prodigal_orf_count" = "ORF Count",
                   "prodigal_coding_density" = "Coding Density",
                   "gibson_overlap2_tm" = "Gibson Overlap Tm"
    )
  ) %>%
  mutate(group = factor(group, levels = names(all_vars_list))) %>%
  group_by(group) %>% arrange(estimate) %>% ungroup() %>%
  mutate(label_ordered = factor(label, levels = unique(label)))

# ==============================================================================
# 5. PLOT
# ==============================================================================
p_forest_clean <- ggplot(all_res, aes(x = estimate, y = label_ordered)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
   geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = is_significant), height = 0.4) +
   geom_point(aes(color = is_significant), size = 3) +
  facet_grid(group ~ ., scales = "free_y", space = "free") +
  scale_color_manual(values = c("TRUE" = "#2E9FDF", "FALSE" = "grey70")) +
  labs(x = NULL, y = NULL) +
   theme_bw(base_family = "Arial") +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 9, color = "white"),
    strip.background = element_rect(fill = "#444444"),
    axis.text = element_text(size = 10, color = "black"),
    panel.grid.major.y = element_blank()
  )

print(p_forest_clean)

ggsave("~/Documents/PhD/Projects/Biosecurity/bacteriophage_paper/figure_1f_univ_clean.png", 
       p_forest_clean, width = 5.5, height = 7,dpi = 600, device = "png")

#1g
baseline_vars <- c(
  "num_syntenic_genes", "average_protein_percent_identity", "architecture_similarity_score", 
  "reference_genome_percent_identity", "median_patristic_dist", 
  "gibson_overlap2_tm", "gc_dist_from_optimum"
)

model_data <- analysis_df %>%
  select(validated, evo2_7b_likelihood, all_of(baseline_vars)) %>%
  drop_na() %>%
  mutate(across(-validated, scale))

model_baseline <- glm(validated ~ . -evo2_7b_likelihood, data = model_data, family = "binomial")
model_combined <- glm(validated ~ ., data = model_data, family = "binomial")

calc_mcfadden <- function(model) { 1 - (logLik(model) / logLik(update(model, . ~ 1))) }
r2_base <- as.numeric(calc_mcfadden(model_baseline))
r2_comb <- as.numeric(calc_mcfadden(model_combined))

# ==============================================================================
# 2. PLOT
# ==============================================================================
lift_data <- data.frame(
  Model = factor(c("Baseline", "Combined"), levels = c("Baseline", "Combined")),
  R2 = c(r2_base, r2_comb)
)

p_lift_clean <- ggplot(lift_data, aes(x = Model, y = R2, fill = Model)) +
  
  # Bars
  geom_col(width = 0.6, color = "black", alpha = 0.9) +
  scale_fill_manual(values = c("grey70", "#FC4E07")) +
  scale_y_continuous(
    breaks = c(0, 0.2, 0.4),           
    labels = c("0", "20", "40"),      
    limits = c(0, 0.45),               
    expand = c(0, 0)                  
  ) +
   labs(x = NULL, y = NULL) +
  theme_bw(base_family = "sans") +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 10, color = "black", face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
     panel.border = element_rect(color = "black", fill = NA)
  )

print(p_lift_clean)
ggsave("~/Documents/PhD/Projects/Biosecurity/bacteriophage_paper/figure_1g_lift.png", p_lift_clean, width = 4, height = 4, device = "png")


# TOY TREE
col_base <- "#4A90A4"    # Blue - base/original lineages
col_new  <- "#D64550"    # Red - new/additional taxa
# Structure: new taxa (n1-n6) are inserted as sisters to various base taxa (b1-b8)
tree1 <- read.tree(text = paste0(
    "(((b1:2,(b2:1,n1:1):1):1,(b3:2,n2:2):1):1,",
    "((b4:2,(n3:1,b5:1):1):1,((b6:1,n4:1):1,(n5:1,(b7:1,(b8:1,n6:1):0.5):0.5):1):1):1);"
))
  

tree1_groups <- data.frame(
  label = tree1$tip.label,
  group = ifelse(grepl("^b", tree1$tip.label), "base", "new")
)

tree2 <- read.tree(text = paste0(
  "(((b1:2,b2:2):1,(b3:2,b4:2):1):1,",
  "((b5:2,(b6:1,(b7:1,b8:1):0.5):1):1,",
  "((n1:1,n2:1):1,(n3:1,(n4:1,(n5:1,n6:1):0.5):0.5):1):1):1);"
))

tree2_groups <- data.frame(
  label = tree2$tip.label,
  group = ifelse(grepl("^b", tree2$tip.label), "base", "new")
)

minimal_theme <- theme_void() +
  theme(
    plot.margin = margin(30, 30, 30, 30)
  )

# -----------------------------------------------------------------------------
# Plot 1: Intermingled
# -----------------------------------------------------------------------------

p1 <- ggtree(tree1, size = 1.1, colour = "grey40") %<+% tree1_groups +
  geom_tippoint(aes(colour = group), size = 6) +
  scale_colour_manual(values = c("base" = col_base, "new" = col_new)) +
  minimal_theme +
  theme(legend.position = "none")
  
# -----------------------------------------------------------------------------
# Plot 2: Clustered / New branches
# -----------------------------------------------------------------------------

p2 <- ggtree(tree2, size = 1.1, colour = "grey40") %<+% tree2_groups +
  geom_tippoint(aes(colour = group), size = 6) +
  scale_colour_manual(values = c("base" = col_base, "new" = col_new)) +
  minimal_theme +
  theme(legend.position = "none")

# -----------------------------------------------------------------------------
# Combine - stacked vertically
# -----------------------------------------------------------------------------

final_plot <- p1 / p2

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------

ggsave("~/Documents/PhD/Projects/Biosecurity/bacteriophage_paper/evo_trees_comparison.png", final_plot, 
       width = 4, height = 7, dpi = 600, device = "png", bg = "white")
 
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
############################# fig 2 now ##################
##########################################################

out_path <- "~/Documents/Phd/Projects/Biosecurity/bacteriophage_paper/"

# Define Standard Theme (Arial, No Titles)
theme_pub <- theme_bw(base_size = 12, base_family = "Arial") +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    axis.title = element_blank(), # Removes X/Y labels
    legend.title = element_blank(), # Removes legend title
    legend.position = "right",
    text = element_text(family = "Arial")
  )

L <- 5800
p_vals <- c(0.10, 0.20, 0.30)
identities <- seq(0.80, 1.00, by = 0.001)

df_fig1 <- expand.grid(identity = identities, p_lethal = p_vals) %>%
  mutate(
    k = (1 - identity) * L,
    viability = (1 - p_lethal)^k,
    p_label = paste0("p_lethal = ", p_lethal)
  )

# Empirical Point
df_point <- data.frame(identity = 0.97, viability = 0.06)

gg1 <- ggplot(df_fig1, aes(x = identity * 100, y = viability, colour = p_label)) +
  geom_line(linewidth = 1.2) +
  
  # Add Empirical Point
  geom_point(data = df_point, aes(x = identity * 100, y = viability),
             colour = "black", size = 4, inherit.aes = FALSE) +
  
  # Log Scale
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  
  # Colors consistent with previous figures (Blue/Orange/Gold style if preferred, or distinct)
  scale_color_manual(values = c("#2E9FDF", "#FC4E07", "grey50")) +
  
  # Clean Labels (None)
  labs(x = NULL, y = NULL, colour = NULL, title = NULL) +
  
  theme_pub +
  theme(legend.position = c(0.2, 0.2),
        legend.background = element_blank())

print(gg1)
ggsave(paste0(out_path, "fig2_viability_vs_identity.png"), gg1, width = 5, height = 4, dpi = 300)

# Parameters
L <- 5800; ident <- 0.97
viable <- 18; total <- 300
p_obs <- viable / total

# CI Calculation
se <- sqrt(p_obs * (1 - p_obs) / total)
p_low <- max(p_obs - 1.96 * se, 1e-10)
p_high <- min(p_obs + 1.96 * se, 1 - 1e-10)

k <- (1 - ident) * L
p_eff <- 1 - p_obs^(1 / k)
p_eff_low <- 1 - p_high^(1 / k)
p_eff_high <- 1 - p_low^(1 / k)

dms_vals <- c(0.10, 0.20, 0.30)

df_fig2 <- tibble::tibble(
  source = c("Model", "DMS", "DMS", "DMS"),
  label  = factor(c(sprintf("Model (%.0f%%)", ident * 100), "DMS: 10%", "DMS: 20%", "DMS: 30%"),
                  levels = c(sprintf("Model (%.0f%%)", ident * 100), "DMS: 10%", "DMS: 20%", "DMS: 30%")),
  p = c(p_eff, dms_vals)
)

df_ci <- tibble::tibble(
  label = factor(sprintf("Model (%.0f%%)", ident * 100)),
  p = p_eff, p_low = p_eff_low, p_high = p_eff_high
)

gg2 <- ggplot(df_fig2, aes(x = label, y = p)) +
  geom_col(fill = "grey40", width = 0.6) +
  
  # Error bar for Model
  geom_errorbar(data = df_ci, aes(ymin = p_low, ymax = p_high), width = 0.2, size = 0.8) +
  
  labs(x = NULL, y = NULL, title = NULL) +
  
  theme_pub +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

print(gg2)
ggsave(paste0(out_path, "fig2_effective_p.png"), gg2, width = 4, height = 4, dpi = 300)

identity_fixed <- 0.97
L_seq <- seq(1000, 15000, by = 250)

df_fig3 <- expand.grid(L = L_seq, p_lethal = p_vals) %>%
  mutate(
    k = (1 - identity_fixed) * L,
    viability = (1 - p_lethal)^k,
    p_label = paste0("p_lethal = ", p_lethal)
  )

gg3 <- ggplot(df_fig3, aes(x = L / 1000, y = viability, colour = p_label)) +
  geom_line(linewidth = 1.2) +
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_color_manual(values = c("#2E9FDF", "#FC4E07", "grey50")) +
  
  labs(x = NULL, y = NULL, colour = NULL, title = NULL) +
  theme_pub +
  theme(legend.position = "none") # Legend redundant if shown in gg1

print(gg3)
ggsave(paste0(out_path, "fig2_viability_vs_length.png"), gg3, width = 5, height = 4, dpi = 300)



set.seed(123)
theta <- seq(0, 2 * pi, length.out = 400)
outer_radius <- 1

df_outer <- data.frame(x = outer_radius * cos(theta), y = outer_radius * sin(theta))

# Inner "manifold"
inner_radius_x <- 0.35
inner_radius_y <- 0.20
df_inner <- data.frame(
  x = inner_radius_x * cos(theta) + 0.25,
  y = inner_radius_y * sin(theta) - 0.10
)

# Random Points
n_rand <- 120
rand_pts <- tibble(x = runif(n_rand*2, -1, 1), y = runif(n_rand*2, -1, 1)) %>%
  filter(x^2 + y^2 <= 1) %>% dplyr::slice(1:n_rand) %>% mutate(generator = "Random")

# Model Points
n_model <- 60
u <- runif(n_model, 0, 2 * pi); r <- sqrt(runif(n_model, 0, 1))
model_pts <- tibble(
  x = 0.25 + inner_radius_x * r * cos(u),
  y = -0.10 + inner_radius_y * r * sin(u),
  generator = "Generative Model"
)

df_pts <- bind_rows(rand_pts, model_pts)

gg_schema <- ggplot() +
  geom_path(data = df_outer, aes(x = x, y = y), linewidth = 0.8, color = "grey30") +
  geom_polygon(data = df_inner, aes(x = x, y = y), fill = "#FC4E07", alpha = 0.1) +
  geom_path(data = df_inner, aes(x = x, y = y), linewidth = 0.5, linetype = "dashed", color = "#FC4E07") +
  
  geom_point(data = df_pts, aes(x = x, y = y, color = generator), size = 2.5, alpha = 0.8) +
  
  scale_color_manual(values = c("Random" = "grey60", "Generative Model" = "#FC4E07")) +
  
  # Text Annotations
  annotate("text", x = 0, y = 1.1, label = "Sequence Space (97% Id)", size = 4, family = "Arial") +
  annotate("text", x = 0.35, y = -0.4, label = "Viability Manifold", size = 4, color = "#FC4E07", family = "Arial", fontface="bold") +
  
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none")

print(gg_schema)
ggsave(paste0(out_path, "fig2_schematic_manifold.png"), gg_schema, width = 5, height = 5, dpi = 300)

L <- 5800; P_model <- 0.06
identities <- c(0.95, 0.96, 0.97, 0.98, 0.99)
p_baseline <- seq(0.001, 0.40, by = 0.005)

df_uplift <- expand.grid(identity = identities, p_lethal = p_baseline) %>%
  mutate(
    k = (1 - identity) * L,
    P_baseline = (1 - p_lethal)^k,
    uplift = P_model / P_baseline,
    identity_label = paste0(identity * 100, "%")
  )

gg_uplift <- ggplot(df_uplift, aes(x = p_lethal, y = uplift, colour = identity_label)) +
  geom_line(linewidth = 1.2) +
  scale_x_log10() + 
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  
  labs(x = NULL, y = NULL, colour = NULL) +
  
  theme_pub +
  theme(legend.position = c(0.8, 0.7),
        legend.background = element_rect(fill="white", color=NA))

print(gg_uplift)
ggsave(paste0(out_path, "fig2_uplift.png"), gg_uplift, width = 5, height = 4, dpi = 300)

# Note: This block assumes 'natural_genomes_aln.fasta' exists.
# If you don't have it, set f_var manually (e.g., f_var <- 0.3) to run the plot logic.

#f_var <- 0.3 # Placeholder if file missing; Replace with calculated value if available

# --- Uncomment to calculate from file ---
aln <- readDNAStringSet("~/Downloads/natural_genomes_aln.fasta")
aln_mat <- do.call(rbind, lapply(as.character(aln), function(s) strsplit(s, "")[[1]]))
summ_cols <- function(col) { bases <- col[!col %in% c("-", "N", "n")]; length(unique(bases)) }
n_unique <- apply(aln_mat, 2, summ_cols)
f_var <- mean(n_unique >= 2)
# ---------------------------------------

G <- 5800; P_model <- 0.06
identities <- c(0.95, 0.96, 0.97, 0.98, 0.99)
p_grid <- seq(0.001, 0.30, by = 0.001)

df_msa <- expand.grid(identity = identities, p_lethal_msa = p_grid) %>%
  mutate(
    m = 1 - identity,
    feasible = m <= f_var,
    k_allowed = ifelse(feasible, (m / f_var) * G, NA),
    P_baseline = ifelse(feasible, (1 - p_lethal_msa)^k_allowed, 0),
    uplift = ifelse(P_baseline > 0, P_model / P_baseline, NA),
    identity_label = paste0(round(identity * 100), "%")
  )

gg_msa <- ggplot(df_msa %>% filter(feasible), aes(x = p_lethal_msa, y = uplift, colour = identity_label)) +
  geom_line(linewidth = 1.2) +
  scale_x_log10() + 
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  
  labs(x = NULL, y = NULL, colour = NULL) +
  theme_pub +
  theme(legend.position = "none") # Consistent style

print(gg_msa)
ggsave(paste0(out_path, "fig2_uplift_msa.png"), gg_msa, width = 5, height = 4, dpi = 300)

# Define Path
out_path <- "~/Documents/PhD/Projects/Biosecurity/bacteriophage_paper/"

# Clean Theme for Data Plots (Box, No Text)
theme_clean <- theme_bw(base_size = 12, base_family = "Arial") +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    legend.background = element_blank()
  )

# Void Theme for Schematic (No Box at all)
theme_void_custom <- theme_void(base_family = "Arial") +
  theme(legend.position = "none")

# --- COLOR PALETTES ---
# 1. Comparison Colors
col_model <- "#D55E00"  # Vermilion (High Vis)
  col_ref   <- "#999999"  # Grey
    
  # 2. P-value Palette (Viability plots)
  # Mapping specific p-values to distinct colors
  cols_p_lethal <- c(
    "p = 0.1" = "#009E73", # Green (Optimistic)
    "p = 0.2" = "#56B4E9", # Light Blue (Moderate)
    "p = 0.3" = "#0072B2"  # Dark Blue (Harsh/Random)
  )
  
  # 3. Identity Palette (Uplift plots)
  # Use Viridis for continuous feel but discrete steps
  cols_identity <- scales::viridis_pal(option = "plasma", end = 0.8)(5)
  names(cols_identity) <- c("95%", "96%", "97%", "98%", "99%")
  
  
  # ==============================================================================
  # 1. VIABILITY vs IDENTITY (gg1)
  # ==============================================================================
  L <- 5800
  p_vals <- c(0.10, 0.20, 0.30)
  identities <- seq(0.80, 1.00, by = 0.001)
  
  df_fig1 <- expand.grid(identity = identities, p_lethal = p_vals) %>%
    mutate(
      k = (1 - identity) * L,
      viability = (1 - p_lethal)^k,
      p_label = paste0("p = ", p_lethal)
    )
  
  df_point <- data.frame(identity = 0.97, viability = 0.06)
  
  # Base Plot
  p_gg1 <- ggplot(df_fig1, aes(x = identity * 100, y = viability, colour = p_label)) +
    geom_line(linewidth = 1.2) +
    geom_point(data = df_point, aes(x = identity * 100, y = viability),
               colour = "black", size = 2.5, inherit.aes = FALSE) +
    scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_color_manual(values = cols_p_lethal, name = "Lethal Rate") +
    theme_clean
  
  # 1a. Save Legend
  legend_gg1 <- get_legend(p_gg1 + theme(legend.box.margin = margin(10, 10, 10, 10)))
  ggsave(paste0(out_path, "legend_viability.png"), legend_gg1, width = 2, height = 2, dpi=300)
  
  # 1b. Save Plot (No Legend)
  gg1_final <- p_gg1 + theme(legend.position = "none")
  ggsave(paste0(out_path, "fig2_viability_vs_identity.png"), gg1_final, width = 3, height = 3, dpi = 300)
  
  
  # ==============================================================================
  # 2. EFFECTIVE P (gg2)
  # ==============================================================================
  L <- 5800; ident <- 0.97
  viable <- 18; total <- 300
  p_obs <- viable / total
  
  se <- sqrt(p_obs * (1 - p_obs) / total)
  p_low <- max(p_obs - 1.96 * se, 1e-10)
  p_high <- min(p_obs + 1.96 * se, 1 - 1e-10)
  
  k <- (1 - ident) * L
  p_eff <- 1 - p_obs^(1 / k)
  p_eff_low <- 1 - p_high^(1 / k)
  p_eff_high <- 1 - p_low^(1 / k)
  
  dms_vals <- c(0.10, 0.20, 0.30)
  
  df_fig2 <- tibble(
    label  = factor(c("Model", "DMS 10%", "DMS 20%", "DMS 30%"),
                    levels = c("Model", "DMS 10%", "DMS 20%", "DMS 30%")),
    p = c(p_eff, dms_vals),
    type = c("Model", "Ref", "Ref", "Ref")
  )
  
  df_ci <- tibble(
    label = factor("Model", levels = levels(df_fig2$label)),
    p_low = p_eff_low, p_high = p_eff_high
  )
  
  gg2 <- ggplot(df_fig2, aes(x = label, y = p, fill = type)) +
    geom_col(width = 0.6) +
    geom_errorbar(data = df_ci, aes(x = label, ymin = p_low, ymax = p_high), 
                  width = 0.2, inherit.aes = FALSE, size = 0.8) +
    scale_fill_manual(values = c("Model" = col_model, "Ref" = col_ref)) +
    theme_clean +
    theme(legend.position = "none") # No legend needed for bar chart usually
  
  ggsave(paste0(out_path, "fig2_effective_p.png"), gg2, width = 2.5, height = 2.5, dpi = 300)
  
  
  # ==============================================================================
  # 3. VIABILITY vs LENGTH (gg3)
  # ==============================================================================
  identity_fixed <- 0.97
  L_seq <- seq(1000, 15000, by = 250)
  
  df_fig3 <- expand.grid(L = L_seq, p_lethal = p_vals) %>%
    mutate(
      k = (1 - identity_fixed) * L,
      viability = (1 - p_lethal)^k,
      p_label = paste0("p = ", p_lethal)
    )
  
  p_gg3 <- ggplot(df_fig3, aes(x = L / 1000, y = viability, colour = p_label)) +
    geom_line(linewidth = 1.2) +
    scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_color_manual(values = cols_p_lethal, name = "Lethal Rate") +
    theme_clean
  
  # 3a. Save Legend (Same as gg1, but saved separately as requested)
  legend_gg3 <- get_legend(p_gg3 + theme(legend.box.margin = margin(10, 10, 10, 10)))
  ggsave(paste0(out_path, "legend_length.png"), legend_gg3, width = 2, height = 2, dpi=300)
  
  # 3b. Save Plot
  gg3_final <- p_gg3 + theme(legend.position = "none")
  ggsave(paste0(out_path, "fig2_viability_vs_length.png"), gg3_final, width = 3, height = 3, dpi = 300)
  
  
  # ==============================================================================
  # 4. SCHEMATIC MANIFOLD (gg_schema) - Removed Box/Grid
  # ==============================================================================
  set.seed(123)
  theta <- seq(0, 2 * pi, length.out = 400)
  outer_radius <- 1
  df_outer <- data.frame(x = outer_radius * cos(theta), y = outer_radius * sin(theta))
  
  inner_radius_x <- 0.30; inner_radius_y <- 0.15
  df_inner <- data.frame(x = inner_radius_x * cos(theta) + 0.25, y = inner_radius_y * sin(theta) - 0.10)
  
  n_rand <- 120
  rand_pts <- tibble(x = runif(n_rand*2, -1, 1), y = runif(n_rand*2, -1, 1)) %>%
    filter(x^2 + y^2 <= 1) %>% dplyr::slice(1:n_rand) %>% mutate(generator = "Random")
  
  n_model <- 60
  u <- runif(n_model, 0, 2 * pi); r <- sqrt(runif(n_model, 0, 1))
  model_pts <- tibble(x = 0.25 + inner_radius_x * r * cos(u), y = -0.10 + inner_radius_y * r * sin(u), generator = "Generative Model")
  
  df_pts <- bind_rows(rand_pts, model_pts)
  
  gg_schema <- ggplot() +
    # Shapes
    geom_path(data = df_outer, aes(x = x, y = y), linewidth = 0.8, color = "grey40") +
    geom_polygon(data = df_inner, aes(x = x, y = y), fill = col_model, alpha = 0.1) +
    geom_path(data = df_inner, aes(x = x, y = y), linewidth = 0.5, linetype = "dashed", color = col_model) +
    
    # Points
    geom_point(data = df_pts, aes(x = x, y = y, color = generator), size = 2.0, alpha = 0.8) +
    
    # Colors
    scale_color_manual(values = c("Random" = "grey60", "Generative Model" = col_model)) +
    
    coord_fixed() +
    theme_void_custom # Uses the void theme to remove all boxes/grids
  
  ggsave(paste0(out_path, "fig2_schematic_manifold.png"), gg_schema, width = 4, height = 4, dpi = 300)
  
  
  # ==============================================================================
  # 5. UPLIFT CURVES (gg_uplift)
  # ==============================================================================
  L <- 5800; P_model <- 0.06
  identities <- c(0.95, 0.96, 0.97, 0.98, 0.99)
  p_baseline <- seq(0.001, 0.40, by = 0.005)
  
  df_uplift <- expand.grid(identity = identities, p_lethal = p_baseline) %>%
    mutate(
      k = (1 - identity) * L,
      P_baseline = (1 - p_lethal)^k,
      uplift = P_model / P_baseline,
      identity_label = paste0(identity * 100, "%")
    )
  
  p_uplift <- ggplot(df_uplift, aes(x = p_lethal, y = uplift, colour = identity_label)) +
    geom_line(linewidth = 1.2) +
    scale_x_log10() + 
    scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_color_manual(values = cols_identity, name = "Identity") +
    theme_clean
  
  # 5a. Save Legend
  legend_uplift <- get_legend(p_uplift + theme(legend.box.margin = margin(10, 10, 10, 10)))
  ggsave(paste0(out_path, "legend_uplift.png"), legend_uplift, width = 2, height = 2, dpi=300)
  
  # 5b. Save Plot
  gg_uplift_final <- p_uplift + theme(legend.position = "none")
  ggsave(paste0(out_path, "fig2_uplift.png"), gg_uplift_final, width = 3, height = 3, dpi = 300)
  
  
  # ==============================================================================
  # 6. MSA-AWARE UPLIFT (gg_msa)
  # ==============================================================================
  if(!exists("f_var")) f_var <- 0.3 
  
  G <- 5800; P_model <- 0.06
  identities <- c(0.95, 0.96, 0.97, 0.98, 0.99)
  p_grid <- seq(0.001, 0.30, by = 0.001)
  
  df_msa <- expand.grid(identity = identities, p_lethal_msa = p_grid) %>%
    mutate(
      m = 1 - identity,
      feasible = m <= f_var,
      k_allowed = ifelse(feasible, (m / f_var) * G, NA),
      P_baseline = ifelse(feasible, (1 - p_lethal_msa)^k_allowed, 0),
      uplift = ifelse(P_baseline > 0, P_model / P_baseline, NA),
      identity_label = paste0(round(identity * 100), "%")
    )
  
  p_msa <- ggplot(df_msa %>% filter(feasible), aes(x = p_lethal_msa, y = uplift, colour = identity_label)) +
    geom_line(linewidth = 1.2) +
    scale_x_log10() + 
    scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_color_manual(values = cols_identity, name = "Identity") +
    theme_clean
  
  # Save Plot (Reuse uplift legend as colors are identical)
  gg_msa_final <- p_msa + theme(legend.position = "none")
  ggsave(paste0(out_path, "fig2_uplift_msa.png"), gg_msa_final, width = 3, height = 3, dpi = 300)
  
  
###########################################
###################### 1c #################
###########################################  
  
  
  x <- seq(0.001, 1.6, length.out = 1000)
  
  # Model each DFE as mixture:
  # Component 1: exponential decay from 0 (lethal/deleterious)
  # Component 2: bell-shaped peak near 1 (neutral/viable)
  
  # Standard DFE
  w1_std <- 0.50  # weight on deleterious component
  rate_std <- 12  # steeper exponential = sharper peak at 0
  mu_std <- 0.95  # peak location for viable
  sd_std <- 0.08  # spread of viable peak
  
  # Supercharged DFE  
  w1_sup <- 0.80  # much more weight on deleterious
  rate_sup <- 15  # even steeper = taller narrower spike at 0
  mu_sup <- 1.0   # peak slightly higher
  sd_sup <- 0.18  # wider spread = fatter tails both sides
  
  # Calculate densities
  y1_exp <- dexp(x, rate = rate_std)
  y1_norm <- dnorm(x, mean = mu_std, sd = sd_std)
  y1 <- w1_std * y1_exp + (1 - w1_std) * y1_norm
  
  y2_exp <- dexp(x, rate = rate_sup)
  y2_norm <- dnorm(x, mean = mu_sup, sd = sd_sup)
  y2 <- w1_sup * y2_exp + (1 - w1_sup) * y2_norm
  
  # Normalize to equal area
  dx <- x[2] - x[1]
  y1 <- y1 / (sum(y1) * dx)
  y2 <- y2 / (sum(y2) * dx)
  
  # Create data frame
  df <- data.frame(
    x = rep(x, 2),
    y = c(y1, y2),
    distribution = rep(c("standard", "supercharged"), each = length(x))
  )
  
  # Create the plot
  p <- ggplot(df, aes(x = x, y = y, fill = distribution)) +
    geom_area(alpha = 0.5, position = "identity") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
 #   geom_vline(xintercept = 1.3, linetype = "dotted", color = "black", linewidth = 0.5) +
    scale_fill_manual(values = c("standard" = "#4477AA", "supercharged" = "#CC6677")) +
    scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 1.5), ylim = c(0, max(c(y1, y2)) * 1.05)) +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      legend.position = "none",
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    )
  
  # Display the plot
  print(p)

ggsave(paste0(out_path, "cartoon_fitness.png"), p, width = 4, height = 2.5, dpi = 300)
  
 