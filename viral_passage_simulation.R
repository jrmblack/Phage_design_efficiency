# ------------------------------------------------------------
# Serial passage simulator (minimal, fast, and interpretable)
# - Asexual haploid population
# - Within each passage: g generations of selection + mutation
# - Between passages: bottleneck of size B
# - Tracks mutation counts in 3 classes: deleterious / neutral / beneficial
# - Reports *sequence identity* using Jukes–Cantor (JC69), not 1 - k/L
#
# Defaults are "ΦX174-like": L = 5386, mu ~ 1e-6 per site per copying round
# ------------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

# JC69: substitutions/site (d) -> expected observed identity
jc_identity_from_d <- function(d) 0.25 + 0.75 * exp(-4 * d / 3)
# (and p-distance would be 1 - identity)

simulate_serial_passage_phiX174 <- function(
    P = 1000,             # number of passages
    g = 5,                # effective generations per passage
    L = 5386,             # ΦX174 genome length
    mu = 1e-6,            # per-site mutation rate per generation
    N = 1e7,            # within-passage effective population size
    B = 200,              # bottleneck founders per passage
    # DFE mixture (toy but useful)
    p_del = 0.70,
    p_neu = 0.299,
    p_ben = 0.001,
    s_del = 0.02,         # fitness penalty per deleterious mutation (multiplicative via exp)
    s_ben = 0.05,         # fitness gain per beneficial mutation (multiplicative via exp)
    seed = 1
) {
  stopifnot(abs((p_del + p_neu + p_ben) - 1) < 1e-12)
  set.seed(seed)
  
  lambda <- L * mu  # per-genome mutation rate per generation (≈0.005386 for defaults)
  
  # We exploit that lambda is small for ΦX174: P(>=2 mutations/genome/gen) is tiny.
  # So we use a "0-or-1 mutation per genome per generation" approximation:
  #   mutated genomes ~ Binomial(n, lambda)
  # This is accurate for small lambda and keeps the simulator fast.
  #
  # If you later want multi-hit handling, tell me and I'll add an efficient version.
  
  parse_key <- function(key) as.integer(strsplit(key, "\\|", fixed = FALSE)[[1]])
  make_key  <- function(kd, kn, kb) paste0(kd, "|", kn, "|", kb)
  fitness   <- function(kd, kb) exp(-s_del * kd + s_ben * kb)
  
  # population as hash map: key "kd|kn|kb" -> count
  pop <- new.env(hash = TRUE, parent = emptyenv())
  pop[["0|0|0"]] <- B  # start with B identical founders
  
  summarize_pop <- function(pop_env) {
    keys <- ls(pop_env)
    counts <- as.integer(mget(keys, pop_env, inherits = FALSE))
    k_mat <- t(vapply(keys, parse_key, integer(3)))
    k_tot <- rowSums(k_mat)
    
    Ntot <- sum(counts)
    mean_k <- sum(counts * k_tot) / Ntot
    
    mode_idx <- which.max(counts)
    mode_k <- k_tot[mode_idx]
    
    # substitutions/site (d)
    mean_d <- mean_k / L
    mode_d <- mode_k / L
    
    # expected observed identity (JC69)
    mean_id_jc <- jc_identity_from_d(mean_d)
    mode_id_jc <- jc_identity_from_d(mode_d)
    
    list(
      N = Ntot,
      mean_k = mean_k,
      mode_k = mode_k,
      mean_d = mean_d,
      mode_d = mode_d,
      mean_identity_jc = mean_id_jc,
      mode_identity_jc = mode_id_jc
    )
  }
  
  out <- data.frame(
    passage = 0:P,
    N = NA_integer_,
    mean_k = NA_real_,
    mode_k = NA_real_,
    mean_d = NA_real_,
    mode_d = NA_real_,
    mean_identity_jc = NA_real_,
    mode_identity_jc = NA_real_
  )
  
  # init summary
  s0 <- summarize_pop(pop)
  out[1, names(s0)] <- s0
  
  for (p in seq_len(P)) {
    
    # Within passage: g generations
    for (gen in seq_len(g)) {
      keys <- ls(pop)
      counts <- as.integer(mget(keys, pop, inherits = FALSE))
      k_mat <- t(vapply(keys, parse_key, integer(3)))
      
      # --- Selection (Wright–Fisher with selection) ---
      w <- fitness(kd = k_mat[, 1], kb = k_mat[, 3])
      probs <- counts * w
      probs <- probs / sum(probs)
      next_counts <- as.integer(rmultinom(1, size = N, prob = probs))
      
      # --- Mutation (0-or-1 per genome approximation) ---
      pop2 <- new.env(hash = TRUE, parent = emptyenv())
      
      for (i in seq_along(keys)) {
        n_i <- next_counts[i]
        if (n_i == 0L) next
        
        kd <- k_mat[i, 1]; kn <- k_mat[i, 2]; kb <- k_mat[i, 3]
        
        # how many of these n_i genomes mutate this generation?
        m_i <- rbinom(1, size = n_i, prob = lambda)
        n_nomut <- n_i - m_i
        
        # keep non-mutated genomes in same state
        if (n_nomut > 0L) {
          key0 <- make_key(kd, kn, kb)
          pop2[[key0]] <- (pop2[[key0]] %||% 0L) + n_nomut
        }
        
        # mutated genomes: exactly 1 mutation each; split by DFE class
        if (m_i > 0L) {
          splits <- as.integer(rmultinom(1, size = m_i, prob = c(p_del, p_neu, p_ben)))
          
          # each class produces a new state with +1 in that class
          if (splits[1] > 0L) {
            keyd <- make_key(kd + 1L, kn, kb)
            pop2[[keyd]] <- (pop2[[keyd]] %||% 0L) + splits[1]
          }
          if (splits[2] > 0L) {
            keyn <- make_key(kd, kn + 1L, kb)
            pop2[[keyn]] <- (pop2[[keyn]] %||% 0L) + splits[2]
          }
          if (splits[3] > 0L) {
            keyb <- make_key(kd, kn, kb + 1L)
            pop2[[keyb]] <- (pop2[[keyb]] %||% 0L) + splits[3]
          }
        }
      }
      
      pop <- pop2
    }
    
    # --- Bottleneck between passages ---
    keys <- ls(pop)
    counts <- as.integer(mget(keys, pop, inherits = FALSE))
    probs <- counts / sum(counts)
    founders <- as.integer(rmultinom(1, size = B, prob = probs))
    
    popB <- new.env(hash = TRUE, parent = emptyenv())
    for (i in seq_along(keys)) {
      if (founders[i] > 0L) popB[[keys[i]]] <- founders[i]
    }
    pop <- popB
    
    sp <- summarize_pop(pop)
    out[p + 1, names(sp)] <- sp
  }
  
  out
}

# ---------------- Example ----------------
res <- simulate_serial_passage_phiX174(
  P = 2000,
  g = 10,
  N = 1e7,
  B = 1000,
  L = 5386,
  mu = 1e-6,
  p_del = 0.7,
  p_neu = 0.29,
  p_ben = 0.01,
  s_del = 0.2,
  s_ben = 0.05,
  seed = 1
)

# First passage where "consensus-like" (mode) expected identity <= 97%
hit <- which(res$mode_identity_jc <= 0.97)
if (length(hit)) {
  cat("First passage with mode identity <= 0.97:", res$passage[min(hit)], "\n")
} else {
  cat("Mode identity never <= 0.97 within", max(res$passage), "passages.\n")
}

# Plot mean vs mode identity (JC-corrected)
plot(res$passage, res$mean_identity_jc, type = "l",
     xlab = "Passage", ylab = "Expected sequence identity (JC69)",
     main = "Serial passage: ΦX174-like parameters")
lines(res$passage, res$mode_identity_jc)
abline(h = 0.97, lty = 2)
legend("topright",
       legend = c("Mean identity (JC69)", "Mode identity (JC69)", "97% threshold"),
       lty = c(1, 1, 2), bty = "n")

df_res <- data.frame(passage = res$passage, mean_identity_jc = res$mean_identity_jc)

gg1 <- df_res %>%
ggplot(aes(x = passage, y = mean_identity_jc)) + geom_line() + 
  labs(x = "Passage", y = "Expected Sequence identity (JC69)", title = "Serial passage: ΦX174-like parameters") + 
  theme_bw() + 
  geom_hline(yintercept = 0.97, lty = 2)

# ------------------------------------------------------------
# Grid runner for serial passage simulator
# - Expands a parameter grid
# - Runs multiple seeds per condition
# - Computes time-to-threshold T97 (first passage where identity <= 0.97)
# - Returns tidy results + optional trajectory summaries
#
# Assumes you have simulate_serial_passage_phiX174() defined (from earlier)
# and that its output has columns: passage, mean_identity_jc, mode_identity_jc
# ------------------------------------------------------------

library(dplyr)
library(tidyr)
library(purrr)

# ---- helpers ----

first_hit_time <- function(df, col = c("mode_identity_jc", "mean_identity_jc"), thr = 0.97) {
  col <- match.arg(col)
  idx <- which(df[[col]] <= thr)
  if (length(idx) == 0) return(NA_integer_)
  df$passage[min(idx)]
}

run_one_condition <- function(params, seeds, thr = 0.97) {
  # params: one-row data.frame/tibble with parameter columns
  # seeds: integer vector of RNG seeds
  
  p <- as.list(params)
  
  out <- purrr::map_dfr(seeds, function(seed) {
    sim <- do.call(simulate_serial_passage_phiX174, c(p, list(seed = seed)))
    
    tibble::tibble(
      seed = seed,
      T97_mode = first_hit_time(sim, "mode_identity_jc", thr),
      T97_mean = first_hit_time(sim, "mean_identity_jc", thr),
      final_mode_id = sim$mode_identity_jc[nrow(sim)],
      final_mean_id = sim$mean_identity_jc[nrow(sim)]
    )
  })
  
  # replicate the 1-row params to match number of rows in out
  params_rep <- params[rep(1, nrow(out)), , drop = FALSE]
  dplyr::bind_cols(out, params_rep)
}


# Optional: summarize trajectories across seeds (median + bands) for plotting
summarise_trajectories <- function(trajs_df, col, probs = c(0.1, 0.5, 0.9)) {
  # trajs_df: tidy with columns: passage, value, plus condition cols + seed
  # col: name of value column (already selected), e.g. "mode_identity_jc"
  stopifnot(col %in% names(trajs_df))
  
  trajs_df %>%
    group_by(across(-c(seed, all_of(col)))) %>%
    summarise(
      qlo = quantile(.data[[col]], probs[1], na.rm = TRUE),
      q50 = quantile(.data[[col]], probs[2], na.rm = TRUE),
      qhi = quantile(.data[[col]], probs[3], na.rm = TRUE),
      .groups = "drop"
    )
}

require(dplyr)
require(tidyr)
require(purrr)
# ---- build parameter grid ----
# ΦX174-ish defaults, varying the main axes:
# B (bottleneck), mu (mutation rate), p_del (constraint density), g (gens/passage)
param_grid <- tidyr::crossing(
  P = 2000L,
  L = 5386L,
  N = 1e7L,
  B = c(5L, 20L, 100L, 1000L),
  g = c(3L, 5L, 10L),
  mu = c(3e-7, 1e-6, 3e-6),
  p_del = c(0.3, 0.5, 0.7),
  # keep beneficial minimal; can also toggle p_ben = c(0, 1e-3) if you want
  p_ben = 1e-3,
  s_del = c(0.1, 0.2, 0.5),
  s_ben = c(0.05, 0.1, 0.2)
) %>%
  mutate(
    p_neu = 1 - p_del - p_ben
  ) %>%
  # basic sanity check
  filter(p_neu >= 0)

# ---- choose seeds ----
n_seeds <- 3
seeds <- 1:n_seeds

# ---- run grid (summary outcomes only) ----
# This returns one row per (condition × seed)
results <- param_grid %>%
  mutate(cond_id = row_number()) %>%
  group_split(cond_id) %>%
  map_dfr(~run_one_condition(.x %>% select(-cond_id), seeds = seeds, thr = 0.97)) %>%
  mutate(
    hit_mode = !is.na(T97_mode),
    hit_mean = !is.na(T97_mean)
  )

saveRDS(results, "simulation_serial_passage_results.rds")

# ---- summarize per condition ----
summary_by_condition <- results %>%
  group_by(P, L, N, B, g, mu, p_del, p_neu, p_ben, s_del, s_ben) %>%
  summarise(
    n = n(),
    hit_rate_mode = mean(hit_mode),
    hit_rate_mean = mean(hit_mean),
    T97_mode_median = median(T97_mode, na.rm = TRUE),
    T97_mode_q10 = quantile(T97_mode, 0.10, na.rm = TRUE),
    T97_mode_q90 = quantile(T97_mode, 0.90, na.rm = TRUE),
    T97_mean_median = median(T97_mean, na.rm = TRUE),
    T97_mean_q10 = quantile(T97_mean, 0.10, na.rm = TRUE),
    T97_mean_q90 = quantile(T97_mean, 0.90, na.rm = TRUE),
    final_mode_id_median = median(final_mode_id, na.rm = TRUE),
    final_mean_id_median = median(final_mean_id, na.rm = TRUE),
    .groups = "drop"
  )

library(dplyr)

P_max <- max(summary_by_condition$P)

plot_df <- summary_by_condition %>%
  mutate(
    T97_mode_plot = ifelse(is.na(T97_mode_median), P_max + 1, T97_mode_median),
    T97_mean_plot = ifelse(is.na(T97_mean_median), P_max + 1, T97_mean_median),
    censored_mode = is.na(T97_mode_median),
    censored_mean = is.na(T97_mean_median)
  )


library(ggplot2)

ggplot(plot_df, aes(x = factor(B), y = T97_mode_plot, group = mu, color = factor(mu))) +
  geom_point() +
  geom_line(aes(group = factor(mu))) +
  geom_errorbar(aes(ymin = T97_mode_q10, ymax = T97_mode_q90), width = 0.2, alpha = 0.6) +
  facet_grid(p_del ~ g, labeller = label_both) +
  scale_y_continuous(
    trans = "log10",
    breaks = c(10, 30, 100, 300, 1000, 3000, 10000),
    minor_breaks = NULL
  ) +
  labs(
    x = "Bottleneck size (B)",
    y = "Passages to reach 97% identity (mode, median; q10–q90)",
    color = expression(mu~"(per-site per generation)"),
    title = "Time to 97% identity depends strongly on bottleneck and mutation rate"
  ) +
  theme_bw()

ggplot(plot_df, aes(x = factor(B), y = factor(mu), fill = T97_mode_plot)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_grid(p_del ~ g, labeller = label_both) +
  scale_fill_continuous(trans = "log10") +
  labs(
    x = "Bottleneck size (B)",
    y = expression(mu),
    fill = "T97 (mode)\n(passages)",
    title = "Heatmap of passages to 97% identity (censored at P+1 if not reached)"
  ) +
  theme_bw()

ggplot(plot_df, aes(x = factor(B), y = factor(mu))) +
  geom_tile(aes(fill = T97_mode_plot), color = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.2f", hit_rate_mode)), size = 3) +
  facet_grid(p_del ~ g, labeller = label_both) +
  scale_fill_continuous(trans = "log10") +
  labs(
    x = "Bottleneck size (B)",
    y = expression(mu),
    fill = "T97 (mode)\n(passages)",
    title = "T97 heatmap (numbers are hit-rate across seeds)"
  ) +
  theme_bw()

ggplot(plot_df, aes(x = factor(B), y = final_mode_id_median, color = factor(mu))) +
  geom_point() +
  geom_line(aes(group = factor(mu))) +
  facet_grid(p_del ~ g, labeller = label_both) +
  geom_hline(yintercept = 0.97, linetype = 2) +
  labs(
    x = "Bottleneck size (B)",
    y = "Median final identity (mode, JC69)",
    color = expression(mu),
    title = "Median identity after 2000 passages"
  ) +
  theme_bw()

summary_by_condition <- results %>%
  group_by(P, L, N, B, g, mu, p_del, p_neu, p_ben, s_del, s_ben) %>%
  summarise(
    final_mode_id_median = median(final_mode_id, na.rm = TRUE),
    final_mode_id_q10    = quantile(final_mode_id, 0.10, na.rm = TRUE),
    final_mode_id_q90    = quantile(final_mode_id, 0.90, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(summary_by_condition,
       aes(x = factor(B), y = final_mode_id_median, color = factor(mu))) +
  geom_point(position = position_dodge(width = 0.4), size = 2.2) +
  facet_grid(p_del ~ g + s_ben, labeller = label_both) +
  geom_hline(yintercept = 0.97, linetype = 2, linewidth = 0.6) +
  labs(
    x = "Bottleneck size (B)",
    y = "Final sequence identity after P passages (median ± 10–90%)",
    color = expression(mu),
    title = "Final identity across bottleneck, mutation rate, and constraint strength"
  ) +
  theme_bw(base_size = 12)

ggplot(summary_by_condition,
       aes(x = factor(B), y = factor(mu), fill = final_mode_id_median)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_grid(p_del ~ g, labeller = label_both) +
  scale_fill_viridis_c(option = "magma", direction = -1) +
  geom_text(aes(label = sprintf("%.3f", final_mode_id_median)), size = 3) +
  geom_hline(yintercept = NA) +
  labs(
    x = "Bottleneck size (B)",
    y = expression(mu),
    fill = "Final identity",
    title = "Final sequence identity after P passages"
  ) +
  theme_bw(base_size = 12)

library(dplyr)
library(ggplot2)

param_cols <- c("B", "g", "mu", "p_del", "s_del", "s_ben", "N")

X <- summary_by_condition %>%
  select(all_of(param_cols)) %>%
  mutate(across(everything(), scale)) %>%
  as.matrix()

pc <- prcomp(X)

pc_df <- summary_by_condition %>%
  mutate(
    PC1 = pc$x[,1],
    PC2 = pc$x[,2]
  )

ggplot(pc_df, aes(PC1, PC2, color = final_mode_id_median)) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_viridis_c(direction = -1) +
  labs(
    title = "Global parameter landscape colored by final identity",
    color = "Final identity"
  ) +
  theme_bw(  )

library(GGally)

summary_by_condition %>%
  select(B, g, mu, p_del, s_del, s_ben, final_mode_id_median) %>%
  ggparcoord(
    columns = 1:6,
    groupColumn = 7,
    scale = "uniminmax",
    alphaLines = 0.3
  ) +
  scale_color_viridis_c(direction = -1) +
  labs(
    title = "Parameter trajectories colored by final identity"
  ) +
  theme_bw()

best_case_df <- summary_by_condition %>%
  group_by(B, g) %>%
  summarise(
    best_identity = min(final_mode_id_median),
    .groups = "drop"
  )

gg2 <- ggplot(best_case_df, aes(x = factor(B), y = best_identity, group = g, color = factor(g))) +
  geom_point(size = 2) +
  geom_line() +
  geom_hline(yintercept = 0.97, linetype = 2) +
  labs(
    x = "Bottleneck size (B)",
    y = "Best-case final identity across all other parameters",
    color = "g (generations / passage)",
    title = "Even best-case parameter combinations rarely reach 97%"
  ) +
  theme_bw()

require(patchwork)
gg1 + gg2


summary_by_condition %>%
  group_by(B, g) %>%
  summarise(
    best_identity = min(final_mode_id_median),   # lowest identity = most divergence
    median_identity = median(final_mode_id_median),
    .groups = "drop"
  ) %>%
  arrange(best_identity)

thr <- 0.97
bg_bar_df <- summary_by_condition %>%
  mutate(
    reached_97 = final_mode_id_median <= thr,
    Bg = paste0("B=", B, ", g=", g)
  ) %>%
  group_by(B, g, Bg, p_del) %>%
  summarise(
    frac_reached = mean(reached_97),
    n = n(),
    .groups = "drop"
  )

ggplot(bg_bar_df, aes(x = factor(B), y = frac_reached)) +
  geom_col(fill = "steelblue") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Bottleneck / generations per passage",
    y = "Fraction reaching ≤97% identity",
    title = "Probability of reaching 97% identity across regimes"
  ) +
  facet_grid(p_del~g) + 
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(ggplot2)
library(scales)

gg3 <- ggplot(bg_bar_df, aes(x = factor(B), y = frac_reached)) +
  geom_col(
    width = 0.7,
    fill = "darkgreen"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  facet_grid(
    p_del ~ g,
    labeller = label_both
  ) +
  labs(
    x = "Bottleneck size (B)",
    y = "Fraction of parameter settings reaching ≤97% identity",
    title = "Regimes that reach 97% identity"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 13, face = "bold"),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )




