# ══════════════════════════════════════════════════════════
#  POST-BROUK REGENERATION CLUSTERING
#  Goal: group 208 CZ field plots into ecologically
#  meaningful regeneration archetypes for iLand initialization - only year 2025, the most recent status
#
#  Input:  env_chars_subplots.gpkg
#  Output: cluster assignments + archetype profiles
#
#  Logic:
#  1. Load & filter Czech Republic, most recent survey year
#  2. Aggregate subplots → plot level
#  3. Derive climate variables from ERA5 (optional, for validation), or from SPEI in czechia?
#  4. Define regeneration categories from field rules
#  5. Cluster remaining plots on composition + structure
#  6. Validate & describe archetypes
# ══════════════════════════════════════════════════════════

library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cluster)      # silhouette
library(factoextra)   # cluster viz
library(RColorBrewer)
library(data.table)
library(terra)
library(ggpubr)
library(ggpmisc)      # run quick model and paste output on teh ggplot
library(mgcv)
library(GGally)

getwd()


#  0. Load data ───────────────────────────────────────── --------------
# read elevation info
plot_v_all_countries <- terra::vect("../../r_post_brouk/outDataShare/Karim_AEF/cleaned/env_chars_plots.gpkg") 
plot_v_cz <- plot_v_all_countries[plot_v_all_countries$country_name == "Czech Republic" & 
                                  plot_v_all_countries$year == 2025, ]
plot_df <- as.data.frame(plot_v_cz) %>% 
  select(-subplot, - plot_id, - country, - region)

# get raw data of stems, species, ... on subplot level
sub_veg <- fread("../../r_post_brouk/outData/dat_full_species_all_countries.csv") %>%
  filter(country_name == "Czechia" & year == 2025) %>%
  mutate(
    plot      = factor(plot)
  )

# export only czechia for Claude analysis
#fwrite(dat_veg, "../../r_post_brouk/outData/dat_full_species_Czechia_2025.csv")
#fwrite(plot_df, "../../r_post_brouk/outData/env_chars_Czechia_2025.csv")



# Any trend with elevation -------------------------------------------------
# is there any trend with elevation? share of spruce, species richness? check again on plot level ----


## low k so it can't overfit; spruce_share is a 0–1 proportion → quasibinomial
g_eff <- gam(effective_numbers ~ s(elevation, k = 4), data = plot_df, method = "REML")
g_spr <- gam(spruce_share ~ s(elevation, k = 4), data = plot_df,
             family = quasibinomial, method = "REML")
summary(g_eff); summary(g_spr)   # read: edf, Deviance explained, p of s(elevation)

## label + plot
gam_lab <- function(m) { s <- summary(m)
sprintf("dev.expl = %.1f%%,  edf = %.2f,  p = %.3f",
        100*s$dev.expl, s$s.table[1,"edf"], s$s.table[1,"p-value"]) }

p1 <- ggplot(plot_df, aes(elevation, effective_numbers)) +
  geom_point(alpha=.5) + geom_smooth(method="gam", formula = y ~ s(x, k=4)) +
  annotate("text", -Inf, Inf, hjust=-0.05, vjust=1.4, label = gam_lab(g_eff))

p2 <- ggplot(plot_df, aes(elevation, spruce_share)) +
  geom_point(alpha=.5) +
  geom_smooth(method="gam", formula = y ~ s(x, k=4),
              method.args = list(family = quasibinomial)) +
  annotate("text", -Inf, Inf, hjust=-0.05, vjust=1.4, label = gam_lab(g_spr))

ggarrange(p1, p2)


# Make extra data - how many subplots are empty? -----------------------------

## get emptiness - how many subplots are empty per plot? ---------------------------------------------------------
## --- emptiness from dat_veg: a subplot is empty if it has zero stems ---
empty_by_plot <- dat_veg |>
  group_by(plot, subplot) |>
  summarise(sub_stems = sum(n, na.rm = TRUE), .groups = "drop") |>
  group_by(plot) |>
  summarise(
    n_subplots     = n_distinct(subplot),
    n_empty_sub    = sum(sub_stems == 0),
    frac_empty_sub = mean(sub_stems == 0),
    all_empty      = all(sub_stems == 0),
    any_empty      = any(sub_stems == 0),
    .groups = "drop"
  ) |>
  mutate(plot = as.character(plot))


### 2. Emptiness summary 
cat("=== Regeneration status ===\n")
cat("ALL subplots empty (truly bare):", sum(empty_by_plot$all_empty),
    sprintf("(%.1f%%)\n", 100 * mean(empty_by_plot$all_empty)))
cat("ANY subplot empty:              ", sum(empty_by_plot$any_empty),
    sprintf("(%.1f%%)\n", 100 * mean(empty_by_plot$any_empty)))
cat("NO empty subplots:              ", sum(!empty_by_plot$any_empty),
    sprintf("(%.1f%%)\n\n", 100 * mean(!empty_by_plot$any_empty)))
cat("Empty-subplot distribution:\n"); print(table(empty_by_plot$n_empty_sub)); cat("\n")


## 1b. Mature-tree presence per plot  (from sub_veg) ---------------------------------------------------------------------
# Mature = dbh >= 10 cm, encoded as vegtype == "mature".
# We flag it at plot level so we can test whether a "residual-overstory"
# archetype emerges (likely tracks lower disturbance severity /
# partial disturbance) rather than assuming it.
mature_by_plot <- sub_veg  %>% 
  filter(!is.na(n), n > 0, vegtype == "mature") %>%
  group_by(plot) %>%
  summarise(n_mature = sum(n), .groups = "drop")

merged_df <- plot_df %>%
  left_join(mature_by_plot, by = join_by(plot)) %>%
  left_join(empty_by_plot, by = join_by(plot)) %>%
  mutate(
    n_mature   = tidyr::replace_na(n_mature, 0),
    has_mature = n_mature > 0
  )

cat("Plots with >=1 residual mature tree:", sum(plots$has_mature), "\n\n")

# I have 10 plots with presence of mature trees

# ---------------------------------------------------------------------
# 3. Rule-based pre-classification: bare / mature-present / to-cluster
# ---------------------------------------------------------------------
# as my groups are rather a gradinet that distinct grousp, i wpould lie to keep them ecologically separated

#  T0_bare   : all subplots empty            -> iLand bare-soil init
#  T_mature  : has_mature & not bare         -> residual overstory archetype
#              (too few plots, n=10, to split further by composition —
#               pooled into one type regardless of understory profile)
#  Everything else -> data-driven clustering (regen only, no legacy trees).
merged_df <- merged_df %>%
  mutate(regen_type = if_else(all_empty, "T0_bare", "to_cluster")) %>% 
  mutate(
    regen_type = case_when(
      all_empty            ~ "T0_bare",
      has_mature            ~ "T_mature",
      TRUE                   ~ "to_cluster"
    )
  )

cat("Pre-classified as T0_bare:  ", sum(merged_df$regen_type == "T0_bare"), "\n")
cat("Pre-classified as T_mature: ", sum(merged_df$regen_type == "T_mature"), "\n")
cat("Remaining for clustering:  ", sum(merged_df$regen_type == "to_cluster"), "\n\n")

#


# PCA: identify teh main drivers 
# get teh main drivers - for PCA ---------------------------
resp_vars <- c("dens_ha", "mean_hgt", "cv_hgt", "frac_empty_sub",
               "sp_richness", #"shannon_sp", "evenness_sp", "effective_numbers",
               "CWM_shade", "CWM_drought", "spruce_share", "IV_coniferous")

pca_df <- merged_df %>%
  filter(regen_type != "T0_bare") %>%     # bare has no meaningful position; join regen_type on first if not present
  select(all_of(resp_vars)) %>%
  drop_na()

# check correlation before PCA — expect the richness/diversity cluster to be tight
print(round(cor(pca_df), 2))

pca_fit <- prcomp(pca_df, scale. = TRUE)
summary(pca_fit)            # variance explained per PC
print(round(pca_fit$rotation[, 1:3], 2))   # loadings on PC1-3

library(ggfortify)
autoplot(pca_fit, data = merged_df  %>% filter(regen_type != "T0_bare") %>% drop_na(all_of(resp_vars)),
         colour = "regen_type",
         loadings = TRUE, loadings.label = TRUE, loadings.colour = "black") +
  scale_colour_brewer(palette = "Set2") +
  theme_bw()






# 4. PAM (k-medoids) clustering of stocked plots
# ---------------------------------------------------------------------
# Variables = what/how-much is growing (response side). NOT management
# (explains how it got there) and NOT site vars (iLand handles internally).
#
# NOTE on redundancy: IV_conif + IV_decid ~= 1 (collinear) -> IV_decid
# dropped by default. spruce_share and IV_conif also overlap; inspect the
# correlation matrix below and prune further if you like.
vars <- c("dens_ha",
          "mean_hgt",
          "spruce_share",
          "sp_richness",
          "frac_empty_sub")

# run clustering on the subset of data - excluded mature trees as this is its own group
clust_data <- merged_df %>%
  filter(regen_type == "to_cluster") %>%
  select(plot,
        all_of(vars)
  ) %>%
  drop_na()

cat("Plots entering clustering:", nrow(clust_data), "\n\n")

X    <- clust_data %>% select(-plot)
cat("=== Correlation of clustering variables (check redundancy) ===\n")
print(round(cor(X), 2)); cat("\n")

X_sc <- scale(X)


### see plots 1:1 to checsk wheather there are any clusters possible ---------------------------------
# the 5 variables that entered PAM — edit names to match your clust_data

pair_df <- merged_df %>%                                    # plots that entered clustering
  select(plot, all_of(vars)) %>% 
  left_join(plots %>% select(plot, has_mature, n_mature), by = "plot") %>%
  mutate(mature = ifelse(has_mature, "mat_present", "mat_abs"))

ggpairs(
  pair_df, columns = vars,
  mapping = aes(colour = mature),
  lower = list(continuous = wrap("points", alpha = 0.6, size = 1.3)),
  upper = list(continuous = wrap("cor", size = 3)),
  diag  = list(continuous = wrap("densityDiag", alpha = 0.4))
) +
  scale_colour_manual(values = c("mat_abs" = "grey80", "mat_present" = "darkgreen")) +
  scale_fill_manual(  values = c("mat_abs" = "grey80", "mat_present" = "darkgreen")) +
  theme_bw()




# ---- 4a. choose k via average silhouette width (PAM = Partitioning Around Medoids ~ alternative to k-menas clustering ----------------
set.seed(42)
sil_scores <- sapply(2:8, function(k) {
  pam(X_sc, k = k, diss = FALSE)$silinfo$avg.width
})
cat("Silhouette (avg width) by k:\n")
for (i in seq_along(sil_scores)) cat(sprintf("  k=%d: %.3f\n", i + 1, sil_scores[i]))
#k_opt <- which.max(sil_scores) + 1
k_opt <- 3 # set manually to two
cat(sprintf("\nOptimal k = %d\n\n", k_opt))

sil_df <- data.frame(k = 2:8, silhouette = sil_scores)

ggplot(sil_df, aes(k, silhouette)) +
    geom_line() + geom_point(size = 3) +
    geom_vline(xintercept = k_opt, linetype = "dashed", colour = "red") +
    labs(title = "Silhouette width by number of clusters (PAM)",
         x = "k", y = "mean silhouette width") +
    theme_bw()


pam_k4 <- pam(X_sc, k = 4, diss = FALSE)
table(pam_k4$clustering)

pam_k3 <- pam(X_sc, k = 3, diss = FALSE)
table(pam_k3$clustering)

pam_k2 <- pam(X_sc, k = 2, diss = FALSE)
table(pam_k2$clustering)


# ---- 4b. final PAM: clusters + medoid (representative real plot) -----
pam_final <- pam(X_sc, k = k_opt, diss = FALSE)
clust_data$cluster_raw <- pam_final$clustering

# medoid plot_id for each cluster -> the plot we extract iLand input from
medoids <- tibble::tibble(
  cluster_raw = seq_len(k_opt),
  medoid_plot = clust_data$plot[pam_final$id.med]   # was clust_data$plot_id
)

cat("=== Medoid (representative) plot per cluster ===\n")
print(as.data.frame(medoids)); cat("\n")

merged_df_cluster <- merged_df %>%
  left_join(clust_data %>% select(plot, cluster_raw), by = "plot") %>%
  mutate(regen_type = case_when(
    regen_type %in% c("T0_bare", "T_mature") ~ regen_type,
    TRUE                                       ~ paste0("T", cluster_raw)
  ))

cat("Final type distribution:\n"); print(table(plots$regen_type)); cat("\n")

### Clusters profiles: get a boxplot -----------------
# do only for regeneration clusters
comparisons <- combn(c("T1", "T2", "T3"), 2, simplify = FALSE)
merged_df_cluster %>%
  select(regen_type, all_of(vars)) %>%
  pivot_longer(all_of(vars), names_to = "variable", values_to = "value") %>%
  ggplot(aes(regen_type, value, fill = regen_type)) +
  geom_boxplot(outlier.alpha = 0.1, size = 0.1) +
#  geom_jitter(alpha = 0.5) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                     p.adjust.method = "BH", label = "p.signif",
                     hide.ns = TRUE, size = 3) +
  facet_wrap(~ variable, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() + theme(legend.position = "none")




# do for environmental gradient as well
env_vars <- c("elevation", "disturbance_severity", "dist_to_edge",
              "planting_intensity", "anti_browsing_intensity")

merged_df_cluster %>%
  select(regen_type, all_of(env_vars)) %>%
  pivot_longer(all_of(env_vars), names_to = "variable", values_to = "value") %>%
  ggplot(aes(regen_type, value, fill = regen_type)) +
  geom_boxplot(outlier.alpha = 0.4) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                     p.adjust.method = "BH", label = "p.signif",
                     hide.ns = TRUE, size = 3) +
  facet_wrap(~ variable, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() + theme(legend.position = "none")


# ---------------------------------------------------------------------
# 5. Describe archetypes  (+ test the "mature present" hypothesis)
# ---------------------------------------------------------------------
summary_vars <- c("dens_ha", "mean_hgt", #"cv_hgt",
                  "spruce_share", "sp_richness",
                  #"effective_numbers", 
                  "CWM_drought",
                  "frac_empty_sub",
                  "clear_intensity", "grndwrk_intensity", "logging_trail_intensity",
                  "planting_intensity", "anti_browsing_intensity",
                  "elevation", "dist_to_edge",
                  "disturbance_severity", "n_mature")

archetype_profiles <- merged_df_cluster %>%
  group_by(regen_type) %>%
  summarise(n = n(),
            across(all_of(summary_vars), ~ round(mean(.x, na.rm = TRUE), 3)),
            .groups = "drop") %>%
  arrange(regen_type)

cat("=== ARCHETYPE PROFILES ===\n\n")
print(as.data.frame(
  archetype_profiles %>%
    select(regen_type, n, dens_ha, mean_hgt, spruce_share,
           sp_richness, CWM_drought,
           frac_empty_sub, disturbance_severity, n_mature)
))

cat("\n=== Does residual overstory align with a type? ===\n")
print(table(regen_type = merged_df_cluster$regen_type,
            has_mature = merged_df_cluster$has_mature))
# If has_mature concentrates in one type -> that's your residual-overstory
# archetype, and its medoid's tree file will carry the survivors.

# ---------------------------------------------------------------------
# 5b. Driver check: does elevation structure the response? (post-hoc)
# ---------------------------------------------------------------------
driver_long <- merged_df_cluster %>%
  select(regen_type, elevation,
         effective_numbers, sp_richness, spruce_share) %>%
  pivot_longer(c(effective_numbers, sp_richness, spruce_share),
               names_to = "response", values_to = "value")

print(
  ggplot(driver_long, aes(elevation, value)) +
    geom_point(aes(colour = regen_type), alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, colour = "black") +
    facet_wrap(~ response, scales = "free_y") +
    scale_colour_brewer(palette = "Set2") +
    labs(title = "Elevation vs regeneration response (plot level)",
         x = "elevation (m)", y = NULL, colour = "type") +
    theme_bw()
)

cat("\n=== Elevation ~ response (linear fits) ===\n")
for (v in c("effective_numbers", "sp_richness", "spruce_share")) {
  m <- lm(reformulate("elevation", v), data = merged_df_cluster)
  cat(sprintf("  %-17s slope=%+.4f  R2=%.3f  p=%.3g\n",
              v, coef(m)[2], summary(m)$r.squared,
              summary(m)$coefficients[2, 4]))
}
cat("\n")


# see the 3 regeneration groups
merged_df_cluster %>%
  filter(regen_type %in% c("T1","T2","T3")) %>%          # drop T0_bare/T_mature: too few/degenerate for a trend line
  select(regen_type, elevation, effective_numbers, sp_richness, spruce_share) %>%
  pivot_longer(c(effective_numbers, sp_richness, spruce_share),
               names_to = "response", values_to = "value") %>%
  ggplot(aes(elevation, value, colour = regen_type)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ response, scales = "free_y") +
  scale_colour_brewer(palette = "Set2") +
  labs(title = "Elevation vs regeneration response, by archetype",
       x = "elevation (m)", y = NULL, colour = "type") +
  theme_bw()



# =====================================================================
# 7. iLand initialization  —  extract each medoid's REAL stems
# =====================================================================
# For each archetype we take its medoid plot and read that plot's actual
# species x height/dbh-class stems from sub_veg. No species averaging:
# the archetype inherits the observed assemblage of one real plot.

# ---- 7a. tuning knobs (CHECK) ---------------------------------------
SUBPLOT_AREA_M2     <- 4       # quadrat area per subplot (verify!)
ILAND_STAND_AREA_M2 <- 10000   # target stand area: 10000 = 1 ha RU,
# 100 = one 10 m stand pixel

# ---- 7b. species-code crosswalk (CHECK / EDIT) ----------------------
# Most 4-letter codes are already valid iLand species. The genus-level
# "*sp" codes and "ots1" are NOT valid iLand species and MUST be mapped.
species_map <- c(
  piab = "piab", fasy = "fasy", pisy = "pisy", lade = "lade",
  psme = "psme", acps = "acps", acpl = "acpl", potr = "potr",
  cabe = "cabe", prav = "prav", soau = "soau", soar = "soar",
  algl = "algl", alin = "alin", saca = "saca",
  # --- genus-level field codes: assign a real iLand species ---
  besp = "bepe",   # Betula sp.  -> B. pendula
  qusp = "quro",   # Quercus sp. -> Q. robur   (or "qupe")
  sasp = "saca",   # Salix sp.   -> S. caprea
  absp = "abal",   # Abies sp.   -> A. alba
  tisp = "tico",   # Tilia sp.   -> T. cordata (or "tipl")
  ots1 = NA_character_   # "other": decide (drop or assign)
)

# ---- 7c. height-class -> iLand (from,to); >4 capped to the 2-4 bin ---
norm_dash <- function(x) gsub("[\u2012\u2013\u2014-]", "-", x)  # any dash -> "-"
hgt_map <- data.frame(
  key         = c("0.2-0.4","0.4-0.6","0.6-0.8","0.8-1.0","1.0-1.3","1.3-2.0","2-4",">4"),
  height_from = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.3, 2.0, 2.0),
  height_to   = c(0.4, 0.6, 0.8, 1.0, 1.3, 2.0, 4.0, 4.0),   # >4 -> 2-4 (cap)
  stringsAsFactors = FALSE
)

# ---- 7d. builders ----------------------------------------------------
# saplings: regeneration < 4 m (vegtype small + advanced)
build_saplings <- function(plot_key, stand_id) {
  d <- sub_veg %>%
    filter(plot == plot_key, !is.na(n), n > 0, vegtype %in% c("small","advanced")) %>%
    mutate(key = norm_dash(hgt)) %>%
    left_join(hgt_map, by = "key")
  
  if (any(is.na(d$height_from)))
    warning("Unmapped height class in plot ", plot_key, ": ",
            paste(unique(d$hgt[is.na(d$height_from)]), collapse = ", "))
  
  d %>%
    mutate(n_sub = n_distinct(subplot)) %>%
    group_by(species, height_from, height_to) %>%
    summarise(stems = sum(n), n_sub = first(n_sub), .groups = "drop") %>%
    mutate(
      iland_species = species_map[species],
      count = round(stems / (n_sub * SUBPLOT_AREA_M2) * ILAND_STAND_AREA_M2),
      stand_id = stand_id
    ) %>%
    filter(!is.na(iland_species), count > 0) %>%
    group_by(stand_id, species = iland_species, height_from, height_to) %>%
    summarise(count = sum(count), .groups = "drop") %>%   # merge mapped duplicates
    select(stand_id, species, count, height_from, height_to)
}

# trees: mature, dbh >= 10 cm (only where present)
build_trees <- function(plot_key, stand_id) {
  d <- sub_veg %>%
    filter(plot == plot_key, !is.na(n), n > 0, vegtype == "mature")
  if (nrow(d) == 0) return(NULL)
  
  d %>%
    mutate(
      dbh_from = case_when(dbh_est < 20 ~ 10, dbh_est < 40 ~ 20,
                           dbh_est < 60 ~ 40, TRUE ~ 60),
      dbh_to   = case_when(dbh_est < 20 ~ 20, dbh_est < 40 ~ 40,
                           dbh_est < 60 ~ 60, TRUE ~ 80),
      hd  = round(hgt_est * 100 / dbh_est),   # height(cm)/dbh(cm)
      age = 0,                                # 0 -> iLand estimates from dbh
      n_sub = n_distinct(subplot),
      iland_species = species_map[species],
      count = round(n / (n_sub * SUBPLOT_AREA_M2) * ILAND_STAND_AREA_M2),
      stand_id = stand_id
    ) %>%
    filter(!is.na(iland_species), count > 0) %>%
    group_by(stand_id, species = iland_species, dbh_from, dbh_to, hd, age) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    select(stand_id, species, count, dbh_from, dbh_to, hd, age)
}

# ---- 7e. run over archetypes ----------------------------------------
# stand_id is a placeholder here (1,2,3...); replace with the actual
# stand-grid value each archetype occupies in your landscape.
# add for bare soil and for mature-Trees present
bare_plot <- merged_df_cluster %>%
  filter(regen_type == "T0_bare") %>%
  pull(plot)

mature_data <- merged_df_cluster %>%
  filter(regen_type == "T_mature") %>%
  select(plot, all_of(vars)) %>%     # vars = dens_ha, mean_hgt, spruce_share, sp_richness, frac_empty_sub
  drop_na()

mature_dist <- as.matrix(dist(scale(mature_data %>% select(-plot))))
mature_medoid_plot <- mature_data$plot[which.min(rowSums(mature_dist))]

# add other groups into medoids
medoids <- bind_rows(
  medoids %>% mutate(regen_type = paste0("T", cluster_raw)),
  tibble(cluster_raw = NA_integer_, medoid_plot = bare_plot,          regen_type = "T0_bare"),
  tibble(cluster_raw = NA_integer_, medoid_plot = mature_medoid_plot, regen_type = "T_mature")
) %>%
  mutate(stand_id = row_number())

medoids

out_dir <- "outData/iland_init"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (i in seq_len(nrow(medoids))) {
  pk  <- medoids$medoid_plot[i]
  sid <- medoids$stand_id[i]
  typ <- medoids$regen_type[i]          # <- also fixes the TNA bug from before
  
  sap <- build_saplings(pk, sid)
  tre <- build_trees(pk, sid)
  
  cat(sprintf("\n--- %s | medoid plot %s | stand_id %d ---\n", typ, pk, sid))
  cat("saplings:\n"); print(as.data.frame(sap))
  write.csv(sap, file.path(out_dir, sprintf("iniplt_%s_sap.txt", typ)),
            row.names = FALSE, quote = TRUE)
  
  if (!is.null(tre)) {
    cat("mature trees:\n"); print(as.data.frame(tre))
    write.csv(tre, file.path(out_dir, sprintf("iniplt_%s_tree.txt", typ)),
              row.names = FALSE, quote = TRUE)
  } else {
    cat("mature trees: none (empty tree layer)\n")
  }
}

cat("\nT0_bare -> initialize as bare soil (no init file).\n")

# ---------------------------------------------------------------------
# 8. Export
# ---------------------------------------------------------------------
write.csv(merged_df_cluster,  "outData/CZ_plots_archetypes.csv",   row.names = FALSE)
write.csv(archetype_profiles, "outData/CZ_archetype_profiles.csv", row.names = FALSE)
write.csv(medoids,            "outData/CZ_archetype_medoids.csv",  row.names = FALSE)

cat("\nDone. iLand init files in ./iland_init/ ; one sapling (+ optional tree)\n",
    "file per archetype, each drawn from its medoid's real stems.\n")
