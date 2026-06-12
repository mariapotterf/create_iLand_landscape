# ══════════════════════════════════════════════════════════
#  POST-BROUK REGENERATION CLUSTERING
#  Goal: group 208 CZ field plots into ecologically
#  meaningful regeneration archetypes for iLand initialization
#
#  Input:  env_chars_subplots.gpkg
#  Output: cluster assignments + archetype profiles
#
#  Logic:
#  1. Load & filter Czech Republic, most recent survey year
#  2. Aggregate subplots → plot level
#  3. Derive climate variables from ERA5 (optional, for validation)
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

# ── 0. Load data ─────────────────────────────────────────

gdf <- st_read("env_chars_subplots.gpkg", quiet = TRUE)

# Filter Czech Republic only
cz <- gdf %>%
  filter(country_name == "Czech Republic") %>%
  # Keep most recent survey year per subplot
  group_by(subplot) %>%
  slice_max(year, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("CZ subplots:", nrow(cz), "\n")
cat("CZ plots:   ", n_distinct(cz$plot_id), "\n")
cat("Years:      ", paste(sort(unique(cz$year)), collapse = ", "), "\n\n")

# ── 1. Aggregate subplots → plot level ───────────────────
#
#  For each plot compute:
#  STRUCTURE:   stem density, height, height CV
#  COMPOSITION: spruce share, conifer IV, deciduous IV,
#               species richness, CWM drought tolerance
#  STATUS:      fraction of subplots with zero regeneration
#  MANAGEMENT:  planting intensity, browsing protection
#  SITE:        distance to edge, disturbance severity,
#               elevation (from field data)

plots <- cz %>%
  st_drop_geometry() %>%
  group_by(plot_id) %>%
  summarise(
    # ── survey info
    year             = first(year),
    n_subplots       = n(),
    disturbance_year = first(disturbance_year),
    time_since_dist  = year - disturbance_year,
    
    # ── structure
    mean_dens_ha     = mean(dens_ha,    na.rm = TRUE),
    sd_dens_ha       = sd(dens_ha,      na.rm = TRUE),
    mean_hgt         = mean(mean_hgt,   na.rm = TRUE),
    mean_cv_hgt      = mean(cv_hgt,     na.rm = TRUE),
    
    # ── composition
    mean_spruce      = mean(spruce_share,   na.rm = TRUE),
    mean_IV_conif    = mean(IV_coniferous,  na.rm = TRUE),
    mean_IV_decid    = mean(IV_deciduous,   na.rm = TRUE),
    mean_sp_rich     = mean(sp_richness,    na.rm = TRUE),
    mean_CWM_drought = mean(CWM_drought,    na.rm = TRUE),
    
    # ── regeneration status (subplot level)
    n_empty_sub      = sum(active_regeneration == 0, na.rm = TRUE),
    frac_empty_sub   = mean(active_regeneration == 0, na.rm = TRUE),
    all_empty        = all(active_regeneration == 0, na.rm = TRUE),
    any_empty        = any(active_regeneration == 0, na.rm = TRUE),
    
    # ── management
    planting         = mean(planting_intensity,       na.rm = TRUE),
    browsing_prot    = mean(anti_browsing_intensity,  na.rm = TRUE),
    
    # ── site
    elevation        = mean(elevation,           na.rm = TRUE),
    dist_to_edge     = mean(dist_to_edge,        na.rm = TRUE),
    disturbance_sev  = mean(disturbance_severity, na.rm = TRUE),
    
    .groups = "drop"
  )

cat("Plot-level dataset:", nrow(plots), "plots\n\n")

# ── 2. Print emptiness summary ────────────────────────────

cat("=== Regeneration status ===\n")
cat("Plots with ALL subplots empty (truly bare):",
    sum(plots$all_empty), sprintf("(%.1f%%)\n", 100*mean(plots$all_empty)))
cat("Plots with ANY subplot empty:",
    sum(plots$any_empty), sprintf("(%.1f%%)\n", 100*mean(plots$any_empty)))
cat("Plots with NO empty subplots:",
    sum(!plots$any_empty), sprintf("(%.1f%%)\n\n", 100*mean(!plots$any_empty)))

cat("Empty subplot distribution (0=none, 5=half, 10=all):\n")
print(table(plots$n_empty_sub))
cat("\n")

# ── 3. Rule-based pre-classification ─────────────────────
#
#  Before clustering, assign rule-based categories that
#  are ecologically unambiguous and should NOT be split
#  by a clustering algorithm:
#
#  TYPE 0 — BARE/DELAYED:
#    All subplots empty (active_regeneration == 0 everywhere)
#    → iLand: NBS initialization, very low/zero stem input
#
#  All remaining plots go into data-driven clustering below.

plots <- plots %>%
  mutate(
    regen_type = case_when(
      all_empty ~ "T0_bare",
      TRUE      ~ "to_cluster"
    )
  )

cat("Pre-classified as T0_bare:", sum(plots$regen_type == "T0_bare"), "\n")
cat("Remaining for clustering: ", sum(plots$regen_type == "to_cluster"), "\n\n")

# ── 4. Clustering of non-bare plots ──────────────────────

clust_data <- plots %>%
  filter(regen_type == "to_cluster") %>%
  # Variables that matter for iLand initialization:
  # what species are present, how dense, how drought-adapted
  # NOT: management (explains HOW it got there, not WHAT is there)
  # NOT: site vars (iLand handles those internally)
  select(plot_id,
         mean_dens_ha,    # how many stems to initialize
         mean_hgt,        # structural stage
         mean_spruce,     # spruce dominance
         mean_IV_conif,   # conifer community
         mean_IV_decid,   # deciduous community
         mean_sp_rich,    # diversity
         mean_CWM_drought # drought tolerance of community
  ) %>%
  drop_na()

cat("Plots entering clustering:", nrow(clust_data), "\n\n")

# Scale
X <- clust_data %>% select(-plot_id)
X_sc <- scale(X)

# ── 4a. Determine optimal k ──────────────────────────────

set.seed(42)
sil_scores <- sapply(2:8, function(k) {
  km  <- kmeans(X_sc, centers = k, nstart = 50, iter.max = 100)
  sil <- silhouette(km$cluster, dist(X_sc))
  mean(sil[, 3])
})

cat("Silhouette scores:\n")
for (i in seq_along(sil_scores)) {
  cat(sprintf("  k=%d: %.3f\n", i+1, sil_scores[i]))
}
k_opt <- which.max(sil_scores) + 1
cat(sprintf("\nOptimal k = %d\n\n", k_opt))

# Visualize silhouette
sil_df <- data.frame(k = 2:8, silhouette = sil_scores)
p_sil <- ggplot(sil_df, aes(x = k, y = silhouette)) +
  geom_line() + geom_point(size = 3) +
  geom_vline(xintercept = k_opt, linetype = "dashed", color = "red") +
  labs(title = "Silhouette score by number of clusters",
       x = "Number of clusters (k)", y = "Mean silhouette width") +
  theme_bw()
print(p_sil)

# ── 4b. Final clustering ──────────────────────────────────

set.seed(42)
km_final <- kmeans(X_sc, centers = k_opt, nstart = 50, iter.max = 100)
clust_data$cluster_raw <- km_final$cluster

# Merge back
plots <- plots %>%
  left_join(clust_data %>% select(plot_id, cluster_raw),
            by = "plot_id") %>%
  mutate(
    regen_type = case_when(
      regen_type == "T0_bare" ~ "T0_bare",
      TRUE ~ paste0("T", cluster_raw)
    )
  )

cat("Final type distribution:\n")
print(table(plots$regen_type))
cat("\n")

# ── 5. Describe archetypes ────────────────────────────────

summary_vars <- c("mean_dens_ha", "mean_hgt", "mean_cv_hgt",
                  "mean_spruce", "mean_IV_conif", "mean_IV_decid",
                  "mean_sp_rich", "mean_CWM_drought",
                  "frac_empty_sub", "planting", "browsing_prot",
                  "elevation", "dist_to_edge", "time_since_dist")

archetype_profiles <- plots %>%
  group_by(regen_type) %>%
  summarise(
    n             = n(),
    across(all_of(summary_vars), ~round(mean(.x, na.rm=TRUE), 3)),
    .groups = "drop"
  ) %>%
  arrange(regen_type)

cat("=== ARCHETYPE PROFILES ===\n\n")
print(archetype_profiles %>%
        select(regen_type, n, mean_dens_ha, mean_hgt,
               mean_spruce, mean_IV_conif, mean_IV_decid,
               mean_sp_rich, mean_CWM_drought,
               frac_empty_sub, planting) %>%
        as.data.frame())

# ── 6. Visualize archetypes ───────────────────────────────

# Radar-style bar chart of key variables
key_vars <- c("mean_dens_ha", "mean_spruce", "mean_IV_conif",
              "mean_IV_decid", "mean_CWM_drought", "frac_empty_sub",
              "planting", "mean_hgt")

plot_df <- archetype_profiles %>%
  select(regen_type, n, all_of(key_vars)) %>%
  # Normalize each variable 0-1 for display
  mutate(across(all_of(key_vars), ~(.x - min(.x)) / (max(.x) - min(.x) + 1e-9))) %>%
  pivot_longer(cols = all_of(key_vars), names_to = "variable", values_to = "value") %>%
  mutate(label = paste0(regen_type, "\n(n=", n, ")"))

p_bar <- ggplot(plot_df, aes(x = variable, y = value, fill = regen_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Regeneration archetype profiles (normalized)",
       subtitle = "CZ post-disturbance plots — for iLand initialization",
       x = NULL, y = "Normalized value (0-1)", fill = "Type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
print(p_bar)

# Map of cluster assignments
plots_sf <- cz %>%
  select(plot_id, geometry) %>%
  distinct(plot_id, .keep_all = TRUE) %>%
  left_join(plots %>% select(plot_id, regen_type), by = "plot_id")

p_map <- ggplot(plots_sf) +
  geom_sf(aes(color = regen_type), size = 2.5) +
  scale_color_brewer(palette = "Set2") +
  labs(title = "Spatial distribution of regeneration archetypes",
       color = "Type") +
  theme_bw()
print(p_map)

# ── 7. Export ─────────────────────────────────────────────

write.csv(plots, "CZ_plots_archetypes.csv", row.names = FALSE)
write.csv(archetype_profiles, "CZ_archetype_profiles.csv", row.names = FALSE)

cat("\n=== iLand initialization summary ===\n")
cat("Use archetype_profiles.csv as input to iLand stand initialization.\n")
cat("Key fields per archetype:\n")
cat("  mean_dens_ha  → initial stem density\n")
cat("  mean_hgt      → initial mean height\n")
cat("  mean_spruce   → Picea abies fraction\n")
cat("  mean_IV_decid → deciduous fraction\n")
cat("  mean_CWM_drought → drought tolerance of community\n")
cat("  T0_bare       → initialize as bare soil (NBS)\n")