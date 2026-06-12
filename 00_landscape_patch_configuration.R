library(terra)

# ══════════════════════════════════════════════════════════
#  PARAMETERS
# ══════════════════════════════════════════════════════════

# Grid dimensions
n      <- 100   # inner grid: 100×100 @ 10 m = 1 km²
belt   <- 20    # belt: 20 cells = 200 m
res_m  <- 10    # stand grid resolution (m)
res_ru <- 100   # RU grid resolution (m)

# C4 2022: 38.6% disturbed (Table D2: 114,550 / 296,860 ha)
target_frac <- 0.386
target_size <- round(target_frac * n * n)  # 3860 cells = 3.86 ha

# Coordinates — from real Kostelec landscape files
# (standgrid_20260330.asc and ru_ids_20260330.asc)
# CRS: Křovák EPSG:5514
XLL <- -823100
YLL <- -1178600

# iLand stand grid values — confirmed from standgrid_20260330.asc:
#   0   = NODATA (truly outside, non-forested) — written as NODATA_value
#  -2   = forested outside belt (provides seed rain into disturbed area)
#  -1   = non-forested outside (no trees, no seeds) — not used here
#   1   = disturbed stand (post-bark-beetle, C4 2022)
#   2   = mature forest stand
VAL_DISTURBED      <-  1L
VAL_MATURE         <-  2L
VAL_BELT           <- -2L   # forested outside — enables seed rain
VAL_NODATA_STAND   <-  0L   # NODATA for stand grid

# RU grid values — confirmed from ru_ids_20260330.asc:
#   0      = outside (NODATA)
#   1,2,.. = sequential RU IDs for active cells
VAL_NODATA_RU      <-  0L

# ══════════════════════════════════════════════════════════
#  HELPER: add forested belt around inner grid
#  Belt = -2 (forested outside) matching real landscape files
# ══════════════════════════════════════════════════════════
add_belt <- function(grid) {
  new_n <- nrow(grid) + 2 * belt
  out   <- matrix(VAL_BELT, nrow = new_n, ncol = new_n)
  out[(belt+1):(belt+n), (belt+1):(belt+n)] <- grid
  out
}

# ══════════════════════════════════════════════════════════
#  HELPER: matrix → terra raster
#  Origin = Kostelec coordinates, CRS = Křovák EPSG:5514
# ══════════════════════════════════════════════════════════
to_raster <- function(mat, cell_res = res_m) {
  new_n    <- nrow(mat)
  extent_m <- new_n * cell_res
  r <- rast(nrows = new_n, ncols = new_n,
            xmin = XLL,
            xmax = XLL + extent_m,
            ymin = YLL,
            ymax = YLL + extent_m)
  values(r) <- as.vector(t(mat))
  crs(r)    <- "EPSG:5514"
  r
}

# ══════════════════════════════════════════════════════════
#  CORE: grow n_blobs disturbed patches
#
#  Algorithm:
#  1. Place n_blobs random seed cells (one per patch)
#  2. Maintain frontier of candidate neighbour cells per blob
#  3. At each step: pick eligible blob, pick random frontier
#     cell, add it, expand frontier
#  4. Stop when total disturbed cells == target_size
# ══════════════════════════════════════════════════════════
make_landscape <- function(n_blobs, seed = 42) {
  set.seed(seed)
  
  grid <- matrix(VAL_MATURE, nrow = n, ncol = n)
  dirs <- matrix(c(1,0,-1,0,0,1,0,-1), ncol = 2, byrow = TRUE)
  
  # Equal cells per blob; remainder goes to first blob
  blob_targets    <- rep(floor(target_size / n_blobs), n_blobs)
  blob_targets[1] <- blob_targets[1] + (target_size - sum(blob_targets))
  blob_counts     <- integer(n_blobs)
  
  # Random seed placement — avoid 5-cell margin from edge
  margin    <- 5L
  all_cells <- expand.grid(r = (margin+1):(n-margin),
                           c = (margin+1):(n-margin))
  seed_idx  <- sample(nrow(all_cells), n_blobs, replace = FALSE)
  
  # Initialise blobs and frontiers
  frontiers <- vector("list", n_blobs)
  for (i in seq_len(n_blobs)) {
    r  <- all_cells$r[seed_idx[i]]
    cc <- all_cells$c[seed_idx[i]]
    grid[r, cc]    <- VAL_DISTURBED
    blob_counts[i] <- 1L
    nbrs <- sweep(dirs, 2, c(r, cc), "+")
    nbrs <- nbrs[nbrs[,1] >= 1 & nbrs[,1] <= n &
                   nbrs[,2] >= 1 & nbrs[,2] <= n, , drop = FALSE]
    frontiers[[i]] <- nbrs
  }
  
  total_disturbed <- n_blobs
  
  while (total_disturbed < target_size) {
    
    remaining <- blob_targets - blob_counts
    eligible  <- which(remaining > 0 &
                         sapply(frontiers, nrow) > 0)
    
    if (length(eligible) == 0) {
      warning(sprintf(
        "Frontiers exhausted: placed %d of %d cells.",
        total_disturbed, target_size))
      break
    }
    
    b     <- eligible[sample(length(eligible), 1)]
    front <- frontiers[[b]]
    idx   <- sample(nrow(front), 1)
    r     <- front[idx, 1]
    cc    <- front[idx, 2]
    
    frontiers[[b]] <- front[-idx, , drop = FALSE]
    
    if (grid[r, cc] == VAL_MATURE) {
      grid[r, cc]     <- VAL_DISTURBED
      blob_counts[b]  <- blob_counts[b] + 1L
      total_disturbed <- total_disturbed + 1L
      
      nbrs <- sweep(dirs, 2, c(r, cc), "+")
      nbrs <- nbrs[nbrs[,1] >= 1 & nbrs[,1] <= n &
                     nbrs[,2] >= 1 & nbrs[,2] <= n, , drop = FALSE]
      mature_nbrs <- nbrs[apply(nbrs, 1, function(x)
        grid[x[1], x[2]] == VAL_MATURE), , drop = FALSE]
      if (nrow(mature_nbrs) > 0)
        frontiers[[b]] <- rbind(frontiers[[b]], mature_nbrs)
    }
  }
  
  n_dist  <- sum(grid == VAL_DISTURBED)
  blob_ha <- (n_dist / n_blobs) * res_m^2 / 10000
  cat(sprintf(
    "[%d blob(s)] disturbed: %d cells = %.2f ha (%.1f%%) | %.3f ha/blob | sizes: %s\n",
    n_blobs, n_dist, n_dist * res_m^2 / 10000,
    100 * n_dist / (n * n), blob_ha,
    paste(blob_counts, collapse = "+")))
  
  to_raster(add_belt(grid))
}

# ══════════════════════════════════════════════════════════
#  RESOURCE UNIT GRID
#  100 m resolution — same origin and CRS as stand grid
#  Values: 0 = outside (NODATA), 1..n = sequential RU IDs
#  Belt = 2 cells @ 100 m = 200 m (matches stand grid belt)
# ══════════════════════════════════════════════════════════
make_ru_grid <- function() {
  new_n_10m <- n + 2 * belt            # 140 cells @ 10 m
  extent_m  <- new_n_10m * res_m       # 1400 m
  n_ru      <- extent_m / res_ru       # 14 cells @ 100 m
  belt_ru   <- (belt * res_m) / res_ru # 2 cells
  
  g  <- matrix(VAL_NODATA_RU, nrow = n_ru, ncol = n_ru)
  id <- 1L
  for (i in (belt_ru+1):(n_ru-belt_ru)) {
    for (j in (belt_ru+1):(n_ru-belt_ru)) {
      g[i, j] <- id
      id <- id + 1L
    }
  }
  
  n_active <- id - 1L
  cat(sprintf(
    "RU grid: %d×%d cells @ %d m | %d active RUs | belt = %d cells\n",
    n_ru, n_ru, res_ru, n_active, belt_ru))
  
  to_raster(g, cell_res = res_ru)
}

# ══════════════════════════════════════════════════════════
#  GENERATE LANDSCAPES
# ══════════════════════════════════════════════════════════

cat("=== iLand landscape generator — C4 2022 post-disturbance ===\n")
cat(sprintf("Inner grid:  %d×%d cells @ %d m = %d×%d m\n",
            n, n, res_m, n*res_m, n*res_m))
cat(sprintf("With belt:   %d×%d cells @ %d m = %d×%d m\n",
            n+2*belt, n+2*belt, res_m,
            (n+2*belt)*res_m, (n+2*belt)*res_m))
cat(sprintf("Origin:      xll=%d, yll=%d (Křovák EPSG:5514)\n",
            XLL, YLL))
cat(sprintf("Disturbed:   %.1f%% = %d cells = %.2f ha [C4 2022]\n\n",
            100*target_frac, target_size,
            target_size*res_m^2/10000))

r1   <- make_landscape(n_blobs = 1, seed = 42)
r2   <- make_landscape(n_blobs = 2, seed = 42)
r4   <- make_landscape(n_blobs = 4, seed = 42)
r8   <- make_landscape(n_blobs = 8, seed = 42)
r_ru <- make_ru_grid()

# ══════════════════════════════════════════════════════════
#  WRITE ASCII FILES
#  NODATA_value in header will be set by terra to match
#  the minimum value in the raster — force it explicitly
# ══════════════════════════════════════════════════════════

writeRaster(r1,   "outData/stand_C4_2022_1blob.asc",
            filetype = "AAIGrid", NAflag = 0, overwrite = TRUE)
writeRaster(r2,   "outData/stand_C4_2022_2blobs.asc",
            filetype = "AAIGrid", NAflag = 0, overwrite = TRUE)
writeRaster(r4,   "outData/stand_C4_2022_4blobs.asc",
            filetype = "AAIGrid", NAflag = 0, overwrite = TRUE)
writeRaster(r8,   "outData/stand_C4_2022_8blobs.asc",
            filetype = "AAIGrid", NAflag = 0, overwrite = TRUE)
writeRaster(r_ru, "outData/ru_grid.asc",
            filetype = "AAIGrid", NAflag = 0, overwrite = TRUE)

cat("\n=== Files written to outData/ ===\n")
cat("Stand grids: NODATA=0 | -2=forested belt | 1=disturbed | 2=mature\n")
cat("RU grid:     NODATA=0 | 1..100=active RU IDs\n")

# ══════════════════════════════════════════════════════════
#  SANITY CHECKS
# ══════════════════════════════════════════════════════════

cat("\n=== Sanity checks (r1) ===\n")
cat(sprintf("Dimensions: %d rows × %d cols (expect 140×140)\n",
            nrow(r1), ncol(r1)))
cat(sprintf("Extent: xmin=%.0f xmax=%.0f ymin=%.0f ymax=%.0f\n",
            xmin(r1), xmax(r1), ymin(r1), ymax(r1)))
cat(sprintf("CRS: %s\n", crs(r1, describe=TRUE)$code))
cat(sprintf("Belt value (cell 1,1): %d (expect -2)\n",
            as.integer(r1[1,1])))
cat(sprintf("Inner cell (21,21): %d (expect 1 or 2)\n",
            as.integer(r1[21,21])))
cat(sprintf("Unique values: %s\n",
            paste(sort(unique(values(r1))), collapse=", ")))

cat("\n=== Sanity checks (r_ru) ===\n")
cat(sprintf("Dimensions: %d rows × %d cols (expect 14×14)\n",
            nrow(r_ru), ncol(r_ru)))
cat(sprintf("Active RUs (>0): %d (expect 100)\n",
            sum(values(r_ru) > 0)))
cat(sprintf("Max RU ID: %d (expect 100)\n",
            max(values(r_ru))))

# ══════════════════════════════════════════════════════════
#  VISUAL CHECK
# ══════════════════════════════════════════════════════════

cols   <- c("gray30", "pink", "darkgreen")
breaks <- c(-2.5, 0.5, 1.5, 2.5)

# 5 panels: 4 stand grids + 1 legend
layout(matrix(1:5, nrow = 1), widths = c(4, 4, 4, 4, 1.5))
par(mar = c(2, 2, 3, 1))

plot(r1, col = cols, breaks = breaks, legend = FALSE,
     main = "1 blob\n3.86 ha")
plot(r2, col = cols, breaks = breaks, legend = FALSE,
     main = "2 blobs\n1.93 ha each")
plot(r4, col = cols, breaks = breaks, legend = FALSE,
     main = "4 blobs\n0.97 ha each")
plot(r8, col = cols, breaks = breaks, legend = FALSE,
     main = "8 blobs\n0.48 ha each")

par(mar = c(2, 0, 3, 1))
plot.new()
legend("center",
       legend = c("belt (-2)", "disturbed (1)", "mature (2)"),
       fill   = cols,
       title  = "Stand type",
       cex    = 0.9,
       bty    = "n")
