library(terra)

# ══════════════════════════════════════════════════════════
#  PARAMETERS
# ══════════════════════════════════════════════════════════

n      <- 100   # inner grid: 100×100 @ 10 m = 1 km²
belt   <- 20    # belt: 20 cells = 200 m
res_m  <- 10    # resolution (m)
res_ru <- 100   # RU resolution (m)

# C4 2022: 38.6% disturbed (Table D2: 114,550 / 296,860 ha)
target_frac <- 0.386
target_size <- round(target_frac * n * n)  # 3860 cells = 3.86 ha

VAL_DISTURBED <- 1L
VAL_MATURE    <- 2L
VAL_EDGE      <- -1L
VAL_RU_FOREST <- 1L
VAL_RU_EDGE   <- 0L

# ══════════════════════════════════════════════════════════
#  HELPER: add belt
# ══════════════════════════════════════════════════════════
add_belt <- function(grid, belt, edge_val) {
  new_n <- nrow(grid) + 2 * belt
  out   <- matrix(edge_val, nrow = new_n, ncol = new_n)
  out[(belt+1):(belt+n), (belt+1):(belt+n)] <- grid
  out
}

# ══════════════════════════════════════════════════════════
#  HELPER: matrix → terra raster
# ══════════════════════════════════════════════════════════
to_raster <- function(mat , crs_str = "EPSG:25833" # crs(r) <- "EPSG:25833"  # ETRS89 / UTM zone 33N — Central Europe standard
                      ) {
  new_n <- nrow(mat)
  r <- rast(nrows = new_n, ncols = new_n,
            xmin = 0, xmax = new_n * res_m,
            ymin = 0, ymax = new_n * res_m)
  values(r) <- as.vector(t(mat))
  crs(r)    <- crs_str
  r
}

# ══════════════════════════════════════════════════════════
#  CORE: grow n_blobs patches until total = target_size
#
#  1. Place n_blobs random seed cells (one per patch)
#  2. Maintain a frontier of candidate neighbour cells
#     for each patch
#  3. At each step: pick a random patch (weighted by how
#     many cells it still needs), pick a random cell from
#     its frontier, add it, expand frontier
#  4. Stop when total disturbed cells == target_size
# ══════════════════════════════════════════════════════════
make_landscape <- function(n_blobs, seed = 42) {
  set.seed(seed)
  
  grid  <- matrix(VAL_MATURE, nrow = n, ncol = n)
  dirs  <- matrix(c(1,0,-1,0,0,1,0,-1), ncol = 2, byrow = TRUE)
  
  # Equal target per blob; remainder to first
  blob_targets        <- rep(floor(target_size / n_blobs), n_blobs)
  blob_targets[1]     <- blob_targets[1] + (target_size - sum(blob_targets))
  blob_counts         <- integer(n_blobs)  # current cell count per blob
  
  # Random seed cells — sample without replacement from inner grid
  margin    <- 5L
  all_cells <- expand.grid(r = (margin+1):(n-margin),
                           c = (margin+1):(n-margin))
  seed_idx  <- sample(nrow(all_cells), n_blobs, replace = FALSE)
  
  # Frontier: list of candidate neighbour cells per blob
  # stored as list of 2-col matrices
  frontiers <- vector("list", n_blobs)
  
  for (i in seq_len(n_blobs)) {
    r <- all_cells$r[seed_idx[i]]
    cc <- all_cells$c[seed_idx[i]]
    grid[r, cc]    <- VAL_DISTURBED
    blob_counts[i] <- 1L
    # Initial frontier: valid neighbours of seed
    nbrs <- sweep(dirs, 2, c(r, cc), "+")
    nbrs <- nbrs[nbrs[,1] >= 1 & nbrs[,1] <= n &
                   nbrs[,2] >= 1 & nbrs[,2] <= n, , drop = FALSE]
    frontiers[[i]] <- nbrs
  }
  
  total_disturbed <- n_blobs  # seeds already placed
  
  while (total_disturbed < target_size) {
    
    # Which blobs still need more cells?
    remaining <- blob_targets - blob_counts
    eligible  <- which(remaining > 0 & 
                         sapply(frontiers, nrow) > 0)
    
    if (length(eligible) == 0) {
      warning("All frontiers exhausted before reaching target_size. ",
              "Placed ", total_disturbed, " of ", target_size, " cells.")
      break
    }
    
    # Pick a blob randomly (uniform among eligible)
    b <- eligible[sample(length(eligible), 1)]
    
    # Pick a random cell from its frontier
    front <- frontiers[[b]]
    idx   <- sample(nrow(front), 1)
    r     <- front[idx, 1]
    cc    <- front[idx, 2]
    
    # Remove from frontier regardless (avoids duplicate checks)
    frontiers[[b]] <- front[-idx, , drop = FALSE]
    
    # Only add if still mature (another blob may have taken it)
    if (grid[r, cc] == VAL_MATURE) {
      grid[r, cc]    <- VAL_DISTURBED
      blob_counts[b] <- blob_counts[b] + 1L
      total_disturbed <- total_disturbed + 1L
      
      # Expand frontier with new neighbours
      nbrs <- sweep(dirs, 2, c(r, cc), "+")
      nbrs <- nbrs[nbrs[,1] >= 1 & nbrs[,1] <= n &
                     nbrs[,2] >= 1 & nbrs[,2] <= n, , drop = FALSE]
      # Only add cells that are still mature
      mature_nbrs <- nbrs[apply(nbrs, 1, function(x)
        grid[x[1], x[2]] == VAL_MATURE), , drop = FALSE]
      if (nrow(mature_nbrs) > 0)
        frontiers[[b]] <- rbind(frontiers[[b]], mature_nbrs)
    }
  }
  
  cat(sprintf(
    "[%d blob(s)] disturbed: %d cells = %.2f ha (%.1f%%) | %.3f ha/blob | blob sizes: %s\n",
    n_blobs,
    sum(grid == VAL_DISTURBED),
    sum(grid == VAL_DISTURBED) * res_m^2 / 10000,
    100 * sum(grid == VAL_DISTURBED) / (n * n),
    (sum(grid == VAL_DISTURBED) / n_blobs) * res_m^2 / 10000,
    paste(blob_counts, collapse = "+")))
  
  to_raster(add_belt(grid, belt, VAL_EDGE))
}

# ══════════════════════════════════════════════════════════
#  RESOURCE UNIT GRID
# ══════════════════════════════════════════════════════════
make_ru_grid <- function() {
  new_n    <- n + 2 * belt
  extent_m <- new_n * res_m
  n_ru     <- extent_m / res_ru
  belt_ru  <- (belt * res_m) / res_ru
  
  g <- matrix(VAL_RU_FOREST, nrow = n_ru, ncol = n_ru)
  g[1:belt_ru, ]              <- VAL_RU_EDGE
  g[(n_ru-belt_ru+1):n_ru, ] <- VAL_RU_EDGE
  g[, 1:belt_ru]              <- VAL_RU_EDGE
  g[, (n_ru-belt_ru+1):n_ru] <- VAL_RU_EDGE
  
  r <- rast(nrows = n_ru, ncols = n_ru,
            xmin = 0, xmax = extent_m,
            ymin = 0, ymax = extent_m)
  values(r) <- as.vector(t(g))
  crs(r) <- "EPSG:25833"
  #crs(r)    <- "EPSG:5514"
  r
}

# ══════════════════════════════════════════════════════════
#  GENERATE
# ══════════════════════════════════════════════════════════

cat("=== iLand landscape generator — C4 2022 ===\n")
cat(sprintf("Total disturbed: %d cells = %.2f ha (%.1f%%)\n\n",
            target_size, target_size*res_m^2/10000, 100*target_frac))

r1   <- make_landscape(n_blobs = 1, seed = 42)
r2   <- make_landscape(n_blobs = 2, seed = 42)
r4   <- make_landscape(n_blobs = 4, seed = 42)
r8   <- make_landscape(n_blobs = 8, seed = 42)
r_ru <- make_ru_grid()

# ══════════════════════════════════════════════════════════
#  WRITE
# ══════════════════════════════════════════════════════════

writeRaster(r1,   "outData/stand_C4_2022_1blob.asc",
            filetype = "AAIGrid", overwrite = TRUE)
writeRaster(r2,   "outData/stand_C4_2022_2blobs.asc",
            filetype = "AAIGrid", overwrite = TRUE)
writeRaster(r4,   "outData/stand_C4_2022_4blobs.asc",
            filetype = "AAIGrid", overwrite = TRUE)
writeRaster(r8,   "outData/stand_C4_2022_8blobs.asc",
            filetype = "AAIGrid", overwrite = TRUE)
writeRaster(r_ru, "outData/climate_ru.asc",
            filetype = "AAIGrid", overwrite = TRUE)


# ══════════════════════════════════════════════════════════
#  VISUAL CHECK
# ══════════════════════════════════════════════════════════

cols   <- c("gray15", "pink", "darkgreen")
breaks <- c(-1.5, 0.5, 1.5, 2.5)

par(mfrow = c(1, 4), mar = c(2, 2, 3, 1))

plot(r1, col = cols, breaks = breaks, legend = FALSE,
     main = "1 blob — compact (3.86 ha)")
legend("topright", legend = c("edge","disturbed","mature"),
       fill = cols, cex = 0.7, bty = "n")

plot(r2, col = cols, breaks = breaks, legend = FALSE,
     main = "2 blobs (1.93 ha each)")

plot(r4, col = cols, breaks = breaks, legend = FALSE,
     main = "4 blobs (0.97 ha each)")

plot(r8, col = cols, breaks = breaks, legend = FALSE,
     main = "8 blobs (0.48 ha each)")


# Should show: 140 rows, 140 cols, extent 0–1400 m
print(r1)

# Belt check: edge cells should be -1, inner starts at row/col 21
values(r1)[1:5]        # should all be -1 (edge)
r1[21, 21]             # should be 1 or 2 (inner grid)
