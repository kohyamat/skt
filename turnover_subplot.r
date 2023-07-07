source("functions.r")

turnover_subplot <- function(d, grid_size, area, dbh_min, n_tree_min) {
  sx <- floor(d$x / grid_size)
  sy <- floor(d$y / grid_size)
  d$sub_id <- sprintf("%02i_%02i", sx, sy)

  # Split data frame by subplot
  d_by_sub <- split(d, d$sub_id)

  # Calculate turnover rates for each species in each subplot
  apply_func <- function(ds, dbh_min, area, n_tree_min) {
    n_surv <- tapply((ds$dbh1 >= dbh_min & ds$dbh2 >= dbh_min), ds$species, sum, na.rm = TRUE)
    rare_sp <- names(n_surv)[n_surv < n_tree_min]
    ds$species[ds$species %in% rare_sp] <- "Others"
    ds_by_sp <- split(ds, ds$species)
    res <- lapply(ds_by_sp, function(x) data.frame(est_turnover_rates(x$dbh1, x$dbh2, x$t, dbh_min = dbh_min)))
    names(res) <- NULL
    res <- do.call(rbind, res)
    res$sub_id <- unique(ds$sub_id)
    res$species <- names(ds_by_sp)
    res$N <- res$N / area # per-ha period mean number of stems
    res$B <- res$B / area # per-ha period mean biomass
    res$Bl <- res$Bl / area # per-ha period mean biomass
    res$R <- res$r * res$N # absolute recruitment rate
    res$M <- res$m * res$N # absolute mortality rate
    res$P <- res$p * res$B # absolute productivity rate
    res$L <- res$l * res$B # absolute loss rate
    toleft <- c("sub_id", "species")
    return(res[, c(toleft, colnames(res)[!(colnames(res) %in% toleft)])])
  }

  out <- lapply(d_by_sub, apply_func, dbh_min, area, n_tree_min)
  names(out) <- NULL
  out <- do.call(rbind, out)

  return(out)
}
