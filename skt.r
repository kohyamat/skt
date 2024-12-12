# Required packages: glmmTMB
library(glmmTMB)

# Load source for functions
source("functions.r")

# Load data
df0 <- read.csv("./data/pfr_sub_100by100.csv.gz")

# data infomation
# -------------------
# subplot: subplot ID, total area was divided into by 100m x 100m grid
# species: species ID, integer

# Variables
dbh_min <- 1.0 # Boundary size in dbh (cm)
n_tree_min <- 2 # minimum popul. size
grid_size <- 100 # one-side length of sub-grid (m)
q_min <- 10 # minimum frequency in terms of subplot number

# Output directory
outdir <- "output"
dir.create(outdir, showWarnings = FALSE)

# Target plot area definition (all 50 ha)
area <- 50.0 #  area of entire plot in ha
area_sub <- 1.0 # subplot area in ha
n_sub <- area / area_sub # number of subplots


# Calculate turnover rates for each species in each subplot
turnover_subplot <- function(ds, dbh_min, area, n_tree_min) {
  n_surv <- tapply((ds$dbh1 >= dbh_min & ds$dbh2 >= dbh_min), ds$species, sum, na.rm = TRUE)
  rare_sp <- names(n_surv)[n_surv < n_tree_min]
  ds$species[as.character(ds$species) %in% rare_sp] <- "Others"
  dss <- split(ds, ds$species)
  res <- lapply(
    dss,
    function(x) data.frame(est_turnover_rates(x$dbh1, x$dbh2, x$t, dbh_min = dbh_min))
  )
  sp_names <- names(res)
  names(res) <- NULL
  res <- do.call(rbind, res)
  res$subplot <- unique(ds$subplot)
  res$species <- sp_names
  res$N <- res$N / area # per-ha period mean number of stems
  res$B <- res$B / area # per-ha period mean biomass
  res$Bl <- res$Bl / area # per-ha period mean biomass
  res$P <- (res$p_g + res$p_f) * res$B # absolute productivity rate
  res$L <- res$l * res$B # absolute loss rate
  to_left <- c("subplot", "species")
  return(res[, c(to_left, colnames(res)[!(colnames(res) %in% to_left)])])
}

# Data preparation
df0_split <- split(df0, df0$subplot) # split data by subplot
df1 <- lapply(df0_split, turnover_subplot, dbh_min, area_sub, n_tree_min)
names(df1) <- NULL
df1 <- do.call(rbind, df1)

write.csv(df1, file = file.path(outdir, "res_turnovers_100by100.csv"), row.names = FALSE, quote = FALSE)

# F ~ B parameter estimation (Eq. 1)
fitF <- glmmTMB(log(Bl) ~ log(B) + (1 + log(B) | species),
  data = df1,
)
bi <- coef(fitF)$cond$species[, 1]
ai <- coef(fitF)$cond$species[, 2]
predF <- data.frame(
  "bi" = exp(bi + sigma(fitF)^2 / 2), # fix the lognormal-distribution bias
  "ai" = ai,
  "species" = rownames(coef(fitF)$cond$species)
)

# Merge predictions to df1
df1 <- merge(df1, predF, by = "species")
df1$x <- df1$B^df1$ai
df1$bx <- df1$bi * df1$B^df1$ai

# Subplot-level aggragation
agg_func_sub <- function(dat) {
  with(dat, data.frame(
    "S" = sum(B > 0), # species richness per subplot
    "Bsum" = sum(B),
    "Fsum" = sum(Bl),
    "xsum" = sum(x),
    "bxsum" = sum(bi * x),
    "Psum" = sum(P)
  ))
}

df_sub <- do.call(rbind, lapply(split(df1, df1$subplot), agg_func_sub))
df_sub$subplot <- levels(as.factor(df1$subplot))

# Species-level aggregation
agg_func_sp <- function(dat, n_sub) {
  with(dat, data.frame(
    "Q" = sum(B > 0), # frequency in subplot number
    "N_obs" = sum(N) / n_sub,
    "F_obs" = sum(Bl) / n_sub,
    "B_obs" = sum(B) / n_sub, # landscape-wide biomass per ha
    "Wmax_obs" = mean(Umax),
    "bxsp" = sum(bx) / n_sub,
    "p_obs" = sum(p * B) / sum(B),
    "l_obs" = sum(l * B) / sum(B),
    "r_obs" = (sum(p * B) - sum(l * B)) / sum(B),
    "P_obs" = sum(p * B) / n_sub
  ))
}

df_sp <- do.call(rbind, lapply(split(df1, df1$species), agg_func_sp, n_sub))
df_sp$species <- levels(as.factor(df1$species))

# Merge aggregated data
df1 <- merge(df1, df_sub, by = "subplot")
df1 <- merge(df1, df_sp, by = "species")

# Calculate heterospecific variables
df1 <- within(df1, {
  xh <- xsum - x # heterospecific x
  bxh <- bxsum - bx # heterospecific bi*x
})

# Species selection
# q_min: minimum species frequency in number of subplots
df1 <- subset(df1, species != "Others" & Q >= q_min)

# -------------------
# Standardization of variables
df1 <- within(df1, {
  y <- x / mean(x)
  by <- bx / mean(bx)
  byh <- bxh / mean(bxh)
  bysum <- bxsum / mean(bxsum)
})

# producrivity_growth model
fitg <- glmmTMB(p_g ~ 1 + by + byh + (1 + by + byh | species),
  REML = TRUE,
  data = df1,
  start = list(theta = rep(0, 6))
) 

# prouctivity_recruit model
fitf <- glmmTMB(p_f ~ 1 + (1 | species),
  data = df1
)

# loss rate model
fitl <- glmmTMB(l ~ 1 + (1 | species),
  data = df1
)

dn <- data.frame(
  "species" = rownames(coef(fitg)$cond$species),
  "cf1" = coef(fitg)$cond$species[, 1],
  "cf2" = coef(fitg)$cond$species[, 2] / mean(df1$bx),
  "cf3" = coef(fitg)$cond$species[, 3] / mean(df1$bxh),
  "cff1" = coef(fitf)$cond$species[, 1],
  "cfl1" = coef(fitl)$cond$species[, 1]
)

dn <- merge(dn, predF, by = "species") # adding bi
dn <- merge(dn, df_sp, by = "species") # adding bi

# Parameters of Eq. 1
dn <- within(dn, {
  ri <- cf1 + cff1 - cfl1
  mi <- cfl1
  ui <- -cf2
  vi <- -cf3
  inv <- ri / vi # species invasibility
  qi <- vi / (ui - vi)
})

# Filter species
dn <- subset(dn, ri > 0) # remove lethal species (3 species)
dn <- subset(dn, ui > 0) # remove conspecific facilitation (one species)
dn <- subset(dn, vi > 0)

N <- nrow(dn) # Species pool, 485 spp.

#  Order species data by invasibility
dn <- dn[order(dn$inv, decreasing = TRUE), ]
dn <- within(dn, Ci <- cumsum(inv * qi) / (1 + cumsum(qi))) # C(i)

#  Equilibrium solution
dn <- within(dn, {
  coex <- sign(inv - Ci)
  bx_con <- ri / ui
  B_con <- (ri / (ui * bi))^(1 / ai)
})

dP <- subset(dn, coex > 0) # positive species
dQ <- subset(dn, coex <= 0) # zero species
s <- nrow(dP) # number of coexisting species
Cs <- dn$Ci[s] # C(s)

dP <- within(dP, {
  bx_eq <- (inv - Cs) * qi
  B_eq <- (bx_eq / bi)^(1 / ai)
  cum_bx_eq <- cumsum(bx_eq)
  cum_B_eq <- cumsum(B_eq)
  P_eq <- cfl1 * B_eq
})

B_eq <- sum(dP$B_eq)
P_eq <- sum(dP$P_eq)


# Output
# -------------------
# Draw Fig. 1 of the main text
pdf(file.path(outdir, "fig1.pdf"))
plot(dn$inv,
  type = "l", ylim = c(0, 25),
  col = "royalblue", lwd = 2,
  xlab = "Species rank",
  ylab = "Leaf biomass (Mg/ha)"
)
lines(dn$Ci, col = "orangered", lwd = 2)
lines(dP$cum_bx_eq, lwd = 2)
lines(c(s, N), c(Cs, Cs), lwd = 2)
points(dP$bx_eq * 100, col = "dimgray", type = "o")
abline(v = s, lty = 2)
lines(c(0, nrow(dn)), c(0, 0), lwd = 0.6)
dev.off()

# Community-level summary
sprintf("N, s: %i, %i", N, s)
sprintf("B_obs, B_equil: %.3f, %.3f", sum(df_sub$Bsum) / n_sub, B_eq)
sprintf("P_obs, P_equil: %.3f, %.3f", sum(df_sub$Psum) / n_sub, P_eq)


dn <- within(dn, {
  F_equil <- 0
  B_equil <- 0
  P_equil <- 0
})
for (i in seq_len(nrow(dP))) {
  dn$F_equil[i] <- dP$bx_eq[i]
  dn$B_equil[i] <- dP$B_eq[i]
  dn$P_equil[i] <- dP$P_eq[i]
}

# List of species parameters to be output
output_cols <- c(
  "species",
  "F_obs",
  "B_obs",
  "Wmax_obs",
  "p_obs",
  "l_obs",
  "r_obs",
  "P_obs",
  "ri",
  "mi",
  "ui",
  "vi",
  "ai",
  "bi",
  "inv",
  "Ci",
  "F_equil",
  "B_equil",
  "P_equil"
)

write.csv(dn[output_cols], file.path(outdir, "skt_list.csv"), quote = FALSE, row.names = FALSE)
