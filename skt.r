library(glmmTMB)

source("turnover_subplot.r")

# Load data
df0 <- read.csv("data/pfr.csv.gz")

# Variables
dbh_min <- 1. # Boundary size in dbh (cm)
n_tree_min <- 2 # minimum popul. size
grid_size <- 100 # one-side length of sub-grid (m)
q_min <- 10 # minimum frequency in terms of subplot number

# Output directory
outdir <- "output"
dir.create(outdir, showWarnings = FALSE)

# Flitering
df0$species <- df0$Cd
df0 <- subset(df0, (dbh1 >= dbh_min) | (dbh2 >= dbh_min))
df0 <- subset(df0, dbh2 - dbh1 < 100)

# Target plot area definition (all 50 ha)
x_min <- 0
x_max <- 1000
y_min <- 0
y_max <- 500
df0 <- subset(df0, (x >= x_min) & (x < x_max) & (y >= y_min) & (y < y_max))
df0$x <- df0$x - x_min
df0$y <- df0$y - y_min
x_len <- x_max - x_min # total plot x-axis length in m
y_len <- y_max - y_min # total plot y-axis length in m
area <- (x_len * y_len) / 10000 #  area of entire plot in ha

# Subplot (grid) definition
area_sub <- grid_size * grid_size / 10000 # subplot area in ha
n_sub <- area / area_sub # number of subplots

# Calculate turnover rates for each species in each subplot
df1 <- turnover_subplot(df0, grid_size, area_sub, dbh_min, n_tree_min)

write.csv(df1, file = file.path(outdir, "res_turnovers_100.csv"), row.names = FALSE, quote = FALSE)


# F ~ B parameter estimation (Eq. 1)
fitF <- glmmTMB(log(Bl) ~ log(B) + (1 | species),
  weights = B^(2 / 3),
  data = df1, family = "gaussian"
)

# Weights for subplot-level match between observed and predicted leaf mass
bi <- coef(fitF)$cond$species[, 1]
k <- mean(coef(fitF)$cond$species[, 2])

predF <- data.frame(
  "bi" = exp(bi + sigma(fitF)^2 / 2), # fixing the log-distribution bias
  "species" = rownames(coef(fitF)$cond$species)
)

df1 <- merge(df1, predF, by = "species")


df1$x <- df1$B^k
df1$bx <- df1$bi * df1$B^k


# Subplot-level aggragation
agg_func_sub <- function(dat) {
  with(dat, data.frame(
    "S" = sum(ifelse(B > 0, 1, 0)), # species richness per subplot
    "Bsum" = sum(B),
    "Fsum" = sum(Bl),
    "xsum" = sum(x),
    "bxsum" = sum(bi * x),
    "Psum" = sum(P)
  ))
}

df1_by_sub <- split(df1, df1$sub_id)
df_sub <- lapply(df1_by_sub, agg_func_sub)
names(df_sub) <- NULL
df_sub <- do.call(rbind, df_sub)
df_sub$sub_id <- names(df1_by_sub)


# Species-level aggregation
agg_func_sp <- function(dat) {
  with(dat, data.frame(
    "Q" = sum(ifelse(B > 0, 1, 0)), # frequency in subplot number
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

df1_by_sp <- split(df1, df1$species)
df_sp <- lapply(df1_by_sp, agg_func_sp)
names(df_sp) <- NULL
df_sp <- do.call(rbind, df_sp)
df_sp$species <- names(df1_by_sp)

df1 <- merge(df1, df_sub, by = "sub_id")
df1 <- merge(df1, df_sp, by = "species")

df1$xh <- df1$xsum - df1$x # heterospecific x
df1$bxh <- df1$bxsum - df1$bx # heterospecific bi*x


# Species selection
df1 <- subset(df1, df1$species != "Others")
df1 <- subset(df1, df1$Q >= q_min) # minimum species frequency in number of subplots


# -------------------
# Standardization of variables
df1$y <- df1$x / mean(df1$x)
df1$by <- df1$bx / mean(df1$bx)
df1$byh <- df1$bxh / mean(df1$bxh)
df1$bysum <- df1$bxsum / mean(df1$bxsum)

# Woody producrivity model
fitp <- glmmTMB(p ~ 1 + by + byh + (1 + by + byh | species),
  REML = TRUE,
  data = df1, family = "gaussian"
) # selected model by AIC

# Woody loss rate model
fitl <- glmmTMB(l ~ 1 + (1 | species),
  data = df1, family = "gaussian"
)

dn <- data.frame(
  "species" = rownames(coef(fitp)$cond$species),
  "cf1" = coef(fitp)$cond$species[, 1],
  "cf2" = coef(fitp)$cond$species[, 2] / mean(df1$bx),
  "cf3" = coef(fitp)$cond$species[, 3] / mean(df1$bxh),
  "cfl1" = coef(fitl)$cond$species[, 1]
)

dn <- merge(dn, predF, by = "species") # adding bi
dn <- merge(dn, df_sp, by = "species") # adding bi


# Parameters of Eq. 1
dn$ri <- dn$cf1 - dn$cfl1
dn$di <- dn$cfl1
dn$ui <- -dn$cf2
dn$vi <- -dn$cf3

dn$inv <- dn$ri / dn$vi # species invasibility
dn$qi <- dn$vi / (dn$ui - dn$vi)

dn <- subset(dn, ri > 0) # remove lethal species (3 species)
dn <- subset(dn, ui > 0) # remove conspecific facilitation (one species)
dn <- subset(dn, vi > 0)

N <- nrow(dn) # Species pool, 485 spp.

#  Order species data by invasibility
dn <- dn[order(dn$inv, decreasing = TRUE), ]
dn <- within(dn, Ci <- cumsum(inv * qi) / (1 + cumsum(qi))) # C(i)

#  Equilibrium solution
dn$coex <- sign(dn$inv - dn$Ci)
dn$bx_con <- dn$ri / dn$ui
dn$B_con <- (dn$ri / (dn$ui * dn$bi))^(1 / k)

dP <- subset(dn, coex > 0) # positive species
dQ <- subset(dn, coex <= 0) # zero species
s <- nrow(dP) # number of coexisting species
Cs <- dn$Ci[s] # C(s)

dP$bx_eq <- (dP$inv - Cs) * dP$qi
dP$B_eq <- (dP$bx_eq / dP$bi)^(1 / k)
dP$cum_bx_eq <- cumsum(dP$bx_eq)
dP$cum_B_eq <- cumsum(dP$B_eq)
B_eq <- sum(dP$B_eq)
dP$P_eq <- dP$cfl1 * dP$B_eq
P_eq <- sum(dP$P_eq)

# Figure 1
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
print(c("N, s: "))
print(c(N, s))
print(c("B_obs, B_equil: "))
print(c(sum(df_sub$Bsum) / n_sub, B_eq))
print(c("P_obs, P_equil: "))
print(c(sum(df_sub$Psum) / n_sub, P_eq))

dn$F_equil <- 0
dn$B_equil <- 0
dn$P_equil <- 0
for (i in seq_len(nrow(dP))) {
  dn$F_equil[i] <- dP$bx_eq[i]
  dn$B_equil[i] <- dP$B_eq[i]
  dn$P_equil[i] <- dP$P_eq[i]
}


# Replacing species code by species name
sp_list <- read.csv("data/Pasoh_species.csv")
sp_list$species <- paste(sp_list$Genus, sp_list$Species, sep = " ")
dn$Cd <- dn$species
dn$species <- NULL

dn <- merge(dn, sp_list, by = "Cd")

# List of species parameters
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
  "di",
  "ui",
  "vi",
  "bi",
  "inv",
  "Ci",
  "F_equil",
  "B_equil",
  "P_equil"
)

write.csv(dn[output_cols], file.path(outdir, "skt_list.csv"), quote = FALSE, row.names = FALSE)
