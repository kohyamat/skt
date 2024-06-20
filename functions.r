#' Calculate biomass of a tree
#'
#' This function calculates the biomass of a tree in Pasoh Forest Reserve based on
#' the diameter at breast height (DBH) using allometric equations developed by
#' Niiyama et al. (2010).
#'
#' @param dbh Diameter at breast height in centimeters.
#' @param component Character string indicating the component of biomass to calculate.
#'                  Possible values are: 'agb', 'all', 'stem', 'leaf', 'root'.
#'                  Default is 'agb'.
#'
#' @details The calculation of biomass is based on the equations described in the
#' following publication:
#'   Niiyama, K., Kajimoto, T., Matsuura, Y., et al. (2010). Estimation of root biomass
#'   based on excavation of individual root systems in a primary dipterocarp forest in
#'   Pasoh Forest Reserve, Peninsular Malaysia. Journal of Tropical Ecology, 26(3),
#'   271-284. doi:10.1017/S0266467410000040
#'
#' @return The calculated biomass based on the specified component (kg in dry mass).
#'
#' @examples
#' biomass_pfr(30, "agb")
#' biomass_pfr(30, "all")
#' biomass_pfr(30, "stem")
#' biomass_pfr(30, "leaf")
#' biomass_pfr(30, "root")
#'
#' @export
biomass_pfr <- function(dbh, component = "agb") {
  h <- 1 / (1 / (1.61 * dbh) + 1 / 69.) # tree height (m)
  w_s <- 0.036 * (dbh * dbh * h)^1.01 # trunk plus branch
  w_l <- 1 / (1 / (0.108 * w_s^0.75) + 1 / 105) # leaves
  w_r <- 0.023 * dbh^2.59 # coarse root
  if (component == "agb") {
    return(w_s + w_l)
  } else if (component == "all") {
    return(w_s + w_l + w_r)
  } else if (component == "stem") {
    return(w_s)
  } else if (component == "leaf") {
    return(w_l)
  } else if (component == "root") {
    return(w_r)
  } else {
    stop(
      sprintf(
        "Invalid component '%s'. Expected one of: 'all', 'agb', 'stem', 'leaf', 'root'",
        component
      )
    )
  }
}


#' Calculate the period mean of two values
#'
#' This function calculates the period mean between two values \code{x1} and \code{x2}.
#' The period mean is calculated as (x2 - x1) / log(x2 / x1).
#'
#' @param x1 Numeric value for the first value.
#' @param x2 Numeric value for the second value.
#' @param omit_na Logical indicating whether to omit NA values when calculating the mean.
#'                Default is \code{FALSE}.
#'
#' @return The calculated period mean.
#'
#' @examples
#' period_mean(2, 10)
#' period_mean(5, 5)
#' period_mean(3, NA)
#' period_mean(NA, 7, omit_na = TRUE)
#'
#' @export
period_mean <- function(x1, x2, omit_na = FALSE) {
  if (is.na(x1) || is.na(x2)) {
    if (omit_na) {
      return(mean(c(x1, x2), na.rm = TRUE))
    } else {
      return(NA)
    }
  } else if (x1 == x2) {
    return(x1)
  } else {
    return((x2 - x1) / log(x2 / x1))
  }
}


#' Calculate turnover rate
#'
#' This function calculates the turnover rate based on the observed values of the
#' survival status or biomass of individual organisms at time 1 and time 2 (z, y)
#' and the time interval (t).
#'
#' @param y A vector representing the observed survival status or biomass at time 2.
#' @param z A vector representing the observed survival status or biomass at time 1.
#' @param t A scalar or vector representing the time interval(s) between the two observations.
#' @param init The initial value for the turnover rate. Default is 0.01.
#' @param precision The desired precision for convergence. Default is 1.0e-12.
#'
#' @return The calculated turnover rate.
#'
#' @examples
#' y <- c(100, 200, 300)
#' z <- c(50, 150, 250)
#' t <- 1
#' turnover(y, z, t)
#'
#' @export
turnover <- function(y, z, t, init = 0.01, precision = 1.0e-12) {
  if (identical(y, z)) {
    return(0.0)
  } else {
    f <- function(rho) {
      sum(y * exp(-rho * t) - z)
    }
    fprime <- function(rho) {
      sum(-t * y * exp(-rho * t))
    }
    # Newton-Raphson iteration
    rho <- init
    change <- precision + 1.0
    while (change > precision) {
      rho2 <- rho - f(rho) / fprime(rho)
      change <- abs(rho2 - rho)
      rho <- rho2
    }
    return(rho)
  }
}


#' Estimate turnover rates based on diameter at breast height (DBH)
#'
#' This function estimates various turnover rates based on the changes in DBH over time.
#'
#' @param dbh1 Numeric vector of initial DBH values.
#' @param dbh2 Numeric vector of final DBH values.
#' @param t Numeric vector representing the time interval between measurements.
#' @param dbh_min Numeric value representing the minimum DBH threshold.
#'
#' @return A list of estimated turnover rates and other metrics.
#'
#' @examples
#' dbh1 <- c(10, 12, 15, 9)
#' dbh2 <- c(13, 14, 17, 10)
#' t <- c(2, 2, 2, 2)
#' est_turnover_rates(dbh1, dbh2, t, dbh_min = 1.0)
#'
#' @importFrom stats quantile
#'
#' @export
est_turnover_rates <- function(dbh1, dbh2, t, dbh_min = 1.0) {
  dbh1[is.na(dbh1)] <- 0
  dbh2[is.na(dbh2)] <- 0

  w1 <- biomass_pfr(dbh1, "agb") / 1000 # in Mg
  w2 <- biomass_pfr(dbh2, "agb") / 1000 # in Mg
  wl1 <- biomass_pfr(dbh1, "leaf") / 1000 # in Mg
  wl2 <- biomass_pfr(dbh2, "leaf") / 1000 # in Mg

  q1 <- as.numeric((dbh1 >= dbh_min) & (dbh2 >= dbh_min))
  q2 <- as.numeric((dbh1 >= dbh_min) & (dbh2 < dbh_min))
  q3 <- as.numeric((dbh1 < dbh_min) & (dbh2 >= dbh_min))

  n_s0 <- sum(q1, na.rm = TRUE)
  n_0 <- n_s0 + sum(q2, na.rm = TRUE)
  n_t <- n_s0 + sum(q3, na.rm = TRUE)

  b_s0 <- sum(q1 * w1, na.rm = TRUE)
  b_st <- sum(q1 * w2, na.rm = TRUE)
  b_0 <- b_s0 + sum(q2 * w1, na.rm = TRUE)
  b_t <- b_st + sum(q3 * w2, na.rm = TRUE)

  bl_s0 <- sum(q1 * wl1, na.rm = TRUE)
  bl_st <- sum(q1 * wl2, na.rm = TRUE)
  bl_0 <- bl_s0 + sum(q2 * wl1, na.rm = TRUE)
  bl_t <- bl_st + sum(q3 * wl2, na.rm = TRUE)

  n_m <- period_mean(n_0, n_t) # period mean number of stems
  b_m <- period_mean(b_0, b_t) # period mean biomass
  bl_m <- period_mean(bl_0, bl_t) # period mean leaf biomass

  # r <- turnover(q1 + q3, q1, t) # relative recruitment rate
  # m <- turnover(q1 + q2, q1, t) # relative mortality rate
  p <- turnover(w2, q1 * w1, t) # relative above ground biomass productivity rate
  l <- turnover(w1, q1 * w1, t) # relative above ground biomass loss rate
  pl <- turnover(wl2, q1 * wl1, t) # relative leaf mass productivity rate
  ll <- turnover(wl1, q1 * wl1, t) # relative leaf mass loss rate

  u_0 <- max(w1)
  u_t <- max(w2)
  umax <- period_mean(u_0, u_t)
  ux <- as.numeric(quantile(w1[!q3], 0.99))
  udmax <- max(dbh1[!q3])
  umaw0 <- max(w1[!q3])

  return(list(
    "t" = mean(t * (w1 + w2)) / mean(w1 + w2),
    "N" = n_m,
    "B" = b_m,
    "Bl" = bl_m,
    # "r" = r,
    # "m" = m,
    "p" = p,
    "l" = l,
    "pl" = pl,
    "ll" = ll,
    "Umax" = umax,
    "Ux" = ux,
    "Udmax" = udmax,
    "Umaw0" = umaw0
  ))
}
