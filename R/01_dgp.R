## =========================================================================
## 01_dgp.R  —  Procesos generadores de datos (DGP)
##
## Cada generador devuelve una lista con:
##   y : serie simulada (numeric, longitud Tn)
##   z : variable de umbral continua exógena (numeric o NULL)
##   s : régimen categórico observable (integer o NULL)
##   info : lista con los parámetros efectivos del DGP
## =========================================================================

## -------------------------------------------------------------------------
## Bloque lineal
## -------------------------------------------------------------------------

gen_ar1 <- function(Tn, phi = 0.5, sigma = 1, c0 = 0, burn = 200) {
  n <- Tn + burn
  e <- rnorm(n, 0, sigma)
  y <- numeric(n)
  for (t in 2:n) y[t] <- c0 + phi * y[t - 1] + e[t]
  list(y = y[(burn + 1):n], z = NULL, s = NULL,
       info = list(dgp = "AR(1)", phi = phi, sigma = sigma))
}

gen_arma11 <- function(Tn, phi = 0.6, theta = 0.4, sigma = 1, c0 = 0, burn = 200) {
  n <- Tn + burn
  e <- rnorm(n, 0, sigma)
  y <- numeric(n)
  for (t in 2:n) y[t] <- c0 + phi * y[t - 1] + e[t] + theta * e[t - 1]
  list(y = y[(burn + 1):n], z = NULL, s = NULL,
       info = list(dgp = "ARMA(1,1)", phi = phi, theta = theta, sigma = sigma))
}

gen_ma1 <- function(Tn, theta = 0.6, sigma = 1, c0 = 0, burn = 200) {
  n <- Tn + burn
  e <- rnorm(n, 0, sigma)
  y <- numeric(n)
  for (t in 2:n) y[t] <- c0 + e[t] + theta * e[t - 1]
  list(y = y[(burn + 1):n], z = NULL, s = NULL,
       info = list(dgp = "ARIMA(0,0,1)", theta = theta, sigma = sigma))
}

gen_arima110 <- function(Tn, phi = 0.4, mu = 0.02, sigma = 1, burn = 200) {
  n <- Tn + burn
  e <- rnorm(n, 0, sigma)
  dy <- numeric(n)
  for (t in 2:n) dy[t] <- mu + phi * dy[t - 1] + e[t]
  y <- cumsum(dy)
  list(y = y[(burn + 1):n], z = NULL, s = NULL,
       info = list(dgp = "ARIMA(1,1,0)", phi = phi, mu = mu, sigma = sigma))
}

gen_arima111 <- function(Tn, phi = 0.4, theta = 0.3, mu = 0.02, sigma = 1, burn = 200) {
  n <- Tn + burn
  e <- rnorm(n, 0, sigma)
  dy <- numeric(n)
  for (t in 2:n) dy[t] <- mu + phi * dy[t - 1] + e[t] + theta * e[t - 1]
  y <- cumsum(dy)
  list(y = y[(burn + 1):n], z = NULL, s = NULL,
       info = list(dgp = "ARIMA(1,1,1)", phi = phi, theta = theta,
                   mu = mu, sigma = sigma))
}

## -------------------------------------------------------------------------
## Bloque de umbrales
## -------------------------------------------------------------------------

## D6 — SETAR(2;1,1): variable de umbral CONTINUA y AUTOEXCITADA (y_{t-d})
gen_setar <- function(Tn, phi1 = 0.85, phi2 = 0.30,
                      c1 = -0.30, c2 = 0.30,
                      r = 0, d = 1, sigma = 1, burn = 200) {
  ## Condición suficiente de estacionariedad ergódica (Petruccelli-Woolford)
  stopifnot(phi1 < 1, phi2 < 1, phi1 * phi2 < 1)
  n <- Tn + burn
  e <- rnorm(n, 0, sigma)
  y <- numeric(n)
  for (t in (d + 1):n) {
    if (y[t - d] <= r) y[t] <- c1 + phi1 * y[t - 1] + e[t]
    else               y[t] <- c2 + phi2 * y[t - 1] + e[t]
  }
  reg <- as.integer(y[seq_len(n)] > r) + 1L
  list(y = y[(burn + 1):n], z = NULL, s = NULL,
       info = list(dgp = "SETAR(2;1,1)", phi1 = phi1, phi2 = phi2,
                   c1 = c1, c2 = c2, r = r, d = d, sigma = sigma,
                   delta = abs(phi1 - phi2),
                   prop_reg2 = mean(reg[(burn + 1):n] == 2L)))
}

## D7 — TAR con variable de umbral CONTINUA EXÓGENA z_{t-1}
gen_tar_zcont <- function(Tn, phi1 = 0.85, phi2 = 0.30,
                          c1 = -0.30, c2 = 0.30,
                          rho_z = 0.7, sigma = 1, burn = 200) {
  n <- Tn + burn
  e <- rnorm(n, 0, sigma)
  u <- rnorm(n, 0, 1)
  z <- numeric(n)
  for (t in 2:n) z[t] <- rho_z * z[t - 1] + u[t]
  r <- stats::median(z)                     # umbral = mediana de z
  y <- numeric(n)
  for (t in 2:n) {
    if (z[t - 1] <= r) y[t] <- c1 + phi1 * y[t - 1] + e[t]
    else               y[t] <- c2 + phi2 * y[t - 1] + e[t]
  }
  list(y = y[(burn + 1):n], z = z[(burn + 1):n], s = NULL,
       info = list(dgp = "TAR-z continua", phi1 = phi1, phi2 = phi2,
                   c1 = c1, c2 = c2, r = r, rho_z = rho_z, sigma = sigma,
                   delta = abs(phi1 - phi2)))
}

## D8 — TAR con variable de umbral CATEGÓRICA EXÓGENA s_t (regímenes observables)
P_BASE <- matrix(c(0.95, 0.04, 0.01,
                   0.03, 0.94, 0.03,
                   0.02, 0.08, 0.90),
                 nrow = 3, byrow = TRUE)

gen_tar_scat <- function(Tn,
                         P      = P_BASE,
                         cs     = c(0.40,  0.00, -0.60),
                         phis   = c(0.30,  0.70,  0.85),
                         sigmas = c(0.80,  1.00,  2.20),
                         escala_sigma = 1,
                         misclass = 0, burn = 200) {
  K <- nrow(P)
  n <- Tn + burn
  s <- integer(n); s[1] <- 2L
  for (t in 2:n) s[t] <- sample.int(K, 1L, prob = P[s[t - 1], ])
  y <- numeric(n)
  for (t in 2:n) {
    k <- s[t]
    y[t] <- cs[k] + phis[k] * y[t - 1] + rnorm(1, 0, sigmas[k] * escala_sigma)
  }
  s_obs <- s
  if (misclass > 0) {                       # indicador de régimen con error
    flip <- runif(n) < misclass
    if (any(flip)) s_obs[flip] <- sample.int(K, sum(flip), replace = TRUE)
  }
  list(y = y[(burn + 1):n], z = NULL, s = s_obs[(burn + 1):n],
       info = list(dgp = "TAR-s categórica", K = K, cs = cs, phis = phis,
                   sigmas = sigmas * escala_sigma, misclass = misclass,
                   frec_reg = as.numeric(table(factor(s[(burn + 1):n],
                                                      levels = 1:K))) / Tn))
}

## -------------------------------------------------------------------------
## Rejilla de DGP del experimento
## -------------------------------------------------------------------------
## Cada entrada: id, etiqueta, familia, función generadora y argumentos fijos.
## `sigma` se sobrescribe con el factor de ruido del diseño factorial.

DGP_GRID <- list(
  D1a = list(id = "D1a", etiqueta = "AR(1) phi=0.5",       familia = "lineal",
             fn = gen_ar1,       args = list(phi = 0.5)),
  D1b = list(id = "D1b", etiqueta = "AR(1) phi=0.9",       familia = "lineal",
             fn = gen_ar1,       args = list(phi = 0.9)),
  D2  = list(id = "D2",  etiqueta = "ARMA(1,1)",           familia = "lineal",
             fn = gen_arma11,    args = list(phi = 0.6, theta = 0.4)),
  D3  = list(id = "D3",  etiqueta = "ARIMA(1,1,1)",        familia = "integrado",
             fn = gen_arima111,  args = list(phi = 0.4, theta = 0.3, mu = 0.02)),
  D4  = list(id = "D4",  etiqueta = "ARIMA(1,1,0)",        familia = "integrado",
             fn = gen_arima110,  args = list(phi = 0.4, mu = 0.02)),
  D5  = list(id = "D5",  etiqueta = "ARIMA(0,0,1)",        familia = "lineal",
             fn = gen_ma1,       args = list(theta = 0.6)),
  D6a = list(id = "D6a", etiqueta = "SETAR delta baja",    familia = "umbral_auto",
             fn = gen_setar,     args = list(phi1 = 0.65, phi2 = 0.45,
                                             c1 = -0.15, c2 = 0.15)),
  D6b = list(id = "D6b", etiqueta = "SETAR delta alta",    familia = "umbral_auto",
             fn = gen_setar,     args = list(phi1 = 0.85, phi2 = 0.30,
                                             c1 = -0.30, c2 = 0.30)),
  D7  = list(id = "D7",  etiqueta = "TAR umbral continuo", familia = "umbral_cont",
             fn = gen_tar_zcont, args = list(phi1 = 0.85, phi2 = 0.30,
                                             c1 = -0.30, c2 = 0.30)),
  D8  = list(id = "D8",  etiqueta = "TAR umbral categórico", familia = "umbral_cat",
             fn = gen_tar_scat,  args = list())
)

## Genera una serie de la celda (dgp_id, Tn, sigma)
simular_dgp <- function(dgp_id, Tn, sigma = 1, burn = 200) {
  esp <- DGP_GRID[[dgp_id]]
  if (is.null(esp)) stop("DGP no reconocido: ", dgp_id)
  args <- esp$args
  args$Tn <- Tn
  args$burn <- burn
  ## El nivel de ruido entra como sigma; en el DGP categórico escala las
  ## varianzas por régimen preservando su heterocedasticidad relativa.
  if (identical(esp$id, "D8")) args$escala_sigma <- sigma else args$sigma <- sigma
  out <- do.call(esp$fn, args)
  out$info$dgp_id   <- esp$id
  out$info$etiqueta <- esp$etiqueta
  out$info$familia  <- esp$familia
  out$info$Tn       <- Tn
  out$info$sigma    <- sigma
  out
}
