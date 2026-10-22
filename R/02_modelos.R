## =========================================================================
## 02_modelos.R  —  Modelos de pronóstico
##
## Interfaz común:
##   m_<modelo>(y, h, z = NULL, s = NULL, ...)
##     y : serie de entrenamiento (numeric)
##     h : horizonte máximo
##     z : variable de umbral continua exógena, alineada con y (o NULL)
##     s : régimen categórico observable, alineado con y (o NULL)
##   devuelve list(mean, lower, upper, sim, ok, msg)
##
## TAR/SETAR, MS-AR y el modelo bayesiano están implementados aquí y no
## dependen de paquetes externos: eso garantiza reproducibilidad, control
## de la densidad predictiva (necesaria para CRPS y cobertura) y soporte
## uniforme para umbral autoexcitado, continuo exógeno y categórico.
## =========================================================================

## -------------------------------------------------------------------------
## Utilidades
## -------------------------------------------------------------------------

salida_prono <- function(mean, sim = NULL, niveles = NIVELES,
                         ok = TRUE, msg = "") {
  h <- length(mean)
  lower <- upper <- matrix(NA_real_, h, length(niveles),
                           dimnames = list(NULL, paste0("q", niveles)))
  if (!is.null(sim) && all(is.finite(sim))) {
    for (j in seq_along(niveles)) {
      a <- (1 - niveles[j] / 100) / 2
      lower[, j] <- apply(sim, 2, quantile, probs = a,     na.rm = TRUE)
      upper[, j] <- apply(sim, 2, quantile, probs = 1 - a, na.rm = TRUE)
    }
  }
  list(mean = mean, lower = lower, upper = upper, sim = sim, ok = ok, msg = msg)
}

prono_fallido <- function(h, msg = "fallo") {
  salida_prono(rep(NA_real_, h), NULL, ok = FALSE, msg = msg)
}

## Matriz de rezagos: devuelve list(Y, X) alineados
matriz_rezagos <- function(y, p) {
  n <- length(y)
  Y <- y[(p + 1):n]
  X <- sapply(seq_len(p), function(i) y[(p + 1 - i):(n - i)])
  X <- matrix(X, ncol = p)
  colnames(X) <- paste0("lag", seq_len(p))
  list(Y = Y, X = X)
}

## Diferenciación e integración de pronósticos --------------------------------
orden_diferenciacion <- function(y, max_d = 2) {
  d <- tryCatch(forecast::ndiffs(y, test = "kpss", max.d = max_d),
                error = function(e) 0L)
  as.integer(min(d, max_d))
}

diferenciar <- function(y, d) if (d == 0) y else diff(y, differences = d)

## Reintegra media y trayectorias simuladas de la escala diferenciada a niveles
integrar_prono <- function(mean_d, sim_d, y_orig, d) {
  if (d == 0) return(list(mean = mean_d, sim = sim_d))
  n <- length(y_orig)
  mean_i <- mean_d
  sim_i  <- sim_d
  for (k in seq_len(d)) {
    ## nivel de anclaje: última observación de la serie diferenciada (d-k) veces
    base_serie <- if (d - k == 0) y_orig else diff(y_orig, differences = d - k)
    ancla <- base_serie[length(base_serie)]
    mean_i <- ancla + cumsum(mean_i)
    if (!is.null(sim_i)) sim_i <- ancla + t(apply(sim_i, 1, cumsum))
  }
  list(mean = mean_i, sim = sim_i)
}

## Envoltura genérica: diferencia si hace falta, ajusta y reintegra
con_diferenciacion <- function(fun, y, h, z = NULL, s = NULL, max_d = 2, ...) {
  d <- orden_diferenciacion(y, max_d = max_d)
  yd <- diferenciar(y, d)
  zd <- if (!is.null(z) && d > 0) z[-seq_len(d)] else z
  sd_ <- if (!is.null(s) && d > 0) s[-seq_len(d)] else s
  res <- fun(yd, h, z = zd, s = sd_, ...)
  if (!isTRUE(res$ok)) return(res)
  ii <- integrar_prono(res$mean, res$sim, y, d)
  out <- salida_prono(ii$mean, ii$sim)
  ## conserva los diagnósticos del modelo (K, umbral, tipo, orden, motor...)
  extra <- setdiff(names(res), c("mean", "lower", "upper", "sim", "ok", "msg"))
  for (nm in extra) out[[nm]] <- res[[nm]]
  out$d_aplicado <- d
  out
}

## -------------------------------------------------------------------------
## M0 — Naive / random walk
## -------------------------------------------------------------------------
m_naive <- function(y, h, z = NULL, s = NULL, B = B_PRED, ...) {
  n <- length(y)
  dy <- diff(y)
  if (length(dy) < 5) return(prono_fallido(h, "serie muy corta"))
  mean_f <- rep(y[n], h)
  ## trayectorias por bootstrap de las primeras diferencias (media cero)
  dy_c <- dy - mean(dy)
  paths <- matrix(sample(dy_c, B * h, replace = TRUE), nrow = B, ncol = h)
  sim <- y[n] + t(apply(paths, 1, cumsum))
  salida_prono(mean_f, sim)
}

## -------------------------------------------------------------------------
## M1 — ARIMA automático (benchmark lineal)
## -------------------------------------------------------------------------
m_arima <- function(y, h, z = NULL, s = NULL, B = B_PRED, ic = "bic", ...) {
  fit <- tryCatch(forecast::auto.arima(y, ic = ic, stepwise = TRUE,
                                       approximation = FALSE, seasonal = FALSE),
                  error = function(e) NULL)
  if (is.null(fit)) return(prono_fallido(h, "auto.arima falló"))
  fc <- tryCatch(forecast::forecast(fit, h = h, level = 95),
                 error = function(e) NULL)
  if (is.null(fc)) return(prono_fallido(h, "forecast falló"))
  mean_f <- as.numeric(fc$mean)
  sdv <- (as.numeric(fc$upper[, 1]) - mean_f) / qnorm(0.975)
  sdv[!is.finite(sdv) | sdv <= 0] <- sd(residuals(fit), na.rm = TRUE)
  sim <- matrix(rnorm(B * h, rep(mean_f, each = B), rep(sdv, each = B)),
                nrow = B, ncol = h)
  out <- salida_prono(mean_f, sim)
  out$orden <- paste(forecast::arimaorder(fit), collapse = ",")
  out
}

## -------------------------------------------------------------------------
## M2 — TAR / SETAR  (umbral autoexcitado, continuo exógeno o categórico)
## -------------------------------------------------------------------------
## tipo: "auto"  -> variable de umbral q_t = y_{t-d}
##       "zcont" -> variable de umbral q_t = z_{t-d}  (continua exógena)
##       "scat"  -> regímenes dados por s_t           (categórica observable)

tar_ajustar <- function(y, z = NULL, s = NULL, p = 1, d = 1,
                        tipo = c("auto", "zcont", "scat"),
                        n_grid = 40, min_prop = 0.15) {
  tipo <- match.arg(tipo)
  n <- length(y)
  ## el arranque debe cubrir tanto los rezagos (p) como el retardo de la
  ## variable de umbral (d), de lo contrario q y Y quedan desalineados
  m0 <- max(p, if (tipo == "scat") 0L else d)
  if (n - m0 < 4 * (p + 2)) return(NULL)
  idx <- (m0 + 1):n                      # índices de y asociados a cada fila
  Y <- y[idx]
  X <- cbind(1, matrix(vapply(seq_len(p), function(i) y[idx - i],
                              numeric(length(idx))), ncol = p))

  ajustar_por_regimen <- function(reg, K) {
    coef <- vector("list", K); resid <- rep(NA_real_, length(Y))
    sig <- rep(NA_real_, K); ssr <- 0
    for (k in seq_len(K)) {
      sel <- which(reg == k)
      if (length(sel) < ncol(X) + 3) return(NULL)
      fitk <- stats::lm.fit(X[sel, , drop = FALSE], Y[sel])
      coef[[k]] <- fitk$coefficients
      resid[sel] <- fitk$residuals
      sig[k] <- sqrt(sum(fitk$residuals^2) / max(1, length(sel) - ncol(X)))
      ssr <- ssr + sum(fitk$residuals^2)
    }
    list(coef = coef, resid = resid, sigma = sig, ssr = ssr, reg = reg, K = K)
  }

  if (tipo == "scat") {
    if (is.null(s)) stop("tipo 'scat' requiere s")
    K <- length(unique(s))
    reg <- as.integer(factor(s[idx], levels = sort(unique(s))))
    fit <- ajustar_por_regimen(reg, K)
    if (is.null(fit)) return(NULL)
    fit$tipo <- "scat"; fit$p <- p; fit$d <- 0; fit$r <- NA
    fit$niveles_s <- sort(unique(s))
    return(fit)
  }

  q <- if (tipo == "auto") y[idx - d] else z[idx - d]
  if (any(!is.finite(q))) return(NULL)
  cand <- unique(quantile(q, probs = seq(min_prop, 1 - min_prop,
                                         length.out = n_grid), na.rm = TRUE))
  mejor <- NULL
  for (r in cand) {
    reg <- as.integer(q > r) + 1L
    f <- ajustar_por_regimen(reg, 2L)
    if (is.null(f)) next
    if (is.null(mejor) || f$ssr < mejor$ssr) { mejor <- f; mejor$r <- r }
  }
  if (is.null(mejor)) return(NULL)
  mejor$tipo <- tipo; mejor$p <- p; mejor$d <- d; mejor$q <- q
  mejor
}

## Selección de p (y d en el caso autoexcitado) por BIC
tar_seleccionar <- function(y, z = NULL, s = NULL,
                            tipo = c("auto", "zcont", "scat"),
                            p_max = 3, d_max = 2) {
  tipo <- match.arg(tipo)
  n <- length(y)
  mejor <- NULL; mejor_bic <- Inf
  ds <- if (tipo == "scat") 0 else seq_len(d_max)
  for (p in seq_len(p_max)) for (d in ds) {
    f <- tryCatch(tar_ajustar(y, z, s, p = p, d = d, tipo = tipo),
                  error = function(e) NULL)
    if (is.null(f)) next
    nobs <- length(f$resid)
    k <- f$K * (p + 1) + (tipo != "scat")
    bic <- nobs * log(f$ssr / nobs) + k * log(nobs)
    if (is.finite(bic) && bic < mejor_bic) { mejor <- f; mejor_bic <- bic }
  }
  if (!is.null(mejor)) mejor$bic <- mejor_bic
  mejor
}

## Pronóstico por simulación de trayectorias
tar_pronosticar <- function(fit, y, h, z = NULL, s = NULL, B = B_PRED) {
  p <- fit$p; d <- fit$d; K <- fit$K
  n <- length(y)
  res_por_reg <- split(fit$resid, fit$reg)
  res_por_reg <- lapply(res_por_reg, function(v) v[is.finite(v)])

  ## Trayectoria futura de la variable de umbral
  z_fut <- NULL; s_fut <- NULL
  if (fit$tipo == "zcont") {
    fz <- tryCatch(stats::ar(z, order.max = 1, aic = FALSE, method = "mle"),
                   error = function(e) NULL)
    rho <- if (!is.null(fz)) as.numeric(fz$ar[1]) else 0
    mz  <- if (!is.null(fz)) fz$x.mean else mean(z)
    sz  <- if (!is.null(fz)) sqrt(fz$var.pred) else sd(z)
    z_fut <- matrix(NA_real_, B, h + d)
    z_fut[, seq_len(d)] <- matrix(rep(z[(length(z) - d + 1):length(z)], each = B),
                                  nrow = B)
    for (j in (d + 1):(h + d))
      z_fut[, j] <- mz + rho * (z_fut[, j - 1] - mz) + rnorm(B, 0, sz)
  }
  if (fit$tipo == "scat") {
    niv <- fit$niveles_s
    sf <- factor(s, levels = niv)
    P <- matrix(0, K, K)
    tab <- table(head(sf, -1), tail(sf, -1))
    P <- sweep(as.matrix(tab) + 1e-6, 1, rowSums(as.matrix(tab) + 1e-6), "/")
    s_fut <- matrix(NA_integer_, B, h)
    actual <- rep(match(s[length(s)], niv), B)
    cumP <- t(apply(P, 1, cumsum))
    for (j in seq_len(h)) {
      actual <- max.col(runif(B) < cumP[actual, , drop = FALSE],
                        ties.method = "first")
      s_fut[, j] <- actual
    }
  }

  hist <- matrix(rep(y[(n - max(p, d) + 1):n], each = B), nrow = B)
  paths <- matrix(NA_real_, B, h)
  buf <- cbind(hist, matrix(NA_real_, B, h))
  off <- ncol(hist)
  for (j in seq_len(h)) {
    reg <- switch(fit$tipo,
      auto  = as.integer(buf[, off + j - d] > fit$r) + 1L,
      zcont = as.integer(z_fut[, j] > fit$r) + 1L,
      scat  = s_fut[, j])
    xj <- cbind(1, buf[, off + j - seq_len(p), drop = FALSE])
    mu <- numeric(B)
    ej <- numeric(B)
    for (k in seq_len(K)) {
      sel <- which(reg == k)
      if (!length(sel)) next
      mu[sel] <- as.numeric(xj[sel, , drop = FALSE] %*% fit$coef[[k]])
      rk <- res_por_reg[[as.character(k)]]
      ej[sel] <- if (length(rk) > 5) sample(rk, length(sel), replace = TRUE)
                 else rnorm(length(sel), 0, fit$sigma[k])
    }
    buf[, off + j] <- mu + ej
    paths[, j] <- buf[, off + j]
  }
  list(mean = colMeans(paths, na.rm = TRUE), sim = paths)
}

m_tar <- function(y, h, z = NULL, s = NULL, tipo = NULL, B = B_PRED, ...) {
  if (is.null(tipo)) tipo <- if (!is.null(s)) "scat" else if (!is.null(z)) "zcont" else "auto"
  fit <- tryCatch(tar_seleccionar(y, z, s, tipo = tipo),
                  error = function(e) NULL)
  if (is.null(fit)) return(prono_fallido(h, "TAR no estimable"))
  pr <- tryCatch(tar_pronosticar(fit, y, h, z, s, B = B),
                 error = function(e) NULL)
  if (is.null(pr)) return(prono_fallido(h, "TAR pronóstico falló"))
  out <- salida_prono(pr$mean, pr$sim)
  out$umbral <- fit$r; out$p <- fit$p; out$d <- fit$d; out$tipo <- fit$tipo
  out
}

## Variantes explícitas usadas en el experimento
m_setar    <- function(y, h, z = NULL, s = NULL, ...) m_tar(y, h, NULL, NULL, tipo = "auto",  ...)
m_tar_z    <- function(y, h, z = NULL, s = NULL, ...) if (is.null(z)) prono_fallido(h, "sin z") else m_tar(y, h, z, NULL, tipo = "zcont", ...)
m_tar_s    <- function(y, h, z = NULL, s = NULL, ...) if (is.null(s)) prono_fallido(h, "sin s") else m_tar(y, h, NULL, s, tipo = "scat",  ...)

## -------------------------------------------------------------------------
## M3 — MS-AR: cambio de régimen markoviano (filtro de Hamilton + ML)
## -------------------------------------------------------------------------
msar_P <- function(a, K) {
  P <- matrix(0, K, K); idx <- 1
  for (i in seq_len(K)) {
    fila <- numeric(K)
    for (j in seq_len(K)) {
      if (j == i) fila[j] <- 0 else { fila[j] <- a[idx]; idx <- idx + 1 }
    }
    ex <- exp(pmin(pmax(fila, -20), 20))
    P[i, ] <- ex / sum(ex)
  }
  P
}

msar_desempacar <- function(theta, K, orden) {
  cs   <- theta[1:K]
  phis <- matrix(theta[(K + 1):(K + K * orden)], nrow = K)
  lsig <- pmin(pmax(theta[(K + K * orden + 1):(2 * K + K * orden)], -6), 3)
  a    <- theta[(2 * K + K * orden + 1):length(theta)]
  list(cs = cs, phis = phis, sig = exp(lsig), P = msar_P(a, K))
}

## NOTA: los argumentos NO se llaman `p` para evitar el emparejamiento parcial
## con el formal `par` de optim(), que rompe silenciosamente la optimización.
msar_loglik <- function(theta, Y, X, K, orden, devolver_filtro = FALSE) {
  ## ruta rápida: todo el cálculo ocurre en C++; es la que ejecuta el
  ## optimizador miles de veces, así que no debe pasar por código R.
  if (isTRUE(HAMILTON_CPP) && !devolver_filtro)
    return(msar_nll_cpp(theta, Y, X, as.integer(K), as.integer(orden)))

  th <- msar_desempacar(theta, K, orden)
  n <- length(Y)
  ## distribución ergódica como estado inicial
  xi0 <- tryCatch({
    A <- rbind(diag(K) - t(th$P), rep(1, K))
    as.numeric(qr.solve(A, c(rep(0, K), 1)))
  }, error = function(e) rep(1 / K, K))
  if (any(!is.finite(xi0)) || any(xi0 < 0)) xi0 <- rep(1 / K, K)

  if (isTRUE(HAMILTON_CPP)) {
    r <- hamilton_filtro(Y, X, th$cs, th$phis, th$sig, th$P, xi0, devolver_filtro)
    if (!isTRUE(r$ok)) return(if (devolver_filtro) NULL else 1e10)
    if (devolver_filtro)
      return(list(ll = r$ll, xi_T = r$xi_T, filt = r$filt, th = th))
    return(-r$ll)
  }

  eta <- matrix(0, n, K)
  for (k in seq_len(K)) {
    mu <- th$cs[k] + as.numeric(X %*% th$phis[k, ])
    eta[, k] <- dnorm(Y, mu, th$sig[k])
  }
  eta[!is.finite(eta)] <- 1e-300
  eta <- pmax(eta, 1e-300)
  xi <- xi0
  ll <- 0; filt <- matrix(0, n, K)
  for (t in seq_len(n)) {
    pred <- as.numeric(t(th$P) %*% xi)
    num <- pred * eta[t, ]
    f <- sum(num)
    if (!is.finite(f) || f <= 0) return(if (devolver_filtro) NULL else 1e10)
    ll <- ll + log(f)
    xi <- num / f
    filt[t, ] <- xi
  }
  if (devolver_filtro) return(list(ll = ll, xi_T = xi, filt = filt, th = th))
  -ll
}

msar_ajustar <- function(y, K = 2, p = 1) {
  lr <- matriz_rezagos(y, p); Y <- lr$Y; X <- lr$X
  ## valores iniciales: partición por cuantiles del nivel
  cortes <- quantile(Y, probs = seq(0, 1, length.out = K + 1))
  grp <- cut(Y, breaks = unique(cortes), include.lowest = TRUE, labels = FALSE)
  cs <- numeric(K); phis <- matrix(0, K, p); lsig <- numeric(K)
  for (k in seq_len(K)) {
    sel <- which(grp == k)
    if (length(sel) < p + 3) sel <- seq_along(Y)
    fk <- stats::lm.fit(cbind(1, X[sel, , drop = FALSE]), Y[sel])
    cs[k] <- fk$coefficients[1]
    phis[k, ] <- fk$coefficients[-1]
    lsig[k] <- log(max(sd(fk$residuals), 1e-3))
  }
  a0 <- rep(log(0.05 / 0.9), K * (K - 1))
  par0 <- c(cs, as.numeric(phis), lsig, a0)
  op <- tryCatch(optim(par = par0, fn = msar_loglik,
                       Y = Y, X = X, K = K, orden = p,
                       method = "Nelder-Mead",
                       control = list(maxit = 4000, reltol = 1e-8)),
                 error = function(e) NULL)
  if (is.null(op) || !is.finite(op$value) || op$value >= 1e9) return(NULL)
  ft <- msar_loglik(op$par, Y, X, K, p, devolver_filtro = TRUE)
  if (is.null(ft)) return(NULL)
  list(par = op$par, ll = ft$ll, xi_T = ft$xi_T, th = ft$th, K = K, p = p,
       npar = length(op$par), nobs = length(Y))
}

msar_pronosticar <- function(fit, y, h, B = B_PRED) {
  K <- fit$K; p <- fit$p; th <- fit$th
  n <- length(y)
  estado <- sample.int(K, B, replace = TRUE, prob = fit$xi_T)
  buf <- cbind(matrix(rep(y[(n - p + 1):n], each = B), nrow = B),
               matrix(NA_real_, B, h))
  off <- p
  paths <- matrix(NA_real_, B, h)
  ## muestreo vectorizado del régimen por CDF inversa (evita B llamadas a
  ## sample.int por paso, que dominaban el costo del pronóstico)
  cumP <- t(apply(th$P, 1, cumsum))
  for (j in seq_len(h)) {
    estado <- max.col(runif(B) < cumP[estado, , drop = FALSE],
                      ties.method = "first")
    xj <- buf[, off + j - seq_len(p), drop = FALSE]
    mu <- th$cs[estado] + rowSums(xj * th$phis[estado, , drop = FALSE])
    buf[, off + j] <- mu + rnorm(B, 0, th$sig[estado])
    paths[, j] <- buf[, off + j]
  }
  list(mean = colMeans(paths), sim = paths)
}

m_msar <- function(y, h, z = NULL, s = NULL, K_grid = c(2, 3), p = 1,
                   B = B_PRED, ...) {
  mejor <- NULL; mejor_bic <- Inf
  for (K in K_grid) {
    f <- tryCatch(msar_ajustar(y, K = K, p = p), error = function(e) NULL)
    if (is.null(f)) next
    bic <- -2 * f$ll + f$npar * log(f$nobs)
    if (is.finite(bic) && bic < mejor_bic) { mejor <- f; mejor_bic <- bic }
  }
  if (is.null(mejor)) return(prono_fallido(h, "MS-AR no estimable"))
  pr <- tryCatch(msar_pronosticar(mejor, y, h, B = B), error = function(e) NULL)
  if (is.null(pr)) return(prono_fallido(h, "MS-AR pronóstico falló"))
  out <- salida_prono(pr$mean, pr$sim)
  out$K <- mejor$K; out$bic <- mejor_bic
  out
}

## -------------------------------------------------------------------------
## M4 — Bayesiano: AR(p) con prior Normal-Gamma Inversa (tipo Minnesota)
## -------------------------------------------------------------------------
## Posterior conjugada -> densidad predictiva por simulación (sin MCMC).
## Es el análogo tratable de un BSTS con componente autorregresiva; el
## componente de nivel local se añade en la aplicación empírica.

m_bayes_ar <- function(y, h, z = NULL, s = NULL, p = 4, B = B_PRED,
                       lambda = 0.2, ...) {
  n <- length(y)
  if (n < p + 10) return(prono_fallido(h, "serie corta para AR(p) bayesiano"))
  lr <- matriz_rezagos(y, p)
  Y <- lr$Y; X <- cbind(1, lr$X)
  k <- ncol(X); nobs <- length(Y)

  ## Prior: media de random walk en el primer rezago, shrinkage tipo Minnesota
  b0 <- c(0, 1, rep(0, p - 1))
  vpr <- c(1e6, (lambda / seq_len(p))^2)
  V0i <- diag(1 / vpr, k)
  a0 <- 3; d0 <- 2 * var(Y)

  XtX <- crossprod(X); XtY <- crossprod(X, Y)
  Vn  <- solve(V0i + XtX)
  bn  <- as.numeric(Vn %*% (V0i %*% b0 + XtY))
  an  <- a0 + nobs / 2
  dn  <- as.numeric(d0 + 0.5 * (sum(Y^2) + t(b0) %*% V0i %*% b0 -
                                  t(bn) %*% (V0i + XtX) %*% bn))
  if (!is.finite(dn) || dn <= 0) dn <- max(1e-6, 0.5 * sum((Y - X %*% bn)^2))

  sig2 <- 1 / rgamma(B, shape = an, rate = dn)
  Ch <- tryCatch(chol(Vn), error = function(e) NULL)
  if (is.null(Ch)) return(prono_fallido(h, "posterior no definida positiva"))
  beta <- matrix(bn, B, k, byrow = TRUE) +
    (matrix(rnorm(B * k), B, k) %*% Ch) * sqrt(sig2)

  buf <- cbind(matrix(rep(y[(n - p + 1):n], each = B), nrow = B),
               matrix(NA_real_, B, h))
  paths <- matrix(NA_real_, B, h)
  for (j in seq_len(h)) {
    xj <- cbind(1, buf[, p + j - seq_len(p), drop = FALSE])
    mu <- rowSums(xj * beta)
    buf[, p + j] <- mu + rnorm(B, 0, sqrt(sig2))
    paths[, j] <- buf[, p + j]
  }
  out <- salida_prono(colMeans(paths), paths)
  out$p <- p
  out
}

## -------------------------------------------------------------------------
## M5 — Machine learning de árboles sobre rezagos
## -------------------------------------------------------------------------
## Usa xgboost si está disponible; si no, gbm; si no, randomForest.
## Pronóstico recursivo; intervalos por bootstrap de residuos in-sample.

m_ml_arboles <- function(y, h, z = NULL, s = NULL, p = 12, B = B_PRED,
                         motor = NULL, ...) {
  n <- length(y)
  if (n < p + 30) return(prono_fallido(h, "serie corta para ML"))
  lr <- matriz_rezagos(y, p); Y <- lr$Y; X <- lr$X
  if (is.null(motor)) {
    motor <- if (PKG_OK[["xgboost"]]) "xgboost" else
             if (PKG_OK[["gbm"]]) "gbm" else "rf"
  }
  pred_fun <- NULL
  if (motor == "xgboost") {
    fit <- tryCatch(xgboost::xgboost(data = X, label = Y, nrounds = 200,
                                     max_depth = 3, eta = 0.05,
                                     subsample = 0.8, colsample_bytree = 0.8,
                                     objective = "reg:squarederror",
                                     verbose = 0),
                    error = function(e) NULL)
    if (is.null(fit)) return(prono_fallido(h, "xgboost falló"))
    pred_fun <- function(newx) as.numeric(predict(fit, newx))
  } else if (motor == "gbm") {
    df <- data.frame(Y = Y, X)
    fit <- tryCatch(gbm::gbm(Y ~ ., data = df, distribution = "gaussian",
                             n.trees = 500, interaction.depth = 3,
                             shrinkage = 0.05, bag.fraction = 0.8,
                             n.minobsinnode = 5, verbose = FALSE),
                    error = function(e) NULL)
    if (is.null(fit)) return(prono_fallido(h, "gbm falló"))
    pred_fun <- function(newx) {
      nd <- as.data.frame(newx); names(nd) <- colnames(X)
      as.numeric(gbm::predict.gbm(fit, nd, n.trees = 500))
    }
  } else {
    fit <- tryCatch(randomForest::randomForest(x = X, y = Y, ntree = 500),
                    error = function(e) NULL)
    if (is.null(fit)) return(prono_fallido(h, "randomForest falló"))
    pred_fun <- function(newx) as.numeric(predict(fit, newx))
  }
  ## Escala de los residuos: los residuos in-sample de un ensamble de árboles
  ## están fuertemente sesgados a la baja y producirían intervalos absurdamente
  ## angostos. Se usa un holdout temporal (último 20 % del entrenamiento) para
  ## estimar la dispersión predictiva de forma honesta.
  resid <- Y - pred_fun(X)
  n_ho <- max(10L, floor(0.2 * length(Y)))
  if (length(Y) - n_ho > p + 20) {
    tr <- seq_len(length(Y) - n_ho); ho <- setdiff(seq_along(Y), tr)
    fit_ho <- tryCatch({
      if (motor == "xgboost")
        xgboost::xgboost(data = X[tr, , drop = FALSE], label = Y[tr],
                         nrounds = 200, max_depth = 3, eta = 0.05,
                         subsample = 0.8, colsample_bytree = 0.8,
                         objective = "reg:squarederror", verbose = 0)
      else if (motor == "gbm")
        gbm::gbm(Y ~ ., data = data.frame(Y = Y[tr], X[tr, , drop = FALSE]),
                 distribution = "gaussian", n.trees = 500,
                 interaction.depth = 3, shrinkage = 0.05,
                 bag.fraction = 0.8, n.minobsinnode = 5, verbose = FALSE)
      else randomForest::randomForest(x = X[tr, , drop = FALSE], y = Y[tr],
                                      ntree = 500)
    }, error = function(e) NULL)
    if (!is.null(fit_ho)) {
      ph <- tryCatch({
        if (motor == "gbm") {
          nd <- as.data.frame(X[ho, , drop = FALSE]); names(nd) <- colnames(X)
          as.numeric(gbm::predict.gbm(fit_ho, nd, n.trees = 500))
        } else as.numeric(predict(fit_ho, X[ho, , drop = FALSE]))
      }, error = function(e) NULL)
      if (!is.null(ph)) {
        e_ho <- Y[ho] - ph
        if (sd(resid) > 0 && is.finite(sd(e_ho)))
          resid <- resid * (sd(e_ho) / sd(resid))   # reescala a dispersión honesta
      }
    }
  }

  ## media: recursión determinista
  hist <- y[(n - p + 1):n]
  buf <- c(hist, rep(NA_real_, h))
  for (j in seq_len(h)) {
    xj <- matrix(buf[p + j - seq_len(p)], nrow = 1)
    colnames(xj) <- colnames(X)
    buf[p + j] <- pred_fun(xj)
  }
  mean_f <- buf[(p + 1):(p + h)]

  ## trayectorias: recursión con residuos remuestreados
  Bp <- min(B, 200L)
  bufB <- cbind(matrix(rep(hist, each = Bp), nrow = Bp),
                matrix(NA_real_, Bp, h))
  for (j in seq_len(h)) {
    xj <- bufB[, p + j - seq_len(p), drop = FALSE]
    colnames(xj) <- colnames(X)
    bufB[, p + j] <- pred_fun(xj) + sample(resid, Bp, replace = TRUE)
  }
  sim <- bufB[, (p + 1):(p + h), drop = FALSE]
  out <- salida_prono(mean_f, sim)
  out$motor <- motor
  out
}

## -------------------------------------------------------------------------
## M6 — Red neuronal autorregresiva (NNAR) como representante de DL
## -------------------------------------------------------------------------
## En el Monte Carlo se usa nnetar (feed-forward autorregresiva, liviana).
## La LSTM/GRU en torch se reserva para la aplicación empírica; si torch
## está disponible, m_lstm la reemplaza aquí sin cambiar la interfaz.

m_nnar <- function(y, h, z = NULL, s = NULL, B = 200L, orden = NULL,
                   size = NULL, ...) {
  ## nnetar usa missing() para p y size: pasar NULL rompe la llamada
  args <- list(y = stats::as.ts(y), repeats = 20)
  if (!is.null(orden)) args$p <- orden
  if (!is.null(size))  args$size <- size
  fit <- tryCatch(do.call(forecast::nnetar, args), error = function(e) NULL)
  if (is.null(fit)) return(prono_fallido(h, "nnetar falló"))
  fc <- tryCatch(forecast::forecast(fit, h = h), error = function(e) NULL)
  if (is.null(fc)) return(prono_fallido(h, "nnetar pronóstico falló"))
  sim <- tryCatch({
    m <- matrix(NA_real_, B, h)
    for (b in seq_len(B)) m[b, ] <- as.numeric(simulate(fit, nsim = h, future = TRUE))
    m
  }, error = function(e) NULL)
  salida_prono(as.numeric(fc$mean), sim)
}

## -------------------------------------------------------------------------
## Registro de modelos del Monte Carlo
## -------------------------------------------------------------------------
## `usa_z` / `usa_s`: el modelo recibe la variable de umbral cuando existe.
## `diferenciar`: aplicar la envoltura de diferenciación automática.

MODELOS <- list(
  M0 = list(id = "M0", etiqueta = "Naive (RW)",        familia = "benchmark",
            fn = m_naive,       diferenciar = FALSE, usa_z = FALSE, usa_s = FALSE),
  M1 = list(id = "M1", etiqueta = "ARIMA (auto)",      familia = "lineal",
            fn = m_arima,       diferenciar = FALSE, usa_z = FALSE, usa_s = FALSE),
  M2 = list(id = "M2", etiqueta = "SETAR autoexcitado", familia = "no lineal clásico",
            fn = m_setar,       diferenciar = TRUE,  usa_z = FALSE, usa_s = FALSE),
  M2z = list(id = "M2z", etiqueta = "TAR umbral continuo (z)", familia = "no lineal clásico",
            fn = m_tar_z,       diferenciar = TRUE,  usa_z = TRUE,  usa_s = FALSE),
  M2s = list(id = "M2s", etiqueta = "TAR umbral categórico (s)", familia = "no lineal clásico",
            fn = m_tar_s,       diferenciar = TRUE,  usa_z = FALSE, usa_s = TRUE),
  M3 = list(id = "M3", etiqueta = "MS-AR",             familia = "no lineal clásico",
            fn = m_msar,        diferenciar = TRUE,  usa_z = FALSE, usa_s = FALSE),
  M4 = list(id = "M4", etiqueta = "AR(p) bayesiano",   familia = "bayesiano",
            fn = m_bayes_ar,    diferenciar = TRUE,  usa_z = FALSE, usa_s = FALSE),
  M5 = list(id = "M5", etiqueta = "ML de árboles",     familia = "machine learning",
            fn = m_ml_arboles,  diferenciar = TRUE,  usa_z = FALSE, usa_s = FALSE),
  M6 = list(id = "M6", etiqueta = "Red neuronal AR",   familia = "deep learning",
            fn = m_nnar,        diferenciar = TRUE,  usa_z = FALSE, usa_s = FALSE)
)

## Ejecuta un modelo respetando su configuración
correr_modelo <- function(mod, y, h, z = NULL, s = NULL) {
  zz <- if (isTRUE(mod$usa_z)) z else NULL
  ss <- if (isTRUE(mod$usa_s)) s else NULL
  t0 <- proc.time()[["elapsed"]]
  out <- tryCatch({
    if (isTRUE(mod$diferenciar))
      con_diferenciacion(mod$fn, y, h, z = zz, s = ss)
    else
      mod$fn(y, h, z = zz, s = ss)
  }, error = function(e) prono_fallido(h, paste("error:", conditionMessage(e))))
  out$tiempo <- proc.time()[["elapsed"]] - t0
  out$modelo <- mod$id
  out
}
