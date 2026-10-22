## =========================================================================
## 03_evaluacion.R  —  Métricas, pruebas de superioridad y evaluación
##                     con origen móvil (expanding window)
## =========================================================================

## -------------------------------------------------------------------------
## Métricas puntuales
## -------------------------------------------------------------------------
rmse <- function(e) sqrt(mean(e^2, na.rm = TRUE))
mae  <- function(e) mean(abs(e), na.rm = TRUE)

## MASE: escala = MAE del naive de un paso en la muestra de entrenamiento
escala_mase <- function(y_train) mean(abs(diff(y_train)), na.rm = TRUE)

mape <- function(y, yhat) {
  ok <- is.finite(y) & is.finite(yhat) & abs(y) > 1e-8
  if (!any(ok)) return(NA_real_)
  mean(abs((y[ok] - yhat[ok]) / y[ok])) * 100
}

smape <- function(y, yhat) {
  den <- (abs(y) + abs(yhat)) / 2
  ok <- is.finite(den) & den > 1e-8
  if (!any(ok)) return(NA_real_)
  mean(abs(y[ok] - yhat[ok]) / den[ok]) * 100
}

## -------------------------------------------------------------------------
## Métricas de densidad / intervalos
## -------------------------------------------------------------------------
## CRPS a partir de una muestra predictiva (fórmula de energía)
crps_muestra <- function(y, x) {
  x <- x[is.finite(x)]
  m <- length(x)
  if (m < 2 || !is.finite(y)) return(NA_real_)
  xs <- sort(x)
  t1 <- mean(abs(xs - y))
  ## E|X - X'| con muestra ordenada: 2/m^2 * sum_i (2i - m - 1) x_(i)
  i <- seq_len(m)
  t2 <- 2 * sum((2 * i - m - 1) * xs) / (m^2)
  t1 - 0.5 * t2
}

## Pérdida pinball en un cuantil
pinball <- function(y, q, tau) {
  if (!is.finite(y) || !is.finite(q)) return(NA_real_)
  if (y >= q) tau * (y - q) else (1 - tau) * (q - y)
}

## Cobertura y ancho del intervalo
cobertura <- function(y, lo, hi) as.numeric(is.finite(lo) & is.finite(hi) & y >= lo & y <= hi)
ancho     <- function(lo, hi) hi - lo

## -------------------------------------------------------------------------
## Pruebas de superioridad predictiva
## -------------------------------------------------------------------------
## Diebold-Mariano con corrección de muestra pequeña (Harvey-Leybourne-Newbold)
dm_hln <- function(e1, e2, h = 1, potencia = 2) {
  ok <- is.finite(e1) & is.finite(e2)
  e1 <- e1[ok]; e2 <- e2[ok]
  n <- length(e1)
  if (n < 8) return(c(estadistico = NA_real_, p = NA_real_))
  d <- abs(e1)^potencia - abs(e2)^potencia
  dbar <- mean(d)
  ## varianza de largo plazo con truncamiento en h-1
  gam <- function(k) sum((d[(k + 1):n] - dbar) * (d[1:(n - k)] - dbar)) / n
  v <- gam(0)
  if (h > 1) for (k in seq_len(h - 1)) v <- v + 2 * gam(k)
  if (!is.finite(v) || v <= 0) return(c(estadistico = NA_real_, p = NA_real_))
  DM <- dbar / sqrt(v / n)
  corr <- sqrt((n + 1 - 2 * h + h * (h - 1) / n) / n)
  DM <- DM * corr
  p <- 2 * pt(-abs(DM), df = n - 1)
  c(estadistico = DM, p = p)
}

## -------------------------------------------------------------------------
## Evaluación con origen móvil
## -------------------------------------------------------------------------
## Devuelve un data.frame con una fila por (modelo, origen, horizonte)

evaluar_serie <- function(y, z = NULL, s = NULL,
                          modelos = MODELOS,
                          horizontes = HORIZONTES,
                          prop_train = 0.70, prop_val = 0.15,
                          n_origenes = 10L,
                          niveles = NIVELES,
                          verbose = FALSE) {
  n <- length(y)
  h_max <- max(horizontes)
  ini_test <- floor(n * (prop_train + prop_val))
  fin_test <- n - h_max
  if (fin_test <= ini_test) stop("serie demasiado corta para el diseño de evaluación")
  origenes <- unique(round(seq(ini_test, fin_test, length.out = n_origenes)))

  filas <- list(); k <- 0
  for (o in origenes) {
    y_tr <- y[1:o]
    z_tr <- if (!is.null(z)) z[1:o] else NULL
    s_tr <- if (!is.null(s)) s[1:o] else NULL
    esc  <- escala_mase(y_tr)
    y_fu <- y[(o + 1):(o + h_max)]

    for (mid in names(modelos)) {
      mod <- modelos[[mid]]
      ## modelos que requieren z o s se omiten si el DGP no los provee
      if (isTRUE(mod$usa_z) && is.null(z)) next
      if (isTRUE(mod$usa_s) && is.null(s)) next
      fc <- correr_modelo(mod, y_tr, h_max, z = z_tr, s = s_tr)
      for (hh in horizontes) {
        yh <- y_fu[hh]
        pm <- if (length(fc$mean) >= hh) fc$mean[hh] else NA_real_
        e  <- yh - pm
        simh <- if (!is.null(fc$sim)) fc$sim[, hh] else NULL
        k <- k + 1
        filas[[k]] <- data.frame(
          modelo   = mid,
          etiqueta = mod$etiqueta,
          familia  = mod$familia,
          origen   = o,
          h        = hh,
          y_real   = yh,
          y_hat    = pm,
          error    = e,
          escala   = esc,
          crps     = if (!is.null(simh)) crps_muestra(yh, simh) else NA_real_,
          cob80    = cobertura(yh, fc$lower[hh, 1], fc$upper[hh, 1]),
          amp80    = ancho(fc$lower[hh, 1], fc$upper[hh, 1]),
          cob95    = cobertura(yh, fc$lower[hh, 2], fc$upper[hh, 2]),
          amp95    = ancho(fc$lower[hh, 2], fc$upper[hh, 2]),
          pin05    = if (!is.null(simh)) pinball(yh, quantile(simh, 0.05, na.rm = TRUE), 0.05) else NA_real_,
          pin95    = if (!is.null(simh)) pinball(yh, quantile(simh, 0.95, na.rm = TRUE), 0.95) else NA_real_,
          tiempo   = fc$tiempo,
          ok       = isTRUE(fc$ok),
          stringsAsFactors = FALSE
        )
      }
      if (verbose) cat(sprintf("  origen %d | %-28s | ok=%s | %.2fs\n",
                               o, mod$etiqueta, isTRUE(fc$ok), fc$tiempo))
    }
  }
  do.call(rbind, filas)
}

## -------------------------------------------------------------------------
## Agregación por réplica: una fila por (modelo, h)
## -------------------------------------------------------------------------
resumir_replica <- function(df) {
  clave <- interaction(df$modelo, df$h, drop = TRUE)
  out <- lapply(split(df, clave), function(g) {
    esc <- mean(g$escala, na.rm = TRUE)
    data.frame(
      modelo   = g$modelo[1],
      etiqueta = g$etiqueta[1],
      familia  = g$familia[1],
      h        = g$h[1],
      rmse     = rmse(g$error),
      mae      = mae(g$error),
      mase     = if (is.finite(esc) && esc > 0) mae(g$error) / esc else NA_real_,
      smape    = smape(g$y_real, g$y_hat),
      crps     = mean(g$crps, na.rm = TRUE),
      picp80   = mean(g$cob80, na.rm = TRUE),
      mpiw80   = mean(g$amp80, na.rm = TRUE),
      picp95   = mean(g$cob95, na.rm = TRUE),
      mpiw95   = mean(g$amp95, na.rm = TRUE),
      pinball  = mean(c(g$pin05, g$pin95), na.rm = TRUE),
      tiempo   = mean(g$tiempo, na.rm = TRUE),
      tasa_ok  = mean(g$ok),
      n_orig   = nrow(g),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

## -------------------------------------------------------------------------
## Diebold-Mariano por réplica contra un modelo pivote
## -------------------------------------------------------------------------
dm_vs_pivote <- function(df, pivote = "M1") {
  res <- list(); k <- 0
  for (hh in unique(df$h)) {
    d <- df[df$h == hh, ]
    piv <- d[d$modelo == pivote, ]
    if (!nrow(piv)) next
    piv <- piv[order(piv$origen), ]
    for (mid in setdiff(unique(d$modelo), pivote)) {
      m <- d[d$modelo == mid, ]
      m <- m[order(m$origen), ]
      comun <- intersect(piv$origen, m$origen)
      if (length(comun) < 8) next
      e1 <- piv$error[match(comun, piv$origen)]
      e2 <- m$error[match(comun, m$origen)]
      tst <- dm_hln(e1, e2, h = hh)
      k <- k + 1
      res[[k]] <- data.frame(modelo = mid, pivote = pivote, h = hh,
                             dm = unname(tst[1]), p = unname(tst[2]),
                             mejor_que_pivote = unname(tst[1]) > 0,
                             stringsAsFactors = FALSE)
    }
  }
  if (!k) return(NULL)
  do.call(rbind, res)
}
