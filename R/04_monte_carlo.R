## =========================================================================
## 04_monte_carlo.R  —  Motor del experimento factorial
##
## Diseño: DGP (10 celdas) x T (2) x sigma (2) x réplicas
## Para cada serie simulada se evalúan todos los modelos del registro con
## origen móvil y se guardan los resultados agregados por réplica.
## =========================================================================

## Rejilla factorial del experimento
construir_diseno <- function(dgps = names(DGP_GRID),
                             Ts = c(250L, 1000L),
                             sigmas = c(0.5, 2.0),
                             n_rep = 500L) {
  d <- expand.grid(dgp = dgps, Tn = Ts, sigma = sigmas,
                   stringsAsFactors = FALSE)
  d$celda <- sprintf("%s_T%d_s%s", d$dgp, d$Tn, sub("\\.", "", format(d$sigma)))
  d$n_rep <- n_rep
  d
}

## -------------------------------------------------------------------------
## Una réplica: simula, evalúa, resume
## -------------------------------------------------------------------------
correr_replica <- function(dgp_id, Tn, sigma, rep_id,
                           modelos = MODELOS,
                           horizontes = HORIZONTES,
                           n_origenes = 10L,
                           semilla_base = SEMILLA_MAESTRA,
                           guardar_detalle = FALSE,
                           verbose = FALSE) {
  set.seed(semilla_base + rep_id * 1009L +
             sum(utf8ToInt(dgp_id)) * 7L + Tn + round(sigma * 100))
  sim <- simular_dgp(dgp_id, Tn, sigma, burn = BURN_IN)
  det <- evaluar_serie(sim$y, z = sim$z, s = sim$s,
                       modelos = modelos, horizontes = horizontes,
                       n_origenes = n_origenes, verbose = verbose)
  res <- resumir_replica(det)
  res$dgp   <- dgp_id
  res$Tn    <- Tn
  res$sigma <- sigma
  res$rep   <- rep_id
  res$etiqueta_dgp <- DGP_GRID[[dgp_id]]$etiqueta
  res$familia_dgp  <- DGP_GRID[[dgp_id]]$familia

  dm <- dm_vs_pivote(det, pivote = "M1")
  if (!is.null(dm)) { dm$dgp <- dgp_id; dm$Tn <- Tn; dm$sigma <- sigma; dm$rep <- rep_id }

  list(resumen = res, dm = dm,
       detalle = if (guardar_detalle) det else NULL,
       info = sim$info)
}

## -------------------------------------------------------------------------
## Corrida completa (secuencial o en paralelo con mclapply)
## -------------------------------------------------------------------------
correr_montecarlo <- function(diseno,
                              modelos = MODELOS,
                              horizontes = HORIZONTES,
                              n_origenes = 10L,
                              n_cores = 1L,
                              etiqueta = "mc",
                              guardar = TRUE,
                              verbose = TRUE) {
  t_ini <- Sys.time()
  tareas <- do.call(rbind, lapply(seq_len(nrow(diseno)), function(i) {
    data.frame(dgp = diseno$dgp[i], Tn = diseno$Tn[i], sigma = diseno$sigma[i],
               rep_id = seq_len(diseno$n_rep[i]), stringsAsFactors = FALSE)
  }))
  if (verbose) cat(sprintf("Tareas totales: %d (%d celdas)\n",
                           nrow(tareas), nrow(diseno)))

  aplicar <- if (n_cores > 1L && requireNamespace("parallel", quietly = TRUE)) {
    function(X, FUN) parallel::mclapply(X, FUN, mc.cores = n_cores)
  } else lapply

  salida <- aplicar(seq_len(nrow(tareas)), function(i) {
    tr <- tareas[i, ]
    if (verbose && (i %% 5 == 0 || i == 1))
      cat(sprintf("[%s] tarea %d/%d  %s T=%d sigma=%.1f rep=%d\n",
                  format(Sys.time(), "%H:%M:%S"), i, nrow(tareas),
                  tr$dgp, tr$Tn, tr$sigma, tr$rep_id))
    tryCatch(correr_replica(tr$dgp, tr$Tn, tr$sigma, tr$rep_id,
                            modelos = modelos, horizontes = horizontes,
                            n_origenes = n_origenes),
             error = function(e) list(resumen = NULL, dm = NULL,
                                      error = conditionMessage(e)))
  })

  resumen <- do.call(rbind, lapply(salida, `[[`, "resumen"))
  dm      <- do.call(rbind, lapply(salida, function(x) x$dm))
  errores <- vapply(salida, function(x) if (is.null(x$resumen)) 1L else 0L, integer(1))

  if (guardar) {
    write.csv(resumen, file.path(RUTA_RESULTADOS, paste0(etiqueta, "_resumen.csv")),
              row.names = FALSE)
    if (!is.null(dm))
      write.csv(dm, file.path(RUTA_RESULTADOS, paste0(etiqueta, "_dm.csv")),
                row.names = FALSE)
    saveRDS(list(resumen = resumen, dm = dm, diseno = diseno,
                 sesion = utils::sessionInfo(), fecha = Sys.time()),
            file.path(RUTA_RESULTADOS, paste0(etiqueta, "_completo.rds")))
  }
  if (verbose) cat(sprintf("Terminado en %.1f min | réplicas fallidas: %d\n",
                           as.numeric(difftime(Sys.time(), t_ini, units = "mins")),
                           sum(errores)))
  list(resumen = resumen, dm = dm, n_fallos = sum(errores))
}

## -------------------------------------------------------------------------
## Agregación entre réplicas: una fila por (celda, modelo, h)
## -------------------------------------------------------------------------
agregar_montecarlo <- function(resumen) {
  clave <- interaction(resumen$dgp, resumen$Tn, resumen$sigma,
                       resumen$modelo, resumen$h, drop = TRUE)
  out <- lapply(split(resumen, clave), function(g) {
    data.frame(
      dgp = g$dgp[1], etiqueta_dgp = g$etiqueta_dgp[1],
      familia_dgp = g$familia_dgp[1],
      Tn = g$Tn[1], sigma = g$sigma[1],
      modelo = g$modelo[1], etiqueta = g$etiqueta[1], familia = g$familia[1],
      h = g$h[1], n_rep = nrow(g),
      rmse_med = median(g$rmse, na.rm = TRUE),
      rmse_ee  = sd(g$rmse, na.rm = TRUE) / sqrt(sum(is.finite(g$rmse))),
      mae_med  = median(g$mae, na.rm = TRUE),
      mase_med = median(g$mase, na.rm = TRUE),
      crps_med = median(g$crps, na.rm = TRUE),
      picp80   = mean(g$picp80, na.rm = TRUE),
      picp95   = mean(g$picp95, na.rm = TRUE),
      mpiw95   = median(g$mpiw95, na.rm = TRUE),
      pinball  = median(g$pinball, na.rm = TRUE),
      tiempo_med = median(g$tiempo, na.rm = TRUE),
      tasa_ok  = mean(g$tasa_ok, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  ag <- do.call(rbind, out)
  ## RMSE relativo al benchmark ARIMA dentro de cada celda-horizonte
  ag$celda_h <- paste(ag$dgp, ag$Tn, ag$sigma, ag$h, sep = "|")
  base <- ag$rmse_med[ag$modelo == "M1"]
  names(base) <- ag$celda_h[ag$modelo == "M1"]
  ag$rmse_rel_arima <- ag$rmse_med / base[ag$celda_h]
  base_n <- ag$rmse_med[ag$modelo == "M0"]
  names(base_n) <- ag$celda_h[ag$modelo == "M0"]
  ag$theil_u <- ag$rmse_med / base_n[ag$celda_h]
  ag[order(ag$dgp, ag$Tn, ag$sigma, ag$h, ag$rmse_med), ]
}

## Proporción de rechazos del DM por celda
agregar_dm <- function(dm, alpha = 0.05) {
  if (is.null(dm)) return(NULL)
  clave <- interaction(dm$dgp, dm$Tn, dm$sigma, dm$modelo, dm$h, drop = TRUE)
  out <- lapply(split(dm, clave), function(g) {
    data.frame(dgp = g$dgp[1], Tn = g$Tn[1], sigma = g$sigma[1],
               modelo = g$modelo[1], h = g$h[1], n = nrow(g),
               prop_mejor_sig = mean(g$p < alpha & g$dm > 0, na.rm = TRUE),
               prop_peor_sig  = mean(g$p < alpha & g$dm < 0, na.rm = TRUE),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}
